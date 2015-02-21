c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c Externals required:
c--------------------------------------------------------------------
c == subroutine evaluateNonlinearResidual(ntot,x,res)
c
c      Evaluates nonlinear residual for the problem of interest
c
c        * ntot is the problem size
c        * x    is the newton state vector (input)
c        * res  is the nonlinear residual (output)
c
c == subroutine setupPC(ntot,xk,nwt_it)
c
c      Initializes preconditioner based on current Newton iterate
c
c == subroutine applyPC(ntot,rhs,sol,iout)
c
c      Applies preconditioner operator by solving P sol = rhs.
c
c        * rhs is the independent term (input)
c        * sol is the solution (output)
c        * iout, if > 0, indicates output level (input)
c
c == subroutine killPC
c
c      Deallocates dynamic storage space used in preconditioner
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c module nk
c ###################################################################
      module nk

      use grid_mpi

      use math

      use fixed_point_accelerator

      integer, parameter :: NOERR=0,EXACT=-1,NOCONV=1,CUTDT=2

      real(8) :: pdt,check

      real(8) :: epsmac=0d-15,roundoff

      integer,private :: jit=0

      logical :: pseudo_dt,prnt,is_residual=.true.

      type :: nk_options  !Explanation below
        integer      :: etak_meth=0
        integer      :: global_meth=0
        integer      :: ksmax=10
        integer      :: gmmax=10
        integer      :: nwt_max_it_acc=5
        integer      :: nwt_max_it_rej=10
        integer      :: gm_it_out
        integer      :: nwt_it_out
        real(8)      :: damp=1d0
        real(8)      :: pdt0=1d30
        real(8)      :: eta0=1d-1
        real(8)      :: atol=0d0
        real(8)      :: rtol=1d-4
        real(8)      :: stol=0d0
        real(8)      :: mf_eps=1d-6
        logical      :: mf_first_order=.true.
        character(2) :: krylov_method='gm'
        logical      :: vol_wgt=.true.

        !FPA setup
        type (fpa_state) :: fpas
        logical      :: fpa=.false.
        integer      :: fpa_it=5

        !Newton vectors
        real(8),pointer,dimension(:) :: xk,rk
      end type nk_options

      type(nk_options),save :: nk_conf

      contains

c     current_Newton_iteration
c     ###############################################################
      function current_Newton_iteration() result(nw_it)

      integer :: nw_it

      nw_it = jit
      
      end function current_Newton_iteration

c     jfnk
c     ###############################################################
      subroutine jfnk(neq,ntot,x,iguess,out,ierr,matvec)
c     ---------------------------------------------------------------
c     Performs Jacobian-free inexact Newton iteration on res(x) = 0, 
c     where res is calculated in 'evaluateNonlinearResidual' (provided
c     by user). 
c
c     Call parameters:
c       * neq: number of equations
c       * ntot: dimension of vector of unknowns.
c       * x: on input, initial guess; on output, solution.
c       * iguess: whether to use preconditioner to provide initial guess
c                 (iguess = 1).
c       * out: level of output:
c              out = 0 --> No output
c              out = 1 --> Newton convergence output
c              out = 2 --> Previous plus GMRES convergence output
c              out = 3 --> Previous plus residual information output
c              out > 4 --> Previous plus preconditioner output (controlled
c                          in preconditioner routine).
c       * ierr: error flag:
c              ierr = EXACT  --> Initial guess is exact solution
c              ierr = NOERR  --> No error
c              ierr = NOCONV --> No Newton convergence to prescribed tolerance in
c                                prescribed number of iterations (ntit_max_rej)
c              ierr = CUTDT  --> Newton converged, but took too many iterations
c                                (>  ntit_max_acc)
c       * matvec: optional matrix-vector multiplication routine
c
c     Configuration parameters (input via nk_options type):
c       * etak_meth: specifies method to determine inexact Newton forcing 
c           parameter. Currently:
c               etak_meth = 0 --> constant
c               etak_meth > 0 --> power law adaptive strategy
c       * damp  : damping parameter for Newton update
c       * global_meth: determines the globalization procedure:
c               global_meth = 0 --> no globalization
c               global_meth = 1 --> linesearch backtracking 
c               global_meth = 2 --> pseudo-transient
c       * pdt0: initial time step for pseudo-transient
c       * eta0: initial inexact Newton parameter (GMRES convergence tolerance)
c       * ksmax: dimension of krylov subspace for restarted GMRES(Ksmax)
c       * gmmax: maximum number of GMRES its. (restart if > ksmax)
c       * rtol: Newton relative convergence tolerance
c       * atol: Newton absolute convergence tolerance
c       * stol: Newton update convergence tolerance
c       * ntit_max_acc: maximum number of Newton its. to accept
c                       solution without subcycling time step
c       * ntit_max_rej: maximum number of Newton its. to reject solution
c       * gmit_out: on output, actual number of gmres its.
c       * ntit_out: on output, actual number of newton its.
c     ---------------------------------------------------------------

      implicit none   !For safe fortran

c     Call variables

      integer :: ntot,neq,out,iguess,ierr

      real(8) :: x(ntot)

      optional :: matvec

      INTERFACE
        subroutine matvec(nn,z,y)
        integer    :: nn
        real(8)    :: z(nn),y(nn)
        end subroutine matvec
      END INTERFACE

c     Local variables

      integer :: itk,i,j,ig,iout,etak_meth,ksmax,gmmax,ntit_max_acc
     .          ,ntit_max_rej,global,ntgl

      real(8) :: dxnorm,check_lim,ddx(ntot)
     .          ,eta0,atol,rtol,stol,damp,pdt0,leps

      real(8) :: f0,fkm,fk,fkp,fm,flimit,etak,etakm,theta,dampm

      logical :: convergence,failure,vol_wgt

c     Begin program

      prnt = (my_rank == 0)

      ntgl = ntot*npx*npy*npz

c     Find machine round-off

      if (epsmac == 0d0) epsmac = findRoundOff()

#if defined(petsc)
      leps = epsmac
      call MPI_Allreduce(leps,epsmac,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif

      roundoff = sqrt(1d0*ntgl)*epsmac

c     Initialize variables

      etak_meth    = nk_conf%etak_meth
      ksmax        = nk_conf%ksmax
      gmmax        = nk_conf%gmmax
      ntit_max_acc = nk_conf%nwt_max_it_acc
      ntit_max_rej = nk_conf%nwt_max_it_rej
      global       = nk_conf%global_meth

      eta0         = nk_conf%eta0
      damp         = nk_conf%damp
      pdt0         = nk_conf%pdt0
      atol         = nk_conf%atol
      rtol         = nk_conf%rtol
      stol         = nk_conf%stol

      vol_wgt      = nk_conf%vol_wgt

      nk_conf%gm_it_out  = 0
      nk_conf%nwt_it_out = 0

      if (ntit_max_acc > ntit_max_rej) then
        call pstop('nk'
     $            ,'Error in Newton input: ntit_max_acc > ntit_max_rej')
      endif

      if (nk_conf%fpa) then
        ntit_max_rej =  ntit_max_rej*2
        ntit_max_acc =  ntit_max_rej/2
      endif

      ierr     = 0
      itk      = 0

      if (.not.(nk_conf%fpa)) then
         if (out.ge.1 .and. prnt) write (*,200)
      else
         if (out.ge.1 .and. prnt) write (*,205)
      endif

      if (global == 2) then   !Activate pseudo-transient
        pseudo_dt = .true.
        call pstop('nk','Pseudo-t option not working')
      else
        pseudo_dt = .false.
cc        pdt0      = 1d30
      endif

      check_lim = 1d30

      pdt = pdt0

      if (atol == 0d0) atol = 1d3*epsmac !Set absolute tolerance to roundoff
      if (stol == 0d0) stol = 1d3*epsmac !Set update tolerance to roundoff

      convergence = .false.
      failure     = .false.

      jit = 0

c     Evaluate rhs and norms

      allocate (nk_conf%xk(ntot),nk_conf%rk(ntot))

      nk_conf%xk = x    !Save previous Newton state

      call evaluateNewtonResidual(ntot,x,nk_conf%rk)

      f0 = rms(nk_conf%rk)

      if (vol_wgt) f0 = ntgl*f0

      !Check if initial guess is exact, and if so exit Newton step
      if (f0.lt.atol) then
        ierr = -1
cc        deallocate(xx0)
        deallocate(nk_conf%xk,nk_conf%rk)
        return
      endif

c     Initial output

      if (.not.(nk_conf%fpa)) then
         if (out.ge.1 .and. prnt)
     $     write (*,210) 0,0d0,f0,1d0,1d0,pdt0,eta0,0
      else
         if (out.ge.1 .and. prnt)
     $     write (*,215) 0,0d0,f0,1d0,1d0,pdt0
      endif

      if (out.ge.3) then
        call dumpResidualNorms(neq,-nk_conf%rk,vol_wgt)
      endif

c     Start Newton iteration

      fk    = f0
      fkm   = f0
      flimit= atol + rtol*f0

      etakm = eta0
      etak  = eta0

      if (nk_conf%fpa) call fpa_create(nk_conf%fpas,nk_conf%rk
     .                                ,maxv=nk_conf%fpa_it)

      do jit = 1,ntit_max_rej

        !Determine inexact Newton constant
        if (.not.(nk_conf%fpa)) then
           call find_etak(fk,fkm,flimit,eta0,etak,etakm,etak_meth)
        endif

        !Setup preconditioner
        call setupPC(ntot,x,jit)

        !Jacobian-free solve
        iout = out - 1

        if (nk_conf%fpa) then
          call applyPC(ntot,-nk_conf%rk,ddx,iout-1)

          call fpa_correction(nk_conf%fpas,ddx,dp=dot2)
cc          ddx = -ddx
        else
          ddx = 0d0
          if (PRESENT(matvec)) then
            call gmresDriver(matvec
     .                    ,ntot,-nk_conf%rk,ddx,itk,etak,ksmax,gmmax
     .                    ,iguess,iout,ierr)
          else
            call gmresDriver(matrixFreeMatVec
     .                    ,ntot,-nk_conf%rk,ddx,itk,etak,ksmax,gmmax
     .                    ,iguess,iout,ierr)
          endif
        endif

        !Kill preconditioner
        call killPC

        !Update counters
        nk_conf%gm_it_out  = nk_conf%gm_it_out + itk
        nk_conf%nwt_it_out = jit

        !Check for error in GMRES
        if (.not.(nk_conf%fpa)) then
          if (ierr.eq.1) then
            if (prnt) write(*,*) 
     .            '   Exceeded maximum GMRES iterations (',gmmax,')'
            ierr = 0  !Use GMRES solution for Newton update regardless
          elseif (ierr.eq.-1) then !Got to atol in the Newton residual
            ierr = 0
            exit !Do not continue Newton iteration
          endif
        endif

        !Globalization procedure: damping coefficient
        if (global >= 1) call findDamping (ntot,x,ddx,etak,fk,damp)
cc        if (global == 1) call findDamping (ntot,x,ddx,etak,fk,damp)

        !Update solution
        x = x + damp*ddx

cc        call updateNewtonSolution(ntot,x,ddx,damp,dxnorm)

        nk_conf%xk = x   !Save current Newton state

        !Evaluate rhs and norms
        call evaluateNewtonResidual(ntot,x,nk_conf%rk)

        fkp = rms(nk_conf%rk)

        if (vol_wgt) fkp = ntgl*fkp

        !Check Newton convergence/failure
        check = fkp/f0

        dxnorm = rms(ddx)

        if (nk_conf%fpa) then
          if (out.ge.1 .and. prnt)
     $       write(*,215)jit,dxnorm,fkp,check,damp,pdt
        else
          if (out.ge.1 .and. prnt)
     $       write(*,210)jit,dxnorm,fkp,check,damp,pdt,etak,itk
        endif

        if (out.ge.3) then
          call dumpResidualNorms(neq,-nk_conf%rk,vol_wgt)
        endif

        convergence = (fkp <= flimit .or.dxnorm <= stol)
        failure     = (check > check_lim .or. jit == ntit_max_rej)
     .               .and.(.not.convergence)
cc        if (fkp <= flimit .or.dxnorm <= stol) then
cc          convergence = .true.
cc        elseif (check > check_lim .or. jit == ntit_max_rej) then
cc          failure = .true.
cc        endif

        if (convergence .or. failure) exit

        !Store magnitude of residual
        fkm   = fk
        fk    = fkp

        etakm = etak

        !Change pseudo-transient time step
cc        if (pseudo_dt) then
cccc          pdt = min(pdt/check,1d30)
cc          pdt = min(pdt0/sqrt(check),1d30)
cc          call evaluateNewtonResidual(ntot,x,rk)
cc
cccc          if (check < 1d-3) pdt = min(pdt0*jit,1d0)
cccc          pdt = min(pdt0/check**2,1d30)
cc        endif

      enddo       !End of Newton loop

      if (nk_conf%fpa) call fpa_destroy(nk_conf%fpas)

c     Check error in Newton convergence

      !No convergence: reject solution
      if (failure) then 
        if (jit == ntit_max_rej .and. prnt) write (*,220) check
        if (check > check_lim .and. prnt) write (*,230) check,check_lim
        if (damp  < 1d-4 .and. prnt)
     $       write (*,*) 'Damping parameter too small'
        ierr = 1
      !Accept solution, but warn user
      elseif (((jit-1).ge.ntit_max_acc).and.(ntit_max_acc > 0)) then 
        if (prnt) write (*,240) ntit_max_acc
        ierr = 2
      endif

c     End program

      deallocate (nk_conf%rk,nk_conf%xk)
cc      deallocate (xx0)

 200  format 
     .   (/,' New_it   Av_updt    Abs_res   Rel_res     Damping'
     .     ,'      pdt_n      eta_k     GMRES')
 205  format 
     .   (/,' FPA_it   Av_updt    Abs_res   Rel_res     Damping'
     .     ,'      pdt_n')
 210  format (i4,3x,1p,6e11.3,3x,i4)
 215  format (i4,3x,1p,5e11.3)
 220  format ('    Max newton its. exceeded; rel. residual: ',1p,1e10.2)
 230  format ('    Relative residual =',f7.2,' >',f7.2)
 240  format ('    Newton converged in too many iterations (>',i2,')')

      end subroutine jfnk

c     dumpResidualNorms
c     ###############################################################
      subroutine dumpResidualNorms(neq,b,vw)

c     ---------------------------------------------------------------
c     Calculates RMS of residuals
c     ---------------------------------------------------------------

      implicit none     !For safe fortran

c     Call variables

      integer    :: neq
      real(8)    :: b(:)
      logical    :: vw

c     Local variables

      integer    :: i,j,ntot,ntgl
      real(8)    :: ravg(neq),lravg(neq)

c     Begin program

      ntot = size(b)

      ntgl = ntot*npx*npy*npz

c     Calculate residuals

      ravg = 0d0

      do i = 1,ntot-neq,neq
        do j = 1,neq
          ravg(j) = ravg(j) + b(i - 1 + j)**2
        enddo
      enddo

c     Calculate magnitude of residuals

#if defined(petsc)
      if (.not.asm) then
        lravg = ravg
        call MPI_Allreduce(lravg,ravg,neq,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      endif
#endif

      ravg = sqrt(ravg/ntgl)

      if (vw) ravg=ravg*ntgl

      if (my_rank == 0) then
        do i = 1,neq
          write (*,320) i,ravg(i)
        enddo
      endif

c     End program

 320  format (' RMS residual eqn. ',i2,': ',1pe9.2)

      end subroutine dumpResidualNorms

ccc     updateNewtonSolution
ccc     ###############################################################
cc      subroutine updateNewtonSolution(ntot,x,ddx,damp,dxnorm)
cc      implicit none         !For safe fortran
ccc     ---------------------------------------------------------------
ccc     Updates solution in Newton iteration
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: ntot
cc      real(8)    :: x(ntot),ddx(ntot),damp,dxnorm
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc      x = x + damp*ddx
cc
cccc      dxnorm = sqrt(dot(ntot,ddx,ddx))
cc      dxnorm = rms(ddx)
cc
ccc     End program
cc
cc      end subroutine updateNewtonSolution

c     find_etak
c     ###############################################################
      subroutine find_etak (fk,fkm,flimit,eta0,etak,etakm,etak_meth)
      implicit none       !For safe fortran
c     ---------------------------------------------------------------
c     Finds inexact Newton forcing parameter, using two methods:
c     constant (etak_meth = 0) or power law adaptive strategy.
c
c     In call sequence:
c        * fk,fkm: residual norms at current (k) and previous (k-1)
c                  nonlinear iteration level
c        * flimit: nonlinear tolerance
c        * eta0: initial linear tolerance 
c        * etakm: linear tolerance at (k-1)
c        * etak (output): current iteration linear tolerance
c        * etak_meth: integer that selects method of determining etak
c                     (=0 -> constant; <>0 -> adaptive)
c     ---------------------------------------------------------------

c     Call variables

      integer    :: etak_meth
      real(8)    :: fk,fkm,flimit,eta0,etak,etakm

c     Local variables

      real(8)    :: gamm,alph

c     Begin program

      if (etak_meth.eq.0) then
        etak = eta0
      else
        gamm = .9
        alph = 1.5

        !For superlinear convergence, alph>1
        etak = gamm*(fk/fkm)**alph

        !First safeguard: avoid sharp decrease of etak
        etak = min(eta0,max(etak,gamm*etakm**alph))

        !Second safeguard: avoid oversolving
        etak = min(eta0,max(etak,gamm*flimit/fk))
      endif

c     End program

      end subroutine find_etak

c     findDamping
c     ###############################################################
      subroutine findDamping (ntot,x,ddx,etak,fk,damp)
      implicit none       !For safe fortran
c     ---------------------------------------------------------------
c     If global is true, uses linesearch backtracking to ensure
c     sufficient reduction in the Newton residual norm.
c     In call:
c       - ntot (input): vector dimension
c       - x (input): current Newton state
c       - ddx (input): Newton update
c       - etak (input/output): forcing for inexact Newton
c       - fk (input): current Newton residual
c       - damp (output): damping parameter
c
c     Taken from C. T. Kelley's "Iterative methods for Linear and
c     Nonlinear Equations", Sec. 8.3.
c     ---------------------------------------------------------------

c     Call variables

      integer    :: ntot
      real(8)    :: x(ntot),ddx(ntot),fk,etak,damp

c     Local variables

      integer    :: idamp
      real(8)    :: dampm,theta,etak0,dxnorm,fkp,fm,df,ddf
      real(8)    :: b(ntot),dummy(ntot)
      real(8)    :: alpha,sigma0,sigma1

c     Begin program

c     Initialize parameters

      damp = 1d0    !Damping
      etak0 = etak  !Inexact Newton forcing

      sigma0 = 0.1  !Lower limit in damping parameter decrease
      sigma1 = 0.5  !Damping parameter decrease for linesearch backtracking

      alpha = 1d-1  !Residual reduction parameter

c     Start quadratic linesearch backtracking

      do idamp = 1,100

        !Find residual norm |F(xk+dxk)| --> fkp
        dummy = x + damp*ddx
        call evaluateNewtonResidual(ntot,dummy,b)
cc        fkp = sqrt(dot(ntot,b,b))
        fkp = rms(b)

        !Check residual minimation condition
        if (fkp.lt.(1.-alpha*damp)*fk) then

          exit  !All is well

        else
          !New damping parameter
          dampm = sigma1*damp

          !Evaluate residual norm |F(xk+lambda*dxk)| --> fm
          dummy = x + dampm*ddx
          call evaluateNewtonResidual(ntot,dummy,b)
cc          fm = sqrt(dot(ntot,b,b))
          fm = rms(b)

          !Check convergence
          if (fm.lt.(1.-alpha*dampm)*fk) then   !Good enough

            damp = dampm
            exit

          else  !Quadratic minimization in lambda

            !Find out curvature of parabola
            ddf = 2./(damp-dampm)*((fkp-fk)/damp-(fm-fk)/dampm)

            if (ddf > 0) then !Parabola has a minimum

              df = 1./(damp-dampm)*(-dampm/damp*(fkp-fk)
     .                              +damp/dampm*(fm -fk))
              theta = -df/ddf

              damp = fmed(sigma0*damp,sigma1*damp,theta)  !Takes middle value

            else  !Try again

              damp = dampm

            endif

          endif
cc          if (damp < 1d-2) then
cc            damp = 1d0
cc            exit
cc          endif
cc          write (*,*) 'Residuals',fk,fm,fkp
cc          write (*,*) 'Damping parameters',dampm,damp
        endif
      enddo

c     End program

      end subroutine findDamping

c     gmresDriver
c     ################################################################
      subroutine gmresDriver(matvec,nn,b,x,itk,etak,ksmax,maxitgm,iguess
     .                      ,iout,ierr)
c     ----------------------------------------------------------------
c     This subroutine solves the linear system
c
c        A x = b
c
c     using the Generalized Minimal Residual (GMRES(K))
c     iteration algorithm with right pre-conditioning
c     Note that this subroutine calls the gmres subroutine from
c     the SPARSKIT package by Y. Saad, modified by L. Chacon.
c     ----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: itk,ierr,nn,ksmax,iguess,iout,maxitgm
      real(8)    :: b(nn),x(nn),etak

      external   :: matvec

c     Local variables

      integer    :: iiout
      real(8)    :: eps

c     Begin program

c     Compute the initial guess x

      if (iguess.eq.1) then
        iiout = iout - 2
        call applyPC(nn,b,x,iiout)
cc      else
cc        x = 0d0
      endif

c     Initialize the GMRES(k) algorithm

      eps = etak               !GMRES convergence tolerance

c     Call the preconditioned GMRES(k) algorithm

      if (nk_conf%krylov_method=='fg') then
        call fgmres(matvec,nn,b,x,eps,ksmax,maxitgm,iout,ierr,itk)
      else
        call  gmres(matvec,nn,b,x,eps,ksmax,maxitgm,iout,ierr,itk)
      endif

c     End program

      end subroutine gmresDriver

c     fgmres
c     #######################################################################
      subroutine fgmres(matvec,ntot,rhs,sol,eps,im,maxits,iout,ierr,its)
c     ----------------------------------------------------------------------*
c                                                                           *
c                    *** Preconditioned Flexible GMRES ***                  *
c                                                                           *
c     ----------------------------------------------------------------------*
c      This is a simple version of the right-preconditioned GMRES algorithm.*
c      The stopping criterion utilized is based simply on reducing the      *
c      residual norm by epsilon.                                            *
c     ----------------------------------------------------------------------*
c      parameters                                                           *
c     -----------                                                           *
c      on entry:                                                            *
c     ==========                                                            *
c                                                                           *
c      matvec== matrix-free routine                                         *
c      ntot  == integer. The dimension of the matrix.                       *
c      rhs   == real vector of length n containing the right hand side.     *
c               Destroyed on return.                                        *
c      sol   == real vector of length n containing an initial guess to the  *
c               solution on input. approximate solution on output           *
c      eps   == tolerance for stopping criterion. process is stopped        *
c               as soon as ( ||.|| is the euclidean norm):                  *
c               || current residual||/||initial residual|| <= eps           *
c      im    == size of krylov subspace.                                    *
c      maxits== maximum number of GMRES iterations allowed                  *
c      iout  == output unit number number for printing intermediate results *
c               if (iout .le. 0) nothing is printed out.                    *
c                                                                           *
c      on return:                                                           *
c     ==========                                                            *
c      sol   == contains an approximate solution (upon successful return).  *
c      ierr  == integer. Error message with the following meaning.          *
c               ierr = 0 --> successful return.                             *
c               ierr = 1 --> convergence not achieved in itmax iterations.  *
c               ierr =-1 --> the initial guess seems to be the exact        *
c                            solution (initial residual computed was zero)  *
c      its   == final number of GMRES iterations                            *
c                                                                           *
c     ----------------------------------------------------------------------*
c                                                                           *
c      work arrays:                                                         *
c     =============                                                         *
c      vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c               basis)                                                      *
c      zz    == work array of length  n x (im+1) (used to store the Arnoli  *
c               basis times the preconditioner operator (FGMRES))
c     ----------------------------------------------------------------------*
c      arnoldi size should not exceed im=50 in this version.                *
c     ----------------------------------------------------------------------*

      implicit none             !For safe fortran

c     Call variables

      integer    :: ntot,im,maxits,iout,ierr,its
      real(8)    :: rhs(ntot),sol(ntot),eps

      external   matvec

c     Local variables

      real(8)    :: hh(im+1,im), c(im), s(im), rs(im+1)
      real(8)    :: vv(ntot,im+1),zz(ntot,im+1),dsol(ntot),dum(ntot)
      real(8)    :: rold,ro,eps1,gam,t,xsvmx(im+1),xsvmn(im+1)
     .             ,smx_est,smn_est,smx_new,smn_new,snmx,snmn,csmx,csmn
      integer    :: i,j,i1,k,k1,ii,jj,rstrt,irstrt,precout

      logical    :: sing_mat=.false.

c     Begin program

      precout = iout - 2

      its = 0

c     Calculate magnitude of nonlinear residual

      rold = sqrt(dot2(rhs,rhs))

      eps1=eps*rold

      if (rold < roundoff) then
        ierr = -1
        return
      endif

c     Compute initial residual vector

      call matvec(ntot,sol,vv(:,1))

      vv(:,1) = rhs(:) - vv(:,1)

c     Calculate number of restarting loops

      rstrt = min(maxits/im + 1,maxits)

cc      if (maxits > 0) then
cc        rstrt = maxits/im + 1
cc      else
cc        rstrt = 0
cc      endif

c     Restarted FGMRES loop

      do irstrt = 1,rstrt

        rs = 0d0

        ro = sqrt(dot2(vv(:,1),vv(:,1)))

        if (iout .gt. 0 .and. prnt) then
          if (its .eq. 0) then
            write (*,199) its, ro, ro/rold
          else
            write (*,*) 'Restarting FGMRES... (',irstrt-1,')'
          endif
        endif

        if (ro.le.eps1) exit

        t = 1d0/ro
        vv(:,1) = vv(:,1)*t

c       Initialize 1-st term  of rhs of hessenberg system

        rs(1) = ro

c       FGMRES iteration

        do i = 1,im

          its = its + 1
          i1  = i + 1

          call applyPC(ntot,vv(:,i),zz(:,i),precout)

          call matvec(ntot,zz(:,i),vv(:,i1))

c       Modified gram - schmidt.

          do j=1,i
            t = dot2(vv(:,j),vv(:,i1))
            hh(j,i) = t
            vv(:,i1) = vv(:,i1)-t*vv(:,j)
          enddo

          t = sqrt(dot2(vv(:,i1),vv(:,i1)))
          hh(i1,i) = t

          if (t.ne.0d0) then
            t = 1d0/t
            vv(:,i1) = vv(:,i1)*t
          endif

c       Done with modified Gram-Schimdt and arnoldi step
c       Now update factorization of hh
c       Perform previous transformations on i-th column of h

          if (i .gt. 1) then
            do k=2,i
              k1 = k-1
              t = hh(k1,i)
              hh(k1,i) =  c(k1)*t + s(k1)*hh(k,i)
              hh(k ,i) = -s(k1)*t + c(k1)*hh(k,i)
            enddo
          endif

          gam = sqrt(hh(i,i)*hh(i,i) + hh(i1,i)*hh(i1,i))

c       If gamma is zero then any small value will do
c       Will affect only residual estimate

          if (gam.eq.0d0) gam = epsmac

c       Get next plane rotation

          c(i)   = hh(i ,i)/gam
          s(i)   = hh(i1,i)/gam

          hh(i,i)= c(i)*hh(i,i) + s(i)*hh(i1,i)
          hh(i1,i) = 0d0

c       Estimate condition number of hh

          if (sing_mat) then
            if (i > 1) then
              !Max Sing. value
              call dlaic1(1,i-1,xsvmx(1:i-1),smx_est,hh(1:i-1,i),hh(i,i)
     .                   ,smx_new,snmx,csmx)

              !Min Sing. value
              call dlaic1(2,i-1,xsvmn(1:i-1),smn_est,hh(1:i-1,i),hh(i,i)
     .                   ,smn_new,snmn,csmn)

cc              write (*,*) 'Mat column=',hh(i,1:i-1)
cc              write (*,*) 'Gamma=',hh(i,i)
cc              write (*,*) 'Max sv ',xsvmx(1:i-1)
cc     .                   ,'norm=',sqrt(sum(xsvmx(1:i-1)**2))
cc     .                   ,'sine=',snmx,'cos=',csmx
cc              write (*,*) 
cc              write (*,*) 'Min sv ',xsvmn(1:i-1)
cc     .                   ,'norm=',sqrt(sum(xsvmn(1:i-1)**2))
cc     .                   ,'sine=',snmn,'cos=',csmn
cc              write (*,*) 

              !Condition number
              write (*,*)
              write (*,*) 'FGMRES diag: Max. sing. val. =',smx_new
              write (*,*) 'FGMRES diag: Min. sing. val. =',smn_new
              if (smn_new /= 0d0) then
                write (*,*) 'FGMRES diag: Condition number='
     .                      ,smx_new/smn_new
              endif

              !Update quantities
              smx_est = smx_new
              smn_est = smn_new

              xsvmx(1:i-1) = snmx*xsvmx(1:i-1)
              xsvmx(i)     = csmx
              xsvmn(1:i-1) = snmn*xsvmn(1:i-1)
              xsvmn(i)     = csmn
            else
              smx_est = hh(1,1)  !Initial SV is single element
              smn_est = hh(1,1)  !Initial SV is single element
              xsvmx(1) = 1d0
              xsvmn(1) = 1d0
            endif
          endif

c       Determine residual norm and test for convergence

          rs(i1) = -s(i)*rs(i)
          rs(i)  =  c(i)*rs(i)

          ro = abs(rs(i1))

          if (iout .gt. 0 .and. prnt) then
            write(*, 199) its, ro, ro/rold
          endif

          !The i=im condition below is necessary because i is used afterwards
          if (ro <= eps1.or.its >= maxits.or.i.eq.im) exit

        enddo

c     Now compute solution. First solve upper triangular system

        rs(i) = rs(i)/hh(i,i)
        do ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do j=k1,i
            t = t-hh(k,j)*rs(j)
          enddo
          rs(k) = t/hh(k,k)
        enddo

c     Form linear combination of zz=P*vv to get solution

        t = rs(1)
        dsol(:) = zz(:,1)*t

        do j=2, i
          t = rs(j)
          dsol(:) = dsol(:) + t*zz(:,j)
        enddo

c     Call preconditioner

        vv(:,im+1) = dsol(:)

c     Update solution

        sol = sol + dsol

c diag *****
cc        if (sing_mat) then
cc          write (*,*)
cc          call matvec(ntot,sol,dum)
cc          write (*,*) 'FGMRES diag: Residual=',sqrt(sum((dum-rhs)**2))
cc          write (*,*) 'FGMRES diag: Sol. update=',sqrt(sum(sol**2))
cccc          write (*,*) 'Solution coeffs=',rs
cc          pause
cc        endif
c diag *****

c     Check convergence and restart outer loop if necessary

        if (ro <= eps1) then
          ierr = 0
          exit
        endif

        if (its >= maxits) then
          ierr = 1
          exit
        endif

c     Else compute residual vector and continue

        do j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj)   =  c(jj-1)*rs(jj)
        enddo

        do j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          vv(:,1) = vv(:,1) + t*vv(:,j)
        enddo

c     Restart outer loop

      enddo

c     End program

 199  format(' FGMRES its =', i4,';  Res =',1p,d10.4
     .      ,';  Ratio =', d10.4)

      end subroutine fgmres

c     gmres
c     ######################################################################
      subroutine gmres(matvec,ntot,rhs,sol,eps,im,maxits,iout,ierr,its)
c     ----------------------------------------------------------------------*
c                                                                           *
c                    *** Preconditioned GMRES ***                  *
c                                                                           *
c     ----------------------------------------------------------------------*
c      This is a simple version of the right-preconditioned GMRES algorithm.*
c      The stopping criterion utilized is based simply on reducing the      *
c      residual norm by epsilon.                                            *
c     ----------------------------------------------------------------------*
c      parameters                                                           *
c     -----------                                                           *
c      on entry:                                                            *
c     ==========                                                            *
c                                                                           *
c      matvec== matrix-free routine                                         *
c      ntot  == integer. The dimension of the matrix.                       *
c      rhs   == real vector of length n containing the right hand side.     *
c               Destroyed on return.                                        *
c      sol   == real vector of length n containing an initial guess to the  *
c               solution on input. approximate solution on output           *
c      eps   == tolerance for stopping criterion. process is stopped        *
c               as soon as ( ||.|| is the euclidean norm):                  *
c               || current residual||/||initial residual|| <= eps           *
c      im    == size of krylov subspace.                                    *
c      maxits== maximum number of GMRES iterations allowed                  *
c      iout  == output unit number number for printing intermediate results *
c               if (iout .le. 0) nothing is printed out.                    *
c                                                                           *
c      on return:                                                           *
c     ==========                                                            *
c      sol   == contains an approximate solution (upon successful return).  *
c      ierr  == integer. Error message with the following meaning.          *
c               ierr = 0 --> successful return.                             *
c               ierr = 1 --> convergence not achieved in itmax iterations.  *
c               ierr =-1 --> the initial guess seems to be the exact        *
c                            solution (initial residual computed was zero)  *
c      its   == final number of GMRES iterations                            *
c                                                                           *
c     ----------------------------------------------------------------------*
c                                                                           *
c      work arrays:                                                         *
c     =============                                                         *
c      vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c               basis)                                                      *
c     ----------------------------------------------------------------------*
c      arnoldi size should not exceed im=50 in this version.                *
c     ----------------------------------------------------------------------*

      implicit none             !For safe fortran

c     Call variables

      integer    :: ntot,im,maxits,iout,ierr,its
      real(8)    :: rhs(ntot),sol(ntot),eps

      external   matvec

c     Local variables

      real(8)    :: hh(im+1,im), c(im), s(im), rs(im+1)
      real(8)    :: vv(ntot,im+1),dsol(ntot)
      real(8)    :: rold,ro,eps1,gam,t,xsvmx(im+1),xsvmn(im+1)
     .             ,smx_est,smn_est,smx_new,smn_new,snmx,snmn,csmx,csmn
      integer    :: i,j,i1,k,k1,ii,jj,rstrt,irstrt,precout

      logical    :: sing_mat=.false.

c     Begin program

      precout = iout - 2

      its = 0

c     Calculate magnitude of nonlinear residual

      rold = sqrt(dot2(rhs,rhs))

      eps1=eps*rold

      if (rold < roundoff) then
        ierr = -1
        return
      endif

c     Compute initial residual vector

      call matvec(ntot,sol,vv(:,1))

      vv(:,1) = rhs(:) - vv(:,1)

c     Calculate number of restarting loops

      if (maxits > 0) then
        rstrt = maxits/im + 1
      else
        rstrt = 0
      endif

c     Restarted FGMRES loop

      do irstrt = 1,rstrt

        ro = sqrt(dot2(vv(:,1),vv(:,1)))

        if (iout .gt. 0 .and. prnt) then
          if (its .eq. 0) then
            write (*,199) its, ro, ro/rold
          else
            write (*,*) 'Restarting GMRES... (',irstrt-1,')'
          endif
        endif

        if (ro.le.eps1) exit

        t = 1d0/ro
        vv(:,1) = vv(:,1)*t

c     Initialize 1-st term  of rhs of hessenberg system

        rs(1) = ro

c     FGMRES iteration

        do i = 1,im

          its = its + 1
          i1  = i + 1

          call applyPC(ntot,vv(1,i),rhs,precout)

          call matvec(ntot,rhs,vv(1,i1))

c       Modified gram - schmidt.

          do j=1,i
            t = dot2(vv(:,j),vv(:,i1))
            hh(j,i) = t
            vv(:,i1) = vv(:,i1)-t*vv(:,j)
          enddo

          t = sqrt(dot2(vv(:,i1),vv(:,i1)))
          hh(i1,i) = t

          if (t.ne.0d0) then
            t = 1d0/t
            vv(:,i1) = vv(:,i1)*t
          endif

c       Done with modified Gram-Schimdt and arnoldi step
c       Now update factorization of hh
c       Perform previous transformations on i-th column of h

          if (i .gt. 1) then
            do k=2,i
              k1 = k-1
              t = hh(k1,i)
              hh(k1,i) =  c(k1)*t + s(k1)*hh(k,i)
              hh(k ,i) = -s(k1)*t + c(k1)*hh(k,i)
            enddo
          endif

          gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)

c       If gamma is zero then any small value will do
c       Will affect only residual estimate

          if (gam.eq.0d0) gam = epsmac

c       Get next plane rotation

          c(i)   = hh(i ,i)/gam
          s(i)   = hh(i1,i)/gam

          hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
          hh(i1,i) = 0d0

c       Estimate condition number of hh

          if (sing_mat) then
            if (i > 1) then
              !Max Sing. value
              call dlaic1(1,i-1,xsvmx(1:i-1),smx_est,hh(1:i-1,i),hh(i,i)
     .                   ,smx_new,snmx,csmx)

              !Min Sing. value
              call dlaic1(2,i-1,xsvmn(1:i-1),smn_est,hh(1:i-1,i),hh(i,i)
     .                   ,smn_new,snmn,csmn)

cc              write (*,*) 'Mat column=',hh(i,1:i-1)
cc              write (*,*) 'Gamma=',hh(i,i)
cc              write (*,*) 'Max sv ',xsvmx(1:i-1)
cc     .                   ,'norm=',sqrt(sum(xsvmx(1:i-1)**2))
cc     .                   ,'sine=',snmx,'cos=',csmx
cc              write (*,*) 
cc              write (*,*) 'Min sv ',xsvmn(1:i-1)
cc     .                   ,'norm=',sqrt(sum(xsvmn(1:i-1)**2))
cc     .                   ,'sine=',snmn,'cos=',csmn
cc              write (*,*) 

              !Condition number
              write (*,*) 'Max. sing. val. =',smx_new
              write (*,*) 'Min. sing. val. =',smn_new
              if (smn_new /= 0d0) then
                write (*,*) 'Condition number=',smx_new/smn_new
              endif

              !Update quantities
              smx_est = smx_new
              smn_est = smn_new

              xsvmx(1:i-1) = snmx*xsvmx(1:i-1)
              xsvmx(i)     = csmx
              xsvmn(1:i-1) = snmn*xsvmn(1:i-1)
              xsvmn(i)     = csmn
            else
              smx_est = hh(1,1)  !Initial SV is single element
              smn_est = hh(1,1)  !Initial SV is single element
              xsvmx(1) = 1d0
              xsvmn(1) = 1d0
            endif
          endif

c       Determine residual norm and test for convergence

          rs(i1) = -s(i)*rs(i)
          rs(i)  =  c(i)*rs(i)
          ro = abs(rs(i1))

          if (iout .gt. 0 .and. prnt) then
            write(*, 199) its, ro, ro/rold
          endif

          !The i=im condition below is necessary because i is used afterwards
          if (ro <= eps1.or.its >= maxits.or.i.eq.im) exit

        enddo

c     Now compute solution. First solve upper triangular system

        rs(i) = rs(i)/hh(i,i)
        do ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do j=k1,i
            t = t-hh(k,j)*rs(j)
          enddo
          rs(k) = t/hh(k,k)
        enddo

c     Form linear combination of zz=P*vv to get solution

        t = rs(1)
        dsol(:) = vv(:,1)*t

        do j=2, i
          t = rs(j)
          dsol(:) = dsol(:) + t*vv(:,j)
        enddo

c     Call preconditioner

       call applyPC(ntot,dsol,vv(1,im+1),precout)

c     Update solution

        sol(:) = sol(:) + vv(:,im+1)

c     Check convergence and restart outer loop if necessary

        if (ro <= eps1) then
          ierr = 0
          exit
        endif

        if (its >= maxits) then
          ierr = 1
          exit
        endif

c     Else compute residual vector and continue

        do j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj)   =  c(jj-1)*rs(jj)
        enddo

        do j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          vv(:,1) = vv(:,1) + t*vv(:,j)
        enddo

c     Restart outer loop

      enddo

c     End program

 199  format(' GMRES its =', i4,';  Res =',1p,d10.2
     .      ,';  Ratio =', d10.2)

      end subroutine gmres

c     matrixFreeMatVec
c     #############################################################
      subroutine matrixFreeMatVec(nn,z,y)

c     -------------------------------------------------------------
c     This subroutine computes
c
c       y = J(xk)*z,
c
c     where J(xk) is the Jacobian matrix at the previous Newton iterate.
c     A finite difference approximation is used to compute this product,
c
c       y = (F(xk+e*z) - F(xk))/e
c
c     where e is some small perturbation constant
c     -------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      integer :: nn
      real(8) :: z(nn),y(nn)

c     Local variables

      integer :: i
      real(8) :: dummy(nn),yk(nn),pert,ipert,modz,modx,xdotz

c     Begin program

      modz = dot2(z,z)

c     Calculate J.x

      if (sqrt(modz) == 0d0) then  !Failsafe for the case z=0

        y = 0d0

      else

c       Calculate difference parameter

        modx  = dot2(nk_conf%xk,nk_conf%xk)  !x is grid based, and hence we compute volume-weighed norm
        xdotz = dot2(z         ,nk_conf%xk)  !x is grid based, and hence we compute volume-weighed norm

        if (modx < roundoff) then
          pert = nk_conf%mf_eps/sqrt(modz)
        else
          pert = nk_conf%mf_eps*sqrt(modx)/sqrt(modz)
        endif
cc        pert=nk_conf%mf_eps*(sqrt(modz)+abs(xdotz))/modz*sign(1d0,xdotz)
cc        if(my_rank == 0) write (*,*) modz,modx,xdotz,pert

c       Perturb state variables x + eps z --> dummy

        dummy = nk_conf%xk + pert*z

c       Nonlinear function evaluation --> y

        call evaluateNewtonResidual(nn,dummy,y)

c       Compute the product J.x

cc        ipert = 1d0/pert

        if (nk_conf%mf_first_order) then  ! using a first-order matrix-free approx
          y = (y-nk_conf%rk)/pert
        else                      ! Second order formula
          dummy = nk_conf%xk - pert*z
          call evaluateNewtonResidual(nn,dummy,yk)
          y = 0.5*(y-yk)/pert
        endif

      endif

c     End

      end subroutine matrixFreeMatVec

c     evaluateNewtonResidual
c     ###############################################################
      subroutine evaluateNewtonResidual(ntot,x,f)

c     ---------------------------------------------------------------
c     Calculates nonlinear residuals. 
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: ntot
      real(8)  :: x(ntot),f(ntot)

c     Local variables

      real(8)  :: invpdt

c     Begin program

c     Evaluate nonlinear residual

      call evaluateNonlinearResidual(ntot,x,f)

c     Add pseudo-transient term

cc      if (pseudo_dt) then
cc        invpdt = 1d0/pdt
cc        if (abs(invpdt) < 1d-5) invpdt = 0d0
cc        f = (x - xk )*invpdt + f  !This doesn't really work for vol-weighed residuals!!!
cccc        f = (x - xx0)*invpdt + f
cc        if (nk_conf%vol_wgt) then
cc          write (*,*) 'WARNING: psedo-transient term only works for'//
cc     .               ' non-vol-weighed residuals'
cc        endif
cc      endif

c     End program

      end subroutine evaluateNewtonResidual

ccc     volumeWeighVec
ccc     ###############################################################
cc      subroutine volumeWeighVec(ntot,x)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates nonlinear residuals. 
ccc     ---------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer  :: ntot
cc      real(8)  :: x(ntot)
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc      do k = 1,nzd
cc        do j = 1,nyd
cc          do i = 1,nxd
cc
cc            call getMGmap(i,j,k,1,1,1,ig,jg,kg)
cc
cc            ii = neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))
cc
cc            do ieq=1,neqd
cc              f(ii+ieq) = (varray%array_var(ieq)%array(i,j,k)
cc     .                    -   u_n%array_var(ieq)%array(i,j,k))
cc     .                    *one_over_dt(ieq)
cc     .                    + (1d0-cnf(ieq))*f   (ii+ieq)
cc     .                    +      cnf(ieq) *fold(ii+ieq)
cc     .                    -                fsrc(ii+ieq)
cc            enddo
cc
cc            !Do not include Jacobian here to allow for moving grid cases
cccc            dvol = grid_params%dxh(ig)
cccc     .            *grid_params%dyh(jg)
cccc     .            *grid_params%dzh(kg)
cccc            f(ii+1:ii+neqd) = f(ii+1:ii+neqd)*dvol
cc
cc            if (vol_wgt) 
cc     .        f(ii+1:ii+neqd) = f(ii+1:ii+neqd)
cc     .                         *gmetric%grid(1)%dvol(i,j,k)
cc
cc          enddo
cc        enddo
cc      enddo
cc
ccc     End program
cc
cc      end subroutine volumeWeighVec
cc
ccc     fmedval
ccc     ###############################################################
cc      real(8) function fmedval(p1,p2,p3)
cc      implicit none                !For safe fortran
ccc     ---------------------------------------------------------------
ccc     This function computes intermediate value of p1, p2, p3.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      real(8) :: p1,p2,p3
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc      fmedval = min( max(p1,p2) , max( p3,min(p1,p2) ) )
cc
ccc     End
cc
cc      end function fmedval

      end module nk

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
c == subroutine setupPreconditioner(ntot,x)
c
c      Initializes preconditioner
c
c == subroutine applyPreconditioner(ntot,rhs,sol,iout)
c
c      Applies preconditioner operator by solving P sol = rhs.
c
c        * rhs is the independent term (input)
c        * sol is the solution (output)
c        * iout, if > 0, indicates output level (input)
c
c == subroutine killPreconditioner
c
c      Deallocates variables used in preconditioner
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


c module newtonGmres
c ###################################################################
      module newton_gmres

      double precision, allocatable, dimension(:):: res

      end module newton_gmres

c newtonGmres
c####################################################################
      subroutine newtonGmres(neq,ntot,x,method,global,tolgm,kmax
     .           ,gmit_max,tolnewt,ntit_max_acc,ntit_max_rej,gmit_out
     .           ,ntit_out,iguess,out,ierr)
c--------------------------------------------------------------------
c     Performs Jacobian-free inexact Newton iteration on b(x) = 0, 
c     where b is calculated in 'evaluateNonlinearResidual' (provided
c     by user). 
c
c     Call parameters:
c       * neq: number of equations
c       * ntot: dimension of vector of unknowns.
c       * x: on input, initial guess; on output, solution.
c       * method: specifies method to determine inexact Newton forcing 
c                 parameter. Currently:
c               method = 0 --> constant
c               method = 1 --> unused
c               method = 2 --> power law adaptive strategy
c       * global: if true, uses linesearch backtracking as globalization
c                 procedure
c       * tolgm: GMRES convergence tolerance
c       * kmax: dimension of krylov subspace for restarted GMRES(Kmax)
c       * gmit_max: maximum number of GMRES its. (restart if > kmax)
c       * tolnewt: Newton convergence tolerance
c       * ntit_max_acc: maximum number of Newton its. to accept
c                       solution without subcycling time step
c       * ntit_max_rej: maximum number of Newton its. to reject solution
c       * gmit_out: on output, actual number of gmres its.
c       * ntit_out: on output, actual number of newton its.
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
c              ierr =-1 --> Initial guess is exact solution
c              ierr = 0 --> No error
c              ierr = 1 --> No Newton convergence to prescribed tolerance in
c                           prescribed number of iterations (ntit_max_rej)
c              ierr = 2 --> Newton converged, but took too many iterations
c                           (>  ntit_max_acc)
c--------------------------------------------------------------------

      use newton_gmres

      implicit none   !For safe fortran

c Call variables

      integer          :: ntot,neq,method,kmax,gmit_max,ntit_max_acc
     .                   ,ntit_max_rej,gmit_out,ntit_out,out,iguess,ierr

      double precision :: x(ntot),tolgm,tolnewt

      logical          :: global

c Local variables

      integer*4   itk,jit,i,j,ig,idamp,iout

      real*8      dxavg,check_lim,check,damp,residuals(neq)

      real*8      b(ntot),ddx(ntot)

      real*8      f0,fkm,fk,fkp,fm,etak,etakm,theta,dampm,gamm,alph

c Begin program

      if (ntit_max_acc > ntit_max_rej) then
        write (*,*) 'Error in Newton input: ntit_max_acc > ntit_max_rej'
        write (*,*) 'Aborting...'
        stop
      endif

c Initialize

      allocate (res(ntot))

      ierr     = 0
      gmit_out = 0
      itk      = 0

      check_lim = 1d1

      if (out.ge.1) write (*,200)

c Evaluate rhs and norms

      call evaluateNonlinearResidual(ntot,x,res)

      b = -res

      f0 = sqrt(sum(b*b))

c Start Newton iteration

      fk    = f0
      fkm   = f0
      etakm = tolgm
      etak  = tolgm

      do jit = 1,ntit_max_rej

c     Setup preconditioner

        call setupPreconditioner(ntot,x)

c     Jacobian-free solve

        iout = out - 1

        call gmresDriver(ntot,x,b,ddx,itk,etak,kmax,gmit_max
     .                  ,iguess,iout,ierr)

c     Kill preconditioner

        call killPreconditioner

c     Update counters

        gmit_out = gmit_out + itk
        ntit_out = jit

c     Check error in gmresDriver

        if (ierr.eq.1) then
          write(*,*) 
     .          '   Exceeded maximum GMRES iterations (',gmit_max,')'
          ierr = 0  !Use GMRES solution for Newton update regardless
        elseif (ierr.eq.-1) then
          write(*,*)
          write(*,*) '   Initial guess seems to be exact solution'
          write(*,*) '   Aborting...'
          exit
        endif

c     Globalization procedure

c       Determine inexact Newton constant

        call find_etak (fk,fkm,tolgm,etak,etakm,method)

c       Calculate damping coefficient

        call findDamping (ntot,x,ddx,etak,fk,damp,global)

c     Update solution (x = x + ddx)

        call updateNewtonSolution(x,ddx,ntot,damp,dxavg)

c     Evaluate rhs and norms

        call evaluateNonlinearResidual(ntot,x,res)

        b = -res

        fkp = sqrt(sum(b*b))

c     Check Newton convergence

        check = fkp/f0

        if (out.ge.1) write (*,210) jit,dxavg,check,damp,itk

        if (out.ge.3) then
          call calculateResidualNorms(neq,ntot,b,residuals)
          do i = 1,neq
            write (*,320) i,residuals(i)
          enddo
        endif

        if(fkp.lt.(1d-11 + tolnewt*f0).or.check.gt.check_lim) exit
cc        if(fkp.lt.(1d-11 + tolnewt*f0)) exit

c     Store magnitude of residual

        fkm   = fk
        fk    = fkp

        etakm = etak

      enddo       !End of Newton loop

c Check error in Newton convergence

cc      if (jit.gt.ntit_max_rej) then          !No convergence: reject solution
cc        write (*,220) check
cc        ierr = 1
      if (jit.gt.ntit_max_rej.or.check.gt.check_lim) then !No convergence: reject solution
        if (jit.gt.ntit_max_rej) write (*,220) check
        if (check.gt.check_lim) write (*,230) check,check_lim
        ierr = 1
      elseif ((jit-1).gt.ntit_max_acc) then  !Accept solution, but warn user
        write (*,240) ntit_max_acc
        ierr = 2
      endif

c End program

      deallocate (res)

      return

 200  format (/,' New_it   Av_updt    Rel_res    Damping  GMRES')
 210  format (i4,3x,1p3e11.3,i4)
 220  format ('    Max newton its. exceeded; rel. residual: ',1p1e10.2)
 230  format ('    Relative residual =',f7.2,' >',f7.2)
 240  format ('    Newton converged in too many iterations (>',i2,')')
 320  format ('Residual  eqn. ',i2,': ',1pe9.2)

      end

c calculateResidualNorms
c####################################################################
      subroutine calculateResidualNorms(neq,ntot,b,ravg)

c--------------------------------------------------------------------
c     Calculates norms of residuals
c--------------------------------------------------------------------

      implicit none     !For safe fortran

c Call variables

      integer ::  neq,ntot
      real*8      b(ntot),ravg(neq)

c Local variables

      integer*4   i,j

c Begin program

c Calculate residuals

      ravg = 0d0

      do i = 1,ntot-neq,neq
        do j = 1,neq
          ravg(j) = ravg(j) + b(i - 1 + j)**2
        enddo
      enddo

c    Calculate magnitude of residuals

      ravg = sqrt(ravg)

c End program

      return
      end

c updateNewtonSolution
c####################################################################
      subroutine updateNewtonSolution(x,ddx,ntot,damp,dxavg)
      implicit none         !For safe fortran
c--------------------------------------------------------------------
c    Updates solution in Newton iteration
c--------------------------------------------------------------------

c Call variables

      integer*4   ntot
      real*8      x(ntot),ddx(ntot),damp,dxavg

c Local variables

c Begin program

      x = x + damp*ddx

      dxavg = damp/dfloat(ntot)*sqrt(sum(ddx*ddx))

c End program

      return
      end

c find_etak
c####################################################################
      subroutine find_etak (fk,fkm,tolgm,etak,etakm,method)
      implicit none       !For safe fortran
c--------------------------------------------------------------------
c     Finds inexact Newton forcing parameter, using two methods:
c     constant (method = 0) or power law adaptive strategy.
c--------------------------------------------------------------------

c Call variables

      integer*4   method
      real*8      fk,fkm,tolgm,etak,etakm

c Local variables

      real*8      gamm,alph

c Begin program

      if (method.eq.0) then
        etak = tolgm
      else
        gamm = .5
        alph = 1.5
        etak = gamm*(fk/fkm)**alph
        etak = min(tolgm,max(etak,gamm*etakm**alph))
      endif

c End program

      return
      end

c findDamping
c####################################################################
      subroutine findDamping (ntot,x,ddx,etak,fk,damp,global)
      implicit none       !For safe fortran
c--------------------------------------------------------------------
c     If global is true, uses linesearch backtracking to ensure
c     sufficient reduction in the Newton residual norm.
c--------------------------------------------------------------------

c Call variables

      logical ::  global
      integer ::  ntot
      double precision ::  x(ntot),ddx(ntot),fk,etak,damp

c Local variables

      integer ::  idamp
      double precision :: dampm,theta,dxavg,b(ntot),dummy(ntot),fkp,fm

c External

      double precision :: fmedval
      external            fmedval

c Begin program

      damp = 1d0

      if (global) then
        do idamp = 1,100
          dummy = x
          call updateNewtonSolution(dummy,ddx,ntot,damp,dxavg)
          call evaluateNonlinearResidual(ntot,dummy,b)
          fkp = sqrt(sum(b*b))
          if (fkp.gt.(1.-1e-4*(1.-etak))*fk) then
            dampm = .5*damp
            dummy = x
            call updateNewtonSolution(dummy,ddx,ntot,dampm,dxavg)
            call evaluateNonlinearResidual(ntot,dummy,b)
            fm = sqrt(sum(b*b))
            theta = .5*(1.5*fk + 0.5*fkp - 2.*fm)/(fk + fkp - 2.*fm)
            theta = fmedval(1d-1,8d-1,theta)
            damp = damp*theta
cc            write (*,*) 'New damping parameter',damp
            etak = 1.- theta*(1.-etak)
          else
            exit
          endif
        enddo
      endif

c End program

      return
      end

c gmresDriver
c ####################################################################
      subroutine gmresDriver(nn,x0,b,x,itk,etak,kmax,maxitgm,iguess
     .                      ,iout,ierr)
c --------------------------------------------------------------------
c   This subroutine solves the linear system
c
c      A(x0) x = b
c
c   using the Generalized Minimal Residual (GMRES(K))
c   iteration algorithm with MG right pre-conditioning
c   Note that this subroutine calls the gmres subroutine from
c   the SPARSKIT package by Y. Saad, modified by L. Chacon.
c --------------------------------------------------------------------

      implicit none    !For safe fortran

c Call variables

      integer*4   itk,ierr,nn,kmax,iguess,iout,maxitgm
      real*8      b(nn),x(nn),x0(nn),etak

c Local variables

      integer*4   iiout
      real*8      vv(nn,kmax+1),zz(nn,kmax+1),eps

c Begin program

c Compute the initial guess

      if(iguess.eq.1) then
        iiout = iout - 2
        call applyPreconditioner(nn,b,x,iiout)
      else
        x = 0d0
      endif

c Initialize the GMRES(k) algorithm

      eps = etak               !GMRES convergence tolerance

c Call the preconditioned GMRES(k) algorithm

      call gmresMeth(nn,x0,b,x,vv,zz,eps,kmax,maxitgm,iout,ierr,itk)

c End program

      return
      end

c gmresMeth
c #######################################################################
      subroutine gmresMeth(ntot,x0,rhs,sol,vv,zz,eps,im,maxits
     .                    ,iout,ierr,its)
      implicit none             !For safe fortran
c----------------------------------------------------------------------*
c                                                                      *
c                     *** Preconditioned GMRES ***                     *
c                                                                      *
c----------------------------------------------------------------------*
c This is a simple version of the ILUT preconditioned GMRES algorithm. *
c The ILUT preconditioner uses a dual strategy for dropping elements   *
c instead  of the usual level of-fill-in approach. See details in ILUT *
c subroutine documentation. PGMRES uses the L and U matrices generated *
c from the subroutine ILUT to precondition the GMRES algorithm.        *
c The preconditioning is applied to the right. The stopping criterion  *
c utilized is based simply on reducing the residual norm by epsilon.   *
c This preconditioning is more reliable than ilu0 but requires more    *
c storage. It seems to be much less prone to difficulties related to   *
c strong nonsymmetries in the matrix. We recommend using a nonzero tol *
c (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
c lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
c more reliable the code is. Efficiency may also be much improved.     *
c Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
c Gaussian elimination without pivoting.                               *
c                                                                      *
c ILU(0) and MILU(0) are also provided for comparison purposes         *
c USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
c then call pgmres.                                                    *
c----------------------------------------------------------------------*
c Coded by Y. Saad - This version dated May, 7, 1990.                  *
c----------------------------------------------------------------------*
c parameters                                                           *
c-----------                                                           *
c on entry:                                                            *
c==========                                                            *
c                                                                      *
c n     == integer. The dimension of the matrix.                       *
c im    == size of krylov subspace.                                    *
c x0    == current newton state to calculate Jacobian                  *
c rhs   == real vector of length n containing the right hand side.     *
c          Destroyed on return.                                        *
c sol   == real vector of length n containing an initial guess to the  *
c          solution on input. approximate solution on output           *
c eps   == tolerance for stopping criterion. process is stopped        *
c          as soon as ( ||.|| is the euclidean norm):                  *
c          || current residual||/||initial residual|| <= eps           *
c maxits== maximum number of iterations allowed                        *
c iout  == output unit number number for printing intermediate results *
c          if (iout .le. 0) nothing is printed out.                    *
c                                                                      *
c aa, ja,                                                              *
c ia    == the input matrix in compressed sparse row format:           *
c          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
c          ja(1:nnz) = corresponding column indices.                   *
c          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
c                                                                      *
c alu,jlu== A matrix stored in Modified Sparse Row format containing   *
c           the L and U factors, as computed by subroutine ilut.       *
c                                                                      *
c ju     == integer array of length n containing the pointers to       *
c           the beginning of each row of U in alu, jlu as computed     *
c           by subroutine ILUT.                                        *
c                                                                      *
c on return:                                                           *
c==========                                                            *
c sol   == contains an approximate solution (upon successful return).  *
c ierr  == integer. Error message with the following meaning.          *
c          ierr = 0 --> successful return.                             *
c          ierr = 1 --> convergence not achieved in itmax iterations.  *
c          ierr =-1 --> the initial guess seems to be the exact        *
c                       solution (initial residual computed was zero)  *
c                                                                      *
c----------------------------------------------------------------------*
c                                                                      *
c work arrays:                                                         *
c=============                                                         *
c vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c          basis)                                                      *
c----------------------------------------------------------------------*
c arnoldi size should not exceed im=50 in this version.                *
c----------------------------------------------------------------------*
c----------------------------------------------------------------------*
c subroutines called :                                                 *
c amux   : SPARSKIT routine to do the matrix by vector multiplication  *
c          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
c lusol0 : combined forward and backward solves (Preconditioning ope.) *
c BLAS1  routines.                                                     *
c----------------------------------------------------------------------*

c Call variables

      integer*4     ntot,im,maxits,iout,ierr,its
      real*8        x0(ntot),rhs(ntot),sol(ntot),eps

c Local variables

      real*8        hh(im+1,im), c(im), s(im), rs(im+1)
      real*8        vv(ntot,im+1),zz(ntot,im+1)
      real*8        epsmac,rold,ro,eps1,gam,t
      integer*4     i,j,i1,k,k1,ii,jj,n,rstrt,irstrt,precout

      data epsmac /1.d-16/

c Begin program

      precout = iout - 2

      n = ntot
      its = 0

c Calculate magnitude of nonlinear residual

      rold = sqrt(sum(rhs*rhs))

      eps1=eps*rold

      if (rold .lt. (n*1d-16)) then
        ierr = -1
        return
      endif

c Compute initial residual vector

      call matrixFreeMatVec(n,x0,sol,vv)

      vv(:,1) = rhs(:) - vv(:,1)

c Calculate number of restarting loops

      rstrt = maxits/im + 1

c Restarted GMRES loop

      do irstrt = 1,rstrt

        ro = sqrt(sum(vv(:,1)*vv(:,1)))

        if (iout .gt. 0 .and. its .eq. 0) then
          write(*, 199) its, ro/rold
        endif

        if (ro.le.eps1) exit

        t = 1.0d0/ro
        vv(:,1) = vv(:,1)*t

cUGH        if (its .eq. 0) eps1=eps*ro

c     Initialize 1-st term  of rhs of hessenberg system

        rs(1) = ro

c     GMRES iteration

        do i = 1,im

          its = its + 1
          i1  = i + 1

cFG           call applyPreconditioner(n,vv(1,i),rhs,precout)
          call applyPreconditioner(n,vv(1,i),zz(1,i),precout)

cFG           call matrixFreeMatVec(n,x0,rhs,vv(1,i1))
          call matrixFreeMatVec(n,x0,zz(1,i),vv(1,i1))

c       Modified gram - schmidt.

          do j=1,i
            t = sum(vv(:,j)*vv(:,i1))
            hh(j,i) = t
            vv(:,i1) = vv(:,i1)-t*vv(:,j)
          enddo

          t = sqrt(sum(vv(:,i1)*vv(:,i1)))
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

          if (gam .eq. 0.0d0) gam = epsmac

c       Get next plane rotation

          c(i)   = hh(i ,i)/gam
          s(i)   = hh(i1,i)/gam
          rs(i1) = -s(i)*rs(i)
          rs(i)  =  c(i)*rs(i)

c       Determine residual norm and test for convergence

          hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
          ro = abs(rs(i1))

          if (iout .gt. 0) then
            write(*, 199) its, ro/rold
          endif

          if (ro.le.eps1.or.i.eq.im) exit

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

c      Form linear combination of v(*,i)'s to get solution

        t = rs(1)
        rhs(:) = zz(:,1)*t
cFG       rhs(:) = vv(:,1)*t

        do j=2, i
          t = rs(j)
          rhs(:) = rhs(:) + t*zz(:,j)
cFG         rhs(:) = rhs(:) + t*vv(:,j)
        enddo

c     Call preconditioner

cFG       call applyPreconditioner(n,rhs,vv(1,im+1),precout)

c     Update solution

        sol = sol + rhs
cFG       sol(:) = sol(:) + vv(:,im+1)

c      Check convergence and restart outer loop if necessary

        if (ro.le.eps1) then
          ierr = 0
          exit
        endif

        if (its.ge.maxits) then
          ierr = 1
          exit
        endif

c      Else compute residual vector and continue

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

c      Restart outer loop

       enddo

c End program

       return

 199   format(' GMRES its =', i4,';  res/rold norm =', 1pd10.2)

       end

c matrixFreeMatVec
c##################################################################
      subroutine matrixFreeMatVec(nn,x,z,y)

c------------------------------------------------------------------
c   This subroutine computes
c
c     y = J(x)*z
c
c   where J is the Jacobian matrix.  A finite difference
c   approximation is used to compute this product,
c
c     y ~ (F(x+ez) - F(x))/e
c
c   where e is some small perturbation constant
c
c   The integer array idiag contains the information
c   that relates the actual column to the main diagonal
c   since only the non-zero diagonals of A are stored, i.e.,
c
c     i = row number
c     j = column number = i + idiag(k)
c     k = non-zero diagonal index numbered from left to right
c------------------------------------------------------------------

      use newton_gmres

      implicit none       !For safe fortran

c Call variables

      integer*4   nn
      real*8      x(nn),z(nn),y(nn)

c Local variables

      real*8      dummy(nn)

      real*8      pert,pertinv,scale,eps
      integer*4   i,j,m,ii

c Begin program

      pert  = 1d0
      scale = 0d0

      scale = sqrt(sum(z*z))
      pert  = sqrt(pert + sum(x*x))

c Calculate J.x

      if (scale .lt. 1.0d-16) then

         do i = 1,nn
            y(i) = 0d0
         enddo

      else

        eps = 1d-6

cc        pert = eps*dsqrt(pert/scale)/nn
        pert = eps*(1.+pert/scale/nn)
cc        pert = eps*(1.+dsqrt(pert/scale))
cc        pert = eps*2d0/dsqrt(scale)

c     Perturb state variables x + eps z --> dummy

        dummy = x + pert*z

c     Nonlinear function evaluation --> y

        call evaluateNonlinearResidual(nn,dummy,y)

c     Compute the product J.x using the matrix-free approx

        y = (y-res)/pert

      endif

c End

      return
      end

c fmedval
c####################################################################
      real*8 function fmedval(p1,p2,p3)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c    This function computes intermediate value of p1, p2, p3.
c--------------------------------------------------------------------

c Call variables

      real*8       p1,p2,p3

c Local variables

c Begin program

      fmedval = min( max(p1,p2) , max( p3,min(p1,p2) ) )

c End

      return
      end

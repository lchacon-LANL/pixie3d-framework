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
c      Deallocates dynamic storage space used in preconditioner
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


c module newtonGmres
c ###################################################################
      module newton_gmres

      double precision, allocatable, dimension(:):: res,x0,xk

      double precision :: dt,check

      integer          :: jit

      end module newton_gmres

c newtonGmres
c####################################################################
      subroutine newtonGmres(neq,ntot,x,etak_meth,damp,global,dt0
     .           ,tolgm,kmax,gmit_max,tolnewt,ntit_max_acc
     .           ,ntit_max_rej,gmit_out,ntit_out,iguess,out,ierr)
c--------------------------------------------------------------------
c     Performs Jacobian-free inexact Newton iteration on b(x) = 0, 
c     where b is calculated in 'evaluateNonlinearResidual' (provided
c     by user). 
c
c     Call parameters:
c       * neq: number of equations
c       * ntot: dimension of vector of unknowns.
c       * x: on input, initial guess; on output, solution.
c       * etak_meth: specifies method to determine inexact Newton forcing 
c           parameter. Currently:
c               etak_meth = 0 --> constant
c               etak_meth > 0 --> power law adaptive strategy
c       * damp  : damping parameter for Newton update
c       * global: determines the globalization procedure:
c               global = 0 --> no globalization
c               global = 1 --> linesearch backtracking 
c               global = 2 --> pseudo-transient
c       * dt0: initial time step for pseudo-transient
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

      integer          :: ntot,neq,etak_meth,kmax,gmit_max,ntit_max_acc
     .                   ,ntit_max_rej,gmit_out,ntit_out,out,iguess,ierr

      double precision :: x(ntot),tolgm,tolnewt,damp

      integer          :: global

c Local variables

      integer*4   itk,i,j,ig,iout

      real*8      dxavg,check_lim,residuals(neq)

      real*8      b(ntot),ddx(ntot),dt0,roundoff

      real*8      f0,fkm,fk,fkp,fm,flimit,etak,etakm,theta,dampm

c Begin program

      if (ntit_max_acc > ntit_max_rej) then
        write (*,*) 'Error in Newton input: ntit_max_acc > ntit_max_rej'
        write (*,*) 'Aborting...'
        stop
      endif

c Initialize

      allocate (res(ntot),x0(ntot),xk(ntot))

      ierr     = 0
      gmit_out = 0
      itk      = 0

      if (out.ge.1) write (*,200)

      check_lim = 1d30
cc      if (global.ne.2) check_lim = 1d1

      if (global.ne.2) dt0 = 1d30
      dt = dt0

      x0 = x    !Save initial guess
      xk = x    !Save previous Newton state

      roundoff = ntot*1d-15

c Evaluate rhs and norms

      call evaluateNewtonResidual(ntot,x,res)

      b = -res

      f0 = sqrt(sum(b*b))

      !Check if initial guess is exact, and if so exit Newton step
      if (f0.lt.roundoff) then
        ierr = -1
        return
      endif

c Initial output

      if (out.ge.1) write (*,210) 0,0d0,f0,1d0,1d0,0

      if (out.ge.3) then
        call calculateResidualNorms(neq,ntot,b,residuals)
        do i = 1,neq
          write (*,320) i,residuals(i)
        enddo
      endif

c Start Newton iteration

      fk    = f0
      fkm   = f0
      flimit= roundoff + tolnewt*f0

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

c     Check for error in GMRES

        if (ierr.eq.1) then
          write(*,*) 
     .          '   Exceeded maximum GMRES iterations (',gmit_max,')'
          ierr = 0  !Use GMRES solution for Newton update regardless
        elseif (ierr.eq.-1) then !Got to roundoff in the Newton residual
          ierr = 0
          exit !Do not continue Newton iteration
        endif

c     Globalization procedure

c       Determine inexact Newton constant

        call find_etak (fk,fkm,flimit,tolgm,etak,etakm,etak_meth)

c       Calculate damping coefficient

cc        if (global.eq.1) call findDamping (ntot,x,ddx,etak,fk,damp)
        if (global.ge.1) call findDamping (ntot,x,ddx,etak,fk,damp)

c     Update solution (x = x + ddx)

        xk = x   !Save previous Newton state

        call updateNewtonSolution(x,ddx,ntot,damp,dxavg)

c     Evaluate rhs and norms

        call evaluateNewtonResidual(ntot,x,res)

        b = -res

        fkp = sqrt(sum(b*b))

c     Check Newton convergence

        check = fkp/f0

        if (out.ge.1) write (*,210) jit,dxavg,fkp,check,damp,dt,itk

        if (out.ge.3) then
          call calculateResidualNorms(neq,ntot,b,residuals)
          do i = 1,neq
            write (*,320) i,residuals(i)
          enddo
        endif

        if(fkp.lt.flimit.or.check.gt.check_lim) exit

c     Store magnitude of residual

        fkm   = fk
        fk    = fkp

        etakm = etak

c     Change pseudo-transient time step

cc        if (global.eq.2) dt = min(dt/check,1d30)
        if (global.eq.2) dt = min(dt0/check,1d30)
cc        if (global.eq.2.and.check.lt.1d-3) dt = min(dt0*jit,1d0)
cc        if (global.eq.2) dt = min(dt0/check**2,1d30)

      enddo       !End of Newton loop

c Check error in Newton convergence

      !No convergence: reject solution
      if (jit.gt.ntit_max_rej.or.check.gt.check_lim) then 
        if (jit.gt.ntit_max_rej) write (*,220) check
        if (check.gt.check_lim) write (*,230) check,check_lim
        if (check.gt.check_lim) write (*,*) check,check_lim
        ierr = 1
      !Accept solution, but warn user
      elseif ((jit-1).gt.ntit_max_acc) then 
        write (*,240) ntit_max_acc
        ierr = 2
      endif

c End program

      deallocate (res,x0,xk)

      return

 200  format 
     .   (/,' New_it   Av_updt    Abs_res   Rel_res     Damping'
     .     ,'      dt_n     GMRES')
 210  format (i4,3x,1p5e11.3,3x,i4)
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
      subroutine find_etak (fk,fkm,flimit,tolgm,etak,etakm,etak_meth)
      implicit none       !For safe fortran
c--------------------------------------------------------------------
c     Finds inexact Newton forcing parameter, using two methods:
c     constant (etak_meth = 0) or power law adaptive strategy.
c--------------------------------------------------------------------

c Call variables

      integer*4   etak_meth
      real*8      fk,fkm,flimit,tolgm,etak,etakm

c Local variables

      real*8      gamm,alph

c Begin program

      if (etak_meth.eq.0) then
        etak = tolgm
      else
        gamm = .9
        alph = 1.5

        !For superlinear convergence, alph>1
        etak = gamm*(fk/fkm)**alph

        !First safeguard: avoid sharp decrease of etak
        etak = min(tolgm,max(etak,gamm*etakm**alph))

        !Second safeguard: avoid oversolving
        etak = min(tolgm,max(etak,gamm*flimit/fk))
      endif

cc      write (*,*) etak

c End program

      return
      end

c findDamping
c####################################################################
      subroutine findDamping (ntot,x,ddx,etak,fk,damp)
      implicit none       !For safe fortran
c--------------------------------------------------------------------
c     If global is true, uses linesearch backtracking to ensure
c     sufficient reduction in the Newton residual norm.
c--------------------------------------------------------------------

c Call variables

      integer ::  ntot
      double precision ::  x(ntot),ddx(ntot),fk,etak,damp

c Local variables

      integer ::  idamp
      double precision :: dampm,theta,etak0,dxavg,fkp,fm
      double precision :: b(ntot),dummy(ntot)

c External

      double precision :: fmedval
      external            fmedval

c Begin program

      damp = 1d0
      etak0 = etak

      do idamp = 1,100
        dummy = x
        call updateNewtonSolution(dummy,ddx,ntot,damp,dxavg)
        call evaluateNewtonResidual(ntot,dummy,b)
        fkp = sqrt(sum(b*b))
        if (fkp.gt.(1.-1e-4*(1.-etak))*fk) then
          dampm = .5*damp
          dummy = x
          call updateNewtonSolution(dummy,ddx,ntot,dampm,dxavg)
          call evaluateNewtonResidual(ntot,dummy,b)
          fm = sqrt(sum(b*b))
          theta = .5*(1.5*fk + 0.5*fkp - 2.*fm)/(fk + fkp - 2.*fm)
          theta = fmedval(1d-1,8d-1,theta)
          damp = damp*theta
cc          write (*,*) 'New damping parameter',damp
          etak = 1.- theta*(1.-etak)
        else
          exit
        endif
      enddo

cc      if (damp.lt.1d-5) etak = etak0

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
c   iteration algorithm with right pre-conditioning
c   Note that this subroutine calls the gmres subroutine from
c   the SPARSKIT package by Y. Saad, modified by L. Chacon.
c --------------------------------------------------------------------

      implicit none    !For safe fortran

c Call variables

      integer*4   itk,ierr,nn,kmax,iguess,iout,maxitgm
      real*8      b(nn),x(nn),x0(nn),etak

c Local variables

      integer*4   iiout
      real*8      eps

c Begin program

c Compute the initial guess x

      if(iguess.eq.1) then
        iiout = iout - 2
        call applyPreconditioner(nn,b,x,iiout)
      else
        x = 0d0
      endif

c Initialize the GMRES(k) algorithm

      eps = etak               !GMRES convergence tolerance

c Call the preconditioned GMRES(k) algorithm

      call fgmresMeth(nn,x0,b,x,eps,kmax,maxitgm,iout,ierr,itk)

c End program

      return
      end

c gmresMeth
c ######################################################################
      subroutine fgmresMeth(ntot,x0,rhs,sol,eps,im,maxits,iout,ierr,its)
      implicit none             !For safe fortran
c----------------------------------------------------------------------*
c                                                                      *
c               *** Preconditioned Flexible GMRES ***                  *
c                                                                      *
c----------------------------------------------------------------------*
c This is a simple version of the right-preconditioned GMRES algorithm.*
c The stopping criterion utilized is based simply on reducing the      *
c residual norm by epsilon.                                            *
c----------------------------------------------------------------------*
c parameters                                                           *
c-----------                                                           *
c on entry:                                                            *
c==========                                                            *
c                                                                      *
c ntot  == integer. The dimension of the matrix.                       *
c x0    == current newton state to calculate Jacobian                  *
c rhs   == real vector of length n containing the right hand side.     *
c          Destroyed on return.                                        *
c sol   == real vector of length n containing an initial guess to the  *
c          solution on input. approximate solution on output           *
c eps   == tolerance for stopping criterion. process is stopped        *
c          as soon as ( ||.|| is the euclidean norm):                  *
c          || current residual||/||initial residual|| <= eps           *
c im    == size of krylov subspace.                                    *
c maxits== maximum number of GMRES iterations allowed                  *
c iout  == output unit number number for printing intermediate results *
c          if (iout .le. 0) nothing is printed out.                    *
c                                                                      *
c on return:                                                           *
c==========                                                            *
c sol   == contains an approximate solution (upon successful return).  *
c ierr  == integer. Error message with the following meaning.          *
c          ierr = 0 --> successful return.                             *
c          ierr = 1 --> convergence not achieved in itmax iterations.  *
c          ierr =-1 --> the initial guess seems to be the exact        *
c                       solution (initial residual computed was zero)  *
c its   == final number of GMRES iterations                            *
c                                                                      *
c----------------------------------------------------------------------*
c                                                                      *
c work arrays:                                                         *
c=============                                                         *
c vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c          basis)                                                      *
c zz    == work array of length  n x (im+1) (used to store the Arnoli  *
c          basis times the preconditioner operator (FGMRES))
c----------------------------------------------------------------------*
c arnoldi size should not exceed im=50 in this version.                *
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

      if (rold .lt. (n*1d-15)) then
        ierr = -1
        return
      endif

c Compute initial residual vector

      call matrixFreeMatVec(n,x0,sol,vv)

      vv(:,1) = rhs(:) - vv(:,1)

c Calculate number of restarting loops

      if (maxits > 0) then
        rstrt = maxits/im + 1
      else
        rstrt = 0
      endif

c Restarted GMRES loop

      do irstrt = 1,rstrt

        ro = sqrt(sum(vv(:,1)*vv(:,1)))

        if (iout .gt. 0 ) then
          if (its .eq. 0) then
            write (*,199) its, ro/rold
          else
            write (*,*) 'Restarting GMRES... (',irstrt-1,')'
          endif
        endif

        if (ro.le.eps1) exit

        t = 1.0d0/ro
        vv(:,1) = vv(:,1)*t

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

c      Form linear combination of zz=P*vv to get solution

        t = rs(1)
        rhs(:) = zz(:,1)*t
cFG       rhs(:) = vv(:,1)*t

        do j=2, i
          t = rs(j)
          rhs(:) = rhs(:) + t*zz(:,j)
cFG         rhs(:) = rhs(:) + t*vv(:,j)
        enddo

c     Call preconditioner

        vv(:,im+1) = rhs(:)
cFG       call applyPreconditioner(n,rhs,vv(1,im+1),precout)

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

c End program

      return

 199  format(' GMRES its =', i4,';  res/rold norm =', 1pd10.2)

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

      integer*4   i
      real*8      dummy(nn),pert,scale,eps

c Begin program

      eps   = 1d-6

c Calculate difference parameter

      scale = sum(z*z)

cc      pert  = (1d0 + sum(x*x))/scale
cc      pert  = eps*sqrt(pert)

      pert  = sum(z*x)
      pert  = eps*(sqrt(scale)+abs(pert))/scale*sign(1d0,pert)

cc      write (*,*) pert,scale

c Calculate J.x

      if (scale .lt. 1.0d-16) then

         do i = 1,nn
            y(i) = 0d0
         enddo

      else

c     Perturb state variables x + eps z --> dummy

        dummy = x + pert*z

c     Nonlinear function evaluation --> y

        call evaluateNewtonResidual(nn,dummy,y)

c     Compute the product J.x using the matrix-free approx

        y = (y-res)/pert

      endif

c End

      end subroutine matrixFreeMatVec

c evaluateNewtonResidual
c####################################################################
      subroutine evaluateNewtonResidual(ntot,x,f)

c--------------------------------------------------------------------
c     Calculates nonlinear residuals. 
c--------------------------------------------------------------------

      use newton_gmres

      implicit none

c Call variables

      integer  :: ntot
      real*8      x(ntot),f(ntot)

c Local variables

c Begin program

c Evaluate nonlinear residual

      call evaluateNonlinearResidual(ntot,x,f)

c Add pseudo-transient term

      f = (x - x0)/dt - f

c End program

      end subroutine

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

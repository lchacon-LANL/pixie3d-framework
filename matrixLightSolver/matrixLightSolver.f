c TO DO list:
c
c 1) Implement boundary conditions information for restriction, prolongation
c    and matrix element determination.
c
c 2) Posibility of linking to user-provided solvers.
c
c 3) Separating 2D dependent routines, to allow more general use in 3D.


c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c THIS ROUTINE WORKS IN COMBINATION WITH THE FOLLOWING ROUTINES:
c     coupledMG.f, mlsolver_mod.f mg_mod.f
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c getSolver
c####################################################################
      recursive subroutine getSolver(ntot,b,x,matvec,igrid,guess
     .                              ,out,depth)
c--------------------------------------------------------------------
c     Solves Ax=b matrix-free using a variety of solvers.
c--------------------------------------------------------------------

      use mlsolverSetup

      implicit none     !For safe fortran

c Call variables

      integer*4      ntot,igrid,iter,guess,out,depth
      real*8         x(ntot),b(ntot)

      external       matvec

c Local variables

      type (solver_unit)   :: solver_def

      type (solver_options):: options
      character*(2)  solver

c Begin program

c Read solver definition from solver hierarchy

      call readSolverHierarchy (solver_def,depth)

      solver  = solver_def%solver
      options = solver_def%options

c Symmetry test

      if (options%sym_test) call symm_test(2,matvec)

c Invoke solver

      call solver_meth(ntot,b,x,matvec,igrid
     .                ,guess,out,solver,options,depth)

c Check convergence

      if (options%tol_out.gt.options%tol.and.out.ge.2) then
        write (*,100) solver,options%tol_out,options%tol
c        stop
      endif

c Output information

      solver_def%options%tol_out  = options%tol_out
      solver_def%options%iter_out = options%iter_out

      call writeSolverHierarchy(solver_def,depth)

cc      write (*,*) options%iter_out

c End program

 100  format (/,' ',a2,' did not converge to prescribed tolerance',
     .        /,' Residual =',1pe10.2,' > Tolerance =',e10.2,/)

      end subroutine

c solver_meth
c####################################################################
      recursive subroutine solver_meth(ntot,b,x,matvec,igrid
     .                               ,guess,out,solver,options,depth)
c--------------------------------------------------------------------
c     Solves Ax=b matrix-free using a variety of solvers.
c--------------------------------------------------------------------

      use mlsolverSetup

      implicit none     !For safe fortran

c Call variables

      integer*4      ntot,igrid,iter,guess,out,depth
      real*8         x(ntot),b(ntot)

      type (solver_options):: options
      character*(2)  solver

      external       matvec

c Local variables

c Begin program

c Solve

      select case (solver)
      case ('cg')

        call cgmeth(ntot,b,x,matvec,options,igrid,guess,out,depth)

      case ('gm')

        call gmres(ntot,b,x,matvec,options,igrid,guess,out,depth)

      case ('jb')

        call jb(ntot,b,x,matvec,options,igrid,guess,out,depth)

      case ('gs')

        call gs(ntot,b,x,matvec,options,igrid,guess,out,depth)

      case ('mg')

        call mg(ntot,b,x,matvec,options,igrid,guess,out,depth)

      case ('id')

        x = b

      case default

        write (*,*) 'Solver specified (',solver,') not valid.'
        write (*,*) 'Aborting...'
        stop

      end select

c End program

      end subroutine

c cgmeth
c####################################################################
      recursive subroutine cgmeth(ntot,b,x,matvec,options,igrid
     .                           ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free Preconditioned Conjugate Gradient routine to solve
c     Ax = b. Call variables:
c       * ntot: grid dimension
c       * b,x: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0->no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      implicit none     !For safe fortran

c Call variables

      integer*4    ntot,igrid,guess,out,depth
      real*8       x(ntot),b(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer*4    niter,stp_test,depth1
      integer*4    i,iter,smit,nn,vcycle,precout
      real*8       alpha,beta,mag,mag1,rr0,rho0,rho1
      real*8       smtol,tol,abstol

      double precision:: rr(ntot),zz(ntot),pp(ntot),qq(ntot)

c Begin program

      nn = ntot

      precout = max(0,out-2)

      depth1 = depth + 1

      abstol = 1d-16*nn

c Extract options

      niter    = options%iter      !Maximum number of iterations
      stp_test = options%stp_test  !Stopping test type (0 -> rhs, 1-> residual)
      tol      = options%tol       !Convergence tolerance

c Compute initial residual

      if (guess.eq.0) then

c     For zero initial guess

        x  = 0d0
        rr = b

      else

c     For arbitrary initial guess (rr = b - Ax)

        call matvec(0,ntot,x,rr,igrid)

        rr = b - rr

      endif

c CG solve

      if (stp_test.eq.1) then
        rr0 = sqrt(sum(rr*rr))
      else
        rr0 = sqrt(sum(b*b))
      endif

      if (rr0.lt.abstol) then
        if (out.ge.1) 
     .        write (*,*) 'Initial guess seems exact solution in CG' 
        return
      endif

cc      iter = 0
cc      if (out.ge.2) write (*,10) iter,rr0,rr0/rr0

c     Preconditioner (zz = P rr)

      call getSolver(ntot,rr,zz,matvec,igrid,0,precout,depth1)

c     CG preparation

      rho1 = sum(rr*zz)

      pp = zz

c     CG step

      do iter= 1,niter

        rho0 = rho1

        call matvec(0,ntot,pp,qq,igrid)

        alpha = rho0/sum(pp*qq)

        x  = x  + alpha*pp

        rr = rr - alpha*qq

c     Check convergence

        mag = sqrt(sum(rr*rr))
        mag1 = mag/rr0
        if (out.ge.2) write (*,10) iter,mag,mag1
cc        if (mag1.lt.tol) exit
        if (mag.lt.(abstol + tol*rr0)) exit

c     Preconditioner (zz = P rr)

        call getSolver(ntot,rr,zz,matvec,igrid,0,precout,depth1)

c     Prepare next step

        rho1 = sum(rr*zz)

        beta = rho1/rho0

        pp = beta*pp + zz

      enddo

      if (out.eq.1) write (*,10) iter,mag,mag1
      if (out.ge.2) write (*,*)

      options%iter_out = iter 
      options%tol_out  = mag1

c End program

      return

 10   format (' CG Iteration:',i4,'; Residual:',1p1e12.4,
     .        '; Ratio:',1p1e12.4)

      end subroutine

c gmres
c#######################################################################
      recursive subroutine gmres(ntot,rhs,sol,matvec,options,igrid
     .                          ,guess,iout,depth)
c--------------------------------------------------------------------
c     Matrix-free Preconditioned Conjugate Gradient routine to solve
c     Ax = b. Call variables:
c       * ntot: grid dimension
c       * rhs,sol: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0->no initial guess; 1 -> initial guess provided.
c       * iout: convergence info output on screen if iout > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      implicit none       !For safe fortran
      
c Call variables

      integer*4     ntot,im,iout,ierr,igrid,its,guess,depth
      real*8        rhs(ntot),sol(ntot)

      external      matvec

      type (solver_options) :: options

c Local variables

      integer*4     kmax,stp_test,depth1
      real*8        eps,epsmac,rold,ro,eps1,gam,tt,mag,abstol

      integer*4     i,j,i1,k,k1,ii,jj
      integer*4     nn,rstrt,irstrt,precout,maxits

      double precision,allocatable,dimension(:)  :: c,s,rs
      double precision,allocatable,dimension(:,:):: hh,vv,zz

      data epsmac /1.d-16/

c Begin program

      its = 0
      nn = ntot

      precout = max(0,iout-2)

      depth1 = depth + 1

      abstol = 1d-16*nn

c Extract options

      maxits   = options%iter             !Maximum number of GMRES its
      stp_test = options%stp_test         !Stopping test type (0 -> rhs, 1-> residual)
      eps      = options%tol              !Convergence tolerance
      kmax     = options%krylov_subspace  !Maximum krylov subspace

c Allocate work arrays

      allocate(hh(kmax+1,kmax),vv(ntot,kmax+1),zz(ntot,kmax+1))
      allocate(c(kmax),s(kmax),rs(kmax+1))

c Compute initial residual vector

      if (guess.eq.0) then

c     For zero initial guess

        vv(:,1) = rhs(:)

      else

c     For arbitrary initial guess (vv(:,1) = b - Ax)

        call matvec(0,ntot,sol,vv,igrid)

        vv(:,1) = rhs(:) - vv(:,1)

      endif

c Calculate restarting loops

      rstrt = maxits/kmax + 1

c Calculate magnitude of initial residual

      if (stp_test.eq.1) then
        rold = sqrt(sum(vv(:,1)*vv(:,1)))
      else
        rold = sqrt(sum(rhs*rhs))
      endif

      eps1=eps*rold

c Restarted GMRES loop

      do irstrt = 1,rstrt

        ro = sqrt(sum(vv(:,1)*vv(:,1)))

cc        if (iout.ge.2.and.its.eq.0) write(*,10) its,ro,ro/rold

cc        if (its .eq. 0) eps1=eps*rold

        if (ro.lt.abstol) then
          ierr = -1
          if (iout.ge.1) 
     .       write (*,*) 'Initial guess seems exact solution in GMRES' 
          exit
        elseif (ro.le.(abstol+eps1)) then
          ierr = 0
          exit
        endif

        tt = 1.0d0/ro
        vv(:,1) = vv(:,1)*tt

c      Initialize 1-st term  of rhs of hessenberg system

        rs(1) = ro

c      GMRES iteration

        do i = 1,kmax

          its = its + 1
          i1  = i + 1

c        Call preconditioner

          call getSolver(ntot,vv(1,i),zz(1,i),matvec,igrid,0,precout
     .                  ,depth1)

          call matvec(0,ntot,zz(1,i),vv(1,i1),igrid)

c        Modified gram - schmidt.

          do j=1, i
            tt = sum(vv(:,j)*vv(:,i1))
            hh(j,i) = tt
            vv(:,i1) = vv(:,i1) - tt*vv(:,j)
          enddo

          tt = sqrt(sum(vv(:,i1)*vv(:,i1)))
          hh(i1,i) = tt

          if ( tt .ne. 0.0d0) then
            tt = 1.0d0/tt
            vv(:,i1) = vv(:,i1)*tt
          endif

c        Done with modified Gram-Schimdt and arnoldi step
c        Now update factorization of hh
c        Perform previous transformations on i-th column of h

          if (i .gt. 1) then
            do k=2,i
              k1 = k-1
              tt = hh(k1,i)
              hh(k1,i) =  c(k1)*tt + s(k1)*hh(k,i)
              hh(k ,i) = -s(k1)*tt + c(k1)*hh(k,i)
            enddo
          endif

          gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)

c        If gamma is zero then any small value will do
c        Will affect only residual estimate

          if (gam .eq. 0.0d0) gam = epsmac

c        Get next plane rotation

          c(i)   = hh(i ,i)/gam
          s(i)   = hh(i1,i)/gam
          rs(i1) = -s(i)*rs(i)
          rs(i)  =  c(i)*rs(i)

c        Determine residual norm and test for convergence

          hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
          ro = abs(rs(i1))

          if (iout.ge.2) write(*,10) its,ro,ro/rold

          if (ro.le.(abstol+eps1).or.i.ge.min(kmax,maxits)) exit

        enddo

c      Now compute solution. First solve upper triangular system

        rs(i) = rs(i)/hh(i,i)
        do ii=2,i
          k=i-ii+1
          k1 = k+1
          tt=rs(k)
          do j=k1,i
            tt = tt-hh(k,j)*rs(j)
          enddo
          rs(k) = tt/hh(k,k)
        enddo

c      Form linear combination of z(*,i)'s to get solution

        tt = rs(1)
        do k=1, nn
          rhs(k) = zz(k,1)*tt
        enddo

        do j=2, i
          tt = rs(j)
          do k=1, nn
            rhs(k) = rhs(k)+tt*zz(k,j)
          enddo
        enddo

c      Form solution

        sol = sol + rhs

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
          tt = rs(j)
          if (j .eq. 1)  tt = tt - 1d0
          vv(:,1) = vv(:,1) + tt*vv(:,j)
        enddo

c      Restart outer loop

      enddo

      if (iout.eq.1) write (*,10) its,ro,ro/rold
      if (iout.ge.2) write (*,*)

      options%iter_out = its
      options%tol_out  = ro/rold

c End program

      deallocate(hh,vv,zz)
      deallocate(c,s,rs)

      return

 10   format (' GMRES Iteration:',i4,'; Residual:',1p1e12.4,
     .        '; Ratio:',1p1e12.4)

      end subroutine

c symm_test
c####################################################################
      recursive subroutine symm_test(igrid,matvec)
c--------------------------------------------------------------------
c     Performs symmetry test on matvec.
c--------------------------------------------------------------------

      use mg_setup

      implicit none     !For safe fortran

c Call variables

      integer*4      igrid

      external       matvec

c Local variables

      double precision,allocatable,dimension(:)::x1,dummy,dummy2

      real*8         dd1,dd2,error
      integer*4      nx,ny,i,j,ii,nn,jj,i1,j1,i2,j2,ix1,iy1,ix2,iy2

c Begin program

      nx = nxvp(igrid)
      ny = nyvp(igrid)

      nn = nx*ny

      write (*,*) 'Performing symmetry test of system matrix ',
     .            'on grid:',nx,'x',ny,'...'

      allocate(x1(nn),dummy(nn),dummy2(nn))

      error = 0d0

      do ii = 1,nn
        x1    (ii) = 0d0
        dummy (ii) = 0d0
        dummy2(ii) = 0d0
      enddo

      do ii = 1,nn

c     Find column vector ii

        call findBaseVector(ii,1,1,nn,igrid,x1,1d0)

        call matvec(0,nn,x1,dummy,igrid)

        call findBaseVector(ii,1,1,nn,igrid,x1,0d0)

c     Compare column vector ii with corresponding row vector (intersect in
c     diagonal)

        do jj = ii,nn

          call findBaseVector(jj,1,1,nn,igrid,x1,1d0)

          call matvec(ii,nn,x1,dummy2,igrid)

          call findBaseVector(jj,1,1,nn,igrid,x1,0d0)

          dd1 = abs(dummy(jj) - dummy2(ii))
          if (abs(dummy(jj)).gt.1d-15.or.abs(dummy2(ii)).gt.1d-15) then
            write(*,15) jj,ii,dummy(jj),ii,jj,dummy2(ii),dd1
     .                ,100*dd1/max(abs(dummy(jj)),abs(dummy2(ii)))
            error = error + dd1
          endif

        enddo
      enddo

      write (*,20) error

      stop

c End program

      deallocate (x1,dummy,dummy2)

      return

 15   format ('(',i3,',',i3,'):',1pe10.2,'; (',i3,',',i3,'):',e10.2,
     .        '  Error:',e10.2,'  %error:',0pf7.2)
 20   format (/,'Total relative error:',1pe10.3)

      end subroutine

c limits
c####################################################################
      subroutine limits(elem,nx,ny,imin,imax,jmin,jmax)
      implicit none
c--------------------------------------------------------------------
c     Finds limits on loops for matvec routines. Used in finding 
c     diagonal from matvec.
c--------------------------------------------------------------------

c Call variables

      integer*4       elem,nx,ny,imin,imax,jmin,jmax

c Local variables

c Begin program

      if (elem.eq.0) then
        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
      else
        imin = mod(elem,nx)
        if (imin.eq.0) imin = nx
        imax = imin
        jmin = 1 + (elem - imin)/nx
        jmax = jmin
      endif

c End program

      return
      end subroutine

c TO DO list:
c
c 1) Posibility of linking to user-provided solvers.

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c THIS ROUTINE WORKS IN COMBINATION WITH THE FOLLOWING ROUTINES:
c     mlsolver_mod.f,coupledMG.f
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c THIS ROUTINE REQUIRES THE FOLLOWING EXTERNAL ROUTINES:
c
c     subroutine matvec(elem,neq,ntot,x,b,igrid,bcnd)
c
c WHERE:
c     * elem: whether matvec is to be applied to the whole vector
c             (elem=0), or to find only the component elem<>0 of
c             the vector.
c     * neq: number of coupled equations
c     * ntot: size of vector
c     * x(ntot): vector to apply matrix operator on.
c     * b(ntot): resulting vector (b=Ax)
c     * igrid (integer): grid level operation is applied at.
c     * bcnd(6,neq): boundary condition information for all
c             dimensions of the problem.
c    
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c getSolver
c####################################################################
      recursive subroutine getSolver(neq,ntot,b,x,matvec,igrid,bcnd
     .                              ,guess,out,depth)
c--------------------------------------------------------------------
c     Solves Ax=b matrix-free using a variety of solvers.
c--------------------------------------------------------------------

      use mlsolverSetup

      implicit none     !For safe fortran

c Call variables

      integer(4) ::  neq,ntot,igrid,iter,guess,out,depth,bcnd(6,neq)
      real(8)    ::  x(ntot),b(ntot)

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

      if (options%sym_test)
     .     call symm_test(neq,options%ngrd_tst,matvec,bcnd)

c Invoke solver

      call solver_meth(neq,ntot,b,x,matvec,igrid,bcnd
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
     .        /,' Relat. residual =',1pe10.2,' > Tolerance =',e10.2,/)

      end subroutine getSolver

c solver_meth
c####################################################################
      recursive subroutine solver_meth(neq,ntot,b,x,matvec,igrid,bcnd
     .                                ,guess,out,solver,options,depth)
c--------------------------------------------------------------------
c     Solves Ax=b matrix-free using a variety of solvers.
c--------------------------------------------------------------------

      use mlsolverSetup

      implicit none     !For safe fortran

c Call variables

      integer(4) ::  neq,ntot,igrid,iter,guess,out,depth,bcnd(6,neq)
      real(8)    ::  x(ntot),b(ntot)

      type (solver_options):: options
      character*(2)  solver

      external       matvec

c Local variables

c Begin program

c Solve

      select case (solver)
      case ('cg')

        call cg(neq,ntot,b,x,matvec,options,igrid,bcnd,guess,out,depth)

      case ('gm')

        call gm(neq,ntot,b,x,matvec,options,igrid,bcnd,guess,out,depth)

      case ('jb')

        call jb(neq,ntot,b,x,matvec,options,igrid,bcnd,guess,out,depth)

      case ('gs')

        call gs(neq,ntot,b,x,matvec,options,igrid,bcnd,guess,out,depth)

      case ('mg')

        call mg(neq,ntot,b,x,matvec,options,igrid,bcnd,guess,out,depth)

      case ('id')

        x = b

      case default

        write (*,*) 'Solver specified (',solver,') not valid.'
        write (*,*) 'Aborting...'
        stop

      end select

c End program

      end subroutine solver_meth

c cg
c####################################################################
      recursive subroutine cg(neq,ntot,b,x,matvec,options,igrid,bcnd
     .                       ,guess,out,depth)
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

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: x(ntot),b(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: niter,stp_test,depth1,ierr
      integer(4) :: i,iter,smit,nn,vcycle,precout
      real(8)    :: alpha,beta,mag,mag1,rr0,rho0,rho1
      real(8)    :: smtol,tol,abstol
      real(8)    :: rr(ntot),zz(ntot),pp(ntot),qq(ntot)

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

        call matvec(0,neq,ntot,x,rr,igrid,bcnd)

        rr = b - rr

      endif

c CG solve

      if (stp_test.eq.1) then
        rr0 = sqrt(sum(rr*rr))
      else
        rr0 = sqrt(sum(b*b))
      endif

      if (rr0.lt.abstol) then
        if (out .ge. 1) then
          write (*,*) 'CG residual',rr0
          write (*,*) 'Initial guess seems exact solution in CG' 
        endif
        return
      endif

cc      iter = 0
cc      if (out.ge.2) write (*,10) iter,rr0,rr0/rr0

c     Preconditioner (zz = P rr)

      call getSolver(neq,ntot,rr,zz,matvec,igrid,bcnd,0,precout,depth1)

c     CG preparation

      rho1 = sum(rr*zz)

      pp = zz

c     CG step

      do iter= 1,niter

        rho0 = rho1

        call matvec(0,neq,ntot,pp,qq,igrid,bcnd)

        alpha = rho0/sum(pp*qq)

        x  = x  + alpha*pp

        rr = rr - alpha*qq

c     Check convergence

        mag = sqrt(sum(rr*rr))
        mag1 = mag/rr0
        if (out.ge.2) write (*,10) iter,mag,mag1
cc        if (mag1.lt.tol) exit

        if (mag < (abstol + tol*rr0)) then
          if (mag > tol*rr0) then
            if (out.gt.1) write (*,*) 'Solution in CG is in round-off.'
            if (out.gt.1) write (*,*) 'Aborting CG...'
          endif
          ierr = 0
          exit
        endif

c     Preconditioner (zz = P rr)

        call getSolver(neq,ntot,rr,zz,matvec,igrid,bcnd,0,precout
     .                ,depth1)

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

 10   format (' CG Iteration:',i4,'; Residual:',1p,1e12.4,
     .        '; Ratio:',1p,1e12.4)

      end subroutine cg

c gm
c#######################################################################
      recursive subroutine gm(neq,ntot,b,x,matvec,options,igrid,bcnd
     .                       ,guess,iout,depth)
c--------------------------------------------------------------------
c     Matrix-free Preconditioned Conjugate Gradient routine to solve
c     Ax = b. Call variables:
c       * ntot: grid dimension
c       * b,x: rhs, solution vectors
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

      integer(4) :: neq,ntot,im,iout,ierr,igrid,its,guess,depth
     .             ,bcnd(6,neq)
      real(8)    :: b(ntot),x(ntot)

      external      matvec

      type (solver_options) :: options

c Local variables

      integer(4) :: k_max,stp_test,depth1
      real(8)    :: eps,epsmac,rold,ro,eps1,gam,tt,mag,abstol,dx(ntot)

      integer(4) :: i,j,i1,k,k1,ii,jj
      integer(4) :: nn,rstrt,irstrt,precout,maxits

      double precision,allocatable,dimension(:)  :: c,s,rs
      double precision,allocatable,dimension(:,:):: hh,vv,zz

      data epsmac /1d-20/

c Begin program

      its = 0
      nn = ntot

      precout = max(0,iout-2)

      depth1 = depth + 1

      abstol = epsmac*nn

c Extract options

      maxits   = options%iter                       !Maximum number of GMRES its
      stp_test = options%stp_test                   !Stopping test type (0 -> rhs, 1-> residual)
      eps      = options%tol                        !Convergence tolerance
      k_max    = min(options%krylov_subspace,ntot)  !Maximum krylov subspace

c Allocate work arrays

      allocate(hh(k_max+1,k_max),vv(ntot,k_max+1),zz(ntot,k_max+1))
      allocate(c(k_max),s(k_max),rs(k_max+1))

c Compute initial residual vector

      if (guess.eq.0) then

c     For zero initial guess

        vv(:,1) = b(:)

        x = 0d0

      else

c     For arbitrary initial guess (vv(:,1) = b - Ax)

        call matvec(0,neq,ntot,x,vv(:,1),igrid,bcnd)

        vv(:,1) = b(:) - vv(:,1)

      endif

c Calculate restarting loops

      rstrt = min(maxits/k_max + 1,maxits)

c Calculate magnitude of initial residual

      if (stp_test.eq.1) then
        rold = sqrt(sum(vv(:,1)*vv(:,1)))
      else
        rold = sqrt(sum(b*b))
      endif

      if (rold.lt.abstol) then
        ierr = -1
        if (iout .ge. 1)
     .       write (*,*) 'Initial guess seems exact solution in GMRES'
        call killgm
        return
      endif

      eps1=eps*rold

      ro = rold

c Restarted GMRES loop

      do irstrt = 1,rstrt

        ro = sqrt(sum(vv(:,1)*vv(:,1)))

        tt = 1.0d0/ro
        vv(:,1) = vv(:,1)*tt

c      Initialize 1-st term  of rhs of hessenberg system

        rs(1) = ro

c      GMRES iteration

        do i = 1,k_max

          its = its + 1
          i1  = i + 1

c        Call preconditioner

          call getSolver(neq,ntot,vv(1,i),zz(1,i),matvec,igrid,bcnd,0
     .                  ,precout,depth1)

          call matvec(0,neq,ntot,zz(1,i),vv(1,i1),igrid,bcnd)

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

          if (    ro .le.(abstol+eps1)
     .        .or.i  .ge.min(k_max,maxits)
     .        .or.its.ge.maxits          ) exit

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
          dx(k) = zz(k,1)*tt
        enddo

        do j=2, i
          tt = rs(j)
          do k=1, nn
            dx(k) = dx(k)+tt*zz(k,j)
          enddo
        enddo

c      Form solution

        x = x + dx

c      Check convergence and restart outer loop if necessary

        if (ro.le.(abstol+eps1)) then
          if (ro > eps1 .and. ro < abstol) then
            if (iout.gt.1)
     $           write (*,*) 'Solution in GMRES is in round-off.'
            if (iout.gt.1) write (*,*) 'Aborting GMRES...'
          endif
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

        if (iout > 1 .and. rstrt > 1) then
          write (*,*)
          write (*,*) 'Restarting GMRES; restart level #',irstrt+1
          write (*,*) 
        endif

      enddo

      if (iout.eq.1) write (*,20) its,ro,ro/rold,min(rstrt,irstrt)-1
      if (iout.ge.2) write (*,*)

      options%iter_out = its
      options%tol_out  = ro/rold

c End program

      call killgm

 10   format (' GMRES Iteration:',i4,'; Residual:',1p,1e12.4,
     .        '; Ratio:',1p,1e12.4)
 20   format (' GMRES Iteration:',i4,'; Residual:',1p,1e12.4,
     .        '; Ratio:',1p,1e12.4,'; # restarts:',i2)

      contains

      subroutine killgm

        deallocate(hh,vv,zz)
        deallocate(c,s,rs)

      end subroutine killgm

      end subroutine gm

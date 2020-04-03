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
c     subroutine matvec(g_def,elem,neq,ntot,x,b,igrid,bcnd)
c
c WHERE:
c     * g_def: grid definition structure
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

      integer :: neq,ntot,igrid,iter,guess,out,depth,bcnd(6,neq)
      real(8) :: x(ntot),b(ntot)

      external    matvec

      INTERFACE
        subroutine matvec(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
          use grid
          type(grid_mg_def),pointer :: g_def
          integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
          real(8) :: x(ntot),y(ntot)
        end subroutine matvec
      END INTERFACE

c Local variables

      type (solver_unit)   :: solver_def

cc      type (solver_options):: options
cc      character*(2)  solver

c Begin program

c Read solver definition from solver hierarchy

      call readSolverHierarchy (solver_def,depth)

cc      solver  = solver_def%solver
cc      options = solver_def%options
cc      call xferSolverOptions(options,solver_def%options)

c Symmetry test

c$$$      if (options%sym_test)
c$$$     .     call symm_test(neq,options%ngrd_tst,matvec,bcnd)

c Invoke solver

      call solver_meth(neq,ntot,b,x,matvec,igrid,bcnd,guess,out
     .                ,solver_def%solver,solver_def%options,depth)

c Check convergence

      if (solver_def%options%tol_out
     .   >solver_def%options%tol.and.out.ge.2.and.my_rank==0) then
           write (*,100) solver_def%solver
     .                  ,solver_def%options%tol_out
     .                  ,solver_def%options%tol
      endif

c Output information

cc      solver_def%options%tol_out  = options%tol_out
cc      solver_def%options%iter_out = options%iter_out

      call writeSolverHierarchy(solver_def,depth)

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

      integer    ::  neq,ntot,igrid,iter,guess,out,depth,bcnd(6,neq)
      real(8)    ::  x(ntot),b(ntot)

      type (solver_options):: options
      character*(2)  solver

      external       matvec

      INTERFACE
        subroutine matvec(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
          use grid
          type(grid_mg_def),pointer :: g_def
          integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
          real(8) :: x(ntot),y(ntot)
        end subroutine matvec
      END INTERFACE

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

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: x(ntot),b(ntot)

      type (solver_options) :: options

      external     matvec

      INTERFACE
        subroutine matvec(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
          use grid
          type(grid_mg_def),pointer :: g_def
          integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
          real(8) :: x(ntot),y(ntot)
        end subroutine matvec
      END INTERFACE

c Local variables

      integer    :: niter,stp_test,depth1,ierr,mpi_comm
      integer    :: i,iter,smit,nn,vcycle,precout
      real(8)    :: alpha,beta,mag,mag1,rr0,rho0,rho1
      real(8)    :: rtol,atol,stol
      real(8)    :: rr(ntot),zz(ntot),pp(ntot),qq(ntot)

      type (solver_unit) :: solver

c Begin program

      nn = ntot

      precout = max(0,out-2)

      depth1 = depth + 1
      !Transfer grid config to PC
      call readSolverHierarchy (solver,depth1)
cc      solver%options%diag => options%diag
      solver%options%mg_grid_def => options%mg_grid_def
      call writeSolverHierarchy(solver,depth1)

c Extract options

      niter    = options%iter      !Maximum number of iterations
      stp_test = options%stp_test  !Stopping test type (0 -> rhs, 1-> residual)
      rtol     = options%tol       !Relative convergence tolerance
      atol     = options%atol      !Absolute convergence tolerance
      stol     = options%stol      !Update convergence tolerance (not used)

      if (atol == 0d0) atol = find_mach_eps()*sqrt(1.*nn)

#if defined(petsc)
      if (associated(options%mg_grid_def)) then
        mpi_comm = options%mg_grid_def%mpi_comm
      else
        mpi_comm = MPI_COMM_WORLD
      endif
#endif

c Compute initial residual

      if (guess.eq.0) then

c     For zero initial guess

        x  = 0d0
        rr = b

      else

c     For arbitrary initial guess (rr = b - Ax)

        call matvec(options%mg_grid_def,0,neq,ntot,x,rr,igrid,bcnd)

        rr = b - rr

      endif

c CG solve

      if (stp_test.eq.1) then
        rr0 = sqrt(dot3(rr,rr,mpi_comm=mpi_comm))
      else
        rr0 = sqrt(dot3(b,b,mpi_comm=mpi_comm))
      endif

      if (rr0.lt.atol) then
        if (out .ge. 1 .and. my_rank == 0) then
          write (*,*) 'CG residual',rr0
          write (*,*) 'Initial guess seems exact solution in CG' 
        endif
        return
      endif

c     Preconditioner (zz = P rr)

      call getSolver(neq,ntot,rr,zz,matvec,igrid,bcnd,0,precout,depth1)

c     CG preparation

      rho1 = dot3(rr,zz,mpi_comm=mpi_comm)

      pp = zz

c     CG step

      do iter= 1,niter

        rho0 = rho1

        call matvec(options%mg_grid_def,0,neq,ntot,pp,qq,igrid,bcnd)

        alpha = rho0/dot3(pp,qq,mpi_comm=mpi_comm)

        x  = x  + alpha*pp

        rr = rr - alpha*qq

c     Check convergence

        mag = sqrt(dot3(rr,rr,mpi_comm=mpi_comm))
        mag1 = mag/rr0
        if (out.ge.2.and. my_rank == 0) write (*,10) iter,mag,mag1

        if (mag < (atol + rtol*rr0) .or. iter == niter) then
          if (mag > rtol*rr0 .and. mag < atol) then
            if (out.gt.1.and. my_rank == 0)
     .           write (*,*) 'Solution in CG is in round-off.'
            if (out.gt.1.and. my_rank == 0) write (*,*) 'Aborting CG...'
          endif
          ierr = 0
          exit
        endif

c     Preconditioner (zz = P rr)

        call getSolver(neq,ntot,rr,zz,matvec,igrid,bcnd,0,precout
     .                ,depth1)

c     Prepare next step

        rho1 = dot3(rr,zz,mpi_comm=mpi_comm)

        beta = rho1/rho0

        pp = beta*pp + zz

      enddo

      if (out.eq.1.and. my_rank == 0) write (*,10) iter,mag,mag1
      if (out.ge.2.and. my_rank == 0) write (*,*)

      options%iter_out = iter 
      options%tol_out  = mag1

c End program

 10   format (' CG Iteration:',i4,'; Residual:',1p,1e12.4,
     .        '; Ratio:',1p,1e12.4)

      end subroutine cg

c gm
c#######################################################################
      recursive subroutine gm(neq,ntot,b,x,matvec,options,igrid,bcnd
     .                       ,guess,iout,depth)
c--------------------------------------------------------------------
c     Matrix-free Preconditioned Flexible restarted GMRES routine.
c     Has capability of dealing with singular matrices.
c
c     Call variables:
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

      integer    :: neq,ntot,im,iout,ierr,igrid,its,guess,depth
     .             ,bcnd(6,neq)
      real(8)    :: b(ntot),x(ntot)

      type (solver_options) :: options

      external      matvec

      INTERFACE
        subroutine matvec(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
          use grid
          type(grid_mg_def),pointer :: g_def
          integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
          real(8) :: x(ntot),y(ntot)
        end subroutine matvec
      END INTERFACE

c Local variables

      integer    :: k_max,stp_test,depth1
      real(8)    :: eps,rold,ro,eps1,gam,tt,mag,dx(ntot)
     .             ,smx_est,smn_est,smx_new,smn_new,snmx,snmn,csmx,csmn
     .             ,tol,atol,rtol,stol

      integer    :: i,j,i1,k,k1,ii,jj,mpi_comm
      integer    :: nn,rstrt,irstrt,precout,maxits

      real(8),allocatable,dimension(:)  :: c,s,rs,xsvmn,xsvmx
      real(8),allocatable,dimension(:,:):: hh,vv,zz

      real(8)    :: kappa,kappa_max,epsmac

      logical    :: sing_mat,sing_mat_stop,roundoff_stop

      type (solver_unit) :: solver

c Begin program

      epsmac = find_mach_eps()

      its = 0
      nn = ntot

      precout = max(0,iout-2)

      depth1 = depth + 1

      !Transfer grid config to PC
      call readSolverHierarchy (solver,depth1)
      solver%options%diag => options%diag
      solver%options%mg_grid_def => options%mg_grid_def
      call writeSolverHierarchy(solver,depth1)

      sing_mat_stop = .false.
      roundoff_stop = .false.

c Extract options

      maxits   = options%iter                       !Maximum number of GMRES its
      stp_test = options%stp_test                   !Stopping test type (0 -> rhs, 1-> residual)
      rtol     = options%tol                        !Relative convergence tolerance
      atol     = options%atol                       !Absolute convergence tolerance
      stol     = options%stol                       !Update convergence tolerance (not used)
      k_max    = options%krylov_subspace            !Maximum krylov subspace
      sing_mat = options%singular_matrix            !Tell GMRES is dealing with a singular matrix

      if (sing_mat .and. iout >= 1.and. my_rank == 0)
     .  write (*,*)
     .     '>>> GMRES SINGULAR MAT: processing singular matrix'

      if (atol == 0d0) atol = epsmac*sqrt(1d0*nn)

#if defined(petsc)
      if (associated(options%mg_grid_def)) then
        mpi_comm = options%mg_grid_def%mpi_comm
      else
        mpi_comm = MPI_COMM_WORLD
      endif
#endif

c Allocate work arrays

      allocate(hh(k_max+1,k_max),vv(ntot,k_max+1),zz(ntot,k_max+1))
      allocate(c(k_max),s(k_max),rs(k_max+1)
     .        ,xsvmn(k_max+1),xsvmx(k_max+1))

c Compute initial residual vector

      if (guess.eq.0) then

c     For zero initial guess

        vv(:,1) = b(:)

        x = 0d0

      else

c     For arbitrary initial guess (vv(:,1) = b - Ax)

        call matvec(options%mg_grid_def,0,neq,ntot,x,vv(:,1),igrid,bcnd)

        vv(:,1) = b(:) - vv(:,1)

      endif

c Calculate restarting loops

      rstrt = min(maxits/k_max + 1,maxits)

c Calculate magnitude of initial residual

      if (stp_test == 1) then
        rold = sqrt(dot3(vv(:,1),vv(:,1),mpi_comm=mpi_comm))
      else
        rold = sqrt(dot3(b,b,mpi_comm=mpi_comm))
      endif

      if (rold < atol) then
        ierr = -1
        if (iout >= 1.and. my_rank == 0)
     .       write (*,*) 'Initial guess seems exact solution in GMRES'
        call killgm
        return
      endif

      rtol = rtol*rold

      tol  = atol + rtol

      ro   = rold

c Restarted GMRES loop

      do irstrt = 1,rstrt

        if (stp_test /= 1)
     .       ro = sqrt(dot3(vv(:,1),vv(:,1),mpi_comm=mpi_comm))

        tt = 1.0d0/ro
        vv(:,1) = vv(:,1)*tt

c      Initialize 1-st term  of rhs of hessenberg system

        rs(1) = ro

c      GMRES iteration

        do i = 1,k_max

          its = its + 1
          i1  = i + 1

c        Call preconditioner

          call getSolver(neq,ntot,vv(:,i),zz(:,i),matvec,igrid,bcnd,0
     .                  ,precout,depth1)

          call matvec(options%mg_grid_def,0,neq,ntot,zz(:,i),vv(:,i1)
     .               ,igrid,bcnd)

c        Modified gram - schmidt.

          do j=1, i
            tt = dot3(vv(:,j),vv(:,i1),mpi_comm=mpi_comm)
            hh(j,i) = tt
            vv(:,i1) = vv(:,i1) - tt*vv(:,j)
          enddo

          tt = sqrt(dot3(vv(:,i1),vv(:,i1),mpi_comm=mpi_comm))
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

          gam = sqrt(hh(i,i)*hh(i,i) + hh(i1,i)*hh(i1,i))

c        If gamma is zero then any small value will do
c        Will affect only residual estimate

          if (gam .eq. 0.0d0) gam = epsmac

c        Get next plane rotation

          c(i)   = hh(i ,i)/gam
          s(i)   = hh(i1,i)/gam

          hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
          hh(i1,i)= 0d0

c        Estimate condition number of hh for SINGULAR MATRICES

          if (sing_mat) then
            if (i > 1) then
              !Max Sing. value
              call dlaic1(1,i-1,xsvmx(1:i-1),smx_est,hh(1:i-1,i),hh(i,i)
     .                   ,smx_new,snmx,csmx)

              !Min Sing. value
              call dlaic1(2,i-1,xsvmn(1:i-1),smn_est,hh(1:i-1,i),hh(i,i)
     .                   ,smn_new,snmn,csmn)

              !Condition number
              if (smn_new /= 0d0) then
                kappa = smx_new/smn_new
              else
                kappa = 0d0
              endif

              kappa_max = sqrt(ro/rold/epsmac)

              if (iout.ge.2.and. my_rank == 0) then
                write(*,*)' >>> GMRES SINGULAR MAT: Max. sing. val. ='
     .                    ,smx_new
                write(*,*)' >>> GMRES SINGULAR MAT: Min. sing. val. ='
     .                    ,smn_new
                if (smn_new /= 0d0) then
                  write(*,*)' >>> GMRES SINGULAR MAT: Condition number='
     .                      ,kappa
                endif
                write(*,*)' >>> GMRES SINGULAR MAT: Cond. number limit='
     .                     ,kappa_max
              endif

              !Update quantities
              smx_est = smx_new
              smn_est = smn_new

              xsvmx(1:i-1) = snmx*xsvmx(1:i-1)
              xsvmx(i)     = csmx
              xsvmn(1:i-1) = snmn*xsvmn(1:i-1)
              xsvmn(i)     = csmn

              sing_mat_stop = (kappa > kappa_max)
            else
              smx_est = hh(1,1)  !Initial SV is single element
              smn_est = hh(1,1)  !Initial SV is single element
              xsvmx(1) = 1d0
              xsvmn(1) = 1d0
            endif

          endif

c        Determine residual norm and test for convergence

          rs(i1) = -s(i)*rs(i)
          rs(i)  =  c(i)*rs(i)

          ro = abs(rs(i1))

          if (iout.ge.2.and. my_rank == 0) write(*,10) its,ro,ro/rold

          if (    ro  <= tol
     .        .or.i   >= min(k_max,maxits)
     .        .or.its >= maxits
     .        .or.sing_mat_stop ) exit

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
        dx = zz(:,1)*tt

        do j=2, i
          tt = rs(j)
          dx = dx + tt*zz(:,j)
        enddo

c      Store last krylov vector

        vv(:,k_max+1) = dx(:)

c      Form solution

        x = x + dx

c      Check convergence and restart outer loop if necessary

        if (ro <= tol) then
          if (ro > rtol .and. ro < atol) then
cc            roundoff_stop = .true.
            if (iout > 1 .and. my_rank == 0) write (*,*)
            if (iout.gt.1.and. my_rank == 0)
     .        write(*,*)'GMRES solution is in round-off. Terminating.'
          endif
          ierr = 0
          exit
        endif

        if (its.ge.maxits) then
          ierr = 1
          exit
        endif

        if (sing_mat_stop) then
          ierr = 2
          if (my_rank == 0) then
            if (iout > 1) write (*,*)
            write (*,*)
     .           '>>> GMRES SINGULAR MAT: found least-squares solution'
          endif
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

        if (iout > 1 .and. rstrt > 1 .and. my_rank == 0) then
          write (*,*)
          write (*,*) 'Restarting GMRES; restart level #',irstrt+1
          write (*,*) 
        endif

      enddo

      if (iout.eq.1.and. my_rank == 0)
     .     write (*,20) its,ro,ro/rold,min(rstrt,irstrt)-1
      if (iout.ge.2.and. my_rank == 0) write (*,*)

      options%iter_out = its

      if (sing_mat_stop) then
        options%tol_out  = kappa_max/kappa*tol
      elseif (roundoff_stop) then
        options%tol_out  = 0d0
      else
        options%tol_out  = ro/rold
      endif

c diag ******
cc        if (sing_mat) then
cc          call matvec(0,neq,ntot,x,dx,igrid,bcnd)
cc          dx = b - dx
cc          ro = sqrt(dot3(dx,dx,mpi_comm=mpi_comm))
cc          write (*,*) 'DIAG -- gm; Real residual=',ro
cc        endif
c diag ******

c End program

      call killgm

 10   format ('   FGMRES Iteration:',i4,'; Residual:',1p,1e12.4,
     .        '; Ratio:',1e12.4)
 20   format ('   FGMRES Iteration:',i4,'; Residual:',1p,1e12.4,
     .        '; Ratio:',1e12.4,'; # restarts:',i2)

      contains

      subroutine killgm

        deallocate(hh,vv,zz)
        deallocate(c,s,rs)

      end subroutine killgm

      end subroutine gm

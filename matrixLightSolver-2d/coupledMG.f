c TODO list:
c
c 1) Implement direct block solver for neq > 2.
c
c 2) Clean mg_setup module (including pertaining routines).
c
c 3) Document i/o variables in all subroutines consistently.
c
c 4) Implement general robin boundary conditions. This requires
c    defining a new structure that contains information of
c    robin constants alpha, beta, and gamma for each boundary:
c
c         alpha*u + beta*u' = gamma
c
c    The following could work:
c
c        structure defBC:
c           logical periodic
c           real(8) alpha
c           real(8) beta
c           real(8) gamma
c
c    Implementation: define structure defBC in module mlsolversetup, 
c    and declare structure array bcond(4,neq) of type defBC where used.

c module mg_internal
c######################################################################
      module mg_internal

        use mg_setup

        double precision, allocatable, dimension(:,:) :: diag

        integer,dimension(:),allocatable :: istart,ntotv

      contains

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,ngrid,fpointers)

c     -----------------------------------------------------------------
c     Initializes MG and creates 2D uniform grid
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer :: neq,ngrid
      logical :: fpointers

c     Local variables

      integer :: i,alloc_stat

c     Begin program

      fpointers = .false.

      allocate(istart(ngrid),ntotv(ngrid), STAT = alloc_stat)

      !Successful memory allocation
      if (alloc_stat == 0) then
        istart(1) = 1
        do i = 2,ngrid
          istart(i) = neq*nxvp(i-1)*nyvp(i-1) + istart (i-1)
          ntotv (i) = neq*nxvp(i)*nyvp(i)
        enddo
        fpointers = .true.
      endif

c     End program

      end subroutine allocPointers

c     deallocPointers
c     #################################################################
      subroutine deallocPointers(fpointers)

c     -----------------------------------------------------------------
c     Initializes MG and creates 2D uniform grid
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers

      if (fpointers) deallocate(istart,ntotv)

      end subroutine deallocPointers

      end module mg_internal

c mg
c#######################################################################
      recursive subroutine mg(neq,ntot,y,x,matvec,options,igrid,bcond
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free coupled MG routine to solve
c     Ax = y. Call variables:
c       * ntot: grid dimension
c       * y,x: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0->no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer*4    neq,ntot,igrid,guess,out,depth,bcond(4,neq)
      real*8       x(ntot),y(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer          :: iter,igridmin,vcyc
      integer          :: orderres,orderprol,alloc_stat

      integer          :: guess2,outc
      integer          :: izero,i,j,ii,ivcyc

      double precision :: xx(2*ntot),yy(2*ntot),wrk(2*ntot)
      double precision :: rr0,rr1,mag,mag1,mgtol
      double precision :: dummy(ntot),rr(ntot)

      logical          :: fdiag,fpointers

c Begin program

      igridmin = options%igridmin
      vcyc     = options%vcyc
      mgtol    = options%tol
      orderres = options%orderres
      orderprol= options%orderprol
      fdiag    = options%fdiag

      if (out.ge.2.and.igridmin.lt.igrid) write (*,5)

      outc = out
      if (igridmin.lt.igrid) outc = out - 2

c Set pointers

      call allocPointers(neq,igrid,fpointers)

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

c Find diagonal for smoothers

      if (fdiag) then

        allocate(diag(neq,2*ntot),STAT = alloc_stat)

        if (alloc_stat.ne.0) then !Diagonal is known (failed alloc)
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          call find_mf_diag_neq(neq,ntot,matvec,igrid,bcond,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Initialize local solution vector (xx) and local r.h.s. (yy)

      xx  = 0d0
      yy  = 0d0
      wrk = 0d0

c memcheck ****
      dummy = 0d0
c memcheck ****

      if (guess.eq.0) then
        x = 0d0
        yy(istart(igrid):istart(igrid) + ntot - 1) = y(:)
      else
        xx(istart(igrid):istart(igrid) + ntot - 1) = x(:)
        yy(istart(igrid):istart(igrid) + ntot - 1) = y(:)
      endif

c Compute initial residual and check convergence

      if (guess.eq.0) then
        rr0 = sqrt(sum(y*y))
      else
        call matvec(0,ntot,x,dummy,igrid,bcond)
        rr = y - dummy
        rr0 = sqrt(sum(rr*rr))
      endif

      if (rr0.lt.1d-16*ntot) then
        if (out.ge.1) write (*,*) 'Initial solution seems exact in MG'
        call killmg
        return
      endif

      rr1 = rr0

c Start V-cycle

      guess2 = guess

      do ivcyc = 1,vcyc

        if (outc.ge.1.and.igridmin.lt.igrid) write (*,7) ivcyc

c     Perform Vcycle recursively

        if (igrid.eq.igridmin) then
          call solve (igrid)
          exit
        else
          call vcycle(igrid)
        endif

c     Check MG convergence

        call matvec(0,ntot,xx(istart(igrid)),dummy,igrid,bcond)

        rr   = y - dummy
        mag  = sqrt(sum(rr*rr))
        mag1 = mag/rr1

        if (out.ge.3) then
          write (*,10) mag,mag1
        elseif (out.ge.2) then
          write (*,20) mag,mag1,ivcyc
        endif

        rr1 = mag

        if (mag/rr0 < mgtol .or. mag < 1d-20*ntot) exit

      enddo

c Map solution from local vector xx to external vector x

      x(:) = xx(istart(igrid):istart(igrid) + ntot - 1)

c MG convergence info

      mag1 = mag/rr0

      if (igridmin.lt.igrid) then
        if (out.eq.1) then
          write (*,20) mag,mag1,min(ivcyc,vcyc)
        elseif (out.ge.2.and.vcyc.gt.1) then
          write (*,*) 
          write (*,*) 'Final MG convergence info:'
          write (*,20) mag,mag1,min(ivcyc,vcyc)
          write (*,*)
        endif
      endif

c End program

      options%tol_out  = mag1
      options%iter_out = min(ivcyc,vcyc)

      call killmg

      return

 5    format (/,' MG method output:')
 7    format (/,' MG V-cycle #:',i3)
 10   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4)
 20   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4,
     .          '; V-cycle #:'i3)

      contains

c     vcycle
c     ###################################################################
        recursive subroutine vcycle(igr)

          implicit none
          integer :: igr,isigm,isig,nn

          real(8) :: norm

c       Begin program

          nn    = ntotv(igr)
          isig  = istart(igr)
          isigm = istart(igr-1)

c       Relax error/solution on grid number igr/igrid (find new xx)

          call solve(igr)

c       Evaluate residual (ie wrk = yy - A xx = yy - wrk )

          call matvec(0,nn,xx(isig),wrk(isig),igr,bcond)

          wrk(isig:isig+nn-1) = yy(isig:isig+nn-1) - wrk(isig:isig+nn-1)

c       Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) down a grid

          call crestrict (neq
     .                   ,yy(isigm),ntotv(igr-1),nxvp(igr-1),nyvp(igr-1)
     .                   ,wrk(isig),ntotv(igr  ),nxvp(igr  ),nyvp(igr  )
     .                   ,orderres,igr,bcond,4d0)

          guess2 = 0

c       If on coarsest grid, solve for error, else descend a grid level

          if (igr-1.eq.igridmin) then
cc            call solve (igr-1)
            call coarseSolve (igr-1)
          else
            call vcycle(igr-1)
          endif

c       Cycle back up to grid igr updating errors (xx)
c       with fixed R.H.S. (yy)

          guess2 = 1

c       Prolong error (wrk = P * xx_1) on grid 1 to grid 2

          call cprolong (neq
     .                  ,wrk(isig),ntotv(igr  ),nxvp(igr  ),nyvp(igr  )
     .                  ,xx(isigm),ntotv(igr-1),nxvp(igr-1),nyvp(igr-1)
     .                  ,orderprol,igr-1,bcond)

c       Update existing error on grid 2 (i.e. xx_2): xx_2 = xx_2 + wrk 

          xx(isig:isig+nn-1) = xx(isig:isig+nn-1) + wrk(isig:isig+nn-1)

cc          norm = sqrt(sum(xx*xx))
cc          norm = sqrt(sum(wrk*wrk))/norm
cc          write (*,*) norm
c       Relax updated error on igr (i.e. xx_igr)

          call solve(igr)

        end subroutine vcycle

c       solve
c       ###################################################################
        recursive subroutine solve(igr)

          implicit none
          integer :: igr,nn,isig,depth1

          nn   = ntotv(igr)
          isig = istart(igr)

          depth1 = depth + 1

          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcond
     .                  ,guess2,outc,depth1)

        end subroutine solve

c       coarseSolve
c       ###################################################################
        recursive subroutine coarseSolve(igr)

          implicit none
          integer :: igr,nn,isig,depth1,depth2

          nn   = ntotv(igr)
          isig = istart(igr)

          call countSolverLevels(depth2)

          depth1 = min(depth+2,depth2)

          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcond
     .                  ,guess2,outc,depth1)

        end subroutine coarseSolve

c       killmg
c       ###################################################################
        subroutine killmg

          implicit none

          if (fdiag)  deallocate(diag)
          call deallocPointers(fpointers)

        end subroutine killmg

      end subroutine

c jb
c#######################################################################
      recursive subroutine jb(neq,ntot,rr,zz,matvec,options,igrid,bcond
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free Jacobi routine to solve Ax = b. Call variables:
c       * ntot: grid dimension
c       * rr,zz: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0 -> no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer          :: neq,ntot,igrid,guess,out,depth,bcond(4,neq)
      double precision :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer          :: iter,alloc_stat,isig
      double precision :: omega0,omega10,omega01,tol
      logical          :: fdiag,fpointers

      integer*4    i,j,itr,nn,ieq,nx,ny
      integer*4    ii,iii,iip,iim,jj,jjp,jjm
      integer*4    ig,ipg,img,jg,jpg,jmg,iig
      real*8       mag0,mag1,mag,yy(ntot),delta
      double precision :: aw,ae,an,as

      double precision,allocatable, dimension(:)  :: dummy,rhs

c Begin program

      iter   = options%iter
      omega0 = options%omega
      omega10= options%omega10
      omega01= options%omega01
      tol    = options%tol
      fdiag  = options%fdiag

      if (out.ge.2) write (*,*)

c Allocate pointers

      call allocPointers(neq,igrid,fpointers)

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

      isig = istart(igrid)

c Find diagonal for smoothers

      if (fdiag) then

        allocate(diag(neq,2*ntot),STAT = alloc_stat)

        if (alloc_stat.ne.0) then !Diagonal is known (failed alloc)
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          call find_mf_diag_neq(neq,ntot,matvec,igrid,bcond,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Allocate internal arrays

      allocate(dummy(neq),rhs(neq))

c Jacobi iteration

      mag = 0d0

      nn = ntot

      if (guess.eq.0) then
        do ii=1,nn
          zz(ii) = 0d0
        enddo
      endif

      do itr=1,iter

        call matvec(0,ntot,zz,yy,igrid,bcond)

        mag1 = mag
        mag  = 0d0

        select case (neq)
        case (1)

          if (omega10.eq.0d0.and.omega01.eq.0d0) then

            do ii = 1,nn
              ig = ii + isig - 1
              rhs(neq) = rr(ii) - yy(ii)
              dummy(neq) = rhs(neq)/diag(neq,ig)
              zz(ii) = zz(ii) + omega0*dummy(neq)
              mag = mag + rhs(neq)**2
            enddo

          else

            nx = nxvp(igrid)
            ny = nyvp(igrid)

            do i = 1,nx
              do j = 1,ny
                call shiftcoeffs (i,j,bcond(1,1))

                call shiftIndices(i,j,ii,iim,iip,jj,jjm,jjp,nx,ny,1
     .                            ,bcond(1,1))
                call shiftIndices(i,j,ig,img,ipg,jg,jmg,jpg,nx,ny,isig
     .                            ,bcond(1,1))

                zz(ii) =zz(ii)
     .                 +omega0*     (rr(ii) -yy(ii) )/diag(neq,ig )
     .                 +omega10*(ae*(rr(iip)-yy(iip))/diag(neq,ipg)
     .                          +aw*(rr(iim)-yy(iim))/diag(neq,img))
     .                 +omega01*(an*(rr(jjp)-yy(jjp))/diag(neq,jpg)
     .                          +as*(rr(jjm)-yy(jjm))/diag(neq,jmg))

                mag = mag + (rr(ii)-yy(ii))**2
              enddo
            enddo

          endif

        case (2)

          if (omega10.ne.0d0.or.omega01.ne.0d0) then
            write (*,*)'Weighed jacobi not available in coupled Jacobi'
            write (*,*)'Aborting...'
            stop
          endif

          do ii = 1,nn/neq

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              rhs(ieq) = rr(iii) - yy(iii)
            enddo

            iig = neq*(ii - 1) + isig - 1
            delta = diag(1,iig+1)*diag(2,iig+2)
     .             -diag(2,iig+1)*diag(1,iig+2)
            dummy(1) = (rhs(1)*diag(2,iig+2)-rhs(2)*diag(2,iig+1))/delta
            dummy(2) = (rhs(2)*diag(1,iig+1)-rhs(1)*diag(1,iig+2))/delta

            iii = neq*(ii - 1)
            do ieq=1,neq
              zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
            enddo

            mag = mag + sum(rhs**2)
          enddo

        case default

          write (*,*) 'Jacobi smoothing for neq > 2 not implemented yet'
          stop

        end select

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) min(itr,iter),mag,mag/mag0

      options%iter_out = min(itr,iter)
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag)  deallocate(diag)
      call deallocPointers(fpointers)

      return

 10   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

      contains

c     shiftcoeffs
c     ###############################################################
      subroutine shiftcoeffs(i,j,bcond)
      implicit none
c     ---------------------------------------------------------------
c     Shifts indices in a 5pt stencil for both arrays (isig = 0) 
c     and MG vectors (isig = pointer to grid), accounting for boundary 
c     conditions as defined in integer array bcond.
c     ---------------------------------------------------------------

c     Call variables

      integer          :: i,j,bcond(4)

c     Local variables

c     Begin program

      an = 1d0
      as = 1d0
      ae = 1d0
      aw = 1d0

      if (i+1 > nx .and. bcond(4) /= 0) ae = 0d0  !Zero out if not periodic
      if (i-1 < 1  .and. bcond(3) /= 0) aw = 0d0
      if (j+1 > ny .and. bcond(2) /= 0) an = 0d0
      if (j-1 < 1  .and. bcond(1) /= 0) as = 0d0

c     End program

      end subroutine

      end subroutine

c gs
c#######################################################################
      recursive subroutine gs(neq,ntot,rr,zz,matvec,options,igrid,bcond
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free Jacobi routine to solve Ax = b. Call variables:
c       * ntot: grid dimension
c       * rr,zz: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0 -> no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer          :: neq,ntot,igrid,guess,out,depth,bcond(4,neq)
      double precision :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer          :: iter,alloc_stat,isig
      double precision :: omega0,tol
      logical          :: fdiag,fpointers,fbcond

      integer*4    i,j,itr,nn,ieq,nx,ny
      integer*4    ii,iii,iip,iim,jjp,jjm
      integer*4    ig,ipg,img,jpg,jmg,iig
      real*8       mag0,mag1,mag,yy(ntot),delta

      double precision,allocatable, dimension(:)  :: dummy,rhs

c Begin program

      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag

      if (out.ge.2) write (*,*)

c Set pointers

      call allocPointers(neq,igrid,fpointers)

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

      isig = istart(igrid)

c Find diagonal for smoothers

      if (fdiag) then

        allocate(diag(neq,2*ntot),STAT = alloc_stat)

        if (alloc_stat.ne.0) then !Diagonal is known (failed alloc)
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          call find_mf_diag_neq(neq,ntot,matvec,igrid,bcond,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Allocate internal arrays

      allocate(dummy(neq),rhs(neq))

c SGS iteration

      nn = ntot

      if (guess.eq.0) then
        do ii=1,nn
          zz(ii) = 0d0
        enddo
      endif

      do itr=1,iter

        mag1 = mag
        mag = 0d0

        select case (neq)
        case (1)

c       Forward pass

          do ii = 1,nn
            ig = ii + isig - 1
            call matvec(ii,ntot,zz,yy,igrid,bcond)
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ig)
            zz(ii) = zz(ii) + omega0*dummy(neq)
          enddo

c       Backward pass

          do ii = nn,1,-1
            ig = ii + isig - 1
            call matvec(ii,ntot,zz,yy,igrid,bcond)
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ig)
            zz(ii) = zz(ii) + omega0*dummy(neq)
            mag = mag + rhs(neq)**2
          enddo

        case (2)

c       Forward pass

          do ii = 1,nn/neq

            call matvec(ii,ntot,zz,yy,igrid,bcond)

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              rhs(ieq) = rr(iii) - yy(iii)
            enddo

            iig = neq*(ii - 1) + isig - 1
            delta = diag(1,iig+1)*diag(2,iig+2)
     .             -diag(2,iig+1)*diag(1,iig+2)
            dummy(1) = (rhs(1)*diag(2,iig+2)-rhs(2)*diag(2,iig+1))/delta
            dummy(2) = (rhs(2)*diag(1,iig+1)-rhs(1)*diag(1,iig+2))/delta

            iii = neq*(ii-1)
            do ieq=1,neq
              zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
            enddo

          enddo

c       Backward pass

          do ii = nn/neq,1,-1

            call matvec(ii,ntot,zz,yy,igrid,bcond)

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              rhs(ieq) = rr(iii) - yy(iii)
            enddo

            iig = neq*(ii - 1) + isig - 1
            delta = diag(1,iig+1)*diag(2,iig+2)
     .             -diag(2,iig+1)*diag(1,iig+2)
            dummy(1) = (rhs(1)*diag(2,iig+2)-rhs(2)*diag(2,iig+1))/delta
            dummy(2) = (rhs(2)*diag(1,iig+1)-rhs(1)*diag(1,iig+2))/delta

            iii = neq*(ii-1)
            do ieq=1,neq
              zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
            enddo

            mag = mag + sum(rhs**2)
cc            mag = mag + rhs(2)**2
          enddo

        case default

          write (*,*) 'SGS smoothing for neq > 2 not implemented yet'
          stop

        end select

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) min(itr,iter),mag,mag/mag0

      options%iter_out = min(itr,iter)
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag)  deallocate(diag)
      call deallocPointers(fpointers)

      return

 10   format (' GS Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' GS Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

      end subroutine

*deck cprolong
c######################################################################
      subroutine cprolong(neq,xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc,order
     .                   ,igc,bcond)
c----------------------------------------------------------------------
c     This is a prolongation routine for MG. 
c
c     If order = 0, it employs simple injection.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none            ! For safe Fortran

c Call variables

      integer :: neq,ntotc,nxc,nyc,ntotf,nxf,nyf,order,igc
     .          ,bcond(4,neq)
      double precision ::  xc(ntotc),xf(ntotf)

c Local variables
 
      real*8       xxc(ntotc/neq),xxf(ntotf/neq)

      integer*4    ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

cc      if (order.eq.0) then

cc      do jc = 1,nyc

cc        do ic = 1,nxc
cc          jf = 2*jc
cc          if = 2*ic

cc          do ieq = 1,neq
cc            iic = ieq + neq*(ic-1 + nxc*(jc-1))
cc            iif = ieq + neq*(if-1 + nxf*(jf-1))

cc            xf(iif                ) = xc(iic) 
cc            xf(iif - neq          ) = xc(iic) 
cc            xf(iif - neq*nxf      ) = xc(iic) 
cc            xf(iif - neq*nxf - neq) = xc(iic) 
cc          enddo
cc        enddo

cc      enddo

cc      else

c Interpolation

        nntotc = ntotc/neq
        nntotf = ntotf/neq

        do ieq = 1,neq

c     Unpack vector

          do i = 1,nntotc
            xxc(i) = xc(neq*(i-1)+ieq)
          enddo

c     Use scalar prolongation

          call prolong(xxf,nntotf,nxf,nyf
     .                ,xxc,nntotc,nxc,nyc
     .                ,order,igc,bcond(1,ieq))

c     Repack prolonged vector

          do i=1,nntotf
            xf(neq*(i-1)+ieq)=xxf(i)
          enddo

        enddo

cc      endif

c End program

      return
      end subroutine

*deck crestrict
c######################################################################
      subroutine crestrict(neq,xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf,order
     .                    ,igf,bcond,volf)
c----------------------------------------------------------------------
c     This is a restriction routine for MG. 
c
c     If order = 0, it employs simple agglomeration.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer ::   ntotc,nxc,nyc,ntotf,nxf,nyf,order,igf,neq
     .            ,bcond(4,neq)

      double precision :: xc(ntotc),xf(ntotf),volf

c Local variables
 
      real*8       xxc(ntotc/neq),xxf(ntotf/neq)

      integer*4    ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

c Agglomeration

cc      if (order.eq.0) then

cc        do jc = 1,nyc

cc          do ic = 1,nxc
cc            jf = 2*jc
cc            if = 2*ic

cc            do ieq = 1,neq
cc              iic = ieq + neq*(ic-1 + nxc*(jc-1))
cc              iif = ieq + neq*(if-1 + nxf*(jf-1))

cc              xc(iic) = (xf(iif)        +xf(iif        -neq)
cc     .                 + xf(iif-nxf*neq)+xf(iif-nxf*neq-neq))/4d0*volf
cc            enddo
cc          enddo

cc        enddo

cc      else

c Interpolation

        nntotc = ntotc/neq
        nntotf = ntotf/neq

        do ieq = 1,neq

c     Unpack vector

          do i = 1,nntotf
            xxf(i) = xf(neq*(i-1)+ieq)
          enddo

c     Use scalar restriction

          call restrict(xxc,nntotc,nxc,nyc
     .                 ,xxf,nntotf,nxf,nyf
     .                 ,order,igf,bcond(1,ieq),volf)

c     Repack restricted vector

          do i=1,nntotc
            xc(neq*(i-1)+ieq)=xxc(i)
          enddo

        enddo

cc      endif

c End program

      return
      end subroutine

*deck prolong
c######################################################################
      subroutine prolong(xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc,order
     .                  ,igc,bbcond)
c----------------------------------------------------------------------
c     This is a prolongation routine for MG. 
c
c     If order = 0, it employs simple injection.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none            ! For safe Fortran

c Call variables

      integer      ntotc,nxc,nyc,ntotf,nxf,nyf,order,igc,bbcond(4)
      real*8       xc(ntotc),xf(ntotf)

c Local variables

      double precision :: xxf,yyf,ff

      double precision :: xx(nxc+2),yy(nyc+2),uc(0:nxc+1,0:nyc+1)

      integer          :: ic,if,jc,jf,ii,i,igf,iic,iif

c Interpolation

      integer*4    kx,ky,nx,ny,dim,flg
      real*8, dimension(:),allocatable:: tx,ty,work
      real*8, dimension(:,:),allocatable:: bcoef

      real*8       db2val
      external     db2val

c Begin program

c Injection

      if (order.eq.0) then

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic
            iic = ic + nxc*(jc-1)
            iif = if + nxf*(jf-1)

            xf(iif)           = xc(iic) 
            xf(iif-1)         = xc(iic) 
            xf(iif - nxf)     = xc(iic) 
            xf(iif - nxf - 1) = xc(iic) 
          enddo
        enddo

      else

c Interpolation

c     Determine current coarse grid

c     Map vectors into arrays

        xx(1:nxc+2) = xl(0:nxc+1,igc)
        yy(1:nyc+2) = yl(0:nyc+1,igc)

        call mapMGVectorToArray(nxc,nyc,xc,uc,1,bbcond)

c     Calculate interpolation

        flg = 0
        kx = order+1
        ky = order+1
        nx = nxc + 2
        ny = nyc + 2
        dim = nx*ny + max(2*kx*(nx+1),2*ky*(ny+1))

        allocate(tx(nx+kx))
        allocate(ty(ny+ky))
        allocate(work(dim))
        allocate(bcoef(nx,ny))

        call db2ink(xx,nx,yy,ny,uc,nx,kx,ky,tx,ty,bcoef,work,flg)

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic
            iif = if + nxf*(jf-1)

            igf = igc+1         !Consider next finer grid

            xxf = xl(if,igf)
            yyf = yl(jf,igf)

            xf(iif) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

            xxf = xl(if-1,igf)
            yyf = yl(jf  ,igf)

            xf(iif-1) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

            xxf = xl(if-1,igf)
            yyf = yl(jf-1,igf)

            xf(iif-1-nxf) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

            xxf = xl(if  ,igf)
            yyf = yl(jf-1,igf)

            xf(iif-nxf) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

          enddo
        enddo

        deallocate(tx)
        deallocate(ty)
        deallocate(work)
        deallocate(bcoef)

      endif

c End program

      return
      end subroutine

*deck restrict
c######################################################################
      subroutine restrict(xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf,order
     .                   ,igf,bbcond,volf)
c----------------------------------------------------------------------
c     This is a restriction routine for MG. 
c
c     If order = 0, it employs simple agglomeration.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf,order,igf,bbcond(4)

      real*8       xc(ntotc),xf(ntotf),volf

c Local variables
 
      real*8       uf(0:nxf+1,0:nyf+1)

      real*8       xxc,yyc,xx(nxf+2),yy(nyf+2)

      integer*4    ic,if,jc,jf,iif,iic,i,igc

c Interpolation

      integer*4    kx,ky,nx,ny,dim,flg
      real*8, dimension(:)  ,allocatable:: tx,ty,work
      real*8, dimension(:,:),allocatable:: bcoef

      real*8       db2val
      external     db2val

c Begin program

c Agglomeration

      if (order.eq.0) then

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic
            iic = ic + nxc*(jc-1)
            iif = if + nxf*(jf-1)

            xc(iic) = (xf(iif)     + xf(iif-1)
     .               + xf(iif-nxf) + xf(iif-nxf-1) )/4d0*volf
          enddo
        enddo

      else

c Interpolation

c     Map vectors into arrays

        xx(1:nxf+2) = xl(0:nxf+1,igf)
        yy(1:nyf+2) = yl(0:nyf+1,igf)

        call mapMGVectorToArray(nxf,nyf,xf,uf,1,bbcond)

c     Calculate interpolation

        flg = 0
        kx = order + 1
        ky = order + 1
        nx = nxf+2
        ny = nyf+2
        dim = nx*ny + max(2*kx*(nx+1),2*ky*(ny+1))

        allocate(tx(nx+kx))
        allocate(ty(nx+ky))
        allocate(work(dim))
        allocate(bcoef(nx,ny))

        call db2ink(xx,nx,yy,ny,uf,nx,kx,ky,tx,ty,bcoef,work,flg)

        do jc = 1,nyc
          do ic = 1,nxc

            iic = ic + nxc*(jc-1)

            igc = igf-1         !Consider next coarser grid

            xxc = xl(ic,igc)
            yyc = yl(jc,igc)
            
            xc(iic) = 
     .           volf*db2val(xxc,yyc,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

          enddo
        enddo

        deallocate(tx)
        deallocate(ty)
        deallocate(work)
        deallocate(bcoef)

      endif

c End program

      return
      end subroutine

c find_mf_diag_neq
c####################################################################
      subroutine find_mf_diag_neq(neq,ntot,matvec,ngrid,bbcnd,diag1)
c--------------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec.
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer*4      neq,ntot,ngrid,bbcnd(4,neq)

      double precision :: diag1(neq,2*ntot)

      external       matvec

c Local variables

      real*8         x1(ntot),dummy(ntot)
      integer*4      ii,jj,nn
      integer*4      nx,ny,igrid,ieq,alloc_stat
      logical        fpointers

c Begin program

c Allocate pointers

      call allocPointers(neq,ngrid,fpointers)

c Form diagonal

      do igrid = ngrid,2,-1

c     Form diagonal terms for smoother

        nx  = nxvp(igrid)
        ny  = nyvp(igrid)

        nn = neq*nx*ny

        x1(1:nn) = 0d0

c Finds block diagonals for neq equations.

        do ii = 1,nn/neq

          jj  = (ii-1)*neq + istart(igrid) - 1

          do ieq = 1,neq

c         Find column vector corresponding to grid node ii and equation ieq

            call findBaseVector(ii,ieq,neq,nn,igrid,bbcnd(1,ieq),x1,1d0)

            call matvec(ii,nn,x1,dummy,igrid,bbcnd)

            call findBaseVector(ii,ieq,neq,nn,igrid,bbcnd(1,ieq),x1,0d0)

c         Fill diagonal

            diag1(ieq,jj+1:jj+neq) = dummy((ii-1)*neq+1:ii*neq)

          enddo

        enddo

      enddo

c Deallocate pointers

      call deallocPointers(fpointers)

c End program

      return
      end subroutine

c symm_test
c####################################################################
      recursive subroutine symm_test(neq,igrid,matvec,bcond)
c--------------------------------------------------------------------
c     Performs symmetry test on matvec.
c--------------------------------------------------------------------

      use mg_setup

      implicit none     !For safe fortran

c Call variables

      integer*4      neq,igrid,bcond(4,*)

      external       matvec

c Local variables

      double precision,allocatable,dimension(:)::x1,dummy,dummy2

      real(8) ::   dd1,dd2,error
      integer ::   nx,ny,ii,iii,jj,jjj,ieq,nn

c Begin program

      nx = nxvp(igrid)
      ny = nyvp(igrid)

      nn = neq*nx*ny

      write (*,*) 'Performing symmetry test of system matrix ',
     .            'on grid:',nx,'x',ny,'...'

      allocate(x1(nn),dummy(nn),dummy2(nn))

      error = 0d0

      do ii = 1,nn
        x1    (ii) = 0d0
        dummy (ii) = 0d0
        dummy2(ii) = 0d0
      enddo

      do iii = 1,nn

        ieq = mod(iii,neq)
        if (ieq.eq.0) ieq = neq

        ii = 1 + (iii - ieq)/neq

c     Find column vector iii = neq*(ii-1) + ieq

        call findBaseVector(ii,ieq,neq,nn,igrid,bcond(1,ieq),x1,1d0)

        call matvec(0,nn,x1,dummy,igrid,bcond)

        call findBaseVector(ii,ieq,neq,nn,igrid,bcond(1,ieq),x1,0d0)

c     Compare column vector iii with corresponding row vector (intersect in
c     diagonal)

        do jjj = iii,nn

          ieq = mod(jjj,neq)
          if (ieq.eq.0) ieq = neq

          jj = 1 + (jjj - ieq)/neq

          call findBaseVector(jj,ieq,neq,nn,igrid,bcond(1,ieq),x1,1d0)

          call matvec(ii,nn,x1,dummy2,igrid,bcond)

          call findBaseVector(jj,ieq,neq,nn,igrid,bcond(1,ieq),x1,0d0)

          dd1 = abs(dummy(jjj) - dummy2(iii))
          if(    abs(dummy (jjj)).gt.1d-15
     .       .or.abs(dummy2(iii)).gt.1d-15) then
            write (*,15) jjj,iii,dummy(jjj)
     .                  ,iii,jjj,dummy2(iii)
     .                  ,dd1
     .                  ,100*dd1/max(abs(dummy(jjj)),abs(dummy2(iii)))
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

c findBaseVector
c####################################################################
      subroutine findBaseVector(ii,ieq,neq,ntot,igrid,bcond,x1,coef)
c--------------------------------------------------------------------
c     Finds base vector corresponding to grid node ii and equation ieq.
c     Assumes that x1 has been initialized to zero elsewhere.
c--------------------------------------------------------------------

      use mg_setup

      implicit none      !For safe fortran

c Call variables

      integer          :: neq,ii,ieq,ntot,igrid,bcond(4)

      double precision :: x1(ntot),coef

c Local variables

      integer          :: nx,ny,iii,imin,imax,jmin,jmax

c Begin program

      nx  = nxvp(igrid)
      ny  = nyvp(igrid)

c Locate equation and grid point

cc      ieq = mod(iii,neq)
cc      if (ieq.eq.0) ieq = neq
cc
cc      ii = 1 + (iii - ieq)/neq

      call limits(ii,nx,ny,imin,imax,jmin,jmax)

      if (bcond(3) == 0 .or. bcond(4) == 0) then
        if (imin.eq.1) then
          iii = neq*(ii-1) + ieq
          x1(iii) = coef
          iii = neq*(ii+nx-2) + ieq
          x1(iii) = coef
        elseif (imin.eq.nx) then
          iii = neq*(ii-1) + ieq
          x1(iii) = coef
          iii = neq*(ii-nx) + ieq
          x1(iii) = coef
        else
          iii = neq*(ii-1) + ieq
          x1(iii) = coef
        endif
      elseif (bcond(1) == 0 .or. bcond(2) == 0) then
        if (jmin.eq.1) then
          iii = neq*(ii-1) + ieq
          x1(iii) = coef
          iii = neq*(ii-1+nx*(ny-1)) + ieq
          x1(iii) = coef
        elseif (jmin.eq.ny) then
          iii = neq*(ii-1) + ieq
          x1(iii) = coef
          iii = neq*(ii-1-nx*(ny-1)) + ieq
          x1(iii) = coef
        else
          iii = neq*(ii-1) + ieq
          x1(iii) = coef
        endif
      else
        iii = neq*(ii-1) + ieq
        x1(iii) = coef
      endif

c End program

      return
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

c shiftIndices
c####################################################################
      subroutine shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,isig,bcond)
      implicit none
c--------------------------------------------------------------------
c     Shifts indices in a 5pt stencil for both arrays (isig = 0) 
c     and MG vectors (isig = pointer to grid), accounting for boundary 
c     conditions as defined in integer array bcond.
c--------------------------------------------------------------------

c Call variables

      integer :: i,j,ii,jj,im,ip,jm,jp,nx,ny,isig,bcond(4)

c Local variables

      integer :: nxx,row

c Functions

      row(i,j) = i + nx*(j-1) + isig - 1

c Begin program

      if (bcond(3) == 0 .or. bcond(4) == 0) then
        call    periodicBC(i,nx,ii,im,ip)
      else
        call nonperiodicBC(i,nx,ii,im,ip,bcond(3),bcond(4))
      endif

      if (bcond(1) == 0 .or. bcond(2) == 0) then
        call    periodicBC(j,ny,jj,jm,jp)
      else
        call nonperiodicBC(j,ny,jj,jm,jp,bcond(1),bcond(2))
      endif

      if (isig.gt.0) then    !Transform to MG vector coordinates
        ip = min(row(ip,jj),row(nx,ny))
        im = max(row(im,jj),row(1 , 1))
        jp = min(row(ii,jp),row(nx,ny))
        jm = max(row(ii,jm),row(1 , 1))
        ii = row(ii,jj)
        jj = ii
      endif

c End program

      contains

      subroutine periodicBC(i,nx,ii,im,ip)

      implicit none

      integer :: i,ii,im,ip,nx

      if (i.lt.1) then
        ii = nx + i - 1 
      elseif (i.gt.nx) then
        ii = i - nx + 1 
      else
        ii = i 
      endif

      if (i.eq.nx) then
        ip = 2 
      else
        ip = ii + 1
      endif

      if (i.eq.1 ) then
        im = nx - 1 
      else
        im = ii - 1
      endif

      end subroutine

      subroutine nonperiodicBC(i,nx,ii,im,ip,bcs1,bcs2)

      implicit none

      integer :: i,ii,im,ip,nx,bcs1,bcs2

      ii = i
      if (i == 1 .and. bcs1 == 2) then
cc        im = i
        im = i+1
      elseif (i == 0) then
        im = i
      else
        im = i-1
      endif

      if (i == nx .and. bcs2 == 2) then
cc        ip = i
        ip = i-1
      elseif (i == nx+1) then
        ip = i
      else
        ip = i+1
      endif

      end subroutine

      end subroutine

c setBoundaryConditions
c####################################################################
      subroutine setBoundaryConditions(array,imin,imax,jmin,jmax
     .                                ,nx,ny,bcs)

c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * array  -> real array with ghost-nodes
c       * nx,ny  -> domain size
c       * bcs    -> integer array of size 4 containing BC setup:
c           + bcs(1) ---> bottom
c           + bcs(2) ---> top
c           + bcs(3) ---> left
c           + bcs(4) ---> right
c         If bcs = 0  --> periodic BC's
c         If bcs = 1  --> dirichlet BC's (always zero because MG deals with
c                                         error updates)
c         If bcs = 2  --> Neumann BC's
c         If bcs = 3  --> second-order derivative = 0 (linear extrapolation)
c--------------------------------------------------------------------

      implicit none       !For safe fortran

c Call variables

      integer*4  nx,ny,imin,imax,jmin,jmax,bcs(4)
      real*8     array(0:nx+1,0:ny+1)

c Local variables

      integer*4  i,j

c Begin program

c Bottom

      if (bcs(1).eq.0) then
        do i=imin,imax
          array(i,0)    = array(i,ny-1)
        enddo
      elseif (bcs(1).eq.1) then
        do i=imin,imax
          array(i,0)    = 0d0
        enddo
      elseif (bcs(1).eq.2) then
        do i=imin,imax
cc          array(i,0)    = array(i,2)
          array(i,0)    = array(i,1)
        enddo
      elseif (bcs(1).eq.3) then
        do i=imin,imax
          array(i,0)    = 2*array(i,1) - array(i,2)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at bottom'
        stop
      endif

c Top

      if (bcs(2).eq.0) then
        do i=imin,imax
          array(i,ny+1) = array(i,2)
        enddo
      elseif (bcs(2).eq.1) then
        do i=imin,imax
          array(i,ny+1) = 0d0
        enddo
      elseif (bcs(2).eq.2) then
        do i=imin,imax
cc          array(i,ny+1) = array(i,ny-1)
          array(i,ny+1) = array(i,ny)
        enddo
      elseif (bcs(2).eq.3) then
        do i=imin,imax
          array(i,ny+1) = 2*array(i,ny) - array(i,ny-1)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at top'
        stop
      endif

c Left

      if (bcs(3).eq.0) then
        do j=jmin-1,jmax+1
          array(0   ,j) = array(nx-1,j)
        enddo
      elseif (bcs(3).eq.1) then
        do j=jmin-1,jmax+1
          array(0   ,j) = 0d0
        enddo
      elseif (bcs(3).eq.2) then
        do j=jmin-1,jmax+1
cc          array(0   ,j) = array(2,j)
          array(0   ,j) = array(1,j)
        enddo
      elseif (bcs(3).eq.3) then
        do j=jmin-1,jmax+1
          array(0   ,j) = 2*array(1,j) - array(2,j)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at left'
        stop
      endif

c Right

      if (bcs(4).eq.0) then
        do j=jmin-1,jmax+1
          array(nx+1,j) = array(2,j)
        enddo
      elseif (bcs(4).eq.1) then
        do j=jmin-1,jmax+1
          array(nx+1,j) = 0d0
        enddo
      elseif (bcs(4).eq.2) then
        do j=jmin-1,jmax+1
cc          array(nx+1,j) = array(nx-1,j)
          array(nx+1,j) = array(nx,j)
        enddo
      elseif (bcs(4).eq.3) then
        do j=jmin-1,jmax+1
          array(nx+1,j) = 2*array(nx,j) - array(nx-1,j)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at right'
        stop
      endif

c Resolve singularity if no BC is Dirichlet

cc      if (      bcs(1) /= 1 .and. bcs(2) /= 1 
cc     .    .and. bcs(3) /= 1 .and. bcs(4) /= 1 ) then
cc        array(0,0) = 0d0
cc      endif

c End

      end subroutine

c setBoundaryConditions
c####################################################################
      subroutine setBoundaryConditions2(array,nx,ny,bcs)

c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * array  -> real array with ghost-nodes
c       * nx,ny  -> domain size
c       * bcs    -> integer array of size 4 containing BC setup:
c           + bcs(1) ---> bottom
c           + bcs(2) ---> top
c           + bcs(3) ---> left
c           + bcs(4) ---> right
c         If bcs = 0  --> periodic BC's
c         If bcs = 1  --> dirichlet BC's (always zero because MG deals with
c                                         error updates)
c         If bcs = 2  --> Neumann BC's
c         If bcs = 3  --> second-order derivative = 0 (linear extrapolation)
c--------------------------------------------------------------------

      implicit none       !For safe fortran

c Call variables

      integer*4  nx,ny,bcs(4)
      real*8     array(0:nx+1,0:ny+1)

c Local variables

      integer*4  i,j

c Begin program

c Bottom

      if (bcs(1).eq.0) then
        do i=1,nx
          array(i,0)    = array(i,ny-1)
        enddo
      elseif (bcs(1).eq.1) then
        do i=1,nx
          array(i,0)    = 0d0
        enddo
      elseif (bcs(1).eq.2) then
        do i=1,nx
          array(i,0)    = array(i,2)
cc          array(i,0)    = array(i,1)
        enddo
      elseif (bcs(1).eq.3) then
        do i=1,nx
          array(i,0)    = 2*array(i,1) - array(i,2)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at bottom'
        stop
      endif

c Top

      if (bcs(2).eq.0) then
        do i=1,nx
          array(i,ny+1) = array(i,2)
        enddo
      elseif (bcs(2).eq.1) then
        do i=1,nx
          array(i,ny+1) = 0d0
        enddo
      elseif (bcs(2).eq.2) then
        do i=1,nx
          array(i,ny+1) = array(i,ny-1)
cc          array(i,ny+1) = array(i,ny)
        enddo
      elseif (bcs(2).eq.3) then
        do i=1,nx
          array(i,ny+1) = 2*array(i,ny) - array(i,ny-1)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at top'
        stop
      endif

c Left

      if (bcs(3).eq.0) then
        do j=0,ny+1
          array(0   ,j) = array(nx-1,j)
        enddo
      elseif (bcs(3).eq.1) then
        do j=0,ny+1
          array(0   ,j) = 0d0
        enddo
      elseif (bcs(3).eq.2) then
        do j=0,ny+1
          array(0   ,j) = array(2,j)
cc          array(0   ,j) = array(1,j)
        enddo
      elseif (bcs(3).eq.3) then
        do j=0,ny+1
          array(0   ,j) = 2*array(1,j) - array(2,j)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at left'
        stop
      endif

c Right

      if (bcs(4).eq.0) then
        do j=0,ny+1
          array(nx+1,j) = array(2,j)
        enddo
      elseif (bcs(4).eq.1) then
        do j=0,ny+1
          array(nx+1,j) = 0d0
        enddo
      elseif (bcs(4).eq.2) then
        do j=0,ny+1
          array(nx+1,j) = array(nx-1,j)
cc          array(nx+1,j) = array(nx,j)
        enddo
      elseif (bcs(4).eq.3) then
        do j=0,ny+1
          array(nx+1,j) = 2*array(nx,j) - array(nx-1,j)
        enddo
      else
        write (*,*) 'Boundary condition type not implemented at right'
        stop
      endif


c End

      end subroutine

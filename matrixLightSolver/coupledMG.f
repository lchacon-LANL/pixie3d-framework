c TODO list:
c
c 1) Implement direct block solver for neq > 2.
c
c 2) Clean mg_setup module (including pertaining routines).
c
c 3) Eliminate hardwired boundary conditions (periodic) in Find_mf 
c    routine.

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

      allocate(istart (ngrid),ntotv (ngrid), STAT = alloc_stat)

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
      recursive subroutine mg(ntot,y,x,matvec,options,igrid
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

      integer*4    ntot,igrid,guess,out,depth
      real*8       x(ntot),y(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer          :: iter,neq,igridmin,vcyc
      integer          :: orderres,orderprol,alloc_stat

      integer          :: guess2,outc
      integer          :: izero,i,j,ii,ivcyc

      double precision :: xx(2*ntot),yy(2*ntot),wrk(2*ntot)
      double precision :: rr0,rr1,mag,mag1,mgtol
      double precision :: dummy(ntot),rr(ntot)

      logical          :: fdiag
      logical          :: fpointers

c Begin program

      neq      = options%neq
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
          call find_mf_diag_neq(neq,ntot,matvec,igrid,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Initialize local solution vector (xx) and local r.h.s. (yy)

      xx  = 0d0
      yy  = 0d0
      wrk = 0d0

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
        call matvec(0,ntot,x,dummy,igrid)
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

        call matvec(0,ntot,xx(istart(igrid)),dummy,igrid)

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
cc        if (mag/rr0 < mgtol) exit

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

      options%tol_out = mag1

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

c       Begin program

          nn    = ntotv(igr)
          isig  = istart(igr)
          isigm = istart(igr-1)

c       Relax error/solution on grid number igr/igrid (find new xx)

          call solve(igr)

c       Evaluate residual (ie wrk = yy - A xx = yy - wrk )

          call matvec(0,nn,xx(isig),wrk(isig),igr)

          wrk(isig:isig+nn-1) = yy(isig:isig+nn-1) - wrk(isig:isig+nn-1)

c       Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) down a grid

          call crestrict (neq
     .                   ,yy(isigm),ntotv(igr-1),nxvp(igr-1),nyvp(igr-1)
     .                   ,wrk(isig),ntotv(igr  ),nxvp(igr  ),nyvp(igr  )
     .                   ,orderres,igr,4d0)

          guess2 = 0

c       If on coarsest grid, solve for error, else descend a grid level

          if (igr-1.eq.igridmin) then
            call solve (igr-1)
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
     .                  ,orderprol,igr-1)

c       Update existing error on grid 2 (i.e. xx_2): xx_2 = xx_2 + wrk 

          xx(isig:isig+nn-1) = xx(isig:isig+nn-1) + wrk(isig:isig+nn-1)

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

          call getSolver(nn,yy(isig),xx(isig),matvec,igr,guess2,outc
     .                  ,depth1)

        end subroutine solve

c       killmg
c       ###################################################################
        subroutine killmg

          implicit none

          if (fdiag) deallocate(diag)

          call deallocPointers(fpointers)

        end subroutine killmg

      end subroutine

c jb
c#######################################################################
      recursive subroutine jb(ntot,rr,zz,matvec,options,igrid
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

      integer          :: ntot,igrid,guess,out,depth
      double precision :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer          :: iter,neq,alloc_stat,isig
      double precision :: omega0,omega10,omega01,tol
      logical          :: fdiag
      logical          :: fpointers

      integer*4    i,j,itr,nn,ieq,nx,ny
      integer*4    ii,iii,iip,iim,jjp,jjm
      integer*4    ig,ipg,img,jpg,jmg,iig
      real*8       mag0,mag1,mag,yy(ntot),delta

      double precision,allocatable, dimension(:)  :: dummy,rhs

c Begin program

      neq    = options%neq
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
          call find_mf_diag_neq(neq,ntot,matvec,igrid,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Allocate internal arrays

      allocate(dummy(neq),rhs(neq))

c Jacobi iteration

      nn = ntot

      if (guess.eq.0) then
        do ii=1,nn
          zz(ii) = 0d0
        enddo
      endif

      do itr=1,iter

        call matvec(0,ntot,zz,yy,igrid)

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

              j = 1

              call shiftIndices(i,j,ii,iim,iip,jjm,jjp,nx,ny,1)
              call shiftIndices(i,j,ig,img,ipg,jmg,jpg,nx,ny,isig)

              zz(ii) = zz(ii)+omega0*  (rr(ii) -yy(ii) )/diag(neq,ig )
     .                       +omega10*((rr(iip)-yy(iip))/diag(neq,ipg)
     .                                +(rr(iim)-yy(iim))/diag(neq,img))
     .                       +omega01* (rr(jjp)-yy(jjp))/diag(neq,jpg)

              mag = mag + (rr(ii)-yy(ii))**2

              do j = 2,ny-1

                call shiftIndices(i,j,ii,iim,iip,jjm,jjp,nx,ny,1)
                call shiftIndices(i,j,ig,img,ipg,jmg,jpg,nx,ny,isig)

                zz(ii) =zz(ii)+omega0*  (rr(ii) -yy(ii) )/diag(neq,ig )
     .                        +omega10*((rr(iip)-yy(iip))/diag(neq,ipg)
     .                                 +(rr(iim)-yy(iim))/diag(neq,img))
     .                        +omega01*((rr(jjp)-yy(jjp))/diag(neq,jpg)
     .                                 +(rr(jjm)-yy(jjm))/diag(neq,jmg))

                mag = mag + (rr(ii)-yy(ii))**2
              enddo

              j = ny

              call shiftIndices(i,j,ii,iim,iip,jjm,jjp,nx,ny,1)
              call shiftIndices(i,j,ig,img,ipg,jmg,jpg,nx,ny,isig)

              zz(ii) = zz(ii)+omega0*  (rr(ii) -yy(ii) )/diag(neq,ig )
     .                       +omega10*((rr(iip)-yy(iip))/diag(neq,ipg)
     .                                +(rr(iim)-yy(iim))/diag(neq,img))
     .                       +omega01* (rr(jjm)-yy(jjm))/diag(neq,jmg)

              mag = mag + (rr(ii)-yy(ii))**2

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
cc          if (mag/mag0.lt.tol) exit
        endif

      enddo

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) itr,mag,mag/mag0

      options%iter_out = itr
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag)     deallocate(diag)
      call deallocPointers(fpointers)

      return

 10   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

      end subroutine

c gs
c#######################################################################
      recursive subroutine gs(ntot,rr,zz,matvec,options,igrid
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

      integer          :: ntot,igrid,guess,out,depth
      double precision :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer          :: iter,neq,alloc_stat,isig
      double precision :: omega0,tol
      logical          :: fdiag
      logical          :: fpointers

      integer*4    i,j,itr,nn,ieq,nx,ny
      integer*4    ii,iii,iip,iim,jjp,jjm
      integer*4    ig,ipg,img,jpg,jmg,iig
      real*8       mag0,mag1,mag,yy(ntot),delta

      double precision,allocatable, dimension(:)  :: dummy,rhs

c Begin program

      neq    = options%neq
      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag

      if (out.ge.2) write (*,*)

      allocate(dummy(neq),rhs(neq))

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
          call find_mf_diag_neq(neq,ntot,matvec,igrid,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

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
            call matvec(ii,ntot,zz,yy,igrid)
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ig)
            zz(ii) = zz(ii) + omega0*dummy(neq)
          enddo

c       Backward pass

          do ii = nn,1,-1
            ig = ii + isig - 1
            call matvec(ii,ntot,zz,yy,igrid)
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ig)
            zz(ii) = zz(ii) + omega0*dummy(neq)
            mag = mag + rhs(neq)**2
          enddo

        case (2)

c       Forward pass

          do ii = 1,nn/neq

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              call matvec(iii,ntot,zz,yy,igrid)
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

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              call matvec(iii,ntot,zz,yy,igrid)
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
cc          if (mag/mag0.lt.tol) exit
        endif

      enddo

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) itr,mag,mag/mag0

      options%iter_out = itr
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag) deallocate(diag)
      call deallocPointers(fpointers)

      return

 10   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

      end subroutine

*deck cprolong
c######################################################################
      subroutine cprolong(neq,xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc,order
     .                   ,igc)
c----------------------------------------------------------------------
c     This is a prolongation routine for MG. 
c
c     If order = 0, it employs simple injection.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none            ! For safe Fortran

c Call variables

      integer      neq,ntotc,nxc,nyc,ntotf,nxf,nyf,order,igc
      real*8       xc(ntotc),xf(ntotf)

c Local variables
 
      real*8       xxc(ntotc/neq),xxf(ntotf/neq)

      integer*4    ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

c Injection

      if (order.eq.0) then

      do jc = 1,nyc

        do ic = 1,nxc
          jf = 2*jc
          if = 2*ic

          do ieq = 1,neq
            iic = ieq + neq*(ic-1 + nxc*(jc-1))
            iif = ieq + neq*(if-1 + nxf*(jf-1))

            xf(iif                ) = xc(iic) 
            xf(iif - neq          ) = xc(iic) 
            xf(iif - neq*nxf      ) = xc(iic) 
            xf(iif - neq*nxf - neq) = xc(iic) 
          enddo
        enddo

      enddo

      else

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
     .                ,xxc,nntotc,nxc,nyc,order,igc)

c     Repack prolonged vector

          do i=1,nntotf
            xf(neq*(i-1)+ieq)=xxf(i)
          enddo

        enddo

      endif

c End program

      return
      end subroutine

*deck crestrict
c######################################################################
      subroutine crestrict(neq,xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf,order
     .                   ,igf,volf)
c----------------------------------------------------------------------
c     This is a restriction routine for MG. 
c
c     If order = 0, it employs simple agglomeration.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf,order,igf,neq

      real*8       xc(ntotc),xf(ntotf),volf

c Local variables
 
      real*8       xxc(ntotc/neq),xxf(ntotf/neq)

      integer*4    ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

c Agglomeration

      if (order.eq.0) then

        do jc = 1,nyc

          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic

            do ieq = 1,neq
              iic = ieq + neq*(ic-1 + nxc*(jc-1))
              iif = ieq + neq*(if-1 + nxf*(jf-1))

              xc(iic) = (xf(iif)        +xf(iif        -neq)
     .                 + xf(iif-nxf*neq)+xf(iif-nxf*neq-neq))/4d0*volf
            enddo
          enddo

        enddo

      else

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
     .                 ,xxf,nntotf,nxf,nyf,order,igf,volf)

c     Repack restricted vector

          do i=1,nntotc
            xc(neq*(i-1)+ieq)=xxc(i)
          enddo

        enddo

      endif

c End program

      return
      end subroutine

*deck prolong
c######################################################################
      subroutine prolong(xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc,order,igc)
c----------------------------------------------------------------------
c     This is a prolongation routine for MG. 
c
c     If order = 0, it employs simple injection.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none            ! For safe Fortran

c Call variables

      integer      ntotc,nxc,nyc,ntotf,nxf,nyf,order,igc
      real*8       xc(ntotc),xf(ntotf)

c Local variables
 
      real*8       uc(nxc,nyc)

      real*8       xxf,yyf,ff,xx(nxc),yy(nyc)

      integer*4    ic,if,jc,jf,ii,i,igf,iic,iif

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

        xx(1:nxc) = xl(1:nxc,igc)
        yy(1:nyc) = yl(1:nyc,igc)

        do ic = 1,nxc
          do jc = 1,nyc
            ii = ic + nxc*(jc-1)
            uc(ic,jc) = xc(ii)
          enddo
        enddo

c     Calculate interpolation

        flg = 0
        kx = order+1
        ky = order+1
        nx = nxc
        ny = nyc
        dim = nx*ny + max(2*kx*(nx+1),2*ky*(ny+1))

        allocate(tx(nx+kx))
        allocate(ty(nx+ky))
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
     .                   ,igf,volf)
c----------------------------------------------------------------------
c     This is a restriction routine for MG. 
c
c     If order = 0, it employs simple agglomeration.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf,order,igf

      real*8       xc(ntotc),xf(ntotf),volf

c Local variables
 
      real*8       uf(nxf,nyf)

      real*8       xxc,yyc,xx(nxf),yy(nyf)

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

        xx(1:nxf) = xl(1:nxf,igf)
        yy(1:nyf) = yl(1:nyf,igf)

        do if = 1,nxf
          do jf = 1,nyf
            iif = if + nxf*(jf-1)
            uf(if,jf) = xf(iif)
          enddo
        enddo

c     Calculate interpolation

        flg = 0
        kx = order + 1
        ky = order + 1
        nx = nxf
        ny = nyf
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
      subroutine find_mf_diag_neq(neq,ntot,matvec,ngrid,diag1)
c--------------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec.
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer*4      neq,ntot,ngrid

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

            call findBaseVector(ii,ieq,neq,nn,igrid,x1,1d0)

            call matvec(ii,nn,x1,dummy,igrid)

            call findBaseVector(ii,ieq,neq,nn,igrid,x1,0d0)

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

c findBaseVector
c####################################################################
      subroutine findBaseVector(ii,ieq,neq,ntot,igrid,x1,coef)
c--------------------------------------------------------------------
c     Finds base vector corresponding to grid node ii and equation ieq.
c     Assumes that x1 has been initialized to zero elsewhere.
c--------------------------------------------------------------------

      use mg_setup

      implicit none      !For safe fortran

c Call variables

      integer          :: neq,ii,ieq,ntot,igrid

      double precision :: x1(ntot),coef

c Local variables

      integer          :: nx,ny,iii

c Begin program

      nx  = nxvp(igrid)
      ny  = nyvp(igrid)

      if (mod(ii,nx).eq.1) then
        iii = neq*(ii-1) + ieq
        x1(iii) = coef
        iii = neq*(ii+nx-2) + ieq
        x1(iii) = coef
      elseif (mod(ii,nx).eq.0) then
        iii = neq*(ii-1) + ieq
        x1(iii) = coef
        iii = neq*(ii-nx) + ieq
        x1(iii) = coef
      else
        iii = neq*(ii-1) + ieq
        x1(iii) = coef
      endif

c End program

      return
      end subroutine

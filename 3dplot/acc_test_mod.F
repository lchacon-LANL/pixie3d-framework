c module acc_tests
c ###########################################################################
      module acc_tests

        use parameters

        use variables

        use var_io

        use ts_setup

        logical :: temp_acc=.false.,spat_acc=.false.

        character(200) :: acc_ref_file=''

      contains

c     readAccTestConfig
c     ######################################################################
      function readAccTestConfig(file) result(ioerr)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ioerr
      character(*) :: file

c     Local variables

      namelist /acc_test/ temp_acc,spat_acc,acc_ref_file

c     Read graphics initialization parameters

      open(unit=25,file=trim(file),status='old')
      read(25,acc_test,iostat=ioerr)
      close(unit=25)

c     End program

      end function readAccTestConfig

c     acc_test_space
c     #######################################################################
      subroutine acc_test_space

c     -----------------------------------------------------------------------
c     Tests spatial convergence on variable number icomp.
c     -----------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer        :: i,j,k,ig,jg,kg,nx,ny,nz,ivar
      integer        :: ierr,system,nplot,igx,igy,igz

      real(8)        :: mag,mag0,jac,tref,tsol,lmax
      real(8),allocatable,dimension(:,:,:,:) :: pos 
      real(8),pointer,dimension(:,:,:) :: array

      type (var_array),pointer :: vref => null(),vsol => null()

c     Interpolation

      real(8)    :: xp,yp,zp,interp
      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer    :: kx,ky,kz,nxs,nys,nzs,dim,flg,order,ufile

      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db2val,db3val
      external      db2val,db3val

c     Debug

      character(20) :: debugfile

c     Begin program

      if (my_rank == 0) then
        write (*,*)
        write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write (*,*) 'Starting spatial convergence test...'
        write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      endif

      igx = 1
      igy = 1
      igz = 1

c     Close current record files

      call closeRecordFiles

c     Read reference grid info

      nullify(grid_params)
      call destroyGrid(gv%gparams)

      ufile=openRestartFileForRead(file=acc_ref_file)

      if (my_rank == 0) then
        write (*,*)
        write (*,*) ' Reading reference file...'
        write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd
      endif

      !Create grid
      call createGrid(nxd,nyd,nzd,gv%gparams)

      call setVectorDimensions

      !Read solutions
      call allocateDerivedType(vref)

      do
        call readTimeStep(ufile,itime,time,dt,vref,ierr)
        if (ierr /= 0) exit
      enddo

      tref = time

      call closeRestartFileForRead(ufile)

c Spline reference solution

      nx = ihig-ilog+1
      ny = jhig-jlog+1
      nz = khig-klog+1

      order = 2  !Order of interpolation of ref solution

      nxs = nx+2
      nys = ny+2
      nzs = nz+2

      allocate(xx(nxs),yy(nys),zz(nzs))

      call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)

      xx(1:nxs) = grid_params%xx(ig-1:ig+nx)
      yy(1:nys) = grid_params%yy(jg-1:jg+ny)
      zz(1:nzs) = grid_params%zz(kg-1:kg+nz)

      flg = 0

      kx = min(order+1,nxs-1)
      ky = min(order+1,nys-1)
      kz = min(order+1,nzs-1)

      dim = nxs*nys*nzs + max(2*kx*(nxs+1),2*ky*(nys+1),2*kz*(nzs+1))

      allocate(tx(nxs+kx))
      allocate(ty(nys+ky))
      allocate(tz(nzs+kz))
      allocate(work(dim))
      allocate(bcoef(nxs,nys,nzs))

c Initialize current grid info and read current grid solution

      nullify(grid_params)
      call destroyGrid(gv%gparams)

      igx = 1
      igy = 1
      igz = 1

      ufile=openRestartFileForRead(file=recordfile)

      if (my_rank == 0) then
        write (*,*)
        write (*,*) ' Reading solution file...'
        write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd
      endif

      call createGrid(nxd,nyd,nzd,gv%gparams)

      call setVectorDimensions

      call allocateDerivedType(vsol)

      do
        call readTimeStep(ufile,itime,time,dt,vsol,ierr)
        if (ierr /= 0) exit
      enddo

      tsol = time

c Calculate difference

      if (abs(tref- tsol)>1d-8) then
        if (my_rank == 0) write (*,*) 'T_ref-T_sol=',tref-tsol
        call pstop('acc_test_space'
     .            ,'Reference and solution are not at same time')
      endif

      if (my_rank == 0) write (*,*)

      nx = ihig-ilog+1
      ny = jhig-jlog+1
      nz = khig-klog+1

      do ivar=1,vref%nvar

        !Spline reference component
        array => vref%array_var(ivar)%array

        call db3ink(xx,nxs,yy,nys,zz,nzs,array
     .             ,nxs,nys,kx,ky,kz,tx,ty,tz,bcoef,work,flg)

        !Interpolate to solution grid and compute L2-norm
        array => vsol%array_var(ivar)%array

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              xp = grid_params%xx(ig)
              yp = grid_params%yy(jg)
              zp = grid_params%zz(kg)

              interp = db3val(xp,yp,zp,0,0,0,tx,ty,tz,nxs,nys,nzs
     .                       ,kx,ky,kz,bcoef,work)

              array(i,j,k) = abs(array(i,j,k) - interp)
            enddo
          enddo
        enddo

c diag ****
        write (*,*) 'DIAG -- dumping error file'
        debugfile = 'debug_var='//trim(int2char(ivar))//'.bin'
        open(unit=110,file=debugfile,form='unformatted'
     .      ,status='replace')

        call contour(array(1:nx,1:ny,1),nx,ny
     .              ,0d0,1d0,0d0,1d0,0,110)
        call contour(log10(array(1:nx,1:ny,1)),nx,ny
     .              ,0d0,1d0,0d0,1d0,1,110)

        close(110)
c diag ****

        !L1 norm
        mag = integral(igx,igy,igz,nx,ny,nz,array,average=.true.)

        if (my_rank == 0)
     .    write (*,10) ' L1-norm error, component:',ivar,'=',mag

        !L2 norm
        mag = integral(igx,igy,igz,nx,ny,nz,array**2,average=.true.)
        mag = sqrt(mag)

        if (my_rank == 0)
     .    write (*,10) ' L2-norm error, component:',ivar,'=',mag

        !L-infty norm
        mag = maxval(array(1:nx,1:ny,1:nz))
#if defined(petsc)
        lmax = mag
        call MPI_Allreduce(lmax,mag,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif

        if (my_rank == 0)
     .    write (*,10) ' Linfty-norm error, component:',ivar,'=',mag

        !Number of points with sizable error
        write (*,*) ' Number of points with error > 0.1 Linfty='
     .              ,count((array(1:nx,1:ny,1:nz)>0.1*mag))
        write (*,*) ' Number of points with error > 0.9 Linfty='
     .              ,count((array(1:nx,1:ny,1:nz)>0.9*mag))

      enddo

 10   format(a,i2,a,1pe12.3)

c End program

      deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

      nullify(array)
      call deallocateDerivedType(vref)
      call deallocateDerivedType(vsol)

      end subroutine acc_test_space

c     acc_test_time
c     #######################################################################
      subroutine acc_test_time

c     -----------------------------------------------------------------------
c     Tests temporal convergence on variable number icomp based on variable
c     stored in recordfile.
c     -----------------------------------------------------------------------

      use xdraw_io

      implicit none

c     Call variables

c     Local variables

      integer :: i,j,k,ig,jg,kg,nx,ny,nz,ivar,nxs,nys,nzs
      integer :: ierr,system,nplot,igx,igy,igz,ufile

      real(8) :: mag,mag0,jac,tref,tsol,lmax
      real(8),allocatable,dimension(:,:,:,:) :: pos 
      real(8),pointer,dimension(:,:,:) :: aref,asol

      type (var_array),pointer :: vref => null(),vsol => null()

c     Debug

      character(20) :: debugfile

c     Begin program

      if (my_rank == 0) then
        write (*,*)
        write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write (*,*) 'Starting temporal convergence test...'
        write (*,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      endif

      igx = 1
      igy = 1
      igz = 1

      nxs = nxd
      nys = nyd
      nzs = nzd

      call closeRecordFiles

      time = 0d0

c     Read reference field

      if (my_rank == 0) then
        write (*,*)
        write (*,*) '>>>>> Reading reference file...'
      endif

      !This redefines nxd,nyd,nzd
      ufile=openRestartFileForRead(file=acc_ref_file)

      !Consistency check
      if (nxd /= nxs .or. nyd /= nys .or. nzd /= nzs) then
         call pstop('acc_test_time'
     .   ,'Reference and solution are at different spatial resolutions')
      endif

      call allocateDerivedType(vref)

      do
        call readTimeStep(ufile,itime,time,dt,vref,ierr)
        if (ierr /= 0
     .     .or. (abs(time-tmax)<1d-10.and.tmax > 0d0)
     .     .or. (numtime == itime.and.numtime > 0)) exit
      enddo

      if (my_rank == 0) then
        write (*,*) 'Grid size =',nxd,nyd,nzd
        write (*,*) 'Final time=',time
        write (*,*)
      endif

      call closeRestartFileForRead(ufile)

      tref = time

c     Read current solution

      if (my_rank == 0) then
        write (*,*)
        write (*,*) '>>>>> Reading solution file...'
      endif

      !This redefines nxd,nyd,nzd
      ufile=openRestartFileForRead(file=recordfile)

      call allocateDerivedType(vsol)

      do
        call readTimeStep(ufile,itime,time,dt,vsol,ierr)
        if (ierr /= 0
     .     .or. (abs(time-tmax)<1d-10.and.tmax > 0d0)
     .     .or. (numtime == itime.and.numtime > 0)) exit
      enddo

      if (my_rank == 0) then
        write (*,*) 'Grid size =',nxd,nyd,nzd
        write (*,*) 'Final time=',time
        write (*,*)
      endif

      tsol = time

c     Consistency checks

cc      if (nx /= nxs .or. ny /= nys .or. nz /= nzs) then
cc         call pstop('acc_test_time'
cc     .   ,'Reference and solution are at different spatial resolutions')
cc      endif

      if (abs(tref-tsol)>1d-10) then
        if (my_rank == 0) write (*,*) 'T_ref-T_sol=',tref-tsol
        call pstop('acc_test_time'
     .            ,'Reference and solution are not at same time')
      endif

c     Calculate solution differences

      nx = ihig-ilog+1
      ny = jhig-jlog+1
      nz = khig-klog+1

      if (my_rank == 0) then
        write (*,*) '>>>>> Temporal errors...'
        write (*,*)
      endif

      do ivar=1,vref%nvar
        aref => vref%array_var(ivar)%array
        asol => vsol%array_var(ivar)%array

        asol = abs(asol-aref)

c diag ****
cc        write (*,*) 'DIAG -- acc_test_time'
cc        debugfile = 'debug_var='//trim(int2char(ivar))//'.bin'
cc        open(unit=110,file=debugfile,form='unformatted'
cc     .      ,status='replace')
cc
cc        call contour(asol(1:nx,1:ny,1),nx,ny
cc     .              ,0d0,1d0,0d0,1d0,0,110)
cc
cc        close(110)
c diag ****

        !L1 norm
        mag = integral(igx,igy,igz,nx,ny,nz,asol,average=.true.)

        if (my_rank == 0)
     .    write (*,10) ' L1-norm time error, component:',ivar,'=',mag

        !L2 norm
        mag = integral(igx,igy,igz,nx,ny,nz,asol**2,average=.true.)
        mag = sqrt(mag)

        if (my_rank == 0)
     .    write (*,10) ' L2-norm time error, component:',ivar,'=',mag

        !L-infty norm
        mag = maxval(asol(1:nx,1:ny,1:nz))
#if defined(petsc)
        lmax = mag
        call MPI_Allreduce(lmax,mag,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif

        if (my_rank == 0)
     .   write (*,10) ' Linfty-norm time error, component:',ivar,'=',mag

      enddo

      if (my_rank == 0) write (*,10)

c     End program

 10   format(a,i2,a,1pe12.3)

      nullify(asol,aref)
      call deallocateDerivedType(vref)
      call deallocateDerivedType(vsol)

      end subroutine acc_test_time

      end module acc_tests

c module acc_tests
c ###########################################################################
      module acc_tests

        use parameters

        use variables

        use iosetup

        use timeStepping

        logical :: temp_acc=.false.,spat_acc=.false.

        character(50) :: acc_ref_file

      contains

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

      real(8)        :: mag,mag0,jac,tref,tsol
      real(8),allocatable,dimension(:,:,:,:) :: pos 
      real(8),pointer,dimension(:,:,:) :: array

      type (var_array),pointer :: vref => null(),vsol => null()

c     Interpolation

      real(8)    :: xp,yp,zp,interp
      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer    :: kx,ky,kz,nxs,nys,nzs,dim,flg,order

      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db2val,db3val
      external      db2val,db3val

c     Begin program

      if(my_rank == 0) then
        write (*,*)
        write (*,*) ' Starting spatial convergence test...'
      endif

      igx = 1
      igy = 1
      igz = 1

c     Read reference grid info

      call destroyGrid(grid_params)

      call openRecordFiles(urecord,acc_ref_file,murecord
     .                    ,ilog,ihig,jlog,jhig,klog,khig)

      nx = ihig-ilog+1
      ny = jhig-jlog+1
      nz = khig-klog+1

      if (my_rank == 0) then
        write (*,*)
        write (*,*) ' Reading reference file...'
        write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd
      endif

      call setVectorDimensions

      !Create grid
      call createGrid(nxd,nyd,nzd,gv%gparams)

      !Read solutions
      call allocateDerivedType(vref)

      do
        call readTimeStep(vref,ierr)
cc        ierr = readRecordFile(urecord,itime,tref,dt,vref)
        if (ierr /= 0) exit
      enddo

      tref = time

      call closeRecordFiles

c Spline reference solution

      order = 2

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

      call destroyGrid(grid_params)

      igx = 1
      igy = 1
      igz = 1

      call openRecordFiles(urecord,recordfile,murecord
     .                    ,ilog,ihig,jlog,jhig,klog,khig)

      nx = ihig-ilog+1
      ny = jhig-jlog+1
      nz = khig-klog+1

      if (my_rank == 0) then
        write (*,*)
        write (*,*) ' Reading solution file...'
        write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd
      endif

      call setVectorDimensions

      call createGrid(nxd,nyd,nzd,gv%gparams)

      call allocateDerivedType(vsol)

      do
        call readTimeStep(vsol,ierr)
cc        ierr = readRecordFile(urecord,itime,tsol,dt,vsol)
        if (ierr /= 0) exit
      enddo

      tsol = time

      call closeRecordFiles

c Calculate difference

      if (abs(tref- tsol)>1d-8) then
        if (my_rank == 0) write (*,*) 'T_ref-T_sol=',tref-tsol
        call pstop('acc_test_space'
     .            ,'Reference and solution are not at same time')
      endif

      if (my_rank == 0) write (*,*)

      do ivar=1,vref%nvar

        !Spline component
        array => vref%array_var(ivar)%array

        call db3ink(xx,nxs,yy,nys,zz,nzs,array(0:nx+1,0:ny+1,0:nz+1)
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

              interp = db3val(xp,yp,zp,0,0,0,tx,ty,tz,nx,ny,nz
     .                       ,kx,ky,kz,bcoef,work)

              array(i,j,k) = (array(i,j,k) - interp)**2
            enddo
          enddo
        enddo

        mag = integral(igx,igy,igz,nx,ny,nz,array,average=.true.)

        mag = sqrt(mag)

        !Write output
        write (*,10) ' L2-norm time error, component:',ivar,'=',mag

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
      integer :: ierr,system,nplot,igx,igy,igz

      real(8) :: mag,mag0,jac,tref,tsol
      real(8),allocatable,dimension(:,:,:,:) :: pos 
      real(8),pointer,dimension(:,:,:) :: aref,asol

      type (var_array),pointer :: vref => null(),vsol => null()

c     Debug

      character(13) :: debugfile

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

      call closeRecordFiles

      time = 0d0

c     Read reference field

      if (my_rank == 0) then
        write (*,*)
        write (*,*) '>>>>> Reading reference file...'
      endif

      call openRecordFiles(urecord,acc_ref_file,murecord
     .                    ,ilog,ihig,jlog,jhig,klog,khig)

      call allocateDerivedType(vref)

      do
        call readTimeStep(vref,ierr)
        if (ierr /= 0 .or. abs(time-tmax)<1d-8 ) exit
      enddo

      nx = ihig-ilog+1
      ny = jhig-jlog+1
      nz = khig-klog+1

      if (my_rank == 0) then
        write (*,*) 'Grid size =',nx,ny,nz
        write (*,*) 'Final time=',time
        write (*,*)
      endif

      call closeRecordFiles

      tref = time

c     Read current solution

      if (my_rank == 0) then
        write (*,*)
        write (*,*) '>>>>> Reading solution file...'
      endif

      call openRecordFiles(urecord,recordfile,murecord
     .                    ,ilog,ihig,jlog,jhig,klog,khig)

      call allocateDerivedType(vsol)

      do
        call readTimeStep(vsol,ierr)
        if (ierr /= 0 .or. abs(time-tmax)<1d-10) exit
      enddo

      nxs = ihig-ilog+1
      nys = jhig-jlog+1
      nzs = khig-klog+1

      if (my_rank == 0) then
        write (*,*) 'Grid size =',nxs,nys,nzs
        write (*,*) 'Final time=',time
        write (*,*)
      endif

      call closeRecordFiles

      tsol = time

c     Consistency checks

      if (nx /= nxs .or. ny /= nys .or. nz /= nzs) then
         call pstop('acc_test_time'
     .   ,'Reference and solution are at different spatial resolutions')
      endif

      if (abs(tref-tsol)>1d-10) then
        if (my_rank == 0) write (*,*) 'T_ref-T_sol=',tref-tsol
        call pstop('acc_test_time'
     .            ,'Reference and solution are not at same time')
      endif

c     Calculate difference

      if (my_rank == 0) then
        write (*,*) '>>>>> Temporal errors...'
        write (*,*)
      endif

      do ivar=1,vref%nvar
        aref => vref%array_var(ivar)%array
        asol => vsol%array_var(ivar)%array

        asol = asol-aref
        mag = integral(igx,igy,igz,nx,ny,nz,asol**2,average=.true.)
        mag = sqrt(mag)

        if (my_rank == 0)
     .    write (*,10) ' L2-norm time error, component:',ivar,'=',mag
      enddo

      if (my_rank == 0) write (*,10)

c diag ****
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      call contour(asol(0:nx+1,0:ny+1,1),nx+2,ny+2
cc     .            ,0d0,1d0,0d0,1d0,0,110)
cc
cc      close(110)
c diag ****

c     End program

 10   format(a,i2,a,1pe12.3)

      nullify(asol,aref)
      call deallocateDerivedType(vref)
      call deallocateDerivedType(vsol)

      end subroutine acc_test_time

      end module acc_tests

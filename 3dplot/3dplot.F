      program plot3d

c##########################################################################
c     Postprocess output information from 3DMHD and creates graphics
c     in various formats.
c##########################################################################

      use parameters

      use variables

      use graphics

      use timeStepping

      use grid

      use iosetup

      implicit none

c Local variables

      integer(4)     :: i,j,k,ig,jg,kg
      integer(4)     :: ierr,system,nplot,igx,igy,igz
      character*(40) :: command

c diag
      integer(4)      :: uref=100
c diag

c Begin program

      igx = 1
      igy = 1
      igz = 1

c Initialize 

#if defined(petsc)
      call MPI_Init(mpierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)
#endif

      call initializeCalculation(igx,igy,igz)

c Initialize counters

      nplot     = 0
      tmplot    = 0d0

c Time loop

      do

c     Store previous record

        u_n = u_np

c     Read next record

        call readRecordFile(urecord,itime,time,dt,u_np,ierr)

c     Check for error

        if (ierr == -1) cycle  !Error, but not EOF
        if (ierr == -2) exit   !EOF

c     Update counters

        nplot  = nplot + 1

c     Time level plots (xdraw)

        plot = .false.
        if (nplot.eq.ndplot.or.(time-tmplot).ge.0.999*dplot) then

          call postProcessSolution(u_np,igx,igy,igz)

          call evaluateDiagnostics(u_np,igx,igy,igz,.false.)

          plot = .true.
          nplot  = 0
          if (itime.gt.0) tmplot = time
          call dumpTimeStepPlots
        endif

c     Output per time step

        call output(plot)

      enddo       !End of time loop

c Close graphics files and create draw*.in files

      call finalizeGraphics

      call finalizeDiagnostics

      close(urecord)

c Finalize MPI

#if defined(petsc)
      call MPI_Finalize(mpierr)
#endif

c Second order test

cc      call order_tst

c Move files to 'plot' directory

cc      ierr=system('[ -d plots ] || mkdir plots')
cc      ierr=system('mv draw*in plots')

c Remove backup copy

cc      command='rm '//trim(recordfile)//'.bak'
cc      ierr=system(trim(command))

c End program

      end

c output
c######################################################################
      subroutine output(plot)

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use timeStepping

      use grid

      implicit none

c Call variables

      logical :: plot

c Local variables

      integer(4)  :: ngrd(3)

c Begin program

#if defined(petsc)
      if (my_rank == 0) then
#endif

      ngrd = (/ grid_params%ngrdx,grid_params%ngrdy,grid_params%ngrdz /)

      if (itime.eq.0) then

        write (*,10) nxd,nyd,nzd,ngrd

        write (*,*) 
        write (*,100)

      endif

      if (plot) then
        write (*,120) itime,time,' --> dumped plots'
      else
        write (*,110) itime,time
      endif

#if defined(petsc)
      endif
#endif

c End program

 10   format (/,'  Grid mesh:................',i4,'x',i3,'x',i3
     .        /,'  Number of grid levels:....',i4,',',i2,',',i2)
 100  format ('   itime      time')
 110  format (i7,2x,1p,e10.2)
 120  format (i7,2x,1p,e10.2,a)

      end subroutine output

c initializeCalculation
c######################################################################
      subroutine initializeCalculation(igx,igy,igz)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use iosetup

      use graphics

      use constants

      use generalPurposeFunctions

      implicit none

c Call variables

      integer(4) :: igx,igy,igz,ierr

c Local variables

c Begin program

c Initialize graphics files

#if defined(petsc)
      urecord    = 25 + my_rank

      if (np > 1) then
        recordfile = 'record_proc'//trim(int2char(my_rank))//'.bin'
      else
        recordfile = 'record.bin'
      endif
#else
      urecord    = 2
      recordfile = 'record.bin'
#endif

c Read user initializations

      call readGraphicsInput(sel_diag,sel_graph,ndplot,dplot)

c Open graphics file

#if !defined(petsc)
      call mergeGraphicsFiles
#endif

      call openGraphicsFile(urecord,recordfile)

c Initialize vector dimensions

      call setVectorDimensions

c Initialize MG and create grid

      call createGrid(ilog,ihig,jlog,jhig,klog,khig,nxd,nyd,nzd)
cc      call checkGrid

c Define application arrays (external)

      call allocateApplicationVariables     !External

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Allocate records

      call allocateDerivedType(u_ic)
      call allocateDerivedType(u_0)
      call allocateDerivedType(u_np)
      call allocateDerivedType(u_graph)

c Read equilibrium u_0 (Do NOT postprocess for BC)

      call readRecordFile(urecord,itime,time,dt,u_0,ierr)

      if (ierr /= 0) then
        write (*,*) 'Unable to read equilibrium'
        write (*,*) 'Aborting...'
        stop
      endif

      u_ic = u_0   !Store initial condition (with BCs)
      call postProcessSolution(u_ic,igx,igy,igz)

c Initialize graphics

      u_np = u_ic  !Prime u_np

      call initializeGraphics(igx,igy,igz,bcond)

      call dumpTimeStepPlots

c Initialize diagnostics

      u_n = u_np

      call initializeDiagnostics(u_np,igx,igy,igz)

c End program

      end subroutine initializeCalculation

c openGraphicsFile
c######################################################################
      subroutine openGraphicsFile(unit,file)

c----------------------------------------------------------------------
c     In parallel runs, it merges all output files 'record_proc*.bin'
c     to a serial file 'record.bin', and then opens it for postprocessing.
c     Finally, it reads equilibrium information and stores it in u_0.
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variables

      use graphics

      use iosetup

      use generalPurposeFunctions

      implicit none

c Call variables

      integer(4)   :: unit
      character(*) :: file

c Local variables

      integer(4) :: nx,ny,nz,ierr,nfiles,ifile
      integer(4),allocatable,dimension(:)     :: murecord
      character*(20),allocatable,dimension(:) :: mrecfile

c Externals

      integer(4) :: system

c Begin program

c Open graphics file

      open(unit,file=file,form='unformatted',status='unknown')

c Consistency check

      read (unit) nx,ilog,ihig
      read (unit) ny,jlog,jhig
      read (unit) nz,klog,khig

      nxl  = ihig-ilog+1
      nyl  = jhig-jlog+1
      nzl  = khig-klog+1

#if defined(petsc)
      if ((nx /= nxl .or. ny /= nyl .or. nz /= nzl)
     .    .and.(my_rank == 0)) then
        write (*,*) 'Grid meshes do not agree; cannot restart'
        write (*,*) 'Aborting...'
        stop
      endif
#else
      if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
        write (*,*) 'Grid size does not agree'
        write (*,*) 'Aborting graphics dump...'
        stop
      endif
#endif

c End programs

      end subroutine openGraphicsFile

c mergeGraphicsFile
c######################################################################
      subroutine mergeGraphicsFiles

c----------------------------------------------------------------------
c     In parallel runs, it merges all output files 'record_proc*.bin'
c     to a serial file 'record.bin'.
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variables

      use graphics

      use iosetup

      use generalPurposeFunctions

      implicit none

c Call variables

c Local variables

      integer(4) :: nx,ny,nz,ierr,nfiles,ifile
      integer(4),allocatable,dimension(:)     :: murecord
      character*(20),allocatable,dimension(:) :: mrecfile

c Externals

      integer(4) :: system

c Begin program

      ierr = system('ls record_proc* >& /dev/null')

      if (ierr == 0) then

        write (*,*) 'Merging graphics files...'

        !Find out number of graphics files
        nfiles = -1
        ierr   = 0
        do while (ierr == 0)
          nfiles=nfiles+1
          ierr=
     .        system('ls record_proc'//trim(int2char(nfiles))
     .               //'.bin >& /dev/null')
        enddo

        allocate(mrecfile(nfiles),murecord(nfiles))

        do ifile=1,nfiles
          murecord(ifile) = 20+ifile
          mrecfile(ifile) =
     .         'record_proc'//trim(int2char(ifile-1))//'.bin'
        enddo

        open(urecord,file=recordfile,form='unformatted'
     .       ,status='unknown')

        !Merge graphics files: INITIALIZATION
        nx = 0
        ny = 0
        nz = 0
        do ifile=1,nfiles

          open(murecord(ifile),file=mrecfile(ifile),form='unformatted'
     .        ,status='old')

          read (murecord(ifile)) nxl,ilog,ihig
          read (murecord(ifile)) nyl,jlog,jhig
          read (murecord(ifile)) nzl,klog,khig

          if (nxl /= nxd) then
            nx = nx + nxl
          else
            nx = nxl
          endif
          if (nyl /= nyd) then
            ny = ny + nyl
          else
            ny = nyl
          endif
          if (nzl /= nzd) then
            nz = nz + nzl
          else
            nz = nzl
          endif

        enddo

        write (urecord) nxd,1,nxd
        write (urecord) nyd,1,nyd
        write (urecord) nzd,1,nzd

        !Consistency check

cc        if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
cc          write (*,*) 'Grid size does not agree:',nx,ny,nz
cc          write (*,*) 'Aborting graphics files merging...'
cc          stop
cc        endif

        ihig = nxd
        ilog = 1
        jhig = nyd
        jlog = 1
        khig = nzd
        klog = 1

        call setVectorDimensions

        call allocateDerivedType(u_0)

        do 

          do ifile=1,nfiles

            call readRecordFile(murecord(ifile),itime,time,dt,u_0,ierr)
            if (ierr == -1) cycle !Error, but not EOF
            if (ierr == -2) cycle !EOF

          enddo

          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) exit  !EOF

          !Reset vector dimensions (readRecordFile alters defs of some values)

          ihig = nxd
          ilog = 1
          jhig = nyd
          jlog = 1
          khig = nzd
          klog = 1

          call setVectorDimensions

          call writeRecordFile(urecord,itime,time,dt,u_0)

        enddo

        !Close files
        close(urecord)
        do ifile=1,nfiles
          close(murecord(ifile))
        enddo

        deallocate(mrecfile,murecord)
      endif

      call deallocateDerivedType(u_0)

c End programs

      end subroutine mergeGraphicsFiles

ccc compare
ccc ########################################################################
cc      subroutine compare
cc
cc      use parameters
cc
cc      use variables
cc
cc      use grid
cc
cc      use iosetup
cc
cc      use timeStepping
cc
cc      use equilibrium
cc
cc      use slatec_splines
cc
cc      implicit none
cc
ccc Local variables
cc
cc      integer(4)     :: i,j,k,ig,jg,kg,nx,ny,nz
cc      integer(4)     :: ierr,system,nplot,igx,igy,igz
cc
cc      real(8)        :: mag,mag0,jac
cc      real(8),allocatable,dimension(:,:,:,:) :: pos 
cc      real(8),pointer,dimension(:,:,:) :: array
cc
cc      type (var_array),target :: vref,vsol
cc
ccc Interpolation
cc
cc      real(8)    :: xp,yp,zp,interp
cc      real(8),allocatable,dimension(:) :: xx,yy,zz
cc
cccc      integer(4) :: kx,ky,kz,nx,ny,nz,dim,flg,order
cccc
cccc      real(8), dimension(:),allocatable:: tx,ty,tz,work
cccc      real(8), dimension(:,:,:),allocatable:: bcoef
cccc
cccc      real(8)    :: db2val,db3val
cccc      external      db2val,db3val
cc
ccc Begin program
cc
cc      write (*,*)
cc      write (*,*) ' Starting reference test...'
cc
ccc Read reference grid info
cc
cc      open(unit=urecord,file='record-ref.bin'
cc     .    ,form='unformatted',status='old')
cc
cc      read (urecord) nxd,ilog,ihig
cc      read (urecord) nyd,jlog,jhig
cc      read (urecord) nzd,klog,khig
cc
cc      write (*,*)
cc      write (*,*) ' Reading reference file...'
cc      write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd
cc
cc      !Read solutions
cc      call allocateDerivedType(vref)
cc
cc      do
cc        call readRecordFile(urecord,itime,time,dt,vref,ierr)
cc
cc        if (ierr /= 0) exit
cc      enddo
cc
cc      close (urecord)
cc
cc      array => vref%array_var(IVX)%array
cc
ccc Spline reference solution
cc
cc      order = 2
cc
cc      nx = nxd+2
cc      ny = nyd+2
cc      nz = nzd+2
cc
cc      allocate(xx(nx),yy(ny),zz(nz))
cc
cc      call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)
cc
cc      xx(1:nx) = grid_params%xx(ig-1:ig+nxd)
cc      yy(1:ny) = grid_params%yy(jg-1:jg+nyd)
cc      zz(1:nz) = grid_params%zz(kg-1:kg+nzd)
cc
cc      flg = 0
cc
cc      kx = min(order+1,nx-1)
cc      ky = min(order+1,ny-1)
cc      kz = min(order+1,nz-1)
cc
cc      dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))
cc
cc      allocate(tx(nx+kx))
cc      allocate(ty(ny+ky))
cc      allocate(tz(nz+kz))
cc      allocate(work(dim))
cc      allocate(bcoef(nx,ny,nz))
cc
cc      call db3ink(xx,nx,yy,ny,zz,nz,array(0:nxd+1,0:nyd+1,0:nzd+1)
cc     .           ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)
cc
cc      call deallocateGridMetric(gmetric)
cc      call deallocateGridStructure(grid_params)
cc
ccc Initialize current grid info and read current grid solution
cc
cc      igx = 1
cc      igy = 1
cc      igz = 1
cc
cc      open(urecord,file=recordfile,form='unformatted',status='unknown')
cc
cc      read (urecord) nxd,ilog,ihig
cc      read (urecord) nyd,jlog,jhig
cc      read (urecord) nzd,klog,khig
cc
cc      write (*,*)
cc      write (*,*) ' Reading solution file...'
cc      write (*,*) ' Grid: ',nxd,'x',nyd,'x',nzd
cc
cc      call setVectorDimensions
cc
cc      call createGrid(1,nxd,1,nyd,1,nzd,nxd,nyd,nzd)
cc
cc      call allocateDerivedType(vsol)
cc
cc      do
cc        call readRecordFile(urecord,itime,time,dt,vsol,ierr)
cc
cc        if (ierr /= 0) exit
cc      enddo
cc
cc      close (urecord)
cc
cc      array => vsol%array_var(IVX)%array
cc
ccc Calculate difference
cc
cc      mag0= 0d0
cc      mag = 0d0
cc
cc      do k = 1,nzd
cc        do j = 1,nyd
cc          do i = 1,nxd
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            xp = grid_params%xx(ig)
cc            yp = grid_params%yy(jg)
cc            zp = grid_params%zz(kg)
cc
cc            jac = gmetric%grid(igx)%jac(i,j,k)
cc
cc            interp = db3val(xp,yp,zp,0,0,0,tx,ty,tz,nx,ny,nz
cc     .                     ,kx,ky,kz,bcoef,work)
cc
cccc            mag0 = mag0 + jac* array(i,j,k)**2
cc            mag0 = mag0 + jac* interp**2
cc            mag  = mag  + jac*(array(i,j,k)-interp)**2
cc
cc          enddo
cc        enddo
cc      enddo
cc
cc      mag = sqrt(mag/mag0)
cc
cc      write (*,*)
cc      write (*,*) ' L2-norm relative error:',mag
cc
ccc End program
cc
cc      deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)
cc
cc      end subroutine compare

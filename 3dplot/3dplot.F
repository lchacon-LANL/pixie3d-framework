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

c Backup record file

cc      write (*,*) 'Backing up record file...'
cc
cc      command='cp '//trim(recordfile)//' '//trim(recordfile)//'.bak'
cc
cc      ierr=system(trim(command))
cc
cc      if (ierr /= 0) then
cc        write (*,*) 'Unable to find record file'
cc        write (*,*) 'Aborting...'
cc        stop
cc      else
cc        write (*,*) 'Done!'
cc      endif

c Initialize counters

      nplot     = 0
      tmplot    = 0d0

c Time loop

      do

c     Store previous record

        u_n = u_np

c     Read next record

        call readRecord(urecord,itime,time,dt,u_np,ierr)

c     Check for error

cc        if (ierr /= 0) exit

cc        write (*,*) ierr

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
        write (*,120) itime,time,' dumped plots'
      else
        write (*,110) itime,time
      endif

#if defined(petsc)
      endif
#endif

c End program

 10   format (/,'  Grid mesh:................',i4,'x',i3,'x',i3
     .        /,'  Number of grid levels:....',i4,',',i2,',',i2)
 100  format ('   itime    time')
 110  format (i7,f9.2)
 120  format (i7,f9.2,a)

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

      implicit none

c Call variables

      integer(4) :: igx,igy,igz,ierr

c Local variables

c Begin program

c Read user initializations

      call readGraphicsInput(sel_diag,sel_graph,ndplot,dplot)

c Open graphics file

#if defined(petsc)
      call openGraphicsFile
#endif

c Initialize vector dimensions

#if !defined(petsc)
      ihig = nxd
      ilog = 1
      jhig = nyd
      jlog = 1
      khig = nzd
      klog = 1
#endif

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

c Open graphics file

#if !defined(petsc)
      call openGraphicsFile
#endif

c Read equilibrium u_0 (Do NOT postprocess for BC)

      call readRecord(urecord,itime,time,dt,u_0,ierr)

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
      subroutine openGraphicsFile

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

c Local variables

      integer(4) :: nx,ny,nz,ierr,nfiles,ifile
      integer(4),allocatable,dimension(:)     :: murecord
      character*(20),allocatable,dimension(:) :: mrecfile

c Externals

      integer(4) :: system

c Begin program

c Initialize graphics

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

c For parallel runs, merge graphics files

#if !defined(petsc)
      ierr = system('ls record_proc* >& /dev/null')

      if (ierr == 0) then

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

        do 

          do ifile=1,nfiles

            call readRecord(murecord(ifile),itime,time,dt,u_0,ierr)
            if (ierr == -1) cycle !Error, but not EOF
            if (ierr == -2) cycle !EOF

          enddo

          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) exit  !EOF

          !Reset vector dimensions (readRecord alters defs of some values)

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
#endif

c Open graphics file

      open(urecord,file=recordfile,form='unformatted',status='unknown')

c Consistency check

      read (urecord) nx,ilog,ihig
      read (urecord) ny,jlog,jhig
      read (urecord) nz,klog,khig

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

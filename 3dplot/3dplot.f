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

c Initialize vector dimensions

      call setVectorDimensions

c Initialize MG and create grid

      call createGrid(1,nxd,1,nyd,1,nzd,nxd,nyd,nzd)
cc      call checkGrid

c Define application arrays (external)

      call allocateApplicationVariables     !External

      allocate(zeros (0:nxdp,0:nydp,0:nzdp))
      allocate(vzeros(0:nxdp,0:nydp,0:nzdp,3))
      allocate(ones  (0:nxdp,0:nydp,0:nzdp))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Allocate records

      call allocateDerivedType(u_0)
      call allocateDerivedType(u_np)
      call allocateDerivedType(u_graph)

c Open graphics file and read equilibrium u_0

      call openGraphicsFile

c Postprocess equilibrium

      call postProcessSolution(u_0,igx,igy,igz)

c Initialize graphics

      u_np    = u_0   !Required to "prime" u_np for pointers

      call initializeGraphics(igx,igy,igz,bcond)

      call dumpTimeStepPlots

c Initialize diagnostics

      call initializeDiagnostics(u_0,igx,igy,igz)

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

      use generalOperators

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

      urecord    = 2
      recordfile = 'record.bin'

c For parallel runs, merge graphics files

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

          read (murecord(ifile)) nxl
          read (murecord(ifile)) nyl
          read (murecord(ifile)) nzl

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

        write (urecord) nx
        write (urecord) ny
        write (urecord) nz

        !Consistency check

        if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
          write (*,*) 'Grid size does not agree'
          write (*,*) 'Aborting graphics files merging...'
          stop
        endif

        do 

          do ifile=1,nfiles

            call readRecord(murecord(ifile),itime,time,dt,u_0,ierr)
            if (ierr == -1) cycle !Error, but not EOF
            if (ierr == -2) cycle !EOF

          enddo

          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) exit  !EOF

          !Reset vector dimensions (readRecord alters defs of some values)
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

c Open graphics file

      open(urecord,file=recordfile,form='unformatted',status='unknown')

c Consistency check

      read (urecord) nx
      read (urecord) ny
      read (urecord) nz

      if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
        write (*,*) 'Grid size does not agree'
        write (*,*) 'Aborting graphics dump...'
        stop
      endif

c Read equilibrium

      call readRecord(urecord,itime,time,dt,u_0,ierr)
      if (ierr /= 0) then
        write (*,*) 'Unable to read equilibrium'
        write (*,*) 'Aborting...'
        stop
      endif

c End programs

      end subroutine openGraphicsFile

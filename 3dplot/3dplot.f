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

      integer(4) :: ierr,system,nplot,igx,igy,igz
      real(8)    :: tmplot
      character*(40) :: command

c Begin program

      igx = 1
      igy = 1
      igz = 1

c Initialize 

      call initializeCalculation(igx,igy,igz)

c Backup record file

      command='cp '//trim(recordfile)//' '//trim(recordfile)//'.bak'

      ierr=system(trim(command))

      if (ierr /= 0) then
        write (*,*) 'Unable to find record file'
        write (*,*) 'Aborting...'
        stop
      endif

c Initialize counters

      nplot     = 0
      tmplot    = 0d0

c Time loop

      do

c     Read next record

        call readRecord(urecord,itime,time,u_np,ierr)

c     Check for error

cc        if (ierr /= 0) exit

cc        write (*,*) ierr

        if (ierr == -1) cycle  !Error, but not EOF
        if (ierr == -2) exit   !EOF

c     Update counters (only if timeStep is successful)

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

c Move files to 'plot' directory

cc      ierr=system('[ -d plots ] || mkdir plots')
cc      ierr=system('mv draw*in plots')

c Remove backup copy

      command='rm '//trim(recordfile)//'.bak'
      ierr=system(trim(command))

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

      implicit none

c Call variables

      integer(4) :: igx,igy,igz

c Local variables

      integer(4) :: nx,ny,nz,ierr

c Begin program

      urecord = 2
      recordfile = 'record.bin'

c Read user initializations

      call readGraphicsInput

c Open graphics file

      open(urecord,file=recordfile,form='unformatted',status='unknown')

      read (urecord) nx
      read (urecord) ny
      read (urecord) nz

c Consistency check

      if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
        write (*,*) 'Grid size does not agree'
        write (*,*) 'Aborting graphics dump...'
        stop
      endif

c Initialize vector dimensions

      call setVectorDimensions

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd)

c Define application arrays (external)

      call allocateApplicationVariables     !External

c Allocate records

      call allocateDerivedType(u_0)
      call allocateDerivedType(u_n)
      call allocateDerivedType(u_np)
      call allocateDerivedType(u_graph)

c Read equilibrium

      call readRecord(urecord,itime,time,u_0,ierr)
      if (ierr /= 0) then
        write (*,*) 'Unable to read equilibrium'
        write (*,*) 'Aborting...'
        stop
      endif

      call postProcessSolution(u_0,igx,igy,igz)

c Initialize graphics

      u_graph = u_0

      u_np = u_0   !Required to "prime" u_np for pointers

      call initializeGraphics(igx,igy,igz,bcond)

c Initialize diagnostics

      call initializeDiagnostics(u_0,igx,igy,igz)

c End program

      end subroutine initializeCalculation

c readGraphicsInput
c######################################################################
      subroutine readGraphicsInput

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use graphics

      implicit none

c Call variables

c Local variables

      namelist /graphdef/ sel_diag,sel_graph,ndplot,dplot

c Begin program

c Read computation initializations (external)

      call readInput

c Graphics defaults

      ndplot = 0
      dplot  = 0d0

c Read graphics initialization parameters

      open(unit=25,file='3dmhd.in',status='old')
      read(25,graphdef)
      close(unit=25)

c End program

      end subroutine readGraphicsInput

ccc createGraphics
ccc######################################################################
cc      subroutine createGraphics(igx,igy,igz)
cc
ccc----------------------------------------------------------------------
ccc     Creates graphics pointers, defines dumping intervals
ccc----------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use constants
cc
cc      use timeStepping
cc
cc      use variables
cc
cc      use graphics
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer(4) :: igx,igy,igz
cc
ccc Local variables
cc
ccc Begin program
cc
ccc Initialize graphics
cc
cc
ccc End programs
cc
cc      end subroutine createGraphics

ccc imposeBoundaryConditions
ccc####################################################################
cc      subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
ccc--------------------------------------------------------------------
ccc     Sets adequate boundary conditions on array structure varray.
ccc--------------------------------------------------------------------
cc
cc      use variables
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer(4) :: iigx,iigy,iigz
cc
cc      type (var_array) :: varray
cc
ccc Local variables
cc
ccc Begin program
cc
ccc End program
cc
cc      end subroutine imposeBoundaryConditions

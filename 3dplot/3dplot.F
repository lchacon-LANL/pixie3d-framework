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

      integer        :: i,j,k,ig,jg,kg
      integer        :: ierr,gierr,system,nplot,igx,igy,igz
      character*(40) :: command
      logical        :: exitloop

c Interfaces

      INTERFACE
        subroutine readTimeStep(varray,ierr)
        use variable_setup
        type(var_array),pointer :: varray
        integer         :: ierr
        end subroutine readTimeStep
      END INTERFACE

      INTERFACE
        subroutine postProcessSolution(varray,iigx,iigy,iigz)
        use variable_setup
        integer    :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine postProcessSolution
      END INTERFACE

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use variable_setup
        integer    :: iigx,iigy,iigz
        logical    :: init
        type(var_array),pointer :: varray
        end subroutine evaluateDiagnostics
      END INTERFACE

c Begin program

      igx = 1
      igy = 1
      igz = 1

c Initialize 

      call initializeCalculation(igx,igy,igz)

c Initialize counters

      nplot  = 0
      tmplot = 0d0

c Time loop

      do

c diag ***** Redefine perturbations for preconditioner debugging
cc        write (*,*) 'Redefined perts for PC debugging in 3dplot'
cc        call equateDerivedType(u_ic,u_n)
cc        call imposeBoundaryConditions(u_ic,igx,igy,igz)
c diag *****

c     Read next record

        call readTimeStep(u_np,ierr)

c     Gather error info

#if defined(petsc)
        call MPI_Allreduce(ierr,gierr,1,MPI_INTEGER,MPI_SUM
     .                    ,MPI_COMM_WORLD,mpierr)
        exitloop = (gierr /= 0)
#else
        exitloop = (ierr  /= 0)
#endif
        
c     Check for error

cc        if (ierr == -1) cycle  !Error, but not EOF
cc        if (ierr == -2) exit   !EOF
        if (exitloop) exit

c     Update counters

        nplot  = nplot + 1

c     Time level plots (xdraw)

        plot = .false.

        if (    nplot        == ndplot
     .     .or.(time-tmplot) >= dplot-0.9*dt
     .     .or. time         == 0d0) then

          !Post-process solution for diagnostics
          call postProcessSolution(u_np,igx,igy,igz)

          call evaluateDiagnostics(u_np,igx,igy,igz,.false.)

          !Dump time step plots
          plot  = .true.
          nplot = 0
          if (itime.gt.0) tmplot = time
          call dumpTimeStepPlots

          !Store previous record (for time derivatives in diagnostics)
          call equateDerivedType(u_n,u_np)

        endif

c     Output per time step

        call output(plot)

c     Check for tmax

        if (time.ge.tmax*0.99999)  exit

      enddo       !End of time loop

c Close graphics files and create draw*.in files

      call finalizeGraphics

      call finalizeDiagnostics

      call closeRecordFiles

      if (clean) call closeCleanFiles

#if defined(petsc)
      call PetscFinalize(mpierr)
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

      use parameters

      use grid

      implicit none

c Call variables

      logical :: plot

c Local variables

      integer     :: ngrd(3)

c Begin program

      if (my_rank == 0) then

        ngrd =(/grid_params%ngrdx,grid_params%ngrdy,grid_params%ngrdz/)

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

      endif

c End program

 10   format (/,'  Grid mesh:................',i4,'x',i4,'x',i4
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

      integer    :: igx,igy,igz

c Local variables

      integer    :: ierr
      logical    :: petsc_flg

c Interfaces

      INTERFACE
        subroutine readTimeStep(varray,ierr)
        use variable_setup
        type(var_array),pointer :: varray
        integer         :: ierr
        end subroutine readTimeStep
      END INTERFACE

      INTERFACE
        subroutine postProcessSolution(varray,iigx,iigy,iigz)
        use variable_setup
        integer    :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine postProcessSolution
      END INTERFACE

      INTERFACE
        subroutine initializeDiagnostics(varray,iigx,iigy,iigz)
        use variable_setup
        integer          :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine initializeDiagnostics
      END INTERFACE

c Begin program

c Allocate global variable

      call allocateGlobalVar(gv)

c Initialize PETSc (parallel version)

#if defined(petsc)
      call PetscInitialize(PETSC_NULL_CHARACTER,mpierr)

      !Read PETSc runtime options
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npx',npx
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npy',npy
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npz',npz
     .                       ,petsc_flg,ierr)
#endif

c Read user initializations

      call readGraphicsInput

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd,g_pack,gv%gparams)

c Open graphics files (this initializes local limits, among other things)

      call openGraphicsFiles

c Initialize vector dimensions

      call setVectorDimensions

c Define application arrays (external)

      call allocateApplicationVariables     !External

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Allocate records

      call allocateStructures

c Read equilibrium u_0 (Do NOT postprocess for BC)

      call readTimeStep(u_0,ierr)

      if (ierr /= 0) call pstop('initializeCalculation',
     .                          'Unable to read equilibrium')

c Initialize graphics

      !Set graphics range
      call setGraphicsRange(nxl,nyl,nzl,bcond)

      !Store initial condition w/ BCs and initialize postprocessing quantities
      call equateDerivedType(u_ic,u_0)
      call equateDerivedType(u_n ,u_0)

      call postProcessSolution(u_ic,igx,igy,igz)

      !Prime u_np and initialize graphics quantities
      call equateDerivedType(u_np,u_ic)

      call initializeGraphics(igx,igy,igz,bcond)

      !Dump source for debugging
      if (chk_src) then
        itime = -1
        call dumpTimeStepPlots
      endif

c Initialize diagnostics

      call equateDerivedType(u_n,u_np)

      call initializeDiagnostics(u_np,igx,igy,igz)

c End program

      end subroutine initializeCalculation

c openGraphicsFiles
c######################################################################
      subroutine openGraphicsFiles

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

      integer    :: nx,ny,nz,ierr

c Externals

      integer    :: system

c Begin program

c Open graphics files

#if !defined(petsc)
      ierr = system('ls '//trim(recordfile)//'_proc* >& /dev/null')

      if (ierr == 0) then

        write (*,*)
        write (*,*) 'Opening multiple graphics files...'

        !Find out number of graphics files
        nfiles = -1
        ierr   = 0
        do while (ierr == 0)
          nfiles=nfiles+1
          ierr = system('ls '//trim(recordfile)//'_proc'
     .                 //trim(int2char(nfiles))//' >& /dev/null')
        enddo

        allocate(mrecfile(nfiles),murecord(nfiles))

        do ifile=1,nfiles
          murecord(ifile) = urecord+ifile
          mrecfile(ifile) =
     .        trim(recordfile)//'_proc'//trim(int2char(ifile-1))
          write (*,*) '  => ',trim(mrecfile(ifile))
     $               ,', unit=',murecord(ifile)
        enddo

        if (clean) then
           write (*,*)
           write (*,*) 'Creating clean-up graphics files...'

           allocate(mrecsave(nfiles),musave(nfiles))

           do ifile=1,nfiles
             musave  (ifile) = usave + ifile
             mrecsave(ifile) = trim(mrecfile(ifile))//'_sav'
             write (*,*) '  => ',trim(mrecsave(ifile))
     $                        ,', unit=',musave(ifile)
          enddo
        endif

        write (*,*) 'Finished!'
        write (*,*)

        nx = 0
        ny = 0
        nz = 0
        do ifile=1,nfiles

          open(murecord(ifile),file=trim(mrecfile(ifile))
     $        ,form='unformatted',status='old')

          read (murecord(ifile)) nxl,ilog,ihig
          read (murecord(ifile)) nyl,jlog,jhig
          read (murecord(ifile)) nzl,klog,khig

          if (clean) then
            open(musave(ifile),file=trim(mrecsave(ifile))
     $          ,form='unformatted',status='unknown')

            write (musave(ifile)) nxl,ilog,ihig
            write (musave(ifile)) nyl,jlog,jhig
            write (musave(ifile)) nzl,klog,khig
          endif

        enddo

        !Set global limits (NOT a parallel run)
        ihig = nxd
        ilog = 1
        jhig = nyd
        jlog = 1
        khig = nzd
        klog = 1

      else

        if (clean) recordsave = trim(recordfile)//'_sav'

#else

        urecord = urecord + my_rank

        if (np > 1)
     .     recordfile=trim(recordfile)//'_proc'//trim(int2char(my_rank))

        if (clean) then
          usave = usave + my_rank
          recordsave = trim(recordfile)//'_sav'
        endif

#endif

        open(urecord,file=trim(recordfile),form='unformatted'
     .      ,status='unknown')

        if (clean) open(usave,file=(recordsave),form='unformatted'
     .      ,status='unknown')

c     Consistency check

        read (urecord) nx,ilog,ihig
        read (urecord) ny,jlog,jhig
        read (urecord) nz,klog,khig

        if (clean) then
          write (usave) nx,ilog,ihig
          write (usave) ny,jlog,jhig
          write (usave) nz,klog,khig
        endif

cc        nxl  = grid_params%nlx
cc        nyl  = grid_params%nly
cc        nzl  = grid_params%nlz
        nxl  = ihig-ilog+1
        nyl  = jhig-jlog+1
        nzl  = khig-klog+1

#if defined(petsc)
        if ((nx /= nxl .or. ny /= nyl .or. nz /= nzl)) then
          call pstop('openGraphicsFiles','Grid meshes do not agree')
        endif
#else
        if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
          call pstop('openGraphicsFiles','Grid size does not agree')
        endif

      endif
#endif
        
c End programs

      end subroutine openGraphicsFiles

c closeCleanFiles
c######################################################################
      subroutine closeCleanFiles

c----------------------------------------------------------------------
c     Close graphics files
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

c Begin program

      if (my_rank == 0) then
        write (*,*)
        write (*,*) 'Closing clean-up files...'
        write (*,*)
      endif

#if !defined(petsc)
      if (allocated(musave)) then
        do ifile=1,nfiles
          close(musave(ifile))
        enddo

        deallocate(mrecsave,musave)
      else
        close (usave)
      endif
#else
      close (usave)
#endif
        
c End programs

      end subroutine closeCleanFiles

c closeRecordFiles
c######################################################################
      subroutine closeRecordFiles

c----------------------------------------------------------------------
c     Close graphics files
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

c Begin program

      if (my_rank == 0) then
        write (*,*)
        write (*,*) 'Closing record files...'
        write (*,*)
      endif

#if !defined(petsc)
      if (allocated(murecord)) then
        do ifile=1,nfiles
          close(murecord(ifile))
        enddo

        deallocate(mrecfile,murecord)
      else
        close(urecord)
      endif
#else
      close(urecord)
#endif

c End programs

      end subroutine closeRecordFiles

c readTimeStep
c######################################################################
      subroutine readTimeStep(varray,ierr)

c----------------------------------------------------------------------
c     Reads time step info from single or multiple files.
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variable_setup

      use graphics

      use iosetup

      use generalPurposeFunctions

      implicit none

c Call variables

      type(var_array),pointer :: varray
      integer         :: ierr,system

c Local variables

c Begin program

#if defined(petsc)

      call readRecordFile(urecord,itime,time,dt,varray,ierr)
      if (clean) call writeRecordFile(usave,itime,time,dt,varray)

#else

      if (allocated(murecord)) then
        do ifile=1,nfiles
          call readRecordFile(murecord(ifile),itime,time,dt,varray,ierr)
          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) cycle !EOF
          if (clean) then
cc            call writeRecordFile(musave(ifile),itime,time,dt,varray)
            write (*,*) 'Recycling a parallel run using the serial code'
            write (*,*) 'does not work yet. Need to fix the serial'
            write (*,*) 'write routine to dump only local domain info'
            call closeCleanFiles
            write (*,*) 'Removing clean-up files...'
            ierr = system('rm *_sav >& /dev/null')
            call pstop('readTimeStep','')
          endif
        enddo
      else
        call readRecordFile(urecord,itime,time,dt,varray,ierr)
        if (clean) call writeRecordFile(usave,itime,time,dt,varray)
      endif

#endif

c End program

      end subroutine readTimeStep

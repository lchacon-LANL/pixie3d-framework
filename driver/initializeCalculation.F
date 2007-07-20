c initializeCalculation
c######################################################################
      subroutine initializeCalculation

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      implicit none

c Call variables

c Local variables

      integer    :: ierr,system

      logical    :: petsc_flg

c Interfaces

      INTERFACE
         subroutine setInitialCondition(varrayn,varraynp)
         use variable_setup
         type(var_array),pointer :: varrayn,varraynp
         end subroutine setInitialCondition
      END INTERFACE

      INTERFACE
        subroutine evaluateNonlinearFunction(varray,fi)
        use parameters
        use variable_setup
        real(8)          :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluateNonlinearFunction
      END INTERFACE

c Begin program

c Assign pointers

      g_pack%dim(:)%pack = .false.

c Read user initializations

      call readInput

c Read PETSc runtime options

#if defined(petsc) && !defined(petsc_c)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npx',npx
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npy',npy
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npz',npz
     .                       ,petsc_flg,ierr)
#endif

c Check for autoinitializing parameters

      pi = acos(-1d0)

      if (maxitnwt.eq.0) 
     .      maxitnwt = max(floor(1.5*log(rtol)/log(tolgm)),10)

      alpha = 1. - cnfactor

      dtbase = dt

      !Time smoothing
      if (sm_flag == 0) then
        sm_pass= 0
        if (cnfactor > 1d0 .or. cnfactor < 0d0) then
          call pstop('initializeCalculation'
     .              ,'cnfactor out of range')
        endif
      endif

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd,g_pack,gv%gparams)

c Initialize global domain limits

      call setVectorDimensions

c Allocate constant arrays

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Create nonlinear solver

      call createNonlinearSolver

c Create nonlinear function

      call createNonlinearFunction

c Create equilibrium u_0

      call createEquilibrium

c Initialize old time solution

      call equateDerivedType(u_n,u_0)

c Set unperturbed forcing fields

      jit = -1   !This informs nlfunction that we are processing 
                 !  n time level info

      !This not only evaluates fsrc, but defines BCs on u_n
      call evaluateNonlinearFunction(u_n,fsrc)

      if (.not.source) fsrc = 0d0

c Set output files

#if defined(petsc) 
      if (.not.restart) then
        if (my_rank == 0)
     .     ierr=system('rm -f '//trim(recordfile)//'* >& /dev/null')
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      endif

      urecord = urecord + my_rank

      if (np > 1) then
        recordfile=trim(recordfile)//'_proc'//trim(int2char(my_rank))
      endif
#else
      if (.not.restart) then
         ierr=system('rm -f '//trim(recordfile)//'* >& /dev/null')
      endif
#endif

c Set initial condition

      call setInitialCondition(u_n,u_np)

c Initialize record file

      call initializeRecordFile

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        if (my_rank == 0) then
          write(*,*)
          write(*,*) 'Tmax is less or equal than restarting time'
          write(*,*) 'Aborting'
        endif
#if defined(petsc)
        call PetscFinalize(mpierr)
#endif
        stop
      endif

c End program

      end subroutine initializeCalculation

c createNonlinearSolver
c######################################################################
      subroutine createNonlinearSolver

c----------------------------------------------------------------------
c     Allocates nonlinear solver arrays
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variables

      implicit none

c Call variables

c Local variables

c Begin program

      call allocateStructures

cc      u_0 => global_var%u_0

      allocate(fold(ntotd),fsrc(ntotd))

      allocate(cnf(neqd),one_over_dt(neqd)
     .        ,bdfp(neqd),bdfn(neqd),bdfnm(neqd))

c Initialization

      bdfp  =  1d0
      bdfn  = -1d0
      bdfnm =  0d0

c End programs

      end subroutine createNonlinearSolver

c createEquilibrium
c######################################################################
      subroutine createEquilibrium

c----------------------------------------------------------------------
c     Creates equilibrium according to user input.
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variable_setup

      use grid

      implicit none

c Call variables

c Local variables

      integer    :: ieq

      real(8),allocatable,dimension(:,:,:,:) :: var

      character*(20),allocatable,dimension(:):: label

      integer   ,allocatable,dimension(:,:)  :: bcs

c Begin program

c Set equilibrium u_0 and define BCs on all variables

      allocate(var(ilom:ihip,jlom:jhip,klom:khip,neqd)
     $            ,label(neqd),bcs(6,neqd))

      !Initialize boundary conditions
      do ieq = 1,neqd
        bcs(:,ieq) = bcond
      enddo

      var = 0d0

      call setEquilibrium(1,1,1,var,bcs,label)

      call packVariables(u_0)

      deallocate(var,label,bcs)

c End programs

      contains

c     packVariables
c     #################################################################
      subroutine packVariables(varray)

c     -----------------------------------------------------------------
c     Packs arrays into variable structure 
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer    :: ieq

c     Begin program

      varray%nvar = neqd

      do ieq = 1,neqd
        call VarPack (var(:,:,:,ieq),bcs(:,ieq),label(ieq),ieq,varray)
      enddo

c     End program

      end subroutine PackVariables

      end subroutine createEquilibrium

c setInitialCondition
c####################################################################
      subroutine setInitialCondition(varrayn,varraynp)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation. Variables
c     are stored in u_np
c--------------------------------------------------------------------

      use icond

      use grid

      use variable_setup

      use timeStepping

      use constants

      use iosetup

      implicit none

c Call variables

      type(var_array),pointer :: varrayn,varraynp

c Local variables

      integer    ::  ieq,nx,ny,nz

c Begin program

c Perturb equilibrium

      if (.not.restart) then

        do ieq = 1,neqd
          call perturbEquilibrium(varrayn%array_var(ieq)%array
     .                           ,varrayn%array_var(ieq)%bconds
     .                           ,pert(ieq),ieq)
        enddo

        time     = 0d0
        inewtime = 1

        varraynp = varrayn

      else

        call readRestartFile (nx,ny,nz,itime,time,varrayn,varraynp)

        inewtime = itime+1

      endif

c End program

      contains

c     readRestartFile
c     #################################################################
      subroutine readRestartFile(nx,ny,nz,itime,time,vn,vnp)

c     -----------------------------------------------------------------
c     Reads restart file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer   ,intent(OUT) :: nx,ny,nz,itime
      real(8),intent(OUT)    :: time

      type(var_array),pointer       :: vn,vnp

c     Local variables

      integer    :: ierr

      type(var_array),pointer       :: vmed

c     Begin program

      call allocateDerivedType(vmed)

      open(unit=urecord,file=recordfile,form='unformatted',status='old')

#if defined(petsc) 
      read (urecord) nx,ilog,ihig
      read (urecord) ny,jlog,jhig
      read (urecord) nz,klog,khig

      if ((nx /= nxl .or. ny /= nyl .or. nz /= nzl)) then
#else
      read (urecord) nx
      read (urecord) ny
      read (urecord) nz

      if (nx /= nxd.or.ny /= nyd.or.nz /= nzd) then
#endif
        call pstop('initializeCalculation'
     $            ,'Grid meshes do not agree; cannot restart')
      endif

      if (my_rank == 0) write (*,*) 'Reading restart file(s)...'

      call readRecordFile(urecord,itime,time,dt,vn,ierr)

      vnp = vn

      do
        call readRecordFile(urecord,itime,time,dt,vmed,ierr)

        if (ierr /= 0) then
          exit
        else
          vn  = vnp
          vnp = vmed
        endif
      enddo

      close (urecord)

      if (my_rank == 0) write (*,*) 'Done!'

c     End

      call deallocateDerivedType(vmed)

      end subroutine readRestartFile

      end subroutine setInitialCondition

c initializeRecordFile
c######################################################################
      subroutine initializeRecordFile

c----------------------------------------------------------------------
c     Creates graphics pointers, defines dumping intervals
c----------------------------------------------------------------------

      use timeStepping

      use variables

      use iosetup

      implicit none

c Call variables

c Local variables

      integer    :: ierr

      integer    :: system
      external      system

c Begin program

c Set data dumping intervals

      dfreq = 8d0

      if (tmax.gt.0d0) then
        if (dstep.eq.0d0) dstep = dt*max(int((tmax-time)/dfreq/dt),1)
        rstep = min(dt*max(int((tmax-time)/dfreq/dt),1),dstep)
        ndstep = -1
        numtime = -1
      else
        if (ndstep.eq.0) ndstep = max(numtime/int(dfreq),1)
        nrstep = min(max(numtime/int(dfreq),1),ndstep)
        dstep = 1e30
        tmax  = 0d0
      endif

c Open record file

      if (.not.restart) then
        !Initially dump u_n instead of u_0 (w/BCs) for comparison w/ preconditioner solution
cc        if (debug) then
cc          write (*,*) 'DIAG -- Dumping graphs for PC testing'
cc          u_graph = u_n
cc        else
        !Impose BC's on u_graph <- u_0 (do not overwrite u_0, since it contains equil. BCs)
           u_graph = u_0
cc        endif

        !Check source (need to uncomment dumpTimeStepPlots call in 3dplot.F)
         if (chk_src) then
cc           u_graph = fsrc
           call mapVectorToStructure(u_graph,fsrc)
cc           u_graph%array_var(1)%array = u_0%array_var(1)%array !Set rho finite
cc           call imposeBoundaryConditions(u_graph,1,1,1)
         endif

        !Open record file
        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='unknown')

#if defined(petsc)
        write (urecord) nxl,ilog,ihig
        write (urecord) nyl,jlog,jhig
        write (urecord) nzl,klog,khig
#else
        write (urecord) nxd,1,nxd
        write (urecord) nyd,1,nyd
        write (urecord) nzd,1,nzd
#endif

cc        call writeRecordFile(urecord,0,0d0,dt,u_n)
        call writeRecordFile(urecord,0,0d0,dt,u_graph)
        call writeRecordFile(urecord,0,0d0,dt,u_np)

      else

        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='old',position='append')
          
      endif

c End programs

      end subroutine initializeRecordFile

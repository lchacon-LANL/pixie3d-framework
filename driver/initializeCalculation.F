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

      use iosetup

      implicit none

c Call variables

c Local variables

      integer    :: ierr

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

      INTERFACE
         subroutine setEquilibrium(igx,igy,igz,varray)
         use variable_setup
         integer :: igx,igy,igz
         type(var_array),pointer :: varray
         end subroutine setEquilibrium
      END INTERFACE

c Begin program

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
      call PetscOptionsGetTruth(PETSC_NULL_CHARACTER,'-test',tst_flg
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

      itime = 0

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd,g_pack,gv%gparams)

c Initialize global domain limits

      call setVectorDimensions

#if defined(petsc)
      if (.not.tst_flg)
     .   write (*,'(a,i4,a,6i5)') 'Processor:',my_rank
     .           ,'; Grid limits:',ilog,ihig,jlog,jhig,klog,khig
#endif

c Allocate structures

      call allocateStructures

c Create nonlinear solver

      call createTS

c Create nonlinear function

      call createNonlinearFunction    !External

c Create equilibrium u_0

      call setEquilibrium(1,1,1,u_0)  !External

c Initialize old time solution

      call equateDerivedType(u_n,u_0)

c Set unperturbed forcing fields

      source_eval = .true.   !This informs nlfunction that we are processing the source
      
      !This not only evaluates fsrc, but defines BCs on u_n
      call evaluateNonlinearFunction(u_n,fsrc)

      !Dump source
      if (chk_src) call dump_src

      if (.not.source) fsrc = 0d0

      source_eval = .false.

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

        call fortranDestroy
        stop
      endif

c End program

      contains

c     dump_src
c     ###############################################################
      subroutine dump_src

        use grid_debug,ONLY:C_contour,cont_conf,dbg,createDrawInCfile

        integer :: ieq,nx,ny,nz,nunit

        character(10) :: srcf
        character(4)  :: dbg_desc(neqd)
        character(3)  :: xlabel,ylabel

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '>>>>>>>>>>>>>>>>>'
          write (*,*) 'Dumping source...'
          write (*,*) '>>>>>>>>>>>>>>>>>'
        endif

        nx = grid_params%nxv(1)
        ny = grid_params%nyv(1)
        nz = grid_params%nzv(1)

        xlabel=cont_conf%label(1)
        ylabel=cont_conf%label(2)

        do ieq=1,neqd
          dbg_desc(ieq) = 'dgb'//trim(int2char(ieq))
        enddo

        srcf='src_p'//trim(int2char(my_rank))//'.bin'

        call createDrawInCfile(neqd,trim(srcf),'Source'
     .        ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
     .        ,'drawsrc_p'//trim(int2char(my_rank))//'.in')

        call mapVectorToStructure(u_graph,fsrc)

        nunit = 110

        open(nunit,file=trim(srcf),form='unformatted'
     .        ,status='replace')

        do ieq=1,neqd
          call C_contour(u_graph%array_var(ieq)%array(1:nx,1:ny,1:nz)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,ieq-1,nunit)
        enddo

        close (nunit)

      end subroutine dump_src

      end subroutine initializeCalculation

c setInitialCondition
c####################################################################
      subroutine setInitialCondition(varrayn,varraynp)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation. Variables
c     are stored in u_np
c--------------------------------------------------------------------

      use grid

      use variable_setup

      use timeStepping

      use iosetup

      implicit none

c Call variables

      type(var_array),pointer :: varrayn,varraynp

c Local variables

      integer    ::  ieq,nx,ny,nz

c Interfaces

      INTERFACE
        subroutine applyBC (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine applyBC
      END INTERFACE

c Begin program

c Perturb equilibrium

      if (.not.restart) then

        !Equate to varraynp
        call equateDerivedType(varraynp,varrayn)

        do ieq = 1,neqd
          call perturbEquilibrium(varraynp%array_var(ieq)%array
     .                           ,varraynp%array_var(ieq)%bconds
     .                           ,ieq)
        enddo

        time     = 0d0
        itime    = 1   !Needed for BC logic
        inewtime = 1

      else

        call readRestartFile(nx,ny,nz,itime,time,varrayn,varraynp)

        inewtime = itime+1

      endif

c Impose BC on varrayn

      call applyBC(varraynp,1,1,1)

c End program

      contains

#if !defined(adios)
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

      type(var_array),pointer :: vn,vnp

c     Local variables

      integer    :: ierr

      type(var_array),pointer :: vmed => null()

c     Begin program

c     Open restart files and check mesh sizes

      call openBinFile(urecord,recordfile,status='old')

      read (urecord) nx,ilog,ihig
      read (urecord) ny,jlog,jhig
      read (urecord) nz,klog,khig

#if defined(petsc) 
      if (nx /= nxl .or. ny /= nyl .or. nz /= nzl) then
#else
      if (nx /= nxd .or .ny /= nyd .or. nz /= nzd) then
#endif
        call pstop('initializeCalculation'
     $            ,'Grid meshes do not agree; cannot restart')
      endif

c     Allocate intermediate variable

      call allocateDerivedType(vmed)

c     Read files

      if (my_rank == 0) write (*,*) ' Reading restart file(s)...'

      ierr = readRecordFile(urecord,itime,time,dt,vn)

      call equateDerivedType(vnp,vn)

      do
        ierr = readRecordFile(urecord,itime,time,dt,vmed)

        if (ierr /= 0) then
          exit
        else
          call equateDerivedType(vn ,vnp)
          call equateDerivedType(vnp,vmed)
        endif
      enddo

      if (my_rank == 0) write (*,*) ' Done!'

c     End

      close(urecord)

      call deallocateDerivedType(vmed)

      end subroutine readRestartFile

#else

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

      type(var_array),pointer :: vn,vnp

c     Local variables

      integer    :: ierr

      type(var_array),pointer :: vmed => null()

c     Begin program

      call allocateDerivedType(vmed)

      call adios_init_noxml (adios_err)

      ! allocate 100MB buffer for ADIOS
      call adios_allocate_buffer (100,adios_err)

      if (my_rank == 0) write (*,*) ' Reading restart file(s)...'

      call openRecordFileAdiosForRead(adios_err)

      if (adios_err /= 0) then
        call pstop('readRestartFile'
     $            ,'Error reading ADIOS restart file')
      endif

      ierr = readRecordFile(urecord,itime,time,dt,vn)

      call equateDerivedType(vnp,vn)

      do
        ierr = readRecordFile(urecord,itime,time,dt,vmed)

        if (ierr /= 0) then
          exit
        else
          call equateDerivedType(vn ,vnp)
          call equateDerivedType(vnp,vmed)
        endif
      enddo

      call closeRecordFileAdiosForRead(adios_err)

      if (my_rank == 0) write (*,*) ' Done!'

c     End

      call deallocateDerivedType(vmed)

      end subroutine readRestartFile
#endif

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

      integer    :: ierr,ieq

c Interfaces

      INTERFACE
        subroutine applyBC (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine applyBC
      END INTERFACE

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

        !PC debugging
cc        if (chk_src) then
cc          call mapVectorToStructure(u_graph,fsrc)
cc        else
cccc         if (debug) then
cccc           write (*,*) 'DIAG -- Dumping graphs for PC testing'
cccc           u_graph = u_n
cccc         endif
cc        endif

        !Impose BC's on u_graph <- u_0
        !(do not overwrite u_0, since it contains equil. BCs;
        ! do nothing at EQU boundaries, to preserve EQU BCs)
        call equateDerivedType(u_graph,u_0)

cc        call applyBC(u_graph,1,1,1,obc=0)
        call applyBC(u_graph,1,1,1)

#if !defined(adios)
        !Initialize POSIX files
        call init_POSIX
#endif

        !Dump initial data
        call writeRecordFile(urecord,0,0d0,dt,u_graph)
        call writeRecordFile(urecord,0,0d0,dt,u_np)

      else

#if defined(adios)
        adios_open_num = 1  !Set ADIOS file to append
#else
        call openBinFile(urecord,recordfile,status='old'
     .                  ,position='append'
     .                  ,check=.true.)
#endif
      endif

c End programs

      contains

c     init_POSIX
c     ##############################################################
      subroutine init_POSIX

        if (my_rank == 0) ierr = rm_files(recordfile)

#if defined(petsc) 
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

        !Initialize record file
        call openBinFile(urecord,recordfile,status='unknown')

        write (urecord) nxl,ilog,ihig
        write (urecord) nyl,jlog,jhig
        write (urecord) nzl,klog,khig

      end subroutine init_POSIX

      end subroutine initializeRecordFile

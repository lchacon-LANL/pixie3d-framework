c initializeCalculation
c######################################################################
      subroutine initializeCalculation

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use timeStepping

      use nk

      use var_io

      use mg_solver

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr

      integer :: mgpcf,mgcs

      logical :: petsc_flg
     .          ,bc_debg=.false.
     .          ,mg_coarse=.true.

c Begin program
     
c Read user initializations

      call readInput

c Account for MK equation

#if defined(MK)
      if (mk_grid) neqd = neqd+1
#endif

c Read PETSc runtime options

#if defined(petsc) && !defined(petsc_c)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nxd',nxd
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nyd',nyd
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nzd',nzd
     .                       ,petsc_flg,ierr)

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npx',npx
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npy',npy
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npz',npz
     .                       ,petsc_flg,ierr)

      call PetscOptionsGetTruth(PETSC_NULL_CHARACTER,'-pc_debug'
     .                         ,pc_debug,petsc_flg,ierr)

      call PetscOptionsGetTruth(PETSC_NULL_CHARACTER,'-bc_debug'
     .                         ,bc_debg,petsc_flg,ierr)
      if (petsc_flg) call set_BC_debug(bc_debg)

      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-rtol',rtol
     .                        ,petsc_flg,ierr)

#if defined(coarse_MG)
      call PetscOptionsGetTruth(PETSC_NULL_CHARACTER,'-mg_coarse'
     .                         ,mg_coarse,petsc_flg,ierr)
      if (petsc_flg) call MG_set_coarse_solve(mg_coarse)

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-mgpcf',mgpcf
     .                       ,petsc_flg,ierr)
      if (petsc_flg) call MG_set_proc_coarsen_factor(mgpcf)

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-mgcs',mgcs
     .                       ,petsc_flg,ierr)
      if (petsc_flg) call MG_set_coarse_size(mgcs)
#endif

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-ilevel',ilevel
     .                       ,petsc_flg,ierr)

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-numtime',numtime
     .                       ,petsc_flg,ierr)
      if (petsc_flg) tmax = 0d0

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-ndstep',ndstep
     .                       ,petsc_flg,ierr)
      if (petsc_flg) dstep = 0d0

      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-tmax',tmax
     .                        ,petsc_flg,ierr)
      if (petsc_flg) numtime = 0

      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-dstep',dstep
     .                        ,petsc_flg,ierr)
      if (petsc_flg) ndstep = 0

      call PetscOptionsGetTruth(PETSC_NULL_CHARACTER,'-restart'
     .                         ,ts_restart,petsc_flg,ierr)

      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-trestart',trestart
     .                        ,petsc_flg,ierr)
      if (petsc_flg) ts_restart = .true.

      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-cpu_quota',quota
     .                        ,petsc_flg,ierr)

      call PetscOptionsGetTruth(PETSC_NULL_CHARACTER,'-test'
     .                         ,tst_flg,petsc_flg,ierr)

#else
!!$  n_args = COMMAND_ARGUMENT_COUNT()
!!$
!!$  if (n_args < 4) then
!!$     call pstop('parareal_advance','Wrong call sequence')
!!$  endif
!!$
!!$  if (debug) write (*,*) "Reading command arguments"
!!$
!!$  call get_command_argument(1, ifile1)
!!$  call get_command_argument(2, ifile2)
!!$  call get_command_argument(3, ifile3)
!!$  call get_command_argument(4, ofile)
#endif

c Allocate global variables

      call allocateGlobalVar(gv)

#if defined(petsc)
c$$$        if (.not.tst_flg)
c$$$     .       write (*,'(a,i4,a,6i5)') 'Processor:',my_rank
c$$$     .           ,'; Grid limits:',ilog,ihig,jlog,jhig,klog,khig
#endif

c Create Time Stepping variables

      call TS_create

c Create nonlinear function

      call createNLF    !External

c Initialize MK variables

#if defined(MK)
      if (mk_grid) call MK_alloc

      relax_init_grid = mk_grid.and.mk_init_relax.and.(.not.ts_restart)

      if (relax_init_grid) init_grid_offset = 1
#endif
     
c Initialize computation

      if (.not.ts_restart) then
        call setInitialState
      else
        call readCurrentState(recordfile)
      endif

c Initialize record file

      call initializeRecordFile(recordfile)

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        call pstop("initializeCalculation"
     $            ,"Tmax is less or equal than restarting time")
      endif

c End program

      end subroutine initializeCalculation

c setInitialState
c####################################################################
      subroutine setInitialState

c--------------------------------------------------------------------
c     Set initial state for initial value calculation.
c--------------------------------------------------------------------

      use grid

      use timeStepping

      use var_io

#if defined(MK)
      use mk
#endif

      implicit none

c Local variables

      integer :: ieq,nx,ny,nz,igrid

c Interfaces

      INTERFACE
        subroutine setEquilibrium(igx,varray)
        use var_setup
        integer :: igx
        type(var_array),pointer :: varray
        end subroutine setEquilibrium
      END INTERFACE

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

      INTERFACE
        subroutine get_src(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine get_src
      END INTERFACE

c Begin program

      igrid = 1

c Create equilibrium u_0

      call setEquilibrium(igrid,u_0)  !External

#if defined(MK)
      if (mk_grid) then
        u_0%array_var(neqd)%descr = 'MK phi'
        u_0%array_var(neqd)%array = 0d0

        call MK_initGrid(igrid,bcs=bcond)  !Change also in 3dplot.F
      endif
#endif

      call equateDerivedType(u_n,u_0)

      call get_src(u_n)  !This also imposes BCs

      call mapStructureToVector(gv%gparams,x0,u_n)

c Write geometry when applicable

cc      if (coords == 'ext') call writeExtMap(gv%gparams)

c Perturb equilibrium

      call equateDerivedType(u_np,u_n)

      if (.not.relax_init_grid) call perturbEquilibrium(u_np)

c Impose BC (to initialize all auxiliary variables)

      call applyBC(igrid,u_np,gv%aux,it=1)

#if defined(MK)
      if (mk_grid) call MK_applyBC(igrid,u_np%array_var(neqd)%array)
#endif

c End program

      end subroutine setInitialState

c readCurrentState
c####################################################################
      subroutine readCurrentState(file)

c--------------------------------------------------------------------
c     Read current state from restart file.
c--------------------------------------------------------------------

      use grid

      use timeStepping

      use var_io

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      character(*) :: file
      
c Local variables

      integer :: ieq,igrid

c Interfaces

      INTERFACE
        subroutine get_src(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine get_src
      END INTERFACE

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

c Begin program

      igrid = 1

c Read geometry when applicable

cc      if (coords == 'ext') call readExtMap(gv%gparams)

c Read restart file

      call readRestartFile(itime,time)

#if defined(MK)
      if(mk_grid) call MK_setup_grid(igrid,u_np%array_var(neqd)%array)

      if(mk_grid.and.mk_init_relax.and.itime==1) itime = 0
#endif

c Introduce perturbation at restart

      inewtime = itime+1

      if (perturb_rstrt) call perturbEquilibrium(u_np)

c Impose BCs (to initialize all auxiliary variables)

      call applyBC(igrid,u_np,gv%aux,it=itime)

#if defined(MK)
      if (mk_grid) call MK_applyBC(igrid,u_np%array_var(neqd)%array)
#endif

c End program

      contains

c     readRestartFile
c     #################################################################
      subroutine readRestartFile(itime,time)

c     -----------------------------------------------------------------
c     Reads restart file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer,intent(OUT) :: itime
      real(8),intent(OUT) :: time

c     Local variables

      integer :: nx,ny,nz,ierr,ufile,itime_old

      type(var_array),pointer :: vmed => null()

c     Begin program

c     Open restart files

      ufile = openRestartFileForRead(file=file)

c     Allocate intermediate variable

      call allocateDerivedType(vmed)

c     Read equilibrium

      if (my_rank == 0) write (*,*) ' Reading restart file(s)...'

      ierr = readRecordFile(ufile,itime,time,dt,gammat,u_0)

      if (ierr /= 0) then
        call pstop('readCurrentState','Unable to read u_0')
      endif

      call TS_update_preCounters(time,exclude_itime=.true.)

c     Read files

      ierr = readRecordFile(ufile,itime,time,dt,gammat,u_n)

      if (ierr /= 0) then
        call pstop('readCurrentState','Unable to read u_n')
      endif

      call equateDerivedType(u_np,u_n)

      do
        itime_old = itime + 1

        ierr = readRecordFile(ufile,itime,time,dt,gammat,vmed)

        if(switch_dumpfile(itime_old,itime,ufile)) then
          ierr = readRecordFile(ufile,itime,time,dt,gammat,vmed)
        endif

        if (ierr /= 0) then !EOF: Rollback counters and exit
          call TS_rollback_counters(exclude_itime=.true.)
          exit
        else
          call TS_update_preCounters(time,exclude_itime=.true.)
          call equateDerivedType(u_n ,u_np)
          call equateDerivedType(u_np,vmed)
          if (   (itime > nrestart.and.nrestart > 0)
     .       .or.(time  > trestart.and.trestart > 0d0)) then
            call setEarlyTimeRestartFlag(.true.)
            exit
          endif
        endif

      enddo

c     Setup source

      call equateDerivedType(vmed,u_0)
      call get_src(vmed)

c     Exit

      if (my_rank == 0) write (*,*) ' Done!'

#if defined(pit)
      dt = dtbase  !Reset time step with parallel-in-time
#endif

c     End

      call deallocateDerivedType(vmed)

      call closeRestartFileForRead(ufile)

      end subroutine readRestartFile

      end subroutine readCurrentState

c initializeRecordFile
c######################################################################
      subroutine initializeRecordFile(file)

c----------------------------------------------------------------------
c     Creates graphics pointers, defines dumping intervals
c----------------------------------------------------------------------

      use timeStepping

      use var_io

      implicit none

c Call variables

      character(*) :: file
      
c Local variables

      integer :: ierr,ieq,igrid

      type(var_array),pointer :: u_graph => null()

c Interfaces

c Begin program

#if defined(petsc)
      !Begin profiling IO
      call switch_profile_stage(iPOST)
#endif

      igrid = 1

#if defined(adios)
      ierr = init_ADIOS_IO()
#endif

#if defined(ADIOS2)
      ierr = init_ADIOS2_IO()
#endif

c Open record file

      if (.not.no_graphics_io) then
        if (.not.ts_restart) then
          !Dump initial data
          call writeRecordFile(file,0,0d0,dt,0d0,u_0)
          call writeRecordFile(file,0,0d0,dt,0d0,u_np)
        elseif (getEarlyTimeRestartFlag()) then
          file = new_dumpfile(itime)
          call writeRecordFile(file,itime,time,dt,gammat,u_np)
        else
          call setAppendIOMode
        endif
      endif

#if defined(petsc)
      !End profiling IO
      call switch_profile_stage(iSTUP)
#endif

c End programs

      end subroutine initializeRecordFile

c get_src
c####################################################################
      subroutine get_src(varray)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation. Variables
c     are stored in u_np
c--------------------------------------------------------------------

      use var_setup

      use timeStepping

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      type(var_array),pointer :: varray

c Local variables

      integer :: ieq,igr

      INTERFACE
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use var_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c Begin program      

#if defined(petsc)
        !Begin profiling NL solver
        call switch_profile_stage(iNL)
#endif

      igr = 1

      source_eval = .true.   !This informs nlfunction that we are processing the source

      !This not only evaluates fsrc, but defines BCs on varray
      call evaluate_NLF(igr,varray,fsrc)

      !Spline unperturbed equilibrium if moving mesh
#if defined(MK)
      if (relax_init_grid) call MK_spline_equ(igr,varray)
#endif

      !Dump source
      if (chk_src) call dump_src

      if (.not.source) fsrc = 0d0

      source_eval = .false.

      !Spline source if moving mesh
#if defined(MK)
      if (mk_grid) call MK_spline_src(igr,fsrc)
#endif

#if defined(petsc)
        !End profiling NL solver
        call switch_profile_stage(iSTUP)
#endif


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

        type(var_array),pointer :: u_graph => null()

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '>>>>>>>>>>>>>>>>>'
          write (*,*) 'Dumping source...'
          write (*,*) '>>>>>>>>>>>>>>>>>'
        endif

        nx = gv%gparams%nxv(1)
        ny = gv%gparams%nyv(1)
        nz = gv%gparams%nzv(1)

        xlabel=cont_conf%label(1)
        ylabel=cont_conf%label(2)

        do ieq=1,neqd
          dbg_desc(ieq) = 'dgb'//trim(int2char(ieq))
        enddo

        srcf='src_p'//trim(int2char(my_rank))//'.bin'

        call createDrawInCfile(neqd,trim(srcf),'Source'
     .        ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
     .        ,'drawsrc_p'//trim(int2char(my_rank))//'.in')

        call mapVectorToStructure(gv%gparams,u_graph,fsrc)

        nunit = 110

        open(nunit,file=trim(srcf),form='unformatted'
     .        ,status='replace')

        do ieq=1,neqd
          call C_contour(u_graph%array_var(ieq)%array(1:nx,1:ny,1:nz)
     .                  ,gv%gparams%gxmin,gv%gparams%gxmax
     .                  ,gv%gparams%gymin,gv%gparams%gymax
     .                  ,gv%gparams%gzmin,gv%gparams%gzmax
     .                  ,ieq-1,nunit)
        enddo

        close (nunit)

        call deallocateDerivedType(u_graph)

      end subroutine dump_src

      end subroutine get_src

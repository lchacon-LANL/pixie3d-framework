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

      integer :: mgcf,mgcs

      logical :: petsc_flg,petsc_flg2,petsc_flg3
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
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .    ,"-nxd",nxd,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-nyd",nyd,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-nzd",nzd,petsc_flg,ierr)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-npx",npx,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-npy",npy,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-npz",npz,petsc_flg,ierr)

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-pc_debug",pc_debug,petsc_flg,ierr)

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-bc_debug",bc_debg,petsc_flg,ierr)
      if (petsc_flg) call set_BC_debug(bc_debg)

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-rtol",rtol,petsc_flg,ierr)

#if defined(coarse_MG)
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-mg_coarse",mg_coarse,petsc_flg,ierr)
      if (petsc_flg) call MG_set_coarse_solve(mg_coarse)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-mgcf",mgcf,petsc_flg,ierr)
      if (petsc_flg) call MG_set_proc_coarsen_factor(mgcf)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-mgcs",mgcs,petsc_flg,ierr)
      if (petsc_flg) call MG_set_coarse_size(mgcs)
#endif

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-ilevel",ilevel,petsc_flg,ierr)

      !Should NOT specify dt at run time, as APP internal checks may use input deck one
c$$$      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
c$$$     .   ,"-dt",dt,petsc_flg,ierr)

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-bdf1",petsc_flg,petsc_flg2,ierr)
      if (petsc_flg.and.petsc_flg2) then
        cnfactor = 0d0 ; sm_flag = 0
      endif

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-bdf2",petsc_flg,petsc_flg2,ierr)
      if (petsc_flg.and.petsc_flg2) then
        cnfactor = 0d0 ; sm_flag = 2
      endif

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-cn",petsc_flg,petsc_flg2,ierr)
      if (petsc_flg.and.petsc_flg2) then
        cnfactor = 0.5d0 ; sm_flag = 0
      endif

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-rannacher",petsc_flg,petsc_flg2,ierr)
      if (petsc_flg.and.petsc_flg2) then
        sm_flag = 1
      endif

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-numtime",numtime,petsc_flg,ierr)
      if (petsc_flg) tmax = 0d0

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-ndstep",ndstep,petsc_flg,ierr)
      if (petsc_flg) dstep = 0d0

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-tmax",tmax,petsc_flg,ierr)
      if (petsc_flg) numtime = 0

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-dstep",dstep,petsc_flg,ierr)
      if (petsc_flg) ndstep = 0

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-chk_src",chk_src,petsc_flg,ierr)

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-adapt_dt",ts_adapt_dt,petsc_flg,ierr)

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-restart",ts_restart,petsc_flg2,ierr)
      
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-trestart",trestart,petsc_flg,ierr)
      if (petsc_flg.and.(.not.petsc_flg2)) ts_restart = .true.

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     $   ,"-dt_ramp_up",ts_ramp_up_dt,petsc_flg,ierr)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     $   ,"-dt_ramp_up_frac",ts_ramp_up_dt_fraction,petsc_flg2,ierr)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     $   ,"-dt_ramp_up_recov",ts_ramp_up_dt_recovery,petsc_flg3,ierr)
      if ((petsc_flg2.or.petsc_flg3).and.(.not.petsc_flg))
     .                         ts_ramp_up_dt = .true.

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-fpa",fpa,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-fpa_it",fpa_it,petsc_flg,ierr)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-nwt_it",maxitnwt,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-ksp_it",maxksp,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-gm_it",maxitgm,petsc_flg,ierr)

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-rtol",rtol,petsc_flg,ierr)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-atol",atol,petsc_flg,ierr)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-stol",stol,petsc_flg,ierr)

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-cpu_quota",quota,petsc_flg,ierr)

      call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-test",tst_flg,petsc_flg,ierr)
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

      call readSimulationState(recordfile)

c Initialize record file

      call initializeRecordFile(recordfile)

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        call pstop("initializeCalculation"
     $            ,"Tmax is less or equal than restarting time")
      endif

c End program

      end subroutine initializeCalculation

c$$$c setInitialState
c$$$c####################################################################
c$$$      subroutine setInitialState
c$$$
c$$$c--------------------------------------------------------------------
c$$$c     Set initial state for initial value calculation.
c$$$c--------------------------------------------------------------------
c$$$
c$$$      use grid
c$$$
c$$$      use timeStepping
c$$$
c$$$      use var_io
c$$$
c$$$#if defined(MK)
c$$$      use mk
c$$$#endif
c$$$
c$$$      implicit none
c$$$
c$$$c Local variables
c$$$
c$$$      integer :: ieq,nx,ny,nz,igrid
c$$$
c$$$c Interfaces
c$$$
c$$$      INTERFACE
c$$$        subroutine setEquilibrium(igx,varray)
c$$$        use var_setup
c$$$        integer :: igx
c$$$        type(var_array),pointer :: varray
c$$$        end subroutine setEquilibrium
c$$$      END INTERFACE
c$$$
c$$$      INTERFACE
c$$$        subroutine perturbEquilibrium(varray)
c$$$        use var_setup
c$$$        type(var_array),pointer :: varray
c$$$        end subroutine perturbEquilibrium
c$$$      END INTERFACE
c$$$
c$$$      INTERFACE
c$$$        subroutine get_src(varray)
c$$$        use var_setup
c$$$        type(var_array),pointer :: varray
c$$$        end subroutine get_src
c$$$      END INTERFACE
c$$$
c$$$c Begin program
c$$$
c$$$      igrid = 1
c$$$
c$$$c Create equilibrium u_0
c$$$
c$$$      call setEquilibrium(igrid,u_0)  !External
c$$$
c$$$#if defined(MK)
c$$$      if (mk_grid) then
c$$$        u_0%array_var(neqd)%descr = 'MK phi'
c$$$        u_0%array_var(neqd)%array = 0d0
c$$$
c$$$        call MK_initGrid(igrid,bcs=bcond)  !Change also in 3dplot.F
c$$$      endif
c$$$#endif
c$$$
c$$$      call equateDerivedType(u_n,u_0)
c$$$
c$$$      call equateDerivedType(u_np,u_n)
c$$$
c$$$c Set up source
c$$$      
c$$$      call get_src(u_n)  !This also imposes BCs
c$$$
c$$$      call mapStructureToVector(gv%gparams,x0,u_n)
c$$$
c$$$c Write geometry when applicable
c$$$
c$$$cc      if (coords == 'ext') call writeExtMap(gv%gparams)
c$$$
c$$$c Perturb equilibrium
c$$$
c$$$      if (.not.relax_init_grid) call perturbEquilibrium(u_np)
c$$$
c$$$c Impose BC (to initialize all auxiliary variables)
c$$$
c$$$      call applyBC(igrid,u_np,gv%aux,it=1)
c$$$
c$$$#if defined(MK)
c$$$      if (mk_grid) call MK_applyBC(igrid,u_np%array_var(neqd)%array)
c$$$#endif
c$$$
c$$$c End program
c$$$
c$$$      end subroutine setInitialState

c readSimulationState
c####################################################################
      subroutine readSimulationState(file)

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

      type(var_array),pointer :: vmed => null()

c Interfaces

      INTERFACE
        subroutine setEquilibrium(igx,varray)
        use var_setup
        integer :: igx
        type(var_array),pointer :: varray
        end subroutine setEquilibrium
      END INTERFACE

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

c Read init data

      if (.not.ts_restart) then
        call readEquilibrium(itime,time)
      else
        call readRestartFile(itime,time)
      endif

c Set up source

      call allocateDerivedType(vmed)

      call equateDerivedType(vmed,u_0)
      call get_src(vmed)

      call deallocateDerivedType(vmed)

c Introduce perturbation at restart

      inewtime = itime+1

      if ((perturb_rstrt).or.((.not.relax_init_grid).and.(itime==0)))
     $     call perturbEquilibrium(u_np)

c Impose BCs (to initialize all auxiliary variables)

      old_time_eval = .true. ; init_time_eval = .true.
      call applyBC(igrid,u_np,gv%aux)
      old_time_eval = .false.; init_time_eval = .false.

#if defined(MK)
      if (mk_grid) call MK_applyBC(igrid,u_np%array_var(neqd)%array)
#endif

c End program

      contains

c     readEquilibrium
c     #################################################################
      subroutine readEquilibrium(itime,time)

c     -----------------------------------------------------------------
c     Reads restart file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer,intent(OUT) :: itime
      real(8),intent(OUT) :: time

c     Local variables

c     Begin program

      call setEquilibrium(igrid,u_0)  !External

#if defined(MK)
      if (mk_grid) then
        u_0%array_var(neqd)%descr = 'MK phi'
        u_0%array_var(neqd)%array = 0d0

        call MK_initGrid(igrid,bcs=bcond)  !Change also in 3dplot.F
      endif
#endif

      call equateDerivedType(u_n,u_0)

      call equateDerivedType(u_np,u_n)

      itime = 0
      time  = 0d0
      
      end subroutine readEquilibrium
      
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

      if (my_rank == 0) write (*,'(a)') ' Reading restart file(s)...'

      ierr = readRecordFile(ufile,itime,time,dt,gammat,u_0)

      if (ierr /= 0) then
        call pstop('readSimulationState','Unable to read u_0')
      endif

      call TS_update_preCounters(time,exclude_itime=.true.)

c     Read files

      ierr = readRecordFile(ufile,itime,time,dt,gammat,u_n)

      if (ierr /= 0) then  !Starting from file-based initial condition

        call equateDerivedType(u_n ,u_0)
        call equateDerivedType(u_np,u_n)
cc        call pstop('readSimulationState','Unable to read u_n')
cc        ts_restart = .false.
      else

        call equateDerivedType(u_np,u_n)

        do
          itime_old = itime + 1

          ierr = readRecordFile(ufile,itime,time,dt,gammat,vmed)

          if(ierr==0.and.switch_dumpfile(itime_old,itime,ufile)) then
            ierr = readRecordFile(ufile,itime,time,dt,gammat,vmed)
          endif

          if (ierr /= 0) then   !EOF: Rollback counters and exit
            call TS_rollback_counters(exclude_itime=.true.)
            exit
          else
            call TS_update_preCounters(time,exclude_itime=.true.)
            call equateDerivedType(u_n ,u_np)
            call equateDerivedType(u_np,vmed)
            if (   (itime > nrestart.and.nrestart > 0)
     .         .or.(time  > trestart.and.trestart > 0d0)) then
              call setEarlyTimeRestartFlag(.true.)
              exit
            endif
          endif

        enddo

      endif
      
c     Exit

      if (my_rank == 0) write (*,'(a)') ' Done!'

#if defined(pit)
      dt = dtbase  !Reset time step with parallel-in-time
#endif

c     End

      call deallocateDerivedType(vmed)

      call closeRestartFileForRead(ufile)

#if defined(MK)
      if(mk_grid) call MK_setup_grid(igrid,u_np%array_var(neqd)%array)

      if(mk_grid.and.mk_init_relax.and.itime==1) itime = 0
#endif

      end subroutine readRestartFile

      end subroutine readSimulationState

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
        if ((.not.ts_restart).or.(ts_restart.and.(itime==0))) then
          !Dump initial data
          call writeRecordFile(file,0,0d0,dt,0d0,u_0)
          call writeRecordFile(file,0,0d0,dt,0d0,u_np)
c$$$        elseif (getEarlyTimeRestartFlag()) then
c$$$          file = new_dumpfile(itime)
c$$$          call writeRecordFile(file,itime,time,dt,gammat,u_np)
c$$$        else
c$$$          call setAppendIOMode
        else
          file = new_dumpfile(itime)
          call writeRecordFile(file,itime,time,dt,gammat,u_np)
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
      use mg_io

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
      if (chk_src) then
cc        call dump_src
        call vec_plot_gl(gv%gparams,neqd,fsrc,igr,'src','Source')
      endif
      
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

c$$$      contains
c$$$
c$$$c     dump_src
c$$$c     ###############################################################
c$$$      subroutine dump_src
c$$$
c$$$        use grid_debug,ONLY:C_contour,cont_conf,dbg,createDrawInCfile
c$$$
c$$$        integer :: ieq,nx,ny,nz,nunit
c$$$
c$$$        character(10) :: srcf
c$$$        character(4)  :: dbg_desc(neqd)
c$$$        character(3)  :: xlabel,ylabel
c$$$
c$$$        type(var_array),pointer :: u_graph => null()
c$$$
c$$$        if (my_rank == 0) then
c$$$          write (*,*)
c$$$          write (*,*) '>>>>>>>>>>>>>>>>>'
c$$$          write (*,*) 'Dumping source...'
c$$$          write (*,*) '>>>>>>>>>>>>>>>>>'
c$$$        endif
c$$$        
c$$$        nx = gv%gparams%nxv(1)
c$$$        ny = gv%gparams%nyv(1)
c$$$        nz = gv%gparams%nzv(1)
c$$$
c$$$        xlabel=cont_conf%label(1)
c$$$        ylabel=cont_conf%label(2)
c$$$
c$$$        do ieq=1,neqd
c$$$          dbg_desc(ieq) = 'dgb'//trim(int2char(ieq))
c$$$        enddo
c$$$
c$$$        srcf='src_p'//trim(int2char(my_rank))//'.bin'
c$$$
c$$$        call createDrawInCfile(neqd,trim(srcf),'Source'
c$$$     .        ,'t',trim(xlabel),trim(ylabel),dbg_desc,'-c -X0 -L57'
c$$$     .        ,'drawsrc_p'//trim(int2char(my_rank))//'.in')
c$$$
c$$$        call mapVectorToStructure(gv%gparams,u_graph,fsrc)
c$$$
c$$$        nunit = 110
c$$$
c$$$        open(nunit,file=trim(srcf),form='unformatted'
c$$$     .        ,status='replace')
c$$$
c$$$        do ieq=1,neqd
c$$$          call C_contour(u_graph%array_var(ieq)%array(1:nx,1:ny,1:nz)
c$$$     .                  ,gv%gparams%gxmin,gv%gparams%gxmax
c$$$     .                  ,gv%gparams%gymin,gv%gparams%gymax
c$$$     .                  ,gv%gparams%gzmin,gv%gparams%gzmax
c$$$     .                  ,ieq-1,nunit)
c$$$        enddo
c$$$
c$$$        close (nunit)
c$$$
c$$$        call deallocateDerivedType(u_graph)
c$$$
c$$$      end subroutine dump_src

      end subroutine get_src

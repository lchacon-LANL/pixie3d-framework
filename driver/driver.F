      program driver_3d

c ******************************************************************
c  This program drives the time stepping of an arbitrary 3D system.
c ******************************************************************

      use variables

      use timeStepping

      use nk

      use var_io

      use grid

#if defined(MK)
      use mk, ONLY: mk_grid,MK_output
#endif

c Common variables

      implicit none

c Local variables

      integer :: ierr=0
      
c Interfaces

      INTERFACE
        subroutine output(varr)
        use var_setup
        type(var_array),pointer :: varr
        end subroutine output
      END INTERFACE

      INTERFACE
        subroutine preProcessSol(vnm,vn)
        use var_setup
        type(var_array),pointer :: vn,vnm
        end subroutine preProcessSol
      END INTERFACE

      INTERFACE
        subroutine postProcessSol(vnm,vn,vnp)
        use var_setup
        type(var_array),pointer :: vnm,vn,vnp
        end subroutine postProcessSol
      END INTERFACE

c Begin program

      call DR_init

c Initialize calculation

      call initializeCalculation

c Initialize timers

      call TS_init_timers

c Initial output

      call output(u_np)  !External

#if defined(MK)
      if (mk_grid) call MK_output
#endif

c Time loop

      do

        if (TS_chk_end()) exit

c     Find new time step and update time counters

        call TS_correct_dt(ierr)

c     Store old time solution

        call equateDerivedType(u_nm,u_n)
        call equateDerivedType(u_n ,u_np)

c     Preprocess old solution (for e.g. splitting)

        if (preprocess) call preProcessSol(u_nm,u_n) !External

c     Time update

        ierr = TS_update(u_np,ilevel)

c     Check for error in time stepping

        select case(ierr)
        case(NOCONV)    !Restart time step

          call equateDerivedType(u_np,u_n )
          call equateDerivedType(u_n ,u_nm)

          call TS_rollback_counters()

          if (ts_adapt_dt) then
            cycle
          else
            call pstop('Main'
     .                ,'ERROR: Solver failed to converge;'
     .                 //' rerun with "ts_adapt_dt = t"')
          endif

        case(EXACT)     !Initial guess is exact

          if (my_rank == 0) then
            write(*,*)
            write(*,*) '   Found steady state solution.'
            write(*,*) '   Aborting...'
          endif

          if (tst_flg) EXIT_STATUS = SUCCESS

          exit

        end select

c     Update other TS counters

        call TS_update_postCounters

c     Postprocess solution

        if (postprocess) call postProcessSol(u_nm,u_n,u_np)  !External

c     Output per time step

#if defined(MK)
        if (mk_grid) call MK_output
#endif

        call output(u_np)

c     Time level data dump

        call TS_dump(u_np,dump_now=TS_chk_end())  !Force restart dump at end

c     Exit if eigenvalue computation

        if (eigen) exit

      enddo       !End of time loop

c End timers

      call TS_end_timers

c Final I/O

      call TS_final_io

c End program

      call DR_destroy

      call exit(EXIT_STATUS)

      end program

c DR_init
c ######################################################################
      subroutine DR_init

      use var_io

      use grid

      use ts_setup

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr,provided

c Begin program

c Initialize PETSc and MPI (parallel version)

#if defined(petsc)
CGMM_OMP      call MPI_init_thread(MPI_THREAD_MULTIPLE,provided,ierr)
CGMM_OMP      if (provided < MPI_THREAD_MULTIPLE) then
CGMM_OMP        call pstop("DR_init","Cannot provide MPI threading requested")
CGMM_OMP      endif

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      
      !Setup profiling
      call init_petsc_profile

      call initMPI(MPI_COMM_WORLD,np,my_rank)

      !Profiling stage 1
      call switch_profile_stage(iSTUP)
#endif

c End program

      end subroutine DR_init

c DR_destroy
c ######################################################################
      subroutine DR_destroy

      use var_io

      use grid

      use ts_setup

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr

c Begin program

c Free memory

      call deallocateStructures

      call TS_destroy

      call destroyNLF   !External

#if defined(MK)
      if (mk_grid) call MK_dealloc
#endif

      call deallocateGlobalVar(gv)

      call finalize_IO

#if defined(petsc)
      call end_petsc_profile

      call PetscFinalize(ierr)
#endif

c End program

      end subroutine DR_destroy

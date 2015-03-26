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

      integer :: ierr=0,offset

c Interfaces

      INTERFACE
        subroutine output(varr)
        use var_setup
        type(var_array),pointer :: varr
        end subroutine output
      END INTERFACE

      INTERFACE
        subroutine preProcessSol(vn,vnm)
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

      INTERFACE
        subroutine setInitialCondition(vn,vnp)
        use var_setup
        type(var_array),pointer :: vn,vnp
        end subroutine setInitialCondition
      END INTERFACE

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

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

      call Init

c Setup profiling

#if defined(petsc)
      call init_petsc_profile

      !Profiling stage 1
      call switch_profile_stage(iSTUP)
#endif

c Initialize calculation

      call initializeCalculation

      call nk_setup_consistency

c Initialize counters

      call TS_init_otherCounters

c Initial output

      call output(u_np)  !External

#if defined(MK)
      if (mk_grid) call MK_output
#endif

c Time loop

#if defined(petsc)
      call PetscGetTime   (prof_time(1),prof_err)
      call PetscGetCPUTime(prof_time(2),prof_err)

      !Begin profiling third stage
      call switch_profile_stage(iPOST)
#endif

      offset = 0
      if (relax_init_grid) offset = 1

      do

        if (TS_chk_end()) exit

c     Find new time step and update time counters

        call TS_correct_dt(offset,ierr)

c     Store old time solution

        call equateDerivedType(u_nm,u_n)
        call equateDerivedType(u_n ,u_np)

c     Preprocess old solution (for e.g. splitting)

        if (preprocess) call preProcessSol(u_nm,u_n) !External

c     Evaluate old-time NLF

#if defined(petsc)
        !Begin profiling NL solver
        call switch_profile_stage(iNL)
#endif

        old_time_eval = .true.  !This informs nlfunction that we are processing
                                !  n time level info

        call evaluate_NLF(1,u_n,fold)

        old_time_eval = .false.

c     Store BDF2 old-time solution (in u_np)

        if (bdf2) call AXPYDerivedType(-cn/cnp,u_n,-cnm/cnp,u_nm,u_np)

c     Time update

        if (eigen) then
          !Eigenvalue solver
          ierr = eigenSolution (u_np,ilevel)
        elseif (prj_int_ts) then
          ierr = prj_int_update(u_np,ilevel)
        else
          !Time step
          ierr = updateSolution(u_np,ilevel)
        endif

c     Check for error in time stepping

        if (ierr == NOCONV) then     !Restart time step

          call equateDerivedType(u_np,u_n )
          call equateDerivedType(u_n ,u_nm)

          call TS_rollback_timecounters()

          itime = itime - 1

cc          time  = tn
cc          tn    = tnm
cc          tnm   = 0d0

          if (timecorr) then
            cycle
          else
            exit
          endif

        elseif (ierr == EXACT) then !Initial guess is exact

          if (my_rank == 0) then
            write(*,*)
            write(*,*) '   Found steady state solution.'
            write(*,*) '   Aborting...'
          endif

          exit

        endif

c     Postprocess solution

#if defined(petsc)
        !Begin profiling third stage
        call switch_profile_stage(iPOST)
#endif

        if (postprocess) call postProcessSol(u_nm,u_n,u_np)  !External

c     Update other counters (only if timeStep is successful)

        if (.not.relax_init_grid) call TS_update_otherCounters

c     Output per time step

#if defined(MK)
        if (mk_grid) call MK_output
#endif

        call output(u_np)

c     Time level data dump

#if defined(RFX)
        if (nrst.eq.ndstep.or.tmrst.ge.dstep
#else
        if (nrst.eq.ndstep.or.tmrst.ge.0.999999*dstep
#endif
     .      .or.relax_init_grid) then
          if (.not.relax_init_grid) nrst  = 0
          if(itime-offset > 0) tmrst=tmrst-dstep
          call writeRecordFile(recordfile,itime,time,dt,u_np)
        endif

        if (eigen) exit

        if (relax_init_grid) then
          relax_init_grid = .false.
cc          if (source) call get_src(u_np,.true.)  !Update source with interpolated equilibrium
          call equateDerivedType (u_n,u_np)     !Save interpolated (unperturbed) equilibrium
          call perturbEquilibrium(u_np)         !Perturb interpolated equilibrium
        endif

      enddo       !End of time loop

c Average explicit time step

      if (cnfactor > 5d-1.and.my_rank == 0.and.(itime > inewtime-1))then
        dtexp = dtexp/(itime - inewtime+1)
        write (*,'(a,1pe10.4)') 'Average explicit time step: ',dtexp
      endif

c Final statistics

#if defined(petsc)
      call PetscGetTime   (prof_time(3),prof_err)
      call PetscGetCPUTime(prof_time(4),prof_err)

      if (my_rank == 0 .and. .not.(tst_flg)) then
        write (*,*) 'CPU statistics'
        write (*,*)
        write (*,*) 'Elapsed time =', prof_time(3)-prof_time(1)
        write (*,*) 'CPU time     =',(prof_time(4)-prof_time(2))*1d4
        write (*,*) 'COMM % time  ='
     .      ,(prof_time(3)-prof_time(1)-(prof_time(4)-prof_time(2))*1d4)
     .      /(prof_time(3)-prof_time(1))*100
      endif
#endif

      if (my_rank == 0) then
cc        prec_tot = gmres_tot + iguess*newt_tot

        write(*,300) 

        if (itime-inewtime+1 > 0 .and. newt_tot > 0) then
          write(*,310) itime,float(newt_tot )/(itime-inewtime+1)
     .                      ,float(gmres_tot)/(itime-inewtime+1)
     .                      ,float(gmres_tot)/newt_tot
        else
          write(*,310) 0,0d0,0d0,0d0
        endif
      endif

c Formats

 300  format (/,'Final statistics',/,/,
     .          ' itime  Newt/ndt  GMRES/ndt  GMRES/Newt')
 310  format (i5,3x,f7.1,4x,f7.1,4x,f7.1)

c End program

      call Destroy

      end program

c Init
c ######################################################################
      subroutine Init

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

c Initialize PETSc and MPI (parallel version)

#if defined(petsc)
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      call initMPI(MPI_COMM_WORLD,np,my_rank)
#endif

c Initialize ADIOS for componentized I/O

#if defined(adios) && defined(petsc)
      if (.not.is_file('adios_config.xml')) then
        call pstop('driver_3d','ADIOS config file not present')
      else
        call MPI_Comm_dup (MPI_COMM_WORLD, adios_world_comm, ierr)
        call ADIOS_init('adios_config.xml', adios_err)
      endif
#endif

c End program

      end subroutine Init

c Destroy
c ######################################################################
      subroutine Destroy

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

      call destroyTS

      call destroyNLF   !External

#if defined(MK)
      if (mk_grid) call MK_dealloc
#endif

      call deallocateGlobalVar(gv)

#if defined(adios)
      call ADIOS_finalize(my_rank, adios_err)
#endif

#if defined(petsc)
      call end_petsc_profile

      call PetscFinalize(ierr)
#endif

c End program

      end subroutine Destroy

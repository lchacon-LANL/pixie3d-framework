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

      integer :: prec_tot,ierr,perr,prof_stage(3),offset,ieq
      real(8) :: time1,time2,time3,time4

c Interfaces

cc      INTERFACE
cc        subroutine timeStep(vn,vnp,ierr)
cc        use var_setup
cc        integer          :: ierr
cc        type(var_array),pointer :: vnp,vn
cc        end subroutine timeStep
cc      END INTERFACE

      INTERFACE
        subroutine output(varr)
        use var_setup
        type(var_array),pointer :: varr
        end subroutine output
      END INTERFACE

      INTERFACE
        subroutine preProcessSol(v_n,v_nm)
        use var_setup
        type(var_array),pointer :: v_n,v_nm
        end subroutine preProcessSol
      END INTERFACE

      INTERFACE
        subroutine postProcessSol(v_arr)
        use var_setup
        type(var_array),pointer :: v_arr
        end subroutine postProcessSol
      END INTERFACE

      INTERFACE
        subroutine setInitialCondition(varrayn,varraynp)
        use var_setup
        type(var_array),pointer :: varrayn,varraynp
        end subroutine setInitialCondition
      END INTERFACE

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

c Begin program

c Initialize PETSc and MPI (parallel version)

#if defined(petsc)
      call PetscInitialize(PETSC_NULL_CHARACTER,perr)

      call initMPI(MPI_COMM_WORLD,np,my_rank)

      !Set profiling stages
      call PetscLogStageRegister('Initial setup'  ,prof_stage(1),perr)
      call PetscLogStageRegister('Nonlinear solve',prof_stage(2),perr)
      call PetscLogStageRegister('Sol. processing',prof_stage(3),perr)  
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

c Initialize calculation

#if defined(petsc)
      !Profiling stage 1
      call PetscLogStagePush(prof_stage(1),perr)
#endif

      call initializeCalculation

      call nk_setup_consistency

c Initialize counters

      gmres_tot = 0
      newt_tot  = 0
      ierr      = 0

      itime     = inewtime - 1

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0

c Initial output

      call output(u_np)

#if defined(MK)
      if (mk_grid) call MK_output
#endif

c Time loop

#if defined(petsc)
      call PetscGetTime(time1,perr)
      call PetscGetCPUTime(time2,perr)

      !Begin profiling third stage
      call PetscLogStagePop(perr)
      call PetscLogStagePush(prof_stage(3),perr)
#endif

      offset = 0
      if (relax_init_grid) offset = 1

      do

        itime = itime + 1

        if (tmax.gt.0d0.and.time.ge.(tmax-1d-9*dt))  exit
cc        if (tmax.gt.0d0.and.time.gt.tmax)  exit
        if (numtime.ge.0.and.itime.ge.(numtime+inewtime)) exit

c     Find new time step

        if (    cnfactor == 1d0                            !Explicit
     .     .or.(timecorr.and.ierr==0.and.itime-offset==1)  !First explicit step if timecorr
     .     ) then
          call findExplicitDt(dt)
          call correctTimeStep(itime-offset,ierr,dt_exp=dt)
cc        elseif (timecorr) then
        else
          call correctTimeStep(itime-offset,ierr)
        endif

c     Assign old time solution

        call equateDerivedType(u_nm,u_n)
        call equateDerivedType(u_n ,u_np)

c     Preprocess old solution (for e.g. splitting)

        if (preprocess) call preProcessSol(u_nm,u_n) !External

c     Store BDF2 old-time solution (in u_np)

        if (bdf2) then
          call AXPYDerivedType(-cn/cnp,u_n,-cnm/cnp,u_nm,u_np)
        endif

c     Time update

#if defined(petsc)
        !Begin profiling second stage
        call PetscLogStagePop(perr)
        call PetscLogStagePush(prof_stage(2),perr)
#endif

        if (.not.relax_init_grid) then
          tnm = tn
          tn  = time
          time = tn + dt
        endif

        if (eigen) then
          !Eigenvalue solver
          ierr = eigenSolution(u_np,ilevel)
        else
          !Time step
          ierr = updateSolution(dt,u_np,ilevel)
cc          call timeStep(u_n,u_np,ierr)
        endif

c     Check for error in time stepping

        if (ierr == NOCONV) then     !Restart time step

          call equateDerivedType(u_np,u_n )
          call equateDerivedType(u_n ,u_nm)

          itime = itime - 1

          time  = tn
          tn    = tnm

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
        call PetscLogStagePop(perr)
        call PetscLogStagePush(prof_stage(3),perr)
#endif

        if (postprocess) call postProcessSol(u_np)  !External

c     Update counters (only if timeStep is successful)

        if (.not.relax_init_grid) then
          tmrst  = tmrst + dt
          nrst   = nrst  + 1
        endif

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

      if (cnfactor.eq.1d0.and.my_rank == 0.and.(itime > inewtime)) then
        dtexp = dtexp/(itime - inewtime)
        write (*,'(a,1pe10.4)') 'Average explicit time step: ',dtexp
      endif

c Final statistics

#if defined(petsc)
      call PetscGetTime(time3,perr)
      call PetscGetCPUTime(time4,perr)

      if (my_rank == 0 .and. .not.(tst_flg)) then
        write (*,*) 'CPU statistics'
        write (*,*)
        write (*,*) 'Elapsed time =', time3-time1
        write (*,*) 'CPU time     =',(time4-time2)*1d4
        write (*,*) 'COMM % time  =',(time3-time1-(time4-time2)*1d4)
     .                              /(time3-time1)*100
      endif
#endif

      if (my_rank == 0) then
        prec_tot = gmres_tot + iguess*newt_tot

        write(*,300) 

        if (itime-inewtime > 0 .and. newt_tot > 0) then
          write(*,310) (itime-1),float(newt_tot )/(itime-inewtime)
     .                          ,float(gmres_tot)/(itime-inewtime)
     .                          ,float(gmres_tot)/newt_tot
        else
          write(*,310) 0,0d0,0d0,0d0
        endif
      endif

c Formats

 300  format (/,'Final statistics',/,/,
     .          ' itime  Newt/ndt  GMRES/ndt  GMRES/Newt')
 310  format (i5,3x,f7.1,4x,f7.1,4x,f7.1)

c End program

      call fortranDestroy

#if defined(petsc)
      call PetscLogStagePop(perr)

      call PetscFinalize(perr)
#endif

#if defined(adios)
      call ADIOS_finalize(my_rank, adios_err)
#endif

      end program

c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

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

      call destroyNonlinearFunction   !External

#if defined(MK)
      if (mk_grid) call MK_dealloc
#endif

      call deallocateGlobalVar(gv)

c End program

      end subroutine fortranDestroy

ccc timeStep
ccc######################################################################
cc      subroutine timeStep(vn,vnp,ierr)
cc
ccc----------------------------------------------------------------------
ccc     Performs time advance of solution from vn to vnp
ccc----------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use var_setup
cc
cc      use timeStepping
cc
cc      use nk
cc
cc      use var_io
cc
cc      use eigen_mod
cc
cc#if defined(MK)
cc      use mk
cc#endif
cc
cc      implicit none
cc
ccc Call variables
cc
cc      integer :: ierr
cc
cc      type(var_array),pointer :: vnp,vn
cc
ccc Local variables
cc
cc      real(8) :: x(ntotd)
cc
ccc Interfaces
cc
cc      INTERFACE
cc        subroutine evaluate_NLF(igr,varray,fi)
cc        use parameters
cc        use var_setup
cc        integer :: igr
cc        real(8) :: fi(ntotd)
cc        type(var_array),pointer :: varray
cc        end subroutine evaluate_NLF
cc      END INTERFACE
cc
cccc      INTERFACE
cccc        subroutine explicit(vnp,iter,iout)
cccc        use var_setup
cccc        integer :: iter,iout
cccc        type(var_array),pointer :: vnp
cccc        end subroutine explicit
cccc      END INTERFACE
cc
ccc Begin program
cc
cc      ierr = 0
cc
ccc Evaluate nonlinear function at old time for theta scheme
cc
cc      old_time_eval = .true.     !This informs nlfunction that we are processing
cc                                 !  n time level info
cc
cc      call evaluate_NLF(1,vn,fold)
cc
cc      old_time_eval = .false.
cc 
ccc Map old solution (use vnp in case we are doing BDF2)
cc
cc      call mapStructureToVector(x_old,vnp)
cc
ccc Implicit update (Newton)
cc
cc      if (cnfactor.lt.1d0) then
cc
ccc     Initial guess
cc
cc        if (predictor.and.(.not.bdf2)) then
cc          call findGuess(vn,vnp)
cc          call mapStructureToVector(x,vnp)
cc        else
cc          x = x_old
cc        endif
cc
ccc     Newton-Krylov iteration
cc
cc        nk_conf%etak_meth      = method
cc        nk_conf%ksmax          = maxksp
cc        nk_conf%gmmax          = maxitgm
cc        nk_conf%nwt_max_it_acc = 0.7*maxitnwt
cc        nk_conf%nwt_max_it_rej = maxitnwt
cc        nk_conf%global_meth    = global
cc
cc        nk_conf%eta0   = tolgm
cc        nk_conf%damp   = damp
cc        nk_conf%pdt0   = dt0
cc        if (atol  > 0d0) nk_conf%atol   = atol
cc        if (rtol  > 0d0) nk_conf%rtol   = rtol
cc        if (stol  > 0d0) nk_conf%stol   = stol
cc        if (mf_eps> 0d0) nk_conf%mf_eps = mf_eps
cc
cc        nk_conf%krylov_method='fg'
cccc        nk_conf%krylov_method='gm'
cc        nk_conf%fpa=fpa
cc        nk_conf%fpa_it=fpa_it
cc
cc        nk_conf%vol_wgt = vol_wgt
cc
cc        call jfnk(neqd,ntotd,x,iguess,ilevel,ierr)
cc
cc        itgmres = nk_conf%gm_it_out
cc        itnewt  = nk_conf%nwt_it_out
cc
ccc     Map Newton solution to vnp (if no error)
cc
cc        if (ierr.eq.0.or.ierr.eq.2) then
cc          call mapVectorToStructure(vnp,x)
cc
cc          call applyBC(1,vnp,it=itime)
cc
cc#if defined(MK)
cc          if (mk_grid) call MK_applyBC(1,vnp%array_var(neqd)%array)
cc#endif
cc
cc          if (predictor) call storeTSinfo(vn,vnp)
cc        endif
cc
cc      else
cc
ccc Explicit update
cc
cc        itgmres = maxitgm
cc
cccc        call explicit(vnp,itgmres,ilevel)
cc
cc        call applyBC(1,vnp,it=itime)
cc
cc#if defined(MK)
cc        if (mk_grid) call MK_applyBC(1,vnp%array_var(neqd)%array)
cc#endif
cc
cc        itnewt = 1
cc
cc      endif
cc
ccc Update counters
cc
cc      gmres_tot = gmres_tot + itgmres
cc      newt_tot  = newt_tot  + itnewt
cc
ccc End program
cc
cc      end subroutine timeStep

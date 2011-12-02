      program driver_3d

c ******************************************************************
c  This program drives the time stepping of an arbitrary 3D system.
c ******************************************************************

      use variables

      use timeStepping

      use newtongm

      use counters

      use iosetup

      use grid

#if defined(MK)
      use mk, ONLY: mk_grid,MK_output
#endif

c Common variables

      implicit none

c Local variables

      integer    :: prec_tot,ierr,perr,prof_stage(3),offset,ieq
      real(8)    :: time1,time2,time3,time4

c Interfaces

      INTERFACE
        subroutine timeStep(vn,vnp,ierr)
        use variable_setup
        integer          :: ierr
        type(var_array),pointer :: vnp,vn
        end subroutine timeStep
      END INTERFACE

      INTERFACE
        subroutine eigen_sol(vn,vnp,ierr)
        use variable_setup
        integer          :: ierr
        type(var_array),pointer :: vnp,vn
        end subroutine eigen_sol
      END INTERFACE

      INTERFACE
        subroutine output(varr)
        use variable_setup
        type(var_array),pointer :: varr
        end subroutine output
      END INTERFACE

      INTERFACE
        subroutine preProcessSol(v_n,v_nm)
        use variable_setup
        type(var_array),pointer :: v_n,v_nm
        end subroutine preProcessSol
      END INTERFACE

      INTERFACE
        subroutine postProcessSol(v_arr)
        use variable_setup
        type(var_array),pointer :: v_arr
        end subroutine postProcessSol
      END INTERFACE

      INTERFACE
        subroutine setInitialCondition(varrayn,varraynp)
        use variable_setup
        type(var_array),pointer :: varrayn,varraynp
        end subroutine setInitialCondition
      END INTERFACE

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use variable_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

c Begin program

c Allocate global variable

      call allocateGlobalVar(gv)
      
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
        call adios_init('adios_config.xml', adios_err)
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

        if (.not.relax_init_grid) time = time + dt

        if (eigen) then
          !Eigenvalue solver
          call eigen_sol(u_n,u_np,ierr)
        else
          !Time step
          call timeStep(u_n,u_np,ierr)
        endif

c     Check for error in time stepping

#if defined(petsc)
        !Begin profiling third stage
        call PetscLogStagePop(perr)
        call PetscLogStagePush(prof_stage(3),perr)
#endif

        if (ierr == 1) then     !Restart time step

          call equateDerivedType(u_np,u_n )
          call equateDerivedType(u_n ,u_nm)

          itime = itime - 1
          time  = time - dt

          if (timecorr) then
            cycle
          else
            exit
          endif

        elseif (ierr == -1) then !Initial guess is exact

          if (my_rank == 0) then
            write(*,*)
            write(*,*) '   Found steady state solution.'
            write(*,*) '   Aborting...'
          endif

          exit

        endif

c     Postprocess solution

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
          call writeRecordFile(urecord,itime,time,dt,u_np)
        endif

        if (eigen) exit

        if (relax_init_grid) then
          relax_init_grid = .false.
cc          if (source) call get_src(u_np,.true.)  !Update source with interpolated equilibrium
          call equateDerivedType (u_n,u_np)     !Save interpolated (unperturbed) equilibrium
          call perturbEquilibrium(u_np)         !Perturb interpolated equilibrium
cc          call applyBC(u_np,1,1,1)
        endif

      enddo       !End of time loop

c Average explicit time step

      if (cnfactor.eq.1d0.and.my_rank == 0.and.(itime > inewtime)) then
        dtexp = dtexp/(itime - inewtime)
        write (*,*) 'Average explicit time step: ',dtexp
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

        if (itime-inewtime > 0) then
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

      end

c timeStep
c######################################################################
      subroutine timeStep(vn,vnp,ierr)

c----------------------------------------------------------------------
c     Performs time advance of solution from vn to vnp
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use timeStepping

      use newtongm

      use counters

      use iosetup

      use eigen_mod

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ierr

      type(var_array),pointer :: vnp,vn

c Local variables

      real(8) :: x(ntotd)

c Interfaces

      INTERFACE
        subroutine evaluateNonlinearFunction(igr,varray,fi)
        use parameters
        use variable_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluateNonlinearFunction
      END INTERFACE

      INTERFACE
        subroutine explicit(vnp,iter,iout)
        use variable_setup
        integer :: iter,iout
        type(var_array),pointer :: vnp
        end subroutine explicit
      END INTERFACE

      INTERFACE
        subroutine applyBC(varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine applyBC
      END INTERFACE

c Begin program

      ierr = 0

c Evaluate nonlinear function at old time for theta scheme

      old_time_eval = .true.     !This informs nlfunction that we are processing
                                 !  n time level info

      call evaluateNonlinearFunction(1,vn,fold)

      old_time_eval = .false.
 
c Map old solution (use vnp in case we are doing BDF2)

      call mapStructureToVector(x_old,vnp)

c Implicit update (Newton)

      if (cnfactor.lt.1d0) then

c     Initial guess

        if (predictor.and.(.not.bdf2)) then
          call findGuess(vn,vnp)
          call mapStructureToVector(x,vnp)
        else
          x = x_old
        endif

c     Newton-Krylov iteration

        nk_conf%etak_meth      = method
        nk_conf%ksmax          = maxksp
        nk_conf%gmmax          = maxitgm
        nk_conf%nwt_max_it_acc = 0.7*maxitnwt
        nk_conf%nwt_max_it_rej = maxitnwt
        nk_conf%global_meth    = global

        nk_conf%eta0   = tolgm
        nk_conf%damp   = damp
        nk_conf%pdt0   = dt0
        if (atol  > 0d0) nk_conf%atol   = atol
        if (rtol  > 0d0) nk_conf%rtol   = rtol
        if (stol  > 0d0) nk_conf%stol   = stol
        if (mf_eps> 0d0) nk_conf%mf_eps = mf_eps

        nk_conf%krylov_method='fg'
cc        nk_conf%krylov_method='gm'
        nk_conf%fpa=fpa

        nk_conf%vol_wgt = vol_wgt

        call jfnk(neqd,ntotd,x,iguess,ilevel,ierr)

        itgmres = nk_conf%gm_it_out
        itnewt  = nk_conf%nwt_it_out

c     Map Newton solution to vnp (if no error)

        if (ierr.eq.0.or.ierr.eq.2) then
          call mapVectorToStructure(vnp,x)

          call applyBC(vnp,1,1,1)

#if defined(MK)
          if (mk_grid) call MK_applyBC(1,vnp%array_var(neqd)%array)
#endif

          if (predictor) call storeTSinfo(vn,vnp)
        endif

      else

c Explicit update

        itgmres = maxitgm

        call explicit(vnp,itgmres,ilevel)

        call applyBC(vnp,1,1,1)

#if defined(MK)
        if (mk_grid) call MK_applyBC(1,vnp%array_var(neqd)%array)
#endif

        itnewt = 1

      endif

c Update counters

      gmres_tot = gmres_tot + itgmres
      newt_tot  = newt_tot  + itnewt

c End program

      end subroutine timeStep

c eigen_sol
c######################################################################
      subroutine eigen_sol(vn,vnp,ierr)

c----------------------------------------------------------------------
c     Performs time advance of solution from vn to vnp
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use timeStepping

      use newtongm

      use counters

      use iosetup

      use eigen_mod

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ierr

      type(var_array),pointer :: vnp,vn

c Local variables

      real(8) :: x(ntotd)

c Interfaces

      INTERFACE
        subroutine evaluateNonlinearFunction(igr,varray,fi)
        use parameters
        use variable_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluateNonlinearFunction
      END INTERFACE

      INTERFACE
        subroutine applyPC(n,x,y,out)
        integer :: n,out
        real(8),dimension(n) :: x,y
        end subroutine
      END INTERFACE

      INTERFACE
        subroutine applyBC(varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine applyBC
      END INTERFACE

c Begin program

      ierr = 0

c Evaluate nonlinear function at old time for theta scheme

      old_time_eval = .true.     !This informs nlfunction that we are processing
                                 !  n time level info

      call evaluateNonlinearFunction(1,vn,fold)

      old_time_eval = .false.

c Eigenvalue solver

      !Setup matrix-free vectors (xk, rk)
      allocate(nk_conf%xk(ntotd),nk_conf%rk(ntotd))

      call mapStructureToVector(nk_conf%xk,u_0)
      call evaluateNewtonResidual(ntotd,nk_conf%xk,nk_conf%rk)

      !Setup preconditioner
      call setupPC(ntotd,nk_conf%xk,1)

      !Setup initial guess for eigenfunction
      call mapStructureToVector(x ,vnp)
      x = x - nk_conf%xk

        !Eigensolve
cc#if defined(arpack)
cc        call arnoldi_eig(neqd,ntotd,x,ilevel)
cc#endif
      eigen_conf%my_rank = my_rank

      eigen_conf%ksmax   = maxksp
      eigen_conf%armax   = maxitgm
      if (rtol > 0d0) eigen_conf%rtol = rtol
      if (atol > 0d0) eigen_conf%atol = atol
      eigen_conf%step    = eig_dt

      iguess = 1
      call lehoucq_eig(ntotd,x,iguess,ilevel,ierr
     .                ,matrixFreeMatVec
     .                ,applyPC
     .                ,ldot=dot2)

      select case(ierr)
      case(1)
        call pstop('lehoucq_eig','Matrix not positive definite')
      end select

      !Output
      if (my_rank == 0) then
         write (*,*)
         write (*,*) '**************************************'
         write (*,'(a,1p,1e10.3,a)') ' ****** Eigenvalue is='
     $           ,1d0/dt-eigen_conf%eigv,' ******' 
         write (*,*) '**************************************'
      endif

      itgmres = eigen_conf%iter_out

      !Map eigenvector for plotting
      x = x + nk_conf%xk
      call mapVectorToStructure(vnp,x)

      !Free memory
      deallocate(nk_conf%xk,nk_conf%rk)
 
c End program

      end subroutine eigen_sol

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

cc      use icond

c Common variables

      implicit none

c Local variables

      integer    :: prec_tot,ierr,perr,prof_stage(3)
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
        subroutine output(varr)
        use variable_setup
        type(var_array),pointer :: varr
        end subroutine output
      END INTERFACE

      INTERFACE
        subroutine postProcessSol(v_arr)
        use variable_setup
        type(var_array),pointer :: v_arr
        end subroutine postProcessSol
      END INTERFACE

c Begin program

c Allocate global variable

      call allocateGlobalVar(gv)
      
c Initialize PETSc and MPI (parallel version)

#if defined(petsc)
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      !Set profiling stages
      call PetscLogStageRegister('Initial setup'  ,prof_stage(1),perr)
      call PetscLogStageRegister('Nonlinear solve',prof_stage(2),perr)
      call PetscLogStageRegister('Sol. processing',prof_stage(3),perr)  
#endif

c Initialize calculation

#if defined(petsc)
      !Profiling stage 1
      call PetscLogStagePush(prof_stage(1),perr)
#endif

      call initializeCalculation

c Initialize counters

      gmres_tot = 0
      newt_tot  = 0
      wh_tot    = 0
      ierr      = 0

      itime     = inewtime - 1

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0

c Initial output

      call output(u_np)

c Time loop

#if defined(petsc)
      call PetscGetTime(time1,perr)
      call PetscGetCPUTime(time2,perr)

      !Begin profiling third stage
      call PetscLogStagePop(perr)
      call PetscLogStagePush(prof_stage(3),perr)
#endif

      do

        itime = itime + 1

        if (tmax.gt.0d0.and.time.ge.(tmax-1d-5*dt))  exit
cc        if (tmax.gt.0d0.and.time.gt.tmax)  exit
        if (numtime.ge.0.and.itime.ge.(numtime+inewtime)) exit

c     Find new time step

        if (    cnfactor == 1d0                          !Explicit
     .     .or.(timecorr.and.ierr.eq.0.and.itime.eq.1)   !First explicit step if timecorr
     .     ) then
          call findExplicitDt
          call correctTimeStep(itime,ierr,dt_exp=dt)
        else
          call correctTimeStep(itime,ierr)
        endif

c     Assign old time solution

        call equateDerivedType(u_nm,u_n)
        call equateDerivedType(u_n ,u_np)

c     Time update

#if defined(petsc)
        !Begin profiling second stage
        call PetscLogStagePop(perr)
        call PetscLogStagePush(prof_stage(2),perr)
#endif

        time   = time  + dt

        call timeStep(u_n,u_np,ierr)

c     Check for error in time stepping

#if defined(petsc)
        !Begin profiling third stage
        call PetscLogStagePop(perr)
        call PetscLogStagePush(prof_stage(3),perr)
#endif

        if (ierr == 1) then     !Restart time step

          call equateDerivedType(u_np,u_n)

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

        if (postprocess) call postProcessSol(u_np)

c     Update counters (only if timeStep is successful)

        tmrst  = tmrst + dt
        nrst   = nrst  + 1

c     Output per time step

        call output(u_np)

c     Time level data dump

#if defined(RFX)
        if (nrst.eq.ndstep.or.tmrst.ge.dstep) then
#else
        if (nrst.eq.ndstep.or.tmrst.ge.0.99*dstep) then
#endif
          nrst  = 0
          if (itime.gt.0) tmrst = tmrst - dstep
          call writeRecordFile(urecord,itime,time,dt,u_np)
        endif

        if (eigen) exit

      enddo       !End of time loop

c Average explicit time step

      if (cnfactor.eq.1d0 .and. my_rank == 0) then
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

#if defined(petsc)
      !Finalize profiling
      call PetscLogStagePop(perr)

      call PetscFinalize(perr)
#endif

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

      implicit none

c Call variables

      integer :: ierr

      type(var_array),pointer :: vnp,vn

c Local variables

      real(8) :: x(ntotd)

c Interfaces

      INTERFACE
        subroutine evaluateNonlinearFunction(varray,fi)
        use parameters
        use variable_setup
        real(8)          :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluateNonlinearFunction
      END INTERFACE

      INTERFACE
        subroutine explicit(vn,vnp,iter,iout)
        use variable_setup
        integer          :: iter,iout
        type(var_array),pointer :: vnp,vn
        end subroutine explicit
      END INTERFACE

      INTERFACE
        subroutine applyPreconditioner(n,x,y,out)
        integer :: n,out
        real(8),dimension(n) :: x,y
        end subroutine
      END INTERFACE

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

c Begin program

      ierr = 0

      itwhistler = 0

c Evaluate nonlinear function at old time for theta scheme

cc      jit = -1
      old_time_eval = .true.     !This informs nlfunction that we are processing
                                 !  n time level info

      call evaluateNonlinearFunction(vnp,fold)

      old_time_eval = .false.

c Eigenvalue solver

      if (eigen) then

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
     .                  ,matrixFreeMatVec
     .                  ,applyPreconditioner
     .                  ,dot=dot2)

        select case(ierr)
        case(1)
          call pstop('lehoucq_eig'
     .              ,'Matrix not positive definite')
        end select

        !Output
        if (my_rank == 0) then
           write (*,*)
           write (*,*) '**************************************'
           write (*,'(a,1p,1e10.3,a)') ' ****** Eigenvalue is='
     $             ,1d0/dt-eigen_conf%eigv,' ******' 
           write (*,*) '**************************************'
        endif

        itgmres = eigen_conf%iter_out

        !Map eigenvector for plotting
        x = x + nk_conf%xk
        call mapVectorToStructure(vnp,x)

        !Free memory
        deallocate(nk_conf%xk,nk_conf%rk)
 
c Implicit update (Newton)

      elseif (cnfactor.lt.1d0) then

c     Initial guess

        if (predictor) then
          call findGuess(vn,vnp)
          call mapStructureToVector(x,vnp)
        else
          if (bdf2) then
            call AXPYDerivedType(-cn/cnp,vn,-cnm/cnp,u_nm,vnp)
            call mapStructureToVector(x,vnp)
          else
            call mapStructureToVector(x,vnp)
          endif
        endif

c     Newton-Krylov iteration

        nk_conf%etak_meth      = method
        nk_conf%ksmax          = maxksp
        nk_conf%gmmax          = maxitgm
        nk_conf%nwt_max_it_acc = maxitnwt/2
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

        call nk(neqd,ntotd,x,iguess,ilevel,ierr)

        itgmres = nk_conf%gm_it_out
        itnewt  = nk_conf%nwt_it_out

c     Map Newton solution to vnp (if no error)

        if (ierr.eq.0.or.ierr.eq.2) then
          call mapVectorToStructure(vnp,x)
          call imposeBoundaryConditions(vnp,1,1,1)
          if (predictor) call storeTSinfo(vn,vnp)
        endif

      else

c Explicit update

        itgmres = maxitgm

        call explicit(vn,vnp,itgmres,ilevel)

        call imposeBoundaryConditions(vnp,1,1,1)

        itnewt = 1

      endif

c Update counters

      gmres_tot = gmres_tot + itgmres
      newt_tot  = newt_tot  + itnewt
      wh_tot    = wh_tot    + itwhistler

c End program

      end subroutine timeStep


c module timeStepping
c ######################################################################
      module timeStepping

#if defined(samrai)
        use ts_io
#else
        use prj_int

        logical,private :: ts_after_restart = .false.

      contains

c     TS_chk_after_restart
c     ###############################################################
      function TS_chk_after_restart(time_slices) result(flag)

c     ---------------------------------------------------------------
c     Check for end of time loop
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      logical :: flag

      integer :: time_slices

c     Local variables

      integer :: counter

c     Begin program

      counter = itime-inewtime + 1

      flag = ts_after_restart.and.(counter <= time_slices)

c     End program

      end function TS_chk_after_restart

c     TS_correct_dt
c     ###############################################################
      subroutine TS_correct_dt(ierr)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer,INTENT(IN) :: ierr

c     Local variables

      integer :: itm

c     Begin program

c     Estimate local growth rate

      if (sm_flag == 1) gammat = TS_growth_rate(u_np)

c     Update counters

      call TS_update_preCounters(time)

c     Find new time step

      itm = itime-init_grid_offset

      call TS_find_new_dt(itm,ierr,dt)

c     Adjust to final time, if needed

      if (time + dt > tmax) then
        dt = tmax - time
        time = tmax
      else
        time = time + dt
      endif

c     End program

      end subroutine TS_correct_dt

c     TS_find_new_dt
c     ####################################################################
      subroutine TS_find_new_dt(itm,ierr,dt)

c     --------------------------------------------------------------------
c     Find new time step
c     --------------------------------------------------------------------

      use nk, ONLY: NOERR,CUTDT

      implicit none

c     Call variables

      real(8) :: dt
      integer,INTENT(IN) ::  ierr,itm

c     Local variables

#if defined(USER_PROVIDED_DT)
      external set_user_dt
#endif
      external set_explicit_dt

c     Begin program

      if (cnfactor > 5d-1) then  !Explicit
        call set_explicit_dt(dt)
  
        dtexp = dtexp + dt
      else
        call adapt_dt

#if defined(USER_PROVIDED_DT)
        call set_user_dt(dt)  !External
#endif

        call calculate_cnfactor
      endif

c     Reset restart flags

      if (ts_restart) ts_after_restart = .true.

      ts_restart = .false.

c     End program

      contains

c     adapt_dt
c     #######################################################################
      subroutine adapt_dt

        real(8) ::    coef1=0.8   !Time subcycling coefficient
     .               ,coef2=1.05  !Time recovery   coefficient

c$$$        if (ts_ramp_up_dt) then
c$$$          coef1 = ts_ramp_up_dt_fraction
c$$$          coef2 = ts_ramp_up_dt_recovery
c$$$        endif

        if (ierr /= NOERR) then   !Error in time stepping

          dt = dt*coef1

          if (dt < 1d-6*dtbase) then
            call pstop('adapt_dt','Time step too small')
          endif

          if (my_rank == 0 .and. ierr == CUTDT) write (*,240)
          if (my_rank == 0) write (*,400) dt

        else                 !No error in time stepping

          select case (sm_flag)
          case(1)  !Rannacher time stepping

            if (ts_ramp_up_dt) then
              if (itm == 1) then
                dt = dtbase*ts_ramp_up_dt_fraction
              else
                dt = min(dtbase,dt*ts_ramp_up_dt_recovery)
                ts_ramp_up_dt = (dt /= dtbase) !Stop dt ramp-up
              endif
            elseif (itm == 1.and.ts_adapt_dt) then
              call set_explicit_dt(dt)
            elseif (   (itm <= sm_pass  .and.(.not.ts_adapt_dt))
     .             .or.(itm <= sm_pass+1.and.      ts_adapt_dt )) then
              dt = 0.5*dtbase
            elseif((itm == sm_pass+1.and.(.not.ts_adapt_dt))
     .         .or.(itm == sm_pass+2.and.      ts_adapt_dt )) then
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
            endif

          case(2)  !BDF time stepping

            if (ts_ramp_up_dt) then
              if (itm == 1) then
                dt = dtbase*ts_ramp_up_dt_fraction
              else
                dt = min(dtbase,dt*ts_ramp_up_dt_recovery)
                ts_ramp_up_dt = (dt /= dtbase) !Stop dt ramp-up
              endif
            elseif (ts_restart) then
              dt = dt*coef2**(-15) !Quiet BDF2 start to ramp up in a few time steps
            elseif (ts_adapt_dt) then
              if (itm == 1) then
                call set_explicit_dt(dt)
              elseif (itm == 2) then
                dt = coef1*dtbase
c$$$              elseif (itm == 3) then
c$$$                dt = dtbase/coef2**15  !Quiet BDF2 start to ramp up in a few time steps
              else
                dt = min(dtbase,dt*coef2)
              endif
            else
              if  (itm <= bdf2_init_offset) then
                dt = 0.1*dtbase
                dtold = dt
              elseif (itm == 1+bdf2_init_offset) then
                dt = dtbase - bdf2_init_offset*dt
              elseif (itm == 2+bdf2_init_offset) then
                dt = dtbase
              else
                dt = min(dtbase,dt*coef2)
              endif
            endif

          case default   !Theta time stepping

            if (ts_ramp_up_dt) then
              if (itm == 1) then
                dt = dtbase*ts_ramp_up_dt_fraction
              else
                dt = min(dtbase,dt*ts_ramp_up_dt_recovery)
                ts_ramp_up_dt = (dt /= dtbase) !Stop dt ramp-up
              endif
            elseif (ts_adapt_dt) then
              if (itm == 1) then
                call set_explicit_dt(dt)
              elseif (itm == 2) then
                if (dt < dtbase) then
                  dt = dtbase - dt
                else
                  dt = dtbase
                endif
              elseif (itm == 3) then
                dt = dtbase
              else
                dt = min(dtbase,dt*coef2)
              endif
            else
              if (itm == 1) then
                dt = dtbase
              else
                dt = min(dtbase,dt*coef2)
              endif
            endif

          end select

        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step: ',1p,e11.4)

      end subroutine adapt_dt

c     calculate_cnfactor
c     #######################################################################
      subroutine calculate_cnfactor

      if (sm_flag == 1) then      !Rannacher

        if (itm == 1.and.ts_adapt_dt) then
          cnfactor = .5
        elseif (    (itm <= sm_pass+1 .and.      ts_adapt_dt )
     .          .or.(itm <= sm_pass   .and.(.not.ts_adapt_dt)))then
          cnfactor = .3
        elseif (    (itm >= sm_pass+1.and.(.not.ts_adapt_dt))
     .          .or.(itm >= sm_pass+2.and.      ts_adapt_dt )) then
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      elseif (sm_flag == 2) then  !BDF2

        bdf2 = (itm > 1 .and. (.not.ts_adapt_dt) )
     .    .or. (itm > 2 .and.       ts_adapt_dt  .and. ierr == NOERR)

        bdf2 = bdf2.and.(.not.ts_restart)

        if (bdf2) then
          !All these are affected by factor 1/dt (done in evaluateNLF)
          cnp = (2*dt+dtold)/(dt+dtold)
          cn  = -(dt+dtold)/dtold
          cnm = -cnp-cn    !Sum should be zero identically

          cnfactor = 0d0
        else
          cnp = 1d0
          cn  =-1d0
          cnm = 0d0

          bdf1 = (itm == 1)

          if (ts_restart) cnfactor = 0.45d0
        endif

      else  !Theta scheme

        bdf1 = (cnfactor == 0d0)

      endif

      alpha = 1d0 - cnfactor

      end subroutine calculate_cnfactor

      end subroutine TS_find_new_dt

c     TS_update
c     ######################################################################
      function TS_update(vnp,iout) result(ierr)

c     ----------------------------------------------------------------------
c     Performs TS update
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr,iout

      type(var_array),pointer :: vnp

c     Local variables
      
      INTERFACE
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use var_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c     Begin program
      
#if defined(petsc)
        !Begin profiling NL solver
        call switch_profile_stage(iNL)
#endif

c     Perform TS update

      if (eigen) then
        !Eigenvalue solver
        ierr = eigenSolution (vnp,iout)
      elseif (prj_int_ts) then
        ierr = prj_int_update(vnp,iout)
      else
        !Time step
        ierr = updateSolution(vnp,iout)
      endif

c     End program

      end function TS_update
#endif

      end module timeStepping

c module ts_base
c ######################################################################
      module ts_base

        use ts_setup

        use var_io

        use nk, ONLY: NOERR,CUTDT

#if defined(MK)
        use mk
#endif

        integer :: EXIT_STATUS=FATAL !Hardwire to fatal termination


        integer :: itnewt
     .            ,itgmres
     .            ,nrstep

        real(8) :: rstep,dfreq

        integer :: nrst             = 0
     .            ,gmres_tot        = 0
     .            ,newt_tot         = 0
     .            ,init_grid_offset = 0
     .            ,bdf2_init_offset = 1

        real(8) :: tmrst = 0d0   !Dump time for restart
     .            ,quota = 1d30  !Max CPU time (in seconds)
     .            ,dtexp = 0d0   !Explicit time accumulator
     .            ,exec_time=0d0 !Execution time (in seconds)

        logical,private :: ts_after_restart = .false.

#if !defined(samrai)

cc        real(8),pointer,dimension(:,:) :: old_f

c$$$#if !defined(petsc_c)
c$$$        type(var_array),pointer :: dn,dnm,dnp
c$$$        real(8) :: l_np,l_n,et1,et2
c$$$
      contains

c$$$c     findGuess
c$$$c     ###############################################################
c$$$      subroutine findGuess(vn,vnp)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Correct time step
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      type(var_array),pointer :: vn,vnp
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$      if (itime <= 2) then
c$$$        !Previous time step
c$$$        vnp = vn
c$$$cc      elseif (itime == 2) then
c$$$cc        !Linear extrapolation
c$$$cc        l_np = dt/dt_n
c$$$cc        call AXPYDerivedType(1d0,vn,l_np,dn,vnp)
c$$$cc        dnp = dn   !Initialize recursion relation
c$$$      else
c$$$        !Quadratic extrapolation
c$$$        l_np = dt/dt_n
c$$$        et1 = l_np*(1d0+l_n*(1+l_np)/(1+l_n))
c$$$        et2 =-l_np*l_n**2*(1+l_np)/(1+l_n)
c$$$        call AXPYDerivedType(et1,dn,et2,dnm,dnp) !Recursion relation
c$$$        call AXPYDerivedType(1d0,vn,1d0,dnp,vnp) !vnp = vn + dnp
c$$$      endif
c$$$
c$$$cc      write (*,*) itime,l_np,l_n,et1,et2
c$$$
c$$$cc      if (itime <= 2) then
c$$$cc        !Previous time step
c$$$cc        vnp = vn
c$$$cccc      elseif (itime == 2) then
c$$$cccc        !FE predictor
c$$$cccc        dn = dt*old_f(:,1)
c$$$cccc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp)
c$$$cc      else
c$$$cc        !AB@ predictor
c$$$cc        old_f(:,2) = 0.5*dt*(3*old_f(:,1)-old_f(:,2))
c$$$cc        dn = old_f(:,2)
c$$$cc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp) !vnp = vn - dt/2*(3*f_n - f_n-1)
c$$$cc      endif
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine findGuess
c$$$
c$$$c     storeTSinfo
c$$$c     ###############################################################
c$$$      subroutine storeTSinfo(vn,vnp)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Advance recursion relation
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      type(var_array),pointer :: vn,vnp
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$      l_n = l_np
c$$$
c$$$      call AXPYDerivedType(1d0,vnp,-1d0,vn,dnp)
c$$$
c$$$      dnm = dn
c$$$      dn  = dnp
c$$$
c$$$cc      old_f(:,2) = old_f(:,1)
c$$$cc      call evaluateNonlinearFunction(vnp,old_f(:,1))
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine storeTSinfo

c     TS_set_dt_ramp_up
c     ###############################################################
      subroutine TS_set_dt_ramp_up(fraction,recovery)

c     ---------------------------------------------------------------
c     Set up initial dt ramp up
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: fraction,recovery

c     Local variables

c     Begin program

      ts_ramp_up_dt = .true.
      ts_ramp_up_dt_fraction = fraction
      ts_ramp_up_dt_recovery = recovery

c     End program

      end subroutine TS_set_dt_ramp_up

c     TS_set_bdf2_init_offset
c     ###############################################################
      subroutine TS_set_bdf2_init_offset(offset)

c     ---------------------------------------------------------------
c     Set initial offset for bdf2 integration
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: offset

c     Local variables

c     Begin program

      bdf2_init_offset = offset

c     End program

      end subroutine TS_set_bdf2_init_offset

c     TS_chk_end
c     ###############################################################
      function TS_chk_end() result(exit_tloop)

c     ---------------------------------------------------------------
c     Check for end of time loop
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      logical :: exit_tloop

c     Local variables

c     Begin program

      exit_tloop = (tmax.gt.0d0.and.time.ge.(tmax-1d-9*dt))
     .         .or.(numtime.ge.0.and.itime+1.ge.(numtime+inewtime))
     .         .or.(exec_time > quota)

      if (exit_tloop .and. exec_time <= quota) EXIT_STATUS = SUCCESS

      if (exec_time > quota) EXIT_STATUS = NON_FATAL
      
c     End program

      end function TS_chk_end

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

      time = time + dt

c     End program

      end subroutine TS_correct_dt

c     TS_update_postCounters
c     ####################################################################
      subroutine TS_update_postCounters

c     --------------------------------------------------------------------
c     Update time stepping counters
c     --------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      if (relax_init_grid) return

      tmrst  = tmrst + dt
      nrst   = nrst  + 1

      gmres_tot = gmres_tot + itgmres
      newt_tot  = newt_tot  + itnewt
      
#if defined(petsc)
      !Begin profiling postprocessing stage
      call switch_profile_stage(iPOST)
#endif

c     End program

      end subroutine TS_update_postCounters

c     TS_find_new_dt
c     ####################################################################
      subroutine TS_find_new_dt(itm,ierr,dt)

c     --------------------------------------------------------------------
c     Find new time step
c     --------------------------------------------------------------------

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

          if (ts_restart) cnfactor = 0.45d0
        endif

      else  !Theta scheme

        bdf1 = (cnfactor == 0d0)

      endif

      alpha = 1d0 - cnfactor

      end subroutine calculate_cnfactor

      end subroutine TS_find_new_dt

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

      integer, save :: counter = 0

c     Begin program

      counter = counter + 1

      flag = ts_after_restart.and.(counter <= time_slices)

c     End program

      end function TS_chk_after_restart

#endif

      end module ts_base

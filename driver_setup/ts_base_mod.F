c module ts_base
c ######################################################################
      module ts_base

        use ts_setup

        use var_io

#if defined(MK)
        use mk
#endif

        integer :: EXIT_STATUS=FATAL !Hardwire to fatal termination

        integer :: itnewt
     .            ,itgmres

        integer :: nrst             = 0
     .            ,gmres_tot        = 0
     .            ,newt_tot         = 0
     .            ,init_grid_offset = 0
     .            ,bdf2_init_offset = 1

        real(8) :: tmrst = 0d0   !Dump time for restart
     .            ,quota = 1d30  !Max CPU time (in seconds)
     .            ,dtexp = 0d0   !Explicit time accumulator

#if !defined(samrai)
      contains

c     findGuess
c     ###############################################################
      subroutine findGuess(vnm,vn,vnp)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: vnm,vn,vnp

c     Local variables

c$$$      real(8),pointer,dimension(:,:) :: old_f
c$$$
c$$$      type(var_array),pointer :: dn,dnm,dnp
c$$$      real(8) :: l_np,l_n,et1,et2

#if defined(EXT_NEWT_GUESS)
      INTERFACE
        function app_newt_guess(vnm,vn,vnp)
        use var_setup
        type(var_array),pointer :: vnm,vn,vnp
        end function app_newt_guess
      END INTERFACE
#endif

c     Begin program

#if defined(EXT_NEWT_GUESS)
      call app_newt_guess(vnm,vn,vnp)
#else
#if !defined(NO_BDF2_XOLD)
      call AXPYDerivedType(-cn/cnp,vn,-cnm/cnp,vnm,vnp)
#endif
 
c$$$      if (itime < 2) then
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

c$$$      write (*,*) itime,l_np,l_n,et1,et2
c$$$
c$$$      if (itime <= 2) then
c$$$        !Previous time step
c$$$        vnp = vn
c$$$cc      elseif (itime == 2) then
c$$$cc        !FE predictor
c$$$cc        dn = dt*old_f(:,1)
c$$$cc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp)
c$$$      else
c$$$        !AB@ predictor
c$$$        old_f(:,2) = 0.5*dt*(3*old_f(:,1)-old_f(:,2))
c$$$        dn = old_f(:,2)
c$$$        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp) !vnp = vn - dt/2*(3*f_n - f_n-1)
c$$$      endif
#endif

c     End program

      end subroutine findGuess

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

      real(8) :: e_time
      
c     Begin program

#if defined(petsc)
      call PetscGetTime   (prof_time(3),prof_err)
      e_time = prof_time(3)-prof_time(1)
cc      if (my_rank == 0) write (*,*) "Time counters",e_time, quota
#else
      e_time = 0d0
#endif

      exit_tloop= (tmax/=0d0.and.time>=(tmax-1d-13*itime*dt).and.dt>0d0)
     .        .or.(              time<=(tmax-1d-13*itime*dt).and.dt<0d0)
     .        .or.(numtime.ge.0.and.itime+1.ge.(numtime+inewtime))
     .        .or.(e_time > quota)

      if (exit_tloop .and. e_time <= quota) EXIT_STATUS = SUCCESS

      if (e_time > quota) EXIT_STATUS = NON_FATAL
      
c     End program

      end function TS_chk_end

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

#endif

      end module ts_base

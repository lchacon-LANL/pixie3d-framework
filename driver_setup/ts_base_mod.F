c module ts_base
c ######################################################################
      module ts_base

        use ts_setup

        use var_io

        use nk, ONLY: NOERR,CUTDT

#if defined(MK)
        use mk
#endif

        integer :: nrstep
     .            ,nrst
        integer :: itnewt
     .            ,itgmres
     .            ,gmres_tot
     .            ,newt_tot
     .            ,init_grid_offset=0
     .            ,bdf2_init_offset=1

        real(8) :: tmrst       !Dump time for restart
     .            ,quota=1d30  !Max CPU time (in seconds)
        real(8) :: rstep,dfreq,dtexp

#if !defined(samrai)

cc        real(8),pointer,dimension(:,:) :: old_f

c$$$#if !defined(petsc_c)
c$$$        type(var_array),pointer :: dn,dnm,dnp
c$$$        real(8) :: l_np,l_n,et1,et2
c$$$
      contains
c$$$
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

c     TS_set_bdf2_init_offset
c     ###############################################################
      subroutine TS_set_bdf2_init_offset(offset)

c     ---------------------------------------------------------------
c     Check for end of time loop
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

      if (sm_flag == 1) call calculate_gammat

c     Update counters

      call TS_update_timecounters(time)

c     Find new time step

      itm = itime-init_grid_offset

      call findNewTimeStep(itm,ierr,dt)

      time = time + dt

c     End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

c     Calculation of local growth rate for CN

        integer :: ieq

        real(8) :: dmag1,dmag2,dpert,mag(neqd)

        real(8) :: array(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

        do ieq=1,neqd

          array = (u_np%array_var(ieq)%array
     .            -u_0 %array_var(ieq)%array )
          array = array*array

          dpert = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)

          dpert = sqrt(dpert)

          array = (u_n%array_var(ieq)%array
     .            -u_0%array_var(ieq)%array )
          array = array*array

          dmag1 = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)

          array = (u_np%array_var(ieq)%array
     .            +u_n %array_var(ieq)%array
     .         -2.*u_0 %array_var(ieq)%array )
          array = array*array

          dmag2 = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)

          if (dpert/=sqrt(dmag1).and.dmag2 /=0d0) then
            mag(ieq) = .5*dt*sqrt(dmag2)/(dpert-sqrt(dmag1))
          else
            mag(ieq) = 1e30
          endif

        enddo

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

      end subroutine TS_correct_dt

c$$$#else
c$$$
c$$$      contains
c$$$
c$$$c     correctTimeStepPetsc
c$$$c     ######################################################################
c$$$      subroutine correctTimeStepPetsc(dn,dnh,dnp,ierr,dt_to_c)
c$$$
c$$$      use parameters
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer     :: ierr
c$$$      real(8)     :: dn(neqd),dnh(neqd),dnp(neqd),dt_to_c
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$c     Estimate local growth rate
c$$$
c$$$      if (sm_flag == 1) call calculate_gammat
c$$$
c$$$c     Find new time step
c$$$
c$$$      call findNewTimeStep(itime+1,ierr)
c$$$
c$$$c     Update counters (only if timeStep is successful)
c$$$
c$$$      itime  = itime + 1
c$$$
c$$$      time   = time  + dt
c$$$      tmrst  = tmrst + dt
c$$$      nrst   = nrst  + 1
c$$$
c$$$      dt_to_c = dt
c$$$
c$$$c     End program
c$$$
c$$$      contains
c$$$
c$$$c     calculate_gammat
c$$$c     #######################################################################
c$$$      subroutine calculate_gammat
c$$$
c$$$c     -----------------------------------------------------------------------
c$$$c     Calculation of local growth rate for CN
c$$$c    ------------------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$c     Local variables
c$$$
c$$$        integer     :: ieq
c$$$
c$$$        real(8)     :: mag(neqd)
c$$$
c$$$c     Begin program
c$$$
c$$$        where (dnp /= dn .and. dnh /= 0d0)
c$$$          mag = dt*dnh/(dnp-dn)
c$$$        elsewhere
c$$$          mag = 1e30
c$$$        end where
c$$$
c$$$        gammat = 1./minval(abs(mag))
c$$$
c$$$      end subroutine calculate_gammat
c$$$
c$$$      end subroutine correctTimeStepPetsc
c$$$
c$$$#endif

c     TS_init_otherCounters
c     ####################################################################
      subroutine TS_init_otherCounters

c     --------------------------------------------------------------------
c     Initialize time stepping counters
c     --------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      gmres_tot = 0
      newt_tot  = 0

      itime     = inewtime - 1

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0
      
      if (relax_init_grid) init_grid_offset = 1

c     PETSc CPU counters and profiling
      
#if defined(petsc)
      call PetscGetTime   (prof_time(1),prof_err)
      call PetscGetCPUTime(prof_time(2),prof_err)

      !Begin profiling postprocessing stage
      call switch_profile_stage(iPOST)
#endif

c     End program

      end subroutine TS_init_otherCounters

c     TS_update_otherCounters
c     ####################################################################
      subroutine TS_update_otherCounters

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

      end subroutine TS_update_otherCounters

c     findNewTimeStep
c     ####################################################################
      subroutine findNewTimeStep(itm,ierr,dt)

c     --------------------------------------------------------------------
c     Find new time step
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dt
      integer,INTENT(IN) ::  ierr,itm

c     Local variables

      external findExplicitDt

c     Begin program

      if (    cnfactor > 5d-1                   !Explicit
     .   .or.(timecorr.and.ierr==0.and.itm==1)  !First explicit step if timecorr
     .   ) then
        call findExplicitDt(dt)
  
        if (cnfactor > 5d-1) dtexp = dtexp + dt
      else
        call adapt_dt

        call calculate_cnfactor
      endif

c     End program

      contains

c     adapt_dt
c     #######################################################################
      subroutine adapt_dt

        real(8) ::    coef1=0.8   !Time subcycling coefficient
     .               ,coef2=1.05  !Time recovery   coefficient

        if (ierr > 0) then   !Error in time stepping

          dt = dt*coef1
          if (dt < 1d-3*dtbase) then
            call pstop('adapt_dt','Time step too small')
          endif
          if (my_rank == 0 .and. ierr == CUTDT) write (*,240)
          if (my_rank == 0) write (*,400) dt

        else                 !No error in time stepping

          select case (sm_flag)
          case(1)  !Rannacher time stepping

            if (   (itm <= sm_pass  .and.(.not.timecorr))
     .         .or.(itm <= sm_pass+1.and.      timecorr )) then
              dt = 0.5*dtbase
            elseif ((itm == sm_pass+1.and.(.not.timecorr))
     .          .or.(itm == sm_pass+2.and.      timecorr )) then
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
            endif

          case(2)  !BDF time stepping

            if (timecorr) then
              if  (itm == 2) then
                dt = 0.1*dtbase
              elseif (itm == 3) then
                dt = dtbase/coef2**15  !Quiet BDF2 start to ramp up in a few time steps
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

          case default   !Default time stepping

            if     ((itm == 1.and.(.not.timecorr))) then
              dt = dtbase
            elseif ((itm == 2.and.      timecorr )) then
              if (dt < dtbase) then
                dt = dtbase - dt
              else
                dt = dtbase
              endif
            elseif ((itm == 3.and.      timecorr )) then
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
            endif

          end select

        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step: ',1p,e10.4)

      end subroutine adapt_dt

c     calculate_cnfactor
c     #######################################################################
      subroutine calculate_cnfactor

      if (sm_flag == 1) then      !Rannacher

        if (itm == 1.and.timecorr) then
          cnfactor = .5
        elseif (    (itm <= sm_pass+1 .and.      timecorr )
     .          .or.(itm <= sm_pass   .and.(.not.timecorr)))then
          cnfactor = .3
        elseif (    (itm >= sm_pass+1.and.(.not.timecorr))
     .          .or.(itm >= sm_pass+2.and.      timecorr )) then
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      elseif (sm_flag == 2) then  !BDF2

        cnfactor = 0d0

        bdf2 = (itm > 1 .and. (.not.timecorr) )
     .    .or. (itm > 2 .and.       timecorr  .and. ierr == NOERR)

        bdf2 = bdf2.and.(.not.restart)

        if (bdf2) then
          !All these are affected by factor 1/dt (done in evaluateNLF)
          cnp = (2*dt+dtold)/(dt+dtold)
          cn  = -(dt+dtold)/dtold
          cnm = -cnp-cn    !Sum should be zero identically
        else
          cnp = 1d0
          cn  =-1d0
          cnm = 0d0

          restart = .false.
        endif

      else  !Theta scheme

        bdf1 = (cnfactor == 0d0)

      endif

      alpha = 1d0 - cnfactor

      end subroutine calculate_cnfactor

      end subroutine findNewTimeStep
#endif

      end module ts_base

c module ts_setup
c ######################################################################
      module ts_setup

        use parameters

        use ts_solver_setup

        use io

        use OMP

        !Default TS variables
        real(8) :: cnfactor=0.5        ! Crank-Nicolson factor (implicit if <= 0.5)
     .            ,dstep=0d0           ! Time interval between plots (if zero,
                                       !        dstep is calculated in code)
     .            ,eig_dt=1d0          ! Time step for eigenvalue dynamical system
     .            ,dt      = 1d0       ! Time step (if zero, dt is calculated in code)
     .            ,tmax    =-1d0       ! Target time, in Alfven times (not active if negative)
     .            ,trstrt  =-1d0       ! Restart time (not active if negative)

        integer :: sm_pass=2           ! Number of initial smoother passes for Rannacher TS
     .            ,sm_flag=0           ! Time smoothing flag:
                                       !   0 -> Theta scheme
                                       !   1 -> Rannacher time stepping
                                       !   2 -> BDF2
     .            ,numtime=-1          ! Number of time steps (not active if negative)
     .            ,ndstep=0            ! # time steps between plots (if zero,
                                       !        ndstep is calculated in code)
     .            ,nrstrt=-1           ! Restart cycle number (not active if negative)

        logical :: ts_adapt_dt=.false. ! Time adaptive algorithm (based on Newton convergence)
     .            ,ts_restart=.false.  ! Restarting flag
     .            ,source=.true.       ! Impose source to provide equilibrium
     .            ,predictor=.false.   ! Whether to extrapolate in time to find first guess (not active)
     .            ,postprocess=.false. ! Whether to run solution postprocessing after time step
     .            ,preprocess=.false.  ! Whether to run solution preprocessing before time step
     .            ,eigen=.false.       ! Whether to run framework in eigensolver mode
     .            ,perturb_rstrt=.false.! Whether to apply perturbation after restart

        logical, private :: early_restart !Whether restart has happened from intermediate time
        logical :: ts_ramp_up_dt = .false.       !Whether to ramp up dt initially
        real(8) :: ts_ramp_up_dt_fraction = 5d-1 !Fraction of initial timestep
     .            ,ts_ramp_up_dt_recovery = 1.05 !Ramp up coefficient

        !TS working variables
        real(8),allocatable,dimension(:) :: cnf,one_over_dt
     .                                     ,bdfp,bdfn,bdfnm

#if !defined(samrai)
        real(8),pointer,dimension(:) :: fold,fsrc,x_old,x0
#endif

        integer :: ntimemax,itime=0,inewtime

        real(8) :: dt_nm,dt_n,dtold,time,dtbase,tn,tnm,tnm2,alpha,gammat

        real(8) :: cnp=1d0,cn=-1d0,cnm=0d0

        logical :: chk_src=.false.
     .            ,source_eval=.false.
     .            ,old_time_eval=.false.
     .            ,bdf2=.false.
     .            ,bdf1=.false.
     .            ,relax_init_grid=.false.

      contains

c     setEarlyTimeRestartFlag
c     ###############################################################
      subroutine setEarlyTimeRestartFlag(flag)

c     ---------------------------------------------------------------
c     Performs update of time counters
c     ---------------------------------------------------------------

      implicit none

      logical :: flag

c     Begin program

      early_restart = flag

c     End program

      end subroutine setEarlyTimeRestartFlag

c     getEarlyTimeRestartFlag
c     ###############################################################
      function getEarlyTimeRestartFlag() result(flag)

c     ---------------------------------------------------------------
c     Performs update of time counters
c     ---------------------------------------------------------------

      implicit none

      logical :: flag

c     Begin program

      flag = early_restart

c     End program

      end function getEarlyTimeRestartFlag

c     TS_init_timecounters
c     ###############################################################
      subroutine TS_init_timecounters()

c     ---------------------------------------------------------------
c     Performs update of time counters
c     ---------------------------------------------------------------

      implicit none

c     Begin program

      itime = 0

      inewtime = 1

      tnm2 = 0d0
      tnm  = 0d0
      tn   = 0d0
      time = 0d0

      dtold = dt
      dt_n  = dt
      dt_nm = dt

c     End program

      end subroutine TS_init_timecounters

c     TS_rollback_timecounters
c     ###############################################################
      subroutine TS_rollback_timecounters(exclude_itime)

c     ---------------------------------------------------------------
c     Performs rollback of time counters when subcycling/restarting
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      logical,optional :: exclude_itime

c     Local variables

      logical :: ex_it

c     Begin program

      if (PRESENT(exclude_itime)) then
        ex_it = exclude_itime
      else
        ex_it = .false.
      endif

      if (.not.ex_it) itime = itime - 1

      dtold = tn - tnm
      dt_n  = dtold
      dt_nm = 0d0

      time  = tn
      tn    = tnm
      tnm   = tnm2
      tnm2  = 0d0

c     End program

      end subroutine TS_rollback_timecounters

c     TS_update_timecounters
c     ###############################################################
      subroutine TS_update_timecounters(time,exclude_itime)

c     ---------------------------------------------------------------
c     Performs update of time counters
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: time

      logical,optional :: exclude_itime

c     Local variables

      logical :: ex_it

c     Begin program

      if (PRESENT(exclude_itime)) then
        ex_it = exclude_itime
      else
        ex_it = .false.
      endif

      if (.not.ex_it) itime = itime + 1

      tnm2= tnm
      tnm = tn
      tn  = time

      dtold = tn - tnm
      if (dtold == 0d0) dtold = dt

      dt_nm = dt_n
      dt_n  = dtold

c     End program

      end subroutine TS_update_timecounters

#if !defined(samrai)
c     TS_create
c     #################################################################
      subroutine TS_create

c     -----------------------------------------------------------------
c     Allocates time-stepping arrays
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      allocate(fold(ntotd),fsrc(ntotd),x_old(ntotd),x0(ntotd))

      allocate(cnf(neqd),one_over_dt(neqd)
     .        ,bdfp(neqd),bdfn(neqd),bdfnm(neqd))

c     Initialization

      bdfp  =  1d0
      bdfn  = -1d0
      bdfnm =  0d0

      fold  = 0d0
      fsrc  = 0d0
      x_old = 0d0
      x0    = 0d0

      call TS_init_timecounters

      alpha = 1. - cnfactor

      dtbase = dt

      !Time smoothing
      if (sm_flag == 0) then
        sm_pass= 0
        if (cnfactor > 1d0 .or. cnfactor < 0d0) then
          call pstop('createTS','cnfactor out of range')
        endif
      endif

      call ts_solver_consistency_chk

c     End program

      end subroutine TS_create

c     TS_destroy
c     #################################################################
      subroutine TS_destroy

c     -----------------------------------------------------------------
c     Deallocates time-stepping arrays
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      deallocate(fold,fsrc,x_old,x0)

      deallocate(cnf,one_over_dt,bdfp,bdfn,bdfnm)

c     End program

      end subroutine TS_destroy
#endif

c     readTSInput
c     #################################################################
      function readTSInput(file,dump) result(ioerr)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ioerr
      character(*) :: file
      logical :: dump

c     Local variables

      integer :: uinput=1234

c     Namelists

      !Temporal discretization
      namelist /ts/ dt,tmax,dstep,numtime,ndstep,source,ts_restart
     .             ,ts_adapt_dt
     .             ,cnfactor,sm_flag,sm_pass,predictor,postprocess
     .             ,preprocess,eigen,eig_dt,perturb_rstrt
     .             ,ts_ramp_up_dt,ts_ramp_up_dt_fraction
     .             ,ts_ramp_up_dt_recovery

c     Begin program

c     Read initialization parameters

      do while(is_open_unit(uinput))
        uinput = uinput + 1
      enddo

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,ts,iostat=ioerr)

      close(unit=uinput)

cc      if (ioerr/=0) call pstop('readGridInput','Problem reading "ts"')
      if (ioerr/=0) return

c     Consistency check

      if (eigen) then
        source  = .false.
        cnfactor = 0d0
        sm_flag  = 0
      endif

#if defined(pit)
      postprocess = .true.
#endif

c     Write initialization parameters to standard output

      if (my_rank == 0.and.dump) then
        write (*,*) '------------------'
        write (*,*) 'TS configuration'
        write (*,*) '------------------'
        write (*,ts)
      endif

c     End program

      end function readTSInput

c     TS_growth_rate
c     #######################################################################
      function TS_growth_rate(varray,d_pert) result(gammat)

        use var_setup

        implicit none

        type(var_array),pointer :: varray

        real(8) :: gammat

        real(8),optional :: d_pert(neqd)

c     Calculation of local growth rate for CN

        integer :: ieq

        real(8) :: dpert,dpert_n,dpert_nh,mag(neqd)

        real(8) :: array(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

        do ieq=1,neqd

          array = (varray%array_var(ieq)%array
     .            -u_0   %array_var(ieq)%array )
          array = array*array

          dpert = integral(gv%gparams,1,array)

          dpert = sqrt(dpert)

          if (PRESENT(d_pert)) d_pert(ieq) = dpert

          array = (u_n%array_var(ieq)%array
     .            -u_0%array_var(ieq)%array )
          array = array*array

          dpert_n = integral(gv%gparams,1,array)

          dpert_n = sqrt(dpert_n)

c$$$          array = (varray%array_var(ieq)%array
c$$$     .            +u_n   %array_var(ieq)%array
c$$$     .         -2.*u_0   %array_var(ieq)%array )
c$$$          array = array*array
c$$$
c$$$          dpert_nh = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)
c$$$
c$$$          dpert_nh = sqrt(dpert_nh)

          if (dpert > 0d0 .and. dpert_n > 0d0 .and. time > tn) then
            mag(ieq) = (log(dpert)-log(dpert_n))/(time-tn)
c$$$          if (dpert /= dpert_n.and.dpert_nh /=0d0.and.time > tn) then
c$$$            mag(ieq) = .5*(time-tn)*dpert_nh/(dpert-dpert_n)
          else
            mag(ieq) = 0d0
          endif

        enddo

c$$$        gammat = 1d0/minval(abs(mag))
        gammat = maxval(abs(mag))

      end function TS_growth_rate

      end module ts_setup

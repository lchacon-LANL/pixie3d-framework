c module ts_setup
c ######################################################################
      module ts_setup

        use parameters

        use nk_setup

        use io

        !Default TS variables
        real(8) :: cnfactor=0.5        ! Crank-Nicolson factor (implicit if <= 0.5)
     .            ,dstep=0d0           ! Time interval between plots (if zero,
                                       !        dstep is calculated in code)
     .            ,eig_dt=1d0          ! Time step for eigenvalue dynamical system
     .            ,dt      = 1d0       ! Time step (if zero, dt is calculated in code)
     .            ,tmax    =-1d0       ! Target time, in Alfven times (not active if negative)

        integer :: sm_pass=2           ! Number of initial smoother passes for Rannacher TS
     .            ,sm_flag=0           ! Time smoothing flag:
                                       !   0 -> Theta scheme
                                       !   1 -> Rannacher time stepping
                                       !   2 -> BDF2
     .            ,numtime=-1          ! Number of time steps
     .            ,ndstep=0            ! # time steps between plots (if zero,
                                       !        ndstep is calculated in code)

        logical :: timecorr=.true.     ! Time adaptive algorithm (based on Newton convergence)
     .            ,restart=.false.     ! Restarting flag
     .            ,source=.true.       ! Impose source to provide equilibrium
     .            ,predictor=.false.   ! Whether to extrapolate in time to find first guess (not active)
     .            ,postprocess=.false. ! Whether to run solution postprocessing after time step
     .            ,preprocess=.false.  ! Whether to run solution preprocessing before time step
     .            ,eigen=.false.       ! Whether to run framework in eigensolver mode
     .            ,perturb_rstrt=.false.! Whether to apply perturbation after restart

        !TS working variables
        real(8),allocatable,dimension(:) :: cnf,one_over_dt
     .                                     ,bdfp,bdfn,bdfnm

#if !defined(samrai)
        real(8),pointer,dimension(:) :: fold,fsrc,x_old,x0
#endif

        integer :: ntimemax,itime=0,inewtime

        real(8) :: dt_n,dtold,time,dtbase,tn,tnm,alpha,gammat

        real(8) :: cnp=1d0,cn=-1d0,cnm=0d0

        logical :: chk_src=.false.
     .            ,source_eval=.false.
     .            ,old_time_eval=.false.
     .            ,bdf2=.false.
     .            ,bdf1=.false.
     .            ,relax_init_grid=.false.


      contains

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

      tnm  = 0d0
      tn   = 0d0
      time = 0d0

      dtold = dt
      dt_n  = dt

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

      time  = tn
      tn    = tnm
      tnm   = 0d0

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

      tnm = tn
      tn  = time

      dtold = tn - tnm
      if (dtold == 0d0) dtold = dt

      dt_n  = dtold

c     End program

      end subroutine TS_update_timecounters

#if !defined(samrai)
c     createTS
c     #################################################################
      subroutine createTS

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

      call TS_init_timecounters()

      alpha = 1. - cnfactor

      dtbase = dt

      !Time smoothing
      if (sm_flag == 0) then
        sm_pass= 0
        if (cnfactor > 1d0 .or. cnfactor < 0d0) then
          call pstop('createTS','cnfactor out of range')
        endif
      endif

c     End program

      end subroutine createTS

c     destroyTS
c     #################################################################
      subroutine destroyTS

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

      end subroutine destroyTS
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
      namelist /ts/ dt,tmax,dstep,numtime,ndstep,source,restart,timecorr
     .                ,cnfactor,sm_flag,sm_pass,predictor,postprocess
     .                ,preprocess,eigen,eig_dt,perturb_rstrt

c     Begin program

c     Set defaults (set in module)

cc      dt       = 5.            ! Time step (if zero, dt is calculated in code)
cc      tmax     = -1d0          ! Target time, in Alfven times.
cc      numtime  = -1            ! Number of time steps
cc      ndstep   = 0             ! # time steps between plots (if zero,
cc                               !        ndstep is calculated in code)
cc      dstep    = 0.            ! Time interval between plots (if zero,
cc                               !        dstep is calculated in code)
cc
cc      restart  = .false.       ! Restarting flag
cc      perturb_rstrt = .false.! Whether to apply perturbation after restart
cc
cc      timecorr = .true.        ! Time adaptive algorithm (based on Newton convergence)
cc      postprocess = .false.    ! Whether to run solution postprocessing after time step
cc      preprocess = .false.     ! Whether to run solution preprocessing before time step
cc      source   = .true.        ! Impose source to provide equilibrium
cc      predictor= .false.       ! Whether to extrapolate in time to find first guess (not active)
cc
cc      eigen    = .false.       ! Whether to run framework in eigensolver mode
cc      eig_dt   = 1d0           ! Time step for eigenvalue dynamical system
cc
cc      cnfactor = 0.5           ! Crank-Nicolson factor (implicit if <= 0.5)
cc      sm_flag  = 0             ! Time smoothing flag:
cc                               !   0 -> Theta scheme
cc                               !   1 -> Rannacher time stepping
cc                               !   2 -> BDF2
cc      sm_pass  = 2             ! Number of initial smoother passes for Rannacher TS

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

      end module ts_setup

c module ts_setup
c ######################################################################
      module ts_setup

        use parameters

        use nk_setup

        real(8),allocatable,dimension(:) :: cnf,one_over_dt
     .                                     ,bdfp,bdfn,bdfnm

#if !defined(samrai)
        real(8),pointer,dimension(:)   :: fold,fsrc,x_old
#endif
        real(8),pointer,dimension(:,:) :: old_f

        integer    :: numtime,ndstep,nrstep,inewtime,nrst,ntimemax

        real(8)    :: tmrst,dt,dtbase,cnfactor,alpha,gammat,dt_n
        real(8)    :: tmax,dstep,rstep,time,tn,tnm
     $               ,dfreq,dtexp,dtold,eig_dt

        !Initialized variables
        integer    :: itime=0
     .               ,sm_pass=2
     .               ,sm_flag=0

        logical    :: timecorr=.false.
     .               ,restart=.false.
     .               ,source=.true.
     .               ,predictor=.false.
     .               ,postprocess=.false.
     .               ,preprocess=.false.
     .               ,chk_src=.false.
     .               ,source_eval=.false.
     .               ,old_time_eval=.false.
     .               ,bdf2=.false.
     .               ,bdf1=.false.
     .               ,eigen=.false.
     .               ,relax_init_grid=.false.

        real(8)    :: cnp=1d0,cn=-1d0,cnm=0d0

      contains

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

      allocate(fold(ntotd),fsrc(ntotd),x_old(ntotd))

      allocate(cnf(neqd),one_over_dt(neqd)
     .        ,bdfp(neqd),bdfn(neqd),bdfnm(neqd))

c     Initialization

      bdfp  =  1d0
      bdfn  = -1d0
      bdfnm =  0d0

      fold = 0d0
      fsrc = 0d0
      x_old = 0d0

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

      deallocate(fold,fsrc,x_old)

      deallocate(cnf,one_over_dt,bdfp,bdfn,bdfnm)

c     End program

      end subroutine destroyTS
#endif

c     readTSInput
c     #################################################################
      subroutine readTSInput(file)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: file

c     Local variables

      integer :: ioerr,uinput=1234

c     Namelists

      !Temporal discretization
      namelist /ts/ dt,tmax,dstep,numtime,ndstep,source,restart,timecorr
     .                ,cnfactor,sm_flag,sm_pass,predictor,postprocess
     .                ,eigen,eig_dt

c     Begin program

c     Set defaults

      dt       = 5.            ! Time step (if zero, dt is calculated in code)
      tmax     = -1d0          ! Target time, in Alfven times.
      numtime  = -1            ! Number of time steps
      ndstep   = 0             ! # time steps between plots (if zero,
                               !        ndstep is calculated in code)
      dstep    = 0.            ! Time interval between plots (if zero,
                               !        dstep is calculated in code)

      restart  = .false.       ! Restarting flag
      timecorr = .true.        ! Time adaptive algorithm (based on Newton convergence)
      postprocess = .false.    ! Whether to run solution postprocessing after time step
      source   = .true.        ! Impose source to provide equilibrium
      predictor= .false.       ! Whether to extrapolate in time to find first guess (not active)

      eigen    = .false.       ! Whether to run PIXIE3D in eigensolver mode
      eig_dt   = 1d0           ! Time step for eigenvalue dynamical system

      cnfactor = 0.5           ! Crank-Nicolson factor (implicit if <= 0.5)
      sm_flag  = 0             ! Time smoothing flag:
                               !   0 -> Theta scheme
                               !   1 -> Rannacher time stepping
                               !   2 -> BDF2
      sm_pass  = 2             ! Number of initial smoother passes for Rannacher TS

c     Read initialization parameters

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,ts,iostat=ioerr)

      if (ioerr/=0) call pstop('readGridInput','Problem reading "ts"')

      close(unit=uinput)

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

      if (my_rank == 0) then
        write (*,*) '------------------'
        write (*,*) 'TS configuration'
        write (*,*) '------------------'
        write (*,ts)
      endif

c     End program

      end subroutine readTSInput

      end module ts_setup

c module nk_setup
c ######################################################################
      module nk_setup

        use grid, ONLY: np

        integer :: maxitnwt,maxitgm,iguess,maxksp,method,global

        real(8) :: tolgm=0d0,rtol=0d0,atol=0d0,stol=0d0,mf_eps=1d-8

        logical :: fpa=.false.

        !PC setup
        integer :: pc_iter        = 1     !Number of SI iterations in preconditioner
     .            ,mg_ores        = 0     !Restriction order for MG
     .            ,mg_oprol       = 2     !Prolongation order for MG
     .            ,mg_coarse_size = 2     !Minimum coarse grid size (mg_ratio^mg_coarse_size)
     .            ,mg_vcyc        = 1     !Maximum number of MG V-cycles
     .            ,mg_zebra_it    = 1     !Number of relaxations in ZEBRA smoothing
     .            ,mg_zebra_dir   = 0     !Preferred smoothing direction in ZEBRA smoothing (all if zero)
     .            ,sm_ncolors     = 8     !Number of colors for diagonal calc., smoothing
     .            ,sm_iter        = 5     !Number of MG smoothing passes

        real(8) :: pc_tol        = 1d-3   !PC solvers convergence tolerance
     .            ,sm_omega      = 0.7    !Damping for MG smoother (only when 'jb' is chosen)
     .            ,mg_zebra_omg  = 1d0    !Damping for ZEBRA smoothing

        character(2) :: pc_type  = 'id'   !Type of preconditioner (application dependent)
     .                 ,sm_type  = 'jb'   !Type of smoother ('jg','gs','gm')

        logical :: mg_glrkin     =.false. !Whether to use Galerking MG
     .            ,mg_gm_coarse  =.false. !Whether to use GMRES as the coarse solve
     .            ,sm_line_relax =.false. !Whether to perform line relaxation in MG
     .            ,sm_zebra_relax=.false. !Whether to perform ZEBRA line relaxation in MG
     .            ,mg_debug      =.false. !Whether to go into debugging mode in MG
     .            ,mg_cvrg_test  =.false. !Whether to perform convergence test of smoother
     .            ,pc_asm        =.false. !Whether to do additive Schwartz PC
     .            ,pc_debug      =.false. !PC debugging flag
     .            ,pc_id         =.false. !Whether we are doing identity PC
     .            ,pc_gm_driver  =.false. !whether GMRES drives PC solves or not

cc        type :: pcsetup
cc          integer :: sm_iter,maxvcyc,pc_max_iter,debug_it
cc
cc          integer :: pc_iter,ores,oprol,mg_coarse_size,ncolors
cc
cc          real(8) :: mgtol,omeg_smth
cc
cc          character(2) :: smooth
cc          character(10):: pc_type
cc
cc          logical :: gm_smooth
cc     .              ,mg_glrkin
cc     .              ,gm_coarse
cc     .              ,line_relax
cc     .              ,zebra_relax
cc     .              ,inv_B
cc     .              ,v_solve
cc     .              ,mg_debug
cc     .              ,cvrg_test
cc     .              ,asm_PC
cc        end type pcsetup
cc
cc        type(pcsetup) :: pc__setup

      contains

c     nk_setup_consistency
c     ##################################################################
      subroutine nk_setup_consistency

      implicit none

      pc_id = (pc_type == 'id')

cc      if (pc_id) iguess = 0

      if (np > 1) mg_gm_coarse = .false.       !Do not perform GMRES coarse solve in parallel

      if (sm_type /= 'jb') sm_omega = 1d0

      end subroutine nk_setup_consistency

      end module nk_setup

c module counters
c ######################################################################
      module counters

        integer    :: itnewt,itgmres
        integer    :: gmres_tot,newt_tot

      end module counters

c module ts_setup
c ######################################################################
      module ts_setup

        use iosetup

        integer    :: numtime,ndstep,nrstep,inewtime,nrst

        real(8)    :: tmrst,dt,dtbase,cnfactor,alpha,gammat,dt_n
        real(8)    :: tmax,dstep,rstep,time,dfreq,dtexp,damp,dt0,dtold
     .               ,eig_dt

        real(8),allocatable,dimension(:) :: cnf,one_over_dt
     .                                     ,bdfp,bdfn,bdfnm

#if !defined(samrai)
        real(8),allocatable,dimension(:)   :: fold,fsrc,x_old
#endif
        real(8),allocatable,dimension(:,:) :: old_f

        !Initialized variables
        integer    :: itime=0,sm_pass=2,sm_flag=0

        logical    :: timecorr=.false.,restart=.false.,source=.true.
     .               ,predictor=.false.,postprocess=.false.
     .               ,preprocess=.false.
     .               ,chk_src=.false.,source_eval=.false.
     .               ,old_time_eval=.false.,vol_wgt=.true.
     .               ,bdf2=.false.,eigen=.false.
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

c     readTimeStep
c     #################################################################
      subroutine readTimeStep(varray,ierr)

c     -----------------------------------------------------------------
c     Reads time step info from single or multiple files.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray
      integer :: ierr

c     Local variables

      integer :: ifile

c     Begin program

      if (associated(murecord)) then
        do ifile=1,size(murecord)
          ierr = readRecordFile(murecord(ifile),itime,time,dt,varray)
          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) cycle !EOF
        enddo
      else
        ierr = readRecordFile(urecord,itime,time,dt,varray)
      endif

c     End program

      end subroutine readTimeStep

      end module ts_setup

c module timeStepping
c ######################################################################
      module timeStepping

        use parameters

        use variables

        use iosetup

        use ts_setup

#if !defined(petsc_c)
        type(var_array),pointer :: dn,dnm,dnp
        real(8) :: l_np,l_n,et1,et2

      contains

c     findGuess
c     ###############################################################
      subroutine findGuess(vn,vnp)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: vn,vnp

c     Local variables

c     Begin program

      if (itime <= 2) then
        !Previous time step
        vnp = vn
cc      elseif (itime == 2) then
cc        !Linear extrapolation
cc        l_np = dt/dt_n
cc        call AXPYDerivedType(1d0,vn,l_np,dn,vnp)
cc        dnp = dn   !Initialize recursion relation
      else
        !Quadratic extrapolation
        l_np = dt/dt_n
        et1 = l_np*(1d0+l_n*(1+l_np)/(1+l_n))
        et2 =-l_np*l_n**2*(1+l_np)/(1+l_n)
        call AXPYDerivedType(et1,dn,et2,dnm,dnp) !Recursion relation
        call AXPYDerivedType(1d0,vn,1d0,dnp,vnp) !vnp = vn + dnp
      endif

cc      write (*,*) itime,l_np,l_n,et1,et2

cc      if (itime <= 2) then
cc        !Previous time step
cc        vnp = vn
cccc      elseif (itime == 2) then
cccc        !FE predictor
cccc        dn = dt*old_f(:,1)
cccc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp)
cc      else
cc        !AB@ predictor
cc        old_f(:,2) = 0.5*dt*(3*old_f(:,1)-old_f(:,2))
cc        dn = old_f(:,2)
cc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp) !vnp = vn - dt/2*(3*f_n - f_n-1)
cc      endif

c     End program

      end subroutine findGuess

c     storeTSinfo
c     ###############################################################
      subroutine storeTSinfo(vn,vnp)

c     ---------------------------------------------------------------
c     Advance recursion relation
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: vn,vnp

c     Local variables

c     Begin program

      l_n = l_np

      call AXPYDerivedType(1d0,vnp,-1d0,vn,dnp)

      dnm = dn
      dn  = dnp

cc      old_f(:,2) = old_f(:,1)
cc      call evaluateNonlinearFunction(vnp,old_f(:,1))

c     End program

      end subroutine storeTSinfo

c     correctTimeStep
c     ###############################################################
      subroutine correctTimeStep(itm,ierr,dt_exp)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer          :: ierr,itm
      real(8),optional :: dt_exp

c     Local variables

c     Begin program

c     Estimate local growth rate

      if (sm_flag == 1) call calculate_gammat

c     Find new time step

      call findNewTimeStep(itm,ierr,dt_exp=dt_exp)

c     End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

c     Calculation of local growth rate for CN

        integer     :: ieq

        real(8)     :: dmag1,dmag2,dpert,mag(neqd)

        real(8)     :: array(ilom:ihip,jlom:jhip,klom:khip)

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

      end subroutine correctTimeStep

#else

      contains

c     correctTimeStepPetsc
c     ######################################################################
      subroutine correctTimeStepPetsc(dn,dnh,dnp,ierr,dt_to_c)

      use counters

      use parameters

      implicit none

c     Call variables

      integer     :: ierr
      real(8)     :: dn(neqd),dnh(neqd),dnp(neqd),dt_to_c

c     Local variables

c     Begin program

c     Estimate local growth rate

      if (sm_flag == 1) call calculate_gammat

c     Find new time step

      call findNewTimeStep(itime+1,ierr)

c     Update counters (only if timeStep is successful)

      itime  = itime + 1

      time   = time  + dt
      tmrst  = tmrst + dt
      nrst   = nrst  + 1

      dt_to_c = dt

c     End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

c     -----------------------------------------------------------------------
c     Calculation of local growth rate for CN
c    ------------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer     :: ieq

        real(8)     :: mag(neqd)

c     Begin program

        where (dnp /= dn .and. dnh /= 0d0)
          mag = dt*dnh/(dnp-dn)
        elsewhere
          mag = 1e30
        end where

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

      end subroutine correctTimeStepPetsc

#endif

c     findNewTimeStep
c     ####################################################################
      subroutine findNewTimeStep(itm,ierr,dt_exp)

c     --------------------------------------------------------------------
c     Find new time step
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer    ::  ierr,itm
      real(8),optional :: dt_exp

c     Local variables

c     Begin program

      dtold = dt
      dt_n  = dt

      if (PRESENT(dt_exp)) then
        dt = dt_exp
      else
        call adapt_dt

        call calculate_cnfactor
      endif

      alpha = 1d0 - cnfactor

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
          if (my_rank == 0 .and. ierr == 2) write (*,240)
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
              if  (itm == 1) then
                dt = 0.1*dtbase
                dtold = dt
              elseif (itm == 2) then
                dt = dtbase - dt
              elseif (itm == 3) then
                dt = dtbase
              else
                dt = min(dtbase,dt*coef2)
              endif
            endif

          case default   !Default time stepping
            if (    (itm == 1.and.(.not.timecorr))
     .          .or.(itm == 2.and.      timecorr )) then
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

        if (itm.eq.1) then
          cnfactor = .5
        elseif (itm.le.sm_pass+1) then
          cnfactor = .3
        else
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      elseif (sm_flag == 2) then  !BDF2

        cnfactor = 0d0

        bdf2 = (itm > 1 .and. (.not.timecorr) )
     .    .or. (itm > 2 .and.       timecorr  )

        if (bdf2) then
          !All these are affected by factor 1/dt (done in evaluateNLF)
          cnp = (2*dt+dtold)/(dt+dtold)
          cn  = -(dt+dtold)/dtold
          cnm = -cnp-cn    !Sum should be zero identically
        endif

      endif

      end subroutine calculate_cnfactor

      end subroutine findNewTimeStep

      end module timeStepping

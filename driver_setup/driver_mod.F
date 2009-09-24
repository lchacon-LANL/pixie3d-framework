c module icond
c ######################################################################
      module icond

        implicit none

        integer :: nh1,nh2,nh3,npw1=1,npw2=1,npw3=1

        logical :: odd,random

        real(8),dimension(10) :: pert

      end module icond

c module nk_setup
c ######################################################################
      module nk_setup

        integer    ::  maxitnwt,maxitgm,iguess,maxksp,method,global

        real(8)    ::  tolgm=0d0,rtol=0d0,atol=0d0,stol=0d0,mf_eps=1d-8

      end module nk_setup

c module counters
c ######################################################################
      module counters

        integer    :: itnewt,itgmres,itwhistler
        integer    :: gmres_tot,newt_tot,wh_tot

      end module counters

c module ts_setup
c ######################################################################
      module ts_setup

        integer    :: numtime,ndstep,nrstep,sm_pass,sm_flag=0,inewtime

        real(8)    :: dt,dtbase,cnfactor,alpha,gammat
        real(8)    :: tmax,dstep,rstep,time,dfreq,dtexp,damp,dt0,dtold
     .               ,eigstep

#if !defined(samrai)
        real(8),allocatable,dimension(:)   :: fold,fsrc
#endif
        real(8),allocatable,dimension(:,:) :: old_f

        real(8)    :: vx_max,vy_max,vz_max,bx_max,by_max,bz_max

        logical    :: timecorr=.false.,restart=.false.,source=.true.
     .               ,predictor=.false.,postprocess=.false.
     .               ,chk_src=.false.

        real(8)    :: tmrst,cnp,cn,cnm,dt_n

        integer    :: itime,nrst

        real(8),allocatable,dimension(:) :: cnf,one_over_dt
     .                                     ,bdfp,bdfn,bdfnm

        logical    :: vol_wgt=.true.,bdf2=.false.,eigen=.false.

      end module ts_setup

c module timeStepping
c ######################################################################
      module timeStepping

        use parameters

        use variables

        use iosetup

        use generalPurposeFunctions

        use ts_setup

#if !defined(petsc_c)
        type(var_array),pointer :: dn,dnm,dnp
        real(8)    :: l_np,l_n,et1,et2

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

          dpert = integral(nxl,nyl,nzl,array,1,1,1,.true.)

          dpert = sqrt(dpert)

          array = (u_n%array_var(ieq)%array
     .            -u_0%array_var(ieq)%array )
          array = array*array

          dmag1 = integral(nxl,nyl,nzl,array,1,1,1,.true.)

          array = (u_np%array_var(ieq)%array
     .            +u_n %array_var(ieq)%array
     .         -2.*u_0 %array_var(ieq)%array )
          array = array*array

          dmag2 = integral(nxl,nyl,nzl,array,1,1,1,.true.)

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
      endif

      call calculate_cnfactor

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

            if (   (itm <= sm_pass+1.and.(.not.timecorr))
     .         .or.(itm <= sm_pass+2.and.      timecorr )) then
              dt = 0.5*dtbase
            elseif ((itm == sm_pass+2.and.(.not.timecorr))
     .          .or.(itm == sm_pass+3.and.      timecorr )) then
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
            endif

          case(2)  !BDF time stepping
            if     (    (itm == 1 .and.(.not.timecorr))
     .              .or.(itm == 2 .and.      timecorr )) then
              dt = 0.1*dtbase
            elseif (    (itm == 2 .and.(.not.timecorr))
     .              .or.(itm == 3 .and.      timecorr )) then
cc              dt = min(dtbase,2*dt)
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
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

      cnp = 1d0
      cn  =-1d0
      cnm = 0d0

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

      alpha = 1d0 - cnfactor

      end subroutine calculate_cnfactor

      end subroutine findNewTimeStep

      end module timeStepping

c module icond
c ######################################################################
      module icond

        implicit none

        integer :: nh1,nh2,nh3

        logical :: odd,random

        real(8),dimension(10) :: pert

      end module icond

c module newtongm
c ######################################################################
      module newtongm

#if !defined(petsc)
        use nk_mod
#endif

        integer(4) ::  maxitnwt,maxitgm,iguess,maxksp,method,global

        real(8)    ::  tolgm=0d0,rtol=0d0,atol=0d0,stol=0d0,mf_eps=1d-8

      end module newtongm


c module precond_setup
c ######################################################################
      module precond_setup

        integer(4)    :: nsweep,maxvcyc,precpass,debug_it
        real(8)       :: mgtol,omega

        character*(10):: precon

        logical       :: asm_PC=.false.

      end module precond_setup

c module counters
c ######################################################################
      module counters

        integer(4) :: itnewt,itgmres,itwhistler
        integer(4) :: gmres_tot,newt_tot,wh_tot

      end module counters

c module timeStepping
c ######################################################################
      module timeStepping

        use parameters

        use variables

        use iosetup

        use generalPurposeFunctions

        integer(4) :: numtime,ndstep,nrstep,sm_pass,sm_flag=0,inewtime

        real(8)    :: dt,dtbase,cnfactor,alpha,gammat
        real(8)    :: tmax,dstep,rstep,time,dfreq,dtexp,damp,dt0,dtold

        real(8),allocatable,dimension(:)   :: fold,fsrc
        real(8),allocatable,dimension(:,:) :: old_f

        real(8)    :: vx_max,vy_max,vz_max,bx_max,by_max,bz_max

        logical    :: timecorr=.false.,restart=.false.,source=.true.
     .               ,predictor=.false.

        real(8)    :: tmrst,cnp,cn,cnm,dt_n

        integer(4) :: itime,nrst

        real(8),allocatable,dimension(:) :: cnf,one_over_dt
     .                                     ,bdfp,bdfn,bdfnm

        logical    :: vol_wgt=.true.,bdf2=.false.

#if !defined(petsc)
        type(var_array) :: dn,dnm,dnp
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

      type (var_array) :: vn,vnp

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

      type (var_array) :: vn,vnp

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
      subroutine correctTimeStep(varray,itm,ierr)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer(4)       ::  ierr,itm

      type (var_array) :: varray

cc      external findExplicitDt

c     Local variables

c     Begin program

      dtold = dt

c     Estimate local growth rate

      call calculate_gammat

c     Find new time step

      call findNewTimeStep(itm,ierr)

c     End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

c     Calculation of local growth rate for CN

        integer (4) :: ieq

        real(8)     :: dmag1,dmag2,dpert,mag(neqd)

        real(8)     :: array(0:nxdp,0:nydp,0:nzdp)

c     Begin program

        do ieq=1,neqd

          array = (varray%array_var(ieq)%array
     .            -u_0   %array_var(ieq)%array )
          array = array*array

          dpert = integral(nxd,nyd,nzd,array,1,1,1,.true.)

          dpert = sqrt(dpert)

          array = (u_n%array_var(ieq)%array
     .            -u_0%array_var(ieq)%array )
          array = array*array

          dmag1 = integral(nxd,nyd,nzd,array,1,1,1,.true.)

          array = (varray%array_var(ieq)%array
     .            +u_n   %array_var(ieq)%array
     .         -2.*u_0   %array_var(ieq)%array )
          array = array*array

          dmag2 = integral(nxd,nyd,nzd,array,1,1,1,.true.)

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

c     correctTimeStep
c     ######################################################################
      subroutine correctTimeStep(dn,dnh,dnp,ierr,dt_to_c)

      use counters

      use parameters

      implicit none

c     Call variables

      integer(4)  :: ierr
      real(8)     :: dn(neqd),dnh(neqd),dnp(neqd),dt_to_c

c     Local variables

c     Begin program

cc      write (*,*) 'ierr',ierr

      dtold = dt

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

        integer (4) :: ieq

        real(8)     :: mag(neqd)

c     Begin program

        where (dnp /= dn .and. dnh /= 0d0)
          mag = dt*dnh/(dnp-dn)
        elsewhere
          mag = 1e30
        end where

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

      end subroutine correctTimeStep

#endif

c     findNewTimeStep
c     ####################################################################
      subroutine findNewTimeStep(itm,ierr)

c     --------------------------------------------------------------------
c     Find new time step
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) ::  ierr,itm

c     Local variables

c     Begin program

      call calculate_dt

      cnp = 1d0
      cn  =-1d0
      cnm = 0d0

      if (sm_flag == 1) then

        call calculate_cnfactor

      elseif (sm_flag == 2) then
        if (itm == 1 ) then

          cnfactor = 0.5

          bdf2 = .false.

        else
          cnfactor = 0d0

          bdf2 = .true.

          !All these are affected by factor 1/dt (done in evaluateNLF)
          cnp = (2*dt+dtold)/(dt+dtold)
          cn  = -(dt+dtold)/dtold
          cnm = dt**2/(dtold*(dt+dtold))

        endif
      endif

      alpha = 1. - cnfactor

c     End program

      contains

c     calculate_cnfactor
c     #######################################################################
      subroutine calculate_cnfactor

        if (itm.eq.1) then
          cnfactor = .5
        elseif (itm.le.sm_pass+1) then
          cnfactor = .3
        else
cc          write (*,*) 'gamma ',gammat,'dt ',dt
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      end subroutine calculate_cnfactor

c     calculate_dt
c     #######################################################################
      subroutine calculate_dt

        dt_n = dt

#if !defined(petsc)
        if (timecorr .or. cnfactor == 1d0) then
          if (ierr.eq.0 .and. (itm.eq.1 .or. cnfactor == 1d0)) then
            call findExplicitDt
          else
            call adapt_dt(dtbase)
          endif
        else
          dt = dtbase
        endif
#else
        if (timecorr) then
          call adapt_dt(dtbase)
        else
          dt = dtbase
        endif
#endif

      end subroutine calculate_dt

c     adapt_dt
c     #######################################################################
      subroutine adapt_dt(dtbase)

        real(8) ::    dtbase
        real(8) ::    coef1,coef2

        coef1 = 0.8             !Time subcycling coefficient
        coef2 = 1.05            !Time recovery   coefficient

        if (ierr.gt.0) then
          dt = dt*coef1
          if (dt < 1d-3*dtbase) then
            write (*,*) 'Time step too small'
            write (*,*) 'Aborting...'
            stop
          endif
          if (ierr.eq.2) write (*,240)
          write (*,400) dt
        else
          if (sm_flag == 1 .and. itm.le.sm_pass+1) then
            dt = dtbase/2.
          elseif (itm.eq.(sm_pass+2).or. sm_flag /= 1) then
            dt = dtbase
          else
            dt = min(dtbase,dt*coef2)
          endif
        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step: ',1p,e10.4)

      end subroutine adapt_dt

      end subroutine findNewTimeStep

      end module timeStepping

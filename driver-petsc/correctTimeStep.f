c correctTimeStep
c ######################################################################
      subroutine correctTimeStep(dn,dnh,dnp,ierr,dt_to_c)

      use timeStepping

      use counters

      implicit none

c Call variables

      integer(4)  :: ierr
      real(8)     :: dn(neqd),dnh(neqd),dnp(neqd),dt_to_c

c Local variables

c Begin program

cc      write (*,*) 'ierr',ierr

c Estimate local growth rate

      call calculate_gammat

c Find new time step

      call findNewTimeStep(itime+1,ierr)

c Update counters (only if timeStep is successful)

      itime  = itime + 1

      time   = time  + dt
      tmrst  = tmrst + dt
      nrst   = nrst  + 1

      dt_to_c = dt

c End program

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

c findNewTimeStep
c####################################################################
      subroutine findNewTimeStep(itm,ierr)

c--------------------------------------------------------------------
c     Find new time step
c--------------------------------------------------------------------

      use timeStepping

      use iosetup

      use variables

      use parameters

      implicit none

c Call variables

      integer(4) ::  ierr,itm

c Local variables

c Begin program

      call calculate_dt

      if (sm_flag.eq.1) call calculate_cnfactor

      alpha = 1. - cnfactor

c End program

      contains

c     calculate_cnfactor
c     #######################################################################
      subroutine calculate_cnfactor

cc        if (itm.eq.1) then
cc          cnfactor = .5
cc        elseif (itm.le.sm_pass+1) then
        if (itm.le.sm_pass) then
          cnfactor = .3
        else
cc          write (*,*) 'gamma ',gammat,'dt ',dt
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      end subroutine calculate_cnfactor

c     calculate_dt
c     #######################################################################
      subroutine calculate_dt

        if (timecorr) then
          call adapt_dt(dtbase)
        else
          dt = dtbase
        endif

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
cc          if (itm.le.sm_pass+1) then
          if (itm.le.sm_pass) then
            dt = dtbase/2.
cc          elseif (itm.eq.(sm_pass+2)) then
          elseif (itm.eq.(sm_pass+1)) then
            dt = dtbase
          else
            dt = min(dtbase,dt*coef2)
          endif
        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step:',f7.4)

      end subroutine adapt_dt

      end subroutine findNewTimeStep

c processOldSolution
c ######################################################################
      subroutine processOldSolution(array,imin,imax,jmin,jmax,kmin,kmax
     .                             ,gmits,nwits)

      use variables

      use timeStepping

      use iosetup

      use counters

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax,nwits,gmits

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      type(petsc_array) :: petscarray

      integer(4)        :: ierr,i,j,k,ieq

      real(8) :: mag

c Begin program

      ierr = 0

      call allocatePetscType(petscarray)

      petscarray%array(imin:imax,jmin:jmax,kmin:kmax)
     .         = array(imin:imax,jmin:jmax,kmin:kmax)

c Map old solution

      u_n = petscarray

c Evaluate nonlinear function at old time for theta scheme

      call evaluateNonlinearFunction(u_n,fold)

cc      u_n = u_n - u_np
cc      write (*,*) imin,imax,jmin,jmax,kmin,kmax
cc
cc      do ieq=1,neqd
cc        mag = sqrt(sum(u_n%array_var(ieq)%array**2))
cc        write (*,*) mag
cc      enddo
cc
cccc      call writeDerivedType(u_n,6,.true.)
cc      stop

c Time level plots (xdraw)

      if (nrst.eq.ndstep.or.tmrst.ge.dstep) then
        nrst  = 0
        if (itime.gt.0) tmrst = tmrst - dstep
        call writeRecordFile(itime,time,dt,u_n)
cc        write (*,*) itime,time
      endif

c Output per time step

      itgmres = gmits
      itnewt  = nwits

      call output

c Find new time step

      call correctTimeStep(u_n,itime+1,ierr)

c Update counters (only if timeStep is successful)

      itime  = itime + 1

      time   = time  + dt
      tmrst  = tmrst + dt
      nrst   = nrst  + 1

c Deallocate memory

      call deallocatePetscType(petscarray)

      end subroutine processOldSolution

c correctTimeStep
c####################################################################
      subroutine correctTimeStep(varray,itm,ierr)

c--------------------------------------------------------------------
c     Correct time step
c--------------------------------------------------------------------

      use timeStepping

      use iosetup

      use variables

      use parameters

      implicit none

c Call variables

      integer(4)       ::  ierr,itm

      type (var_array) :: varray

c Local variables

c Begin program

      call calculate_gammat

      call calculate_dt

      if (sm_flag.eq.1) call calculate_cnfactor

      alpha = 1. - cnfactor

c End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

        use generalOperators

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

          if (dpert /= 0d0.and.dmag2 /= 0d0) then
            mag(ieq) = .5*dt*sqrt(dmag2)/(dpert-sqrt(dmag1))
          else
            mag(ieq) = 1e30
          endif

        enddo

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

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


        if (timecorr) then
          if (ierr.eq.0 .and. (itm.eq.1 .or. cnfactor .eq. 1d0)) then
            call findExplicitDt
          else
            call adapt_dt(dtbase)
          endif
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
          if (itm.le.sm_pass+1) then
            dt = dtbase/2.
          elseif (itm.eq.(sm_pass+2)) then
            dt = dtbase
          else
            dt = min(dtbase,dt*coef2)
          endif
        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step:',f7.4)

      end subroutine adapt_dt

      end subroutine correctTimeStep

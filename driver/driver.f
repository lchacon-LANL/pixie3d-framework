      program driver_3d

c ******************************************************************
c  This program drives the time stepping of an arbitrary 3D system.
c ******************************************************************

      use variables

      use timeStepping

      use newtongm

      use counters

      use iosetup

      use grid

      use icond

c Common variables

      implicit none

c Local variables

      integer(4) :: prec_tot,ierr

c Begin program

c Initialize calculation

      call initializeCalculation

c Initialize counters

      gmres_tot = 0
      newt_tot  = 0
      wh_tot    = 0
      ierr      = 0

      itime     = inewtime - 1

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0

c Initial output

      call output

c Time loop

      do

        itime = itime + 1

        if (tmax.gt.0d0.and.time.ge.(tmax-1d-5))  exit
        if (numtime.ge.0.and.itime.ge.(numtime+inewtime)) exit

c     Find new time step

        call correctTimeStep(u_np,itime,ierr)

c     Assign old time solution

        u_n = u_np       !Overloaded assignment

c     Time update

        time   = time  + dt

        call timeStep(u_n,u_np,ierr)

c     Check for error in time stepping

        if (ierr.eq.1) then     !Restart time step

          itime = itime - 1
          time  = time - dt

          if (timecorr) then
            cycle
          else
            exit
          endif

        elseif (ierr.eq.-1) then !Initial guess is exact

          write(*,*)
          write(*,*) '   Found steady state solution.'
          write(*,*) '   Aborting...'

          exit

        endif
 
c     Update counters (only if timeStep is successful)

        tmrst  = tmrst + dt
        nrst   = nrst  + 1

c     Output per time step

        call output

c     Time level data dump

        if (nrst.eq.ndstep.or.tmrst.ge.0.99*dstep) then
cc          write (*,*) 'Dumped here'
          nrst  = 0
          if (itime.gt.0) tmrst = tmrst - dstep
          call imposeBoundaryConditions(u_np,1,1,1)
          call writeRecordFile(urecord,itime,time,dt,u_np)
        endif

      enddo       !End of time loop

c Average explicit time step

      if (cnfactor.eq.1d0) then
        dtexp = dtexp/(itime - inewtime)
        write (*,*) 'Average explicit time step: ',dtexp
      endif

c Final statistics

      prec_tot = gmres_tot + iguess*newt_tot

      write(*,300) 
      write(*,310) (itime-1),float(newt_tot )/(itime-inewtime)
     .                      ,float(gmres_tot)/(itime-inewtime)
     .                      ,float(gmres_tot)/newt_tot
     .                      ,float(wh_tot)   /prec_tot

c Formats

 300  format (/,'Final statistics',/,/,
     .          ' itime  Newt/ndt  GMRES/ndt  GMRES/Newt  Whst/GMRES')
 310  format (i5,3x,f7.1,4x,f7.1,4x,f7.1,4x,f7.1)

c End program

      close (urecord)

      end

c timeStep
c######################################################################
      subroutine timeStep(vn,vnp,ierr)

c----------------------------------------------------------------------
c     Performs time advance of solution from vn to vnp
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use timeStepping

      use newtongm

      use counters

      use iosetup

      implicit none

c Call variables

      integer          :: ierr

      type (var_array) :: vnp,vn

c Local variables

      double precision :: x(ntotd)

c Diagnostics

cc      double precision :: dummy(ntotd),dummy2(nxd),norm
cc      double precision, pointer, dimension (:,:):: arr
cc      integer          :: ieq

c Begin program

      ierr = 0

      itwhistler = 0

c Evaluate nonlinear function at old time for theta scheme

      call evaluateNonlinearFunction(vn,fold)

c Update BC for resistive wall BC using old time quantities

cc      if (reswall) then
cc        !Update BC in phi and store in vn
cc        call resistiveWallFlow(nxd,dummy2)
cc
cc        vn%array_var(1)%array(1:nxd,nyd+1) = dummy2
cc
cc        !Update BC in psi and store in vn
cc        call resistiveWallBC  (nxd,dummy2)
cc
cc        vn%array_var(3)%array(1:nxd,nyd+1) = dummy2
cc
cc      endif

c Implicit update (Newton)

      if (cnfactor.lt.1d0) then

c     Map previous time step solution into Newton vector for initial guess

        x = vn  !Overloaded assignment

c     Newton iteration

        call newtonGmres(neqd,ntotd,x,method,damp,global,dt0
     .                  ,tolgm,maxksp,maxitgm,rtol,atol,0.5*maxitnwt
     .                  ,maxitnwt,itgmres,itnewt,iguess,ilevel,ierr)

c     If no error, map Newton solution to vnp

        if (ierr.eq.0.or.ierr.eq.2) vnp = x   !Overloaded assignment

      else

c Explicit update

        itgmres = maxitgm

        call explicit(vn,vnp,itgmres,ilevel)

        itnewt = 1

      endif

c Update counters

      gmres_tot = gmres_tot + itgmres
      newt_tot  = newt_tot  + itnewt
      wh_tot    = wh_tot    + itwhistler

c End program

      end subroutine timeStep

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

        use generalPurposeFunctions

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

        if (timecorr .or. cnfactor == 1d0) then
          if (ierr.eq.0 .and. (itm.eq.1 .or. cnfactor == 1d0)) then
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
 400    format ('    Subcycling time step... New time step: ',1p,e10.4)

      end subroutine adapt_dt

      end subroutine correctTimeStep

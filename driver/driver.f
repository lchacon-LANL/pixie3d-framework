      program threed_mhd_driver

c ******************************************************************
c  This program computes the implicit solution for the 3-D XMHD
c  system in a generalized coordinate system.
c ******************************************************************

      use variables

      use timeStepping

      use newtongm

      use counters

      use iosetup

      use grid

      use graphics

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

      nplot     = 0
      tmplot    = 0d0

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

        call correctTimeStep(itime,ierr)

c     Assign old time solution

        u_n = u_np       !Overloaded assignment

c     Time update

        call timeStep(u_n,u_np,ierr)

c     Check for error in time stepping

        if (ierr.eq.1) then     !Restart time step

          itime = itime - 1

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

        time   = time   + dt
        tmplot = tmplot + dt
        tmrst  = tmrst  + dt

        nplot  = nplot + 1
        nrst   = nrst  + 1

c     Time level plots (xdraw)

        if (plot) then

          call evaluateDiagnostics(u_np,.false.)

          if (nplot.eq.ndstep.or.tmplot.ge.dstep) then
            nplot  = 0
            if (itime.gt.0) tmplot = tmplot - dstep
            call dumpTimeStepPlots
          endif

        endif

c     Output per time step

        call output

c     Periodic dump to restart

        if (nrst.eq.nrstep.or.tmrst.ge.rstep) then
          nrst  = 0
          tmrst = tmrst - rstep
          call writeRestartFile(nxd,nyd,nzd,u_np)
        endif

      enddo       !End of time loop

c Dump to restart

      call writeRestartFile(nxd,nyd,nzd,u_np)

c Average explicit time step

      if (cnfactor.eq.1d0) then
        dtexp = dtexp/(itime - inewtime)
        write (*,*) 'Average explicit time step: ',dtexp
      endif

c Final statistics

      prec_tot = gmres_tot + iguess*newt_tot

      write(*,300) 
      write(*,310) (itime-1),dfloat(newt_tot )/(itime-inewtime)
     .                      ,dfloat(gmres_tot)/(itime-inewtime)
     .                      ,dfloat(gmres_tot)/newt_tot
     .                      ,dfloat(wh_tot)   /prec_tot

c Close graphics files and create draw*.in files

      if (plot) then
        call finalizeGraphics
        call finalizeDiagnostics
      endif

c Formats

 300  format (/,'Final statistics',/,/,
     .          'itime  Newt/ndt  GMRES/ndt  GMRES/Newt  Whst/GMRES')
 310  format (i4,3x,f7.1,4x,f7.1,4x,f7.1,4x,f7.1)

c End program

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

      type (var_array),target :: vnp,vn

c Local variables

      double precision :: x(ntotd)

c Diagnostics

      double precision :: dummy(ntotd),dummy2(nxd),norm
      double precision, pointer, dimension (:,:):: arr

      integer          :: ieq

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
     .                  ,tolgm,maxksp,maxitgm,rtol,atol,maxitnwt
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

      end subroutine

c readRestartFile
c######################################################################
      subroutine readRestartFile(nx,ny,nz,varray)

c----------------------------------------------------------------------
c     Reads restart file
c----------------------------------------------------------------------

      use timeStepping

      use variables

      use graphics

      implicit none

c Call variables

      integer*4   nx,ny,nz

      type (var_array):: varray

c Local variables

      integer(4) :: ieq

c Begin program

      open(2,file='restart.bin',form='unformatted',status='unknown')

      read (2) time
      read (2) itime

      read (2) nx
      read (2) ny
      read (2) nz

      call readDerivedType(varray,2,.false.)

      read (2) vx_max,vy_max,vy_max
      read (2) bx_max,by_max,by_max
      read (2) diagnostics
      read (2) gammat

      close (2)

c End

      end subroutine readRestartFile

c writeRestartFile
c######################################################################
      subroutine writeRestartFile(nx,ny,nz,varray)

c----------------------------------------------------------------------
c     Writes restart file
c----------------------------------------------------------------------

      use timeStepping

      use variables

      use graphics

      implicit none

c Call variables

      integer*4   nx,ny,nz

      type (var_array):: varray

c Local variables

      integer(4) :: ieq

c Begin program

      open(2,file='restart.bin',form='unformatted',status='unknown')

      write (2) time
      write (2) itime

      write (2) nx
      write (2) ny
      write (2) nz

      call writeDerivedType(varray,2,.false.)

      write (2) vx_max,vy_max,vy_max
      write (2) bx_max,by_max,by_max
      write (2) diagnostics
      write (2) gammat

      close (2)

c End

      end subroutine writeRestartFile

c correctTimeStep
c####################################################################
      subroutine correctTimeStep(itm,ierr)

c--------------------------------------------------------------------
c     Correct time step
c--------------------------------------------------------------------

      use timeStepping

      use iosetup

      implicit none

c Call variables

      integer*4     ierr,itm

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

        real(8) ::    coef1,coef2

        if (itm.eq.1 .or. cnfactor .eq. 1d0) then
          call findExplicitDt(dt)
        else
          call adapt_dt(dtbase)
        endif

      end subroutine calculate_dt

c     adapt_dt
c     #######################################################################
      subroutine adapt_dt(dtbase)

        real(8) ::    dtbase
        real(8) ::    coef1,coef2

        coef1 = 0.8             !Time subcycling coefficient
        coef2 = 1.05            !Time recovery   coefficient

        if (timecorr) then
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
        else
          dt = dtbase
        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step:',f7.4)

      end subroutine adapt_dt

      end subroutine correctTimeStep

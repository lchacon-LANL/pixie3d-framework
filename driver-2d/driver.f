      program twod_mhd_driver

c ******************************************************************
c  This program drives the time stepping of a 2D MHD model. 
c ******************************************************************

      use variables

      use timeStepping

      use newtongm

      use resistiveWall

      use counters

      use iosetup

c Common variables

      implicit none

c Local variables

      integer(4) :: i,j,itime,nrst,nplot,ierr,corr,prec_tot

      real(8)    :: tmrst,tmplot,error,tu0

c Begin program

c Define i/o

      call defineFiles

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

      tu0       = 0d0

      dtexp     = 0d0

c Initialize graphics

cc      if (restart) call evaluateDiagnostics

      if (plot) then
        call initializeGraphics
        call dumpTimeStepPlots(nplot,tmplot,itime)
      endif

c Initial output

      call output(itime)

c Setup variable u00 computation

      if (tu0_step.gt.0d0) then
        tmax = ( (u01-u00)/du0 + 1 )*tu0_step
cc        du0 = (u01-u00)/(tmax/tu0_step-1)
cc        u00 = u00 + du0
        write (*,'(a)') 'Variable u00 computation'
        write (*,'(a,1pe10.2)') 'Initial u00 set to ',u00
      endif

c Time loop

      do

        itime = itime + 1

        if (tmax.gt.0d0.and.time.ge.(tmax-1d-5))  exit
        if (numtime.ge.0.and.itime.ge.(numtime+inewtime))   exit

c     Find new time step

        call correctTimeStep(itime,ierr)

c     Store old time solution

        call equateDerivedType(u_np,u_n)

c     Time update

        call timeStep(u_n,u_np,ierr)

c     Find total updates (for debugging)

        do i = 1,neqd
          utmp%array_var(i)%array = u_np%array_var(i)%array
     .                             -u_n %array_var(i)%array
        enddo

c     Check for error in time stepping

        if (ierr.eq.1) then     !Restart time step

          itime = itime - 1

          if (timecorr) then
            cycle
          else
            write(*,*)
            write(*,*) '   Newton failed to converge'
            write(*,*) '   Aborting...'
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
        tu0    = tu0    + dt

        nplot  = nplot + 1
        nrst   = nrst  + 1

c     Output per time step

        call output(itime)

c     Time step diagnostics

        call evaluateDiagnostics

c     Time level plots (xdraw)

        if (plot) call dumpTimeStepPlots(nplot,tmplot,itime)

c     Periodic dump to restart

        if (nrst.eq.nrstep.or.tmrst.ge.rstep) then
          nrst  = 0
          tmrst = tmrst - rstep
          call writeRestartFile(itime,nxd,nyd,u_np)
        endif

c     Periodic change of u_0 in resistive wall

        if (reswall.and.tu0_step.gt.0d0.and.tu0.ge.tu0_step) then
          tu0 = tu0 - tu0_step
          u00 = u00 + du0
          write (*,'(a,1pe10.2)') 'Changing u00 to ',u00
        endif

      enddo       !End of time loop

c Dump to restart

      call writeRestartFile(itime,nxd,nyd,u_np)

c Average explicit time step

      if (cnfactor.eq.1d0) then
        dtexp = dtexp/itime
        write (*,*) 'Average explicit time step: ',dtexp
      endif

c Second order test

c     Write gauge solution

cc      open(2,file='sndord-test.bin',form='unformatted',status='unknown')

cc      write (2) u_np

cc      close (2)

c     Read gauge solution and compare psi

cc      open(2,file='sndord-test.bin',form='unformatted',status='unknown')

cc      read (2) u_n

cc      close (2)

cc      error = sqrt(sum( (u_n %array_var(3)%array
cc     .                  -u_np%array_var(3)%array)**2))

cc      write (*,*) error

c Final statistics

      prec_tot = gmres_tot + iguess*newt_tot

      write(*,300) 
      write(*,310) (itime-1),dfloat(newt_tot )/(itime-inewtime)
     .                      ,dfloat(gmres_tot)/(itime-inewtime)
     .                      ,dfloat(gmres_tot)/newt_tot
     .                      ,dfloat(wh_tot)   /prec_tot

c Close graphics files and create draw*.in files

      if (plot) call finalizeGraphics

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

      use resistiveWall

      implicit none

c Call variables

      integer          :: ierr

      type (var_array),target :: vnp,vn

c Local variables

      double precision :: x(ntotd),damp,dt0,rtol,atol

c Diagnostics

      double precision :: dummy(ntotd),dummy2(nxd)
      double precision, pointer, dimension (:,:):: arr

      integer          :: ieq

c Begin program

      ierr = 0

      itwhistler = 0

c Save old time step solution before imposing res-wall BCs

      call equateDerivedType(vn,utmp)

c Evaluate nonlinear function at old time for theta scheme

      call evaluateNonlinearFunction(vn,fold)

c Update BC for resistive wall BC using old time quantities

      if (reswall) then
        !Update BC in phi and store in vn
        call resistiveWallFlow(nxd,dummy2)

        vn%array_var(1)%array(1:nxd,nyd+1) = dummy2

        !Update BC in psi and store in vn
        call resistiveWallBC  (nxd,dummy2)

        vn%array_var(3)%array(1:nxd,nyd+1) = dummy2

      endif

c Implicit update (Newton)

      if (cnfactor.lt.1d0) then

c     Map previous time step solution into Newton vector for initial guess

        call mapStructureToVector(vn,x)

c     Newton iteration

        damp = 1d0
        dt0 = 1d30
        atol = 0d0
        rtol = tolnewt

        call newtonGmres(neqd,ntotd,x,method,damp,global,dt0
     .                  ,tolgm,kmax,maxitgm,rtol,atol,maxitnwt
     .                  ,maxitnwt,itgmres,itnewt,iguess,ilevel,ierr)

c     If no error, map Newton solution to vnp

        if (ierr.eq.0.or.ierr.eq.2) then
          call mapVectorToStructure(x,vnp,vn)
        else
          call equateDerivedType(utmp,vn)  !Recover unchanged old time step
        endif

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
      subroutine readRestartFile(itime,nx,ny,varray)

c----------------------------------------------------------------------
c     Reads restart file
c----------------------------------------------------------------------

      use timeStepping

      use variables

      use graphics

      implicit none

c Call variables

      integer*4   itime,nx,ny

      type (var_array):: varray

c Local variables

      integer(4) :: ieq

c Begin program

      open(2,file='restart.bin',form='unformatted',status='unknown')

      read (2) nx
      read (2) ny

      call readDerivedType(varray,2,.false.)

      read (2) time
      read (2) itime
      read (2) u_max
      read (2) v_max
      read (2) diagnostics
      read (2) gamma

      close (2)

c End

      return
      end

c writeRestartFile
c######################################################################
      subroutine writeRestartFile(itime,nx,ny,varray)

c----------------------------------------------------------------------
c     Writes restart file
c----------------------------------------------------------------------

      use timeStepping

      use variables

      use graphics

      implicit none

c Call variables

      integer*4   itime,nx,ny

      type (var_array):: varray

c Local variables

      integer(4) :: ieq

c Begin program

      open(2,file='restart.bin',form='unformatted',status='unknown')

      write (2) nx
      write (2) ny

      call writeDerivedType(varray,2,.false.)

      write (2) time
      write (2) itime
      write (2) u_max
      write (2) v_max
      write (2) diagnostics
      write (2) gamma

      close (2)

c End

      return
      end

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

      alpha = 1d0 - cnfactor

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
cc          write (*,*) 'gamma ',gamma,'dt ',dt
          cnfactor = max(.5 - gamma/12.*dt,0d0)
        endif

      end subroutine calculate_cnfactor

c     calculate_dt
c     #######################################################################
      subroutine calculate_dt

        real(8) ::    coef1,coef2

        coef1 = 0.8             !Time subcycling coefficient
        coef2 = 1.05            !Time recovery   coefficient

        if (timecorr.or. cnfactor .eq. 1d0) then
          if (itm.eq.1 .or. cnfactor .eq. 1d0) then
            call findExplicitDt(dt)
cc            if (debug) dt = dtbase
cc            dt = dtbase
          elseif (itm.le.sm_pass+1) then
            dt = dtbase/2.
          elseif (itm.eq.(sm_pass+2)) then
            dt = dtbase
          else
            if (ierr.gt.0) then
              dt = dt*coef1
              if (ierr.eq.2) write (*,240)
              write (*,400) dt
            else
              dt = min(dtbase,dt*coef2)
            endif
          endif
        else
          dt = dtbase
        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step:',f7.4)

      end subroutine calculate_dt

      end subroutine correctTimeStep

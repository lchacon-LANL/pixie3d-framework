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

        u_nm = u_n   !Overloaded assignment
        u_n  = u_np  !Overloaded assignment

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

      if (itime-inewtime > 0d0) then
        write(*,310) (itime-1),float(newt_tot )/(itime-inewtime)
     .                        ,float(gmres_tot)/(itime-inewtime)
     .                        ,float(gmres_tot)/newt_tot
     .                        ,float(wh_tot)   /prec_tot
      else
        write(*,310) 0,0d0,0d0,0d0,0d0
      endif

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

c Begin program

      ierr = 0

      itwhistler = 0

      jit = -1     !This informs nlfunction that we are processing
                   !  n time level info

c Evaluate nonlinear function at old time for theta scheme

      call evaluateNonlinearFunction(vn,fold)

c Implicit update (Newton)

      if (cnfactor.lt.1d0) then

c     Initial guess

        if (predictor) then
          call findGuess(vn,vnp)
          x = vnp               !Overloaded assignment
        else
          if (bdf2) then
            call AXPYDerivedType(-cn/cnp,vn,-cnm/cnp,u_nm,vnp)
            x = vnp
          else
            x = vn
          endif
        endif
        
c     Newton iteration

        nk_conf%etak_meth      = method
        nk_conf%ksmax          = maxksp
        nk_conf%gmmax          = maxitgm
        nk_conf%nwt_max_it_acc = maxitnwt/2
        nk_conf%nwt_max_it_rej = maxitnwt
        nk_conf%global_meth    = global

        nk_conf%eta0   = tolgm
        nk_conf%damp   = damp
        nk_conf%pdt0   = dt0
        if (atol  > 0d0) nk_conf%atol   = atol
        if (rtol  > 0d0) nk_conf%rtol   = rtol
        if (stol  > 0d0) nk_conf%stol   = stol
        if (mf_eps> 0d0) nk_conf%mf_eps = mf_eps

        nk_conf%krylov_method='fg'

        call nk(neqd,ntotd,x,iguess,ilevel,ierr)

        itgmres = nk_conf%gm_it_out
        itnewt  = nk_conf%nwt_it_out

c     Map Newton solution to vnp (if no error)

        if (ierr.eq.0.or.ierr.eq.2) then
          vnp = x               !Overloaded assignment
          if (predictor) call storeTSinfo(vn,vnp)
        elseif (ierr == 1) then
          vnp = vn
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

      end subroutine timeStep

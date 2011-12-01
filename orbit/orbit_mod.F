
c module orbit
c ######################################################################
      module orbit

        use xdraw_io

        use spline_field

        use grid, ONLY: my_rank,pstop,messg

        !Diagnostics
        integer,private :: nunitd(4)
        logical,private :: flux = .false.

        !B field diagnostics private variables
        logical,private :: bdiag=.false.
        integer,private :: inx,nnx,nny
        integer,private :: icross,ncross=100000,nbins=100
        real(8),private,allocatable,dimension(:,:) :: fpos,qdif

        !Integrator public variables
        real(8) :: xold,yold,zold,lsode_out(10)
        logical :: sing_point=.false.,sln,exit_loop,orb_init

        integer :: iiout

        !LSODE private variables
        INTEGER,private :: LIW, LRW, NEQ
        PARAMETER  (NEQ=20,LIW = 20,LRW = LIW + 16*NEQ)

        INTEGER,private :: IOPT,IOUT,ISTATE=0,ITASK,ITOL
     .                    ,IWORK(LIW),MF,NEQD(1)
        REAL(8),private :: LS_ATOL(1),LS_RTOL(1),RWORK(LRW)
     .                    ,T,TOUT,YY(NEQ)

        !ORBIT private variables
        type,private :: ode_solver
          integer :: niter     !Number of nonlinear iterations
          integer :: order     !Order of integration (1,2)
          real(8) :: atol      !Absolute convergence tolerance
          real(8) :: rtol      !Relative convergence tolerance
          real(8) :: dtime     !Time step
          real(8) :: tmax      !Maximum time
          real(8) :: smax      !Maximum arc length
          logical :: solen     !Whether use solenoidal version or LSODE
          logical :: picard    !Whether to use picard or Newton
          logical :: ds_is_Bds !Whether we use B.ds instead of ds as independent variable
          logical :: lsode_dt_lim  !Whether to limit lsode time step
          logical :: lsode_s_indep !Whether to arc-length is indep. var. in LSODE
          logical :: lsode_to_smax !Whether to integrate LSODE exactly to target arc-length
          integer :: lsode_neq !Number of equations for LSODE
        end type ode_solver

        type,private :: splines
          integer :: order  !Order of splines
          integer :: a_gauge!Gauge for vector potential (1,2)
        end type splines

        type,private :: orbit_def
          type(ode_solver) :: ode_conf
          type(splines)    :: spline_conf
        end type orbit_def

        type(orbit_def),private  :: orb_default
        type(ode_solver),private :: ode_conf

        integer,private :: itime,iterror,ssorder

        integer,private :: nunitx,nunity,nunitz,nunitrz
        integer,private :: npointx,npointy,npointz

        real(8),private :: tvold,told,time,tmax
     .                    ,svold,sold,ss  ,smax
        real(8),private :: x,y,z,xcar,ycar,zcar,xss0,yss0,zss0
     .                    ,xvold,yvold,zvold,x00,y00,z00
        real(8),private :: xss,yss,zss,tss,rmin

        real(8),private :: xmin,xmax,ymin,ymax,zmin,zmax

        character(50),private :: ssx,ssy,ssz,ssrz

        logical,private :: cartesian=.false.,cross=.false.
     .                    ,cartesian_flag=.false.
        logical,private :: lssx,lssy,lssz,rzss,ss_io=.false.

      contains

c     B_diag_setup
c     #################################################################
      subroutine B_diag_setup(b_diag,ninitx,ninity,ninitz)

c     ----------------------------------------------------------------------
c     Set up B diagnostics: set size of non-periodic direction (nnx)
c     and remaining periodic direction (nny).
c     ----------------------------------------------------------------------

      implicit none

      logical :: b_diag
      integer :: ninitx,ninity,ninitz

c     Begin program

      bdiag = b_diag
      if (.not.bdiag) return

      if ( (lssx.and.lssy).or.(lssx.and.lssz).or.(lssy.and.lssz)) then
        write (*,*)
     .       'Cannot do B diagnostics for more than one Surf-section'
        stop
      endif

c     Initialize array sizes and update orbit crossings counter

      if (lssx) then
        nnx = ninity
        nny = ninitz
      endif

      if (lssy) then
        nnx = ninitx
        nny = ninitz
      endif

      if (lssz) then
        nnx = ninitx
        nny = ninity
      endif

c     Setup diag B diagnostic arrays

      allocate(fpos(-nbins/2:nbins/2,ncross))
      fpos = 0d0

      allocate(qdif(nnx,ncross))
      qdif = 0d0

c     End program

      end subroutine B_diag_setup

c     B_diag_init
c     #################################################################
      subroutine B_diag_init(initx,inity,initz,x0,y0,z0)

c     ----------------------------------------------------------------------
c     Initialize counters for B diagnostic at the beginning of each orbit.
c     Counter is a binning index in non-periodic direction (i.e., radius in
c     polar coordinates).
c     ----------------------------------------------------------------------

      implicit none

      integer :: initx,inity,initz
      real(8) :: x0,y0,z0

c     Begin program

      if (.not.bdiag) return

      !Initialize orbit-crossings counter
      icross = 0

      if (lssx) then
        inx = inity             !Only makes sense if Y is not periodic...

        if (x0 == xss0) then
          icross = icross + 1
          fpos(0,icross) = fpos(0,icross) + 1
        endif
      endif

      if (lssy) then
        inx = initx

        if (y0 == yss0) then
          icross = icross + 1
          fpos(0,icross) = fpos(0,icross) + 1
        endif
      endif

      if (lssz) then
        inx = initx

        if (z0 == zss0) then
          icross = icross + 1
          fpos(0,icross) = fpos(0,icross) + 1
        endif
      endif

c     End program

      end subroutine B_diag_init

c     B_diag_update
c     #################################################################
      subroutine B_diag_update(xs,x0,ts)

c     ----------------------------------------------------------------------
c     Accumulates B diagnostics counters.
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: xs,x0,ts

c     Local variables

      integer :: bin

c     Begin program

      if (.not.bdiag) return

      icross = icross + 1

      if (icross > ncross) then
        write (*,*) 'Exceeded max number of Poincare crossings'
      else
        write (*,*)
        write (*,*) 'B diag: Poincare crossing ',icross

        !PDF diagnostic
        bin = sign(1d0,xs-x0)
     .       *floor(0.5*abs(x0**2-xs**2)*nbins+0.5)
        fpos(bin,icross) = fpos(bin,icross)+1

        !Add quasilinear diffusion diagnostic
        qdif(inx,icross) = qdif(inx,icross)
     .                   + (0.5*(xs**2-x0**2))**2/ts
      endif

c     End program

      end subroutine B_diag_update

c     B_diag_finalize
c     #################################################################
      subroutine B_diag_finalize

c     -----------------------------------------------------------------
c     Finalizes B diagnostics, including calculation of statistics
c     and dumping of data to file "b_stat.txt".
c     -----------------------------------------------------------------

        implicit none

c     Local variables

        integer :: slice
        real(8) :: rr

        real(8),allocatable,dimension(:) ::  mean,sigma,vol,diff

c     Begin program

        if (.not.bdiag) return

        open(unit=50,file='b_stat.txt',status='unknown')

        allocate(mean(icross),sigma(icross),vol(icross),diff(icross))

        mean  = 0d0
        sigma = 0d0
        vol   = 0d0

        do slice=-nbins/2,nbins/2
          rr = (1.*slice)/(nbins+1)
          vol   = vol   +       fpos(slice,1:icross)
          mean  = mean  + rr   *fpos(slice,1:icross)
          sigma = sigma + rr**2*fpos(slice,1:icross)
        enddo

        mean  = mean/vol
        sigma = sqrt(sigma/vol - mean**2)

        do slice = 1,icross
          diff(slice) = sigma(slice)**2/(zss0+(slice-1)*(zmax-zmin))
        enddo

        do slice=-nbins/2,nbins/2
          write(50,*) (1.*slice)/(nbins+1),fpos(slice,1:icross)/vol
        enddo

        write (*,*)
        write (*,*) 'Mean =',mean
        write (*,*)
        write (*,*) 'Sigma=',sigma
        write (*,*)
        write (*,*) 'Diff PDF =',diff

        do slice = 1,icross
          diff(slice) = sum(qdif(:,slice))/size(qdif,1)/nny
        enddo

        write (*,*)
        write (*,*) 'Diff MC =',diff

        close(50)

        open(unit=50,file='b_dif.txt',status='unknown')
        do slice=1,size(qdif,1)
          write(50,*) slice,qdif(slice,1:icross)/nny
        enddo
        close(50)

        deallocate(mean,sigma,vol,diff,fpos,qdif)

c     End program

      end subroutine B_diag_finalize

c     orbit_conf_def
c     ######################################################################
      subroutine orbit_conf_def

c     ----------------------------------------------------------------------
c     Default configuration for orbit integrator
c     ----------------------------------------------------------------------

      !ODE solver
      orb_default%ode_conf%tmax   = 0d0    !Maximum orbit pseudotime
      orb_default%ode_conf%smax   = 0d0    !Maximum orbit arc length
      orb_default%ode_conf%dtime  = 0.01   !Maximum time step

      orb_default%ode_conf%solen  = .true. !Whether to use solenoidal ODE int.
      orb_default%ode_conf%order  = 1      !Order of integration in solenoidal int. (1,2)
      orb_default%ode_conf%picard = .true. !Whether to use functional iteration
      orb_default%ode_conf%niter  = 20     !Number of nonlinear iterations
      orb_default%ode_conf%rtol   = 0d-8   !Relative tolerance
      orb_default%ode_conf%atol   = 1d-8   !Absolute tolerance

      orb_default%ode_conf%ds_is_Bds = .false.  !Whether we use B.ds instead of ds as independent variable

      !LSODE default setup
      orb_default%ode_conf%lsode_dt_lim = .false. !Do not limit LSODE time step
      orb_default%ode_conf%lsode_s_indep  = .false. !Find arc-length as dependent variable
      if (orb_default%ode_conf%lsode_s_indep) then
        orb_default%ode_conf%lsode_neq = 3 !Number of LSODE equations
      else
        orb_default%ode_conf%lsode_neq = 5 !Number of LSODE equations
      endif
      orb_default%ode_conf%lsode_to_smax  = .false. !Return after every single LSODE step

      !Splines
      orb_default%spline_conf%order = 2    !Spline interpolation order
      orb_default%spline_conf%a_gauge = 2  !Gauge for vector potential

      end subroutine orbit_conf_def

c     orbit_setup
c     ######################################################################
      subroutine orbit_setup(nx,ny,nz,xx,yy,zz
     .                      ,bx,by,bz,bxcar,bycar,bzcar,jac,xmap
     .                      ,tmax,smax,dtime
     .                      ,niter,atol,rtol,order,solen,picard
     .                      ,spline_order,gauge,B_input_is_A,r_min
     .                      ,lsode_neq,lsode_dt_lim,lsode_to_smax
     .                      ,lsode_s_indep,ds_is_Bds)

c     ----------------------------------------------------------------------
c     Sets up orbit integrator.
c     ----------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: nx,ny,nz
        real(8) :: xx(nx),yy(ny),zz(nz)
        real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
        real(8) :: bxcar(nx,ny,nz),bycar(nx,ny,nz),bzcar(nx,ny,nz)
        real(8) :: xmap(nx,ny,nz,3),jac(nx,ny,nz)

        !Optional parameters
        real(8),optional :: tmax,smax,dtime,atol,rtol,r_min
        integer,optional :: niter,order,spline_order,gauge,lsode_neq
        logical,optional :: solen,picard,B_input_is_A,lsode_dt_lim
     .                     ,lsode_to_smax,lsode_s_indep,ds_is_Bds

c     Begin program

        call orbit_conf_def   !Defaults are defined here.

        if (PRESENT(tmax)) then
          orb_default%ode_conf%tmax = tmax
        endif

        if (PRESENT(smax)) then
          orb_default%ode_conf%smax = smax
        endif

        if (PRESENT(dtime)) then
          orb_default%ode_conf%dtime = dtime
        endif

        if (PRESENT(niter)) then
          orb_default%ode_conf%niter = niter
        endif

        if (PRESENT(atol)) then
          orb_default%ode_conf%atol = atol
        endif

        if (PRESENT(rtol)) then
          orb_default%ode_conf%rtol = rtol
        endif

        if (PRESENT(order)) then
          orb_default%ode_conf%order = order
        endif

        if (PRESENT(solen)) then
          orb_default%ode_conf%solen = solen
        endif

        if (PRESENT(picard)) then
          orb_default%ode_conf%picard = picard
        endif

        if (PRESENT(lsode_neq)) then
          orb_default%ode_conf%lsode_neq = lsode_neq
        endif

        if (PRESENT(lsode_dt_lim)) then
          orb_default%ode_conf%lsode_dt_lim = lsode_dt_lim
        endif

        if (PRESENT(lsode_s_indep)) then
          orb_default%ode_conf%lsode_s_indep = lsode_s_indep
        endif

        if (PRESENT(lsode_to_smax)) then
          orb_default%ode_conf%lsode_to_smax = lsode_to_smax
        endif

        if (PRESENT(ds_is_Bds)) then
          orb_default%ode_conf%ds_is_Bds = ds_is_Bds
        endif

        if (PRESENT(spline_order)) then
          orb_default%spline_conf%order = spline_order
        endif

        if (PRESENT(gauge)) then
          orb_default%spline_conf%a_gauge = gauge
        endif

        if (PRESENT(r_min)) then
          rmin = r_min
        else
          rmin = 1d-3
        endif

        ode_conf = orb_default%ode_conf

c     Consistency checks

        if (ode_conf%lsode_to_smax) then
          ode_conf%solen = .false.
          ode_conf%lsode_s_indep = .true.
        endif

cc        if (ode_conf%lsode_s_indep) then
cc          ode_conf%lsode_neq = 3 !Number of LSODE equations for field line integration
cc        else
cc          ode_conf%lsode_neq = 5 !Number of LSODE equations for field line integration
cc        endif

c     Do splines

        call setupSplines(nx,ny,nz,xx,yy,zz
     .                   ,orb_default%spline_conf%order
     .                   ,xmin,xmax,ymin,ymax,zmin,zmax)

        if (ode_conf%solen) then
          call splineA(bx,by,bz,orb_default%spline_conf%a_gauge
     .                ,input_is_A=B_input_is_A)
        else
          call splineB(bx,by,bz)
        endif

        call splineBcar(bxcar,bycar,bzcar)

        call splineJ(jac )
        call splineX(xmap)

      end subroutine orbit_setup

c     orbit_destroy
c     ######################################################################
      subroutine orbit_destroy

c     ----------------------------------------------------------------------
c     Frees up memory from orbit integration
c     ----------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        call killSplines

c     Finalize B diagnostics

        call B_diag_finalize

      end subroutine orbit_destroy

c     orbit_find
c     ######################################################################
      subroutine orbit_find(iout,x0,y0,z0,fwd,line_int,lsode_f,ss_order
     .                     ,t_max,s_max)

c     ----------------------------------------------------------------------
c     Integrates orbit up to specified limit in arc length or pseudotime
c     (specified in module). On input:
c       * iout: level of diagnostic output
c       * x0,y0,z0: starting point of orbit
c       * fwd: whether to go forward or backward
c     Optional input:
c       * line_int: external routine to perform line integrals
c       * ode_line_int: external routine to perform line integrals using LSODE
c       * ss_check: check for surface of section (default: not when
c                   performing integrals along line)
c       * ss_order: order of interpolation for SS (default: 1)
c     ----------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: iout
        real(8),intent(IN) :: x0,y0,z0
        logical :: fwd

        integer,optional :: ss_order
        real(8),optional :: t_max,s_max

        optional :: line_int

        !Interface for generic line integral routine (optional)
        INTERFACE
          function line_int(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: line_int
          logical :: fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function line_int
        END INTERFACE

        optional :: lsode_f

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Local variables

        real(8) :: temp_s,xh,yh,zh,sh
        logical :: chk_ss,integral

c     Begin program

        iiout = iout

        integral = PRESENT(line_int)

        chk_ss = (.not.integral).and.ss_io  !Do SS only if not integrating and if IO has been set up

        !Override maximum limits
        if (PRESENT(t_max)) then
          tmax = t_max
        else
          tmax = ode_conf%tmax 
        endif

        if (PRESENT(s_max)) then
          smax = s_max
        else
          smax = ode_conf%smax
        endif

c     Check ranges

        if((x0.le.xmin.or.x0.ge.xmax)
     .     .and.(bcond(1)/=PER))then
          write(*,*) 'x0(init)=',x0,' out of range:',xmin,xmax
          return
        endif

        if((y0.le.ymin.or.y0.ge.ymax)
     .     .and.(bcond(3)/=PER))then
          write(*,*) 'y0(init)=',y0,' out of range:',ymin,ymax
          return
        endif

        if((z0.le.zmin.or.z0.ge.zmax)
     .     .and.(bcond(5)/=PER))then
          write(*,*) 'z0(init)=',z0,' out of range:',zmin,zmax
          return
        endif

c     Initialize orbit

        call orbit_init(x0,y0,z0)

c     Daniele,23/03/2009: dump the initial point of the orbit

        if (chk_ss) then

          if(lssx .and. x0 == xss0) then
            write(nunitx) real(y0),real(z0)
          endif

          if(lssy .and. y0 == yss0) then
            if (flux) then
              write(nunity) real(0.5*x0**2),real(z0)
            else
              write(nunity) real(x0),real(z0)
            endif
          endif

          if(lssz .and. z0 == zss0) then
            if (flux) then
              write(nunitz) real(0.5*x0**2),real(y0) !In flux
            else
              write(nunitz) real(x0),real(y0)        !In radius
            endif

          endif

        endif

c     Orbit integral

        do

          itime = itime + 1

          !Save old times and positions (for SS interpolation below)
          tvold = told
          told  = time

          svold = sold
          sold  = ss

          xvold = xold
          yvold = yold
          zvold = zold

          xold  = x
          yold  = y
          zold  = z

          !Find new LOGICAL positions, new time, new arc length
          call orbit_step(x,y,z,iterror,forward=fwd,lsode_f=lsode_f)

          !Check for orbit end and interpolate quantities
          if (smax /= 0d0) then
            exit_loop = (abs(ss) >= smax)

            if (exit_loop) then
              x = xold + (smax-abs(sold))/abs(ss-sold)*(x-xold)
              y = yold + (smax-abs(sold))/abs(ss-sold)*(y-yold)
              z = zold + (smax-abs(sold))/abs(ss-sold)*(z-zold)
              ss = sign(1d0,ss)*smax
            endif
          elseif (tmax /= 0d0) then
            exit_loop = (abs(time) >= tmax)

            if (exit_loop) then
              x = xold + (tmax-abs(told))/abs(time-told)*(x-xold)
              y = yold + (tmax-abs(told))/abs(time-told)*(y-yold)
              z = zold + (tmax-abs(told))/abs(time-told)*(z-zold)
              time = sign(1d0,time)*tmax
            endif
          else
            exit_loop = .false.
          endif

          !Line-orbit integration
          if (integral) then
            exit_loop = (line_int(fwd,xold,x,yold,y,zold,z,sold,ss))
     .                  .or.(exit_loop)
          endif

          !Check for intersection with surface of section
          if (iterror /= 2 .and. itime > 2 .and. chk_ss) then
            if (smax /= 0d0) then
              call sschk(ssorder,svold,sold,ss
     .                  ,xvold,xold,x,yvold,yold,y
     .                  ,zvold,zold,z,xss,yss,zss,tss)
            else
              call sschk(ssorder,tvold,told,time
     .                  ,xvold,xold,x,yvold,yold,y
     .                  ,zvold,zold,z,xss,yss,zss,tss)
            endif
          endif

          !Show progress
          if (smax /= 0d0) then
            call orbit_progress(ss  ,smax,fwd,exit_loop)
          elseif (tmax /= 0d0) then
            call orbit_progress(time,tmax,fwd,exit_loop)
cc          else
cc            exit_loop = .false.
          endif

          if (exit_loop .or. orbit_error(iterror)) then
cc          if (exit_loop) then
            if (iiout >=0) write (*,*)
            exit
          endif

        enddo

c       Write orbit separator in XDRAW file

        if (chk_ss) then
          if (lssx) write(nunitx)
          if (lssy) write(nunity)
          if (lssz) write(nunitz)
          if (rzss) write(nunitrz)
        endif

      end subroutine orbit_find

c     orbit_init
c     ######################################################################
      subroutine orbit_init(x0,y0,z0)

c     ----------------------------------------------------------------------
c     Initializes orbit integrator
c     ----------------------------------------------------------------------

        implicit none

        real(8) :: x0,y0,z0

c     Begin program

        x00  = x0
        y00  = y0
        z00  = z0

        x    = x0
        y    = y0
        z    = z0

        xold = x0
        yold = y0
        zold = z0

        time = 0d0  !Initialize pseudo-time along line
        told = 0d0

        ss   = 0d0  !Initialize arc-length along line
        sold = 0d0

        lsode_out = 0d0   !Initialize additional LSODE output

        npointx=0
        npointy=0
        npointz=0

        itime = 0

        exit_loop = .false.

        orb_init = .true.

        if (iiout == 0) then
          write (*,FMT="(a,3f9.4,a)",ADVANCE="NO")
     .         ' Orbit (x,y,z)=',x0,y0,z0,', Progress= '
        elseif (iiout > 0) then
          write (*,*) '**********************************************'
          write (*,FMT="(a,3f9.4)") ' Orbit (x,y,z)=',x0,y0,z0
          write (*,*) '**********************************************'
          write (*,*)
        endif

      end subroutine orbit_init

c     orbit_step
c     ##################################################################
      subroutine orbit_step(x,y,z,ierror,forward,lsode_f)

c     ------------------------------------------------------------------
c     Crank-Nicolson (Picard/Newton) leapfrog, symmetrized (or not).
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierror
      real(8) :: x,y,z

      logical,optional :: forward

c     Local Variables

      real(8) :: x0,y0,z0,x1,y1,z1,half,jach,dt=10.
      real(8) :: vx,vy,vz,error(ode_conf%niter),atol,rtol,dtime
      real(8) :: xh,yh,zh,xnew,ynew,znew,xc,yc,zc,ds,cc
      integer :: iter,order,niter
      logical :: restart,nt,fwd,lsode_3eq,lsode_to_smax

      save dt

      optional :: lsode_f

      !Interface for optional LSODE function evaluation
      INTERFACE
        subroutine lsode_f(neq,t,y,ydot)
        integer :: neq
        real(8) :: t,y(neq),ydot(neq)
        end subroutine lsode_f
      END INTERFACE

c     Begin program

c     ODE configuration

      dtime = ode_conf%dtime

      niter = ode_conf%niter
      rtol  = ode_conf%rtol
      atol  = ode_conf%atol
      order = ode_conf%order
      sln   = ode_conf%solen
      nt    = .not.ode_conf%picard
      lsode_to_smax = ode_conf%lsode_to_smax

      if (lsode_to_smax) then
        lsode_3eq=.true.
      else
        lsode_3eq = ode_conf%lsode_s_indep
      endif

      if (PRESENT(forward)) then
        fwd = forward
      else
        fwd = .true.
      endif

      ierror = 0    !Initialize domain bounds error flag

      restart = .false. !Initialize restarting status flag

      if (fwd) then
        cc = 1d0
      else
        cc =-1d0
      endif

c     Check for SP

      if (.not.cartesian) cartesian = (bcond(1)==SP).and.(x<rmin)

c     Time integration

      if (sln .and. (.not.cartesian)) then  !Solenoidal

        call vp_solver(atol,rtol,x,y,z,ierror)

      else                                  !LSODE (non-solenoidal)

        if (lsode_3eq) then
          call lsode_step_3eq(atol,rtol,x,y,z,ierror)
        else
          call lsode_step    (atol,rtol,x,y,z,ierror)
        endif
      endif

c     Output

      if (iiout > 0) then
        write (*,'(a,e10.3,a,e10.3,a,e10.3,a,e10.3,a,3es12.3,a,l3)')
     .         ' Time=',time
     .        ,' ;  dt=',dt
     .        ,' ;  Arc length=',ss
     .        ,' ;  ds=',ds
     .        ,' ;  Position=',x,y,z
     .        ,' ;  Forward=',fwd
      endif

c     End program

      contains

c     lsode_step
c     ####################################################################
      subroutine lsode_step(aatol,rrtol,x,y,z,iterr)

c     --------------------------------------------------------------------
c     Performs orbit integration step using LSODE (arbitrary order < 12th,
c     implicit Adams-Moulton, adaptive time stepping).
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: iterr
      real(8) :: x,y,z,aatol,rrtol

c     Local variables

      real(8) :: xc,yc,zc
      real(8) :: xh,yh,zh,jach

      logical  :: external_f

c     Interfaces

      INTERFACE
        subroutine ff(neq,t,y,ydot)
        integer :: neq
        real(8) :: t,y(neq),ydot(neq)
        end subroutine ff
      END INTERFACE

c     Begin program

      if (cartesian .and. (.not.sing_point)) ISTATE=1 !Restart LSODE for cartesian step
      external_f = PRESENT(lsode_f)                   !LSODE function is external

      sing_point = cartesian

      NEQD = ode_conf%lsode_neq    !Number of equations

c     Initial conditions

      if (cartesian) then
        !Find Cartesian coordinates
        call evalX(x,y,z,xc,yc,zc)

        YY(1) = xc
        YY(2) = yc
        YY(3) = zc
        YY(4) = time
        YY(5) = ss
        if (NEQD(1) > 5) YY(6:NEQD(1))=lsode_out(1:NEQD(1)-5)

        TOUT = 0.1*cc*dtime !Output time value (ignored after first call, since ITASK=2)
      else
        YY(1) = x
        YY(2) = y
        YY(3) = z
        YY(4) = time
        YY(5) = ss
        if (NEQD(1) > 5) YY(6:NEQD(1))=lsode_out(1:NEQD(1)-5)

        TOUT  = cc*dtime    !Output time value (ignored after first call, since ITASK=2)
      endif

c     LSODE configuration and call

      T    = 0d0   !Independent variable

      ITOL = 1      !Whether atol is array or not (2 => true)

      if (ISTATE /= 3) then
        LS_RTOL = rrtol            !Relative tolerance
        LS_ATOL = aatol            !Absolute tolerance
      endif

      ITASK = 2     !Mode of operation: take one step and return

      !Whether it is first call (ISTATE=1) or subsequent call (ISTATE=2)
      if ((time == 0d0).or. (sln .and. (.not.cartesian))) ISTATE = 1

      if (ode_conf%lsode_dt_lim) then
cc        write (*,*) 'LSODE -- limiting time step'
        IOPT = 1                !Whether optional inputs are needed
        RWORK(6) = dtime        !Maximum time step allowed (optional input)
      else
        IOPT = 0                !Whether optional inputs are needed
      endif

      MF = 10       !Specifies non-stiff problem: implicit Adams-Moulton method

cc      write (*,*)'Number of equations',NEQD
cc      write (*,*)'Input state=',YY(1:NEQD(1))
cc      write (*,*)'Current time=',T
cc      write (*,*)'Next time=',TOUT
cc      write (*,*)'Tolerances=',ITOL,LS_RTOL,LS_ATOL
cc      write (*,*)'LSODE task=',ITASK
cc      write (*,*)'LSODE state=',ISTATE
cc      write (*,*)'Optional inputs required?',IOPT
cc      write (*,*)'Real work array size=',size(RWORK),' declared=',LRW
cc      write (*,*)'Integer work array size=',size(IWORK),' declared=',LIW
cc      write (*,*)'LSODE method=',MF

      if (external_f) then
        CALL LSODE (lsode_f,NEQD(1),YY(1:NEQD(1)),T,TOUT
     .             ,ITOL,LS_RTOL,LS_ATOL
     .             ,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,lsode_f, MF)
      else
        CALL LSODE (ff     ,NEQD(1),YY(1:NEQD(1)),T,TOUT
     .             ,ITOL,LS_RTOL,LS_ATOL
     .             ,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,ff      ,MF)
      endif

      select case(ISTATE)
      case(2)       !Uneventful
        continue
      case(-2)      !Rescale tolerances
        ISTATE = 3
        RTOL=RTOL*RWORK(14)
        ATOL=ATOL*RWORK(14)
      case default  !Exit LSODE
        write (*,*) 'Error code ISTATE=',ISTATE
        write (*,*) 'See lsode.f for instructions'
        write (*,*) 'Aborting...'
        stop
      end select

c     Postprocessing to find logical coordinates

      if (cartesian) then
        x = xold ; y = yold ; z = zold
        call evalXi(iiout,YY(1),YY(2),YY(3),x,y,z,iterr)

        if (iterr /= 0) then
          iterr = 6
          return
        endif

        if (x >= rmin) then
          iterr = -2
          ISTATE = 1   !Instruct LSODE to start anew for next time step
          cartesian = .false.
          sing_point= .false.
        endif
      else
        x = YY(1)
        y = YY(2)
        z = YY(3)
      endif

c     Advance time and arc length

      dt   = cc*YY(4)-time !For output
      ds   = cc*YY(5)-ss   !For output

      time = cc*YY(4)
      ss   = cc*YY(5)

      if (NEQD(1) > 5) lsode_out(1:NEQD(1)-5)=cc*YY(6:NEQD(1))

c     End program

      end subroutine lsode_step

c     lsode_step_3eq
c     ####################################################################
      subroutine lsode_step_3eq(aatol,rrtol,x,y,z,iterr)

c     --------------------------------------------------------------------
c     Performs orbit integration step using LSODE (arbitrary order < 12th,
c     implicit Adams-Moulton, adaptive time stepping).
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: iterr
      real(8) :: x,y,z,aatol,rrtol

c     Local variables

      real(8) :: xc,yc,zc
      real(8) :: xh,yh,zh,jach

      logical  :: external_f

c     Interfaces

      INTERFACE
        subroutine ff_3eq(neq,t,y,ydot)
        integer :: neq
        real(8) :: t,y(neq),ydot(neq)
        end subroutine ff_3eq
      END INTERFACE

c     Begin program

      if (cartesian .and. (.not.sing_point)) ISTATE=1 !Restart LSODE for cartesian step
      external_f = PRESENT(lsode_f)                   !LSODE function is external

      sing_point = cartesian

      NEQD = ode_conf%lsode_neq    !Number of equations

c     Initial conditions

      if (cartesian) then
        !Find Cartesian coordinates
        call evalX(x,y,z,xc,yc,zc)

        YY(1) = xc
        YY(2) = yc
        YY(3) = zc
        if (NEQD(1) > 3) YY(4:NEQD(1))=lsode_out(1:NEQD(1)-3)

        TOUT = 0.1*cc*dtime !Output time value (ignored after first call, since ITASK=2)
      else
        YY(1) = x
        YY(2) = y
        YY(3) = z
        if (NEQD(1) > 3) YY(4:NEQD(1))=lsode_out(1:NEQD(1)-3)

        TOUT  = cc*dtime    !Output time value (ignored after first call, since ITASK=2)
      endif

c     LSODE configuration and call

      T = time   !Independent variable

      ITOL = 1      !Whether atol is array or not (2 => true)

      if (ISTATE /= 3) then
        LS_RTOL = rrtol            !Relative tolerance
        LS_ATOL = aatol            !Absolute tolerance
      endif

      if (.not.lsode_to_smax) then
        ITASK = 2   !Mode of operation: take one step and return
      else
        TOUT = smax !Set smax for end of integration
cc        ITASK = 1   !Mode of operation: integrate to end and return at TOUT=smax by interpolation
        ITASK = 4   !Mode of operation: integrate to end and return at TOUT=smax w/o overshooting
      endif

      !Whether it is first call (ISTATE=1) or subsequent call (ISTATE=2)
      if ((time == 0d0).or. (sln .and. (.not.cartesian))) ISTATE = 1

      IOPT = 0                  !Whether optional inputs are needed
      if (ode_conf%lsode_dt_lim) then
cc        write (*,*) 'LSODE -- limiting time step'
        IOPT = 1                !Whether optional inputs are needed
        RWORK(6) = dtime        !Maximum time step allowed (optional input)
      endif
      if (ITASK==4) then
        IOPT = 1                !Whether optional inputs are needed
        RWORK(1) = TOUT         !TCRIT: not to be overshot
        IWORK(6) = 10000        !Maximum number of steps
      endif

      MF = 10       !Specifies non-stiff problem: implicit Adams-Moulton method

cc      write (*,*)'Number of equations',NEQD
cc      write (*,*)'Input state=',YY(1:NEQD(1))
cc      write (*,*)'Current time=',T
cc      write (*,*)'Next time=',TOUT
cc      write (*,*)'Tolerances=',ITOL,LS_RTOL,LS_ATOL
cc      write (*,*)'LSODE task=',ITASK
cc      write (*,*)'LSODE state=',ISTATE
cc      write (*,*)'Optional inputs required?',IOPT
cc      write (*,*)'Real work array size=',size(RWORK),' declared=',LRW
cc      write (*,*)'Integer work array size=',size(IWORK),' declared=',LIW
cc      write (*,*)'LSODE method=',MF

      if (external_f) then
        CALL LSODE (lsode_f,NEQD(1),YY(1:NEQD(1)),T,TOUT
     .             ,ITOL,LS_RTOL,LS_ATOL
     .             ,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,lsode_f, MF)
      else
        CALL LSODE (ff_3eq ,NEQD(1),YY(1:NEQD(1)),T,TOUT
     .             ,ITOL,LS_RTOL,LS_ATOL
     .             ,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,ff_3eq  ,MF)
      endif

      select case(ISTATE)
      case(2)       !Uneventful
        continue
      case(-1)      !Too much work: call again
        ISTATE = 3
      case(-2)      !Rescale tolerances
        ISTATE = 3
        RTOL=RTOL*RWORK(14)
        ATOL=ATOL*RWORK(14)
      case default  !Exit LSODE
        write (*,*) 'Error code ISTATE=',ISTATE
        write (*,*) 'See lsode.f for instructions'
        write (*,*) 'Aborting...'
        stop
      end select

c     Postprocessing to find logical coordinates

      if (cartesian) then
        x = xold ; y = yold ; z = zold
        call evalXi(iiout,YY(1),YY(2),YY(3),x,y,z,iterr)

        if (iterr /= 0) then
          iterr = 6
          return
        endif

        if (x >= rmin) then
          iterr = -2
          ISTATE = 1   !Instruct LSODE to start anew for next time step
          cartesian = .false.
          sing_point= .false.
        endif
      else
        x = YY(1)
        y = YY(2)
        z = YY(3)
      endif

c     Advance time and arc length

      dt   = (abs(T)-time)   !For output
      ds   = (abs(T)-time)   !For output

      time = abs(T)
      ss   = abs(T)

      if (NEQD(1) > 3) lsode_out(1:NEQD(1)-3)=cc*YY(4:NEQD(1))

c     End program

      end subroutine lsode_step_3eq

c     vp_solver
c     #########################################################################
      subroutine vp_solver(atol,rtol,x,y,z,ierror)

c     ------------------------------------------------------------------------
c     Performs home-grown Crank-Nicolson integration
c     ------------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: ierror
        real(8) :: atol,rtol,x,y,z

c     Local variables

        integer :: iterr

c     Begin program

        dt = min(dtime,1.2*dt)

 1      iterr = 0

        !Cartesian step (near SP)
        if (cartesian) then

          call vp_step(nt,niter,atol,rtol,cc*dt,x,y,z,0,iterr,error)

          if (iterr == 0) iterr = -2  !Notify of Cartesian step

        !Second order splitting
        elseif (order==2) then

          half = 0.5*dt

c         Begin first half-step split

          call vp_step(nt,niter,atol,rtol,cc*half,x,y,z,1,iterr,error)

c         Begin whole step split

          call vp_step(nt,niter,atol,rtol,cc*dt  ,x,y,z,2,iterr,error)

c         Begin second half-step split

          call vp_step(nt,niter,atol,rtol,cc*half,x,y,z,1,iterr,error)

        !First-order integration (30% cheaper!)
        else

          call vp_step(nt,niter,atol,rtol,cc*dt,x,y,z,1,iterr,error)

          call vp_step(nt,niter,atol,rtol,cc*dt,x,y,z,2,iterr,error)

        endif

c       Find jacobian at mid-point (for time/arc-length integral)

        if (cartesian) then
          jach = 1d0
        else
          xh = 0.5*(xold+x)
          yh = 0.5*(yold+y)
          zh = 0.5*(zold+z)
          jach = evalJ(xh,yh,zh,iterr)

          if (iterr /= 0) iterr = 4
        endif

c       Check for errors

        ierror = iterr

        if (iterr /= 0 .and. iterr /= -2) then

          if (iterr < 4 .and. iterr > 0) then
            if (iiout == 1) then
              write (*,*) 'VP -- Step',iterr
     .                ,' in CN did not converge to prescribed accuracy'
            elseif (iiout == 2) then
              write(*,*)'VP --Initial error=',error(1)
              write(*,*)'VP --Final error  =',error(niter),'> tol=',atol
            elseif (iiout > 2) then
              write (*,*)'VP --Error history=',error
            endif
          elseif (iterr == 4) then
            call pstop('vp_solver','Negative Jacobian')
          else
            write (messg,*) 'Unknown error ',iterr
            call pstop('orbit_step',messg)
          endif

          !Return to initial conditions and try again with smaller time step
          x = xold
          y = yold
          z = zold

          dt = 0.8*dt

          if (dt < 1d-8) call pstop('vp_solver','Time step too small')

          if (iiout > 1)
     .         write(*,*)'VP -- Trying smaller time step:',dt,'...'

          restart = .true.

          go to 1

        elseif (restart) then

          restart = .false.

          ierror = -1

          if (iiout > 1) then
            write (*,*) 'VP -- CN succeeded! Proceeding.'
            write (*,*)
          endif

        endif

c       Advance time and arc length

        call updt_measures(xh,yh,zh,jach*dt)

      end subroutine vp_solver

c     updt_measures
c     #########################################################################
      subroutine updt_measures(x,y,z,dt)

      real(8),intent(IN) :: x,y,z,dt

      real(8) :: vx,vy,vz

      !Update time step
      time=time+dt  !Include Jacobian factor to advance time

      !Update arc length
      call evalBcar(x,y,z,vx,vy,vz)

      ds = dt*sqrt(vx**2+vy**2+vz**2)
      ss = ss + ds

      end subroutine updt_measures

c     vp_step
c     ####################################################################
      subroutine vp_step(newton,niter,atol,rtol,dt,x,y,z,flag,iterr
     .                  ,error)

c     --------------------------------------------------------------------
c     Performs orbit integration step in a 3D simplectic manner (i.e.,
c     using only one component of vector potential at a time; see Finn &
c     Chacon, PoP 2005).  Currently allows both a Picard and a Newton
c     iteration for nonlinear convergence.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: flag,iterr,niter
      real(8) :: x,y,z,dt,error(niter),rtol,atol
      logical :: newton

c     Local variables

      integer :: iter

      integer,parameter :: size=3,icol=1

      real(8) :: tol
      real(8) :: JJ(size,size),res(size,icol),dx(size,icol)

c     Begin program

      !Check whether we are doing correct step (based on vector potential components)
      if (iterr /= 0 .or. (istep /= 0 .and. flag /= istep)) return

      !Initialize position
      x0=x
      y0=y
      z0=z

      if (iiout > 1) write (*,*) 'vp_step ',flag,' -- X0=',x0,y0,z0

cc      if (newton) then
cc        !Form Jacobian
cc        JJ = formJacobian(size,flag,dt,x0,y0,z0)
cc        call blockSolve(size,JJ,icol,res,dx,ret_inv=.true.)
cc      endif

      !Form initial residual
      res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z)

      tol = atol + rtol*sqrt(sum(res**2))

      !Iterations for first v1 half-step split
      do iter=1,niter

         !Find update
         if (newton) then

           if (iiout > 2) write (*,*) 'vp_step ',flag,' -- res=',res

           xh = 0.5*(x+x0)
           yh = 0.5*(y+y0)
           zh = 0.5*(z+z0)
           JJ = formJacobian(size,flag,dt,xh,yh,zh)
           call blockSolve(size,JJ,icol,res,dx)
cc           dx = matmul(JJ,res)

         else
           dx = res
         endif

         if (iiout > 2) write (*,*) 'vp_step ',flag,' -- dx=',dx

         !Update solution
         x = x + dx(1,1)
         y = y + dx(2,1)
         z = z + dx(3,1)

         if (iiout > 2) write (*,*) 'vp_step ',flag,' -- X =',x,y,z

         !Check convergence
         error(iter) = sqrt(sum(dx**2))

         if (iiout > 2) write (*,*) 'vp_step ',flag,' -- error='
     .                               ,error(iter)

         if (error(iter) < tol) exit

         !Form new residual
         res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z)

      enddo

      if (iiout == 2) then
        write (*,*) 'vp_step ',flag,' -- X =',x,y,z
     .             ,'; error history=',error(1:min(iter,niter))
      endif

      if (iter > niter) iterr = flag

c     End program

      end subroutine vp_step

c     formResidual
c     ###################################################################
      function formResidual(flag,dt,x0,y0,z0,x1,y1,z1) result(res)

      implicit none

      integer :: flag
      real(8) :: dt,x0,y0,z0,x1,y1,z1,res(3)

      real(8) :: xh,yh,zh,vx,vy,vz,xc0,yc0,zc0,x11,y11,z11

      !Find mid point
      xh = 0.5*(x1+x0)
      yh = 0.5*(y1+y0)
      zh = 0.5*(z1+z0)

      !Find vector field at midpoint
      call getB(xh,yh,zh,vx,vy,vz,.true.,cartesian,flag=flag)

      !Evaluate residual
      if (cartesian) then

        !Find Cartesian coordinates
        call evalX(x0,y0,z0,xc0,yc0,zc0)

        !Update Cartesian coordinates
        xc0 = xc0 +dt*vx
        yc0 = yc0 +dt*vy
        zc0 = zc0 +dt*vz

        !Find logical coordinates
        x11 = x0 ; y11 = y0 ; z11 = z0
        call evalXi(iiout,xc0,yc0,zc0,x11,y11,z11,ierror)

        if (ierror /= 0) then
          call pstop('formResidual','Could not find logical coords')
        endif

        res(1) = x1 - x11
        res(2) = y1 - y11
        res(3) = z1 - z11

      else
        res(1) = x1 - x0 -dt*vx
        res(2) = y1 - y0 -dt*vy
        res(3) = z1 - z0 -dt*vz
      endif

      end function formResidual

c     formJacobian
c     ###################################################################
      function formJacobian(size,flag,dt,x1,y1,z1) result(JJ)

      implicit none

      integer :: size,flag
      real(8) :: dt,x1,y1,z1,JJ(size,size)

      call evalCurlA(x1,y1,z1,JJ(1,1),JJ(2,1),JJ(3,1),flag=flag
     .              ,idx=1,idy=0,idz=0)

      call evalCurlA(x1,y1,z1,JJ(1,2),JJ(2,2),JJ(3,2),flag=flag
     .              ,idx=0,idy=1,idz=0)

      call evalCurlA(x1,y1,z1,JJ(1,3),JJ(2,3),JJ(3,3),flag=flag
     .              ,idx=0,idy=0,idz=1)

      JJ =-0.5*dt*JJ

      JJ(1,1) = JJ(1,1) + 1d0
      JJ(2,2) = JJ(2,2) + 1d0
      JJ(3,3) = JJ(3,3) + 1d0

      end function formJacobian

      end subroutine orbit_step

c     lsode_orbit_f_3eq
c     #####################################################################
      subroutine lsode_orbit_f_3eq(neq,t,y,ydot,xx1,yy1,zz1,B2)

c     ---------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = lsode_orbit_f(y,t), using
c     either arc-length or B-weighed arc-length (determined at orbit setup).
c
c     xx1,yy1,zz1 are optional, and provide logical coordinates on output.
c     B2 is optional, and provides B.B
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)
      real(8),optional,INTENT(OUT) :: xx1,yy1,zz1,B2

c     Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,vx,vy,vz,xc0,yc0,zc0,bnorm,jac
      logical :: bwgt

c     Begin program

      if (neq /= 3) then
        write (*,*) 'Wrong number of Eqs. in lsode_orbit_f_3eq:',neq
        stop
      endif

      bwgt = ode_conf%ds_is_Bds

      !Find logical coordinates
      if (sing_point) then
        x1 = xold ; y1 = yold ; z1 = zold
        call evalXi(iiout,y(1),y(2),y(3),x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) ierror
          write (*,*)
     .         'Could not find logical coords in lsode_orbit_f_3eq'
          stop
        endif
      else
        x1 = y(1)
        y1 = y(2)
        z1 = y(3)
      endif

      !Get vector field (Cartesian if sing_point is .true.)
      call getB(x1,y1,z1,vx,vy,vz,sln,sing_point)

      ydot(1) = vx
      ydot(2) = vy
      ydot(3) = vz

      !Time and arc-length ODEs
      if (.not.sing_point) then
        call evalBcar(x1,y1,z1,vx,vy,vz)
        bnorm = vx**2+vy**2+vz**2
        if (.not.bwgt) bnorm = sqrt(bnorm)

        jac = evalJ(x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) 'Negative Jacobian in lsode_orbit_f_3eq'
          stop
        endif
      else
        jac = 1d0
        bnorm = vx**2+vy**2+vz**2
        if (.not.bwgt) bnorm = sqrt(bnorm)
      endif

      ydot = ydot/(bnorm*jac)

      if (PRESENT(xx1)) xx1 = x1
      if (PRESENT(yy1)) yy1 = y1
      if (PRESENT(zz1)) zz1 = z1
      if (PRESENT(B2)) then
        if (bwgt) then
          B2 = bnorm
        else
          B2 = bnorm*bnorm
        endif
      endif

c     End program

      end subroutine lsode_orbit_f_3eq

c     lsode_orbit_f
c     #####################################################################
      subroutine lsode_orbit_f(neq,t,y,ydot,xx1,yy1,zz1,B2)

c     ---------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = lsode_orbit_f(y,t), using
c     pseudo-time (dt=ds/J/B). The fourth equation integrates geometric
c     pseudo-time (dt'=ds/B), and the fifth either arc-length
c     or B-weighed arc-length (determined at orbit setup).
c
c     xx1,yy1,zz1 are optional, and provide logical coordinates on output.
c     B2 is optional, and provides B.B
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)
      real(8),optional,INTENT(OUT) :: xx1,yy1,zz1,B2

c     Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,vx,vy,vz,xc0,yc0,zc0,bnorm
      logical :: bwgt

c     Begin program

      if (neq /= 5) then
        write (*,*) 'Wrong number of Eqs. for lsode_orbit',neq
        stop
      endif

      bwgt = ode_conf%ds_is_Bds

      !Find logical coordinates
      if (sing_point) then
        x1 = xold ; y1 = yold ; z1 = zold
        call evalXi(iiout,y(1),y(2),y(3),x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) ierror
          write (*,*) 'Could not find logical coords in lsode_orbit_f'
          stop
        endif
      else
        x1 = y(1)
        y1 = y(2)
        z1 = y(3)
      endif

      !Get vector field (Cartesian if sing_point is .true.)
      call getB(x1,y1,z1,vx,vy,vz,sln,sing_point)

      ydot(1) = vx
      ydot(2) = vy
      ydot(3) = vz

      !Time and arc-length ODEs
      if (.not.sing_point) then
        ydot(4) = evalJ(x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) 'Position=',x1,y1,z1
          write (*,*) 'Jacobian=',ydot(4)
          write (*,*) 'Negative Jacobian in lsode_orbit_f'
          stop
        endif

        call evalBcar(x1,y1,z1,vx,vy,vz)
        bnorm = vx**2+vy**2+vz**2
        if (.not.bwgt) bnorm = sqrt(bnorm)

        ydot(5) = ydot(4)*bnorm
      else
        ydot(4) = 1d0

        bnorm = vx**2+vy**2+vz**2
        if (.not.bwgt) bnorm = sqrt(bnorm)

        ydot(5) = bnorm
      endif

      if (PRESENT(xx1)) xx1 = x1
      if (PRESENT(yy1)) yy1 = y1
      if (PRESENT(zz1)) zz1 = z1
      if (PRESENT(B2)) then
        if (bwgt) then
          B2 = bnorm
        else
          B2 = bnorm*bnorm
        endif
      endif

c     End program

      end subroutine lsode_orbit_f

c     orbit_progress
c     #################################################################
      subroutine orbit_progress(time,tmax,fwd,exit_loop)

c     -----------------------------------------------------------------
c     Shows progress of orbit integration in various ways, depending
c     on level of output
c     -----------------------------------------------------------------

        implicit none

        real(8) :: time,tmax
        logical :: fwd,exit_loop

        real(8) :: tprev=0d0
        save :: tprev

        if (iiout == 0) then
          if (abs(time-tprev) > 0.1*tmax .or. exit_loop) then
            if (cartesian_flag) then
              write (*,FMT="(a1)",ADVANCE="NO") 'C'
              cartesian_flag = .false.
            elseif (cross) then
              write (*,FMT="(a1)",ADVANCE="NO") 'X'
              cross = .false.
            else
              if (fwd) then
                write (*,FMT="(a1)",ADVANCE="NO") '+'
              else
                write (*,FMT="(a1)",ADVANCE="NO") '-'
              endif
            endif
            tprev = time
          endif
        endif

      end subroutine orbit_progress

c     orbit_error
c     #####################################################################
      function orbit_error(iterror) result(exitloop)

c     -----------------------------------------------------------------
c     Checks error in ODE integration
c     -----------------------------------------------------------------

        implicit none

        integer :: iterror
        logical :: exitloop

        exitloop = .false.

        select case(iterror)
        case(0)    !No error
          continue
        case(-1)   !Time step cut performed (not fatal)
          cross = .true.
        case(-2)   !Cartesian integration performed (not fatal)
          cartesian_flag = .true.   
        case(1)
          write (*,*) 'Time step too small. Skipping orbit'
          exitloop = .true.
        case(5)
          write (*,*) 'Negative Jacobian in Cartesian integration'
          stop
        case(6)
          write (*,*) 'Error in Cartesian integration around SP'
          exitloop = .true.
        case default
          write (*,*) 'Error in orbit_step routine=',iterror
          stop
        end select

      end function orbit_error

c     ss_init
c     #################################################################
      subroutine ss_init(ss_order,lxss,xsurf,lyss,ysurf,lzss,zsurf,rflx)

c     -----------------------------------------------------------------
c     Initializes surface of section (SS) internal variables
c     -----------------------------------------------------------------

        implicit none

        integer :: ss_order
        real(8) :: xsurf,ysurf,zsurf
        logical :: lxss,lyss,lzss
        logical,optional :: rflx

c     Begin program

cc        if (my_rank == 0 .and. iiout >= 0) write (*,*)

c     Initialize private variables

        ssorder = ss_order

        lssx = lxss
        lssy = lyss
        lssz = lzss

        if (PRESENT(rflx)) then
          flux = rflx
        else
          flux = .false.
        endif

c     Set position of SS

        xss0 = xsurf
        yss0 = ysurf
        zss0 = zsurf

c     Consistency check

        if (lssx .and. lssy .and. lssz) then
          call pstop('ss_init'
     .              ,'Cannot select all 3 surfaces of section!!')
        endif

      end subroutine ss_init

c     ss_io_init
c     #################################################################
      subroutine ss_io_init(file,rz_ss)

c     -----------------------------------------------------------------
c     Sets up XDRAW IO for surfaces of section
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        character(*) :: file
        logical, optional :: rz_ss

c     Local variables

        integer :: cpos
        character(50) :: drawin

c     Begin program

        ss_io = .true.   !Turn on SS XDRAW IO

        if (PRESENT(rz_ss)) then
          rzss = rz_ss
        else
          rzss = .false.
        endif

cc        cpos = len_trim(ssx)-4
        cpos = index(trim(file),'.bin') - 1

        ssx = 'ssx'//trim(file)
        ssy = 'ssy'//trim(file)
        ssz = 'ssz'//trim(file)
        ssrz= 'ssrz'//trim(file)

        if (lssx) then
          drawin = 'draw'//ssx(1:cpos+3)//'.in'
          call createDrawInGfile(1,1,trim(ssx),'X surface of section'
     .                       ,'y',(/1/),(/'z'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunitx=27
          open(unit=nunitx,file=trim(ssx),form='unformatted')

        endif

        if (lssy) then
          drawin = 'draw'//ssy(1:cpos+3)//'.in'
          call createDrawInGfile(1,1,trim(ssy),'Y surface of section'
     .                       ,'x',(/1/),(/'z'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunity=28
          open(unit=nunity,file=trim(ssy),form='unformatted')
        endif

        if (lssz) then
          drawin = 'draw'//ssz(1:cpos+3)//'.in'
          call createDrawInGfile(1,1,trim(ssz),'Z surface of section'
     .                       ,'x',(/1/),(/'y'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunitz=29
          open(unit=nunitz,file=trim(ssz),form='unformatted')
        endif

        if (rzss) then
          drawin = 'draw'//ssrz(1:cpos+4)//'.in'
          call createDrawInGfile(1,1,trim(ssrz),'RZ surface of section'
     .                       ,'r',(/1/),(/'z'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunitrz=30
          open(unit=nunitrz,file=trim(ssrz),form='unformatted')
        endif

      end subroutine ss_io_init

c     ss_io_end
c     #################################################################
      subroutine ss_io_end

c     -----------------------------------------------------------------
c     Closes XDRAW surface of section output files
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        if (lssx) close(nunitx)
        if (lssy) close(nunity)
        if (lssz) close(nunitz)
        if (rzss) close(nunitrz)

      end subroutine ss_io_end

c     sschk
c     #################################################################
      subroutine sschk(order,tv,to,tt,xvold,xold,x,yvold,yold,y
     .                ,zvold,zold,z,xss,yss,zss,tss)

c     -----------------------------------------------------------------
c     Checks whether orbit intersects specified surfaces of section
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer,intent(IN) :: order
        real(8),intent(IN) :: tv,to,tt
     .                       ,xvold,xold,x
     .                       ,yvold,yold,y
     .                       ,zvold,zold,z
        real(8),intent(OUT):: xss,yss,zss,tss

c     Local variables

        real(8) :: rr
        integer :: iwhich

c     Begin program

        xss = 0d0 ; yss = 0d0 ; zss = 0d0 ; tss = 0d0

c     Lower order of interpolation for SS if derivatives are too different

        iwhich = order

        if (to > tv) then
          if (tt > to) then
            if(    (   abs(xold-xvold)/(to-tv)
     .              >5*abs(x    -xold)/(tt-to))
     .         .or.(   abs(yold-yvold)/(to-tv)
     .              >5*abs(y    -yold)/(tt-to))
     .         .or.(   abs(zold-zvold)/(to-tv)
     .              >5*abs(z    -zold)/(tt-to)) ) iwhich = 1
          else
            iwhich = 1
          endif
        else
          iwhich = 1
        endif

c     Check for intersection with X surface of section

        if (lssx
     .      .and.chk_ss_pos(xss0,xold,x,bcond(1),xmin,xmax,xss)) then

          call sursec(yvold,zvold,xvold,tv
     .               ,yold ,zold ,xold ,to
     .               ,y    ,z    ,x    ,tt
     .               ,xss,npointx,yss,zss,tss,iwhich)

          if (ss_io) then
            call chk_pos(xss,yss,zss) !this resets logical coords
            write(nunitx) real(yss),real(zss)
          endif

          !Accumulate B diagnostics
          call B_diag_update(yss,y00,tss)
        endif

c     Check for intersection with Y surface of section

        if (lssy
     .      .and.chk_ss_pos(yss0,yold,y,bcond(3),ymin,ymax,yss)) then

          call sursec(zvold,xvold,yvold,tv
     .               ,zold ,xold ,yold ,to
     .               ,z    ,x    ,y    ,tt
     .               ,yss,npointy,zss,xss,tss,iwhich)

          if (ss_io) then
            call chk_pos(xss,yss,zss) !this resets logical coords
            
            if (flux) then
              write(nunity) real(0.5*xss**2),real(zss)
            else
              write(nunity) real(xss),real(zss)
            endif
          endif

          !Accumulate B diagnostics
          call B_diag_update(xss,x00,tss)
        endif

c     Check for intersection with Z surface of section

        if ((lssz.or.rzss)
     .      .and.chk_ss_pos(zss0,zold,z,bcond(5),zmin,zmax,zss)) then

          call sursec(xvold,yvold,zvold,tv
     .               ,xold ,yold ,zold ,to
     .               ,x    ,y    ,z    ,tt
     .               ,zss,npointz,xss,yss,tss,iwhich)

          if (ss_io) then
            if (rzss) then   !Find R-Z plane positions
              call evalX(xss,yss,zss,xcar,ycar,zcar) !Eval Cartesian coords; also resets logical coords
              rr = sqrt(xcar**2+ycar**2)
              write(nunitrz) real(rr),real(zcar)
            endif

            if (lssz) then
              call chk_pos(xss,yss,zss) !this resets logical coords

              if (flux) then
                write(nunitz) real(0.5*xss**2),real(yss) !In flux
              else
                write(nunitz) real(xss)       ,real(yss) !In radius
              endif
            endif
          endif

          !Accumulate B diagnostics
          call B_diag_update(xss,x00,tss)
        endif

      end subroutine sschk

c     chk_ss_pos
c     ##################################################################
      function chk_ss_pos(x0,x1,x2,bc,xmin,xmax,xss) result(found)

c     ------------------------------------------------------------------
c     Determines whether surface of section at x0 is between x1 and x2.
c     It automatically accounts for periodic BCs.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      logical :: found
      integer :: bc
      real(8) :: x0,x1,x2,xmin,xmax,xss

c     Local variables

      integer :: pp,p1,p2
      real(8) :: x11,x22,x00

c     Begin program

      if (bc == PER) then
        p1 = floor((x1-xmin)/(xmax-xmin))
        p2 = floor((x2-xmin)/(xmax-xmin))

        x00 = x0

        if (p1 == p2) then
          pp = p1
          if (pp >= 0) x00 = x0 +(xmax-xmin)
          if (pp <  0) x00 = x0 -(xmax-xmin)
        else
          if     (p2 > p1.and.p1 >= 0) then  !x2 > x1, and cross xmax periodic boundary
            pp = min(p1,p2)
            x00 = x0 +(xmax-xmin)
          elseif (p2 > p1.and.p1 <  0) then  !x2 > x1, and cross xmin periodic boundary
            pp = max(p1,p2)
            x00 = x0 -(xmax-xmin)
          elseif (p2 < p1.and.p1 >= 0) then  !x1 > x2, and cross xmin periodic boundary
            pp = max(p1,p2)
            x00 = x0 -(xmax-xmin)
          elseif (p2 < p1.and.p1 <  0) then  !x1 < x2, and cross xmax periodic boundary
            pp = min(p1,p2)
            x00 = x0 +(xmax-xmin)
          endif
        endif

        x11 = x1 - pp*(xmax-xmin)
        x22 = x2 - pp*(xmax-xmin)
      else
        x00 = x0
        x11 = x1
        x22 = x2
        pp = 0
      endif

      found =    ((x11-x00)*(x22-x00) <= 0d0)
     .       .or.((x11-x0 )*(x22-x0 ) <= 0d0)

      if ((x11-x00)*(x22-x00) <= 0d0) xss = x00+pp*(xmax-xmin)
      if ((x11-x0 )*(x22-x0 ) <= 0d0) xss = x0 +pp*(xmax-xmin)

      end function chk_ss_pos

c     sursec
c     ##################################################################
      subroutine sursec(xvold,yvold,zvold,tvold,xold,yold,zold,told
     $                 ,x,y,z,time,zsurf,npoints,xss,yss,tss,order)

c     ------------------------------------------------------------------
c     Surface of section routine. Allows linear (order=1) and quadratic
c     (order=2) interpolation of time orbits.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: xvold,yvold,zvold,tvold,xold,yold,zold,told,x,y,z
      real(8) :: time,zsurf,xss,yss,tss
      integer :: npoints,order

c     Local variables

      real(8) :: a11,a12,a21,a22,b11,b12,b21,b22,aaa,bbb,ccc,det

c     Begin program

      select case(order)
      case(1)

c       Linear interpolation. xvold,yvold,zvold not used.

        tss=told+(zsurf-zold)*(time-told)/(z-zold)
        xss=xold+(x-xold)*(tss-told)/(time-told)
        yss=yold+(y-yold)*(tss-told)/(time-told)

      case(2)

c       Quadratic interpolation. Expand about zold.  
c       Vandermonde matrix with zvold, zold, z
c       Careful if time stepping is made nonuniform!

        a11=tvold-told
        a12=(tvold-told)**2
        a21=time-told
        a22=(time-told)**2

        det=a11*a22-a12*a21

        b11=a22/det
        b12=-a12/det
        b21=-a21/det
        b22=a11/det

        bbb=b11*(zvold-zold)+b12*(z-zold)
        ccc=b21*(zvold-zold)+b22*(z-zold)
        aaa=zold

        tss=(zsurf-zold)/bbb
        tss=told+tss           !Linear time
        if (ccc /= 0d0) then   !Calculate quadratic correction
          if(zsurf.ge.zold) then
            tss=(-bbb+sqrt(bbb**2-4.*ccc*(zold-zsurf)))/(2.*ccc)
            tss=told+tss
          else
            tss=(-bbb-sqrt(bbb**2-4.*ccc*(zold-zsurf)))/(2.*ccc)
            tss=told+tss
          endif
        endif

        bbb=b11*(xvold-xold)+b12*(x-xold)
        ccc=b21*(xvold-xold)+b22*(x-xold)
        aaa=xold
        xss=aaa+bbb*(tss-told)+ccc*(tss-told)**2

        bbb=b11*(yvold-yold)+b12*(y-yold)
        ccc=b21*(yvold-yold)+b22*(y-yold)
        aaa=yold
        yss=aaa+bbb*(tss-told)+ccc*(tss-told)**2

      case default

        write (*,*)
     .     'Order of accuracy in SS interpolation not implemented'
        stop

      end select

      npoints=npoints+1

c     End program

      end subroutine sursec

      end module orbit

c ff
c #########################################################################
      subroutine ff(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Interface driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use orbit

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

c Begin program

      call lsode_orbit_f(neq,t,y,ydot)

c End program

      end subroutine ff

c ff_3eq
c #########################################################################
      subroutine ff_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Interface driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use orbit

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

c Begin program

      call lsode_orbit_f_3eq(neq,t,y,ydot)

c End program

      end subroutine ff_3eq

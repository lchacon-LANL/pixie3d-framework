
c module orbit
c ######################################################################
      module orbit

        use xdraw_io

        use spline_field

        use grid, ONLY: my_rank,pstop,messg

        !ODE solver
        type :: ode_solver
          integer :: niter  !Number of nonlinear iterations
          integer :: order  !Order of integration (1,2)
          real(8) :: tol    !Convergence tolerance
          real(8) :: dtime  !Time step
          real(8) :: tmax   !Maximum time
          real(8) :: smax   !Maximum arc length
          logical :: solen  !Whether use solenoidal version or LSODE
          logical :: picard !Whether to use picard or Newton
        end type ode_solver

        !Splines
        type :: splines
          integer :: order  !Order of splines
          integer :: a_gauge!Gauge for vector potential (1,2)
        end type splines

        !Orbit config
        type :: orbit_def
          type(ode_solver) :: ode_conf
          type(splines)    :: spline_conf
        end type orbit_def

        !LSODE public variables
        integer :: odeflag
        real(8) :: xold,yold,zold
        logical :: sing_point=.false.,sln,odepack=.false.

        integer :: iiout

        !Private variables

        type(orbit_def),private  :: orb_default
        type(ode_solver),private :: ode_conf

        integer,private :: itime,iterror,ssorder

        integer,private :: nunitx,nunity,nunitz,nunitrz
        integer,private :: npointx,npointy,npointz

        real(8),private :: tvold,told,tprev,time,tmax
     .                    ,svold,sold,sprev,ss  ,smax
        real(8),private :: x,y,z,xcar,ycar,zcar,xss0,yss0,zss0
     .                    ,xvold,yvold,zvold
        real(8),private :: xss,yss,zss,tss

        real(8),private :: xmin,xmax,ymin,ymax,zmin,zmax

        character(50),private :: ssx,ssy,ssz,ssrz

        logical,private :: exit_loop,cartesian=.false.,cross=.false.
     .                    ,cartesian_flag=.false.
        logical,private :: lssx,lssy,lssz,rzss,ss_io=.false.

cc        INTERFACE orbit_setup
cc          module procedure orbit_setup_file, orbit_setup_gmetric
cc        end INTERFACE

      contains

c     orbit_conf_def
c     ######################################################################
      subroutine orbit_conf_def

c     ----------------------------------------------------------------------
c     Default configuration for orbit integrator
c     ----------------------------------------------------------------------

      !ODE solver
      orb_default%ode_conf%tmax   = 0d0
      orb_default%ode_conf%smax   = 20000.
      orb_default%ode_conf%dtime  = 0.01
      orb_default%ode_conf%niter  = 20
      orb_default%ode_conf%tol    = 1d-8
      orb_default%ode_conf%order  = 1
      orb_default%ode_conf%solen  = .true.
      orb_default%ode_conf%picard = .true.

      !Splines
      orb_default%spline_conf%order = 2
      orb_default%spline_conf%a_gauge = 2

      end subroutine orbit_conf_def

c     orbit_setup
c     ######################################################################
      subroutine orbit_setup(nx,ny,nz,xx,yy,zz
     .                      ,bx,by,bz,bxcar,bycar,bzcar
     .                      ,jac,xmap,orb_conf)

c     ----------------------------------------------------------------------
c     Sets up orbit splines for orbit integrator.
c     ----------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: nx,ny,nz
        real(8) :: xx(nx),yy(ny),zz(nz)
        real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
        real(8) :: bxcar(nx,ny,nz),bycar(nx,ny,nz),bzcar(nx,ny,nz)
        real(8) :: xmap(nx,ny,nz,3),jac(nx,ny,nz)

        type(orbit_def),optional :: orb_conf

c     Begin program

        if (PRESENT(orb_conf)) then
          orb_default = orb_conf
        else
          call orbit_conf_def
        endif

        ode_conf = orb_default%ode_conf

c     Setup domain limits

        if (bcond(1) == PER) then
          xmin = xx(2)
          xmax = xx(nx)
        else
          xmin = xx(1)
          xmax = xx(nx)
        endif

        if (bcond(3) == PER) then
          ymin = yy(2)
          ymax = yy(ny)
        else
          ymin = yy(1)
          ymax = yy(ny)
        endif

        if (bcond(5) == PER) then
          zmin = zz(2)
          zmax = zz(nz)
        else
          zmin = zz(1)
          zmax = zz(nz)
        endif

c     Do splines

        call setupSplines(nx,ny,nz,xx,yy,zz
     .                   ,orb_default%spline_conf%order)

        if (ode_conf%solen) then
          call splineA(bx,by,bz,orb_default%spline_conf%a_gauge)
        else
          call splineB(bx,by,bz)
        endif

        call splineBcar(bxcar,bycar,bzcar)

        call splineJ(jac )
        call splineX(xmap)

      end subroutine orbit_setup

ccc     orbit_setup_gmetric
ccc     ######################################################################
cc      subroutine orbit_setup_gmetric(igrid,nx,ny,nz,bx,by,bz,orb_conf)
cc
ccc     ----------------------------------------------------------------------
ccc     Sets up orbit splines for orbit integrator.
ccc     ----------------------------------------------------------------------
cc
cc        use grid
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer :: igrid,nx,ny,nz
cc        real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
cc
cc        type(orbit_def),optional :: orb_conf
cc
ccc     Local variables
cc
cc        integer :: igx,igy,igz,ig,jg,kg
cc        real(8) :: xx(nx),yy(ny),zz(nz)
cc        real(8) :: bxcar(nx,ny,nz),bycar(nx,ny,nz),bzcar(nx,ny,nz)
cc
ccc     Begin program
cc
cc        igx = igrid
cc        igy = igrid
cc        igz = igrid
cc
cc        if (     (nx /= grid_params%nxv(igx)+2)
cc     .      .or. (ny /= grid_params%nyv(igy)+2)
cc     .      .or. (nz /= grid_params%nzv(igz)+2)) then
cc          call pstop('orbit_setup','Inconsistent grid level')
cc        endif
cc
cc        call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)
cc
cc        xx = grid_params%xx(ig:ig+nx-1)
cc        yy = grid_params%yy(jg:jg+ny-1)
cc        zz = grid_params%zz(kg:kg+nz-1)
cc
ccc     Setup domain limits
cc
cc        if (bcond(1) == PER) then
cc          xmin = xx(2)
cc          xmax = xx(nx)
cc        else
cc          xmin = xx(1)
cc          xmax = xx(nx)
cc        endif
cc
cc        if (bcond(3) == PER) then
cc          ymin = yy(2)
cc          ymax = yy(ny)
cc        else
cc          ymin = yy(1)
cc          ymax = yy(ny)
cc        endif
cc
cc        if (bcond(5) == PER) then
cc          zmin = zz(2)
cc          zmax = zz(nz)
cc        else
cc          zmin = zz(1)
cc          zmax = zz(nz)
cc        endif
cc
cc        bxcar = bx
cc        bycar = by
cc        bzcar = bz
cc        call transformVector(igx,igy,igz
cc     .                      ,0,nx-1,0,ny-1,0,nz-1
cc     .                      ,bxcar,bycar,bzcar,'cnv','car')
cc
cc
cc        !ODE configuration
cc        if (PRESENT(orb_conf)) then
cc          orb_default = orb_conf
cc        else
cc          call orbit_conf_def
cc        endif
cc
cc        ode_conf = orb_default%ode_conf
cc
ccc     Do splines
cc
cc        call setupSplines(nx,ny,nz,xx,yy,zz
cc     .                   ,orb_default%spline_conf%order)
cc
cc        if (ode_conf%solen) then
cc          call splineA(bx,by,bz,orb_default%spline_conf%a_gauge)
cc        else
cc          call splineB(bx,by,bz)
cc        endif
cc
cc        call splineBcar(bxcar,bycar,bzcar)
cc
cc        call splineJ(gmetric%grid(igrid)%jac)
cc        call splineX(gmetric%grid(igrid)%car)
cc
cc        if (my_rank == 0) then
cc          write (*,*) 'Orbit_setup: Set up splines.'
cc          write (*,*)
cc        endif
cc
cc      end subroutine orbit_setup_gmetric

c     orbit_destroy
c     ######################################################################
      subroutine orbit_destroy

c     ----------------------------------------------------------------------
c     Frees up memory from orbit integration
c     ----------------------------------------------------------------------

c     Call variables

c     Begin program

        call killSplines

      end subroutine orbit_destroy

c     orbit_find
c     ######################################################################
      subroutine orbit_find(iout,x0,y0,z0,fwd,line_int,ss_order)

c     ----------------------------------------------------------------------
c     Integrates orbit up to specified limit in arc length or pseudotime
c     (specified in module). On input:
c       * iout: level of diagnostic output
c       * x0,y0,z0: starting point of orbit
c       * fwd: whether to go forward or backward
c     Optional input:
c       * line_int: external routine to perform line integrals
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

        optional :: line_int

        !Interface for generic line integral routine (optional)
        INTERFACE
          function line_int(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: line_int
cc          logical, optional :: line_int
          logical :: fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          logical :: ex_loop
          end function line_int
        END INTERFACE

c     Local variables

        real(8) :: temp_s,xh,yh,zh,sh
        logical :: chk_ss,integral

c     Begin program

        iiout = iout

        integral = PRESENT(line_int)

        chk_ss = (.not.integral).and.ss_io  !Do SS only if not integrating and if IO has been set up

        tmax = ode_conf%tmax
        smax = ode_conf%smax

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

          if(lssx) then
            write(nunitx) real(y0),real(z0)
          endif
          if(lssy) then
            write(nunity) real(x0),real(z0)
          endif
          if(lssz) then
            write(nunitz) real(x0),real(y0)
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
          call orbit_step(x,y,z,time,ss,iterror,forward=fwd)

          exit_loop = orbit_error(iterror)

          !Line-orbit integration
          if (integral) then
            exit_loop = line_int(fwd,xold,x,yold,y,zold,z,sold,ss)

            if (exit_loop) then
              if (iiout >=0) write (*,*)
              exit
            endif
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

          !Reset logical coords according to BCs
          call chk_pos(x,y,z)

          !Show progress/terminate
          if (smax /= 0d0) then
            call orbit_progress(ss  ,sprev,smax,fwd,exit_loop)
          else
            call orbit_progress(time,tprev,tmax,fwd,exit_loop)
          endif

          if (exit_loop) then
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

        x    = x0
        y    = y0
        z    = z0

        xold = x0
        yold = y0
        zold = z0

        time = 0d0
        told = 0d0
        tprev= 0d0

        ss   = 0d0  !Initialize arc-length along line
        sold = 0d0
        sprev= 0d0

        npointx=0
        npointy=0
        npointz=0

        itime = 0

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
      subroutine orbit_step(x,y,z,time,ss,ierror,forward)

c     ------------------------------------------------------------------
c     Crank-Nicolson (Picard/Newton) leapfrog, symmetrized (or not).
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierror
      real(8) :: x,y,z,time,ss

      logical,optional :: forward

c     Local Variables

      real(8) :: x0,y0,z0,x1,y1,z1,half,jach,dt=10.
      real(8) :: vx,vy,vz,error(ode_conf%niter),tol,dtime
      real(8) :: xh,yh,zh,xnew,ynew,znew,xc,yc,zc,ds,cc
      integer :: iter,iterr,order,niter
      logical :: restart,nt,fwd

      real(8) :: rmin=1d-3

      save dt

c     Begin program

c     ODE configuration

      dtime = ode_conf%dtime

      niter = ode_conf%niter
      tol   = ode_conf%tol
      order = ode_conf%order
      sln   = ode_conf%solen
      nt    = .not.ode_conf%picard

      if (PRESENT(forward)) then
        fwd = forward
      else
        fwd = .true.
      endif

      ierror = 0    !Initialize domain bounds error flag
      iterr  = 0    !Initialize time stepping error flag

      restart = .false. !Initialize restarting status flag

      if (fwd) then
        cc = 1d0
      else
        cc =-1d0
      endif

c     Check for SP

      cartesian = (bcond(1)==SP).and.(x<rmin)

c     Time integration

      if (sln .and. (.not.cartesian)) then  !Solenoidal

        if (odepack) then
          call lsode_step(tol,x,y,z,1,iterr)
          call lsode_step(tol,x,y,z,2,iterr)
        else
          call cn_solver
        endif

      else                                  !LSODE (non-solenoidal)

        call lsode_step(tol,x,y,z,0,iterr)

      endif

c     Output

      if (iiout > 0) then
        write (*,'(a,e10.3,a,e10.3,a,e10.3,a,e10.3,a,1p3e12.3)')
     .         ' Time=',time
     .        ,' ;  dt=',dt
     .        ,' ;  Arc length=',ss
     .        ,' ;  ds=',ds
     .        ,' ;  Position=',x,y,z
      endif

c     End program

      contains

c     cn_solver
c     #########################################################################
      subroutine cn_solver

c     ------------------------------------------------------------------------
c     Performs home-grown Crank-Nicolson integration
c     ------------------------------------------------------------------------

        implicit none

c     Begin program

        dt = min(dtime,1.2*dt)

        !Cartesian step (near SP)
 1      if (cartesian) then

          call simp_step(nt,niter,tol,dt,x,y,z,0,iterr,error)

          if (iterr == 0) iterr = -2  !Notify of Cartesian step

        !Second order splitting
        elseif (order==2) then

          half = 0.5*dt

c         Begin first half-step split

          call simp_step(nt,niter,tol,half,x,y,z,1,iterr,error)

c         Begin whole step split

          call simp_step(nt,niter,tol,dt  ,x,y,z,2,iterr,error)

c         Begin second half-step split

          call simp_step(nt,niter,tol,half,x,y,z,1,iterr,error)

        !First-order integration (30% cheaper!)
        else

          call simp_step(nt,niter,tol,dt,x,y,z,1,iterr,error)

          call simp_step(nt,niter,tol,dt,x,y,z,2,iterr,error)

        endif

c       Find jacobian at mid-point (for time/arc-length integral)

        if (cartesian) then
          jach = 1d0
        else
          xh = 0.5*(xold+x)
          yh = 0.5*(yold+y)
          zh = 0.5*(zold+z)
          jach = evalJ(xh,yh,zh,ierror)

          if (ierror /= 0) iterr = 4
        endif

c       Check for errors

        ierror = iterr

        if (iterr /= 0 .and. iterr /= -2) then

          if (iterr < 4 .and. iterr > 0) then
            if (iiout == 1) then
              write (*,*) 'Step',iterr
     .                ,' in CN did not converge to prescribed accuracy'
            elseif (iiout == 2) then
              write(*,*) 'Initial error=',error(1)
              write(*,*) 'Final error  =',error(niter),'> tol=',tol
            elseif (iiout > 2) then
              write (*,*) 'Error history=',error
            endif
          elseif (iterr == 4) then
            if (iiout > 1) write (*,*) 'Negative Jacobian'
          else
            write (messg,*) 'Unknown error ',iterr
            call pstop('orbit_step',messg)
          endif

          !Return to initial conditions and try again with smaller time step
          x = xold
          y = yold
          z = zold

          dt = 0.8*dt

          if (dt < 1d-8) then
            ierror = 1
            return
          endif

          if (iiout>1) write(*,*)'Trying smaller time step:',dt,'...'

          restart = .true.

          iterr = 0

          go to 1

        elseif (restart) then

          restart = .false.

          ierror = -1

          if (iiout > 1) then
            write (*,*) 'CN succeeded! Proceeding.'
            write (*,*)
          endif

        endif

c       Advance time and arc length

        call updt_measures(xh,yh,zh,jach*dt,time,ds,ss)

      end subroutine cn_solver

c     updt_measures
c     #########################################################################
      subroutine updt_measures(x,y,z,dt,time,ds,ss)

      real(8),intent(IN) :: x,y,z,dt
      real(8),intent(OUT):: time,ss,ds

      real(8) :: vx,vy,vz

      !Update time step
      time=time+dt  !Include Jacobian factor to advance time

      !Update arc length
      call evalBcar(x,y,z,vx,vy,vz)

      ds = dt*sqrt(vx**2+vy**2+vz**2)
      ss = ss + ds

      end subroutine updt_measures

c     lsode_step
c     ####################################################################
      subroutine lsode_step(tol,x,y,z,flag,iterr)

c     --------------------------------------------------------------------
c     Performs orbit integration step in a 3D simplectic manner (i.e.,
c     using only one component of vector potential at a time; see Finn &
c     Chacon, PoP 2005). Currently employs a Picard iteration for
c     nonlinear convergence. This should be improved to use Newton.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: flag,iterr
      real(8) :: x,y,z,tol

c     Local variables

      real(8) :: xc,yc,zc

      !LSODE setup
      INTEGER :: IOPT, IOUT, ISTATE=0, ITASK, ITOL, IWORK(20), LIW, LRW,
     .           MF, NEQ, NEQD(1)
      PARAMETER  (NEQ=5,LRW = 20 + 16*NEQ,LIW = 20)
      REAL(8) :: ATOL(NEQ),RTOL(1),RWORK(LRW),T,TOUT,YY(NEQ)

      EXTERNAL :: FF

      save ISTATE

c     Begin program

      if (cartesian .and. (.not.sing_point) ) ISTATE=1 !Restart LSODE for cartesian step

      sing_point = cartesian

c     Initial conditions

      if (cartesian) then
        !Find Cartesian coordinates
        call evalX(x,y,z,xc,yc,zc)

        YY(1) = xc
        YY(2) = yc
        YY(3) = zc
        YY(4) = time
        YY(5) = ss

        TOUT = 0.1*cc*dtime !Output time value (ignored after first call, since ITASK=2)
      else
        YY(1) = x
        YY(2) = y
        YY(3) = z
        YY(4) = time
        YY(5) = ss
        TOUT  = cc*dtime    !Output time value (ignored after first call, since ITASK=2)
      endif

c     LSODE configuration and call

      NEQD = NEQ

      T    = 0d0    !Independent variable

      ITOL = 2      !Whether atol is array or not (2 => true)
      RTOL = 1d-8   !Relative tolerance
      ATOL = tol    !Absolute tolerance
      ITASK = 2     !Mode of operation: take one step and return

      !Whether it is first call (ISTATE=1) or subsequent call (ISTATE=2)
      if ((time == 0d0).or. (sln .and. (.not.cartesian))) ISTATE = 1
cc      if ((time == 0d0)) ISTATE = 1
cc      if (.not.cartesian) ISTATE=1

      IOPT = 0      !Whether optional inputs are needed

      MF = 10       !Specifies non-stiff problem: implicit Adams method

      odeflag = flag   !Whether we use split step (1,2) or full step (0)

      CALL LSODE (ff, NEQD, YY, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     .            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, FF, MF)

      if (ISTATE /= 2) then
        write (*,*) 'Problems in LSODE: ISTATE=', ISTATE
        stop
      endif

c     Postprocessing to find logical coordinates

      if (cartesian) then
        x = xold ; y = yold ; z = zold
        call evalXi(iiout,YY(1),YY(2),YY(3),x,y,z,ierror)

        if (ierror /= 0) then
          ierror = 6
          return
        endif

        if (x > rmin) then
          ierror = -2
          ISTATE = 1   !Instruct LSODE to start anew for next time step
        endif
      else
        x = YY(1)
        y = YY(2)
        z = YY(3)

        call chk_pos(x,y,z)
      endif

c     Advance time and arc length

      if (flag == 0 .or. flag == 2) then
        ds   = YY(5)-ss
        dt   = YY(4)-time

        time = cc*YY(4)
        ss   = cc*YY(5)
      endif

c     End program

      end subroutine lsode_step

c     simp_step
c     ####################################################################
      subroutine simp_step(newton,niter,tol,dt,x,y,z,flag,iterr,error)

c     --------------------------------------------------------------------
c     Performs orbit integration step in a 3D simplectic manner (i.e.,
c     using only one component of vector potential at a time; see Finn &
c     Chacon, PoP 2005). Currently employs a Picard iteration for
c     nonlinear convergence. This should be improved to use Newton.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: flag,iterr,niter
      real(8) :: x,y,z,dt,error(niter),tol
      logical :: newton

c     Local variables

      integer :: iter

      integer,parameter :: size=3,icol=1

      real(8) :: JJ(size,size),res(size,icol),dx(size,icol)

c     Begin program

      if (iterr /= 0) return

      x0=x
      y0=y
      z0=z

      if (iiout > 1) write (*,*) 'simp_step ',flag,' -- X0=',x0,y0,z0

cc      if (newton) then
cc        !Form Jacobian
cc        JJ = formJacobian(size,flag,dt,x0,y0,z0)
cc        call blockSolve(size,JJ,icol,res,dx,ret_inv=.true.)
cc      endif

      !Iterations for first v1 half-step split
      do iter=1,niter

         !Form residual
         res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z)

         if (iiout > 2) write (*,*) 'simp_step ',flag,' -- res=',res

         !Find update
         if (newton) then
           xh = 0.5*(x+x0)
           yh = 0.5*(y+y0)
           zh = 0.5*(z+z0)
           JJ = formJacobian(size,flag,dt,xh,yh,zh)
           call blockSolve(size,JJ,icol,res,dx)
cc           dx = matmul(JJ,res)
         else
           dx = res
         endif

         if (iiout > 2) write (*,*) 'simp_step ',flag,' -- dx=',dx

         !Check convergence
         error(iter) = sqrt(sum(dx**2))

         if (iiout > 2) write (*,*) 'simp_step ',flag,' -- error='
     .                               ,error(iter)

         if (error(iter) < tol) exit

         !Update solution
         x = x + dx(1,1)
         y = y + dx(2,1)
         z = z + dx(3,1)

         if (iiout > 2) write (*,*) 'simp_step ',flag,' -- X =',x,y,z

      enddo

      if (iiout > 1) then
        write (*,*) 'simp_step ',flag,' -- X =',x,y,z
     .             ,'; error=',error(min(iter,niter))
cc        pause
      endif

      if (iiout > 2) write (*,*) 'simp_step ',flag,' -- error history'
     .                            ,error(1:min(iter,niter))

      if (iter > niter) iterr = flag

c     End program

      end subroutine simp_step

c     formResidual
c     ###################################################################
      function formResidual(flag,dt,x0,y0,z0,x1,y1,z1) result(res)

      implicit none

      integer :: flag
      real(8) :: dt,x0,y0,z0,x1,y1,z1,res(3)

      real(8) :: cdt,xh,yh,zh,vx,vy,vz,xc0,yc0,zc0,x11,y11,z11

      cdt = cc*dt

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
        xc0 = xc0 +cdt*vx
        yc0 = yc0 +cdt*vy
        zc0 = zc0 +cdt*vz

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
        res(1) = x1 - x0 -cdt*vx
        res(2) = y1 - y0 -cdt*vy
        res(3) = z1 - z0 -cdt*vz
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

      JJ =-0.5*cc*dt*JJ

      JJ(1,1) = JJ(1,1) + 1d0
      JJ(2,2) = JJ(2,2) + 1d0
      JJ(3,3) = JJ(3,3) + 1d0

      end function formJacobian

      end subroutine orbit_step

c     orbit_progress
c     #################################################################
      subroutine orbit_progress(time,tprev,tmax,fwd,exit_loop)

c     -----------------------------------------------------------------
c     Shows progress of orbit integration in various ways, depending
c     on level of output
c     -----------------------------------------------------------------

        implicit none

        real(8) :: time,tprev,tmax
        logical :: fwd,exit_loop

        exit_loop = (abs(time) > tmax)

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
c     Checks error in time integration
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
      subroutine ss_init(ss_order,lxss,xsurf,lyss,ysurf,lzss,zsurf)

c     -----------------------------------------------------------------
c     Initializes surface of section (SS) internal variables
c     -----------------------------------------------------------------

        implicit none

        integer :: ss_order
        real(8) :: xsurf,ysurf,zsurf
        logical :: lxss,lyss,lzss

c     Begin program

        if (my_rank == 0 .and. iiout >= 0) write (*,*)

c     Initialize private variables

        ssorder = ss_order

        lssx = lxss
        lssy = lyss
        lssz = lzss

c     Set position of SS

        xss0 = xsurf
        yss0 = ysurf
        zss0 = zsurf

c     Consistency check

        if ((.not.lssx) .and. (.not.lssy) .and. (.not.lssz)) then
          call pstop('ss_init'
     .              ,'No Surface of section selected. Aborting.')
        elseif (lssx .and. lssy .and. lssz) then
          call pstop('ss_init'
     .              ,'Cannot select all 3 surfaces of section!!')
        endif

      end subroutine ss_init

c     ss_io_init
c     #################################################################
      subroutine ss_io_init(file,rz_ss)

c     -----------------------------------------------------------------
c     Sets up XDRAW surface of section output files
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        character(*) :: file
        logical, optional :: rz_ss

c     Local variables

        character(50) :: drawin

c     Begin program

        ss_io = .true.   !Turn on SS XDRAW IO

        if (PRESENT(rz_ss)) then
          rzss = rz_ss
        else
          rzss = .false.
        endif

cc        trace = 'trc_'//trim(file)
        ssx = 'ssx'//trim(file)
        ssy = 'ssy'//trim(file)
        ssz = 'ssz'//trim(file)
        ssrz= 'ssrz'//trim(file)

        if (lssx) then
          drawin = 'draw'//ssx(1:len_trim(ssx)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssx),'X surface of section'
     .                       ,'y',(/1/),(/'z'/),(/0/),(/0/)
     .                       ,drawin,.false.,.false.)
          nunitx=27
          open(unit=nunitx,file=trim(ssx),form='unformatted')
        endif

        if (lssy) then
          drawin = 'draw'//ssy(1:len_trim(ssy)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssy),'Y surface of section'
     .                       ,'x',(/1/),(/'z'/),(/0/),(/0/)
     .                       ,drawin,.false.,.false.)
          nunity=28
          open(unit=nunity,file=trim(ssy),form='unformatted')
        endif

        if (lssz) then
          drawin = 'draw'//ssz(1:len_trim(ssz)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssz),'Z surface of section'
     .                       ,'x',(/1/),(/'y'/),(/0/),(/0/)
     .                       ,drawin,.false.,.false.)
          nunitz=29
          open(unit=nunitz,file=trim(ssz),form='unformatted')
        endif

        if (rzss) then
          drawin = 'draw'//ssrz(1:len_trim(ssrz)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssrz),'RZ surface of section'
     .                       ,'r',(/1/),(/'z'/),(/0/),(/0/)
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

        iwhich = order

c     Lower order of interpolation for SS if derivatives are too different

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
        endif

c     Check for intersection with X surface of section

        if (lssx) then
          !Check for periodic BCs (suggested by Daniele Bonfiglio, 4/8/09)
          if ((x-xss0)*(xold-xss0) <= 0.0) then
            xss = xss0
          elseif ((x-xmax)*(xold-xmax) <= 0.0
     .           .and. bcond(1)==PER
     .           .and. xss0==xmin) then
            xss = xmax
          elseif ((x-xmin)*(xold-xmin) <= 0.0
     .           .and. bcond(1)==PER
     .           .and. xss0==xmax) then
            xss = xmin
          else
            goto 2
          endif

          call sursec(yvold,zvold,xvold,tv
     .               ,yold ,zold ,xold ,to
     .               ,y    ,z    ,x    ,tt
     .               ,xss,npointx,yss,zss,tss,iwhich)

          call chk_pos(xss,yss,zss) !this resets logical coords

          if (ss_io) write(nunitx) real(yss),real(zss)
        endif

c     Check for intersection with Y surface of section

 2      if (lssy) then
          !Check for periodic BCs (suggested by Daniele Bonfiglio, 4/8/09)
          if((y-yss0)*(yold-yss0) <= 0.0) then
            yss = yss0
          elseif ((y-ymax)*(yold-ymax) <= 0.0
     .           .and. bcond(3)==PER
     .           .and. yss0==ymin) then
            yss = ymax
          elseif ((y-ymin)*(yold-ymin) <= 0.0
     .           .and. bcond(3)==PER
     .           .and. yss0==ymax) then
            yss = ymin
          else
            goto 3
          endif

          call sursec(zvold,xvold,yvold,tv
     .               ,zold ,xold ,yold ,to
     .               ,z    ,x    ,y    ,tt
     .               ,yss,npointy,zss,xss,tss,iwhich)

          call chk_pos(xss,yss,zss) !this resets logical coords

          if (ss_io) write(nunity) real(xss),real(zss)
        endif

c     Check for intersection with Z surface of section

 3      if (lssz) then
          !Check for periodic BCs (suggested by Daniele Bonfiglio, 4/8/09)
          if((z-zss0)*(zold-zss0) <= 0.0) then
            zss = zss0
          elseif ((z-zmax)*(zold-zmax) <= 0.0
     .           .and. bcond(5)==PER
     .           .and. zss0==zmin) then
            zss = zmax
          elseif ((z-zmin)*(zold-zmin) <= 0.0
     .           .and. bcond(5)==PER
     .           .and. zss0==zmax) then
            zss = zmin
          else
            return
          endif

          call sursec(xvold,yvold,zvold,tv
     .               ,xold ,yold ,zold ,to
     .               ,x    ,y    ,z    ,tt
     .               ,zss,npointz,xss,yss,tss,iwhich)

          if (rzss.and.ss_io) then   !Find R-Z plane positions
            call evalX(xss,yss,zss,xcar,ycar,zcar) !Eval Cartesian coords; also resets logical coords
            rr = sqrt(xcar**2+ycar**2)
            write(nunitrz) real(rr),real(zcar)
cc            write(*,*) real(rr),real(zcar)   !diag
          else
            call chk_pos(xss,yss,zss) !this resets logical coords
          endif

          if (ss_io) write(nunitz) real(xss),real(yss)
        endif

      end subroutine sschk

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
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use orbit, ONLY: spt=>sing_point,sln,odeflag,iiout,xold,yold,zold
      use spline_field

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,vx,vy,vz,xc0,yc0,zc0

c Begin program

      !Find logical coordinates
      if (spt) then
        x1 = xold ; y1 = yold ; z1 = zold
        call evalXi(iiout,y(1),y(2),y(3),x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) ierror
          write (*,*) 'Could not find logical coords in ff'
          stop
        endif
      else
        x1 = y(1)
        y1 = y(2)
        z1 = y(3)
      endif

      !Get vector field (Cartesian if spt is .true.)
      call getB(x1,y1,z1,vx,vy,vz,sln,spt,flag=odeflag)

      ydot(1) = vx
      ydot(2) = vy
      ydot(3) = vz

      !Time and arc-length ODEs
      if (.not.spt) then
        ydot(4) = evalJ(x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) 'Negative Jacobian in ff'
          stop
        endif

        call evalBcar(x1,y1,z1,vx,vy,vz)
        ydot(5) = ydot(4)*sqrt(vx**2+vy**2+vz**2)
      else
        ydot(4) = 1d0
        ydot(5) = sqrt(vx**2+vy**2+vz**2)
      endif

c End program

      end subroutine ff

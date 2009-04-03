
c module orbit
c ######################################################################
      module orbit

        use xdraw_io

        use spline_field

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

        logical :: sing_point=.false.

        !Private variables

        type(ode_solver),private :: ode_conf

        integer,private :: ilevel,itime,iterror,ssorder

        integer,private :: nunitx,nunity,nunitz,nunitrz
        integer,private :: npointx,npointy,npointz

        real(8),private :: tvold,told,tprev,time,tmax
     .                    ,svold,sold,sprev,ss  ,smax
        real(8),private :: x,y,z,xcar,ycar,zcar,xss0,yss0,zss0
     .                    ,xvold,yvold,zvold,xold,yold,zold
        real(8),private :: xss,yss,zss,tss

        character(50),private :: ssx,ssy,ssz,ssrz

        logical,private :: exit_loop,cartesian=.false.,cross=.false.
        logical,private :: lssx,lssy,lssz,rzss

cc        INTERFACE orbit_setup
cc          module procedure orbit_setup_file
cc        end INTERFACE

      contains

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
     .            ,ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
        real(8) :: bxcar(nx,ny,nz),bycar(nx,ny,nz),bzcar(nx,ny,nz)
        real(8) :: xmap(nx,ny,nz,3),jac(nx,ny,nz)

        type(orbit_def) :: orb_conf

c     Begin program

        ode_conf = orb_conf%ode_conf

c     Do splines

        call setupSplines(nx,ny,nz,xx,yy,zz,orb_conf%spline_conf%order)

        if (ode_conf%solen) then
          call splineA(bx,by,bz,orb_conf%spline_conf%a_gauge)
        else
          call splineB(bx,by,bz)
        endif

        call splineBcar(bxcar,bycar,bzcar)

        call splineJ(jac )
        call splineX(xmap)

        write (*,*) 'Set up splines.'
        write (*,*)

      end subroutine orbit_setup

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

        ilevel = iout

        integral = PRESENT(line_int)

        chk_ss = .not.integral

        tmax = ode_conf%tmax
        smax = ode_conf%smax

c     Initialize orbit

        call orbit_init(x0,y0,z0)

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

          call orbit_error(iterror)

          !Line-orbit integration
          if (integral) then
            exit_loop = line_int(fwd,xold,x,yold,y,zold,z,sold,ss)

            if (exit_loop) then
              if (.not.fwd) write (*,*)
              exit
            endif
          endif

          !Check for intersection with surface of section
          if (iterror /= 2 .and. itime > 2 .and. chk_ss) then
            if (smax /= 0d0) then
              call sschk(ssorder,svold,sold,ss  ,xss,yss,zss,tss)
            else
              call sschk(ssorder,tvold,told,time,xss,yss,zss,tss)
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
            write (*,*)
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

      end subroutine orbit_init

c     orbit_step
c     ##################################################################
      subroutine orbit_step(x,y,z,time,ss,ierror,forward)

c     ------------------------------------------------------------------
c     Crank-Nicolson (Picard/Newton) leapfrog, symmetrized (or not).
c     ------------------------------------------------------------------

      implicit none

c     Call variables

cc      type(ode_solver) :: ode_conf

      integer :: ierror
      real(8) :: x,y,z,time,ss

      logical,optional :: forward

c     Local Variables

      real(8) :: x0,y0,z0,x1,y1,z1,half,jach,dt=10.
      real(8) :: vx,vy,vz,error(ode_conf%niter),tol,dtime
      real(8) :: xh,yh,zh,xnew,ynew,znew,xc,yc,zc,ds,cc
      integer :: iter,iterr,order,niter
      logical :: restart,cartesian,nt,fwd,sln

      real(8) :: rmin=1d-3

      save dt

      !LSODE setup
      INTEGER :: IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,
     .           MF, NEQ, NEQD(1)
      PARAMETER  (NEQ=5,LRW = 20 + 16*NEQ,LIW = 20)
      REAL(8) :: ATOL(NEQ),RTOL(1),RWORK(20 + 16*NEQ),T,TOUT,YY(NEQ)

      EXTERNAL :: FF

c     Begin program

c     ODE configuration

      dtime  = ode_conf%dtime

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

cc      sing_point = cartesian
cc
cc      if (cartesian) then
cc
cc        !Find local magnetic field
cc        call evalCurlA(x,y,z,vx,vy,vz) 
cc
cccc        sing_point = (vx < 0d0)      !Perform SP integration only if vx < 0
cc
cc      endif

      if (cartesian) then

        if (ilevel > 0) then
          write (*,*)
          write (*,*) 'Performing Cartesian integration'
          write (*,*) 'time step -- initial xi=',x,y,z
          write (*,*)
        endif

        !Estimate time step
        jach = evalJ(x,y,z,ierror)

        if (ierror /= 0) then
          ierror = 5
          return
        endif

        call evalCurlA(x,y,z,vx,vy,vz) 

        if (cc*vx < 0d0) then    !Perform integration only if cc*vx < 0
          half = 2*abs(rmin*jach/vx)

          !Start time cycle
          do

            !Find Cartesian coordinates
            call evalX(x,y,z,xc,yc,zc)

            !Find Cartesian vector components
            call evalBcar(x,y,z,vx,vy,vz)

            !Update cartesian coordinates (Forward Euler, single step)
            xc=xc+cc*half*vx
            yc=yc+cc*half*vy
            zc=zc+cc*half*vz

            !Find logical coordinates
            call evalXi(xc,yc,zc,x,y,z,ierror)

            if (ierror /= 0) then
              ierror = 6
              return
            endif

cc            call updt_measures(half,time,ds,ss)

            !Update time step
            time=time+half  !Include Jacobian factor to advance time

            !Update arc length
            ds = half*sqrt(vx**2+vy**2+vz**2)
            ss = ss + ds

            !Check we are out of inner circle or radius rmin
            if (x > rmin) exit

          enddo

          if (ilevel > 0) then
            write (*,*)
            write (*,*) 'time step -- final xi=',x,y,z
            write (*,*) 'End Cartesian integration'
            write (*,*)
          endif

          ierror = -2

          return
        endif
      endif

c     Solenoidal time integration

      if (sln) then

        dt = min(dtime,1.2*dt)

        !Second order splitting
 1      if (order==2) then

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

        xh = 0.5*(xold+x)
        yh = 0.5*(yold+y)
        zh = 0.5*(zold+z)
        jach = evalJ(xh,yh,zh,ierror)

        if (ierror /= 0) iterr = 4

c       Check for errors

        if (iterr /= 0) then

          if (ilevel > 1) then

            if (iterr < 4) then
              write (*,*) 'Step',iterr
     .                ,' in CN did not converge to prescribed accuracy'

              if (ilevel > 2) then
                write (*,*) 'Error history=',error
              else
                write(*,*) 'Initial error=',error(1)
                write(*,*) 'Final error  =',error(niter),'> tol=',tol
              endif
            elseif (iterr == 4) then
              write (*,*) 'Negative Jacobian'
            else
              ierror = iterr
              return
            endif

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

          if (ilevel>0) write(*,*)'Trying smaller time step:',dt,'...'

          restart = .true.

          iterr = 0

          go to 1

        elseif (restart) then

          restart = .false.

          ierror = -1

          if (ilevel > 0) then
            write (*,*) 'CN succeeded! Proceeding.'
            write (*,*)
          endif

        endif

c       Advance time and arc length

        call updt_measures(xh,yh,zh,jach*dt,time,ds,ss)

      else

c       Standard time integration

        if (sing_point) then

          write (*,*) sing_point
          stop

          if (ilevel > 0) then
            write (*,*)
            write (*,*) 'Performing Cartesian integration'
            write (*,*) 'time step -- initial xi=',x,y,z
            write (*,*)
          endif

          !Find Cartesian coordinates
          call evalX(x,y,z,xc,yc,zc)

          YY(1) = xc
          YY(2) = yc
          YY(3) = zc
          YY(4) = time
          YY(5) = ss

          !Estimate time step
          jach = evalJ(x,y,z,ierror)

          if (ierror /= 0) then
            ierror = 5
            return
          endif

          TOUT = cc*2*abs(rmin*jach/vx) !Output time value (ignored after first call, since ITASK=2)

        else
          YY(1) = x
          YY(2) = y
          YY(3) = z
          YY(4) = time
          YY(5) = ss
          TOUT  = cc*dtime             !Output time value (ignored after first call, since ITASK=2)
        endif

        NEQD = NEQ

        T    = 0d0    !Independent variable

        ITOL = 2      !Whether atol is array or not (2 => true)
        RTOL = 1.D-8  !Relative tolerance
        ATOL = tol    !Absolute tolerance
        ITASK = 2     !Mode of operation: take one step and return
cc        ITASK = 1     !Mode of operation: take one step and return

        ISTATE = 1    !Whether it is first call (1) or subsequent call (2)
cc        if (time > 0d0) ISTATE = 2
        if (ss > 0d0) ISTATE = 2

        IOPT = 0      !Whether optional inputs are needed

        MF = 10       !Specifies non-stiff problem

        CALL LSODE (ff, NEQD, YY, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     .              ISTATE, IOPT, RWORK, LRW, IWORK, LIW, FF, MF)

        if (ISTATE /= 2) then
          write (*,*) 'Problems in LSODE: ISTATE=', ISTATE
          stop
        endif

        if (sing_point) then
          !Find logical coordinates
          call evalXi(YY(1),YY(2),YY(3),x,y,z,ierror)

          if (ierror /= 0) then
            ierror = 6
            return
          endif

          if (ilevel > 0 .and. x > rmin) then
            write (*,*)
            write (*,*) 'time step -- final xi=',x,y,z
            write (*,*) 'End Cartesian integration'
            write (*,*)

            ierror = -2
          endif

        else
          x = YY(1)
          y = YY(2)
          z = YY(3)

cc          call chk_pos(x,y,z)

        endif

        if (ilevel>=1) write(*,'(a,e10.3,a,e10.3,a,e10.3
     .                            ,a,f10.3,a,1p3e12.3)')
     .                       ' dt=',YY(4)-time
     .                      ,' ;  Time=',time
     .                      ,' ;  ds=',YY(5)-ss
     .                      ,' ;  Arc length=',ss
     .                      ,' ;  Position=',x,y,z

        time = cc*YY(4)
        ss   = cc*YY(5)

      endif

c     End program

      contains

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

      if (ilevel>=1) write(*,'(a,e10.3,a,e10.3,a,e10.3
     .                            ,a,f10.3,a,1p3e12.3)')
     .                       ' dt=',dt
     .                      ,' ;  Time=',time
     .                      ,' ;  ds=',ds
     .                      ,' ;  Arc length=',ss
     .                      ,' ;  Position=',x,y,z

      end subroutine updt_measures

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

      if (ilevel > 1) write (*,*) 'simp_step -- x0=',x0,y0,z0

cc      if (newton) then
cc        !Form Jacobian
cc        JJ = formJacobian(size,flag,dt,x0,y0,z0)
cc        call blockSolve(size,JJ,icol,res,dx,ret_inv=.true.)
cc      endif

        !Iterations for first v1 half-step split
      do iter=1,niter

         !Form residual
         res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z)

         if (ilevel > 1) write (*,*) 'simp_step -- res=',res

cc         !Check convergence
cc         error(iter) = sqrt(sum(res**2))
cc
cc         if (ilevel > 1) write (*,*) 'simp_step -- error=',error(iter)
cc
cc         if (error(iter) < tol) exit

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

         if (ilevel > 1) write (*,*) 'simp_step -- dx=',dx

         !Check convergence
         error(iter) = sqrt(sum(dx**2))

         if (ilevel > 1) write (*,*) 'simp_step -- error=',error(iter)

         if (error(iter) < tol) exit

         !Update solution
         x = x + dx(1,1)
         y = y + dx(2,1)
         z = z + dx(3,1)

         if (ilevel > 1) write (*,*) 'simp_step -- pos=',x,y,z

      enddo

cc      else
cc
cc        do iter=1,niter
cc           !Find mid point (logical)
cc           xh = 0.5*(x0+x)
cc           yh = 0.5*(y0+y)
cc           zh = 0.5*(z0+z)
cc
cc           !Find vector field at midpoint
cc           call evalCurlA(xh,yh,zh,vx,vy,vz,flag=flag)
cc
cc           !Find new position
cc           xnew=x0+cc*dt*vx
cc           ynew=y0+cc*dt*vy         !Leave out y eq if A_x=0 gauge.
cc           znew=z0+cc*dt*vz
cc
cc           !Check error
cc           error(iter) = sqrt((xnew-x)**2+(ynew-y)**2+(znew-z)**2)
cc
cc           !Update current position
cc           x = xnew
cc           y = ynew
cc           z = znew
cc
cc           if (error(iter) < tol) exit
cc        enddo
cc
cc      endif

      if (ilevel > 0) write (*,*) 'simp_step -- error history'
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

      real(8) :: xh,yh,zh,vx,vy,vz

      !Find mid point
      xh = 0.5*(x1+x0)
      yh = 0.5*(y1+y0)
      zh = 0.5*(z1+z0)

      !Find vector field at midpoint
      call evalCurlA(xh,yh,zh,vx,vy,vz,flag=flag)

      !Evaluate residual
      res(1) = x1 - x0 -cc*dt*vx
      res(2) = y1 - y0 -cc*dt*vy
      res(3) = z1 - z0 -cc*dt*vz

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

cc        exit_loop = .false.
cc
cc        if (abs(time) > tmax) then
cc          if (ilevel == 0) then
cc            if (cartesian) then
cc              write (*,FMT="(a1)") 'C'
cc              cartesian = .false.
cc            elseif (cross) then
cc              write (*,FMT="(a1)") 'X'
cc              cross = .false.
cc            else
cc              if (fwd) then
cc                write (*,FMT="(a1)") '+'
cc              else
cc                write (*,FMT="(a1)") '-'
cc              endif
cc            endif
cc          endif
cc          exit_loop = .true.
cc        elseif (ilevel == 0) then
cc          if (abs(time-tprev) > 0.1*tmax) then
cc            if (cartesian) then
cc              write (*,FMT="(a1)",ADVANCE="NO") 'C'
cc              cartesian = .false.
cc            elseif (cross) then
cc              write (*,FMT="(a1)",ADVANCE="NO") 'X'
cc              cross = .false.
cc            else
cc              if (fwd) then
cc                write (*,FMT="(a1)",ADVANCE="NO") '+'
cc              else
cc                write (*,FMT="(a1)",ADVANCE="NO") '-'
cc              endif
cc            endif
cc            tprev = time
cc          endif
cc        endif

        exit_loop = (abs(time) > tmax)

        if (ilevel == 0) then
          if (abs(time-tprev) > 0.1*tmax) then
            if (cartesian) then
              write (*,FMT="(a1)",ADVANCE="NO") 'C'
              cartesian = .false.
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
      subroutine orbit_error(iterror)

c     -----------------------------------------------------------------
c     Checks error in time integration
c     -----------------------------------------------------------------

        implicit none

        integer :: iterror

        select case(iterror)
        case(0)    !No error
          continue
        case(-1)   !Time step cut performed (not fatal)
          cross = .true.
        case(-2)   !Cartesian integration performed (not fatal)
          cartesian = .true.
        case(1)
          write (*,*) 'Time step too small. Skipping orbit'
          stop
        case(5)
          write (*,*) 'Negative Jacobian in Cartesian integration'
          stop
        case(6)
          write (*,*) 'Error in Cartesian integration around SP'
          stop
        case default
          write (*,*) 'Error in orbit_step routine=',iterror
          stop
        end select

      end subroutine orbit_error

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

        write (*,*)

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
          write (*,*) 'No Surface of section selected. Aborting.'
          stop
        elseif (lssx .and. lssy .and. lssz) then
          write (*,*) 'Cannot select all 3 surfaces of section!!'
          stop
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
      subroutine sschk(order,tv,to,tt,xss,yss,zss,tss)

c     -----------------------------------------------------------------
c     Checks whether orbit intersects specified surfaces of section
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer,intent(IN) :: order
        real(8),intent(IN) :: tv,to,tt
        real(8),intent(OUT):: xss,yss,zss,tss

c     Local variables

        real(8) :: rr
        integer :: iwhich

c     Begin program

        iwhich = order

        if (to > tv) then
          if(   abs(xold-xvold)/(to-tv)
     .       >5*abs(x    -xold)/(tt-to))
     .      iwhich = 1      !Lower order of interpolation for SS
        endif

        !Check for intersection with X surface of section
        if((x-xss0)*(xold-xss0) <= 0.0.and.lssx) then
          xss = xss0
          call sursec(yvold,zvold,xvold,tv
     .               ,yold ,zold ,xold ,to
     .               ,y    ,z    ,x    ,tt
     .               ,xss,npointx,yss,zss,tss,iwhich)

          call chk_pos(xss,yss,zss) !this resets logical coords

          write(nunitx) real(yss),real(zss)
        endif

        !Check for intersection with Y surface of section
        if((y-yss0)*(yold-yss0) <= 0.0.and.lssy) then
          yss = yss0

          call sursec(zvold,xvold,yvold,tv
     .               ,zold ,xold ,yold ,to
     .               ,z    ,x    ,y    ,tt
     .               ,yss,npointy,zss,xss,tss,iwhich)

          call chk_pos(xss,yss,zss) !this resets logical coords

          write(nunity) real(xss),real(zss)
        endif

        !Check for intersection with Z surface of section
        if((z-zss0)*(zold-zss0) <= 0.0.and.lssz) then
          zss = zss0

          call sursec(xvold,yvold,zvold,tv
     .               ,xold ,yold ,zold ,to
     .               ,x    ,y    ,z    ,tt
     .               ,zss,npointz,xss,yss,tss,iwhich)

          if (rzss) then   !Find R-Z plane positions
            call evalX(xss,yss,zss,xcar,ycar,zcar) !Eval Cartesian coords; also resets logical coords
            rr = sqrt(xcar**2+ycar**2)
            write(nunitrz) real(rr),real(zcar)
          else
            call chk_pos(xss,yss,zss) !this resets logical coords
          endif

          write(nunitz) real(xss),real(yss)
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

      use orbit, ONLY: spt => sing_point
      use spline_field

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,vx,vy,vz

c Begin program

      !Find vector field
      if (spt) then
        !Find logical coords
        call evalXi(y(1),y(2),y(3),x1,y1,z1,ierror)
        if (ierror /= 0) then
          write (*,*) 'Could not find logical coords in ff'
          stop
        endif

        !Get vector field (Cartesian comps)
        call getB(x1  ,y1  ,z1  ,ydot(1),ydot(2),ydot(3)
     .           ,.false.,spt)
      else
        call getB(y(1),y(2),y(3),ydot(1),ydot(2),ydot(3)
     .           ,.false.,spt)
      endif

      !Find jacobian for time step update
      if (.not.spt) then
        call evalBcar(y(1),y(2),y(3),vx,vy,vz)

        ydot(4) = evalJ(y(1),y(2),y(3),ierror)

        ydot(5) = ydot(4)
     .           *sqrt(vx**2+vy**2+vz**2)
        if (ierror /= 0) then
          write (*,*) 'Negative Jacobian in ff'
          stop
        endif
      else
        call evalBcar(x1,y1,z1,vx,vy,vz)

        ydot(4) = 1d0
        ydot(5) = sqrt(vx**2+vy**2+vz**2)
      endif

c End program

      end subroutine ff

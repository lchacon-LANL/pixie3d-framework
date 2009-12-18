
c module orbit
c ######################################################################
      module orbit

        use xdraw_io

        use spline_field

        use grid, ONLY: my_rank,pstop,messg

        !Diagnostics
        integer,private :: nunitd=12345

        logical :: flux = .false.

        !LSODE public variables
        real(8) :: xold,yold,zold,lsode_out(10)
        logical :: sing_point=.false.,sln,exit_loop

        integer :: iiout

        !Private variables
        type,private :: ode_solver
          integer :: niter  !Number of nonlinear iterations
          integer :: order  !Order of integration (1,2)
          real(8) :: atol   !Absolute convergence tolerance
          real(8) :: rtol   !Relative convergence tolerance
          real(8) :: dtime  !Time step
          real(8) :: tmax   !Maximum time
          real(8) :: smax   !Maximum arc length
          logical :: solen  !Whether use solenoidal version or LSODE
          logical :: picard !Whether to use picard or Newton
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

        !B field diagnostics
        integer :: icross,ncross=100000,nbins=100,inx,nnx,nny
        real(8),allocatable,dimension(:,:) :: fpos,qdif
        logical :: bdiag=.false.

        !Private variables
        integer,private :: itime,iterror,ssorder

        integer,private :: nunitx,nunity,nunitz,nunitrz
        integer,private :: npointx,npointy,npointz

        real(8),private :: tvold,told,tprev,time,tmax
     .                    ,svold,sold,sprev,ss  ,smax
        real(8),private :: x,y,z,xcar,ycar,zcar,xss0,yss0,zss0
     .                    ,xvold,yvold,zvold,x00,y00,z00
        real(8),private :: xss,yss,zss,tss

        real(8),private :: xmin,xmax,ymin,ymax,zmin,zmax,rmin

        character(50),private :: ssx,ssy,ssz,ssrz

        logical,private :: cartesian=.false.,cross=.false.
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
      orb_default%ode_conf%tmax   = 0d0    !Maximum orbit pseudotime
      orb_default%ode_conf%smax   = 20000. !Maximum orbit arc length
      orb_default%ode_conf%dtime  = 0.01   !Maximum time step

      orb_default%ode_conf%solen  = .true. !Whether to use solenoidal ODE int.
      orb_default%ode_conf%order  = 1      !Order of integration in solenoidal int. (1,2)
      orb_default%ode_conf%picard = .true. !Whether to use functional iteration
      orb_default%ode_conf%niter  = 20     !Number of nonlinear iterations
      orb_default%ode_conf%rtol   = 0d-8   !Relative tolerance
      orb_default%ode_conf%atol   = 1d-8   !Absolute tolerance
c3-eq      orb_default%ode_conf%lsode_neq = 3   !Number of LSODE equations
      orb_default%ode_conf%lsode_neq = 5   !Number of LSODE equations

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
     .                      ,lsode_neq)

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
        logical,optional :: solen,picard,B_input_is_A

c     Begin program

        call orbit_conf_def

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
          call splineA(bx,by,bz,orb_default%spline_conf%a_gauge
     .                ,input_is_A=B_input_is_A)
        else
          call splineB(bx,by,bz)
        endif

        call splineBcar(bxcar,bycar,bzcar)

        call splineJ(jac )
        call splineX(xmap)

c     Setup diag B diagnostic arrays

        if (bdiag) then
          allocate(fpos(-nbins/2:nbins/2,ncross))
          fpos = 0d0
          allocate(qdif(nnx,ncross))
          qdif = 0d0
        endif

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

        implicit none

c     Call variables

c     Local variables

        integer :: slice
        real(8) :: rr

        real(8),allocatable,dimension(:) ::  mean,sigma,vol,diff

c     Begin program

        call killSplines

c     Finalize B diagnostics

        if (bdiag) then
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

          write (*,*) 'Mean =',mean
          write (*,*) 'Sigma=',sigma
          write (*,*) 'Diff PDF =',diff

          do slice = 1,icross
            diff(slice) = sum(qdif(:,slice))/size(qdif,1)/nny
          enddo

          write (*,*) 'Diff MC =',diff

          close(50)

          open(unit=50,file='b_dif.txt',status='unknown')
          do slice=1,size(qdif,1)
            write(50,*) slice,qdif(slice,1:icross)/nny
          enddo
          close(50)

          deallocate(mean,sigma,vol,diff,fpos,qdif)
        endif

      end subroutine orbit_destroy

c     orbit_find
c     ######################################################################
      subroutine orbit_find(iout,x0,y0,z0,fwd,line_int,lsode_f,ss_order)

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
cc        real(8),optional :: xss,yss,zss,tss

        optional :: line_int

        !Interface for generic line integral routine (optional)
        INTERFACE
          function line_int(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: line_int
          logical :: fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          logical :: ex_loop
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

cc        chk_ss = (.not.integral).and.ss_io  !Do SS only if not integrating and if IO has been set up
        chk_ss = ss_io  !Do SS only if IO has been set up

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

          if (bdiag) icross = 0  !Initialize for B diagnostics

          if(lssx .and. x0 == xss0) then
            write(nunitx) real(y0),real(z0)
          endif

c delta-diag
cc          if (lssx) then
cc            write(nunitd,FMT="(2e14.5)",ADVANCE="NO") 0.5*x0**2,y0  !Initial (Psi,theta) of orbit
cc          endif
c delta-diag

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

            !Initialize B diagnostics
            if (bdiag) then
              icross = icross + 1
              fpos(0,icross) = fpos(0,icross) + 1
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

          !Line-orbit integration
          if (integral) then
            exit_loop = line_int(fwd,xold,x,yold,y,zold,z,sold,ss)
          endif

          if (exit_loop .or. orbit_error(iterror)) then
            if (iiout >=0) write (*,*)
            exit
          endif

          !Check for intersection with surface of section
          if (iterror /= 2 .and. itime > 2 .and. chk_ss) then
            if (smax /= 0d0) then
              call sschk(ssorder,svold,sold,ss
     .                  ,xvold,xold,x,yvold,yold,y
     .                  ,zvold,zold,z,xss,yss,zss,tss)

c delta-diag
cc              if (lssx .and. tss /= 0d0) then
cc                if (fwd) then
cc                  write(nunitd,FMT="(e14.5)",ADVANCE="NO")  tss
cc                else
cc                  write(nunitd,FMT="(e14.5)",ADVANCE="NO") -tss
cc                endif
cc              endif
cccc          write(*,FMT="(e10.3)") tss
c delta-diag

            else
              call sschk(ssorder,tvold,told,time
     .                  ,xvold,xold,x,yvold,yold,y
     .                  ,zvold,zold,z,xss,yss,zss,tss)
            endif
          endif

          !Show progress/terminate
          if (smax /= 0d0) then
            call orbit_progress(ss  ,sprev,smax,fwd,exit_loop)
          elseif (tmax /= 0d0) then
            call orbit_progress(time,tprev,tmax,fwd,exit_loop)
          else
            exit_loop = .false.
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

c delta-diag
cc          if (lssx) then
cc            write(nunitd,*) !Carriage return
cc          endif
c delta-diag

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
        tprev= 0d0

        ss   = 0d0  !Initialize arc-length along line
        sold = 0d0
        sprev= 0d0

        lsode_out = 0d0   !Initialize additional LSODE output

        npointx=0
        npointy=0
        npointz=0

        itime = 0

        exit_loop = .false.

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
      logical :: restart,nt,fwd

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

        call cn_solver(atol,rtol,x,y,z,ierror)

      else                                  !LSODE (non-solenoidal)

        call lsode_step(atol,rtol,x,y,z,ierror)

      endif

c     Output

      if (iiout > 0) then
        write (*,'(a,e10.3,a,e10.3,a,e10.3,a,e10.3,a,1p3e12.3,a,l)')
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

      !LSODE setup
      INTEGER :: LIW, LRW, NEQ
      PARAMETER  (NEQ=10,LIW = 20,LRW = LIW + 16*NEQ)

      INTEGER :: IOPT, IOUT, ISTATE=0, ITASK, ITOL, IWORK(LIW)
     .          ,MF, NEQD(1)
      REAL(8) :: ATOL(NEQ),RTOL(1),RWORK(LRW),T,TOUT,YY(NEQ)

      save ISTATE

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
c3-eq        if (NEQD(1) > 3) YY(4:NEQD(1))=lsode_out(1:NEQD(1)-3)

        TOUT = 0.1*cc*dtime !Output time value (ignored after first call, since ITASK=2)
      else
        YY(1) = x
        YY(2) = y
        YY(3) = z
        YY(4) = time
        YY(5) = ss
        if (NEQD(1) > 5) YY(6:NEQD(1))=lsode_out(1:NEQD(1)-5)
c3-eq        if (NEQD(1) > 3) YY(4:NEQD(1))=lsode_out(1:NEQD(1)-3)

        TOUT  = cc*dtime    !Output time value (ignored after first call, since ITASK=2)
      endif

c     LSODE configuration and call

c3-eq      T    = time   !Independent variable
      T    = 0d0   !Independent variable

      ITOL = 2      !Whether atol is array or not (2 => true)

      if (ISTATE /= 3) then
        RTOL = rrtol            !Relative tolerance
        ATOL = aatol            !Absolute tolerance
      endif

      ITASK = 2     !Mode of operation: take one step and return

      !Whether it is first call (ISTATE=1) or subsequent call (ISTATE=2)
      if ((time == 0d0).or. (sln .and. (.not.cartesian))) ISTATE = 1

      !Optional inputs
cc      if (chk_ss) then
cc        IOPT = 1                !Optional inputs are needed
cc        RWORK(6) = minval( (/(xmax-xmin)/30.
cc     .                      ,(ymax-ymin)/30.
cc     .                      ,(zmax-zmin)/30. /) !Maximum time step allowed (optional input)
cc      else
        IOPT = 0                !Optional inputs not needed
cc      endif

      MF = 10       !Specifies non-stiff problem: implicit Adams-Moulton method

      if (external_f) then
        CALL LSODE (lsode_f,NEQD,YY,T,TOUT,ITOL,RTOL,ATOL,ITASK,
     .              ISTATE,IOPT,RWORK,LRW,IWORK,LIW,lsode_f, MF)
      else
        CALL LSODE (ff     ,NEQD,YY,T,TOUT,ITOL,RTOL,ATOL,ITASK,
     .              ISTATE,IOPT,RWORK,LRW,IWORK,LIW,ff      ,MF)
      endif

      select case(ISTATE)
      case(2)       !Uneventful
        continue
      case(-2)      !Rescale tolerances
        ISTATE = 3
        RTOL=RTOL*RWORK(14)
        ATOL=ATOL*RWORK(14)
      case default  !Exit LSODE
        write (*,*) 'Problems in LSODE: ISTATE=', ISTATE
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

c3-eq      dt   = (abs(T)-time)   !For output
c3-eq      ds   = (abs(T)-time)   !For output
c3-eq
c3-eq      time = abs(T)
c3-eq      ss   = abs(T)
c3-eq
c3-eq      if (NEQD(1) > 3) lsode_out(1:NEQD(1)-3)=cc*YY(4:NEQD(1))

c     End program

      end subroutine lsode_step

c     cn_solver
c     #########################################################################
      subroutine cn_solver(atol,rtol,x,y,z,ierror)

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

          call simp_step(nt,niter,atol,rtol,cc*dt,x,y,z,0,iterr,error)

          if (iterr == 0) iterr = -2  !Notify of Cartesian step

        !Second order splitting
        elseif (order==2) then

          half = 0.5*dt

c         Begin first half-step split

          call simp_step(nt,niter,atol,rtol,cc*half,x,y,z,1,iterr,error)

c         Begin whole step split

          call simp_step(nt,niter,atol,rtol,cc*dt  ,x,y,z,2,iterr,error)

c         Begin second half-step split

          call simp_step(nt,niter,atol,rtol,cc*half,x,y,z,1,iterr,error)

        !First-order integration (30% cheaper!)
        else

          call simp_step(nt,niter,atol,rtol,cc*dt,x,y,z,1,iterr,error)

          call simp_step(nt,niter,atol,rtol,cc*dt,x,y,z,2,iterr,error)

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
              write (*,*) 'Step',iterr
     .                ,' in CN did not converge to prescribed accuracy'
            elseif (iiout == 2) then
              write(*,*) 'Initial error=',error(1)
              write(*,*) 'Final error  =',error(niter),'> tol=',atol
            elseif (iiout > 2) then
              write (*,*) 'Error history=',error
            endif
          elseif (iterr == 4) then
            call pstop('cn_solver','Negative Jacobian')
          else
            write (messg,*) 'Unknown error ',iterr
            call pstop('orbit_step',messg)
          endif

          !Return to initial conditions and try again with smaller time step
          x = xold
          y = yold
          z = zold

          dt = 0.8*dt

          if (dt < 1d-8) call pstop('cn_solver','Time step too small')

          if (iiout > 1) write(*,*)'Trying smaller time step:',dt,'...'

          restart = .true.

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

        call updt_measures(xh,yh,zh,jach*dt)

      end subroutine cn_solver

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

c     simp_step
c     ####################################################################
      subroutine simp_step(newton,niter,atol,rtol,dt,x,y,z,flag,iterr
     .                    ,error)

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

      if (iterr /= 0 .or. (istep /= 0 .and. flag /= istep)) return

      x0=x
      y0=y
      z0=z

      if (iiout > 1) write (*,*) 'simp_step ',flag,' -- X0=',x0,y0,z0

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

           if (iiout > 2) write (*,*) 'simp_step ',flag,' -- res=',res

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

         !Update solution
         x = x + dx(1,1)
         y = y + dx(2,1)
         z = z + dx(3,1)

         if (iiout > 2) write (*,*) 'simp_step ',flag,' -- X =',x,y,z

         !Check convergence
         error(iter) = sqrt(sum(dx**2))

         if (iiout > 2) write (*,*) 'simp_step ',flag,' -- error='
     .                               ,error(iter)

         if (error(iter) < tol) exit

         !Form new residual
         res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z)

      enddo

      if (iiout == 2) then
        write (*,*) 'simp_step ',flag,' -- X =',x,y,z
     .             ,'; error history=',error(1:min(iter,niter))
      endif

      if (iter > niter) iterr = flag

c     End program

      end subroutine simp_step

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

c3-eqc     lsode_orbit_f
c3-eqc     #####################################################################
c3-eq      subroutine lsode_orbit_f(neq,t,y,ydot,xx1,yy1,zz1)
c3-eq
c3-eqc     ---------------------------------------------------------------------
c3-eqc     Driver routine for LSODE integration: ydot = lsode_orbit_f(y,t)
c3-eqc     ---------------------------------------------------------------------
c3-eq
c3-eq      implicit none
c3-eq
c3-eqc     Call variables
c3-eq
c3-eq      integer :: neq
c3-eq      real(8) :: t,y(neq),ydot(neq)
c3-eq      real(8),optional :: xx1,yy1,zz1
c3-eq
c3-eqc     Local variables
c3-eq
c3-eq      integer :: ierror
c3-eq      real(8) :: x1,y1,z1,vx,vy,vz,xc0,yc0,zc0,coeff
c3-eq
c3-eqc     Begin program
c3-eq
c3-eq      if (neq /= 3) then
c3-eq        write (*,*) 'Wrong number of Eqs. in lsode_orbit_f:',neq
c3-eq        stop
c3-eq      endif
c3-eq
c3-eq      !Find logical coordinates
c3-eq      if (sing_point) then
c3-eq        x1 = xold ; y1 = yold ; z1 = zold
c3-eq        call evalXi(iiout,y(1),y(2),y(3),x1,y1,z1,ierror)
c3-eq        if (ierror /= 0) then
c3-eq          write (*,*) ierror
c3-eq          write (*,*) 'Could not find logical coords in lsode_orbit_f'
c3-eq          stop
c3-eq        endif
c3-eq      else
c3-eq        x1 = y(1)
c3-eq        y1 = y(2)
c3-eq        z1 = y(3)
c3-eq      endif
c3-eq
c3-eq      !Get vector field (Cartesian if sing_point is .true.)
c3-eq      call getB(x1,y1,z1,vx,vy,vz,sln,sing_point)
c3-eq
c3-eq      ydot(1) = vx
c3-eq      ydot(2) = vy
c3-eq      ydot(3) = vz
c3-eq
c3-eq      !Time and arc-length ODEs
c3-eq      if (.not.sing_point) then
c3-eq        coeff = evalJ(x1,y1,z1,ierror)
c3-eq        if (ierror /= 0) then
c3-eq          write (*,*) 'Negative Jacobian in lsode_orbit_f'
c3-eq          stop
c3-eq        endif
c3-eq
c3-eq        call evalBcar(x1,y1,z1,vx,vy,vz)
c3-eq        coeff = coeff*sqrt(vx**2+vy**2+vz**2)
c3-eq      else
c3-eq        coeff = sqrt(vx**2+vy**2+vz**2)
c3-eq      endif
c3-eq
c3-eq      ydot = ydot/coeff
c3-eq
c3-eq      if (PRESENT(xx1)) xx1 = x1
c3-eq      if (PRESENT(yy1)) yy1 = y1
c3-eq      if (PRESENT(zz1)) zz1 = z1
c3-eq
c3-eqc     End program
c3-eq
c3-eq      end subroutine lsode_orbit_f

c     lsode_orbit_f
c     #####################################################################
      subroutine lsode_orbit_f(neq,t,y,ydot,xx1,yy1,zz1)

c     ---------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = lsode_orbit_f(y,t)
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)
      real(8),optional :: xx1,yy1,zz1

c     Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,vx,vy,vz,xc0,yc0,zc0

c     Begin program

      if (neq /= 5) then
        write (*,*) 'Wrong number of Eqs. for lsode_orbit',neq
        stop
      endif

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
          write (*,*) 'Negative Jacobian in lsode_orbit_f'
          stop
        endif

        call evalBcar(x1,y1,z1,vx,vy,vz)
        ydot(5) = ydot(4)*sqrt(vx**2+vy**2+vz**2)
      else
        ydot(4) = 1d0
        ydot(5) = sqrt(vx**2+vy**2+vz**2)
      endif

      if (PRESENT(xx1)) xx1 = x1
      if (PRESENT(yy1)) yy1 = y1
      if (PRESENT(zz1)) zz1 = z1

c     End program

      end subroutine lsode_orbit_f

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

cc        if (my_rank == 0 .and. iiout >= 0) write (*,*)

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

cc        if ((.not.lssx) .and. (.not.lssy) .and. (.not.lssz)) then
cc          call pstop('ss_init'
cc     .              ,'No Surface of section selected. Aborting.')
cc        elseif (lssx .and. lssy .and. lssz) then
        if (lssx .and. lssy .and. lssz) then
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

        ssx = 'ssx'//trim(file)
        ssy = 'ssy'//trim(file)
        ssz = 'ssz'//trim(file)
        ssrz= 'ssrz'//trim(file)

        if (lssx) then
          drawin = 'draw'//ssx(1:len_trim(ssx)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssx),'X surface of section'
     .                       ,'y',(/1/),(/'z'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunitx=27
          open(unit=nunitx,file=trim(ssx),form='unformatted')

c delta-diag
cc          nunitd=1234
cc          open(unit=nunitd,file='delta-test.txt')
c delta-diag
        endif

        if (lssy) then
          drawin = 'draw'//ssy(1:len_trim(ssy)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssy),'Y surface of section'
     .                       ,'x',(/1/),(/'z'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunity=28
          open(unit=nunity,file=trim(ssy),form='unformatted')
        endif

        if (lssz) then
          drawin = 'draw'//ssz(1:len_trim(ssz)-4)//'.in'
          call createDrawInGfile(1,1,trim(ssz),'Z surface of section'
     .                       ,'x',(/1/),(/'y'/),(/0/),(/0,0/)
     .                       ,drawin,.false.,.false.)
          nunitz=29
          open(unit=nunitz,file=trim(ssz),form='unformatted')
        endif

        if (rzss) then
          drawin = 'draw'//ssrz(1:len_trim(ssrz)-4)//'.in'
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
        integer :: iwhich,bin

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
        endif

c     Check for intersection with X surface of section

        if (lssx
     .      .and.chk_ss_pos(xss0,xold,x,bcond(1),xmin,xmax,xss)) then

          call sursec(yvold,zvold,xvold,tv
     .               ,yold ,zold ,xold ,to
     .               ,y    ,z    ,x    ,tt
     .               ,xss,npointx,yss,zss,tss,iwhich)

          call chk_pos(xss,yss,zss) !this resets logical coords
          write(nunitx) real(yss),real(zss)
cc          write(*,*) real(yss),real(zss)

        endif

c     Check for intersection with Y surface of section

        if (lssy
     .      .and.chk_ss_pos(yss0,yold,y,bcond(3),ymin,ymax,yss)) then

          call sursec(zvold,xvold,yvold,tv
     .               ,zold ,xold ,yold ,to
     .               ,z    ,x    ,y    ,tt
     .               ,yss,npointy,zss,xss,tss,iwhich)

          call chk_pos(xss,yss,zss) !this resets logical coords

          if (flux) then
            write(nunity) real(0.5*xss**2),real(zss)
          else
            write(nunity) real(xss),real(zss)
          endif
        endif

c     Check for intersection with Z surface of section

        if (lssz
     .      .and.chk_ss_pos(zss0,zold,z,bcond(5),zmin,zmax,zss)) then

          call sursec(xvold,yvold,zvold,tv
     .               ,xold ,yold ,zold ,to
     .               ,x    ,y    ,z    ,tt
     .               ,zss,npointz,xss,yss,tss,iwhich)

          if (rzss) then   !Find R-Z plane positions
            call evalX(xss,yss,zss,xcar,ycar,zcar) !Eval Cartesian coords; also resets logical coords
            rr = sqrt(xcar**2+ycar**2)
            write(nunitrz) real(rr),real(zcar)
cc            write(*,*) 'diag',z,zold,real(rr),real(zcar)   !diag
          else
            call chk_pos(xss,yss,zss) !this resets logical coords
          endif

          if (flux) then
            write(nunitz) real(0.5*xss**2),real(yss) !In flux
          else
            write(nunitz) real(xss),real(yss)        !In radius
          endif

          !Accumulate B diagnostics
          if (bdiag) then
            icross = icross + 1
            if (icross > ncross) then
              write (*,*) 'Exceeded max number of Poincare crossings'
            else
              write (*,*)
              write (*,*) 'B diag: Poincare crossing ',icross

              !PDF diagnostic
              bin = sign(1d0,xss-x00)
     .             *floor(0.5*abs(x00**2-xss**2)*nbins+0.5)
              fpos(bin,icross) = fpos(bin,icross)+1

              !Add quasilinear diffusion diagnostic
              qdif(inx,icross) = qdif(inx,icross)
     .                         + (0.5*(xss**2-x00**2))**2/tss

            endif
          endif

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

      integer :: pp
      real(8) :: x11,x22,x00

c     Begin program

      if (bc == PER) then

c       Shift so that both x1, x2 are positive

        pp =        floor((x1-xmin)/(xmax-xmin))
        pp = min(pp,floor((x2-xmin)/(xmax-xmin)))

        x11 = x1 - pp*(xmax-xmin)
        x22 = x2 - pp*(xmax-xmin)

c       Redefine SS position

        x00 = x0

        if (x0 == xmin) x00 = xmax

        !If x0 is close to domain limit, switch to most convenient one
        if (   (pp <  0 .and.(x11-x0 >xmax-xmin))
     .     .or.(pp >= 0 .and.(x22-x0 >xmax-xmin))) x00 =x0+(xmax-xmin)

c       Check for SS

        found = ((x11-x00)*(x22-x00) <= 0d0)

c       Shift SS position for sursec routine

        if (found) xss = x00+pp*(xmax-xmin)

      else

c       Check for SS

        found = ((x1-x0)*(x2-x0) <= 0d0)

cc        write (*,*) x0,x1,x2
        if (found) then
          xss = x0
        endif

      endif

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
c     Driver routine for LSODE integration: ydot = ff(y,t)
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

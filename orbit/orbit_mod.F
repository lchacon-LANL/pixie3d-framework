
c module orbit
c ######################################################################
      module orbit

        use xdraw_io

        use spline_field

        use grid, ONLY: np,my_rank,pstop,messg,UNDEF_DP

        use file_io
     
        use lyapn
   
        !Diagnostics
        integer,private :: nunitd(4)
        logical,private :: flux = .false.

        !B field diagnostics private variables
        logical,private :: bdiag=.false.
        integer,private :: inx,nnx,nny
        integer,private :: icross,ncross=100000,nbins=100
        real(8),private,allocatable,dimension(:,:) :: fpos,qdif

        !Integrator public variables
        logical :: sln
!$OMP THREADPRIVATE(sln)

        integer :: iiout
!$OMP THREADPRIVATE(iiout)

        !ORBIT private variables
        type :: ode_solver
          integer :: niter     !Number of nonlinear iterations
          integer :: order     !Order of integration (1,2)
          real(8) :: atol      !Absolute convergence tolerance
          real(8) :: rtol      !Relative convergence tolerance
          real(8) :: dtime     !Time step
          real(8) :: tmax      !Maximum pseudotime
          real(8) :: smax      !Maximum arc length
          logical :: solen     !Whether use solenoidal version or LSODE
          logical :: picard    !Whether to use picard or Newton
          logical :: ds_is_Bds !Whether we use B.ds instead of ds as independent variable
          logical :: lsode_dt_lim  !Whether to limit lsode time step
          logical :: lsode_s_indep !Whether to arc-length is indep. var. in LSODE
          logical :: lsode_to_smax !Whether to integrate LSODE exactly to target arc-length
          integer :: lsode_neq !Number of equations for LSODE
        end type ode_solver

        type :: splines
          integer :: order  !Order of splines
          integer :: a_gauge!Gauge for vector potential (1,2)
          logical :: input_is_A !Whether input is vector potential instead of magnetic field
        end type splines

        type :: orbit_def
          type(ode_solver) :: ode_conf
          type(splines)    :: spline_conf
        end type orbit_def

        type(orbit_def),target  :: orbit_conf
        type(ode_solver),pointer :: ode_conf

        integer,private :: ssorder=0,obcnd(6),HIST_LIM=20!,iterror,itime

        real(8),private :: tmax,smax
!$OMP THREADPRIVATE(tmax,smax)
        real(8),private :: xss0,yss0,zss0
        real(8),private :: rmin

        real(8),private :: xmin,xmax,ymin,ymax,zmin,zmax

	!Global variables >>> THREADPRIVATE
        integer, dimension(:), allocatable :: nunitx,nunity,
     .                                        nunitz,nunitrz
!!!!$OMP THREADPRIVATE(nunitx,nunity,nunitz,nunitrz)
	real(8) :: xold,yold,zold,xvold,yvold,zvold,svold,sold
!$OMP THREADPRIVATE(xold,yold,zold,xvold,yvold,zvold,svold,sold)
        real(8),private :: tvold,told,time,ss
!$OMP THREADPRIVATE(tvold,told,time,ss)
        real(8),private,dimension(:),allocatable :: xhist,yhist,zhist
     .                                             ,shist,thist
!$OMP THREADPRIVATE(xhist,yhist,zhist,shist,thist)

        logical :: sing_point=.false.,exit_loop=.false.,orb_init
!$OMP THREADPRIVATE(sing_point,exit_loop,orb_init)
	integer,private :: LSODE_IERR=ORB_OK
!$OMP THREADPRIVATE(LSODE_IERR)
        logical,private :: cartesian_flag = .false.,cross = .false.,
     .             cartesian = .false.
!$OMP THREADPRIVATE(cartesian_flag,cross,cartesian)
        real(8),private :: x,y,z,xcar,ycar,zcar
!$OMP THREADPRIVATE(x,y,z,xcar,ycar,zcar)
        real(8),private :: x00,y00,z00,xss,yss,zss,tss
!$OMP THREADPRIVATE(x00,y00,z00,xss,yss,zss,tss)
        integer,private :: npointx,npointy,npointz
!$OMP THREADPRIVATE(npointx,npointy,npointz)
!#!Marco, number of intersection counter (to deal with ridge=.true.)
        integer,private :: mx,my,mz
!$OMP THREADPRIVATE(mx,my,mz)

!LSODE private variables

        INTEGER :: LSODE_ACC_DIM
        PARAMETER (LSODE_ACC_DIM=30)
        INTEGER,PRIVATE :: LIW, LRW, NEQ
        PARAMETER  (NEQ=20,LIW = 20,LRW = LIW + 16*NEQ)

        INTEGER,PRIVATE :: IOPT,IOUT,ISTATE=0,ITASK,ITOL
     .                  ,IWORK(LIW),MF,NEQD(1)
!$OMP THREADPRIVATE(IOPT,IOUT,ISTATE,ITASK,ITOL,IWORK,MF,NEQD)
        REAL(8),PRIVATE :: LS_ATOL(1),LS_RTOL(1),RWORK(LRW)=0d0
     .                    ,T,TOUT,YY(NEQ),LSODE_ACC(LSODE_ACC_DIM)
!$OMP THREADPRIVATE(LS_ATOL,LS_RTOL,RWORK,T,TOUT,YY,LSODE_ACC)
!     Store of LSODE common block variables
        real(8),dimension(218),private :: RSAV
!$OMP THREADPRIVATE(RSAV)
        integer,dimension(41),private :: ISAV
!$OMP THREADPRIVATE(ISAV)

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

      icross = icross + 1

      if (icross > ncross) then
        write (*,*) 'Exceeded max number of Poincare crossings'
      else
        write (*,*)
        write (*,*) 'B diag: Poincare crossing ',icross

        !PDF diagnostic
        bin = sign(1d0,xs-x0)
     .       *floor(0.5d0*abs(x0**2-xs**2)*nbins+0.5)
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

        open(unit=50,file='b_stat.txt',status='unknown')

        allocate(mean(icross),sigma(icross),vol(icross),diff(icross))

        mean  = 0d0
        sigma = 0d0
        vol   = 0d0

        do slice=-nbins/2,nbins/2
          rr = (1d0*slice)/(nbins+1)
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
      orbit_conf%ode_conf%tmax   = 0d0    !Maximum orbit pseudotime
      orbit_conf%ode_conf%smax   = 0d0    !Maximum orbit arc length
      orbit_conf%ode_conf%dtime  = 1d-2   !Maximum time step

      orbit_conf%ode_conf%solen  = .true. !Whether to use solenoidal ODE int.
      orbit_conf%ode_conf%order  = 1      !Order of integration in solenoidal int. (1,2)
      orbit_conf%ode_conf%picard = .true. !Whether to use functional iteration
      orbit_conf%ode_conf%niter  = 20     !Number of nonlinear iterations
      orbit_conf%ode_conf%rtol   = 0d-8   !Relative tolerance
      orbit_conf%ode_conf%atol   = 1d-8   !Absolute tolerance

      orbit_conf%ode_conf%ds_is_Bds = .false.  !Whether we use B.ds instead of ds as independent variable

      !LSODE default setup
      orbit_conf%ode_conf%lsode_dt_lim = .false. !Do not limit LSODE time step
      orbit_conf%ode_conf%lsode_s_indep  = .false. !Find arc-length as dependent variable
      if (orbit_conf%ode_conf%lsode_s_indep) then
        orbit_conf%ode_conf%lsode_neq = 3 !Number of LSODE equations
      else
        orbit_conf%ode_conf%lsode_neq = 5 !Number of LSODE equations
      endif
      orbit_conf%ode_conf%lsode_to_smax  = .false. !Return after every single LSODE step

      !Splines
      orbit_conf%spline_conf%order = 2    !Spline interpolation order
      orbit_conf%spline_conf%a_gauge = 2  !Gauge for vector potential

      end subroutine orbit_conf_def

c     orbit_params
c     ######################################################################
      subroutine orbit_params(tmax,smax,dtime
     .                       ,niter,atol,rtol,order,solen,picard
     .                       ,spline_order,gauge,B_input_is_A,r_min
     .                       ,lsode_neq,lsode_dt_lim,lsode_to_smax
     .                       ,lsode_s_indep,ds_is_Bds)

c     ----------------------------------------------------------------------
c     Sets up parameters for orbit integrator.
c     ----------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),optional :: tmax,smax,dtime,atol,rtol,r_min
        integer,optional :: niter,order,spline_order,gauge,lsode_neq
        logical,optional :: solen,picard,B_input_is_A,lsode_dt_lim
     .                     ,lsode_to_smax,lsode_s_indep,ds_is_Bds

c     Local variables

c     Begin program

        call orbit_conf_def   !Defaults are defined here.

        if (PRESENT(tmax)) then
          orbit_conf%ode_conf%tmax = tmax
        endif

        if (PRESENT(smax)) then
          orbit_conf%ode_conf%smax = smax
        endif

        if (PRESENT(dtime)) then
          orbit_conf%ode_conf%dtime = dtime
        endif

        if (PRESENT(niter)) then
          orbit_conf%ode_conf%niter = niter
        endif

        if (PRESENT(atol)) then
          orbit_conf%ode_conf%atol = atol
        endif

        if (PRESENT(rtol)) then
          orbit_conf%ode_conf%rtol = rtol
        endif

        if (PRESENT(order)) then
          orbit_conf%ode_conf%order = order
        endif

        if (PRESENT(solen)) then
          orbit_conf%ode_conf%solen = solen
        endif

        if (PRESENT(picard)) then
          orbit_conf%ode_conf%picard = picard
        endif

        if (PRESENT(lsode_neq)) then
          orbit_conf%ode_conf%lsode_neq = lsode_neq
        endif

        if (PRESENT(lsode_dt_lim)) then
          orbit_conf%ode_conf%lsode_dt_lim = lsode_dt_lim
        endif

        if (PRESENT(lsode_s_indep)) then
          orbit_conf%ode_conf%lsode_s_indep = lsode_s_indep
        endif

        if (PRESENT(lsode_to_smax)) then
          orbit_conf%ode_conf%lsode_to_smax = lsode_to_smax
        endif

        if (PRESENT(ds_is_Bds)) then
          orbit_conf%ode_conf%ds_is_Bds = ds_is_Bds
        endif

        if (PRESENT(spline_order)) then
          orbit_conf%spline_conf%order = spline_order
        endif

        if (PRESENT(gauge)) then
          orbit_conf%spline_conf%a_gauge = gauge
        endif

        if (PRESENT(B_input_is_A)) then
          orbit_conf%spline_conf%input_is_A = B_input_is_A
        endif

        if (PRESENT(r_min)) then
          rmin = r_min
        else
          rmin = 1d-3
        endif

        ode_conf => orbit_conf%ode_conf

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

      end subroutine orbit_params

c     orbit_setup
c     ######################################################################
      subroutine orbit_setup(nx,ny,nz,xx,yy,zz
     .                      ,bx,by,bz,bxcar,bycar,bzcar,jac,xmap,bcnd
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

        integer :: nx,ny,nz,bcnd(6)
        real(8) :: xx(nx),yy(ny),zz(nz)
        real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
        real(8) :: bxcar(nx,ny,nz),bycar(nx,ny,nz),bzcar(nx,ny,nz)
        real(8) :: xmap(nx,ny,nz,3),jac(nx,ny,nz)

        !Optional parameters
        real(8),optional :: tmax,smax,dtime,atol,rtol,r_min
        integer,optional :: niter,order,spline_order,gauge,lsode_neq
        logical,optional :: solen,picard,B_input_is_A,lsode_dt_lim
     .                     ,lsode_to_smax,lsode_s_indep,ds_is_Bds

c     Local variables

        integer :: nn

c     Begin program

        call orbit_params(
     .        tmax          =tmax
     .       ,smax          =smax         
     .       ,dtime         =dtime        
     .       ,niter         =niter        
     .       ,atol          =atol         
     .       ,rtol          =rtol         
     .       ,order         =order        
     .       ,solen         =solen        
     .       ,picard        =picard       
     .       ,spline_order  =spline_order 
     .       ,gauge         =gauge        
     .       ,B_input_is_A  =B_input_is_A 
     .       ,r_min         =r_min        
     .       ,lsode_neq     =lsode_neq    
     .       ,lsode_dt_lim  =lsode_dt_lim 
     .       ,lsode_to_smax =lsode_to_smax
     .       ,lsode_s_indep =lsode_s_indep
     .       ,ds_is_Bds     =ds_is_Bds    
     .       )

        obcnd = bcnd

c     Do splines

        call setupSplines(nx,ny,nz,xx,yy,zz
     .                   ,orbit_conf%spline_conf%order,obcnd,solen)

        call get_sp_domain_limits(xmin,xmax,ymin,ymax,zmin,zmax)

        if (    ode_conf%solen
     .      .or.orbit_conf%spline_conf%input_is_A) then
          call splineA(bx,by,bz,  orbit_conf%spline_conf%a_gauge
     .                ,input_is_A=orbit_conf%spline_conf%input_is_A)
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

c     Close SS output
!       if (ss_io) call ss_io_end
       if (ss_io .and. (.not. ridge) 
     .     .and. (.not. auto_corr)) call ss_io_end

c     Kill spline info

        call killSplines

c     Finalize B diagnostics

        if (bdiag) call B_diag_finalize

      end subroutine orbit_destroy

c     orbit_find
c     ######################################################################
      subroutine orbit_find(iout,x0,y0,z0,fwd,line_int,lsode_f
     .                     ,t_max,s_max,ierror,lsode_sum)

c     ----------------------------------------------------------------------
c     Integrates orbit up to specified limit in arc length or pseudotime
c     (specified in module). On input:
c       * iout: level of diagnostic output
c       * x0,y0,z0: starting point of orbit
c       * fwd: whether to go forward or backward
c     Optional input:
c       * line_int: external routine to perform line integrals
c       * lsode_f: external routine to perform line integrals using LSODE
c       * ss_check: check for surface of section (default: not when
c                   performing integrals along line)
ccc       * ss_order: order of interpolation for SS (default: 1)
c       * t_max,s_max: override maximum orbit lengths, measured in pseudotime
c                      units (t_max, for solenoidal int), or arclength (s_max)
c     Optional output:
c       * ierror: gives error code (see spline_mod.f for codes)
c       * lsode_sum: gives LSODE accumulators (real)
c     ----------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: iout
        real(8),intent(IN) :: x0,y0,z0
        logical :: fwd

        real(8),optional,INTENT(IN) :: t_max,s_max

        integer,optional,INTENT(OUT) :: ierror
        real(8),optional,INTENT(OUT) :: lsode_sum(LSODE_ACC_DIM)

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

	integer :: iterror,itime,np
        real(8) :: temp_s,xh,yh,zh,sh
        logical :: chk_ss,integral,non_fatal_err
!#!     Marco, to deal with connection length in the backward case
        logical :: cnc_outside,cnc_cross

        integer,save :: thrtot,thrnum
!$OMP THREADPRIVATE(thrtot,thrnum)

c     Begin program

        call get_omp_thread_id(thrnum,thrtot)

        iiout = iout

        iterror = 0

        integral = PRESENT(line_int)
        
!#!Marco
        if (auto_corr) ss_io = .true.
        chk_ss = (.not.integral).and.ss_io  !Do SS only if not integrating and if IO has been set up

        !Override maximum limits
        if (PRESENT(t_max)) then
          tmax = t_max
          ode_conf%tmax = tmax
        else
          tmax = ode_conf%tmax 
        endif

        if (PRESENT(s_max)) then
          smax = s_max
          ode_conf%smax = smax
        else
          smax = ode_conf%smax
        endif

        if (PRESENT(lsode_f).and.(.not.PRESENT(lsode_sum))
     .     .and.my_rank==0) then
          write(*,*)'WARNING: no LSODE accumulators given to orbit_find'
        endif
          
c     Check ranges

        if((x0.le.xmin.or.x0.ge.xmax)
     .     .and.(obcnd(1)/=PER))then
          write(*,*) 'x0(init)=',x0,' out of range:',xmin,xmax
          return
        endif

        if((y0.le.ymin.or.y0.ge.ymax)
     .     .and.(obcnd(3)/=PER))then
          write(*,*) 'y0(init)=',y0,' out of range:',ymin,ymax
          return
        endif

        if((z0.le.zmin.or.z0.ge.zmax)
     .     .and.(obcnd(5)/=PER))then
          write(*,*) 'z0(init)=',z0,' out of range:',zmin,zmax
          return
        endif

c     Initialize orbit

        call orbit_init(x0,y0,z0)

c     Daniele,23/03/2009: dump the initial point of the orbit

        if (chk_ss) then

          if(lssx .and. x0 == xss0) then
           if (.not. ridge .and. .not. auto_corr) then
            write(nunitx(thrnum)) real(y0),real(z0)
           endif
          endif

          if(lssy .and. y0 == yss0) then
           if (.not. ridge .and. .not. auto_corr) then
            if (flux) then
              write(nunity(thrnum)) real(0.5*x0**2),real(z0)
            else
              write(nunity(thrnum)) real(x0),real(z0)
            endif
           endif
          endif

          if (lssy .and. ridge) then
!!#! Marco, write initial position in the y=yss0 plane,ridgeridge
           call ridge_write(x0,z0,thrnum)
          endif

          if(lssz .and. z0 == zss0) then
           if (.not. ridge .and. .not. auto_corr) then
            if (flux) then
              write(nunitz(thrnum)) real(0.5*x0**2),real(y0) !In flux
            else
              write(nunitz(thrnum)) real(x0),real(y0)        !In radius
            endif
           endif
          endif

          if (lssz .and. ridge) then
!!#! Marco, write initial position in the z=zss0 plane
!            write(*,*) 'WRITE_ci',x0,y0,thrnum
            call ridge_write(x0,y0,thrnum)
          endif

        endif

!#! lyapunov
        if (lyap) call lyap_io_init

!#! auto_corr
        if (auto_corr) then
          call auto_corr_io_init(deltass, ss_min, x0,y0,z0,ss_check)
          ac_rtrn = .false.
        endif

!     Orbit integral

!#! Marco, initialize intersections counters
        mx = 0 ; my = 0 ; mz = 0

!#! Marco, initialize cnc_outside
        if (connect) cnc_outside = .not.cnc_inside_bdry(x0,y0,z0)

        itime = 0 

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
        
          if (lyap) call lyap_write(x,ss)

          if (auto_corr) then
           call auto_corr_write(x,y,z,atc_tol,deltass
     .                         ,ss_check,ac_rtrn)
           if (ac_rtrn) then
            call auto_corr_io_end(lssx,lssy,lssz,rzss
     .                           ,x0,y0,z0,ss)
             return
           endif
          endif
          
          xhist(3:HIST_LIM) = xhist(2:HIST_LIM-1)
          xhist(2) = x

          yhist(3:HIST_LIM) = yhist(2:HIST_LIM-1)
          yhist(2) = y

          zhist(3:HIST_LIM) = zhist(2:HIST_LIM-1)
          zhist(2) = z

          thist(3:HIST_LIM) = thist(2:HIST_LIM-1)
          thist(2) = time

          shist(3:HIST_LIM) = shist(2:HIST_LIM-1)
          shist(2) = ss

          !Find new LOGICAL positions, new time, new arc length
          call orbit_step(x,y,z,iterror,forward=fwd,lsode_f=lsode_f)

          !#! Marco, january 15th 2019
          if (connect) then
            cnc_cross =
     .          ((.not.cnc_outside).and.(.not.cnc_inside_bdry(x,y,z)))
     .           .or. (cnc_outside .and.      cnc_inside_bdry(x,y,z))
            if (cnc_cross) then
              call connect_write(x0,y0,z0,ss,thrnum)
              !#! Marco, no need to continue after threshold is reached
              exit
            endif
          endif

          exit_loop = (orbit_error(iterror).or.get_orbit_exit())
          non_fatal_err = (iterror<=0)

          xhist(1) = x
          yhist(1) = y
          zhist(1) = z
          thist(1) = time
          shist(1) = ss

          !Check for orbit end and interpolate quantities
          if (smax /= 0d0) then
            exit_loop = exit_loop.or.(abs(ss) >= smax)
            !#!Marco, qui per smax in connect
            if (exit_loop .and. connect) then
              call connect_write(x0,y0,z0,big_n,thrnum)
            endif

            if (abs(ss) > smax) then !4th order interpolation of orbit end points
              ss = sign(1d0,ss)*smax
              call sursec_ho(ss,shist,xhist,yhist,zhist,np,x,y,z,4)
            endif
          elseif (tmax /= 0d0) then
            exit_loop = exit_loop.or.(abs(time) >= tmax)
            !#!Marco, qui per smax in connect
            if (exit_loop .and. connect) then
              call connect_write(x0,y0,z0,big_n,thrnum)
            endif

            if (abs(time) > tmax) then !4th order interpolation of orbit end points
              time = sign(1d0,time)*tmax
              call sursec_ho(time,thist,xhist,yhist,zhist,np,x,y,z,4)
            endif
          endif

          !Line-orbit integration
          if (integral.and.non_fatal_err) then
            exit_loop = exit_loop
     .              .or.(line_int(fwd,xold,x,yold,y,zold,z,sold,ss))
          endif

          !Check for intersection with surface of section
          if (non_fatal_err .and. itime > 2 .and. chk_ss) then
            call sschk_ho(xss,yss,zss,tss)
cc            if (smax /= 0d0) then
cc              call sschk(svold,sold,ss
cc     .                  ,xvold,xold,x,yvold,yold,y
cc     .                  ,zvold,zold,z,xss,yss,zss,tss)
cc            else
cc              call sschk(tvold,told,time
cc     .                  ,xvold,xold,x,yvold,yold,y
cc     .                  ,zvold,zold,z,xss,yss,zss,tss)
cc            endif
          endif

          !Show progress
          if (smax /= 0d0) then
            call orbit_progress(ss  ,smax,fwd,exit_loop)
          elseif (tmax /= 0d0) then
            call orbit_progress(time,tmax,fwd,exit_loop)
          endif

          if (exit_loop) then
            if (iiout >=0) write (*,*)
            exit
          endif

          !stop orbit integration if # of intersections > prescribed number
          if (auto_corr.and.auto_corr_inters_stop(lssx,lssy,lssz)) exit
!          if (auto_corr) then
!           write(*,*) 'DIAGac2',intersect_x,intersect_y,intersect_z
!          endif

!#!Marco, check maximum intersection number
          if (ridge) then

           if (my .gt. max_inter_y) exit

           if (mz .gt. max_inter_z) exit

          endif

        enddo

        if (auto_corr) call auto_corr_io_end(lssx,lssy,lssz,rzss
     .                                      ,x0,y0,z0,ss)

        if (lyap) call lyap_io_end

        if (PRESENT(ierror)) ierror = iterror

        if (PRESENT(lsode_sum)) lsode_sum = LSODE_ACC

c     Write orbit separator in XDRAW file

       if ((.not. ridge) .and. (.not. auto_corr) .and. chk_ss) then
          if (lssx) write(nunitx(thrnum))
          if (lssy) write(nunity(thrnum))
          if (lssz) write(nunitz(thrnum))
          if (rzss) write(nunitrz(thrnum))
       endif

c     End orbit

        call orbit_end

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

        allocate(xhist(HIST_LIM)
     .          ,yhist(HIST_LIM)
     .          ,zhist(HIST_LIM)
     .          ,thist(HIST_LIM)
     .          ,shist(HIST_LIM))

        xhist = UNDEF_DP  !Initialize orbit history accumulators
        yhist = UNDEF_DP
        zhist = UNDEF_DP
        thist = UNDEF_DP
        shist = UNDEF_DP

        LSODE_ACC = 0d0   !Initialize additional LSODE output

        npointx=0   !Poincare points accumulator for X surface of section
        npointy=0   !Poincare points accumulator for Y surface of section
        npointz=0   !Poincare points accumulator for Z surface of section

        if (auto_corr) call auto_corr_counters_init

        call set_orbit_exit(.false.)

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

c     orbit_end
c     ######################################################################
      subroutine orbit_end
c     ----------------------------------------------------------------------
c     Finalizes orbit integrator
c     ----------------------------------------------------------------------

      integer :: alloc_stat

      deallocate(xhist,yhist,zhist,thist,shist,STAT=alloc_stat)

      end subroutine orbit_end

c     set_orbit_exit
c     ######################################################################
      subroutine set_orbit_exit(exitl)

c     ----------------------------------------------------------------------
c     Sets orbit exit control variable
c     ----------------------------------------------------------------------

        implicit none

        logical :: exitl

        exit_loop = exitl

      end subroutine set_orbit_exit

c     get_orbit_exit
c     ######################################################################
      function get_orbit_exit() result(exitl)
c     ----------------------------------------------------------------------
c     Gets orbit exit control variable
c     ----------------------------------------------------------------------

        implicit none

        logical :: exitl

        exitl = exit_loop

      end function get_orbit_exit

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

      !Interface for optional LSODE function evaluation
      INTERFACE
        subroutine lsode_f(neq,t,y,ydot)
        integer :: neq
        real(8) :: t,y(neq),ydot(neq)
        end subroutine lsode_f
      END INTERFACE

      optional :: lsode_f

c     Local Variables

      real(8) :: x0,y0,z0,x1,y1,z1,half,jach,dt=10d0
      real(8) :: vx,vy,vz,error(ode_conf%niter),atol,rtol,dtime
      real(8) :: xh,yh,zh,xnew,ynew,znew,xc,yc,zc,ds,cc
      integer :: iter,order,niter
      logical :: restart,nt,fwd,lsode_3eq,lsode_to_smax

      character(5) :: solver

      logical  :: external_f

      save dt
!$OMP THREADPRIVATE(dt)

c     Begin program

      external_f = PRESENT(lsode_f) !LSODE function is external

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

      ierror = ORB_OK   !Initialize domain bounds error flag

      restart = .false. !Initialize restarting status flag

      if (fwd) then
        cc = 1d0
      else
        cc =-1d0
      endif

c     Check for SP

      if (.not.cartesian) cartesian = (obcnd(1)==SP).and.(x<rmin)

c     Time integration

      if (sln .and. (.not.cartesian)) then  !Solenoidal

        call vp_solver(atol,rtol,x,y,z,ierror)

        solver = 'VP'

      else                                  !LSODE (non-solenoidal)

        sln = orbit_conf%spline_conf%input_is_A !Whether input is A

        if (lsode_3eq) then
          call lsode_step_3eq(atol,rtol,x,y,z,ierror)
        else
          call lsode_step    (atol,rtol,x,y,z,ierror)
        endif

        solver = 'LSODE'

      endif

c     Output

      if (iiout > 0) then
        write (*,'(a,e10.3,a,e10.3,a,e10.3,a,e10.3,a,3es12.3,a,l3)')
     .         solver//'; Time=',time
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

c     Interfaces

      INTERFACE
        subroutine ff(neq,t,y,ydot)
        integer :: neq
        real(8) :: t,y(neq),ydot(neq)
        end subroutine ff
      END INTERFACE

c     Begin program

      if (cartesian .and. (.not.sing_point)) ISTATE=1 !Restart LSODE for cartesian step

      sing_point = cartesian

      NEQD = ode_conf%lsode_neq    !Number of equations

c     Initial conditions

      if (cartesian) then
        !Find Cartesian coordinates
        call evalX(x,y,z,xc,yc,zc,ierr=iterr)
        if (iterr /= ORB_OK) return

        YY(1) = xc
        YY(2) = yc
        YY(3) = zc
        YY(4) = time
        YY(5) = ss
        if (NEQD(1) > 5) YY(6:NEQD(1))=LSODE_ACC(1:NEQD(1)-5)

        TOUT = 0.1*cc*dtime !Output time value (ignored after first call, since ITASK=2)
      else
        YY(1) = x
        YY(2) = y
        YY(3) = z
        YY(4) = time
        YY(5) = ss
        if (NEQD(1) > 5) YY(6:NEQD(1))=LSODE_ACC(1:NEQD(1)-5)

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

!      write (*,*)'Number of equations',NEQD
!      write (*,*)'Input state=',YY(1:NEQD(1))
!      write (*,*)'Current time=',T
!      write (*,*)'Next time=',TOUT
!      write (*,*)'Tolerances=',ITOL,LS_RTOL,LS_ATOL
!      write (*,*)'LSODE task=',ITASK
!      write (*,*)'LSODE state=',ISTATE
!      write (*,*)'Optional inputs required?',IOPT
!      write (*,*)'Real work array size=',size(RWORK),' declared=',LRW
!      write (*,*)'Integer work array size=',size(IWORK),' declared=',LIW
!      write (*,*)'LSODE method=',MF

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
      case(-2)      !Rescale tolerances and continue
        ISTATE = 3
        RTOL=RTOL*RWORK(14)
        ATOL=ATOL*RWORK(14)
        write (*,*) 'Rescaling tolerances and continuing...'
      case default  !Exit LSODE
        write (*,*) 'Error code ISTATE=',ISTATE
        write (*,*) 'See lsode.f for instructions'
        write (*,*) 'Aborting...'
        stop
      end select

      iterr = LSODE_IERR !Recover errors from LSODE residual call
      if (iterr /= ORB_OK) return

c     Postprocessing to find logical coordinates

      if (cartesian) then
        x = xold ; y = yold ; z = zold
        call evalXi(iiout,YY(1),YY(2),YY(3),x,y,z,iterr)

        if (iterr /= 0) then
          iterr = ORB_SP_ERR  !Error around SP
          return
        endif

        if (x >= rmin) then
          iterr = ORB_CAR_ST   !End of Cartesian step
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

      if (NEQD(1) > 5) LSODE_ACC(1:NEQD(1)-5)=cc*YY(6:NEQD(1))

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

c     Interfaces

      INTERFACE
        subroutine ff_3eq(neq,t,y,ydot)
        integer :: neq
        real(8) :: t,y(neq),ydot(neq)
        end subroutine ff_3eq
      END INTERFACE

c     Begin program

      if (cartesian .and. (.not.sing_point)) ISTATE=1 !Restart LSODE for cartesian step

      sing_point = cartesian

      NEQD = ode_conf%lsode_neq    !Number of equations

c     Initial conditions

      if (cartesian) then
        !Find Cartesian coordinates
        call evalX(x,y,z,xc,yc,zc,ierr=iterr)
        if (iterr /= ORB_OK) return

        YY(1) = xc
        YY(2) = yc
        YY(3) = zc
        if (NEQD(1) > 3) YY(4:NEQD(1))=LSODE_ACC(1:NEQD(1)-3)

        TOUT = 0.1*cc*dtime !Output time value (ignored after first call, since ITASK=2)
      else
        YY(1) = x
        YY(2) = y
        YY(3) = z
        if (NEQD(1) > 3) YY(4:NEQD(1))=LSODE_ACC(1:NEQD(1)-3)

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
        IWORK(6) = 100000        !Maximum number of steps
      endif

      MF = 10       !Specifies non-stiff problem: implicit Adams-Moulton method

c$$$      if (external_f) then
c$$$      write (*,*)'Number of equations',NEQD
c$$$      write (*,*)'Input state=',YY(1:NEQD(1))
c$$$      write (*,*)'Current time=',T
c$$$      write (*,*)'Next time=',TOUT
c$$$      write (*,*)'Tolerances=',ITOL,LS_RTOL,LS_ATOL
c$$$      write (*,*)'LSODE task=',ITASK
c$$$      write (*,*)'LSODE state=',ISTATE
c$$$      write (*,*)'Optional inputs required?',IOPT
c$$$      write (*,*)'Real work array size=',size(RWORK),' declared=',LRW
c$$$      write (*,*)'Integer work array size=',size(IWORK),' declared=',LIW
c$$$      write (*,*)'LSODE method=',MF
c$$$      write (*,*)'External function=',external_f
c$$$      write (*,*)'====================================================='
c$$$      endif

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
        write (*,*) 'Continuing with LSODE integral...'
      case(-2)      !Rescale tolerances and call again
        ISTATE = 3
        RTOL=RTOL*RWORK(14)
        ATOL=ATOL*RWORK(14)
        write (*,*) 'Rescaling tolerances and continuing...'
      case default  !Exit LSODE
        write (*,*) 'Error code ISTATE=',ISTATE
        write (*,*) 'See lsode.f for instructions'
        write (*,*) 'Aborting...'
        stop
      end select

      iterr = LSODE_IERR  !Set in LSODE integration routine lsode_f or ff_3eq
      if (iterr /= ORB_OK) return

c     Postprocessing to find logical coordinates

      if (cartesian) then
        x = xold ; y = yold ; z = zold
        call evalXi(iiout,YY(1),YY(2),YY(3),x,y,z,iterr)

        if (iterr /= 0) then
          iterr = ORB_SP_ERR !Error around SP
          return
        endif

        if (x >= rmin) then
          iterr = ORB_CAR_ST   !End of Cartesian integration
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

      if (NEQD(1) > 3) LSODE_ACC(1:NEQD(1)-3)=cc*YY(4:NEQD(1))

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

 1      iterr = ORB_OK

        !Cartesian step (near SP)
        if (cartesian) then

          call vp_step(nt,niter,atol,rtol,cc*dt,x,y,z,0,iterr,error)

          if (iterr == ORB_OK) iterr = ORB_CAR_ST  !Notify of Cartesian step

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

c       Find center-point

        xh = 0.5*(xold+x)
        yh = 0.5*(yold+y)
        zh = 0.5*(zold+z)

c       Find jacobian at mid-point (for time/arc-length integral)

        if (cartesian) then
          jach = 1d0
        else
          jach = evalJ(xh,yh,zh,ierror)
          if (ierror /= ORB_OK) return
        endif

c       Check for errors

        if (iterr /= ORB_OK .and. iterr /= ORB_CAR_ST) then

          if (iterr == 1 .or. iterr == 2 .or. iterr == 3) then  !Convergence failure in VP stage "iterr"

            if (iiout >= 2) then
              write (*,*) 'VP -- Step',iterr
     .                 ,' in CN did not converge to prescribed accuracy'
            endif

            if (iiout == 2) then
              write(*,*)'VP --Initial error=',error(1)
              write(*,*)'VP --Final error  =',error(niter),'> tol=',atol
            elseif (iiout > 2) then
              write(*,*)'VP --Error history=',error
            endif

            !Return to initial conditions and try again with smaller time step
            x = xold
            y = yold
            z = zold

            dt = 0.8*dt

            if (dt < 1d-8) then
              ierror = ORB_SML_DT  !Time step too small
              return
            else
              if (iiout > 1)
     .           write(*,*)'VP -- Trying smaller time step:',dt,'...'

              restart = .true.

              go to 1
            endif

          else

            ierror = iterr !Return with error
            return

          endif

        elseif (restart) then

          restart = .false.

          ierror = ORB_ADJ_DT  !Indicates time step has been decreased

          if (iiout > 1) then
            write (*,*) 'VP -- CN succeeded! Proceeding.'
            write (*,*)
          endif

        endif

c       Advance time and arc length

        call updt_measures(xh,yh,zh,jach*dt,iterr)
        if (iterr /= ORB_OK) ierror = iterr
       
      end subroutine vp_solver

c     updt_measures
c     #########################################################################
      subroutine updt_measures(x,y,z,dt,ierr)

      implicit none

      integer :: ierr
      real(8) :: x,y,z,dt

      real(8) :: bx,by,bz

      !Update time
      time=time+dt  !Include Jacobian factor to advance time

      !Update arc length
      call evalBcar(x,y,z,bx,by,bz,ierr)

      ds = dt*sqrt(bx*bx+by*by+bz*bz)

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
      if (iterr /= ORB_OK.or.(istep /= 0 .and. flag /= istep)) return

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
      res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z,iterr)
      if (iterr /= ORB_OK) return

      tol = atol + rtol*sqrt(sum(res**2))

      !Iterations for half-step split
      do iter=1,niter

         !Find update
         if (newton) then

           if (iiout > 2) write (*,*) 'vp_step ',flag,' -- res=',res

           xh = 0.5*(x+x0)
           yh = 0.5*(y+y0)
           zh = 0.5*(z+z0)
           JJ = formJacobian(size,flag,dt,xh,yh,zh,iterr)
           if (iterr /= ORB_OK) return

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
         res(:,1) = -formResidual(flag,dt,x0,y0,z0,x,y,z,iterr)
         if (iterr /= ORB_OK) return

      enddo

      if (iiout == 2) then
        write (*,*) 'vp_step ',flag,' -- X =',x,y,z
     .             ,'; error history=',error(1:min(iter,niter))
      endif

      if (iter > niter .and. iterr == ORB_OK) iterr = flag

c     End program

      end subroutine vp_step

c     formResidual
c     ###################################################################
      function formResidual(flag,dt,x0,y0,z0,x1,y1,z1,ierr) result(res)

      implicit none

      integer :: flag,ierr
      real(8) :: dt,x0,y0,z0,x1,y1,z1,res(3)

      real(8) :: xh,yh,zh,vx,vy,vz,xc0,yc0,zc0,x11,y11,z11

      ierr = ORB_OK

      !Find mid point
      xh = 0.5*(x1+x0)
      yh = 0.5*(y1+y0)
      zh = 0.5*(z1+z0)

      !Find vector field at midpoint
      call getB(xh,yh,zh,vx,vy,vz,.true.,cartesian,ierr,flag=flag)
      if (ierr /= ORB_OK) return

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
        call evalXi(iiout,xc0,yc0,zc0,x11,y11,z11,ierr)
        if (ierr /= ORB_OK) ierr = ORB_SP_ERR !Error around SP

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
      function formJacobian(size,flag,dt,x1,y1,z1,ierr) result(JJ)

      implicit none

      integer :: size,flag,ierr
      real(8) :: dt,x1,y1,z1,JJ(size,size)

      call evalCurlA(x1,y1,z1,JJ(1,1),JJ(2,1),JJ(3,1),ierr,flag=flag
     .              ,idx=1,idy=0,idz=0)

      call evalCurlA(x1,y1,z1,JJ(1,2),JJ(2,2),JJ(3,2),ierr,flag=flag
     .              ,idx=0,idy=1,idz=0)

      call evalCurlA(x1,y1,z1,JJ(1,3),JJ(2,3),JJ(3,3),ierr,flag=flag
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

cc      integer :: ierror
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
        call evalXi(iiout,y(1),y(2),y(3),x1,y1,z1,LSODE_IERR)
        if (LSODE_IERR /= ORB_OK) then
          LSODE_IERR = ORB_SP_ERR
          return
        endif
      else
        x1 = y(1)
        y1 = y(2)
        z1 = y(3)
      endif

      !Get vector field (Cartesian if sing_point is .true.)
      call getB(x1,y1,z1,vx,vy,vz,sln,sing_point,LSODE_IERR)
      if (LSODE_IERR /= ORB_OK) return

      ydot(1) = vx
      ydot(2) = vy
      ydot(3) = vz

      !Time and arc-length ODEs
      if (.not.sing_point) then
        jac = evalJ(x1,y1,z1,LSODE_IERR)
        if (LSODE_IERR /= ORB_OK) return

        call evalBcar(x1,y1,z1,vx,vy,vz,LSODE_IERR)
        if (LSODE_IERR /= ORB_OK) return
      else
        jac = 1d0
      endif

      bnorm = vx*vx+vy*vy+vz*vz
      if (.not.bwgt) bnorm = sqrt(bnorm)

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
        call evalXi(iiout,y(1),y(2),y(3),x1,y1,z1,LSODE_IERR)
        if (LSODE_IERR /= ORB_OK) then
          LSODE_IERR = ORB_SP_ERR
          return
        endif
      else
        x1 = y(1)
        y1 = y(2)
        z1 = y(3)
      endif

      !Get vector field (Cartesian if sing_point is .true.)
      call getB(x1,y1,z1,vx,vy,vz,sln,sing_point,LSODE_IERR)
      if (LSODE_IERR /= ORB_OK) return

      ydot(1) = vx
      ydot(2) = vy
      ydot(3) = vz

      !Time and arc-length ODEs
      if (.not.sing_point) then
        ydot(4) = evalJ(x1,y1,z1,LSODE_IERR)
        if (LSODE_IERR /= ORB_OK) return

        call evalBcar(x1,y1,z1,vx,vy,vz,LSODE_IERR)
        if (LSODE_IERR /= ORB_OK) return

        bnorm = vx*vx+vy*vy+vz*vz
        if (.not.bwgt) bnorm = sqrt(bnorm)

        ydot(5) = ydot(4)*bnorm
      else
        ydot(4) = 1d0

        bnorm = vx*vx+vy*vy+vz*vz
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
          if (time-tprev >= 0.1*tmax .or. exit_loop) then
            if (cartesian_flag) then
              write (*,FMT="(a1)",ADVANCE="NO") 'C'
cc              call flush(6)
              cartesian_flag = .false.
            elseif (cross) then
              write (*,FMT="(a1)",ADVANCE="NO") 'X'
cc              call flush(6)
              cross = .false.
            else
              if (fwd) then
                write (*,FMT="(a1)",ADVANCE="NO") '+'
cc                call flush(6)
              else
                write (*,FMT="(a1)",ADVANCE="NO") '-'
cc                call flush(6)
              endif
            endif
            tprev = time
          endif

          if (time >= tmax) tprev = 0d0
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
        case(ORB_OK)       !No error
          continue
        case(ORB_ADJ_DT)   !Time step cut performed (not fatal)
          cross = .true.
        case(ORB_CAR_ST)   !Cartesian integration performed (not fatal)
          cartesian_flag = .true.   
        case(ORB_SML_DT)
          if (iiout > 0)
     .      write (*,*) 'Time step too small. Skipping orbit'
          exitloop = .true.
        case(ORB_OUT_DOM)
          if (iiout > 0)
     .      write (*,*) 'Out of domain. Skipping orbit'
          exitloop = .true.
        case(ORB_NEG_J)
          write (*,*) 'Negative Jacobian. Terminating'
          stop
        case(ORB_MAP_INV)
          if (iiout > 0)
     .      write (*,*) 'Could not invert map. Skipping orbit'
          exitloop = .true.
        case(ORB_SP_ERR)
          if (iiout > 0)
     .      write (*,*) 'Error around SP. Skipping orbit'
          exitloop = .true.
        case default
          write (*,*) 'Unknown error in orbit_step routine=',iterror
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

        if ((.not. lssx) .and. (.not. lssy) .and. (.not. lssz)
     $      .and.(.not.(lyap .or. auto_corr .or. connect))) then
          call pstop('ss_init'
     .              ,'Did not select any surface of section!!')
        endif

      end subroutine ss_init

cc!     ss_io_init
cc!     #################################################################
cc      subroutine ss_io_init(file,rz_ss)
cc
cc!     -----------------------------------------------------------------
cc!     Sets up XDRAW IO for surfaces of section
cc!     -----------------------------------------------------------------
cc
cc        implicit none
cc
cc!     Call variables
cc
cc        character(*) :: file
cc        logical, optional :: rz_ss
cc
cc!     Local variables
cc
cc        integer :: cpos
cc        character(50) :: drawin
cc
cc!     Begin program
cc
cc        ss_io = .true.   !Turn on SS XDRAW IO
cc
cc        if (PRESENT(rz_ss)) then
cc          rzss = rz_ss
cc        else
cc          rzss = .false.
cc        endif
cc
cc!        cpos = len_trim(ssx)-4
cc        cpos = index(trim(file),'.bin') - 1
cc
cc        ssx = 'ssx'//trim(file)
cc        ssy = 'ssy'//trim(file)
cc        ssz = 'ssz'//trim(file)
cc        ssrz= 'ssrz'//trim(file)
cc
cc        if (lssx) then
cc          drawin = 'draw'//ssx(1:cpos+3)//'.in'
cc          call createDrawInGfile(1,1,trim(ssx),'X surface of section'
cc     .                       ,'y',(/1/),(/'z'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc          nunitx=6027
cc          open(unit=nunitx,file=trim(ssx),form='unformatted')
cc        endif
cc
cc        if (lssy) then
cc          drawin = 'draw'//ssy(1:cpos+3)//'.in'
cc          call createDrawInGfile(1,1,trim(ssy),'Y surface of section'
cc     .                       ,'x',(/1/),(/'z'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc          nunity=6028
cc          open(unit=nunity,file=trim(ssy),form='unformatted')
cc        endif
cc
cc        if (lssz) then
cc          drawin = 'draw'//ssz(1:cpos+3)//'.in'
cc          call createDrawInGfile(1,1,trim(ssz),'Z surface of section'
cc     .                       ,'x',(/1/),(/'y'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc          nunitz=6029
cc          open(unit=nunitz,file=trim(ssz),form='unformatted')
cc        endif
cc
cc        if (rzss) then
cc          drawin = 'draw'//ssrz(1:cpos+4)//'.in'
cc          call createDrawInGfile(1,1,trim(ssrz),'RZ surface of section'
cc     .                       ,'r',(/1/),(/'z'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc          nunitrz=6030
cc          open(unit=nunitrz,file=trim(ssrz),form='unformatted')
cc        endif
cc
cc      end subroutine ss_io_init

cc!     ss_io_init
cc!     #################################################################
cc      subroutine ss_io_init(file,rz_ss)
cc
cc!     -----------------------------------------------------------------
cc!     Sets up XDRAW IO for surfaces of section
cc!     -----------------------------------------------------------------
cc
cc        implicit none
cc
cc!     Call variables
cc
cc        character(*) :: file
cc        logical, optional :: rz_ss
cc
cc!     Local variables
cc
cc        integer :: cpos
cc        character(50) :: drawin
cc
cc        character(50) :: ssx,ssy,ssz,ssrz
cc
cc!     Begin program
cc        ss_io = .true.   !Turn on SS XDRAW IO
cc
cc        if (PRESENT(rz_ss)) then
cc          rzss = rz_ss
cc        else
cc          rzss = .false.
cc        endif
cc
cc        cpos = index(trim(file),'.bin') - 1
cc	
cc        ssx = 'ssx'//trim(file)
cc        ssy = 'ssy'//trim(file)
cc        ssz = 'ssz'//trim(file)
cc        ssrz= 'ssrz'//trim(file)
cc	
cc        if (lssx) then
cc          drawin = 'draw'//ssx(1:cpos+3)//'.in'
cc          call createDrawInGfile(1,1,trim(ssx),'X surface of section'
cc     .                       ,'y',(/1/),(/'z'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc        endif
cc
cc        if (lssy) then
cc          drawin = 'draw'//ssy(1:cpos+3)//'.in'
cc          call createDrawInGfile(1,1,trim(ssy),'Y surface of section'
cc     .                       ,'x',(/1/),(/'z'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc        endif
cc
cc        if (lssz) then
cc          drawin = 'draw'//ssz(1:cpos+3)//'.in'
cc          call createDrawInGfile(1,1,trim(ssz),'Z surface of section'
cc     .                       ,'x',(/1/),(/'y'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc        endif
cc
cc        if (rzss) then
cc          drawin = 'draw'//ssrz(1:cpos+4)//'.in'
cc          call createDrawInGfile(1,1,trim(ssrz),'RZ surface of section'
cc     .                       ,'r',(/1/),(/'z'/),(/0/),(/0,0/)
cc     .                       ,drawin,.false.,.false.)
cc        endif
cc
cc        if (lssx) then
cc          nunitx = 27
cc          open(unit=nunitx,file=trim(ssx),form='unformatted')
cccc          call openBinFile(nunitx,trim(ssx))
cc        endif
cc
cc        if (lssy) then
cc          nunity = 1027
cc          open(unit=nunity,file=trim(ssy),form='unformatted')
cccc          call openBinFile(nunity,trim(ssy))
cc        endif
cc
cc        if (lssz) then
cc          nunitz = 2027
cc          open(unit=nunitz,file=trim(ssz),form='unformatted')
cccc          call openBinFile(nunitz,trim(ssz))
cc        endif
cc
cc        if (rzss) then
cc          nunitrz = 3027
cc          open(unit=nunitrz,file=trim(ssrz),form='unformatted')
cccc          call openBinFile(nunitrz,trim(ssrz))
cc        endif
cc
cc      end subroutine ss_io_init

!!     ss_io_init
!!     #################################################################
!      subroutine ss_io_init(file,rz_ss)
!
!!     -----------------------------------------------------------------
!!     Sets up XDRAW IO for surfaces of section
!!     -----------------------------------------------------------------
!
!        implicit none
!
!!     Call variables
!
!        character(*) :: file
!        logical, optional :: rz_ss
!
!!     Local variables
!
!        integer :: cpos,nfiles,ifile,ith,iproc,ierr
!        character(50) :: drawin,dummy
!
!        character(50),dimension(:),allocatable :: ssx,ssy,ssz,ssrz
!
!!     Begin program
!
!        cpos = index(trim(file),'.bin') - 1
!
!        nfiles = np*thr_tot
!        allocate(ssx(0:nfiles-1),ssy (0:nfiles-1)
!     .          ,ssz(0:nfiles-1),ssrz(0:nfiles-1))
!
!cc        write (*,*) 'DIAG -- ss_io_init',thr_tot,np
!
!        do ith=0,thr_tot-1
!          do iproc=0,np-1
!            ifile = ith+thr_tot*iproc
!            ssx(ifile) = 'ssx'//trim(file)//'_proc'//int2char(ifile)
!            ssy(ifile) = 'ssy'//trim(file)//'_proc'//int2char(ifile)
!            ssz(ifile) = 'ssz'//trim(file)//'_proc'//int2char(ifile)
!            ssrz(ifile)= 'ssrz'//trim(file)//'_proc'//int2char(ifile)
!          enddo
!        enddo
!
!        ss_io = .true.   !Turn on SS XDRAW IO
!  
!        if (PRESENT(rz_ss)) then
!          rzss = rz_ss
!        else
!          rzss = .false.
!        endif
!
!        if (my_rank == 0.and.thr_num==0) then
!          if (lssx) then
!            dummy = 'ssx'//trim(file)
!cc            ierr = rm_files(dummy)
!            drawin = 'draw'//dummy(1:cpos+3)//'.in'
!            call createDrawInGfile(1,1,ssx
!     .                         ,'X surface of section'
!     .                         ,'y',(/1/),(/'z'/),(/0/),(/0,0/)
!     .                         ,drawin,.false.,.false.)
!          endif
!          if (lssy) then
!            dummy = 'ssy'//trim(file)
!cc            ierr = rm_files(dummy)
!            drawin = 'draw'//dummy(1:cpos+3)//'.in'
!            call createDrawInGfile(1,1,ssy
!     .                         ,'Y surface of section'
!     .                         ,'x',(/1/),(/'z'/),(/0/),(/0,0/)
!     .                         ,trim(drawin),.false.,.false.)
!          endif
!          if (lssz) then
!            dummy = 'ssz'//trim(file)
!cc            ierr = rm_files(dummy)
!            drawin = 'draw'//dummy(1:cpos+3)//'.in'
!            call createDrawInGfile(1,1,ssz
!     .                         ,'Z surface of section'
!     .                         ,'x',(/1/),(/'y'/),(/0/),(/0,0/)
!     .                         ,drawin,.false.,.false.)
!          endif
!          if (rzss) then
!            dummy = 'ssrz'//trim(file)
!cc            ierr = rm_files(dummy)
!            drawin = 'draw'//dummy(1:cpos+4)//'.in'
!            call createDrawInGfile(1,1,ssrz
!     .                         ,'RZ surface of section'
!     .                         ,'r',(/1/),(/'z'/),(/0/),(/0,0/)
!     .                         ,drawin,.false.,.false.)
!          endif
!
!        endif
!
!c!$OMP CRITICAL
!        ifile = thr_num + thr_tot*my_rank
!        if (lssx) then
!          nunitx=6027+thr_num
!          open(unit=nunitx,file=trim(ssx(ifile))
!     .        ,form='unformatted')
!        endif
!
!        if (lssy) then
!          nunity=7027+thr_num
!          open(unit=nunity,file=trim(ssy(ifile))
!     .        ,form='unformatted')
!        endif
!
!        if (lssz) then
!          nunitz=8027+thr_num
!          open(unit=nunitz,file=trim(ssz(ifile))
!     .        ,form='unformatted')
!        endif
!
!        if (rzss) then
!          nunitrz=9027+thr_num
!          open(unit=nunitrz,file=trim(ssrz(ifile))
!     .        ,form='unformatted')
!        endif
!c!$OMP END CRITICAL
!
!!$OMP BARRIER
!
!        deallocate(ssx,ssy,ssz,ssrz)
!
!      end subroutine ss_io_init

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

        integer :: cpos,nfiles,ifile,ith,iproc,ierr,ic,thrtot

        character(50) :: drawin,dummy

        character(50),dimension(:),allocatable :: ssx,ssy,ssz,ssrz

c     Begin program

        thrtot = 1

#ifdef _OPENMP
!$OMP PARALLEL
        thrtot = omp_get_num_threads()
!$OMP END PARALLEL
#endif

        ss_io = .true.   !Turn on SS XDRAW IO
  
        if (PRESENT(rz_ss)) then
          rzss = rz_ss
        else
          rzss = .false.
        endif

c     Define surface of section files (one per OMP thread and MPI job)

        cpos = index(trim(file),'.bin') - 1

        nfiles = np*thrtot
        allocate(ssx(0:nfiles-1),ssy (0:nfiles-1)
     .          ,ssz(0:nfiles-1),ssrz(0:nfiles-1))

        do ith=0,thrtot-1
          do iproc=0,np-1
            ifile = ith+thrtot*iproc
            ssx(ifile) = 'ssx'//trim(file)//'_proc'//int2char(ifile)
            ssy(ifile) = 'ssy'//trim(file)//'_proc'//int2char(ifile)
            ssz(ifile) = 'ssz'//trim(file)//'_proc'//int2char(ifile)
            ssrz(ifile)= 'ssrz'//trim(file)//'_proc'//int2char(ifile)
          enddo
        enddo

c     Create draw.in files (inside MPI region, outside OMP threading)

        if (my_rank == 0) then
          if (lssx) then
            dummy = 'ssx'//trim(file)
cc            ierr = rm_files(dummy)
            drawin = 'draw'//dummy(1:cpos+3)//'.in'
            call createDrawInGfile(1,1,ssx
     .                         ,'X surface of section'
     .                         ,'y',(/1/),(/'z'/),(/0/),(/0,0/)
     .                         ,drawin,.false.,.false.)
          endif
          if (lssy) then
            dummy = 'ssy'//trim(file)
cc            ierr = rm_files(dummy)
            drawin = 'draw'//dummy(1:cpos+3)//'.in'
            call createDrawInGfile(1,1,ssy
     .                         ,'Y surface of section'
     .                         ,'x',(/1/),(/'z'/),(/0/),(/0,0/)
     .                         ,trim(drawin),.false.,.false.)
          endif
          if (lssz) then
            dummy = 'ssz'//trim(file)
cc            ierr = rm_files(dummy)
            drawin = 'draw'//dummy(1:cpos+3)//'.in'
            call createDrawInGfile(1,1,ssz
     .                         ,'Z surface of section'
     .                         ,'x',(/1/),(/'y'/),(/0/),(/0,0/)
     .                         ,drawin,.false.,.false.)
          endif
          if (rzss) then
            dummy = 'ssrz'//trim(file)
cc            ierr = rm_files(dummy)
            drawin = 'draw'//dummy(1:cpos+4)//'.in'
            call createDrawInGfile(1,1,ssrz
     .                         ,'RZ surface of section'
     .                         ,'R',(/1/),(/'Z'/),(/0/),(/0,0/)
     .                         ,drawin,.false.,.false.)
          endif

        endif

c     Open output files (one per thread and MPI process)

        allocate(nunitx(0:thrtot-1),nunity (0:thrtot-1))
        allocate(nunitz(0:thrtot-1),nunitrz(0:thrtot-1))

        if (lssx) then
          do ic = 0, thrtot - 1
            ifile = ic + thrtot*my_rank 
            nunitx(ic) = 6027+ic
            open(unit=nunitx(ic),file=trim(ssx(ifile))
     .           ,form='unformatted')
          enddo
        endif

        if (lssy) then
          do ic = 0, thrtot - 1
            ifile = ic + thrtot*my_rank 
            nunity(ic) = 7027+ic
            open(unit=nunity(ic),file=trim(ssy(ifile))
     .           ,form='unformatted')
          enddo
        endif

        if (lssz) then
          do ic = 0, thrtot - 1
            ifile = ic + thrtot*my_rank 
            nunitz(ic) = 8027+ic
            open(unit=nunitz(ic),file=trim(ssz(ifile))
     .           ,form='unformatted')
          enddo
        endif

        if (rzss) then
          do ic = 0, thrtot - 1
            ifile = ic + thrtot*my_rank 
            nunitrz(ic) = 9027+ic
            open(unit=nunitrz(ic),file=trim(ssrz(ifile))
     .           ,form='unformatted')
          enddo
        endif

        deallocate(ssx,ssy,ssz,ssrz)

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

        integer :: ith,thrtot

c     Begin program

        thrtot = 1

#ifdef _OPENMP
!$OMP PARALLEL
        thrtot = omp_get_num_threads()
!$OMP END PARALLEL
#endif

        do ith=0,thrtot-1
          if (lssx) close(nunitx(ith))
          if (lssy) close(nunity(ith))
          if (lssz) close(nunitz(ith))
          if (rzss) close(nunitrz(ith))
        enddo

        deallocate(nunitx,nunity,nunitz,nunitrz)

      end subroutine ss_io_end

c     sschk
c     #################################################################
      subroutine sschk(tv,to,tt,xvold,xold,x,yvold,yold,y
     .                ,zvold,zold,z,xss,yss,zss,tss)

c     -----------------------------------------------------------------
c     Checks whether orbit intersects specified surfaces of section
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN) :: tv,to,tt
     .                       ,xvold,xold,x
     .                       ,yvold,yold,y
     .                       ,zvold,zold,z
        real(8),intent(OUT):: xss,yss,zss,tss

c     Local variables

        real(8) :: rr
        integer :: iwhich
        integer :: thrnum,thrtot

c     Begin program

        call get_omp_thread_id(thrnum,thrtot)

        xss = 0d0 ; yss = 0d0 ; zss = 0d0 ; tss = 0d0

c     Lower order of interpolation for SS if derivatives are too different

        if (ssorder == 0) then
          call pstop('sschk','Need to call ss_init!')
        endif

        iwhich = ssorder

        if (to > tv) then
          if (tt > to) then
            if(    (     abs(xold-xvold)/(to-tv)
     .              >5d0*abs(x    -xold)/(tt-to))
     .         .or.(     abs(yold-yvold)/(to-tv)
     .              >5d0*abs(y    -yold)/(tt-to))
     .         .or.(     abs(zold-zvold)/(to-tv)
     .              >5d0*abs(z    -zold)/(tt-to)) ) iwhich = 1
          else
            iwhich = 1
          endif
        else
          iwhich = 1
        endif

c     Check for intersection with X surface of section

        if (lssx
     .      .and.chk_ss_pos(xss0,xold,x,obcnd(1),xmin,xmax,xss)) then

          call sursec(yvold,zvold,xvold,tv
     .               ,yold ,zold ,xold ,to
     .               ,y    ,z    ,x    ,tt
     .               ,xss,npointx,yss,zss,tss,iwhich)

          if (ss_io) then
            call chk_pos(xss,yss,zss) !this resets logical coords
            write(nunitx(thrnum)) real(yss),real(zss)
          endif

          !Accumulate B diagnostics
          if (bdiag) call B_diag_update(yss,y00,tss)
        endif

c     Check for intersection with Y surface of section

        if (lssy
     .      .and.chk_ss_pos(yss0,yold,y,obcnd(3),ymin,ymax,yss)) then

          call sursec(zvold,xvold,yvold,tv
     .               ,zold ,xold ,yold ,to
     .               ,z    ,x    ,y    ,tt
     .               ,yss,npointy,zss,xss,tss,iwhich)

          if (ss_io) then
            call chk_pos(xss,yss,zss) !this resets logical coords
            
            if (flux) then
              write(nunity(thrnum)) real(0.5*xss**2),real(zss)
            else
              write(nunity(thrnum)) real(xss),real(zss)
            endif
          endif

          !Accumulate B diagnostics
          if (bdiag) call B_diag_update(xss,x00,tss)
        endif

c     Check for intersection with Z surface of section

        if ((lssz.or.rzss)
     .      .and.chk_ss_pos(zss0,zold,z,obcnd(5),zmin,zmax,zss)) then

          call sursec(xvold,yvold,zvold,tv
     .               ,xold ,yold ,zold ,to
     .               ,x    ,y    ,z    ,tt
     .               ,zss,npointz,xss,yss,tss,iwhich)

          if (ss_io) then
            if (rzss) then   !Find R-Z plane positions
              call evalX(xss,yss,zss,xcar,ycar,zcar) !Eval Cartesian coords; also resets logical coords
              rr = sqrt(xcar**2+ycar**2)
              write(nunitrz(thrnum)) real(rr),real(zcar)
            endif

            if (lssz) then
              call chk_pos(xss,yss,zss) !this resets logical coords

              if (flux) then
                write(nunitz(thrnum)) real(0.5*xss**2),real(yss) !In flux
              else
                write(nunitz(thrnum)) real(xss)       ,real(yss) !In radius
              endif
            endif
          endif

          !Accumulate B diagnostics
          if (bdiag) call B_diag_update(xss,x00,tss)
        endif

      end subroutine sschk

c     sursec
c     ##################################################################
      subroutine sursec(x0,y0,z0,t0
     .                 ,x1,y1,z1,t1
     .                 ,x2,y2,z2,t2
     .                 ,zsurf,npoints,xs,ys,ts,order)

c     ------------------------------------------------------------------
c     Surface of section routine in Z direction. Allows linear (order=1)
c     and quadratic (order=2) interpolation of time orbits. Returns (x,y)
c     and arc-length at surface of section in xs,ys,ts respectively.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: x0,y0,z0,t0,x1,y1,z1,t1,x2,y2,z2,t2
      real(8) :: zsurf,xs,ys,ts
      integer :: npoints,order

c     Local variables

      real(8) :: a11,a12,a21,a22,det
     .          ,b(order,order),vec(order)
     .          ,xcf(order),ycf(order),tcf(order)

c     Begin program

      select case(order)
      case(1)

c       Linear interpolation. x0,y0,z0 not used.

        ts=t1+(zsurf-z1)*(t2-t1)/(z2-z1)
        xs=x1+(zsurf-z1)*(x2-x1)/(z2-z1)
        ys=y1+(zsurf-z1)*(y2-y1)/(z2-z1)
cc        xs=x1+(x-x1)*(ts-t1)/(tt-t1)
cc        ys=y1+(y-y1)*(ts-t1)/(tt-t1)

      case(2)

c       Quadratic interpolation. Expand about z1.  
c       Vandermonde matrix with z0, z1, z
c       Careful if time stepping is made nonuniform!

        a11= z0-z1
        a12=(z0-z1)**2
        a21= z2-z1
        a22=(z2-z1)**2

        det=a11*a22-a12*a21

        if (abs(det) < 1d-15) then
          write (*,*) 'z0,z1,z2=',z0,z1,z2
          write (*,*) 'Singular Vandermonde matrix'
          write (*,*) 'Aborting in routine sursec'
          stop
        endif

        b(1,1)= a22/det
        b(1,2)=-a12/det
        b(2,1)=-a21/det
        b(2,2)= a11/det

        vec = (/ x0-x1, x2-x1 /)
        xcf = matmul(b,vec)

        vec = (/ y0-y1, y2-y1 /)
        ycf = matmul(b,vec)

        vec = (/ t0-t1, t2-t1 /)
        tcf = matmul(b,vec)

        vec = (/ zsurf-z1, (zsurf-z1)**2 /)
        ts = t1 + dot_product(tcf,vec)
        xs = x1 + dot_product(xcf,vec)
        ys = y1 + dot_product(ycf,vec)

c       Quadratic interpolation. Expand about z1.  
c       Vandermonde matrix with z0, z1, z
c       Careful if time stepping is made nonuniform!

cc        a11=t0-t1
cc        a12=(t0-t1)**2
cc        a21=tt-t1
cc        a22=(tt-t1)**2
cc
cc        det=a11*a22-a12*a21
cc
cc        b11=a22/det
cc        b12=-a12/det
cc        b21=-a21/det
cc        b22=a11/det
cc
cc        bbb=b11*(z0-z1)+b12*(z-z1)
cc        ccc=b21*(z0-z1)+b22*(z-z1)
cc        aaa=z1
cc
cc        ts=(zsurf-z1)/bbb
cc        ts=t1+ts           !Linear time
cc        if (ccc /= 0d0) then   !Calculate quadratic correction
cc          if(zsurf.ge.z1) then
cc            ts=(-bbb+sqrt(bbb**2-4.*ccc*(z1-zsurf)))/(2.*ccc)
cc            ts=t1+ts
cc          else
cc            ts=(-bbb-sqrt(bbb**2-4.*ccc*(z1-zsurf)))/(2.*ccc)
cc            ts=t1+ts
cc          endif
cc        endif
cc
cc        bbb=b11*(x0-x1)+b12*(x-x1)
cc        ccc=b21*(x0-x1)+b22*(x-x1)
cc        aaa=x1
cc        xs=aaa+bbb*(ts-t1)+ccc*(ts-t1)**2
cc
cc        bbb=b11*(y0-y1)+b12*(y-y1)
cc        ccc=b21*(y0-y1)+b22*(y-y1)
cc        aaa=y1
cc        ys=aaa+bbb*(ts-t1)+ccc*(ts-t1)**2

      case default

        write (*,*)
     .     'Order of accuracy in SS interpolation not implemented'
        stop

      end select

      npoints=npoints+1

c     End program

      end subroutine sursec

c     sschk_ho
c     #################################################################
      subroutine sschk_ho(xss,yss,zss,tss)

c     -----------------------------------------------------------------
c     Checks whether orbit intersects specified surfaces of section.
c     Returns (x,y,z,t) at surface of section
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(OUT):: xss,yss,zss,tss

c     Local variables

        real(8) :: rr
        integer :: iwhich,iierr
        integer :: thrnum,thrtot

c     Begin program

        call get_omp_thread_id(thrnum,thrtot)

        xss = 0d0 ; yss = 0d0 ; zss = 0d0 ; tss = 0d0

c     Lower order of interpolation for SS if derivatives are too different

        if (ssorder == 0) then
          call pstop('sschk','Need to call ss_init!')
        endif

        iwhich = ssorder

c     Check for intersection with X surface of section

        if (lssx
     .      .and.chk_ss_pos(xss0,xold,x,obcnd(1),xmin,xmax,xss)) then

          call sursec_ho(xss,xhist,yhist,zhist,shist
     .                  ,npointx,yss,zss,tss,iwhich)

          if (ss_io) then
           if (auto_corr) then
            !update number of intersections
            call auto_corr_intersx_update(intersect_x)
           endif

            mx = mx + 1
            if (.not. ridge .and. .not. auto_corr) then
             call chk_pos(xss,yss,zss) !this resets logical coords
!#! Marco, november 2017, insert the ierr in chk_pos for consistency
!with t3d which needs to stop when orbit exits domain
!             call chk_pos(xss,yss,zss) !this resets logical coords
             call chk_pos(xss,yss,zss,ierr=iierr) !this resets logical coords
             if (iierr .ne. 0) then
              return
             endif
             write(nunitx(thrnum)) real(yss),real(zss)
            endif
          endif

          !Accumulate B diagnostics
          if (bdiag) call B_diag_update(yss,y00,tss)
        endif

c     Check for intersection with Y surface of section

        if (lssy
     .      .and.chk_ss_pos(yss0,yold,y,obcnd(3),ymin,ymax,yss)) then

          call sursec_ho(yss,yhist,zhist,xhist,shist
     .                  ,npointy,zss,xss,tss,iwhich)

         if (ss_io) then
           if (auto_corr) then
            !update number of intersections
            call auto_corr_intersy_update(intersect_y)
           endif

!#!Marco, update intersection counter
          my = my + 1

          if (ridge) then 
           if ((my .lt. max_inter_y) 
     .         .and. (mod(my,delta_y) .eq. 0)) then
!            call chk_pos(xss,yss,zss,ierr=iierr) !this resets logical coords
            if (xss .gt. xmax) then
!             write(*,*) 'DIAGyridge',my,delta_y,xss,zss
             call ridge_write(1.d0,0.d0,thrnum)
            else
!#!Marco, ridge write
            call ridge_write(xss,zss,thrnum)
            endif
           endif
          endif
 
          if (.not. ridge .and. .not. auto_corr) then
           call chk_pos(xss,yss,zss) !this resets logical coords
!#! Marco, november 2017, insert the ierr in chk_pos for consistency
!with t3d which needs to stop when orbit exits domain
!             call chk_pos(xss,yss,zss) !this resets logical coords
             call chk_pos(xss,yss,zss,ierr=iierr) !this resets logical coords
             if (iierr .ne. 0) then
              return
             endif
           if (flux) then
             write(nunity(thrnum)) real(0.5*xss**2),real(zss)
           else
             write(nunity(thrnum)) real(xss),real(zss)
           endif
          endif
         endif

         !Accumulate B diagnostics
         if (bdiag) call B_diag_update(xss,x00,tss)
        endif

c     Check for intersection with Z surface of section

        if ((lssz.or.rzss)
     .      .and.chk_ss_pos(zss0,zold,z,obcnd(5),zmin,zmax,zss)) then

          call sursec_ho(zss,zhist,xhist,yhist,shist
     .                  ,npointz,xss,yss,tss,iwhich)

          if (ss_io) then
            if (rzss) then   !Find R-Z plane positions
             if (auto_corr) then
              !update number of intersections
              call auto_corr_interstor_update(intersect_tor)
             endif

             if (.not. ridge .and. .not. auto_corr) then
              call evalX(xss,yss,zss,xcar,ycar,zcar) !Eval Cartesian coords; also resets logical coords
              rr = sqrt(xcar**2+ycar**2)
              write(nunitrz(thrnum)) real(rr),real(zcar)
             endif
            endif

           if (lssz) then
            if (auto_corr) then
             !update number of intersections
             call auto_corr_intersz_update(intersect_z)
!             write(*,*) 'DIAGac0',thrnum,intersect_z
            endif

!#!Marco, update intersection counter
            mz = mz + 1
!            write(*,*) 'DIAGz0',thrnum

            if (ridge) then 
             if ((mz .lt. max_inter_z)
     .           .and. (mod(mz,delta_z) .eq. 0)) then
                if (thrnum .eq. 0) then
!                 write(*,*) 'DIAGridge0',mz,delta_z,xss,yss
                endif
              if (xss .gt. xmax) then
!               write(*,*) 'DIAGridge',xss,xmax
               call ridge_write(1.d0,0.d0,thrnum)
              else
!#!Marco, ridge write
!               write(*,*) 'DIAGridge',thrnum
               call ridge_write(xss,yss,thrnum)
              endif
             endif

            endif
 
            if (.not. ridge .and. .not. auto_corr) then
!#! Marco, november 2017, insert the ierr in chk_pos for consistency
!with t3d which needs to stop when orbit exits domain
!             call chk_pos(xss,yss,zss) !this resets logical coords
             call chk_pos(xss,yss,zss,ierr=iierr) !this resets logical coords
             if (iierr .ne. 0) then
              return
             endif
             if (flux) then
               write(nunitz(thrnum)) real(0.5*xss**2),real(yss) !In flux
             else
               write(nunitz(thrnum)) real(xss)       ,real(yss) !In radius
             endif

            endif
           endif
          endif

          !Accumulate B diagnostics
          if (bdiag) call B_diag_update(xss,x00,tss)
        endif

      end subroutine sschk_ho

c     sursec_ho
c     ##################################################################
      subroutine sursec_ho(zsurf,z,x,y,t,npoints,xs,ys,ts,order)

c     ------------------------------------------------------------------
c     Surface of section routine in Z direction.  Returns (x,y) and
c     arc-length/time at surface of section in xs,ys,ts
c     respectively. Allows linear (order=1) and quadratic (order=2)
c     interpolation of time orbits. Arrays are ordered with most recent
c     point first (first position), and older points later.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: npoints,order
      real(8) :: x(:),y(:),z(:),t(:)
      real(8) :: zsurf,xs,ys,ts

c     Local variables

      integer :: i,j,ii,lorder
      real(8),allocatable,dimension(:,:) :: a,rhs,sol
      real(8),allocatable,dimension(:)   :: vec

c     Begin program

      lorder = min(order,count(z /= UNDEF_DP)-1)

cc      write (*,*) count(t /= UNDEF_DP),lorder

      select case(lorder)
      case(1)

c       Linear interpolation. x0,y0,z0 not used.

        ts=t(2)+(zsurf-z(2))*(t(1)-t(2))/(z(1)-z(2))
        xs=x(2)+(zsurf-z(2))*(x(1)-x(2))/(z(1)-z(2))
        ys=y(2)+(zsurf-z(2))*(y(1)-y(2))/(z(1)-z(2))

      case default

        allocate(a  (lorder,lorder)
     .          ,rhs(lorder,3)
     .          ,sol(lorder,3)
     .          ,vec(lorder))

c       Quadratic interpolation. Expand about z(2)
c       Vandermonde matrix
c       Careful if time stepping is made nonuniform!

        do i=1,lorder
          do j=1,lorder
            ii = i
            if (i>=2) ii=i+1
            a(i,j) = (z(ii)-z(2))**j
          enddo
        enddo

        rhs(1,1) = x(1)-x(2)
        rhs(1,2) = y(1)-y(2)
        rhs(1,3) = t(1)-t(2)
        rhs(2:lorder,1) = x(3:lorder+1)-x(2)
        rhs(2:lorder,2) = y(3:lorder+1)-y(2)
        rhs(2:lorder,3) = t(3:lorder+1)-t(2)

        call blockSolve(lorder,a,3,rhs,sol)

        do i=1,lorder
          vec(i) = (zsurf-z(2))**i
        enddo

        xs = x(2) + dot_product(sol(:,1),vec)
        ys = y(2) + dot_product(sol(:,2),vec)
        ts = t(2) + dot_product(sol(:,3),vec)

        deallocate(a,rhs,sol,vec)
cc      case default
cc
cc        write (*,*)
cc     .     'Order of accuracy in SS interpolation not implemented'
cc        stop

      end select

      npoints=npoints+1

c     End program

      end subroutine sursec_ho

c     chk_ss_pos
c     ##################################################################
      function chk_ss_pos(x0,x1,x2,bc,xmin,xmax,xss) result(found)

c     ------------------------------------------------------------------
c     Determines whether surface of section at x0 is between x1 and x2.
c     It automatically accounts for periodic BCs. If periodic BCs, it
c     returns shifted surface of section location, xss, such that
c                  x1 < xss < x2
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
cc        p1 = floor((x1-xmin)/(xmax-xmin))
cc        p2 = floor((x2-xmin)/(xmax-xmin))
        p1 = floor((x1-x0)/(xmax-xmin))
        p2 = floor((x2-x0)/(xmax-xmin))

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

cc      write (*,*)'chk_ss_pos',x0,x00,x1,x2,x11,x22,p1,p2,found
cc      pause

      if ((x11-x00)*(x22-x00) <= 0d0) xss = x00+pp*(xmax-xmin)
      if ((x11-x0 )*(x22-x0 ) <= 0d0) xss = x0 +pp*(xmax-xmin)

      end function chk_ss_pos

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

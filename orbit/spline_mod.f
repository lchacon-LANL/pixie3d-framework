c module spline_field
c #####################################################################
      module spline_field

        use xdraw_io

!#! Marco
        use lyapn

        use grid, ONLY:pstop,my_rank,grid_mg_def

        use bc_def, bcond2 => bcond

        real(8)  :: db3val
        external :: db3val

        integer :: istep=0

        !Error codes
        integer, parameter :: ORB_OK=0,ORB_ADJ_DT=-1,ORB_CAR_ST=-2  !Non fatal
     .                       ,ORB_SML_DT=1,ORB_OUT_DOM=2,ORB_NEG_J=4!Fatal
     .                       ,ORB_MAP_INV=5,ORB_SP_ERR=6            !Fatal
     

        !Private variables
        integer,private :: nx,ny,nz,ag,sbcnd(6)

        real(8),private :: xsmin,xsmax,ysmin,ysmax,zsmin,zsmax
        integer,private :: kx,ky,kz,dim,dime,flg
        real(8),private,dimension(:),allocatable :: tx,ty,tz,work
     .                                             ,xs,ys,zs,worke

        real(8),private,dimension(:,:,:),pointer:: jcoef   => null()
     .                                            ,acoefx  => null()
     .                                            ,acoefy  => null()
     .                                            ,acoefz  => null()
     .                                            ,fldcoef => null()

        real(8),private,dimension(:,:,:,:),pointer :: xcoef => null()
     .                                               ,bcoef => null()
     .                                               ,bcarcoef => null()

        logical, private :: x_is_log,y_is_log,z_is_log

        INTERFACE setupSplines
          module procedure setupSplines_gst,setupSplines_nogst
        END INTERFACE

      contains

c     set_sp_domain_limits
c     ##################################################################
      subroutine set_sp_domain_limits(x,bc,xmin,xmax)

c     ------------------------------------------------------------------
c     Set domain limits according to dimension vector x and boundary
c     condition bc.
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: bc
        real(8) :: x(:),xmin,xmax

c     Local Variables

        integer :: nn

c     Begin program

        nn = size(x)

        if (bc == PER) then  !Ghost cell is at x(1)
cc          xmin = 5d-1*(x(1)+x(2))
cc          xmax = 5d-1*(x(nn-1)+x(nn))
          x = x - x(2)  !Set angle = 0 at x(2)
          xmin = x(2)
          xmax = x(nn)
        else
          xmin = x(1)
          xmax = x(nn)
        endif

c     End program

      end subroutine set_sp_domain_limits

c     get_sp_domain_limits
c     ##################################################################
      subroutine get_sp_domain_limits(xmin,xmax,ymin,ymax,zmin,zmax)

c     ------------------------------------------------------------------
c     Get spline domain limits.
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: xmin,xmax,ymin,ymax,zmin,zmax

c     Local Variables

c     Begin program

        xmin = xsmin
        ymin = ysmin
        zmin = zsmin

        xmax = xsmax
        ymax = ysmax
        zmax = zsmax

c     End program

      end subroutine get_sp_domain_limits

c     chk_pos
c     ##################################################################
      subroutine chk_pos(x,y,z,no_per_bc,warning,ierr)

c     ------------------------------------------------------------------
c     Check whether we are within LOGICAL domain, and if not:
c       * If periodic, return position within first period
c       * Otherwise:
c          + If ierr is not present, terminate with error
c          + Else, return ierr /= 0.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: x,y,z
      logical,optional :: no_per_bc,warning
      integer,optional :: ierr

c     Local Variables

      integer :: ierror
      logical :: per_bc,message

c     Begin program

      if (PRESENT(no_per_bc)) then
        per_bc = .not.no_per_bc
      else
        per_bc = .true.
      endif

      if (PRESENT(warning)) then
        message = warning
      else
        message = .true.
      endif

      ierror = ORB_OK

      if (x > xsmax) then
        if (sbcnd(2) == PER) then   !Periodic BC
          if (per_bc) x = xsmin + mod(x-xsmin,(xsmax-xsmin))
        else
          ierror = ORB_OUT_DOM
        endif
      endif

      if (x < 0d0 .and. sbcnd(1) == SP) then !Singular point BC
        x =-x
        y = y + acos(-1d0)
      endif

      if (x < xsmin) then
        if (sbcnd(1) == PER) then   !Periodic BC
          if (per_bc) x = xsmax - mod(xsmin-x,xsmax-xsmin)
        else
          ierror = ORB_OUT_DOM
        endif
      endif

      if (y > ysmax) then
        if (sbcnd(4) == PER) then   !Periodic BC
          if (per_bc) y = ysmin + mod(y-ysmin,(ysmax-ysmin))
        else
          ierror = ORB_OUT_DOM
        endif
      endif

      if (y < ysmin) then
        if (sbcnd(3) == PER) then   !Periodic BC
          if (per_bc) y = ysmax - mod(ysmin-y,ysmax-ysmin)
        else
          ierror = ORB_OUT_DOM
        endif
      endif

      if (z > zsmax) then
        if (sbcnd(6) == PER) then   !Periodic BC
          if (per_bc) z = zsmin + mod(z-zsmin,(zsmax-zsmin))
        else
          ierror = ORB_OUT_DOM
        endif
      endif

      if (z < zsmin) then
        if (sbcnd(5) == PER) then   !Periodic BC
          if (per_bc) z = zsmax - mod(zsmin-z,zsmax-zsmin)
        else
          ierror = ORB_OUT_DOM
        endif
      endif

cc      if (.not.PRESENT(ierr)) then
        if (ierror /= 0.and.message) then
          write (*,*) 
          write (*,*) 'Error in chk_pos: out of domain!'
          write (*,*)
          write (*,*) 'Current position',x,y,z
          write (*,*)
          write (*,*) 'X domain limits=',xsmin,xsmax
          write (*,*) 'Y domain limits=',ysmin,ysmax
          write (*,*) 'Z domain limits=',zsmin,zsmax
          write (*,*)
          write (*,*) 'BCs',sbcnd
cc          stop
        endif
cc      else
cc        ierr = ierror
cc      endif

        if (PRESENT(ierr)) then
          ierr = ierror
        else
          if (ierror /= 0) stop
        endif

c     End program

      end subroutine chk_pos

c     setupSplines_gst
c     #################################################################
      subroutine setupSplines_gst(g_def,igrid,order,bcnd)
c     -----------------------------------------------------------------
c     This routine sets up 3D splines, including allocation of memory
c     space.
c     -----------------------------------------------------------------

        implicit none            ! For safe Fortran

c     Call variables

        type(grid_mg_def) :: g_def

        integer :: igrid,order
        integer :: bcnd(6)

c     Local variables

        integer :: alloc_stat,i

c     Begin program

        sbcnd = bcnd

c     Initialize private variables

        nx = g_def%nxgl(igrid) + 2
        ny = g_def%nygl(igrid) + 2
        nz = g_def%nzgl(igrid) + 2

        allocate(xs(nx),ys(ny),zs(nz))

c     Initialize spline domain arrays

        xs = g_def%xg
        ys = g_def%yg
        zs = g_def%zg

c     Define domain limits

        xsmin = g_def%gxmin
        xsmax = g_def%gxmax
          
        ysmin = g_def%gymin
        ysmax = g_def%gymax
          
        zsmin = g_def%gzmin
        zsmax = g_def%gzmax

c     Prepare 3d spline interpolation

        flg = 0 !Let spline routine find interpolation knots
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim  = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(work(dim),stat=alloc_stat)
        allocate(tx(nx+kx),stat=alloc_stat)
        allocate(ty(ny+ky),stat=alloc_stat)
        allocate(tz(nz+kz),stat=alloc_stat)

!$omp parallel
        call set_omp_thread_id()
!$omp end parallel

	dime = ky*kz + 3*max(kx,ky,kz) + kz
	allocate(worke(dime*thr_tot),stat=alloc_stat)

c     End program

      end subroutine setupSplines_gst

c     setupSplines_nogst
c     #################################################################
      subroutine setupSplines_nogst(nnx,nny,nnz,xx,yy,zz,order,bcnd)
c     -----------------------------------------------------------------
c     This routine sets up 3D splines, including allocation of memory
c     space.
c     -----------------------------------------------------------------

        implicit none            ! For safe Fortran

c     Call variables

        integer :: nnx,nny,nnz,order
        real(8) :: xx(nnx),yy(nny),zz(nnz)
        integer :: bcnd(6)

c     Local variables

        integer :: alloc_stat,i!,get_omp_numthreads

c     Begin program

        sbcnd = bcnd

c     Initialize private variables

        nx = nnx
        ny = nny
        nz = nnz

        allocate(xs(nx),ys(ny),zs(nz))

c     Initialize spline domain arrays

        xs = xx
        ys = yy
        zs = zz

c     Define spline domain limits

        call set_sp_domain_limits(xs,sbcnd(1),xsmin,xsmax)
        call set_sp_domain_limits(ys,sbcnd(3),ysmin,ysmax)
        call set_sp_domain_limits(zs,sbcnd(5),zsmin,zsmax)

c     Prepare 3d spline interpolation

        flg = 0  !Let spline routine find interpolation knots

        if (order == 0) call pstop('setupSplines',"spline order = 0")
       
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim  = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(work(dim),stat=alloc_stat)
        allocate(tx(nx+kx),stat=alloc_stat)
        allocate(ty(ny+ky),stat=alloc_stat)
        allocate(tz(nz+kz),stat=alloc_stat)

!$omp parallel
        call set_omp_thread_id()
!$omp end parallel                                                             

        dime = ky*kz + 3*max(kx,ky,kz) + kz
        allocate(worke(dime*thr_tot),stat=alloc_stat)

c     End program

      end subroutine setupSplines_nogst
      
c     splineA
c     #################################################################
      subroutine splineA(bx,by,bz,a_gauge,input_is_A)
c     -----------------------------------------------------------------
c     This routine splines up vector components (ax,ay,az) on a mesh
c     of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer    :: a_gauge
      real(8)    :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
      logical,optional :: input_is_A

c     Local variables

      integer :: alloc_stat,i
      logical :: inputA

      real(8) :: ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
     .          ,sumax,sumay,sumaz

c     Begin program

      if (PRESENT(input_is_A)) then
        inputA = input_is_A
      else
        inputA = .false.
      endif

      ag = a_gauge

c     Get vector potential

      if (inputA) then
        ax = bx
        ay = by
        az = bz
      else
        call getA_on_mesh(nx,ny,nz,xs,ys,zs,bx,by,bz,ax,ay,az,ag)
      endif

c     Determine integration steps

      sumax = sum(ax**2)
      sumay = sum(ay**2)
      sumaz = sum(az**2)

      select case(ag)
      case(1)
        if (sumay == 0d0) istep=2
        if (sumaz == 0d0) istep=1
      case(2)
        if (sumax == 0d0) istep=2
        if (sumaz == 0d0) istep=1
      end select

c     Spline vector potential

      select case(ag)
      case(1) !Ax=0
        !ay
        allocate(acoefy(nx,ny,nz),stat=alloc_stat)

        call db3ink(xs,nx,ys,ny,zs,nz
     .             ,ay,nx,ny,kx,ky,kz,tx,ty,tz,acoefy,work,flg)
      case(2) !Ay=0
        !ax
        allocate(acoefx(nx,ny,nz),stat=alloc_stat)

        call db3ink(xs,nx,ys,ny,zs,nz
     .             ,ax,nx,ny,kx,ky,kz,tx,ty,tz,acoefx,work,flg)
      end select

      !az
      allocate(acoefz(nx,ny,nz),stat=alloc_stat)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,az,nx,ny,kx,ky,kz,tx,ty,tz,acoefz,work,flg)

c     End program

      end subroutine splineA

c     getA_on_mesh
c     ###############################################################
      subroutine getA_on_mesh(nx,ny,nz,xx,yy,zz,bx,by,bz,ax,ay,az,gauge)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential "a" from CONTRAVARIANT vector
c     field components "b". Employs either gauge A_1=0 or A_2=0d0,
c     as specified in gauge.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,gauge
      real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
      real(8) :: ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
      real(8) :: xx(0:nx-1),yy(0:ny-1),zz(0:nz-1)

c     Local variables

      integer    :: i,j,k,ii,jj,ig,jg,kg,ivar,nxl,nyl,nzl
      real(8)    :: cm,c0,cp,dxx,dyy,dzz,dvol,cov(3)
      real(8)    :: a(0:nx-1,0:ny-1,0:nz-1,3)
     .             ,b(0:nx-1,0:ny-1,0:nz-1,3)

      logical    :: spoint

c     Begin program

      nxl = nx-2
      nyl = ny-2
      nzl = nz-2

      b(:,:,:,1) = bx
      b(:,:,:,2) = by
      b(:,:,:,3) = bz

      !Account for collapsed dimensions
      if (nxl == 1) gauge = 2
      if (nyl == 1) gauge = 3
      if (nzl == 1) gauge = 1

      spoint = bcSP()

      a=0d0

c     Accommodate different gauges

      select case(gauge)
      case(1) !A1=0 code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if (spoint) then
          !Start at x=0 to enforce regularity conditions in A
          call line_int_xz(1,nxl,1,nzl
     .                 ,b(0:nxl+1,0:nyl+1,0:nzl+1,1:3)
     .                 ,a(0:nxl+1,0:nyl+1,0:nzl+1,1:3))
        else
          !Lower-left quadrant
          call line_int_xz(nxl/2-1,1,nzl/2-1,1
     .                 ,b(0:nxl/2,0:nyl+1,0:nzl/2,1:3)
     .                 ,a(0:nxl/2,0:nyl+1,0:nzl/2,1:3))

          !Lower-right quadrant
          call line_int_xz(nxl/2+1,nxl,nzl/2-1,1
     .                 ,b(nxl/2:nxl+1,0:nyl+1,0:nzl/2,1:3)
     .                 ,a(nxl/2:nxl+1,0:nyl+1,0:nzl/2,1:3))

          !Upper-left quadrant
          call line_int_xz(nxl/2-1,1,nzl/2+1,nzl
     .                 ,b(0:nxl/2,0:nyl+1,nzl/2:nzl+1,1:3)
     .                 ,a(0:nxl/2,0:nyl+1,nzl/2:nzl+1,1:3))

          !Upper-right quadrant
          call line_int_xz(nxl/2+1,nxl,nzl/2+1,nzl
     .                 ,b(nxl/2:nxl+1,0:nyl+1,nzl/2:nzl+1,1:3)
     .                 ,a(nxl/2:nxl+1,0:nyl+1,nzl/2:nzl+1,1:3))
        endif

      case(2) !A2=0 code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c       Quadrant-based code

        !Lower-left quadrant
        call line_int_yz(nyl/2-1,1,nzl/2-1,1
     .               ,b(0:nxl+1,0:nyl/2,0:nzl/2,1:3)
     .               ,a(0:nxl+1,0:nyl/2,0:nzl/2,1:3))

        !Lower-right quadrant
        call line_int_yz(nyl/2+1,nyl,nzl/2-1,1
     .               ,b(0:nxl+1,nyl/2:nyl+1,0:nzl/2,1:3)
     .               ,a(0:nxl+1,nyl/2:nyl+1,0:nzl/2,1:3))

        !Upper-left quadrant
        call line_int_yz(nyl/2-1,1,nzl/2+1,nzl
     .               ,b(0:nxl+1,0:nyl/2,nzl/2:nzl+1,1:3)
     .               ,a(0:nxl+1,0:nyl/2,nzl/2:nzl+1,1:3))

        !Upper-right quadrant
        call line_int_yz(nyl/2+1,nyl,nzl/2+1,nzl
     .               ,b(0:nxl+1,nyl/2:nyl+1,nzl/2:nzl+1,1:3)
     .               ,a(0:nxl+1,nyl/2:nyl+1,nzl/2:nzl+1,1:3))

      case default
        call pstop('getA','Gauge not implemented for ny=1')
      end select

c diag ****
cc      write (*,*) 'DIAG -- getA_on_mesh',nx,ny,nz
cc      open(unit=110,file='vecpot.bin',form='unformatted'
cc     .      ,status='replace')
cc      call contour(a(1:nxl+1,2,1:nzl+1,1)
cc     .            ,nx-1,nz-1,0d0,xmax,0d0,zmax,0,110)
cc      call contour(a(1:nxl+1,2,1:nzl+1,2)
cc     .            ,nx-1,nz-1,0d0,xmax,0d0,zmax,1,110)
cc      call contour(a(1:nxl+1,2,1:nzl+1,3)
cc     .            ,nx-1,nz-1,0d0,xmax,0d0,zmax,1,110)
cc
cccc      call contour(a(1:nxl+1,1:nyl+1,1,1)
cccc     .            ,nx-1,ny-1,0d0,xmax,0d0,ymax,0,110)
cccc      call contour(a(1:nxl+1,1:nyl+1,1,2)
cccc     .            ,nx-1,ny-1,0d0,xmax,0d0,ymax,1,110)
cccc      call contour(a(1:nxl+1,1:nyl+1,1,3)
cccc     .            ,nx-1,ny-1,0d0,xmax,0d0,ymax,1,110)
cc      close(110)
cccc      stop
c diag ****

c     Take care of collapsed dimensions (no motion along those)

      if (nxl == 1) a(:,:,:,2:3) = 0d0

      if (nyl == 1) then
        a(:,:,:,1) = 0d0
        a(:,:,:,3) = 0d0
      endif

      if (nzl == 1) a(:,:,:,1:2) = 0d0

c     Return A components

      ax = a(:,:,:,1)
      ay = a(:,:,:,2)
      az = a(:,:,:,3)

c     End program

      contains

ccc     line_int
ccc     #############################################################
cc      subroutine line_int(ilo,ihi,jlo,jhi,bb,aa)
cc
cc      implicit none
cc
ccc     -------------------------------------------------------------
ccc     Performs line integral to find vector potential from magnetic
ccc     field from starting point (ilo,jlo) in any given quadrant
ccc     defined by [ilo,ihi]x[jlo,jhi]. We do NOT assume ilo<ihi or
ccc     jlo<jhi.
ccc     -------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer :: ilo,ihi,jlo,jhi
cc        real(8) :: bb(min(ilo,ihi)-1:max(ilo,ihi)+1
cc     .               ,min(jlo,jhi)-1:max(jlo,jhi)+1,0:nzl+1,1:3)
cc        real(8) :: aa(min(ilo,ihi)-1:max(ilo,ihi)+1
cc     .               ,min(jlo,jhi)-1:max(jlo,jhi)+1,0:nzl+1,1:3)
cc
ccc     Local variables
cc
cc        integer :: i,im,ip,imin,imax,istep
cc     .            ,j,jm,jp,jmin,jmax,jstep
cc
ccc     Begin program
cc
cc        istep = -1
cc        if (ilo < ihi) istep = 1
cc
cc        jstep = -1
cc        if (jlo < jhi) jstep = 1
cc
cc        a(ilo-istep,jlo-jstep,:,3)=0d0
cc
cc        i = ilo-istep
cc
cc        jmin = jlo-min(jstep,0)
cc        jmax = jhi+max(jstep,0)
cc
cc        do j=jmin,jmax,jstep
cc          jm = j-max(jstep,0)
cc          jp = j+min(jstep,0)
cc
cc          dyy = yy(jp)-yy(jm)
cc
cccc          if (spoint) then  !r=0 face
cccc            aa(i,j,:,3)=aa(i,j-1,:,3) + dyy*2.5d-1*(bb(i  ,j-1,:,1)
cccc     .                                             +bb(i  ,j  ,:,1)
cccc     .                                             +bb(i+1,j-1,:,1)
cccc     .                                             +bb(i+1,j  ,:,1))
cccc          else
cc            aa(i,jp,:,3)=aa(i,jm,:,3) + dyy*5d-1*(bb(i,jm,:,1)
cc     .                                           +bb(i,jp,:,1))
cccc          endif
cc        enddo
cc
cc        aa(ilo-istep,:,:,2)=0d0
cc
cc        imin = ilo-min(istep,0)
cc        imax = ihi+max(istep,0)
cc
cc        do i=imin,imax,istep
cc
cc          im = i-max(istep,0)
cc          ip = i+min(istep,0)
cc
cccc          if (.not.glbl) call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cccc
cccc          if (spoint) then
cccc            if (.not.glbl) then
cccc              dxx = g_def%dxh(ig-1)
cccc            else
cccccc              dxx = 0.5*(g_def%xg(i+1)-g_def%xg(i-1))
cccc              dxx = g_def%xg(i)-g_def%xg(i-1)
cccc            endif
cccc
cccc            a(i,:,:,3) = a(i-1,:,:,3) - dxx*b(i,:,:,2)
cccc
cccc            if (isSP(i,1,1,igx,igy,igz)) then
cccc              a(i,:,:,2) = a(i-1,:,:,2) + 5d-1*dxx*b(i,:,:,3)  !Factor of 1/2 due to geometry
cccc            else
cccc              a(i,:,:,2) = a(i-1,:,:,2) +      dxx*b(i,:,:,3)
cccc            endif
cccc
cccc          else
cc
cc            dxx = xx(ip)-xx(im)
cc
cc            aa(ip,:,:,3) = aa(im,:,:,3) - dxx*5d-1*(bb(im,:,:,2)
cc     .                                             +bb(ip,:,:,2))
cc            aa(ip,:,:,2) = aa(im,:,:,2) + dxx*5d-1*(bb(im,:,:,3)
cc     .                                             +bb(ip,:,:,3))
cc
cccc          endif
cc        enddo
cc
cc      !Average from radial faces to nodes in SP coordinate systems
cccc        if (spoint) aa(1:nx+1,:,:,2:3) = 0.5*(aa(1:nx+1,:,:,2:3)
cccc     .                                       +aa(0:nx  ,:,:,2:3))
cc
cc      end subroutine line_int

c     line_int_xz
c     #############################################################
      subroutine line_int_xz(ilo,ihi,klo,khi,bb,aa)

c     -------------------------------------------------------------
c     Performs line integral to find vector potential from magnetic
c     field. Vector potential assumes gauge A_1 = 0. integration is
c     performed from starting point (ilo,jlo) in any given quadrant
c     defined by [ilo,ihi]x[klo,khi]. We do NOT assume ilo<ihi or
c     jlo<jhi.
c     -------------------------------------------------------------

      implicit none

c     Call variables

        integer :: ilo,ihi,klo,khi
        real(8) :: bb(min(ilo,ihi)-1:max(ilo,ihi)+1,0:nyl+1
     .               ,min(klo,khi)-1:max(klo,khi)+1,1:3)
        real(8) :: aa(min(ilo,ihi)-1:max(ilo,ihi)+1,0:nyl+1
     .               ,min(klo,khi)-1:max(klo,khi)+1,1:3)

c     Local variables

        integer :: i,im,ip,imin,imax,istep
     .            ,k,km,kp,kmin,kmax,kstep

c     Begin program

        istep = -1
        if (ilo < ihi) istep = 1

        kstep = -1
        if (klo < khi) kstep = 1

        imin = ilo-min(istep,0)
        imax = ihi+max(istep,0)

        kmin = klo-min(kstep,0)
        kmax = khi+max(kstep,0)

        i = imin

        do k=kmin,kmax,kstep
          km = k-max(kstep,0)
          kp = k+min(kstep,0)

          dzz = zz(kp)-zz(km)

          aa(i,:,kp,2)=aa(i,:,km,2) - dzz*5d-1*(bb(i,:,km,1)
     .                                         +bb(i,:,kp,1))
        enddo

        do i=imin,imax,istep
          im = i-max(istep,0)
          ip = i+min(istep,0)

          dxx = xx(ip)-xx(im)

          aa(ip,:,:,3) = aa(im,:,:,3) - dxx*5d-1*(bb(im,:,:,2)
     .                                           +bb(ip,:,:,2))
          aa(ip,:,:,2) = aa(im,:,:,2) + dxx*5d-1*(bb(im,:,:,3)
     .                                           +bb(ip,:,:,3))
        enddo

      end subroutine line_int_xz

c     #############################################################
      subroutine line_int_yz(jlo,jhi,klo,khi,bb,aa)

      implicit none

c     -------------------------------------------------------------
c     Performs line integral to find vector potential from magnetic
c     field from starting point (ilo,jlo) in any given quadrant
c     defined by [ilo,ihi]x[jlo,jhi]. We do NOT assume ilo<ihi or
c     jlo<jhi.
c     -------------------------------------------------------------

c     Call variables

        integer :: jlo,jhi,klo,khi
        real(8) :: bb(0:nxl+1,min(jlo,jhi)-1:max(jlo,jhi)+1
     .                       ,min(klo,khi)-1:max(klo,khi)+1,1:3)
        real(8) :: aa(0:nxl+1,min(jlo,jhi)-1:max(jlo,jhi)+1
     .                       ,min(klo,khi)-1:max(klo,khi)+1,1:3)

c     Local variables

        integer :: j,jm,jp,jmin,jmax,jstep
     .            ,k,km,kp,kmin,kmax,kstep

c     Begin program

        jstep = -1
        if (jlo < jhi) jstep = 1

        kstep = -1
        if (klo < khi) kstep = 1

        jmin = jlo-min(jstep,0)
        jmax = jhi+max(jstep,0)

        kmin = klo-min(kstep,0)
        kmax = khi+max(kstep,0)

        j = jmin

        do k=kmin,kmax,kstep
          km = k-max(kstep,0)
          kp = k+min(kstep,0)

          dzz = zz(kp)-zz(km)

          aa(:,j,kp,1)=aa(:,j,km,1) + dzz*5d-1*(bb(:,j,km,2)
     .                                         +bb(:,j,kp,2))
        enddo

        do j=jmin,jmax,jstep

          jm = j-max(jstep,0)
          jp = j+min(jstep,0)

          dyy = yy(jp)-yy(jm)

          aa(:,jp,:,1) = aa(:,jm,:,1) - dyy*5d-1*(bb(:,jm,:,3)
     .                                           +bb(:,jp,:,3))
          aa(:,jp,:,3) = aa(:,jm,:,3) + dyy*5d-1*(bb(:,jm,:,1)
     .                                           +bb(:,jp,:,1))
        enddo

      end subroutine line_int_yz

      end subroutine getA_on_mesh

c     splineB
c     #################################################################
      subroutine splineB(bx,by,bz)
c     -----------------------------------------------------------------
c     This routine splines up vector components (bx,by,abz) on a mesh
c     of size (nx,ny,nz) with positions given in xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

c     Local variables

      integer :: alloc_stat,i
      real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)

c     Begin program

      allocate(bcoef(nx,ny,nz,3),stat=alloc_stat)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bx,nx,ny,kx,ky,kz,tx,ty,tz,bcoef(:,:,:,1),work,flg)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,by,nx,ny,kx,ky,kz,tx,ty,tz,bcoef(:,:,:,2),work,flg)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bz,nx,ny,kx,ky,kz,tx,ty,tz,bcoef(:,:,:,3),work,flg)

      end subroutine splineB

c     splineBcar
c     #################################################################
      subroutine splineBcar(bx,by,bz)
c     -----------------------------------------------------------------
c     This routine splines up the Cartesian components of B on a uniform
c     grid of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)

c     Local variables

      integer :: alloc_stat,i

c     Begin program

      allocate(bcarcoef(nx,ny,nz,3),stat=alloc_stat)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bx,nx,ny,kx,ky,kz,tx,ty,tz,bcarcoef(:,:,:,1),work,flg)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,by,nx,ny,kx,ky,kz,tx,ty,tz,bcarcoef(:,:,:,2),work,flg)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bz,nx,ny,kx,ky,kz,tx,ty,tz,bcarcoef(:,:,:,3),work,flg)

      end subroutine splineBcar

c     splineJ
c     #################################################################
      subroutine splineJ(jac)
c     -----------------------------------------------------------------
c     This routine splines up the jacobian on a mesh
c     of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: jac(nx,ny,nz)

c     Local variables

      integer :: alloc_stat,i

c     Begin program

      allocate(jcoef(nx,ny,nz),stat=alloc_stat)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,jac,nx,ny,kx,ky,kz,tx,ty,tz,jcoef,work,flg)

      end subroutine splineJ

c     splineFlds
c     #################################################################
      subroutine splineFlds(flds,fldcoef)
c     -----------------------------------------------------------------
c     This routine splines up a set of fields on a
c     mesh of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: flds(:,:,:,:)
      real(8),optional,INTENT(OUT) :: fldcoef(:,:,:,:)

c     Local variables

      integer :: alloc_stat,i
      real(8),pointer,dimension(:,:,:,:) :: ext_xcoef,lxcoef

c     Begin program

      if (PRESENT(fldcoef)) then
        allocate(ext_xcoef(size(flds,1),size(flds,2)
     .                    ,size(flds,3),size(flds,4)),stat=alloc_stat)
        lxcoef => ext_xcoef
      else
        allocate(xcoef(size(flds,1),size(flds,2)
     .                ,size(flds,3),size(flds,4)),stat=alloc_stat)
        lxcoef => xcoef
      endif

      do i=1,size(flds,4)
        call db3ink(xs,nx,ys,ny,zs,nz,flds(:,:,:,i)
     .             ,nx,ny,kx,ky,kz,tx,ty,tz,lxcoef(:,:,:,i),work,flg)
      enddo

      if (PRESENT(fldcoef)) then
        fldcoef = lxcoef
        deallocate(ext_xcoef)
      endif

      nullify(lxcoef)

c     End program

      end subroutine splineFlds

c     splineX
c     #################################################################
      subroutine splineX(xcar,xcoef_ext)
c     -----------------------------------------------------------------
c     This routine splines up the Cartesian map on a
c     grid of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: xcar(:,:,:,:)
      real(8),optional,INTENT(OUT) :: xcoef_ext(:,:,:,:)

c     Local variables

      call splineFlds(xcar,fldcoef=xcoef_ext)

c     Detect whether physical coord is logical coord (needed for evalXi)

      x_is_log = sqrt(sum((xcar(nx/2,:,:,1)-xcar(nx/2,1,1,1))**2))<1d-10
      y_is_log = sqrt(sum((xcar(:,ny/2,:,2)-xcar(1,ny/2,1,2))**2))<1d-10
      z_is_log = sqrt(sum((xcar(:,:,nz/2,3)-xcar(1,1,nz/2,3))**2))<1d-10

c     End program

      end subroutine splineX

c     evalA
c     #################################################################
      subroutine evalA(x1,x2,x3,ax,ay,az,ierr)
c     -----------------------------------------------------------------
c     This evaluates vector potential components at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierr
      real(8) :: x1,x2,x3,ax,ay,az

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3,ierr=ierr)

      if (ierr /= ORB_OK) return

      select case(ag)
      case(1)
        ax = 0d0
        ay = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .             ,kx,ky,kz,acoefy
     .             ,worke(1+thr_num*dime:(thr_num+1)*dime))
      case(2)
        ax = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .             ,kx,ky,kz,acoefx
     .             ,worke(1+thr_num*dime:(thr_num+1)*dime))
        ay = 0d0
      end select

      az = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,acoefz
     .           ,worke(1+thr_num*dime:(thr_num+1)*dime))

      end subroutine evalA

c     evalCurlA
c     #################################################################
      subroutine evalCurlA(x1,x2,x3,bx,by,bz,ierr,flag,idx,idy,idz)
c     -----------------------------------------------------------------
c     This evaluates B=curl(A) at logical position (x1,x2,x3). Other
c     variables are:
c       * flag: indicates whether two vector potential components
c               are nonzero (0), or two are zero (1,2). Which are
c               zero depends on flag and the gauge chosen (passed to
c               this routine in variable "ag" from module). This
c               capability is needed for the VP integrator.
c       * idx,idy,idz: derivative index in x,y,z respectively, in
c               addition to what is needed to evaluate the curl. This
c               is needed to evaluate the Jacobian matrix of the implicit
c               orbit equations.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierr
      real(8) :: x1,x2,x3,bx,by,bz
      integer,optional :: flag,idx,idy,idz

c     Local variables

      integer :: flg,iidx,iidy,iidz
      real(8) :: ax_y,ax_z,ay_x,ay_z,az_x,az_y

c     Begin program

      if (PRESENT(idx)) then
        iidx = idx
      else
        iidx = 0
      endif

      if (PRESENT(idy)) then
        iidy = idy
      else
        iidy = 0
      endif

      if (PRESENT(idz)) then
        iidz = idz
      else
        iidz = 0
      endif

      if (PRESENT(flag)) then
        flg = flag
      else
        flg = 0
      endif

c     Calculate derivatives

      call chk_pos(x1,x2,x3,ierr=ierr)

      if (ierr /= ORB_OK) return

      select case(ag)
      case(1) !Ax=0

        ax_y = 0d0
        ax_z = 0d0

        select case(flg)
        case(0)
          ay_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          ay_z = db3val(x1,x2,x3,iidx,iidy,1+iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          az_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          az_y = db3val(x1,x2,x3,iidx,1+iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))
        case(1)
          az_x = 0d0
          az_y = 0d0
          ay_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          ay_z = db3val(x1,x2,x3,iidx,iidy,1+iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))
        case(2)
          ay_x = 0d0
          ay_z = 0d0
          az_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          az_y = db3val(x1,x2,x3,iidx,1+iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))
        case default

          write (*,*) 'Flag ',flg,' not implemented in evalCurlA'
          stop

        end select

      case(2) !Ay=0

        ay_x = 0d0
        ay_z = 0d0

        select case(flg)
        case(0)
          ax_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          ax_z = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          az_x = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          az_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))
        case(1)
          az_x = 0d0
          az_y = 0d0
          ax_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          ax_z = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))
        case(2)
          ax_y = 0d0
          ax_z = 0d0
          az_x = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))

          az_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz
     .               ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .               ,work(1+thr_num*dim:(thr_num+1)*dim))
        case default

          write (*,*) 'Flag not implemented in evalCurlA'
          stop

        end select

      end select

c     Find curl

      bx = az_y - ay_z
      by = ax_z - az_x
      bz = ay_x - ax_y

      end subroutine evalCurlA

c     evalB
c     #################################################################
      subroutine evalB(x1,x2,x3,bx,by,bz,ierr)
c     -----------------------------------------------------------------
c     This evaluates magnetic field components at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierr
      real(8) :: x1,x2,x3,bx,by,bz

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3,ierr=ierr)

      if (ierr /= ORB_OK) return

      bx = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,bcoef(:,:,:,1)
     .           ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .           ,work(1+thr_num*dim:(thr_num+1)*dim))
      by = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,bcoef(:,:,:,2)
     .           ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .           ,work(1+thr_num*dim:(thr_num+1)*dim))
      bz = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,bcoef(:,:,:,3)
     .           ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .           ,work(1+thr_num*dim:(thr_num+1)*dim))

      end subroutine evalB

c     evalBcar
c     #################################################################
      subroutine evalBcar(x1,x2,x3,bx,by,bz,ierr)
c     -----------------------------------------------------------------
c     This evaluates Cartesian components of magnetic field at
c     logical position (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierr
      real(8) :: x1,x2,x3,bx,by,bz

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3,ierr=ierr)

      if (ierr /= ORB_OK) return

      bx = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,bcarcoef(:,:,:,1)
     .            ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .            ,work(1+thr_num*dim:(thr_num+1)*dim))

      by = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,bcarcoef(:,:,:,2)
     .            ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .            ,work(1+thr_num*dim:(thr_num+1)*dim))

      bz = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,bcarcoef(:,:,:,3)
     .            ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .            ,work(1+thr_num*dim:(thr_num+1)*dim))

      end subroutine evalBcar

c     getB
c     ##################################################################
      subroutine getB(x1,x2,x3,b1,b2,b3,solen,car,ierr,flag)

c     ------------------------------------------------------------------
c     Wrapper routine to find magnetic field components on specified
c     location in LOGICAL space (x1,x2,x3) using various options:
c      * solen: using curl(A). If so, variable flag indicates which
c               components of A are nonzero (see evalCurlA)>
c      * car: cartesian components
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr
      real(8) :: x1,x2,x3,b1,b2,b3
      logical,INTENT(IN) :: solen,car
      integer,optional :: flag

c     Local variables

c     Begin program

      if (car) then

        !Find Cartesian vector components
        call evalBcar(x1,x2,x3,b1,b2,b3,ierr)

      else
        if (solen) then
          call evalCurlA(x1,x2,x3,b1,b2,b3,ierr,flag=flag)
        else
          call evalB    (x1,x2,x3,b1,b2,b3,ierr)
        endif
      endif

c     End program

      end subroutine getB

c     evalFlds
c     #################################################################
      function evalFlds(x1,x2,x3,spcoef) result(vec)
c     -----------------------------------------------------------------
c     This evaluates cartesian position (x,y,z) at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: spcoef(:,:,:,:)
      real(8) :: x1,x2,x3,vec(size(spcoef,4))

c     Local variables

      integer :: i

c     Begin program

      do i = 1,size(spcoef,4)
        vec(i) = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .                 ,kx,ky,kz,spcoef(:,:,:,i)
     .                 ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .                 ,work(1+thr_num*dim:(thr_num+1)*dim))
      enddo

      end function evalFlds

c     evalX
c     #################################################################
      subroutine evalX(x1,x2,x3,x,y,z,ierr,xcoef_ext,xder,yder,zder)
c     -----------------------------------------------------------------
c     This evaluates cartesian position (x,y,z) at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: x,y,z,x1,x2,x3
      integer,optional :: ierr,xder,yder,zder
      real(8),optional,pointer :: xcoef_ext(:,:,:,:)

c     Local variables

      integer :: ierror,xd,yd,zd
      real(8),pointer :: lxcoef(:,:,:,:)

c     Begin program

      if (PRESENT(xcoef_ext)) then
        lxcoef => xcoef_ext
      else
        lxcoef => xcoef
      endif

      if (PRESENT(xder)) then
        xd = xder
      else
        xd = 0
      endif

      if (PRESENT(xder)) then
        yd = yder
      else
        yd = 0
      endif

      if (PRESENT(xder)) then
        zd = zder
      else
        zd = 0
      endif

      call chk_pos(x1,x2,x3,ierr=ierror)

      if (PRESENT(ierr)) ierr = ierror

      if (ierror /= ORB_OK) return

      x = db3val(x1,x2,x3,xd,yd,zd,tx,ty,tz,nx,ny,nz
     .          ,kx,ky,kz,lxcoef(:,:,:,1)
     .          ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .          ,work(1+thr_num*dim:(thr_num+1)*dim))

      y = db3val(x1,x2,x3,xd,yd,zd,tx,ty,tz,nx,ny,nz
     .          ,kx,ky,kz,lxcoef(:,:,:,2)
     .          ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .          ,work(1+thr_num*dim:(thr_num+1)*dim))

      z = db3val(x1,x2,x3,xd,yd,zd,tx,ty,tz,nx,ny,nz
     .          ,kx,ky,kz,lxcoef(:,:,:,3)
     .          ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .          ,work(1+thr_num*dim:(thr_num+1)*dim))

      nullify(lxcoef)

      end subroutine evalX

c     evalXi
c     #################################################################
      subroutine evalXi(ilevel,x,y,z,x1,x2,x3,ierror,xcoef_ext)
c     -----------------------------------------------------------------
c     This evaluates logical position (x1,x2,x3) at Cartesian position
c     (x,y,z). Requires Newton inversion of map x(xi). On input,
c     (x1,x2,x3) contains initial guess. 
c     -----------------------------------------------------------------

      implicit none            !For safe Fortran

c     Call variables

      integer :: ilevel,ierror
      real(8) :: x,y,z,x1,x2,x3
      real(8),optional,pointer :: xcoef_ext(:,:,:,:)

c     Local variables

      integer,parameter :: maxit=100,size=3,icol=1
      real(8),parameter :: tol=1d-12

      integer :: iter
      real(8) :: JJ(size,size),res(size,icol)
     .          ,dxi(size,icol),error(maxit),damp(maxit)

      logical :: prnt

      real(8),pointer :: lxcoef(:,:,:,:)

c     Begin program

      if (PRESENT(xcoef_ext)) then
        lxcoef => xcoef_ext
      else
        lxcoef => xcoef
      endif

      prnt = (ilevel > 1)

      ierror = ORB_OK

c$$$      x1 = sqrt(x*x+y*y)  !Radius
c$$$      x2 = acos(x/x1)     !Theta
c$$$      x3 = z              !Z

      if (prnt) then
        write (*,*)
        write (*,*) 'evalXi -- x =',x,y,z
        write (*,*) 'evalXi -- xi=',x1,x2,x3
      endif

c     Perform iteration

      do iter=1,maxit

        if (prnt) write (*,*) '>>>>Iteration=',iter

        !Form residual and check convergence
        res(:,1) = formResidual(x1,x2,x3,ierror)

        if (ierror /= ORB_OK) exit

        if (prnt) write (*,*) 'evalXi -- res=',res

        error(iter) = sqrt(sum(res*res))

        if (prnt) write (*,*) 'evalXi -- error=',error(iter)

        if (error(iter) < tol) exit

        !Form Jacobian
        JJ = formJacobian(x1,x2,x3,ierror)

        if (ierror /= ORB_OK) exit

        if (prnt) write (*,*) 'evalXi -- J(1,:)=',JJ(1,:)
        if (prnt) write (*,*) 'evalXi -- J(2,:)=',JJ(2,:)
        if (prnt) write (*,*) 'evalXi -- J(3,:)=',JJ(3,:)

        !Find update dxi =-JJ^-1(res)
        call blockSolve(size,JJ,icol,res,dxi)

        if (prnt) write (*,*) 'evalXi -- dxi=',dxi

        !Linesearch
c$$$        if (iter > 1) then
c$$$          call linesearch(x1,x2,x3,dxi(:,1),error(iter),damp(iter)
c$$$     .                   ,ierror)
c$$$        else
          damp(iter) = 1d0
          x1 = x1 + dxi(1,1)
          x2 = x2 + dxi(2,1)
          x3 = x3 + dxi(3,1)
c$$$        endif

        if (ierror /= ORB_OK) exit

        !Singular point BC
        call chk_pos(x1,x2,x3,ierr=ierror,no_per_bc=.true.)

        if (ierror /= ORB_OK) exit

        if (prnt) write (*,*) 'evalXi -- xi=',x1,x2,x3

      enddo

      if (iter > maxit.or.ierror /= ORB_OK) ierror = ORB_MAP_INV

      if (ilevel == 1 .or. ierror == ORB_MAP_INV) then
        write (*,*)
        write (*,*) 'evalXi -- x =',x,y,z
        write (*,*) 'evalXi -- xi=',x1,x2,x3
        write (*,*) 'evalXi -- domain=',xsmin,xsmax
     .                                 ,ysmin,ysmax
     .                                 ,zsmin,zsmax
        write (*,*) 'evalXi -- Convergence history: '
     .              ,error(1:min(iter,maxit))
        write (*,*) 'evalXi -- Damp parameter history: '
     .              ,damp (1:min(iter,maxit))
        stop
      endif

      contains

c     formResidual
c     ###################################################################
      function formResidual(x1,x2,x3,ierr) result(res)

      implicit none

      integer :: ierr
      real(8) :: x1,x2,x3,res(3)

      real(8) :: x11,x22,x33

      x11 = x1
      x22 = x2
      x33 = x3

      call chk_pos(x11,x22,x33,ierr=ierr)

      if (ierr /= ORB_OK) return

      res(1) = x - db3val(x11,x22,x33,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,lxcoef(:,:,:,1)
     .                   ,worke(1+thr_num*dime:(thr_num+1)*dime))

      res(2) = y - db3val(x11,x22,x33,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,lxcoef(:,:,:,2)
     .                   ,worke(1+thr_num*dime:(thr_num+1)*dime))

      res(3) = z - db3val(x11,x22,x33,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,lxcoef(:,:,:,3)
     .                   ,worke(1+thr_num*dime:(thr_num+1)*dime))

      !Shift residual if physical coord is logical coord
      if (x_is_log) res(1) = res(1) - (x1-x11)
      if (y_is_log) res(2) = res(2) - (x2-x22)
      if (z_is_log) res(3) = res(3) - (x3-x33)

cc      if (prnt) then
cc        write (*,*)
cc        write (*,*) 'evalXi -- sbcnd',sbcnd
cc        write (*,*) 'evalXi -- limits',xsmin,xsmax
cc     .                                ,ysmin,ysmax
cc     .                                ,zsmin,zsmax,tol
cc        write (*,*) 'evalXi -- input log coords',x1 ,x2 ,x3
cc        write (*,*) 'evalXi -- shifted log coords',x11,x22,x33
cc        write (*,*) 'evalXi -- Shifts=',x1-x11,x2-x22,x3-x33
cc        write (*,*) 'evalXi -- res before mod=',res
cc      endif

cc      if (sbcnd(1) == PER) res(1) = mod(res(1),xsmax-xsmin-0.1*tol)
cc      if (sbcnd(3) == PER) res(2) = mod(res(2),ysmax-ysmin-0.1*tol)
cc      if (sbcnd(5) == PER) res(3) = mod(res(3),zsmax-zsmin-0.1*tol)

      end function formResidual

c     formJacobian
c     ###################################################################
      function formJacobian(x1,x2,x3,ierror) result(JJ)

      implicit none

      real(8) :: x1,x2,x3,JJ(3,3)
      integer :: ierror

      real(8) :: x11,x22,x33
      integer :: ieq

      x11 = x1
      x22 = x2
      x33 = x3

      call chk_pos(x11,x22,x33,ierr=ierror)

      if (ierror /= ORB_OK) return

      do ieq=1,3
        JJ(ieq,1) = db3val(x11,x22,x33,1,0,0,tx,ty,tz,nx,ny,nz
     .                    ,kx,ky,kz,lxcoef(:,:,:,ieq)
     .                    ,worke(1+thr_num*dime:(thr_num+1)*dime))
        JJ(ieq,2) = db3val(x11,x22,x33,0,1,0,tx,ty,tz,nx,ny,nz
     .                    ,kx,ky,kz,lxcoef(:,:,:,ieq)
     .                    ,worke(1+thr_num*dime:(thr_num+1)*dime))
        JJ(ieq,3) = db3val(x11,x22,x33,0,0,1,tx,ty,tz,nx,ny,nz
     .                    ,kx,ky,kz,lxcoef(:,:,:,ieq)
     .                    ,worke(1+thr_num*dime:(thr_num+1)*dime))
      enddo

      end function formJacobian

c     linesearch
c     ###############################################################
      subroutine linesearch (x1,x2,x3,ddx,fk,damp,ierror)
      implicit none       !For safe fortran
c     ---------------------------------------------------------------
c     If global is true, uses linesearch backtracking to ensure
c     sufficient reduction in the Newton residual norm.
c     In call:
c       - x1,x2,x3 (input/output): current/next Newton state
c       - ddx (input): Newton update
c       - fk (input): current Newton residual
c       - damp (output): damping parameter
c       - ierror (output): error code
c
c     Taken from C. T. Kelley's "Iterative methods for Linear and
c     Nonlinear Equations", Sec. 8.3.
c     ---------------------------------------------------------------

c     Call variables

      integer :: ierror
      real(8) :: x1,x2,x3,ddx(size),fk,damp

c     Local variables

      integer :: idamp
      real(8) :: dampm,theta,dxnorm,fkp,fm,df,ddf
      real(8) :: b(size),dummy(size),x(size)
      real(8) :: alpha,sigma0,sigma1

c     Begin program

c     Initialize parameters

      damp = 1d0    !Damping

      sigma0 = 0.1  !Lower limit in damping parameter decrease
      sigma1 = 0.5  !Damping parameter decrease for linesearch backtracking

      alpha = 1d-1  !Residual reduction parameter

c     Start quadratic linesearch backtracking

      x = (/x1,x2,x3/)

      do idamp = 1,100

        !Find residual norm |F(xk+dxk)| --> fkp
        dummy = x + damp*ddx
        b = formResidual(dummy(1),dummy(2),dummy(3),ierror)
        fkp = sqrt(sum(b*b))

        !Check residual minimation condition
        if (fkp.lt.(1.-alpha*damp)*fk) then

          exit  !All is well

        else
          !New damping parameter
          dampm = sigma1*damp

          !Evaluate residual norm |F(xk+lambda*dxk)| --> fm
          dummy = x + dampm*ddx
          b = formResidual(dummy(1),dummy(2),dummy(3),ierror)
          fm = sqrt(sum(b*b))

          !Check convergence
          if (fm.lt.(1.-alpha*dampm)*fk) then   !Good enough

            damp = dampm
            exit

          else  !Quadratic minimization in lambda

            !Find out curvature of parabola
            ddf = 2./(damp-dampm)*((fkp-fk)/damp-(fm-fk)/dampm)

            if (ddf > 0) then !Parabola has a minimum

              df = 1./(damp-dampm)*(-dampm/damp*(fkp-fk)
     .                              +damp/dampm*(fm -fk))
              theta = -df/ddf

              damp = fmed(sigma0*damp,sigma1*damp,theta)  !Takes middle value

            else  !Try again

              damp = dampm

            endif

          endif
cc          if (damp < 1d-2) then
cc            damp = 1d0
cc            exit
cc          endif
cc          write (*,*) 'Residuals',fk,fm,fkp
cc          write (*,*) 'Damping parameters',dampm,damp
        endif
      enddo

      x1 = dummy(1)
      x2 = dummy(2)
      x3 = dummy(3)

      if (damp < 1d-8) ierror = ORB_SP_ERR

      if (prnt.and.damp<1d0) write (*,*) 'evalXi -- damping=',damp

c     End program

      end subroutine linesearch

      end subroutine evalXi

c     blockSolve
c     #################################################################
      subroutine blockSolve(size,mat,icol,rhs,x,ret_inv)

c     -----------------------------------------------------------------
c     Solves block systems using a direct solve approach. Requires
c     linking with LAPACK.
c
c     In the call sequence, we have:
c       * size (integer): block size
c       * mat  (real array): block matrix
c       * icol (integer) : number of columns of rhs
c       * rhs  (real array): rhs of equation
c       * x    ( "     "   ): solution (output)
c       * ret_inv (logical): return inverse in mat
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer  :: size,icol

        real(8)  :: mat(size,size),rhs(size,icol),x(size,icol)

        logical,optional :: ret_inv

c     Local variables

        integer  :: ipiv(size),info

        real(8)  :: mat2(size,size)

        external dgesv

c     Begin program

        x = rhs

        if (PRESENT(ret_inv)) then

          call dgesv(size,icol,mat ,size,ipiv,x,size,info) !LAPACK routine

        else

          mat2 = mat            !Avoid overwritting mat

          call dgesv(size,icol,mat2,size,ipiv,x,size,info) !LAPACK routine

        endif

        call LAPACK_error(info)

      end subroutine blockSolve

c     LAPACK_error
c     #################################################################
      subroutine LAPACK_error(info)

c     -----------------------------------------------------------------
c     Error routine for LAPACK calls
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: info

c     Local variables

c     Begin program

        if (info /= 0) then
          if (info < 0) then
            write (*,*) 'Problem in factorization in argument '
     .                 ,-info
            stop
          else
            write (*,*) 'Matrix in blocksolve is singular'
            stop
          endif
        endif

c     End program

      end subroutine LAPACK_error

c     evalJ
c     #################################################################
      function evalJ(x1,x2,x3,ierror) result(jac)
c     -----------------------------------------------------------------
c     This evaluates jacobian at logical position (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierror
      real(8) :: x1,x2,x3,jac

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3,ierr=ierror)

      if (ierror /= ORB_OK) return

      jac = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,jcoef
     .            ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .            ,work(1+thr_num*dim:(thr_num+1)*dim))

      if (jac < 0d0) ierror = ORB_NEG_J !Negative Jacobian

      end function evalJ

c     splineFld
c     #################################################################
      subroutine splineFld(fld,fcoef_ext)
c     -----------------------------------------------------------------
c     This routine splines up a user-input field on a mesh
c     of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: fld(nx,ny,nz)
      real(8),optional,pointer :: fcoef_ext(:,:,:)

c     Local variables

      real(8),pointer,dimension(:,:,:) :: lfcoef => null()

c     Begin program

      if (PRESENT(fcoef_ext)) then
        lfcoef => fcoef_ext
      else
        if (.not.associated(fldcoef)) allocate(fldcoef(nx,ny,nz))
        lfcoef => fldcoef
      endif

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,fld,nx,ny,kx,ky,kz,tx,ty,tz,lfcoef,work,flg)

      nullify(lfcoef)

      end subroutine splineFld

c     evalFld
c     #################################################################
      function evalFld(x1,x2,x3,ierror,fcoef_ext) result(ff)
c     -----------------------------------------------------------------
c     This evaluates a user-input positive-definite field at logical
c     position (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierror
      real(8) :: x1,x2,x3,ff
      real(8),optional,pointer :: fcoef_ext(:,:,:)

c     Local variables

      real(8),pointer :: lfcoef(:,:,:)

c     Begin program

      if (PRESENT(fcoef_ext)) then
        lfcoef => fcoef_ext
      else
        lfcoef => fldcoef
      endif

      call chk_pos(x1,x2,x3,ierr=ierror)

      if (ierror /= ORB_OK) return

      ff = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,lfcoef
     .           ,worke(1+thr_num*dime:(thr_num+1)*dime))
!     .           ,work(1+thr_num*dim:(thr_num+1)*dim))

      nullify(lfcoef)

      end function evalFld

c     killSplines
c     #################################################################
      subroutine killSplines
c     -----------------------------------------------------------------
c     This deallocates memory space for spline routines.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

      integer :: alloc_stat

c     Begin program

      deallocate(worke,work,tx,ty,tz,xs,ys,zs,stat=alloc_stat)
      deallocate(acoefx,acoefy,acoefz,stat=alloc_stat)
      deallocate(jcoef,xcoef,stat=alloc_stat)
      deallocate(bcoef,bcarcoef,stat=alloc_stat)
      deallocate(fldcoef,stat=alloc_stat)

c     End programs

      end subroutine killSplines

      end module spline_field

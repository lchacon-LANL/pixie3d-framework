c module spline_field
c #####################################################################
      module spline_field

        use bc_def

        real(8)  :: db3val
        external :: db3val

        !Private variables
        integer,private :: kx,ky,kz,nx,ny,nz,dim,flg,ag

        real(8),private :: xmin,xmax,ymin,ymax,zmin,zmax

        real(8),private,dimension(:),allocatable    :: tx,ty,tz,work
     .                                                ,xs,ys,zs

        real(8),private,dimension(:,:,:),allocatable:: jcoef
     .                                        ,acoefx,acoefy,acoefz
     .                                        ,bcoefx,bcoefy,bcoefz
     .                                        ,xcoefx,xcoefy,xcoefz
     .                                        ,b2coefx,b2coefz
     .                                        ,b3coefx,b3coefy
     .                                        ,bxcarcoef,bycarcoef
     .                                        ,bzcarcoef

      contains

c     chk_pos
c     ##################################################################
      subroutine chk_pos(x,y,z)

c     ------------------------------------------------------------------
c     Check we are still within LOGICAL domain
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: x,y,z

c     Local Variables

      integer :: ierror

c     Begin program

      ierror = 0

      if(x.gt.xmax) then
        if (bcond(2) == PER) then   !Periodic BC
          x = xmin + mod(x,(xmax-xmin))
        else
          ierror = 1
        endif
      endif

      if(x.lt.xmin) then
        if (bcond(1) == PER) then   !Periodic BC
          x = xmax - mod(xmin-x,xmax-xmin)
        else
          ierror = 1
        endif
      endif

      if(y.gt.ymax) then
        if (bcond(4) == PER) then   !Periodic BC
          y = ymin + mod(y,(ymax-ymin))
        else
          ierror = 1
        endif
      endif

      if(y.lt.ymin) then
        if (bcond(3) == PER) then   !Periodic BC
cc          y = ymax - (ymin-y)
          y = ymax - mod(ymin-y,ymax-ymin)
        else
          ierror = 1
        endif
      endif

      if(z.gt.zmax) then
        if (bcond(6) == PER) then   !Periodic BC
          z = zmin + mod(z,(zmax-zmin))
        else
          ierror = 1
        endif
      endif

      if(z.lt.zmin) then
        if (bcond(5) == PER) then   !Periodic BC
          z = zmax - mod(zmin-z,zmax-zmin)
        else
          ierror = 1
        endif
      endif

      if (ierror /= 0) then
        write (*,*) 'OOPS; out of domain!'
        write (*,*) 'Domain limits=',xmin,xmax,ymin,ymax,zmin,zmax
        write (*,*) 'Current position',x,y,z
cc        write (*,*) 'Old     position',xold,yold,zold
        stop
      endif

c     End program

      end subroutine chk_pos

c     setupSplines
c     #################################################################
      subroutine setupSplines(nnx,nny,nnz,xx,yy,zz,order)
c     -----------------------------------------------------------------
c     This routine splines up vector components (ax,ay,az) on a uniform
c     grid of size (nx,ny,nz) with positions xx,yy,zz. The spline order
c     is set in variable 'order'.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: nnx,nny,nnz,order
      real(8)    :: xx(nnx),yy(nny),zz(nnz)

c     Local variables

      integer :: alloc_stat,i

c     Begin program

c     Initialize private variables

        nx = nnx
        ny = nny
        nz = nnz

        allocate(xs(nx),ys(ny),zs(nz))

c     Initialize spline domain arrays

        xs = xx
        ys = yy
        zs = zz

        if (bcond(1) == PER) then
          xmin = xs(2)
          xmax = xs(nx)
        else
          xmin = xs(1)
          xmax = xs(nx)
        endif

        if (bcond(3) == PER) then
          ymin = ys(2)
          ymax = ys(ny)
        else
          ymin = ys(1)
          ymax = ys(ny)
        endif

        if (bcond(5) == PER) then
          zmin = zs(2)
          zmax = zs(nz)
        else
          zmin = zs(1)
          zmax = zs(nz)
        endif

c     Prepare 3d spline interpolation

        flg = 0 !Let spline routine find interpolation knots
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(work(dim),stat=alloc_stat)
        allocate(tx(nx+kx),stat=alloc_stat)
        allocate(ty(ny+ky),stat=alloc_stat)
        allocate(tz(nz+kz),stat=alloc_stat)

c     End program

      end subroutine setupSplines

c     splineA
c     #################################################################
      subroutine splineA(bx,by,bz,ax,ay,az,a_gauge)
c     -----------------------------------------------------------------
c     This routine splines up vector components (ax,ay,az) on a uniform
c     grid of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer    :: a_gauge
      real(8)    :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
     .             ,ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)

c     Local variables

      integer :: alloc_stat,i

c     Begin program

      ag = a_gauge

c     Get vector potential

      call getA_on_mesh(nx,ny,nz,xs,ys,zs,bx,by,bz,ax,ay,az,a_gauge)

      write (*,*) 'Got vector potential.'
      write (*,*)

c     Spline vector potential

      select case(a_gauge)
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

      end subroutine splineA


c     getA_on_mesh
c     ###############################################################
      subroutine getA_on_mesh(nx,ny,nz,xx,yy,zz,bx,by,bz,ax,ay,az,gauge)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential "a" from CONTRAVARIANT vector
c     field components "b". Employs gauge A_2 = 0d0. In SP systems,
c     it integrates at faces in r, and then averages to nodes.
c     
c     Integrals done by Numerical Recipes p. 128, 4.1.12 -- actually
c     one bit better -- extrapolates quadratically to f(-1) and uses
c     the internal formula, rather than using midpoint for the first
c     point. Thus even the first point has third order accuracy.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,gauge
      real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
      real(8) :: ax(nx,ny,nz),ay(nx,ny,nz),az(nx,ny,nz)
      real(8) :: xx(0:nx-1),yy(0:ny-1),zz(0:nz-1)

c     Local variables

      integer    :: i,j,k,ii,ig,jg,kg,ivar,nxl,nyl,nzl
      real(8)    :: cm,c0,cp,dxx,dyy,dzz,dvol,cov(3)
      real(8)    :: a(0:nx-1,0:ny-1,0:nz-1,3)
     .             ,b(0:nx-1,0:ny-1,0:nz-1,3)

c     Begin program

      nxl = nx-2
      nyl = ny-2
      nzl = nz-2

      b(:,:,:,1) = bx
      b(:,:,:,2) = by
      b(:,:,:,3) = bz

c     Accommodate different gauges

      select case(gauge)
      case(1) !A1=0 code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c       Gauge

cc        a = 0d0
        a(:,:,:,1)=0d0

c       Quadrant-based code (SP ready, parallel-ready)

        !Lower-left quadrant
        call line_int(nxl/2-1,1,nyl/2-1,1
     .               ,b(0:nxl/2,0:nyl/2,0:nzl+1,1:3)
     .               ,a(0:nxl/2,0:nyl/2,0:nzl+1,1:3))

        !Lower-right quadrant
        call line_int(nxl/2+1,nxl,nyl/2-1,1
     .               ,b(nxl/2:nxl+1,0:nyl/2,0:nzl+1,1:3)
     .               ,a(nxl/2:nxl+1,0:nyl/2,0:nzl+1,1:3))

        !Upper-left quadrant
        call line_int(nxl/2-1,1,nyl/2+1,nyl
     .               ,b(0:nxl/2,nyl/2:nyl+1,0:nzl+1,1:3)
     .               ,a(0:nxl/2,nyl/2:nyl+1,0:nzl+1,1:3))

        !Upper-right quadrant
        call line_int(nxl/2+1,nxl,nyl/2+1,nyl
     .               ,b(nxl/2:nxl+1,nyl/2:nyl+1,0:nzl+1,1:3)
     .               ,a(nxl/2:nxl+1,nyl/2:nyl+1,0:nzl+1,1:3))

      case(2) !A2=0 code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c       Gauge

cc        a = 0d0
        a(:,:,:,2)=0d0

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
        write (*,*) 'Gauge not implemented'
        stop
      end select

c     Transform A components to Ay=0 gauge

      ax = a(:,:,:,1)
      ay = a(:,:,:,2)
      az = a(:,:,:,3)

c     End program

      contains

c     line_int
c     #############################################################
      subroutine line_int(ilo,ihi,jlo,jhi,bb,aa)

c     -------------------------------------------------------------
c     Performs line integral to find vector potential from magnetic
c     field from starting point (ilo,jlo) in any given quadrant
c     defined by [ilo,ihi]x[jlo,jhi]. We do NOT assume ilo<ihi or
c     jlo<jhi.
c     -------------------------------------------------------------

c     Call variables

        integer :: ilo,ihi,jlo,jhi
        real(8) :: bb(min(ilo,ihi)-1:max(ilo,ihi)+1
     .               ,min(jlo,jhi)-1:max(jlo,jhi)+1,0:nzl+1,1:3)
        real(8) :: aa(min(ilo,ihi)-1:max(ilo,ihi)+1
     .               ,min(jlo,jhi)-1:max(jlo,jhi)+1,0:nzl+1,1:3)

c     Local variables

        integer :: i,im,ip,imin,imax,istep
     .            ,j,jm,jp,jmin,jmax,jstep

c     Begin program

        istep = -1
        if (ilo < ihi) istep = 1

        jstep = -1
        if (jlo < jhi) jstep = 1

        a(ilo-istep,jlo-jstep,:,3)=0d0

        i = ilo-istep

        jmin = jlo-min(jstep,0)
        jmax = jhi+max(jstep,0)

        do j=jmin,jmax,jstep
          jm = j-max(jstep,0)
          jp = j+min(jstep,0)

          dyy = yy(jp)-yy(jm)

cc          if (spoint) then  !r=0 face
cc            aa(i,j,:,3)=aa(i,j-1,:,3) + dyy*2.5d-1*(bb(i  ,j-1,:,1)
cc     .                                             +bb(i  ,j  ,:,1)
cc     .                                             +bb(i+1,j-1,:,1)
cc     .                                             +bb(i+1,j  ,:,1))
cc          else
            aa(i,jp,:,3)=aa(i,jm,:,3) + dyy*5d-1*(bb(i,jm,:,1)
     .                                           +bb(i,jp,:,1))
cc          endif
        enddo

        aa(ilo-istep,:,:,2)=0d0

        imin = ilo-min(istep,0)
        imax = ihi+max(istep,0)

        do i=imin,imax,istep

          im = i-max(istep,0)
          ip = i+min(istep,0)

cc          if (.not.glbl) call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc          if (spoint) then
cc            if (.not.glbl) then
cc              dxx = grid_params%dxh(ig-1)
cc            else
cccc              dxx = 0.5*(grid_params%xg(i+1)-grid_params%xg(i-1))
cc              dxx = grid_params%xg(i)-grid_params%xg(i-1)
cc            endif
cc
cc            a(i,:,:,3) = a(i-1,:,:,3) - dxx*b(i,:,:,2)
cc
cc            if (isSP(i,1,1,igx,igy,igz)) then
cc              a(i,:,:,2) = a(i-1,:,:,2) + 5d-1*dxx*b(i,:,:,3)  !Factor of 1/2 due to geometry
cc            else
cc              a(i,:,:,2) = a(i-1,:,:,2) +      dxx*b(i,:,:,3)
cc            endif
cc
cc          else

            dxx = xx(ip)-xx(im)

            aa(ip,:,:,3) = aa(im,:,:,3) - dxx*5d-1*(bb(im,:,:,2)
     .                                             +bb(ip,:,:,2))
            aa(ip,:,:,2) = aa(im,:,:,2) + dxx*5d-1*(bb(im,:,:,3)
     .                                             +bb(ip,:,:,3))

cc          endif
        enddo

      !Average from radial faces to nodes in SP coordinate systems
cc        if (spoint) aa(1:nx+1,:,:,2:3) = 0.5*(aa(1:nx+1,:,:,2:3)
cc     .                                       +aa(0:nx  ,:,:,2:3))

      end subroutine line_int

c     #############################################################
      subroutine line_int_yz(jlo,jhi,klo,khi,bb,aa)

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

        a(:,jlo-jstep,klo-kstep,3)=0d0

        j = jlo-jstep

        kmin = klo-min(kstep,0)
        kmax = khi+max(kstep,0)

        do k=kmin,kmax,kstep
          km = k-max(kstep,0)
          kp = k+min(kstep,0)

          dzz = zz(kp)-zz(km)

          aa(:,j,kp,1)=aa(:,j,km,1) + dzz*5d-1*(bb(:,j,km,2)
     .                                         +bb(:,j,kp,2))
        enddo

        aa(:,jlo-jstep,:,2)=0d0

        jmin = jlo-min(jstep,0)
        jmax = jhi+max(jstep,0)

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


ccc getA_on_mesh
ccc#######################################################################
cc      subroutine getA_on_mesh(nx,ny,nz,xx,yy,zz,bx,by,bz,ax,az,order)
cc
cc      implicit none
cc
ccc-----------------------------------------------------------------------
ccc     Variable-order inverse-curl operation, performed by choosing 
ccc     gauge Ay=0d0 and integrating resulting line integrals. The
ccc     order of integration is determined by order:
ccc       * order = 2: Integrals done by trapezoidal rule.
ccc       * order = 3: Integrals done by Numerical Recipes p. 128, 4.1.12 
ccc                    -- actually one bit better -- extrapolates quadratically
ccc                    to f(-1) and uses the internal formula, rather than
ccc                    using midpoint for the first point. Thus even the
ccc                    first point has third order accuracy.
ccc-----------------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer :: nx,ny,nz,order
cc      real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
cc      real(8) :: ax(nx,ny,nz),az(nx,ny,nz)
cc      real(8) :: xx(nx),yy(ny),zz(nz)
cc
ccc     Local variables
cc
cc      integer :: i,j,k
cc      real(8) :: tm,t0,cm,c0,cp,halfx,halfy,halfz
cc      real(8) :: dx,dy,dz
cc
ccc     Begin program
cc
cc      ax(1,:,1) = 0d0
cc
cc      select case(order)
cc      case(1,2)
cc
cc        do k=2,nz
cc          halfz = 0.5*(zz(k)-zz(k-1))
cc          ax(:,1,k)=ax(:,1,k-1)+halfz*(by(:,1,k-1)+by(:,1,k))
cc        enddo
cc
cc        az(:,1,:) = 0d0
cc
cc        do j=2,ny
cc          halfy = 0.5*(yy(j)-yy(j-1))
cc          ax(:,j,:)=ax(:,j-1,:)-halfy*(bz(:,j-1,:)+bz(:,j,:))
cc          az(:,j,:)=az(:,j-1,:)+halfy*(bx(:,j-1,:)+bx(:,j,:))
cc        enddo
cc
cc      case(3)
cc
cc        tm=0.5
cc        t0=0.5
cc        cm=-0.0833333333333
cc        c0= 0.6666666666667
cc        cp= 0.4166666666667
cc
cc        do k=2,nz
cc          dz = (zz(k)-zz(k-1))
cc          if(k.eq.2) then
cccc            ax(:,1,k)=ax(:,1,k-1)+ dz*(tm*by(:,1,k-1)+t0*by(:,1,k))
cc            ax(:,1,k)=ax(:,1,k-1)+ dz*(cp*by(:,1,k-1)
cc     .                                +c0*by(:,1,k  )
cc     .                                +cm*by(:,1,k+1))
cc          else
cc            ax(:,1,k)=ax(:,1,k-1)+ dz*(cm*by(:,1,k-2)
cc     .                                +c0*by(:,1,k-1)
cc     .                                +cp*by(:,1,k  ))
cc          endif
cc        enddo
cc
cc        az(:,1,:) = 0d0
cc
cc        do j=2,ny
cc          dy = (yy(j)-yy(j-1))
cc          if (j == 2) then
cccc            ax(:,j,:)=ax(:,j-1,:)-dy*(tm*bz(:,j-1,:)+t0*bz(:,j,:))
cccc            az(:,j,:)=az(:,j-1,:)+dy*(tm*bx(:,j-1,:)+t0*bx(:,j,:))
cc            ax(:,j,:)=ax(:,j-1,:)- dy*(cp*bz(:,j-1,:)
cc     .                                +c0*bz(:,j  ,:)
cc     .                                +cm*bz(:,j+1,:))
cc            az(:,j,:)=az(:,j-1,:)+ dy*(cp*bx(:,j-1,:)
cc     .                                +c0*bx(:,j  ,:)
cc     .                                +cm*bx(:,j+1,:))
cc          else
cc            ax(:,j,:)=ax(:,j-1,:)- dy*(cm*bz(:,j-2,:)
cc     .                                +c0*bz(:,j-1,:)
cc     .                                +cp*bz(:,j  ,:))
cc            az(:,j,:)=az(:,j-1,:)+ dy*(cm*bx(:,j-2,:)
cc     .                                +c0*bx(:,j-1,:)
cc     .                                +cp*bx(:,j  ,:))
cc
cc          endif
cc        enddo
cc
cc      case default
cc
cc        write (*,*) 'Order not implemented yet in getA'
cc        stop
cc
cc      end select
cc
cc      end subroutine getA_on_mesh

c     splineB
c     #################################################################
      subroutine splineB(bx,by,bz)
c     -----------------------------------------------------------------
c     This routine splines up vector components (ax,ay,az) on a uniform
c     grid of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

c     Local variables

      integer :: alloc_stat,i,j,k
cc      real(8) :: b1(nx,ny,nz),b2(nx,ny,nz),bx1,by1,bz1
      real(8) :: bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)

c     Begin program

cc      !B2 magnetic field
cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc            call evalCurlA(xs(i),ys(j),zs(k),bx1,by1,bz1,flag=1)
cc            b1(i,j,k)=bx1
cc            b2(i,j,k)=bz1
cc          enddo
cc        enddo
cc      enddo
cc
cc      allocate(b2coefx(nx,ny,nz),stat=alloc_stat)
cc      call db3ink(xs,nx,ys,ny,zs,nz
cc     .           ,b1,nx,ny,kx,ky,kz,tx,ty,tz,b2coefx,work,flg)
cc
cc      allocate(b2coefz(nx,ny,nz),stat=alloc_stat)
cc      call db3ink(xs,nx,ys,ny,zs,nz
cc     .           ,b2,nx,ny,kx,ky,kz,tx,ty,tz,b2coefz,work,flg)
cc
cc      !B3 magnetic field
cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc            call evalCurlA(xs(i),ys(j),zs(k),bx1,by1,bz1,flag=2)
cc            b1(i,j,k)=bx1
cc            b2(i,j,k)=by1
cc          enddo
cc        enddo
cc      enddo
cc
cc      allocate(b3coefx(nx,ny,nz),stat=alloc_stat)
cc      call db3ink(xs,nx,ys,ny,zs,nz
cc     .           ,b1,nx,ny,kx,ky,kz,tx,ty,tz,b3coefx,work,flg)
cc
cc      allocate(b3coefy(nx,ny,nz),stat=alloc_stat)
cc      call db3ink(xs,nx,ys,ny,zs,nz
cc     .           ,b2,nx,ny,kx,ky,kz,tx,ty,tz,b3coefy,work,flg)

      allocate(bcoefx(nx,ny,nz),stat=alloc_stat)
      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bx,nx,ny,kx,ky,kz,tx,ty,tz,bcoefx,work,flg)

      allocate(bcoefy(nx,ny,nz),stat=alloc_stat)
      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,by,nx,ny,kx,ky,kz,tx,ty,tz,bcoefy,work,flg)

      allocate(bcoefz(nx,ny,nz),stat=alloc_stat)
      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bz,nx,ny,kx,ky,kz,tx,ty,tz,bcoefz,work,flg)

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

      allocate(bxcarcoef(nx,ny,nz)
     .        ,bycarcoef(nx,ny,nz)
     .        ,bzcarcoef(nx,ny,nz),stat=alloc_stat)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bx,nx,ny,kx,ky,kz,tx,ty,tz,bxcarcoef,work,flg)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,by,nx,ny,kx,ky,kz,tx,ty,tz,bycarcoef,work,flg)

      call db3ink(xs,nx,ys,ny,zs,nz
     .           ,bz,nx,ny,kx,ky,kz,tx,ty,tz,bzcarcoef,work,flg)

      end subroutine splineBcar

c     splineJ
c     #################################################################
      subroutine splineJ(jac)
c     -----------------------------------------------------------------
c     This routine splines up the jacobian on a uniform
c     grid of size (nx,ny,nz) with positions xs,ys,zs.
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

c     splineX
c     #################################################################
      subroutine splineX(xcar)
c     -----------------------------------------------------------------
c     This routine splines up the Cartesian map on a uniform
c     grid of size (nx,ny,nz) with positions xs,ys,zs.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8)    :: xcar(nx,ny,nz,3)

c     Local variables

      integer :: alloc_stat

c     Begin program

      allocate(xcoefx(nx,ny,nz),stat=alloc_stat)
      allocate(xcoefy(nx,ny,nz),stat=alloc_stat)
      allocate(xcoefz(nx,ny,nz),stat=alloc_stat)

      call db3ink(xs,nx,ys,ny,zs,nz,xcar(:,:,:,1)
     .           ,nx,ny,kx,ky,kz,tx,ty,tz,xcoefx,work,flg)
      call db3ink(xs,nx,ys,ny,zs,nz,xcar(:,:,:,2)
     .           ,nx,ny,kx,ky,kz,tx,ty,tz,xcoefy,work,flg)
      call db3ink(xs,nx,ys,ny,zs,nz,xcar(:,:,:,3)
     .           ,nx,ny,kx,ky,kz,tx,ty,tz,xcoefz,work,flg)

c     End program

      end subroutine splineX

c     evalA
c     #################################################################
      subroutine evalA(x1,x2,x3,ax,ay,az)
c     -----------------------------------------------------------------
c     This evaluates vector potential components at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: x1,x2,x3,ax,ay,az

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3)

      select case(ag)
      case(1)
        ax = 0d0
        ay = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .             ,kx,ky,kz,acoefy,work)
      case(2)
        ax = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .             ,kx,ky,kz,acoefx,work)
        ay = 0d0
      end select

      az = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,acoefz,work)

      end subroutine evalA

c     evalCurlA
c     #################################################################
      subroutine evalCurlA(x1,x2,x3,bx,by,bz,flag,idx,idy,idz)
c     -----------------------------------------------------------------
c     This evaluates curl(A) at logical position (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

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

      call chk_pos(x1,x2,x3)

      if (PRESENT(flag)) then
        flg = flag
      else
        flg = 0
      endif

c     Calculate derivatives

      select case(ag)
      case(1) !Ax=0

        ax_y = 0d0
        ax_z = 0d0

        select case(flg)
        case(0)
          ay_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy,work)

          ay_z = db3val(x1,x2,x3,iidx,iidy,1+iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy,work)

          az_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)

          az_y = db3val(x1,x2,x3,iidx,1+iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)
        case(1)
          az_x = 0d0
          az_y = 0d0
          ay_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy,work)

          ay_z = db3val(x1,x2,x3,iidx,iidy,1+iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefy,work)
        case(2)
          ay_x = 0d0
          ay_z = 0d0
          az_x = db3val(x1,x2,x3,1+iidx,iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)

          az_y = db3val(x1,x2,x3,iidx,1+iidy,iidz,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)
        case default

          write (*,*) 'Flag not implemented in evalCurlA'
          stop

        end select

      case(2) !Ay=0

        ay_x = 0d0
        ay_z = 0d0

        select case(flg)
        case(0)
          ax_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx,work)

          ax_z = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx,work)

          az_x = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)

          az_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)
        case(1)
          az_x = 0d0
          az_y = 0d0
          ax_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx,work)

          ax_z = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefx,work)
        case(2)
          ax_y = 0d0
          ax_z = 0d0
          az_x = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)

          az_y = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .               ,kx,ky,kz,acoefz,work)
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
      subroutine evalB(x1,x2,x3,bx,by,bz)
c     -----------------------------------------------------------------
c     This evaluates vector potential components at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: x1,x2,x3,bx,by,bz

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3)

      bx = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,bcoefx,work)
      by = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,bcoefy,work)
      bz = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .           ,kx,ky,kz,bcoefz,work)

      end subroutine evalB

c     evalBcar
c     #################################################################
      subroutine evalBcar(x1,x2,x3,bx,by,bz)
c     -----------------------------------------------------------------
c     This evaluates Cartesian components of magnetic field at
c     logical position (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: x1,x2,x3,bx,by,bz

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3)

      bx = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,bxcarcoef,work)

      by = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,bycarcoef,work)

      bz = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,bzcarcoef,work)

      end subroutine evalBcar

c     getB
c     ##################################################################
      subroutine getB(x1,x2,x3,b1,b2,b3,solen,car)

c     ------------------------------------------------------------------
c     Finds magnetic field components on specified location in LOGICAL
c     space (x1,x2,x3)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: x1,x2,x3,b1,b2,b3
      logical :: solen,car

c     Local variables

      real(8) :: bx1,by1,bz1,bx2,by2,bz2

c     Begin program

      if (car) then

        !Find Cartesian vector components
        call evalBcar(x1,x2,x3,b1,b2,b3)

      else
        if (solen) then
          call evalCurlA(x1,x2,x3,b1,b2,b3)
        else
          call evalB    (x1,x2,x3,b1,b2,b3)
        endif
      endif

c     End program

      end subroutine getB

c     evalX
c     #################################################################
      subroutine evalX(x1,x2,x3,x,y,z)
c     -----------------------------------------------------------------
c     This evaluates cartesian position (x,y,z) at logical position
c     (x1,x2,x3).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      real(8) :: x,y,z,x1,x2,x3

c     Local variables

c     Begin program

      call chk_pos(x1,x2,x3)

      x = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .          ,kx,ky,kz,xcoefx,work)

      y = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .          ,kx,ky,kz,xcoefy,work)

      z = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .          ,kx,ky,kz,xcoefz,work)

      end subroutine evalX

c     evalXi
c     #################################################################
      subroutine evalXi(x,y,z,x1,x2,x3,ierror)
c     -----------------------------------------------------------------
c     This evaluates logical position (x1,x2,x3) at Cartesian position
c     (x,y,z). Requires Newton inversion of map x(xi). On input,
c     (x1,x2,x3) contains initial guess. 
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      integer :: ierror
      real(8) :: x,y,z,x1,x2,x3

c     Local variables

      integer,parameter :: maxit=100,size=3,icol=1

      integer :: iter
      real(8) :: JJ(size,size),res(size,icol)
     .          ,dxi(size,icol),error(maxit)

c     Begin program

cc      if (ilevel > 0) write (*,*) 'evalXi -- xi=',x1,x2,x3

      do iter=1,maxit

        !Form residual and check convergence
        res(:,1) = formResidual(x1,x2,x3)

        error(iter) = sqrt(sum(res**2))

cc        if (ilevel > 0) write (*,*) 'evalXi -- error=',error(iter)

        if (error(iter) < 1d-10) exit

        !Form Jacobian
        JJ = formJacobian(x1,x2,x3)

        !Find update dxi = JJ^-1(res)
        call blockSolve(size,JJ,icol,res,dxi)

cc        if (ilevel > 0) write (*,*) 'evalXi -- dxi=',dxi

        !Update solution
        x1 = x1 + dxi(1,1)
        x2 = x2 + dxi(2,1)
        x3 = x3 + dxi(3,1)

        if (x1 < 0d0) then  !Apply regularity
          x1 = -x1
          x2 = x2 + acos(-1d0)
        endif

        call chk_pos(x1,x2,x3)

cc        if (ilevel > 0) write (*,*) 'evalXi -- xi=',x1,x2,x3

      enddo

      if (iter > maxit) then
        ierror = 1
        write (*,*) 'No convergence in Newton'
        write (*,*) 'Convergence history: ',error
        return
      endif

      contains

c     formResidual
c     ###################################################################
      function formResidual(x1,x2,x3) result(res)

      implicit none

      real(8) :: x1,x2,x3,res(3)

      res(1) = x - db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,xcoefx,work)

      res(2) = y - db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,xcoefy,work)

      res(3) = z - db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .                   ,kx,ky,kz,xcoefz,work)

      end function formResidual

c     formJacobian
c     ###################################################################
      function formJacobian(x1,x2,x3) result(JJ)

      implicit none

      real(8) :: x1,x2,x3,JJ(3,3)

      JJ(1,1) = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefx,work)
      JJ(1,2) = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefx,work)
      JJ(1,3) = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefx,work)

      JJ(2,1) = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefy,work)
      JJ(2,2) = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefy,work)
      JJ(2,3) = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefy,work)

      JJ(3,1) = db3val(x1,x2,x3,1,0,0,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefz,work)
      JJ(3,2) = db3val(x1,x2,x3,0,1,0,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefz,work)
      JJ(3,3) = db3val(x1,x2,x3,0,0,1,tx,ty,tz,nx,ny,nz
     .                ,kx,ky,kz,xcoefz,work)

      end function formJacobian

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

      call chk_pos(x1,x2,x3)

      jac = db3val(x1,x2,x3,0,0,0,tx,ty,tz,nx,ny,nz
     .            ,kx,ky,kz,jcoef,work)

      ierror = 0
      if (jac < 0d0) ierror = 1

      end function evalJ

c     killSplines
c     #################################################################
      subroutine killSplines
c     -----------------------------------------------------------------
c     This evaluates curl(A) at (x,y,z).
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

      integer :: alloc_stat

c     Begin program

      deallocate(work,tx,ty,tz,xs,ys,zs,stat=alloc_stat)
      deallocate(acoefx,acoefy,acoefz,stat=alloc_stat)
      deallocate(jcoef,xcoefx,xcoefy,xcoefz,stat=alloc_stat)
      deallocate(b2coefx,b2coefz,b3coefx,b3coefy,stat=alloc_stat)
      deallocate(bcoefx,bcoefy,bcoefz,stat=alloc_stat)

c     End programs

      end subroutine killSplines

      end module spline_field

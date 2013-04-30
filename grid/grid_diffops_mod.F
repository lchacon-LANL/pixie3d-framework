c grid_vertex_ops
c ######################################################################
      module grid_vertex_ops

        use grid_create

      contains

c     grad_vrtx
c     ##################################################################
      function grad_vrtx(i,j,k,igx,igy,igz,phi,vrtx) result(grd)

c     ------------------------------------------------------------------
c     Calculates gradient of phi from nearest neighbors.
c     If vrtx is true, it finds grad at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds grad at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: phi(2,2,2),grd(3)
      logical :: vrtx

c     Local variables

      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0,nx,ny,nz
      real(8) :: idhx,idhy,idhz,fp,fm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      nx = grid_params%nxgl(igx)
      ny = grid_params%nygl(igy)
      nz = grid_params%nzgl(igz)

      if (vrtx) then
        idhx = 1./grid_params%dx(ig)
        idhy = 1./grid_params%dy(jg)
        idhz = 1./grid_params%dz(kg)
      else
        idhx = 1./grid_params%dxh(ig)
        idhy = 1./grid_params%dyh(jg)
        idhz = 1./grid_params%dzh(kg)
      endif

      ip = 2
      i0 = 1
      jp = 2
      j0 = 1
      kp = 2
      k0 = 1

      !X component
      fp = 0.25*(phi(ip,j0,k0)+phi(ip,jp,k0)
     .          +phi(ip,j0,kp)+phi(ip,jp,kp))
      fm = 0.25*(phi(i0,j0,k0)+phi(i0,jp,k0)
     .          +phi(i0,j0,kp)+phi(i0,jp,kp))

      grd(1) = (fp-fm)*idhx

      !Y component
      fp = 0.25*(phi(i0,jp,k0)+phi(ip,jp,k0)
     .          +phi(i0,jp,kp)+phi(ip,jp,kp))
      fm = 0.25*(phi(i0,j0,k0)+phi(ip,j0,k0)
     .          +phi(i0,j0,kp)+phi(ip,j0,kp))

      grd(2) = (fp-fm)*idhy

      !Z component
      fp = 0.25*(phi(ip,j0,kp)+phi(i0,j0,kp)
     .          +phi(ip,jp,kp)+phi(i0,jp,kp))
      fm = 0.25*(phi(ip,j0,k0)+phi(i0,j0,k0)
     .          +phi(ip,jp,k0)+phi(i0,jp,k0))

      grd(3) = (fp-fm)*idhz

c     End program

      end function grad_vrtx

c     curl_vrtx
c     ##################################################################
      function curl_vrtx(i,j,k,igx,igy,igz,vv,vrtx) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curl of vv (covariant) from nearest neighbors. Returns
c     contravariant components.
c
c     If vrtx is true, it finds curl at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds curl at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: vv(2,2,2,3),cnv(3)
      logical :: vrtx

c     Local variables

      real(8) :: gvx(3),gvy(3),gvz(3)

c     Begin program

      gvx = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,1),vrtx)
      gvy = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,2),vrtx)
      gvz = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,3),vrtx)

      cnv(1) = gvz(2)-gvy(3)
      cnv(2) = gvx(3)-gvz(1)
      cnv(3) = gvy(1)-gvx(2)

c     End program

      end function curl_vrtx

c     div_vrtx
c     ##################################################################
      function div_vrtx(i,j,k,igx,igy,igz,v,vrtx) result(div)

c     ------------------------------------------------------------------
c     Calculates div of cnv vector defined at nearest neighbors.
c     If vrtx is true, it finds div at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds div at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: v(2,2,2,3),div
      logical :: vrtx

c     Local variables

      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0
      real(8) :: dhx,dhy,dhz,fp,fm,jac

c     Begin program

      if (isSP(i+1,j,k,igx,igy,igz).and.vrtx) then
        div = 0d0
        return
      endif

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      if (vrtx) then
        dhx = grid_params%dx(ig)
        dhy = grid_params%dy(jg)
        dhz = grid_params%dz(kg)
      else
        dhx = grid_params%dxh(ig)
        dhy = grid_params%dyh(jg)
        dhz = grid_params%dzh(kg)
      endif

      ip = 2
      i0 = 1
      jp = 2
      j0 = 1
      kp = 2
      k0 = 1

      !X component
      fp = 0.25*(v(ip,j0,k0,1)+v(ip,jp,k0,1)
     .          +v(ip,j0,kp,1)+v(ip,jp,kp,1))
      fm = 0.25*(v(i0,j0,k0,1)+v(i0,jp,k0,1)
     .          +v(i0,j0,kp,1)+v(i0,jp,kp,1))

      div = (fp-fm)/dhx

      !Y component
      fp = 0.25*(v(i0,jp,k0,2)+v(ip,jp,k0,2)
     .          +v(i0,jp,kp,2)+v(ip,jp,kp,2))
      fm = 0.25*(v(i0,j0,k0,2)+v(ip,j0,k0,2)
     .          +v(i0,j0,kp,2)+v(ip,j0,kp,2))

      div = div + (fp-fm)/dhy

      !Z component
      fp = 0.25*(v(ip,j0,kp,3)+v(i0,j0,kp,3)
     .          +v(ip,jp,kp,3)+v(i0,jp,kp,3))
      fm = 0.25*(v(ip,j0,k0,3)+v(i0,j0,k0,3)
     .          +v(ip,jp,k0,3)+v(i0,jp,k0,3))

      div = div + (fp-fm)/dhz

      !Jacobian factor
      if (vrtx) then
        jac = gmetric%grid(igx)%jac_v(i,j,k)
      else
        jac = gmetric%grid(igx)%jac(i,j,k)
      endif

      div = div/jac

c     End program

      end function div_vrtx

c     graddiv_vrtx
c     ###############################################################
      function graddiv_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
     .         result(cov)

c     ---------------------------------------------------------------
c     Calculates grad(div(A)) at cell centers in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cov(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: div(2,2,2)

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

c     Divergence at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2

            div(il,jl,kl) = div_vrtx(ii,jj,kk,igx,igy,igz
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,:)
     .                         ,.true.)
          enddo
        enddo
      enddo

c     Gradient at cell centers

      cov = grad_vrtx(i,j,k,igx,igy,igz,div,.false.)

      if (vol_wgt) cov=cov*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function graddiv_vrtx

c     lap_vrtx
c     ###############################################################
      function lap_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
     .         result(lap)

c     ---------------------------------------------------------------
c     Calculates lap(A) at cell centers in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1),lap

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: cov(3),cnv(2,2,2,3)

      logical    :: vol_wgt,enfbc

c     Begin program

      if (PRESENT(vol)) then
        vol_wgt = vol
      else
        vol_wgt = .false.
      endif

c     Gradient at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2

            cov = grad_vrtx(ii,jj,kk,igx,igy,igz
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1)
     .                     ,.true.)

            !Transform to contravariant at vertex
            cnv(il,jl,kl,:) =
     .           matmul(gmetric%grid(igx)%gsup_v(ii,jj,kk,:,:),cov)
          enddo
        enddo
      enddo

c     Divergence at cell centers

      lap = div_vrtx(i,j,k,igx,igy,igz,cnv,.false.)

      if (vol_wgt) lap=lap*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_vrtx

c$$$c     curlcurl_vrtx_ndiff
c$$$c     ###############################################################
c$$$      function curlcurl_vrtx_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
c$$$     .         result(cnv)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates curl(curl(A)) in general non-orthogonal
c$$$c     coordinates, preserving the SPD property. The vector A is
c$$$c     covariant, and returns contravariant components.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)
c$$$
c$$$      logical,optional,intent(IN) :: vol
c$$$
c$$$c     Local variables
c$$$
c$$$      integer    :: il,jl,kl,ii,jj,kk
c$$$
c$$$      real(8)    :: crl(2,2,2,3),cov(3)
c$$$
c$$$      logical    :: vol_wgt
c$$$
c$$$c     Begin program
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$c     Curl at vertices
c$$$
c$$$      do kl=1,2
c$$$        do jl=1,2
c$$$          do il=1,2
c$$$            ii = i+il-2
c$$$            jj = j+jl-2
c$$$            kk = k+kl-2
c$$$            cnv = curl_vrtx(ii,jj,kk,igx,igy,igz
c$$$     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1,:)
c$$$     .                     ,.true.)
c$$$
c$$$            call transformFromCurvToCurv(ii,jj,kk,igx,igy,igz
c$$$     .                                  ,cov(1),cov(2),cov(3)
c$$$     .                                  ,cnv(1),cnv(2),cnv(3)
c$$$     .                                  ,.false.
c$$$     .                                  ,hex=1,hey=1,hez=1)
c$$$
c$$$            crl(il,jl,kl,:) = cov
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     Curl at cell center
c$$$
c$$$      cnv = curl_vrtx(i,j,k,igx,igy,igz,crl,.false.)
c$$$
c$$$      if (vol_wgt) cnv=cnv*gmetric%grid(igx)%dvol(i,j,k)
c$$$
c$$$c     End program
c$$$
c$$$      end function curlcurl_vrtx_ndiff

c     curlcurl_vrtx
c     ###############################################################
      function curlcurl_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,aa,diff,vol)
     $         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(D*curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)
      real(8),optional,target :: diff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: il,jl,kl,ii,jj,kk

      real(8) :: crl(2,2,2,3),cov(3),dff

      logical :: vol_wgt,have_coef

      real(8),pointer,dimension(:,:,:) :: coef => null()

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      have_coef = PRESENT(diff)
      if (have_coef) coef => diff

c     Curl at vertices

      dff = 1d0

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2
            cnv = curl_vrtx(ii,jj,kk,igx,igy,igz
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1,:)
     .                     ,.true.)

            call transformFromCurvToCurv(ii,jj,kk,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.
     .                                  ,hex=1,hey=1,hez=1)

            if (have_coef) dff=0.125*sum(coef(ii:ii+1,jj:jj+1,kk:kk+1))

            crl(il,jl,kl,:) = cov*dff
          enddo
        enddo
      enddo

      nullify(coef)

c     Curl at cell center

      cnv = curl_vrtx(i,j,k,igx,igy,igz,crl,.false.)

      if (vol_wgt) cnv=cnv*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function curlcurl_vrtx

      end module grid_vertex_ops

c grid_nabla_ops
c ######################################################################
      module grid_nabla_ops

        use grid_create

        INTERFACE grad
          module procedure grad_ijk,grad_mesh
        end INTERFACE

        INTERFACE curl
          module procedure curl_ijk,curl_mesh
        end INTERFACE

        INTERFACE div
          module procedure div_ijk,div_mesh
        end INTERFACE

      contains

c     grad_ijk
c     ###############################################################
      function grad_ijk(i,j,k,nx,ny,nz,igx,igy,igz,arr,hex,hey,hez
     .                 ,order) result(cov)

c     ---------------------------------------------------------------
c     Calculates grad(A) in general non-orthogonal coordinates.
c     The vector grad(A) is covariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1),cov(3)

      integer,optional :: order,hex,hey,hez

c     Local variables

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg,hx,hy,hz
      real(8) :: dhx,dhy,dhz,dx1,dx2,dy1,dy2,dz1,dz2
      real(8) :: gx,gy,gz

c     Begin program

      hx = 0; hy = 0; hz = 0

      if (PRESENT(hex)) hx = hex
      if (PRESENT(hey)) hy = hey
      if (PRESENT(hez)) hz = hez

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c     Central-difference formulas

cc      dhx = 2.*grid_params%dxh(ig)
cc      dhy = 2.*grid_params%dyh(jg)
cc      dhz = 2.*grid_params%dzh(kg)
cccc      dhx = grid_params%dx(ig-1)+grid_params%dx(ig)
cccc      dhy = grid_params%dy(jg-1)+grid_params%dy(jg)
cccc      dhz = grid_params%dz(kg-1)+grid_params%dz(kg)
cc
ccc     Exceptions
cc
cc      if (hx == 1) then
cc        im = i
cc        dhx = grid_params%dx(ig)
cc      elseif (hx == -1) then
cc        ip = i
cc        dhx = grid_params%dx(ig-1)
cc      endif
cc
cc      if (hy == 1) then
cc        jm = j
cc        dhy = grid_params%dy(jg)
cc      elseif (hy == -1) then
cc        jp = j
cc        dhy = grid_params%dy(jg-1)
cc      endif
cc
cc      if (hz == 1) then
cc        km = k
cc        dhz = grid_params%dz(kg)
cc      elseif (hz == -1) then
cc        kp = k
cc        dhz = grid_params%dz(kg-1)
cc      endif
cc
cc      gx = (arr(ip,j,k)-arr(im,j,k))/dhx
cc      gy = (arr(i,jp,k)-arr(i,jm,k))/dhy
cc      gz = (arr(i,j,kp)-arr(i,j,km))/dhz
cc
cc      cov = (/ gx,gy,gz /)
cc
cc      return

c     2nd order formulas

      !X comp
      select case(hx)
      case(0)
        if (i==0) then
          dx1=grid_params%dx(ig)
          dx2=grid_params%dx(ig+1)
          gx = (-arr(i+2,j,k)*dx1/dx2/(dx1+dx2)
     .          +arr(i+1,j,k)*(dx1+dx2)/dx1/dx2
     .          -arr(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
        elseif (i==nx+1) then
          dx1=grid_params%dx(ig-1)
          dx2=grid_params%dx(ig-2)
          gx = -(-arr(i-2,j,k)*dx1/dx2/(dx1+dx2)
     .           +arr(i-1,j,k)*(dx1+dx2)/dx1/dx2
     .           -arr(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
        else
          dx1=grid_params%dx(ig-1)
          dx2=grid_params%dx(ig)

cc          gx = ( arr(i+1,j,k)*dx1/dx2/(dx1+dx2)
cc     .          +arr(i  ,j,k)*(dx2-dx1)/dx2/dx1
cc     .          -arr(i-1,j,k)*dx2/dx1/(dx1+dx2))
          gx = (arr(i+1,j,k)-arr(i-1,j,k))/(dx1+dx2)
        endif
      case(1)
        gx = (arr(i+1,j,k)-arr(i  ,j,k))/grid_params%dx(ig  )
      case(-1)
        gx = (arr(i  ,j,k)-arr(i-1,j,k))/grid_params%dx(ig-1)
      end select

      !Y comp
      select case (hy)
      case(0)
        if (j==0) then
          dy1=grid_params%dy(jg)
          dy2=grid_params%dy(jg+1)
          gy = (-arr(i,j+2,k)*dy1/dy2/(dy1+dy2)
     .          +arr(i,j+1,k)*(dy1+dy2)/dy1/dy2
     .          -arr(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
        elseif (j==ny+1) then
          dy1=grid_params%dy(jg-1)
          dy2=grid_params%dy(jg-2)
          gy = -(-arr(i,j-2,k)*dy1/dy2/(dy1+dy2)
     .           +arr(i,j-1,k)*(dy1+dy2)/dy1/dy2
     .           -arr(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
        else
          dy1=grid_params%dy(jg-1)
          dy2=grid_params%dy(jg)

cc          gy = ( arr(i,j+1,k)*dy1/dy2/(dy1+dy2)
cc     .          +arr(i,j  ,k)*(dy2-dy1)/dy2/dy1
cc     .          -arr(i,j-1,k)*dy2/dy1/(dy1+dy2))
          gy = (arr(i,j+1,k)-arr(i,j-1,k))/(dy1+dy2)
        endif
      case(1)
        gy = (arr(i,j+1,k)-arr(i,j  ,k))/grid_params%dy(jg  )
      case(-1)
        gy = (arr(i,j  ,k)-arr(i,j-1,k))/grid_params%dy(jg-1)
      end select

      !Z comp
      select case(hz)
      case(0)
        if (k==0) then
          dz1=grid_params%dz(kg)
          dz2=grid_params%dz(kg+1)
          gz = (-arr(i,j,k+2)*dz1/dz2/(dz1+dz2)
     .          +arr(i,j,k+1)*(dz1+dz2)/dz1/dz2
     .          -arr(i,j,k  )*(1./dz1+1./(dz1+dz2)))
        elseif (k==nz+1) then
          dz1=grid_params%dz(kg-1)
          dz2=grid_params%dz(kg-2)
          gz = -(-arr(i,j,k-2)*dz1/dz2/(dz1+dz2)
     .           +arr(i,j,k-1)*(dz1+dz2)/dz1/dz2
     .           -arr(i,j,k  )*(1./dz1+1./(dz1+dz2)))
        else
          dz1=grid_params%dz(kg-1)
          dz2=grid_params%dz(kg)

cc          gz = ( arr(i,j,k+1)*dz1/dz2/(dz1+dz2)
cc     .          +arr(i,j,k  )*(dz2-dz1)/dz2/dz1
cc     .          -arr(i,j,k-1)*dz2/dz1/(dz1+dz2))
          gz = (arr(i,j,k+1)-arr(i,j,k-1))/(dz1+dz2)
        endif
      case(1)
        gz = (arr(i,j,k+1)-arr(i,j,k  ))/grid_params%dz(kg  )
      case(-1)
        gz = (arr(i,j,k  )-arr(i,j,k-1))/grid_params%dz(kg-1)
      end select

      cov = (/ gx,gy,gz /)

c     End program

      end function grad_ijk

c     grad_mesh
c     ###############################################################
      function grad_mesh(igrid,arr) result(cov)

c     ---------------------------------------------------------------
c     Calculates grad(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: arr(0:,0:,0:)
     .          ,cov(0:size(arr,1)-1,0:size(arr,2)-1,0:size(arr,3)-1,3)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            cov(i,j,k,:)=grad_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,arr)
          enddo
        enddo
      enddo

c     End program

      end function grad_mesh

c     curl_ijk
c     ###############################################################
      function curl_ijk(i,j,k,nx,ny,nz,igx,igy,igz,a,hex,hey,hez)
     .         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,comp,nx,ny,nz,igx,igy,igz

      real(8)    :: cnv(3)

      real(8)    :: a(0:nx+1,0:ny+1,0:nz+1,3)

      integer,optional :: hex,hey,hez

c     Local variables

      integer :: icomp
      real(8) :: gv(3,3)

c     Begin program

      do icomp=1,3
        gv(:,icomp) = grad_ijk(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,a(:,:,:,icomp),hex=hex,hey=hey,hez=hez)
      enddo

      cnv(1) = gv(2,3)-gv(3,2)
      cnv(2) = gv(3,1)-gv(1,3)
      cnv(3) = gv(1,2)-gv(2,1)

c     End program

      end function curl_ijk

c     curl_mesh
c     ###############################################################
      function curl_mesh(igrid,cov) result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cnv(0:size(cov,1)-1,0:size(cov,2)-1,0:size(cov,3)-1,3)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            cnv(i,j,k,:)=curl_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,cov)
          enddo
        enddo
      enddo

c     End program

      end function curl_mesh

c     div_ijk
c     ###############################################################
      function div_ijk(i,j,k,nx,ny,nz,igx,igy,igz,cnv,hex,hey,hez,sp)
     .         result(div)

      implicit none

c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cnv(0:nx+1,0:ny+1,0:nz+1,3),div
      integer,optional :: hex,hey,hez
      logical,optional :: sp

c     Local variables

      integer :: ig,jg,kg,igrid,ip,im,jp,jm,kp,km,hx,hy,hz
      real(8) :: dxx,dyy,dzz,x0,y0,z0,jacp,jacm,jac0,jach
      logical :: spoint

c     Begin program

      hx = 0; hy = 0; hz = 0

      if (PRESENT(hex)) hx = hex
      if (PRESENT(hey)) hy = hey
      if (PRESENT(hez)) hz = hez

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = .false.
      endif

      igrid = igx

      if (i < 1 ) hx =  1
      if (i > nx) hx = -1

      if (j < 1 ) hy =  1
      if (j > ny) hy = -1

      if (k < 1 ) hz =  1
      if (k > nz) hz = -1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      dxx = 2*grid_params%dxh(ig)
      dyy = 2*grid_params%dyh(jg)
      dzz = 2*grid_params%dzh(kg)

      if (hx == 1) then
        im = i
        dxx = grid_params%dx(ig)
      elseif (hx == -1) then
        ip = i
        dxx = grid_params%dx(ig-1)
      endif

      if (hy == 1) then
        jm = j
        dyy = grid_params%dy(jg)
      elseif (hy == -1) then
        jp = j
        dyy = grid_params%dy(jg-1)
      endif

      if (hz == 1) then
        km = k
        dzz = grid_params%dz(kg)
      elseif (hz == -1) then
        kp = k
        dzz = grid_params%dz(kg-1)
      endif

      if (isSP2(i,igx).and.hx == 0.and.spoint) then
        jacp = gmetric%grid(igrid)%jac(i+1,j,k)
        jach = 0.5*(jacp+jac0)   !Only good for cylindrical-like geom.

        div = ((cnv(i+1,j  ,k  ,1)/jacp
     .         +cnv(i  ,j  ,k  ,1)/jac0)*jach        )/dxx
     .        +(cnv(i  ,j+1,k  ,2)-cnv(i  ,j-1,k  ,2))/dyy
     .        +(cnv(i  ,j  ,k+1,3)-cnv(i  ,j  ,k-1,3))/dzz
      else
        div =  (cnv(ip,j ,k ,1)-cnv(im,j ,k ,1))/dxx
     .        +(cnv(i ,jp,k ,2)-cnv(i ,jm,k ,2))/dyy
     .        +(cnv(i ,j ,kp,3)-cnv(i ,j ,km,3))/dzz
      endif

c     Affect calculation by inverse jacobian

cc      select case(half_elem)
cc      case(1)
cccc        im = i
cccc        dxx = grid_params%dx(ig)
cc        jac0 = 0.5*(gmetric%grid(igrid)%jac(ip,j,k)
cc     .             +gmetric%grid(igrid)%jac(i ,j,k))
cc      case(2)
cccc        jm = j
cccc        dyy = grid_params%dy(jg)
cc        jac0 = 0.5*(gmetric%grid(igrid)%jac(i,jp,k)
cc     .             +gmetric%grid(igrid)%jac(i,j ,k))
cc      case(3)
cccc        km = k
cccc        dzz = grid_params%dz(kg)
cc        jac0 = 0.5*(gmetric%grid(igrid)%jac(i,j,kp)
cc     .             +gmetric%grid(igrid)%jac(i,j,k ))
cc      case default
        jac0 = gmetric%grid(igrid)%jac(i,j,k)
cc      end select

      div = div/jac0
      
c     End 

      end function div_ijk

c     div_mesh
c     ###############################################################
      function div_mesh(igrid,cnv) result(scl)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cnv(0:,0:,0:,:)
     .          ,scl(0:size(cnv,1)-1,0:size(cnv,2)-1,0:size(cnv,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      scl = 0d0

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            scl(i,j,k) = div_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,cnv)
          enddo
        enddo
      enddo

c     End program

      end function div_mesh

      end module grid_nabla_ops

c grid_nabla2_ops
c ######################################################################
      module grid_nabla2_ops

        use grid_create

      contains

c     gradDiv
c     #####################################################################
      function gradDiv(i,j,k,nx,ny,nz,igx,igy,igz,da,vol,cnv)
     .         result(vec)

c     ---------------------------------------------------------------
c     Calculates grad(div(A)) in general non-orthogonal
c     coordinates in a compact form. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: vec(3)
      real(8)    :: da(0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional,intent(IN) :: vol,cnv

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: divip,divim,divjp,divjm,divkp,divkm
     .             ,jacip,jacim,jacjp,jacjm,jackp,jackm
     .             ,jac

      logical    :: vol_wgt,cnvt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      cnvt = .false.
      if (PRESENT(cnv)) cnvt = cnv

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      jacip  = gmetric%grid(igx)%jac(ip,j,k)
      jacim  = gmetric%grid(igx)%jac(im,j,k)
      jacjp  = gmetric%grid(igx)%jac(i,jp,k)
      jacjm  = gmetric%grid(igx)%jac(i,jm,k)
      jackp  = gmetric%grid(igx)%jac(i,j,kp)
      jackm  = gmetric%grid(igx)%jac(i,j,km)
      jac    = gmetric%grid(igx)%jac(i,j,k )

c     Fluxes at faces for calculation of grad[div(da)]

      !!Divergence at faces i+-1/2, etc.
      divip = (da(ip,j ,k,1)-da(i ,j ,k,1))/grid_params%dx(ig)
     .       +(da(i ,jp,k,2)-da(i ,jm,k,2)
     .        +da(ip,jp,k,2)-da(ip,jm,k,2))/grid_params%dyh(jg)*0.25
     .       +(da(i ,j,kp,3)-da(i ,j,km,3)
     .        +da(ip,j,kp,3)-da(ip,j,km,3))/grid_params%dzh(kg)*0.25
      divim = (da(i ,j ,k,1)-da(im,j ,k,1))/grid_params%dx(ig-1)
     .       +(da(i ,jp,k,2)-da(i ,jm,k,2)
     .        +da(im,jp,k,2)-da(im,jm,k,2))/grid_params%dyh(jg)*0.25
     .       +(da(i ,j,kp,3)-da(i ,j,km,3)
     .        +da(im,j,kp,3)-da(im,j,km,3))/grid_params%dzh(kg)*0.25

      divjp = (da(ip,j ,k,1)-da(im,j ,k,1)
     .        +da(ip,jp,k,1)-da(im,jp,k,1))/grid_params%dxh(ig)*0.25
     .       +(da(i ,jp,k,2)-da(i ,j ,k,2))/grid_params%dy(jg)
     .       +(da(i,j ,kp,3)-da(i,j ,km,3)
     .        +da(i,jp,kp,3)-da(i,jp,km,3))/grid_params%dzh(kg)*0.25
      divjm = (da(ip,j ,k,1)-da(im,j ,k,1)
     .        +da(ip,jm,k,1)-da(im,jm,k,1))/grid_params%dxh(ig)*0.25
     .       +(da(i ,j ,k,2)-da(i ,jm,k,2))/grid_params%dy(jg-1)
     .       +(da(i,j ,kp,3)-da(i,j ,km,3)
     .        +da(i,jm,kp,3)-da(i,jm,km,3))/grid_params%dzh(kg)*0.25

      divkp = (da(ip,j,k ,1)-da(im,j,k ,1)
     .        +da(ip,j,kp,1)-da(im,j,kp,1))/grid_params%dxh(ig)*0.25
     .       +(da(i,jp,k ,2)-da(i,jm,k ,2)
     .        +da(i,jp,kp,2)-da(i,jm,kp,2))/grid_params%dyh(jg)*0.25
     .       +(da(i,j ,kp,3)-da(i,j ,k ,3))/grid_params%dz(kg)
      divkm = (da(ip,j,k ,1)-da(im,j,k ,1)
     .        +da(ip,j,km,1)-da(im,j,km,1))/grid_params%dxh(ig)*0.25
     .       +(da(i,jp,k ,2)-da(i,jm,k ,2)
     .        +da(i,jp,km,2)-da(i,jm,km,2))/grid_params%dyh(jg)*0.25
     .       +(da(i,j ,k ,3)-da(i,j ,km,3))/grid_params%dz(kg-1)

      divip = 2*divip/(jac+jacip)
      if (.not.isSP(i,j,k,igx,igy,igz)) then
        divim = 2*divim/(jac+jacim)
      else
        divim = 0d0
      endif

      divjp = 2*divjp/(jac+jacjp)
      divjm = 2*divjm/(jac+jacjm)

      divkp = 2*divkp/(jac+jackp)
      divkm = 2*divkm/(jac+jackm)

c     Calculate gradient

      vec(1) = (divip - divim)/grid_params%dxh(ig)
      vec(2) = (divjp - divjm)/grid_params%dyh(jg)
      vec(3) = (divkp - divkm)/grid_params%dzh(kg)

      if (vol_wgt) vec=vec*gmetric%grid(igx)%dvol(i,j,k)

      if (cnvt) vec = matmul(gmetric%grid(igx)%gsup(i,j,k,:,:),vec)

      end function gradDiv

c     lap_flx_x_ijk
c     ###############################################################
      function lap_flx_x_ijk(i,j,k,igrid,arr,diff,vertex) result(flx)
c     ---------------------------------------------------------------
c     Calculates X-flux for scalar laplacian at faces (i+1/2,j,k)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx

      real(8),optional,target :: diff(0:,0:,0:)
      logical,optional,intent(IN) :: vertex

c     Local variables

      integer    :: im,jm,km,ip,jp,kp,ig,jg,kg

      real(8)    :: d_xx_ip
     .             ,d_xy_ipjp,d_xy_ipjm
     .             ,d_xz_ipkp,d_xz_ipkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      real(8),dimension(:,:,:),pointer :: dff

      logical  :: vrtx_met

c     Begin program

      if (PRESENT(vertex)) then
        vrtx_met = vertex
      else
        vrtx_met = .true.
      endif

      if (PRESENT(diff)) dff => diff

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (vrtx_met) then
        !Vertex metrics
        gsup => gmetric%grid(igrid)%gsup_v

        d_xx_ip = 0.25*(gsup(i,j ,k ,1,1)
     .                 +gsup(i,jm,k ,1,1)
     .                 +gsup(i,j ,km,1,1)
     .                 +gsup(i,jm,km,1,1))

        d_xy_ipjp = 0.5*(gsup(i,j ,k ,1,2)
     .                  +gsup(i,j ,km,1,2))
        d_xy_ipjm = 0.5*(gsup(i,jm,k ,1,2)
     .                  +gsup(i,jm,km,1,2))

        d_xz_ipkp = 0.5*(gsup(i,j ,k ,1,3)
     .                  +gsup(i,jm,k ,1,3))
        d_xz_ipkm = 0.5*(gsup(i,j ,km,1,3)
     .                  +gsup(i,jm,km,1,3))

        if (PRESENT(diff)) then
          d_xx_ip = 0.5*(dff(i ,j,k)+dff(ip,j,k))*d_xx_ip
          d_xy_ipjp = 0.25*(dff(i ,j ,k)          
     .                     +dff(ip,j ,k)          
     .                     +dff(i ,jp,k)          
     .                     +dff(ip,jp,k))*d_xy_ipjp
          d_xy_ipjm = 0.25*(dff(i ,j ,k)          
     .                     +dff(ip,j ,k)          
     .                     +dff(i ,jm,k)          
     .                     +dff(ip,jm,k))*d_xy_ipjm
          d_xz_ipkp = 0.25*(dff(i ,j,k )           
     .                     +dff(ip,j,k )           
     .                     +dff(i ,j,kp)           
     .                     +dff(ip,j,kp))*d_xz_ipkp
          d_xz_ipkm = 0.25*(dff(i ,j,k )           
     .                     +dff(ip,j,k )           
     .                     +dff(i ,j,km)           
     .                     +dff(ip,j,km))*d_xz_ipkm
        endif

      else
        !Cell-centered metrics
        gsup => gmetric%grid(igrid)%gsup

        if (PRESENT(diff)) then
          d_xx_ip = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
     .                  +dff(ip,j,k)*gsup(ip,j,k,1,1))
          d_xy_ipjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                     +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                     +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                     +dff(ip,jp,k)*gsup(ip,jp,k,1,2))
          d_xy_ipjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                     +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                     +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
     .                     +dff(ip,jm,k)*gsup(ip,jm,k,1,2))
          d_xz_ipkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                     +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                     +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                     +dff(ip,j,kp)*gsup(ip,j,kp,1,3))
          d_xz_ipkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                     +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                     +dff(i ,j,km)*gsup(i ,j,km,1,3)
     .                     +dff(ip,j,km)*gsup(ip,j,km,1,3))
        else
          d_xx_ip = 0.5*(gsup(i ,j,k,1,1)
     .                  +gsup(ip,j,k,1,1))

          d_xy_ipjp = 0.25*(gsup(i ,j ,k,1,2)
     .                     +gsup(ip,j ,k,1,2)
     .                     +gsup(i ,jp,k,1,2)
     .                     +gsup(ip,jp,k,1,2))
          d_xy_ipjm = 0.25*(gsup(i ,j ,k,1,2)
     .                     +gsup(ip,j ,k,1,2)
     .                     +gsup(i ,jm,k,1,2)
     .                     +gsup(ip,jm,k,1,2))

          d_xz_ipkp = 0.25*(gsup(i ,j,k ,1,3)
     .                     +gsup(ip,j,k ,1,3)
     .                     +gsup(i ,j,kp,1,3)
     .                     +gsup(ip,j,kp,1,3))
          d_xz_ipkm = 0.25*(gsup(i ,j,k ,1,3)
     .                     +gsup(ip,j,k ,1,3)
     .                     +gsup(i ,j,km,1,3)
     .                     +gsup(ip,j,km,1,3))
        endif
      endif
                                                                       
      flx   = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))

c     End program

      end function lap_flx_x_ijk

c     lap_flx_y_ijk
c     ###############################################################
      function lap_flx_y_ijk(i,j,k,igrid,arr,diff,vertex) result(flx)
c     ---------------------------------------------------------------
c     Calculates Y-flux for scalar laplacian at faces (i,j+1/2,k)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx

      real(8),optional,target :: diff(0:,0:,0:)
      logical,optional,intent(IN) :: vertex

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_yy_jp
     .             ,d_xy_ipjp,d_xy_imjp
     .             ,d_yz_jpkp,d_yz_jpkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup
      real(8),dimension(:,:,:),pointer :: dff

      logical  :: vrtx_met

c     Begin program

      if (PRESENT(vertex)) then
        vrtx_met = vertex
      else
        vrtx_met = .true.
      endif

      if (PRESENT(diff)) dff => diff

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (vrtx_met) then
        !Vertex metrics
        gsup => gmetric%grid(igrid)%gsup_v

        d_yy_jp = 0.25*(gsup(i ,j,k ,2,2)
     .                 +gsup(im,j,k ,2,2)
     .                 +gsup(i ,j,km,2,2)
     .                 +gsup(im,j,km,2,2))

        d_xy_ipjp = 0.5*(gsup(i ,j,k ,1,2)
     .                  +gsup(i ,j,km,1,2))
        d_xy_imjp = 0.5*(gsup(im,j,k ,1,2)
     .                  +gsup(im,j,km,1,2))

        d_yz_jpkp = 0.5*(gsup(i ,j,k ,2,3)
     .                  +gsup(im,j,k ,2,3))
        d_yz_jpkm = 0.5*(gsup(i ,j,km,2,3)
     .                  +gsup(im,j,km,2,3))

        !To fix singularity with gsup(2,2)
        if (isSP2(i,igrid)) then
          gsup => gmetric%grid(igrid)%gsup
          d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
     .                  +gsup(i,jp,k,2,2))
        endif

        if (PRESENT(diff)) then
          d_yy_jp = 0.5*(dff(i,j ,k)+dff(i,jp,k))*d_yy_jp
          d_xy_ipjp = 0.25*(dff(i ,j ,k)          
     .                     +dff(ip,j ,k)          
     .                     +dff(i ,jp,k)          
     .                     +dff(ip,jp,k))*d_xy_ipjp
          d_xy_imjp = 0.25*(dff(i ,j ,k)          
     .                     +dff(im,j ,k)          
     .                     +dff(i ,jp,k)          
     .                     +dff(im,jp,k))*d_xy_imjp
          d_yz_jpkp = 0.25*(dff(i,j ,k )           
     .                     +dff(i,jp,k )           
     .                     +dff(i,j ,kp)           
     .                     +dff(i,jp,kp))*d_yz_jpkp
          d_yz_jpkm = 0.25*(dff(i,j ,k )           
     .                     +dff(i,jp,k )           
     .                     +dff(i,j ,km)           
     .                     +dff(i,jp,km))*d_yz_jpkm
        endif
      else
        !Cell-centered metrics
        gsup => gmetric%grid(igrid)%gsup

        if (PRESENT(diff)) then
          d_yy_jp = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
     .                  +dff(i,jp,k)*gsup(i,jp,k,2,2))
          d_xy_ipjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                     +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                     +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                     +dff(ip,jp,k)*gsup(ip,jp,k,1,2))
          d_xy_imjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                     +dff(im,j ,k)*gsup(im,j ,k,1,2)
     .                     +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                     +dff(im,jp,k)*gsup(im,jp,k,1,2))
          d_yz_jpkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                     +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                     +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                     +dff(i,jp,kp)*gsup(i,jp,kp,2,3))
          d_yz_jpkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                     +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                     +dff(i,j ,km)*gsup(i,j ,km,2,3)
     .                     +dff(i,jp,km)*gsup(i,jp,km,2,3))
        else
          d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
     .                  +gsup(i,jp,k,2,2))

          d_xy_ipjp = 0.25*(gsup(i ,j ,k,1,2)
     .                     +gsup(ip,j ,k,1,2)
     .                     +gsup(i ,jp,k,1,2)
     .                     +gsup(ip,jp,k,1,2))
          d_xy_imjp = 0.25*(gsup(i ,j ,k,1,2)
     .                     +gsup(im,j ,k,1,2)
     .                     +gsup(i ,jp,k,1,2)
     .                     +gsup(im,jp,k,1,2))

          d_yz_jpkp = 0.25*(gsup(i,j ,k ,2,3)
     .                     +gsup(i,jp,k ,2,3)
     .                     +gsup(i,j ,kp,2,3)
     .                     +gsup(i,jp,kp,2,3))
          d_yz_jpkm = 0.25*(gsup(i,j ,k ,2,3)
     .                     +gsup(i,jp,k ,2,3)
     .                     +gsup(i,j ,km,2,3)
     .                     +gsup(i,jp,km,2,3))
        endif
      endif
                                                                       
      flx =   d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
     .                  -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))

c     End program

      end function lap_flx_y_ijk

c     lap_flx_z_ijk
c     ###############################################################
      function lap_flx_z_ijk(i,j,k,igrid,arr,diff,vertex) result(flx)
c     ---------------------------------------------------------------
c     Calculates Z-flux for scalar laplacian at faces (i,j,k-1/2)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx

      real(8),optional,target :: diff(0:,0:,0:)
      logical,optional,intent(IN) :: vertex

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_zz_kp
     .             ,d_xz_ipkp,d_xz_imkp
     .             ,d_yz_jpkp,d_yz_jmkp

      real(8),dimension(:,:,:,:,:),pointer :: gsup
      real(8),dimension(:,:,:),pointer :: dff

      logical  :: vrtx_met

c     Begin program

      if (PRESENT(vertex)) then
        vrtx_met = vertex
      else
        vrtx_met = .true.
      endif

      if (PRESENT(diff)) dff => diff

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (vrtx_met) then
        !Vertex metrics
        gsup => gmetric%grid(igrid)%gsup_v

        d_zz_kp = 0.25*(gsup(i ,j ,k,3,3)
     .                 +gsup(im,j ,k,3,3)
     .                 +gsup(i ,jm,k,3,3)
     .                 +gsup(im,jm,k,3,3))

        d_xz_ipkp = 0.5*(gsup(i ,j ,k,1,3)
     .                  +gsup(i ,jm,k,1,3))
        d_xz_imkp = 0.5*(gsup(im,j ,k,1,3)
     .                  +gsup(im,jm,k,1,3))

        d_yz_jpkp = 0.5*(gsup(i ,j ,k,2,3)
     .                  +gsup(im,j ,k,2,3))
        d_yz_jmkp = 0.5*(gsup(i ,jm,k,2,3)
     .                  +gsup(im,jm,k,2,3))

        if (PRESENT(diff)) then
          d_zz_kp = 0.5*(dff(i,j,k )+dff(i,j,kp))*d_zz_kp
          d_xz_ipkp = 0.25*(dff(i ,j,k )           
     .                     +dff(ip,j,k )           
     .                     +dff(i ,j,kp)           
     .                     +dff(ip,j,kp))*d_xz_ipkp
          d_xz_imkp = 0.25*(dff(i ,j,k )           
     .                     +dff(im,j,k )           
     .                     +dff(i ,j,kp)           
     .                     +dff(im,j,kp))*d_xz_imkp
          d_yz_jpkp = 0.25*(dff(i,j ,k )           
     .                     +dff(i,jp,k )           
     .                     +dff(i,j ,kp)           
     .                     +dff(i,jp,kp))*d_yz_jpkp
          d_yz_jmkp = 0.25*(dff(i,j ,k )           
     .                     +dff(i,jm,k )           
     .                     +dff(i,j ,kp)           
     .                     +dff(i,jm,kp))*d_yz_jmkp
        endif
      else
        !Cell-centered metrics
        gsup => gmetric%grid(igrid)%gsup

        if (PRESENT(diff)) then
          d_zz_kp = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
     .                  +dff(i,j,kp)*gsup(i,j,kp,3,3))
          d_xz_ipkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                     +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                     +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                     +dff(ip,j,kp)*gsup(ip,j,kp,1,3))
          d_xz_imkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                     +dff(im,j,k )*gsup(im,j,k ,1,3)
     .                     +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                     +dff(im,j,kp)*gsup(im,j,kp,1,3))
          d_yz_jpkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                     +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                     +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                     +dff(i,jp,kp)*gsup(i,jp,kp,2,3))
          d_yz_jmkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                     +dff(i,jm,k )*gsup(i,jm,k ,2,3)
     .                     +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                     +dff(i,jm,kp)*gsup(i,jm,kp,2,3))
        else
          d_zz_kp = 0.5*(gsup(i,j,k ,3,3)
     .                  +gsup(i,j,kp,3,3))

          d_xz_ipkp = 0.25*(gsup(i ,j,k ,1,3)
     .                     +gsup(ip,j,k ,1,3)
     .                     +gsup(i ,j,kp,1,3)
     .                     +gsup(ip,j,kp,1,3))
          d_xz_imkp = 0.25*(gsup(i ,j,k ,1,3)
     .                     +gsup(im,j,k ,1,3)
     .                     +gsup(i ,j,kp,1,3)
     .                     +gsup(im,j,kp,1,3))

          d_yz_jpkp = 0.25*(gsup(i,j ,k ,2,3)
     .                     +gsup(i,jp,k ,2,3)
     .                     +gsup(i,j ,kp,2,3)
     .                     +gsup(i,jp,kp,2,3))
          d_yz_jmkp = 0.25*(gsup(i,j ,k ,2,3)
     .                     +gsup(i,jm,k ,2,3)
     .                     +gsup(i,j ,kp,2,3)
     .                     +gsup(i,jm,kp,2,3))
        endif
      endif
                                                                       
      flx   = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))

c     End program

      end function lap_flx_z_ijk

ccc     lap_ndiff
ccc     ###############################################################
cc      function lap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol
cc     .                  ,vertex_metrics)
ccc     ---------------------------------------------------------------
ccc     Calculates lap(arr) at cell centers in general non-orthog.
ccc     coordinates, preserving the SPD property.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: arr (0:nx+1,0:ny+1,0:nz+1),lap_ndiff
cc
cc      logical,optional,intent(IN) :: vol,vertex_metrics
cc
ccc     Local variables
cc
cc      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
cc
cc      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
cc     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
cc     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
cc     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
cc     .             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc      real(8),dimension(:,:,:,:,:),pointer :: gsup
cc
cc      logical  :: vol_wgt,vrtx_met
cc
ccc     Begin program
cc
cc      if (PRESENT(vol)) then
cc        vol_wgt = vol
cc      else
cc        vol_wgt = .false.
cc      endif
cc
cc      if (PRESENT(vertex_metrics)) then
cc        vrtx_met = vertex_metrics
cc      else
cc        vrtx_met = .true.
cc      endif
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      if (vrtx_met) then
cc        !Vertex metrics
cc        gsup => gmetric%grid(igx)%gsup_v
cc
cc        d_xx_ip = 0.25*(gsup(i ,j ,k ,1,1)
cc     .                 +gsup(i ,jm,k ,1,1)
cc     .                 +gsup(i ,j ,km,1,1)
cc     .                 +gsup(i ,jm,km,1,1))
cc        d_xx_im = 0.25*(gsup(im,j ,k ,1,1)
cc     .                 +gsup(im,jm,k ,1,1)
cc     .                 +gsup(im,j ,km,1,1)
cc     .                 +gsup(im,jm,km,1,1))
cc
cc        d_yy_jp = 0.25*(gsup(i ,j ,k ,2,2)
cc     .                 +gsup(im,j ,k ,2,2)
cc     .                 +gsup(i ,j ,km,2,2)
cc     .                 +gsup(im,j ,km,2,2))
cc        d_yy_jm = 0.25*(gsup(i ,jm,k ,2,2)
cc     .                 +gsup(im,jm,k ,2,2)
cc     .                 +gsup(i ,jm,km,2,2)
cc     .                 +gsup(im,jm,km,2,2))
cc
cc        d_zz_kp = 0.25*(gsup(i ,j ,k ,3,3)
cc     .                 +gsup(im,j ,k ,3,3)
cc     .                 +gsup(i ,jm,k ,3,3)
cc     .                 +gsup(im,jm,k ,3,3))
cc        d_zz_km = 0.25*(gsup(i ,j ,km,3,3)
cc     .                 +gsup(im,j ,km,3,3)
cc     .                 +gsup(i ,jm,km,3,3)
cc     .                 +gsup(im,jm,km,3,3))
cc
cc        d_xy_ipjp = 0.5*(gsup(i ,j ,k ,1,2)
cc     .                  +gsup(i ,j ,km,1,2))
cc        d_xy_ipjm = 0.5*(gsup(i ,jm,k ,1,2)
cc     .                  +gsup(i ,jm,km,1,2))
cc        d_xy_imjp = 0.5*(gsup(im,j ,k ,1,2)
cc     .                  +gsup(im,j ,km,1,2))
cc        d_xy_imjm = 0.5*(gsup(im,jm,k ,1,2)
cc     .                  +gsup(im,jm,km,1,2))
cc
cc        d_xz_ipkp = 0.5*(gsup(i ,j ,k ,1,3)
cc     .                  +gsup(i ,jm,k ,1,3))
cc        d_xz_ipkm = 0.5*(gsup(i ,j ,km,1,3)
cc     .                  +gsup(i ,jm,km,1,3))
cc        d_xz_imkp = 0.5*(gsup(im,j ,k ,1,3)
cc     .                  +gsup(im,jm,k ,1,3))
cc        d_xz_imkm = 0.5*(gsup(im,j ,km,1,3)
cc     .                  +gsup(im,jm,km,1,3))
cc
cc        d_yz_jpkp = 0.5*(gsup(i ,j ,k ,2,3)
cc     .                  +gsup(im,j ,k ,2,3))
cc        d_yz_jpkm = 0.5*(gsup(i ,j ,km,2,3)
cc     .                  +gsup(im,j ,km,2,3))
cc        d_yz_jmkp = 0.5*(gsup(i ,jm,k ,2,3)
cc     .                  +gsup(im,jm,k ,2,3))
cc        d_yz_jmkm = 0.5*(gsup(i ,jm,km,2,3)
cc     .                  +gsup(im,jm,km,2,3))
cc
cc        !To fix singularity with gsup(2,2)
cc        if (isSP2(i,igx)) then
cc          gsup => gmetric%grid(igx)%gsup
cc          d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
cc     .                  +gsup(i,jp,k,2,2))
cc          d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
cc     .                  +gsup(i,jm,k,2,2))
cc        endif
cc      else
cc        !Cell-centered metrics
cc        gsup => gmetric%grid(igx)%gsup
cc
cc        d_xx_ip = 0.5*(gsup(i ,j,k,1,1)
cc     .                +gsup(ip,j,k,1,1))
cc        d_xx_im = 0.5*(gsup(i ,j,k,1,1)
cc     .                +gsup(im,j,k,1,1))
cc
cc        d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
cc     .                +gsup(i,jp,k,2,2))
cc        d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
cc     .                +gsup(i,jm,k,2,2))
cc
cc        d_zz_kp = 0.5*(gsup(i,j,k ,3,3)
cc     .                +gsup(i,j,kp,3,3))
cc        d_zz_km = 0.5*(gsup(i,j,k ,3,3)
cc     .                +gsup(i,j,km,3,3))
cc
cc        d_xy_ipjp = 0.25*(gsup(i ,j ,k,1,2)
cc     .                   +gsup(ip,j ,k,1,2)
cc     .                   +gsup(i ,jp,k,1,2)
cc     .                   +gsup(ip,jp,k,1,2))
cc        d_xy_ipjm = 0.25*(gsup(i ,j ,k,1,2)
cc     .                   +gsup(ip,j ,k,1,2)
cc     .                   +gsup(i ,jm,k,1,2)
cc     .                   +gsup(ip,jm,k,1,2))
cc        d_xy_imjp = 0.25*(gsup(i ,j ,k,1,2)
cc     .                   +gsup(im,j ,k,1,2)
cc     .                   +gsup(i ,jp,k,1,2)
cc     .                   +gsup(im,jp,k,1,2))
cc        d_xy_imjm = 0.25*(gsup(i ,j ,k,1,2)
cc     .                   +gsup(im,j ,k,1,2)
cc     .                   +gsup(i ,jm,k,1,2)
cc     .                   +gsup(im,jm,k,1,2))
cc
cc        d_xz_ipkp = 0.25*(gsup(i ,j,k ,1,3)
cc     .                   +gsup(ip,j,k ,1,3)
cc     .                   +gsup(i ,j,kp,1,3)
cc     .                   +gsup(ip,j,kp,1,3))
cc        d_xz_ipkm = 0.25*(gsup(i ,j,k ,1,3)
cc     .                   +gsup(ip,j,k ,1,3)
cc     .                   +gsup(i ,j,km,1,3)
cc     .                   +gsup(ip,j,km,1,3))
cc        d_xz_imkp = 0.25*(gsup(i ,j,k ,1,3)
cc     .                   +gsup(im,j,k ,1,3)
cc     .                   +gsup(i ,j,kp,1,3)
cc     .                   +gsup(im,j,kp,1,3))
cc        d_xz_imkm = 0.25*(gsup(i ,j,k ,1,3)
cc     .                   +gsup(im,j,k ,1,3)
cc     .                   +gsup(i ,j,km,1,3)
cc     .                   +gsup(im,j,km,1,3))
cc
cc        d_yz_jpkp = 0.25*(gsup(i,j ,k ,2,3)
cc     .                   +gsup(i,jp,k ,2,3)
cc     .                   +gsup(i,j ,kp,2,3)
cc     .                   +gsup(i,jp,kp,2,3))
cc        d_yz_jpkm = 0.25*(gsup(i,j ,k ,2,3)
cc     .                   +gsup(i,jp,k ,2,3)
cc     .                   +gsup(i,j ,km,2,3)
cc     .                   +gsup(i,jp,km,2,3))
cc        d_yz_jmkp = 0.25*(gsup(i,j ,k ,2,3)
cc     .                   +gsup(i,jm,k ,2,3)
cc     .                   +gsup(i,j ,kp,2,3)
cc     .                   +gsup(i,jm,kp,2,3))
cc        d_yz_jmkm = 0.25*(gsup(i,j ,k ,2,3)
cc     .                   +gsup(i,jm,k ,2,3)
cc     .                   +gsup(i,j ,km,2,3)
cc     .                   +gsup(i,jm,km,2,3))
cc      endif
cc
cc      flxip = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
cc     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
cc     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
cc     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
cc     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
cc     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
cc     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                       
cc      flxim =-d_xx_im  *(arr(im,j ,k)-arr(i,j ,k))/grid_params%dx(ig-1)
cc     .  +.25*(d_xy_imjp*(arr(im,jp,k)+arr(i,jp,k)               
cc     .                  -arr(im,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_xy_imjm*(arr(im,j ,k)+arr(i,j ,k)               
cc     .                  -arr(im,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
cc     .  +     d_xz_imkp*(arr(im,j,kp)+arr(i,j,kp)               
cc     .                  -arr(im,j,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_xz_imkm*(arr(im,j,k )+arr(i,j,k )               
cc     .                  -arr(im,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                       
cc      flxjp = d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
cc     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
cc     .                  -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
cc     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
cc     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
cc     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                       
cc      flxjm =-d_yy_jm  *(arr(i ,jm,k)-arr(i ,j,k))/grid_params%dy(jg-1)
cc     .  +.25*(d_xy_ipjm*(arr(ip,jm,k)+arr(ip,j,k)               
cc     .                  -arr(i ,jm,k)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xy_imjm*(arr(i ,jm,k)+arr(i ,j,k)               
cc     .                  -arr(im,jm,k)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jmkp*(arr(i,jm,kp)+arr(i,j,kp)               
cc     .                  -arr(i,jm,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_yz_jmkm*(arr(i,jm,k )+arr(i,j,k )               
cc     .                  -arr(i,jm,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                       
cc      flxkp = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
cc     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
cc     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
cc     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
cc     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
cc     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))
cc                                                                       
cc      flxkm =-d_zz_km  *(arr(i ,j,km)-arr(i ,j,k))/grid_params%dz(kg-1)
cc     .  +.25*(d_xz_ipkm*(arr(ip,j,km)+arr(ip,j,k)               
cc     .                  -arr(i ,j,km)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xz_imkm*(arr(i ,j,km)+arr(i ,j,k)               
cc     .                  -arr(im,j,km)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkm*(arr(i,jp,km)+arr(i,jp,k)               
cc     .                  -arr(i,j ,km)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_yz_jmkm*(arr(i,j ,km)+arr(i,j ,k)               
cc     .                  -arr(i,jm,km)-arr(i,jm,k))/grid_params%dy(jg-1))
cc
cc      lap_ndiff= ((flxip -flxim)/grid_params%dxh(ig)
cc     .           +(flxjp -flxjm)/grid_params%dyh(jg)
cc     .           +(flxkp -flxkm)/grid_params%dzh(kg))
cc     .           /gmetric%grid(igx)%jac(i,j,k)
cc
cc      if (vol_wgt) lap_ndiff=lap_ndiff*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function lap_ndiff
c$$$
c$$$c     lap_ndiff
c$$$c     ###############################################################
c$$$      function lap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol
c$$$     .                  ,vertex_metrics)
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates lap(arr) at cell centers in general non-orthog.
c$$$c     coordinates, preserving the SPD property.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8)    :: arr (0:nx+1,0:ny+1,0:nz+1),lap_ndiff
c$$$
c$$$      logical,optional,intent(IN) :: vol,vertex_metrics
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: ip,jp,kp,im,jm,km,ig,jg,kg
c$$$      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
c$$$
c$$$      real(8),dimension(:,:,:,:,:),pointer :: gsup
c$$$
c$$$      logical  :: vol_wgt,vrtx_met
c$$$
c$$$c     Begin program
c$$$
c$$$      if (PRESENT(vol)) then
c$$$        vol_wgt = vol
c$$$      else
c$$$        vol_wgt = .false.
c$$$      endif
c$$$
c$$$      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$      im = i-1
c$$$      jm = j-1
c$$$      km = k-1
c$$$
c$$$      flxip = lap_flx_x_ijk(i ,j,k,igx,arr,vertex=vertex_metrics)
c$$$      flxim = lap_flx_x_ijk(im,j,k,igx,arr,vertex=vertex_metrics)
c$$$
c$$$      flxjp = lap_flx_y_ijk(i,j ,k,igx,arr,vertex=vertex_metrics)
c$$$      flxjm = lap_flx_y_ijk(i,jm,k,igx,arr,vertex=vertex_metrics)
c$$$
c$$$      flxkp = lap_flx_z_ijk(i,j,k ,igx,arr,vertex=vertex_metrics)
c$$$      flxkm = lap_flx_z_ijk(i,j,km,igx,arr,vertex=vertex_metrics)
c$$$
c$$$      lap_ndiff= ((flxip -flxim)/grid_params%dxh(ig)
c$$$     .           +(flxjp -flxjm)/grid_params%dyh(jg)
c$$$     .           +(flxkp -flxkm)/grid_params%dzh(kg))
c$$$     .           /gmetric%grid(igx)%jac(i,j,k)
c$$$
c$$$      if (vol_wgt) lap_ndiff=lap_ndiff*gmetric%grid(igx)%dvol(i,j,k)
c$$$
c$$$c     End program
c$$$
c$$$      end function lap_ndiff

ccc     lap_diff
ccc     ###############################################################
cc      function lap_diff(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol
cc     .                 ,vertex_metrics) result (laplacian)
ccc     ---------------------------------------------------------------
ccc     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
ccc     coordinates, preserving the SPD property.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian
cc
cc      real(8)    :: dff(0:nx+1,0:ny+1,0:nz+1)
cc
cc      logical,optional,intent(IN) :: vol,vertex_metrics
cc
ccc     Local variables
cc
cc      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
cc
cc      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
cc     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
cc     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
cc     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
cc     .             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc      real(8),dimension(:,:,:,:,:),pointer :: gsup
cc
cc      logical    :: vol_wgt,vrtx_met
cc
ccc     Begin program
cc
cc      vol_wgt = .false.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
cc      if (PRESENT(vertex_metrics)) then
cc        vrtx_met = vertex_metrics
cc      else
cc        vrtx_met = .true.
cc      endif
cc
cc      igrid = igx
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc      
cc      if (vrtx_met) then
cc        !Average metrics (defined at vertices)
cc        gsup => gmetric%grid(igrid)%gsup_v
cc
cc        d_xx_ip = 0.25*(gsup(i ,j ,k ,1,1)
cc     .                 +gsup(i ,jm,k ,1,1)
cc     .                 +gsup(i ,j ,km,1,1)
cc     .                 +gsup(i ,jm,km,1,1))
cc        d_xx_im = 0.25*(gsup(im,j ,k ,1,1)
cc     .                 +gsup(im,jm,k ,1,1)
cc     .                 +gsup(im,j ,km,1,1)
cc     .                 +gsup(im,jm,km,1,1))
cc
cc        d_yy_jp = 0.25*(gsup(i ,j ,k ,2,2)
cc     .                 +gsup(im,j ,k ,2,2)
cc     .                 +gsup(i ,j ,km,2,2)
cc     .                 +gsup(im,j ,km,2,2))
cc        d_yy_jm = 0.25*(gsup(i ,jm,k ,2,2)
cc     .                 +gsup(im,jm,k ,2,2)
cc     .                 +gsup(i ,jm,km,2,2)
cc     .                 +gsup(im,jm,km,2,2))
cc
cc        d_zz_kp = 0.25*(gsup(i ,j ,k ,3,3)
cc     .                 +gsup(im,j ,k ,3,3)
cc     .                 +gsup(i ,jm,k ,3,3)
cc     .                 +gsup(im,jm,k ,3,3))
cc        d_zz_km = 0.25*(gsup(i ,j ,km,3,3)
cc     .                 +gsup(im,j ,km,3,3)
cc     .                 +gsup(i ,jm,km,3,3)
cc     .                 +gsup(im,jm,km,3,3))
cc
cc        d_xy_ipjp = 0.5*(gsup(i ,j ,k ,1,2)
cc     .                  +gsup(i ,j ,km,1,2))
cc        d_xy_ipjm = 0.5*(gsup(i ,jm,k ,1,2)
cc     .                  +gsup(i ,jm,km,1,2))
cc        d_xy_imjp = 0.5*(gsup(im,j ,k ,1,2)
cc     .                  +gsup(im,j ,km,1,2))
cc        d_xy_imjm = 0.5*(gsup(im,jm,k ,1,2)
cc     .                  +gsup(im,jm,km,1,2))
cc                         
cc        d_xz_ipkp = 0.5*(gsup(i ,j ,k ,1,3)
cc     .                  +gsup(i ,jm,k ,1,3))
cc        d_xz_ipkm = 0.5*(gsup(i ,j ,km,1,3)
cc     .                  +gsup(i ,jm,km,1,3))
cc        d_xz_imkp = 0.5*(gsup(im,j ,k ,1,3)
cc     .                  +gsup(im,jm,k ,1,3))
cc        d_xz_imkm = 0.5*(gsup(im,j ,km,1,3)
cc     .                  +gsup(im,jm,km,1,3))
cc
cc        d_yz_jpkp = 0.5*(gsup(i ,j ,k ,2,3)
cc     .                  +gsup(im,j ,k ,2,3))
cc        d_yz_jpkm = 0.5*(gsup(i ,j ,km,2,3)
cc     .                  +gsup(im,j ,km,2,3))
cc        d_yz_jmkp = 0.5*(gsup(i ,jm,k ,2,3)
cc     .                  +gsup(im,jm,k ,2,3))
cc        d_yz_jmkm = 0.5*(gsup(i ,jm,km,2,3)
cc     .                  +gsup(im,jm,km,2,3))
cc
cc        !Fix singularity with gsup(2,2)
cc        if (isSP2(i,igx)) then
cc          gsup => gmetric%grid(igrid)%gsup
cc          d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
cc     .                  +gsup(i,jp,k,2,2))
cc          d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
cc     .                  +gsup(i,jm,k,2,2))
cc        endif
cc
cc        !Add diffusion coefficient (defined at cell centers)
cc        d_xx_ip = 0.5*(dff(i ,j,k)+dff(ip,j,k))*d_xx_ip
cc        d_xx_im = 0.5*(dff(i ,j,k)+dff(im,j,k))*d_xx_im
cc        d_yy_jp = 0.5*(dff(i,j ,k)+dff(i,jp,k))*d_yy_jp
cc        d_yy_jm = 0.5*(dff(i,j ,k)+dff(i,jm,k))*d_yy_jm
cc        d_zz_kp = 0.5*(dff(i,j,k )+dff(i,j,kp))*d_zz_kp
cc        d_zz_km = 0.5*(dff(i,j,k )+dff(i,j,km))*d_zz_km
cc
cc        d_xy_ipjp = 0.25*(dff(i ,j ,k)
cc     .                   +dff(ip,j ,k)
cc     .                   +dff(i ,jp,k)
cc     .                   +dff(ip,jp,k))*d_xy_ipjp
cc        d_xy_ipjm = 0.25*(dff(i ,j ,k)          
cc     .                   +dff(ip,j ,k)          
cc     .                   +dff(i ,jm,k)          
cc     .                   +dff(ip,jm,k))*d_xy_ipjm
cc        d_xy_imjp = 0.25*(dff(i ,j ,k)          
cc     .                   +dff(im,j ,k)          
cc     .                   +dff(i ,jp,k)          
cc     .                   +dff(im,jp,k))*d_xy_imjp
cc        d_xy_imjm = 0.25*(dff(i ,j ,k)          
cc     .                   +dff(im,j ,k)          
cc     .                   +dff(i ,jm,k)          
cc     .                   +dff(im,jm,k))*d_xy_imjm
cc                                                
cc        d_xz_ipkp = 0.25*(dff(i ,j,k )          
cc     .                   +dff(ip,j,k )          
cc     .                   +dff(i ,j,kp)          
cc     .                   +dff(ip,j,kp))*d_xz_ipkp
cc        d_xz_ipkm = 0.25*(dff(i ,j,k )           
cc     .                   +dff(ip,j,k )           
cc     .                   +dff(i ,j,km)           
cc     .                   +dff(ip,j,km))*d_xz_ipkm
cc        d_xz_imkp = 0.25*(dff(i ,j,k )           
cc     .                   +dff(im,j,k )           
cc     .                   +dff(i ,j,kp)           
cc     .                   +dff(im,j,kp))*d_xz_imkp
cc        d_xz_imkm = 0.25*(dff(i ,j,k )           
cc     .                   +dff(im,j,k )           
cc     .                   +dff(i ,j,km)           
cc     .                   +dff(im,j,km))*d_xz_imkm
cc                                                
cc        d_yz_jpkp = 0.25*(dff(i,j ,k )          
cc     .                   +dff(i,jp,k )          
cc     .                   +dff(i,j ,kp)          
cc     .                   +dff(i,jp,kp))*d_yz_jpkp
cc        d_yz_jpkm = 0.25*(dff(i,j ,k )           
cc     .                   +dff(i,jp,k )           
cc     .                   +dff(i,j ,km)           
cc     .                   +dff(i,jp,km))*d_yz_jpkm
cc        d_yz_jmkp = 0.25*(dff(i,j ,k )           
cc     .                   +dff(i,jm,k )           
cc     .                   +dff(i,j ,kp)           
cc     .                   +dff(i,jm,kp))*d_yz_jmkp
cc        d_yz_jmkm = 0.25*(dff(i,j ,k )           
cc     .                   +dff(i,jm,k )           
cc     .                   +dff(i,j ,km)           
cc     .                   +dff(i,jm,km))*d_yz_jmkm
cc      else
cc        !Cell-centered metrics
cc        gsup => gmetric%grid(igrid)%gsup
cc
cc        d_xx_ip = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
cc     .                +dff(ip,j,k)*gsup(ip,j,k,1,1))
cc        d_xx_im = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
cc     .                +dff(im,j,k)*gsup(im,j,k,1,1))
cc        d_yy_jp = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
cc     .                +dff(i,jp,k)*gsup(i,jp,k,2,2))
cc        d_yy_jm = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
cc     .                +dff(i,jm,k)*gsup(i,jm,k,2,2))
cc        d_zz_kp = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
cc     .                +dff(i,j,kp)*gsup(i,j,kp,3,3))
cc        d_zz_km = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
cc     .                +dff(i,j,km)*gsup(i,j,km,3,3))
cc
cc        d_xy_ipjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
cc     .                   +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
cc     .                   +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
cc     .                   +dff(ip,jp,k)*gsup(ip,jp,k,1,2))
cc        d_xy_ipjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
cc     .                   +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
cc     .                   +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
cc     .                   +dff(ip,jm,k)*gsup(ip,jm,k,1,2))
cc        d_xy_imjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
cc     .                   +dff(im,j ,k)*gsup(im,j ,k,1,2)
cc     .                   +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
cc     .                   +dff(im,jp,k)*gsup(im,jp,k,1,2))
cc        d_xy_imjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
cc     .                   +dff(im,j ,k)*gsup(im,j ,k,1,2)
cc     .                   +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
cc     .                   +dff(im,jm,k)*gsup(im,jm,k,1,2))
cc                                                        
cc        d_xz_ipkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
cc     .                   +dff(ip,j,k )*gsup(ip,j,k ,1,3)
cc     .                   +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
cc     .                   +dff(ip,j,kp)*gsup(ip,j,kp,1,3))
cc        d_xz_ipkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
cc     .                   +dff(ip,j,k )*gsup(ip,j,k ,1,3)
cc     .                   +dff(i ,j,km)*gsup(i ,j,km,1,3)
cc     .                   +dff(ip,j,km)*gsup(ip,j,km,1,3))
cc        d_xz_imkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
cc     .                   +dff(im,j,k )*gsup(im,j,k ,1,3)
cc     .                   +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
cc     .                   +dff(im,j,kp)*gsup(im,j,kp,1,3))
cc        d_xz_imkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
cc     .                   +dff(im,j,k )*gsup(im,j,k ,1,3)
cc     .                   +dff(i ,j,km)*gsup(i ,j,km,1,3)
cc     .                   +dff(im,j,km)*gsup(im,j,km,1,3))
cc                                                        
cc        d_yz_jpkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
cc     .                   +dff(i,jp,k )*gsup(i,jp,k ,2,3)
cc     .                   +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
cc     .                   +dff(i,jp,kp)*gsup(i,jp,kp,2,3))
cc        d_yz_jpkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
cc     .                   +dff(i,jp,k )*gsup(i,jp,k ,2,3)
cc     .                   +dff(i,j ,km)*gsup(i,j ,km,2,3)
cc     .                   +dff(i,jp,km)*gsup(i,jp,km,2,3))
cc        d_yz_jmkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
cc     .                   +dff(i,jm,k )*gsup(i,jm,k ,2,3)
cc     .                   +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
cc     .                   +dff(i,jm,kp)*gsup(i,jm,kp,2,3))
cc        d_yz_jmkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
cc     .                   +dff(i,jm,k )*gsup(i,jm,k ,2,3)
cc     .                   +dff(i,j ,km)*gsup(i,j ,km,2,3)
cc     .                   +dff(i,jm,km)*gsup(i,jm,km,2,3))
cc      endif
cc
cc      flxip = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
cc     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
cc     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
cc     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
cc     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
cc     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
cc     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxim =-d_xx_im  *(arr(im,j ,k)-arr(i,j ,k))/grid_params%dx(ig-1)
cc     .  +.25*(d_xy_imjp*(arr(im,jp,k)+arr(i,jp,k)               
cc     .                  -arr(im,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_xy_imjm*(arr(im,j ,k)+arr(i,j ,k)               
cc     .                  -arr(im,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
cc     .  +     d_xz_imkp*(arr(im,j,kp)+arr(i,j,kp)               
cc     .                  -arr(im,j,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_xz_imkm*(arr(im,j,k )+arr(i,j,k )               
cc     .                  -arr(im,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxjp = d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
cc     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
cc     .                  -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
cc     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
cc     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
cc     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxjm =-d_yy_jm  *(arr(i ,jm,k)-arr(i ,j,k))/grid_params%dy(jg-1)
cc     .  +.25*(d_xy_ipjm*(arr(ip,jm,k)+arr(ip,j,k)               
cc     .                  -arr(i ,jm,k)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xy_imjm*(arr(i ,jm,k)+arr(i ,j,k)               
cc     .                  -arr(im,jm,k)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jmkp*(arr(i,jm,kp)+arr(i,j,kp)               
cc     .                  -arr(i,jm,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_yz_jmkm*(arr(i,jm,k )+arr(i,j,k )               
cc     .                  -arr(i,jm,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxkp = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
cc     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
cc     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
cc     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
cc     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
cc     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))
cc                                                                     
cc      flxkm =-d_zz_km  *(arr(i ,j,km)-arr(i ,j,k))/grid_params%dz(kg-1)
cc     .  +.25*(d_xz_ipkm*(arr(ip,j,km)+arr(ip,j,k)               
cc     .                  -arr(i ,j,km)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xz_imkm*(arr(i ,j,km)+arr(i ,j,k)               
cc     .                  -arr(im,j,km)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkm*(arr(i,jp,km)+arr(i,jp,k)               
cc     .                  -arr(i,j ,km)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_yz_jmkm*(arr(i,j ,km)+arr(i,j ,k)               
cc     .                  -arr(i,jm,km)-arr(i,jm,k))/grid_params%dy(jg-1))
cc
cc      laplacian= ((flxip -flxim)/grid_params%dxh(ig)
cc     .           +(flxjp -flxjm)/grid_params%dyh(jg)
cc     .           +(flxkp -flxkm)/grid_params%dzh(kg))
cc     .           /gmetric%grid(igx)%jac(i,j,k)
cc
cc      if (vol_wgt) laplacian=laplacian*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function lap_diff

c$$$c     lap_diff
c$$$c     ###############################################################
c$$$      function lap_diff(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol
c$$$     .                 ,vertex_metrics) result (laplacian)
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
c$$$c     coordinates, preserving the SPD property.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian
c$$$
c$$$      real(8)    :: dff(0:nx+1,0:ny+1,0:nz+1)
c$$$
c$$$      logical,optional,intent(IN) :: vol,vertex_metrics
c$$$
c$$$c     Local variables
c$$$
c$$$      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
c$$$
c$$$      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
c$$$
c$$$      logical    :: vol_wgt,vrtx_met
c$$$
c$$$c     Begin program
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$      im = i-1
c$$$      jm = j-1
c$$$      km = k-1
c$$$
c$$$      flxip=lap_flx_x_ijk(i ,j,k,igx,arr,diff=dff,vertex=vertex_metrics)
c$$$      flxim=lap_flx_x_ijk(im,j,k,igx,arr,diff=dff,vertex=vertex_metrics)
c$$$
c$$$      flxjp=lap_flx_y_ijk(i,j ,k,igx,arr,diff=dff,vertex=vertex_metrics)
c$$$      flxjm=lap_flx_y_ijk(i,jm,k,igx,arr,diff=dff,vertex=vertex_metrics)
c$$$
c$$$      flxkp=lap_flx_z_ijk(i,j,k ,igx,arr,diff=dff,vertex=vertex_metrics)
c$$$      flxkm=lap_flx_z_ijk(i,j,km,igx,arr,diff=dff,vertex=vertex_metrics)
c$$$
c$$$      laplacian = ((flxip -flxim)/grid_params%dxh(ig)
c$$$     .            +(flxjp -flxjm)/grid_params%dyh(jg)
c$$$     .            +(flxkp -flxkm)/grid_params%dzh(kg))
c$$$     .            /gmetric%grid(igx)%jac(i,j,k)
c$$$
c$$$      if (vol_wgt) laplacian=laplacian*gmetric%grid(igx)%dvol(i,j,k)
c$$$
c$$$c     End program
c$$$
c$$$      end function lap_diff

c     lap
c     ###############################################################
      function lap(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol
     .            ,vertex_metrics) result (laplacian)
c     ---------------------------------------------------------------
c     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      real(8),optional :: dff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol,vertex_metrics

c     Local variables

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      logical :: vol_wgt,vrtx_met

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      im = i-1
      jm = j-1
      km = k-1

      flxip=lap_flx_x_ijk(i ,j,k,igx,arr,diff=dff,vertex=vertex_metrics)
      flxim=lap_flx_x_ijk(im,j,k,igx,arr,diff=dff,vertex=vertex_metrics)

      flxjp=lap_flx_y_ijk(i,j ,k,igx,arr,diff=dff,vertex=vertex_metrics)
      flxjm=lap_flx_y_ijk(i,jm,k,igx,arr,diff=dff,vertex=vertex_metrics)

      flxkp=lap_flx_z_ijk(i,j,k ,igx,arr,diff=dff,vertex=vertex_metrics)
      flxkm=lap_flx_z_ijk(i,j,km,igx,arr,diff=dff,vertex=vertex_metrics)

      laplacian = ((flxip -flxim)/grid_params%dxh(ig)
     .            +(flxjp -flxjm)/grid_params%dyh(jg)
     .            +(flxkp -flxkm)/grid_params%dzh(kg))
     .            /gmetric%grid(igx)%jac(i,j,k)

      if (vol_wgt) laplacian=laplacian*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap

c     lap_tnsr_flx_x_ijk
c     ###############################################################
      function lap_tnsr_flx_x_ijk(i,j,k,igrid,arr,tnsr) result(flx)
c     ---------------------------------------------------------------
c     Calculates X-flux for scalar laplacian at faces (i+1/2,j,k)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx
      real(8),target :: tnsr(0:,0:,0:,:)

c     Local variables

      integer :: im,jm,km,ip,jp,kp,ig,jg,kg

      real(8) :: d_xx_ip
     .          ,d_xy_ipjp,d_xy_ipjm
     .          ,d_xz_ipkp,d_xz_ipkm

c     Begin program

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      d_xx_ip = 0.5*(tnsr(i ,j,k,1)
     .              +tnsr(ip,j,k,1))

      d_xy_ipjp = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(ip,j ,k,2)
     .                 +tnsr(i ,jp,k,2)
     .                 +tnsr(ip,jp,k,2))
      d_xy_ipjm = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(ip,j ,k,2)
     .                 +tnsr(i ,jm,k,2)
     .                 +tnsr(ip,jm,k,2))

      d_xz_ipkp = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(ip,j,k ,3)
     .                 +tnsr(i ,j,kp,3)
     .                 +tnsr(ip,j,kp,3))
      d_xz_ipkm = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(ip,j,k ,3)
     .                 +tnsr(i ,j,km,3)
     .                 +tnsr(ip,j,km,3))

      flx   = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))

c     End program

      end function lap_tnsr_flx_x_ijk

c     lap_tnsr_flx_y_ijk
c     ###############################################################
      function lap_tnsr_flx_y_ijk(i,j,k,igrid,arr,tnsr) result(flx)
c     ---------------------------------------------------------------
c     Calculates Y-flux for scalar laplacian at faces (i,j+1/2,k)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx
      real(8),target :: tnsr(0:,0:,0:,:)

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_yy_jp
     .             ,d_xy_ipjp,d_xy_imjp
     .             ,d_yz_jpkp,d_yz_jpkm

c     Begin program

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      d_yy_jp = 0.5*(tnsr(i,j ,k,4)
     .              +tnsr(i,jp,k,4))

      d_xy_ipjp = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(ip,j ,k,2)
     .                 +tnsr(i ,jp,k,2)
     .                 +tnsr(ip,jp,k,2))
      d_xy_imjp = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(im,j ,k,2)
     .                 +tnsr(i ,jp,k,2)
     .                 +tnsr(im,jp,k,2))

      d_yz_jpkp = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jp,k ,5)
     .                 +tnsr(i,j ,kp,5)
     .                 +tnsr(i,jp,kp,5))
      d_yz_jpkm = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jp,k ,5)
     .                 +tnsr(i,j ,km,5)
     .                 +tnsr(i,jp,km,5))
                                                                       
      flx =   d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
     .                  -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))

c     End program

      end function lap_tnsr_flx_y_ijk

c     lap_tnsr_flx_z_ijk
c     ###############################################################
      function lap_tnsr_flx_z_ijk(i,j,k,igrid,arr,tnsr) result(flx)
c     ---------------------------------------------------------------
c     Calculates Z-flux for scalar laplacian at faces (i,j,k+1/2)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx
      real(8),target :: tnsr(0:,0:,0:,:)

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_zz_kp
     .             ,d_xz_ipkp,d_xz_imkp
     .             ,d_yz_jpkp,d_yz_jmkp

c     Begin program

      call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      d_zz_kp = 0.5*(tnsr(i,j,k ,6)
     .              +tnsr(i,j,kp,6))

      d_xz_ipkp = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(ip,j,k ,3)
     .                 +tnsr(i ,j,kp,3)
     .                 +tnsr(ip,j,kp,3))
      d_xz_imkp = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(im,j,k ,3)
     .                 +tnsr(i ,j,kp,3)
     .                 +tnsr(im,j,kp,3))

      d_yz_jpkp = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jp,k ,5)
     .                 +tnsr(i,j ,kp,5)
     .                 +tnsr(i,jp,kp,5))
      d_yz_jmkp = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jm,k ,5)
     .                 +tnsr(i,j ,kp,5)
     .                 +tnsr(i,jm,kp,5))
                                                                       
      flx   = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))

c     End program

      end function lap_tnsr_flx_z_ijk

ccc     lap_tnsr
ccc     ###############################################################
cc      function lap_tnsr(i,j,k,nx,ny,nz,igx,igy,igz,tnsr,arr,vol)
cc     .         result (laplacian)
ccc     ---------------------------------------------------------------
ccc     Calculates div(tensor*grad(arr)) at cell centers in general
ccc     non-orthog. coordinates, preserving the SPD property. Tensor
ccc     "tnsr" must be SPD. Only six components of "tnsr" are stored,
ccc     in a lexicographic way. Convention:
ccc
ccc                 ( 1  2  3 )
ccc                 ( -  4  5 )
ccc                 ( -  -  6 )
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)        :: arr (0:nx+1,0:ny+1,0:nz+1),laplacian
cc      real(8),target :: tnsr(0:nx+1,0:ny+1,0:nz+1,6)
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
cc
cc      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
cc     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
cc     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
cc     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
cc     .             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      if (PRESENT(vol)) then
cc        vol_wgt = vol
cc      else
cc        vol_wgt = .false.
cc      endif
cc
cc      igrid = igx
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      d_xx_ip = 0.5*(tnsr(i ,j,k,1)
cc     .              +tnsr(ip,j,k,1))
cc      d_xx_im = 0.5*(tnsr(i ,j,k,1)
cc     .              +tnsr(im,j,k,1))
cc
cc      d_yy_jp = 0.5*(tnsr(i,j ,k,4)
cc     .              +tnsr(i,jp,k,4))
cc      d_yy_jm = 0.5*(tnsr(i,j ,k,4)
cc     .              +tnsr(i,jm,k,4))
cc
cc      d_zz_kp = 0.5*(tnsr(i,j,k ,6)
cc     .              +tnsr(i,j,kp,6))
cc      d_zz_km = 0.5*(tnsr(i,j,k ,6)
cc     .              +tnsr(i,j,km,6))
cc
cc      d_xy_ipjp = 0.25*(tnsr(i ,j ,k,2)
cc     .                 +tnsr(ip,j ,k,2)
cc     .                 +tnsr(i ,jp,k,2)
cc     .                 +tnsr(ip,jp,k,2))
cc      d_xy_ipjm = 0.25*(tnsr(i ,j ,k,2)
cc     .                 +tnsr(ip,j ,k,2)
cc     .                 +tnsr(i ,jm,k,2)
cc     .                 +tnsr(ip,jm,k,2))
cc      d_xy_imjp = 0.25*(tnsr(i ,j ,k,2)
cc     .                 +tnsr(im,j ,k,2)
cc     .                 +tnsr(i ,jp,k,2)
cc     .                 +tnsr(im,jp,k,2))
cc      d_xy_imjm = 0.25*(tnsr(i ,j ,k,2)
cc     .                 +tnsr(im,j ,k,2)
cc     .                 +tnsr(i ,jm,k,2)
cc     .                 +tnsr(im,jm,k,2))
cc
cc      d_xz_ipkp = 0.25*(tnsr(i ,j,k ,3)
cc     .                 +tnsr(ip,j,k ,3)
cc     .                 +tnsr(i ,j,kp,3)
cc     .                 +tnsr(ip,j,kp,3))
cc      d_xz_ipkm = 0.25*(tnsr(i ,j,k ,3)
cc     .                 +tnsr(ip,j,k ,3)
cc     .                 +tnsr(i ,j,km,3)
cc     .                 +tnsr(ip,j,km,3))
cc      d_xz_imkp = 0.25*(tnsr(i ,j,k ,3)
cc     .                 +tnsr(im,j,k ,3)
cc     .                 +tnsr(i ,j,kp,3)
cc     .                 +tnsr(im,j,kp,3))
cc      d_xz_imkm = 0.25*(tnsr(i ,j,k ,3)
cc     .                 +tnsr(im,j,k ,3)
cc     .                 +tnsr(i ,j,km,3)
cc     .                 +tnsr(im,j,km,3))
cc
cc      d_yz_jpkp = 0.25*(tnsr(i,j ,k ,5)
cc     .                 +tnsr(i,jp,k ,5)
cc     .                 +tnsr(i,j ,kp,5)
cc     .                 +tnsr(i,jp,kp,5))
cc      d_yz_jpkm = 0.25*(tnsr(i,j ,k ,5)
cc     .                 +tnsr(i,jp,k ,5)
cc     .                 +tnsr(i,j ,km,5)
cc     .                 +tnsr(i,jp,km,5))
cc      d_yz_jmkp = 0.25*(tnsr(i,j ,k ,5)
cc     .                 +tnsr(i,jm,k ,5)
cc     .                 +tnsr(i,j ,kp,5)
cc     .                 +tnsr(i,jm,kp,5))
cc      d_yz_jmkm = 0.25*(tnsr(i,j ,k ,5)
cc     .                 +tnsr(i,jm,k ,5)
cc     .                 +tnsr(i,j ,km,5)
cc     .                 +tnsr(i,jm,km,5))
cc
cc      flxip = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
cc     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
cc     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
cc     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
cc     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
cc     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
cc     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxim =-d_xx_im  *(arr(im,j ,k)-arr(i,j ,k))/grid_params%dx(ig-1)
cc     .  +.25*(d_xy_imjp*(arr(im,jp,k)+arr(i,jp,k)               
cc     .                  -arr(im,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_xy_imjm*(arr(im,j ,k)+arr(i,j ,k)               
cc     .                  -arr(im,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
cc     .  +     d_xz_imkp*(arr(im,j,kp)+arr(i,j,kp)               
cc     .                  -arr(im,j,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_xz_imkm*(arr(im,j,k )+arr(i,j,k )               
cc     .                  -arr(im,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxjp = d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
cc     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
cc     .                  -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
cc     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
cc     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
cc     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxjm =-d_yy_jm  *(arr(i ,jm,k)-arr(i ,j,k))/grid_params%dy(jg-1)
cc     .  +.25*(d_xy_ipjm*(arr(ip,jm,k)+arr(ip,j,k)               
cc     .                  -arr(i ,jm,k)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xy_imjm*(arr(i ,jm,k)+arr(i ,j,k)               
cc     .                  -arr(im,jm,k)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jmkp*(arr(i,jm,kp)+arr(i,j,kp)               
cc     .                  -arr(i,jm,k )-arr(i,j,k ))/grid_params%dz(kg)
cc     .  +     d_yz_jmkm*(arr(i,jm,k )+arr(i,j,k )               
cc     .                  -arr(i,jm,km)-arr(i,j,km))/grid_params%dz(kg-1))
cc                                                                     
cc      flxkp = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
cc     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
cc     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
cc     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
cc     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
cc     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))
cc                                                                     
cc      flxkm =-d_zz_km  *(arr(i ,j,km)-arr(i ,j,k))/grid_params%dz(kg)
cc     .  +.25*(d_xz_ipkm*(arr(ip,j,km)+arr(ip,j,k)               
cc     .                  -arr(i ,j,km)-arr(i ,j,k))/grid_params%dx(ig)
cc     .  +     d_xz_imkm*(arr(i ,j,km)+arr(i ,j,k)               
cc     .                  -arr(im,j,km)-arr(im,j,k))/grid_params%dx(ig-1)
cc     .  +     d_yz_jpkm*(arr(i,jp,km)+arr(i,jp,k)               
cc     .                  -arr(i,j ,km)-arr(i,j ,k))/grid_params%dy(jg)
cc     .  +     d_yz_jmkm*(arr(i,j ,km)+arr(i,j ,k)               
cc     .                  -arr(i,jm,km)-arr(i,jm,k))/grid_params%dy(jg-1))
cc
cc      laplacian=((flxip -flxim)/grid_params%dxh(ig)
cc     .          +(flxjp -flxjm)/grid_params%dyh(jg)
cc     .          +(flxkp -flxkm)/grid_params%dzh(kg))
cc     .          /gmetric%grid(igx)%jac(i,j,k)
cc
cc      if (vol_wgt) laplacian=laplacian*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function lap_tnsr

c     lap_tnsr
c     ###############################################################
      function lap_tnsr(i,j,k,nx,ny,nz,igx,igy,igz,tnsr,arr,vol)
     .         result (laplacian)
c     ---------------------------------------------------------------
c     Calculates div(tensor*grad(arr)) at cell centers in general
c     non-orthog. coordinates, preserving the SPD property. Tensor
c     "tnsr" must be SPD. Only six components of "tnsr" are stored,
c     in a lexicographic way. Convention:
c
c                 ( 1  2  3 )
c                 ( -  4  5 )
c                 ( -  -  6 )
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)        :: arr (0:nx+1,0:ny+1,0:nz+1),laplacian
      real(8),target :: tnsr(0:nx+1,0:ny+1,0:nz+1,6)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      im = i-1
      jm = j-1
      km = k-1

      flxip=lap_tnsr_flx_x_ijk(i ,j,k,igx,arr,tnsr)
      flxim=lap_tnsr_flx_x_ijk(im,j,k,igx,arr,tnsr)

      flxjp=lap_tnsr_flx_y_ijk(i,j ,k,igx,arr,tnsr)
      flxjm=lap_tnsr_flx_y_ijk(i,jm,k,igx,arr,tnsr)

      flxkp=lap_tnsr_flx_z_ijk(i,j,k ,igx,arr,tnsr)
      flxkm=lap_tnsr_flx_z_ijk(i,j,km,igx,arr,tnsr)

      laplacian = ((flxip -flxim)/grid_params%dxh(ig)
     .            +(flxjp -flxjm)/grid_params%dyh(jg)
     .            +(flxkp -flxkm)/grid_params%dzh(kg))
     .            /gmetric%grid(igx)%jac(i,j,k)

      if (vol_wgt) laplacian=laplacian*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_tnsr

      end module grid_nabla2_ops

c grid_tensor_ops
c ######################################################################
      module grid_tensor_ops

        use grid_create

        use grid_nabla2_ops

c$$$        INTERFACE veclaplacian
c$$$          module procedure veclap!,veclap_diff,veclap_ndiff
c$$$        end INTERFACE

c$$$        INTERFACE veclap_cov
c$$$          module procedure veclap_cov_diff,veclap_cov_ndiff
c$$$        end INTERFACE

        INTERFACE curlcurl
          module procedure curlcurl_ijk,curlcurl_mesh
        end INTERFACE

        INTERFACE curlcurl_cov
          module procedure curlcurl_cov_ijk,curlcurl_cov_mesh
        end INTERFACE

        real(8),pointer,dimension(:,:,:,:) :: vec => null()
     .                                       ,vec1=> null()
     .                                       ,vec2=> null()
        real(8),pointer,dimension(:,:,:)   :: coef=> null()
     .                                       ,sc1 => null()
     .                                       ,sc2 => null()

      contains

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                 ,half_elem,cnv)
     .         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_l^m (if cnv = .false.) or T^lm
c     (if cnv = .true.) of nabla(vec v), with v contravariant, and fills
c     tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

        logical,optional :: cnv

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: tsrc(3,3),gsuper(3,3)
     .               ,jacip,jacjp,jacjm,jackp,jackm
     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac

        logical    :: cntv

c     Begin program

        igrid = igx

        !Defaults

        if (PRESENT(cnv)) then
          cntv = cnv
        else
          cntv = .false.
        endif

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*grid_params%dxh(ig)
        dhy = 2.*grid_params%dyh(jg)
        dhz = 2.*grid_params%dzh(kg)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        !Exceptions

        select case(half_elem)
        case (1)
          dhx = grid_params%dx(ig)

          if (cntv) then
            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,:,:)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,:,:))
          endif

          vxip = ax(ip,j,k)
          vxim = ax(i ,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(i ,j,k)
          vzip = az(ip,j,k)
          vzim = az(i ,j,k)

          vxjp = 0.5*(ax(i,jp,k)+ax(ip,jp,k))
          vxjm = 0.5*(ax(i,jm,k)+ax(ip,jm,k))
          vyjp = 0.5*(ay(i,jp,k)+ay(ip,jp,k))
          vyjm = 0.5*(ay(i,jm,k)+ay(ip,jm,k))
          vzjp = 0.5*(az(i,jp,k)+az(ip,jp,k))
          vzjm = 0.5*(az(i,jm,k)+az(ip,jm,k))
                                      
          vxkp = 0.5*(ax(i,j,kp)+ax(ip,j,kp))
          vxkm = 0.5*(ax(i,j,km)+ax(ip,j,km))
          vykp = 0.5*(ay(i,j,kp)+ay(ip,j,kp))
          vykm = 0.5*(ay(i,j,km)+ay(ip,j,km))
          vzkp = 0.5*(az(i,j,kp)+az(ip,j,kp))
          vzkm = 0.5*(az(i,j,km)+az(ip,j,km))

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

        case (2)
          dhy = grid_params%dy(jg)

          if (cntv)
     .         gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,:,:)
     .                      +gmetric%grid(igrid)%gsup(i,jp,k,:,:))

          vxip = 0.5*(ax(ip,j,k)+ax(ip,jp,k))
          vxim = 0.5*(ax(im,j,k)+ax(im,jp,k))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,jp,k))
          vyim = 0.5*(ay(im,j,k)+ay(im,jp,k))
          vzip = 0.5*(az(ip,j,k)+az(ip,jp,k))
          vzim = 0.5*(az(im,j,k)+az(im,jp,k))

          vxjp = ax(i,jp,k)
          vxjm = ax(i,j ,k)
          vyjp = ay(i,jp,k)
          vyjm = ay(i,j ,k)
          vzjp = az(i,jp,k)
          vzjm = az(i,j ,k)

          vxkp = 0.5*(ax(i,j,kp)+ax(i,jp,kp))
          vxkm = 0.5*(ax(i,j,km)+ax(i,jp,km))
          vykp = 0.5*(ay(i,j,kp)+ay(i,jp,kp))
          vykm = 0.5*(ay(i,j,km)+ay(i,jp,km))
          vzkp = 0.5*(az(i,j,kp)+az(i,jp,kp))
          vzkm = 0.5*(az(i,j,km)+az(i,jp,km))

          tsrc = 0.5*(nabla_v_src(i,jp,k)
     .               +nabla_v_src(i,j ,k))

        case (3)
          dhz = grid_params%dz(kg)

          if (cntv)
     .         gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,:,:)
     .                      +gmetric%grid(igrid)%gsup(i,j,kp,:,:))

          vxip = 0.5*(ax(ip,j,k)+ax(ip,j,kp))
          vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,j,kp))
          vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
          vzip = 0.5*(az(ip,j,k)+az(ip,j,kp))
          vzim = 0.5*(az(im,j,k)+az(im,j,kp))

          vxjp = 0.5*(ax(i,jp,k)+ax(i,jp,kp))
          vxjm = 0.5*(ax(i,jm,k)+ax(i,jm,kp))
          vyjp = 0.5*(ay(i,jp,k)+ay(i,jp,kp))
          vyjm = 0.5*(ay(i,jm,k)+ay(i,jm,kp))
          vzjp = 0.5*(az(i,jp,k)+az(i,jp,kp))
          vzjm = 0.5*(az(i,jm,k)+az(i,jm,kp))

          vxkp = ax(i,j,kp)
          vxkm = ax(i,j,k )
          vykp = ay(i,j,kp)
          vykm = ay(i,j,k )
          vzkp = az(i,j,kp)
          vzkm = az(i,j,k )

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

        case default

          if(cntv) gsuper=gmetric%grid(igrid)%gsup(i,j,k ,:,:)

          vxip = ax(ip,j,k)
          vxim = ax(im,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(im,j,k)
          vzip = az(ip,j,k)
          vzim = az(im,j,k)

          vxjp = ax(i,jp,k)
          vxjm = ax(i,jm,k)
          vyjp = ay(i,jp,k)
          vyjm = ay(i,jm,k)
          vzjp = az(i,jp,k)
          vzjm = az(i,jm,k)

          vxkp = ax(i,j,kp)
          vxkm = ax(i,j,km)
          vykp = ay(i,j,kp)
          vykm = ay(i,j,km)
          vzkp = az(i,j,kp)
          vzkm = az(i,j,km)

          tsrc = nabla_v_src(i,j,k)

        end select

        idhx = 1d0/dhx
        idhy = 1d0/dhy
        idhz = 1d0/dhz

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)*idhx

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)*idhx

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)*idhx

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)*idhy

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)*idhy

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)*idhy

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)*idhz

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)*idhz

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)*idhz

      ! Add geometric source

        tensor = tensor + tsrc

        if (cntv) tensor = matmul(gsuper,tensor)

c     End program

      contains

c     nabla_v_src
c     #############################################################
      function nabla_v_src(i,j,k) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k
        real(8)    :: tensor(3,3)

c     Local variables

        integer    :: l,m
        real(8)    :: cnv(3)

c     Begin program

        cnv = (/ ax(i,j,k),ay(i,j,k),az(i,j,k) /)

        do l=1,3
          do m=1,3
            tensor(l,m) =
     .           cnv(1)* gmetric%grid(igrid)%Gamma(i,j,k,m,1,l)
     .          +cnv(2)* gmetric%grid(igrid)%Gamma(i,j,k,m,2,l)
     .          +cnv(3)* gmetric%grid(igrid)%Gamma(i,j,k,m,3,l)
     .          -cnv(m)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,l)
     .                  +gmetric%grid(igrid)%Gamma(i,j,k,2,2,l)
     .                  +gmetric%grid(igrid)%Gamma(i,j,k,3,3,l))
          enddo
        enddo

c     End program

      end function nabla_v_src

      end function fnabla_v

c     fnabla_v_upwd
c     #############################################################
      function fnabla_v_upwd(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                      ,hex,hey,hez) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + hex,hey,hez = 0 => i,j,k
c       + hex=+-1 --> i+-1/2
c       + hey=+-1 --> j+-1/2
c       + hez=+-1 --> k+-1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,hex,hey,hez,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
        real(8)    :: dhx,dhy,dhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: hessian(3,3,3)
        logical    :: cartsn

c     Begin program

        igrid = igx

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*grid_params%dxh(ig)
        dhy = 2.*grid_params%dyh(jg)
        dhz = 2.*grid_params%dzh(kg)

        hessian = -gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

c     Exceptions

        if (hex == 1) then
          im = i
          dhx = grid_params%dx(ig)
        elseif (hex == -1) then
          ip = i
          dhx = grid_params%dx(ig-1)
        endif

        if (hey == 1) then
          jm = j
          dhy = grid_params%dy(jg)
        elseif (hey == -1) then
          jp = j
          dhy = grid_params%dy(jg-1)
        endif

        if (hez == 1) then
          km = k
          dhz = grid_params%dz(kg)
        elseif (hez == -1) then
          kp = k
          dhz = grid_params%dz(kg-1)
        endif

c     Vectors

        vxx = ax(i,j,k)
        vyy = ay(i,j,k)
        vzz = az(i,j,k)

cc        if (sing_point) then
cc          vxip = ax(ip,j,k)+ax(i,j,k)
cc          vxim = 2.*ax(im,j,k)
cc          vyip = ay(ip,j,k)+ay(i,j,k)
cc          vyim = 2.*ay(im,j,k)
cc          vzip = az(ip,j,k)+az(i,j,k)
cc          vzim = 2.*az(im,j,k)
cc        else
          vxip = ax(ip,j,k)
          vxim = ax(im,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(im,j,k)
          vzip = az(ip,j,k)
          vzim = az(im,j,k)
cc        endif

        vxjp = ax(i,jp,k)
        vxjm = ax(i,jm,k)
        vyjp = ay(i,jp,k)
        vyjm = ay(i,jm,k)
        vzjp = az(i,jp,k)
        vzjm = az(i,jm,k)

        vxkp = ax(i,j,kp)
        vxkm = ax(i,j,km)
        vykp = ay(i,j,kp)
        vykm = ay(i,j,km)
        vzkp = az(i,j,kp)
        vzkm = az(i,j,km)

c     Calculate nabla_v tensor

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)/dhx
     .               + vxx*(hessian(1,1,1)
     .                    + hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(1,1,1)
     .               - vyy* hessian(1,2,1)
     .               - vzz* hessian(1,3,1)

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)/dhx
     .               + vyy*(hessian(1,1,1)
     .                    + hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(2,1,1)
     .               - vyy* hessian(2,2,1)
     .               - vzz* hessian(2,3,1)

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)/dhx
     .               + vzz*(hessian(1,1,1)
     .                    + hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(3,1,1)
     .               - vyy* hessian(3,2,1)
     .               - vzz* hessian(3,3,1)

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)/dhy
     .               + vxx*(hessian(1,1,2)
     .                    + hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(1,1,2)
     .               - vyy* hessian(1,2,2)
     .               - vzz* hessian(1,3,2)

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)/dhy
     .               + vyy*(hessian(1,1,2)
     .                    + hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(2,1,2)
     .               - vyy* hessian(2,2,2)
     .               - vzz* hessian(2,3,2)

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)/dhy
     .               + vzz*(hessian(1,1,2)
     .                    + hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(3,1,2)
     .               - vyy* hessian(3,2,2)
     .               - vzz* hessian(3,3,2)

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)/dhz
     .               + vxx*(hessian(1,1,3)
     .                    + hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(1,1,3)
     .               - vyy* hessian(1,2,3)
     .               - vzz* hessian(1,3,3)

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)/dhz
     .               + vyy*(hessian(1,1,3)
     .                    + hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(2,1,3)
     .               - vyy* hessian(2,2,3)
     .               - vzz* hessian(2,3,3)

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)/dhz
     .               + vzz*(hessian(1,1,3)
     .                    + hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(3,1,3)
     .               - vyy* hessian(3,2,3)
     .               - vzz* hessian(3,3,3)

c     End program

      end function fnabla_v_upwd

c     fnabla_v_cov
c     #############################################################
      function fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                     ,half_elem) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_lm of nabla(vec v), with
c     v covariant, and fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid,jpp,jmm
        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: tsrc(3,3),jacip,jacjp,jacjm,jackp,jackm
     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac
        real(8)    :: cf1,cf2,cf3,cf4

        logical    :: spoint

c     Begin program

        igrid = igx

        !Defaults

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*grid_params%dxh(ig)
        dhy = 2.*grid_params%dyh(jg)
        dhz = 2.*grid_params%dzh(kg)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        spoint = bcSP()
cc        spoint = isSP(i,j,k,igx,igy,igz)
cc        spoint = .false.

cc        !!Coeffs for 4th order first derivative
cc        cf1 = 1./8.
cc        cf2 = 9./8.
cc
cc        !!Coeffs for 4th order second derivative
cc        cf3 =  1./12.
cc        cf4 = 15./12.

        select case(half_elem)
        case (1)
          dhx = grid_params%dx(ig)

          vxip = ax(ip,j,k)
          vxim = ax(i ,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(i ,j,k)
          vzip = az(ip,j,k)
          vzim = az(i ,j,k)

cc          if (spoint) then
ccc FIX PARALLEL for PER BC
cc            !4th order
cc            jpp = jp+1
cc            jmm = jm-1
cc            if (j == ny) jpp = 2
cc            if (j == 1 ) jmm = ny-1
cc
cc            vxjp = 0.5*(-cf1*(ax(ip,jpp,k)+ax(i,jpp,k))
cc     .                  +cf2*(ax(ip,jp ,k)+ax(i,jp ,k))
cc     .                  +cf2*(ax(ip,j  ,k)+ax(i,j  ,k))
cc     .                  -cf1*(ax(ip,jm ,k)+ax(i,jm ,k)))
cc            vxjm = 0.5*(-cf1*(ax(ip,jp ,k)+ax(i,jp ,k))
cc     .                  +cf2*(ax(ip,j  ,k)+ax(i,j  ,k))
cc     .                  +cf2*(ax(ip,jm ,k)+ax(i,jm ,k))
cc     .                  -cf1*(ax(ip,jmm,k)+ax(i,jmm,k)))
cc
cc            vyjp = 0.5*(-cf1*(ay(ip,jpp,k)+ay(i,jpp,k))
cc     .                  +cf2*(ay(ip,jp ,k)+ay(i,jp ,k))
cc     .                  +cf2*(ay(ip,j  ,k)+ay(i,j  ,k))
cc     .                  -cf1*(ay(ip,jm ,k)+ay(i,jm ,k)))
cc            vyjm = 0.5*(-cf1*(ay(ip,jp ,k)+ay(i,jp ,k))
cc     .                  +cf2*(ay(ip,j  ,k)+ay(i,j  ,k))
cc     .                  +cf2*(ay(ip,jm ,k)+ay(i,jm ,k))
cc     .                  -cf1*(ay(ip,jmm,k)+ay(i,jmm,k)))
cc
cc            vzjp = 0.5*(-cf1*(az(ip,jpp,k)+az(i,jpp,k))
cc     .                  +cf2*(az(ip,jp ,k)+az(i,jp ,k))
cc     .                  +cf2*(az(ip,j  ,k)+az(i,j  ,k))
cc     .                  -cf1*(az(ip,jm ,k)+az(i,jm ,k)))
cc            vzjm = 0.5*(-cf1*(az(ip,jp ,k)+az(i,jp ,k))
cc     .                  +cf2*(az(ip,j  ,k)+az(i,j  ,k))
cc     .                  +cf2*(az(ip,jm ,k)+az(i,jm ,k))
cc     .                  -cf1*(az(ip,jmm,k)+az(i,jmm,k)))
cc          else
cc          !2nd order
            vxjp = 0.5*(ax(ip,jp,k)+ax(i,jp,k))
            vxjm = 0.5*(ax(ip,jm,k)+ax(i,jm,k))
            vyjp = 0.5*(ay(ip,jp,k)+ay(i,jp,k))
            vyjm = 0.5*(ay(ip,jm,k)+ay(i,jm,k))
            vzjp = 0.5*(az(ip,jp,k)+az(i,jp,k))
            vzjm = 0.5*(az(ip,jm,k)+az(i,jm,k))
cc          endif

          vxkp = 0.5*(ax(ip,j,kp)+ax(i,j,kp))
          vxkm = 0.5*(ax(ip,j,km)+ax(i,j,km))
          vykp = 0.5*(ay(ip,j,kp)+ay(i,j,kp))
          vykm = 0.5*(ay(ip,j,km)+ay(i,j,km))
          vzkp = 0.5*(az(ip,j,kp)+az(i,j,kp))
          vzkm = 0.5*(az(ip,j,km)+az(i,j,km))

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

        case (2)
          dhy = grid_params%dy(jg)

          vxip = 0.5*(ax(ip,j,k)+ax(ip,jp,k))
          vxim = 0.5*(ax(im,j,k)+ax(im,jp,k))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,jp,k))
          vyim = 0.5*(ay(im,j,k)+ay(im,jp,k))
          vzip = 0.5*(az(ip,j,k)+az(ip,jp,k))
          vzim = 0.5*(az(im,j,k)+az(im,jp,k))

cc          if (spoint) then
cc            !4th order
ccc FIX PARALLEL for PER BC
cc            jpp = jp+1
cc            if (j == ny) jpp = 2
cc            if (j == 0 ) jm  = ny-1
cc
cc            vxjp = cf4*ax(i,jp,k)-cf3*ax(i,jpp,k)
cc            vxjm = cf4*ax(i,j ,k)-cf3*ax(i,jm ,k)
cc            vyjp = cf4*ay(i,jp,k)-cf3*ay(i,jpp,k)
cc            vyjm = cf4*ay(i,j ,k)-cf3*ay(i,jm ,k)
cc            vzjp = cf4*az(i,jp,k)-cf3*az(i,jpp,k)
cc            vzjm = cf4*az(i,j ,k)-cf3*az(i,jm ,k)
cc
cc            tsrc = 0.5*(-cf1*nabla_v_src(i,jpp,k)
cc     .                  +cf2*nabla_v_src(i,jp ,k)
cc     .                  +cf2*nabla_v_src(i,j  ,k)
cc     .                  -cf1*nabla_v_src(i,jm ,k))
cc          else
            !2nd order
            vxjp = ax(i,jp,k)
            vxjm = ax(i,j ,k)
            vyjp = ay(i,jp,k)
            vyjm = ay(i,j ,k)
            vzjp = az(i,jp,k)
            vzjm = az(i,j ,k)

            tsrc = 0.5*(nabla_v_src(i,jp,k)
     .                 +nabla_v_src(i,j ,k))
cc          endif

          vxkp = 0.5*(ax(i,j,kp)+ax(i,jp,kp))
          vxkm = 0.5*(ax(i,j,km)+ax(i,jp,km))
          vykp = 0.5*(ay(i,j,kp)+ay(i,jp,kp))
          vykm = 0.5*(ay(i,j,km)+ay(i,jp,km))
          vzkp = 0.5*(az(i,j,kp)+az(i,jp,kp))
          vzkm = 0.5*(az(i,j,km)+az(i,jp,km))

        case (3)
          dhz = grid_params%dz(kg)

          vxip = 0.5*(ax(ip,j,k)+ax(ip,j,kp))
          vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,j,kp))
          vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
          vzip = 0.5*(az(ip,j,k)+az(ip,j,kp))
          vzim = 0.5*(az(im,j,k)+az(im,j,kp))

cc          if (spoint) then
ccc FIX PARALLEL for PER BC
cc          !4th order
cc            jpp = jp+1
cc            jmm = jm-1
cc            if (j == ny) jpp = 2
cc            if (j == 1 ) jmm = ny-1
cc
cc            vxjp = 0.5*(-cf1*(ax(i,jpp,kp)+ax(i,jpp,k))
cc     .                  +cf2*(ax(i,jp ,kp)+ax(i,jp ,k))
cc     .                  +cf2*(ax(i,j  ,kp)+ax(i,j  ,k))
cc     .                  -cf1*(ax(i,jm ,kp)+ax(i,jm ,k)))
cc            vxjm = 0.5*(-cf1*(ax(i,jp ,kp)+ax(i,jp ,k))
cc     .                  +cf2*(ax(i,j  ,kp)+ax(i,j  ,k))
cc     .                  +cf2*(ax(i,jm ,kp)+ax(i,jm ,k))
cc     .                  -cf1*(ax(i,jmm,kp)+ax(i,jmm,k)))
cc
cc            vyjp = 0.5*(-cf1*(ay(i,jpp,kp)+ay(i,jpp,k))
cc     .                  +cf2*(ay(i,jp ,kp)+ay(i,jp ,k))
cc     .                  +cf2*(ay(i,j  ,kp)+ay(i,j  ,k))
cc     .                  -cf1*(ay(i,jm ,kp)+ay(i,jm ,k)))
cc            vyjm = 0.5*(-cf1*(ay(i,jp ,kp)+ay(i,jp ,k))
cc     .                  +cf2*(ay(i,j  ,kp)+ay(i,j  ,k))
cc     .                  +cf2*(ay(i,jm ,kp)+ay(i,jm ,k))
cc     .                  -cf1*(ay(i,jmm,kp)+ay(i,jmm,k)))
cc
cc            vzjp = 0.5*(-cf1*(az(i,jpp,kp)+az(i,jpp,k))
cc     .                  +cf2*(az(i,jp ,kp)+az(i,jp ,k))
cc     .                  +cf2*(az(i,j  ,kp)+az(i,j  ,k))
cc     .                  -cf1*(az(i,jm ,kp)+az(i,jm ,k)))
cc            vzjm = 0.5*(-cf1*(az(i,jp ,kp)+az(i,jp ,k))
cc     .                  +cf2*(az(i,j  ,kp)+az(i,j  ,k))
cc     .                  +cf2*(az(i,jm ,kp)+az(i,jm ,k))
cc     .                  -cf1*(az(i,jmm,kp)+az(i,jmm,k)))
cc          else
cc          !2nd order
            vxjp = 0.5*(ax(i,jp,k)+ax(i,jp,kp))
            vxjm = 0.5*(ax(i,jm,k)+ax(i,jm,kp))
            vyjp = 0.5*(ay(i,jp,k)+ay(i,jp,kp))
            vyjm = 0.5*(ay(i,jm,k)+ay(i,jm,kp))
            vzjp = 0.5*(az(i,jp,k)+az(i,jp,kp))
            vzjm = 0.5*(az(i,jm,k)+az(i,jm,kp))
cc          endif

          vxkp = ax(i,j,kp)
          vxkm = ax(i,j,k )
          vykp = ay(i,j,kp)
          vykm = ay(i,j,k )
          vzkp = az(i,j,kp)
          vzkm = az(i,j,k )

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

        case default

          vxip = ax(ip,j,k)
          vxim = ax(im,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(im,j,k)
          vzip = az(ip,j,k)
          vzim = az(im,j,k)

cc          if (spoint) then
ccc FIX PARALLEL for PER BC
cc          !4th order
cc            jpp = jp+1
cc            jmm = jm-1
cc            if (j == ny) jpp = 2
cc            if (j == 1 ) jmm = ny-1
cc
cc            vxjp = (-cf1*(ax(i,jpp,k))
cc     .              +cf2*(ax(i,jp ,k))
cc     .              +cf2*(ax(i,j  ,k))
cc     .              -cf1*(ax(i,jm ,k)))
cc            vxjm = (-cf1*(ax(i,jp ,k))
cc     .              +cf2*(ax(i,j  ,k))
cc     .              +cf2*(ax(i,jm ,k))
cc     .              -cf1*(ax(i,jmm,k)))
cc
cc            vyjp = (-cf1*(ay(i,jpp,k))
cc     .              +cf2*(ay(i,jp ,k))
cc     .              +cf2*(ay(i,j  ,k))
cc     .              -cf1*(ay(i,jm ,k)))
cc            vyjm = (-cf1*(ay(i,jp ,k))
cc     .              +cf2*(ay(i,j  ,k))
cc     .              +cf2*(ay(i,jm ,k))
cc     .              -cf1*(ay(i,jmm,k)))
cc
cc            vzjp = (-cf1*(az(i,jpp,k))
cc     .              +cf2*(az(i,jp ,k))
cc     .              +cf2*(az(i,j  ,k))
cc     .              -cf1*(az(i,jm ,k)))
cc            vzjm = (-cf1*(az(i,jp ,k))
cc     .              +cf2*(az(i,j  ,k))
cc     .              +cf2*(az(i,jm ,k))
cc     .              -cf1*(az(i,jmm,k)))
cc          else
          !2nd order
            vxjp = ax(i,jp,k)
            vxjm = ax(i,jm,k)
            vyjp = ay(i,jp,k)
            vyjm = ay(i,jm,k)
            vzjp = az(i,jp,k)
            vzjm = az(i,jm,k)
cc          endif

          vxkp = ax(i,j,kp)
          vxkm = ax(i,j,km)
          vykp = ay(i,j,kp)
          vykm = ay(i,j,km)
          vzkp = az(i,j,kp)
          vzkm = az(i,j,km)

          tsrc = nabla_v_src(i,j,k)

        end select

        idhx = 1./dhx
        idhy = 1./dhy
        idhz = 1./dhz

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)*idhx

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)*idhx

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)*idhx

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)*idhy

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)*idhy

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)*idhy

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)*idhz

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)*idhz

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)*idhz

      ! Add geometric source

        tensor = tensor - tsrc

c     End program

      contains

c     nabla_v_src
c     #############################################################
      function nabla_v_src(i,j,k) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

        integer    :: i,j,k
        real(8)    :: tensor(3,3)

        real(8)    :: hessian(3,3,3)

c     Begin program

        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

      ! l = 1, m = 1
        tensor(1,1) =+ ax(i,j,k)* hessian(1,1,1)
     .               + ay(i,j,k)* hessian(2,1,1)
     .               + az(i,j,k)* hessian(3,1,1)

      ! l = 1, m = 2
        tensor(1,2) =+ ax(i,j,k)* hessian(1,1,2)
     .               + ay(i,j,k)* hessian(2,1,2)
     .               + az(i,j,k)* hessian(3,1,2)

      ! l = 1, m = 3
        tensor(1,3) =+ ax(i,j,k)* hessian(1,1,3)
     .               + ay(i,j,k)* hessian(2,1,3)
     .               + az(i,j,k)* hessian(3,1,3)

      ! l = 2, m = 1
        tensor(2,1) =+ ax(i,j,k)* hessian(1,2,1)
     .               + ay(i,j,k)* hessian(2,2,1)
     .               + az(i,j,k)* hessian(3,2,1)

      ! l = 2, m = 2
        tensor(2,2) =+ ax(i,j,k)* hessian(1,2,2)
     .               + ay(i,j,k)* hessian(2,2,2)
     .               + az(i,j,k)* hessian(3,2,2)

      ! l = 2, m = 3
        tensor(2,3) =+ ax(i,j,k)* hessian(1,2,3)
     .               + ay(i,j,k)* hessian(2,2,3)
     .               + az(i,j,k)* hessian(3,2,3)

      ! l = 3, m = 1
        tensor(3,1) =+ ax(i,j,k)* hessian(1,3,1)
     .               + ay(i,j,k)* hessian(2,3,1)
     .               + az(i,j,k)* hessian(3,3,1)

      ! l = 3, m = 2
        tensor(3,2) =+ ax(i,j,k)* hessian(1,3,2)
     .               + ay(i,j,k)* hessian(2,3,2)
     .               + az(i,j,k)* hessian(3,3,2)

      ! l = 3, m = 3
        tensor(3,3) =+ ax(i,j,k)* hessian(1,3,3)
     .               + ay(i,j,k)* hessian(2,3,3)
     .               + az(i,j,k)* hessian(3,3,3)

c     End program

      end function nabla_v_src

      end function fnabla_v_cov

c     dtnsr_flx_ijk
c     ###############################################################
      function dtnsr_flx_ijk(i,j,k,nx,ny,nz,igrid,alt_eom,flux)
     .         result(flx)
c     ---------------------------------------------------------------
c     Calculates fluxes for contravariant components of div(tensor) in
c     general non-orthogonal coordinates at faces. The routine flux must
c     provide the contravariant tersor components T^lm/jac at faces
c     (i+1/2,j,k),(i,j+1/2,k),(i,j,k+1/2) (which face is computed
c     depends on which flux routine is called).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8) :: flx(3)

      integer :: i,j,k,igrid,nx,ny,nz

      logical :: alt_eom

      external   flux

c     Local variables

      real(8) :: jac

c     Begin program

      call flux(i,j,k,nx,ny,nz,igrid,igrid,igrid,alt_eom
     .         ,flx(1),flx(2),flx(3),1)

c     End program

      end function dtnsr_flx_ijk

c     dtnsr_src_ijk
c     ###############################################################
      function dtnsr_src_ijk(i,j,k,nx,ny,nz,igr,alt_eom,flxx,flxy,flxz)
     .         result(src)
c     ---------------------------------------------------------------
c     Calculates source for contravariant components of div(tensor) in
c     general non-orthogonal coordinates at (i,j,k). The routines
c     flxx,flxy,flxz must provide the contravariant tersor components
c     T^lm (if flag /= 0, they actually provide T^lm/jac).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8) :: src(3)

      integer :: i,j,k,igr,nx,ny,nz

      logical :: alt_eom

      external   flxx,flxy,flxz

c     Local variables

      real(8) :: hess(3,3,3),to(3,3),jac,ijac

c     Begin program

      if (.not.no_map) then
        hess = gmetric%grid(igr)%Gamma(i,j,k,:,:,:)

        jac = gmetric%grid(igr)%jac(i,j,k)
        ijac = 1d0/jac

        call flxx(i,j,k,nx,ny,nz,igr,igr,igr,alt_eom
     .           ,to(1,1),to(1,2),to(1,3),0)
        call flxy(i,j,k,nx,ny,nz,igr,igr,igr,alt_eom
     .           ,to(2,1),to(2,2),to(2,3),0)
        call flxz(i,j,k,nx,ny,nz,igr,igr,igr,alt_eom
     .           ,to(3,1),to(3,2),to(3,3),0)

        src(1) = sum(to*hess(1,:,:))*ijac

        if (alt_eom) then
          src(2) =(sum(to*hess(2,:,:))
     .            -dot_product(to(:,2),hess(1,:,1))
     .            -dot_product(to(:,2),hess(2,:,2))
     .            -dot_product(to(:,2),hess(3,:,3)))
        else
          src(2)=sum(to*hess(2,:,:))*ijac
        endif

        src(3)=sum(to*hess(3,:,:))*ijac
      else
        src = 0d0
      endif

c     End program

      end function dtnsr_src_ijk

c     div_tensor
c     ###############################################################
      function div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,flxx,flxy,flxz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculates the contravariant components of div(tensor) at cell
c     centers in general non-orthogonal coordinates, given by:
c
c         div(T)^i = partial_j(T^ji/jac) + T^lm/jac Gamma(i)_lm)
c
c     The routines flxx,flxy,flxz must provide the contravariant
c     tensor fluxes and sources (if flag /= 0, they actually provide
c     T^lm/jac; otherwise, T^lm).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8) :: divt(3)

      integer :: i,j,k,igx,igy,igz,nx,ny,nz

      logical :: alt_eom

      external   flxx,flxy,flxz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

      real(8) :: dvol,dxx,dyy,dzz,jac,ijac

      real(8) :: tp(3,3),tm(3,3),to(3,3)

      real(8) :: msource(3)

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      igrid = igx

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

c     Fluxes

      tp(1,:) = dtnsr_flx_ijk(i ,j,k,nx,ny,nz,igx,alt_eom,flxx)
      tm(1,:) = dtnsr_flx_ijk(im,j,k,nx,ny,nz,igx,alt_eom,flxx)

      tp(2,:) = dtnsr_flx_ijk(i,j ,k,nx,ny,nz,igx,alt_eom,flxy)
      tm(2,:) = dtnsr_flx_ijk(i,jm,k,nx,ny,nz,igx,alt_eom,flxy)

      tp(3,:) = dtnsr_flx_ijk(i,j,k ,nx,ny,nz,igx,alt_eom,flxz)
      tm(3,:) = dtnsr_flx_ijk(i,j,km,nx,ny,nz,igx,alt_eom,flxz)

c     Source

      msource = dtnsr_src_ijk(i,j,k,nx,ny,nz,igx,alt_eom,flxx,flxy,flxz)

c     Assemble divergence of tensor

      divt = (tp(1,:)-tm(1,:))/dxx
     .      +(tp(2,:)-tm(2,:))/dyy
     .      +(tp(3,:)-tm(3,:))/dzz + msource

      if (alt_eom) then
        divt(2) = divt(2)/gmetric%grid(igx)%jac(i,j,k)
      endif

      !Volume factor
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dlvol(i,j,k)  !Logical volume

c     End program

      end function div_tensor

ccc     div_tensor
ccc     ###############################################################
cc      function div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                   ,tsrx,tsry,tsrz,vol) result (divt)
cc
ccc     ---------------------------------------------------------------
ccc     Calculatesthe contravariant components of div(tensor) at cell
ccc     centers in general non-orthogonal coordinates, given by:
ccc
ccc         div(T)^i = partial_j(T^ji/jac) + T^lm/jac Gamma(i)_lm)
ccc
ccc     The routines tsrx,tsry,tsrz must provide the contravariant
ccc     tersor components T^lm (if flag /= 0, they actually provide
ccc     T^lm/jac).
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      real(8)    :: divt(3)
cc
cc      integer    :: i,j,k,igx,igy,igz,nx,ny,nz
cc
cc      logical    :: alt_eom
cc
cc      external   tsrx,tsry,tsrz
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
cc
cc      real(8)    :: jac,ijac,dvol,dS1,dS2,dS3,dxx,dyy,dzz
cc
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
cc     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
cc     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o
cc
cc      real(8)    :: hess(3,3,3),msource,dx1,dx2,ll
cc
cc      logical    :: vol_wgt
cc
cc      real(8)    :: dum1,dum2,coeff
cc
ccc     Begin program
cc
cc      vol_wgt = .false.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
cc      igrid = igx
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      jac = gmetric%grid(igrid)%jac(i,j,k)
cc      ijac = 1./jac
cc
cc      dxx = grid_params%dxh(ig)
cc      dyy = grid_params%dyh(jg)
cc      dzz = grid_params%dzh(kg)
cc
cc      if (.not.no_map) then
cc        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc      endif
cc
cc      call tsrx(i ,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11p,t12p,t13p, 1)
cc      call tsrx(im,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11m,t12m,t13m,-1)
cc      if (.not.no_map)
cc     .    call tsrx(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11o,t12o,t13o,0)
cc
cc      call tsry(i,j ,k,nx,ny,nz,igx,igy,igz,alt_eom,t21p,t22p,t23p, 1)
cc      call tsry(i,jm,k,nx,ny,nz,igx,igy,igz,alt_eom,t21m,t22m,t23m,-1)
cc      if (.not.no_map)
cc     .    call tsry(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t21o,t22o,t23o,0)
cc
cc      call tsrz(i,j,k ,nx,ny,nz,igx,igy,igz,alt_eom,t31p,t32p,t33p, 1)
cc      call tsrz(i,j,km,nx,ny,nz,igx,igy,igz,alt_eom,t31m,t32m,t33m,-1)
cc      if (.not.no_map)
cc     .    call tsrz(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t31o,t32o,t33o,0)
cc
cc      msource = 0d0
cc
cc      !Component 1
cc      if (.not.no_map) then
cc        msource =  (t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
cc     .             +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
cc     .             +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
cc     .             *ijac
cc      endif
cc
cc      divt(1) =  (t11p - t11m)/dxx
cc     .          +(t21p - t21m)/dyy
cc     .          +(t31p - t31m)/dzz + msource
cc
cc      !Component 2
cc      coeff = 1d0
cc
cc      if (.not.no_map) then
cc        if (alt_eom) then
cc
cc          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
cc     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
cc     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
cc     .              -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
cc     .              -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
cc     .              -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))*ijac
cc
cc          coeff = ijac
cc
cc        else
cc          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
cc     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
cc     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
cc     .              *ijac
cc        endif
cc
cc      endif
cc
cc      divt(2) = ( (t12p - t12m)/dxx
cc     .           +(t22p - t22m)/dyy
cc     .           +(t32p - t32m)/dzz)*coeff  + msource
cc
cc      !Component 3
cc      if (.not.no_map) then
cc        msource = (t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
cc     .            +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
cc     .            +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
cc     .            *ijac
cc      endif
cc
cc      divt(3) =  (t13p - t13m)/dxx
cc     .          +(t23p - t23m)/dyy
cc     .          +(t33p - t33m)/dzz + msource
cc
cc      !Volume factor
cc      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function div_tensor

c     dtnsr_cov_flx_ijk
c     ###############################################################
      function dtnsr_cov_flx_ijk(i,j,k,nx,ny,nz,igrid,flux)
     .         result(flx)
c     ---------------------------------------------------------------
c     Calculates fluxes for contravariant components of div(tensor) in
c     general non-orthogonal coordinates at faces. The routine flux must
c     provide the contravariant tersor components T^lm/jac at faces
c     (i+1/2,j,k),(i,j+1/2,k),(i,j,k+1/2) (which face is computed
c     depends on which flux routine is called).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8) :: flx(3)

      integer :: i,j,k,igrid,nx,ny,nz

      external   flux

c     Local variables

      call flux(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .         ,flx(1),flx(2),flx(3),1)

c     End program

      end function dtnsr_cov_flx_ijk

c     dtnsr_cov_src_ijk
c     ###############################################################
      function dtnsr_cov_src_ijk(i,j,k,nx,ny,nz,igr,flxx,flxy,flxz)
     .         result(src)
c     ---------------------------------------------------------------
c     Calculates source for contravariant components of div(tensor) in
c     general non-orthogonal coordinates at (i,j,k). The routines
c     flxx,flxy,flxz must provide the contravariant tersor components
c     T^lm (if flag /= 0, they actually provide T^lm/jac).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8) :: src(3)

      integer :: i,j,k,igr,nx,ny,nz

      external   flxx,flxy,flxz

c     Local variables

      real(8) :: hess(3,3,3),to(3,3)
      real(8) :: t11o,t12o,t13o
     .          ,t21o,t22o,t23o
     .          ,t31o,t32o,t33o

c     Begin program

      if (.not.no_map) then
        hess = gmetric%grid(igr)%Gamma(i,j,k,:,:,:)

cc        call flxx(i,j,k,nx,ny,nz,igr,igr,igr,t11o,t12o,t13o,0)
cc        call flxy(i,j,k,nx,ny,nz,igr,igr,igr,t21o,t22o,t23o,0)
cc        call flxz(i,j,k,nx,ny,nz,igr,igr,igr,t31o,t32o,t33o,0)
cc
cc        src(1) =-(t11o*hess(1,1,1)+t12o*hess(2,1,1)+t13o*hess(3,1,1)
cc     .           +t21o*hess(1,2,1)+t22o*hess(2,2,1)+t23o*hess(3,2,1)
cc     .           +t31o*hess(1,3,1)+t32o*hess(2,3,1)+t33o*hess(3,3,1))
cc        src(2) =-(t11o*hess(1,1,2)+t12o*hess(2,1,2)+t13o*hess(3,1,2)
cc     .           +t21o*hess(1,2,2)+t22o*hess(2,2,2)+t23o*hess(3,2,2)
cc     .           +t31o*hess(1,3,2)+t32o*hess(2,3,2)+t33o*hess(3,3,2))
cc        src(3) =-(t11o*hess(1,1,3)+t12o*hess(2,1,3)+t13o*hess(3,1,3)
cc     .           +t21o*hess(1,2,3)+t22o*hess(2,2,3)+t23o*hess(3,2,3)
cc     .           +t31o*hess(1,3,3)+t32o*hess(2,3,3)+t33o*hess(3,3,3))

        call flxx(i,j,k,nx,ny,nz,igr,igr,igr,to(1,1),to(1,2),to(1,3),0)
        call flxy(i,j,k,nx,ny,nz,igr,igr,igr,to(2,1),to(2,2),to(2,3),0)
        call flxz(i,j,k,nx,ny,nz,igr,igr,igr,to(3,1),to(3,2),to(3,3),0)

        to = transpose(to)
        src(1) =-sum(to*hess(:,:,1))
        src(2) =-sum(to*hess(:,:,2))
        src(3) =-sum(to*hess(:,:,3))
      else
        src = 0d0
      endif

c     End program

      end function dtnsr_cov_src_ijk

c     div_tensor_cov
c     ###############################################################
      function div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                       ,flxx,flxy,flxz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculates covariant components of div(tensor) at cell centers in
c     general non-orthogonal coordinates, which is given by:
c
c         div(T)_i = jac^(-1)*[ partial_j(T^j_i)-T^l_m Gamma(m)_il) ]
c
c     The routines flxx,flxy,flxz must provide mixed representation 
c     tensor components T^l_m.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8) :: divt(3)

      integer :: i,j,k,igx,igy,igz,nx,ny,nz

      external   flxx,flxy,flxz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

      real(8) :: dxx,dyy,dzz

      real(8) :: tp(3,3),tm(3,3),to(3,3)

      real(8) :: hess(3,3,3),msource(3)

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      igrid = igx

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

c     Fluxes

      tp(1,:) = dtnsr_cov_flx_ijk(i ,j,k,nx,ny,nz,igx,flxx)
      tm(1,:) = dtnsr_cov_flx_ijk(im,j,k,nx,ny,nz,igx,flxx)

      tp(2,:) = dtnsr_cov_flx_ijk(i,j ,k,nx,ny,nz,igx,flxy)
      tm(2,:) = dtnsr_cov_flx_ijk(i,jm,k,nx,ny,nz,igx,flxy)

      tp(3,:) = dtnsr_cov_flx_ijk(i,j,k ,nx,ny,nz,igx,flxz)
      tm(3,:) = dtnsr_cov_flx_ijk(i,j,km,nx,ny,nz,igx,flxz)

c     Source

      msource = dtnsr_cov_src_ijk(i,j,k,nx,ny,nz,igx,flxx,flxy,flxz)

c     Assemble divergence of tensor

      divt = (tp(1,:)-tm(1,:))/dxx
     .      +(tp(2,:)-tm(2,:))/dyy
     .      +(tp(3,:)-tm(3,:))/dzz + msource

      divt = divt/gmetric%grid(igrid)%jac(i,j,k)

      !Volume factor
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)

c     End program

      end function div_tensor_cov

c     veclap_mesh
c     ###############################################################
      function veclap_mesh(igrid,cnv,diff) result(cnv_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cnv(0:,0:,0:,:)
     .          ,cnv_out(0:size(cnv,1)-1
     $                  ,0:size(cnv,2)-1
     $                  ,0:size(cnv,3)-1,3)

      real(8),optional :: diff(0:,0:,0:)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cnv_out(i,j,k,:)=veclap(i,j,k,nx,ny,nz
     .                             ,igrid,igrid,igrid,cnv
     .                             ,diff=diff,alteom=alt__eom())
          enddo
        enddo
      enddo

c     End program

      end function veclap_mesh

c     veclap
c     ###############################################################
      function veclap(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .               ,diff,alteom,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*div(diff grad(vector)) at cell centers in general
c     non-orthog.  coordinates. Vector is assumed in contravariant
c     representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target,optional :: diff(0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol,alteom

c     Local variables

      integer    :: icomp

      logical    :: vol_wgt,alt_eom

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      alt_eom = .false.
      if (PRESENT(alteom)) alt_eom = alteom

      if (no_map) then
        do icomp=1,3
          vlap(icomp)=lap(i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vfield(:,:,:,icomp),dff=diff
     .                   ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        if (PRESENT(diff)) coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
        if (PRESENT(diff)) nullify(coef)
      endif

c     End program

      end function veclap

c$$$c     veclap_diff
c$$$c     ###############################################################
c$$$      function veclap_diff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
c$$$     .                    ,diff,alteom,vol) result (vlap)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates dvol*div(diff grad(vector)) at cell centers in general
c$$$c     non-orthog.  coordinates. Vector is assumed in contravariant
c$$$c     representation.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
c$$$     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)
c$$$
c$$$      real(8)    :: vlap(3)
c$$$
c$$$      logical,optional,intent(IN) :: vol,alteom
c$$$
c$$$c     Local variables
c$$$
c$$$      integer    :: icomp
c$$$
c$$$      logical    :: vol_wgt,alt_eom
c$$$
c$$$c     Begin program
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$      alt_eom = .false.
c$$$      if (PRESENT(alteom)) alt_eom = alteom
c$$$
c$$$      if (no_map) then
c$$$        do icomp=1,3
c$$$          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                           ,vfield(:,:,:,icomp),diff
c$$$     .                           ,vol=vol_wgt)
c$$$        enddo
c$$$      else
c$$$        vec  => vfield !Pointer passed to nabtensor routines
c$$$        coef => diff   !Pointer passed to nabtensor routines
c$$$        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
c$$$     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
c$$$     $                   ,vol=vol_wgt)
c$$$        nullify(coef)
c$$$      endif
c$$$
c$$$c     End program
c$$$
c$$$      end function veclap_diff
c$$$
c$$$c     veclap_ndiff
c$$$c     ###############################################################
c$$$      function veclap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
c$$$     .                     ,alteom,vol) result (vlap)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates dvol*lap(vector) at cell centers in contravariant
c$$$c     representation. Vector vfield is contravariant.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
c$$$
c$$$      real(8)    :: vlap(3)
c$$$
c$$$      logical,optional,intent(IN) :: vol,alteom
c$$$
c$$$c     Local variables
c$$$
c$$$      integer    :: icomp
c$$$
c$$$      logical    :: vol_wgt,alt_eom
c$$$
c$$$c     Begin program
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$      alt_eom = .false.
c$$$      if (PRESENT(alteom)) alt_eom = alteom
c$$$
c$$$      if (no_map) then
c$$$        do icomp=1,3
c$$$          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                           ,vfield(:,:,:,icomp)
c$$$     .                           ,vol=vol_wgt)
c$$$        enddo
c$$$      else
c$$$        vec  => vfield !Pointer passed to nabtensor routines
c$$$        nullify(coef)
c$$$        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
c$$$     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
c$$$     .                   ,vol=vol_wgt)
c$$$      endif
c$$$
c$$$c     End program
c$$$
c$$$      end function veclap_ndiff

c     nabtensor_x
c     #############################################################
      subroutine nabtensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for nabla(vec)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: x,y,z,jac,jac0,jacp,ijac
        real(8)    :: nabla_v_cnv(3,3),dd

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
        ijac = 1d0/jac
cc        if (isSP2(i+1,igx)) ijac = 0d0

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),1,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*nabla_v_cnv(1,1)
        t12 = dd*nabla_v_cnv(1,2)
        t13 = dd*nabla_v_cnv(1,3)

        if (flag /= 0) then
          t11 = t11*ijac
          if (.not.alteom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine nabtensor_x

c     nabtensor_y
c     #############################################################
      subroutine nabtensor_y(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v_cnv(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        ijac = 1d0/jac

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),2,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j ,k))
        else
          dd = 1d0
        endif

        t21 = dd*nabla_v_cnv(2,1)
        t22 = dd*nabla_v_cnv(2,2)
        t23 = dd*nabla_v_cnv(2,3)

        if (flag /= 0) then
          t21 = t21*ijac
          if (.not.alteom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine nabtensor_y

c     nabtensor_z
c     #############################################################
      subroutine nabtensor_z(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v_cnv(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        ijac = 1d0/jac

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),3,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k ))
        else
          dd = 1d0
        endif

        t31 = dd*nabla_v_cnv(3,1)
        t32 = dd*nabla_v_cnv(3,2)
        t33 = dd*nabla_v_cnv(3,3)

        if (flag /= 0) then
          t31 = t31*ijac
          if (.not.alteom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine nabtensor_z

c$$$c     veclap_cov_ndiff_mesh
c$$$c     ###############################################################
c$$$      function veclap_cov_ndiff_mesh(igrid,cov) result(cov_out)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates curl(A) in general non-orthogonal coordinates. Is
c$$$c     boundary-ready (can extrapolate). The vector A is both covariant
c$$$c     and contravariant, and returns the covariant component "comp".
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: igrid
c$$$
c$$$      real(8) :: cov(0:,0:,0:,:)
c$$$     .          ,cov_out(0:size(cov,1)-1
c$$$     $                  ,0:size(cov,2)-1
c$$$     $                  ,0:size(cov,3)-1,3)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: nx,ny,nz,i,j,k
c$$$
c$$$c     Begin program
c$$$
c$$$      nx = grid_params%nxv(igrid)
c$$$      ny = grid_params%nyv(igrid)
c$$$      nz = grid_params%nzv(igrid)
c$$$
c$$$      do k=1,nz
c$$$        do j=1,ny
c$$$          do i=1,nx
c$$$            cov_out(i,j,k,:)=veclap_cov_ndiff(i,j,k,nx,ny,nz
c$$$     $                                       ,igrid,igrid,igrid,cov)
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     End program
c$$$
c$$$      end function veclap_cov_ndiff_mesh
c$$$
c$$$c     veclap_cov_diff_mesh
c$$$c     ###############################################################
c$$$      function veclap_cov_diff_mesh(igrid,diff,cov) result(cov_out)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates curl(A) in general non-orthogonal coordinates. Is
c$$$c     boundary-ready (can extrapolate). The vector A is both covariant
c$$$c     and contravariant, and returns the covariant component "comp".
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: igrid
c$$$
c$$$      real(8) :: cov(0:,0:,0:,:)
c$$$     .          ,diff(0:,0:,0:)
c$$$     .          ,cov_out(0:size(cov,1)-1
c$$$     $                  ,0:size(cov,2)-1
c$$$     $                  ,0:size(cov,3)-1,3)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: nx,ny,nz,i,j,k
c$$$
c$$$c     Begin program
c$$$
c$$$      nx = grid_params%nxv(igrid)
c$$$      ny = grid_params%nyv(igrid)
c$$$      nz = grid_params%nzv(igrid)
c$$$
c$$$      do k=1,nz
c$$$        do j=1,ny
c$$$          do i=1,nx
c$$$            cov_out(i,j,k,:)=veclap_cov_diff(i,j,k,nx,ny,nz
c$$$     $                                      ,igrid,igrid,igrid,cov,diff)
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     End program
c$$$
c$$$      end function veclap_cov_diff_mesh

c     veclap_cov_mesh
c     ###############################################################
      function veclap_cov_mesh(igrid,cov,diff) result(cov_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is both covariant
c     and contravariant, and returns the covariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cov_out(0:size(cov,1)-1
     $                  ,0:size(cov,2)-1
     $                  ,0:size(cov,3)-1,3)

      real(8),optional :: diff(0:,0:,0:)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cov_out(i,j,k,:)=veclap_cov(i,j,k,nx,ny,nz
     $                                 ,igrid,igrid,igrid,cov,diff=diff)
          enddo
        enddo
      enddo

c     End program

      end function veclap_cov_mesh

c     veclap_cov
c     ###############################################################
      function veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,vfield,diff,vol)
     $         result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of dvol*lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target,optional :: diff(0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      if (no_map) then
        do icomp=1,3
          vlap(icomp)=lap(i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vfield(:,:,:,icomp),dff=diff
     .                   ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        if (PRESENT(diff)) coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
        if (PRESENT(diff)) nullify(coef)
      endif

c     End program

      end function veclap_cov

c$$$c     veclap_cov_ndiff
c$$$c     ###############################################################
c$$$      function veclap_cov_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
c$$$     .                         ,vol) result (vlap)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates covariant components of dvol*lap(vfield) at cell
c$$$c     centers in general non-orthogonal coordinates. Vector vfield
c$$$c     in covariant representation.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
c$$$
c$$$      real(8)    :: vlap(3)
c$$$
c$$$      logical,optional,intent(IN) :: vol
c$$$
c$$$c     Local variables
c$$$
c$$$      integer    :: icomp
c$$$
c$$$      logical    :: vol_wgt
c$$$
c$$$c     Begin program
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$      if (no_map) then
c$$$        do icomp=1,3
c$$$          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                         ,vfield(:,:,:,icomp),vol=vol_wgt)
c$$$        enddo
c$$$      else
c$$$        vec => vfield !Pointer passed to nabtensor routines
c$$$        nullify(coef)
c$$$        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
c$$$     $             ,vol=vol_wgt)
c$$$      endif
c$$$
c$$$c     End program
c$$$
c$$$      end function veclap_cov_ndiff

c$$$c     veclap_cov_diff
c$$$c     ###############################################################
c$$$      function veclap_cov_diff(i,j,k,nx,ny,nz,igx,igy,igz,vfield,diff
c$$$     .                         ,vol) result (vlap)
c$$$
c$$$c     ---------------------------------------------------------------
c$$$c     Calculates covariant components of dvol*lap(vfield) at cell
c$$$c     centers in general non-orthogonal coordinates. Vector vfield
c$$$c     in covariant representation.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none           !For safe fortran
c$$$
c$$$c     Call variables
c$$$
c$$$      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$
c$$$      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
c$$$     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)
c$$$
c$$$      real(8)    :: vlap(3)
c$$$
c$$$      logical,optional,intent(IN) :: vol
c$$$
c$$$c     Local variables
c$$$
c$$$      integer    :: icomp
c$$$
c$$$      logical    :: vol_wgt
c$$$
c$$$c     Begin program
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$      if (no_map) then
c$$$        do icomp=1,3
c$$$          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                         ,vfield(:,:,:,icomp),diff
c$$$     .                         ,vol=vol_wgt)
c$$$        enddo
c$$$      else
c$$$        vec  => vfield !Pointer passed to nabtensor routines
c$$$        coef => diff   !Pointer passed to nabtensor routines
c$$$        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
c$$$     $             ,vol=vol_wgt)
c$$$        nullify(coef)
c$$$      endif
c$$$
c$$$c     End program
c$$$
c$$$      end function veclap_cov_diff

c     nabtensor_cov_x
c     #############################################################
      subroutine nabtensor_cov_x(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t11-t13 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: x,y,z,jac,jac0,jacp
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),1)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*( gsuper(1,1)*nabla_v(1,1)
     .            +gsuper(1,2)*nabla_v(2,1)
     .            +gsuper(1,3)*nabla_v(3,1) )

        t12 = dd*( gsuper(1,1)*nabla_v(1,2)
     .            +gsuper(1,2)*nabla_v(2,2)
     .            +gsuper(1,3)*nabla_v(3,2) )

        t13 = dd*( gsuper(1,1)*nabla_v(1,3)
     .            +gsuper(1,2)*nabla_v(2,3)
     .            +gsuper(1,3)*nabla_v(3,3) )

c     End program

      end subroutine nabtensor_cov_x

c     nabtensor_cov_y
c     #############################################################
      subroutine nabtensor_cov_y(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t21-t23 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: x,y,z,jac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),2)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t21 = dd*( gsuper(2,1)*nabla_v(1,1)
     .            +gsuper(2,2)*nabla_v(2,1)
     .            +gsuper(2,3)*nabla_v(3,1) )

        t22 = dd*( gsuper(2,1)*nabla_v(1,2)
     .            +gsuper(2,2)*nabla_v(2,2)
     .            +gsuper(2,3)*nabla_v(3,2) )

        t23 = dd*( gsuper(2,1)*nabla_v(1,3)
     .            +gsuper(2,2)*nabla_v(2,3)
     .            +gsuper(2,3)*nabla_v(3,3) )

c     End program

      end subroutine nabtensor_cov_y

c     nabtensor_cov_z
c     #############################################################
      subroutine nabtensor_cov_z(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t31-t33 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: x,y,z,jac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),3)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t31 = dd*( gsuper(3,1)*nabla_v(1,1)
     .            +gsuper(3,2)*nabla_v(2,1)
     .            +gsuper(3,3)*nabla_v(3,1) )

        t32 = dd*( gsuper(3,1)*nabla_v(1,2)
     .            +gsuper(3,2)*nabla_v(2,2)
     .            +gsuper(3,3)*nabla_v(3,2) )

        t33 = dd*( gsuper(3,1)*nabla_v(1,3)
     .            +gsuper(3,2)*nabla_v(2,3)
     .            +gsuper(3,3)*nabla_v(3,3) )

c     End program

      end subroutine nabtensor_cov_z

c     curlcurl_cov_mesh
c     ###############################################################
      function curlcurl_cov_mesh(igrid,cov) result(cov_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is both covariant
c     and contravariant, and returns the covariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cov_out(0:size(cov,1)-1
     $                  ,0:size(cov,2)-1
     $                  ,0:size(cov,3)-1,3)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cov_out(i,j,k,:)=curlcurl_cov_ijk(i,j,k,nx,ny,nz
     $                                       ,igrid,igrid,igrid,cov)
          enddo
        enddo
      enddo

c     End program

      end function curlcurl_cov_mesh

c     curlcurl_cov_ijk
c     ###############################################################
      function curlcurl_cov_ijk(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol)
     .         result(cov)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal coordinates,
c     preserving the SPD property. The vector A is both in covariant and
c     contravariant representations, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cov(3)

      real(8),target :: acov(0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional,intent(IN) :: vol

c     Local variables

c     Begin program

      vec => acov

      cov = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                ,curlcurl_cov_x,curlcurl_cov_y,curlcurl_cov_z
     .                ,vol=vol)

      nullify(vec)

c     End program

      end function curlcurl_cov_ijk

c     curlcurl_cov_x
c     #############################################################
      subroutine curlcurl_cov_x(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t11-t13 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),1)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        nabla_v =-nabla_v + transpose(nabla_v)

        t11 =( gsuper(1,1)*nabla_v(1,1)
     .        +gsuper(1,2)*nabla_v(2,1)
     .        +gsuper(1,3)*nabla_v(3,1) )

        t12 =( gsuper(1,1)*nabla_v(1,2)
     .        +gsuper(1,2)*nabla_v(2,2)
     .        +gsuper(1,3)*nabla_v(3,2) )

        t13 =( gsuper(1,1)*nabla_v(1,3)
     .        +gsuper(1,2)*nabla_v(2,3)
     .        +gsuper(1,3)*nabla_v(3,3) )

c     End program

      end subroutine curlcurl_cov_x

c     curlcurl_cov_y
c     #############################################################
      subroutine curlcurl_cov_y(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t21-t23 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),2)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        nabla_v =-nabla_v + transpose(nabla_v)

        t21 =( gsuper(2,1)*nabla_v(1,1)
     .        +gsuper(2,2)*nabla_v(2,1)
     .        +gsuper(2,3)*nabla_v(3,1) )

        t22 =( gsuper(2,1)*nabla_v(1,2)
     .        +gsuper(2,2)*nabla_v(2,2)
     .        +gsuper(2,3)*nabla_v(3,2) )

        t23 =( gsuper(2,1)*nabla_v(1,3)
     .        +gsuper(2,2)*nabla_v(2,3)
     .        +gsuper(2,3)*nabla_v(3,3) )

c     End program

      end subroutine curlcurl_cov_y

c     curlcurl_cov_z
c     #############################################################
      subroutine curlcurl_cov_z(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t31-t33 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),3)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        nabla_v = -nabla_v + transpose(nabla_v)

        t31 =( gsuper(3,1)*nabla_v(1,1)
     .        +gsuper(3,2)*nabla_v(2,1)
     .        +gsuper(3,3)*nabla_v(3,1) )

        t32 =( gsuper(3,1)*nabla_v(1,2)
     .        +gsuper(3,2)*nabla_v(2,2)
     .        +gsuper(3,3)*nabla_v(3,2) )

        t33 =( gsuper(3,1)*nabla_v(1,3)
     .        +gsuper(3,2)*nabla_v(2,3)
     .        +gsuper(3,3)*nabla_v(3,3) )

c     End program

      end subroutine curlcurl_cov_z

c     curlcurl_mesh
c     ###############################################################
      function curlcurl_mesh(igrid,cnv,diff) result(cnv_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cnv(0:,0:,0:,:)
     .          ,cnv_out(0:size(cnv,1)-1
     $                  ,0:size(cnv,2)-1
     $                  ,0:size(cnv,3)-1,3)

      real(8), optional :: diff(0:size(cnv,1)-1
     $                         ,0:size(cnv,2)-1
     $                         ,0:size(cnv,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cnv_out(i,j,k,:)=curlcurl_ijk(i,j,k,nx,ny,nz
     $                                   ,igrid,igrid,igrid,cnv
     $                                   ,diff=diff)
          enddo
        enddo
      enddo

c     End program

      end function curlcurl_mesh

c     curlcurl_ijk
c     ###############################################################
      function curlcurl_ijk(i,j,k,nx,ny,nz,igx,igy,igz,acnv,diff,vol)
     $         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     contravariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cnv(3)

      real(8),target :: acnv(0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target,optional :: diff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

c     Begin program

c$$$      cnv=-veclap (i,j,k,nx,ny,nz,igx,igy,igz,acnv,vol=vol)
c$$$     .    +gradDiv(i,j,k,nx,ny,nz,igx,igy,igz,acnv,vol=vol,cnv=.true.)

      vec1 => acnv
      if (PRESENT(diff)) then
         coef => diff
      endif

      cnv = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt__eom()
     .                ,curlcurl_x,curlcurl_y,curlcurl_z
     .                ,vol=vol)

      if (associated(coef)) nullify(coef)
      nullify(vec1)

c     End program

      end function curlcurl_ijk

c     curlcurl_x
c     #############################################################
      subroutine curlcurl_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                     ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for curl(curl), computed
c     as div(grad - transpose(grad)). In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate curl(curl) in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
        real(8)    :: nabla_v_cnv(3,3)
        real(8)    :: jac,ijac

c     Begin program

        if (.not.associated(vec1)) then
          call pstop('curlcurl_x','Pointers not associated')
        endif

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),1,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        nabla_v_cnv = nabla_v_cnv - transpose(nabla_v_cnv)

        if (associated(coef)) nabla_v_cnv = nabla_v_cnv
     $                                     *0.5*(coef(ip,j,k)
     $                                          +coef(i ,j,k))

        t11 = -nabla_v_cnv(1,1)
        t12 = -nabla_v_cnv(1,2)
        t13 = -nabla_v_cnv(1,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine curlcurl_x

c     curlcurl_y
c     #############################################################
      subroutine curlcurl_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                     ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for curl(curl). In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate curl(curl) in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
        real(8)    :: nabla_v_cnv(3,3)
        real(8)    :: jac,ijac

c     Begin program

        if (.not.associated(vec1)) then
          call pstop('curlcurl_y','Pointers not associated')
        endif

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))

c       Diffusive part

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),2,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        nabla_v_cnv = nabla_v_cnv - transpose(nabla_v_cnv)

        if (associated(coef)) nabla_v_cnv = nabla_v_cnv
     $                                     *0.5*(coef(i,jp,k)
     $                                          +coef(i,j ,k))

        t21 = -nabla_v_cnv(2,1)
        t22 = -nabla_v_cnv(2,2)
        t23 = -nabla_v_cnv(2,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine curlcurl_y

c     curlcurl_z
c     #############################################################
      subroutine curlcurl_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                     ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for curl(curl). In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate curl(curl) in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
        real(8)    :: nabla_v_cnv(3,3)
        real(8)    :: jac,ijac

c     Begin program

        if (.not.associated(vec1)) then
          call pstop('curlcurl_z','Pointers not associated')
        endif

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))

c       Diffusive part

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),3,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        nabla_v_cnv = nabla_v_cnv - transpose(nabla_v_cnv)

        if (associated(coef)) nabla_v_cnv = nabla_v_cnv
     $                                     *0.5*(coef(i,j,kp)
     $                                          +coef(i,j,k ))

        t31 = -nabla_v_cnv(3,1)
        t32 = -nabla_v_cnv(3,2)
        t33 = -nabla_v_cnv(3,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine curlcurl_z

      end module grid_tensor_ops

c grid_advec_ops
c ######################################################################
      module grid_advec_ops

        use grid_create

        use grid_tensor_ops

        integer :: adv_meth

        INTERFACE inter
          module procedure in2face_scl,in2face_vec
        end INTERFACE

      contains

ccc     advec
ccc     ###############################################################
cc      function advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,phi,vol
cc     .              ,upwind)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of vector field (v.phi) at cell centers in
ccc     general non-orthogonal geometry.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),advec
cc
cc      logical,optional :: vol,upwind
cc
ccc     Local variables
cc
cc      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
cc      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jac
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc      logical    :: vw,upw
cc
ccc     Begin program
cc
cc      if (PRESENT(vol)) then
cc        vw = vol
cc      else
cc        vw = .false.
cc      endif
cc
cc      if (PRESENT(upwind)) then
cc        upw = upwind
cc      else
cc        upw = .false.
cc      endif
cc
cc      igrid = igx
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      dxx = dxh(ig)
cc      dyy = dyh(jg)
cc      dzz = dzh(kg)
cc
cc      jac   = gmetric%grid(igrid)%jac(i ,j,k)
cccc      jacip = gmetric%grid(igrid)%jac(ip,j,k)
cccc      jacim = gmetric%grid(igrid)%jac(im,j,k)
cc
cc      if (upw) then
cc        flxip = 0.25*(
cc     .         (    (v1(i,j,k)+v1(ip,j,k))
cc     .          +abs(v1(i,j,k)+v1(ip,j,k)) )*phi(i ,j,k)
cc     .        +(    (v1(i,j,k)+v1(ip,j,k))          
cc     .          -abs(v1(i,j,k)+v1(ip,j,k)) )*phi(ip,j,k)) 
cc
cc        if (isSP(i,j,k,igx,igy,igz)) then
cc          flxim = 0d0
cc        else
cc          flxim = 0.25*(
cc     .         (    (v1(i,j,k)+v1(im,j,k))
cc     .          +abs(v1(i,j,k)+v1(im,j,k)) )*phi(im,j,k)
cc     .        +(    (v1(i,j,k)+v1(im,j,k))          
cc     .          -abs(v1(i,j,k)+v1(im,j,k)) )*phi(i ,j,k))
cc        endif
cc
cc        flxjp = 0.25*(
cc     .         (    (v2(i,j,k)+v2(i,jp,k))
cc     .          +abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,j ,k)
cc     .        +(    (v2(i,j,k)+v2(i,jp,k))          
cc     .          -abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,jp,k))
cc        flxjm = 0.25*(
cc     .         (    (v2(i,j,k)+v2(i,jm,k))
cc     .          +abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,jm,k)
cc     .        +(    (v2(i,j,k)+v2(i,jm,k))          
cc     .          -abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,j ,k))
cc
cc        flxkp = 0.25*(
cc     .         (    (v3(i,j,k)+v3(i,j,kp))
cc     .          +abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,k )
cc     .        +(    (v3(i,j,k)+v3(i,j,kp))             
cc     .          -abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,kp))
cc        flxkm = 0.25*(
cc     .         (    (v3(i,j,k)+v3(i,j,km))
cc     .          +abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,km)
cc     .        +(    (v3(i,j,k)+v3(i,j,km))             
cc     .          -abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,k ))
cc
cc      else
cccc        flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cccc        if (isSP(i,j,k,igx,igy,igz)) then
cccc          flxim = 0d0
cccc        else
cccc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cccc        endif
cc
cccc        if (isSP(i,j,k,igx,igy,igz)) then
cccc          flxip = 0.5*(v1(ip,j,k)*phi(ip,j,k) + v1(i,j,k)*phi(i,j,k))
cccc          flxim = 0d0
cccc        elseif (isSP(i-1,j,k,igx,igy,igz)) then
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k) + v1(i,j,k)*phi(ip,j,k))
cccc          flxim = 0.5*(v1(im,j,k)*phi(im,j,k) + v1(i,j,k)*phi(i ,j,k))
cccc        else
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cccc        endif
cc
cc        flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
cc        flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
cc
cc        flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
cc        flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
cc      endif
cc
cc      advec =( (flxip - flxim)/dxx
cc     .       + (flxjp - flxjm)/dyy
cc     .       + (flxkp - flxkm)/dzz )/jac
cc      
cc      if (vw) advec = advec*gmetric%grid(igrid)%dvol(i,j,k)
cc
ccc     End 
cc
cc      end function advec

c     scl_advc_flx_x
c     ###############################################################
      function scl_advc_flx_x(i,j,k,nx,ny,nz,igr,v,arr
     .                  ,method,sp,zip_vel) result(flxip)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_1 of scalar arr at face
c     (i+1/2,j,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see inter)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxip

      logical,optional :: sp,zip_vel

c     Local variables

      integer :: im,ip,ipp,igm,igp,igpp,ig,jg,kg

      real(8) :: jacip,jacim,jac,qip,vxp

      logical :: spoint,bc_safe,zv,at_bdry

c     Begin program

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = .false.
      endif

      if (PRESENT(zip_vel)) then
        zv = zip_vel
      else
        zv = .true.
      endif

      bc_safe = .not.(method > 2 .and. method <= 6)

c     Update grid indices

      call getMGmap(i,j,k,igr,igr,igr,ig,jg,kg)

      ip  = i+1
      im  = max(i-1,0)
      ipp = min(ip+1,nx+1)

      igp  = ig + (ip-i)
      igm  = ig + (im-i)
      igpp = ig + (ipp-i)

c     Grid info

      jac = gmetric%grid(igr)%jac(i,j,k)

      if (spoint) then
        jacip = gmetric%grid(igr)%jac(ip,j,k)
      else
        jacip = jac
      endif

c     Find X ZIP flux (used even if method /= ZIP)

      flxip = 0d0

      at_bdry = isBdry(i,igr,2) .or. isBdry(i+1,igr,1)

      if (.not.at_bdry) then               !Not at boundary
        flxip = 0.25*(v(ip,j,k)*arr(i ,j,k)/jacip
     .              + v(i ,j,k)*arr(ip,j,k)/jac  )*(jac+jacip)
      else
        if (   (v(ip,j,k) == -v(i,j,k))
     .     .or.(isSP2(i+1,igr,ibc=1).and.spoint)) then
          flxip = 0d0    !Shut off flux for vn=0 BC
        else
          flxip = 0.5*(v(ip,j,k)*arr(i,j,k)+v(i,j,k)*arr(ip,j,k))
        endif
      endif

c     Other advective methods

      if (method /= 2) then

        !Find velocities at faces from ZIP fluxes
        if (zv) then
          vxp = 2.*flxip/(arr(ip,j,k)+arr(i,j,k)+1d-12)
        else
          vxp = 0.5*(v(ip,j,k)+v(i,j,k))
        endif

        !Find advected quantity at faces
        if ((i == nx .or. i == 0).and.(.not.bc_safe)) then
          qip=in2face_scl(grid_params%xx(igpp),arr(ipp,j,k)
     .                   ,grid_params%xx(igp ),arr(ip ,j,k)
     .                   ,grid_params%xx(ig  ),arr(i  ,j,k)
     .                   ,grid_params%xx(igm ),arr(im ,j,k)
     .                   ,vxp,1)
        else
          qip=in2face_scl(grid_params%xx(igpp),arr(ipp,j,k)
     .                   ,grid_params%xx(igp ),arr(ip ,j,k)
     .                   ,grid_params%xx(ig  ),arr(i  ,j,k)
     .                   ,grid_params%xx(igm ),arr(im ,j,k)
     .                   ,vxp,method)
        endif

        !Setup fluxes (ZIP at boundaries by default for conservation)
        if (.not.(at_bdry).or.bc_safe) flxip=vxp*qip

      endif

c     End program

      end function scl_advc_flx_x

c     scl_advc_flx_y
c     ###############################################################
      function scl_advc_flx_y(i,j,k,nx,ny,nz,igr,v,arr
     .                       ,method,zip_vel) result(flxjp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_2 of scalar arr at face
c     (i,j+1/2,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see inter)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxjp

      logical,optional :: zip_vel

c     Local variables

      integer :: jm,jp,jpp,jgm,jgp,jgpp,ig,jg,kg

      real(8) :: jac,qjp,vyp

      logical :: bc_safe,zv,at_bdry

c     Begin program

      if (PRESENT(zip_vel)) then
        zv = zip_vel
      else
        zv = .true.
      endif

      bc_safe = .not.(method > 2 .and. method <= 6)

c     Update grid indices

      call getMGmap(i,j,k,igr,igr,igr,ig,jg,kg)

      jp  = j+1
      jm  = max(j-1,0)
      jpp = min(jp+1,ny+1)

      jgp  = jg + (jp-j)
      jgm  = jg + (jm-j)
      jgpp = jg + (jpp-j)

c     Grid info

      jac = gmetric%grid(igr)%jac(i,j,k)

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = isBdry(j,igr,4) .or. isBdry(j+1,igr,3)

      flxjp = 0.5*(v(i,jp,k)*arr(i,j,k) + v(i,j,k)*arr(i,jp,k))

c     Other advective methods

      if (method /= 2) then

        !Find velocities at faces from ZIP fluxes
        if (zv) then
          vyp = 2.*flxjp/(arr(i,jp,k)+arr(i,j,k)+1d-12)
        else
          vyp = 0.5*(v(i,jp,k)+v(i,j,k))
        endif

        !Find advected quantity at faces
        if ((j == ny .or. j == 0).and.(.not.bc_safe)) then
          qjp=in2face_scl(grid_params%yy(jgpp),arr(i,jpp,k)
     .                   ,grid_params%yy(jgp ),arr(i,jp ,k)
     .                   ,grid_params%yy(jg  ),arr(i,j  ,k)
     .                   ,grid_params%yy(jgm ),arr(i,jm ,k)
     .                   ,vyp,1)
        else
          qjp=in2face_scl(grid_params%yy(jgpp),arr(i,jpp,k)
     .                   ,grid_params%yy(jgp ),arr(i,jp ,k)
     .                   ,grid_params%yy(jg  ),arr(i,j  ,k)
     .                   ,grid_params%yy(jgm ),arr(i,jm ,k)
     .                   ,vyp,method)
        endif

        !Setup fluxes (ZIP at boundaries by default for conservation)
        if (.not.(at_bdry).or.bc_safe) flxjp=vyp*qjp

      endif

c     End program

      end function scl_advc_flx_y

c     scl_advc_flx_z
c     ###############################################################
      function scl_advc_flx_z(i,j,k,nx,ny,nz,igr,v,arr
     .                       ,method,zip_vel) result(flxkp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_3 of scalar arr at face
c     (i,j,k+1/2) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see inter)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxkp

      logical,optional :: zip_vel

c     Local variables

      integer :: km,kp,kpp,kgm,kgp,kgpp,ig,jg,kg

      real(8) :: jac,qkp,vzp

      logical :: bc_safe,zv,at_bdry

c     Begin program

      if (PRESENT(zip_vel)) then
        zv = zip_vel
      else
        zv = .true.
      endif

      bc_safe = .not.(method > 2 .and. method <= 6)

c     Update grid indices

      call getMGmap(i,j,k,igr,igr,igr,ig,jg,kg)

      kp  = k+1
      km  = max(k-1,0)
      kpp = min(kp+1,nz+1)

      kgp  = kg + (kp-k)
      kgm  = kg + (km-k)
      kgpp = kg + (kpp-k)

c     Grid info

      jac = gmetric%grid(igr)%jac(i,j,k)

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = isBdry(k,igr,6) .or. isBdry(k+1,igr,5)

      flxkp = 0.5*(v(i,j,kp)*arr(i,j,k) + v(i,j,k)*arr(i,j,kp))

c     Other advective methods

      if (method /= 2) then

        !Find velocities at faces from ZIP fluxes
        if (zv) then
          vzp = 2.*flxkp/(arr(i,j,kp)+arr(i,j,k)+1d-12)
        else
          vzp = 0.5*(v(i,j,kp)+v(i,j,k))
        endif

        !Find advected quantity at faces
        if ((k == nz .or. k == 0).and.(.not.bc_safe)) then
          qkp=in2face_scl(grid_params%zz(kgpp),arr(i,j,kpp)
     .                   ,grid_params%zz(kgp ),arr(i,j,kp )
     .                   ,grid_params%zz(kg  ),arr(i,j,k  )
     .                   ,grid_params%zz(kgm ),arr(i,j,km )
     .                   ,vzp,1)
        else
          qkp=in2face_scl(grid_params%zz(kgpp),arr(i,j,kpp)
     .                   ,grid_params%zz(kgp ),arr(i,j,kp )
     .                   ,grid_params%zz(kg  ),arr(i,j,k  )
     .                   ,grid_params%zz(kgm ),arr(i,j,km )
     .                   ,vzp,method)
        endif

        !Setup fluxes (ZIP at boundaries by default for conservation)
        if (.not.(at_bdry).or.bc_safe) flxkp=vzp*qkp

      endif

c     End program

      end function scl_advc_flx_z

c     flx_advec
c     ###############################################################
      function flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                  ,method,sp,vol,zip_vel) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see inter)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1,3),advec

      logical,optional :: sp,vol,zip_vel

c     Local variables

      integer :: ig,jg,kg

      real(8) :: dxx,dyy,dzz,jac
      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical :: vw

c     Begin program

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .false.
      endif

c     Grid info

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

      jac = gmetric%grid(igx)%jac(i,j,k)

c     Find fluxes

      flxip = scl_advc_flx_x(i  ,j,k,nx,ny,nz,igx,v(:,:,:,1),arr
     .                      ,method,sp=sp,zip_vel=zip_vel)
      flxim = scl_advc_flx_x(i-1,j,k,nx,ny,nz,igx,v(:,:,:,1),arr
     .                      ,method,sp=sp,zip_vel=zip_vel)

      flxjp = scl_advc_flx_y(i,j  ,k,nx,ny,nz,igx,v(:,:,:,2),arr
     .                      ,method,zip_vel=zip_vel)
      flxjm = scl_advc_flx_y(i,j-1,k,nx,ny,nz,igx,v(:,:,:,2),arr
     .                      ,method,zip_vel=zip_vel)

      flxkp = scl_advc_flx_z(i,j,k  ,nx,ny,nz,igx,v(:,:,:,3),arr
     .                      ,method,zip_vel=zip_vel)
      flxkm = scl_advc_flx_z(i,j,k-1,nx,ny,nz,igx,v(:,:,:,3),arr
     .                      ,method,zip_vel=zip_vel)

c     Find advective term

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac

      if (vw) advec = advec*gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function flx_advec

ccc     c_advec
ccc     ###############################################################
cc      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,phi,sp,upwind,vol)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of vector field (v phi) at cell centers in
ccc     general non-orthogonal geometry.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v  (0:nx+1,0:ny+1,0:nz+1,3),c_advec
cc
cc      logical,optional :: sp,vol,upwind
cc
ccc     Local variables
cc
cc      logical :: upw
cc
ccc     Begin program
cc
cc      if (PRESENT(upwind)) then
cc        upw = upwind
cc      else
cc        upw = .false.
cc      endif
cc
cc      if (upw) then
cc        c_advec = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,phi
cc     .                  ,1,sp=sp,vol=vol,zip_vel=.false.)
cc      else
cc        c_advec = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,phi
cc     .                  ,2,sp=sp,vol=vol,zip_vel=.true.)
cc      endif
cc
ccc     End 
cc
cc      end function c_advec
cc
ccc     c_advec
ccc     ###############################################################
cc      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,phi
cc     .                ,sp,upwind,vol)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of vector field (v.phi) at cell centers in
ccc     general non-orthogonal geometry.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),c_advec
cc
cc      logical,optional :: sp,vol,upwind
cc
ccc     Local variables
cc
cc      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
cc      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc      logical    :: vw,spoint,upw
cc
ccc     Begin program
cc
cc      if (PRESENT(vol)) then
cc        vw = vol
cc      else
cc        vw = .false.
cc      endif
cc
cc      if (PRESENT(upwind)) then
cc        upw = upwind
cc      else
cc        upw = .false.
cc      endif
cc
cc      if (PRESENT(sp)) then
cc        spoint = sp
cc      else
cc        spoint = .false.
cc      endif
cc
cc      igrid = igx
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      dxx = grid_params%dxh(ig)
cc      dyy = grid_params%dyh(jg)
cc      dzz = grid_params%dzh(kg)
cc
cc      jac   = gmetric%grid(igrid)%jac(i,j,k)
cc
cc      if (spoint) then
cc        jacip = gmetric%grid(igrid)%jac(ip,j,k)
cc        jacim = gmetric%grid(igrid)%jac(im,j,k)
cc      else
cc        jacip = jac
cc        jacim = jac
cc      endif
cc
cc      if (upw) then
cc        !X flux
cc        if (     i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .or. bcond(2) == PER) then
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.25*jach
cc     .      *( (    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)
cc     .          +abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(i ,j,k)
cc     .        +(    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)          
cc     .          -abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(ip,j,k) )
cc        else
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc        endif
cc
cc        if (     i+grid_params%ilo(igx)-1 > 1
cc     .      .or. spoint
cc     .      .or. bcond(1) == PER) then
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.25*jach
cc     .       *( (    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
cc     .           +abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(im,j,k)
cc     .         +(    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
cc     .           -abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(i ,j,k) )
cc        else
cc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cc        endif
cc
cc        if (     j+grid_params%jlo(igy)-1 < grid_params%nygl(igy)
cc     .      .or. bcond(4) == PER) then
cc          flxjp = 0.25*(
cc     .         (    (v2(i,j,k)+v2(i,jp,k))
cc     .          +abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,j ,k)
cc     .        +(    (v2(i,j,k)+v2(i,jp,k))          
cc     .          -abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,jp,k))
cc        else
cc          flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
cc        endif
cc
cc        if (     j+grid_params%jlo(igy)-1 > 1
cc     .      .or. bcond(3) == PER) then
cc          flxjm = 0.25*(
cc     .         (    (v2(i,j,k)+v2(i,jm,k))
cc     .          +abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,jm,k)
cc     .        +(    (v2(i,j,k)+v2(i,jm,k))          
cc     .          -abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,j ,k))
cc        else
cc          flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
cc        endif
cc
cc        if (     k+grid_params%klo(igz)-1 < grid_params%nzgl(igz)
cc     .      .or. bcond(6) == PER) then
cc          flxkp = 0.25*(
cc     .         (    (v3(i,j,k)+v3(i,j,kp))
cc     .          +abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,k )
cc     .        +(    (v3(i,j,k)+v3(i,j,kp))             
cc     .          -abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,kp))
cc        else
cc          flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
cc        endif
cc
cc        if (     k+grid_params%klo(igz)-1 > 1
cc     .      .or. bcond(5) == PER) then
cc          flxkm = 0.25*(
cc     .         (    (v3(i,j,k)+v3(i,j,km))
cc     .          +abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,km)
cc     .        +(    (v3(i,j,k)+v3(i,j,km))             
cc     .          -abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,k ))
cc        else
cc          flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
cc        endif
cc
cc      else
cc        !X flux
cc        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
cc     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
cc
cc          jach = 0.5*(jac+jacim)
cc          if (isSP(i,j,k,igx,igy,igz)) jach = 0d0
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        endif
cc
cccc        if (i+grid_params%ilo(igx)-1 == 1) then
cccc          jach = 0.5*(jac+jacip)
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
cccc     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
cccc          jach = 0.5*(jac+jacim)
cccc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cccc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cccc        elseif (i+grid_params%ilo(igx)-1 == 2) then
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cccc
cccc          jach = 0.5*(jac+jacim)
cccc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cccc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cccc        else
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cccc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cccc        endif
cc
cc        !Y flux
cc        flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
cc        flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
cc
cc        !Z flux
cc        flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
cc        flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
cc
cc      endif
cc
cc      c_advec =( (flxip - flxim)/dxx
cc     .         + (flxjp - flxjm)/dyy
cc     .         + (flxkp - flxkm)/dzz )/jac
cc      
cc      if (vw) c_advec = c_advec*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End 
cc
cc      end function c_advec

ccc     flx_advec
ccc     ###############################################################
cc      function flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
cc     .                  ,method,sp,vol,zip_vel,conserv) result(advec)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of vector field (v.arr) at cell centers in
ccc     general non-orthogonal geometry.
ccc
ccc     method: specifies the advective scheme (see inter)
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer :: i,j,k,nx,ny,nz,igx,igy,igz,method
cc
cc      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
cc     .          ,v  (0:nx+1,0:ny+1,0:nz+1,3),advec
cc
cc      logical,optional :: sp,vol,zip_vel,conserv
cc
ccc     Local variables
cc
cc      integer :: igrid,half_elem,su
cc
cc      integer :: im,ip,ipp,imm,igm,igp,igpp,igmm
cc     .          ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
cc     .          ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
cc     .          ,ig,jg,kg,igl,jgl,kgl
cc
cc      real(8) :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
cc      real(8) :: vxp,vxm,vyp,vym,vzp,vzm
cc     .          ,qip,qim,qjp,qjm,qkp,qkm
cc      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc      logical :: vw,spoint,bc_safe,zv,cv
cc
ccc     Begin program
cc
cc      if (PRESENT(conserv)) then
cc        cv = conserv
cc      else
cc        cv = .true.
cc      endif
cc
cc      if (.not.cv) then
cc        advec = flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
cc     .                      ,method,sp=sp,vol=vol)
cc        return
cc      endif
cc
ccc     Set defaults
cc
cc      if (PRESENT(vol)) then
cc        vw = vol
cc      else
cc        vw = .false.
cc      endif
cc
cc      if (PRESENT(sp)) then
cc        spoint = sp
cc      else
cc        spoint = .false.
cc      endif
cc
cc      if (PRESENT(zip_vel)) then
cc        zv = zip_vel
cc      else
cc        zv = .true.
cc      endif
cc
cccc      bc_safe = (method <=2)
cc      bc_safe = .not.(method > 2 .and. method <= 6)
ccc     General setup
cc
cc      igrid = igx
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      ipp = ip+1
cc      imm = im-1
cc      ipp = min(ipp,nx+1)
cc      imm = max(imm,0)
cc
cc      igp  = ig+1
cc      igm  = ig-1
cc      igpp = ig + (ipp - i)
cc      igmm = ig + (imm - i)
cc
cc      jpp = jp+1
cc      jmm = jm-1
cc      jpp = min(jpp,ny+1)
cc      jmm = max(jmm,0)
cc
cc      jgp  = jg+1
cc      jgm  = jg-1
cc      jgpp = jg + (jpp-j)
cc      jgmm = jg + (jmm-j)
cc
cc      kpp = kp+1
cc      kmm = km-1
cc      kpp = min(kpp,nz+1)
cc      kmm = max(kmm,0)
cc
cc      kgp  = kg+1
cc      kgm  = kg-1
cc      kgpp = kg + (kpp-k)
cc      kgmm = kg + (kmm-k)
cc
cc      igl = i+grid_params%ilo(igx)-1
cc      jgl = j+grid_params%jlo(igy)-1
cc      kgl = k+grid_params%klo(igz)-1
cc
ccc     Grid info
cc
cc      dxx = grid_params%dxh(ig)
cc      dyy = grid_params%dyh(jg)
cc      dzz = grid_params%dzh(kg)
cc
cc      jac = gmetric%grid(igrid)%jac(i,j,k)
cc
cc      if (spoint) then
cc        jacip = gmetric%grid(igrid)%jac(ip,j,k)
cc        jacim = gmetric%grid(igrid)%jac(im,j,k)
cc      else
cc        jacip = jac
cc        jacim = jac
cc      endif
cc
ccc     Find ZIP fluxes (used even if method /= ZIP)
cc
cc      flxip = 0d0
cc      flxim = 0d0
cc      flxjp = 0d0
cc      flxjm = 0d0
cc      flxkp = 0d0
cc      flxkm = 0d0
cc
cc      !X flux
cc      if (.not.isBdry(i,igx,2)) then  !Not at r=1 boundary
cc        jach  = 0.5*(jac+jacip)
cc        flxip = 0.5*(v(ip,j,k,1)*arr(i ,j,k)/jacip
cc     .             + v(i ,j,k,1)*arr(ip,j,k)/jac  )*jach
cc
cc        if (isSP2(i,igx,ibc=1).and.spoint) then
cc          flxim = 0d0
cc        else
cc          jach  = 0.5*(jac+jacim)
cc          flxim = 0.5*(v(im,j,k,1)*arr(i ,j,k)/jacim
cc     .               + v(i ,j,k,1)*arr(im,j,k)/jac  )*jach
cc        endif
cc      else
cc        if (v(ip,j,k,1) == -v(i,j,k,1)) then
cc          flxip = 0d0    !Shut off flux for vn=0 BC
cc        else
cc          flxip = 0.5*(v(ip,j,k,1)*arr(i,j,k)+v(i,j,k,1)*arr(ip,j,k))
cc        endif
cc
cc        jach  = 0.5*(jac+jacim)
cc        flxim = 0.5*(v(im,j,k,1)*arr(i ,j,k)/jacim
cc     .             + v(i ,j,k,1)*arr(im,j,k)/jac  )*jach
cc      endif
cc
cc      !Y flux
cc      flxjp = 0.5*(v(i,jp,k,2)*arr(i,j,k) + v(i,j,k,2)*arr(i,jp,k))
cc      flxjm = 0.5*(v(i,jm,k,2)*arr(i,j,k) + v(i,j,k,2)*arr(i,jm,k))
cc
cc      !Z flux
cc      flxkp = 0.5*(v(i,j,kp,3)*arr(i,j,k) + v(i,j,k,3)*arr(i,j,kp))
cc      flxkm = 0.5*(v(i,j,km,3)*arr(i,j,k) + v(i,j,k,3)*arr(i,j,km))
cc
ccc     Other advective methods
cc
cc      if (method /= 2) then
cc
cc        !Find velocities at faces from ZIP fluxes
cc        if (zv) then
cc          vxp = 2.*flxip/(arr(ip,j,k)+arr(i,j,k)+1d-12)
cc          vxm = 2.*flxim/(arr(im,j,k)+arr(i,j,k)+1d-12)
cc                                                    
cc          vyp = 2.*flxjp/(arr(i,jp,k)+arr(i,j,k)+1d-12)
cc          vym = 2.*flxjm/(arr(i,jm,k)+arr(i,j,k)+1d-12)
cc                                                    
cc          vzp = 2.*flxkp/(arr(i,j,kp)+arr(i,j,k)+1d-12)
cc          vzm = 2.*flxkm/(arr(i,j,km)+arr(i,j,k)+1d-12)
cc        else
cc          vxp = 0.5*(v(ip,j,k,1)+v(i,j,k,1))
cc          vxm = 0.5*(v(im,j,k,1)+v(i,j,k,1))
cc                                                    
cc          vyp = 0.5*(v(i,jp,k,2)+v(i,j,k,2))
cc          vym = 0.5*(v(i,jm,k,2)+v(i,j,k,2))
cc                                                    
cc          vzp = 0.5*(v(i,j,kp,3)+v(i,j,k,3))
cc          vzm = 0.5*(v(i,j,km,3)+v(i,j,k,3))
cc        endif
cc
cc        !Find advected quantity at faces
cc        if ((i == nx).and.(.not.bc_safe)) then
cc          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
cc     .             ,grid_params%xx(igp ),arr(ip ,j,k)
cc     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
cc     .             ,grid_params%xx(igm ),arr(im ,j,k)
cc     .             ,vxp,1)
cc        else
cc          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
cc     .             ,grid_params%xx(igp ),arr(ip ,j,k)
cc     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
cc     .             ,grid_params%xx(igm ),arr(im ,j,k)
cc     .             ,vxp,method)
cc        endif
cc
cc        if ((i == 1).and.(.not.bc_safe)) then
cc          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
cc     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
cc     .             ,grid_params%xx(igm ),arr(im ,j,k)
cc     .             ,grid_params%xx(igmm),arr(imm,j,k)
cc     .             ,vxm,1)
cc        else
cc          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
cc     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
cc     .             ,grid_params%xx(igm ),arr(im ,j,k)
cc     .             ,grid_params%xx(igmm),arr(imm,j,k)
cc     .             ,vxm,method)
cc        endif
cc
cc        if ((j == ny).and.(.not.bc_safe)) then
cc          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
cc     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
cc     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
cc     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
cc     .             ,vyp,1)
cc        else
cc          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
cc     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
cc     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
cc     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
cc     .             ,vyp,method)
cc        endif
cc
cc        if ((j == 1).and.(.not.bc_safe)) then
cc          qjm=inter(grid_params%yy(jgp ),arr(i,jp ,k)
cc     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
cc     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
cc     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
cc     .             ,vym,1)
cc        else
cc          qjm=inter(grid_params%yy(jgp ),arr(i,jp,k )
cc     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
cc     .             ,grid_params%yy(jgm ),arr(i,jm,k )
cc     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
cc     .             ,vym,method)
cc        endif
cc
cc        if ((k == nz).and.(.not.bc_safe)) then
cc          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
cc     .             ,grid_params%zz(kgp ),arr(i,j,kp )
cc     .             ,grid_params%zz(kg  ),arr(i,j,k  )
cc     .             ,grid_params%zz(kgm ),arr(i,j,km )
cc     .             ,vzp,1)
cc        else
cc          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
cc     .             ,grid_params%zz(kgp ),arr(i,j,kp )
cc     .             ,grid_params%zz(kg  ),arr(i,j,k  )
cc     .             ,grid_params%zz(kgm ),arr(i,j,km )
cc     .             ,vzp,method)
cc        endif
cc
cc        if ((k == 1).and.(.not.bc_safe)) then
cc          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
cc     .             ,grid_params%zz(kg  ),arr(i,j,k  )
cc     .             ,grid_params%zz(kgm ),arr(i,j,km )
cc     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
cc     .             ,vzm,1)
cc        else
cc          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
cc     .             ,grid_params%zz(kg  ),arr(i,j,k  )
cc     .             ,grid_params%zz(kgm ),arr(i,j,km )
cc     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
cc     .             ,vzm,method)
cc        endif
cc
cc        !Setup fluxes (ZIP at boundaries by default for conservation)
cc        if (.not.isBdry(i,igx,2).or.bc_safe) flxip=vxp*qip
cccc        bc_safe = (method<=2).and.(.not.(isSP2(i,igx,ibc=1).and.spoint))
cc        if (.not.isBdry(i,igx,1).or.bc_safe) flxim=vxm*qim
cc
cc        if (.not.isBdry(j,igy,4).or.bc_safe) flxjp=vyp*qjp
cc        if (.not.isBdry(j,igy,3).or.bc_safe) flxjm=vym*qjm
cc
cc        if (.not.isBdry(k,igz,6).or.bc_safe) flxkp=vzp*qkp
cc        if (.not.isBdry(k,igz,5).or.bc_safe) flxkm=vzm*qkm
cc
cc      endif
cc
ccc     Find advective term
cc
cc      advec =( (flxip - flxim)/dxx
cc     .       + (flxjp - flxjm)/dyy
cc     .       + (flxkp - flxkm)/dzz )/jac
cc
cc      if (vw) advec = advec*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End 
cc
cc      end function flx_advec

ccc     c_advec
ccc     ###############################################################
cc      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,phi,sp,upwind,vol)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of vector field (v phi) at cell centers in
ccc     general non-orthogonal geometry.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v  (0:nx+1,0:ny+1,0:nz+1,3),c_advec
cc
cc      logical,optional :: sp,vol,upwind
cc
ccc     Local variables
cc
cc      logical :: upw
cc
ccc     Begin program
cc
cc      if (PRESENT(upwind)) then
cc        upw = upwind
cc      else
cc        upw = .false.
cc      endif
cc
cc      if (upw) then
cc        c_advec = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,phi
cc     .                  ,1,sp=sp,vol=vol,zip_vel=.false.)
cc      else
cc        c_advec = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,phi
cc     .                  ,2,sp=sp,vol=vol,zip_vel=.true.)
cc      endif
cc
ccc     End 
cc
cc      end function c_advec
cc
ccc     c_advec
ccc     ###############################################################
cc      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,phi
cc     .                ,sp,upwind,vol)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of vector field (v.phi) at cell centers in
ccc     general non-orthogonal geometry.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),c_advec
cc
cc      logical,optional :: sp,vol,upwind
cc
ccc     Local variables
cc
cc      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
cc      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc      logical    :: vw,spoint,upw
cc
ccc     Begin program
cc
cc      if (PRESENT(vol)) then
cc        vw = vol
cc      else
cc        vw = .false.
cc      endif
cc
cc      if (PRESENT(upwind)) then
cc        upw = upwind
cc      else
cc        upw = .false.
cc      endif
cc
cc      if (PRESENT(sp)) then
cc        spoint = sp
cc      else
cc        spoint = .false.
cc      endif
cc
cc      igrid = igx
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      dxx = grid_params%dxh(ig)
cc      dyy = grid_params%dyh(jg)
cc      dzz = grid_params%dzh(kg)
cc
cc      jac   = gmetric%grid(igrid)%jac(i,j,k)
cc
cc      if (spoint) then
cc        jacip = gmetric%grid(igrid)%jac(ip,j,k)
cc        jacim = gmetric%grid(igrid)%jac(im,j,k)
cc      else
cc        jacip = jac
cc        jacim = jac
cc      endif
cc
cc      if (upw) then
cc        !X flux
cc        if (     i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .or. bcond(2) == PER) then
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.25*jach
cc     .      *( (    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)
cc     .          +abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(i ,j,k)
cc     .        +(    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)          
cc     .          -abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(ip,j,k) )
cc        else
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc        endif
cc
cc        if (     i+grid_params%ilo(igx)-1 > 1
cc     .      .or. spoint
cc     .      .or. bcond(1) == PER) then
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.25*jach
cc     .       *( (    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
cc     .           +abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(im,j,k)
cc     .         +(    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
cc     .           -abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(i ,j,k) )
cc        else
cc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cc        endif
cc
cc        if (     j+grid_params%jlo(igy)-1 < grid_params%nygl(igy)
cc     .      .or. bcond(4) == PER) then
cc          flxjp = 0.25*(
cc     .         (    (v2(i,j,k)+v2(i,jp,k))
cc     .          +abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,j ,k)
cc     .        +(    (v2(i,j,k)+v2(i,jp,k))          
cc     .          -abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,jp,k))
cc        else
cc          flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
cc        endif
cc
cc        if (     j+grid_params%jlo(igy)-1 > 1
cc     .      .or. bcond(3) == PER) then
cc          flxjm = 0.25*(
cc     .         (    (v2(i,j,k)+v2(i,jm,k))
cc     .          +abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,jm,k)
cc     .        +(    (v2(i,j,k)+v2(i,jm,k))          
cc     .          -abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,j ,k))
cc        else
cc          flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
cc        endif
cc
cc        if (     k+grid_params%klo(igz)-1 < grid_params%nzgl(igz)
cc     .      .or. bcond(6) == PER) then
cc          flxkp = 0.25*(
cc     .         (    (v3(i,j,k)+v3(i,j,kp))
cc     .          +abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,k )
cc     .        +(    (v3(i,j,k)+v3(i,j,kp))             
cc     .          -abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,kp))
cc        else
cc          flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
cc        endif
cc
cc        if (     k+grid_params%klo(igz)-1 > 1
cc     .      .or. bcond(5) == PER) then
cc          flxkm = 0.25*(
cc     .         (    (v3(i,j,k)+v3(i,j,km))
cc     .          +abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,km)
cc     .        +(    (v3(i,j,k)+v3(i,j,km))             
cc     .          -abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,k ))
cc        else
cc          flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
cc        endif
cc
cc      else
cc        !X flux
cc        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
cc     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
cc
cc          jach = 0.5*(jac+jacim)
cc          if (isSP(i,j,k,igx,igy,igz)) jach = 0d0
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        endif
cc
cccc        if (i+grid_params%ilo(igx)-1 == 1) then
cccc          jach = 0.5*(jac+jacip)
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
cccc     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
cccc          jach = 0.5*(jac+jacim)
cccc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cccc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cccc        elseif (i+grid_params%ilo(igx)-1 == 2) then
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cccc
cccc          jach = 0.5*(jac+jacim)
cccc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cccc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cccc        else
cccc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cccc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cccc        endif
cc
cc        !Y flux
cc        flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
cc        flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
cc
cc        !Z flux
cc        flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
cc        flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
cc
cc      endif
cc
cc      c_advec =( (flxip - flxim)/dxx
cc     .         + (flxjp - flxjm)/dyy
cc     .         + (flxkp - flxkm)/dzz )/jac
cc      
cc      if (vw) c_advec = c_advec*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End 
cc
cc      end function c_advec

c     flx_nc_advec
c     ###############################################################
      function flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                     ,method,sp,vol,reverse) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see inter)
c     ---------------------------------------------------------------

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz,method

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1)
     .             ,v  (0:nx+1,0:ny+1,0:nz+1,3),advec

      logical,optional :: sp,vol,reverse

c     Local variables

      integer    :: igrid,half_elem,su

      integer    :: im,ip,ipp,imm,igm,igp,igpp,igmm
     .             ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
     .             ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
     .             ,ig,jg,kg,igl,jgl,kgl

      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jac,sgn
      real(8)    :: vxp,vxm,vyp,vym,vzp,vzm
     .             ,qip,qim,qjp,qjm,qkp,qkm
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw,spoint,bc_safe,rev_sgn

c     Begin program

c     Set defaults

      bc_safe = .not.(method > 2 .and. method <= 6)

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .false.
      endif

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = .false.
      endif

      if (PRESENT(reverse)) then
        rev_sgn = reverse
      else
        rev_sgn = .false.
      endif

      sgn = 1d0
      if (rev_sgn) sgn=-1d0  !Takes (-v) instead of v

c     General setup

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      ipp = ip+1
      imm = im-1
      ipp = min(ipp,nx+1)
      imm = max(imm,0)

      igp  = ig+1
      igm  = ig-1
      igpp = ig + (ipp - i)
      igmm = ig + (imm - i)

      jpp = jp+1
      jmm = jm-1
      jpp = min(jpp,ny+1)
      jmm = max(jmm,0)

      jgp  = jg+1
      jgm  = jg-1
      jgpp = jg + (jpp-j)
      jgmm = jg + (jmm-j)

      kpp = kp+1
      kmm = km-1
      kpp = min(kpp,nz+1)
      kmm = max(kmm,0)

      kgp  = kg+1
      kgm  = kg-1
      kgpp = kg + (kpp-k)
      kgmm = kg + (kmm-k)

      igl = i+grid_params%ilo(igx)-1
      jgl = j+grid_params%jlo(igy)-1
      kgl = k+grid_params%klo(igz)-1

c     Grid info

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

      jac = gmetric%grid(igrid)%jac(i,j,k)

c     Find DEFAULT non-conservative fluxes (central differences)

      !X flux
      flxip = 0.5*(arr(i,j,k) + arr(ip,j,k))*v(i,j,k,1)
      flxim = 0.5*(arr(i,j,k) + arr(im,j,k))*v(i,j,k,1)

      !Y flux
      flxjp = 0.5*(arr(i,j,k) + arr(i,jp,k))*v(i,j,k,2)
      flxjm = 0.5*(arr(i,j,k) + arr(i,jm,k))*v(i,j,k,2)

      !Z flux
      flxkp = 0.5*(arr(i,j,k) + arr(i,j,kp))*v(i,j,k,3)
      flxkm = 0.5*(arr(i,j,k) + arr(i,j,km))*v(i,j,k,3)

c     Find specialized non-conservative fluxes (reverts to CD at logical domain boundaries)

      if (method /= 2) then

        vxp = sgn*v(i,j,k,1)
        vxm = sgn*v(i,j,k,1)
                                                    
        vyp = sgn*v(i,j,k,2)
        vym = sgn*v(i,j,k,2)
                                                    
        vzp = sgn*v(i,j,k,3)
        vzm = sgn*v(i,j,k,3)

        !Find advected quantity at faces 
        if ((i == nx).and.(.not.bc_safe)) then
          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
     .             ,grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,vxp,2)
        else
          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
     .             ,grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,vxp,method)
        endif

        if ((i == 1).and.(.not.bc_safe)) then
          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,grid_params%xx(igmm),arr(imm,j,k)
     .             ,vxm,2)
        else
          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,grid_params%xx(igmm),arr(imm,j,k)
     .             ,vxm,method)
        endif

        if ((j == ny).and.(.not.bc_safe)) then
          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,2)
        else
          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,method)
        endif

        if ((j == 1).and.(.not.bc_safe)) then
          qjm=inter(grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
     .             ,vym,2)
        else
          qjm=inter(grid_params%yy(jgp ),arr(i,jp,k )
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm,k )
     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
     .             ,vym,method)
        endif

        if ((k == nz).and.(.not.bc_safe)) then
          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
     .             ,grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,vzp,2)
        else
          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
     .             ,grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,vzp,method)
        endif

        if ((k == 1).and.(.not.bc_safe)) then
          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,2)
        else
          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,method)
        endif

        !Setup fluxes
        if (.not.isBdry(i,igx,2).or.bc_safe) flxip=vxp*qip
        if (.not.isBdry(i,igx,1).or.bc_safe) flxim=vxm*qim
                                            
        if (.not.isBdry(j,igy,4).or.bc_safe) flxjp=vyp*qjp
        if (.not.isBdry(j,igy,3).or.bc_safe) flxjm=vym*qjm
                                            
        if (.not.isBdry(k,igz,6).or.bc_safe) flxkp=vzp*qkp
        if (.not.isBdry(k,igz,5).or.bc_safe) flxkm=vzm*qkm

      endif

c     Find advective term

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac
      
      if (vw) advec = advec*gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function flx_nc_advec

c     tensor_nc_advc
c     #####################################################################
      function tensor_nc_advc(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,v1,meth)
     .         result(cnv)

c     ---------------------------------------------------------------------
c     Calculates v0.nabla(v1), with both v0, v1 vector fields.
c
c     meth: specifies the advective scheme (only upwind=1 and central=2)
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz,meth
      real(8) :: cnv(3)
      real(8),target :: v1(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .                 ,v0(0:nxx+1,0:nyy+1,0:nzz+1,3)

c     Local variables

      integer    :: ieq,hex,hey,hez
      real(8)    :: nabla_v(3,3),jac

c     Begin program

      jac = gmetric%grid(igx)%jac(i,j,k)

      select case(meth)
      case(1) !Upwinding
        hex = int(sign(1d0,-v0(i,j,k,1)))
        hey = int(sign(1d0,-v0(i,j,k,2)))
        hez = int(sign(1d0,-v0(i,j,k,3)))

        nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                         ,v1(:,:,:,1)
     .                         ,v1(:,:,:,2)
     .                         ,v1(:,:,:,3),hex,hey,hez)
      case(2) !Central differences
        nabla_v = fnabla_v(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                    ,v1(:,:,:,1)
     .                    ,v1(:,:,:,2)
     .                    ,v1(:,:,:,3),0)
      case default
        call pstop('tensor_nc_advc','Advection option not implemented')
      end select

      cnv = matmul(v0(i,j,k,:),nabla_v)/jac

      end function tensor_nc_advc

c     tensor_advc
c     #####################################################################
      function tensor_advc(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,v1,method)
     .         result(cnv)

c     ---------------------------------------------------------------------
c     Calculates divergence of tensor, div(v0 v1), at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see inter)
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz,method
      real(8) :: cnv(3)
      real(8),target :: v1(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .                 ,v0(0:nxx+1,0:nyy+1,0:nzz+1,3)

c     Local variables

c     Begin program

      adv_meth = method

      vec1 => v0
      vec2 => v1
      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt__eom()
     .                ,tnsr_advc_x,tnsr_advc_y,tnsr_advc_z)
      nullify(vec1,vec2)

      end function tensor_advc

c     tnsr_advc_x
c     #############################################################
      subroutine tnsr_advc_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for tensor adv. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate tensor adv in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm,icomp,adv__meth
        real(8)    :: vxp,qp(3),t1(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jac = 0.5*(gmetric%grid(igx)%jac(ip,j,k)
     .            +gmetric%grid(igx)%jac(i ,j,k))

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

c     ZIP advective fluxes

        t1(1) = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,1)
     .                +vec2(i ,j,k,1)*vec1(ip,j,k,1)
     .                +vec1(i ,j,k,1)*vec2(ip,j,k,1)
     .                +vec1(ip,j,k,1)*vec2(i ,j,k,1))
     .              *ijacp*ijac0*jac**2

        t1(2) = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,2)*ijacp
     .                +vec2(i ,j,k,1)*vec1(ip,j,k,2)*ijac0
     .                +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
     .                +vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)
     .              *jac

        t1(3) = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,3)
     .                +vec2(i ,j,k,1)*vec1(ip,j,k,3)
     .                +vec1(i ,j,k,1)*vec2(ip,j,k,3)
     .                +vec1(ip,j,k,1)*vec2(i ,j,k,3))
     .              *ijacp*ijac0*jac**2

c     Flux limiting

        if (adv_meth /= 2 .and. flag /= 0) then

          !Revert extended support methods to upwinding at processor boundaries
          if ((i == nx .or. i == 0)
     .       .and.(adv_meth > 2 .and. adv_meth <= 6)) then
            adv__meth = 1
          else
            adv__meth = adv_meth
          endif

          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

          ipp  = min(ip+1,nx+1)
          im   = max(i -1,0   )

          igp  = ig + (ip  - i)
          igm  = ig + (im  - i)
          igpp = ig + (ipp - i)

          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac

          !Limit flux based on velocity
          t1 = vxp*inter(grid_params%xx(igpp),vec2(ipp,j,k,:)
     .                  ,grid_params%xx(igp ),vec2(ip ,j,k,:)
     .                  ,grid_params%xx(ig  ),vec2(i  ,j,k,:)
     .                  ,grid_params%xx(igm ),vec2(im ,j,k,:)
     .                  ,vxp,adv__meth)
        endif

c     Postprocessing

        t11 = t1(1)
        t12 = t1(2)
        t13 = t1(3)

        if (flag /= 0) then
          ijac  = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine tnsr_advc_x

c     tnsr_advc_y
c     #############################################################
      subroutine tnsr_advc_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for div(v1 v2). In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate tensor adv in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8) :: t21,t22,t23
        logical :: alt_eom

c     Local variables

        integer :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm,adv__meth
        real(8) :: vyp,qp(3),jac,ijac,t2(3)

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jac = 0.5*(gmetric%grid(igx)%jac(i,jp,k)
     .            +gmetric%grid(igx)%jac(i,j ,k))

c     ZIP advective fluxes

        t2 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,:)
     .             +vec2(i,j ,k,2)*vec1(i,jp,k,:)
     .             +vec1(i,j ,k,2)*vec2(i,jp,k,:)
     .             +vec1(i,jp,k,2)*vec2(i,j ,k,:))

c     Flux limiting

        if (adv_meth /= 2 .and. flag /= 0) then

          !Revert extended support methods to upwinding at processor boundaries
          if ((j == ny .or. j == 0)
     .       .and.(adv_meth > 2 .and. adv_meth <= 6)) then
            adv__meth = 1
          else
            adv__meth = adv_meth
          endif

          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

          jpp = min(jp+1,ny+1)
          jm  = max(j -1,0   )

          jgp  = jg + (jp  - j)
          jgm  = jg + (jm  - j)
          jgpp = jg + (jpp - j)

          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))

          !Limit flux based on velocity
          t2 = vyp*inter(grid_params%yy(jgpp),vec2(i,jpp,k,:)
     .                  ,grid_params%yy(jgp ),vec2(i,jp ,k,:)
     .                  ,grid_params%yy(jg  ),vec2(i,j  ,k,:)
     .                  ,grid_params%yy(jgm ),vec2(i,jm ,k,:)
     .                  ,vyp,adv__meth)

        endif

c     Postprocessing

        t21 = t2(1)
        t22 = t2(2)
        t23 = t2(3)

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine tnsr_advc_y

c     tnsr_advc_z
c     #############################################################
      subroutine tnsr_advc_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for tensor adv. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate tensor adv in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8) :: t31,t32,t33
        logical :: alt_eom

c     Local variables

        integer :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm,adv__meth
        real(8) :: vzp,qp(3),jac,ijac,t3(3)

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac= 0.5*(gmetric%grid(igx)%jac(i,j,kp)
     .           +gmetric%grid(igx)%jac(i,j,k ))

c     ZIP advective fluxes

        t3 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,:)
     .             +vec2(i,j,k ,3)*vec1(i,j,kp,:)
     .             +vec1(i,j,k ,3)*vec2(i,j,kp,:)
     .             +vec1(i,j,kp,3)*vec2(i,j,k ,:))

c     Flux limiting

        if (adv_meth /= 2 .and. flag /= 0) then

         !Revert extended support methods to upwinding at processor boundaries
          if ((k == nz .or. k == 0)
     .       .and.(adv_meth > 2 .and. adv_meth <= 6)) then
            adv__meth = 1
          else
            adv__meth = adv_meth
          endif

          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

          kpp = min(kp+1,nz+1)
          km  = max(k -1,0   )

          kgp  = kg + (kp  - k)
          kgm  = kg + (km  - k)
          kgpp = kg + (kpp - k)

          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))

            !Limit flux based on velocity
          t3=vzp*inter(grid_params%zz(kgpp),vec2(i,j,kpp,:)
     .                ,grid_params%zz(kgp ),vec2(i,j,kp ,:)
     .                ,grid_params%zz(kg  ),vec2(i,j,k  ,:)
     .                ,grid_params%zz(kgm ),vec2(i,j,km ,:)
     .                ,vzp,adv__meth)

        endif

c     Postprocessing

        t31 = t3(1)
        t32 = t3(2)
        t33 = t3(3)

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine tnsr_advc_z

c     in2face_scl
c     ###############################################################
      function in2face_scl(x4,q4,x3,q3,x2,q2,x1,q1,vel,advect)
     .         result(inter)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation at control
c     volume face based of face velocity "vel". Options:
c       * advect = 1 => First-order upwind
c       * advect = 2 => Centered
c       * advect = 3 => High-order upwind (QUICK)
c       * advect = 4 => Monotone high-order upwind (SMART)
c       * advect = 5 => SMART with smooth transition
c       * advect = 6 => Centered, high-order
c       * advect = 7 => Gamma (modification of SMART; not operational)
c
c     Convention:
c
c       ---x-----x--o--x-----x---->
c          q1    q2 ^  q3    q4
c                   ^
c           Location of face
c     ---------------------------------------------------------------

      use oned_int

      implicit none                !For safe fortran

c     Call variables

      real(8)    :: x4,x3,x2,x1,q4,q3,q2,q1,vel,inter
      integer    :: advect

c     Local variables

      real(8)    :: qt1,qt2,qt3,qt4,slp1,slp2,a,b,c,xp1,xp2,curv

c     Begin program

      select case(advect)
      case(1)

c     Upwind

        if (vel.gt.0d0) then
          inter = q2
        else
          inter = q3
        endif

      case(2)

c     Centered

        inter = 0.5*(q3+q2)

      case(3)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
          curv = ((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/0.5/(x3-x1)
        else
          curv = ((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/0.5/(x4-x2)
        endif

        inter = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2

      case(4)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = 1.5
        slp2 = .5
        if (vel.gt.0d0) then
          curv = ((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/0.5/(x3-x1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q2 + (1.-slp1)*q1
          qt3 = slp2*q2 + (1.-slp2)*q3
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
          curv = ((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/0.5/(x4-x2)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q3 + (1.-slp1)*q4
          qt3 = slp2*q3 + (1.-slp2)*q2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(5)

c     Smooth SMART

        xp1 = .3
        xp2 = .5
        if (vel.gt.0d0) then
          qt1  = (q2-q1)/(q3-q1)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1.)
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 3.*q3/8. + 3.*q2/4.- q1/8.
          else
            inter = q2
          endif
        else
          qt1  = (q3-q4)/(q2-q4)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            c = 1.
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 3.*q2/8. + 3.*q3/4.- q4/8.
          else
            inter = q3
          endif
        endif

      case(6)

c     Centered high-order

        inter = (-q1 + 9.*q2 + 9.*q3 - q4)/16.

cc      case(7)
cc
ccc     Gamma [Jasak, Weller, and Gosman, Int. J. Numer. Meth. Fluids, 32, p431 (1999)]
cc
cc        slp1 = 1.5
cc        slp2 = .5
cc        if (vel.gt.0d0) then
cc          curv = ((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/0.5/(x3-x1)
cc          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
cccc          qt1 = 3.*q3/8. + 3.*q2/4.- q1/8.
cc          qt2 = slp1*q2 + (1.-slp1)*q1
cc          qt3 = slp2*q2 + (1.-slp2)*q3
cc          qt4 = fmed(q2,qt2,qt3)
cc          inter = fmed(q2,qt4,qt1)
cc        else
cc          curv = ((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/0.5/(x4-x2)
cc          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
cccc          qt1 = 3.*q2/8. + 3.*q3/4.- q4/8.
cc          qt2 = slp1*q3 + (1.-slp1)*q4
cc          qt3 = slp2*q3 + (1.-slp2)*q2
cc          qt4 = fmed(q3,qt2,qt3)
cc          inter = fmed(q3,qt4,qt1)
cc        endif

      case default

        call pstop('in2face_scl','Advection scheme not available')

      end select

c     End

      end function in2face_scl

c     in2face_vec
c     ###############################################################
      function in2face_vec(x4,q4,x3,q3,x2,q2,x1,q1,vel,advect)
     .         result(inter)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation of vector
c     quantity at control volume face based of face velocity "vel".
c     Options:
c       * advect = 1 => First-order upwind
c       * advect = 2 => Centered
c       * advect = 3 => High-order upwind (QUICK)
c       * advect = 4 => Monotone high-order upwind (SMART)
c       * advect = 5 => SMART with smooth transition
c       * advect = 6 => Centered, high-order
c       * advect = 7 => Gamma (modification of SMART; not operational)
c
c     Convention:
c
c       ---x-----x--o--x-----x---->
c          q1    q2 ^  q3    q4
c                   ^
c           Location of face
c     ---------------------------------------------------------------

      use oned_int

      implicit none                !For safe fortran

c     Call variables

      real(8) :: x4,x3,x2,x1,q4(:),q3(:),q2(:),q1(:),vel
      real(8) :: inter(size(q1))
      integer :: advect

c     Local variables

      real(8) :: qt1(size(q1)),qt2(size(q1))
     .          ,qt3(size(q1)),qt4(size(q1))
     .          ,slp1,slp2,a(size(q1)),b(size(q1)),c(size(q1))
     .          ,xp1,xp2,curv(size(q1))

c     Begin program

      select case(advect)
      case(1)

c     Upwind

        if (vel.gt.0d0) then
          inter = q2
        else
          inter = q3
        endif

      case(2)

c     Centered

        inter = 0.5*(q3+q2)

      case(3)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
          curv = ((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/0.5/(x3-x1)
        else
          curv = ((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/0.5/(x4-x2)
        endif

        inter = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2

      case(4)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = 1.5
        slp2 = .5
        if (vel.gt.0d0) then
          curv = ((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/0.5/(x3-x1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q2 + (1.-slp1)*q1
          qt3 = slp2*q2 + (1.-slp2)*q3
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
          curv = ((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/0.5/(x4-x2)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q3 + (1.-slp1)*q4
          qt3 = slp2*q3 + (1.-slp2)*q2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(5)

c     Smooth SMART

        xp1 = .3
        xp2 = .5
        if (vel.gt.0d0) then
          qt1  = (q2-q1)/(q3-q1)
cc          if (qt1.lt.xp1.and.qt1.ge.0.) then
cc            a = .25/xp1**3*(xp1-3.)
cc            b = 1./xp1**2*(9./8.-xp1/2.)
cc            qt2 = qt1*(a*qt1**2 + b*qt1 + 1.)
cc            inter = q1 + qt2*(q3-q1)
cc          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
cc            a = (xp2/4.-.5)/(xp2-1.)**3
cc            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
cc            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
cc            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
cc            inter = q1 + qt2*(q3-q1)
cc          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
cc            inter = 3.*q3/8. + 3.*q2/4.- q1/8.
cc          else
cc            inter = q2
cc          endif
          where (qt1.lt.xp1.and.qt1.ge.0.)
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1.)
            inter = q1 + qt2*(q3-q1)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1.)
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q1 + qt2*(q3-q1)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
            inter = 3.*q3/8. + 3.*q2/4.- q1/8.
          elsewhere
            inter = q2
          end where
        else
          qt1  = (q3-q4)/(q2-q4)
cc          if (qt1.lt.xp1.and.qt1.ge.0.) then
cc            a = .25/xp1**3*(xp1-3.)
cc            b = 1./xp1**2*(9./8.-xp1/2.)
cc            c = 1.
cc            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
cc            inter = q4 + qt2*(q2-q4)
cc          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
cc            a = (xp2/4.-.5)/(xp2-1.)**3
cc            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
cc            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
cc            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
cc            inter = q4 + qt2*(q2-q4)
cc          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
cc            inter = 3.*q2/8. + 3.*q3/4.- q4/8.
cc          else
cc            inter = q3
cc          endif
          where (qt1.lt.xp1.and.qt1.ge.0.)
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            c = 1.
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
            inter = q4 + qt2*(q2-q4)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1.)
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q4 + qt2*(q2-q4)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
            inter = 3.*q2/8. + 3.*q3/4.- q4/8.
          elsewhere
            inter = q3
          end where
        endif

      case(6)

c     Centered high-order

        inter = (-q1 + 9.*q2 + 9.*q3 - q4)/16.

cc      case(7)
cc
ccc     Gamma [Jasak, Weller, and Gosman, Int. J. Numer. Meth. Fluids, 32, p431 (1999)]
cc
cc        slp1 = 1.5
cc        slp2 = .5
cc        if (vel.gt.0d0) then
cc          curv = ((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/0.5/(x3-x1)
cc          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
cccc          qt1 = 3.*q3/8. + 3.*q2/4.- q1/8.
cc          qt2 = slp1*q2 + (1.-slp1)*q1
cc          qt3 = slp2*q2 + (1.-slp2)*q3
cc          qt4 = fmed(q2,qt2,qt3)
cc          inter = fmed(q2,qt4,qt1)
cc        else
cc          curv = ((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/0.5/(x4-x2)
cc          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
cccc          qt1 = 3.*q2/8. + 3.*q3/4.- q4/8.
cc          qt2 = slp1*q3 + (1.-slp1)*q4
cc          qt3 = slp2*q3 + (1.-slp2)*q2
cc          qt4 = fmed(q3,qt2,qt3)
cc          inter = fmed(q3,qt4,qt1)
cc        endif

      case default

        call pstop('in2face_vec','Advection scheme not available')

      end select

c     End

      end function in2face_vec

      end module grid_advec_ops

c grid_diff_ops
c ######################################################################
      module grid_diff_ops

       use grid_vertex_ops

       use grid_nabla_ops

       use grid_nabla2_ops

       use grid_tensor_ops

       use grid_advec_ops

      contains

c     integral
c     ################################################################
      function integral(igx,igy,igz,nx,ny,nz,array,average,logical_dom)

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on domain (nx)x(ny)x(nz).
c     ---------------------------------------------------------------

      implicit none

c    Call variables

      integer :: igx,igy,igz,nx,ny,nz
      real(8) :: array(0:nx+1,0:ny+1,0:nz+1),integral
      logical,optional :: average,logical_dom

c     Local variables

      integer :: i,j,k

      real(8) :: volume,dvol,send_buf(2),rec_buf(2)

      logical :: avg,log_dom

c     Begin program

      avg = .false.
      if (PRESENT(average)) avg = average

      log_dom = .false.
      if (PRESENT(logical_dom)) log_dom = logical_dom

c     Integrate

      integral = 0d0
      volume   = 0d0

      if (log_dom) then
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              dvol = gmetric%grid(igx)%dlvol(i,j,k)

              if (isSYM(i,igx,dim=1,loc=0)
     .        .or.isSYM(i,igx,dim=1,loc=1)
     .        .or.isSYM(j,igy,dim=2,loc=0)
     .        .or.isSYM(j,igy,dim=2,loc=1)
     .        .or.isSYM(k,igz,dim=3,loc=0)
     .        .or.isSYM(k,igz,dim=3,loc=1)) dvol = 0.5*dvol

              integral = integral + array(i,j,k)*dvol
              volume = volume + dvol
            enddo
          enddo
        enddo
      else
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              dvol = gmetric%grid(igx)%dvol(i,j,k)

              if (isSYM(i,igx,dim=1,loc=0)
     .        .or.isSYM(i,igx,dim=1,loc=1)
     .        .or.isSYM(j,igy,dim=2,loc=0)
     .        .or.isSYM(j,igy,dim=2,loc=1)
     .        .or.isSYM(k,igz,dim=3,loc=0)
     .        .or.isSYM(k,igz,dim=3,loc=1)) dvol = 0.5*dvol

              integral = integral + array(i,j,k)*dvol
              volume = volume + dvol
            enddo
          enddo
        enddo
      endif

#if defined(petsc)
      send_buf = (/ integral, volume /)
      call MPI_Allreduce(send_buf,rec_buf,2,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      integral = rec_buf(1)
      volume = rec_buf(2)
#endif

      if (avg) integral = integral/volume

c     End 

      end function integral

      end module grid_diff_ops


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

c     curl_vrtx_mesh
c     ##################################################################
      function curl_vrtx_mesh(igrid,vv,vrtx) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curl of vv (covariant) from nearest neighbors. Returns
c     contravariant components.
c
c     If vrtx is true, it finds curl at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds curl at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid
      logical :: vrtx
      real(8) :: vv(0:,0:,0:,1:)
     $          ,cnv(0:size(vv,1)-1,0:size(vv,2)-1,0:size(vv,3)-1,3)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(vv,1)-2
      ny = size(vv,2)-2
      nz = size(vv,3)-2

      do k=0,nz
        do j=0,ny
          do i=0,nx
            cnv(i,j,k,:) = curl_vrtx(i,j,k,igrid,igrid,igrid
     .                              ,vv(i:i+1,j:j+1,k:k+1,:),vrtx)
          enddo
        enddo
      enddo

c     End program

      end function curl_vrtx_mesh

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

        INTERFACE div_stg
          module procedure div_stg_ijk,div_stg_mesh
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

      gx = 0d0 ; gy = 0d0 ; gz = 0d0

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
        if (i<nx+1) gx =(arr(i+1,j,k)-arr(i  ,j,k))/grid_params%dx(ig)
      case(-1)
        if (i>0   ) gx =(arr(i  ,j,k)-arr(i-1,j,k))/grid_params%dx(ig-1)
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
        if (j<ny+1) gy =(arr(i,j+1,k)-arr(i,j  ,k))/grid_params%dy(jg  )
      case(-1)
        if (j>0   ) gy =(arr(i,j  ,k)-arr(i,j-1,k))/grid_params%dy(jg-1)
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
        if (k<nz+1) gz =(arr(i,j,k+1)-arr(i,j,k  ))/grid_params%dz(kg  )
      case(-1)
        if (k>0   ) gz =(arr(i,j,k  )-arr(i,j,k-1))/grid_params%dz(kg-1)
      end select

      cov = (/ gx,gy,gz /)

c     End program

      end function grad_ijk

c     grad_mesh
c     ###############################################################
      function grad_mesh(igrid,arr,hex,hey,hez) result(cov)

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

      integer,optional :: hex,hey,hez

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            cov(i,j,k,:)=grad_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,arr
     .                           ,hex=hex,hey=hey,hez=hez)
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
      function curl_mesh(igrid,cov,hex,hey,hez) result(cnv)

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

      integer,optional :: hex,hey,hez

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            cnv(i,j,k,:)=curl_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,cov
     .                           ,hex=hex,hey=hey,hez=hez)
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

c     div_stg_ijk
c     ###############################################################
      function div_stg_ijk(i,j,k,nx,ny,nz,igx,igy,igz,cnv)
     .         result(div)

      implicit none

c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cnv(0:nx+1,0:ny+1,0:nz+1,3),div

c     Local variables

      integer :: ig,jg,kg,igrid,ip,im,jp,jm,kp,km,hx,hy,hz
      real(8) :: dxx,dyy,dzz

c     Begin program

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      im = i-1
      jm = j-1
      km = k-1

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

      div =  (cnv(i,j,k,1)-cnv(im,j ,k ,1))/dxx
     .      +(cnv(i,j,k,2)-cnv(i ,jm,k ,2))/dyy
     .      +(cnv(i,j,k,3)-cnv(i ,j ,km,3))/dzz

      div = div/gmetric%grid(igrid)%jac(i,j,k)
      
c     End 

      end function div_stg_ijk

c     div_stg_mesh
c     ###############################################################
      function div_stg_mesh(igrid,curv,iscnv) result(scl)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: curv(0:,0:,0:,:)
     .          ,scl(0:size(curv,1)-1,0:size(curv,2)-1,0:size(curv,3)-1)

      logical,optional :: iscnv

c     Local variables

      integer :: nx,ny,nz,i,j,k
      real(8) :: cnv(0:size(curv,1)-1,0:size(curv,2)-1,0:size(curv,3)-1
     .              ,size(curv,4))
      logical :: is_cnv

c     Begin program

      if (PRESENT(iscnv)) then
        is_cnv = iscnv
      else
        is_cnv = .true.
      endif

      scl = 0d0

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      if (is_cnv) then

        do k=1,nz
          do j=1,ny
            do i=1,nx
              scl(i,j,k)=div_stg_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                              ,curv)
            enddo
          enddo
        enddo

      else

        cnv = XformToCnv_stg_mesh(igrid,curv)

        do k=1,nz
          do j=1,ny
            do i=1,nx
              scl(i,j,k)=div_stg_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                              ,cnv)
            enddo
          enddo
        enddo

      endif

c     End program

      end function div_stg_mesh

c     curl_stg_mesh
c     ###############################################################
      function curl_stg_mesh(igr,cov) result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates for a
c     staggered vector A.  The vector A is covariant, and returns a
c     contravariant vector, also staggered.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igr

      real(8) :: cov(0:,0:,0:,:)
     .          ,cnv(0:size(cov,1)-1
     .              ,0:size(cov,2)-1
     .              ,0:size(cov,3)-1
     .                ,size(cov,4))

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      cnv = 0d0

      nx = size(cov,1)-2
      ny = size(cov,2)-2
      nz = size(cov,3)-2

      do k=0,nz
        do j=0,ny
          do i=0,nx
            cnv(i,j,k,:)=curl_stg_ijk(i,j,k,nx,ny,nz,igr,igr,igr,cov)
          enddo
        enddo
      enddo

c     End program

      end function curl_stg_mesh

c     curl_stg_ijk
c     ###############################################################
      function curl_stg_ijk(i,j,k,nx,ny,nz,igx,igy,igz,cov)
     .         result(curl)

      implicit none

c     ---------------------------------------------------------------
c     Calculates curl of staggered covariant vector field to find
c     staggered contravariant field in general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cov(0:nx+1,0:ny+1,0:nz+1,3),curl(3)

c     Local variables

      integer :: ig,jg,kg,igrid,ip,jp,kp,km,hx,hy,hz
      real(8) :: dxx,dyy,dzz

c     Begin program

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = min(i+1,nx+1)
      jp = min(j+1,ny+1)
      kp = min(k+1,nz+1)

      dxx = grid_params%dx(ig)
      dyy = grid_params%dy(jg)
      dzz = grid_params%dz(kg)

      curl(1) = (cov(i,jp,k ,3)-cov(i,j,k,3))/dyy
     .         -(cov(i,j ,kp,2)-cov(i,j,k,2))/dzz !@(i,j+1/2,k+1/2)
      curl(2) = (cov(i ,j,kp,1)-cov(i,j,k,1))/dzz
     .         -(cov(ip,j,k ,3)-cov(i,j,k,3))/dxx !@(i+1/2,j,k+1/2)
      curl(3) = (cov(ip,j ,k,2)-cov(i,j,k,2))/dxx
     .         -(cov(i ,jp,k,1)-cov(i,j,k,1))/dyy !@(i+1/2,j+1/2,k)
      
c     End 

      end function curl_stg_ijk

      end module grid_nabla_ops

c grid_nabla2_ops
c ######################################################################
      module grid_nabla2_ops

        use grid_create

        INTERFACE lap
          module procedure lap_mesh,lap_ijk
        end INTERFACE

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

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: vec(3)
      real(8) :: da(0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional,intent(IN) :: vol,cnv

c     Local variables

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8) :: divip,divim,divjp,divjm,divkp,divkm
     .          ,jacip,jacim,jacjp,jacjm,jackp,jackm
     .          ,jac

      logical :: vol_wgt,cnvt

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

      integer :: im,jm,km,ip,jp,kp,ig,jg,kg

      real(8) :: d_xx_ip
     .          ,d_xy_ipjp,d_xy_ipjm
     .          ,d_xz_ipkp,d_xz_ipkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      real(8),dimension(:,:,:),pointer :: dff

      logical  :: vrtx_met,have_diff

c     Begin program

      if (PRESENT(vertex)) then
        vrtx_met = vertex
      else
        vrtx_met = .true.
      endif

      have_diff = .false.
      if (PRESENT(diff)) then
        dff => diff
        have_diff = .true.
      endif

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

        if (have_diff) then
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

        if (have_diff) then
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

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8) :: d_yy_jp
     .          ,d_xy_ipjp,d_xy_imjp
     .          ,d_yz_jpkp,d_yz_jpkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup
      real(8),dimension(:,:,:),pointer :: dff

      logical  :: vrtx_met,have_diff

c     Begin program

      if (PRESENT(vertex)) then
        vrtx_met = vertex
      else
        vrtx_met = .true.
      endif

      have_diff = .false.
      if (PRESENT(diff)) then
        dff => diff
        have_diff = .true.
      endif

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
        else
          d_yy_jp = 0.25*(gsup(i ,j,k ,2,2)
     .                   +gsup(im,j,k ,2,2)
     .                   +gsup(i ,j,km,2,2)
     .                   +gsup(im,j,km,2,2))
        endif

        if (have_diff) then
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

        if (have_diff) then
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

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8) :: d_zz_kp
     .          ,d_xz_ipkp,d_xz_imkp
     .          ,d_yz_jpkp,d_yz_jmkp

      real(8),dimension(:,:,:,:,:),pointer :: gsup
      real(8),dimension(:,:,:),pointer :: dff

      logical  :: vrtx_met,have_diff

c     Begin program

      if (PRESENT(vertex)) then
        vrtx_met = vertex
      else
        vrtx_met = .true.
      endif

      have_diff = .false.
      if (PRESENT(diff)) then
        dff => diff
        have_diff = .true.
      endif

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

        if (have_diff) then
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

        if (have_diff) then
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

c     lap_ijk
c     ###############################################################
      function lap_ijk(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol
     .                ,vertex_metrics) result (laplacian)
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

      end function lap_ijk

c     lap_mesh
c     ###############################################################
      function lap_mesh(igr,arr,dff,vol,vertex_metrics) result (lap)
c     ---------------------------------------------------------------
c     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igr

      real(8) :: arr(0:,0:,0:)
     $          ,lap(0:size(arr,1)-1
     $              ,0:size(arr,2)-1
     $              ,0:size(arr,3)-1)

      real(8),optional :: dff(0:size(arr,1)-1
     $                       ,0:size(arr,2)-1
     $                       ,0:size(arr,3)-1)

      logical,optional,intent(IN) :: vol,vertex_metrics

c     Local variables

      integer :: i,j,k,nx,ny,nz

c     Begin program

      nx = size(arr,1)-2
      ny = size(arr,2)-2
      nz = size(arr,3)-2

      do k=1,nz
        do j=1,ny
          do i=1,nx
            lap(i,j,k) = lap_ijk(i,j,k,nx,ny,nz,igr,igr,igr,arr
     $              ,dff=dff,vol=vol,vertex_metrics=vertex_metrics)
          enddo
        enddo
      enddo

c     End program

      end function lap_mesh

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

        INTERFACE veclap
          module procedure veclap_mesh,veclap_ijk
        end INTERFACE

        INTERFACE veclap_cov
          module procedure veclap_cov_mesh,veclap_cov_ijk
        end INTERFACE

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
        real(8)    :: cnv(3),dum(3)

c     Begin program

        cnv = (/ ax(i,j,k),ay(i,j,k),az(i,j,k) /)

        dum = gmetric%grid(igrid)%Gamma(i,j,k,1,1,:)
     .       +gmetric%grid(igrid)%Gamma(i,j,k,2,2,:)
     .       +gmetric%grid(igrid)%Gamma(i,j,k,3,3,:)

        do l=1,3
          do m=1,3
            tensor(l,m) =
     .           dot_product(cnv,gmetric%grid(igrid)%Gamma(i,j,k,m,:,l))
     .          -cnv(m)*dum(l)
cc     .           cnv(1)* gmetric%grid(igrid)%Gamma(i,j,k,m,1,l)
cc     .          +cnv(2)* gmetric%grid(igrid)%Gamma(i,j,k,m,2,l)
cc     .          +cnv(3)* gmetric%grid(igrid)%Gamma(i,j,k,m,3,l)
cc     .          -cnv(m)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,l)
cc     .                  +gmetric%grid(igrid)%Gamma(i,j,k,2,2,l)
cc     .                  +gmetric%grid(igrid)%Gamma(i,j,k,3,3,l))
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

c     veclap_ijk
c     ###############################################################
      function veclap_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                   ,diff,alteom,vol) result (vlap)

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

      end function veclap_ijk

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

cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
        if (isSP2(i+1,igx).and. flag /= 0) then
          ijac = 0d0
        else
          ijac = 1d0/jac
        endif

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

c     veclap_cov_mesh
c     ###############################################################
      function veclap_cov_mesh(igrid,cov,diff) result(cov_out)

c     ---------------------------------------------------------------
c     Calculates covariant components of lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cov_out(0:size(cov,1)-1
     $                  ,0:size(cov,2)-1
     $                  ,0:size(cov,3)-1,3)

      real(8),target,optional,dimension(:,:,:) :: diff

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

cc#if defined(samrai)
cc      real(8),dimension(:,:,:),pointer :: dff
cc      if (PRESENT(diff)) then
cc        dff => diff
cc      endif
cc#endif

      nx = size(cov,1)-2
      ny = size(cov,2)-2
      nz = size(cov,3)-2

cc#if defined(samrai)
cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc            cov_out(i,j,k,:)=veclap_cov_ijk(i,j,k,nx,ny,nz
cc     .                                     ,igrid,igrid,igrid
cc     .                                     ,cov,dff=dff)
cc          enddo
cc        enddo
cc      enddo
cc#else
      if (PRESENT(diff)) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              cov_out(i,j,k,:)=veclap_cov_ijk(i,j,k,nx,ny,nz
     .                                       ,igrid,igrid,igrid
     .                                       ,cov,dff=diff)
            enddo
          enddo
        enddo
      else
        do k=1,nz
          do j=1,ny
            do i=1,nx
              cov_out(i,j,k,:)=veclap_cov_ijk(i,j,k,nx,ny,nz
     .                                       ,igrid,igrid,igrid,cov)
            enddo
          enddo
        enddo
      endif
cc#endif

c     End program

      end function veclap_cov_mesh

c     veclap_cov_ijk
c     ###############################################################
      function veclap_cov_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vfield,dff,vol)
     $         result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target,optional :: dff(0:nx+1,0:ny+1,0:nz+1)

      real(8) :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: icomp

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      if (no_map) then
        do icomp=1,3
          vlap(icomp)=lap(i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vfield(:,:,:,icomp),dff=dff
     .                   ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        if (PRESENT(dff)) coef => dff   !Pointer passed to nabtensor routines
        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
        if (PRESENT(dff)) nullify(coef)
      endif

c     End program

      end function veclap_cov_ijk

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

      cnv = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
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

        integer,private :: bc_adv_meth=2

        !Smooth SMART config
        real(8),private :: ssmart_xp1=0.25 !0.3
     .                    ,ssmart_xp2=0.52 !0.5

        !SMART config
        real(8),private :: smart_slp1= 1.5  !3.0
     .                    ,smart_slp2= 0.5  !0.3

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
     .                       ,method,sp,zip_vel) result(flxip)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_1 of scalar arr at face
c     (i+1/2,j,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxip

      logical,optional :: sp,zip_vel

c     Local variables

      integer :: im,ip,ipp,igm,igp,igpp,ig,jg,kg,admth

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
          admth = bc_adv_meth
        else
          admth = method
        endif

        qip=in2face_scl(grid_params%xx(igpp),arr(ipp,j,k)
     .                 ,grid_params%xx(igp ),arr(ip ,j,k)
     .                 ,grid_params%xx(ig  ),arr(i  ,j,k)
     .                 ,grid_params%xx(igm ),arr(im ,j,k)
     .                 ,vxp,admth)

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
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxjp

      logical,optional :: zip_vel

c     Local variables

      integer :: jm,jp,jpp,jgm,jgp,jgpp,ig,jg,kg,admth

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
          admth = bc_adv_meth
        else
          admth = method
        endif

        qjp=in2face_scl(grid_params%yy(jgpp),arr(i,jpp,k)
     .                 ,grid_params%yy(jgp ),arr(i,jp ,k)
     .                 ,grid_params%yy(jg  ),arr(i,j  ,k)
     .                 ,grid_params%yy(jgm ),arr(i,jm ,k)
     .                 ,vyp,admth)

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
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxkp

      logical,optional :: zip_vel

c     Local variables

      integer :: km,kp,kpp,kgm,kgp,kgpp,ig,jg,kg,admth

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
          admth = bc_adv_meth
        else
          admth = method
        endif

        qkp=in2face_scl(grid_params%zz(kgpp),arr(i,j,kpp)
     .                 ,grid_params%zz(kgp ),arr(i,j,kp )
     .                 ,grid_params%zz(kg  ),arr(i,j,k  )
     .                 ,grid_params%zz(kgm ),arr(i,j,km )
     .                 ,vzp,admth)

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
c     method: specifies the advective scheme (see int2face)
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

c     flx_nc_advec
c     ###############################################################
      function flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                     ,method,sp,vol,reverse) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1,3),advec

      logical,optional :: sp,vol,reverse

c     Local variables

      integer :: igrid,half_elem,su,bc_meth

      integer :: im,ip,ipp,imm,igm,igp,igpp,igmm
     .          ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
     .          ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
     .          ,ig,jg,kg,igl,jgl,kgl

      real(8) :: dxx,dyy,dzz,x0,y0,z0,jac,sgn,vel
      real(8) :: vxp,vxm,vyp,vym,vzp,vzm
     .          ,qip,qim,qjp,qjm,qkp,qkm
      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical :: vw,spoint,bc_safe,rev_sgn

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

      bc_meth = 2  !Advective method at boundaries

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
     .             ,vxp,bc_meth)
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
     .             ,vxm,bc_meth)
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
     .             ,vyp,bc_meth)
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
     .             ,vym,bc_meth)
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
     .             ,vzp,bc_meth)
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
     .             ,vzm,bc_meth)
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
c     method: specifies the advective scheme (see int2face)
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

        integer :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8) :: t11,t12,t13
        logical :: alt_eom

c     Local variables

        integer :: ig,jg,kg,ip,im,ipp,igpp,igp,igm,icomp,adv__meth
        real(8) :: vxp,qp(3),t1(3)
        real(8) :: jac,jac0,jacp,ijac,ijac0,ijacp

        logical :: bc_safe

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

        bc_safe = .not.(adv_meth > 2 .and. adv_meth <= 6)

        if (adv_meth /= 2 .and. flag /= 0) then

          !Revert extended support methods to upwinding at processor boundaries
          if ((i == nx .or. i == 0).and.(.not.bc_safe))then
            adv__meth = bc_adv_meth
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

        logical :: bc_safe

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

        bc_safe = .not.(adv_meth > 2 .and. adv_meth <= 6)

        if (adv_meth /= 2 .and. flag /= 0) then

          !Revert extended support methods to upwinding at processor boundaries
          if ((j == ny .or. j == 0).and.(.not.bc_safe))then
            adv__meth = bc_adv_meth
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

        logical :: bc_safe

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

        bc_safe = .not.(adv_meth > 2 .and. adv_meth <= 6)

        if (adv_meth /= 2 .and. flag /= 0) then

         !Revert extended support methods to upwinding at processor boundaries
          if ((k == nz .or. k == 0).and.(.not.bc_safe))then
            adv__meth = bc_adv_meth
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
          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
        else
          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
        endif

        inter = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2

      case(4)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = smart_slp1
        slp2 = smart_slp2
        if (vel.gt.0d0) then
          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q2 + (1.-slp1)*q1
          qt3 = slp2*q2 + (1.-slp2)*q3
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q3 + (1.-slp1)*q4
          qt3 = slp2*q3 + (1.-slp2)*q2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(5)

c     Smooth SMART

        xp1 = ssmart_xp1
        xp2 = ssmart_xp2

        if (vel.gt.0d0) then
          qt1  = (q2-q1)/(q3-q1)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1d0)
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp2.and.qt1.lt.1d0) then
            a = (0.25d0*xp2-.5d0)/(xp2-1d0)**3
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 0.375d0*q3 + 0.75d0*q2 - 0.125d0*q1
          else
            inter = q2
          endif
        else
          qt1  = (q3-q4)/(q2-q4)
          if (qt1.lt.xp1.and.qt1.ge.0d0) then
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            c = 1d0
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp2.and.qt1.lt.1d0) then
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 0.375d0*q2 + 0.75d0*q3- 0.125d0*q4
          else
            inter = q3
          endif
        endif

      case(6)

c     Centered high-order

        inter = 0.0625d0*(-q1 + 9.*q2 + 9.*q3 - q4)

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

        slp1 = smart_slp1
        slp2 = smart_slp2

        if (vel.gt.0d0) then
          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q2 + (1.-slp1)*q1
          qt3 = slp2*q2 + (1.-slp2)*q3
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
          qt2 = slp1*q3 + (1.-slp1)*q4
          qt3 = slp2*q3 + (1.-slp2)*q2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(5)

c     Smooth SMART

        xp1 = ssmart_xp1
        xp2 = ssmart_xp2

        if (vel.gt.0d0) then
          qt1  = (q2-q1)/(q3-q1)
          where (qt1.lt.xp1.and.qt1.ge.0d0)
            a = .25*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-.5d0*xp1)*xp1**(-2)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1d0)
            inter = q1 + qt2*(q3-q1)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1d0)
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
            inter = q1 + qt2*(q3-q1)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
            inter = 0.375d0*q3 + 0.75d0*q2- 0.125d0*q1
          elsewhere
            inter = q2
          end where
        else
          qt1  = (q3-q4)/(q2-q4)
          where (qt1.lt.xp1.and.qt1.ge.0d0)
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            c = 1d0
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
            inter = q4 + qt2*(q2-q4)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1d0)
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
            inter = q4 + qt2*(q2-q4)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
            inter = 0.375d0*q2 + 0.75d0*q3 - 0.125d0*q4
          elsewhere
            inter = q3
          end where
        endif

      case(6)

c     Centered high-order

        inter = 0.0625d0*(-q1 + 9.*q2 + 9.*q3 - q4)

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

c grid_pardiff_ops
c ######################################################################
      module grid_pardiff_ops

        use grid_create

      contains

c     par_diff
c     ################################################################
      function par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,par_order,symmetric,limited,vol) result(par2)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 2nd or 4th order
c     accuracy.
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz,par_order
      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1),par2
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)

      logical,optional :: symmetric,limited,vol

c     Local variables

      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      real(8)  :: bx2,by2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
      integer  :: ii,jj,kk,ig,jg,kg
      logical  :: limit,symm,vol_wgt

      real(8),dimension(:,:,:),pointer :: jac

c     Begin program

      Fluxx=0
      Fluxy=0
      Fluxz=0

      if (PRESENT(limited)) then
        limit = limited
      else
        limit = .false.
      endif

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      if (PRESENT(symmetric)) then
        symm = symmetric
      else
        symm = .false.
      endif

      jac => gmetric%grid(igx)%jac

c     Find fluxes

      select case(par_order)
      case(2) !second order
        call par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,Fluxx,Fluxy,Fluxz,symmetric=symm,limited=limit)
      case(4) !fourth order, 5-pt stencil in co-directed dervs.
        call par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,Fluxx,Fluxy,Fluxz,limited=limit,stencil=5)
      case(-4) !fourth order, 7-pt stencil in co-directed dervs.
        call par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,Fluxx,Fluxy,Fluxz,limited=limit,stencil=7)
      case default
        write (*,*) 'Order not implemented'
        stop
      end select

c     Assemble fluxes

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      par2 = (Fluxx(i,j,k)-Fluxx(i-1,j,k))/grid_params%dxh(ig)
     .      +(Fluxy(i,j,k)-Fluxy(i,j-1,k))/grid_params%dyh(jg)
     .      +(Fluxz(i,j,k)-Fluxz(i,j,k-1))/grid_params%dzh(kg)

      par2 = par2/jac(i,j,k)

      if (vol_wgt) par2=par2*gmetric%grid(igx)%dvol(i,j,k)

      end function par_diff

c     par_flux_2nd
c     ################################################################
      subroutine par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                       ,Fluxx,Fluxy,Fluxz,symmetric,limited)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 2nd order accuracy
c     ----------------------------------------------------------------

cc      use transport_params

      implicit none

c     Call variables

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      logical,optional :: limited,symmetric

c     Local variables

      real(8)  :: bx2,by2,bz2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
     .           ,jac_avg
      integer  :: ig,jg,kg,ii,jj,kk
      logical  :: limit,symm

      real(8),dimension(:,:,:),pointer :: jac
      real(8),dimension(:,:,:,:,:),pointer :: gsup

c     Begin program

      Fluxx=0
      Fluxy=0
      Fluxz=0

      if (PRESENT(limited)) then
        limit = limited
      else
        limit = .false.
      endif

      if (PRESENT(symmetric)) then
        symm = symmetric
      else
        symm = .false.
      endif

      jac  => gmetric%grid(igx)%jac
      gsup => gmetric%grid(igx)%gsup

c     Symmetric discretization

      if (symm) then

        do kk = k-1,k
          do jj = j-1,j
            do ii = i-1,i
              call getMGmap(ii,jj,kk,igx,igy,igz,ig,jg,kg)

              jac_avg = 0.25*(jac(ii  ,jj  ,kk)
     .                       +jac(ii+1,jj  ,kk)
     .                       +jac(ii  ,jj+1,kk)
     .                       +jac(ii+1,jj+1,kk))

              bx2 = 0d0
              by2 = 0d0
              bxby = 0d0

              if (jac_avg /= 0d0) then
                bx2 = (0.25*(bx(ii  ,jj  ,kk)
     .                      +bx(ii+1,jj  ,kk)
     .                      +bx(ii  ,jj+1,kk)
     .                      +bx(ii+1,jj+1,kk)))**2/jac_avg

                bxby= (.25*(bx(ii  ,jj  ,kk)
     .                     +bx(ii+1,jj  ,kk)
     .                     +bx(ii  ,jj+1,kk)
     .                     +bx(ii+1,jj+1,kk)))
     .               *(.25*(by(ii  ,jj  ,kk)
     .                     +by(ii+1,jj  ,kk)
     .                     +by(ii  ,jj+1,kk)
     .                     +by(ii+1,jj+1,kk)))/jac_avg

                by2 = (0.25*(by(ii  ,jj  ,kk)
     .                      +by(ii+1,jj  ,kk)
     .                      +by(ii  ,jj+1,kk)
     .                      +by(ii+1,jj+1,kk)))**2/jac_avg
              endif

              dtx = .5*(tmp(ii+1,jj  ,kk)-tmp(ii,jj  ,kk)
     .                 +tmp(ii+1,jj+1,kk)-tmp(ii,jj+1,kk))
     .                 /grid_params%dx(ig)
              dty = .5*(tmp(ii  ,jj+1,kk)-tmp(ii  ,jj,kk)
     .                 +tmp(ii+1,jj+1,kk)-tmp(ii+1,jj,kk))
     .                 /grid_params%dy(jg)

              Fluxx(ii,jj,kk) = bx2*dtx + bxby*dty
              Fluxy(ii,jj,kk) = by2*dty + bxby*dtx
              Fluxz(ii,jj,kk) = 0d0
            enddo
          enddo
        enddo

        !Average fluxes from vertices to faces
        Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
     .                      +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))

        Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
     .                      +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))

        Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
     .                      +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))

c     Asymmetric discretization

      else

        if (limit) then !Limited differences (Sharma, Hammett, JCP 227 (2007))

          !X-component cross flux
          do ii = i-1,i
            call getMGmap(ii,j,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)
     .                 *by(ii  ,j,k)
     .                 +bx(ii+1,j,k)/jac(ii+1,j,k)
     .                 *by(ii+1,j,k))

            flx0 = mc(tmp(ii  ,j  ,k)-tmp(ii  ,j-1,k)
     .               ,tmp(ii  ,j+1,k)-tmp(ii  ,j  ,k))
            flxp = mc(tmp(ii+1,j  ,k)-tmp(ii+1,j-1,k)
     .               ,tmp(ii+1,j+1,k)-tmp(ii+1,j  ,k))

            Fluxx(ii,j,k)= bxby/grid_params%dyh(jg)*mc(flx0,flxp)
          enddo

          !Y-component cross flux 
          do jj = j-1,j
            call getMGmap(i,jj,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)
     .                 *by(i,jj  ,k)
     .                 +bx(i,jj+1,k)/jac(i,jj+1,k)
     .                 *by(i,jj+1,k))

            flx0 = mc(tmp(i  ,jj  ,k)-tmp(i-1,jj  ,k)
     .               ,tmp(i+1,jj  ,k)-tmp(i  ,jj  ,k))
            flxp = mc(tmp(i  ,jj+1,k)-tmp(i-1,jj+1,k)
     .               ,tmp(i+1,jj+1,k)-tmp(i  ,jj+1,k))

            Fluxy(i,jj,k)= bxby/grid_params%dxh(ig)*mc(flx0,flxp)
          enddo

        else !Standard central average

          !X-component cross flux
          do ii = i-1,i
            call getMGmap(ii,j,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)*by(ii  ,j,k)
     .                 +bx(ii+1,j,k)/jac(ii+1,j,k)*by(ii+1,j,k))

            Fluxx(ii,j,k)= 0.25/grid_params%dyh(jg)*bxby
     .             *(tmp(ii  ,j+1,k)-tmp(ii  ,j-1,k)
     .              +tmp(ii+1,j+1,k)-tmp(ii+1,j-1,k))
          enddo

          !Y-component cross flux 
          do jj = j-1,j
            call getMGmap(i,jj,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)*by(i,jj  ,k)
     .                 +bx(i,jj+1,k)/jac(i,jj+1,k)*by(i,jj+1,k))

            Fluxy(i,jj,k)= 0.25/grid_params%dxh(ig)*bxby
     .                   *(tmp(i+1,jj  ,k)-tmp(i-1,jj  ,k)
     .                    +tmp(i+1,jj+1,k)-tmp(i-1,jj+1,k))
          enddo

        endif

        !Add co-derivative terms: standard centered
        do ii = i-1,i
          if (isSP2(ii+1,igx)) cycle

          call getMGmap(ii,j,k,igx,igy,igz,ig,jg,kg)

cc          bx2 = (0.5*(bx(ii+1,j,k)/jac(ii+1,j,k)
cc     .               +bx(ii  ,j,k)/jac(ii  ,j,k)))**2

          jac_avg = 0.5*(jac(ii+1,j,k)+jac(ii,j,k))
          bx2 = (0.5*(bx(ii+1,j,k)+bx(ii  ,j,k)))**2/jac_avg

          Fluxx(ii,j,k) = Fluxx(ii,j,k)
     .                  + bx2*(tmp(ii+1,j,k)-tmp(ii,j,k))
     .                       /grid_params%dx(ig)
        enddo

        do jj = j-1,j
          call getMGmap(i,jj,k,igx,igy,igz,ig,jg,kg)
 
cc          by2 = (0.5*(by(i,jj+1,k)/jac(i,jj+1,k)
cc     .               +by(i,jj  ,k)/jac(i,jj  ,k)))**2

          jac_avg = 0.5*(jac(i,jj+1,k)+jac(i,jj,k))
          by2 = (0.5*(by(i,jj+1,k)+by(i,jj,k)))**2/jac_avg

          Fluxy(i,jj,k) = Fluxy(i,jj,k)
     .                  + by2*(tmp(i,jj+1,k)-tmp(i,jj,k))
     .                       /grid_params%dy(jg)
        enddo

        do kk = k-1,k
          call getMGmap(i,j,kk,igx,igy,igz,ig,jg,kg)

cc          bz2 = (0.5*(bz(i,j,kk+1)/jac(i,j,kk+1)
cc     .               +bz(i,j,kk  )/jac(i,j,kk  )))**2

          jac_avg = 0.5*(jac(i,j,kk+1)+jac(i,j,kk))
          bz2 = (0.5*(bz(i,j,kk+1)+bz(i,j,kk)))**2/jac_avg

          Fluxz(i,j,kk) = Fluxz(i,j,kk)
     .                  + bz2*(tmp(i,j,kk+1)-tmp(i,j,kk))
     .                       /grid_params%dz(kg)
        enddo

cc      !Add co-derivative terms: ZIP
cc      do i = 0,nx
cc        call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz) 
cc     .              +bx(i+1,1:ny,1:nz)*bx(i,1:ny,1:nz)
cc     .                   *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cc      enddo
cc
cc      do j = 0,ny
cc        call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz) 
cc     .              +by(1:nx,j+1,1:nz)*by(1:nx,j,1:nz)
cc     .                   *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cc      enddo
cc
cc      do k = 0,nz
cc        call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
cc     .              +bz(1:nx,1:ny,k+1)*bz(1:nx,1:ny,k)
cc     .                   *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cc      enddo

      endif

c     End program

      end subroutine par_flux_2nd

ccc     iso_flux_2nd
ccc     ################################################################
cc      subroutine iso_flux_2nd(nx,ny,nz,igx,igy,igz,tmp
cc     .                       ,Fluxx,Fluxy,Fluxz,symmetric)
ccc     ----------------------------------------------------------------
ccc     Calculates div.(grad) in 3D with 2nd order accuracy
ccc     ----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer  :: nx,ny,nz,igx,igy,igz
cc      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1)
cc      real(8)  :: Fluxx(0:nx,0:ny,0:nz)
cc     .           ,Fluxy(0:nx,0:ny,0:nz)
cc     .           ,Fluxz(0:nx,0:ny,0:nz)
cc
cc      logical,optional :: symmetric
cc
ccc     Local variables
cc
cc      real(8)  :: dtx,dty,dtz,a1,a2,sgn,bb,flx0,flxp
cc      integer  :: i,j,k,ig,jg,kg
cc      logical  :: limit,symm
cc
ccc     Begin program
cc
cc      Fluxx=0
cc      Fluxy=0
cc      Fluxz=0
cc
cc      if (PRESENT(symmetric)) then
cc        symm = symmetric
cc      else
cc        symm = .false.
cc      endif
cc
ccc     Symmetric discretization
cc
cc      if (symm) then
cc
cc        do k = 0,nz
cc          do j = 0,ny
cc            do i = 0,nx
cc
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              dtx = .5*(tmp(i+1,j  ,k)-tmp(i,j  ,k)
cc     .                 +tmp(i+1,j+1,k)-tmp(i,j+1,k))/dx(ig)
cc              dty = .5*(tmp(i  ,j+1,k)-tmp(i  ,j,k)
cc     .                 +tmp(i+1,j+1,k)-tmp(i+1,j,k))/dy(jg)
cc              dtz = 0d0
cc
cc              Fluxx(i,j,k)= dtx
cc              Fluxy(i,j,k)= dty
cc              Fluxz(i,j,k)= dtz
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Average fluxes to faces
cc        do k = nz,1,-1
cc          do j = ny,1,-1
cc            Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
cc     .                          +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))
cc          enddo
cc        enddo
cc
cc        do k = nz,1,-1
cc          do i = nx,1,-1
cc            Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
cc     .                          +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))
cc          enddo
cc        enddo
cc
cc        do j = ny,1,-1
cc          do i = nx,1,-1
cc            Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
cc     .                          +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))
cc          enddo
cc        enddo
cc
ccc     Asymmetric discretization
cc
cc      else
cc
cc        !Add co-derivative terms: standard centered
cc        do i = 0,nx
cc          call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxx(i,1:ny,1:nz) =
cc     .         (tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cc        enddo
cc
cc        do j = 0,ny
cc          call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxy(1:nx,j,1:nz) =
cc     .         (tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cc        enddo
cc
cc        do k = 0,nz
cc          call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxz(1:nx,1:ny,k) =
cc     .         (tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cc        enddo
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine iso_flux_2nd

c     par_flux-4th
c     ################################################################
      subroutine par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,F,Bx,By,Bz
     .                       ,Fluxx,Fluxy,Fluxz,limited,stencil)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 4th order accuracy
c     Can deal with 1D, 2D, and 3D.
c     Checks if there are enough points for the stencil in all directions
c     CANNOT handle different dx, dy, dz, (except for reduced dimensions)
c     Need information about MG to do that, hence for now dx=dy=dz=h
c
c     Coded by Natalia Krashenninikova, 8/2008
c     Modified by L. Chacon, 9/2008
c     ----------------------------------------------------------------
      implicit none

c     Call variables

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)  :: F (0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      integer,optional :: stencil
      logical,optional :: limited

c     Local variables

      real(8)  :: dx,dy,dz
      integer  :: tgx,tgy,tgz,xi,xf,yi,yf,zi,zf,ig,jg,kg,ii,jj,kk,stncl
      logical  :: limit

      real(8),dimension(:,:,:),pointer :: jac

c     Begin program

      if (PRESENT(limited)) then
        limit = limited
      else
        limit = .false.
      endif

      if (PRESENT(stencil)) then
        stncl = stencil
      else
        stncl = 5
      endif

      Fluxx=0
      Fluxy=0
      Fluxz=0

c     Get (constant) mesh spacings

      call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)

      dx=grid_params%dxh(ig)
      dy=grid_params%dyh(jg)
      dz=grid_params%dzh(kg)

      jac => gmetric%grid(igx)%jac

c     Calculate fluxes

      !X flux
      if (nx > 1) then
        do ii=i-1,i
          tgx=IOR(min(ii-2,0),max(ii+2-nx,0))
          if (ii == 0) tgx = -2
          xi=max(ii-2,0)-max(ii+3-(nx+1),0)
          xf=min(ii+3,nx+1)+max(2-ii,0)
          Fluxx(ii,j,k)=Get_Flux(dx,tgx,F(xi:xf,j,k)
     .              ,Bx(xi:xf,j,k)*Bx(xi:xf,j,k)/jac(xi:xf,j,k))

          if (ny > 1) then
            tgx=IOR(min(ii-2,0),max(ii+2-(nx+1),0))
            tgy=IOR(min(j-2,0),max(j+2-(ny+1),0))
            if (ii == 0) tgx = -2
            xi=max(ii-2,0)-max(ii+2-(nx+1),0)
            xf=min(ii+2,nx+1)+max(2-ii,0)
            yi=max(j-2,0)-max(j+2-(ny+1),0)
            yf=min(j+2,ny+1)+max(2-j,0)
            Fluxx(ii,j,k)=Fluxx(ii,j,k)
     .                   +Get_Flux_Cross(dy,tgx,tgy,F(xi:xf,yi:yf,k)
     .                   ,Bx(xi:xf,j,k)*By(xi:xf,j,k)/jac(xi:xf,j,k))
          endif

          if (nz > 1) then
            tgx=IOR(min(ii-2,0),max(ii+2-(nx+1),0))
            if (ii == 0) tgx = -2
            tgz=IOR(min(k-2,0),max(k+2-(nz+1),0))
            xi=max(ii-2,0)-max(ii+2-(nx+1),0)
            xf=min(ii+2,nx+1)+max(2-ii,0)
            zi=max(k-2,0)-max(k+2-(nz+1),0)
            zf=min(k+2,nz+1)+max(2-k,0)

            Fluxx(ii,j,k)=Fluxx(ii,j,k)
     .                  +Get_Flux_Cross(dz,tgx,tgz,F(xi:xf,j,zi:zf)
     .                  ,Bx(xi:xf,j,k)*Bz(xi:xf,j,k)/jac(xi:xf,j,k))
         endif

        enddo
      endif

      !Y flux
      if (ny > 1) then
        do jj=j-1,j
          tgy=IOR(min(jj-2,0),max(jj+2-ny,0))
          if (jj == 0) tgy = -2
          yi=max(jj-2,0)-max(jj+3-(ny+1),0)
          yf=min(jj+3,ny+1)+max(2-jj,0)
          Fluxy(i,jj,k)=Get_Flux(dy,tgy,F(i,yi:yf,k)
     .                    ,By(i,yi:yf,k)*By(i,yi:yf,k)/jac(i,yi:yf,k))

          if (nx > 1) then
            tgx=IOR(min(i-2,0),max(i+2-(nx+1),0))
            tgy=IOR(min(jj-2,0),max(jj+2-(ny+1),0))
            if (jj == 0) tgy = -2
            xi=max(i-2,0)-max(i+2-(nx+1),0)
            xf=min(i+2,nx+1)+max(2-i,0)
            yi=max(jj-2,0)-max(jj+2-(ny+1),0)
            yf=min(jj+2,ny+1)+max(2-jj,0)
            Fluxy(i,jj,k)=Fluxy(i,jj,k)
     .                  +Get_Flux_Cross(dx,tgy,tgx
     .                     ,Transpose(F(xi:xf,yi:yf,k))
     .                     ,By(i,yi:yf,k)*Bx(i,yi:yf,k)/jac(i,yi:yf,k))
          endif

          if (nz > 1) then
            tgy=IOR(min(jj-2,0),max(jj+2-(ny+1),0))
            if (jj == 0) tgy = -2
            tgz=IOR(min(k-2,0),max(k+2-(nz+1),0))
            yi=max(jj-2,0)-max(jj+2-(ny+1),0)
            yf=min(jj+2,ny+1)+max(2-jj,0)
            zi=max(k-2,0)-max(k+2-(nz+1),0)
            zf=min(k+2,nz+1)+max(2-k,0)
            Fluxy(i,jj,k)=Fluxy(i,jj,k)
     .                   +Get_Flux_Cross(dz,tgy,tgz,F(i,yi:yf,zi:zf)
     .                     ,By(i,yi:yf,k)*Bz(i,yi:yf,k)/jac(i,yi:yf,k))
          endif
        enddo
      endif

      !Z flux
      if (nz > 1) then
        do kk=k-1,k
          tgz=IOR(min(kk-2,0),max(kk+2-nz,0))
          if (kk == 0) tgz = -2
          zi=max(kk-2,0)-max(kk+2-(nz+1),0)
          zf=min(kk+2,nz+1)+max(2-kk,0)
          Fluxz(i,j,kk)=Get_Flux(dz,tgz,F(i,j,zi:zf)
     .                    ,Bz(i,j,zi:zf)*Bz(i,j,zi:zf)/jac(i,j,zi:zf))

          if (ny > 1) then
            tgy=IOR(min(j-2,0),max(j+2-(ny+1),0))
            tgz=IOR(min(kk-2,0),max(kk+2-(nz+1),0))
            if (kk == 0) tgz = -2
            yi=max(j-2,0)-max(j+2-(ny+1),0)
            yf=min(j+2,ny+1)+max(2-j,0)
            zi=max(kk-2,0)-max(kk+2-(nz+1),0)
            zf=min(kk+2,nz+1)+max(2-kk,0)
            Fluxz(i,j,kk)=Fluxz(i,j,kk)
     .                  +Get_Flux_Cross(dy,tgz,tgy
     .                     ,Transpose(F(i,yi:yf,zi:zf))
     .                     ,Bz(i,j,zi:zf)*By(i,j,zi:zf)/jac(i,j,zi:zf))
          endif

          if (nx > 1) then
            tgx=IOR(min(i-2,0),max(i+2-(nx+1),0))
            tgz=IOR(min(kk-2,0),max(kk+2-(nz+1),0))
            if (kk == 0) tgz = -2
            xi=max(i-2,0)-max(i+2-(nx+1),0)
            xf=min(i+2,nx+1)+max(2-i,0)
            zi=max(kk-2,0)-max(kk+2-(nz+1),0)
            zf=min(kk+2,nz+1)+max(2-kk,0)
            Fluxz(i,j,kk)=Fluxz(i,j,kk)
     .                  +Get_Flux_Cross(dx,tgz,tgx
     .                     ,Transpose(F(xi:xf,j,zi:zf))
     .                     ,Bz(i,j,zi:zf)*Bx(i,j,zi:zf)/jac(i,j,zi:zf))
          endif
        enddo
      endif

      contains

c     ################################################################
      function Get_Flux_Cross(h,tagx,tagy,F,B2) result(val)
c     ----------------------------------------------------------------
c     Calculates the flux for d^2/dx dy part of div.(BB.grad) operator 
c     with 4th order accuracy at a point.
c     Can handle boundary points
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: tagx,tagy
      real(8)  :: F(5,5),B2(5),val,h

c     Local variables

      integer  :: a(-1:1,5),cm(4),cl(5),cr(5)
      real(8)  :: twelve=12d0,bxby,flx0,flx1,flx
      logical  :: chk0,chk1,chk2

      val=0

      a(-1,:)=(/-3,-10,18,-6,1/)
      a(0,:)=(/1,-8,0,8,-1/)
      a(1,:)=(/-1,6,-18,10,3/)

      cm=(/-1,7,7,-1/)
      cl=(/2,17,-11,5,-1/)
      cr=(/-1,5,-11,17,2/)

      if((tagx+2)*(tagx+1)*tagx*(tagx-1)+(tagy+1)*tagy*(tagy-1)==0)then
         select case(tagx)
         case (-2)
            val=sum(cl*B2*matmul(F,a(tagy,:))/twelve/h)/twelve
         case (-1)
            val=sum(cm*B2(1:4)*matmul(F(1:4,:),a(tagy,:))
     .           /twelve/h)/twelve
         case (0)
            val=sum(cm*B2(2:5)*matmul(F(2:5,:),a(tagy,:))
     .           /twelve/h)/twelve
         case (1)
            val=sum(cr*B2*matmul(F,a(tagy,:))/twelve/h)/twelve
         end select
      else
         write(*,*)"Error in Get_Flux_Cross, wrong tags"
         write(*,*)"tagx=",tagx,"tagy=",tagy
         write(*,*)"Choices are -1, 0, or 1 for left boundary, "
         write(*,*)"interior or right boundary points"
         stop
      endif

      !Check for monotonicity (Sharma, Hammett, JCP 227 (2007))
      if (limit) then
        bxby = 0.5*(B2(4)+B2(3))

        flx0 = mc(F(3,3)-F(3,2),F(3,4)-F(3,3),chk=chk0)
        flx1 = mc(F(4,3)-F(4,2),F(4,4)-F(4,3),chk=chk1)
        flx  = mc(flx0         ,flx1         ,chk=chk2)

cc        write (*,*) 'here',bxby,chk0,chk1,chk2
      
        if (chk0 .or. chk1 .or. chk2) then
          val = bxby/h*flx
cc          write (*,*) 'here'
cc          val = 0d0
        endif
      endif

      end function Get_Flux_Cross

c     ################################################################
      function Get_Flux(h,tag,F,B2) result (val)
c     ----------------------------------------------------------------
c     Calculates the flux for d^2/dx^2 part of div.(BB.grad) operator 
c     with 4th order accuracy at a point.
c     Can handle boundary points
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: tag
      real(8)  :: F(6),B2(6),val,h

c     Local variables
cc      integer  :: a(5,6)
      real(8)  :: sv2=72d0,shty=720d0
      integer  :: a(0:5,6),cm(4),cl(5),cr(5)
      real(8)  :: twelve=12d0,sixty=60d0

      select case(stncl)
      case(5)

        val=0
        a=0
        a(1,2:5)=(/0,1,4,-5/)
        a(2,2:5)=(/11,-54,33,10/)
        a(3,2:5)=(/-10,-33,54,-11/)
        a(4,2:5)=(/5,-4,-1,0/)

        select case(tag)
        case (-2)
           a(1,:)=(/-36, -185, 400, -230, 60, -9/)
           a(2,:)=(/-610, -290, 1770, -1340, 560, -90/)
           a(3,:)=(/116, 840, -1620, 970, -360, 54/)
           a(4,:)=(/-22, -430, 710, -360, 120, -18/)
           a(5,:)=(/12, 65, -120, 60, -20, 3/)
           val=sum(B2(1:5)*matmul(a(1:5,:),F)/h)/shty
           return
        case (-1)
           val=sum(B2(1:4)*matmul(a(1:4,2:5),F(1:4))/h)/sv2
           return
        case (0)
           val=sum(B2(2:5)*matmul(a(1:4,2:5),F(2:5))/h)/sv2
           return
        case (1)
           val=sum(B2(3:6)*matmul(a(1:4,2:5),F(3:6))/h)/sv2
           return
        case (2)
           a(1,:)=(/-3, 20, -60, 120, -65, -12/)
           a(2,:)=(/18, -120, 360, -710, 430, 22/)
           a(3,:)=(/-54, 360, -970, 1620, -840, -116/)
           a(4,:)=(/90, -560, 1340, -1770, 290, 610/)
           a(5,:)=(/9, -60, 230, -400, 185, 36/)
           val=sum(B2(2:6)*matmul(a(1:5,:),F)/h)/shty
           return
        case default
           write(*,*)"Error in Get_Flux, wrong tag=",tag
           write(*,*)"Choices are -2, -1, 0, 1 or 2 for left boundary, "
           write(*,*)"interior or right boundary points"
           return
        end select
      case(7)
        a(0,:)=(/-137,300,-300,200,-75,12/)
        a(1,:)=(/-12,-65,120,-60,20,-3/)
        a(2,:)=(/3,-30,-20,60,-15,2/)
        a(3,:)=(/-2,15,-60,20,30,-3/)
        a(4,:)=(/3,-20,60,-120,65,12/)
        a(5,:)=(/-12,75,-200,300,-300,137/)

        cm=(/-1,7,7,-1/)
        cl=(/2,17,-11,5,-1/)
        cr=(/-1,5,-11,17,2/)

        select case(tag)
        case (-2)
           val=sum(cl*B2(1:5)*matmul(a(0:4,:),F)/sixty/h)/twelve
           return
        case (-1)
           val=sum(cm*B2(1:4)*matmul(a(0:3,:),F)/sixty/h)/twelve
           return
        case (0)
           val=sum(cm*B2(2:5)*matmul(a(1:4,:),F)/sixty/h)/twelve
           return
        case (1)
           val=sum(cm*B2(3:6)*matmul(a(2:5,:),F)/sixty/h)/twelve
           return
        case (2)
           val=sum(cr*B2(2:6)*matmul(a(1:5,:),F)/sixty/h)/twelve
           return
        case default
           write(*,*)"Error in Get_Flux, wrong tag=",tag
           write(*,*)"Choices are -2, -1, 0, 1 or 2 for left boundary, "
           write(*,*)"interior or right boundary points"
           return
        end select
      case default
        call pstop('Get_Flux in par_flux_4th','Stencil not available')
      end select

      end function Get_Flux

      end subroutine par_flux_4th

c     minmod
c     #################################################################
      function minmod(a,b)

      real(8) :: a,b,minmod

      if (a*b <= 0d0) then
        minmod = 0d0
        return
      endif

      if (a > 0d0) then
        minmod = min(a,b)
      else
        minmod = max(a,b)
      endif

      end function minmod

c     mc
c     #################################################################
      function mc(a,b,chk)

      real(8) :: a,b,mc,ho_flx

      logical,optional :: chk

      ho_flx = 0.5*(a+b)

      mc = minmod(2*minmod(a,b),ho_flx)

      chk = (mc /= ho_flx)

      end function mc

      end module grid_pardiff_ops

c grid_diff_ops
c ######################################################################
      module grid_diff_ops

       use grid_vertex_ops

       use grid_nabla_ops

       use grid_nabla2_ops

       use grid_tensor_ops

       use grid_advec_ops

       use grid_pardiff_ops

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
              !Averages in logical volume don't make sense
cc              dvol = gmetric%grid(igx)%dvol(i,j,k)
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


c grid_vertex_ops
c ######################################################################
      module grid_vertex_ops

        use grid_create

        INTERFACE curlcurl_vrtx
          module procedure curlcurl_vrtx_ijk,curlcurl_vrtx_mesh
        end INTERFACE

        INTERFACE gradDiv_vrtx
          module procedure gradDiv_vrtx_ijk,gradDiv_vrtx_mesh
        end INTERFACE

      contains

c     grad_vrtx
c     ##################################################################
      function grad_vrtx(g_def,i,j,k,igx,igy,igz,phi,vrtx) result(grd)

c     ------------------------------------------------------------------
c     Calculates gradient of phi from nearest neighbors.
c     If vrtx is true, it finds grad at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds grad at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igx,igy,igz
      real(8) :: phi(2,2,2),grd(3)
      logical :: vrtx

c     Local variables

      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0,nx,ny,nz
      real(8) :: idhx,idhy,idhz,fp,fm

c     Begin program

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      nx = g_def%nxgl(igx)
      ny = g_def%nygl(igy)
      nz = g_def%nzgl(igz)

      if (vrtx) then
        idhx = 1./g_def%dx(ig)
        idhy = 1./g_def%dy(jg)
        idhz = 1./g_def%dz(kg)
      else
        idhx = 1./g_def%dxh(ig)
        idhy = 1./g_def%dyh(jg)
        idhz = 1./g_def%dzh(kg)
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
      function curl_vrtx(g_def,i,j,k,igx,igy,igz,vv,vrtx) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curl of vv (covariant) from nearest neighbors. Returns
c     contravariant components.
c
c     If vrtx is true, it finds curl at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds curl at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igx,igy,igz
      real(8) :: vv(2,2,2,3),cnv(3)
      logical :: vrtx

c     Local variables

      real(8) :: gvx(3),gvy(3),gvz(3)

c     Begin program

      gvx = grad_vrtx(g_def,i,j,k,igx,igy,igz,vv(:,:,:,1),vrtx)
      gvy = grad_vrtx(g_def,i,j,k,igx,igy,igz,vv(:,:,:,2),vrtx)
      gvz = grad_vrtx(g_def,i,j,k,igx,igy,igz,vv(:,:,:,3),vrtx)

      cnv(1) = gvz(2)-gvy(3)
      cnv(2) = gvx(3)-gvz(1)
      cnv(3) = gvy(1)-gvx(2)

c     End program

      end function curl_vrtx

c     curl_cell2vrtx_mesh
c     ##################################################################
      function curl_cell2vrtx_mesh(g_def,igrid,vv) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curl of vv (covariant) at vertices from cell centers.
c     Returns contravariant components.
c
c     If vrtx is true, it finds curl at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds curl at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

      integer :: igrid
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
            cnv(i,j,k,:) = curl_vrtx(g_def,i,j,k,igrid,igrid,igrid
     .                              ,vv(i:i+1,j:j+1,k:k+1,:),.true.)
          enddo
        enddo
      enddo

c     End program

      end function curl_cell2vrtx_mesh

c     curl_vrtx2cell_mesh
c     ##################################################################
      function curl_vrtx2cell_mesh(g_def,igrid,vv) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curl of vv (covariant) from nearest neighbors. Returns
c     contravariant components.
c
c     If vrtx is true, it finds curl at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds curl at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: vv(0:,0:,0:,1:)
     $          ,cnv(0:size(vv,1)-1,0:size(vv,2)-1,0:size(vv,3)-1,3)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(vv,1)-2
      ny = size(vv,2)-2
      nz = size(vv,3)-2

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cnv(i,j,k,:) = curl_vrtx(g_def,i-1,j-1,k-1,igrid,igrid,igrid
     .                              ,vv(i-1:i,j-1:j,k-1:k,:),.false.)
          enddo
        enddo
      enddo

c     End program

      end function curl_vrtx2cell_mesh

c     div_vrtx
c     ##################################################################
      function div_vrtx(g_def,i,j,k,igx,igy,igz,v,vrtx) result(div)

c     ------------------------------------------------------------------
c     Calculates div of cnv vector defined at nearest neighbors.
c     If vrtx is true, it finds div at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds div at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igx,igy,igz
      real(8) :: v(2,2,2,3),div
      logical :: vrtx

c     Local variables

      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0
      real(8) :: dhx,dhy,dhz,fp,fm,jac

c     Begin program

      if (isSP(g_def,i+1,j,k,igx,igy,igz).and.vrtx) then
        div = 0d0
        return
      endif

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      if (vrtx) then
        dhx = g_def%dx(ig)
        dhy = g_def%dy(jg)
        dhz = g_def%dz(kg)
      else
        dhx = g_def%dxh(ig)
        dhy = g_def%dyh(jg)
        dhz = g_def%dzh(kg)
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
        jac = g_def%gmetric%grid(igx)%jac_v(i,j,k)
      else
        jac = g_def%gmetric%grid(igx)%jac(i,j,k)
      endif

      div = div/jac

c     End program

      end function div_vrtx

c     graddiv_vrtx_ijk
c     ###############################################################
      function graddiv_vrtx_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cnv,vol
     .                         ,ret_cnv) result(gdiv)

c     ---------------------------------------------------------------
c     Calculates grad(div(A)) at cell centers in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: cnv(0:nx+1,0:ny+1,0:nz+1,3),gdiv(3)

      logical,optional,intent(IN) :: vol,ret_cnv

c     Local variables

      integer :: il,jl,kl,ii,jj,kk

      real(8) :: div(2,2,2)

      logical :: vol_wgt,tocnv

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      tocnv = .false.
      if (PRESENT(ret_cnv)) tocnv = ret_cnv

c     Divergence at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2

            div(il,jl,kl) = div_vrtx(g_def,ii,jj,kk,igx,igy,igz
     .                         ,cnv(ii:ii+1,jj:jj+1,kk:kk+1,:)
     .                         ,.true.)
          enddo
        enddo
      enddo

c     Gradient at cell centers

      gdiv = grad_vrtx(g_def,i,j,k,igx,igy,igz,div,.false.)

      if (vol_wgt) gdiv=gdiv*g_def%gmetric%grid(igx)%dvol(i,j,k)

      if (tocnv) gdiv = XformToCnv(g_def,i,j,k,igx,gdiv)

c     End program

      end function graddiv_vrtx_ijk

c     graddiv_vrtx_mesh
c     ###############################################################
      function graddiv_vrtx_mesh(g_def,igr,cnv,vol,ret_cnv) result(gdiv)

c     ---------------------------------------------------------------
c     Calculates grad(div(A)) at cell centers in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A is
c     contravariant, and returns covariant (or contravariant if
c     ret_cnv=.true) components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr

      real(8) :: cnv(0:,0:,0:,:)

      real(8) :: gdiv(0:size(cnv,1)-1
     .               ,0:size(cnv,2)-1
     .               ,0:size(cnv,3)-1
     .               ,3)

      logical,optional,intent(IN) :: vol,ret_cnv

c     Local variables

      integer :: i,j,k,nx,ny,nz

c     Begin program

      nx = size(cnv,1)-2
      ny = size(cnv,2)-2
      nz = size(cnv,3)-2

      do k=1,nz
        do j=1,ny
          do i=1,nx
            gdiv(i,j,k,:) = gradDiv_vrtx_ijk(g_def,i,j,k,nx,ny,nz
     .                                      ,igr,igr,igr,cnv
     .                                      ,vol=vol,ret_cnv=ret_cnv)
          enddo
        enddo
      enddo

c     End program

      end function gradDiv_vrtx_mesh

c     lap_vrtx
c     ###############################################################
      function lap_vrtx(g_def,i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
     .         result(lap)

c     ---------------------------------------------------------------
c     Calculates lap(A) at cell centers in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

            cov = grad_vrtx(g_def,ii,jj,kk,igx,igy,igz
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1)
     .                     ,.true.)

            !Transform to contravariant at vertex
            cnv(il,jl,kl,:) =
     .          matmul(g_def%gmetric%grid(igx)%gsup_v(ii,jj,kk,:,:),cov)
          enddo
        enddo
      enddo

c     Divergence at cell centers

      lap = div_vrtx(g_def,i,j,k,igx,igy,igz,cnv,.false.)

      if (vol_wgt) lap=lap*g_def%gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_vrtx

c     curlcurl_vrtx_ijk
c     ###############################################################
      function curlcurl_vrtx_ijk(g_def,i,j,k,igr,aa,diff,vol)
     .         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(D*curl(A)) at cell centers in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igr

      real(8) :: aa(0:,0:,0:,1:),cnv(3)

      real(8),optional,target :: diff(0:size(aa,1)-1
     .                               ,0:size(aa,2)-1
     .                               ,0:size(aa,3)-1)

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
            cnv = curl_vrtx(g_def,ii,jj,kk,igr,igr,igr
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1,:)
     .                     ,.true.)

            cov = XformToCov_ijk_vrtx(g_def,ii,jj,kk,igr,cnv)

            if (have_coef) dff=0.125*sum(coef(ii:ii+1,jj:jj+1,kk:kk+1))

            crl(il,jl,kl,:) = cov*dff
          enddo
        enddo
      enddo

      nullify(coef)

c     Curl at cell center

      cnv = curl_vrtx(g_def,i,j,k,igr,igr,igr,crl,.false.)

      if (vol_wgt) cnv=cnv*g_def%gmetric%grid(igr)%dvol(i,j,k)

c     End program

      end function curlcurl_vrtx_ijk

c     curlcurl_vrtx_mesh
c     ##################################################################
      function curlcurl_vrtx_mesh(g_def,igrid,vv,diff) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curlcurl of vv (covariant) from nearest neighbors. 
c     Returns contravariant components at cell centers.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      logical :: vrtx
      real(8) :: vv(0:,0:,0:,1:)
     $          ,cnv(0:size(vv,1)-1,0:size(vv,2)-1,0:size(vv,3)-1,3)

      real(8),optional,target :: diff(0:size(vv,1)-1
     .                               ,0:size(vv,2)-1
     .                               ,0:size(vv,3)-1)

c     Local variables

      integer :: nx,ny,nz,ivar

      logical :: have_coef
      
      real(8),pointer,dimension(:,:,:) :: coef => null()

c     Begin program

      nx = size(vv,1)-2
      ny = size(vv,2)-2
      nz = size(vv,3)-2

      have_coef = PRESENT(diff)
      if (have_coef) coef => diff

      cnv = curl_cell2vrtx_mesh(g_def,igrid,vv)
      cnv = XformToCov_vrtx    (g_def,igrid,cnv)

      if (have_coef) then
        do ivar=1,3
          cnv(0:nx,0:ny,0:nz,ivar)
     .      = 0.125*(coef(0:nx  ,0:ny  ,0:nz  )
     .              +coef(0:nx  ,1:ny+1,0:nz  )
     .              +coef(0:nx  ,0:ny  ,1:nz+1)
     .              +coef(0:nx  ,1:ny+1,1:nz+1)
     .              +coef(1:nx+1,0:ny  ,0:nz  )
     .              +coef(1:nx+1,1:ny+1,0:nz  )
     .              +coef(1:nx+1,0:ny  ,1:nz+1)
     .              +coef(1:nx+1,1:ny+1,1:nz+1))
     .         *cnv(0:nx,0:ny,0:nz,ivar)
        enddo
      endif
      
      cnv = curl_vrtx2cell_mesh(g_def,igrid,cnv)

c     End program

      end function curlcurl_vrtx_mesh

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
      function grad_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,arr,hex,hey,hez
     .                 ,order) result(cov)

c     ---------------------------------------------------------------
c     Calculates grad(A) in general non-orthogonal coordinates.
c     The vector grad(A) is covariant.
c     ---------------------------------------------------------------

      use oned_int
      
      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1),cov(3)

      integer,optional :: hex,hey,hez,order

c     Local variables

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg,hx,hy,hz,ord,ordmax
     .          ,imin,imax,igmin,igmax
      real(8) :: dhx,dhy,dhz,dx1,dx2,dy1,dy2,dz1,dz2

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

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

c     Central-difference formulas

cc      dhx = 2.*g_def%dxh(ig)
cc      dhy = 2.*g_def%dyh(jg)
cc      dhz = 2.*g_def%dzh(kg)
cccc      dhx = g_def%dx(ig-1)+g_def%dx(ig)
cccc      dhy = g_def%dy(jg-1)+g_def%dy(jg)
cccc      dhz = g_def%dz(kg-1)+g_def%dz(kg)
cc
ccc     Exceptions
cc
cc      if (hx == 1) then
cc        im = i
cc        dhx = g_def%dx(ig)
cc      elseif (hx == -1) then
cc        ip = i
cc        dhx = g_def%dx(ig-1)
cc      endif
cc
cc      if (hy == 1) then
cc        jm = j
cc        dhy = g_def%dy(jg)
cc      elseif (hy == -1) then
cc        jp = j
cc        dhy = g_def%dy(jg-1)
cc      endif
cc
cc      if (hz == 1) then
cc        km = k
cc        dhz = g_def%dz(kg)
cc      elseif (hz == -1) then
cc        kp = k
cc        dhz = g_def%dz(kg-1)
cc      endif
cc
cc      gx = (arr(ip,j,k)-arr(im,j,k))/dhx
cc      gy = (arr(i,jp,k)-arr(i,jm,k))/dhy
cc      gz = (arr(i,j,kp)-arr(i,j,km))/dhz
cc
cc      cov = (/ gx,gy,gz /)
cc
cc      return

c     Arbitrary order formulas

      cov = 0d0

      if (PRESENT(order)) then
        ordmax = order
      else
        ordmax = 3
      endif
 
      !X comp
      ord = min(ordmax,nx+1)
      select case(hx)
      case(0)
        if (i==0) then
          call IntDriver1d(ord+1,g_def%xx(ig:ig+ord)
     .                    ,arr(i:i+ord,j,k)
     .                    ,1,g_def%xx(ig:ig),cov(1:1),ord,deriv=1)
        elseif (i==nx+1) then
          call IntDriver1d(ord+1,g_def%xx(ig-ord:ig),arr(i-ord:i,j,k)
     .                    ,1,g_def%xx(ig:ig),cov(1:1),ord,deriv=1)
        else
          dx1=g_def%dx(ig-1)
          dx2=g_def%dx(ig)
          cov(1) = (arr(i+1,j,k)-arr(i-1,j,k))/(dx1+dx2)
c$$$          cov(1) = q_int(3,arr(i-1:i+1,j,k),g_def%xx(ig-1:ig+1),g_def%xx(ig),1)
        endif
        
c$$$cc        if (i<(ord+1)/2+1) then
c$$$        if (i<2) then
c$$$          ord   = min(3,ord)
c$$$          imin  = max(i-ord/2-mod(ord,2),0)
c$$$     .          - max(i+ord/2-(nx+1),0) 
c$$$          imax  = min(i+ord/2 ,nx+1)
c$$$     .          - min(i-ord/2-mod(ord,2),0)
c$$$          igmin = ig - (i-imin)
c$$$          igmax = ig + (imax-i)
c$$$          call IntDriver1d(ord+1
c$$$     .                    ,g_def%xx(igmin:igmax),arr(imin:imax,j,k)
c$$$     .                    ,1,g_def%xx(ig:ig),cov(1:1),ord,deriv=1)
c$$$        elseif (i>=nx+1-(ord+1)/2) then
c$$$          ord   = min(3,ord)
c$$$          imin  = max(i-ord/2-mod(ord,2),0)
c$$$     .          - max(i+ord/2-(nx+1),0) 
c$$$          imax  = min(i+ord/2 ,nx+1)
c$$$     .          - min(i-ord/2-mod(ord,2),0)
c$$$          igmin = ig - (i-imin)
c$$$          igmax = ig + (imax-i)
c$$$          call IntDriver1d(ord+1
c$$$     .                    ,g_def%xx(igmin:igmax),arr(imin:imax,j,k)
c$$$     .                    ,1,g_def%xx(ig:ig),cov(1:1),ord,deriv=1)
c$$$        else
c$$$c$$$          imin  = max(i-ord/2-mod(ord,2),0)
c$$$c$$$     .          - max(i+ord/2-(nx+1),0) 
c$$$c$$$          imax  = min(i+ord/2 ,nx+1)
c$$$c$$$     .          - min(i-ord/2-mod(ord,2),0)
c$$$c$$$c$$$          imin  = max(i-ord/2,0)
c$$$c$$$c$$$     .          - max(i+ord/2+mod(ord,2)-(nx+1),0) 
c$$$c$$$c$$$          imax  = min(i+ord/2+mod(ord,2) ,nx+1)
c$$$c$$$c$$$     .          - min(i-ord/2,0)
c$$$c$$$          
c$$$c$$$          igmin = ig - (i-imin)
c$$$c$$$          igmax = ig + (imax-i)
c$$$c$$$
c$$$c$$$          write (*,*) "grad tst",i,imin,imax
c$$$c$$$
c$$$c$$$          call IntDriver1d(ord+1
c$$$c$$$     .                    ,g_def%xx(igmin:igmax),arr(imin:imax,j,k)
c$$$c$$$     .                    ,1,g_def%xx(ig:ig),cov(1:1),ord,deriv=1)
c$$$          write (*,*) i,ord
c$$$          ord = 4
c$$$          call IntDriver1d(ord+1,g_def%xx(ig-(ord+1)/2:ig+ord/2)
c$$$     .                    ,arr(i-(ord+1)/2:i+ord/2,j,k)
c$$$     .                    ,1,g_def%xx(ig:ig),cov(1:1),ord,deriv=1)
c$$$        endif
      case(1)
        if (i<nx+1) cov(1) =(arr(i+1,j,k)-arr(i  ,j,k))/g_def%dx(ig)
      case(-1)
        if (i>0   ) cov(1) =(arr(i  ,j,k)-arr(i-1,j,k))/g_def%dx(ig-1)
      end select

      !Y comp
      ord = min(ordmax,ny+1)
      select case (hy)
      case(0)
        if (j==0) then
          call IntDriver1d (ord+1,g_def%yy(jg:jg+ord),arr(i,j:j+ord,k)
     .                     ,1,g_def%yy(jg:jg),cov(2:2),ord,deriv=1)
        elseif (j==ny+1) then
          call IntDriver1d (ord+1,g_def%yy(jg-ord:jg),arr(i,j-ord:j,k)
     .                     ,1,g_def%yy(jg:jg),cov(2:2),ord,deriv=1)
        else
          dy1=g_def%dy(jg-1)
          dy2=g_def%dy(jg)
          cov(2) = (arr(i,j+1,k)-arr(i,j-1,k))/(dy1+dy2)
c$$$          cov(2) = q_int(3,arr(i,j-1:j+1,k),g_def%yy(jg-1:jg+1),g_def%yy(jg),1)
        endif
      case(1)
        if (j<ny+1) cov(2) =(arr(i,j+1,k)-arr(i,j  ,k))/g_def%dy(jg  )
      case(-1)
        if (j>0   ) cov(2) =(arr(i,j  ,k)-arr(i,j-1,k))/g_def%dy(jg-1)
      end select

      !Z comp
      ord = min(ordmax,nz+1)
      select case(hz)
      case(0)
        if (k==0) then
          call IntDriver1d (ord+1,g_def%zz(kg:kg+ord),arr(i,j,k:k+ord)
     .                     ,1,g_def%zz(kg:kg),cov(3:3),ord,deriv=1)
        elseif (k==nz+1) then
          call IntDriver1d (ord+1,g_def%zz(kg-ord:kg),arr(i,j,k-ord:k)
     .                     ,1,g_def%zz(kg:kg),cov(3:3),ord,deriv=1)
        else
          dz1=g_def%dz(kg-1)
          dz2=g_def%dz(kg)
          cov(3) = (arr(i,j,k+1)-arr(i,j,k-1))/(dz1+dz2)
c$$$          cov(3) = q_int(3,arr(i,j,k-1:k+1),g_def%zz(kg-1:kg+1),g_def%zz(kg),1)
        endif
      case(1)
        if (k<nz+1) cov(3) =(arr(i,j,k+1)-arr(i,j,k  ))/g_def%dz(kg  )
      case(-1)
        if (k>0   ) cov(3) =(arr(i,j,k  )-arr(i,j,k-1))/g_def%dz(kg-1)
      end select

c     End program

      end function grad_ijk

c     grad_mesh
c     ###############################################################
      function grad_mesh(g_def,igrid,arr,hex,hey,hez,order) result(cov)

c     ---------------------------------------------------------------
c     Calculates grad(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: arr(0:,0:,0:)
     .          ,cov(0:size(arr,1)-1,0:size(arr,2)-1,0:size(arr,3)-1,3)

      integer,optional :: hex,hey,hez,order

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            cov(i,j,k,:)=grad_ijk(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                         ,arr,hex=hex,hey=hey,hez=hez,order=order)
          enddo
        enddo
      enddo

c     End program

      end function grad_mesh

c     curl_ijk
c     ###############################################################
      function curl_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,a,hex,hey,hez
     .                 ,order)
     .         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: i,j,k,comp,nx,ny,nz,igx,igy,igz

      real(8)    :: cnv(3)

      real(8)    :: a(0:nx+1,0:ny+1,0:nz+1,3)

      integer,optional :: hex,hey,hez,order

c     Local variables

      integer :: icomp
      real(8) :: gv(3,3)

c     Begin program

      do icomp=1,3
        gv(:,icomp) = grad_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,a(:,:,:,icomp),hex=hex,hey=hey,hez=hez
     .                    ,order=order)
      enddo

      cnv(1) = gv(2,3)-gv(3,2)
      cnv(2) = gv(3,1)-gv(1,3)
      cnv(3) = gv(1,2)-gv(2,1)

c     End program

      end function curl_ijk

c     curl_mesh
c     ###############################################################
      function curl_mesh(g_def,igrid,cov,hex,hey,hez,order) result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cnv(0:size(cov,1)-1,0:size(cov,2)-1,0:size(cov,3)-1,3)

      integer,optional :: hex,hey,hez,order

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            cnv(i,j,k,:)=curl_ijk(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                         ,cov,hex=hex,hey=hey,hez=hez,order=order)
          enddo
        enddo
      enddo

c     End program

      end function curl_mesh

c     div_ijk
c     ###############################################################
      function div_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cnv
     .                ,hex,hey,hez,sp) result(div)

      implicit none

c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      dxx = 2*g_def%dxh(ig)
      dyy = 2*g_def%dyh(jg)
      dzz = 2*g_def%dzh(kg)

      if (hx == 1) then
        im = i
        dxx = g_def%dx(ig)
      elseif (hx == -1) then
        ip = i
        dxx = g_def%dx(ig-1)
      endif

      if (hy == 1) then
        jm = j
        dyy = g_def%dy(jg)
      elseif (hy == -1) then
        jp = j
        dyy = g_def%dy(jg-1)
      endif

      if (hz == 1) then
        km = k
        dzz = g_def%dz(kg)
      elseif (hz == -1) then
        kp = k
        dzz = g_def%dz(kg-1)
      endif

      if (isSP2(g_def,i,igx).and.hx == 0.and.spoint) then
        jacp = g_def%gmetric%grid(igrid)%jac(i+1,j,k)
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
cccc        dxx = g_def%dx(ig)
cc        jac0 = 0.5*(g_def%gmetric%grid(igrid)%jac(ip,j,k)
cc     .             +g_def%gmetric%grid(igrid)%jac(i ,j,k))
cc      case(2)
cccc        jm = j
cccc        dyy = g_def%dy(jg)
cc        jac0 = 0.5*(g_def%gmetric%grid(igrid)%jac(i,jp,k)
cc     .             +g_def%gmetric%grid(igrid)%jac(i,j ,k))
cc      case(3)
cccc        km = k
cccc        dzz = g_def%dz(kg)
cc        jac0 = 0.5*(g_def%gmetric%grid(igrid)%jac(i,j,kp)
cc     .             +g_def%gmetric%grid(igrid)%jac(i,j,k ))
cc      case default
        jac0 = g_def%gmetric%grid(igrid)%jac(i,j,k)
cc      end select

      div = div/jac0
      
c     End 

      end function div_ijk

c     div_mesh
c     ###############################################################
      function div_mesh(g_def,igrid,cnv) result(scl)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cnv(0:,0:,0:,:)
     .          ,scl(0:size(cnv,1)-1,0:size(cnv,2)-1,0:size(cnv,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      scl = 0d0

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            scl(i,j,k) = div_ijk(g_def,i,j,k,nx,ny,nz
     .                          ,igrid,igrid,igrid,cnv)
          enddo
        enddo
      enddo

c     End program

      end function div_mesh

c     div_stg_ijk
c     ###############################################################
      function div_stg_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cnv)
     .         result(div)

      implicit none

c     ---------------------------------------------------------------
c     Calculates divergence of a staggered vector field at cell centers
c     in general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cnv(0:nx+1,0:ny+1,0:nz+1,3),div

c     Local variables

      integer :: ig,jg,kg,igrid,ip,im,jp,jm,kp,km,hx,hy,hz
      real(8) :: idx,idy,idz

c     Begin program

      igrid = igx

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      im = i-1
      jm = j-1
      km = k-1

      idx = 1d0/g_def%dxh(ig)
      idy = 1d0/g_def%dyh(jg)
      idz = 1d0/g_def%dzh(kg)
      if (nx == 1) idx = 0d0
      if (ny == 1) idy = 0d0
      if (nz == 1) idz = 0d0

      div =  (cnv(i,j,k,1)-cnv(im,j ,k ,1))*idx
     .      +(cnv(i,j,k,2)-cnv(i ,jm,k ,2))*idy
     .      +(cnv(i,j,k,3)-cnv(i ,j ,km,3))*idz

      div = div/g_def%gmetric%grid(igrid)%jac(i,j,k)
      
c     End 

      end function div_stg_ijk

c     div_stg_mesh
c     ###############################################################
      function div_stg_mesh(g_def,igrid,curv,iscnv,bc,bc_order)
     .         result(scl)

c     ---------------------------------------------------------------
c     Calculates div of a staggered vector field in general
c     non-orthogonal coordinates. The vector A can be covariant or
c     contravariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: curv(0:,0:,0:,:)
     .          ,scl(0:size(curv,1)-1,0:size(curv,2)-1,0:size(curv,3)-1)

      logical,optional :: iscnv,bc
      integer,optional :: bc_order

c     Local variables

      integer :: bcs(6),nx,ny,nz,i,j,k
      real(8) :: cnv(0:size(curv,1)-1
     .              ,0:size(curv,2)-1
     .              ,0:size(curv,3)-1
     .              ,size(curv,4))

      real(8) :: zero(0:size(curv,1)-1
     .               ,0:size(curv,2)-1
     .               ,0:size(curv,3)-1)

      logical :: is_cnv,do_bc
      integer :: bc_ord

c     Begin program

      if (PRESENT(iscnv)) then
        is_cnv = iscnv
      else
        is_cnv = .true.
      endif

      if (PRESENT(bc)) then
        do_bc = bc
      else
        do_bc = PRESENT(bc_order)
      endif

      if (PRESENT(bc_order)) then
        bc_ord = bc_order
      else
        bc_ord = 1
      endif

c     Compute divergence

      scl = 0d0

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      if (is_cnv) then

        do k=1,nz
          do j=1,ny
            do i=1,nx
              scl(i,j,k)=div_stg_ijk(g_def,i,j,k,nx,ny,nz
     .                              ,igrid,igrid,igrid,curv)
            enddo
          enddo
        enddo

      else

        cnv = XformToCnv_stg_mesh(g_def,igrid,curv)

        do k=1,nz
          do j=1,ny
            do i=1,nx
              scl(i,j,k)=div_stg_ijk(g_def,i,j,k,nx,ny,nz
     .                              ,igrid,igrid,igrid,cnv)
            enddo
          enddo
        enddo

      endif

c     Apply topological BCs to divergence

      if (do_bc) then
        call get_topol_bcs(bcs)
        where(bcs == DEF) bcs = DIR
        zero = 0d0
        call setBC(g_def,1,nx,ny,nz,scl,zero,bcs
     .            ,igrid,igrid,igrid,iorder=bc_ord)
      endif

c     End program

      end function div_stg_mesh

c     curl_stg_mesh
c     ###############################################################
      function curl_stg_mesh(g_def,igr,cov,vrtx) result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates for a
c     staggered vector A.  The vector A is covariant, and returns a
c     contravariant vector, also staggered on the dual mesh.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr

      real(8) :: cov(0:,0:,0:,:)
     .          ,cnv(0:size(cov,1)-1
     .              ,0:size(cov,2)-1
     .              ,0:size(cov,3)-1
     .                ,size(cov,4))

      logical,optional :: vrtx

c     Local variables

      integer :: nx,ny,nz,i,j,k
      logical :: vertex

c     Begin program

      if (PRESENT(vrtx)) then
        vertex = vrtx
      else
        vertex = .false.
      endif

      cnv = 0d0

      nx = size(cov,1)-2
      ny = size(cov,2)-2
      nz = size(cov,3)-2

      if (vertex) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              cnv(i,j,k,:)=
     .           curl_stg_vrtx_ijk(g_def,i,j,k,nx,ny,nz,igr,igr,igr,cov)
            enddo
          enddo
        enddo
      else
        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              cnv(i,j,k,:)=
     .           curl_stg_ijk(g_def,i,j,k,nx,ny,nz,igr,igr,igr,cov)
            enddo
          enddo
        enddo
      endif

c     End program

      end function curl_stg_mesh

c     curl_stg_ijk
c     ###############################################################
      function curl_stg_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cov)
     .         result(curl)

      implicit none

c     ---------------------------------------------------------------
c     Calculates curl of CELL staggered covariant vector field to find
c     VERTEX staggered contravariant field in general non-orthogonal
c     geometry.
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cov(0:nx+1,0:ny+1,0:nz+1,3),curl(3)

c     Local variables

      integer :: ig,jg,kg,igrid,ip,jp,kp,km,hx,hy,hz
      real(8) :: dxx,dyy,dzz

c     Begin program

      igrid = igx

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      jp = j+1
      kp = k+1

      dxx = g_def%dx(ig)
      dyy = g_def%dy(jg)
      dzz = g_def%dz(kg)

      curl = 0d0

      if (jp <= ny+1.and.kp <= nz+1) then
        curl(1) = (cov(i,jp,k ,3)-cov(i,j,k,3))/dyy
     .           -(cov(i,j ,kp,2)-cov(i,j,k,2))/dzz !@(i,j+1/2,k+1/2)
      endif
      if (ip <= nx+1.and.kp <= nz+1) then
        curl(2) = (cov(i ,j,kp,1)-cov(i,j,k,1))/dzz
     .           -(cov(ip,j,k ,3)-cov(i,j,k,3))/dxx !@(i+1/2,j,k+1/2)
      endif
      if (ip <= nx+1.and.jp <= ny+1) then
        curl(3) = (cov(ip,j ,k,2)-cov(i,j,k,2))/dxx
     .           -(cov(i ,jp,k,1)-cov(i,j,k,1))/dyy !@(i+1/2,j+1/2,k)
      endif
      
c     End 

      end function curl_stg_ijk

c     curl_stg_vrtx_ijk
c     ###############################################################
      function curl_stg_vrtx_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cov)
     .         result(curl)

      implicit none

c     ---------------------------------------------------------------
c     Calculates curl of VERTEX staggered covariant vector field to find
c     CELL staggered contravariant field in general non-orthogonal
c     geometry.
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cov(0:nx+1,0:ny+1,0:nz+1,3),curl(3)

c     Local variables

      integer :: ig,jg,kg,igrid,im,jm,km
      real(8) :: dxx,dyy,dzz

c     Begin program

      igrid = igx

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      im = max(i-1,0)
      jm = max(j-1,0)
      km = max(k-1,0)

      dxx = g_def%dx(ig)
      dyy = g_def%dy(jg)
      dzz = g_def%dz(kg)

      curl(1) = (cov(i,j,k,3)-cov(i,jm,k ,3))/dyy
     .         -(cov(i,j,k,2)-cov(i,j ,km,2))/dzz !@(i+1/2,j,k)
      curl(2) = (cov(i,j,k,1)-cov(i ,j,km,1))/dzz
     .         -(cov(i,j,k,3)-cov(im,j,k ,3))/dxx !@(i,j+1/2,k)
      curl(3) = (cov(i,j,k,2)-cov(im,j ,k,2))/dxx
     .         -(cov(i,j,k,1)-cov(i ,jm,k,1))/dyy !@(i,j,k+1/2)
      
c     End 

      end function curl_stg_vrtx_ijk

      end module grid_nabla_ops

c grid_nabla2_ops
c ######################################################################
      module grid_nabla2_ops

        use grid_create

        use grid_nabla_ops

        INTERFACE lap
          module procedure lap_mesh,lap_ijk
        end INTERFACE

        INTERFACE gradDiv
          module procedure gradDiv_mesh,gradDiv_ijk
        end INTERFACE

      contains

c     gradDiv_ijk
c     ##################################################################
      function gradDiv_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,da,vol,cnv)
     .         result(vec)

c     -----------------------------------------------------------------
c     Calculates grad(div(A)) in general non-orthogonal coordinates in a
c     compact form. The vector A is contravariant, and returns covariant
c     components unless cnv=true.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      jacip  = g_def%gmetric%grid(igx)%jac(ip,j,k)
      jacim  = g_def%gmetric%grid(igx)%jac(im,j,k)
      jacjp  = g_def%gmetric%grid(igx)%jac(i,jp,k)
      jacjm  = g_def%gmetric%grid(igx)%jac(i,jm,k)
      jackp  = g_def%gmetric%grid(igx)%jac(i,j,kp)
      jackm  = g_def%gmetric%grid(igx)%jac(i,j,km)
      jac    = g_def%gmetric%grid(igx)%jac(i,j,k )

c     Fluxes at faces for calculation of grad[div(da)]

      !!Divergence at faces i+-1/2, etc.
      divip = (da(ip,j ,k,1)-da(i ,j ,k,1))/g_def%dx(ig)
     .       +(da(i ,jp,k,2)-da(i ,jm,k,2)
     .        +da(ip,jp,k,2)-da(ip,jm,k,2))/g_def%dyh(jg)*0.25
     .       +(da(i ,j,kp,3)-da(i ,j,km,3)
     .        +da(ip,j,kp,3)-da(ip,j,km,3))/g_def%dzh(kg)*0.25
      divim = (da(i ,j ,k,1)-da(im,j ,k,1))/g_def%dx(ig-1)
     .       +(da(i ,jp,k,2)-da(i ,jm,k,2)
     .        +da(im,jp,k,2)-da(im,jm,k,2))/g_def%dyh(jg)*0.25
     .       +(da(i ,j,kp,3)-da(i ,j,km,3)
     .        +da(im,j,kp,3)-da(im,j,km,3))/g_def%dzh(kg)*0.25

      divjp = (da(ip,j ,k,1)-da(im,j ,k,1)
     .        +da(ip,jp,k,1)-da(im,jp,k,1))/g_def%dxh(ig)*0.25
     .       +(da(i ,jp,k,2)-da(i ,j ,k,2))/g_def%dy(jg)
     .       +(da(i,j ,kp,3)-da(i,j ,km,3)
     .        +da(i,jp,kp,3)-da(i,jp,km,3))/g_def%dzh(kg)*0.25
      divjm = (da(ip,j ,k,1)-da(im,j ,k,1)
     .        +da(ip,jm,k,1)-da(im,jm,k,1))/g_def%dxh(ig)*0.25
     .       +(da(i ,j ,k,2)-da(i ,jm,k,2))/g_def%dy(jg-1)
     .       +(da(i,j ,kp,3)-da(i,j ,km,3)
     .        +da(i,jm,kp,3)-da(i,jm,km,3))/g_def%dzh(kg)*0.25

      divkp = (da(ip,j,k ,1)-da(im,j,k ,1)
     .        +da(ip,j,kp,1)-da(im,j,kp,1))/g_def%dxh(ig)*0.25
     .       +(da(i,jp,k ,2)-da(i,jm,k ,2)
     .        +da(i,jp,kp,2)-da(i,jm,kp,2))/g_def%dyh(jg)*0.25
     .       +(da(i,j ,kp,3)-da(i,j ,k ,3))/g_def%dz(kg)
      divkm = (da(ip,j,k ,1)-da(im,j,k ,1)
     .        +da(ip,j,km,1)-da(im,j,km,1))/g_def%dxh(ig)*0.25
     .       +(da(i,jp,k ,2)-da(i,jm,k ,2)
     .        +da(i,jp,km,2)-da(i,jm,km,2))/g_def%dyh(jg)*0.25
     .       +(da(i,j ,k ,3)-da(i,j ,km,3))/g_def%dz(kg-1)

      divip = 2*divip/(jac+jacip)
      if (.not.isSP(g_def,i,j,k,igx,igy,igz)) then
        divim = 2*divim/(jac+jacim)
      else
        divim = 0d0
      endif

      divjp = 2*divjp/(jac+jacjp)
      divjm = 2*divjm/(jac+jacjm)

      divkp = 2*divkp/(jac+jackp)
      divkm = 2*divkm/(jac+jackm)

c     Calculate gradient

      vec(1) = (divip - divim)/g_def%dxh(ig)
      vec(2) = (divjp - divjm)/g_def%dyh(jg)
      vec(3) = (divkp - divkm)/g_def%dzh(kg)

      if (vol_wgt) vec=vec*g_def%gmetric%grid(igx)%dvol(i,j,k)

      if (cnvt) vec = XformToCnv(g_def,i,j,k,igx,vec)

      end function gradDiv_ijk

c     gradDiv_mesh
c     ###############################################################
      function gradDiv_mesh(g_def,igr,cnv,vol,ret_cnv) result(cov)
c     ---------------------------------------------------------------
c     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr

      real(8) :: cnv(0:,0:,0:,:)
     $          ,cov(0:size(cnv,1)-1
     $              ,0:size(cnv,2)-1
     $              ,0:size(cnv,3)-1,3)

      logical,optional,intent(IN) :: vol,ret_cnv

c     Local variables

      integer :: i,j,k,nx,ny,nz

c     Begin program

      nx = size(cnv,1)-2
      ny = size(cnv,2)-2
      nz = size(cnv,3)-2

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cov(i,j,k,:) = gradDiv_ijk(g_def,i,j,k,nx,ny,nz,igr,igr,igr
     .                                ,cnv,vol=vol,cnv=ret_cnv)
          enddo
        enddo
      enddo

c     End program

      end function gradDiv_mesh

c     lap_flx_x_ijk
c     ###############################################################
      function lap_flx_x_ijk(g_def,i,j,k,igrid,arr,diff,vertex)
     .         result(flx)
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

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (vrtx_met) then
        !Vertex metrics
        gsup => g_def%gmetric%grid(igrid)%gsup_v

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
        gsup => g_def%gmetric%grid(igrid)%gsup

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

      flx   = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/g_def%dx(ig)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
     .                  -arr(ip,j ,k)-arr(i,j ,k))/g_def%dy(jg)
     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
     .                  -arr(ip,jm,k)-arr(i,jm,k))/g_def%dy(jg-1)
     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
     .                  -arr(ip,j,k )-arr(i,j,k ))/g_def%dz(kg)
     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
     .                  -arr(ip,j,km)-arr(i,j,km))/g_def%dz(kg-1))

c     End program

      end function lap_flx_x_ijk

c     lap_flx_y_ijk
c     ###############################################################
      function lap_flx_y_ijk(g_def,i,j,k,igrid,arr,diff,vertex)
     .         result(flx)
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

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (vrtx_met) then
        !Vertex metrics
        gsup => g_def%gmetric%grid(igrid)%gsup_v

        d_xy_ipjp = 0.5*(gsup(i ,j,k ,1,2)
     .                  +gsup(i ,j,km,1,2))
        d_xy_imjp = 0.5*(gsup(im,j,k ,1,2)
     .                  +gsup(im,j,km,1,2))

        d_yz_jpkp = 0.5*(gsup(i ,j,k ,2,3)
     .                  +gsup(im,j,k ,2,3))
        d_yz_jpkm = 0.5*(gsup(i ,j,km,2,3)
     .                  +gsup(im,j,km,2,3))

        !To fix singularity with gsup(2,2)
        if (isSP2(g_def,i,igrid)) then
          gsup => g_def%gmetric%grid(igrid)%gsup
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
        gsup => g_def%gmetric%grid(igrid)%gsup

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
                                                                       
      flx =   d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/g_def%dy(jg)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
     .                  -arr(i ,jp,k)-arr(i ,j,k))/g_def%dx(ig)
     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
     .                  -arr(im,jp,k)-arr(im,j,k))/g_def%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
     .                  -arr(i,jp,k )-arr(i,j,k ))/g_def%dz(kg)
     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
     .                  -arr(i,jp,km)-arr(i,j,km))/g_def%dz(kg-1))

c     End program

      end function lap_flx_y_ijk

c     lap_flx_z_ijk
c     ###############################################################
      function lap_flx_z_ijk(g_def,i,j,k,igrid,arr,diff,vertex)
     .         result(flx)
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

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (vrtx_met) then
        !Vertex metrics
        gsup => g_def%gmetric%grid(igrid)%gsup_v

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
        gsup => g_def%gmetric%grid(igrid)%gsup

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
                                                                       
      flx   = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/g_def%dz(kg)
     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
     .                  -arr(i ,j,kp)-arr(i ,j,k))/g_def%dx(ig)
     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
     .                  -arr(im,j,kp)-arr(im,j,k))/g_def%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
     .                  -arr(i,j ,kp)-arr(i,j ,k))/g_def%dy(jg)
     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
     .                  -arr(i,jm,kp)-arr(i,jm,k))/g_def%dy(jg-1))

c     End program

      end function lap_flx_z_ijk

c     lap_ijk
c     ###############################################################
      function lap_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol
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

      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      logical :: vol_wgt,vrtx_met

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      im = i-1
      jm = j-1
      km = k-1

      flxip=lap_flx_x_ijk(g_def,i ,j,k,igx,arr,diff=dff
     .                   ,vertex=vertex_metrics)
      flxim=lap_flx_x_ijk(g_def,im,j,k,igx,arr,diff=dff
     .                   ,vertex=vertex_metrics)

      flxjp=lap_flx_y_ijk(g_def,i,j ,k,igy,arr,diff=dff
     .                   ,vertex=vertex_metrics)
      flxjm=lap_flx_y_ijk(g_def,i,jm,k,igy,arr,diff=dff
     .                   ,vertex=vertex_metrics)

      flxkp=lap_flx_z_ijk(g_def,i,j,k ,igz,arr,diff=dff
     .                   ,vertex=vertex_metrics)
      flxkm=lap_flx_z_ijk(g_def,i,j,km,igz,arr,diff=dff
     .                   ,vertex=vertex_metrics)

      laplacian = ((flxip -flxim)/g_def%dxh(ig)
     .            +(flxjp -flxjm)/g_def%dyh(jg)
     .            +(flxkp -flxkm)/g_def%dzh(kg))
     .            /g_def%gmetric%grid(igx)%jac(i,j,k)

      if (vol_wgt) laplacian=laplacian
     .                      *g_def%gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_ijk

c     lap_mesh
c     ###############################################################
      function lap_mesh(g_def,igr,arr,dff,vol,vertex_metrics)result(lap)
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

      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer :: i,j,k,nx,ny,nz

c     Begin program

      nx = size(arr,1)-2
      ny = size(arr,2)-2
      nz = size(arr,3)-2

      lap = 0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx
            lap(i,j,k) = lap_ijk(g_def,i,j,k,nx,ny,nz,igr,igr,igr,arr
     $              ,dff=dff,vol=vol,vertex_metrics=vertex_metrics)
          enddo
        enddo
      enddo

c     End program

      end function lap_mesh

c     lap_tnsr_flx_x_ijk
c     ###############################################################
      function lap_tnsr_flx_x_ijk(g_def,i,j,k,igrid,arr,tnsr)result(flx)
c     ---------------------------------------------------------------
c     Calculates X-flux for scalar laplacian at faces (i+1/2,j,k)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx
      real(8),target :: tnsr(0:,0:,0:,:)

c     Local variables

      integer :: im,jm,km,ip,jp,kp,ig,jg,kg

      real(8) :: d_xx_ip
     .          ,d_xy_ipjp,d_xy_ipjm
     .          ,d_xz_ipkp,d_xz_ipkm

c     Begin program

      call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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

      flx   = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/g_def%dx(ig)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
     .                  -arr(ip,j ,k)-arr(i,j ,k))/g_def%dy(jg)
     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
     .                  -arr(ip,jm,k)-arr(i,jm,k))/g_def%dy(jg-1)
     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
     .                  -arr(ip,j,k )-arr(i,j,k ))/g_def%dz(kg)
     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
     .                  -arr(ip,j,km)-arr(i,j,km))/g_def%dz(kg-1))

c     End program

      end function lap_tnsr_flx_x_ijk

c     lap_tnsr_flx_y_ijk
c     ###############################################################
      function lap_tnsr_flx_y_ijk(g_def,i,j,k,igrid,arr,tnsr)result(flx)
c     ---------------------------------------------------------------
c     Calculates Y-flux for scalar laplacian at faces (i,j+1/2,k)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx
      real(8),target :: tnsr(0:,0:,0:,:)

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_yy_jp
     .             ,d_xy_ipjp,d_xy_imjp
     .             ,d_yz_jpkp,d_yz_jpkm

c     Begin program

      call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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
                                                                       
      flx =   d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/g_def%dy(jg)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
     .                  -arr(i ,jp,k)-arr(i ,j,k))/g_def%dx(ig)
     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
     .                  -arr(im,jp,k)-arr(im,j,k))/g_def%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
     .                  -arr(i,jp,k )-arr(i,j,k ))/g_def%dz(kg)
     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
     .                  -arr(i,jp,km)-arr(i,j,km))/g_def%dz(kg-1))

c     End program

      end function lap_tnsr_flx_y_ijk

c     lap_tnsr_flx_z_ijk
c     ###############################################################
      function lap_tnsr_flx_z_ijk(g_def,i,j,k,igrid,arr,tnsr)result(flx)
c     ---------------------------------------------------------------
c     Calculates Z-flux for scalar laplacian at faces (i,j,k+1/2)
c     in general non-orthogonal coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igrid

      real(8) :: arr (0:,0:,0:),flx
      real(8),target :: tnsr(0:,0:,0:,:)

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_zz_kp
     .             ,d_xz_ipkp,d_xz_imkp
     .             ,d_yz_jpkp,d_yz_jmkp

c     Begin program

      call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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
                                                                       
      flx   = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/g_def%dz(kg)
     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
     .                  -arr(i ,j,kp)-arr(i ,j,k))/g_def%dx(ig)
     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
     .                  -arr(im,j,kp)-arr(im,j,k))/g_def%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
     .                  -arr(i,j ,kp)-arr(i,j ,k))/g_def%dy(jg)
     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
     .                  -arr(i,jm,kp)-arr(i,jm,k))/g_def%dy(jg-1))

c     End program

      end function lap_tnsr_flx_z_ijk

c     lap_tnsr
c     ###############################################################
      function lap_tnsr(g_def,i,j,k,nx,ny,nz,igx,igy,igz,tnsr,arr,vol)
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

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      im = i-1
      jm = j-1
      km = k-1

      flxip=lap_tnsr_flx_x_ijk(g_def,i ,j,k,igx,arr,tnsr)
      flxim=lap_tnsr_flx_x_ijk(g_def,im,j,k,igx,arr,tnsr)

      flxjp=lap_tnsr_flx_y_ijk(g_def,i,j ,k,igx,arr,tnsr)
      flxjm=lap_tnsr_flx_y_ijk(g_def,i,jm,k,igx,arr,tnsr)

      flxkp=lap_tnsr_flx_z_ijk(g_def,i,j,k ,igx,arr,tnsr)
      flxkm=lap_tnsr_flx_z_ijk(g_def,i,j,km,igx,arr,tnsr)

      laplacian = ((flxip -flxim)/g_def%dxh(ig)
     .            +(flxjp -flxjm)/g_def%dyh(jg)
     .            +(flxkp -flxkm)/g_def%dzh(kg))
     .            /g_def%gmetric%grid(igx)%jac(i,j,k)

      if (vol_wgt)
     .     laplacian=laplacian*g_def%gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_tnsr

      end module grid_nabla2_ops

c grid_tensor_ops
c ######################################################################
      module grid_tensor_ops

        use grid_create

        use grid_nabla2_ops

        use grid_vertex_ops

        INTERFACE veclap
          module procedure veclap_mesh,veclap_ijk
        end INTERFACE
        
        INTERFACE veclap_cnv
          module procedure veclap_cnv_mesh,veclap_cnv_ijk
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

c$$$c     fnabla_v
c$$$c     #############################################################
c$$$      function fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
c$$$     .                 ,half_elem,cnv)
c$$$     .         result(tensor)
c$$$c     -------------------------------------------------------------
c$$$c     Calculates the tensor components T_l^m (if cnv = .false.) or T^lm
c$$$c     (if cnv = .true.) of nabla(vec v), with v contravariant, and fills
c$$$c     tensor(l,m) at the following positions:
c$$$c       + half_elem=0 --> i,j,k
c$$$c       + half_elem=1 --> i+1/2,j,k
c$$$c       + half_elem=2 --> i,j+1/2,k
c$$$c       + half_elem=3 --> i,j,k+1/2
c$$$c     -------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
c$$$        real(8)    :: tensor(3,3)
c$$$        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
c$$$     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
c$$$     .               ,az(0:nx+1,0:ny+1,0:nz+1)
c$$$
c$$$        logical,optional :: cnv
c$$$
c$$$c     Local variables
c$$$
c$$$        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
c$$$        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
c$$$        real(8)    :: vxx,vyy,vzz
c$$$     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
c$$$     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
c$$$     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
c$$$        real(8)    :: tsrc(3,3),gsuper(3,3),nblv(3,3)
c$$$     .               ,jacip,jacjp,jacjm,jackp,jackm
c$$$     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac
c$$$
c$$$        logical    :: cntv
c$$$
c$$$c     Begin program
c$$$
c$$$        igrid = igx
c$$$
c$$$        !Defaults
c$$$
c$$$        if (PRESENT(cnv)) then
c$$$          cntv = cnv
c$$$        else
c$$$          cntv = .false.
c$$$        endif
c$$$
c$$$        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$        dhx = 2.*g_def%dxh(ig)
c$$$        dhy = 2.*g_def%dyh(jg)
c$$$        dhz = 2.*g_def%dzh(kg)
c$$$
c$$$        ip = i+1
c$$$        im = i-1
c$$$        jp = j+1
c$$$        jm = j-1
c$$$        kp = k+1
c$$$        km = k-1
c$$$
c$$$        nblv = nabla_v_src(i,j,k)
c$$$
c$$$        !Exceptions
c$$$
c$$$        select case(half_elem)
c$$$        case (1)
c$$$          dhx = g_def%dx(ig)
c$$$
c$$$          if (cntv) 
c$$$c$$$            gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:)
c$$$c$$$     .                   +g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:))
c$$$     .      gsuper=.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k-1,:,:))
c$$$
c$$$          vxip = ax(ip,j,k)
c$$$          vxim = ax(i ,j,k)
c$$$          vyip = ay(ip,j,k)
c$$$          vyim = ay(i ,j,k)
c$$$          vzip = az(ip,j,k)
c$$$          vzim = az(i ,j,k)
c$$$
c$$$          vxjp = 0.5*(ax(i,jp,k)+ax(ip,jp,k))
c$$$          vxjm = 0.5*(ax(i,jm,k)+ax(ip,jm,k))
c$$$          vyjp = 0.5*(ay(i,jp,k)+ay(ip,jp,k))
c$$$          vyjm = 0.5*(ay(i,jm,k)+ay(ip,jm,k))
c$$$          vzjp = 0.5*(az(i,jp,k)+az(ip,jp,k))
c$$$          vzjm = 0.5*(az(i,jm,k)+az(ip,jm,k))
c$$$                                      
c$$$          vxkp = 0.5*(ax(i,j,kp)+ax(ip,j,kp))
c$$$          vxkm = 0.5*(ax(i,j,km)+ax(ip,j,km))
c$$$          vykp = 0.5*(ay(i,j,kp)+ay(ip,j,kp))
c$$$          vykm = 0.5*(ay(i,j,km)+ay(ip,j,km))
c$$$          vzkp = 0.5*(az(i,j,kp)+az(ip,j,kp))
c$$$          vzkm = 0.5*(az(i,j,km)+az(ip,j,km))
c$$$
c$$$          tsrc = 0.5*(nabla_v_src(ip,j,k)+nblv)
c$$$
c$$$        case (2)
c$$$          dhy = g_def%dy(jg)
c$$$
c$$$          if (cntv)
c$$$c$$$     .         gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:)
c$$$c$$$     .                      +g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:))
c$$$     .      gsuper=.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k-1,:,:))
c$$$
c$$$          vxip = 0.5*(ax(ip,j,k)+ax(ip,jp,k))
c$$$          vxim = 0.5*(ax(im,j,k)+ax(im,jp,k))
c$$$          vyip = 0.5*(ay(ip,j,k)+ay(ip,jp,k))
c$$$          vyim = 0.5*(ay(im,j,k)+ay(im,jp,k))
c$$$          vzip = 0.5*(az(ip,j,k)+az(ip,jp,k))
c$$$          vzim = 0.5*(az(im,j,k)+az(im,jp,k))
c$$$
c$$$          vxjp = ax(i,jp,k)
c$$$          vxjm = ax(i,j ,k)
c$$$          vyjp = ay(i,jp,k)
c$$$          vyjm = ay(i,j ,k)
c$$$          vzjp = az(i,jp,k)
c$$$          vzjm = az(i,j ,k)
c$$$
c$$$          vxkp = 0.5*(ax(i,j,kp)+ax(i,jp,kp))
c$$$          vxkm = 0.5*(ax(i,j,km)+ax(i,jp,km))
c$$$          vykp = 0.5*(ay(i,j,kp)+ay(i,jp,kp))
c$$$          vykm = 0.5*(ay(i,j,km)+ay(i,jp,km))
c$$$          vzkp = 0.5*(az(i,j,kp)+az(i,jp,kp))
c$$$          vzkm = 0.5*(az(i,j,km)+az(i,jp,km))
c$$$
c$$$          tsrc = 0.5*(nabla_v_src(i,jp,k)+nblv)
c$$$
c$$$        case (3)
c$$$          dhz = g_def%dz(kg)
c$$$
c$$$          if (cntv)
c$$$c$$$     .         gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:)
c$$$c$$$     .                      +g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:))
c$$$     .      gsuper=.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j-1,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j-1,k,:,:))
c$$$
c$$$          vxip = 0.5*(ax(ip,j,k)+ax(ip,j,kp))
c$$$          vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
c$$$          vyip = 0.5*(ay(ip,j,k)+ay(ip,j,kp))
c$$$          vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
c$$$          vzip = 0.5*(az(ip,j,k)+az(ip,j,kp))
c$$$          vzim = 0.5*(az(im,j,k)+az(im,j,kp))
c$$$
c$$$          vxjp = 0.5*(ax(i,jp,k)+ax(i,jp,kp))
c$$$          vxjm = 0.5*(ax(i,jm,k)+ax(i,jm,kp))
c$$$          vyjp = 0.5*(ay(i,jp,k)+ay(i,jp,kp))
c$$$          vyjm = 0.5*(ay(i,jm,k)+ay(i,jm,kp))
c$$$          vzjp = 0.5*(az(i,jp,k)+az(i,jp,kp))
c$$$          vzjm = 0.5*(az(i,jm,k)+az(i,jm,kp))
c$$$
c$$$          vxkp = ax(i,j,kp)
c$$$          vxkm = ax(i,j,k )
c$$$          vykp = ay(i,j,kp)
c$$$          vykm = ay(i,j,k )
c$$$          vzkp = az(i,j,kp)
c$$$          vzkm = az(i,j,k )
c$$$
c$$$          tsrc = 0.5*(nabla_v_src(i,j,kp)+nblv)
c$$$
c$$$        case default
c$$$
c$$$          if (cntv) gsuper=g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:)
c$$$
c$$$          vxip = ax(ip,j,k)
c$$$          vxim = ax(im,j,k)
c$$$          vyip = ay(ip,j,k)
c$$$          vyim = ay(im,j,k)
c$$$          vzip = az(ip,j,k)
c$$$          vzim = az(im,j,k)
c$$$
c$$$          vxjp = ax(i,jp,k)
c$$$          vxjm = ax(i,jm,k)
c$$$          vyjp = ay(i,jp,k)
c$$$          vyjm = ay(i,jm,k)
c$$$          vzjp = az(i,jp,k)
c$$$          vzjm = az(i,jm,k)
c$$$
c$$$          vxkp = ax(i,j,kp)
c$$$          vxkm = ax(i,j,km)
c$$$          vykp = ay(i,j,kp)
c$$$          vykm = ay(i,j,km)
c$$$          vzkp = az(i,j,kp)
c$$$          vzkm = az(i,j,km)
c$$$
c$$$          tsrc = nblv
c$$$
c$$$        end select
c$$$
c$$$        idhx = 1d0/dhx
c$$$        idhy = 1d0/dhy
c$$$        idhz = 1d0/dhz
c$$$
c$$$      ! l = 1, m = 1
c$$$        tensor(1,1) = (vxip-vxim)*idhx
c$$$
c$$$      ! l = 1, m = 2
c$$$        tensor(1,2) = (vyip-vyim)*idhx
c$$$
c$$$      ! l = 1, m = 3
c$$$        tensor(1,3) = (vzip-vzim)*idhx
c$$$
c$$$      ! l = 2, m = 1
c$$$        tensor(2,1) = (vxjp-vxjm)*idhy
c$$$
c$$$      ! l = 2, m = 2
c$$$        tensor(2,2) = (vyjp-vyjm)*idhy
c$$$
c$$$      ! l = 2, m = 3
c$$$        tensor(2,3) = (vzjp-vzjm)*idhy
c$$$
c$$$      ! l = 3, m = 1
c$$$        tensor(3,1) = (vxkp-vxkm)*idhz
c$$$
c$$$      ! l = 3, m = 2
c$$$        tensor(3,2) = (vykp-vykm)*idhz
c$$$
c$$$      ! l = 3, m = 3
c$$$        tensor(3,3) = (vzkp-vzkm)*idhz
c$$$
c$$$      ! Add geometric source
c$$$
c$$$        tensor = tensor + tsrc
c$$$
c$$$        if (cntv) tensor = matmul(gsuper,tensor)
c$$$
c$$$c     End program
c$$$
c$$$      contains
c$$$
c$$$c     nabla_v_src
c$$$c     #############################################################
c$$$      function nabla_v_src(i,j,k) result(tensor)
c$$$
c$$$c     -------------------------------------------------------------
c$$$c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c$$$c     -------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        integer    :: i,j,k
c$$$        real(8)    :: tensor(3,3)
c$$$
c$$$c     Local variables
c$$$
c$$$        integer    :: l,m
c$$$        real(8)    :: cnv(3),dum(3)
c$$$
c$$$c     Begin program
c$$$
c$$$        cnv = (/ ax(i,j,k),ay(i,j,k),az(i,j,k) /)
c$$$
c$$$        dum = g_def%gmetric%grid(igrid)%Gamma(i,j,k,1,1,:)
c$$$     .       +g_def%gmetric%grid(igrid)%Gamma(i,j,k,2,2,:)
c$$$     .       +g_def%gmetric%grid(igrid)%Gamma(i,j,k,3,3,:)
c$$$
c$$$        do l=1,3
c$$$          do m=1,3
c$$$            tensor(l,m) =
c$$$     .              dot_product(cnv
c$$$     .                   ,g_def%gmetric%grid(igrid)%Gamma(i,j,k,m,:,l))
c$$$     .             -cnv(m)*dum(l)
c$$$cc     .           cnv(1)* g_def%gmetric%grid(igrid)%Gamma(i,j,k,m,1,l)
c$$$cc     .          +cnv(2)* g_def%gmetric%grid(igrid)%Gamma(i,j,k,m,2,l)
c$$$cc     .          +cnv(3)* g_def%gmetric%grid(igrid)%Gamma(i,j,k,m,3,l)
c$$$cc     .          -cnv(m)*(g_def%gmetric%grid(igrid)%Gamma(i,j,k,1,1,l)
c$$$cc     .                  +g_def%gmetric%grid(igrid)%Gamma(i,j,k,2,2,l)
c$$$cc     .                  +g_def%gmetric%grid(igrid)%Gamma(i,j,k,3,3,l))
c$$$          enddo
c$$$        enddo
c$$$
c$$$c     End program
c$$$
c$$$      end function nabla_v_src
c$$$
c$$$      end function fnabla_v

c     fnabla_v
c     #############################################################
      function fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz,aa
     .                 ,half_elem,cnv) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_l^m of nabla(vec v), with
c     v covariant, and fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
        real(8) :: tensor(3,3)
        real(8) :: aa(0:nx+1,0:ny+1,0:nz+1,3)

        logical,optional :: cnv

c     Local variables

        integer :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
        real(8) :: dhx,dhy,dhz,idh(3)
        real(8) :: v(3),vip(3),vim(3),vjp(3),vjm(3),vkp(3),vkm(3)
        real(8) :: tsrc(3,3),gsuper(3,3)
        logical :: cntv
        
c     Begin program

        igrid = igx

        if (PRESENT(cnv)) then
          cntv = cnv
        else
          cntv = .false.
        endif

        !Defaults

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*g_def%dxh(ig)
        dhy = 2.*g_def%dyh(jg)
        dhz = 2.*g_def%dzh(kg)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        select case(half_elem)
        case (1)
          dhx = g_def%dx(ig)

          vip(:) = aa(ip,j,k,:)
          vim(:) = aa(i ,j,k,:)

          vjp(:) = 0.5*(aa(ip,jp,k,:)+aa(i,jp,k,:))
          vjm(:) = 0.5*(aa(ip,jm,k,:)+aa(i,jm,k,:))

          vkp(:) = 0.5*(aa(ip,j,kp,:)+aa(i,j,kp,:))
          vkm(:) = 0.5*(aa(ip,j,km,:)+aa(i,j,km,:))

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

          if (cntv)
c$$$     .    gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k-1,:,:))
     .      gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .                   +g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        case (2)
          dhy = g_def%dy(jg)

          vip(:) = 0.5*(aa(ip,j,k,:)+aa(ip,jp,k,:))
          vim(:) = 0.5*(aa(im,j,k,:)+aa(im,jp,k,:))

          vjp(:) = aa(i,jp,k,:)
          vjm(:) = aa(i,j ,k,:)
                
          vkp(:) = 0.5*(aa(i,j,kp,:)+aa(i,jp,kp,:))
          vkm(:) = 0.5*(aa(i,j,km,:)+aa(i,jp,km,:))

          tsrc = 0.5*(nabla_v_src(i,jp,k)
     .               +nabla_v_src(i,j ,k))

          if (cntv)
c$$$     .    gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k-1,:,:))
     .    gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        case (3)
          dhz = g_def%dz(kg)

          vip(:) = 0.5*(aa(ip,j,k,:)+aa(ip,j,kp,:))
          vim(:) = 0.5*(aa(im,j,k,:)+aa(im,j,kp,:))

          vjp(:) = 0.5*(aa(i,jp,k,:)+aa(i,jp,kp,:))
          vjm(:) = 0.5*(aa(i,jm,k,:)+aa(i,jm,kp,:))
                
          vkp(:) = aa(i,j,kp,:)
          vkm(:) = aa(i,j,k ,:)

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

          if (cntv)
c$$$     .    gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j-1,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j-1,k,:,:))
     .    gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        case default

          vip(:) = aa(ip,j,k,:)
          vim(:) = aa(im,j,k,:)

          vjp(:) = aa(i,jp,k,:)
          vjm(:) = aa(i,jm,k,:)
                
          vkp(:) = aa(i,j,kp,:)
          vkm(:) = aa(i,j,km,:)

          tsrc = nabla_v_src(i,j,k)

          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

        end select

        idh(1) = 1d0/dhx
        idh(2) = 1d0/dhy
        idh(3) = 1d0/dhz

        tensor(1,:) = (vip(:)-vim(:))*idh(1)
        tensor(2,:) = (vjp(:)-vjm(:))*idh(2)
        tensor(3,:) = (vkp(:)-vkm(:))*idh(3)

c$$$        if (cntv) tensor = matmul(gsuper,tensor)

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

        cnv = aa(i,j,k,:)

        dum = g_def%gmetric%grid(igrid)%Gamma(i,j,k,1,1,:)
     .       +g_def%gmetric%grid(igrid)%Gamma(i,j,k,2,2,:)
     .       +g_def%gmetric%grid(igrid)%Gamma(i,j,k,3,3,:)

        do l=1,3
          do m=1,3
            tensor(l,m) =
     .              dot_product(cnv
     .                   ,g_def%gmetric%grid(igrid)%Gamma(i,j,k,m,:,l))
     .             -cnv(m)*dum(l)
          enddo
        enddo

c$$$        if (cntv)
c$$$     .   tensor=matmul(g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:),tensor)

c     End program

      end function nabla_v_src
      
      end function fnabla_v

c     fnabla_v_upwd
c     #############################################################
      function fnabla_v_upwd(g_def,i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
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

        type(grid_mg_def),pointer :: g_def

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

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*g_def%dxh(ig)
        dhy = 2.*g_def%dyh(jg)
        dhz = 2.*g_def%dzh(kg)

        hessian = -g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

c     Exceptions

        if (hex == 1) then
          im = i
          dhx = g_def%dx(ig)
        elseif (hex == -1) then
          ip = i
          dhx = g_def%dx(ig-1)
        endif

        if (hey == 1) then
          jm = j
          dhy = g_def%dy(jg)
        elseif (hey == -1) then
          jp = j
          dhy = g_def%dy(jg-1)
        endif

        if (hez == 1) then
          km = k
          dhz = g_def%dz(kg)
        elseif (hez == -1) then
          kp = k
          dhz = g_def%dz(kg-1)
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
      
c     fnabla_v_cnv
c     #################################################################
      function fnabla_v_cnv(g_def,i,j,k,igrid,vcnv,half_elem)
     .         result(tensor)
c     -----------------------------------------------------------------
c     Calculates the tensor components T^lm of nabla(vec v), with
c     v contravariant, and fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,half_elem,igrid
        real(8) :: tensor(3,3)
        real(8) :: vcnv(0:,0:,0:,:)

c     Local variables

        integer :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8) :: dhx,dhy,dhz,idh(3)
        real(8) :: v(3),vip(3),vim(3),vjp(3),vjm(3),vkp(3),vkm(3)
        real(8) :: tsrc(3,3),gsuper(3,3),jach

        real(8),pointer,dimension(:,:,:) :: jac

c     Begin program

        !Defaults

        call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

        dhx = 2.*g_def%dxh(ig)
        dhy = 2.*g_def%dyh(jg)
        dhz = 2.*g_def%dzh(kg)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        jac => g_def%gmetric%grid(igrid)%jac

        select case(half_elem)
        case (1)
          dhx = g_def%dx(ig)

          vip(:) = vcnv(ip,j,k,:)/jac(ip,j,k)
          vim(:) = vcnv(i ,j,k,:)/jac(i ,j,k)

          vjp(:) = 0.5*(vcnv(ip,jp,k,:)/jac(ip,jp,k)
     .                 +vcnv(i ,jp,k,:)/jac(i ,jp,k))
          vjm(:) = 0.5*(vcnv(ip,jm,k,:)/jac(ip,jm,k)
     .                 +vcnv(i ,jm,k,:)/jac(i ,jm,k))

          vkp(:) = 0.5*(vcnv(ip,j,kp,:)/jac(ip,j,kp)
     .                 +vcnv(i ,j,kp,:)/jac(i ,j,kp))
          vkm(:) = 0.5*(vcnv(ip,j,km,:)/jac(ip,j,km)
     .                 +vcnv(i ,j,km,:)/jac(i ,j,km))

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k-1,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:))
c$$$          jach   =0.25*(g_def%gmetric%grid(igrid)%jac_v(i,j  ,k  )
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i,j-1,k  )
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i,j  ,k-1)
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i,j-1,k-1))
          jach  = 0.5*(g_def%gmetric%grid(igrid)%jac (ip,j,k)
     .                +g_def%gmetric%grid(igrid)%jac (i ,j,k))

        case (2)
          dhy = g_def%dy(jg)

          vip(:) = 0.5*(vcnv(ip,j ,k,:)/jac(ip,j ,k)
     .                 +vcnv(ip,jp,k,:)/jac(ip,jp,k))
          vim(:) = 0.5*(vcnv(im,j ,k,:)/jac(im,j ,k)
     .                 +vcnv(im,jp,k,:)/jac(im,jp,k))

          vjp(:) = vcnv(i,jp,k,:)/jac(i,jp,k)
          vjm(:) = vcnv(i,j ,k,:)/jac(i,j ,k)
                
          vkp(:) = 0.5*(vcnv(i,j ,kp,:)/jac(i,j ,kp)
     .                 +vcnv(i,jp,kp,:)/jac(i,jp,kp))
          vkm(:) = 0.5*(vcnv(i,j ,km,:)/jac(i,j ,km)
     .                 +vcnv(i,jp,km,:)/jac(i,jp,km))

          tsrc = 0.5*(nabla_v_src(i,jp,k)
     .               +nabla_v_src(i,j ,k))

c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k-1,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:))
c$$$          jach   =0.25*(g_def%gmetric%grid(igrid)%jac_v(i  ,j,k  )
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i-1,j,k  )
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i  ,j,k-1)
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i-1,j,k-1))
          jach  = 0.5*(g_def%gmetric%grid(igrid)%jac (i,jp,k)
     .                +g_def%gmetric%grid(igrid)%jac (i,j ,k))

        case (3)
          dhz = g_def%dz(kg)

          vip(:) = 0.5*(vcnv(ip,j,k ,:)/jac(ip,j,k )
     .                 +vcnv(ip,j,kp,:)/jac(ip,j,kp))
          vim(:) = 0.5*(vcnv(im,j,k ,:)/jac(im,j,k )
     .                 +vcnv(im,j,kp,:)/jac(im,j,kp))

          vjp(:) = 0.5*(vcnv(i,jp,k ,:)/jac(i,jp,k )
     .                 +vcnv(i,jp,kp,:)/jac(i,jp,kp))
          vjm(:) = 0.5*(vcnv(i,jm,k ,:)/jac(i,jm,k )
     .                 +vcnv(i,jm,kp,:)/jac(i,jm,kp))
                
          vkp(:) = vcnv(i,j,kp,:)/jac(i,j,kp)
          vkm(:) = vcnv(i,j,k ,:)/jac(i,j,k )

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j-1,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j-1,k,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:))
c$$$          jach   =0.25*(g_def%gmetric%grid(igrid)%jac_v(i  ,j  ,k)
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i-1,j  ,k)
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i  ,j-1,k)
c$$$     .                 +g_def%gmetric%grid(igrid)%jac_v(i-1,j-1,k))
          jach  = 0.5*(g_def%gmetric%grid(igrid)%jac (i,j,kp)
     .                +g_def%gmetric%grid(igrid)%jac (i,j,k ))

        case default

          vip(:) = vcnv(ip,j,k,:)/jac(ip,j,k)
          vim(:) = vcnv(im,j,k,:)/jac(im,j,k)
                                             
          vjp(:) = vcnv(i,jp,k,:)/jac(i,jp,k)
          vjm(:) = vcnv(i,jm,k,:)/jac(i,jm,k)
                                             
          vkp(:) = vcnv(i,j,kp,:)/jac(i,j,kp)
          vkm(:) = vcnv(i,j,km,:)/jac(i,j,km)

          tsrc = nabla_v_src(i,j,k)

          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

          jach   = g_def%gmetric%grid(igrid)%jac (i,j,k)

        end select

        !Compute gradients
        
        idh(1) = 1d0/dhx
        idh(2) = 1d0/dhy
        idh(3) = 1d0/dhz

        if (jach == 0d0) jach = SP_flsv
        
        tensor(1,:) = (vip(:)-vim(:))*idh(1)*jach
        tensor(2,:) = (vjp(:)-vjm(:))*idh(2)*jach
        tensor(3,:) = (vkp(:)-vkm(:))*idh(3)*jach

cc        tensor = matmul(gsuper,tensor)

        !Add geometric source

        tensor = tensor + tsrc

        tensor = matmul(gsuper,tensor)

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

        hessian = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

        tensor(:,1) = matmul(vcnv(i,j,k,:),hessian(1,:,:))
        tensor(:,2) = matmul(vcnv(i,j,k,:),hessian(2,:,:))
        tensor(:,3) = matmul(vcnv(i,j,k,:),hessian(3,:,:))

cc        tensor=matmul(g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:),tensor)
        
c     End program

      end function nabla_v_src

      end function fnabla_v_cnv

c     fnabla_v_mixed_cov
c     #############################################################
      function fnabla_v_mixed_cov(g_def,i,j,k,igrid,vcov,half_elem)
     .         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T^l_m of nabla(vec v), with
c     v covariant, and fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,half_elem,nx,ny,nz,igrid
        real(8) :: tensor(3,3)
        real(8) :: vcov(0:,0:,0:,:)

c     Local variables

        integer :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8) :: dhx,dhy,dhz,idh(3)
        real(8) :: v(3),vip(3),vim(3),vjp(3),vjm(3),vkp(3),vkm(3)
        real(8) :: tsrc(3,3),gsuper(3,3)

c     Begin program

        !Defaults

        call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

        dhx = 2.*g_def%dxh(ig)
        dhy = 2.*g_def%dyh(jg)
        dhz = 2.*g_def%dzh(kg)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        select case(half_elem)
        case (1)
          dhx = g_def%dx(ig)

          vip(:) = vcov(ip,j,k,:)
          vim(:) = vcov(i ,j,k,:)

          vjp(:) = 0.5*(vcov(ip,jp,k,:)+vcov(i,jp,k,:))
          vjm(:) = 0.5*(vcov(ip,jm,k,:)+vcov(i,jm,k,:))

          vkp(:) = 0.5*(vcov(ip,j,kp,:)+vcov(i,j,kp,:))
          vkm(:) = 0.5*(vcov(ip,j,km,:)+vcov(i,j,km,:))

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k-1,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        case (2)
          dhy = g_def%dy(jg)

          vip(:) = 0.5*(vcov(ip,j,k,:)+vcov(ip,jp,k,:))
          vim(:) = 0.5*(vcov(im,j,k,:)+vcov(im,jp,k,:))

          vjp(:) = vcov(i,jp,k,:)
          vjm(:) = vcov(i,j ,k,:)
                
          vkp(:) = 0.5*(vcov(i,j,kp,:)+vcov(i,jp,kp,:))
          vkm(:) = 0.5*(vcov(i,j,km,:)+vcov(i,jp,km,:))

          tsrc = 0.5*(nabla_v_src(i,jp,k)
     .               +nabla_v_src(i,j ,k))

c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k-1,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        case (3)
          dhz = g_def%dz(kg)

          vip(:) = 0.5*(vcov(ip,j,k,:)+vcov(ip,j,kp,:))
          vim(:) = 0.5*(vcov(im,j,k,:)+vcov(im,j,kp,:))

          vjp(:) = 0.5*(vcov(i,jp,k,:)+vcov(i,jp,kp,:))
          vjm(:) = 0.5*(vcov(i,jm,k,:)+vcov(i,jm,kp,:))
                
          vkp(:) = vcov(i,j,kp,:)
          vkm(:) = vcov(i,j,k ,:)

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j-1,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j-1,k,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        case default

          vip(:) = vcov(ip,j,k,:)
          vim(:) = vcov(im,j,k,:)

          vjp(:) = vcov(i,jp,k,:)
          vjm(:) = vcov(i,jm,k,:)
                
          vkp(:) = vcov(i,j,kp,:)
          vkm(:) = vcov(i,j,km,:)

          tsrc = nabla_v_src(i,j,k)

          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

        end select

        idh(1) = 1d0/dhx
        idh(2) = 1d0/dhy
        idh(3) = 1d0/dhz

        tensor(1,:) = (vip(:)-vim(:))*idh(1)
        tensor(2,:) = (vjp(:)-vjm(:))*idh(2)
        tensor(3,:) = (vkp(:)-vkm(:))*idh(3)

cc        tensor = matmul(gsuper,tensor)

      ! Add geometric source

cc        tsrc = nabla_v_src2(i,j,k,half_elem)

        tensor = tensor + tsrc

        tensor = matmul(gsuper,tensor)

c     End program

      contains

c     nabla_v_src
c     #############################################################
      function nabla_v_src(i,j,k) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

        integer :: i,j,k
        real(8) :: tensor(3,3)

        real(8) :: hessian(3,3,3),gsuper(3,3)

c     Begin program

        hessian = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

        tensor(:,1) = -matmul(vcov(i,j,k,:),hessian(:,:,1))
        tensor(:,2) = -matmul(vcov(i,j,k,:),hessian(:,:,2))
        tensor(:,3) = -matmul(vcov(i,j,k,:),hessian(:,:,3))

cc        gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)
cc        tensor = matmul(gsuper,tensor)

c     End program

      end function nabla_v_src

c     nabla_v_src2
c     #############################################################
      function nabla_v_src2(i,j,k,he) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

        integer :: i,j,k,he
        real(8) :: tensor(3,3)

        real(8) :: hessian(3,3,3),gsuper(3,3),cov(3)

c     Begin program

        select case(he)
        case (1)

          hessian = 0.5*(g_def%gmetric%grid(igrid)%Gamma(i ,j,k,:,:,:)
     .                  +g_def%gmetric%grid(igrid)%Gamma(ip,j,k,:,:,:))

          cov = 0.5*(vcov(i ,j,k,:)
     .              +vcov(ip,j,k,:))
          
        case (2)

          hessian = 0.5*(g_def%gmetric%grid(igrid)%Gamma(i,j ,k,:,:,:)
     .                  +g_def%gmetric%grid(igrid)%Gamma(i,jp,k,:,:,:))

          cov = 0.5*(vcov(i,j ,k,:)
     .              +vcov(i,jp,k,:))

        case (3)
          
          hessian = 0.5*(g_def%gmetric%grid(igrid)%Gamma(i,j,k ,:,:,:)
     .                  +g_def%gmetric%grid(igrid)%Gamma(i,j,kp,:,:,:))

          cov = 0.5*(vcov(i,j,k ,:)
     .              +vcov(i,j,kp,:))

        case default
          
          hessian = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

          cov = vcov(i,j,k,:)
          
        end select

        tensor(:,1) = -matmul(cov,hessian(:,:,1))
        tensor(:,2) = -matmul(cov,hessian(:,:,2))
        tensor(:,3) = -matmul(cov,hessian(:,:,3))

cc        gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)
cc        tensor = matmul(gsuper,tensor)

c     End program

      end function nabla_v_src2

      end function fnabla_v_mixed_cov

c     fnabla_v_cov
c     #############################################################
      function fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
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

        type(grid_mg_def),pointer :: g_def

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

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*g_def%dxh(ig)
        dhy = 2.*g_def%dyh(jg)
        dhz = 2.*g_def%dzh(kg)

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
          dhx = g_def%dx(ig)

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
          dhy = g_def%dy(jg)

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
          dhz = g_def%dz(kg)

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

        hessian = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

      ! l = 1, m = 1
        tensor(1,1) =- ax(i,j,k)* hessian(1,1,1)
     .               - ay(i,j,k)* hessian(2,1,1)
     .               - az(i,j,k)* hessian(3,1,1)

      ! l = 1, m = 2
        tensor(1,2) =- ax(i,j,k)* hessian(1,1,2)
     .               - ay(i,j,k)* hessian(2,1,2)
     .               - az(i,j,k)* hessian(3,1,2)

      ! l = 1, m = 3
        tensor(1,3) =- ax(i,j,k)* hessian(1,1,3)
     .               - ay(i,j,k)* hessian(2,1,3)
     .               - az(i,j,k)* hessian(3,1,3)

      ! l = 2, m = 1
        tensor(2,1) =- ax(i,j,k)* hessian(1,2,1)
     .               - ay(i,j,k)* hessian(2,2,1)
     .               - az(i,j,k)* hessian(3,2,1)

      ! l = 2, m = 2
        tensor(2,2) =- ax(i,j,k)* hessian(1,2,2)
     .               - ay(i,j,k)* hessian(2,2,2)
     .               - az(i,j,k)* hessian(3,2,2)

      ! l = 2, m = 3
        tensor(2,3) =- ax(i,j,k)* hessian(1,2,3)
     .               - ay(i,j,k)* hessian(2,2,3)
     .               - az(i,j,k)* hessian(3,2,3)

      ! l = 3, m = 1
        tensor(3,1) =- ax(i,j,k)* hessian(1,3,1)
     .               - ay(i,j,k)* hessian(2,3,1)
     .               - az(i,j,k)* hessian(3,3,1)

      ! l = 3, m = 2
        tensor(3,2) =- ax(i,j,k)* hessian(1,3,2)
     .               - ay(i,j,k)* hessian(2,3,2)
     .               - az(i,j,k)* hessian(3,3,2)

      ! l = 3, m = 3
        tensor(3,3) =- ax(i,j,k)* hessian(1,3,3)
     .               - ay(i,j,k)* hessian(2,3,3)
     .               - az(i,j,k)* hessian(3,3,3)

c     End program

      end function nabla_v_src

      end function fnabla_v_cov

c     fnabla_v_upwd_cov
c     #############################################################
      function fnabla_v_upwd_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,ax,ay,az,hex,hey,hez) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_lm of nabla(vec v), with
c     v covariant, and fills tensor(l,m) at the following positions:
c       + hex,hey,hez = 0 => i,j,k
c       + hex=+-1 --> i+-1/2
c       + hey=+-1 --> j+-1/2
c       + hez=+-1 --> k+-1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,hex,hey,hez,nx,ny,nz,igx,igy,igz
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
        real(8)    :: cf1,cf2,cf3,cf4,hessian(3,3,3)

        logical    :: spoint

c     Begin program

        igrid = igx

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*g_def%dxh(ig)
        dhy = 2.*g_def%dyh(jg)
        dhz = 2.*g_def%dzh(kg)

        hessian = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

c     Exceptions

        if (hex == 1) then
          im = i
          dhx = g_def%dx(ig)
        elseif (hex == -1) then
          ip = i
          dhx = g_def%dx(ig-1)
        endif

        if (hey == 1) then
          jm = j
          dhy = g_def%dy(jg)
        elseif (hey == -1) then
          jp = j
          dhy = g_def%dy(jg-1)
        endif

        if (hez == 1) then
          km = k
          dhz = g_def%dz(kg)
        elseif (hez == -1) then
          kp = k
          dhz = g_def%dz(kg-1)
        endif

c     Vectors

        vxx = ax(i,j,k)
        vyy = ay(i,j,k)
        vzz = az(i,j,k)

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

c     Calculate nabla_v tensor
        
        idhx = 1d0/dhx
        idhy = 1d0/dhy
        idhz = 1d0/dhz

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)*idhx
     .               -(vxx*hessian(1,1,1)
     .               + vyy*hessian(2,1,1)
     .               + vzz*hessian(3,1,1))

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)*idhx
     .               -(vxx*hessian(1,1,2)
     .               + vyy*hessian(2,1,2)
     .               + vzz*hessian(3,1,2))

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)*idhx
     .               -(vxx*hessian(1,1,3)
     .               + vyy*hessian(2,1,3)
     .               + vzz*hessian(3,1,3))

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)*idhy
     .               -(vxx*hessian(1,2,1)
     .               + vyy*hessian(2,2,1)
     .               + vzz*hessian(3,2,1))

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)*idhy
     .               -(vxx*hessian(1,2,2)
     .               + vyy*hessian(2,2,2)
     .               + vzz*hessian(3,2,2))

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)*idhy
     .               -(vxx*hessian(1,2,3)
     .               + vyy*hessian(2,2,3)
     .               + vzz*hessian(3,2,3))

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)*idhz
     .               -(vxx*hessian(1,3,1)
     .               + vyy*hessian(2,3,1)
     .               + vzz*hessian(3,3,1))

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)*idhz
     .               -(vxx*hessian(1,3,2)
     .               + vyy*hessian(2,3,2)
     .               + vzz*hessian(3,3,2))

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)*idhz
     .               -(vxx*hessian(1,3,3)
     .               + vyy*hessian(2,3,3)
     .               + vzz*hessian(3,3,3))

c     End program

      end function fnabla_v_upwd_cov

c     dtnsr_flx_ijk
c     ###############################################################
      function dtnsr_flx_ijk(g_def,i,j,k,nx,ny,nz,igrid,alt_eom,flux)
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

      type(grid_mg_def),pointer :: g_def

      real(8) :: flx(3)

      integer :: i,j,k,igrid,nx,ny,nz

      logical :: alt_eom

      INTERFACE
        subroutine flux(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flux
      END INTERFACE

c     Local variables

      real(8) :: jac

c     Begin program

      call flux(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid,alt_eom
     .         ,flx(1),flx(2),flx(3),1)

c     End program

      end function dtnsr_flx_ijk

c     dtnsr_src_ijk
c     ###############################################################
      function dtnsr_src_ijk(g_def,i,j,k,nx,ny,nz,igr,alt_eom
     .                      ,flxx,flxy,flxz)
     .         result(src)
c     ---------------------------------------------------------------
c     Calculates source for contravariant components of div(tensor) in
c     general non-orthogonal coordinates at (i,j,k). The routines
c     flxx,flxy,flxz must provide the contravariant tersor components
c     T^lm (if flag /= 0, they actually provide T^lm/jac).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8) :: src(3)

      integer :: i,j,k,igr,nx,ny,nz

      logical :: alt_eom

      INTERFACE
        subroutine flxx(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flxx

        subroutine flxy(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flxy

        subroutine flxz(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flxz
      END INTERFACE

c     Local variables

      real(8) :: hess(3,3,3),to(3,3),jac,ijac

c     Begin program

      if (have_map()) then
        hess = g_def%gmetric%grid(igr)%Gamma(i,j,k,:,:,:)

        jac = g_def%gmetric%grid(igr)%jac(i,j,k)
        ijac = 1d0/jac

        call flxx(g_def,i,j,k,nx,ny,nz,igr,igr,igr,alt_eom
     .           ,to(1,1),to(1,2),to(1,3),0)
        call flxy(g_def,i,j,k,nx,ny,nz,igr,igr,igr,alt_eom
     .           ,to(2,1),to(2,2),to(2,3),0)
        call flxz(g_def,i,j,k,nx,ny,nz,igr,igr,igr,alt_eom
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
      function div_tensor(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

      type(grid_mg_def),pointer :: g_def

      real(8) :: divt(3)

      integer :: i,j,k,igx,igy,igz,nx,ny,nz

      logical :: alt_eom

      logical,optional,intent(IN) :: vol

      INTERFACE
        subroutine flxx(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flxx

        subroutine flxy(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flxy

        subroutine flxz(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
          logical    :: alt_eom
        end subroutine flxz
      END INTERFACE

c     Local variables

      integer :: ig,jg,kg,im,jm,km

      real(8) :: dvol,dxx,dyy,dzz,jac,ijac

      real(8) :: tp(3,3),tm(3,3),to(3,3)

      real(8) :: msource(3)

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      im = i-1
      jm = j-1
      km = k-1

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = g_def%dxh(ig)
      dyy = g_def%dyh(jg)
      dzz = g_def%dzh(kg)

c     Fluxes

      tp(1,:) = dtnsr_flx_ijk(g_def,i ,j,k,nx,ny,nz,igx,alt_eom,flxx)
      tm(1,:) = dtnsr_flx_ijk(g_def,im,j,k,nx,ny,nz,igx,alt_eom,flxx)

      tp(2,:) = dtnsr_flx_ijk(g_def,i,j ,k,nx,ny,nz,igx,alt_eom,flxy)
      tm(2,:) = dtnsr_flx_ijk(g_def,i,jm,k,nx,ny,nz,igx,alt_eom,flxy)

      tp(3,:) = dtnsr_flx_ijk(g_def,i,j,k ,nx,ny,nz,igx,alt_eom,flxz)
      tm(3,:) = dtnsr_flx_ijk(g_def,i,j,km,nx,ny,nz,igx,alt_eom,flxz)

c     Source

      msource = dtnsr_src_ijk(g_def,i,j,k,nx,ny,nz,igx,alt_eom
     .                       ,flxx,flxy,flxz)

c     Assemble divergence of tensor

      divt = (tp(1,:)-tm(1,:))/dxx
     .      +(tp(2,:)-tm(2,:))/dyy
     .      +(tp(3,:)-tm(3,:))/dzz + msource

      if (alt_eom) then
        divt(2) = divt(2)/g_def%gmetric%grid(igx)%jac(i,j,k)
      endif

      !Volume factor
      if (vol_wgt) divt=divt*g_def%gmetric%grid(igx)%dlvol(i,j,k)  !Logical volume

c     End program

      end function div_tensor

c     dtnsr_cov_flx_ijk
c     ###############################################################
      function dtnsr_cov_flx_ijk(g_def,i,j,k,nx,ny,nz,igrid,flux)
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

      type(grid_mg_def),pointer :: g_def

      real(8) :: flx(3)

      integer :: i,j,k,igrid,nx,ny,nz

      INTERFACE
        subroutine flux(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flux
      END INTERFACE

c     Local variables

      call flux(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .         ,flx(1),flx(2),flx(3),1)

c     End program

      end function dtnsr_cov_flx_ijk

c     dtnsr_cov_src_ijk
c     ###############################################################
      function dtnsr_cov_src_ijk(g_def,i,j,k,nx,ny,nz,igr
     .                          ,flxx,flxy,flxz)
     .         result(src)
c     ---------------------------------------------------------------
c     Calculates source for contravariant components of div(tensor) in
c     general non-orthogonal coordinates at (i,j,k). The routines
c     flxx,flxy,flxz must provide the contravariant tersor components
c     T^lm (if flag /= 0, they actually provide T^lm/jac).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8) :: src(3)

      integer :: i,j,k,igr,nx,ny,nz

      INTERFACE
        subroutine flxx(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flxx

        subroutine flxy(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flxy

        subroutine flxz(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flxz
      END INTERFACE

c     Local variables

      real(8) :: hess(3,3,3),to(3,3)
      real(8) :: t11o,t12o,t13o
     .          ,t21o,t22o,t23o
     .          ,t31o,t32o,t33o

c     Begin program

      if (.not.have_map()) then
        src = 0d0
      else
        hess = g_def%gmetric%grid(igr)%Gamma(i,j,k,:,:,:)

        call flxx(g_def,i,j,k,nx,ny,nz,igr,igr,igr,to(1,1),to(1,2)
     .           ,to(1,3),0)
        call flxy(g_def,i,j,k,nx,ny,nz,igr,igr,igr,to(2,1),to(2,2)
     .           ,to(2,3),0)
        call flxz(g_def,i,j,k,nx,ny,nz,igr,igr,igr,to(3,1),to(3,2)
     .           ,to(3,3),0)

        to = transpose(to)
        src(1) =-sum(to*hess(:,:,1))
        src(2) =-sum(to*hess(:,:,2))
        src(3) =-sum(to*hess(:,:,3))
      endif

c     End program

      end function dtnsr_cov_src_ijk

c     div_tensor_cov
c     ###############################################################
      function div_tensor_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
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

      type(grid_mg_def),pointer :: g_def

      real(8) :: divt(3)

      integer :: i,j,k,igx,igy,igz,nx,ny,nz

      logical,optional,intent(IN) :: vol

      INTERFACE
        subroutine flxx(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flxx

        subroutine flxy(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flxy

        subroutine flxz(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,t11,t12,t13,flag)
          use grid_def
          type(grid_mg_def),pointer :: g_def
          integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
          real(8)    :: t11,t12,t13
        end subroutine flxz
      END INTERFACE

c     Local variables

      integer :: ig,jg,kg,im,jm,km

      real(8) :: dxx,dyy,dzz

      real(8) :: tp(3,3),tm(3,3),to(3,3)

      real(8) :: hess(3,3,3),msource(3)

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      im = i-1
      jm = j-1
      km = k-1

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = g_def%dxh(ig)
      dyy = g_def%dyh(jg)
      dzz = g_def%dzh(kg)

c     Fluxes

      tp(1,:) = dtnsr_cov_flx_ijk(g_def,i ,j,k,nx,ny,nz,igx,flxx)
      tm(1,:) = dtnsr_cov_flx_ijk(g_def,im,j,k,nx,ny,nz,igx,flxx)

      tp(2,:) = dtnsr_cov_flx_ijk(g_def,i,j ,k,nx,ny,nz,igx,flxy)
      tm(2,:) = dtnsr_cov_flx_ijk(g_def,i,jm,k,nx,ny,nz,igx,flxy)

      tp(3,:) = dtnsr_cov_flx_ijk(g_def,i,j,k ,nx,ny,nz,igx,flxz)
      tm(3,:) = dtnsr_cov_flx_ijk(g_def,i,j,km,nx,ny,nz,igx,flxz)

c     Source

      msource=dtnsr_cov_src_ijk(g_def,i,j,k,nx,ny,nz,igx,flxx,flxy,flxz)

c     Assemble divergence of tensor

      divt = (tp(1,:)-tm(1,:))/dxx
     .      +(tp(2,:)-tm(2,:))/dyy
     .      +(tp(3,:)-tm(3,:))/dzz + msource

      divt = divt/g_def%gmetric%grid(igx)%jac(i,j,k)

      !Volume factor
      if (vol_wgt) divt=divt*g_def%gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function div_tensor_cov

c     veclap_mesh
c     ###############################################################
      function veclap_mesh(g_def,igrid,cnv,diff) result(cnv_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cnv(0:,0:,0:,:)
     .          ,cnv_out(0:size(cnv,1)-1
     $                  ,0:size(cnv,2)-1
     $                  ,0:size(cnv,3)-1,3)

      real(8),optional :: diff(0:size(cnv,1)-1
     $                        ,0:size(cnv,2)-1
     $                        ,0:size(cnv,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cnv_out(i,j,k,:)=veclap_ijk(g_def,i,j,k,nx,ny,nz
     .                                 ,igrid,igrid,igrid,cnv
     .                                 ,diff=diff)
          enddo
        enddo
      enddo

c     End program

      end function veclap_mesh

c     veclap_ijk
c     ###############################################################
      function veclap_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                   ,diff,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*div(diff grad(vector)) at cell centers in general
c     non-orthog.  coordinates. Vector is assumed in contravariant
c     representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target,optional :: diff(0:nx+1,0:ny+1,0:nz+1)

      real(8) :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: icomp

      logical :: alt_eom

c     Begin program

      alt_eom = alt__eom()

      if (.not.have_map()) then
        do icomp=1,3
          vlap(icomp)=lap(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vfield(:,:,:,icomp),dff=diff
     .                   ,vol=vol)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        if (PRESENT(diff)) coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol)
        if (PRESENT(diff)) nullify(coef)
      endif

c     End program

      end function veclap_ijk

c     nabtensor_x
c     #############################################################
      subroutine nabtensor_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for nabla(vec)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: x,y,z,jac,jac0,jacp,ijac
        real(8)    :: nabla_v_cnv(3,3),dd

c     Begin program

        igrid = igx

        if (flag == 0) then
          ip = i
          jac = g_def%gmetric%grid(igrid)%jac(i,j,k)

          nabla_v_cnv = fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec,0,cnv=.true.)
        else
          ip = i+1
c$$$          jac = 0.25*(g_def%gmetric%grid(igrid)%jac_v(i,j  ,k  )
c$$$     .               +g_def%gmetric%grid(igrid)%jac_v(i,j-1,k  )
c$$$     .               +g_def%gmetric%grid(igrid)%jac_v(i,j  ,k-1)
c$$$     .               +g_def%gmetric%grid(igrid)%jac_v(i,j-1,k-1))
          jac  = 0.5*(g_def%gmetric%grid(igrid)%jac (ip,j,k)
     .               +g_def%gmetric%grid(igrid)%jac (i ,j,k))

          nabla_v_cnv = fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec,1,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*nabla_v_cnv(1,1)
        t12 = dd*nabla_v_cnv(1,2)
        t13 = dd*nabla_v_cnv(1,3)

c       Postprocessing

        if (flag /= 0) then
          if (isSP2(g_def,i+1,igx)) then
            ijac = 0d0
          else
            ijac = 1d0/jac
          endif

          t11 = t11*ijac
          if (.not.alteom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine nabtensor_x

c     nabtensor_y
c     #############################################################
      subroutine nabtensor_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v_cnv(3,3),dd

c     Begin program

        igrid = igx

        if (flag == 0) then
          jp = j
          jac= g_def%gmetric%grid(igrid)%jac(i,j,k)

          nabla_v_cnv = fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec,0,cnv=.true.)
        else
          jp = j+1
c$$$          jac= 0.25*(g_def%gmetric%grid(igrid)%jac_v(i  ,j,k  )
c$$$     .              +g_def%gmetric%grid(igrid)%jac_v(i-1,j,k  )
c$$$     .              +g_def%gmetric%grid(igrid)%jac_v(i  ,j,k-1)
c$$$     .              +g_def%gmetric%grid(igrid)%jac_v(i-1,j,k-1))
          jac  = 0.5*(g_def%gmetric%grid(igrid)%jac (i,jp,k)
     .               +g_def%gmetric%grid(igrid)%jac (i,j ,k))

          nabla_v_cnv = fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec,2,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j ,k))
        else
          dd = 1d0
        endif

        t21 = dd*nabla_v_cnv(2,1)
        t22 = dd*nabla_v_cnv(2,2)
        t23 = dd*nabla_v_cnv(2,3)

c     Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alteom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine nabtensor_y

c     nabtensor_z
c     #############################################################
      subroutine nabtensor_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v_cnv(3,3),dd

c     Begin program

        igrid = igx

        if (flag == 0) then
          kp = k
          jac= g_def%gmetric%grid(igrid)%jac(i,j,k)

          nabla_v_cnv = fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec,0,cnv=.true.)
        else
          kp = k+1
c$$$          jac= 0.25*(g_def%gmetric%grid(igrid)%jac_v(i  ,j  ,k)
c$$$     .              +g_def%gmetric%grid(igrid)%jac_v(i-1,j  ,k)
c$$$     .              +g_def%gmetric%grid(igrid)%jac_v(i  ,j-1,k)
c$$$     .              +g_def%gmetric%grid(igrid)%jac_v(i-1,j-1,k))
          jac  = 0.5*(g_def%gmetric%grid(igrid)%jac (i,j,kp)
     .               +g_def%gmetric%grid(igrid)%jac (i,j,k ))

          nabla_v_cnv = fnabla_v(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec,3,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k ))
        else
          dd = 1d0
        endif

        t31 = dd*nabla_v_cnv(3,1)
        t32 = dd*nabla_v_cnv(3,2)
        t33 = dd*nabla_v_cnv(3,3)

c     Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alteom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine nabtensor_z

c     veclap_cnv_mesh
c     ###############################################################
      function veclap_cnv_mesh(g_def,igrid,cnv,diff) result(cnv_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cnv(0:,0:,0:,:)
     .          ,cnv_out(0:size(cnv,1)-1
     $                  ,0:size(cnv,2)-1
     $                  ,0:size(cnv,3)-1,3)

      real(8),optional :: diff(0:size(cnv,1)-1
     $                        ,0:size(cnv,2)-1
     $                        ,0:size(cnv,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cnv_out(i,j,k,:)=veclap_cnv_ijk(g_def,i,j,k,nx,ny,nz
     .                                     ,igrid,igrid,igrid,cnv
     .                                     ,diff=diff)
          enddo
        enddo
      enddo

c     End program

      end function veclap_cnv_mesh

c     veclap_cnv_ijk
c     ###############################################################
      function veclap_cnv_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                      ,diff,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*div(diff grad(vector)) at cell centers in general
c     non-orthog.  coordinates. Vector is assumed in contravariant
c     representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target,optional :: diff(0:nx+1,0:ny+1,0:nz+1)

      real(8) :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer :: icomp

      logical :: alt_eom

c     Begin program

      alt_eom = alt__eom()

      if (.not.have_map()) then
        do icomp=1,3
          vlap(icomp)=lap(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vfield(:,:,:,icomp),dff=diff
     .                   ,vol=vol)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        if (PRESENT(diff)) coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     $                  ,nabtensor_cnv_x,nabtensor_cnv_y,nabtensor_cnv_z
     $                  ,vol=vol)
        if (PRESENT(diff)) nullify(coef)
      endif

c     End program

      end function veclap_cnv_ijk

c     nabtensor_cnv_x
c     #############################################################
      subroutine nabtensor_cnv_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,alt_eom,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t11-t13 for
c     nabla(vec), [nab(vec))^l_m, with vec in contravariant rep.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8) :: t11,t12,t13

        logical :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip
        real(8)    :: x,y,z,ijac,jac
        real(8)    :: nabla_v(3,3),dd

c     Begin program

        if (flag == 0) then
          ip = i

          jac = g_def%gmetric%grid(igx)%jac(i,j,k)

          nabla_v = fnabla_v_cnv(g_def,i,j,k,igx,vec,0)
        else
          ip = i+1

c$$$          jac = 0.25*(g_def%gmetric%grid(igx)%jac_v(i,j  ,k  )
c$$$     .               +g_def%gmetric%grid(igx)%jac_v(i,j-1,k  )
c$$$     .               +g_def%gmetric%grid(igx)%jac_v(i,j  ,k-1)
c$$$     .               +g_def%gmetric%grid(igx)%jac_v(i,j-1,k-1))
          jac  = 0.5*(g_def%gmetric%grid(igx)%jac (ip,j,k)
     .               +g_def%gmetric%grid(igx)%jac (i ,j,k))

          nabla_v = fnabla_v_cnv(g_def,i,j,k,igx,vec,1)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*nabla_v(1,1)
        t12 = dd*nabla_v(1,2)
        t13 = dd*nabla_v(1,3)

c       Postprocessing

        if (flag /= 0) then
          if (isSP2(g_def,i+1,igx)) then
            ijac = 1d0/SP_flsv
          else
            ijac = 1d0/jac
          endif

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine nabtensor_cnv_x

c     nabtensor_cnv_y
c     #############################################################
      subroutine nabtensor_cnv_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,alt_eom,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t21-t23 for
c     nabla(vec), [nab(vec))^l_m, with vec in contravariant rep.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8) :: t21,t22,t23

        logical :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp
        real(8)    :: x,y,z,ijac,jac
        real(8)    :: nabla_v(3,3),dd

c     Begin program

        if (flag == 0) then
          jp = j
          jac= g_def%gmetric%grid(igx)%jac(i,j,k)

          nabla_v = fnabla_v_cnv(g_def,i,j,k,igx,vec,0)
        else
          jp = j+1
c$$$          jac= 0.25*(g_def%gmetric%grid(igx)%jac_v(i  ,j,k  )
c$$$     .              +g_def%gmetric%grid(igx)%jac_v(i-1,j,k  )
c$$$     .              +g_def%gmetric%grid(igx)%jac_v(i  ,j,k-1)
c$$$     .              +g_def%gmetric%grid(igx)%jac_v(i-1,j,k-1))
          jac  = 0.5*(g_def%gmetric%grid(igx)%jac (i,jp,k)
     .               +g_def%gmetric%grid(igx)%jac (i,j ,k))

          nabla_v = fnabla_v_cnv(g_def,i,j,k,igx,vec,2)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t21 = dd*nabla_v(2,1)
        t22 = dd*nabla_v(2,2)
        t23 = dd*nabla_v(2,3)

c     Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine nabtensor_cnv_y

c     nabtensor_cnv_z
c     #############################################################
      subroutine nabtensor_cnv_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,alt_eom,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t31-t33 for
c     nabla(vec), [nab(vec))^l_m, with vec in contravariant rep.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8) :: t31,t32,t33

        logical :: alt_eom

c     Local variables

        integer :: ig,jg,kg,kp
        real(8) :: x,y,z,ijac,jac
        real(8) :: nabla_v(3,3),dd

c     Begin program

        if (flag == 0) then
          kp = k
          jac= g_def%gmetric%grid(igx)%jac(i,j,k)

          nabla_v = fnabla_v_cnv(g_def,i,j,k,igx,vec,0)
        else
          kp = k+1
c$$$          jac= 0.25*(g_def%gmetric%grid(igx)%jac_v(i  ,j  ,k)
c$$$     .              +g_def%gmetric%grid(igx)%jac_v(i-1,j  ,k)
c$$$     .              +g_def%gmetric%grid(igx)%jac_v(i  ,j-1,k)
c$$$     .              +g_def%gmetric%grid(igx)%jac_v(i-1,j-1,k))
          jac  = 0.5*(g_def%gmetric%grid(igx)%jac (i,j,kp)
     .               +g_def%gmetric%grid(igx)%jac (i,j,k ))

          nabla_v = fnabla_v_cnv(g_def,i,j,k,igx,vec,3)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t31 = dd*nabla_v(3,1)
        t32 = dd*nabla_v(3,2)
        t33 = dd*nabla_v(3,3)

c     Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine nabtensor_cnv_z

c     veclap_cov_mesh
c     ###############################################################
      function veclap_cov_mesh(g_def,igrid,cov,diff) result(cov_out)

c     ---------------------------------------------------------------
c     Calculates covariant components of lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cov_out(0:size(cov,1)-1
     $                  ,0:size(cov,2)-1
     $                  ,0:size(cov,3)-1,3)

      real(8),optional :: diff(0:size(cov,1)-1
     $                        ,0:size(cov,2)-1
     $                        ,0:size(cov,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(cov,1)-2
      ny = size(cov,2)-2
      nz = size(cov,3)-2

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cov_out(i,j,k,:)=veclap_cov_ijk(g_def,i,j,k,nx,ny,nz
     .                                     ,igrid,igrid,igrid
     .                                     ,cov,dff=diff)
          enddo
        enddo
      enddo

c     End program

      end function veclap_cov_mesh

c     veclap_cov_ijk
c     ###############################################################
      function veclap_cov_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                       ,dff,vol)  result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     is in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

      if (.not.have_map()) then
        do icomp=1,3
          vlap(icomp)=lap(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vfield(:,:,:,icomp),dff=dff
     .                   ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        if (PRESENT(dff)) coef => dff   !Pointer passed to nabtensor routines
        vlap = div_tensor_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
        if (PRESENT(dff)) nullify(coef)
      endif

c     End program

      end function veclap_cov_ijk

c     nabtensor_cov_x
c     #############################################################
      subroutine nabtensor_cov_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t11-t13 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13

c     Local variables

        integer    :: ig,jg,kg,ip
        real(8)    :: x,y,z
        real(8)    :: nabla_v(3,3),dd

c     Begin program

        if (flag == 0) then
          ip = i

          nabla_v = fnabla_v_mixed_cov(g_def,i,j,k,igx,vec,0)
        else
          ip = i+1
          
          nabla_v = fnabla_v_mixed_cov(g_def,i,j,k,igx,vec,1)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*nabla_v(1,1)
        t12 = dd*nabla_v(1,2)
        t13 = dd*nabla_v(1,3)

c     End program

      end subroutine nabtensor_cov_x

c     nabtensor_cov_y
c     #############################################################
      subroutine nabtensor_cov_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t21-t23 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23

c     Local variables

        integer    :: ig,jg,kg,jp
        real(8)    :: x,y,z
        real(8)    :: nabla_v(3,3),dd

c     Begin program

        if (flag == 0) then
          jp = j

          nabla_v = fnabla_v_mixed_cov(g_def,i,j,k,igx,vec,0)
        else
          jp = j+1

          nabla_v = fnabla_v_mixed_cov(g_def,i,j,k,igx,vec,2)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t21 = dd*nabla_v(2,1)
        t22 = dd*nabla_v(2,2)
        t23 = dd*nabla_v(2,3)

c     End program

      end subroutine nabtensor_cov_y

c     nabtensor_cov_z
c     #############################################################
      subroutine nabtensor_cov_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t31-t33 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33

c     Local variables

        integer    :: ig,jg,kg,kp
        real(8)    :: x,y,z
        real(8)    :: nabla_v(3,3),dd

c     Begin program

        if (flag == 0) then
          kp = k

          nabla_v = fnabla_v_mixed_cov(g_def,i,j,k,igx,vec,0)
        else
          kp = k+1

          nabla_v = fnabla_v_mixed_cov(g_def,i,j,k,igx,vec,3)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t31 = dd*nabla_v(3,1)
        t32 = dd*nabla_v(3,2)
        t33 = dd*nabla_v(3,3)

c     End program

      end subroutine nabtensor_cov_z

c$$$c     nabtensor_cov_x
c$$$c     #############################################################
c$$$      subroutine nabtensor_cov_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                      ,t11,t12,t13,flag)
c$$$c     -------------------------------------------------------------
c$$$c     Calculates mixed representation tensor components t11-t13 for
c$$$c     nabla(vec), [nab(vec))^l_m.
c$$$c     -------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
c$$$        real(8)    :: t11,t12,t13
c$$$
c$$$c     Local variables
c$$$
c$$$        integer    :: ig,jg,kg,ip,igrid
c$$$        real(8)    :: x,y,z
c$$$        real(8)    :: nabla_v(3,3),dd,gsuper(3,3)
c$$$
c$$$c     Begin program
c$$$
c$$$        igrid = igx
c$$$
c$$$        if (flag == 0) then
c$$$          ip = i
c$$$          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)
c$$$
c$$$          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                          ,vec(:,:,:,1)
c$$$     .                          ,vec(:,:,:,2)
c$$$     .                          ,vec(:,:,:,3),0)
c$$$        else
c$$$          ip = i+1
c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k-1,:,:))
c$$$c$$$          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:)
c$$$c$$$     .                 +g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:))
c$$$
c$$$          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                          ,vec(:,:,:,1)
c$$$     .                          ,vec(:,:,:,2)
c$$$     .                          ,vec(:,:,:,3),1)
c$$$        endif
c$$$
c$$$        if (associated(coef)) then
c$$$          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
c$$$        else
c$$$          dd = 1d0
c$$$        endif
c$$$
c$$$        t11 = dd*( gsuper(1,1)*nabla_v(1,1)
c$$$     .            +gsuper(1,2)*nabla_v(2,1)
c$$$     .            +gsuper(1,3)*nabla_v(3,1) )
c$$$
c$$$        t12 = dd*( gsuper(1,1)*nabla_v(1,2)
c$$$     .            +gsuper(1,2)*nabla_v(2,2)
c$$$     .            +gsuper(1,3)*nabla_v(3,2) )
c$$$
c$$$        t13 = dd*( gsuper(1,1)*nabla_v(1,3)
c$$$     .            +gsuper(1,2)*nabla_v(2,3)
c$$$     .            +gsuper(1,3)*nabla_v(3,3) )
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine nabtensor_cov_x
c$$$
c$$$c     nabtensor_cov_y
c$$$c     #############################################################
c$$$      subroutine nabtensor_cov_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                      ,t21,t22,t23,flag)
c$$$c     -------------------------------------------------------------
c$$$c     Calculates mixed representation tensor components t21-t23 for
c$$$c     nabla(vec), [nab(vec))^l_m.
c$$$c     -------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
c$$$        real(8)    :: t21,t22,t23
c$$$
c$$$c     Local variables
c$$$
c$$$        integer    :: ig,jg,kg,jp,igrid
c$$$        real(8)    :: x,y,z
c$$$        real(8)    :: nabla_v(3,3),gsuper(3,3),dd
c$$$
c$$$c     Begin program
c$$$
c$$$        igrid = igx
c$$$
c$$$        if (flag == 0) then
c$$$          jp = j
c$$$          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)
c$$$
c$$$          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                          ,vec(:,:,:,1)
c$$$     .                          ,vec(:,:,:,2)
c$$$     .                          ,vec(:,:,:,3),0)
c$$$        else
c$$$          jp = j+1
c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k-1,:,:))
c$$$c$$$          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:)
c$$$c$$$     .                 +g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:))
c$$$
c$$$          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                          ,vec(:,:,:,1)
c$$$     .                          ,vec(:,:,:,2)
c$$$     .                          ,vec(:,:,:,3),2)
c$$$        endif
c$$$
c$$$        if (associated(coef)) then
c$$$          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j,k))
c$$$        else
c$$$          dd = 1d0
c$$$        endif
c$$$
c$$$        t21 = dd*( gsuper(2,1)*nabla_v(1,1)
c$$$     .            +gsuper(2,2)*nabla_v(2,1)
c$$$     .            +gsuper(2,3)*nabla_v(3,1) )
c$$$
c$$$        t22 = dd*( gsuper(2,1)*nabla_v(1,2)
c$$$     .            +gsuper(2,2)*nabla_v(2,2)
c$$$     .            +gsuper(2,3)*nabla_v(3,2) )
c$$$
c$$$        t23 = dd*( gsuper(2,1)*nabla_v(1,3)
c$$$     .            +gsuper(2,2)*nabla_v(2,3)
c$$$     .            +gsuper(2,3)*nabla_v(3,3) )
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine nabtensor_cov_y
c$$$
c$$$c     nabtensor_cov_z
c$$$c     #############################################################
c$$$      subroutine nabtensor_cov_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                      ,t31,t32,t33,flag)
c$$$c     -------------------------------------------------------------
c$$$c     Calculates mixed representation tensor components t31-t33 for
c$$$c     nabla(vec), [nab(vec))^l_m.
c$$$c     -------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
c$$$        real(8)    :: t31,t32,t33
c$$$
c$$$c     Local variables
c$$$
c$$$        integer    :: ig,jg,kg,kp,igrid
c$$$        real(8)    :: x,y,z
c$$$        real(8)    :: nabla_v(3,3),gsuper(3,3),dd
c$$$
c$$$c     Begin program
c$$$
c$$$        igrid = igx
c$$$
c$$$        if (flag == 0) then
c$$$          kp = k
c$$$          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)
c$$$
c$$$          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                          ,vec(:,:,:,1)
c$$$     .                          ,vec(:,:,:,2)
c$$$     .                          ,vec(:,:,:,3),0)
c$$$        else
c$$$          kp = k+1
c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j-1,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j-1,k,:,:))
c$$$c$$$          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:)
c$$$c$$$     .                 +g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:))
c$$$
c$$$          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                          ,vec(:,:,:,1)
c$$$     .                          ,vec(:,:,:,2)
c$$$     .                          ,vec(:,:,:,3),3)
c$$$        endif
c$$$
c$$$        if (associated(coef)) then
c$$$          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k))
c$$$        else
c$$$          dd = 1d0
c$$$        endif
c$$$
c$$$        t31 = dd*( gsuper(3,1)*nabla_v(1,1)
c$$$     .            +gsuper(3,2)*nabla_v(2,1)
c$$$     .            +gsuper(3,3)*nabla_v(3,1) )
c$$$
c$$$        t32 = dd*( gsuper(3,1)*nabla_v(1,2)
c$$$     .            +gsuper(3,2)*nabla_v(2,2)
c$$$     .            +gsuper(3,3)*nabla_v(3,2) )
c$$$
c$$$        t33 = dd*( gsuper(3,1)*nabla_v(1,3)
c$$$     .            +gsuper(3,2)*nabla_v(2,3)
c$$$     .            +gsuper(3,3)*nabla_v(3,3) )
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine nabtensor_cov_z

c     veclap_vrtx
c     ###############################################################
      function veclap_vrtx(g_def,igrid,cov,ret_cnv) result(veclap)
c     ---------------------------------------------------------------
c     Calculates covariant components of vector Laplacian for
c     covariant vector representation satisfying div(lap(A)) =
c     lap(div(A)) at faces everywhere on the mesh.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,veclap(0:size(cov,1)-1
     .                 ,0:size(cov,2)-1
     .                 ,0:size(cov,3)-1
     .                 ,3)

      logical,optional :: ret_cnv

c     Local variables

      real(8) :: cnv(0:size(cov,1)-1
     .              ,0:size(cov,2)-1
     .              ,0:size(cov,3)-1,3)

      logical :: tocov

c     Begin program

      if (PRESENT(ret_cnv)) then
        tocov = .not.ret_cnv
      else
        tocov = .false.
      endif

      cnv = XformToCnv(g_def,igrid,cov)

c     Assemble curl(curl) on the mesh

      veclap =-curlcurl_vrtx_mesh(g_def,igrid,cov)

      if (tocov) veclap = XformToCov(g_def,igrid,veclap)

c     Assemble vector Laplacian on the mesh

      veclap = veclap
     .        +gradDiv(g_def,igrid,cnv,ret_cnv=.not.tocov)
cc     .        +gradDiv_vrtx(g_def,igrid,cnv)
      
c     End program

      end function veclap_vrtx

c     curlcurl_cov_mesh
c     ###############################################################
      function curlcurl_cov_mesh(g_def,igrid,cov) result(cov_out)

c     ---------------------------------------------------------------
c     Calculates curl(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is both covariant
c     and contravariant, and returns the covariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,cov_out(0:size(cov,1)-1
     $                  ,0:size(cov,2)-1
     $                  ,0:size(cov,3)-1,3)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cov_out(i,j,k,:)=curlcurl_cov_ijk(g_def,i,j,k,nx,ny,nz
     $                                       ,igrid,igrid,igrid,cov)
          enddo
        enddo
      enddo

c     End program

      end function curlcurl_cov_mesh

c     curlcurl_cov_ijk
c     ###############################################################
      function curlcurl_cov_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,acov,vol)  result(cov)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal coordinates,
c     preserving the SPD property. The vector A is both in covariant
c     representation, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: cov(3)

      real(8),target :: acov(0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional,intent(IN) :: vol

c     Local variables

c     Begin program

      vec => acov

      cov = div_tensor_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                ,curlcurl_cov_x,curlcurl_cov_y,curlcurl_cov_z
     .                ,vol=vol)

      nullify(vec)

c     End program

      end function curlcurl_cov_ijk

c     curlcurl_cov_x
c     #############################################################
      subroutine curlcurl_cov_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 of curl(curl) in covariant
c     representation.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        if (flag == 0) then
          ip = i
          gsuper = g_def%gmetric%grid(igrid)%gsup  (i,j  ,k  ,:,:)

          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        else
          ip = i+1
c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j  ,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i,j-1,k-1,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i ,j,k,:,:))

          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),1)
        endif

        nabla_v = -nabla_v + transpose(nabla_v)

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
      subroutine curlcurl_cov_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 of curl(curl) in covariant
c     representation.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        if (flag == 0) then
          jp = j
          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        else
          jp = j+1
c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k  ,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j,k-1,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j,k-1,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j ,k,:,:))

          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),2)
        endif

        nabla_v = -nabla_v + transpose(nabla_v)

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
      subroutine curlcurl_cov_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 of curl(curl) in covariant
c     representation.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        if (flag == 0) then
          kp = k
          gsuper = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        else
          kp = k+1
c$$$          gsuper =0.25*(g_def%gmetric%grid(igrid)%gsup_v(i  ,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j  ,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i  ,j-1,k,:,:)
c$$$     .                 +g_def%gmetric%grid(igrid)%gsup_v(i-1,j-1,k,:,:))
          gsuper = 0.5*(g_def%gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .                 +g_def%gmetric%grid(igrid)%gsup(i,j,k ,:,:))

          nabla_v = fnabla_v_cov(g_def,i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),3)
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
      function curlcurl_mesh(g_def,igrid,vec,iscnv) result(cnv_out)

c     ---------------------------------------------------------------
c     Calculates curlcurl(A) in general non-orthogonal coordinates. The
c     vector A can be covariant or contravariant, and returns
c     contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: vec(0:,0:,0:,:)
     .          ,cnv_out(0:size(vec,1)-1
     $                  ,0:size(vec,2)-1
     $                  ,0:size(vec,3)-1,3)

c$$$      real(8), optional :: diff(0:size(vec,1)-1
c$$$     $                         ,0:size(vec,2)-1
c$$$     $                         ,0:size(vec,3)-1)

      logical :: iscnv

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            cnv_out(i,j,k,:)=curlcurl_ijk(g_def,i,j,k,nx,ny,nz
     $                                   ,igrid,igrid,igrid,vec,iscnv)
          enddo
        enddo
      enddo

c     End program

      end function curlcurl_mesh

c     curlcurl_ijk
c     ###############################################################
      function curlcurl_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vec,iscnv
     .                     ,vol)  result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A can be
c     covariant or contravariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vec(0:nx+1,0:ny+1,0:nz+1,3)
cc      real(8),target,optional :: diff(0:nx+1,0:ny+1,0:nz+1)

      logical :: iscnv
      
      logical,optional,intent(IN) :: vol

c     Local variables

      real(8) :: cnv(3),cov(3)

c     Begin program

      if (iscnv) then
        cnv=-veclap (g_def,i,j,k,nx,ny,nz,igx,igy,igz,vec,vol=vol)
     .      +gradDiv(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vec,vol=vol
     .              ,cnv=.true.)
      else
        cov = curlcurl_cov_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,vec)
        cnv = XformToCnv(g_def,i,j,k,igx,cov)
      endif
      
c     End program

      end function curlcurl_ijk

c     veclap_stg_mesh
c     ###############################################################
      function veclap_stg_mesh(g_def,igrid,cov,vol,bc_order,apply_bc)
     .         result(veclap)
c     ---------------------------------------------------------------
c     Calculates contravariant components of vector Laplacian for
c     covariant vector representation satisfying div(lap(A)) =
c     lap(div(A)) at faces everywhere on the mesh.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: cov(0:,0:,0:,:)
     .          ,veclap(0:size(cov,1)-1
     .                 ,0:size(cov,2)-1
     .                 ,0:size(cov,3)-1
     .                 ,3)

      logical,optional,intent(IN) :: vol,apply_bc
      integer,optional :: bc_order

c     Local variables

      integer :: bcs(6,3),nx,ny,nz,i,j,k
      logical :: vol_wgt,do_bc
      real(8) :: divA(0:size(cov,1)-1
     .               ,0:size(cov,2)-1
     .               ,0:size(cov,3)-1)

      real(8),pointer,dimension(:,:,:) :: div

c     Begin program

      nx = size(cov,1)-2
      ny = size(cov,2)-2
      nz = size(cov,3)-2

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      do_bc = .false.
      if (PRESENT(apply_bc)) do_bc = apply_bc
  
c     Compute div(cov)

      divA = div_stg(g_def,igrid,cov,iscnv=.false.,bc=.true.
     .              ,bc_order=bc_order)

c     Assemble curl(curl) on the mesh

      veclap =-curlcurl_stg_mesh(g_def,igrid,cov)

c$$$      !Old way
c$$$      veclap = 0d0
c$$$      do k=0,nz
c$$$        do j=0,ny
c$$$          do i=0,nx
c$$$            veclap(i,j,k,:)= -curlcurl_stg_ijk(g_def,i,j,k,nx,ny,nz
c$$$     .                                        ,igrid,igrid,igrid,cov)
c$$$          enddo
c$$$        enddo
c$$$      enddo

c     Assemble vector Laplacian on the mesh

      do k=1,nz
        do j=1,ny
          do i=1,nx
            !Add grad(div) term
            veclap(i,j,k,1) = lap_flx_x_ijk(g_def,i,j,k,igrid,divA)
     .                       +veclap(i,j,k,1)
            veclap(i,j,k,2) = lap_flx_y_ijk(g_def,i,j,k,igrid,divA)
     .                       +veclap(i,j,k,2)
            veclap(i,j,k,3) = lap_flx_z_ijk(g_def,i,j,k,igrid,divA)
     .                       +veclap(i,j,k,3)
            if (vol_wgt) veclap(i,j,k,:)
     .                  =veclap(i,j,k,:)
     .                  *g_def%gmetric%grid(igrid)%dlvol(i,j,k)
          enddo
        enddo
      enddo

      if (do_bc) then
        call get_topol_bcs(bcs(:,1))
        call get_topol_bcs(bcs(:,2))
        call get_topol_bcs(bcs(:,3))
        where(bcs == DEF) bcs = EXT
        call setMGBC_stg(g_def,0,3,nx,ny,nz,igrid,veclap
     .                  ,bcs,is_cnv=.true.,iorder=bc_order)
      endif

c     End program

      end function veclap_stg_mesh

c     curlcurl_stg_ijk
c     ###############################################################
      function curlcurl_stg_ijk(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cov
     .                         ,vol)  result (curlcurl)
c     ---------------------------------------------------------------
c     Calculates contravariant components of curl(curl) operator for
c     covariant vector representation at faces on a per-cell basis.
c
c     NOTE: ONLY WORKS IN 2D FOR NOW
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: cov(0:nx+1,0:ny+1,0:nz+1,3),curlcurl(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      real(8) :: Bz_ipjp,Bz_imjp,Bz_ipjm
      integer :: ig,jg,kg
      logical :: vol_wgt
      real(8),pointer,dimension(:,:,:,:,:) :: gsupv
      real(8),pointer,dimension(:,:,:)     :: jacv,jac

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      !Code for tensor-product mesh
      gsupv => g_def%gmetric%grid(igx)%gsup_v
      jacv  => g_def%gmetric%grid(igx)%jac_v
      jac   => g_def%gmetric%grid(igx)%jac

      Bz_ipjp = 1d0/jacv(i,j  ,k)
     .          *((cov(i+1,j  ,k,2)-cov(i,j  ,k,2))/g_def%dx(ig)
     .           -(cov(i  ,j+1,k,1)-cov(i,j  ,k,1))/g_def%dy(jg))
      Bz_ipjm = 1d0/jacv(i,j-1,k)
     .          *((cov(i+1,j-1,k,2)-cov(i,j-1,k,2))/g_def%dx(ig)
     .           -(cov(i  ,j  ,k,1)-cov(i,j-1,k,1))/g_def%dy(jg))
      if (isSP2(g_def,i,igx)) then
        Bz_imjp = (cov(i  ,j,k,2)/jac(i  ,j,k)
     .            -cov(i-1,j,k,2)/jac(i-1,j,k))/g_def%dx(ig)
      else
        Bz_imjp = 1d0/jacv(i-1,j,k)
     .            *((cov(i  ,j  ,k,2)-cov(i-1,j,k,2))/g_def%dx(ig)
     .             -(cov(i-1,j+1,k,1)-cov(i-1,j,k,1))/g_def%dy(jg))
      endif

      curlcurl(1) = (Bz_ipjp - Bz_ipjm)/g_def%dy(jg)

      curlcurl(2) =-(Bz_ipjp - Bz_imjp)/g_def%dx(ig)

      curlcurl(3) =-lap(g_def,i,j,k,nx,ny,nz,igx,igy,igz,cov(:,:,:,3))
     .             *jac(i,j,k)

      if(vol_wgt) curlcurl=curlcurl*g_def%gmetric%grid(igx)%dlvol(i,j,k)

c     End program

      end function curlcurl_stg_ijk

c     curlcurl_stg_mesh
c     ###############################################################
      function curlcurl_stg_mesh(g_def,igr,cov,vol) result(curlcurl)
c     ---------------------------------------------------------------
c     Calculates contravariant components of curl(curl) operator for
c     covariant vector representation at faces on a per-cell basis.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr,bcs(6)

      real(8) :: cov(0:,0:,0:,:)
     .          ,curlcurl(0:size(cov,1)-1
     .                   ,0:size(cov,2)-1
     .                   ,0:size(cov,3)-1,3)

      logical,optional,intent(IN) :: vol

c     Local variables

      logical :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      !Find first curl (vertex-based, contravariant)
      curlcurl = curl_stg_mesh(g_def,igr,cov)

      !Transform to covariant
      curlcurl = XformToCov_stg_mesh(g_def,igr,curlcurl
     .                              ,vertex_in =.true.
     .                              ,vertex_out=.true.)

      !Compute second curl (cell-based, contravariant)
      curlcurl = curl_stg_mesh(g_def,igr,curlcurl,vrtx=.true.)

      if (vol_wgt) then
       curlcurl(:,:,:,1)=curlcurl(:,:,:,1)*g_def%gmetric%grid(igr)%dlvol
       curlcurl(:,:,:,2)=curlcurl(:,:,:,2)*g_def%gmetric%grid(igr)%dlvol
       curlcurl(:,:,:,3)=curlcurl(:,:,:,3)*g_def%gmetric%grid(igr)%dlvol
      endif

c     End program

      end function curlcurl_stg_mesh

      end module grid_tensor_ops

c grid_advec_ops
c ######################################################################
      module grid_advec_ops

        use grid_create

        use grid_tensor_ops

c     Options:
c       * advect = 1 => First-order upwind
c       * advect = 2 => Centered
c       * advect = 3 => High-order upwind (QUICK)
c       * advect = 4 => Monotone high-order upwind (SMART)
c       * advect = 5 => SMART with smooth transition
c       * advect = 6 => Centered, high-order
c       * advect = 7 => Gamma (modification of SMART; not operational)

        integer,parameter :: UPWD=1,ZIP=2,QUICK=3,SMART=4,SSMART=5,CHO=6

        integer :: adv_meth

        integer,private :: bc_adv_meth=ZIP   !Advective method at boundaries

        !Smooth SMART config
        real(8),private :: ssmart_xp1=0.25 !0.3
     .                    ,ssmart_xp2=0.52 !0.5

        !SMART config
        real(8),private :: smart_slp1= 1.5  !3.0
     .                    ,smart_slp2= 0.5  !0.3

        !ADVECTION positivity threshold for ZIP fluxes
        real(8),private :: zip_threshold=1d-8
        
        INTERFACE inter
          module procedure in2face_scl,in2face_vec
cc     .                    ,in2face_scl_new,in2face_vec_new
        end INTERFACE

      contains

c     set_bc_advc_meth
c     ###############################################################
      subroutine set_bc_advc_meth(bcmth)

      implicit none

c     ---------------------------------------------------------------
c     Sets advective method at boundaries
c     ---------------------------------------------------------------

      integer :: bcmth

      bc_adv_meth = bcmth

      end subroutine set_bc_advc_meth

c     scl_advc_flx_x
c     ###############################################################
      function scl_advc_flx_x(g_def,i,j,k,nx,ny,nz,igr,v,arr
     .                       ,method,sp,zip_vel) result(flxip)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_1 of scalar arr at face
c     (i+1/2,j,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1),flxip

      logical,optional :: sp,zip_vel

c     Local variables

      integer :: im,ip,ipp,igm,igp,igpp,ig,jg,kg,admth

      real(8) :: jacip,jacim,jac,qip,vxp,vxh

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

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      ip  = i+1
      im  = max(i-1,0)
      ipp = min(ip+1,nx+1)

      igp  = ig + (ip-i)
      igm  = ig + (im-i)
      igpp = ig + (ipp-i)

c     Grid info

      jac = g_def%gmetric%grid(igr)%jac(i,j,k)

      if (spoint) then
        jacip = g_def%gmetric%grid(igr)%jac(ip,j,k)
      else
        jacip = jac
      endif

c     Find X ZIP flux (used even if method /= ZIP)

      flxip = 0d0

      at_bdry = isBdry(g_def,i,igr,2) .or. isBdry(g_def,i+1,igr,1)

      if (.not.at_bdry) then               !Not at boundary
        flxip = 0.25*(v(ip,j,k)*arr(i ,j,k)/jacip
     .              + v(i ,j,k)*arr(ip,j,k)/jac  )*(jac+jacip)
      else
        if (   (v(ip,j,k) == -v(i,j,k))
     .     .or.(isSP2(g_def,i+1,igr,ibc=1).and.spoint)) then
          flxip = 0d0    !Shut off flux for vn=0 BC
        else
          flxip = 0.5*(v(ip,j,k)*arr(i,j,k)+v(i,j,k)*arr(ip,j,k))
        endif
      endif

c     Other advective methods

      if (method /= 2.and.(.not.(at_bdry).or.bc_safe)) then

        !Find velocities at faces from ZIP fluxes
        if (zv
     .     .and.(arr(ip,j,k)*arr(i,j,k)>0d0)
     .     .and.(abs(arr(ip,j,k)+arr(i,j,k))>zip_threshold)) then
          vxp =2.*flxip/(arr(ip,j,k)+arr(i,j,k))
        else
          vxp = 0.5*(v(ip,j,k)+v(i,j,k))
        endif

        !Find advected quantity at faces
        if ((i == nx .or. i == 0).and.(.not.bc_safe)) then
          admth = bc_adv_meth
        else
          admth = method
        endif

        qip=in2face_scl(g_def%xx(igpp),arr(ipp,j,k)
     .                 ,g_def%xx(igp ),arr(ip ,j,k)
     .                 ,g_def%xx(ig  ),arr(i  ,j,k)
     .                 ,g_def%xx(igm ),arr(im ,j,k)
     .                 ,vxp,admth)

        !Setup fluxes (ZIP at boundaries by default for conservation)
        flxip=vxp*qip

      endif

c     End program

      end function scl_advc_flx_x

c     scl_advc_flx_y
c     ###############################################################
      function scl_advc_flx_y(g_def,i,j,k,nx,ny,nz,igr,v,arr
     .                       ,method,zip_vel) result(flxjp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_2 of scalar arr at face
c     (i,j+1/2,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      jp  = j+1
      jm  = max(j-1,0)
      jpp = min(jp+1,ny+1)

      jgp  = jg + (jp-j)
      jgm  = jg + (jm-j)
      jgpp = jg + (jpp-j)

c     Grid info

      jac = g_def%gmetric%grid(igr)%jac(i,j,k)

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = isBdry(g_def,j,igr,4) .or. isBdry(g_def,j+1,igr,3)

      flxjp = 0.5*(v(i,jp,k)*arr(i,j,k) + v(i,j,k)*arr(i,jp,k))

c     Other advective methods

      if (method /= 2.and.(.not.(at_bdry).or.bc_safe)) then

        !Find velocities at faces from ZIP fluxes
        if (zv
     .     .and.(arr(i,jp,k)*arr(i,j,k)>0d0)
     .     .and.(abs(arr(i,jp,k)+arr(i,j,k))>zip_threshold)) then
          vyp =2.*flxjp/(arr(i,jp,k)+arr(i,j,k))
        else
          vyp = 0.5*(v(i,jp,k)+v(i,j,k))
        endif

        !Find advected quantity at faces
        if ((j == ny .or. j == 0).and.(.not.bc_safe)) then
          admth = bc_adv_meth
        else
          admth = method
        endif

        qjp=in2face_scl(g_def%yy(jgpp),arr(i,jpp,k)
     .                 ,g_def%yy(jgp ),arr(i,jp ,k)
     .                 ,g_def%yy(jg  ),arr(i,j  ,k)
     .                 ,g_def%yy(jgm ),arr(i,jm ,k)
     .                 ,vyp,admth)

        !Setup fluxes (ZIP at boundaries by default for conservation)
        flxjp=vyp*qjp

      endif

c     End program

      end function scl_advc_flx_y

c     scl_advc_flx_z
c     ###############################################################
      function scl_advc_flx_z(g_def,i,j,k,nx,ny,nz,igr,v,arr
     .                       ,method,zip_vel) result(flxkp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_3 of scalar arr at face
c     (i,j,k+1/2) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      kp  = k+1
      km  = max(k-1,0)
      kpp = min(kp+1,nz+1)

      kgp  = kg + (kp-k)
      kgm  = kg + (km-k)
      kgpp = kg + (kpp-k)

c     Grid info

      jac = g_def%gmetric%grid(igr)%jac(i,j,k)

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = isBdry(g_def,k,igr,6) .or. isBdry(g_def,k+1,igr,5)

      flxkp = 0.5*(v(i,j,kp)*arr(i,j,k) + v(i,j,k)*arr(i,j,kp))

c     Other advective methods

      if (method /= 2.and.(.not.(at_bdry).or.bc_safe)) then

        !Find velocities at faces from ZIP fluxes
        if (zv
     .      .and.(arr(i,j,kp)*arr(i,j,k)>0d0)
     .      .and.abs(arr(i,j,kp)+arr(i,j,k))>zip_threshold) then
          vzp =2.*flxkp/(arr(i,j,kp)+arr(i,j,k))
        else
          vzp = 0.5*(v(i,j,kp)+v(i,j,k))
        endif

        !Find advected quantity at faces
        if ((k == nz .or. k == 0).and.(.not.bc_safe)) then
          admth = bc_adv_meth
        else
          admth = method
        endif

        qkp=in2face_scl(g_def%zz(kgpp),arr(i,j,kpp)
     .                 ,g_def%zz(kgp ),arr(i,j,kp )
     .                 ,g_def%zz(kg  ),arr(i,j,k  )
     .                 ,g_def%zz(kgm ),arr(i,j,km )
     .                 ,vzp,admth)

        !Setup fluxes (ZIP at boundaries by default for conservation)
        flxkp=vzp*qkp

      endif

c     End program

      end function scl_advc_flx_z

c     flx_advec
c     ###############################################################
      function flx_advec(g_def,i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                  ,method,sp,vol,zip_vel) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = g_def%dxh(ig)
      dyy = g_def%dyh(jg)
      dzz = g_def%dzh(kg)

      jac = g_def%gmetric%grid(igx)%jac(i,j,k)

c     Find fluxes

      flxip = scl_advc_flx_x(g_def,i  ,j,k,nx,ny,nz,igx,v(:,:,:,1),arr
     .                      ,method,sp=sp,zip_vel=zip_vel)
      flxim = scl_advc_flx_x(g_def,i-1,j,k,nx,ny,nz,igx,v(:,:,:,1),arr
     .                      ,method,sp=sp,zip_vel=zip_vel)

      flxjp = scl_advc_flx_y(g_def,i,j  ,k,nx,ny,nz,igx,v(:,:,:,2),arr
     .                      ,method,zip_vel=zip_vel)
      flxjm = scl_advc_flx_y(g_def,i,j-1,k,nx,ny,nz,igx,v(:,:,:,2),arr
     .                      ,method,zip_vel=zip_vel)

      flxkp = scl_advc_flx_z(g_def,i,j,k  ,nx,ny,nz,igx,v(:,:,:,3),arr
     .                      ,method,zip_vel=zip_vel)
      flxkm = scl_advc_flx_z(g_def,i,j,k-1,nx,ny,nz,igx,v(:,:,:,3),arr
     .                      ,method,zip_vel=zip_vel)

c     Find advective term

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac

      if (vw) advec = advec*g_def%gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function flx_advec

c     flx_nc_advec
c     ###############################################################
      function flx_nc_advec(g_def,i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                     ,method,sp,vol,reverse) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igx,igy,igz,method

      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
     .          ,v  (0:nx+1,0:ny+1,0:nz+1,3),advec

      logical,optional :: sp,vol,reverse

c     Local variables

      integer :: igrid,half_elem,su

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

      bc_safe = .not.(method > ZIP .and. method <= 6)

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

      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

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

      igl = i+g_def%ilo(igx)-1
      jgl = j+g_def%jlo(igy)-1
      kgl = k+g_def%klo(igz)-1

c     Grid info

      dxx = g_def%dxh(ig)
      dyy = g_def%dyh(jg)
      dzz = g_def%dzh(kg)

      jac = g_def%gmetric%grid(igrid)%jac(i,j,k)

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

      if (method /= ZIP) then

        vxp = sgn*v(i,j,k,1)
        vxm = sgn*v(i,j,k,1)
                                                    
        vyp = sgn*v(i,j,k,2)
        vym = sgn*v(i,j,k,2)
                                                    
        vzp = sgn*v(i,j,k,3)
        vzm = sgn*v(i,j,k,3)

        !Find advected quantity at faces 
        if ((i == nx).and.(.not.bc_safe)) then
          qip=inter(g_def%xx(igpp),arr(ipp,j,k)
     .             ,g_def%xx(igp ),arr(ip ,j,k)
     .             ,g_def%xx(ig  ),arr(i  ,j,k)
     .             ,g_def%xx(igm ),arr(im ,j,k)
     .             ,vxp,bc_adv_meth)
        else
          qip=inter(g_def%xx(igpp),arr(ipp,j,k)
     .             ,g_def%xx(igp ),arr(ip ,j,k)
     .             ,g_def%xx(ig  ),arr(i  ,j,k)
     .             ,g_def%xx(igm ),arr(im ,j,k)
     .             ,vxp,method)
        endif

        if ((i == 1).and.(.not.bc_safe)) then
          qim=inter(g_def%xx(igp ),arr(ip ,j,k)
     .             ,g_def%xx(ig  ),arr(i  ,j,k)
     .             ,g_def%xx(igm ),arr(im ,j,k)
     .             ,g_def%xx(igmm),arr(imm,j,k)
     .             ,vxm,bc_adv_meth)
        else
          qim=inter(g_def%xx(igp ),arr(ip ,j,k)
     .             ,g_def%xx(ig  ),arr(i  ,j,k)
     .             ,g_def%xx(igm ),arr(im ,j,k)
     .             ,g_def%xx(igmm),arr(imm,j,k)
     .             ,vxm,method)
        endif

        if ((j == ny).and.(.not.bc_safe)) then
          qjp=inter(g_def%yy(jgpp),arr(i,jpp,k)
     .             ,g_def%yy(jgp ),arr(i,jp ,k)
     .             ,g_def%yy(jg  ),arr(i,j  ,k)
     .             ,g_def%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,bc_adv_meth)
        else
          qjp=inter(g_def%yy(jgpp),arr(i,jpp,k)
     .             ,g_def%yy(jgp ),arr(i,jp ,k)
     .             ,g_def%yy(jg  ),arr(i,j  ,k)
     .             ,g_def%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,method)
        endif

        if ((j == 1).and.(.not.bc_safe)) then
          qjm=inter(g_def%yy(jgp ),arr(i,jp ,k)
     .             ,g_def%yy(jg  ),arr(i,j  ,k)
     .             ,g_def%yy(jgm ),arr(i,jm ,k)
     .             ,g_def%yy(jgmm),arr(i,jmm,k)
     .             ,vym,bc_adv_meth)
        else
          qjm=inter(g_def%yy(jgp ),arr(i,jp,k )
     .             ,g_def%yy(jg  ),arr(i,j  ,k)
     .             ,g_def%yy(jgm ),arr(i,jm,k )
     .             ,g_def%yy(jgmm),arr(i,jmm,k)
     .             ,vym,method)
        endif

        if ((k == nz).and.(.not.bc_safe)) then
          qkp=inter(g_def%zz(kgpp),arr(i,j,kpp)
     .             ,g_def%zz(kgp ),arr(i,j,kp )
     .             ,g_def%zz(kg  ),arr(i,j,k  )
     .             ,g_def%zz(kgm ),arr(i,j,km )
     .             ,vzp,bc_adv_meth)
        else
          qkp=inter(g_def%zz(kgpp),arr(i,j,kpp)
     .             ,g_def%zz(kgp ),arr(i,j,kp )
     .             ,g_def%zz(kg  ),arr(i,j,k  )
     .             ,g_def%zz(kgm ),arr(i,j,km )
     .             ,vzp,method)
        endif

        if ((k == 1).and.(.not.bc_safe)) then
          qkm=inter(g_def%zz(kgp ),arr(i,j,kp )
     .             ,g_def%zz(kg  ),arr(i,j,k  )
     .             ,g_def%zz(kgm ),arr(i,j,km )
     .             ,g_def%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,bc_adv_meth)
        else
          qkm=inter(g_def%zz(kgp ),arr(i,j,kp )
     .             ,g_def%zz(kg  ),arr(i,j,k  )
     .             ,g_def%zz(kgm ),arr(i,j,km )
     .             ,g_def%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,method)
        endif

        !Setup fluxes
        if (.not.isBdry(g_def,i,igx,2).or.bc_safe) flxip=vxp*qip
        if (.not.isBdry(g_def,i,igx,1).or.bc_safe) flxim=vxm*qim
                                            
        if (.not.isBdry(g_def,j,igy,4).or.bc_safe) flxjp=vyp*qjp
        if (.not.isBdry(g_def,j,igy,3).or.bc_safe) flxjm=vym*qjm
                                            
        if (.not.isBdry(g_def,k,igz,6).or.bc_safe) flxkp=vzp*qkp
        if (.not.isBdry(g_def,k,igz,5).or.bc_safe) flxkm=vzm*qkm

      endif

c     Find advective term

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac
      
      if (vw) advec = advec*g_def%gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function flx_nc_advec

c     scl_advc_flx_x_new
c     ###############################################################
      function scl_advc_flx_x_new(g_def,i,j,k,nx,ny,nz,igr,v,arr,darr
     .                           ,method,sp,zip_vel) result(flxip)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_1 of scalar arr at face
c     (i+1/2,j,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr (0:,0:,0:)
     .          ,darr(0:,0:,0:)
     .          ,v   (0:,0:,0:),flxip

      logical,optional :: sp,zip_vel

c     Local variables

      integer :: ip,igp,ig,jg,kg

      real(8) :: jacip,jacim,jac,qip,vxp

      logical :: spoint,zv,at_bdry

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

c     Update grid indices

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      ip  = i+1
      igp = ig + (ip-i)

c     Grid info

      jac = g_def%gmetric%grid(igr)%jac(i,j,k)

      if (spoint) then
        jacip = g_def%gmetric%grid(igr)%jac(ip,j,k)
      else
        jacip = jac
      endif

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = (isBdry(g_def,i ,igr,2).and.bcond(2)/=PER)
     .      .or.(isBdry(g_def,ip,igr,1).and.bcond(1)/=PER)

c$$$      if ((v(ip,j,k) == -v(i,j,k)).and.at_bdry) then      !Impenetrable-wall bdry
c$$$        flxip = 0d0
c$$$      elseif(isSP2(g_def,i+1,igr,ibc=1).and.spoint) then  !SP boundary
      if(isSP2(g_def,i+1,igr,ibc=1).and.spoint) then !SP boundary
        flxip = 0d0
      elseif(at_bdry) then                                !Other boundaries
        flxip = 0.5*(v(ip,j,k)*arr(i,j,k)+v(i,j,k)*arr(ip,j,k))
      else                                                !Not at boundary
        flxip = 0.25*(v(ip,j,k)*arr(i ,j,k)/jacip
     .              + v(i ,j,k)*arr(ip,j,k)/jac  )*(jac+jacip)
      endif

c     Other advective methods

      if (method /= 2.and.(.not.at_bdry)) then

        !Find velocities at faces from ZIP fluxes
        if (zv
     .     .and.(arr(ip,j,k)*arr(i,j,k)>0d0)
     .     .and.(abs(arr(ip,j,k)+arr(i,j,k))>zip_threshold)) then
          vxp =2.*flxip/(arr(ip,j,k)+arr(i,j,k))
        else
          vxp = 0.5*(v(ip,j,k)+v(i,j,k))
        endif

        qip=in2face_scl_new(g_def%xx(ig ),arr(i ,j,k)
     .                     ,g_def%xx(igp),arr(ip,j,k)
     .                     ,2*g_def%dxh(ig ),darr(i ,j,k)
     .                     ,2*g_def%dxh(igp),darr(ip,j,k)
     .                     ,vxp,method)

        !Setup fluxes (ZIP at boundaries by default for conservation)
        flxip=vxp*qip

      endif

c     End program

      end function scl_advc_flx_x_new

c     scl_advc_flx_y_new
c     ###############################################################
      function scl_advc_flx_y_new(g_def,i,j,k,nx,ny,nz,igr,v,arr,darr
     .                           ,method,zip_vel) result(flxjp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_2 of scalar arr at face
c     (i,j+1/2,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr (0:,0:,0:)
     .          ,darr(0:,0:,0:)
     .          ,v   (0:,0:,0:),flxjp

      logical,optional :: zip_vel

c     Local variables

      integer :: jp,jgp,ig,jg,kg

      real(8) :: jac,qjp,vyp

      logical :: zv,at_bdry

c     Begin program

      if (PRESENT(zip_vel)) then
        zv = zip_vel
      else
        zv = .true.
      endif

c     Update grid indices

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      jp  = j+1
      jgp = jg + (jp-j)

c     Grid info

      jac = g_def%gmetric%grid(igr)%jac(i,j,k)

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = (isBdry(g_def,j ,igr,4).and.bcond(4)/=PER)
     .      .or.(isBdry(g_def,jp,igr,3).and.bcond(3)/=PER)

c$$$      if ((v(i,j,k) == -v(i,jp,k)).and.at_bdry) then
c$$$        flxjp = 0d0             !Shut off flux for vn=0 BC
c$$$      else
        flxjp = 0.5*(v(i,jp,k)*arr(i,j,k) + v(i,j,k)*arr(i,jp,k))
c$$$      endif

c     Other advective methods

      if (method /= 2.and.(.not.at_bdry)) then

        !Find velocities at faces from ZIP fluxes
        if (zv
     .     .and.(    arr(i,jp,k)*arr(i,j,k) >0d0)
     .     .and.(abs(arr(i,jp,k)+arr(i,j,k))>zip_threshold)) then
          vyp =2.*flxjp/(arr(i,jp,k)+arr(i,j,k))
        else
          vyp = 0.5*(v(i,jp,k)+v(i,j,k))
        endif

        qjp=in2face_scl_new(g_def%yy(jg ),arr(i,j ,k)
     .                     ,g_def%yy(jgp),arr(i,jp,k)
     .                     ,2*g_def%dyh(jg ),darr(i,j ,k)
     .                     ,2*g_def%dyh(jgp),darr(i,jp,k)
     .                     ,vyp,method)

        !Setup fluxes (ZIP at boundaries by default for conservation)
        flxjp=vyp*qjp

      endif

c     End program

      end function scl_advc_flx_y_new

c     scl_advc_flx_z_new
c     ###############################################################
      function scl_advc_flx_z_new(g_def,i,j,k,nx,ny,nz,igr,v,arr,darr
     .                           ,method,zip_vel) result(flxkp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_3 of scalar arr at face
c     (i,j,k+1/2) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr (0:,0:,0:)
     .          ,darr(0:,0:,0:)
     .          ,v   (0:,0:,0:),flxkp

      logical,optional :: zip_vel

c     Local variables

      integer :: kp,kgp,ig,jg,kg

      real(8) :: jac,qkp,vzp

      logical :: zv,at_bdry

c     Begin program

      if (PRESENT(zip_vel)) then
        zv = zip_vel
      else
        zv = .true.
      endif

c     Update grid indices

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      kp  = k+1
      kgp = kg + (kp-k)

c     Grid info

      jac = g_def%gmetric%grid(igr)%jac(i,j,k)

c     Find X ZIP flux (used even if method /= ZIP)

      at_bdry = (isBdry(g_def,k ,igr,6).and.bcond(6)/=PER)
     .      .or.(isBdry(g_def,kp,igr,5).and.bcond(5)/=PER)

c$$$      if ((v(i,j,k) == -v(i,j,kp)).and.at_bdry) then
c$$$        flxkp = 0d0             !Shut off flux for vn=0 BC
c$$$      else
        flxkp = 0.5*(v(i,j,kp)*arr(i,j,k) + v(i,j,k)*arr(i,j,kp))
c$$$      endif

c     Other advective methods

      if (method /= 2.and.(.not.at_bdry)) then

        !Find velocities at faces from ZIP fluxes
        if (zv
     .      .and.(arr(i,j,kp)*arr(i,j,k)>0d0)
     .      .and.abs(arr(i,j,kp)+arr(i,j,k))>zip_threshold) then
          vzp =2.*flxkp/(arr(i,j,kp)+arr(i,j,k))
        else
          vzp = 0.5*(v(i,j,kp)+v(i,j,k))
        endif

        !Find advected quantity at faces
        qkp=in2face_scl_new(g_def%zz(kg ),arr(i,j,k )
     .                     ,g_def%zz(kgp),arr(i,j,kp)
     .                     ,2*g_def%dzh(kg ),darr(i,j,k )
     .                     ,2*g_def%dzh(kgp),darr(i,j,kp)
     .                     ,vzp,method)

        !Setup fluxes (ZIP at boundaries by default for conservation)
        flxkp=vzp*qkp

      endif

c     End program

      end function scl_advc_flx_z_new

c     flx_advec_mesh
c     ###############################################################
      function flx_advec_mesh(g_def,igrid,v,arr,method,bcs,sp,vol,no_jac
     .                       ,zip_vel) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid,method,bcs(6)

      real(8) :: arr(0:,0:,0:)
     .          ,v  (0:,0:,0:,:)
     $          ,advec(0:size(v,1)-1
     $                ,0:size(v,2)-1
     $                ,0:size(v,3)-1)

      logical,optional :: sp,vol,no_jac,zip_vel

c     Local variables

      integer :: i,j,k,nx,ny,nz,ig,jg,kg,bcnd(6,3),ivar

      real(8) :: dxx,dyy,dzz
      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      real(8),pointer,dimension(:,:,:,:) :: darr
      logical :: vw,njac,zv

c     Begin program

      nx = size(arr,1)-2
      ny = size(arr,2)-2
      nz = size(arr,3)-2

      allocate(darr(0:nx+1,0:ny+1,0:nz+1,3))

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .false.
      endif

      if (PRESENT(no_jac)) then
        njac = no_jac
      else
        njac = .false.
      endif

c     Compute difference arrays and communicate BCs

      darr(1:nx,:,:,1) = arr(2:nx+1,:,:)-arr(0:nx-1,:,:)
      darr(:,1:ny,:,2) = arr(:,2:ny+1,:)-arr(:,0:ny-1,:)
      darr(:,:,1:nz,3) = arr(:,:,2:nz+1)-arr(:,:,0:nz-1)

      bcnd(:,1) = abs(bcs)
c$$$      bcnd(:,1) = (/bcs(1),bcs(2),DEF,DEF,DEF,DEF/)
      !Change parity (darr is gradient)
      where     (bcnd(1:2,1) == NEU)
        bcnd(1:2,1) = DIR
      elsewhere (bcnd(1:2,1) == DIR)
        bcnd(1:2,1) = NEU
      end where
      
      bcnd(:,2) = abs(bcs)
c$$$      bcnd(:,2) = (/DEF,DEF,bcs(3),bcs(4),DEF,DEF/)
      !Change parity (darr is gradient)
      where     (bcnd(3:4,2) == DIR)
        bcnd(3:4,2) = NEU
      elsewhere (bcnd(3:4,2) == NEU)
        bcnd(3:4,2) = DIR
      end where

      bcnd(:,3) = abs(bcs)
c$$$      bcnd(:,3) = (/DEF,DEF,DEF,DEF,bcs(5),bcs(6)/)
      !Change parity (darr is gradient)
      where     (bcnd(5:6,3) == DIR)
        bcnd(5:6,3) = NEU
      elsewhere (bcnd(5:6,3) == NEU)
        bcnd(5:6,3) = DIR
      end where

      do ivar=1,3
        call setMGBC(g_def,0,1,nx,ny,nz,igrid,darr(:,:,:,ivar)
     $              ,bcnd(:,ivar),icomp=(/1/),iorder=1)
      enddo

c$$$      open(unit=110,file='debug.bin',form='unformatted'
c$$$     .    ,status='replace')
c$$$      call contour(darr(:,:,1,1),nx+2,ny+2,0d0,1d0,0d0,1d0,0,110)
c$$$      call contour(darr(:,:,1,2),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
c$$$      call contour(darr(:,:,1,3),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
c$$$      close(110)
      
c     Compute advection

      advec = 0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx
            
            !Grid info
            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            dxx = g_def%dxh(ig)
            dyy = g_def%dyh(jg)
            dzz = g_def%dzh(kg)

            !Find fluxes
            flxip = scl_advc_flx_x_new(g_def,i  ,j,k,nx,ny,nz,igrid
     .                                ,v(:,:,:,1),arr,darr(:,:,:,1)
     .                                ,method,sp=sp,zip_vel=zip_vel)
            flxim = scl_advc_flx_x_new(g_def,i-1,j,k,nx,ny,nz,igrid
     .                                ,v(:,:,:,1),arr,darr(:,:,:,1)
     .                                ,method,sp=sp,zip_vel=zip_vel)

            flxjp = scl_advc_flx_y_new(g_def,i,j  ,k,nx,ny,nz,igrid
     .                                ,v(:,:,:,2),arr,darr(:,:,:,2)
     .                                ,method,zip_vel=zip_vel)
            flxjm = scl_advc_flx_y_new(g_def,i,j-1,k,nx,ny,nz,igrid
     .                                ,v(:,:,:,2),arr,darr(:,:,:,2)
     .                                ,method,zip_vel=zip_vel)

            flxkp = scl_advc_flx_z_new(g_def,i,j,k  ,nx,ny,nz,igrid
     .                                ,v(:,:,:,3),arr,darr(:,:,:,3)
     .                                ,method,zip_vel=zip_vel)
            flxkm = scl_advc_flx_z_new(g_def,i,j,k-1,nx,ny,nz,igrid
     .                                ,v(:,:,:,3),arr,darr(:,:,:,3)
     .                                ,method,zip_vel=zip_vel)

            !Find advective term

            advec(i,j,k) =( (flxip - flxim)/dxx
     .                    + (flxjp - flxjm)/dyy
     .                    + (flxkp - flxkm)/dzz )

c$$$            !Old advective code
c$$$            advec(i,j,k) = flx_advec(g_def,i,j,k,nx,ny,nz,igrid,igrid
c$$$     .                            ,igrid,v,arr,method,zip_vel=zv,vol=vw)
c$$$     .                    *g_def%gmetric%grid(igrid)%jac(i,j,k)

          enddo
        enddo
      enddo

      if (vw) then
        advec = advec*g_def%gmetric%grid(igrid)%dlvol
      elseif (.not.njac) then
        advec = advec/g_def%gmetric%grid(igrid)%jac
      endif

      deallocate(darr)

c     End 

      end function flx_advec_mesh

c     scl_nc_advc_flx_x_new
c     ###############################################################
      function scl_nc_advc_flx_x_new(g_def,i,j,k,nx,ny,nz,igr,vxp,arr
     .                              ,darr,method) result(flxip)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_1 of scalar arr at face
c     (i+1/2,j,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr (0:nx+1,0:ny+1,0:nz+1)
     .          ,darr(0:nx+1,0:ny+1,0:nz+1)
     .          ,vxp,flxip

c     Local variables

      integer :: ip,igp,ig,jg,kg,meth

      real(8) :: qip

      logical :: at_bdry

c     Begin program

c     Update grid indices

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      ip  = i+1
      igp = ig + (ip-i)

c     Figure out advective method

      at_bdry = (isBdry(g_def,i  ,igr,2).and.bcond(2)/=PER)
     .      .or.(isBdry(g_def,i+1,igr,1).and.bcond(1)/=PER)

      if (at_bdry.and.method > ZIP) then
         meth = UPWD  !Revert to upwinding
      else
         meth = method
      endif

c     Find advected quantity at faces

      qip=in2face_scl_new(g_def%xx(ig ),arr(i ,j,k)
     .                   ,g_def%xx(igp),arr(ip,j,k)
     .                   ,2*g_def%dxh(ig ),darr(i ,j,k)
     .                   ,2*g_def%dxh(igp),darr(ip,j,k)
     .                   ,vxp,meth)

      !Setup fluxes
      flxip=vxp*qip

c     End program

      end function scl_nc_advc_flx_x_new

c     scl_nc_advc_flx_y_new
c     ###############################################################
      function scl_nc_advc_flx_y_new(g_def,i,j,k,nx,ny,nz,igr,vyp,arr
     .                              ,darr,method) result(flxjp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_2 of scalar arr at face
c     (i,j+1/2,k) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr (0:nx+1,0:ny+1,0:nz+1)
     .          ,darr(0:nx+1,0:ny+1,0:nz+1)
     .          ,vyp,flxjp

c     Local variables

      integer :: jp,jgp,ig,jg,kg,meth

      real(8) :: qjp

      logical :: at_bdry

c     Begin program

c     Update grid indices

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      jp  = j+1
      jgp = jg + (jp-j)

c     Figure out advective method

      at_bdry = (isBdry(g_def,j  ,igr,4).and.bcond(4)/=PER)
     .      .or.(isBdry(g_def,j+1,igr,3).and.bcond(3)/=PER)

      if (at_bdry.and.method > ZIP) then
         meth = UPWD  !Revert to upwinding
      else
         meth = method
      endif

c     Find advected quantity at faces

      qjp=in2face_scl_new(g_def%yy(jg ),arr(i,j ,k)
     .                   ,g_def%yy(jgp),arr(i,jp,k)
     .                   ,2*g_def%dyh(jg ),darr(i,j ,k)
     .                   ,2*g_def%dyh(jgp),darr(i,jp,k)
     .                   ,vyp,method)

      !Setup fluxes
      flxjp=vyp*qjp

c     End program

      end function scl_nc_advc_flx_y_new

c     scl_nc_advc_flx_z_new
c     ###############################################################
      function scl_nc_advc_flx_z_new(g_def,i,j,k,nx,ny,nz,igr,vzp,arr
     .                              ,darr,method) result(flxkp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates advective flux (v.arr)_3 of scalar arr at face
c     (i,j,k+1/2) in general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr,method

      real(8) :: arr (0:nx+1,0:ny+1,0:nz+1)
     .          ,darr(0:nx+1,0:ny+1,0:nz+1)
     .          ,vzp,flxkp

c     Local variables

      integer :: kp,kgp,ig,jg,kg,meth

      real(8) :: qkp

      logical :: at_bdry

c     Begin program

c     Update grid indices

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      kp  = k+1
      kgp = kg + (kp-k)

c     Figure out advective method

      at_bdry = (isBdry(g_def,k  ,igr,6).and.bcond(6)/=PER)
     .      .or.(isBdry(g_def,k+1,igr,5).and.bcond(5)/=PER)

      if (at_bdry.and.method > ZIP) then
         meth = UPWD  !Revert to upwinding
      else
         meth = method
      endif

c     Find advected quantity at faces

      qkp=in2face_scl_new(g_def%zz(kg ),arr(i,j,k )
     .                   ,g_def%zz(kgp),arr(i,j,kp)
     .                   ,2*g_def%dzh(kg ),darr(i,j,k )
     .                   ,2*g_def%dzh(kgp),darr(i,j,kp)
     .                   ,vzp,method)

      !Setup fluxes
      flxkp=vzp*qkp

c     End program

      end function scl_nc_advc_flx_z_new

c     flx_nc_advec_mesh
c     ###############################################################
      function flx_nc_advec_mesh(g_def,igrid,v,arr,method,bcs,vol)
     $         result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid,method,bcs(6)

      real(8) :: arr(0:,0:,0:)
     .          ,v  (0:,0:,0:,:)
     $          ,advec(0:size(v,1)-1
     $                ,0:size(v,2)-1
     $                ,0:size(v,3)-1)

      logical,optional :: vol

c     Local variables

      integer :: i,j,k,nx,ny,nz,ig,jg,kg,bcnd(6,3),ivar

      real(8) :: dxx,dyy,dzz
      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      real(8),pointer,dimension(:,:,:,:) :: darr
      logical :: vw

c     Debug

      real(8),allocatable,dimension(:,:,:) :: debug2
      character(13) :: debugfile

c     Begin program

      nx = size(arr,1)-2
      ny = size(arr,2)-2
      nz = size(arr,3)-2

      allocate(darr(0:nx+1,0:ny+1,0:nz+1,3))

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .false.
      endif

c     Compute difference arrays and communicate BCs

      darr = 0d0

      darr(1:nx,:,:,1) = arr(2:nx+1,:,:)-arr(0:nx-1,:,:)
      darr(:,1:ny,:,2) = arr(:,2:ny+1,:)-arr(:,0:ny-1,:)
      darr(:,:,1:nz,3) = arr(:,:,2:nz+1)-arr(:,:,0:nz-1)

c$$$      bcnd(:,1) = bcs
c$$$c$$$      bcnd(:,1) = (/bcs(1),bcs(2),DEF,DEF,DEF,DEF/)
c$$$c$$$      bcnd(:,1) = (/bcond(1),bcond(2),DEF,DEF,DEF,DEF/)
c$$$      !Change parity (darr is gradient)
c$$$      where (bcnd(1:2,1) == NEU) bcnd(1:2,1) = DIR
c$$$      where (bcnd(1:2,1) == DIR) bcnd(1:2,1) = NEU
c$$$
c$$$      bcnd(:,2) = bcs
c$$$c$$$      bcnd(:,2) = (/DEF,DEF,bcs(3),bcs(4),DEF,DEF/)
c$$$c$$$      bcnd(:,2) = (/DEF,DEF,bcond(3),bcond(4),DEF,DEF/)
c$$$      !Change parity (darr is gradient)
c$$$      where (bcnd(3:4,2) == DIR) bcnd(3:4,2) = NEU
c$$$      where (bcnd(3:4,2) == NEU) bcnd(3:4,2) = DIR
c$$$
c$$$      bcnd(:,3) = bcs
c$$$c$$$      bcnd(:,3) = (/DEF,DEF,DEF,DEF,bcs(5),bcs(6)/)
c$$$c$$$      bcnd(:,3) = (/DEF,DEF,DEF,DEF,bcond(5),bcond(6)/)
c$$$      !Change parity (darr is gradient)
c$$$      where (bcnd(5:6,3) == DIR) bcnd(5:6,3) = NEU
c$$$      where (bcnd(5:6,3) == NEU) bcnd(5:6,3) = DIR
c$$$
c$$$      do ivar=1,3
c$$$        call setMGBC(0,1,nx,ny,nz,igrid,darr(:,:,:,ivar)
c$$$     $              ,bcnd(:,ivar),icomp=(/1/),iorder=1)
c$$$      enddo

      !Setup BCs (do nothing for known BC values)
      bcnd(:,1) = (/bcs(1),bcs(2),DEF,DEF,DEF,DEF/)
      bcnd(:,2) = (/DEF,DEF,bcs(3),bcs(4),DEF,DEF/)
      bcnd(:,3) = (/DEF,DEF,DEF,DEF,bcs(5),bcs(6)/)

      !Check for special boundaries
      where (bcond == SP ) bcnd(:,1) = SP
      where (bcond == SP ) bcnd(:,2) = SP
      where (bcond == SP ) bcnd(:,3) = SP
      where (bcond == PER) bcnd(:,1) = PER
      where (bcond == PER) bcnd(:,2) = PER
      where (bcond == PER) bcnd(:,3) = PER
cc      where (bcnd /= DEF .and. bcnd /= PER .and. bcnd /= SP) bcnd = EXT
      where (bcnd /= PER .and. bcnd /= SP) bcnd = EXT

      call setMGBC(g_def,0,3,nx,ny,nz,igrid,darr,bcnd,icomp=(/1/)
     .            ,is_cnv=.false.,is_vec=.true.,iorder=1)
      
c diag *****
c$$$      if (my_rank == 0) then
c$$$        write (*,*) 'comp 1',darr(:,1,1,1)
c$$$        write (*,*) 'comp 2',darr(:,1,1,1)
c$$$        write (*,*) 'comp 3',darr(:,1,1,1)
c$$$      endif
c$$$
c$$$      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
c$$$      open(unit=110,file=debugfile,form='unformatted'
c$$$     .    ,status='replace')
c$$$
c$$$      k   = 1
c$$$
c$$$      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
c$$$
c$$$      debug2 = arr
c$$$      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2,0d0
c$$$     .     ,xmax,0d0,ymax,0,110)
c$$$      debug2 = darr(:,:,:,1)
c$$$      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2,0d0
c$$$     .     ,xmax,0d0,ymax,1,110)
c$$$      debug2 = darr(:,:,:,2)
c$$$      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2,0d0
c$$$     .     ,xmax,0d0,ymax,1,110)
c$$$      debug2 = darr(:,:,:,3)
c$$$      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2,0d0
c$$$     .     ,xmax,0d0,ymax,1,110)
c$$$
c$$$      deallocate(debug2)
c$$$      close(110)
c$$$
c$$$cc      call pstop('flx_nc_advec_mesh','DIAG')
c$$$cc      stop
c$$$cc      pause
c diag *****
c     Compute advection

      advec =0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx
            
            !Grid info
            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            dxx = g_def%dxh(ig)
            dyy = g_def%dyh(jg)
            dzz = g_def%dzh(kg)

            !Find fluxes
            flxip = scl_nc_advc_flx_x_new(g_def,i  ,j,k,nx,ny,nz,igrid
     $                                   ,v(i,j,k,1),arr,darr(:,:,:,1)
     .                                   ,method)
            flxim = scl_nc_advc_flx_x_new(g_def,i-1,j,k,nx,ny,nz,igrid
     $                                   ,v(i,j,k,1),arr,darr(:,:,:,1)
     .                                   ,method)

            flxjp = scl_nc_advc_flx_y_new(g_def,i,j  ,k,nx,ny,nz,igrid
     $                                   ,v(i,j,k,2),arr,darr(:,:,:,2)
     .                                   ,method)
            flxjm = scl_nc_advc_flx_y_new(g_def,i,j-1,k,nx,ny,nz,igrid
     $                                   ,v(i,j,k,2),arr,darr(:,:,:,2)
     .                                   ,method)

            flxkp = scl_nc_advc_flx_z_new(g_def,i,j,k  ,nx,ny,nz,igrid
     $                                   ,v(i,j,k,3),arr,darr(:,:,:,3)
     .                                   ,method)
            flxkm = scl_nc_advc_flx_z_new(g_def,i,j,k-1,nx,ny,nz,igrid
     $                                   ,v(i,j,k,3),arr,darr(:,:,:,3)
     .                                   ,method)

            !Find advective term

            advec(i,j,k) =( (flxip - flxim)/dxx
     .                    + (flxjp - flxjm)/dyy
     .                    + (flxkp - flxkm)/dzz )

          enddo
        enddo
      enddo

      if (vw) then
        advec = advec*g_def%gmetric%grid(igrid)%dlvol
      else
        advec = advec/g_def%gmetric%grid(igrid)%jac
      endif

      deallocate(darr)

c     End 

      end function flx_nc_advec_mesh

c     tensor_nc_advc
c     #####################################################################
      function tensor_nc_advc(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,v1
     .                       ,meth,transpose,v1_cov)  result(vec)

c     ---------------------------------------------------------------------
c     Calculates v0.nabla(v1), with both v0, v1 vector fields.  Vector
c     v0 is always contravariant; v1 can be covariant or
c     contravariant. Result is in same representation as v1.
c
c     meth: specifies the advective scheme (only upwind=1 and central=2)
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz,meth
      real(8) :: vec(3)
      real(8),target :: v1(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .                 ,v0(0:nxx+1,0:nyy+1,0:nzz+1,3)

      logical, optional :: transpose,v1_cov

c     Local variables

      integer    :: ieq,hex,hey,hez
      real(8)    :: nabla_v(3,3),jac

      logical :: transp,v1cov

c     Begin program

      if (PRESENT(transpose)) then
        transp = transpose
      else
        transp = .false.
      endif

      if (PRESENT(v1_cov)) then
        v1cov = v1_cov
      else
        v1cov = .false.
      endif

      jac = g_def%gmetric%grid(igx)%jac(i,j,k)

      select case(meth)
      case(1) !Upwinding
        hex = int(sign(1d0,-v0(i,j,k,1)))
        hey = int(sign(1d0,-v0(i,j,k,2)))
        hez = int(sign(1d0,-v0(i,j,k,3)))

        if (v1cov) then
          nabla_v = fnabla_v_upwd_cov(g_def,i,j,k,nxx,nyy,nzz
     .                           ,igx,igy,igz
     .                           ,v1(:,:,:,1)
     .                           ,v1(:,:,:,2)
     .                           ,v1(:,:,:,3),hex,hey,hez)
        else
          nabla_v = fnabla_v_upwd(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                           ,v1(:,:,:,1)
     .                           ,v1(:,:,:,2)
     .                           ,v1(:,:,:,3),hex,hey,hez)
        endif
      case(2) !Central differences
        if (v1cov) then
          nabla_v = fnabla_v_cov(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                          ,v1(:,:,:,1)
     .                          ,v1(:,:,:,2)
     .                          ,v1(:,:,:,3),0)
        else
          nabla_v = fnabla_v(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,v1,0)
        endif
      case default
        call pstop('tensor_nc_advc','Advection option not implemented')
      end select

      if (transp) then
        vec = matmul(nabla_v,v0(i,j,k,:))/jac
      else
        vec = matmul(v0(i,j,k,:),nabla_v)/jac
      endif

      end function tensor_nc_advc

c     tensor_nc_advc_mesh2
c     #####################################################################
      function tensor_nc_advc_mesh2(g_def,igr,v0,v1,meth,v1_cov
     .                             ,transpose)
     .         result(vec)

c     ---------------------------------------------------------------------
c     Calculates v0.nabla(v1), with both v0, v1 vector fields, on the whole
c     mesh. Vector v0 is always contravariant; v1 can be covariant or
c     contravariant. Result is in same representation as v1.
c
c     meth: specifies the advective scheme
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr,meth
      real(8) :: v1(0:,0:,0:,:)
     .          ,v0(0:,0:,0:,:)
      real(8) :: vec(0:size(v1,1)-1
     .              ,0:size(v1,2)-1
     .              ,0:size(v1,3)-1,size(v1,4))

      logical :: v1_cov
      logical,optional :: transpose
      
c     Local variables

      integer :: i,j,k,nx,ny,nz
    
c     Begin program

      nx = size(v1,1)-2
      ny = size(v1,2)-2
      nz = size(v1,3)-2

      vec = 0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx
            vec(i,j,k,:)=tensor_nc_advc(g_def,i,j,k,nx,ny,nz
     .                                 ,igr,igr,igr,v0,v1,meth
     .                                 ,transpose=transpose
     .                                 ,v1_cov=v1_cov)
          enddo
        enddo
      enddo

c     End program

      end function tensor_nc_advc_mesh2

c     tensor_nc_advc_mesh
c     #####################################################################
      function tensor_nc_advc_mesh(g_def,igr,v0,v1,bcs,meth,v1_cov)
     .         result(vec)

c     ---------------------------------------------------------------------
c     Calculates v0.nabla(v1), with both v0, v1 vector fields, on the whole
c     mesh. Vector v0 is always contravariant; v1 can be covariant or
c     contravariant. Result is in same representation as v1.
c
c     meth: specifies the advective scheme
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr,meth,bcs(6,3)
      real(8) :: v1(0:,0:,0:,:)
     .          ,v0(0:,0:,0:,:)
      real(8) :: vec(0:size(v1,1)-1
     .              ,0:size(v1,2)-1
     .              ,0:size(v1,3)-1,size(v1,4))

      logical :: v1_cov
      
c     Local variables

      integer :: i,j,k,nx,ny,nz
      real(8) :: jac,tensor(3,3)

c     Begin program

      nx = size(v1,1)-2
      ny = size(v1,2)-2
      nz = size(v1,3)-2

      !Perform advective step
      vec(:,:,:,1)
     .     = flx_nc_advec_mesh(g_def,igr,v0,v1(:,:,:,1),meth,bcs(:,1))
      vec(:,:,:,2)
     .     = flx_nc_advec_mesh(g_def,igr,v0,v1(:,:,:,2),meth,bcs(:,2))
      vec(:,:,:,3)
     .     = flx_nc_advec_mesh(g_def,igr,v0,v1(:,:,:,3),meth,bcs(:,3))

      !Add geometric source
      if (v1_cov) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              jac = g_def%gmetric%grid(igr)%jac(i,j,k)
              tensor = nabla_v_src_cov(i,j,k,v1(i,j,k,:))
              vec(i,j,k,:)=vec(i,j,k,:)+matmul(v0(i,j,k,:),tensor)/jac
            enddo
          enddo
        enddo
      else
        do k=1,nz
          do j=1,ny
            do i=1,nx
              jac = g_def%gmetric%grid(igr)%jac(i,j,k)

              tensor = nabla_v_src(i,j,k,v1(i,j,k,:))
              vec(i,j,k,:)=vec(i,j,k,:)+matmul(v0(i,j,k,:),tensor)/jac
            enddo
          enddo
        enddo
      endif

c     End program

      contains

c     nabla_v_src
c     #############################################################
      function nabla_v_src(i,j,k,cnv) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k
        real(8)    :: tensor(3,3),cnv(3)

c     Local variables

        integer    :: l,m
        real(8)    :: dum(3)

c     Begin program

        dum = g_def%gmetric%grid(igr)%Gamma(i,j,k,1,1,:)
     .       +g_def%gmetric%grid(igr)%Gamma(i,j,k,2,2,:)
     .       +g_def%gmetric%grid(igr)%Gamma(i,j,k,3,3,:)

        do l=1,3
          do m=1,3
            tensor(l,m) =
     .       dot_product(cnv,g_def%gmetric%grid(igr)%Gamma(i,j,k,m,:,l))
     .      -cnv(m)*dum(l)
          enddo
        enddo

c     End program

      end function nabla_v_src

c     nabla_v_src_cov
c     #############################################################
      function nabla_v_src_cov(i,j,k,cov) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

        integer    :: i,j,k
        real(8)    :: tensor(3,3),cov(3)

        integer    :: l,m
        real(8)    :: hessian(3,3,3)

c     Begin program

        hessian = g_def%gmetric%grid(igr)%Gamma(i,j,k,:,:,:)

        do l=1,3
          do m=1,3
            tensor(l,m) = -dot_product(cov,hessian(:,l,m))
          enddo
        enddo

c     End program

      end function nabla_v_src_cov

      end function tensor_nc_advc_mesh

c     tensor_advc
c     #####################################################################
      function tensor_advc(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,v1
     .                    ,method) result(cnv)

c     ---------------------------------------------------------------------
c     Calculates divergence of tensor, div(v0 v1), at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme (see int2face)
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz,method
      real(8) :: cnv(3)
      real(8),target :: v1(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .                 ,v0(0:nxx+1,0:nyy+1,0:nzz+1,3)

c     Local variables

c     Begin program

      adv_meth = method

      vec1 => v0
      vec2 => v1
      cnv = div_tensor(g_def,i,j,k,nxx,nyy,nzz,igx,igy,igz,alt__eom()
     .                ,tnsr_advc_x,tnsr_advc_y,tnsr_advc_z)
      nullify(vec1,vec2)

      end function tensor_advc

c     tnsr_advc_x
c     #############################################################
      subroutine tnsr_advc_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                      ,t11,t12,t13,flag)
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

        type(grid_mg_def),pointer :: g_def

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

        jac = 0.5*(g_def%gmetric%grid(igx)%jac(ip,j,k)
     .            +g_def%gmetric%grid(igx)%jac(i ,j,k))

        if ( i + g_def%ilo(igx)-1 < g_def%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = g_def%gmetric%grid(igx)%jac(ip,j,k)
          jac0 = g_def%gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (isSP(g_def,i+1,j,k,igx,igy,igz)) jac = SP_flsv

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

          call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

          ipp  = min(ip+1,nx+1)
          im   = max(i -1,0   )

          igp  = ig + (ip  - i)
          igm  = ig + (im  - i)
          igpp = ig + (ipp - i)

          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac

          !Limit flux based on velocity
          t1 = vxp*inter(g_def%xx(igpp),vec2(ipp,j,k,:)
     .                  ,g_def%xx(igp ),vec2(ip ,j,k,:)
     .                  ,g_def%xx(ig  ),vec2(i  ,j,k,:)
     .                  ,g_def%xx(igm ),vec2(im ,j,k,:)
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
      subroutine tnsr_advc_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                      ,t21,t22,t23,flag)
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

        type(grid_mg_def),pointer :: g_def

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

        jac = 0.5*(g_def%gmetric%grid(igx)%jac(i,jp,k)
     .            +g_def%gmetric%grid(igx)%jac(i,j ,k))

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

          call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

          jpp = min(jp+1,ny+1)
          jm  = max(j -1,0   )

          jgp  = jg + (jp  - j)
          jgm  = jg + (jm  - j)
          jgpp = jg + (jpp - j)

          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))

          !Limit flux based on velocity
          t2 = vyp*inter(g_def%yy(jgpp),vec2(i,jpp,k,:)
     .                  ,g_def%yy(jgp ),vec2(i,jp ,k,:)
     .                  ,g_def%yy(jg  ),vec2(i,j  ,k,:)
     .                  ,g_def%yy(jgm ),vec2(i,jm ,k,:)
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
      subroutine tnsr_advc_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                      ,t31,t32,t33,flag)
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

        type(grid_mg_def),pointer :: g_def

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

        jac= 0.5*(g_def%gmetric%grid(igx)%jac(i,j,kp)
     .           +g_def%gmetric%grid(igx)%jac(i,j,k ))

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

          call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

          kpp = min(kp+1,nz+1)
          km  = max(k -1,0   )

          kgp  = kg + (kp  - k)
          kgm  = kg + (km  - k)
          kgpp = kg + (kpp - k)

          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))

            !Limit flux based on velocity
          t3=vzp*inter(g_def%zz(kgpp),vec2(i,j,kpp,:)
     .                ,g_def%zz(kgp ),vec2(i,j,kp ,:)
     .                ,g_def%zz(kg  ),vec2(i,j,k  ,:)
     .                ,g_def%zz(kgm ),vec2(i,j,km ,:)
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

c     tensor_advc_mesh
c     #####################################################################
      function tensor_advc_mesh(g_def,igr,v0,v1,bcs,meth) result(vec)

c     ---------------------------------------------------------------------
c     Calculates v0.nabla(v1), with both v0, v1 vector fields, on the whole
c     mesh. Vector v0 is always contravariant; v1 can be covariant or
c     contravariant. Result is in same representation as v1.
c
c     meth: specifies the advective scheme
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr,meth,bcs(6,3)
      real(8),target :: v1(0:,0:,0:,:)
     .                 ,v0(0:,0:,0:,:)
      real(8) :: vec(0:size(v1,1)-1
     .              ,0:size(v1,2)-1
     .              ,0:size(v1,3)-1,size(v1,4))

c     Local variables

      integer :: i,j,k,nx,ny,nz
      real(8) :: msource(3)
     
c     Begin program

      nx = size(v1,1)-2
      ny = size(v1,2)-2
      nz = size(v1,3)-2

      !Perform advective step
      vec(:,:,:,1) = v1(:,:,:,1)/g_def%gmetric%grid(igr)%jac
      if (.not.alt__eom()) then
        vec(:,:,:,2) = v1(:,:,:,2)/g_def%gmetric%grid(igr)%jac
      else
        vec(:,:,:,2) = v1(:,:,:,2)
      endif
      vec(:,:,:,3) = v1(:,:,:,3)/g_def%gmetric%grid(igr)%jac

      vec(:,:,:,1) = flx_advec_mesh(g_def,igr,v0,vec(:,:,:,1),meth
     .                          ,bcs(:,1),no_jac=.true.,zip_vel=.false.)
      vec(:,:,:,2) = flx_advec_mesh(g_def,igr,v0,vec(:,:,:,2),meth
     .                          ,bcs(:,2),no_jac=.true.,zip_vel=.false.)
      vec(:,:,:,3) = flx_advec_mesh(g_def,igr,v0,vec(:,:,:,3),meth
     .                          ,bcs(:,3),no_jac=.true.,zip_vel=.false.)

      !Add geometric source
      vec1 => v0
      vec2 => v1
      if (.not.alt__eom()) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              msource =dtnsr_src_ijk(g_def,i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,tnsr_advc_x,tnsr_advc_y,tnsr_advc_z)
              vec(i,j,k,:) = vec(i,j,k,:) + msource
            enddo
          enddo
        enddo
      else
        do k=1,nz
          do j=1,ny
            do i=1,nx
              msource =dtnsr_src_ijk(g_def,i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,tnsr_advc_x,tnsr_advc_y,tnsr_advc_z)
              vec(i,j,k,:) = vec(i,j,k,:) + msource
              vec(i,j,k,2) = vec(i,j,k,2)
     .                      /g_def%gmetric%grid(igr)%jac(i,j,k)
            enddo
          enddo
        enddo
      endif

c     End program

      end function tensor_advc_mesh

c     in2face_scl
c     ###############################################################
      function in2face_scl(x4,q4,x3,q3,x2,q2,x1,q1,vel,advect)
     .         result(inter)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation at control
c     volume face based of face velocity "vel".
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
     $             ,dq1,dq2,dq3,dx1,dx2,dx3

c     Begin program

      dq1 = q2-q1
      dq2 = q3-q2
      dq3 = q4-q3

      dx1 = x2-x1
      dx2 = x3-x2
      dx3 = x4-x3

      select case(advect)
      case(UPWD)

c     Upwind

        if (vel.gt.0d0) then
          inter = q2
        else
          inter = q3
        endif

      case(ZIP)

c     Centered

        inter = 0.5*(q3+q2)

      case(QUICK)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
cc          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
        else
cc          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
        endif

cc        inter = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
        inter = 0.5*(q3+q2) - 0.125*curv*(dx2)**2

      case(SMART)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = smart_slp1
        slp2 = smart_slp2
        if (vel.gt.0d0) then
c$$$          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
c$$$          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q2 + (1.-slp1)*q1
c$$$          qt3 = slp2*q2 + (1.-slp2)*q3
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q2 + (slp1-1.)*dq1
          qt3 = q2 - (slp2-1.)*dq2
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
c$$$          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
c$$$          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q3 + (1.-slp1)*q4
c$$$          qt3 = slp2*q3 + (1.-slp2)*q2
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q3 - (slp1-1.)*dq3
          qt3 = q3 + (slp2-1.)*dq2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(SSMART)

c     Smooth SMART (uniform mesh only)

        xp1 = ssmart_xp1
        xp2 = ssmart_xp2

        if (vel.gt.0d0) then
cc          qt1  = (q2-q1)/(q3-q1)
          qt1  = dq1/(dq2+dq1)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1d0)
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elseif (qt1.ge.xp2.and.qt1.lt.1d0) then
            a = (0.25d0*xp2-.5d0)/(xp2-1d0)**3
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
cc            inter = 0.375d0*q3 + 0.75d0*q2 - 0.125d0*q1
            inter = (3*q3 + 5*q2 + dq1)/8.
          else
            inter = q2
          endif
        else
cc          qt1  = (q3-q4)/(q2-q4)
          qt1  = dq3/(dq3+dq2)
          if (qt1.lt.xp1.and.qt1.ge.0d0) then
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            c = 1d0
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq2+dq3)
          elseif (qt1.ge.xp2.and.qt1.lt.1d0) then
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq3+dq2)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
cc            inter = 0.375d0*q2 + 0.75d0*q3- 0.125d0*q4
            inter = (3*q2 + 5*q3 - dq3)/8.
          else
            inter = q3
          endif
        endif

      case(CHO)

c     Centered high-order

cc        inter = 0.0625d0*(-q1 + 9.*q2 + 9.*q3 - q4)
        inter = 0.0625d0*(dq1-dq3) + 0.5*(q2 + q3)

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

c     in2face_scl_new
c     ###############################################################
      function in2face_scl_new(x2,q2,x3,q3,dcx2,dcq2,dcx3,dcq3,vel,advc)
     .         result(inter)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation at control
c     volume face based of face velocity "vel".
c
c     Convention:
c                  dcq3,dcx3
c                /-----^-----\
c            dcq2,dcx2
c          /-----^-----\
c       ---x-----x--o--x-----x---->
c                q2 ^  q3
c                   ^
c           Location of face
c     ---------------------------------------------------------------

      use oned_int

      implicit none                !For safe fortran

c     Call variables

      real(8)    :: x2,x3,q2,q3,dcx2,dcq2,dcx3,dcq3,vel,inter
      integer    :: advc

c     Local variables

      real(8)    :: qt1,qt2,qt3,qt4,slp1,slp2,a,b,c,xp1,xp2,curv
     $             ,dq1,dq2,dq3,dx1,dx2,dx3

c     Begin program

      dx2 = x3-x2
      dx1 = dcx2 - dx2
      dx3 = dcx3 - dx2

      if (dx3 == 0d0) dx3 = dx2
      if (dx1 == 0d0) dx1 = dx2

      dq2 = q3-q2
      dq1 = dcq2 - dq2
      dq3 = dcq3 - dq2
      
      select case(advc)
      case(UPWD)

c     Upwind

        if (vel.gt.0d0) then
          inter = q2
        else
          inter = q3
        endif

      case(ZIP)

c     Centered

        inter = 0.5*(q3+q2)

      case(QUICK)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
cc          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
        else
cc          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
        endif

cc        inter = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
        inter = 0.5*(q3+q2) - 0.125*curv*(dx2)**2

      case(SMART)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = smart_slp1
        slp2 = smart_slp2
        if (vel.gt.0d0) then
c$$$          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
c$$$          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q2 + (1.-slp1)*q1
c$$$          qt3 = slp2*q2 + (1.-slp2)*q3
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q2 + (slp1-1.)*dq1
          qt3 = q2 - (slp2-1.)*dq2
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
c$$$          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
c$$$          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q3 + (1.-slp1)*q4
c$$$          qt3 = slp2*q3 + (1.-slp2)*q2
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q3 - (slp1-1.)*dq3
          qt3 = q3 + (slp2-1.)*dq2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(SSMART)

c     Smooth SMART (uniform mesh only)

        xp1 = ssmart_xp1
        xp2 = ssmart_xp2

        if (vel.gt.0d0) then
cc          qt1  = (q2-q1)/(q3-q1)
          qt1  = dq1/(dq2+dq1)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1d0)
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elseif (qt1.ge.xp2.and.qt1.lt.1d0) then
            a = (0.25d0*xp2-.5d0)/(xp2-1d0)**3
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
cc            inter = 0.375d0*q3 + 0.75d0*q2 - 0.125d0*q1
            inter = (3*q3 + 5*q2 + dq1)/8.
          else
            inter = q2
          endif
        else
cc          qt1  = (q3-q4)/(q2-q4)
          qt1  = dq3/(dq3+dq2)
          if (qt1.lt.xp1.and.qt1.ge.0d0) then
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            c = 1d0
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq2+dq3)
          elseif (qt1.ge.xp2.and.qt1.lt.1d0) then
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq3+dq2)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
cc            inter = 0.375d0*q2 + 0.75d0*q3- 0.125d0*q4
            inter = (3*q2 + 5*q3 - dq3)/8.
          else
            inter = q3
          endif
        endif

      case(CHO)

c     Centered high-order

cc        inter = 0.0625d0*(-q1 + 9.*q2 + 9.*q3 - q4)
        inter = 0.0625d0*(dq1-dq3) + 0.5*(q2 + q3)

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

        call pstop('in2face_scl_new','Advection scheme not available')

      end select

c     End

      end function in2face_scl_new

c     in2face_vec
c     ###############################################################
      function in2face_vec(x4,q4,x3,q3,x2,q2,x1,q1,vel,advect)
     .         result(inter)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation of vector
c     quantity at control volume face based of face velocity "vel".
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
     $          ,dq1(size(q1)),dq2(size(q1)),dq3(size(q1))
     .          ,dx1,dx2,dx3

c     Begin program

      dq1 = q2-q1
      dq2 = q3-q2
      dq3 = q4-q3

      dx1 = x2-x1
      dx2 = x3-x2
      dx3 = x4-x3

      select case(advect)
      case(UPWD)

c     Upwind

        if (vel.gt.0d0) then
          inter = q2
        else
          inter = q3
        endif

      case(ZIP)

c     Centered

        inter = 0.5*(q3+q2)

      case(QUICK)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
cc          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
        else
cc          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
        endif

        inter = 0.5*(q3+q2) - 0.125*curv*(dx2)**2

      case(SMART)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = smart_slp1
        slp2 = smart_slp2

        if (vel.gt.0d0) then
c$$$          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
c$$$          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q2 + (1.-slp1)*q1
c$$$          qt3 = slp2*q2 + (1.-slp2)*q3
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q2 + (slp1-1.)*dq1
          qt3 = q2 - (slp2-1.)*dq2
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
c$$$          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
c$$$          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q3 + (1.-slp1)*q4
c$$$          qt3 = slp2*q3 + (1.-slp2)*q2
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q3 - (slp1-1.)*dq3
          qt3 = q3 + (slp2-1.)*dq2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(SSMART)

c     Smooth SMART

        xp1 = ssmart_xp1
        xp2 = ssmart_xp2

        if (vel.gt.0d0) then
cc          qt1  = (q2-q1)/(q3-q1)
          qt1  = dq1/(dq2+dq1)
          where (qt1.lt.xp1.and.qt1.ge.0d0)
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1d0)
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1d0)
            a = (0.25d0*xp2-.5d0)/(xp2-1d0)**3
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
cc            inter = 0.375d0*q3 + 0.75d0*q2 - 0.125d0*q1
            inter = (3*q3 + 5*q2 + dq1)/8.
          elsewhere
            inter = q2
          end where
        else
cc          qt1  = (q3-q4)/(q2-q4)
          qt1  = dq3/(dq3+dq2)
          where (qt1.lt.xp1.and.qt1.ge.0d0)
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            c = 1d0
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq2+dq3)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1d0)
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq3+dq2)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
cc            inter = 0.375d0*q2 + 0.75d0*q3- 0.125d0*q4
            inter = (3*q2 + 5*q3 - dq3)/8.
          elsewhere
            inter = q3
          end where
        endif

      case(CHO)

c     Centered high-order

cc        inter = 0.0625d0*(-q1 + 9.*q2 + 9.*q3 - q4)
        inter = 0.0625d0*(dq1-dq3) + 0.5*(q2 + q3)

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

c     in2face_vec_new
c     ###############################################################
      function in2face_vec_new(q2,q3,dx1,dq1,dx2,dq2,dx3,dq3,vel,advect)
     .         result(inter)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation at control
c     volume face based of face velocity "vel".
c
c     Convention:
c
c            dq1   dq2   dq3
c       ---x-----x--o--x-----x---->
c                q2 ^  q3
c                   ^
c           Location of face
c     ---------------------------------------------------------------

      use oned_int

      implicit none                !For safe fortran

c     Call variables

      real(8) :: q2(:),q3(:),vel
     $          ,dq1(:),dq2(:),dq3(:),dx1,dx2,dx3
      real(8) :: inter(size(q2))
      integer :: advect

c     Local variables

      real(8) :: qt1(size(q2)),qt2(size(q2))
     .          ,qt3(size(q2)),qt4(size(q2))
     .          ,slp1,slp2,a(size(q2)),b(size(q2)),c(size(q2))
     .          ,xp1,xp2,curv(size(q2))

c     Begin program

      select case(advect)
      case(UPWD)

c     Upwind

        if (vel.gt.0d0) then
          inter = q2
        else
          inter = q3
        endif

      case(ZIP)

c     Centered

        inter = 0.5*(q3+q2)

      case(QUICK)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
cc          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
        else
cc          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
        endif

        inter = 0.5*(q3+q2) - 0.125*curv*(dx2)**2

      case(SMART)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = smart_slp1
        slp2 = smart_slp2

        if (vel.gt.0d0) then
c$$$          curv = 2d0*((q3-q2)/(x3-x2)-(q2-q1)/(x2-x1))/(x3-x1)
c$$$          qt1 = 0.5*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q2 + (1.-slp1)*q1
c$$$          qt3 = slp2*q2 + (1.-slp2)*q3
          curv = 2d0*(dq2/dx2-dq1/dx1)/(dx2+dx1)
          qt1 = 0.5*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q2 + (slp1-1.)*dq1
          qt3 = q2 - (slp2-1.)*dq2
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        else
c$$$          curv = 2d0*((q4-q3)/(x4-x3)-(q3-q2)/(x3-x2))/(x4-x2)
c$$$          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(x3-x2)**2
c$$$          qt2 = slp1*q3 + (1.-slp1)*q4
c$$$          qt3 = slp2*q3 + (1.-slp2)*q2
          curv = 2d0*(dq3/dx3-dq2/dx2)/(dx3+dx2)
          qt1 = 0.5d0*(q3+q2) - 0.125*curv*(dx2)**2
          qt2 = q3 - (slp1-1.)*dq3
          qt3 = q3 + (slp2-1.)*dq2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        endif

      case(SSMART)

c     Smooth SMART

        xp1 = ssmart_xp1
        xp2 = ssmart_xp2

        if (vel.gt.0d0) then
cc          qt1  = (q2-q1)/(q3-q1)
          qt1  = dq1/(dq2+dq1)
          where (qt1.lt.xp1.and.qt1.ge.0d0)
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1d0)
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1d0)
            a = (0.25d0*xp2-.5d0)/(xp2-1d0)**3
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q1 + qt2*(q3-q1)
            inter = q3 + (qt2-1.)*(dq2+dq1)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
cc            inter = 0.375d0*q3 + 0.75d0*q2 - 0.125d0*q1
            inter = (3*q3 + 5*q2 + dq1)/8.
          elsewhere
            inter = q2
          end where
        else
cc          qt1  = (q3-q4)/(q2-q4)
          qt1  = dq3/(dq3+dq2)
          where (qt1.lt.xp1.and.qt1.ge.0d0)
            a = .25d0*(xp1-3d0)*xp1**(-3)
            b = (1.125d0-0.5d0*xp1)*xp1**(-2)
            c = 1d0
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq2+dq3)
          elsewhere (qt1.ge.xp2.and.qt1.lt.1d0)
            a = (0.25d0*xp2-.5d0)*(xp2-1d0)**(-3)
            b = 0.125d0*(-4d0*xp2**2+7d0*xp2+1d0)*(xp2-1d0)**(-3)
            c = 1d0 - (-.5d0*xp2**2-0.375d0+1.125d0*xp2)*(xp2-1d0)**(-3)
            qt2 = (qt1 - 1d0)*(a*qt1**2 + b*qt1 + c) + 1d0
cc            inter = q4 + qt2*(q2-q4)
            inter = q2 - (qt2-1.)*(dq3+dq2)
          elsewhere (qt1.ge.xp1.and.qt1.lt.xp2)
cc            inter = 0.375d0*q2 + 0.75d0*q3- 0.125d0*q4
            inter = (3*q2 + 5*q3 - dq3)/8.
          elsewhere
            inter = q3
          end where
        endif

      case(CHO)

c     Centered high-order

cc        inter = 0.0625d0*(-q1 + 9.*q2 + 9.*q3 - q4)
        inter = 0.0625d0*(dq1-dq3) + 0.5*(q2 + q3)

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

      end function in2face_vec_new

      end module grid_advec_ops

c grid_pardiff_ops
c ######################################################################
      module grid_pardiff_ops

        use grid_create

        INTERFACE par_diff
          module procedure par_diff_ijk,par_diff_mesh
        end INTERFACE

      contains

c     par_diff_mesh
c     ###############################################################
      function par_diff_mesh(g_def,igrid,arr,BB,par_order,method,vol
     .                      ,tmp0) result(lap)

c     ---------------------------------------------------------------
c     Calculates grad(A) in general non-orthogonal coordinates. Is
c     boundary-ready (can extrapolate). The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid,par_order

      real(8) :: arr(0:,0:,0:)
     .          ,lap(0:size(arr,1)-1,0:size(arr,2)-1,0:size(arr,3)-1)
      real(8) :: BB(0:,0:,0:,:)

      character(*) :: method

      logical,optional :: vol
      real(8),optional :: tmp0(0:size(arr,1)-1
     .                        ,0:size(arr,2)-1
     .                        ,0:size(arr,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = g_def%nxv(igrid)
      ny = g_def%nyv(igrid)
      nz = g_def%nzv(igrid)

      do k=1,nz
        do j=1,ny
          do i=1,nx
            lap(i,j,k)=par_diff_ijk(g_def,i,j,k,igrid,arr,BB,par_order
     .                             ,method,vol=vol,tmp0=tmp0)
          enddo
        enddo
      enddo

c     End program

      end function par_diff_mesh

c     par_diff_ijk
c     ################################################################
      function par_diff_ijk(g_def,i,j,k,igr,tmp,BB,par_order,method,vol
     .                     ,tmp0) result(par2)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 2nd or 4th order
c     accuracy. Here:
c       * par_order=2,4 (discretization order)
c       * method='symm','sharma','mimetic','central'
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer  :: i,j,k,igr,par_order
      real(8)  :: tmp(0:,0:,0:),par2
      real(8)  :: BB(0:,0:,0:,:)

      character(*) :: method

      logical,optional :: vol
      real(8),optional :: tmp0(0:size(tmp,1)-1
     .                        ,0:size(tmp,2)-1
     .                        ,0:size(tmp,3)-1)

c     Local variables

      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      integer  :: nx,ny,nz,ig,jg,kg
      logical  :: vol_wgt,limited

c     Begin program

      nx = g_def%nxv(igr)
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      Fluxx=0d0
      Fluxy=0d0
      Fluxz=0d0

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

c     Find fluxes

      select case(par_order)
      case(2)                   !second order
        call par_flux_2nd(g_def,i,j,k,nx,ny,nz,igr,igr,igr,tmp
     .                   ,BB(:,:,:,1),BB(:,:,:,2),BB(:,:,:,3)
     .                   ,Fluxx,Fluxy,Fluxz,method,tmp0=tmp0)
      case(4)                   !fourth order, 5-pt stencil in co-directed dervs.
        limited = (trim(method) == 'sharma')
        call par_flux_4th(g_def,i,j,k,nx,ny,nz,igr,igr,igr,tmp
     .                   ,BB(:,:,:,1),BB(:,:,:,2),BB(:,:,:,3)
     .                   ,Fluxx,Fluxy,Fluxz,limited=limited,stencil=5)
      case(-4) !fourth order, 7-pt stencil in co-directed dervs.
        limited = (trim(method) == 'sharma')
        call par_flux_4th(g_def,i,j,k,nx,ny,nz,igr,igr,igr,tmp
     .                   ,BB(:,:,:,1),BB(:,:,:,2),BB(:,:,:,3)
     .                   ,Fluxx,Fluxy,Fluxz,limited=limited,stencil=7)
      case default
        write (*,*) 'Order not implemented'
        stop
      end select

c     Assemble fluxes

      call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

      par2 = (Fluxx(i,j,k)-Fluxx(i-1,j,k))/g_def%dxh(ig)
     .      +(Fluxy(i,j,k)-Fluxy(i,j-1,k))/g_def%dyh(jg)
     .      +(Fluxz(i,j,k)-Fluxz(i,j,k-1))/g_def%dzh(kg)

      par2 = par2/g_def%gmetric%grid(igr)%jac(i,j,k)

      if (vol_wgt) par2=par2*g_def%gmetric%grid(igr)%dvol(i,j,k)

      end function par_diff_ijk

c     par_flux_2nd
c     ################################################################
      subroutine par_flux_2nd(g_def,i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .               ,Bx,By,Bz,Fluxx,Fluxy,Fluxz,method,tmp0)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 2nd order accuracy
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8),target  :: tmp(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      character(*) :: method

      real(8),optional,target :: tmp0(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      real(8)  :: bx2,by2,bz2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
     .           ,jac_avg
      integer  :: ig,jg,kg,ii,jj,kk
      logical  :: limit,symm

      real(8),dimension(:,:,:),pointer :: jac
      real(8),dimension(:,:,:,:,:),pointer :: gsup

      real(8),pointer,dimension(:,:,:) :: tmp_0

c     Begin program

      Fluxx=0d0
      Fluxy=0d0
      Fluxz=0d0

      jac  => g_def%gmetric%grid(igx)%jac
      gsup => g_def%gmetric%grid(igx)%gsup

      select case(trim(method))
      case("symm")  !Symmetric discretization

        do kk = k-1,k
          do jj = j-1,j
            do ii = i-1,i
              call getMGmap(g_def,ii,jj,kk,igx,igy,igz,ig,jg,kg)

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
     .                 /g_def%dx(ig)
              dty = .5*(tmp(ii  ,jj+1,kk)-tmp(ii  ,jj,kk)
     .                 +tmp(ii+1,jj+1,kk)-tmp(ii+1,jj,kk))
     .                 /g_def%dy(jg)

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

      case ("sharma") !Limited differences (Sharma, Hammett, JCP 227 (2007))
        
        !X-component cross flux
        do ii = i-1,i
          call getMGmap(g_def,ii,j,k,igx,igy,igz,ig,jg,kg)

          bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)
     .               *by(ii  ,j,k)
     .               +bx(ii+1,j,k)/jac(ii+1,j,k)
     .               *by(ii+1,j,k))

          flx0 = mc(tmp(ii  ,j  ,k)-tmp(ii  ,j-1,k)
     .             ,tmp(ii  ,j+1,k)-tmp(ii  ,j  ,k))
          flxp = mc(tmp(ii+1,j  ,k)-tmp(ii+1,j-1,k)
     .             ,tmp(ii+1,j+1,k)-tmp(ii+1,j  ,k))

          Fluxx(ii,j,k)= bxby/g_def%dyh(jg)*mc(flx0,flxp)
        enddo

        !Y-component cross flux 
        do jj = j-1,j
          call getMGmap(g_def,i,jj,k,igx,igy,igz,ig,jg,kg)

          bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)
     .               *by(i,jj  ,k)
     .               +bx(i,jj+1,k)/jac(i,jj+1,k)
     .               *by(i,jj+1,k))

          flx0 = mc(tmp(i  ,jj  ,k)-tmp(i-1,jj  ,k)
     .             ,tmp(i+1,jj  ,k)-tmp(i  ,jj  ,k))
          flxp = mc(tmp(i  ,jj+1,k)-tmp(i-1,jj+1,k)
     .             ,tmp(i+1,jj+1,k)-tmp(i  ,jj+1,k))

          Fluxy(i,jj,k)= bxby/g_def%dxh(ig)*mc(flx0,flxp)
        enddo

        !Add co-derivative terms: standard centered
        do ii = i-1,i
          if (isSP2(g_def,ii+1,igx)) cycle

          call getMGmap(g_def,ii,j,k,igx,igy,igz,ig,jg,kg)

cc          bx2 = (0.5*(bx(ii+1,j,k)/jac(ii+1,j,k)
cc     .               +bx(ii  ,j,k)/jac(ii  ,j,k)))**2

          jac_avg = 0.5*(jac(ii+1,j,k)+jac(ii,j,k))
          bx2 = (0.5*(bx(ii+1,j,k)+bx(ii  ,j,k)))**2/jac_avg

          Fluxx(ii,j,k) = Fluxx(ii,j,k)
     .                  + bx2*(tmp(ii+1,j,k)-tmp(ii,j,k))
     .                       /g_def%dx(ig)
        enddo

        do jj = j-1,j
          call getMGmap(g_def,i,jj,k,igx,igy,igz,ig,jg,kg)
 
cc          by2 = (0.5*(by(i,jj+1,k)/jac(i,jj+1,k)
cc     .               +by(i,jj  ,k)/jac(i,jj  ,k)))**2

          jac_avg = 0.5*(jac(i,jj+1,k)+jac(i,jj,k))
          by2 = (0.5*(by(i,jj+1,k)+by(i,jj,k)))**2/jac_avg

          Fluxy(i,jj,k) = Fluxy(i,jj,k)
     .                  + by2*(tmp(i,jj+1,k)-tmp(i,jj,k))
     .                       /g_def%dy(jg)
        enddo

        do kk = k-1,k
          call getMGmap(g_def,i,j,kk,igx,igy,igz,ig,jg,kg)

cc          bz2 = (0.5*(bz(i,j,kk+1)/jac(i,j,kk+1)
cc     .               +bz(i,j,kk  )/jac(i,j,kk  )))**2

          jac_avg = 0.5*(jac(i,j,kk+1)+jac(i,j,kk))
          bz2 = (0.5*(bz(i,j,kk+1)+bz(i,j,kk)))**2/jac_avg

          Fluxz(i,j,kk) = Fluxz(i,j,kk)
     .                  + bz2*(tmp(i,j,kk+1)-tmp(i,j,kk))
     .                       /g_def%dz(kg)
        enddo

      case ("mimetic") !Lypnikov et al. 2012 Russ J. Numer. Anal. Math. Modeling

        if (PRESENT(tmp0)) then
          tmp_0 => tmp0
        else
          tmp_0 => tmp
        endif

        do ii = i-1,i
          Fluxx(ii,j,k) = Mimetic_D_flux_x(g_def,ii,j,k,igx,igy,igz
     .                 ,bx(ii:ii+1,j,k)*bx(ii:ii+1,j,k)/jac(ii:ii+1,j,k)
     .                 ,bx(ii:ii+1,j,k)*by(ii:ii+1,j,k)/jac(ii:ii+1,j,k)
     .                 ,by(ii:ii+1,j,k)*by(ii:ii+1,j,k)/jac(ii:ii+1,j,k)
     .                 ,tmp,tmp_0)
        enddo

        do jj = j-1,j
          Fluxy(i,jj,k) = Mimetic_D_flux_y(g_def,i,jj,k,igx,igy,igz
     .                 ,bx(i,jj:jj+1,k)*bx(i,jj:jj+1,k)/jac(i,jj:jj+1,k)
     .                 ,bx(i,jj:jj+1,k)*by(i,jj:jj+1,k)/jac(i,jj:jj+1,k)
     .                 ,by(i,jj:jj+1,k)*by(i,jj:jj+1,k)/jac(i,jj:jj+1,k)
     .                 ,tmp,tmp_0)
        enddo

        do kk = k-1,k
          Fluxz(i,j,kk) = 0d0
        enddo

        nullify(tmp_0)

      case ("central") !Standard central average

        !X-component cross flux
        do ii = i-1,i
          call getMGmap(g_def,ii,j,k,igx,igy,igz,ig,jg,kg)

          bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)*by(ii  ,j,k)
     .               +bx(ii+1,j,k)/jac(ii+1,j,k)*by(ii+1,j,k))

          Fluxx(ii,j,k)= 0.25/g_def%dyh(jg)*bxby
     .           *(tmp(ii  ,j+1,k)-tmp(ii  ,j-1,k)
     .            +tmp(ii+1,j+1,k)-tmp(ii+1,j-1,k))
        enddo

        !Y-component cross flux 
        do jj = j-1,j
          call getMGmap(g_def,i,jj,k,igx,igy,igz,ig,jg,kg)

          bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)*by(i,jj  ,k)
     .               +bx(i,jj+1,k)/jac(i,jj+1,k)*by(i,jj+1,k))

          Fluxy(i,jj,k)= 0.25/g_def%dxh(ig)*bxby
     .                 *(tmp(i+1,jj  ,k)-tmp(i-1,jj  ,k)
     .                  +tmp(i+1,jj+1,k)-tmp(i-1,jj+1,k))
        enddo

        !Add co-derivative terms: standard centered
        do ii = i-1,i
          if (isSP2(g_def,ii+1,igx)) cycle

          call getMGmap(g_def,ii,j,k,igx,igy,igz,ig,jg,kg)

cc          bx2 = (0.5*(bx(ii+1,j,k)/jac(ii+1,j,k)
cc     .               +bx(ii  ,j,k)/jac(ii  ,j,k)))**2

          jac_avg = 0.5*(jac(ii+1,j,k)+jac(ii,j,k))
          bx2 = (0.5*(bx(ii+1,j,k)+bx(ii  ,j,k)))**2/jac_avg

          Fluxx(ii,j,k) = Fluxx(ii,j,k)
     .                  + bx2*(tmp(ii+1,j,k)-tmp(ii,j,k))
     .                       /g_def%dx(ig)
        enddo

        do jj = j-1,j
          call getMGmap(g_def,i,jj,k,igx,igy,igz,ig,jg,kg)
 
cc          by2 = (0.5*(by(i,jj+1,k)/jac(i,jj+1,k)
cc     .               +by(i,jj  ,k)/jac(i,jj  ,k)))**2

          jac_avg = 0.5*(jac(i,jj+1,k)+jac(i,jj,k))
          by2 = (0.5*(by(i,jj+1,k)+by(i,jj,k)))**2/jac_avg

          Fluxy(i,jj,k) = Fluxy(i,jj,k)
     .                  + by2*(tmp(i,jj+1,k)-tmp(i,jj,k))
     .                       /g_def%dy(jg)
        enddo

        do kk = k-1,k
          call getMGmap(g_def,i,j,kk,igx,igy,igz,ig,jg,kg)

cc          bz2 = (0.5*(bz(i,j,kk+1)/jac(i,j,kk+1)
cc     .               +bz(i,j,kk  )/jac(i,j,kk  )))**2

          jac_avg = 0.5*(jac(i,j,kk+1)+jac(i,j,kk))
          bz2 = (0.5*(bz(i,j,kk+1)+bz(i,j,kk)))**2/jac_avg

          Fluxz(i,j,kk) = Fluxz(i,j,kk)
     .                  + bz2*(tmp(i,j,kk+1)-tmp(i,j,kk))
     .                       /g_def%dz(kg)
        enddo

      case default

        call pstop("par_flux_2nd"
     .       ,"Wrong flux option; use symm,sharma,mimetic,central")

      end select

cc      !Add co-derivative terms: ZIP
cc      do i = 0,nx
cc        call getMGmap(g_def,i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz) 
cc     .              +bx(i+1,1:ny,1:nz)*bx(i,1:ny,1:nz)
cc     .                   *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cc      enddo
cc
cc      do j = 0,ny
cc        call getMGmap(g_def,1,j,1,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz) 
cc     .              +by(1:nx,j+1,1:nz)*by(1:nx,j,1:nz)
cc     .                   *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cc      enddo
cc
cc      do k = 0,nz
cc        call getMGmap(g_def,1,1,k,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
cc     .              +bz(1:nx,1:ny,k+1)*bz(1:nx,1:ny,k)
cc     .                   *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cc      enddo

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
cc              call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)
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
cc          call getMGmap(g_def,i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxx(i,1:ny,1:nz) =
cc     .         (tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cc        enddo
cc
cc        do j = 0,ny
cc          call getMGmap(g_def,1,j,1,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxy(1:nx,j,1:nz) =
cc     .         (tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cc        enddo
cc
cc        do k = 0,nz
cc          call getMGmap(g_def,1,1,k,igx,igy,igz,ig,jg,kg)
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
      subroutine par_flux_4th(g_def,i,j,k,nx,ny,nz,igx,igy,igz,F
     .                      ,Bx,By,Bz,Fluxx,Fluxy,Fluxz,limited,stencil)
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

      type(grid_mg_def),pointer :: g_def

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

      call getMGmap(g_def,1,1,1,igx,igy,igz,ig,jg,kg)

      dx=g_def%dxh(ig)
      dy=g_def%dyh(jg)
      dz=g_def%dzh(kg)

      jac => g_def%gmetric%grid(igx)%jac

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

      if (PRESENT(chk)) chk = (mc /= ho_flx)

      end function mc

c     G_covector_calc
c     ##################################################################
      function G_covector_calc(Dxx,Dyy,Dxy,nr,np) 
     .result(lvec)
c     ------------------------------------------------------------------
c     This function will compute the co-vector between the 
c     cell centered stress tensor and the cell face in which the
c     flux is being computed
c     ------------------------------------------------------------------

      implicit none
      
c     Call variables
      real(8) :: Dxx, Dyy, Dxy, nr, np
      real(8) :: lvec(2)

c     Local variables

c     Begin program

      lvec(1) = Dxx*nr + Dxy*np
      lvec(2) = Dxy*nr + Dyy*np

c     End program

      end function G_covector_calc

c     solver_2x2_local
c     ##################################################################
      function solver_2x2_local(A,b) result(x)
c     ------------------------------------------------------------------
c     This simple routine will perform an exact 2x2 solve. 
c     I'm doing this because I don't want to call Lapack since it's
c     a general subroutine and not optimal for a 2x2 (which gets called
c     awefully a lot for a 2V system in our app)
c
c     Inputs: 
c     A:      A 2x2 matrix that the user inputs
c     b:      The RHS that the user inputs
c     x:      The solution vector that gets outputted
c
c     ------------------------------------------------------------------

      implicit none

c     Call variables
      real(8) :: A(2,2)
      real(8) :: b(2), x(2)

c     Local variables
      real(8) :: A_inv(2,2), A_det

c     Begin program

      ! Compute the determinant
      A_det        = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      ! Compute the inverse of matrix
      A_inv(1,1)   = A(2,2)/A_det
      A_inv(1,2)   =-A(1,2)/A_det
      A_inv(2,1)   =-A(2,1)/A_det
      A_inv(2,2)   = A(1,1)/A_det
      ! Solve
      x(1)         = A_inv(1,1)*b(1) + A_inv(1,2)*b(2)
      x(2)         = A_inv(2,1)*b(1) + A_inv(2,2)*b(2)
      
c     End Program

      end function solver_2x2_local

c     Mimetic_D_flux_x
c     ##################################################################
      function Mimetic_D_flux_x(g_def,i,j,k,igx,igy,igz
     .                         ,Dxx,Dxy,Dyy,f_in,f_pc)
     .result(flux_x)
c     ------------------------------------------------------------------
c     This function will compute the X diffusion flux @ (i+1/2,j,k)
c     ------------------------------------------------------------------

      implicit none
      
c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8) :: f_in(0:,0:,0:), f_pc(0:,0:,0:)
      real(8) :: Dxx(i:i+1), Dxy(i:i+1), Dyy(i:i+1)
      integer :: i,j,k,igx,igy,igz,ig,jg,kg
      real(8) :: flux_x

c     Local variables

      real(8) :: A(2,2), b(2)

      real(8) :: lvec_ij    (2), lvec_ipj  (2)
      real(8) :: tvec_ij    (2), tvec_ipj  (2)
     .          ,tvec_ijpm  (2), tvec_ipjpm(2)

      real(8) :: alpha_ij   (2), alpha_ipj (2)
      real(8) :: mu_R          , mu_L
      real(8) :: flux_R        , flux_L
      real(8) :: flux_R_pc     , flux_L_pc

c     Begin program

      ! Get MGmap
      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      ! Compute covectors
      lvec_ij = G_covector_calc(Dxx(i  ), Dyy(i  ) 
     .                         ,Dxy(i  ), 1d0,0d0)
      lvec_ipj= G_covector_calc(Dxx(i+1), Dyy(i+1)
     .                         ,Dxy(i+1),-1d0,0d0)

      !===============================================================
      ! Below is for the homogeneous isotropic medium case
      !===============================================================
      ! Compute the displacement vector between cell centers
      ! Aligned diagonal term
      tvec_ij  = (/ g_def%dxh(ig), 0d0/)
      tvec_ipj = (/-g_def%dxh(ig), 0d0/)
      ! Off-diagonal term
      if (lvec_ij(2)  > 0d0) then
         tvec_ijpm  = (/0d0, g_def%dyh(jg)/)
      else
         tvec_ijpm  = (/0d0,-g_def%dyh(jg)/)
      endif
      if (lvec_ipj(2) > 0d0) then
         tvec_ipjpm = (/0d0, g_def%dyh(jg)/)
      else
         tvec_ipjpm = (/0d0,-g_def%dyh(jg)/)
      endif

      ! Compute the rotated diffusion coefficients
      ! For cell j  ,k
      A(1,1)   = tvec_ij(1); A(1,2) = tvec_ijpm(1)
      A(2,1)   = tvec_ij(2); A(2,2) = tvec_ijpm(2)
      b        = lvec_ij
      alpha_ij = solver_2x2_local(A,b)
      ! For cell j+1,k
      A(1,1)   = tvec_ipj(1); A(1,2) = tvec_ipjpm(1)
      A(2,1)   = tvec_ipj(2); A(2,2) = tvec_ipjpm(2)
      b        = lvec_ipj
      alpha_ipj= solver_2x2_local(A,b)

c      if (alpha_jk(1) <= 0d0 .or. alpha_jk(2) <= 0d0 .or.
c     .    alpha_jpk(1)<= 0d0 .or. alpha_jpk(2)<= 0d0) then
c         write(*,*) 'alpha_jk = ', alpha_jk, 'alpha_jpk = ', alpha_jpk
c      endif

      ! Construct the flux from two sides
      if (lvec_ij(2)  > 0d0) then
         flux_L   = - alpha_ij (2)*(f_in(i,j+1,k) - f_in(i,j,k) )
         flux_L_pc= - alpha_ij (2)*(f_pc(i,j+1,k) - f_pc(i,j,k) )
      else
         flux_L   = - alpha_ij (2)*(f_in(i,j-1,k) - f_in(i,j,k) )
         flux_L_pc= - alpha_ij (2)*(f_pc(i,j-1,k) - f_pc(i,j,k) )
      endif

      if (lvec_ipj(2) > 0d0) then
         flux_R   = - alpha_ipj(2)*(f_in(i+1,j+1,k) - f_in(i+1,j,k))
         flux_R_pc= - alpha_ipj(2)*(f_pc(i+1,j+1,k) - f_pc(i+1,j,k))
      else
         flux_R   = - alpha_ipj(2)*(f_in(i+1,j-1,k) - f_in(i+1,j,k))
         flux_R_pc= - alpha_ipj(2)*(f_pc(i+1,j-1,k) - f_pc(i+1,j,k))
      endif

      !===============================================================      
      ! Compute the nonlinear weights by checking the 
      ! conditions of flux_L*flux_R
      !===============================================================
      if (abs(flux_R_pc)+abs(flux_L_pc) == 0d0) then           ! Check if the solution is not unique
         mu_L = 0.5d0
         mu_R = 0.5d0
      else                                                     ! Check product
         mu_L = abs(flux_R_pc)/(abs(flux_L_pc)+abs(flux_R_pc))
         mu_R = abs(flux_L_pc)/(abs(flux_L_pc)+abs(flux_R_pc))
      endif

      ! Assemble the flux
      flux_x  =  mu_L*flux_L     -mu_R*flux_R 
     .         +(mu_L*alpha_ij(1)+mu_R*alpha_ipj(1))
     .         *(f_in(i,j,k)-f_in(i+1,j,k))

      flux_x =-flux_x

c     End program

      end function Mimetic_D_flux_x

c     Mimetic_D_flux_y
c     ##################################################################
      function Mimetic_D_flux_y(g_def,i,j,k,igx,igy,igz
     .                         ,Dxx,Dxy,Dyy,f_in,f_pc)
     .result(flux_y)
c     ------------------------------------------------------------------
c     This function will compute the Y diffusion flux @ (i,j+1/2,k)
c     ------------------------------------------------------------------

      implicit none
      
c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8) :: f_in(0:,0:,0:), f_pc(0:,0:,0:)
      real(8) :: Dxx(j:j+1), Dxy(j:j+1), Dyy(j:j+1)
      integer :: i,j,k,igx,igy,igz,ig,jg,kg
      real(8) :: flux_y

c     Local variables

      real(8) :: A(2,2), b(2)

      real(8) :: lvec_ij    (2), lvec_ijp  (2)
      real(8) :: tvec_ij    (2), tvec_ijp  (2)
     .          ,tvec_ipmj  (2), tvec_ipmjp(2)

      real(8) :: alpha_ij   (2), alpha_ijp (2)
      real(8) :: mu_T          , mu_B
      real(8) :: flux_T        , flux_B
      real(8) :: flux_T_pc     , flux_B_pc

c     Begin program

      ! Get MGmap
      call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

      ! Compute covectors
      lvec_ij = G_covector_calc(Dxx(j)  ,Dyy(j)  ,Dxy(j)  , 0d0, 1d0)
      lvec_ijp= G_covector_calc(Dxx(j+1),Dyy(j+1),Dxy(j+1), 0d0,-1d0)
      
      !===============================================================
      ! Below is for the homogeneous isotropic medium case
      !===============================================================
      ! Compute the displacement vector between cell centers
      ! Aligned diagonal term
      tvec_ij  = (/ 0d0, g_def%dyh(jg)/)
      tvec_ijp = (/ 0d0,-g_def%dyh(jg)/)
      ! Off-diagonal term
      if (lvec_ij(1)  > 0d0) then
         tvec_ipmj  = (/ g_def%dxh(ig), 0d0/)
      else
         tvec_ipmj  = (/-g_def%dxh(ig), 0d0/)
      endif
      if (lvec_ijp(1) > 0d0) then
         tvec_ipmjp = (/ g_def%dxh(ig), 0d0/)
      else
         tvec_ipmjp = (/-g_def%dxh(ig), 0d0/)
      endif

      ! Compute the rotated diffusion coefficients
      ! For cell j  ,k
      A(1,1)   = tvec_ij(1); A(1,2) = tvec_ipmj(1)
      A(2,1)   = tvec_ij(2); A(2,2) = tvec_ipmj(2)
      b        = lvec_ij
      alpha_ij = solver_2x2_local(A,b)
      ! For cell j  ,k+1
      A(1,1)   = tvec_ijp(1); A(1,2) = tvec_ipmjp(1)
      A(2,1)   = tvec_ijp(2); A(2,2) = tvec_ipmjp(2)
      b        = lvec_ijp
      alpha_ijp= solver_2x2_local(A,b)

c$$$      if (alpha_ij(1) <= 0d0 .or. alpha_ij(2) <= 0d0 .or.
c$$$     .    alpha_ijp(1)<= 0d0 .or. alpha_ijp(2)<= 0d0) then
c$$$         write(*,*) 'j,k=',j,k,'alpha_ij = ', alpha_ij
c$$$     ., 'alpha_ijp = ', alpha_ijp
c$$$     .,'lvec_ijp = ',lvec_ijp, 'Drr_ijp=',Dxx(j,k+1)
c$$$     .,'Dpp_ijp =',Dyy(j,k+1), 'Drp_ijp = ',Dxy(j,k+1)
c$$$      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! I don't understand why I have to use alpha_ij(2) 
      ! and alpha_ijp(2) instead of alpha_ij(1) and alpha_ijp(1).
      ! If I don't, I get a completely wrong flux. Must understand
      ! why...
      !!!!!!!!!!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Construct the flux from two sides
      if (lvec_ij(1)  > 0d0) then
         flux_B   = - alpha_ij (2)*(f_in(i+1,j,k) - f_in(i,j,k))
         flux_B_pc= - alpha_ij (2)*(f_pc(i+1,j,k) - f_pc(i,j,k))
      else
         flux_B   = - alpha_ij (2)*(f_in(i-1,j,k) - f_in(i,j,k))
         flux_B_pc= - alpha_ij (2)*(f_pc(i-1,j,k) - f_pc(i,j,k))
      endif

      if (lvec_ijp(1) > 0d0) then
         flux_T   = - alpha_ijp(2)*(f_in(i+1,j+1,k) - f_in(i,j+1,k))
         flux_T_pc= - alpha_ijp(2)*(f_pc(i+1,j+1,k) - f_pc(i,j+1,k))
      else
         flux_T   = - alpha_ijp(2)*(f_in(i-1,j+1,k) - f_in(i,j+1,k))
         flux_T_pc= - alpha_ijp(2)*(f_pc(i-1,j+1,k) - f_pc(i,j+1,k))
      endif

      !===============================================================      
      ! Compute the nonlinear weights by checking the 
      ! conditions of flux_L*flux_R
      !===============================================================
      if (abs(flux_T_pc)+abs(flux_B_pc) == 0d0) then           ! Check if the solution is not unique
         mu_B = 0.5d0;  mu_T = 0.5d0
      else                                               ! Check product
         mu_B = abs(flux_T_pc)/(abs(flux_B_pc)+abs(flux_T_pc))
         mu_T = abs(flux_B_pc)/(abs(flux_B_pc)+abs(flux_T_pc))
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! I don't understand why I have to use alpha_ij(1) 
      ! and alpha_ijp(1) instead of alpha_ij(2) and alpha_ijp(2).
      ! If I don't, I get a completely wrong flux. Must understand
      ! why...
      !!!!!!!!!!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Assemble the flux
      flux_y =  mu_B*flux_B     -mu_T*flux_T 
     .        +(mu_B*alpha_ij(1)+mu_T*alpha_ijp(1))
     .        *(f_in(i,j,k)-f_in(i,j+1,k))

      flux_y =-flux_y

c      write(*,*) 'j,k=',j,k,'flux_vp = ', flux_vp

c     End program

      end function Mimetic_D_flux_y

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
      function integral(g_def,igr,array,average,logical_dom)

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on domain (nx)x(ny)x(nz).
c     ---------------------------------------------------------------

      implicit none

c    Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: array(0:,0:,0:),integral
      logical,optional :: average,logical_dom

c     Local variables

      integer :: i,j,k,nx,ny,nz

      real(8) :: volume,dvol,send_buf(2),rec_buf(2)

      logical :: avg,log_dom

c     Begin program

      nx = size(array,1)-2
      ny = size(array,2)-2
      nz = size(array,3)-2

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
              dvol = g_def%gmetric%grid(igr)%dlvol(i,j,k)

              if (isSYM(g_def,i,igr,dim=1,loc=0)
     .        .or.isSYM(g_def,i,igr,dim=1,loc=1)
     .        .or.isSYM(g_def,j,igr,dim=2,loc=0)
     .        .or.isSYM(g_def,j,igr,dim=2,loc=1)
     .        .or.isSYM(g_def,k,igr,dim=3,loc=0)
     .        .or.isSYM(g_def,k,igr,dim=3,loc=1)) dvol = 0.5*dvol

              integral = integral + array(i,j,k)*dvol
              !Averages in logical volume don't make sense
cc              dvol = g_def%gmetric%grid(igx)%dvol(i,j,k)
              volume = volume + dvol
            enddo
          enddo
        enddo
      else
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              dvol = g_def%gmetric%grid(igr)%dvol(i,j,k)

              if (isSYM(g_def,i,igr,dim=1,loc=0)
     .        .or.isSYM(g_def,i,igr,dim=1,loc=1)
     .        .or.isSYM(g_def,j,igr,dim=2,loc=0)
     .        .or.isSYM(g_def,j,igr,dim=2,loc=1)
     .        .or.isSYM(g_def,k,igr,dim=3,loc=0)
     .        .or.isSYM(g_def,k,igr,dim=3,loc=1)) dvol = 0.5*dvol

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

c     integral_stg
c     ################################################################
      function integral_stg(g_def,igr,array,average,logical_dom)
     .         result(integral)

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on domain (nx)x(ny)x(nz).
c     ---------------------------------------------------------------

      implicit none

c    Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: array(0:,0:,0:),integral
      logical,optional :: average,logical_dom

c     Local variables

      integer :: i,j,k,nx,ny,nz

      real(8) :: volume,dvol,send_buf(2),rec_buf(2)

      logical :: avg,log_dom

c     Begin program

      nx = size(array,1)-2
      ny = size(array,2)-2
      nz = size(array,3)-2

      avg = .false.
      if (PRESENT(average)) avg = average

      log_dom = .false.
      if (PRESENT(logical_dom)) log_dom = logical_dom

c     Integrate

      integral = 0d0
      volume   = 0d0

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            if (log_dom) then
              dvol = g_def%gmetric%grid(igr)%dlvol(i,j,k)
            else
              dvol = g_def%gmetric%grid(igr)%dvol (i,j,k)
            endif

            if (isSYM(g_def,i,igr,dim=1,loc=0)
     .      .or.isSYM(g_def,i,igr,dim=1,loc=1)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=0)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=1)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=0)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=1)) dvol = 0.5*dvol

            integral = integral + array(i,j,k)*dvol
            !Averages in logical volume don't make sense
            volume = volume + dvol
          enddo
        enddo
      enddo

      if (isBdry(g_def,1,igr,1).and.(.not.bcPER(1))) then
        i = 0
        do k = 1,nz
          do j = 1,ny
            if (log_dom) then
              dvol = g_def%gmetric%grid(igr)%dlvol(i,j,k)
            else
              dvol = g_def%gmetric%grid(igr)%dvol (i,j,k)
            endif

            if (isSYM(g_def,i,igr,dim=1,loc=0)
     .      .or.isSYM(g_def,i,igr,dim=1,loc=1)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=0)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=1)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=0)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=1)) dvol = 0.5*dvol

            integral = integral + array(i,j,k)*dvol
            !Averages in logical volume don't make sense
            volume = volume + dvol
          enddo
        enddo
      endif

      if (isBdry(g_def,1,igr,3).and.(.not.bcPER(2))) then
        j = 0
        do k = 1,nz
          do i = 1,nx
            if (log_dom) then
              dvol = g_def%gmetric%grid(igr)%dlvol(i,j,k)
            else
              dvol = g_def%gmetric%grid(igr)%dvol (i,j,k)
            endif

            if (isSYM(g_def,i,igr,dim=1,loc=0)
     .      .or.isSYM(g_def,i,igr,dim=1,loc=1)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=0)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=1)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=0)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=1)) dvol = 0.5*dvol

            integral = integral + array(i,j,k)*dvol
            !Averages in logical volume don't make sense
            volume = volume + dvol
          enddo
        enddo
      endif

      if (isBdry(g_def,1,igr,5).and.(.not.bcPER(3))) then
        k = 0
        do j = 1,ny
          do i = 1,nx
            if (log_dom) then
              dvol = g_def%gmetric%grid(igr)%dlvol(i,j,k)
            else
              dvol = g_def%gmetric%grid(igr)%dvol (i,j,k)
            endif

            if (isSYM(g_def,i,igr,dim=1,loc=0)
     .      .or.isSYM(g_def,i,igr,dim=1,loc=1)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=0)
     .      .or.isSYM(g_def,j,igr,dim=2,loc=1)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=0)
     .      .or.isSYM(g_def,k,igr,dim=3,loc=1)) dvol = 0.5*dvol

            integral = integral + array(i,j,k)*dvol
            !Averages in logical volume don't make sense
            volume = volume + dvol
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

      end function integral_stg

      end module grid_diff_ops


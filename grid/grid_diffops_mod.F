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

c     curlcurl_vrtx_ndiff
c     ###############################################################
      function curlcurl_vrtx_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
     .         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: crl(2,2,2,3),cov(3)

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

c     Curl at vertices

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

            crl(il,jl,kl,:) = cov
          enddo
        enddo
      enddo

c     Curl at cell center

      cnv = curl_vrtx(i,j,k,igx,igy,igz,crl,.false.)

      if (vol_wgt) cnv=cnv*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function curlcurl_vrtx_ndiff

c     curlcurl_vrtx_diff
c     ###############################################################
      function curlcurl_vrtx_diff(i,j,k,nx,ny,nz,igx,igy,igz,aa,diff
     .                           ,vol) result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(D*curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)
     .             ,diff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: crl(2,2,2,3),cov(3),dff

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

c     Curl at vertices

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

            dff = 0.125*sum(diff(ii:ii+1,jj:jj+1,kk:kk+1))

            crl(il,jl,kl,:) = cov*dff
          enddo
        enddo
      enddo

c     Curl at cell center

      cnv = curl_vrtx(i,j,k,igx,igy,igz,crl,.false.)

      if (vol_wgt) cnv=cnv*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function curlcurl_vrtx_diff

ccc     div_vrtx
ccc     ##################################################################
cc      function div_vrtx(i,j,k,igx,igy,igz,vx,vy,vz,vrtx) result(div)
cc
ccc     ------------------------------------------------------------------
ccc     Calculates div of cnv vector defined at nearest neighbors.
ccc     If vrtx is true, it finds div at (i+1/2,j+1/2,k+1/2) from (i,j,k)
ccc     If vrtx is false, it finds div at (i,j,k) from (i+1/2,j+1/2,k+1/2)
ccc     ------------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: i,j,k,igx,igy,igz
cc      real(8) :: vx(2,2,2),vy(2,2,2),vz(2,2,2),div
cc      logical :: vrtx
cc
ccc     Local variables
cc
cc      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0
cc      real(8) :: dhx,dhy,dhz,fp,fm,jac
cc
ccc     Begin program
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      if (vrtx) then
cc        dhx = dx(ig)
cc        dhy = dy(jg)
cc        dhz = dz(kg)
cc      else
cc        dhx = dxh(ig)
cc        dhy = dyh(jg)
cc        dhz = dzh(kg)
cc      endif
cc
cc      ip = 2
cc      i0 = 1
cc      jp = 2
cc      j0 = 1
cc      kp = 2
cc      k0 = 1
cc
cc      !X component
cc      fp = 0.25*(vx(ip,j0,k0)+vx(ip,jp,k0)
cc     .          +vx(ip,j0,kp)+vx(ip,jp,kp))
cc      fm = 0.25*(vx(i0,j0,k0)+vx(i0,jp,k0)
cc     .          +vx(i0,j0,kp)+vx(i0,jp,kp))
cc
cc      div = (fp-fm)/dhx
cc
cc      !Y component
cc      fp = 0.25*(vy(i0,jp,k0)+vy(ip,jp,k0)
cc     .          +vy(i0,jp,kp)+vy(ip,jp,kp))
cc      fm = 0.25*(vy(i0,j0,k0)+vy(ip,j0,k0)
cc     .          +vy(i0,j0,kp)+vy(ip,j0,kp))
cc
cc      div = div + (fp-fm)/dhy
cc
cc      !Z component
cc      fp = 0.25*(vz(ip,j0,kp)+vz(i0,j0,kp)
cc     .          +vz(ip,jp,kp)+vz(i0,jp,kp))
cc      fm = 0.25*(vz(ip,j0,k0)+vz(i0,j0,k0)
cc     .          +vz(ip,jp,k0)+vz(i0,jp,k0))
cc
cc      div = div + (fp-fm)/dhz
cc
cc      !Jacobian factor
cc      if (vrtx) then
cc        jac = 0.125*(gmetric%grid(igx)%jac(i  ,j  ,k  )
cc     .              +gmetric%grid(igx)%jac(i+1,j  ,k  )
cc     .              +gmetric%grid(igx)%jac(i  ,j+1,k  )
cc     .              +gmetric%grid(igx)%jac(i+1,j+1,k  ) 
cc     .              +gmetric%grid(igx)%jac(i  ,j  ,k+1)
cc     .              +gmetric%grid(igx)%jac(i+1,j  ,k+1)
cc     .              +gmetric%grid(igx)%jac(i  ,j+1,k+1)
cc     .              +gmetric%grid(igx)%jac(i+1,j+1,k+1))
cc      else
cc        jac = gmetric%grid(igx)%jac(i,j,k)
cc      endif
cc
cc      if (isSP(i+1,j,k,igx,igy,igz).and.vrtx) then
cc        div = 0d0
cc      else
cc        div = div/jac
cc      endif
cc
ccc     End program
cc
cc      end function div_vrtx
cc
ccc     grad_vrtx
ccc     ##################################################################
cc      function grad_vrtx(i,j,k,igx,igy,igz,phi,vrtx) result(grd)
cc
ccc     ------------------------------------------------------------------
ccc     Calculates gradient of phi from nearest neighbors.
ccc     If vrtx is true, it finds grad at (i+1/2,j+1/2,k+1/2) from (i,j,k)
ccc     If vrtx is false, it finds grad at (i,j,k) from (i+1/2,j+1/2,k+1/2)
ccc     ------------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: i,j,k,igx,igy,igz
cc      real(8) :: phi(2,2,2),grd(3)
cc      logical :: vrtx
cc
ccc     Local variables
cc
cc      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0,nx,ny,nz
cc      real(8) :: idhx,idhy,idhz,fp,fm
cc
ccc     Begin program
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      nx = grid_params%nxgl(igx)
cc      ny = grid_params%nygl(igy)
cc      nz = grid_params%nzgl(igz)
cc
cc      idhx = 0d0
cc      idhy = 0d0
cc      idhz = 0d0
cc
cc      if (vrtx) then
cc        if (nx > 1) idhx = 1./dx(ig)
cc        if (ny > 1) idhy = 1./dy(jg)
cc        if (nz > 1) idhz = 1./dz(kg)
cc      else
cc        if (nx > 1) idhx = 1./dxh(ig)
cc        if (ny > 1) idhy = 1./dyh(jg)
cc        if (nz > 1) idhz = 1./dzh(kg)
cc      endif
cc
cccc      if (vrtx) then
cccc        idhx = 1./dx(ig)
cccc        idhy = 1./dy(jg)
cccc        idhz = 1./dz(kg)
cccc      else
cccc        idhx = 1./dxh(ig)
cccc        idhy = 1./dyh(jg)
cccc        idhz = 1./dzh(kg)
cccc      endif
cccc
cccc      if (vrtx) then
cccc        ip = 1 + min(grid_params%nxgl(igx)-(i+grid_params%ilo(igx)-1),1)
cccc        i0 = 2 - min(i+grid_params%ilo(igx)-1,1)
cccc        jp = 1 + min(grid_params%nygl(igy)-(j+grid_params%jlo(igy)-1),1)
cccc        j0 = 2 - min(j+grid_params%jlo(igy)-1,1)
cccc        kp = 1 + min(grid_params%nzgl(igz)-(k+grid_params%klo(igz)-1),1)
cccc        k0 = 2 - min(k+grid_params%klo(igz)-1,1)
cccc
cccc        write (*,*) i,j,k,"  ;  ",i0,ip,j0,jp,k0,kp
cccc
cccc      else
cc        ip = 2
cc        i0 = 1
cc        jp = 2
cc        j0 = 1
cc        kp = 2
cc        k0 = 1
cccc      endif
cc
cc      !X component
cc      fp = 0.25*(phi(ip,j0,k0)+phi(ip,jp,k0)
cc     .          +phi(ip,j0,kp)+phi(ip,jp,kp))
cc      fm = 0.25*(phi(i0,j0,k0)+phi(i0,jp,k0)
cc     .          +phi(i0,j0,kp)+phi(i0,jp,kp))
cc
cc      grd(1) = (fp-fm)*idhx
cc
cc      !Y component
cc      fp = 0.25*(phi(i0,jp,k0)+phi(ip,jp,k0)
cc     .          +phi(i0,jp,kp)+phi(ip,jp,kp))
cc      fm = 0.25*(phi(i0,j0,k0)+phi(ip,j0,k0)
cc     .          +phi(i0,j0,kp)+phi(ip,j0,kp))
cc
cc      grd(2) = (fp-fm)*idhy
cc
cc      !Z component
cc      fp = 0.25*(phi(ip,j0,kp)+phi(i0,j0,kp)
cc     .          +phi(ip,jp,kp)+phi(i0,jp,kp))
cc      fm = 0.25*(phi(ip,j0,k0)+phi(i0,j0,k0)
cc     .          +phi(ip,jp,k0)+phi(i0,jp,k0))
cc
cc      grd(3) = (fp-fm)*idhz
cc
ccc     End program
cc
cc      end function grad_vrtx
cc

ccc     grad_vrtx
ccc     ##################################################################
cc      function grad_vrtx(i,j,k,igx,igy,igz,phi) result(grd)
cc
ccc     ------------------------------------------------------------------
ccc     Calculates gradient at vertex (i+1/2,j+1/2,k+1/2) from scalar
ccc     defined at surrounding cell centers
ccc     ------------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: i,j,k,igx,igy,igz
cc      real(8) :: phi(2,2,2),grd(3)
cc
ccc     Local variables
cc
cc      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0,nx,ny,nz
cc      real(8) :: idhx,idhy,idhz,fp,fm
cc
ccc     Begin program
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      idhx = 0d0
cc      idhy = 0d0
cc      idhz = 0d0
cc
cc      nx = grid_params%nxgl(igx)
cc      ny = grid_params%nygl(igy)
cc      nz = grid_params%nzgl(igz)
cc
cc      if (nx > 1) idhx = 1d0/grid_params%dx(ig)
cc      if (ny > 1) idhy = 1d0/grid_params%dy(jg)
cc      if (nz > 1) idhz = 1d0/grid_params%dz(kg)
cc
cc      ip = 2
cc      i0 = 1
cc      jp = 2
cc      j0 = 1
cc      kp = 2
cc      k0 = 1
cc
cc      !X component
cc      fp = 0.25*(phi(ip,j0,k0)+phi(ip,jp,k0)
cc     .          +phi(ip,j0,kp)+phi(ip,jp,kp))
cc      fm = 0.25*(phi(i0,j0,k0)+phi(i0,jp,k0)
cc     .          +phi(i0,j0,kp)+phi(i0,jp,kp))
cc
cc      grd(1) = (fp-fm)*idhx
cc
cc      !Y component
cc      fp = 0.25*(phi(i0,jp,k0)+phi(ip,jp,k0)
cc     .          +phi(i0,jp,kp)+phi(ip,jp,kp))
cc      fm = 0.25*(phi(i0,j0,k0)+phi(ip,j0,k0)
cc     .          +phi(i0,j0,kp)+phi(ip,j0,kp))
cc
cc      grd(2) = (fp-fm)*idhy
cc
cc      !Z component
cc      fp = 0.25*(phi(ip,j0,kp)+phi(i0,j0,kp)
cc     .          +phi(ip,jp,kp)+phi(i0,jp,kp))
cc      fm = 0.25*(phi(ip,j0,k0)+phi(i0,j0,k0)
cc     .          +phi(ip,jp,k0)+phi(i0,jp,k0))
cc
cc      grd(3) = (fp-fm)*idhz
cc
ccc     End program
cc
cc      end function grad_vrtx
cc
ccc     curl_vrtx
ccc     ##################################################################
cc      function curl_vrtx(i,j,k,igx,igy,igz,vv,vrtx) result(cnv)
cc
ccc     ------------------------------------------------------------------
ccc     Calculates curl of vv (covariant) from nearest neighbors. Returns
ccc     contravariant components.
ccc
ccc     If vrtx is true, it finds curl at (i+1/2,j+1/2,k+1/2) from (i,j,k)
ccc     If vrtx is false, it finds curl at (i,j,k) from (i+1/2,j+1/2,k+1/2)
ccc     ------------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: i,j,k,igx,igy,igz
cc      real(8) :: vv(2,2,2,3),cnv(3)
cc      logical :: vrtx
cc
ccc     Local variables
cc
cc      real(8) :: gvx(3),gvy(3),gvz(3)
cc
ccc     Begin program
cc
cc      gvx = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,1),vrtx)
cc      gvy = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,2),vrtx)
cc      gvz = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,3),vrtx)
cc
cc      cnv(1) = gvz(2)-gvy(3)
cc      cnv(2) = gvx(3)-gvz(1)
cc      cnv(3) = gvy(1)-gvx(2)
cc
ccc     End program
cc
cc      end function curl_vrtx
cc
ccc     graddiv_vrtx
ccc     ###############################################################
cc      function graddiv_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
cc     .         result(cov)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates grad(div(A)) at cell centers in general non-orthogonal
ccc     coordinates, preserving the self-adjoint property. The vector A
ccc     is contravariant, and returns covariant components.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cov(3)
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: il,jl,kl,ii,jj,kk
cc
cc      real(8)    :: div(2,2,2)
cc
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
ccc     Divergence at vertices
cc
cc      do kl=1,2
cc        do jl=1,2
cc          do il=1,2
cc            ii = i+il-2
cc            jj = j+jl-2
cc            kk = k+kl-2
cc
cc            div(il,jl,kl) = div_vrtx(ii,jj,kk,igx,igy,igz
cc     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,1)
cc     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,2)
cc     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,3)
cc     .                         ,.true.)
cc          enddo
cc        enddo
cc      enddo
cc
ccc     Gradient at cell centers
cc
cc      cov = grad_vrtx(i,j,k,igx,igy,igz,div,.false.)
cc
cc      if (vol_wgt) cov=cov*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function graddiv_vrtx

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
c     Calculates curl(A) in general non-orthogonal coordinates. Is
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

cc      if (PRESENT(he)) then
cc        half_elem = he
cc
cc        select case(he)
cc        case(1)
cc          hx = 1
cc        case(2)
cc          hy = 1
cc        case(3)
cc          hz = 1
cc        end select
cc      else
cc        half_elem = 0

        if (PRESENT(hex)) hx = hex
        if (PRESENT(hey)) hy = hey
        if (PRESENT(hez)) hz = hez
cc      endif

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

cc      if (isSP2(i,igx).and.half_elem == 0.and.spoint) then
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

        INTERFACE laplacian
cc          module procedure lap_diff_old,lap_ndiff_old,lap_tnsr
          module procedure lap_diff,lap_ndiff,lap_tnsr
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

c     lap_diff_old
c     ###############################################################
      function lap_diff_old(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol)
     .         result (laplacian)
c     ---------------------------------------------------------------
c     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      real(8)    :: dff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1
      
      if (i == 0 .or. j == 0 .or. k == 0 ) then
        write (*,*) 'Error in laplace; i,j,k=0'
      elseif (i == nx+1 .or. j == ny+1 .or. k == nz+1) then
        write (*,*) 'Error in laplace; i,j,k=nmax+1'
      endif

      gsup => gmetric%grid(igrid)%gsup

      d_xx_ip = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
     .              +dff(ip,j,k)*gsup(ip,j,k,1,1))
      d_xx_im = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
     .              +dff(im,j,k)*gsup(im,j,k,1,1))
      d_yy_jp = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
     .              +dff(i,jp,k)*gsup(i,jp,k,2,2))
      d_yy_jm = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
     .              +dff(i,jm,k)*gsup(i,jm,k,2,2))
      d_zz_kp = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
     .              +dff(i,j,kp)*gsup(i,j,kp,3,3))
      d_zz_km = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
     .              +dff(i,j,km)*gsup(i,j,km,3,3))

      d_xy_ipjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                 +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                 +dff(ip,jp,k)*gsup(ip,jp,k,1,2))
      d_xy_ipjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                 +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
     .                 +dff(ip,jm,k)*gsup(ip,jm,k,1,2))
      d_xy_imjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(im,j ,k)*gsup(im,j ,k,1,2)
     .                 +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                 +dff(im,jp,k)*gsup(im,jp,k,1,2))
      d_xy_imjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(im,j ,k)*gsup(im,j ,k,1,2)
     .                 +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
     .                 +dff(im,jm,k)*gsup(im,jm,k,1,2))
                                                      
      d_xz_ipkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                 +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                 +dff(ip,j,kp)*gsup(ip,j,kp,1,3))
      d_xz_ipkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                 +dff(i ,j,km)*gsup(i ,j,km,1,3)
     .                 +dff(ip,j,km)*gsup(ip,j,km,1,3))
      d_xz_imkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(im,j,k )*gsup(im,j,k ,1,3)
     .                 +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                 +dff(im,j,kp)*gsup(im,j,kp,1,3))
      d_xz_imkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(im,j,k )*gsup(im,j,k ,1,3)
     .                 +dff(i ,j,km)*gsup(i ,j,km,1,3)
     .                 +dff(im,j,km)*gsup(im,j,km,1,3))
                                                      
      d_yz_jpkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                 +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                 +dff(i,jp,kp)*gsup(i,jp,kp,2,3))
      d_yz_jpkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                 +dff(i,j ,km)*gsup(i,j ,km,2,3)
     .                 +dff(i,jp,km)*gsup(i,jp,km,2,3))
      d_yz_jmkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jm,k )*gsup(i,jm,k ,2,3)
     .                 +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                 +dff(i,jm,kp)*gsup(i,jm,kp,2,3))
      d_yz_jmkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jm,k )*gsup(i,jm,k ,2,3)
     .                 +dff(i,j ,km)*gsup(i,j ,km,2,3)
     .                 +dff(i,jm,km)*gsup(i,jm,km,2,3))

      laplacian = 
     .     grid_params%dyh(jg)*grid_params%dzh(kg)
     .       *( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/grid_params%dx(ig)
     .         +d_xx_im*(arr(im,j,k)-arr(i,j,k))/grid_params%dx(ig-1) )
     .    +grid_params%dxh(ig)*grid_params%dzh(kg)
     .       *( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/grid_params%dy(jg)
     .         +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/grid_params%dy(jg-1) )
     .    +grid_params%dxh(ig)*grid_params%dyh(jg)*(
     .          d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/grid_params%dz(kg)
     .         +d_zz_km*(arr(i,j,km)-arr(i,j,k))/grid_params%dz(kg-1) )
     .    +0.5*( grid_params%dzh(kg)*(
     .                     d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
     .                    +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
     .                    -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
     .                    -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
     .          +grid_params%dyh(jg)*(
     .                     d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
     .                    +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
     .                    -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
     .                    -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
     .          +grid_params%dxh(ig)*(
     .                     d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
     .                    +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
     .                    -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
     .                    -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )
cc     .     ( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/dx(ig)
cc     .      +d_xx_im*(arr(im,j,k)-arr(i,j,k))/dx(ig-1) )/dxh(ig)
cc     .    +( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/dy(jg)
cc     .      +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/dy(jg-1) )/dyh(jg)
cc     .    +( d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/dz(kg)
cc     .      +d_zz_km*(arr(i,j,km)-arr(i,j,k))/dz(kg-1) )/dzh(kg)
cc     .    +0.5*(( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
cc     .           +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
cc     .           -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
cc     .           -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )/dxh(ig)/dyh(jg)
cc     .         +( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
cc     .           +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
cc     .           -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
cc     .           -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )/dxh(ig)/dzh(kg)
cc     .         +( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
cc     .           +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
cc     .           -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
cc     .           -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) )/dyh(jg)/dzh(kg))
cc     .    /gmetric%grid(igx)%jac(i,j,k)

      if (.not.vol_wgt) laplacian=laplacian
     .                           /gmetric%grid(igx)%dvol(i,j,k)
cc      if (vol_wgt) laplacian=laplacian
cc     .                      *gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_diff_old

c     lap_diff
c     ###############################################################
      function lap_diff(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol)
c     ---------------------------------------------------------------
c     Calculates div(dff grad(arr)) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),lap_diff

      real(8)    :: dff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
     .             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      !Average metrics (defined at vertices)
      gsup => gmetric%grid(igrid)%gsup_v

      d_xx_ip = 0.25*(gsup(i ,j ,k ,1,1)
     .               +gsup(i ,jm,k ,1,1)
     .               +gsup(i ,j ,km,1,1)
     .               +gsup(i ,jm,km,1,1))
      d_xx_im = 0.25*(gsup(im,j ,k ,1,1)
     .               +gsup(im,jm,k ,1,1)
     .               +gsup(im,j ,km,1,1)
     .               +gsup(im,jm,km,1,1))

      d_yy_jp = 0.25*(gsup(i ,j ,k ,2,2)
     .               +gsup(im,j ,k ,2,2)
     .               +gsup(i ,j ,km,2,2)
     .               +gsup(im,j ,km,2,2))
      d_yy_jm = 0.25*(gsup(i ,jm,k ,2,2)
     .               +gsup(im,jm,k ,2,2)
     .               +gsup(i ,jm,km,2,2)
     .               +gsup(im,jm,km,2,2))

      d_zz_kp = 0.25*(gsup(i ,j ,k ,3,3)
     .               +gsup(im,j ,k ,3,3)
     .               +gsup(i ,jm,k ,3,3)
     .               +gsup(im,jm,k ,3,3))
      d_zz_km = 0.25*(gsup(i ,j ,km,3,3)
     .               +gsup(im,j ,km,3,3)
     .               +gsup(i ,jm,km,3,3)
     .               +gsup(im,jm,km,3,3))

      d_xy_ipjp = 0.5*(gsup(i ,j ,k ,1,2)
     .                +gsup(i ,j ,km,1,2))
      d_xy_ipjm = 0.5*(gsup(i ,jm,k ,1,2)
     .                +gsup(i ,jm,km,1,2))
      d_xy_imjp = 0.5*(gsup(im,j ,k ,1,2)
     .                +gsup(im,j ,km,1,2))
      d_xy_imjm = 0.5*(gsup(im,jm,k ,1,2)
     .                +gsup(im,jm,km,1,2))
                       
      d_xz_ipkp = 0.5*(gsup(i ,j ,k ,1,3)
     .                +gsup(i ,jm,k ,1,3))
      d_xz_ipkm = 0.5*(gsup(i ,j ,km,1,3)
     .                +gsup(i ,jm,km,1,3))
      d_xz_imkp = 0.5*(gsup(im,j ,k ,1,3)
     .                +gsup(im,jm,k ,1,3))
      d_xz_imkm = 0.5*(gsup(im,j ,km,1,3)
     .                +gsup(im,jm,km,1,3))

      d_yz_jpkp = 0.5*(gsup(i ,j ,k ,2,3)
     .                +gsup(im,j ,k ,2,3))
      d_yz_jpkm = 0.5*(gsup(i ,j ,km,2,3)
     .                +gsup(im,j ,km,2,3))
      d_yz_jmkp = 0.5*(gsup(i ,jm,k ,2,3)
     .                +gsup(im,jm,k ,2,3))
      d_yz_jmkm = 0.5*(gsup(i ,jm,km,2,3)
     .                +gsup(im,jm,km,2,3))

cc      if (isSP(i,j,k,igx,igy,igz).and.(.not.numerical_grid)) then !To fix singularity with gsup(2,2)
cc        gsup => gmetric%grid(igx)%gsup
cc
cc        d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
cc     .                +gsup(i,jp,k,2,2))
cc        d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
cc     .                +gsup(i,jm,k,2,2))
cc      endif

      !Add diffusion coefficient (defined at cell centers)
      d_xx_ip = 0.5*(dff(i ,j,k)+dff(ip,j,k))*d_xx_ip
      d_xx_im = 0.5*(dff(i ,j,k)+dff(im,j,k))*d_xx_im
      d_yy_jp = 0.5*(dff(i,j ,k)+dff(i,jp,k))*d_yy_jp
      d_yy_jm = 0.5*(dff(i,j ,k)+dff(i,jm,k))*d_yy_jm
      d_zz_kp = 0.5*(dff(i,j,k )+dff(i,j,kp))*d_zz_kp
      d_zz_km = 0.5*(dff(i,j,k )+dff(i,j,km))*d_zz_km

      d_xy_ipjp = 0.25*(dff(i ,j ,k)
     .                 +dff(ip,j ,k)
     .                 +dff(i ,jp,k)
     .                 +dff(ip,jp,k))*d_xy_ipjp
      d_xy_ipjm = 0.25*(dff(i ,j ,k)          
     .                 +dff(ip,j ,k)          
     .                 +dff(i ,jm,k)          
     .                 +dff(ip,jm,k))*d_xy_ipjm
      d_xy_imjp = 0.25*(dff(i ,j ,k)          
     .                 +dff(im,j ,k)          
     .                 +dff(i ,jp,k)          
     .                 +dff(im,jp,k))*d_xy_imjp
      d_xy_imjm = 0.25*(dff(i ,j ,k)          
     .                 +dff(im,j ,k)          
     .                 +dff(i ,jm,k)          
     .                 +dff(im,jm,k))*d_xy_imjm
                                              
      d_xz_ipkp = 0.25*(dff(i ,j,k )          
     .                 +dff(ip,j,k )          
     .                 +dff(i ,j,kp)          
     .                 +dff(ip,j,kp))*d_xz_ipkp
      d_xz_ipkm = 0.25*(dff(i ,j,k )           
     .                 +dff(ip,j,k )           
     .                 +dff(i ,j,km)           
     .                 +dff(ip,j,km))*d_xz_ipkm
      d_xz_imkp = 0.25*(dff(i ,j,k )           
     .                 +dff(im,j,k )           
     .                 +dff(i ,j,kp)           
     .                 +dff(im,j,kp))*d_xz_imkp
      d_xz_imkm = 0.25*(dff(i ,j,k )           
     .                 +dff(im,j,k )           
     .                 +dff(i ,j,km)           
     .                 +dff(im,j,km))*d_xz_imkm
                                              
      d_yz_jpkp = 0.25*(dff(i,j ,k )          
     .                 +dff(i,jp,k )          
     .                 +dff(i,j ,kp)          
     .                 +dff(i,jp,kp))*d_yz_jpkp
      d_yz_jpkm = 0.25*(dff(i,j ,k )           
     .                 +dff(i,jp,k )           
     .                 +dff(i,j ,km)           
     .                 +dff(i,jp,km))*d_yz_jpkm
      d_yz_jmkp = 0.25*(dff(i,j ,k )           
     .                 +dff(i,jm,k )           
     .                 +dff(i,j ,kp)           
     .                 +dff(i,jm,kp))*d_yz_jmkp
      d_yz_jmkm = 0.25*(dff(i,j ,k )           
     .                 +dff(i,jm,k )           
     .                 +dff(i,j ,km)           
     .                 +dff(i,jm,km))*d_yz_jmkm

      flxip = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                     
      flxim =-d_xx_im  *(arr(im,j ,k)-arr(i,j ,k))/grid_params%dx(ig-1)
     .  +.25*(d_xy_imjp*(arr(im,jp,k)+arr(i,jp,k)               
     .                  -arr(im,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_xy_imjm*(arr(im,j ,k)+arr(i,j ,k)               
     .                  -arr(im,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
     .  +     d_xz_imkp*(arr(im,j,kp)+arr(i,j,kp)               
     .                  -arr(im,j,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_xz_imkm*(arr(im,j,k )+arr(i,j,k )               
     .                  -arr(im,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                     
      flxjp = d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
     .                 -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                     
      flxjm =-d_yy_jm  *(arr(i ,jm,k)-arr(i ,j,k))/grid_params%dy(jg-1)
     .  +.25*(d_xy_ipjm*(arr(ip,jm,k)+arr(ip,j,k)               
     .                  -arr(i ,jm,k)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xy_imjm*(arr(i ,jm,k)+arr(i ,j,k)               
     .                  -arr(im,jm,k)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jmkp*(arr(i,jm,kp)+arr(i,j,kp)               
     .                  -arr(i,jm,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_yz_jmkm*(arr(i,jm,k )+arr(i,j,k )               
     .                  -arr(i,jm,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                     
      flxkp = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))
                                                                     
      flxkm =-d_zz_km  *(arr(i ,j,km)-arr(i ,j,k))/grid_params%dz(kg)
     .  +.25*(d_xz_ipkm*(arr(ip,j,km)+arr(ip,j,k)               
     .                 -arr(i ,j,km)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xz_imkm*(arr(i ,j,km)+arr(i ,j,k)               
     .                  -arr(im,j,km)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkm*(arr(i,jp,km)+arr(i,jp,k)               
     .                  -arr(i,j ,km)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_yz_jmkm*(arr(i,j ,km)+arr(i,j ,k)               
     .                  -arr(i,jm,km)-arr(i,jm,k))/grid_params%dy(jg-1))

      lap_diff= ((flxip -flxim)/grid_params%dxh(ig)
     .          +(flxjp -flxjm)/grid_params%dyh(jg)
     .          +(flxkp -flxkm)/grid_params%dzh(kg))
     .          /gmetric%grid(igx)%jac(i,j,k)

      if (vol_wgt) lap_diff=lap_diff*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_diff

c     lap_ndiff_old
c     ###############################################################
      function lap_ndiff_old(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol)
     .         result (laplacian)
c     ---------------------------------------------------------------
c     Calculates lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      logical    :: vol_wgt

c     Begin program

      if (PRESENT(vol)) then
        vol_wgt = vol
      else
        vol_wgt = .false.
      endif

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1
      
      if (i == 0 .or. j == 0 .or. k == 0 ) then
        write (*,*) 'Error in laplace; i,j,k=0'
      elseif (i == nx+1 .or. j == ny+1 .or. k == nz+1) then
        write (*,*) 'Error in laplace; i,j,k=nmax+1'
      endif

      if (coords == 'car') then
        laplacian = 
     .     grid_params%dyh(jg)*grid_params%dzh(kg)
     .       *( (arr(ip,j,k)-arr(i,j,k))/grid_params%dx(ig)
     .         +(arr(im,j,k)-arr(i,j,k))/grid_params%dx(ig-1) )
     .    +grid_params%dxh(ig)*grid_params%dzh(kg)
     .       *( (arr(i,jp,k)-arr(i,j,k))/grid_params%dy(jg)
     .         +(arr(i,jm,k)-arr(i,j,k))/grid_params%dy(jg-1) )
     .    +grid_params%dxh(ig)*grid_params%dyh(jg)
     .       *( (arr(i,j,kp)-arr(i,j,k))/grid_params%dz(kg)
     .         +(arr(i,j,km)-arr(i,j,k))/grid_params%dz(kg-1) )
      else

        gsup => gmetric%grid(igrid)%gsup

        d_xx_ip = 0.5*(gsup(i ,j,k,1,1)
     .                +gsup(ip,j,k,1,1))
        d_xx_im = 0.5*(gsup(i ,j,k,1,1)
     .                +gsup(im,j,k,1,1))

        d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
     .                +gsup(i,jp,k,2,2))
        d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
     .                +gsup(i,jm,k,2,2))

        d_zz_kp = 0.5*(gsup(i,j,k ,3,3)
     .                +gsup(i,j,kp,3,3))
        d_zz_km = 0.5*(gsup(i,j,k ,3,3)
     .                +gsup(i,j,km,3,3))

        d_xy_ipjp = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(ip,j ,k,1,2)
     .                   +gsup(i ,jp,k,1,2)
     .                   +gsup(ip,jp,k,1,2))
        d_xy_ipjm = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(ip,j ,k,1,2)
     .                   +gsup(i ,jm,k,1,2)
     .                   +gsup(ip,jm,k,1,2))
        d_xy_imjp = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(im,j ,k,1,2)
     .                   +gsup(i ,jp,k,1,2)
     .                   +gsup(im,jp,k,1,2))
        d_xy_imjm = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(im,j ,k,1,2)
     .                   +gsup(i ,jm,k,1,2)
     .                   +gsup(im,jm,k,1,2))

        d_xz_ipkp = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(ip,j,k ,1,3)
     .                   +gsup(i ,j,kp,1,3)
     .                   +gsup(ip,j,kp,1,3))
        d_xz_ipkm = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(ip,j,k ,1,3)
     .                   +gsup(i ,j,km,1,3)
     .                   +gsup(ip,j,km,1,3))
        d_xz_imkp = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(im,j,k ,1,3)
     .                   +gsup(i ,j,kp,1,3)
     .                   +gsup(im,j,kp,1,3))
        d_xz_imkm = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(im,j,k ,1,3)
     .                   +gsup(i ,j,km,1,3)
     .                   +gsup(im,j,km,1,3))

        d_yz_jpkp = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jp,k ,2,3)
     .                   +gsup(i,j ,kp,2,3)
     .                   +gsup(i,jp,kp,2,3))
        d_yz_jpkm = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jp,k ,2,3)
     .                   +gsup(i,j ,km,2,3)
     .                   +gsup(i,jp,km,2,3))
        d_yz_jmkp = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jm,k ,2,3)
     .                   +gsup(i,j ,kp,2,3)
     .                   +gsup(i,jm,kp,2,3))
        d_yz_jmkm = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jm,k ,2,3)
     .                   +gsup(i,j ,km,2,3)
     .                   +gsup(i,jm,km,2,3))

cc        if (isSP(i,j,k,igx,igy,igz)) then
cc          d_xx_im   = 0d0
cc          d_xy_imjp = 0d0
cc          d_xy_imjm = 0d0
cc          d_xz_imkp = 0d0
cc          d_xz_imkm = 0d0
cc        endif

        laplacian = 
     .     grid_params%dyh(jg)*grid_params%dzh(kg)
     .       *( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/grid_params%dx(ig)
     .         +d_xx_im*(arr(im,j,k)-arr(i,j,k))/grid_params%dx(ig-1) )
     .    +grid_params%dxh(ig)*grid_params%dzh(kg)
     .       *( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/grid_params%dy(jg)
     .         +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/grid_params%dy(jg-1) )
     .    +grid_params%dxh(ig)*grid_params%dyh(jg)
     .       *( d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/grid_params%dz(kg)
     .         +d_zz_km*(arr(i,j,km)-arr(i,j,k))/grid_params%dz(kg-1) )
     .    +0.5*( grid_params%dzh(kg)
     .                  *( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
     .                    +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
     .                    -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
     .                    -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
     .          +grid_params%dyh(jg)
     .                  *( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
     .                    +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
     .                    -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
     .                    -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
     .          +grid_params%dxh(ig)
     .                  *( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
     .                    +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
     .                    -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
     .                    -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )

      endif

      if (.not.vol_wgt) laplacian=laplacian
     .                           /gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_ndiff_old

c     lap_ndiff
c     ###############################################################
      function lap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol)
c     ---------------------------------------------------------------
c     Calculates lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr (0:nx+1,0:ny+1,0:nz+1),lap_ndiff

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
     .             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      logical    :: vol_wgt

c     Begin program

      if (PRESENT(vol)) then
        vol_wgt = vol
      else
        vol_wgt = .false.
      endif

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (coords == 'car') then
        lap_ndiff = 
     .     grid_params%dyh(jg)*grid_params%dzh(kg)
     .       *( (arr(ip,j,k)-arr(i,j,k))/grid_params%dx(ig)
     .         +(arr(im,j,k)-arr(i,j,k))/grid_params%dx(ig-1) )
     .    +grid_params%dxh(ig)*grid_params%dzh(kg)
     .       *( (arr(i,jp,k)-arr(i,j,k))/grid_params%dy(jg)
     .         +(arr(i,jm,k)-arr(i,j,k))/grid_params%dy(jg-1) )
     .    +grid_params%dxh(ig)*grid_params%dyh(jg)
     .       *( (arr(i,j,kp)-arr(i,j,k))/grid_params%dz(kg)
     .         +(arr(i,j,km)-arr(i,j,k))/grid_params%dz(kg-1) )

        if (.not.vol_wgt) lap_ndiff=lap_ndiff
     .                           /gmetric%grid(igx)%dvol(i,j,k)
      else

        gsup => gmetric%grid(igx)%gsup_v

        d_xx_ip = 0.25*(gsup(i ,j ,k ,1,1)
     .                 +gsup(i ,jm,k ,1,1)
     .                 +gsup(i ,j ,km,1,1)
     .                 +gsup(i ,jm,km,1,1))
        d_xx_im = 0.25*(gsup(im,j ,k ,1,1)
     .                 +gsup(im,jm,k ,1,1)
     .                 +gsup(im,j ,km,1,1)
     .                 +gsup(im,jm,km,1,1))

        d_yy_jp = 0.25*(gsup(i ,j ,k ,2,2)
     .                 +gsup(im,j ,k ,2,2)
     .                 +gsup(i ,j ,km,2,2)
     .                 +gsup(im,j ,km,2,2))
        d_yy_jm = 0.25*(gsup(i ,jm,k ,2,2)
     .                 +gsup(im,jm,k ,2,2)
     .                 +gsup(i ,jm,km,2,2)
     .                 +gsup(im,jm,km,2,2))

        d_zz_kp = 0.25*(gsup(i ,j ,k ,3,3)
     .                 +gsup(im,j ,k ,3,3)
     .                 +gsup(i ,jm,k ,3,3)
     .                 +gsup(im,jm,k ,3,3))
        d_zz_km = 0.25*(gsup(i ,j ,km,3,3)
     .                 +gsup(im,j ,km,3,3)
     .                 +gsup(i ,jm,km,3,3)
     .                 +gsup(im,jm,km,3,3))

        d_xy_ipjp = 0.5*(gsup(i ,j ,k ,1,2)
     .                  +gsup(i ,j ,km,1,2))
        d_xy_ipjm = 0.5*(gsup(i ,jm,k ,1,2)
     .                  +gsup(i ,jm,km,1,2))
        d_xy_imjp = 0.5*(gsup(im,j ,k ,1,2)
     .                  +gsup(im,j ,km,1,2))
        d_xy_imjm = 0.5*(gsup(im,jm,k ,1,2)
     .                  +gsup(im,jm,km,1,2))

        d_xz_ipkp = 0.5*(gsup(i ,j ,k ,1,3)
     .                  +gsup(i ,jm,k ,1,3))
        d_xz_ipkm = 0.5*(gsup(i ,j ,km,1,3)
     .                  +gsup(i ,jm,km,1,3))
        d_xz_imkp = 0.5*(gsup(im,j ,k ,1,3)
     .                  +gsup(im,jm,k ,1,3))
        d_xz_imkm = 0.5*(gsup(im,j ,km,1,3)
     .                  +gsup(im,jm,km,1,3))

        d_yz_jpkp = 0.5*(gsup(i ,j ,k ,2,3)
     .                  +gsup(im,j ,k ,2,3))
        d_yz_jpkm = 0.5*(gsup(i ,j ,km,2,3)
     .                  +gsup(im,j ,km,2,3))
        d_yz_jmkp = 0.5*(gsup(i ,jm,k ,2,3)
     .                  +gsup(im,jm,k ,2,3))
        d_yz_jmkm = 0.5*(gsup(i ,jm,km,2,3)
     .                  +gsup(im,jm,km,2,3))

cc        if (isSP(i,j,k,igx,igy,igz).and.(.not.numerical_grid)) then  !To fix singularity with gsup(2,2)
cc          gsup => gmetric%grid(igx)%gsup
cc
cc          d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
cc     .                  +gsup(i,jp,k,2,2))
cc          d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
cc     .                  +gsup(i,jm,k,2,2))
cc        endif

        flxip = d_xx_ip*(arr(ip,j ,k)-arr(i,j ,k))/grid_params%dx(ig)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)               
     .                  -arr(ip,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)               
     .                  -arr(ip,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
     .  +     d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)               
     .                  -arr(ip,j,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )               
     .                  -arr(ip,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                       
        flxim =-d_xx_im*(arr(im,j ,k)-arr(i,j ,k))/grid_params%dx(ig-1)
     .  +.25*(d_xy_imjp*(arr(im,jp,k)+arr(i,jp,k)               
     .                  -arr(im,j ,k)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_xy_imjm*(arr(im,j ,k)+arr(i,j ,k)               
     .                  -arr(im,jm,k)-arr(i,jm,k))/grid_params%dy(jg-1)
     .  +     d_xz_imkp*(arr(im,j,kp)+arr(i,j,kp)               
     .                  -arr(im,j,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_xz_imkm*(arr(im,j,k )+arr(i,j,k )               
     .                  -arr(im,j,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                       
        flxjp = d_yy_jp*(arr(i ,jp,k)-arr(i ,j,k))/grid_params%dy(jg)
     .  +.25*(d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)               
     .                 -arr(i ,jp,k)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)               
     .                  -arr(im,jp,k)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)               
     .                  -arr(i,jp,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )               
     .                  -arr(i,jp,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                       
        flxjm =-d_yy_jm*(arr(i ,jm,k)-arr(i ,j,k))/grid_params%dy(jg-1)
     .  +.25*(d_xy_ipjm*(arr(ip,jm,k)+arr(ip,j,k)               
     .                  -arr(i ,jm,k)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xy_imjm*(arr(i ,jm,k)+arr(i ,j,k)               
     .                  -arr(im,jm,k)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jmkp*(arr(i,jm,kp)+arr(i,j,kp)               
     .                  -arr(i,jm,k )-arr(i,j,k ))/grid_params%dz(kg)
     .  +     d_yz_jmkm*(arr(i,jm,k )+arr(i,j,k )               
     .                  -arr(i,jm,km)-arr(i,j,km))/grid_params%dz(kg-1))
                                                                       
        flxkp = d_zz_kp*(arr(i ,j,kp)-arr(i ,j,k))/grid_params%dz(kg)
     .  +.25*(d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)               
     .                  -arr(i ,j,kp)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)               
     .                  -arr(im,j,kp)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)               
     .                  -arr(i,j ,kp)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)               
     .                  -arr(i,jm,kp)-arr(i,jm,k))/grid_params%dy(jg-1))
                                                                       
        flxkm =-d_zz_km*(arr(i ,j,km)-arr(i ,j,k))/grid_params%dz(kg)
     .  +.25*(d_xz_ipkm*(arr(ip,j,km)+arr(ip,j,k)               
     .                  -arr(i ,j,km)-arr(i ,j,k))/grid_params%dx(ig)
     .  +     d_xz_imkm*(arr(i ,j,km)+arr(i ,j,k)               
     .                  -arr(im,j,km)-arr(im,j,k))/grid_params%dx(ig-1)
     .  +     d_yz_jpkm*(arr(i,jp,km)+arr(i,jp,k)               
     .                  -arr(i,j ,km)-arr(i,j ,k))/grid_params%dy(jg)
     .  +     d_yz_jmkm*(arr(i,j ,km)+arr(i,j ,k)               
     .                  -arr(i,jm,km)-arr(i,jm,k))/grid_params%dy(jg-1))

        lap_ndiff= ((flxip -flxim)/grid_params%dxh(ig)
     .             +(flxjp -flxjm)/grid_params%dyh(jg)
     .             +(flxkp -flxkm)/grid_params%dzh(kg))
     .             /gmetric%grid(igx)%jac(i,j,k)

        if (vol_wgt) lap_ndiff=lap_ndiff*gmetric%grid(igx)%dvol(i,j,k)

      endif

c     End program

      end function lap_ndiff

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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)        :: arr (0:nx+1,0:ny+1,0:nz+1),laplacian
      real(8),target :: tnsr(0:nx+1,0:ny+1,0:nz+1,6)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      logical    :: vol_wgt

c     Begin program

      if (PRESENT(vol)) then
        vol_wgt = vol
      else
        vol_wgt = .false.
      endif

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (i == 0 .or. j == 0 .or. k == 0 ) then
        write (*,*) 'Error in laplace; i,j,k=0'
      elseif (i == nx+1 .or. j == ny+1 .or. k == nz+1) then
        write (*,*) 'Error in laplace; i,j,k=nmax+1'
      endif

      d_xx_ip = 0.5*(tnsr(i ,j,k,1)
     .              +tnsr(ip,j,k,1))
      d_xx_im = 0.5*(tnsr(i ,j,k,1)
     .              +tnsr(im,j,k,1))

      d_yy_jp = 0.5*(tnsr(i,j ,k,4)
     .              +tnsr(i,jp,k,4))
      d_yy_jm = 0.5*(tnsr(i,j ,k,4)
     .              +tnsr(i,jm,k,4))

      d_zz_kp = 0.5*(tnsr(i,j,k ,6)
     .              +tnsr(i,j,kp,6))
      d_zz_km = 0.5*(tnsr(i,j,k ,6)
     .              +tnsr(i,j,km,6))

      d_xy_ipjp = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(ip,j ,k,2)
     .                 +tnsr(i ,jp,k,2)
     .                 +tnsr(ip,jp,k,2))
      d_xy_ipjm = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(ip,j ,k,2)
     .                 +tnsr(i ,jm,k,2)
     .                 +tnsr(ip,jm,k,2))
      d_xy_imjp = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(im,j ,k,2)
     .                 +tnsr(i ,jp,k,2)
     .                 +tnsr(im,jp,k,2))
      d_xy_imjm = 0.25*(tnsr(i ,j ,k,2)
     .                 +tnsr(im,j ,k,2)
     .                 +tnsr(i ,jm,k,2)
     .                 +tnsr(im,jm,k,2))

      d_xz_ipkp = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(ip,j,k ,3)
     .                 +tnsr(i ,j,kp,3)
     .                 +tnsr(ip,j,kp,3))
      d_xz_ipkm = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(ip,j,k ,3)
     .                 +tnsr(i ,j,km,3)
     .                 +tnsr(ip,j,km,3))
      d_xz_imkp = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(im,j,k ,3)
     .                 +tnsr(i ,j,kp,3)
     .                 +tnsr(im,j,kp,3))
      d_xz_imkm = 0.25*(tnsr(i ,j,k ,3)
     .                 +tnsr(im,j,k ,3)
     .                 +tnsr(i ,j,km,3)
     .                 +tnsr(im,j,km,3))

      d_yz_jpkp = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jp,k ,5)
     .                 +tnsr(i,j ,kp,5)
     .                 +tnsr(i,jp,kp,5))
      d_yz_jpkm = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jp,k ,5)
     .                 +tnsr(i,j ,km,5)
     .                 +tnsr(i,jp,km,5))
      d_yz_jmkp = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jm,k ,5)
     .                 +tnsr(i,j ,kp,5)
     .                 +tnsr(i,jm,kp,5))
      d_yz_jmkm = 0.25*(tnsr(i,j ,k ,5)
     .                 +tnsr(i,jm,k ,5)
     .                 +tnsr(i,j ,km,5)
     .                 +tnsr(i,jm,km,5))

      laplacian = 
     .     grid_params%dyh(jg)*grid_params%dzh(kg)
     .       *( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/grid_params%dx(ig)
     .         +d_xx_im*(arr(im,j,k)-arr(i,j,k))/grid_params%dx(ig-1) )
     .    +grid_params%dxh(ig)*grid_params%dzh(kg)
     .       *( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/grid_params%dy(jg)
     .         +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/grid_params%dy(jg-1) )
     .    +grid_params%dxh(ig)*grid_params%dyh(jg)
     .       *( d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/grid_params%dz(kg)
     .         +d_zz_km*(arr(i,j,km)-arr(i,j,k))/grid_params%dz(kg-1) )
     .    +0.5*( grid_params%dzh(kg)
     .                  *( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
     .                    +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
     .                    -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
     .                    -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
     .          +grid_params%dyh(jg)
     .                  *( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
     .                    +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
     .                    -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
     .                    -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
     .          +grid_params%dxh(ig)
     .                  *( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
     .                    +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
     .                    -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
     .                    -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )

      if (.not.vol_wgt) laplacian=laplacian
     .                           /gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_tnsr

      end module grid_nabla2_ops

c grid_tensor_ops
c ######################################################################
      module grid_tensor_ops

        use grid_create

        use grid_nabla2_ops

        INTERFACE veclaplacian
          module procedure veclap_diff,veclap_ndiff
        end INTERFACE

        INTERFACE veclap_cov
          module procedure veclap_cov_diff,veclap_cov_ndiff
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

ccc     fnabla_v
ccc     #############################################################
cc      function fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
cc     .                 ,half_elem,cnv)
cc     .         result(tensor)
ccc     -------------------------------------------------------------
ccc     Calculates the tensor components T_l^m of nabla(vec v) and
ccc     fills tensor(l,m) at the following positions:
ccc       + half_elem=0 --> i,j,k
ccc       + half_elem=1 --> i+1/2,j,k
ccc       + half_elem=2 --> i,j+1/2,k
ccc       + half_elem=3 --> i,j,k+1/2
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
cc        real(8)    :: tensor(3,3)
cc        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
cc     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
cc     .               ,az(0:nx+1,0:ny+1,0:nz+1)
cc
cc        logical,optional :: cnv
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
cc        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
cc        real(8)    :: vxx,vyy,vzz
cc     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc        real(8)    :: tsrc(3,3),jacip,jacjp,jacjm,jackp,jackm
cc     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac,jach
cc     .               ,gsuper(3,3)
cc
cc        real(8), pointer, dimension(:,:,:) :: jb
cc
cc        logical    :: cntv
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        !Defaults
cc
cc        if (PRESENT(cnv)) then
cc          cntv = cnv
cc        else
cc          cntv = .false.
cc        endif
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        dhx = 2.*dxh(ig)
cc        dhy = 2.*dyh(jg)
cc        dhz = 2.*dzh(kg)
cc
cc        ip = i+1
cc        im = i-1
cc        jp = j+1
cc        jm = j-1
cc        kp = k+1
cc        km = k-1
cc
cc        jb => gmetric%grid(igrid)%jac
cc
cc        jach = jb(i,j,k)
cc
cc        if(cntv) gsuper=gmetric%grid(igrid)%gsup(i,j,k ,:,:)
cc
cc        !Exceptions
cc
cc        select case(half_elem)
cc        case (1)
cc          dhx = dx(ig)
cc
cc          if (isSP(ip,j,k,igx,igy,igz)) then
cc            jach = 0d0
cc          else
cc            jach = 0.5*(jb(i,j,k)+jb(ip,j,k))
cc          endif
cc
cc          if (cntv) then
cc            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,:,:)
cc     .                   +gmetric%grid(igrid)%gsup(ip,j,k,:,:))
cc
cc            if (isSP(ip,j,k,igx,igy,igz)) then
cc              gsuper(1,1) = 0d0
cc              gsuper(2,2) = 0d0
cc              gsuper(3,3) = 0d0
cccc              write (*,*) 'DIAG -- fnabla_v',gsuper
cc            endif
cc          endif
cc  
cc          vxip = ax(ip,j,k)/jb(ip,j,k)
cc          vxim = ax(i ,j,k)/jb(i ,j,k)
cc          vyip = ay(ip,j,k)/jb(ip,j,k)
cc          vyim = ay(i ,j,k)/jb(i ,j,k)
cc          vzip = az(ip,j,k)/jb(ip,j,k)
cc          vzim = az(i ,j,k)/jb(i ,j,k)
cc
cc          vxjp = 0.5*(ax(i,jp,k)/jb(i,jp,k)+ax(ip,jp,k)/jb(ip,jp,k))
cc          vxjm = 0.5*(ax(i,jm,k)/jb(i,jm,k)+ax(ip,jm,k)/jb(ip,jm,k))
cc          vyjp = 0.5*(ay(i,jp,k)/jb(i,jp,k)+ay(ip,jp,k)/jb(ip,jp,k))
cc          vyjm = 0.5*(ay(i,jm,k)/jb(i,jm,k)+ay(ip,jm,k)/jb(ip,jm,k))
cc          vzjp = 0.5*(az(i,jp,k)/jb(i,jp,k)+az(ip,jp,k)/jb(ip,jp,k))
cc          vzjm = 0.5*(az(i,jm,k)/jb(i,jm,k)+az(ip,jm,k)/jb(ip,jm,k))
cc                                      
cc          vxkp = 0.5*(ax(i,j,kp)/jb(i,j,kp)+ax(ip,j,kp)/jb(ip,j,kp))
cc          vxkm = 0.5*(ax(i,j,km)/jb(i,j,km)+ax(ip,j,km)/jb(ip,j,km))
cc          vykp = 0.5*(ay(i,j,kp)/jb(i,j,kp)+ay(ip,j,kp)/jb(ip,j,kp))
cc          vykm = 0.5*(ay(i,j,km)/jb(i,j,km)+ay(ip,j,km)/jb(ip,j,km))
cc          vzkp = 0.5*(az(i,j,kp)/jb(i,j,kp)+az(ip,j,kp)/jb(ip,j,kp))
cc          vzkm = 0.5*(az(i,j,km)/jb(i,j,km)+az(ip,j,km)/jb(ip,j,km))
cc
cc          tsrc = 0.5*(nabla_v_src(ip,j,k)
cc     .               +nabla_v_src(i ,j,k))
cc
cc        case (2)
cc          dhy = dy(jg)
cc
cc          jach = 0.5*(jb(i,j,k)+jb(i,jp,k))
cc
cc          if (cntv) gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,:,:)
cc     .                           +gmetric%grid(igrid)%gsup(i,jp,k,:,:))
cc
cc          vxip = 0.5*(ax(ip,j,k)/jb(ip,j,k)+ax(ip,jp,k)/jb(ip,jp,k))
cc          vxim = 0.5*(ax(im,j,k)/jb(im,j,k)+ax(im,jp,k)/jb(im,jp,k))
cc          vyip = 0.5*(ay(ip,j,k)/jb(ip,j,k)+ay(ip,jp,k)/jb(ip,jp,k))
cc          vyim = 0.5*(ay(im,j,k)/jb(im,j,k)+ay(im,jp,k)/jb(im,jp,k))
cc          vzip = 0.5*(az(ip,j,k)/jb(ip,j,k)+az(ip,jp,k)/jb(ip,jp,k))
cc          vzim = 0.5*(az(im,j,k)/jb(im,j,k)+az(im,jp,k)/jb(im,jp,k))
cc
cc          vxjp = ax(i,jp,k)/jb(i,jp,k)
cc          vxjm = ax(i,j ,k)/jb(i,j ,k)
cc          vyjp = ay(i,jp,k)/jb(i,jp,k)
cc          vyjm = ay(i,j ,k)/jb(i,j ,k)
cc          vzjp = az(i,jp,k)/jb(i,jp,k)
cc          vzjm = az(i,j ,k)/jb(i,j ,k)
cc
cc          vxkp = 0.5*(ax(i,j,kp)/jb(i,j,kp)+ax(i,jp,kp)/jb(i,jp,kp))
cc          vxkm = 0.5*(ax(i,j,km)/jb(i,j,km)+ax(i,jp,km)/jb(i,jp,km))
cc          vykp = 0.5*(ay(i,j,kp)/jb(i,j,kp)+ay(i,jp,kp)/jb(i,jp,kp))
cc          vykm = 0.5*(ay(i,j,km)/jb(i,j,km)+ay(i,jp,km)/jb(i,jp,km))
cc          vzkp = 0.5*(az(i,j,kp)/jb(i,j,kp)+az(i,jp,kp)/jb(i,jp,kp))
cc          vzkm = 0.5*(az(i,j,km)/jb(i,j,km)+az(i,jp,km)/jb(i,jp,km))
cc
cc          tsrc = 0.5*(nabla_v_src(i,jp,k)
cc     .               +nabla_v_src(i,j ,k))
cc
cc        case (3)
cc          dhz = dz(kg)
cc
cc          jach = 0.5*(jb(i,j,k)+jb(i,j,kp))
cc
cc          if (cntv) gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,:,:)
cc     .                           +gmetric%grid(igrid)%gsup(i,j,kp,:,:))
cc
cc          vxip = 0.5*(ax(ip,j,k)/jb(ip,j,k)+ax(ip,j,kp)/jb(ip,j,kp))
cc          vxim = 0.5*(ax(im,j,k)/jb(im,j,k)+ax(im,j,kp)/jb(im,j,kp))
cc          vyip = 0.5*(ay(ip,j,k)/jb(ip,j,k)+ay(ip,j,kp)/jb(ip,j,kp))
cc          vyim = 0.5*(ay(im,j,k)/jb(im,j,k)+ay(im,j,kp)/jb(im,j,kp))
cc          vzip = 0.5*(az(ip,j,k)/jb(ip,j,k)+az(ip,j,kp)/jb(ip,j,kp))
cc          vzim = 0.5*(az(im,j,k)/jb(im,j,k)+az(im,j,kp)/jb(im,j,kp))
cc
cc          vxjp = 0.5*(ax(i,jp,k)/jb(i,jp,k)+ax(i,jp,kp)/jb(i,jp,kp))
cc          vxjm = 0.5*(ax(i,jm,k)/jb(i,jm,k)+ax(i,jm,kp)/jb(i,jm,kp))
cc          vyjp = 0.5*(ay(i,jp,k)/jb(i,jp,k)+ay(i,jp,kp)/jb(i,jp,kp))
cc          vyjm = 0.5*(ay(i,jm,k)/jb(i,jm,k)+ay(i,jm,kp)/jb(i,jm,kp))
cc          vzjp = 0.5*(az(i,jp,k)/jb(i,jp,k)+az(i,jp,kp)/jb(i,jp,kp))
cc          vzjm = 0.5*(az(i,jm,k)/jb(i,jm,k)+az(i,jm,kp)/jb(i,jm,kp))
cc
cc          vxkp = ax(i,j,kp)/jb(i,j,kp)
cc          vxkm = ax(i,j,k )/jb(i,j,k )
cc          vykp = ay(i,j,kp)/jb(i,j,kp)
cc          vykm = ay(i,j,k )/jb(i,j,k )
cc          vzkp = az(i,j,kp)/jb(i,j,kp)
cc          vzkm = az(i,j,k )/jb(i,j,k )
cc
cc          tsrc = 0.5*(nabla_v_src(i,j,kp)
cc     .               +nabla_v_src(i,j,k ))
cc
cc        case default
cc
cc          vxip = ax(ip,j,k)/jb(ip,j,k)
cc          vxim = ax(im,j,k)/jb(im,j,k)
cc          vyip = ay(ip,j,k)/jb(ip,j,k)
cc          vyim = ay(im,j,k)/jb(im,j,k)
cc          vzip = az(ip,j,k)/jb(ip,j,k)
cc          vzim = az(im,j,k)/jb(im,j,k)
cc                                      
cc          vxjp = ax(i,jp,k)/jb(i,jp,k)
cc          vxjm = ax(i,jm,k)/jb(i,jm,k)
cc          vyjp = ay(i,jp,k)/jb(i,jp,k)
cc          vyjm = ay(i,jm,k)/jb(i,jm,k)
cc          vzjp = az(i,jp,k)/jb(i,jp,k)
cc          vzjm = az(i,jm,k)/jb(i,jm,k)
cc                                      
cc          vxkp = ax(i,j,kp)/jb(i,j,kp)
cc          vxkm = ax(i,j,km)/jb(i,j,km)
cc          vykp = ay(i,j,kp)/jb(i,j,kp)
cc          vykm = ay(i,j,km)/jb(i,j,km)
cc          vzkp = az(i,j,kp)/jb(i,j,kp)
cc          vzkm = az(i,j,km)/jb(i,j,km)
cc
cc          tsrc = nabla_v_src(i,j,k)
cc
cc        end select
cc
cc        idhx = 1d0/dhx
cc        idhy = 1d0/dhy
cc        idhz = 1d0/dhz
cc
cc      ! l = 1, m = 1
cc        tensor(1,1) = (vxip-vxim)*idhx
cc
cc      ! l = 1, m = 2
cc        tensor(1,2) = (vyip-vyim)*idhx
cc
cc      ! l = 1, m = 3
cc        tensor(1,3) = (vzip-vzim)*idhx
cc
cc      ! l = 2, m = 1
cc        tensor(2,1) = (vxjp-vxjm)*idhy
cc
cc      ! l = 2, m = 2
cc        tensor(2,2) = (vyjp-vyjm)*idhy
cc
cc      ! l = 2, m = 3
cc        tensor(2,3) = (vzjp-vzjm)*idhy
cc
cc      ! l = 3, m = 1
cc        tensor(3,1) = (vxkp-vxkm)*idhz
cc
cc      ! l = 3, m = 2
cc        tensor(3,2) = (vykp-vykm)*idhz
cc
cc      ! l = 3, m = 3
cc        tensor(3,3) = (vzkp-vzkm)*idhz
cc
cc      ! Version I (requires uncommenting code in nabla_v_src)
cc
cccc        tensor = jach*tensor
cccc        if (cntv) tensor = matmul(gsuper,tensor)
cccc
cccc        tensor = tensor + tsrc
cc
cc      ! Version II
cc
cc        tensor = jach*tensor + tsrc
cc
cc        if (cntv) tensor = matmul(gsuper,tensor)
cc
ccc     End program
cc
cc      contains
cc
ccc     nabla_v_src
ccc     #############################################################
cc      function nabla_v_src(i,j,k) result(tensor)
cc
ccc     -------------------------------------------------------------
ccc     Finds geometric source of tensor nabla(v) at cell (i,j,k)
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k
cc        real(8)    :: tensor(3,3)
cc
ccc     Local variables
cc
cc        integer    :: l,m,n
cc        real(8)    :: hessian(3,3,3),gsuper(3,3)
cc
ccc     Begin program
cc
cc        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cc        !Code for Version I
cccc        gsuper  = gmetric%grid(igrid)%gsup (i,j,k,:,:)
cccc        if (cntv) then
cccc          do m=1,3
cccc            do n=1,3
cccc              do l=1,3
cccc                hessian(l,n,m)=dot_product(gsuper(n,:),hessian(l,:,m))
cccc              enddo
cccc            enddo
cccc          enddo
cccc        endif
cc
cc        do l=1,3
cc          do m=1,3
cc            tensor(l,m) = ax(i,j,k)*hessian(m,l,1)
cc     .                   +ay(i,j,k)*hessian(m,l,2)
cc     .                   +az(i,j,k)*hessian(m,l,3)
cc          enddo
cc        enddo
cc
ccc     End program
cc
cc      end function nabla_v_src
cc
ccccc     nabla_v_src
ccccc     #############################################################
cccc      function nabla_v_src(i,j,k) result(tensor)
cccc
ccccc     -------------------------------------------------------------
ccccc     Finds geometric source of tensor nabla(v) at cell (i,j,k)
ccccc     -------------------------------------------------------------
cccc
cccc        implicit none
cccc
ccccc     Call variables
cccc
cccc        integer    :: i,j,k
cccc        real(8)    :: tensor(3,3)
cccc
ccccc     Local variables
cccc
cccc        integer    :: l,m,vxh,vyh,vzh
cccc        real(8)    :: hessian(3,3,3)
cccc
ccccc     Begin program
cccc
cccc        select case(half_elem)
cccc        case (1)
cccc
cccc          hessian  = 0.5*(gmetric%grid(igrid)%Gamma(i ,j,k,:,:,:)
cccc     .                   +gmetric%grid(igrid)%Gamma(ip,j,k,:,:,:))
cccc
cccc          vxh = 0.5*(ax(ip,j,k)+ax(i,j,k))
cccc          vyh = 0.5*(ay(ip,j,k)+ay(i,j,k))
cccc          vzh = 0.5*(az(ip,j,k)+az(i,j,k))
cccc
cccc        case (2)
cccc
cccc          hessian  = 0.5*(gmetric%grid(igrid)%Gamma(i,j ,k,:,:,:)
cccc     .                   +gmetric%grid(igrid)%Gamma(i,jp,k,:,:,:))
cccc
cccc          vxh = 0.5*(ax(i,jp,k)+ax(i,j,k))
cccc          vyh = 0.5*(ay(i,jp,k)+ay(i,j,k))
cccc          vzh = 0.5*(az(i,jp,k)+az(i,j,k))
cccc
cccc        case (3)
cccc
cccc          hessian  = 0.5*(gmetric%grid(igrid)%Gamma(i,j,k ,:,:,:)
cccc     .                   +gmetric%grid(igrid)%Gamma(i,j,kp,:,:,:))
cccc
cccc          vxh = 0.5*(ax(i,j,kp)+ax(i,j,k))
cccc          vyh = 0.5*(ay(i,j,kp)+ay(i,j,k))
cccc          vzh = 0.5*(az(i,j,kp)+az(i,j,k))
cccc
cccc        case default
cccc
cccc          hessian  = gmetric%grid(igrid)%Gamma(i,j,k ,:,:,:)
cccc
cccc          vxh = ax(i,j,k)
cccc          vyh = ay(i,j,k)
cccc          vzh = az(i,j,k)
cccc
cccc        end select
cccc
cccc        do l=1,3
cccc          do m=1,3
cccc            tensor(l,m) = vxh*hessian(m,1,l)
cccc     .                   +vyh*hessian(m,2,l)
cccc     .                   +vzh*hessian(m,3,l)
cccc          enddo
cccc        enddo
cccc
ccccc     End program
cccc
cccc      end function nabla_v_src
cc
cc      end function fnabla_v

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

c     div_tensor
c     ###############################################################
      function div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,tsrx,tsry,tsrz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculatesthe contravariant components of div(tensor) at cell
c     centers in general non-orthogonal coordinates, given by:
c
c         div(T)^i = partial_j(T^ji/jac) + T^lm/jac Gamma(i)_lm)
c
c     The routines tsrx,tsry,tsrz must provide the contravariant
c     tersor components T^lm (if flag /= 0, they actually provide
c     T^lm/jac).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: divt(3)

      integer    :: i,j,k,igx,igy,igz,nx,ny,nz

      logical    :: alt_eom

      external   tsrx,tsry,tsrz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

      real(8)    :: jac,ijac,dvol,dS1,dS2,dS3,dxx,dyy,dzz

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: hess(3,3,3),msource,dx1,dx2,ll

      logical    :: vol_wgt

      real(8)    :: dum1,dum2,coeff

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

      jac = gmetric%grid(igrid)%jac(i,j,k)
      ijac = 1./jac

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

      if (coords /= 'car') then
        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
      endif

      call tsrx(i ,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11p,t12p,t13p, 1)
      call tsrx(im,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11m,t12m,t13m,-1)
      if (coords /= 'car')
     .    call tsrx(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11o,t12o,t13o,0)

      call tsry(i,j ,k,nx,ny,nz,igx,igy,igz,alt_eom,t21p,t22p,t23p, 1)
      call tsry(i,jm,k,nx,ny,nz,igx,igy,igz,alt_eom,t21m,t22m,t23m,-1)
      if (coords /= 'car')
     .    call tsry(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t21o,t22o,t23o,0)

      call tsrz(i,j,k ,nx,ny,nz,igx,igy,igz,alt_eom,t31p,t32p,t33p, 1)
      call tsrz(i,j,km,nx,ny,nz,igx,igy,igz,alt_eom,t31m,t32m,t33m,-1)
      if (coords /= 'car')
     .    call tsrz(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t31o,t32o,t33o,0)

      msource = 0d0

      !Component 1
      if (coords /= 'car') then
        msource =  (t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
     .             +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
     .             +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
     .             *ijac
      endif

      divt(1) =  (t11p - t11m)/dxx
     .          +(t21p - t21m)/dyy
     .          +(t31p - t31m)/dzz + msource

      !Component 2
      coeff = 1d0

      if (coords /= 'car') then
        if (alt_eom) then

          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
     .              -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
     .              -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
     .              -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))*ijac

          coeff = ijac

        else
          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
     .              *ijac
        endif

      endif

      divt(2) = ( (t12p - t12m)/dxx
     .           +(t22p - t22m)/dyy
     .           +(t32p - t32m)/dzz)*coeff  + msource

      !Component 3
      if (coords /= 'car') then
        msource = (t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
     .            +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
     .            +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
     .            *ijac
      endif

      divt(3) =  (t13p - t13m)/dxx
     .          +(t23p - t23m)/dyy
     .          +(t33p - t33m)/dzz + msource

      !Volume factor
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)

c     End program

      end function div_tensor

c     div_tensor_cov
c     ###############################################################
      function div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                       ,tsrx,tsry,tsrz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculates covariant components of div(tensor) at cell centers in
c     general non-orthogonal coordinates, which is given by:
c
c         div(T)_i = jac^(-1)*[ partial_j(T^j_i)-T^l_m Gamma(m)_il) ]
c
c     The routinex tsrx,tsry,tsrz must provide mixed representation 
c     tensor components T^l_m.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: divt(3)

      integer    :: i,j,k,igx,igy,igz,nx,ny,nz

      external   tsrx,tsry,tsrz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

      real(8)    :: jac,dvol,dS1,dS2,dS3,dxx,dyy,dzz

cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: hess(3,3,3),msource,dx1,dx2,ll

      logical    :: vol_wgt

      real(8)    :: dum1,dum2

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

      jac = gmetric%grid(igrid)%jac(i,j,k)

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

      dvol = dxx*dyy*dzz

      if (coords /= 'car') then
        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
      endif

      call tsrx(i ,j,k,nx,ny,nz,igx,igy,igz,t11p,t12p,t13p, 1)
      call tsrx(im,j,k,nx,ny,nz,igx,igy,igz,t11m,t12m,t13m,-1)
      if (coords /= 'car')
     .    call tsrx(i,j,k,nx,ny,nz,igx,igy,igz,t11o,t12o,t13o,0)

      call tsry(i,j ,k,nx,ny,nz,igx,igy,igz,t21p,t22p,t23p, 1)
      call tsry(i,jm,k,nx,ny,nz,igx,igy,igz,t21m,t22m,t23m,-1)
      if (coords /= 'car')
     .    call tsry(i,j,k,nx,ny,nz,igx,igy,igz,t21o,t22o,t23o,0)

      call tsrz(i,j,k ,nx,ny,nz,igx,igy,igz,t31p,t32p,t33p, 1)
      call tsrz(i,j,km,nx,ny,nz,igx,igy,igz,t31m,t32m,t33m,-1)
      if (coords /= 'car')
     .    call tsrz(i,j,k,nx,ny,nz,igx,igy,igz,t31o,t32o,t33o,0)

      msource = 0d0

      !Component 1
      if (coords /= 'car') then
        msource =  (t11o*hess(1,1,1)+t12o*hess(2,1,1)+t13o*hess(3,1,1)
     .             +t21o*hess(1,2,1)+t22o*hess(2,2,1)+t23o*hess(3,2,1)
     .             +t31o*hess(1,3,1)+t32o*hess(2,3,1)+t33o*hess(3,3,1))
      endif

      divt(1) =  (t11p - t11m)/dxx
     .          +(t21p - t21m)/dyy
     .          +(t31p - t31m)/dzz - msource

      !Component 2
      if (coords /= 'car') then
        msource=  (t11o*hess(1,1,2)+t12o*hess(2,1,2)+t13o*hess(3,1,2)
     .            +t21o*hess(1,2,2)+t22o*hess(2,2,2)+t23o*hess(3,2,2)
     .            +t31o*hess(1,3,2)+t32o*hess(2,3,2)+t33o*hess(3,3,2))
      endif

      divt(2) =  (t12p - t12m)/dxx
     .          +(t22p - t22m)/dyy
     .          +(t32p - t32m)/dzz  - msource

      !Component 3
      if (coords /= 'car') then
        msource = (t11o*hess(1,1,3)+t12o*hess(2,1,3)+t13o*hess(3,1,3)
     .            +t21o*hess(1,2,3)+t22o*hess(2,2,3)+t23o*hess(3,2,3)
     .            +t31o*hess(1,3,3)+t32o*hess(2,3,3)+t33o*hess(3,3,3))
      endif

      divt(3) =  (t13p - t13m)/dxx
     .          +(t23p - t23m)/dyy
     .          +(t33p - t33m)/dzz - msource

      divt = divt/jac

      !Volume factor
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)

c     End program

      end function div_tensor_cov

ccc     veclap2
ccc     ###############################################################
cc      function veclap2(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcov,vol)
cc     .         result (vlap)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates dvol*lap(vector) at cell centers in general non-orthog.
ccc     coordinates using
ccc                  lap(vector)=grad(div(v)-curl(curl(v)). 
ccc     Returns lap(vector) in contravariant components.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8),target :: vcnv(0:nx+1,0:ny+1,0:nz+1,3)
cc     .                 ,vcov(0:nx+1,0:ny+1,0:nz+1,3)
cc
cc      real(8)    :: vlap(3)
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: icomp
cc
cc      real(8)    :: cov(3)
cc
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
cc      if (coords == 'car') then
cc        do icomp=1,3
cc          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                           ,vcnv(:,:,:,icomp)
cc     .                           ,vol=vol_wgt)
cc        enddo
cc      else
cc        cov = graddiv(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vol=vol_wgt)
cc        vlap = matmul(gmetric%grid(igx)%gsup(i,j,k,:,:),cov) !Xforms to cnv
cc
cc        vlap = vlap - curlcurl(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                        ,vcov,vol=vol_wgt)
cc      endif
cc
ccc     End program
cc
cc      end function veclap2

c     veclap_diff
c     ###############################################################
      function veclap_diff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                     ,alteom,diff,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*div(diff grad(vector)) at cell centers in general
c     non-orthog.  coordinates. Vector is assumed in contravariant
c     representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)

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

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,vfield(:,:,:,icomp),diff
     .                           ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
        nullify(coef)
      endif

c     End program

      end function veclap_diff

c     veclap_ndiff
c     ###############################################################
      function veclap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                     ,alteom,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(vector) at cell centers in contravariant
c     representation. Vector vfield is contravariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)

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

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,vfield(:,:,:,icomp)
     .                           ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        nullify(coef)
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     .                   ,vol=vol_wgt)
      endif

c     End program

      end function veclap_ndiff

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
        real(8)    :: nabla_v_cnv(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

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

c     veclap_cov_ndiff
c     ###############################################################
      function veclap_cov_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                         ,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of dvol*lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,vfield(:,:,:,icomp),vol=vol_wgt)
        enddo
      else
        vec => vfield !Pointer passed to nabtensor routines
        nullify(coef)
        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
      endif

c     End program

      end function veclap_cov_ndiff

c     veclap_cov_diff
c     ###############################################################
      function veclap_cov_diff(i,j,k,nx,ny,nz,igx,igy,igz,vfield,diff
     .                         ,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of dvol*lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,vfield(:,:,:,icomp),diff
     .                         ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
        nullify(coef)
      endif

c     End program

      end function veclap_cov_diff

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

c     curlcurl
c     ###############################################################
      function curlcurl(i,j,k,nx,ny,nz,igx,igy,igz,a,vol) result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     contravariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: a(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)

      logical,optional,intent(IN) :: vol

c     Local variables

c     Begin program

      cnv=-veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,a,vol=vol)
     .    +gradDiv(i,j,k,nx,ny,nz,igx,igy,igz,a,vol=vol,cnv=.true.)

c     End program

      end function curlcurl

      end module grid_tensor_ops

c grid_advec_ops
c ######################################################################
      module grid_advec_ops

        use grid_create

        use grid_tensor_ops

        integer,private :: adv_meth

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
cc
c     c_advec
c     ###############################################################
      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,phi
     .                ,sp,upwind,vol)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.phi) at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),c_advec

      logical,optional :: sp,vol,upwind

c     Local variables

      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw,spoint,upw

c     Begin program

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .false.
      endif

      if (PRESENT(upwind)) then
        upw = upwind
      else
        upw = .false.
      endif

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = .false.
      endif

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

      jac   = gmetric%grid(igrid)%jac(i,j,k)

      if (spoint) then
        jacip = gmetric%grid(igrid)%jac(ip,j,k)
        jacim = gmetric%grid(igrid)%jac(im,j,k)
      else
        jacip = jac
        jacim = jac
      endif

      if (upw) then
        !X flux
        if (     i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .or. bcond(2) == PER) then
          jach = 0.5*(jac+jacip)
          flxip = 0.25*jach
     .      *( (    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)
     .          +abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(i ,j,k)
     .        +(    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)          
     .          -abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(ip,j,k) )
        else
          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
        endif

        if (     i+grid_params%ilo(igx)-1 > 1
     .      .or. spoint
     .      .or. bcond(1) == PER) then
          jach = 0.5*(jac+jacim)
          flxim = 0.25*jach
     .       *( (    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
     .           +abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(im,j,k)
     .         +(    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
     .           -abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(i ,j,k) )
        else
          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
        endif

        if (     j+grid_params%jlo(igy)-1 < grid_params%nygl(igy)
     .      .or. bcond(4) == PER) then
          flxjp = 0.25*(
     .         (    (v2(i,j,k)+v2(i,jp,k))
     .          +abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,j ,k)
     .        +(    (v2(i,j,k)+v2(i,jp,k))          
     .          -abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,jp,k))
        else
          flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
        endif

        if (     j+grid_params%jlo(igy)-1 > 1
     .      .or. bcond(3) == PER) then
          flxjm = 0.25*(
     .         (    (v2(i,j,k)+v2(i,jm,k))
     .          +abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,jm,k)
     .        +(    (v2(i,j,k)+v2(i,jm,k))          
     .          -abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,j ,k))
        else
          flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
        endif

        if (     k+grid_params%klo(igz)-1 < grid_params%nzgl(igz)
     .      .or. bcond(6) == PER) then
          flxkp = 0.25*(
     .         (    (v3(i,j,k)+v3(i,j,kp))
     .          +abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,k )
     .        +(    (v3(i,j,k)+v3(i,j,kp))             
     .          -abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,kp))
        else
          flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
        endif

        if (     k+grid_params%klo(igz)-1 > 1
     .      .or. bcond(5) == PER) then
          flxkm = 0.25*(
     .         (    (v3(i,j,k)+v3(i,j,km))
     .          +abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,km)
     .        +(    (v3(i,j,k)+v3(i,j,km))             
     .          -abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,k ))
        else
          flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
        endif

      else
        !X flux
        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
          jach = 0.5*(jac+jacip)
          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach

          jach = 0.5*(jac+jacim)
          if (isSP(i,j,k,igx,igy,igz)) jach = 0d0
          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))

          jach = 0.5*(jac+jacim)
          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
        endif

cc        if (i+grid_params%ilo(igx)-1 == 1) then
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
cc     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        elseif (i+grid_params%ilo(igx)-1 == 2) then
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        else
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cc        endif

        !Y flux
        flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
        flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))

        !Z flux
        flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
        flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))

      endif

      c_advec =( (flxip - flxim)/dxx
     .         + (flxjp - flxjm)/dyy
     .         + (flxkp - flxkm)/dzz )/jac
      
      if (vw) c_advec = c_advec*gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function c_advec

c     flx_advec
c     ###############################################################
      function flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                  ,method,sp,vol,zip_vel,conserv) result(advec)
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

      logical,optional :: sp,vol,zip_vel,conserv

c     Local variables

      integer :: igrid,half_elem,su

      integer :: im,ip,ipp,imm,igm,igp,igpp,igmm
     .          ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
     .          ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
     .          ,ig,jg,kg,igl,jgl,kgl

      real(8) :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
      real(8) :: vxp,vxm,vyp,vym,vzp,vzm
     .          ,qip,qim,qjp,qjm,qkp,qkm
      real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical :: vw,spoint,bc_safe,zv,cv

c     Begin program

c     Set defaults

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

      if (PRESENT(zip_vel)) then
        zv = zip_vel
      else
        zv = .true.
      endif

      if (PRESENT(conserv)) then
        cv = conserv
      else
        cv = .true.
      endif

      bc_safe = (method <=2)

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

      if (spoint) then
        jacip = gmetric%grid(igrid)%jac(ip,j,k)
        jacim = gmetric%grid(igrid)%jac(im,j,k)
      else
        jacip = jac
        jacim = jac
      endif

c     Find ZIP fluxes (used even if method /= ZIP)

      flxip = 0d0
      flxim = 0d0
      flxjp = 0d0
      flxjm = 0d0
      flxkp = 0d0
      flxkm = 0d0

      !X flux
      if (.not.isBdry(i,igx,2)) then  !Not at r=1 boundary
        jach  = 0.5*(jac+jacip)
        flxip = 0.5*(v(ip,j,k,1)*arr(i ,j,k)/jacip
     .             + v(i ,j,k,1)*arr(ip,j,k)/jac  )*jach

        if (isSP2(i,igx,ibc=1).and.spoint) then
          flxim = 0d0
        else
          jach  = 0.5*(jac+jacim)
          flxim = 0.5*(v(im,j,k,1)*arr(i ,j,k)/jacim
     .               + v(i ,j,k,1)*arr(im,j,k)/jac  )*jach
        endif
      else
        if (v(ip,j,k,1) == -v(i,j,k,1)) then
          flxip = 0d0    !Shut off flux for vn=0 BC
        else
          flxip = 0.5*(v(ip,j,k,1)*arr(i,j,k)+v(i,j,k,1)*arr(ip,j,k))
        endif

        jach  = 0.5*(jac+jacim)
        flxim = 0.5*(v(im,j,k,1)*arr(i ,j,k)/jacim
     .             + v(i ,j,k,1)*arr(im,j,k)/jac  )*jach
      endif

      !Y flux
      flxjp = 0.5*(v(i,jp,k,2)*arr(i,j,k) + v(i,j,k,2)*arr(i,jp,k))
      flxjm = 0.5*(v(i,jm,k,2)*arr(i,j,k) + v(i,j,k,2)*arr(i,jm,k))

      !Z flux
      flxkp = 0.5*(v(i,j,kp,3)*arr(i,j,k) + v(i,j,k,3)*arr(i,j,kp))
      flxkm = 0.5*(v(i,j,km,3)*arr(i,j,k) + v(i,j,k,3)*arr(i,j,km))

c     Other advective methods

      if (method /= 2 .or. (.not.cv)) then

        !Find velocities at faces from ZIP fluxes
        if (zv) then
          vxp = 2.*flxip/(arr(ip,j,k)+arr(i,j,k)+1d-12)
          vxm = 2.*flxim/(arr(im,j,k)+arr(i,j,k)+1d-12)
                                                    
          vyp = 2.*flxjp/(arr(i,jp,k)+arr(i,j,k)+1d-12)
          vym = 2.*flxjm/(arr(i,jm,k)+arr(i,j,k)+1d-12)
                                                    
          vzp = 2.*flxkp/(arr(i,j,kp)+arr(i,j,k)+1d-12)
          vzm = 2.*flxkm/(arr(i,j,km)+arr(i,j,k)+1d-12)
        else
          if (cv) then
            vxp = 0.5*(v(ip,j,k,1)+v(i,j,k,1))
            vxm = 0.5*(v(im,j,k,1)+v(i,j,k,1))
                                                    
            vyp = 0.5*(v(i,jp,k,2)+v(i,j,k,2))
            vym = 0.5*(v(i,jm,k,2)+v(i,j,k,2))
                                                    
            vzp = 0.5*(v(i,j,kp,3)+v(i,j,k,3))
            vzm = 0.5*(v(i,j,km,3)+v(i,j,k,3))
          else
            vxp = v(i,j,k,1)
            vxm = v(i,j,k,1)
                                                    
            vyp = v(i,j,k,2)
            vym = v(i,j,k,2)
                                                    
            vzp = v(i,j,k,3)
            vzm = v(i,j,k,3)
           endif
        endif

        !Find advected quantity at faces
        if ((i == nx).and.(method > 2 .and. method <= 6)) then
          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
     .             ,grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,vxp,1)
        else
          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
     .             ,grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,vxp,method)
        endif

        if ((i == 1).and.(method > 2 .and. method <= 6)) then
          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,grid_params%xx(igmm),arr(imm,j,k)
     .             ,vxm,1)
        else
          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,grid_params%xx(igmm),arr(imm,j,k)
     .             ,vxm,method)
        endif

        if ((j == ny).and.(method > 2 .and. method <= 6)) then
          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,1)
        else
          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,method)
        endif

        if ((j == 1).and.(method > 2 .and. method <= 6)) then
          qjm=inter(grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
     .             ,vym,1)
        else
          qjm=inter(grid_params%yy(jgp ),arr(i,jp,k )
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm,k )
     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
     .             ,vym,method)
        endif

        if ((k == nz).and.(method > 2 .and. method <= 6)) then
          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
     .             ,grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,vzp,1)
        else
          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
     .             ,grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,vzp,method)
        endif

        if ((k == 1).and.(method > 2 .and. method <= 6)) then
          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,1)
        else
          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,method)
        endif

        !Setup fluxes (ZIP at boundaries by default for conservation)
        if (.not.isBdry(i,igx,2).or.bc_safe.or.(.not.cv)) flxip=vxp*qip
        if (.not.isBdry(i,igx,1).or.bc_safe.or.(.not.cv)) flxim=vxm*qim

        if (.not.isBdry(j,igy,4).or.bc_safe.or.(.not.cv)) flxjp=vyp*qjp
        if (.not.isBdry(j,igy,3).or.bc_safe.or.(.not.cv)) flxjm=vym*qjm

        if (.not.isBdry(k,igz,6).or.bc_safe.or.(.not.cv)) flxkp=vzp*qkp
        if (.not.isBdry(k,igz,5).or.bc_safe.or.(.not.cv)) flxkm=vzm*qkm

      endif

c     Find advective term

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac

      if (vw) advec = advec*gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function flx_advec

ccc     flx_advec
ccc     ###############################################################
cc      function flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
cc     .                     ,method,sp,vol) result(advec)
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
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz,method
cc
cc      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,v  (0:nx+1,0:ny+1,0:nz+1,3),advec
cc
cc      logical,optional :: sp,vol
cc
ccc     Local variables
cc
cc      integer    :: igrid,half_elem,su
cc
cc      integer    :: im,ip,ipp,imm,igm,igp,igpp,igmm
cc     .             ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
cc     .             ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
cc     .             ,ig,jg,kg,igl,jgl,kgl
cc
cc      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
cc      real(8)    :: vxp,vxm,vyp,vym,vzp,vzm
cc     .             ,qip,qim,qjp,qjm,qkp,qkm
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc     .             ,flxhip,flxhim,flxhjp,flxhjm,flxhkp,flxhkm
cc     .             ,flxlip,flxlim,flxljp,flxljm,flxlkp,flxlkm
cc      logical    :: vw,spoint,bc_safe
cc
ccc     Begin program
cc
ccc     Set defaults
cc
cc      bc_safe = (method <=2)
cc
cc      if (PRESENT(vol)) then
cc        vw = vol
cc      else
cc        vw = .true.
cc      endif
cc
cc      if (PRESENT(sp)) then
cc        spoint = sp
cc      else
cc        spoint = .false.
cc      endif
cc
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
ccc     Find low-order (upwind) fluxes
cc
cc      jach = 0.5*(jac+jacip)
cc      flxlip = 0.25*jach
cc     .      *( (    (v(i,j,k,1)/jac+v(ip,j,k,1)/jacip)
cc     .          +abs(v(i,j,k,1)/jac+v(ip,j,k,1)/jacip) ) *phi(i ,j,k)
cc     .        +(    (v(i,j,k,1)/jac+v(ip,j,k,1)/jacip)          
cc     .          -abs(v(i,j,k,1)/jac+v(ip,j,k,1)/jacip) ) *phi(ip,j,k) )
cc
cc      jach = 0.5*(jac+jacim)
cc      flxlim = 0.25*jach
cc     .       *( (    (v(i,j,k,1)/jac+v(im,j,k,1)/jacim)          
cc     .           +abs(v(i,j,k,1)/jac+v(im,j,k,1)/jacim) ) *phi(im,j,k)
cc     .         +(    (v(i,j,k,1)/jac+v(im,j,k,1)/jacim)          
cc     .           -abs(v(i,j,k,1)/jac+v(im,j,k,1)/jacim) ) *phi(i ,j,k) )
cc
cc      flxljp = 0.25*(
cc     .         (    (v(i,j,k,2)+v(i,jp,k,2))
cc     .          +abs(v(i,j,k,2)+v(i,jp,k,2)) )*phi(i,j ,k)
cc     .        +(    (v(i,j,k,2)+v(i,jp,k,2))          
cc     .          -abs(v(i,j,k,2)+v(i,jp,k,2)) )*phi(i,jp,k))
cc
cc      flxljm = 0.25*(
cc     .         (    (v(i,j,k,2)+v(i,jm,k,2))
cc     .          +abs(v(i,j,k,2)+v(i,jm,k,2)) )*phi(i,jm,k)
cc     .        +(    (v(i,j,k,2)+v(i,jm,k,2))          
cc     .          -abs(v(i,j,k,2)+v(i,jm,k,2)) )*phi(i,j ,k))
cc
cc      flxlkp = 0.25*(
cc     .         (    (v(i,j,k,3)+v(i,j,kp,3))
cc     .          +abs(v(i,j,k,3)+v(i,j,kp,3)) )*phi(i,j,k )
cc     .        +(    (v(i,j,k,3)+v(i,j,kp,3))             
cc     .          -abs(v(i,j,k,3)+v(i,j,kp,3)) )*phi(i,j,kp))
cc
cc      flxlkm = 0.25*(
cc     .         (    (v(i,j,k,3)+v(i,j,km,3))
cc     .          +abs(v(i,j,k,3)+v(i,j,km,3)) )*phi(i,j,km)
cc     .        +(    (v(i,j,k,3)+v(i,j,km,3))             
cc     .          -abs(v(i,j,k,3)+v(i,j,km,3)) )*phi(i,j,k ))
cc
ccc     Find high-order (ZIP) fluxes
cc
cc      !X flux
cc      if (.not.isBdry(i,igx,2)) then  !Not at r=1 boundary
cccc      if     (igl < grid_params%nxgl(igx)) then
cc        jach  = 0.5*(jac+jacip)
cc        flxhip = 0.5*(v(ip,j,k,1)*arr(i ,j,k)/jacip
cc     .             + v(i ,j,k,1)*arr(ip,j,k)/jac  )*jach
cc
cc        if (isSP2(i,igx,ibc=1).and.spoint) then
cc          flxhim = 0d0
cc        else
cc          jach  = 0.5*(jac+jacim)
cc          flxhim = 0.5*(v(im,j,k,1)*arr(i ,j,k)/jacim
cc     .               + v(i ,j,k,1)*arr(im,j,k)/jac  )*jach
cc        endif
cccc      elseif (igl == grid_params%nxgl(igx)) then
cc      else
cc        if (v(ip,j,k,1) == -v(i,j,k,1)) then
cc          flxhip = 0d0    !Shut off flux for vn=0 BC
cc        else
cc          flxhip = 0.5*(v(ip,j,k,1)*arr(i,j,k) + v(i,j,k,1)*arr(ip,j,k))
cc        endif
cc
cc        jach  = 0.5*(jac+jacim)
cc        flxhim = 0.5*(v(im,j,k,1)*arr(i ,j,k)/jacim
cc     .             + v(i ,j,k,1)*arr(im,j,k)/jac  )*jach
cc      endif
cc
cc      !Y flux
cc      flxhjp = 0.5*(v(i,jp,k,2)*arr(i,j,k) + v(i,j,k,2)*arr(i,jp,k))
cc      flxhjm = 0.5*(v(i,jm,k,2)*arr(i,j,k) + v(i,j,k,2)*arr(i,jm,k))
cc
cc      !Z flux
cc      flxhkp = 0.5*(v(i,j,kp,3)*arr(i,j,k) + v(i,j,k,3)*arr(i,j,kp))
cc      flxhkm = 0.5*(v(i,j,km,3)*arr(i,j,k) + v(i,j,k,3)*arr(i,j,km))
cc
ccc     Compute flux limiters
cc
cc      rr = (phi(i,j,k)-phi(im,j,k))/(phi(ip,j,k)-phi(i,j,k))
cc      flx_lim_x = flux_lim(rr,method)
cc
cc      rr = (phi(i,j,k)-phi(i,jm,k))/(phi(i,jp,k)-phi(i,j,k))
cc      flx_lim_y = flux_lim(rr,method)
cc
cc      rr = (phi(i,j,k)-phi(i,j,km))/(phi(i,j,kp)-phi(i,j,k))
cc      flx_lim_z = flux_lim(rr,method)
cc
ccc     Find fluxes
cc
ccc     Find velocities at faces from ZIP fluxes
cc
cc      if (method /= 2) then
cc
cccc        vxp = 2.*flxip/(arr(ip,j,k)+arr(i,j,k)+1d-10)
cccc        vxm = 2.*flxim/(arr(im,j,k)+arr(i,j,k)+1d-10)
cccc
cccc        vyp = 2.*flxjp/(arr(i,jp,k)+arr(i,j,k)+1d-10)
cccc        vym = 2.*flxjm/(arr(i,jm,k)+arr(i,j,k)+1d-10)
cccc
cccc        vzp = 2.*flxkp/(arr(i,j,kp)+arr(i,j,k)+1d-10)
cccc        vzm = 2.*flxkm/(arr(i,j,km)+arr(i,j,k)+1d-10)
cc        vxp = 2.*flxip/(arr(ip,j,k)+arr(i,j,k)+1d-12)
cc        vxm = 2.*flxim/(arr(im,j,k)+arr(i,j,k)+1d-12)
cc                                                    
cc        vyp = 2.*flxjp/(arr(i,jp,k)+arr(i,j,k)+1d-12)
cc        vym = 2.*flxjm/(arr(i,jm,k)+arr(i,j,k)+1d-12)
cc                                                    
cc        vzp = 2.*flxkp/(arr(i,j,kp)+arr(i,j,k)+1d-12)
cc        vzm = 2.*flxkm/(arr(i,j,km)+arr(i,j,k)+1d-12)
cc
ccc     Find advected quantity at faces
cc
cc        if ((i == nx).and.(method > 2 .and. method <= 6)) then
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
cc        if ((i == 1).and.(method > 2 .and. method <= 6)) then
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
cc        if ((j == ny).and.(method > 2 .and. method <= 6)) then
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
cc        if ((j == 1).and.(method > 2 .and. method <= 6)) then
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
cc        if ((k == nz).and.(method > 2 .and. method <= 6)) then
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
cc        if ((k == 1).and.(method > 2 .and. method <= 6)) then
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
ccc     Setup fluxes (ZIP at boundaries by default for conservation)
cc
cccc        if (igl<grid_params%nxgl(igx).or.bcond(2)==PER) flxip=vxp*qip
cccc        if (igl>1                    .or.bcond(1)==PER) flxim=vxm*qim
cccc
cccc        if (jgl<grid_params%nygl(igy).or.bcond(4)==PER) flxjp=vyp*qjp
cccc        if (jgl>1                    .or.bcond(3)==PER) flxjm=vym*qjm
cccc
cccc        if (kgl<grid_params%nzgl(igz).or.bcond(6)==PER) flxkp=vzp*qkp
cccc        if (kgl>1                    .or.bcond(5)==PER) flxkm=vzm*qkm
cc        if (.not.isBdry(i,igx,2).or.bc_safe) flxip=vxp*qip
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

c     flx_nc_advec
c     ###############################################################
      function flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz,v,arr
     .                     ,method,sp,vol) result(advec)
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

      logical,optional :: sp,vol

c     Local variables

      integer    :: igrid,half_elem,su

      integer    :: im,ip,ipp,imm,igm,igp,igpp,igmm
     .             ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
     .             ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
     .             ,ig,jg,kg,igl,jgl,kgl

      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jac
      real(8)    :: vxp,vxm,vyp,vym,vzp,vzm
     .             ,qip,qim,qjp,qjm,qkp,qkm
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw,spoint,bc_safe

c     Begin program

c     Set defaults

      bc_safe = (method <=2)

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

c     Find non-conservative fluxes

      !X flux
      flxip = 0.5*(arr(i,j,k) + arr(ip,j,k))*v(i,j,k,1)
      flxim = 0.5*(arr(i,j,k) + arr(im,j,k))*v(i,j,k,1)

      !Y flux
      flxjp = 0.5*(arr(i,j,k) + arr(i,jp,k))*v(i,j,k,2)
      flxjm = 0.5*(arr(i,j,k) + arr(i,jm,k))*v(i,j,k,2)

      !Z flux
      flxkp = 0.5*(arr(i,j,k) + arr(i,j,kp))*v(i,j,k,3)
      flxkm = 0.5*(arr(i,j,k) + arr(i,j,km))*v(i,j,k,3)

c     Find velocities at faces from ZIP fluxes

      if (method /= 2) then

        vxp = v(i,j,k,1)
        vxm = v(i,j,k,1)
                                                    
        vyp = v(i,j,k,2)
        vym = v(i,j,k,2)
                                                    
        vzp = v(i,j,k,3)
        vzm = v(i,j,k,3)

c     Find advected quantity at faces

        if ((i == nx).and.(method > 2 .and. method <= 6)) then
          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
     .             ,grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,vxp,1)
        else
          qip=inter(grid_params%xx(igpp),arr(ipp,j,k)
     .             ,grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,vxp,method)
        endif

        if ((i == 1).and.(method > 2 .and. method <= 6)) then
          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,grid_params%xx(igmm),arr(imm,j,k)
     .             ,vxm,1)
        else
          qim=inter(grid_params%xx(igp ),arr(ip ,j,k)
     .             ,grid_params%xx(ig  ),arr(i  ,j,k)
     .             ,grid_params%xx(igm ),arr(im ,j,k)
     .             ,grid_params%xx(igmm),arr(imm,j,k)
     .             ,vxm,method)
        endif

        if ((j == ny).and.(method > 2 .and. method <= 6)) then
          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,1)
        else
          qjp=inter(grid_params%yy(jgpp),arr(i,jpp,k)
     .             ,grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,vyp,method)
        endif

        if ((j == 1).and.(method > 2 .and. method <= 6)) then
          qjm=inter(grid_params%yy(jgp ),arr(i,jp ,k)
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm ,k)
     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
     .             ,vym,1)
        else
          qjm=inter(grid_params%yy(jgp ),arr(i,jp,k )
     .             ,grid_params%yy(jg  ),arr(i,j  ,k)
     .             ,grid_params%yy(jgm ),arr(i,jm,k )
     .             ,grid_params%yy(jgmm),arr(i,jmm,k)
     .             ,vym,method)
        endif

        if ((k == nz).and.(method > 2 .and. method <= 6)) then
          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
     .             ,grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,vzp,1)
        else
          qkp=inter(grid_params%zz(kgpp),arr(i,j,kpp)
     .             ,grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,vzp,method)
        endif

        if ((k == 1).and.(method > 2 .and. method <= 6)) then
          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,1)
        else
          qkm=inter(grid_params%zz(kgp ),arr(i,j,kp )
     .             ,grid_params%zz(kg  ),arr(i,j,k  )
     .             ,grid_params%zz(kgm ),arr(i,j,km )
     .             ,grid_params%zz(kgmm),arr(i,j,kmm)
     .             ,vzm,method)
        endif

c     Setup fluxes

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

c     tensor_advc
c     #####################################################################
      function tensor_advc(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,v1,method
     .                    ,cons) result(cnv)

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
      logical,optional :: cons

c     Local variables

      integer    :: ieq,hex,hey,hez
      real(8)    :: nabla_v(3,3),jac
      logical    :: conserv

c     Begin program

      adv_meth = method

      if (PRESENT(cons)) then
        conserv = cons
      else
        conserv = .true.
      endif

      if (conserv) then

        vec1 => v0
        vec2 => v1
        cnv = div_tensor(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt__eom()
     .                  ,tnsr_advc_x,tnsr_advc_y,tnsr_advc_z
     .                  ,vol=.false.)
        nullify(vec1,vec2)

      else

        jac = gmetric%grid(igx)%jac(i,j,k)

        select case(method)
        case(1) !Upwinding
          hex = int(sign(1d0,-v0(i,j,k,1)))
          hey = int(sign(1d0,-v0(i,j,k,2)))
          hez = int(sign(1d0,-v0(i,j,k,3)))

          nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                           ,v1(:,:,:,1)
     .                           ,v1(:,:,:,2)
     .                           ,v1(:,:,:,3),hex,hey,hez)
        case(2) !Central differences
          nabla_v = fnabla_v(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                      ,v1(:,:,:,1)
     .                      ,v1(:,:,:,2)
     .                      ,v1(:,:,:,3),0)
        case default
          call pstop('tensor_advc','Advection option not implemented')
        end select

        cnv = matmul(v0(i,j,k,:),nabla_v)/jac

      endif

      end function tensor_advc

c     tnsr_advc_x
c     #############################################################
      subroutine tnsr_advc_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
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
c       * alt_eom: whether to use alternate EOM in singular coord.
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
c     Calculates tensor components t31-t33 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
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

ccc     tnsr_advc_x
ccc     #############################################################
cc      subroutine tnsr_advc_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t11,t12,t13: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine assumes vec1->v (velocity), vec2->p (momentum).
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
cc        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,vis
cc
ccc     Begin program
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
cc     .               +gmetric%grid(igx)%jac (i ,j,k))
cc
cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
cc          jacp = gmetric%grid(igx)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igx)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc
cc        ijac0 = 1d0/jac0
cc        ijacp = 1d0/jacp
cc
ccc       Advective part
cc
cc        t11 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,1)
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,1)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,1)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,1))
cc     .            *ijacp*ijac0*jac**2
cc
cc        t12 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,2)*ijacp
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,2)*ijac0
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)
cc     .            *jac
cc
cc        t13 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,3)
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,3)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,3)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,3))
cc     .            *ijacp*ijac0*jac**2
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac  = 1d0/jac
cc
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine tnsr_advc_x
cc
ccc     tnsr_advc_y
ccc     #############################################################
cc      subroutine tnsr_advc_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t21-t23 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t21,t22,t23: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
cc        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,ijac,vis
cc
ccc     Begin program
cc
cc        jp = j+1
cc        if (flag == 0) jp = j
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
cc     .               +gmetric%grid(igx)%jac (i,j ,k))
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc       Advective part
cc
cc        t21 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,1)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,1)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,1)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,1))
cc
cc        t22 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,2)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,2)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,2)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,2))
cc
cc        t23 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,3)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,3)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,3)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,3))
cc
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine tnsr_advc_y
cc
ccc     tnsr_advc_z
ccc     #############################################################
cc      subroutine tnsr_advc_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t31-t33 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t31,t32,t33: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
cc        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,ijac,vis
cc
ccc     Begin program
cc
cc        kp = k+1
cc        if (flag == 0) kp = k
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
cc     .               +gmetric%grid(igx)%jac (i,j,k ))
cccc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
cccc     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc       Advective part
cc
cc        t31 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,1)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,1)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,1)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,1))
cc
cc        t32 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,2)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,2)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,2)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,2))
cc
cc        t33 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,3)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,3)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,3)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,3))
cc
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine tnsr_advc_z

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

ccc     EOM_advc
ccc     #####################################################################
cc      function EOM_advc(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,dv,cons,upwind)
cc     .                 result(cnv)
cc
ccc     ---------------------------------------------------------------------
ccc     Find advection/diffusion of EOM. Return contravariant components.
ccc     ---------------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz
cc      real(8) :: cnv(3)
cc      real(8),target :: dv  (0:nxx+1,0:nyy+1,0:nzz+1,3)
cc     .                 ,v0  (0:nxx+1,0:nyy+1,0:nzz+1,3)
cc      logical,optional :: upwind,cons
cc
ccc     Local variables
cc
cc      integer    :: ieq,hex,hey,hez
cc      real(8)    :: nabla_v(3,3),jac
cc      logical    :: upwd,conserv
cc
ccc     Begin program
cc
cc      if (PRESENT(upwind)) then
cc        upwd = upwind
cc      else
cc        upwd = .false.
cc      endif
cc
cc      if (PRESENT(cons)) then
cc        conserv = cons
cc      else
cc        conserv = .not.nc_eom_v
cc      endif
cc
cc      if (conserv) then
cc
cc        vec1 => v0
cc        vec2 => dv
cc        cnv = div_tensor(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .                  ,eom_advc_x,eom_advc_y,eom_advc_z
cc     .                  ,vol=.false.)
cc        nullify(vec1,vec2)
cc
cc      else
cc
cc        jac = gmetric%grid(igx)%jac(i,j,k)
cc
cc        if (upwd) then !Upwinding
cccc          hex = 1
cccc          hey = 1
cccc          hez = 1
cccc          if (v0(i,j,k,1) > 0d0) hex = -1
cccc          if (v0(i,j,k,2) > 0d0) hey = -1
cccc          if (v0(i,j,k,3) > 0d0) hez = -1
cc          hex = sign(1d0,-v0(i,j,k,1))
cc          hey = sign(1d0,-v0(i,j,k,2))
cc          hez = sign(1d0,-v0(i,j,k,3))
cc
cc          nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
cc     .                           ,dv(:,:,:,1)
cc     .                           ,dv(:,:,:,2)
cc     .                           ,dv(:,:,:,3),hex,hey,hez)
cc        else
cc          nabla_v = fnabla_v(i,j,k,nxx,nyy,nzz,igx,igy,igz
cc     .                      ,dv(:,:,:,1)
cc     .                      ,dv(:,:,:,2)
cc     .                      ,dv(:,:,:,3),0)
cc        endif
cc
cc        cnv = matmul(v0(i,j,k,:),nabla_v)/jac
cccc        do ieq=1,3
cccc          cnv(ieq) = ( v0(i,j,k,1)*nabla_v(1,ieq)
cccc     .                +v0(i,j,k,2)*nabla_v(2,ieq)
cccc     .                +v0(i,j,k,3)*nabla_v(3,ieq))/jac
cccc        enddo
cc
cc      endif
cc
cc      end function EOM_advc
cc
ccc     eom_advc_x
ccc     #############################################################
cc      subroutine eom_advc_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t11,t12,t13: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine assumes vec1->v (velocity), vec2->p (momentum).
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
cc        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,vis
cc
ccc     Begin program
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
cc     .               +gmetric%grid(igx)%jac (i ,j,k))
cccc        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
cccc     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))
cc
cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
cc          jacp = gmetric%grid(igx)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igx)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc
cc        ijac0 = 1d0/jac0
cc        ijacp = 1d0/jacp
cc
ccc       Advective part
cc
cc        t11 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,1)
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,1)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,1)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,1))
cc     .            *ijacp*ijac0*jac**2
cc
cc        t12 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,2)*ijacp
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,2)*ijac0
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)
cc     .            *jac
cc
cc        t13 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,3)
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,3)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,3)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,3))
cc     .            *ijacp*ijac0*jac**2
cc
cccc        if (flag /= 0 .and. advect /= 2) then
cccc          ipp = min(ip+1,nx+1)
cccc          im  = max(i -1,0   )
cccc
cccc          igp  = ig + (ip  - i)
cccc          igm  = ig + (im  - i)
cccc          igpp = ig + (ipp - i)
cccc        endif
cccc
cccc        !Advection of px
cccc        t11 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
cccc     .             +vec1(i ,j,k,1)*vec2(ip,j,k,1))
cccc     .       *ijacp*ijac0*jac**2
cccc
cccc        ph = 0.5*(vec2(i,j,k,1)+vec2(ip,j,k,1))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t11    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vxp = t11/ph
cccc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc
cccc          if (i == nx .or. i == 0) then
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,1),xx(igp),vec2(ip,j,k,1)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,1),xx(igm),vec2(im,j,k,1)
cccc     .               ,vxp,min(advect,2))
cccc          else
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,1),xx(igp),vec2(ip,j,k,1)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,1),xx(igm),vec2(im,j,k,1)
cccc     .               ,vxp,advect)
cccc          endif
cccc
cccc          t11 = vxp*qip
cccc        endif
cccc
cccc        !Advection of py
cccc        t12 = 0.5*( vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
cccc     .            + vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)*jac
cccc        ph = 0.5*(vec2(i,j,k,2)+vec2(ip,j,k,2))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t12    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vxp = t12/ph
cccc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc
cccc          if (i == nx .or. i == 0) then
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,2),xx(igp),vec2(ip,j,k,2)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,2),xx(igm),vec2(im,j,k,2)
cccc     .               ,vxp,min(advect,2))
cccc          else
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,2),xx(igp),vec2(ip,j,k,2)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,2),xx(igm),vec2(im,j,k,2)
cccc     .               ,vxp,advect)
cccc          endif
cccc
cccc          t12 = vxp*qip
cccc        endif
cccc
cccc        !Advection of pz
cccc        t13 = 0.5*( vec1(i ,j,k,1)*vec2(ip,j,k,3)
cccc     .            + vec1(ip,j,k,1)*vec2(i ,j,k,3))
cccc     .           *ijacp*ijac0*jac**2
cccc        ph = 0.5*(vec2(i,j,k,3)+vec2(ip,j,k,3))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t13    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vxp = t13/ph
cccc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc
cccc          if (i == nx .or. i == 0) then
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,3),xx(igp),vec2(ip,j,k,3)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,3),xx(igm),vec2(im,j,k,3)
cccc     .               ,vxp,min(advect,2))
cccc          else
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,3),xx(igp),vec2(ip,j,k,3)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,3),xx(igm),vec2(im,j,k,3)
cccc     .               ,vxp,advect)
cccc          endif
cccc
cccc          t13 = vxp*qip
cccc        endif
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac  = 1d0/jac
cc
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine eom_advc_x
cc
ccc     eom_advc_y
ccc     #############################################################
cc      subroutine eom_advc_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t21-t23 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t21,t22,t23: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
cc        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,ijac,vis
cc
ccc     Begin program
cc
cc        jp = j+1
cc        if (flag == 0) jp = j
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
cc     .               +gmetric%grid(igx)%jac (i,j ,k))
cccc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
cccc     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc       Advective part
cc
cc        t21 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,1)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,1)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,1)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,1))
cc
cc        t22 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,2)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,2)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,2)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,2))
cc
cc        t23 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,3)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,3)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,3)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,3))
cc
cccc        if (flag /= 0 .and. advect /= 2) then
cccc          jpp = min(jp+1,ny+1)
cccc          jm  = max(j-1 ,0   )
cccc
cccc          jgp  = jg + (jp  - j)
cccc          jgm  = jg + (jm  - j)
cccc          jgpp = jg + (jpp - j)
cccc        endif
cccc
cccc        t21 = 0.5*(vec1(i,j ,k,2)*vec1(i,jp,k,1)
cccc     .            +vec1(i,jp,k,2)*vec1(i,j ,k,1))
cccc
cccc        ph = 0.5*(vec2(i,j,k,1)+vec2(i,jp,k,1))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t21    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vyp = t21/ph
cccc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc
cccc          if (j == ny .or. j == 0) then
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,1),yy(jgp),vec2(i,jp,k,1)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,1),yy(jgm),vec2(i,jm,k,1)
cccc     .               ,vyp,min(advect,2))
cccc          else
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,1),yy(jgp),vec2(i,jp,k,1)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,1),yy(jgm),vec2(i,jm,k,1)
cccc     .               ,vyp,advect)
cccc          endif
cccc
cccc          t21 = vyp*qjp
cccc        endif
cccc
cccc        t22 = 0.5*( vec1(i,j ,k,2)*vec2(i,jp,k,2)
cccc     .            + vec1(i,jp,k,2)*vec2(i,j ,k,2))
cccc
cccc        ph = 0.5*(vec2(i,j,k,2)+vec2(i,jp,k,2))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t22    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vyp = t22/ph
cccc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc
cccc          if (j == ny .or. j == 0) then
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,2),yy(jgp),vec2(i,jp,k,2)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,2),yy(jgm),vec2(i,jm,k,2)
cccc     .               ,vyp,min(advect,2))
cccc          else
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,2),yy(jgp),vec2(i,jp,k,2)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,2),yy(jgm),vec2(i,jm,k,2)
cccc     .               ,vyp,advect)
cccc          endif
cccc
cccc          t22 = vyp*qjp
cccc        endif
cccc
cccc        t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,jp,k,3)
cccc     .            + vec1(i,jp,k,2)*vec2(i,j ,k,3))
cccc
cccc        ph = 0.5*(vec2(i,j,k,3)+vec2(i,jp,k,3))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t23    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vyp = t23/ph
cccc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc
cccc          if (j == ny .or. j == 0) then
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,3),yy(jgp),vec2(i,jp,k,3)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,3),yy(jgm),vec2(i,jm,k,3)
cccc     .               ,vyp,min(advect,2))
cccc          else
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,3),yy(jgp),vec2(i,jp,k,3)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,3),yy(jgm),vec2(i,jm,k,3)
cccc     .               ,vyp,advect)
cccc          endif
cccc
cccc          t23 = vyp*qjp
cccc        endif
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine eom_advc_y
cc
ccc     eom_advc_z
ccc     #############################################################
cc      subroutine eom_advc_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t31-t33 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t31,t32,t33: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
cc        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,ijac,vis
cc
ccc     Begin program
cc
cc        kp = k+1
cc        if (flag == 0) kp = k
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
cc     .               +gmetric%grid(igx)%jac (i,j,k ))
cccc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
cccc     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc       Advective part
cc
cc        t31 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,1)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,1)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,1)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,1))
cc
cc        t32 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,2)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,2)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,2)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,2))
cc
cc        t33 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,3)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,3)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,3)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,3))
cc
cccc        if (flag /= 0 .and. advect /= 2) then
cccc          kpp = min(kp+1,nx+1)
cccc          km  = max(k -1,0   )
cccc
cccc          kgp  = kg + (kp  - k)
cccc          kgm  = kg + (km  - k)
cccc          kgpp = kg + (kpp - k)
cccc        endif
cccc
cccc        t31 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,kp,1)
cccc     .            + vec1(i,j,kp,3)*vec2(i,j,k ,1) )
cccc
cccc        ph = 0.5*(vec2(i,j,k,1)+vec2(i,j,kp,1))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t31    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vzp = t31/ph
cccc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc
cccc          if (k == nz .or. k == 0) then
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,1),zz(kgp),vec2(i,j,kp,1)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,1),zz(kgm),vec2(i,j,km,1)
cccc     .               ,vzp,min(advect,2))
cccc          else
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,1),zz(kgp),vec2(i,j,kp,1)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,1),zz(kgm),vec2(i,j,km,1)
cccc     .               ,vzp,advect)
cccc          endif
cccc
cccc          t31 = vzp*qkp
cccc        endif
cccc
cccc        t32 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,kp,2)
cccc     .            + vec1(i,j,kp,3)*vec2(i,j,k ,2) )
cccc
cccc        ph = 0.5*(vec2(i,j,k,2)+vec2(i,j,kp,2))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t32    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vzp = t32/ph
cccc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc
cccc          if (k == nz .or. k == 0) then
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,2),zz(kgp),vec2(i,j,kp,2)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,2),zz(kgm),vec2(i,j,km,2)
cccc     .               ,vzp,min(advect,2))
cccc          else
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,2),zz(kgp),vec2(i,j,kp,2)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,2),zz(kgm),vec2(i,j,km,2)
cccc     .               ,vzp,advect)
cccc          endif
cccc
cccc          t32 = vzp*qkp
cccc        endif
cccc
cccc        t33 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
cccc     .            + vec1(i,j,k ,3)*vec2(i,j,kp,3) )
cccc
cccc        ph = 0.5*(vec2(i,j,k,3)+vec2(i,j,kp,3))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t33    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vzp = t33/ph
cccc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc
cccc          if (k == nz .or. k == 0) then
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,3),zz(kgp),vec2(i,j,kp,3)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,3),zz(kgm),vec2(i,j,km,3)
cccc     .               ,vzp,min(advect,2))
cccc          else
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,3),zz(kgp),vec2(i,j,kp,3)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,3),zz(kgm),vec2(i,j,km,3)
cccc     .               ,vzp,advect)
cccc          endif
cccc
cccc          t33 = vzp*qkp
cccc        endif
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine eom_advc_z

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
              dvol = gmetric%grid(igx)%dvol(i,j,k)
     .              /gmetric%grid(igx)%jac (i,j,k)

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


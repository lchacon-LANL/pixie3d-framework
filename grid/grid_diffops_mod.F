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

cc      idhx = 0d0
cc      idhy = 0d0
cc      idhz = 0d0
cc
cc      if (vrtx) then
cc        if (nx > 1) idhx = 1./grid_params%dx(ig)
cc        if (ny > 1) idhy = 1./grid_params%dy(jg)
cc        if (nz > 1) idhz = 1./grid_params%dz(kg)
cc      else
cc        if (nx > 1) idhx = 1./grid_params%dxh(ig)
cc        if (ny > 1) idhy = 1./grid_params%dyh(jg)
cc        if (nz > 1) idhz = 1./grid_params%dzh(kg)
cc      endif

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
      function div_vrtx(i,j,k,igx,igy,igz,vx,vy,vz,vrtx) result(div)

c     ------------------------------------------------------------------
c     Calculates div of cnv vector defined at nearest neighbors.
c     If vrtx is true, it finds div at (i+1/2,j+1/2,k+1/2) from (i,j,k)
c     If vrtx is false, it finds div at (i,j,k) from (i+1/2,j+1/2,k+1/2)
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: vx(2,2,2),vy(2,2,2),vz(2,2,2),div
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
      fp = 0.25*(vx(ip,j0,k0)+vx(ip,jp,k0)
     .          +vx(ip,j0,kp)+vx(ip,jp,kp))
      fm = 0.25*(vx(i0,j0,k0)+vx(i0,jp,k0)
     .          +vx(i0,j0,kp)+vx(i0,jp,kp))

      div = (fp-fm)/dhx

      !Y component
      fp = 0.25*(vy(i0,jp,k0)+vy(ip,jp,k0)
     .          +vy(i0,jp,kp)+vy(ip,jp,kp))
      fm = 0.25*(vy(i0,j0,k0)+vy(ip,j0,k0)
     .          +vy(i0,j0,kp)+vy(ip,j0,kp))

      div = div + (fp-fm)/dhy

      !Z component
      fp = 0.25*(vz(ip,j0,kp)+vz(i0,j0,kp)
     .          +vz(ip,jp,kp)+vz(i0,jp,kp))
      fm = 0.25*(vz(ip,j0,k0)+vz(i0,j0,k0)
     .          +vz(ip,jp,k0)+vz(i0,jp,k0))

      div = div + (fp-fm)/dhz

      !Jacobian factor
      if (vrtx) then
cc        jac = 0.125*(gmetric%grid(igx)%jac(i  ,j  ,k  )
cc     .              +gmetric%grid(igx)%jac(i+1,j  ,k  )
cc     .              +gmetric%grid(igx)%jac(i  ,j+1,k  )
cc     .              +gmetric%grid(igx)%jac(i+1,j+1,k  ) 
cc     .              +gmetric%grid(igx)%jac(i  ,j  ,k+1)
cc     .              +gmetric%grid(igx)%jac(i+1,j  ,k+1)
cc     .              +gmetric%grid(igx)%jac(i  ,j+1,k+1)
cc     .              +gmetric%grid(igx)%jac(i+1,j+1,k+1))
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
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,1)
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,2)
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,3)
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

      real(8)    :: cov(3),grad(2,2,2,3)

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .false.
      if (PRESENT(vol)) vol_wgt = vol

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
            grad(il,jl,kl,:) =
     .           matmul(gmetric%grid(igx)%gsup_v(i,j,k,:,:),cov)
          enddo
        enddo
      enddo

c     Divergence at cell centers

      lap = div_vrtx(ii,jj,kk,igx,igy,igz
     .              ,grad(:,:,:,1)
     .              ,grad(:,:,:,2)
     .              ,grad(:,:,:,3)
     .              ,.false.)

      if (vol_wgt) lap=lap*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_vrtx

      end module grid_vertex_ops

c grid_diff_ops
c ######################################################################
      module grid_diff_ops

       use grid_vertex_ops

      end module grid_diff_ops

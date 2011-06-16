c module MK
c ######################################################################
      module mk

        use variables

        use grid

        use io

        use ts_setup

        use mg_solver

        use mg_io

        use imposeBCinterface

        use spline_field

        integer :: sp_order=2             !Spline order for interpolation
     .            ,mk_advect=2            !Grid velocity advection method (see common/grid/grid_diffops_mod.F)

        real(8) :: tau=0d0                !Time constant for MMPDE1
     .            ,g_ratio=1d0            !Cell volume limiting factorx

        logical :: mk_grid=.true.         !Whether we turn on MK mesh motion
     .            ,mk_nc=.true.           !Whether we are doing non-conservative discret.
     .            ,lagged_err_mon=.true.  !Whether we are lagging the monitor function
cc     .            ,limit_err_mon =.false. !Whether we are limiting the monitor function
     .            ,mk_init_relax =.true.  !Whether to relax initial grid according to monitor function

        real(8),target,allocatable,dimension(:,:,:,:) :: gvel

        real(8),target,allocatable,dimension(:,:,:) :: jacn,jacnm

        !MK auxiliary variables
        integer :: IMK=1,mk_bcs(6),mk_map_bcs(6,3)

        real(8),target,allocatable
     .                ,dimension(:,:,:,:) :: xnew,xold,xnm,x0
        real(8),target,allocatable,dimension(:,:,:) :: det,tjac

        type(grid_mg_def),pointer :: mk_gr_params => null()

        type(MG_grid_metrics),pointer :: mk_gmetric => null()

        real(8),private :: mk_det_int

        logical :: mk_res_eval=.true.
     .            ,mk_3d=.true.

        !Spline variables
        real(8),pointer,dimension(:,:,:,:),private :: src_coef
     .                                               ,equ_coef
     .                                               ,x0_coef
     .                                               ,xold_coef

        !Interfaces
        INTERFACE mk_det_3d
          procedure mk_det_op1,mk_det_op2
        END INTERFACE

      contains

c     mk_div_vrtx
c     ##################################################################
      function mk_div_vrtx(i,j,k,igx,igy,igz,v,vrtx) result(div)

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
        dhx = mk_gr_params%dx(ig)
        dhy = mk_gr_params%dy(jg)
        dhz = mk_gr_params%dz(kg)
      else
        dhx = mk_gr_params%dxh(ig)
        dhy = mk_gr_params%dyh(jg)
        dhz = mk_gr_params%dzh(kg)
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
        jac = mk_gmetric%grid(igx)%jac_v(i,j,k)
      else
        jac = mk_gmetric%grid(igx)%jac(i,j,k)
      endif

      div = div/jac

c     End program

      end function mk_div_vrtx

c     mk_grad_vrtx
c     ##################################################################
      function mk_grad_vrtx(i,j,k,igx,igy,igz,phi,vrtx) result(grd)

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

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

      if (vrtx) then
        idhx = 1./mk_gr_params%dx(ig)
        idhy = 1./mk_gr_params%dy(jg)
        idhz = 1./mk_gr_params%dz(kg)
      else
        idhx = 1./mk_gr_params%dxh(ig)
        idhy = 1./mk_gr_params%dyh(jg)
        idhz = 1./mk_gr_params%dzh(kg)
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

      end function mk_grad_vrtx

c     mk_lap
c     ###############################################################
      function mk_lap(i,j,k,nx,ny,nz,igx,igy,igz,arr,enf_bc)
     .         result(laplacian)
c     ---------------------------------------------------------------
c     Calculates J*lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr (0:nx+1,0:ny+1,0:nz+1),laplacian

      logical,optional,intent(IN) :: enf_bc

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
     .             ,flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      real(8)    :: xxold(3),dr(3,3),cov(3,3),grd(3),jac,ijac

      logical    :: vol_wgt,enfbc

c     Begin program

      if (PRESENT(enf_bc)) then
        enfbc = enf_bc
      else
        enfbc = .false.
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

      gsup => mk_gmetric%grid(igrid)%gsup_v

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

      flxip = d_xx_ip  *(arr(ip,j ,k)-arr(i,j ,k))/mk_gr_params%dx(ig)
     .  +0.25*d_xy_ipjp*(arr(ip,jp,k)+arr(i,jp,k)
     .                  -arr(ip,j ,k)-arr(i,j ,k))/mk_gr_params%dy(jg)
     .  +0.25*d_xy_ipjm*(arr(ip,j ,k)+arr(i,j ,k)
     .                  -arr(ip,jm,k)-arr(i,jm,k))/mk_gr_params%dy(jg-1)
     .  +0.25*d_xz_ipkp*(arr(ip,j,kp)+arr(i,j,kp)
     .                  -arr(ip,j,k )-arr(i,j,k ))/mk_gr_params%dz(kg)
     .  +0.25*d_xz_ipkm*(arr(ip,j,k )+arr(i,j,k )
     .                  -arr(ip,j,km)-arr(i,j,km))/mk_gr_params%dz(kg-1)

      flxim =-d_xx_im  *(arr(im,j ,k)-arr(i,j ,k))/mk_gr_params%dx(ig-1)
     .  +0.25*d_xy_imjp*(arr(im,jp,k)+arr(i,jp,k)
     .                  -arr(im,j ,k)-arr(i,j ,k))/mk_gr_params%dy(jg)
     .  +0.25*d_xy_imjm*(arr(im,j ,k)+arr(i,j ,k)
     .                  -arr(im,jm,k)-arr(i,jm,k))/mk_gr_params%dy(jg-1)
     .  +0.25*d_xz_imkp*(arr(im,j,kp)+arr(i,j,kp)
     .                  -arr(im,j,k )-arr(i,j,k ))/mk_gr_params%dz(kg)
     .  +0.25*d_xz_imkm*(arr(im,j,k )+arr(i,j,k )
     .                  -arr(im,j,km)-arr(i,j,km))/mk_gr_params%dz(kg-1)

      flxjp = d_yy_jp  *(arr(i ,jp,k)-arr(i ,j,k))/mk_gr_params%dy(jg)
     .  +0.25*d_xy_ipjp*(arr(ip,jp,k)+arr(ip,j,k)
     .                  -arr(i ,jp,k)-arr(i ,j,k))/mk_gr_params%dx(ig)
     .  +0.25*d_xy_imjp*(arr(i ,jp,k)+arr(i ,j,k)
     .                  -arr(im,jp,k)-arr(im,j,k))/mk_gr_params%dx(ig-1)
     .  +0.25*d_yz_jpkp*(arr(i,jp,kp)+arr(i,j,kp)
     .                  -arr(i,jp,k )-arr(i,j,k ))/mk_gr_params%dz(kg)
     .  +0.25*d_yz_jpkm*(arr(i,jp,k )+arr(i,j,k )
     .                  -arr(i,jp,km)-arr(i,j,km))/mk_gr_params%dz(kg-1)

      flxjm =-d_yy_jm  *(arr(i ,jm,k)-arr(i ,j,k))/mk_gr_params%dy(jg-1)
     .  +0.25*d_xy_ipjm*(arr(ip,jm,k)+arr(ip,j,k)              
     .                  -arr(i ,jm,k)-arr(i ,j,k))/mk_gr_params%dx(ig)
     .  +0.25*d_xy_imjm*(arr(i ,jm,k)+arr(i ,j,k)              
     .                  -arr(im,jm,k)-arr(im,j,k))/mk_gr_params%dx(ig-1)
     .  +0.25*d_yz_jmkp*(arr(i,jm,kp)+arr(i,j,kp)              
     .                  -arr(i,jm,k )-arr(i,j,k ))/mk_gr_params%dz(kg)
     .  +0.25*d_yz_jmkm*(arr(i,jm,k )+arr(i,j,k )              
     .                  -arr(i,jm,km)-arr(i,j,km))/mk_gr_params%dz(kg-1)

      flxkp = d_zz_kp  *(arr(i ,j,kp)-arr(i ,j,k))/mk_gr_params%dz(kg)
     .  +0.25*d_xz_ipkp*(arr(ip,j,kp)+arr(ip,j,k)              
     .                  -arr(i ,j,kp)-arr(i ,j,k))/mk_gr_params%dx(ig)
     .  +0.25*d_xz_imkp*(arr(i ,j,kp)+arr(i ,j,k)              
     .                  -arr(im,j,kp)-arr(im,j,k))/mk_gr_params%dx(ig-1)
     .  +0.25*d_yz_jpkp*(arr(i,jp,kp)+arr(i,jp,k)              
     .                  -arr(i,j ,kp)-arr(i,j ,k))/mk_gr_params%dy(jg)
     .  +0.25*d_yz_jmkp*(arr(i,j ,kp)+arr(i,j ,k)              
     .                  -arr(i,jm,kp)-arr(i,jm,k))/mk_gr_params%dy(jg-1)

      flxkm =-d_zz_km  *(arr(i ,j,km)-arr(i ,j,k))/mk_gr_params%dz(kg)
     .  +0.25*d_xz_ipkm*(arr(ip,j,km)+arr(ip,j,k)              
     .                  -arr(i ,j,km)-arr(i ,j,k))/mk_gr_params%dx(ig)
     .  +0.25*d_xz_imkm*(arr(i ,j,km)+arr(i ,j,k)              
     .                  -arr(im,j,km)-arr(im,j,k))/mk_gr_params%dx(ig-1)
     .  +0.25*d_yz_jpkm*(arr(i,jp,km)+arr(i,jp,k)              
     .                  -arr(i,j ,km)-arr(i,j ,k))/mk_gr_params%dy(jg)
     .  +0.25*d_yz_jmkm*(arr(i,j ,km)+arr(i,j ,k)              
     .                  -arr(i,jm,km)-arr(i,jm,k))/mk_gr_params%dy(jg-1)

      if (enfbc) then
        if ((isBdry(i,igx,1)).and.(mk_bcs(1) /= PER)) flxim=0d0
        if ((isBdry(i,igx,2)).and.(mk_bcs(2) /= PER)) flxip=0d0

        if ((isBdry(j,igy,3)).and.(mk_bcs(3) /= PER)) flxjm=0d0
        if ((isBdry(j,igy,4)).and.(mk_bcs(4) /= PER)) flxjp=0d0

        if ((isBdry(k,igz,5)).and.(mk_bcs(5) /= PER)) flxkm=0d0
        if ((isBdry(k,igz,6)).and.(mk_bcs(6) /= PER)) flxkp=0d0
      endif

      laplacian = (flxip -flxim)/mk_gr_params%dxh(ig)
     .           +(flxjp -flxjm)/mk_gr_params%dyh(jg)
     .           +(flxkp -flxkm)/mk_gr_params%dzh(kg)

c     End program

      end function mk_lap

c     mk_lap_vrtx
c     ###############################################################
      function mk_lap_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,aa,enf_bc)
     .         result(lap)

c     ---------------------------------------------------------------
c     Calculates J*lap(A) at cell centers in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1),lap

      logical,optional,intent(IN) :: enf_bc

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: cov(3),cnv(2,2,2,3)

      logical    :: vol_wgt,enfbc

c     Begin program

      if (PRESENT(enf_bc)) then
        enfbc = enf_bc
      else
        enfbc = .false.
      endif

c     Gradient at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2

            cov = mk_grad_vrtx(ii,jj,kk,igx,igy,igz
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1)
     .                     ,.true.)

            !Transform to contravariant at vertex
            cnv(il,jl,kl,:) =
     .           matmul(mk_gmetric%grid(igx)%gsup_v(ii,jj,kk,:,:),cov)
          enddo
        enddo
      enddo

      if (enfbc) then
        if (isBdry(i,igx,1).and.(mk_bcs(1) /= PER)) cnv(1,:,:,1)=0d0
        if (isBdry(i,igx,2).and.(mk_bcs(2) /= PER)) cnv(2,:,:,1)=0d0

        if (isBdry(j,igy,3).and.(mk_bcs(3) /= PER)) cnv(:,1,:,2)=0d0
        if (isBdry(j,igy,4).and.(mk_bcs(4) /= PER)) cnv(:,2,:,2)=0d0

        if (isBdry(k,igz,5).and.(mk_bcs(5) /= PER)) cnv(:,:,1,3)=0d0
        if (isBdry(k,igz,6).and.(mk_bcs(6) /= PER)) cnv(:,:,2,3)=0d0
      endif

c     Divergence at cell centers

      lap = mk_div_vrtx(i,j,k,igx,igy,igz,cnv,.false.)
     .     *mk_gmetric%grid(igx)%jac(i,j,k)

c     End program

      end function mk_lap_vrtx

c     MK_alloc
c     ###################################################################
      subroutine MK_alloc

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: alloc_stat,ntotd2p

c     Begin program

        allocate (x0   (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,xnew (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,xold (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,xnm  (ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,det  (ilom:ihip,jlom:jhip,klom:khip)
     .           ,tjac (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jacn (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jacnm(ilom:ihip,jlom:jhip,klom:khip)
     .           ,gvel (ilom:ihip,jlom:jhip,klom:khip,3))

c     End program

      end subroutine MK_alloc

c     MK_dealloc
c     ###################################################################
      subroutine MK_dealloc

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat

c     Begin program

        deallocate (x0
     .             ,xold
     .             ,xnm
     .             ,xnew
     .             ,tjac
     .             ,jacn
     .             ,jacnm
     .             ,det
     .             ,gvel,STAT=alloc_stat)

        call killSplines

        call destroyGrid(mk_gr_params)

        deallocate(src_coef,x0_coef,xold_coef,equ_coef
     .            ,STAT=alloc_stat)

c     End program

      end subroutine MK_dealloc

c     MK_initGrid
c     ###############################################################
      subroutine MK_initGrid(igrid,bcs)

c     ---------------------------------------------------------------
c     Initializes grid in unit square.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: igrid
      integer,optional :: bcs(6)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,exp,ugrid,ierr,igx,igy,igz,nx,ny,nz
      real(8)    :: eps,x1,y1,z1,xmin,xmax,ymin,ymax,zmin,zmax
      real(8),allocatable :: xf(:,:,:,:)
      logical    :: is_restart

c     Begin program

c     Hardwire solver configuration

      !Force volume weighing for conservative form
      if (.not.mk_nc) vol_wgt = .true.

      !Fix solver to BDF to ensure proper time centering of metrics
      cnfactor = 0d0
cc      sm_flag  = 2

c     Setup boundary conditions

      if(PRESENT(bcs)) then
        mk_bcs = bcs
        where (mk_bcs /= PER) mk_bcs = NEU
      else
        mk_bcs = NEU !Homogeneous Neumann BCs by default for MK potential
      endif

      !Map BCs (enforces PER bcs at boundary, because it acts on map updates)
      mk_map_bcs(:,1) = mk_bcs
      mk_map_bcs(:,2) = mk_bcs
      mk_map_bcs(:,3) = mk_bcs

      where (mk_map_bcs /= PER) mk_map_bcs = EXT !Extrapolate map by default

      where (mk_map_bcs(1:2,1) /= PER) mk_map_bcs(1:2,1) = DIR
      where (mk_map_bcs(3:4,2) /= PER) mk_map_bcs(3:4,2) = DIR
      where (mk_map_bcs(5:6,3) /= PER) mk_map_bcs(5:6,3) = DIR

c     Use currently provided map

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      x0 = gmetric%grid(igrid)%car

cc      if (nz == 1) then
cc        mk_3d = .false.
cc      else
cc        mk_3d = .true.
cc      endif

c     Define metrics for initial grid

      call defineGridMetric(grid_params,xcar=x0)

      call equateGridStructure(mk_gr_params,grid_params)

      mk_gmetric => mk_gr_params%gmetric

      jacn  = gmetric%grid(igx)%jac
      jacnm = jacn

c     Set up interpolation by splines

      call getMGmap(0,0,0,igx,igy,igz,ig,jg,kg)

      call setupSplines(nx+2,ny+2,nz+2
     .                 ,mk_gr_params%xx(ig:ig+nx+1)
     .                 ,mk_gr_params%yy(jg:jg+ny+1)
     .                 ,mk_gr_params%zz(kg:kg+nz+1),sp_order
     .                 ,xmin,xmax,ymin,ymax,zmin,zmax !-->Outputs logical mesh dimensions
     .                 ,bcnd=mk_bcs)

cc      write (*,*) 'MK domain limits',xmin,xmax,ymin,ymax,zmin,zmax

      allocate(x0_coef  (nx+2,ny+2,nz+2,3))
      allocate(xold_coef(nx+2,ny+2,nz+2,3))

      call splineX(x0,xcoef_ext=x0_coef)

      xnew = x0
      xold = x0
      xnm  = x0

c     End program

      end subroutine MK_initGrid

c     mk_setup_NLF
c     ################################################################
      subroutine mk_setup_NLF(igrid,phi,xold)

cc      use xdraw_io

      implicit none

c     ----------------------------------------------------------------
c     Sets up MK PDE by computing Jacobian of new map
c        xnew = xold + grad(phi)
c     ----------------------------------------------------------------

c     Call variables

      integer :: igrid
      real(8) :: phi(:,:,:),xold(0:,0:,0:,:)

c     Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz
      real(8) :: mag,xnew(0:size(xold,1)-1
     .                   ,0:size(xold,2)-1
     .                   ,0:size(xold,3)-1
     .                   ,size(xold,4))

ccc     Debug
cc
cc      integer :: ieq
cc      real(8),allocatable,dimension(:,:,:) :: debug2
cc      character(13) :: debugfile

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

c     Find vertex positions (xold+grad(phi))

      if (mk_3d) then
       !3D
        do k=0,nz
          do j=0,ny
            do i=0,nx
              xnew(i,j,k,:)=findPos(i,j,k,igx,igy,igz,nx,ny,nz,xold
     .                             ,phi,1)
            enddo
          enddo
        enddo

        !Compute determinant at cell centers
        det = 0d0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              det(i,j,k) = mk_det_3d(i,j,k,nx,ny,nz,igx,igy,igz,xnew)
            enddo
          enddo
        enddo

      else
        !2D
        do k=0,nz
          do j=0,ny
            do i=0,nx
              xnew(i,j,k,:)=findPos(i,j,k,igx,igy,igz,nx,ny,nz,xold
cc     .                             ,phi,1,abs=.false.)
     .                             ,phi,1)
            enddo
          enddo
        enddo

        !Compute Poisson bracket at cell centers
        det = 0d0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              det(i,j,k) = mk_det_2d(i,j,k,nx,ny,nz,igx,igy,igz
     .                             ,xnew(:,:,:,1),xnew(:,:,:,2))
cc     .            +mk_gmetric%grid(igx)%jac(i,j,k)
cc     .            +mk_lap(i,j,k,nx,ny,nz,igx,igy,igz,phi,enf_bc=.true.)
            enddo
          enddo
        enddo

      endif

c diag ****
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      k = 1
cc
cc      nx = grid_params%nxv(igx)
cc      ny = grid_params%nyv(igx)
cc      nz = grid_params%nzv(igx)
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      debug2 = det
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,0d0,1d0,0d0,1d0,0,110)
cccc      debug2 = gmetric%grid(igx)%jac
cccc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cccc     .              ,0d0,1d0,0d0,1d0,0,110)
cccc
cccc      debug2 = gmetric%grid(igx)%gsup(:,:,:,1,1)
cccc      call contour(debug2(1:nx,1:ny,k),nx,ny
cccc     .              ,0d0,1d0,0d0,1d0,1,110)
cccc
cccc      debug2 = gmetric%grid(igx)%gsup(:,:,:,2,2)
cccc      call contour(debug2(1:nx,1:ny,k),nx,ny
cccc     .              ,0d0,1d0,0d0,1d0,1,110)
cccc
cccc      debug2 = gmetric%grid(igx)%gsup(:,:,:,2,1)
cccc      call contour(debug2(1:nx,1:ny,k),nx,ny
cccc     .              ,0d0,1d0,0d0,1d0,1,110)
cccc
cccc      do ieq=1,3
cccc        debug2 = xold(0:nx+1,0:ny+1,0:nz+1,ieq)
cccc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cccc     .              ,0d0,1d0,0d0,1d0,1,110)
cccc      enddo
cccc
cccc      do ieq=1,3
cccc        debug2 = xnew(0:nx+1,0:ny+1,0:nz+1,ieq)
cccc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cccc     .              ,0d0,1d0,0d0,1d0,1,110)
cccc      enddo
cccc
cc      deallocate(debug2)
cc      close(110)
cc
cc#if defined(petsc)
cc      call PetscEnd(mpierr)
cc#endif
cc      stop
c diag ****

      end subroutine mk_setup_NLF

c     MK_residual
c     ###############################################################
      function MK_residual(i,j,k,igrid,phi) result(mk_res)

c     ---------------------------------------------------------------
c     Computes MK nonlinear residual
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igrid

      real(8) :: phi(0:,0:,0:),mk_res

c     Local variables

      integer :: nx,ny,nz,igx,igy,igz

c     Begin program
 
      if (.not.mk_res_eval) then
        mk_res = 0d0
        return
      endif

      igx = igrid
      igy = igrid
      igz = igrid

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

      mk_res = det(i,j,k)
     .        -mk_lap_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,phi,enf_bc=.true.)
     .        +mk_lap     (i,j,k,nx,ny,nz,igx,igy,igz,phi,enf_bc=.true.)

      mk_res = -(mk_res - tjac(i,j,k))

c     End program

      end function MK_residual

c     mk_det_op1
c     ###############################################################
      function mk_det_op1(i,j,k,nx,ny,nz,igx,igy,igz,map_v)
     .         result(mk_det_op)
c     ---------------------------------------------------------------
c     Calculates triple product [map_v1,map_v2,map_v3]
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: map_v(0:nx+1,0:ny+1,0:nz+1,3),mk_det_op

c     Local variables

      integer :: l

      real(8) :: gr(3,3)

c     Begin program

c     Find gradients at cell center from vertices

      do l=1,3
        gr(:,l) = grad_vrtx(i,j,k,igx,igy,igz
     .                     ,map_v(i-1:i,j-1:j,k-1:k,l),.false.)
      enddo

c     Find triple product

      mk_det_op = determ3(gr)

c     End program

      end function mk_det_op1

c     mk_det_op2
c     ###############################################################
      function mk_det_op2(i,j,k,nx,ny,nz,igx,igy,igz,mapx,mapy,mapz)
     .         result(mk_det_op)
c     ---------------------------------------------------------------
c     Calculates triple product [map_v1,map_v2,map_v3]
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8) :: mapx(0:nx+1,0:ny+1,0:nz+1)
     .          ,mapy(0:nx+1,0:ny+1,0:nz+1)
     .          ,mapz(0:nx+1,0:ny+1,0:nz+1),mk_det_op

c     Local variables

      integer :: im,jm,km,l

      real(8) :: gr(3,3)

c     Begin program

      im = i-1
      jm = j-1
      km = k-1

c     Find gradients at cell center from vertices

      gr(:,1) = grad_vrtx(i,j,k,igx,igy,igz,mapx(im:i,jm:j,km:k)
     .                   ,.false.)
      gr(:,2) = grad_vrtx(i,j,k,igx,igy,igz,mapy(im:i,jm:j,km:k)
     .                   ,.false.)
      gr(:,3) = grad_vrtx(i,j,k,igx,igy,igz,mapz(im:i,jm:j,km:k)
     .                   ,.false.)

c     Find triple product

      mk_det_op = determ3(gr)

c     End program

      end function mk_det_op2

c     mk_det_2d
c     ###############################################################
      function mk_det_2d(i,j,k,nx,ny,nz,igx,igy,igz,grdx,grdy)

c     ---------------------------------------------------------------
c     Calculates the elliptic operator.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: grdx(0:nx+1,0:ny+1,0:nz+1)
     .             ,grdy(0:nx+1,0:ny+1,0:nz+1),mk_det_2d

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: dxx,dyy,dxy,dyx

c     Begin program

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

c     Fing grid metrics

      dxx = (grdx(i ,jm,k)+grdx(i ,j,k)
     .      -grdx(im,jm,k)-grdx(im,j,k))/2d0/mk_gr_params%dxh(ig)
      dyy = (grdy(im,j ,k)+grdy(i,j ,k)
     .      -grdy(im,jm,k)-grdy(i,jm,k))/2d0/mk_gr_params%dyh(jg)
      dxy = (grdy(i ,jm,k)+grdy(i ,j,k)
     .      -grdy(im,jm,k)-grdy(im,j,k))/2d0/mk_gr_params%dxh(ig)
      dyx = (grdx(im,j ,k)+grdx(i,j ,k)
     .      -grdx(im,jm,k)-grdx(i,jm,k))/2d0/mk_gr_params%dyh(jg)

      mk_det_2d = (dxx*dyy - dxy*dyx)

c     End program

      end function mk_det_2d

c     MK_get_mon
c     ###############################################################
      subroutine MK_get_mon(igrid)

c     ---------------------------------------------------------------
c     Updatess target Jacobian, jac', on current mesh.
c     ---------------------------------------------------------------

      implicit none             !For safe fortran

c     Call variables

      integer :: igrid

c     Local variables

      integer :: nx,ny,nz
      real(8),allocatable,dimension(:,:,:)  :: zeros
      logical :: fgc

      external app_mon_fun

c     Begin program

      nx = mk_gr_params%nxv(igrid)
      ny = mk_gr_params%nyv(igrid)
      nz = mk_gr_params%nzv(igrid)

      if (lagged_err_mon) then
        if (old_time_eval.or.source_eval) then
          !Set target jacobian for old mesh x'=x+grad(phi_n)
          call app_mon_func(igrid,nx,ny,nz,tjac)  !External

          call splineFld(tjac)  !Spline jacobian on xold mesh
        else
          !Interpolate old-time target jacobian on current mesh
          call mk_interp_J(igrid,nx,ny,nz,tjac,fill_gc=.true.)
cc          call app_mon_func(igrid,nx,ny,nz,tjac) !External
        endif
      else
        call app_mon_func(igrid,nx,ny,nz,tjac)   !External
      endif

      call mk_norm_J(igrid,nx,ny,nz,tjac,int_det=.true.)

c     End program

      end subroutine MK_get_mon

c     mk_interp_J
c     ###############################################################
      subroutine mk_interp_J(igrid,nx,ny,nz,jac,fill_gc)

c     ---------------------------------------------------------------
c     Interpolates target Jacobian jac' at new map x=xold+grad(phi).
c     ---------------------------------------------------------------

      implicit none             !For safe fortran

c     Call variables

      integer :: nx,ny,nz,igrid
      real(8) :: jac(0:nx+1,0:ny+1,0:nz+1)
      logical,optional :: fill_gc

c     Local variables

      integer :: igx,igy,igz,i,j,k,ig,jg,kg,bcs(6)
      real(8) :: dummy(0:nx+1,0:ny+1,0:nz+1),x1,x2,x3
      logical :: fgc

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      if (PRESENT(fill_gc)) then
        fgc = fill_gc
      else
        fgc = .false.
      endif

      jac = 0d0
      dummy = 0d0

      if ((igrid > 1).and.fgc) then
        call pstop('MK_def_targetj'
     .            ,'Cannot fill ghost cells in coarse grids')
      endif

      do k=1,nz
        do j=1,ny
          do i=1,nx
            call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
            jac(i,j,k) = evalJac(i,j,k,igrid,x1,x2,x3)
cc            jac(i,j,k) = 
cc     .           mk_interp_J_ijk(i,j,k,igx,igy,igz,nx,ny,nz,xold,phi)
          enddo
        enddo
      enddo

      !Impose *physics* topological/logical boundary conditions
      if (fgc) then
        bcs = bcond
        where (bcond == DEF) bcs = EXT
        call setBC(IMK,nx,ny,nz,jac,dummy,bcs,igx,igy,igz
cc     .            ,iorder=sp_order)
     .            ,iorder=0)
      endif

c     End program

      end subroutine mk_interp_J

ccc     mk_interp_J_ijk
ccc     ##############################################################
cc      function mk_interp_J_ijk(i,j,k,igx,igy,igz,nxx,nyy,nzz,xold,phi
cc     .                    ,smooth)
cc     .         result(jac)
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,igx,igy,igz,nxx,nyy,nzz
cc      real(8)    :: jac
cc     .             ,phi  (0:nxx+1,0:nyy+1,0:nzz+1)
cc     .             ,xold (0:nxx+1,0:nyy+1,0:nzz+1,3)
cc      logical,optional :: smooth
cc
ccc     Local variables
cc
cc      integer    :: npts,nptsx,nptsy,nptsz
cc     .             ,ig,jg,kg,igp,igm,jgp,jgm,kgp,kgm
cc     .             ,iloc,jloc,kloc,im,ip,jm,jp,km,kp
cc      real(8)    :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
cc     .             ,xi1,xi2,xi3,dlx,dly,dlz,x1,x2,x3,pos(3),dv
cc     .             ,x1ipjpkp,x1imjpkp,x1ipjmkp,x1imjmkp
cc     .             ,x2ipjpkp,x2imjpkp,x2ipjmkp,x2imjmkp
cc     .             ,x3ipjpkp,x3imjpkp,x3ipjmkp,x3imjmkp
cc     .             ,x1ipjpkm,x1imjpkm,x1ipjmkm,x1imjmkm
cc     .             ,x2ipjpkm,x2imjpkm,x2ipjmkm,x2imjmkm
cc     .             ,x3ipjpkm,x3imjpkm,x3ipjmkm,x3imjmkm
cc
cc      logical    :: smth
cc
ccc     Begin program
cc
cc      if (PRESENT(smooth)) then
cc        smth = smooth
cc      else
cc        smth =     (i >= 1).and.(i <= nxx)
cc     .        .and.(j >= 1).and.(j <= nyy)
cc     .        .and.(k >= 1).and.(k <= nzz)
cc      endif
cc
cc      ip=i+1
cc      im=i-1
cc      jp=j+1
cc      jm=j-1
cc      kp=k+1
cc      km=k-1
cc
ccc     Number of points according to grid position
cc
cc      if (smth) then
cc        npts = 3*igx
cccc        npts = igx
cc      else
cc        npts = igx
cc      endif
cc
cc      npts = 1
cc
ccc     Define adaptivity measure
cc
cc      if (npts == 1 .and. (.not.smth)) then
cc
cc        pos = findPos(i,j,k,igx,igy,igz,nxx,nyy,nzz,xold,phi,0)
cc        jac = evalJac(i,j,k,igx,pos(1),pos(2),pos(3))
cc
cc      else
cc
cc        call getMGmap(i ,j,k,igx,igy,igz,ig ,jg,kg)
cc        call getMGmap(ip,j,k,igx,igy,igz,igp,jg,kg)
cc        call getMGmap(im,j,k,igx,igy,igz,igm,jg,kg)
cc
cc        call getMGmap(i,j ,k,igx,igy,igz,ig,jg ,kg)
cc        call getMGmap(i,jp,k,igx,igy,igz,ig,jgp,kg)
cc        call getMGmap(i,jm,k,igx,igy,igz,ig,jgm,kg)
cc
cc        call getMGmap(i,j,k ,igx,igy,igz,ig,jg,kg )
cc        call getMGmap(i,j,kp,igx,igy,igz,ig,jg,kgp)
cc        call getMGmap(i,j,km,igx,igy,igz,ig,jg,kgm)
cc
cc        !Define integration domain
cc        lxmin = 0.5*(mk_gr_params%xx(igm)+mk_gr_params%xx(ig))
cc        lxmax = 0.5*(mk_gr_params%xx(igp)+mk_gr_params%xx(ig))
cc        lymin = 0.5*(mk_gr_params%yy(jgm)+mk_gr_params%yy(jg))
cc        lymax = 0.5*(mk_gr_params%yy(jgp)+mk_gr_params%yy(jg))
cc        lzmin = 0.5*(mk_gr_params%zz(kgm)+mk_gr_params%zz(kg))
cc        lzmax = 0.5*(mk_gr_params%zz(kgp)+mk_gr_params%zz(kg))
cc
cc        !Define integration subintervals
cc        if (nxx == 1) then
cc          nptsx = 1
cc        else
cc          nptsx = npts
cc        endif
cc
cc        if (nyy == 1) then
cc          nptsy = 1
cc        else
cc          nptsy = npts
cc        endif
cc
cc        if (nzz == 1) then
cc          nptsz = 1
cc        else
cc          nptsz = npts
cc        endif
cc
cc        dlx = (lxmax-lxmin)/nptsx
cc        dly = (lymax-lymin)/nptsy
cc        dlz = (lzmax-lzmin)/nptsz
cc
cc        !Limits in Cartesian space
cc        pos=findPos(i ,j ,k,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1ipjpkp=pos(1)
cc        x2ipjpkp=pos(2)
cc        x3ipjpkp=pos(3)
cc
cc        pos=findPos(i ,jm,k,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1ipjmkp=pos(1)
cc        x2ipjmkp=pos(2)
cc        x3ipjmkp=pos(3)
cc
cc        pos=findPos(im,j ,k,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1imjpkp=pos(1)
cc        x2imjpkp=pos(2)
cc        x3imjpkp=pos(3)
cc
cc        pos=findPos(im,jm,k,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1imjmkp=pos(1)
cc        x2imjmkp=pos(2)
cc        x3imjmkp=pos(3)
cc
cc        pos=findPos(i ,j ,km,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1ipjpkm=pos(1)
cc        x2ipjpkm=pos(2)
cc        x3ipjpkm=pos(3)
cc
cc        pos=findPos(i ,jm,km,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1ipjmkm=pos(1)
cc        x2ipjmkm=pos(2)
cc        x3ipjmkm=pos(3)
cc
cc        pos=findPos(im,j ,km,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1imjpkm=pos(1)
cc        x2imjpkm=pos(2)
cc        x3imjpkm=pos(3)
cc
cc        pos=findPos(im,jm,km,igx,igy,igz,nxx,nyy,nzz,xold,phi,1)
cc        x1imjmkm=pos(1)
cc        x2imjmkm=pos(2)
cc        x3imjmkm=pos(3)
cc
cc        !Perform integral in logical space
cc        jac = 0d0
cc        dv  = 0d0
cc        do kloc=1,nptsz
cc          do jloc=1,nptsy
cc            do iloc=1,nptsx
cc            !Logical coordinates
cc            xi1 = lxmin-0.5*dlx + iloc*dlx
cc            xi2 = lymin-0.5*dly + jloc*dly
cc            xi3 = lzmin-0.5*dlz + kloc*dlz
cc            !3D trilinear interpolation
cc            x1 = x1ipjpkp*(xi1  -lxmin)*(xi2  -lymin)*(xi3  -lzmin)
cc     .                  /((lxmax-lxmin)*(lymax-lymin)*(lzmax-lzmin))
cc     .          +x1ipjmkp*(xi1  -lxmin)*(xi2  -lymax)*(xi3  -lzmin)
cc     .                  /((lxmax-lxmin)*(lymin-lymax)*(lzmax-lzmin))
cc     .          +x1imjpkp*(xi1  -lxmax)*(xi2  -lymin)*(xi3  -lzmin)
cc     .                  /((lxmin-lxmax)*(lymax-lymin)*(lzmax-lzmin))
cc     .          +x1imjmkp*(xi1  -lxmax)*(xi2  -lymax)*(xi3  -lzmin)
cc     .                  /((lxmin-lxmax)*(lymin-lymax)*(lzmax-lzmin))
cc     .          +x1ipjpkm*(xi1  -lxmin)*(xi2  -lymin)*(xi3  -lzmax)
cc     .                  /((lxmax-lxmin)*(lymax-lymin)*(lzmin-lzmax))
cc     .          +x1ipjmkm*(xi1  -lxmin)*(xi2  -lymax)*(xi3  -lzmax)
cc     .                  /((lxmax-lxmin)*(lymin-lymax)*(lzmin-lzmax))
cc     .          +x1imjpkm*(xi1  -lxmax)*(xi2  -lymin)*(xi3  -lzmax)
cc     .                  /((lxmin-lxmax)*(lymax-lymin)*(lzmin-lzmax))
cc     .          +x1imjmkm*(xi1  -lxmax)*(xi2  -lymax)*(xi3  -lzmax)
cc     .                  /((lxmin-lxmax)*(lymin-lymax)*(lzmin-lzmax))
cc
cc            x2 = x2ipjpkp*(xi1  -lxmin)*(xi2  -lymin)*(xi3  -lzmin)
cc     .                  /((lxmax-lxmin)*(lymax-lymin)*(lzmax-lzmin))
cc     .          +x2ipjmkp*(xi1  -lxmin)*(xi2  -lymax)*(xi3  -lzmin)
cc     .                  /((lxmax-lxmin)*(lymin-lymax)*(lzmax-lzmin))
cc     .          +x2imjpkp*(xi1  -lxmax)*(xi2  -lymin)*(xi3  -lzmin)
cc     .                  /((lxmin-lxmax)*(lymax-lymin)*(lzmax-lzmin))
cc     .          +x2imjmkp*(xi1  -lxmax)*(xi2  -lymax)*(xi3  -lzmin)
cc     .                  /((lxmin-lxmax)*(lymin-lymax)*(lzmax-lzmin))
cc     .          +x2ipjpkm*(xi1  -lxmin)*(xi2  -lymin)*(xi3  -lzmax)
cc     .                  /((lxmax-lxmin)*(lymax-lymin)*(lzmin-lzmax))
cc     .          +x2ipjmkm*(xi1  -lxmin)*(xi2  -lymax)*(xi3  -lzmax)
cc     .                  /((lxmax-lxmin)*(lymin-lymax)*(lzmin-lzmax))
cc     .          +x2imjpkm*(xi1  -lxmax)*(xi2  -lymin)*(xi3  -lzmax)
cc     .                  /((lxmin-lxmax)*(lymax-lymin)*(lzmin-lzmax))
cc     .          +x2imjmkm*(xi1  -lxmax)*(xi2  -lymax)*(xi3  -lzmax)
cc     .                  /((lxmin-lxmax)*(lymin-lymax)*(lzmin-lzmax))
cc
cc            x3 = x3ipjpkp*(xi1  -lxmin)*(xi2  -lymin)*(xi3  -lzmin)
cc     .                  /((lxmax-lxmin)*(lymax-lymin)*(lzmax-lzmin))
cc     .          +x3ipjmkp*(xi1  -lxmin)*(xi2  -lymax)*(xi3  -lzmin)
cc     .                  /((lxmax-lxmin)*(lymin-lymax)*(lzmax-lzmin))
cc     .          +x3imjpkp*(xi1  -lxmax)*(xi2  -lymin)*(xi3  -lzmin)
cc     .                  /((lxmin-lxmax)*(lymax-lymin)*(lzmax-lzmin))
cc     .          +x3imjmkp*(xi1  -lxmax)*(xi2  -lymax)*(xi3  -lzmin)
cc     .                  /((lxmin-lxmax)*(lymin-lymax)*(lzmax-lzmin))
cc     .          +x3ipjpkm*(xi1  -lxmin)*(xi2  -lymin)*(xi3  -lzmax)
cc     .                  /((lxmax-lxmin)*(lymax-lymin)*(lzmin-lzmax))
cc     .          +x3ipjmkm*(xi1  -lxmin)*(xi2  -lymax)*(xi3  -lzmax)
cc     .                  /((lxmax-lxmin)*(lymin-lymax)*(lzmin-lzmax))
cc     .          +x3imjpkm*(xi1  -lxmax)*(xi2  -lymin)*(xi3  -lzmax)
cc     .                  /((lxmin-lxmax)*(lymax-lymin)*(lzmin-lzmax))
cc     .          +x3imjmkm*(xi1  -lxmax)*(xi2  -lymax)*(xi3  -lzmax)
cc     .                  /((lxmin-lxmax)*(lymin-lymax)*(lzmin-lzmax))
cc
cc            !Function evaluation
cc            jac = jac + dlx*dly*dlz*evalJac(i,j,k,igx,x1,x2,x3)
cc            dv  = dv  + dlx*dly*dlz
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        jac = jac/dv
cc      endif
cc
cc      end function mk_interp_J_ijk

c     evalJac
c     ##############################################################
      function evalJac(i,j,k,igrid,x1,x2,x3) result(jac)

c     --------------------------------------------------------------
c     Evaluate jacobian at physical coordinates x1,x2,x3
c     --------------------------------------------------------------
        implicit none

c     Call variables

        integer :: i,j,k,igrid
        real(8) :: jac,x1,x2,x3
      
c     Local variables

        integer :: ilevel,ierror,ig,jg,kg
        real(8) :: xi1,xi2,xi3

c     Begin program

c     Initial guess for logical coords

        call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid,ig,jg,kg
     .                                ,xi1,xi2,xi3)

c     Invert old map to find logical coords (this function only called when tjac is lagged)

        ilevel = 0 !For debugging level in evalXi
        call evalXi(ilevel,x1,x2,x3,xi1,xi2,xi3,ierror
     .             ,xcoef_ext=xold_coef)!,no_perbc=.true.)

        if (ierror /= 0) then
          call pstop('evalJac in MK','Map inversion failed')
        endif

c     Check for error

        jac = evalFld(xi1,xi2,xi3,ierror)

      end function evalJac

c     mk_norm_J
c     ###############################################################
      subroutine mk_norm_J(igrid,nx,ny,nz,jac,int_det)

c     ---------------------------------------------------------------
c     Normalizes target Jacobian jac'.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8) :: jac(0:nx+1,0:ny+1,0:nz+1)
      logical,optional :: int_det

c     Local variables

      integer :: igx,igy,igz,i,j,k,ig,jg,kg
      real(8) :: cc_new,sendbuf(1),recbuf(1),jac_norm
      logical :: fgc,int__det

      save jac_norm

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      if (PRESENT(int_det)) then
        int__det = int_det
      else
        int__det = .false.
      endif

c     Find integral of det(d x'/dx)

      if (int__det) then
        mk_det_int = integral(igx,igy,igz,nx,ny,nz,det
     .                       ,average=.false.
     .                       ,logical_dom=.true.)
      endif

c     Rescale target jacobian according to max-to-min ratio

      call mk_limit_J(jac)

c     Find integral of error monitor function, Int(J' dxi)

      cc_new = integral(igx,igy,igz,nx,ny,nz,jac
     .                 ,average=.false.
     .                 ,logical_dom=.true.)

      jac_norm = cc_new/mk_det_int  !Int(J' dxi)/Int(dxi)

      jac = jac/jac_norm   !jac_norm defined in module

c     End program

      end subroutine mk_norm_J

c     mk_limit_J
c     ##############################################################
      subroutine mk_limit_J(mon)

c     --------------------------------------------------------------
c     Sets max-to-min cell volume ratio.
c     --------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: mon(:,:,:)

c     Local variables

      real(8) :: maxv,minv,alpha,gext(2),lext(2)

c     Begin program

      if (g_ratio > 1d0) then
        maxv = maxval(mon)
        minv = minval(mon)

#if defined(petsc)
        lext = (/ maxv, minv/)
        call MPI_AllReduce(lext,gext,2,MPI_DOUBLE_PRECISION
     .                    ,MPI_MAX,MPI_COMM_WORLD,mpierr)
        maxv = gext(1)
        minv = gext(2)
#endif

        if (maxv > g_ratio*minv) then
          alpha = (maxv - g_ratio*minv)/(g_ratio - 1d0)
          mon = mon + alpha
        endif
      endif

      end subroutine mk_limit_J

c     findPos
c     ##############################################################
      function findPos(i,j,k,igx,igy,igz,nxx,nyy,nzz,xold,phi,he,abs)
     .         result(pos)

c     --------------------------------------------------------------
c     Finds grid vertex and cell-center CARTESIAN positions from phi.
c     On input:
c        * xold: old Cartesian cell-center positions
c        * he: determines whether we find vertex positions (he=1)
c              or cell center positions (he=0). If vertex, computes
c              position at (i+1/2,j+1/2,k+1/2).
c        * abs: determines whether we find relative (grad(phi)) or
c               absolute (xold+grad(phi)) coordinates
c     --------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igx,igy,igz,nxx,nyy,nzz,he
      real(8)    :: phi (0:nxx+1,0:nyy+1,0:nzz+1),pos(3)
     .             ,xold(0:nxx+1,0:nyy+1,0:nzz+1,3)
      logical,optional :: abs

c     Local variables

      integer    :: ig,jg,kg,hex,hey,hez
      real(8)    :: xxold(3),dr(3,3),cov(3,3),grd(3)
      logical    :: absolute

c     Begin program

      if (PRESENT(abs)) then
        absolute = abs
      else
        absolute = .true.
      endif

      select case(he)
      case(1)
        grd = grad_vrtx(i,j,k,igx,igy,igz
     .                 ,phi(i:i+1,j:j+1,k:k+1),.true.)
      
        cov = mk_gmetric%grid(igx)%cov_v(i,j,k,:,:)
        grd = matmul(transpose(cov),grd)

        !Enforce n.grad(phi) = 0 at boundaries
        if ((isBdry(i+1,igx,1).or.isBdry(i,igx,2))
     .      .and.(mk_bcs(1) /= PER)) grd(1) = 0d0
        if ((isBdry(j+1,igy,3).or.isBdry(j,igy,4))
     .      .and.(mk_bcs(3) /= PER)) grd(2) = 0d0
        if ((isBdry(k+1,igz,5).or.isBdry(k,igz,6))
     .      .and.(mk_bcs(5) /= PER)) grd(3) = 0d0

        if (absolute) then
          xxold = 0.125*(xold(i+1,j,k  ,:)+xold(i+1,j+1,k  ,:)
     .                  +xold(i  ,j,k  ,:)+xold(i  ,j+1,k  ,:)
     .                  +xold(i+1,j,k+1,:)+xold(i+1,j+1,k+1,:)
     .                  +xold(i  ,j,k+1,:)+xold(i  ,j+1,k+1,:))
        endif

      case default
        grd = grad(i,j,k,nxx,nyy,nzz,igx,igy,igz,phi)

        cov = mk_gmetric%grid(igx)%cov(i,j,k,:,:)
        grd = matmul(transpose(cov),grd)
      
        if (absolute) xxold = xold(i,j,k,:)
      end select

      if (absolute) then
        pos = xxold + grd
      else
        pos = grd
      endif

      end function findPos

ccc     mk_xfer_mesh
ccc     ##############################################################
cc      subroutine mk_xfer_mesh(igrid,nx,ny,nz,phi,xold,xnew)
ccc     --------------------------------------------------------------
ccc     Transfers x'=x+grad(phi) -> x at cell centers
ccc     --------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer :: igrid,nx,ny,nz
cc      real(8),INTENT(IN) :: phi (0:nx+1,0:ny+1,0:nz+1)
cc     .                     ,xold(0:nx+1,0:ny+1,0:nz+1,3)
cc      real(8) :: xnew(0:nx+1,0:ny+1,0:nz+1,3)
cc
ccc     Local variables
cc
cc      integer :: i,j,k,igx,igy,igz
cc      real(8) :: xx0 (0:nx+1,0:ny+1,0:nz+1,3)
cc
cc      real(8) :: sendbuf(1),recbuf(1)
cc
ccc     Begin program
cc
cc      igx = igrid
cc      igy = igrid
cc      igz = igrid
cc
ccc     Find vertex positions and average to cell centers
cc
cc      xx0 = 0d0
cc
cc      !Find vertex positions
cc      do k=0,nz
cc        do j=0,ny
cc          do i=0,nx
cc            xx0(i,j,k,:)=findPos(i,j,k,igx,igy,igz,nx,ny,nz,xold
cc     .                           ,phi,1,abs=.false.)
cc          enddo
cc        enddo
cc      enddo
cc
cc      !Average to cell centers
cc      xnew(1:nx,1:ny,1:nz,:) = xold(1:nx,1:ny,1:nz,:)
cc     .  +0.125*(xx0(0:nx-1,1:ny,0:nz-1,:)+xx0(0:nx-1,0:ny-1,0:nz-1,:)
cc     .         +xx0(1:nx  ,1:ny,0:nz-1,:)+xx0(1:nx  ,0:ny-1,0:nz-1,:)
cc     .         +xx0(0:nx-1,1:ny,1:nz  ,:)+xx0(0:nx-1,0:ny-1,1:nz  ,:)
cc     .         +xx0(1:nx  ,1:ny,1:nz  ,:)+xx0(1:nx  ,0:ny-1,1:nz  ,:))
cc
ccc     Prepare to fill ghost cells by extrapolation
cc
cccc      call map__bc(1)
cccc      call map__bc(2)
cccc      call map__bc(3)
cc
cc      !X component BC
cc      if (isBdry(1,igx,1)) then
cc        i = 0
cc        xx0(i,1:ny,1:nz,1) = 0.5*(xold(i  ,1:ny,1:nz,1)
cc     .                           +xold(i+1,1:ny,1:nz,1))
cc     .         +0.25*(xx0(i,1:ny,0:nz-1,1)+xx0(i,0:ny-1,0:nz-1,1)
cc     .               +xx0(i,1:ny,1:nz  ,1)+xx0(i,0:ny-1,1:nz  ,1))
cc      endif
cc
cc      if (isBdry(nx,igx,2)) then
cc        i = nx+1
cc        xx0(i,1:ny,1:nz,1) = 0.5*(xold(i  ,1:ny,1:nz,1)
cc     .                           +xold(i-1,1:ny,1:nz,1))
cc     .         +0.25*(xx0(i-1,1:ny,0:nz-1,1)+xx0(i-1,0:ny-1,0:nz-1,1)
cc     .               +xx0(i-1,1:ny,1:nz  ,1)+xx0(i-1,0:ny-1,1:nz  ,1))
cc      endif
cc
cc      call map_bc(1,igx,nx,ny,nz,xnew(:,:,:,1),xx0(:,:,:,1),iorder=1)
cc
cc      !Y component BC
cc      if (isBdry(1,igy,3)) then
cc        j = 0
cc        xx0(1:nx,j,1:nz,2) = 0.5*(xold(1:nx,j  ,1:nz,2)
cc     .                           +xold(1:nx,j+1,1:nz,2))
cc     .         +0.25*(xx0(1:nx,j,0:nz-1,2)+xx0(0:nx-1,j,0:nz-1,2)
cc     .               +xx0(1:nx,j,1:nz  ,2)+xx0(0:nx-1,j,1:nz  ,2))
cc      endif
cc
cc      if (isBdry(ny,igy,4)) then
cc        j = ny+1
cc        xx0(1:nx,j,1:nz,2) = 0.5*(xold(1:nx,j  ,1:nz,2)
cc     .                           +xold(1:nx,j-1,1:nz,2))
cc     .         +0.25*(xx0(1:nx,j-1,0:nz-1,2)+xx0(0:nx-1,j-1,0:nz-1,2)
cc     .               +xx0(1:nx,j-1,1:nz  ,2)+xx0(0:nx-1,j-1,1:nz  ,2))
cc      endif
cc
cc      call map_bc(2,igx,nx,ny,nz,xnew(:,:,:,2),xx0(:,:,:,2),iorder=1)
cc
cc      !Z component BC
cc      if (isBdry(1,igz,5)) then
cc        k = 0
cc        xx0(1:nx,1:ny,k,3) = 0.5*(xold(1:nx,1:ny,k  ,3)
cc     .                           +xold(1:nx,1:ny,k+1,3))
cc     .         +0.25*(xx0(1:nx,0:ny-1,k,3)+xx0(0:nx-1,0:ny-1,k,3)
cc     .               +xx0(1:nx,1:ny  ,k,3)+xx0(0:nx-1,1:ny  ,k,3))
cc      endif
cc
cc      if (isBdry(nz,igz,6)) then
cc        k = nz+1
cc        xx0(1:nx,1:ny,k,3) = 0.5*(xold(1:nx,1:ny,k  ,3)
cc     .                           +xold(1:nx,1:ny,k-1,3))
cc     .         +0.25*(xx0(1:nx,0:ny-1,k-1,3)+xx0(0:nx-1,0:ny-1,k-1,3)
cc     .               +xx0(1:nx,1:ny  ,k-1,3)+xx0(0:nx-1,1:ny  ,k-1,3))
cc      endif
cc
cc      call map_bc(3,igx,nx,ny,nz,xnew(:,:,:,3),xx0(:,:,:,3),iorder=1)
cc
cccc        write (*,*)sum((xnew(1:nx,1:ny,1:nz,:)-xold(1:nx,1:ny,1:nz,:))**2)
cccc        write (*,*)sum((xnew(:,1:ny,1:nz,:)-xold(:,1:ny,1:nz,:))**2)
cccc        write (*,*)sum((xnew(1:nx,:,1:nz,:)-xold(1:nx,:,1:nz,:))**2)
cccc        write (*,*)sum((xnew(1:nx,1:ny,:,:)-xold(1:nx,1:ny,:,:))**2)
cccc        stop
cc
cc      contains
cc
ccc     map__bc
ccc     ###############################################################
cc      subroutine map__bc(ieq)
cc
cc      integer :: ieq,map_bcnd(6)
cc      real(8) :: x00 (0:nx+1,0:ny+1,0:nz+1,3)
cc
cc      map_bcnd = EQU
cc
cc      !X component BC
cc      if (isBdry(1,igx,1)) then
cc        i = 0
cc        x00(i,1:ny,1:nz,ieq) = 0.5*(xold(i  ,1:ny,1:nz,ieq)
cc     .                             +xold(i+1,1:ny,1:nz,ieq))
cc     .       +0.25*(xx0(i,1:ny,0:nz-1,ieq)+xx0(i,0:ny-1,0:nz-1,ieq)
cc     .             +xx0(i,1:ny,1:nz  ,ieq)+xx0(i,0:ny-1,1:nz  ,ieq))
cc      endif
cc
cc      if (isBdry(nx,igx,2)) then
cc        i = nx+1
cc        x00(i,1:ny,1:nz,ieq) = 0.5*(xold(i  ,1:ny,1:nz,ieq)
cc     .                             +xold(i-1,1:ny,1:nz,ieq))
cc     .       +0.25*(xx0(i-1,1:ny,0:nz-1,ieq)+xx0(i-1,0:ny-1,0:nz-1,ieq)
cc     .             +xx0(i-1,1:ny,1:nz  ,ieq)+xx0(i-1,0:ny-1,1:nz  ,ieq))
cc      endif
cc
cc      if (isBdry(1,igy,3)) then
cc        j = 0
cc        x00(1:nx,j,1:nz,ieq) = 0.5*(xold(1:nx,j  ,1:nz,ieq)
cc     .                             +xold(1:nx,j+1,1:nz,ieq))
cc     .       +0.25*(xx0(1:nx,j,0:nz-1,ieq)+xx0(0:nx-1,j,0:nz-1,ieq)
cc     .             +xx0(1:nx,j,1:nz  ,ieq)+xx0(0:nx-1,j,1:nz  ,ieq))
cc      endif
cc
cc      if (isBdry(ny,igy,4)) then
cc        j = ny+1
cc        x00(1:nx,j,1:nz,ieq) = 0.5*(xold(1:nx,j  ,1:nz,ieq)
cc     .                             +xold(1:nx,j-1,1:nz,ieq))
cc     .       +0.25*(xx0(1:nx,j-1,0:nz-1,ieq)+xx0(0:nx-1,j-1,0:nz-1,ieq)
cc     .             +xx0(1:nx,j-1,1:nz  ,ieq)+xx0(0:nx-1,j-1,1:nz  ,ieq))
cc      endif
cc
cc      if (isBdry(1,igz,5)) then
cc        k = 0
cc        x00(1:nx,1:ny,k,ieq) = 0.5*(xold(1:nx,1:ny,k  ,ieq)
cc     .                             +xold(1:nx,1:ny,k+1,ieq))
cc     .       +0.25*(xx0(1:nx,0:ny-1,k,ieq)+xx0(0:nx-1,0:ny-1,k,ieq)
cc     .             +xx0(1:nx,1:ny  ,k,ieq)+xx0(0:nx-1,1:ny  ,k,ieq))
cc      endif
cc
cc      if (isBdry(nz,igz,6)) then
cc        k = nz+1
cc        x00(1:nx,1:ny,k,ieq) = 0.5*(xold(1:nx,1:ny,k  ,ieq)
cc     .                             +xold(1:nx,1:ny,k-1,ieq))
cc     .       +0.25*(xx0(1:nx,0:ny-1,k-1,ieq)+xx0(0:nx-1,0:ny-1,k-1,ieq)
cc     .             +xx0(1:nx,1:ny  ,k-1,ieq)+xx0(0:nx-1,1:ny  ,k-1,ieq))
cc      endif
cc
cc      call map_bc(ieq,igx,nx,ny,nz,xnew(:,:,:,ieq),x00(:,:,:,ieq)
cc     .           ,iorder=1)!,bconds=map_bcnd)
cc
cc      end subroutine map__bc
cc
cc      end subroutine mk_xfer_mesh

c     mk_xfer_mesh
c     ##############################################################
      subroutine mk_xfer_mesh(igrid,nx,ny,nz,phi,xold,xnew)
c     --------------------------------------------------------------
c     Transfers x'=x+grad(phi) -> x at cell centers
c     --------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,nx,ny,nz
      real(8),INTENT(IN) :: phi (0:nx+1,0:ny+1,0:nz+1)
     .                     ,xold(0:nx+1,0:ny+1,0:nz+1,3)
      real(8) :: xnew(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer :: i,j,k,igx,igy,igz!,map_bcs(6,3)
      real(8) :: xx0 (0:nx+1,0:ny+1,0:nz+1,3)

      real(8) :: sendbuf(1),recbuf(1)

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

c     Find grad(phi) in Cartesian geom. at vertices

      xx0 = 0d0

      do k=0,nz
        do j=0,ny
          do i=0,nx
            xx0(i,j,k,:)=findPos(i,j,k,igx,igy,igz,nx,ny,nz,xold
     .                           ,phi,1,abs=.false.)
          enddo
        enddo
      enddo

c     Average to cell centers

      xnew(1:nx,1:ny,1:nz,:) = 
     .   0.125*(xx0(0:nx-1,1:ny,0:nz-1,:)+xx0(0:nx-1,0:ny-1,0:nz-1,:)
     .         +xx0(1:nx  ,1:ny,0:nz-1,:)+xx0(1:nx  ,0:ny-1,0:nz-1,:)
     .         +xx0(0:nx-1,1:ny,1:nz  ,:)+xx0(0:nx-1,0:ny-1,1:nz  ,:)
     .         +xx0(1:nx  ,1:ny,1:nz  ,:)+xx0(1:nx  ,0:ny-1,1:nz  ,:))

c     Fill ghost cells of xnew=grad(phi)

cc      map_bcs(:,1) = mk_bcs
cc      map_bcs(:,2) = mk_bcs
cc      map_bcs(:,3) = mk_bcs
cc
cc      where (map_bcs /= PER) map_bcs = EXT !Extrapolate map by default
cc
cc      where (mk_bcs(1:2) /= PER) map_bcs(1:2,1) = EQU
cc      where (mk_bcs(3:4) /= PER) map_bcs(3:4,2) = EQU
cc      where (mk_bcs(5:6) /= PER) map_bcs(5:6,3) = EQU

      do i=1,3
        call map_bc(i,igrid,nx,ny,nz,xnew(:,:,:,i),xx0(:,:,:,i),iorder=1
     .             ,bcs=mk_map_bcs(:,i))
      enddo

c     Add old map at cell centers

      xnew = xnew + xold

      end subroutine mk_xfer_mesh

c     MK_setup_grid
c     ###############################################################
      subroutine MK_setup_grid(igrid,phi)

c     ---------------------------------------------------------------
c     Sets up new MK quantities for nonlinear residual computation
c     and finds new grid metrics based on phi.
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid
        real(8) :: phi(:,:,:)

c     Local variables

        integer :: nx,ny,nz,igx,igy,igz
        real(8) :: dummy(size(phi,1),size(phi,2),size(phi,3))

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        nx = mk_gr_params%nxv(igx)
        ny = mk_gr_params%nyv(igy)
        nz = mk_gr_params%nzv(igz)

        !Impose BCs for MK potential
        call mk_applyBC(igrid,phi)

        !Setup MK residual evaluation
        call mk_setup_NLF(igrid,phi,x0)

        !Find new mesh at cell centers
        call mk_xfer_mesh(igrid,nx,ny,nz,phi,x0,xnew)

        !Calculate grid velocity (cnv)
        if (old_time_eval.or.source_eval.or.itime==0) then
          xnm  = xold
          xold = xnew

          jacnm= jacn
          jacn = gmetric%grid(igrid)%jac  !Old jacobian, about to be rewritten

          gvel = 0d0

cc          phi = 0d0

          !Spline map (for lagged monitor function only)
          if (lagged_err_mon) call splineX(xold,xcoef_ext=xold_coef)
        else

          !In Cartesian geometry
          if (bdf2) then
            gvel =-(cnp*xnew+cn*xold+cnm*xnm)/dt
          else
            gvel =-(xnew-xold)/dt
          endif
          gvel = XformToCurv(igrid,gvel,.true.)

        endif

        !Define metrics in current grid
        call defineGridMetric(grid_params,xcar=xnew,igr=igrid)
        gmetric => grid_params%gmetric !For backward compatibility

c     End program

      end subroutine MK_setup_grid

c     mk_applyBC
c     ##############################################################
      subroutine mk_applyBC(igrid,phi)

c     --------------------------------------------------------------
c     Apply boundary conditions to MK potential.
c     --------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid
      real(8) :: phi(:,:,:)

c     Local variables

      integer :: nx,ny,nz
      real(8) :: dummy(size(phi,1),size(phi,2),size(phi,3))

c     Begin program

      nx = mk_gr_params%nxv(igrid)
      ny = mk_gr_params%nyv(igrid)
      nz = mk_gr_params%nzv(igrid)

      dummy = 0d0
      call setBC(1,nx,ny,nz,phi,dummy,mk_bcs
     .          ,igrid,igrid,igrid,iorder=0)

c     End program

      end subroutine mk_applyBC

c     MK_spline_src
c     ##############################################################
      subroutine MK_spline_src(igrid,fsrc)

c     --------------------------------------------------------------
c     Spline source for source evaluation with moving mesh.
c     --------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid
      real(8) :: fsrc(:)

c     Local variables

      integer :: i,j,k,ig,jg,kg,igx,igy,igz,nx,ny,nz,neq,ieq

      type (var_array),pointer :: vsrc => null()

c     Spline-related

      real(8), dimension(:,:,:),pointer :: var,tmp_coef

      integer :: flg,alloc_stat

ccc     Debug
cc
cccc      integer :: ieq
cc      real(8),allocatable,dimension(:,:,:) :: debug2
cc      character(13) :: debugfile

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

c     Unpack vector fsrc

      call mapVectorToStructure(vsrc,fsrc)

c     Fill ghost cells (extrapolate at physical boundaries)

      neq = vsrc%nvar-1

      do ieq=1,neq
        vsrc%array_var(ieq)%bconds = bcond  !Topological boundary conditions
        where (bcond == DEF)  vsrc%array_var(ieq)%bconds = EXT

        call setBC(ieq,nx,ny,nz
     .            ,vsrc%array_var(ieq)%array
     .            ,u_0 %array_var(ieq)%array
     .            ,vsrc%array_var(ieq)%bconds
     .            ,igx,igy,igz,iorder=sp_order)
      enddo

c     Allocate spline coefficient arrays

      allocate(src_coef(nx+2,ny+2,nz+2,neq)
     .        ,tmp_coef(nx+2,ny+2,nz+2),STAT=alloc_stat)

c     Spline variables and current map for later source interpolation

      do ieq=1,neq
        var => vsrc%array_var(ieq)%array
        call splineFld(var,fcoef_ext=tmp_coef)
        src_coef(:,:,:,ieq) = tmp_coef
      enddo

cc      call splineX(xnew,xcoef_ext=x0_coef)

c diag ****
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      k = 1
cc
cc      nx = grid_params%nxv(igx)
cc      ny = grid_params%nyv(igx)
cc      nz = grid_params%nzv(igx)
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      do ieq=1,neqd
cc        debug2 = vsrc%array_var(ieq)%array
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,0d0,1d0,0d0,1d0,ieq-1,110)
cc      enddo
cc
cc      deallocate(debug2)
cc      close(110)
cc
cc#if defined(petsc)
cc      call PetscEnd(mpierr)
cc#endif
cc      stop
c diag ****

c     Free memory

      call deallocateDerivedType(vsrc)
      deallocate(tmp_coef)
      nullify(var)

c     End program

      end subroutine MK_spline_src

c     MK_eval_src
c     ##############################################################
      function MK_eval_src(igrid,i,j,k) result(src)

        implicit none

c     Call variables

        integer :: igrid,i,j,k
        real(8) :: src(neqd)

c     Local variables

        integer :: ilevel,ierror,ieq,ig,jg,kg
        real(8) :: x1,x2,x3,xi1,xi2,xi3

        real(8), dimension(:,:,:),pointer :: tmp_coef

c     Begin program

c     Find current physical coordinates

        x1 = gmetric%grid(igrid)%car(i,j,k,1)
        x2 = gmetric%grid(igrid)%car(i,j,k,2)
        x3 = gmetric%grid(igrid)%car(i,j,k,3)

c     Invert map to find logical coordinates

        !Initial guess for logical coords
        call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid,ig,jg,kg
     .                                ,xi1,xi2,xi3)

        ilevel = 0 !For debugging output level in evalXi
        call evalXi(ilevel,x1,x2,x3,xi1,xi2,xi3,ierror
     .             ,xcoef_ext=x0_coef)

        if (ierror /= 0) then
          call pstop('MK_eval_src','Map inversion failed')
        endif

c     Interpolate at logical coordinates

        allocate(tmp_coef(size(src_coef,1)
     .                   ,size(src_coef,2)
     .                   ,size(src_coef,3)))

        do ieq=1,neqd-1
          tmp_coef = src_coef(:,:,:,ieq)
          src(ieq) = evalFld(xi1,xi2,xi3,ierror,fcoef_ext=tmp_coef)
        enddo

        src(neqd) = 0d0

        deallocate(tmp_coef)

c     End program

      end function MK_eval_src

c     MK_spline_equ
c     ##############################################################
      subroutine MK_spline_equ(igrid,vequ)

c     --------------------------------------------------------------
c     Spline source for source evaluation with moving mesh.
c     --------------------------------------------------------------

cc      use xdraw_io

      implicit none

c     Call variables

      integer :: igrid

      type (var_array),pointer :: vequ

c     Local variables

      integer :: i,j,k,ig,jg,kg,igx,igy,igz,nx,ny,nz,neq,ieq

c     Spline-related

      real(8), dimension(:,:,:),pointer :: var,tmp_coef

      integer :: flg,alloc_stat

c     Debug

cc      integer :: ierror
cc      real(8) :: xi1,xi2,xi3
cc      real(8),allocatable,dimension(:,:,:) :: debug2
cc      character(13) :: debugfile

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

c     Fill ghost cells (extrapolate at physical boundaries)

      neq = vequ%nvar-1

      do ieq=1,neq
        vequ%array_var(ieq)%bconds = u_0%array_var(ieq)%bconds

        call setBC(ieq,nx,ny,nz
     .            ,vequ%array_var(ieq)%array
     .            ,u_0 %array_var(ieq)%array
     .            ,vequ%array_var(ieq)%bconds
     .            ,igx,igy,igz,iorder=sp_order)
      enddo

c diag ****
cc      write (*,*) 'DIAG -- MK_spline_equ'
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      k = 1
cc
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      do ieq=1,neq
cc        debug2 = vequ%array_var(ieq)%array
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,0d0,1d0,0d0,1d0,ieq-1,110)
cc      enddo
cc
cc      deallocate(debug2)
cc      close(110)
cc
cc#if defined(petsc)
cc      call PetscEnd(mpierr)
cc#endif
cc      stop
c diag ****

c     Allocate spline coefficient arrays

      allocate(equ_coef (nx+2,ny+2,nz+2,neq)
     .        ,tmp_coef (nx+2,ny+2,nz+2),STAT=alloc_stat)

c     Spline gradient function

      do ieq=1,neq
        var => vequ%array_var(ieq)%array
        call splineFld(var,fcoef_ext=tmp_coef)
        equ_coef(:,:,:,ieq) = tmp_coef
      enddo

c diag ****
cc      write (*,*) 'DIAG -- MK_spline_equ'
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      do ieq=1,neq
cccc        debug2 = vequ%array_var(ieq)%array
cc        do k=0,nz+1
cc          do j=0,ny+1
cc            do i=0,nx+1
cc              call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
cc     .                                      ,ig,jg,kg,xi1,xi2,xi3)
cc              tmp_coef = equ_coef(:,:,:,ieq)
cc              debug2(i,j,k) = evalFld(xi1,xi2,xi3,ierror
cc     .                               ,fcoef_ext=tmp_coef)
cc            enddo
cc          enddo
cc        enddo
cc
cc        k = 2
cc
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,0d0,1d0,0d0,1d0,ieq-1,110)
cc      enddo
cc
cc      deallocate(debug2)
cc      close(110)
cc
cc#if defined(petsc)
cc      call PetscEnd(mpierr)
cc#endif
cc      stop
c diag ****

c     Free memory

      deallocate(tmp_coef)
      nullify(var)

c     End program

      end subroutine MK_spline_equ

c     MK_eval_equ
c     ##############################################################
      function MK_eval_equ(igrid,i,j,k) result(equ)

        implicit none

c     Call variables

        integer :: igrid,i,j,k
        real(8) :: equ(neqd)

c     Local variables

        integer :: ilevel,ierror,ieq,ig,jg,kg
        real(8) :: x1,x2,x3,xi1,xi2,xi3

        real(8), dimension(:,:,:),pointer :: tmp_coef

c     Begin program

c     Get current physical coordinates

        x1 = gmetric%grid(igrid)%car(i,j,k,1)
        x2 = gmetric%grid(igrid)%car(i,j,k,2)
        x3 = gmetric%grid(igrid)%car(i,j,k,3)

c     Invert map to find logical coordinates based on map x0

        !Initial guess for logical coords
        call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid,ig,jg,kg
     .                                ,xi1,xi2,xi3)

        ilevel = 0 !For debugging output level in evalXi
        call evalXi(ilevel,x1,x2,x3,xi1,xi2,xi3,ierror
     .             ,xcoef_ext=x0_coef)

        if (ierror /= 0) then
          call pstop('MK_eval_equ','Map inversion failed')
        endif

c     Interpolate at logical coordinates

        allocate(tmp_coef(size(equ_coef,1)
     .                   ,size(equ_coef,2)
     .                   ,size(equ_coef,3)))

        do ieq=1,size(equ_coef,4)
          tmp_coef = equ_coef(:,:,:,ieq)
          equ(ieq) = evalFld(xi1,xi2,xi3,ierror,fcoef_ext=tmp_coef)
        enddo

        equ(neqd) = 0d0

        deallocate(tmp_coef)

c     End program

      end function MK_eval_equ

c     MK_xfer_equ
c     ##############################################################
      subroutine MK_xfer_equ(igrid,vequ)

        implicit none

c     Call variables

        integer :: igrid

        type (var_array),pointer :: vequ

c     Local variables

        integer :: i,j,k,igx,igy,igz,nx,ny,nz,ieq

        real(8) :: equ_val(neqd)

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        nx = mk_gr_params%nxv(igx)
        ny = mk_gr_params%nyv(igy)
        nz = mk_gr_params%nzv(igz)

        do k=1,nz
          do j=1,ny
            do i=1,nx
              equ_val = MK_eval_equ(igrid,i,j,k)
              do ieq=1,neqd-1
                vequ%array_var(ieq)%array(i,j,k) = equ_val(ieq)
              enddo
            enddo
          enddo
        enddo

        do ieq=1,neqd-1
          call setBC(ieq,nx,ny,nz
     .            ,vequ%array_var(ieq)%array
     .            ,u_0 %array_var(ieq)%array
     .            ,vequ%array_var(ieq)%bconds
     .            ,igx,igy,igz)
        enddo

c     End program

      end subroutine MK_xfer_equ

c     MK_output
c     ###################################################################
      subroutine MK_output

c     -------------------------------------------------------------------
c     MK output
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        real(8) :: gvmax(3),lgvmax(3)

c     Begin program

        if (ilevel > 0.and.my_rank == 0) write (*,*)

        if (itime.lt.inewtime) then
          if (relax_init_grid.and.my_rank == 0)
     .           write (*,*) '     >>> MK --- Relaxing initial grid...'
        else
          gvmax(1) = maxval(abs(gvel(:,:,:,1)))
          gvmax(2) = maxval(abs(gvel(:,:,:,2)))
          gvmax(3) = maxval(abs(gvel(:,:,:,3)))

#if defined(petsc)
          lgvmax = gvmax
          call MPI_Reduce(lgvmax,gvmax,3,MPI_DOUBLE_PRECISION
     .                   ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#endif

          if (my_rank == 0) 
     .      write (*,10) '      >>> MK --- Max grid velocities:',gvmax
        endif

c     End program

 10    format (a,1p3e12.3)

      end subroutine MK_output

      end module mk

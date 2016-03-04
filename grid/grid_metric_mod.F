
c module grid_metric_init
c #####################################################################
      module grid_metric_init

c ---------------------------------------------------------------------
c     This module allocates and fills the grid metric structure gmetric.
c     It contains the following routines: 
c        * allocateGridMetric
c        * deallocateGridMetric
c        * findGridMetrics
c           -> gammaZeroForce
c        * coarsenGridMetrics
c           -> coarsen
c        * topol_bc
c        * map_bc
c        * defineGridMetric
c ---------------------------------------------------------------------

        use grid_anal_map

        use grid_mpi

        use setMGBC_interface

        implicit none

        logical :: numerical_grid,check_grid

c     PRIVATE VARIABLES & PROCEDURES

        integer,private :: ilo ,ihi ,jlo ,jhi ,klo ,khi
     .                    ,ilom,ihip,jlom,jhip,klom,khip

        integer,private ::
     .             i,j,k,igx,igy,igz
     .            ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .            ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp

        logical,private :: map_present,vmap_present

        INTERFACE find_dr
          module procedure find_dr_1st, find_dr_2nd
        END INTERFACE

        private :: topol_bc,findGridmetrics,find_dr

      contains

c     alt__eom
c     #################################################################
      function alt__eom()

        logical :: alt__eom

        alt__eom = bcSP()

      end function alt__eom


c     checkMapDatabase
c     #################################################################
      function checkMapDatabase() result(is_map)

c     -----------------------------------------------------------------
c     Checks grid database for analytical mappings
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        logical :: is_map

c     Begin program

        is_map = checkAnalMapDatabase().or.(coords=='ext')

      end function checkMapDatabase

c     find_dr_1st
c     #####################################################################
      function find_dr_1st(i,j,k,igrid,hex,hey,hez,g_def,x) result(dr)

c     ---------------------------------------------------------------------
c     This routine finds dx/d_xi tensor at the point (i,j,k)
c     from the position vector x. It can do this for:
c        * Cell centers: (hex,hey,hez)=(0,0,0)
c        * Faces, e.g. (hex,hey,hez)=(1,0,0) for i+1/2 face.
c        * Corners, e.g. (hex,hey,hez)=(1,1,1) for i+1/2,j+1/2,k+1/2 corner.
c     on an arbitrary grid level.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,igrid,hex,hey,hez
        real(8)    :: dr(3,3),x(0:,0:,0:,:)

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
     .               ,igp,igm,jgp,jgm,kgp,kgm,hel,nx,ny,nz

        real(8)    :: carp(3),carm(3),dh(3),coeff

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        ip=min(i+1,nx+1)
        im=max(i-1,0)
        jp=min(j+1,ny+1)
        jm=max(j-1,0)
        kp=min(k+1,nz+1)
        km=max(k-1,0)

        call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg,g_def=g_def)

c     Exceptions

        if (hex == 1) im = i
        if (hey == 1) jm = j
        if (hez == 1) km = k

        !X
cc        if (nx == 1) then
cc          dr(1,:) = (/ 1d0, 0d0, 0d0 /)
cc        else
          igp = ig+(ip-i)
          igm = ig-(i-im)

          dh(1)= (g_def%xx(igp)-g_def%xx(igm))

          carp = x(ip,j,k,:)
          carm = x(im,j,k,:)
          coeff = 1d0

          if (hey == 1) then
            carp = carp + x(ip,jp,k,:)
            carm = carm + x(im,jp,k,:)
            coeff = 0.5
          endif

          if (hez == 1) then
            carp = carp + x(ip,j,kp,:)
            carm = carm + x(im,j,kp,:)
            coeff = 0.5
          endif

          if (hey == 1 .and. hez == 1) then
            carp = carp + x(ip,jp,kp,:)
            carm = carm + x(im,jp,kp,:)
            coeff = 0.25
          endif

          dr(:,1) = coeff*(carp-carm)/dh(1)
cc        endif

        !Y
cc        if (ny == 1) then
cc          dr(2,:) = (/ 0d0, 1d0, 0d0 /)
cc        else
          jgp = jg+(jp-j)
          jgm = jg-(j-jm)

          dh(2)= (g_def%yy(jgp)-g_def%yy(jgm))

          carp = x(i,jp,k,:)
          carm = x(i,jm,k,:)
          coeff = 1d0

          if (hex == 1) then
            carp = carp + x(ip,jp,k,:)
            carm = carm + x(ip,jm,k,:)
            coeff = 0.5
          endif

          if (hez == 1) then
            carp = carp + x(i,jp,kp,:)
            carm = carm + x(i,jm,kp,:)
            coeff = 0.5
          endif

          if (hex == 1 .and. hez == 1) then
            carp = carp + x(ip,jp,kp,:)
            carm = carm + x(ip,jm,kp,:)
            coeff = 0.25
          endif

          dr(:,2) = coeff*(carp-carm)/dh(2)
cc        endif

        !Z
cc        if (nz == 1) then
cc          dr(3,:) = (/ 0d0, 0d0, 1d0 /)
cc        else
          kgp = kg+(kp-k)
          kgm = kg-(k-km)

          dh(3)= (g_def%zz(kgp)-g_def%zz(kgm))

          carp = x(i,j,kp,:)
          carm = x(i,j,km,:)
          coeff = 1d0

          if (hex == 1) then
            carp = carp + x(ip,j,kp,:)
            carm = carm + x(ip,j,km,:)
            coeff = 0.5
          endif

          if (hey == 1) then
            carp = carp + x(i,jp,kp,:)
            carm = carm + x(i,jp,km,:)
            coeff = 0.5
          endif

          if (hex == 1 .and. hey == 1) then
            carp = carp + x(ip,jp,kp,:)
            carm = carm + x(ip,jp,km,:)
            coeff = 0.25
          endif

          dr(:,3) = coeff*(carp-carm)/dh(3)
cc        endif

      end function find_dr_1st

c     find_dr_2nd
c     #####################################################################
      function find_dr_2nd(i,j,k,igrid,g_def,x) result(dr)

c     ---------------------------------------------------------------------
c     This routine finds dx/d_xi tensor at the point (i,j,k)
c     from the position vector x. It can do this for:
c        * Cell centers: (hex,hey,hez)=(0,0,0)
c        * Faces, e.g. (hex,hey,hez)=(1,0,0) for i+1/2 face.
c        * Corners, e.g. (hex,hey,hez)=(1,1,1) for i+1/2,j+1/2,k+1/2 corner.
c     on an arbitrary grid level.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,j,k,igrid
        real(8) :: dr(3,3),x(0:,0:,0:,:)

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: ip,im,i0,jp,jm,j0,kp,km,k0,ig,jg,kg
     .            ,igp,igm,jgp,jgm,kgp,kgm,nx,ny,nz

        real(8) :: carp(3),carm(3),car0(3),coeff

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

c     Exceptions

        !X
        i0 = min(max(i,1),nx)
        ip = i0+1
        im = i0-1

        call getMGmap(i0,j,k,igx,igy,igz,ig,jg,kg,g_def=g_def)

        igp = ig+1
        igm = ig-1

        carp = x(ip,j,k,:)
        car0 = x(i0,j,k,:)
        carm = x(im,j,k,:)

        dr(:,1) = dr_dxi(g_def%xx(igm),carm
     .                  ,g_def%xx(ig ),car0
     .                  ,g_def%xx(igp),carp
     .                  ,i-i0)

        !Y
        j0 = min(max(j,1),ny)
        jp = j0+1
        jm = j0-1

        call getMGmap(i,j0,k,igx,igy,igz,ig,jg,kg,g_def=g_def)

        jgp = jg+1
        jgm = jg-1

        carp = x(i,jp,k,:)
        car0 = x(i,j0,k,:)
        carm = x(i,jm,k,:)

        dr(:,2) = dr_dxi(g_def%yy(jgm),carm
     .                  ,g_def%yy(jg ),car0
     .                  ,g_def%yy(jgp),carp
     .                  ,j-j0)

        !Z
        k0 = min(max(k,1),nz)
        kp = k0+1
        km = k0-1

        call getMGmap(i,j,k0,igx,igy,igz,ig,jg,kg,g_def=g_def)

        kgp = kg+1
        kgm = kg-1

        carp = x(i,j,kp,:)
        car0 = x(i,j,k0,:)
        carm = x(i,j,km,:)

        dr(:,3) = dr_dxi(g_def%zz(kgm),carm
     .                  ,g_def%zz(kg ),car0
     .                  ,g_def%zz(kgp),carp
     .                  ,k-k0)

      contains

c     dr_dxi
c     #################################################################
      function dr_dxi(xm,rm,x0,r0,xp,rp,pos) result(dr)

        implicit none

        integer :: pos
        real(8) :: xm,x0,xp,rm(3),r0(3),rp(3),dr(3)

        real(8) :: dx1,dx2

c     Begin program

        dx1=x0-xm
        dx2=xp-x0

        select case(pos)
        case(-1)  !Left boundary
          dr = (-rp* dx1/dx2/(dx1+dx2)
     .          +r0*(dx1+dx2)/dx1/dx2
     .          -rm*(1./dx1+1./(dx1+dx2)))
        case(1)  !Right boundary
          dr = -(-rm* dx1/dx2/(dx1+dx2)
     .           +r0*(dx1+dx2)/dx1/dx2
     .           -rp*(1./dx1+1./(dx1+dx2)))
        case(0)
          dr = ( rp* dx1/dx2/(dx1+dx2)
     .          +r0*(dx2-dx1)/dx2/dx1
     .          -rm* dx2/dx1/(dx1+dx2))
        case default
          write (*,*) 'Error in dr_dxi'
          stop
        end select

      end function dr_dxi

      end function find_dr_2nd

c     coarsen_map
c     #################################################################
      subroutine coarsen_map(igc,g_def,xf,xc,bcs)

c     -----------------------------------------------------------------
c     Coarsen fine map xf(xi) to coarse map xc(xi). This is done
c     dimension by dimension, using first-order extrapolation.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igc
        real(8) :: xf(0:,0:,0:,:)
     .            ,xc(0:,0:,0:,:)

        type(grid_mg_def),pointer :: g_def

        integer,optional :: bcs(6,size(xf,4))

c     Local variables

        integer :: ii,igf,nxf,nyf,nzf,nxc,nyc,nzc,neq

c     Begin program

        igf = igc - 1

        nxf = g_def%nxv(igf)
        nyf = g_def%nyv(igf)
        nzf = g_def%nzv(igf)

        nxc = g_def%nxv(igc)
        nyc = g_def%nyv(igc)
        nzc = g_def%nzv(igc)

        neq = size(xf,4)

        if (PRESENT(bcs)) then
          do ii=1,neq
            call crse_map_comp(ii,xf(:,:,:,ii)
     .                           ,xc(:,:,:,ii)
     .                           ,bcs=bcs(:,ii))
          enddo
        else
          do ii=1,neq
            call crse_map_comp(ii,xf(:,:,:,ii)
     .                           ,xc(:,:,:,ii))
          enddo
        endif

      contains

c     crse_map_comp
c     #################################################################
      subroutine crse_map_comp(dim,xf_j,xc_j,bcs)

c     -----------------------------------------------------------------
c     Coarsen component "dim" of xf_j to xc_j one grid level.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: dim
        real(8)    :: xf_j(0:nxf+1,0:nyf+1,0:nzf+1)
     .               ,xc_j(0:nxc+1,0:nyc+1,0:nzc+1)

        integer,optional :: bcs(6)

c     Local variables

        integer    :: igf,ic,jc,kc,if,jf,kf,mgx,mgy,mgz
     .               ,i_off,j_off,k_off
        real(8)    :: x0_j(0:nxc+1,0:nyc+1,0:nzc+1)

        real(8)    :: volc,ac

c     Begin program

        igf = igc - 1

        mgx = g_def%mg_ratio_x(igf)
        mgy = g_def%mg_ratio_y(igf)
        mgz = g_def%mg_ratio_z(igf)

        do kc = 0,nzc+1
          do jc = 0,nyc+1
            do ic = 0,nxc+1

              i_off = min(max(ic,1),nxc)-ic
              j_off = min(max(jc,1),nyc)-jc
              k_off = min(max(kc,1),nzc)-kc

              !Fill inner domain
              volc   = 0d0
              ac     = 0d0

              do kf=mgz*((kc+k_off)-1)+1-k_off,mgz*(kc+k_off)-k_off
                do jf=mgy*((jc+j_off)-1)+1-j_off,mgy*(jc+j_off)-j_off
                  do if=mgx*((ic+i_off)-1)+1-i_off,mgx*(ic+i_off)-i_off
                    ac   = ac   + xf_j(if,jf,kf)
                    volc = volc + 1d0
                  enddo
                enddo
              enddo

              xc_j(ic,jc,kc) = ac/volc

            enddo
          enddo
        enddo

c       Prepare boundaries

        x0_j = xc_j

        if ((nxc == 1).and.(dim == 1)) then
          i=0
          x0_j(i,:,:) = 0.5*(x0_j(i,:,:)+x0_j(i+1,:,:))
          i=nxc+1
          x0_j(i,:,:) = 0.5*(x0_j(i,:,:)+x0_j(i-1,:,:))
        endif

        if ((nyc == 1).and.(dim == 2)) then
          j=0
          x0_j(:,j,:) = 0.5*(x0_j(:,j,:)+x0_j(:,j+1,:))
          j=nyc+1
          x0_j(:,j,:) = 0.5*(x0_j(:,j,:)+x0_j(:,j-1,:))
        endif

        if ((nzc == 1).and.(dim == 3)) then
          k=0
          x0_j(:,:,k) = 0.5*(x0_j(:,:,k)+x0_j(:,:,k+1))
          k=nzc+1
          x0_j(:,:,k) = 0.5*(x0_j(:,:,k)+x0_j(:,:,k-1))
        endif

c     Fill ghost cells

        call map_bc(dim,igc,nxc,nyc,nzc,xc_j,x0_j,bcs=bcs
     .             ,g_def=g_def)

      end subroutine crse_map_comp

      end subroutine coarsen_map

c     map_bc
c     #################################################################
      subroutine map_bc(dim,igrid,nxl,nyl,nzl,x_j,x0_j,iorder,bcs,g_def)

c     -----------------------------------------------------------------
c     Impose BCs on map component x_j using boundary of x0_j.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer :: dim,igrid,nxl,nyl,nzl

        real(8) :: x_j (0:nxl+1,0:nyl+1,0:nzl+1)
     .            ,x0_j(0:nxl+1,0:nyl+1,0:nzl+1)

        integer,optional,intent(IN) :: iorder
        integer,optional :: bcs(6)

        type(grid_mg_def),optional,pointer :: g_def

c     Local variables

        integer :: bcnd(6),iibc,order

c     Begin program

        if (PRESENT(iorder)) then
          order = iorder
        else
          order = 1    !Default to first order for predictable results at boundaries
        endif

c     Setup boundary conditions

        if (PRESENT(bcs)) then
          bcnd = bcs
        else
          bcnd = bcond

          !Force boundaries along dimension
          select case(dim)
          case(1)
            bcnd(1:2) = EQU
          case(2)
            bcnd(3:4) = EQU
          case(3)
            bcnd(5:6) = EQU
          end select

          where(bcnd == DEF) bcnd = EXT
        endif

c     Impose BCs

        call setMGBC(0,1,nxl,nyl,nzl,igrid,x_j,bcnd
     .              ,arr0=x0_j,iorder=order,g_def=g_def)

c     End program

      end subroutine map_bc

c     findGridMetrics
c     #################################################################
      subroutine findGridMetrics(g_def,igrid,coarsen,ierr)

c     -----------------------------------------------------------------
c     Find grid metrics at grid level igrid.
c     -----------------------------------------------------------------

        use grid_debug

        implicit none

c     Input variables

        integer :: igrid
        type(grid_mg_def),pointer :: g_def

        logical :: coarsen

        integer,optional :: ierr

c     Local variables

        integer :: ii,nxf,nyf,nzf,nxc,nyc,nzc,nx,ny,nz,ig,jg,kg,m,l
     .            ,lerr

        logical :: return_error

c     Begin program

        return_error = PRESENT(ierr)

        if (return_error) lerr = 0

        igx = igrid
        igy = igrid
        igz = igrid

        ilo = 1
        jlo = 1
        klo = 1
        ihi = g_def%nxv(igrid)
        jhi = g_def%nyv(igrid)
        khi = g_def%nzv(igrid)

        ilom = ilo-1
        jlom = jlo-1
        klom = klo-1
        ihip = ihi+1
        jhip = jhi+1
        khip = khi+1

        nx = ihi-ilo+1
        ny = jhi-jlo+1
        nz = khi-klo+1

c       ANALYTICAL QUANTITIES **********************************

        if ((.not.numerical_grid).and.anal_map) then

          !Cell-centered metrics
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                g_def%gmetric%grid(igrid)%car  (i,j,k,:)
     .                    = map     (i,j,k,igx,igy,igz
     .                              ,g_def=g_def)
                g_def%gmetric%grid(igrid)%jac  (i,j,k)
     .                    = jacobian(i,j,k,igx,igy,igz
     .                              ,g_def=g_def)
                g_def%gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                    = g_sub   (i,j,k,igx,igy,igz
     .                              ,g_def=g_def)
                g_def%gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                    = g_sup   (i,j,k,igx,igy,igz
     .                              ,g_def=g_def)
                g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .                    = christ_2knd(i,j,k,igx,igy,igz
     .                                 ,g_def=g_def)
                do l=1,3
                  g_def%gmetric%grid(igrid)%cov(i,j,k,l,:)
     .                    = covariantVector    (l,i,j,k,igx,igy,igz
     .                                         ,g_def=g_def)
                  g_def%gmetric%grid(igrid)%cnv(i,j,k,l,:)
     .                    = contravariantVector(l,i,j,k,igx,igy,igz
     .                                         ,g_def=g_def)
                enddo

                !Cell volume
                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg,g_def=g_def)

                g_def%gmetric%grid(igrid)%dlvol(i,j,k) = 1d0

                if (g_def%nxgl(igrid) > 1)
     .               g_def%gmetric%grid(igrid)%dlvol(i,j,k) =
     .               g_def%gmetric%grid(igrid)%dlvol(i,j,k) *
     .                           g_def%dxh(ig)

                if (g_def%nygl(igrid) > 1)
     .               g_def%gmetric%grid(igrid)%dlvol(i,j,k) =
     .               g_def%gmetric%grid(igrid)%dlvol(i,j,k) *
     .                           g_def%dyh(jg)

                if (g_def%nzgl(igrid) > 1)
     .               g_def%gmetric%grid(igrid)%dlvol(i,j,k) =
     .               g_def%gmetric%grid(igrid)%dlvol(i,j,k) *
     .                           g_def%dzh(kg)

                g_def%gmetric%grid(igrid)%dvol(i,j,k)=
     .                          g_def%gmetric%grid(igrid)%jac(i,j,k)
     .                         *g_def%gmetric%grid(igrid)%dlvol(i,j,k)
              enddo
            enddo
          enddo

          !Vertex metrics
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                g_def%gmetric%grid(igrid)%jac_v(i,j,k)
     .                    = jacobian(i,j,k,igx,igy,igz,vertex=.true.
     .                              ,g_def=g_def)
                g_def%gmetric%grid(igrid)%gsub_v(i,j,k,:,:)
     .                    = g_sub   (i,j,k,igx,igy,igz,vertex=.true.
     .                              ,g_def=g_def)
                g_def%gmetric%grid(igrid)%gsup_v (i,j,k,:,:)
     .                    = g_sup   (i,j,k,igx,igy,igz,vertex=.true.
     .                              ,g_def=g_def)
                do l=1,3
                  g_def%gmetric%grid(igrid)%cov_v(i,j,k,l,:)
     .                    = covariantVector    (l,i,j,k,igx,igy,igz
     .                                         ,vertex=.true.
     .                                         ,g_def=g_def)
                  g_def%gmetric%grid(igrid)%cnv_v(i,j,k,l,:)
     .                    = contravariantVector(l,i,j,k,igx,igy,igz
     .                                         ,vertex=.true.
     .                                         ,g_def=g_def)
                enddo

              enddo
            enddo
          enddo

          !Zero force condition on Christoffel symbols (only on finest grid)
          if (igrid == 1) call hessianZeroForceCorr(igrid)

c       NUMERICAL QUANTITIES ***********************************

        else

          numerical_grid = .true.

          call findCellMap(g_def)

          call findVertexMap(g_def)

          call findVertexMetrics(g_def)

          call findCellMetrics(g_def)

        endif

        if (return_error) ierr = lerr

c     End program

      contains

c     findCellMap
c     #################################################################
      subroutine findCellMap(g_def)

c     -----------------------------------------------------------------
c     Finds Cartesian map at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: xf,xc

c     Begin program

        xc => g_def%gmetric%grid(igrid)%car

        if (map_present) then

          if (igrid > 1 .and. coarsen) then  !Coarsen map if at coarse grid level
            xf => g_def%gmetric%grid(igrid-1)%car

            call coarsen_map(igrid,g_def,xf,xc)

            nullify(xf)
          endif

        else

          if (checkAnalMapDatabase()) then

            do k = klom,khip
              do j = jlom,jhip
                do i = ilom,ihip
                  xc(i,j,k,:) = map(i,j,k,igx,igy,igz,g_def=g_def)
                enddo
              enddo
            enddo

          else
            call pstop('findCellMap','Map not available')
          endif

        endif

c     End program

      end subroutine findCellMap

c     findVertexMap
c     #################################################################
      subroutine findVertexMap(g_def)

c     -----------------------------------------------------------------
c     Finds Cartesian map at vertices.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: xc,xv

c     Begin program

        xv => g_def%gmetric%grid(igrid)%car_v

        if ((.not.vmap_present).or.(igrid > 1)) then

          xc => g_def%gmetric%grid(igrid)%car !Cell center map

          xv(0:nx,0:ny,0:nz,:) = 0.125*(xc(0:nx  ,0:ny  ,0:nz  ,:)
     .                                 +xc(0:nx  ,1:ny+1,0:nz  ,:)
     .                                 +xc(0:nx  ,0:ny  ,1:nz+1,:)
     .                                 +xc(0:nx  ,1:ny+1,1:nz+1,:)
     .                                 +xc(1:nx+1,0:ny  ,0:nz  ,:)
     .                                 +xc(1:nx+1,1:ny+1,0:nz  ,:)
     .                                 +xc(1:nx+1,0:ny  ,1:nz+1,:)
     .                                 +xc(1:nx+1,1:ny+1,1:nz+1,:))

        endif

c     End program

        nullify(xv,xc)

      end subroutine findVertexMap

c     findVertexMetrics
c     #################################################################
      subroutine findVertexMetrics(g_def)

c     -----------------------------------------------------------------
c     Finds metrics at vertices (i+1/2,j+1/2,k+1/2) (stored at (i,j,k))
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: he,imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz,bcs(6)
     .            ,bcond_sv(6)

        real(8) :: r(3,3),cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),jac,ijac
     .            ,gamma(3,3,3),d2r(3,3)

        real(8),pointer,dimension(:,:,:,:,:) :: dr
        real(8),pointer,dimension(:,:,:,:)   :: xc

c     Begin program

        imin = ilom
        imax = ihip

        jmin = jlom
        jmax = jhip

        kmin = klom
        kmax = khip

        he = 1

        allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))

c     Local grid size

        nx = ihi-ilo+1
        ny = jhi-jlo+1
        nz = khi-klo+1

c     Evaluate dr=dx/dxi vectors at vertices

        xc => g_def%gmetric%grid(igrid)%car

        dr = 0d0

        do k = klom,khi
          do j = jlom,jhi
            do i = ilom,ihi
              dr(i,j,k,:,:)
     .             = find_dr(i,j,k,igrid,he,he,he,g_def,xc)
            enddo
          enddo
        enddo

        !Modify BC setup to do staggered
        bcond_sv = bcond
        bcs = bcond
        !Do no overwrite i=j=k=0
        if (bcs(2) == DEF) bcs(2) = EXT
        if (bcs(4) == DEF) bcs(4) = EXT 
        if (bcs(6) == DEF) bcs(6) = EXT 
        if (bcs(2) == FSYM .or. bcs(2) == SYM) then
          bcs  (2) = NEU
          bcond(2) = SYM  !Set BCs to do NEU across face
        endif
        if (bcs(4) == FSYM .or. bcs(4) == SYM) then
          bcs  (4) = NEU 
          bcond(4) = SYM  !Set BCs to do NEU across face
        endif
        if (bcs(6) == FSYM .or. bcs(6) == SYM) then
          bcs  (6) = NEU
          bcond(6) = SYM  !Set BCs to do NEU across face
        endif

        call set_force_corners(.true.)

        do i=1,3
          do j=1,3
            call setMGBC(0,1,nx,ny,nz,igrid,dr(:,:,:,i,j),bcs
     .                  ,iorder=2,g_def=g_def)
          enddo
        enddo

        !Recover original setup
        call set_force_corners(.false.)
        bcond = bcond_sv

c     Evaluate grid quantities (at vertices)

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(:,1),r(:,2),r(:,3))
              g_def%gmetric%grid(igrid)%jac_v(i,j,k) = jac

              if (jac <= 0d0 .and. (.not.isSP2(i+1,igrid,ibc=1))) then
                if (.not.check_grid .and. (.not.return_error)) then
                  write (*,*) 'Negative Jacobian in half mesh'
                  write (*,*) 'Grid level',igrid
                  write (*,*) 'Processor rank=',my_rank
                  write (*,*) 'Grid position=',i,j,k
                  write (*,*) 'dr/dxi_1=',r(:,1)
                  write (*,*) 'dr/dxi_2=',r(:,2)
                  write (*,*) 'dr/dxi_3=',r(:,3)
                  write (*,*) 'Jacobian=',jac

                  call sstop('findVertexMetrics','Inverted Jacobian')
                elseif (return_error) then
                  lerr = 1 ; exit
                endif
              elseif (isSP2(i+1,igrid,ibc=1,g_def=g_def)) then
cc              if (isSP2(i+1,igrid,ibc=1)) then
                ijac = 0d0
              else
                ijac = 1d0/jac
              endif

              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
              cnv(:,:) = r(:,:)*ijac

              g_def%gmetric%grid(igrid)%cnv_v(i,j,k,:,:)
     .                 =transpose(cnv)

              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
              cov(:,1) = cross_product(r(:,2),r(:,3))
              cov(:,2) = cross_product(r(:,3),r(:,1))
              cov(:,3) = cross_product(r(:,1),r(:,2))

              cov = cov*ijac

              g_def%gmetric%grid(igrid)%cov_v(i,j,k,:,:)
     .               = transpose(cov)

              !Metric tensors <<<<<<<<<<<<<<<<<<<<<<<<
              do m=1,3
                do l=m,3
                  gsub(l,m) = jac*dot_product(cnv(:,l),cnv(:,m))
                  gsub(m,l) = gsub(l,m) !Symmetry
                  gsup(l,m) = jac*dot_product(cov(:,l),cov(:,m))
                  gsup(m,l) = gsup(l,m) !Symmetry
                enddo
              enddo
              g_def%gmetric%grid(igrid)%gsub_v(i,j,k,:,:) = gsub
              g_def%gmetric%grid(igrid)%gsup_v(i,j,k,:,:) = gsup

            enddo
          enddo
        enddo

c     Christoffel symbols (at cell centers; BCs enforced in findCellMetrics)

        do k = kmin+1,kmax
          do j = jmin+1,jmax
            do i = imin+1,imax

              cov = 0.125*(
     .             g_def%gmetric%grid(igrid)%cov_v(i-1,j-1,k-1,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i  ,j-1,k-1,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i-1,j  ,k-1,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i  ,j  ,k-1,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i-1,j-1,k  ,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i  ,j-1,k  ,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i-1,j  ,k  ,:,:)
     .            +g_def%gmetric%grid(igrid)%cov_v(i-1,j-1,k  ,:,:))

              do l=1,3
                do m=1,3
                  do n=1,3
                    d2r=find_dr(i-1,j-1,k-1,igrid,1,1,1,g_def
     .                         ,dr(:,:,:,:,m))
                    gamma(l,m,n) = dot_product(d2r(:,n),cov(l,:))
                  enddo
                enddo
                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
              enddo

              g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma

            enddo
          enddo
        enddo

c     End program

        deallocate(dr)
        nullify(xc)

      end subroutine findVertexMetrics

c     findCellMetrics
c     #################################################################
      subroutine findCellMetrics(g_def)  !OLD version

c     -----------------------------------------------------------------
c     Finds metrics at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: he,imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz

        real(8) :: r(3,3),car0(3),carp(3),carm(3),jac,ijac
     .            ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .            ,gamma(3,3,3),mag,dhh,dhm,dhp,d2r(3,3)

        real(8),pointer,dimension(:,:,:,:,:) :: dr
        real(8),pointer,dimension(:,:,:,:)   :: xc

c     Begin program

        imin = ilom
        imax = ihip

        jmin = jlom
        jmax = jhip

        kmin = klom
        kmax = khip

c     Local grid size

        nx = ihi-ilo+1
        ny = jhi-jlo+1
        nz = khi-klo+1

c     Evaluate dr=dx/dxi vectors

        xc => g_def%gmetric%grid(igrid)%car

        allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))

        dr = 0d0

        if (igrid == 1) then  !Second order for finest grid
          do k = kmin,kmax
            do j = jmin,jmax
              do i = imin,imax
                dr(i,j,k,:,:) = find_dr(i,j,k,igrid,g_def,xc)
              enddo
            enddo
          enddo
        else                  !First order in coarser grids (to avoid negative jacobians)
          do k = kmin,kmax
            do j = jmin,jmax
              do i = imin,imax
                dr(i,j,k,:,:) = find_dr(i,j,k,igrid,0,0,0,g_def,xc)
              enddo
            enddo
          enddo
        endif

c     Enforce topological constraints on dr (including SP)

        do i=1,3
          do j=1,3
            call topol_bc(g_def,igrid,nx,ny,nz,dr(:,:,:,i,j))
            if (bcSP()) then
              igxbc = igrid
              igybc = igrid
              igzbc = igrid

              nxbc = nx
              nybc = ny
              nzbc = nz
              !dr behaves like a covariant vector
              call vectorSingularBC(j,dr(:,:,:,i,j),.true.)
            endif
          enddo
        enddo

c     Evaluate grid quantities

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(:,1),r(:,2),r(:,3))
              g_def%gmetric%grid(igrid)%jac(i,j,k) = jac

              if (jac <= 0d0 .and.(.not.isSP2(i+1,igrid,ibc=1))) then
                if (.not.check_grid .and. (.not.return_error)) then
                  write (*,*) 'Negative Jacobian in integer mesh'
                  write (*,*) 'Grid level',igrid
                  write (*,*) 'Grid position=',i,j,k
                  write (*,*) 'dr/dxi_1=',r(:,1)
                  write (*,*) 'dr/dxi_2=',r(:,2)
                  write (*,*) 'dr/dxi_3=',r(:,3)
                  write (*,*) 'Jacobian=',jac

                  call sstop('findCellMetrics','Inverted jacobian')
                elseif (return_error) then
                  lerr = 2 ; exit
                endif
              endif

              ijac = 1d0/jac

              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
              cnv(:,:) = r(:,:)*ijac

              g_def%gmetric%grid(igrid)%cnv(i,j,k,:,:) =transpose(cnv)

              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
              cov(:,1) = cross_product(r(:,2),r(:,3))
              cov(:,2) = cross_product(r(:,3),r(:,1))
              cov(:,3) = cross_product(r(:,1),r(:,2))

              cov = cov*ijac

              g_def%gmetric%grid(igrid)%cov(i,j,k,:,:) =transpose(cov)

              !Metric tensors <<<<<<<<<<<<<<<<<<<<<<<<
              do m=1,3
                do l=m,3
                  gsub(l,m) = jac*dot_product(cnv(:,l),cnv(:,m))
                  gsub(m,l) = gsub(l,m) !Symmetry
                  gsup(l,m) = jac*dot_product(cov(:,l),cov(:,m))
                  gsup(m,l) = gsup(l,m) !Symmetry
                enddo
              enddo
              g_def%gmetric%grid(igrid)%gsub(i,j,k,:,:) = gsub
              g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:) = gsup

              !Cell volume <<<<<<<<<<<<<<<<<<
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg,g_def=g_def)

              g_def%gmetric%grid(igrid)%dlvol(i,j,k) = 1d0

              if (g_def%nxgl(igrid) > 1)
     .             g_def%gmetric%grid(igrid)%dlvol(i,j,k) =
     .             g_def%gmetric%grid(igrid)%dlvol(i,j,k) *
     .                         g_def%dxh(ig)

              if (g_def%nygl(igrid) > 1)
     .             g_def%gmetric%grid(igrid)%dlvol(i,j,k) =
     .             g_def%gmetric%grid(igrid)%dlvol(i,j,k) *
     .                         g_def%dyh(jg)

              if (g_def%nzgl(igrid) > 1)
     .             g_def%gmetric%grid(igrid)%dlvol(i,j,k) =
     .             g_def%gmetric%grid(igrid)%dlvol(i,j,k) *
     .                         g_def%dzh(kg)

              g_def%gmetric%grid(igrid)%dvol(i,j,k)=
     .                         g_def%gmetric%grid(igrid)%jac(i,j,k)
     .                        *g_def%gmetric%grid(igrid)%dlvol(i,j,k)

              !Christoffel symbols <<<<<<<<<<<<<<<<<<
              do l=1,3
                do m=1,3
                  do n=1,3
                    d2r = find_dr(i,j,k,igrid,g_def,dr(:,:,:,:,m))
                    gamma(l,m,n) = dot_product(d2r(:,n),cov(:,l))
                  enddo
                enddo
                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
              enddo

              g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma

            enddo
          enddo
        enddo

        !Zero-force condition on Christoffle symbols 
        if (igrid == 1) call hessianZeroForceCorr(igrid)

        !Enforce topological constraints on Christoffel symbols
        do i=1,3
          do j=1,3
            do k=1,3
              call topol_bc(g_def,igrid,nx,ny,nz
     .           ,g_def%gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
            enddo
          enddo
        enddo

        deallocate(dr)
        nullify(xc)

c     End program

      end subroutine findCellMetrics

c     hessianZeroForceCorr
c     #################################################################
      subroutine hessianZeroForceCorr(igrid)

c     -----------------------------------------------------------------
c     Postprocess Christoffel symbol to satisfy zero-force condition
c     to machine accuracy
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer    :: igrid

c     Local variables

        integer :: i,j,k,ii,ll,kk,mm,jj,ig,jg,kg
        real(8) :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8) :: hess(3,3,3),table(3,3,3),gamma(3,3,3)
     .            ,gsub(3,3),gsup(3,3),gsupp(3,3),gsupm(3,3)
        logical :: alt_eom

c     Begin program

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

            gamma = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
            gsub  = g_def%gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsup  = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg,g_def=g_def)

            do ii=1,3 !cycle through Chistoffel symbols

              alt_eom = (ii == 2 .and. alt__eom())

              if (alt_eom) then
                jac  = g_def%gmetric%grid(igrid)%jac(i,j,k)
                const= 2d0
              else
                jac  = 1d0
                const= 1d0
              endif

              !Find rhs of Gamma def.
              do ll=1,3
                do kk=1,3
                  select case (kk)
                  case(1)
                    dh = 2*g_def%dxh(ig)
                    if (alt_eom) then
                      jacp = g_def%gmetric%grid(igrid)%jac(i+1,j,k)
                      jacm = g_def%gmetric%grid(igrid)%jac(i-1,j,k)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp
     .                   *g_def%gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
                    gsupm = jacm
     .                   *g_def%gmetric%grid(igrid)%gsup(i-1,j,k,:,:)
                  case(2)
                    dh = 2*g_def%dyh(jg)
                    if (alt_eom) then
                      jacp = g_def%gmetric%grid(igrid)%jac(i,j+1,k)
                      jacm = g_def%gmetric%grid(igrid)%jac(i,j-1,k)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp
     .                   *g_def%gmetric%grid(igrid)%gsup(i,j+1,k,:,:)
                    gsupm = jacm
     .                   *g_def%gmetric%grid(igrid)%gsup(i,j-1,k,:,:)
                  case(3)
                    dh = 2*g_def%dzh(kg)
                    if (alt_eom) then
                      jacp = g_def%gmetric%grid(igrid)%jac(i,j,k+1)
                      jacm = g_def%gmetric%grid(igrid)%jac(i,j,k-1)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp
     .                   *g_def%gmetric%grid(igrid)%gsup(i,j,k+1,:,:)
                    gsupm = jacm
     .                   *g_def%gmetric%grid(igrid)%gsup(i,j,k-1,:,:)
                  end select
                    
                  table(ii,ll,kk) =
     .                  -gsub(ll,1)*(gsupp(1,ii)-gsupm(1,ii))/dh/jac
     .                  -gsub(ll,2)*(gsupp(2,ii)-gsupm(2,ii))/dh/jac
     .                  -gsub(ll,3)*(gsupp(3,ii)-gsupm(3,ii))/dh/jac
     .                  +const*delta(ii,ll)*(gamma(1,kk,1)
     .                                      +gamma(2,kk,2)
     .                                      +gamma(3,kk,3))

                  summ=0d0
                  do mm=1,3
                    do jj=1,3
                      summ=summ+gsub(ll,mm)*gsup(ii,jj)*gamma(mm,jj,kk)
                    enddo
                  enddo

                  table(ii,ll,kk) = table(ii,ll,kk) - summ
                enddo
              enddo

              !Postprocess alt_eom case (ii == 2 .and. bcond(1) == SP)
cc              if (alt_eom) then
cc                do kk=1,3
cc                  table(ii,ii,kk)
cc     $                 = 0.5*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk))
cc                enddo
cc              endif

            enddo

            !Symmetrize
            do kk=1,3
              do ll=1,3
                do ii=1,3
                  g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
     $                        = 0.5*( table(ii,ll,kk)+table(ii,kk,ll) )
                enddo
              enddo
            enddo

            enddo
          enddo
        enddo

cc        !Postprocessing for alt_eom case (ii == 2 .and. bcond(1) == SP) 
cc        if (bcSP()) then
cc
cc          gamma = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cc          jac   = g_def%gmetric%grid(igrid)%jac(i,j,k)
cc
cc          do ii=1,3             !cycle through Chistoffel symbols
cc            do ll=1,3
cc              do kk=1,3
cc                select case (kk)
cc                case(1)
cc                  dh = g_def%dxh(ig)
cc                  jacp = hlf*(jac+g_def%gmetric%grid(igrid)%jac(i+1,j,k))
cc                  jacm = hlf*(jac+g_def%gmetric%grid(igrid)%jac(i-1,j,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+g_def%gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+g_def%gmetric%grid(igrid)%gsup(i-1,j,k,:,:))
cc                case(2)
cc                  dh = g_def%dyh(jg)
cc                  jacp = hlf*(jac+g_def%gmetric%grid(igrid)%jac(i,j+1,k))
cc                  jacm = hlf*(jac+g_def%gmetric%grid(igrid)%jac(i,j-1,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+g_def%gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+g_def%gmetric%grid(igrid)%gsup(i,j-1,k,:,:))
cc                case(3)
cc                  dh = g_def%dzh(kg)
cc                  jacp = hlf*(jac+g_def%gmetric%grid(igrid)%jac(i,j,k+1))
cc                  jacm = hlf*(jac+g_def%gmetric%grid(igrid)%jac(i,j,k-1))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+g_def%gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+g_def%gmetric%grid(igrid)%gsup(i,j,k-1,:,:))
cc                end select
cc
cc                table(ii,kk,ll) =
cc     $                -gsub(ll,1)*(gsupp(1,ii)-gsupm(1,ii))/dh/jac
cc     .                -gsub(ll,2)*(gsupp(2,ii)-gsupm(2,ii))/dh/jac
cc     .                -gsub(ll,3)*(gsupp(3,ii)-gsupm(3,ii))/dh/jac
cc     .                +delta(ii,ll)*(gamma(1,kk,1)
cc     .                              +gamma(2,kk,2)
cc     .                              +gamma(3,kk,3))
cc
cc                summ=0d0
cc                do mm=1,3
cc                  do jj=1,3
cc                    summ= summ + gsub(ll,mm)*gsup(ii,jj)*gamma(mm,jj,kk)
cc                  enddo
cc                enddo
cc
cc                table(ii,kk,ll) = table(ii,kk,ll) - summ
cc
cc              enddo
cc            enddo
cc          enddo
cc
cc          !Symmetrize
cc          do ii=1,3             !cycle through Chistoffel symbols
cc            do ll=1,3
cc              do kk=ll,3
cc                table(ii,kk,ll) = hlf*(table(ii,kk,ll)+table(ii,ll,kk))
cc                table(ii,ll,kk) = table(ii,kk,ll)
cc              enddo
cc            enddo
cc          enddo
cc
cc          !Invert SI systems
cc          do ii=1,3
cc            do ll=1,3
cc              do kk=1,3
cc                if (ii==ll.and.ii/=kk) then
cc                  select case(ii)
cc                  case(1)
cc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,kk,2)-table(3,kk,3)
cc                  case(2)
cc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(3,kk,3)
cc                  case(3)
cc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(2,kk,2)
cc                  end select
cc                elseif(ii==kk .and.ii/=ll) then
cc                  select case(ii)
cc                  case(1)
cc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,2,ll)-table(3,3,ll)
cc                  case(2)
cc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(3,3,ll)
cc                  case(3)
cc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(2,2,ll)
cc                  end select
cc                elseif(ii==kk .and.ii==ll) then
cc                  g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = gamma(ii,kk,ll)
cc                else
cc                  g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = table(ii,kk,ll)
cc                endif
cccc                if (ii == ll) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     .                 = matmul(inv,table(ii,:,ii))
cccccc     $                 = hlf*(table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccccc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cccccc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccc                  case(2)
cccc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccccc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cccccc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccc                  case(3)
cccc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccccc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cccccc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccc                  end select
cccc                elseif (ii == kk) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccccc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cccccc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccc                  case(2)
cccc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccccc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cccccc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccc                  case(3)
cccc                    g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccccc     $                 =0.25*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll)
cccccc     $                       +table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccc                  end select
cccc                else
cccc                  g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cccc                endif
cc              enddo
cc            enddo
cc
cc          enddo

cc                if (ii /= ll .and. ii /= kk) then
cc                  g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = table(ii,ll,kk)
cc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cc                else
cc                  if (ii == ll) then
cc                    select case(ii)
cc                    case(1)
cc                      g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk))
cc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cc                    case(2)
cc                      g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk))
cc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cc                    case(3)
cc                      g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk))
cc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cc                    end select
cc                  else
cc                    select case(ii)
cc                    case(1)
cc                      g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll))
cc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cc                    case(2)
cc                      g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll))
cc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cc                    case(3)
cc                      g_def%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
ccc$$$     $                 = hlf*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll))
cc     $                 =0.25*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll)
cc     $                       +table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cc                    end select
cc
cc                  endif
cc                endif
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif

c     End program

      end subroutine hessianZeroForceCorr

      end subroutine findGridMetrics

c     topol_bc
c     #################################################################
      subroutine topol_bc(gdef,igrid,nxl,nyl,nzl,array)

c     -----------------------------------------------------------------
c     Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: gdef
        integer :: igrid,nxl,nyl,nzl
        real(8) :: array(0:nxl+1,0:nyl+1,0:nzl+1)

c     Local variables

#if defined(petsc)
        type(petsc_da_ctx),pointer :: da_ctx
#endif

c     External

c     Begin program

#if !defined(petsc)
        if (bcond(1) == PER .or. bcond(2) == PER) then
          array(0    ,:,:) = array(nxl,:,:)
          array(nxl+1,:,:) = array(1  ,:,:)
        endif

        if (bcond(3) == PER .or. bcond(4) == PER) then
          array(:,0    ,:) = array(:,nyl,:)
          array(:,nyl+1,:) = array(:,1  ,:)
        endif
        
        if (bcond(5) == PER .or. bcond(6) == PER) then
          array(:,:,0    ) = array(:,:,nzl)
          array(:,:,nzl+1) = array(:,:,1  )
        endif
#else
#if !defined(samrai)
        da_ctx => gdef%dactx(igrid)
        call fillPetscGhostCells(da_ctx,array)
#else
        write(*,*) 'WARNING: metrics in SAMRAI are not being'
     .            ,' communicated yet.'
#endif
#endif

c     End program

      end subroutine topol_bc

ccc     topol_bc_vrtx
ccc     #################################################################
cc      subroutine topol_bc_vrtx(igrid,nxl,nyl,nzl,array)
cc
ccc     -----------------------------------------------------------------
ccc     Enforce topological constraints on array
ccc     -----------------------------------------------------------------
cc
cc        implicit  none
cc
ccc     Call variables
cc
cc        integer    :: igrid,nxl,nyl,nzl
cc        real(8)    :: array(0:nxl+1,0:nyl+1,0:nzl+1)
cc
ccc     Local variables
cc
cc#if defined(petsc)
cc        type(petsc_da_ctx),pointer :: da_ctx
cc#endif
cc
ccc     External
cc
ccc     Begin program
cc
cc#if !defined(petsc)
cc        if (bcond(1) == PER .or. bcond(2) == PER) then
cccc          array(0    ,:,:) = array(nxl,:,:)
cccc          array(nxl+1,:,:) = array(1,:,:)
cc          array(1,:,:) = array(nxl+1,:,:)
cc        endif
cc
cc        if (bcond(3) == PER .or. bcond(4) == PER) then
cccc          array(:,0    ,:) = array(:,nyl,:)
cccc          array(:,nyl+1,:) = array(:,1,:)
cc          array(:,1,:) = array(:,nyl+1,:)
cc        endif
cc        
cc        if (bcond(5) == PER .or. bcond(6) == PER) then
cccc          array(:,:,0    ) = array(:,:,nzl)
cccc          array(:,:,nzl+1) = array(:,:,1)
cc          array(:,:,1) = array(:,:,nzl+1)
cc        endif
cc#else
cc#if !defined(samrai)
cc          write(*,*) 'WARNING: vertex metrics in PETSc are not being'
cc     .              ,' communicated yet.'
cccccc        if (igrid == 1) then
cccc          da_ctx => grid_params%dactx(igrid)
cccc          call fillPetscGhostCells(da_ctx,array)
cccccc        endif
cc#else
cc          write(*,*) 'WARNING: metrics in SAMRAI are not being'
cc     .              ,' communicated yet.'
cc#endif
cc#endif
cc
ccc     End program
cc
cc      end subroutine topol_bc_vrtx

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric(g_def,xcar,xcarv,igr,ierr)

c     -----------------------------------------------------------------
c     This routine calculates all grid metric quantities required for
c     the curvilinear representation of a set of PDE's: jacobian,
c     metric tensors, covariant and contravariant vectors, Christoffel
c     symbols of the second kind. All quantities are stored in
c     structure g_metric. There are two modes of computation:
c        * Analytical (numerical_grid=.false.)
c        * Numerical  (numerical_grid=.true.)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def
        real(8),optional,intent(IN),dimension(:,:,:,:) :: xcar,xcarv
        integer,optional :: igr,ierr

c     Local variables

        integer :: igrmin,igrmax,igrid,istat,lerr

c     Begin program

c     Checks

        map_present = PRESENT(xcar)
        vmap_present= PRESENT(xcarv)

        anal_map = checkAnalMapDatabase().and.(.not.map_present)

c     Define grid levels

        if (PRESENT(igr)) then
          igrmin = igr
          igrmax = igr
        else
          igrmin = 1
          igrmax = g_def%ngrid
        endif

c     Allocate arrays

        call allocateGridMetric(g_def%gmetric,g_def)

        if ( map_present) g_def%gmetric%grid(igrmin)%car   = xcar
        if (vmap_present) g_def%gmetric%grid(igrmin)%car_v = xcarv

c     Find grid metrics

        do igrid=igrmin,igrmax
          call findGridMetrics(g_def,igrid,.not.PRESENT(igr),ierr=ierr)
        enddo

c     Gather error flags

#if defined(petsc)
        if (PRESENT(ierr)) then
          lerr = ierr
          call MPI_Allreduce(lerr,ierr,1,MPI_INTEGER
     .                      ,MPI_MAX,g_def%mpi_comm,mpierr)
        endif
#endif

c     End program

      end subroutine defineGridMetric

      end module grid_metric_init

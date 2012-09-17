
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

        alt__eom=.false.
        if (bcSP()) alt__eom=.true.

      end function alt__eom

c     cross_product
c     #################################################################
      function cross_product(vec1,vec2) result(vec3)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: vec1(3),vec2(3),vec3(3)

c     Local variables

c     Begin program

        vec3(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
        vec3(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
        vec3(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

      end function cross_product

c     triple_product
c     #################################################################
      function triple_product(vec1,vec2,vec3) result(scalar)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: vec1(3),vec2(3),vec3(3),scalar

c     Local variables

c     Begin program

        scalar = dot_product(vec1,cross_product(vec2,vec3))

      end function triple_product

c     delta
c     #################################################################
      function delta(i,j)

        integer    :: i,j
        real(8)    :: delta

        if (i == j) then
          delta = 1d0
        else
          delta = 0d0
        endif

      end function delta

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
      function find_dr_1st(i,j,k,igrid,hex,hey,hez,gparams,x) result(dr)

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

        type(grid_mg_def),pointer :: gparams

c     Local variables

        integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
     .               ,igp,igm,jgp,jgm,kgp,kgm,hel,nx,ny,nz

        real(8)    :: carp(3),carm(3),dh(3),coeff

c     Begin program

        nx = gparams%nxv(igrid)
        ny = gparams%nyv(igrid)
        nz = gparams%nzv(igrid)

        ip=min(i+1,nx+1)
        im=max(i-1,0)
        jp=min(j+1,ny+1)
        jm=max(j-1,0)
        kp=min(k+1,nz+1)
        km=max(k-1,0)

        call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

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

          dh(1)= (gparams%xx(igp)-gparams%xx(igm))

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

          dh(2)= (gparams%yy(jgp)-gparams%yy(jgm))

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

          dh(3)= (gparams%zz(kgp)-gparams%zz(kgm))

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
      function find_dr_2nd(i,j,k,igrid,gparams,x) result(dr)

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

        type(grid_mg_def),pointer :: gparams

c     Local variables

        integer :: ip,im,i0,jp,jm,j0,kp,km,k0,ig,jg,kg
     .            ,igp,igm,jgp,jgm,kgp,kgm,nx,ny,nz

        real(8) :: carp(3),carm(3),car0(3),coeff

c     Begin program

        nx = gparams%nxv(igrid)
        ny = gparams%nyv(igrid)
        nz = gparams%nzv(igrid)

c     Exceptions

        !X
        i0 = min(max(i,1),nx)
        ip = i0+1
        im = i0-1

        call getMGmap(i0,j,k,igx,igy,igz,ig,jg,kg)

        igp = ig+1
        igm = ig-1

        carp = x(ip,j,k,:)
        car0 = x(i0,j,k,:)
        carm = x(im,j,k,:)

        dr(:,1) = dr_dxi(gparams%xx(igm),carm
     .                  ,gparams%xx(ig ),car0
     .                  ,gparams%xx(igp),carp
     .                  ,i-i0)

        !Y
        j0 = min(max(j,1),ny)
        jp = j0+1
        jm = j0-1

        call getMGmap(i,j0,k,igx,igy,igz,ig,jg,kg)

        jgp = jg+1
        jgm = jg-1

        carp = x(i,jp,k,:)
        car0 = x(i,j0,k,:)
        carm = x(i,jm,k,:)

        dr(:,2) = dr_dxi(gparams%yy(jgm),carm
     .                  ,gparams%yy(jg ),car0
     .                  ,gparams%yy(jgp),carp
     .                  ,j-j0)

        !Z
        k0 = min(max(k,1),nz)
        kp = k0+1
        km = k0-1

        call getMGmap(i,j,k0,igx,igy,igz,ig,jg,kg)

        kgp = kg+1
        kgm = kg-1

        carp = x(i,j,kp,:)
        car0 = x(i,j,k0,:)
        carm = x(i,j,km,:)

        dr(:,3) = dr_dxi(gparams%zz(kgm),carm
     .                  ,gparams%zz(kg ),car0
     .                  ,gparams%zz(kgp),carp
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
      subroutine coarsen_map(igc,gparams,xf,xc,bcs)

c     -----------------------------------------------------------------
c     Coarsen fine map xf(xi) to coarse map xc(xi). This is done
c     dimension by dimension, using first-order extrapolation.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igc
        real(8) :: xf(0:,0:,0:,:)
     .            ,xc(0:,0:,0:,:)

        type(grid_mg_def),pointer :: gparams

        integer,optional :: bcs(6,3)

c     Local variables

        integer :: ii,igf,nxf,nyf,nzf,nxc,nyc,nzc,neq

c     Begin program

        igf = igc - 1

        nxf = gparams%nxv(igf)
        nyf = gparams%nyv(igf)
        nzf = gparams%nzv(igf)

        nxc = gparams%nxv(igc)
        nyc = gparams%nyv(igc)
        nzc = gparams%nzv(igc)

        neq = size(xf,4)

        if (PRESENT(bcs)) then
          do ii=1,neq
            call crse_map_comp(ii,xf(:,:,:,ii)
     .                           ,xc(:,:,:,ii)
     .                           ,bcs=bcs(:,ii))
          enddo
        else
          do ii=1,3
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

        mgx = gparams%mg_ratio_x(igf)
        mgy = gparams%mg_ratio_y(igf)
        mgz = gparams%mg_ratio_z(igf)

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

        call map_bc(dim,igc,nxc,nyc,nzc,xc_j,x0_j,bcs=bcs)

      end subroutine crse_map_comp

      end subroutine coarsen_map

c     map_bc
c     #################################################################
      subroutine map_bc(dim,igrid,nxl,nyl,nzl,x_j,x0_j,iorder,bcs)

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
     .              ,arr0=x0_j,iorder=order)

c     End program

      end subroutine map_bc

c     findGridMetrics
c     #################################################################
      subroutine findGridMetrics(gparams,igrid,coarsen)

c     -----------------------------------------------------------------
c     Find grid metrics at grid level igrid.
c     -----------------------------------------------------------------

        use grid_debug

        implicit none

c     Input variables

        integer :: igrid
        type(grid_mg_def),pointer :: gparams

        logical :: coarsen

c     Local variables

        integer :: ii,nxf,nyf,nzf,nxc,nyc,nzc,nx,ny,nz,ig,jg,kg,m,l

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        ilo = 1
        jlo = 1
        klo = 1
        ihi = gparams%nxv(igrid)
        jhi = gparams%nyv(igrid)
        khi = gparams%nzv(igrid)

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
                gparams%gmetric%grid(igrid)%car  (i,j,k,:)
     .                    = map     (i,j,k,igx,igy,igz)
                gparams%gmetric%grid(igrid)%jac  (i,j,k)
     .                    = jacobian(i,j,k,igx,igy,igz)
                gparams%gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                    = g_sub   (i,j,k,igx,igy,igz)
                gparams%gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                    = g_sup   (i,j,k,igx,igy,igz)
                gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .                    = christ_2knd(i,j,k,igx,igy,igz)
                do l=1,3
                  gparams%gmetric%grid(igrid)%cov(i,j,k,l,:)
     .                    = covariantVector    (l,i,j,k,igx,igy,igz)
                  gparams%gmetric%grid(igrid)%cnv(i,j,k,l,:)
     .                    = contravariantVector(l,i,j,k,igx,igy,igz)
                enddo

                !Cell volume
                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

                if (gparams%nxgl(igrid) > 1)
     .               gparams%gmetric%grid(igrid)%dvol(i,j,k) =
     .                           gparams%dxh(ig)

                if (gparams%nygl(igrid) > 1)
     .               gparams%gmetric%grid(igrid)%dvol(i,j,k) =
     .               gparams%gmetric%grid(igrid)%dvol(i,j,k) *
     .                           gparams%dyh(jg)

                if (gparams%nzgl(igrid) > 1)
     .               gparams%gmetric%grid(igrid)%dvol(i,j,k) =
     .               gparams%gmetric%grid(igrid)%dvol(i,j,k) *
     .                           gparams%dzh(kg)

                gparams%gmetric%grid(igrid)%dvol(i,j,k)=
     .                           gparams%gmetric%grid(igrid)%jac(i,j,k)
     .                          *gparams%gmetric%grid(igrid)%dvol(i,j,k)
              enddo
            enddo
          enddo

          !Vertex metrics
          do k = klom,khi
            do j = jlom,jhi
              do i = ilom,ihi
                gparams%gmetric%grid(igrid)%jac_v(i,j,k)
     .                    = jacobian(i,j,k,igx,igy,igz,vertex=.true.)
                gparams%gmetric%grid(igrid)%gsub_v(i,j,k,:,:)
     .                    = g_sub   (i,j,k,igx,igy,igz,vertex=.true.)
                gparams%gmetric%grid(igrid)%gsup_v (i,j,k,:,:)
     .                    = g_sup   (i,j,k,igx,igy,igz,vertex=.true.)
                do l=1,3
                  gparams%gmetric%grid(igrid)%cov_v(i,j,k,l,:)
     .                    = covariantVector    (l,i,j,k,igx,igy,igz
     .                                         ,vertex=.true.)
                  gparams%gmetric%grid(igrid)%cnv_v(i,j,k,l,:)
     .                    = contravariantVector(l,i,j,k,igx,igy,igz
     .                                         ,vertex=.true.)
                enddo

              enddo
            enddo
          enddo

          !Zero force condition on Christoffel symbols (only on finest grid)
          if (igrid == 1) call hessianZeroForceCorr(igrid)

c       NUMERICAL QUANTITIES ***********************************

        else

          numerical_grid = .true.

          call findCellMap

          call findVertexMap

          call findVertexMetrics

          call findCellMetrics

        endif

c     End program

      contains

c     findCellMap
c     #################################################################
      subroutine findCellMap

c     -----------------------------------------------------------------
c     Finds Cartesian map at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: xf,xc
        real(8),pointer,dimension(:) :: xg,yg,zg

        integer,save :: c_rat_x,c_rat_y,c_rat_z

c     Begin program

        xc => gparams%gmetric%grid(igrid)%car

        if (map_present) then

          if (igrid > 1 .and. coarsen) then  !Coarsen map if at coarse grid level
            xf => gparams%gmetric%grid(igrid-1)%car

            call coarsen_map(igrid,gparams,xf,xc)

            nullify(xf)
          endif

        else

          if (checkAnalMapDatabase()) then

            do k = klom,khip
              do j = jlom,jhip
                do i = ilom,ihip
                  xc(i,j,k,:) = map(i,j,k,igx,igy,igz)
                enddo
              enddo
            enddo

cc          elseif (coords == 'pck') then
cc
cc            !X coordinate
cc            allocate(xg(0:gparams%nxgl(igx)+1))
cc
cc            if (igx == 1) then
cc              c_rat_x = 1
cc            else
cc              c_rat_x = c_rat_x*gparams%mg_ratio_x(igx-1)
cc            endif
cc
cc            call findGlobalGrid(gparams%gxmin,gparams%gxmax
cc     .                         ,bcond(1),bcond(2),xg
cc     .                         ,gparams%g_pack(1),c_rat_x)
cc
cc            do i = ilom,ihip
cc              xc(i,:,:,1) = xg(gparams%ilo(igx)+i-1)
cc            enddo
cc
cc            deallocate(xg)
cc
cc            !Y coordinate
cc            allocate(yg(0:gparams%nygl(igy)+1))
cc
cc            if (igy == 1) then
cc              c_rat_y = 1
cc            else
cc              c_rat_y = c_rat_y*gparams%mg_ratio_y(igy-1)
cc            endif
cc
cc            call findGlobalGrid(gparams%gymin,gparams%gymax
cc     .                         ,bcond(3),bcond(4),yg
cc     .                         ,gparams%g_pack(2),c_rat_y)
cc
cc            do j = jlom,jhip
cc              xc(:,j,:,2) = yg(gparams%jlo(igy)+j-1)
cc            enddo
cc
cc            deallocate(yg)
cc
cc            !Z coordinate
cc            allocate(zg(0:gparams%nzgl(igz)+1))
cc
cc            if (igz == 1) then
cc              c_rat_z = 1
cc            else
cc              c_rat_z = c_rat_z*gparams%mg_ratio_z(igz-1)
cc            endif
cc
cc            call findGlobalGrid(gparams%gzmin,gparams%gzmax
cc     .                         ,bcond(5),bcond(6),zg
cc     .                         ,gparams%g_pack(3),c_rat_z)
cc
cc            do k = klom,khip
cc              xc(:,:,k,3) = zg(gparams%klo(igz)+k-1)
cc            enddo
cc
cc            deallocate(zg)

cc#if defined(vmec)
cc          elseif (coords == 'vmc') then
cc
cc            vmec_met = .false.   !Whether to use VMEC metric quantities or not
cc            vmec_mets = vmec_met .and. (igrid == 1)
cc
cc            nxf=ihi-ilo+1
cc            nyf=jhi-jlo+1
cc            nzf=khi-klo+1
cc
cc            call vmec_map(igrid,nxf,nyf,nzf,xc,vmec_mets) !External
cc
cc            if (vmec_mets) then
cc              call vmec_metrics(igrid,nxf,nyf,nzf         !External
cc     .                         ,gparams%gmetric%grid(igrid)%jac
cc     .                         ,gparams%gmetric%grid(igrid)%gsub
cc     .                         ,gparams%gmetric%grid(igrid)%gsup)
cc            endif
cc#endif
          
          else
            call pstop('findCellMap','Map not available')
          endif

        endif

c     End program

      end subroutine findCellMap

c     findVertexMap
c     #################################################################
      subroutine findVertexMap

c     -----------------------------------------------------------------
c     Finds Cartesian map at vertices.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: xc,xv

c     Begin program

        xv => gparams%gmetric%grid(igrid)%car_v

        if ((.not.vmap_present).or.(igrid > 1)) then

          xc => gparams%gmetric%grid(igrid)%car !Cell center map

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
      subroutine findVertexMetrics

c     -----------------------------------------------------------------
c     Finds metrics at vertices (i+1/2,j+1/2,k+1/2) (stored at (i,j,k))
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        integer :: he,imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz

        real(8) :: r(3,3),cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),jac,ijac
     .            ,gamma(3,3,3),d2r(3,3)

        real(8),pointer,dimension(:,:,:,:,:) :: dr
        real(8),pointer,dimension(:,:,:,:)   :: xc

c     Begin program

        imin = ilom
        imax = ihi

        jmin = jlom
        jmax = jhi

        kmin = klom
        kmax = khi

        he = 1

c     Local grid size

        nx = ihi-ilo+1
        ny = jhi-jlo+1
        nz = khi-klo+1

c     Evaluate dr=dx/dxi vectors at vertices

        xc => gparams%gmetric%grid(igrid)%car

        allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))

        dr = 0d0

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              dr(i,j,k,:,:)
     .             = find_dr(i,j,k,igrid,he,he,he,gparams,xc)
            enddo
          enddo
        enddo

c     Enforce topological constraints

cc        do i=1,3
cc          do j=1,3
cc            call topol_bc_vrtx(igrid,nx,ny,nz,dr(:,:,:,i,j))
cc          enddo
cc        enddo

c     Evaluate grid quantities (at vertices)

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(:,1),r(:,2),r(:,3))
              gparams%gmetric%grid(igrid)%jac_v(i,j,k) = jac

              if (jac <= 0d0 .and. (.not.isSP2(i+1,igrid,ibc=1))) then
                write (*,*) 'Negative Jacobian in half mesh'
                write (*,*) 'Grid level',igrid
                write (*,*) 'Grid position=',i,j,k
                write (*,*) 'dr/dxi_1=',r(:,1)
                write (*,*) 'dr/dxi_2=',r(:,2)
                write (*,*) 'dr/dxi_3=',r(:,3)
                write (*,*) 'Jacobian=',jac
                if (.not.check_grid)
     .             call pstop('findVertexMetrics','Inverted jacobian')
              elseif (isSP2(i+1,igrid,ibc=1)) then
cc              if (isSP2(i+1,igrid,ibc=1)) then
                ijac = 0d0
              else
                ijac = 1d0/jac
              endif

              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
              cnv(:,:) = r(:,:)*ijac

              gparams%gmetric%grid(igrid)%cnv_v(i,j,k,:,:)
     .                 =transpose(cnv)

              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
              cov(:,1) = cross_product(r(:,2),r(:,3))
              cov(:,2) = cross_product(r(:,3),r(:,1))
              cov(:,3) = cross_product(r(:,1),r(:,2))

              cov = cov*ijac

              gparams%gmetric%grid(igrid)%cov_v(i,j,k,:,:)
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
              gparams%gmetric%grid(igrid)%gsub_v(i,j,k,:,:) = gsub
              gparams%gmetric%grid(igrid)%gsup_v(i,j,k,:,:) = gsup

            enddo
          enddo
        enddo

c     Christoffel symbols (at cell centers; BCs enforced in findCellMetrics)

        do k = kmin+1,kmax
          do j = jmin+1,jmax
            do i = imin+1,imax

              cov = 0.125*(
     .             gparams%gmetric%grid(igrid)%cov_v(i-1,j-1,k-1,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i  ,j-1,k-1,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i-1,j  ,k-1,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i  ,j  ,k-1,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i-1,j-1,k  ,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i  ,j-1,k  ,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i-1,j  ,k  ,:,:)
     .            +gparams%gmetric%grid(igrid)%cov_v(i-1,j-1,k  ,:,:))

              do l=1,3
                do m=1,3
                  do n=1,3
                    d2r=find_dr(i-1,j-1,k-1,igrid,1,1,1,gparams
     .                         ,dr(:,:,:,:,m))
                    gamma(l,m,n) = dot_product(d2r(:,n),cov(l,:))
                  enddo
                enddo
                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
              enddo

              gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma

            enddo
          enddo
        enddo

c     End program

        deallocate(dr)
        nullify(xc)

      end subroutine findVertexMetrics

ccc     findCellMetrics
ccc     #################################################################
cc      subroutine findCellMetrics
cc
ccc     -----------------------------------------------------------------
ccc     Finds metrics at cell centers.
ccc     -----------------------------------------------------------------
cc
cc        implicit  none
cc
ccc     Call variables
cc
ccc     Local variables
cc
cc        integer :: he,imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
cc
cc        real(8) :: r(3,3),car0(3),carp(3),carm(3),jac,ijac
cc     .            ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
cc     .            ,gamma(3,3,3),mag,dhh,dhm,dhp,d2r(3,3)
cc
cc        real(8),pointer,dimension(:,:,:,:,:) :: dr
cc        real(8),pointer,dimension(:,:,:,:)   :: xc
cc
ccc     Begin program
cc
cc        imin = ilom
cc        imax = ihip
cc
cc        jmin = jlom
cc        jmax = jhip
cc
cc        kmin = klom
cc        kmax = khip
cc
ccc     Local grid size
cc
cc        nx = ihi-ilo+1
cc        ny = jhi-jlo+1
cc        nz = khi-klo+1
cc
ccc     Evaluate dr=dx/dxi vectors
cc
cc        xc => gparams%gmetric%grid(igrid)%car
cc
cc        allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))
cc
cc        dr = 0d0
cc
cc        if (igrid == 1) then  !Second order for finest grid
cc          do k = kmin,kmax
cc            do j = jmin,jmax
cc              do i = imin,imax
cc                dr(i,j,k,:,:) = find_dr(i,j,k,igrid,gparams,xc)
cc              enddo
cc            enddo
cc          enddo
cc        else                  !First order in coarser grids (to avoid negative jacobians)
cc          do k = kmin,kmax
cc            do j = jmin,jmax
cc              do i = imin,imax
cc                dr(i,j,k,:,:) = find_dr(i,j,k,igrid,0,0,0,gparams,xc)
cc              enddo
cc            enddo
cc          enddo
cc        endif
cc
ccc     Enforce topological constraints on dr (including SP)
cc
cc        do i=1,3
cc          do j=1,3
cc            call topol_bc(igrid,nx,ny,nz,dr(:,:,:,i,j))
cc            if (bcSP()) then
cc              igxbc = igrid
cc              igybc = igrid
cc              igzbc = igrid
cc
cc              nxbc = nx
cc              nybc = ny
cc              nzbc = nz
cc              !dr behaves like a covariant vector
cc              call vectorSingularBC(j,dr(:,:,:,i,j),.true.)
cc            endif
cc          enddo
cc        enddo
cc
ccc     Evaluate grid quantities
cc
cc        do k = kmin,kmax
cc          do j = jmin,jmax
cc            do i = imin,imax
cc              r = dr(i,j,k,:,:)
cc
cc              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cc              jac = triple_product(r(:,1),r(:,2),r(:,3))
cc              gparams%gmetric%grid(igrid)%jac(i,j,k) = jac
cc
cc              if (jac <= 0d0 .and.(.not.isSP2(i+1,igrid,ibc=1))) then
cc                write (*,*) 'Grid level',igrid
cc                write (*,*) 'Grid position=',i,j,k
cc                write (*,*) 'dr/dxi_1=',r(:,1)
cc                write (*,*) 'dr/dxi_2=',r(:,2)
cc                write (*,*) 'dr/dxi_3=',r(:,3)
cc                write (*,*) 'Jacobian=',jac
cc                if (.not.check_grid)
cc     .             call pstop('findCellMetrics','Inverted jacobian')
cc              endif
cc
cc              ijac = 1d0/jac
cc
cc              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
cc              cnv(:,:) = r(:,:)*ijac
cc
cc              gparams%gmetric%grid(igrid)%cnv(i,j,k,:,:) =transpose(cnv)
cc
cc              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
cc              cov(:,1) = cross_product(r(:,2),r(:,3))
cc              cov(:,2) = cross_product(r(:,3),r(:,1))
cc              cov(:,3) = cross_product(r(:,1),r(:,2))
cc
cc              cov = cov*ijac
cc
cc              gparams%gmetric%grid(igrid)%cov(i,j,k,:,:) =transpose(cov)
cc
cc              !Metric tensors <<<<<<<<<<<<<<<<<<<<<<<<
cc              do m=1,3
cc                do l=m,3
cc                  gsub(l,m) = jac*dot_product(cnv(:,l),cnv(:,m))
cc                  gsub(m,l) = gsub(l,m) !Symmetry
cc                  gsup(l,m) = jac*dot_product(cov(:,l),cov(:,m))
cc                  gsup(m,l) = gsup(l,m) !Symmetry
cc                enddo
cc              enddo
cc              gparams%gmetric%grid(igrid)%gsub(i,j,k,:,:) = gsub
cc              gparams%gmetric%grid(igrid)%gsup(i,j,k,:,:) = gsup
cc
cc              !Cell volume <<<<<<<<<<<<<<<<<<
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              if (gparams%nxgl(igrid) > 1)
cc     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) =
cc     .                         gparams%dxh(ig)
cc
cc              if (gparams%nygl(igrid) > 1)
cc     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) =
cc     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) *
cc     .                         gparams%dyh(jg)
cc
cc              if (gparams%nzgl(igrid) > 1)
cc     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) =
cc     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) *
cc     .                         gparams%dzh(kg)
cc
cc              gparams%gmetric%grid(igrid)%dvol(i,j,k)=
cc     .                           gparams%gmetric%grid(igrid)%jac(i,j,k)
cc     .                          *gparams%gmetric%grid(igrid)%dvol(i,j,k)
cccc     .                          *gparams%dxh(ig)
cccc     .                          *gparams%dyh(jg)
cccc     .                          *gparams%dzh(kg)
cc
cc              !Christoffel symbols <<<<<<<<<<<<<<<<<<
cc              do l=1,3
cc                do m=1,3
cc                  do n=1,3
cc                    d2r = find_dr(i,j,k,igrid,gparams,dr(:,:,:,:,m))
cc                    gamma(l,m,n) = dot_product(d2r(:,n),cov(:,l))
cc                  enddo
cc                enddo
cc                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
cc              enddo
cc
cc              gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Zero-force condition on Christoffle symbols 
cc        if (igrid == 1) call hessianZeroForceCorr(igrid)
cc
cc        !Enforce topological constraints on Christoffel symbols
cc        do i=1,3
cc          do j=1,3
cc            do k=1,3
cc              call topol_bc(igrid,nx,ny,nz
cc     .           ,gparams%gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
cc            enddo
cc          enddo
cc        enddo
cc
cc        deallocate(dr)
cc        nullify(xc)
cc
ccc     End program
cc
cc      end subroutine findCellMetrics

c     findCellMetrics
c     #################################################################
      subroutine findCellMetrics

c     -----------------------------------------------------------------
c     Finds metrics at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        integer :: he,imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
     .            ,bc_dr(6)

        real(8) :: r(3,3),car0(3),carp(3),carm(3),jac,ijac
     .            ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .            ,gamma(3,3,3),mag,dhh,dhm,dhp,d2r(3,3)

        real(8),pointer,dimension(:,:,:,:,:) :: dr
        real(8),pointer,dimension(:,:,:,:)   :: xv!,xc
        real(8),pointer,dimension(:,:,:)     :: zeros

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

        xv => gparams%gmetric%grid(igrid)%car_v

c     Evaluate dr=dx/dxi vectors at cell centers from vertex map

        allocate(dr   (ilom:ihip,jlom:jhip,klom:khip,3,3)
     .          ,zeros(ilom:ihip,jlom:jhip,klom:khip))

        dr = 0d0
        zeros = 0d0

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi
              dr(i,j,k,:,:)=find_dr(i-1,j-1,k-1,igrid,1,1,1,gparams,xv)
            enddo
          enddo
        enddo

cc        !OLD
cc        xc => gparams%gmetric%grid(igrid)%car
cc
cc        do k = kmin,kmax
cc          do j = jmin,jmax
cc            do i = imin,imax
cc              dr(i,j,k,:,:) = find_dr(i,j,k,igrid,gparams,xc)
cc            enddo
cc          enddo
cc        enddo
cc
cc        nullify(xc)

c     Enforce topological constraints on dr (including SP) and extrapolate otherwise

        bc_dr = bcond
        where (bc_dr == DEF) bc_dr = EXT
        where (bc_dr == SP ) bc_dr = DEF

        do j=1,3
          do i=1,3
            call setBC(0,nx,ny,nz,dr(:,:,:,i,j),zeros,bc_dr
     .                ,igrid,igrid,igrid,iorder=1)

cc            !Periodic BCs (OLD)
cc            call topol_bc(igrid,nx,ny,nz,dr(:,:,:,i,j))

            !Singular point BCs
            if (bcSP()) then
              !dr behaves like a covariant vector
              call vectorSingularBC(j,dr(:,:,:,i,j),.true.)
            endif

          enddo
        enddo

c     Evaluate grid quantities

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(:,1),r(:,2),r(:,3))
              gparams%gmetric%grid(igrid)%jac(i,j,k) = jac

              if (jac <= 0d0 .and.(.not.isSP2(i+1,igrid,ibc=1))) then
                write (*,*) 'Grid level',igrid
                write (*,*) 'Grid position=',i,j,k
                write (*,*) 'dr/dxi_1=',r(:,1)
                write (*,*) 'dr/dxi_2=',r(:,2)
                write (*,*) 'dr/dxi_3=',r(:,3)
                write (*,*) 'Jacobian=',jac
                if (.not.check_grid)
     .             call pstop('findCellMetrics','Inverted jacobian')
              endif

              ijac = 1d0/jac

              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
              cnv(:,:) = r(:,:)*ijac

              gparams%gmetric%grid(igrid)%cnv(i,j,k,:,:) =transpose(cnv)

              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
              cov(:,1) = cross_product(r(:,2),r(:,3))
              cov(:,2) = cross_product(r(:,3),r(:,1))
              cov(:,3) = cross_product(r(:,1),r(:,2))

              cov = cov*ijac

              gparams%gmetric%grid(igrid)%cov(i,j,k,:,:) =transpose(cov)

              !Metric tensors <<<<<<<<<<<<<<<<<<<<<<<<
              do m=1,3
                do l=m,3
                  gsub(l,m) = jac*dot_product(cnv(:,l),cnv(:,m))
                  gsub(m,l) = gsub(l,m) !Symmetry
                  gsup(l,m) = jac*dot_product(cov(:,l),cov(:,m))
                  gsup(m,l) = gsup(l,m) !Symmetry
                enddo
              enddo
              gparams%gmetric%grid(igrid)%gsub(i,j,k,:,:) = gsub
              gparams%gmetric%grid(igrid)%gsup(i,j,k,:,:) = gsup

              !Cell volume <<<<<<<<<<<<<<<<<<
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              if (gparams%nxgl(igrid) > 1)
     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) =
     .                         gparams%dxh(ig)

              if (gparams%nygl(igrid) > 1)
     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) =
     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) *
     .                         gparams%dyh(jg)

              if (gparams%nzgl(igrid) > 1)
     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) =
     .             gparams%gmetric%grid(igrid)%dvol(i,j,k) *
     .                         gparams%dzh(kg)

              gparams%gmetric%grid(igrid)%dvol(i,j,k)=
     .                           gparams%gmetric%grid(igrid)%jac(i,j,k)
     .                          *gparams%gmetric%grid(igrid)%dvol(i,j,k)

cc              !Christoffel symbols (OLD) <<<<<<<<<<<<<<<<<<
cc              do l=1,3
cc                do m=1,3
cc                  do n=1,3
cc                    d2r = find_dr(i,j,k,igrid,gparams,dr(:,:,:,:,m))
cc                    gamma(l,m,n) = dot_product(d2r(:,n),cov(:,l))
cc                  enddo
cc                enddo
cc                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
cc              enddo
cc
cc              gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma

            enddo
          enddo
        enddo

c     Enforce zero-force condition and topological constraints on Christoffle symbols 

        if (igrid == 1) call hessianZeroForceCorr(igrid)

        do i=1,3
          do j=1,3
            do k=1,3
cc              !OLD
cc              call topol_bc(igrid,nx,ny,nz
cc     .           ,gparams%gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
              call setBC(0,nx,ny,nz
     .                  ,gparams%gmetric%grid(igrid)%Gamma(:,:,:,i,j,k)
     .                  ,zeros,bc_dr,igrid,igrid,igrid,iorder=1)
            enddo
          enddo
        enddo

c     End program

        deallocate(dr,zeros)

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

        integer    :: i,j,k,ii,ll,kk,mm,jj,ig,jg,kg
        real(8)    :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),table(3,3,3),gsub(3,3)
     .               ,gsup(3,3),gsupp(3,3),gsupm(3,3),gamma(3,3,3)
        logical    :: alt_eom

c     Begin program

        do k = klo,khi
          do j = jlo,jhi
            do i = ilo,ihi

            gamma = gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
            gsub  = gparams%gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsup  = gparams%gmetric%grid(igrid)%gsup(i,j,k,:,:)

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            do ii=1,3 !cycle through Chistoffel symbols

              alt_eom = (ii == 2 .and. alt__eom())

              if (alt_eom) then
                jac  = gparams%gmetric%grid(igrid)%jac(i,j,k)
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
                    dh = 2*gparams%dxh(ig)
                    if (alt_eom) then
                      jacp = gparams%gmetric%grid(igrid)%jac(i+1,j,k)
                      jacm = gparams%gmetric%grid(igrid)%jac(i-1,j,k)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp
     .                   *gparams%gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
                    gsupm = jacm
     .                   *gparams%gmetric%grid(igrid)%gsup(i-1,j,k,:,:)
                  case(2)
                    dh = 2*gparams%dyh(jg)
                    if (alt_eom) then
                      jacp = gparams%gmetric%grid(igrid)%jac(i,j+1,k)
                      jacm = gparams%gmetric%grid(igrid)%jac(i,j-1,k)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp
     .                   *gparams%gmetric%grid(igrid)%gsup(i,j+1,k,:,:)
                    gsupm = jacm
     .                   *gparams%gmetric%grid(igrid)%gsup(i,j-1,k,:,:)
                  case(3)
                    dh = 2*gparams%dzh(kg)
                    if (alt_eom) then
                      jacp = gparams%gmetric%grid(igrid)%jac(i,j,k+1)
                      jacm = gparams%gmetric%grid(igrid)%jac(i,j,k-1)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp
     .                   *gparams%gmetric%grid(igrid)%gsup(i,j,k+1,:,:)
                    gsupm = jacm
     .                   *gparams%gmetric%grid(igrid)%gsup(i,j,k-1,:,:)
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
            do ii=1,3
              do ll=1,3             !cycle through tensor components
                do kk=1,3
                  gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
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
cc          gamma = gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cc          jac   = gparams%gmetric%grid(igrid)%jac(i,j,k)
cc
cc          do ii=1,3             !cycle through Chistoffel symbols
cc            do ll=1,3
cc              do kk=1,3
cc                select case (kk)
cc                case(1)
cc                  dh = gparams%dxh(ig)
cc                  jacp = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i+1,j,k))
cc                  jacm = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i-1,j,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i-1,j,k,:,:))
cc                case(2)
cc                  dh = gparams%dyh(jg)
cc                  jacp = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j+1,k))
cc                  jacm = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j-1,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j-1,k,:,:))
cc                case(3)
cc                  dh = gparams%dzh(kg)
cc                  jacp = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j,k+1))
cc                  jacm = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j,k-1))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j,k-1,:,:))
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
cc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,kk,2)-table(3,kk,3)
cc                  case(2)
cc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(3,kk,3)
cc                  case(3)
cc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(2,kk,2)
cc                  end select
cc                elseif(ii==kk .and.ii/=ll) then
cc                  select case(ii)
cc                  case(1)
cc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,2,ll)-table(3,3,ll)
cc                  case(2)
cc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(3,3,ll)
cc                  case(3)
cc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(2,2,ll)
cc                  end select
cc                elseif(ii==kk .and.ii==ll) then
cc                  gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = gamma(ii,kk,ll)
cc                else
cc                  gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = table(ii,kk,ll)
cc                endif
cccc                if (ii == ll) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     .                 = matmul(inv,table(ii,:,ii))
cccccc     $                 = hlf*(table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccccc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cccccc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccc                  case(2)
cccc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccccc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cccccc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccc                  case(3)
cccc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccccc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cccccc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccc                  end select
cccc                elseif (ii == kk) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccccc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cccccc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccc                  case(2)
cccc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccccc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cccccc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccc                  case(3)
cccc                    gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccccc     $                 =0.25*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll)
cccccc     $                       +table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccc                  end select
cccc                else
cccc                  gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cccc                endif
cc              enddo
cc            enddo
cc
cc          enddo

cc                if (ii /= ll .and. ii /= kk) then
cc                  gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = table(ii,ll,kk)
cc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cc                else
cc                  if (ii == ll) then
cc                    select case(ii)
cc                    case(1)
cc                      gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk))
cc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cc                    case(2)
cc                      gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk))
cc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cc                    case(3)
cc                      gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk))
cc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cc                    end select
cc                  else
cc                    select case(ii)
cc                    case(1)
cc                      gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll))
cc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cc                    case(2)
cc                      gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll))
cc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cc                    case(3)
cc                      gparams%gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
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
      subroutine topol_bc(igrid,nxl,nyl,nzl,array)

c     -----------------------------------------------------------------
c     Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

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
cc        if (igrid == 1) then
          da_ctx => grid_params%dactx(igrid)
          call fillPetscGhostCells(da_ctx,array)
cc        endif
#else
          write(*,*) 'WARNING: metrics in SAMRAI are not being'
     .              ,' communicated yet.'
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
      subroutine defineGridMetric(gparams,xcar,xcarv,igr)

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

        type(grid_mg_def),pointer :: gparams
        real(8),optional,intent(IN),dimension(:,:,:,:) :: xcar,xcarv
        integer,optional :: igr

c     Local variables

        integer :: igrmin,igrmax,igrid,istat

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
          igrmax = gparams%ngrid
        endif

c     Allocate arrays

        call allocateGridMetric(gparams%gmetric,gparams)

        if ( map_present) gparams%gmetric%grid(igrmin)%car   = xcar
        if (vmap_present) gparams%gmetric%grid(igrmin)%car_v = xcarv

c     Find grid metrics

        do igrid=igrmin,igrmax
          call findGridMetrics(gparams,igrid,.not.PRESENT(igr))
        enddo

c     End program

      end subroutine defineGridMetric

      end module grid_metric_init

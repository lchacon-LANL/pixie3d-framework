
c module grid_metric_init
c #####################################################################
      module grid_metric_init

c ---------------------------------------------------------------------
c     This module allocates and fills the grid metric structure gmetric.
c     It contains the following routines: 
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

        integer,private :: i,j,k,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .                    ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp

        logical,private :: map_present,vmap_present

        INTERFACE find_dr
          module procedure find_dr_1st,find_dr_2nd !4th
        END INTERFACE

        private :: findGridmetrics,find_dr

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

c     find_d2r
c     #####################################################################
      function find_d2r(g_def,igrid,x) result(d2r)

c     ---------------------------------------------------------------------
c     This routine finds d2x/d_xi/d_xj tensor everywhere on the mesh
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: x(0:,0:,0:,:)
     .            ,d2r(0:size(x,1)-1
     .                ,0:size(x,2)-1
     .                ,0:size(x,3)-1,3,3,3)

c     Local variables

        integer :: i,j,k,nx,ny,nz,bcs(6,3),m
        real(8) :: dr(0:size(x,1)-1
     .               ,0:size(x,2)-1
     .               ,0:size(x,3)-1,3,3)
        
c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

c     Evaluate dr=dx/dxi vectors at vertices

        do k = 0,nz
          do j = 0,ny
            do i = 0,nx
              dr(i,j,k,:,:) = find_dr(g_def,i,j,k,igrid,1,1,1,x)
            enddo
          enddo
        enddo

c     Evaluate d2r=d2x/dxi/dxj hessian at cell centers

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              do m= 1,3
                d2r(i,j,k,:,:,m)=find_dr(g_def,i-1,j-1,k-1,igrid,1,1,1
     .                                  ,dr(:,:,:,:,m))
              enddo
            enddo
          enddo
        enddo

      end function find_d2r

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

c     Fill ghost cells

        call map_bc(g_def,dim,igc,nxc,nyc,nzc,xc_j,bcs=bcs)

      end subroutine crse_map_comp

      end subroutine coarsen_map

c     map_bc
c     #################################################################
      subroutine map_bc(g_def,dim,igrid,nxl,nyl,nzl,x_j,iorder,bcs)

c     -----------------------------------------------------------------
c     Impose BCs on map component x_j using boundary of x0_j.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: dim,igrid,nxl,nyl,nzl

        real(8) :: x_j (0:nxl+1,0:nyl+1,0:nzl+1)

        integer,optional,intent(IN) :: iorder
        integer,optional :: bcs(6)

c     Local variables

        integer :: bcnd(6),order
        real(8) :: x0_j(0:nxl+1,0:nyl+1,0:nzl+1)

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

c     Prepare boundaries

        x0_j = x_j

        if ((nxl == 1).and.(dim == 1)) then
          i=0
          x0_j(i,:,:) = 0.5*(x0_j(i,:,:)+x0_j(i+1,:,:))
          i=nxl+1
          x0_j(i,:,:) = 0.5*(x0_j(i,:,:)+x0_j(i-1,:,:))
        endif

        if ((nyl == 1).and.(dim == 2)) then
          j=0
          x0_j(:,j,:) = 0.5*(x0_j(:,j,:)+x0_j(:,j+1,:))
          j=nyl+1
          x0_j(:,j,:) = 0.5*(x0_j(:,j,:)+x0_j(:,j-1,:))
        endif

        if ((nzl == 1).and.(dim == 3)) then
          k=0
          x0_j(:,:,k) = 0.5*(x0_j(:,:,k)+x0_j(:,:,k+1))
          k=nzl+1
          x0_j(:,:,k) = 0.5*(x0_j(:,:,k)+x0_j(:,:,k-1))
        endif

c     Impose BCs

        call setMGBC(g_def,0,1,nxl,nyl,nzl,igrid,x_j,bcnd
     .              ,arr0=x0_j,iorder=order)

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

        lerr = 0

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

c       ANALYTICAL QUANTITIES **********************************

        if ((.not.numerical_grid).and.anal_map) then

          !Cell-centered metrics
          do k = 0,nz+1
            do j = 0,ny+1
              do i = 0,nx+1
                g_def%gmetric%grid(igrid)%car  (i,j,k,:)
     .                    = map     (g_def,i,j,k,igrid,igrid,igrid)
                g_def%gmetric%grid(igrid)%jac  (i,j,k)
     .                    = jacobian(g_def,i,j,k,igrid,igrid,igrid)
                g_def%gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                    = g_sub   (g_def,i,j,k,igrid,igrid,igrid)
                g_def%gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                    = g_sup   (g_def,i,j,k,igrid,igrid,igrid)
                g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .                    = christ_2knd(g_def,i,j,k,igrid,igrid,igrid)
                do l=1,3
                  g_def%gmetric%grid(igrid)%cov(i,j,k,l,:)
     .                = covariantVector(g_def,l,i,j,k,igrid,igrid,igrid)
                  g_def%gmetric%grid(igrid)%cnv(i,j,k,l,:)
     .            = contravariantVector(g_def,l,i,j,k,igrid,igrid,igrid)
                enddo

                !Cell volume
                call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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
          do k = 0,nz+1
            do j = 0,ny+1
              do i = 0,nx+1
                g_def%gmetric%grid(igrid)%jac_v(i,j,k)
     .           = jacobian(g_def,i,j,k,igrid,igrid,igrid,vertex=.true.)
                g_def%gmetric%grid(igrid)%gsub_v(i,j,k,:,:)
     .           = g_sub   (g_def,i,j,k,igrid,igrid,igrid,vertex=.true.)
                g_def%gmetric%grid(igrid)%gsup_v (i,j,k,:,:)
     .           = g_sup   (g_def,i,j,k,igrid,igrid,igrid,vertex=.true.)
                do l=1,3
                  g_def%gmetric%grid(igrid)%cov_v(i,j,k,l,:)
     .                 = covariantVector(g_def,l,i,j,k,igrid,igrid,igrid
     .                                  ,vertex=.true.)
                  g_def%gmetric%grid(igrid)%cnv_v(i,j,k,l,:)
     .             = contravariantVector(g_def,l,i,j,k,igrid,igrid,igrid
     .                                   ,vertex=.true.)
                enddo

              enddo
            enddo
          enddo

          !Zero force condition on Christoffel symbols (only on finest grid)
          if (igrid == 1) call hessianZeroForceCorr(g_def,igrid)

c       NUMERICAL QUANTITIES ***********************************

        else

          numerical_grid = .true.

          call findCellMap(g_def,igrid,coarsen)

          call findVertexMap(g_def,igrid)

          call findCellMetrics(g_def,igrid)

          call findVertexMetrics(g_def,igrid)

          call Christoffel(g_def,igrid)

        endif

c     Gather error flags

        if (return_error) then
#if defined(petsc)
          ierr = pmax(lerr,mpi_comm=g_def%mpi_comm)
#else
          ierr = lerr
#endif
        endif

c     End program

      contains

c     findVertexMetrics
c     #################################################################
      subroutine findVertexMetrics(g_def,igrid)

c     -----------------------------------------------------------------
c     Finds metrics at vertices (i+1/2,j+1/2,k+1/2) (stored at (i,j,k))
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        
c     Local variables

        integer :: he,nx,ny,nz

        real(8) :: r(3,3),cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),jac,ijac

        real(8),pointer,dimension(:,:,:,:,:) :: dr
        real(8),pointer,dimension(:,:,:,:)   :: xc

c     Begin program

        he = 1

c     Local grid size

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

c     Evaluate dr=dx/dxi vectors at vertices

        allocate(dr(0:nx+1,0:ny+1,0:nz+1,3,3))

        xc => g_def%gmetric%grid(igrid)%car

        dr = 0d0

        do k = 0,nz
          do j = 0,ny
            do i = 0,nx
              dr(i,j,k,:,:) = find_dr(g_def,i,j,k,igrid,he,he,he,xc)
            enddo
          enddo
        enddo

c     Impose BCs on dr

        call dr_bc_vrtx(g_def,igrid,dr)

c     Evaluate grid quantities (at vertices)

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nx+1
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(:,1),r(:,2),r(:,3))
              g_def%gmetric%grid(igrid)%jac_v(i,j,k) = jac

              if (jac <= 0d0
     .           .and. (.not.isSP2(g_def,i+1,igrid))) then
                if (.not.(check_grid.or.return_error)) then
                  write (*,*) 'Negative or zero Jacobian in half mesh'
                  write (*,*) 'Grid level',igrid
                  write (*,*) 'Processor rank=',my_rank
                  write (*,*) 'Grid position=',i,j,k
                  write (*,*) 'dr/dxi_1=',r(:,1)
                  write (*,*) 'dr/dxi_2=',r(:,2)
                  write (*,*) 'dr/dxi_3=',r(:,3)
                  write (*,*) 'Jacobian=',jac

                  call sstop(0,'findVertexMetrics','Inverted Jacobian')
                elseif (return_error) then
                  lerr = 1 ;  exit
                endif
              elseif (isSP2(g_def,i+1,igrid)) then
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

c     End program

        deallocate(dr)
        nullify(xc)

      end subroutine findVertexMetrics

c     findCellMetrics
c     #################################################################
      subroutine findCellMetrics(g_def,igrid)

c     -----------------------------------------------------------------
c     Finds metrics at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        
c     Local variables

        integer :: he,nx,ny,nz,bcs(6)

        real(8) :: r(3,3),car0(3),carp(3),carm(3),jac,ijac
     .            ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .            ,gamma(3,3,3),mag,dhh,dhm,dhp,d2r(3,3)

        real(8),pointer,dimension(:,:,:,:,:) :: dr
        real(8),pointer,dimension(:,:,:,:)   :: xc

c     Begin program

c     Local grid size

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

c     Evaluate dr=dx/dxi vectors

        xc => g_def%gmetric%grid(igrid)%car

        allocate(dr(0:nx+1,0:ny+1,0:nz+1,3,3))

        dr = 0d0

        if (igrid == 1) then  !Second order for finest grid
          do k = 0,nz+1
            do j = 0,ny+1
              do i = 0,nx+1
                dr(i,j,k,:,:) = find_dr(g_def,i,j,k,igrid,xc)
              enddo
            enddo
          enddo
        else                  !First order in coarser grids (to avoid negative jacobians)
          do k = 0,nz+1
            do j = 0,ny+1
              do i = 0,nx+1
                dr(i,j,k,:,:) = find_dr(g_def,i,j,k,igrid,0,0,0,xc)
              enddo
            enddo
          enddo
        endif

c     Apply BCs on dr (including SP)

        call dr_bc(g_def,igrid,dr)

c     Evaluate rest of grid quantities

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nx+1
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(:,1),r(:,2),r(:,3))
              g_def%gmetric%grid(igrid)%jac(i,j,k) = jac

              if (jac <= 0d0
     .           .and.(.not.isSP2(g_def,i+1,igrid))) then
                if (.not.check_grid .and. (.not.return_error)) then
                  write (*,*) 'Negative Jacobian in integer mesh'
                  write (*,*) 'Grid level',igrid
                  write (*,*) 'Grid position=',i,j,k
                  write (*,*) 'dr/dxi_1=',r(:,1)
                  write (*,*) 'dr/dxi_2=',r(:,2)
                  write (*,*) 'dr/dxi_3=',r(:,3)
                  write (*,*) 'Jacobian=',jac

                  call sstop(0,'findCellMetrics','Inverted jacobian')
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
              call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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

            enddo
          enddo
        enddo

c     End program

        deallocate(dr)
        nullify(xc)

      end subroutine findCellMetrics

c     Christoffel
c     #################################################################
      subroutine Christoffel(g_def,igrid)

c     -----------------------------------------------------------------
c     Compute Christoffel symbols of the second kind at cell centers
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        
c     Local variables

        integer :: nx,ny,nz

        integer :: i,j,k,l,m,n,he,bcs(6,3)
        real(8) :: r(3,3),cov(3,3),gamma(3,3,3)

        real(8),pointer,dimension(:,:,:,:,:,:) :: d2r
        real(8),pointer,dimension(:,:,:,:,:)   :: dr
        real(8),pointer,dimension(:,:,:,:)     :: xc

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        xc => g_def%gmetric%grid(igrid)%car

c     NEW WAY: compute directly second-order derivatives

        allocate(d2r(0:nx+1,0:ny+1,0:nz+1,3,3,3))

        d2r = find_d2r(g_def,igrid,xc)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx

              cov = g_def%gmetric%grid(igrid)%cov(i,j,k,:,:)
              
              do l=1,3
                do n=1,3
                  do m=1,3
                    gamma(l,m,n)=dot_product(d2r(i,j,k,:,m,n),cov(l,:))
                  enddo
                enddo
                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
              enddo

              g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma

            enddo
          enddo
        enddo

        deallocate(d2r)

c$$$c     LEGACY WAY: cell-centered derivatives
c$$$
c$$$        allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))
c$$$
c$$$        dr = 0d0
c$$$
c$$$        do k = klom,khip
c$$$          do j = jlom,jhip
c$$$            do i = ilom,ihip
c$$$              dr(i,j,k,:,:) = find_dr(g_def,i,j,k,igrid,xc)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$
c$$$        call dr_bc(g_def,igrid,dr)
c$$$
c$$$        do k = klom,khip
c$$$          do j = jlom,jhip
c$$$            do i = ilom,ihip
c$$$
c$$$              cov = g_def%gmetric%grid(igrid)%cov(i,j,k,:,:)
c$$$              
c$$$              do l=1,3
c$$$                do m=1,3
c$$$                  d2_r(:,:,m) = find_dr(g_def,i,j,k,igrid,dr(:,:,:,:,m))
c$$$                  do n=1,3
c$$$                    gamma(l,n,m) = dot_product(d2_r(:,n,m),cov(l,:))
c$$$                  enddo
c$$$                enddo
c$$$                gamma(l,:,:) =0.5*(gamma(l,:,:)+transpose(gamma(l,:,:))) !Symmetrize
c$$$              enddo
c$$$
c$$$              g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
c$$$
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$        
c$$$        deallocate(dr)
        
        nullify(xc)

c     Enforce zero-force condition on Christoffle symbols
        
        if (igrid == 1) call hessianZeroForceCorr(g_def,igrid)

c     Apply BCs to Christoffel symbols

        !Regularize
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
               g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) =
     .              g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .             *g_def%gmetric%grid(igrid)%jac  (i,j,k)
             enddo
          enddo
        enddo

        !Apply BCs
        bcs(:,1) = bcond
        bcs(:,2) = bcond
        bcs(:,3) = bcond
        where (bcs == DEF.or.bcs == SP) bcs = EXT
        do k=1,3
          do j=1,3
            do i=1,3
              call setMGBC(g_def,0,1,nx,ny,nz,igrid
     .                   ,g_def%gmetric%grid(igrid)%Gamma(:,:,:,i:i,j,k)
     .                   ,bcs(:,1:1),iorder=2)
            enddo
          enddo
        enddo

        !Transform back
        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nx+1
               g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) =
     .              g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .             /g_def%gmetric%grid(igrid)%jac  (i,j,k)
             enddo
          enddo
        enddo

c     End program

      end subroutine Christoffel
      
c     hessianZeroForceCorr
c     #################################################################
      subroutine hessianZeroForceCorr(g_def,igrid)

c     -----------------------------------------------------------------
c     Postprocess Christoffel symbol to satisfy zero-force condition
c     to machine accuracy
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        
c     Local variables

        integer :: i,j,k,ii,ll,kk,mm,jj,ig,jg,kg,nx,ny,nz
        real(8) :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8) :: hess(3,3,3),table(3,3,3),gamma(3,3,3)
     .            ,gsub(3,3),gsup(3,3),gsupp(3,3),gsupm(3,3)
        logical :: alt_eom

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx

            gamma = g_def%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
            gsub  = g_def%gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsup  = g_def%gmetric%grid(igrid)%gsup(i,j,k,:,:)

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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

c     findCellMap
c     #################################################################
      subroutine findCellMap(g_def,igrid,coarsen)

c     -----------------------------------------------------------------
c     Finds Cartesian map at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer :: igrid
        type(grid_mg_def),pointer :: g_def
        logical :: coarsen

c     Local variables

        integer :: nx,ny,nz
        real(8),pointer,dimension(:,:,:,:) :: xf,xc

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        xc => g_def%gmetric%grid(igrid)%car

        if (map_present) then

          if (igrid > 1 .and. coarsen) then  !Coarsen map if at coarse grid level
            xf => g_def%gmetric%grid(igrid-1)%car

            call coarsen_map(igrid,g_def,xf,xc)

            nullify(xf)
          endif

        else

          if (checkAnalMapDatabase()) then

            do k = 0,nz+1
              do j = 0,ny+1
                do i = 0,nx+1
                  xc(i,j,k,:) = map(g_def,i,j,k,igrid,igrid,igrid)
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
      subroutine findVertexMap(g_def,igrid)

c     -----------------------------------------------------------------
c     Finds Cartesian map at vertices.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer :: igrid
        type(grid_mg_def),pointer :: g_def

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: xc,xv

        integer :: ig1,jg1,kg1,nx,ny,nz

        real(8) :: x1,y1,z1

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        xv => g_def%gmetric%grid(igrid)%car_v

        if ((.not.vmap_present).or.(igrid > 1)) then

          if (checkAnalMapDatabase()) then

            do k = 0,nz
              do j = 0,ny
                do i = 0,nx
                  xv(i,j,k,:) = map(g_def,i,j,k,igrid,igrid,igrid
     .                             ,vertex=.true.)
                enddo
              enddo
            enddo

          else

            xc => g_def%gmetric%grid(igrid)%car !Cell center map

c$$$            call setupSplines(g_def,igrid,2,bcond)
c$$$
c$$$            call splineX(xc)
c$$$
c$$$            do k = 0,nz
c$$$              do j = 0,ny
c$$$                do i = 0,nx
c$$$                  call getCurvilinearCoordinates(g_def,i,j,k,igrid,igrid,igrid
c$$$     .                                          ,ig1,jg1,kg1,x1,y1,z1
c$$$     .                                          ,vertex=.true.)
c$$$
c$$$                  call evalX(x1,y1,z1,xv(i,j,k,1)
c$$$     .                               ,xv(i,j,k,2)
c$$$     .                               ,xv(i,j,k,3))
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$
c$$$            call killSplines

            xv(0:nx,0:ny,0:nz,:) = 0.125*(xc(0:nx  ,0:ny  ,0:nz  ,:)
     .                                   +xc(0:nx  ,1:ny+1,0:nz  ,:)
     .                                   +xc(0:nx  ,0:ny  ,1:nz+1,:)
     .                                   +xc(0:nx  ,1:ny+1,1:nz+1,:)
     .                                   +xc(1:nx+1,0:ny  ,0:nz  ,:)
     .                                   +xc(1:nx+1,1:ny+1,0:nz  ,:)
     .                                   +xc(1:nx+1,0:ny  ,1:nz+1,:)
     .                                   +xc(1:nx+1,1:ny+1,1:nz+1,:))
          endif

        endif

c     End program

        nullify(xv,xc)

      end subroutine findVertexMap

c     dr_bc
c     #################################################################
      subroutine dr_bc(g_def,igrid,dr)

c     -----------------------------------------------------------------
c     Apply BCs on dr/dxi vector
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: dr(:,:,:,:,:)

c     Local variables

        integer :: bcond_sv(6),i,j,bcs(6,3),nx,ny,nz

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

c     Save topological BC config

        bcond_sv = bcond

c     Modify BC setup

        bcs(:,1) = bcond
        bcs(:,2) = bcond
        bcs(:,3) = bcond

c$$$        where (bcs == SP) bcs = DEF

c$$$        if (g_def%nglx == 1 .and. bcPER(1)) bcs(1:2,:) = DEF
c$$$        if (g_def%ngly == 1 .and. bcPER(2)) bcs(3:4,:) = DEF
c$$$        if (g_def%nglz == 1 .and. bcPER(3)) bcs(5:6,:) = DEF

        !Special case for helical geometry when limiting angular domain
        if(coords=='hel'.and.g_def%params(1)/=1) then
          bcond(3:4) = DEF
          bcs(3:4,:) = EXT !DEF
          bcs(1  ,:) = DEF
        endif

c     Impose BCs

        do i=1,3  !Cycle over Cartesian coordinates
          call setBC_ez(g_def,igrid,dr(:,:,:,i,:),.false.,bcs,2)
        enddo

c     Recover topological BC config

        bcond = bcond_sv

c     End program

      end subroutine dr_bc

c     dr_bc_vrtx
c     #################################################################
      subroutine dr_bc_vrtx(g_def,igrid,dr)

c     -----------------------------------------------------------------
c     Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: dr(0:,0:,0:,:,:)

c     Local variables

        integer :: bcond_sv(6),bcs(6,3),i

c     Begin program

c     Save topological BC config

        bcond_sv = bcond

c     Start with topological BCs

        bcs(:,1) = bcond
        bcs(:,2) = bcond
        bcs(:,3) = bcond

c     Set up i=0, j=0, k=0 boundaries (do NOTHING: dr is trustworthy there)

        where (bcs(1,:) /= PER) bcs(1,:) = DEF
        where (bcs(3,:) /= PER) bcs(3,:) = DEF
        where (bcs(5,:) /= PER) bcs(5,:) = DEF

c     Set up i=imax, j=jmax, k=kmax boundaries

        where (bcs(2,:) == DEF) bcs(2,:) = EXT
        where (bcs(4,:) == DEF) bcs(4,:) = EXT 
        where (bcs(6,:) == DEF) bcs(6,:) = EXT 

        !Set BCs to do NEU across face
        if (bcond_sv(2) == FSYM .or. bcond_sv(2) == SYM) bcond(2) = SYM
        if (bcond_sv(4) == FSYM .or. bcond_sv(4) == SYM) bcond(4) = SYM
        if (bcond_sv(6) == FSYM .or. bcond_sv(6) == SYM) bcond(6) = SYM

        where (bcs(2,:) == FSYM .or. bcs(2,:) == SYM) bcs(2,:) = NEU
        where (bcs(4,:) == FSYM .or. bcs(4,:) == SYM) bcs(4,:) = NEU
        where (bcs(6,:) == FSYM .or. bcs(6,:) == SYM) bcs(6,:) = NEU

c     Impose BCs

        do i=1,3  !Cycle over Cartesian coordinates
          call imposeBConVector_stg_ez(g_def,igrid,dr(:,:,:,i,:)
     .                                ,.false.,bcs,2)
        enddo

c     Recover topological BC config

        bcond = bcond_sv

c     End program

      end subroutine dr_bc_vrtx

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric(g_def,xcar,xcarv,igr,ierr)

c     -----------------------------------------------------------------
c     This routine calculates all grid metric quantities required for
c     the curvilinear representation of a set of PDE's: jacobian,
c     metric tensors, covariant and contravariant vectors, Christoffel
c     symbols of the second kind. All quantities are stored in
c     structure g_metric.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def
        real(8),optional,intent(IN),dimension(:,:,:,:) :: xcar,xcarv
        integer,optional :: igr,ierr

c     Local variables

        integer :: igrmin,igrmax,igrid,istat,lerr

c     Begin program

        lerr = 0
        
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

c     Find grid metrics at all grid levels

        do igrid=igrmin,igrmax
          call findGridMetrics(g_def,igrid,.not.PRESENT(igr),ierr=lerr)
c$$$          if (my_rank == 0) write (*,*) 'igrid=',igrid,
c$$$     $                            'findGridMetrics error',lerr
          if (lerr /= 0) exit
        enddo

        !Reset grid level to last usable
        if (lerr /= 0 .and. (.not.PRESENT(igr))) g_def%ngrid = igrid - 1

        if (g_def%ngrid < 1.and.(.not.PRESENT(ierr))) then
          call pstop("defineGridMetric"
     .              ,"Cannot evaluate metrics on top mesh")
        endif

c     Gather error flags

        if (PRESENT(ierr)) ierr = lerr

c     End program

      end subroutine defineGridMetric

c     coarsenGridMetrics
c     #################################################################
      subroutine coarsenGridMetrics(g_def,igc,ierr)

c     -----------------------------------------------------------------
c     Coasens grid metrics to grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igc
        type(grid_mg_def),pointer :: g_def

        integer, optional :: ierr

c     Local variables

        integer :: igf,ic,jc,kc,if,jf,kf,icg,jcg,kcg
     .            ,nxc,nyc,nzc,nxf,nyf,nzf,l,m,bcs(6)

c     Begin program

        !From analytical mapping
        if ((.not.numerical_grid).and.anal_map) then

          call findGridMetrics(g_def,igc,.false.,ierr=ierr)

        !From numerical mapping
        else

          call findCellMap(g_def,igc,.true.)

          call findVertexMap(g_def,igc)

          call get_topol_bcs(bcs)
          where(bcs == DEF) bcs = EXT

          igf = igc - 1

          nxf = g_def%nxv(igf)
          nyf = g_def%nyv(igf)
          nzf = g_def%nzv(igf)

          nxc = g_def%nxv(igc)
          nyc = g_def%nyv(igc)
          nzc = g_def%nzv(igc)

c         Coarsen volumes and jacobians (needed for coarsening all else)

          !Coarsen cell volumes (with ghost cells)
          call coarsen(igc,g_def%gmetric%grid(igf)%dvol
     .                    ,g_def%gmetric%grid(igc)%dvol
     .                ,.false.)

          !Calculate jacobian
          do kc = 1,nzc
            do jc = 1,nyc
              do ic = 1,nxc
                call getMGmap(g_def,ic,jc,kc,igc,igc,igc,icg,jcg,kcg)

                g_def%gmetric%grid(igc)%jac(ic,jc,kc) =
     .                      g_def%gmetric%grid(igc)%dvol(ic,jc,kc)
     .                                    /(g_def%dxh(icg)
     .                                     *g_def%dyh(jcg)
     .                                     *g_def%dzh(kcg))
              enddo
            enddo
          enddo

          call imposeBConScalar_ez(g_def,igc
     .              ,g_def%gmetric%grid(igc)%jac,bcs
     .              ,iorder=1)

c         Coarsen everything else (with ghost cells)

          !Gsup
          do l=1,3
            do m=l,3
              call coarsen(igc
     .                    ,g_def%gmetric%grid(igf)%gsup(:,:,:,l,m)
     .                    *g_def%gmetric%grid(igf)%jac
     .                    ,g_def%gmetric%grid(igc)%gsup(:,:,:,l,m)
     .                    ,.true.)
              g_def%gmetric%grid(igc)%gsup(:,:,:,l,m) =
     .             g_def%gmetric%grid(igc)%gsup(:,:,:,l,m)
     .            /g_def%gmetric%grid(igc)%jac
              g_def%gmetric%grid(igc)%gsup(:,:,:,m,l)
     .             = g_def%gmetric%grid(igc)%gsup(:,:,:,l,m) !Symmetry
            enddo
          enddo

          !Gsub
          do l=1,3
            do m=l,3
              call coarsen(igc
     .                    ,g_def%gmetric%grid(igf)%gsub(:,:,:,l,m)
     .                    *g_def%gmetric%grid(igf)%jac
     .                    ,g_def%gmetric%grid(igc)%gsub(:,:,:,l,m)
     .                    ,.true.)
              g_def%gmetric%grid(igc)%gsub(:,:,:,l,m) =
     .             g_def%gmetric%grid(igc)%gsub(:,:,:,l,m)
     .            /g_def%gmetric%grid(igc)%jac
              g_def%gmetric%grid(igc)%gsub(:,:,:,m,l)
     .             = g_def%gmetric%grid(igc)%gsub(:,:,:,l,m) !Symmetry
            enddo
          enddo

          !Gamma
          do m=1,3
            do n=m,3
              do l=1,3
                call coarsen(igc
     .                      ,g_def%gmetric%grid(igf)%Gamma(:,:,:,l,m,n)
     .                      *g_def%gmetric%grid(igf)%jac
     .                      ,g_def%gmetric%grid(igc)%Gamma(:,:,:,l,m,n)
     .                      ,.true.)
                g_def%gmetric%grid(igc)%Gamma(:,:,:,l,m,n) =
     .             g_def%gmetric%grid(igc)%Gamma(:,:,:,l,m,n)
     .            /g_def%gmetric%grid(igc)%jac
              enddo
              g_def%gmetric%grid(igc)%Gamma(:,:,:,:,n,m)
     .             = g_def%gmetric%grid(igc)%Gamma(:,:,:,:,m,n) !Symmetry
            enddo
          enddo

          !Cov
          do l=1,3
            do m=1,3
              call coarsen(igc
     .                    ,g_def%gmetric%grid(igf)%cov(:,:,:,l,m)
     .                    *g_def%gmetric%grid(igf)%jac
     .                    ,g_def%gmetric%grid(igc)%cov(:,:,:,l,m)
     .                    ,.true.)
              g_def%gmetric%grid(igf)%cov(:,:,:,l,m) =
     .              g_def%gmetric%grid(igf)%cov(:,:,:,l,m)
     .             /g_def%gmetric%grid(igf)%jac
            enddo
          enddo

          !Cnv
          do l=1,3
            do m=1,3
              call coarsen(igc
     .                    ,g_def%gmetric%grid(igf)%cnv(:,:,:,l,m)
     .                    *g_def%gmetric%grid(igf)%jac
     .                    ,g_def%gmetric%grid(igc)%cnv(:,:,:,l,m)
     .                    ,.true.)
              g_def%gmetric%grid(igf)%cnv(:,:,:,l,m) =
     .              g_def%gmetric%grid(igf)%cnv(:,:,:,l,m)
     .             /g_def%gmetric%grid(igf)%jac
            enddo
          enddo

        endif

c     End program

      contains

c     coarsen
c     #################################################################
      subroutine coarsen(igc,arrayf,arrayc,volf)

c     -----------------------------------------------------------------
c     Coarsen arrayf to arrayc at grid level igc.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igc
        real(8) :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1)
     .            ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1)
        logical :: volf

c     Local variables

        integer :: igf,ic,jc,kc,if,jf,kf,mgx,mgy,mgz
     .            ,i_off,j_off,k_off

        real(8) :: volc,ac

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

              volc   = 0d0
              ac     = 0d0

              do kf=mgz*((kc+k_off)-1)+1-k_off,mgz*(kc+k_off)-k_off
                do jf=mgy*((jc+j_off)-1)+1-j_off,mgy*(jc+j_off)-j_off
                  do if=mgx*((ic+i_off)-1)+1-i_off,mgx*(ic+i_off)-i_off

                    if (volf) then
                      ac   = ac
     .                     + g_def%gmetric%grid(igf)%dvol(if,jf,kf)
     .                      *arrayf(if,jf,kf)
                      volc = volc
     .                     + g_def%gmetric%grid(igf)%dvol(if,jf,kf)
                    else
                      ac   = ac   + arrayf(if,jf,kf)
                    endif

                  enddo
                enddo
              enddo

              if (volf) ac = ac/volc

              arrayc(ic,jc,kc) = ac
            enddo
          enddo
        enddo

c     Fill ghost cells by 2nd-order extrapolation

        call imposeBConScalar_ez(g_def,igc,arrayc,bcs,iorder=1)

      end subroutine coarsen

      end subroutine coarsenGridMetrics

c     find_dr_1st
c     #####################################################################
      function find_dr_1st(g_def,i,j,k,igrid,hex,hey,hez,x) result(dr)

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

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid,hex,hey,hez
        real(8)    :: dr(3,3),x(0:,0:,0:,:)

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

        call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

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
      function find_dr_2nd(g_def,i,j,k,igrid,x) result(dr)

c     ---------------------------------------------------------------------
c     This routine finds dx/d_xi tensor at the point (i,j,k)
c     from the position vector x. 
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: dr(3,3),x(0:,0:,0:,:)

c     Local variables

        integer :: ip,im,i0,jp,jm,j0,kp,km,k0,ig,jg,kg
     .            ,igp,igm,jgp,jgm,kgp,kgm,nx,ny,nz

        real(8) :: carp(3),carm(3),car0(3),coeff

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        !X
        i0 = min(max(i,1),nx)
        ip = i0+1
        im = i0-1

        call getMGmap(g_def,i0,j,k,igrid,igrid,igrid,ig,jg,kg)

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

        call getMGmap(g_def,i,j0,k,igrid,igrid,igrid,ig,jg,kg)

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

        call getMGmap(g_def,i,j,k0,igrid,igrid,igrid,ig,jg,kg)

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
          dr = (-rp* dx1/(dx2*(dx1+dx2))
     .          +r0*(1d0/dx1+1d0/dx2)
     .          -rm*(1d0/dx1+1d0/(dx1+dx2)))
        case(1)  !Right boundary
          dr = ( rm* dx2/(dx1*(dx1+dx2))
     .          -r0*(1d0/dx1+1d0/dx2)
     .          +rp*(1d0/dx2+1d0/(dx1+dx2)))
        case(0)
          dr = ( rp* dx1/(dx2*(dx1+dx2))
     .          +r0*(1d0/dx1-1d0/dx2)
     .          -rm* dx2/(dx1*(dx1+dx2)))
        case default
          write (*,*) 'Error in dr_dxi'
          stop
        end select

      end function dr_dxi

      end function find_dr_2nd

c     find_dr_4th
c     #####################################################################
      function find_dr_4th(g_def,i,j,k,igrid,x) result(dr)

c     ---------------------------------------------------------------------
c     This routine finds dx/d_xi tensor at the point (i,j,k)
c     from the position vector x. 
c     ---------------------------------------------------------------------

        use oned_int
      
        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: dr(3,3),x(0:,0:,0:,:)

c     Local variables

        integer :: ipp,ip,i0,im,imm
     .            ,jpp,jp,j0,jm,jmm
     .            ,kpp,kp,k0,km,kmm
     .            ,igpp,igp,ig0,igm,igmm
     .            ,jgpp,jgp,jg0,jgm,jgmm
     .            ,kgpp,kgp,kg0,kgm,kgmm
     .            ,ig,jg,kg,nx,ny,nz

        real(8) :: carpp(3),carp(3),car0(3),carm(3),carmm(3),coeff

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

c     Exceptions

        !X
        if (nx > 3) then
          
          i0 = min(max(i,2),nx-1)
          ipp = i0+2
          ip  = i0+1
          im  = i0-1
          imm = i0-2

          call getMGmap(g_def,i0,j,k,igrid,igrid,igrid,ig0,jg,kg)

          igpp = ig0+2
          igp  = ig0+1
          igm  = ig0-1
          igmm = ig0-2

          carpp = x(ipp,j,k,:)
          carp  = x(ip ,j,k,:)
          car0  = x(i0 ,j,k,:)
          carm  = x(im ,j,k,:)
          carmm = x(imm,j,k,:)
        
          dr(:,1) = quartic_interp(g_def%xx(igmm),carmm
     .                            ,g_def%xx(igm ),carm
     .                            ,g_def%xx(ig0 ),car0
     .                            ,g_def%xx(igp ),carp
     .                            ,g_def%xx(igpp),carpp
     .                            ,g_def%xx(ig),1)

        else
          i0 = min(max(i,1),nx)
          ip = i0+1
          im = i0-1

          call getMGmap(g_def,i0,j,k,igrid,igrid,igrid,ig0,jg,kg)

          igp = ig0+1
          igm = ig0-1

          carp = x(ip,j,k,:)
          car0 = x(i0,j,k,:)
          carm = x(im,j,k,:)

          dr(:,1) = dr_dxi(g_def%xx(igm),carm
     .                    ,g_def%xx(ig0),car0
     .                    ,g_def%xx(igp),carp
     .                    ,i-i0)
        endif
          
        !Y
        if (ny > 3) then
          j0 = min(max(j,2),ny-1)
          jpp = j0+2
          jp  = j0+1
          jm  = j0-1
          jmm = j0-2

          call getMGmap(g_def,i,j0,k,igrid,igrid,igrid,ig,jg0,kg)

          jgpp = jg0+2
          jgp  = jg0+1
          jgm  = jg0-1
          jgmm = jg0-2

          carpp = x(i,jpp,k,:)
          carp  = x(i,jp ,k,:)
          car0  = x(i,j0 ,k,:)
          carm  = x(i,jm ,k,:)
          carmm = x(i,jmm,k,:)
        
          dr(:,2) = quartic_interp(g_def%yy(jgmm),carmm
     .                            ,g_def%yy(jgm ),carm
     .                            ,g_def%yy(jg0 ),car0
     .                            ,g_def%yy(jgp ),carp
     .                            ,g_def%yy(jgpp),carpp
     .                            ,g_def%yy(jg),1)

        else
          j0 = min(max(j,1),ny)
          jp = j0+1
          jm = j0-1

          call getMGmap(g_def,i,j0,k,igrid,igrid,igrid,ig,jg0,kg)

          jgp = jg0+1
          jgm = jg0-1

          carp = x(i,jp,k,:)
          car0 = x(i,j0,k,:)
          carm = x(i,jm,k,:)

          dr(:,2) = dr_dxi(g_def%yy(jgm),carm
     .                    ,g_def%yy(jg0),car0
     .                    ,g_def%yy(jgp),carp
     .                    ,j-j0)
        endif
        
        !Z
        if (nz > 3) then
          k0 = min(max(k,2),nz-1)
          kpp = k0+2
          kp  = k0+1
          km  = k0-1
          kmm = k0-2

          call getMGmap(g_def,i,j,k0,igrid,igrid,igrid,ig,jg,kg0)

          kgpp = kg0+2
          kgp  = kg0+1
          kgm  = kg0-1
          kgmm = kg0-2

          carpp = x(i,j,kpp,:)
          carp  = x(i,j,kp ,:)
          car0  = x(i,j,k0 ,:)
          carm  = x(i,j,km ,:)
          carmm = x(i,j,kmm,:)
        
          dr(:,2) = quartic_interp(g_def%zz(kgmm),carmm
     .                            ,g_def%zz(kgm ),carm
     .                            ,g_def%zz(kg0 ),car0
     .                            ,g_def%zz(kgp ),carp
     .                            ,g_def%zz(kgpp),carpp
     .                            ,g_def%zz(kg),1)
        else
          k0 = min(max(k,1),nz)
          kp = k0+1
          km = k0-1

          call getMGmap(g_def,i,j,k0,igrid,igrid,igrid,ig,jg,kg0)

          kgp = kg0+1
          kgm = kg0-1

          carp = x(i,j,kp,:)
          car0 = x(i,j,k0,:)
          carm = x(i,j,km,:)

          dr(:,3) = dr_dxi(g_def%zz(kgm),carm
     .                    ,g_def%zz(kg0),car0
     .                    ,g_def%zz(kgp),carp
     .                    ,k-k0)
        endif
        
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
          dr = (-rp* dx1/(dx2*(dx1+dx2))
     .          +r0*(1d0/dx1+1d0/dx2)
     .          -rm*(1d0/dx1+1d0/(dx1+dx2)))
        case(1)  !Right boundary
          dr = ( rm* dx2/(dx1*(dx1+dx2))
     .          -r0*(1d0/dx1+1d0/dx2)
     .          +rp*(1d0/dx2+1d0/(dx1+dx2)))
        case(0)
          dr = ( rp* dx1/(dx2*(dx1+dx2))
     .          +r0*(1d0/dx1-1d0/dx2)
     .          -rm* dx2/(dx1*(dx1+dx2)))
        case default
          write (*,*) 'Error in dr_dxi'
          stop
        end select

      end function dr_dxi

      end function find_dr_4th

      end module grid_metric_init

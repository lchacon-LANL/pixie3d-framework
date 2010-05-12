
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
c           -> zeroForceCheck
c        * coarsenGridMetrics
c           -> coarsen
c        * topol_bc
c        * map_bc
c        * defineGridMetric
c           -> numericalGridCoeffsCheck
c           -> graphicsCheck
c ---------------------------------------------------------------------

        use grid_anal_map

        use grid_mpi

        use setMGBC_interface

        implicit none

        integer :: funit
        logical :: check_grid,numerical_grid

c     PRIVATE VARIABLES

        integer   ,private :: ilo ,ihi ,jlo ,jhi ,klo ,khi
     .                       ,ilom,ihip,jlom,jhip,klom,khip

        integer   ,private ::
     .                i,j,k,igx,igy,igz
     .               ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .               ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp
cc        real(8),private :: car0(3),carp(3),carm(3),jac,ijac
cc     .               ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
cc     .               ,gamma(3,3,3),gamm1(3,3,3),mag,dhp,dhm,dhh

        logical,private    :: map_present

        private :: topol_bc,findGridmetrics

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

        is_map = checkAnalMapDatabase()

        !Check external contributed maps
cc        if (.not.is_map) then
cc          select case(coords)
cc          case('vmc')
cc            is_map = .true.
cc          end select
cc        endif

      end function checkMapDatabase

c     find_dr
c     #####################################################################
      function find_dr(i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz
     .                ,gparams,x)
     .         result(dr)

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

        integer    :: i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz
        real(8)    :: dr(3,3),x(0:nx+1,0:ny+1,0:nz+1,3)

        type(grid_mg_def),pointer :: gparams

c     Local variables

        integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
     .               ,igp,igm,jgp,jgm,kgp,kgm,hel

        real(8)    :: carp(3),carm(3),dh(3),coeff

c     Begin program

        dr = 1d0*reshape((/ 1,0,0
     .                     ,0,1,0
     .                     ,0,0,1 /)
     .                   ,(/ 3,3 /))

        ip=min(i+1,nx+1)
        im=max(i-1,0)
        jp=min(j+1,ny+1)
        jm=max(j-1,0)
        kp=min(k+1,nz+1)
        km=max(k-1,0)

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

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

          dr(1,:) = coeff*(carp-carm)/dh(1)
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

          dr(2,:) = coeff*(carp-carm)/dh(2)
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

          dr(3,:) = coeff*(carp-carm)/dh(3)
cc        endif

      end function find_dr

c     coarsen_map
c     #################################################################
      subroutine coarsen_map(igc,gparams,nxf,nyf,nzf,xf,nxc,nyc,nzc,xc)

c     -----------------------------------------------------------------
c     Coarsen fine map xf(xi) to coarse map xc(xi). This is done
c     dimension by dimension, using first-order extrapolation.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: igc,nxf,nyf,nzf,nxc,nyc,nzc
        real(8)    :: xf(0:nxf+1,0:nyf+1,0:nzf+1,3)
     .               ,xc(0:nxc+1,0:nyc+1,0:nzc+1,3)

        type(grid_mg_def),pointer :: gparams

c     Local variables

        integer    :: ii,igf

c     Begin program

        do ii=1,3
          call crse_map_comp(ii,igc,nxf,nyf,nzf,xf(:,:,:,ii)
     .                             ,nxc,nyc,nzc,xc(:,:,:,ii))
        enddo

      contains

c     crse_map_comp
c     #################################################################
      subroutine crse_map_comp(dim,igc,nxf,nyf,nzf,xf_j
     .                                ,nxc,nyc,nzc,xc_j)

c     -----------------------------------------------------------------
c     Coarsen xf_j to xc_j at grid level igc.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: igc,nxf,nyf,nzf,nxc,nyc,nzc,dim
        real(8)    :: xf_j(0:nxf+1,0:nyf+1,0:nzf+1)
     .               ,xc_j(0:nxc+1,0:nyc+1,0:nzc+1)

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

        x0_j = 0d0

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

              volc   = 0d0
              ac     = 0d0

              do kf = mgz*(kc-1)+1,mgz*kc
                do jf = mgy*(jc-1)+1,mgy*jc
                  do if = mgx*(ic-1)+1,mgx*ic
                    ac   = ac   + xf_j(if,jf,kf)
                    volc = volc + 1d0
                  enddo
                enddo
              enddo

              xc_j(ic,jc,kc) = ac/volc

              i_off = 0
              j_off = 0
              k_off = 0

              ac = 0d0

              select case(dim)
              case (1)
                if (ic == 1  .and.nxc /= 1) i_off =-1
                if (ic == nxc.and.nxc /= 1) i_off = 1

                if (i_off /= 0) then
                  do kf = mgz*(kc-1)+1,mgz*kc
                    do jf = mgy*(jc-1)+1,mgy*jc
                      do if = mgx*(ic-1)+1,mgx*ic
                        ac = ac + xf_j(if+i_off,jf+j_off,kf+k_off)
                      enddo
                    enddo
                  enddo
                endif
              case (2)
                if (jc == 1  .and.nyc /= 1) j_off =-1
                if (jc == nyc.and.nyc /= 1) j_off = 1

                if (j_off /= 0) then
                  do kf = mgz*(kc-1)+1,mgz*kc
                    do jf = mgy*(jc-1)+1,mgy*jc
                      do if = mgx*(ic-1)+1,mgx*ic
                        ac = ac + xf_j(if+i_off,jf+j_off,kf+k_off)
                      enddo
                    enddo
                  enddo
                endif
              case (3)
                if (kc == 1  .and.nzc /= 1) k_off =-1
                if (kc == nzc.and.nzc /= 1) k_off = 1

                if (k_off /= 0) then
                  do kf = mgz*(kc-1)+1,mgz*kc
                    do jf = mgy*(jc-1)+1,mgy*jc
                      do if = mgx*(ic-1)+1,mgx*ic
                        ac = ac + xf_j(if+i_off,jf+j_off,kf+k_off)
                      enddo
                    enddo
                  enddo
                endif
              end select

              x0_j(ic+i_off,jc+j_off,kc+k_off) = ac/volc
            enddo
          enddo
        enddo

c     Fill ghost cells

        call map_bc(dim,igc,nxc,nyc,nzc,xc_j,x0_j)

      end subroutine crse_map_comp

      end subroutine coarsen_map

c     map_bc
c     #################################################################
      subroutine map_bc(dim,igrid,nxl,nyl,nzl,x_j,x0_j,iorder)

c     -----------------------------------------------------------------
c     Impose BCs on map component x_j using boundary of x0_j.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer    :: dim,igrid,nxl,nyl,nzl

        real(8)    :: x_j (0:nxl+1,0:nyl+1,0:nzl+1)
     .               ,x0_j(0:nxl+1,0:nyl+1,0:nzl+1)

        integer   ,optional,intent(IN) :: iorder

c     Local variables

        integer    :: bcnd(6),iibc,order

c     Begin program

        if (PRESENT(iorder)) then
          order = iorder
        else
          order = 1    !Default to first order for predictable results at boundaries
        endif

c     Setup boundary conditions (topological conditions are meaningless for maps!)

        bcnd = EXT

        !Force boundaries along dimension
        select case(dim)
        case(1)
          if (nxl > 1) then
            bcnd(1:2) = EQU
          else
            bcnd(1:2) = PER
          endif
        case(2)
          if (nyl > 1) then
            bcnd(3:4) = EQU
          else
            bcnd(3:4) = PER
          endif
        case(3)
          if (nzl > 1) then
            bcnd(5:6) = EQU
          else
            bcnd(5:6) = PER
          endif
        end select

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

        integer    :: ii,nxf,nyf,nzf,nxc,nyc,nzc,nx,ny,nz
     .               ,ig,jg,kg

        real(8),pointer,dimension(:,:,:,:)   :: xc
        real(8),allocatable,dimension(:,:,:,:,:) :: dr

        logical :: vmec_met=.false.,vmec_mets=.false.

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

c       ANALYTICAL QUANTITIES **********************************

        if ((.not.numerical_grid).and.anal_map) then

          numerical_grid = .false.

          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                gparams%gmetric%grid(igrid)%car  (i,j,k,:)
     .                    = map(i,j,k,igx,igy,igz)
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

                gparams%gmetric%grid(igrid)%dvol(i,j,k)=
     .                           gparams%gmetric%grid(igrid)%jac(i,j,k)
     .                          *gparams%dxh(ig)
     .                          *gparams%dyh(jg)
     .                          *gparams%dzh(kg)
              enddo
            enddo
          enddo

          !Zero force condition on Christoffel symbols (only on finest grid)
          if (igrid == 1) then
            call gammaZeroForce(igx,igy,igz)
          endif

c       NUMERICAL QUANTITIES ***********************************

        else

          numerical_grid = .true.

          allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))

          call findMap

          call findCellMetrics

#if !defined(samrai)
          call findVertexMetrics
#endif
          deallocate(dr)
          nullify(xc)

        endif

c     End program

      contains

c     findMap
c     #################################################################
      subroutine findMap

c     -----------------------------------------------------------------
c     Finds Cartesian map.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: xf

c     Begin program

        xc => gparams%gmetric%grid(igrid)%car

        if (map_present) then

          if (igrid > 1 .and. coarsen) then  !Coarsen map if at coarse grid level
            xf => gparams%gmetric%grid(igrid-1)%car

            nxf=size(xf,1)-2
            nyf=size(xf,2)-2
            nzf=size(xf,3)-2

            nxc=size(xc,1)-2
            nyc=size(xc,2)-2
            nzc=size(xc,3)-2

            call coarsen_map(igrid,gparams,nxf,nyf,nzf,xf
     .                                    ,nxc,nyc,nzc,xc)

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
            call pstop('findMap','Map not available')
          endif

        endif

c     End program

      end subroutine findMap

c     findCellMetrics
c     #################################################################
      subroutine findCellMetrics

c     -----------------------------------------------------------------
c     Finds metrics at cell centers.
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        integer :: he,imin,imax,jmin,jmax,kmin,kmax

        real(8) :: r(3,3),car0(3),carp(3),carm(3),jac,ijac
     .            ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .            ,gamma(3,3,3),mag,dhh,dhm,dhp

c     Begin program

        imin = ilom
        imax = ihip

        jmin = jlom
        jmax = jhip

        kmin = klom
        kmax = khip

        he = 0

c     Local grid size

        nx = ihi-ilo+1
        ny = jhi-jlo+1
        nz = khi-klo+1

c     Evaluate dr=dx/dxi vectors

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              dr(i,j,k,:,:)
     .             = find_dr(i,j,k,igx,igy,igz,he,he,he,nx,ny,nz
     .                      ,gparams,xc)
            enddo
          enddo
        enddo

c     Enforce topological constraints on dr (including SP)

        do i=1,3
          do j=1,3
            call topol_bc(igrid,nx,ny,nz,dr(:,:,:,i,j))
            if (bcSP()) then
              igxbc = igrid
              igybc = igrid
              igzbc = igrid

              nxbc = nx
              nybc = ny
              nzbc = nz
              !dr behaves like a covariant vector
              call vectorSingularBC(i,dr(:,:,:,i,j),.true.)
            endif
          enddo
        enddo

        !Evaluate grid quantities
        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cc              if (.not.vmec_mets) then
                jac = triple_product(r(1,:),r(2,:),r(3,:))
                gparams%gmetric%grid(igrid)%jac(i,j,k) = jac
cc              else
cc                jac = gparams%gmetric%grid(igrid)%jac(i,j,k)
cc              endif

              if (jac <= 0d0 .and. (.not.isSP2(i+1,igrid,ibc=1))) then
                write (*,*) 'Grid level',igrid
                write (*,*) 'Grid position=',i,j,k
                write (*,*) 'dr/dxi_1=',r(1,:)
                write (*,*) 'dr/dxi_2=',r(2,:)
                write (*,*) 'dr/dxi_3=',r(3,:)
                write (*,*) 'Jacobian=',jac
                call pstop('findCellMetrics','Inverted jacobian')
              endif

              ijac = 1d0/jac

              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
              cnv(:,:) = r(:,:)*ijac

              gparams%gmetric%grid(igrid)%cnv(i,j,k,:,:) = cnv

              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
              cov(1,:) = cross_product(r(2,:),r(3,:))
              cov(2,:) = cross_product(r(3,:),r(1,:))
              cov(3,:) = cross_product(r(1,:),r(2,:))

              cov = cov*ijac

              gparams%gmetric%grid(igrid)%cov(i,j,k,:,:) = cov

              !Metric tensors <<<<<<<<<<<<<<<<<<<<<<<<
cc              if (.not.vmec_mets) then
                do m=1,3
                  do l=m,3
                    gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
                    gsub(m,l) = gsub(l,m) !Symmetry
                    gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
                    gsup(m,l) = gsup(l,m) !Symmetry
                  enddo
                enddo
                gparams%gmetric%grid(igrid)%gsub(i,j,k,:,:) = gsub
                gparams%gmetric%grid(igrid)%gsup(i,j,k,:,:) = gsup
cc              endif

              !Cell volume <<<<<<<<<<<<<<<<<<
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              gparams%gmetric%grid(igrid)%dvol(i,j,k)=
     .                           gparams%gmetric%grid(igrid)%jac(i,j,k)
     .                          *gparams%dxh(ig)
     .                          *gparams%dyh(jg)
     .                          *gparams%dzh(kg)

              !Christoffel symbols <<<<<<<<<<<<<<<<<<
              if (i==imin) then
                i0=ilo
              elseif (i==imax) then
                i0=ihi
              else
                i0=i
              endif

              if (j==jmin) then
                j0=jlo
              elseif (j==jmax) then
                j0=jhi
              else
                j0=j
              endif

              if (k==kmin) then
                k0=klo
              elseif (k==kmax) then
                k0=khi
              else
                k0=k
              endif

              ip=i0+1
              im=i0-1
              jp=j0+1
              jm=j0-1
              kp=k0+1
              km=k0-1

              do l=1,3
                do m=1,3
                  do n=m,3
                    if (m == n) then
                      select case(m)
                      case(1)
                        call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
                        call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
                        call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
                        car0 = xc(i0,j0,k0,:)
                        carp = xc(ip,j0,k0,:)
                        carm = xc(im,j0,k0,:)
                        dhp = (gparams%xx(igp)
     .                        -gparams%xx(ig0))
                        dhm = (gparams%xx(ig0)
     .                        -gparams%xx(igm))
                      case(2)
                        call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                        call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                        call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                        car0 = xc(i0,j0,k0,:)
                        carp = xc(i0,jp,k0,:)
                        carm = xc(i0,jm,k0,:)
                        dhp = (gparams%yy(jgp)
     .                        -gparams%yy(jg0))
                        dhm = (gparams%yy(jg0)
     .                        -gparams%yy(jgm))
                      case(3)
                        call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                        call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                        call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                        car0 = xc(i0,j0,k0,:)
                        carp = xc(i0,j0,kp,:)
                        carm = xc(i0,j0,km,:)
                        dhp = (gparams%zz(kgp)
     .                        -gparams%zz(kg0))
                        dhm = (gparams%zz(kg0)
     .                        -gparams%zz(kgm))
                      end select
                      carp = (carp-car0)/dhp
                      carm = (car0-carm)/dhm
                      dhh  = 0.5*(dhp+dhm)
                    else
                      select case(m)
                      case(1)
                       call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
                       call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
                       call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
                       carp= 0.5*(dr(ip,j0,k0,n,:)+dr(i0,j0,k0,n,:))
                       carm= 0.5*(dr(im,j0,k0,n,:)+dr(i0,j0,k0,n,:))
                       dhh = 0.5*(gparams%xx(igp)
     .                           -gparams%xx(igm))
                      case(2)
                       call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                       call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                       call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                       carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
                       carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
                       dhh = 0.5*(gparams%yy(jgp)
     .                           -gparams%yy(jgm))
                      case(3)
                       call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                       call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                       call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                       carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
                       carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
                       dhh = 0.5*(gparams%zz(kgp)
     .                           -gparams%zz(kgm))
                      end select
                    endif
                    vec = (carp-carm)/dhh

                    gamma(l,m,n) = dot_product(vec,cov(l,:))
                    gamma(l,n,m) = gamma(l,m,n) !Symmetry
                  enddo
                enddo
              enddo

              gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma


            enddo
          enddo
        enddo

        !Zero-force condition on Christoffle symbols 
        if (igrid == 1) then
          call gammaZeroForce(igx,igy,igz)
        endif

        !Enforce topological constraints on Christoffel symbols
        do i=1,3
          do j=1,3
            do k=1,3
              call topol_bc(igrid,nx,ny,nz
     .           ,gparams%gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
            enddo
          enddo
        enddo

c     End program

      end subroutine findCellMetrics

#if !defined(samrai)
c     findVertexMetrics
c     #################################################################
      subroutine findVertexMetrics

c     -----------------------------------------------------------------
c     Finds metrics at vertices (i+1/2,j+1/2,k+1/2) (stored at (i,j,k))
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

c     Local variables

        integer :: he,imin,imax,jmin,jmax,kmin,kmax

        real(8) :: r(3,3),cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),jac,ijac

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

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              dr(i,j,k,:,:)
     .             = find_dr(i,j,k,igx,igy,igz,he,he,he,nx,ny,nz
     .                      ,gparams,xc)
            enddo
          enddo
        enddo

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
cc              call vectorSingularBC(i,dr(:,:,:,i,j),.true.)
cc            endif
cc          enddo
cc        enddo

        !Evaluate grid quantities
        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              r = dr(i,j,k,:,:)

              !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              jac = triple_product(r(1,:),r(2,:),r(3,:))
              gparams%gmetric%grid(igrid)%jac_v(i,j,k) = jac

              if (jac <= 0d0 .and. (.not.isSP2(i+1,igrid,ibc=1))) then
                write (*,*) 'Negative Jacobian in half mesh'
                write (*,*) 'Grid level',igrid
                write (*,*) 'Grid position=',i,j,k
                write (*,*) 'dr/dxi_1=',dr(i,j,k,1,:)
                write (*,*) 'dr/dxi_2=',dr(i,j,k,2,:)
                write (*,*) 'dr/dxi_3=',dr(i,j,k,3,:)
                write (*,*) 'Jacobian=',jac
                call pstop('findVertexMetrics','Inverted jacobian')
              elseif (isSP2(i+1,igrid,ibc=1)) then
                ijac = 0d0
              else
                ijac = 1d0/jac
              endif

              !Contravariant vectors <<<<<<<<<<<<<<<<<<<
              cnv(:,:) = r(:,:)*ijac

              gparams%gmetric%grid(igrid)%cnv_v(i,j,k,:,:) = cnv

              !Covariant vectors <<<<<<<<<<<<<<<<<<<<<<
              cov(1,:) = cross_product(r(2,:),r(3,:))
              cov(2,:) = cross_product(r(3,:),r(1,:))
              cov(3,:) = cross_product(r(1,:),r(2,:))

              cov = cov*ijac

              gparams%gmetric%grid(igrid)%cov_v(i,j,k,:,:) = cov

              !Metric tensors <<<<<<<<<<<<<<<<<<<<<<<<
cc              if (.not.vmec_mets) then
                do m=1,3
                  do l=m,3
                    gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
                    gsub(m,l) = gsub(l,m) !Symmetry
                    gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
                    gsup(m,l) = gsup(l,m) !Symmetry
                  enddo
                enddo
                gparams%gmetric%grid(igrid)%gsub_v(i,j,k,:,:) = gsub
                gparams%gmetric%grid(igrid)%gsup_v(i,j,k,:,:) = gsup
cc              endif

            enddo
          enddo
        enddo

c     End program

      end subroutine findVertexMetrics
#endif
      
c     gammaZeroForce
c     #################################################################
      subroutine gammaZeroForce(igx,igy,igz)

c     -----------------------------------------------------------------
c     Postprocess Christoffel symbol to satisfy zero-force condition
c     to machine accuracy
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer    :: igx,igy,igz

c     Local variables

        integer    :: i,j,k,ii,ll,kk,mm,jj,ig,jg,kg,igrid
        real(8)    :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),table(3,3,3),gsub(3,3)
     .               ,gsup(3,3),gsupp(3,3),gsupm(3,3),gamma(3,3,3)
        logical    :: alt_eom

c     Begin program

        igrid = igx

        if (check_grid) then
          write (funit,*)
          write (funit,*) '********************************'
          write (funit,*) 'Checking zero-force condition...'
          write (funit,*) '********************************'
          write (funit,*)
        endif

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

c           START CHECKS

            if (check_grid)
     .         call zeroForceCheck(i,j,k,igrid,gamma
     $                 ,gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:))

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

      end subroutine gammaZeroForce

c     zeroForceCheck
c     #################################################################
      subroutine zeroForceCheck(i,j,k,igrid,gold,gnew)

c     -----------------------------------------------------------------
c     Check zero-force condition of Christoffle symbols
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer    :: i,j,k,igrid
        real(8)    :: gold(3,3,3),gnew(3,3,3)
        
c     Local variables

        integer    :: ii,ll,kk,mm,jj,ig,jg,kg
        real(8)    :: summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),gsup(3,3),gsub(3,3)
     $               ,gsupp(3,3),gsupm(3,3)
        logical    :: alt_eom

c     Begin program

        write(funit,*)
        write(funit,*) 'Grid:',igrid,'  Grid node:',i,j,k

        call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

        gsup  = gparams%gmetric%grid(igrid)%gsup(i,j,k,:,:)
        gsub  = gparams%gmetric%grid(igrid)%gsub(i,j,k,:,:)

        hlf = 0.5

c      Check difference between correction and prediction

        hess = gnew - gold

        summ2= sqrt(sum(gold*gold))
        if (summ2 /= 0d0) then
          summ = sqrt(sum(hess*hess))/summ2
        else
          summ = 0d0
        endif

        write(funit,*) '  Correction size:',summ

c      Check cancellation property g^lk(d_il G^j_kj-g_lm g^ij G^m_jk) = 0

        do ii=1,3 !cycle through Chistoffel symbols
          do ll=1,3
            do kk=1,3
              hess(ii,ll,kk) =
     .               delta(ii,ll)*(gnew(1,kk,1)
     .                            +gnew(2,kk,2)
     .                            +gnew(3,kk,3))

              summ=0d0
              do mm=1,3
                do jj=1,3
                  summ = summ + gsub(ll,mm)*gsup(ii,jj)*gnew(mm,jj,kk)
                enddo
              enddo

              hess(ii,ll,kk) = hess(ii,ll,kk) - summ
            enddo
          enddo

          summ = 0d0
          do ll=1,3
            do kk=1,3
              summ = summ + gsup(ll,kk)*hess(ii,ll,kk)
            enddo
          enddo

          write(funit,'(a,i1,a,e10.3)')
     .       '  Cancel. of corrector term for Gamma^',ii,':',summ
        enddo

c      Check zero-force property: g^lk G^i_lk + nabla_m(g^mi) = 0

        do ii=1,3 !cycle through Chistoffel symbols

          alt_eom = (ii == 2 .and. alt__eom())

          if (alt_eom) then
            jac  = gparams%gmetric%grid(igrid)%jac(i,j,k)
            const= 1d0
          else
            jac  = 1d0
            const= 0d0
          endif

          summ = 0d0

cc          dh = gparams%dxh(ig)
cc          if (alt_eom) then
cc            jacp = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i+1,j,k))
cc            jacm = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i-1,j,k))
cc          else
cc            jacp=1d0
cc            jacm=1d0
cc          endif
cc          gsupp = jacp
cc     $         *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
cc          gsupm = jacm
cc     $         *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i-1,j,k,:,:))
          dh = 2*gparams%dxh(ig)
          if (alt_eom) then
            jacp = gparams%gmetric%grid(igrid)%jac(i+1,j,k)
            jacm = gparams%gmetric%grid(igrid)%jac(i-1,j,k)
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp*gparams%gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
          gsupm = jacm*gparams%gmetric%grid(igrid)%gsup(i-1,j,k,:,:)

          summ = summ +(gsupp(1,ii)-gsupm(1,ii))/dh/jac

cc          dh = gparams%dyh(jg)
cc          if (alt_eom) then
cc            jacp = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j+1,k))
cc            jacm = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j-1,k))
cc          else
cc            jacp=1d0
cc            jacm=1d0
cc          endif
cc          gsupp = jacp
cc     $         *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
cc          gsupm = jacm
cc     $         *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j-1,k,:,:))
          dh = 2*gparams%dyh(jg)
          if (alt_eom) then
            jacp = gparams%gmetric%grid(igrid)%jac(i,j+1,k)
            jacm = gparams%gmetric%grid(igrid)%jac(i,j-1,k)
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp*gparams%gmetric%grid(igrid)%gsup(i,j+1,k,:,:)
          gsupm = jacm*gparams%gmetric%grid(igrid)%gsup(i,j-1,k,:,:)

          summ = summ +(gsupp(2,ii)-gsupm(2,ii))/dh/jac

cc          dh = gparams%dzh(kg)
cc          if (alt_eom) then
cc            jacp = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j,k+1))
cc            jacm = hlf*(jac+gparams%gmetric%grid(igrid)%jac(i,j,k-1))
cc          else
cc            jacp=1d0
cc            jacm=1d0
cc          endif
cc          gsupp = jacp
cc     $         *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
cc          gsupm = jacm
cc     $         *hlf*(gsup+gparams%gmetric%grid(igrid)%gsup(i,j,k-1,:,:))
          dh = 2*gparams%dzh(kg)
          if (alt_eom) then
            jacp = gparams%gmetric%grid(igrid)%jac(i,j,k+1)
            jacm = gparams%gmetric%grid(igrid)%jac(i,j,k-1)
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp*gparams%gmetric%grid(igrid)%gsup(i,j,k+1,:,:)
          gsupm = jacm*gparams%gmetric%grid(igrid)%gsup(i,j,k-1,:,:)

          summ = summ +(gsupp(3,ii)-gsupm(3,ii))/dh/jac

          summ2 = 0d0
          do ll=1,3
            do kk=1,3
              summ2 = summ2 + gsup(ll,kk)*gnew(ii,ll,kk)
     $                - const*gsup(ll,ii)*gnew(kk,ll,kk)
            enddo
          enddo

cc          write(*,*) '  Zero force cancellation for i=',ii,':',summ,summ2,summ+summ2
          write(funit,*) '  Zero force cancellation for i=',ii,':'
     .               ,summ+summ2
        enddo

c     End program

      end subroutine zeroForceCheck

      end subroutine findGridMetrics

c     topol_bc
c     #################################################################
      subroutine topol_bc(igrid,nxl,nyl,nzl,array)

c     -----------------------------------------------------------------
c     !Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer    :: igrid,nxl,nyl,nzl
        real(8)    :: array(0:nxl+1,0:nyl+1,0:nzl+1)

c     Local variables

        integer    :: nxg,nyg,nzg,dim

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

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric(gparams,xcar,igr)

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
        real(8),optional,intent(IN) :: xcar(:,:,:,:)
        integer,optional :: igr

c     Local variables

        integer    :: igrmin,igrmax,igrid,istat

c     Begin program

c     Checks

        map_present = PRESENT(xcar)

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

        if (map_present) then
          gparams%gmetric%grid(igrmin)%car = xcar
        endif

c     Find grid metrics

        do igrid=igrmin,igrmax
          call findGridMetrics(gparams,igrid,.not.PRESENT(igr))
        enddo
        
c     Check grid: test of numerical calculation vs. analytical calculation

        if (check_grid.and.(.not.PRESENT(igr))) then
          call numericalGridCoeffsCheck
          call graphicsCheck
        endif

c     End program

      contains

c     numericalGridCoeffsCheck
c     #################################################################
      subroutine numericalGridCoeffsCheck

c     -----------------------------------------------------------------
c     Test of numerical calculation vs. analytical calculation of
c     grid quantities.
c     -----------------------------------------------------------------

      implicit none

      real(8)    :: mag,mag1,mag2(3,3),mag3(3,3),mag4(3,3,3)

c     Begin program

      if (anal_map) then
        write(funit,*)
        write(funit,*)'***********************************************'
        write(funit,*)'Checking accuracy of numerical grid quantities.'
        write(funit,*)'***********************************************'
        write(funit,*)
      else
        write(funit,*)
        write(funit,*)'***********************************************'
        write(funit,*)'Cannot check numerical grid quantities:'
        write(funit,*)'   analytical mapping unavailable'
        write(funit,*)'***********************************************'
        write(funit,*)
        return
      endif

      do igrid=1,gparams%ngrid

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

        mag1 = 0d0
        mag2 = 0d0
        mag3 = 0d0
        mag4 = 0d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              mag1 = mag1
     .             + (gparams%gmetric%grid(igrid)%jac  (i,j,k)
     .               -jacobian(i,j,k,igx,igy,igz))**2
              mag2 = mag2
     .             + (gparams%gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .               -g_sub   (i,j,k,igx,igy,igz))**2
              mag3 = mag3
     .             + (gparams%gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .               -g_sup   (i,j,k,igx,igy,igz))**2
              mag4 = mag4
     .             + (gparams%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .               -christ_2knd(i,j,k,igx,igy,igz))**2
            enddo
          enddo
        enddo

        write(funit,*) 'Grid level:',igrid
          
        mag = sqrt(    mag1/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1))
        write(funit,*)'   Jacobian tst=',mag
        mag = sum(sqrt(mag2/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1)))/9.
        write(funit,*)'   Gsub tst    =',mag
        mag = sum(sqrt(mag3/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1)))/9.
        write(funit,*)'   Gsup tst    =',mag
        mag = sum(sqrt(mag4/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1)))/27.
        write(funit,*)'   Gamma tst   =',mag
        write(funit,*)

      enddo

      end subroutine numericalGridCoeffsCheck

c     graphicsCheck
c     #################################################################
      subroutine graphicsCheck

c     -----------------------------------------------------------------
c     Test of numerical calculation vs. analytical calculation of
c     grid quantities.
c     -----------------------------------------------------------------

      use grid_debug

      implicit none

      integer    :: nx,ny,nz,i,j,k,l,m,n,igrid,ivar
      real(8)    :: tensor(3,3,3),vec(3)
      character(30) :: debugfile,drawfile

      integer   ,parameter :: nvar=50
      character(30) :: descr(nvar)

      logical :: dump_diff

c     Begin program

c     Output

      write(funit,*)
      write(funit,*) '**********************************'
      write(funit,*) 'Dumping grid debugging graphics...'
      write(funit,*) '**********************************'
      write(funit,*)

cc      dump_diff = .false.
      dump_diff = anal_map.and.numerical_grid

      if (dump_diff) then
        write(funit,*)'Analytical mapping available; dumping diffs...'
      else
        if (anal_map) then
          write(funit,*)'Dumping analytical qtys...'
        else
          write(funit,*)
     .              'Analytical map unavailable; dumping full qtys...'
        endif
      endif
      write(funit,*)

c     Dump quantities

      do igrid=1,gparams%ngrid

        nx = gparams%nxv(igrid)
        ny = gparams%nyv(igrid)
        nz = gparams%nzv(igrid)

        do k=0,nz+1

        write(funit,*)
        write(funit,*) '   ******************************'
        write(funit,*) '   Grid level',igrid,'on slice k=',k
        write(funit,*) '   ******************************'
        write(funit,*)

        !Open .bin file
        debugfile = 'grid_gl'//trim(int2char(igrid))
     .                 //'k='//trim(int2char(k))
     .                 //'p' //trim(int2char(my_rank))//'.bin'
        drawfile = 'drawgl'//trim(int2char(igrid))
     .              //'_k='//trim(int2char(k))
     .               //'_p'//trim(int2char(my_rank))//'.in'
        open(unit=110,file=trim(debugfile),form='unformatted'
     .    ,status='replace')

        !Allocate dummy array
        allocate(dbg(0:nx+1,0:ny+1,1,1))

        !Jacobian
        do j=0,ny+1
          do i=0,nx+1
            dbg(i,j,1,1) = gparams%gmetric%grid(igrid)%jac(i,j,k)
            if (dump_diff) dbg(i,j,1,1) = dbg(i,j,1,1)
     .                    -jacobian(i,j,k,igrid,igrid,igrid)
          enddo
        enddo
        write(funit,*) '   Dumping Jacobian'
        ivar = 1
        descr(ivar) = 'jac'
        call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,0,110)

        !Metric tensors
        do l=1,3
          do m=l,3
            do j=0,ny+1
              do i=0,nx+1
                dbg(i,j,1,1) =
     .               gparams%gmetric%grid(igrid)%gsub(i,j,k,l,m)
     .              *gparams%gmetric%grid(igrid)%jac(i,j,k)
                if (dump_diff) then
                  tensor(:,:,1)  = g_sub(i,j,k,igrid,igrid,igrid)
                  dbg(i,j,1,1) = dbg(i,j,1,1)
     .                    -tensor(l,m,1)
     .                    *jacobian(i,j,k,igrid,igrid,igrid)
                endif
              enddo
            enddo
            write(funit,*) '   Dumping Gsub component',l,m
            ivar = ivar+1
            descr(ivar)='gsub('//trim(int2char(l))//','
     .                         //trim(int2char(m))//')'
            call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1,110)
          enddo
        enddo

        do l=1,3
          do m=l,3
            do j=0,ny+1
              do i=0,nx+1
                dbg(i,j,1,1) =
     .               gparams%gmetric%grid(igrid)%gsup(i,j,k,l,m)
     .              *gparams%gmetric%grid(igrid)%jac(i,j,k)
                if (dump_diff) then
                  tensor(:,:,1)  = g_sup(i,j,k,igrid,igrid,igrid)
                  dbg(i,j,1,1) = dbg(i,j,1,1)
     .                        -tensor(l,m,1)
     .                    *jacobian(i,j,k,igrid,igrid,igrid)
                endif
              enddo
            enddo
            write(funit,*) '   Dumping Gsup component',l,m
            ivar = ivar+1
            descr(ivar)='gsup('//trim(int2char(l))//','
     .                         //trim(int2char(m))//')'
            call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1,110)
          enddo
        enddo

        !Cnv vectors
        do l=1,3
          do m=1,3
            do j=0,ny+1
              do i=0,nx+1
                dbg(i,j,1,1) =
     .               gparams%gmetric%grid(igrid)%cnv(i,j,k,l,m)
     .              *gparams%gmetric%grid(igrid)%jac(i,j,k)
                if (dump_diff) then
                  vec = contravariantVector(l,i,j,k,igrid,igrid,igrid)
     .                 *jacobian(i,j,k,igrid,igrid,igrid)
                  dbg(i,j,1,1) = dbg(i,j,1,1)-vec(m)
                endif
              enddo
            enddo
            write(funit,*) '   Dumping Cnv vector',l,' component',m
            ivar = ivar+1
            descr(ivar)='Cnv_'//trim(int2char(l))//'('
     .                        //trim(int2char(m))//')'
            call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1
     .                    ,110)
          enddo
        enddo

        !Christoffle symbols
        do l=1,3
          do m=1,3
            do n=m,3
              do j=0,ny+1
                do i=0,nx+1
                  dbg(i,j,1,1) =
     .                 gparams%gmetric%grid(igrid)%Gamma(i,j,k,l,m,n)
     .                *gparams%gmetric%grid(igrid)%jac(i,j,k)
                  if (dump_diff) then
                    tensor  = christ_2knd(i,j,k,igrid,igrid,igrid)
                    dbg(i,j,1,1) = dbg(i,j,1,1)
     .                    -tensor(l,m,n)
     .                    *jacobian(i,j,k,igrid,igrid,igrid)
                  endif
                enddo
              enddo
              write(funit,*) '   Dumping Gamma component',l,m,n
              ivar = ivar+1
              descr(ivar)='Gamma'//trim(int2char(l))//'('
     .                           //trim(int2char(m))//','
     .                           //trim(int2char(n))//')'
              call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1
     .                    ,110)
            enddo
          enddo
        enddo

        deallocate(dbg)
        close(110)

        !Create draw*.in file
        call createDrawInCfile(ivar,debugfile,'Grid debug'
     .                        ,'t','x','y',descr,'-c -X0 -L57'
     .                        ,drawfile)

        enddo
      enddo

      write(funit,*)
      write(funit,*) 'Done!'
      write(funit,*) 

      end subroutine graphicsCheck

      end subroutine defineGridMetric

      end module grid_metric_init


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

        logical :: check_grid,numerical_grid

c     PRIVATE VARIABLES

        integer   ,private :: ilo ,ihi ,jlo ,jhi ,klo ,khi
     .                       ,ilom,ihip,jlom,jhip,klom,khip

        integer   ,private ::
     .                i,j,k,igx,igy,igz
     .               ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .               ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp
        real(8),private ::
     .                r(3,3),car0(3),carp(3),carm(3),dh(3),jac,ijac
     .               ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .               ,gamma(3,3,3),gamm1(3,3,3),mag,dhp,dhm,dhh

        real(8),allocatable,dimension(:,:,:,:),private :: xcr

        logical,private    :: map_present

        private :: topol_bc,map_bc,findGridmetrics,coarsenGridMetrics

      contains

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
        if (.not.is_map) then
          select case(coords)
          case('vmc')
            is_map = .true.
          end select
        endif

      end function checkMapDatabase
cc
ccc     find_dr
ccc     #####################################################################
cc      function find_dr(i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq,x)
cc     .         result(dr)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine finds dx/d_xi tensor at the point (i,j,k)
ccc     from the position vector x. It can do this for:
ccc        * Cell centers: (hex,hey,hez)=(0,0,0)
ccc        * Faces, e.g. (hex,hey,hez)=(1,0,0) for i+1/2 face.
ccc        * Corners, e.g. (hex,hey,hez)=(1,1,1) for i+1/2,j+1/2,k+1/2 corner.
ccc     on an arbitrary grid level.
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq
cc        real(8)    :: dr(3,3),x(0:nx+1,0:ny+1,0:nz+1,neq)
cc
ccc     Local variables
cc
cc        integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
cc     .               ,igp,igm,jgp,jgm,kgp,kgm,hel
cc
cc        real(8)    :: carp(neq),carm(neq),dh(3),coeff
cc
ccc     Begin program
cc
cc        dr = 1d0*reshape((/ 1,0,0
cc     .                     ,0,1,0
cc     .                     ,0,0,1 /)
cc     .                   ,(/ 3,3 /))
cc
cc        ip=min(i+1,nx+1)
cc        im=max(i-1,0)
cc        jp=min(j+1,ny+1)
cc        jm=max(j-1,0)
cc        kp=min(k+1,nz+1)
cc        km=max(k-1,0)
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc     Exceptions
cc
cc        if (hex == 1) im = i
cc        if (hey == 1) jm = j
cc        if (hez == 1) km = k
cc
cc        !X
cc        igp = ig+(ip-i)
cc        igm = ig-(i-im)
cc
cc        dh(1)= (gparams%xx(igp)-gparams%xx(igm))
cc
cc        carp = x(ip,j,k,:)
cc        carm = x(im,j,k,:)
cc        coeff = 1d0
cc
cc        if (hey == 1) then
cc          carp = carp + x(ip,jp,k,:)
cc          carm = carm + x(im,jp,k,:)
cc          coeff = 0.5
cc        endif
cc
cc        if (hez == 1) then
cc          carp = carp + x(ip,j,kp,:)
cc          carm = carm + x(im,j,kp,:)
cc          coeff = 0.5
cc        endif
cc
cc        if (hey == 1 .and. hez == 1) then
cc          carp = carp + x(ip,jp,kp,:)
cc          carm = carm + x(im,jp,kp,:)
cc          coeff = 0.25
cc        endif
cc
cc        dr(1,1:neq) = coeff*(carp-carm)/dh(1)
cc
cc        !Y
cc        jgp = jg+(jp-j)
cc        jgm = jg-(j-jm)
cc
cc        dh(2)= (gparams%yy(jgp)-gparams%yy(jgm))
cc
cc        carp = x(i,jp,k,:)
cc        carm = x(i,jm,k,:)
cc        coeff = 1d0
cc
cc        if (hex == 1) then
cc          carp = carp + x(ip,jp,k,:)
cc          carm = carm + x(ip,jm,k,:)
cc          coeff = 0.5
cc        endif
cc
cc        if (hez == 1) then
cc          carp = carp + x(i,jp,kp,:)
cc          carm = carm + x(i,jm,kp,:)
cc          coeff = 0.5
cc        endif
cc
cc        if (hex == 1 .and. hez == 1) then
cc          carp = carp + x(ip,jp,kp,:)
cc          carm = carm + x(ip,jm,kp,:)
cc          coeff = 0.25
cc        endif
cc
cc        dr(2,1:neq) = coeff*(carp-carm)/dh(2)
cc
cc        !Z
cc        kgp = kg+(kp-k)
cc        kgm = kg-(k-km)
cc
cc        dh(3)= (gparams%zz(kgp)-gparams%zz(kgm))
cc
cc        carp = x(i,j,kp,:)
cc        carm = x(i,j,km,:)
cc        coeff = 1d0
cc
cc        if (hex == 1) then
cc          carp = carp + x(ip,j,kp,:)
cc          carm = carm + x(ip,j,km,:)
cc          coeff = 0.5
cc        endif
cc
cc        if (hey == 1) then
cc          carp = carp + x(i,jp,kp,:)
cc          carm = carm + x(i,jp,km,:)
cc          coeff = 0.5
cc        endif
cc
cc        if (hex == 1 .and. hey == 1) then
cc          carp = carp + x(ip,jp,kp,:)
cc          carm = carm + x(ip,jp,km,:)
cc          coeff = 0.25
cc        endif
cc
cc        dr(3,1:neq) = coeff*(carp-carm)/dh(3)
cc
cc      end function find_dr

ccc     find_gmetric
ccc     #####################################################################
cc      subroutine find_gmetric(i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq
cc     .                       ,x,jac,vol,ccov,ccnv,ggsup,ggsub,dtensor)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine finds relevant metric quantities at the point (i,j,k)
ccc     from the position vector x. It can do this for:
ccc        * Cell centers: (hex,hey,hez)=(0,0,0)
ccc        * Faces, e.g. (hex,hey,hez)=(1,0,0) for i+1/2 face.
ccc        * Corners, e.g. (hex,hey,hez)=(1,1,1) for i+1/2,j+1/2,k+1/2 corner.
ccc     on an arbitrary grid level.
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq
cc        real(8)    :: x(0:nx+1,0:ny+1,0:nz+1,neq),jac
cc
cc        real(8),optional,intent(OUT) :: ggsup(3,3),ggsub(3,3),vol
cc     .                                 ,ccov(3,3),ccnv(3,3)
cc        real(8),optional,intent(IN)  :: dtensor(3,3)
cc
ccc     Local variables
cc
cc        integer    :: m,l,ig,jg,kg
cc        real(8)    :: r(3,3),cnv(3,3),cov(3,3),x1,x2
cc     .               ,g_sup(3,3),g_sub(3,3),tensor(3,3)
cc        logical    :: wtnsr
cc
ccc     Begin program
cc
cc        !Set tensor
cc        if (PRESENT(dtensor)) then
cc          wtnsr = .true.
cc        else
cc          wtnsr = .false.
cc        endif
cc
cc        if (wtnsr) tensor = dtensor
cc
cc        !Find dx/dxi
cc        r = find_dr(i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq,x)
cc
cc        !Evaluate Jacobian
cc        jac = triple_product(r(1,:),r(2,:),r(3,:))
cc
cc        !Volume
cc        if (PRESENT(vol)) then
cc          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc          vol = jac*gparams%dxh(ig)
cc     .             *gparams%dyh(jg)
cc     .             *gparams%dzh(kg)
cc        endif
cc
cc        !Vectors
cc        if (PRESENT(ggsup).or.PRESENT(ccov)) then
cc          !Covariant vectors
cc          cov(1,:) = cross_product(r(2,:),r(3,:))
cc          cov(2,:) = cross_product(r(3,:),r(1,:))
cc          cov(3,:) = cross_product(r(1,:),r(2,:))
cc
cc          cov = cov/jac
cc
cc          if (PRESENT(ccov)) ccov = cov
cc        endif
cc
cc        if (PRESENT(ggsub).or.PRESENT(ccnv)) then
cc          !Contravariant vectors
cc          cnv(:,:) = r(:,:)/jac
cc
cc          if (PRESENT(ccov)) ccnv = cnv
cc        endif
cc
cc        !Metric tensors
cc        if (PRESENT(ggsup)) then
cc
cc          if (wtnsr) then
cc            do m=1,3
cc              do l=m,3
cc                g_sup(l,m) = jac
cc     .                *dot_product(cov(l,:),matmul(tensor,cov(m,:)))
cc                g_sup(m,l) = g_sup(l,m) !Symmetry
cc              enddo
cc            enddo
cc          else
cc            do m=1,3
cc              do l=m,3
cc                g_sup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
cc                g_sup(m,l) = g_sup(l,m) !Symmetry
cc              enddo
cc            enddo
cc          endif
cc
cc          ggsup = g_sup
cc
cc        endif
cc
cc        if (PRESENT(ggsub)) then
cc
cc          if (wtnsr) then
cc            do m=1,3
cc              do l=m,3
cc                g_sub(l,m) = jac*dot_product(cnv(l,:)
cc     .                                    ,matmul(tensor,cnv(m,:)))
cc                g_sub(m,l) = g_sub(l,m) !Symmetry
cc              enddo
cc            enddo
cc          else
cc            do m=1,3
cc              do l=m,3
cc                g_sub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
cc                g_sub(m,l) = g_sub(l,m) !Symmetry
cc              enddo
cc            enddo
cc          endif
cc
cc          ggsub = g_sub
cc
cc        endif
cc
cc      end subroutine find_gmetric

c     findGridMetrics
c     #################################################################
      subroutine findGridMetrics(gparams,igrid)

c     -----------------------------------------------------------------
c     Find grid metrics at grid level igrid.
c     -----------------------------------------------------------------

        use grid_debug

        implicit none

c     Input variables

        integer    :: igrid
        type(grid_mg_def),pointer :: gparams

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
              enddo
            enddo
          enddo

cc          write (*,*) 'DIAG -- findGridMetrics'
          !Zero force condition on Christoffel symbols (only on finest grid)
          if (igrid == 1) then
            call gammaZeroForce(igx,igy,igz)
          endif

c       NUMERICAL QUANTITIES ***********************************

        else

          !Rename map
          xc => gparams%gmetric%grid(igrid)%car

          numerical_grid = .true.

          allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))

          !Find cartesian map
          if (map_present) then

            if (igrid == 1) then
              nxf=ihi-ilo+1
              nyf=jhi-jlo+1
              nzf=khi-klo+1

              xc = xcr
            else
              nxf=size(xcr,1)
              nyf=size(xcr,2)
              nzf=size(xcr,3)

              nzc=ihi-ilo+1
              nyc=jhi-jlo+1
              nzc=khi-klo+1

              do ii=1,3
                call coarsen_map(igrid,nxf,nyf,nzf,xcr(:,:,:,ii)
     .                                ,nxc,nyc,nzc,xc (:,:,:,ii))
              enddo

              !Prepare next coarsening
              deallocate(xcr)
              allocate(xcr(0:nxc+1,0:nyc+1,0:nzc+1,3))

              xcr = xc
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

#if defined(vmec)
            elseif (coords == 'vmc') then

                vmec_met = .false.   !Whether to use VMEC metric quantities or not
                vmec_mets = vmec_met .and. (igrid == 1)

                nxf=ihi-ilo+1
                nyf=jhi-jlo+1
                nzf=khi-klo+1

                call vmec_map(igrid,nxf,nyf,nzf,xc,vmec_mets) !External

                if (vmec_mets) then
                  call vmec_metrics(igrid,nxf,nyf,nzf         !External
     .                             ,gparams%gmetric%grid(igrid)%jac
     .                             ,gparams%gmetric%grid(igrid)%gsub
     .                             ,gparams%gmetric%grid(igrid)%gsup)
                endif
#endif
            else
              call pstop('findGridMetrics','Mapping not available')
            endif
          endif

          nx = ihi-ilo+1
          ny = jhi-jlo+1
          nz = khi-klo+1

          !Evaluate dr=dx/dxi vectors
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                dr(i,j,k,:,:)
     .               = find_dr(i,j,k,igx,igy,igz,0,0,0,nx,ny,nz,3,xc)
              enddo
            enddo
          enddo

          !Enforce topological constraints on dr (including SP)
          do i=1,3
            do j=1,3
              call topol_bc(igrid,nx,ny,nz,dr(:,:,:,i,j))
              if (bcSP()) then
cc                write (*,*) 'DIAG -- fingGridMetrics'
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
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                r = dr(i,j,k,:,:)

                !Jacobian <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                if (.not.vmec_mets) then
                  jac = triple_product(r(1,:),r(2,:),r(3,:))
                  gparams%gmetric%grid(igrid)%jac(i,j,k) = jac
                else
                  jac = gparams%gmetric%grid(igrid)%jac(i,j,k)
                endif

                if (jac == 0d0) then
                  write (*,*) i,j,k,r,jac
                  call pstop('findGridMetrics','Zero jacobian')
                  stop
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
                if (.not.vmec_mets) then
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
                endif

                !Christoffel symbols <<<<<<<<<<<<<<<<<<
                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
                dh(1) = gparams%dxh(ig)
                dh(2) = gparams%dyh(jg)
                dh(3) = gparams%dzh(kg)

                if (i==ilom) then
                  i0=ilo
                elseif (i==ihip) then
                  i0=ihi
                else
                  i0=i
                endif

                if (j==jlom) then
                  j0=jlo
                elseif (j==jhip) then
                  j0=jhi
                else
                  j0=j
                endif

                if (k==klom) then
                  k0=klo
                elseif (k==khip) then
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
     .                          -gparams%xx(ig0))
                          dhm = (gparams%xx(ig0)
     .                          -gparams%xx(igm))
                        case(2)
                          call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                          call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                          call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                          car0 = xc(i0,j0,k0,:)
                          carp = xc(i0,jp,k0,:)
                          carm = xc(i0,jm,k0,:)
                          dhp = (gparams%yy(jgp)
     .                          -gparams%yy(jg0))
                          dhm = (gparams%yy(jg0)
     .                          -gparams%yy(jgm))
                        case(3)
                          call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                          call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                          call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                          car0 = xc(i0,j0,k0,:)
                          carp = xc(i0,j0,kp,:)
                          carm = xc(i0,j0,km,:)
                          dhp = (gparams%zz(kgp)
     .                          -gparams%zz(kg0))
                          dhm = (gparams%zz(kg0)
     .                          -gparams%zz(kgm))
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
     .                             -gparams%xx(igm))
                        case(2)
                         call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                         call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                         call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                         carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
                         carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
                         dhh = 0.5*(gparams%yy(jgp)
     .                             -gparams%yy(jgm))
                        case(3)
                         call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                         call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                         call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                         carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
                         carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
                         dhh = 0.5*(gparams%zz(kgp)
     .                             -gparams%zz(kgm))
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

          deallocate(dr)
          nullify(xc)

          !Zero-force condition on Christoffle symbols 
          if (igrid == 1) then
            call gammaZeroForce(igx,igy,igz)
          endif

          !Enforce topological constraints on Christoffel symbols
          do i=1,3
            do j=1,3
              do k=1,3
                call topol_bc(igrid,nx,ny,nz
     .             ,gparams%gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
              enddo
            enddo
          enddo

        endif

c       Store cell volume

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              gparams%gmetric%grid(igrid)%dvol(i,j,k)=
     .                           gparams%gmetric%grid(igrid)%jac(i,j,k)
     .                          *gparams%dxh(ig)
     .                          *gparams%dyh(jg)
     .                          *gparams%dzh(kg)
cc              gparams%gmetric%grid(igrid)%dvol(i,j,k)=volume(i,j,k,igx,igy,igz)
            enddo
          enddo
        enddo

c     End program

      contains


c     find_dr
c     #####################################################################
      function find_dr(i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq,x)
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

        integer    :: i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq
        real(8)    :: dr(3,3),x(0:nx+1,0:ny+1,0:nz+1,neq)

c     Local variables

        integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
     .               ,igp,igm,jgp,jgm,kgp,kgm,hel

        real(8)    :: carp(neq),carm(neq),dh(3),coeff

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

        dr(1,1:neq) = coeff*(carp-carm)/dh(1)

        !Y
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

        dr(2,1:neq) = coeff*(carp-carm)/dh(2)

        !Z
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

        dr(3,1:neq) = coeff*(carp-carm)/dh(3)

      end function find_dr

c     coarsen_map
c     #################################################################
      subroutine coarsen_map(igc,nxf,nyf,nzf,arrayf,nxc,nyc,nzc,arrayc)

c     -----------------------------------------------------------------
c     Coarsen arrayf to arrayc at grid level igc.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: igc,nxf,nyf,nzf,nxc,nyc,nzc
        real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1)
     .               ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1)
cc        logical    :: volf

c     Local variables

        integer    :: igf,ic,jc,kc,if,jf,kf,mgx,mgy,mgz

        real(8)    :: volc,ac

c     Begin program

        igf = igc - 1

        mgx = gparams%mg_ratio_x(igf)
        mgy = gparams%mg_ratio_y(igf)
        mgz = gparams%mg_ratio_z(igf)

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

              volc   = 0d0
              ac     = 0d0

              do kf = mgz*(kc-1)+1,mgz*kc
                do jf = mgy*(jc-1)+1,mgy*jc
                  do if = mgx*(ic-1)+1,mgx*ic

cc                    if (volf) then
cc                      ac   = ac   + gparams%gmetric%grid(igf)%dvol(if,jf,kf)
cc     .                             *arrayf(if,jf,kf)
cc                      volc = volc + gparams%gmetric%grid(igf)%dvol(if,jf,kf)
cc                    else
                      ac   = ac   + arrayf(if,jf,kf)
                      volc = volc + 1d0
cc                    endif

                  enddo
                enddo
              enddo

              arrayc(ic,jc,kc) = ac/volc
            enddo
          enddo
        enddo

c     Fill ghost cells by 2nd-order extrapolation

        call map_bc(igc,nxc,nyc,nzc,arrayc,IBC=EXT,iorder=1)

      end subroutine coarsen_map

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
          write (*,*)
          write (*,*) '********************************'
          write (*,*) 'Checking zero-force condition...'
          write (*,*) '********************************'
          write (*,*)
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

            if (check_grid) call zeroForceCheck(i,j,k,igrid,gamma
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

        write (*,*)
        write (*,*) 'Grid:',igrid,'  Grid node:',i,j,k

        call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

        gsup  = gparams%gmetric%grid(igrid)%gsup(i,j,k,:,:)
        gsub  = gparams%gmetric%grid(igrid)%gsub(i,j,k,:,:)

        hlf = 0.5

c      Check difference between correction and prediction

        hess = gnew - gold

        summ2= sqrt(sum(gold*gold))
        summ = sqrt(sum(hess*hess))/summ2

        write (*,*) '  Correction size:',summ

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

          write(*,*)'  Cancel. of corrector term for Gamma^',ii,':',summ
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
          write (*,*) '  Zero force cancellation for i=',ii,':'
     .               ,summ+summ2
        enddo

c     End program

      end subroutine zeroForceCheck

      end subroutine findGridMetrics

c     coarsenGridMetrics
c     #################################################################
      subroutine coarsenGridMetrics(gparams,igc)

c     -----------------------------------------------------------------
c     Coasens grid metrics to grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: igc
cc        type(MG_grid_metrics),pointer :: g_metric
        type(grid_mg_def),pointer :: gparams

c     Local variables

        integer    :: igf,ic,jc,kc,if,jf,kf,icg,jcg,kcg
     .               ,nxc,nyc,nzc,nxf,nyf,nzf,l,m

c     Begin program

        call findGridMetrics(gparams,igc)
        return

        !From analytical mapping
        if (.not.numerical_grid) then

          call findGridMetrics(gparams,igc)

        !From numerical mapping
        else

          igf = igc - 1

          nxf = gparams%nxv(igf)
          nyf = gparams%nyv(igf)
          nzf = gparams%nzv(igf)

          nxc = gparams%nxv(igc)
          nyc = gparams%nyv(igc)
          nzc = gparams%nzv(igc)

c         Coarsen volumes and jacobians (needed for coarsening all else)

          !Coarsen cell volumes (with ghost cells)
          call coarsen(igc,nxf,nyf,nzf,gparams%gmetric%grid(igf)%dvol
     .                    ,nxc,nyc,nzc,gparams%gmetric%grid(igc)%dvol
     .                ,.false.)

          !Calculate jacobian
cc          do kc = 0,nzc+1
cc            do jc = 0,nyc+1
cc              do ic = 0,nxc+1
          do kc = 1,nzc
            do jc = 1,nyc
              do ic = 1,nxc
                call getMGmap(ic,jc,kc,igc,igc,igc,icg,jcg,kcg)

                gparams%gmetric%grid(igc)%jac(ic,jc,kc) =
     .                      gparams%gmetric%grid(igc)%dvol(ic,jc,kc)
     .                                    /(gparams%dxh(icg)
     .                                     *gparams%dyh(jcg)
     .                                     *gparams%dzh(kcg))
cc                write (*,*) 'Grid=',igc,'i,j,k=',ic,jc,kc
cc     .               ,' jac=',gparams%gmetric%grid(igc)%jac (ic,jc,kc)
cc     .               ,' vol=',gparams%gmetric%grid(igc)%dvol(ic,jc,kc)
cc     .               ,' dxh=',gparams%dxh(icg)
              enddo
            enddo
          enddo

          call map_bc(igc,nxc,nyc,nzc,gparams%gmetric%grid(igc)%jac
     .               ,IBC=EXT,iorder=1)

c         Coarsen everything else (with ghost cells)

          !Gsup
          do l=1,3
            do m=l,3
              call coarsen(igc
     .           ,nxf,nyf,nzf,gparams%gmetric%grid(igf)%gsup(:,:,:,l,m)
     .           ,nxc,nyc,nzc,gparams%gmetric%grid(igc)%gsup(:,:,:,l,m)
     .           ,.true.)
              gparams%gmetric%grid(igc)%gsup(:,:,:,m,l)
     .             = gparams%gmetric%grid(igc)%gsup(:,:,:,l,m) !Symmetry
            enddo
          enddo

          !Gsub
          do l=1,3
            do m=l,3
              call coarsen(igc
     .           ,nxf,nyf,nzf,gparams%gmetric%grid(igf)%gsub(:,:,:,l,m)
     .           ,nxc,nyc,nzc,gparams%gmetric%grid(igc)%gsub(:,:,:,l,m)
     .           ,.true.)
              gparams%gmetric%grid(igc)%gsub(:,:,:,m,l)
     .             = gparams%gmetric%grid(igc)%gsub(:,:,:,l,m) !Symmetry
            enddo
          enddo

          !Gamma
          do m=1,3
            do n=m,3
              do l=1,3
                call coarsen(igc
     .         ,nxf,nyf,nzf,gparams%gmetric%grid(igf)%Gamma(:,:,:,l,m,n)
     .         ,nxc,nyc,nzc,gparams%gmetric%grid(igc)%Gamma(:,:,:,l,m,n)
     .         ,.true.)
              enddo
              gparams%gmetric%grid(igc)%Gamma(:,:,:,:,n,m)
     .             = gparams%gmetric%grid(igc)%Gamma(:,:,:,:,m,n) !Symmetry
            enddo
          enddo

          !Cov
          do l=1,3
            do m=1,3
              call coarsen(igc
     .             ,nxf,nyf,nzf,gparams%gmetric%grid(igf)%cov(:,:,:,l,m)
     .             ,nxc,nyc,nzc,gparams%gmetric%grid(igc)%cov(:,:,:,l,m)
     .             ,.true.)
            enddo
          enddo

          !Cnv
          do l=1,3
            do m=1,3
              call coarsen(igc
     .             ,nxf,nyf,nzf,gparams%gmetric%grid(igf)%cnv(:,:,:,l,m)
     .             ,nxc,nyc,nzc,gparams%gmetric%grid(igc)%cnv(:,:,:,l,m)
     .             ,.true.)
            enddo
          enddo

        endif

c     End program

      contains

c     coarsen
c     #################################################################
      subroutine coarsen(igc,nxf,nyf,nzf,arrayf,nxc,nyc,nzc,arrayc,volf)

c     -----------------------------------------------------------------
c     Coarsen arrayf to arrayc at grid level igc.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: igc,nxf,nyf,nzf,nxc,nyc,nzc
        real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1)
     .               ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1)
        logical    :: volf

c     Local variables

        integer    :: igf,ic,jc,kc,if,jf,kf,mgx,mgy,mgz

        real(8)    :: volc,ac

c     Begin program

        igf = igc - 1

        mgx = gparams%mg_ratio_x(igf)
        mgy = gparams%mg_ratio_y(igf)
        mgz = gparams%mg_ratio_z(igf)

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

              volc   = 0d0
              ac     = 0d0

              do kf = mgz*(kc-1)+1,mgz*kc
                do jf = mgy*(jc-1)+1,mgy*jc
                  do if = mgx*(ic-1)+1,mgx*ic

                    if (volf) then
                      ac   = ac
     .                     + gparams%gmetric%grid(igf)%dvol(if,jf,kf)
     .                      *arrayf(if,jf,kf)
                      volc = volc
     .                     + gparams%gmetric%grid(igf)%dvol(if,jf,kf)
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

        call map_bc(igc,nxc,nyc,nzc,arrayc,IBC=EXT,iorder=1)

      end subroutine coarsen

      end subroutine coarsenGridMetrics

c     topol_bc
c     #################################################################
      subroutine topol_bc(igrid,nxl,nyl,nzl,array)

c     -----------------------------------------------------------------
c     !Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        real(8)    :: array(0:nxl+1,0:nyl+1,0:nzl+1)

        integer    :: igrid

c     Local variables

        integer    :: nxl,nyl,nzl,nxg,nyg,nzg,dim

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
          call fillPetscGhostCells(array,igrid)
cc        endif
#else
          write(*,*) 'WARNING: metrics in SAMRAI are not being'
     .              ,' communicated yet.'
#endif
#endif

c     End program

      end subroutine topol_bc

c     map_bc
c     #################################################################
      subroutine map_bc(igrid,nxl,nyl,nzl,array,IBC,iorder)

c     -----------------------------------------------------------------
c     Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer    :: igrid,nxl,nyl,nzl

        real(8)    :: array(0:nxl+1,0:nyl+1,0:nzl+1)

        integer   ,optional,intent(IN) :: IBC,iorder

c     Local variables

        integer    :: bcnd(6),iibc,order
        real(8)    :: array0(0:nxl+1,0:nyl+1,0:nzl+1)

c     Degug

        integer    :: i,j,k

c     External

c     Begin program

        if (PRESENT(IBC)) then
          iibc = IBC
        else
          iibc = EQU
        endif

        if (PRESENT(iorder)) then
          order = iorder
        else
          order = 1
        endif

        array0 = array

c     Set boundary conditions to IBC

        bcnd = bcond

cc        if (iibc == EQU) then
          where(bcond /= PER) 
            bcnd = iibc
          elsewhere
            bcnd = DEF
          endwhere
cc        else
cc          where(bcond == DEF) bcnd = iibc
cc        endif

c     Impose BCs

        call setMGBC(0,1,nxl,nyl,nzl,igrid,array,bcnd
     .              ,arr0=array0,iorder=order)

c diag ******
cc        write (*,*) nxl,nyl,nzl
cc        write (*,*) 'bcnd',bcnd
cc        do j=1,nyl
cc          write (*,*) 'before, j=',j,array0(:,j,1)
cc          write (*,*) 'after , j=',j,array (:,j,1)
cc          write (*,*)
cc        enddo
cc        write (*,*) '==============================='
c diag ******

c     End program

      end subroutine map_bc

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric(gparams,xcar)

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
cc        type(MG_grid_metrics),pointer :: g_metric
        real(8),optional,intent(IN) :: xcar(:,:,:,:)

c     Local variables

        integer    :: igrmin,igrmax,igrid,istat

c     Begin program

c     Checks

        map_present = PRESENT(xcar)

        anal_map = checkAnalMapDatabase().and.(.not.map_present)

c     Allocate arrays

        call allocateGridMetric(gparams%gmetric,gparams)

        if (map_present) then
          allocate(xcr(1:gparams%nxv(1)
     .                ,1:gparams%nyv(1)
     .                ,1:gparams%nzv(1),3))
          xcr = xcar
        endif

c     Define grid levels

cc        if (PRESENT(igr)) then
cc          igrmin = igr
cc          igrmax = igr
cc        else
          igrmin = 1
          igrmax = gparams%ngrid
cc        endif
cc
ccc     Find grid metrics
cc
cc        do igrid=igrmin,igrmax
cc          call findGridMetrics(igrid)
cc        enddo

c     Find grid metrics in fine grid

        call findGridMetrics(gparams,igrmin)

c     Coarsen grid metrics (for MG)

        do igrid=igrmin+1,igrmax
          call coarsenGridMetrics(gparams,igrid)
        enddo
        
c     Check grid: test of numerical calculation vs. analytical calculation

cc        if (check_grid .and. numerical_grid) then
        if (check_grid) then
          call numericalGridCoeffsCheck
          call graphicsCheck
        endif

c     End program

        deallocate(xcr,stat=istat)

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
        write (*,*)
        write (*,*) '*************************************************'
        write (*,*) 'Checking accuracy of numerical grid quantities...'
        write (*,*) '*************************************************'
        write (*,*)
      else
        write (*,*)
        write (*,*) '*************************************************'
        write (*,*) 'Cannot check numerical grid quantities:'
        write (*,*) '   analytical mapping unavailable'
        write (*,*) '*************************************************'
        write (*,*)
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

        write (*,*) 'Grid level:',igrid
          
        mag = sqrt(    mag1 /(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1))
        write (*,*)'   Jacobian tst=',mag
        mag = sum(sqrt(mag2)/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1))/9.
        write (*,*)'   Gsub tst    =',mag
        mag = sum(sqrt(mag3)/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1))/9.
        write (*,*)'   Gsup tst    =',mag
        mag = sum(sqrt(mag4)/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1))/27.
        write (*,*)'   Gamma tst   =',mag
        write (*,*)

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

      write (*,*)
      write (*,*) '**********************************'
      write (*,*) 'Dumping grid debugging graphics...'
      write (*,*) '**********************************'
      write (*,*)

cc      dump_diff = anal_map.and.(.not.numerical_grid)
      dump_diff = anal_map.and.numerical_grid

      if (dump_diff) then
        write(*,*)'Analytical mapping available; dumping diffs...'
      else
        if (anal_map) then
          write(*,*)'Dumping analytical qtys...'
        else
          write(*,*)'Analytical map unavailable; dumping full qtys...'
        endif
      endif
      write (*,*)

c     Dump quantities

      do igrid=1,gparams%ngrid

        nx = gparams%nxv(igrid)
        ny = gparams%nyv(igrid)
        nz = gparams%nzv(igrid)

        do k=0,nz+1

        write (*,*)
        write (*,*) '   ******************************'
        write (*,*) '   Grid level',igrid,'on slice k=',k
        write (*,*) '   ******************************'
        write (*,*)

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
        write (*,*) '   Dumping Jacobian'
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
            write (*,*) '   Dumping Gsub component',l,m
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
            write (*,*) '   Dumping Gsup component',l,m
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
            write (*,*) '   Dumping Cnv vector',l,' component',m
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
              write (*,*) '   Dumping Gamma component',l,m,n
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

      write (*,*)
      write (*,*) 'Done!'
      write (*,*) 

      end subroutine graphicsCheck

      end subroutine defineGridMetric

      end module grid_metric_init

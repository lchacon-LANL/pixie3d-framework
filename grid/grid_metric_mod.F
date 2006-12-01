
c module grid_metric
c #####################################################################
      module grid_metric

c ---------------------------------------------------------------------
c     This module allocates and fills the grid metric structure gmetric.
c     It contains the following routines: 
c        * allocateGridMetric
c        * deallocateGridMetric
c        * defineGridMetric
c           -> map
c           -> topol_bc
c           -> gammaZeroForce
c           -> zeroForceCheck
c        * cross_product
c        * triple_product
c        * delta
c ---------------------------------------------------------------------

        use grid_definition

        use local_mpi

        use bc_def

        implicit none

        type :: grid_metrics
cc          real(8),pointer,dimension(:,:,:,:)     :: car   !Node positions in cartesian space
          real(8),pointer,dimension(:,:,:)       :: dvol  !Local cell volume
          real(8),pointer,dimension(:,:,:)       :: jac   !Jacobian factor at grid cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsub  !Covariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsup  !Contravariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:,:) :: Gamma !Christoffel symbol at cell centers
          real(8),pointer,dimension(:,:,:,:,:)   :: cov   !Covariant vectors
          real(8),pointer,dimension(:,:,:,:,:)   :: cnv   !Contravariant vectors
        end type grid_metrics

        type :: MG_grid_metrics
          type(grid_metrics),pointer,dimension(:) :: grid
        end type MG_grid_metrics

        type(MG_grid_metrics) :: gmetric

        logical :: check_grid

      contains

c     allocateGridMetric
c     #################################################################
      subroutine allocateGridMetric(gmetric)

        implicit none

c     Call variables

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4)      :: igrid,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        if (.not.associated(gmetric%grid)) then
          allocate(gmetric%grid(grid_params%ngrid))
          do igrid=1,grid_params%ngrid
            ilom = 0
            jlom = 0
            klom = 0
            ihip = grid_params%nxv(igrid)+1
            jhip = grid_params%nyv(igrid)+1
            khip = grid_params%nzv(igrid)+1
cc            allocate(gmetric%grid(igrid)
cc     $           %car  (ilom:ihip,jlom:jhip,klom:khip,3))
            allocate(gmetric%grid(igrid)
     $           %dvol (ilom:ihip,jlom:jhip,klom:khip))
            allocate(gmetric%grid(igrid)
     $           %jac  (ilom:ihip,jlom:jhip,klom:khip))
            allocate(gmetric%grid(igrid)
     $           %gsub (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %gsup (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %cov  (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %cnv  (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %Gamma(ilom:ihip,jlom:jhip,klom:khip,3,3,3))

          enddo
        endif

c     End program

      end subroutine allocateGridMetric

c     deallocateGridMetric
c     #################################################################
      subroutine deallocateGridMetric(gmetric)

        implicit none

c     Call variables

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4)      :: igrid,nxp,nyp,nzp

c     Begin program

        if (associated(gmetric%grid)) then
          do igrid=1,grid_params%ngrid
cc            deallocate(gmetric%grid(igrid)%car  )
            deallocate(gmetric%grid(igrid)%dvol )
            deallocate(gmetric%grid(igrid)%jac  )
            deallocate(gmetric%grid(igrid)%gsub )
            deallocate(gmetric%grid(igrid)%gsup )
            deallocate(gmetric%grid(igrid)%cov  )
            deallocate(gmetric%grid(igrid)%cnv  )
            deallocate(gmetric%grid(igrid)%Gamma)
          enddo
          deallocate(gmetric%grid)
        endif

c     End program

      end subroutine deallocateGridMetric

c     defineGridMetric
c     #################################################################
cc      subroutine defineGridMetric(gmetric,igr)
      subroutine defineGridMetric(gmetric,xcar)

c     -----------------------------------------------------------------
c     This routine calculates all grid metric quantities required for
c     the curvilinear representation of a set of PDE's: jacobian,
c     metric tensors, covariant and contravariant vectors, Christoffel
c     symbols of the second kind. All quantities are stored in
c     structure gmetric. There are two modes of computation:
c        * Analytical (numerical_grid=.false.)
c        * Numerical  (numerical_grid=.true.)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

cc        integer(4),optional,intent(IN) :: igr

        type(MG_grid_metrics) :: gmetric
        real(8),optional,intent(IN) :: xcar(:,:,:,:)

c     Local variables

        integer(4) :: igrmin,igrmax,igrid,ilo,ihi,jlo,jhi,klo,khi
     $               ,ilom,ihip,jlom,jhip,klom,khip,istat

        integer(4) :: i,j,k,igx,igy,igz
     .               ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .               ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp
        real(8)    :: r(3,3),car0(3),carp(3),carm(3),dh(3),jac,ijac
     .               ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .               ,gamma(3,3,3),gamm1(3,3,3),mag,dhp,dhm,dhh
        real(8),allocatable,dimension(:,:,:,:,:) :: dr
        real(8),allocatable,dimension(:,:,:,:)   :: xcr

        logical    :: map_present

c     Begin program

c     Checks

        map_present = PRESENT(xcar)

        anal_map = checkAnalMapDatabase().and.(.not.map_present)

c     Allocate arrays

        call allocateGridMetric(gmetric)

        if (map_present) then
          allocate(xcr(1:grid_params%nxv(1)
     .                ,1:grid_params%nyv(1)
     .                ,1:grid_params%nzv(1),3))
          xcr = xcar
        endif

c     Define grid levels

cc        if (PRESENT(igr)) then
cc          igrmin = igr
cc          igrmax = igr
cc        else
          igrmin = 1
          igrmax = grid_params%ngrid
cc        endif
cc
ccc     Find grid metrics
cc
cc        do igrid=igrmin,igrmax
cc          call findGridMetrics(igrid)
cc        enddo

c     Find grid metrics in fine grid

        call findGridMetrics(igrmin)

c     Coarsen grid metrics (for MG)

        do igrid=igrmin+1,igrmax
          call coarsenGridMetrics(igrid)
        enddo
        
c     Check grid: test of numerical calculation vs. analytical calculation

        if (check_grid .and. numerical_grid) then
          call numericalGridCoeffsCheck
          call graphicsCheck
        endif

c     End program

        deallocate(xcr,stat=istat)

      contains

c     findGridMetrics
c     #################################################################
      subroutine findGridMetrics(igrid)

c     -----------------------------------------------------------------
c     Find grid metrics at grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: igrid

c     Local variables

        integer(4) :: ii

        real(8),allocatable,dimension(:,:,:,:)   :: xc

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        ilo = 1
        jlo = 1
        klo = 1
        ihi = grid_params%nxv(igrid)
        jhi = grid_params%nyv(igrid)
        khi = grid_params%nzv(igrid)

        ilom = ilo-1
        jlom = jlo-1
        klom = klo-1
        ihip = ihi+1
        jhip = jhi+1
        khip = khi+1

c       ANALYTICAL QUANTITIES **********************************

        if ((.not.numerical_grid).and.anal_map) then

          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                gmetric%grid(igrid)%jac  (i,j,k)
     .                    = jacobian(i,j,k,igx,igy,igz)
                gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                    = g_sub   (i,j,k,igx,igy,igz)
                gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                    = g_sup   (i,j,k,igx,igy,igz)
                gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .                    = christ_2knd(i,j,k,igx,igy,igz)
                do l=1,3
                  gmetric%grid(igrid)%cov(i,j,k,l,:)
     .                    = covariantVector    (l,i,j,k,igx,igy,igz)
                  gmetric%grid(igrid)%cnv(i,j,k,l,:)
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

          numerical_grid = .true.

          allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3)
     .            ,xc(ilom:ihip,jlom:jhip,klom:khip,3))

          !Find cartesian map
          if (map_present) then

            if (igrid == 1) then
              xc = xcr
            else
              do ii=1,3
                call coarsen(igrid,xcr(:,:,:,ii),xc(:,:,:,ii),.false.)
              enddo
              deallocate(xcr)
              allocate(xcr(ilom:ihip,jlom:jhip,klom:khip,3))
              xcr = xc    !For coarsening later
            endif

          else

            if (checkAnalMapDatabase()) then

              do k = klom,khip
                do j = jlom,jhip
                  do i = ilom,ihip
                    xc(i,j,k,:) = map(i,j,k,igx,igy,igz,ig,jg,kg)
                  enddo
                enddo
              enddo

#if defined(vmec)
            elseif (coords == 'vmc') then

              if (igrid == 1) then
                call vmec_map(igrid,ihi-ilo+1,jhi-jlo+1,khi-klo+1,xc) !External
                allocate(xcr(ilom:ihip,jlom:jhip,klom:khip,3))
                xcr = xc   !For coarsening later
              else
                do ii=1,3
                  call coarsen(igrid,xcr(:,:,:,ii),xc(:,:,:,ii),.false.)
                enddo
                deallocate(xcr)
                allocate(xcr(ilom:ihip,jlom:jhip,klom:khip,3))
                xcr = xc   !For coarsening later
              endif
#endif
            else
              call pstop('findGridMetrics'
     .                  ,'Mapping not available')
            endif
          endif

          if (.not.checkAnalMapDatabase()) then
            !Enforce topological constraints on the map xc
            do i=1,3
              call topol_bc(xc(:,:,:,i),igrid)
            enddo
          endif

          do i=1,3
            do k=klom,khip
              do j=jlom,jhip
                write (*,*) 'comp=',i,'  ',j,k,xc(:,j,k,i)
              enddo
            enddo
            write (*,*)
          enddo
          stop

          !Evaluate dr=dx/dxi vectors
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip

                dr(i,j,k,:,:)=find_dr(i,j,k,igx,igy,igz,0,0,0
     .                            ,ihi-ilo+1,jhi-jlo+1,khi-klo+1,3,xc)

cc                ip=min(i+1,ihip)
cc                im=max(i-1,ilom)
cc                jp=min(j+1,jhip)
cc                jm=max(j-1,jlom)
cc                kp=min(k+1,khip)
cc                km=max(k-1,klom)
cc
cc                carp = map(ip,j,k,igx,igy,igz,igp,jg,kg)
cc                carm = map(im,j,k,igx,igy,igz,igm,jg,kg)
cc                dh(1)= (grid_params%xx(igp)-grid_params%xx(igm))
cc                dr(i,j,k,1,:) = (carp-carm)/dh(1)
cc
cc                carp = map(i,jp,k,igx,igy,igz,ig,jgp,kg)
cc                carm = map(i,jm,k,igx,igy,igz,ig,jgm,kg)
cc                dh(2)= (grid_params%yy(jgp)-grid_params%yy(jgm))
cc                dr(i,j,k,2,:) = (carp-carm)/dh(2)
cc
cc                carp = map(i,j,kp,igx,igy,igz,ig,jg,kgp)
cc                carm = map(i,j,km,igx,igy,igz,ig,jg,kgm)
cc                dh(3)= (grid_params%zz(kgp)-grid_params%zz(kgm))
cc                dr(i,j,k,3,:) = (carp-carm)/dh(3)

              enddo
            enddo
          enddo

          !Enforce topological constraints on dr
          do i=1,3
            do j=1,3
              call topol_bc(dr(:,:,:,i,j),igrid)
            enddo
          enddo

          !Evaluate grid quantities
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                r = dr(i,j,k,:,:)

                !Evaluate Jacobian
                jac = triple_product(r(1,:),r(2,:),r(3,:))
                ijac = 1d0/jac

                !Contravariant vectors
                cnv(:,:) = r(:,:)*ijac

                !Covariant vectors
                cov(1,:) = cross_product(r(2,:),r(3,:))
                cov(2,:) = cross_product(r(3,:),r(1,:))
                cov(3,:) = cross_product(r(1,:),r(2,:))

                cov = cov*ijac

                !Metric tensors
                do m=1,3
                  do l=m,3
                    gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
                    gsub(m,l) = gsub(l,m) !Symmetry
                    gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
                    gsup(m,l) = gsup(l,m) !Symmetry
                  enddo
                enddo

                !Grid spacings
                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
                dh(1) = grid_params%dxh(ig)
                dh(2) = grid_params%dyh(jg)
                dh(3) = grid_params%dzh(kg)

                !Christoffel symbols

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
cc                          car0 = map(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
cc                          carp = map(ip,j0,k0,igx,igy,igz,igp,jg,kg)
cc                          carm = map(im,j0,k0,igx,igy,igz,igm,jg,kg)
                          call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
                          call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
                          call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
                          car0 = xc(i0,j0,k0,:)
                          carp = xc(ip,j0,k0,:)
                          carm = xc(im,j0,k0,:)
                          dhp = (grid_params%xx(igp)
     .                          -grid_params%xx(ig0))
                          dhm = (grid_params%xx(ig0)
     .                          -grid_params%xx(igm))
                        case(2)
cc                          car0 = map(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
cc                          carp = map(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
cc                          carm = map(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                          call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                          call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                          call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                          car0 = xc(i0,j0,k0,:)
                          carp = xc(i0,jp,k0,:)
                          carm = xc(i0,jm,k0,:)
                          dhp = (grid_params%yy(jgp)
     .                          -grid_params%yy(jg0))
                          dhm = (grid_params%yy(jg0)
     .                          -grid_params%yy(jgm))
                        case(3)
cc                          car0 = map(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
cc                          carp = map(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
cc                          carm = map(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                          call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                          call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                          call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                          car0 = xc(i0,j0,k0,:)
                          carp = xc(i0,j0,kp,:)
                          carm = xc(i0,j0,km,:)
                          dhp = (grid_params%zz(kgp)
     .                          -grid_params%zz(kg0))
                          dhm = (grid_params%zz(kg0)
     .                          -grid_params%zz(kgm))
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
                         dhh = 0.5*(grid_params%xx(igp)
     .                             -grid_params%xx(igm))
                        case(2)
                         call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                         call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                         call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                         carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
                         carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
                         dhh = 0.5*(grid_params%yy(jgp)
     .                             -grid_params%yy(jgm))
                        case(3)
                         call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                         call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                         call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                         carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
                         carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
                         dhh = 0.5*(grid_params%zz(kgp)
     .                             -grid_params%zz(kgm))
                        end select
                      endif
                      vec = (carp-carm)/dhh

                      gamma(l,m,n) = dot_product(vec,cov(l,:))
                      gamma(l,n,m) = gamma(l,m,n) !Symmetry
                    enddo
                  enddo
                enddo

                !Store grid quantities
                gmetric%grid(igrid)%jac  (i,j,k)       = jac
                gmetric%grid(igrid)%gsub (i,j,k,:,:)   = gsub
                gmetric%grid(igrid)%gsup (i,j,k,:,:)   = gsup
                gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
                gmetric%grid(igrid)%cov  (i,j,k,:,:)   = cov
                gmetric%grid(igrid)%cnv  (i,j,k,:,:)   = cnv

              enddo
            enddo
          enddo

          deallocate(dr,xc)

          !Zero-force condition on Christoffle symbols 
          if (igrid == 1) then
            call gammaZeroForce(igx,igy,igz)
          endif

          !Enforce topological constraints on Christoffel symbols
          do i=1,3
            do j=1,3
              do k=1,3
                call topol_bc(gmetric%grid(igrid)%Gamma(:,:,:,i,j,k)
     .                       ,igrid)
              enddo
            enddo
          enddo

        endif

c       Store cell volume

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
             gmetric%grid(igrid)%dvol(i,j,k)=volume(i,j,k,igx,igy,igz)
            enddo
          enddo
        enddo

c     End program

      end subroutine findGridMetrics

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

        integer(4) :: i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq
        real(8)    :: dr(3,3),x(0:nx+1,0:ny+1,0:nz+1,neq)

c     Local variables

        integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
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

        dh(1)= (grid_params%xx(igp)-grid_params%xx(igm))

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

        dh(2)= (grid_params%yy(jgp)-grid_params%yy(jgm))

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

        dh(3)= (grid_params%zz(kgp)-grid_params%zz(kgm))

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
cc        integer(4) :: i,j,k,igx,igy,igz,hex,hey,hez,nx,ny,nz,neq
cc        real(8)    :: x(0:nx+1,0:ny+1,0:nz+1,neq),jac
cc
cc        real(8),optional,intent(OUT) :: ggsup(3,3),ggsub(3,3),vol
cc     .                                 ,ccov(3,3),ccnv(3,3)
cc        real(8),optional,intent(IN)  :: dtensor(3,3)
cc
ccc     Local variables
cc
cc        integer(4) :: m,l,ig,jg,kg
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
cc          vol = jac*grid_params%dxh(ig)
cc     .             *grid_params%dyh(jg)
cc     .             *grid_params%dzh(kg)
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

c     coarsenGridMetrics
c     #################################################################
      subroutine coarsenGridMetrics(igc)

c     -----------------------------------------------------------------
c     Coasens grid metrics to grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: igc

c     Local variables

        integer(4) :: igf,ic,jc,kc,if,jf,kf,icg,jcg,kcg
     .               ,nxc,nyc,nzc,nxf,nyf,nzf
     .               ,mgx,mgy,mgz,l,m

        real(8)    :: volc

c     Begin program

        !From coarse grid (good to find ghost cells)
        call findGridMetrics(igc)
        return

c       Coarsen volumes and jacobians (needed for coarsening all else)

        igf = igc - 1

        nxc = grid_params%nxv(igc)
        nyc = grid_params%nyv(igc)
        nzc = grid_params%nzv(igc)

        mgx = grid_params%mg_ratio_x(igf)
        mgy = grid_params%mg_ratio_y(igf)
        mgz = grid_params%mg_ratio_z(igf)

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

              volc   = 0d0
              do kf = mgz*(kc-1)+1,mgz*kc
                do jf = mgy*(jc-1)+1,mgy*jc
                  do if = mgx*(ic-1)+1,mgx*ic
                    volc = volc + gmetric%grid(igf)%dvol(if,jf,kf)
                  enddo
                enddo
              enddo

              gmetric%grid(igc)%dvol(ic,jc,kc) = volc

              call getMGmap(ic,jc,kc,igc,igc,igc,icg,jcg,kcg)

              gmetric%grid(igc)%jac (ic,jc,kc) = volc
     .                                          /(grid_params%dxh(icg)
     .                                           *grid_params%dyh(jcg)
     .                                           *grid_params%dzh(kcg))

            enddo
          enddo
        enddo

c       Coarsen everything else

        !Gsup
        do l=1,3
          do m=l,3
            call coarsen(igc,gmetric%grid(igf)%gsup(:,:,:,l,m)
     .                      ,gmetric%grid(igc)%gsup(:,:,:,l,m),.true.)
            gmetric%grid(igc)%gsup(:,:,:,m,l)
     .           = gmetric%grid(igc)%gsup(:,:,:,l,m) !Symmetry
          enddo
        enddo

        !Gsub
        do l=1,3
          do m=l,3
            call coarsen(igc,gmetric%grid(igf)%gsub(:,:,:,l,m)
     .                      ,gmetric%grid(igc)%gsub(:,:,:,l,m),.true.)
            gmetric%grid(igc)%gsub(:,:,:,m,l)
     .           = gmetric%grid(igc)%gsub(:,:,:,l,m) !Symmetry
          enddo
        enddo

        !Gamma
        do m=1,3
          do n=m,3
            do l=1,3
              call coarsen(igc,gmetric%grid(igf)%Gamma(:,:,:,l,m,n)
     .                        ,gmetric%grid(igc)%Gamma(:,:,:,l,m,n)
     .                    ,.true.)
            enddo
            gmetric%grid(igc)%Gamma(:,:,:,:,n,m)
     .           = gmetric%grid(igc)%Gamma(:,:,:,:,m,n) !Symmetry
          enddo
        enddo

        !Cov
        do l=1,3
          do m=1,3
            call coarsen(igc,gmetric%grid(igf)%cov(:,:,:,l,m)
     .                      ,gmetric%grid(igc)%cov(:,:,:,l,m),.true.)
          enddo
        enddo

        !Cnv
        do l=1,3
          do m=1,3
            call coarsen(igc,gmetric%grid(igf)%cnv(:,:,:,l,m)
     .                      ,gmetric%grid(igc)%cnv(:,:,:,l,m),.true.)
          enddo
        enddo

c     End program

      end subroutine coarsenGridMetrics

c     coarsen
c     #################################################################
      subroutine coarsen(igc,arrayf,arrayc,volf)

c     -----------------------------------------------------------------
c     Coarsen arrayf to arrayc at grid level igc.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: igc
        real(8)    :: arrayf(:,:,:),arrayc(:,:,:)
        logical    :: volf

c     Local variables

        integer(4) :: igf,ic,jc,kc,if,jf,kf
     .               ,nxc,nyc,nzc,mgx,mgy,mgz

        real(8)    :: volc,ac

c     Begin program

        igf = igc - 1

        nxc = grid_params%nxv(igc)
        nyc = grid_params%nyv(igc)
        nzc = grid_params%nzv(igc)

        mgx = grid_params%mg_ratio_x(igf)
        mgy = grid_params%mg_ratio_y(igf)
        mgz = grid_params%mg_ratio_z(igf)

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

              volc   = 0d0
              ac     = 0d0

              do kf = mgz*(kc-1)+1,mgz*kc
                do jf = mgy*(jc-1)+1,mgy*jc
                  do if = mgx*(ic-1)+1,mgx*ic

                    if (volf) then
                      ac   = ac   + gmetric%grid(igf)%dvol(if,jf,kf)
     .                             *arrayf(if,jf,kf)
                      volc = volc + gmetric%grid(igf)%dvol(if,jf,kf)
                    else
                      ac   = ac   + arrayf(if,jf,kf)
                      volc = volc + 1d0
                    endif

                  enddo
                enddo
              enddo

cc              if (volf) ac = ac/volc

              arrayc(ic,jc,kc) = ac/volc
            enddo
          enddo
        enddo

        call pstop('coarsen'
     .            ,'Routine does not fill ghost cells')

      end subroutine coarsen

c     topol_bc
c     #################################################################
      subroutine topol_bc(array,igrid)

c     -----------------------------------------------------------------
c     !Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        real(8)    :: array(ilom:ihip,jlom:jhip,klom:khip)

        integer(4) :: igrid

c     Local variables

        integer(4) :: nxl,nyl,nzl,nxg,nyg,nzg,dim

c     External

c     Begin program

#if !defined(petsc)

        nxl = grid_params%nxv(igrid)
        nyl = grid_params%nyv(igrid)
        nzl = grid_params%nzv(igrid)

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
cc        if (igrid == 1) then
          call fillPetscGhostCells(array,igrid)
cc        endif
#endif

c     End program

      end subroutine topol_bc

c     gammaZeroForce
c     #################################################################
      subroutine gammaZeroForce(igx,igy,igz)

c     -----------------------------------------------------------------
c     Postprocess Christoffel symbol to satisfy zero-force condition
c     to machine accuracy
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer(4) :: igx,igy,igz

c     Local variables

        integer(4) :: i,j,k,ii,ll,kk,mm,jj,ig,jg,kg,igrid
        real(8)    :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),table(3,3,3),gsub(3,3)
     .               ,gsup(3,3),gsupp(3,3),gsupm(3,3),gamma(3,3,3)
        logical    :: alt_eom

c     Begin program

cc        if  (    (i == ilom .or. i == ihip)       !Avoid x-boundaries
cc     .       .or.(j == jlom .or. j == jhip)       !Avoid y-boundaries
cc     .       .or.(k == klom .or. k == khip)       !Avoid z-boundaries
cccc     .       .or. isSP(i,j,k,igx,igy,igz)         !Avoid singular point
cc     .       ) return

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

            gamma = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
            gsub  = gmetric%grid(igrid)%gsub(i,j,k,:,:)
            gsup  = gmetric%grid(igrid)%gsup(i,j,k,:,:)

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            do ii=1,3 !cycle through Chistoffel symbols

              alt_eom = (ii == 2 .and. alt__eom())

              if (alt_eom) then
                jac  = gmetric%grid(igrid)%jac(i,j,k)
                const= 2d0
              else
                jac  = 1d0
                const= 1d0
              endif

              do ll=1,3
                do kk=1,3
                  select case (kk)
                  case(1)
                    dh = 2*grid_params%dxh(ig)
                    if (alt_eom) then
                      jacp = gmetric%grid(igrid)%jac(i+1,j,k)
                      jacm = gmetric%grid(igrid)%jac(i-1,j,k)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp*gmetric%grid(igrid)%gsup (i+1,j,k,:,:)
                    gsupm = jacm*gmetric%grid(igrid)%gsup (i-1,j,k,:,:)
                  case(2)
                    dh = 2*grid_params%dyh(jg)
                    if (alt_eom) then
                      jacp = gmetric%grid(igrid)%jac(i,j+1,k)
                      jacm = gmetric%grid(igrid)%jac(i,j-1,k)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp*gmetric%grid(igrid)%gsup (i,j+1,k,:,:)
                    gsupm = jacm*gmetric%grid(igrid)%gsup (i,j-1,k,:,:)
                  case(3)
                    dh = 2*grid_params%dzh(kg)
                    if (alt_eom) then
                      jacp = gmetric%grid(igrid)%jac(i,j,k+1)
                      jacm = gmetric%grid(igrid)%jac(i,j,k-1)
                    else
                      jacp=1d0
                      jacm=1d0
                    endif
                    gsupp = jacp*gmetric%grid(igrid)%gsup (i,j,k+1,:,:)
                    gsupm = jacm*gmetric%grid(igrid)%gsup (i,j,k-1,:,:)
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

            enddo

            !Symmetrize
            do ii=1,3
              do ll=1,3             !cycle through tensor components
                do kk=1,3
                  gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
     $                        = 0.5*( table(ii,ll,kk)+table(ii,kk,ll) )
                enddo
              enddo
            enddo

c           START CHECKS

            if (check_grid) call zeroForceCheck(i,j,k,igrid,gamma
     $                     ,gmetric%grid(igrid)%Gamma(i,j,k,:,:,:))

            enddo
          enddo
        enddo

cc        !Postprocessing for alt_eom case (ii == 2 .and. bcond(1) == SP) 
cc        if (bcSP()) then
cc
cc          gamma = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cc          jac   = gmetric%grid(igrid)%jac(i,j,k)
cc
cc          do ii=1,3             !cycle through Chistoffel symbols
cc            do ll=1,3
cc              do kk=1,3
cc                select case (kk)
cc                case(1)
cc                  dh = grid_params%dxh(ig)
cc                  jacp = hlf*(jac+gmetric%grid(igrid)%jac(i+1,j,k))
cc                  jacm = hlf*(jac+gmetric%grid(igrid)%jac(i-1,j,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i-1,j,k,:,:))
cc                case(2)
cc                  dh = grid_params%dyh(jg)
cc                  jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j+1,k))
cc                  jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j-1,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j-1,k,:,:))
cc                case(3)
cc                  dh = grid_params%dzh(kg)
cc                  jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k+1))
cc                  jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k-1))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k-1,:,:))
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
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,kk,2)-table(3,kk,3)
cc                  case(2)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(3,kk,3)
cc                  case(3)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(2,kk,2)
cc                  end select
cc                elseif(ii==kk .and.ii/=ll) then
cc                  select case(ii)
cc                  case(1)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,2,ll)-table(3,3,ll)
cc                  case(2)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(3,3,ll)
cc                  case(3)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(2,2,ll)
cc                  end select
cc                elseif(ii==kk .and.ii==ll) then
cc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = gamma(ii,kk,ll)
cc                else
cc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = table(ii,kk,ll)
cc                endif
cccc                if (ii == ll) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     .                 = matmul(inv,table(ii,:,ii))
cccccc     $                 = hlf*(table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccccc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cccccc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccc                  case(2)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccccc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cccccc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccc                  case(3)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccccc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cccccc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccc                  end select
cccc                elseif (ii == kk) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccccc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cccccc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccc                  case(2)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccccc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cccccc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccc                  case(3)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccccc     $                 =0.25*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll)
cccccc     $                       +table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccc                  end select
cccc                else
cccc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cccc                endif
cc              enddo
cc            enddo
cc
cc          enddo

cc                if (ii /= ll .and. ii /= kk) then
cc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = table(ii,ll,kk)
cc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cc                else
cc                  if (ii == ll) then
cc                    select case(ii)
cc                    case(1)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk))
cc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cc                    case(2)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk))
cc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cc                    case(3)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk))
cc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cc                    end select
cc                  else
cc                    select case(ii)
cc                    case(1)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll))
cc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cc                    case(2)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll))
cc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cc                    case(3)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
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

        integer(4) :: i,j,k,igrid
        real(8)    :: gold(3,3,3),gnew(3,3,3)
        
c     Local variables

        integer(4) :: ii,ll,kk,mm,jj,ig,jg,kg
        real(8)    :: summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),gsup(3,3),gsub(3,3)
     $               ,gsupp(3,3),gsupm(3,3)
        logical    :: alt_eom

c     Begin program

        write (*,*)
        write (*,*) 'Grid:',igrid,'  Grid node:',i,j,k

        call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

        gsup  = gmetric%grid(igrid)%gsup(i,j,k,:,:)
        gsub  = gmetric%grid(igrid)%gsub(i,j,k,:,:)

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
            jac  = gmetric%grid(igrid)%jac(i,j,k)
            const= 1d0
          else
            jac  = 1d0
            const= 0d0
          endif

          summ = 0d0

          dh = grid_params%dxh(ig)
          if (alt_eom) then
            jacp = hlf*(jac+gmetric%grid(igrid)%jac(i+1,j,k))
            jacm = hlf*(jac+gmetric%grid(igrid)%jac(i-1,j,k))
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
          gsupm = jacm
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i-1,j,k,:,:))

          summ = summ +(gsupp(1,ii)-gsupm(1,ii))/dh/jac

          dh = grid_params%dyh(jg)
          if (alt_eom) then
            jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j+1,k))
            jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j-1,k))
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
          gsupm = jacm
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j-1,k,:,:))

          summ = summ +(gsupp(2,ii)-gsupm(2,ii))/dh/jac

          dh = grid_params%dzh(kg)
          if (alt_eom) then
            jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k+1))
            jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k-1))
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
          gsupm = jacm
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k-1,:,:))

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
        write (*,*) 'Cannot check numerical grid quantities: analytical'
        write (*,*) '  mapping unavailable'
      endif

      do igrid=1,grid_params%ngrid

        igx = igrid
        igy = igrid
        igz = igrid

        mag1 = 0d0
        mag2 = 0d0
        mag3 = 0d0
        mag4 = 0d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              mag1 = mag1 + (gmetric%grid(igrid)%jac  (i,j,k)
     .                      -jacobian(i,j,k,igx,igy,igz))**2
              mag2 = mag2 + (gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                      -g_sub   (i,j,k,igx,igy,igz))**2
              mag3 = mag3 + (gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                      -g_sup   (i,j,k,igx,igy,igz))**2
              mag4 = mag4 + (gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .                      -christ_2knd(i,j,k,igx,igy,igz))**2
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

      use debug_mod

      implicit none

      integer(4) :: nx,ny,nz,i,j,k,l,m,n,igrid
      real(8)    :: tensor(3,3,3)
      character(17) :: debugfile

c     Begin program

      if (anal_map) then
        write (*,*)
        write (*,*) '**********************************'
        write (*,*) 'Dumping grid debugging graphics...'
        write (*,*) '**********************************'
        write (*,*)
      else
        write (*,*) 'Cannot dump debugging graphics: analytical'
        write (*,*) '  mapping unavailable'
        return
      endif

      igrid = 1

      debugfile = 'grid_debug'//trim(int2char(my_rank))//'.bin'
      open(unit=110,file=trim(debugfile),form='unformatted'
     .    ,status='replace')

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      allocate(dbg(0:nx+1,0:ny+1,0:nz+1,1))

      !Jacobian
      k = 1
      do j=0,ny+1
        do i=0,nx+1
          dbg(i,j,1,1) = gmetric%grid(igrid)%jac(i,j,k)
     .                  -jacobian(i,j,k,igrid,igrid,igrid)
        enddo
      enddo
      call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,0,110)

      !Metric tensors
      do l=1,3
        do m=1,3
          k = 1
          do j=0,ny+1
            do i=0,nx+1
              tensor(:,:,1)  = g_sub(i,j,k,igrid,igrid,igrid)
              dbg(i,j,k,1) = gmetric%grid(igrid)%gsub(i,j,k,l,m)
     .                  *gmetric%grid(igrid)%jac(i,j,k)
     .                      -tensor(l,m,1)
     .                  *jacobian(i,j,k,igrid,igrid,igrid)
            enddo
          enddo
          call contour(dbg(:,:,k,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1,110)
        enddo
      enddo

      !Christoffle symbols
      do l=1,3
        do m=1,3
          do n=1,3
            k = 1
            do j=0,ny+1
              do i=0,nx+1
                tensor  = christ_2knd(i,j,k,igrid,igrid,igrid)
                dbg(i,j,k,1) = gmetric%grid(igrid)%Gamma(i,j,k,l,m,n)
     .                  *gmetric%grid(igrid)%jac(i,j,k)
     .                        -tensor(l,m,n)
     .                  *jacobian(i,j,k,igrid,igrid,igrid)
              enddo
            enddo
            call contour(dbg(:,:,k,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1,110)
          enddo
        enddo
      enddo

      deallocate(dbg)
      close(110)

      write (*,*) 'Done...'
      write (*,*) 

      end subroutine graphicsCheck

      end subroutine defineGridMetric

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

        integer(4) :: i,j
        real(8)    :: delta

        if (i == j) then
          delta = 1d0
        else
          delta = 0d0
        endif

      end function delta

c     volume
c     #################################################################
      function volume(i,j,k,igx,igy,igz) result(vol)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: vol

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,dx1,dx2,dx3,jac

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dx1 = grid_params%dxh(ig)
        dx2 = grid_params%dyh(jg)
        dx3 = grid_params%dzh(kg)

        jac = gmetric%grid(igx)%jac(i,j,k)
        
        vol = jac*dx1*dx2*dx3

      end function volume

      end module grid_metric

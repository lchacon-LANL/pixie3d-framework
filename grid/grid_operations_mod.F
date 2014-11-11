
c module grid_operations
c #####################################################################
      module grid_operations

c ---------------------------------------------------------------------
c     This module packs routines that perform operations on grid
c     quantities, such as coordinate transformation of vector components,
c     vector norms and scalar products. It contains the following
c     routines:
c          * transformVectorToCartesian 
c          * transformVectorToCurvilinear 
c          * transformFromCurvToCurv 
c          * vectorNorm
c          * scalarProduct
c     It is assumed that the grid metric structure gmetric is allocated 
c     and filled.
c ---------------------------------------------------------------------

        use grid_def

        implicit none

        INTERFACE crossProduct
          module procedure crossProduct_ijk,crossProduct_mesh
        end INTERFACE

        INTERFACE vectorNorm
          module procedure vectorNorm_ijk,vectorNorm_mesh
        end INTERFACE

        INTERFACE scalarProduct
          module procedure scalarProduct_ijk,scalarProduct_mesh
        end INTERFACE

        INTERFACE XformToCar
          module procedure XformToCar_ijk,XformToCar_mesh
        end INTERFACE

        INTERFACE XformToCurv
          module procedure XformToCurv_ijk,XformToCurv_mesh
        end INTERFACE

        INTERFACE XformVector
          module procedure XformVector_mesh!,XformVector_mesh2
        end INTERFACE

        INTERFACE XformVector_stg
          module procedure XformVector_stg_mesh
        end INTERFACE

cc        private :: transformVectorToCartesian
cc     .            ,transformVectorToCurvilinear
cc     .            ,transformFromCurvToCurv

      contains

c     getMGmap
c     #################################################################
      subroutine getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c     -----------------------------------------------------------------
c     Gets MG vector components (ig,jg,kg) for grid quantities
c     corresponding to node position (i,j,k) in grid levels igx,igy,igz
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz,ig,jg,kg

c     Local variables

c     Begin program

        ig = i + grid_params%istartx(igx)
        jg = j + grid_params%istarty(igy)
        kg = k + grid_params%istartz(igz)

      end subroutine getMGmap

c     getCurvilinearCoordinates
c     #################################################################
      subroutine getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Finds curvilinear coordinates for position (i,j,k)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        integer    :: ii,jj,ny

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x1 = grid_params%xx(ig)
        y1 = grid_params%yy(jg)
        z1 = grid_params%zz(kg)

      end subroutine getCurvilinearCoordinates

c     getCartesianCoordinates
c     #################################################################
      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Returns Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        real(8)    :: car(3)

c     Begin program

        !Find coordinates
        x1 = gmetric%grid(igx)%car(i,j,k,1)
        y1 = gmetric%grid(igy)%car(i,j,k,2)
        z1 = gmetric%grid(igz)%car(i,j,k,3)

      end subroutine getCartesianCoordinates

c     vecPos
c     #################################################################
      function vecPos(neq,i,j,k,igx,igy,igz) result (ijk)

        implicit none

        integer    :: ijk,i,j,k,neq,igx,igy,igz

        integer    :: nxl,nyl,nzl,ilo,jlo,klo

        ilo = 1
        jlo = 1
        klo = 1

        nxl = grid_params%nxv(igx)
        nyl = grid_params%nyv(igy)
        nzl = grid_params%nzv(igz)

        ijk = neq*(i-ilo + nxl*(j-jlo) + nxl*nyl*(k-klo))
      
      end function vecPos

c     fromGlobalToLocalLimits
c     #################################################################
      subroutine fromGlobalToLocalLimits(ig,jg,kg,il,jl,kl,igx,igy,igz)

        implicit none

        integer    :: ig,jg,kg,il,jl,kl,igx,igy,igz

        il = ig - grid_params%ilo(igx) + 1
        jl = jg - grid_params%jlo(igy) + 1
        kl = kg - grid_params%klo(igz) + 1

      end subroutine fromGlobalToLocalLimits

c     fromLocalToGlobalLimits
c     #################################################################
      subroutine fromLocalToGlobalLimits(il,jl,kl,ig,jg,kg,igx,igy,igz)

        implicit none

        integer    :: ig,jg,kg,il,jl,kl,igx,igy,igz

        ig = il + grid_params%ilo(igx) - 1
        jg = jl + grid_params%jlo(igy) - 1
        kg = kl + grid_params%klo(igz) - 1

      end subroutine fromLocalToGlobalLimits

c     localIndex
c     #################################################################
      function localIndex(igrid,gindex) result(lindex)

        implicit none

        integer  :: igrid,lindex(3),gindex(3)

        lindex(1) = gindex(1) - grid_params%ilo(igrid) + 1
        lindex(2) = gindex(2) - grid_params%jlo(igrid) + 1
        lindex(3) = gindex(3) - grid_params%klo(igrid) + 1

      end function localIndex

c     globalIndex
c     #################################################################
      function globalIndex(igrid,lindex) result(gindex)

        implicit none

        integer  :: igrid,lindex(3),gindex(3)

        gindex(1) = lindex(1) + grid_params%ilo(igrid) - 1
        gindex(2) = lindex(2) + grid_params%jlo(igrid) - 1
        gindex(3) = lindex(3) + grid_params%klo(igrid) - 1

      end function globalIndex

c     XformToCar_ijk
c     #################################################################
      function XformToCar_ijk(i,j,k,igrid,curv,covariant) result(car)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igrid
        real(8)    :: curv(3),car(3)
        logical    :: covariant

c     Local variables

c     Begin program

        if (covariant) then
          car = matmul(curv,gmetric%grid(igrid)%cov(i,j,k,:,:))
        else
          car = matmul(curv,gmetric%grid(igrid)%cnv(i,j,k,:,:))
        endif

      end function XformToCar_ijk

c     XformToCar_mesh
c     #################################################################
      function XformToCar_mesh(igrid,curv,covariant) result(car)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector to Cartesian 
c     over the whole mesh. Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igrid
        real(8) :: curv(0:,0:,0:,:)
     .            ,car(0:size(curv,1)-1
     .                ,0:size(curv,2)-1
     .                ,0:size(curv,3)-1,size(curv,4))
        logical :: covariant

c     Local variables

        integer :: nxl,nyl,nzl

c     Begin program

        nxl = grid_params%nxv(igrid)
        nyl = grid_params%nyv(igrid)
        nzl = grid_params%nzv(igrid)

        do k=0,nzl+1
          do j=0,nyl+1
            do i=0,nxl+1
              car(i,j,k,:) = XformToCar_ijk(i,j,k,igrid,curv(i,j,k,:)
     .                                     ,covariant)
            enddo
          enddo
        enddo

      end function XformToCar_mesh

c     XformToCurv_ijk
c     #################################################################
      function XformToCurv_ijk(i,j,k,igrid,car,tocnv,nojac) result(curv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector to curvilinear representation at
c     grid coordinates (i,j,k). Curvilinear vector is covariant if
c     tocnv=.false., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igrid
        real(8) :: car(3),curv(3)
        logical :: tocnv

        logical,optional :: nojac

c     Local variables

        real(8) :: jac
        logical :: nj

c     Begin program

        if (PRESENT(nojac)) then
          nj = nojac
        else
          nj = .false.
        endif

        if (nj) then
          jac = 1d0
        else
          jac = gmetric%grid(igrid)%jac(i,j,k)
        endif

        if (.not.tocnv) then
          curv = jac*matmul(gmetric%grid(igrid)%cnv(i,j,k,:,:),car)
        else
          curv = jac*matmul(gmetric%grid(igrid)%cov(i,j,k,:,:),car)
        endif

      end function XformToCurv_ijk

c     XformToCurv_mesh
c     #################################################################
      function XformToCurv_mesh(igrid,car,tocnv) result(curv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector to curvilinear
c     over the whole mesh. Curvilinear vector is contravariant if
c     tocnv=.true., and covariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igrid
        real(8) :: car(0:,0:,0:,:)
     .            ,curv(0:size(car,1)-1
     .                 ,0:size(car,2)-1
     .                 ,0:size(car,3)-1,size(car,4))
        logical :: tocnv

c     Local variables

        integer :: nxl,nyl,nzl

c     Begin program

        nxl = grid_params%nxv(igrid)
        nyl = grid_params%nyv(igrid)
        nzl = grid_params%nzv(igrid)

        do k=0,nzl+1
          do j=0,nyl+1
            do i=0,nxl+1
              curv(i,j,k,:) = XformToCurv_ijk(i,j,k,igrid,car(i,j,k,:)
     .                                       ,tocnv)
            enddo
          enddo
        enddo

      end function XformToCurv_mesh

c     XformToCurv_mesh2
c     #################################################################
      subroutine XformToCurv_mesh2(igrid,vec,tocnv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector vec to curvilinear
c     over the whole mesh. Curvilinear vector is contravariant if
c     tocnv=.true., and covariant otherwise. Result is returned in vec.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igrid
        real(8) :: vec(0:,0:,0:,:)
        logical :: tocnv

c     Local variables

        integer :: nxl,nyl,nzl
        real(8) :: dum(size(vec,4))

c     Begin program

        nxl = grid_params%nxv(igrid)
        nyl = grid_params%nyv(igrid)
        nzl = grid_params%nzv(igrid)

        do k=0,nzl+1
          do j=0,nyl+1
            do i=0,nxl+1
              dum = vec(i,j,k,:)
              vec(i,j,k,:) = XformToCurv_ijk(i,j,k,igrid,dum,tocnv)
            enddo
          enddo
        enddo

      end subroutine XformToCurv_mesh2

c     transformVectorToCurvilinear
c     #################################################################
      subroutine transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .                                       ,cx,cy,cz,tocnv
     .                                       ,c1,c2,c3)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector (cx,cy,cz) to curvilinear (c1,c2,c3)
c     (contravariant if tocnv=.true., contravariant otherwise)
c     at grid coordinates (i,j,k). 
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: tocnv

c     Local variables

        integer    :: ic
        real(8)    :: T_to_curv(3,3),vec(3),jac

c     Begin program

        jac = gmetric%grid(igx)%jac(i,j,k)

        vec = (/ cx,cy,cz /)

        if (.not.tocnv) then

          c1 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,1,:))
          c2 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,2,:))
          c3 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,3,:))
            
        else

          c1 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,1,:))
          c2 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,2,:))
          c3 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,3,:))

        endif

      end subroutine transformVectorToCurvilinear

c     transformFromCurvToCurv
c     #################################################################
      subroutine transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv,half_elem
     .             ,hex,hey,hez)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at cell centers, faces, and vertices
c     of grid coordinate (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        integer   ,optional,intent(IN) :: half_elem,hex,hey,hez

        real(8)    :: cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical    :: tocnv

c     Local variables

        integer    :: he,he_x,he_y,he_z,ip,jp,kp
        real(8)    :: cov(3),cnv(3),tensor(3,3),coeff
        real(8),pointer,dimension(:,:,:,:,:) :: gtnsr

c     Begin program

c     Special case of cartesian coords

        if (no_map) then

          if (tocnv) then
            cnv1 = cov1
            cnv2 = cov2
            cnv3 = cov3
          else
            cov1 = cnv1
            cov2 = cnv2
            cov3 = cnv3
          endif

          return
        endif

c     Optional parameters

        if (PRESENT(half_elem)) then
          he = half_elem
        else
          he = 0
        endif

        if (PRESENT(hex)) then
          he_x = hex
        else
          he_x = 0
        endif

        if (PRESENT(hey)) then
          he_y = hey
        else
          he_y = 0
        endif

        if (PRESENT(hez)) then
          he_z = hez
        else
          he_z = 0
        endif

c     Consistency check

        if (he == 1) then
          he_x = 1; he_y = 0 ; he_z = 0
        endif

        if (he == 2) then
          he_x = 0; he_y = 1 ; he_z = 0
        endif

        if (he == 3) then
          he_x = 0; he_y = 0 ; he_z = 1
        endif

c     Find metric tensor

        ip = i+1
        jp = j+1
        kp = k+1

        if (tocnv) then
          gtnsr => gmetric%grid(igx)%gsup
        else
          gtnsr => gmetric%grid(igx)%gsub
        endif

        tensor = gtnsr(i,j,k,:,:)
        coeff = 1d0

        if (he_x == 1) then

          tensor = tensor + gtnsr(ip,j,k,:,:)
          coeff = 0.5

          if (he_y == 1) then
            tensor = tensor + gtnsr(i ,jp,k,:,:)
     .                      + gtnsr(ip,jp,k,:,:)
            coeff = 0.25
          endif

          if (he_z == 1) then
            tensor = tensor + gtnsr(i ,j ,kp,:,:)
     .                      + gtnsr(ip,j ,kp,:,:)
     .                      + gtnsr(i ,jp,kp,:,:)
     .                      + gtnsr(ip,jp,kp,:,:)
            coeff = 0.125
          endif

        elseif (he_y == 1) then

          tensor = tensor + gtnsr(i,jp,k,:,:)
          coeff = 0.5

          if (he_z == 1) then
            tensor = tensor + gtnsr(i,j ,kp,:,:)
     .                      + gtnsr(i,jp,kp,:,:)
            coeff = 0.25
          endif

        elseif (he_z == 1) then

          tensor = tensor + gtnsr(i,j,kp,:,:)
          coeff = 0.5

        endif

        tensor = coeff*tensor

c     Transform vector

        if (tocnv) then
          cnv1 = tensor(1,1)*cov1+tensor(1,2)*cov2+tensor(1,3)*cov3
          cnv2 = tensor(2,1)*cov1+tensor(2,2)*cov2+tensor(2,3)*cov3
          cnv3 = tensor(3,1)*cov1+tensor(3,2)*cov2+tensor(3,3)*cov3
        else
          cov1 = tensor(1,1)*cnv1+tensor(1,2)*cnv2+tensor(1,3)*cnv3
          cov2 = tensor(2,1)*cnv1+tensor(2,2)*cnv2+tensor(2,3)*cnv3
          cov3 = tensor(3,1)*cnv1+tensor(3,2)*cnv2+tensor(3,3)*cnv3
        endif

      end subroutine transformFromCurvToCurv

c     transformVector
c     ######################################################################
      subroutine transformVector(igx,igy,igz
     .                          ,imin,imax,jmin,jmax,kmin,kmax
     .                          ,arr1,arr2,arr3,in_repr,out_repr)

c     ----------------------------------------------------------------------
c     Transforms vectors components in arrays arr1,arr2,arr3 from
c     Cartesian to curvilinear and viceversa, or from curvilinear to
c     curvilinear, depending on the value of in_repr (which indicates the
c     representation of the input vector) and out_repr (which indicates
c     representation of the input vector). Possible representations are:
c     'car','cnv',and 'cov'.
c     ----------------------------------------------------------------------

      implicit none

c     Input variables

        integer    :: imin,imax,jmin,jmax,kmin,kmax
        integer    :: igx,igy,igz
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

        character(3) :: in_repr,out_repr

c     Local variables

        integer    :: i,j,k
        real(8)    :: vec(3)

        logical    :: covariant,cnvariant,cartesian

c     Begin program

        cartesian = (in_repr == 'car')
        covariant = (in_repr == 'cov')
        cnvariant = (.not. cartesian) .and. (.not. covariant)

        select case(out_repr)
        case('car')

          if (cartesian) return

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                vec = XformToCar(i,j,k,igx
     .               ,(/arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)/)
     .               ,covariant)

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        case('cnv')

          if (cnvariant) return

          if (cartesian) then
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformVectorToCurvilinear
     .                 (i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,.true.
     .                 ,vec(1),vec(2),vec(3))

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          else
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformFromCurvToCurv
     .                 (i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,vec(1),vec(2),vec(3)
     .                 ,.true.,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        case('cov')

          if (covariant) return

          if (cartesian) then
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformVectorToCurvilinear
     .                 (i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,.false.
     .                 ,vec(1),vec(2),vec(3))

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          else
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformFromCurvToCurv
     .                 (i,j,k,igx,igy,igz
     .                 ,vec(1),vec(2),vec(3)
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,.false.,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        end select

      end subroutine transformVector

c     XformVector_mesh2
c     #################################################################
      subroutine XformVector_mesh2(igrid,vec,in_repr,out_repr)

c     -----------------------------------------------------------------
c     Transforms vector vec from Cartesian to curvilinear and viceversa,
c     or from curvilinear to curvilinear, depending on the value of
c     in_repr (which indicates the representation of the input vector)
c     and out_repr (which indicates representation of the input
c     vector). Possible representations are: 'car','cnv',and 'cov'.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igrid
        real(8) :: vec(:,:,:,:)

        character(3) :: in_repr,out_repr

c     Local variables

        integer :: imin,imax,jmin,jmax,kmin,kmax

c     Begin program

        imin = 0
        imax = grid_params%nxv(igrid)+1

        jmin = 0
        jmax = grid_params%nyv(igrid)+1

        kmin = 0
        kmax = grid_params%nzv(igrid)+1

        call  transformVector(igrid,igrid,igrid
     .                       ,imin,imax,jmin,jmax,kmin,kmax
     .                       ,vec(:,:,:,1)
     .                       ,vec(:,:,:,2)
     .                       ,vec(:,:,:,3),in_repr,out_repr)

c     End program

      end subroutine XformVector_mesh2

c     XformVector_mesh
c     #################################################################
      function XformVector_mesh(igrid,vec1,tocnv) result(vec2)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at cell centers on all mesh.
c     -----------------------------------------------------------------

      implicit none

c     Input variables

      integer    :: igrid
      real(8)    :: vec1(0:,0:,0:,:)
      real(8)    :: vec2(0:size(vec1,1)-1
     .                  ,0:size(vec1,2)-1
     .                  ,0:size(vec1,3)-1
     .                  ,size(vec1,4))

      logical    :: tocnv

c     Local variables

      integer    :: i,j,k,nxl,nyl,nzl
      real(8),pointer,dimension(:,:,:,:,:) :: tensor

c     Begin program

      nxl = grid_params%nxv(igrid)
      nyl = grid_params%nyv(igrid)
      nzl = grid_params%nzv(igrid)

      if (tocnv) then
        tensor => gmetric%grid(igrid)%gsup
      else
        tensor => gmetric%grid(igrid)%gsub
      endif

      do k=0,nzl+1
        do j=0,nyl+1
          do i=0,nxl+1
            vec2(i,j,k,:) = matmul(tensor(i,j,k,:,:),vec1(i,j,k,:))
          enddo
        enddo
      enddo

      nullify(tensor)

c     End program

      end function XformVector_mesh

c     XformVector_BC
c     ###################################################################
      subroutine XformVector_BC(igrid,vec1,vec2,tocnv)

      implicit none

c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form vec1 to vec2
c     -------------------------------------------------------------------

c     Input variables

      integer :: igrid
      real(8) :: vec1(0:,0:,0:,:)
      real(8) :: vec2(0:,0:,0:,:)
      logical :: tocnv

c     Local variables

      integer :: nxp,nyp,nzp,i,j,k
      real(8),pointer,dimension(:,:,:,:,:) :: tensor

c     Begin program

      nxp = size(vec1,1)-1
      nyp = size(vec1,2)-1
      nzp = size(vec1,3)-1

      if (tocnv) then
        tensor => gmetric%grid(igrid)%gsup
      else
        tensor => gmetric%grid(igrid)%gsub
      endif

      !X boundaries
      do k = 0,nzp
        do j= 0,nyp
          vec2(0  ,j,k,:) = matmul(tensor(0  ,j,k,:,:),vec1(0  ,j,k,:))
          vec2(nxp,j,k,:) = matmul(tensor(nxp,j,k,:,:),vec1(nxp,j,k,:))
        enddo
      enddo

      !Y boundaries
      do k= 0,nzp
        do i=0,nxp
          vec2(i,0  ,k,:) = matmul(tensor(i,0  ,k,:,:),vec1(i,0  ,k,:))
          vec2(i,nyp,k,:) = matmul(tensor(i,nyp,k,:,:),vec1(i,nyp,k,:))
        enddo
      enddo

      !Z boundaries
      do j= 0,nyp
        do i= 0,nxp
          vec2(i,j,0  ,:) = matmul(tensor(i,j,0  ,:,:),vec1(i,j,0  ,:))
          vec2(i,j,nzp,:) = matmul(tensor(i,j,nzp,:,:),vec1(i,j,nzp,:))
        enddo
      enddo

c     End program

      end subroutine XformVector_BC

c     XformToCnv_stg_ijk
c     #################################################################
      function XformToCnv_stg_ijk(i,j,k,igrid,cov) result(cnv)

c     -----------------------------------------------------------------
c     Transforms a staggered covariant vector to another contravariant
c     representation at faces (i,j,k).
c
c     CURRENTLY WORKS ONLY IN 2D
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igrid
        real(8) :: cov(0:,0:,0:,:)
     .            ,cnv(3)

c     Local variables

        integer :: nx,ny,nz,ip,im,jp,jm,kp,km
        real(8) :: gsup11,gsup22
        real(8),pointer,dimension(:,:,:,:,:) :: gsup_v

c     Begin program

        ip = min(i+1,size(cov,1)-1)
        im = max(i-1,0)

        jp = min(j+1,size(cov,2)-1)
        jm = max(j-1,0)

        kp = min(k+1,size(cov,3)-1)
        km = max(k-1,0)

        gsup_v => gmetric%grid(igrid)%gsup_v

        gsup11 = 0.5*(gsup_v(i,jm,k,1,1)
     .               +gsup_v(i,j ,k,1,1))

        cnv(1) = gsup11*cov(i,j,k,1)
     .                +0.25*(gsup_v(i,j ,k,1,2)*(cov(i ,j ,k,2)
     .                                          +cov(ip,j ,k,2))
     .                      +gsup_v(i,jm,k,1,2)*(cov(i ,jm,k,2)
     .                                          +cov(ip,jm,k,2)))

        gsup22 = 0.5*(gsup_v(im,j,k,2,2)
     .               +gsup_v(i ,j,k,2,2))

        cnv(2) = gsup22*cov(i,j,k,2)
     .                +0.25*(gsup_v(i ,j,k,2,1)*(cov(i ,j ,k,1)
     .                                          +cov(i ,jp,k,1))
     .                      +gsup_v(im,j,k,2,1)*(cov(im,j ,k,1)
     .                                          +cov(im,jp,k,1)))

        cnv(3) = cov(i,j,k,3)

      end function XformToCnv_stg_ijk

c     XformToCov_stg_ijk
c     #################################################################
      function XformToCov_stg_ijk(i,j,k,igrid,cnv) result(cov)

c     -----------------------------------------------------------------
c     Transforms a staggered contravariant vector to a covariant
c     representation at faces for cell (i,j,k).
c
c     CURRENTLY WORKS ONLY IN 2D
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igrid
        real(8) :: cnv(0:,0:,0:,:)
     .            ,cov(3)

c     Local variables

        integer :: nx,ny,nz,ip,im,jp,jm,kp,km
        real(8) :: gsub11,gsub22
        real(8),pointer,dimension(:,:,:,:,:) :: gsub_v

c     Begin program

        ip = min(i+1,size(cnv,1)-1)
        im = max(i-1,0)

        jp = min(j+1,size(cnv,2)-1)
        jm = max(j-1,0)

        kp = min(k+1,size(cnv,3)-1)
        km = max(k-1,0)

        gsub_v => gmetric%grid(igrid)%gsub_v

        gsub11 = 0.5*(gsub_v(i,jm,k,1,1)
     .               +gsub_v(i,j ,k,1,1))

        cov(1) = gsub11*cnv(i,j,k,1)
     .                +0.25*(gsub_v(i,j ,k,1,2)*(cnv(i ,j ,k,2)
     .                                          +cnv(ip,j ,k,2))
     .                      +gsub_v(i,jm,k,1,2)*(cnv(i ,jm,k,2)
     .                                          +cnv(ip,jm,k,2)))

        gsub22 = 0.5*(gsub_v(im,j,k,2,2)
     .               +gsub_v(i ,j,k,2,2))

        cov(2) = gsub22*cnv(i,j,k,2)
     .                +0.25*(gsub_v(i ,j,k,2,1)*(cnv(i ,j ,k,1)
     .                                          +cnv(i ,jp,k,1))
     .                      +gsub_v(im,j,k,2,1)*(cnv(im,j ,k,1)
     .                                          +cnv(im,jp,k,1)))

        cov(3) = cnv(i,j,k,3)

      end function XformToCov_stg_ijk

c     Xform_From_VrtxCnv_To_CellCov_stg_ijk
c     #################################################################
      function Xform_From_VrtxCnv_To_CellCov_stg_ijk(i,j,k,igr,cnv_vrtx)
     .         result(cov_cell)

c     -----------------------------------------------------------------
c     Transforms a staggered contravariant vector to a covariant
c     representation at faces for cell (i,j,k).
c
c     CURRENTLY WORKS ONLY IN 2D
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igr
        real(8) :: cnv_vrtx(0:,0:,0:,:)
     .            ,cov_cell(3)

c     Local variables

        integer :: nx,ny,nz,ip,im,jp,jm,kp,km
        real(8) :: gsub12,gsub21
        real(8),pointer,dimension(:,:,:,:,:) :: gsub_v

c     Begin program

        ip = min(i+1,size(cnv_vrtx,1)-1)
        im = max(i-1,0)

        jp = min(j+1,size(cnv_vrtx,2)-1)
        jm = max(j-1,0)

        kp = min(k+1,size(cnv_vrtx,3)-1)
        km = max(k-1,0)

        gsub_v => gmetric%grid(igr)%gsub_v

        gsub12 = 0.5*(gsub_v(i,jm,k,1,2)
     .               +gsub_v(i,j ,k,1,2))

        cov_cell(1) = gsub12*cnv_vrtx(i,j,k,2)
     .          +0.25*(gsub_v(i,j ,k,1,1)*(cnv_vrtx(i ,j ,k,1)
     .                                    +cnv_vrtx(ip,j ,k,1))
     .                +gsub_v(i,jm,k,1,1)*(cnv_vrtx(i ,jm,k,1)
     .                                    +cnv_vrtx(ip,jm,k,1)))

        gsub21 = 0.5*(gsub_v(im,j,k,2,1)
     .               +gsub_v(i ,j,k,2,1))

        cov_cell(2) = gsub21*cnv_vrtx(i,j,k,1)
     .          +0.25*(gsub_v(i ,j,k,2,2)*(cnv_vrtx(i ,j ,k,2)
     .                                    +cnv_vrtx(i ,jp,k,2))
     .                +gsub_v(im,j,k,2,2)*(cnv_vrtx(im,j ,k,2)
     .                                    +cnv_vrtx(im,jp,k,2)))

        cov_cell(3) = 0.25*(cnv_vrtx(im,jm,k,3)
     .                     +cnv_vrtx(im,j ,k,3)
     .                     +cnv_vrtx(i ,j ,k,3)
     .                     +cnv_vrtx(i ,jm,k,3))

      end function Xform_From_VrtxCnv_To_CellCov_stg_ijk

c     XformToCnv_stg_mesh
c     #################################################################
      function XformToCnv_stg_mesh(igrid,cov) result(cnv)

c     -----------------------------------------------------------------
c     Transforms a staggered covariant vector to a contravariant
c     representation at faces everywhere on the mesh.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igrid
        real(8) :: cov(0:,0:,0:,:)
     .            ,cnv(0:size(cov,1)-1
     .                ,0:size(cov,2)-1
     .                ,0:size(cov,3)-1,3)

c     Local variables

        integer :: i,j,k,nx,ny,nz,ip,im,jp,jm,kp,km
        real(8) :: gsup11,gsup22
        real(8),pointer,dimension(:,:,:,:,:) :: gsup_v

c     Begin program

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        cnv = 0d0

        do k=0,nz
          do j=0,ny
            do i=0,nx
              cnv(i,j,k,:) = XformToCnv_stg_ijk(i,j,k,igrid,cov)
            enddo
          enddo
        enddo

      end function XformToCnv_stg_mesh

c     XformToCov_stg_mesh
c     #################################################################
      function XformToCov_stg_mesh(igrid,cnv,vertex) result(cov)

c     -----------------------------------------------------------------
c     Transforms a staggered contravariant vector (vertex-based if
c     vertex=.true.) to a cell-based covariant staggered representation
c     everywhere on the mesh.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igrid
        real(8),INTENT(IN) :: cnv(0:,0:,0:,:)
        real(8) :: cov(0:size(cnv,1)-1
     .                ,0:size(cnv,2)-1
     .                ,0:size(cnv,3)-1,3)
        logical, optional :: vertex

c     Local variables

        integer :: i,j,k,nx,ny,nz
        logical :: vrtx

c     Begin program

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        nx = size(cnv,1)-2
        ny = size(cnv,2)-2
        nz = size(cnv,3)-2

        cov = 0d0

        if (vrtx) then
          do k=0,nz
            do j=0,ny
              do i=0,nx
                cov(i,j,k,:) =
     .            Xform_From_VrtxCnv_To_CellCov_stg_ijk(i,j,k,igrid,cnv)
              enddo
            enddo
          enddo
        else
          do k=0,nz
            do j=0,ny
              do i=0,nx
                cov(i,j,k,:) = XformToCov_stg_ijk(i,j,k,igrid,cnv)
              enddo
            enddo
          enddo
        endif

      end function XformToCov_stg_mesh

ccc     XformVector_stg_BC
ccc     ###################################################################
cc      subroutine XformVector_stg_BC(igrid,vec1,vec2,tocnv)
cc
cc      implicit none
cc
ccc     -------------------------------------------------------------------
ccc     Routine to transfer ghost cell values form vec1 to vec2
ccc     -------------------------------------------------------------------
cc
ccc     Input variables
cc
cc      integer :: igrid
cc      real(8) :: vec1(0:,0:,0:,:)
cc      real(8) :: vec2(0:,0:,0:,:)
cc      logical :: tocnv
cc
ccc     Local variables
cc
cc      integer :: nxp,nyp,nzp,i,j,k
cc
ccc     Begin program
cc
cc      nxp = size(vec1,1)-1
cc      nyp = size(vec1,2)-1
cc      nzp = size(vec1,3)-1
cc
cc      if (tocnv) then
cc
cc        !X boundaries
cc        do k= 0,nzp
cc          do j= 0,nyp
cc            vec2(0  ,j,k,:) = XformToCnv_stg_ijk(0  ,j,k,igrid,vec1)
cc            vec2(nxp,j,k,:) = XformToCnv_stg_ijk(nxp,j,k,igrid,vec1)
cc          enddo
cc        enddo
cc
cc        !Y boundaries
cc        do k= 0,nzp
cc          do i=0,nxp
cc            vec2(i,0  ,k,:) = XformToCnv_stg_ijk(i,0  ,k,igrid,vec1)
cc            vec2(i,nyp,k,:) = XformToCnv_stg_ijk(i,nyp,k,igrid,vec1)
cc          enddo
cc        enddo
cc
cc        !Z boundaries
cc        do j= 0,nyp
cc          do i= 0,nxp
cc            vec2(i,j,0  ,:) = XformToCnv_stg_ijk(i,j,0  ,igrid,vec1)
cc            vec2(i,j,nzp,:) = XformToCnv_stg_ijk(i,j,nzp,igrid,vec1)
cc          enddo
cc        enddo
cc
cc      else
cc        call pstop("XformVector_stg_BC","Option not implemented")
cc      endif
cc
ccc     End program
cc
cc      end subroutine XformVector_stg_BC

c     XformVector_stg_mesh
c     #################################################################
      function XformVector_stg_mesh(igrid,vec1,tocnv) result(vec2)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at cell centers on all mesh.
c     -----------------------------------------------------------------

      implicit none

c     Input variables

      integer    :: igrid
      real(8)    :: vec1(0:,0:,0:,:)
      real(8)    :: vec2(0:size(vec1,1)-1
     .                  ,0:size(vec1,2)-1
     .                  ,0:size(vec1,3)-1
     .                  ,size(vec1,4))

      logical    :: tocnv

c     Local variables

c     Begin program

      if (tocnv) then
        vec1 = XformToCnv_stg_mesh(igrid,vec1)
      else
        vec1 = XformToCov_stg_mesh(igrid,vec1)
      endif

c     End program

      end function XformVector_stg_mesh

c     Xfer_stg_to_cell_mesh
c     #################################################################
      function Xfer_stg_to_cell_mesh(stg,vrtx) result(cell)

c     -----------------------------------------------------------------
c     Averages staggered vector to cell centers. The variable vrtx
c     determines whether vector is cell-based (components orthogonal to
c     cell-center volumes) or vertex-based (components orthogonal to
c     vertex-centered volumes)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: igrid
        real(8) :: stg(0:,0:,0:,:)
     .            ,cell(0:size(stg,1)-1
     .                 ,0:size(stg,2)-1
     .                 ,0:size(stg,3)-1
     .                 ,  size(stg,4)  )

        logical,optional :: vrtx

c     Local variables

        integer :: i,j,k,nx,ny,nz,neq
        logical :: vertex

c     Begin program

        if (PRESENT(vrtx)) then
          vertex = vrtx
        else
          vertex = .false.
        endif

        nx = size(stg,1)-2
        ny = size(stg,2)-2
        nz = size(stg,3)-2
        neq= size(stg,4)

        cell = 0d0

        if (.not.vertex) then
          cell(1:nx,:,:,1) = 0.5*(stg(0:nx-1,:,:,1)+stg(1:nx,:,:,1))
          cell(:,1:ny,:,2) = 0.5*(stg(:,0:ny-1,:,2)+stg(:,1:ny,:,2))
          cell(:,:,1:nz,3) = 0.5*(stg(:,:,0:nz-1,3)+stg(:,:,1:nz,3))
        else
          cell(1:nx,1:ny,1:nz,1) = 0.25*(stg(1:nx,0:ny-1,0:nz-1,1)
     .                                  +stg(1:nx,0:ny-1,1:nz  ,1)
     .                                  +stg(1:nx,1:ny  ,1:nz  ,1)
     .                                  +stg(1:nx,1:ny  ,0:nz-1,1))
          cell(1:nx,1:ny,1:nz,2) = 0.25*(stg(0:nx-1,1:ny,0:nz-1,2)
     .                                  +stg(0:nx-1,1:ny,1:nz  ,2)
     .                                  +stg(1:nx  ,1:ny,1:nz  ,2)
     .                                  +stg(1:nx  ,1:ny,0:nz-1,2))
          cell(1:nx,1:ny,1:nz,3) = 0.25*(stg(0:nx-1,0:ny-1,1:nz,3)
     .                                  +stg(0:nx-1,1:ny  ,1:nz,3)
     .                                  +stg(1:nx  ,1:ny  ,1:nz,3)
     .                                  +stg(1:nx  ,0:ny-1,1:nz,3))
        endif

      end function Xfer_stg_to_cell_mesh

c     vectorNorm_ijk
c     ################################################################
      function vectorNorm_ijk(i,j,k,igrid,vec,covar) result(norm)

c     ---------------------------------------------------------------
c     Finds L2-norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igrid
      real(8)    :: vec(3),norm
      logical    :: covar

c     Local variables

      real(8)    :: tensor(3,3)

c     Begin program

      if (covar) then
        tensor = gmetric%grid(igrid)%gsup(i,j,k,:,:)
      else
        tensor = gmetric%grid(igrid)%gsub(i,j,k,:,:)
      endif

      norm = dot_product(vec,matmul(tensor,vec))
     .      /gmetric%grid(igrid)%jac(i,j,k)

c     End 

      end function vectorNorm_ijk

c     vectorNorm_mesh
c     ################################################################
      function vectorNorm_mesh(igrid,vec,covar) result(norm)

c     ---------------------------------------------------------------
c     Finds L2-norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid
      real(8) :: vec(0:,0:,0:,:)
     .          ,norm(0:size(vec,1)-1,0:size(vec,2)-1,0:size(vec,3)-1)
      logical :: covar

c     Local variables

      integer :: i,j,k,nxl,nyl,nzl

c     Begin program

      nxl = grid_params%nxv(igrid)
      nyl = grid_params%nyv(igrid)
      nzl = grid_params%nzv(igrid)

      do k=0,nzl+1
        do j=0,nyl+1
          do i=0,nxl+1
            norm(i,j,k) = vectorNorm_ijk(i,j,k,igrid,vec(i,j,k,:),covar)
          enddo
        enddo
      enddo

c     End 

      end function vectorNorm_mesh

c     scalarProduct_ijk
c     ################################################################
      function scalarProduct_ijk(i,j,k,igrid,cov,cnv) result (dot)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igrid
      real(8) :: dot,cnv(3),cov(3)

c     Local variables

c     Begin program

      dot = dot_product(cov,cnv)/gmetric%grid(igrid)%jac(i,j,k)

c     End 

      end function scalarProduct_ijk

c     scalarProduct_mesh
c     ################################################################
      function scalarProduct_mesh(igrid,vec1,vec2) result (sp)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid
      real(8) :: vec1(:,:,:,:),vec2(:,:,:,:)
      real(8) :: sp(size(vec1,1),size(vec1,2),size(vec1,3))

c     Local variables

c     Begin program

      sp = (vec1(:,:,:,1)*vec2(:,:,:,1)
     .     +vec1(:,:,:,2)*vec2(:,:,:,2)
     .     +vec1(:,:,:,3)*vec2(:,:,:,3))/gmetric%grid(igrid)%jac

c     End 

      end function scalarProduct_mesh

c     tensorScalarProduct
c     ################################################################
      function tensorScalarProduct(i,j,k,igrid,tnsr1,tnsr2)
     .         result (dot)

c     ---------------------------------------------------------------
c     Finds tensor scalar product of two tensors, one covariant
c     and the other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igrid
      real(8)    :: dot,tnsr1(3,3),tnsr2(3,3)

c     Local variables

c     Begin program

      dot = sum(tnsr1*tnsr2)/gmetric%grid(igrid)%jac(i,j,k)**2

c     End 

      end function tensorScalarProduct

c     crossProduct_ijk
c     ################################################################
      function crossProduct_ijk(i,j,k,igrid,vec1,vec2,cnv_in)
     .         result (vec3)

c     ---------------------------------------------------------------
c     Finds cross product of two vectors, either both covariant
c     (cnv_in=.false.) or both contravariant (cnv_in=.true.) Result
c     is covariant if cnv_in=.true., and contravariant otherwise.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igrid
      real(8)    :: vec1(3),vec2(3),vec3(3)
      logical    :: cnv_in

c     Local variables

c     Begin program

      if (cnv_in) then  !Result is covariant
        vec3(1) = (vec1(2)*vec2(3)
     .            -vec1(3)*vec2(2))
     .            /gmetric%grid(igrid)%jac(i,j,k)
        vec3(2) = (vec1(3)*vec2(1)
     .            -vec1(1)*vec2(3))
     .            /gmetric%grid(igrid)%jac(i,j,k)
        vec3(3) = (vec1(1)*vec2(2)
     .            -vec1(2)*vec2(1))
     .            /gmetric%grid(igrid)%jac(i,j,k)
      else           !Result is contravariant
        vec3(1) = (vec1(2)*vec2(3)
     .            -vec1(3)*vec2(2))
        vec3(2) = (vec1(3)*vec2(1)
     .            -vec1(1)*vec2(3))
        vec3(3) = (vec1(1)*vec2(2)
     .            -vec1(2)*vec2(1))
      endif

c     End 

      end function crossProduct_ijk

c     crossProduct_mesh
c     ################################################################
      function crossProduct_mesh(igrid,vec1,vec2,cnv_in)
     .         result (vec3)

c     ---------------------------------------------------------------
c     Finds cross product of two vectors, either both covariant
c     (cnv_in=.false.) or both contravariant (cnv_in=.true.) Result
c     is covariant if cnv_in=.true., and contravariant otherwise.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: igrid
      real(8)    :: vec1(:,:,:,:),vec2(:,:,:,:)
      real(8)    :: vec3(size(vec1,1),size(vec1,2)
     .                  ,size(vec1,3),size(vec1,4))
      logical    :: cnv_in

c     Local variables

c     Begin program

      if (cnv_in) then  !Result is covariant
        vec3(:,:,:,1) = (vec1(:,:,:,2)*vec2(:,:,:,3)
     .                  -vec1(:,:,:,3)*vec2(:,:,:,2))
     .                  /gmetric%grid(igrid)%jac
        vec3(:,:,:,2) = (vec1(:,:,:,3)*vec2(:,:,:,1)
     .                  -vec1(:,:,:,1)*vec2(:,:,:,3))
     .                  /gmetric%grid(igrid)%jac
        vec3(:,:,:,3) = (vec1(:,:,:,1)*vec2(:,:,:,2)
     .                  -vec1(:,:,:,2)*vec2(:,:,:,1))
     .                  /gmetric%grid(igrid)%jac
      else           !Result is contravariant
        vec3(:,:,:,1) = (vec1(:,:,:,2)*vec2(:,:,:,3)
     .                  -vec1(:,:,:,3)*vec2(:,:,:,2))
        vec3(:,:,:,2) = (vec1(:,:,:,3)*vec2(:,:,:,1)
     .                  -vec1(:,:,:,1)*vec2(:,:,:,3))
        vec3(:,:,:,3) = (vec1(:,:,:,1)*vec2(:,:,:,2)
     .                  -vec1(:,:,:,2)*vec2(:,:,:,1))
      endif

c     End 

      end function crossProduct_mesh

      end module grid_operations

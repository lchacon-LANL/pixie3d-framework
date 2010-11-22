
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

        INTERFACE XformVector
          module procedure XformVector_mesh
        end INTERFACE

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

c     transformVectorToCartesian
c     #################################################################
      subroutine transformVectorToCartesian(i,j,k,igx,igy,igz
     .                                     ,c1,c2,c3,covariant
     .                                     ,cx,cy,cz)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: covariant
        

c     Local variables

        integer    :: ic
        real(8)    :: T_to_car(3,3),vec(3)

c     Begin program

        vec = (/ c1,c2,c3 /)

        if (covariant) then
          cx = dot_product(vec,gmetric%grid(igx)%cov(i,j,k,:,1))
          cy = dot_product(vec,gmetric%grid(igx)%cov(i,j,k,:,2))
          cz = dot_product(vec,gmetric%grid(igx)%cov(i,j,k,:,3))
        else
          cx = dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,:,1))
          cy = dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,:,2))
          cz = dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,:,3))
        endif

      end subroutine transformVectorToCartesian

c     transformVectorToCurvilinear
c     #################################################################
      subroutine transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .                                       ,cx,cy,cz,covariant
     .                                       ,c1,c2,c3)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector (cx,cy,cz) to curvilinear (c1,c2,c3)
c     (covariant if covariant=.true., contravariant otherwise)
c     at grid coordinates (i,j,k). 
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: covariant

c     Local variables

        integer    :: ic
        real(8)    :: T_to_curv(3,3),vec(3),jac

c     Begin program

        jac = gmetric%grid(igx)%jac(i,j,k)

        vec = (/ cx,cy,cz /)

        if (covariant) then

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

        if (coords == 'car') then

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
        cnvariant = (.not. cartesian .and. .not. covariant)

        select case(out_repr)
        case('car')

          if (cartesian) return

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCartesian
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

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
     .                 ,covariant
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
     .                 ,covariant,0)

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
     .                 ,covariant
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
     .                 ,covariant,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        end select

      end subroutine transformVector

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

c     vectorNorm_ijk
c     ################################################################
      function vectorNorm_ijk(i,j,k,igx,igy,igz,ax,ay,az,covar)
     .         result(norm)

c     ---------------------------------------------------------------
c     Finds L2-norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igx,igy,igz
      real(8)    :: ax,ay,az,norm
      logical    :: covar

c     Local variables

      real(8)    :: tensor(3,3),cnv(3),cov(3),jac

c     Begin program

      jac = gmetric%grid(igx)%jac(i,j,k)

      if (covar) then
        tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
        cov = (/ ax,ay,az /)
        cnv = matmul(tensor,cov)
      else
        tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
        cnv = (/ ax,ay,az /)
        cov = matmul(tensor,cnv)
      endif

      norm = dot_product(cov,cnv)/jac

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
            norm(i,j,k) = vectorNorm_ijk(i,j,k,igrid,igrid,igrid
     .                                  ,vec(i,j,k,1)
     .                                  ,vec(i,j,k,2)
     .                                  ,vec(i,j,k,3),covar)
          enddo
        enddo
      enddo

c     End 

      end function vectorNorm_mesh

c     scalarProduct_ijk
c     ################################################################
      function scalarProduct_ijk(i,j,k,igx,igy,igz,cov1,cov2,cov3
     .                       ,cnv1,cnv2,cnv3) result (dot)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igx,igy,igz
      real(8)    :: dot,cov1,cov2,cov3,cnv1,cnv2,cnv3

c     Local variables

      real(8)    :: cnv(3),cov(3)

c     Begin program

      cnv = (/ cnv1,cnv2,cnv3 /)
      cov = (/ cov1,cov2,cov3 /)

      dot = dot_product(cov,cnv)/gmetric%grid(igx)%jac(i,j,k)

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
cc      real(8) :: vec1(0:,0:,0:,:),vec2(0:,0:,0:,:)
cc      real(8) :: sp(0:size(vec1,1)-1,0:size(vec1,2)-1,0:size(vec1,3)-1)

c     Local variables

cc      integer :: i,j,k,nxl,nyl,nzl

c     Begin program

      sp = (vec1(:,:,:,1)*vec2(:,:,:,1)
     .     +vec1(:,:,:,2)*vec2(:,:,:,2)
     .     +vec1(:,:,:,3)*vec2(:,:,:,3))/gmetric%grid(igrid)%jac

cc      nxl = grid_params%nxv(igrid)
cc      nyl = grid_params%nyv(igrid)
cc      nzl = grid_params%nzv(igrid)
cc
cc      do k=0,nzl+1
cc        do j=0,nyl+1
cc          do i=0,nxl+1
cc            sp(i,j,k) = scalarProduct_ijk(i,j,k,igrid,igrid,igrid
cc     .                 ,vec1(i,j,k,1),vec1(i,j,k,2),vec1(i,j,k,3)
cc     .                 ,vec2(i,j,k,1),vec2(i,j,k,2),vec2(i,j,k,3))
cc          enddo
cc        enddo
cc      enddo

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

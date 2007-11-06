
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

      contains

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        integer    :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point = (i+grid_params%ilo(igx)-1==1).and.bcSP()

      end function isSP

c     isSYM
c     #################################################################
      function isSYM(i,igr,dim,loc)

        integer    :: i,igr,dim,loc
        logical    :: isSYM

        integer    :: ibc,ilog,itst

        ibc = (1+loc)+2*(dim-1)

        isSYM=.false.

        select case(dim)
        case(1)
          ilog = grid_params%ilo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nxgl(igr)
          end select
        case(2)
          ilog = grid_params%jlo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nygl(igr)
          end select
        case(3)
          ilog = grid_params%klo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nzgl(igr)
          end select
        end select

        if (i+ilog-1 == itst .and. bcond(ibc) == SYM)  isSYM=.true.

      end function isSYM

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

cc        if (tocnv) then
cc
cc          select case(he)
cc          case(1)
cc            ip = i+1
cccc            if (isSP(ip,j,k,igx,igy,igz)) then !We are right @ SP ghost cell
cccc              tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
cccc            elseif (isSP(i,j,k,igx,igy,igz)) then !We are right SP
cccc              tensor = (gmetric%grid(igx)%gsup(ip,j,k,:,:)
cccc     .                 *gmetric%grid(igx)%jac (ip,j,k)
cccc     .                 +gmetric%grid(igx)%gsup(i ,j,k,:,:)
cccc     .                 *gmetric%grid(igx)%jac (i ,j,k))
cccc     .                 /(gmetric%grid(igx)%jac(ip,j,k)
cccc     .                  +gmetric%grid(igx)%jac(i ,j,k))
cccc            else
cc              tensor = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
cc     .                     +gmetric%grid(igx)%gsup(i ,j,k,:,:))
cccc            endif
cc          case(2)
cc            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j+1,k,:,:)
cc     .                   +gmetric%grid(igx)%gsup(i,j  ,k,:,:))
cc          case(3)
cc            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j,k+1,:,:)
cc     .                   +gmetric%grid(igx)%gsup(i,j,k  ,:,:))
cc          case default
cc            tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
cc          end select
cc
cc          cnv1 = tensor(1,1)*cov1+tensor(1,2)*cov2+tensor(1,3)*cov3
cc          cnv2 = tensor(2,1)*cov1+tensor(2,2)*cov2+tensor(2,3)*cov3
cc          cnv3 = tensor(3,1)*cov1+tensor(3,2)*cov2+tensor(3,3)*cov3
cc        else
cc          select case(he)
cc          case(1)
cc            ip = i+1
cccc            if (isSP(ip,j,k,igx,igy,igz)) then !We are right @ SP ghost cell
cccc              tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
cccc            elseif (isSP(i,j,k,igx,igy,igz)) then !We are right SP
cccc              tensor = (gmetric%grid(igx)%gsub(ip,j,k,:,:)
cccc     .                 *gmetric%grid(igx)%jac (ip,j,k)
cccc     .                 +gmetric%grid(igx)%gsub(i ,j,k,:,:)
cccc     .                 *gmetric%grid(igx)%jac (i ,j,k))
cccc     .                 /(gmetric%grid(igx)%jac(ip,j,k)
cccc     .                  +gmetric%grid(igx)%jac(i ,j,k))
cccc            else
cc              tensor = 0.5*(gmetric%grid(igx)%gsub(ip,j,k,:,:)
cc     .                     +gmetric%grid(igx)%gsub(i ,j,k,:,:))
cccc            endif
cc          case(2)
cc            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j+1,k,:,:)
cc     .                   +gmetric%grid(igx)%gsub(i,j  ,k,:,:))
cc          case(3)
cc            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j,k+1,:,:)
cc     .                   +gmetric%grid(igx)%gsub(i,j,k  ,:,:))
cc          case default
cc            tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
cc          end select
cc
cc          cov1 = tensor(1,1)*cnv1+tensor(1,2)*cnv2+tensor(1,3)*cnv3
cc          cov2 = tensor(2,1)*cnv1+tensor(2,2)*cnv2+tensor(2,3)*cnv3
cc          cov3 = tensor(3,1)*cnv1+tensor(3,2)*cnv2+tensor(3,3)*cnv3
cc        endif

      end subroutine transformFromCurvToCurv

c     vectorNorm
c     ################################################################
      function vectorNorm(i,j,k,igx,igy,igz,ax,ay,az,covar)

c     ---------------------------------------------------------------
c     Finds L2-norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: i,j,k,igx,igy,igz
      real(8)    :: ax,ay,az,vectorNorm
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

      vectorNorm = dot_product(cov,cnv)/jac

c     End 

      end function vectorNorm

c     scalarProduct
c     ################################################################
      function scalarProduct(i,j,k,igx,igy,igz,cov1,cov2,cov3
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

      end function scalarProduct

      end module grid_operations

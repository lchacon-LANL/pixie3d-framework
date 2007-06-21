
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
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv,half_elem)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at grid coordinates (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        integer   ,optional,intent(IN) :: half_elem

        real(8)    :: cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical    :: tocnv

c     Local variables

        integer    :: he,ip
        real(8)    :: cov(3),cnv(3),tensor(3,3)

c     Begin program

        if (PRESENT(half_elem)) then
          he = half_elem
        else
          he = 0
        endif

        if (tocnv) then

          select case(he)
          case(1)
            ip = i+1
cc            if (isSP(ip,j,k,igx,igy,igz)) then !We are right @ SP ghost cell
cc              tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
cc            elseif (isSP(i,j,k,igx,igy,igz)) then !We are right SP
cc              tensor = (gmetric%grid(igx)%gsup(ip,j,k,:,:)
cc     .                 *gmetric%grid(igx)%jac (ip,j,k)
cc     .                 +gmetric%grid(igx)%gsup(i ,j,k,:,:)
cc     .                 *gmetric%grid(igx)%jac (i ,j,k))
cc     .                 /(gmetric%grid(igx)%jac(ip,j,k)
cc     .                  +gmetric%grid(igx)%jac(i ,j,k))
cc            else
              tensor = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .                     +gmetric%grid(igx)%gsup(i ,j,k,:,:))
cc            endif
          case(2)
            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j+1,k,:,:)
     .                   +gmetric%grid(igx)%gsup(i,j  ,k,:,:))
          case(3)
            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j,k+1,:,:)
     .                   +gmetric%grid(igx)%gsup(i,j,k  ,:,:))
          case default
            tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
          end select

          cnv1 = tensor(1,1)*cov1+tensor(1,2)*cov2+tensor(1,3)*cov3
          cnv2 = tensor(2,1)*cov1+tensor(2,2)*cov2+tensor(2,3)*cov3
          cnv3 = tensor(3,1)*cov1+tensor(3,2)*cov2+tensor(3,3)*cov3
        else
          select case(he)
          case(1)
            ip = i+1
cc            if (isSP(ip,j,k,igx,igy,igz)) then !We are right @ SP ghost cell
cc              tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
cc            elseif (isSP(i,j,k,igx,igy,igz)) then !We are right SP
cc              tensor = (gmetric%grid(igx)%gsub(ip,j,k,:,:)
cc     .                 *gmetric%grid(igx)%jac (ip,j,k)
cc     .                 +gmetric%grid(igx)%gsub(i ,j,k,:,:)
cc     .                 *gmetric%grid(igx)%jac (i ,j,k))
cc     .                 /(gmetric%grid(igx)%jac(ip,j,k)
cc     .                  +gmetric%grid(igx)%jac(i ,j,k))
cc            else
              tensor = 0.5*(gmetric%grid(igx)%gsub(ip,j,k,:,:)
     .                     +gmetric%grid(igx)%gsub(i ,j,k,:,:))
cc            endif
          case(2)
            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j+1,k,:,:)
     .                   +gmetric%grid(igx)%gsub(i,j  ,k,:,:))
          case(3)
            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j,k+1,:,:)
     .                   +gmetric%grid(igx)%gsub(i,j,k  ,:,:))
          case default
            tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
          end select

          cov1 = tensor(1,1)*cnv1+tensor(1,2)*cnv2+tensor(1,3)*cnv3
          cov2 = tensor(2,1)*cnv1+tensor(2,2)*cnv2+tensor(2,3)*cnv3
          cov3 = tensor(3,1)*cnv1+tensor(3,2)*cnv2+tensor(3,3)*cnv3
        endif

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

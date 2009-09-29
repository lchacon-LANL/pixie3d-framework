      module grid_mgmt

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
      implicit none

      contains

      subroutine setup_var_array(fgrid,i,c_array_ptr,xs,ys,zs,xe,ye,ze)
        implicit none
        type(var_array) :: fgrid
        integer    :: i,xs,ys,zs,xe,ye,ze
        REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze)
        fgrid%array_var(i+1)%array => c_array_ptr
      end subroutine setup_var_array

      subroutine setup_aux_var(fgrid,i,c_array_ptr,xs,ys,zs,xe,ye,ze)
        implicit none
        type(aux_array) :: fgrid
        integer    :: i,xs,ys,zs,xe,ye,ze
        REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze)
        fgrid%var_list(i+1)%array => c_array_ptr
      end subroutine setup_aux_var

      subroutine setup_aux_vec(fgrid,i,c_array_ptr,xs,ys,zs,xe,ye,ze)
        implicit none
        type(aux_array) :: fgrid
        integer    :: i,xs,ys,zs,xe,ye,ze
        REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze,1:3)
        fgrid%vec_list(i+1)%vec => c_array_ptr
      end subroutine setup_aux_vec

      subroutine setup_patch(fgrid,u_0,u_n,aux,gparams)
        implicit none
        type(patch), pointer :: fgrid
        type(var_array), TARGET :: u_0, u_n
        type(aux_array), TARGET :: aux
        type(grid_mg_def), TARGET :: gparams
        fgrid%u_0 => u_0
        fgrid%u_n => u_n
        fgrid%aux => aux
        fgrid%gparams => gparams
      end subroutine setup_patch

      end module grid_mgmt

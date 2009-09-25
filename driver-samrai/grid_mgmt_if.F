c f_create_var_array
c######################################################################
      subroutine f_create_var_array(patch_obj,nvar)

      use variable_setup, only: var_array
      implicit none

      integer*8 :: patch_obj
      integer :: nvar
      type(var_array), pointer :: grid_data
      allocate(grid_data)
      grid_data%nvar = nvar
      allocate(grid_data%array_var(1:nvar))
      call pflow_loc(patch_obj, grid_data)
      
      end subroutine f_create_var_array


c f_create_aux_array
c######################################################################
      subroutine f_create_aux_array(patch_obj,nvar,nvec)

      use variable_setup, only: aux_array
 
      integer*8 :: patch_obj
      integer :: nvar, nvec
      type(aux_array), pointer :: grid_data
      allocate(grid_data)
      grid_data%nvar = nvar
      grid_data%nvec = nvec
      allocate(grid_data%var_list(1:nvar))
      allocate(grid_data%vec_list(1:nvec))
      call pflow_loc(patch_obj, grid_data)

      end subroutine f_create_aux_array


c f_create_grid_mg_def
c######################################################################
      subroutine f_create_grid_mg_def(patch_obj,nglx,ngly,nglz
     .      ,ilog,jlog,klog,ihig,jhig,khig,gcw)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
c      use grid_mgmt, only: setup_patch
 
      integer*8 :: patch_obj
      integer :: nglx, ngly, nglz
      integer :: ilog, jlog, klog, ihig, jhig, khig, gcw
      type(grid_mg_def), pointer :: grid_data

      allocate(grid_data)

      grid_data%ngrdx = 1
      grid_data%ngrdy = 1
      grid_data%ngrdz = 1
      grid_data%ngrid = 1
      grid_data%nglx  = nglx
      grid_data%ngly  = ngly
      grid_data%nglz  = nglz
      grid_data%nlx   = ihig-ilog+1
      grid_data%nly   = jhig-jlog+1
      grid_data%nlz   = khig-klog+1
      grid_data%ilog  = ilog
      grid_data%ihig  = ihig
      grid_data%jlog  = jlog
      grid_data%jhig  = jhig
      grid_data%klog  = klog
      grid_data%khig  = khig
      grid_data%gcw   = 1
cc#if defined(samrai)
cc      nullify(grid_data%zeros)
cc      nullify(grid_data%ones)
cc      nullify(grid_data%vzeros)
cc#endif

      call pflow_loc(patch_obj, grid_data)

      end subroutine f_create_grid_mg_def


c f_create_patch_data
c######################################################################
      subroutine f_create_patch_data(patch_obj,u_0,u_n,aux,gparams)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
c      use grid_mgmt, only: setup_patch
 
      integer*8 :: patch_obj
      type(var_array),INTENT(INOUT) :: u_0, u_n
      type(aux_array),INTENT(INOUT) :: aux
      type(grid_mg_def),INTENT(INOUT) :: gparams
      type(patch), pointer :: grid_data

      INTERFACE
         subroutine setup_patch(fgrid,u_0,u_n,aux,gparams)
            use grid_def_st, only: grid_mg_def
            use variable_setup, only: var_array, aux_array, patch
            implicit none
            type(patch), pointer :: fgrid
            type(var_array), TARGET :: u_0, u_n
            type(aux_array), TARGET :: aux
            type(grid_mg_def), TARGET :: gparams
         end subroutine setup_patch
      END INTERFACE

      allocate(grid_data)
      nullify(grid_data%u_np)
      nullify(grid_data%u_n)
      nullify(grid_data%u_nm)
      nullify(grid_data%u_0)
      nullify(grid_data%u_graph)
      nullify(grid_data%gparams)
      nullify(grid_data%aux)
      nullify(grid_data%bc_seq)

      call setup_patch(grid_data,u_0,u_n,aux,gparams)
      call pflow_loc(patch_obj, grid_data)

      end subroutine f_create_patch_data


c fill_var_arrayf
c######################################################################
      subroutine fill_var_array(fgrid,i,c_array,xs,ys,zs,xe,ye,ze)

      use variable_setup
c      use grid_mgmt
      implicit none

      type(var_array),INTENT(INOUT) :: fgrid
      integer,INTENT(IN)       :: i,xs,ys,zs,xe,ye,ze
      REAL*8, INTENT(IN)       :: c_array(xs:xe,ys:ye,zs:ze)

      call setup_var_array(fgrid,i,c_array,xs,ys,zs,xe,ye,ze)

      end subroutine fill_var_array


c fill_aux_array_var
c######################################################################
      subroutine fill_aux_array_var(fgrid,i,c_array,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
c      use grid_mgmt

      type(aux_array),INTENT(INOUT) :: fgrid
      integer,INTENT(IN)       :: i,xs,ys,zs,xe,ye,ze
      REAL*8, INTENT(IN)       :: c_array(xs:xe,ys:ye,zs:ze)

      call setup_aux_var(fgrid,i,c_array,xs,ys,zs,xe,ye,ze)

      end subroutine fill_aux_array_var


c ill_aux_array_vec
c######################################################################
      subroutine fill_aux_array_vec(fgrid,i,c_array,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
c      use grid_mgmt

      type(aux_array),INTENT(INOUT) :: fgrid
      integer,INTENT(IN)       :: i,xs,ys,zs,xe,ye,ze
      REAL*8, INTENT(IN)       :: c_array(xs:xe,ys:ye,zs:ze,1:3)

      call setup_aux_vec(fgrid,i,c_array,xs,ys,zs,xe,ye,ze)

      end subroutine fill_aux_array_vec


c f_delete_patch_data
c######################################################################
      subroutine f_delete_patch_data(fgrid)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
      use parameters

      type(patch), target :: fgrid
      type(patch), pointer :: grid_data

      grid_data => fgrid;
c     deallocate(grid_data%u_0%array_var(1:grid_data%u_0%nvar)) 
c     deallocate(grid_data%u_n%array_var(1:grid_data%u_n%nvar)) 
      deallocate(grid_data%u_0)
      deallocate(grid_data%u_n) 
c     deallocate(grid_data%aux%var_list(1:grid_data%aux%nvar)) 
c     deallocate(grid_data%aux%vec_list(1:grid_data%aux%nvec)) 
      deallocate(grid_data%aux)
      deallocate(grid_data%gparams)
      deallocate(grid_data)

      end subroutine f_delete_patch_data


c f_delete_var_array
c######################################################################
      subroutine f_delete_var_array(fgrid)

      use variable_setup, only: var_array

      type(var_array), target :: fgrid
      type(var_array), pointer :: grid_data

      grid_data => fgrid;
c     deallocate(grid_data%array_var(1:grid_data%nvar)) 
      deallocate(grid_data)

      end subroutine f_delete_var_array


c setup_var_array
c######################################################################
      subroutine setup_var_array(fgrid,i,c_array_ptr,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: var_array
      implicit none

      type(var_array) :: fgrid
      integer    :: i,xs,ys,zs,xe,ye,ze
      REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze)

      fgrid%array_var(i+1)%array => c_array_ptr

      end subroutine setup_var_array


c setup_aux_var
c######################################################################
      subroutine setup_aux_var(fgrid,i,c_array_ptr,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
      implicit none

      type(aux_array) :: fgrid
      integer    :: i,xs,ys,zs,xe,ye,ze
      REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze)
      fgrid%var_list(i+1)%array => c_array_ptr

      end subroutine setup_aux_var


c setup_aux_vec
c######################################################################
      subroutine setup_aux_vec(fgrid,i,c_array_ptr,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
      implicit none

      type(aux_array) :: fgrid
      integer    :: i,xs,ys,zs,xe,ye,ze
      REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze,1:3)

      fgrid%vec_list(i+1)%vec => c_array_ptr

      end subroutine setup_aux_vec


c setup_patch
c######################################################################
      subroutine setup_patch(fgrid,u_0,u_n,aux,gparams)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
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



c f_create_var_array
c######################################################################
      subroutine f_create_var_array(c_var_array_ptr,nvar)

      use variable_setup, only: var_array
      implicit none

      integer*8 :: c_var_array_ptr
      integer :: nvar
      integer :: i
      type(var_array), pointer :: var_array_data
      allocate(var_array_data)
      var_array_data%nvar = nvar
      allocate(var_array_data%array_var(1:nvar))
      do i=1,nvar
           nullify(var_array_data%array_var(i)%array)         
      enddo

      call setup_f90c_ptrlinks(c_var_array_ptr, var_array_data)
      
      end subroutine f_create_var_array

c f_create_aux_array
c######################################################################
      subroutine f_create_aux_array(c_aux_array_ptr,nvar,nvec)

      use variable_setup, only: aux_array
 
      integer*8 :: c_aux_array_ptr
      integer :: nvar, nvec
      type(aux_array), pointer :: aux_array_data
      allocate(aux_array_data)
      aux_array_data%nvar = nvar
      aux_array_data%nvec = nvec
      allocate(aux_array_data%var_list(1:nvar))
      allocate(aux_array_data%vec_list(1:nvec))
      call setup_f90c_ptrlinks(c_aux_array_ptr, aux_array_data)

      end subroutine f_create_aux_array


c f_create_grid_mg_def
c######################################################################
      subroutine f_create_grid_mg_def(c_grid_mg_def_ptr, cmin, cmax,
     & nglx, ngly,nglz,ilog,jlog,klog,ihig,jhig,khig,gcw)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
c      use grid_mgmt, only: setup_patch
 
      integer*8 :: c_grid_mg_def_ptr
      double precision ::cmin(3), cmax(3)
      integer :: nglx, ngly, nglz
      integer :: ilog, jlog, klog, ihig, jhig, khig, gcw
      type(grid_mg_def), pointer :: grid_data

      allocate(grid_data)
      nullify(grid_data%ilo)
      nullify(grid_data%jlo)
      nullify(grid_data%klo)
      nullify(grid_data%ihi)
      nullify(grid_data%jhi)
      nullify(grid_data%khi)

      nullify(grid_data%xg)
      nullify(grid_data%yg)
      nullify(grid_data%zg)
      nullify(grid_data%xx)
      nullify(grid_data%yy)
      nullify(grid_data%zz)

      nullify(grid_data%dx)
      nullify(grid_data%dy)
      nullify(grid_data%dz)
      nullify(grid_data%dxh)
      nullify(grid_data%dyh)
      nullify(grid_data%dzh)

      nullify(grid_data%nxv)
      nullify(grid_data%nyv)
      nullify(grid_data%nzv)
      nullify(grid_data%nxgl)
      nullify(grid_data%nygl)
      nullify(grid_data%nzgl)

      nullify(grid_data%istartx)
      nullify(grid_data%istarty)
      nullify(grid_data%istartz)
      nullify(grid_data%istartp)

      nullify(grid_data%mg_ratio_x)
      nullify(grid_data%mg_ratio_y)
      nullify(grid_data%mg_ratio_z)
      nullify(grid_data%iline)
      nullify(grid_data%jline)
      nullify(grid_data%kline)

      nullify(grid_data%gmetric)
      nullify(grid_data%g_pack)
      
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

      grid_data%gxmin = cmin(1)
      grid_data%gymin = cmin(2)
      grid_data%gzmin = cmin(3)
      grid_data%gxmax = cmax(1)
      grid_data%gymax = cmax(2)
      grid_data%gzmax = cmax(3)

      grid_data%gcw   = 1
      grid_data%aligner   = 1
      
      call setup_f90c_ptrlinks(c_grid_mg_def_ptr, grid_data)

      end subroutine f_create_grid_mg_def


c f_create_patch_data
c######################################################################
      subroutine f_create_patch_data(c_patch_ptr,u_0,u_n,aux,gparams)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, bcgrp, patch
 
      integer*8 :: c_patch_ptr
      type(var_array),INTENT(INOUT) :: u_0, u_n
      type(aux_array),INTENT(INOUT) :: aux
      type(grid_mg_def),INTENT(INOUT) :: gparams
      type(patch), pointer :: patch_data

      INTERFACE
         subroutine setup_patch(fpatch,u_0,u_n,aux,gparams)
            use grid_def_st, only: grid_mg_def
            use variable_setup, only: var_array, aux_array, patch
            implicit none
            type(patch), pointer :: fpatch
            type(var_array), TARGET :: u_0, u_n
            type(aux_array), TARGET :: aux
            type(grid_mg_def), TARGET :: gparams
         end subroutine setup_patch
      END INTERFACE

      allocate(patch_data)
      nullify(patch_data%bc_grp)
      nullify(patch_data%u_0)
      nullify(patch_data%u_n)
      nullify(patch_data%u_np)
      nullify(patch_data%u_nm)
      nullify(patch_data%u_graph)
      nullify(patch_data%gparams)
      nullify(patch_data%aux)

      call setup_patch(patch_data,u_0,u_n,aux,gparams)
      call setup_f90c_ptrlinks(c_patch_ptr, patch_data)

      end subroutine f_create_patch_data


c fill_var_arrayf
c######################################################################
      subroutine fill_var_array(fvararray,i,c_array,xs,ys,zs,xe,ye,ze)

      use variable_setup
c      use grid_mgmt
      implicit none

      type(var_array),INTENT(INOUT) :: fvararray
      integer,INTENT(IN)       :: i,xs,ys,zs,xe,ye,ze
      REAL*8, INTENT(IN)       :: c_array(xs:xe,ys:ye,zs:ze)

      call setup_var_array(fvararray,i,c_array,xs,ys,zs,xe,ye,ze)

      end subroutine fill_var_array


c fill_aux_array_var
c######################################################################
      subroutine fill_aux_array_var(fauxvar,i,c_array,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
c      use grid_mgmt

      type(aux_array),INTENT(INOUT) :: fauxvar
      integer,INTENT(IN)       :: i,xs,ys,zs,xe,ye,ze
      REAL*8, INTENT(IN)       :: c_array(xs:xe,ys:ye,zs:ze)

      call setup_aux_var(fauxvar,i,c_array,xs,ys,zs,xe,ye,ze)

      end subroutine fill_aux_array_var


c ill_aux_array_vec
c######################################################################
      subroutine fill_aux_array_vec(fauxvec,i,c_array,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
c      use grid_mgmt

      type(aux_array),INTENT(INOUT) :: fauxvec
      integer,INTENT(IN)       :: i,xs,ys,zs,xe,ye,ze
      REAL*8, INTENT(IN)       :: c_array(xs:xe,ys:ye,zs:ze,1:3)

      call setup_aux_vec(fauxvec,i,c_array,xs,ys,zs,xe,ye,ze)

      end subroutine fill_aux_array_vec


c f_delete_patch_data
c######################################################################
      subroutine f_delete_patch_data(fpatch)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
      use parameters

      type(patch), target :: fpatch
      type(patch), pointer :: patch_data

      patch_data => fpatch;
      deallocate(patch_data%u_0)
      deallocate(patch_data%u_n) 
      deallocate(patch_data%aux)
      deallocate(patch_data%gparams)
      deallocate(patch_data)

      end subroutine f_delete_patch_data


c f_delete_var_array
c######################################################################
      subroutine f_delete_var_array(farray)

      use variable_setup, only: var_array

      type(var_array), target :: farray
      type(var_array), pointer :: array_data

      array_data => farray;
      deallocate(array_data)

      end subroutine f_delete_var_array


c setup_var_array
c######################################################################
      subroutine setup_var_array(farray,i,c_array_ptr,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: var_array
      implicit none

      type(var_array) :: farray
      integer    :: i,xs,ys,zs,xe,ye,ze
      REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze)

      farray%array_var(i+1)%array => c_array_ptr

      end subroutine setup_var_array


c setup_aux_var
c######################################################################
      subroutine setup_aux_var(fauxvar,i,c_array_ptr,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
      implicit none

      type(aux_array) :: fauxvar
      integer    :: i,xs,ys,zs,xe,ye,ze
      REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze)
      fauxvar%var_list(i+1)%array => c_array_ptr

      end subroutine setup_aux_var


c setup_aux_vec
c######################################################################
      subroutine setup_aux_vec(fauxvec,i,c_array_ptr,xs,ys,zs,xe,ye,ze)

      use variable_setup, only: aux_array
      implicit none

      type(aux_array) :: fauxvec
      integer    :: i,xs,ys,zs,xe,ye,ze
      REAL*8, TARGET :: c_array_ptr(xs:xe,ys:ye,zs:ze,1:3)

      fauxvec%vec_list(i+1)%vec => c_array_ptr

      end subroutine setup_aux_vec


c setup_patch
c######################################################################
      subroutine setup_patch(fpatch,u_0,u_n,aux,gparams)

      use grid_def_st, only: grid_mg_def
      use variable_setup, only: var_array, aux_array, patch
      implicit none

      type(patch), pointer :: fpatch
      type(var_array), TARGET :: u_0, u_n
      type(aux_array), TARGET :: aux
      type(grid_mg_def), TARGET :: gparams

      fpatch%u_0 => u_0
      fpatch%u_n => u_n
      fpatch%aux => aux
      fpatch%gparams => gparams

      end subroutine setup_patch

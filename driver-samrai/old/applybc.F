c formInitialCondition
c######################################################################
      subroutine applybc(i,patch_var)

c----------------------------------------------------------------------
c     Applies the boundary conditions
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use generalpurposefunctions

      use auxiliaryVariables

      use grid_aliases

      implicit none

c Call variables

      integer :: i
      type(patch), TARGET :: patch_var

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c Interfaces
 
      INTERFACE
        subroutine applyDepVarBC(bc_inf,varray,vaux,iigx,iigy,iigz)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine applyDepVarBC
      END INTERFACE
 
      INTERFACE
        subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine applyAuxVarBC
      END INTERFACE

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions
      varray => gv%u_n
      vaux => gv%aux
      u_0 => gv%u_0
      xx  => gv%gparams%xx
      yy  => gv%gparams%yy
      zz  => gv%gparams%zz
      dxh => gv%gparams%dxh
      dyh => gv%gparams%dyh
      dzh => gv%gparams%dzh
      dx  => gv%gparams%dx
      dy  => gv%gparams%dy
      dz  => gv%gparams%dz

c Set pointers for data
      pcnv => gv%aux_spec%pcnv
      vecnv => gv%aux_spec%vecnv
      vdummy => gv%aux_spec%vdummy
      vcnv => gv%aux_spec%vcnv
      vcov => gv%aux_spec%vcov
      grd => gv%aux_spec%grd
#ifdef defined(vec_pot)
      acnv => gv%aux_spec%acnv
      acov => gv%aux_spec%acov
      vlap_a => gv%aux_spec%vlap_a
#else 
      bcnv => gv%aux_spec%bcnv
      bcov => gv%aux_spec%bcov
#endif
      b_n => gv%aux_spec%b_n
      dv_dt => gv%aux_spec%dv_dt
      p_n => gv%aux_spec%p_n
      eeta => gv%aux_spec%eeta
      nuu => gv%aux_spec%nuu
      divrgV => gv%aux_spec%divrgV

c Set boundary conditions

      if (gv%bc_seq(i+1,1)>0) then  !Process dep. variable
        call applyDepVarBC(gv%bc_seq(i+1,:),varray,vaux,1,1,1)
      else                          !Process aux. variable
        call applyAuxVarBC(gv%bc_seq(i+1,:),varray,vaux,1,1,1)
      endif

c End program

      end subroutine applybc



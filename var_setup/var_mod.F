
c module variables
c #####################################################################
      module variables

        use var_bc

#if defined(samrai)
        integer :: bc_order = 1
#endif

      end module variables

#if defined(samrai)
c     applybc
c     #################################################################
      subroutine applybc(patch_var,nbc_seq,bc_seq)

c     -----------------------------------------------------------------
c     Applies the boundary conditions for SAMRAI
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer :: nbc_seq, bc_seq(nbc_seq,3)
      type(patch), TARGET :: patch_var

c     Local variables

      integer :: i
      type(var_array),pointer :: varray => null()

c     Interfaces

      INTERFACE
        subroutine applyDepVarBC(igr,bc_inf,varray,order)
          use var_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          integer,optional :: order
        end subroutine applyDepVarBC
      END INTERFACE

      INTERFACE
        subroutine applyAuxVarBC(igr,bc_inf,varray,vaux,order)
          use var_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyAuxVarBC
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use var_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

c     Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

      call setAppBCs(varray,gv%aux)
cc      call allocAuxVariables(gv%aux)

c     Set boundary conditions

      !Impose BCs on group
      do i = 1,nbc_seq
        if (bc_seq(i,1)*AUX < 0) then !Process dep. variable
          call applyDepVarBC(1,bc_seq(i,:),varray,order=bc_order)
        elseif (bc_seq(i,1)*AUX > 0) then !Process aux. variable
          call applyAuxVarBC(1,bc_seq(i,:),varray,gv%aux
     .                      ,order=bc_order)
        endif
      enddo

c     Nullify aux pointers

cc      call deallocLocalAuxVars
      call destroyAppBCs

c     End program

      end subroutine applybc

c     postproc
c     #################################################################
      subroutine postproc(patch_var,nbc_seq,bc_seq)

c     -----------------------------------------------------------------
c     Applies the postprocessor for SAMRAI
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer :: nbc_seq, bc_seq(nbc_seq,3)
      type(patch), TARGET :: patch_var

c     Local variables

      integer :: i
      type(var_array),pointer :: varray => null()

c     Interfaces

      INTERFACE
        subroutine postproc_AuxVar(igr,bc_inf,varray,vaux)
          use var_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine postproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use var_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

c     Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions
      varray => gv%u_n
      u_0    => gv%u_0

      call setAppBCs(varray,gv%aux)
cc      call allocAuxVariables(gv%aux)

c     Postprocess auxiliary quantities

      do i = 1,nbc_seq
        if (bc_seq(i,1)*AUX > 0) then
          call postproc_AuxVar(1,bc_seq(i,:),varray,gv%aux)
        endif
      enddo

c     Nullify aux pointers

cc      call deallocLocalAuxVars
      call destroyAppBCs

c     End program

      end subroutine postproc

c     initializeAuxVar
c     #################################################################
      subroutine initializeAuxVar(patch_var,igrp)

c     -----------------------------------------------------------------
c     Fills auxiliary variable arrays
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      type(patch), TARGET :: patch_var
      integer:: igrp,i

c     Local variables

      type(var_array),pointer :: varray => null()

      INTERFACE
        subroutine preproc_AuxVar(igr,bc_inf,varray,vaux,order)
          use var_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine preproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use var_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

c     Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

      call setAppBCs(varray,gv%aux)
cc      call allocAuxVariables(gv%aux)

c     Fill auxiliary variable arrays

      do i = 1,gv%bc_grp(igrp)%nbc_seq
        if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
          call preproc_AuxVar(1,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                       ,gv%aux,order=bc_order)
        endif
      enddo

c     Nullify pointers

cc      call deallocAuxVariables(gv%aux)
      call destroyAppBCs

      end subroutine initializeAuxVar

c     setupVarInitSeq
c     #################################################################
      subroutine setupVarInitSeq(patch_var,it,inherit)
c     -----------------------------------------------------------------
c     Sets up application BC for SAMRAI
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer :: it
      type(patch), TARGET :: patch_var
      logical*1 :: inherit

c     Local variables

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use var_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

      gv => patch_var

c     Set application BCs (external)

      call setAppBCs(gv%u_n,gv%aux)

c     Set BC schedule

      call create_bc_schedule(gv%u_n,gv%aux,gv%nbc_grp,gv%bc_grp
     .            ,logical(inherit),.false.,.false.,.false.)

c     End program

      end subroutine setupVarInitSeq

#endif

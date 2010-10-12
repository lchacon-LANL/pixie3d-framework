c getNumberOfBCGroups
c ######################################################################
      subroutine getNumberOfBCGroups(patch_var, nbcgrps)
cc      use variable_setup
      implicit none

      type(patch), TARGET :: patch_var
      integer :: nbcgrps

cc      gv => patch_var
cc      grid_params => gv%gparams
cc      gmetric => gv%gparams%gmetric
      
c     Setup SAMRAI f90 pointers
      
      nbcgrps = patch_var%nbc_grp

      end subroutine getNumberOfBCGroups

cc  L. Chacon, 10/12/2010: I have commented out the routine below and added this
cc                         on the application side in applyBoundaryConditions.F.
ccc setupVariableInitializationSequence
ccc ######################################################################
cc      subroutine setupVariableInitializationSequence(patch_var, it)
cc      use parameters
cc      
cc      use variable_setup
cc      
cc      use grid
cc
cc      implicit none
cc      
ccc     Call variables
cc      
cc      type(patch), TARGET :: patch_var
cc      integer :: it
cc
cccc      gv => patch_var
cccc      grid_params => gv%gparams
cccc      gmetric => gv%gparams%gmetric
cc
cc      call setup_app_BC(patch_var,it)     !"it" indicates time level, with it=0 indicating
cc                                          !equilibrium setup. It needs to come via routine call.
cc
cc      end subroutine setupVariableInitializationSequence

c getBC
c ######################################################################
      subroutine getBC(patch_var, igrp, nbc_seq,bc_seq)

c ----------------------------------------------------------------------
c     Initializes BC pointers
c ----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      implicit none

c Call variables

      type(patch), TARGET :: patch_var
      integer :: igrp, nbc_seq, bc_seq

c Begin program

cc      gv => patch_var
cc      grid_params => gv%gparams
cc      gmetric => gv%gparams%gmetric
cc      u_n => gv%u_n
cc      u_0 => gv%u_0
cc      call setVectorDimensions

c Setup SAMRAI f90 pointers
      
      nbc_seq = patch_var%bc_grp(igrp)%nbc_seq
      call setup_f90c_ptrlinks(bc_seq,patch_var%bc_grp(igrp)%bc_seq)

c End program
 
      end subroutine getBC



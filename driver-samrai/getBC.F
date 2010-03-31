c     getNumberOfBCGroups
      subroutine getNumberOfBCGroups(patch_var, nbcgrps)
      use variable_setup
      implicit none

      type(patch), TARGET :: patch_var
      integer :: nbcgrps

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      
c     Setup SAMRAI f90 pointers
      
      nbcgrps = patch_var%nbc_grp

      end subroutine getNumberOfBCGroups

c     setupVariableInitializationSequence
      subroutine setupVariableInitializationSequence(patch_var, it)
      use parameters
      
      use variable_setup
      
      use grid

      implicit none
      
c     Call variables
      
      type(patch), TARGET :: patch_var
      integer :: it

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric

      call setup_app_BC(it)     !Bobby, "it" indicates time level, with it=0 indicating
                             !equilibrium setup. It needs to come via routine call.
                             !Left it undefined so that it breaks at compilation time so
                                !that you fix it on your end.

      end subroutine setupVariableInitializationSequence

c     getBC
c######################################################################
      subroutine getBC(patch_var, igrp, nbc_seq,bc_seq)

c----------------------------------------------------------------------
c     Initializes BC pointers
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      implicit none

c Call variables

      type(patch), TARGET :: patch_var
      integer :: igrp, nbc_seq, bc_seq

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
cc      u_n => gv%u_n
cc      u_0 => gv%u_0
cc      call setVectorDimensions

c Setup SAMRAI f90 pointers
      
      nbc_seq = patch_var%bc_grp(igrp)%nbc_seq
      call setup_f90c_ptrlinks(bc_seq,patch_var%bc_grp(igrp)%bc_seq)

c End program
 
      end subroutine getBC



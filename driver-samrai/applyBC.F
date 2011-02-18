c getNumberOfBCGroups
c ######################################################################
      subroutine getNumberOfBCGroups(patch_var, nbcgrps)
      use variable_setup
      implicit none

      type(patch), TARGET :: patch_var
      integer :: nbcgrps
      
      nbcgrps = patch_var%nbc_grp

      end subroutine getNumberOfBCGroups


c getBC
c ######################################################################
      subroutine getBCSequence(patch_var, igrp, nbc_seq, bc_seq)
      use parameters
      use variable_setup
      use grid
      implicit none

      type(patch), TARGET :: patch_var
      integer :: igrp, nbc_seq, bc_seq

      nbc_seq = patch_var%bc_grp(igrp)%nbc_seq
      call setup_f90c_ptrlinks(bc_seq,patch_var%bc_grp(igrp)%bc_seq)
 
      end subroutine getBCSequence



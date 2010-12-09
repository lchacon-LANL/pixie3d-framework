c getNumberOfBCGroups
c ######################################################################
      subroutine getNumberOfBCGroups(patch_var, nbcgrps)
      use variable_setup
      implicit none

      type(patch), TARGET :: patch_var
      integer :: nbcgrps

c     Setup SAMRAI f90 pointers
      
      nbcgrps = patch_var%nbc_grp

      end subroutine getNumberOfBCGroups

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

c Setup SAMRAI f90 pointers
      
      nbc_seq = patch_var%bc_grp(igrp)%nbc_seq
      call setup_f90c_ptrlinks(bc_seq,patch_var%bc_grp(igrp)%bc_seq)

c End program
 
      end subroutine getBC



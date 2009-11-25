c getBC
c######################################################################
      subroutine getBC(patch_var,nbc_seq,bc_seq)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      implicit none

c Call variables

c      integer         :: imin,imax,jmin,jmax,kmin,kmax,ieq

      integer :: nbc_seq, bc_seq
      type(patch), TARGET :: patch_var
      nbc_seq = patch_var%nbc_seq
      call setup_f90c_ptrlinks(bc_seq,patch_var%bc_seq)


c End program
 
      end subroutine getBC



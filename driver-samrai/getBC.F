c getBC
c######################################################################
      subroutine getBC(patch_var,nbc_seq,bc_seq, it)

c----------------------------------------------------------------------
c     Initializes BC pointers
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      implicit none

c Call variables

      integer :: nbc_seq, bc_seq, it
      type(patch), TARGET :: patch_var

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
cc      u_n => gv%u_n
cc      u_0 => gv%u_0
cc      call setVectorDimensions

c Setup application BCs (external)

      call setup_app_BC(it)  !Bobby, "it" indicates time level, with it=0 indicating
                             !equilibrium setup. It needs to come via routine call.
                             !Left it undefined so that it breaks at compilation time so
                             !that you fix it on your end.

c Setup SAMRAI f90 pointers

      nbc_seq = patch_var%nbc_seq
      call setup_f90c_ptrlinks(bc_seq,patch_var%bc_seq)

c End program
 
      end subroutine getBC



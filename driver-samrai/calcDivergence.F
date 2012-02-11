c calcDivergence
c####################################################################
      subroutine calcDivergence(patch_var,div)
c--------------------------------------------------------------------
c     Calculates the divergence of a quantity
c--------------------------------------------------------------------

      use parameters
      use grid
      use variables
      implicit none

c Call variables

      type(patch), TARGET :: patch_var
      integer :: nx, ny, nz
      real(8),dimension(:,:,:) :: div

c Interfaces


c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      nx = patch_var%gparams%nlx
      ny = patch_var%gparams%nly
      nz = patch_var%gparams%nlz      


c Deallocate structures

c End program

      end subroutine calcDivergence


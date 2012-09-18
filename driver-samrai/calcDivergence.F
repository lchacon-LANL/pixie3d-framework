c calcDivergence
c####################################################################
      subroutine calcDivergence(patch_var,div,nx,ny,nz)
c--------------------------------------------------------------------
c     Calculates the divergence of a quantity
c--------------------------------------------------------------------

      use parameters
      use grid
      use variables
      use grid_vertex_ops
      implicit none

c Call variables

      type(patch), TARGET :: patch_var
      integer :: nx, ny, nz, i, j, k, igrid
      real(8) :: div(nx,ny,nz)
      real(8),pointer :: cnv(:,:,:,:)

c Interfaces


c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      cnv => patch_var%aux%vec_list(5)%vec
      igrid = 1    

c Calculate the divergence of B
      do k=1,nz
        do j=1,ny
          do i=1,nx
            div(i,j,k) = div_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,cnv)
          enddo
        enddo
      enddo   

c End program

      end subroutine calcDivergence


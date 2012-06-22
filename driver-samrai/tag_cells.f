      subroutine tag_cells(
     & lo, hi, gcw, gcw_tag, J0, J_tol, tag_array )

c  Tag cells for refinement

      implicit none

      integer lo(0:2)                               ! Low index
      integer hi(0:2)                               ! High index
      integer gcw(0:2)                              ! Ghost cell width to use for cell-centered variables
      integer gcw_tag(0:2)                          ! Ghost cell width to use for tagging
      double precision J0(lo(0)-gcw(0):hi(0)+gcw(0),
     &   lo(1)-gcw(1):hi(1)+gcw(1),
     &   lo(2)-gcw(2):hi(2)+gcw(2), 3 )             ! Cell-centered current density
      double precision J_tol                        ! Cell-centered electron energy density
      integer tag_array(lo(0)-gcw_tag(0):hi(0)+gcw_tag(0),
     &   lo(1)-gcw_tag(1):hi(1)+gcw_tag(1),
     &   lo(2)-gcw_tag(2):hi(2)+gcw_tag(2))         ! Cell-centered array indicating which cells are tagged (output)
      
      integer i, j, k
      double precision J_mag, J_tol2

c  Initialize the tag_array
      tag_array(:,:,:) = 0
      
c  Set the tag_array
      J_tol2 = J_tol*J_tol
      do i = lo(0), hi(0)
         do j = lo(1), hi(1)
            do k = lo(2), hi(2)
               J_mag = J0(i,j,k,1)*J0(i,j,k,1) + 
     &                 J0(i,j,k,2)*J0(i,j,k,2) + 
     &                 J0(i,j,k,3)*J0(i,j,k,3) 
               if ( J_mag >= J_tol2 ) then
                  tag_array(i,j,k) = 1
               end if
            end do
         end do
      end do

      return
      end


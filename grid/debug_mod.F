c module debug
c ######################################################################
      module debug_mod

        real(8) :: lxmax,lxmin,lymax,lymin
        real(8),allocatable,dimension(:,:,:,:) :: dbg

      contains

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c     ---------------------------------------------------------------------
c     Contours 2D array in xdraw format. In call:
c       * arr: 2D array to be plotted
c       * nx,ny: dimensions of array
c       * xmin,xmax,ymin,ymax: 2D domain limits
c       * iopt: whether to initialize xdraw plot (iopt=0) or not.
c       * nunit: integer file identifier.
c     ---------------------------------------------------------------------

c     Call variables

      integer(4) :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer(4) :: i,j

c     Begin program

      if(iopt == 0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4)
     .              ,real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*)

c     End program

      end subroutine contour

      end module debug_mod

c module mg_setup
c######################################################################
      module mg_setup

        integer :: ngrd
        double precision,dimension(:,:),allocatable:: xl,yl
        double precision,dimension(:)  ,allocatable:: dx,dy
        integer         ,dimension(:)  ,allocatable:: istartp,
     .                                                nxvp,nyvp,ntotvp

        double precision :: xlength,ylength  !2D grid dimension

      contains

c     setupMG
c     #################################################################
      subroutine setupMG(nx,ny,bcond)

c     -----------------------------------------------------------------
c     Initializes MG and creates 2D uniform grid.
c
c     Input variables:
c         * nx : number of mesh points in X
c         * ny : number of mesh points in Y
c         * bcond : vector with definition of boundary conditions
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer          :: nx,ny,bcond(4)
cc      double precision :: xlength,ylength
cc      logical          :: xper,yper

c     Local variables

      integer          :: n1,n2,i,j,ig,if,jf,nx1,ny1

c     Begin program

c     Find number of grid levels (ngrd)

      if (ngrd.eq.0) then
        n1 = int(dlog(1d0*nx)/dlog(2d0)+0.001)
        n2 = int(dlog(1d0*ny)/dlog(2d0)+0.001)
        ngrd = min(n1,n2)-1
cc        ngrd = min(n1,n2)-2
        do i = ngrd,1,-1
          nx1 = nx/2**(i-1)
          ny1 = ny/2**(i-1)
          if (nx1*2**(i-1).eq.nx.and.ny1*2**(i-1).eq.ny) exit
        enddo
        ngrd = i
      else
        nx1 = nx/2**(ngrd-1)
        ny1 = ny/2**(ngrd-1)
      endif

c     Allocate MG arrays

      allocate(nxvp(ngrd),nyvp(ngrd))
      allocate(dx(ngrd),dy(ngrd))
      allocate(xl(0:nx+1,ngrd),yl(0:ny+1,ngrd))
      allocate(istartp(ngrd),ntotvp(ngrd))

c     Initialize MG arrays

      nxvp(1) = nx1
      nyvp(1) = ny1
      istartp(1) = 1
      ntotvp(1)  = nxvp(1)*nyvp(1)
      do i = 2,ngrd
        nxvp(i) = nxvp(i-1) * 2
        nyvp(i) = nyvp(i-1) * 2
        istartp(i) = nxvp(i-1)*nyvp(i-1) + istartp(i-1)
        ntotvp(i)  = nxvp(i)*nyvp(i)
      enddo

c     Error check

      if (nxvp(ngrd).ne.nx.or.nyvp(ngrd).ne.ny) then
         write (*,*) 'MG grids do not match original grid. Aborting.'
         stop
      endif

c     Define uniform grid at all grid levels

      call defineGrid(nx,xl,dx,nxvp,ngrd,xlength,bcond(3),bcond(4))

      call defineGrid(ny,yl,dy,nyvp,ngrd,ylength,bcond(1),bcond(2))

c     End program

      contains

      subroutine defineGrid (nn,xx,dx,nx,ngrid,length,bc1,bc2)

        implicit none

        integer          :: nn,ngrid,nx(ngrid),bc1,bc2
        double precision :: xx(0:nn+1,ngrid),dx(ngrid),length
        
        integer          :: ig,i

        if (bc1 == 0 .or. bc2 == 0) then

          ig = ngrid

          dx(ig) = length/dfloat(nx(ig)-1)

          xx(1,ig) = 0d0
          do i = 2,nx(ig)+1
            xx(i,ig) = xx(i-1,ig) + dx(ig)
          enddo
          xx(0,ig) = xx(1,ig) - dx(ig)

          do ig = ngrid-1,1,-1

cc            dx(ig) = dx(ig+1)*2d0
cc
cc            do i = 1,nx(ig)
cc              if = 2*i
cc              xx(i,ig) = .5*(xx(if,ig+1)+xx(if-1,ig+1))
cc            enddo
cc            xx(0       ,ig) = xx(1     ,ig) - dx(ig)/2.
cc            xx(nx(ig)+1,ig) = xx(nx(ig),ig) + dx(ig)/2.

            dx(ig) = length/dfloat(nx(ig)-1)

            xx(1,ig) = 0d0
            do i = 2,nx(ig)+1
              xx(i,ig) = xx(i-1,ig) + dx(ig)
            enddo
            xx(0,ig) = xx(1,ig) - dx(ig)

          enddo

        elseif (bc1 == 1 .and. bc2 == 1) then

          ig = ngrid

          dx(ig) = length/dfloat(nx(ig)+1)

          xx(0,ig) = 0d0
          do j = 1,nx(ig)+1
            xx(j,ig) = xx(j-1,ig) + dx(ig)
          enddo

          do ig = ngrid-1,1,-1

            dx(ig) = length/dfloat(nx(ig)+1)

            xx(0,ig) = 0d0
            do j = 1,nx(ig)+1
              xx(j,ig) = xx(j-1,ig) + dx(ig)
            enddo

          enddo

        elseif (bc1 == 2 .and. bc2 == 1) then

          ig = ngrid

          dx(ig) = length/dfloat(nx(ig))

          xx(0,ig) = -dx(ig)
          do j = 1,nx(ig)+1
            xx(j,ig) = xx(j-1,ig) + dx(ig)
          enddo

          do ig = ngrid-1,1,-1

            dx(ig) = length/dfloat(nx(ig))

            xx(0,ig) = -dx(ig)
            do j = 1,nx(ig)+1
              xx(j,ig) = xx(j-1,ig) + dx(ig)
            enddo

          enddo

        elseif (bc1 == 1 .and. bc2 == 2) then

          ig = ngrid

          dx(ig) = length/dfloat(nx(ig))

          xx(0,ig) = 0d0
          do j = 1,nx(ig)+1
            xx(j,ig) = xx(j-1,ig) + dx(ig)
          enddo

          do ig = ngrid-1,1,-1

            dx(ig) = length/dfloat(nx(ig))

            xx(0,ig) = 0d0
            do j = 1,nx(ig)+1
              xx(j,ig) = xx(j-1,ig) + dx(ig)
            enddo

          enddo

        elseif (bc1 == 2 .and. bc2 == 2) then

          ig = ngrid

          dx(ig) = length/dfloat(nx(ig)-1)

          xx(0,ig) = -dx(ig)
          do j = 1,nx(ig)+1
            xx(j,ig) = xx(j-1,ig) + dx(ig)
          enddo

          do ig = ngrid-1,1,-1

            dx(ig) = length/dfloat(nx(ig)-1)

            xx(0,ig) = -dx(ig)
            do j = 1,nx(ig)+1
              xx(j,ig) = xx(j-1,ig) + dx(ig)
            enddo

          enddo

        endif

      end subroutine defineGrid

      end subroutine setupMG

c     restrictArray
c     ###############################################################
      subroutine restrictArray(nx,ny,array,mgvector,order,bcond)
c     --------------------------------------------------------------
c     Restricts array to mgvector in all grids (without ghost nodes)
c     --------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer         nx,ny,order,bcond(4)
      real*8          array(0:nx+1,0:ny+1)
      real*8          mgvector(2*nx*ny)

c     Local variables

      integer*4       nxf,nyf,nxc,nyc,igrid,igridf

      real*8          dx1,dy1

      double precision, allocatable,dimension(:,:)::arrayf,arrayc

c     Begin program

c     Map array in finest grid onto MG vector

      dx1 = dx(ngrd)
      dy1 = dy(ngrd)

      nxf  = nxvp(ngrd)
      nyf  = nyvp(ngrd)

      call mapArrayToMGVector(nxf,nyf,array,mgvector,ngrd)

      allocate(arrayf(0:nxf+1,0:nyf+1))

      arrayf = array

      igridf = ngrd
 
c     Restrict array to coarser grids

      do igrid = ngrd-1,2,-1

c       Characterize coarse grid and define arrays

        dx1 = dx(igrid)
        dy1 = dy(igrid)

        nxc  = nxvp(igrid)
        nyc  = nyvp(igrid)

        allocate(arrayc(0:nxc+1,0:nyc+1))

c       Restrict arrayf -> arrayc

        call restrictArraytoArray(nxf,nyf,arrayf,igridf
     .                           ,nxc,nyc,arrayc,igrid ,order,bcond)

c       Map arrayc onto MG vector

        call mapArrayToMGVector(nxc,nyc,arrayc,mgvector,igrid)

c       Transfer grid information

        if (order.eq.0) then
          igridf = igrid
          nxf = nxc
          nyf = nyc

          deallocate(arrayf)
          allocate(arrayf(0:nxf+1,0:nyf+1))

          arrayf = arrayc
        endif

c       Deallocate variables

        deallocate(arrayc)

      enddo

c     End program

      deallocate(arrayf)

      return
      end subroutine

c     mapArrayToMGVector
c     ###############################################################
      subroutine mapArrayToMGVector(nx,ny,array,mgvector,igrid)
c     --------------------------------------------------------------
c     Maps array into a MG vector excluding ghost nodes.
c     --------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer ::      igrid,nx,ny
      real*8          mgvector(2*nx*ny),array(0:nx+1,0:ny+1)

c     Local variables

      integer*4       i,j,ii

c     Begin program

      do i = 1,nx
        do j = 1,ny
          ii = i + nx*(j-1) + istartp(igrid)-1
          mgvector(ii) = array(i,j)
        enddo
      enddo

c     End program

      end subroutine

c     mapMGVectorToArray
c     ###############################################################
      subroutine mapMGVectorToArray(nx,ny,mgvector,array,igrid,bcond)
c     --------------------------------------------------------------
c     Maps a MG vector into an array, filling ghost cells.
c     --------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer ::      igrid,nx,ny,bcond(4)
      real*8          mgvector(2*nx*ny),array(0:nx+1,0:ny+1)

c     Local variables

      integer*4       i,j,ii

c     Begin program

      do i = 1,nx
        do j = 1,ny
          ii = i + nx*(j-1) + istartp(igrid)-1
          array(i,j) = mgvector(ii)
        enddo
      enddo

      call setBoundaryConditions(array,1,nx,1,ny,nx,ny,bcond)

c     End program

      end subroutine

c     restrictArraytoArray
c     ###############################################################
      subroutine restrictArraytoArray(nxf,nyf,arrayf,igridf
     .                               ,nxc,nyc,arrayc,igridc,order,bcond)
c     --------------------------------------------------------------
c     Restricts full array (including ghost nodes) from igridf to
c     igridc with interpolation order 'order'.
c     --------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer ::      igridc,igridf,nxc,nyc,nxf,nyf,order,bcond(4)
      real*8          arrayf(0:nxf+1,0:nyf+1),arrayc(0:nxc+1,0:nyc+1)

c     Local variables

      integer*4       ic,jc,if,jf

      real*8          xxc,yyc,xx(nxf+2),yy(nyf+2)

c     Interpolation

      integer*4    kx,ky,nx,ny,dim,flg
      real*8, dimension(:),allocatable:: tx,ty,work
      real*8, dimension(:,:),allocatable:: bcoef

      real*8       db2val
      external     db2val

c     Begin program

c     Agglomeration

      if (order.eq.0) then

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic

            arrayc(ic,jc) = (arrayf(if,jf  ) + arrayf(if-1,jf  )
     .                      +arrayf(if,jf-1) + arrayf(if-1,jf-1))/4d0
          enddo
        enddo

        call setBoundaryConditions(arrayc,1,nxc,1,nyc,nxc,nyc,bcond)

      else

c     Interpolate fine grid array using splines

        xx(1:nxf+2) = xl(0:nxf+1,igridf)
        yy(1:nyf+2) = yl(0:nyf+1,igridf)

        flg = 0
        kx = order+1
        ky = order+1
        nx = nxf+2
        ny = nyf+2
        dim = nx*ny + max(2*kx*(nx+1),2*ky*(ny+1))

        allocate(tx(nx+kx))
        allocate(ty(nx+ky))
        allocate(work(dim))
        allocate(bcoef(nx,ny))

        call db2ink(xx,nx,yy,ny,arrayf,nx,kx,ky,tx,ty,bcoef,work,flg)

c     Map array into vector vv for all grids

        do jc = 1,nyc
          do ic = 1,nxc

            xxc = xl(ic,igridc)
            yyc = yl(jc,igridc)

            arrayc(ic,jc) = 
     .           db2val(xxc,yyc,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

          enddo
        enddo

        call setBoundaryConditions(arrayc,1,nxc,1,nyc,nxc,nyc,bcond)

        deallocate(tx)
        deallocate(ty)
        deallocate(work)
        deallocate(bcoef)

      endif

c     End program

      end subroutine

      end module mg_setup

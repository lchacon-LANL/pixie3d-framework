c module mg_setup
c######################################################################
      module mg_setup

        integer :: ngrd
        double precision,dimension(:,:),allocatable:: xl,yl
        double precision,dimension(:)  ,allocatable:: dx,dy
        integer         ,dimension(:)  ,allocatable:: istartp,
     .           nxvp,nyvp,ntotvp
        integer         ,dimension(:)  ,allocatable:: istart,
     .           ntotv

      contains

c     setupMG
c     #################################################################
      subroutine setupMG(nx,ny,xlength,ylength)

c     -----------------------------------------------------------------
c     Initializes MG and creates 2D uniform grid
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer         :: nx,ny
      double precision:: xlength,ylength

c     Local variables

      integer*4      n1,n2,i,j,ig,if,jf,nx1,ny1

c     Begin program

c     Find ngrd

      if (ngrd.eq.0) then
        n1 = int(dlog(1d0*nx)/dlog(2d0)+0.001)
        n2 = int(dlog(1d0*ny)/dlog(2d0)+0.001)
        ngrd = min(n1,n2)-1
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

c     Allocate arrays

      allocate(nxvp(ngrd),nyvp(ngrd))
      allocate(dx(ngrd),dy(ngrd))
      allocate(xl(0:nx+1,ngrd),yl(0:ny+1,ngrd))
      allocate(istartp(ngrd),ntotvp(ngrd))
      allocate(istart (ngrd),ntotv (ngrd))

c     Initialize MG

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

c     Define the uniform grid at all levels

      dx(ngrd) = xlength/dfloat(nx-1)
      dy(ngrd) = ylength/dfloat(ny+1)

      ig = ngrd

      xl(1,ig) = 0d0
      do i = 2,nxvp(ig)+1
        xl(i,ig) = xl(i-1,ig) + dx(ig)
      enddo
      xl(0,ig) = xl(1,ig) - dx(ig)

      yl(0,ig) = 0d0
      do j = 1,nyvp(ig)+1
        yl(j,ig) = yl(j-1,ig) + dy(ig)
      enddo

      do ig = ngrd-1,1,-1

cc        dx(ig) = xlength/dfloat(nxvp(ig)-1)
        dx(ig) = dx(ig+1)*2d0

        do i = 1,nxvp(ig)
          if = 2*i
          xl(i,ig) = .5*(xl(if,ig+1)+xl(if-1,ig+1))
        enddo
        xl(0         ,ig) = xl(1       ,ig) - dx(ig)/2.
        xl(nxvp(ig)+1,ig) = xl(nxvp(ig),ig) + dx(ig)/2.

cc        dy(ig) = dy(ig+1)*2d0
cc        yl(0,ig) = 0d0
cc        do j = 1,nyvp(ig)
cc          jf = 2*j
cc          yl(j,ig) = .5*(yl(jf,ig+1)+yl(jf-1,ig+1))
cc        enddo
cc        yl(nyvp(ig)+1,ig) = 1d0

        dy(ig) = 1d0/dfloat(nyvp(ig)+1)
        yl(0,ig) = 0d0
        do j = 1,nyvp(ig)+1
          yl(j,ig) = yl(j-1,ig) + dy(ig)
        enddo

cc        write (*,*) 'next grid: ',ig
cc        write (*,*) xl(0:nxvp(ig)+1,ig)
cc        write (*,*)
cc        write (*,*) yl(0:nyvp(ig)+1,ig)
      enddo

c     End program

      return
      end subroutine

      end module mg_setup

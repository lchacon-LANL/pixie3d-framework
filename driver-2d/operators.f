c curl
c####################################################################
      subroutine curl(i,j,nx,ny,dxl,dyl,sf,ax,ay)

c--------------------------------------------------------------------
c     This subroutine finds (ax,ay)=-curl(sf*ez), where sf is a stream 
c     function.
c--------------------------------------------------------------------

      implicit none

c Call variables

      integer*4    nx,ny,bcond(4)
      real*8       sf(0:nx+1,0:ny+1)
      real*8       ax,ay
      real*8       dxl,dyl

c Local variables

      integer*4    i,j,ii,ip,im,jj,jp,jm

c Begin program
 
      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      if (j.eq.0) then
        ax =  (3*sf(ii,jj)-4*sf(ii,jp)+sf(ii,jp+1))/dyl/2d0
        ay =  (sf(ip,jj)-sf(im,jj))/dxl/2d0
      elseif (j.eq.ny+1) then
        ax = -(3*sf(ii,jj)-4*sf(ii,jm)+sf(ii,jm-1))/dyl/2d0
        ay =  (sf(ip,jj)-sf(im,jj))/dxl/2d0
      else
        ax = -(sf(ii,jp)-sf(ii,jm))/dyl/2d0
        ay =  (sf(ip,jj)-sf(im,jj))/dxl/2d0
      endif

c End

      return
      end

c grad
c####################################################################
      subroutine grad(i,j,nx,ny,dxl,dyl,sf,ax,ay)

c--------------------------------------------------------------------
c     This subroutine finds (ax,ay)=grad(sf), where sf is a stream 
c     function.
c--------------------------------------------------------------------

      implicit none

c Call variables

      integer*4    nx,ny,bcond(4)
      real*8       sf(0:nx+1,0:ny+1)
      real*8       ax,ay
      real*8       dxl,dyl

c Local variables

      integer*4    i,j,ii,ip,im,jj,jp,jm

c Begin program
 
      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      if (j.eq.0) then
        ay = -(3*sf(ii,jj)-4*sf(ii,jp)+sf(ii,jp+1))/dyl/2d0
        ax =  (sf(ip,jj)-sf(im,jj))/dxl/2d0
      elseif (j.eq.ny+1) then
        ay =  (3*sf(ii,jj)-4*sf(ii,jm)+sf(ii,jm-1))/dyl/2d0
        ax =  (sf(ip,jj)-sf(im,jj))/dxl/2d0
      else
        ay =  (sf(ii,jp)-sf(ii,jm))/dyl/2d0
        ax =  (sf(ip,jj)-sf(im,jj))/dxl/2d0
      endif

c End

      return
      end

c div
c####################################################################
      subroutine div(i,j,nx,ny,dx,dy,sf,array,dv)

c--------------------------------------------------------------------
c     This subroutine finds div(-curl(sf)*array), where sf is a stream 
c     function, on a co-located mesh using ZIP differencing.
c--------------------------------------------------------------------

      implicit none

c Call variables

      integer*4    nx,ny,bcond(4)
      real*8       dx,dy,dv
      real*8       sf(0:nx+1,0:ny+1),array(0:nx+1,0:ny+1)

c Local variables

      integer*4    i,j,ii,ip,im,jj,jp,jm

      real(8) ::   coeff,dd,ax,ay,ax_e,ax_w,ay_n,ay_s
      real(8) ::   arr,arr_e,arr_w,arr_n,arr_s
      real(8) ::   cflxe,cflxw,cflxn,cflxs

c Begin program
 
      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      coeff = 2d0
      if (j.eq.0) then
cc        jm = jj
        coeff = 1d0
      elseif (j.eq.ny+1) then
cc        jp = jj
        coeff = 1d0
      endif

      arr   = array(ii,jj)
      arr_e = array(ip,jj)
      arr_w = array(im,jj)
      arr_n = array(ii,jp)
      arr_s = array(ii,jm)

cc      call curl(ii,jj,nx,ny,dx,dy,sf,ax,ay)
cc      call curl(ii,jp,nx,ny,dx,dy,sf,dd,ay_n)
cc      call curl(ii,jm,nx,ny,dx,dy,sf,dd,ay_s)
cc      call curl(ip,jj,nx,ny,dx,dy,sf,ax_e,dd)
cc      call curl(im,jj,nx,ny,dx,dy,sf,ax_w,dd)

      ax    = -(sf(ii,jp)-sf(ii,jm))/dy/coeff
      ax_e  = -(sf(ip,jp)-sf(ip,jm))/dy/coeff
      ax_w  = -(sf(im,jp)-sf(im,jm))/dy/coeff

      ay_n  =  (sf(ip,jp)-sf(im,jp))/dx/2d0
      ay_s  =  (sf(ip,jm)-sf(im,jm))/dx/2d0
      ay    =  (sf(ip,jj)-sf(im,jj))/dx/2d0

      cflxe = (ax*arr_e + ax_e*arr)/2d0
      cflxw = (ax*arr_w + ax_w*arr)/2d0
      cflxn = (ay*arr_n + ay_n*arr)/2d0
      cflxs = (ay*arr_s + ay_s*arr)/2d0

      dv = ( (cflxe - cflxw)/dx + (cflxn - cflxs)/dy*2d0/coeff )

c End

      return
      end

c laplace
c####################################################################
      real*8 function laplace(i,j,nx,ny,dxx,dyy,arr)
      implicit none           !For safe fortran
c--------------------------------------------------------------------
c     Calculates lap(arr).
c--------------------------------------------------------------------

c Call variables

      integer*4  i,j,nx,ny,bcond(4)
      real*8     dxx,dyy,arr(0:nx+1,0:ny+1)

c Local variables

      integer*4  ip,im,jp,jm,ii,jj

c Begin program

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      if (j.eq.0) then
        laplace = ( (arr(ip,jj)-arr(ii,jj))/dxx
     .             -(arr(ii,jj)-arr(im,jj))/dxx )/dxx
     .          + ( 2*arr(ii,jj  )-5*arr(ii,jp  )
     .             +4*arr(ii,jp+1)-  arr(ii,jp+2) )/dyy**2
      elseif (j.eq.ny+1) then
        laplace = ( (arr(ip,jj)-arr(ii,jj))/dxx
     .             -(arr(ii,jj)-arr(im,jj))/dxx )/dxx
     .          + ( 2*arr(ii,jj  )-5*arr(ii,jm  )
     .             +4*arr(ii,jm-1)-  arr(ii,jm-2) )/dyy**2
      else
        laplace = ( (arr(ip,jj)-arr(ii,jj))/dxx
     .             -(arr(ii,jj)-arr(im,jj))/dxx )/dxx
     .          + ( (arr(ii,jp)-arr(ii,jj))/dyy
     .             -(arr(ii,jj)-arr(ii,jm))/dyy )/dyy
      endif

c End program

      return
      end

c laplace_new
c####################################################################
      function laplace_new(i,j,nx,ny,igr,arr) result(laplace)
c--------------------------------------------------------------------
c     Calculates lap(arr).
c--------------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c Call variables

      integer*4  i,j,nx,ny,bcond(4),igr
      real*8     arr(0:nx+1,0:ny+1),laplace

c Local variables

      integer*4  ip,im,jp,jm,ii,jj,ig,jg
      real(8)  :: dxx,dyy

c Begin program

      call getMGmap(i,j,igr,ig,jg)

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      
      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      if (j.eq.0) then
        laplace = ( (arr(ip,jj)-arr(ii,jj))/dxx
     .             -(arr(ii,jj)-arr(im,jj))/dxx )/dxx
     .          + ( 2*arr(ii,jj  )-5*arr(ii,jp  )
     .             +4*arr(ii,jp+1)-  arr(ii,jp+2) )/dyy**2
      elseif (j.eq.ny+1) then
        laplace = ( (arr(ip,jj)-arr(ii,jj))/dxx
     .             -(arr(ii,jj)-arr(im,jj))/dxx )/dxx
     .          + ( 2*arr(ii,jj  )-5*arr(ii,jm  )
     .             +4*arr(ii,jm-1)-  arr(ii,jm-2) )/dyy**2
      else
        laplace = ( (arr(ip,jj)-arr(ii,jj))/grid_params%dx(ig  )
     .             -(arr(ii,jj)-arr(im,jj))/grid_params%dx(ig-1) )/dxx
     .          + ( (arr(ii,jp)-arr(ii,jj))/grid_params%dy(jg  )      
     .             -(arr(ii,jj)-arr(ii,jm))/grid_params%dy(jg-1) )/dyy
      endif

c End program

      return
      end

c divDgrad
c####################################################################
      real*8 function divDgrad(i,j,nx,ny,dxx,dyy,arr,diff)
      implicit none           !For safe fortran
c--------------------------------------------------------------------
c     Calculates div(Diff*grad(arr)).
c--------------------------------------------------------------------

c Call variables

      integer*4  i,j,nx,ny,bcond(4)
      real*8     dxx,dyy,arr(0:nx+1,0:ny+1),diff(0:nx+1,0:ny+1)

c Local variables

      integer*4  ip,im,jp,jm,ii,jj
      real*8     diffip,diffim,diffjp,diffjm

c Begin program

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      if (j.eq.0) then
        write (*,*) 'Error in divDgrad; j=0'
      elseif (j.eq.ny+1) then
        write (*,*) 'Error in divDgrad; j=ny+1'
      else
        diffip = 2./(1./diff(ip,jj)+1./diff(ii,jj))
        diffim = 2./(1./diff(im,jj)+1./diff(ii,jj))
        diffjp = 2./(1./diff(ii,jp)+1./diff(ii,jj))
        diffjm = 2./(1./diff(ii,jm)+1./diff(ii,jj))
        divDgrad =( diffip*(arr(ip,jj)-arr(ii,jj))/dxx
     .             -diffim*(arr(ii,jj)-arr(im,jj))/dxx )/dxx
     .          + ( diffjp*(arr(ii,jp)-arr(ii,jj))/dyy
     .             -diffjm*(arr(ii,jj)-arr(ii,jm))/dyy )/dyy
      endif

c End program

      return
      end

c lap_vec
c####################################################################
      real*8 function lap_vec(i,j,nx,ny,dx1,dy1,zz,array,igrid)

c--------------------------------------------------------------------
c     Calculates lap(zz), where zz is a vector, taking boundary
c     conditions from array.
c--------------------------------------------------------------------

      use parameters

      use mg_setup

      implicit none

c Call variables

      integer*4  i,j,nx,ny,igrid,bcond(4)
      real*8     dx1,dy1,zz(nx*ny),array(0:nxd+1,0:nyd+1)

c Local variables

      integer*4  ii,iip,iim,jj,jjp,jjm,nxx
      integer*4  if,imf,ipf,jf,jmf,jpf

c Begin program

      if (j.lt.1.or.j.gt.ny) then
        lap_vec = 0d0
        return
      endif

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0
      call shiftIndices(i,j,ii,iim,iip,jj,jjm,jjp,nx,ny,1,bcond)

      if (j.eq.1) then
        if = i*2**(ngrd-igrid)
        lap_vec = (zz(iip) -2.*zz(ii) + zz(iim)     )/dx1**2
     .           +(zz(jjp) -2.*zz(ii) + array(if,0) )/dy1**2
      elseif (j.eq.ny) then
        if = i*2**(ngrd-igrid)
        lap_vec = (zz(iip)         -2.*zz(ii) + zz(iim))/dx1**2
     .           +(array(if,nyd+1) -2.*zz(ii) + zz(jjm))/dy1**2
      else
        lap_vec = (zz(iip) -2.*zz(ii) + zz(iim))/dx1**2
     .           +(zz(jjp) -2.*zz(ii) + zz(jjm))/dy1**2
      endif

c End program

      return
      end

c curr_vec
c####################################################################
      real*8 function curr_vec(i,j,nx,ny,dx1,dy1,zz,array,igrid)

c--------------------------------------------------------------------
c     Calculates curr(zz), where zz is a vector, taking BC's from
c     array
c--------------------------------------------------------------------

      use parameters

      use mg_setup

      implicit none      !For safe fortran

c Call variables

      integer*4  i,j,nx,ny,igrid,bcond(4)
      real*8     dx1,dy1,zz(nx*ny),array(0:nx+1,0:ny+1)

c Local variables

      integer*4  ii,iip,iim,jj,jjp,jjm,nxx
      integer*4  if,imf,ipf,jf,jmf,jpf

      real*8     current0,current1,current2

c Externals

      real*8     lap_vec
      external   lap_vec

c Begin program

      if (j.eq.0) then
        current0 = lap_vec(i,1,nx,ny,dx1,dy1,zz,array,igrid)
        current1 = lap_vec(i,2,nx,ny,dx1,dy1,zz,array,igrid)
        current2 = lap_vec(i,3,nx,ny,dx1,dy1,zz,array,igrid)
        curr_vec = 2*current0 - current1
cc        curr_vec = current2 -3*current1 + 3*current0
      elseif (j.eq.ny+1) then
        current0 = lap_vec(i,ny  ,nx,ny,dx1,dy1,zz,array,igrid)
        current1 = lap_vec(i,ny-1,nx,ny,dx1,dy1,zz,array,igrid)
        current2 = lap_vec(i,ny-2,nx,ny,dx1,dy1,zz,array,igrid)
        curr_vec = 2*current0 - current1
cc        curr_vec = current2 -3*current1 + 3*current0
      else
        curr_vec = lap_vec(i,j,nx,ny,dx1,dy1,zz,array,igrid)
      endif

      return
      end

c bgrad
c####################################################################
      real*8 function bgrad(i,j,nx,ny,dx,dy,tb,array,cons)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c     This function computes div[B array].
c--------------------------------------------------------------------

c Call variables

      integer*4  i,j,nx,ny,bcond(4)
      real*8     dx,dy,tb(0:nx+1,0:ny+1)
      real*8     array(0:nx+1,0:ny+1)

      logical    cons

c Local variables

      integer*4  ip,im,jp,jm,ii,jj
      real*8     arr,arr_e,arr_w,arr_n,arr_s
      real*8     bx,bx_e,bx_w,by,by_n,by_s
      real*8     cflxe,cflxw,cflxn,cflxs
      real*8     coeff,ax,ay,dd

c Begin program

      if (cons) then

        call div (i,j,nx,ny,dx,dy,tb,array,bgrad)

      else

        call curl(i,j,nx,ny,dx,dy,tb   ,bx,by)

        call grad(i,j,nx,ny,dx,dy,array,ax,ay)

        bgrad =  bx*ax + by*ay

      endif

c End

      return
      end

c bgrad2
c####################################################################
      real*8 function bgrad2(i,j,nx,ny,dx,dy,bxx,byy,arr1)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c     This function computes (B.grad)^2(arr1*arr2).
c--------------------------------------------------------------------

c Call variables

      integer*4  i,j,nx,ny,bcond(4)
      real*8     dx,dy,bxx(nx*ny),byy(nx*ny)
      real*8     arr1(0:nx+1,0:ny+1)

c Local variables

      integer*4  ij,ipj,imj,ijp,ijm,imjp,ipjm
      integer :: ii,ip,im,jj,jp,jm,dum
      real*8     arr,arr_e,arr_n
      real*8     bx,bx_e,bx_w,bx_n,bx_nw
      real*8     by,by_e,by_n,by_s,by_se

c Begin program

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0

      call shiftIndices(i,j,ij,imj,ipj,dum,ijm,ijp,nx,ny,1,bcond)

      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      imjp = imj + nx
      ijp  = ij  + nx

      ipjm = ipj - nx
      ijm  = ij  - nx

      bx    = bxx(ij)
      bx_e  = bxx(ipj)
      bx_w  = bxx(imj)
      if (j.eq.ny) then
        bx_nw = 0d0
        bx_n  = 0d0
      else
        bx_nw = bxx(imjp)
        bx_n  = bxx(ijp)
      endif

      by    = byy(ij)
      by_n  = byy(ijp)
      by_e  = byy(ipj)
      if (j.eq.1) then
        by_se = 0d0
        by_s  = 0d0
      else
        by_se = byy(ipjm)
        by_s  = byy(ijm)
      endif

      arr   = .5*(bx   + bx_w )*(arr1(ii,jj) - arr1(im,jj) )/dx
     .       +.5*(by   + by_s )*(arr1(ii,jj) - arr1(ii,jm) )/dy

      arr_e = .5*(bx   + bx_e )*(arr1(ip,jj) - arr1(ii,jj) )/dx
     .       +.5*(by_e + by_se)*(arr1(ip,jj) - arr1(ip,jm) )/dy

      arr_n = .5*(bx_n + bx_nw)*(arr1(ii,jp) - arr1(im,jp) )/dx
     .       +.5*(by   + by_n )*(arr1(ii,jp) - arr1(ii,jj) )/dy

      bgrad2 = (.5*(bx_e + bx)*arr_e
     .         -.5*(bx_w + bx)*arr   )/dx
     .       + (.5*(by_n + by)*arr_n
     .         -.5*(by_s + by)*arr   )/dy

cc      bx    = -(tb(i  ,j+1)-tb(i  ,j-1))/dy/2d0
cc      bx_e  = -(tb(i+1,j+1)-tb(i+1,j-1))/dy/2d0
cc      bx_w  = -(tb(i-1,j+1)-tb(i-1,j-1))/dy/2d0
cc      if (j.eq.ny) then
cc        bx_nw = 0d0
cc        bx_n  = 0d0
cc      else
cc        bx_nw = -(tb(i-1,j+2)-tb(i-1,j))/dy/2d0
cc        bx_n  = -(tb(i  ,j+2)-tb(i  ,j))/dy/2d0
cc      endif

cc      by    =  (tb(i+1,j  )-tb(i-1,j  ))/dx/2d0
cc      by_n  =  (tb(i+1,j+1)-tb(i-1,j+1))/dx/2d0
cc      by_s  =  (tb(i+1,j-1)-tb(i-1,j-1))/dx/2d0
cc      if (i.eq.nx) then
cc        by_e  =  (tb(3  ,j  )-tb(i  ,j  ))/dx/2d0
cc        by_se =  (tb(3  ,j-1)-tb(i  ,j-1))/dx/2d0
cc      else
cc        by_e  =  (tb(i+2,j  )-tb(i  ,j  ))/dx/2d0
cc        by_se =  (tb(i+2,j-1)-tb(i  ,j-1))/dx/2d0
cc      endif

cc      arr   = .5*(bx + bx_w)*(arr1(i,j) - arr1(i-1,j) )/dx
cc     .       +.5*(by + by_s)*(arr1(i,j) - arr1(i,j-1) )/dy

cc      arr_e = .5*(bx   + bx_e )*(arr1(i+1,j) - arr1(i  ,j  ) )/dx
cc     .       +.5*(by_e + by_se)*(arr1(i+1,j) - arr1(i+1,j-1) )/dy

cc      arr_n = .5*(bx_n + bx_nw)*(arr1(i,j+1) - arr1(i-1,j+1) )/dx
cc     .       +.5*(by   + by_n )*(arr1(i,j+1) - arr1(i,j    ) )/dy

cc      bgrad2 = (.5*(bx_e + bx)*arr_e
cc     .         -.5*(bx_w + bx)*arr   )/dx
cc     .       + (.5*(by_n + by)*arr_n
cc     .         -.5*(by_s + by)*arr   )/dy

c End

      return
      end

c bgrad_vec
c####################################################################
      real*8 function bgrad_vec(i,j,nx,ny,dx1,dy1,bxx,byy,zz,array,ig)

c--------------------------------------------------------------------
c     This subroutine calculates B.grad(zz), where zz is a vector,
c     taking boundary conditions from array.
c--------------------------------------------------------------------

      use parameters

      use mg_setup

      implicit none

c Call variables

      integer*4      i,j,nx,ny,ig,bcond(4)
      real*8         dx1,dy1,bxx(nx*ny),byy(nx*ny)
      real*8         zz(nx*ny),array(0:nxd+1,0:nyd+1)

c Local variables

      integer*4      ii,iip,iim,jj,jjp,jjm,nxx,iff

c Externals

c Begin program

      if (j.lt.1.or.j.gt.ny) then
        bgrad_vec = 0d0
        return
      endif

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0

      call shiftIndices(i,j,ii,iim,iip,jj,jjm,jjp,nx,ny,1,bcond)

      if (j.eq.1) then
        iff = i*2**(ngrd-ig)
        bgrad_vec = bxx(ii)*(zz(iip) - zz(iim)      )/2d0/dx1
     .            + byy(ii)*(zz(jjp) - array(iff,0) )/2d0/dy1
      elseif (j.eq.ny) then
        iff = i*2**(ngrd-ig)
        bgrad_vec = bxx(ii)*(zz(iip)          - zz(iim))/2d0/dx1
     .            + byy(ii)*(array(iff,nyd+1) - zz(jjm))/2d0/dy1
      else
        bgrad_vec = bxx(ii)*(zz(iip) - zz(iim))/2d0/dx1
     .            + byy(ii)*(zz(jjp) - zz(jjm))/2d0/dy1
      endif

c End program

      return
      end

c bgradj2
c####################################################################
      real*8 function bgradj2(i,j,nx,ny,dx,dy,tb,tj,cons)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c     This function computes div(B.grad(J)).
c--------------------------------------------------------------------

c Call variables

      integer*4  i,j,nx,ny,bcond(4)
      real*8     dx,dy,tj(0:nx+1,0:ny+1),tb(0:nx+1,0:ny+1)

      logical    cons

c Local variables

      integer*4  ip,im,jp,jm,ii,jj
      real*8     arr,arr_e,arr_w,arr_n,arr_s
      real*8     bx,bx_e,bx_w,by,by_n,by_s
      real*8     cflxe,cflxw,cflxn,cflxs
      real*8     coeff

c Externals

      real*8     laplace
      external   laplace

c Begin program

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0

      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      coeff = 2d0
      if (j.eq.0) then
        jm = jj
        coeff = 1d0
      elseif (j.eq.ny+1) then
        jp = jj
        coeff = 1d0
      endif

      arr   = laplace(i ,j,nx,ny,dx,dy,tj)
      arr_e = laplace(ip,j,nx,ny,dx,dy,tj)
      arr_w = laplace(im,j,nx,ny,dx,dy,tj)
      arr_n = laplace(i,jp,nx,ny,dx,dy,tj)
      arr_s = laplace(i,jm,nx,ny,dx,dy,tj)

      if (cons) then

        bx    = -(tb(ii,jp)-tb(ii,jm))/dy/coeff
        bx_e  = -(tb(ip,jp)-tb(ip,jm))/dy/coeff
        bx_w  = -(tb(im,jp)-tb(im,jm))/dy/coeff

        by_n  =  (tb(ip,jp)-tb(im,jp))/dx/2d0
        by_s  =  (tb(ip,jm)-tb(im,jm))/dx/2d0
        by    =  (tb(ip,jj)-tb(im,jj))/dx/2d0

        cflxe = (bx*arr_e + bx_e*arr)/2d0
        cflxw = (bx*arr_w + bx_w*arr)/2d0
        cflxn = (by*arr_n + by_n*arr)/2d0
        cflxs = (by*arr_s + by_s*arr)/2d0

        bgradj2 = ( (cflxe - cflxw)/dx + (cflxn - cflxs)/dy*2d0/coeff )

      else

        bx    = -(tb(ii,jp)-tb(ii,jm))/dy/coeff
        by    =  (tb(ip,jj)-tb(im,jj))/dx/2d0

        bgradj2 =  bx*(arr_e - arr_w)/dx/2d0
     .            +by*(arr_n - arr_s)/dy/coeff

      endif

      return
      end

c bgradj
c####################################################################
      real*8 function bgradj(i,j,nx,ny,dx,dy,tb,tj,cons)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c     This function computes div(B.grad(J)).
c--------------------------------------------------------------------

c Call variables

      integer*4  i,j,nx,ny,bcond(4)
      real*8     dx,dy,tj(0:nx+1,0:ny+1),tb(0:nx+1,0:ny+1)

      logical    cons

c Local variables

      integer*4  ip,im,jp,jm,ii,jj
      real*8     arr(0:nx+1,0:ny+1)

c Externals

      real*8     laplace,bgrad
      external   laplace,bgrad

c Begin program

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0

      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      arr(ii,jj) = laplace(i ,j,nx,ny,dx,dy,tj)
      arr(ip,jj) = laplace(ip,j,nx,ny,dx,dy,tj)
      arr(im,jj) = laplace(im,j,nx,ny,dx,dy,tj)
      arr(ii,jp) = laplace(i,jp,nx,ny,dx,dy,tj)
      arr(ii,jm) = laplace(i,jm,nx,ny,dx,dy,tj)

      bgradj = bgrad(i,j,nx,ny,dx,dy,tb,arr,cons)

      return
      end

c advec
c######################################################################
      real*8 function advec(i,j,nx,ny,dxx,dyy,vxx,vyy,arr,cons,method)
      implicit none        !For safe fortran
c----------------------------------------------------------------------
c     This function computes the advection term, v.grad(arr) at
c     the (i,j) cell, where v is given by its components vxx, vyy.
c
c     cons: logical; decices between conserv. and non-conserv. scheme.
c
c     method: integer; type of advective scheme: 1 -> upwind,
c             2 -> centered, 3-> QUICK, 4 -> SMART, 5 -> smooth SMART
c             6 -> centered high-order, 7 -> explicit TVD
c----------------------------------------------------------------------

c Call variables

      integer*4    i,j,nx,ny,method,bcond(4)

      real*8       dxx,dyy,diff
      real*8       vxx(0:nx+1,0:ny+1),vyy(0:nx+1,0:ny+1)
      real*8       arr(0:nx+1,0:ny+1)

      logical      cons

c Local variables

      real*8       cflxe, cflxw, cflxn, cflxs
     .            ,vveln, vvels, uvele, uvelw
     .            ,phi_e, phi_w, phi_n, phi_s

      real*8       qe,qw,qn,qs,vvel,uvel

      integer ::   ii,im,ip,jj,jm,jp,ipp,imm,jpp,jmm

      integer*4    advx,advy,adv1,adv2

c Externals

      real*8       inter
      external     inter

c Begin program

      bcond(1) = 1
      bcond(2) = 1
      bcond(3) = 0
      bcond(4) = 0

      call shiftIndices(i,j,ii,im,ip,jj,jm,jp,nx,ny,0,bcond)

      ipp = min(ip+1,nx+1)
      imm = max(im-1,0)

      jpp = min(jp+1,nx+1)
      jmm = max(jm-1,0)

      if (cons) then

c Conservative/ZIP differencing

        cflxe = (vxx(ii,jj)*arr(ip,jj)+vxx(ip,jj)*arr(ii,jj))/2d0
        cflxw = (vxx(ii,jj)*arr(im,jj)+vxx(im,jj)*arr(ii,jj))/2d0
        cflxn = (vyy(ii,jj)*arr(ii,jp)+vyy(ii,jp)*arr(ii,jj))/2d0
        cflxs = (vyy(ii,jj)*arr(ii,jm)+vyy(ii,jm)*arr(ii,jj))/2d0

      else

c Nonconservative differencing, in various flavors

        uvele = vxx(i,j)
        vveln = vyy(i,j)

        uvelw = uvele
        vvels = vveln

        phi_e = 0d0
        phi_w = 0d0
        phi_n = 0d0
        phi_s = 0d0

        if (method.eq.7) then
         call tvd(i,j,nx,ny,dxx,dyy,vxx,vyy,arr,phi_e,phi_w,phi_n,phi_s)
        endif

        advx = method
        advy = method

        qe=inter(arr(ipp,jj),arr(ip,jj),arr(i,jj),arr(im,jj)
     .          ,uvele,advx)

        qw=inter(arr(ip,jj),arr(i,jj),arr(im,jj),arr(imm,jj)
     .          ,uvelw,advx)

        if (j.eq.ny) then
          qn=inter(arr(i,jpp),arr(ii,jp),arr(ii,jj),arr(ii,jm)
     .            ,vveln,min(advy,2))
        else
          qn=inter(arr(i,jpp),arr(ii,jp),arr(ii,jj),arr(ii,jm)
     .            ,vveln,advy)
        endif

        if (j.eq.1) then
          qs=inter(arr(ii,jp),arr(ii,jj),arr(ii,jm),arr(ii,jmm)
     .            ,vvels,min(advy,2))
        else
          qs=inter(arr(ii,jp),arr(ii,jj),arr(ii,jm),arr(ii,jmm)
     .            ,vvels,advy)
        endif

cc        if (i.eq.nx) then
cc          qe=inter(arr(3,j),arr(i+1,j),arr(i,j),arr(i-1,j),uvele,advx)
cc        else
cc          qe=inter(arr(i+2,j),arr(i+1,j),arr(i,j),arr(i-1,j),uvele,advx)
cc        endif

cc        if (i.eq.1) then
cc         qw=inter(arr(i+1,j),arr(i,j),arr(i-1,j),arr(nx-2,j),uvelw,advx)
cc        else
cc          qw=inter(arr(i+1,j),arr(i,j),arr(i-1,j),arr(i-2,j),uvelw,advx)
cc        endif

cc        if (j.eq.ny) then
cc          if (vveln.gt.0d0.or.advy.eq.1) then
cc            qn=inter(0d0,arr(ii,jp),arr(ii,jj),arr(ii,jm),vveln,advy)
cc          else
cc            qn=inter(0d0,arr(ii,jp),arr(ii,jj),arr(ii,jm),vveln,2)
cc          endif
cc        else
cc          qn=inter(arr(i,jp+1),arr(ii,jp),arr(ii,jj),arr(ii,jm)
cc     .            ,vveln,advy)
cc        endif

cc        if (j.eq.1) then
cc          if (vvels.lt.0d0.or.advy.eq.1) then
cc            qs=inter(arr(ii,jp),arr(ii,jj),arr(ii,jm),0d0,vvels,advy)
cc          else
cc            qs=inter(arr(ii,jp),arr(ii,jj),arr(ii,jm),0d0,vvels,2)
cc          endif
cc        else
cc          qs=inter(arr(ii,jp),arr(ii,jj),arr(ii,jm),arr(ii,jm-1)
cc     .            ,vvels,advy)
cc        endif

c     Calculate advection fluxes at faces

        cflxe = uvele*qe + phi_e*(arr(ip,jj)-arr(ii,jj))/dxx
        cflxw = uvelw*qw + phi_w*(arr(ii,jj)-arr(im,jj))/dxx
        cflxn = vveln*qn + phi_n*(arr(ii,jp)-arr(ii,jj))/dyy
        cflxs = vvels*qs + phi_s*(arr(ii,jj)-arr(ii,jm))/dyy

      endif

c Calculate advection term

      advec = ( (cflxe - cflxw)/dxx + (cflxn - cflxs)/dyy )

c End

      return
      end

c inter
c####################################################################
      real*8 function inter(q1,q2,q3,q4,vel,advect)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c    This function computes the advection interpolation at control
c    volume face. Options:
c      * advect = 1 or 7 => First-order upwind
c      * advect = 2 => Centered
c      * advect = 3 => High-order upwind (QUICK)
c      * advect = 4 => Monotone high-order upwind (SMART)
c      * advect = 5 => SMART with smooth transition
c      * advect = 6 => Centered, high-order
c
c    Convention:
c
c       ---x-----x--o--x-----x---->
c          q4    q3 ^  q2    q1
c                   ^
c           Location of face
c--------------------------------------------------------------------

c Call variables

      real*8       q1,q2,q3,q4,vel
      integer*4    advect

c Local variables

      real*8       qt1,qt2,qt3,qt4,slp1,slp2,a,b,c,xp1,xp2

c Externals

      real*8       fmed
      external     fmed

c Begin program

c     Upwind

      if     (advect.eq.1.or.advect.eq.7) then
        if (vel.gt.0d0) then
          inter = q3
        else
          inter = q2
        endif

c     Centered

      elseif (advect.eq.2) then
        inter = (q2+q3)/2d0

c     QUICK

      elseif (advect.eq.3) then
        if (vel.gt.0d0) then
          inter = 3.*q2/8. + 3.*q3/4.- q4/8.
        else
          inter = 3.*q3/8. + 3.*q2/4.- q1/8.
        endif

c     SMART

      elseif (advect.eq.4) then
        slp1 = 1.5
        slp2 = .5
        if (vel.gt.0d0) then
          qt1 = 3.*q2/8. + 3.*q3/4.- q4/8.
          qt2 = slp1*q3 + (1.-slp1)*q4
          qt3 = slp2*q3 + (1.-slp2)*q2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        else
          qt1 = 3.*q3/8. + 3.*q2/4.- q1/8.
          qt2 = slp1*q2 + (1.-slp1)*q1
          qt3 = slp2*q2 + (1.-slp2)*q3
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        endif

c     Smooth SMART

      elseif (advect.eq.5) then
        xp1 = .3
        xp2 = .5
        if (vel.gt.0d0) then
          qt1  = (q3-q4)/(q2-q4)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1.)
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 3.*q2/8. + 3.*q3/4.- q4/8.
          else
            inter = q3
          endif
        else
          qt1  = (q2-q1)/(q3-q1)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            c = 1.
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 3.*q3/8. + 3.*q2/4.- q1/8.
          else
            inter = q2
          endif
        endif

c     Centered high-order

      elseif (advect.eq.6) then
        inter = (-q4 + 9.*q3 + 9.*q2 - q1)/16.

      endif

c End

      return
      end

c fmed
c####################################################################
      real*8 function fmed(p1,p2,p3)
      implicit none                !For safe fortran
c--------------------------------------------------------------------
c    This function computes intermediate value of p1, p2, p3.
c--------------------------------------------------------------------

c Call variables

      real*8       p1,p2,p3

c Local variables

c Begin program

      fmed = min( max(p1,p2) , max( p3,min(p1,p2) ) )

c End

      return
      end

c tvd
c######################################################################
      subroutine tvd(i,j,nx,ny,dx1,dy1,vxx,vyy,arr
     .              ,phi_e,phi_w,phi_n,phi_s)
cc      implicit none        !For safe fortran
c----------------------------------------------------------------------
c     This function computes Van Leer's TVD flux limiter coefficients.
c----------------------------------------------------------------------

      use timeStepping

      implicit none

c Call variables

      integer*4    nx,ny,i,j
      real*8       dx1,dy1,arr(0:nx+1,0:ny+1)
      real*8       vxx(0:nx+1,0:ny+1),vyy(0:nx+1,0:ny+1)
      real*8       phi_e,phi_w,phi_n,phi_s

c Local variables

      real*8       vveln, vvels, uvele, uvelw
     .            ,phip_e, phip_w, phip_n, phip_s
     .            ,phim_e, phim_w, phim_n, phim_s
     .            ,k_e  , k_w  , k_n  , k_s

      real*8       qe,qw,qn,qs,vvel,uvel,dup,dum

c Begin program

c Calculate face velocities

      uvele = (vxx(i+1,j) + vxx(i  ,j))/2.
      uvelw = (vxx(i  ,j) + vxx(i-1,j))/2.
      vveln = (vyy(i,j+1) + vyy(i,j  ))/2.
      vvels = (vyy(i,j  ) + vyy(i,j-1))/2.

c Calculate coefficients

      k_n = .5*( dy1 - dt*abs(vveln) )*abs(vveln)
      k_s = .5*( dy1 - dt*abs(vvels) )*abs(vvels)
      k_e = .5*( dx1 - dt*abs(uvele) )*abs(uvele)
      k_w = .5*( dx1 - dt*abs(uvelw) )*abs(uvelw)

c Check error

      if (    (k_n.lt.0d0).or.(k_s.lt.0d0)
     .    .or.(k_e.lt.0d0).or.(k_w.lt.0d0) ) then
        write (*,*) 'Explicit method is unstable'
        stop
      endif

c East flux

      dup = (arr(i+1,j) - arr(i  ,j))
      dum = (arr(i  ,j) - arr(i-1,j))
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phip_e=k_e*(abs(dum) + dum*sign(1d0,dup))/(abs(dup)+abs(dum))

      dup = (arr(i+2,j) - arr(i+1,j))
      dum = (arr(i+1,j) - arr(i  ,j))
      if (i.eq.nxd) dup = (arr(3,j) - arr(i+1,j))
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phim_e=k_e*(abs(dup) + dup*sign(1d0,dum))/(abs(dup)+abs(dum))

      if (uvele.ne.0d0)
     .     phi_e = (uvele + abs(uvele))/2./abs(uvele)*phip_e
     .           - (uvele - abs(uvele))/2./abs(uvele)*phim_e

c West flux

      dup = (arr(i  ,j) - arr(i-1,j))
      dum = (arr(i-1,j) - arr(i-2,j))
      if (i.eq.1) dum = (arr(i-1,j) - arr(nx-2,j))
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phip_w=k_w*(abs(dum) + dum*sign(1d0,dup))/(abs(dup)+abs(dum))

      dup = (arr(i+1,j) - arr(i  ,j))
      dum = (arr(i  ,j) - arr(i-1,j))
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phim_w=k_w*(abs(dup) + dup*sign(1d0,dum))/(abs(dup)+abs(dum))

      if (uvelw.ne.0d0)
     .     phi_w = (uvelw + abs(uvelw))/2./abs(uvelw)*phip_w
     .           - (uvelw - abs(uvelw))/2./abs(uvelw)*phim_w

c North flux

      dup = (arr(i,j+1) - arr(i,j  ))
      dum = (arr(i  ,j) - arr(i,j-1))
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phip_n=k_n*(abs(dum) + dum*sign(1d0,dup))/(abs(dup)+abs(dum))

      dup = (arr(i,j+2) - arr(i,j+1))
      dum = (arr(i,j+1) - arr(i,j  ))
      if (j.eq.nyd) dup = dum
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phim_n=k_n*(abs(dup) + dup*sign(1d0,dum))/(abs(dup)+abs(dum))

      if (vveln.ne.0d0)
     .     phi_n = (vveln + abs(vveln))/2./abs(vveln)*phip_n
     .           - (vveln - abs(vveln))/2./abs(vveln)*phim_n

c South flux

      dup = (arr(i,j  ) - arr(i,j-1))
      dum = (arr(i,j-1) - arr(i,j-2))
      if (j.eq.1) dum = dup
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phip_s=k_s*(abs(dum) + dum*sign(1d0,dup))/(abs(dup)+abs(dum))

      dup = (arr(i,j+1) - arr(i,j  ))
      dum = (arr(i,j  ) - arr(i,j-1))
cc      if (j.eq.1) dum = dup
      if (dup.ne.0d0.or.dum.ne.0d0) 
     .     phim_s=k_s*(abs(dup) + dup*sign(1d0,dum))/(abs(dup)+abs(dum))

      if (vvels.ne.0d0)
     .     phi_s = (vvels + abs(vvels))/2./abs(vvels)*phip_s
     .           - (vvels - abs(vvels))/2./abs(vvels)*phim_s

c Check error

      if (    (phi_n.lt.0d0).or.(phi_s.lt.0d0)
     .    .or.(phi_e.lt.0d0).or.(phi_w.lt.0d0) ) then
        write (*,10) phi_n,phi_s,phi_e,phi_w,vveln,vvels,uvele,uvelw
cc        write (*,*) vveln,vvels,uvele,uvelw
        stop
      endif

c End program

      return

 10   format ('Flux correction is anti-diffusive',/,
     .        '  phi_n =',e10.3,'  phi_s =',e10.3,
     .        '  phi_e =',e10.3,'  phi_w =',e10.3,/,
     .        '  vel_n =',e10.3,'  vel_s =',e10.3,
     .        '  vel_e =',e10.3,'  vel_w =',e10.3 )
      end

c vecad2
c###################################################################
      subroutine vecad2(nx,cc1,x,cc2,y)
      implicit none             !For safe fortran
c-------------------------------------------------------------------
c     This subroutine computes
c
c     x = cc1*x + cc2*y
c
c     where x, y are vectors
c-------------------------------------------------------------------

c Call variables

      integer*4   nx
      real*8      x(nx),y(nx),cc1,cc2

c Local variables

      integer*4   i

c Begin program

      do i = 1,nx

        x(i) = cc1*x(i) + cc2*y(i)

      enddo

      return
      end subroutine

c vecadd
c###############################################################
      subroutine vecadd(nx,cc1,x,cc2,y,z)
      implicit none           !For safe fortran
c---------------------------------------------------------------
c     This subroutine computes
c
c     z = cc1*x + cc2*y
c
c     where x, y, and z are vectors
c---------------------------------------------------------------

c Call variables

      integer*4   nx
      real*8      x(nx),y(nx),z(nx),cc1,cc2

c Local variables

      integer*4   i

c Begin program

      do i = 1,nx
c
        z(i) = cc1*x(i) + cc2*y(i)

      enddo

      return
      end subroutine

c dot
c####################################################################
      double precision function dot(n,x,y)
      implicit none           !For safe fortran
c-------------------------------------------------------------------
c     Calculates dot product of vectors x, y
c-------------------------------------------------------------------

c Call variables

      integer*4   n           
      double precision :: x(n), y(n)

c Local variables

      integer*4   i
      double precision :: val

c Begin program

      val = 0d0
      do i=1,n
        val = val + x(i) * y(i)
      enddo

      dot = val

c End program

      return
      end function

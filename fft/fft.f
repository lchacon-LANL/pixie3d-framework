c perf_fft
c#######################################################################
      subroutine perf_fft(nn,psi,x,inter,ism)

c***********************************************************************
c     Preprocesses and calculates FFT of input 'psi', and gives output
c     in psir (real part) and psii (imaginary part).
c
c     Preprocessor assumes real signal; odd entries in data() are the 
c     actual time series, i.e. real part, even entries in data() are the 
c     imaginary part.
c
c     Taken from J. M. Finn.
c
c     In call:
c        + 'nn' is dimension of psi
c        + 'psi' is signal in real space
c        + 'ntot' is dimension of x
c        + 'x' is signal independent variable (i.e.,time)
c        + 'int' decides if interpolation is required (1).
c        + 'ism' decides if Hanning smoothing is performed (1).
c***********************************************************************

      implicit none

c Call variables

      integer*4     nn,inter,ism
      real*8        psi(nn),x(nn)

c Local variables

      real*8        dxx,fk_r(nn),fk_i(nn),omega(nn)
     .             ,ps(nn)
      integer*4     nx,i

c Begin program

c Initialize arrays

      do i = 1,nn
        fk_r(i) = 0d0
        fk_i(i) = 0d0
      enddo

c Calculate new length of vector (power of 2 for fft)

      nx = int(dlog(1d0*(nn-1))/dlog(2d0)+0.001) + 1
      if (    (abs(nn-2**nx).gt.abs(nn-2**(nx-1)))
     .    .or.(2**nx.gt.nn)) nx = nx -1

cc      nx = int(dlog(1d0*(nn-1))/dlog(2d0)+0.001)

      nx = 2**nx 

      write (*,*) 'Number of points provided   for FFT: ',nn
      write (*,*) 'Number of points considered for FFT: ',nx

c Copy and supplement time series up to nx

      do i = 1,nn
        ps(i) = psi(i)
      enddo

      dxx = x(nn)-x(nn-1)
      do i = nn + 1,nx
        x(i) = x(i-1) + dxx
        ps(i) = 0d0
      enddo

c Perform fft 

      call fft(nx,ps,x,fk_r,fk_i,omega,inter,ism)

c Find power spectrum

      do i = 1,nx
        ps(i) = log(fk_r(i)**2 + fk_i(i)**2)
      enddo

c Plot results

      open(unit = 20,file='fft.bin',form='unformatted')

cc      do i = 1,nx
cc        if (i.eq.nx/2+1) cycle
cc        write (20) real(omega(i)),real(fk_r(i)),real(fk_i(i))
cc     .            ,real(ps(i))
cc      enddo

      do i = nx/2+2,nx
        write (20) real(omega(i)),real(fk_r(i)),real(fk_i(i))
     .            ,real(ps(i))
      enddo
      do i = 1,nx/2
        write (20) real(omega(i)),real(fk_r(i)),real(fk_i(i))
     .            ,real(ps(i))
      enddo
      write (20)

      close (20)

c End program

      return
      end

c fft
c#######################################################################
      subroutine fft(nx,psi,x,nv,psir,psii,ak,inter,ism)
      implicit none     !For safe Fortran
c***********************************************************************
c     Preprocesses and calculates FFT of input 'psi', and gives output
c     in psir (real part) and psii (imaginary part).
c
c     Preprocessor assumes real signal; odd entries in data() are the 
c     actual time series, i.e. real part, even entries in data() are the 
c     imaginary part.
c
c     Taken from J. M. Finn.
c
c     In call:
c        + 'nx' is signal vector dimension
c        + 'psi' is signal in real space
c        + 'x' is signal independent variable
c        + 'nv' is fft vector dimension
c        + 'psir' is real part of transform
c        + 'psii' is imaginary part of transform
c        + 'ak' is new abcissae (k or omega)
c        + 'inter' decides if interpolation is required (1)
c        + 'ism' decides if smoothing is performed (1).
c***********************************************************************

c Call variables

      integer*4   nx,nv,inter,ism
      real*8      psi(nx),x(nx)
      real*8      psir(nv),psii(nv),ak(nv)

c Local variables

      integer*4   i
      real*8      pi,twopi,xmed,lx
      real*8      psi_wrk(nv),x_wrk(nv),data(2*nv)

c Begin program

      pi   =acos(-1d0)
      twopi=2*pi

      lx   =  x(nx)-x(1)
      xmed = (x(nx)+x(1))/2.

c Smooth initial data (Hanning)

      if (ism.eq.1) psi = cos(pi*(x-xmed)/lx)**2*psi

c Check for power of 2 in signal

cc      if (   ceiling(log(1d0*nx)/log(2d0)) 
cc     .    /= int    (log(1d0*nx)/log(2d0))) then

cc        nv = 2**ceiling(log(1d0*nx)/log(2d0))

cc        inter = 1
cc      else
cc        nv = nx
cc      endif

      if (   ceiling(log(1d0*nv)/log(2d0)) 
     .    /= int    (log(1d0*nv)/log(2d0))) then

        write (*,*) 'FFT vectors sizes are not a power of 2.'
        write (*,*) 'Aborting...'
        stop
      endif

c Interpolate time plots to provide uniform spacing

      if (inter.eq.1 .or. nx /= nv) then
        call IntDriver1d (nx,x,psi,nv,x_wrk,psi_wrk,3)
      else
        psi_wrk = psi
        x_wrk   = x
      endif

c Preprocessor for FFT

      do i=1,nv
         data(2*i-1)=psi_wrk(i)
         data(2*i)  =0d0
         if(i.le.nv/2+1) then
           ak(i)= twopi*(   i-1)/lx
         else
           ak(i)=-twopi*(nv+1-i)/lx
         endif
      enddo

c FFT

      call four1(data,nv,1)
      
      do i=1,nv
         psir(i)=data(2*i-1)/nv
         psii(i)=data(2*i  )/nv
      enddo

c End program

      return
      end

c ifft
c#######################################################################
      subroutine ifft(nx,psi,x,nv,psir,psii,ak,inter)
cc      implicit none     !For safe Fortran
c***********************************************************************
c     Performs inverse Fourier transform. Odd entries in data() are real 
c     parts of the FT; even entries are the imaginary parts of the FT.  
c     four1(.,.,-1) does the inverse transform giving the original signal 
c     in psi(). zero() should be identically zero, but I compute it to 
c     check (during debugging). 
c
c     Taken from J. M. Finn.
c
c     In call:
c        + 'nx' is signal vector dimension
c        + 'psi' is signal in real space
c        + 'x' is signal independent variable
c        + 'nv' is fft vector dimension
c        + 'psir' is real part of transform
c        + 'psii' is imaginary part of transform
c        + 'ak' is the abcissae (k or omega)
c        + 'inter' is dummy
c***********************************************************************

c Call variables

      implicit real*8(a-h,o-z)

      integer*4     nx,nv,inter

      real*8        psi(nx),x(nx)
      real*8        psir(nv),psii(nv),ak(nv)

c Local variables

      integer*4     im
      real*8        zero(nv),data(2*nv),psi_wrk(nv),x_wrk(nv),lx

c Begin program

      pi   =acos(-1d0)
      twopi=2*pi

      do i=1,nv
         data(2*i-1)=psir(i)
         data(2*i)  =psii(i)
      enddo
      
      call four1(data,nv,-1)

      do i=1,nv
         psi_wrk(i)=data(2*i-1)
         zero   (i)=data(2*i  )
      enddo

c Reconstruct uniform grid from ak

      lx = x(nx)-x(1)

      x_wrk(1) = x(1)
      do i=2,nv
         if(i.le.nv/2+1) then
           im = int(lx*ak(i)/twopi+1d-3)
         else
           im = nv+int(lx*ak(i)/twopi-1d-3)
         endif
        x_wrk(i) = x_wrk(1) + lx/(nv-1)*im
      enddo

c Interpolate ifft to original grid

      if (nx /= nv) then
        call IntDriver1d (nv,x_wrk,psi_wrk,nx,x,psi,3)
      else
        psi = psi_wrk
        x   = x_wrk
      endif

c End program

      return
      end

c four1
c########################################################################
      subroutine four1(data,nx,isign)
cc      implicit none       !For safe Fortran
c************************************************************************
c     Performs direct and inverse Fourier transform of 'data'.
c
c     Taken from J. M. Finn (Copied from Numerical Recipes.  
c     It's good for something.)
c
c     In call:
c        + 'data' is postprocessed signal
c        + 'nx' is signal vector length (must be power of 2)
c        + 'isign' determines if direct (1) or inverse (-1) FFT is to be
c          performed. 
c************************************************************************

c Variables

      implicit real*8(a-h,o-z)

      integer*4   nx
      real*8 data(2*nx)

c Begin program

      pi=4*atan(1.0)
      twopi=2*pi
      nn=nx
      n=2*nn
      j=1

      do i=1,n,2
         if(j.gt.i) then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=n/2
 1       if((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
            go to 1
         endif
         j=j+m
      enddo

      mmax=2
 2    if(n.gt.mmax) then
         istep=2*mmax
         theta=twopi/(isign*mmax)
         wpr=-2*sin(0.5*theta)**2
         wpi=sin(theta)
         wr=1
         wi=0
         do m=1,mmax,2
            do i=m,n,istep
               j=i+mmax
               tempr=wr*data(j)-wi*data(j+1)
               tempi=wr*data(j+1)+wi*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
            enddo 
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         enddo 
         mmax=istep
      go to 2
      endif

      return
      end

c resistiveWallBC
c###################################################################
      subroutine resistiveWallBC(nn,array)

c-------------------------------------------------------------------
c     Integrates vacuum equations to update resistive wall boundary.
c
c     Integrals are done in Fourier space, using fft's. Since the
c     periodic BC are imposed at nodes 1 and nn (which is a multiple
c     of 2), the fft requires to exclude one of the end points, which 
c     results in array sizes not multiple of 2. Hence, we interpolate 
c     the signal to an array of 1:nn+1 and take only points (1:nn).
c-------------------------------------------------------------------

      use mg_setup

      use variables

      use resistiveWall

cc      use geometry

      use timeStepping

      implicit none

c Call variables

      integer          :: nn
      double precision :: array(nn)

c Local variables

      integer          :: i,itp_order
      double precision :: x(nn),signal(nn)
      double precision :: x_wrk(nn+1),s_wrk(nn+1)
      double precision :: kn(nn),fk_r(nn),fk_i(nn)
      double complex   :: ank(nn),bnk(nn),anp(nn),anp2k(nn)
      double complex   :: ankp(nn),anpk(nn),sn(nn),cn(nn)

      double precision :: twopi,dyy,tauinv,ko

c Begin program

      twopi = 2*acos(-1d0)

      itp_order = 3

c Setup independent variable

      x = xl(1:nn,ngrd)

      dyy = dy(ngrd)

c Interpolation test

cc      signal = u_n%array_var(3)%array(1:nn,nyd+1)
cccc      signal = 100*x**4 -x + 1d0
cc
cc      call IntForFFT (nn,x,signal,nn+1,x_wrk,s_wrk,itp_order)
cc
cc      call IntForFFT (nn+1,x_wrk,s_wrk,nn,x,array,itp_order)
cc
cccc      where (array /= 0d0)
cccc        array = (array - signal)/array*100
cccc      elsewhere
cccc        array = 0d0
cccc      end where
cc
cc      write (*,*) signal
cc      write (*,*) array
cc      stop

c fft old psi at resistive wall

      signal = u_n%array_var(3)%array(1:nn,nyd+1)

      call IntForFFT (nn,x,signal,nn+1,x_wrk,s_wrk,itp_order)

      call fft(nn,s_wrk(1:nn),x_wrk(1:nn),nn,fk_r,fk_i,kn,0,0)

      ank = cmplx(fk_r,fk_i)

c Define fundamental mode

      ko = abs(kn(2))

c fft old psi at inner grid point in plasma

      signal = u_n%array_var(3)%array(1:nn,nyd)

      call IntForFFT (nn,x,signal,nn+1,x_wrk,s_wrk,itp_order)

      call fft(nn,s_wrk(1:nn),x_wrk(1:nn),nn,fk_r,fk_i,kn,0,0)

      anpk = cmplx(fk_r,fk_i)

c fft old psi at second inner grid point in plasma

      signal = u_n%array_var(3)%array(1:nn,nyd-1)

      call IntForFFT (nn,x,signal,nn+1,x_wrk,s_wrk,itp_order)

      call fft(nn,s_wrk(1:nn),x_wrk(1:nn),nn,fk_r,fk_i,kn,0,0)

      anp2k = cmplx(fk_r,fk_i)

c fft old psi at vacuum wall

      signal = vacuumWallSignal(x,time)

      call IntForFFT (nn,x,signal,nn+1,x_wrk,s_wrk,itp_order)

      call fft(nn,s_wrk(1:nn),x_wrk(1:nn),nn,fk_r,fk_i,kn,0,0)

      bnk = cmplx(fk_r,fk_i)

c Proportional feedback

      !Define gain

      if ( (gi /= 0d0 .or. gr /= 0d0) .and. gain == 0d0) then
        gain = sqrt(gi**2 + gr**2)
        shift = atan(gi/gr)/ko
      endif

      if (allmodes) then

      !Feedback on all modes

        where (kn /= 0d0)
          bnk = bnk - gain*cmplx(cos(kn*shift),sin(kn*shift))*ank
        end where

      else

      !Feedback on longest-wavelength mode ko

        i = 2
        bnk(i) = bnk(i)
     .         - gain*cmplx(cos(kn(i)*shift),sin(kn(i)*shift))*ank(i)

        i = nn
        bnk(i) = bnk(i)
     .         - gain*cmplx(cos(kn(i)*shift),sin(kn(i)*shift))*ank(i)
      
      endif

c Perform time integral of resistive wall condition

      tauinv = 1d0/tau_wall

      where (kn == 0d0)
        cn = tauinv/(vw-1)
        sn = tauinv*(bnk/(vw-1) - (3*ank-4*anpk+anp2k)/2d0/dyy)
cc        sn = tauinv*(bnk/(vw-1) - (ank-anpk)/dyy)
      elsewhere
cc        cn = tauinv*kn*(cosh(kn*(vw-1))
cc     .                 +gain*cmplx(cos(kn*shift),sin(kn*shift)))
cc     .                 /sinh(kn*(vw-1))
        cn = tauinv*kn*cosh(kn*(vw-1))/sinh(kn*(vw-1))
        sn = tauinv*(kn*bnk/sinh(kn*(vw-1))
     .              - (3*ank-4*anpk+anp2k)/2d0/dyy )
cc        sn = tauinv*(kn*bnk/sinh(kn*(vw-1)) - (ank-anpk)/dyy)
      end where

      where (cn == 0d0)
        ankp = ank + sn*dt
      elsewhere
        ankp = exp(-cn*dt)*ank + sn*(1.-exp(-cn*dt))/cn
      end where

c Perform inverse fft

      fk_r = real(ankp)
      fk_i = imag(ankp)

      call ifft(nn,s_wrk(1:nn),x_wrk(1:nn),nn,fk_r,fk_i,kn,0)

      s_wrk(nn+1) = s_wrk(1)  !Impose periodic boundary condition

      call IntForFFT (nn+1,x_wrk,s_wrk,nn,x,array,itp_order)

c End program

      contains

c     vacuumWallSignal
c     ##############################################################
      function vacuumWallSignal(x,time) result (signal)

        implicit none

c     Call variables

        double precision :: time
        double precision,dimension(:),intent(in) :: x
        double precision,dimension (size(x))     :: signal

c     Local variables

        double precision :: kk,psiw

c     Begin program

c     Calculate flux at vacuum wall to guarantee an equilibrium

cc        psiw = (u_0%array_var(3)%array(nxd/2,nyd+1)
cc     .         -u_0%array_var(3)%array(nxd/2,nyd  ))/dyy*(vw-1)

        psiw = (3*u_0%array_var(3)%array(nxd/2,nyd+1)
     .         -4*u_0%array_var(3)%array(nxd/2,nyd  )
     .         +  u_0%array_var(3)%array(nxd/2,nyd-1))/2d0/dyy*(vw-1)

c     Initialize vacuum wall signal

        kk    = ko*nhvw

        phase = w_vw*time

c     Calculate signal

        signal = psiw + eps*cos(kk*x-phase)

c     End program

      end function vacuumWallSignal
      
      end subroutine resistiveWallBC

c resistiveWallFlow
c###################################################################
      subroutine resistiveWallFlow(nn,array)

c-------------------------------------------------------------------
c     Integrates vacuum equations to update resistive wall boundary
c     Integrals are done in Fourier space, using fft's.
c-------------------------------------------------------------------

      use mg_setup

      use variables

      use resistiveWall

      use timeStepping

      implicit none

c Call variables

      integer          :: nn
      double precision :: array(nn)

c Local variables

      integer          :: i

      double precision :: dxx,dyy,ss

c Begin program

      dxx = dx(ngrd)
      dyy = dy(ngrd)

c Set Maxwell stress source

cc      ss = sum(bx(1:nxd-1,nyd+1)*by(1:nxd-1,nyd+1),1)/(nxd-1)

      ss = maxstrs

cc      write (*,*) 'Maxwell stress tensor ',ss

c Perform time integral of <ux>

      uxk  = -u_n%array_var(1)%array(nxd/2,nyd+1)

      if (drag.ne.0d0) then
        uxkp = u00 + exp(-drag*dt)*(uxk-u00) +ss*(1.-exp(-drag*dt))/drag
cc        uxkp = uxk + (ss-drag*(uxk-u00))*dt
cc        uxkp = (uxk + (ss+drag*u00)*dt)/(1 + drag*dt)
      else
        uxkp = uxk + ss*dt
      endif

c Set output array

      array = - uxkp

c End program

      end subroutine resistiveWallFlow

c IntForFFT
c#######################################################################
      subroutine IntForFFT (nx,x,vec,nv,x1,vec1,order)
c***********************************************************************
c     Interpolates vector vec along coordinate x to vector vec1 along x1
c     using uniform intervals for FFT, with interpolation order "order". 
c     Currently, we consider first order (order=1), second order (2)
c     and cubic splines (3). 
c***********************************************************************

      use oned_int

      implicit none       !For safe Fortran

c Call variables

      integer*4      nx,nv,order
      real*8         x(nx),vec(nx),x1(nv),vec1(nv)

c Local variables

      real*8         dxx
      integer*4      i,j

c Begin program

c Define new uniform coordinates for FFT

      dxx = (x(nx)-x(1))/(nv-1)

      x1(1) = x(1)
      do i = 2,nv-1
        x1(i) = x1(i-1) + dxx
      enddo
      x1(nv) = x(nx)

cc      if ( abs(x1(nv)-x(nx)).gt.1d-3) then
cc        write (*,*) 'Something is wrong with x1; x1(nv) - x(nx) ='
cc     .              ,x1(nv)-x(nx)
cc        stop
cc      endif

c Interpolate vector on new coordinates

      call IntDriver1d (nx,x,vec,nv,x1,vec1,order)

c End program

      return
      end

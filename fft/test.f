      program FastFourierTransform

      call test(60,64)

      end program

c test
c#######################################################################
      subroutine test(ftime,ftime1)
      implicit none         !For safe Fortran
c***********************************************************************
c     Test for FFT and surface of section using analytical functions.
c***********************************************************************

c Local variables

      integer*4       ftime,ftime1,i
      real*8          tmax,tmin,dxx,time(ftime),tst(ftime),dummy(ftime)
      real*8          fk_r(ftime1),fk_i(ftime1),pi,omega(ftime1),f0

c Begin program

      pi = acos(-1d0)

c Define sampling interval

      !This should be dxx=(tmax-tmin)/(ftime-1), but then it turns out
      !that the first and last value are equal if T=32 is the exact
      !period of a sinusoidal, and that screws up the fft. The following
      !definition states that there should be ftime values sampling the
      !sinusoidal without repeating the first and the last.
      tmax  = 16d0
      tmin  = -16d0
      dxx = (tmax-tmin)/ftime

c Reconstruct independent variable

      time(1) = tmin
      do i = 2,ftime,2
cc        time(i  ) = time(i-1) + 0.5*dxx
cc        time(i+1) = time(i  ) + 1.5*dxx

        time(i) = time(i-1) + dxx
        time(i+1) = time(i) + dxx
      enddo

c Create signal test

      f0 = 1./9.143
      do i = 1,ftime
cc        tst(i) = 2.*f0*sin(2*pi*f0*time(i))/(2*pi*f0*time(i))
cc        if (time(i).eq.0d0) tst(i) = 2*f0

        tst(i) = sin(2.*pi*time(i)/8.)

cc        tst(i) = f0*exp(-f0*abs(time(i)))
cc        tst(1) = .5*f0

cc        tst(i) = 0d0
cc        if (abs(time(i)).lt.5.) tst(i) = 1d0

cc        tst(i) = cos(2*pi*f0*time(i))
      enddo

c FFT test

      dummy = tst

      call fft(ftime,dummy,time,ftime1,fk_r,fk_i,omega,0,0)

      open(unit = 20,file='fft.bin',form='unformatted')

      do i = ftime1/2+2,ftime1
        write (20) real(omega(i)),real(fk_r(i)),real(fk_i(i))
cc        write (*,*) real(omega(i)),real(fk_r(i)),real(fk_i(i))
      enddo
      do i = 1,ftime1/2+1
        write (20) real(omega(i)),real(fk_r(i)),real(fk_i(i))
cc        write (*,*) real(omega(i)),real(fk_r(i)),real(fk_i(i))
      enddo

      close (20)

c Inverse FFT test

      dummy = 0d0

      call ifft(ftime,dummy,time,ftime1,fk_r,fk_i,omega,0)

      open(unit = 20,file='ifft.bin',form='unformatted')

      do i = 1,ftime
        write (20) real(time(i)),real(tst(i)),real(dummy(i))
cc        write (*,*) real(time(i)),real(tst(i)),real(dummy(i))
      enddo

      close (20)

c End program

      return
      end


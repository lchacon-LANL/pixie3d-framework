c-----------------------------------------------------------------------
c demonstration program for the lsode package.
c this is the version of 13 august, 1981.
c
c this version is in double precision.
c
c for computer systems requiring a program card, the following (with
c the c in column 1 removed) may be used..
c     program lsdem(lsout,tape6=lsout)
c
c the package is used to solve two simple problems,
c one with a full jacobian, the other with a banded jacobian,
c with all 8 of the appropriate values of mf in each case.
c if the errors are too large, or other difficulty occurs,
c a warning message is printed.  all output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1, f2, jac2
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter, miter1,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, er, erm, ero, hu, rtol, rwork, t,
     1   tout, tout1, y
      dimension y(25), rwork(697), iwork(45)
      data lout/6/, tout1/1.39283880203d0/, dtout/2.214773875d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-6
      lrw = 697
      liw = 45
      iopt = 0
c
c first problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(1h1/1x,40h demonstration program for lsode package////
     1  1x,39h problem 1..   van der pol oscillator../
     2  1x,38h  xdotdot - 3*(1 - x**2)*xdot + x = 0, ,
     3  24h   x(0) = 2, xdot(0) = 0/
     4  1x,6h neq =,i2/
     5  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1//)
c
      do 195 meth = 1,2
      do 190 miter1 = 1,4
      miter = miter1 - 1
      mf = 10*meth + miter
      write (lout,120) mf
 120  format(////1x,5h mf =,i3///
     1  6x,1ht,15x,1hx,15x,4hxdot,7x,2hnq,6x,1hh//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      tout = tout1
      ero = 0.0d0
      do 170 iout = 1,nout
        call lsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac1,mf)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,140) t,y(1),y(2),nqu,hu
 140    format(1x,e15.5,e16.5,e14.3,i5,e14.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 170
        er = dabs(y(1))/atol
        ero = dmax1(ero,er)
        if (er .lt. 1000.0d0) go to 170
        write (lout,150)
 150    format(//1x,41h warning.. error exceeds 1000 * tolerance//)
        nerr = nerr + 1
 170    tout = tout + dtout
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - neq*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//1x,32h final statistics for this run../
     1  1x,13h rwork size =,i4,15h   iwork size =,i4/
     2  1x,18h number of steps =,i5/
     3  1x,18h number of f-s   =,i5/
     4  1x,18h (excluding j-s) =,i5/
     5  1x,18h number of j-s   =,i5/
     6  1x,16h error overrun =,e10.2)
 190  continue
 195  continue
c
c second problem
c
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      nout = 5
      write (lout,210) neq,ml,mu,itol,rtol,atol
 210  format(1h1/1x,40h demonstration program for lsode package////
     1  1x,33h problem 2.. ydot = a * y , where,
     2  39h  a is a banded lower triangular matrix/
     3  1x,32h  derived from 2-d advection pde/
     4  1x,6h neq =,i3,7h   ml =,i2,7h   mu =,i2/
     5  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1//)
      do 295 meth = 1,2
      do 290 miter1 = 1,6
      miter = miter1 - 1
      if (miter .eq. 1 .or. miter .eq. 2) go to 290
      mf = 10*meth + miter
      write (lout,220) mf
 220  format(////1x,5h mf =,i3///
     1  6x,1ht,13x,8hmax.err.,5x,2hnq,6x,1hh//)
      t = 0.0d0
      do 230 i = 2,neq
 230    y(i) = 0.0d0
      y(1) = 1.0d0
      itask = 1
      istate = 1
      tout = 0.01d0
      ero = 0.0d0
      do 270 iout = 1,nout
        call lsode(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,mf)
        call edit2(y,t,erm)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,240) t,erm,nqu,hu
 240    format(1x,e15.5,e14.3,i5,e14.3)
        if (istate .lt. 0) go to 275
        er = erm/atol
        ero = dmax1(ero,er)
        if (er .le. 1000.0d0) go to 270
        write (lout,150)
        nerr = nerr + 1
 270    tout = tout*10.0d0
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 5) nfea = nfe - mband*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 290  continue
 295  continue
      write (lout,300) nerr
 300  format(////1x,31h number of errors encountered =,i3)
      stop
      end
      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(2), ydot(2)
      ydot(1) = y(2)
      ydot(2) = 3.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end
      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(2), pd(nrowpd,2)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
      return
      end
      subroutine f2 (neq, t, y, ydot)
      integer neq, i, j, k, ng
      double precision t, y, ydot, alph1, alph2, d
      dimension y(1), ydot(1)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end
      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(1), pd(nrowpd,1)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end
      subroutine edit2 (y, t, erm)
      integer i, j, k, ng
      double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
      dimension y(25)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      erm = 0.0d0
      if (t .eq. 0.0d0) return
      ex = 0.0d0
      if (t .le. 30.0d0) ex = dexp(-2.0d0*t)
      a2 = 1.0d0
      do 60 j = 1,ng
        a1 = 1.0d0
        do 50 i = 1,ng
          k = i + (j - 1)*ng
          yt = t**(i+j-2)*ex*a1*a2
          er = dabs(y(k)-yt)
          erm = dmax1(erm,er)
          a1 = a1*alph1/dfloat(i)
 50       continue
        a2 = a2*alph2/dfloat(j)
 60     continue
      return
      end

c-----------------------------------------------------------------------
c example problem.
c
c the following is a simple example problem, with the coding
c needed for its solution by lsode.  the problem is from chemical
c kinetics, and consists of the following three rate equations..
c     dy1/dt = -.04*y1 + 1.e4*y2*y3
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c     dy3/dt = 3.e7*y2**2
c on the interval from t = 0.0 to t = 4.e10, with initial conditions
c y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
c
c the following coding solves this problem with lsode, using mf = 21
c and printing results at t = .4, 4., ..., 4.e10.  it uses
c itol = 2 and atol much smaller for y2 than y1 or y3 because
c y2 has much smaller values.
c at the end of the run, statistical quantities of interest are
c printed (see optional outputs in the full description below).
c
      external fex, jex
      double precision atol, rtol, rwork, t, tout, y
      dimension y(3), atol(3), rwork(58), iwork(23)
      neq = 3
      y(1) = 1.d0
      y(2) = 0.d0
      y(3) = 0.d0
      t = 0.d0
      tout = .4d0
      itol = 2
      rtol = 1.d-4
      atol(1) = 1.d-6
      atol(2) = 1.d-10
      atol(3) = 1.d-6
      itask = 1
      istate = 1
      iopt = 0
      lrw = 58
      liw = 23
      mf = 21
      do 40 iout = 1,12
        call lsode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jex,mf)
        write(6,20)t,y(1),y(2),y(3)
  20    format(7h at t =,1pe12.4,6h   y =,1p3e14.6)
        if (istate .lt. 0) go to 80
  40    tout = tout*10.d0
      write(6,60)iwork(11),iwork(12),iwork(13)
  60  format(/12h no. steps =,i4,11h  no. f-s =,i4,11h  no. j-s =,i4)
      stop
  80  write(6,90)istate
  90  format(///22h error halt.. istate =,i3)
      stop
      end
 
      subroutine fex (neq, t, y, ydot)
      double precision t, y, ydot
      dimension y(3), ydot(3)
      ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
      ydot(3) = 3.d7*y(2)*y(2)
      ydot(2) = -ydot(1) - ydot(3)
      return
      end
 
      subroutine jex (neq, t, y, ml, mu, pd, nrpd)
      double precision pd, t, y
      dimension y(3), pd(nrpd,3)
      pd(1,1) = -.04d0
      pd(1,2) = 1.d4*y(3)
      pd(1,3) = 1.d4*y(2)
      pd(2,1) = .04d0
      pd(2,3) = -pd(1,3)
      pd(3,2) = 6.d7*y(2)
      pd(2,2) = -pd(1,2) - pd(3,2)
      return
      end   

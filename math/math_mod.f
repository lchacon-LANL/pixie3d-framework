      module math

      INTERFACE  solve_quadratic
        module procedure solve_quadratic_real,solve_quadratic_cmplx
      END INTERFACE

      INTERFACE  solve_cubic
        module procedure solve_cubic_real,solve_cubic_cmplx
      END INTERFACE

      real(8) :: pi

cc      INTERFACE determ
cc        module procedure determ3
cc      END INTERFACE

      INTERFACE fmed
        module procedure fmed_scl,fmed_vec
      END INTERFACE

      contains

c     factorial
c     ###############################################################
      recursive function factorial(n,m) result(fac)

      implicit none

c     ---------------------------------------------------------------
c     Finds n!/m!
c     ---------------------------------------------------------------

c     Call variables

      real(8) :: fac
      integer :: n,m

c     Local variables

c     Begin program

      if (n == m) then
        fac = 1d0
        return
      elseif (n > m) then
        fac = n*factorial(n-1,m)
      else
        fac = factorial(n,m-1)/m
      endif

      end function factorial

c     findRoundOff
c     ###############################################################
      subroutine findRoundOff (eps)

      implicit none

c     ---------------------------------------------------------------
c     Finds machine round-off constant epsmac
c     ---------------------------------------------------------------

c     Call variables

      real(8) :: eps

c     Local variables

      real(8) :: mag,mag2

      mag = 1d0
      do
        eps = mag
        mag = mag/2
        mag2 = 1d0 + mag
        if (.not.(mag2 > 1d0)) exit
      enddo

      end subroutine findRoundOff

c     determ3
c     #################################################################
      function determ3(tensor)

      real(8) :: tensor(3,3)

      real(8) :: determ3

      determ3 = tensor(1,1)*tensor(2,2)*tensor(3,3)
     .         +tensor(3,2)*tensor(2,1)*tensor(1,3)
     .         +tensor(1,2)*tensor(2,3)*tensor(3,1)
     .         -tensor(1,3)*tensor(2,2)*tensor(3,1)
     .         -tensor(1,1)*tensor(2,3)*tensor(3,2)
     .         -tensor(3,3)*tensor(1,2)*tensor(2,1)

      end function determ3

c     atanh
c     #################################################################
      real(8) function atanh(x)

        real(8) :: x

        atanh = 0.5*(log( (1+x)/(1-x) ) )

      end function atanh

c     acosh
c     #################################################################
      real(8) function acosh(x)

        real(8) :: x

        acosh = log(x+sqrt(x**2-1))

      end function acosh

c     bessel
c     ###################################################################
      SUBROUTINE BESSEL (N,X,BJ,BY)
c     -------------------------------------------------------------------
c     !!!!!!!!!!!!!!!!!!!!!!   Program 5.5   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     !                                                                 !
c     !Please Note:                                                     !
c     !                                                                 !
c     !(1) This computer program is written by Tao Pang in conjunction  ! 
c     !    with his book, "An Introduction to Computational Physics,"   !
c     !    published by Cambridge University Press in 1997.             !
c     !                                                                 !
c     !(2) No warranties, express or implied, are made for this program.!
c     !                                                                 !
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     Subroutine to generate J_n(x) and Y_n(x) with given x and
c     up to N=NMAX-NTEL.
c
c     -------------------------------------------------------------------

c     Call variables

      INTEGER, INTENT (IN) :: N
      REAL(8), INTENT (IN) :: X
      REAL(8), INTENT (OUT), DIMENSION (0:N) :: BJ,BY

c     Local variables

      INTEGER, PARAMETER :: NMAX=30,NTEL=5
      INTEGER :: I,J,K
      REAL(8) :: PI,GAMMA,SUM,SUM1
      REAL(8), DIMENSION (0:NMAX) :: B1

c     Begin program

      PI = 4d0*ATAN(1d0)
      GAMMA = 0.577215664901532860606512

      B1(NMAX)   = 0d0
      B1(NMAX-1) = 1d0

c     Generating J_n(x)

      SUM = 0.0
      DO I = NMAX-1, 1, -1 
        B1(I-1) = 2*I*B1(I)/X-B1(I+1)
        IF (MOD(I,2).EQ.0) SUM = SUM+2d0*B1(I)
      END DO
      SUM = SUM+B1(0)

      DO I = 0, N
        BJ(I) = B1(I)/SUM
      END DO

c     Generating Y_n(x) starts here

      SUM1 = 0.0
      DO K = 1, NMAX/2
        SUM1 = SUM1+(-1)**K*B1(2*K)/K
      END DO

      SUM1 = -4.0*SUM1/(PI*SUM)
      BY(0) = 2.0*(LOG(0.5d0*X)+GAMMA)*BJ(0)/PI+SUM1
      BY(1) = (BJ(1)*BY(0)-2d0/(PI*X))/BJ(0)
      DO I  = 1, N-1
        BY(I+1) = 2*I*BY(I)/X-BY(I-1)
      END DO

      END SUBROUTINE BESSEL

c     solve_quadratic_real
c     ##########################################################
      function solve_quadratic_real(a,b,c) result(root)

c     ----------------------------------------------------------
c     Solves for maximum of roots of quadratic polynomial:
c        a x^2 + b x + c = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: a,b,c
      complex(8) :: root(2)

c     Local variables

c     Begin program

      if (a == 0d0) then

        root(1) = -c/b
        root(2) = 0d0

      else

        root(1) = 0.5*(-b + sqrt(b**2-4*a*c))/a
        root(2) = 0.5*(-b - sqrt(b**2-4*a*c))/a

      endif

c     End program

      end function solve_quadratic_real

c     solve_quadratic_cmplx
c     ##########################################################
      function solve_quadratic_cmplx(a,b,c) result(root)

c     ----------------------------------------------------------
c     Solves for maximum of roots of quadratic polynomial:
c        a x^2 + b x + c = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      complex(8) :: root(2),a,b,c

c     Local variables

c     Begin program

      if (a == cmplx(0d0,0d0)) then

        root(1) = -c/b
        root(2) = 0d0

      else

        root(1) = (-b + sqrt(b**2-4*a*c))/2./a
        root(2) = (-b - sqrt(b**2-4*a*c))/2./a

      endif

c     End program

      end function solve_quadratic_cmplx

c     solve_cubic_real
c     ##########################################################
      function solve_cubic_real(a,b,c,d) result(root)

c     ----------------------------------------------------------
c     Solves for  roots of cubic polynomial:
c        a x^3 + b x^2 + c x + d = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: a,b,c,d
      complex(8) :: root(3)

c     Local variables

      integer :: iroot,inewt
      complex(8) :: S,T,R,Q,Det,a2,a1,a0
     .             ,res,res0,x,dx

c     Begin program

c     Special case of a=0

      if (a == 0d0) then
        root(1:2) = solve_quadratic(b,c,d)
        root(3) = 0d0
        return
      endif

c     Define auxiliary quantities
      
      a2=b/a
      a1=c/a
      a0=d/a
      Q=(3d0*a1-a2**2)/9d0
      R=(9d0*a2*a1-27d0*a0-2d0*a2**3)/54d0
      Det=exp(2d0*log(R)-3d0*log(Q))
      Det=exp(0.5d0*(3d0*log(Q)+log(1+Det)))
      S=(R+Det)**(1d0/3d0)
      T=(R-Det)**(1d0/3d0)
cc      S=S**(1d0/3d0)
cc      T=T**(1d0/3d0)

c     Find roots

      root(1)=-a2/3d0+(S+T)
      root(2)=-a2/3d0-(S+T-sqrt(3d0)*(0d0,1d0)*(S-T))/2d0
      root(3)=-a2/3d0-(S+T+sqrt(3d0)*(0d0,1d0)*(S-T))/2d0

c     Converge with Newton for accuracy

      do iroot=1,3
        x = root(iroot)

        res0 = poly(x,a,b,c,d,0)
        res = res0

cc        write (*,*) 'Residual root',iroot,'=',res/res0

        do inewt=1,100
          dx = -res/poly(x,a,b,c,d,1)

          !Globalization
          do while (abs(poly(x+dx,a,b,c,d,0)) > abs(res)) 
            dx = 0.5*dx
cc            write (*,*) 'Update root=',dx
          enddo

          x = x + dx

          res = poly(x,a,b,c,d,0)

cc          write (*,*) 'Residual root',iroot,'=',res/res0
cc          write (*,*) 'Update root  ',iroot,'=',dx/x

          if (abs(res/res0) < 1d-8 .or. abs(dx/x) < 1d-10) exit
        enddo
        root(iroot) = x
cc        write (*,*)
      enddo

cc      write (*,*) 'Solution roots=',root
cc      write (*,*)

c     End program

      contains

      function poly(x,a,b,c,d,order)

      integer    :: order
      real(8)    :: a,b,c,d
      complex(8) :: x,poly

      select case(order)
      case(0)
        poly = a*x**3 + b*x**2 + c*x + d
      case(1)
        poly = 3*a*x**2 + 2*b*x + c
      end select
      
      end function poly

      end function solve_cubic_real

c     solve_cubic_cmplx
c     ##########################################################
      function solve_cubic_cmplx(a,b,c,d) result(root)

c     ----------------------------------------------------------
c     Solves for  roots of cubic polynomial:
c        a x^3 + b x^2 + c x + d = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      complex(8) :: root(3),a,b,c,d

c     Local variables

      integer :: iroot,inewt
      complex(8) :: S,T,R,Q,Det,a2,a1,a0,res,res0,x,dx

c     Begin program

c     Special case of a=0

      if (a == cmplx(0d0,0d0)) then
        root(1:2) = solve_quadratic(b,c,d)
        root(3) = 0d0
        return
      endif

c     Define auxiliary quantities

      a2=b/a
      a1=c/a
      a0=d/a
      Q=(3d0*a1-a2**2)/9d0
      R=(9d0*a2*a1-27d0*a0-2d0*a2**3)/54d0
      Det=exp(2d0*log(R)-3d0*log(Q))
      Det=exp(0.5d0*(3d0*log(Q)+log(1+Det)))
      S=(R+Det)**(1d0/3d0)
      T=(R-Det)**(1d0/3d0)

cc      a2=b/a
cc      a1=c/a
cc      a0=d/a
cc      Q=(3*a1-a2**2)/9d0
cc      R=(9*a2*a1-27*a0-2*a2**3)/54d0
cc      Det=Q**3+R**2
cc      S=R+sqrt(Det)
cc      T=R-sqrt(Det)
cc
ccc     Need to make sure that (-8)^(1/3)=-2
cc
cc      if (aImag(S)==0 .and. Real(S)<0) then
cc         S=-(abs(S))**(1d0/3d0)
cc      else
cc         S=S**(1d0/3d0)
cc      endif
cc      if (aImag(T)==0 .and. Real(T)<0) then
cc         T=-(abs(T))**(1d0/3d0)
cc      else
cc         T=T**(1d0/3d0)
cc      endif

c     Find roots

      root(1)=-a2/3d0+(S+T)
      root(2)=-a2/3d0-(S+T-sqrt(3d0)*(0d0,1d0)*(S-T))/2d0
      root(3)=-a2/3d0-(S+T+sqrt(3d0)*(0d0,1d0)*(S-T))/2d0

c     Converge with Newton for accuracy

cc      do iroot=1,3
cc        x = root(iroot)
cc        res0 = a*x**3 + b*x**2 + c*x + d
cc        res = res0
cc        do inewt=1,10
cc          dx = -res/(3*a*x**2 + 2*b*x + c)
cc          x = x + dx
cc          write (*,*) 'Residual root',iroot,'=',res/res0
cc          write (*,*) 'Update root',iroot,'=',dx/x
cc          if (abs(res/res0) < 1d-8 .or. abs(dx/x) < 1d-10) exit
cc          res = a*x**3 + b*x**2 + c*x + d
cc        enddo
cc        root(iroot) = x
cc      enddo

      do iroot=1,3
        x = root(iroot)

        res0 = poly(x,a,b,c,d,0)
        res = res0

cc        write (*,*) 'Residual root',iroot,'=',res/res0

        do inewt=1,100
          dx = -res/poly(x,a,b,c,d,1)

          !Globalization
          do while (abs(poly(x+dx,a,b,c,d,0)) > abs(res)) 
            dx = 0.5*dx
cc            write (*,*) 'Update root=',dx
          enddo

          x = x + dx

          res = poly(x,a,b,c,d,0)

cc          write (*,*) 'Residual root',iroot,'=',res/res0
cc          write (*,*) 'Update root  ',iroot,'=',dx/x

          if (abs(res/res0) < 1d-8 .or. abs(dx/x) < 1d-10) exit
        enddo
        root(iroot) = x
cc        write (*,*)
      enddo

cc      write (*,*)

c     End program

      contains

      function poly(x,a,b,c,d,order)

      integer    :: order
      complex(8) :: a,b,c,d
      complex(8) :: x,poly

      select case(order)
      case(0)
        poly = a*x**3 + b*x**2 + c*x + d
      case(1)
        poly = 3*a*x**2 + 2*b*x + c
      end select
      
      end function poly

      end function solve_cubic_cmplx

c     int2char
c     ################################################################
      function int2char(n) result (chr)

      implicit none

      integer      :: n
      character(10):: chr

      integer      :: i,exp,k,j
      character(3) :: c

      if (abs(n) > 0) then
         exp = int(log(float(n))/log(1d1))
      else
         exp = 0
      endif

      if (n >= 0) then
        chr=''
      else
        chr='-'
      endif

      k = abs(n)
      do i=exp,0,-1
         j = k/10**i
         c = achar(48+j)
         chr = trim(chr)//trim(c)
         if (i > 0 .and. j /= 0) k = mod(k,j*10**i)
      enddo

      end function int2char

c     fmed_scl
c     ###############################################################
      function fmed_scl(p1,p2,p3) result(fmed)
      implicit none                !For safe fortran
c     ---------------------------------------------------------------
c     This function computes intermediate value of p1, p2, p3.
c     ---------------------------------------------------------------

c     Call variables

      real(8) :: p1,p2,p3,fmed

c     Local variables

c     Begin program

      fmed = min( max(p1,p2) , max( p3,min(p1,p2) ) )

c     End

      end function fmed_scl

c     fmed_vec
c     ###############################################################
      function fmed_vec(p1,p2,p3) result(fmed)
      implicit none                !For safe fortran
c     ---------------------------------------------------------------
c     This function computes intermediate value of p1, p2, p3.
c     ---------------------------------------------------------------

c     Call variables

      real(8) :: p1(:),p2(:),p3(:),fmed(size(p1))

c     Local variables

c     Begin program

      fmed = min( max(p1,p2) , max( p3,min(p1,p2) ) )

c     End

      end function fmed_vec

      end module math

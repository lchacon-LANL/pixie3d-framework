      module math

      contains

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

c     solve_quadratic
c     ##########################################################
      function solve_quadratic(a,b,c) result(root)

c     ----------------------------------------------------------
c     Solves for maximum of roots of cubic polynomial:
c        a + b x +c x^2 = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)   :: root(2),a,b,c

c     Local variables

c     Begin program

      root(1) = (-b + sqrt(b**2-4*a*c))/2./a

      root(2) = (-b - sqrt(b**2-4*a*c))/2./a

c     End program

      end function solve_quadratic

c     solve_cubic
c     ##########################################################
      function solve_cubic(a,b,c,d) result(root)

c     ----------------------------------------------------------
c     Solves for maximum of roots of cubic polynomial:
c        a + b x +c x^2 + d x^3 = 0
c     ----------------------------------------------------------

      implicit none

c     Call variables

      real(8)   :: root(3),a,b,c,d

c     Local variables

c     Begin program

      root(1) =
     -     -c/(3.*d) - (2**0.3333333333333333*(-c**2 + 3*b*d))/
     -   (3.*d*(-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333) + 
     -  (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -      Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -        (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -    0.3333333333333333/(3.*2**0.3333333333333333*d)

      root(2) =
     -     -c/(3.*d) + ((1 + (0,1)*Sqrt(3.))*(-c**2 + 3*b*d))/
     -   (3.*2**0.6666666666666666*d*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333) - 
     -  ((1 - (0,1)*Sqrt(3.))*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333)/(6.*2**0.3333333333333333*d)

      root(3) =
     .     -c/(3.*d) + ((1 - (0,1)*Sqrt(3.))*(-c**2 + 3*b*d))/
     -   (3.*2**0.6666666666666666*d*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333) - 
     -  ((1 + (0,1)*Sqrt(3.))*
     -     (-2*c**3 + 9*b*c*d - 27*a*d**2 + 
     -        Sqrt(4*(-c**2 + 3*b*d)**3 + 
     -          (-2*c**3 + 9*b*c*d - 27*a*d**2)**2))**
     -      0.3333333333333333)/(6.*2**0.3333333333333333*d)

c     End program

      end function solve_cubic

c     int2char
c     ################################################################
      function int2char(n) result (chr)

      implicit none

      integer(4)   :: n
      character(10):: chr

      integer(4)   :: i,exp,k,j
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

c     fmed
c     ###############################################################
      function fmed(p1,p2,p3)
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

      end function fmed

      end module math

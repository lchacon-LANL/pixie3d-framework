      module math

      use io

      INTERFACE  solve_quadratic
        module procedure solve_quadratic_real,solve_quadratic_cmplx
      END INTERFACE

      INTERFACE  solve_cubic
        module procedure solve_cubic_real,solve_cubic_cmplx
      END INTERFACE

      real(8),parameter :: pi=acos(-1d0)

      INTERFACE det
        module procedure determ3_dp,determ3_sp
      END INTERFACE

      INTERFACE fmed
        module procedure fmed_scl,fmed_vec
      END INTERFACE

      INTERFACE ellipticK
        module procedure ellipticK_scl,ellipticK_vec
      END INTERFACE

      INTERFACE ellipticE
        module procedure ellipticE_scl,ellipticE_vec
      END INTERFACE

      INTERFACE chkRoundOff
        module procedure chkRoundOff_dp,chkRoundOff_sp
      END INTERFACE

#ifdef _SPRNG
#include "sprng_f.h"
      SPRNG_POINTER stream
      logical :: is_sprng_init=.false.
!$OMP THREADPRIVATE(is_sprng_init)
#endif

      integer :: lcg_seed = 123456789
!$OMP THREADPRIVATE(lcg_seed)

      contains

c     is_NAN
c     #################################################################
      function is_NAN(a)

c     -----------------------------------------------------------------
c     Detects whether number a is NaN
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        logical :: is_NAN
        real(8) :: a

c     Local variables

c     Begin program

        is_NAN = (a /= a)

      end function is_NAN

c     is_INF
c     #################################################################
      function is_INF(a)

c     -----------------------------------------------------------------
c     Detects whether number a is INF
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        logical :: is_INF
        real(8) :: a

c     Local variables

c     Begin program

        is_INF = (abs(a) > HUGE(a))

      end function is_INF
      
c     cross_product
c     #################################################################
      function cross_product(vec1,vec2) result(vec3)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: vec1(3),vec2(3),vec3(3)

c     Local variables

c     Begin program

        vec3(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
        vec3(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
        vec3(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

      end function cross_product

c     triple_product
c     #################################################################
      function triple_product(vec1,vec2,vec3) result(scalar)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: vec1(3),vec2(3),vec3(3),scalar

c     Local variables

c     Begin program

        scalar = dot_product(vec1,cross_product(vec2,vec3))

      end function triple_product

c     delta
c     #################################################################
      function delta(i,j)

        integer    :: i,j
        real(8)    :: delta

        if (i == j) then
          delta = 1d0
        else
          delta = 0d0
        endif

      end function delta

c     findRoundOff_sp
c     ###############################################################
      function findRoundOff_sp() result(eps)

      implicit none

c     ---------------------------------------------------------------
c     Finds machine round-off constant epsmac
c     ---------------------------------------------------------------

c     Call variables

      real(4) :: eps

c     Local variables

      real(4) :: mag,mag2

      mag = 1d0
      do
        eps = mag
        mag = mag/2
        mag2 = 1d0 + mag
        if (.not.(mag2 > 1d0)) exit
      enddo

      end function findRoundOff_sp

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
      function findRoundOff() result(eps)

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

      end function findRoundOff

c     chkRoundOff_dp
c     ###############################################################
      function chkRoundOff_dp(test,var) result(eps)

      implicit none

c     ---------------------------------------------------------------
c     Finds machine round-off constant epsmac
c     ---------------------------------------------------------------

c     Call variables

      logical :: test
      real(8) :: var,eps

c     Local variables

      real(8),save :: roundoff=0d0

c     Begin program

      eps = var
      if (.not.test) return
      
      if (roundoff == 0d0) roundoff=findRoundOff()

      if (abs(var) < 1000*roundoff) eps = 0d0

      end function chkRoundOff_dp

c     chkRoundOff_sp
c     ###############################################################
      function chkRoundOff_sp(test,var) result(eps)

      implicit none

c     ---------------------------------------------------------------
c     Finds machine round-off constant epsmac
c     ---------------------------------------------------------------

c     Call variables

      logical :: test
      real(4) :: var,eps

c     Local variables

      real(4),save :: roundoff=0e0

c     Begin program

      eps = var
      if (.not.test) return
      
      if (roundoff == 0d0) roundoff=findRoundOff_sp()

      if (abs(var) < 1000*roundoff) eps = 0d0

      end function chkRoundOff_sp

c     determ3_dp
c     #################################################################
      function determ3_dp(tensor)

      real(8) :: tensor(3,3)

      real(8) :: determ3_dp

      determ3_dp = tensor(1,1)*tensor(2,2)*tensor(3,3)
     .            +tensor(3,2)*tensor(2,1)*tensor(1,3)
     .            +tensor(1,2)*tensor(2,3)*tensor(3,1)
     .            -tensor(1,3)*tensor(2,2)*tensor(3,1)
     .            -tensor(1,1)*tensor(2,3)*tensor(3,2)
     .            -tensor(3,3)*tensor(1,2)*tensor(2,1)

      end function determ3_dp

c     determ3_sp
c     #################################################################
      function determ3_sp(tensor)

      real(4) :: tensor(3,3)

      real(4) :: determ3_sp

      determ3_sp = tensor(1,1)*tensor(2,2)*tensor(3,3)
     .            +tensor(3,2)*tensor(2,1)*tensor(1,3)
     .            +tensor(1,2)*tensor(2,3)*tensor(3,1)
     .            -tensor(1,3)*tensor(2,2)*tensor(3,1)
     .            -tensor(1,1)*tensor(2,3)*tensor(3,2)
     .            -tensor(3,3)*tensor(1,2)*tensor(2,1)

      end function determ3_sp

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

!************************************************************************
!*                                                                      *
!*    Program to calculate the first kind modified Bessel function of   *
!*    integer order N, for any REAL X, using the function BESSI(N,X).   *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*                                                                      *
!*    SAMPLE RUN:                                                       *
!*                                                                      *
!*    (Calculate Bessel function for N=2, X=0.75).                      *
!*                                                                      *
!*    Bessel function of order  2 for X =  0.7500:                      *
!*                                                                      *
!*         Y =  0.73666878E-01                                          *
!*                                                                      *
!* -------------------------------------------------------------------- *
!*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!*                                                                      *
!*                               F90 Release 1.1 By J-P Moreau, Paris.  *
!*                                        (www.jpmoreau.fr)             *
!*                                                                      *
!*   Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not  *
!*                1.2067429) Aug. 2011.                                 *
!************************************************************************
c$$$PROGRAM TBESSI
c$$$
c$$$  REAL*8  BESSI, X, Y
c$$$  INTEGER N
c$$$
c$$$  N=2
c$$$  X=0.75D0
c$$$
c$$$  Y = BESSI(N,X)
c$$$
c$$$  write(*,10) N, X
c$$$  write(*,20) Y
c$$$
c$$$  stop
c$$$
c$$$ 10    format (/' Bessel Function of order ',I2,' for X=',F8.4,':')
c$$$ 20     format(/'      Y = ',E15.8/)
c$$$
c$$$END

! ----------------------------------------------------------------------
      FUNCTION BESSI(N,X)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL(8) :: X,BESSI,TOX,BIM,BI,BIP
      IF (N.EQ.0) THEN
      BESSI = BESSI0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSI = BESSI1(X)
      RETURN
      ENDIF
      IF(X.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO 12 J = M,1,-1
      BIM = BIP+DFLOAT(J)*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (J.EQ.N) BESSI = BIP
 12    CONTINUE
      BESSI = BESSI*BESSI0(X)/BI
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(X)
      REAL(8) :: X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,
     .     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,
     .     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     .     0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     .     0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------
      FUNCTION BESSI1(X)
      REAL(8) :: X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,
     .     Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     .     0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     .     -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     .     -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END FUNCTION
! ----------------------------------------------------------------------

      FUNCTION BESSK(N,X)
      IMPLICIT NONE
      INTEGER N,J
      REAL *8 X,BESSK,TOX,BK,BKM,BKP
! ------------------------------------------------------------------------
!     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
!     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
!     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
!
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE SECOND
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
! ------------------------------------------------------------------------ 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
! ------------------------------------------------------------------------
      IF (N.EQ.0) THEN
      BESSK = BESSK0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSK = BESSK1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.D0) THEN
      BESSK = 1.D30
      RETURN
      ENDIF
      TOX = 2.D0/X
      BK  = BESSK1(X)
      BKM = BESSK0(X)
      DO 11 J=1,N-1
      BKP = BKM+DFLOAT(J)*TOX*BK
      BKM = BK
      BK  = BKP
   11 CONTINUE
      BESSK = BK
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSK0(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     .     0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     .     -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF(X.EQ.0.D0) THEN
      BESSK0=1.D30
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=-LOG(X/2.D0)*BESSI0(X)
      BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/DSQRT(X)
      BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

! ----------------------------------------------------------------------
      FUNCTION BESSK1(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,
     .        -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     .      0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF(X.EQ.0.D0) THEN
      BESSK1=1.D32
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=LOG(X/2.D0)*BESSI1(X)
      BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/DSQRT(X)
      BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

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

        if (abs(res0) < 1d2*findRoundOff()) cycle

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

          if (abs(res)/abs(res0)<1d-8.or.abs(dx)/abs(x)<1d-10) exit
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

        if (abs(res0) < 1d2*findRoundOff()) cycle

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

          if (abs(res)/abs(res0)<1d-8.or.abs(dx)/abs(x)<1d-10) exit
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

c     ellipticK_vec
c     ####################################################################
      function ellipticK_vec(eta,pole) result(ellK)

      implicit none

c     Call variables

      real(8) :: eta(:),ellK(size(eta))
      logical,optional :: pole(size(eta))

c     Local variables

      real(8) :: alc0,alc1,alc2,alc3,alc4
     .          ,blc0,blc1,blc2,blc3,blc4
      logical :: ple(size(eta))

c     Begin program

      if (PRESENT(pole)) then
        ple = pole
      else
        ple = .false.
      endif

      alc0 = 1.38629436112d0
      alc1 = 0.09666344259d0
      alc2 = 0.03590092383d0
      alc3 = 0.03742563713d0
      alc4 = 0.01451196212d0
                                                 
      blc0 = 0.5d0          
      blc1 = 0.12498593597d0
      blc2 = 0.06880248576d0
      blc3 = 0.03328355346d0
      blc4 = 0.00441787012d0

      where (ple)  !Substract logarithmic singularity
        ellK = alc0+eta*(alc1+eta*(alc2+eta*(alc3+eta*alc4)))
     .       -(    +eta*(blc1+eta*(blc2+eta*(blc3+eta*blc4))))
     .       *log(eta)
      elsewhere
        ellK = alc0+eta*(alc1+eta*(alc2+eta*(alc3+eta*alc4)))
     .       -(blc0+eta*(blc1+eta*(blc2+eta*(blc3+eta*blc4))))
     .       *log(eta)
      end where

c     End program

      end function ellipticK_vec

c     ellipticE_vec
c     ####################################################################
      function ellipticE_vec(eta) result(ellE)

      implicit none

c     Call variables

      real(8) :: eta(:),ellE(size(eta))

c     Local variables

      real(8) :: alk1,alk2,alk3,alk4
     .          ,blk1,blk2,blk3,blk4

c     Begin program
                                    
      alk1 = 0.44325141463d0
      alk2 = 0.06260601220d0
      alk3 = 0.04757383546d0
      alk4 = 0.01736506451d0
                                                 
      blk1 = 0.24998368310d0
      blk2 = 0.09200180037d0
      blk3 = 0.04069697526d0
      blk4 = 0.00526449639d0

      ellE = 1d0 
     .      +eta*(alk1+eta*(alk2+eta*(alk3+eta*alk4)))
     .      -eta*(blk1+eta*(blk2+eta*(blk3+eta*blk4)))*log(eta)

c     End program

      end function ellipticE_vec

c     ellipticK_scl
c     ####################################################################
      function ellipticK_scl(eta) result(ellK)

      implicit none

c     Call variables

      real(8) :: eta,ellK

c     Local variables

      real(8) :: alc0,alc1,alc2,alc3,alc4
     .          ,blc0,blc1,blc2,blc3,blc4

c     Begin program

      alc0 = 1.38629436112d0
      alc1 = 0.09666344259d0
      alc2 = 0.03590092383d0
      alc3 = 0.03742563713d0
      alc4 = 0.01451196212d0
                                                 
      blc0 = 0.5d0          
      blc1 = 0.12498593597d0
      blc2 = 0.06880248576d0
      blc3 = 0.03328355346d0
      blc4 = 0.00441787012d0

      ellK = alc0+eta*(alc1+eta*(alc2+eta*(alc3+eta*alc4)))
     .     -(blc0+eta*(blc1+eta*(blc2+eta*(blc3+eta*blc4))))*log(eta)

c     End program

      end function ellipticK_scl

c     ellipticE_scl
c     ####################################################################
      function ellipticE_scl(eta) result(ellE)

      implicit none

c     Call variables

      real(8) :: eta,ellE

c     Local variables

      real(8) :: alk1,alk2,alk3,alk4
     .          ,blk1,blk2,blk3,blk4

c     Begin program
                                    
      alk1 = 0.44325141463d0
      alk2 = 0.06260601220d0
      alk3 = 0.04757383546d0
      alk4 = 0.01736506451d0
                                                 
      blk1 = 0.24998368310d0
      blk2 = 0.09200180037d0
      blk3 = 0.04069697526d0
      blk4 = 0.00526449639d0

      ellE = 1d0 
     .      +eta*(alk1+eta*(alk2+eta*(alk3+eta*alk4)))
     .      -eta*(blk1+eta*(blk2+eta*(blk3+eta*blk4)))*log(eta)

c     End program

      end function ellipticE_scl

c     simpson_scl
c     ####################################################################
      recursive function simpson_scl(f,a,b,epsilon,level,level_max,iout)
     .          result (simpson_result)

      implicit none

c     --------------------------------------------------------------------
c     Computes integral of f to tolerance epsilon by adaptive simpson rule
c     --------------------------------------------------------------------

c     Call variables

      interface
        function f(x)   result (f_result)
        real(8), intent (in) :: x
        real(8) :: f_result
        end function f
      end interface

      real(8), intent(in) :: a, b, epsilon
      real(8) :: simpson_result
      integer :: level,level_max,iout

c     Local variables

      real(8) :: h, c, d, e
      real(8) :: one_simpson, two_simpson
      real(8) :: left_simpson, right_simpson

c     Begin program

      level = level + 1

      h = b - a
      c = (a+b)/2.0

c     Perform two-level integral

      if (iout > 0) then
        print *,"a, b", a, b      
        print *, "level", level
      endif

      one_simpson = h*(f(a) + 4d0*f(c) + f(b))/6d0
      d = 0.5d0*(a+c)
      e = 0.5d0*(c+b)		
      two_simpson = h*(f(a) + 4d0*f(d) + 2d0*f(c)
     .                      + 4d0*f(e) + f(b))/12d0

c     Check convergence

      if (level >= level_max) then
         simpson_result = two_simpson
      else   
         if ( abs(two_simpson - one_simpson) <= 15.0*epsilon ) then
            simpson_result = two_simpson
            if (iout > 0) then
              print *, "simpson NO SPLIT; a, b=", a, b,'level=',level
     .               ,"Result=", simpson_result
            endif
         else
            left_simpson
     .          = simpson_scl(f,a,c,0.5*epsilon,level,level_max,iout)
            right_simpson
     .          = simpson_scl(f,c,b,0.5*epsilon,level,level_max,iout)
            simpson_result = left_simpson + right_simpson
            if (iout > 0) then
               print *, "Simpson SPLIT; a, b", a, b,'level=',level
     .               ,"Result=", simpson_result
            endif
         end if
      end if   

      end function simpson_scl

ccc     simpson_vec
ccc     ####################################################################
cc      recursive function simpson_vec(n,f,a,b,eps,level,level_max,iout)
cc     .          result (simpson_result)
cc
cc      implicit none
cc
ccc     --------------------------------------------------------------------
ccc     Computes integral of f to tolerance eps by adaptive simpson rule
ccc     --------------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer :: level,level_max,iout,n
cc
cc      interface
cc        function f(n,x) result (f_result)
cc          integer :: n
cc          real(8), intent (in) :: x
cc          real(8) :: f_result(n)
cc        end function f
cc      end interface
cc
cc      real(8), intent(in) :: a, b, eps
cc      real(8) :: simpson_result(n)
cc
ccc     Local variables
cc
cc      real(8) :: h, c, d, e
cc      real(8) :: one_simpson(n), two_simpson(n)
cc      real(8) :: left_simpson(n), right_simpson(n)
cc
ccc     Begin program
cc
cc      level = level + 1
cc
cc      h = b - a
cc      c = (a+b)/2.0
cc
ccc     Perform two-level integral
cc
cc      if (iout > 0) then
cc        print *,"a, b", a, b      
cc        print *, "level", level
cc      endif
cc
cc      one_simpson = h*(f(n,a) + 4d0*f(n,c) + f(n,b))/6d0
cc      d = 0.5d0*(a+c)
cc      e = 0.5d0*(c+b)		
cc      two_simpson = h*(f(n,a) + 4d0*f(n,d) + 2d0*f(n,c)
cc     .                      + 4d0*f(n,e) + f(n,b))/12d0
cc
ccc     Check convergence
cc
cc      if (level >= level_max) then
cc         simpson_result = two_simpson
cc      else   
cc         if ( abs(two_simpson - one_simpson) <= 15.0*eps ) then
cc            simpson_result = two_simpson
cc            if (iout > 0) then
cc              print *, "simpson NO SPLIT; a, b=", a, b,'level=',level
cc     .               ,"Result=", simpson_result
cc            endif
cc         else
cc            left_simpson  = simpson(f,a,c,0.5*eps,level,level_max)
cc            right_simpson = simpson(f,c,b,0.5*eps,level,level_max)
cc            simpson_result = left_simpson + right_simpson
cc            if (iout > 0) then
cc               print *, "Simpson SPLIT; a, b", a, b,'level=',level
cc     .               ,"Result=", simpson_result
cc            endif
cc         end if
cc      end if   
cc
cc      end function simpson_vec

c     gauss2 (Box-Muller transformation)
c     #################################################################
      SUBROUTINE gauss2(y1)

      implicit none

c     arguments      

      real(8), intent(out) :: y1
      
c     local variables

      real(8) :: x1,x2,w,r

c     Begin program

      call unif_rng(r)
      x1 = 2d0*r - 1d0      

      call unif_rng(r)
      x2 = 2d0*r - 1d0
      w = x1*x1 + x2*x2

      do while(w>1d0)
        call unif_rng(r)
        x1 = 2d0*r - 1d0

        call unif_rng(r)
        x2 = 2d0*r - 1d0
        w = x1*x1 + x2*x2
      end do
      
      w = sqrt((-2d0*log(w) ) / w )
      
      y1 = x2*w

c     End program
      
      end subroutine gauss2

c     unif_rng
c     #################################################################
      SUBROUTINE unif_rng(r)

      implicit none

c     arguments      

      real(8), intent(out) :: r
      
c     local variables

c     Begin program

#ifdef _SPRNG
      r = sprng(stream)         
#else
#ifdef _OPENMP
      r = randlcg(lcg_seed)         !careful
#else
      call random_number(r)         !Problematic with OpenMP
#endif
#endif

c     End program
      
      end subroutine unif_rng

c     randlcg
c     ##################################################################
      function randlcg(seed)

      implicit none

      integer :: seed
      real(8) :: randlcg

      seed = int(mod(16807d0*dble(seed),2147483647d0))
      randlcg = dble(seed)/2147483647d0

      end function randlcg

c     init_random_seed
c     #########################################################################
      SUBROUTINE init_random_seed(nstreams,streamnum,output,seed)

c     -------------------------------------------------------------------------
c     Initializes random seed for stream number "streamnum" of "nstreams"
c     -------------------------------------------------------------------------
#ifdef _OPENMP
        use omp_lib
#endif

        implicit none
        integer,INTENT(IN) :: streamnum,nstreams
        logical,INTENT(IN) :: output
        integer,INTENT(IN),optional :: seed

#ifdef _SPRNG
        integer :: junk,sd
        integer :: gtype
        real(8) :: rn

        if (is_sprng_init) return
        is_sprng_init = .true.
        
        gtype = 1       !lcg
        if (PRESENT(seed)) then
          sd = 985456376 + 37*seed
        else
          sd = 985456376
        endif
        
        stream = init_sprng(gtype,streamnum,nstreams,sd,SPRNG_DEFAULT)   

        if (output) junk = print_sprng(stream)
#else
#ifdef _OPENMP
        if (PRESENT(seed)) then
          lcg_seed = 985456376 + seed + streamnum*37
        else
          lcg_seed = 985456376        + streamnum*37
        endif
#else
        INTEGER :: i, n, clock, sd
        INTEGER, DIMENSION(:), ALLOCATABLE :: sd2

        CALL RANDOM_SEED(size = n)
        ALLOCATE(sd2(n))

        if (PRESENT(seed)) then
          sd = seed
        else
          sd = streamnum
        endif
          
c$$$        CALL SYSTEM_CLOCK(COUNT=clock)
c$$$          
c$$$        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        sd2 = sd + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = sd2)
          
        DEALLOCATE(sd2)
#endif
#endif
      end subroutine init_random_seed
  
c     free_random_seed
c     #########################################################################
      SUBROUTINE free_random_seed()

        implicit none

#ifdef _SPRNG
        integer :: junk

        junk = free_sprng(stream)
#endif

      END SUBROUTINE free_random_seed

c     erfcc
c     #########################################################################
      FUNCTION erfcc(x)

      implicit none

c     ------------------------------------------------------------------------
c     Returns the complementary error function erfc(x)
c     with fractional error everywhere less than 1.2e-7
c     ------------------------------------------------------------------------

c     Call variables

      REAL(8) :: erfcc,x

c     Local variables

      REAL(8) :: t,z

c     Begin program

      z=abs(x)

      t=1d0/(1.+0.5*z)

      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     $     t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     $     t*(1.48851587+t*(-.82215223+t*.17087277)))))))))

      if (x.lt.0d0) erfcc=2d0-erfcc

      END FUNCTION erfcc

c     CALERF
c     #######################################################################
      SUBROUTINE CALERF(ARG,RESULT,JINT)
C------------------------------------------------------------------
C
C This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
C   for a real argument  x.  It contains three FUNCTION type
C   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
C   and one SUBROUTINE type subprogram, CALERF.  The calling
C   statements for the primary entries are:
C
C                   Y=ERF(X)     (or   Y=DERF(X)),
C
C                   Y=ERFC(X)    (or   Y=DERFC(X)),
C   and
C                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
C
C   The routine  CALERF  is intended for internal packet use only,
C   all computations within the packet being concentrated in this
C   routine.  The function subprograms invoke  CALERF  with the
C   statement
C
C          CALL CALERF(ARG,RESULT,JINT)
C
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALERF
C       call              ARG                  Result          JINT
C
C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
C     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
C     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
C
C   The main computation evaluates near-minimax approximations
C   from "Rational Chebyshev approximations for the error function"
C   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
C   transportable program uses rational functions that theoretically
C   approximate  erf(x)  and  erfc(x)  to at least 18 significant
C   decimal digits.  The accuracy achieved depends on the arithmetic
C   system, the compiler, the intrinsic functions, and proper
C   selection of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XMIN   = the smallest positive floating-point number.
C   XINF   = the largest positive finite floating-point number.
C   XNEG   = the largest negative argument acceptable to ERFCX;
C            the negative of the solution to the equation
C            2*exp(x*x) = XINF.
C   XSMALL = argument below which erf(x) may be represented by
C            2*x/sqrt(pi)  and above which  x*x  will not underflow.
C            A conservative value is the largest machine number X
C            such that   1.0 + X = 1.0   to machine precision.
C   XBIG   = largest argument acceptable to ERFC;  solution to
C            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
C            W(x) = exp(-x*x)/[x*sqrt(pi)].
C   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
C            machine precision.  A conservative value is
C            1/[2*sqrt(XSMALL)]
C   XMAX   = largest acceptable argument to ERFCX; the minimum
C            of XINF and 1/[sqrt(pi)*XMIN].
C
C   Approximate values for some important machines are:
C
C                          XMIN       XINF        XNEG     XSMALL
C
C  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
C  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
C  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
C  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
C  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
C  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
C
C
C                          XBIG       XHUGE       XMAX
C
C  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
C  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
C  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
C  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
C  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
C  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ERFC = 0      for  ARG .GE. XBIG;
C
C                       ERFCX = XINF  for  ARG .LT. XNEG;
C      and
C                       ERFCX = 0     for  ARG .GE. XMAX.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 19, 1990
C
C------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1     A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI,
     2     TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL,
     3     Y,YSQ,ZERO,IY
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
C------------------------------------------------------------------
C  Mathematical constants
C------------------------------------------------------------------
CS    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
CS   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
CS   2     SIXTEN/16.0E0/
      DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/,
     1     SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/,
     2     SIXTEEN/16.0D0/,ISIXTEEN/0.0625D0/
C------------------------------------------------------------------
C  Machine-dependent constants
C------------------------------------------------------------------
CS    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
CS   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
      DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/,
     1     XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
C------------------------------------------------------------------
C  Coefficients for approximation to  erf  in first interval
C------------------------------------------------------------------
CS    DATA A/3.16112374387056560E00,1.13864154151050156E02,
CS   1       3.77485237685302021E02,3.20937758913846947E03,
CS   2       1.85777706184603153E-1/
CS    DATA B/2.36012909523441209E01,2.44024637934444173E02,
CS   1       1.28261652607737228E03,2.84423683343917062E03/
      DATA A/3.16112374387056560D00,1.13864154151050156D02,
     1       3.77485237685302021D02,3.20937758913846947D03,
     2       1.85777706184603153D-1/
      DATA B/2.36012909523441209D01,2.44024637934444173D02,
     1       1.28261652607737228D03,2.84423683343917062D03/
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in second interval
C------------------------------------------------------------------
CS    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
CS   1       6.61191906371416295E01,2.98635138197400131E02,
CS   2       8.81952221241769090E02,1.71204761263407058E03,
CS   3       2.05107837782607147E03,1.23033935479799725E03,
CS   4       2.15311535474403846E-8/
CS    DATA D/1.57449261107098347E01,1.17693950891312499E02,
CS   1       5.37181101862009858E02,1.62138957456669019E03,
CS   2       3.29079923573345963E03,4.36261909014324716E03,
CS   3       3.43936767414372164E03,1.23033935480374942E03/
      DATA C/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/
      DATA D/1.57449261107098347D01,1.17693950891312499D02,
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in third interval
C------------------------------------------------------------------
CS    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
CS   1       1.25781726111229246E-1,1.60837851487422766E-2,
CS   2       6.58749161529837803E-4,1.63153871373020978E-2/
CS    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
CS   1       5.27905102951428412E-1,6.05183413124413191E-2,
CS   2       2.33520497626869185E-3/
      DATA P/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/
      DATA Q/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/
C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      if (Y > 0d0) IY = ONE/Y

      IF (Y .LE. THRESH) THEN
C------------------------------------------------------------------
C  Evaluate  erf  for  |X| <= 0.46875
C------------------------------------------------------------------
            YSQ = ZERO
            IF (Y .GT. XSMALL) YSQ = Y * Y
            XNUM = A(5)*YSQ
            XDEN = YSQ
            DO 20 I = 1, 3
               XNUM = (XNUM + A(I)) * YSQ
               XDEN = (XDEN + B(I)) * YSQ
   20       CONTINUE
            RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
            IF (JINT .NE. 0) RESULT = ONE - RESULT
            IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
            GO TO 800
C------------------------------------------------------------------
C  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
C------------------------------------------------------------------
      ELSE IF (Y .LE. FOUR) THEN
            XNUM = C(9)*Y
            XDEN = Y
            DO 120 I = 1, 7
               XNUM = (XNUM + C(I)) * Y
               XDEN = (XDEN + D(I)) * Y
  120       CONTINUE
            RESULT = (XNUM + C(8)) / (XDEN + D(8))
            IF (JINT .NE. 2) THEN
               YSQ = AINT(Y*SIXTEEN)*ISIXTEEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ-DEL)*RESULT
            END IF
C------------------------------------------------------------------
C  Evaluate  erfc  for |X| > 4.0
C------------------------------------------------------------------
      ELSE
            RESULT = ZERO
            IF (Y .GE. XBIG) THEN
               IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
               IF (Y .GE. XHUGE) THEN
                  RESULT = SQRPI * IY !/ Y
                  GO TO 300
               END IF
            END IF
            YSQ = IY * IY !ONE / (Y * Y)
            XNUM = P(6)*YSQ
            XDEN = YSQ
            DO 240 I = 1, 4
               XNUM = (XNUM + P(I)) * YSQ
               XDEN = (XDEN + Q(I)) * YSQ
  240       CONTINUE
            RESULT = (SQRPI -  YSQ *(XNUM + P(5)) / (XDEN + Q(5)))*IY ! / Y
            IF (JINT .NE. 2) THEN
               YSQ = AINT(Y*SIXTEEN)*ISIXTEEN
               DEL = (Y-YSQ)*(Y+YSQ)
               RESULT = EXP(-YSQ*YSQ -DEL) * RESULT
            END IF
      END IF
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
  300 IF (JINT .EQ. 0) THEN
            RESULT = (HALF - RESULT) + HALF
            IF (X .LT. ZERO) RESULT = -RESULT
         ELSE IF (JINT .EQ. 1) THEN
            IF (X .LT. ZERO) RESULT = TWO - RESULT
         ELSE
            IF (X .LT. ZERO) THEN
               IF (X .LT. XNEG) THEN
                     RESULT = XINF
                  ELSE
                     YSQ = AINT(X*SIXTEEN)*ISIXTEEN
                     DEL = (X-YSQ)*(X+YSQ)
                     Y = EXP(YSQ*YSQ + DEL)
                     RESULT = (Y+Y) - RESULT
               END IF
            END IF
      END IF
  800 RETURN
C---------- Last card of CALERF ----------
      END
CS    REAL FUNCTION ERF(X)
      DOUBLE PRECISION FUNCTION DERF(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for erf(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 0
      CALL CALERF(X,RESULT,JINT)
CS    ERF = RESULT
      DERF = RESULT
      RETURN
C---------- Last card of DERF ----------
      END
CS    REAL FUNCTION ERFC(X)
      DOUBLE PRECISION FUNCTION DERFC(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 1
      CALL CALERF(X,RESULT,JINT)
CS    ERFC = RESULT
      DERFC = RESULT
      RETURN
C---------- Last card of DERFC ----------
      END
CS    REAL FUNCTION ERFCX(X)
      DOUBLE PRECISION FUNCTION DERFCX(X)
C------------------------------------------------------------------
C
C This subprogram computes approximate values for exp(x*x) * erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, March 30, 1987
C
C------------------------------------------------------------------
      INTEGER JINT
CS    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
C------------------------------------------------------------------
      JINT = 2
      CALL CALERF(X,RESULT,JINT)
CS    ERFCX = RESULT
      DERFCX = RESULT
      RETURN
C---------- Last card of DERFCX ----------
      END
      
!**********************************************************************************
!  Copyright (C) 1985 Numerical Recipes Software -- GAMMLN
!**********************************************************************************
      FUNCTION GAMMLN(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX
!c      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
!c    *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA COF,STP/76.18009173,-86.50532033,24.01409822,
     *    -1.231739516,.120858003E-2,-.536382E-5,2.50662827465/
      DATA HALF,ONE,FPF/0.5,1.0,5.5/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END


!c     SYMINV
!c     ########################################################################
      SUBROUTINE SYMINV(A,N,C,W,NULLTY,IFAULT,RMAX,CMY,TOL)
!C
!C       Algorithm AS7, Applied Statistics, Vol.17, 1968.
!C
!C       ARGUMENTS:-
!C       A()     = input, the symmetric matrix to be inverted, stored in
!C                 lower triangular form
!C       N       = input, order of the matrix
!C       C()     = output, the inverse of a (a generalized inverse if c is
!C                 singular), also stored in lower triangular.
!C                 c and a may occupy the same locations.
!C       W()     = workspace, dimension at least n.
!C       NULLTY  = output, the rank deficiency of a.
!C       IFAULT  = output, error indicator
!C                       = 1 if n < 1
!C                       = 2 if a is not +ve semi-definite
!C                       = 0 otherwise
!C       RMAX    = output, approximate bound on the accuracy of the diagonal
!C                 elements of c.  e.g. if rmax = 1.e-04 then the diagonal
!C                 elements of c will be accurate to about 4 dec. digits.
!C
!C       Latest revision - 18 April 1981

!C        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT NONE             ! GC 2016      
!C     call variables
       integer :: N,NULLTY,IFAULT
       real(8) :: A(N*N),C(N*N),W(N),CMY(N*(N+1)/2),RMAX,TOL

!C     local varialbe
       integer :: NROW,NN,IMY,IROW,ICOL,JCOL,MDIAG,L,NDIAG,K,I,J
       real(8) :: X

        NROW=N
        IFAULT=1
        IF(NROW.LE.0) GO TO 100
        IFAULT=0
!C       Cholesky factorization of A, result in C
        CALL CHOLA(A,NROW,C,NULLTY,IFAULT,RMAX,W,TOL)
        IF(IFAULT.NE.0) GO TO 100
!C       Invert C & form the product (CINV)'*CINV, where CINV is the inverse
!C       of C, row by row starting with the last row.
!C       IROW = the row number, NDIAG = location of last element in the row.
        NN=NROW*(NROW+1)/2
        DO 200 IMY=1,NN
200       CMY(IMY)=C(IMY)
        IROW=NROW
        NDIAG=NN
10      IF(C(NDIAG).EQ.0.0) GO TO 60
        L=NDIAG
        DO 20 I=IROW,NROW
          W(I)=C(L)
          L=L+I
20      CONTINUE
        ICOL=NROW
        JCOL=NN
        MDIAG=NN
30      L=JCOL
        X=0.0
        IF(ICOL.EQ.IROW) X=1.0/W(IROW)
        K=NROW
40      IF(K.EQ.IROW) GO TO 50
        X=X-W(K)*C(L)
        K=K-1
        L=L-1
        IF(L.GT.MDIAG) L=L-K+1
        GO TO 40
50      C(L)=X/W(IROW)
        IF(ICOL.EQ.IROW) GO TO 80
        MDIAG=MDIAG-ICOL
        ICOL=ICOL-1
        JCOL=JCOL-1
        GO TO 30
60      L=NDIAG
        DO 70 J=IROW,NROW
          C(L)=0.0
          L=L+J
70      CONTINUE
80      NDIAG=NDIAG-IROW
        IROW=IROW-1
        IF(IROW.NE.0) GO TO 10
100     RETURN
        END SUBROUTINE SYMINV

!c     CHOLA
!c     ########################################################################
        SUBROUTINE CHOLA(A,N,U,NULLTY,IFAULT,RMAX,R,TOL)
!C
!C       Algorithm AS6, Applied Statistics, Vol.17, 1968, with
!C       modifications by A.J.Miller
!C
!C       ARGUMENTS:-
!C       A()     = input, a +ve definite matrix stored in lower-triangular
!C                 form.
!C       N       = input, the order of A
!C       U()     = output, a lower triangular matrix such that U*U' = A.
!C                 A & U may occupy the same locations.
!C       NULLTY  = output, the rank deficiency of A.
!C       IFAULT  = output, error indicator
!C                       = 1 if N < 1
!C                       = 2 if A is not +ve semi-definite
!C                       = 0 otherwise
!C       RMAX    = output, an estimate of the relative accuracy of the
!C                 diagonal elements of U.
!C       R()     = output, array containing bounds on the relative accuracy
!C                 of each diagonal element of U.
!C
!C       Latest revision - 18 April 1981
!C
!C        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT NONE ! GC 2016
!C     call variables
       integer :: N,NULLTY,IFAULT
       real(8) :: A(N*N),U(N*N),R(N),RMAX,TOL
       
!C     Local variables
       integer :: I,J,K,ICOL,IROW,L,M
       real(8) :: ETA,W,RSQ

!C       ETA should be set equal to the smallest +ve value such that
!C       1.0 + eta is calculated as being greater than 1.0 in the accuracy
!C       being used.
!C        DATA ETA/1.e-07/
        ETA=TOL
        IFAULT=1
        IF(N.LE.0) GO TO 100
        IFAULT=2
        NULLTY=0
        RMAX=ETA
        R(1)=ETA
        J=1
        K=0
!C       Factorize column by column, ICOL = Column No.
        DO 80 ICOL=1,N
          L=0
!C       IROW = row number within column ICOL
          DO 40 IROW=1,ICOL
            K=K+1
            W=A(K)
            IF(IROW.EQ.ICOL) RSQ=(W*ETA)**2
            M=J
            DO 10 I=1,IROW
              L=L+1
              IF(I.EQ.IROW) GO TO 20
              W=W-U(L)*U(M)
              IF(IROW.EQ.ICOL) RSQ=RSQ+(U(L)**2*R(I))**2
              M=M+1
10          CONTINUE
20          IF(IROW.EQ.ICOL) GO TO 50
            IF(U(L).EQ.0.0) GO TO 30
            U(K)=W/U(L)
            GO TO 40
30          U(K)=0.0
            IF(ABS(W).GT.ABS(RMAX*A(K))) GO TO 100
40        CONTINUE
!C         End of row, estimate relative accuracy of diagonal element.
50        RSQ=SQRT(RSQ)
          IF(ABS(W).LE.5.*RSQ) GO TO 60
          IF(W.LT.0.0) GO TO 100
          U(K)=SQRT(W)
          R(I)=RSQ/W
          IF(R(I).GT.RMAX) RMAX=R(I)
          GO TO 70
60        U(K)=0.0
          NULLTY=NULLTY+1
70        J=J+ICOL
80      CONTINUE
        IFAULT=0.0
100     RETURN
        END SUBROUTINE CHOLA

c     blockSolve
c     #################################################################
      subroutine blockSolve(size,mat,icol,rhs,x)

c     -----------------------------------------------------------------
c     Solves block systems using a direct solve approach. Requires
c     linking with LAPACK.
c
c     In the call sequence, we have:
c       * size (integer): block size
c       * mat  (real array): block matrix
c       * icol (integer) : number of columns of rhs
c       * rhs  (real array): rhs of equation
c       * x    ( "     "   ): solution (output)
c       * smoother (character): identifies calling subroutine (JB, GS)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: size,icol

        real(8)    :: mat(size,size),rhs(size,icol),x(size,icol)

c     Local variables

        integer    :: ipiv(size),info

        real(8)    :: mat2(size,size)

        external dgesv

c     Begin program

        mat2 = mat !Avoid overwritting mat
        x = rhs

        call dgesv(size,icol,mat2,size,ipiv,x,size,info) !LAPACK routine

        call LAPACK_error(info)

      end subroutine blockSolve

c     blockInv
c     #################################################################
      subroutine blockInv(size,mat)

c     -----------------------------------------------------------------
c     Inverts block systems using a direct solve approach. Requires
c     linking with LAPACK.
c
c     In the call sequence, we have:
c       * size (integer): block size
c       * mat  (real array): block matrix on input, inverse on output
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: size

        real(8)    :: mat(size,size)

c     Local variables

        integer    :: ipiv(size),info,lwork

        real(8),allocatable,dimension(:) :: work

        external   dgetri,dgetrf

c     Begin program

c     Find LU decomposition

        call dgetrf(size,size,mat,size,ipiv,info) !LAPACK routine

        call LAPACK_error(info)

c     Invert matrix

        !Workspace query
        lwork=-1
cc        lwork=size
        allocate(work(1))
        call dgetri(size,mat,size,ipiv,work,lwork,info) !LAPACK routine
        lwork=floor(work(1))
        deallocate(work)

        !Matrix inversion
        allocate(work(lwork))
        call dgetri(size,mat,size,ipiv,work,lwork,info) !LAPACK routine
        deallocate(work)

        call LAPACK_error(info)

      end subroutine blockInv

c     LAPACK_error
c     #################################################################
      subroutine LAPACK_error(info)

c     -----------------------------------------------------------------
c     Error routine for LAPACK calls
c     -----------------------------------------------------------------

        use io
      
        implicit none

c     Call variables

        integer    :: info

c     Local variables

        character(50) :: messg

c     Begin program

        if (info /= 0) then
          if (info < 0) then
            messg = 'Problem in factorization in argument '
     .             //int2char(-info)
            call sstop(0,'blockSolve',messg)
          else
            call sstop(0,'blockSolve','Matrix is singular')
          endif
        endif

c     End program

      end subroutine LAPACK_error

      end module math

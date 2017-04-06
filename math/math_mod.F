      module math

      INTERFACE  solve_quadratic
        module procedure solve_quadratic_real,solve_quadratic_cmplx
      END INTERFACE

      INTERFACE  solve_cubic
        module procedure solve_cubic_real,solve_cubic_cmplx
      END INTERFACE

      real(8),parameter :: pi=acos(-1d0)

      INTERFACE det
        module procedure determ3
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

#ifdef _SPRNG
#include "sprng_f.h"

        integer ::junk
        SPRNG_POINTER stream
#endif

      contains

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

c     chkRoundOff
c     ###############################################################
      function chkRoundOff(test,var) result(eps)

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

      end function chkRoundOff

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
      integer, save :: seed = 1 !careful

c     Begin program

#ifdef _SPRNG
      r = sprng(stream)         
#else
      CALL RANDOM_NUMBER(r)
cc      r = randlcg(seed)         !careful
#endif
      x1 = 2d0*r - 1d0      
#ifdef _SPRNG
      r = sprng(stream)         
#else
      CALL RANDOM_NUMBER(r)
cc      r = randlcg(seed)         !careful
#endif
      x2 = 2d0*r - 1d0
      w = x1*x1 + x2*x2

      do while(w>1d0)
#ifdef _SPRNG
        r = sprng(stream)         
#else
        CALL RANDOM_NUMBER(r)
cc        r = randlcg(seed)       !careful
#endif
        x1 = 2d0*r - 1d0
#ifdef _SPRNG
        r = sprng(stream)         
#else
        CALL RANDOM_NUMBER(r) 
cc        r = randlcg(seed)
#endif
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
cc      r = randlcg(seed)         !careful
      call random_number(r)
#endif

c     End program
      
      end subroutine unif_rng

c     randlcg
c     ##################################################################
      function randlcg(seed)

      implicit none

      integer :: seed
      real(8) :: randlcg

      seed = int(mod(16807d0*dble(seed), 2147483647d0))
      randlcg = dble(seed)/2147483647d0      

      end function randlcg

c     init_random_seed
c     #########################################################################
      SUBROUTINE init_random_seed(nproc,streamnum,output)

        implicit none
        integer,INTENT(IN) :: streamnum,nproc
        logical,INTENT(IN) :: output

#ifdef _SPRNG
        integer :: nstreams, junk,seed
        integer :: gtype
        real(8) :: rn

        gtype = 1       !lcg
        seed = 985456376      
        nstreams = nproc
        
        stream = init_sprng(gtype,streamnum,nstreams,seed,SPRNG_DEFAULT)   

        if (output) junk = print_sprng(stream)
#else
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
          
c$$$        CALL SYSTEM_CLOCK(COUNT=clock)
c$$$          
c$$$        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        seed = streamnum + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
          
        DEALLOCATE(seed)
#endif
        end subroutine init_random_seed
  
c     free_random_seed
c     #########################################################################
      SUBROUTINE free_random_seed()

        implicit none

#ifdef _SPRNG
        integer ::junk

        junk = free_sprng(stream)
#endif

      END SUBROUTINE free_random_seed

      end module math

      module math

      contains

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
      
      end module math

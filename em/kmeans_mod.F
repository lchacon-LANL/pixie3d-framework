C     module KMEANS_MOD
C     #######################################################################
      MODULE KMEANS_MOD
C     This group of subroutines implements the K-means Clustering algorithm.
C     Implemented by David Peel May 1994

      use gm_utils

      integer,private :: mpierr

      CONTAINS

!######################################################################
!   GMM_KMEANS_INIT
!######################################################################
      
!c     GMM_INIT_KMEANS
!c     ##################################################################
      SUBROUTINE GMM_INIT_KMEANS(NIND,NATT,NG,X,XMU,XVAR,T,WTOT,IDT)
!C     This subroutine uses the KMEANS CLUSTERING to initialize
!C     Gaussian parameters for EM. Also adaptively 'kill' cluster(s)
!C     that have small number of points

!c     This will take initial values of XMU from INIT_GMMEM
!c     subroutine as initial values for centroids in K-MEANS
!c     clustering algorithm to provide 'better' initial values
!c     for XMU. On outputs, number of elements & XMU's initial
!c     values are given for clusters 1->NG. Using these info,
!c     init values for proportions T & for XVAR can be computed.
!c     See GMM_KMEANS_INIT subroutine

      IMPLICIT NONE

!c     Call variables
      
      INTEGER, INTENT(IN) :: NIND,NATT
      INTEGER, INTENT(INOUT) :: NG
      INTEGER, INTENT(INOUT) :: IDT(NIND)
      REAL(8), INTENT(IN) :: X(NIND,NATT), WTOT
      REAL(8), INTENT(INOUT):: T(NG)
     .                      ,XMU(NG,NATT)
     .                      ,XVAR(NG,NATT,NATT)

!c     Local variables
      
      REAL(8) :: RAN            !Random number (0,1)
     $          ,TEMP
     $          ,KM_TOL,KLdist
     $          ,XVAR_K(1:NATT,1:NATT)
     $		,LSUMX(1:NATT,1:NATT)
      INTEGER :: I,J,K,II,JJ
     $          ,RANI,MAXIT,ITOUT,IER
      INTEGER, PARAMETER :: DP = KIND(1.0D0)    
      INTEGER :: GRP_COUNT(NG)

!c     Begin program

      MAXIT = 200 	! set max iter = 200 for K-MEANS

      if (my_rank_gmm==0.and.std_io) then
!$OMP CRITICAL
        write (*,*) "Performing KMEANS GMM initialization"
!$OMP END CRITICAL
      endif
      
      IF (DIAG_IO) THEN
        WRITE(FYLENO,*) '--------------- K-MEANS OUTPUT ---------------'
        WRITE(FYLENO,977) NG
 977    FORMAT (/2X,'Before K-MEANS, NG = ',2I2)
      END IF
      
      GRP_COUNT = 0
      if (NIND > 0) then
        KM_TOL = 1D-4
!ccc	Using K-means clustering to obtain 'good guess' for means
!ccc	call KMEANS(NIND,NATT,NG,X,XMU,IDT,KM_TOL,MAXIT,ITOUT,IER)
!c	Using Adaptive K-Means to obtain 'good guesses' for number
!c     of clusters to be used for EM and those clusters' respective means

!c      TODO: add variable particle weight
        call ADAPT_KMEANS(NIND,NATT,NG,X(1:NIND,1:NATT),XMU(1:NG,1:NATT)
     .                   ,IDT(1:NIND),GRP_COUNT(1:NG),KM_TOL,MAXIT,ITOUT
     .                   ,IER)

!c	NG changed after the above call!
!c	Calculate groups' proportions
        DO K = 1,NG
          T(K) = real(GRP_COUNT(K)/WTOT,DP) ! proportion of group K
!C      Compute new estimate of covariance matrix for each group
          DO J=1,NATT
           DO I=1,J
             XVAR_K(I,J)=0.0D0
             DO JJ=1,NIND
              IF ( IDT(JJ) .EQ. K ) THEN
                XVAR_K(I,J)=XVAR_K(I,J)
     .             +(X(JJ,I)-XMU(K,I))*(X(JJ,J)-XMU(K,J))
              ENDIF
            END DO
            XVAR_K(J,I)=XVAR_K(I,J)
          END DO
         END DO
#if defined(petsc)
         lsumx = XVAR_K
         call MPI_Allreduce(lsumx,XVAR_K,NATT*NATT
     .              ,MPI_DOUBLE_PRECISION,MPI_SUM
     .              ,MPI_COMM_WORLD,mpierr)
#endif

         XVAR_K = XVAR_K/dble(GRP_COUNT(K)) + (1.11d-16)
         XVAR(K,1:NATT,1:NATT) = XVAR_K(1:NATT,1:NATT)
        END DO  ! end do K = 1,NG
      endif	! end if (NIND > 0)

!c	Ttot = sum(T(1:NG))
!c	XMUtot = 0D0
!c	DO K=1,NG
!c	  XMUtot(:) = XMUtot(:) + T(K)*XMU(K,:)
!c	END DO
	
!c	WRITE(*,*) 'WTOT after K-Means:',Ttot
!c	WRITE(*,*) 'XMUTOT after K-Means:',(XMUtot(I),I=1,NATT)

      IF (DIAG_IO) THEN
        WRITE(FYLENO,1077) NG
 1077   FORMAT (/2X,'After K-MEANS, NG = ',2I2)
        WRITE(FYLENO,*) 'K-MEANS grouping of data'
        WRITE (FYLENO,1177) (IDT(I),I=1,NIND)
 1177   FORMAT (2X,10I4)
        WRITE(FYLENO,1277) (GRP_COUNT(K),K=1,NG)
 1277   FORMAT (/2X,'Group Counts:',20I4)
        WRITE (FYLENO,*) 'MAX K-MEANS ITERS = ',ITOUT
        WRITE(FYLENO,1377) NG        
 1377   FORMAT (/2X,'Group-means for ',I2,' groups')
        DO K=1,NG
          WRITE (FYLENO,1477) (XMU(K,J),J=1,NATT)
 1477     FORMAT (2X,5G13.5)
        END DO
        WRITE(FYLENO,*) '------------- END K-MEANS OUTPUT -------------'
      ENDIF
      
      END SUBROUTINE GMM_INIT_KMEANS

      SUBROUTINE ADAPT_KMEANS(NIND,NATT,NG,X,XMU,IDT,GRP_CT
     $			     ,EPSL,TT,T,IER)
C      Main subroutine
      
       IMPLICIT NONE

C      Call variables

       INTEGER, INTENT(IN) :: NIND,NATT,TT
       INTEGER, INTENT(INOUT) :: NG
       INTEGER, INTENT(OUT) :: IER,T
       REAL(8), INTENT(IN) :: X(NIND,NATT),EPSL
       REAL(8), INTENT(INOUT) :: XMU(NG,NATT)
       INTEGER, INTENT(INOUT) :: GRP_CT(NG),IDT(NIND)

C      Local variables

       INTEGER I,J,II,KK,LL,LGRP(NG)
       REAL(8) :: XK(NG,NATT),XKOLD(NG,NATT)
     &           ,XSTAN(NIND,NATT),ET
       

      IER=0
c      XMU(1:NG,1:NATT) = 0D0
c      CALL KSTAND(NIND,NATT,X,XSTAN) 		! normalize sample
c      CALL KSEED(NIND,NATT,NG,XSTAN,XK,IER)	! pick centroids randomly from normalized sample.
      XSTAN(1:NIND,1:NATT) = X(1:NIND,1:NATT)	! do not normalize sample X
c      CALL KSEED(NIND,NATT,NG,XSTAN,XK,IER)	! pick centroids randomly from original sample.
      XK(1:NG,1:NATT) = XMU(1:NG,1:NATT)	! pick centroids = pseudo-random init of means
      DO T=1,TT				! iter loop
        GRP_CT = 0	
        XKOLD = XK

	DO KK=1,NIND		! particle loop & assign particle to a cluster
          CALL WINNER1(NIND,NATT,NG,XSTAN(KK,1:NATT)
     $		      ,XK(1:NG,1:NATT),IDT(KK),GRP_CT(1:NG),IER)
        END DO

#if defined(petsc)
        LGRP(1:NG) = GRP_CT(1:NG)
        call MPI_Allreduce(lgrp(1:NG),GRP_CT(1:NG),NG,MPI_INTEGER
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

	! adaptive (killing small cluster goes here)
	! First, in the above WINNER call, we know how many pts
	! are there in each cluster i.e., GRP_CT(1:NG)
	I=1	! group index
	DO WHILE (I .LE. NG .AND. NG .GT. 1)
	  ! Now, if the GRP_CT(I) < 2**NATT, we kill grp I by
	  ! copying info from grp NG to grp I
	  IF ( GRP_CT(I) .LT. 2**NATT) THEN
!	  IF ( DBLE(GRP_CT(I)) .LT. 0.02D0*DBLE(NIND) ) THEN
	    GRP_CT(I) = GRP_CT(NG); GRP_CT(NG) = 0
	    XK(I,1:NATT) = XK(NG,1:NATT)
            XKOLD(I,1:NATT) = XKOLD(NG,1:NATT)
	    ! Here, reassign pts in killed grp to 'remaining' grps
	    DO KK=1,NIND
	      IF (IDT(KK).EQ.I) THEN
		CALL WINNER1(NIND,NATT,NG-1,XSTAN(KK,1:NATT)
     $			    ,XK(1:NG-1,1:NATT),IDT(KK),GRP_CT(1:NG-1),IER)
	      ENDIF
	      ! assign pts in grp NG to grp I
	      IF (IDT(KK).EQ.NG) IDT(KK)=I
            END DO
	    NG = NG-1
#if defined(petsc)
            LGRP(1:NG) = GRP_CT(1:NG)
            call MPI_Allreduce(lgrp(1:NG),GRP_CT(1:NG),NG,MPI_INTEGER
     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
	  ELSE
	    I = I+1	! Consider the next grp/cluster
	  ENDIF
	END DO	! end while loop I

	! Update means using pts in cluster:
        CALL UPDATE(NIND,NATT,NG,XSTAN,XK(1:NG,1:NATT),IDT,IER)
        ET=RULE(NG,NATT,XKOLD(1:NG,1:NATT),XK(1:NG,1:NATT))
        XMU = 0d0
	XMU(1:NG,1:NATT) = XK(1:NG,1:NATT)
        IF (ET .LE. EPSL) RETURN
      END DO	! end iter loop T
      
      IER=-41

      END SUBROUTINE ADAPT_KMEANS

c$$$      SUBROUTINE KMEANS(NIND,NATT,NG,X,XMU,IDT,EPSILON,TT,T,IER)
c$$$C      Main subroutine
c$$$      
c$$$       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$       INTEGER T,TT
c$$$       EXTERNAL RANDNUM
c$$$       INTEGER FLAGS(40),FYLENO
c$$$       COMMON /STORE2/ FLAGS,FYLENO
c$$$       DOUBLE PRECISION RANDNUM
c$$$       DIMENSION X(NIND,NATT),XK(NG,NATT),
c$$$     &           XKOLD(NG,NATT),IDT(NIND),
c$$$     &           XSTAN(NIND,NATT),XMU(NG,NATT)
c$$$      IER=0
c$$$c      XMU(1:NG,1:NATT) = 0D0
c$$$c      CALL KSTAND(NIND,NATT,X,XSTAN) 		! normalize sample
c$$$c      CALL KSEED(NIND,NATT,NG,XSTAN,XK,IER)	! pick centroids randomly from normalized sample.
c$$$      XSTAN(1:NIND,1:NATT) = X(1:NIND,1:NATT)	! do not normalize sample X
c$$$c      CALL KSEED(NIND,NATT,NG,XSTAN,XK,IER)	! pick centroids randomly from original sample.
c$$$      XK(1:NG,1:NATT) = XMU(1:NG,1:NATT)	! pick centroids = pseudo-random init of means
c$$$      DO 30 T=1,TT
c$$$      DO 430 KK=1,NG
c$$$        DO 420 LL=1,NATT 
c$$$            XKOLD(KK,LL)=XK(KK,LL)
c$$$420   	  CONTINUE
c$$$430     CONTINUE
c$$$	DO 20 KK=1,NIND
c$$$          CALL WINNER(NIND,NATT,NG,XSTAN,KK,XK,IDT,IER)
c$$$20      CONTINUE
c$$$          CALL UPDATE(NIND,NATT,NG,XSTAN,XK,IDT,IER)
c$$$          ET=RULE(NG,NATT,XKOLD,XK)
c$$$	  XMU(1:NG,1:NATT) = XK(1:NG,1:NATT)
c$$$          IF (ET.LE.EPSILON) GO TO 99
c$$$30    CONTINUE
c$$$      
c$$$      !WRITE (FYLENO,*) 'REACHED MAXIMUM NUMBER OF ',TT,' ITERATIONS'
c$$$      IER=-41
c$$$99    RETURN
c$$$      END

c$$$      SUBROUTINE KSTAND(NIND,NATT,X,XNEW) 
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      DIMENSION X(NIND,NATT),XNEW(NIND,NATT),
c$$$     &          XVAR(NATT),XMU(NATT)
c$$$      DO 200 J=1,NATT
c$$$       XMU(J)=0
c$$$       DO 200 I=1,NIND
c$$$        XMU(J)=XMU(J)+X(I,J)/NIND 
c$$$200   CONTINUE       
c$$$      DO 210 J=1,NATT
c$$$       XVAR(J)=0
c$$$       DO 210 I=1,NIND
c$$$        XVAR(J)=XVAR(J)+(X(I,J)-XMU(J))*
c$$$     &            (X(I,J)-XMU(J))/(NIND-1)
c$$$210   CONTINUE
c$$$      DO 220 J=1,NATT
c$$$       DO 220 I=1,NIND
c$$$         XNEW(I,J)=(X(I,J)-XMU(J))/XVAR(J)
c$$$220   CONTINUE
c$$$      RETURN
c$$$      END

c$$$      SUBROUTINE KSEED(NIND,NATT,NG,XSTAN,XK,IER)               
c$$$c     This Subroutine chooses the initial K seeds (Means of clusters)
c$$$c     for the algorithm. At present they are chosen from data set at 
c$$$c     random.
c$$$
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      INTEGER CHOICE
c$$$      EXTERNAL RANDNUM
c$$$      DOUBLE PRECISION RANDNUM
c$$$      DIMENSION XSTAN(NIND,NATT),XK(NG,NATT)
c$$$      DO 210 I=1,NG
c$$$c        R=RANDNUM()
c$$$c	!CALL RANDOM_NUMBER(R)
c$$$c        R=R*NIND
c$$$c       Convert CHOICE to integer
c$$$c	CHOICE=INT(R)+1
c$$$	CHOICE = max(int(NIND*(dble(I))/dble(NG)),1)
c$$$	DO 200 J=1,NATT
c$$$          XK(I,J)=XSTAN(CHOICE,J) 
c$$$200	CONTINUE
c$$$210   CONTINUE
c$$$      RETURN
c$$$      END

c$$$      SUBROUTINE WINNER(NIND,NATT,NG,XSTAN,KK,XK,IDT,IER)
c$$$c     This subroutine determines the allocation of the KKth point 
c$$$c     ie which mean is closest to the given data point (Euclidean).
c$$$
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      DIMENSION XSTAN(NIND,NATT),XK(NG,NATT),
c$$$     &          IDT(NIND)
c$$$      DO 310 I=1,NG
c$$$        DIST=0
c$$$  	DO 300 J=1,NATT
c$$$          DIST=DIST+(XSTAN(KK,J)-XK(I,J))**2
c$$$300 	CONTINUE
c$$$        IF (I.EQ.1) DISTB=DIST 
c$$$        IF (DIST.LE.DISTB) THEN
c$$$	  IDT(KK)=I
c$$$          DISTB=DIST
c$$$        ENDIF
c$$$310   CONTINUE
c$$$      RETURN
c$$$      END

      SUBROUTINE WINNER1(NIND,NATT,NG,XSTAN,XK,IDT,GRP_CT,IER)
c     This subroutine is exactly same as WINNER subroutine, except
c     it also compute element counts for groups, GRP_CT(1:NG)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XSTAN(NATT),XK(NG,NATT)
      INTEGER :: GRP_CT(NG),IB,IDT

      IB = 1
      DISTB = 0
      
      DO 310 I=1,NG
        DIST=0
  	DO 300 J=1,NATT
          DIST=DIST+(XSTAN(J)-XK(I,J))**2
300 	CONTINUE
        IF (DIST.LE.DISTB) THEN
	  IDT=I
	  IB=I
          DISTB=DIST
        ENDIF
310   CONTINUE
      GRP_CT(IB)=GRP_CT(IB)+1

      RETURN
      END 
          
      SUBROUTINE UPDATE(NIND,NATT,NG,XSTAN,XK,IDT,IER)
      IMPLICIT NONE

      INTEGER :: NIND,NATT,NG,IER
      REAL(8) :: XK(NG,NATT),XSTAN(NIND,NATT),SUMXK(NG,NATT)
      INTEGER :: N(NG),SUMN(NG),II,I,LL,IDT(NIND)
      
        DO 410 II=1,NG
          N(II)=0
          DO 410 LL=1,NATT
            XK(II,LL)=0 
410     CONTINUE
        DO 450 I=1,NIND
          II=IDT(I)
          N(II)=N(II)+1
c         Update rules
          DO 440 LL=1,NATT
            XK(II,LL)=XK(II,LL)+XSTAN(I,LL)
440       CONTINUE
450     CONTINUE

#if defined(petsc)
         SUMXK = XK
         call MPI_Allreduce(SUMXK,XK,NG*NATT
     .              ,MPI_DOUBLE_PRECISION,MPI_SUM
     .              ,MPI_COMM_WORLD,mpierr)

         SUMN = N
         call MPI_Allreduce(SUMN,N,NG
     .              ,MPI_INTEGER,MPI_SUM
     .              ,MPI_COMM_WORLD,mpierr)
#endif

        DO 499 II=1,NG
          DO 499 LL=1,NATT
           IF (N(II).NE.0) THEN
            XK(II,LL)=XK(II,LL)/N(II)       
           ELSE
             RETURN
           ENDIF
499     CONTINUE
      RETURN
      END


      FUNCTION RULE(NG,NATT,XKOLD,XK)
c     This function returns the value used to determine if the algorithm
c     has converged it is a measure of the change in the nodes from iteration 
c     to iteration. 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER R
      DIMENSION XK(NG,NATT),XKOLD(NG,NATT)
      RULE=0.0
      DO 510 KK=1,NATT
        DO 500 R=1,NG
           RULE=RULE+abs(XK(R,KK)-XKOLD(R,KK))
500     CONTINUE
510   CONTINUE
      RETURN
      END


c           FUNCTION RANDNUM()
cC          This is the function called by the program NMM. If you
cC          wish to use your own portable random number generator 
cC          then it should be used in place of this function.

c	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	   COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
c	   IF (RANDTYPE.EQ.1) THEN
c             RANDNUM=RANDOM(IX,IY,IZ)
c	   ELSE
c      WRITE(*,*)'ERROR: As previously described due to random number'
c      WRITE(*,*)'       generator problems features utilising'
c      WRITE(*,*)'       random numbers are unavailable'
c	   STOP
c	   ENDIF
c          RETURN
c	   END

      END MODULE KMEANS_MOD

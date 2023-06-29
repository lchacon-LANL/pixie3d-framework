      MODULE GS_KMEANS_MOD

      USE GM_IO

      USE grid_mpi

      USE MATH, ONLY: GAUSS2,UNIF_RNG

      use OMP
      
      PRIVATE  ! Make everything private
      PUBLIC :: GMM_KMEANS_MULTI_INIT
      PUBLIC :: OPTIMAL_K
      
!********************************************************************************      
!********************************************************************************
! Authors: L. Chacon, G. Chen, T. Nguyen (2020)
! Note (as of 08/24/2020):
! - Main subroutines of K-means module is taken EMMIX.F by David Peel May 1994.
! - The module also includes: implementation of adaptive K-means using
!     similar idea to adaptive EM (kill a cluster if too small) -- NOT GOOD.
! - Determine the number of cluster by SILHOUETTE VALUE -- NOT EFFICIENT, O(N^2).
! - Determine the number of cluster by GAP-STATISTICS -- RELIABLE & EFFICIENT!
! See [Tibshirani2001] for the gap statistics method
!
!
! To do list (as of 10/20/2020): 
!   - Need to perform careful parallelization for routines that support GS-Kmeans.
!   - Need to add particles' weights.
!********************************************************************************
!********************************************************************************

      CONTAINS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    Collection of subroutines for K-MEANS CLUSTERING
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!######################################################################
!   GMM_KMEANS_MULTI_INIT
!######################################################################
!   This subroutines perform K-means (with fixed number of clusters)
!   for 'KM_RUNS' number of times and choose the best run using the 
!   smallest inertia (SSE) value.
      SUBROUTINE GMM_KMEANS_MULTI_INIT(NIND,NATT,NG,MAXIT,KM_TOL
     $      ,WX,X,XMU,XVAR,T,IDT,GRP_COUNT,WTOT,ITOUT,KM_RUNS)
!     This subroutine uses the KMEANS to initialize G-parameters for EM.
!     Approach: Run K-Means clustering for 'KM_RUNS' times & select best cluster
      IMPLICIT NONE

!     Call variables

      INTEGER, INTENT(IN) :: NIND,NATT,MAXIT,KM_RUNS
      INTEGER, INTENT(INOUT) :: NG
      INTEGER, INTENT(OUT) :: ITOUT
      INTEGER, INTENT(INOUT) :: IDT(NIND)
      REAL(8), INTENT(IN) :: WX(NIND),X(NIND,NATT),WTOT,KM_TOL
      REAL(8), INTENT(INOUT):: T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT)
      REAL(8), INTENT(OUT) :: GRP_COUNT(NG)

!     Local variables

      REAL(8) :: RAN,TEMP,KLdist
     $     ,XVAR_K(1:NATT,1:NATT),LSUMX(1:NATT,1:NATT)
     $     ,INERTIA,MIN_INERTIA
      INTEGER :: I,J,K,II,JJ,RANI,IER,NINDT                     !RANI=Random integer
      REAL(8) :: XMU_KM(NG,NATT),LGRP(NG)
      INTEGER :: IDT_KM(NIND)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)    

      REAL(8) :: WS(NIND),XMUS(NATT),XVARS(NATT,NATT)
     $          ,WSUM,LWTOT

      REAL(8) :: XSTAN(NIND,NATT),XMIN(NATT),XMAX(NATT)
     $          ,XRANGE(NATT)

!     Begin program
      
      NINDT = psum((/NIND/))
      
      ! Standardize the particles' velocities X-->XSTAN
      DO I=1,NATT
        XMIN(I)=pmin(MINVAL(X(:,I)))
        XMAX(I)=pmax(MAXVAL(X(:,I)))
        XRANGE(I)=XMAX(I)-XMIN(I)
        XSTAN(:,I)=X(:,I)/XRANGE(I)
cc        XSTAN(:,I)=X(:,I)
      END DO

      MIN_INERTIA = 1.0D8   ! initialize to be large
      GRP_COUNT = 0D0

      IF (NINDT > 0) then
        ! Using K-means clustering to obtain 'good guess' for means
        ! Standard K-Means (performed for 'KM_RUNS' times)
        DO JJ = 1,KM_RUNS
          CALL KMEANS(NIND,NATT,NG,WX,XSTAN,XMU_KM,IDT_KM
     $               ,KM_TOL,MAXIT,ITOUT,IER,INERTIA)
          if (STD_IO) WRITE(*,*) "ITER =",JJ,"SSE =",INERTIA
          if (DIAG_IO) WRITE(FYLENO,*) "ITER =",JJ,"SSE =",INERTIA
          ! Select best run using 'smallest' value for inertia (SSE)
          IF (INERTIA .LE. MIN_INERTIA) THEN
            XMU = XMU_KM
            IDT = IDT_KM
            MIN_INERTIA = INERTIA
          END IF
        END DO
        if (STD_IO) THEN
          WRITE(*,*) "Final clustering, min_SSE =",MIN_INERTIA
          WRITE(*,*)
          WRITE(*,*)
        ENDIF
        if (DIAG_IO) THEN
          WRITE(FYLENO,*) "Final clustering, min_SSE =",MIN_INERTIA
          WRITE(FYLENO,*)
          WRITE(FYLENO,*)
        ENDIF

        ! Once done, calculate number of elements in each cluster/grp
        ! taking into account particle weights
        GRP_COUNT(1:NG) = 0D0
        DO K = 1,NG
          DO II=1,NIND
            IF (IDT(II) .EQ. K) GRP_COUNT(K)=GRP_COUNT(K)+WX(II)
          END DO
        END DO
#if defined(petsc)
        LGRP(1:NG) = GRP_COUNT(1:NG)
        call MPI_Allreduce(LGRP(1:NG),GRP_COUNT(1:NG),NG
     .                    ,MPI_DOUBLE_PRECISION,MPI_SUM
     .                    ,MPI_COMM_WORLD,mpierr)
#endif
        
        ! Here, we calculate total-mean & total-variance
        ! aka sample mean & sample variance matrix
        WSUM=PSUM(WX(1:NIND))
c$$$#if defined(petsc)
c$$$        lwtot = WSUM
c$$$        call MPI_Allreduce(lwtot,WSUM,1,MPI_DOUBLE_PRECISION
c$$$     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
c$$$#endif
        CALL LCALK(NIND,NATT,X,WX,WSUM,XMUS(1:NATT)
     $            ,XVARS(1:NATT,1:NATT)) 

        DO K = 1,NG
          ! Calculate groups' proportions
          T(K) = GRP_COUNT(K)/WTOT     ! Gaussian weight K := proportion of group K
          ! Compute new estimate of 'with-in' covariance matrix for each cluster
          ! This is initial guess for Gaussian cov-matrix K

          ! Rescale the means back
          DO I=1,NATT
            XMU(K,I)=XMU(K,I)*XRANGE(I)
          END DO

          IF ( INT(T(K)*NINDT) .GT. 2**NATT ) THEN
            DO J=1,NATT
              DO I=1,J
                XVAR_K(I,J)=0.0D0
                DO JJ=1,NIND
                  IF ( IDT(JJ) .EQ. K ) THEN
                    XVAR_K(I,J)=XVAR_K(I,J) +
     $                WX(JJ)*(X(JJ,I)-XMU(K,I))*(X(JJ,J)-XMU(K,J))
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
            XVAR(K,:,:) = XVAR_K(:,:)
          ELSE
            XVAR(K,:,:) = XVARS(:,:)/DBLE(NG) ! smaller variance may help converging faster?
          END IF                ! end if ( INT(T(K)*NIND) .GT. 2**NATT ) 
        END DO                  ! end do K = 1,NG
      END IF                    ! end if (NIND > 0)

      END SUBROUTINE GMM_KMEANS_MULTI_INIT
      

!######################################################################
!   This subroutines calculates KMeans optimal K using 
!       Gap-Statistic from Tibshirani, Walther, Hastie 2001
!   We choose the optimal number of components K to be the
!   smallest k such that Gap(k) > Gap(k+1) - 0.5*sd_{k+1} for k=KMIN,...,KMAX
!   This condition is proposed in the paper by Tibshirani, Walther, Hastie 2001
      SUBROUTINE OPTIMAL_K(NIND,NATT,WX,X,NREFS,KMIN,KMAX,KM_TOL,OPT_K)

      IMPLICIT NONE

      !   Call variables
      INTEGER, INTENT(IN) :: NIND,NATT,NREFS,KMIN,KMAX
      REAL(8), INTENT(IN) :: WX(NIND),X(NIND,NATT),KM_TOL
      INTEGER, INTENT(OUT) :: OPT_K
    
      !   Local variables
      INTEGER :: I,J,K,GAP_IDX,REF_SET_METHOD
      REAL(8) :: GAPS(KMAX-KMIN+1),refDisps(NREFS),origDisp
     $     ,mean_refDisps,XREF(NIND,NATT),SD(KMAX-KMIN+1)
     $     ,SQDEV(NREFS),WTOT,WXREF(NIND),LWTOT

      INTEGER :: IDT(NIND),MAXIT,ITOUT,IER,K_VECTOR(KMAX-KMIN+1)
      REAL(8) :: XMU(1:KMAX,1:NATT),INERTIA,A1,A2
    
      ! Code begins
    
      ! Set some parameters
      REF_SET_METHOD = 2    ! 1. Generate each reference feature uniformly over the range of 
                            !     the observed values for that feature.
                            ! 2. Generate the reference features from a uniform distribution 
                            !     over a box aligned with the principal components of the data.
                            ! 3. Generate the reference data set from a normal distribution
                            !     with parameters := samples' mean & samples' covariance 
#if defined(petsc)
      REF_SET_METHOD = 3    ! Method 2 is not MPI-ready
#endif
      
      WTOT = PSUM(WX)

      MAXIT = 200           ! Max iter for K-means
      A2 = 0.0D0            ! Scaling factor for the squared standard dev.
                            ! We find (heuristically) that it is good to set A2=1 or 2
                            ! to capture the right number of components for our synthetic & PIC data sets
      A1 = 1D0 + 1D0/DBLE(NREFS)
    
      WXREF(1:NIND) = 1D0   ! REF. SET. will be i.i.d (equal weights)

      ! Begin the loop to determine best/optimal number of components
      K_VECTOR = 0
      GAP_IDX = 1
      DO K=KMIN,KMAX
        K_VECTOR(GAP_IDX) = K
        refDisps(1:NREFS) = 0D0
        DO I=1,NREFS
          CALL GENERATE_REFERENCE_SET(NIND,NATT,X,XREF,REF_SET_METHOD)
          CALL KMEANS(NIND,NATT,K,WXREF,XREF,XMU(1:K,:),IDT
     $               ,KM_TOL,MAXIT,ITOUT,IER,INERTIA)
          refDisps(I) = INERTIA
        END DO
        CALL KMEANS(NIND,NATT,K,WX,X,XMU(1:K,:),IDT
     $           ,KM_TOL,MAXIT,ITOUT,IER,INERTIA)
        origDisp = INERTIA
        !   Calculate gap-statistics value associated with K-component
        mean_refDisps = SUM(LOG(refDisps))/NREFS
        GAPS(GAP_IDX) = mean_refDisps - LOG(origDisp)
        !   Begin the process of selecting optimal K
        DO I=1,NREFS
          SQDEV(I)=(LOG(refDisps(I)) - mean_refDisps)**2
        END DO
        SD(GAP_IDX)=SQRT(A1*SUM(SQDEV)/NREFS)
        IF (GAP_IDX .GT. 1) THEN
          IF (GAPS(GAP_IDX-1) .GE. GAPS(GAP_IDX)+A2*SD(GAP_IDX)) THEN
            OPT_K = K_VECTOR(GAP_IDX-1)
            EXIT
          END IF
        END IF
        ! increase GAP_IDX by 1 and move on to the next K
        GAP_IDX = GAP_IDX+1
      END DO

      END SUBROUTINE OPTIMAL_K
    
    
!######################################################################
!   GENERATE_REFERENCE_SET
!######################################################################
      SUBROUTINE GENERATE_REFERENCE_SET(NIND,NATT,X,XREF,METHOD)
      ! This subroutine generates the reference distribution for the 
      ! Monte-Carlo GS method. There are 3 ways of doing it:
      ! 1. Generate each reference feature uniformly over the range of 
      !     the observed values for that feature.
      ! 2. Generate the reference features from a uniform distribution 
      !     over a box aligned with the principal components of the data.
      ! 3. Generate the reference data set from a normal distribution
      !     with parameters := samples' mean & samples' covariance

      IMPLICIT NONE

      !   Call variables
      INTEGER, INTENT(IN) :: NIND,NATT,METHOD
      REAL(8), INTENT(IN) :: X(NIND,NATT)
      REAL(8), INTENT(OUT) :: XREF(NIND,NATT)

      !   Local variables
      INTEGER :: I,J,K
      REAL(8) :: XMIN,XMAX,XP(NIND,NATT),DIFF
      REAL(8) :: XMU(NATT),XVAR(NATT,NATT)

      !   Variable for the SVD of X
      REAL(8) :: U(NIND,NATT),S(NATT),VT(NATT,NATT),V(NATT,NATT)
     $   ,ZP(NIND,NATT)
      INTEGER, PARAMETER :: LWORK = 10000   ! set to be large integer?
      INTEGER :: INFO
      REAL(8) :: WORK(LWORK)

      ! Code begins
      IF (METHOD .EQ. 1) THEN
        ! Generate each reference feature uniformly over the range of 
        ! the observed values for that feature.
        DO I=1,NATT
          XMIN = pmin(MINVAL(X(:,I)))
          XMAX = pmax(MAXVAL(X(:,I)))
          DO J=1,NIND
            CALL UNIF_RNG(XREF(J,I))
            XREF(J,I) = DIFF*XREF(J,I)+XMIN
          END DO
        END DO
      ELSE IF (METHOD .EQ. 2) THEN
#if defined(petsc)
        call pstop("generate_reference_set","Broken in parallel")
#endif
        ! Generate the reference features from a uniform distribution 
        ! over a box aligned with the principal components of the data.
        ! Compute SVD of X(1:NIND,1:NATT), note that the DGESVD returns V^T not V
        ! and that on output, the incoming data will be destroyed.
        ! Hence, we make a copy:
        XP = X
        ! Then do SVD on XP
        CALL DGESVD('N','A',NIND,NATT,XP,NIND,S,U,1,VT,NATT
     $          ,WORK,LWORK,INFO)
        V = TRANSPOSE(VT) ! obtain V from transpose of V
        ! Do the transformation: XP=X*V
        XP = 0D0
        CALL DGEMM('N','N',NIND,NATT,NATT,1.0D0,X,NIND,V,NATT
     $          ,0.0D0,XP,NIND)
        DO I=1,NATT
          XMIN = pmin(MINVAL(XP(:,I)))
          XMAX = pmax(MAXVAL(XP(:,I)))
          DIFF = XMAX-XMIN
          !DIFF = 1.0D0
          !XMIN = 0D0
CC          CALL RANDOM_NUMBER(ZP(1:NIND,I))  !Problematic with OpenMP
CC          ZP(1:NIND,I) = DIFF*ZP(1:NIND,I)+XMIN
          DO J=1,NIND
            CALL UNIF_RNG(ZP(J,I))
            ZP(J,I) = DIFF*ZP(J,I)+XMIN
          END DO
        END DO
        ! Finally, transform back via XREF = Z = ZP*V^T
        XREF = 0D0
        CALL DGEMM('N','N',NIND,NATT,NATT,1.0D0,ZP,NIND,VT,NATT
     $          ,0.0D0,XREF,NIND)    
      ELSE IF (METHOD .EQ. 3) THEN
        ! first we compute the 'global' sample's mean & covariance
        CALL GET_SAMPLE_MEAN_COV(NIND,NATT,X,XMU,XVAR)
        ! then use Box-Muller transformation to sample from Normal(mean,covariance)
        CALL BOX_MULLER_GAUSSIAN_SAMPLING(NIND,NATT,XMU,XVAR,XREF)
      END IF  ! End if in METHOD   

      END SUBROUTINE GENERATE_REFERENCE_SET
    
!######################################################################
!   BOX_MULLER_GAUSSIAN_SAMPLING
!######################################################################
      SUBROUTINE BOX_MULLER_GAUSSIAN_SAMPLING(NIND,NATT,XMU,XVAR,XREF)
      ! This subroutine takes parameters mean (XMU) & covariance (XVAR)
      ! and sample 'NIND' points from a Gaussian distributions with
      ! these parameters.... Output in XREF(NIND,NATT)

      USE AAMEM_UTILITIES_MOD, ONLY: MCHOLESKY

      IMPLICIT NONE

      !   Call variables
      INTEGER, INTENT(IN) :: NIND,NATT
      REAL(8), INTENT(IN) :: XMU(NATT),XVAR(NATT,NATT)
      REAL(8), INTENT(OUT) :: XREF(NIND,NATT)

      !   Local variables
      REAL(8), PARAMETER :: PI=3.14159265358979D0
      INTEGER :: I,J,K,II,JJ
      REAL(8) :: L(NATT,NATT),LX(NATT)
      CHARACTER(300) :: FILENAME

      ! Code begins
      DO I=1,NATT
        DO J=1,NIND
          CALL GAUSS2(XREF(J,I))
        END DO
      END DO
        
      ! Compute the lower-triangular matrix L in the Cholesky decomp. of XVAR: 
      !     XVAR = L * L^T
      CALL MCHOLESKY(XVAR,NATT,L)    ! Truong: Here, because of the implementation
                                     ! of MCHOLESKY, the returned L is symmetric
                                     ! so we need to zero-out the upper triangular part of L
      ! Zero out the upper triangular part of L:
      DO II=1,NATT
        DO JJ=II+1,NATT
          L(II,JJ)=0D0
        END DO
      END DO

      !For each data point XREF(I,:), do XREF(I,:) = L*XREF(I,:) + XMU(:)
      DO I=1,NIND
        LX = XREF(I,:)
        CALL DGEMV('N',NATT,NATT,1.0D0,L,NATT,LX,1,0D0,XREF(I,1:NATT),1)
        XREF(I,1:NATT) = XREF(I,1:NATT) + XMU(1:NATT)
      END DO
      ! On output, XREF ~ Normal(XMU,XVAR)


      END SUBROUTINE BOX_MULLER_GAUSSIAN_SAMPLING
    
!######################################################################
!   GET_SAMPLE_MEAN_COV
!######################################################################
      SUBROUTINE GET_SAMPLE_MEAN_COV(NIND,NATT,X,XMU,XVAR)

      IMPLICIT NONE

        !   Call variables
      INTEGER, INTENT(IN) :: NIND,NATT
      REAL(8), INTENT(IN) :: X(NIND,NATT)
      REAL(8), INTENT(OUT) :: XMU(NATT),XVAR(NATT,NATT)

        !   Local variables
      INTEGER :: I,J,K,JJ
      REAL(8) :: WTOT,LSUMX(NATT,NATT),LSUM(NATT)

      ! Code begins
      WTOT = DBLE(NIND)
      XMU = 0D0
      XVAR = 0D0
    
      DO I=1,NIND
        XMU(1:NATT) = XMU(1:NATT)+X(I,1:NATT)
      END DO
#if defined(petsc)
      LSUM = XMU
      call MPI_Allreduce(lsum,XMU,NATT,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      XMU = XMU/WTOT + (1.11d-16)

      DO J=1,NATT
        DO I=1,J
          XVAR(I,J)=0.0D0
          DO JJ=1,NIND
            XVAR(I,J)=XVAR(I,J)+(X(JJ,I)-XMU(I))*(X(JJ,J)-XMU(J)) 
          END DO
          XVAR(J,I)=XVAR(I,J)
        END DO
      END DO
#if defined(petsc)
      lsumx = XVAR
      call MPI_Allreduce(lsumx,XVAR,NATT*NATT,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      XVAR = XVAR/WTOT + (1.11d-16)    

      END SUBROUTINE GET_SAMPLE_MEAN_COV

      
!   #######################################################################
!   This group of subroutines implements the K-means Clustering algorithm.
!   Implemented by David Peel May 1994

!######################################################################
!   SEEDING_INIT_FOR_KMEANS
!######################################################################
      SUBROUTINE SEEDING_INIT_FOR_KMEANS(NIND,NATT,NG,WX,X,C)
      !This subroutine is taken from:
      !https://rosettacode.org/wiki/K-means%2B%2B_clustering#Fortran
      !Implemented based on the paper "kmeans++: the advantages of careful seeding"
      !by David Arthur and Sergei Vassilvitskii (2007)
      !It is tested with gfortran and works well

      IMPLICIT NONE

      !Call variables
      INTEGER, INTENT(IN) :: NIND,NATT,NG
      REAL(8), INTENT(IN) :: WX(NIND),X(NIND,NATT)
      REAL(8), INTENT(OUT) :: C(NG,NATT)

      !     Local variables
      INTEGER :: I,J,K,CHOICE
      REAL(8) :: DIST(NIND,NG),MIN_DIST
     $  ,BEST,D2,TOT,WORK(NIND),W

      ! Code begins
      DO I=1,NIND
        WORK(I) = 1D12
      END DO

      CALL UNIF_RNG(W)
      CHOICE=MIN(INT(W*DBLE(NIND))+1, NIND)
      C(1,1:NATT) = X(CHOICE,1:NATT)

      IF (NG .GT. 1) THEN
        DO K=2,NG
          TOT = 0D0   ! initialize cummulative squared distance
          DO I=1,NIND
            BEST = WORK(I)
            D2 = 0D0
            DO J=1,NATT
              D2=D2+WX(I)*(X(I,J)-C(K-1,J))**2     ! (weighted) Euclidean distance between X_I and C_J
              IF (D2 .GE. BEST) GO TO 10
            END DO
            IF (D2 .LT. BEST) BEST = D2
            WORK(I) = BEST
10          TOT = TOT + BEST
          END DO
          CALL UNIF_RNG(W)
          W = W * TOT
          TOT = 0D0
          DO I=1,NIND
            CHOICE=I
            TOT = TOT + WORK(I)
            IF (TOT .GT. W) GO TO 20
          END DO
20        CONTINUE
          C(K,1:NATT) = X(CHOICE,1:NATT)      
        END DO  ! END DO K
      END IF    ! END IF (NG .GT. 1)

      END SUBROUTINE SEEDING_INIT_FOR_KMEANS
      
      


!######################################################################
!   KMEANS
!######################################################################
      SUBROUTINE KMEANS(NIND,NATT,NG,WX,X,XMU,IDT,EPSILON
     $                 ,TT,T,IER,INERTIA_SUM)
!     Main subroutine - standard K-means for fixed number of cluster
!     K-means can be initialized randomly or with improved seedings (preferred!)

      IMPLICIT NONE

!     Call variables

      INTEGER :: NIND,NATT,NG,T,TT,IER,IDT(NIND)
      REAL(8) :: WX(NIND),X(NIND,NATT),XSTAN(NIND,NATT)
     $          ,XKOLD(NG,NATT),XK(NG,NATT)
     $          ,XMU(NG,NATT),EPSILON

!     Local variables
      
      INTEGER :: LL,KK
      REAL(8) :: INERTIA(NIND),INERTIA_SUM,LSUM,DISTB
    
!     Begin program
      
      IER=0

      XSTAN = X ! do not normalize sample X

      !Pick centroids randomly
CC      CALL KSEED(NIND,NATT,NG,XSTAN,XK,IER)
      !Pick centroids using kmeans++ style with improved seedings
      CALL SEEDING_INIT_FOR_KMEANS(NIND,NATT,NG,WX,XSTAN,XK)

      DO T=1,TT                 ! Begin iterative loop in T

        XKOLD = XK

        DO KK=1,NIND
          CALL WINNER(NIND,NATT,NG,WX,XSTAN,KK,XK,IDT,IER,DISTB)
          INERTIA(KK) = DISTB   ! squared distance of point kk^{th} to its closest cluster
        ENDDO
        
        CALL UPDATE_LOC(NIND,NATT,NG,WX,XSTAN,XK,IDT,IER)
        
        ! Convergence check: also compute the total inertia
        IF (sum(abs(XK-XKOLD)) <= EPSILON) THEN
          ! Regular inertia
          INERTIA_SUM = PSUM(INERTIA)
          XMU = XK
          RETURN
        END IF
      ENDDO

      !WRITE (FYLENO,*) 'REACHED MAXIMUM NUMBER OF ',TT,' ITERATIONS'
      IER=-41
      
      END SUBROUTINE KMEANS

c$$$!######################################################################
c$$$!   KSTAND
c$$$!######################################################################
c$$$!     Normalizing the data to Z-score
c$$$!     For our application of using A-AMEM, we may NOT need to use this
c$$$      SUBROUTINE KSTAND(NIND,NATT,WX,X,XNEW) 
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER :: J,I,NIND,NATT
c$$$      REAL(8) :: X(NIND,NATT),XNEW(NIND,NATT),XVAR(NATT),XMU(NATT)
c$$$     $          ,WX(NIND),WTOT
c$$$      WTOT = SUM(WX(1:NIND))
c$$$      DO 200 J=1,NATT
c$$$        XMU(J)=0
c$$$        DO 200 I=1,NIND
c$$$          XMU(J)=XMU(J)+WX(I)*X(I,J)/WTOT 
c$$$200     CONTINUE       
c$$$      DO 210 J=1,NATT
c$$$        XVAR(J)=0
c$$$        DO 210 I=1,NIND
c$$$          XVAR(J)=XVAR(J)+WX(I)*
c$$$     $            (X(I,J)-XMU(J))*(X(I,J)-XMU(J))/(WTOT)
c$$$210     CONTINUE
c$$$      DO 220 J=1,NATT
c$$$        DO 220 I=1,NIND
c$$$          XNEW(I,J)=(X(I,J)-XMU(J))/XVAR(J)
c$$$220     CONTINUE
c$$$      RETURN
c$$$      END SUBROUTINE KSTAND
    
!######################################################################
!   KSEED
!######################################################################
      SUBROUTINE KSEED(NIND,NATT,NG,XSTAN,XK,IER)               
!   This Subroutine chooses the initial K seeds (Means of clusters)
!   for the algorithm. At present they are chosen from data set at 
!   random. Not recommend to use. Should use improved seedings
!   See subroutine SEEDING_INIT_FOR_KMEANS

      IMPLICIT NONE
      INTEGER :: CHOICE,NIND,NATT,NG,IER,I

      REAL(8) :: XSTAN(NIND,NATT),XK(NG,NATT),R

      DO I=1,NG
        CALL UNIF_RNG(R)
        ! Convert CHOICE to integer
        R=R*NIND
        CHOICE=MAX(INT(R),1)
        XK(I,:)=XSTAN(CHOICE,:) 
      ENDDO
      RETURN
      END
    
!######################################################################
!   WINNER
!######################################################################
      SUBROUTINE WINNER(NIND,NATT,NG,WX,XSTAN,KK,XK,IDT,IER,DISTB)
!     This subroutine determines the allocation of the KKth point 
!     ie which mean is closest to the given data point (Euclidean).

      IMPLICIT NONE
      INTEGER :: NIND,NATT,NG,IDT(NIND),IER,I,J,KK
      
      REAL(8) :: XSTAN(NIND,NATT),XK(NG,NATT),WX(NIND),DISTB,DIST

      DISTB=DIST
      
      DO 310 I=1,NG
        DIST=0D0
        DO 300 J=1,NATT
          DIST=DIST+WX(KK)*(XSTAN(KK,J)-XK(I,J))**2
300     CONTINUE
cc        IF (I.EQ.1) DISTB=DIST 
        IF (DIST.LE.DISTB) THEN
          IDT(KK)=I
          DISTB=DIST
        END IF
310   CONTINUE
      RETURN
      END

!######################################################################
!   UPDATE_LOC
!######################################################################
!     This subroutine update the new location for means in K-means alg.
      SUBROUTINE UPDATE_LOC(NIND,NATT,NG,WX,XSTAN,XK,IDT,IER)
      IMPLICIT NONE

      INTEGER :: NIND,NATT,NG,IDT(NIND),II,I,LL,IER
      REAL(8) :: XK(NG,NATT),XSTAN(NIND,NATT),SUMN(NG),N(NG)
     $          ,WX(NIND),SUMXK(NG,NATT)

      N = 0d0
      XK = 0d0
  
      DO I=1,NIND
        II=IDT(I)
        N(II)=N(II)+WX(I)
        XK(II,:)=XK(II,:)+WX(I)*XSTAN(I,:)
      ENDDO

#if defined(petsc)
      SUMXK = XK
      call MPI_Allreduce(SUMXK,XK,NG*NATT
     .              ,MPI_DOUBLE_PRECISION,MPI_SUM
     .              ,MPI_COMM_WORLD,mpierr)

      SUMN = N
      call MPI_Allreduce(SUMN,N,NG
     .              ,MPI_DOUBLE_PRECISION,MPI_SUM
     .              ,MPI_COMM_WORLD,mpierr)
#endif

      DO II=1,NG
        IF (N(II) > 0d0) THEN
          XK(II,:)=XK(II,:)/N(II)       
        ELSE
          RETURN
        ENDIF
      ENDDO
      
      RETURN
      END

c$$$!######################################################################
c$$$!   RULE
c$$$!######################################################################
c$$$      FUNCTION RULE(NG,NATT,XKOLD,XK)
c$$$!     This function returns the value used to determine if the algorithm
c$$$!     has converged it is a measure of the change in the nodes from iteration 
c$$$!     to iteration. 
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER :: NG,NATT,R,KK
c$$$      REAL(8) :: XK(NG,NATT),XKOLD(NG,NATT),RULE
c$$$
c$$$      RULE= sum(abs(XK-XKOLD))
c$$$
c$$$      RETURN
c$$$      END

c$$$!######################################################################
c$$$!   FUNCTION RANDNUM()
c$$$!######################################################################
c$$$           FUNCTION RANDNUM()
c$$$!C          This is the function called by the program NMM. If you
c$$$!C          wish to use your own portable random number generator 
c$$$!C          then it should be used in place of this function.
c$$$
c$$$	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$	   COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
c$$$	   IF (RANDTYPE.EQ.1) THEN
c$$$             RANDNUM=RANDOM(IX,IY,IZ)
c$$$	   ELSE
c$$$      WRITE(*,*)'ERROR: As previously described due to random number'
c$$$      WRITE(*,*)'       generator problems features utilising'
c$$$      WRITE(*,*)'       random numbers are unavailable'
c$$$	   STOP
c$$$	   ENDIF
c$$$          RETURN
c$$$	   END
c$$$!C   End of subroutine group for KMEANS-CLUSTERING
c$$$!C   #######################################################################

!c     LCALK
!c     A copy yof LCALK from SEMGM.F
!##################################################################
      SUBROUTINE LCALK(NIND,NATT,X,W,WSUM,XMU,XVAR)
!C     Compute estimate of mean, covariance matrix and mixing proportion
!C     for a given Gaussian

      IMPLICIT NONE

!c     Call variables

      INTEGER, INTENT(IN)    :: NIND,NATT
      REAL(8), INTENT(IN)    :: X(NIND,NATT),W(NIND),WSUM
      REAL(8), INTENT(INOUT) :: XMU(NATT), XVAR(NATT,NATT)

!c     Local variables

      INTEGER :: I,J,JJ
      REAL(8) :: SUM(NATT),LSUM(NATT),LSUMX(NATT,NATT)

!C     Compute new estimates of group means (XMU)

      SUM =0.0
      DO J=1,NATT
         DO JJ=1,NIND
            SUM(J)=SUM(J)+X(JJ,J)*W(JJ)
         END DO
      END DO
#if defined(petsc)
      LSUM = SUM
      call MPI_Allreduce(lsum,SUM,NATT,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      XMU=SUM/WSUM

!C     Compute new estimate of covariance matrix for each group

      XVAR=0.0
      DO J=1,NATT
        DO I=1,J
          DO JJ=1,NIND
            XVAR(I,J)=XVAR(I,J)
     .               +(X(JJ,I)-XMU(I))*(X(JJ,J)-XMU(J))*W(JJ)
          END DO
          XVAR(J,I)=XVAR(I,J)
        END DO
      END DO
#if defined(petsc)
      lsumx = XVAR
      call MPI_Allreduce(lsumx,XVAR,NATT*NATT,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      XVAR = XVAR/WSUM + (1.11d-16)

      END SUBROUTINE LCALK

    
    
      !************************************************************************************************
      ! END OF GS--K-MEANS MODULE
      !************************************************************************************************
      END MODULE GS_KMEANS_MOD

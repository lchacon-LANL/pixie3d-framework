c     module gm_init
c     ##################################################################
      module gm_init

      use kmeans_mod

      use gs_kmeans_mod

      integer,private :: mpierr

      REAL(8)                              :: TS
      REAL(8), DIMENSION(:  ), ALLOCATABLE :: XMUS
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: XVARS
CGMM_OMP!$OMP THREADPRIVATE(TS,XMUS,XVARS)

      INTEGER,PRIVATE :: NINDT  ! Global number of observations
CGMM_OMP!$OMP THREADPRIVATE(NINDT)

      LOGICAL :: GMM_RANDOM_INIT=.false.           !Whether we use RANDOM GMM init (else GRP8)
     &          ,GMM_KMEANS_INIT=.false.           !Whether we use KMEANS GMM (on top of RND or GRP8)
     &          ,GMM_GS_KMEANS_MULTI_INIT=.false.  !Whether we use GS-MULTI-KMEANS GMM (on top of RND or GRP8)

      INTEGER :: KM_RUNS    = 10  ! number of K-means trials after optimal NG is obtained
     &          ,KM_REF_SET = 3   ! number of ref. sets to be used in gap-statistics
     &          ,KM_NG_MIN  = 2   ! min # of Gaussians
     &          ,KM_NG_MAX  = 10  ! max # of Gaussians
     &          ,KM_NG_ADJ  = 2   ! to avoid underestimate # of Gaussians

c$$$      namelist /gmm_init/ gmm_random_init,gmm_kmeans_init
c$$$     .                   ,gmm_gs_kmeans_multi_init,km_runs,km_ref_set
c$$$     .                   ,km_ng_min,km_ng_max,km_ng_adj

      contains

c     DESTROY_GM
c     ##################################################################
      SUBROUTINE DESTROY_GM()

      IMPLICIT NONE

      DEALLOCATE(XMUS,XVARS)

      END SUBROUTINE DESTROY_GM

c     INIT_GM
!     #########################################################################
      SUBROUTINE INIT_GM(NIND,NATT,NG,X,WX,IDT,XMU,XVAR,T,WTOT,DV,V)

      IMPLICIT NONE

c    Call variables

      INTEGER, INTENT(IN)    :: NIND,NATT
      INTEGER, INTENT(INOUT) :: NG,IDT(NIND)
      REAL(8), INTENT(IN)    :: X(NIND,NATT),WX(NIND)
      REAL(8), INTENT(OUT)   :: XMU(NG,NATT),XVAR(NG,NATT,NATT),T(NG)
     .                         ,WTOT,DV(NG),V(NG,NATT,NATT)

c    Local variables

      REAL(8) :: RAN            !Random number (0,1)
     $          ,TEMP
     $          ,WSUM(NG)    !Sum of weights of a group
     $          ,KOK,K0,Np,Nd
     $          ,GDV,GV(NATT,NATT),LWTOT,W(NIND,NG)

      INTEGER :: I,J,K,KK,II,RANI
      real(8) :: r(NATT,NATT),sgn(3),NPARBY2
      INTEGER :: info,IER

      ! Variables for GS-kmeans
      INTEGER :: OPT_NG,KMITER,GRP_COUNT(NG)
      REAL(8) :: NK(NG)

c     Begin program

#if defined(petsc)
      call MPI_Allreduce(NIND,NINDT,1,MPI_INTEGER
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
      NINDT = NIND
#endif

      ALLOCATE(XMUS(NATT),XVARS(NATT,NATT))

      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)

      if (NINDT < NPARBY2) then
        write (*,*) "Too few particles per cell. Aborting..."
        STOP
      endif

      T    = 1d0/NG
      XMU  = 0d0
      XVAR = 0d0

      W(1:NIND,1) = WX(1:NIND)

      WTOT = sum(WX) !adding all particle weights for calculating the total variance 
#if defined(petsc)
      lwtot = WTOT
      call MPI_Allreduce(lwtot,WTOT,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      WSUM = WTOT

      if (GMM_RANDOM_INIT) then

CGMM_OMP!$OMP MASTER
        if (my_rank_gmm==0.and.std_io)
     .       PRINT*,"Performing RANDOM init"
CGMM_OMP!$OMP END MASTER

        call GMM_INIT_RANDOM(NIND,NATT,NG,W,WSUM,X,XMU,XVAR,T)

      elseif (GMM_KMEANS_INIT) THEN

CGMM_OMP!$OMP MASTER
        if (my_rank_gmm==0.and.std_io)
     .       PRINT*,"Performing GS-Kmeans init"
CGMM_OMP!$OMP END MASTER

c      Perform KMEANS adaptive refinement
c      Option to use GS--Multi-Kmeans initialization
        ! Perform Gap-statistics to obtain optimal number of
        ! components
        if (GMM_GS_KMEANS_MULTI_INIT) THEN
          CALL OPTIMAL_K(NIND,NATT,WX,X,KM_REF_SET
     $                ,KM_NG_MIN,KM_NG_MAX,OPT_NG)

          NG = MIN(OPT_NG + KM_NG_ADJ,NG)
          if (STD_IO) write (*,*) "Optimal Kmeans # gaussians",NG
          if (DIAG_IO) write (FYLENO,*) "Optimal Kmeans # gaussians",NG
c$$$        else
c$$$          NG = NG
        end if

        ! Given the 'optimal' no. of component, do multi-Kmeans & select the best:
        CALL GMM_KMEANS_MULTI_INIT(NIND,NATT,NG,500,WX,X
     $          ,XMU(1:NG,:),XVAR(1:NG,:,:)
     $          ,T(1:NG),IDT,NK(1:NG),WTOT,KMITER,KM_RUNS)

        if (NG<size(T)) then
          T(NG+1:) = 0d0
          XMU(NG+1,:) = 0d0
          XVAR(NG+1,:,:) = 0d0
        endif
        
c$$$      ELSE
c$$$        CALL GMM_INIT_KMEANS(NIND,NATT,NG,X,XMU(1:NG,1:NATT)
c$$$     $          ,XVAR(1:NG,1:NATT,1:NATT),T(1:NG),WTOT,IDT)

      else

CGMM_OMP!$OMP MASTER
        if (my_rank_gmm==0.and.std_io) PRINT*,"Performing GRP8 init"
CGMM_OMP!$OMP END MASTER

        !Default: group of eight init
        call GMM_INIT_GRP8  (NIND,NATT,NG,W,WSUM,X,XMU,XVAR,T)

      endif

C    !Invert covariance matrices, returned by GS-KMeans init

      IER = 0

      CALL GDET(NATT,NG,XVAR(1:NG,:,:),V(1:NG,:,:),DV(1:NG),IER)

      IER=0

C     Calculate total data size when data is weighted ie weighted
C     likelihood and sum of weights does not equal one.

      WTOT=sum(WX)
#if defined(petsc)
      lwtot = WTOT
      call MPI_Allreduce(lwtot,WTOT,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      ! If (adaptive/GS)--K-means is used then NG is changed.
      ! Here, 
      !   T(1:NG)         stores the the initial weights
      !   XMU(1:NG,:)     stores the the initial means
      !   XVAR(1:NG,:,:)  stores the the initial cov-matrices

      END SUBROUTINE INIT_GM

c     GMM_INIT_RANDOM
c     ##################################################################
      SUBROUTINE GMM_INIT_RANDOM(NIND,NATT,NG,W,WSUM,X,XMU,XVAR,T)

      IMPLICIT NONE

c     Call variables

      INTEGER, INTENT(IN)  :: NIND,NATT,NG
      REAL(8), INTENT(IN)  :: X(NIND,NATT),WSUM(NG),W(NIND,NG)
      REAL(8), INTENT(OUT) :: XMU(NG,NATT),XVAR(NG,NATT,NATT),T(NG)
      
c     Local variables

      REAL(8) :: RAN            !Random number (0,1)
     $          ,TEMP
     $          ,KOK,K0,Np,Nd
     $          ,GDV,GV(1:NATT,1:NATT),LWTOT
      INTEGER :: I,J,K
     $          ,RANI,RCORE
      
c     Begin program

      if (my_rank_gmm == 0.and.std_io) then
        write (*,*) "Performing RANDOM GMM initialization"
      endif
      
c     Compute initial single-Gaussian moments
      
      K = 1
      CALL LCALK(NIND,NATT,X,W(1:NIND,K),WSUM(K),XMU(K,1:NATT)
     .          ,XVAR(K,1:NATT,1:NATT))

      DO 450 K=1,1 !NG (ASK Gy)
        IF (STD_IO) WRITE (*      ,445) K
        IF (DIAG_IO) WRITE (FYLENO,445) K
 445    FORMAT (/2X,'Global mean, covariance matrix  ',I2)
        IF (STD_IO) WRITE (*     ,415) 'XMU',(XMU(K,J),J=1,NATT)
        IF (DIAG_IO)WRITE (FYLENO,415) 'XMU',(XMU(K,J),J=1,NATT)
        DO 450 J=1,NATT
          IF (STD_IO) WRITE (*     ,415) 'XVAR',(XVAR(K,I,J),I=1,J)
          IF (DIAG_IO)WRITE (FYLENO,415) 'XVAR',(XVAR(K,I,J),I=1,J)
 450  CONTINUE

 415  FORMAT (2X,5G30.12)

c     Make copies of global prop, mean & covariance for conservation check
      
      TS = 1D0
      XMUS(1:NATT) = XMU(1,1:NATT)
      XVARS(1:NATT,1:NATT) = XVAR(1,1:NATT,1:NATT)

c     Initialize weights:

      T(1:NG) = 1.0/DBLE(NG)

c     Initialize variances according to global variance

c$$$      !Compute inverse covariance matrix
c$$$      K = 1
c$$$      call MDETINV(NATT,XVAR(K,1:NATT,1:NATT),GDV,GV(1:NATT,1:NATT))
c$$$
c$$$      KOK = KURTOSIS(NIND,NATT,X(1:NIND,1:NATT),W(1:NIND,K)
c$$$     $              ,XMU(K,1:NATT),WSUM(K),GV(1:NATT,1:NATT))
c$$$
c$$$      Np = DBLE(NINDT)
c$$$      Nd = DBLE(NATT)
c$$$      K0 = EKurt(Np,Nd)

c$$$      IF(KOK>K0) then
c$$$c     common mu
c$$$         DO K = 2,NG
c$$$            XMU(K,1:NATT) = XMU(1,1:NATT)
c$$$         END DO
c$$$         CALL RANDOM_NUMBER(RAN)
c$$$         XVAR(2,:,:) = XVAR(1,:,:)*2d0*RAN
c$$$         XVAR(1,:,:) = XVAR(1,:,:)*2d0*(1d0-RAN)
c$$$      ELSE
c$$$         XVAR(1,:,:) = XVAR(1,:,:)*T(1)
         DO K = 2,NG
            XVAR(K,1:NATT,1:NATT) = XVAR(1,1:NATT,1:NATT)
         END DO
c$$$       END IF

c     Initialize means randomly

      if (NINDT > 0) then

#if defined(petsc)
        DO K = 1,NG
          !Pick random processor
          call unif_rng(RAN)    !Works in parallel with SPRNG
          if (my_rank_gmm == 0) RCORE = floor(nproc_gmm*RAN)
          call MPI_Bcast(RCORE,1,MPI_INTEGER,0
     .                  ,MPI_COMM_WORLD,mpierr)
          if (my_rank_gmm == RCORE) then
           !Pick specified particle locations (not parallel)
            RANI = max(int(NIND*(dble(K))/dble(NG)),1)              
            XMU(K,1:NATT) = X(RANI,1:NATT)
          endif
          call MPI_Bcast(XMU(K,1:NATT),NATT,MPI_DOUBLE_PRECISION,RCORE
     .                  ,MPI_COMM_WORLD,mpierr)
        END DO
#else
        DO K = 1,NG
c$$$          !Pick random particle locations (not parallel)
c$$$          call unif_rng(RAN)  !Works in parallel
c$$$          RANI = ceiling(NIND*RAN)
          !Pick specified particle locations (not parallel)
          RANI = max(int(NIND*(dble(K))/dble(NG)),1)              
          XMU(K,1:NATT) = X(RANI,1:NATT)
        END DO
#endif
      endif

      END SUBROUTINE GMM_INIT_RANDOM

c     GMM_INIT_GRP8
c     ##################################################################
      SUBROUTINE GMM_INIT_GRP8(NIND,NATT,NG,W,WSUM,X,XMU,XVAR,T)

      IMPLICIT NONE

c     Call variables

      INTEGER, INTENT(IN)  :: NIND,NATT,NG
      REAL(8), INTENT(IN)  :: X(NIND,NATT),WSUM(NG),W(NIND,NG)
      REAL(8), INTENT(OUT) :: XMU(NG,NATT)
     $                       ,XVAR(NG,NATT,NATT),T(NG)
      
c     Local variables

      INTEGER :: K,KK,I,J
      real(8) :: r(NATT,NATT),sgn(3)
      integer  info
      
c     Begin program

      if(mod(NG,8)/=0) then
        print *, "Error Msg: NG=",NG," is not multiple of 8. Quit."
        stop
      end if
       
      if (my_rank_gmm == 0.and.std_io) then
CGMM_OMP!$OMP CRITICAL
        write (*,*) "Performing GR8 GMM initialization"
CGMM_OMP!$OMP END CRITICAL
      endif
      
c     Calculate initial Gaussian moments
      
      K = 1
      CALL LCALK(NIND,NATT,X,W(1:NIND,K),WSUM(K),XMU(K,1:NATT)
     .          ,XVAR(K,1:NATT,1:NATT))

      DO K = 2,NG
        XMU(K,:) = XMU(1,:)
        XVAR(K,1:NATT,1:NATT) = XVAR(1,1:NATT,1:NATT)/8d0 !smaller initial variance may perform better.
      END DO

c     Make copies of global prop, mean & covariance for conservation checks
      
      TS = 1D0
      XMUS(1:NATT) = XMU(1,1:NATT)
      XVARS(1:NATT,1:NATT) = XVAR(1,1:NATT,1:NATT)

c     Assign Gaussians in groups of 8 according to principal axes
      
      if (NIND > 0) then
        r(1:NATT,1:NATT) = XVAR(1,1:NATT,1:NATT)
      
        call r8po_fa( NATT, r, info)
        sgn = -1
        KK = 1
        do i=1,2
          sgn(1) = -sgn(1);
          do j=1,2
            sgn(2) = -sgn(2);
            do k=1,2
              sgn(3) = -sgn(3);      
              XMU(KK,1:NATT)= sgn(1:NATT)
              KK = KK+1
            end do
          end do
        end do
cc        print *,"r=",r
        
        DO K = 1,NG
          XMU(K,1:NATT) = XMUS(1:NATT)
     .                  + matmul(XMU(K,1:NATT),r(1:NATT,1:NATT))
cc          print *,"mean=",K,XMU(K,1:NATT)
        END DO
      endif

c     Initialize Gaussian weights:

      T(1:NG) = 1d0/DBLE(NG)

c     IO

      DO K=1,1 !NG
        IF (STD_IO) WRITE (*,445) K
        IF (DIAG_IO) WRITE (FYLENO,445) K
        IF (STD_IO) WRITE (*,415) (XMU(K,J),J=1,NATT)
        IF (DIAG_IO)WRITE (FYLENO,415) (XMU(K,J),J=1,NATT)
        DO J=1,NATT
          IF (STD_IO) WRITE (*,415) (XVAR(K,I,J),I=1,J)
          IF (DIAG_IO)WRITE (FYLENO,415) (XVAR(K,I,J),I=1,J)
        ENDDO
      ENDDO

 445  FORMAT (/2X,'Global mean, covariance matrix  ',I2)
 415  FORMAT (2X,5G30.12)
      
      END SUBROUTINE GMM_INIT_GRP8

      end module gm_init

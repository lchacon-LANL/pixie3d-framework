c     module gm_io
c     ##################################################################
      module gm_io

#if defined(m_p_i) || defined(petsc)
      use mpi
#endif

      integer,private :: mpierr

      INTEGER :: my_rank_gmm=0,nproc_gmm=1

      LOGICAL :: STD_IO=.true.                     !Whether we output to STD output
     &          ,DIAG_IO=.false.                   !Whether we output diagnostics to file
CGMM_OMP!$OMP THREADPRIVATE(STD_IO,DIAG_IO)

      INTEGER :: FYLENO                            !File id for GMM output file
     &          ,ML_IO                             !File id for ML output file
CGMM_OMP!$OMP THREADPRIVATE(FYLENO,ML_IO)

      contains

c     INIT_EM_IO
c     ##################################################################
      SUBROUTINE INIT_EM_IO(diagio,stdio,id)
      IMPLICIT NONE

      character(*) :: id
      logical :: diagio,stdio
     
C     MPI init

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank_gmm, mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc_gmm  , mpierr)
#endif
      
c     Set diagnostic files
 
      diag_io = diagio.and.(my_rank_gmm==0)
      std_io  = stdio.and.(my_rank_gmm==0)
      
      if (diag_io) CALL set_diag_files(id)

      END SUBROUTINE INIT_EM_IO

c     CLOSE_EM_FILES
c     ##################################################################
      SUBROUTINE CLOSE_EM_FILES()
      IMPLICIT NONE
       
      CALL close_diag_files

      END SUBROUTINE CLOSE_EM_FILES

c     set_diag_files
c     ########################################################################
      subroutine set_diag_files(id)

      Implicit none

      character(*) :: id

      character(100) :: file

      file = "GMMEM-diag-"//trim(id)//".log"
      FYLENO=22
      OPEN (UNIT=FYLENO,FILE=trim(file),STATUS = 'UNKNOWN') 

      file = "LF-diag-"//trim(id)//".log"
      ML_IO = 31
      OPEN (UNIT=ML_IO,FILE=trim(file),STATUS = 'UNKNOWN')
      write (ML_IO,*) "#ITER  Log-likelihood  Message-length  T(1:K)"

      end subroutine set_diag_files

c     close_diag_files
c     ########################################################################
      Subroutine close_diag_files

      Implicit none

      CLOSE(FYLENO)
      CLOSE(ML_IO)

      end subroutine close_diag_files

c     FINOUT
c     ##################################################################
      SUBROUTINE FINOUT(NIND,NATT,NCOV,NG,XUU,XMU,XVAR,T,IDT,XCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FLAGS(40)
      EXTERNAL AUTOPARTITION
      EXTERNAL RANDNUM,DETERRANDOM
      DOUBLE PRECISION RANDNUM
      DIMENSION XVAR(NG,NATT,NATT),T(NG),
     &         XMU(NG,NATT),IDT(NIND),XUU(NG),N(NG),
     &         XCC(NG),ISUM(NG)
      
CC   local arrays
      real(8) :: xvar_l(3,3),xmu_l(3) !limited to 3 dims

C     FINAL OUTPUT Section

C     Output parameter estimates to file
      IF (FLAGS(7).GT.0) THEN
      WRITE (FYLENO,26)
26    FORMAT (/2X,'Estimated Nu for each component')
        WRITE (FYLENO,35) (XUU(K),K=1,NG)
      ENDIF
      IF (FLAGS(20).EQ.1) THEN
       DO 1030 K=1,NG
        WRITE (FYLENO,35) (XMU(K,J),J=1,NATT)
        IF (NCOV.EQ.1) THEN
          DO 1040 J=1,NATT
            WRITE (FYLENO,45) (XVAR(1,I,J),I=1,J)
1040      CONTINUE
        ELSE
          DO 1041 J=1,NATT
            WRITE (FYLENO,45) (XVAR(K,I,J),I=1,J)
1041      CONTINUE
        ENDIF
1030   CONTINUE

1031   CONTINUE
       WRITE (FYLENO,45) (T(I),I=1,NG)
      ENDIF
      DO 1180 K=1,NG
1180    N(K)=0
      DO 1185 I=1,NIND
        K=IDT(I)
        IF (K.EQ.0) GO TO 1185
          N(K)=N(K)+1
1185  CONTINUE
#if defined(petsc)
      ISUM = N
      call MPI_Allreduce(ISUM,N,NG,MPI_INTEGER
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      IF (NG.GT.1) THEN
      WRITE (FYLENO,1187)
1187  FORMAT (/2X,'Number assigned to each component')
      WRITE (FYLENO,1189) (N(K),K=1,NG)
1189   FORMAT (2X,10I6)
      WRITE (FYLENO,1191)
1191   FORMAT (/2X,'Estimate of mixing proportion for each component')
      WRITE (FYLENO,1193) (T(K),K=1,NG)
1193   FORMAT (2X,10F10.6)
c$$$C     Compute estimates of correct allocation rates
c$$$      CC=0.0
c$$$      DO 1195 K=1,NG
c$$$        XCC(K)=XCC(K)/(NIND*T(K))
c$$$        CC=CC+T(K)*XCC(K)
c$$$1195  CONTINUE
c$$$      WRITE (FYLENO,1197)
c$$$1197   FORMAT (/2X,'Estimates of correct allocation rates for ',
c$$$     &        'each component')
c$$$      WRITE (FYLENO,1193) (XCC(K),K=1,NG)
c$$$      WRITE (FYLENO,1198) CC
c$$$1198   FORMAT (/2X,'Estimate of overall correct allocation rate ',
c$$$     &        F7.3)
      ENDIF

!     accumulate local arrays
      xmu_l = 0
      xvar_l= 0
      WRITE (FYLENO,25)
25    FORMAT (/2X,'Estimated mean (as a row vector) for each component')
      DO 30 K=1,NG
         do J=1,NATT
            xmu_l(J) = xmu_l(J) + xmu(K,J)*T(K)
         end do
30      WRITE (FYLENO,35) (XMU(K,J),J=1,NATT)
35      FORMAT (2X,5G13.5)
C     Test if a common covariance matrix is specified (NCOV = 1)
      IF ((NCOV.EQ.1).OR.(NCOV.EQ.3)) THEN
       IF (NCOV.EQ.1) THEN
        WRITE (FYLENO,37)
37      FORMAT (/2X,'Estimated common covariance matrix ')
       ELSE
        WRITE (FYLENO,38)
38      FORMAT (/2X,'Estimated common diagonal covariance matrix ')
       ENDIF
        DO 40 J=1,NATT
40        WRITE (FYLENO,45) (XVAR(1,I,J),I=1,J)
45      FORMAT (5X,5G14.6)
      ELSEIF ((NCOV.EQ.2).OR.(NCOV.EQ.4)) THEN
        DO 50 K=1,NG
         IF (NCOV.EQ.2) THEN
          WRITE (FYLENO,46) K
46        FORMAT (/2X,'Estimated covariance matrix for component ',I2)
         ELSE
          WRITE (FYLENO,48) K
48    FORMAT(/2X,'Estimated diagonal covariance matrix for component '
     &,I2)
      ENDIF

          DO 50 J=1,NATT
            WRITE (FYLENO,47) (XVAR(K,I,J),I=1,J)
47      FORMAT (5X,5G20.12)
 50   CONTINUE

!     first term
      DO K=1,NG
         DO J=1,NATT      
            DO I=1,J
               xvar_l(I,J) = xvar_l(I,J) + XVAR(K,I,J)*T(K)
            END DO
         END DO
      END DO
      
!     2nd term
      DO K=1,NG
         DO J=1,NATT
            DO I=1,J
               xvar_l(I,J) = xvar_l(I,J)+xmu(K,I)*xmu(K,J)*T(K)
            END DO
         END DO
      END DO

!     3nd term
      DO J=1,NATT
         DO I=1,J
            xvar_l(I,J) = xvar_l(I,J)-xmu_l(I)*xmu_l(J)
         end do
      end do

      ENDIF

!     write the local arrays
      WRITE (FYLENO,55) (xmu_l(J),J=1,NATT)
 55   FORMAT (2X,"total mu and var:",3G20.12)
      DO J=1,NATT
         WRITE (FYLENO,47) (xvar_l(I,J),I=1,J)         
      END DO
cc      if(NG>1) stop
99    RETURN
      END SUBROUTINE FINOUT

c     OUTLOOP
c     ##################################################################
      SUBROUTINE OUTLOOP(NIND,NATT,NG,XMU,DV,T,NCOV,IOUNT,XLOGL,MLENG
     .                  ,W,IDT,X)
c      This subroutine displays all the relevant information from the
c      EM algorithm applied to the best partition.
      implicit double precision (a-h,o-z)

      DIMENSION DV(:),T(NG),W(NIND,NG),
     &          XCC(NG),XLOGL(:),XMU(NG,NATT),
     &          IDT(NIND),X(NIND,NATT)

      REAL(8) :: MLENG(:)

      IF (IOUNT>2) then
        write(FYLENO,*)
        write(FYLENO,*)

        WRITE (FYLENO,9550) IOUNT-1
 9550   FORMAT(2X,"Final iteration:",I5)
        WRITE (FYLENO,9551) XLOGL(1)
 9551   FORMAT(2X,"Initial LOG-LIKELIHOOD:",2F24.16)
        WRITE (FYLENO,9552) XLOGL(IOUNT-1),XLOGL(IOUNT-2)
 9552   FORMAT(2X,"Last two LOG-LIKELIHOOD:",2F24.16)
      
        WRITE (FYLENO,9553) MLENG(1)
 9553   FORMAT(2X,"Initial Message-length:",2F24.16)
        WRITE (FYLENO,9554) MLENG(IOUNT-1),MLENG(IOUNT-2)
 9554   FORMAT(2X,"Last two Message-Length:",2F24.16)
      ENDIF

      WRITE(FYLENO,*)

      WRITE(FYLENO,*) 'Final NG = ', NG

      WRITE (FYLENO,1105) IOUNT,XLOGL(IOUNT)
1105    FORMAT (/2X,'In loop ',I3,' log likelihood is ',F15.3)

      WRITE (FYLENO,*) 'Observation | mixture log density|',
     &'Grp1 Grp2 ..etc...|'

      DO 1160 I=1,NIND
        WRITE (FYLENO,1170) I,(W(I,K),K=1,NG)
 1160 CONTINUE

 1170 FORMAT (2X,I6,2X,10G13.5,'*',I3)
cc 1170 FORMAT (2X,I6,2X,G13.5,' ',10F7.4,'*',I3)
 1171 FORMAT (2X,I6,2X,7G12.5)

      WRITE (FYLENO,1175) NG
 1175 FORMAT (/2X,'Implied grouping of the entities into ',
     &     I3,' groups')

      WRITE (FYLENO,1177) (IDT(I),I=1,NIND)
 1177 FORMAT (2X,10I4)

      END SUBROUTINE OUTLOOP

C     CHECK_CONSERVATION
C     ######################################################################
      SUBROUTINE CHECK_CONSERVATION(NG,NATT,TS,XMUS,XVARS,T,XMU,XVAR)
C     Input: TS = (should be 1): sample's proportion
C	     XMUS = sample's mean
C	     XVARS = sample's covariance
      IMPLICIT NONE

C     Call variables
      INTEGER, INTENT(IN) :: NG,NATT
      REAL(8), INTENT(IN) :: TS,XMUS(NATT)
     $		,XVARS(NATT,NATT)
     $		,T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT)

C     Local variables
      REAL(8) :: T_TOT,XMU_TOT(NATT),XVAR_TOT(NATT,NATT)
      INTEGER :: I,J,K,CONSERVED_FLAG(3)
      REAL(8) :: DIFF,TOL

      CONSERVED_FLAG = 1
      TOL = 1D-8
      DIFF = 0D0

C     Check conservation of mass
      DIFF = ABS( SUM(T(1:NG)) - TS )
      IF ( DIFF .LT. TOL ) THEN
        WRITE(FYLENO,*) 'TOTAL MASS CONSERVED!'
      ELSE
        WRITE(FYLENO,*) 'TOTAL MASS NOT-CONSERVED!'
        CONSERVED_FLAG(1) = 0
      ENDIF
      WRITE(FYLENO,*) '		DIFF=',DIFF

C     Check conservation of momentum
      XMU_TOT = 0D0
      DO K=1,NG
        XMU_TOT(1:NATT) = XMU_TOT(1:NATT) + T(K)*XMU(K,1:NATT)
      END DO
      DIFF = NORM2(XMU_TOT(1:NATT) - XMUS(1:NATT))
      IF ( ABS(DIFF) .LT. TOL ) THEN
        WRITE(FYLENO,*) 'TOTAL MOMENTUM CONSERVED!'
      ELSE
        WRITE(FYLENO,*) 'TOTAL MOMENTUM NOT-CONSERVED!'
        CONSERVED_FLAG(2) = 0
      ENDIF
      WRITE(FYLENO,*) '		DIFF=',DIFF

C     Check conservation of covariance
      XVAR_TOT = 0D0
      DO K=1,NG
        DO J=1,NATT      
          DO I=1,J
            XVAR_TOT(I,J) = XVAR_TOT(I,J) + XVAR(K,I,J)*T(K)
            XVAR_TOT(I,J) = XVAR_TOT(I,J) +
     $                      XMU(K,I)*XMU(K,J)*T(K)
          END DO
        END DO
      END DO      
      DO J=1,NATT
        DO I=1,J
          XVAR_TOT(I,J) = XVAR_TOT(I,J)-XMU_TOT(I)*XMU_TOT(J)
          XVAR_TOT(J,I) = XVAR_TOT(I,J)
        END DO
      END DO
      DIFF = NORM2(XVAR_TOT(1:NATT,1:NATT)-XVARS(1:NATT,1:NATT))
      IF ( ABS(DIFF) .LT. TOL ) THEN
        WRITE(FYLENO,*) 'TOTAL COVAR CONSERVED!'
      ELSE
        WRITE(FYLENO,*) 'TOTAL COVAR NOT-CONSERVED!'
        CONSERVED_FLAG(3) = 0
      ENDIF
      WRITE(FYLENO,*) '		DIFF=',DIFF

      DO J=1,3
        IF (CONSERVED_FLAG(J) .EQ. 0) STOP
      END DO

      END SUBROUTINE CHECK_CONSERVATION

      end module gm_io

c     module gm_utils
c     ##################################################################
      module gm_utils

      use math

      use mv_normal
      
      use gm_io

      integer,private :: mpierr

      contains

c     Kullback-Leibler divergence (distance)
c     ##################################################################
      SUBROUTINE KLDIV(NATT,NG,I,J,TI,TJ,XMUI,XMUJ,XVARI,XVARJ,KLdist)
c     This subroutine will compute the Kullback-Leibler divergence FROM
c     Gaussian J ~ Normal(XMUJ, XVARJ) TO Gaussian I ~ Normal(XMUI, XVARI)
c     i.e., KLdist(GI|GJ) = integral log(GI/GJ) * GI
c     VI, VJ are inverses of XVARI, XVARJ respectively
c     DVI, DVJ are determinant of XVARI, XVARJ respectively

      IMPLICIT NONE

c     Call variables      
      INTEGER, INTENT(IN)    :: NATT,NG,I,J
      REAL(8), INTENT(IN)    :: TI, TJ	! proportions of Gaussians I,J
     $	   ,XMUI(1,NATT), XMUJ(1,NATT)
     $     ,XVARI(NATT,NATT), XVARJ(NATT,NATT)
      REAL(8), INTENT(OUT)   :: KLdist

c     Local variable
      integer :: ii,jj,NULL,IERR
      ! storage for inverses & determinants of XVARI & XVARJ
      real(8) :: VI(1:NATT,1:NATT), VJ(1:NATT,1:NATT), DVI, DVJ
      real(8) :: tmp_mat(1:NATT,1:NATT)
      real(8) :: Mdist, Tval

c     Compute inverses & determinants of cov-matrices XVARJ & XVARJ
c     VI = inverse(XVARI), VJ = inverse(XVARJ)
c     DVI = det(XVARI), DVJ = det(XVARJ)
      call GDETK(NATT,I,XVARI(1:NATT,1:NATT),VI(1:NATT,1:NATT),DVI,IERR
     .          ,NULL)
      call GDETK(NATT,J,XVARJ(1:NATT,1:NATT),VJ(1:NATT,1:NATT),DVJ,IERR
     .          ,NULL)

cc     handle the situation where cov-matrix could be singular
c      if (DVI .LT. 1D-16) DVI = 1D-13
c      if (DVJ .LT. 1D-16) DVJ = 1D-13

c     compute trace(VJ*XVARI)
      tmp_mat = MATMUL(VJ(1:NATT,1:NATT),XVARI(1:NATT,1:NATT))
      Tval = 0.0D0
      do ii = 1,NATT
        Tval = Tval + tmp_mat(ii,ii)
      end do
c     compute Mahalanobis distance of XMUJ & XMUI under XVARJ
      call Mahalanobis(NATT,XMUJ,XMUI,VJ,Mdist)
      KLdist = 0.5D0*(LOG(DVJ/DVI) - dble(NATT) + Mdist + Tval)      
      KLdist = KLdist + LOG(TI/TJ)

      END SUBROUTINE KLDIV

c     EKURT
c     ##################################################################
      FUNCTION EKURT(Np,Nd) result(EK)
C     Compute expcted multivariate kurtosis of particles of one group K

      IMPLICIT NONE
      real(8) :: Np,Nd
      real(8) :: EK
      EK = Nd*(Nd+2d0)*(Np-1d0)*(Np-1d0)*(Np-1d0)/(Np*Np*(Np+1d0))

      END FUNCTION 

c     KURTOSIS
c     ##################################################################
      FUNCTION KURTOSIS(NIND,NATT,X,W,XMU,WSUM,V) result(KURT)
C     Compute multivariate kurtosis of particles of one group K

      IMPLICIT NONE

c     Call variables

      INTEGER, INTENT(IN)    :: NIND,NATT
      REAL(8), INTENT(IN)    :: X(NIND,NATT),W(NIND),WSUM
      REAL(8), INTENT(INOUT) :: XMU(NATT), V(NATT,NATT)

c     Local variables

      INTEGER :: I,J,JJ
      REAL(8) :: SUM,AL,KURT,LSUM

C     loop over all points
      SUM=0d0
      if (NIND > 0) then
         DO JJ=1,NIND
            CALL Mahalanobis(NATT,X(JJ,:),XMU,V,AL)
            SUM=SUM+AL*AL*W(JJ) 
         END DO
      endif
#if defined(petsc)
      LSUM = SUM
      call MPI_Allreduce(LSUM,SUM,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      KURT=SUM/WSUM
      
      RETURN
      END FUNCTION KURTOSIS

c     Mahalanobis
c     ##################################################################
      SUBROUTINE Mahalanobis(NATT,X,XMU,V,MD)
!     Mahalanobis distance for one point of one group
      IMPLICIT NONE

c     call variables
      INTEGER,INTENT(IN) :: NATT
      REAL(8),INTENT(IN) :: X(NATT),XMU(NATT),V(NATT,NATT)
      REAL(8),INTENT(OUT):: MD

c     Local variables
      INTEGER :: I,J

      MD = 0d0
      !!$OMP SIMD COLLAPSE(2) REDUCTION(+:MD)
      DO J = 1,NATT
         DO I = 1,NATT
            MD = MD + (X(I)-XMU(I))*V(I,J)*(X(J)-XMU(J))
         END DO
      END DO
      
      END SUBROUTINE Mahalanobis

c     splitK
c     ##################################################################
      Subroutine splitK(Np,kspl,NG,NATT,XMU,XVAR,T,DV,V)
!     Split the candidate (kspl) following the chi-squared distribution 
      implicit none
      
c     Call variables      
      INTEGER, INTENT(IN)    :: kspl,NATT
      INTEGER, INTENT(INOUT) :: NG
      REAL(8), INTENT(INOUT)    :: Np,XMU(NG+1,NATT), XVAR(NG+1,NATT
     $     ,NATT),T(NG+1),DV(NG+1),V(NG+1,NATT,NATT)

c     Local variable
      real(4) :: ran_chisq(2*NATT) !single-precision
      integer :: Nu,ISEED,i,j
      real(8) :: diag1(1:NATT),diag2(1:NATT)

!     half the prior
      T(NG+1) = 0.5d0*T(kspl)
      T(kspl) = T(NG+1)

!     using a common mean     
      XMU(NG+1,1:NATT) = XMU(kspl,1:NATT)

!     chi-squred random diagonal variance (Ververidis 2008)
      Nu = int(Np) - 1

      CALL SYSTEM_CLOCK(COUNT=iseed)
      
cc      ISEED = 1
      Call CHSRAN(2*NATT,Nu,ISEED,ran_chisq) 

      do i = 1,NATT
         diag1(i) = sqrt(ran_chisq(i     )*DV(kspl)/(2d0*NATT*Nu))
         diag2(i) = sqrt(ran_chisq(i+NATT)*DV(kspl)/(2d0*NATT*Nu))
      end do
      Do i = 1,NATT
         DO j = 1,NATT
            if(i==j) then
               XVAR(kspl,i,j) = diag1(i)
               XVAR(NG+1,i,j) = diag2(i)
            else
               XVAR(kspl,i,j) = 0d0
               XVAR(NG+1,i,j) = 0d0
            end if
         END DO
      END DO

!     find the inverse and determinant
      call MDETINV(NATT,XVAR(kspl,1:NATT,1:NATT),DV(kspl)
     $            ,V(kspl,1:NATT,1:NATT))
      call MDETINV(NATT,XVAR(NG+1,1:NATT,1:NATT),DV(NG+1)
     $            ,V(NG+1,1:NATT,1:NATT))

      NG = NG+1

      End subroutine splitK

c     merge_two_Gaussians
c     ##################################################################
      SUBROUTINE merge_two_Gaussians(NATT,NG,K,L,T,XMU,XVAR
     &				    ,V,DV)
c     This subroutine merges 2 Gaussians together:
c     Gaussian K = Gaussian K merged with Gaussian L

      IMPLICIT NONE

c     Call variables      
      INTEGER, INTENT(IN)    :: NATT,NG,K,L
      REAL(8), INTENT(INOUT) :: T(NG),XMU(NG,NATT)
     $     ,XVAR(NG,NATT,NATT),V(NG,NATT,NATT),DV(NG)

c     Local variables
      INTEGER :: IERR,NUL,I,J
      REAL(8) :: T_tmp,xmu_tmp(1:NATT)


C     Start merging: Gaussian K <- Gaussian K + Gaussian L
      T_tmp = T(K) + T(L)
      do I=1,NATT
	xmu_tmp(I) = T(K)*XMU(K,I) + T(L)*XMU(L,I)
	xmu_tmp(I) = xmu_tmp(I)/T_tmp
      end do
      do I=1,NATT
	do J=1,I
	  XVAR(K,I,J) = T(K)*(XVAR(K,I,J)+XMU(K,I)*XMU(K,J))
     .		      + T(L)*(XVAR(L,I,J)+XMU(L,I)*XMU(L,J))
          XVAR(K,I,J)=XVAR(K,I,J)/T_tmp - xmu_tmp(I)*xmu_tmp(J)
	  XVAR(K,J,I)=XVAR(K,I,J)
	end do
      end do
      T(K) = T_tmp
      XMU(K,1:NATT) = xmu_tmp(1:NATT)

C     Update det(XVAR(K,:,:)) & inv(XVAR(K,:,:)) after merging
      call GDETK(NATT,K,XVAR(K,1:NATT,1:NATT),V(K,1:NATT,1:NATT),DV(K)
     .          ,IERR,NUL)      	

      END SUBROUTINE merge_two_Gaussians

c     find_split_cand
c     ##################################################################
      Subroutine find_split_cand(NG,NATT,XMU,XVAR,T,kspl)
!     Find the split candidate by the maximum minimum upper-bound KL
!     distance between component pairs
      implicit none
      
c     Call variables
      INTEGER, INTENT(IN)    :: NG,NATT
      REAL(8), INTENT(IN)    :: XMU(NG,NATT), XVAR(NG,NATT,NATT),T(NG)
      INTEGER, INTENT(OUT)   :: kspl

c     Local variables
      integer :: i,j,K,imin,imaxmin
      real(8) :: B,Bmin,Bmaxmin,determ,LDPK,wki,wk,wi,dMU(NATT,1),P(NATT
     $     ,NATT)
c      real(8) :: MDET

c     Scan all component pairs
      select case (NG)
      case (1)
         kspl = 1
      case (2)
         kspl = 1
         if(T(2)>T(1)) kspl = 2
      case default

      Bmin = Huge(1d0)
      Bmaxmin = 0d0

      Do K = 1,NG
         determ = MDET(NATT,XVAR(K,1:NATT,1:NATT))
         LDPK = LOG(determ)
         DO i = 1,NG     
            if(i==K) CYCLE
            wki = T(K)+T(i)
            wk  = T(K)/wki
            wi  = T(i)/wki
            dMU(1:NATT,1) = XMU(K,1:NATT)-XMU(i,1:NATT)
            P = wk*XVAR(K,1:NATT,1:NATT) + wi*XVAR(i,1:NATT,1:NATT) +
     $           wk*wi*MATMUL(dMU(1:NATT,1:1),TRANSPOSE(dMU(1:NATT
     $           ,1:1)))
            determ = MDET(NATT,P)
            B = wki*Log(determ) - wk*LDPK 
            determ = MDET(NATT,XVAR(i ,1:NATT,1:NATT))
            B = B - wi*Log(determ)
            if(B<Bmin) then
               Bmin = B
               imin = i
            end if                        
         END DO
         if(Bmin>Bmaxmin) then
            Bmaxmin = Bmin
            imaxmin = imin
         end if         
      END DO
      kspl = imaxmin
      end select

      End subroutine find_split_cand

c     find_merge_cand
c     ##################################################################
      subroutine find_merge_cand(NATT,NG,K,T,XMU,XVAR,criterion,L)
c     find merged candidate to Gaussian K according to:
c	KL distance
c	smallest/largest weights
c	closest Euclidean distance of means

      implicit none

c     Call variables
      INTEGER, INTENT(IN) :: NG,NATT,K,criterion
      INTEGER, INTENT(OUT) :: L
      REAL(8), INTENT(IN) :: XMU(NG,NATT), XVAR(NG,NATT,NATT)
     $     ,T(NG)

c     Local variables
      INTEGER :: LL,I,J
      REAL(8) :: minKLdist,KLdist,minT,maxT
     $	   ,xmu_dist,mindist



      if (criterion .eq. 1) then
C       Using KL distance. I.e., find Gaussian L st it is closest to
C	Gaussian K under KL distance
	minKLdist = HUGE(1D0)
	do LL=1,NG
	  KLdist = HUGE(1D0)
	  if (LL.NE.K) then
	    ! get KL-div from Gauss K to Gauss LL
	    call KLDIV(NATT,NG,LL,K,T(LL),T(K),XMU(LL,1:NATT)
     $		,XMU(K,1:NATT),XVAR(LL,1:NATT,1:NATT),XVAR(K,1:NATT,1:NATT)
     .          ,KLdist)
	  endif
	  if (minKLdist .GT. KLdist) then
	    minKLdist = KLdist
	    L = LL
	  endif
        end do
      else if (criterion .eq. 2) then
C	Find Gaussian with smallest weight & merge to Gaussian K
	minT = HUGE(1D0)
	do LL=1,NG
	  if (T(LL) .lt. minT .and. LL .ne. K) then
	    minT = T(LL)
	    L = LL
	  endif
	end do
      else if (criterion .eq. 3) then
C	Find Gaussian with largest weight & merge to Gaussian K
	maxT = 0D0
	do LL=1,NG
	  if (T(LL) .gt. maxT .and. LL .ne. K) then
	    maxT = T(LL)
	    L = LL
	  endif
	end do
      else
C	Find Gaussian whose mean is closest & merge to Gaussian K
	mindist = HUGE(1D0)
	do LL=1,NG
	  if (LL .ne. K) then
	    xmu_dist = 0.0D0
	    do I=1,NATT
	      xmu_dist = xmu_dist + (XMU(K,I)-XMU(LL,I))**2
	    end do
	  else
	    xmu_dist = HUGE(1.0D0)
	  endif
	  if (xmu_dist .lt. mindist) then
	    mindist = xmu_dist
	    L = LL
	  endif
	end do
      endif

      end subroutine find_merge_cand

c     MDETINV
c     ########################################################################
      SUBROUTINE MDETINV (dim,A,DET,invA)

      IMPLICIT NONE

      INTEGER :: dim
      DOUBLE PRECISION, DIMENSION(dim,dim), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(dim,dim), INTENT(OUT) :: invA

      DOUBLE PRECISION :: DET


      select case (dim)
      case (1)
         DET = A(1,1)
         invA(1,1) = 1d0/A(1,1)
      case (2)
         DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
         invA(1,1) = A(2,2)/DET
         invA(2,1) =-A(2,1)/DET
         invA(1,2) =-A(1,2)/DET
         invA(2,2) = A(1,1)/DET
         
      case(3)
         DET =  A(1,1)*A(2,2)*A(3,3)  
     &        - A(1,1)*A(2,3)*A(3,2)  
     &        - A(1,2)*A(2,1)*A(3,3)  
     &        + A(1,2)*A(2,3)*A(3,1)  
     &        + A(1,3)*A(2,1)*A(3,2)  
     &        - A(1,3)*A(2,2)*A(3,1)

         INVA(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
         INVA(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/DET
         INVA(3,1) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))/DET
         INVA(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/DET
         INVA(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
         INVA(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/DET
         INVA(1,3) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))/DET
         INVA(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/DET
         INVA(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET            

      case default
         print *, "MDETINV can only do rand<=3."
         stop
      end select
      RETURN

      END SUBROUTINE MDETINV

c     MDET
c     ########################################################################
      FUNCTION MDET (dim,A) RESULT (DET)

      IMPLICIT NONE

      INTEGER :: dim
      DOUBLE PRECISION, DIMENSION(dim,dim), INTENT(IN)  :: A

      DOUBLE PRECISION :: DET

      select case (dim)
      case (1)
         DET = A(1,1)
      case (2)
         DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      case(3)
         DET =  A(1,1)*A(2,2)*A(3,3)  
     &        - A(1,1)*A(2,3)*A(3,2)  
     &        - A(1,2)*A(2,1)*A(3,3)  
     &        + A(1,2)*A(2,3)*A(3,1)  
     &        + A(1,3)*A(2,1)*A(3,2)  
     &        - A(1,3)*A(2,2)*A(3,1)
      case default
         print *, "MDET can only do rand<=3."
         stop
      end select
      RETURN

      END FUNCTION MDET


!**********************************************************************************
!  M33DET  -  Compute the determinant of a 3x3 matrix.
!**********************************************************************************

      FUNCTION M33DET (A) RESULT (DET)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A

      DOUBLE PRECISION :: DET


       DET =   A(1,1)*A(2,2)*A(3,3)  
     &       - A(1,1)*A(2,3)*A(3,2)  
     &       - A(1,2)*A(2,1)*A(3,3)  
     &       + A(1,2)*A(2,3)*A(3,1)  
     &       + A(1,3)*A(2,1)*A(3,2)  
     &       - A(1,3)*A(2,2)*A(3,1)

       RETURN

      END FUNCTION M33DET

c     GDETK
c     ########################################################################
      SUBROUTINE GDETK(NATT,K,XVAR,V,DV,IER,NULL)
C     This subroutine reads one covariance matrx, then calls
C     SYMINV, which inverts a matrix and calculates its determinant.
      IMPLICIT NONE

c     Call variables

      INTEGER, INTENT(IN)    :: NATT,K
      INTEGER, INTENT(OUT)   :: NULL,IER
      REAL(8), INTENT(IN)    :: XVAR(NATT,NATT)
      REAL(8), INTENT(INOUT) :: V(NATT,NATT),DV

c     Local variables

      INTEGER :: I,J,IT,NNULL,JJ
      REAL(8) :: TRIANG(NATT*NATT),ATRIANG(NATT*NATT)
     &          ,TEMP(NATT*NATT),CMY(NATT*NATT)
      REAL(8) :: TOL,PROD,RMAX

c     Begin program

      NNULL=0
      DO I=1,NATT*NATT
         ATRIANG(I)=0
      END DO

      IT=0
      TOL=0.0
      DO J=1,NATT
         TOL=TOL+SQRT(XVAR(J,J))
         DO I=1,J
            IT=IT+1
            TRIANG(IT)=XVAR(I,J)
         END DO
      END DO

      ! TN: Some of the cov-matrices will have very small determinant
      ! Should use small tolerance:
      TOL=(TOL/NATT)*1.0D-10
      !TOL=(TOL/NATT)*1.0D-6    ! original TOL
      CALL SYMINV (TRIANG,NATT,ATRIANG,TEMP,NULL,IER,RMAX,CMY,TOL)
      IF (IER.GT.0) GO TO 899
      IF (NULL.NE.0) THEN
        WRITE (*,815)  K, NULL
         IF (DIAG_IO) WRITE (FYLENO,815)  K, NULL
 815     FORMAT (/2X,'Rank deficiency of covariance matrix ',I3,
     &        ' is ',I3)
         NNULL=NULL+1
      ENDIF
      IT=0
      PROD=1.0
      DO J=1,NATT
         JJ=J*(J+1)/2
         PROD=PROD*CMY(JJ)*CMY(JJ)
         DO I=1,J
            IT=IT+1
            V(I,J)=ATRIANG(IT)
            V(J,I)=V(I,J)
         END DO
      END DO

      DV=PROD
      IF (DV.EQ.0) THEN
         IF (DIAG_IO) WRITE (FYLENO,*)
     .       '  Determinant is equal to zero for grp ',K
         IER=5
         GOTO 900
      ENDIF

      NULL=NNULL
      GOTO 900

 895  FORMAT (/2X,'Terminal error in SYMINV for matrix ',I3,
     &        'as IFAULT is ',I3)
 899  IF (DIAG_IO) WRITE (FYLENO,895) K, IER
      IF (STD_IO) WRITE (*,895)  K, IER

 900  RETURN
      END SUBROUTINE GDETK

c     GDET
c     ###################################################################
      SUBROUTINE GDET(NATT,NG,XVAR,V,DV,IER,NULL)
      IMPLICIT NONE
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.

C     Call variables

      INTEGER :: NATT,NG,IER,NULL
      
      REAL(8) :: XVAR(NG,NATT,NATT),V(NG,NATT,NATT),DV(NG)

c     Local variables
      
      INTEGER :: K,I

c     Begin program
      
      DO K=1,NG
        call GDETK(NATT,K,XVAR(K,1:NATT,1:NATT),V(K,1:NATT,1:NATT)
     .            ,DV(K),IER,NULL)
        if (IER == 5) then
          write (*,*) "Determinant equal to zero in group",K
          do I=1,NATT
            write (*,*) xvar(K,I,:)
          enddo
          stop
        endif
      ENDDO

      END SUBROUTINE GDET

c$$$c     GDET
c$$$c     ########################################################################
c$$$      SUBROUTINE GDET(NATT,NG,XVAR,V,DV,IER,NULL)
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$C     This subroutine reads all covariance matrices, then calls
c$$$C     SYMINV, which inverts a matrix and calculates its determinant,
c$$$C     for each covariance matrix in turn.
c$$$
c$$$      DIMENSION XVAR(NG,NATT,NATT),V(NG,NATT,NATT),
c$$$     &          DV(NG),
c$$$     &          TRIANG(NATT*NATT),ATRIANG(NATT*NATT),
c$$$     &          TEMP(NATT*NATT),CMY(NATT*NATT)
c$$$      NNULL=0
c$$$      DO 805 I=1,NATT*NATT
c$$$805      ATRIANG(I)=0
c$$$      DO 830 K=1,NG
c$$$        IT=0
c$$$        TOL=0.0
c$$$        DO 810 J=1,NATT
c$$$          TOL=TOL+SQRT(XVAR(K,J,J))
c$$$          DO 810 I=1,J
c$$$            IT=IT+1
c$$$            TRIANG(IT)=XVAR(K,I,J)
c$$$810     CONTINUE
c$$$        TOL=(TOL/NATT)*1.0d-10          ! Similar to GDETK(), we lower tolerance
c$$$        !TOL=(TOL/NATT)*0.000001        ! original TOL
c$$$        CALL SYMINV (TRIANG,NATT,ATRIANG,TEMP,NULL,IER,RMAX,CMY,TOL)
c$$$        IF (IER.GT.0) GO TO 899
c$$$        IF (NULL.NE.0) THEN
c$$$          IF (DIAG_IO) WRITE (FYLENO,815)  K, NULL
c$$$815       FORMAT (/2X,'Rank deficiency of covariance matrix ',I3,
c$$$     &                ' is ',I3)
c$$$
c$$$cc          WRITE (*,815) K, NULL
c$$$          NNULL=NULL+1
c$$$        ENDIF
c$$$        IT=0
c$$$        PROD=1.0
c$$$        DO 820 J=1,NATT
c$$$          JJ=J*(J+1)/2
c$$$          PROD=PROD*CMY(JJ)*CMY(JJ)
c$$$          DO 820 I=1,J
c$$$            IT=IT+1
c$$$            V(K,I,J)=ATRIANG(IT)
c$$$            V(K,J,I)=V(K,I,J)
c$$$820     CONTINUE
c$$$        DV(K)=PROD
c$$$        IF (DV(K).EQ.0) THEN
c$$$          IF (DIAG_IO) WRITE (FYLENO,*)
c$$$     .         '  Determinant is equal to zero for grp ',K
c$$$          IER=5
c$$$          RETURN
c$$$        ENDIF
c$$$830   CONTINUE
c$$$      NULL=NNULL
c$$$      RETURN
c$$$
c$$$ 895  FORMAT (/2X,'Terminal error in SYMINV for matrix ',I3,
c$$$     &        'as IFAULT is ',I3)
c$$$ 899  IF (DIAG_IO) WRITE (FYLENO,895) K, IER
c$$$      IF (STD_IO) WRITE (*,895)  K, IER
c$$$      RETURN
c$$$      END SUBROUTINE GDET

c     LCALK
c     ##################################################################
      SUBROUTINE LCALK(NIND,NATT,X,W,WSUM,XMU,XVAR) 
C     Compute estimate of mean, covariance matrix and mixing proportion
C     for a given Gaussian

      IMPLICIT NONE

c     Call variables

      INTEGER, INTENT(IN)    :: NIND,NATT
      REAL(8), INTENT(IN)    :: X(NIND,NATT),W(NIND),WSUM
      REAL(8), INTENT(INOUT) :: XMU(NATT), XVAR(NATT,NATT)

c     Local variables

      INTEGER :: I,J,JJ
      REAL(8) :: SUM(NATT),LSUM(NATT),LSUMX(NATT,NATT)

C     Compute new estimates of group means (XMU)

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

C     Compute new estimate of covariance matrix for each group

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

c     INVRT
c     ########################################################################
      SUBROUTINE INVRT(NATT,NCOV,NG,XVAR,V,DV,IER,NULL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAR(NG,NATT,NATT),DV(NG),V(NG,NATT,NATT)

      IF (NCOV.EQ.1) THEN
        CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
        DO 4401 K=2,NG
          DV(K)=DV(1)
          DO 4401 J=1,NATT
            DO 4401 I=1,NATT
              XVAR(K,I,J)=XVAR(1,I,J)
              V(K,I,J)=V(1,I,J)
 4401       CONTINUE

      ELSE
        CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
      ENDIF
      
      IF (IER.NE.0.AND.DIAG_IO) THEN
        WRITE (FYLENO,*) '  Problem in INVRT:'
        WRITE (FYLENO,*) '    User defined covariance matrix'
        WRITE (FYLENO,*) '    is singular.'
      ENDIF
      IF (IER.NE.0) THEN
        WRITE (*,*) '  Problem in INVRT:'
        WRITE (*,*) '    User defined covariance matrix'
        WRITE (*,*) '    is singular. Aborting...'
        stop
      ENDIF
      RETURN
      END SUBROUTINE INVRT

      end module gm_utils

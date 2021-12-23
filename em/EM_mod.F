      
c     module EM_mod
c     ##################################################################
      module EM_mod

      use em_utils

      PRIVATE ! Make everything in module private by default 
      PUBLIC :: EM

      INTEGER :: NCOV

      integer,private :: mpierr
      
      CONTAINS

c     EM
c     ##################################################################
       SUBROUTINE EM(INFO,FILEID,NIND,NATT,NG,X,WX,TOL
     $              ,TOUT,XMUOUT,XVAROUT,MITER,PLLH_OUT)
C     ------------------------------------------------------------------
C     This subroutine uses the EM algorithm from a specified starting
C     value to find a solution of the likelihood equation
C     (Component-wise [Celeux2001]), and automatically selects the
C     number of components via minimum message length criterion [Figueiredo2002].
C     On input:
C     * INFO: level of IO information (>2 -> STD output, >3 -> IO to file)
C     * FILEID: character specifying output file (if INFO>3)
C     * NIND: number of observations
C     * NATT: dimension of space of observations
C     * NG: number of initial gaussians
C     * X: array with observations
C     * WX: array with weights of observations
C     * TOL: tolerance used for convergence of PML EM iteration
C     * MITER: max # iterations used for convergence of PML EM iteration
C     On output:
C     * TOUT,XMUOUT,XVAROUT: Gaussian weights, means, variances
C     * MITER: iteration count
C     * PLLH_OUT: converged penalized likelihood value
C     ------------------------------------------------------------------
   
      IMPLICIT NONE

c    Call variables

      INTEGER, INTENT(IN) :: INFO
      INTEGER, INTENT(IN) :: NIND,NATT
      INTEGER, INTENT(INOUT) :: NG,MITER
      REAL(8), INTENT(IN) :: X(NIND,NATT),WX(NIND),TOL  ! X is the particles' vel., WX is the particle weights.
      REAL(8), INTENT(OUT):: PLLH_OUT,TOUT(NG)
     .                      ,XMUOUT(NG,NATT)
     .                      ,XVAROUT(NG,NATT,NATT)
      CHARACTER(*) :: FILEID
      
c    Local variables
      
      INTEGER :: NG_ORG ! original starting # of Gaussians
      REAL(8) :: XLOGL(MITER),MLENG(MITER)
     $     ,XCC(NG),XLA(MITER),W(NIND,NG)
     $     ,ALIM,DIFF,TEMP,CKOUT,GUM(NIND)
     $     ,WSUM(NG),XLOGL_OLD,WTOT,LWTOT,XLOGL_STAR
     $     ,AL(NIND,NG)
     $     ,RESP_SUM(NG)
      INTEGER :: I,J,K,LAST,KMAX,II,JJ,NULL,KK,IOUNT
     $     ,NPARAMS,NG_OLD,ALIVE_GRP(NG)
      
      real(8) :: CN, tmpsum,MD
      logical :: LSPLIT,AEM_CONV_FLAG
      integer :: ISPLIT,imin,lsum
      real(8) :: MML
      real(8) :: lwsum,lxlogl

      real(8) :: Ttot,XMUtot(1:NATT)
      real(8) :: min_T_old,min_T,min_T_diff,diff_perc
      real(8) :: p_ik,AL_JJ(NG) ! parameters needed only for computing XLOGL at iteration 1

      REAL(8),DIMENSION(:,:),ALLOCATABLE :: U     !(may not need)

      INTEGER,DIMENSION(:),ALLOCATABLE :: IDT     !Group index

      REAL(8),DIMENSION(:,:),ALLOCATABLE::  XMU   !Group means

      REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: V    !Inverse Variance
     &                                       ,XVAR !Group covariance matrices

      REAL(8),DIMENSION(:),ALLOCATABLE :: T !Group mixing proportions
     &                                   ,XUU     !(many not need)
     &                                   ,DV      !Determinant of V

c    Begin program

      IF (MITER > MITERMAX) THEN
        write (*,*) MITER,MITERMAX
         stop "Iteration number too large... stop."
      END IF

      IF (NG.EQ.1) THEN
         stop "At least two Gaussian mixtures are needed... stop."
      END IF

      MAXSPLIT = 0
      NCOV = 2

c$$$      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)
      
C    Set EM IO

      CALL INIT_EM_IO((info>3),(info>2),trim(FILEID))

      if (std_io) PRINT*,"Start EM: with NG=",NG

C    Setup GM EM algorithm

      ALLOCATE(U(NIND,NG),IDT(NIND),T(NG),XUU(NG),DV(NG)
     .        ,XMU(NG,NATT),V(NG,NATT,NATT),XVAR(NG,NATT,NATT))

      XUU=0

      CALL INIT_GM(NIND,NATT,NG,X,WX,IDT,XMU,XVAR,T,WTOT,DV,V)

C    Initialize output G-params to be 0D0

      TOUT(1:NG)                  = 0D0
      XMUOUT(1:NG,1:NATT)         = 0D0
      XVAROUT(1:NG,1:NATT,1:NATT) = 0D0
      
C    Parameters initialization for EM

      NG_ORG  = NG              ! original (starting) number of components
      XLOGL   = 0D0
      WTOT = SUM(WX(1:NIND))    !DBLE(NIND)
#if defined(petsc)
      LWTOT = WTOT
      call MPI_Allreduce(LWTOT,WTOT,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      ALIVE_GRP(1:NG_ORG) = 1   ! we start with all alive Gaussians; 0: killed

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C	     REINIT GAUSS. PARAMS WITH ADAPTIVE K-MEANS
C	     (comment this block out to not use kmeans)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      CN=1d0/(2d0*PI)**(NATT*0.5d0)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        SPLIT LOOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ISPLIT = 1 !the outer loop does the split operation
      MML = HUGE(1d0)
      LSPLIT = .true.

      DO WHILE (LSPLIT)
        if (DIAG_IO) write(ML_IO,*)
c         write(ML_IO,*) "#SPLIT",ISPLIT

        !initialise weights to zero
        W     = 0d0
        DIFF  = 1d0
        ALIM  = 0d0
        IOUNT = 1

        !Calculate XLOGL & MLENG at iteration 1
        DO JJ=1,NIND            !particle
          DO KK=1,NG            !group
            !Calculate component density
            CALL Mahalanobis(NATT,X(JJ,:),XMU(KK,:),V(KK,:,:),p_ik)

            !Compute particle probability per Gaussian
            p_ik=(-0.5d0)*MIN(p_ik,1d0*DENMAX)
            AL(JJ,KK)=T(KK)*EXP(p_ik)/SQRT(DV(KK))*CN
          END DO                !group
          !Calculate mixture density
          GUM(JJ) = SUM(AL(JJ,1:NG))
        END DO                  !particle

        CALL compute_pllh(NATT,T,GUM,WX,ALIVE_GRP
     .                   ,XLOGL(IOUNT),MLENG(IOUNT))
        if (DIAG_IO) write(ML_IO,*) IOUNT,XLOGL(IOUNT),MLENG(IOUNT)
     .                           ,MINVAL(T,mask=(ALIVE_GRP==1))

        min_T_old = 1D0
        
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !C       MAIN EM ITERATIVE LOOP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
        DO WHILE (IOUNT < MITER)
          !Adaptive component-wise EM
          call CW_AEM_step(NIND,NATT,NG,X,WX,W(:,1:NG),AL(:,1:NG),GUM
     .         ,DV(1:NG),T(1:NG),XMU(1:NG,:),XVAR(1:NG,:,:),V(1:NG,:,:)
     .         ,WSUM(1:NG),WTOT,ALIVE_GRP(1:NG))

          !Update iteration counter
          IOUNT = IOUNT + 1

          !Compute log-likelihood
          CALL compute_pllh(NATT,T,GUM,WX,ALIVE_GRP
     .                     ,XLOGL(IOUNT),MLENG(IOUNT))
          if (DIAG_IO) write(ML_IO,*) IOUNT,XLOGL(IOUNT),MLENG(IOUNT)
     .                               ,MINVAL(T,mask=(ALIVE_GRP==1))

          !Standard convergence criteria      
          IF (IOUNT .GT. 1) THEN
            LAST = IOUNT-1
            ALIM=TOL*XLOGL(LAST)
            DIFF=XLOGL(IOUNT)-XLOGL(LAST)         
          END IF

          if(ABS(DIFF).le.ABS(ALIM)) EXIT               

        END DO !END LOOP OF WHILE (IOUNT<MITER ...)
          
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !C      END MAIN EM ITERATIVE LOOP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
C      Here, we have NG Gaussian available, others are gone!

        T   (NG+1:NG_ORG)               = 0D0
        XMU (NG+1:NG_ORG,1:NATT)        = 0D0
        XVAR(NG+1:NG_ORG,1:NATT,1:NATT) = 0D0

C      Perform one 'final' standard (non-adaptive) EM step to 
C      enforce conservation

        call std_EM_step(NIND,NATT,NG,X,WX,W(:,1:NG),GUM,DV(1:NG)
     .      ,T(1:NG),XMU(1:NG,:),XVAR(1:NG,:,:),V(1:NG,:,:),WSUM(1:NG)
     .      ,WTOT)

        if (STD_IO) print *,"End EM-loop with: ",NG," components"

        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !C     END OF EM LOOP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        LSPLIT = split_gaussian(isplit,TOL,NIND,NATT,NG,X,WX,W(:,1:NG)
     .       ,GUM,DV(1:NG),T(1:NG),XMU(1:NG,:),XVAR(1:NG,:,:)
     .       ,V(1:NG,:,:),WSUM(1:NG),XLOGL(IOUNT-1),MLENG(IOUNT-1))
       
      END DO !SPLIT-LOOP

c     Final work after convergence
      
c$$$      IF (IER.EQ.-111) THEN
c$$$        WRITE(*,*) 'WARNING : Some points have zero Likelihood'
c$$$        WRITE(*,*) '         (will denote with 0 in grouping)'
c$$$        IER=0
c$$$      ENDIF

      CALL CAPART(NIND,NATT,NG,W,IDT,XCC)

C     Enforce conservation for moments

c$$$      CALL ENFORCE_MOMS_CONSERVATION(NG,NATT,TS,XMUS,XVARS
c$$$     $	,T(1:NG),XMU(1:NG,:),XVAR(1:NG,:,:))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     IO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (STD_IO.and.IOUNT>2) then
        write (*,*) "Final # groups:", NG
        write (*,*) "Final iteration:",IOUNT-1
        write (*,*) "Relative Tolerance:",TOL
        write (*,*) "Initial LOG-LIKELIHOOD:",XLOGL(1)
        write (*,*) "Last two LOG-LIKELIHOOD:",XLOGL(IOUNT-1:IOUNT)
        write (*,*) "Initial Message-length:",MLENG(1)
        write (*,*) "Last two Message-length:",MLENG(IOUNT-1:IOUNT)

        write (*,*)
        write (*,*) "Mixture groups:"
        write (*,*)
        
        DO K=1,NG
          WRITE (*,*) "Weight=",T(K)
          WRITE (*,*) "Mean  =",(XMU(K,J),J=1,NATT)
          DO J=1,NATT
            WRITE (*,*) "Variance=",(XVAR(K,I,J),I=1,J)
          ENDDO
          WRITE (*,*) "Det chk=",trace(NATT,XVAR(K,:,:))
     $                         /(NATT*(det(XVAR(K,:,:)))**(1d0/NATT))
          WRITE (*,*)
        ENDDO
      ENDIF

      IF (DIAG_IO) THEN
        CALL OUTLOOP(NIND,NATT,NG,XMU(1:NG,1:NATT),DV(1:NG),T(1:NG)
     &              ,NCOV,IOUNT,XLOGL,MLENG,W(1:NIND,1:NG),IDT,X)

        CALL FINOUT(NIND,NATT,NCOV,NG,XUU(1:NG),XMU(1:NG,1:NATT)
     &             ,XVAR(1:NG,1:NATT,1:NATT),T(1:NG),IDT,XCC(1:NG))
      ENDIF

C    Transfer Gaussian weights, mu, and sigma and PLLH to subroutine outputs      

      TOUT   (1:NG)     = T   (1:NG)
      XMUOUT (1:NG,:)   = XMU (1:NG,:)
      XVAROUT(1:NG,:,:) = XVAR(1:NG,:,:)

      PLLH_OUT = XLOGL(IOUNT-1)

      MITER = IOUNT
      
C     Destroy temporary storage      

      DEALLOCATE(U,IDT,XMU,V,XVAR,T,XUU,DV)

      CALL DESTROY_GM()

      if (diag_io) Call CLOSE_EM_FILES()

      END SUBROUTINE EM

c$$$C######################################################################
c$$$	subroutine check_correct_spd(A,n,tol)
c$$$	implicit none
c$$$C	Call variables
c$$$	integer, intent(in) :: n
c$$$	real(8), intent(in) :: tol
c$$$	real(8), intent(inout) :: A(n,n)
c$$$C	Local variables
c$$$	real(8) :: w(n),evec(n,n),evals(n,n)
c$$$	real(8) :: tempmat(n,n)
c$$$	integer :: i,j,k
c$$$
c$$$	evals = A
c$$$	evec = 0D0
c$$$	call EIG_Jacobi(evals,evec,1D-5,n)
c$$$	do i=1,n
c$$$	  w(i) = evals(i,i)
c$$$	end do
c$$$	if (MINVAL(w) .lt. 0D0) then
c$$$	  write(FYLENO,*) 'XVAR(K,:,:) NOT SPD'
c$$$	  do i=1,n
c$$$	    if (w(i) .lt. tol) w(i) = tol
c$$$	    tempmat(:,i) = w(i)*evec(:,i)
c$$$	  end do
c$$$	  A = 0.0D0
c$$$	  call dgemm('N','T',n,n,n,1D0,tempmat,n,evec,n,0D0,A,n)
c$$$	  A = 0.5D0*(A+transpose(A))
c$$$	endif
c$$$
c$$$	end subroutine check_correct_spd
c$$$

      end module EM_mod

c     module AMEM_mod
c     ##################################################################
      module AMEM_mod

      use em_utils
            
      USE AAMEM_UTILITIES_MOD
      USE REG_AA_MOD

      PRIVATE ! Make everything in module private by default 
      PUBLIC :: AMEM

      INTEGER :: NCOV
CTOP_OMP!$OMP THREADPRIVATE(NCOV)

      integer,private :: mpierr

      CONTAINS

c     AMEM
c     ##################################################################
       SUBROUTINE AMEM(INFO,FILEID,NIND,NATT,NG,X,WX,TOL
     $                ,TOUT,XMUOUT,XVAROUT,MITER,PLLH_OUT)
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
C     * MITER: # iterations used for convergence of PML EM iteration
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
      REAL(8), INTENT(IN) :: X(NIND,NATT)    ! X is the particles' vel.
     .                      ,WX(NIND)        ! WX is the particle weight
     .                      ,TOL
      REAL(8), INTENT(OUT):: PLLH_OUT,TOUT(NG)
     .                      ,XMUOUT(NG,NATT)
     .                      ,XVAROUT(NG,NATT,NATT)
      CHARACTER(*) :: FILEID
      
c    Local variables

      INTEGER :: G_DOF,NPARBY2,IOUNT
      INTEGER :: NG_ORG ! original starting # of Gaussians
      REAL(8) :: XLOGL(MITER),MLENG(MITER)
     $     ,XCC(NG),W(NIND,NG)
     $     ,ALIM,DIFF,GUM(NIND)
     $     ,WSUM(NG),WTOT,LWTOT,XLOGL_STAR
     $     ,AL(NIND,NG),RESP_SUM(NG)
      INTEGER :: I,J,K,LAST,II,JJ,NULL,KK,IERR
     $     ,NPARAMS,ALIVE_GRP_OLD(NG),NG_OLD,ALIVE_GRP(NG)

      real(8) :: CN,MML
      logical :: LSPLIT
      integer :: ISPLIT

      REAL(8),DIMENSION(:,:),ALLOCATABLE :: U      !(may not need)

      INTEGER,DIMENSION(:),ALLOCATABLE :: IDT      !Group index

      REAL(8),DIMENSION(:,:),ALLOCATABLE::  XMU    !Group means

      REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: V    !Inverse Variance
     &                                       ,XVAR !Group covariance matrices

      REAL(8),DIMENSION(:),ALLOCATABLE :: T        !Group mixing proportions
     &                                   ,XUU      !(many not need)
     &                                   ,DV       !Determinant of V

c     Variables for A-AMEM

      INTEGER :: ITER
      TYPE(AA_STATE) :: AASTATE
      REAL(8) :: S_NEW,RAA_VARS(5),SK0

c     Variables for monotonicity control using the gradient of LLH w.r.t parameters
      LOGICAL :: AA_RELAX_FLAG,AA_MON_FLAG

      REAL(8),ALLOCATABLE,DIMENSION(:) :: THETA,THETA_OLD,GVAL,DTHETA
      REAL(8) :: time_start,time_end

!c     Variable for debugging
!      CHARACTER(300) :: debug_io

c      Variables for plotting output (xfer-amem-paper)
       INTEGER :: kill_gauss_hist(NG), aa_mono_hist(MITER)
       LOGICAL :: amem_diag

      
c    Begin program


      amem_diag = DIAG_AMEM     ! make a copy of DIAG_AMEM (from input deck)

      IF (MITER > MITERMAX) THEN
         stop "Iteration number too large... stop."
      END IF

      IF (NG.EQ.1) THEN
         stop "At least two Gaussian mixtures are needed... stop."
      END IF

      XLOGL   = 0D0
      MLENG   = 0D0
      XCC     = 0D0
      RESP_SUM= 0D0
      GUM     = 0D0
      WSUM    = 0D0
      WTOT    = 0D0
      AL      = 0D0
      XLOGL   = 0d0

      S_NEW   = 0D0

      NULL    = 0
      NG_ORG  = NG              ! original (starting) number of components

      MAXSPLIT = 0
      NCOV = 2

      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)

      ! G_DOF = degree of freedom of 1 anisotropic Gaussian:
      ! - 1 for weight
      ! - NVDIM for mean
      ! - NVDIM*(NVDIM+1)/2 for covariance
      ! Useful for creating AA-vector consisting of all Gaussians' parameters
      G_DOF = 1+NATT+NATT*(NATT+1)/2
      NPARAMS = G_DOF*NG_ORG    ! starting total number of Gaussian params

      ALLOCATE(THETA(NPARAMS),THETA_OLD(NPARAMS)
     .        ,GVAL(NPARAMS),DTHETA(NPARAMS))

      THETA = 0d0 ; THETA_OLD = 0d0 ; GVAL = 0d0 ; DTHETA = 0d0
      
C    Set CWEM IO

      CALL INIT_EM_IO((info>3),(info>2),trim(FILEID))

      if (std_io) PRINT*,"Start AMEM: with NG=",NG

C    Setup GM EM algorithm

      ALLOCATE(U(NIND,NG),IDT(NIND),T(NG),XUU(NG),DV(NG)
     .        ,XMU(NG,NATT),V(NG,NATT,NATT),XVAR(NG,NATT,NATT))

      U = 0 ; IDT = 0 ; XUU=0 ; T = 0 ; XUU = 0
      DV = 0 ; XMU = 0 ; V = 0 ; XVAR = 0

      CALL INIT_GM(NIND,NATT,NG,X,WX,IDT,XMU,XVAR,T,WTOT,DV,V)  !This can change NG

C     Initialize output G-params to be 0D0

      TOUT    = 0D0
      XMUOUT  = 0D0
      XVAROUT = 0D0


      ! timing the AMEM routine (without/with AA)
      call cpu_time(time_start)
      
C    Parameters initialization for CWEM

      WTOT = SUM(WX)
#if defined(petsc)
      LWTOT = WTOT
      call MPI_Allreduce(LWTOT,WTOT,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      ALIVE_GRP = 0
      ALIVE_GRP(1:NG) = 1   ! we start with all alive Gaussians; 0: killed
      ALIVE_GRP_OLD = ALIVE_GRP

      CN=(2d0*PI)**(-NATT*0.5d0)  !Notice negative sign in exponent

      CALL FORM_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $       ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),THETA)

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
        ITER = 0

        !Calculate XLOGL & MLENG at iteration 1
        CALL COMPUTE_PLLH_VER2(ITER,NIND,G_DOF,NG_ORG,NG,NATT
     $         ,WX,X,THETA,XLOGL(IOUNT),MLENG(IOUNT),AL,GUM,ALIVE_GRP)
        
        XLOGL_STAR = XLOGL(IOUNT)
        DIFF = 1.0D0

        !Initialize (regularized) Anderson acceleration state
        RAA_VARS = (/ LAMDA,RVAL,ALPHA,KAPPA,SK /) ! must be in this order
        CALL AA_CREATE(AASTATE,MVEC_AA,AASTART,NPARAMS,RAA_VARS) ! create AA state
        SK0 = SK

        kill_gauss_hist = 0
        aa_mono_hist = 0
        II=1
        
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !C       MAIN EM ITERATIVE LOOP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
        !--------------------------------------------------------------------------------
        ! Import from subroutine ADAPTIVE_DAAREM (in GMMEM_MOD of A-AMEM stand-alone code)
        !--------------------------------------------------------------------------------
        ! Main iterative loop
        DO WHILE (ITER < MITER-1)
          AA_RELAX_FLAG = .FALSE.
          AA_MON_FLAG   = .TRUE.

          !----------------------------------------------------------------------------
          ! STEP (I): PICARD FOR A-EM ALGORITHM & CONVERGENCE CHECK
          !----------------------------------------------------------------------------
          ! GVAL = G(THETA), GUM will be needed for computing loglikelihood at current iteration
          ALIVE_GRP_OLD = ALIVE_GRP
          NG_OLD = NG
          THETA_OLD = THETA    ! Gaussian params at 'current' iteration, ITER
          ! Perform one A-EM step: GVAL = EM(THETA)
          CALL AEM_STEP(G_DOF,NG_ORG,NG,NIND,NATT,WTOT
     $         ,WX,X,THETA,GVAL,AL,GUM,RESP_SUM,ALIVE_GRP)

          !Check convergence of nonlinear residual
          if (sqrt(sum((theta-gval)*(theta-gval))/NPARAMS) < 1d-13) exit
          
          ! Restart AA if at least one group is killed
          IF (NG .LT. NG_OLD) THEN
            CALL RESTART_AA(AASTATE)

            if (STD_IO) write (*,*) "Gaussian killed...restarting AA"

            if (amem_diag) then 
              kill_gauss_hist(II)=ITER
              II=II+1
            end if

            ! Sanity/safety check: if after 1 Picard update, all Gaussian(s) are killed, 
            ! then we use the previous solution & terminate the iteration!
            IF (NG .EQ. 0) THEN
              ALIVE_GRP = ALIVE_GRP_OLD
              NG = NG_OLD
              THETA = THETA_OLD
              EXIT
            END IF
          END IF
          ! residual (shall be used in monotonicity check)
          DTHETA = GVAL - THETA

          ! Update iteration counter
          IOUNT = IOUNT + 1

          !------------------------------------------------------------------------------
          ! STEP (II): ANDERSON ACCELERATION W/ OR W/OUT (TIKKONOV)
          ! REGULARIZATION TERM
          !------------------------------------------------------------------------------
          ! IF MVEC or MVEC_AA == 0, then APPLY_AA simply choose THETA =
          ! GVAL and move on to the next iteration. No AA's least-sqrt
          ! problem is solved.
          CALL APPLY_AA(AA_RELAX_FLAG)


          ! Handle the case when ITER = 0 or AASTATE%MAA = 0.
          !   In these cases, haven't done AA yet, simply compute PLLH(G(THETA)) & move on
c$$$          IF ( (AASTATE%MVEC .EQ. 0) .OR. (ITER .LT. AASTATE%AASTART)
c$$$     $                               .OR. (AASTATE%MAA .EQ. 0) ) THEN
          IF (    (AASTATE%MVEC == 0)
     .        .OR.(ITER < AASTATE%AASTART)
     .        .OR.(AASTATE%MAA == 0)
     .       ) THEN
            CALL COMPUTE_PLLH_VER2(IOUNT,NIND,G_DOF,NG_ORG,NG,NATT
     $               ,WX,X,THETA,XLOGL(IOUNT),MLENG(IOUNT),AL,GUM
     $               ,ALIVE_GRP)
            SK0 = S_NEW
          ELSE
cc            !------------------------------------------------------------------------------
cc            ! STEP (II): ANDERSON ACCELERATION W/ OR W/OUT (TIKKONOV) REGULARIZATION TERM
cc            !------------------------------------------------------------------------------
cc            !CALL APPLY_AA(AA_RELAX_FLAG)


            !------------------------------------------------------------------------------
            ! STEP (III): MONOTONICITY CONTROL STEP FOR PENALIZED LOG-LIKELIHOOD
            ! For this step, we follow exactly Ref. [2] but using score function
            ! (PLLH's 1st deriv) checking the penalized log-likelihood of AA & EM iterates
            !------------------------------------------------------------------------------
            CALL PLLH_MONOTONICITY_CTRL(AA_RELAX_FLAG,AA_MON_FLAG)
            IF(AA_MON_FLAG) aa_mono_hist(IOUNT)=1

            !----------------------------------------------------------------------
            ! STEP (IV): AA RESTART?
            !----------------------------------------------------------------------
            CALL AA_RESTART(AA_RELAX_FLAG,AA_MON_FLAG)
            
          END IF                ! end-if (...) OR (...) OR (...)

          ! Update value of SK in AA state variable (bug-fixed 01/21/2021)
          AASTATE%RAA_VARS%SK = SK0
          
          !Check convergence
          IF (IOUNT > 1) THEN
            LAST = IOUNT-1
            ALIM = TOL*XLOGL(LAST)
            DIFF = XLOGL(IOUNT)-XLOGL(LAST)
          END IF

          if(ABS(DIFF) <= ABS(ALIM)) EXIT
          
          ITER = ITER + 1

        END DO                  ! enddo while (ITER .LT. MITER)
!--------------------------------------------------------------------------------

        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !C      END MAIN EM ITERATIVE LOOP
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!C      If AA is turned on, the converged solution is in vector THETA(:)
!C      we need to extract it. To do so, use the ALIVE_GRP vector to extract
!C      the remaining alive Gaussians

        CALL AA_DESTROY(AASTATE)
          
        ! Extract converged params
        CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $       ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),THETA)
        ! Here at this point, in T,XMU,XVAR, there are dead Gaussians,
        ! i.e., T(K) = 0 for some K=1,...,NG_ORG
        ! and ALIVE_GRP(K) = 0 for some K=1,...,NG_ORG
        ! Extract alive Gaussians: we keep alive Gaussians in Gaussians 1,...,NG
        KK=0
        DO K=1,NG_ORG
          IF (ALIVE_GRP(K) .EQ. 1) THEN
            KK=KK+1
            TOUT(KK)                 =T(K)
            XMUOUT(KK,1:NATT)        =XMU(K,1:NATT)
            XVAROUT(KK,1:NATT,1:NATT)=XVAR(K,1:NATT,1:NATT)
          END IF
        END DO
        ! Sanity/safety check: (hope we don't run into this problem)
        IF (KK .NE. NG) THEN
          PRINT*,"No. of alive Gaussians mismatch!STOP!"
          STOP
        END IF
        T(1:NG)                 =TOUT(1:NG)
        XMU(1:NG,1:NATT)        =XMUOUT(1:NG,1:NATT)
        XVAR(1:NG,1:NATT,1:NATT)=XVAROUT(1:NG,1:NATT,1:NATT)
        ! Update inverses & determinants of XVAR of alive Gaussians
        DO K=1,NG
          CALL GDETK(NATT,K,XVAR(K,1:NATT,1:NATT),V(K,1:NATT,1:NATT)
     .              ,DV(K),IERR,NULL)
        END DO
        ALIVE_GRP(1:NG) = 1
        ALIVE_GRP(NG+1:NG_ORG) = 0
        
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


      ! Finish timing AMEM routine (without/with AA)
      call cpu_time(time_end)
      if (amem_diag) print*,'amem runtime =',time_end-time_start

c     Final work after convergence
      
c$$$      IF (IERR.EQ.-111) THEN
c$$$        WRITE(*,*) 'WARNING : Some points have zero Likelihood'
c$$$        WRITE(*,*) '         (will denote with 0 in grouping)'
c$$$        IERR=0
c$$$      ENDIF

      CALL CAPART(NIND,NATT,NG,W,IDT,XCC)

cC     Enforce conservation for moments
c      CALL ENFORCE_MOMS_CONSERVATION(NG,NATT,TS,XMUS,XVARS
c     $	,T(1:NG),XMU(1:NG,:),XVAR(1:NG,:,:))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     IO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (STD_IO.and.IOUNT>2) then
        print *,"Final # groups:", NG
        print *,"Final iteration:",IOUNT-1
        print *,"Relative Tolerance:",TOL
        print *,"Initial LOG-LIKELIHOOD:",XLOGL(1)
        print *,"Last two LOG-LIKELIHOOD:",XLOGL(IOUNT-1)
     .         ,XLOGL(IOUNT-2)
        print *,"Initial Message-length:",MLENG(1)
        print *,"Last two Message-length:",MLENG(IOUNT-1)
     .         ,MLENG(IOUNT-2)
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
      DEALLOCATE(THETA,THETA_OLD,GVAL,DTHETA)
      
      CALL DESTROY_GM()

      if (diag_io) Call CLOSE_EM_FILES()



C     Extract the pllh history for plotting (xfer-aamem-paper)
      if (amem_diag) then
      open(1,file='pllh_amem.dat',action='write')
      DO I=1,IOUNT-1
        WRITE(1,*) I,XLOGL(I)
      END DO
      close(1)

      open(1,file='kill_gauss_hist_amem.dat',action='write')
      DO I=1,NG_ORG
        IF(kill_gauss_hist(I).GT.0) WRITE(1,*) kill_gauss_hist(I)
      END DO
      close(1)

      open(1,file='aa_mono_hist.dat',action='write')
      DO I=1,IOUNT
        WRITE(1,*) I,aa_mono_hist(I)
      END DO
      WRITE(1,*) I+1,0  ! the in the very last iteration, we do standard EM to conserve moments
      end if

      CONTAINS

!c     APPLY_AA
!c     ###############################################################
      SUBROUTINE APPLY_AA(AA_RELAX)

      IMPLICIT NONE

      LOGICAL :: AA_RELAX
      
      !Local variables
      REAL(8) :: GVAL_L(AASTATE%N),THETA_L(AASTATE%N)
      
      !Begin program

      GVAL_L = 0d0 ; THETA_L = 0d0
      
      ! Convert cov-matrices in THETA,GVAL into lower-triangular
      ! matrices L in THETA_L,GVALL
      CALL XVAR_CHOLESKY_TRANSFORM(G_DOF,NG_ORG,NATT
     $        ,THETA,GVAL,THETA_L,GVAL_L,ALIVE_GRP)
      ! GVALL & THETA_L stores the lower-triangular matrices instead of cov-matrices for
      ! GVAL & THETA, respectively. This is to avoid the violation of positivity of cov-matrices in AA
      !FVAL = GVALL - THETA_L
      IF (REGULARIZED_FLAG) THEN
        !On output THETA or THETA_L is accelerated solution
!        CALL REGULARIZED_AA(AASTATE,ITER,THETA,GVAL)
        CALL REGULARIZED_AA(AASTATE,ITER,THETA_L,GVAL_L)
      ELSE
!        CALL STANDARD_AA(AASTATE,ITER,THETA,GVAL)
        CALL STANDARD_AA(AASTATE,ITER,THETA_L,GVAL_L)
      END IF


      ! Here, we re-obtain the cov-matrices from lower-triangular matrices:
      !   In particular, THETA now holds Gaussian weights, means & sigmas
      !   where THETA_L holds weights, means, and lower-triangular matrices.
      CALL CHOLESKY_XVAR_TRANSFORM(G_DOF,NG_ORG,NATT
     $                    ,THETA_L,THETA,ALIVE_GRP)

      ! Relax AA solution: if a weight is negative, roll back to EM solution
      !   and move on to the next iteration (without checking monotonicity)
      ! Output AA_RELAX_FLAG of RELAX_ACCEL_SOLN() will help us moving on to the next
      !   iteration
      CALL RELAX_ACCEL_SOLN(G_DOF,NG_ORG,NATT,ITER
     $            ,THETA,GVAL,AA_RELAX,ALIVE_GRP)
      ! If we restart AA when a group is killed, the following call is not necessary.
      !   However, still keep it for sanity check - cheap task!
      CALL ENFORCE_DEAD_GAUSSIAN(G_DOF,NG_ORG,NATT,THETA,ALIVE_GRP)
      ! Normalize Gaussian weights (???)
      THETA(1:NG_ORG) = THETA(1:NG_ORG)/SUM(THETA(1:NG_ORG))

      END SUBROUTINE APPLY_AA

!c     PLLH_MONOTONICITY_CTRL
!c     ###############################################################
      SUBROUTINE PLLH_MONOTONICITY_CTRL(AA_RELAX,MON_CHK)

      IMPLICIT NONE

      LOGICAL :: AA_RELAX,MON_CHK
      
      ! If AA solution is fine (AA_RELAX_FLAG = .FALSE.), we compare LLH(AA) vs LLH(EM)
      ! using Taylor expansion idea (see subroutine CHECK_MONOTONICITY). This avoid the
      ! re-evaluation of LLH(AA) which is very expensive in EM-GMM context.
      IF (.NOT.AA_RELAX) THEN
        ! Smart LLH comparison using Taylor expansion
        MON_CHK = CHECK_MONOTONICITY(G_DOF,NIND,NG_ORG
     .       ,NG,NATT,THETA,THETA_OLD,EPS_VAL,RESP_SUM,DTHETA
     .       ,ALIVE_GRP)

 
        ! Naive LLH comparison -- keep this code (dont delete), useful for efficiency check.
c        CALL COMPUTE_PLLH_VER2(ITER,NIND,G_DOF,NG_ORG,NG,NATT,WX
c     $          ,X,THETA,XLOGL(IOUNT),MLENG(IOUNT),AL,GUM,ALIVE_GRP)
c        IF (XLOGL(IOUNT)-XLOGL(IOUNT-1) .GT. -1d0*EPS_VAL) THEN
c          MON_CHK = .TRUE.
c        ELSE
c          MON_CHK = .FALSE.
c        END IF

        if (STD_IO) write (*,*) "PLLH monotonicity chk passes?",mon_chk

        ! If LLH-monotonicity chk passes, compute PLLH(AA-soln) and move on
        IF (MON_CHK) THEN
          CALL COMPUTE_PLLH_VER2(ITER,NIND,G_DOF,NG_ORG,NG,NATT,WX
     $          ,X,THETA,XLOGL(IOUNT),MLENG(IOUNT),AL,GUM,ALIVE_GRP)
          S_NEW = SK0 + 1D0
        ELSE
        ! If LLH-monotonicity chk does not pass, revert back to EM-soln & compute PLLH(EM-soln)
          THETA = GVAL
          S_NEW = SK0
          CALL COMPUTE_PLLH_VER2(ITER,NIND,G_DOF,NG_ORG,NG,NATT,WX
     $          ,X,THETA,XLOGL(IOUNT),MLENG(IOUNT),AL,GUM,ALIVE_GRP)
        END IF
      ! Else-if (AA_RELAX) then 
      !   we just evaluate the LLH of EM solution & move on
      ! No monotonicity check is required.
      ELSE
        MON_CHK = .TRUE.
        CALL COMPUTE_PLLH_VER2(ITER,NIND,G_DOF,NG_ORG,NG,NATT,WX
     $        , X,GVAL,XLOGL(IOUNT),MLENG(IOUNT),AL,GUM,ALIVE_GRP)
        S_NEW = SK0 + 1D0
      END IF              ! end-if (AA_RELAX .EQV. .FALSE.)
      SK0 = S_NEW
      END SUBROUTINE PLLH_MONOTONICITY_CTRL

!c     AA_RESTART
!c     ###############################################################
      SUBROUTINE AA_RESTART(AA_RELAX,MON_CHK)

      IMPLICIT NONE

      LOGICAL :: AA_RELAX,MON_CHK

      LOGICAL :: AA_RSTRT

c      AA_RSTRT=
c     .     (MOD(AASTATE%MAA,MVEC_AA).EQ.0).OR.AA_RELAX.OR.(.NOT.MON_CHK)
      
      ! TN: I am not sure if we need to restart AA when we relax AA
      ! solution back to EM solution or if MON-CHCK not passed...
      ! We can visit this condition later.

      AA_RSTRT=(MOD(AASTATE%MAA,MVEC_AA).EQ.0)

      !IF (MOD(ITER,MVEC_AA) .EQ. 0) THEN
      IF (AA_RSTRT) THEN
        if (STD_IO)  write (*,*) "Periodic restart of AA..."

        IF (XLOGL(IOUNT) .GE. XLOGL_STAR - EPS_C) THEN
          CALL RESTART_AA(AASTATE)
          XLOGL_STAR = XLOGL(IOUNT)
        ELSE
          S_NEW = MAX((S_NEW-MVEC_AA)*1D0,-1.0D0*DVAL)
          CALL RESTART_AA(AASTATE)
          XLOGL_STAR = XLOGL(IOUNT)
        END IF
      END IF  ! end-if (MOD(AASTATE%MAA,MVEC_AA) .EQ. 0)
      
      SK0 = S_NEW
      
      END SUBROUTINE AA_RESTART

      END SUBROUTINE AMEM
      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   Relevant subroutines for A-AMEM algorithm will be added here
!   In particular, we need to include the following (so far, could be more)
!   - 'simultaneous' adaptive_EM_step
!   - compute_pllh (using GUM which is output from adaptive_EM_step)
!   - check_monotonicity subroutine
!   - helper subroutines such as cholesky_xvar, enforce_dead_gaussian, relax_aa_solution, etc...
!   - anything else??
!   Note: Subroutines for doing AA will be employed from REG_AA_MOD (driver for regularized AA)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!c     AEM_STEP
!c     ###############################################################
      SUBROUTINE AEM_STEP(G_DOF,NG_ORG,NG,NIND,NATT,WTOT
     $            ,WX,X,THETA,GVAL,AL,GUM,WSUM,ALIVE_GRP)
        ! Unlike adaptive_EM_step() which utilizes component-wise EM (Celeux2001),
        ! This subroutine performs 1 adaptive-em step in simultaneous manner
        !   (given a good guess for starting number of component)
        ! This subroutine follows the implementation of subroutine ADAPTIVE_EM_STEP_OPT
        ! in GMMEM_MOD (from A-AMEM's stand-alone code)
        ! Small update: considering vector of particles' weights, WX(1:NIND), in the E-step

      IMPLICIT NONE

      ! Call variables
      INTEGER, INTENT(IN) :: G_DOF,NG_ORG,NIND,NATT
      INTEGER, INTENT(INOUT) :: NG    ! This variable keeps track of current number of alive Gaussian
      INTEGER, INTENT(INOUT) :: ALIVE_GRP(NG_ORG)
      REAL(8), INTENT(IN) :: X(NIND,NATT),WX(NIND),WTOT  ! X = particle velocities, WX = particle weights
      REAL(8), INTENT(IN) :: THETA(G_DOF*NG_ORG)
     $                     ,GUM(NIND),AL(NIND,NG_ORG)
      REAL(8), INTENT(OUT) :: GVAL(G_DOF*NG_ORG),WSUM(NG_ORG)

      ! Local variables
      REAL(8) :: T(NG_ORG),XMU(NG_ORG,NATT),XVAR(NG_ORG,NATT,NATT)
     $          ,V(NG_ORG,NATT,NATT),DV(NG_ORG),SQRT_DV(NG_ORG)
     $          ,W(NIND,NG_ORG)

      ! Begin program

      ! First we unpack the solution vector U into T(1:NG_ORG), XMU(1:NG_ORG,:)
      !   XVAR(1:NG_ORG,:,:)
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $      ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),THETA)

      ! Perform AEM step
c      call SIM_AEM_step(NIND,NATT,NG,X,WX,W(:,1:NG),AL(:,1:NG),GUM
c     .         ,DV(1:NG),T(1:NG),XMU(1:NG,:),XVAR(1:NG,:,:),V(1:NG,:,:)
c     .         ,WSUM(1:NG),WTOT,ALIVE_GRP(1:NG))

      call SIM_AEM_step_2(NIND,NATT,NG_ORG,NG,X,WX,W(:,1:NG_ORG)
     .         ,AL(:,1:NG_ORG),GUM,DV(1:NG_ORG),T(1:NG_ORG)
     .         ,XMU(1:NG_ORG,:),XVAR(1:NG_ORG,:,:),V(1:NG_ORG,:,:)
     .         ,WSUM(1:NG_ORG),WTOT,ALIVE_GRP(1:NG_ORG))

      ! Finally, use the updated params to form the output soln vector: GVAL = Picard_EM(THETA)
      CALL FORM_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $      ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),GVAL)
      !--------------------------------------------------------------------------------------   
      END SUBROUTINE AEM_STEP

!c     CHECK_MONOTONICITY
!c     ###############################################################
!c     This subroutine uses score function (1st deriv. of LLH function)
!c     to do the checking monotonoticity of AAEM solutions 
!c     I.e., we are using Taylor's expansion idea to estimate LHS of
!c     LLH(THETA_AA) - LLH(THETA_{it}) > - epsilon_val = -0.01
!c     <=> D LLH(THETA_{it}) / D THETA \DOT (THETA_AA - THETA_{it}) > - 0.01 ??
!######################################################################
      FUNCTION CHECK_MONOTONICITY(G_DOF,NIND,NG_ORG,NG,NATT,THETA_AA
     $                           ,THETA,EPS_VAL,WSUM,FVAL,ALIVE_GRP)
     .         RESULT(MONOTONICITY)

      IMPLICIT NONE
    
      !   Call variables
      INTEGER, INTENT(IN) :: G_DOF,NIND,NG_ORG,NG,NATT
      INTEGER, INTENT(IN) :: ALIVE_GRP(NG_ORG)
      REAL(8), INTENT(IN) :: WSUM(NG_ORG)
      REAL(8), INTENT(IN) :: THETA_AA(G_DOF*NG_ORG),THETA(G_DOF*NG_ORG)
     $                     ,EPS_VAL
      REAL(8), INTENT(IN) :: FVAL(G_DOF*NG_ORG)   ! residual of AA: GVAL - U = G(U^{it}) - U^{it}
      LOGICAL :: MONOTONICITY        ! if TRUE, use AA solution
                                     ! if FALSE, fall back to PICARD solution
      !   Local variables
      INTEGER :: I,J,K,JJ,IER,NULL,NINDT
      REAL(8) :: GRAD_LLH(G_DOF*NG_ORG),NPARBY2
      REAL(8) :: T(NG_ORG),XMU(NG_ORG,NATT),XVAR(NG_ORG,NATT,NATT)
     $       ,DT(NG_ORG),DXMU(NG_ORG,NATT),DXVAR(NG_ORG,NATT,NATT)
     $       ,IXVAR(NG_ORG,NATT,NATT),DETXVAR(NG_ORG)
     $       ,DLLH_DT(NG_ORG),DLLH_DXMU(NG_ORG,NATT)
     $       ,DLLH_DXVAR(NG_ORG,NATT,NATT)
      ! Temporary variables:
      REAL(8) :: VV(NATT),MM(NATT,NATT),LHS_VAL

      ! Code begins
      GRAD_LLH= 0d0
      NPARBY2 = 0.5D0*(NATT+NATT*(NATT+1)/2)

#if defined(petsc)
      call MPI_Allreduce(NIND,NINDT,1,MPI_INTEGER
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
      NINDT = NIND
#endif

      ! First, Extract vector solution THETA to G-parameters
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $        ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),THETA)
      ! Extract DTHETA = FVAL = GVAL - THETA solution
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,DT(1:NG_ORG)
     $        ,DXMU(1:NG_ORG,1:NATT),DXVAR(1:NG_ORG,1:NATT,1:NATT),FVAL)
      ! Compute the gradient of LLH w.r.t \THETA = (T,XMU,XVAR) at iteration it
      DLLH_DT = 0D0
      DLLH_DXMU = 0D0
      DLLH_DXVAR = 0D0
      DO K = 1,NG_ORG
        IF (ALIVE_GRP(K) .EQ. 1) THEN
          ! First, we compute the inverse of Gaussian cov-matrix, XVAR(K,:,:)
          CALL GDETK(NATT,K,XVAR(K,1:NATT,1:NATT),IXVAR(K,1:NATT,1:NATT)
     $              ,DETXVAR(K),IER,NULL)
          ! Compute gradient LLH w.r.t Gaussian weights
          DLLH_DT(K) = WSUM(K)/T(K)-NPARBY2/T(K)-1D0*NINDT+NPARBY2*NG
          ! Compute gradient LLH w.r.t Gaussian means
          DLLH_DXMU(K,1:NATT)
     .         = MATMUL( IXVAR(K,1:NATT,1:NATT),WSUM(K)*DXMU(K,1:NATT) )
          ! Compute gradient LLH w.r.t. Gaussian covariance-matrices
          MM(1:NATT,1:NATT)
     .          =MATMUL( IXVAR(K,1:NATT,1:NATT)
     .                 ,WSUM(K)*DXVAR(K,1:NATT,1:NATT) )
          DLLH_DXVAR(K,1:NATT,1:NATT)=0.5D0
     .          *MATMUL( MM(1:NATT,1:NATT),IXVAR(K,1:NATT,1:NATT) )
        END IF    ! end if ALIVE_GRP
      END DO
      CALL FORM_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,DLLH_DT(1:NG_ORG)
     $         ,DLLH_DXMU(1:NG_ORG,1:NATT)
     .         ,DLLH_DXVAR(1:NG_ORG,1:NATT,1:NATT),GRAD_LLH)
      ! Compute the dot-product of GRAD_LLH & THETA_AA - THETA_{it}
      LHS_VAL = DOT_PRODUCT( GRAD_LLH,THETA_AA-THETA )

      MONOTONICITY = (LHS_VAL .GE. -1D0*EPS_VAL)

      END FUNCTION CHECK_MONOTONICITY
    
!c     COMPUTE_PLLH_VER2
!c     ###############################################################
!c     This function compute the penalized LLH (same as COMPUTE_LLH subroutine)
!c     but it takes the GMM parameters and form GMM PDF (GUM) and then compute
!c     the value of penalized LLH
!c     ###############################################################
      SUBROUTINE COMPUTE_PLLH_VER2(ITER,NIND,G_DOF,NG_ORG,NG,NATT
     $            ,WX,X,THETA,XLOGL_OUT,MLENG_OUT,AL,GUM,ALIVE_GRP)
!   Given the Gaussian parameters vector U & the data points X(1:NIND,1:NATT), we
!	compute the XLOGL, the LOGLIKELIHOOD value, of Gaussian mixture at iteration ITER
!   During the process of computing XLOGL, we actually computed responsibilities, 
!   pdf of X(J,:) in Gaussian(K), AL(J,K), and pdf of GMM, GUM(J) J=1,...,NIND & K=1,...,NG_ORG.
!   So we output them (AL & GUM) as well

      IMPLICIT NONE

      !	Call variables
      INTEGER, INTENT(IN) :: NIND,G_DOF,NG_ORG,NG,NATT,ITER
      INTEGER, INTENT(IN) :: ALIVE_GRP(NG_ORG)
      REAL(8), INTENT(IN) :: THETA(G_DOF*NG_ORG),X(NIND,NATT),WX(NIND)
      REAL(8), INTENT(OUT) :: AL(NIND,NG_ORG),GUM(NIND),XLOGL_OUT
     $       ,MLENG_OUT
      !	Local variables
      INTEGER :: I,J,K,II,JJ,KK,IERR,NUL,NG_ALIVE
      REAL(8) :: CN,detcov,p_ik,PENALIZED_TERM,LXLOGL
      INTEGER, PARAMETER :: DENMAX=175
      REAL(8), PARAMETER :: PI=3.14159265358979D0 
      REAL(8) :: NPARBY2,DVAL,T(NG_ORG),XMU(NG_ORG,NATT)
     $       ,XVAR(NG_ORG,NATT,NATT),V(NG_ORG,NATT,NATT)
     $       ,DV(NG_ORG),EVAL(NG_ORG,NATT,NATT),EVEC(NATT,NATT)
     $       ,SQRT_DV(NG_ORG),NPARBY2_LOGT
     $       ,XMUK(NATT),VK(NATT,NATT),ISQDV,TK

      ! Code begins
      CN=1d0/((2d0*PI)**(NATT*0.5d0))
      GUM(1:NIND) = 0D0
      AL(1:NIND,1:NG_ORG) = 0D0
      XLOGL_OUT = 0d0

      ! Extract Gaussian parameters (of GMM)
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $         ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),THETA)

      ! Obtain the inverse (V) and determinant (DV) of each group's cov-matrix
      DO K=1,NG_ORG               ! group
        IF (ALIVE_GRP(K) .EQ. 1) THEN
          CALL GDETK(NATT,K,XVAR(K,1:NATT,1:NATT),V(K,1:NATT,1:NATT)
     $                  ,DV(K),IERR,NUL)
          !------------------------------------------------------

          XMUK = XMU(K,:)
          VK   = V(K,:,:)
          ISQDV= CN/SQRT(DV(K))
          TK   = T(K)
          DO JJ=1,NIND          ! particle
            CALL Mahalanobis(NATT,X(JJ,:),XMUK,VK,p_ik)
            p_ik = (-0.5d0)*MIN(p_ik,1d0*DENMAX)
            AL(JJ,K) = ISQDV*TK*EXP(p_ik)
            GUM(JJ) = GUM(JJ) + AL(JJ,K)
          END DO                ! end particle
        !------------------------------------------------------
        END IF
      END DO                      ! end group

      ! Compute penalized LLH (taking into account particles' weights)
      CALL compute_pllh(NATT,T,GUM,WX,ALIVE_GRP,XLOGL_OUT,MLENG_OUT)

      if (DIAG_IO) write(ML_IO,*) ITER,XLOGL_OUT,MLENG_OUT
     .                           ,MINVAL(T,mask=(ALIVE_GRP==1))

      END SUBROUTINE COMPUTE_PLLH_VER2

!c     RELAX_ACCEL_SOLN
!c     ###############################################################
      SUBROUTINE RELAX_ACCEL_SOLN(G_DOF,NG,NATT,CUR_ITER,THETA,GVAL
     $                           ,RELAX_FLAG,ALIVE_GRP)
!     Perform relaxation on accelerated solutions if positivity (of
!     Gaussian weights) is not preserved.  Since we accelerate XVAR's
!     Cholesky decomp. XVAR(K,:,:) = L(K,:,:) * L(K,:,:)^T , we do not
!     care about positivity of XVAR.

      IMPLICIT NONE

      !	Call variables
      INTEGER, INTENT(IN) :: G_DOF,NG,NATT,CUR_ITER
      INTEGER, INTENT(IN) :: ALIVE_GRP(NG)
      REAL(8), INTENT(IN) :: GVAL(G_DOF*NG)
      REAL(8), INTENT(INOUT) :: THETA(G_DOF*NG)
      LOGICAL, INTENT(OUT) :: RELAX_FLAG

      !	Some controllable Variables
      LOGICAL, PARAMETER :: PRINT_MODE = .FALSE.

      !   Local Variables    
      INTEGER :: I,J,K,KK
      REAL(8) :: TOL,TOL1,TOL2,DXVAR(NG),BETA,BETA0,DET_XVAR
      REAL(8) :: T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT)
      REAL(8) :: T_GVAL(NG),XMU_GVAL(NG,NATT),XVAR_GVAL(NG,NATT,NATT)
      REAL(8) :: Y(G_DOF*NG),EVEC(NATT,NATT),EVAL(NATT),MIN_EVAL
      REAL(8) :: XVAR_TEMP(NATT,NATT)

      !	Variables for DSYEV
      INTEGER :: LWORK, DSYEV_INFO
      REAL(8) :: WORK(3*NATT-1)

      !	Code begins
      LWORK = 100!3*NATT-1
      RELAX_FLAG = .FALSE.

      !First, we extract Gauss params from sol vectors GVAL & THETA
      !CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,T_GVAL,XMU_GVAL,XVAR_GVAL,GVAL)
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,T,XMU,XVAR,THETA)
      ! Sanity checks:
      !   If a Gaussian weight < 0, or
      !   If there is a non-pd cov-matrix/2nd-moment, or
      !   If there is a cov-matrix/2nd-moment that is close to rank-deficient, or
      !   Then reroll to EM solution: GVAL --> THETA
      TOL1 = 1.0D-9  ! used for checking Gaussian weights & eigenvalues of cov-matrices/2nd-moments
      TOL2 = 1.0D-13  ! used for checking determinant of cov-matrices/2nd-moments

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ! First, we test Gaussian weights or zeroth moments
      DO K = 1,NG
        IF ((ALIVE_GRP(K) .EQ. 1) .AND. (T(K) .LE. TOL1)) THEN
          IF (PRINT_MODE) THEN
            PRINT*, "ITER =",CUR_ITER,". T(K) =",T(K),
     $            ". RELAX_ACCEL_SOLN: NEGATIVE WEIGHT!"
          END IF
          THETA = GVAL
          RELAX_FLAG = .TRUE.
          RETURN
        END IF
      END DO

!c      ! Second, we test if covariance matrices have small determinant:
!c      DO K = 1,NG
!c        IF ((ALIVE_GRP(K) .EQ. 1)) THEN
!c          DET_XVAR = M33DET(XVAR(K,:,:))
!c          IF (DET_XVAR .LT. TOL2) THEN
!c            THETA = GVAL
!c            RELAX_FLAG = .TRUE.
!c            RETURN
!c          END IF 
!c        END IF
!c      END DO

      END SUBROUTINE RELAX_ACCEL_SOLN
    
!c     ENFORCE_DEAD_GAUSSIAN
!c     ###############################################################   
      SUBROUTINE ENFORCE_DEAD_GAUSSIAN(G_DOF,NG,NATT,THETA,ALIVE_GRP)
      ! This subroutine will enforce dead Gaussians after AA since AA can bring
      !   a group that has weight = 0 to having very small weights or negative weights

      IMPLICIT NONE

      ! Call variables
      INTEGER, INTENT(IN) :: G_DOF,NG,NATT
      INTEGER, INTENT(IN) :: ALIVE_GRP(NG)
      REAL(8), INTENT(INOUT) :: THETA(G_DOF*NG)
      ! Local variables
      INTEGER :: K
      REAL(8) :: T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT)

      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,T,XMU,XVAR,THETA)

      DO K = 1,NG
        IF (ALIVE_GRP(K) .EQ. 0) THEN
          T(K) = 0D0
          XMU(K,1:NATT) = 0D0
          XVAR(K,1:NATT,1:NATT) = 0D0
        END IF
      END DO
      T(1:NG) = T(1:NG)/SUM(T(1:NG))

      CALL FORM_GPARAMS_VECTOR(G_DOF,NG,NATT,T,XMU,XVAR,THETA)

      END SUBROUTINE ENFORCE_DEAD_GAUSSIAN
    
!---------------------------------------------------------------------------------------
! 2 subroutines that will perform XVAR(K,:,:) -> L(K,:,:) & L(K,:,:) = XVAR(K,:,:)
!   where L(K,:,:) is Cholesky fact. of XVAR(K,:,:)
!---------------------------------------------------------------------------------------

!c    CHOLESKY_XVAR_TRANSFORM
!c    ###############################################################   
      SUBROUTINE CHOLESKY_XVAR_TRANSFORM(G_DOF,NG_ORG,NATT
     $                  ,UL,UOUT,ALIVE_GRP)
      ! Input: vectors UL stores weights, means & L (chol of cov-matrices)
      ! Ouput: vectors UOUT stores weights, means & cov-matrices

      IMPLICIT NONE

      !   Call variables
      INTEGER, INTENT(IN) :: G_DOF,NG_ORG,NATT
      INTEGER, INTENT(IN) :: ALIVE_GRP(NG_ORG)
      REAL(8), INTENT(IN) :: UL(G_DOF*NG_ORG)
      REAL(8), INTENT(OUT) :: UOUT(G_DOF*NG_ORG)

      !   Local variables
      INTEGER :: I,J,K
      REAL(8) :: T(NG_ORG),XVAR(NG_ORG,NATT,NATT)
     $      ,XMU(NG_ORG,NATT),LOWER_XVAR(NG_ORG,NATT,NATT)


      ! Code begins
      LOWER_XVAR = 0D0
      ! First extract UL to Gaussian weights, means, and lower-triangular matrices:
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $      ,XMU(1:NG_ORG,1:NATT),LOWER_XVAR(1:NG_ORG,1:NATT,1:NATT),UL)
      ! Compute Cholesky fact. of XVAR(1:NG_ORG,:,:) & GXVAR(1:NG_ORG,:,:):
      DO K=1,NG_ORG
        IF (ALIVE_GRP(K) .EQ. 1) THEN
          !----------------------------------------------------
          ! Since LOWER_XVAR is stored as symmetric matrices, before doing anything with it
          ! we enforce the 'strictly-upper-triangular' portion of it to be zero
          DO I=1,NATT
            DO J=I+1,NATT
              LOWER_XVAR(K,I,J)=0D0
            END DO
          END DO
          !----------------------------------------------------
          XVAR(K,1:NATT,1:NATT) = 0D0
          CALL DGEMM('N','T',NATT,NATT,NATT,1.0D0
     .         ,LOWER_XVAR(K,1:NATT,1:NATT)
     .         ,NATT,LOWER_XVAR(K,1:NATT,1:NATT),NATT,0.0D0
     .         ,XVAR(K,1:NATT,1:NATT),NATT)
        END IF
      END DO
      ! Form UOUT = (T,XMU,LOWER_XVAR) & GVALOUT = (GT,GXMU,LOWER_GXVAR)
      UOUT = 0D0
      CALL FORM_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $    ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),UOUT)
    
      END SUBROUTINE CHOLESKY_XVAR_TRANSFORM

    
    
!c     XVAR_CHOLESKY_TRANSFORM
!c     ###############################################################   
      SUBROUTINE XVAR_CHOLESKY_TRANSFORM(G_DOF,NG_ORG,NATT
     $                  ,THETA,GVAL,UOUT,GVALOUT,ALIVE_GRP)
      ! Input: vectors THETA & GVAL that stores Gaussians parameters at iteration IT & IT+1, respectively
      ! Ouput: vectors UOUT & GVALOUT that stores weights, means & L (chol of cov-matrices)
      !        DF = GVALOUT-UOUT := residual vector of GVALOUT, UOUT

      IMPLICIT NONE
    
      !   Call variables
      INTEGER, INTENT(IN) :: G_DOF,NG_ORG,NATT
      INTEGER, INTENT(IN) :: ALIVE_GRP(NG_ORG)
      REAL(8), INTENT(IN) :: THETA(G_DOF*NG_ORG),GVAL(G_DOF*NG_ORG)
      REAL(8), INTENT(OUT) :: UOUT(G_DOF*NG_ORG),GVALOUT(G_DOF*NG_ORG)
    
      !   Local variables
      INTEGER :: I,J,K
      REAL(8) :: T(NG_ORG),XMU(NG_ORG,NATT),XVAR(NG_ORG,NATT,NATT)
     $    ,GT(NG_ORG),GXMU(NG_ORG,NATT),GXVAR(NG_ORG,NATT,NATT)
     $    ,LOWER_XVAR(NG_ORG,NATT,NATT),LOWER_GXVAR(NG_ORG,NATT,NATT)


      ! Code begins
      LOWER_XVAR = 0D0
      LOWER_GXVAR = 0D0
      ! First extract U & GVAL to Gaussian parameters:
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $    ,XMU(1:NG_ORG,1:NATT),XVAR(1:NG_ORG,1:NATT,1:NATT),THETA)
      CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,GT(1:NG_ORG)
     $    ,GXMU(1:NG_ORG,1:NATT),GXVAR(1:NG_ORG,1:NATT,1:NATT),GVAL)
      ! Compute Cholesky fact. of XVAR(1:NG_ORG,:,:) & GXVAR(1:NG_ORG,:,:):
      DO K=1,NG_ORG
        IF (ALIVE_GRP(K) .EQ. 1) THEN
          !CALL CHOLESKY_FACT(XVAR(K,:,:),NATT,LOWER_XVAR(K,:,:))
          !CALL CHOLESKY_FACT(GXVAR(K,:,:),NATT,LOWER_GXVAR(K,:,:))
          CALL M33CHOLESKLY(XVAR(K,1:NATT,1:NATT),NATT
     .                    ,LOWER_XVAR(K,1:NATT,1:NATT))
          CALL M33CHOLESKLY(GXVAR(K,1:NATT,1:NATT),NATT
     .                     ,LOWER_GXVAR(K,1:NATT,1:NATT))
          ! Truong: note that the L matrices LOWER_XVAR & LOWER_GXVAR 
          !   returned here are symmetrix of the lower-triangular matrices
          !   but in the following, we only store the lower-triangular part
          !   so we are all good here.
        END IF
      END DO
      ! Form UOUT = (T,XMU,LOWER_XVAR) & GVALOUT = (GT,GXMU,LOWER_GXVAR)
      ! Truong: the routine FORM_GPARAMS_VECTOR will only store the lower-triangular
      !   part of the (covariance/lower) matrices
      CALL FORM_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,T(1:NG_ORG)
     $       ,XMU(1:NG_ORG,1:NATT),LOWER_XVAR(1:NG_ORG,1:NATT,1:NATT)
     .       ,UOUT)
      CALL FORM_GPARAMS_VECTOR(G_DOF,NG_ORG,NATT,GT(1:NG_ORG)
     $       ,GXMU(1:NG_ORG,1:NATT),LOWER_GXVAR(1:NG_ORG,1:NATT,1:NATT)
     .       ,GVALOUT)
    
      END SUBROUTINE XVAR_CHOLESKY_TRANSFORM

!---------------------------------------------------------------------------------------

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   End of 'necessary' subroutines needed for A-AMEM.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      end module AMEM_mod

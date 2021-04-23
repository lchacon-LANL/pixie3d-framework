       ! Module for AA with regularization
       MODULE REG_AA_MOD

       PRIVATE  ! Make everything private
       PUBLIC :: REGULARIZED_AA,STANDARD_AA,RESTART_AA,AA_CREATE
     .          ,AA_DESTROY

       PUBLIC :: TEST_SVD_LAPACK

       !******************************************************************************************************************!
       !******************************************************************************************************************!
       ! Anderson acceleration driver with/without regularization & restart
       ! Authors: L. Chacon, G. Chen, T. Nguyen (07/2020)
       ! This module contains implementation for 'standard' Anderson acceleration and Anderson
       !     acceleration with regularization, periodical restart, and monotonicity control for EM
       ! The Anderson acceleration implementation follows closely the paper written by Prof. Homer Walker in 2011
       ! The regularized Anderson acceleration follows the paper written by Nicholas Henderson & Ravi Varadhan in 2018
       ! Refs:
       ! [1] H. F. Walker and P. Ni. Anderson acceleration for fixed-point iterations. SIAM J. Numer. Anal., 49:1715–1735, 2011
       ! [2] Nicholas C. Henderson & Ravi Varadhan (2019) Damped Anderson Acceleration With Restarts and Monotonicity Control 
       ! for Accelerating EM and EM-like Algorithms, Journal of Computational and Graphical Statistics, 28:4, 834-846,
       ! DOI: 10.1080/10618600.2019.1594835

       ! Updated (08/2020):
       ! The code for AA follows the implementation style of FPA module by N.N.Carlson.
       ! I.e., we create a data type for AA variables. We call it AA_STATE variable.
       !
       !---------------------- Usage -----------------------
       !       
       ! Given MVEC    = max number of past residuals,
       !       AASTART = delayed parameter for AA,
       !       N       = nrow of solution vector X
       !
       ! CALL AA_CREATE(STATE,MVEC,AASTART,N) 
       ! DO ITER=1,2,... UNTIL CONVERGED:
       !   GVAL = PICARD(X)                                    ! Fixed point/Picard iteration
       !   CALL REGULARIZED_AA(STATE,ITER,N,X,GVAL)            ! CUR_IT = current iteration index
       !   IF MOD(ITER,MVEC) == 0 CALL RESTART_AA(STATE)       ! (optional)
       ! END DO
       ! CALL AA_DESTROY(STATE)
       !
       !----------------------------------------------------
       !******************************************************************************************************************!
       !******************************************************************************************************************!
    
       !######################################################################
       ! Create AA_STATE
       !######################################################################
       TYPE, PRIVATE :: REGULARIZED_VARS
         REAL(8) :: LAMDA=0d0,RVAL=0d0,ALPHA=0d0,KAPPA=0d0,SK=0d0      ! parameters to obtain regularized factor in AA's LSQ problem
       END TYPE REGULARIZED_VARS

       TYPE, PUBLIC :: AA_STATE
         INTEGER :: MVEC=0                         ! maximum number of past residuals used.
         INTEGER :: AASTART=0                      ! delayed parameter for starting AA.
         INTEGER :: MAA=0                          ! AA iteration counter, if AA is reset, MAA = 0
         INTEGER :: N=0                            ! length of solution vector
         REAL(8), POINTER :: FK(:,:) => null()     ! storage of past residuals (upto MAA)
         REAL(8), POINTER :: DG(:,:) => null()     ! storage of G(x^{it}) - G(x^{it-1})
         REAL(8), POINTER :: F_OLD(:) => null()    ! residual vector from 'previous' iteration 
         REAL(8), POINTER :: G_OLD(:) => null()    ! picard solution from 'previous' iteration
         TYPE(REGULARIZED_VARS) :: RAA_VARS        ! parameters to compute reg-factor in AA's LSQ problem (see Ref. [2])
       END TYPE AA_STATE

!c    Variables for A-AMEM configuration via namelist

       INTEGER,PUBLIC :: MVEC_AA =0          !# maximum number of past residualis ~5 or ~10
     .                  ,AASTART =0          !# delayed AA? - should use 0
       REAL(8),PUBLIC :: LAMDA   = 1.0D7     !# start with large values
     .                  ,KAPPA   = 25.0D0    !# value recommended by [Henderson&Varadhan2019]
     .                  ,ALPHA   = 1.5D0     !# value recommended by [Henderson&Varadhan2019]
     .                  ,SK      = 0D0       !# a parameter for regularized-AA
     .                  ,EPS_VAL = 3.0D-2    !# apparently for DPIC, should choose [0.1, 0.5] instead of ~0.01
     .                  ,EPS_C   = 1.0D-3    !# [Henderson&Varadhan2019] suggested to set to small
     .                  ,RVAL    = 1.0D0
     .                  ,DVAL    = 0.0D0
       LOGICAL,PUBLIC :: REGULARIZED_FLAG=.false.        !# regularization flag
     .                  ,DIAG_AMEM=.false.               !# diagnostics flag

!     List of variables to be read from aamem.in
      namelist /aamem_params/ MVEC_AA,AASTART,LAMDA,KAPPA,ALPHA
     .              ,SK,RVAL,DVAL,EPS_C,EPS_VAL,REGULARIZED_FLAG
     .              ,DIAG_AMEM

      PUBLIC aamem_params
      
      CONTAINS


       !######################################################################
       !	AA_CREATE
       !   Create state for AA
       !######################################################################
       SUBROUTINE AA_CREATE(STATE,MVEC,AASTART,N,REG_VARS)
       ! Inputs: MVEC = maximum number of residuals (specified by user)
       !         AASTART = delayed parameter for AA. AA won't start after 'AASTART' iterations
       !         N = size (number of rows) of solution vector
       !         REG_VARS = vector of 5 elements (specified by user) for initial values of 
       !                    the parameters that are needed to compute regularization-factor \lambda
       !                    in solving the AA's LSQ problem.
       !                    After 1st iter, at every iter, these parameters will be updated. (see Ref. [2])

       IMPLICIT NONE

       ! Call variables
       INTEGER, INTENT(IN) :: MVEC,AASTART,N
       REAL(8), INTENT(IN) :: REG_VARS(5)
       TYPE(AA_STATE), INTENT(OUT) :: STATE

       ! Code begins:
       STATE%N = N
       STATE%MAA = 0     ! initial number of residual = 0 (i.e., at iter = 0, we have 0 residual)
       STATE%MVEC = MVEC
       STATE%AASTART = AASTART
       ALLOCATE( STATE%FK(N,MVEC), STATE%DG(N,MVEC) )
       ALLOCATE( STATE%F_OLD(N), STATE%G_OLD(N) )

       ! parameters for regularized AA (again, see Ref. [2])
       ! one should store the regularization parameters in the following orders:
       STATE%RAA_VARS%LAMDA = REG_VARS(1)
       STATE%RAA_VARS%RVAL  = REG_VARS(2)
       STATE%RAA_VARS%ALPHA = REG_VARS(3)
       STATE%RAA_VARS%KAPPA = REG_VARS(4)
       STATE%RAA_VARS%SK    = REG_VARS(5)

       END SUBROUTINE AA_CREATE
    
       !######################################################################
       !	AA_DESTROY
       !   Create state for AA
       !######################################################################
       SUBROUTINE AA_DESTROY(STATE)

       IMPLICIT NONE

       ! Call variables
       TYPE(AA_STATE), INTENT(INOUT) :: STATE

       ! Code begins:
       DEALLOCATE( STATE%FK,STATE%DG,STATE%F_OLD,STATE%G_OLD )

       END SUBROUTINE AA_DESTROY
      
       !######################################################################
       !	REGULARIZED_AA
       !   Main implementation of regularized-AA
       !   See Ref. [2] for the implementation details
       !######################################################################
       SUBROUTINE REGULARIZED_AA(STATE,CUR_IT,X,GVAL)

       ! This implementation of regularized-AA follows the implementation of H.Walker's AA
       ! The ONLY difference is the Tikkonov-regularization term.

       ! Here is how it goes:
       !   If cur_it == 0 or cur_it < AASTART, we simply set xAA = Picard iterate, G(x)
       !   If MAA == 0, also set xAA = G(x)
       !   Otherwise, solve the LSQ problem with regularized term (coming from EXT_RAA_VARS)
       !       EXT_RAA_VARS consists of variables that needed for regularization.
       !       EXT_RAA_VARS should be of length 5, EXT_RAA_VARS = [LAMDA,RVAL,ALPHA,KAPPA,SK]

       ! Some 'other' variables explaination
       ! STATE = AA_STATE variables, see AA_CREATE subroutine
       ! CUR_IT = current iteration counter in the main iterative function
       !          i.e., this subroutine needs to be called at every iteration!
       ! N = first dimension of FK,DG or dimension of FVAL,F_OLD,GVAL,G_OLD
       !          (In GMMEM context, it is the total number of parameters in GMM)
       ! X = solution at current iteration, on output, X = updated (AA) solution
       ! GVAL = G(X) = picard update of X


       IMPLICIT NONE

       ! Call variables
       TYPE(AA_STATE), INTENT(INOUT) :: STATE  ! AA state variable
       INTEGER, INTENT(IN) :: CUR_IT
       REAL(8), INTENT(IN) :: GVAL(STATE%N)
       REAL(8), INTENT(INOUT) :: X(STATE%N)    ! Input is solution from previous iteration
                                               ! Output can either be GVAL or AA solution

       ! Local variables
       INTEGER :: I,J,K,IT,N,MAA
       REAL(8) :: LAMDA,DELTA_K
       REAL(8) :: FVAL(STATE%N)
       REAL(8) :: DF(STATE%N),MATRIX_A(STATE%MVEC,STATE%MVEC)
     $   ,VECTOR_B(STATE%MVEC),TEMPVEC(STATE%N)
     $   ,GAMMA(STATE%MVEC),ID_MAT(STATE%MVEC,STATE%MVEC)

       ! Variables for DGESV (LAPACK routine)
       INTEGER :: IPIV(STATE%MVEC),INFO

       ! Code begins

       IT = CUR_IT
       N = STATE%N                  ! Do this to shorten the code
       MAA = STATE%MAA              ! Do this to shorten the code

       ! Compute residual vector: G(x) - x, GVAL - X
       FVAL = GVAL - X 
              
       ! Anderson acceleration in H.Walker style
       IF ( (IT .LT. STATE%AASTART) .OR. (STATE%MVEC .EQ. 0) ) THEN
         X = GVAL
       ELSE
         ! Apply AA. Compute the DF vector := FVAL-F_OLD
         IF (IT .GT. STATE%AASTART) THEN
           DF = FVAL - STATE%F_OLD
           IF (MAA .LT. STATE%MVEC) THEN
             STATE%DG(1:N,MAA+1) = GVAL - STATE%G_OLD
           ELSE
             STATE%DG(1:N,1:MAA-1) = STATE%DG(1:N,2:MAA)
             STATE%DG(1:N,MAA) = GVAL - STATE%G_OLD
           END IF
           MAA = MAA + 1
         END IF
         ! Updating F_OLD & G_OLD for the next iteration
         STATE%F_OLD = FVAL
         STATE%G_OLD = GVAL

         IF (MAA .EQ. 0) THEN
           X = GVAL
         ELSE
           IF (MAA .GT. STATE%MVEC) THEN
             STATE%FK(1:N,1:MAA-2) = STATE%FK(1:N,2:MAA-1)
             MAA = MAA - 1
           END IF
           STATE%FK(1:N,MAA) = DF
           ! Now, we are ready to solve LSQ problem (with regularization)
           DELTA_K = 1D0/(1D0 +
     $  STATE%RAA_VARS%ALPHA**(STATE%RAA_VARS%KAPPA-STATE%RAA_VARS%SK))
           CALL COMPUTE_REGULARIZED_FACTOR(N,MAA,STATE%RAA_VARS%LAMDA
     $      ,STATE%RAA_VARS%RVAL,STATE%FK(:,1:MAA),FVAL,DELTA_K
     $      ,STATE%RAA_VARS%ALPHA,STATE%RAA_VARS%KAPPA,STATE%RAA_VARS%SK
     $      ,MATRIX_A(1:MAA,1:MAA),VECTOR_B(1:MAA))
           ! From the call of subroutine COMPUTE_REGULARIZED_FACTOR(),
           ! the MATRIX_A := FK^T * Fk and the VECTOR_B := Fk^T * FVAL are also computed

           ! Form identity matrix ID_MAT of size (MAA x MAA)
           ID_MAT(1:MAA,1:MAA) = 0D0
           DO I = 1,MAA
             ID_MAT(I,I) = 1D0
           END DO

           ! MATRIX_A & VECTOR_B, output from subroutine COMPUTE_REGULARIZED_FACTOR, are defined to be
           !   MATRIX_A = FK^T * FK & VECTOR_B = FK^T * FVAL
           ! Now, form the matrix FK^T * FK + \lambda*I ==> MATRIX_A
           MATRIX_A(1:MAA,1:MAA) = MATRIX_A(1:MAA,1:MAA)
     $                       + STATE%RAA_VARS%LAMDA*ID_MAT(1:MAA,1:MAA)

           ! Then solve for gamma = MATRIX_A \ VECTOR_B in the system
           !   (FK^T * FK + \lambda*I) * \gamma = FK^T * FVAL
           GAMMA(1:MAA) = VECTOR_B(1:MAA)
           CALL DGESV(MAA,1,MATRIX_A(1:MAA,1:MAA),MAA,IPIV(1:MAA)
     $      ,GAMMA(1:MAA),MAA,INFO)
           ! Now, obtain X --- the accelerated solution using the values of gamma
           TEMPVEC(1:N) = 0D0
           CALL DGEMV('N',N,MAA,1.0D0,STATE%DG(1:N,1:MAA),N
     $      ,GAMMA(1:MAA),1,0D0,TEMPVEC(1:N),1)
           X = GVAL - TEMPVEC
         END IF    
       END IF
       STATE%MAA = MAA ! Update the AA state variable MAA
       END SUBROUTINE REGULARIZED_AA
    
    
       !######################################################################
       !	STANDARD_AA
       !   Main implementation of regularized-AA
       !######################################################################
       SUBROUTINE STANDARD_AA(STATE,CUR_IT,X,GVAL)

       ! This implementation of regularized-AA follows the implementation of H.Walker AA

       ! Here is how it goes:
       !   If cur_it == 0 or cur_it < AASTART, we simply set xAA = Picard iterate, G(x)
       !   If MAA == 0, also set xAA = G(x)
       !   Otherwise, solve the LSQ problem WITHOUT regularized term

       ! Some 'other' variables explaination
       ! CUR_IT = current iteration counter in the main iterative function
       !          i.e., this subroutine needs to be called at every iteration!
       ! N = first dimension of FK,DG or dimension of FVAL,F_OLD,GVAL,G_OLD
       !          in GMMEM context, it is the total number of parameters in GMM
       ! MAA = current number of residuals
       ! MVEC = maximum number of residuals AA can use
       ! AASTART = parameter for delaying AA
       ! FK (see paper by H.Walker AA in 2011 i.e., ref[1])
       ! DG stores history of G(x) - x in its columns
       ! FVAL = G(x) - x
       ! F_OLD = previous G(x) - x
       ! GVAL = G(x)
       ! G_OLD = previous G(x)


       IMPLICIT NONE

       ! Call variables
       TYPE(AA_STATE), INTENT(INOUT) :: STATE      ! AA state variable
       INTEGER, INTENT(IN) :: CUR_IT
       REAL(8), INTENT(IN) :: GVAL(STATE%N)
       REAL(8), INTENT(INOUT) :: X(STATE%N) ! The AA solution (presummably a better solution than GVAL)

       ! Local variables
       INTEGER :: I,J,K,IT,N,MAA
       REAL(8) :: DF(STATE%N),MATRIX_A(STATE%MVEC,STATE%MVEC)
     $   ,VECTOR_B(STATE%MVEC),TEMPVEC(STATE%N)
     $   ,GAMMA(STATE%MVEC),FVAL(STATE%N)

       ! Variables for DGESV (LAPACK routine)
       INTEGER :: IPIV(STATE%MVEC),INFO


       ! Code begins

       N = STATE%N         ! Do this to shorten the code
       FVAL = GVAL - X     ! compute residual vector
       IT = CUR_IT
       MAA = STATE%MAA     ! Do this to shorten the code

       ! Anderson acceleration in H.Walker style
       IF ( (IT .LT. STATE%AASTART) .OR. (STATE%MVEC .EQ. 0) ) THEN
         X = GVAL
       ELSE
         ! Apply AA. Compute the DF vector := FVAL-F_OLD
         IF (IT .GT. STATE%AASTART) THEN
           DF = FVAL - STATE%F_OLD
           IF (MAA .LT. STATE%MVEC) THEN
             STATE%DG(1:N,MAA+1) = GVAL - STATE%G_OLD
           ELSE
             STATE%DG(1:N,1:MAA-1) = STATE%DG(1:N,2:MAA)
             STATE%DG(1:N,MAA) = GVAL - STATE%G_OLD
           END IF
           MAA = MAA + 1
         END IF
         ! Updating F_OLD & G_OLD for the next iteration
         STATE%F_OLD = FVAL
         STATE%G_OLD = GVAL

         IF (MAA .EQ. 0) THEN
           X = GVAL
         ELSE
           IF (MAA .GT. STATE%MVEC) THEN
             STATE%FK(1:N,1:MAA-2) = STATE%FK(1:N,2:MAA-1)
             MAA = MAA - 1
           END IF
           STATE%FK(1:N,MAA) = DF

           ! Now, we are ready to solve LSQ problem (no regularization)
           ! See netlib.org documentation for DGEMM & DGEMV
           ! Form the MATRIX_A := (Fk^T)*Fk
           MATRIX_A(1:MAA,1:MAA) = 0D0
           VECTOR_B(1:MAA) = 0D0
           CALL DGEMM('T','N',MAA,MAA,N,1.0D0,STATE%FK(1:N,1:MAA),N
     $      ,STATE%FK(1:N,1:MAA),N,0.0D0,MATRIX_A(1:MAA,1:MAA),MAA)
           ! Form the vector VECTOR_B := (Fk^T)*fval
           CALL DGEMV('T',N,MAA,1.0D0,STATE%FK(1:N,1:MAA),N,FVAL,1
     $      ,0.0D0,VECTOR_B(1:MAA),1)

           ! Then solve for gamma = MATRIX_A \ VECTOR_B in the system
           !   (FK^T * FK) * \gamma = FK^T * FVAL
           GAMMA(1:MAA) = VECTOR_B(1:MAA)
           CALL DGESV(MAA,1,MATRIX_A(1:MAA,1:MAA),MAA,IPIV(1:MAA)
     $                  ,GAMMA(1:MAA),MAA,INFO)
           ! Now, obtain XAA --- the accelerated solution using the values of gamma
           TEMPVEC(1:N) = 0D0
           CALL DGEMV('N',N,MAA,1.0D0,STATE%DG(1:N,1:MAA),N,GAMMA(1:MAA)
     $      ,1,0D0,TEMPVEC(1:N),1)
           X = GVAL - TEMPVEC
         END IF    
       END IF
       STATE%MAA = MAA     ! update state variable's MAA
       END SUBROUTINE STANDARD_AA
    
    
    
    
    
    
    
    
       !!----------------------------------------------------------------------------------------!!
       !! The following subroutines are necessary for DAAREM subroutine
       !!----------------------------------------------------------------------------------------!!



       !######################################################################
       !	FUNCTION_S
       !######################################################################
       SUBROUTINE FUNCTION_S(N,MAA,LAMDA,FK,FVAL,SVAL,A,VECB)
       !   Input: FK(1:N,1:MAA) = storage for fval - fval_old
       !          (See H.Walker AA paper (ref [1]) for the definition of FK)
       !          N is the 1st dimension of FK, in adaptive EM context, N = G_DOF*NG_ORG
       !          LAMDA = regularized factor in the LSQ problem for computing AA coefficients
       !          FVAL = G(X) - G(X_OLD) = residual
       !          A = FK^T * FK
       !          VECB = FK^T * FVAL
       !   Output: SVAL: = S(LAMDA)

       IMPLICIT NONE

       ! Call variables
       INTEGER, INTENT(IN) :: N,MAA
       REAL(8), INTENT(IN) :: LAMDA,FK(N,MAA),FVAL(N)
     $       ,A(MAA,MAA),VECB(MAA)
       REAL(8), INTENT(OUT) :: SVAL(MAA)

       ! Local variables
       INTEGER :: I,J
       REAL(8) :: ID_MAT(MAA,MAA),MAT_A(MAA,MAA)

       ! Variables needed for LAPACK routines
       INTEGER :: IPIV(MAA),INFO

       ! Code begins
       ! First, build identity matrix of size MAA-by-MAA
       ID_MAT = 0D0
       DO I = 1,MAA
         ID_MAT(I,I) = 1D0
       END DO

       ! Since A = transpose(FK) * FK
       ! MAT_A = transpose(FK)*FK + \lambda*I
       MAT_A(1:MAA,1:MAA) = A(1:MAA,1:MAA) + LAMDA*ID_MAT(1:MAA,1:MAA)

       ! Compute S(LAMDA) = SVAL = (A^-1) * VECB
       SVAL(1:MAA) = VECB(1:MAA)
       CALL DGESV(MAA,1,MAT_A,MAA,IPIV,SVAL,MAA,INFO)
       ! On output SVAL is the solution!

       END SUBROUTINE FUNCTION_S
    

       !######################################################################
       !	FUNCTION_PHI
       !######################################################################
       SUBROUTINE FUNCTION_PHI(MAA,LAMDA,DELTA_K,NORM_S_LAMDA
     $                  ,NORM_S_ZERO,PHI)
       ! This function takes outputs computed by subroutine FUNCTION_S()
       !   to compute PHI = NORM(S_LAMDA) - V_K
       !   where V_K = sqrt(DELTA_K) * NORM(S_ZERO)
       ! NOTE: S_LAMDA = FUNCTION_S(...,LAMDA,...)
       !       S_ZERO  = FUNCTION_S(...,0D0,...)

       IMPLICIT NONE

       ! Call variables
       INTEGER, INTENT(IN) :: MAA
       REAL(8), INTENT(IN) :: LAMDA,DELTA_K,NORM_S_LAMDA,NORM_S_ZERO
       REAL(8), INTENT(OUT) :: PHI ! return value for function phi(lamda)

       ! Local variables
       REAL(8) :: V_K

       ! Code begins
       V_K = SQRT(DELTA_K)*NORM_S_ZERO
       PHI = NORM_S_LAMDA - V_K

       END SUBROUTINE FUNCTION_PHI

    
       !######################################################################
       !	FUNCTION_DPHI
       !######################################################################
       SUBROUTINE FUNCTION_DPHI(N,MAA,LAMDA,DK,UK_DOT_FVAL,NORM_S_LAMDA
     $                          ,DPHI)
       ! This subroutine is derivative of function PHI above w.r.t LAMDA

       IMPLICIT NONE

       ! Call variables
       INTEGER, INTENT(IN) :: N,MAA
       REAL(8), INTENT(IN) :: LAMDA,NORM_S_LAMDA
     $   ,DK(MAA),UK_DOT_FVAL(MAA) 
       ! DK = singular values of FK, UK_DOT_FVAL stores the values of dot_product of (columns of UK) & (FVAL)
       REAL(8), INTENT(OUT) :: DPHI    ! return value for function dphi(lamda)

       ! Local variables
       INTEGER :: K,L
       REAL(8) :: CC,XX,YY,FK_COPY(N,MAA)

       ! Code begins
       ! Compute derivative of PHI, DPHI, w.r.t lamda
       CC = -1.0D0/NORM_S_LAMDA
       DPHI = 0D0
       DO L = 1,MAA
         XX = DK(L)*UK_DOT_FVAL(L)
         YY = DK(L)*DK(L) + LAMDA
         DPHI = DPHI + (XX*XX)/(YY*YY*YY)
       END DO
       DPHI = CC*DPHI

       END SUBROUTINE FUNCTION_DPHI
    
    
    
    
       !######################################################################
       !	COMPUTE_REGULARIZED_FACTOR
       !######################################################################
       SUBROUTINE COMPUTE_REGULARIZED_FACTOR(N,MAA,LAMDA,RVAL,FK,FVAL
     $               ,DELTA_K,ALPHA,KAPPA,SK,MATRIX_A,VECTOR_B)
       ! This subroutine will compute the Tikkonov-regularized factor for the LSQ problem
       ! that is needed to compute the AA-coeffs.
       ! Again, in EM/adaptive EM context, N = G_DOF*NG_ORG = first dimension of FK
       !                                   MAA = second dimension of Fk
       ! LAMDA is the regularized factor, it is input & output
       ! RVAL also is input and output
       ! FK (See H.Walker AA paper (ref [1]) for the definition of FK)
       ! FVAL = G(THETA) - THETA = residual vector at current iteration
       ! DELTA_K, ALPHA, KAPPA, SK are parameters for regularized-restart-monotonicity-AAEM
       ! Output: LAMDA & RVAL & the precomputed FK^T * FK, FK^T * FVAL

       IMPLICIT NONE

       ! Call variables
       INTEGER, INTENT(IN) :: N,MAA
       REAL(8), INTENT(IN) :: FK(N,MAA),FVAL(N)
     $   ,DELTA_K,ALPHA,KAPPA,SK
       REAL(8), INTENT(INOUT) :: LAMDA,RVAL
       REAL(8),INTENT(OUT) :: MATRIX_A(MAA,MAA),VECTOR_B(MAA)

       ! Local variables
       INTEGER :: I,J,K,MAXIT
       REAL(8) :: LAMDA_OLD,R_OLD,PHI_OF_LAMDA
     $   ,DPHI_OF_LAMDA,NORM_BETA_LS,V_K,LVAL,UVAL
     $   ,LSTOP,USTOP,PHI_OF_ZERO,DPHI_OF_ZERO,VECB(MAA)
     $   ,LowerBound,UpperBound
       REAL(8) :: S_ZERO(MAA),S_LAMDA(MAA),PHI_OVER_DPHI
     $   ,NORM_S_ZERO,NORM_S_LAMDA
       REAL(8) :: FK_COPY(N,MAA),UK(N,N),DK(MAA),VKT(MAA,MAA)
     $   ,UK_DOT_FVAL(MAA)

       ! Variables for SVD-LAPACK subroutine
       INTEGER, PARAMETER :: LWORK = 10000         ! set to be large integer?
       INTEGER :: INFO
         REAL(8) :: WORK(LWORK)


       ! Code begins

       ! To optimize the code, we evaluate some precomputed variables here:
       !   MATRIX_A = transpose(FK) * FK
       !   VECTOR_B = transpose(FK) * FVAL
       ! These variables will be wisely used in FUNCTION_S:
       MATRIX_A(1:MAA,1:MAA) = 0D0
       VECTOR_B(1:MAA) = 0D0
       ! See netlib.org documentation for DGEMM & DGEMV
       ! Form the MATRIX_A := (Fk^T)*Fk
       CALL DGEMM('T','N',MAA,MAA,N,1.0D0,FK,N,FK,N,0.0D0,MATRIX_A,MAA)
       ! Form the vector VECTOR_B := (Fk^T)*fval
       CALL DGEMV('T',N,MAA,1.0D0,FK,N,FVAL,1,0.0D0,VECTOR_B,1)

       ! Also, compute SVD of FK for only ONE time!
       FK_COPY = FK
       CALL DGESVD('A','N',N,MAA,FK_COPY,N,DK,UK,N,VKT,MAA
     $                  ,WORK,LWORK,INFO)
       ! Then compute the dot_product(UK(:,L),FVAL) for each L=1,...,MAA
       ! to be used in FUNCTION_DPHI()
       DO I = 1,MAA
         UK_DOT_FVAL(I) = DOT_PRODUCT(UK(:,I),FVAL)
       END DO

       ! Step 0: re-compute r_old using value of lamda_old
       LAMDA_OLD = LAMDA
       R_OLD = RVAL
       !   Compute S(0D0) & S(LAMDA_OLD)
       CALL FUNCTION_S(N,MAA,0D0,FK,FVAL,S_ZERO,MATRIX_A,VECTOR_B)
       CALL FUNCTION_S(N,MAA,LAMDA_OLD,FK,FVAL,S_LAMDA
     $                  ,MATRIX_A,VECTOR_B)
       NORM_S_LAMDA = NORM2(S_LAMDA)
       NORM_S_ZERO  = NORM2(S_ZERO)
       !   Compute PHI(LAMDA) & DPHI(LAMDA)
       CALL FUNCTION_PHI(MAA,LAMDA_OLD,DELTA_K,NORM_S_LAMDA
     $                  ,NORM_S_ZERO,PHI_OF_LAMDA)
       CALL FUNCTION_DPHI(N,MAA,LAMDA_OLD,DK,UK_DOT_FVAL
     $                  ,NORM_S_LAMDA,DPHI_OF_LAMDA)
       R_OLD = NORM2(S_LAMDA)*PHI_OF_LAMDA/DPHI_OF_LAMDA

       ! Step 1: initialization of lamda_0:
       V_K = SQRT(DELTA_K)*NORM_S_ZERO
       LAMDA = LAMDA_OLD - R_OLD/V_K
       !   Initialization of L_0 & U_0:
       CALL FUNCTION_PHI(MAA,0D0,DELTA_K,NORM_S_ZERO
     $                  ,NORM_S_ZERO,PHI_OF_ZERO)
       CALL FUNCTION_DPHI(N,MAA,0D0,DK,UK_DOT_FVAL
     $                  ,NORM_S_ZERO,DPHI_OF_ZERO)
       LVAL = -1.0D0*PHI_OF_ZERO/DPHI_OF_ZERO
       !   Compute transpose(FK) * FVAL = VECB
       CALL DGEMV('T',N,MAA,1.0D0,FK,N,FVAL,1,0.0D0,VECB,1)
       UVAL = NORM2(VECB)/V_K
       !   Compute stopping parameters for L & U
       LSTOP = 1.0D0/SQRT( 1.0D0 + ALPHA**(KAPPA-SK+0.5) )
       USTOP = 1.0D0/SQRT( 1.0D0 + ALPHA**(KAPPA-SK-0.5) )

       ! Step 2: iterate until convergence to find LAMDA
       MAXIT = 100 
       DO I = 1,MAXIT
         ! Update value of LAMDA
         IF ((LAMDA .LE. LVAL) .OR. (LAMDA .GE. UVAL)) THEN
           LAMDA = MAX(0.001D0*UVAL,SQRT(UVAL*LVAL))
         END IF
         ! Evaluate S(LAMDA_T), PHI(LAMDA_T), DPHI(LAMDA_T)
         CALL FUNCTION_S(N,MAA,LAMDA,FK,FVAL,S_LAMDA,MATRIX_A,VECTOR_B)
         NORM_S_LAMDA = NORM2(S_LAMDA)
         CALL FUNCTION_PHI(MAA,LAMDA,DELTA_K,NORM_S_LAMDA
     $                          ,NORM_S_ZERO,PHI_OF_LAMDA)
         CALL FUNCTION_DPHI(N,MAA,LAMDA,DK,UK_DOT_FVAL
     $                          ,NORM_S_LAMDA,DPHI_OF_LAMDA)
         LowerBound = LSTOP*NORM_S_ZERO
         UpperBound = USTOP*NORM_S_ZERO
         PHI_OVER_DPHI = PHI_OF_LAMDA/DPHI_OF_LAMDA
         ! Terminating condition
         IF ((NORM_S_LAMDA .GE. LowerBound) .AND. 
     $       (NORM_S_LAMDA .LE. UpperBound)) THEN
           RVAL = NORM_S_LAMDA*PHI_OVER_DPHI
           EXIT
         END IF
         ! Update LVAL & UVAL for iteration I+1
         IF ( PHI_OF_LAMDA .LT. 0D0 ) UVAL = LAMDA
         LVAL = MAX(LVAL,LAMDA - PHI_OVER_DPHI)
         ! Update LAMDA for iteration I+1
         LAMDA = LAMDA - PHI_OVER_DPHI*NORM_S_LAMDA/V_K
       END DO

       !Step 3: If LAMDA doesn't converge, then compute RVAL using current value of LAMDA
       CALL FUNCTION_S(N,MAA,LAMDA,FK,FVAL,S_LAMDA,MATRIX_A,VECTOR_B)
       NORM_S_LAMDA = NORM2(S_LAMDA)
       CALL FUNCTION_PHI(MAA,LAMDA,DELTA_K,NORM_S_LAMDA
     $                  ,NORM_S_ZERO,PHI_OF_LAMDA)
       CALL FUNCTION_DPHI(N,MAA,LAMDA,DK,UK_DOT_FVAL
     $                  ,NORM_S_LAMDA,DPHI_OF_LAMDA)
       RVAL = NORM2(S_LAMDA)*PHI_OF_LAMDA/DPHI_OF_LAMDA

       ! At this point, LAMDA & RVAL should be output
       END SUBROUTINE COMPUTE_REGULARIZED_FACTOR
    
    
    
       !######################################################################
       !	RESTART_AA
       !######################################################################
       SUBROUTINE RESTART_AA(STATE)
       ! This should routine sets everything that is needed for AA to zero
       ! for restarting AA from 0 residual
       ! DG(1:N,1:MVEC) = pre-allocated storage for history of G(theta) - theta

       IMPLICIT NONE

       ! Call variables
       TYPE(AA_STATE), INTENT(INOUT) :: STATE

       STATE%MAA = 0
       STATE%DG(:,:) = 0D0
       STATE%FK(:,:) = 0D0

       !STATE%MAA = 1
       !STATE%DG(:,1) = STATE%DG(:,STATE%MVEC)
       !STATE%FK(:,1) = STATE%FK(:,STATE%MVEC)
       !STATE%DG(:,2:STATE%MVEC) = 0D0
       !STATE%FK(:,2:STATE%MVEC) = 0D0

       END SUBROUTINE RESTART_AA
    
    
    
       !!----------------------------------------------------------------------------------------!!
       !! The following subroutines are for testing necessary routines
       !!----------------------------------------------------------------------------------------!!

       !######################################################################
       !	TEST_SVD_LAPACK
       !######################################################################
       SUBROUTINE TEST_SVD_LAPACK(M,N,A,U,D,VT)
       IMPLICIT NONE
       ! Call variables
       INTEGER, INTENT(IN) :: M,N
       REAL(8), INTENT(IN) :: A(M,N)
       REAL(8), INTENT(OUT) :: U(M,M),D(N),VT(N,N)
       ! Local variables
       INTEGER :: I,J
       ! Variables for LAPACK subroutine
       INTEGER, PARAMETER :: LWORK = 10000
       INTEGER :: INFO
       REAL(8) :: WORK(LWORK)

       ! See netlib.org documentation for DGESVD
       CALL DGESVD('A','N',M,N,A,M,D,U,M,VT,N,WORK,LWORK,INFO)
       WRITE(*,*) "PRINT U(:,:)"
       DO I = 1,M
         WRITE(*,*) (U(I,J),J=1,M)
       END DO
       WRITE(*,*) "PRINT D(:)"
       WRITE(*,*) (D(I),I=1,N)
       WRITE(*,*) "PRINT VT(:,:)"
       DO I = 1,N
         WRITE(*,*) (VT(I,J),J=1,N)
       END DO
       END SUBROUTINE TEST_SVD_LAPACK


       !************************************************************************************************
       ! END OF REG_AA_MOD MODULE
       !************************************************************************************************
       END MODULE REG_AA_MOD

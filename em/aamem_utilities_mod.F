       MODULE AAMEM_UTILITIES_MOD

       PRIVATE
       PUBLIC :: FORM_GPARAMS_VECTOR,EXTRACT_GPARAMS_VECTOR,
     $    GPARAMS_TO_MOMENTS,GMOMENTS_TO_PARAMS,
     $    PRINT_GPARAMS,LOGLIKELIHOOD_IO,
     $    MONOTONICITY_DIAGNOSTIC_IO,REGULARIZED_FACTOR_IO,
     $    KILL_GAUSS_HISTORY_IO,MCHOLESKY,
     $    WRITE_SAMPLES_IO,WRITE_GPARAMS_IO

       CONTAINS




       !######################################################################
       !     FORM_GPARAMS_VECTOR
       !######################################################################
       SUBROUTINE FORM_GPARAMS_VECTOR(G_DOF,NG,NATT,T,XMU,XVAR,U)
       !     This subroutine takes parameters from NG Gaussian and put them in a
       !     column vector U = [T(1:NG)'; XMU(1:NG,:); XVAR(1:NG,1:6)]
       !     NOTE: XVAR in vector form is a .5*NATT*(NATT+1)-by-1 vector, NATT=3

       IMPLICIT NONE

       !   Call variables
       INTEGER, INTENT(IN) :: NG,NATT,G_DOF
       REAL(8), INTENT(IN) :: T(NG),XMU(NG,NATT)
       REAL(8), INTENT(IN) :: XVAR(NG,NATT,NATT)
       REAL(8), INTENT(OUT) :: U(G_DOF*NG)

       !   Local variables
       INTEGER :: I,J,K,IDX

       IDX=0
       !   Insert NG Gaussian weights to solution vector:
       DO K=1,NG
        IDX=IDX+1
        U(K) = T(K)
       END DO

       !   Insert NG Gaussian means to solution vector:
       DO K=1,NG
        DO I=1,NATT
          IDX=IDX+1
          U(IDX) = XMU(K,I)
        END DO
       END DO   

       !   Insert NG Gaussian cov to solution vector (only using 6 components):
       DO K=1,NG
        DO I=1,NATT
          DO J=1,I
            IDX=IDX+1
            U(IDX) = XVAR(K,I,J)
          END DO
        END DO
       END DO

       IF (IDX .NE. G_DOF*NG) THEN
        WRITE(*,*) 'FORM_GPARAMS_VECTOR: MISSING ELT IN SOLN VECTOR!'
        STOP
       ENDIF

       ! NOTE: The first NG entries of U will contain the weights of Gaussians, U(1:NG)
       !		The next  NG*NATT entries of U will contain the means of Gaussians, U(NG+1:NG*NATT)
       !		The next  NG*( NATT*(NATT+1)/2 ) of U will contain the sigmas of Gaussians, U(NG*NATT+1:END)
       !		If U holds the moments of Gaussians, the order is the same

       END SUBROUTINE FORM_GPARAMS_VECTOR



       !######################################################################
       !     EXTRACT_GPARAMS_VECTOR
       !######################################################################
       SUBROUTINE EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT
     $			,TOUT,XMUOUT,XVAROUT,U)
        !    This subroutine extract parameters from a vector U to Gaussians'
        !    parameters: T, XMU, XVAR
        !    NOTE: column vector U = [T(1:NG)'; XMU(1:NG,:); XVAR(1:NG,1:6)]
        !    NOTE: XVAR in vector form is a .5*NATT*(NATT+1)-by-1 vector, NATT=3
        !    In particular, U(1:NG)          -> T(1:NG)
        !                   U(NG+1:NG*NATT)  -> XMU(1:NG,1:NATT)
        !                   U(NG*NATT+1:END) -> XVAR(1:NG,1:NATT,1:NATT)

       IMPLICIT NONE

       !   Call variables
       INTEGER, INTENT(IN) :: NG,NATT,G_DOF
       REAL(8), INTENT(IN) :: U(G_DOF*NG)
       REAL(8), INTENT(OUT) :: TOUT(NG),XMUOUT(NG,NATT)
     $			,XVAROUT(NG,NATT,NATT)

       !   Local variables
       INTEGER :: I,J,K,IDX

       IDX=0
       !   Extract NG Gaussians' weights from solution vector:
       DO K=1,NG
        IDX=IDX+1
        TOUT(K) = U(K)
       END DO

       !   Extract NG Gaussians' means from solution vector:
       DO K=1,NG
        DO I=1,NATT
          IDX=IDX+1
          XMUOUT(K,I) = U(IDX)
        END DO
       END DO   

       !   Extract NG Gaussians' cov matrices from solution vector:
       !   NOTE: cov matrices are symmetric
       DO K=1,NG
        DO I=1,NATT
          DO J=1,I
            IDX=IDX+1
            XVAROUT(K,I,J) = U(IDX)
            XVAROUT(K,J,I) = XVAROUT(K,I,J)
          END DO
        END DO
       END DO

       IF (IDX .NE. G_DOF*NG) THEN
        WRITE(*,*) 'EXTRACT_GPARAMS_VECTOR: MISSING ELT IN SOLN VECTOR!'
        STOP
       ENDIF

       END SUBROUTINE EXTRACT_GPARAMS_VECTOR
    
    
       !######################################################################
       !     GPARAMS_TO_MOMENTS
       !######################################################################
       SUBROUTINE GPARAMS_TO_MOMENTS(G_DOF,NG,NATT,U,MOMS)
!    This subroutine is used to convert Gauss parameters (weights,means,sigmas)
!    into Gaussian moments
!    Note (08/24/2020): Decided not to accelerate Gaussian moments for adaptive AMEM
!                   So this subroutine is not used
    
       IMPLICIT NONE

       !   Call variables
       INTEGER, INTENT(IN) :: G_DOF,NG,NATT
       REAL(8), INTENT(IN) :: U(G_DOF*NG)
       REAL(8), INTENT(OUT) :: MOMS(G_DOF*NG)

       !   Local variables
       INTEGER :: I,J,K
       REAL(8) :: T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT)
     $          ,M0(NG),M1(NG,NATT),M2(NG,NATT,NATT)

       MOMS = 0D0

       ! First we extract the Gaussian params from the vector solution U
       CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,T,XMU,XVAR,U)

       M0 = 0D0
       M1 = 0D0
       M2 = 0D0
       ! Computing the moments from gauss params:
       ! For each Gaussian, we know that
       !   Zeroth moment = weight
       M0(1:NG) = T(1:NG)
       !   First moment = weight*mu
       !   Second moment = weight*( sigma + mu*mu^T )
       DO K = 1,NG
        DO I = 1,NATT
          M1(K,I) = T(K)*XMU(K,I)
          DO J = 1,I
            M2(K,I,J) = T(K)*( XVAR(K,I,J) + XMU(K,I)*XMU(K,J) )
            M2(K,J,I) = M2(K,I,J)
          END DO
        END DO
       END DO

       ! Finally, packing M0, M1, M2 into vector MOMS for output
       CALL FORM_GPARAMS_VECTOR(G_DOF,NG,NATT,M0,M1,M2,MOMS)

       END SUBROUTINE GPARAMS_TO_MOMENTS
    
    
   
       !######################################################################
       !     GMOMENTS_TO_PARAMS
       !######################################################################
       SUBROUTINE GMOMENTS_TO_PARAMS(G_DOF,NG,NATT,MOMS,U)
!    This subroutine is used to convert Gauss moments (M0,M1,M2)
!    into Gaussian parameters (weights,means,sigmas)
!    Note (08/24/2020): Decided not to accelerate Gaussian moments for adaptive AMEM
!                      So this subroutine is not used
    
       IMPLICIT NONE

       !   Call variables
       INTEGER, INTENT(IN) :: G_DOF,NG,NATT
       REAL(8), INTENT(IN) :: MOMS(G_DOF*NG)
       REAL(8), INTENT(OUT) :: U(G_DOF*NG)

       !   Local variables
       INTEGER :: I,J,K
       REAL(8) :: T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT),
     $          M0(NG),M1(NG,NATT),M2(NG,NATT,NATT)


       U = 0D0

       ! First, we extract the moments from solution vector MOMS:
       CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,M0,M1,M2,MOMS)

       ! We then compute the G-params from the moments:
       T(1:NG) = M0(1:NG)
       DO K = 1,NG
        DO I = 1,NATT
          XMU(K,I) = M1(K,I)/T(K)
          DO J = 1,I
            XVAR(K,I,J) = M2(K,I,J)/T(K) - XMU(K,I)*XMU(K,J)
            XVAR(K,J,I) = XVAR(K,I,J)
          END DO
        END DO
       END DO

       ! Form the solution vector U consisting of T,XMU,XVAR
       CALL FORM_GPARAMS_VECTOR(G_DOF,NG,NATT,T,XMU,XVAR,U)

       END SUBROUTINE GMOMENTS_TO_PARAMS
    
    
!######################################################################
!       MCHOLESKY
!       Cholesky factorization for square matrix using CHOLA
!######################################################################        
       SUBROUTINE MCHOLESKY(A,N,L)

       IMPLICIT NONE

       !       Call variables
       INTEGER, INTENT(IN) :: N
       REAL(8), INTENT(IN) :: A(N,N)   ! input matrix for computing cholesky
       REAL(8), INTENT(OUT) :: L(N,N)  ! output lower-triangular matrix
                                     ! BUT in the form of symmetric matrix

       !       Local variables
       INTEGER :: I,J,K

       select case(N)
       case(3)
         L(:,:) = 0D0
         L(1,1) = SQRT(A(1,1))
         L(2,1) = A(2,1)/L(1,1)
         L(3,1) = A(3,1)/L(1,1)

         L(1,2) = L(2,1)
         L(2,2) = SQRT( A(2,2)-L(2,1)*L(2,1) )
         L(3,2) = (A(3,2)-L(3,1)*L(2,1))/L(2,2)

         L(1,3) = L(3,1)
         L(2,3) = L(3,2)
         L(3,3) = SQRT( A(3,3)-L(3,1)*L(3,1)-L(3,2)*L(3,2) )
       case(2)
         L(:,:) = 0D0
         L(1,1) = SQRT(A(1,1))
         L(2,1) = A(2,1)/L(1,1)

         L(1,2) = L(2,1)
         L(2,2) = SQRT( A(2,2)-L(2,1)*L(2,1) )
       case(1)
         L(1,1) = SQRT(A(1,1))
       case default
         WRITE(*,*) 'Cannot factorize A in M33Cholesky'
         STOP
       end select
       
       END SUBROUTINE MCHOLESKY
    
    
       !######################################################################
       !     PRINT_GPARAMS
       !######################################################################
       SUBROUTINE PRINT_GPARAMS(G_DOF,NG,NATT,EM_PARAMS,AAEM_PARAMS)
       !   This subroutine acts as I/O function to double check & compare the results
       !   returned by adaptive EM & adaptive AMEM algorithms

       IMPLICIT NONE

       !   Call variables
       INTEGER, INTENT(IN) :: G_DOF,NG,NATT
       REAL(8), INTENT(IN) :: EM_PARAMS(G_DOF*NG),AAEM_PARAMS(G_DOF*NG)

       !   Local variables
       INTEGER :: I,J,K
       REAL(8) :: T_em(NG),XMU_em(NG,NATT),XVAR_em(NG,NATT,NATT)
       REAL(8) :: T_aa(NG),XMU_aa(NG,NATT),XVAR_aa(NG,NATT,NATT)



       ! Here, for print out, we extract standard EM solution vector EM_PARAMS
       !   to (T_em,XMU_em,XVAR_em) & extract AAEM solution vector AAEM_PARAMS
       !   to (T_aa,XMU_aa,XVAR_aa) so we can print them out and compare.
       CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,T_em,XMU_em
     $          ,XVAR_em,EM_PARAMS)
       CALL EXTRACT_GPARAMS_VECTOR(G_DOF,NG,NATT,T_aa,XMU_aa
     $          ,XVAR_aa,AAEM_PARAMS)

       PRINT*,"  "
       PRINT*,"------------------------------------------------"
       PRINT*, "Final weights EM (left col) & AAEM (right col) methods"
       DO K = 1,NG
          WRITE(*,*) T_em(K), T_aa(K)
       END DO
       PRINT*,"------------------------------------------------"
       PRINT*, "Final means from EM(left col) & AAEM(right col) methods"
       DO K = 1,NG
        DO I = 1,NATT
          WRITE(*,*) XMU_em(K,I), XMU_aa(K,I)
        END DO
        WRITE(*,*) " "
       END DO
       PRINT*,"------------------------------------------------"
       PRINT*, "Final cov-matrices from EM(left col) & 
     $                  AAEM(right col) methods"
       DO K = 1,NG
        DO I = 1,NATT
        DO J = 1,I
            WRITE(*,*) XVAR_em(K,I,J), XVAR_aa(K,I,J)
        END DO
        END DO
        WRITE(*,*) " "
       END DO


       END SUBROUTINE PRINT_GPARAMS
    
    
    
    
       !######################################################################
       !     WRITE_GPARAMS_IO
       !######################################################################
       SUBROUTINE WRITE_GPARAMS_IO(FILENAME,NIND,G_DOF
     $                  ,NG,NATT,EM_PARAMS)
       IMPLICIT NONE

       !   Call variables
       INTEGER, INTENT(IN) :: NIND,G_DOF,NG,NATT
       REAL(8), INTENT(IN) :: EM_PARAMS(G_DOF*NG)
       CHARACTER(*), INTENT(IN) :: FILENAME

       !   Local variables
       INTEGER :: I,J,K    

       OPEN(UNIT=99,FILE=FILENAME,ACTION='WRITE')
       WRITE(99,*) NIND
       WRITE(99,*) G_DOF
       WRITE(99,*) NG
       WRITE(99,*) NATT
       DO I=1,G_DOF*NG
        WRITE(99,*) EM_PARAMS(I)
       END DO
       CLOSE(99) 


       END SUBROUTINE WRITE_GPARAMS_IO
    
    
    
       !-----------------------------------------------------------------------------
       ! I/O SUBROUTINES FOR WRITING STUFF OUT TO FILES FOR THE DIAGNOSTIC PURPOSES
       !-----------------------------------------------------------------------------

       !######################################################################
       !     LOGLIKELIHOOD_IO
       !######################################################################
       SUBROUTINE LOGLIKELIHOOD_IO(MITER,XLOGL,ITER,FILENAME)
       !   This subroutine write the history of (penalized) log-likelihood as 
       !   a function of iterations to a file

       IMPLICIT NONE

       !	Call variables
       INTEGER, INTENT(IN) :: MITER,ITER
       REAL(8), INTENT(IN) :: XLOGL(MITER)
       CHARACTER(*), INTENT(IN) :: FILENAME
       !	Local variables
       INTEGER :: I,J

       OPEN(UNIT=99,FILE=FILENAME,ACTION='WRITE')
       DO I=1,ITER
          WRITE(99,*) I, XLOGL(I)
       END DO
       CLOSE(99)

       END SUBROUTINE LOGLIKELIHOOD_IO
    
    
    
    
       !######################################################################
       !     MONOTONICITY_DIAGNOSTIC_IO
       !######################################################################
       SUBROUTINE MONOTONICITY_DIAGNOSTIC_IO(MITER,MONOTONICITY_HIST
     $                  ,ITER,FILENAME)
       !   This subroutine write the history of A-AMEM's monotonicity control (binary values of 0 or 1) 
       !   as a function of iterations to a file

       IMPLICIT NONE

       !	Call variables
       INTEGER, INTENT(IN) :: MITER,ITER
       INTEGER, INTENT(IN) :: MONOTONICITY_HIST(MITER)
       CHARACTER(*), INTENT(IN) :: FILENAME
       !	Local variables
       INTEGER :: I,J

       OPEN(UNIT=99,FILE=FILENAME,ACTION='WRITE')
       DO I=1,ITER
          WRITE(99,*) I, MONOTONICITY_HIST(I)
       END DO
       CLOSE(99)

       END SUBROUTINE MONOTONICITY_DIAGNOSTIC_IO
    
    
    
    
       !######################################################################
       !     REGULARIZED_FACTOR_IO
       !######################################################################
       SUBROUTINE REGULARIZED_FACTOR_IO(MITER,REG_FACTOR,ITER,FILENAME)
       !   This subroutine write the history of A-AMEM's penalization factor \lambda
       !   as a function of iterations to a file

       IMPLICIT NONE

       !	Call variables
       INTEGER, INTENT(IN) :: MITER,ITER
       REAL(8), INTENT(IN) :: REG_FACTOR(MITER)
       CHARACTER(*), INTENT(IN) :: FILENAME
       !	Local variables
       INTEGER :: I,J

       OPEN(UNIT=99,FILE=FILENAME,ACTION='WRITE')
       DO I=1,ITER
          WRITE(99,*) I, REG_FACTOR(I)
       END DO
       CLOSE(99)

       END SUBROUTINE REGULARIZED_FACTOR_IO



       !######################################################################
       !     KILL_GAUSS_HISTORY_IO
       !######################################################################
       SUBROUTINE KILL_GAUSS_HISTORY_IO(NG_ORG,KILL_GAUSS_HIST,FILENAME)
       !   This subroutine write the iteration^{th} where a Gaussian is killed to a file
       !   for both A-EM and A-AMEM method.

       IMPLICIT NONE

       !	Call variables
       INTEGER, INTENT(IN) :: NG_ORG
       INTEGER, INTENT(IN) :: KILL_GAUSS_HIST(NG_ORG)
       CHARACTER(*), INTENT(IN) :: FILENAME
       !   Local variables
       INTEGER :: I


       OPEN(UNIT=99,FILE=FILENAME,ACTION='WRITE')
       DO I=1,NG_ORG
        IF (KILL_GAUSS_HIST(I) .GT. 0) WRITE(99,*) KILL_GAUSS_HIST(I)
       END DO
       CLOSE(99)

       END SUBROUTINE KILL_GAUSS_HISTORY_IO
    
    
    
       !######################################################################
       !     WRITE_SAMPLES_IO
       !######################################################################
       SUBROUTINE WRITE_SAMPLES_IO(FILENAME,NIND,NATT,X)
       IMPLICIT NONE

       ! Call variables
       INTEGER, INTENT(IN) :: NIND,NATT
       REAL(8), INTENT(IN) :: X(NIND,NATT)
       CHARACTER(*), INTENT(IN) :: FILENAME
       ! Local variables
       INTEGER :: I,J


       OPEN(UNIT=99,FILE=FILENAME,ACTION='WRITE')
       DO I=1,NIND
        WRITE(99,*) X(I,:)
       END DO
       CLOSE(99)   


       END SUBROUTINE WRITE_SAMPLES_IO


       END MODULE AAMEM_UTILITIES_MOD

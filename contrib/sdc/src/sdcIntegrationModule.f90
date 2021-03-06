MODULE sdcIntegrationModule

  USE realKindModule, ONLY: &
    RealKind, &
    Zero
  USE sdcModule, ONLY: &
    sdcType, &
    createSDC, &
    destroySDC, &
    setSubNodes, &
    setIntegrationMatrix
!!$  USE sdcInterfaceModule, ONLY: &
!!$    computeProvisionalSolution, &
!!$    evaluateOdeRhs, &
!!$    computeCorrection

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
    stepSDC

CONTAINS


  SUBROUTINE stepSDC(my_rank,U,timeBegin,timeEnd,nSteps,nSubNodes,nCorrections,iout &
                    ,SDC_init_solution,SDC_eval_correction,SDC_eval_rhs,SDC_eval_residual)

    ! -- Input:
    !
    !  my_rank      : parallel rank
    !  U            : vector of evolved fields
    !  timeBegin    : time corresponding to the input state of U
    !  timeEnd      : time state corresponding to the output state of U
    !  nSteps       : number of steps to update U from timeBegin to timeEnd
    !  nSubNodes    : number of SDC sub-nodes (including both end points; >= 2)
    !  nCorrections : number of SDC corrections applied
    !  iout         : indicates level of output (0 -> no output, >0 -> increasing levels of output)
    !
    !  SDC_init_solution: routine to evaluate initial solution
    !  SDC_eval_correction: routine to compute SDC correction
    !  SDC_eval_rhs (optional) : routine to evaluate ODE RHS
    !  SDC_eval_residual (optional): routine to evaluate SDC residual
    !
    ! -- Output:
    !
    !  U            : vector of evolved fields at endTime

    !Call variables

    INTEGER, INTENT(in) :: &
      my_rank, &
      nSteps, &
      nSubNodes, &
      nCorrections, &
      iout
    REAL(kind=RealKind), INTENT(in) :: &
      timeBegin, &
      timeEnd
    REAL(kind=RealKind), DIMENSION(:), INTENT(inout) :: &
      U

    OPTIONAL :: SDC_eval_residual,SDC_eval_RHS

    INTERFACE
       SUBROUTINE SDC_init_solution(iout,nn,nt,t,U)
         USE realKindModule, ONLY: RealKind
         INTEGER :: iout,nn,nt
         REAL(kind=RealKind), DIMENSION(nt) :: t
         REAL(kind=RealKind), DIMENSION(nn,nt) :: U
       END SUBROUTINE SDC_init_solution

       SUBROUTINE SDC_eval_residual(iout,nn,nt,t,U,M,Residual)
         USE realKindModule, ONLY: RealKind
         INTEGER :: iout,nn,nt
         REAL(kind=RealKind), DIMENSION(nt) :: t
         REAL(kind=RealKind), DIMENSION(nt,nt) :: M
         REAL(kind=RealKind), DIMENSION(nn,nt) :: U
         REAL(kind=RealKind), DIMENSION(nn,nt), INTENT(out) :: Residual
       END SUBROUTINE SDC_eval_residual

       SUBROUTINE SDC_eval_RHS(iout,nn,nt,t,U,RHS)
         USE realKindModule, ONLY: RealKind
         INTEGER :: iout,nn,nt
         REAL(kind=RealKind), DIMENSION(nt) :: t
         REAL(kind=RealKind), DIMENSION(nn,nt) :: U, RHS
       END SUBROUTINE SDC_eval_RHS

       SUBROUTINE SDC_eval_correction(iout,nn,nt,t,U,Residual,Correction)
         USE realKindModule, ONLY: RealKind
         INTEGER :: iout,nn,nt
         REAL(kind=RealKind), DIMENSION(nt) :: t
         REAL(kind=RealKind), DIMENSION(nn,nt) :: U, Residual, Correction
       END SUBROUTINE SDC_eval_correction
    END INTERFACE

    !Local variables

    INTEGER :: &
      iStep, &
      iSubNode, &
      iCorrection
    REAL(kind=RealKind) :: &
      dt, &
      time
    TYPE(sdcType), POINTER :: &
      SDC

    LOGICAL :: have_eval_res,have_eval_rhs

    !Begin program

    IF (nSubNodes < 2) then
       WRITE (*,*) ' nSubNodes in stepSDC must be >=2'
       WRITE (*,*) ' Aborting...'
       STOP
    ENDIF

    have_eval_res = PRESENT(SDC_eval_residual)
    have_eval_rhs = PRESENT(SDC_eval_rhs)

    IF (.not.(have_eval_res.or.have_eval_rhs)) then
       WRITE (*,*) ' Do not have routine to evaluate RHS or RESIDUAL'
       WRITE (*,*) ' Aborting...'
       STOP
    ENDIF

    dt   = (timeEnd-timeBegin) / REAL(nSteps)
    time = timeBegin

    CALL createSDC(U, nSubNodes, SDC)

    DO iStep = 1, nSteps

      IF (iout > 0.and. my_rank == 0) THEN
         PRINT*
         PRINT*, "  stepSDC, iStep=",iStep," time =", time
         PRINT*
      ENDIF

      CALL setSubNodes(time, time+dt, SDC%t)

      CALL setIntegrationMatrix(SDC%t, SDC%M)

      SDC%U (:,1) = U(:)

      IF (iout > 0 .and. my_rank == 0) THEN
        IF (iout > 1) write (*,*) '  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write (*,*) '  SDC >>> Computing provisional solution...'
        IF (iout > 1) write (*,*) '  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      ENDIF

      CALL SDC_init_solution(iout,SIZE(SDC%U,1),SIZE(SDC%U,2),SDC%t, SDC%U)

      DO iCorrection = 1, nCorrections

        IF (iout > 0.and. my_rank == 0) THEN
           IF (iout > 1) PRINT*, "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "
           PRINT*, "  SDC >>> iCorrection  =", iCorrection
           IF (iout > 1) PRINT*, "  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> "
        ENDIF

        IF (iout > 0.and. my_rank == 0) THEN
           IF (iout > 1) write (*,*)
           write (*,*) '   SDC >>> Computing residual...'
           IF (iout > 1) write (*,*)
        ENDIF

        if (have_eval_rhs) then
           CALL computeResidual(iout,SDC%t, SDC%U, SDC%M, SDC%Residual,SDC_eval_rhs)
        else
           CALL SDC_eval_residual(iout,SIZE(SDC%U,1), SIZE(SDC%U,2),SDC%t, SDC%U, SDC%M, SDC%Residual)
        endif

        IF (iout > 1 .and. my_rank == 0) THEN
           write (*,*)
           do iSubnode=2,nSubNodes
              PRINT*, "    SDC >>>  Max |Residual|/Max|Solution|, time   =", iSubNode,"is ", &
                MAXVAL(ABS(SDC%Residual(:,iSubNode)))/MAXVAL(ABS(SDC%U(:,iSubNode)))
           enddo
           write (*,*)
!!$           do iSubnode=2,nSubNodes
!!$              PRINT*, "    SDC >>>  Max |Residual|, time   =", iSubNode,"is ", &
!!$                MAXVAL(ABS(SDC%Residual(:,iSubNode)))
!!$           enddo
!!$           write (*,*)
        ENDIF


        IF (iout > 0 .and. my_rank == 0) THEN
           write (*,*) '   SDC >>> Computing correction...'
        ENDIF

        CALL SDC_eval_correction(iout,SIZE(SDC%U,1), SIZE(SDC%U,2),SDC%t, SDC%U, SDC%Residual, SDC%Correction)

        IF (iout > 1 .and. my_rank == 0) THEN
           PRINT*
           do iSubnode=2,nSubNodes
              PRINT*, "    SDC >>> Max|Correction|/Max|Solution|, time   =", iSubNode,"is ", &
                MAXVAL(ABS(SDC%Correction(:,iSubNode)))/MAXVAL((SDC%U(:,iSubnode)))
           enddo
           PRINT*
        ENDIF

!!$! diag ****
!!$!!        do iSubnode=2,nSubNodes
!!$           isubnode = nsubnodes
!!$           if (my_rank == 0) write (*,*) 'Plot Qtys itime=',isubnode
!!$           call plotQtys('debug.bin',size(SDC%U,1),time+dt &
!!$                     ,SDC%U(:,1),SDC%U(:,iSubNode)     &
!!$                     ,SDC%Residual(:,iSubNode)         &
!!$                     ,SDC%Correction(:,iSubNode))
!!$!!        enddo
!!$! diag ****

        SDC%U = SDC%U + SDC%Correction

      END DO

      time = time + dt

      U(:) = SDC%U(:, nSubNodes)

!!$      PRINT*
!!$      PRINT*, "    iStep, time =", iStep, time

    END DO

    CALL destroySDC(SDC)

  END SUBROUTINE stepSDC

  SUBROUTINE computeResidual(iout,t,U,M,Residual,evalRHS)

    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      t
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(in) :: &
      U, M
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(out) :: &
      Residual

    INTERFACE
       SUBROUTINE evalRHS(iout,nn,nt,t,U,RHS)
         USE realKindModule, ONLY: RealKind
         INTEGER :: iout,nn,nt
         REAL(kind=RealKind), DIMENSION(nt) :: t
         REAL(kind=RealKind), DIMENSION(nn,nt) :: U, RHS
       END SUBROUTINE evalRHS
    END INTERFACE

    INTEGER :: &
      i, j, iout
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      Integral
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      RHS

!   Begin program

    ALLOCATE(Integral(SIZE(U,DIM=1)),RHS(SIZE(U,DIM=1),SIZE(U,DIM=2)))

!   Evaluate ODE RHS

    CALL evalRHS(iout,SIZE(U,1),SIZE(U,2),t, U, RHS)

!   Compute residual at SDC nodes

    Residual(:,1) = Zero

    DO i = 2, SIZE(M,DIM=1)
      Integral(:) = Zero
      DO j = 1, SIZE(M,DIM=2)
        Integral(:) = Integral(:) + M(i,j) * RHS(:,j)
      END DO
      Residual(:,i) = U(:,1) + Integral(:) - U(:,i)
    END DO

    DEALLOCATE(Integral,RHS)

  END SUBROUTINE computeResidual

END MODULE sdcIntegrationModule

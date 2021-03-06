MODULE thermalWaveSolutionModule

  USE realKindModule, ONLY: &
    RealKind, &
    Zero, &
    Half, &
    One, &
    Two
  USE sdcModule, ONLY: &
    polynomialIntegral
  USE thermalWaveModule, ONLY: &
    nReactionSubCycles, &
    MeshScale, &
    Delta, &
    U_s

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
    updateThermalWaveBackwardEuler, &
    updateThermalWaveSplit, &
    updateThermalWaveMultiImplicit, &
    evaluateOdeRhsThermalWave, &
    computeCorrectionThermalWave, &
    computeCorrectionThermalWaveSplit, &
    computeCorrectionThermalWaveMI

  INTEGER, PARAMETER :: &
    OLD    = 1, &
    NEW    = 2, &
    CHANGE = 3

  TYPE :: thermalWaveImplicitType

    REAL(kind=RealKind) :: &
      timeStep, &
      meshScale
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      rhs
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      Temperature, &
      jacobian

  END TYPE thermalWaveImplicitType

CONTAINS


  SUBROUTINE updateThermalWaveBackwardEuler(t0, t1, Told, Tnew)

    REAL(kind=RealKind), INTENT(in) :: &
      t0, t1
    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      Told
    REAL(kind=RealKind), DIMENSION(:), INTENT(out) :: &
      Tnew

    LOGICAL :: &
      converged
    INTEGER :: &
      nZones, &
      iteration
    REAL(kind=RealKind) :: &
      newTemperature
    TYPE(thermalWaveImplicitType), POINTER :: &
      TWBE

    nZones = SIZE(Told)

    ALLOCATE(TWBE)
    CALL createThermalWaveImplicitType(TWBE, nZones)

    TWBE%meshScale          = MeshScale
    TWBE%timeStep           = (t1-t0)
    TWBE%temperature(:,OLD) = Told(:)

    ! --- Initial guess:
    TWBE%temperature(:,NEW) = TWBE%temperature(:,OLD)

    iteration = 0
    converged = .FALSE.
    DO WHILE(.NOT.converged)

      iteration = iteration + 1

      CALL populateLinearizedSystemRHS( &
             TWBE%temperature(:,OLD), TWBE%temperature(:,NEW), t0, &
             TWBE%timeStep, TWBE%meshScale, nZones, TWBE%rhs)

      CALL populateJacobianMatrix( &
             TWBE%temperature(:,NEW), TWBE%timeStep, TWBE%meshScale, &
             nZones, TWBE%jacobian)

      CALL solvePentaDiagonalLinearSystem( &
             TWBE%jacobian, TWBE%rhs, nZones, TWBE%temperature(:,CHANGE))

      !Update new temperature:
      TWBE%temperature(:,NEW) &
        = TWBE%temperature(:,NEW) + TWBE%temperature(:,CHANGE)

      WHERE(TWBE%temperature(:,NEW) < Zero) &
        TWBE%temperature(:,NEW) = 1.0E-16

      WHERE(TWBE%temperature(:,NEW) > One) &
        TWBE%temperature(:,NEW) = One - 1.0E-16

      CALL checkNonlinearConvergence(converged, iteration, &
             TWBE%temperature(:,CHANGE), TWBE%temperature(:,NEW), TWBE%rhs)

    END DO

    Tnew(:) = TWBE%temperature(:,NEW)

    CALL destroyThermalWaveImplicitType(TWBE)
    DEALLOCATE(TWBE)

  END SUBROUTINE updateThermalWaveBackwardEuler


  SUBROUTINE updateThermalWaveSplit(t0, t1, Told, Tnew, reverseOption, &
                                    semiAnalyticReactionUpdateOption)

    REAL(kind=RealKind), INTENT(in) :: &
      t0, t1
    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      Told
    REAL(kind=RealKind), DIMENSION(:), INTENT(out) :: &
      Tnew
    LOGICAL, OPTIONAL :: &
      reverseOption, &
      semiAnalyticReactionUpdateOption

    LOGICAL :: &
      reverse, &
      semiAnalyticReactionUpdate, &
      converged
    INTEGER :: &
      nZones, &
      iteration
    TYPE(thermalWaveImplicitType), POINTER :: &
      TWS

    reverse = .FALSE.
    IF(PRESENT(reverseOption)) &
      reverse = reverseOption

    semiAnalyticReactionUpdate = .FALSE.
    IF(PRESENT(semiAnalyticReactionUpdateOption)) &
      semiAnalyticReactionUpdate &
        = semiAnalyticReactionUpdateOption

    nZones = SIZE(Told)

    ALLOCATE(TWS)
    CALL createThermalWaveImplicitType(TWS, nZones)

    TWS%meshScale          = MeshScale
    TWS%timeStep           = (t1-t0)
    TWS%temperature(:,OLD) = Told(:)

    IF(reverse)THEN

      ! --- Update temperature with diffusion operator

      CALL updateWithDiffusionOperator(TWS, t0, nZones)

      TWS%temperature(:,OLD) = TWS%temperature(:,NEW)

    END IF

    ! --- Update temperature with reaction operator

    ! --- Initial guess:
    TWS%temperature(:,NEW) = TWS%temperature(:,OLD)

    iteration = 0
    converged = .FALSE.
    DO WHILE(.NOT.converged)

      iteration = iteration + 1

      IF(semiAnalyticReactionUpdate)THEN

        CALL populateRhsAndJacobianSmAnltc( &
               TWS%temperature(:,OLD), TWS%temperature(:,NEW), &
               TWS%timeStep, nZones, TWS%rhs, TWS%jacobian)

      ELSE

        CALL populateLinearizedSystemRHS( &
               TWS%temperature(:,OLD), TWS%temperature(:,NEW), t0, &
               TWS%timeStep, TWS%meshScale, nZones, TWS%rhs, &
               ignoreDiffusionOption = .TRUE.)

        CALL populateJacobianMatrix( &
               TWS%temperature(:,NEW), TWS%timeStep, TWS%meshScale, &
               nZones, TWS%jacobian, ignoreDiffusionOption = .TRUE.)

      END IF

      !Jacobian is diagonal, linearized system trivial to invert:
      TWS%temperature(:,CHANGE) &
        = TWS%rhs(:) / TWS%jacobian(:,0)
      !Update new temperature:
      TWS%temperature(:,NEW) &
        = TWS%temperature(:,NEW) &
            + TWS%temperature(:,CHANGE)

      WHERE(TWS%temperature(:,NEW) < Zero) &
        TWS%temperature(:,NEW) = 1.0E-16

      WHERE(TWS%temperature(:,NEW) > One) &
        TWS%temperature(:,NEW) = One - 1.0E-16

      CALL checkNonlinearConvergence(converged, iteration, &
             TWS%temperature(:,CHANGE), TWS%temperature(:,NEW), TWS%rhs)

    END DO

    IF(.NOT.reverse)THEN

      ! --- Update temperature with diffusion operator

      TWS%temperature(:,OLD) = TWS%temperature(:,NEW)

      CALL updateWithDiffusionOperator(TWS, t0, nZones)

    END IF

    Tnew(:) = TWS%temperature(:,NEW)

    CALL destroyThermalWaveImplicitType(TWS)
    DEALLOCATE(TWS)

  END SUBROUTINE updateThermalWaveSplit


  SUBROUTINE updateThermalWaveMultiImplicit(t0, t1, Told, Tnew, iSubNode)

    REAL(kind=RealKind), INTENT(in) :: &
      t0, t1
    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      Told
    REAL(kind=RealKind), DIMENSION(:), INTENT(out) :: &
      Tnew
    INTEGER, INTENT(in) :: &
      iSubNode

    LOGICAL :: &
      converged
    INTEGER :: &
      nZones, &
      iReactionSubCycle, &
      iteration
    REAL(kind=RealKind) :: &
      subTimeStep, &
      subTime0, &
      subTime1
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      diffusionOperator
    TYPE(thermalWaveImplicitType), POINTER :: &
      TWS

    nZones = SIZE(Told)

    ALLOCATE(TWS)
    CALL createThermalWaveImplicitType(TWS, nZones)

    TWS%meshScale          = MeshScale
    TWS%timeStep           = (t1-t0)
    TWS%temperature(:,OLD) = Told(:)

    ! --- Update temperature with diffusion operator

    CALL updateWithDiffusionOperator(TWS, t0, nZones)

    ALLOCATE(diffusionOperator(nZones))
    CALL computeDiffusionOperator( &
           TWS%temperature(:,NEW), t1, MeshScale, nZones, diffusionOperator)

    U_s(:,0,iSubNode) = TWS%temperature(:,NEW)

    ! --- Update temperature with reaction operator

    subTimeStep = TWS%timeStep / nReactionSubCycles
    subTime0    = t0
    subTime1    = t0 + subTimeStep

    DO iReactionSubCycle = 1, nReactionSubCycles

      IF(nReactionSubCycles > 1)THEN
        PRINT*
        PRINT*, '  INFO: Subcycling Reaction Update', iReactionSubCycle
        PRINT*, '    Evolving from time ', subTime0
        PRINT*, '               to time ', subTime1, ' in this cycle'
        PRINT*, '    Start Time  = ', t0
        PRINT*, '    Target time = ', t1
        PRINT*
      END IF

      iteration = 0
      converged = .FALSE.
      DO WHILE(.NOT.converged)

        iteration = iteration + 1

        CALL populateLinearizedSystemRHS( &
               TWS%temperature(:,OLD), TWS%temperature(:,NEW), subTime0, &
               subTimeStep, TWS%meshScale, nZones, TWS%rhs, &
               ignoreDiffusionOption = .TRUE.)

        !Add the diffusion operator computed above:
        TWS%rhs(:) = TWS%rhs(:) + subTimestep * diffusionOperator(:)

        CALL populateJacobianMatrix( &
               TWS%temperature(:,NEW), subTimeStep, TWS%meshScale, &
               nZones, TWS%jacobian, ignoreDiffusionOption = .TRUE.)

        !Jacobian is diagonal, linearized system trivial to invert:
        TWS%temperature(:,CHANGE) &
          = TWS%rhs(:) / TWS%jacobian(:,0)
        !Update new temperature:
        TWS%temperature(:,NEW) &
          = TWS%temperature(:,NEW) + TWS%temperature(:,CHANGE)

        WHERE(TWS%temperature(:,NEW) < Zero) &
          TWS%temperature(:,NEW) = 1.0E-16

        WHERE(TWS%temperature(:,NEW) > One) &
          TWS%temperature(:,NEW) = One - 1.0E-16

        CALL checkNonlinearConvergence(converged, iteration, &
               TWS%temperature(:,CHANGE), TWS%temperature(:,NEW), TWS%rhs)

      END DO

      U_s(:,iReactionSubCycle,iSubNode) = TWS%temperature(:,NEW)

      subTime0 = subTime1
      subTime1 = subTime1 + subTimeSTep

      TWS%temperature(:,OLD) = TWS%temperature(:,NEW)

    END DO

    DEALLOCATE(diffusionOperator)

    Tnew(:) = TWS%temperature(:,NEW)
    
    CALL destroyThermalWaveImplicitType(TWS)
    DEALLOCATE(TWS)

  END SUBROUTINE updateThermalWaveMultiImplicit


  SUBROUTINE updateWithDiffusionOperator(TW,t0,nZones)

    TYPE(thermalWaveImplicitType), POINTER :: &
      TW
    REAL(kind=RealKind), INTENT(in) :: &
      t0
    INTEGER, INTENT(in) :: &
      nZones

    CALL populateLinearizedSystemRHS( &
           TW%temperature(:,OLD), TW%temperature(:,NEW), t0, &
           TW%timeStep, TW%meshScale, nZones, TW%rhs, &
           ignoreReactionOption = .TRUE.)

    CALL populateJacobianMatrix( &
           TW%temperature(:,NEW), TW%timeStep, TW%meshScale, &
           nZones, TW%jacobian, ignoreReactionOption = .TRUE.)

    CALL solvePentaDiagonalLinearSystem( &
           TW%jacobian, TW%rhs, nZones, TW%temperature(:,CHANGE))

    !Update new temperature:
    TW%temperature(:,NEW) &
      = TW%temperature(:,NEW) + TW%temperature(:,CHANGE)

  END SUBROUTINE updateWithDiffusionOperator

  SUBROUTINE populateLinearizedSystemRHS(Told, Tnew, t, dt, dx, nZones, rhs, &
               ignoreDiffusionOption, ignoreReactionOption)

    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(in) :: &
      Told, &
      Tnew
    REAL(kind=RealKind), INTENT(in) :: &
      t, &
      dt, &
      dx
    INTEGER, INTENT(in) :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(out) :: &
      rhs
    LOGICAL, INTENT(in), OPTIONAL :: &
      ignoreDiffusionOption, &
      ignoreReactionOption

    LOGICAL :: &
      ignoreDiffusion, &
      ignoreReaction
    REAL(kind=RealKind), DIMENSION(1:nZones) :: &
      diffusionOperator, &
      reactionOperator

    ignoreDiffusion = .FALSE.
    IF(PRESENT(ignoreDiffusionOption)) &
      ignoreDiffusion = ignoreDiffusionOption

    ignoreReaction = .FALSE.
    IF(PRESENT(ignoreReactionOption)) &
      ignoreReaction = ignoreReactionOption

    diffusionOperator(:) = Zero
    reactionOperator(:)  = Zero

    IF(.NOT.ignoreDiffusion) &
      CALL computeDiffusionOperator( &
             Tnew, t+dt, dx, nZones, diffusionOperator)

    IF(.NOT.ignoreReaction) &
      CALL computeReactionOperator( &
             Tnew, nZones, reactionOperator)

    rhs(:) &
      = - ( (Tnew(:)-Told(:)) &
            - dt * (diffusionOperator(:) + reactionOperator(:)) )

  END SUBROUTINE populateLinearizedSystemRHS


  SUBROUTINE computeDiffusionOperator(T,time,dx,nZones,diffusionOperator)

    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(in) :: &
      T
    REAL(kind=RealKind), INTENT(in) :: &
      time, &
      dx
    INTEGER :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(out) :: &
      diffusionOperator

    INTEGER :: &
      iZone
    REAL(kind=RealKind), DIMENSION(-1:nZones+2) :: &
      TMP

    TMP(1:nZones) = T(1:nZones)
    !Dirichlet boundary condition (using analytic solution):  
    TMP(-1) &
      = Half*(One-TANH(-(10.0_RealKind+1.5_RealKind*dx+Two*time/Delta)/Delta))
    TMP( 0) &
      = Half*(One-TANH(-(10.0_RealKind+0.5_RealKind*dx+Two*time/Delta)/Delta))
    TMP(nZones+1) &
      = Half*(One-TANH( (10.0_RealKind+0.5_RealKind*dx-Two*time/Delta)/Delta))
    TMP(nZones+2) &
      = Half*(One-TANH( (10.0_RealKind+1.5_RealKind*dx-Two*time/Delta)/Delta))

    DO iZone = 1, nZones
      diffusionOperator(iZone) &
        = ( - TMP(iZone+2) &
            + 16.0_RealKind * TMP(iZone+1) &
            - 30.0_RealKind * TMP(iZone) &
            + 16.0_RealKind * TMP(iZone-1) &
            - TMP(iZone-2) ) &
          / ( 12.0_RealKind * dx**2 )
    END DO

  END SUBROUTINE computeDiffusionOperator


  SUBROUTINE computeReactionOperator(T, nZones, reactionOperator)

    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(in) :: &
      T
    INTEGER :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(out) :: &
      reactionOperator

    reactionOperator(:) &
      = 8.0_RealKind * T(:)**2 * (One - T(:)) / Delta**2

  END SUBROUTINE computeReactionOperator


  SUBROUTINE populateJacobianMatrix(Tnew, dt, dx, nZones, jacobian, &
               ignoreDiffusionOption, ignoreReactionOption)

    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(in) :: &
      Tnew
    REAL(kind=RealKind), INTENT(in) :: &
      dt, &
      dx
    INTEGER :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones,-2:2), INTENT(out) :: &
      jacobian
    LOGICAL, OPTIONAL :: &
      ignoreDiffusionOption, &
      ignoreReactionOption

    LOGICAL :: &
      ignoreDiffusion, &
      ignoreReaction
    REAL(kind=RealKind), DIMENSION(1:nZones) :: &
      jacobianReaction
    REAL(kind=RealKind), DIMENSION(1:nZones,-2:2) :: &
      jacobianDiffusion

    ignoreDiffusion = .FALSE.
    IF(PRESENT(ignoreDiffusionOption)) &
      ignoreDiffusion = ignoreDiffusionOption

    ignoreReaction = .FALSE.
    IF(PRESENT(ignoreReactionOption)) &
      ignoreReaction = ignoreReactionOption

    jacobian(:,:)          = Zero
    jacobianDiffusion(:,:) = Zero
    jacobianReaction(:)    = Zero

    jacobian(:,0) = One

    IF(.NOT.ignoreDiffusion)THEN

      CALL jacobianFromDiffusionOperator(dx, nZones, jacobianDiffusion)

      jacobian(:,:) = jacobian(:,:) - dt * jacobianDiffusion(:,:)

    END IF

    IF(.NOT.ignoreReaction)THEN

      CALL jacobianFromReactionOperator(Tnew, nZones, jacobianReaction)

      jacobian(:,0) = jacobian(:,0) - dt * jacobianReaction(:)

    END IF

  END SUBROUTINE populateJacobianMatrix


  SUBROUTINE jacobianFromDiffusionOperator(dx, nZones, jacobian)

    REAL(kind=RealKind), INTENT(in) :: &
      dx
    INTEGER :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones,-2:2), INTENT(out) :: &
      jacobian

    INTEGER :: &
      iZone
    REAL(kind=RealKind) :: &
      inverseDenominator

    inverseDenominator = One / ( 12.0_RealKind * dx**2 )

    iZone = 1
    jacobian(iZone, 0) &
      = - 30.0_RealKind * inverseDenominator
    jacobian(iZone,+1) &
      =   16.0_RealKind * inverseDenominator
    jacobian(iZone,+2) &
      = - One * inverseDenominator

    iZone = 2
    jacobian(iZone,-1) &
      =   16.0_RealKind * inverseDenominator
    jacobian(iZone, 0) &
      = - 30.0_RealKind * inverseDenominator
    jacobian(iZone,+1) &
      =   16.0_RealKind * inverseDenominator
    jacobian(iZone,+2) &
      = - One * inverseDenominator

    DO iZone = 3, nZones - 2
      jacobian(iZone,-2) &
        = - One * inverseDenominator
      jacobian(iZone,-1) &
        =   16.0_RealKind * inverseDenominator
      jacobian(iZone, 0) &
        = - 30.0_RealKind * inverseDenominator
      jacobian(iZone,+1) &
        =   16.0_RealKind * inverseDenominator
      jacobian(iZone,+2) &
        = - One * inverseDenominator
    END DO

    iZone = nZones - 1
    jacobian(iZone,-2) &
      = - One * inverseDenominator
    jacobian(iZone,-1) &
      =   16.0_RealKind * inverseDenominator
    jacobian(iZone, 0) &
      = - 30.0_RealKind * inverseDenominator
    jacobian(iZone,+1) &
      =   16.0_RealKind * inverseDenominator

    iZone = nZones
    jacobian(iZone,-2) &
      = - One * inverseDenominator
    jacobian(iZone,-1) &
      =   16.0_RealKind * inverseDenominator
    jacobian(iZone, 0) &
      = - 30.0_RealKind * inverseDenominator

  END SUBROUTINE jacobianFromDiffusionOperator


  SUBROUTINE jacobianFromReactionOperator(T, nZones, jacobian)

    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(in) :: &
      T
    INTEGER :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(out) :: &
      jacobian

    jacobian(:) &
      = 16.0_RealKind * T(:) * ( One - 1.5_RealKind * T(:) ) / Delta**2

  END SUBROUTINE jacobianFromReactionOperator


  SUBROUTINE populateRhsAndJacobianSmAnltc(Told,Tnew,dt,nZones,rhs,jacobian)

    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(in) :: &
      Told, &
      Tnew
    REAL(kind=RealKind), INTENT(in) :: &
      dt
    INTEGER, INTENT(in) :: &
      nZones
    REAL(kind=RealKind), DIMENSION(1:nZones), INTENT(out) :: &
      rhs
    REAL(kind=RealKind), DIMENSION(1:nZones,-2:2), INTENT(out) :: &
      jacobian

    rhs(:) &
      = - ( LOG(ABS((Tnew(:)*(One-Told(:)))/(Told(:)*(One-Tnew(:))))) &
            + (Tnew(:)-Told(:))/(Tnew(:)*Told(:)) &
            - 8.0_RealKind*dt/Delta**2 )

    jacobian(:,:) = Zero
    jacobian(:,0) = One/(Tnew(:)**2*(One-Tnew(:)))

  END SUBROUTINE populateRhsAndJacobianSmAnltc


  SUBROUTINE createThermalWaveImplicitType(TW, nZones)

    TYPE(thermalWaveImplicitType), POINTER :: &
      TW
    INTEGER, INTENT(in) :: &
      nZones

    TW%timeStep  = Zero
    TW%meshScale = Zero

    ALLOCATE(TW%rhs(nZones))
    TW%rhs(:) = Zero

    ALLOCATE(TW%jacobian(nZones,-2:2))
    TW%jacobian(:,:) = Zero

    ALLOCATE(TW%temperature(nZones,3))
    TW%temperature(:,:) = Zero

  END SUBROUTINE createThermalWaveImplicitType


  SUBROUTINE destroyThermalWaveImplicitType(TW)

    TYPE(thermalWaveImplicitType), POINTER :: &
      TW

    DEALLOCATE(TW%rhs)
    DEALLOCATE(TW%jacobian)
    DEALLOCATE(TW%temperature)

  END SUBROUTINE destroyThermalWaveImplicitType


  SUBROUTINE solvePentaDiagonalLinearSystem(Matrix, RHS, N, x)

    !
    ! --- Subroutine solves pentadiagonal linear system with lapack (DGBSV)
    !

    REAL(kind=RealKind), DIMENSION(1:N,-2:2), INTENT(in) :: &
      Matrix
    REAL(kind=RealKind), DIMENSION(1:N), INTENT(in) :: &
      RHS
    INTEGER :: &
      N
    REAL(kind=RealKind), DIMENSION(1:N), INTENT(out) :: &
      x

    INTEGER, PARAMETER :: &
      KL   = 2, &
      KU   = 2, &
      NRHS = 1
    INTEGER :: &
      LDAB, &
      LDB, &
      i, j, ii, &
      INFO
    INTEGER, DIMENSION(:), ALLOCATABLE :: &
      IPIV
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      B
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      AB

    LDAB = 2*KU+KL+1
    LDB  = N

    ALLOCATE(IPIV(N))
    ALLOCATE(AB(LDAB,N))
    ALLOCATE(B(LDB))

    DO j = 1, N
      DO i = MAX(1,j-KU), MIN(N,j+KL)
        !Insert matrix elements in AB:
        ii = KL+KU+1+i-j
        AB(ii,j) = Matrix(j,i-j)
      END DO
      !Right-hand side:
      B(j) = RHS(j)
    END DO

    CALL DGBSV(N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)

    x(:) = B(:)

    DEALLOCATE(IPIV)
    DEALLOCATE(AB)
    DEALLOCATE(B)

  END SUBROUTINE solvePentaDiagonalLinearSystem


  SUBROUTINE checkNonlinearConvergence(converged,iteration,dT,T,rhs)

    LOGICAL, INTENT(out) :: &
      converged
    INTEGER, INTENT(in) :: &
      iteration
    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      dT, T, rhs

    INTEGER :: &
      nZones, &
      iZone
    INTEGER, PARAMETER :: &
      MaxIterations = 150
    REAL(kind=RealKind) :: &
      relativeChange, &
      relativeChangeNorm, &
      maxRelativeChange, &
      nonlinearResidualNorm, &
      relativeNonlinearResidualNorm
    REAL(kind=RealKind), SAVE :: &
      initialNonlinearResidualNorm
    REAL(kind=RealKind), PARAMETER :: &
      RelativeChangeNormTolerance   = 1.0E-10_RealKind, &
      MaxRelativeChangeTolerance    = 1.0E-09_RealKind, &
      RelativeResidualNormTolerance = 1.0E-02_RealKind

    nZones = SIZE(dT)

    relativeChangeNorm = Zero
    maxRelativeChange  = Zero
    DO iZone = 1, nZones
      relativeChange &
        = SQRT( ( dT(iZone) / T(iZone) )**2 )
      relativeChangeNorm &
        = relativeChangeNorm + relativeChange**2
      maxRelativeChange &
        = MAX(relativeChange, maxRelativeChange)
    END DO
    relativeChangeNorm &
      = SQRT( relativeChangeNorm )

    IF(iteration == 1)THEN
      initialNonlinearResidualNorm &
        = SQRT( DOT_PRODUCT(rhs(:), rhs(:)) )
      nonlinearResidualNorm &
        = initialNonlinearResidualNorm
    ELSE
      nonlinearResidualNorm &
        = SQRT( DOT_PRODUCT(rhs(:), rhs(:)) )
    END IF
    relativeNonlinearResidualNorm &
      = nonlinearResidualNorm &
          / MAX(initialNonlinearResidualNorm, 1.0E-16_RealKind)

    PRINT*
    PRINT*, '    INFO: Iteration ', iteration
    PRINT*, '                  Relative Change Norm =', &
      relativeChangeNorm
    PRINT*, '                   Max Relative Change =', &
      maxRelativeChange
    PRINT*, '               Nonlinear Residual Norm =', &
      nonlinearResidualNorm
    PRINT*, '      Relative Nonlinear Residual Norm =', &
      relativeNonlinearResidualNorm

    IF(iteration > MaxIterations)THEN
      PRINT*, '    INFO: Iteration Exceeds MaxIterations.  Exiting.'
      PRINT*, '          Iteration =', iteration
      PRINT*, '      MaxIterations =', MaxIterations
      converged = .TRUE.
    ELSE
      IF(relativeChangeNorm <= RelativeChangeNormTolerance &
         .AND. maxRelativeChange <= MaxRelativeChangeTolerance) &
      THEN
        converged = .TRUE.
      END IF
    END IF

    IF(converged)THEN
      PRINT*
      PRINT*, '  INFO: Converged after ', iteration, ' iterations'
      PRINT*
      PRINT*, '                Relative Change Norm =', &
        relativeChangeNorm
      PRINT*, '                 Max Relative Change =', &
        maxRelativeChange
      PRINT*, '     Initial Nonlinear Residual Norm =', &
        initialNonlinearResidualNorm
      PRINT*, '             Nonlinear Residual Norm =', &
        nonlinearResidualNorm
      PRINT*, '    Relative Nonlinear Residual Norm =', &
        relativeNonlinearResidualNorm
      PRINT*
    END IF

  END SUBROUTINE checkNonlinearConvergence


  SUBROUTINE evaluateOdeRhsThermalWave(t,U,RHS)

    REAL(kind=RealKind), INTENT(in) :: &
      t
    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      U
    REAL(kind=RealKind), DIMENSION(:), INTENT(out) :: &
      RHS

    INTEGER :: &
      nZones
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      diffusionOperator, &
      reactionOperator

    nZones = SIZE(U)

    ALLOCATE(diffusionOperator(nZones))
    ALLOCATE(reactionOperator(nZones))

    CALL computeDiffusionOperator( &
           U, t, MeshScale, nZones, diffusionOperator)

    CALL computeReactionOperator( &
           U, nZones, reactionOperator)

    RHS(:) &
      = diffusionOperator(:) + reactionOperator(:)

    DEALLOCATE(diffusionOperator)
    DEALLOCATE(reactionOperator)

  END SUBROUTINE evaluateOdeRhsThermalWave


  SUBROUTINE computeCorrectionThermalWave(t,U,Residual,Correction)

    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      t
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(in) :: &
      U, &
      Residual
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(out) :: &
      Correction

    INTEGER :: &
      nZones, &
      iSubNode
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      RHS
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      Matrix

    nZones = SIZE(U,DIM=1)

    ALLOCATE(RHS(nZones))
    ALLOCATE(Matrix(nZones,-2:2))

    Correction(:,:) = Zero
    DO iSubNode = 2, SIZE(t)

      RHS(:) = Correction(:,iSubNode-1) &
                 + (Residual(:,iSubNode)-Residual(:,iSubNode-1))

      CALL populateJacobianMatrix( &
             U(:,iSubNode), t(iSubNode)-t(iSubNode-1), MeshScale, nZones, &
             Matrix)

      CALL solvePentaDiagonalLinearSystem( &
             Matrix, RHS, nZones, Correction(:,iSubNode))

    END DO

    DEALLOCATE(RHS)
    DEALLOCATE(Matrix)

  END SUBROUTINE computeCorrectionThermalWave


  SUBROUTINE computeCorrectionThermalWaveSplit(t,U,Residual,Correction, &
               reverseOption)

    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      t
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(in) :: &
      U, &
      Residual
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(out) :: &
      Correction
    LOGICAL, INTENT(in), OPTIONAL :: &
      reverseOption

    LOGICAL :: &
      reverse
    INTEGER :: &
      nZones, &
      iSubNode, &
      iSubCycle
    REAL(kind=RealKind) :: &
      dt
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      RHS
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      Matrix

    reverse = .FALSE.
    IF(PRESENT(reverseOption)) &
      reverse = reverseOption

    nZones = SIZE(U,DIM=1)

    ALLOCATE(RHS(nZones))
    ALLOCATE(Matrix(nZones,-2:2))

    Correction(:,:) = Zero
    DO iSubNode = 2, SIZE(t)

      IF(.NOT.reverse)THEN

        ! --- Compute correction from reaction operator

        RHS(:) = Correction(:,iSubNode-1) &
                   + (Residual(:,iSubNode)-Residual(:,iSubNode-1))

        CALL populateJacobianMatrix( &
               U(:,iSubNode), t(iSubNode)-t(iSubNode-1), MeshScale, nZones, &
               Matrix, ignoreDiffusionOption = .TRUE.)

        Correction(:,iSubNode) = RHS(:) / Matrix(:,0)

      END IF

      ! --- Update correction from diffusion operator

      IF(.NOT.reverse)THEN
        RHS(:) = Correction(:,iSubNode)
      ELSE
        RHS(:) = Correction(:,iSubNode-1) &
                   + (Residual(:,iSubNode)-Residual(:,iSubNode-1))
      END IF

      CALL populateJacobianMatrix( &
             U(:,iSubNode), t(iSubNode)-t(iSubNode-1), MeshScale, nZones, &
             Matrix, ignoreReactionOption = .TRUE.)

      CALL solvePentaDiagonalLinearSystem( &
             Matrix, RHS, nZones, Correction(:,iSubNode))

      IF(reverse)THEN

        dt = (t(iSubNode)-t(iSubNode-1)) / REAL(nReactionSubCycles)
        DO iSubCycle = 1, nReactionSubCycles

          RHS(:) = Correction(:,iSubNode)

          CALL populateJacobianMatrix( &
                 U(:,iSubNode), dt, MeshScale, nZones, &
                 Matrix, ignoreDiffusionOption = .TRUE.)

          Correction(:,iSubNode) = RHS(:) / Matrix(:,0)

        END DO

      END IF

    END DO

    DEALLOCATE(RHS)
    DEALLOCATE(Matrix)

  END SUBROUTINE computeCorrectionThermalWaveSplit


  SUBROUTINE computeCorrectionThermalWaveMI(t,U,F,Residual,Correction)

    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      t
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(in) :: &
      U, &
      F, &
      Residual
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(out) :: &
      Correction

    INTEGER :: &
      nZones, &
      iSubNode, &
      iSubSubNode
    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      RHS, &
      F_D0, F_D1, &
      t_s
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      Matrix, &
      Correction_s

    nZones = SIZE(U,DIM=1)

    ALLOCATE(RHS(nZones))
    ALLOCATE(F_D0(nZones))
    ALLOCATE(F_D1(nZones))
    ALLOCATE(Matrix(nZones,-2:2))
    ALLOCATE(t_s(0:nReactionSubCycles))
    ALLOCATE(Correction_s(1:nZones,0:nReactionSubCycles))

    Correction(:,:) = Zero
    DO iSubNode = 2, SIZE(t)

      RHS(:) = Correction(:,iSubNode-1) &
                 + (Residual(:,iSubNode)-Residual(:,iSubNode-1))

      CALL populateJacobianMatrix( &
             U(:,iSubNode), t(iSubNode)-t(iSubNode-1), MeshScale, nZones, &
             Matrix, ignoreReactionOption = .TRUE.)

      CALL solvePentaDiagonalLinearSystem( &
             Matrix, RHS, nZones, Correction(:,iSubNode))

      CALL computeDiffusionOperator( &
             U(:,iSubNode), t(iSubNode), &
             MeshScale, nZones, F_D0(:))
      CALL computeDiffusionOperator( &
             U(:,iSubNode)+Correction(:,iSubNode), t(iSubNode), &
             MeshScale, nZones, F_D1(:))

      t_s(0) = t(iSubNode-1)
      DO iSubSubNode = 1, nReactionSubCycles
        t_s(iSubSubNode) &
          = t_s(iSubSubNode-1) &
              + (t(iSubNode)-t(iSubNode-1))/REAL(nReactionSubCycles)
      END DO

      Correction_s(:,0) = Correction(:,iSubNode-1)
      DO iSubSubNode = 1, nReactionSubCycles

        RHS(:) = Correction_s(:,iSubSubNode-1) &
                   + (t_s(iSubSubNode)-t_s(iSubSubNode-1)) &
                     * (F_D1(:)-F_D0(:)) &
                   + polynomialIntegral( &
                       t(:), F(:,:), t_s(iSubSubNode-1), t_s(iSubSubNode)) &
                   - (U_s(:,iSubSubNode,iSubNode) &
                      - U_s(:,iSubSubNode-1,iSubNode))

        CALL populateJacobianMatrix( &
               U_s(:,iSubSubNode,iSubNode), &
               (t_s(iSubSubNode)-t_s(iSubSubNode-1)), MeshScale, nZones, &
               Matrix, ignoreDiffusionOption = .TRUE.)

        Correction_s(:,iSubSubNode) = RHS(:) / Matrix(:,0)

      END DO

      Correction(:,iSubNode) = Correction_s(:,nReactionSubCycles)

    END DO

    DEALLOCATE(RHS)
    DEALLOCATE(F_D0)
    DEALLOCATE(F_D1)
    DEALLOCATE(Matrix)
    DEALLOCATE(t_s)
    DEALLOCATE(Correction_s)

  END SUBROUTINE computeCorrectionThermalWaveMI


END MODULE thermalWaveSolutionModule

PROGRAM thermalWave

  USE realKindModule, ONLY: &
    RealKind, &
    Zero, &
    Half, &
    One, &
    Two
  USE meshModule, ONLY: &
    Center
  USE fieldsModule, ONLY: &
    Temperature, &
    OldState
  USE fieldsInputOutputModule, ONLY: &
    writeFields1D
  USE universeModule, ONLY: &
    universeType
  USE universeCreationModule, ONLY: &
    createUniverse, &
    destroyUniverse
  USE sdcIntegrationModule, ONLY: &
    stepSDC
  USE thermalWaveInitializationModule, ONLY: &
    initializeThermalWave
  USE thermalWaveModule, ONLY: &
    Delta, &
    nReactionSubCycles, &
    U_s

  IMPLICIT NONE

  TYPE(universeType), POINTER :: &
    U => NULL()

  LOGICAL :: &
    writeData, &
    done
  INTEGER :: &
    iCycle, &
    nZones, &
    iZone
  INTEGER, PARAMETER :: &
    nSubNodes    = 2, &
    nCorrections = 0
  REAL(kind=RealKind) :: &
    writeTime, &
    timeStep, &
    L1ErrorNorm
  REAL(kind=RealKind), PARAMETER :: &
    StartTime         = Zero, &
    EndTime           = 1.000E-1_RealKind, &
    WriteTimeInterwal = 1.280E-0_RealKind
  REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
    AnalyticSolution

  CALL createUniverse(U, name = 'thermalWave', &
         innerBoundaries = (/-10.0_RealKind/), &
         outerBoundaries = (/+10.0_RealKind/), &
         nZones = (/8192,1,1/))

  ALLOCATE(U_s(1:SIZE(U%fields%data(:,1,1,Temperature,OldState)), &
               0:nReactionSubCycles, 1:nSubNodes))

  PRINT*
  PRINT*, 'INFO:  thermalWave.f90'
  PRINT*, '  StartTime    = ', StartTime
  PRINT*, '  EndTime      = ', EndTime
  PRINT*, '  nSubNodes    = ', nSubNodes
  PRINT*, '  nCorrections = ', nCorrections

  CALL initializeThermalWave(U)

  U%time    = StartTime
  writeTime = WriteTimeInterwal
  writeData = .FALSE.

  done   = .FALSE.
  iCycle = 0
  DO WHILE(.NOT.done)

    iCycle = iCycle + 1

    PRINT*
    PRINT*, '  cycle    =', iCycle

    timeStep = MIN(1.0E-1_RealKind, EndTime-U%time)

    IF(writeTime < U%time + timeStep + 1.0E-15_RealKind)THEN
      PRINT*
      PRINT*, '  INFO: Reducing time step to hit target write time interwal'
      timeStep  = writeTime - U%time
      writeData = .TRUE.
    END IF

    PRINT*
    PRINT*, '  time     =', U%time
    PRINT*, '  timeStep =', timeStep

    CALL stepSDC(U%fields%data(:, 1, 1, Temperature, OldState), &
                 timeBegin = U%time, timeEnd = U%time+timeStep, &
                 nSteps = 1, nSubNodes = nSubNodes, &
                 nCorrections = nCorrections, name = U%name)

    U%time = U%time + timeStep

    IF(writeData)THEN
      CALL writeFields1D(U)
      writeTime = U%time + WriteTimeInterwal
      writeData = .FALSE.
    END IF

    IF(ABS(U%time-EndTime)/EndTime <= 1.0E-15_RealKind) &
      done = .TRUE.

  END DO

  CALL writeFields1D(U)

  nZones = SIZE(U%fields%data(:,1,1,Temperature,OldState))

  ALLOCATE(AnalyticSolution(nZones))

  AnalyticSolution(:) &
    = Half * ( One - TANH( (U%mesh%positions(1)%data(:,Center) &
                            - Two*U%time/Delta) / Delta ) )

  L1ErrorNorm = Zero
  DO iZone = 1, nZones
    L1ErrorNorm &
      = L1ErrorNorm &
          + ABS(AnalyticSolution(iZone) &
                - U%fields%data(iZone,1,1,Temperature,OldState))
  END DO
  L1ErrorNorm = L1ErrorNorm / REAL(nZones)

  DEALLOCATE(AnalyticSolution)

  PRINT*
  PRINT*, 'INFO: Finished thermalWave.f90'
  PRINT*
  PRINT*, '  time     =', U%time
  PRINT*, '  L1-error =', L1ErrorNorm
  PRINT*

  DEALLOCATE(U_s)

  CALL destroyUniverse(U)

END PROGRAM thermalWave

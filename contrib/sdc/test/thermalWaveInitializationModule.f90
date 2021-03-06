MODULE thermalWaveInitializationModule

  USE realKindModule, ONLY: &
    RealKind, &
    Zero, &
    Half, &
    One
  USE meshModule, ONLY: &
    Center
  USE fieldsModule, ONLY: &
    Temperature, &
    OldState
  USE fieldsInputOutputModule, ONLY: &
    writeFields1D
  USE universeModule, ONLY: &
    universeType
  USE thermalWaveModule, ONLY: &
    MeshScale, &
    Delta

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
    initializeThermalWave

CONTAINS

  SUBROUTINE initializeThermalWave(U)

    TYPE(universeType), POINTER :: &
      U

    INTEGER :: &
      InitialState, &
      iZone, jZone, kZone
    REAL(kind=RealKind), DIMENSION(U%nZones(1)) :: &
      x

    !Initialize MeshScale defined in thermalWaveModule
    MeshScale = U%mesh%meshScale(1)

    InitialState = OldState
    x = U%mesh%positions(1)%data(:,Center)

    DO kZone = 1, U%nZones(3)
      DO jZone = 1, U%nZones(2)
        DO iZone = 1, U%nZones(1)
          U%fields%data(iZone,jZone,kZone,Temperature,InitialState) &
            = Half * (One - TANH(x(iZone)/Delta))
        END DO
      END DO
    END DO

    CALL writeFields1D(U,stateOption=InitialState)

  END SUBROUTINE initializeThermalWave

END MODULE thermalWaveInitializationModule

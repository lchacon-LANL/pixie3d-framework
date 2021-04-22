MODULE thermalWaveModule

  USE realKindModule, ONLY: &
    RealKind

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: &
    nReactionSubCycles = 1
  REAL(kind=RealKind), PUBLIC :: &
    MeshScale
  REAL(kind=RealKind), PUBLIC, PARAMETER :: &
    Delta = 1.00E-0_RealKind

  REAL(kind=RealKind), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: &
    U_s ! Used to store subcyele updates in Multi Implicit Scheme

END MODULE thermalWaveModule

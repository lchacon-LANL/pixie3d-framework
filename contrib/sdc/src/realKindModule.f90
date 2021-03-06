MODULE realKindModule

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PARAMETER :: &
    RealKind = KIND(1.d0)

  REAL(kind=RealKind), PUBLIC, PARAMETER :: &
    Zero  = 0.0_RealKind, &
    Half  = 0.5_RealKind, &
    One   = 1.0_RealKind, &
    Two   = 2.0_RealKind, &
    Three = 3.0_RealKind, &
    Four  = 4.0_RealKind

END MODULE realKindModule

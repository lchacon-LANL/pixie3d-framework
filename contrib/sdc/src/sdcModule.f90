MODULE sdcModule

  USE realKindModule, ONLY: &
    RealKind, &
    Zero, &
    Half, &
    Two, &
    Three, &
    Four

  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: sdcType

    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
      t
    REAL(kind=RealKind), DIMENSION(:,:), ALLOCATABLE :: &
      U, M, Residual, Correction!, OdeRhs

  END TYPE sdcType

  PUBLIC :: &
    createSDC, &
    destroySDC, &
    setSubNodes, &
    setIntegrationMatrix
!!!    polynomialIntegral

CONTAINS


  SUBROUTINE createSDC(U,nSubNodes,SDC)

    INTEGER, INTENT(in) :: &
      nSubNodes
    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      U
    TYPE(sdcType), POINTER, INTENT(inout) :: &
      SDC

    ALLOCATE(SDC)
    ALLOCATE(SDC%t(nSubNodes))
    ALLOCATE(SDC%U (SIZE(U), nSubNodes))
    ALLOCATE(SDC%M(nSubNodes, nSubNodes))
!!    ALLOCATE(SDC%OdeRhs(SIZE(U), nSubNodes))
    ALLOCATE(SDC%Residual(SIZE(U), nSubNodes))
    ALLOCATE(SDC%Correction(SIZE(U), nSubNodes))

    SDC%t(:)            = Zero
    SDC%U(:,:)          = Zero
    SDC%M(:,:)          = Zero
!!    SDC%OdeRhs(:,:)     = Zero
    SDC%Residual(:,:)   = Zero
    SDC%Correction(:,:) = Zero

  END SUBROUTINE createSDC


  SUBROUTINE destroySDC(SDC)

    TYPE(sdcType), POINTER, INTENT(inout) :: &
      SDC

    DEALLOCATE(SDC%Correction)
    DEALLOCATE(SDC%Residual)
!!    DEALLOCATE(SDC%OdeRhs)
    DEALLOCATE(SDC%M)
    DEALLOCATE(SDC%U)
    DEALLOCATE(SDC%t)
    DEALLOCATE(SDC)

  END SUBROUTINE destroySDC


  SUBROUTINE setSubNodes(t0,t1,t)

    REAL(kind=RealKind), INTENT(in) :: &
      t0, t1
    REAL(kind=RealKind), DIMENSION(:), INTENT(inout) :: &
      t

    INTEGER :: &
      i

    !Uniform sub-node spacing:
    DO i = 1, SIZE(t)
       t(i) = t0 + REAL(i-1)*(t1-t0)/REAL(SIZE(t)-1)
    END DO

  END SUBROUTINE setSubNodes

  SUBROUTINE setIntegrationMatrix(t, M)

    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
      t
    REAL(kind=RealKind), DIMENSION(:,:), INTENT(out) :: &
      M

    REAL(kind=RealKind) :: &
      t0, t1, t2, t3

    SELECT CASE (SIZE(t))

      CASE(2)
        t0 = t(1)-t(1)
        t1 = t(2)-t(1)
        ! Row 1:
        M(1,:) = Zero
        ! Row 2:
        M(2,1) = ( Half*(t1-t0)**2-t1*(t1-t0) )/( t0-t1 )
        M(2,2) = ( Half*(t1-t0)**2-t0*(t1-t0) )/( t1-t0 )

      CASE(3)
        t0 = t(1)-t(1)
        t1 = t(2)-t(1)
        t2 = t(3)-t(1)
        ! Row 1:
        M(1,:) = Zero
        ! Row 2:
        M(2,1) = ( (t1-t0)**3 / Three - (t1+t2)*(t1-t0)**2 / Two &
                   + t1*t2*(t1-t0) ) / ( (t0-t1)*(t0-t2) )
        M(2,2) = ( (t1-t0)**3 / Three - (t0+t2)*(t1-t0)**2 / Two &
                   + t0*t2*(t1-t0) ) / ( (t1-t0)*(t1-t2) )
        M(2,3) = ( (t1-t0)**3 / Three - (t0+t1)*(t1-t0)**2 / Two &
                   + t0*t1*(t1-t0) ) / ( (t2-t0)*(t2-t1) )
        ! Row 3:
        M(3,1) = ( (t2-t0)**3 / Three - (t1+t2)*(t2-t0)**2 / Two &
                   + t1*t2*(t2-t0) ) / ( (t0-t1)*(t0-t2) )
        M(3,2) = ( (t2-t0)**3 / Three - (t0+t2)*(t2-t0)**2 / Two &
                   + t0*t2*(t2-t0) ) / ( (t1-t0)*(t1-t2) )
        M(3,3) = ( (t2-t0)**3 / Three - (t0+t1)*(t2-t0)**2 / Two &
                   + t0*t1*(t2-t0) ) / ( (t2-t0)*(t2-t1) )

      CASE(4)
        t0 = t(1)-t(1)
        t1 = t(2)-t(1)
        t2 = t(3)-t(1)
        t3 = t(4)-t(1)
        ! Row 1:
        M(1,:) = Zero
        ! Row 2:
        M(2,1) = ( (t1-t0)**4 / Four - (t1+t2+t3)*(t1-t0)**3 / Three &
                   + (t1*t2+t1*t3+t2*t3)*(t1-t0)**2 / Two &
                   - t1*t2*t3*(t1-t0) ) / ( (t0-t1)*(t0-t2)*(t0-t3) )
        M(2,2) = ( (t1-t0)**4 / Four - (t0+t2+t3)*(t1-t0)**3 / Three &
                   + (t0*t2+t0*t3+t2*t3)*(t1-t0)**2 / Two &
                   - t0*t2*t3*(t1-t0) ) / ( (t1-t0)*(t1-t2)*(t1-t3) )
        M(2,3) = ( (t1-t0)**4 / Four - (t0+t1+t3)*(t1-t0)**3 / Three &
                   + (t0*t1+t0*t3+t1*t3)*(t1-t0)**2 / Two &
                   - t0*t1*t3*(t1-t0) ) / ( (t2-t0)*(t2-t1)*(t2-t3) )
        M(2,4) = ( (t1-t0)**4 / Four - (t0+t1+t2)*(t1-t0)**3 / Three &
                   + (t0*t1+t0*t2+t1*t2)*(t1-t0)**2 / Two &
                   - t0*t1*t2*(t1-t0) ) / ( (t3-t0)*(t3-t1)*(t3-t2) )
        ! Row 3:
        M(3,1) = ( (t2-t0)**4 / Four - (t1+t2+t3)*(t2-t0)**3 / Three &
                   + (t1*t2+t1*t3+t2*t3)*(t2-t0)**2 / Two &
                   - t1*t2*t3*(t2-t0) ) / ( (t0-t1)*(t0-t2)*(t0-t3) )
        M(3,2) = ( (t2-t0)**4 / Four - (t0+t2+t3)*(t2-t0)**3 / Three &
                   + (t0*t2+t0*t3+t2*t3)*(t2-t0)**2 / Two &
                   - t0*t2*t3*(t2-t0) ) / ( (t1-t0)*(t1-t2)*(t1-t3) )
        M(3,3) = ( (t2-t0)**4 / Four - (t0+t1+t3)*(t2-t0)**3 / Three &
                   + (t0*t1+t0*t3+t1*t3)*(t2-t0)**2 / Two &
                   - t0*t1*t3*(t2-t0) ) / ( (t2-t0)*(t2-t1)*(t2-t3) )
        M(3,4) = ( (t2-t0)**4 / Four - (t0+t1+t2)*(t2-t0)**3 / Three &
                   + (t0*t1+t0*t2+t1*t2)*(t2-t0)**2 / Two &
                   - t0*t1*t2*(t2-t0) ) / ( (t3-t0)*(t3-t1)*(t3-t2) )
        ! Row 4:
        M(4,1) = ( (t3-t0)**4 / Four - (t1+t2+t3)*(t3-t0)**3 / Three &
                   + (t1*t2+t1*t3+t2*t3)*(t3-t0)**2 / Two &
                   - t1*t2*t3*(t3-t0) ) / ( (t0-t1)*(t0-t2)*(t0-t3) )
        M(4,2) = ( (t3-t0)**4 / Four - (t0+t2+t3)*(t3-t0)**3 / Three &
                   + (t0*t2+t0*t3+t2*t3)*(t3-t0)**2 / Two &
                   - t0*t2*t3*(t3-t0) ) / ((t1-t0)*(t1-t2)*(t1-t3))
        M(4,3) = ( (t3-t0)**4 / Four - (t0+t1+t3)*(t3-t0)**3 / Three &
                   + (t0*t1+t0*t3+t1*t3)*(t3-t0)**2 / Two &
                   - t0*t1*t3*(t3-t0) ) / ( (t2-t0)*(t2-t1)*(t2-t3) )
        M(4,4) = ( (t3-t0)**4 / Four - (t0+t1+t2)*(t3-t0)**3 / Three &
                   + (t0*t1+t0*t2+t1*t2)*(t3-t0)**2 / Two &
                   - t0*t1*t2*(t3-t0) ) / ( (t3-t0)*(t3-t1)*(t3-t2) )

      CASE DEFAULT

!!        M(:,:) = Zero
        PRINT*
        PRINT*, '  INFO: Integration matrix only implemented for n=2,3,4'
        PRINT*, '    n =', SIZE(t)
        STOP

    END SELECT

  END SUBROUTINE setIntegrationMatrix


!!$  FUNCTION polynomialIntegral(t,f,t1,t2) RESULT(integral)
!!$
!!$    !  Function integrates polynomial f(t) from t1 to t2.  
!!$    !    Polynomial constructed with t(1:n) and f(1:n).  
!!$    !    Note: t1,t2 in [t(1),t(n)].  
!!$
!!$    REAL(kind=RealKind), DIMENSION(:), INTENT(in) :: &
!!$      t
!!$    REAL(kind=RealKind), DIMENSION(:,:), INTENT(in) :: &
!!$      f
!!$    REAL(kind=RealKind), INTENT(in) :: &
!!$      t1, t2
!!$    REAL(kind=RealKind), DIMENSION(:) :: &
!!$      integral
!!$
!!$    INTEGER :: &
!!$      n, i
!!$    REAL(kind=RealKind), DIMENSION(:), ALLOCATABLE :: &
!!$      cardinalIntegral
!!$
!!$    n = SIZE(t)
!!$
!!$    ALLOCATE(cardinalIntegral(n))
!!$
!!$    SELECT CASE (n)
!!$
!!$      CASE(2)
!!$
!!$        cardinalIntegral(1) &
!!$          = ((t(2)-Half*t2)*t2-(t(2)-Half*t1)*t1)/(t(2)-t(1))
!!$        cardinalIntegral(2) &
!!$          = ((Half*t2-t(1))*t2-(Half*t1-t(1))*t1)/(t(2)-t(1))
!!$
!!$      CASE(3)
!!$
!!$      CASE(4)
!!$
!!$      CASE DEFAULT
!!$
!!$        PRINT*
!!$        PRINT*, '  INFO: polynomialIntegral only implemented for n=2,3,4'
!!$        PRINT*, '    n =', n
!!$        STOP
!!$
!!$    END SELECT
!!$
!!$    integral(:) = Zero
!!$    DO i = 1, n
!!$      integral(:) &
!!$        = integral(:) + cardinalIntegral(i) * f(:,i)
!!$    END DO
!!$
!!$    DEALLOCATE(cardinalIntegral)
!!$
!!$  END FUNCTION polynomialIntegral


END MODULE sdcModule

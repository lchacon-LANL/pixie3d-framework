c module Block_TriDiagonal_Solver
c ##############################################################################
      Module Block_TriDiagonal_Solver

      real(8), pointer :: B_a(:,:)=>null(),B_b(:,:)=>null()
     .                   ,B_c(:,:)=>null()
      real(8), pointer :: invbm(:,:)=>null(),lm(:,:)=>null()
     .       ,dm(:,:)=>null(),rm(:,:)=>null()
      real(8), pointer, private :: rr(:)=>null()
      real(8), pointer :: Pm(:,:)
      integer :: CODE
!     for lapack
      integer :: NRHS,LDA,LDB,INFO
      integer, pointer :: IPIV(:)

      CONTAINS  

c     BTDS
c     ##########################################################################
      SUBROUTINE BTDS(nx,neq,nr,aa,bb,cc,r,sol)

      implicit none
      integer, intent(in) :: nx, neq, nr
      real(8), intent(in), target :: aa(nx*neq,neq)
     .                              ,bb(nx*neq,neq)
     .                              ,cc(nx*neq,neq)
     .                              ,r (nx*neq)
      real(8), intent(out) :: sol(nx*neq)

c     Begin program

      B_a => aa(neq+1:(nx-1)*neq,1:neq)
      B_b => bb(neq+1:(nx-1)*neq,1:neq)
      B_c => cc(neq+1:(nx-1)*neq,1:neq)

      allocate(IPIV(2*neq))
      allocate(Pm(2*neq,2*neq))

      allocate(rr(2*(neq)))
      allocate(invbm(1:(nx-2)*neq,1:nr)) !1st col:B^-1*rm(:,1); 2nd,3rd col: B^-1*U

c     Fill

      rr(1:neq) = r(1:neq)
      rr(neq+1:2*neq) = r((nx-1)*neq+1:nx*neq)

      allocate(rm(1:(nx-2)*neq,1:nr)) !U
      rm(1:(nx-2)*neq,1) = r(neq+1:(nx-1)*neq)
      rm(1:neq,2:1+neq)  = aa(neq+1:2*neq,1:neq)
      rm(neq+1:(nx-2)*neq,2:1+neq) = 0d0
      rm(1:(nx-3)*neq,1+neq+1:1+neq*2) = 0d0
      rm((nx-3)*neq+1:(nx-2)*neq,1+neq+1:1+neq*2) = 
     .     cc((nx-2)*neq+1:(nx-1)*neq,1:neq)

      allocate(lm(2*neq,(nx-2)*neq)) !L
      lm(1:neq,1:neq) = cc(1:neq,1:neq)
      lm(1:neq,neq+1:(nx-2)*neq)   = 0d0
      lm(neq+1:2*neq,1:(nx-3)*neq) = 0d0
      lm(neq+1:2*neq,(nx-3)*neq+1:(nx-2)*neq) = 
     .     aa((nx-1)*neq+1:nx*neq,1:neq)
      
      allocate(dm(2*neq,2*neq)) !D
      dm(1:neq,1:neq) = bb(1:neq,1:neq)
      dm(1:neq,neq+1:2*neq) = aa(1:neq,1:neq)
      dm(neq+1:2*neq,1:neq) = cc((nx-1)*neq+1:nx*neq,1:neq)
      dm(neq+1:2*neq,neq+1:2*neq) = bb((nx-1)*neq+1:nx*neq,1:neq)

c     Solve

      call BLOCK_TRIDIAG(B_a,B_b,B_c,rm,invbm,nx-2,neq,nr,CODE)

c      dm = dm-matmul(lm,invbm(:,2:5))
      dm = dm-matmul(lm,invbm(:,2:2*neq+1)) !careful   
      
!     !$  call Invert_matrix(dm,Pm,4,CODE)
!     !$  rr = matmul(Pm, rr-matmul(lm,invbm(:,1)))

!     or call lapack
      rr = rr-matmul(lm,invbm(:,1))
      
      NRHS = 1
      LDA  = 2*neq
      LDB  = LDA
      CALL DGESV( 2*neq, NRHS, dm, LDA, IPIV, rr, LDB, INFO )
      
      sol(1:neq) = rr(1:neq)
      sol((nx-1)*neq+1:nx*neq)= rr(neq+1:2*neq)

      sol(neq+1:(nx-1)*neq) = invbm(:,1) - matmul(invbm(:,2:2*neq+1),rr)
      
      deallocate(IPIV,Pm,rm,rr,invbm,lm,dm)

c     End program

      end SUBROUTINE BTDS

c     CLEAN_BTDS
c     ###############################################################
      SUBROUTINE CLEAN_BTDS()

      nullify(B_a,B_b,B_c)
c      deallocate(rm,rr,invbm,lm,dm,Pm,IPIV)

      END SUBROUTINE CLEAN_BTDS

c     BLOCK_TRIDIAG
c     #################################################################
      SUBROUTINE BLOCK_TRIDIAG(AA,BB,CC,RR,UU,N,NEQ,NR,CODE)
!     *****************************************************************
!     Solves for a vector U of length N the tridiagonal linear set
!     M U = R, where A, B and C are the three main diagonals of matrix
!     M(N,N), the other terms are 0. R is the right side vector.
!     *****************************************************************
      implicit none
      integer :: N, NEQ, NR
      real(8) :: AA(N*NEQ,NEQ),BB(N*NEQ,NEQ),CC(N*NEQ,NEQ),RR(N*NEQ,NR)
     .     ,UU(N*NEQ,NR)
      REAL(8) :: BET(NEQ,NEQ),GAM(N*NEQ,NEQ), inv2(NEQ,NEQ)
      integer :: J,J1,J2,Jm1,Jm2,Jp1,Jp2
      INTEGER :: CODE

      call Invert_matrix(BB(1:NEQ,1:NEQ),BET,NEQ,CODE)
      UU(1:NEQ,:)=matmul(BET,RR(1:NEQ,:))

      DO J=2,N                  !Decomposition and forward substitution
         J1 = (J-1)*NEQ+1
         J2 = J*NEQ
         Jm1 = (J-2)*NEQ+1
         Jm2 = (J-1)*NEQ
         GAM(J1:J2,:)=matmul(BET,CC(Jm1:Jm2, :))
         BET=BB(J1:J2,:)-matmul(AA(J1:J2,:),GAM(J1:J2,:))

         call Invert_matrix(BET,inv2,NEQ,CODE)
         BET = inv2
         UU(J1:J2,:)=matmul(BET,RR(J1:J2,:)
     .              -matmul(AA(J1:J2,:),UU(Jm1:Jm2,:)))
      END DO
      
      DO J=N-1,1,-1             !Back substitution
         J1 = (J-1)*NEQ+1
         J2 = J*NEQ
         Jp1= J*NEQ+1
         Jp2= (J+1)*NEQ
         UU(J1:J2,:)=UU(J1:J2,:)-matmul(GAM(Jp1:Jp2,:),UU(Jp1:Jp2,:))
      END DO
      
      CODE=0

      END SUBROUTINE BLOCK_TRIDIAG

c     Invert_matrix
c     #################################################################
      SUBROUTINE Invert_matrix(A,invA,N,CODE)
      implicit none

      integer, intent(in) :: N 
      real(8), intent(in) :: A(N,N)
      real(8), intent(out):: invA(N,N)
      integer, intent(out) :: CODE

!     for lapack
      integer :: LDA,INFO,LWORK
      integer, pointer :: IPIV(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK

!     call lapack
      allocate(IPIV(N))
      LWORK = N*N
      allocate(WORK(LWORK))

!     LU factorization
      invA = A
      LDA = N
      CALL DGETRF(N, N, invA, LDA, IPIV, INFO )

      IF(INFO.LT.0)THEN
         PRINT '(" LU decomposition:  illegal value ")'
         STOP
      ENDIF
      IF(INFO.GT.0)THEN
         WRITE(*,35)INFO,INFO
 35             FORMAT( 'LU decomposition: U(',I4,',',I4,') = 0 ')
      ENDIF

!     inverse
      CALL DGETRI(N, invA, N, IPIV, WORK, LWORK, INFO)

      IF (info.NE.0) THEN
         stop 'Matrix inversion failed!'
      ENDIF

      deallocate(IPIV,WORK)

      CODE = 1

      end SUBROUTINE Invert_matrix

c$$$      SUBROUTINE Invert_matrix(M,invM,N,CODE)
c$$$      implicit none
c$$$      integer, intent(in) :: N
c$$$      real(8), intent(in) :: M(N,N)
c$$$      real(8), intent(out):: invM(N,N)
c$$$      integer, intent(out) :: CODE
c$$$      real(8) ::det
c$$$      
c$$$      select case(N)
c$$$      case (2)
c$$$         det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
c$$$         invM(1,1) = M(2,2)
c$$$         invM(2,2) = M(1,1)
c$$$         invM(1,2) =-M(1,2)
c$$$         invM(2,1) =-M(2,1)
c$$$         invM = invM/det
c$$$         CODE = 0
c$$$      case (4)
c$$$         det= M(1,4)*M(2,3)*M(3,2)*M(4,1) - M(1,3)*M(2,4)*M(3,2)*M(4,1) 
c$$$     .        - M(1,4)*M(2,2)*M(3,3)*M(4,1) + 
c$$$     .        M(1,2)*M(2,4)*M(3,3)*M(4,1) + M(1,3)*M(2,2)*M(3,4)*M(4,1) 
c$$$     .        - M(1,2)*M(2,3)*M(3,4)*M(4,1) - 
c$$$     .        M(1,4)*M(2,3)*M(3,1)*M(4,2) + M(1,3)*M(2,4)*M(3,1)*M(4,2) 
c$$$     .        + M(1,4)*M(2,1)*M(3,3)*M(4,2) - 
c$$$     .        M(1,1)*M(2,4)*M(3,3)*M(4,2) - M(1,3)*M(2,1)*M(3,4)*M(4,2) 
c$$$     .        + M(1,1)*M(2,3)*M(3,4)*M(4,2) + 
c$$$     .        M(1,4)*M(2,2)*M(3,1)*M(4,3) - M(1,2)*M(2,4)*M(3,1)*M(4,3) 
c$$$     .        - M(1,4)*M(2,1)*M(3,2)*M(4,3) + 
c$$$     .        M(1,1)*M(2,4)*M(3,2)*M(4,3) + M(1,2)*M(2,1)*M(3,4)*M(4,3) 
c$$$     .        - M(1,1)*M(2,2)*M(3,4)*M(4,3) - 
c$$$     .        M(1,3)*M(2,2)*M(3,1)*M(4,4) + M(1,2)*M(2,3)*M(3,1)*M(4,4) 
c$$$     .        + M(1,3)*M(2,1)*M(3,2)*M(4,4) - 
c$$$     .        M(1,1)*M(2,3)*M(3,2)*M(4,4) - M(1,2)*M(2,1)*M(3,3)*M(4,4) 
c$$$     .        + M(1,1)*M(2,2)*M(3,3)*M(4,4)
c$$$
c$$$         invM(1,1) = -M(2,4)*M(3,3)*M(4,2) + M(2,3)*M(3,4)*M(4,2) 
c$$$     .        + M(2,4)*M(3,2)*M(4,3) - 
c$$$     .        M(2,2)*M(3,4)*M(4,3) - M(2,3)*M(3,2)*M(4,4) 
c$$$     .        + M(2,2)*M(3,3)*M(4,4)
c$$$         invM(1,2) =  M(1,4)*M(3,3)*M(4,2) - M(1,3)*M(3,4)*M(4,2) 
c$$$     .        - M(1,4)*M(3,2)*M(4,3) + 
c$$$     .        M(1,2)*M(3,4)*M(4,3) + M(1,3)*M(3,2)*M(4,4) 
c$$$     .        - M(1,2)*M(3,3)*M(4,4)
c$$$         invM(1,3) = -M(1,4)*M(2,3)*M(4,2) + M(1,3)*M(2,4)*M(4,2) 
c$$$     .        + M(1,4)*M(2,2)*M(4,3) - 
c$$$     .        M(1,2)*M(2,4)*M(4,3) - M(1,3)*M(2,2)*M(4,4)
c$$$     .        + M(1,2)*M(2,3)*M(4,4)
c$$$         invM(1,4) =  M(1,4)*M(2,3)*M(3,2) - M(1,3)*M(2,4)*M(3,2) 
c$$$     .        - M(1,4)*M(2,2)*M(3,3) + 
c$$$     .        M(1,2)*M(2,4)*M(3,3) + M(1,3)*M(2,2)*M(3,4) 
c$$$     .        - M(1,2)*M(2,3)*M(3,4)
c$$$         invM(2,1) =  M(2,4)*M(3,3)*M(4,1) - M(2,3)*M(3,4)*M(4,1) 
c$$$     .        - M(2,4)*M(3,1)*M(4,3) + 
c$$$     .        M(2,1)*M(3,4)*M(4,3) + M(2,3)*M(3,1)*M(4,4) 
c$$$     .        - M(2,1)*M(3,3)*M(4,4)
c$$$         invM(2,2) = -M(1,4)*M(3,3)*M(4,1) + M(1,3)*M(3,4)*M(4,1) 
c$$$     .        + M(1,4)*M(3,1)*M(4,3) - 
c$$$     .        M(1,1)*M(3,4)*M(4,3) - M(1,3)*M(3,1)*M(4,4) 
c$$$     .        + M(1,1)*M(3,3)*M(4,4)
c$$$         invM(2,3) =  M(1,4)*M(2,3)*M(4,1) - M(1,3)*M(2,4)*M(4,1) 
c$$$     .        - M(1,4)*M(2,1)*M(4,3) + 
c$$$     .        M(1,1)*M(2,4)*M(4,3) + M(1,3)*M(2,1)*M(4,4) 
c$$$     .        - M(1,1)*M(2,3)*M(4,4)
c$$$         invM(2,4) = -M(1,4)*M(2,3)*M(3,1) + M(1,3)*M(2,4)*M(3,1) 
c$$$     .        + M(1,4)*M(2,1)*M(3,3) - 
c$$$     .        M(1,1)*M(2,4)*M(3,3) - M(1,3)*M(2,1)*M(3,4) 
c$$$     .        + M(1,1)*M(2,3)*M(3,4)
c$$$         invM(3,1) = -M(2,4)*M(3,2)*M(4,1) + M(2,2)*M(3,4)*M(4,1) 
c$$$     .        + M(2,4)*M(3,1)*M(4,2) - 
c$$$     .        M(2,1)*M(3,4)*M(4,2) - M(2,2)*M(3,1)*M(4,4) 
c$$$     .        + M(2,1)*M(3,2)*M(4,4)
c$$$         invM(3,2) =  M(1,4)*M(3,2)*M(4,1) - M(1,2)*M(3,4)*M(4,1)
c$$$     .        - M(1,4)*M(3,1)*M(4,2) + 
c$$$     .        M(1,1)*M(3,4)*M(4,2) + M(1,2)*M(3,1)*M(4,4) 
c$$$     .        - M(1,1)*M(3,2)*M(4,4)         
c$$$         invM(3,3) = -M(1,4)*M(2,2)*M(4,1) + M(1,2)*M(2,4)*M(4,1) 
c$$$     .        + M(1,4)*M(2,1)*M(4,2) - 
c$$$     .        M(1,1)*M(2,4)*M(4,2) - M(1,2)*M(2,1)*M(4,4) 
c$$$     .        + M(1,1)*M(2,2)*M(4,4)
c$$$         invM(3,4) =  M(1,4)*M(2,2)*M(3,1) - M(1,2)*M(2,4)*M(3,1) 
c$$$     .        - M(1,4)*M(2,1)*M(3,2) + 
c$$$     .        M(1,1)*M(2,4)*M(3,2) + M(1,2)*M(2,1)*M(3,4) 
c$$$     .        - M(1,1)*M(2,2)*M(3,4)
c$$$         invM(4,1) =  M(2,3)*M(3,2)*M(4,1) - M(2,2)*M(3,3)*M(4,1) 
c$$$     .        - M(2,3)*M(3,1)*M(4,2) + 
c$$$     .        M(2,1)*M(3,3)*M(4,2) + M(2,2)*M(3,1)*M(4,3) 
c$$$     .        - M(2,1)*M(3,2)*M(4,3)
c$$$         invM(4,2) = -M(1,3)*M(3,2)*M(4,1) + M(1,2)*M(3,3)*M(4,1) 
c$$$     .        + M(1,3)*M(3,1)*M(4,2) - 
c$$$     .        M(1,1)*M(3,3)*M(4,2) - M(1,2)*M(3,1)*M(4,3) 
c$$$     .        + M(1,1)*M(3,2)*M(4,3)
c$$$         invM(4,3) =  M(1,3)*M(2,2)*M(4,1) - M(1,2)*M(2,3)*M(4,1) 
c$$$     .        - M(1,3)*M(2,1)*M(4,2) + 
c$$$     .    M(1,1)*M(2,3)*M(4,2) + M(1,2)*M(2,1)*M(4,3) 
c$$$     .        - M(1,1)*M(2,2)*M(4,3)
c$$$         invM(4,4) = -M(1,3)*M(2,2)*M(3,1) + M(1,2)*M(2,3)*M(3,1) 
c$$$     .        + M(1,3)*M(2,1)*M(3,2) - 
c$$$     .    M(1,1)*M(2,3)*M(3,2) - M(1,2)*M(2,1)*M(3,3) 
c$$$     .        + M(1,1)*M(2,2)*M(3,3)
c$$$
c$$$         invM = invM/det
c$$$
c$$$      case default
c$$$         print *, "could not invert Matrix with dimension =",N
c$$$         print *, "Quit SUBROUTINE Invert_matrix()"
c$$$         CODE = 1
c$$$         stop
c$$$      end select
c$$$      END SUBROUTINE Invert_matrix

      end Module block_tridiagonal_solver
      

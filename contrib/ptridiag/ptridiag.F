cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Daniel R. Reynolds
C UC San Diego, Mathematics
C
C Direct parallel tridiagonal solvers:
C 
C ptridiag_parallel implements an MPI-based parallel solution 
C    approach for distributed, possibly periodic, tridiagonal 
C    systems.  These systems are defined through the three vectors 
C    a, b, and c on each processor, corresponding to the first 
C    subdiagonal, diagonal, and first superdiagonals of the 
C    tridiagonal matrix.  The two assumptions on the system are 
C       (1) the system is full rank
C       (2) in the 1D communicator, the processors with rank one 
C           lower and one higher correspond to the portions of the
C           matrix just to the left and right of a given processor.
C    This solver treats one tridiagonal matrix at a time, through 
C    a three-phase approach: 
C       (a) local solution on each processor, 
C       (b) communication of boundary information to form a 
C           global interface system on each processor, and 
C       (c) solution of the global interface systems and update 
C           of local processor solutions.
C 
C ptridiag_parallel_phase1 performs only the first local solution 
C    phase (a) of the above algorithm.  This should be called once 
C    for each tridiagonal system, with the resulting temporary 
C    data stored appropriately by the calling program.
C
C ptridiag_parallel_phase2 performs only the communication phase 
C    (b) of the above algorithm.  This may be called once for each 
C    tridiagonal system, or temporary data may be stored up to 
C    allow a single call on a larger group of tridiagonal systems.
C
C ptridiag_parallel_phase3 performs only the global correction 
C    solve and local updates (c) for the above algorithm.  This 
C    should be called once for each tridiagonal system under 
C    consideration, using the appropriate data resulting from 
C    the communication phase above.
C
C ptridiag_mrhs solves periodic tridiagonal systems with multiple 
C    right-hand sides in serial (single processor only).
C
C
C NOTE: the two routines ptridiag_parallel and 
C    ptridiag_parallel_phase3 both must be linked with the LAPACK 
C    routines dgetrf and dgetrs at compile time, in order to allow 
C    for solution of the dense global interface systems from phase 
C    (c) of the algorithm.
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine ptridiag_parallel(comm1d,iprocs,iloc,nloc,a,b,c,x,r)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Solves the parallel periodic tridiagonal system
C        R_in(i) = A_in(i)*x(i-1) + B_in(i)*x(i) + C_in(i)*x(i+1)
C     The period is n, the same as the dimension of the matrix.  
C     Assumes that in the 1-dimensional communicator, the processor 
C     one above is the right neighbor and the processor one below 
C     is the left neighbor of this processor in the direction.
C
C     NOTE: this routine expects to be linked with the LAPACK 
C     routines dgetrf and dgetrs, used in solving the small dense 
C     global coupling matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mpi
      implicit none
      integer:: comm1d,iprocs,iloc,nloc,i,j,ipiv(2*iprocs),Glen,ierr
      double precision:: a(nloc), b(nloc), c(nloc), r(nloc), x(nloc)
      double precision:: mrhs(nloc,3), msol(nloc,3), Gb(2*iprocs)
      double precision:: GM(2*iprocs,2*iprocs)
      double precision:: alo, chi, GxL, GxLL, GxR, GxRR
      double precision:: sbuf(6), rbuf(6*iprocs)

c     shift local rank by one
      iloc = iloc+1

c     set global coupling matrix size
      Glen = 2*iprocs

c     set relevant internal arrays
      msol = 0.d0
      mrhs = 0.d0
      Gb = 0.d0
      GM = 0.d0

c     store periodicity components and reset local tridiagonal 
c     matrix arrays without periodicity
      alo = a(1)
      chi = c(nloc)
      a(1) = 0.d0
      c(nloc) = 0.d0

c     set up multiple rhs array with periodicity stuff and original rhs
c       column: 1=rhs, 2=(0,...,0,chi), 3=(alo,0,...,0)
      mrhs(1,1) = r(1)
      mrhs(1,3) = alo
      do i=2,nloc-1
         mrhs(i,1) = r(i)
      enddo
      mrhs(nloc,1) = r(nloc)
      mrhs(nloc,2) = chi

c     solve the local tridiagonal systems to get 
c       msol(:,1) = Tloc^{-1}rhs
c       msol(:,2) = Tloc^{-1}(0,...,0,chi)
c       msol(:,3) = Tloc^{-1}(alo,0,...,0)
      call ptridiag_mrhs(nloc, a, b, c, 3, msol, mrhs)

c     set up the global 1D system on each processor
      if (iprocs==1) then

c        if 1 proc only, then 'global' solution handles periodicity
         if ((alo /= 0.d0) .or. (chi /= 0.d0)) then

c           global matrix and solution vector (directly invert)
            GM(1,1) = 1.d0+msol(1,2)
            GM(2,2) = 1.d0+msol(nloc,3)
            GM(1,2) = msol(1,3)
            GM(2,1) = msol(nloc,2)
            Gb(1) = (GM(2,2)*msol(1,1) - GM(1,2)*msol(nloc,1)) 
     &            / (GM(1,1)*GM(2,2) - GM(1,2)*GM(2,1))
            Gb(2) = (GM(1,1)*msol(nloc,1) - GM(2,1)*msol(1,1)) 
     &            / (GM(1,1)*GM(2,2) - GM(1,2)*GM(2,1))
            GxLL = Gb(2)
            GxL  = Gb(1)
            GxR  = Gb(2)
            GxRR = Gb(1)

c        if 1 proc and non-periodic, copy sol into global solution
         else
            GxLL = 0.d0
            GxL  = msol(1,1)
            GxR  = msol(nloc,1)
            GxRR = 0.d0
         endif

c     multiple processor case remains -- communicate to set up 
c     global matrix and vector, each processor then solves the 
c     same small dense system to obtain their local solution update
      else

c        send buffer contains matrix and rhs entries for boundary
         sbuf(1) = msol(1,3)     ! lower diag, left end
         sbuf(2) = msol(nloc,3)  ! lower diag, right end
         sbuf(3) = msol(1,2)     ! upper diag, left end
         sbuf(4) = msol(nloc,2)  ! upper diag, right end
         sbuf(5) = msol(1,1)     ! rhs, left end
         sbuf(6) = msol(nloc,1)  ! rhs, right end

c        adjust for non-periodicity
         if ((alo == 0.0) .and. (iloc==1)) then
            sbuf(1) = 0.d0
            sbuf(2) = 0.d0
         endif
         if ((chi == 0.0) .and. (iloc==iprocs)) then
            sbuf(3) = 0.d0
            sbuf(4) = 0.d0
         endif

c        communicate matrix information to all procs
         call MPI_AllGather(sbuf, 6, MPI_DOUBLE_PRECISION, 
     &        rbuf, 6, MPI_DOUBLE_PRECISION, comm1d, ierr)
         if (ierr /= 0) then
            write(0,*) 'ptridiag_parallel error in MPI_AllGather =',ierr
            return
         endif

c        set up block matrix and rhs
         GM(1,1)        = 1.d0
         GM(1,3)        = rbuf(3)
         GM(1,2*iprocs) = rbuf(1)
         GM(2,2)        = 1.d0
         GM(2,3)        = rbuf(4)
         GM(2,2*iprocs) = rbuf(2)
         Gb(1)          = rbuf(5)
         Gb(2)          = rbuf(6)
         do i=2,iprocs-1
            GM(2*i-1,2*i-2) = rbuf((i-1)*6+1)
            GM(2*i-1,2*i-1) = 1.d0
            GM(2*i-1,2*i+1) = rbuf((i-1)*6+3)
            GM(2*i,2*i-2)   = rbuf((i-1)*6+2)
            GM(2*i,2*i)     = 1.d0
            GM(2*i,2*i+1)   = rbuf((i-1)*6+4)
            Gb(2*i-1)       = rbuf((i-1)*6+5)
            Gb(2*i)         = rbuf((i-1)*6+6)
         enddo
         GM(2*iprocs-1,1)          = rbuf((iprocs-1)*6+3)
         GM(2*iprocs-1,2*iprocs-2) = rbuf((iprocs-1)*6+1)
         GM(2*iprocs-1,2*iprocs-1) = 1.d0
         GM(2*iprocs,1)            = rbuf((iprocs-1)*6+4)
         GM(2*iprocs,2*iprocs-2)   = rbuf((iprocs-1)*6+2)
         GM(2*iprocs,2*iprocs)     = 1.d0
         Gb(2*iprocs-1)            = rbuf((iprocs-1)*6+5)
         Gb(2*iprocs)              = rbuf((iprocs-1)*6+6)

c        solve matrix system and extract relevant global solution values
         call dgetrf(Glen,Glen,GM,Glen,ipiv,ierr)
         if (ierr /= 0) then
            write(0,*) 'ptridiag_parallel error in dgetrf =',ierr
            return
         endif

c           solve factored system
         call dgetrs('N',Glen,1,GM,Glen,ipiv,Gb,Glen,ierr)
         if (ierr /= 0) then
            write(0,*) 'ptridiag_parallel error in dgetrs =',ierr
            return
         endif

c           handle contributions if at left end of domain
         if (iloc==1) then
            if (alo == 0.d0) then
               GxLL = 0.d0
            else
               GxLL = Gb(iprocs*2)
            endif
         else
            GxLL = Gb((iloc-1)*2)
         endif

c           handle contributions if at right end of domain
         if (iloc==iprocs) then
            if (chi==0.d0) then
               GxRR = 0.d0
            else
               GxRR = Gb(1)
            endif
         else
            GxRR = Gb(iloc*2+1)
         endif

c           handle local contributions
         GxL = Gb((iloc-1)*2+1)
         GxR = Gb((iloc-1)*2+2)

      endif

c     update local solution using global solution components
      x(1) = GxL
      do i=2,nloc-1,1
         x(i) = msol(i,1) - msol(i,3)*GxLL - msol(i,2)*GxRR
      enddo
      x(nloc) = GxR

      return

      end subroutine ptridiag_parallel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






      subroutine ptridiag_parallel_phase1(iprocs, iloc, nloc, a, b, 
     &                                    c, r, x, xL, xR, sbuf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Solves the parallel periodic tridiagonal system
C        R_in(i) = A_in(i)*x(i-1) + B_in(i)*x(i) + C_in(i)*x(i+1)
C     The period is n, the same as the dimension of the matrix.  
C     Assumes that in the 1-dimensional communicator, the processor 
C     one above is the right neighbor and the processor one below 
C     is the left neighbor of this processor in the direction.
C
C     NOTE: these should only be called in lieu of ptridiag_parallel 
C       if in fact there are multiple processers involved in the 
C       parallel tridiagonal solve.
C
C     Phase 1: the local solve and packing the local buffer with 
C          information for the global system.
C
C     Inputs:
C       iprocs - number of processors in this direction
C       iloc   - location in 1D topology i.e. iloc in (1:iprocs)
C       nloc   - length of local portion of overall problem
C       a      - subdiagonal of tridiagonal system
C       b      - diagonal of tridiagonal system
C       c      - superdiagonal of tridiagonal system
C       r      - rhs of tridiagonal system
C
C     Outputs:
C       x      - local solution vector (first pass)
C       xL     - vector corresp. to how local components depend on 
C                left neighbor
C       xR     - vector corresp. to how local components depend on 
C                right neighbor
C       sbuf   - buffer containing local information for global 
C                correction systems over entire face
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mpi
      implicit none
      integer :: iprocs, iloc, nloc, i, ierr
      double precision :: a(nloc), b(nloc), c(nloc), r(nloc)
      double precision :: x(nloc), xL(nloc), xR(nloc)
      double precision :: mrhs(nloc,3), msol(nloc,3)
      double precision :: alo, chi, sbuf(6)

c     set relevant internal arrays
      msol = 0.d0
      mrhs = 0.d0

c     store periodicity components and reset local tridiagonal 
c     matrix arrays without periodicity
      alo = a(1)
      chi = c(nloc)
      a(1) = 0.d0
      c(nloc) = 0.d0

c     set up multiple rhs array with periodicity stuff and original rhs
c       column: 1=rhs, 2=(0,...,0,chi), 3=(alo,0,...,0)
      mrhs(1,1) = r(1)
      mrhs(1,3) = alo
      do i=2,nloc-1
         mrhs(i,1) = r(i)
      enddo
      mrhs(nloc,1) = r(nloc)
      mrhs(nloc,2) = chi

c     solve the local tridiagonal systems to get 
c       msol(:,1) = Tloc^{-1}rhs
c       msol(:,2) = Tloc^{-1}(0,...,0,chi)
c       msol(:,3) = Tloc^{-1}(alo,0,...,0)
      call ptridiag_mrhs(nloc, a, b, c, 3, msol, mrhs)

c     insert msol components into various solution pieces:
c        msol(:,1) -> x (current solution values)
c        msol(:,2) -> xR (how right nbor affects interior)
c        msol(:,3) -> xL (how left  nbor affects interior)
      do i=1,nloc
         x  = msol(:,1)
      enddo
      do i=1,nloc
         xR = msol(:,2)
      enddo
      do i=1,nloc
         xL = msol(:,3)
      enddo

c     pack send buffer with matrix and rhs entries for boundary
      sbuf(1) = msol(1,3)       ! lower diag, left end
      sbuf(2) = msol(nloc,3)    ! lower diag, right end
      sbuf(3) = msol(1,2)       ! upper diag, left end
      sbuf(4) = msol(nloc,2)    ! upper diag, right end
      sbuf(5) = msol(1,1)       ! rhs, left end
      sbuf(6) = msol(nloc,1)    ! rhs, right end

c     adjust for non-periodicity
      if ((alo == 0.0) .and. (iloc==1)) then
         sbuf(1) = 0.d0
         sbuf(2) = 0.d0
      endif
      if ((chi == 0.0) .and. (iloc==iprocs)) then
         sbuf(3) = 0.d0
         sbuf(4) = 0.d0
      endif

      return

      end subroutine ptridiag_parallel_phase1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






      subroutine ptridiag_parallel_phase2(comm1d,iprocs,nynz,sbuf,rbuf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Solves the parallel periodic tridiagonal system
C        R_in(i) = A_in(i)*x(i-1) + B_in(i)*x(i) + C_in(i)*x(i+1)
C     The period is n, the same as the dimension of the matrix.  
C     Assumes that in the 1-dimensional communicator, the processor 
C     one above is the right neighbor and the processor one below 
C     is the left neighbor of this processor in the direction.
C
C     NOTE: these should only be called in lieu of ptridiag_parallel if 
C       in fact there are multiple processers involved in the 
C       parallel tridiagonal solve.
C
C     Phase 2: communication of all the local send buffers over the 
C          2D neighbor interfaces to obtain the relevant 
C          information for the global correction systems.
C
C     Inputs:
C       comm1d - 1-dimensional communicator
C       iprocs - number of processors in this direction
C       nynz   - size of 2D face orthogonal to this direction
C       sbuf   - buffer containing local information for global 
C                correction systems over entire face
C
C     Outputs:
C       rbuf   - buffer containing information for global correction 
C                systems over entire face
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mpi
      implicit none
      integer :: comm1d, iprocs, nynz, ierr
      double precision :: sbuf(6*nynz), rbuf(6*nynz*iprocs)

c     communicate matrix information to all procs
      call MPI_AllGather(sbuf, 6*nynz, MPI_DOUBLE_PRECISION, 
     &     rbuf, 6*nynz, MPI_DOUBLE_PRECISION, comm1d, ierr)
      if (ierr /= 0) then
         write(0,*) 'ptridiag_parallel_phase2 error in MPI_AllGather =',
     &        ierr
         return
      endif

      return

      end subroutine ptridiag_parallel_phase2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






      subroutine ptridiag_parallel_phase3(iprocs, iloc, nloc, x, 
     &                                    xL, xR, perdc, rbuf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Solves the parallel periodic tridiagonal system
C         R_in(i) = A_in(i)*x(i-1) + B_in(i)*x(i) + C_in(i)*x(i+1)
C     The period is n, the same as the dimension of the matrix.  
C     Assumes that in the 1-dimensional communicator, the processor 
C     one above is the right neighbor and the processor one below is 
C     the left neighbor of this processor in the direction.
C     
C     NOTE: these should only be called in lieu of ptridiag_parallel if 
C       in fact there are multiple processers involved in the 
C       parallel tridiagonal solve.
C
C     Phase 3: construction and solution of the global interface 
C          system to correct the local tridiagonal solution.
C
C     Inputs:
C       iprocs - number of processors in this direction
C       iloc   - location in 1D topology i.e. iloc in (1:iprocs)
C       nloc   - length of local portion of overall problem
C       x      - local solution vector (after 1st phase)
C       xL     - vector corresp. to how local components depend on 
C                left neighbor
C       xR     - vector corresp. to how local components depend on 
C                right neighbor
C       perdc  - flag denoting whether dimension is periodic or not
C                  1 => periodic,  0 => otherwise
C       rbuf   - buffer containing information for global correction 
C                system
C
C     Outputs:
C       x      - local solution vector (actual solution)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mpi
      implicit none
      integer :: iprocs,iloc,nloc,i,ipiv(2*iprocs),Glen,ierr,perdc
      double precision :: x(nloc), xL(nloc), xR(nloc)
      double precision :: Gb(2*iprocs), GM(2*iprocs,2*iprocs)
      double precision :: GxL, GxLL, GxR, GxRR
      double precision :: rbuf(6*iprocs)

c     set global coupling matrix size
      Glen = 2*iprocs

c     set up block matrix and rhs
      GM             = 0.d0
      Gb             = 0.d0
      GM(1,1)        = 1.d0
      GM(1,3)        = rbuf(3)
      GM(1,2*iprocs) = rbuf(1)
      GM(2,2)        = 1.d0
      GM(2,3)        = rbuf(4)
      GM(2,2*iprocs) = rbuf(2)
      Gb(1)          = rbuf(5)
      Gb(2)          = rbuf(6)
      do i=2,iprocs-1
         GM(2*i-1,2*i-2) = rbuf((i-1)*6+1)
         GM(2*i-1,2*i-1) = 1.d0
         GM(2*i-1,2*i+1) = rbuf((i-1)*6+3)
         GM(2*i,2*i-2)   = rbuf((i-1)*6+2)
         GM(2*i,2*i)     = 1.d0
         GM(2*i,2*i+1)   = rbuf((i-1)*6+4)
         Gb(2*i-1)       = rbuf((i-1)*6+5)
         Gb(2*i)         = rbuf((i-1)*6+6)
      enddo
      GM(2*iprocs-1,1)          = rbuf((iprocs-1)*6+3)
      GM(2*iprocs-1,2*iprocs-2) = rbuf((iprocs-1)*6+1)
      GM(2*iprocs-1,2*iprocs-1) = 1.d0
      GM(2*iprocs,1)            = rbuf((iprocs-1)*6+4)
      GM(2*iprocs,2*iprocs-2)   = rbuf((iprocs-1)*6+2)
      GM(2*iprocs,2*iprocs)     = 1.d0
      Gb(2*iprocs-1)            = rbuf((iprocs-1)*6+5)
      Gb(2*iprocs)              = rbuf((iprocs-1)*6+6)

c     solve matrix system and extract relevant global solution values
      call dgetrf(Glen,Glen,GM,Glen,ipiv,ierr)
      if (ierr /= 0) then
         write(0,*) 'ptridiag_parallel_phase3 error in dgetrf =',ierr
         return
      endif

c     solve factored system
      call dgetrs('N',Glen,1,GM,Glen,ipiv,Gb,Glen,ierr)
      if (ierr /= 0) then
         write(0,*) 'ptridiag_parallel_phase3 error in dgetrs =',ierr
         return
      endif

c     handle contributions if at left end of domain
      if (iloc==1) then
         if (perdc == 1) then
            GxLL = Gb(iprocs*2)
         else
            GxLL = 0.d0
         endif
      else
         GxLL = Gb((iloc-1)*2)
      endif

c     handle contributions if at right end of domain
      if (iloc==iprocs) then
         if (perdc == 1) then
            GxRR = Gb(1)
         else
            GxRR = 0.d0
         endif
      else
         GxRR = Gb(iloc*2+1)
      endif

c     handle local contributions
      GxL = Gb((iloc-1)*2+1)
      GxR = Gb((iloc-1)*2+2)

c     update local solution using global solution components
      x(1) = GxL
      do i=2,nloc-1,1
         x(i) = x(i) - xL(i)*GxLL - xR(i)*GxRR
      enddo
      x(nloc) = GxR

      return

      end subroutine ptridiag_parallel_phase3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





      subroutine ptridiag_mrhs(n, a, b, c, m, x, r)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Solves the periodic tridiagonal system
C         R_in(i,j) = A_in(i)*x(i-1,j) + B_in(i)*x(i,j) 
C                     + C_in(i)*x(i+1,j)
C     where R_in(:,j), j=1:m, corresponds to one of possibly 
C     multiple right-hand sides.
C     The period is n, the same as the dimension of the matrix.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer:: n, m, i, j
      double precision:: a(n), b(n), c(n), r(n,m), x(n,m)
      
C     Eliminate a on the row below
      a(1) = a(1)/b(1)
      c(1) = c(1)/b(1)
      do j = 1, m
         r(1,j) = r(1,j)/b(1)
      enddo

      do i = 2, n-1
         b(i) = b(i)-a(i)*c(i-1)
         do j = 1, m
            r(i,j) = r(i,j)-a(i)*r(i-1,j)
         enddo
         
         a(i) = -a(i)*a(i-1)/b(i)
         c(i) = c(i)/b(i)
         do j = 1, m
            r(i,j) = r(i,j)/b(i)
         end do
      enddo
      
      c(n-1) = c(n-1)-a(n-1)
      a(n-1) = c(n-1)     

      b(n) = b(n)-a(n)*c(n-1)
      do j = 1, m
         r(n,j) = (r(n,j)-a(n)*r(n-1,j))/b(n)
      enddo

      a(n) = c(n)/b(n)
       
C     Eliminate nth column:
      do i = n-1, 2, -1
         do j = 1, m
            r(i,j) = r(i,j)-r(n,j)*a(i)
         enddo
         a(i) = -a(n)*a(i) 
      end do
      
C     Eliminate c on the rows above:
      do i = n-1, 3, -1
         do j = 1, m
            r(i-1,j) = r(i-1,j)-r(i,j)*c(i-1)
         enddo
         a(i-1) = a(i-1)-a(i)*c(i-1)
      end do
      
C     Back substitution:
      do j = 1, m
         x(1,j) = (r(1,j) - r(n,j)*a(1) - r(2,j)*c(1)) 
     &        / (1.d0 - a(n)*a(1) - a(2)*c(1))
      enddo
      do i = 2, n
         do j = 1, m
            x(i,j) = r(i,j)-a(i)*x(1,j)
         end do
      enddo


      return

      end subroutine ptridiag_mrhs
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

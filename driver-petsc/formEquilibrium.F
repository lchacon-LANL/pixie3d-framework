c formEquilibrium
c######################################################################
      subroutine formEquilibrium(array,imin,imax,jmin,jmax,kmin,kmax)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use mpi

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c FIX PARALLEL
      integer(4)      :: my_rank,ierr
      character*(5)   :: ext
c FIX PARALLEL

c Local variables

      type(petsc_array) :: petscarray

c Begin program

      urecord    = 25

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
      recordfile = 'record.bin'//trim(ext)//'.bin'
#else
      recordfile = 'record.bin'
#endif

c Check for autoinitializing parameters

      pi = acos(-1d0)

      if (maxitnwt.eq.0) 
     .      maxitnwt = max(floor(1.5*log(rtol)/log(tolgm)),10)

      alpha = 1. - cnfactor

      dtbase = dt   

      sm_flag= 0
      if (cnfactor.lt.0d0) then
        sm_flag= 1
      else
        sm_pass= 0
      endif

      cnf_d = 1d0

c Initialize vector dimensions for global and local problems

      ihig = imax
      ilog = imin
      jhig = jmax
      jlog = jmin
      khig = kmax
      klog = kmin

      call setVectorDimensions

c Allocate constant arrays

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Initialize MG and create grid

      call createGrid(ilog,ihig,jlog,jhig,klog,khig,nxd,nyd,nzd)

cc      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
cc      call MPI_Finalize(ierr)
cc      if (my_rank == 3) then
cc        write (*,*) 'Grid in processor',my_rank
cc        call checkGrid
cc      endif
cc      stop

c Create nonlinear solver

      call createNonlinearSolver

c Create nonlinear function

      call createNonlinearFunction

c Create equilibrium u_0

      call createEquilibrium

c Initialize old time solution

      u_n = u_0   !Overloaded assignment

c Transfer to Petsc format

      petscarray = u_n

      array = petscarray%array(ilo:ihi,jlo:jhi,klo:khi)

      call deallocatePetscType(petscarray)

c End program

      contains

      function int2char(n) result (chr)

      integer(4)   :: n
      character(10):: chr

      integer(4)   :: i,exp,k
      character(3) :: c

      if (abs(n) > 0) then
         exp = int(log(float(n))/log(1d1))
      else
         exp = 0
      endif

      if (n >= 0) then
        chr=''
      else
        chr='-'
      endif

      k = abs(n)
      do i=exp,0,-1
         j = k/10**i
         c = achar(48+j)
         chr = trim(chr)//trim(c)
         if (i > 0 ) k = mod(k,j*10**i)
      enddo

      end function int2char

      end subroutine formEquilibrium


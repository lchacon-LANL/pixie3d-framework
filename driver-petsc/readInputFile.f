c readInputFile
c######################################################################
      subroutine readInputFile(inputdata)

c----------------------------------------------------------------------
c     Reads input file
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use precond_setup

      implicit none

c Call variables

      type :: indata
        real(8)    :: tolgm
        real(8)    :: rtol
        real(8)    :: atol
        real(8)    :: damp
        real(8)    :: dt
        real(8)    :: tmax
        integer(4) :: ilevel
        integer(4) :: nxd
        integer(4) :: nyd
        integer(4) :: nzd
        integer(4) :: npx
        integer(4) :: npy
        integer(4) :: npz
        integer(4) :: numtime
        integer(4) :: maxitnwt
        integer(4) :: maxksp
        integer(4) :: maxitgm
        integer(4) :: method
        integer(4) :: global
        integer(4) :: iguess
        integer(4) :: bcsi(6)
        logical    :: user_PC
      end type indata

      type(indata) :: inputdata

c Local variables

c Begin program

      g_pack%dim(:)%pack = .false.

      call readInput

c Determine processor allocation

      call processorAlloc

c Define structure components

      inputdata%ilevel   = ilevel
      inputdata%nxd      = nxd
      inputdata%nyd      = nyd
      inputdata%nzd      = nzd
      inputdata%npx      = npx
      inputdata%npy      = npy
      inputdata%npz      = npz
      inputdata%numtime  = numtime
      inputdata%maxitnwt = maxitnwt
      inputdata%maxksp   = maxksp 
      inputdata%maxitgm  = maxitgm
      inputdata%method   = method 
      inputdata%global   = global 
      inputdata%iguess   = iguess

      where (bcond == PER)
        inputdata%bcsi = 1
      elsewhere
        inputdata%bcsi = 0
      end where

      inputdata%tolgm    = tolgm  
      inputdata%rtol     = rtol   
      inputdata%atol     = atol   
      inputdata%damp     = damp   
      inputdata%dt       = dt   
      inputdata%tmax     = tmax

      if (precon == 'id') then
        inputdata%user_PC = .false.
      else
        inputdata%user_PC = .true.
      endif

c Write structure

cc      write (*,*) inputdata

c End program

      contains

c     processorAlloc
c     ######################################################################
      subroutine processorAlloc

      use error

      use generalPurposeFunctions

      implicit none

c     Call variables

c     Local variables

      integer(4) :: sum_exp,nd,navg,exp(3),npt(3)

c     Begin program

      !Defaults
      npx = 0
      npy = 0
      npz = 0

      exp = 0

c     Find number of processors

      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)

c     Check MG is an option

      sum_exp = floor(log(1d0*np)/log(2d0))

      if (sum_exp < 1) return  !Only one processor

      if (2**sum_exp /= np) then
        messg = 'Number of processors ('//trim(int2char(np))//
     .          ') unsuitable for MG'
        call pstop('processorAlloc',messg)
      endif

      !Find dimensionality
      nd = 3
      if (nxd == 1) nd = nd-1
      if (nyd == 1) nd = nd-1
      if (nzd == 1) nd = nd-1

      if (nd == 0) then
        messg = 'No available dimensions!'
        call pstop('processorAlloc',messg)
      endif

      !Find exponents
      navg = floor((nxd*nyd*nzd/np)**(1./nd))

      npt = (/nxd,nyd,nzd/)
      exp = floor(log(max(1d0*npt/navg,1d0))/log(2d0))

cc      if (my_rank == 0) write (*,*) npt
cc      if (my_rank == 0) write (*,*) exp

      !Consistency check
      if (sum(exp) /= sum_exp) then
        if (nxd >= nyd .and. nxd >= nzd) then
          exp(1) = sum_exp - exp(2) - exp(3)
        elseif (nyd >= nxd .and. nyd >= nzd) then
          exp(2) = sum_exp - exp(1) - exp(3)
        else
          exp(3) = sum_exp - exp(1) - exp(2)
        endif
      endif

      !Find processor distribution
      npx = 2**exp(1)
      npy = 2**exp(2)
      npz = 2**exp(3)

cc      if (my_rank == 0) write (*,*) sum_exp,nd,navg
cc      if (my_rank == 0) write (*,*) npx,npy,npz
cc      call PetscEnd(mpierr)
cc      stop

c     End program

      end subroutine processorAlloc

      end subroutine readInputFile

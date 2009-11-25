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

      use constants

      use iosetup

      use icond

      implicit none

c Call variables

      type :: indata
        real(8)    :: tolgm
        real(8)    :: rtol
        real(8)    :: atol
        real(8)    :: damp
        real(8)    :: dt
        real(8)    :: tmax
        real(8)    :: mf_eps
        integer    :: nvar
        integer    :: nauxs
        integer    :: nauxv
        integer    :: ilevel
        integer    :: nxd
        integer    :: nyd
        integer    :: nzd
        integer    :: npx
        integer    :: npy
        integer    :: npz
        integer    :: numtime
        integer    :: maxitnwt
        integer    :: maxksp
        integer    :: maxitgm
        integer    :: method
        integer    :: global
        integer    :: iguess
        integer    :: precpass
        integer    :: sm_flag
        integer    :: bcsi(6)
        logical    :: asm_PC
      end type indata

      type(indata) :: inputdata

c Local variables

c Begin program

      g_pack%dim(:)%pack = .false.

c Set PETSc defaults

      npx = inputdata%npx
      npy = inputdata%npy
      npz = inputdata%npz

c Read fortran input file (defines also neqd, NAUXS, NAUXV)

      call readInput

c Initialize MPI

cc      call initMPI(nxd,nyd,nzd)

c Define structure components

      inputdata%nvar     = neqd
      inputdata%nauxs    = NAUXS
      inputdata%nauxv    = NAUXV

      inputdata%ilevel   = ilevel
      inputdata%nxd      = nxd
      inputdata%nyd      = nyd
      inputdata%nzd      = nzd
      inputdata%npx      = npx
      inputdata%npy      = npy
      inputdata%npz      = npz
      inputdata%numtime  = numtime
c      inputdata%maxitnwt = maxitnwt
c      inputdata%maxksp   = maxksp 
c      inputdata%maxitgm  = maxitgm
c      inputdata%method   = method 
c      inputdata%global   = global 
c      inputdata%iguess   = iguess
      inputdata%sm_flag  = sm_flag
c      inputdata%precpass = precpass
c      inputdata%asm_PC   = asm_PC

      where (bcond == PER)
        inputdata%bcsi = 1
      elsewhere
        inputdata%bcsi = 0
      end where

c      inputdata%tolgm    = tolgm  
c      inputdata%rtol     = rtol   
c      inputdata%atol     = atol   
      inputdata%damp     = damp   
      inputdata%dt       = dt   
      inputdata%tmax     = tmax
c      inputdata%mf_eps   = mf_eps

c End program

      end subroutine readInputFile

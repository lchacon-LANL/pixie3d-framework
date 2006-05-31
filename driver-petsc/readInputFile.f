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
        integer(4) :: precpass
        integer(4) :: sm_flag
        integer(4) :: bcsi(6)
      end type indata

      type(indata) :: inputdata

c Local variables

c Begin program

      g_pack%dim(:)%pack = .false.

c Set PETSc defaults

      npx = inputdata%npx
      npy = inputdata%npy
      npz = inputdata%npz

c Read fortran input file

      call readInput

c Initialize MPI

      call initMPI(nxd,nyd,nzd)

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
      inputdata%sm_flag  = sm_flag
      inputdata%precpass = precpass

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

c End program

      end subroutine readInputFile

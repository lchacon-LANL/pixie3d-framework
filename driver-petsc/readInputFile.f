c readInputFile
c######################################################################
      subroutine readInputFile(inputdata)

c----------------------------------------------------------------------
c     Reads input file
c----------------------------------------------------------------------

      use parameters

      use grid

      use grid_mpi

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
        real(8)    :: mf_eps
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
        logical    :: tst_flg
      end type indata

      type(indata) :: inputdata

c Local variables

c Begin program

      g_pack%dim(:)%pack = .false.

c Set PETSc defaults

      npx = inputdata%npx
      npy = inputdata%npy
      npz = inputdata%npz

      tst_flg = inputdata%tst_flg

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
      inputdata%asm_PC   = asm_PC

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
      inputdata%mf_eps   = mf_eps

c End program

      end subroutine readInputFile

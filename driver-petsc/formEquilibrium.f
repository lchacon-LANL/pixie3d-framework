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

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c FIX PARALLEL
cc      character*(*)   :: ext
c FIX PARALLEL

c Local variables

      type(petsc_array) :: petscarray

c Begin program

      urecord    = 25
c FIX PARALLEL
      !In the future, "ext" should be passed via subroutine call from C
cc      recordfile = 'record'//trim(ext)//'.bin'
      recordfile = 'record.bin'
c FIX PARALLEL

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
cc      call checkGrid

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

      end subroutine formEquilibrium


c formEquilibrium
c######################################################################
      subroutine formEquilibrium(array,imin,imax,jmin,jmax
     .                                ,kmin,kmax)

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

      integer(4)      ::imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) ::array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      type(petsc_array) :: petscarray

c Begin program

      urecord    = 25
      recordfile = 'record.bin'

c Read user initializations

cc      call readInput

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

c Initialize vector dimensions

      call setVectorDimensions

c Allocate constant arrays

      allocate(zeros (0:nxdp,0:nydp,0:nzdp))
      allocate(vzeros(0:nxdp,0:nydp,0:nzdp,3))
      allocate(ones  (0:nxdp,0:nydp,0:nzdp))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd)
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

      array = petscarray%array(imin:imax,jmin:jmax,kmin:kmax)

      call deallocatePetscType(petscarray)


c End program

      end subroutine formEquilibrium


c createGridStructures
c######################################################################
      subroutine createGridStructures(patch_var)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      use variables

      use timeStepping

      use newtongm

c      use constants

      use iosetup

c      use generalPurposeFunctions

      implicit none

c Call variables

      type(patch), TARGET :: patch_var

c Local variables

      integer ieq

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      u_0 => gv%u_0

      call setVectorDimensions

c Check for autoinitializing parameters

      pi = acos(-1d0)

cc      if (maxitnwt.eq.0) 
cc     .      maxitnwt = max(floor(1.5*log(rtol)/log(tolgm)),10)
cc
cc      alpha = 1. - cnfactor
cc
cc      dtbase = dt   

c Allocate global variable

      call allocateGlobalVar(gv)

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd,g_pack,gv%gparams)

c Initialize vector dimensions for global and local problems

      call setVectorDimensions

cc      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
cc      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
cc      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))
cc
cc      zeros  = 0d0
cc      vzeros = 0d0
cc      ones   = 1d0

c Create nonlinear solver

      call allocateStructures

c Create nonlinear function

cc      call createNonlinearFunction

      end subroutine createGridStructures

c     formEquilibrium
c######################################################################
      subroutine formEquilibrium(patch_var)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      use variables

      use timeStepping

      use newtongm

c      use constants

      use iosetup

c      use generalPurposeFunctions

      implicit none

      type(patch), TARGET :: patch_var

      INTERFACE
         subroutine setEquilibrium(igx,igy,igz,varray)
         use variable_setup
         integer :: igx,igy,igz
         type(var_array),pointer :: varray
         end subroutine setEquilibrium
      END INTERFACE

c Local variables

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      u_0 => gv%u_0

      call setVectorDimensions

c    Create equilibrium u_0

      call setEquilibrium(1,1,1,u_0)  !External

c End program
 
      end subroutine formEquilibrium

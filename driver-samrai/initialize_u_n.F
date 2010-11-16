c initialize_u_n
c######################################################################
      subroutine initialize_u_n(patch_var)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

c      use constants

      use iosetup

c      use generalpurposefunctions

      implicit none

c Call variables

      real(8)         :: time_to_c
      type(patch), TARGET :: patch_var
      
c Local variables

      integer         :: iminl,imaxl,jminl,jmaxl,kminl,kmaxl,ieq
     .                  ,ierr

c External

c      integer    :: system
c      external      system

c Interfaces

      INTERFACE
         subroutine setInitialCondition(varrayn,varraynp)
         use variable_setup
         type(var_array),pointer :: varrayn,varraynp
         end subroutine setInitialCondition
      END INTERFACE

c Initailize gv

      gv => patch_var

c Initialize old time solution

      call initializeDerivedType(gv%u_n)

c      call equateDerivedType(gv%u_n,gv%u_0)

c End program

      end subroutine initialize_u_n



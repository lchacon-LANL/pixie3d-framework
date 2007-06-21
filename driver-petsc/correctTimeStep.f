c correctTimeStep
c ######################################################################
      subroutine petscCorrectTimeStep(dn,dnh,dnp,ierr,dt_to_c)

      use timeStepping

      use counters

      implicit none

c Call variables

      integer     :: ierr
      real(8)     :: dn(neqd),dnh(neqd),dnp(neqd),dt_to_c

c Local variables

c Begin program

      call correctTimeStep(dn,dnh,dnp,ierr,dt_to_c)

c End program

      end subroutine petscCorrectTimeStep

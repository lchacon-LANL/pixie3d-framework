c defineTSParametersRoutine
c####################################################################
      subroutine defineTSParametersRoutine(cnf,one_over_dt)

c--------------------------------------------------------------------
c     Calculates nonlinear residuals, of the form:
c             dt Ui + Fi(Uj) = 0
c--------------------------------------------------------------------

      use grid

      use variables

      use timeStepping

      implicit none

c Call variables

      real(8)    :: cnf(neqd),one_over_dt(neqd)

c Local variables


c Begin program


c Calculate residuals

      call defineTSParameters(cnf,one_over_dt)

c End program

      end subroutine defineTSParametersRoutine

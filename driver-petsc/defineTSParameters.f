c defineTSParametersRoutine
c####################################################################
      subroutine defineTSParametersRoutine(cnf1,one_over_dt1)

c--------------------------------------------------------------------
c     Calculates nonlinear residuals, of the form:
c             dt Ui + Fi(Uj) = 0
c--------------------------------------------------------------------

      use grid

      use variables

      use timeStepping

      implicit none

c Call variables

      real(8)    :: cnf1(neqd),one_over_dt1(neqd)

c Local variables

c Begin program

      call defineTSParameters

      cnf1 = cnf
      one_over_dt1 = one_over_dt

c End program

      end subroutine defineTSParametersRoutine

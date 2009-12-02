c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      implicit none

c Call variables

c Local variables

      integer :: ierr

c Begin program

      call closeRecordFile

      deallocate(zeros,vzeros,ones)

      call destroyNonlinearSolver

      call destroyNonlinearFunction   !External

      call destroyGrid(gv%gparams)

      deallocate(gv%bc_seq,STAT=ierr)

      deallocate(gv,STAT=ierr)

c End program

      end subroutine fortranDestroy

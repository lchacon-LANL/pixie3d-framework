c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      use graphics

      implicit none

c Call variables

c Local variables

      integer :: ierr

c Begin program

      deallocate(zeros,vzeros,ones)

      call deallocateStructures

      call deallocateApplicationVariables   !External

      call destroyGrid(gv%gparams)

cc      deallocate(gv%bc_seq)

      deallocate(gv,STAT=ierr)

c End program

      end subroutine fortranDestroy

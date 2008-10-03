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

      call deallocateDerivedType(u_ic)

      call deallocateApplicationVariables   !External

      call destroyGrid(gv%gparams)

      deallocate(gv%bc_seq,STAT=ierr)

      deallocate(gv,STAT=ierr)

c End program

      end subroutine fortranDestroy

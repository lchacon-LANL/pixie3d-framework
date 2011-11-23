c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      use graphics

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr,igr

c Begin program

      call deallocateStructures

      call deallocateDerivedType(u_ic)

      call deallocateApplicationVariables   !External

#if defined(MK)
      if (mk_grid) call MK_dealloc
#endif

      call deallocateGlobalVar(gv)

c End program

      end subroutine fortranDestroy

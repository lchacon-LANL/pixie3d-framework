c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      use ts_setup

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr

c Begin program

c Free memory

      call deallocateStructures

      call destroyTS

      call destroyNonlinearFunction   !External

#if defined(MK)
      if (mk_grid) call MK_dealloc
#endif

      call deallocateGlobalVar(gv)

c Close record files

#if !defined(adios)
      close(urecord)
#else
      call ADIOS_finalize(my_rank, adios_err)
#endif

c Destroy PETSc

#if defined(petsc)
      call PetscLogStagePop(ierr)

      call PetscFinalize(ierr)
#endif

c End program

      end subroutine fortranDestroy

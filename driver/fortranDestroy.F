c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      use ts_setup

      implicit none

c Call variables

c Local variables

      integer :: ierr,igr

c Begin program

c Free memory

      call deallocateStructures

      call destroyTS

      call destroyNonlinearFunction   !External

      call deallocateGlobalVar(gv)

c Close record files

#if !defined(adios)
      close(urecord)
#else
      call adios_finalize(my_rank, adios_err)
#endif

c Destroy PETSc

#if defined(petsc)
      call PetscLogStagePop(ierr)

      call PetscFinalize(ierr)
#endif

c End program

      end subroutine fortranDestroy

c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      implicit none

c Call variables

c Local variables

      integer :: ierr,igr

c Begin program

c Free memory

      deallocate(zeros,vzeros,ones)

      call destroyNonlinearSolver

      call destroyNonlinearFunction   !External

      call destroyGrid(gv%gparams)

      if (associated(gv%bc_grp)) then
        do igr=1,size(gv%bc_grp)
          deallocate(gv%bc_grp(igr)%bc_seq,STAT=ierr)
        enddo
        deallocate(gv%bc_grp,STAT=ierr)
      endif

      deallocate(gv,STAT=ierr)

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

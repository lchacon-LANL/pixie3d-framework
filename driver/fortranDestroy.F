c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid_mpi

      implicit none

c Call variables

c Local variables

      integer :: ierr

c Begin program

      call closeRecordFile

      deallocate(gv,STAT=ierr)

cc      call deallocateDerivedType(u_0)
cc      call deallocateDerivedType(u_n)
cc      call deallocateDerivedType(u_nm)
cc      call deallocateDerivedType(u_np)
cc      call deallocateDerivedType(u_graph)

#if defined(petsc)
      call destroyFortranMPI
#endif

c End program

      end subroutine fortranDestroy

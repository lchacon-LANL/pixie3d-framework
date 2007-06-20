c fortrandestroy
c ######################################################################
      subroutine fortrandestroy

      use variables

      use iosetup

      use grid_mpi

      implicit none

c Call variables

c Local variables

c Begin program

      call closeRecordFile

      call deallocateDerivedType(u_0)
      call deallocateDerivedType(u_n)
      call deallocateDerivedType(u_nm)
      call deallocateDerivedType(u_np)
      call deallocateDerivedType(u_graph)

#if defined(petsc)
      call destroyFortranMPI
#endif

c End program

      end subroutine fortrandestroy

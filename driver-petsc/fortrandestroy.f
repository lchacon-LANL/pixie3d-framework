c fortrandestroy
c ######################################################################
      subroutine fortrandestroy

      use variables

      use iosetup

      implicit none

c Call variables

c Local variables

c Begin program

      call closeRecordFile

      call deallocateDerivedType(u_0)
      call deallocateDerivedType(u_n)
      call deallocateDerivedType(u_graph)

      call destroyDA

c End program

      end subroutine fortrandestroy

c fortrandestroy
c ######################################################################
      subroutine fortrandestroy

      use variables

      use iosetup

      implicit none

c Call variables

c Local variables

c Begin program

cc      call close(urecord)

      call deallocateStructures

cc      call destroyFortranMPI

c End program

      end subroutine fortrandestroy

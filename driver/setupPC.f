c setupPC
c####################################################################
      subroutine setupPC(ntot,x)

c--------------------------------------------------------------------
c     Interface routine for setupPreconditioner
c--------------------------------------------------------------------

      use grid

      use variables

      implicit none

c Call variables

      integer(4) :: ntot
      real(8)    :: x(ntot)

c Local variables

      type(var_array) :: varray

c Begin program

c Unpack vector x

      varray = x  !This allocates varray; overloaded assignment

      call setupPreconditioner(varray)

      call deallocateDerivedType(varray)

      end subroutine setupPC

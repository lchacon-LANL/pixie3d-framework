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

      integer    :: ntot
      real(8)    :: x(ntot)

c Local variables

      type(var_array),pointer :: varray

c Interface

      INTERFACE
        subroutine setupPreconditioner(varray)
        use variable_setup
        type(var_array),pointer :: varray
        end subroutine setupPreconditioner
      END INTERFACE

c Begin program

c Unpack vector x

cc      varray = x  !This allocates varray; overloaded assignment
      call mapVectorToStructure(varray,x)

      call setupPreconditioner(varray)

      call deallocateDerivedType(varray)

      end subroutine setupPC

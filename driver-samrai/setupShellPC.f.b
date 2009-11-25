c setupShellPC
c######################################################################
      subroutine setupShellPC(array,imin,imax,jmin,jmax,kmin,kmax,nl_it)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use newtongm

      implicit none

c Call variables

      integer         :: imin,imax,jmin,jmax,kmin,kmax,nl_it

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer         :: iminl,imaxl,jminl,jmaxl,kminl,kmaxl,ieq
      type(var_array),pointer :: varray

c Interface

      INTERFACE
        subroutine setupPreconditioner(varray)
        use variable_setup
        type(var_array),pointer :: varray
        end subroutine setupPreconditioner
      END INTERFACE

c Begin program

      call fromGlobalToLocalLimits(imin,jmin,kmin,iminl,jminl,kminl
     $                            ,1,1,1)
      call fromGlobalToLocalLimits(imax,jmax,kmax,imaxl,jmaxl,kmaxl
     $                            ,1,1,1)

c Map petsc array

      call initializeDerivedType(varray)

      do ieq=1,neqd
        varray%array_var(ieq)
     .       %array(iminl:imaxl,jminl:jmaxl,kminl:kmaxl)
     .      = array(imin :imax ,jmin :jmax ,kmin :kmax )%var(ieq)
      enddo

c Store nonlinear iteration (for PC preprocessing and diagnostics)

      jit = nl_it + 1

c Setup parallel BC flags to indicate PETSc provides BCs for varray

      call setup_petsc_BC

c Call PC fortran setup routine

      call setupPreconditioner(varray)

c Deallocate memory

      call deallocateDerivedType(varray)

c End program

      end subroutine setupShellPC

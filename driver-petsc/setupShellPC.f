c setupShellPC
c######################################################################
      subroutine setupShellPC(array,imin,imax,jmin,jmax,kmin,kmax)

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

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer(4)      :: iminl,imaxl,jminl,jmaxl,kminl,kmaxl,ieq
      type(var_array) :: varray

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

c Call PC fortran setup routine

      call setupPreconditioner(varray)

c Deallocate memory

      call deallocateDerivedType(varray)

c End program

      end subroutine setupShellPC

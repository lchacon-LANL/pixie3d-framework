c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(x,f,ilo,ihi,jlo,jhi,klo,khi
     .                           ,ilogc,ihigc,jlogc,jhigc,klogc,khigc)
c--------------------------------------------------------------------
c     Calculates nonlinear residuals, of the form:
c             dt Ui + Fi(Uj) = 0
c--------------------------------------------------------------------

      use grid

      use variables

      use timeStepping

      implicit none

c Call variables

      integer(4)  ::  ilo,ihi,jlo,jhi,klo,khi
     .               ,ilogc,ihigc,jlogc,jhigc,klogc,khigc

      type(petsc_var) :: x(ilogc:ihigc,jlogc:jhigc,klogc:khigc)

      type(petsc_var) :: f(ilo:ihi,jlo:jhi,klo:khi)

c Local variables

      integer(4) :: i,j,k,ieq,ii

      real(8)    :: dudt(neqd),cnf(neqd),one_over_dt(neqd),ff(ntotd)

      type (var_array)  :: varray

      type(petsc_array) :: petscarray

c Begin program

cc      write (*,*) ilo,ihi,jlo,jhi,klo,khi
cc     .               ,ilogc,ihigc,jlogc,jhigc,klogc,khigc
cc      stop
      call allocatePetscType(petscarray)

      petscarray%array(ilogc:ihigc,jlogc:jhigc,klogc:khigc)
     .             = x(ilogc:ihigc,jlogc:jhigc,klogc:khigc)

c Unpack petsc array x

      varray = petscarray   !Overloaded assignment

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,ff)

c Assign ff to f

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))
            do ieq=1,neqd

              f(i,j,k)%var(ieq) = ff  (ii+ieq)

            enddo
          enddo
        enddo
      enddo
      
c Deallocate structures

      call deallocateDerivedType(varray)
      call deallocatePetscType(petscarray)

c End program

      end subroutine

c evaluateNonlinearFunction
c####################################################################
      subroutine evaluateNonlinearFunction(varray,fi)
c--------------------------------------------------------------------
c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c--------------------------------------------------------------------

      use parameters

      use variable_setup

      implicit none

c Call variables

      real(8)          :: fi(ntotd)
      type (var_array) :: varray

c Local variables

      integer(4)       :: i,j,k,ii

c Interfaces

      INTERFACE
         subroutine setupNonlinearFunction(varray)
           use variable_setup
           type (var_array),target :: varray
         end subroutine setupNonlinearFunction
      END INTERFACE

c Begin program

c Prepare auxiliar quantities

      call setupNonlinearFunction(varray)

c Store function evaluation

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd
            ii = neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))
            call nonlinearRHS(i,j,k,varray,fi(ii+1:ii+neqd))
          enddo
        enddo
      enddo

c Deallocate variables

      call killNonlinearFunction

c End program

      end subroutine

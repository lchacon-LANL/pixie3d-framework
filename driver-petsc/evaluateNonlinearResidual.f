c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(x,f
     $                    ,imin  ,imax  ,jmin  ,jmax  ,kmin  ,kmax
     .                    ,imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc)
c--------------------------------------------------------------------
c     Calculates Fi(Uj) in equations of the form: dt Ui + Fi(Uj) = 0
c--------------------------------------------------------------------

      use grid

      use variables

      use timeStepping

      implicit none

c Call variables

      integer(4)  ::  imin,imax,jmin,jmax,kmin,kmax
     .               ,imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc

      type(petsc_var) :: x(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)

      type(petsc_var) :: f(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer(4) :: i,j,k,il,jl,kl,ieq,ii
     .             ,imingcl,imaxgcl,jmingcl,jmaxgcl,kmingcl,kmaxgcl

      real(8)    :: dudt(neqd),ff(ntotd)

      type(var_array)   :: varray

c Begin program

c Find local limits

      call fromGlobalToLocalLimits(imingc ,jmingc ,kmingc
     $                            ,imingcl,jmingcl,kmingcl,1,1,1)
      call fromGlobalToLocalLimits(imaxgc ,jmaxgc ,kmaxgc
     $                            ,imaxgcl,jmaxgcl,kmaxgcl,1,1,1)

c Unpack petsc array

      call allocateDerivedType(varray)

      do ieq=1,neqd
        varray%array_var(ieq)
     .       %array(imingcl:imaxgcl,jmingcl:jmaxgcl,kmingcl:kmaxgcl)
     .       = x(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)%var(ieq)
      enddo

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,ff)

c Assign ff (vector) to f (PETSc array)

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            call fromGlobalToLocalLimits(i,j,k,il,jl,kl,1,1,1)
            ii = vecPos(neqd,il,jl,kl,1,1,1)
            do ieq=1,neqd
              f(i,j,k)%var(ieq) = ff(ii+ieq)
            enddo
          enddo
        enddo
      enddo
      
c Deallocate structures

      call deallocateDerivedType(varray)

c End program

      end subroutine evaluateNonlinearResidual

c evaluateNonlinearFunction
c####################################################################
      subroutine evaluateNonlinearFunction(varray,fi)
c--------------------------------------------------------------------
c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c--------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      implicit none

c Call variables

      real(8)          :: fi(ntotd)
      type (var_array) :: varray

c Local variables

      integer(4)       :: i,j,k,ii

c Diag

      integer(4)       :: ieq

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

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k,1,1,1)
            call nonlinearRHS(i,j,k,varray,fi(ii+1:ii+neqd))
          enddo
        enddo
      enddo

c Deallocate variables

      call killNonlinearFunction

c End program

      end subroutine

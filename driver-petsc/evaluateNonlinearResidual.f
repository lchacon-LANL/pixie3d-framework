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

      real(8)    :: dudt(neqd),cnf(neqd),one_over_dt(neqd),ff(ntotd)

      type(var_array)   :: varray

      type(petsc_array) :: petscarray

c Begin program

cc      write (*,*) 'proc',my_rank,'lim',imin,imax,jmin,jmax,kmin,kmax
cc      write (*,*) my_rank,imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc

      call allocatePetscType(petscarray)

      call fromGlobalToLocalLimits(imingc ,jmingc ,kmingc
     $                            ,imingcl,jmingcl,kmingcl)
      call fromGlobalToLocalLimits(imaxgc ,jmaxgc ,kmaxgc
     $                            ,imaxgcl,jmaxgcl,kmaxgcl)

      petscarray%array(imingcl:imaxgcl,jmingcl:jmaxgcl,kmingcl:kmaxgcl)
     .             = x(imingc :imaxgc ,jmingc :jmaxgc ,kmingc :kmaxgc )

c Unpack petsc array

      varray = petscarray   !Overloaded assignment

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,ff)

c Assign ff (vector) to f (PETSc array)

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            call fromGlobalToLocalLimits(i,j,k,il,jl,kl)
            ii = vecPos(neqd,il,jl,kl)
            do ieq=1,neqd
              f(i,j,k)%var(ieq) = ff(ii+ieq)
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

      use grid

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

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k)
            call nonlinearRHS(i,j,k,varray,fi(ii+1:ii+neqd))
cc            fi(ii+1:ii+neqd) = cos(grid_params%xx(i+1)
cc     .                            *grid_params%yy(j+1))
          enddo
        enddo
      enddo

c Deallocate variables

      call killNonlinearFunction

c End program

      end subroutine

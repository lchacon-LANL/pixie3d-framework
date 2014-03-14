c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(patch_var,n_elem,fsrc,f)
c--------------------------------------------------------------------
c     Calculates Fi(Uj) in equations of the form: dt Ui = Fi(Uj)
c     Note: this function assumes the fluxes have already been set
c--------------------------------------------------------------------

      use parameters
      use grid
      use variables
      use timeStepping
      implicit none

c Call variables
      integer :: n_elem
      type(patch), TARGET :: patch_var
      real(8),dimension(n_elem) :: fsrc
      type(var_array) :: f

c Local variables
      integer    :: i,j,k,il,jl,kl,ieq,ii,ig,jg,kg
      real(8)    :: dudt(neqd),ff(n_elem),dvol
      type(var_array),pointer :: varray

c Interfaces
      INTERFACE
c        Initialize the fluxes
         subroutine evaluateNonlinearResidualFlux(patch_var)
           use parameters
           use grid
           use variables
           implicit none
      type(patch), TARGET :: patch_var
         end subroutine evaluateNonlinearResidualFlux
      END INTERFACE

      INTERFACE
c        evaluate the nonlinearRHS
         subroutine evaluateNonlinearFunction(varray,ff)
           use parameters
           use variable_setup
           use grid
           implicit none
           type (var_array),pointer :: varray
           real(8)          :: ff(ntotd)
         end subroutine evaluateNonlinearFunction
      END INTERFACE

c Begin program
      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions
      varray => patch_var%u_n

c Prepare auxiliar quantities (fluxes)
c      call evaluateNonlinearResidualFlux(patch_var)

c Store function evaluation
      call evaluateNonlinearFunction(varray,ff)

c Evaluate nonlinear function Fi(Uj) at time level (n+1)
      do k = klog,khig
        do j = jlog,jhig
          do i = ilog,ihig
            call fromGlobalToLocalLimits(i,j,k,il,jl,kl,1,1,1)
            ii = vecPos(neqd,il,jl,kl,1,1,1)
            do ieq=1,neqd
              f%array_var(ieq)%array(i,j,k)=fsrc(ii+ieq)-ff(ii+ieq)
            enddo
          enddo
        enddo
      enddo

c Deallocate structures

c End program
      end subroutine evaluateNonlinearResidual




c evaluateNonlinearFunction
c####################################################################
      subroutine evaluateNonlinearFunction(varray,fi)
c--------------------------------------------------------------------
c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c     Note: this function assumes the fluxes have already been set
c--------------------------------------------------------------------

      use parameters
      use variable_setup
      use grid
      implicit none

c Call variables
      type (var_array),pointer :: varray
      real(8)          :: fi(ntotd)

c Local variables
      integer          :: i,j,k,ii

c Interfaces
      INTERFACE
         subroutine nonlinearRHS(i,j,k,igx,igy,igz,varray,ff)
           use variable_setup
           use timeStepping
           implicit none
           
           real(8) :: ff(neqd)
           integer :: i,j,k,igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine nonlinearRHS
      END INTERFACE

c Store function evaluation
      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k,1,1,1)
            call nonlinearRHS(i,j,k,1,1,1,varray,fi(ii+1:ii+neqd))
          enddo
        enddo
      enddo

c End program
      end subroutine evaluateNonlinearFunction




c evaluateNonlinearResidualFlux
c####################################################################
      subroutine evaluateNonlinearResidualFlux(patch_var)
c--------------------------------------------------------------------
c     Calculates the fluxes needed for the evaluateNonlinearResidual
c--------------------------------------------------------------------

      use parameters
      use grid
      use variables

      implicit none

c Call variables
      type(patch), TARGET :: patch_var

c Interfaces
      INTERFACE
         subroutine setupNonlinearFunction(igx,igy,igz,varray)
           use variable_setup
           integer :: igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine setupNonlinearFunction
      END INTERFACE

c Begin program
      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

c Prepare auxiliar quantities (flux and flux_src)
      call setupNonlinearFunction(1,1,1,patch_var%u_n)

c Deallocate variables
      call killNonlinearFunction

c End program
      end subroutine evaluateNonlinearResidualFlux



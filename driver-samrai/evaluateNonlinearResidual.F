c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(patch_var,n_elem,fsrc,f)
c--------------------------------------------------------------------
c     Calculates Fi(Uj) in equations of the form: dt Ui = Fi(Uj)
c--------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

c      use auxiliaryVariables

c      use grid_aliases

      implicit none

c Call variables

      integer :: n_elem
      type(patch), TARGET :: patch_var
      real(8),dimension(n_elem) :: fsrc
      type(var_array) :: f

c Local variables

      integer     ::  imin,imax,jmin,jmax,kmin,kmax
     .               ,imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc

      integer    :: i,j,k,il,jl,kl,ieq,ii,ig,jg,kg
     .             ,imingcl,imaxgcl,jmingcl,jmaxgcl,kmingcl,kmaxgcl

      real(8)    :: dudt(neqd),ff(n_elem),dvol

      type(var_array),pointer :: varray

c Interfaces

      INTERFACE
         subroutine evaluateNonlinearFunction(varray,fi)
           use parameters
           use variable_setup
           use grid
           real(8) :: fi(ntotd)
           type (var_array),pointer :: varray
         end subroutine evaluateNonlinearFunction
      END INTERFACE

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => patch_var%u_n

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,ff)

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

c Diag

      integer          :: ieq

c Interfaces

      INTERFACE
         subroutine setupNonlinearFunction(igx,igy,igz,varray)
           use variable_setup
           integer :: igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine setupNonlinearFunction
      END INTERFACE

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

c Begin program

c     Prepare auxiliar quantities

      call setupNonlinearFunction(1,1,1,varray)

c Store function evaluation

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k,1,1,1)
            call nonlinearRHS(i,j,k,1,1,1,varray,fi(ii+1:ii+neqd))
          enddo
        enddo
      enddo

c Deallocate variables

      call killNonlinearFunction

c End program

      end subroutine

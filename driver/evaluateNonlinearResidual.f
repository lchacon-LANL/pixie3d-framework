c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(ntot,x,f)
c--------------------------------------------------------------------
c     Calculates nonlinear residuals, of the form:
c             dt Ui + Fi(Uj) = 0
c--------------------------------------------------------------------

      use grid

      use variables

      use timeStepping

      implicit none

c Call variables

      integer    :: ntot
      real(8)    :: x(ntot),f(ntot)

c Local variables

      integer    :: i,j,k,ieq,ii,ig,jg,kg
      real(8)    :: dvol

      type(var_array),pointer :: varray => null()

c Interfaces

      INTERFACE
        subroutine evaluateNonlinearFunction(varray,fi)
        use parameters
        use variable_setup
        real(8)          :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluateNonlinearFunction
      END INTERFACE

c Begin program

c Unpack vector x

      call mapVectorToStructure(varray,x)

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,f)

c Calculate residuals

      call defineTSParameters

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi

            call getMGmap(i,j,k,1,1,1,ig,jg,kg)

            ii = vecPos(neqd,i,j,k,1,1,1)

            do ieq=1,neqd
              f(ii+ieq) =
     .              (bdfp (ieq)*varray%array_var(ieq)%array(i,j,k)
     .              +bdfn (ieq)*u_n   %array_var(ieq)%array(i,j,k)
     .              +bdfnm(ieq)*u_nm  %array_var(ieq)%array(i,j,k))
     .              *one_over_dt(ieq)
     .              + ((1d0-cnf(ieq))*f   (ii+ieq)
     .              +       cnf(ieq) *fold(ii+ieq)
     .              -                 fsrc(ii+ieq))
            enddo

            if (vol_wgt) then
              if (.not.checkMapDatabase()) then
              !Do not include Jacobian here to allow for moving grid cases
                dvol = grid_params%dxh(ig)
     .                *grid_params%dyh(jg)
     .                *grid_params%dzh(kg)
                f(ii+1:ii+neqd) = f(ii+1:ii+neqd)*dvol
              else
              !Fixed grid case
                f(ii+1:ii+neqd) = f(ii+1:ii+neqd)
     .                           *gmetric%grid(1)%dvol(i,j,k)
              endif
            endif

          enddo
        enddo
      enddo

c Deallocate structure

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

      implicit none

c Call variables

      real(8)          :: fi(ntotd)
      type(var_array),pointer :: varray

c Local variables

      integer          :: i,j,k,ii

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
           real(8) :: ff(neqd)
           integer :: i,j,k,igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine nonlinearRHS
      END INTERFACE

c Begin program

c Setup parallel BC flags to indicate BCs require communication

      call setup_petsc_BC

c Prepare auxiliar quantities

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

      end subroutine evaluateNonlinearFunction

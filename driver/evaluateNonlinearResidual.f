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

      integer  :: ntot
      real*8      x(ntot),f(ntot)

c Local variables

      integer(4) :: i,j,k,ieq,ii

      real(8)    :: dudt(neqd),cnf(neqd),one_over_dt(neqd)

      type (var_array) :: varray

c Begin program

c Unpack vector x

      varray = x         !Overloaded assignment

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,f)

c Calculate residuals

      call defineTSParameters(cnf,one_over_dt)
      
      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd

            ii = neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))

            do ieq=1,neqd
cc              dudt(ieq) = volume(i,j,k,1,1,1)*one_over_dt(ieq)
cc     .                   *( varray%array_var(ieq)%array(i,j,k)
cc     .                     -   u_n%array_var(ieq)%array(i,j,k) )
cc              f(ii+ieq) = dudt(ieq) + (1.-cnf(ieq))*f   (ii+ieq)
cc     .                              +     cnf(ieq) *fold(ii+ieq)
cc     .                              -               fsrc(ii+ieq)

              f(ii+ieq) = (varray%array_var(ieq)%array(i,j,k)
     .                    -   u_n%array_var(ieq)%array(i,j,k))
     .                    *one_over_dt(ieq)
     .                    + ((1.-cnf(ieq))*f   (ii+ieq)
     .                    +      cnf(ieq) *fold(ii+ieq)
     .                    -                fsrc(ii+ieq))
            enddo

          enddo
        enddo
      enddo

c Deallocate structure

      call deallocateDerivedType(varray)

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

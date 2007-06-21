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

      integer     ::  imin,imax,jmin,jmax,kmin,kmax
     .               ,imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc

      type(petsc_var) :: x(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)

      type(petsc_var) :: f(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer    :: i,j,k,il,jl,kl,ieq,ii,ig,jg,kg
     .             ,imingcl,imaxgcl,jmingcl,jmaxgcl,kmingcl,kmaxgcl

      real(8)    :: dudt(neqd),ff(ntotd),dvol

      type(var_array)   :: varray

c Begin program

c Find local limits

      call fromGlobalToLocalLimits(imingc ,jmingc ,kmingc
     $                            ,imingcl,jmingcl,kmingcl,1,1,1)
      call fromGlobalToLocalLimits(imaxgc ,jmaxgc ,kmaxgc
     $                            ,imaxgcl,jmaxgcl,kmaxgcl,1,1,1)

c Unpack petsc array

      call initializeDerivedType(varray)

      do ieq=1,neqd
        varray%array_var(ieq)
     .       %array(imingcl:imaxgcl,jmingcl:jmaxgcl,kmingcl:kmaxgcl)
     .       = x(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)%var(ieq)
      enddo

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(varray,ff)

c Assign ff (vector) to f (PETSc array)

      call defineTSParameters

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            call fromGlobalToLocalLimits(i,j,k,il,jl,kl,1,1,1)
            ii = vecPos(neqd,il,jl,kl,1,1,1)

            do ieq=1,neqd
              f(i,j,k)%var(ieq)=
     .             (bdfp (ieq)*varray%array_var(ieq)%array(il,jl,kl)
     .             +bdfn (ieq)*u_n   %array_var(ieq)%array(il,jl,kl)
     .             +bdfnm(ieq)*u_nm  %array_var(ieq)%array(il,jl,kl))
     .             *one_over_dt(ieq)
     .             + ((1d0-cnf(ieq))*ff  (ii+ieq)
     .             +       cnf(ieq) *fold(ii+ieq)
     .             -                 fsrc(ii+ieq))
            enddo

            if (vol_wgt) then
              if (.not.checkAnalMapDatabase()) then
              !Do not include Jacobian here to allow for moving grid cases

                call getMGmap(il,jl,kl,1,1,1,ig,jg,kg)

                dvol = grid_params%dxh(ig)
     .                *grid_params%dyh(jg)
     .                *grid_params%dzh(kg)

                f(i,j,k)%var(:) = f(i,j,k)%var(:)*dvol
              else

              !Fixed grid case
                f(i,j,k)%var(:) = f(i,j,k)%var(:)
     .                           *gmetric%grid(1)%dvol(il,jl,kl)
              endif
            endif

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

      integer          :: i,j,k,ii

c Diag

      integer          :: ieq

c Interfaces

      INTERFACE
         subroutine setupNonlinearFunction(varray)
           use variable_setup
           type (var_array),target :: varray
         end subroutine setupNonlinearFunction
      END INTERFACE

c Begin program

c Setup parallel BC flags to indicate PETSc provides BCs

      call setup_petsc_BC

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

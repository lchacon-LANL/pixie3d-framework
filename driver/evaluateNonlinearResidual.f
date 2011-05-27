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

      use mk

      implicit none

c Call variables

      integer    :: ntot
      real(8)    :: x(ntot),f(ntot)

c Local variables

      integer    :: i,j,k,ieq,ii,ig,jg,kg
      real(8)    :: dvol,src(neqd),dfdt(neqd)

      type(var_array),pointer :: varray => null()

      real(8),pointer,dimension(:,:,:) :: jac

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

c Define time-step parameters

      call defineTSParameters

      if (mk_grid) then
        if (mk_relax_init_grid) then
          cnf = 0d0
          bdfp = 0d0
          bdfn = 0d0
          bdfnm = 0d0
          one_over_dt = 0d0
        else
          cnf  (neqd) = -tau/dt
          bdfp (neqd) = 0d0
          bdfn (neqd) = 0d0
          bdfnm(neqd) = 0d0
          one_over_dt(neqd) = 0d0
        endif
      endif

c Calculate residuals

      jac => gmetric%grid(1)%jac

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi

            call getMGmap(i,j,k,1,1,1,ig,jg,kg)

            ii = vecPos(neqd,i,j,k,1,1,1)

            if (source.and.mk_grid) then
              if (mk_relax_init_grid) then
                src = 0d0
              else
                src = MK_eval_src(1,i,j,k) !Interpolate source on current MK mesh
              endif
            else
              src = fsrc(ii+1:ii+neqd)
            endif

            if (mk_grid.and.(.not.mk_nc)) then
              do ieq=1,neqd
                dfdt(ieq) = one_over_dt(ieq)*
     .        (bdfp (ieq)*jac(i,j,k)*varray%array_var(ieq)%array(i,j,k)
     .        +bdfn (ieq)*jacn (i,j,k)*u_n %array_var(ieq)%array(i,j,k)
     .        +bdfnm(ieq)*jacnm(i,j,k)*u_nm%array_var(ieq)%array(i,j,k))
              enddo
            else
              do ieq=1,neqd
                dfdt(ieq) = one_over_dt(ieq)*
     .              (bdfp (ieq)*varray%array_var(ieq)%array(i,j,k)
     .              +bdfn (ieq)*u_n   %array_var(ieq)%array(i,j,k)
     .              +bdfnm(ieq)*u_nm  %array_var(ieq)%array(i,j,k))
              enddo
            endif

            do ieq=1,neqd
              f(ii+ieq) = dfdt(ieq)
     .                  + ((1d0-cnf(ieq))*f   (ii+ieq)
     .                  +       cnf(ieq) *fold(ii+ieq)
     .                  -                  src(   ieq))
            enddo

            if (vol_wgt) then
              if (mk_grid) then
              !Moving grid case
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

      nullify(jac)

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

      use mk

      implicit none

c Call variables

      real(8) :: fi(ntotd)
      type(var_array),pointer :: varray

c Local variables

      integer :: i,j,k,ii,ieq,nx,ny,nz

c Interfaces

      INTERFACE
         subroutine setup_NLF(igx,varray)
           use variable_setup
           integer :: igx
           type(var_array),pointer :: varray
         end subroutine setup_NLF
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

      call setup_NLF(1,varray)

c Store function evaluation

      nx = ihi-ilo+1
      ny = jhi-jlo+1
      nz = khi-klo+1

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k,1,1,1)
            if (mk_relax_init_grid) then
              do ieq=1,neqd-1
                fi(ii+ieq) = varray%array_var(ieq)%array(i,j,k)
              enddo
              fi(ii+1:ii+neqd) = fi(ii+1:ii+neqd) - MK_eval_equ(1,i,j,k)
            else
              call nonlinearRHS(i,j,k,1,1,1,varray,fi(ii+1:ii+neqd))
            endif

            !Compute MK residual, add grid velocity,
            !and multiply by jacobian factor (if conservative)
            if (mk_grid) then
              if (.not.mk_relax_init_grid) then
                do ieq=1,neqd-1
                  if (one_over_dt(ieq) == 0d0) cycle
                  fi(ii+ieq) = fi(ii+ieq)
     .                       + flx_advec(i,j,k,nx,ny,nz,1,1,1,gvel
     .                                  ,varray%array_var(ieq)%array
     .                                  ,mk_advect,zip_vel=.false.
     .                                  ,conserv=.not.mk_nc)
                enddo

                if (.not.mk_nc) then  !Non-conservative MK
                  fi(ii+1:ii+neqd-1) = fi(ii+1:ii+neqd-1)
     .                                *gmetric%grid(1)%jac(i,j,k)
                endif
              endif

              fi(ii+neqd) = MK_residual(i,j,k,1,1,1
     .                                 ,varray%array_var(neqd)%array)
            endif
          enddo
        enddo
      enddo

c Deallocate variables

      call killNonlinearFunction

c End program

      end subroutine evaluateNonlinearFunction

c setup_NLF
c#################################################################
      subroutine setup_NLF(igrid,varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variables

      use mk

      implicit none

c Call variables

      integer :: igrid
      type (var_array),pointer :: varray

c Local variables

c Interfaces

      INTERFACE
         subroutine setupNonlinearFunction(igx,igy,igz,varray)
           use variable_setup
           integer :: igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine setupNonlinearFunction
      END INTERFACE

c Begin program

c Set up MK grid (MK variable is last equation)

      if(mk_grid) call MK_setup_grid(igrid,varray%array_var(neqd)%array)

c Call application setup

      call setupNonlinearFunction(igrid,igrid,igrid,varray)

c Setup MK monitor function

      if(mk_grid) call MK_get_mon(igrid)

c End program

      end subroutine setup_NLF

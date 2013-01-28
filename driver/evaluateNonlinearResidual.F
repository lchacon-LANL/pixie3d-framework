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

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ntot
      real(8) :: x(ntot),f(ntot)

c Local variables

      integer :: i,j,k,ieq,ii,ig,jg,kg,igrid
      real(8) :: dvol,src(neqd),dfdt(neqd)

      real(8),pointer,dimension(:,:,:) :: jac

c Begin program

      igrid = 1

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(ntot,x,f)

c Calculate residuals

      jac => gmetric%grid(igrid)%jac

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ii = vecPos(neqd,i,j,k,igrid,igrid,igrid)

#if defined(MK)
            if (source.and.mk_grid) then
              if (relax_init_grid) then
                src = 0d0
              else
                src = MK_eval_src(igrid,i,j,k) !Interpolate source on current MK mesh
              endif
            else
              src = fsrc(ii+1:ii+neqd)
            endif
#else
            src = fsrc(ii+1:ii+neqd)
#endif

#if defined(MK)
            if (mk_grid.and.(.not.mk_nc)) then
              do ieq=1,neqd
                dfdt(ieq) = one_over_dt(ieq)*
     .        (bdfp (ieq)*jac(i,j,k)*varray%array_var(ieq)%array(i,j,k)
     .        +bdfn (ieq)*jacn (i,j,k)*u_n %array_var(ieq)%array(i,j,k)
     .        +bdfnm(ieq)*jacnm(i,j,k)*u_nm%array_var(ieq)%array(i,j,k))
              enddo
            else
              dfdt=cnp*(x(ii+1:ii+neqd)-x_old(ii+1:ii+neqd))*one_over_dt
            endif
#else
            dfdt = cnp*(x(ii+1:ii+neqd)-x_old(ii+1:ii+neqd))*one_over_dt
#endif

            f(ii+1:ii+neqd) = dfdt + (1d0-cnf)*f   (ii+1:ii+neqd)
     .                             +      cnf *fold(ii+1:ii+neqd) - src

            if (vol_wgt) then
#if defined(MK)
              if (mk_grid.and.(.not.mk_nc)) then
              !Moving grid case
                dvol = gmetric%grid(igrid)%dlvol(i,j,k) !Logical volume
              else
              !Fixed grid case
                dvol = gmetric%grid(igrid)%dvol(i,j,k)
              endif
#else
              dvol = gmetric%grid(igrid)%dlvol(i,j,k)  !Logical volume
#endif
              f(ii+1:ii+neqd) = f(ii+1:ii+neqd)*dvol
            endif

          enddo
        enddo
      enddo

      nullify(jac)

c End program

      end subroutine evaluateNonlinearResidual

c evaluateNonlinearFunction
c####################################################################
      subroutine evaluateNonlinearFunction(ntot,xi,fi)
c--------------------------------------------------------------------
c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c--------------------------------------------------------------------

      use parameters

      use variable_setup

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ntot
      real(8) :: xi(ntot),fi(ntot)

c Local variables

      integer :: igr

      type(var_array),pointer :: varray => null()

c Interfaces

      INTERFACE
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use variable_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c Begin program

      igr = 1

c Prepare auxiliar quantities

      call mapVectorToStructure(varray,xi,igr=igr)

c Evaluate nonlinear function

      call evaluate_NLF(igr,varray,fi)

c Deallocate variables

      call deallocateDerivedType(varray)

c End program

      end subroutine evaluateNonlinearFunction

c evaluate_NLF
c####################################################################
      subroutine evaluate_NLF(igr,varray,fi)
c--------------------------------------------------------------------
c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c--------------------------------------------------------------------

      use parameters

      use variable_setup

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: igr
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

c Define time-step parameters

      call defineTSParameters

#if defined(MK)
      if (mk_grid) then
        if (relax_init_grid) then
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
#endif

c Setup parallel BC flags to indicate BCs require communication

      call setup_petsc_BC

c Prepare auxiliar quantities

      call setup_NLF(igr,varray)

c Store function evaluation

      nx = grid_params%nxv(igr)
      ny = grid_params%nyv(igr)
      nz = grid_params%nzv(igr)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii = vecPos(neqd,i,j,k,igr,igr,igr)

#if defined(MK)
            !Compute MK residual, add grid velocity,
            !and multiply by jacobian factor (if conservative)
            if (mk_grid) then

              if (relax_init_grid) then

                !Enforce equilibrium
                do ieq=1,neqd-1
                  fi(ii+ieq) = varray%array_var(ieq)%array(i,j,k)
                enddo
                fi(ii+1:ii+neqd)=fi(ii+1:ii+neqd)-MK_eval_equ(igr,i,j,k)

              else
                call nonlinearRHS(i,j,k,igr,igr,igr,varray
     .                           ,fi(ii+1:ii+neqd))

cc                do ieq=1,neqd-1
cc                  if (one_over_dt(ieq) == 0d0) cycle
cc                  fi(ii+ieq) = fi(ii+ieq)
cc     .                       + flx_advec(i,j,k,nx,ny,nz,igr,igr,igr,gvel
cc     .                                 ,varray%array_var(ieq)%array
cccc     .                                 ,mk_advect,zip_vel=.true.
cc     .                                 ,mk_advect,zip_vel=(mk_advect==2)
cc     .                                 ,conserv=.not.mk_nc)
cc                enddo

                if (.not.mk_nc) then  !Conservative MK
                  fi(ii+1:ii+neqd-1) = fi(ii+1:ii+neqd-1)
     .                                *gmetric%grid(igr)%jac(i,j,k)
                endif
              endif

              fi(ii+neqd) = MK_residual(i,j,k,igr
     .                                 ,varray%array_var(neqd)%array)
            else
              call nonlinearRHS(i,j,k,igr,igr,igr,varray
     .                         ,fi(ii+1:ii+neqd))
            endif
#else
            call nonlinearRHS(i,j,k,igr,igr,igr,varray,fi(ii+1:ii+neqd))
#endif

          enddo
        enddo
      enddo

c Deallocate variables

      call killNonlinearFunction

c End program

      end subroutine evaluate_NLF

c setup_NLF
c#################################################################
      subroutine setup_NLF(igrid,varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variables

#if defined(MK)
      use mk
#endif

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

#if defined(MK)
      if(mk_grid) call MK_setup_grid(igrid,varray%array_var(neqd)%array)
#endif

c Call application setup

      call setupNonlinearFunction(igrid,igrid,igrid,varray)

c Setup MK monitor function

#if defined(MK)
      if(mk_grid) call MK_get_mon(igrid)
#endif

c End program

      end subroutine setup_NLF

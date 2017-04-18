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

      use nk, ONLY: pseudo_dt,pdt,nk_conf

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ntot
      real(8) :: x(ntot),f(ntot)

c Local variables

      integer :: i,j,k,ieq,ii,ig,jg,kg,igrid

      real(8),pointer,dimension(:) :: src,xdot

c Begin program

      igrid = 1

      allocate(src(ntot),xdot(ntot))

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      if (pseudo_dt) then
        xdot = (x-nk_conf%xk)/pdt
      else
        xdot = cnp*(x-x_old)/dt
      endif

      call evaluateNonlinearFunction(igrid,ntot,x,f)

c Compute source

      call evaluateSource(igrid,ntot,src)

c Temporal update

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = neqd*vPos(gv%gparams,igrid,(/i,j,k/))

            f(ii+1:ii+neqd) = xdot(ii+1:ii+neqd)*dt*one_over_dt
     .                      + (1d0-cnf)*f   (ii+1:ii+neqd)
     .                      +      cnf *fold(ii+1:ii+neqd)
     .                      - src(ii+1:ii+neqd)
          enddo
        enddo
      enddo

c Scale by volume

      if (vol_wgt)
     .     call scaleVector_by_vol1(gv%gparams,igrid,neqd,f,'mul')

c End program

      deallocate(src,xdot)

      end subroutine evaluateNonlinearResidual

c evaluateNonlinearFunction
c####################################################################
      subroutine evaluateNonlinearFunction(igr,ntot,xi,fi)
c--------------------------------------------------------------------
c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c--------------------------------------------------------------------

      use parameters

      use var_setup

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ntot,igr
      real(8) :: xi(ntot),fi(ntot)

c Local variables

      type(var_array),pointer :: varray => null()

c Interfaces

      INTERFACE
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use var_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c Begin program

c Prepare auxiliar quantities

      call mapVectorToStructure(gv%gparams,varray,xi,igr=igr)

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

      use var_setup

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
           use var_setup
           integer :: igx
           type(var_array),pointer :: varray
         end subroutine setup_NLF
      END INTERFACE

      INTERFACE
         subroutine nonlinearRHS(i,j,k,igx,igy,igz,varray,ff)
           use var_setup
           real(8) :: ff(neqd)
           integer :: i,j,k,igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine nonlinearRHS
      END INTERFACE

c Begin program

c Setup parallel BC flags to indicate BCs require communication

      call setup_petsc_BC(gv%gparams)

c Prepare auxiliary quantities

      call setup_NLF(igr,varray)

c Store function evaluation

      nx = gv%gparams%nxv(igr)
      ny = gv%gparams%nyv(igr)
      nz = gv%gparams%nzv(igr)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
cc            ii = vecPos(neqd,i,j,k,igr,igr,igr)
            ii = neqd*vPos(gv%gparams,igr,(/i,j,k/))

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
           use var_setup
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

c evaluateSource
c####################################################################
      subroutine evaluateSource(igrid,ntot,src)
c--------------------------------------------------------------------
c     Calculates residual source, including future IMEX terms.
c--------------------------------------------------------------------

      use grid

      use variables

      use timeStepping

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer :: ntot,igrid
      real(8) :: src(ntot)

c Local variables

      integer :: i,j,k,ieq,ii,ig,jg,kg

c Begin program

#if defined(MK)
      if (source.and.mk_grid) then

        if (relax_init_grid) then
          src = 0d0
        else
          do k = klo,khi
            do j = jlo,jhi
              do i = ilo,ihi
cc                ii = vecPos(neqd,i,j,k,igrid,igrid,igrid)
                ii = neqd*vPos(gv%gparams,igrid,(/i,j,k/))
                src(ii+1:ii+neqd) = MK_eval_src(igrid,i,j,k) !Interpolate source on current MK mesh
              enddo
            enddo
          enddo
        endif
      else
        src = fsrc
      endif
#else
      src = fsrc
#endif

c End program

      end subroutine evaluateSource
      
c$$$c evaluate_IMEX_NLF
c$$$c####################################################################
c$$$      subroutine evaluate_IMEX_NLF(igr,varray,fi)
c$$$c--------------------------------------------------------------------
c$$$c     Stores evaluation of nonlinear function Fi(Uj) in vector fi.
c$$$c     The Uj's are given in varray.
c$$$c--------------------------------------------------------------------
c$$$
c$$$      use parameters
c$$$
c$$$      use var_setup
c$$$
c$$$#if defined(MK)
c$$$      use mk
c$$$#endif
c$$$
c$$$      implicit none
c$$$
c$$$c Call variables
c$$$
c$$$      integer :: igr
c$$$      real(8) :: fi(ntotd)
c$$$      type(var_array),pointer :: varray
c$$$
c$$$c Local variables
c$$$
c$$$      integer :: i,j,k,ii,ieq,nx,ny,nz
c$$$
c$$$c Interfaces
c$$$
c$$$      INTERFACE
c$$$         subroutine setup_IMEX_NLF(igx,varray)
c$$$           use var_setup
c$$$           integer :: igx
c$$$           type(var_array),pointer :: varray
c$$$         end subroutine setup_IMEX_NLF
c$$$      END INTERFACE
c$$$
c$$$      INTERFACE
c$$$         subroutine IMEX_RHS(i,j,k,igx,igy,igz,varray,ff)
c$$$           use var_setup
c$$$           real(8) :: ff(neqd)
c$$$           integer :: i,j,k,igx,igy,igz
c$$$           type(var_array),pointer :: varray
c$$$         end subroutine IMEX_RHS
c$$$      END INTERFACE
c$$$
c$$$c Begin program
c$$$
c$$$c Setup parallel BC flags to indicate BCs require communication
c$$$
c$$$      call setup_petsc_BC
c$$$
c$$$c Prepare auxiliar quantities
c$$$
c$$$      call setup_IMEX_NLF(igr,varray)
c$$$
c$$$c Store function evaluation
c$$$
c$$$      nx = gv%gparams%nxv(igr)
c$$$      ny = gv%gparams%nyv(igr)
c$$$      nz = gv%gparams%nzv(igr)
c$$$
c$$$      do k = 1,nz
c$$$        do j = 1,ny
c$$$          do i = 1,nx
c$$$            ii = vecPos(neqd,i,j,k,igr,igr,igr)
c$$$
c$$$            call IMEX_RHS(i,j,k,igr,igr,igr,varray,fi(ii+1:ii+neqd))
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c Deallocate variables
c$$$
c$$$      call kill_IMEX_NLF
c$$$
c$$$c End program
c$$$
c$$$      end subroutine evaluate_IMEX_NLF

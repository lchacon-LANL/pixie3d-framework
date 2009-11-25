c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(patch_var,n_elem,fsrc,f)
c--------------------------------------------------------------------
c     Calculates Fi(Uj) in equations of the form: dt Ui + Fi(Uj) = 0
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

c      INTERFACE
c         subroutine evaluateNonlinearFunction(varray,fi)
c           use parameters
c           use variable_setup
c           use grid
c           real(8) :: fi(n_elem)
c           type (var_array),pointer :: varray
c         end subroutine evaluateNonlinearFunction
c      END INTERFACE

c Begin program

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions
cc      xx  => gv%gparams%xx
cc      yy  => gv%gparams%yy
cc      zz  => gv%gparams%zz
cc      dxh => gv%gparams%dxh
cc      dyh => gv%gparams%dyh
cc      dzh => gv%gparams%dzh
cc      dx  => gv%gparams%dx
cc      dy  => gv%gparams%dy
cc      dz  => gv%gparams%dz
cc
cc      jcnv => gv%aux%vec_list(IJCNV)%vec
cc      jcov => gv%aux%vec_list(IJCOV)%vec
cc#if defined(vec_pot)
cc      bcnv => gv%aux%vec_list(IBCNV)%vec
cc      bcov => gv%aux%vec_list(IBCOV)%vec
cc#endif
cc
cc      imin = patch_var%gparams%ilog
cc      imax = patch_var%gparams%ihig
cc      jmin = patch_var%gparams%jlog
cc      jmax = patch_var%gparams%jhig
cc      kmin = patch_var%gparams%klog
cc      kmax = patch_var%gparams%khig
cc      varray => patch_var%u_n
cc
ccc Set pointers for data
cc      pcnv => gv%aux_spec%pcnv
cc      vecnv => gv%aux_spec%vecnv
cc      vdummy => gv%aux_spec%vdummy
cc      vcnv => gv%aux_spec%vcnv
cc      vcov => gv%aux_spec%vcov
cc      grd => gv%aux_spec%grd
cc#ifdef defined(vec_pot)
cc      acnv => gv%aux_spec%acnv
cc      acov => gv%aux_spec%acov
cc      vlap_a => gv%aux_spec%vlap_a
cc#else 
cc      bcnv => gv%aux_spec%bcnv
cc      bcov => gv%aux_spec%bcov
cc#endif
cc      b_n => gv%aux_spec%b_n
cc      dv_dt => gv%aux_spec%dv_dt
cc      p_n => gv%aux_spec%p_n
cc      eeta => gv%aux_spec%eeta
cc      nuu => gv%aux_spec%nuu
cc      divrgV => gv%aux_spec%divrgV

c Find local limits

      call fromGlobalToLocalLimits(imingc ,jmingc ,kmingc
     $                            ,imingcl,jmingcl,kmingcl,1,1,1)
      call fromGlobalToLocalLimits(imaxgc ,jmaxgc ,kmaxgc
     $                            ,imaxgcl,jmaxgcl,kmaxgcl,1,1,1)

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
              f%array_var(ieq)%array(i,j,k)=fsrc(ii+ieq)-ff(ii+ieq)
            enddo
          enddo
        enddo
      enddo


c Deallocate structures

c      call deallocateDerivedType(varray)

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
      type (var_array),pointer :: varray

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
           real(8) :: ff(neqd)
           integer :: i,j,k,igx,igy,igz
           type(var_array),pointer :: varray
         end subroutine nonlinearRHS
      END INTERFACE

c Begin program

c Setup parallel BC flags to indicate PETSc provides BCs

#if defined(petsc_c)
      call setup_petsc_BC
#endif

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

      end subroutine

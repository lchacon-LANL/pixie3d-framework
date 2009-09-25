c formEquilibrium
c######################################################################
      subroutine formEquilibrium(patch_var)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use variable_setup

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use generalPurposeFunctions

      use auxiliaryVariables

      implicit none

c Call variables

c      integer         :: imin,imax,jmin,jmax,kmin,kmax,ieq

      type(patch), TARGET :: patch_var

c Local variables

c Begin program

cc      allocate(patch_var%aux_spec)

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

c Check for autoinitializing parameters

      pi = acos(-1d0)

      if (maxitnwt.eq.0) 
     .      maxitnwt = max(floor(1.5*log(rtol)/log(tolgm)),10)

      alpha = 1. - cnfactor

      dtbase = dt   

c Allocate global variable

      call allocateGlobalVar(gv)

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd,g_pack,gv%gparams)

c Initialize vector dimensions for global and local problems

cc      ihig = imax
cc      ilog = imin
cc      jhig = jmax
cc      jlog = jmin
cc      khig = kmax
cc      klog = kmin

      call setVectorDimensions

cc      write (*,*) 'Processor:',my_rank
cc     .           ,'; Grid limits:',ilog,ihig,jlog,jhig,klog,khig

c Allocate constant arrays

cc      allocate(gv%gparams%zeros (ilom:ihip,jlom:jhip,klom:khip))
cc      allocate(gv%gparams%vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
cc      allocate(gv%gparams%ones  (ilom:ihip,jlom:jhip,klom:khip))
cc
cc      zeros  => gv%gparams%zeros
cc      vzeros => gv%gparams%vzeros
cc      ones   => gv%gparams%ones

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Create nonlinear solver

      call createNonlinearSolver

c Create nonlinear function

      call createNonlinearFunction

c Create equilibrium u_0

      call createEquilibrium

c Transfer to Petsc format
c
c      do ieq=1,neqd
c        array(ilog:ihig,jlog:jhig,klog:khig)%var(ieq)
c     .      = u_0%array_var(ieq)%array(ilo:ihi,jlo:jhi,klo:khi)
c      enddo

c Get pointers to internal data

cc      gv%aux_spec%pcnv => pcnv
cc      gv%aux_spec%vecnv => vecnv
cc      gv%aux_spec%vdummy => vdummy
cc      gv%aux_spec%vcnv => vcnv
cc      gv%aux_spec%vcov => vcov
cc      gv%aux_spec%grd => grd
cc#ifdef defined(vec_pot)
cc      gv%aux_spec%acnv => acnv
cc      gv%aux_spec%acov => acov
cc      gv%aux_spec%vlap_a => vlap_a
cc#else 
cc      gv%aux_spec%bcnv =>bcnv 
cc      gv%aux_spec%bcov => bcov
cc#endif
cc      gv%aux_spec%b_n => b_n
cc      gv%aux_spec%dv_dt => dv_dt           
cc      gv%aux_spec%p_n => p_n
cc      gv%aux_spec%eeta => eeta
cc      gv%aux_spec%nuu => nuu
cc      gv%aux_spec%divrgV => divrgV
      
c Nullify pointers to internal data

cc      nullify(pcnv)
cc      nullify(vecnv)
cc      nullify(vdummy)
cc      nullify(vcnv)
cc      nullify(vcov)
cc      nullify(grd)
cc#ifdef defined(vec_pot)
cc      nullify(acnv)
cc      nullify(acov)
cc      nullify(vlap_a)
cc#else 
cc      nullify(bcnv)
cc      nullify(bcov)
cc#endif
cc      nullify(b_n)
cc      nullify(dv_dt)          
cc      nullify(p_n)
cc      nullify(eeta)
cc      nullify(nuu)
cc      nullify(divrgV)

c End program
 
      end subroutine formEquilibrium

c createNonlinearSolver
c######################################################################
      subroutine createNonlinearSolver

c----------------------------------------------------------------------
c     Allocates nonlinear solver arrays
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variables

      implicit none

c Call variables

c Local variables

c Begin program

      call allocateStructures

cc#if !defined(samrai)
cc      allocate(cnf(neqd),one_over_dt(neqd)
cc     .        ,bdfp(neqd),bdfn(neqd),bdfnm(neqd))
cc      allocate(fold(ntotd))
cc      allocate(fsrc(ntotd))
cc#endif

c End programs

      end subroutine createNonlinearSolver

c createEquilibrium
c######################################################################
      subroutine createEquilibrium

c----------------------------------------------------------------------
c     Creates equilibrium according to user input.
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variable_setup

      use grid

      implicit none

c Call variables

c Local variables

      integer    :: ieq,igx,igy,igz

      real(8),allocatable,dimension(:,:,:,:) :: var

      character*(20),allocatable,dimension(:):: label

      integer   ,allocatable,dimension(:,:)  :: bcs

c Begin program

      igx = 1
      igy = 1
      igz = 1

c Set equilibrium u_0 and define BCs on all variables

      allocate(var(ilom:ihip,jlom:jhip,klom:khip,neqd)
     $            ,label(neqd),bcs(6,neqd))

      !Initialize boundary conditions
      do ieq = 1,neqd
        bcs(:,ieq) = bcond
      enddo

      var = 0d0

      call setEquilibrium(igx,igy,igz,var,bcs,label)

      call packVariables(u_0)

      deallocate(var,label,bcs)

c End programs

      contains

c     packVariables
c     #################################################################
      subroutine packVariables(varray)

c     -----------------------------------------------------------------
c     Packs arrays into variable structure 
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      integer    :: ieq

c     Begin program

      varray%nvar = neqd

      do ieq = 1,neqd
        call varPack (var(:,:,:,ieq),bcs(:,ieq),label(ieq)
     .               ,ieq,varray)
      enddo

c     End program

      end subroutine PackVariables

      end subroutine createEquilibrium


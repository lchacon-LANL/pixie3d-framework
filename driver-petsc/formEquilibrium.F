c formEquilibrium
c######################################################################
      subroutine formEquilibrium(array,imin,imax,jmin,jmax,kmin,kmax)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use generalPurposeFunctions

      implicit none

c Call variables

      integer         :: imin,imax,jmin,jmax,kmin,kmax,ieq

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

c Begin program

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

      call setVectorDimensions

      if (.not.tst_flg)
     .   write (*,'(a,i4,a,6i5)') 'Processor:',my_rank
     .           ,'; Grid limits:',ilog,ihig,jlog,jhig,klog,khig

c Allocate constant arrays

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

      do ieq=1,neqd
        array(ilog:ihig,jlog:jhig,klog:khig)%var(ieq)
     .      = u_0%array_var(ieq)%array(ilo:ihi,jlo:jhi,klo:khi)
      enddo

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
      allocate(cnf(neqd),one_over_dt(neqd)
     .        ,bdfp(neqd),bdfn(neqd),bdfnm(neqd))
      allocate(fold(ntotd),fsrc(ntotd))

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


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

      use debug

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax,ieq

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

c Debug

      integer(4) :: nt,ii,i,j,k,igx,nx,ny,nz
      real(8) :: mag
      real(8),allocatable,dimension(:,:) :: v_mat,debug2

c Begin program

c Check for autoinitializing parameters

      pi = acos(-1d0)

      if (maxitnwt.eq.0) 
     .      maxitnwt = max(floor(1.5*log(rtol)/log(tolgm)),10)

      alpha = 1. - cnfactor

      dtbase = dt   

      sm_flag= 0
      if (cnfactor.lt.0d0) then
        sm_flag= 1
      else
        sm_pass= 0
      endif

      cnf_d = 1d0

c Initialize vector dimensions for global and local problems

      ihig = imax
      ilog = imin
      jhig = jmax
      jlog = jmin
      khig = kmax
      klog = kmin

      call setVectorDimensions

c Allocate constant arrays

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Initialize MG, MPI, fortran DA, and create grid

      call createGrid(ilog,ihig,jlog,jhig,klog,khig,nxd,nyd,nzd)

c Check grid

c diag ******* NUMERICAL GRID QUANTITIES
cc      if (my_rank == 2) then
cc        write (*,*) 'proc',my_rank,'lim',ilog,ihig,jlog,jhig,klog,khig
cc        write (*,*) 'Grid in processor',my_rank
cc        call checkGrid
cc      endif
cc
cc      if (my_rank == 0) then
cc        open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc        igx = 1
cc        nx = grid_params%nxv(igx)
cc        ny = grid_params%nyv(igx)
cc        nz = grid_params%nzv(igx)
cc        allocate(debug2(0:nx+1,0:ny+1))
cccc        debug2 = gmetric%grid(igx)%jac(0:nx+1,0:ny+1,1)
cccc        call contour(debug2(0:nx+1,0:ny+1),nx+2,ny+2,0d0
cccc     .            ,xmax,0d0,ymax,0,110)
cc        do i=1,3
cc          do j=1,3
cccc            debug2 = gmetric%grid(igx)%cnv(0:nx+1,0:ny+1,1,i,j)
cccc     .              *gmetric%grid(igx)%jac(0:nx+1,0:ny+1,1)
cc          debug2 = gmetric%grid(igx)%Gamma(0:nx+1,0:ny+1,1,1,i,j)
cc     .            *gmetric%grid(igx)%jac(0:nx+1,0:ny+1,1)
cc            call contour(debug2(0:nx+1,0:ny+1),nx+2,ny+2,0d0
cccc     .            ,xmax,0d0,ymax,1,110)
cc     .           ,xmax,0d0,ymax,i+j-2,110)
cc          enddo
cc        enddo
cc        deallocate(debug2)
cc        close(110)
cc      endif
cc
cc      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc      call MPI_Finalize(mpierr)
cc      stop
c diag *******

c Set output file

      urecord    = 25 + my_rank

      if (np > 1) then
        recordfile = 'record_proc'//trim(int2char(my_rank))//'.bin'
      else
        recordfile = 'record.bin'
      endif

      hdf5file='record.h5'

c Create nonlinear solver

      call createNonlinearSolver

c Create nonlinear function

      call createNonlinearFunction

c Create equilibrium u_0

      call createEquilibrium

c Initialize old time solution

      u_n = u_0   !Overloaded assignment

c Transfer to Petsc format

      do ieq=1,neqd
        array(ilog:ihig,jlog:jhig,klog:khig)%var(ieq)
     .      = u_n%array_var(ieq)%array(ilo:ihi,jlo:jhi,klo:khi)
      enddo

      call deallocateDerivedType(u_n)       !DO NOT deallocate u_0; needed in imposeBC routine

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

      integer(4) :: ieq,igx,igy,igz

      real(8),allocatable,dimension(:,:,:,:) :: var

      character*(20),allocatable,dimension(:):: label

      integer(4),allocatable,dimension(:,:)  :: bcs

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

      type (var_array) :: varray

c     Local variables

      integer(4) :: ieq

c     Begin program

      varray%nvar = neqd

      do ieq = 1,neqd
        call varPack (var(:,:,:,ieq),bcs(:,ieq),label(ieq)
     .               ,ieq,varray)
      enddo

c     End program

      end subroutine PackVariables

      end subroutine createEquilibrium


c initializeCalculation
c######################################################################
      subroutine initializeCalculation

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use mg_setup

      use variables

      use resistiveWall

      use timeStepping

      use newtongm

      use precond_setup

      use iosetup

      use constants

      use discret_params

      use graphics

      use icond

      use grid

      implicit none

c Call variables

c Local variables

      double precision,allocatable,dimension(:,:,:) :: var

      character*(5),allocatable,dimension(:) :: label

      integer(4),allocatable,dimension(:,:)  :: bcs

cc      integer :: igrid

c Begin program

      pi = acos(-1d0)

c Initializations (application specific)

      call readInput

c Initialize time stepping

      if (maxitnwt.eq.0) 
     .      maxitnwt = max(floor(1.5*log(tolnewt)/log(tolgm)),10)

      alpha = 1. - cnfactor

      dtbase = dt   

      sm_flag= 0
      if (cnfactor.lt.0d0) then
        sm_flag= 1
      else
        sm_pass= 0
      endif

      cnf_d = 1d0

cc      if (cnfactor.eq.1d0) then
cc        conserv = .false.
cc        advect  = 3
cc      endif

c Initialize vector dimensions

      call setVectorDimensions

c Allocate variables

      call allocatePhysicsVariables

      allocate(fold(ntotd),fsrc(ntotd))

      call allocateSpecificVariables

      allocate(zeros(0:nxdp,0:nydp))

      zeros = 0d0

c Initialize MG and create grid

      call setupMG(nxd,nyd,bcond)

cc      do igrid=ngrd,1,-1
cc        write (*,*) xl(:,igrid)
cc      enddo
cc      do igrid=ngrd,1,-1
cc        write (*,*) yl(:,igrid)
cc      enddo

      xmin = 0d0
      xmax = xlength
      ymin = 0d0
      ymax = ylength
      zmin = 0d0
      zmax = 1d0
      call createGrid(nxd,nyd,1,bcond)
cc      call checkGrid

c Set equilibrium u_0 and define BCs on all variables

      allocate(var(0:nxdp,0:nydp,neqd),label(neqd),bcs(4,neqd))

      call setEquilibrium(var,bcs,label)

      call packVariables(var,bcs,label,u_0)

      deallocate(var,label,bcs)

c Set unperturbed forcing fields

      if (source) then
        call evaluateNonlinearFunction(u_0,fsrc)
      else
        fsrc = 0d0
      endif

c Set initial condition u_np

      call setInitialCondition(u_0,u_np)

c Calculate derived quantities for plots

      call setupNonlinearFunction(u_np)

c Initialize variable storage template

      call equateDerivedType(u_0,utmp)

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        write(*,*)
        write(*,*) 'Tmax is less or equal than restarting time'
        write(*,*) 'Aborting'
        stop
      endif

      allocate(cnf(neqd),one_over_dt(neqd))

      call defineTSParameters

c End program

      end subroutine initializeCalculation

c packVariables
c######################################################################
      subroutine packVariables(var,bcs,label,varray)

c----------------------------------------------------------------------
c     Packs arrays into variable structure 
c----------------------------------------------------------------------

      use parameters

      use variable_setup
      
      implicit none

c     Call variables

      type (var_array) :: varray

      double precision :: var(0:nxdp,0:nydp,neqd)

      character*(5) :: label(neqd)

      integer(4)    :: bcs(4,neqd)

c     Local variables

      integer(4) :: ieq

c     Begin program

      varray%nvar = neqd

      do ieq = 1,neqd
        call VarPack (nxd,nyd,var(:,:,ieq),bcs(:,ieq),label(ieq)
     .               ,ieq,varray)
      enddo

c     End program

      end subroutine PackVariables

c setInitialCondition
c####################################################################
      subroutine setInitialCondition(varray0,varray)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation. Variables
c     are stored in u_np
c--------------------------------------------------------------------

      use icond

      use mg_setup

      use variable_setup

      use timeStepping

      implicit none

c Call variables

      type (var_array) :: varray0,varray

c Local variables

      integer(4) :: i,j,bcs(4),ieq,nx,ny

c Begin program

c Initialize new time step solution

      call equateDerivedType(varray0,varray)

c Perturb equilibrium

      if (.not.restart) then

        do ieq = 1,neqd
          call setperturb(varray %array_var(ieq)%array
     .                   ,varray0%array_var(ieq)%array,pert(ieq))
        enddo

        time     = 0d0
        inewtime = 1

      else

        call readRestartFile (inewtime,nx,ny,varray)

        if (nx.ne.nxd.or.ny.ne.nyd) then
           write (*,*) 'Grid meshes do not agree; cannot restart'
           stop
        endif

      endif

c End program

      end subroutine

c setperturb
c######################################################################
      subroutine setperturb(array,array0,perturb)

c----------------------------------------------------------------------
c     Introduces perturbations in initial equilibrium
c----------------------------------------------------------------------

      use mg_setup

      use parameters

      use constants

      use resistiveWall

      use icond

      implicit none

c Call variables

      double precision :: perturb
      double precision :: array(0:nxdp,0:nydp),array0(0:nxdp,0:nydp)

c Local variables

      integer*4   i,j
      real*8      rand1,rand2,rand3,rand4,rand5,rand6,rand7,rand8
      real*8      xx(0:nxdp),yy(0:nydp),fy

c Externals

      real*8      laplace
      external    laplace

c Begin program

      xx(:) = xl(:,ngrd)
      yy(:) = yl(:,ngrd)

      if (nh.eq.0) then

        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        call random_number(rand4)
        call random_number(rand5)
        call random_number(rand6)
        call random_number(rand7)
        call random_number(rand8)

        do j = 0,nyd+1
          do i = 0,nxd+1
            array(i,j) = array0(i,j)
     .          + perturb*dsin(pi*yy(j))
     .          *( dsin(2 *pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand1*2*pi)
     .            +dsin(4 *pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand2*2*pi)
     .            +dsin(6 *pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand3*2*pi)
     .            +dsin(8 *pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand4*2*pi)
     .            +dsin(10*pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand5*2*pi)
     .            +dsin(12*pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand6*2*pi) 
     .            +dsin(14*pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand7*2*pi)
     .            +dsin(16*pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-rand8*2*pi))

          enddo
        enddo

      else

        do j = 0,nyd+1
          do i = 0,nxd+1
            if (reswall .and. init == 'slow') then
              fy = yy(j)
            else
cc              fy = dsin(2*pi*yy(j))
              fy = dsin(pi*yy(j))
            endif
            array(i,j) = array0(i,j)
     .            + perturb*fy
     .            * sin(nh*2*pi*(xx(i)-xx(1))/(xx(nxd)-xx(1))-pi/2.)
          enddo
        enddo

      endif

c End program

      return
      end

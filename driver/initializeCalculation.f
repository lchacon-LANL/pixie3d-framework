c initializeCalculation
c######################################################################
      subroutine initializeCalculation

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use precond_setup

      use iosetup

      use constants

      use graphics

      use icond

      implicit none

c Call variables

c Local variables

c Begin program

c Read user initializations

      call readInput

c Check for autoinitializing parameters

      pi = acos(-1d0)

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

c Initialize vector dimensions

      call setVectorDimensions

c Initialize MG and create grid

      call createGrid(nxd,nyd,nzd)
cc      call checkGrid

c Create nonlinear solver

      call createNonlinearSolver

c Create nonlinear function

      call createNonlinearFunction

c Create equilibrium u_0

      call createEquilibrium

c Initialize old time solution

      u_n = u_0   !Overloaded assignment

c Set unperturbed forcing fields

      if (source) then
        call evaluateNonlinearFunction(u_n,fsrc)
      else
        fsrc = 0d0
      endif

c Perturb initial condition u_n

      call setInitialCondition(u_n)

c Initialize new time solution

      u_np = u_n

c Create graphics

      call createGraphics

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        write(*,*)
        write(*,*) 'Tmax is less or equal than restarting time'
        write(*,*) 'Aborting'
        stop
      endif

c End program

      end subroutine initializeCalculation

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

      allocate(fold(ntotd),fsrc(ntotd))

      allocate(zeros(0:nxdp,0:nydp,0:nzdp))

      zeros = 0d0

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

      integer(4) :: ieq

      real(8),allocatable,dimension(:,:,:,:) :: var

      character*(20),allocatable,dimension(:):: label

      integer(4),allocatable,dimension(:,:)  :: bcs

c Begin program

c Set equilibrium u_0 and define BCs on all variables

      allocate(var(0:nxdp,0:nydp,0:nzdp,neqd),label(neqd),bcs(6,neqd))

      !Initialize boundary conditions
      do ieq = 1,neqd
        bcs(:,ieq) = bcond
      enddo

      var = 0d0

      call setEquilibrium(var,bcs,label)

      call packVariables(u_0)

cc      call imposeBoundaryConditions(u_0)

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
        call VarPack (nxd,nyd,nzd,var(:,:,:,ieq),bcs(:,ieq),label(ieq)
     .               ,ieq,varray)
      enddo

c     End program

      end subroutine PackVariables

      end subroutine createEquilibrium

c setInitialCondition
c####################################################################
      subroutine setInitialCondition(varray)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation. Variables
c     are stored in u_np
c--------------------------------------------------------------------

      use icond

      use grid

      use variable_setup

      use timeStepping

      use constants

      implicit none

c Call variables

      type (var_array) :: varray

c Local variables

      integer(4) ::  ieq,nx,ny,nz

c Begin program

c Perturb equilibrium

      if (.not.restart) then

        do ieq = 1,neqd
          call perturbEquilibrium(varray %array_var(ieq)%array
     .                           ,varray %array_var(ieq)%bconds
     .                           ,pert(ieq))
        enddo

cc        call imposeBoundaryConditions(varray)

        time     = 0d0
        inewtime = 1

      else

        call readRestartFile (inewtime,nx,ny,nz,varray)

        if (nx.ne.nxd.or.ny.ne.nyd.or.nz.ne.nzd) then
           write (*,*) 'Grid meshes do not agree; cannot restart'
           write (*,*) 'Aborting.'
           stop
        endif

      endif

c End program

      contains

c     perturbEquilibrium
c     #################################################################
      subroutine perturbEquilibrium(array,bcs,perturb)

c     -----------------------------------------------------------------
c     Perturbs equilibrium quantity in array0 with a sinusoidal
c     perturbation of magnitude perturb, and introduces it in array.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcs(6)
      real(8)    :: perturb
      real(8)    :: array (0:nxdp,0:nydp,0:nzdp)

c     Local variables

      integer*4 :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)   :: x1,y1,z1,car(3)
      real(8)   :: fx(0:nxdp),fy(0:nydp),fz(0:nzdp) 

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      do i = 1,nxd
        call getCurvilinearCoordinates(i,1,1,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fx(i) = factor(xmin,xmax,x1,bcs(1:2),nh1)
      enddo

      do j = 1,nyd
        call getCurvilinearCoordinates(1,j,1,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fy(j) = factor(ymin,ymax,y1,bcs(3:4),nh2)
      enddo

      do k = 1,nzd
        call getCurvilinearCoordinates(1,1,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fz(k) = factor(zmin,zmax,z1,bcs(5:6),nh3)
      enddo

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd
            array(i,j,k) = array(i,j,k) + perturb*fx(i)*fy(j)*fz(k)
          enddo
        enddo
      enddo

c     End program

      end subroutine perturbEquilibrium

c     factor
c     ####################################################################
      real(8) function factor(xmin,xmax,x,bcs,nh) result(ff)

        implicit none

        real(8)    :: xmin,xmax,x,period
        integer(4) :: bcs(2),nh

        period = pi
        if (odd) period = 2*pi

        if (bcs(1) == PER) then
          ff = cos(nh*2*pi*(x-xmin)/(xmax-xmin))
        elseif (random) then
          call random_number(ff)
        elseif (bcs(1) == NEU .and. bcs(2) == NEU) then
          ff = cos(period*(x-xmin)/(xmax-xmin))
        elseif (bcs(1) == NEU .and. bcs(2) == DIR) then
          period = 3*period/4.
          if (.not.odd) period = period/2.
          ff = cos(period*(x-xmin)/(xmax-xmin))
        elseif (bcs(1) == DIR .and. bcs(2) == NEU) then
          period = 3*period/4.
          if (.not.odd) period = period/2.
          ff = sin(period*(x-xmin)/(xmax-xmin))
        elseif (bcs(1) == SP .and. bcs(2) == DIR) then
          ff = x**(nh+1)*sin(period*(x-xmin)/(xmax-xmin)) !To satisfy regularity at r=0 (r^m)
        elseif (bcs(1) == SP .and. bcs(2) == NEU) then
          period = 3*period/4.
          if (.not.odd) period = period/2.
          ff = x**(nh+1)*cos(period*(x-xmin)/(xmax-xmin)) !To satisfy regularity at r=0 (r^m)
        else
          ff = sin(period*(x-xmin)/(xmax-xmin))
        endif

      end function factor

      end subroutine setInitialCondition

c createGraphics
c######################################################################
      subroutine createGraphics

c----------------------------------------------------------------------
c     Creates graphics pointers, defines dumping intervals
c----------------------------------------------------------------------

      use parameters

      use constants

      use timeStepping

      use variables

      use graphics

      implicit none

c Call variables

c Local variables

      integer(4) :: i

c Begin program

c Set data dumping intervals

      dfreq = 8d0

      if (tmax.gt.0d0) then
        if (dstep.eq.0d0) dstep = dt*max(int((tmax-time)/dfreq/dt),1)
        rstep = min(dt*max(int((tmax-time)/dfreq/dt),1),dstep)
        ndstep = -1
      else
        if (ndstep.eq.0) ndstep = max(numtime/int(dfreq),1)
        nrstep = min(max(numtime/int(dfreq),1),ndstep)
        dstep = 1e30
      endif

c Initialize graphics

      u_graph = u_np - u_0
cc      call imposeBoundaryConditions(u_graph)

      if (plot) then
        call initializeDiagnostics
        call initializeGraphics(1,1,1,bcond,restart)
      endif

c End programs

      end subroutine createGraphics


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

      use constants

      use iosetup

      implicit none

c Call variables

c Local variables

c Begin program

      urecord    = 25
      recordfile = 'record.bin'

c Read user initializations

      call readInput

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

c Initialize vector dimensions

      ihig = nxd
      ilog = 1
      jhig = nyd
      jlog = 1
      khig = nzd
      klog = 1

      call setVectorDimensions

c Allocate constant arrays

      allocate(zeros (0:nxdp,0:nydp,0:nzdp))
      allocate(vzeros(0:nxdp,0:nydp,0:nzdp,3))
      allocate(ones  (0:nxdp,0:nydp,0:nzdp))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Initialize MG and create grid

      call createGrid(1,nxd,1,nyd,1,nzd,nxd,nyd,nzd)
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

      !This not only evaluates fsrc, but defines BCs on u_n
      call evaluateNonlinearFunction(u_n,fsrc)

      if (.not.source) fsrc = 0d0

c Set initial condition

      call setInitialCondition(u_n,u_np)

c Initialize record file

      call initializeRecordFile

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

cc      call imposeBoundaryConditions(u_0,igx,igy,igz)

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
        call VarPack (var(:,:,:,ieq),bcs(:,ieq),label(ieq),ieq,varray)
      enddo

c     End program

      end subroutine PackVariables

      end subroutine createEquilibrium

c setInitialCondition
c####################################################################
      subroutine setInitialCondition(varrayn,varraynp)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation. Variables
c     are stored in u_np
c--------------------------------------------------------------------

      use icond

      use grid

      use variable_setup

      use timeStepping

      use constants

      use iosetup

      implicit none

c Call variables

      type (var_array) :: varrayn,varraynp

c Local variables

      integer(4) ::  ieq,nx,ny,nz

c Begin program

c Perturb equilibrium

      if (.not.restart) then

        do ieq = 1,neqd
          call perturbEquilibrium(varrayn%array_var(ieq)%array
     .                           ,abs(varrayn%array_var(ieq)%bconds)
     .                           ,pert(ieq))
        enddo

        time     = 0d0
        inewtime = 1

        varraynp = varrayn

      else

        call readRestartFile (nx,ny,nz,itime,time,varrayn,varraynp)

        inewtime = itime

        if (nx.ne.nxd.or.ny.ne.nyd.or.nz.ne.nzd) then
           write (*,*) 'Grid meshes do not agree; cannot restart'
           write (*,*) 'Aborting.'
           stop
        endif

      endif

c End program

      contains

c     readRestartFile
c     #################################################################
      subroutine readRestartFile(nx,ny,nz,itime,time,vn,vnp)

c     -----------------------------------------------------------------
c     Reads restart file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer(4),intent(OUT) :: nx,ny,nz,itime
      real(8),intent(OUT)    :: time

      type (var_array)       :: vn,vnp

c     Local variables

      integer(4) :: ierr

      type (var_array)       :: vmed

c     Begin program

      call allocateDerivedType(vmed)

      open(unit=urecord,file=recordfile,form='unformatted',status='old')

      read (urecord) nx
      read (urecord) ny
      read (urecord) nz

      write (*,*) ' Reading restart file...'

      call readRecord(urecord,itime,time,dt,vn,ierr)

      vnp = vn

      do
        call readRecord(urecord,itime,time,dt,vmed,ierr)

        if (ierr /= 0) then
          exit
        else
          vn  = vnp
          vnp = vmed
        endif
      enddo

      close (urecord)

      write (*,*) 'Done!'

c     End

      call deallocateDerivedType(vmed)

      end subroutine readRestartFile

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
      real(8)   :: x1,y1,z1
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
c For MSW test case
cc            call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,y1,z1)
cc            jac = gmetric%grid(igx)%jac(i,j,k)
cc            if (ieq == 1) jac = 1d0
cc            array(i,j,k) = array(i,j,k)
cc     .                   + jac*perturb*cos(2*pi*(x1-xmin)/(xmax-xmin)
cc     .                                    +2*pi*(y1-ymin)/(ymax-ymin))
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

c initializeRecordFile
c######################################################################
      subroutine initializeRecordFile

c----------------------------------------------------------------------
c     Creates graphics pointers, defines dumping intervals
c----------------------------------------------------------------------

      use timeStepping

      use variables

      use iosetup

      implicit none

c Call variables

c Local variables

c Begin program

c Set data dumping intervals

      dfreq = 8d0

      if (tmax.gt.0d0) then
        if (dstep.eq.0d0) then
          dstep = dt*max(int((tmax-time)/dfreq/dt),1)
        else
          dstep = max(dstep,dt)
        endif
        rstep = min(dt*max(int((tmax-time)/dfreq/dt),1),dstep)
        ndstep = -1
        numtime = -1
      else
        if (ndstep.eq.0) ndstep = max(numtime/int(dfreq),1)
        nrstep = min(max(numtime/int(dfreq),1),ndstep)
        dstep = 1e30
        tmax  = 0d0
      endif

c Open record file

      if (.not.restart) then

        !Impose BC's on u_graph <- u_0 (do not overwrite u_0, since it contains equil. BCs)
        u_graph = u_0
        !initially dump u_n instead of u_0 (w/BCs) for comparison w/ preconditioner solution
cc        u_graph = u_n
        !Check source
cc        u_graph = fsrc
cc        u_graph%array_var(1)%array = u_0%array_var(1)%array  !Set rho finite
cc        call imposeBoundaryConditions(u_graph,1,1,1)

        !Open record file
        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='replace')

        write (urecord) nxd,1,nxd
        write (urecord) nyd,1,nyd
        write (urecord) nzd,1,nzd

        call writeRecordFile(urecord,0,0d0,dt,u_graph)
        call writeRecordFile(urecord,0,0d0,dt,u_np)

      else

        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='old',position='append')
          
      endif

c End programs

      end subroutine initializeRecordFile


c Subroutines below original "subroutine initializeCalculation()"

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

      integer(4) :: ieq

      real(8),allocatable,dimension(:,:,:,:) :: var

      character*(20),allocatable,dimension(:):: label

      integer(4),allocatable,dimension(:,:)  :: bcs

c Begin program

c Set equilibrium u_0 and define BCs on all variables

      allocate(var(ilom:ihip,jlom:jhip,klom:khip,neqd)
     $            ,label(neqd),bcs(6,neqd))

      !Initialize boundary conditions
      do ieq = 1,neqd
        bcs(:,ieq) = bcond
      enddo

      var = 0d0

      call setEquilibrium(var,bcs,label)

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

        call readRestartFile (itime,time,varrayn,varraynp)

        inewtime = itime + 1

      endif

c End program

      contains

c     readRestartFile
c     #################################################################
      subroutine readRestartFile(itime,time,vn,vnp)

c     -----------------------------------------------------------------
c     Reads restart file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer(4),intent(OUT) :: itime
      real(8),intent(OUT)    :: time

      type (var_array)       :: vn,vnp

c     Local variables

      integer(4)             :: ierr,nx,ny,nz,il,jl,kl,ih,jh,kh

      type (var_array)       :: vmed

c     Begin program

      call allocateDerivedType(vmed)

      open(unit=urecord,file=recordfile,form='unformatted',status='old')

      read (urecord) nx
      read (urecord) ny
      read (urecord) nz

      if (nx /= nxl .or. ny /= nyl .or. nz /= nzl) then
        write (*,*) 'Grid meshes do not agree; cannot restart'
        write (*,*) 'Aborting.'
        stop
      endif

      write (*,*) ' Reading restart file...'

      call readRecord(itime,time,dt,vn,ierr)

      vnp = vn

      do
        call readRecord(itime,time,dt,vmed,ierr)

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
      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)

c     Local variables

      integer*4 :: i,j,k,ig,jg,kg,igx,igy,igz
      real(8)   :: x1,y1,z1,car(3),jac
      real(8)   :: fx(ilo:ihi),fy(jlo:jhi),fz(klo:khi) 

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      do i = ilo,ihi
        call getCurvilinearCoordinates(i,jlo,klo,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fx(i) = factor(xmin,xmax,x1,bcs(1:2),nh1)
      enddo

      do j = jlo,jhi
        call getCurvilinearCoordinates(ilo,j,klo,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fy(j) = factor(ymin,ymax,y1,bcs(3:4),nh2)
      enddo

      do k = klo,khi
        call getCurvilinearCoordinates(ilo,jlo,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fz(k) = factor(zmin,zmax,z1,bcs(5:6),nh3)
      enddo

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
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

      integer(4) :: i

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

        call imposeBoundaryConditions(u_graph,1,1,1)

        !Open record file
        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='replace')

        write (urecord) nxl
        write (urecord) nyl
        write (urecord) nzl

        call writeRecordFile(0,0d0,dt,u_graph)
        call writeRecordFile(0,0d0,dt,u_np)

      else

        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='old',position='append')
          
      endif

c End programs

      end subroutine initializeRecordFile


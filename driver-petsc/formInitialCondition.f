c formInitialCondition
c######################################################################
      subroutine formInitialCondition(array,imin,imax,jmin,jmax
     .                               ,kmin,kmax,time_to_c)

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

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax
      real(8)         :: time_to_c

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer(4)      :: iminl,imaxl,jminl,jmaxl,kminl,kmaxl,ieq

c Begin program

      call fromGlobalToLocalLimits(imin,jmin,kmin,iminl,jminl,kminl
     $                            ,1,1,1)
      call fromGlobalToLocalLimits(imax,jmax,kmax,imaxl,jmaxl,kmaxl
     $                            ,1,1,1)

c Map petsc array (to get BCs from Petsc)

      do ieq=1,neqd
        u_0%array_var(ieq)
     .       %array(iminl:imaxl,jminl:jmaxl,kminl:kmaxl)
     .      = array(imin :imax ,jmin :jmax ,kmin :kmax )%var(ieq)
      enddo

c Initialize old time solution

      call initializeDerivedType(u_n)

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
        if (my_rank == 0) then
          write(*,*)
          write(*,*) 'Tmax is less or equal than restarting time'
          write(*,*) 'Aborting'
        endif
        call PetscFinalize(mpierr)
        stop
      endif

c Initialize counters

      itime     = inewtime-1

      time_to_c = time

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0

c Transfer to Petsc format

      do ieq=1,neqd
        array(imin :imax ,jmin :jmax ,kmin :kmax )%var(ieq)
     .      = u_np%array_var(ieq)
     .            %array(iminl:imaxl,jminl:jmaxl,kminl:kmaxl)
      enddo

      call deallocateDerivedType(u_np)

c End program

      end subroutine formInitialCondition

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

      read (urecord) nx,ilog,ihig
      read (urecord) ny,jlog,jhig
      read (urecord) nz,klog,khig

      if ((nx /= nxl .or. ny /= nyl .or. nz /= nzl)
     .    .and.(my_rank == 0)) then
        write (*,*) 'Grid meshes do not agree; cannot restart'
        write (*,*) 'Aborting...'
        stop
      endif

      if (my_rank == 0) write (*,*) 'Reading restart file...'

      call readRecordFile(urecord,itime,time,dt,vn,ierr)

      vnp = vn

      do
        call readRecordFile(urecord,itime,time,dt,vmed,ierr)

        if (ierr /= 0) then
          exit
        else
          vn  = vnp
          vnp = vmed
        endif
      enddo

      close (urecord)

      if (my_rank == 0) write (*,*) 'Done!'

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
      real(8)   :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 

c     Begin program

      igx = 1
      igy = 1
      igz = 1

      do i = ilom,ihip
        call getCurvilinearCoordinates(i,jlo,klo,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fx(i) = factor(xmin,xmax,x1,bcs(1:2),nh1)
      enddo

      do j = jlom,jhip
        call getCurvilinearCoordinates(ilo,j,klo,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fy(j) = factor(ymin,ymax,y1,bcs(3:4),nh2)
      enddo

      do k = klom,khip
        call getCurvilinearCoordinates(ilo,jlo,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
        fz(k) = factor(zmin,zmax,z1,bcs(5:6),nh3)
      enddo

      do k = klom,khip
        do j = jlom,jhip
          do i = ilom,ihip
            array(i,j,k) = array(i,j,k) + perturb*fx(i)*fy(j)*fz(k)
          enddo
        enddo
      enddo

c     End program

      end subroutine perturbEquilibrium

c     factor
c     ####################################################################
      function factor(xmin,xmax,x,bcs,nh) result(ff)

        implicit none

        real(8)    :: xmin,xmax,x,period,ff
        integer(4) :: bcs(2),nh
        logical    :: neumann(2),dirichlet(2)

        neumann   = (bcs == NEU) .or. (bcs == SYM)
        dirichlet = (bcs == DIR) .or. (bcs ==-SYM) .or. (bcs == EQU)

        period = pi
        if (odd) period = 2*pi

        if (bcs(1) == PER) then
          ff = cos(nh*2*pi*(x-xmin)/(xmax-xmin))
        elseif (random) then
          call random_number(ff)
        elseif (neumann(1) .and. neumann(2)) then
          ff = cos(period*(x-xmin)/(xmax-xmin))
        elseif (neumann(1) .and. dirichlet(2)) then
          if (.not.odd) then
            period = period/2.
          else
            period = 3*period/4.
          endif
          ff = cos(period*(x-xmin)/(xmax-xmin))
        elseif (dirichlet(1) .and. neumann(2)) then
          if (.not.odd) then
            period = period/2.
          else
            period = 3*period/4.
          endif
          ff = sin(period*(x-xmin)/(xmax-xmin))
        elseif (bcs(1) == SP .and. dirichlet(2)) then
          ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
     .         *sign(1d0,sin(period*(x-xmin)/(xmax-xmin)))
        elseif (bcs(1) == SP .and. neumann(2)) then
          if (.not.odd) then
            period = period/2.
            ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
          else
            period = 3*period/4.
            ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
     .        *sign(1d0,sin(period*(x-xmin)/(xmax-xmin)))
          endif
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

      integer(4) :: i,ierr

c External

      integer(4) :: system
      external   :: system

c Begin program

c Set data dumping intervals

      dfreq = 8d0

      if (tmax.gt.0d0) then
        if (dstep.eq.0d0) then
          dstep = dt*max(int((tmax-time)/dfreq/dt),1)
cc        else
cc          dstep = max(dstep,dt)
        endif
        rstep = min(dt*max(int((tmax-time)/dfreq/dt),1),dstep)
        ndstep  = -1
        numtime = -1
      else
        if (ndstep.eq.0) ndstep = max(numtime/int(dfreq),1)
        nrstep = min(max(numtime/int(dfreq),1),ndstep)
        dstep = 1e30
        tmax  = 0d0
      endif

c Open record file

      if (.not.restart) then

        if (my_rank == 0) ierr=system('rm -f record*.bin >& /dev/null')

        call MPI_Barrier(MPI_COMM_WORLD,mpierr)

        u_graph = u_0    !Needed for BCs
cc        u_graph = u_n  !For debugging to compare agains PC solution
cc        u_graph = fsrc !For debugging source

        !Open record file
        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='unknown')

        write (urecord) nxl,ilog,ihig
        write (urecord) nyl,jlog,jhig
        write (urecord) nzl,klog,khig

cc        call writeRecordFile(urecord,0,0d0,dt,u_n)
        call writeRecordFile(urecord,0,0d0,dt,u_graph)
        call writeRecordFile(urecord,0,0d0,dt,u_np)

      else

        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='old',position='append')
          
      endif

c End program

      end subroutine initializeRecordFile


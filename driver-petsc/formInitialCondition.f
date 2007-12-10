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

      use generalpurposefunctions

      implicit none

c Call variables

      integer         :: imin,imax,jmin,jmax,kmin,kmax
      real(8)         :: time_to_c

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer         :: iminl,imaxl,jminl,jmaxl,kminl,kmaxl,ieq
     .                  ,ierr

c External

      integer    :: system
      external      system

c Interfaces

      INTERFACE
         subroutine setInitialCondition(varrayn,varraynp)
         use variable_setup
         type(var_array),pointer :: varrayn,varraynp
         end subroutine setInitialCondition
      END INTERFACE

      INTERFACE
        subroutine evaluateNonlinearFunction(varray,fi)
        use parameters
        use variable_setup
        real(8)          :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluateNonlinearFunction
      END INTERFACE

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

      call equateDerivedType(u_n,u_0)

c Set unperturbed forcing fields

      !This not only evaluates fsrc, but defines BCs on u_n
      call evaluateNonlinearFunction(u_n,fsrc)

      if (.not.source) fsrc = 0d0

c Set output file

      if (.not.restart) then
        if (my_rank == 0)
     .     ierr=system('rm -f '//trim(recordfile)//'* >& /dev/null')
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      endif

      urecord = urecord + my_rank

      if (np > 1) then
        recordfile=trim(recordfile)//'_proc'//trim(int2char(my_rank))
      endif

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

      type(var_array),pointer :: varrayn,varraynp

c Local variables

      integer    ::  ieq,nx,ny,nz

c Interfaces

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
        use variable_setup
        integer    :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

c Begin program

c Perturb equilibrium

      if (.not.restart) then

        !Equate to varraynp
        call equateDerivedType(varraynp,varrayn)

        do ieq = 1,neqd
          call perturbEquilibrium(varraynp%array_var(ieq)%array
     .                           ,varraynp%array_var(ieq)%bconds
     .                           ,pert(ieq),ieq)
        enddo

        time     = 0d0
        inewtime = 1

cc        call equateDerivedType(varraynp,varrayn)

        !Impose BC on varrayn
        call imposeBoundaryConditions(varraynp,1,1,1)

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

      integer   ,intent(OUT) :: itime
      real(8),intent(OUT)    :: time

      type(var_array),pointer       :: vn,vnp

c     Local variables

      integer                :: ierr,nx,ny,nz,il,jl,kl,ih,jh,kh

      type(var_array),pointer       :: vmed

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

      if (my_rank == 0) write (*,*) 'Reading restart file(s)...'

      call readRecordFile(urecord,itime,time,dt,vn,ierr)

      call equateDerivedType(vnp,vn)

      do
        call readRecordFile(urecord,itime,time,dt,vmed,ierr)

        if (ierr /= 0) then
          exit
        else
          call equateDerivedType(vn ,vnp )
          call equateDerivedType(vnp,vmed)
        endif
      enddo

      close (urecord)

      if (my_rank == 0) write (*,*) 'Done!'

c     End

      call deallocateDerivedType(vmed)

      end subroutine readRestartFile

ccc     perturbEquilibrium
ccc     #################################################################
cc      subroutine perturbEquilibrium(array,bcs,perturb,ieq)
cc
ccc     -----------------------------------------------------------------
ccc     Perturbs equilibrium quantity in array0 with a sinusoidal
ccc     perturbation of magnitude perturb, and introduces it in array.
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: bcs(6),ieq
cc      real(8)    :: perturb
cc      real(8)    :: array (ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Local variables
cc
cc      integer    :: i,j,k,ig,jg,kg,igx,igy,igz
cc      real(8)    :: x1,y1,z1,car(3),jac
cc      real(8)    :: fx(ilom:ihip),fy(jlom:jhip),fz(klom:khip) 
cc
ccc     Begin program
cc
cc      igx = 1
cc      igy = 1
cc      igz = 1
cc
cc      do i = ilom,ihip
cc        call getCurvilinearCoordinates(i,jlo,klo,igx,igy,igz,ig,jg,kg
cc     .                                ,x1,y1,z1)
cc        fx(i) = factor(xmin,xmax,x1,bcs(1:2),nh1)
cc      enddo
cc
cc      do j = jlom,jhip
cc        call getCurvilinearCoordinates(ilo,j,klo,igx,igy,igz,ig,jg,kg
cc     .                                ,x1,y1,z1)
cc        fy(j) = factor(ymin,ymax,y1,bcs(3:4),nh2)
cc      enddo
cc
cc      do k = klom,khip
cc        call getCurvilinearCoordinates(ilo,jlo,k,igx,igy,igz,ig,jg,kg
cc     .                                ,x1,y1,z1)
cc        fz(k) = factor(zmin,zmax,z1,bcs(5:6),nh3)
cc      enddo
cc
cc      do k = klom,khip
cc        do j = jlom,jhip
cc          do i = ilom,ihip
cc            array(i,j,k) = array(i,j,k) + perturb*fx(i)*fy(j)*fz(k)
cc          enddo
cc        enddo
cc      enddo
cc
ccc     End program
cc
cc      end subroutine perturbEquilibrium
cc
ccc     factor
ccc     ####################################################################
cc      function factor(xmin,xmax,x,bcs,nh) result(ff)
cc
cc        implicit none
cc
cc        real(8)    :: xmin,xmax,x,period,ff
cc        integer    :: bcs(2),nh
cc
cc        logical    :: neumann(2),dirichlet(2),spoint(2)
cc
cc        spoint    = (bcs == SP )
cc        neumann   = (bcs == NEU) .or. (bcs == SYM)
cc        dirichlet = (bcs == DIR) .or. (bcs ==-SYM) .or. (bcs == EQU)
cc
cc        period = pi
cc        if (odd) period = 2*pi
cc
cc        if (bcs(1) == PER) then
cc          ff = cos(nh*2*pi*(x-xmin)/(xmax-xmin))
cc        elseif (random) then
cc          call random_number(ff)
cc        elseif (neumann(1) .and. neumann(2)) then
cc          ff = cos(period*(x-xmin)/(xmax-xmin))
cc        elseif (neumann(1) .and. dirichlet(2)) then
cc          if (.not.odd) then
cc            period = period/2.
cc          else
cc            period = 3*period/4.
cc          endif
cc          ff = cos(period*(x-xmin)/(xmax-xmin))
cc        elseif (dirichlet(1) .and. neumann(2)) then
cc          if (.not.odd) then
cc            period = period/2.
cc          else
cc            period = 3*period/4.
cc          endif
cc          ff = sin(period*(x-xmin)/(xmax-xmin))
cc        elseif (spoint(1) .and. dirichlet(2)) then
cc          ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
cc     .         *sign(1d0,sin(period*(x-xmin)/(xmax-xmin)))**(nh+1)
cc        elseif (spoint(1) .and. neumann(2)) then
cc          if (.not.odd) then
cc            period = period/2.
cc            ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
cc          else
cc            period = 3*period/4.
cc            ff = (sin(period*(x-xmin)/(xmax-xmin)))**(nh+2) !To satisfy regularity at r=0 (r^m)
cc     .        *sign(1d0,sin(period*(x-xmin)/(xmax-xmin)))
cc          endif
cc        else
cc          ff = sin(period*(x-xmin)/(xmax-xmin))
cc        endif
cc
cc      end function factor

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

      integer    :: i,ierr

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

cc        call equateDerivedType(u_graph,u_0)
cc        u_graph = u_0    !Needed for BCs
cccc        u_graph = u_n  !For debugging to compare agains PC solution
cccc        u_graph = fsrc !For debugging source

        !Open record file
        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='unknown')

        write (urecord) nxl,ilog,ihig
        write (urecord) nyl,jlog,jhig
        write (urecord) nzl,klog,khig

        call writeRecordFile(urecord,0,0d0,dt,u_n)
cc        call writeRecordFile(urecord,0,0d0,dt,u_graph)
        call writeRecordFile(urecord,0,0d0,dt,u_np)

      else

        open(unit=urecord,file=recordfile
     .      ,form='unformatted',status='old',position='append')
          
      endif

c End program

      end subroutine initializeRecordFile


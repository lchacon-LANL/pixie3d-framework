c formInitialCondition
c######################################################################
      subroutine formInitialCondition(patch_var,n_elem,fsrc)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

c      use constants

      use iosetup

c      use generalpurposefunctions

c      use auxiliaryVariables

c      use grid_aliases

      implicit none

c Call variables

      integer :: n_elem
      real(8) :: time_to_c
      type(patch), TARGET :: patch_var
      real(8),dimension(n_elem) :: fsrc

c Local variables

      integer :: iminl,imaxl,jminl,jmaxl,kminl,kmaxl,ieq,ierr

c Interfaces

      INTERFACE
         subroutine setInitialCondition(varrayn,varraynp)
         use variable_setup
         type(var_array),pointer :: varrayn,varraynp
         end subroutine setInitialCondition
      END INTERFACE

      INTERFACE
c        Initialize the fluxes
         subroutine evaluateNonlinearResidualFlux(patch_var)
           use parameters
           use grid
           use variables
           implicit none
           type(patch), TARGET :: patch_var
         end subroutine evaluateNonlinearResidualFlux
      END INTERFACE

      INTERFACE
c        evaluate the nonlinearRHS
         subroutine evaluateNonlinearFunction(varray,ff)
           use parameters
           use variable_setup
           use grid
           implicit none
           type (var_array),pointer :: varray
           real(8)          :: ff(ntotd)
         end subroutine evaluateNonlinearFunction
      END INTERFACE

c Initialize gv

      gv => patch_var
      grid_params => patch_var%gparams
      gmetric => patch_var%gparams%gmetric
      u_n => patch_var%u_n
      u_0 => patch_var%u_0

      call setVectorDimensions

      time_to_c = 0.0

c Set unperturbed forcing fields

      !This not only evaluates fsrc, but defines BCs on u_n
      call evaluateNonlinearResidualFlux(patch_var)
      call evaluateNonlinearFunction(u_n,fsrc)

      if (.not.source) fsrc = 0d0

ccc Set output file
cc
cc      if (.not.restart) then
cc        if (my_rank == 0)
cc     .     ierr=shellcmd('rm -f '//trim(recordfile)//'* >& /dev/null')
ccc        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc      endif
cc
cc      urecord = urecord + my_rank
cc
cc      if (np > 1) then
cc        recordfile=trim(recordfile)//'_proc'//trim(int2char(my_rank))
cc      endif

c Set initial condition

      call setInitialCondition(u_n,u_np)

cc Initialize record file
c
c      call initializeRecordFile

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        if (my_rank == 0) then
          write(*,*)
          write(*,*) 'Tmax is less or equal than restarting time'
          write(*,*) 'Aborting'
        endif
c        call PetscFinalize(mpierr)
        stop
      endif

c Initialize counters

      itime     = inewtime-1

      time_to_c = time

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0
 
c Transfer data

c      do ieq=1,neqd
c        patch_var%u_n%array_var(ieq)%array
c     .      = patch_var%u_np%array_var(ieq)%array
c      enddo

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

      use grid

      use variable_setup

      use timeStepping

c      use constants

      use iosetup

      implicit none

c Call variables

      type(var_array),pointer :: varrayn,varraynp

c Local variables

      integer    ::  ieq,nx,ny,nz

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use variable_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

c Begin program

c Perturb equilibrium

      if (.not.restart) then

        call perturbEquilibrium(varrayn)

        time     = 0d0
        inewtime = 1

        call equateDerivedType(varraynp,varrayn)

      else

        call pstop('setInitialCondition'
     .            ,'No restart capability implemented in SAMRAI')

cc        call readRestartFile (itime,time,varrayn,varraynp)
cc
cc        inewtime = itime + 1

      endif

c End program

      contains

ccc     readRestartFile
ccc     #################################################################
cc      subroutine readRestartFile(itime,time,vn,vnp)
cc
ccc     -----------------------------------------------------------------
ccc     Reads restart file
ccc     -----------------------------------------------------------------
cc
cc      use variables
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer   ,intent(OUT) :: itime
cc      real(8),intent(OUT)    :: time
cc
cc      type(var_array),pointer       :: vn,vnp
cc
ccc     Local variables
cc
cc      integer                :: ierr,nx,ny,nz,il,jl,kl,ih,jh,kh
cc
cc      type(var_array),pointer       :: vmed
cc
ccc     Begin program
cc
cc      call allocateDerivedType(vmed)
cc
cc      open(unit=urecord,file=recordfile,form='unformatted',status='old')
cc
cc      read (urecord) nx,ilog,ihig
cc      read (urecord) ny,jlog,jhig
cc      read (urecord) nz,klog,khig
cc
cc      if ((nx /= nxl .or. ny /= nyl .or. nz /= nzl)
cc     .    .and.(my_rank == 0)) then
cc        write (*,*) 'Grid meshes do not agree; cannot restart'
cc        write (*,*) 'Aborting...'
cc        stop
cc      endif
cc
cc      if (my_rank == 0) write (*,*) 'Reading restart file(s)...'
cc
cc      ierr = readRecordFile(urecord,itime,time,dt,vn)
cc
cc      vnp = vn
cc
cc      do
cc        ierr = readRecordFile(urecord,itime,time,dt,vmed)
cc
cc        if (ierr /= 0) then
cc          exit
cc        else
cc          vn  = vnp
cc          vnp = vmed
cc        endif
cc      enddo
cc
cc      close (urecord)
cc
cc      if (my_rank == 0) write (*,*) 'Done!'
cc
ccc     End
cc
cc      call deallocateDerivedType(vmed)
cc
cc      end subroutine readRestartFile

      end subroutine setInitialCondition

ccc initializeRecordFile
ccc######################################################################
cc      subroutine initializeRecordFile
cc
ccc----------------------------------------------------------------------
ccc     Creates graphics pointers, defines dumping intervals
ccc----------------------------------------------------------------------
cc
cc      use timeStepping
cc
cc      use variables
cc
cc      use iosetup
cc
cc      implicit none
cc
ccc Call variables
cc
ccc Local variables
cc
cc      integer    :: i,ierr
cc
ccc Begin program
cc
ccc Set data dumping intervals
cc
cc      dfreq = 8d0
cc
cc      if (tmax.gt.0d0) then
cc        if (dstep.eq.0d0) then
cc          dstep = dt*max(int((tmax-time)/dfreq/dt),1)
cccc        else
cccc          dstep = max(dstep,dt)
cc        endif
cc        rstep = min(dt*max(int((tmax-time)/dfreq/dt),1),dstep)
cc        ndstep  = -1
cc        numtime = -1
cc      else
cc        if (ndstep.eq.0) ndstep = max(numtime/int(dfreq),1)
cc        nrstep = min(max(numtime/int(dfreq),1),ndstep)
cc        dstep = 1e30
cc        tmax  = 0d0
cc      endif
cc
ccc Open record file
cc
cc      if (.not.restart) then
cc
cc        call equateDerivedType(u_graph,u_0)
cccc        u_graph = u_0    !Needed for BCs
cccc        u_graph = u_n  !For debugging to compare agains PC solution
cccc        u_graph = fsrc !For debugging source
cc
cc        !Open record file
cc        open(unit=urecord,file=recordfile
cc     .      ,form='unformatted',status='unknown')
cc
cc        write (urecord) nxl,ilog,ihig
cc        write (urecord) nyl,jlog,jhig
cc        write (urecord) nzl,klog,khig
cc
cccc        call writeRecordFile(urecord,0,0d0,dt,u_n)
cc        call writeRecordFile(urecord,0,0d0,dt,u_graph)
cc        call writeRecordFile(urecord,0,0d0,dt,u_np)
cc
cc      else
cc
cc        open(unit=urecord,file=recordfile
cc     .      ,form='unformatted',status='old',position='append')
cc          
cc      endif
cc
ccc End program
cc
cc      end subroutine initializeRecordFile
cc

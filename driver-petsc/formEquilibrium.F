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

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      type(petsc_array) :: petscarray

c Begin program

c Initialize global domain limit variables

      ihig = imax
      ilog = imin
      jhig = jmax
      jlog = jmin
      khig = kmax
      klog = kmin

c Initialize parallel MPI info

      call initMPI

c Set output file

      urecord    = 25 + my_rank

      if (np > 1) then
        recordfile = 'record_proc'//trim(int2char(my_rank))//'.bin'
      else
        recordfile = 'record.bin'
      endif

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

      call setVectorDimensions

c Allocate constant arrays

      allocate(zeros (ilom:ihip,jlom:jhip,klom:khip))
      allocate(vzeros(ilom:ihip,jlom:jhip,klom:khip,3))
      allocate(ones  (ilom:ihip,jlom:jhip,klom:khip))

      zeros  = 0d0
      vzeros = 0d0
      ones   = 1d0

c Initialize MG and create grid

      call createGrid(ilog,ihig,jlog,jhig,klog,khig,nxd,nyd,nzd)

cc      if (my_rank == 3) then
cc        write (*,*) 'Grid in processor',my_rank
cc        call checkGrid
cc      endif
cc      stop

c Create nonlinear solver

      call createNonlinearSolver

c Create nonlinear function

      call createNonlinearFunction

c Create equilibrium u_0

      call createEquilibrium

c Initialize old time solution

      u_n = u_0   !Overloaded assignment

c Transfer to Petsc format

      petscarray = u_n

      array = petscarray%array(ilo:ihi,jlo:jhi,klo:khi)

      call deallocatePetscType(petscarray)
      call deallocateDerivedType(u_n)

c End program

      contains

c     initMPI
c     #####################################################################
      subroutine initMPI

c     ---------------------------------------------------------------------
c     Initialize needed MPI constructs
c     ---------------------------------------------------------------------

c     Call variables

c     Local variables

        integer(4),allocatable,dimension(:) :: process_ranks
        integer(4)        :: lnp_sp,gnp_sp,inp,rilog

c     Rank and number of processors

        call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
        call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)

cc      write (*,*) 'proc',my_rank,'coords',imin,imax,jmin,jmax,kmin,kmax
cc
cc      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc      call MPI_Finalize(mpierr)
cc      stop

c     i=1 boundary communicator (for singular point calculations)

        if (bcond(1) /= SP) return

        !Number of processes (-> gnp_sp, in ALL processors)
        lnp_sp = 0
        if (ilog == 1) lnp_sp = 1

        call MPI_Allreduce(lnp_sp,gnp_sp,1,MPI_INTEGER
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)

        allocate(process_ranks(0:gnp_sp-1))

        !Select processes
        if (my_rank == 0) then
          lnp_sp = 0
          if (ilog == 1) then
            process_ranks(lnp_sp)=0
            lnp_sp = lnp_sp+1
          endif
          do inp=1,np-1
            tag = 0
            call MPI_Recv(rilog,1,MPI_INTEGER,inp,tag
     .                   ,MPI_COMM_WORLD,status,mpierr)
            if (rilog == 1) then 
              process_ranks(lnp_sp) = inp
              lnp_sp = lnp_sp+1
            endif
          enddo
        else
          tag  = 0
          dest = 0
          call MPI_Send(ilog,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD
     .                 ,mpierr)
        endif

        !Synchronize process_ranks in ALL processes
        call MPI_Bcast(process_ranks,gnp_sp,MPI_INTEGER,root
     .                ,MPI_COMM_WORLD,mpierr)

        !Get group underlying MPI_COMM_WORLD
        call MPI_Comm_group(MPI_COMM_WORLD,group_world,mpierr)

        !Create new group using process_ranks (-> group_sp)
        call MPI_Group_incl(group_world,gnp_sp,process_ranks,group_sp
     .                     ,mpierr)

        !Create new communicator (-> MPI_COMM_SP)
        call MPI_Comm_create(MPI_COMM_WORLD,group_sp,MPI_COMM_SP,mpierr)

cc      write (*,*) 'sp procs',process_ranks
cc
cc      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc      call MPI_Finalize(mpierr)
cc      stop

        deallocate(process_ranks)

      end subroutine initMPI

      end subroutine formEquilibrium


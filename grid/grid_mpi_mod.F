c module grid_mpi
c ######################################################################
      module grid_mpi

        use io, ONLY: FATAL,sstop

        use grid_operations

        use xdraw_io, ONLY: contour

#if defined(use_pic_mpi)
        use mpi
#endif

        implicit none

c$$$        INTERFACE dot
c$$$          module procedure dot1,dot2,dot3
c$$$        END INTERFACE

        INTERFACE find_global
          module procedure find_global_field_vec,find_global_field_scl
        END INTERFACE

        INTERFACE find_global_nobc
          module procedure find_global_field_vec_nobc
     .                    ,find_global_field_scl_nobc
        END INTERFACE

        INTERFACE find_local
          module procedure find_local_field_vec,find_local_field_scl
        END INTERFACE

        integer :: mpierr

#if defined(petsc)

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdadef.h"

        integer    :: np=1
     .               ,my_rank=0
     .               ,group_world
     .               ,group_sp,tag=0,dest=0,root=0
cc     .               ,status(MPI_STATUS_SIZE)
     .               ,tag_send,tag_recv,request

        integer    :: l_lim(3,0:1)
     .               ,g_lim(3,0:1)
     .               ,rem_l_lim(3,0:1)
     .               ,rp_l_lim(3,0:1)

        character(80) :: messg

        type(grid_mg_def),pointer,private :: g_mpi_def

      contains

ccc     initFortranMPI
ccc     #####################################################################
cc      subroutine initFortranMPI(mpi_comm,nxg,nyg,nzg,mx,my,mz)
cc
ccc     ---------------------------------------------------------------------
ccc     Initializes MPI variables.
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer,INTENT(IN) :: nxg,nyg,nzg,mpi_comm
cc        integer,INTENT(OUT):: mx,my,mz
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        call initMPI(mpi_comm,np,my_rank)
cc
ccc     Allocate processors
cc
cc        call processorAlloc(nxg,nyg,nzg,mx,my,mz)
cc
cc      end subroutine initFortranMPI

c     create_MPI_comms
c     #################################################################
      subroutine create_MPI_comms(g_def)

c     -----------------------------------------------------------------
c     Creates required MPI communicators
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer,INTENT(IN)  :: g_def

c     Local variables

        integer :: my_sp_rank,np_sp

        logical :: chk_mpi_comm=.false.

c     Begin program

        g_mpi_def => g_def

        if (g_mpi_def%mpi_comm == MPI_COMM_SELF) return

c     Singular-point boundary MPI communicator

        if (bcSP()) then
          call splitMPIComm(g_mpi_def%MPI_COMM_SP,selectRank_SP
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
          if (chk_mpi_comm)
     .         call check_MPI_COMM(g_mpi_def%MPI_COMM_SP,'MPI_COMM_SP')
        endif

c     X MPI communicator (domains that share the same Y, Z  slices)

        call splitMPIComm(g_mpi_def%MPI_COMM_X,selectRank_X
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)

        if (bcSP()) g_mpi_def%MPI_COMM_RAD = g_mpi_def%MPI_COMM_X
        if (chk_mpi_comm)
     .       call check_MPI_COMM(g_mpi_def%MPI_COMM_X,'MPI_COMM_X')

c     Y MPI communicator (domains that share same X, Z slices)

        call splitMPIComm(g_mpi_def%MPI_COMM_Y,selectRank_Y
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)

        if (bcSP()) g_mpi_def%MPI_COMM_POL = g_mpi_def%MPI_COMM_Y
        if (chk_mpi_comm)
     .       call check_MPI_COMM(g_mpi_def%MPI_COMM_Y,'MPI_COMM_Y')

c     Z MPI communicator (domains that share same X, Z slices)

        call splitMPIComm(g_mpi_def%MPI_COMM_Z,selectRank_Z
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
        if (chk_mpi_comm)
     .       call check_MPI_COMM(g_mpi_def%MPI_COMM_Z,'MPI_COMM_Z')

c     XY MPI communicator (domains that share same Z slices)

        call splitMPIComm(g_mpi_def%MPI_COMM_XY,selectRank_XY
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
        if (chk_mpi_comm)
     .       call check_MPI_COMM(g_mpi_def%MPI_COMM_XY,'MPI_COMM_XY')

c     XZ MPI communicator (domains that share same Y slices)

        call splitMPIComm(g_mpi_def%MPI_COMM_XZ,selectRank_XZ
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
        if (chk_mpi_comm)
     .       call check_MPI_COMM(g_mpi_def%MPI_COMM_XZ,'MPI_COMM_XZ')

c     YZ MPI communicator (domains that share same X slices)

        call splitMPIComm(g_mpi_def%MPI_COMM_YZ,selectRank_YZ
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
        if (chk_mpi_comm)
     .       call check_MPI_COMM(g_mpi_def%MPI_COMM_YZ,'MPI_COMM_YZ')

c     Periodic BCs communicators

        if (bcond(1) == PER) then
          call splitMPIComm(g_mpi_def%MPI_COMM_PER(1),selectRank_PER_X
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
          if (chk_mpi_comm)
     .    call check_MPI_COMM(g_mpi_def%MPI_COMM_PER(1),'MPI_COMM_PERx')
        endif

        if (bcond(3) == PER) then
          call splitMPIComm(g_mpi_def%MPI_COMM_PER(2),selectRank_PER_Y
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
          if (chk_mpi_comm)
     .    call check_MPI_COMM(g_mpi_def%MPI_COMM_PER(2),'MPI_COMM_PERy')
        endif

        if (bcond(5) == PER) then
          call splitMPIComm(g_mpi_def%MPI_COMM_PER(3),selectRank_PER_Z
     .                     ,MPI_COMM_ORIG=g_mpi_def%mpi_comm)
          if (chk_mpi_comm)
     .    call check_MPI_COMM(g_mpi_def%MPI_COMM_PER(3),'MPI_COMM_PERz')
        endif

        if (chk_mpi_comm) write(*,*) 'Finished checking MPI comms'

c     End program

      end subroutine create_MPI_comms

c     check_MPI_COMM
c     ###################################################################
      subroutine check_MPI_COMM(mpi_comm,mpi_comm_id)

c     -------------------------------------------------------------------
c     Checks MPI communicator
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: mpi_comm
      character(*) :: mpi_comm_id

c     Local variables

      integer :: np_l,my_rank_l

c     Begin program

      if (mpi_comm /= MPI_COMM_NULL) then
        call MPI_Comm_size(MPI_COMM,np_l     ,mpierr)
        call MPI_Comm_rank(MPI_COMM,my_rank_l,mpierr)
        write (*,'(a,i4,a,i4,a,i4,a,i4)')
     .       trim(mpi_comm_id)//'=',MPI_COMM
     .                ,'; nproc_local=',np_l
     .                ,'; my_rank_local=',my_rank_l
     .                ,'; my_rank_global=',my_rank
      endif

c$$$      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
c$$$      call MPI_Finalize(mpierr)
c$$$      stop

c     End program

      end subroutine check_MPI_COMM

c     fillPetscGhostCells
c     #################################################################
      subroutine fillPetscGhostCells(dactx,array)

c     -----------------------------------------------------------------
c     Performs communication to fill ghost cells in array using PETSc
c     distributed arrays in DA context dactx.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(petsc_da_ctx),pointer :: dactx

        real(8) :: array(dactx%ilm:dactx%ihp
     .                  ,dactx%jlm:dactx%jhp 
     .                  ,dactx%klm:dactx%khp)

c     Local variables

        Vec :: localX

        PetscScalar,pointer :: lx_v(:)

        integer    :: ierr

c     Begin program

c     Fill known values of dactx%Xg

        !Get pointer to vector data
        call DAGetLocalVector(dactx%da,localX,ierr)
cc        call fchkerrq(ierr,'DAGetLocalVector')

        call VecGetArrayF90  (localX,lx_v,ierr)

        !Assignment: array -> lx_v
        call fillLocalVec(dactx,array,lx_v)

        !Restore vector
        call VecRestoreArrayF90(localX,lx_v,ierr)

        !Insert values into global vector
        call DALocalToGlobal(dactx%da,localX
     .                      ,INSERT_VALUES,dactx%Xg
     $                      ,ierr)
cc        call DARestoreLocalVector(dactx%da,localX,ierr)

c     Fill ghost cells

        !Get pointer to vector data w/ghost cells
cc        call DAGetLocalVector    (dactx%da,localX,ierr)
        call DAGlobalToLocalBegin(dactx%da,dactx%Xg
     .                           ,INSERT_VALUES,localX
     $                           ,ierr)

        call DAGlobalToLocalEnd  (dactx%da,dactx%Xg
     .                           ,INSERT_VALUES,localX
     $                           ,ierr)


        call VecGetArrayF90(localX,lx_v,ierr)

        !Assignment: lx_v -> array
        call emptyLocalVec(dactx,array,lx_v)

        !Restore vector
        call VecRestoreArrayF90(localX,lx_v,ierr)

c     Deallocate memory

        call DARestoreLocalVector(dactx%da,localX,ierr)

c     End program

      end subroutine fillPetscGhostCells

c     initMPI
c     #####################################################################
      subroutine initMPI(mpi_comm,np,my_rank)

c     ---------------------------------------------------------------------
c     Initializes MPI variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: mpi_comm,np,my_rank

c     Local variables

c     Begin program

c     Rank and number of processors

        call MPI_Comm_rank(mpi_comm,my_rank,mpierr)
        call MPI_Comm_size(mpi_comm,np     ,mpierr)

      end subroutine initMPI

c     splitMPIComm
c     #####################################################################
      subroutine splitMPIComm(MPI_COMM,selProc,MPI_COMM_ORIG)

c     ---------------------------------------------------------------------
c     Creates MPI communicator MPI_COMM according to routine selProc.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: MPI_COMM
        integer, optional :: MPI_COMM_ORIG

        external   selProc

c     Local variables

        integer :: split_key,my_rank_l,mpi_orig

c     Begin program

        if (PRESENT(MPI_COMM_ORIG)) then
          mpi_orig = MPI_COMM_ORIG
        else
          mpi_orig = MPI_COMM_WORLD
        endif

c     Find split key

        split_key = MPI_UNDEFINED

        call selProc(split_key)

c     Create communicator

        call MPI_Comm_rank(mpi_orig,my_rank_l,mpierr)

        call MPI_Comm_split(mpi_orig,split_key,my_rank_l
     .                     ,MPI_COMM,mpierr)

c     End program

      end subroutine splitMPIComm

c     selectRank_SP
c     #####################################################################
      subroutine selectRank_SP(split_key)

c     ---------------------------------------------------------------------
c     Groups processors surrounding singular point
c     (SP boundary condition)
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: split_key

c     Local variables

        integer :: dim,loc,BC,npz

        logical :: cp(0:1)

c     Begin program

        BC = SP
        dim = 1

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npz = g_mpi_def%mz

        split_key = MPI_UNDEFINED

        do loc=0,1
          cp(loc) = checkProc(BC,dim,loc,l_lim,g_lim)
          if (cp(loc)) split_key = BC
        enddo

        !Group communicators based on location in z direction
        if (bcond(dim) == BC .and. cp(0)) then
          split_key = split_key + npz*l_lim(3,1)/g_lim(3,1)
        endif

cc        !Create self-communicator for single-point dimensions
cc        if (      (l_lim(dim,0) == l_lim(dim,1))
cc     .      .and. (cp(0)       .eqv. cp(1)     )  ) then
cc          split_key=split_key+my_rank+1
cc        endif

c     End program

      end subroutine selectRank_SP

c     selectRank_PER_X
c     #####################################################################
      subroutine selectRank_PER_X(split_key)

c     ---------------------------------------------------------------------
c     Groups processors for periodic BC transfer along dimension dim.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(OUT) :: split_key

c     Local variables

        integer :: dim,loc,BC

        logical :: cp(0:1)

c     Begin program

        BC = PER
        dim = 1

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Check whether proc is at boundary of appropriate type

        split_key = MPI_UNDEFINED

        do loc=0,1
          cp(loc) = checkProc(BC,dim,loc,l_lim,g_lim)
          if (cp(loc)) split_key = BC
        enddo

c     Group communicators based on location in face

        if (split_key /= MPI_UNDEFINED) then
          call selectRank_X(split_key)
        endif

cc        !Create self-communicator for single-point dimensions
cc        if (      (l_lim(dim,0) == l_lim(dim,1))
cc     .      .and. (cp(0)       .eqv. cp(1)     )  ) then
cc          split_key=split_key+my_rank+1
cc        endif

c     End program

      end subroutine selectRank_PER_X

c     selectRank_PER_Y
c     #####################################################################
      subroutine selectRank_PER_Y(split_key)

c     ---------------------------------------------------------------------
c     Groups processors for periodic BC transfer along dimension dim.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(OUT) :: split_key

c     Local variables

        integer :: dim,loc,BC

        logical :: cp(0:1)

c     Begin program

        BC = PER
        dim = 2

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Check whether proc is at boundary of appropriate type

        split_key = MPI_UNDEFINED

        do loc=0,1
          cp(loc) = checkProc(BC,dim,loc,l_lim,g_lim)
          if (cp(loc)) split_key =BC
        enddo

c     Group communicators based on location in face

        if (split_key /= MPI_UNDEFINED) then
          call selectRank_Y(split_key)
        endif

cc        !Create self-communicator for single-point dimensions
cc        if (      (l_lim(dim,0) == l_lim(dim,1))
cc     .      .and. (cp(0)       .eqv. cp(1)     )  ) then
cc          split_key=split_key+my_rank+1
cc        endif

c     End program

      end subroutine selectRank_PER_Y

c     selectRank_PER_Z
c     #####################################################################
      subroutine selectRank_PER_Z(split_key)

c     ---------------------------------------------------------------------
c     Groups processors for periodic BC transfer along dimension dim.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(OUT) :: split_key

c     Local variables

        integer :: dim,loc,BC

        logical :: cp(0:1)

c     Begin program

        BC = PER
        dim = 3

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        split_key = MPI_UNDEFINED

        !Check whether proc is at boundary of appropriate type
        do loc=0,1
          cp(loc) = checkProc(BC,dim,loc,l_lim,g_lim)
          if (cp(loc)) split_key = BC
        enddo

        !Group communicators based on location in face
        if (split_key /= MPI_UNDEFINED) then
          call selectRank_Z(split_key)
        endif

        !Create self-communicator for single-point dimensions
cc        if (      (l_lim(dim,0) == l_lim(dim,1))
cc     .      .and. (cp(0)       .eqv. cp(1)     )  ) then
cc          split_key=split_key+my_rank+1
cc        endif

c     End program

      end subroutine selectRank_PER_Z

c     selectRank_X
c     #####################################################################
      subroutine selectRank_X(split_key)

c     ---------------------------------------------------------------------
c     Groups processors in slices in X direction (i.e., order them according
c     to Y and Z positions).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(INOUT) :: split_key

c     Local variables

        integer :: npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npx = g_mpi_def%mx
        npy = g_mpi_def%my
        npz = g_mpi_def%mz

        !Group processors based on location in y and z directions
        !(x direction is shared; 
        ! organize in 2D processor mesh: k + nz*(j-1))
        split_key = split_key +     npz*l_lim(3,1)/g_lim(3,1)
     .                        + npz*npy*l_lim(2,1)/g_lim(2,1)

c     End program

      end subroutine selectRank_X

c     selectRank_Y
c     #####################################################################
      subroutine selectRank_Y(split_key)

c     ---------------------------------------------------------------------
c     Groups processors in slices in Y direction (i.e., order them according
c     to X and Z positions).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(INOUT) :: split_key

c     Local variables

        integer :: npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npx = g_mpi_def%mx
        npy = g_mpi_def%my
        npz = g_mpi_def%mz
        
        !Group processors based on location in x and z directions
        !(y direction is shared; 
        ! organize in 2D processor mesh: k + nz*(i-1))
        split_key = split_key +     npz*l_lim(3,1)/g_lim(3,1)
     .                        + npz*npx*l_lim(1,1)/g_lim(1,1)

c     End program

      end subroutine selectRank_Y

c     selectRank_Z
c     #####################################################################
      subroutine selectRank_Z(split_key)

c     ---------------------------------------------------------------------
c     Groups processors in slices in Z direction (i.e., order them according
c     to X and Y positions).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(INOUT) :: split_key

c     Local variables

        integer :: npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npx = g_mpi_def%mx
        npy = g_mpi_def%my
        npz = g_mpi_def%mz

        !Group processors based on location in x and y directions
        !(z direction is shared; 
        ! organize in 2D processor mesh: j + ny*(i-1))
        split_key = split_key +     npy*l_lim(2,1)/g_lim(2,1)
     .                        + npy*npx*l_lim(1,1)/g_lim(1,1)

c     End program

      end subroutine selectRank_Z

c     selectRank_XY
c     #####################################################################
      subroutine selectRank_XY(split_key)

c     ---------------------------------------------------------------------
c     Groups processors based on shared X position in processor mesh.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: split_key

c     Local variables

        integer :: npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npx = g_mpi_def%mx
        npy = g_mpi_def%my
        npz = g_mpi_def%mz

        !Group processors based on location in x direction
        !(y and z directions are shared; 
        ! organize in 1D processor mesh: i
        split_key = split_key + npz*l_lim(3,1)/g_lim(3,1)

c     End program

      end subroutine selectRank_XY

c     selectRank_XZ
c     #####################################################################
      subroutine selectRank_XZ(split_key)

c     ---------------------------------------------------------------------
c     Groups processors based on shared Y position in processor mesh.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: split_key

c     Local variables

        integer :: npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npx = g_mpi_def%mx
        npy = g_mpi_def%my
        npz = g_mpi_def%mz

        !Group processors based on location in x direction
        !(x and z directions are shared; 
        ! organize in 1D processor mesh: j
        split_key = split_key + npy*l_lim(2,1)/g_lim(2,1)

c     End program

      end subroutine selectRank_XZ

c     selectRank_YZ
c     #####################################################################
      subroutine selectRank_YZ(split_key)

c     ---------------------------------------------------------------------
c     Groups processors based on shared X position in processor mesh.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(INOUT) :: split_key

c     Local variables

        integer :: npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mpi_def%ilo(1)
        l_lim(2,0) = g_mpi_def%jlo(1)
        l_lim(3,0) = g_mpi_def%klo(1)
        l_lim(1,1) = g_mpi_def%ihi(1)
        l_lim(2,1) = g_mpi_def%jhi(1)
        l_lim(3,1) = g_mpi_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mpi_def%nxgl(1)
        g_lim(2,1) = g_mpi_def%nygl(1)
        g_lim(3,1) = g_mpi_def%nzgl(1)

c     Select rank

        npx = g_mpi_def%mx
        npy = g_mpi_def%my
        npz = g_mpi_def%mz

        !Group processors based on location in x direction
        !(y and z directions are shared; 
        ! organize in 1D processor mesh: i
        split_key = split_key + npx*l_lim(1,1)/g_lim(1,1)

c     End program

      end subroutine selectRank_YZ

c     checkProc
c     #####################################################################
      function checkProc(BC,dim,loc,llim,glim) result(include_proc)

c     ---------------------------------------------------------------------
c     Checks whether a processor should be included in an MPI_COMM according
c     to the boundary type BC, and the local vs. global limits
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: BC,dim,loc,llim(3,0:1),glim(3,0:1)
        logical    :: include_proc

c     Local variables

        integer    :: ibc

c     Begin program

        ibc(dim,loc) = (1+loc)+2*(dim-1)

        include_proc= (llim(dim,loc)==glim(dim,loc))
     .                 .and.bcond(ibc(dim,loc))==BC

      end function checkProc

c     fillLocalVec
c     #################################################################
      subroutine fillLocalVec(dactx,array,x)

c     -----------------------------------------------------------------
c     Transfers outer layer of values from array to vector x.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(petsc_da_ctx),pointer :: dactx
        real(8) :: array(dactx%ilm:dactx%ihp
     .                  ,dactx%jlm:dactx%jhp 
     .                  ,dactx%klm:dactx%khp)
     .            ,x    (dactx%gxs:dactx%gxe 
     .                  ,dactx%gys:dactx%gye 
     .                  ,dactx%gzs:dactx%gze)

c     Local variables

c     Begin program

        x(dactx%xs
     .   ,dactx%ys:dactx%ye 
     .   ,dactx%zs:dactx%ze)
     .                       = array(dactx%il
     .                              ,dactx%jl:dactx%jh 
     .                              ,dactx%kl:dactx%kh)

        if (dactx%xe > dactx%xs)
     .    x(dactx%xe 
     .     ,dactx%ys:dactx%ye 
     .     ,dactx%zs:dactx%ze) =
     .                         array(dactx%ih
     .                              ,dactx%jl:dactx%jh 
     .                              ,dactx%kl:dactx%kh)

        x(dactx%xs:dactx%xe 
     .   ,dactx%ys
     .   ,dactx%zs:dactx%ze) =
     .                         array(dactx%il:dactx%ih
     .                              ,dactx%jl
     .                              ,dactx%kl:dactx%kh)

        if (dactx%ye > dactx%ys)
     .    x(dactx%xs:dactx%xe 
     .     ,dactx%ye 
     .     ,dactx%zs:dactx%ze) =
     .                         array(dactx%il:dactx%ih
     .                              ,dactx%jh 
     .                              ,dactx%kl:dactx%kh)

        x(dactx%xs:dactx%xe 
     .   ,dactx%ys:dactx%ye 
     .   ,dactx%zs              ) =
     .                         array(dactx%il:dactx%ih
     .                              ,dactx%jl:dactx%jh 
     .                              ,dactx%kl              )

        if (dactx%ze > dactx%zs)
     .    x(dactx%xs:dactx%xe 
     .     ,dactx%ys:dactx%ye 
     .     ,dactx%ze              ) =
     .                         array(dactx%il:dactx%ih
     .                              ,dactx%jl:dactx%jh 
     .                              ,dactx%kh              )

cc        x(dactx%xs:dactx%xe 
cc     .   ,dactx%ys:dactx%ye 
cc     .   ,dactx%zs:dactx%ze) =
cc     .                         array(dactx%il:dactx%ih
cc     .                              ,dactx%jl:dactx%jh 
cc     .                              ,dactx%kl:dactx%kh)

c     End program

      end subroutine fillLocalVec

c     emptyLocalVec
c     #################################################################
      subroutine emptyLocalVec(dactx,array,x)

c     -----------------------------------------------------------------
c     Transfers ghost cells from x to array.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(petsc_da_ctx),pointer :: dactx
 
        real(8) :: array(dactx%ilm:dactx%ihp
     .                  ,dactx%jlm:dactx%jhp 
     .                  ,dactx%klm:dactx%khp)
     .            ,x    (dactx%gxs:dactx%gxe 
     .                  ,dactx%gys:dactx%gye 
     .                  ,dactx%gzs:dactx%gze)

c     Local variables

c     Begin program

        array(dactx%lgxs
     .       ,dactx%lgys:dactx%lgye
     .       ,dactx%lgzs:dactx%lgze) =
     .                                 x(dactx%gxs
     .                                  ,dactx%gys:dactx%gye 
     .                                  ,dactx%gzs:dactx%gze)

        if (dactx%lgxe > dactx%lgxs) 
     .    array(dactx%lgxe
     .         ,dactx%lgys:dactx%lgye
     .         ,dactx%lgzs:dactx%lgze) =
     .                                 x(dactx%gxe 
     .                                  ,dactx%gys:dactx%gye 
     .                                  ,dactx%gzs:dactx%gze)

        array(dactx%lgxs:dactx%lgxe
     .       ,dactx%lgys
     .       ,dactx%lgzs:dactx%lgze) =
     .                                 x(dactx%gxs:dactx%gxe 
     .                                  ,dactx%gys
     .                                  ,dactx%gzs:dactx%gze)

        if (dactx%lgye > dactx%lgys) 
     .    array(dactx%lgxs:dactx%lgxe
     .         ,dactx%lgye
     .         ,dactx%lgzs:dactx%lgze) =
     .                                 x(dactx%gxs:dactx%gxe 
     .                                  ,dactx%gye 
     .                                  ,dactx%gzs:dactx%gze)

        array(dactx%lgxs:dactx%lgxe
     .       ,dactx%lgys:dactx%lgye
     .       ,dactx%lgzs                ) =
     .                                 x(dactx%gxs:dactx%gxe 
     .                                  ,dactx%gys:dactx%gye 
     .                                  ,dactx%gzs               )

        if (dactx%lgze > dactx%lgzs) 
     .    array(dactx%lgxs:dactx%lgxe
     .         ,dactx%lgys:dactx%lgye
     .         ,dactx%lgze                ) =
     .                                 x(dactx%gxs:dactx%gxe 
     .                                  ,dactx%gys:dactx%gye 
     .                                  ,dactx%gze               )

cc        array(dactx%lgxs:dactx%lgxe
cc     .       ,dactx%lgys:dactx%lgye
cc     .       ,dactx%lgzs:dactx%lgze) =
cc     .                                 x(dactx%gxs:dactx%gxe 
cc     .                                  ,dactx%gys:dactx%gye 
cc     .                                  ,dactx%gzs:dactx%gze)

c     End program                     

      end subroutine emptyLocalVec

c     fchkerrq
c     #################################################################
      subroutine fchkerrq(ierr,routine)

        implicit none

c     Call variables

        integer    :: ierr
        character(*) :: routine

c     Begin program

        if (ierr /= 0) call pstop(routine,'PETSc error')

c     End program

      end subroutine fchkerrq

#else

        integer       :: np=1,my_rank=0
        character(80) :: messg

      contains

#endif

c     pstop
c     ################################################################
      subroutine pstop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        character(*)  :: routine, message

c     Local variables

        integer :: ierr
        
c     Begin program

#if defined(samrai)
c     Call SAMRAI's error handler
        interface
           subroutine pstop_samrai(fun,msg) 
     .                bind(C, name="pstop_samrai")
              use iso_c_binding, only: c_char
              character(kind=c_char) :: fun(*), msg(*)
           end subroutine pstop_samrai
        end interface
        call pstop_samrai(trim(routine),trim(message))
#else
c     Call local error handler

c$$$      call sstop(my_rank,routine,message)

        if (my_rank == 0) then
          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)
        endif

#if defined(petsc)
        call PetscFinalize(ierr)
#endif

        call exit(FATAL)  !Fatal termination

#endif

      end subroutine pstop

c     psum
c     ###################################################################
      function psum(vec1,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel sum of elements of a 1d-array using communicator
c     mpi_comm. On input:
c       * vec1: vector
c       * mpi_comm: MPI communicator
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer,optional :: mpi_comm
      real(8)    :: vec1(:),psum

c     Local variables

      real(8)    :: lsum
      integer    :: mpic

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpic = mpi_comm
      else
        mpic = MPI_COMM_WORLD
      endif

      lsum = sum(vec1)
      call MPI_Allreduce(lsum,psum,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,mpic,mpierr)
#else
      psum = sum(vec1)
#endif

      end function psum

c     ipsum
c     ###################################################################
      function ipsum(vec1,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel sum of elements of a 1d-array using communicator
c     mpi_comm. On input:
c       * vec1: vector
c       * mpi_comm: MPI communicator
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer,optional :: mpi_comm
      integer :: vec1(:),ipsum

c     Local variables

      integer    :: mpic,lsum

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpic = mpi_comm
      else
        mpic = MPI_COMM_WORLD
      endif

      lsum = sum(vec1)
      call MPI_Allreduce(lsum,ipsum,1,MPI_INT
     .                  ,MPI_SUM,mpic,mpierr)
#else
      ipsum = sum(vec1)
#endif

      end function ipsum

c     dot
c     ###################################################################
      function dot1(ntot,vec1,vec2) result(dot)

c     -------------------------------------------------------------------
c     Performs parallel scalar product (vec1,vec2). In call:
c       * ntot: vector dimension
c       * vec1: vector, first term of scalar product.
c       * vec2: vector, second term of scalar product.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ntot
      real(8)    :: vec1(ntot),vec2(ntot),dot

c     Local variables

      real(8)    :: ldot

c     Begin program

      dot = dot_product(vec1,vec2)

#if defined(petsc)
      ldot = dot
      call MPI_Allreduce(ldot,dot,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      end function dot1

c     dot2 (Needed for FPA)
c     ###################################################################
      function dot2(vec1,vec2)

c     -------------------------------------------------------------------
c     Performs parallel scalar product (vec1,vec2). In call:
c       * vec1: vector, first term of scalar product.
c       * vec2: vector, second term of scalar product.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dot2
      real(8),INTENT(IN) :: vec1(:),vec2(:)

c     Local variables

      integer :: nn,lnn
      real(8) :: ldot

c     Begin program

      dot2 = dot_product(vec1,vec2)

#if defined(petsc)
      ldot = dot2
      call MPI_Allreduce(ldot,dot2,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      end function dot2

c     dot3
c     ###################################################################
      function dot3(vec1,vec2,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel scalar product (vec1,vec2). In call:
c       * vec1: vector, first term of scalar product.
c       * vec2: vector, second term of scalar product.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dot3
      real(8),INTENT(IN) :: vec1(:),vec2(:)
      integer,optional :: mpi_comm

c     Local variables

      integer :: mpicomm
      real(8) :: ldot

c     Begin program

      dot3 = dot_product(vec1,vec2)

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif

      ldot = dot3
      call MPI_Allreduce(ldot,dot3,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,mpicomm,mpierr)
#endif

      end function dot3

c     rms
c     ###################################################################
      function rms(vec1,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel root mean square. In call:
c       * vec1: vector
c       * mpi_comm: MPI communicator
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: rms
      real(8),INTENT(IN) :: vec1(:)
      integer,optional :: mpi_comm

c     Local variables

      integer :: nn,lnn,mpicomm
      real(8) :: lrms

c     Begin program

      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
#if defined(petsc)
        mpicomm = MPI_COMM_WORLD
#endif
      endif

      nn = size(vec1)
      rms = dot_product(vec1,vec1)

#if defined(petsc)
      if (.not.asm) then
        lrms = rms
        call MPI_Allreduce(lrms,rms,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,mpicomm,mpierr)

        lnn = nn
        call MPI_Allreduce(lnn,nn,1,MPI_INTEGER
     .                    ,MPI_SUM,mpicomm,mpierr)
      endif
#endif

      rms=sqrt(rms/nn)

      end function rms

c     pmax
c     ###################################################################
      function pmax(val,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel maximum value.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: pmax,val
      integer,optional :: mpi_comm

c     Local variables

      integer :: mpicomm

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif

      call MPI_Allreduce(val,pmax,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_MAX,mpicomm,mpierr)
#else
      pmax = val
#endif

      end function pmax

c     pmin
c     ###################################################################
      function pmin(val,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel maximum value.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: pmin,val
      integer,optional :: mpi_comm

c     Local variables

      integer :: mpicomm

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif

      call MPI_Allreduce(val,pmin,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_MIN,mpicomm,mpierr)
#else
      pmin = val
#endif

      end function pmin

c     ipmin
c     ###################################################################
      function ipmin(val,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel maximum value.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ipmin,val
      integer,optional :: mpi_comm

c     Local variables

      integer :: mpicomm

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif

      call MPI_Allreduce(val,pmin,1,MPI_INT,MPI_MIN,mpicomm,mpierr)
#else
      pmin = val
#endif

      end function ipmin
      
c     isProc
c     #################################################################
      function isProc(proc_id)

        implicit none

c     Call variables

        integer    :: proc_id
        logical    :: isProc

c     Local variables

c     Begin program

        isProc = (my_rank == proc_id)

c     End program

      end function isProc

c     inProc
c     #################################################################
      function inProc(gl_def,igl,jgl,kgl,igx,igy,igz)

        implicit none

c     Call variables

        integer :: igl,jgl,kgl,igx,igy,igz
        logical :: inProc

        type(grid_mg_def),pointer  :: gl_def

c     Local variables

c     Begin program

        inProc =
     .       (igl<0.or.(igl>=gl_def%ilo(igx).and.igl<=gl_def%ihi(igx)))
     .  .and.(jgl<0.or.(jgl>=gl_def%jlo(igy).and.jgl<=gl_def%jhi(igy)))
     .  .and.(kgl<0.or.(kgl>=gl_def%klo(igz).and.kgl<=gl_def%khi(igz)))

c     End program

      end function inProc

c     find_global_field_vec
c     ##############################################################
      subroutine find_global_field_vec(lf,gf,dump,map_bc,no_bc,mpi_comm)

      use xdraw_io

      implicit none

c     --------------------------------------------------------------
c     Performs parallel reduction to find global field INCLUDING
c     ghost cells.
c     --------------------------------------------------------------

c     Call variables

      real(8),dimension(0:,0:,0:,:) :: lf,gf

      logical,optional :: dump,map_bc,no_bc

      integer,optional :: mpi_comm

c     Local variables

      integer :: i,j,k,dl(4),dg(4),lsize,gsize,ibc,ieq,dim,loc
      real(8),pointer,dimension(:) :: local_x,global_x

      logical :: dmp,mapbc,nobc

      integer :: mpicomm,npxl,npyl,npzl,npl

c     Begin program

      if (PRESENT(dump)) then
        dmp = dump
      else
        dmp = .false.
      endif

      if (PRESENT(map_bc)) then
        mapbc = map_bc
      else
        mapbc = .false.
      endif

      if (PRESENT(no_bc)) then
        nobc = no_bc
      else
        nobc = .true.
      endif

      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
#if defined(petsc)
        mpicomm = MPI_COMM_WORLD
#endif
      endif

c     Calculate array sizes

      dl = shape(lf)

      dg = shape(gf)

      lsize = dl(1)*dl(2)*dl(3)*dl(4)

      npxl = (dg(1)-2)/(dl(1)-2)
      npyl = (dg(2)-2)/(dl(2)-2)
      npzl = (dg(3)-2)/(dl(3)-2)

      gsize = dl(1)*npxl*dl(2)*npyl*dl(3)*npzl*dl(4)

#if defined(petsc)
      call MPI_Comm_size(mpicomm,npl,mpierr)

      if (npl /= npxl*npyl*npzl) then
        call pstop('find_global_field_vec'
     .            ,'Processor number does not agree')
      endif
#endif

c     Communicate

#if defined(petsc)
      allocate(local_x(lsize),global_x(gsize))

      local_x = reshape(lf, (/ lsize /))

      call MPI_Allgather(local_x ,lsize,MPI_DOUBLE_PRECISION
     .                  ,global_x,lsize,MPI_DOUBLE_PRECISION
     .                  ,mpicomm,mpierr)

      do k=1,npzl
        do j=1,npyl
          do i=1,npxl

            gf((dl(1)-2)*(i-1):(dl(1)-2)*i+1
     .        ,(dl(2)-2)*(j-1):(dl(2)-2)*j+1
     .        ,(dl(3)-2)*(k-1):(dl(3)-2)*k+1,1:dg(4))
     .        = reshape(global_x(1+lsize*(          (i-1)
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1))
     .                          :  lsize*(           i
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1)) )
     .                 ,dl )
          enddo
        enddo
      enddo

      deallocate(local_x,global_x)
#else
      gf = lf
#endif

c     Synchronize global periodic boundaries

      if (mapbc) then
        do dim =1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (bcond(ibc) == PER) then
              do ieq=1,dg(4)
                if (ieq == dim) cycle
                call per_bc(ibc,gf(:,:,:,ieq))
              enddo
            endif
          enddo
        enddo
      elseif (.not.nobc) then
        do ibc =1,6
          if (bcond(ibc) == PER) then
            do ieq=1,dg(4)
              call per_bc(ibc,gf(:,:,:,ieq))
            enddo
          endif
        enddo
      endif

c     Dump global array

      if (dmp) then
        if (my_rank == 0) then
          open(unit=110,file='debug.bin'
     .        ,form='unformatted',status='replace')
          do i=1,dg(4)
            call contour(gf(0:dg(1)-1,0:dg(2)-1,1,i),dg(1),dg(2)
     .                ,0d0,1d0,0d0,1d0,i-1,110)
          enddo
cc          do i=1,dl(4)
cc            call contour(lf(0:dl(1)-1,0:dl(2)-1,1,i),dl(1),dl(2)
cc     .                ,0d0,xmax,0d0,ymax,i-1,110)
cc          enddo
          close(110)
        endif

        call pstop('find_global_field_vec','Dumping global arrays...')
      endif

c     End program

      end subroutine find_global_field_vec

c     find_global_field_vec_nobc
c     ##############################################################
      subroutine find_global_field_vec_nobc(lf,gf,mpi_comm)

      implicit none

c     --------------------------------------------------------------
c     Performs parallel reduction to find global field EXCLUDING
c     ghost cells.
c     --------------------------------------------------------------

c     Call variables

      integer,optional :: mpi_comm
      real(8),dimension(:,:,:,:) :: lf,gf

c     Local variables

      integer :: i,j,k,dl(4),dg(4),lsize,gsize,ibc,ieq,dim,loc
      real(8),pointer,dimension(:) :: local_x,global_x

      logical :: dmp,mapbc

      integer :: mpicomm,npxl,npyl,npzl,npl

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif
#endif

c     Calculate array sizes

      dl = shape(lf)

      dg = shape(gf)

      lsize = dl(1)*dl(2)*dl(3)*dl(4)

cc      gsize = dl(1)*npx*dl(2)*npy*dl(3)*npz*dl(4)
      gsize = dg(1)*dg(2)*dg(3)*dg(4)

      npxl = dg(1)/dl(1)
      npyl = dg(2)/dl(2)
      npzl = dg(3)/dl(3)

#if defined(petsc)
      call MPI_Comm_size(mpicomm,npl,mpierr)

      if (npl /= npxl*npyl*npzl) then
        call pstop('find_global_field_vec_nobc'
     .            ,'Processor number does not agree')
      endif
#endif

c     Communicate

#if defined(petsc)
      allocate(local_x(lsize),global_x(gsize))

      local_x = reshape(lf, (/ lsize /))

      call MPI_Allgather(local_x ,lsize,MPI_DOUBLE_PRECISION
     .                  ,global_x,lsize,MPI_DOUBLE_PRECISION
     .                  ,mpicomm,mpierr)

      do k=1,npzl
        do j=1,npyl
          do i=1,npxl

            gf(dl(1)*(i-1)+1:dl(1)*i
     .        ,dl(2)*(j-1)+1:dl(2)*j
     .        ,dl(3)*(k-1)+1:dl(3)*k,1:dg(4))
     .        = reshape(global_x(1+lsize*(          (i-1)
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1))
     .                          :  lsize*(           i
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1)) )
     .                 ,dl )
          enddo
        enddo
      enddo

      deallocate(local_x,global_x)
#else
      gf = lf
#endif

c     End program

      end subroutine find_global_field_vec_nobc

c     find_global_field_scl
c     ##############################################################
      subroutine find_global_field_scl(lf,gf,dump,mpi_comm,no_bc)

      use xdraw_io

      implicit none

c     --------------------------------------------------------------
c     Performs parallel reduction to find global field INCLUDIN
c     ghost cells.
c     --------------------------------------------------------------

c     Call variables

      real(8),dimension(0:,0:,0:) :: lf,gf

      logical,optional :: dump,no_bc
      integer,optional :: mpi_comm

c     Local variables

      integer :: i,j,k,dl(3),dg(3),lsize,gsize,ibc
      real(8),pointer,dimension(:) :: local_x,global_x

      logical :: dmp,nobc

      integer :: mpicomm,npxl,npyl,npzl,npl

c     Begin program

      if (PRESENT(dump)) then
        dmp=dump
      else
        dmp=.false.
      endif

      if (PRESENT(no_bc)) then
        nobc=no_bc
      else
        nobc=.false.
      endif

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif
#endif

c     Calculate array sizes

      dl = shape(lf)

      dg = shape(gf)

      lsize = dl(1)*dl(2)*dl(3)

      npxl = (dg(1)-2)/(dl(1)-2)
      npyl = (dg(2)-2)/(dl(2)-2)
      npzl = (dg(3)-2)/(dl(3)-2)

      gsize = dl(1)*npxl*dl(2)*npyl*dl(3)*npzl

#if defined(petsc)
      call MPI_Comm_size(mpicomm,npl,mpierr)

      if (npl /= npxl*npyl*npzl) then
        call pstop('find_global_field_scl'
     .            ,'Processor number does not agree')
      endif
#endif

c     Communicate

#if defined(petsc)
      allocate(local_x(lsize),global_x(gsize))

      local_x = reshape(lf, (/ lsize /))

      call MPI_Allgather(local_x ,lsize,MPI_DOUBLE_PRECISION
     .                  ,global_x,lsize,MPI_DOUBLE_PRECISION
     .                  ,mpicomm,mpierr)

      do k=1,npzl
        do j=1,npyl
          do i=1,npxl
            gf((dl(1)-2)*(i-1):(dl(1)-2)*i+1
     .        ,(dl(2)-2)*(j-1):(dl(2)-2)*j+1
     .        ,(dl(3)-2)*(k-1):(dl(3)-2)*k+1)
     .        = reshape(global_x(1+lsize*(          (i-1)
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1))
     .                          :  lsize*(           i
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1)) )
     .                 ,dl )
          enddo
        enddo
      enddo

      deallocate(local_x,global_x)

      !Synchronize global periodic boundaries
      do ibc=1,6
        if (bcond(ibc)==PER.and.(.not.nobc)) call per_bc(ibc,gf)
      enddo
#else
      gf = lf
#endif

c     Dump global array

      if (dmp) then
        if (my_rank == 0) then
          open(unit=110,file='debug.bin'
     .        ,form='unformatted',status='replace')
          call contour(gf(0:dg(1)-1,0:dg(2)-1,1),dg(1),dg(2)
     .                ,0d0,1d0,0d0,1d0,0,110)
          close(110)
        endif

        call pstop('find_global_field_scl','Dumping global arrays...')
      endif

c     End program

      end subroutine find_global_field_scl

c     find_global_field_scl_nobc
c     ##############################################################
      subroutine find_global_field_scl_nobc(lf,gf,mpi_comm)

      implicit none

c     --------------------------------------------------------------
c     Performs parallel reduction to find global field EXCLUDING
c     ghost cells.
c     --------------------------------------------------------------

c     Call variables

      real(8),dimension(:,:,:) :: lf,gf
      integer,optional :: mpi_comm

c     Local variables

      integer :: i,j,k,dl(3),dg(3),lsize,gsize,ibc
      real(8),pointer,dimension(:) :: local_x,global_x
      integer :: mpicomm,npxl,npyl,npzl,npl

c     Begin program

#if defined(petsc)
      if (PRESENT(mpi_comm)) then
        mpicomm = mpi_comm
      else
        mpicomm = MPI_COMM_WORLD
      endif
#endif

c     Calculate array sizes

      dl = shape(lf)

      dg = shape(gf)

      lsize = dl(1)*dl(2)*dl(3)

cc      gsize = dl(1)*npx*dl(2)*npy*dl(3)*npz
      gsize = dg(1)*dg(2)*dg(3)

      npxl = dg(1)/dl(1)
      npyl = dg(2)/dl(2)
      npzl = dg(3)/dl(3)

#if defined(petsc)
      call MPI_Comm_size(mpicomm,npl,mpierr)

      if (npl /= npxl*npyl*npzl) then
        call pstop('find_global_field_scl_nobc'
     .            ,'Processor number does not agree')
      endif
#endif

c     Communicate

#if defined(petsc)
      allocate(local_x(lsize),global_x(gsize))

      local_x = reshape(lf, (/ lsize /))

      call MPI_Allgather(local_x ,lsize,MPI_DOUBLE_PRECISION
     .                  ,global_x,lsize,MPI_DOUBLE_PRECISION
     .                  ,mpicomm,mpierr)

      do k=1,npzl
        do j=1,npyl
          do i=1,npxl
            gf(dl(1)*(i-1)+1:dl(1)*i
     .        ,dl(2)*(j-1)+1:dl(2)*j
     .        ,dl(3)*(k-1)+1:dl(3)*k)
     .        = reshape(global_x(1+lsize*(          (i-1)
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1))
     .                          :  lsize*(           i
     .                                   +npxl     *(j-1)
     .                                   +npxl*npyl*(k-1)) )
     .                 ,dl )
          enddo
        enddo
      enddo

      deallocate(local_x,global_x)
#else
      gf = lf
#endif

c     End program

      end subroutine find_global_field_scl_nobc

c     find_local_field_vec
c     ##############################################################
      subroutine find_local_field_vec(g_def,igrid,gf,lf)

      implicit none

c     --------------------------------------------------------------
c     Performs restriction from global to local field.
c     --------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8),dimension(0:,0:,0:,:) :: lf,gf

c     Local variables

      integer :: i,j,k,dl(4),igl(3),ilc(3)

c     Begin program

c     Calculate array sizes

      dl = shape(lf)

c     Transfer field to local mesh

      do k=0,dl(3)-1
        do j=0,dl(2)-1
          do i=0,dl(1)-1

            !Find global limits
            ilc = (/i,j,k/)
            igl = globalIndex(g_def,igrid,ilc)

            !Map field
            lf(i,j,k,:) = gf(igl(1),igl(2),igl(3),:)
          enddo
        enddo
      enddo

c     End program

      end subroutine find_local_field_vec

c     find_local_field_scl
c     ##############################################################
      subroutine find_local_field_scl(g_def,igrid,gf,lf)

      implicit none

c     --------------------------------------------------------------
c     Performs parallel reduction to find global field.
c     --------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8),dimension(0:,0:,0:) :: lf,gf

c     Local variables

      integer :: i,j,k,dl(3),igl(3)

c     Begin program

c     Calculate array sizes

      dl = shape(lf)

c     Transfer field to local mesh

      do k=0,dl(3)-1
        do j=0,dl(2)-1
          do i=0,dl(1)-1

            !Find global limits
            igl = globalIndex(g_def,igrid,(/i,j,k/))

            !Map field
            lf(i,j,k) = gf(igl(1),igl(2),igl(3))
          enddo
        enddo
      enddo

c     End program

      end subroutine find_local_field_scl

c$$$c     find_global_vec
c$$$c     ##############################################################
c$$$      subroutine find_global_vec(neq,lv,gv,mpi_comm)
c$$$
c$$$      implicit none
c$$$
c$$$c     --------------------------------------------------------------
c$$$c     Performs parallel reduction to find global field.
c$$$c     --------------------------------------------------------------
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: neq
c$$$      integer,optional :: mpi_comm
c$$$      real(8),dimension(:) :: lv,gv
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: i,j,k,dl(4),dg(4),lsize,gsize,ibc,ieq,dim,loc
c$$$      real(8),pointer,dimension(:) :: local_x,global_x
c$$$
c$$$      logical :: dmp,mapbc
c$$$
c$$$      integer :: mpicomm,npxl,npyl,npzl,npl
c$$$
c$$$c     Begin program
c$$$
c$$$#if defined(petsc)
c$$$      if (PRESENT(mpi_comm)) then
c$$$        mpicomm = mpi_comm
c$$$      else
c$$$        mpicomm = MPI_COMM_WORLD
c$$$      endif
c$$$#endif
c$$$
c$$$c     Calculate array sizes
c$$$
c$$$      dl = shape(lf)
c$$$
c$$$      dg = shape(gf)
c$$$
c$$$      lsize = dl(1)*dl(2)*dl(3)*dl(4)
c$$$
c$$$cc      gsize = dl(1)*npx*dl(2)*npy*dl(3)*npz*dl(4)
c$$$      gsize = dg(1)*dg(2)*dg(3)*dg(4)
c$$$
c$$$      npxl = dg(1)/dl(1)
c$$$      npyl = dg(2)/dl(2)
c$$$      npzl = dg(3)/dl(3)
c$$$
c$$$#if defined(petsc)
c$$$      call MPI_Comm_size(mpicomm,npl,mpierr)
c$$$
c$$$      if (npl /= npxl*npyl*npzl) then
c$$$        call pstop('find_global_field_vec_nobc'
c$$$     .            ,'Processor number does not agree')
c$$$      endif
c$$$#endif
c$$$
c$$$c     Communicate
c$$$
c$$$#if defined(petsc)
c$$$      allocate(local_x(lsize),global_x(gsize))
c$$$
c$$$      local_x = reshape(lf, (/ lsize /))
c$$$
c$$$      call MPI_Allgather(local_x ,lsize,MPI_DOUBLE_PRECISION
c$$$     .                  ,global_x,lsize,MPI_DOUBLE_PRECISION
c$$$     .                  ,mpicomm,mpierr)
c$$$
c$$$
c$$$      do k=1,npzl
c$$$        do j=1,npyl
c$$$          do i=1,npxl
c$$$
c$$$            gf(dl(1)*(i-1)+1:dl(1)*i
c$$$     .        ,dl(2)*(j-1)+1:dl(2)*j
c$$$     .        ,dl(3)*(k-1)+1:dl(3)*k,1:dg(4))
c$$$     .        = reshape(global_x(1+lsize*(          (i-1)
c$$$     .                                   +npxl     *(j-1)
c$$$     .                                   +npxl*npyl*(k-1))
c$$$     .                          :  lsize*(           i
c$$$     .                                   +npxl     *(j-1)
c$$$     .                                   +npxl*npyl*(k-1)) )
c$$$     .                 ,dl )
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      deallocate(local_x,global_x)
c$$$#else
c$$$      gf = lf
c$$$#endif
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine find_global_vec

c     dumpGlobalField
c     #################################################################
      subroutine dumpGlobalField(g_def,file,nx,ny,nz,nxg,nyg,nzg,arr
     .                          ,flg,no_bc)

      use io, ONLY: find_unit

      implicit none

c     -----------------------------------------------------------------
c     Dump XDRAW contour file of global field "arr". Variable "flg"
c     controls IO mode:
c       * flg = 0: init mode
c       * flg = 1: append mode
c     -----------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: nx,ny,nz,nxg,nyg,nzg,flg,mpi_comm
      real(8) :: arr(0:nx+1,0:ny+1,0:nz+1)
      character(*) :: file
      logical, optional :: no_bc

c     Local variables

      integer :: unit
      real(8),allocatable,dimension(:,:,:) :: dbg
      real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
      logical :: nobc

c     Begin program

      if (PRESENT(no_bc)) then
        nobc = no_bc
      else
        nobc = .false.
      endif

      xmin = g_def%gxmin
      xmax = g_def%gxmax

      ymin = g_def%gymin
      ymax = g_def%gymax

      zmin = g_def%gzmin
      zmax = g_def%gzmax

      !Open .bin file
      if (my_rank == 0) then
        if (flg == 0) then
          unit = find_unit(12345)
          open(unit=unit,file=trim(file),form='unformatted'
     .        ,status='unknown')
        elseif (flg == 1) then
          unit = find_unit(12345)
          open(unit=unit,file=trim(file),form='unformatted'
     .        ,access='append',status='old')
        endif
      endif
        
      allocate(dbg(0:nxg+1,0:nyg+1,0:nzg+1))

      if (nobc) then
        call find_global_nobc(arr(1:nx ,1:ny ,1:nz )
     .                       ,dbg(1:nxg,1:nyg,1:nzg)
#if defined(petsc)
     .                       ,mpi_comm=g_def%mpi_comm
#endif
     .                       )
        if (my_rank == 0) then
          call contour(dbg(1:nxg,1:nyg,1),nxg,nyg
     .                ,xmin,xmax,ymin,ymax,flg,unit)
        endif
      else
        call find_global(arr,dbg
#if defined(petsc)
     .                  ,mpi_comm=g_def%mpi_comm
#endif
     .                  ,no_bc=.true.)
        if (my_rank == 0) then
          call contour(dbg(:,:,1),nxg+2,nyg+2
     .                ,xmin,xmax,ymin,ymax,flg,unit)
        endif
      endif

      deallocate(dbg)

      if (my_rank == 0) close(unit)

      end subroutine dumpGlobalField

c     per_BC
c     #################################################################
      subroutine per_BC(ibc,array)
c     -----------------------------------------------------------------
c     Imposes periodic BC. On input:
c        * ibc -> face identifier (1 to 6)
c        * array -> 3D array to impose BC's on.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ibc
      real(8) :: array(0:,0:,0:)

c     Local variables

      integer :: nx,ny,nz

c     Begin program

      nx = size(array,1)-2
      ny = size(array,2)-2
      nz = size(array,3)-2

      select case (ibc)
      case (1)                  !x0
        array(0   ,:,:)=array(nx,:,:)
      case (2)                  !x1
        array(nx+1,:,:)=array(1 ,:,:)
      case (3)                  !y0
        array(:,0   ,:)=array(:,ny,:)
      case (4)                  !y1
        array(:,ny+1,:)=array(:,1 ,:)
      case (5)                  !z0
        array(:,:,0   )=array(:,:,nz)
      case (6)                  !z1
        array(:,:,nz+1)=array(:,:,1 )
      case default
        call pstop('per_BC'
     .            ,'Boundary '//int2char(ibc)//' non existent')
      end select

c     End program

      end subroutine per_BC

#if defined(petsc)
c     processorAlloc
c     ######################################################################
      subroutine processorAlloc(mpi_comm,nx1,ny1,nz1,npx,npy,npz)

c     ----------------------------------------------------------------------
c     Finds processor allocation (npx,npy,npz; defined in module) in the
c     three logical directions based on global dimensions
c     (nx1,nx2,nx3). Only returns processor allocation for directions
c     not specified earlier by user (i.e., when np* is zero).
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer,INTENT(IN) :: nx1,ny1,nz1,mpi_comm
      integer,INTENT(INOUT) :: npx,npy,npz

c     Local variables

      integer :: nx,ny,nz,sum_exp,nd,navg,exp(3),npt(3),nprocs

      logical :: proc_debug=.false.

c     Begin program

      nx = nx1
      ny = ny1
      nz = nz1

      exp = 0

      call initMPI(mpi_comm,nprocs,my_rank)

      proc_debug = proc_debug.and.(my_rank == 0)

c     Eliminate specified directions

      if (npx /= 0) then
        nprocs = max(nprocs/npx,1)
        nx = 1
      endif

      if (npy /= 0) then
        nprocs = max(nprocs/npy,1)
        ny = 1
      endif

      if (npz /= 0) then
        nprocs = max(nprocs/npz,1)
        nz = 1
      endif

      if (nx == 1 .and. ny == 1 .and. nz == 1) then
        if (nprocs == 1) then  !Compatible specification of procs. per direction
          npx = max(npx,1) ; npy = max(npy,1) ; npz = max(npz,1)
          return
        else
          messg = 'np='//trim(int2char(np))//' /= npx*npy*npz='
     $           //trim(int2char(npx*npy*npz))//' ... Aborting'
          call pstop('processorAlloc',messg)
        endif
      endif

c     Check for unallocated procs

      sum_exp = floor(log(1d0*nprocs)/log(2d0))

      if (sum_exp < 1) then   !No processors left to allocate
        if (npx == 0) npx = 1
        if (npy == 0) npy = 1
        if (npz == 0) npz = 1
        return
      endif

      if (2**sum_exp /= nprocs) then
        messg = 'Number of processors unsuitable for MG'
        call pstop('processorAlloc',messg)
      endif

c     Find dimensionality

      nd = 3
      if (nx == 1) nd = nd-1
      if (ny == 1) nd = nd-1
      if (nz == 1) nd = nd-1

      if (nd == 0) then
        messg = 'No available dimensions!'
        call pstop('processorAlloc',messg)
      endif

      !Find exponents
      navg = int(((1d0*nx)*(1d0*ny)*(1d0*nz)/nprocs)**(1d0/nd)) !Average number of points per proc per dim
      navg = 2**(floor(log(1d0*navg))/log(2d0))     !Find nearest power of 2 from below

      npt = (/nx,ny,nz/)

      if (proc_debug) 
     .   write (*,*) 'DIAG -- procAlloc, npt=',npt
     .            ,' nprocs=',nprocs,' navg=',navg

      exp = log(max(1d0*npt/navg,1d0))/log(2d0)     !Find exponents of power of 2 per dim
      exp = max(exp - (sum(exp)-sum_exp)/nd,0)      !Subtract excess procs

      if (proc_debug) write (*,*) 'DIAG -- procAlloc, exp_0=',exp

      !Consistency check (1st layer)
      if (sum(exp) /= sum_exp) then
        if     (exp(1) /= maxval(exp).and. nx > 1) then
          exp(1) = max(sum_exp - exp(2) - exp(3),0)
        elseif (exp(2) /= maxval(exp).and. ny > 1) then
          exp(2) = max(sum_exp - exp(1) - exp(3),0)
        elseif (exp(3) /= maxval(exp).and. nz > 1) then
          exp(3) = max(sum_exp - exp(1) - exp(2),0)
        endif
      endif

      if (proc_debug) write (*,*) 'DIAG -- procAlloc, exp_1=',exp

      !Consistency check (2nd layer)
      if (sum(exp) /= sum_exp) then
        if     (nx >= ny .and. nx >= nz) then
          exp(1) = max(sum_exp - exp(2) - exp(3),0)
        elseif (ny >= nx .and. ny >= nz) then
          exp(2) = max(sum_exp - exp(1) - exp(3),0)
        else
          exp(3) = max(sum_exp - exp(1) - exp(2),0)
        endif
      endif

c     Find processor distribution

      if (npx == 0) npx = 2**exp(1)
      if (npy == 0) npy = 2**exp(2)
      if (npz == 0) npz = 2**exp(3)

      if (proc_debug) then
        write (*,*) 'DIAG -- procAlloc, exp_2=',exp
        write (*,*) 'DIAG -- procAlloc, Procs=',npx,npy,npz
cc        stop
      endif

c     End program

      end subroutine processorAlloc
#endif

      end module grid_mpi

c module error
c ######################################################################
      module error

        use grid_mpi

      end module error

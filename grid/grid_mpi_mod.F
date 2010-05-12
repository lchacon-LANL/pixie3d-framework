c module grid_mpi
c ######################################################################
      module grid_mpi

        use grid_def

        implicit none

#if defined(petsc)

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdadef.h"

        integer    :: np,inp,my_rank,mpierr,group_world
     .               ,group_sp,tag=0,dest=0,root=0
     .               ,status(MPI_STATUS_SIZE)
     .               ,tag_send,tag_recv,request

        integer    :: l_lim(3,0:1),g_lim(3,0:1),rem_l_lim(3,0:1)
     .               ,rp_l_lim(3,0:1)

        integer    :: MPI_COMM_SP,MPI_COMM_POL,MPI_COMM_RAD
     .               ,MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z
     .               ,MPI_COMM_YZ
        integer   ,allocatable,dimension(:) :: MPI_COMM_PER
        integer   ,allocatable,dimension(:) :: MPI_COMM_NBRS

        character(80) :: messg

cc        private :: createMPIComm,selectRank_SP,selectRank_POL
cc     .            ,checkProc,processorAlloc
cc     $            ,fillLocalVec,emptyLocalVec
cc     $            ,createPETScFortranDA,destroyPETScFortranDA

      contains

c     initFortranMPI
c     #####################################################################
      subroutine initFortranMPI(mpi_comm,nxg,nyg,nzg)

c     ---------------------------------------------------------------------
c     Initializes MPI variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nxg,nyg,nzg,mpi_comm

c     Local variables

c     Begin program

        call initMPI(mpi_comm,np,my_rank)

c     Allocate processors

        call processorAlloc(nxg,nyg,nzg)

      end subroutine initFortranMPI

c     create_MPI_comms
c     #################################################################
      subroutine create_MPI_comms

c     -----------------------------------------------------------------
c     Creates required MPI communicators
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: my_sp_rank,np_sp

c     Begin program

c     Singular-point boundary MPI communicator

        if (bcSP()) call createMPIComm(MPI_COMM_SP,selectRank_SP)

cc        call check_MPI_COMM(MPI_COMM_SP,'MPI_COMM_SP')

c     X MPI communicator (domains that share the same Y, Z  slices)

        call createMPIComm(MPI_COMM_X,selectRank_X)

        if (bcSP()) MPI_COMM_RAD = MPI_COMM_X

cc        call check_MPI_COMM(MPI_COMM_X,'MPI_COMM_X')

c     Y MPI communicator (domains that share same X, Z slices)

        call createMPIComm(MPI_COMM_Y,selectRank_Y)

        if (bcSP()) MPI_COMM_POL = MPI_COMM_Y

cc        call check_MPI_COMM(MPI_COMM_Y,'MPI_COMM_Y')

c     YZ MPI communicator (domains that share same X slices)

        call createMPIComm(MPI_COMM_YZ,selectRank_YZ)

cc        call check_MPI_COMM(MPI_COMM_YZ,'MPI_COMM_YZ')

c     End program

      contains

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
        write (*,*) trim(mpi_comm_id)//'=',MPI_COMM
     .                ,'; nproc_local=',np_l
     .                ,'; my_rank_local=',my_rank_l
     .                ,'; my_rank_global=',my_rank
      endif

cc      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc      call MPI_Finalize(mpierr)
cc      stop

c     End program

      end subroutine check_MPI_COMM

      end subroutine create_MPI_comms

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

c     pstop
c     ################################################################
      subroutine pstop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        character(*)  :: routine, message

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)
        endif

        call PetscEnd(mpierr)

        stop

      end subroutine pstop

c     dot2
c     ###################################################################
      function dot2(vec1,vec2)

c     -------------------------------------------------------------------
c     Performs parallel scalar product (vec1,vec2). In call:
c       * vec1: vector, first term of scalar product.
c       * vec2: vector, second term of scalar product.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ntot
      real(8)    :: dot2
      real(8),INTENT(IN) :: vec1(:),vec2(:)

c     Local variables

      real(8)    :: ldot

c     Begin program

      dot2 = dot_product(vec1,vec2)

      if (.not.asm) then
        ldot = dot2
        call MPI_Allreduce(ldot,dot2,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      endif

      end function dot2

c     dot
c     ###################################################################
      function dot(ntot,vec1,vec2)

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

      integer    :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq
      real(8)    :: ldot

c     Begin program

      dot = dot_product(vec1,vec2)

      if (.not.asm) then
        ldot = dot
        call MPI_Allreduce(ldot,dot,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      endif

      end function dot

c     psum
c     ###################################################################
      function psum(ntot,vec1,mpi_comm)

c     -------------------------------------------------------------------
c     Performs parallel sum of element of an array using communicator
c     mpi_comm. On input:
c       * ntot: array dimension
c       * vec1: array
c       * mpi_comm: MPI communicator
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ntot
      integer,optional :: mpi_comm
      real(8)    :: vec1(ntot),psum

c     Local variables

      real(8)    :: lsum
      integer    :: mpic

c     Begin program

      if (PRESENT(mpi_comm)) then
        mpic = mpi_comm
      else
        mpic = MPI_COMM_WORLD
      endif

      psum = sum(vec1)

      if (.not.asm) then
        lsum = psum
        call MPI_Allreduce(lsum,psum,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,mpic,mpierr)
      endif

      end function psum

c     ######################## PRIVATE ROUTINES ###########################

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

ccc     createPETScFortranDA
ccc     #################################################################
cc      subroutine createPETScFortranDA(igrid,nxgl,nygl,nzgl)
cc
ccc     -----------------------------------------------------------------
ccc     Creates PETSc distributed array (DA) for boundary communication
ccc     in fortran at grid level igrid.
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: igrid
cc        integer   ,optional :: nxgl,nygl,nzgl
cc
ccc     Local variables
cc
cc        integer    :: BC,ierr,nxg,nyg,nzg,igx,igy,igz
cc
ccc     Begin program
cc
cc        igx = igrid
cc        igy = igrid
cc        igz = igrid
cc
cc        if (PRESENT(nxgl)) then
cc          nxg = nxgl
cc        else
cc          nxg = grid_params%nxgl(igx)
cc        endif
cc
cc        if (PRESENT(nygl)) then
cc          nyg = nygl
cc        else
cc          nyg = grid_params%nygl(igy)
cc        endif
cc
cc        if (PRESENT(nzgl)) then
cc          nzg = nzgl
cc        else
cc          nzg = grid_params%nzgl(igz)
cc        endif
cc
cc        if (bcond(1)==PER.and.bcond(3)==PER.and.bcond(5)==PER) then
cc          BC = DA_XYZPERIODIC
cc        elseif (bcond(1)==PER.and.bcond(3)==PER.and.bcond(5)/=PER) then
cc          BC = DA_XYPERIODIC
cc        elseif (bcond(1)==PER.and.bcond(3)/=PER.and.bcond(5)==PER) then
cc          BC = DA_XZPERIODIC
cc        elseif (bcond(1)/=PER.and.bcond(3)==PER.and.bcond(5)==PER) then
cc          BC = DA_YZPERIODIC
cc        elseif (bcond(1)==PER.and.bcond(3)/=PER.and.bcond(5)/=PER) then
cc          BC = DA_XPERIODIC
cc        elseif (bcond(1)/=PER.and.bcond(3)==PER.and.bcond(5)/=PER) then
cc          BC = DA_YPERIODIC
cc        elseif (bcond(1)/=PER.and.bcond(3)/=PER.and.bcond(5)==PER) then
cc          BC = DA_ZPERIODIC
cc        elseif (bcond(1)/=PER.and.bcond(3)/=PER.and.bcond(5)/=PER) then
cc          BC = DA_NONPERIODIC
cc        endif
cc
cc        if (npx == 0) npx = PETSC_DECIDE
cc        if (npy == 0) npy = PETSC_DECIDE
cc        if (npz == 0) npz = PETSC_DECIDE
cc
cc        call DACreate3d(PETSC_COMM_WORLD,BC,DA_STENCIL_BOX,nxg,nyg,nzg
cc     &                 ,npx,npy,npz,1,1
cc     &                 ,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER
cc     &                 ,PETSC_NULL_INTEGER,dactx(igrid)%da,ierr)
cc
cc        call DACreateGlobalVector(dactx(igrid)%da,dactx(igrid)%Xg,ierr)
cc
ccc       Get local grid boundaries in global grid
cc
cc        call DAGetCorners(dactx(igrid)%da
cc     .         ,dactx(igrid)%xs,dactx(igrid)%ys,dactx(igrid)%zs
cc     .         ,dactx(igrid)%xm,dactx(igrid)%ym,dactx(igrid)%zm,ierr)
cc        call DAGetGhostCorners(dactx(igrid)%da
cc     .         ,dactx(igrid)%gxs,dactx(igrid)%gys,dactx(igrid)%gzs
cc     .         ,dactx(igrid)%gxm,dactx(igrid)%gym,dactx(igrid)%gzm,ierr)
cc
cc        dactx(igrid)%xs  = dactx(igrid)%xs+1
cc        dactx(igrid)%ys  = dactx(igrid)%ys+1
cc        dactx(igrid)%zs  = dactx(igrid)%zs+1
cc        dactx(igrid)%gxs = dactx(igrid)%gxs+1
cc        dactx(igrid)%gys = dactx(igrid)%gys+1
cc        dactx(igrid)%gzs = dactx(igrid)%gzs+1
cc
cc        dactx(igrid)%ye  = dactx(igrid)%ys+dactx(igrid)%ym-1
cc        dactx(igrid)%xe  = dactx(igrid)%xs+dactx(igrid)%xm-1
cc        dactx(igrid)%ze  = dactx(igrid)%zs+dactx(igrid)%zm-1
cc        dactx(igrid)%gye = dactx(igrid)%gys+dactx(igrid)%gym-1
cc        dactx(igrid)%gxe = dactx(igrid)%gxs+dactx(igrid)%gxm-1
cc        dactx(igrid)%gze = dactx(igrid)%gzs+dactx(igrid)%gzm-1
cc
cc        !With ghost cells (only those that PETSc includes)
cc        dactx(igrid)%lgxs = dactx(igrid)%gxs - dactx(igrid)%xs + 1
cc        dactx(igrid)%lgys = dactx(igrid)%gys - dactx(igrid)%ys + 1
cc        dactx(igrid)%lgzs = dactx(igrid)%gzs - dactx(igrid)%zs + 1
cc        dactx(igrid)%lgxe = dactx(igrid)%gxe - dactx(igrid)%xs + 1
cc        dactx(igrid)%lgye = dactx(igrid)%gye - dactx(igrid)%ys + 1
cc        dactx(igrid)%lgze = dactx(igrid)%gze - dactx(igrid)%zs + 1
cc
cc        !Domain limits (without ghost cells)
cc        dactx(igrid)%lxs = 1
cc        dactx(igrid)%lys = 1
cc        dactx(igrid)%lzs = 1
cc        dactx(igrid)%lxe = dactx(igrid)%xm
cc        dactx(igrid)%lye = dactx(igrid)%ym
cc        dactx(igrid)%lze = dactx(igrid)%zm
cc
cc        !Array limits (with all ghost cells)
cc        dactx(igrid)%il = dactx(igrid)%lxs
cc        dactx(igrid)%jl = dactx(igrid)%lys
cc        dactx(igrid)%kl = dactx(igrid)%lzs
cc        dactx(igrid)%ih = dactx(igrid)%lxe
cc        dactx(igrid)%jh = dactx(igrid)%lye
cc        dactx(igrid)%kh = dactx(igrid)%lze
cc
cc        dactx(igrid)%ilm = dactx(igrid)%il-1
cc        dactx(igrid)%ihp = dactx(igrid)%ih+1
cc        dactx(igrid)%jlm = dactx(igrid)%jl-1
cc        dactx(igrid)%jhp = dactx(igrid)%jh+1
cc        dactx(igrid)%klm = dactx(igrid)%kl-1
cc        dactx(igrid)%khp = dactx(igrid)%kh+1
cc
ccc     Store grid level
cc
cc        dactx(igrid)%igx = igx
cc        dactx(igrid)%igy = igy
cc        dactx(igrid)%igz = igz
cc
ccc     End program
cc
cc      end subroutine createPETScFortranDA
cc
ccc     destroyPETScFortranDA
ccc     #################################################################
cc      subroutine destroyPETScFortranDA(dactx)
cc
ccc     -----------------------------------------------------------------
ccc     Deallocates PETSc DA
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        type(petsc_da_ctx),pointer :: dactx
cc
ccc     Local variables
cc
cc        integer    :: ierr,igr
cc
ccc     Begin program
cc
cccc        do igr=1,ngrid
cccc          call VecDestroy(dactx(igr)%Xg,ierr)
cccc          call DADestroy (dactx(igr)%da,ierr)
cccc        enddo
cc
cc        call VecDestroy(dactx%Xg,ierr)
cc        call DADestroy (dactx%da,ierr)
cc
ccc     End program
cc
cc      end subroutine destroyPETScFortranDA

c     inProc
c     #################################################################
      function inProc(igl,jgl,kgl,igx,igy,igz)

        implicit none

c     Call variables

        integer    :: igl,jgl,kgl,igx,igy,igz
        logical    :: inProc

c     Local variables

c     Begin program

        inProc =igl>=grid_params%ilo(igx).and.igl<=grid_params%ihi(igx)
     .     .and.jgl>=grid_params%jlo(igy).and.jgl<=grid_params%jhi(igy)
     .     .and.kgl>=grid_params%klo(igz).and.kgl<=grid_params%khi(igz)

c     End program

      end function inProc

c     processorAlloc
c     ######################################################################
      subroutine processorAlloc(nx1,ny1,nz1)

c     ----------------------------------------------------------------------
c     Finds processor allocation (npx,npy,npz; defined in module) in the
c     three logical directions based on global dimensions
c     (nx1,nx2,nx3). Only returns processor allocation for directions
c     not specified earlier by user (i.e., when np* is zero).
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: nx1,ny1,nz1

c     Local variables

      integer    :: nx,ny,nz,sum_exp,nd,navg,exp(3),npt(3),nprocs

c     Begin program

      nx = nx1
      ny = ny1
      nz = nz1

      nprocs = np

      exp = 0

c     Eliminate specified directions

      if (npx /= 0) then
        nprocs = nprocs/npx
        nx = 1
cc        asm_dir(1) = (npx > 1)
      endif

      if (npy /= 0) then
        nprocs = nprocs/npy
        ny = 1
cc        asm_dir(2) = (npy > 1)
      endif

      if (npz /= 0) then
        nprocs = nprocs/npz
        nz = 1
cc        asm_dir(3) = (npz > 1)
      endif

      if (nx == 1 .and. ny == 1 .and. nz == 1) then
        if (nprocs == 1) then  !Compatible specification of procs. per direction
          npx = max(npx,1) ; npy = max(npy,1) ; npz = max(npz,1)
          return
        else
          messg = 'np /= npx*npy*npz ... Aborting'
          call pstop('processorAlloc',messg)
        endif
      endif

c     Check MG is an option

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

c     Exclude directions based on topological constraints (PER, SP)

      !Priority: SP, PER
cc      if (bcond(1) == SP  .or.  bcond(3) == PER ) ny = 1
cc      if (bcond(1) == PER .and. ny /= 1         ) nx = 1
cc      if (bcond(5) == PER .and.(nx/=1.and.ny/=1)) nz = 1

      !Priority: SP
cc      if (bcond(1) == SP) ny = 1

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
cc      navg = nint((nx*ny*nz/nprocs)**(1./nd))
      navg = int((nx*ny*nz/nprocs)**(1./nd))

      npt = (/nx,ny,nz/)
cc      exp = nint(log(max(1d0*npt/navg,1d0))/log(2d0) + 0.3)
      exp = int(log(max(1d0*npt/navg,1d0))/log(2d0) + 0.3)

cc      write (*,*) 'DIAG -- procAlloc',nx,ny,nz,nprocs,navg,exp,sum_exp

      !Consistency check
      if (sum(exp) /= sum_exp) then
        if (nx >= ny .and. nx >= nz) then
          exp(1) = sum_exp - exp(2) - exp(3)
        elseif (ny >= nx .and. ny >= nz) then
          exp(2) = sum_exp - exp(1) - exp(3)
        else
          exp(3) = sum_exp - exp(1) - exp(2)
        endif
      endif

c     Find processor distribution

      if (npx == 0) npx = 2**exp(1)
      if (npy == 0) npy = 2**exp(2)
      if (npz == 0) npz = 2**exp(3)

c     Store processor topology for ASM PC treatment

cc      asm_dir(1) = (npx > 1)
cc      asm_dir(2) = (npy > 1)
cc      asm_dir(3) = (npz > 1)

c     End program

      end subroutine processorAlloc

c     createMPIComm
c     #####################################################################
      subroutine createMPIComm(MPI_COMM,selProc)

c     ---------------------------------------------------------------------
c     Creates MPI communicator MPI_COMM according to routine selProc.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: MPI_COMM

c     Local variables

        integer    :: split_key,my_sp_rank,np_sp

        external   selProc

c     Begin program

c     Find split key

        split_key = MPI_UNDEFINED

        call selProc(split_key)

c     Create communicator

        call MPI_Comm_split(MPI_COMM_WORLD,split_key,my_rank
     .                     ,MPI_COMM,mpierr)

c     End program

      end subroutine createMPIComm

c     selectRank_SP
c     #####################################################################
      subroutine selectRank_SP(split_key)

c     ---------------------------------------------------------------------
c     Groups processors surrounding singular point
c     (SP boundary condition)
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: split_key

c     Local variables

        real(8)    :: random
        integer    :: dim,loc,lrank,BC

        logical    :: cp(0:1)

c     Begin program

        BC  = SP
        dim = 1

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = grid_params%ilo(1)
        l_lim(2,0) = grid_params%jlo(1)
        l_lim(3,0) = grid_params%klo(1)
        l_lim(1,1) = grid_params%ihi(1)
        l_lim(2,1) = grid_params%jhi(1)
        l_lim(3,1) = grid_params%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = grid_params%nxgl(1)
        g_lim(2,1) = grid_params%nygl(1)
        g_lim(3,1) = grid_params%nzgl(1)

c     Select rank

        do loc=0,1
          cp(loc) = checkProc(BC,dim,loc,l_lim,g_lim)
          if (cp(loc)) split_key = BC
        enddo

        !Group communicators based on location in z direction
        if (bcond(dim) == SP .and. cp(0)) then
          split_key = split_key + npz*l_lim(3,1)/g_lim(3,1)
        endif

        !Create self-communicator for single-point dimensions
        if (      (l_lim(dim,0) == l_lim(dim,1))
     .      .and. (cp(0)       .eqv. cp(1)     )  ) then
          call random_seed()
          call random_number(random)
          split_key=split_key+my_rank+1
        endif

c     End program

      end subroutine selectRank_SP

ccc     selectRank_SP
ccc     #####################################################################
cc      subroutine selectRank_SP(split_key)
cc
ccc     ---------------------------------------------------------------------
ccc     Groups processors surrounding singular point
ccc     (SP boundary condition)
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: split_key
cc
ccc     Local variables
cc
ccc     Begin program
cc
ccc     Set global limits of local domain l_lim(dim,loc)
cc
cc        l_lim(1,0) = grid_params%ilo(1)
cc        l_lim(2,0) = grid_params%jlo(1)
cc        l_lim(3,0) = grid_params%klo(1)
cc        l_lim(1,1) = grid_params%ihi(1)
cc        l_lim(2,1) = grid_params%jhi(1)
cc        l_lim(3,1) = grid_params%khi(1)
cc
ccc     Set global limits of global domain g_lim(dim,loc)
cc
cc        g_lim(1,0) = 1
cc        g_lim(2,0) = 1
cc        g_lim(3,0) = 1
cc        g_lim(1,1) = grid_params%nxgl(1)
cc        g_lim(2,1) = grid_params%nygl(1)
cc        g_lim(3,1) = grid_params%nzgl(1)
cc
ccc     Select rank: group communicators based on:
ccc        1) SP boundary
ccc        2) location in z direction
cc
cc        if (isSP2(1,1,dim=1,loc=0)) then
cc          split_key = 1
cc          split_key = split_key + npz*l_lim(3,1)/g_lim(3,1)
cc        endif
cc
ccc     End program
cc
cc      end subroutine selectRank_SP

c     selectRank_X
c     #####################################################################
      subroutine selectRank_X(split_key)

c     ---------------------------------------------------------------------
c     Groups processors in slices in X direction (i.e., order them according
c     to Y and Z positions).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: split_key

c     Local variables

        real(8)    :: random
        integer    :: loc,lrank

        logical    :: cp(0:1)

c     Begin program

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = grid_params%ilo(1)
        l_lim(2,0) = grid_params%jlo(1)
        l_lim(3,0) = grid_params%klo(1)
        l_lim(1,1) = grid_params%ihi(1)
        l_lim(2,1) = grid_params%jhi(1)
        l_lim(3,1) = grid_params%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = grid_params%nxgl(1)
        g_lim(2,1) = grid_params%nygl(1)
        g_lim(3,1) = grid_params%nzgl(1)

c     Select rank

        split_key = 1

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

        integer    :: split_key

c     Local variables

        real(8)    :: random
        integer    :: loc,lrank

        logical    :: cp(0:1)

c     Begin program

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = grid_params%ilo(1)
        l_lim(2,0) = grid_params%jlo(1)
        l_lim(3,0) = grid_params%klo(1)
        l_lim(1,1) = grid_params%ihi(1)
        l_lim(2,1) = grid_params%jhi(1)
        l_lim(3,1) = grid_params%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = grid_params%nxgl(1)
        g_lim(2,1) = grid_params%nygl(1)
        g_lim(3,1) = grid_params%nzgl(1)

c     Select rank

        split_key = 1

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

        integer    :: split_key

c     Local variables

        real(8)    :: random
        integer    :: loc,lrank

        logical    :: cp(0:1)

c     Begin program

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = grid_params%ilo(1)
        l_lim(2,0) = grid_params%jlo(1)
        l_lim(3,0) = grid_params%klo(1)
        l_lim(1,1) = grid_params%ihi(1)
        l_lim(2,1) = grid_params%jhi(1)
        l_lim(3,1) = grid_params%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = grid_params%nxgl(1)
        g_lim(2,1) = grid_params%nygl(1)
        g_lim(3,1) = grid_params%nzgl(1)

c     Select rank

        split_key = 1

        !Group processors based on location in x and y directions
        !(z direction is shared; 
        ! organize in 2D processor mesh: j + ny*(i-1))
        split_key = split_key +     npy*l_lim(2,1)/g_lim(2,1)
     .                        + npy*npx*l_lim(1,1)/g_lim(1,1)

c     End program

      end subroutine selectRank_Z

c     selectRank_YZ
c     #####################################################################
      subroutine selectRank_YZ(split_key)

c     ---------------------------------------------------------------------
c     Groups processors based on shared X position in processor mesh.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: split_key

c     Local variables

        real(8)    :: random
        integer    :: loc,lrank

        logical    :: cp(0:1)

c     Begin program

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = grid_params%ilo(1)
        l_lim(2,0) = grid_params%jlo(1)
        l_lim(3,0) = grid_params%klo(1)
        l_lim(1,1) = grid_params%ihi(1)
        l_lim(2,1) = grid_params%jhi(1)
        l_lim(3,1) = grid_params%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = grid_params%nxgl(1)
        g_lim(2,1) = grid_params%nygl(1)
        g_lim(3,1) = grid_params%nzgl(1)

c     Select rank

        split_key = 1

        !Group processors based on location in x direction
        !(y and z directions are shared; 
        ! organize in 1D processor mesh: i
        split_key = split_key + npx*l_lim(1,1)/g_lim(1,1)

c     End program

      end subroutine selectRank_YZ

ccc     createSPComm
ccc     #####################################################################
cc      subroutine createSPComm(dim,MPI_COMM)
cc
ccc     ---------------------------------------------------------------------
ccc     Creates MPI communicator for domains surrounding singular point
ccc     (SP boundary condition)
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: dim,MPI_COMM
cc
ccc     Local variables
cc
cc        integer    :: split_key,my_sp_rank,np_sp
cc
ccc     Begin program
cc
ccc     Find split key
cc
cc        split_key = MPI_UNDEFINED
cc
cc        call selectRank(dim,SP,split_key)
cc
ccc     Create communicator
cc
cc        call MPI_Comm_split(MPI_COMM_WORLD,split_key,my_rank
cc     .                     ,MPI_COMM,mpierr)
cc
ccc diag ****
cccccc        write (*,*) 'Proc',my_rank,', SP COMM' ,MPI_COMM_SP
cccc        if (MPI_COMM_SP /= MPI_COMM_NULL) then
cccc          call MPI_Comm_size(MPI_COMM_SP,np_sp,mpierr)
cccc          call MPI_Comm_rank(MPI_COMM_SP,my_sp_rank,mpierr)
cccc          write (*,*) 'SP_COMM ',split_key,'; nproc=',np_sp
cccc     .                ,'; my_sp_rank=',my_sp_rank
cccc     .                ,'; my_rank=',my_rank
cccc        endif
cccc        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cccc        call MPI_Finalize(mpierr)
cccc        stop
ccc diag ****
cc
ccc     End program
cc
cc        contains
cc
ccc       selectRank
ccc       #####################################################################
cc        subroutine selectRank(dim,BC,split_key)
cc
ccc       ---------------------------------------------------------------------
ccc       Groups processors within groups based on split_key.
ccc       ---------------------------------------------------------------------
cc
cc          implicit none
cc
ccc       Call variables
cc
cc          integer    :: dim,BC,split_key
cc
ccc       Local variables
cc
cc          real(8)    :: random
cc          integer    :: loc,lrank
cc
cc          logical    :: cp(0:1)
cc
ccc       Begin program
cc
ccc       Set global limits of local domain l_lim(dim,loc)
cc
cc          l_lim(1,0) = grid_params%ilo(1)
cc          l_lim(2,0) = grid_params%jlo(1)
cc          l_lim(3,0) = grid_params%klo(1)
cc          l_lim(1,1) = grid_params%ihi(1)
cc          l_lim(2,1) = grid_params%jhi(1)
cc          l_lim(3,1) = grid_params%khi(1)
cc
ccc       Set global limits of global domain g_lim(dim,loc)
cc
cc          g_lim(1,0) = 1
cc          g_lim(2,0) = 1
cc          g_lim(3,0) = 1
cc          g_lim(1,1) = grid_params%nxgl(1)
cc          g_lim(2,1) = grid_params%nygl(1)
cc          g_lim(3,1) = grid_params%nzgl(1)
cc
ccc       Select rank
cc
cc          do loc=0,1
cc            cp(loc) = checkProc(BC,dim,loc,l_lim,g_lim)
cc            if (cp(loc)) split_key = BC
cc          enddo
cc
cc          !Group communicators based on location in z direction
cc          if (bcond(dim) == SP .and. cp(0)) then
cc            split_key = split_key + npz*l_lim(3,1)/g_lim(3,1)
cc          endif
cc
cc          !Create self-communicator for single-point dimensions
cc          if (      l_lim(dim,0) == l_lim(dim,1)
cc     .        .and. cp(0)        == cp(1)       ) then
cc            call random_seed()
cc            call random_number(random)
cc            split_key=split_key+my_rank+1
cc          endif
cc
ccc       End program
cc
cc        end subroutine selectRank
cc
cc      end subroutine createSPComm

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

ccc     setupPETScLocalLimits
ccc     #################################################################
cc      subroutine setupPETScLocalLimits(igrid)
cc
ccc     -----------------------------------------------------------------
ccc     Defines local limits of current domain, and stores them in dactx.
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: igrid
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        !With ghost cells (only those that PETSc includes)
cc        call fromGlobalToLocalLimits
cc     .       (dactx(igrid)%gxs ,dactx(igrid)%gys ,dactx(igrid)%gzs
cc     $       ,dactx(igrid)%lgxs,dactx(igrid)%lgys,dactx(igrid)%lgzs
cc     $       ,dactx(igrid)%igx ,dactx(igrid)%igy ,dactx(igrid)%igz)
cc
cc        call fromGlobalToLocalLimits
cc     .       (dactx(igrid)%gxe ,dactx(igrid)%gye ,dactx(igrid)%gze
cc     $       ,dactx(igrid)%lgxe,dactx(igrid)%lgye,dactx(igrid)%lgze
cc     $       ,dactx(igrid)%igx ,dactx(igrid)%igy ,dactx(igrid)%igz)
cc
cc        !Domain limits (without ghost cells)
cc        call fromGlobalToLocalLimits
cc     .       (dactx(igrid)%xs,dactx(igrid)%ys,dactx(igrid)%zs
cc     $       ,dactx(igrid)%il,dactx(igrid)%jl,dactx(igrid)%kl
cc     $       ,dactx(igrid)%igx,dactx(igrid)%igy,dactx(igrid)%igz)
cc        call fromGlobalToLocalLimits
cc     .       (dactx(igrid)%xe,dactx(igrid)%ye,dactx(igrid)%ze
cc     $       ,dactx(igrid)%ih,dactx(igrid)%jh,dactx(igrid)%kh
cc     $       ,dactx(igrid)%igx,dactx(igrid)%igy,dactx(igrid)%igz)
cc
cc        !Array limits (with all ghost cells)
ccc$$$        dactx(igrid)%il = 1
ccc$$$        dactx(igrid)%jl = 1
ccc$$$        dactx(igrid)%kl = 1
ccc$$$        dactx(igrid)%ih = grid_params%nxv(igx)
ccc$$$        dactx(igrid)%jh = grid_params%nyv(igy)
ccc$$$        dactx(igrid)%kh = grid_params%nzv(igz)
cc
cc        dactx(igrid)%ilm = dactx(igrid)%il-1
cc        dactx(igrid)%ihp = dactx(igrid)%ih+1
cc        dactx(igrid)%jlm = dactx(igrid)%jl-1
cc        dactx(igrid)%jhp = dactx(igrid)%jh+1
cc        dactx(igrid)%klm = dactx(igrid)%kl-1
cc        dactx(igrid)%khp = dactx(igrid)%kh+1
cc
ccc     End program
cc
cc      end subroutine setupPETScLocalLimits

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

c     isProc
c     #################################################################
      function isProc(proc_id)

        implicit none

c     Call variables

        integer    :: proc_id
        logical    :: isProc

c     Local variables

c     Begin program

      isProc = .true.

c     End program

      end function isProc

c     inProc
c     #################################################################
      function inProc(igl,jgl,kgl,igx,igy,igz)

        implicit none

c     Call variables

        integer    :: igl,jgl,kgl,igx,igy,igz
        logical    :: inProc

c     Local variables

c     Begin program

        inProc = .true.

c     End program

      end function inProc

c     pstop
c     ################################################################
      subroutine pstop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        character(*)  :: routine, message

c     Begin program

        write (*,*)
        write (*,*) trim(message)
        write (*,*) 'Program stopped at routine ',trim(routine)

        stop

      end subroutine pstop

c     dot2
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
      real(8), INTENT(IN) :: vec1(:),vec2(:)

c     Local variables

c     Begin program

      dot2 = dot_product(vec1,vec2)

      end function dot2

c     dot
c     ###################################################################
      function dot(ntot,vec1,vec2)

c     -------------------------------------------------------------------
c     Performs scalar product (vec1,vec2). In call:
c       * ntot: vector dimension
c       * vec1: vector, first term of scalar product.
c       * vec2: vector, second term of scalar product.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ntot
      real(8)    :: vec1(ntot),vec2(ntot),dot

c     Local variables

c     Begin program

      dot = dot_product(vec1,vec2)

      end function dot

#endif

      end module grid_mpi

c module error
c ######################################################################
      module error

        use grid_mpi

      end module error

c module grid_mpi
c ######################################################################
      module grid_mpi

        use grid_def

#if defined(petsc)

        implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscda.h"
#include "finclude/petscvec.h90"

        integer    :: np,inp,my_rank,mpierr,group_world
     .               ,group_sp,tag=0,dest=0,root=0
     .               ,status(MPI_STATUS_SIZE)
     .               ,tag_send,tag_recv,request
     .               ,npx=0,npy=0,npz=0

        integer    :: l_lim(3,0:1),g_lim(3,0:1),rem_l_lim(3,0:1)
     .               ,rp_l_lim(3,0:1)

        integer    :: MPI_COMM_SP,MPI_COMM_POL
        integer   ,allocatable,dimension(:) :: MPI_COMM_PER
        integer   ,allocatable,dimension(:) :: MPI_COMM_NBRS

        type :: da_ctx
          DA         :: da
          Vec        :: Xg
          integer    :: xs,xe,xm,gxs,gxe,gxm
          integer    :: ys,ye,ym,gys,gye,gym
          integer    :: zs,ze,zm,gzs,gze,gzm
          integer    :: lxs,lxe,lgxs,lgxe
          integer    :: lys,lye,lgys,lgye
          integer    :: lzs,lze,lgzs,lgze
          integer    :: il ,ih ,jl ,jh ,kl ,kh
          integer    :: ilm,ihp,jlm,jhp,klm,khp
          integer    :: mx,my,mz,rank
          integer    :: igx,igy,igz
        end type da_ctx

        type (da_ctx),dimension(20) :: dactx

        character(80) :: messg

        private :: createMPIComm,selectRank_SP,selectRank_POL
     .            ,checkProc,processorAlloc
     $            ,fillLocalVec,emptyLocalVec
     $            ,createPETScFortranDA,destroyPETScFortranDA

      contains

c     initFortranMPI
c     #####################################################################
      subroutine initFortranMPI(nxg,nyg,nzg)

c     ---------------------------------------------------------------------
c     Initializes MPI variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nxg,nyg,nzg

c     Local variables

c     Begin program

#if !defined(petsc_c)
cc        call PetscInitialize(PETSC_NULL_CHARACTER,mpierr)
#endif
        call initMPI(nxg,nyg,nzg)

        call createPETScFortranDA(1,nxgl=nxg,nygl=nyg,nzgl=nzg)

      end subroutine initFortranMPI

c     destroyFortranMPI
c     #####################################################################
      subroutine destroyFortranMPI

c     ---------------------------------------------------------------------
c     Initializes MPI variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        call destroyPETScFortranDA

      end subroutine destroyFortranMPI

c     createPETScGrid
c     #################################################################
      subroutine createPETScGrid

c     -----------------------------------------------------------------
c     Creates PETSc constructs for boundary communication in fortran.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: igr

        integer    :: my_sp_rank,np_sp

c     Begin program

c     Create PETSc DA MG hierachy

#if defined(petsc_c)
        call createPETScFortranDA(1)
#endif

        do igr=2,grid_params%ngrid
          call createPETScFortranDA(igr)
        enddo

c     Singular-point boundary MPI communicator

cc        if (bcSP()) call createSPComm(1,MPI_COMM_SP)
        if (bcSP()) call createMPIComm(MPI_COMM_SP,selectRank_SP)

c diag ****
cc        if (MPI_COMM_SP /= MPI_COMM_NULL) then
cc          call MPI_Comm_size(MPI_COMM_SP,np_sp,mpierr)
cc          call MPI_Comm_rank(MPI_COMM_SP,my_sp_rank,mpierr)
cc          write (*,*) 'SP_COMM ; nproc=',np_sp
cc     .                ,'; my_sp_rank=',my_sp_rank
cc     .                ,'; my_rank=',my_rank
cc        endif
cc        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc        call MPI_Finalize(mpierr)
cc        stop
c diag ****

c     Poloidal MPI communicator (domains that share same toroidal and radial slices)

#if defined(vmec)
        call createMPIComm(MPI_COMM_POL,selectRank_POL)

c diag ****
c$$$        if (MPI_COMM_POL /= MPI_COMM_NULL) then
c$$$          call MPI_Comm_size(MPI_COMM_POL,np_sp,mpierr)
c$$$          call MPI_Comm_rank(MPI_COMM_POL,my_sp_rank,mpierr)
c$$$          write (*,*) 'POL_COMM ; nproc=',np_sp
c$$$     .                ,'; my_rank_pol=',my_sp_rank
c$$$     .                ,'; my_rank=',my_rank
c$$$        endif
c$$$        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
c$$$        call MPI_Finalize(mpierr)
c$$$        stop
c diag ****
#endif

c     End program

      end subroutine createPETScGrid

c     fillPetscGhostCells
c     #################################################################
      subroutine fillPetscGhostCells(array,igr)

c     -----------------------------------------------------------------
c     Performs communication to fill ghost cells in array using PETSc
c     distributed arrays.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igr
        real(8) :: array(dactx(igr)%ilm:dactx(igr)%ihp
     .                  ,dactx(igr)%jlm:dactx(igr)%jhp 
     .                  ,dactx(igr)%klm:dactx(igr)%khp)

c     Local variables

        Vec :: localX

        PetscScalar,pointer :: lx_v(:)

        integer    :: ierr

c     Begin program

c     Fill known values of dactx%Xg

        !Get pointer to vector data
        call DAGetLocalVector(dactx(igr)%da,localX,ierr)
cc        call fchkerrq(ierr,'DAGetLocalVector')

        call VecGetArrayF90  (localX,lx_v,ierr)

        !Assignment: array -> lx_v
        call fillLocalVec(array,lx_v,igr)

        !Restore vector
        call VecRestoreArrayF90(localX,lx_v,ierr)

        !Insert values into global vector
        call DALocalToGlobal(dactx(igr)%da,localX
     .                      ,INSERT_VALUES,dactx(igr)%Xg
     $                      ,ierr)
cc        call DARestoreLocalVector(dactx(igr)%da,localX,ierr)


c     Fill ghost cells

        !Get pointer to vector data w/ghost cells
cc        call DAGetLocalVector    (dactx(igr)%da,localX,ierr)
        call DAGlobalToLocalBegin(dactx(igr)%da,dactx(igr)%Xg
     .                           ,INSERT_VALUES,localX
     $                           ,ierr)

        call DAGlobalToLocalEnd  (dactx(igr)%da,dactx(igr)%Xg
     .                           ,INSERT_VALUES,localX
     $                           ,ierr)


        call VecGetArrayF90(localX,lx_v,ierr)

        !Assignment: lx_v -> array
        call emptyLocalVec(array,lx_v,igr)

        !Restore vector
        call VecRestoreArrayF90(localX,lx_v,ierr)

c     Deallocate memory

        call DARestoreLocalVector(dactx(igr)%da,localX,ierr)

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

c     ######################## PRIVATE ROUTINES ###########################

c     initMPI
c     #####################################################################
      subroutine initMPI(nxg,nyg,nzg)

c     ---------------------------------------------------------------------
c     Initializes MPI variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nxg,nyg,nzg

c     Local variables

c     Begin program

c     Rank and number of processors

        call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
        call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)

c     Allocate processors

        call processorAlloc(nxg,nyg,nzg)

      end subroutine initMPI

c     createPETScFortranDA
c     #################################################################
      subroutine createPETScFortranDA(igrid,nxgl,nygl,nzgl)

c     -----------------------------------------------------------------
c     Creates PETSc distributed array (DA) for boundary communication
c     in fortran at grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igrid
        integer   ,optional :: nxgl,nygl,nzgl

c     Local variables

        integer    :: BC,ierr,nxg,nyg,nzg,igx,igy,igz

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        if (PRESENT(nxgl)) then
          nxg = nxgl
        else
          nxg = grid_params%nxgl(igx)
        endif

        if (PRESENT(nygl)) then
          nyg = nygl
        else
          nyg = grid_params%nygl(igy)
        endif

        if (PRESENT(nzgl)) then
          nzg = nzgl
        else
          nzg = grid_params%nzgl(igz)
        endif

        if (bcond(1)==PER.and.bcond(3)==PER.and.bcond(5)==PER) then
          BC = DA_XYZPERIODIC
        elseif (bcond(1)==PER.and.bcond(3)==PER.and.bcond(5)/=PER) then
          BC = DA_XYPERIODIC
        elseif (bcond(1)==PER.and.bcond(3)/=PER.and.bcond(5)==PER) then
          BC = DA_XZPERIODIC
        elseif (bcond(1)/=PER.and.bcond(3)==PER.and.bcond(5)==PER) then
          BC = DA_YZPERIODIC
        elseif (bcond(1)==PER.and.bcond(3)/=PER.and.bcond(5)/=PER) then
          BC = DA_XPERIODIC
        elseif (bcond(1)/=PER.and.bcond(3)==PER.and.bcond(5)/=PER) then
          BC = DA_YPERIODIC
        elseif (bcond(1)/=PER.and.bcond(3)/=PER.and.bcond(5)==PER) then
          BC = DA_ZPERIODIC
        elseif (bcond(1)/=PER.and.bcond(3)/=PER.and.bcond(5)/=PER) then
          BC = DA_NONPERIODIC
        endif

        if (npx == 0) npx = PETSC_DECIDE
        if (npy == 0) npy = PETSC_DECIDE
        if (npz == 0) npz = PETSC_DECIDE

        call DACreate3d(PETSC_COMM_WORLD,BC,DA_STENCIL_BOX,nxg,nyg,nzg
     &                 ,npx,npy,npz,1,1
     &                 ,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER
     &                 ,PETSC_NULL_INTEGER,dactx(igrid)%da,ierr)

        call DACreateGlobalVector(dactx(igrid)%da,dactx(igrid)%Xg,ierr)

c       Get local grid boundaries in global grid

        call DAGetCorners(dactx(igrid)%da
     .         ,dactx(igrid)%xs,dactx(igrid)%ys,dactx(igrid)%zs
     .         ,dactx(igrid)%xm,dactx(igrid)%ym,dactx(igrid)%zm,ierr)
        call DAGetGhostCorners(dactx(igrid)%da
     .         ,dactx(igrid)%gxs,dactx(igrid)%gys,dactx(igrid)%gzs
     .         ,dactx(igrid)%gxm,dactx(igrid)%gym,dactx(igrid)%gzm,ierr)

        dactx(igrid)%xs  = dactx(igrid)%xs+1
        dactx(igrid)%ys  = dactx(igrid)%ys+1
        dactx(igrid)%zs  = dactx(igrid)%zs+1
        dactx(igrid)%gxs = dactx(igrid)%gxs+1
        dactx(igrid)%gys = dactx(igrid)%gys+1
        dactx(igrid)%gzs = dactx(igrid)%gzs+1

        dactx(igrid)%ye  = dactx(igrid)%ys+dactx(igrid)%ym-1
        dactx(igrid)%xe  = dactx(igrid)%xs+dactx(igrid)%xm-1
        dactx(igrid)%ze  = dactx(igrid)%zs+dactx(igrid)%zm-1
        dactx(igrid)%gye = dactx(igrid)%gys+dactx(igrid)%gym-1
        dactx(igrid)%gxe = dactx(igrid)%gxs+dactx(igrid)%gxm-1
        dactx(igrid)%gze = dactx(igrid)%gzs+dactx(igrid)%gzm-1

        !With ghost cells (only those that PETSc includes)
        dactx(igrid)%lgxs = dactx(igrid)%gxs - dactx(igrid)%xs + 1
        dactx(igrid)%lgys = dactx(igrid)%gys - dactx(igrid)%ys + 1
        dactx(igrid)%lgzs = dactx(igrid)%gzs - dactx(igrid)%zs + 1
        dactx(igrid)%lgxe = dactx(igrid)%gxe - dactx(igrid)%xs + 1
        dactx(igrid)%lgye = dactx(igrid)%gye - dactx(igrid)%ys + 1
        dactx(igrid)%lgze = dactx(igrid)%gze - dactx(igrid)%zs + 1

        !Domain limits (without ghost cells)
        dactx(igrid)%lxs = 1
        dactx(igrid)%lys = 1
        dactx(igrid)%lzs = 1
        dactx(igrid)%lxe = dactx(igrid)%xm
        dactx(igrid)%lye = dactx(igrid)%ym
        dactx(igrid)%lze = dactx(igrid)%zm

        !Array limits (with all ghost cells)
        dactx(igrid)%il = dactx(igrid)%lxs
        dactx(igrid)%jl = dactx(igrid)%lys
        dactx(igrid)%kl = dactx(igrid)%lzs
        dactx(igrid)%ih = dactx(igrid)%lxe
        dactx(igrid)%jh = dactx(igrid)%lye
        dactx(igrid)%kh = dactx(igrid)%lze

        dactx(igrid)%ilm = dactx(igrid)%il-1
        dactx(igrid)%ihp = dactx(igrid)%ih+1
        dactx(igrid)%jlm = dactx(igrid)%jl-1
        dactx(igrid)%jhp = dactx(igrid)%jh+1
        dactx(igrid)%klm = dactx(igrid)%kl-1
        dactx(igrid)%khp = dactx(igrid)%kh+1

c     Store grid level

        dactx(igrid)%igx = igx
        dactx(igrid)%igy = igy
        dactx(igrid)%igz = igz

c     End program

      end subroutine createPETScFortranDA

c     destroyPETScFortranDA
c     #################################################################
      subroutine destroyPETScFortranDA

c     -----------------------------------------------------------------
c     Deallocates PETSc DA
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer    :: ierr,igr

c     Begin program

        do igr=1,grid_params%ngrid
          call VecDestroy(dactx(igr)%Xg,ierr)
          call DADestroy (dactx(igr)%da,ierr)
        enddo

c     End program

      end subroutine destroyPETScFortranDA

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
        asm_dir(1) = (npx > 1)
      endif

      if (npy /= 0) then
        nprocs = nprocs/npy
        ny = 1
        asm_dir(2) = (npy > 1)
      endif

      if (npz /= 0) then
        nprocs = nprocs/npz
        nz = 1
        asm_dir(3) = (npz > 1)
      endif

      if (nx == 1 .and. ny == 1 .and. nz == 1) then
        if (nprocs == 1) then  !Compatible specification of procs. per direction
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
      navg = nint((nx*ny*nz/nprocs)**(1./nd))

      npt = (/nx,ny,nz/)
      exp = nint(log(max(1d0*npt/navg,1d0))/log(2d0))

cc      write (*,*) 'DIAG -- procAlloc',nx,ny,nz,nprocs,navg,exp

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

      asm_dir(1) = (npx > 1)
      asm_dir(2) = (npy > 1)
      asm_dir(3) = (npz > 1)

cc      write (*,*) 'DIAG -- processorAlloc',npx,npy,npz
cc      call pstop('','')

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
        if (      l_lim(dim,0) == l_lim(dim,1)
     .      .and. cp(0)       .eqv. cp(1)       ) then
          call random_seed()
          call random_number(random)
          split_key=split_key+my_rank+1
        endif

c     End program

      end subroutine selectRank_SP

c     selectRank_POL
c     #####################################################################
      subroutine selectRank_POL(split_key)

c     ---------------------------------------------------------------------
c     Groups processors surrounding singular point
c     (SP boundary condition)
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: split_key

c     Local variables

        real(8)    :: random
        integer    :: dim,loc,lrank

        logical    :: cp(0:1)

c     Begin program

        dim = 3

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
        !(organize in 2D processor mesh: i + nx*(j-1))
        split_key = split_key + npz*l_lim(3,1)/g_lim(3,1)
     .                        + npz*npx*l_lim(1,1)/g_lim(1,1)

c     End program

      end subroutine selectRank_POL

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
      subroutine fillLocalVec(array,x,igr)

c     -----------------------------------------------------------------
c     Transfers outer layer of values from array to vector x.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igr
        real(8) :: array(dactx(igr)%ilm:dactx(igr)%ihp
     .                  ,dactx(igr)%jlm:dactx(igr)%jhp 
     .                  ,dactx(igr)%klm:dactx(igr)%khp)
     .            ,x    (dactx(igr)%gxs:dactx(igr)%gxe 
     .                  ,dactx(igr)%gys:dactx(igr)%gye 
     .                  ,dactx(igr)%gzs:dactx(igr)%gze)

c     Local variables

c     Begin program

        x(dactx(igr)%xs
     .   ,dactx(igr)%ys:dactx(igr)%ye 
     .   ,dactx(igr)%zs:dactx(igr)%ze)
     .                       = array(dactx(igr)%il
     .                              ,dactx(igr)%jl:dactx(igr)%jh 
     .                              ,dactx(igr)%kl:dactx(igr)%kh)

        if (dactx(igr)%xe > dactx(igr)%xs)
     .    x(dactx(igr)%xe 
     .     ,dactx(igr)%ys:dactx(igr)%ye 
     .     ,dactx(igr)%zs:dactx(igr)%ze) =
     .                         array(dactx(igr)%ih
     .                              ,dactx(igr)%jl:dactx(igr)%jh 
     .                              ,dactx(igr)%kl:dactx(igr)%kh)

        x(dactx(igr)%xs:dactx(igr)%xe 
     .   ,dactx(igr)%ys
     .   ,dactx(igr)%zs:dactx(igr)%ze) =
     .                         array(dactx(igr)%il:dactx(igr)%ih
     .                              ,dactx(igr)%jl
     .                              ,dactx(igr)%kl:dactx(igr)%kh)

        if (dactx(igr)%ye > dactx(igr)%ys)
     .    x(dactx(igr)%xs:dactx(igr)%xe 
     .     ,dactx(igr)%ye 
     .     ,dactx(igr)%zs:dactx(igr)%ze) =
     .                         array(dactx(igr)%il:dactx(igr)%ih
     .                              ,dactx(igr)%jh 
     .                              ,dactx(igr)%kl:dactx(igr)%kh)

        x(dactx(igr)%xs:dactx(igr)%xe 
     .   ,dactx(igr)%ys:dactx(igr)%ye 
     .   ,dactx(igr)%zs              ) =
     .                         array(dactx(igr)%il:dactx(igr)%ih
     .                              ,dactx(igr)%jl:dactx(igr)%jh 
     .                              ,dactx(igr)%kl              )

        if (dactx(igr)%ze > dactx(igr)%zs)
     .    x(dactx(igr)%xs:dactx(igr)%xe 
     .     ,dactx(igr)%ys:dactx(igr)%ye 
     .     ,dactx(igr)%ze              ) =
     .                         array(dactx(igr)%il:dactx(igr)%ih
     .                              ,dactx(igr)%jl:dactx(igr)%jh 
     .                              ,dactx(igr)%kh              )

cc        x(dactx(igr)%xs:dactx(igr)%xe 
cc     .   ,dactx(igr)%ys:dactx(igr)%ye 
cc     .   ,dactx(igr)%zs:dactx(igr)%ze) =
cc     .                         array(dactx(igr)%il:dactx(igr)%ih
cc     .                              ,dactx(igr)%jl:dactx(igr)%jh 
cc     .                              ,dactx(igr)%kl:dactx(igr)%kh)

c     End program

      end subroutine fillLocalVec

c     emptyLocalVec
c     #################################################################
      subroutine emptyLocalVec(array,x,igr)

c     -----------------------------------------------------------------
c     Transfers ghost cells from x to array.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igr

        real(8) :: array(dactx(igr)%ilm:dactx(igr)%ihp
     .                  ,dactx(igr)%jlm:dactx(igr)%jhp 
     .                  ,dactx(igr)%klm:dactx(igr)%khp)
     .            ,x    (dactx(igr)%gxs:dactx(igr)%gxe 
     .                  ,dactx(igr)%gys:dactx(igr)%gye 
     .                  ,dactx(igr)%gzs:dactx(igr)%gze)

c     Local variables

c     Begin program

        array(dactx(igr)%lgxs
     .       ,dactx(igr)%lgys:dactx(igr)%lgye
     .       ,dactx(igr)%lgzs:dactx(igr)%lgze) =
     .                                 x(dactx(igr)%gxs
     .                                  ,dactx(igr)%gys:dactx(igr)%gye 
     .                                  ,dactx(igr)%gzs:dactx(igr)%gze)

        if (dactx(igr)%lgxe > dactx(igr)%lgxs) 
     .    array(dactx(igr)%lgxe
     .         ,dactx(igr)%lgys:dactx(igr)%lgye
     .         ,dactx(igr)%lgzs:dactx(igr)%lgze) =
     .                                 x(dactx(igr)%gxe 
     .                                  ,dactx(igr)%gys:dactx(igr)%gye 
     .                                  ,dactx(igr)%gzs:dactx(igr)%gze)

        array(dactx(igr)%lgxs:dactx(igr)%lgxe
     .       ,dactx(igr)%lgys
     .       ,dactx(igr)%lgzs:dactx(igr)%lgze) =
     .                                 x(dactx(igr)%gxs:dactx(igr)%gxe 
     .                                  ,dactx(igr)%gys
     .                                  ,dactx(igr)%gzs:dactx(igr)%gze)

        if (dactx(igr)%lgye > dactx(igr)%lgys) 
     .    array(dactx(igr)%lgxs:dactx(igr)%lgxe
     .         ,dactx(igr)%lgye
     .         ,dactx(igr)%lgzs:dactx(igr)%lgze) =
     .                                 x(dactx(igr)%gxs:dactx(igr)%gxe 
     .                                  ,dactx(igr)%gye 
     .                                  ,dactx(igr)%gzs:dactx(igr)%gze)

        array(dactx(igr)%lgxs:dactx(igr)%lgxe
     .       ,dactx(igr)%lgys:dactx(igr)%lgye
     .       ,dactx(igr)%lgzs                ) =
     .                                 x(dactx(igr)%gxs:dactx(igr)%gxe 
     .                                  ,dactx(igr)%gys:dactx(igr)%gye 
     .                                  ,dactx(igr)%gzs               )

        if (dactx(igr)%lgze > dactx(igr)%lgzs) 
     .    array(dactx(igr)%lgxs:dactx(igr)%lgxe
     .         ,dactx(igr)%lgys:dactx(igr)%lgye
     .         ,dactx(igr)%lgze                ) =
     .                                 x(dactx(igr)%gxs:dactx(igr)%gxe 
     .                                  ,dactx(igr)%gys:dactx(igr)%gye 
     .                                  ,dactx(igr)%gze               )

cc        array(dactx(igr)%lgxs:dactx(igr)%lgxe
cc     .       ,dactx(igr)%lgys:dactx(igr)%lgye
cc     .       ,dactx(igr)%lgzs:dactx(igr)%lgze) =
cc     .                                 x(dactx(igr)%gxs:dactx(igr)%gxe 
cc     .                                  ,dactx(igr)%gys:dactx(igr)%gye 
cc     .                                  ,dactx(igr)%gzs:dactx(igr)%gze)

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

        integer       :: np=1,my_rank=0,npx=1,npy=1,npz=1
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

      integer    :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq
      real(8)    :: ldot

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

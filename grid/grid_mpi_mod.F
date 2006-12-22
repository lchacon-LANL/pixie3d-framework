c module grid_mpi
c ######################################################################
      module grid_mpi

        use grid_def

#if defined(petsc)

cc        use grid_anal_map

        implicit none

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscda.h"
#include "finclude/petscvec.h90"

        integer(4) :: np,inp,my_rank,mpierr,group_world
     .               ,group_sp,tag=0,dest=0,root=0
     .               ,status(MPI_STATUS_SIZE)
     .               ,tag_send,tag_recv,request
     .               ,npx=0,npy=0,npz=0

        integer(4) :: l_lim(3,0:1),g_lim(3,0:1),rem_l_lim(3,0:1)
     .               ,rp_l_lim(3,0:1)

        integer(4) :: MPI_COMM_SP
        integer(4),allocatable,dimension(:) :: MPI_COMM_PER
        integer(4),allocatable,dimension(:) :: MPI_COMM_NBRS

        type :: da_ctx
cc#define PETSC_AVOID_DECLARATIONS
cc#include "finclude/petsc.h"
cc#include "finclude/petscvec.h"
cc#include "finclude/petscda.h"
cc#undef PETSC_AVOID_DECLARATIONS
          DA         :: da
          Vec        :: Xg
          integer(4) :: xs,xe,xm,gxs,gxe,gxm
          integer(4) :: ys,ye,ym,gys,gye,gym
          integer(4) :: zs,ze,zm,gzs,gze,gzm
          integer(4) :: lxs,lxe,lgxs,lgxe
          integer(4) :: lys,lye,lgys,lgye
          integer(4) :: lzs,lze,lgzs,lgze
          integer(4) :: il ,ih ,jl ,jh ,kl ,kh
          integer(4) :: ilm,ihp,jlm,jhp,klm,khp
          integer(4) :: mx,my,mz,rank
          integer(4) :: igx,igy,igz
        end type da_ctx

        type (da_ctx),dimension(20) :: dactx

        character(80) :: messg

      contains

c     initMPI
c     #####################################################################
      subroutine initMPI(nxg,nyg,nzg)

c     ---------------------------------------------------------------------
c     Initialize needed MPI constructs
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nxg,nyg,nzg

c     Local variables

        integer(4) :: dim

c     Begin program

c     Rank and number of processors

        call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
        call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)

c     Allocate processors

        call processorAlloc(nxg,nyg,nzg)

      end subroutine initMPI

c     processorAlloc
c     ######################################################################
      subroutine processorAlloc(nx1,ny1,nz1)

      implicit none

c     Call variables

      integer(4) :: nx1,ny1,nz1

c     Local variables

      integer(4) :: nx,ny,nz,sum_exp,nd,navg,exp(3),npt(3),nprocs

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

      if (sum_exp < 1) return  !Only one processor

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
      if (bcond(1) == SP) ny = 1

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
      navg = floor((nx*ny*nz/nprocs)**(1./nd))

      npt = (/nx,ny,nz/)
      exp = floor(log(max(1d0*npt/navg,1d0))/log(2d0))

cc      if (my_rank == 0) write (*,*) npt
cc      if (my_rank == 0) write (*,*) exp

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

cc      if (my_rank == 0) write (*,*) sum_exp,nd,navg
cc      if (my_rank == 0) write (*,*) npx,npy,npz
cc      call PetscEnd(mpierr)
cc      stop

c     Store processor topology for ASM PC treatment

      asm_dir(1) = (npx > 1)
      asm_dir(2) = (npy > 1)
      asm_dir(3) = (npz > 1)

c     End program

      end subroutine processorAlloc

c     createBoundaryComm
c     #####################################################################
      subroutine createBoundaryComm(dim,BC,MPI_COMM)

c     ---------------------------------------------------------------------
c     Creates MPI communicator for a given boundary condition BC at the
c     calling processor. If BC=0, if builds a communicator for the neighboring
c     domains to a given processor.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: dim,MPI_COMM,BC

c     Local variables

        integer(4) :: split_key

c     Begin program

c     Find split key

        split_key = MPI_UNDEFINED

        call selectRank(dim,BC,split_key)

c     Create communicator

        call MPI_Comm_split(MPI_COMM_WORLD,split_key,my_rank
     .                     ,MPI_COMM,mpierr)

c     End programs

      end subroutine createBoundaryComm

c     selectRank
c     #####################################################################
      subroutine selectRank(dim,BC,split_key)

c     ---------------------------------------------------------------------
c     Finds number of processors in MPI communicator
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: dim,BC,split_key

c     Local variables

        real(8)    :: random
        integer(4) :: loc,lrank

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

        do loc=0,1
          if (checkBCLim(BC,dim,loc,l_lim,g_lim)) split_key = BC
        enddo

        !Create self-communicator for single-point dimensions
        if (l_lim(dim,0) == l_lim(dim,1) .and. split_key == BC) then
          call random_seed()
          call random_number(random)
          split_key=split_key+my_rank+1
        endif

c     End program

      end subroutine selectRank

c     checkBCLim
c     #####################################################################
      function checkBCLim(BC,dim,loc,llim,glim) result(include_proc)

c     ---------------------------------------------------------------------
c     Checks whether a processor should be included in an MPI_COMM according
c     to the boundary type BC, and the local vs. global limits
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: BC,dim,loc,llim(3,0:1),glim(3,0:1)
        logical    :: include_proc

c     Local variables

        integer(4) :: ibc

c     Begin program

        ibc(dim,loc) = (1+loc)+2*(dim-1)

        include_proc= (llim(dim,loc)==glim(dim,loc))
     .                 .and.bcond(ibc(dim,loc))==BC

      end function checkBCLim

c     createPETScGrid
c     #################################################################
      subroutine createPETScGrid

        implicit none

c     -----------------------------------------------------------------
c     Creates PETSc distributed array for boundary communication in
c     fortran.
c     -----------------------------------------------------------------

c     Call variables

c     Local variables

        integer(4) :: igr

c     Begin program

c     Create PETSc DA hierachy

        do igr=1,grid_params%ngrid
          call createFortranDA(igr)
        enddo

c     Singular-point boundary MPI communicator

        if (bcSP()) call createBoundaryComm(1,bcond(1),MPI_COMM_SP)

c diag ****
cc        write (*,*) 'Proc',my_rank,', SP COMM' ,MPI_COMM_SP
cc        if (MPI_COMM_SP /= MPI_COMM_NULL) then
cc          call MPI_Comm_size(MPI_COMM_SP,my_rank,mpierr)
cc          write (*,*) '# ranks SP_COMM ',my_rank
cc        endif
cc        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
cc        call MPI_Finalize(mpierr)
cc        stop
c diag ****

c     End program

      end subroutine createPETScGrid

c     createFortranDA
c     #################################################################
      subroutine createFortranDA(igrid)

c     -----------------------------------------------------------------
c     Creates PETSc distributed array for boundary communication in
c     fortran at grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igrid

c     Local variables

        integer(4) :: BC,ierr,nxg,nyg,nzg,igx,igy,igz

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        nxg = grid_params%nxgl(igx)
        nyg = grid_params%nygl(igy)
        nzg = grid_params%nzgl(igz)

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

        dactx(igrid)%igx = igx
        dactx(igrid)%igy = igy
        dactx(igrid)%igz = igz

c     Define local limits

        !With ghost cells (only those that PETSc includes)
        call fromGlobalToLocalLimits
     .       (dactx(igrid)%gxs ,dactx(igrid)%gys ,dactx(igrid)%gzs
     $       ,dactx(igrid)%lgxs,dactx(igrid)%lgys,dactx(igrid)%lgzs
     $       ,dactx(igrid)%igx ,dactx(igrid)%igy ,dactx(igrid)%igz)
        call fromGlobalToLocalLimits
     .       (dactx(igrid)%gxe ,dactx(igrid)%gye ,dactx(igrid)%gze
     $       ,dactx(igrid)%lgxe,dactx(igrid)%lgye,dactx(igrid)%lgze
     $       ,dactx(igrid)%igx ,dactx(igrid)%igy ,dactx(igrid)%igz)

        !Array limits
        dactx(igrid)%il = 1
        dactx(igrid)%jl = 1
        dactx(igrid)%kl = 1
        dactx(igrid)%ih = grid_params%nxv(igx)
        dactx(igrid)%jh = grid_params%nyv(igy)
        dactx(igrid)%kh = grid_params%nzv(igz)

        dactx(igrid)%ilm = dactx(igrid)%il-1
        dactx(igrid)%ihp = dactx(igrid)%ih+1
        dactx(igrid)%jlm = dactx(igrid)%jl-1
        dactx(igrid)%jhp = dactx(igrid)%jh+1
        dactx(igrid)%klm = dactx(igrid)%kl-1
        dactx(igrid)%khp = dactx(igrid)%kh+1

c     End program

      end subroutine createFortranDA

c     destroyDA
c     #################################################################
      subroutine destroyDA

        implicit none

c     Call variables

c     Local variables

        integer(4) :: ierr,igr

c     Begin program

        do igr=1,grid_params%ngrid
          call VecDestroy(dactx(igr)%Xg,ierr)
          call DADestroy (dactx(igr)%da,ierr)
        enddo

c     End program

      end subroutine destroyDA

c     fillLocalVec
c     #################################################################
      subroutine fillLocalVec(array,x,igr)

        implicit none

c     Call variables

        integer(4) :: igr
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

        implicit none

c     Call variables

        integer(4) :: igr

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

c     fillPetscGhostCells
c     #################################################################
      subroutine fillPetscGhostCells(array,igr)

        implicit none

c     Call variables

        integer(4) :: igr
        real(8) :: array(dactx(igr)%ilm:dactx(igr)%ihp
     .                  ,dactx(igr)%jlm:dactx(igr)%jhp 
     .                  ,dactx(igr)%klm:dactx(igr)%khp)

c     Local variables

        Vec :: localX

        PetscScalar,pointer :: lx_v(:)

        integer(4) :: ierr

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

c     fchkerrq
c     #################################################################
      subroutine fchkerrq(ierr,routine)

        implicit none

c     Call variables

        integer(4) :: ierr
        character(*) :: routine

c     Begin program

        if (ierr /= 0) call pstop(routine,'PETSc error')

c     End program

      end subroutine fchkerrq

c     isProc
c     #################################################################
      function isProc(proc_id)

        implicit none

c     Call variables

        integer(4) :: proc_id
        logical    :: isProc

c     Local variables

c     Begin program

cc        call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)

        isProc = (my_rank == proc_id)

c     End program

      end function isProc

c     inProc
c     #################################################################
      function inProc(igl,jgl,kgl,igx,igy,igz)

        implicit none

c     Call variables

        integer(4) :: igl,jgl,kgl,igx,igy,igz
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

      integer(4) :: ntot
      real(8)    :: vec1(ntot),vec2(ntot),dot

c     Local variables

      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
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

#else

        integer(4)    :: np=1,my_rank=0,npx=0,npy=0,npz=0
        character(80) :: messg

      contains

c     isProc
c     #################################################################
      function isProc(proc_id)

        implicit none

c     Call variables

        integer(4) :: proc_id
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

        integer(4) :: igl,jgl,kgl,igx,igy,igz
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

      integer(4) :: ntot
      real(8)    :: vec1(ntot),vec2(ntot),dot

c     Local variables

      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
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

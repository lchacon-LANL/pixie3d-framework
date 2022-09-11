c module bc_def
c #####################################################################
      module bc_def

c     -----------------------------------------------------------------
c     This module defines the types of BCs possible:
c     Topological BCs (set topology of computational domain):
c       PER: Periodic
c       SP : Singular point (as in cylindrical, toroidal)
c       SYM: Symmetry (> 0 -> homogeneous Neumann
c                     ,< 0 -> homogeneous Dirichlet)
c
c     Standard BCs:
c       EQU: Inhomogeneous dirichlet (generally imposed by equilibrium)
c       DIR: Homogeneous dirichlet
c       NEU: Neumann (homogeneous/inhomogeneous)
c       SYM: Symmetry at cell centers
c       FSYM: Symmetry at faces
c       EXT: By extrapolation
c       IFC: Interpolate to face (and store in ghost cell)
c       DEF: Default (does nothing)
c     The actual integers chosen indicate a particular order that the BC
c     scheduler follows to impose the boundary conditions.
c     -----------------------------------------------------------------

        implicit none

        integer,parameter :: DEF=0   !Default (do nothing)
     .                      ,PER=1   !Periodic  (do not change: NEMATO depends on it)
     .                      ,EQU=2   !Equilibrium (inhomogeneous Dirichlet)
     .                      ,EXT=3   !Extrapolation
     .                      ,NEU=4   !Neumann or Robin
     .                      ,DIR=5   !Homogeneous Dirichlet
     .                      ,IFC=6   !Interpolate-to-face
     .                      ,SP =7   !Polar singular point (do not change: NEMATO depends on it)
     .                      ,SYM=8   !Symmetry (even > 0, odd < 0)
     .                      ,FSYM=9  !Symmetry (at faces; even > 0, odd < 0)
     .                      ,OPN=10  !Open boundary conditions
     .                      ,SYNC=11 !Parallel sync

        integer,parameter :: BCLIM=7 !Defines # of BC cases (excludes SYNC,SYM,FSYM,OPN)
     .                      ,BCTOT=10!Total number of BC defs

        logical :: asm_dir(3)=.false.,asm=.false.

        integer :: bcond(6)          !Contains topological BC configuration

        real(8) :: SP_flsv=1d-5      !Singular point failsafe

      contains

c     initBCs
c     #################################################################
      function initBCs(bcs) result(ierr)

c     -----------------------------------------------------------------
c     Initializes BCs according to character codes in bcs variable
c     On error, it returns ierr =/ 0. Otherwise, ierr == 0.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: ierr
        character*(3) :: bcs(6)

c     Local variables

c     Begin program

        ierr = 0

        bcond = -1

        where (bcs == 'def')
          bcond = DEF
        elsewhere (bcs == 'per')
          bcond = PER
        elsewhere (bcs == 'spt')
          bcond = SP
        elsewhere (bcs == 'sym')
          bcond = SYM
        elsewhere (bcs == 'fsm')
          bcond = FSYM
        elsewhere (bcs == 'equ')
          bcond = EQU
        elsewhere (bcs == 'opn')
          bcond = OPN
        end where

        if (minval(bcond) < 0) ierr = 1

      end function initBCs

c     bcChk
c     #################################################################
      function bcChk(ibc,bctype,bcnd) 

c     -----------------------------------------------------------------
c     Detects whether boundary is of type bctype. Boundary is identified
c     by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        logical :: bcChk
        integer :: ibc,bctype
        integer :: bcnd(6)

c     Local variables

c     Begin program

        bcChk = (ibc>=1.and.ibc<=6)
        if (bcChk) bcChk = (abs(bcnd(ibc)) == bctype)

      end function bcChk

c     bcSP
c     #################################################################
      function bcSP(dim,ibc) result(sing_bc)

c     -----------------------------------------------------------------
c     Detects whether boundary is of SP type. Boundary is identified
c     by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer,optional :: dim,ibc
        logical :: sing_bc

c     Local variables

        integer :: iibc

c     Begin program

        if (PRESENT(ibc)) then
          iibc = ibc
        elseif (PRESENT(dim)) then
          iibc = 1+2*(dim-1)
        else
          iibc = 1
        endif

        sing_bc = bcChk(iibc,SP,bcond)

      end function bcSP

c     bcPER
c     #################################################################
      function bcPER(dim) result(per_bc)

c     -----------------------------------------------------------------
c     Detects whether dimension dim is of PER type.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: dim 
        logical :: per_bc

c     Local variables

        integer :: bc(6),iibc

c     Begin program

        iibc = 1+2*(dim-1)

        per_bc = bcChk(iibc,PER,bcond)

      end function bcPER

c     bcSYM
c     #################################################################
      function bcSYM(ibc) result(sym_bc)

c     -----------------------------------------------------------------
c     Detects whether boundary is of SYM type. Boundary is
c     identified by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: ibc
        logical :: sym_bc

c     Local variables

c     Begin program

        sym_bc = bcChk(ibc,SYM,bcond)

      end function bcSYM

c     bcFSYM
c     #################################################################
      function bcFSYM(ibc) result(sym_bc)

c     -----------------------------------------------------------------
c     Detects whether boundary is of SYM type. Boundary is
c     identified by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: ibc
        logical :: sym_bc

c     Local variables

c     Begin program

        sym_bc = bcChk(ibc,FSYM,bcond)

      end function bcFSYM

c     get_topol_bcs
c     #################################################################
      subroutine get_topol_bcs(bcs)

      integer :: bcs(6)

      bcs = bcond

      end subroutine get_topol_bcs

      end module bc_def

c module profile_petsc
c ######################################################################
      module profile_petsc

       integer,parameter:: iSTUP=1,iNL=2,iPCSTP=3,iPCSLV=4,iPOST=5,iMG=6

       integer :: prof_err,prof_stage(10),iprofile
       real(8) :: prof_time(20)
       logical,private :: prof_init=.true.

#if defined(petsc)
      contains

c     init_petsc_profile
c     ##################################################################
      subroutine init_petsc_profile

      implicit none

      !Set profiling stages
      call PetscLogStageRegister('Setup   ',prof_stage(iSTUP) ,prof_err)
      call PetscLogStageRegister('NL solve',prof_stage(iNL)   ,prof_err)
      call PetscLogStageRegister('PC setup',prof_stage(iPCSTP),prof_err)  
      call PetscLogStageRegister('PC main ',prof_stage(iPCSLV),prof_err)  
      call PetscLogStageRegister('Post+IO ',prof_stage(iPOST) ,prof_err)  
      call PetscLogStageRegister('MG solve',prof_stage(iMG)   ,prof_err)  

      iprofile = 6

      end subroutine init_petsc_profile

c     add_profile_stage
c     ##################################################################
      function add_profile_stage(name) result(handle)

      implicit none

      character(*) :: name
      integer :: handle

      !Set profiling stages
      iprofile = iprofile + 1

      call PetscLogStageRegister(trim(name),prof_stage(iprofile)
     .                          ,prof_err)

      handle = iprofile

      end function add_profile_stage

c     switch_profile_stage
c     ##################################################################
      subroutine switch_profile_stage(handle)

      implicit none

      integer :: handle

      if (.not.prof_init) then
        call PetscLogStagePop(prof_err)
      else
        prof_init = .false.
      endif

      call PetscLogStagePush(handle,prof_err)

      end subroutine switch_profile_stage

c     end_petsc_profile
c     ##################################################################
      subroutine end_petsc_profile

      implicit none

      call PetscLogStagePop(prof_err)

      end subroutine end_petsc_profile
#endif

      end module profile_petsc

c module grid_petsc
c ######################################################################
      module grid_petsc

        use profile_petsc

        use bc_def

#if defined(petsc)

        use mpi

cc#define PETSC_AVOID_MPIF_H
        
#include <petsc/finclude/petscvec.h>
        use Petscvec
        
#include <petsc/finclude/petscdmda.h>
        use Petscdmda

        type :: petsc_da_ctx
          DM         :: da
          Vec        :: Xg
          integer    :: xs,xe,xm,gxs,gxe,gxm
          integer    :: ys,ye,ym,gys,gye,gym
          integer    :: zs,ze,zm,gzs,gze,gzm
          integer    :: lxs,lxe,lgxs,lgxe
          integer    :: lys,lye,lgys,lgye
          integer    :: lzs,lze,lgzs,lgze
          integer    :: il ,ih ,jl ,jh ,kl ,kh
          integer    :: ilm,ihp,jlm,jhp,klm,khp
          integer    :: mx,my,mz        !Number of processors
          integer    :: nxgl,nygl,nzgl  !Global number of mesh points
          integer    :: igrid
          integer    :: MPI_COMM_DA,np,rank
          integer    :: proc_neighbors(-1:1,-1:1,-1:1)
        end type petsc_da_ctx

        integer :: debug_par=0        !For parallel debugging

      contains

c     createPETScFortranDA
c     #################################################################
      subroutine createPETScFortranDA(dactx)

c     -----------------------------------------------------------------
c     Creates PETSc distributed array (DA) for boundary communication
c     in fortran at grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(petsc_da_ctx),pointer :: dactx

c     Local variables

        integer :: BC1,BC2,BC3,ierr,mpierr,i,ix,iy,iz,neighbors(27)
        
c     Begin program

c$$$        if (dactx%mx == 0) dactx%mx = PETSC_DECIDE
c$$$        if (dactx%my == 0) dactx%my = PETSC_DECIDE
c$$$        if (dactx%mz == 0) dactx%mz = PETSC_DECIDE
        if (dactx%mx == 0 .or. dactx%my == 0 .or. dactx%mz == 0) then
          STOP "Unspecified proc. alloc. in createPETScFortranDA"
        endif
          
        !Serial periodic BCs no longer done by PETSc
        BC1 = DM_BOUNDARY_NONE
        if (bcond(1)==PER.and.dactx%mx > 1) BC1 = DM_BOUNDARY_PERIODIC
        BC2 = DM_BOUNDARY_NONE
        if (bcond(3)==PER.and.dactx%my > 1) BC2 = DM_BOUNDARY_PERIODIC
        BC3 = DM_BOUNDARY_NONE
        if (bcond(5)==PER.and.dactx%mz > 1) BC3 = DM_BOUNDARY_PERIODIC

        call DMDACreate3d(dactx%MPI_COMM_DA,BC1,BC2,BC3,DMDA_STENCIL_BOX
     $                   ,dactx%nxgl,dactx%nygl,dactx%nzgl
     &                   ,dactx%mx,dactx%my,dactx%mz,1,1
     &                   ,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER
     &                   ,PETSC_NULL_INTEGER,dactx%da,ierr)
        call DMSetUp(dactx%da,ierr)

        call DMCreateGlobalVector(dactx%da,dactx%Xg,ierr)

c       Get local grid boundaries in global grid

        call DMDAGetCorners(dactx%da
     .         ,dactx%xs,dactx%ys,dactx%zs
     .         ,dactx%xm,dactx%ym,dactx%zm,ierr)
        call DMDAGetGhostCorners(dactx%da
     .         ,dactx%gxs,dactx%gys,dactx%gzs
     .         ,dactx%gxm,dactx%gym,dactx%gzm,ierr)

        dactx%xs  = dactx%xs+1
        dactx%ys  = dactx%ys+1
        dactx%zs  = dactx%zs+1
        dactx%gxs = dactx%gxs+1
        dactx%gys = dactx%gys+1
        dactx%gzs = dactx%gzs+1

        dactx%ye  = dactx%ys+dactx%ym-1
        dactx%xe  = dactx%xs+dactx%xm-1
        dactx%ze  = dactx%zs+dactx%zm-1
        dactx%gye = dactx%gys+dactx%gym-1
        dactx%gxe = dactx%gxs+dactx%gxm-1
        dactx%gze = dactx%gzs+dactx%gzm-1

        !With ghost cells (only those that PETSc includes)
        dactx%lgxs = dactx%gxs - dactx%xs + 1
        dactx%lgys = dactx%gys - dactx%ys + 1
        dactx%lgzs = dactx%gzs - dactx%zs + 1
        dactx%lgxe = dactx%gxe - dactx%xs + 1
        dactx%lgye = dactx%gye - dactx%ys + 1
        dactx%lgze = dactx%gze - dactx%zs + 1

        !Domain limits (without ghost cells)
        dactx%lxs = 1
        dactx%lys = 1
        dactx%lzs = 1
        dactx%lxe = dactx%xm
        dactx%lye = dactx%ym
        dactx%lze = dactx%zm

        !Array limits (with all ghost cells)
        dactx%il = dactx%lxs
        dactx%jl = dactx%lys
        dactx%kl = dactx%lzs
        dactx%ih = dactx%lxe
        dactx%jh = dactx%lye
        dactx%kh = dactx%lze

        dactx%ilm = dactx%il-1
        dactx%ihp = dactx%ih+1
        dactx%jlm = dactx%jl-1
        dactx%jhp = dactx%jh+1
        dactx%klm = dactx%kl-1
        dactx%khp = dactx%kh+1

c     Store MPI info

        call MPI_Comm_rank(dactx%MPI_COMM_DA,dactx%rank,mpierr)
        call MPI_Comm_size(dactx%MPI_COMM_DA,dactx%np  ,mpierr)

        call DMDAGetNeighbors(dactx%da,neighbors,mpierr)

        where (neighbors < 0) neighbors = MPI_PROC_NULL
        
        dactx%proc_neighbors = reshape(neighbors,(/3,3,3/))

        !Make "self" a fake neighbor in serial directions
        if (bcond(1)==PER.and.dactx%mx == 1) then
          dactx%proc_neighbors(-1,0,0) = dactx%proc_neighbors(0,0,0)
          dactx%proc_neighbors( 1,0,0) = dactx%proc_neighbors(0,0,0)
        endif
        if (bcond(3)==PER.and.dactx%my == 1) then
          dactx%proc_neighbors(0,-1,0) = dactx%proc_neighbors(0,0,0)
          dactx%proc_neighbors(0, 1,0) = dactx%proc_neighbors(0,0,0)
        endif
        if (bcond(5)==PER.and.dactx%mz == 1) then
          dactx%proc_neighbors(0,0,-1) = dactx%proc_neighbors(0,0,0)
          dactx%proc_neighbors(0,0, 1) = dactx%proc_neighbors(0,0,0)
        endif

c$$$        write (*,*) dactx%rank,MPI_PROC_NULL
c$$$        write (*,*) dactx%rank,neighbors
c$$$        write (*,*) dactx%rank,dactx%proc_neighbors( 0  , 1,0)
c$$$        write (*,*) dactx%rank,dactx%proc_neighbors(-1:1, 0,0)
c$$$        write (*,*) dactx%rank,dactx%proc_neighbors( 0  ,-1,0)
c$$$        stop

c     End program

      end subroutine createPETScFortranDA

c     destroyPETScFortranDA
c     #################################################################
      subroutine destroyPETScFortranDA(dactx)

c     -----------------------------------------------------------------
c     Deallocates PETSc DA
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(petsc_da_ctx),pointer :: dactx

c     Local variables

        integer    :: ierr

c     Begin program

        call VecDestroy(dactx%Xg,ierr)
        call DMDestroy (dactx%da,ierr)

c     End program

      end subroutine destroyPETScFortranDA

#endif

      end module grid_petsc

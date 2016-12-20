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

        integer,parameter :: BCLIM=7 !Defines # of BC cases (excludes DEF,SYM,FSYM,OPN)
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
c     Detects whether boundary is of SP type. Boundary is identified
c     by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        logical :: bcChk
        integer :: ibc,bctype
        integer,optional :: bcnd(6)

c     Local variables

        integer :: bc(6)

c     Begin program

        if (PRESENT(bcnd)) then
          bc = bcnd
        else
          bc = bcond
        endif

        if (ibc < 1 .or. ibc > 6) then
          bcChk = .false.
        else
          bcChk = (abs(bc(ibc)) == bctype)
        endif

      end function bcChk

c     bcSP
c     #################################################################
      function bcSP(ibc,bcnd) result(sing_bc)

c     -----------------------------------------------------------------
c     Detects whether boundary is of SP type. Boundary is identified
c     by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer,optional :: ibc,bcnd(6)
        logical :: sing_bc

c     Local variables

        integer :: iibc

c     Begin program

        if (PRESENT(ibc)) then
          iibc = ibc
        else
          iibc = 1
        endif

        sing_bc = bcChk(iibc,SP,bcnd=bcnd)

      end function bcSP

c     bcPER
c     #################################################################
      function bcPER(dim,bcnd) result(per_bc)

c     -----------------------------------------------------------------
c     Detects whether dimension dim is of PER type.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: dim 
        logical :: per_bc
        integer,optional :: bcnd(6)

c     Local variables

        integer :: bc(6),iibc

c     Begin program

        iibc = 1+2*(dim-1)

        per_bc = bcChk(iibc,PER,bcnd=bcnd)

      end function bcPER

c     bcSYM
c     #################################################################
      function bcSYM(ibc,bcnd) result(sym_bc)

c     -----------------------------------------------------------------
c     Detects whether boundary is of SYM or FSYM type. Boundary is
c     identified by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: ibc
        logical :: sym_bc
        integer,optional :: bcnd(6)

c     Local variables

        integer :: bc(6)

c     Begin program

        sym_bc = bcChk(ibc,SYM,bcnd=bcnd)

cc        if (PRESENT(bcnd)) then
cc          bc = bcnd
cc        else
cc          bc = bcond
cc        endif
cc
cc        if (ibc < 1 .or. ibc > 6) then
cc          sym_bc = .false.
cc        else
cc          sym_bc =    (abs(bc(ibc)) ==  SYM)
cccc     .            .or.(abs(bc(ibc)) == FSYM)
cc        endif

      end function bcSYM

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
      call PetscLogStageRegister('Setup'   ,prof_stage(iSTUP) ,prof_err)
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

        implicit none

#if defined(petsc)

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscda.h"
#include "finclude/petscvec.h90"

        type :: petsc_da_ctx
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
          integer    :: mx,my,mz        !Number of processors
          integer    :: nxgl,nygl,nzgl  !Global number of mesh points
          integer    :: igrid
          integer    :: MPI_COMM_DA,np,rank
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

        integer :: BC,ierr,mpierr

c     Begin program

        if     (bcond(1)==PER.and.bcond(3)==PER.and.bcond(5)==PER) then
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

        if (dactx%mx == 0) dactx%mx = PETSC_DECIDE
        if (dactx%my == 0) dactx%my = PETSC_DECIDE
        if (dactx%mz == 0) dactx%mz = PETSC_DECIDE

        call DACreate3d(dactx%MPI_COMM_DA,BC,DA_STENCIL_BOX
     $                 ,dactx%nxgl,dactx%nygl,dactx%nzgl
     &                 ,dactx%mx,dactx%my,dactx%mz,1,1
     &                 ,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER
     &                 ,PETSC_NULL_INTEGER,dactx%da,ierr)

        call DACreateGlobalVector(dactx%da,dactx%Xg,ierr)

c       Get local grid boundaries in global grid

        call DAGetCorners(dactx%da
     .         ,dactx%xs,dactx%ys,dactx%zs
     .         ,dactx%xm,dactx%ym,dactx%zm,ierr)
        call DAGetGhostCorners(dactx%da
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
        call DADestroy (dactx%da,ierr)

c     End program

      end subroutine destroyPETScFortranDA

#endif

      end module grid_petsc

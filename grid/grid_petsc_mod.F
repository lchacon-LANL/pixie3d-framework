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
c       SYM: Symmetry (odd,even depend on sign)
c       EXT: By extrapolation
c       IFC: Interpolate to face (and store in ghost cell)
c       DEF: Default (does nothing)
c     The actual integers chosen indicate a particular order that the BC
c     scheduler follows to impose the boundary conditions.
c     -----------------------------------------------------------------

        implicit none

        integer :: PER,DIR,NEU,SP,EQU,DEF,EXT,SYM,IFC
        parameter (DEF=0,PER=1,EQU=2,EXT=3,NEU=4,SYM=5,DIR=6,SP=7,IFC=8)

        integer,parameter :: BCLIM=8 !Defines # of BC cases (excludes DEF)

        logical :: asm_dir(3)=.false.,asm=.false.,par_bc=.false.


        integer :: bcond(6)      !Contains topological BC configuration

        real(8)    :: SP_flsv=1d-5         !SP failsafe

      contains

c     bcSP
c     #################################################################
      function bcSP(ibc) result(sing_bc)

        implicit none

        integer   ,optional :: ibc
        logical    :: sing_bc

        integer    :: iibc

        if (PRESENT(ibc)) then
          iibc = ibc
        else
          iibc = 1
        endif

        sing_bc = (bcond(iibc) == SP)

      end function bcSP

      end module bc_def

c module grid_petsc
c ######################################################################
      module grid_petsc

        use bc_def

        implicit none

#if !defined(petsc)

        integer :: npx=1,npy=1,npz=1  !Serial processor allocation

#else

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
          integer    :: mx,my,mz
          integer    :: igrid
          integer    :: MPI_COMM_DA,np,rank
        end type petsc_da_ctx

        integer :: npx=0,npy=0,npz=0  !Processor allocation

        integer :: debug_par=0        !For parallel debuggin

      contains

c     createPETScFortranDA
c     #################################################################
      subroutine createPETScFortranDA(igrid,nxgl,nygl,nzgl,dactx)

c     -----------------------------------------------------------------
c     Creates PETSc distributed array (DA) for boundary communication
c     in fortran at grid level igrid.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,nxgl,nygl,nzgl
        type(petsc_da_ctx),pointer :: dactx

c     Local variables

        integer    :: BC,ierr,nxg,nyg,nzg,mpierr

c     Begin program

        nxg = nxgl
        nyg = nygl
        nzg = nzgl

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

        if (npx == 0) npx = PETSC_DECIDE
        if (npy == 0) npy = PETSC_DECIDE
        if (npz == 0) npz = PETSC_DECIDE

        call DACreate3d(dactx%MPI_COMM_DA,BC,DA_STENCIL_BOX,nxg,nyg,nzg
     &                 ,npx,npy,npz,1,1
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

c     Store grid level

        dactx%igrid = igrid

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

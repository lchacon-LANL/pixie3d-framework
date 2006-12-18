
c module grid_structures
c #####################################################################
      module grid_structures

        implicit none

        integer(4) :: nxgl,nygl,nzgl

        integer(4),parameter :: ngparams = 6

        type :: grid_mg_def
          integer(4) :: ngrdx                             !# meshes in X
          integer(4) :: ngrdy                             !# meshes in Y
          integer(4) :: ngrdz                             !# meshes in Y
          integer(4) :: ngrid                             !# meshes for MG
          integer(4),pointer,dimension(:)  :: iline       !Restrict ops. to i=iline in MG
          integer(4),pointer,dimension(:)  :: jline       !Restrict ops. to j=jline in MG
          integer(4),pointer,dimension(:)  :: kline       !Restrict ops. to k=kline in MG
          integer(4),pointer,dimension(:)  :: ilo         !Global lower limit in X
          integer(4),pointer,dimension(:)  :: jlo         !Global lower limit in Y
          integer(4),pointer,dimension(:)  :: klo         !Global lower limit in Z
          integer(4),pointer,dimension(:)  :: ihi         !Global higher limit in X
          integer(4),pointer,dimension(:)  :: jhi         !Global higher limit in Y
          integer(4),pointer,dimension(:)  :: khi         !Global higher limit in Z
          real(8)   ,pointer,dimension(:)  :: xg          !Global grid nodes in X (top grid)
          real(8)   ,pointer,dimension(:)  :: yg          !Global grid nodes in Y (")
          real(8)   ,pointer,dimension(:)  :: zg          !Global grid nodes in Z (")
          real(8)   ,pointer,dimension(:)  :: xx          !Local grid nodes in X (all grids)
          real(8)   ,pointer,dimension(:)  :: yy          !Local grid nodes in Y (")
          real(8)   ,pointer,dimension(:)  :: zz          !Local grid nodes in Z (")
          real(8)   ,pointer,dimension(:)  :: dx          !Grid spacings in X for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dy          !Grid spacings in Y for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dz          !Grid spacings in Z for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dxh         !Grid spacings in X for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dyh         !Grid spacings in Y for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dzh         !Grid spacings in Z for half mesh (")
          integer(4),pointer,dimension(:)  :: nxv         !Local # of grid nodes in X  (")
          integer(4),pointer,dimension(:)  :: nyv         !Local # of grid nodes in Y  (")
          integer(4),pointer,dimension(:)  :: nzv         !Local # of grid nodes in Z  (")
          integer(4),pointer,dimension(:)  :: nxgl        !Global # of grid nodes in X  (")
          integer(4),pointer,dimension(:)  :: nygl        !Global # of grid nodes in Y  (")
          integer(4),pointer,dimension(:)  :: nzgl        !Global # of grid nodes in Z  (")
cc          integer(4),pointer,dimension(:)  :: ntotv       !Total # of grid nodes (")
          integer(4),pointer,dimension(:)  :: istartx     !Pointer for MG vectors in X
          integer(4),pointer,dimension(:)  :: istarty     !Pointer for MG vectors in Y
          integer(4),pointer,dimension(:)  :: istartz     !Pointer for MG vectors in Z
          integer(4),pointer,dimension(:)  :: istartp     !Pointer for global MG vectors
          integer(4),pointer,dimension(:)  :: mg_ratio_x  !MG coarsening ratio in X
          integer(4),pointer,dimension(:)  :: mg_ratio_y  !MG coarsening ratio in Y
          integer(4),pointer,dimension(:)  :: mg_ratio_z  !MG coarsening ratio in Z
          real(8)                          :: params(ngparams)   !Grid configuration parameters
        end type grid_mg_def

        type (grid_mg_def) :: grid_params

        INTERFACE ASSIGNMENT (=)
          module procedure equateGridStructure
        END INTERFACE

      contains

c     allocateGridStructure
c     #################################################################
      subroutine allocateGridStructure(nx,ny,nz,ngridx,ngridy,ngridz
     .                                ,grid_st)
c     -----------------------------------------------------------------
c     Allocates grid structure
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4)     :: nx,ny,nz,ngridx,ngridy,ngridz
        type(grid_mg_def) :: grid_st

c     Local variables

        integer(4) :: ngrid,nxmg,nymg,nzmg

c     Begin program

        ngrid = max(ngridx,ngridy,ngridz)

        grid_st%ngrdx = ngridx
        grid_st%ngrdy = ngridy
        grid_st%ngrdz = ngridz
        grid_st%ngrid = ngrid

        nxmg=findMGsize(nx,ngridx,ngrid)
        nymg=findMGsize(ny,ngridy,ngrid)
        nzmg=findMGsize(nz,ngridz,ngrid)

        if (.not.associated(grid_st%xx)) then
          allocate(grid_st%xg(0:nxgl+1))
          allocate(grid_st%yg(0:nygl+1))
          allocate(grid_st%zg(0:nzgl+1))
          allocate(grid_st%xx(nxmg+2*ngrid))
          allocate(grid_st%yy(nymg+2*ngrid))
          allocate(grid_st%zz(nzmg+2*ngrid))
          allocate(grid_st%dx(nxmg+2*ngrid))
          allocate(grid_st%dy(nymg+2*ngrid))
          allocate(grid_st%dz(nzmg+2*ngrid))
          allocate(grid_st%dxh(nxmg+2*ngrid))
          allocate(grid_st%dyh(nymg+2*ngrid))
          allocate(grid_st%dzh(nzmg+2*ngrid))
          allocate(grid_st%nxv(ngrid))
          allocate(grid_st%nyv(ngrid))
          allocate(grid_st%nzv(ngrid))
          allocate(grid_st%nxgl(ngrid))
          allocate(grid_st%nygl(ngrid))
          allocate(grid_st%nzgl(ngrid))
          allocate(grid_st%ilo(ngrid))
          allocate(grid_st%jlo(ngrid))
          allocate(grid_st%klo(ngrid))
          allocate(grid_st%ihi(ngrid))
          allocate(grid_st%jhi(ngrid))
          allocate(grid_st%khi(ngrid))
cc          allocate(grid_st%ntotv(ngrid))
          allocate(grid_st%istartx(ngrid))
          allocate(grid_st%istarty(ngrid))
          allocate(grid_st%istartz(ngrid))
          allocate(grid_st%istartp(ngrid))
          allocate(grid_st%mg_ratio_x(ngrid))
          allocate(grid_st%mg_ratio_y(ngrid))
          allocate(grid_st%mg_ratio_z(ngrid))
          allocate(grid_st%iline(ngrid))
          allocate(grid_st%jline(ngrid))
          allocate(grid_st%kline(ngrid))
        endif

        grid_st%xg = 0.
        grid_st%yg = 0.
        grid_st%zg = 0.
        grid_st%xx = 0.
        grid_st%yy = 0.
        grid_st%zz = 0.
        grid_st%dx = 0.
        grid_st%dy = 0.
        grid_st%dz = 0.
        grid_st%dxh = 0.
        grid_st%dyh = 0.
        grid_st%dzh = 0.
        grid_st%nxv = 0
        grid_st%nyv = 0
        grid_st%nzv = 0
        grid_st%nxgl = 0
        grid_st%nygl = 0
        grid_st%nzgl = 0
        grid_st%ilo = 0
        grid_st%jlo = 0
        grid_st%klo = 0
        grid_st%ihi = 0
        grid_st%jhi = 0
        grid_st%khi = 0
        grid_st%istartx = 0
        grid_st%istarty = 0
        grid_st%istartz = 0
        grid_st%istartp = 0
        grid_st%mg_ratio_x = 0
        grid_st%mg_ratio_y = 0
        grid_st%mg_ratio_z = 0
        grid_st%iline = 0
        grid_st%jline = 0
        grid_st%kline = 0

c     End program

      contains

c     findMGsize
c     #################################################################
      function findMGsize(nn,ngrd,ngrdt) result (nnmg)
      implicit none
c     -----------------------------------------------------------------
c     Finds size for MG vectors, taking into account total grid levels
c     ngrdt, grid levels in the relevant direction ngrd, and the 
c     number of mesh points in the finest grid nn. The formula ensures
c     enough space even if nn=1 in the finest grid. It does NOT include
c     ghost cells (this requires an additional term of 2*ngrdt).
c     -----------------------------------------------------------------

        integer(4) :: nn,nnmg,ngrd,ngrdt

cc        nnmg = 2*nn + ngrdt - ngrd -1
        nnmg = (2+ngrdt-ngrd)*nn

      end function findMGsize

      end subroutine allocateGridStructure

c     deallocateGridStructure
c     #################################################################
      subroutine deallocateGridStructure(grid_st)
c     -----------------------------------------------------------------
c     Allocates grid structure
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def) :: grid_st

c     Begin program

        grid_st%ngrdx = 0
        grid_st%ngrdy = 0
        grid_st%ngrdz = 0
        grid_st%ngrid = 0

        if (associated(grid_st%xx)) then
          deallocate(grid_st%xg)
          deallocate(grid_st%yg)
          deallocate(grid_st%zg)
          deallocate(grid_st%xx)
          deallocate(grid_st%yy)
          deallocate(grid_st%zz)
          deallocate(grid_st%dx)
          deallocate(grid_st%dy)
          deallocate(grid_st%dz)
          deallocate(grid_st%dxh)
          deallocate(grid_st%dyh)
          deallocate(grid_st%dzh)
          deallocate(grid_st%nxv)
          deallocate(grid_st%nyv)
          deallocate(grid_st%nzv)
          deallocate(grid_st%nxgl)
          deallocate(grid_st%nygl)
          deallocate(grid_st%nzgl)
          deallocate(grid_st%ilo)
          deallocate(grid_st%jlo)
          deallocate(grid_st%klo)
          deallocate(grid_st%ihi)
          deallocate(grid_st%jhi)
          deallocate(grid_st%khi)
cc          deallocate(grid_st%ntotv)
          deallocate(grid_st%istartx)
          deallocate(grid_st%istarty)
          deallocate(grid_st%istartz)
          deallocate(grid_st%istartp)
          deallocate(grid_st%mg_ratio_x)
          deallocate(grid_st%mg_ratio_y)
          deallocate(grid_st%mg_ratio_z)
          deallocate(grid_st%iline)
          deallocate(grid_st%jline)
          deallocate(grid_st%kline)
        endif

c     End program

      end subroutine deallocateGridStructure

c     equateGridStructure
c     #################################################################
      subroutine equateGridStructure(grid_st2,grid_st1)

c     -----------------------------------------------------------------
c     Performs grid_st2 = grid_st1
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),intent(IN ) :: grid_st1
        type(grid_mg_def),intent(OUT) :: grid_st2

c     Local variables

        integer(4)     :: ngridx,ngridy,ngridz,nx,ny,nz

c     Begin program

        ngridx = grid_st1%ngrdx
        ngridy = grid_st1%ngrdy
        ngridz = grid_st1%ngrdz

        nx = grid_st1%nxv(1)
        ny = grid_st1%nyv(1)
        nz = grid_st1%nzv(1)

        if (grid_st2%ngrdx /= ngridx) nullify(grid_st2%xx)

        call allocateGridStructure(nx,ny,nz,ngridx,ngridy,ngridz
     .                            ,grid_st2)

        grid_st2%xg         = grid_st1%xg        
        grid_st2%yg         = grid_st1%yg        
        grid_st2%zg         = grid_st1%zg        
        grid_st2%xx         = grid_st1%xx        
        grid_st2%yy         = grid_st1%yy        
        grid_st2%zz         = grid_st1%zz        
        grid_st2%dx         = grid_st1%dx        
        grid_st2%dy         = grid_st1%dy        
        grid_st2%dz         = grid_st1%dz        
        grid_st2%dxh        = grid_st1%dxh       
        grid_st2%dyh        = grid_st1%dyh       
        grid_st2%dzh        = grid_st1%dzh       
        grid_st2%nxv        = grid_st1%nxv       
        grid_st2%nyv        = grid_st1%nyv       
        grid_st2%nzv        = grid_st1%nzv       
        grid_st2%nxgl       = grid_st1%nxgl       
        grid_st2%nygl       = grid_st1%nygl       
        grid_st2%nzgl       = grid_st1%nzgl       
        grid_st2%ilo        = grid_st1%ilo
        grid_st2%jlo        = grid_st1%jlo
        grid_st2%klo        = grid_st1%klo
        grid_st2%ihi        = grid_st1%ihi
        grid_st2%jhi        = grid_st1%jhi
        grid_st2%khi        = grid_st1%khi
cc        grid_st2%ntotv      = grid_st1%ntotv     
        grid_st2%istartx    = grid_st1%istartx   
        grid_st2%istarty    = grid_st1%istarty   
        grid_st2%istartz    = grid_st1%istartz   
        grid_st2%istartp    = grid_st1%istartp   
        grid_st2%mg_ratio_x = grid_st1%mg_ratio_x
        grid_st2%mg_ratio_y = grid_st1%mg_ratio_y
        grid_st2%mg_ratio_z = grid_st1%mg_ratio_z
        grid_st2%iline      = grid_st1%iline
        grid_st2%jline      = grid_st1%jline
        grid_st2%kline      = grid_st1%kline

c     End program

      end subroutine equateGridStructure

c     writeGridStructure
c     #################################################################
      subroutine writeGridStructure(grid_st)

c     -----------------------------------------------------------------
c     Performs grid_st2 = grid_st1
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def) :: grid_st

c     Local variables

c     Begin program

        write (*,*) 'ngrdx',grid_st%ngrdx
        write (*,*) 'ngrdy',grid_st%ngrdy
        write (*,*) 'ngrdz',grid_st%ngrdz
        write (*,*) 'xg',grid_st%xg
        write (*,*) 'yg',grid_st%yg        
        write (*,*) 'zg',grid_st%zg        
        write (*,*) 'xx',grid_st%xx        
        write (*,*) 'yy',grid_st%yy        
        write (*,*) 'zz',grid_st%zz        
        write (*,*) 'dx',grid_st%dx        
        write (*,*) 'dy',grid_st%dy        
        write (*,*) 'dz',grid_st%dz        
        write (*,*) 'dxh',grid_st%dxh       
        write (*,*) 'dyh',grid_st%dyh       
        write (*,*) 'dzh',grid_st%dzh       
        write (*,*) 'nxv',grid_st%nxv       
        write (*,*) 'nyv',grid_st%nyv       
        write (*,*) 'nzv',grid_st%nzv
        write (*,*) 'nxgl',grid_st%nxgl       
        write (*,*) 'nygl',grid_st%nygl       
        write (*,*) 'nzgl',grid_st%nzgl
        write (*,*) 'ilo',grid_st%ilo       
        write (*,*) 'jlo',grid_st%jlo       
        write (*,*) 'klo',grid_st%klo
        write (*,*) 'ihi',grid_st%ihi       
        write (*,*) 'jhi',grid_st%jhi       
        write (*,*) 'khi',grid_st%khi       
cc        write (*,*) 'ntotv',grid_st%ntotv     
        write (*,*) 'istartx',grid_st%istartx   
        write (*,*) 'istarty',grid_st%istarty   
        write (*,*) 'istartz',grid_st%istartz   
        write (*,*) 'istartp',grid_st%istartp   
        write (*,*) 'mg_ratio_x',grid_st%mg_ratio_x
        write (*,*) 'mg_ratio_y',grid_st%mg_ratio_y
        write (*,*) 'mg_ratio_z',grid_st%mg_ratio_z
        write (*,*) 'iline',grid_st%iline
        write (*,*) 'jline',grid_st%jline
        write (*,*) 'kline',grid_st%kline

c     End program

      end subroutine writeGridStructure

c     getMGmap
c     #################################################################
      subroutine getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c     -----------------------------------------------------------------
c     Gets MG vector components (ig,jg,kg) for grid quantities
c     corresponding to node position (i,j,k) in grid levels igx,igy,igz
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg

c     Local variables

c     Begin program

        ig = i + grid_params%istartx(igx)
        jg = j + grid_params%istarty(igy)
        kg = k + grid_params%istartz(igz)

      end subroutine getMGmap

c     getCurvilinearCoordinates
c     #################################################################
      subroutine getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Finds curvilinear coordinates for position (i,j,k)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        integer(4) :: ii,jj,ny

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x1 = grid_params%xx(ig)
        y1 = grid_params%yy(jg)
        z1 = grid_params%zz(kg)

      end subroutine getCurvilinearCoordinates

c     vecPos
c     #################################################################
      function vecPos(neq,i,j,k,igx,igy,igz) result (ijk)

        implicit none

        integer(4) :: ijk,i,j,k,neq,igx,igy,igz

        integer(4) :: nxl,nyl,nzl,ilo,jlo,klo

        ilo = 1
        jlo = 1
        klo = 1

        nxl = grid_params%nxv(igx)
        nyl = grid_params%nyv(igy)
        nzl = grid_params%nzv(igz)

        ijk = neq*(i-ilo + nxl*(j-jlo) + nxl*nyl*(k-klo))
      
      end function vecPos

c     fromGlobalToLocalLimits
c     #################################################################
      subroutine fromGlobalToLocalLimits(ig,jg,kg,il,jl,kl,igx,igy,igz)

        implicit none

        integer(4) :: ig,jg,kg,il,jl,kl,igx,igy,igz

        il = ig - grid_params%ilo(igx) + 1
        jl = jg - grid_params%jlo(igy) + 1
        kl = kg - grid_params%klo(igz) + 1

      end subroutine fromGlobalToLocalLimits

c     fromLocalToGlobalLimits
c     #################################################################
      subroutine fromLocalToGlobalLimits(il,jl,kl,ig,jg,kg,igx,igy,igz)

        implicit none

        integer(4) :: ig,jg,kg,il,jl,kl,igx,igy,igz

        ig = il + grid_params%ilo(igx) - 1
        jg = jl + grid_params%jlo(igy) - 1
        kg = kl + grid_params%klo(igz) - 1

      end subroutine fromLocalToGlobalLimits

      end module grid_structures

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
c     Standard BCs:
c       EQU: Inhomogeneous dirichlet (generally imposed by equilibrium)
c       DIR: Homogeneous dirichlet
c       NEU: Neumann (homogeneous/inhomogeneous)
c       SYM: Symmetry (odd,even depend on sign)
c       EXT: By extrapolation
c       DEF: Default (does nothing)
c     The actual integers chosen indicate a particular order that the BC
c     scheduler follows to impose the boundary conditions.
c     -----------------------------------------------------------------

        use grid_structures

        implicit none

        integer(4) :: PER,DIR,NEU,SP,EQU,DEF,EXT,SYM,IFC
        parameter (DEF=0,PER=1,EQU=2,EXT=3,NEU=4,SYM=5,DIR=6,SP=7,IFC=8)

        integer(4),parameter :: BCLIM=8    !Defines # of BC cases

        integer(4) :: bcond(6)

        logical    :: asm_dir(3)=.false.,asm=.false.,par_bc=.false.

        real(8)    :: SP_flsv=1d-5         !SP failsafe
        
      contains

c     bcSP
c     #################################################################
      function bcSP(ibc) result(sing_bc)

        integer(4),optional :: ibc
        logical    :: sing_bc

        integer(4) :: iibc

        if (PRESENT(ibc)) then
          iibc = ibc
        else
          iibc = 1
        endif

        sing_bc = (bcond(iibc) == SP)

      end function bcSP

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        integer(4) :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point = (i+grid_params%ilo(igx)-1==1).and.bcSP()

      end function isSP

c     isSYM
c     #################################################################
      function isSYM(i,igr,dim,loc)

        integer(4) :: i,igr,dim,loc
        logical    :: isSYM

        integer(4) :: ibc,ilog,itst

        ibc = (1+loc)+2*(dim-1)

        isSYM=.false.

        select case(dim)
        case(1)
          ilog = grid_params%ilo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nxgl(igr)
          end select
        case(2)
          ilog = grid_params%jlo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nygl(igr)
          end select
        case(3)
          ilog = grid_params%klo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nzgl(igr)
          end select
        end select

        if (i+ilog-1 == itst .and. bcond(ibc) == SYM)  isSYM=.true.

      end function isSYM

c     alt__eom
c     #################################################################
      function alt__eom()

        logical :: alt__eom

        alt__eom=.false.
        if (bcSP()) alt__eom=.true.

      end function alt__eom

      end module bc_def

c module grid_metric_def
c #####################################################################
      module grid_metric_def

c ---------------------------------------------------------------------
c     This module allocates and fills the grid metric structure gmetric.
c     It contains the following routines: 
c        * allocateGridMetric
c        * deallocateGridMetric
c        * defineGridMetric
c           -> map
c           -> topol_bc
c           -> gammaZeroForce
c           -> zeroForceCheck
c ---------------------------------------------------------------------

        use grid_structures

        implicit none

        type :: grid_metrics
          real(8),pointer,dimension(:,:,:,:)     :: car   !Node positions in cartesian space
          real(8),pointer,dimension(:,:,:)       :: dvol  !Local cell volume
          real(8),pointer,dimension(:,:,:)       :: jac   !Jacobian factor at grid cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsub  !Covariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsup  !Contravariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:,:) :: Gamma !Christoffel symbol at cell centers
          real(8),pointer,dimension(:,:,:,:,:)   :: cov   !Covariant vectors
          real(8),pointer,dimension(:,:,:,:,:)   :: cnv   !Contravariant vectors
        end type grid_metrics

        type :: MG_grid_metrics
          type(grid_metrics),pointer,dimension(:) :: grid
        end type MG_grid_metrics

        type(MG_grid_metrics) :: gmetric

      contains

c     allocateGridMetric
c     #################################################################
      subroutine allocateGridMetric(gmetric)

        implicit none

c     Call variables

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4)      :: igrid,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        if (.not.associated(gmetric%grid)) then
          allocate(gmetric%grid(grid_params%ngrid))
          do igrid=1,grid_params%ngrid
            ilom = 0
            jlom = 0
            klom = 0
            ihip = grid_params%nxv(igrid)+1
            jhip = grid_params%nyv(igrid)+1
            khip = grid_params%nzv(igrid)+1
            allocate(gmetric%grid(igrid)
     $           %car  (ilom:ihip,jlom:jhip,klom:khip,3))
            allocate(gmetric%grid(igrid)
     $           %dvol (ilom:ihip,jlom:jhip,klom:khip))
            allocate(gmetric%grid(igrid)
     $           %jac  (ilom:ihip,jlom:jhip,klom:khip))
            allocate(gmetric%grid(igrid)
     $           %gsub (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %gsup (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %cov  (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %cnv  (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %Gamma(ilom:ihip,jlom:jhip,klom:khip,3,3,3))

          enddo
        endif

c     End program

      end subroutine allocateGridMetric

c     deallocateGridMetric
c     #################################################################
      subroutine deallocateGridMetric(gmetric)

        implicit none

c     Call variables

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4)      :: igrid,nxp,nyp,nzp

c     Begin program

        if (associated(gmetric%grid)) then
          do igrid=1,grid_params%ngrid
            deallocate(gmetric%grid(igrid)%car  )
            deallocate(gmetric%grid(igrid)%dvol )
            deallocate(gmetric%grid(igrid)%jac  )
            deallocate(gmetric%grid(igrid)%gsub )
            deallocate(gmetric%grid(igrid)%gsup )
            deallocate(gmetric%grid(igrid)%cov  )
            deallocate(gmetric%grid(igrid)%cnv  )
            deallocate(gmetric%grid(igrid)%Gamma)
          enddo
          deallocate(gmetric%grid)
        endif

c     End program

      end subroutine deallocateGridMetric

c     volume
c     #################################################################
      function volume(i,j,k,igx,igy,igz) result(vol)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: vol

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,dx1,dx2,dx3,jac

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dx1 = grid_params%dxh(ig)
        dx2 = grid_params%dyh(jg)
        dx3 = grid_params%dzh(kg)

        jac = gmetric%grid(igx)%jac(i,j,k)
        
        vol = jac*dx1*dx2*dx3

      end function volume

c     getCartesianCoordinates
c     #################################################################
      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Returns Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        real(8)    :: car(3)

c     Begin program

        !Find ig,jg,kg
        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        !Find coordinates
        x1 = gmetric%grid(igx)%car(i,j,k,1)
        y1 = gmetric%grid(igx)%car(i,j,k,2)
        z1 = gmetric%grid(igx)%car(i,j,k,3)

      end subroutine getCartesianCoordinates

      end module grid_metric_def

c module grid_def
c ####################################################################
      module grid_def

        use bc_def

        use grid_metric_def

        real(8)         :: xmax,ymax,zmax,xmin,ymin,zmin  !3D domain dimensions

        character*(3)   :: coords

      end module grid_def

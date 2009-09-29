module grid_module
  implicit none

  type :: var_def
     integer :: bconds(6)                            !Boundary conditions
     real(8),pointer,dimension(:,:,:) :: array
     character(20) :: descr
  end type var_def

  type :: vec_def
     integer :: bconds(6,3)                          !Boundary conditions
     real(8),pointer,dimension(:,:,:,:) :: vec
     character(20) :: descr                          !Description
     logical       :: cnv                            !Whether it is contravariant or not
  end type vec_def

  type :: var_array
     integer :: nvar                                 !Number of variables
     type (var_def),pointer,dimension(:) :: array_var
  end type var_array

  type :: aux_array
     integer :: nvar,nvec                            !Number of variables, vectors
     type (vec_def),pointer,dimension(:) :: vec_list
     type (var_def),pointer,dimension(:) :: var_list
  end type aux_array

  integer   ,parameter :: ngparams = 6,UNDEF_INT=-1234

  type :: dim_pack
     logical    :: pack
     real(8)    :: xp
     real(8)    :: dx0
  end type dim_pack

  type :: grid_pack
     type(dim_pack),dimension(3) :: dim
  end type grid_pack

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

  type :: grid_mg_def
     integer :: ngrdx                                !Grid levels in X
     integer :: ngrdy                                !Grid levels in Y
     integer :: ngrdz                                !Grid levels in Z
     integer :: ngrid                                !Max grid levels for MG

     integer :: nglx                                 !Global number of mesh points in X
     integer :: ngly                                 !Global number of mesh points in Y
     integer :: nglz                                 !Global number of mesh points in Z

     integer :: nlx                                  !Local number of mesh points in X
     integer :: nly                                  !Local number of mesh points in Y
     integer :: nlz                                  !Local number of mesh points in Z

     integer :: ilog                                 !Global lower  limit in X
     integer :: ihig                                 !Global higher limit in X
     integer :: jlog                                 !Global lower  limit in Y
     integer :: jhig                                 !Global higher limit in Y
     integer :: klog                                 !Global lower  limit in Z
     integer :: khig                                 !Global higher limit in Z

     integer :: gcw                                  !Ghost cell width

     real(8) :: lxmin,lxmax                          !Local patch physical dimensions in X
     real(8) :: lymin,lymax                          !Local patch physical dimensions in Y
     real(8) :: lzmin,lzmax                          !Local patch physical dimensions in Z

     integer   ,pointer,dimension(:)  :: ilo         !Global lower limit in X
     integer   ,pointer,dimension(:)  :: jlo         !Global lower limit in Y
     integer   ,pointer,dimension(:)  :: klo         !Global lower limit in Z
     integer   ,pointer,dimension(:)  :: ihi         !Global higher limit in X
     integer   ,pointer,dimension(:)  :: jhi         !Global higher limit in Y
     integer   ,pointer,dimension(:)  :: khi         !Global higher limit in Z

     real(8)   ,pointer,dimension(:)  :: xg          !Global grid nodes in X (top grid only)
     real(8)   ,pointer,dimension(:)  :: yg          !Global grid nodes in Y (")
     real(8)   ,pointer,dimension(:)  :: zg          !Global grid nodes in Z (")
     real(8)   ,pointer,dimension(:)  :: xx          !Local grid nodes in X (all grids)
     real(8)   ,pointer,dimension(:)  :: yy          !Local grid nodes in Y (")
     real(8)   ,pointer,dimension(:)  :: zz          !Local grid nodes in Z (")

     real(8)   ,pointer,dimension(:)  :: dx          !Grid spacings in X for integer mesh (all grids)
     real(8)   ,pointer,dimension(:)  :: dy          !Grid spacings in Y for integer mesh (")
     real(8)   ,pointer,dimension(:)  :: dz          !Grid spacings in Z for integer mesh (")
     real(8)   ,pointer,dimension(:)  :: dxh         !Grid spacings in X for half mesh (")
     real(8)   ,pointer,dimension(:)  :: dyh         !Grid spacings in Y for half mesh (")
     real(8)   ,pointer,dimension(:)  :: dzh         !Grid spacings in Z for half mesh (")

     integer   ,pointer,dimension(:)  :: nxv         !Local # of grid nodes in X  (all grids)
     integer   ,pointer,dimension(:)  :: nyv         !Local # of grid nodes in Y  (")
     integer   ,pointer,dimension(:)  :: nzv         !Local # of grid nodes in Z  (")
     integer   ,pointer,dimension(:)  :: nxgl        !Global # of grid nodes in X  (")
     integer   ,pointer,dimension(:)  :: nygl        !Global # of grid nodes in Y  (")
     integer   ,pointer,dimension(:)  :: nzgl        !Global # of grid nodes in Z  (")

     integer   ,pointer,dimension(:)  :: istartx     !Pointer for MG vectors in X
     integer   ,pointer,dimension(:)  :: istarty     !Pointer for MG vectors in Y
     integer   ,pointer,dimension(:)  :: istartz     !Pointer for MG vectors in Z
     integer   ,pointer,dimension(:)  :: istartp     !Pointer for global MG vectors
     integer   ,pointer,dimension(:)  :: mg_ratio_x  !MG coarsening ratio in X
     integer   ,pointer,dimension(:)  :: mg_ratio_y  !MG coarsening ratio in Y
     integer   ,pointer,dimension(:)  :: mg_ratio_z  !MG coarsening ratio in Z
     integer   ,pointer,dimension(:)  :: iline       !Restrict ops. to i=iline in MG
     integer   ,pointer,dimension(:)  :: jline       !Restrict ops. to j=jline in MG
     integer   ,pointer,dimension(:)  :: kline       !Restrict ops. to k=kline in MG

     real(8)                          :: params(ngparams)   !Grid configuration parameters

     type(MG_grid_metrics),pointer    :: gmetric     !Grid metrics at all levels

     type(grid_pack),pointer          :: g_pack      !Grid packing information

  end type grid_mg_def

  type :: patch
     integer  :: nbc_seq
     integer,pointer,dimension(:,:):: bc_seq
     type(var_array)  ,pointer     :: u_0
     type(var_array)  ,pointer     :: u_n
     type(var_array)  ,pointer     :: u_np
     type(var_array)  ,pointer     :: u_nm
     type(var_array)  ,pointer     :: u_graph
     type(grid_mg_def),pointer     :: gparams
     type(aux_array)  ,pointer     :: aux
  end type patch

end module grid_module

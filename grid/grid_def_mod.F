
c module grid_def_st
c #####################################################################
      module grid_def_st

        use grid_petsc

        use io, ONLY: sstop

        implicit none

        integer   ,parameter :: ngparams = 50,UNDEF_INT=-1234
        real(8)   ,parameter :: UNDEF_DP=-1234d0

        type :: grid_pack
          sequence
          logical    :: pack = .false.
          real(8)    :: xp=0d0
          real(8)    :: dx0=0d0
          logical    :: pack_both_bdrys=.false.
          integer    :: modes=1
          real(8)    :: unif_dom_left_bdry=0d0
          real(8)    :: unif_dom_rite_bdry=0d0
          integer    :: unif_nx=0
          !Put these last
          real(8)    :: slope=0d0
          real(8)    :: dxl=0d0
          real(8)    :: dxr=0d0
          logical    :: processed=.false.
          integer    :: left_nx=0
          integer    :: rite_nx=0
        end type grid_pack

        type :: grid_metrics
#if defined(samrai)
          sequence
#endif
          real(8),pointer,dimension(:,:,:,:)     :: car   => null()!Node positions in cartesian space
          real(8),pointer,dimension(:,:,:)       :: dvol  => null()!Local cell physical volume
          real(8),pointer,dimension(:,:,:)       :: dlvol => null()!Local cell logical volume
          real(8),pointer,dimension(:,:,:)       :: jac   => null()!Jacobian factor at grid cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsub  => null()!Covariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsup  => null()!Contravariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:,:) :: Gamma => null()!Christoffel symbol at cell centers
          real(8),pointer,dimension(:,:,:,:,:)   :: cov   => null()!Covariant vectors
          real(8),pointer,dimension(:,:,:,:,:)   :: cnv   => null()!Contravariant vectors

          real(8),pointer,dimension(:,:,:,:)     :: car_v => null()!Vertex positions in cartesian space
          real(8),pointer,dimension(:,:,:)       :: jac_v => null()!Jacobian factor at grid vertices
          real(8),pointer,dimension(:,:,:,:,:)   :: gsub_v=> null()!Covariant metric tensor at vertices
          real(8),pointer,dimension(:,:,:,:,:)   :: gsup_v=> null()!Contravariant metric tensor at vertices
          real(8),pointer,dimension(:,:,:,:,:)   :: cov_v => null()!Covariant vectors at vertices
          real(8),pointer,dimension(:,:,:,:,:)   :: cnv_v => null()!Contravariant vectors at vertices
        end type grid_metrics

        type :: MG_grid_metrics
#if defined(samrai)
          sequence 
#endif
          type(grid_metrics),pointer,dimension(:) :: grid => null()
        end type MG_grid_metrics

        type :: grid_mg_def
#if defined(samrai)
          sequence 
!!          Variables below have changed from earlier PIXIE3D version!!!
#endif
          !!!!!!!!!External grid manager (e.g. SAMRAI) needs to define these!!!!!!!!
          integer :: nglx    !Global number of mesh points in X, over global domain, not local patch
          integer :: ngly    !Global number of mesh points in Y, over global domain, not local patch
          integer :: nglz    !Global number of mesh points in Z, over global domain, not local patch

          integer :: ilog    !Global lower  limit in X, patch indices relative to global index space, no ghosts
          integer :: ihig    !Global higher limit in X, patch indices relative to global index space
          integer :: jlog    !Global lower  limit in Y, patch indices relative to global index space
          integer :: jhig    !Global higher limit in Y, patch indices relative to global index space
          integer :: klog    !Global lower  limit in Z, patch indices relative to global index space
          integer :: khig    !Global higher limit in Z, patch indices relative to global index space

          real(8) :: gxmin,gxmax        !Physical dim. of global domain in X
          real(8) :: gymin,gymax        !Physical dim. of global domain in Y
          real(8) :: gzmin,gzmax        !Physical dim. of global domain in Z

          real(8) :: lxmin,lxmax        !Physical dim. of LOCAL domain in X
          real(8) :: lymin,lymax        !Physical dim. of LOCAL domain in Y
          real(8) :: lzmin,lzmax        !Physical dim. of LOCAL domain in Z

          !MPI info
          integer :: np,rank,mx,my,mz

          integer :: aligner            !Ghost cell width
          integer :: gcw                !Ghost cell width
          !!!!!!!!!External grid manager (e.g. SAMRAI) needs to define these!!!!!!!!
#if defined(petsc)
          integer :: mpi_comm       =MPI_COMM_NULL  !Current communicator
          !Various group communicators
          integer :: MPI_COMM_SP    =MPI_COMM_NULL
     .              ,MPI_COMM_POL   =MPI_COMM_NULL
     .              ,MPI_COMM_RAD   =MPI_COMM_NULL
     .              ,MPI_COMM_X     =MPI_COMM_NULL
     .              ,MPI_COMM_Y     =MPI_COMM_NULL
     .              ,MPI_COMM_Z     =MPI_COMM_NULL
     .              ,MPI_COMM_XY    =MPI_COMM_NULL
     .              ,MPI_COMM_XZ    =MPI_COMM_NULL
     .              ,MPI_COMM_YZ    =MPI_COMM_NULL
     .              ,MPI_COMM_PER(3)=MPI_COMM_NULL
     .              ,MPI_COMM_CRSE  =MPI_COMM_NULL  !Groups procs for coarse MG solve
#endif
          integer :: ngrdx                                         !Grid levels in X, MG specific
          integer :: ngrdy                                         !Grid levels in Y, MG specific
          integer :: ngrdz                                         !Grid levels in Z, MG specific
          integer :: ngrid                                         !Max grid levels for MG

          integer :: nlx                                           !Local number of mesh points in X
          integer :: nly                                           !Local number of mesh points in Y
          integer :: nlz                                           !Local number of mesh points in Z

          integer   ,pointer,dimension(:)  :: ilo => null()        !Global lower limit in X in all grids
          integer   ,pointer,dimension(:)  :: jlo => null()        !Global lower limit in Y in all grids
          integer   ,pointer,dimension(:)  :: klo => null()        !Global lower limit in Z in all grids
          integer   ,pointer,dimension(:)  :: ihi => null()        !Global higher limit in X in all grids
          integer   ,pointer,dimension(:)  :: jhi => null()        !Global higher limit in Y in all grids
          integer   ,pointer,dimension(:)  :: khi => null()        !Global higher limit in Z in all grids

          real(8)   ,pointer,dimension(:)  :: xg => null()         !Global grid nodes in X (top grid only), x co-ordinates of global domain useful for a packed tensor grid
          real(8)   ,pointer,dimension(:)  :: yg => null()         !Global grid nodes in Y (")
          real(8)   ,pointer,dimension(:)  :: zg => null()         !Global grid nodes in Z (")
          real(8)   ,pointer,dimension(:)  :: xx => null()         !Local grid nodes in X (all grids)
          real(8)   ,pointer,dimension(:)  :: yy => null()         !Local grid nodes in Y (")
          real(8)   ,pointer,dimension(:)  :: zz => null()         !Local grid nodes in Z (")

          real(8)   ,pointer,dimension(:)  :: dx => null()         !Grid spacings in X for integer mesh (all grids)
          real(8)   ,pointer,dimension(:)  :: dy => null()         !Grid spacings in Y for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dz => null()         !Grid spacings in Z for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dxh=> null()         !Grid spacings in X for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dyh=> null()         !Grid spacings in Y for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dzh=> null()         !Grid spacings in Z for half mesh (")

          integer   ,pointer,dimension(:)  :: nxv => null()        !Local # of grid nodes in X  (all MG grids)
          integer   ,pointer,dimension(:)  :: nyv => null()        !Local # of grid nodes in Y  (")
          integer   ,pointer,dimension(:)  :: nzv => null()        !Local # of grid nodes in Z  (")
          integer   ,pointer,dimension(:)  :: nxgl=> null()        !Global # of grid nodes in X  (")
          integer   ,pointer,dimension(:)  :: nygl=> null()        !Global # of grid nodes in Y  (")
          integer   ,pointer,dimension(:)  :: nzgl=> null()        !Global # of grid nodes in Z  (")

          integer   ,pointer,dimension(:)  :: istartx   => null()  !Pointer for MG vectors in X
          integer   ,pointer,dimension(:)  :: istarty   => null()  !Pointer for MG vectors in Y
          integer   ,pointer,dimension(:)  :: istartz   => null()  !Pointer for MG vectors in Z
          integer   ,pointer,dimension(:)  :: istartp   => null()  !Pointer for global MG vectors
          integer   ,pointer,dimension(:)  :: mg_ratio_x=> null()  !MG coarsening ratio in X
          integer   ,pointer,dimension(:)  :: mg_ratio_y=> null()  !MG coarsening ratio in Y
          integer   ,pointer,dimension(:)  :: mg_ratio_z=> null()  !MG coarsening ratio in Z
          integer   ,pointer,dimension(:)  :: iline     => null()  !Restrict ops. to i=iline in MG
          integer   ,pointer,dimension(:)  :: jline     => null()  !Restrict ops. to j=jline in MG
          integer   ,pointer,dimension(:)  :: kline     => null()  !Restrict ops. to k=kline in MG

          real(8)                          :: params(ngparams)   !Grid configuration parameters

          type(MG_grid_metrics),pointer    :: gmetric => null()  !Grid metrics at all levels

          type(grid_pack),dimension(3)     :: g_pack      !Grid packing information

#if defined(petsc) && !defined(samrai)
          type(petsc_da_ctx),pointer,dimension(:) :: dactx => null()  !PETSc DA context
#endif

          type(grid_mg_def),pointer        :: g_crse_def => null()

        end type grid_mg_def

      contains

c     allocateGridStructure
c     #################################################################
      subroutine allocateGridStructure(ngridx,ngridy,ngridz,mx,my,mz
     .                                ,nglx,ngly,nglz,grid_st,mpi_comm)
c     -----------------------------------------------------------------
c     Allocates grid structure
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: ngridx,ngridy,ngridz,nglx,ngly,nglz,mx,my,mz
        type(grid_mg_def),pointer :: grid_st
        integer,optional :: mpi_comm

c     Local variables

        integer :: nx,ny,nz,ngrid,nxmg,nymg,nzmg
     .            ,ilog,ihig,jlog,jhig,klog,khig,mpierr

c     Begin program

        ngrid = max(ngridx,ngridy,ngridz)

        if (.not.associated(grid_st)) then
          allocate(grid_st)
          nullify(grid_st%xx)
          nullify(grid_st%gmetric)

          grid_st%nglx = nglx
          grid_st%ngly = ngly
          grid_st%nglz = nglz

          grid_st%gxmin = UNDEF_DP    !Indicates it is provided by input file
          grid_st%gxmax = UNDEF_DP
          grid_st%gymin = UNDEF_DP
          grid_st%gymax = UNDEF_DP
          grid_st%gzmin = UNDEF_DP
          grid_st%gzmax = UNDEF_DP

          grid_st%lxmin = UNDEF_DP    !Indicates it is provided by input file
          grid_st%lxmax = UNDEF_DP
          grid_st%lymin = UNDEF_DP
          grid_st%lymax = UNDEF_DP
          grid_st%lzmin = UNDEF_DP
          grid_st%lzmax = UNDEF_DP

#if defined(petsc) && !defined(samrai)
          if (PRESENT(mpi_comm)) then
            grid_st%mpi_comm = mpi_comm
          else
            grid_st%mpi_comm = MPI_COMM_NULL
          endif

          grid_st%mx = mx
          grid_st%my = my
          grid_st%mz = mz

          call MPI_Comm_rank(grid_st%mpi_comm,grid_st%rank,mpierr)
          call MPI_Comm_size(grid_st%mpi_comm,grid_st%np  ,mpierr)

          if (mx*my*mz /= grid_st%np) then
            call sstop(grid_st%rank,"allocateGridStructure"
     .          ,"Processor distribution does not match total # procs")
          endif

          allocate(grid_st%dactx(ngrid))

          call createPETScGrid(grid_st,1)

          grid_st%ilog  = grid_st%dactx(1)%xs
          grid_st%ihig  = grid_st%dactx(1)%xe
          grid_st%jlog  = grid_st%dactx(1)%ys
          grid_st%jhig  = grid_st%dactx(1)%ye
          grid_st%klog  = grid_st%dactx(1)%zs
          grid_st%khig  = grid_st%dactx(1)%ze
#else
          grid_st%ilog  = 1
          grid_st%ihig  = nglx
          grid_st%jlog  = 1
          grid_st%jhig  = ngly
          grid_st%klog  = 1
          grid_st%khig  = nglz

          grid_st%mx = 1
          grid_st%my = 1
          grid_st%mz = 1

          grid_st%np    = 1
          grid_st%rank  = 0
#endif
        endif

        grid_st%ngrdx = ngridx
        grid_st%ngrdy = ngridy
        grid_st%ngrdz = ngridz
        grid_st%ngrid = ngrid

        ilog = grid_st%ilog
        ihig = grid_st%ihig
        jlog = grid_st%jlog
        jhig = grid_st%jhig
        klog = grid_st%klog
        khig = grid_st%khig

        grid_st%nlx = ihig-ilog+1
        grid_st%nly = jhig-jlog+1
        grid_st%nlz = khig-klog+1

        nx = grid_st%nlx
        ny = grid_st%nly
        nz = grid_st%nlz

        if (.not.associated(grid_st%xx)) then
          nxmg=findMGsize(nx,grid_st%ngrdx,grid_st%ngrid)
          nymg=findMGsize(ny,grid_st%ngrdy,grid_st%ngrid)
          nzmg=findMGsize(nz,grid_st%ngrdz,grid_st%ngrid)

          allocate(grid_st%xg(0:grid_st%nglx+1))
          allocate(grid_st%yg(0:grid_st%ngly+1))
          allocate(grid_st%zg(0:grid_st%nglz+1))
          allocate(grid_st%xx (nxmg))
          allocate(grid_st%yy (nymg))
          allocate(grid_st%zz (nzmg))
          allocate(grid_st%dx (nxmg))
          allocate(grid_st%dy (nymg))
          allocate(grid_st%dz (nzmg))
          allocate(grid_st%dxh(nxmg))
          allocate(grid_st%dyh(nymg))
          allocate(grid_st%dzh(nzmg))
          allocate(grid_st%nxv (grid_st%ngrid))
          allocate(grid_st%nyv (grid_st%ngrid))
          allocate(grid_st%nzv (grid_st%ngrid))
          allocate(grid_st%nxgl(grid_st%ngrid))
          allocate(grid_st%nygl(grid_st%ngrid))
          allocate(grid_st%nzgl(grid_st%ngrid))
          allocate(grid_st%ilo (grid_st%ngrid))
          allocate(grid_st%jlo (grid_st%ngrid))
          allocate(grid_st%klo (grid_st%ngrid))
          allocate(grid_st%ihi (grid_st%ngrid))
          allocate(grid_st%jhi (grid_st%ngrid))
          allocate(grid_st%khi (grid_st%ngrid))
          allocate(grid_st%istartx(grid_st%ngrid))
          allocate(grid_st%istarty(grid_st%ngrid))
          allocate(grid_st%istartz(grid_st%ngrid))
          allocate(grid_st%istartp(grid_st%ngrid))
          allocate(grid_st%mg_ratio_x(grid_st%ngrid))
          allocate(grid_st%mg_ratio_y(grid_st%ngrid))
          allocate(grid_st%mg_ratio_z(grid_st%ngrid))
          allocate(grid_st%iline(grid_st%ngrid))
          allocate(grid_st%jline(grid_st%ngrid))
          allocate(grid_st%kline(grid_st%ngrid))

          !Initialize arrays
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
          grid_st%gcw = 1     !(ONLY gcw=1 works right now!)
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
          grid_st%params = 0d0

          grid_st%g_pack(1)%pack = .false.
          grid_st%g_pack(2)%pack = .false.
          grid_st%g_pack(3)%pack = .false.
        endif

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
c     enough space even if nn=1 in the finest grid. It includes
c     ghost cells (which requires an additional term of 2*ngrdt).
c     -----------------------------------------------------------------

        integer    :: nn,nnmg,ngrd,ngrdt

cc        nnmg = (2+ngrdt-ngrd)*nn
        nnmg = (2+ngrdt-ngrd)*nn+2*ngrdt

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

        type(grid_mg_def),pointer :: grid_st

c     Local variables

        integer :: igr

#if defined(petsc) && !defined(samrai)
        type(petsc_da_ctx),pointer :: da_ctx => null()
#endif

c     Begin program
        
        if (.not.associated(grid_st)) return

c$$$c     Deallocate coarser levels
c$$$
c$$$        call deallocateGridStructure(grid_st%g_crse_def)

c     Deallocate gmetric

        call deallocateGridMetric(grid_st%gmetric)

c     Deallocate PETSC DA context

#if defined(petsc)  && !defined(samrai)
        if (associated(grid_st%dactx)) then
          do igr=1,size(grid_st%dactx)
            da_ctx => grid_st%dactx(igr)
            call destroyPETScFortranDA(da_ctx)
          enddo

          deallocate(grid_st%dactx)
        endif
#endif

c     Deallocate MG storage variables

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
        
        deallocate(grid_st)

c     End program

      end subroutine deallocateGridStructure

c$$$c     equateGridStructure
c$$$c     #################################################################
c$$$      subroutine equateGridStructure(grid_st2,grid_st1)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Performs grid_st2 = grid_st1, where MG storage arrays are 
c$$$c     strictly copied, PETSC DA are regenerated, and grid metrics
c$$$c     are simply pointed to.
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(grid_mg_def),pointer :: grid_st1
c$$$        type(grid_mg_def),pointer :: grid_st2
c$$$
c$$$c     Local variables
c$$$
c$$$        integer :: ngridx,ngridy,ngridz,nx,ny,nz
c$$$     .            ,nglx,ngly,nglz,mx,my,mz
c$$$     .            ,ilog,ihig,jlog,jhig,klog,khig
c$$$
c$$$c     Begin program
c$$$
c$$$        ngridx = grid_st1%ngrdx
c$$$        ngridy = grid_st1%ngrdy
c$$$        ngridz = grid_st1%ngrdz
c$$$
c$$$        nx = grid_st1%nxv(1)
c$$$        ny = grid_st1%nyv(1)
c$$$        nz = grid_st1%nzv(1)
c$$$
c$$$        nglx = grid_st1%nglx
c$$$        ngly = grid_st1%ngly
c$$$        nglz = grid_st1%nglz
c$$$
c$$$        mx = grid_st1%mx
c$$$        my = grid_st1%my
c$$$        mz = grid_st1%mz
c$$$
c$$$        ilog = grid_st1%ilog
c$$$        ihig = grid_st1%ihig
c$$$        jlog = grid_st1%jlog
c$$$        jhig = grid_st1%jhig
c$$$        klog = grid_st1%klog
c$$$        khig = grid_st1%khig
c$$$
c$$$        if (associated(grid_st2)) then
c$$$          if (    grid_st2%ngrdx /= ngridx
c$$$     .        .or.grid_st2%ngrdy /= ngridy
c$$$     .        .or.grid_st2%ngrdz /= ngridz
c$$$     .        .or.grid_st2%nxv(1)/= nx
c$$$     .        .or.grid_st2%nyv(1)/= ny
c$$$     .        .or.grid_st2%nzv(1)/= nz
c$$$     .        .or.grid_st2%nglx  /= nglx
c$$$     .        .or.grid_st2%ngly  /= ngly
c$$$     .        .or.grid_st2%nglz  /= nglz)
c$$$     .      call deallocateGridStructure(grid_st2)
c$$$        endif
c$$$
c$$$c       Create basic structure
c$$$
c$$$        call allocateGridStructure(ngridx,ngridy,ngridz,mx,my,mz
c$$$     .                            ,nglx,ngly,nglz,grid_st2)
c$$$
c$$$c       Copy values
c$$$
c$$$        grid_st2%xg         = grid_st1%xg        
c$$$        grid_st2%yg         = grid_st1%yg        
c$$$        grid_st2%zg         = grid_st1%zg        
c$$$        grid_st2%xx         = grid_st1%xx        
c$$$        grid_st2%yy         = grid_st1%yy        
c$$$        grid_st2%zz         = grid_st1%zz        
c$$$        grid_st2%dx         = grid_st1%dx        
c$$$        grid_st2%dy         = grid_st1%dy        
c$$$        grid_st2%dz         = grid_st1%dz        
c$$$        grid_st2%dxh        = grid_st1%dxh       
c$$$        grid_st2%dyh        = grid_st1%dyh       
c$$$        grid_st2%dzh        = grid_st1%dzh       
c$$$        grid_st2%nxv        = grid_st1%nxv       
c$$$        grid_st2%nyv        = grid_st1%nyv       
c$$$        grid_st2%nzv        = grid_st1%nzv       
c$$$        grid_st2%nxgl       = grid_st1%nxgl       
c$$$        grid_st2%nygl       = grid_st1%nygl       
c$$$        grid_st2%nzgl       = grid_st1%nzgl       
c$$$        grid_st2%ilo        = grid_st1%ilo
c$$$        grid_st2%jlo        = grid_st1%jlo
c$$$        grid_st2%klo        = grid_st1%klo
c$$$        grid_st2%ihi        = grid_st1%ihi
c$$$        grid_st2%jhi        = grid_st1%jhi
c$$$        grid_st2%khi        = grid_st1%khi
c$$$        grid_st2%gcw        = grid_st1%gcw
c$$$        grid_st2%istartx    = grid_st1%istartx   
c$$$        grid_st2%istarty    = grid_st1%istarty   
c$$$        grid_st2%istartz    = grid_st1%istartz   
c$$$        grid_st2%istartp    = grid_st1%istartp   
c$$$        grid_st2%mg_ratio_x = grid_st1%mg_ratio_x
c$$$        grid_st2%mg_ratio_y = grid_st1%mg_ratio_y
c$$$        grid_st2%mg_ratio_z = grid_st1%mg_ratio_z
c$$$        grid_st2%iline      = grid_st1%iline
c$$$        grid_st2%jline      = grid_st1%jline
c$$$        grid_st2%kline      = grid_st1%kline
c$$$
c$$$        !Grid info
c$$$        grid_st2%gxmin      = grid_st1%gxmin
c$$$        grid_st2%gxmax      = grid_st1%gxmax
c$$$        grid_st2%gymin      = grid_st1%gymin
c$$$        grid_st2%gymax      = grid_st1%gymax
c$$$        grid_st2%gzmin      = grid_st1%gzmin
c$$$        grid_st2%gzmax      = grid_st1%gzmax
c$$$        grid_st2%lxmin      = grid_st1%lxmin
c$$$        grid_st2%lxmax      = grid_st1%lxmax
c$$$        grid_st2%lymin      = grid_st1%lymin
c$$$        grid_st2%lymax      = grid_st1%lymax
c$$$        grid_st2%lzmin      = grid_st1%lzmin
c$$$        grid_st2%lzmax      = grid_st1%lzmax
c$$$
c$$$        grid_st2%g_pack     = grid_st1%g_pack
c$$$        grid_st2%params     = grid_st1%params
c$$$
c$$$        !MPI info
c$$$#if defined(petsc)
c$$$        grid_st2%mpi_comm   = grid_st1%mpi_comm
c$$$        grid_st2%np         = grid_st1%np
c$$$        grid_st2%rank       = grid_st1%rank
c$$$        grid_st2%mx         = grid_st1%mx
c$$$        grid_st2%my         = grid_st1%my
c$$$        grid_st2%mz         = grid_st1%mz
c$$$        grid_st2%MPI_COMM_SP  = grid_st1%MPI_COMM_SP    
c$$$        grid_st2%MPI_COMM_POL = grid_st1%MPI_COMM_POL   
c$$$        grid_st2%MPI_COMM_RAD = grid_st1%MPI_COMM_RAD   
c$$$        grid_st2%MPI_COMM_X   = grid_st1%MPI_COMM_X     
c$$$        grid_st2%MPI_COMM_Y   = grid_st1%MPI_COMM_Y     
c$$$        grid_st2%MPI_COMM_Z   = grid_st1%MPI_COMM_Z     
c$$$        grid_st2%MPI_COMM_XY  = grid_st1%MPI_COMM_XY    
c$$$        grid_st2%MPI_COMM_XZ  = grid_st1%MPI_COMM_XZ    
c$$$        grid_st2%MPI_COMM_YZ  = grid_st1%MPI_COMM_YZ    
c$$$        grid_st2%MPI_COMM_PER = grid_st2%MPI_COMM_PER
c$$$        grid_st2%MPI_COMM_CRSE= grid_st2%MPI_COMM_CRSE
c$$$#endif
c$$$
c$$$c       Create and copy grid metric
c$$$
c$$$        if (associated(grid_st1%gmetric)) then
c$$$          call equateGridMetric(grid_st2%gmetric
c$$$     .                         ,grid_st1%gmetric,grid_st1)
c$$$        else
c$$$          nullify(grid_st2%gmetric)
c$$$        endif
c$$$
c$$$c       Create and copy coarse grid metric
c$$$
c$$$        if (associated(grid_st1%g_crse_def)) then
c$$$          grid_st2%g_crse_def => grid_st1%g_crse_def
c$$$        else
c$$$          nullify(grid_st2%g_crse_def)
c$$$        endif
c$$$
c$$$c       Create fortran PETSc parallel grid hierarchy (parallel version)
c$$$
c$$$#if defined(petsc)
c$$$        call createPETScGrid(grid_st2,0)
c$$$#endif
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine equateGridStructure

c     writeGridStructure
c     #################################################################
      subroutine writeGridStructure(grid_st)

c     -----------------------------------------------------------------
c     Dumps info on grid structure.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: grid_st

c     Local variables

c     Begin program

        write (*,*) 'ngrdx',grid_st%ngrdx
        write (*,*) 'ngrdy',grid_st%ngrdy
        write (*,*) 'ngrdz',grid_st%ngrdz
        write (*,*) 'nglx',grid_st%nglx
        write (*,*) 'ngly',grid_st%ngly
        write (*,*) 'nglz',grid_st%nglz
        write (*,*) 'mx',grid_st%mx
        write (*,*) 'my',grid_st%my
        write (*,*) 'mz',grid_st%mz
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
        write (*,*) 'gcw',grid_st%gcw
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

c     allocateGridMetric
c     #################################################################
      subroutine allocateGridMetric(gmetric,gdef)

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gdef
        type(MG_grid_metrics),pointer :: gmetric

c     Local variables

        integer :: igrid,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        if (.not.associated(gmetric)) then
          allocate(gmetric)
          allocate(gmetric%grid(gdef%ngrid))
          do igrid=1,gdef%ngrid
            ilom = 0
            jlom = 0
            klom = 0
            ihip = gdef%nxv(igrid)+1
            jhip = gdef%nyv(igrid)+1
            khip = gdef%nzv(igrid)+1
            allocate(gmetric%grid(igrid)
     $           %car  (ilom:ihip,jlom:jhip,klom:khip,3))
            allocate(gmetric%grid(igrid)
     $           %dvol (ilom:ihip,jlom:jhip,klom:khip))
            allocate(gmetric%grid(igrid)
     $           %dlvol(ilom:ihip,jlom:jhip,klom:khip))
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

            allocate(gmetric%grid(igrid)
     $           %car_v (ilom:ihip,jlom:jhip,klom:khip,3))
            allocate(gmetric%grid(igrid)
     $           %jac_v (ilom:ihip,jlom:jhip,klom:khip))
            allocate(gmetric%grid(igrid)
     $           %gsub_v(ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %gsup_v(ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %cov_v (ilom:ihip,jlom:jhip,klom:khip,3,3))
            allocate(gmetric%grid(igrid)
     $           %cnv_v (ilom:ihip,jlom:jhip,klom:khip,3,3))

            gmetric%grid(igrid)%car   = 0d0
            gmetric%grid(igrid)%dvol  = 0d0
            gmetric%grid(igrid)%dlvol = 0d0
            gmetric%grid(igrid)%jac   = 0d0
            gmetric%grid(igrid)%gsub  = 0d0
            gmetric%grid(igrid)%gsup  = 0d0
            gmetric%grid(igrid)%cov   = 0d0
            gmetric%grid(igrid)%cnv   = 0d0
            gmetric%grid(igrid)%Gamma = 0d0
            gmetric%grid(igrid)%car_v = 0d0
            gmetric%grid(igrid)%jac_v = 0d0
            gmetric%grid(igrid)%gsub_v= 0d0
            gmetric%grid(igrid)%gsup_v= 0d0
            gmetric%grid(igrid)%cov_v = 0d0
            gmetric%grid(igrid)%cnv_v = 0d0
          enddo
        endif

c     End program

      end subroutine allocateGridMetric

c     equateGridMetric
c     #################################################################
      subroutine equateGridMetric(gmetric2,gmetric1,gdef)

c     -----------------------------------------------------------------
c     Performs gmetric2 = gmetric1
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(MG_grid_metrics),pointer :: gmetric2,gmetric1
        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer :: igrid,nxp,nyp,nzp

c     Begin program

        if (.not.associated(gmetric2))
     .       call allocateGridMetric(gmetric2,gdef)

        do igrid=1,gdef%ngrid
          gmetric2%grid(igrid)%car  = gmetric1%grid(igrid)%car  
          gmetric2%grid(igrid)%dvol = gmetric1%grid(igrid)%dvol 
          gmetric2%grid(igrid)%dlvol= gmetric1%grid(igrid)%dlvol 
          gmetric2%grid(igrid)%jac  = gmetric1%grid(igrid)%jac  
          gmetric2%grid(igrid)%gsub = gmetric1%grid(igrid)%gsub 
          gmetric2%grid(igrid)%gsup = gmetric1%grid(igrid)%gsup 
          gmetric2%grid(igrid)%cov  = gmetric1%grid(igrid)%cov  
          gmetric2%grid(igrid)%cnv  = gmetric1%grid(igrid)%cnv  
          gmetric2%grid(igrid)%Gamma= gmetric1%grid(igrid)%Gamma

          gmetric2%grid(igrid)%car_v  = gmetric1%grid(igrid)%car_v
          gmetric2%grid(igrid)%jac_v  = gmetric1%grid(igrid)%jac_v
          gmetric2%grid(igrid)%gsub_v = gmetric1%grid(igrid)%gsub_v 
          gmetric2%grid(igrid)%gsup_v = gmetric1%grid(igrid)%gsup_v 
          gmetric2%grid(igrid)%cov_v  = gmetric1%grid(igrid)%cov_v  
          gmetric2%grid(igrid)%cnv_v  = gmetric1%grid(igrid)%cnv_v  
        enddo

c     End program

      end subroutine equateGridMetric

c     deallocateGridMetric
c     #################################################################
      subroutine deallocateGridMetric(gmetric)

        implicit none

c     Call variables

        type(MG_grid_metrics),pointer :: gmetric

c     Local variables

        integer         :: igrid,nxp,nyp,nzp

c     Begin program

        if (associated(gmetric)) then
          do igrid=1,size(gmetric%grid)
            deallocate(gmetric%grid(igrid)%car)
            deallocate(gmetric%grid(igrid)%dvol)
            deallocate(gmetric%grid(igrid)%dlvol)
            deallocate(gmetric%grid(igrid)%jac)
            deallocate(gmetric%grid(igrid)%gsub)
            deallocate(gmetric%grid(igrid)%gsup)
            deallocate(gmetric%grid(igrid)%cov)
            deallocate(gmetric%grid(igrid)%cnv)
            deallocate(gmetric%grid(igrid)%Gamma)

            deallocate(gmetric%grid(igrid)%car_v)
            deallocate(gmetric%grid(igrid)%jac_v)
            deallocate(gmetric%grid(igrid)%gsub_v)
            deallocate(gmetric%grid(igrid)%gsup_v)
            deallocate(gmetric%grid(igrid)%cov_v)
            deallocate(gmetric%grid(igrid)%cnv_v)
          enddo
          deallocate(gmetric%grid)
          deallocate(gmetric)
        endif

c     End program

      end subroutine deallocateGridMetric

#if defined(petsc) && !defined(samrai)

c     createPETScGrid
c     #################################################################
      subroutine createPETScGrid(gdef,igr)

c     -----------------------------------------------------------------
c     Creates PETSc constructs for boundary communication in fortran.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igr
        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer :: nxgl,nygl,nzgl,igrid,mpi_comm
        type(petsc_da_ctx),pointer :: da_ctx => null()

c     Begin program

        mpi_comm = gdef%mpi_comm

c     Create PETSc DA MG hierachy

        select case(igr)
        case(1)

          da_ctx => gdef%dactx(igr)

          da_ctx%MPI_COMM_DA = mpi_comm

          da_ctx%igrid = igr

          da_ctx%mx = gdef%mx
          da_ctx%my = gdef%my
          da_ctx%mz = gdef%mz

          da_ctx%nxgl = gdef%nglx
          da_ctx%nygl = gdef%ngly
          da_ctx%nzgl = gdef%nglz

          call createPETScFortranDA(da_ctx)

        case(0)

          do igrid=2,gdef%ngrid

            da_ctx => gdef%dactx(igrid)

            da_ctx%MPI_COMM_DA = mpi_comm

            da_ctx%igrid = igrid

            da_ctx%mx = gdef%mx
            da_ctx%my = gdef%my
            da_ctx%mz = gdef%mz

            da_ctx%nxgl = gdef%nxgl(igrid)
            da_ctx%nygl = gdef%nygl(igrid)
            da_ctx%nzgl = gdef%nzgl(igrid)

            call createPETScFortranDA(da_ctx)

          enddo

        case default

          da_ctx => gdef%dactx(igr)

          da_ctx%MPI_COMM_DA = mpi_comm

          da_ctx%igrid = igr

          da_ctx%mx = gdef%mx
          da_ctx%my = gdef%my
          da_ctx%mz = gdef%mz

          da_ctx%nxgl = gdef%nxgl(igr)
          da_ctx%nygl = gdef%nygl(igr)
          da_ctx%nzgl = gdef%nzgl(igr)

          call createPETScFortranDA(da_ctx)

        end select

        nullify(da_ctx)
        
c     End program

      end subroutine createPETScGrid
#endif

      end module grid_def_st

c module grid_def
c ####################################################################
      module grid_def

        use bc_def

        use grid_def_st

        character*(3) :: coords

      end module grid_def

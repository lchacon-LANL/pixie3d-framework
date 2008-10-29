
c module grid_def_st
c #####################################################################
      module grid_def_st

        implicit none

        integer   ,parameter :: ngparams = 6,UNDEF_INT=-1234
        real(8)   ,parameter :: UNDEF_DP=-1234.

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
          !SAMRAI needs to define these!
          integer :: ngrdx                                !Grid levels in X
          integer :: ngrdy                                !Grid levels in Y
          integer :: ngrdz                                !Grid levels in Z
          integer :: ngrid                                !Max grid levels for MG

          integer :: nglx                                 !Global number of mesh points in X
          integer :: ngly                                 !Global number of mesh points in Y
          integer :: nglz                                 !Global number of mesh points in Z

          integer :: ilog                                 !Global lower  limit in X
          integer :: ihig                                 !Global higher limit in X
          integer :: jlog                                 !Global lower  limit in Y
          integer :: jhig                                 !Global higher limit in Y
          integer :: klog                                 !Global lower  limit in Z
          integer :: khig                                 !Global higher limit in Z
          !SAMRAI needs to define these!

          integer :: nlx                                  !Local number of mesh points in X
          integer :: nly                                  !Local number of mesh points in Y
          integer :: nlz                                  !Local number of mesh points in Z

          integer :: gcw                                  !Ghost cell width

          real(8) :: lxmin,lxmax                          !Global physical dimensions in X
          real(8) :: lymin,lymax                          !Global physical dimensions in Y
          real(8) :: lzmin,lzmax                          !Global physical dimensions in Z

          integer   ,pointer,dimension(:)  :: ilo         !Global lower limit in X in all grids
          integer   ,pointer,dimension(:)  :: jlo         !Global lower limit in Y in all grids
          integer   ,pointer,dimension(:)  :: klo         !Global lower limit in Z in all grids
          integer   ,pointer,dimension(:)  :: ihi         !Global higher limit in X in all grids
          integer   ,pointer,dimension(:)  :: jhi         !Global higher limit in Y in all grids
          integer   ,pointer,dimension(:)  :: khi         !Global higher limit in Z in all grids

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

        type(grid_mg_def),pointer :: grid_params

        type(MG_grid_metrics),pointer :: gmetric

      contains

c     allocateGridStructure
c     #################################################################
      subroutine allocateGridStructure(ngridx,ngridy,ngridz
     .                                ,ilog,ihig,jlog,jhig,klog,khig
     .                                ,nglx,ngly,nglz,grid_st)
c     -----------------------------------------------------------------
c     Allocates grid structure
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: ngridx,ngridy,ngridz,nglx,ngly,nglz
     .            ,ilog,ihig,jlog,jhig,klog,khig
        type(grid_mg_def),pointer :: grid_st

c     Local variables

        integer :: nx,ny,nz,ngrid,nxmg,nymg,nzmg

c     Begin program

        ngrid = max(ngridx,ngridy,ngridz)

        if (.not.associated(grid_st)) then  !Take grid defs from input if existing.
          allocate(grid_st)
          nullify(grid_st%xx)
          nullify(grid_st%gmetric)

          grid_st%ngrdx = ngridx
          grid_st%ngrdy = ngridy
          grid_st%ngrdz = ngridz
          grid_st%ngrid = ngrid

          grid_st%nglx = nglx
          grid_st%ngly = ngly
          grid_st%nglz = nglz

          grid_st%ilog  = ilog
          grid_st%ihig  = ihig
          grid_st%jlog  = jlog
          grid_st%jhig  = jhig
          grid_st%klog  = klog
          grid_st%khig  = khig
        endif

        grid_st%lxmin = UNDEF_DP
        grid_st%lxmax = UNDEF_DP
        grid_st%lymin = UNDEF_DP
        grid_st%lymax = UNDEF_DP
        grid_st%lzmin = UNDEF_DP
        grid_st%lzmax = UNDEF_DP

        grid_st%nlx  = ihig-ilog+1
        grid_st%nly  = jhig-jlog+1
        grid_st%nlz  = khig-klog+1

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
          allocate(grid_st%xx(nxmg+2*grid_st%ngrid))
          allocate(grid_st%yy(nymg+2*grid_st%ngrid))
          allocate(grid_st%zz(nzmg+2*grid_st%ngrid))
          allocate(grid_st%dx(nxmg+2*grid_st%ngrid))
          allocate(grid_st%dy(nymg+2*grid_st%ngrid))
          allocate(grid_st%dz(nzmg+2*grid_st%ngrid))
          allocate(grid_st%dxh(nxmg+2*grid_st%ngrid))
          allocate(grid_st%dyh(nymg+2*grid_st%ngrid))
          allocate(grid_st%dzh(nzmg+2*grid_st%ngrid))
          allocate(grid_st%nxv(grid_st%ngrid))
          allocate(grid_st%nyv(grid_st%ngrid))
          allocate(grid_st%nzv(grid_st%ngrid))
          allocate(grid_st%nxgl(grid_st%ngrid))
          allocate(grid_st%nygl(grid_st%ngrid))
          allocate(grid_st%nzgl(grid_st%ngrid))
          allocate(grid_st%ilo(grid_st%ngrid))
          allocate(grid_st%jlo(grid_st%ngrid))
          allocate(grid_st%klo(grid_st%ngrid))
          allocate(grid_st%ihi(grid_st%ngrid))
          allocate(grid_st%jhi(grid_st%ngrid))
          allocate(grid_st%khi(grid_st%ngrid))
cc          allocate(grid_st%ntotv(grid_st%ngrid))
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

          allocate(grid_st%g_pack)

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
c     enough space even if nn=1 in the finest grid. It does NOT include
c     ghost cells (this requires an additional term of 2*ngrdt).
c     -----------------------------------------------------------------

        integer    :: nn,nnmg,ngrd,ngrdt

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

        type(grid_mg_def),pointer :: grid_st

c     Begin program

        call deallocateGridMetric(grid_st%gmetric)

        deallocate(grid_st%g_pack)

        grid_st%ngrdx = 0
        grid_st%ngrdy = 0
        grid_st%ngrdz = 0
        grid_st%ngrid = 0

        grid_st%nglx = 0
        grid_st%ngly = 0
        grid_st%nglz = 0

        grid_st%ilog = 0
        grid_st%ihig = 0
        grid_st%jlog = 0
        grid_st%jhig = 0
        grid_st%klog = 0
        grid_st%khig = 0

        grid_st%lxmin = UNDEF_DP
        grid_st%lxmax = UNDEF_DP
        grid_st%lymin = UNDEF_DP
        grid_st%lymax = UNDEF_DP
        grid_st%lzmin = UNDEF_DP
        grid_st%lzmax = UNDEF_DP

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

        deallocate(grid_st)

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

        type(grid_mg_def),pointer :: grid_st1
        type(grid_mg_def),pointer :: grid_st2

c     Local variables

        integer        :: ngridx,ngridy,ngridz,nx,ny,nz
     .                   ,nglx,ngly,nglz
     .                   ,ilog,ihig,jlog,jhig,klog,khig

c     Begin program

        ngridx = grid_st1%ngrdx
        ngridy = grid_st1%ngrdy
        ngridz = grid_st1%ngrdz

        nx = grid_st1%nxv(1)
        ny = grid_st1%nyv(1)
        nz = grid_st1%nzv(1)

        nglx = grid_st1%nglx
        ngly = grid_st1%ngly
        nglz = grid_st1%nglz

        ilog = grid_st1%ilog
        ihig = grid_st1%ihig
        jlog = grid_st1%jlog
        jhig = grid_st1%jhig
        klog = grid_st1%klog
        khig = grid_st1%khig

        if (associated(grid_st2)) then
          if (    grid_st2%ngrdx /= ngridx
     .        .or.grid_st2%ngrdy /= ngridy
     .        .or.grid_st2%ngrdz /= ngridz
     .        .or.grid_st2%nxv(1)/= nx
     .        .or.grid_st2%nyv(1)/= ny
     .        .or.grid_st2%nzv(1)/= nz
     .        .or.grid_st2%nglx  /= nglx
     .        .or.grid_st2%ngly  /= ngly
     .        .or.grid_st2%nglz  /= nglz)
     .      call deallocateGridStructure(grid_st2)
        endif

c       Create basic structure

        call allocateGridStructure(ngridx,ngridy,ngridz
     .                            ,ilog,ihig,jlog,jhig,klog,khig
     .                            ,nglx,ngly,nglz,grid_st2)

c       Copy values

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
        grid_st2%gcw        = grid_st1%gcw
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
        grid_st2%lxmin      = grid_st1%lxmin
        grid_st2%lxmax      = grid_st1%lxmax
        grid_st2%lymin      = grid_st1%lymin
        grid_st2%lymax      = grid_st1%lymax
        grid_st2%lzmin      = grid_st1%lzmin
        grid_st2%lzmax      = grid_st1%lzmax

        grid_st2%g_pack     = grid_st1%g_pack

c       Create and copy grid metric

        call allocateGridMetric(grid_st2%gmetric,grid_st1)

        call equateGridMetric(grid_st2%gmetric
     .                       ,grid_st1%gmetric,grid_st1)

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

        type(grid_mg_def),pointer :: grid_st

c     Local variables

c     Begin program

        write (*,*) 'ngrdx',grid_st%ngrdx
        write (*,*) 'ngrdy',grid_st%ngrdy
        write (*,*) 'ngrdz',grid_st%ngrdz
        write (*,*) 'nglx',grid_st%nglx
        write (*,*) 'ngly',grid_st%ngly
        write (*,*) 'nglz',grid_st%nglz
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
      subroutine allocateGridMetric(gmetric,gparams)

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gparams
        type(MG_grid_metrics),pointer :: gmetric

c     Local variables

        integer :: igrid,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        if (.not.associated(gmetric)) then
          allocate(gmetric)
          allocate(gmetric%grid(gparams%ngrid))
          do igrid=1,gparams%ngrid
            ilom = 0
            jlom = 0
            klom = 0
            ihip = gparams%nxv(igrid)+1
            jhip = gparams%nyv(igrid)+1
            khip = gparams%nzv(igrid)+1
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

c     equateGridMetric
c     #################################################################
      subroutine equateGridMetric(gmetric2,gmetric1,gparams)

c     -----------------------------------------------------------------
c     Performs gmetric2 = gmetric1
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(MG_grid_metrics),pointer :: gmetric2,gmetric1
        type(grid_mg_def),pointer :: gparams

c     Local variables

        integer :: igrid,nxp,nyp,nzp

c     Begin program

        if (.not.associated(gmetric2))
     .       call allocateGridMetric(gmetric2,gparams)

        do igrid=1,gparams%ngrid
          gmetric2%grid(igrid)%car  = gmetric1%grid(igrid)%car  
          gmetric2%grid(igrid)%dvol = gmetric1%grid(igrid)%dvol 
          gmetric2%grid(igrid)%jac  = gmetric1%grid(igrid)%jac  
          gmetric2%grid(igrid)%gsub = gmetric1%grid(igrid)%gsub 
          gmetric2%grid(igrid)%gsup = gmetric1%grid(igrid)%gsup 
          gmetric2%grid(igrid)%cov  = gmetric1%grid(igrid)%cov  
          gmetric2%grid(igrid)%cnv  = gmetric1%grid(igrid)%cnv  
          gmetric2%grid(igrid)%Gamma= gmetric1%grid(igrid)%Gamma
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
            deallocate(gmetric%grid(igrid)%jac)
            deallocate(gmetric%grid(igrid)%gsub)
            deallocate(gmetric%grid(igrid)%gsup)
            deallocate(gmetric%grid(igrid)%cov)
            deallocate(gmetric%grid(igrid)%cnv)
            deallocate(gmetric%grid(igrid)%Gamma)
          enddo
          deallocate(gmetric%grid)
          deallocate(gmetric)
        endif

c     End program

      end subroutine deallocateGridMetric

      end module grid_def_st

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

        use grid_def_st

        implicit none

        integer    :: PER,DIR,NEU,SP,EQU,DEF,EXT,SYM,IFC
        parameter (DEF=0,PER=1,EQU=2,EXT=3,NEU=4,SYM=5,DIR=6,SP=7,IFC=8)

        integer   ,parameter :: BCLIM=8    !Defines # of BC cases

        integer    :: bcond(6)

        logical    :: asm_dir(3)=.false.,asm=.false.,par_bc=.false.

        real(8)    :: SP_flsv=1d-5         !SP failsafe
        
      contains

c     bcSP
c     #################################################################
      function bcSP(ibc) result(sing_bc)

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

c     alt__eom
c     #################################################################
      function alt__eom()

        logical :: alt__eom

        alt__eom=.false.
        if (bcSP()) alt__eom=.true.

      end function alt__eom

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        integer    :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point = (i+grid_params%ilo(igx)-1==1).and.bcSP()

      end function isSP

c     isSYM
c     #################################################################
      function isSYM(i,igr,dim,loc,ibc)

        integer    :: i,igr
        logical    :: isSYM
        integer,optional :: dim,loc,ibc

        integer    :: ibcl,ilog,itst,diml,locl

        if (PRESENT(dim) .and.PRESENT(loc)) then
          diml = dim
          locl = loc
          ibcl = (1+locl)+2*(diml-1)
        elseif (PRESENT(ibc)) then
          ibcl = ibc
          locl = abs(mod(ibcl,2)-2)-1
          diml = 1 + (ibcl-1-locl)/2
        else
          write (*,*) 'Error in isSYM'
          write (*,*) 'Aborting...'
          stop
        endif

        isSYM=.false.

        select case(diml)
        case(1)
          ilog = grid_params%ilo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nxgl(igr)
          end select
        case(2)
          ilog = grid_params%jlo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nygl(igr)
          end select
        case(3)
          ilog = grid_params%klo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nzgl(igr)
          end select
        end select

        if (i+ilog-1 == itst .and. bcond(ibcl) == SYM)  isSYM=.true.

      end function isSYM

      end module bc_def

ccc module grid_metric_def
ccc #####################################################################
cc      module grid_metric_def
cc
ccc ---------------------------------------------------------------------
ccc     This module allocates and fills the grid metric structure gmetric.
ccc     It contains the following routines: 
ccc        * allocateGridMetric
ccc        * deallocateGridMetric
ccc        * defineGridMetric
ccc           -> map
ccc           -> topol_bc
ccc           -> gammaZeroForce
ccc           -> zeroForceCheck
ccc ---------------------------------------------------------------------
cc
cc        use grid_structures
cc
cc        implicit none
cc
cc        type :: grid_metrics
cc          real(8),pointer,dimension(:,:,:,:)     :: car   !Node positions in cartesian space
cc          real(8),pointer,dimension(:,:,:)       :: dvol  !Local cell volume
cc          real(8),pointer,dimension(:,:,:)       :: jac   !Jacobian factor at grid cells
cc          real(8),pointer,dimension(:,:,:,:,:)   :: gsub  !Covariant metric tensor at cells
cc          real(8),pointer,dimension(:,:,:,:,:)   :: gsup  !Contravariant metric tensor at cells
cc          real(8),pointer,dimension(:,:,:,:,:,:) :: Gamma !Christoffel symbol at cell centers
cc          real(8),pointer,dimension(:,:,:,:,:)   :: cov   !Covariant vectors
cc          real(8),pointer,dimension(:,:,:,:,:)   :: cnv   !Contravariant vectors
cc        end type grid_metrics
cc
cc        type :: MG_grid_metrics
cc          type(grid_metrics),pointer,dimension(:) :: grid
cc        end type MG_grid_metrics
cc
cccc        type(MG_grid_metrics),pointer :: gmetric
cc
cc      contains
cc
ccc     allocateGridMetric
ccc     #################################################################
cc      subroutine allocateGridMetric(gmetric)
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        type(MG_grid_metrics),pointer :: gmetric
cc
ccc     Local variables
cc
cc        integer :: igrid,ilom,ihip,jlom,jhip,klom,khip
cc
ccc     Begin program
cc
cc        if (.not.associated(gmetric)) then
cc          allocate(gmetric)
cc          allocate(gmetric%grid(grid_params%ngrid))
cc          do igrid=1,grid_params%ngrid
cc            ilom = 0
cc            jlom = 0
cc            klom = 0
cc            ihip = grid_params%nxv(igrid)+1
cc            jhip = grid_params%nyv(igrid)+1
cc            khip = grid_params%nzv(igrid)+1
cc            allocate(gmetric%grid(igrid)
cc     $           %car  (ilom:ihip,jlom:jhip,klom:khip,3))
cc            allocate(gmetric%grid(igrid)
cc     $           %dvol (ilom:ihip,jlom:jhip,klom:khip))
cc            allocate(gmetric%grid(igrid)
cc     $           %jac  (ilom:ihip,jlom:jhip,klom:khip))
cc            allocate(gmetric%grid(igrid)
cc     $           %gsub (ilom:ihip,jlom:jhip,klom:khip,3,3))
cc            allocate(gmetric%grid(igrid)
cc     $           %gsup (ilom:ihip,jlom:jhip,klom:khip,3,3))
cc            allocate(gmetric%grid(igrid)
cc     $           %cov  (ilom:ihip,jlom:jhip,klom:khip,3,3))
cc            allocate(gmetric%grid(igrid)
cc     $           %cnv  (ilom:ihip,jlom:jhip,klom:khip,3,3))
cc            allocate(gmetric%grid(igrid)
cc     $           %Gamma(ilom:ihip,jlom:jhip,klom:khip,3,3,3))
cc
cc          enddo
cc        endif
cc
ccc     End program
cc
cc      end subroutine allocateGridMetric
cc
ccc     deallocateGridMetric
ccc     #################################################################
cc      subroutine deallocateGridMetric(gmetric)
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        type(MG_grid_metrics),pointer :: gmetric
cc
ccc     Local variables
cc
cc        integer         :: igrid,nxp,nyp,nzp
cc
ccc     Begin program
cc
cc        deallocate(gmetric)
cc
cccc        if (associated(gmetric%grid)) then
cccc          do igrid=1,grid_params%ngrid
cccc            deallocate(gmetric%grid(igrid)%car  )
cccc            deallocate(gmetric%grid(igrid)%dvol )
cccc            deallocate(gmetric%grid(igrid)%jac  )
cccc            deallocate(gmetric%grid(igrid)%gsub )
cccc            deallocate(gmetric%grid(igrid)%gsup )
cccc            deallocate(gmetric%grid(igrid)%cov  )
cccc            deallocate(gmetric%grid(igrid)%cnv  )
cccc            deallocate(gmetric%grid(igrid)%Gamma)
cccc          enddo
cccc          deallocate(gmetric%grid)
cccc        endif
cc
ccc     End program
cc
cc      end subroutine deallocateGridMetric
cc
ccccc     volume
ccccc     #################################################################
cccc      function volume(i,j,k,igx,igy,igz) result(vol)
cccc
ccccc     -----------------------------------------------------------------
ccccc     Calculates Jacobian of curvilinear coordinate system
ccccc     -----------------------------------------------------------------
cccc
cccc        implicit none
cccc
ccccc     Input variables
cccc
cccc        integer    :: i,j,k,igx,igy,igz
cccc        real(8)    :: vol
cccc
ccccc     Local variables
cccc
cccc        integer    :: ig,jg,kg
cccc        real(8)    :: x1,x2,x3,dx1,dx2,dx3,jac
cccc
ccccc     Begin program
cccc
cccc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc        dx1 = grid_params%dxh(ig)
cccc        dx2 = grid_params%dyh(jg)
cccc        dx3 = grid_params%dzh(kg)
cccc
cccc        jac = gmetric%grid(igx)%jac(i,j,k)
cccc        
cccc        vol = jac*dx1*dx2*dx3
cccc
cccc      end function volume
cc
ccccc     getCartesianCoordinates
ccccc     #################################################################
cccc      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                  ,x1,y1,z1)
cccc
ccccc     -----------------------------------------------------------------
ccccc     Returns Cartesian coordinates.
ccccc     -----------------------------------------------------------------
cccc
cccc        implicit none
cccc
ccccc     Input variables
cccc
cccc        integer    :: i,j,k,igx,igy,igz,ig,jg,kg
cccc        real(8)    :: x1,y1,z1
cccc
ccccc     Local variables
cccc
cccc        real(8)    :: car(3)
cccc
ccccc     Begin program
cccc
cccc        !Find coordinates
cccc        x1 = gmetric%grid(igx)%car(i,j,k,1)
cccc        y1 = gmetric%grid(igy)%car(i,j,k,2)
cccc        z1 = gmetric%grid(igz)%car(i,j,k,3)
cccc
cccc      end subroutine getCartesianCoordinates
cc
cc      end module grid_metric_def

c module grid_def
c ####################################################################
      module grid_def

        use bc_def

        use grid_def_st

        real(8)         :: xmax,ymax,zmax,xmin,ymin,zmin  !3D global domain dimensions

        character*(3)   :: coords

      end module grid_def

c module grid_structures
c #####################################################################
      module grid_structures

        implicit none

        type :: grid_def
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
          real(8)   ,pointer,dimension(:)  :: xx          !Grid node positions in X (all grids)
          real(8)   ,pointer,dimension(:)  :: yy          !Grid node positions in Y (")
          real(8)   ,pointer,dimension(:)  :: zz          !Grid node positions in Z (")
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
          real(8)                          :: params(5)   !Grid configuration parameters
        end type grid_def

        type (grid_def) :: grid_params

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
        type(grid_def) :: grid_st

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

        nnmg = 2*nn + ngrdt - ngrd -1

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

        type(grid_def) :: grid_st

c     Begin program

        grid_st%ngrdx = 0
        grid_st%ngrdy = 0
        grid_st%ngrdz = 0
        grid_st%ngrid = 0

        if (associated(grid_st%xx)) then
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

        type(grid_def),intent(IN ) :: grid_st1
        type(grid_def),intent(OUT) :: grid_st2

c     Local variables

        integer(4)     :: ngridx,ngridy,ngridz,nx,ny,nz

c     Begin program

        ngridx = grid_st1%ngrdx
        ngridy = grid_st1%ngrdy
        ngridz = grid_st1%ngrdz

        nx = grid_st1%nxv(1)
        ny = grid_st1%nyv(1)
        nz = grid_st1%nzv(1)

        call allocateGridStructure(nx,ny,nz,ngridx,ngridy,ngridz
     .                            ,grid_st2)

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

        type(grid_def) :: grid_st

c     Local variables

c     Begin program

        write (*,*) 'ngrdx',grid_st%ngrdx
        write (*,*) 'ngrdy',grid_st%ngrdy
        write (*,*) 'ngrdz',grid_st%ngrdz
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
      subroutine getMGmap(i,j,igr,ig,jg)

c     -----------------------------------------------------------------
c     Gets MG vector components (ig,jg,kg) for grid quantities
c     corresponding to node position (i,j,k) in grid levels igx,igy,igz
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,igr,ig,jg

c     Local variables

        integer(4) :: ligr

c     Begin program

c     Transform MG grid order from igr=1 finest to igr=1 coarsest

        ligr = grid_params%ngrid - igr + 1

c     Give MG vector coordinates

        ig = i + grid_params%istartx(ligr)
        jg = j + grid_params%istarty(ligr)

      end subroutine getMGmap

c     getCoordinates
c     #################################################################
      subroutine getCoordinates(i,j,igr,ig,jg,x1,y1)

c     -----------------------------------------------------------------
c     Finds curvilinear coordinates for position (i,j,k)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,igr,ig,jg
        real(8)    :: x1,y1

c     Local variables

c     Begin program

        call getMGmap(i,j,igr,ig,jg)

        x1 = grid_params%xx(ig)
        y1 = grid_params%yy(jg)

      end subroutine getCoordinates

      end module grid_structures

c module grid_definition
c #####################################################################
      module grid_definition

        use grid_structures

        implicit none

        real(8)         :: gparams(5)

        character*(3)   :: coords

        real(8)         :: xmax,ymax,zmax,xmin,ymin,zmin  !3D domain dimensions

        real(8),private :: pi,lambda,cc,ypp,eps,mm,kk,aa,phi,major_r

        logical         :: numerical_grid,anal_map

cc      contains
cc
ccc     checkGridDatabase
ccc     #################################################################
cc      function checkGridDatabase() result(anal_map)
cc
ccc     -----------------------------------------------------------------
ccc     Checks grid database for analytical mappings
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        logical :: anal_map
cc
ccc     Begin program
cc
cc        select case(coords)
cc        case('car','scl','cyl','hel','tor','sin')
cc
cc          anal_map = .true.
cc
cc        case default
cc
cc          anal_map = .false.
cc
cc        end select
cc
cc      end function checkGridDatabase
cc
ccc     getCartesianCoordinates
ccc     #################################################################
cc      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,y1,z1)
cc
ccc     -----------------------------------------------------------------
ccc     Inverts curvilinear coordinates to give Cartesian coordinates.
ccc     Requires external routine 'map', with call sequence:
ccc
ccc             map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1)
ccc
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
cc        real(8)    :: x1,y1,z1
cc
ccc     Local variables
cc
cc        real(8)    :: car(3)
cc
ccc     Externals
cc
cc        external   map
cc
ccc     Begin program
cc
cc        if (checkGridDatabase()) then
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,y1,z1)
cc
cc          car = x_xi(x1,y1,z1)
cc
cc          x1 = car(1)
cc          y1 = car(2)
cc          z1 = car(3)
cc
cc        else
cc
cc          call map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1) !External
cc
cc        endif
cc
cc      end subroutine getCartesianCoordinates
cc
ccc     x_xi
ccc     #################################################################
cc      function x_xi(x1,x2,x3) result(car)
cc
ccc     -----------------------------------------------------------------
ccc     Gives Cartesian coordinates from curvilinear coordinates
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        real(8)    :: x1,x2,x3,car(3)
cc
ccc     Local variables
cc
cc        integer(4) :: inewt,ic
cc        real(8)    :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          xx = x1
cc          yy = x2
cc          zz = x3
cc        case ('scl')
cc          lambda = gparams(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*x2/ymax-1.)
cc
cc          xx = x1
cc          yy = 0.5+lambda*atanh(ypp/cc)
cc          zz = x3
cc        case ('cyl')
cc          xx = x1*cos(x2)
cc          yy = x1*sin(x2)
cc          zz = x3
cc        case ('hel')
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm
cc          phi = (x2-aa*x3)
cc
cc          xx = x1*cos(phi)
cc          yy = x1*sin(phi)
cc          zz = x3
cc        case ('tor')
cc          major_r = gparams(1)
cc
cc          xx = (major_r + x1*sin(x2))*cos(x3)
cc          yy = (major_r + x1*sin(x2))*sin(x3)
cc          zz = x1*cos(x2)
cc        case ('sin')
cc          pi = acos(-1d0)
cc          eps = gparams(1)
cc
cc          xx = x1 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
cc          yy = x2 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
cc          zz = x3
cc        case default
cc          write (*,*) 'Grid not implemented in x_xi'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc        car = (/ xx,yy,zz /)
cc
cc      contains
cc
ccc     atanh
ccc     #################################################################
cc      real(8) function atanh(x)
cc
cc        real(8) :: x
cc
cc        atanh = 0.5*(log( (1+x)/(1-x) ) )
cc
cc      end function atanh
cc
cc      end function x_xi
cc
ccc     jacobian
ccc     #################################################################
cc      function jacobian(i,j,k,igx,igy,igz) result(jac)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates Jacobian of curvilinear coordinate system
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        real(8)    :: jac
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: x1,x2,x3,car(3),curv(3)
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          jac = 1d0
cc        case ('scl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gparams(1)
cc
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc        case ('cyl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          jac = x1
cc        case ('hel')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          jac = x1
cc        case ('tor')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          major_r = gparams(1)
cc
cc          jac = x1*(major_r + x1*sin(x2))
cc        case ('sin')
cc          pi = acos(-1d0)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          eps = gparams(1)
cc          jac = (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -          Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -          eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     -          /(xmax*ymax)
cc        case default
cc          write (*,*) 'Grid not implemented in jacobian'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc      end function jacobian
cc
ccc     covariantVector2
ccc     #################################################################
cc      function covariantVector(comp,i,j,k,igx,igy,igz) result (vec)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates covariant vectors of curvilinear coordinate system
ccc     in Cartesian coordinates
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: comp,i,j,k,igx,igy,igz
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: x1,x2,x3,vec(3)
cc
cc        real(8)    :: car(3),curv(3),jac
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          select case (comp)
cc            case (1)
cc              vec = (/ 1d0,0d0,0d0 /)
cc            case (2)
cc              vec = (/ 0d0,1d0,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
cc        case ('scl')
cc          lambda = gparams(1)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          curv = (/ x1,x2,x3 /)
cc
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          jac = 1./jac
cc          select case (comp)
cc            case (1)
cc              vec = (/ 1d0,0d0,0d0 /)
cc            case (2)
cc              vec = (/ 0d0,jac,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
cc        case ('cyl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)
cc            case (2)
cc              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
cc        case ('hel')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm
cc          phi = (curv(2)-aa*curv(3))
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ cos(phi),sin(phi),0d0 /)
cc            case (2)
cc              vec = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
cc        case ('tor')
cc          major_r = gparams(1)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ sin(curv(2))*cos(curv(3))
cc     .                ,sin(curv(2))*sin(curv(3))
cc     .                ,cos(curv(2)) /)
cc            case (2)
cc              vec = (/ cos(curv(2))*cos(curv(3))
cc     .                ,cos(curv(2))*sin(curv(3))
cc     .                ,-sin(curv(2))/)/curv(1)
cc            case (3)
cc              vec = (/ -sin(curv(3)),cos(curv(3)),0d0 /)
cc     .              /(major_r + curv(1)*sin(curv(2)))
cc          end select
cc        case ('sin')
cc
cc          pi = acos(-1d0)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          eps = gparams(1)
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ (xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)*
cc     -                                        Sin((2*Pi*x1)/xmax)))/
cc     -                 (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                            + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              , (-2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
cc     .                               *Sin((2*Pi*x1)/xmax))/
cc     -                (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))  
cc     -                           + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              , 0d0 /)
cc            case (2)
cc              vec = (/ (-2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
cc     .                                *Sin((2*Pi*x2)/ymax))/
cc     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                            + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              , (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)
cc     .                                      *Sin((2*Pi*x2)/ymax)))/
cc     -                (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                           + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              ,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
cc        case default
cc          write (*,*) 'Grid not implemented in covariantVector'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc      end function covariantVector
cc
ccc     contravariantVector2
ccc     #################################################################
cc      function contravariantVector(comp,i,j,k,igx,igy,igz) result (vec)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates contravariant vectors of curvilinear coordinate system
ccc     in Cartesian coordinates
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: comp,i,j,k,igx,igy,igz
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: car(3),curv(3),jac
cc        real(8)    :: x1,x2,x3,vec(3)
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          select case (comp)
cc            case (1)
cc              vec = (/ 1d0,0d0,0d0 /)
cc            case (2)
cc              vec = (/ 0d0,1d0,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
cc       case ('scl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gparams(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          jac = 1./jac
cc          select case (comp)
cc            case (1)
cc              vec = (/ jac,0d0,0d0 /)
cc            case (2)
cc              vec = (/ 0d0,1d0,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,jac /)
cc          end select
cc        case ('cyl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
cc            case (2)
cc              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)/curv(1)
cc          end select
cc        case ('hel')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm
cc          phi = (curv(2)-aa*curv(3))
cc          select case (comp)
cc            case (1)
cc              vec = (/ cos(phi),sin(phi),0d0 /)/curv(1)
cc            case (2)
cc              vec = (/-sin(phi),cos(phi),0d0 /)
cc            case (3)
cc              vec = (/ aa*sin(phi),-aa*cos(phi),1d0/curv(1) /)
cc          end select
cc        case ('tor')
cc          major_r = gparams(1)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ sin(curv(2))*cos(curv(3))
cc     .                ,sin(curv(2))*sin(curv(3))
cc     .                ,cos(curv(2)) /)
cc     .              /curv(1)/(major_r + curv(1)*sin(curv(2)))
cc            case (2)
cc              vec = (/ cos(curv(2))*cos(curv(3))
cc     .                ,cos(curv(2))*sin(curv(3))
cc     .                ,-sin(curv(2))/)
cc     .              /(major_r + curv(1)*sin(curv(2)))
cc            case (3)
cc              vec = (/ -sin(curv(3)),cos(curv(3)),0d0 /)/curv(1)
cc          end select
cc        case ('sin')
cc
cc          pi = acos(-1d0)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          eps = gparams(1)
cc
cc          select case (comp)
cc            case (1)
cc              vec = (/ (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)*
cc     -                                        Sin((2*Pi*x2)/ymax)))/
cc     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                            + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              , (2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
cc     .                              *Sin((2*Pi*x2)/ymax))/
cc     -                (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                           + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              ,0d0 /)
cc            case (2)
cc              vec = (/ (2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
cc     .                               *Sin((2*Pi*x1)/xmax))/
cc     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                            + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              ,(xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)
cc     .                                     *Sin((2*Pi*x1)/xmax)))/
cc     -               (xmax*ymax + eps*Pi*(xmax - ymax)
cc     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                          + eps*Pi*(xmax + ymax)
cc     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .              ,0d0 /)
cc            case (3)
cc              vec = (/ 0d0
cc     .                ,0d0
cc     .                ,(xmax*ymax)/(xmax*ymax
cc     .                  + eps*Pi*(xmax - ymax)
cc     -                       *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                  + eps*Pi*(xmax + ymax)
cc     .                       *Sin(2*Pi*(x1/xmax + x2/ymax))) /)
cc          end select
cc        case default
cc          write (*,*) 'Grid not implemented in contravariantVector'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc      end function contravariantVector
cc
ccc     hessian22
ccc     #################################################################
cc      function hessian22(l,i,j,k,igx,igy,igz) result (tensor)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates hessian elements of curvilinear coordinate system in
ccc     the covariant basis, i.e., 
ccc            hessian[l](i,j) = J^2 <cnv(i)|grad(cov[l])|cnv(j)>
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: l,i,j,k,igx,igy,igz
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: car(3),curv(3),vec(9)
cc        real(8)    :: x1,x2,x3,tensor(3,3)
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
cc        case ('scl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gparams(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc
cc          select case (l)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 2.*ypp/(ypp**2-cc**2)
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc          end select
cc
cc          vec = -vec
cc
cc        case ('cyl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          select case (l)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = curv(1)
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = -1./curv(1)
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc          end select
cc
cc          vec = -vec
cc
cc        case ('hel')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm
cc          phi = (curv(2)-aa*curv(3))
cc          select case (l)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = curv(1)
cc              vec(6) = -aa*curv(1)
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = aa**2*curv(1)
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = -1./curv(1)
cc              vec(3) = aa/curv(1)
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc          end select
cc
cc          vec = -vec
cc
cc        case ('tor')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          major_r = gparams(1)
cc          select case (l)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = curv(1)
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = sin(curv(2))*(major_r + curv(1)*sin(curv(2)))
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = -1./curv(1)
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = (major_r + curv(1)*sin(curv(2)))*cos(curv(2))
cc     .                 /curv(1)
cc            case (3)
cc              vec(1) =  0d0
cc              vec(2) =  0d0
cc              vec(3) =  -sin(curv(2))/(major_r + curv(1)*Sin(curv(2)))
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = -curv(1)*cos(curv(2))
cc     .                 /(major_r + curv(1)*Sin(curv(2)))
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc          end select
cc
cc          vec = -vec
cc
cc        case ('sin')
cc
cc          pi = acos(-1d0)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          eps = gparams(1)
cc
cc          select case (l)
cc            case (1)
cc              vec(1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)*
cc     -                                    Sin((2*Pi*x2)/ymax))/
cc     -          (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -           Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -           eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
cc              vec(2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
cc     .                             *Cos((2*Pi*x2)/ymax))/
cc     -            (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -              Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
cc     .                                   *Sin((2*Pi*x2)/ymax))/
cc     -            (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -            Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)
cc     -                                   *Sin((2*Pi*x2)/ymax))/
cc     -                 (xmax*(xmax*ymax
cc     .                 + eps*Pi*(xmax - ymax)
cc     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                 + eps*Pi*(xmax + ymax)
cc     .                     *Sin(2*Pi*(x1/xmax + x2/ymax))))
cc              vec(2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
cc     .                             *Cos((2*Pi*x2)/ymax))/
cc     -                 (xmax*ymax
cc     .                 + eps*Pi*(xmax - ymax)
cc     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                 + eps*Pi*(xmax + ymax)
cc     .                     *Sin(2*Pi*(x1/xmax + x2/ymax)))
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
cc     .                                   *Sin((2*Pi*x2)/ymax))/
cc     -                  (ymax*(xmax*ymax
cc     .                  + eps*Pi*(xmax - ymax)
cc     -                      *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
cc     -                  + eps*Pi*(xmax + ymax)
cc     .                      *Sin(2*Pi*(x1/xmax + x2/ymax))))
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc          end select
cc
cc        case default
cc          write (*,*) 'Grid not implemented in hessian'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc        tensor = reshape(vec,(/3,3/))
cc
cc      end function hessian22
cc
ccc     christ_2knd
ccc     #################################################################
cc      function christ_2knd(i,j,k,igx,igy,igz) result (tensor)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates elements of Christoffel symbol of the second kind,
ccc     Gamma[l](i,j) defined as: 
ccc            gamma[l](i,j) = -J^2 <cnv(i)|grad(cov[l])|cnv(j)>
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: l,i,j,k,igx,igy,igz
cc        real(8)    :: tensor(3,3,3)
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        tensor(1,:,:) = hessian22(1,i,j,k,igx,igy,igz)
cc        tensor(2,:,:) = hessian22(2,i,j,k,igx,igy,igz)
cc        tensor(3,:,:) = hessian22(3,i,j,k,igx,igy,igz)
cc
ccc     End program
cc
cc      end function christ_2knd
cc
ccc     hessian_cnv
ccc     #################################################################
cc      function hessian_cnv(k,x1,x2,x3) result (tensor)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates elements of tensor grad(cnv) in a mixed coordinate 
ccc     system:
ccc              hessian_cnv[k](i,j) = J^2 <cnv(i)|grad(cnv[k])|cov(j)>
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: k
cc        real(8)    :: x1,x2,x3,tensor(3,3)
cc
ccc     Local variables
cc
cc        real(8)    :: vec(9),car(3),eps
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
cc        case ('scl')
cc          lambda = gparams(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*x2/ymax-1.)
cc          select case (k)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 2.*ypp/(ypp**2-cc**2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 0d0
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 0d0
cc              vec(5) = 0d0
cc              vec(6) = 2.*ypp/(ypp**2-cc**2)
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc          end select
cc        case ('cyl')
cc          select case (k)
cc            case (1)
cc              vec(1) = -1./x1
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 0d0
cc              vec(5) = 1./x1
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = -x1
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = -1./x1
cc              vec(4) = 0d0
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc          end select
cc        case ('hel')
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm 
cc          select case (k)
cc          case (1)
cc            vec(1) = -1./x1
cc            vec(2) = 0d0
cc            vec(3) = 0d0
cc            vec(4) = 0d0
cc            vec(5) = 1./x1
cc            vec(6) = 0d0
cc            vec(7) = 0d0
cc            vec(8) = -aa/x1
cc            vec(9) = 0d0
cc          case (2)
cc            vec(1) = 0d0
cc            vec(2) = 0d0
cc            vec(3) = 0d0
cc            vec(4) = -x1
cc            vec(5) = 0d0
cc            vec(6) = 0d0
cc            vec(7) = aa*x1
cc            vec(8) = 0d0
cc            vec(9) = 0d0
cc          case (3)
cc            vec(1) = 0d0
cc            vec(2) = -aa/x1
cc            vec(3) = -1./x1
cc            vec(4) = aa*x1
cc            vec(5) = 0d0
cc            vec(6) = 0d0
cc            vec(7) = -aa**2*x1
cc            vec(8) = 0d0
cc            vec(9) = 0d0
cc          end select
cc        case ('tor')
cc          major_r = gparams(1)
cc          select case (k)
cc            case (1)
cc              vec(1) = -(major_r + 2*x1*sin(x2))
cc     .                  /x1/(major_r + x1*sin(x2))
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = -x1*cos(x2)/(major_r + x1*sin(x2))
cc              vec(5) = 1./x1
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = sin(x2)/(major_r + x1*sin(x2))
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = -sin(x2)/(major_r + x1*Sin(x2))
cc              vec(3) = 0d0
cc              vec(4) = -x1
cc              vec(5) = -x1*cos(x2)/(major_r + x1*Sin(x2))
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) =  x1*cos(x2)/(major_r + x1*Sin(x2))
cc            case (3)
cc              vec(1) =  0d0
cc              vec(2) =  0d0
cc              vec(3) =  -1./x1
cc              vec(4) = 0d0
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = -sin(x2)*(major_r + x1*Sin(x2))
cc              vec(8) = -cos(x2)*(major_r + x1*Sin(x2))/x1
cc              vec(9) = 0d0
cc          end select
cc
cc        case default
cc          write (*,*) 'Grid not implemented in hessian_cnv'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc        tensor = transpose(reshape(vec, (/3,3/)))
cc
cc      end function hessian_cnv
cc
ccc     g_sub
ccc     #################################################################
cc      function g_sub(i,j,k,igx,igy,igz) result (tensor)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates contravariant metric tensor of curvilinear coordinate 
ccc     system
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        real(8)    :: tensor(3,3)
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: vec(9),jac,car(3),curv(3)
cc        real(8)    :: x1,x2,x3
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          vec = (/ 1d0, 0d0, 0d0
cc     .            ,0d0, 1d0, 0d0
cc     .            ,0d0, 0d0, 1d0 /)
cc        case ('scl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gparams(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          vec = (/ 1d0/jac, 0d0, 0d0
cc     .            ,0d0    , jac, 0d0
cc     .            ,0d0    , 0d0, 1d0/jac /)
cc        case ('cyl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          vec = (/ 1d0/curv(1), 0d0     , 0d0
cc     .            ,0d0        , curv(1) , 0d0
cc     .            ,0d0        , 0d0     , 1d0/curv(1) /)
cc        case ('hel')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm
cc          vec = (/ 1./curv(1),0d0, 0d0
cc     .            ,0d0, curv(1) ,-aa*curv(1)
cc     .            ,0d0,-aa*curv(1), 1./curv(1) + aa**2*curv(1) /)
cc        case ('tor')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          major_r = gparams(1)
cc          vec = (/ 1d0/curv(1)/(major_r + curv(1)*sin(curv(2))),0d0,0d0
cc     .            ,0d0, curv(1)/(major_r + curv(1)*sin(curv(2))), 0d0
cc     .            ,0d0, 0d0, (major_r/curv(1) + sin(curv(2))) /)
cc        case ('sin')
cc
cc          pi = acos(-1d0)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          eps = gparams(1)
cc
cc          vec(1) =
cc     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
cc     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
cc     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
cc     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
cc     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
cc
cc          vec(2) =
cc     .         (eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
cc     -        eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
cc     -        xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
cc     -        ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
cc     -        ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
cc     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
cc
cc          vec(3) = 0d0
cc          vec(4) = vec(2)
cc
cc          vec(5) =
cc     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
cc     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
cc     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
cc     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
cc     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))!
cc
cc          vec(6) = 0d0
cc          vec(7) = vec(3)
cc          vec(8) = vec(6)
cc
cc          vec(9) = (xmax*ymax)/
cc     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
cc
cc        case default
cc          write (*,*) 'Grid not implemented in g_sub'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc        tensor = reshape(vec, (/3,3/))
cc
cc      end function g_sub
cc
ccc     g_sup
ccc     #################################################################
cc      function g_sup(i,j,k,igx,igy,igz) result (tensor)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates covariant metric tensor of curvilinear coordinate 
ccc     system
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        real(8)    :: tensor(3,3)
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: vec(9),jac,car(3),curv(3)
cc        real(8)    :: x1,x2,x3
cc
ccc     Begin program
cc
cc        select case (coords)
cc        case ('car')
cc          vec = (/ 1d0, 0d0, 0d0
cc     .            ,0d0, 1d0, 0d0
cc     .            ,0d0, 0d0, 1d0 /)
cc       case ('scl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gparams(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          vec = (/ jac, 0d0    , 0d0
cc     .            ,0d0, 1d0/jac, 0d0
cc     .            ,0d0, 0d0    , jac /)
cc        case ('cyl')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          vec = (/ curv(1) , 0d0   , 0d0
cc     .            ,0d0, 1d0/curv(1), 0d0
cc     .            ,0d0, 0d0   , curv(1) /)
cc        case ('hel')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          mm = gparams(1)
cc          kk = gparams(2)
cc          aa = kk/mm
cc          vec = (/ curv(1),0d0,0d0
cc     .            ,0d0,1./curv(1) + aa**2*curv(1), aa*curv(1)
cc     .            ,0d0, aa*curv(1), curv(1) /)
cc        case ('tor')
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          major_r = gparams(1)
cc          vec = (/ curv(1)*(major_r + curv(1)*sin(curv(2))), 0d0, 0d0
cc     .            ,0d0, (major_r/curv(1) + sin(curv(2))), 0d0
cc     .            ,0d0, 0d0, curv(1)/(major_r + curv(1)*sin(curv(2))) /)
cc        case ('sin')
cc
cc          pi = acos(-1d0)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          eps = gparams(1)
cc
cc          vec(1) =
cc     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
cc     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
cc     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
cc     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
cc     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
cc
cc          vec(2) =
cc     .         -((eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
cc     -          eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
cc     -          xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
cc     -          ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -          xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
cc     -          ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
cc     -      (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
cc
cc          vec(3) = 0d0
cc          vec(4) = vec(2)
cc
cc          vec(5) =
cc     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
cc     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
cc     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
cc     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
cc     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
cc     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
cc
cc          vec(6) = 0d0
cc          vec(7) = vec(3)
cc          vec(8) = vec(6)
cc
cc          vec(9) = (xmax*ymax + eps*Pi*(xmax - ymax)*
cc     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
cc     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
cc     .         /(xmax*ymax)
cc
cc        case default
cc          write (*,*) 'Grid not implemented in g_sup'
cc          write (*,*) 'Aborting...'
cc          stop
cc        end select
cc
cc        tensor = reshape(vec, (/3,3/))
cc
cc      end function g_sup

      end module grid_definition

c module bc_def
c #####################################################################
      module bc_def

c     -----------------------------------------------------------------
c     This module defines the types of BCs possible:
c     Topological BCs (set topology of computational domain):
c       PER: periodic
c       SP : singular point (as in cylindrical, toroidal)
c     Standard BCs:
c       EQU: Inhomogeneous dirichlet (generally imposed by equilibrium)
c       DIR: Homogeneous dirichlet
c       NEU: Neumann (homogeneous/inhomogeneous)
c       DEF: Default (set in application code)
c       EXT: By extrapolation
c     The actual integers chosen indicate a particular order that the BC
c     scheduler follows to impose the boundary conditions.
c     -----------------------------------------------------------------

        use grid_definition

        implicit none

        integer(4) :: PER,DIR,NEU,SP,EQU,DEF,EXT
        parameter (EQU=1,EXT=3,PER=0,NEU=2,DIR=1,SP=6,DEF=7)

        integer(4),parameter :: BM=1,TP=2,LT=3,RT=4,DN=5,UP=6

        integer(4),parameter :: BCLIM=5

cc        integer(4) :: bcond(6)
        
cc      contains
cc
ccc     isSP
ccc     #################################################################
cc      function isSP(i,j,k,igx,igy,igz) result(sing_point)
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        logical    :: sing_point
cc
cc
cc        sing_point=.false.
cc        if (i+grid_params%ilo(igx)-1 == 1 .and. bcond(1) == SP)
cc     .       sing_point=.true.
cc
cc      end function isSP

      end module bc_def

ccc module grid_metric
ccc #####################################################################
cc      module grid_metric
cc
ccc ---------------------------------------------------------------------
ccc     This module packs routines that perform operations on grid
ccc     quantities, such as coordinate transformation of vector components,
ccc     vector norms and scalar products. It contains the following
ccc     routines:
ccc          * transformVectorToCartesian 
ccc          * transformVectorToCurvilinear 
ccc          * transformFromCurvToCurv 
ccc          * volume
ccc          * vectorNorm
ccc          * scalarProduct
ccc     It is assumed that the grid metric structure gmetric is allocated 
ccc     and filled.
ccc ---------------------------------------------------------------------
cc
cc        use grid_definition
cc
cc        use bc_def
cc
cc        implicit none
cc
cc        type :: grid_metrics
cccc          real(8),pointer,dimension(:,:,:,:)     :: car   !Node positions in cartesian space
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
cc        type(MG_grid_metrics) :: gmetric
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
cc        type(MG_grid_metrics) :: gmetric
cc
ccc     Local variables
cc
cc        integer(4)      :: igrid,ilom,ihip,jlom,jhip,klom,khip
cc
ccc     Begin program
cc
cc        if (.not.associated(gmetric%grid)) then
cc          allocate(gmetric%grid(grid_params%ngrid))
cc          do igrid=1,grid_params%ngrid
cc            ilom = 0
cc            jlom = 0
cc            klom = 0
cc            ihip = grid_params%nxv(igrid)+1
cc            jhip = grid_params%nyv(igrid)+1
cc            khip = grid_params%nzv(igrid)+1
cccc            allocate(gmetric%grid(igrid)
cccc     $           %car  (ilom:ihip,jlom:jhip,klom:khip,3))
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
cc        type(MG_grid_metrics) :: gmetric
cc
ccc     Local variables
cc
cc        integer(4)      :: igrid,nxp,nyp,nzp
cc
ccc     Begin program
cc
cc        if (associated(gmetric%grid)) then
cc          do igrid=1,grid_params%ngrid
cccc            deallocate(gmetric%grid(igrid)%car  )
cc            deallocate(gmetric%grid(igrid)%jac  )
cc            deallocate(gmetric%grid(igrid)%gsub )
cc            deallocate(gmetric%grid(igrid)%gsup )
cc            deallocate(gmetric%grid(igrid)%cov  )
cc            deallocate(gmetric%grid(igrid)%cnv  )
cc            deallocate(gmetric%grid(igrid)%Gamma)
cc          enddo
cc          deallocate(gmetric%grid)
cc        endif
cc
ccc     End program
cc
cc      end subroutine deallocateGridMetric
cc
ccc     defineGridMetric
ccc     #################################################################
cc      subroutine defineGridMetric(gmetric,igr)
cc
ccc     -----------------------------------------------------------------
ccc     This routine calculates all grid metric quantities required for
ccc     the curvilinear representation of a set of PDE's: jacobian,
ccc     metric tensors, covariant and contravariant vectors, Christoffel
ccc     symbols of the second kind. All quantities are stored in
ccc     structure gmetric. There are two modes of computation:
ccc        * Analytical (numerical_grid=.false.)
ccc        * Numerical  (numerical_grid=.true.)
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4),optional,intent(IN) :: igr
cc
cc        type(MG_grid_metrics) :: gmetric
cc
ccc     Local variables
cc
cc        integer(4) :: igrmin,igrmax,igrid,ilo,ihi,jlo,jhi,klo,khi
cc     $               ,ilom,ihip,jlom,jhip,klom,khip
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc     .               ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
cc     .               ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp
cc        real(8)    :: r(3,3),car0(3),carp(3),carm(3),dh(3),jac,ijac
cc     .               ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
cc     .               ,gamma(3,3,3),gamm1(3,3,3),mag,dhp,dhm,dhh
cc        real(8),allocatable,dimension(:,:,:,:,:) :: dr
cc
ccc     Interpolation
cc
cc        integer(4) :: kx,ky,kz,nnx,nny,nnz,dim,flg,order
cc
cc        real(8)    :: xp,yp,zp
cc        real(8),allocatable,dimension(:) :: sx,sy,sz
cc        real(8), dimension(:),allocatable:: tx,ty,tz,work
cc        real(8), dimension(:,:,:,:,:),allocatable:: drbcoef
cc        real(8), dimension(:,:,:,:,:,:),allocatable:: gambcoef
cc
cc        real(8)    :: db3val
cc        external      db3val
cc
ccc     Begin program
cc
cc        anal_map = checkGridDatabase()
cc
cc        call allocateGridMetric(gmetric)
cc
cc        if (PRESENT(igr)) then
cc          igrmin = igr
cc          igrmax = igr
cc        else
cc          igrmin = 1
cc          igrmax = grid_params%ngrid
cc        endif
cc
ccc       ANALYTICAL QUANTITIES
cc        if ((.not.numerical_grid).and.anal_map) then
cc
cc          do igrid=igrmin,igrmax
cc
cc            igx = igrid
cc            igy = igrid
cc            igz = igrid
cc
cc            ilo = 1
cc            jlo = 1
cc            klo = 1
cc            ihi = grid_params%nxv(igrid)
cc            jhi = grid_params%nyv(igrid)
cc            khi = grid_params%nzv(igrid)
cc
cc            ilom = ilo-1
cc            jlom = jlo-1
cc            klom = klo-1
cc            ihip = ihi+1
cc            jhip = jhi+1
cc            khip = khi+1
cc
cc            do k = klom,khip
cc              do j = jlom,jhip
cc                do i = ilom,ihip
cc                  gmetric%grid(igrid)%jac  (i,j,k)
cc     .                      = jacobian(i,j,k,igx,igy,igz)
cc                  gmetric%grid(igrid)%gsub (i,j,k,:,:)
cc     .                      = g_sub   (i,j,k,igx,igy,igz)
cc                  gmetric%grid(igrid)%gsup (i,j,k,:,:)
cc     .                      = g_sup   (i,j,k,igx,igy,igz)
cc                  gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc     .                      = christ_2knd(i,j,k,igx,igy,igz)
cc                  do l=1,3
cc                    gmetric%grid(igrid)%cov(i,j,k,l,:)
cc     .                      = covariantVector    (l,i,j,k,igx,igy,igz)
cc                    gmetric%grid(igrid)%cnv(i,j,k,l,:)
cc     .                      = contravariantVector(l,i,j,k,igx,igy,igz)
cc                  enddo
cc                enddo
cc              enddo
cc            enddo
cc
cc            !Zero force condition on Christoffel symbols (only on finest grid)
cc            if (igrid == 1) then
cc              do k = klo,khi
cc                do j = jlo,jhi
cc                  do i = ilo,ihi
cc                    call gammaZeroForce(i,j,k,igx,igy,igz)
cc                  enddo
cc                enddo
cc              enddo
cc            endif
cc
cc          enddo
cc
ccc       NUMERICAL QUANTITIES
cc        else
cc
cc          do igrid=igrmin,igrmax
cc
cc            igx = igrid
cc            igy = igrid
cc            igz = igrid
cc
cc            ilo = 1
cc            jlo = 1
cc            klo = 1
cc            ihi = grid_params%nxv(igrid)
cc            jhi = grid_params%nyv(igrid)
cc            khi = grid_params%nzv(igrid)
cc
cc            ilom = ilo-1
cc            jlom = jlo-1
cc            klom = klo-1
cc            ihip = ihi+1
cc            jhip = jhi+1
cc            khip = khi+1
cc
cc            allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))
cc
cc            !Evaluate dx/dxi vectors
cc            do k = klom,khip
cc              do j = jlom,jhip
cc                do i = ilom,ihip
cc                  ip=min(i+1,ihip)
cc                  im=max(i-1,ilom)
cc                  jp=min(j+1,jhip)
cc                  jm=max(j-1,jlom)
cc                  kp=min(k+1,khip)
cc                  km=max(k-1,klom)
cc
cc                  carp = map(ip,j,k,igx,igy,igz,igp,jg,kg)
cc                  carm = map(im,j,k,igx,igy,igz,igm,jg,kg)
cc                  dh(1)= (grid_params%xx(igp)-grid_params%xx(igm))
cc                  dr(i,j,k,1,:) = (carp-carm)/dh(1)
cc
cc                  carp = map(i,jp,k,igx,igy,igz,ig,jgp,kg)
cc                  carm = map(i,jm,k,igx,igy,igz,ig,jgm,kg)
cc                  dh(2)= (grid_params%yy(jgp)-grid_params%yy(jgm))
cc                  dr(i,j,k,2,:) = (carp-carm)/dh(2)
cc
cc                  carp = map(i,j,kp,igx,igy,igz,ig,jg,kgp)
cc                  carm = map(i,j,km,igx,igy,igz,ig,jg,kgm)
cc                  dh(3)= (grid_params%zz(kgp)-grid_params%zz(kgm))
cc                  dr(i,j,k,3,:) = (carp-carm)/dh(3)
cc
cc                enddo
cc              enddo
cc            enddo
cc
cc            !Enforce topological constraints on dr
cc            do i=1,3
cc              do j=1,3
cc                call topol_bc(dr(:,:,:,i,j))
cc              enddo
cc            enddo
cc
cccc      if (my_rank == 0) then
cccc        open(unit=110,file='debug0.bin'
cccc     .           ,form='unformatted',status='replace')
cccc      else
cccc        open(unit=110,file='debug1.bin'
cccc     .           ,form='unformatted',status='replace')
cccc      endif
cccc      igx = 1
cccc      do i=1,3
cccc        do j=1,3
cccc          nnx = grid_params%nxv(igx)
cccc          nny = grid_params%nyv(igx)
cccc          nnz = grid_params%nzv(igx)
cccccc          debug = gmetric%grid(igx)%cnv(0:nx+1,0:ny+1,0:nz+1,i,j)
cccccc          debug = gmetric%grid(igx)%Gamma(0:nx+1,0:ny+1,0:nz+1,1,i,j)
cccccc     .           *gmetric%grid(igx)%jac
cccc
cccc          call contour(dr(0:nnx+1,0:nny+1,2,i,j),nnx+2,nny+2,0d0
cccc     .                ,xmax,0d0,ymax,i+j-2,110)
cccc        enddo
cccc      enddo
cccc      close(110)
cccc      call MPI_Finalize(mpierr)
cccc      stop
cc
cc          !Spline dr for interpolation on coarser grids
cccc          order = 2
cccc
cccc          nnx = nxp+1
cccc          nny = nyp+1
cccc          nnz = nzp+1
cccc
cccc          allocate(sx(nnx),sy(nny),sz(nnz))
cccc
cccc          call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)
cccc
cccc          sx(1:nnx) = grid_params%xx(ig-1:ig+nxp-1)
cccc          sy(1:nny) = grid_params%yy(jg-1:jg+nyp-1)
cccc          sz(1:nnz) = grid_params%zz(kg-1:kg+nzp-1)
cccc
cccc          flg = 0
cccc
cccc          kx = min(order+1,nnx-1)
cccc          ky = min(order+1,nny-1)
cccc          kz = min(order+1,nnz-1)
cccc
cccc          dim = nnx*nny*nnz+max(2*kx*(nnx+1),2*ky*(nny+1),2*kz*(nnz+1))
cccc
cccc          allocate(tx(nnx+kx))
cccc          allocate(ty(nny+ky))
cccc          allocate(tz(nnz+kz))
cccc          allocate(work(dim))
cccc          allocate(drbcoef(nnx,nny,nnz,3,3))
cccc
cccc          do i=1,3
cccc            do j=1,3
cccc              call db3ink(sx,nnx,sy,nny,sz,nnz,dr(:,:,:,i,j),nnx,nny
cccc     .                   ,kx,ky,kz,tx,ty,tz,drbcoef(:,:,:,i,j)
cccc     .                   ,work,flg)
cccc            enddo
cccc          enddo
cc
cc          !Evaluate grid quantities
cc            do k = klom,khip
cc              do j = jlom,jhip
cc                do i = ilom,ihip
cc                  r = dr(i,j,k,:,:)
cc
cc                  !Evaluate Jacobian
cc                  jac = triple_product(r(1,:),r(2,:),r(3,:))
cc                  ijac = 1d0/jac
cc
cc                  !Contravariant vectors
cc                  cnv(:,:) = r(:,:)*ijac
cc
cc                  !Covariant vectors
cc                  cov(1,:) = cross_product(r(2,:),r(3,:))*ijac
cc                  cov(2,:) = cross_product(r(3,:),r(1,:))*ijac
cc                  cov(3,:) = cross_product(r(1,:),r(2,:))*ijac
cc
cc                  !Metric tensors
cc                  do m=1,3
cc                    do l=m,3
cc                      gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
cc                      gsub(m,l) = gsub(l,m) !Symmetry
cc                      gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
cc                      gsup(m,l) = gsup(l,m) !Symmetry
cc                    enddo
cc                  enddo
cc
cc                  !Grid spacings
cc                  call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc                  dh(1) = grid_params%dxh(ig)
cc                  dh(2) = grid_params%dyh(jg)
cc                  dh(3) = grid_params%dzh(kg)
cc
cc                  !Christoffel symbols
cc
cc                  if (i==ilom) then
cc                    i0=ilo
cc                  elseif (i==ihip) then
cc                    i0=ihi
cc                  else
cc                    i0=i
cc                  endif
cc
cc                  if (j==jlom) then
cc                    j0=jlo
cc                  elseif (j==jhip) then
cc                    j0=jhi
cc                  else
cc                    j0=j
cc                  endif
cc
cc                  if (k==klom) then
cc                    k0=klo
cc                  elseif (k==khip) then
cc                    k0=khi
cc                  else
cc                    k0=k
cc                  endif
cc
cc                  ip=i0+1
cc                  im=i0-1
cc                  jp=j0+1
cc                  jm=j0-1
cc                  kp=k0+1
cc                  km=k0-1
cc
cc                  xp = grid_params%xx(ig)
cc                  yp = grid_params%yy(jg)
cc                  zp = grid_params%zz(kg)
cc
cc                  do l=1,3
cc                    do m=1,3
cc                      do n=m,3
cc                        if (m == n) then
cc                          select case(m)
cc                          case(1)
cc                            car0 = map(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
cc                            carp = map(ip,j0,k0,igx,igy,igz,igp,jg,kg)
cc                            carm = map(im,j0,k0,igx,igy,igz,igm,jg,kg)
cc                            dhp = (grid_params%xx(igp)
cc     .                            -grid_params%xx(ig0))
cc                            dhm = (grid_params%xx(ig0)
cc     .                            -grid_params%xx(igm))
cc                          case(2)
cc                            car0 = map(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
cc                            carp = map(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
cc                            carm = map(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
cc                            dhp = (grid_params%yy(jgp)
cc     .                            -grid_params%yy(jg0))
cc                            dhm = (grid_params%yy(jg0)
cc     .                            -grid_params%yy(jgm))
cc                          case(3)
cc                            car0 = map(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
cc                            carp = map(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
cc                            carm = map(i0,j0,km,igx,igy,igz,ig,jg,kgm)
cc                            dhp = (grid_params%zz(kgp)
cc     .                            -grid_params%zz(kg0))
cc                            dhm = (grid_params%zz(kg0)
cc     .                            -grid_params%zz(kgm))
cc                          end select
cc                          carp = (carp-car0)/dhp
cc                          carm = (car0-carm)/dhm
cc                          dhh  = 0.5*(dhp+dhm)
cc                        else
cc                          select case(m)
cc                          case(1)
cc                           call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
cc                           call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
cc                           call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
cc                           carp= 0.5*(dr(ip,j0,k0,n,:)+dr(i0,j0,k0,n,:))
cc                           carm= 0.5*(dr(im,j0,k0,n,:)+dr(i0,j0,k0,n,:))
cc                           dhh = 0.5*(grid_params%xx(igp)
cc     .                               -grid_params%xx(igm))
cc                          case(2)
cc                           call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
cc                           call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
cc                           call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
cc                           carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
cc                           carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
cc                           dhh = 0.5*(grid_params%yy(jgp)
cc     .                               -grid_params%yy(jgm))
cc                          case(3)
cc                           call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
cc                           call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
cc                           call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
cc                           carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
cc                           carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
cc                           dhh = 0.5*(grid_params%zz(kgp)
cc     .                               -grid_params%zz(kgm))
cc                          end select
cc                        endif
cc                        vec = (carp-carm)/dhh
cc
cccc                      select case(m)
cccc                      case(1)
cccc                        do p=1,3
cccc                          vec(p)=db3val(xp,yp,zp,1,0,0,tx,ty,tz
cccc     .                          ,nnx,nny,nnz,kx,ky,kz,drbcoef(:,:,:,n,p)
cccc     .                          ,work)
cccc                        enddo
cccc                      case(2)
cccc                        do p=1,3
cccc                          vec(p)=db3val(xp,yp,zp,0,1,0,tx,ty,tz
cccc     .                          ,nnx,nny,nnz,kx,ky,kz,drbcoef(:,:,:,n,p)
cccc     .                          ,work)
cccc                        enddo
cccc                      case(3)
cccc                        do p=1,3
cccc                          vec(p)=db3val(xp,yp,zp,0,0,1,tx,ty,tz
cccc     .                          ,nnx,nny,nnz,kx,ky,kz,drbcoef(:,:,:,n,p)
cccc     .                          ,work)
cccc                        enddo
cccc                      end select
cc
cc                        gamma(l,m,n) = dot_product(vec,cov(l,:))
cc                        gamma(l,n,m) = gamma(l,m,n) !Symmetry
cc                      enddo
cc                    enddo
cc                  enddo
cc
cccc                  gamma = christ_2knd(i,j,k,igx,igy,igz)
cc
cc                  !Store grid quantities
cc                  gmetric%grid(igrid)%jac  (i,j,k)       = jac
cc                  gmetric%grid(igrid)%gsub (i,j,k,:,:)   = gsub
cc                  gmetric%grid(igrid)%gsup (i,j,k,:,:)   = gsup
cc                  gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
cc                  gmetric%grid(igrid)%cov  (i,j,k,:,:)   = cov
cc                  gmetric%grid(igrid)%cnv  (i,j,k,:,:)   = cnv
cc
cc                enddo
cc              enddo
cc            enddo
cc
cc            deallocate(dr)
cc
cc            !Zero-force condition on Christoffle symbols 
cc            if (igrid == 1) then
cc              do k = klo,khi
cc                do j = jlo,jhi
cc                  do i = ilo,ihi
cc                    call gammaZeroForce(i,j,k,igx,igy,igz)
cc                  enddo
cc                enddo
cc              enddo
cc            endif
cc
cc            !Enforce topological constraints on Christoffel symbols
cc            do i=1,3
cc              do j=1,3
cc                do k=1,3
cc                  call topol_bc(gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
cc                enddo
cc              enddo
cc            enddo
cc
cccc      if (my_rank == 0) then
cccc        open(unit=110,file='debug0.bin'
cccc     .           ,form='unformatted',status='replace')
cccc      else
cccc        open(unit=110,file='debug1.bin'
cccc     .           ,form='unformatted',status='replace')
cccc      endif
cccc      igx = 1
cccc      do i=1,3
cccc        do j=1,3
cccc          nnx = grid_params%nxv(igx)
cccc          nny = grid_params%nyv(igx)
cccc          nnz = grid_params%nzv(igx)
cccccc          debug = gmetric%grid(igx)%cnv(0:nx+1,0:ny+1,0:nz+1,i,j)
cccccc          debug = gmetric%grid(igx)%Gamma(0:nx+1,0:ny+1,0:nz+1,1,i,j)
cccccc     .           *gmetric%grid(igx)%jac
cccc
cccc          call contour(gmetric%grid(igrid)%Gamma(:,:,1,i,j,1)
cccc     .                ,nnx+2,nny+2
cccc     .                ,0d0,xmax,0d0,ymax,i+j-2,110)
cccc        enddo
cccc      enddo
cccc      close(110)
cccc      call MPI_Finalize(mpierr)
cccc      stop
cc
cccc        if (my_rank == 0) write (*,*) 'grid_level',igrid
cc
cc          !Spline Christoffel symbols for interpolation on coarser grids
cccc          allocate(gambcoef(nnx,nny,nnz,3,3,3))
cccc          flg = 0
cccc
cccc          do k=1,3
cccc            do j=1,3
cccc              do i=1,3
cccc                call db3ink(sx,nnx,sy,nny,sz,nnz
cccc     .                     ,gmetric%grid(igrid)%Gamma(:,:,:,i,j,k)
cccc     .                     ,nnx,nny,kx,ky,kz,tx,ty,tz
cccc     .                     ,gambcoef(:,:,:,i,j,k),work,flg)
cccc              enddo
cccc            enddo
cccc          enddo
cc
ccc       Restrict grid metric info to coarser grids
cc
cccc          call restrictGridMetrics
cc
ccc       Deallocate spline auxiliary variables
cc
cccc          deallocate(sx,sy,sz,tx,ty,tz,work,drbcoef,gambcoef)
cc
cc          enddo
cc
ccc       Test of numerical calculation vs. analytical calculation
cc
cccc          do igrid=1,grid_params%ngrid
cccc
cccc            igx = igrid
cccc            igy = igrid
cccc            igz = igrid
cccc
cccc            nxp = grid_params%nxv(igrid)+1
cccc            nyp = grid_params%nyv(igrid)+1
cccc            nzp = grid_params%nzv(igrid)+1
cccc
cccc            do k = 0,nzp
cccc              do j = 0,nyp
cccc                do i = 0,nxp
cccc                  gmetric%grid(igrid)%jac  (i,j,k)
cccc     .                      = gmetric%grid(igrid)%jac  (i,j,k)
cccc     .                       -jacobian(i,j,k,igx,igy,igz)
cccc                  gmetric%grid(igrid)%gsub (i,j,k,:,:)
cccc     .                      = gmetric%grid(igrid)%gsub (i,j,k,:,:)
cccc     .                       -g_sub   (i,j,k,igx,igy,igz)
cccc                  gmetric%grid(igrid)%gsup (i,j,k,:,:)
cccc     .                      = gmetric%grid(igrid)%gsup (i,j,k,:,:)
cccc     .                       -g_sup   (i,j,k,igx,igy,igz)
cccc                  gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cccc     .                      = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cccc     .                       -christ_2knd(i,j,k,igx,igy,igz)
cccc                enddo
cccc              enddo
cccc            enddo
cccc
cccc            write (*,*) 'Grid level:',igrid
cccc          
cccc            mag =
cccc     .     sum(gmetric%grid(igrid)%jac(1:nxp-1,1:nyp-1,1:nzp-1)**2)
cccc            write (*,*)'Jacobian tst=',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cccc            mag =
cccc     .     sum(gmetric%grid(igrid)%gsub(2:nxp-1,1:nyp-1,1:nzp-1,:,:)**2)
cccc            write (*,*)'Gsub tst    =',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cccc            mag =
cccc     .     sum(gmetric%grid(igrid)%gsup(2:nxp-1,1:nyp-1,1:nzp-1,:,:)**2)
cccc            write (*,*)'Gsup tst    =',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cccc            mag =
cccc     .  sum(gmetric%grid(igrid)%Gamma(2:nxp-1,1:nyp-1,1:nzp-1,:,:,:)**2)
cccc            write (*,*)'Gamma tst   =',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cccc
cccc          enddo
cccc          stop
cc
cc        endif
cc
ccc     End program
cc
cc      contains
cc
ccc     contour
ccc     #####################################################################
cc      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
cc      implicit none               !For safe fortran
ccc     ---------------------------------------------------------------------
ccc     Contours 2D array in xdraw format. In call:
ccc       * arr: 2D array to be plotted
ccc       * nx,ny: dimensions of array
ccc       * xmin,xmax,ymin,ymax: 2D domain limits
ccc       * iopt: whether to initialize xdraw plot (iopt=0) or not.
ccc       * nunit: integer file identifier.
ccc     ---------------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer(4) :: nx,ny,iopt,nunit
cc      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax
cc
ccc     Local variables
cc
cc      integer(4) :: i,j
cc
ccc     Begin program
cc
cc      if(iopt == 0) then
cc        write(nunit) nx-1,ny-1,0
cc        write(nunit) real(xmin,4),real(xmax,4)
cc     .              ,real(ymin,4),real(ymax,4) 
cc      endif
cc      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cccc      write (*,*) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cccc      write (*,*)
cc
ccc     End program
cc
cc      end subroutine contour
cc
ccc     map
ccc     #################################################################
cc      function map(i,j,k,igx,igy,igz,ig,jg,kg) result(car)
cc
ccc     -----------------------------------------------------------------
ccc     Give Cartesian coordinates corresponding to node (i,j,k) at grid
ccc     level (igx,igy,igz).
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
cc        real(8)    :: car(3)
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                              ,car(1),car(2),car(3))
cc
cc      end function map
cc
ccc$$$c     restrictGridMetrics
ccc$$$c     #################################################################
ccc$$$      subroutine restrictGridMetrics
ccc$$$c     -----------------------------------------------------------------
ccc$$$c     Restricts arrays with geometric grid info.
ccc$$$c     -----------------------------------------------------------------
ccc$$$
ccc$$$      implicit none    !For safe fortran
ccc$$$
ccc$$$c     Call variables
ccc$$$
ccc$$$c     Local variables
ccc$$$
ccc$$$      real (8) :: xp,yp,zp
ccc$$$
ccc$$$c     Begin program
ccc$$$
ccc$$$      do igrid=2,grid_params%ngrid
ccc$$$
ccc$$$        igx = igrid
ccc$$$        igy = igrid
ccc$$$        igz = igrid
ccc$$$
ccc$$$        nxp = grid_params%nxv(igrid)+1
ccc$$$        nyp = grid_params%nyv(igrid)+1
ccc$$$        nzp = grid_params%nzv(igrid)+1
ccc$$$
ccc$$$        allocate(dr(0:nxp,0:nyp,0:nzp,3,3))
ccc$$$
ccc$$$        !Evaluate dx/dxi vectors (splined)
ccc$$$        do k = 0,nzp
ccc$$$          do j = 0,nyp
ccc$$$            do i = 0,nxp
ccc$$$              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
ccc$$$
ccc$$$              xp = grid_params%xx(ig)
ccc$$$              yp = grid_params%yy(jg)
ccc$$$              zp = grid_params%zz(kg)
ccc$$$
ccc$$$              do l=1,3
ccc$$$                do m=1,3
ccc$$$                  dr(i,j,k,l,m) =
ccc$$$     .                db3val(xp,yp,zp,0,0,0,tx,ty,tz,nnx,nny,nnz
ccc$$$     .                      ,kx,ky,kz,drbcoef(:,:,:,l,m),work)
ccc$$$                enddo
ccc$$$              enddo
ccc$$$
ccc$$$            enddo
ccc$$$          enddo
ccc$$$        enddo
ccc$$$
ccc$$$        !Enforce topological constraints on dr
ccc$$$        do i=1,3
ccc$$$          do j=1,3
ccc$$$            call topol_bc(dr(:,:,:,i,j))
ccc$$$          enddo
ccc$$$        enddo
ccc$$$
ccc$$$        !Evaluate grid quantities
ccc$$$        do k = 0,nzp
ccc$$$          do j = 0,nyp
ccc$$$            do i = 0,nxp
ccc$$$
ccc$$$              r = dr(i,j,k,:,:)
ccc$$$
ccc$$$              !Evaluate Jacobian
ccc$$$              jac = triple_product(r(1,:),r(2,:),r(3,:))
ccc$$$              ijac = 1d0/jac
ccc$$$
ccc$$$              !Contravariant vectors
ccc$$$              cnv(:,:) = r(:,:)*ijac
ccc$$$
ccc$$$              !Covariant vectors
ccc$$$              cov(1,:) = cross_product(r(2,:),r(3,:))*ijac
ccc$$$              cov(2,:) = cross_product(r(3,:),r(1,:))*ijac
ccc$$$              cov(3,:) = cross_product(r(1,:),r(2,:))*ijac
ccc$$$
ccc$$$              !Metric tensors
ccc$$$              do m=1,3
ccc$$$                do l=m,3
ccc$$$                  gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
ccc$$$                  gsub(m,l) = gsub(l,m) !Symmetry
ccc$$$                  gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
ccc$$$                  gsup(m,l) = gsup(l,m) !Symmetry
ccc$$$                enddo
ccc$$$              enddo
ccc$$$
ccc$$$              !Grid spacings
ccc$$$              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
ccc$$$              dh(1) = grid_params%dxh(ig)
ccc$$$              dh(2) = grid_params%dyh(jg)
ccc$$$              dh(3) = grid_params%dzh(kg)
ccc$$$
ccc$$$              !Christoffel symbols (splined)
ccc$$$
ccc$$$                if (i==0) then
ccc$$$                  i0=1
ccc$$$                elseif (i==nxp) then
ccc$$$                  i0=nxp-1
ccc$$$                else
ccc$$$                  i0=i
ccc$$$                endif
ccc$$$
ccc$$$                if (j==0) then
ccc$$$                  j0=0
ccc$$$                elseif (j==nyp) then
ccc$$$                  j0=nyp
ccc$$$                else
ccc$$$                  j0=j
ccc$$$                endif
ccc$$$
ccc$$$                if (k==0) then
ccc$$$                  k0=1
ccc$$$                elseif (k==nzp) then
ccc$$$                  k0=nzp-1
ccc$$$                else
ccc$$$                  k0=k
ccc$$$                endif
ccc$$$
ccc$$$                ip=i0+1
ccc$$$                im=i0-1
ccc$$$                jp=j0+1
ccc$$$                jm=j0-1
ccc$$$                kp=k0+1
ccc$$$                km=k0-1
ccc$$$
ccc$$$              xp = grid_params%xx(ig)
ccc$$$              yp = grid_params%yy(jg)
ccc$$$              zp = grid_params%zz(kg)
ccc$$$
ccc$$$              do l=1,3
ccc$$$                do m=1,3
ccc$$$                  do n=m,3
ccc$$$cc                    gamma(l,m,n) =
ccc$$$cc     .                db3val(xp,yp,zp,0,0,0,tx,ty,tz,nnx,nny,nnz
ccc$$$cc     .                      ,kx,ky,kz,gambcoef(:,:,:,l,m,n)
ccc$$$cc     .                      ,work)
ccc$$$
ccc$$$                    select case(m)
ccc$$$                    case(1)
ccc$$$                      call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
ccc$$$                      call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
ccc$$$                      call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
ccc$$$                      carp= 0.5*(dr(ip,j0,k0,n,:)+dr(i0,j0,k0,n,:))
ccc$$$                      carm= 0.5*(dr(im,j0,k0,n,:)+dr(i0,j0,k0,n,:))
ccc$$$                      dhh = 0.5*(grid_params%xx(igp)
ccc$$$     .                              -grid_params%xx(igm))
ccc$$$                    case(2)
ccc$$$                      call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
ccc$$$                      call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
ccc$$$                      call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
ccc$$$                      carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
ccc$$$                      carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
ccc$$$                      dhh = 0.5*(grid_params%yy(jgp)
ccc$$$     .                              -grid_params%yy(jgm))
ccc$$$                    case(3)
ccc$$$                      call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
ccc$$$                      call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
ccc$$$                      call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
ccc$$$                      carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
ccc$$$                      carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
ccc$$$                      dhh = 0.5*(grid_params%zz(kgp)
ccc$$$     .                              -grid_params%zz(kgm))
ccc$$$                    end select
ccc$$$                    vec = (carp-carm)/dhh
ccc$$$
ccc$$$                    gamma(l,m,n) = dot_product(vec,cov(l,:))
ccc$$$                    gamma(l,n,m) = gamma(l,m,n) !Symmetry
ccc$$$                  enddo
ccc$$$                enddo
ccc$$$              enddo
ccc$$$
ccc$$$              !Store grid quantities
ccc$$$              gmetric%grid(igrid)%jac(i,j,k) = jac
ccc$$$              gmetric%grid(igrid)%gsub (i,j,k,:,:)   = gsub
ccc$$$              gmetric%grid(igrid)%gsup (i,j,k,:,:)   = gsup
ccc$$$              gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
ccc$$$              gmetric%grid(igrid)%cov  (i,j,k,:,:)   = cov
ccc$$$              gmetric%grid(igrid)%cnv  (i,j,k,:,:)   = cnv
ccc$$$
ccc$$$            enddo
ccc$$$          enddo
ccc$$$        enddo
ccc$$$
ccc$$$        !Enforce topological constraints on Christoffel symbols
ccc$$$        do i=1,3
ccc$$$          do j=1,3
ccc$$$            do k=1,3
ccc$$$              call topol_bc(gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
ccc$$$            enddo
ccc$$$          enddo
ccc$$$        enddo
ccc$$$
ccc$$$        !Deallocate auxiliary variables
ccc$$$        deallocate(dr)
ccc$$$
ccc$$$      enddo
ccc$$$
ccc$$$c     End program
ccc$$$
ccc$$$      end subroutine restrictGridMetrics
cc
ccc     topol_bc
ccc     #################################################################
cc      subroutine topol_bc(array)
cc
ccc     -----------------------------------------------------------------
ccc     !Enforce topological constraints on array
ccc     -----------------------------------------------------------------
cc
cc        implicit  none
cc
ccc     Call variables
cc
cc        real(8)    :: array(ilom:ihip,jlom:jhip,klom:khip)
cc
ccc     Local variables
cc
cc        integer(4) :: nxl,nyl,nzl,nxg,nyg,nzg,dim
cc        integer(4) :: my_per_rank,dest_rank,count,request1,request2
cc
cc        real(8),allocatable,dimension(:,:) :: put,get
cc
ccc     Begin program
cc
cc        nxl = grid_params%nxv(igrid)
cc        nyl = grid_params%nyv(igrid)
cc        nzl = grid_params%nzv(igrid)
cc
cc#if !defined(petsc)
cc
cc        if (bcond(1) == PER .or. bcond(2) == PER) then
cc          array(0    ,:,:) = array(nxl,:,:)
cc          array(nxl+1,:,:) = array(1  ,:,:)
cc        endif
cc
cc        if (bcond(3) == PER .or. bcond(4) == PER) then
cc          array(:,0    ,:) = array(:,nyl,:)
cc          array(:,nyl+1,:) = array(:,1  ,:)
cc        endif
cc        
cc        if (bcond(5) == PER .or. bcond(6) == PER) then
cc          array(:,:,0    ) = array(:,:,nzl)
cc          array(:,:,nzl+1) = array(:,:,1  )
cc        endif
cc
cc#else
cc
cc        nxg = grid_params%nxgl(1)
cc        nyg = grid_params%nygl(1)
cc        nzg = grid_params%nzgl(1)
cc
ccc       X-direction
cc
cc        dim = 1
cc
cc        call findDestProc(dim,my_per_rank,dest_rank)
cc
cc        if (dest_rank >= 0) then
cc
cc          !X-direction: SEND
cc          if (bcond(1) == PER .and. grid_params%ilo(1)==1) then
cc            tag_send=2
cc            count = size(array(1,:,:))
cc            call MPI_ISend
cc     .         (array(1,:,:),count,MPI_REAL8,dest_rank,tag_send
cc     .         ,MPI_COMM_PER(dim),request1,mpierr)
cc          endif
cc
cc          if (bcond(1) == PER .and. grid_params%ihi(1)==nxg) then
cc            tag_send=1
cc            count = size(array(1,:,:))
cc            call MPI_ISend
cc     .         (array(nxl,:,:),count,MPI_REAL8,dest_rank,tag_send
cc     .         ,MPI_COMM_PER(dim),request2,mpierr)
cc          endif
cc
cc          !X-direction: RECV
cc          if (bcond(1) == PER .and. grid_params%ilo(1)==1) then
cc            tag_recv=1
cc            count = size(array(1,:,:))
cc            call MPI_Recv
cc     .         (array(0,:,:),count,MPI_REAL8,dest_rank,tag_recv
cc     .         ,MPI_COMM_PER(dim),status,mpierr)
cccc            array(0,:,:) = array(nxl,:,:)
cc          endif
cc
cc          if (bcond(1) == PER .and. grid_params%ihi(1)==nxg) then
cc            tag_recv=2
cc            count = size(array(1,:,:))
cc            call MPI_Recv
cc     .         (array(nxl+1,:,:),count,MPI_REAL8,dest_rank,tag_recv
cc     .         ,MPI_COMM_PER(dim),status,mpierr)
cccc            array(nxl+1,:,:) = array(1,:,:)
cc          endif
cc
cc        endif
cc
ccc       Y-direction
cc
cc        dim = 2
cc
cc        call findDestProc(dim,my_per_rank,dest_rank)
cc
cc        if (dest_rank >= 0) then
cc
cc          !Y-direction: SEND
cc          if  (bcond(3) == PER .and. grid_params%jlo(1)==1) then
cc            tag_send=4
cc            count = size(array(:,1,:))
cc            call MPI_Isend
cc     .        (array(:,1,:),count,MPI_REAL8,dest_rank,tag_send
cc     .        ,MPI_COMM_PER(dim),request1,mpierr)
cccc            write (*,*) dest_rank,tag_send,MPI_COMM_PER(dim),count
cccc            write (*,*) array(0:nxl+1,1  ,0:nzl+1)
cc          endif
cc
cc          if (bcond(3) == PER .and. grid_params%jhi(1)==nyg) then
cc            tag_send=3
cc            count = size(array(:,1,:))
cc            call MPI_Isend
cc     .        (array(:,nyl  ,:),count,MPI_REAL8,dest_rank,tag_send
cc     .        ,MPI_COMM_PER(dim),request2,mpierr)
cccc            write (*,*) 'send',dest_rank,tag_send,MPI_COMM_PER(dim),count
cccc            write (*,*) array(0:nxl+1,nyl,0:nzl+1)
cc          endif
cc
cc          !Y-direction: RECV
cc          if  (bcond(3) == PER .and. grid_params%jlo(1)==1) then
cc            tag_recv=3
cc            count = size(array(:,1,:))
cc            call MPI_Recv
cc     .        (array(:,0,:),count,MPI_REAL8,dest_rank,tag_recv
cc     .        ,MPI_COMM_PER(dim),status,mpierr)
cccc            write (*,*)'recv',dest_rank,tag_recv,MPI_COMM_PER(dim),count
cccc            write (*,*) array(0:nxl+1,0,0:nzl+1)
cccc            array(:,0,:) = array(:,nyl,:)
cc          endif
cc
cc          if (bcond(3) == PER .and. grid_params%jhi(1)==nyg) then
cc            tag_recv=4
cc            count = size(array(:,1,:))
cc            call MPI_Recv
cc     .        (array(:,nyl+1,:),count,MPI_REAL8,dest_rank,tag_recv
cc     .        ,MPI_COMM_PER(dim),status,mpierr)
cccc          write (*,*) dest_rank,tag_recv,MPI_COMM_PER(dim),count
cccc          write (*,*) array(0:nxl+1,nyl+1,0:nzl+1)
cccc          array(:,nyl+1,:) = array(:,1,:)
cc          endif
cc
cc        endif
cc
ccc diag ****
cccc        if (igrid > 1) then
cccc          call MPI_Finalize(mpierr)
cccc          stop
cccc        endif
ccc diag ****
cc
ccc       Z-direction
cc
cc        dim = 3
cc
cc        call findDestProc(dim,my_per_rank,dest_rank)
cc
cc        if (dest_rank >= 0) then
cc
cc          !Z-direction: SEND
cc          if (bcond(5) == PER .and. grid_params%klo(1)==1) then
cc            tag_send=6
cc            count=size(array(:,:,1))
cc            call MPI_ISend
cc     .         (array(:,:,1),count,MPI_REAL8,dest_rank,tag_send
cc     .         ,MPI_COMM_PER(dim),request1,mpierr)
cc          endif
cc
cc          if (bcond(5) == PER .and. grid_params%khi(1)==nzg) then
cc            tag_send=5
cc            count=size(array(:,:,nzl))
cc            call MPI_ISend
cc     .         (array(:,:,nzl),count,MPI_REAL8,dest_rank,tag_send
cc     .         ,MPI_COMM_PER(dim),request2,mpierr)
cc          endif
cc
cc          !Z-direction: RECV
cc          if (bcond(5) == PER .and. grid_params%klo(1)==1) then
cc            tag_recv=5
cc            count=size(array(:,:,1))
cc            call MPI_Recv
cc     .         (array(:,:,0),count,MPI_REAL8,dest_rank,tag_recv
cc     .         ,MPI_COMM_PER(dim),status,mpierr)
cccc          array(:,:,0  ) = array(:,:,nzl)
cc          endif
cc
cc          if (bcond(5) == PER .and. grid_params%khi(1)==nzg) then
cc            tag_recv=6
cc            count=size(array(:,:,nzl))
cc            call MPI_Recv
cc     .         (array(:,:,nzl+1),count,MPI_REAL8,dest_rank,tag_recv
cc     .         ,MPI_COMM_PER(dim),status,mpierr)
cccc          array(:,:,nzl+1) = array(:,:,1    )
cc          endif
cc
cc        endif
cc
cc#endif
cc
ccc     End program
cc
cc      end subroutine topol_bc
cc
cc#if defined(petsc)
ccc     findDestProc
ccc     ##################################################################
cc      subroutine findDestProc(dim,my_per_rank,dest_rank)
cc
cc        implicit none
cc
cc        integer(4) :: dest_rank,dim,my_per_size,my_per_rank
cc
cc        if (MPI_COMM_PER(dim) == MPI_COMM_NULL) then
cc          dest_rank = -1
cc          return
cc        endif
cc
cc        call MPI_Comm_rank(MPI_COMM_PER(dim),my_per_rank,mpierr)
cc        call MPI_Comm_size(MPI_COMM_PER(dim),my_per_size,mpierr)
cc
cc        if (my_per_size > 2) then
cc          write (*,*) 'There are more than 2 processors in this'
cc          write (*,*) 'periodic communicator'
cc          write (*,*) 'Aborting topol_bc...'
cc          call MPI_Finalize(mpierr)
cc          stop
cc        else
cc          dest_rank = my_per_size -1 -my_per_rank
cc        endif
cc
cc      end subroutine findDestProc
cc#endif
cc
ccc     gammaZeroForce
ccc     #################################################################
cc      subroutine gammaZeroForce(i,j,k,igx,igy,igz)
cc
ccc     -----------------------------------------------------------------
ccc     Postprocess Christoffel symbol to satisfy zero-force condition
ccc     to machine accuracy
ccc     -----------------------------------------------------------------
cc
cc        implicit  none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc
ccc     Local variables
cc
cc        integer(4) :: ii,ll,kk,mm,jj,ig,jg,kg
cc        real(8)    :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const
cc        real(8)    :: hess(3,3,3),table(3,3,3),gsub(3,3)
cc     .               ,gsup(3,3),gsupp(3,3),gsupm(3,3),gamma(3,3,3)
cc
ccc     Begin program
cc
cc        if  (    (i == ilom .or. i == ihip)       !Avoid x-boundaries
cc     .       .or.(j == jlom .or. j == jhip)       !Avoid y-boundaries
cc     .       .or.(k == klom .or. k == khip)       !Avoid z-boundaries
cc     .       .or. isSP(i,j,k,igx,igy,igz)         !Avoid singular point
cccc     .       .or.(grid_params%ilo(igx) == 1
cccc     .            .and. bcond(1) == SP     )      !Avoid singular point
cc     .       ) return
cc
cc        gamma = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc        gsub  = gmetric%grid(igrid)%gsub(i,j,k,:,:)
cc        gsup  = gmetric%grid(igrid)%gsup(i,j,k,:,:)
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        do ii=1,3 !cycle through Chistoffel symbols
cc
cc          if (ii == 2 .and. bcond(1) == SP) then
cc            jac  = gmetric%grid(igrid)%jac(i,j,k)
cc            const= 2d0
cc          else
cc            jac=1d0
cc            const=1d0
cc          endif
cc
cc          do ll=1,3
cc            do kk=1,3
cc              select case (kk)
cc              case(1)
cc                dh = 2*grid_params%dxh(ig)
cc                if (ii == 2 .and. bcond(1) == SP) then
cc                  jacp = gmetric%grid(igrid)%jac(i+1,j,k)
cc                  jacm = gmetric%grid(igrid)%jac(i-1,j,k)
cc                else
cc                  jacp=1d0
cc                  jacm=1d0
cc                endif
cc                gsupp = jacp*gmetric%grid(igrid)%gsup (i+1,j,k,:,:)
cc                gsupm = jacm*gmetric%grid(igrid)%gsup (i-1,j,k,:,:)
cc              case(2)
cc                dh = 2*grid_params%dyh(jg)
cc                if (ii == 2 .and. bcond(1) == SP) then
cc                  jacp = gmetric%grid(igrid)%jac(i,j+1,k)
cc                  jacm = gmetric%grid(igrid)%jac(i,j-1,k)
cc                else
cc                  jacp=1d0
cc                  jacm=1d0
cc                endif
cc                gsupp = jacp*gmetric%grid(igrid)%gsup (i,j+1,k,:,:)
cc                gsupm = jacm*gmetric%grid(igrid)%gsup (i,j-1,k,:,:)
cc              case(3)
cc                dh = 2*grid_params%dzh(kg)
cc                if (ii == 2 .and. bcond(1) == SP) then
cc                  jacp = gmetric%grid(igrid)%jac(i,j,k+1)
cc                  jacm = gmetric%grid(igrid)%jac(i,j,k-1)
cc                else
cc                  jacp=1d0
cc                  jacm=1d0
cc                endif
cc                gsupp = jacp*gmetric%grid(igrid)%gsup (i,j,k+1,:,:)
cc                gsupm = jacm*gmetric%grid(igrid)%gsup (i,j,k-1,:,:)
cc              end select
cc                
cc              table(ii,ll,kk) =
cc     .              -gsub(ll,1)*(gsupp(1,ii)-gsupm(1,ii))/dh/jac
cc     .              -gsub(ll,2)*(gsupp(2,ii)-gsupm(2,ii))/dh/jac
cc     .              -gsub(ll,3)*(gsupp(3,ii)-gsupm(3,ii))/dh/jac
cc     .              +const*delta(ii,ll)*(gamma(1,kk,1)
cc     .                                  +gamma(2,kk,2)
cc     .                                  +gamma(3,kk,3))
cc
cc              summ=0d0
cc              do mm=1,3
cc                do jj=1,3
cc                  summ = summ + gsub(ll,mm)*gsup(ii,jj)*gamma(mm,jj,kk)
cc                enddo
cc              enddo
cc
cc              table(ii,ll,kk) = table(ii,ll,kk) - summ
cc            enddo
cc          enddo
cc
cc        enddo
cc
cc        gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = table
cc
ccc      START CHECKS
cc
cccc        write (*,*)
cccc        write (*,*) 'Grid:',igrid,'  Grid node:',i,j,k
cc
ccc      Check difference between correction and prediction
cc
cccc        table = table - gamma
cccc
cccc        summ2= sqrt(sum(gamma*gamma))
cccc        summ = sqrt(sum(table*table))/summ2
cccc
cccc        write (*,*) summ
cc
ccc      Check zero-force property: g^lk G^i_lk = -nabla_m(g^mi)
cc
cccc        do ii=1,3 !cycle through Chistoffel symbols
cccc
cccc          summ = 0d0
cccc
cccc          dh = 2*grid_params%dxh(ig)
cccc          gsupp = gmetric%grid(igrid)%gsup (i+1,j,k,:,:)
cccc          gsupm = gmetric%grid(igrid)%gsup (i-1,j,k,:,:)
cccc
cccc          summ = summ -(gsupp(1,ii)-gsupm(1,ii))/dh
cccc
cccc          dh = 2*grid_params%dyh(jg)
cccc          gsupp = gmetric%grid(igrid)%gsup (i,j+1,k,:,:)
cccc          gsupm = gmetric%grid(igrid)%gsup (i,j-1,k,:,:)
cccc
cccc          summ = summ -(gsupp(2,ii)-gsupm(2,ii))/dh
cccc
cccc          dh = 2*grid_params%dzh(kg)
cccc          gsupp = gmetric%grid(igrid)%gsup (i,j,k+1,:,:)
cccc          gsupm = gmetric%grid(igrid)%gsup (i,j,k-1,:,:)
cccc
cccc          summ = summ -(gsupp(3,ii)-gsupm(3,ii))/dh
cccc
cccc          summ2 = 0d0
cccc          do ll=1,3
cccc            do kk=1,3
cccc              summ2 = summ2 + gsup(ll,kk)*gamma(ii,ll,kk)
cccc            enddo
cccc          enddo
cccc
cccc          write (*,*) summ-summ2
cccc        enddo
cc
ccc      Check cancellation property g^lk(d_il G^j_kj-g_lm g^ij G^m_jk) = 0
cc
cccc        do ii=1,3 !cycle through Chistoffel symbols
cccc          do ll=1,3
cccc            do kk=1,3
cccc              table(ii,ll,kk) =
cccc     .               delta(ii,ll)*(gamma(1,kk,1)
cccc     .                            +gamma(2,kk,2)
cccc     .                            +gamma(3,kk,3))
cccc
cccc              summ=0d0
cccc              do mm=1,3
cccc                do jj=1,3
cccc                  summ = summ + gsub(ll,mm)*gsup(ii,jj)*gamma(mm,jj,kk)
cccc                enddo
cccc              enddo
cccc
cccc              table(ii,ll,kk) = table(ii,ll,kk) - summ
cccc            enddo
cccc          enddo
cccc
cccc          summ = 0d0
cccc          do ll=1,3
cccc            do kk=1,3
cccc              summ = summ + gsup(ll,kk)*table(ii,ll,kk)
cccc            enddo
cccc          enddo
cccc
cccc          write (*,*) summ
cccc        enddo
cc
ccc     End program
cc
cc      end subroutine gammaZeroForce
cc
cc      end subroutine defineGridMetric
cc
ccc     cross_product
ccc     #################################################################
cc      function cross_product(vec1,vec2) result(vec3)
cc
ccc     -----------------------------------------------------------------
ccc     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
ccc     where the vectors are in Cartesian coordinates.
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        real(8)    :: vec1(3),vec2(3),vec3(3)
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        vec3(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
cc        vec3(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
cc        vec3(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)
cc
cc      end function cross_product
cc
ccc     triple_product
ccc     #################################################################
cc      function triple_product(vec1,vec2,vec3) result(scalar)
cc
ccc     -----------------------------------------------------------------
ccc     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
ccc     where the vectors are in Cartesian coordinates.
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        real(8)    :: vec1(3),vec2(3),vec3(3),scalar
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        scalar = dot_product(vec1,cross_product(vec2,vec3))
cc
cc      end function triple_product
cc
ccc     delta
ccc     #################################################################
cc      real(8) function delta(i,j)
cc
cc        integer(4) :: i,j
cc
cc        delta = 0d0
cc        if (i == j) delta = 1d0
cc
cc      end function delta
cc
cc      end module grid_metric
cc
ccc module grid_operations
ccc #####################################################################
cc      module grid_operations
cc
ccc ---------------------------------------------------------------------
ccc     This module packs routines that perform operations on grid
ccc     quantities, such as coordinate transformation of vector components,
ccc     vector norms and scalar products. It contains the following
ccc     routines:
ccc          * transformVectorToCartesian 
ccc          * transformVectorToCurvilinear 
ccc          * transformFromCurvToCurv 
ccc          * volume
ccc          * vectorNorm
ccc          * scalarProduct
ccc     It is assumed that the grid metric structure gmetric is allocated 
ccc     and filled.
ccc ---------------------------------------------------------------------
cc
cc        use grid_metric
cc
cc        implicit none
cc
cc      contains
cc
ccc     transformVectorToCartesian
ccc     #################################################################
cc      subroutine transformVectorToCartesian(i,j,k,igx,igy,igz
cc     .                                     ,c1,c2,c3,covariant
cc     .                                     ,cx,cy,cz)
cc
ccc     -----------------------------------------------------------------
ccc     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
ccc     at grid coordinates (i,j,k). Curvilinear vector is covariant if
ccc     covariant=.true., and contravariant otherwise.
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        real(8)    :: c1,c2,c3,cx,cy,cz
cc        logical    :: covariant
cc        
cc
ccc     Local variables
cc
cc        integer(4) :: ic
cc        real(8)    :: T_to_car(3,3),vec(3)
cc
ccc     Begin program
cc
cc        if (covariant) then
cc
cc          do ic =1,3
cc            T_to_car(:,ic) = gmetric%grid(igx)%cov(i,j,k,ic,:)
cc          enddo
cc
cc          vec = (/ c1,c2,c3 /)
cc
cc          vec = matmul(T_to_car,vec)
cc
cc          cx = vec(1)
cc          cy = vec(2)
cc          cz = vec(3)
cc            
cc        else
cc
cc          do ic =1,3
cc            T_to_car(:,ic) = gmetric%grid(igx)%cnv(i,j,k,ic,:)
cc          enddo
cc
cc          vec = (/ c1,c2,c3 /)
cc
cc          vec = matmul(T_to_car,vec)
cc
cc          cx = vec(1)
cc          cy = vec(2)
cc          cz = vec(3)
cc
cc        endif
cc
cc      end subroutine transformVectorToCartesian
cc
ccc     transformVectorToCurvilinear
ccc     #################################################################
cc      subroutine transformVectorToCurvilinear(i,j,k,igx,igy,igz
cc     .                                       ,cx,cy,cz,covariant
cc     .                                       ,c1,c2,c3)
cc
ccc     -----------------------------------------------------------------
ccc     Transforms a Cartesian vector (cx,cy,cz) to curvilinear (c1,c2,c3)
ccc     (covariant if covariant=.true., contravariant otherwise)
ccc     at grid coordinates (i,j,k). 
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        real(8)    :: c1,c2,c3,cx,cy,cz
cc        logical    :: covariant
cc
ccc     Local variables
cc
cc        integer(4) :: ic
cc        real(8)    :: T_to_curv(3,3),vec(3),jac
cc
ccc     Begin program
cc
cc        jac = gmetric%grid(igx)%jac(i,j,k)
cc
cc        if (covariant) then
cc
cc          do ic =1,3
cc            T_to_curv(:,ic) = gmetric%grid(igx)%cnv(i,j,k,ic,:)
cc          enddo
cc
cc          T_to_curv = jac*transpose(T_to_curv)
cc
cc          vec = (/ cx,cy,cz /)
cc
cc          vec = matmul( T_to_curv,vec)
cc
cc          c1 = vec(1)
cc          c2 = vec(2)
cc          c3 = vec(3)
cc            
cc        else
cc
cc          do ic =1,3
cc            T_to_curv(:,ic) = gmetric%grid(igx)%cov(i,j,k,ic,:)
cc          enddo
cc
cc          T_to_curv = jac*transpose(T_to_curv)
cc
cc          vec = (/ cx,cy,cz /)
cc
cc          vec = matmul(T_to_curv,vec)
cc
cc          c1 = vec(1)
cc          c2 = vec(2)
cc          c3 = vec(3)
cc
cc        endif
cc
cc      end subroutine transformVectorToCurvilinear
cc
ccc     transformFromCurvToCurv
ccc     #################################################################
cc      subroutine transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv,half_elem)
ccc     -----------------------------------------------------------------
ccc     Transforms a curvilinear vector from covariant to contravariant 
ccc     (tocnv=.true.) or viceversa at grid coordinates (i,j,k).
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        integer(4),optional,intent(IN) :: half_elem
cc
cc        real(8)    :: cov1,cov2,cov3,cnv1,cnv2,cnv3
cc        logical    :: tocnv
cc
ccc     Local variables
cc
cc        integer(4) :: he,ip
cc        real(8)    :: cov(3),cnv(3),tensor(3,3)
cc
ccc     Begin program
cc
cc        if (PRESENT(half_elem)) then
cc          he = half_elem
cc        else
cc          he = 0
cc        endif
cc
cc        select case(he)
cc        case(1)
cc          ip = i+1
cc          if (isSP(ip,j,k,igx,igy,igz)) ip = i  !SP is located at i=-1 ghost cell
cc          if (tocnv) then
cc            tensor = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
cc     .                   +gmetric%grid(igx)%gsup(i ,j,k,:,:))
cc          else
cc            tensor = 0.5*(gmetric%grid(igx)%gsub(ip,j,k,:,:)
cc     .                   +gmetric%grid(igx)%gsub(i ,j,k,:,:))
cc          endif
cc        case(2)
cc          if (tocnv) then
cc            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j+1,k,:,:)
cc     .                   +gmetric%grid(igx)%gsup(i,j  ,k,:,:))
cc          else                                        
cc            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j+1,k,:,:)
cc     .                   +gmetric%grid(igx)%gsub(i,j  ,k,:,:))
cc          endif
cc        case(3)
cc          if (tocnv) then
cc            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j,k+1,:,:)
cc     .                   +gmetric%grid(igx)%gsup(i,j,k  ,:,:))
cc          else                                          
cc            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j,k+1,:,:)
cc     .                   +gmetric%grid(igx)%gsub(i,j,k  ,:,:))
cc          endif
cc        case(0)
cc          if (tocnv) then
cc            tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
cc          else
cc            tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
cc          endif
cc        end select
cc
cc        if (tocnv) then
cc          cnv  = matmul(tensor,(/ cov1,cov2,cov3 /))
cc          cnv1 = cnv(1)
cc          cnv2 = cnv(2)
cc          cnv3 = cnv(3)
cc        else
cc          cov  = matmul(tensor,(/ cnv1,cnv2,cnv3 /))
cc          cov1 = cov(1)
cc          cov2 = cov(2)
cc          cov3 = cov(3)
cc        endif
cc
cc      end subroutine transformFromCurvToCurv
cc
ccc     volume
ccc     #################################################################
cc      function volume(i,j,k,igx,igy,igz) result(vol)
cc
ccc     -----------------------------------------------------------------
ccc     Calculates Jacobian of curvilinear coordinate system
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Input variables
cc
cc        integer(4) :: i,j,k,igx,igy,igz
cc        real(8)    :: vol
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: x1,x2,x3,dx1,dx2,dx3,jac
cc
ccc     Begin program
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        dx1 = grid_params%dxh(ig)
cc        dx2 = grid_params%dyh(jg)
cc        dx3 = grid_params%dzh(kg)
cc
cc        jac = gmetric%grid(igx)%jac(i,j,k)
cc        
cc        vol = jac*dx1*dx2*dx3
cc
cc      end function volume
cc
ccc     vectorNorm
ccc     ################################################################
cc      real(8) function vectorNorm(i,j,k,igx,igy,igz,ax,ay,az,covar)
cc
ccc     ---------------------------------------------------------------
ccc     Finds norm of vector A given its curvilinear components.
ccc     ---------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k,igx,igy,igz
cc      real(8)    :: ax,ay,az
cc      logical    :: covar
cc
ccc     Local variables
cc
cc      real(8)    :: tensor(3,3),cnv(3),cov(3),jac
cc
ccc     Begin program
cc
cc      jac = gmetric%grid(igx)%jac(i,j,k)
cc
cc      if (covar) then
cc        tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
cc        cov = (/ ax,ay,az /)
cc        cnv = matmul(tensor,cov)
cc      else
cc        tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
cc        cnv = (/ ax,ay,az /)
cc        cov = matmul(tensor,cnv)
cc      endif
cc
cc      vectorNorm = dot_product(cov,cnv)/jac
cc
ccc     End 
cc
cc      end function vectorNorm
cc
ccc     scalarProduct
ccc     ################################################################
cc      function scalarProduct(i,j,k,igx,igy,igz,cov1,cov2,cov3
cc     .                       ,cnv1,cnv2,cnv3) result (dot)
cc
ccc     ---------------------------------------------------------------
ccc     Finds scalar product of two vectors, one covariant and the
ccc     other contravariant.
ccc     ---------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k,igx,igy,igz
cc      real(8)    :: dot,cov1,cov2,cov3,cnv1,cnv2,cnv3
cc
ccc     Local variables
cc
cc      real(8)    :: cnv(3),cov(3)
cc
ccc     Begin program
cc
cc      cnv = (/ cnv1,cnv2,cnv3 /)
cc      cov = (/ cov1,cov2,cov3 /)
cc
cc      dot = dot_product(cov,cnv)/gmetric%grid(igx)%jac(i,j,k)
cc
ccc     End 
cc
cc      end function scalarProduct
cc
cc      end module grid_operations

c module grid
c #####################################################################
      module grid

        use grid_definition

        use bc_def

        implicit none

        integer(4),private :: nxx,nyy,nzz,bcond(6)

        integer(4) :: mg_ratio

        real(8),private :: pi

      contains

c     createGrid
c     #################################################################
      subroutine createGrid(nx,ny,nz,bcs)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities from:
c       * nx,ny: global grid size
c       * bcond: specifies boundary conditions
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nx,ny,nz,bcs(6)

c     Local variables

        integer(4) :: i,n1,n2,n3,imin,imax,jmin,jmax,kmin,kmax
        integer(4) :: ngrdx,ngrdy,ngrdz

c     Begin program

        mg_ratio = 2

        nxx = nx
        nyy = ny
        nzz = nz

        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
        kmin = 1
        kmax = nz

        bcond = bcs

        pi = acos(-1d0)

c     Find adequate number of grid levels (for MG)

        n1 = int(dlog(1d0*nxx)/dlog(1d0*mg_ratio)+0.001)
        n2 = int(dlog(1d0*nyy)/dlog(1d0*mg_ratio)+0.001)
        n3 = int(dlog(1d0*nzz)/dlog(1d0*mg_ratio)+0.001)

        ngrdx = max(n1-1,1)
        do i = ngrdx,1,-1
          n1 = nxx/mg_ratio**(i-1)
          if (n1*mg_ratio**(i-1).eq.nxx) exit
        enddo
        ngrdx = i

        ngrdy = max(n2-1,1)
        do i = ngrdy,1,-1
          n2 = nyy/mg_ratio**(i-1)
          if (n2*mg_ratio**(i-1).eq.nyy) exit
        enddo
        ngrdy = i

        ngrdz = max(n3-1,1)
        do i = ngrdz,1,-1
          n3 = nzz/mg_ratio**(i-1)
          if (n3*mg_ratio**(i-1).eq.nzz) exit
        enddo
        ngrdz = i

c     Allocate grid storage structure

        call allocateGridStructure(nxx,nyy,nzz,ngrdx,ngrdy,ngrdz
     .                            ,grid_params)

c     Initialize MG grid-size arrays

        grid_params%mg_ratio_x = 1
        grid_params%mg_ratio_y = 1
        grid_params%mg_ratio_z = 1

        grid_params%nxgl(1) = nx
        grid_params%nxv (1) = nxx
        grid_params%ilo (1) = imin
        grid_params%ihi (1) = imax
        do i = 2,ngrdx
          grid_params%nxgl(i) = grid_params%nxgl(i-1) / mg_ratio
          grid_params%nxv (i) = grid_params%nxv (i-1) / mg_ratio
cc          grid_params%ilo(i) = max(grid_params%ilo(i-1) / mg_ratio,1)
          grid_params%ilo (i) = imin
          grid_params%ihi (i) = grid_params%ilo(i)+grid_params%nxv(i)-1
          grid_params%mg_ratio_x(i-1) = mg_ratio
        enddo
        do i = ngrdx+1,grid_params%ngrid
          grid_params%nxgl(i) = grid_params%nxgl(i-1)
          grid_params%nxv (i) = grid_params%nxv (i-1)
          grid_params%ilo (i) = grid_params%ilo (i-1)
          grid_params%ihi (i) = grid_params%ihi (i-1)
        enddo

        grid_params%nygl(1) = ny
        grid_params%nyv (1) = nyy
        grid_params%jlo (1) = jmin
        grid_params%jhi (1) = jmax
        do i = 2,ngrdy
          grid_params%nygl(i) = grid_params%nygl(i-1) / mg_ratio
          grid_params%nyv (i) = grid_params%nyv(i-1) / mg_ratio
cc          grid_params%jlo(i) = max(grid_params%jlo(i-1) / mg_ratio,1)
          grid_params%jlo (i) = jmin
          grid_params%jhi (i) = grid_params%jlo(i)+grid_params%nyv(i)-1
          grid_params%mg_ratio_y(i-1) = mg_ratio
        enddo
        do i = ngrdy+1,grid_params%ngrid
          grid_params%nygl(i) = grid_params%nygl(i-1)
          grid_params%nyv (i) = grid_params%nyv (i-1)
          grid_params%jlo (i) = grid_params%jlo (i-1)
          grid_params%jhi (i) = grid_params%jhi (i-1)
        enddo

        grid_params%nzgl(1) = nz
        grid_params%nzv (1) = nzz
        grid_params%klo (1) = kmin
        grid_params%khi (1) = kmax
        do i = 2,ngrdz
          grid_params%nzgl(i) = grid_params%nzgl(i-1) / mg_ratio
          grid_params%nzv (i) = grid_params%nzv (i-1) / mg_ratio
cc          grid_params%klo(i) = max(grid_params%klo(i-1) / mg_ratio,1)
          grid_params%klo (i) = kmin
          grid_params%khi (i) = grid_params%klo(i)+grid_params%nzv(i)-1
          grid_params%mg_ratio_z(i-1) = mg_ratio
        enddo
        do i = ngrdz+1,grid_params%ngrid
          grid_params%nzgl(i) = grid_params%nzgl(i-1)
          grid_params%nzv (i) = grid_params%nzv (i-1)
          grid_params%klo (i) = grid_params%klo (i-1)
          grid_params%khi (i) = grid_params%khi (i-1)
        enddo

c     Initialize MG pointer arrays

        grid_params%istartx(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istartx(i) = grid_params%istartx(i-1)
     .                            +(grid_params%nxv(i-1)+2)
        enddo

        grid_params%istarty(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istarty(i) = grid_params%istarty(i-1)
     .                            +(grid_params%nyv(i-1)+2)
        enddo

        grid_params%istartz(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istartz(i) = grid_params%istartz(i-1)
     .                            +(grid_params%nzv(i-1)+2)
        enddo

        grid_params%istartp(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istartp(i) = grid_params%istartp(i-1)
     .                            +grid_params%nxv(i-1)
     .                            *grid_params%nyv(i-1)
     .                            *grid_params%nzv(i-1)
        enddo

c     Set grid parameters

        grid_params%params = gparams

c     Consistency checks

        call consistencyCheck

c     Define uniform logical LOCAL grid on ALL grid levels

        call createLogicalGrid(nx,grid_params%xx,grid_params%dx
     .                        ,grid_params%dxh,grid_params%nxv
     .                        ,grid_params%ngrid,grid_params%istartx
     .                        ,grid_params%ilo,grid_params%ihi
     .                        ,xmin,xmax,bcond(LT),bcond(RT))

        call createLogicalGrid(ny,grid_params%yy,grid_params%dy
     .                        ,grid_params%dyh,grid_params%nyv
     .                        ,grid_params%ngrid,grid_params%istarty
     .                        ,grid_params%jlo,grid_params%jhi
     .                        ,ymin,ymax,bcond(BM),bcond(TP))

        call createLogicalGrid(nz,grid_params%zz,grid_params%dz
     .                        ,grid_params%dzh,grid_params%nzv
     .                        ,grid_params%ngrid,grid_params%istartz
     .                        ,grid_params%klo,grid_params%khi
     .                        ,zmin,zmax,bcond(DN),bcond(UP))

      end subroutine createGrid

c     createLogicalGrid
c     #################################################################
      subroutine createLogicalGrid (nglobal,xx,dx,dxh,nx,ngrid,istart
     .                             ,ilo,ihi,lmin,lmax,bcs1,bcs2)

        implicit none

c     Call variables

        integer(4) :: nglobal,ngrid,nx(ngrid),istart(ngrid)
     $               ,ilo(ngrid),ihi(ngrid),bcs1,bcs2
        real(8)    :: xx(*),dx(*),dxh(*),lmin,lmax

c     Local variables
        
        integer(4) :: ig,i,isig
        real(8)    :: dh,length,lmax_loc,lmin_loc,lstart

c     Begin program

        length = lmax-lmin

c     Find grid quantities

        do ig = 1,ngrid

          isig = istart(ig)

          !Find cell centers
          if (bcs1 == PER) then
            dh = length/dfloat(nx(ig)-1)
            lstart = 0d0
          elseif(bcs1 == DIR .and. bcs2 == DIR) then
            dh = length/dfloat(nx(ig)+1)
            lstart = dh
          elseif(bcs1 == NEU .and. bcs2 == DIR) then
            dh = length/dfloat(nx(ig))
            lstart = 0d0
          elseif(bcs1 == DIR .and. bcs2 == NEU) then
            dh = length/dfloat(nx(ig))
            lstart = dh
          elseif(bcs1 == NEU .and. bcs2 == NEU) then
            dh = length/dfloat(nx(ig)-1)
            lstart = 0d0
          endif

          xx(1 + isig) = lstart
          do i = 2,nx(ig)+1
            xx(i + isig) = xx(i-1 + isig) + dh
          enddo
          xx(0 + isig) = xx(1+isig) - dh

          !Find integer mesh spacings
          do i = 0,nx(ig)
            dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
          enddo

          !Find half mesh spacings
          dxh(1 + isig) = (dx(1 + isig) + dx(0   + isig))/2.
          do i = 2,nx(ig)
            dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
          enddo

          if (bcs1 == PER) then
            dxh(0        + isig) = dxh(nx(ig) + isig)
            dxh(nx(ig)+1 + isig) = dxh(1      + isig)
          else
            dxh(0        + isig) = dx(0      + isig)/2.
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)/2.
          endif

        enddo

      end subroutine createLogicalGrid

c     consistencyCheck
c     #################################################################
      subroutine consistencyCheck

c     -----------------------------------------------------------------
c     Checks consistency of grid parameters
c     -----------------------------------------------------------------

        implicit none

c     Input variables

c     Local variables

        real(8) :: major_r

c     Begin program

        !Default sizes
        if (xmax == 0d0) then
          xmax = 2*pi
          xmin = 0d0
        endif
        if (ymax == 0d0) then
          ymax = 2*pi
          ymin = 0d0
        endif
        if (zmax == 0d0) then
          zmax = 2*pi
          zmin = 0d0
        endif

        !Consistency
        select case (coords)
        case ('car')
        case ('scl')
        case ('cyl')
          if (xmin /= 0d0 .and. bcond(1) == SP) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif
        case ('hel')
          if (xmin /= 0d0 .and. bcond(1) == SP) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif
        case ('tor')
          if (xmin /= 0d0 .and. bcond(1) == SP) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif

          major_r = grid_params%params(1)

          if (major_r < xmax) then
            write (*,*) 'Ill-defined toroidal coordinate system'
            write (*,*) 'Major radius < minor radius'
            write (*,*) 'Aborting'
            stop
          endif

        end select

        !Ensure ignorable directions are small for numerical computation
        !  of grid parameters
        if (numerical_grid) then
          if (nxx == 1) then
            xmin = 0d0
            xmax = 1d-3
          endif
          if (nyy == 1) then
            ymin = 0d0
            ymax = 1d-3
          endif
          if (nzz == 1) then
            zmin = 0d0
            zmax = 1d-3
          endif
        endif

      end subroutine consistencyCheck

c     checkGrid
c     #################################################################
      subroutine checkGrid

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igx,igy,igz

c     Local variables

        integer(4) :: i,j,n1,n2,n3,isig

c     Begin program

c     Multigrid parameters

        igx = grid_params%ngrid
        igy = grid_params%ngrid
        igz = grid_params%ngrid

        write (*,*)
        write (*,*) 'Coordinate system: ',coords
        write (*,*)
        write (*,*) 'Number of grid levels:'
        write (*,*) 'nx: ',igx,'   ny: ',igy,'   nz: ',igz

        call gridInfo('X',igx,grid_params%istartx,grid_params%nxv
     .               ,grid_params%xx,grid_params%dx,grid_params%dxh)

        call gridInfo('Y',igy,grid_params%istarty,grid_params%nyv
     .               ,grid_params%yy,grid_params%dy,grid_params%dyh)

        call gridInfo('Z',igz,grid_params%istartz,grid_params%nzv
     .               ,grid_params%zz,grid_params%dz,grid_params%dzh)

        stop

c     End program

      contains

c     gridInfo
c     #################################################################  
      subroutine gridInfo(char,ig,istart,nxv,xx,dx,dxh)

        implicit  none

c     Call variables

        character*(1) :: char
        integer(4) :: ig,istart(*),nxv(*)
        real(8)    :: xx(*),dx(*),dxh(*)

c     Local variables

        integer(4) :: i,j,isig

c     Begin program

        write (*,*)
        write (*,*) '***************************'
        write (*,*) 'MG grid in ',char,'-axis'
        write (*,*) '***************************'
        do i = 1,ig
          isig = istart(i)
          write (*,*)
          write (*,*) '************* Grid ',i,' **************'
          write (*,*)
          write (*,*) 'Size ',nxv(i)
          write (*,*) 'MG pointer: ',isig
          write (*,*) 'Grid nodes'
          do j = isig,isig+nxv(i)+1
            write (*,10) 'node :',j-isig
     .                  ,'   position: ',xx(j)
     .                  ,'   int dh: ',dx(j)
     .                  ,'   half dh: ',dxh(j)
 10         format (a,i3,a,f6.3,a,f6.3,a,f6.3)
          enddo
        enddo

c     End program

      end subroutine gridInfo

      end subroutine checkGrid

      end module grid


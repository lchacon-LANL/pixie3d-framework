c module bc_def
c #####################################################################
      module bc_def

        implicit none


        integer(4) :: PER,DIR,NEU,SP,EQU,DEF
        parameter (EQU=1,PER=2,NEU=3,DIR=4,SP=5,DEF=6)

        integer(4) :: bcond(6)

      end module bc_def

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
          real(8)   ,pointer,dimension(:)  :: xx          !Grid node positions in X (all grids)
          real(8)   ,pointer,dimension(:)  :: yy          !Grid node positions in Y (")
          real(8)   ,pointer,dimension(:)  :: zz          !Grid node positions in Z (")
          real(8)   ,pointer,dimension(:)  :: dx          !Grid spacings in X for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dy          !Grid spacings in Y for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dz          !Grid spacings in Z for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dxh         !Grid spacings in X for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dyh         !Grid spacings in Y for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dzh         !Grid spacings in Z for half mesh (")
          integer(4),pointer,dimension(:)  :: nxv         !# of grid nodes in X  (")
          integer(4),pointer,dimension(:)  :: nyv         !# of grid nodes in Y  (")
          integer(4),pointer,dimension(:)  :: nzv         !# of grid nodes in Z  (")
          integer(4),pointer,dimension(:)  :: ntotv       !Total # of grid nodes (")
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
          allocate(grid_st%ntotv(ngrid))
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
          deallocate(grid_st%ntotv)
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

        type(grid_def),intent(in)  :: grid_st1
        type(grid_def),intent(out) :: grid_st2

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

        grid_st2%iline      = grid_st1%iline
        grid_st2%jline      = grid_st1%jline
        grid_st2%kline      = grid_st1%kline
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
        grid_st2%ntotv      = grid_st1%ntotv     
        grid_st2%istartx    = grid_st1%istartx   
        grid_st2%istarty    = grid_st1%istarty   
        grid_st2%istartz    = grid_st1%istartz   
        grid_st2%istartp    = grid_st1%istartp   
        grid_st2%mg_ratio_x = grid_st1%mg_ratio_x
        grid_st2%mg_ratio_y = grid_st1%mg_ratio_y
        grid_st2%mg_ratio_z = grid_st1%mg_ratio_z

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
        write (*,*) 'ntotv',grid_st%ntotv     
        write (*,*) 'istartx',grid_st%istartx   
        write (*,*) 'istarty',grid_st%istarty   
        write (*,*) 'istartz',grid_st%istartz   
        write (*,*) 'istartp',grid_st%istartp   
        write (*,*) 'mg_ratio_x',grid_st%mg_ratio_x
        write (*,*) 'mg_ratio_y',grid_st%mg_ratio_y
        write (*,*) 'mg_ratio_z',grid_st%mg_ratio_z

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

      end module grid_structures

c module grid_definition
c #####################################################################
      module grid_definition

        use grid_structures

        implicit none

        real(8)         :: gparams(5)

        character*(3)   :: coords

        real(8)         :: xmax,ymax,zmax,xmin,ymin,zmin  !3D grid dimension

        real(8),private :: pi,lambda,cc,ypp,eps,mm,kk,aa,phi,major_r

        logical         :: numerical_grid,anal_map

      contains

c     checkGridDatabase
c     #################################################################
      function checkGridDatabase() result(anal_map)

c     -----------------------------------------------------------------
c     Checks grid database for analytical mappings
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        logical :: anal_map

c     Begin program

        select case(coords)
        case('car','scl','cyl','hel','tor','sin')

          anal_map = .true.

        case default

          anal_map = .false.

        end select

      end function checkGridDatabase

c     getCartesianCoordinates
c     #################################################################
      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Inverts curvilinear coordinates to give Cartesian coordinates.
c     Requires external routine 'map', with call sequence:
c
c             map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1)
c
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        real(8)    :: car(3)

c     Externals

        external   :: map

c     Begin program

        if (checkGridDatabase()) then

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

          car = x_xi(x1,y1,z1)

          x1 = car(1)
          y1 = car(2)
          z1 = car(3)

        else

          call map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1)

        endif

      end subroutine getCartesianCoordinates

c     x_xi
c     #################################################################
      function x_xi(x1,x2,x3) result(car)

c     -----------------------------------------------------------------
c     Gives Cartesian coordinates from curvilinear coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x1,x2,x3,car(3)

c     Local variables

        integer(4) :: inewt,ic
        real(8)    :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0

c     Begin program

        select case (coords)
        case ('car')
          xx = x1
          yy = x2
          zz = x3
        case ('scl')
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*x2/ymax-1.)

          xx = x1
          yy = 0.5+lambda*atanh(ypp/cc)
          zz = x3
        case ('cyl')
          xx = x1*cos(x2)
          yy = x1*sin(x2)
          zz = x3
        case ('hel')
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (x2-aa*x3)

          xx = x1*cos(phi)
          yy = x1*sin(phi)
          zz = x3
        case ('tor')
          major_r = gparams(1)

          xx = (major_r + x1*sin(x2))*cos(x3)
          yy = (major_r + x1*sin(x2))*sin(x3)
          zz = x1*cos(x2)
        case ('sin')
          pi = acos(-1d0)
          eps = gparams(1)

          xx = x1 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
          yy = x2 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
          zz = x3
        case default
          write (*,*) 'Grid not implemented in x_xi'
          write (*,*) 'Aborting...'
          stop
        end select

        car = (/ xx,yy,zz /)

      contains

c     atanh
c     #################################################################
      real(8) function atanh(x)

        real(8) :: x

        atanh = 0.5*(log( (1+x)/(1-x) ) )

      end function atanh

      end function x_xi

c     jacobian
c     #################################################################
      function jacobian(i,j,k,igx,igy,igz) result(jac)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: jac

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,car(3),curv(3)

c     Begin program

        select case (coords)
        case ('car')
          jac = 1d0
        case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)

          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          jac = x1
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          jac = x1
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          major_r = gparams(1)

          jac = x1*(major_r + x1*sin(x2))
        case ('sin')
          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)
          jac = (xmax*ymax + eps*Pi*(xmax - ymax)*
     -          Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -          eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
     -          /(xmax*ymax)
        case default
          write (*,*) 'Grid not implemented in jacobian'
          write (*,*) 'Aborting...'
          stop
        end select

      end function jacobian

c     covariantVector2
c     #################################################################
      function covariantVector(comp,i,j,k,igx,igy,igz) result (vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: comp,i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,vec(3)

        real(8)    :: car(3),curv(3),jac

c     Begin program

        select case (coords)
        case ('car')
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('scl')
          lambda = gparams(1)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          curv = (/ x1,x2,x3 /)

          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          jac = 1./jac
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,jac,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))

          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)
            case (2)
              vec = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('tor')
          major_r = gparams(1)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ sin(curv(2))*cos(curv(3))
     .                ,sin(curv(2))*sin(curv(3))
     .                ,cos(curv(2)) /)
            case (2)
              vec = (/ cos(curv(2))*cos(curv(3))
     .                ,cos(curv(2))*sin(curv(3))
     .                ,-sin(curv(2))/)/curv(1)
            case (3)
              vec = (/ -sin(curv(3)),cos(curv(3)),0d0 /)
     .              /(major_r + curv(1)*sin(curv(2)))
          end select
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

          select case (comp)
            case (1)
              vec = (/ (xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)*
     -                                        Sin((2*Pi*x1)/xmax)))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)*
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , (-2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
     .                               *Sin((2*Pi*x1)/xmax))/
     -                (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))  
     -                           + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , 0d0 /)
            case (2)
              vec = (/ (-2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
     .                                *Sin((2*Pi*x2)/ymax))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)
     .                                      *Sin((2*Pi*x2)/ymax)))/
     -                (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                           + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case default
          write (*,*) 'Grid not implemented in covariantVector'
          write (*,*) 'Aborting...'
          stop
        end select

      end function covariantVector

c     contravariantVector2
c     #################################################################
      function contravariantVector(comp,i,j,k,igx,igy,igz) result (vec)

c     -----------------------------------------------------------------
c     Calculates contravariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: comp,i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: car(3),curv(3),jac
        real(8)    :: x1,x2,x3,vec(3)

c     Begin program

        select case (coords)
        case ('car')
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
       case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          jac = 1./jac
          select case (comp)
            case (1)
              vec = (/ jac,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,jac /)
          end select
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)/curv(1)
          end select
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(phi),cos(phi),0d0 /)
            case (3)
              vec = (/ aa*sin(phi),-aa*cos(phi),1d0/curv(1) /)
          end select
        case ('tor')
          major_r = gparams(1)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ sin(curv(2))*cos(curv(3))
     .                ,sin(curv(2))*sin(curv(3))
     .                ,cos(curv(2)) /)
     .              /curv(1)/(major_r + curv(1)*sin(curv(2)))
            case (2)
              vec = (/ cos(curv(2))*cos(curv(3))
     .                ,cos(curv(2))*sin(curv(3))
     .                ,-sin(curv(2))/)
     .              /(major_r + curv(1)*sin(curv(2)))
            case (3)
              vec = (/ -sin(curv(3)),cos(curv(3)),0d0 /)/curv(1)
          end select
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

          select case (comp)
            case (1)
              vec = (/ (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)*
     -                                        Sin((2*Pi*x2)/ymax)))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , (2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
     .                              *Sin((2*Pi*x2)/ymax))/
     -                (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                           + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,0d0 /)
            case (2)
              vec = (/ (2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
     .                               *Sin((2*Pi*x1)/xmax))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,(xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)
     .                                     *Sin((2*Pi*x1)/xmax)))/
     -               (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                          + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,0d0 /)
            case (3)
              vec = (/ 0d0
     .                ,0d0
     .                ,(xmax*ymax)/(xmax*ymax
     .                  + eps*Pi*(xmax - ymax)
     -                       *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                  + eps*Pi*(xmax + ymax)
     .                       *Sin(2*Pi*(x1/xmax + x2/ymax))) /)
          end select
        case default
          write (*,*) 'Grid not implemented in contravariantVector'
          write (*,*) 'Aborting...'
          stop
        end select

      end function contravariantVector

c     hessian22
c     #################################################################
      function hessian22(l,i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates hessian elements of curvilinear coordinate system in
c     the covariant basis, i.e., 
c            hessian[l](i,j) = J^2 <cnv(i)|grad(cov[l])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: l,i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: car(3),curv(3),vec(9)
        real(8)    :: x1,x2,x3,tensor(3,3)

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)

          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 2.*ypp/(ypp**2-cc**2)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)
              vec(6) = -aa*curv(1)
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = aa**2*curv(1)
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = aa/curv(1)
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gparams(1)
          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = sin(curv(2))*(major_r + curv(1)*sin(curv(2)))
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = (major_r + curv(1)*sin(curv(2)))*cos(curv(2))
     .                 /curv(1)
            case (3)
              vec(1) =  0d0
              vec(2) =  0d0
              vec(3) =  -sin(curv(2))/(major_r + curv(1)*Sin(curv(2)))
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = -curv(1)*cos(curv(2))
     .                 /(major_r + curv(1)*Sin(curv(2)))
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          eps = gparams(1)

          select case (l)
            case (1)
              vec(1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)*
     -                                    Sin((2*Pi*x2)/ymax))/
     -          (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -           Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -           eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
     .                             *Cos((2*Pi*x2)/ymax))/
     -            (xmax*ymax + eps*Pi*(xmax - ymax)*
     -              Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
     .                                   *Sin((2*Pi*x2)/ymax))/
     -            (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -            Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)
     -                                   *Sin((2*Pi*x2)/ymax))/
     -                 (xmax*(xmax*ymax
     .                 + eps*Pi*(xmax - ymax)
     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                 + eps*Pi*(xmax + ymax)
     .                     *Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
     .                             *Cos((2*Pi*x2)/ymax))/
     -                 (xmax*ymax
     .                 + eps*Pi*(xmax - ymax)
     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                 + eps*Pi*(xmax + ymax)
     .                     *Sin(2*Pi*(x1/xmax + x2/ymax)))
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
     .                                   *Sin((2*Pi*x2)/ymax))/
     -                  (ymax*(xmax*ymax
     .                  + eps*Pi*(xmax - ymax)
     -                      *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                  + eps*Pi*(xmax + ymax)
     .                      *Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

        case default
          write (*,*) 'Grid not implemented in hessian'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec,(/3,3/))

      end function hessian22

c     christ_2knd
c     #################################################################
      function christ_2knd(i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates elements of Christoffel symbol of the second kind,
c     Gamma[l](i,j) defined as: 
c            gamma[l](i,j) = -J^2 <cnv(i)|grad(cov[l])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: l,i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3,3)

c     Local variables

c     Begin program

        tensor(1,:,:) = hessian22(1,i,j,k,igx,igy,igz)
        tensor(2,:,:) = hessian22(2,i,j,k,igx,igy,igz)
        tensor(3,:,:) = hessian22(3,i,j,k,igx,igy,igz)

c     End program

      end function christ_2knd

c     hessian_cnv
c     #################################################################
      function hessian_cnv(k,x1,x2,x3) result (tensor)

c     -----------------------------------------------------------------
c     Calculates elements of tensor grad(cnv) in a mixed coordinate 
c     system:
c              hessian_cnv[k](i,j) = J^2 <cnv(i)|grad(cnv[k])|cov(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: k
        real(8)    :: x1,x2,x3,tensor(3,3)

c     Local variables

        real(8)    :: vec(9),car(3),eps

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*x2/ymax-1.)
          select case (k)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 2.*ypp/(ypp**2-cc**2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 2.*ypp/(ypp**2-cc**2)
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
          end select
        case ('cyl')
          select case (k)
            case (1)
              vec(1) = -1./x1
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 0d0
              vec(5) = 1./x1
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = -x1
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = -1./x1
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
          end select
        case ('hel')
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm 
          select case (k)
          case (1)
            vec(1) = -1./x1
            vec(2) = 0d0
            vec(3) = 0d0
            vec(4) = 0d0
            vec(5) = 1./x1
            vec(6) = 0d0
            vec(7) = 0d0
            vec(8) = -aa/x1
            vec(9) = 0d0
          case (2)
            vec(1) = 0d0
            vec(2) = 0d0
            vec(3) = 0d0
            vec(4) = -x1
            vec(5) = 0d0
            vec(6) = 0d0
            vec(7) = aa*x1
            vec(8) = 0d0
            vec(9) = 0d0
          case (3)
            vec(1) = 0d0
            vec(2) = -aa/x1
            vec(3) = -1./x1
            vec(4) = aa*x1
            vec(5) = 0d0
            vec(6) = 0d0
            vec(7) = -aa**2*x1
            vec(8) = 0d0
            vec(9) = 0d0
          end select
        case ('tor')
          major_r = gparams(1)
          select case (k)
            case (1)
              vec(1) = -(major_r + 2*x1*sin(x2))
     .                  /x1/(major_r + x1*sin(x2))
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = -x1*cos(x2)/(major_r + x1*sin(x2))
              vec(5) = 1./x1
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = sin(x2)/(major_r + x1*sin(x2))
            case (2)
              vec(1) = 0d0
              vec(2) = -sin(x2)/(major_r + x1*Sin(x2))
              vec(3) = 0d0
              vec(4) = -x1
              vec(5) = -x1*cos(x2)/(major_r + x1*Sin(x2))
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) =  x1*cos(x2)/(major_r + x1*Sin(x2))
            case (3)
              vec(1) =  0d0
              vec(2) =  0d0
              vec(3) =  -1./x1
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = -sin(x2)*(major_r + x1*Sin(x2))
              vec(8) = -cos(x2)*(major_r + x1*Sin(x2))/x1
              vec(9) = 0d0
          end select

        case default
          write (*,*) 'Grid not implemented in hessian_cnv'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = transpose(reshape(vec, (/3,3/)))

      end function hessian_cnv

c     g_sub
c     #################################################################
      function g_sub(i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3)

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: vec(9),jac,car(3),curv(3)
        real(8)    :: x1,x2,x3

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
        case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          vec = (/ 1d0/jac, 0d0, 0d0
     .            ,0d0    , jac, 0d0
     .            ,0d0    , 0d0, 1d0/jac /)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          vec = (/ 1d0/curv(1), 0d0     , 0d0
     .            ,0d0        , curv(1) , 0d0
     .            ,0d0        , 0d0     , 1d0/curv(1) /)
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1./curv(1) + aa**2*curv(1) /)
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gparams(1)
          vec = (/ 1d0/curv(1)/(major_r + curv(1)*sin(curv(2))),0d0,0d0
     .            ,0d0, curv(1)/(major_r + curv(1)*sin(curv(2))), 0d0
     .            ,0d0, 0d0, (major_r/curv(1) + sin(curv(2))) /)
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

          vec(1) =
     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(2) =
     .         (eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -        xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
     -        ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
     -        ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))

          vec(3) = 0d0
          vec(4) = vec(2)

          vec(5) =
     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))!

          vec(6) = 0d0
          vec(7) = vec(3)
          vec(8) = vec(6)

          vec(9) = (xmax*ymax)/
     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))

        case default
          write (*,*) 'Grid not implemented in g_sub'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec, (/3,3/))

      end function g_sub

c     g_sup
c     #################################################################
      function g_sup(i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates covariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3)

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: vec(9),jac,car(3),curv(3)
        real(8)    :: x1,x2,x3

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
       case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          vec = (/ jac, 0d0    , 0d0
     .            ,0d0, 1d0/jac, 0d0
     .            ,0d0, 0d0    , jac /)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          vec = (/ curv(1) , 0d0   , 0d0
     .            ,0d0, 1d0/curv(1), 0d0
     .            ,0d0, 0d0   , curv(1) /)
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2*curv(1), aa*curv(1)
     .            ,0d0, aa*curv(1), curv(1) /)
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gparams(1)
          vec = (/ curv(1)*(major_r + curv(1)*sin(curv(2))), 0d0, 0d0
     .            ,0d0, (major_r/curv(1) + sin(curv(2))), 0d0
     .            ,0d0, 0d0, curv(1)/(major_r + curv(1)*sin(curv(2))) /)
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

          vec(1) =
     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(2) =
     .         -((eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -          eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -          xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
     -          ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -          xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
     -          ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -      (xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(3) = 0d0
          vec(4) = vec(2)

          vec(5) =
     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(6) = 0d0
          vec(7) = vec(3)
          vec(8) = vec(6)

          vec(9) = (xmax*ymax + eps*Pi*(xmax - ymax)*
     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
     .         /(xmax*ymax)

        case default
          write (*,*) 'Grid not implemented in g_sup'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec, (/3,3/))

      end function g_sup

      end module grid_definition

c module grid_metric
c #####################################################################
      module grid_metric

c ---------------------------------------------------------------------
c     This module packs routines that perform operations on grid
c     quantities, such as coordinate transformation of vector components,
c     vector norms and scalar products. It contains the following
c     routines:
c          * transformVectorToCartesian 
c          * transformVectorToCurvilinear 
c          * transformFromCurvToCurv 
c          * volume
c          * vectorNorm
c          * scalarProduct
c     It is assumed that the grid metric structure gmetric is allocated 
c     and filled.
c ---------------------------------------------------------------------

        use grid_definition

        use bc_def

        implicit none

        type :: grid_metrics
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

        integer(4)      :: igrid,nxp,nyp,nzp

c     Begin program

        if (.not.associated(gmetric%grid)) then
          allocate(gmetric%grid(grid_params%ngrid))
        endif

        do igrid=1,grid_params%ngrid
          if (.not.associated(gmetric%grid(igrid)%jac)) then
            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1
            allocate(gmetric%grid(igrid)%jac  (0:nxp,0:nyp,0:nzp))
            allocate(gmetric%grid(igrid)%gsub (0:nxp,0:nyp,0:nzp,3,3))
            allocate(gmetric%grid(igrid)%gsup (0:nxp,0:nyp,0:nzp,3,3))
            allocate(gmetric%grid(igrid)%cov  (0:nxp,0:nyp,0:nzp,3,3))
            allocate(gmetric%grid(igrid)%cnv  (0:nxp,0:nyp,0:nzp,3,3))
            allocate(gmetric%grid(igrid)%Gamma(0:nxp,0:nyp,0:nzp,3,3,3))
          endif
        enddo

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

        do igrid=1,grid_params%ngrid
          if (associated(gmetric%grid(igrid)%jac)) then
            deallocate(gmetric%grid(igrid)%jac  )
            deallocate(gmetric%grid(igrid)%gsub )
            deallocate(gmetric%grid(igrid)%gsup )
            deallocate(gmetric%grid(igrid)%cov  )
            deallocate(gmetric%grid(igrid)%cnv  )
            deallocate(gmetric%grid(igrid)%Gamma)
          endif
        enddo

        if (associated(gmetric%grid)) then
          deallocate(gmetric%grid)
        endif

c     End program

      end subroutine deallocateGridMetric

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric(gmetric)

c     -----------------------------------------------------------------
c     This routine calculates all grid metric quantities required for
c     the curvilinear representation of a set of PDE's: jacobian,
c     metric tensors, covariant and contravariant vectors, Christoffel
c     symbols of the second kind. All quantities are stored in
c     structure gmetric. There are two modes of computation:
c        * Analytical (numerical_grid=.false.)
c        * Numerical  (numerical_grid=.true.)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4) :: igrid,nxp,nyp,nzp,i,j,k,igx,igy,igz
     .               ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .               ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp
        real(8)    :: r(3,3),car0(3),carp(3),carm(3),dh(3),jac,ijac
     .               ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .               ,gamma(3,3,3),gamm1(3,3,3),mag,dhp,dhm,dhh
        real(8),allocatable,dimension(:,:,:,:,:) :: dr

c     Interpolation

        integer(4) :: kx,ky,kz,nnx,nny,nnz,dim,flg,order

        real(8)    :: xp,yp,zp
        real(8),allocatable,dimension(:) :: sx,sy,sz
        real(8), dimension(:),allocatable:: tx,ty,tz,work
        real(8), dimension(:,:,:,:,:),allocatable:: drbcoef
        real(8), dimension(:,:,:,:,:,:),allocatable:: gambcoef

        real(8)    :: db3val
        external   :: db3val

c     Begin program

        anal_map = checkGridDatabase()

        call allocateGridMetric(gmetric)

        if ((.not.numerical_grid).and.anal_map) then

c       Find analytical geometric quantitites on grid

          do igrid=1,grid_params%ngrid

            igx = igrid
            igy = igrid
            igz = igrid

            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1

            do k = 0,nzp
              do j = 0,nyp
                do i = 0,nxp
                  gmetric%grid(igrid)%jac  (i,j,k)
     .                      = jacobian(i,j,k,igx,igy,igz)
                  gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                      = g_sub   (i,j,k,igx,igy,igz)
                  gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                      = g_sup   (i,j,k,igx,igy,igz)
                  gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .                      = christ_2knd(i,j,k,igx,igy,igz)
                  do l=1,3
                    gmetric%grid(igrid)%cov(i,j,k,l,:)
     .                      = covariantVector    (l,i,j,k,igx,igy,igz)
                    gmetric%grid(igrid)%cnv(i,j,k,l,:)
     .                      = contravariantVector(l,i,j,k,igx,igy,igz)
                  enddo
                enddo
              enddo
            enddo

            !Zero force condition on Christoffel symbols (only on finest grid)
cc            if (igrid == 1) then
cc              do k = 1,grid_params%nzv(igrid)
cc                do j = 1,grid_params%nyv(igrid)
cc                  do i = 1,grid_params%nxv(igrid)
cc                    call gammaZeroForce(i,j,k,igx,igy,igz)
cc                  enddo
cc                enddo
cc              enddo
cc            endif

          enddo

        else

c       Find numerical geometric quantities on finest grid

cc          igrid=1

          do igrid=1,grid_params%ngrid

            igx = igrid
            igy = igrid
            igz = igrid

            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1

            allocate(dr(0:nxp,0:nyp,0:nzp,3,3))

            !Evaluate dx/dxi vectors
            do k = 0,nzp
              do j = 0,nyp
                do i = 0,nxp

                  ip=min(i+1,nxp)
                  im=max(i-1,0)
                  jp=min(j+1,nyp)
                  jm=max(j-1,0)
                  kp=min(k+1,nzp)
                  km=max(k-1,0)

                  carp = map(ip,j,k,igx,igy,igz,igp,jg,kg)
                  carm = map(im,j,k,igx,igy,igz,igm,jg,kg)
                  dh(1)= (grid_params%xx(igp)-grid_params%xx(igm))
                  dr(i,j,k,1,:) = (carp-carm)/dh(1)

                  carp = map(i,jp,k,igx,igy,igz,ig,jgp,kg)
                  carm = map(i,jm,k,igx,igy,igz,ig,jgm,kg)
                  dh(2)= (grid_params%yy(jgp)-grid_params%yy(jgm))
                  dr(i,j,k,2,:) = (carp-carm)/dh(2)

                  carp = map(i,j,kp,igx,igy,igz,ig,jg,kgp)
                  carm = map(i,j,km,igx,igy,igz,ig,jg,kgm)
                  dh(3)= (grid_params%zz(kgp)-grid_params%zz(kgm))
                  dr(i,j,k,3,:) = (carp-carm)/dh(3)

                enddo
              enddo
            enddo

            !Enforce topological constraints on dr
            do i=1,3
              do j=1,3
                call topol_bc(dr(:,:,:,i,j))
              enddo
            enddo

          !Spline dr for interpolation on coarser grids
cc          order = 2
cc
cc          nnx = nxp+1
cc          nny = nyp+1
cc          nnz = nzp+1
cc
cc          allocate(sx(nnx),sy(nny),sz(nnz))
cc
cc          call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)
cc
cc          sx(1:nnx) = grid_params%xx(ig-1:ig+nxp-1)
cc          sy(1:nny) = grid_params%yy(jg-1:jg+nyp-1)
cc          sz(1:nnz) = grid_params%zz(kg-1:kg+nzp-1)
cc
cc          flg = 0
cc
cc          kx = min(order+1,nnx-1)
cc          ky = min(order+1,nny-1)
cc          kz = min(order+1,nnz-1)
cc
cc          dim = nnx*nny*nnz+max(2*kx*(nnx+1),2*ky*(nny+1),2*kz*(nnz+1))
cc
cc          allocate(tx(nnx+kx))
cc          allocate(ty(nny+ky))
cc          allocate(tz(nnz+kz))
cc          allocate(work(dim))
cc          allocate(drbcoef(nnx,nny,nnz,3,3))
cc
cc          do i=1,3
cc            do j=1,3
cc              call db3ink(sx,nnx,sy,nny,sz,nnz,dr(:,:,:,i,j),nnx,nny
cc     .                   ,kx,ky,kz,tx,ty,tz,drbcoef(:,:,:,i,j)
cc     .                   ,work,flg)
cc            enddo
cc          enddo

          !Evaluate grid quantities
            do k = 0,nzp
              do j = 0,nyp
                do i = 0,nxp

                  r = dr(i,j,k,:,:)

                  !Evaluate Jacobian
                  jac = triple_product(r(1,:),r(2,:),r(3,:))
                  ijac = 1d0/jac

                  !Contravariant vectors
                  cnv(:,:) = r(:,:)*ijac

                  !Covariant vectors
                  cov(1,:) = cross_product(r(2,:),r(3,:))*ijac
                  cov(2,:) = cross_product(r(3,:),r(1,:))*ijac
                  cov(3,:) = cross_product(r(1,:),r(2,:))*ijac

                  !Metric tensors
                  do m=1,3
                    do l=m,3
                      gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
                      gsub(m,l) = gsub(l,m) !Symmetry
                      gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
                      gsup(m,l) = gsup(l,m) !Symmetry
                    enddo
                  enddo

                  !Grid spacings
                  call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
                  dh(1) = grid_params%dxh(ig)
                  dh(2) = grid_params%dyh(jg)
                  dh(3) = grid_params%dzh(kg)

                  !Christoffel symbols

                  if (i==0) then
                    i0=1
                  elseif (i==nxp) then
                    i0=nxp-1
                  else
                    i0=i
                  endif

                  if (j==0) then
                    j0=0
                  elseif (j==nyp) then
                    j0=nyp
                  else
                    j0=j
                  endif

                  if (k==0) then
                    k0=1
                  elseif (k==nzp) then
                    k0=nzp-1
                  else
                    k0=k
                  endif

                  ip=i0+1
                  im=i0-1
                  jp=j0+1
                  jm=j0-1
                  kp=k0+1
                  km=k0-1

                  xp = grid_params%xx(ig)
                  yp = grid_params%yy(jg)
                  zp = grid_params%zz(kg)

                  do l=1,3
                    do m=1,3
                      do n=m,3
                        if (m == n) then
                          select case(m)
                          case(1)
                            car0 = map(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
                            carp = map(ip,j0,k0,igx,igy,igz,igp,jg,kg)
                            carm = map(im,j0,k0,igx,igy,igz,igm,jg,kg)
                            dhp = (grid_params%xx(igp)
     .                            -grid_params%xx(ig0))
                            dhm = (grid_params%xx(ig0)
     .                            -grid_params%xx(igm))
                          case(2)
                            car0 = map(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                            carp = map(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                            carm = map(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                            dhp = (grid_params%yy(jgp)
     .                            -grid_params%yy(jg0))
                            dhm = (grid_params%yy(jg0)
     .                            -grid_params%yy(jgm))
                          case(3)
                            car0 = map(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                            carp = map(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                            carm = map(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                            dhp = (grid_params%zz(kgp)
     .                            -grid_params%zz(kg0))
                            dhm = (grid_params%zz(kg0)
     .                            -grid_params%zz(kgm))
                          end select
                          carp = (carp-car0)/dhp
                          carm = (car0-carm)/dhm
                          dhh  = 0.5*(dhp+dhm)
                        else
                          select case(m)
                          case(1)
                           call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
                           call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
                           call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
                           carp= 0.5*(dr(ip,j0,k0,n,:)+dr(i0,j0,k0,n,:))
                           carm= 0.5*(dr(im,j0,k0,n,:)+dr(i0,j0,k0,n,:))
                           dhh = 0.5*(grid_params%xx(igp)
     .                               -grid_params%xx(igm))
                          case(2)
                           call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                           call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                           call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                           carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
                           carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
                           dhh = 0.5*(grid_params%yy(jgp)
     .                               -grid_params%yy(jgm))
                          case(3)
                           call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                           call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                           call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                           carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
                           carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
                           dhh = 0.5*(grid_params%zz(kgp)
     .                               -grid_params%zz(kgm))
                          end select
                        endif
                        vec = (carp-carm)/dhh

cc                      select case(m)
cc                      case(1)
cc                        do p=1,3
cc                          vec(p)=db3val(xp,yp,zp,1,0,0,tx,ty,tz
cc     .                          ,nnx,nny,nnz,kx,ky,kz,drbcoef(:,:,:,n,p)
cc     .                          ,work)
cc                        enddo
cc                      case(2)
cc                        do p=1,3
cc                          vec(p)=db3val(xp,yp,zp,0,1,0,tx,ty,tz
cc     .                          ,nnx,nny,nnz,kx,ky,kz,drbcoef(:,:,:,n,p)
cc     .                          ,work)
cc                        enddo
cc                      case(3)
cc                        do p=1,3
cc                          vec(p)=db3val(xp,yp,zp,0,0,1,tx,ty,tz
cc     .                          ,nnx,nny,nnz,kx,ky,kz,drbcoef(:,:,:,n,p)
cc     .                          ,work)
cc                        enddo
cc                      end select

                        gamma(l,m,n) = dot_product(vec,cov(l,:))
                        gamma(l,n,m) = gamma(l,m,n) !Symmetry
                      enddo
                    enddo
                  enddo

cc                  gamma = christ_2knd(i,j,k,igx,igy,igz)

                  !Store grid quantities
                  gmetric%grid(igrid)%jac  (i,j,k)       = jac
                  gmetric%grid(igrid)%gsub (i,j,k,:,:)   = gsub
                  gmetric%grid(igrid)%gsup (i,j,k,:,:)   = gsup
                  gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
                  gmetric%grid(igrid)%cov  (i,j,k,:,:)   = cov
                  gmetric%grid(igrid)%cnv  (i,j,k,:,:)   = cnv

                enddo
              enddo
            enddo

            deallocate(dr)

            !Zero-force condition on Christoffle symbols 
            if (igrid == 1) then
              do k = 1,grid_params%nzv(igrid)
                do j = 1,grid_params%nyv(igrid)
                  do i = 1,grid_params%nxv(igrid)
                    call gammaZeroForce(i,j,k,igx,igy,igz)
                  enddo
                enddo
              enddo
            endif

            !Enforce topological constraints on Christoffel symbols
            do i=1,3
              do j=1,3
                do k=1,3
                  call topol_bc(gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
                enddo
              enddo
            enddo

          !Spline Christoffel symbols for interpolation on coarser grids
cc          allocate(gambcoef(nnx,nny,nnz,3,3,3))
cc          flg = 0
cc
cc          do k=1,3
cc            do j=1,3
cc              do i=1,3
cc                call db3ink(sx,nnx,sy,nny,sz,nnz
cc     .                     ,gmetric%grid(igrid)%Gamma(:,:,:,i,j,k)
cc     .                     ,nnx,nny,kx,ky,kz,tx,ty,tz
cc     .                     ,gambcoef(:,:,:,i,j,k),work,flg)
cc              enddo
cc            enddo
cc          enddo

c       Restrict grid metric info to coarser grids

cc          call restrictGridMetrics

c       Deallocate spline auxiliary variables

cc          deallocate(sx,sy,sz,tx,ty,tz,work,drbcoef,gambcoef)

          enddo

c       Test of numerical calculation vs. analytical calculation

cc          do igrid=1,grid_params%ngrid
cc
cc            igx = igrid
cc            igy = igrid
cc            igz = igrid
cc
cc            nxp = grid_params%nxv(igrid)+1
cc            nyp = grid_params%nyv(igrid)+1
cc            nzp = grid_params%nzv(igrid)+1
cc
cc            do k = 0,nzp
cc              do j = 0,nyp
cc                do i = 0,nxp
cc                  gmetric%grid(igrid)%jac  (i,j,k)
cc     .                      = gmetric%grid(igrid)%jac  (i,j,k)
cc     .                       -jacobian(i,j,k,igx,igy,igz)
cc                  gmetric%grid(igrid)%gsub (i,j,k,:,:)
cc     .                      = gmetric%grid(igrid)%gsub (i,j,k,:,:)
cc     .                       -g_sub   (i,j,k,igx,igy,igz)
cc                  gmetric%grid(igrid)%gsup (i,j,k,:,:)
cc     .                      = gmetric%grid(igrid)%gsup (i,j,k,:,:)
cc     .                       -g_sup   (i,j,k,igx,igy,igz)
cc                  gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc     .                      = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc     .                       -christ_2knd(i,j,k,igx,igy,igz)
cc                enddo
cc              enddo
cc            enddo
cc
cc            write (*,*) 'Grid level:',igrid
cc          
cc            mag =
cc     .     sum(gmetric%grid(igrid)%jac(1:nxp-1,1:nyp-1,1:nzp-1)**2)
cc            write (*,*)'Jacobian tst=',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cc            mag =
cc     .     sum(gmetric%grid(igrid)%gsub(2:nxp-1,1:nyp-1,1:nzp-1,:,:)**2)
cc            write (*,*)'Gsub tst    =',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cc            mag =
cc     .     sum(gmetric%grid(igrid)%gsup(2:nxp-1,1:nyp-1,1:nzp-1,:,:)**2)
cc            write (*,*)'Gsup tst    =',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cc            mag =
cc     .  sum(gmetric%grid(igrid)%Gamma(2:nxp-1,1:nyp-1,1:nzp-1,:,:,:)**2)
cc            write (*,*)'Gamma tst   =',sqrt(mag/(nxp-1)/(nyp-1)/(nzp-1))
cc
cc          enddo
cc          stop

        endif

c     End program

      contains

c     map
c     #################################################################
      function map(i,j,k,igx,igy,igz,ig,jg,kg) result(car)

c     -----------------------------------------------------------------
c     Give Cartesian coordinates corresponding to node (i,j,k) at grid
c     level (igx,igy,igz).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: car(3)

c     Local variables

        real(8)    :: x1,y1,z1

c     Begin program

        call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                              ,x1,y1,z1)

        car = (/ x1,y1,z1 /)

      end function map

c     restrictGridMetrics
c     #################################################################
      subroutine restrictGridMetrics
c     -----------------------------------------------------------------
c     Restricts arrays with geometric grid info.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

c     Local variables

      real (8) :: xp,yp,zp

c     Begin program

      do igrid=2,grid_params%ngrid

        igx = igrid
        igy = igrid
        igz = igrid

        nxp = grid_params%nxv(igrid)+1
        nyp = grid_params%nyv(igrid)+1
        nzp = grid_params%nzv(igrid)+1

        allocate(dr(0:nxp,0:nyp,0:nzp,3,3))

        !Evaluate dx/dxi vectors (splined)
        do k = 0,nzp
          do j = 0,nyp
            do i = 0,nxp
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              xp = grid_params%xx(ig)
              yp = grid_params%yy(jg)
              zp = grid_params%zz(kg)

              do l=1,3
                do m=1,3
                  dr(i,j,k,l,m) =
     .                db3val(xp,yp,zp,0,0,0,tx,ty,tz,nnx,nny,nnz
     .                      ,kx,ky,kz,drbcoef(:,:,:,l,m),work)
                enddo
              enddo

            enddo
          enddo
        enddo

        !Enforce topological constraints on dr
        do i=1,3
          do j=1,3
            call topol_bc(dr(:,:,:,i,j))
          enddo
        enddo

        !Evaluate grid quantities
        do k = 0,nzp
          do j = 0,nyp
            do i = 0,nxp

              r = dr(i,j,k,:,:)

              !Evaluate Jacobian
              jac = triple_product(r(1,:),r(2,:),r(3,:))
              ijac = 1d0/jac

              !Contravariant vectors
              cnv(:,:) = r(:,:)*ijac

              !Covariant vectors
              cov(1,:) = cross_product(r(2,:),r(3,:))*ijac
              cov(2,:) = cross_product(r(3,:),r(1,:))*ijac
              cov(3,:) = cross_product(r(1,:),r(2,:))*ijac

              !Metric tensors
              do m=1,3
                do l=m,3
                  gsub(l,m) = jac*dot_product(cnv(l,:),cnv(m,:))
                  gsub(m,l) = gsub(l,m) !Symmetry
                  gsup(l,m) = jac*dot_product(cov(l,:),cov(m,:))
                  gsup(m,l) = gsup(l,m) !Symmetry
                enddo
              enddo

              !Grid spacings
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
              dh(1) = grid_params%dxh(ig)
              dh(2) = grid_params%dyh(jg)
              dh(3) = grid_params%dzh(kg)

              !Christoffel symbols (splined)

                if (i==0) then
                  i0=1
                elseif (i==nxp) then
                  i0=nxp-1
                else
                  i0=i
                endif

                if (j==0) then
                  j0=0
                elseif (j==nyp) then
                  j0=nyp
                else
                  j0=j
                endif

                if (k==0) then
                  k0=1
                elseif (k==nzp) then
                  k0=nzp-1
                else
                  k0=k
                endif

                ip=i0+1
                im=i0-1
                jp=j0+1
                jm=j0-1
                kp=k0+1
                km=k0-1

              xp = grid_params%xx(ig)
              yp = grid_params%yy(jg)
              zp = grid_params%zz(kg)

              do l=1,3
                do m=1,3
                  do n=m,3
cc                    gamma(l,m,n) =
cc     .                db3val(xp,yp,zp,0,0,0,tx,ty,tz,nnx,nny,nnz
cc     .                      ,kx,ky,kz,gambcoef(:,:,:,l,m,n)
cc     .                      ,work)

                    select case(m)
                    case(1)
                      call getMGmap(i0,j0,k0,igx,igy,igz,ig0,jg,kg)
                      call getMGmap(im,j0,k0,igx,igy,igz,igm,jg,kg)
                      call getMGmap(ip,j0,k0,igx,igy,igz,igp,jg,kg)
                      carp= 0.5*(dr(ip,j0,k0,n,:)+dr(i0,j0,k0,n,:))
                      carm= 0.5*(dr(im,j0,k0,n,:)+dr(i0,j0,k0,n,:))
                      dhh = 0.5*(grid_params%xx(igp)
     .                              -grid_params%xx(igm))
                    case(2)
                      call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg0,kg)
                      call getMGmap(i0,jm,k0,igx,igy,igz,ig,jgm,kg)
                      call getMGmap(i0,jp,k0,igx,igy,igz,ig,jgp,kg)
                      carp= 0.5*(dr(i0,jp,k0,n,:)+dr(i0,j0,k0,n,:))
                      carm= 0.5*(dr(i0,jm,k0,n,:)+dr(i0,j0,k0,n,:))
                      dhh = 0.5*(grid_params%yy(jgp)
     .                              -grid_params%yy(jgm))
                    case(3)
                      call getMGmap(i0,j0,k0,igx,igy,igz,ig,jg,kg0)
                      call getMGmap(i0,j0,km,igx,igy,igz,ig,jg,kgm)
                      call getMGmap(i0,j0,kp,igx,igy,igz,ig,jg,kgp)
                      carp= 0.5*(dr(i0,j0,kp,n,:)+dr(i0,j0,k0,n,:))
                      carm= 0.5*(dr(i0,j0,km,n,:)+dr(i0,j0,k0,n,:))
                      dhh = 0.5*(grid_params%zz(kgp)
     .                              -grid_params%zz(kgm))
                    end select
                    vec = (carp-carm)/dhh

                    gamma(l,m,n) = dot_product(vec,cov(l,:))
                    gamma(l,n,m) = gamma(l,m,n) !Symmetry
                  enddo
                enddo
              enddo

              !Store grid quantities
              gmetric%grid(igrid)%jac(i,j,k) = jac
              gmetric%grid(igrid)%gsub (i,j,k,:,:)   = gsub
              gmetric%grid(igrid)%gsup (i,j,k,:,:)   = gsup
              gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = gamma
              gmetric%grid(igrid)%cov  (i,j,k,:,:)   = cov
              gmetric%grid(igrid)%cnv  (i,j,k,:,:)   = cnv

            enddo
          enddo
        enddo

        !Enforce topological constraints on Christoffel symbols
        do i=1,3
          do j=1,3
            do k=1,3
              call topol_bc(gmetric%grid(igrid)%Gamma(:,:,:,i,j,k))
            enddo
          enddo
        enddo

        !Deallocate auxiliary variables
        deallocate(dr)

      enddo

c     End program

      end subroutine restrictGridMetrics

c     topol_bc
c     #################################################################
      subroutine topol_bc(array)

c     -----------------------------------------------------------------
c     !Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        real(8) :: array(0:nxp,0:nyp,0:nzp)

c     Begin program

        if (bcond(1) == PER .or. bcond(2) == PER) then
          array(0  ,:,:) = array(nxp-1,:,:)
          array(nxp,:,:) = array(1    ,:,:)
        endif

        if (bcond(3) == PER .or. bcond(4) == PER) then
          array(:,0  ,:) = array(:,nyp-1,:)
          array(:,nyp,:) = array(:,1    ,:)
        endif
        
        if (bcond(5) == PER .or. bcond(6) == PER) then
          array(:,:,0  ) = array(:,:,nzp-1)
          array(:,:,nzp) = array(:,:,1    )
        endif

c     End program

      end subroutine topol_bc

c     gammaZeroForce
c     #################################################################
      subroutine gammaZeroForce(i,j,k,igx,igy,igz)

c     -----------------------------------------------------------------
c     Postprocess Christoffel symbol to satisfy zero-force condition
c     to machine accuracy
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer(4) :: i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ii,ll,kk,mm,jj,ig,jg,kg
        real(8)    :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const
        real(8)    :: hess(3,3,3),table(3,3,3),gsub(3,3)
     .               ,gsup(3,3),gsupp(3,3),gsupm(3,3),gamma(3,3,3)

c     Begin program

        if  (    (i == 0 .or. i == grid_params%nxv(igrid)+1) !Avoid x-boundaries
     .       .or.(j == 0 .or. j == grid_params%nyv(igrid)+1) !Avoid y-boundaries
     .       .or.(k == 0 .or. k == grid_params%nzv(igrid)+1) !Avoid z-boundaries
     .       .or.(i == 1 .and. bcond(1) == SP)               !Avoid singular point
     .       ) return

        gamma = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
        gsub  = gmetric%grid(igrid)%gsub(i,j,k,:,:)
        gsup  = gmetric%grid(igrid)%gsup(i,j,k,:,:)

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        do ii=1,3 !cycle through Chistoffel symbols

          if (ii == 2 .and. bcond(1) == SP) then
            jac  = gmetric%grid(igrid)%jac(i,j,k)
            const= 2d0
          else
            jac=1d0
            const=1d0
          endif

          do ll=1,3
            do kk=1,3
              select case (kk)
              case(1)
                dh = 2*grid_params%dxh(ig)
                if (ii == 2 .and. bcond(1) == SP) then
                  jacp = gmetric%grid(igrid)%jac(i+1,j,k)
                  jacm = gmetric%grid(igrid)%jac(i-1,j,k)
                else
                  jacp=1d0
                  jacm=1d0
                endif
                gsupp = jacp*gmetric%grid(igrid)%gsup (i+1,j,k,:,:)
                gsupm = jacm*gmetric%grid(igrid)%gsup (i-1,j,k,:,:)
              case(2)
                dh = 2*grid_params%dyh(jg)
                if (ii == 2 .and. bcond(1) == SP) then
                  jacp = gmetric%grid(igrid)%jac(i,j+1,k)
                  jacm = gmetric%grid(igrid)%jac(i,j-1,k)
                else
                  jacp=1d0
                  jacm=1d0
                endif
                gsupp = jacp*gmetric%grid(igrid)%gsup (i,j+1,k,:,:)
                gsupm = jacm*gmetric%grid(igrid)%gsup (i,j-1,k,:,:)
              case(3)
                dh = 2*grid_params%dzh(kg)
                if (ii == 2 .and. bcond(1) == SP) then
                  jacp = gmetric%grid(igrid)%jac(i,j,k+1)
                  jacm = gmetric%grid(igrid)%jac(i,j,k-1)
                else
                  jacp=1d0
                  jacm=1d0
                endif
                gsupp = jacp*gmetric%grid(igrid)%gsup (i,j,k+1,:,:)
                gsupm = jacm*gmetric%grid(igrid)%gsup (i,j,k-1,:,:)
              end select
                
              table(ii,ll,kk) =
     .              -gsub(ll,1)*(gsupp(1,ii)-gsupm(1,ii))/dh/jac
     .              -gsub(ll,2)*(gsupp(2,ii)-gsupm(2,ii))/dh/jac
     .              -gsub(ll,3)*(gsupp(3,ii)-gsupm(3,ii))/dh/jac
     .              +const*delta(ii,ll)*(gamma(1,kk,1)
     .                                  +gamma(2,kk,2)
     .                                  +gamma(3,kk,3))

              summ=0d0
              do mm=1,3
                do jj=1,3
                  summ = summ + gsub(ll,mm)*gsup(ii,jj)*gamma(mm,jj,kk)
                enddo
              enddo

              table(ii,ll,kk) = table(ii,ll,kk) - summ
            enddo
          enddo

        enddo

        gmetric%grid(igrid)%Gamma(i,j,k,:,:,:) = table

c      START CHECKS

cc        write (*,*)
cc        write (*,*) 'Grid:',igrid,'  Grid node:',i,j,k

c      Check difference between correction and prediction

cc        table = table - gamma
cc
cc        summ2= sqrt(sum(gamma*gamma))
cc        summ = sqrt(sum(table*table))/summ2
cc
cc        write (*,*) summ

c      Check zero-force property: g^lk G^i_lk = -nabla_m(g^mi)

cc        do ii=1,3 !cycle through Chistoffel symbols
cc
cc          summ = 0d0
cc
cc          dh = 2*grid_params%dxh(ig)
cc          gsupp = gmetric%grid(igrid)%gsup (i+1,j,k,:,:)
cc          gsupm = gmetric%grid(igrid)%gsup (i-1,j,k,:,:)
cc
cc          summ = summ -(gsupp(1,ii)-gsupm(1,ii))/dh
cc
cc          dh = 2*grid_params%dyh(jg)
cc          gsupp = gmetric%grid(igrid)%gsup (i,j+1,k,:,:)
cc          gsupm = gmetric%grid(igrid)%gsup (i,j-1,k,:,:)
cc
cc          summ = summ -(gsupp(2,ii)-gsupm(2,ii))/dh
cc
cc          dh = 2*grid_params%dzh(kg)
cc          gsupp = gmetric%grid(igrid)%gsup (i,j,k+1,:,:)
cc          gsupm = gmetric%grid(igrid)%gsup (i,j,k-1,:,:)
cc
cc          summ = summ -(gsupp(3,ii)-gsupm(3,ii))/dh
cc
cc          summ2 = 0d0
cc          do ll=1,3
cc            do kk=1,3
cc              summ2 = summ2 + gsup(ll,kk)*gamma(ii,ll,kk)
cc            enddo
cc          enddo
cc
cc          write (*,*) summ-summ2
cc        enddo

c      Check cancellation property g^lk(d_il G^j_kj-g_lm g^ij G^m_jk) = 0

cc        do ii=1,3 !cycle through Chistoffel symbols
cc          do ll=1,3
cc            do kk=1,3
cc              table(ii,ll,kk) =
cc     .               delta(ii,ll)*(gamma(1,kk,1)
cc     .                            +gamma(2,kk,2)
cc     .                            +gamma(3,kk,3))
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
cc          summ = 0d0
cc          do ll=1,3
cc            do kk=1,3
cc              summ = summ + gsup(ll,kk)*table(ii,ll,kk)
cc            enddo
cc          enddo
cc
cc          write (*,*) summ
cc        enddo

c     End program

      end subroutine gammaZeroForce

      end subroutine defineGridMetric

c     cross_product
c     #################################################################
      function cross_product(vec1,vec2) result(vec3)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: vec1(3),vec2(3),vec3(3)

c     Local variables

c     Begin program

        vec3(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
        vec3(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
        vec3(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)

      end function cross_product

c     triple_product
c     #################################################################
      function triple_product(vec1,vec2,vec3) result(scalar)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: vec1(3),vec2(3),vec3(3),scalar

c     Local variables

c     Begin program

        scalar = dot_product(vec1,cross_product(vec2,vec3))

      end function triple_product

c     delta
c     #################################################################
      real(8) function delta(i,j)

        integer(4) :: i,j

        delta = 0d0
        if (i == j) delta = 1d0

      end function delta

      end module grid_metric

c module grid_operations
c #####################################################################
      module grid_operations

c ---------------------------------------------------------------------
c     This module packs routines that perform operations on grid
c     quantities, such as coordinate transformation of vector components,
c     vector norms and scalar products. It contains the following
c     routines:
c          * transformVectorToCartesian 
c          * transformVectorToCurvilinear 
c          * transformFromCurvToCurv 
c          * volume
c          * vectorNorm
c          * scalarProduct
c     It is assumed that the grid metric structure gmetric is allocated 
c     and filled.
c ---------------------------------------------------------------------

        use grid_metric

        implicit none

      contains

c     transformVectorToCartesian
c     #################################################################
      subroutine transformVectorToCartesian(i,j,k,igx,igy,igz
     .                                     ,c1,c2,c3,covariant
     .                                     ,cx,cy,cz)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: covariant
        

c     Local variables

        integer(4) :: ic
        real(8)    :: T_to_car(3,3),vec(3)

c     Begin program

        if (covariant) then

          do ic =1,3
            T_to_car(:,ic) = gmetric%grid(igx)%cov(i,j,k,ic,:)
          enddo

          vec = (/ c1,c2,c3 /)

          vec = matmul(T_to_car,vec)

          cx = vec(1)
          cy = vec(2)
          cz = vec(3)
            
        else

          do ic =1,3
            T_to_car(:,ic) = gmetric%grid(igx)%cnv(i,j,k,ic,:)
          enddo

          vec = (/ c1,c2,c3 /)

          vec = matmul(T_to_car,vec)

          cx = vec(1)
          cy = vec(2)
          cz = vec(3)

        endif

      end subroutine transformVectorToCartesian

c     transformVectorToCurvilinear
c     #################################################################
      subroutine transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .                                       ,cx,cy,cz,covariant
     .                                       ,c1,c2,c3)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector (cx,cy,cz) to curvilinear (c1,c2,c3)
c     (covariant if covariant=.true., contravariant otherwise)
c     at grid coordinates (i,j,k). 
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: covariant

c     Local variables

        integer(4) :: ic
        real(8)    :: T_to_curv(3,3),vec(3),jac

c     Begin program

        jac = gmetric%grid(igx)%jac(i,j,k)

        if (covariant) then

          do ic =1,3
            T_to_curv(:,ic) = gmetric%grid(igx)%cnv(i,j,k,ic,:)
          enddo

          T_to_curv = jac*transpose(T_to_curv)

          vec = (/ cx,cy,cz /)

          vec = matmul( T_to_curv,vec)

          c1 = vec(1)
          c2 = vec(2)
          c3 = vec(3)
            
        else

          do ic =1,3
            T_to_curv(:,ic) = gmetric%grid(igx)%cov(i,j,k,ic,:)
          enddo

          T_to_curv = jac*transpose(T_to_curv)

          vec = (/ cx,cy,cz /)

          vec = matmul(T_to_curv,vec)

          c1 = vec(1)
          c2 = vec(2)
          c3 = vec(3)

        endif

      end subroutine transformVectorToCurvilinear

c     transformFromCurvToCurv
c     #################################################################
      subroutine transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at grid coordinates (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical    :: tocnv

c     Local variables

        real(8)    :: cov(3),cnv(3)

c     Begin program

        if (tocnv) then
          cnv  = matmul(gmetric%grid(igx)%gsup(i,j,k,:,:)
     .                 , (/ cov1,cov2,cov3 /))
          cnv1 = cnv(1)
          cnv2 = cnv(2)
          cnv3 = cnv(3)
        else
          cov  = matmul(gmetric%grid(igx)%gsub(i,j,k,:,:)
     .                 , (/ cnv1,cnv2,cnv3 /))
          cov1 = cov(1)
          cov2 = cov(2)
          cov3 = cov(3)
        endif

      end subroutine transformFromCurvToCurv

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

c     vectorNorm
c     ################################################################
      real(8) function vectorNorm(i,j,k,igx,igy,igz,ax,ay,az,covar)

c     ---------------------------------------------------------------
c     Finds norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: i,j,k,igx,igy,igz
      real(8)    :: ax,ay,az
      logical    :: covar

c     Local variables

      real(8)    :: tensor(3,3),cnv(3),cov(3),jac

c     Begin program

      jac = gmetric%grid(igx)%jac(i,j,k)

      if (covar) then
        tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
        cov = (/ ax,ay,az /)
        cnv = matmul(tensor,cov)
      else
        tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
        cnv = (/ ax,ay,az /)
        cov = matmul(tensor,cnv)
      endif

      vectorNorm = dot_product(cov,cnv)/jac

c     End 

      end function vectorNorm

c     scalarProduct
c     ################################################################
      function scalarProduct(i,j,k,igx,igy,igz,cov1,cov2,cov3
     .                       ,cnv1,cnv2,cnv3) result (dot)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: i,j,k,igx,igy,igz
      real(8)    :: dot,cov1,cov2,cov3,cnv1,cnv2,cnv3

c     Local variables

      real(8)    :: cnv(3),cov(3)

c     Begin program

      cnv = (/ cnv1,cnv2,cnv3 /)
      cov = (/ cov1,cov2,cov3 /)

      dot = dot_product(cov,cnv)/gmetric%grid(igx)%jac(i,j,k)

c     End 

      end function scalarProduct

      end module grid_operations

c module grid
c #####################################################################
      module grid

        use grid_operations

        implicit none

        integer(4),private :: nxx,nyy,nzz

        integer(4) :: mg_ratio

        real(8),private :: pi

      contains

c     createGrid
c     #################################################################
      subroutine createGrid(nx,ny,nz)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nx,ny,nz

c     Local variables

        integer(4) :: i,n1,n2,n3
        integer(4) :: ngrdx,ngrdy,ngrdz

c     Begin program

        nxx = nx
        nyy = ny
        nzz = nz

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

c     Initialize MG arrays

        grid_params%mg_ratio_x = 1
        grid_params%mg_ratio_y = 1
        grid_params%mg_ratio_z = 1

        grid_params%nxv(1) = nxx
        do i = 2,ngrdx
          grid_params%nxv(i) = grid_params%nxv(i-1) / mg_ratio
          grid_params%mg_ratio_x(i-1) = mg_ratio
        enddo
        do i = ngrdx+1,grid_params%ngrid
          grid_params%nxv(i) = grid_params%nxv(i-1)
        enddo

        grid_params%nyv(1) = nyy
        do i = 2,ngrdy
          grid_params%nyv(i) = grid_params%nyv(i-1) / mg_ratio
          grid_params%mg_ratio_y(i-1) = mg_ratio
        enddo
        do i = ngrdy+1,grid_params%ngrid
          grid_params%nyv(i) = grid_params%nyv(i-1)
        enddo

        grid_params%nzv(1) = nzz
        do i = 2,ngrdz
          grid_params%nzv(i) = grid_params%nzv(i-1) / mg_ratio
          grid_params%mg_ratio_z(i-1) = mg_ratio
        enddo
        do i = ngrdz+1,grid_params%ngrid
          grid_params%nzv(i) = grid_params%nzv(i-1)
        enddo

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

c     Define uniform logical grid on ALL grid levels

        call createLogicalGrid(nxx,grid_params%xx,grid_params%dx
     .                        ,grid_params%dxh,grid_params%nxv
cc     .                        ,grid_params%ngrdx,grid_params%istartx
     .                        ,grid_params%ngrid,grid_params%istartx
     .                        ,xmin,xmax,bcond(1),bcond(2))

        call createLogicalGrid(nyy,grid_params%yy,grid_params%dy
     .                        ,grid_params%dyh,grid_params%nyv
cc     .                        ,grid_params%ngrdy,grid_params%istarty
     .                        ,grid_params%ngrid,grid_params%istarty
     .                        ,ymin,ymax,bcond(3),bcond(4))

        call createLogicalGrid(nzz,grid_params%zz,grid_params%dz
     .                        ,grid_params%dzh,grid_params%nzv
cc     .                        ,grid_params%ngrdz,grid_params%istartz
     .                        ,grid_params%ngrid,grid_params%istartz
     .                        ,zmin,zmax,bcond(5),bcond(6))

c     Store grid metric parameters in grid metric structure

        if (checkGridDatabase()) call defineGridMetric(gmetric)

      end subroutine createGrid

c     createLogicalGrid
c     #################################################################
      subroutine createLogicalGrid (nn,xx,dx,dxh,nx,ngrid,istart
     .                             ,lmin,lmax,bcs1,bcs2)

        implicit none

c     Call variables

        integer(4) :: nn,ngrid,nx(ngrid),istart(ngrid),bcs1,bcs2
        real(8)    :: xx(*),dx(*),dxh(*),lmin,lmax

c     Local variables
        
        integer(4) :: ig,i,isig
        real(8)    :: dh,length

c     Begin program

        length = lmax-lmin

c     Periodic

        if (bcs1 == PER .or. bcs2 == PER) then

          do ig = 1,ngrid

            isig = istart(ig)

          !Find cell centers
            dh = length/dfloat(nx(ig))

            xx(1 + isig) = lmin
            do i = 2,nx(ig)+1
              xx(i + isig) = xx(i-1 + isig) + dh
            enddo
            xx(0 + isig) = xx(1+isig) - dh
cc            xx(0 + isig) = xx(nx(ig)+isig)
cc            xx(nx(ig)+1 +isig) = xx(1+isig)

          !Find integer mesh spacings
            do i = 1,nx(ig)
              dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
            enddo
            dx(0       +isig) = dx(nx(ig)+isig)
            dx(nx(ig)+1+isig) = dx(1     +isig)

          !Find half mesh spacings
            do i = 1,nx(ig)+1
              dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
            enddo
            dxh(0 + isig) = dxh(nx(ig) + isig)

          enddo

c     Radial (singular point at r=0)

        elseif (bcs1 == SP) then

          do ig = 1,ngrid

            isig = istart(ig)

          !Find cell centers
            dh = length/dfloat(nx(ig))

            xx(0 + isig) = 1d-8*dh
            xx(1 + isig) = dh/2.
            do i = 2,nx(ig)+1
cc            xx(0 + isig) = -dh/2.
cc            do i = 1,nx(ig)+1
              xx(i + isig) = xx(i-1 + isig) + dh
            enddo

          !Find integer mesh spacings
cc            dx(0 + isig) = dh/2.
cc            do i = 1,nx(ig)
            do i = 0,nx(ig)
              dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
            enddo

          !Find half mesh spacings
            dxh(1 + isig) = dh
            do i = 2,nx(ig)
cc            do i = 1,nx(ig)
              dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
            enddo
            dxh(0        + isig) = dx(0      + isig)/2.
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)/2.

cc            xx(0 + isig) = 1d-8
          enddo

c     Other

        else

          do ig = 1,ngrid

            isig = istart(ig)

          !Find cell centers
c-ncv            dh = length/dfloat(nx(ig)+1)
            dh = length/dfloat(nx(ig))

c-ncv            xx(0 + isig) = lmin
            xx(0 + isig) = lmin-dh/2.
            do i = 1,nx(ig)+1
              xx(i + isig) = xx(i-1 + isig) + dh
            enddo

          !Find integer mesh spacings
            do i = 0,nx(ig)
              dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
            enddo

          !Find half mesh spacings
            do i = 1,nx(ig)
              dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
            enddo
            dxh(0        + isig) = dx(0      + isig)/2.
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)/2.
          enddo

        endif

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
cc        if (nxx == 1) then
cc          xmin = 0d0
cc          xmax = 1d-3
cc        endif
cc        if (nyy == 1) then
cc          ymin = 0d0
cc          ymax = 1d-3
cc        endif
cc        if (nzz == 1) then
cc          zmin = 0d0
cc          zmax = 1d-3
cc        endif

      end subroutine consistencyCheck

c     checkGrid
c     #################################################################
      subroutine checkGrid

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: i,j,n1,n2,n3
        integer(4) :: igx,igy,igz,isig

c     Begin program

c     Multigrid parameters

        igx = grid_params%ngrid
        igy = grid_params%ngrid
        igz = grid_params%ngrid

cc        igx = 1
cc        igy = 1 
cc        igz = 1

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

        call metricTensorCheck(igx,igy,igz
     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)

cc        call hessianCheck(igx,igy,igz
cc     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)

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

c     hessianCheck
c     #################################################################
      subroutine hessianCheck(igx,igy,igz,nxv,nyv,nzv)

        implicit  none

c     Call variables

        integer(4) :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

c     Local variables

        integer(4) :: i,j,k,i1,j1,k1,ig,jg,kg,ih
        real(8)    :: x1,y1,z1
        real(8)    :: hess(3,3,3),hess_cnv(3,3,3),table(3,3,3)

c     Begin program

        write (*,*) 
        write (*,*) '*************************'
        write (*,*) '     Hessian check'
        write (*,*) '*************************'
        write (*,*)

        do i = 2,nxv(igx)   !Start at 2 to avoid singular points
          do j = 1,nyv(igy)
            do k = 1,nzv(igz)

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = grid_params%xx(ig)
              y1 = grid_params%yy(jg)
              z1 = grid_params%zz(kg)

              hess = -christ_2knd(i,j,k,igx,igy,igz)

              hess_cnv(1,:,:) = hessian_cnv(1,x1,y1,z1)
              hess_cnv(2,:,:) = hessian_cnv(2,x1,y1,z1)
              hess_cnv(3,:,:) = hessian_cnv(3,x1,y1,z1)

              do i1=1,3
                do j1=1,3
                  do k1=1,3
                    table(i1,j1,k1) = hess_cnv(i1,j1,k1)
     .                               -delta(i1,k1)*(hess(1,j1,1)
     .                                             +hess(2,j1,2)
     .                                             +hess(3,j1,3))
     .                               +hess(k1,i1,j1)
                  enddo
                enddo
              enddo

              write (*,5) 'Grid point: (',x1,',',y1,',',z1,')'
 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)

              write (*,*) 
              write (*,*) 'Hessian relation'

              do ih=1,3
                write (*,*) 
                write (*,10) table(ih,1,1:3)
                write (*,10) table(ih,2,1:3)
                write (*,10) table(ih,3,1:3)
 10             format (3f10.3)
              enddo

            enddo
          enddo
        enddo

c     End program

      end subroutine hessianCheck

c     metricTensorCheck
c     #################################################################
      subroutine metricTensorCheck(igx,igy,igz,nxv,nyv,nzv)

        implicit  none

c     Call variables

        integer(4) :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

c     Local variables

        integer(4) :: i,j,k,i1,j1,k1,ig,jg,kg,ih
        real(8)    :: x1,y1,z1,check
        real(8)    :: gup(3,3),gdown(3,3),tensor(3,3)
        logical    :: cartesian

c     Begin program

        check = 0d0

        write (*,*) 
        write (*,*) '*************************'
        write (*,*) '   Metric Tensor check'
        write (*,*) '*************************'
        write (*,*)

        do i = 2,nxv(igx)   !Start at 2 to avoid singular points
          do j = 1,nyv(igy)
            do k = 1,nzv(igz)

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = grid_params%xx(ig)
              y1 = grid_params%yy(jg)
              z1 = grid_params%zz(kg)

cc              gup   = g_sup(i,j,k,igx,igy,igz)
cc              gdown = g_sub(i,j,k,igx,igy,igz)
              gup   = gmetric%grid(igx)%gsup(i,j,k,:,:)
              gdown = gmetric%grid(igx)%gsub(i,j,k,:,:)

              tensor = matmul(gup,gdown)

              write (*,5) 'Grid point: (',x1,',',y1,',',z1,')'
 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)

              write (*,*) 
              write (*,*) 'Metric tensor product'
              write (*,10) tensor(1,1:3)
              write (*,10) tensor(2,1:3)
              write (*,10) tensor(3,1:3)
 10           format (3f10.3)

            enddo
          enddo
        enddo

c     End program

      end subroutine metricTensorCheck

      end subroutine checkGrid

      end module grid


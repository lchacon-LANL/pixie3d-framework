c module debug
c ######################################################################
      module debug

        real(8) :: lxmax,lxmin,lymax,lymin
        real(8),allocatable,dimension(:,:,:,:) :: dbg

      contains

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c     ---------------------------------------------------------------------
c     Contours 2D array in xdraw format. In call:
c       * arr: 2D array to be plotted
c       * nx,ny: dimensions of array
c       * xmin,xmax,ymin,ymax: 2D domain limits
c       * iopt: whether to initialize xdraw plot (iopt=0) or not.
c       * nunit: integer file identifier.
c     ---------------------------------------------------------------------

c     Call variables

      integer(4) :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer(4) :: i,j

c     Begin program

      if(iopt == 0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4)
     .              ,real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*)

c     End program

      end subroutine contour

      end module debug

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

        if (grid_st2%ngrdx /= ngridx) nullify(grid_st2%xx)

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

c$$$        ig = i -grid_params%ilo(igx)+1 + grid_params%istartx(igx)
c$$$        jg = j -grid_params%jlo(igy)+1 + grid_params%istarty(igy)
c$$$        kg = k -grid_params%klo(igz)+1 + grid_params%istartz(igz)
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
        case('car','scl','cyl','hel','hl2','tor','sin')

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

        external   map

c     Begin program

        if (checkGridDatabase()) then

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

          car = x_xi(x1,y1,z1)

          x1 = car(1)
          y1 = car(2)
          z1 = car(3)

        else

          call map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1) !External

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
        case ('hl2')
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (x2-aa*x3)

          xx = sqrt(2*x1)*cos(phi)
          yy = sqrt(2*x1)*sin(phi)
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
        case ('hl2')
          jac = 1d0
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

c     covariantVector
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

        real(8)    :: car(3),curv(3),jac,r

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
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ sqrt(2*x1),x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))

          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)*curv(1)
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
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
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

c     contravariantVector
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
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ sqrt(2*x1),x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(phi),cos(phi),0d0 /)*curv(1)
            case (3)
              vec = (/ aa*sin(phi)*curv(1),-aa*cos(phi)*curv(1),1d0 /)
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

        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ 2*x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (l)
            case (1)
              vec(1) = -1d0/curv(1)
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = -curv(1)
              vec(6) = aa*curv(1)
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = -aa**2*curv(1)
            case (2)
              vec(1) = 0d0
              vec(2) = 1./curv(1)
              vec(3) = -aa/curv(1)
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
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ 2*x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1 + aa**2*curv(1) /)
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
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ 2*x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2, aa
     .            ,0d0, aa, 1d0 /)
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

c module bc_def
c #####################################################################
      module bc_def

c     -----------------------------------------------------------------
c     This module defines the types of BCs possible:
c     Topological BCs (set topology of computational domain):
c       PER: Periodic
c       SP : Singular point (as in cylindrical, toroidal)
c       SP2: Singular point with second-order variable grid
c       SYM: Symmetry (homogeneous Neumann)
c     Standard BCs:
c       EQU: Inhomogeneous dirichlet (generally imposed by equilibrium)
c       DIR: Homogeneous dirichlet
c       NEU: Neumann (homogeneous/inhomogeneous)
c       SYM: Symmetry (odd,even depend on sign)
c       DEF: Default (set in application code)
c       EXT: By extrapolation
c     The actual integers chosen indicate a particular order that the BC
c     scheduler follows to impose the boundary conditions.
c     -----------------------------------------------------------------

        use grid_structures

        implicit none

        integer(4) :: PER,DIR,NEU,SP,SP2,EQU,DEF,EXT,SYM
        parameter (DEF=0,PER=1,EQU=2,EXT=3,NEU=4,SYM=5,DIR=6,SP=7
     .            ,SP2=8)

        integer(4),parameter :: BCLIM=6

        integer(4) :: bcond(6)
        
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

        sing_bc=.false.
        if (bcond(iibc) == SP .or. bcond(iibc) == SP2) sing_bc=.true.

      end function bcSP

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        integer(4) :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point=.false.
        if (i+grid_params%ilo(igx)-1==1.and.bcSP()) sing_point=.true.

      end function isSP

c     isSP1
c     #################################################################
      function isSP1(i,j,k,igx,igy,igz) result(sing_point)

        integer(4) :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point=.false.
        if (i+grid_params%ilo(igx)-1==1.and.bcond(1)==SP)
     .       sing_point=.true.

      end function isSP1

c     isSP2
c     #################################################################
      function isSP2(i,j,k,igx,igy,igz) result(sing_point)

        integer(4) :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point=.false.
        if (i+grid_params%ilo(igx)-1==1.and.bcond(1)==SP2)
     .       sing_point=.true.

      end function isSP2

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

#if defined(petsc)
c module local_mpi
c ######################################################################
      module local_mpi

        use grid_definition

        use bc_def

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

        integer(4),private :: il ,ih ,jl ,jh ,kl ,kh
     .                       ,ilm,ihp,jlm,jhp,klm,khp

        integer(4) :: l_lim(3,0:1),g_lim(3,0:1),rem_l_lim(3,0:1)
     .               ,rp_l_lim(3,0:1)

        integer(4) :: MPI_COMM_SP
        integer(4),allocatable,dimension(:) :: MPI_COMM_PER
        integer(4),allocatable,dimension(:) :: MPI_COMM_NBRS

        type :: userctx
cc#define PETSC_AVOID_DECLARATIONS
cc#include "finclude/petsc.h"
cc#include "finclude/petscvec.h"
cc#include "finclude/petscda.h"
cc#undef PETSC_AVOID_DECLARATIONS
          DA      :: da
          Vec     :: Xg
          integer :: xs,xe,xm,gxs,gxe,gxm
          integer :: ys,ye,ym,gys,gye,gym
          integer :: zs,ze,zm,gzs,gze,gzm
          integer :: lgxs,lgxe
          integer :: lgys,lgye
          integer :: lgzs,lgze
          integer :: mx,my,mz,rank
          integer :: igx,igy,igz
        end type userctx

        type (userctx) :: dactx

      contains

c     initMPI
c     #####################################################################
      subroutine initMPI

c     ---------------------------------------------------------------------
c     Initialize needed MPI constructs
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: dim

c     Begin program

c     Rank and number of processors

        call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
        call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)

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

c     Neighbors MPI communicators (PETSc)

        call createDA

      end subroutine initMPI

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

c     createDA
c     #################################################################
      subroutine createDA

c     -----------------------------------------------------------------
c     Creates PETSc distributed array for boundary communication in
c     fortran.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: BC,ierr,nxg,nyg,nzg,igr

c     Begin program

        igr = 1

        nxg = grid_params%nxgl(igr)
        nyg = grid_params%nygl(igr)
        nzg = grid_params%nzgl(igr)

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
     &                 ,PETSC_NULL_INTEGER,dactx%da,ierr)


        call DACreateGlobalVector(dactx%da,dactx%Xg,ierr)

c       Get local grid boundaries in global grid

        call DAGetCorners(dactx%da,dactx%xs,dactx%ys,dactx%zs
     &     ,dactx%xm,dactx%ym,dactx%zm,ierr)
        call DAGetGhostCorners(dactx%da,dactx%gxs,dactx%gys
     &     ,dactx%gzs,dactx%gxm,dactx%gym,dactx%gzm,ierr)

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

       dactx%igx = igr
       dactx%igy = igr
       dactx%igz = igr

c     Define local limits

        call fromGlobalToLocalLimits(dactx%gxs ,dactx%gys ,dactx%gzs
     $                              ,dactx%lgxs,dactx%lgys,dactx%lgzs
     $                              ,dactx%igx,dactx%igy,dactx%igz)
        call fromGlobalToLocalLimits(dactx%gxe ,dactx%gye ,dactx%gze
     $                              ,dactx%lgxe,dactx%lgye,dactx%lgze
     $                              ,dactx%igx,dactx%igy,dactx%igz)

       il = 1
       jl = 1
       kl = 1
       ih = grid_params%nxv(igr)
       jh = grid_params%nyv(igr)
       kh = grid_params%nzv(igr)

       ilm = il-1
       ihp = ih+1
       jlm = jl-1
       jhp = jh+1
       klm = kl-1
       khp = kh+1

c     End program

      end subroutine createDA

c     destroyDA
c     #################################################################
      subroutine destroyDA

        implicit none

c     Call variables

c     Local variables

        integer(4) :: ierr

c     Begin program

        call VecDestroy(dactx%Xg,ierr)
        call DADestroy (dactx%da,ierr)

c     End program

      end subroutine destroyDA

c     fillLocalVec
c     #################################################################
      subroutine fillLocalVec(array,x)

        implicit none

c     Call variables

        real(8) :: array(ilm:ihp,jlm:jhp,klm:khp)
     .            ,x    (dactx%gxs:dactx%gxe 
     .                  ,dactx%gys:dactx%gye 
     .                  ,dactx%gzs:dactx%gze)

c     Local variables

c     Begin program

        x(dactx%xs:dactx%xe 
     .   ,dactx%ys:dactx%ye 
     .   ,dactx%zs:dactx%ze) = array(il:ih,jl:jh,kl:kh)

c     End program

      end subroutine fillLocalVec

c     emptyLocalVec
c     #################################################################
      subroutine emptyLocalVec(array,x)

        implicit none

c     Call variables

        real(8) :: array(ilm:ihp,jlm:jhp,klm:khp)
     .            ,x    (dactx%gxs:dactx%gxe 
     .                  ,dactx%gys:dactx%gye 
     .                  ,dactx%gzs:dactx%gze)

c     Local variables

c     Begin program

        array(dactx%lgxs:dactx%lgxe
     .       ,dactx%lgys:dactx%lgye
     .       ,dactx%lgzs:dactx%lgze) = x(dactx%gxs:dactx%gxe 
     .                                  ,dactx%gys:dactx%gye 
     .                                  ,dactx%gzs:dactx%gze)

c     End program                     

      end subroutine emptyLocalVec

c     fillPetscGhostCells
c     #################################################################
      subroutine fillPetscGhostCells(array)

        implicit none

c     Call variables

        real(8) :: array(ilm:ihp,jlm:jhp,klm:khp)

c     Local variables

        Vec :: localX

        PetscScalar,pointer :: lx_v(:)

        integer(4) :: ierr

c     Begin program

c     Fill known values of dactx%Xg

        !Get pointer to vector data
        call DAGetLocalVector(dactx%da,localX,ierr)
        call VecGetArrayF90  (localX,lx_v,ierr)

        !Assignment: array -> lx_v
        call fillLocalVec(array,lx_v)

        !Restore vector
        call VecRestoreArrayF90(localX,lx_v,ierr)

        !Insert values into global vector
        call DALocalToGlobal(dactx%da,localX,INSERT_VALUES,dactx%Xg
     $                      ,ierr)
cc        call DARestoreLocalVector(dactx%da,localX,ierr)
  
c     Fill ghost cells

        !Get pointer to vector data w/ghost cells
cc        call DAGetLocalVector    (dactx%da,localX,ierr)
        call DAGlobalToLocalBegin(dactx%da,dactx%Xg,INSERT_VALUES,localX
     $                           ,ierr)
        call DAGlobalToLocalEnd  (dactx%da,dactx%Xg,INSERT_VALUES,localX
     $                           ,ierr)

        call VecGetArrayF90(localX,lx_v,ierr)

        !Assignment: lx_v -> array
        call emptyLocalVec(array,lx_v)

        !Restore vector
        call VecRestoreArrayF90(localX,lx_v,ierr)

c     Deallocate memory

        call DARestoreLocalVector(dactx%da,localX,ierr)

c     End program

      end subroutine fillPetscGhostCells

      end module local_mpi

#endif

c module error
c ######################################################################
      module error

#if defined(petsc)
        use local_mpi
#endif

        character(80) :: messg

      contains

c     pstop
c     ################################################################
      subroutine pstop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        character(*)  :: routine, message

c     Begin program

#if defined(petsc)
        if (my_rank == 0) then
          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)
        endif

        call PetscEnd(mpierr)
#else
        write (*,*)
        write (*,*) trim(message)
        write (*,*) 'Program stopped at routine ',trim(routine)
#endif

        stop

      end subroutine pstop

      end module error

c module grid_metric
c #####################################################################
      module grid_metric

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
c        * cross_product
c        * triple_product
c        * delta
c ---------------------------------------------------------------------

        use grid_definition

#if defined(petsc)
        use local_mpi
#else
        use bc_def
#endif

        implicit none

        type :: grid_metrics
cc          real(8),pointer,dimension(:,:,:,:)     :: car   !Node positions in cartesian space
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
cc            allocate(gmetric%grid(igrid)
cc     $           %car  (ilom:ihip,jlom:jhip,klom:khip,3))
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
cc            deallocate(gmetric%grid(igrid)%car  )
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

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric(gmetric,igr)

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

        integer(4),optional,intent(IN) :: igr

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4) :: igrmin,igrmax,igrid,ilo,ihi,jlo,jhi,klo,khi
     $               ,ilom,ihip,jlom,jhip,klom,khip

        integer(4) :: i,j,k,igx,igy,igz
     .               ,i0,ip,im,j0,jp,jm,k0,kp,km,l,m,n,p
     .               ,ig,ig0,igm,igp,jg,jg0,jgm,jgp,kg,kg0,kgm,kgp
        real(8)    :: r(3,3),car0(3),carp(3),carm(3),dh(3),jac,ijac
     .               ,cnv(3,3),cov(3,3),gsub(3,3),gsup(3,3),vec(3)
     .               ,gamma(3,3,3),gamm1(3,3,3),mag,dhp,dhm,dhh
        real(8),allocatable,dimension(:,:,:,:,:) :: dr

c     Begin program

        anal_map = checkGridDatabase()

        call allocateGridMetric(gmetric)

        if (PRESENT(igr)) then
          igrmin = igr
          igrmax = igr
        else
          igrmin = 1
          igrmax = grid_params%ngrid
        endif

c       ANALYTICAL QUANTITIES **********************************

        if ((.not.numerical_grid).and.anal_map) then

          do igrid=igrmin,igrmax

            igx = igrid
            igy = igrid
            igz = igrid

            ilo = 1
            jlo = 1
            klo = 1
            ihi = grid_params%nxv(igrid)
            jhi = grid_params%nyv(igrid)
            khi = grid_params%nzv(igrid)

            ilom = ilo-1
            jlom = jlo-1
            klom = klo-1
            ihip = ihi+1
            jhip = jhi+1
            khip = khi+1

            do k = klom,khip
              do j = jlom,jhip
                do i = ilom,ihip
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
            if (igrid == 1) then
              do k = klo,khi
                do j = jlo,jhi
                  do i = ilo,ihi
                    call gammaZeroForce(i,j,k,igx,igy,igz)
                  enddo
                enddo
              enddo
            endif

          enddo

c       NUMERICAL QUANTITIES ***********************************

        else

          do igrid=igrmin,igrmax

            igx = igrid
            igy = igrid
            igz = igrid

            ilo = 1
            jlo = 1
            klo = 1
            ihi = grid_params%nxv(igrid)
            jhi = grid_params%nyv(igrid)
            khi = grid_params%nzv(igrid)

            ilom = ilo-1
            jlom = jlo-1
            klom = klo-1
            ihip = ihi+1
            jhip = jhi+1
            khip = khi+1

            allocate(dr(ilom:ihip,jlom:jhip,klom:khip,3,3))

            !Evaluate dx/dxi vectors
            do k = klom,khip
              do j = jlom,jhip
                do i = ilom,ihip
                  ip=min(i+1,ihip)
                  im=max(i-1,ilom)
                  jp=min(j+1,jhip)
                  jm=max(j-1,jlom)
                  kp=min(k+1,khip)
                  km=max(k-1,klom)

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

          !Evaluate grid quantities
            do k = klom,khip
              do j = jlom,jhip
                do i = ilom,ihip
                  r = dr(i,j,k,:,:)

                  !Evaluate Jacobian
                  jac = triple_product(r(1,:),r(2,:),r(3,:))
                  ijac = 1d0/jac

                  !Contravariant vectors
                  cnv(:,:) = r(:,:)*ijac

                  !Covariant vectors
                  cov(1,:) = cross_product(r(2,:),r(3,:))
                  cov(2,:) = cross_product(r(3,:),r(1,:))
                  cov(3,:) = cross_product(r(1,:),r(2,:))

                  cov = cov*ijac

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

                  if (i==ilom) then
                    i0=ilo
                  elseif (i==ihip) then
                    i0=ihi
                  else
                    i0=i
                  endif

                  if (j==jlom) then
                    j0=jlo
                  elseif (j==jhip) then
                    j0=jhi
                  else
                    j0=j
                  endif

                  if (k==klom) then
                    k0=klo
                  elseif (k==khip) then
                    k0=khi
                  else
                    k0=k
                  endif

                  ip=i0+1
                  im=i0-1
                  jp=j0+1
                  jm=j0-1
                  kp=k0+1
                  km=k0-1

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

                        gamma(l,m,n) = dot_product(vec,cov(l,:))
                        gamma(l,n,m) = gamma(l,m,n) !Symmetry
                      enddo
                    enddo
                  enddo

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
              do k = klo,khi
                do j = jlo,jhi
                  do i = ilo,ihi
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

cc      if (my_rank == 0) then
cc        open(unit=110,file='debug0.bin'
cc     .           ,form='unformatted',status='replace')
cc      else
cc        open(unit=110,file='debug1.bin'
cc     .           ,form='unformatted',status='replace')
cc      endif
cc      igx = 1
cc      do i=1,3
cc        do j=1,3
cc          nnx = grid_params%nxv(igx)
cc          nny = grid_params%nyv(igx)
cc          nnz = grid_params%nzv(igx)
cccc          debug = gmetric%grid(igx)%cnv(0:nx+1,0:ny+1,0:nz+1,i,j)
cccc          debug = gmetric%grid(igx)%Gamma(0:nx+1,0:ny+1,0:nz+1,1,i,j)
cccc     .           *gmetric%grid(igx)%jac
cc
cc          call contour(gmetric%grid(igrid)%Gamma(:,:,1,i,j,1)
cc     .                ,nnx+2,nny+2
cc     .                ,0d0,xmax,0d0,ymax,i+j-2,110)
cc        enddo
cc      enddo
cc      close(110)
cc      call MPI_Finalize(mpierr)
cc      stop

cc        if (my_rank == 0) write (*,*) 'grid_level',igrid

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

c     Begin program

        call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                              ,car(1),car(2),car(3))

      end function map

c     topol_bc
c     #################################################################
      subroutine topol_bc(array)

c     -----------------------------------------------------------------
c     !Enforce topological constraints on array
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        real(8)    :: array(ilom:ihip,jlom:jhip,klom:khip)

c     Local variables

        integer(4) :: nxl,nyl,nzl,nxg,nyg,nzg,dim

c     External

c     Begin program

#if !defined(petsc)

        nxl = grid_params%nxv(igrid)
        nyl = grid_params%nyv(igrid)
        nzl = grid_params%nzv(igrid)

        if (bcond(1) == PER .or. bcond(2) == PER) then
          array(0    ,:,:) = array(nxl,:,:)
          array(nxl+1,:,:) = array(1  ,:,:)
        endif

        if (bcond(3) == PER .or. bcond(4) == PER) then
          array(:,0    ,:) = array(:,nyl,:)
          array(:,nyl+1,:) = array(:,1  ,:)
        endif
        
        if (bcond(5) == PER .or. bcond(6) == PER) then
          array(:,:,0    ) = array(:,:,nzl)
          array(:,:,nzl+1) = array(:,:,1  )
        endif

#else
        if (igrid == 1) then
          call fillPetscGhostCells(array)
        endif
#endif

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

        integer(4) :: ii,ll,kk,mm,jj,ig,jg,kg,igrid
        real(8)    :: x1,y1,z1,summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),table(3,3,3),gsub(3,3)
     .               ,gsup(3,3),gsupp(3,3),gsupm(3,3),gamma(3,3,3)
        logical    :: alt_eom

c     Begin program

        if  (    (i == ilom .or. i == ihip)       !Avoid x-boundaries
     .       .or.(j == jlom .or. j == jhip)       !Avoid y-boundaries
     .       .or.(k == klom .or. k == khip)       !Avoid z-boundaries
     .       .or. isSP(i,j,k,igx,igy,igz)         !Avoid singular point
     .       ) return

        igrid = igx

        gamma = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
        gsub  = gmetric%grid(igrid)%gsub(i,j,k,:,:)
        gsup  = gmetric%grid(igrid)%gsup(i,j,k,:,:)

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        do ii=1,3 !cycle through Chistoffel symbols

          alt_eom = (ii == 2 .and. alt__eom())

          if (alt_eom) then
            jac  = gmetric%grid(igrid)%jac(i,j,k)
            const= 2d0
          else
            jac  = 1d0
            const= 1d0
          endif

          do ll=1,3
            do kk=1,3
              select case (kk)
              case(1)
                dh = 2*grid_params%dxh(ig)
                if (alt_eom) then
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
                if (alt_eom) then
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
                if (alt_eom) then
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

        !Symmetrize
        do ii=1,3
          do ll=1,3             !cycle through tensor components
            do kk=1,3
              gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
     $                    = 0.5*( table(ii,ll,kk)+table(ii,kk,ll) )
            enddo
          enddo
        enddo

cc        !Postprocessing for alt_eom case (ii == 2 .and. bcond(1) == SP) 
cc        if (bcSP()) then
cc
cc          gamma = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cc          jac   = gmetric%grid(igrid)%jac(i,j,k)
cc
cc          do ii=1,3             !cycle through Chistoffel symbols
cc            do ll=1,3
cc              do kk=1,3
cc                select case (kk)
cc                case(1)
cc                  dh = grid_params%dxh(ig)
cc                  jacp = hlf*(jac+gmetric%grid(igrid)%jac(i+1,j,k))
cc                  jacm = hlf*(jac+gmetric%grid(igrid)%jac(i-1,j,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i-1,j,k,:,:))
cc                case(2)
cc                  dh = grid_params%dyh(jg)
cc                  jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j+1,k))
cc                  jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j-1,k))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j-1,k,:,:))
cc                case(3)
cc                  dh = grid_params%dzh(kg)
cc                  jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k+1))
cc                  jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k-1))
cc                  gsupp = jacp
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
cc                  gsupm = jacm
cc     $                 *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k-1,:,:))
cc                end select
cc
cc                table(ii,kk,ll) =
cc     $                -gsub(ll,1)*(gsupp(1,ii)-gsupm(1,ii))/dh/jac
cc     .                -gsub(ll,2)*(gsupp(2,ii)-gsupm(2,ii))/dh/jac
cc     .                -gsub(ll,3)*(gsupp(3,ii)-gsupm(3,ii))/dh/jac
cc     .                +delta(ii,ll)*(gamma(1,kk,1)
cc     .                              +gamma(2,kk,2)
cc     .                              +gamma(3,kk,3))
cc
cc                summ=0d0
cc                do mm=1,3
cc                  do jj=1,3
cc                    summ= summ + gsub(ll,mm)*gsup(ii,jj)*gamma(mm,jj,kk)
cc                  enddo
cc                enddo
cc
cc                table(ii,kk,ll) = table(ii,kk,ll) - summ
cc
cc              enddo
cc            enddo
cc          enddo
cc
cc          !Symmetrize
cc          do ii=1,3             !cycle through Chistoffel symbols
cc            do ll=1,3
cc              do kk=ll,3
cc                table(ii,kk,ll) = hlf*(table(ii,kk,ll)+table(ii,ll,kk))
cc                table(ii,ll,kk) = table(ii,kk,ll)
cc              enddo
cc            enddo
cc          enddo
cc
cc          !Invert SI systems
cc          do ii=1,3
cc            do ll=1,3
cc              do kk=1,3
cc                if (ii==ll.and.ii/=kk) then
cc                  select case(ii)
cc                  case(1)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,kk,2)-table(3,kk,3)
cc                  case(2)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(3,kk,3)
cc                  case(3)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,kk,1)-table(2,kk,2)
cc                  end select
cc                elseif(ii==kk .and.ii/=ll) then
cc                  select case(ii)
cc                  case(1)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = -table(2,2,ll)-table(3,3,ll)
cc                  case(2)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(3,3,ll)
cc                  case(3)
cc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = -table(1,1,ll)-table(2,2,ll)
cc                  end select
cc                elseif(ii==kk .and.ii==ll) then
cc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     .                 = gamma(ii,kk,ll)
cc                else
cc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cc     $                 = table(ii,kk,ll)
cc                endif
cccc                if (ii == ll) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     .                 = matmul(inv,table(ii,:,ii))
cccccc     $                 = hlf*(table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccccc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cccccc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cccc                  case(2)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccccc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cccccc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cccc                  case(3)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccccc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cccccc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cccc                  end select
cccc                elseif (ii == kk) then
cccc                  select case(ii)
cccc                  case(1)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccccc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cccccc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cccc                  case(2)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccccc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cccccc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cccc                  case(3)
cccc                    gmetric%grid(igrid)%Gamma(i,j,k,ii,kk,ll)
cccc     $                 = hlf*(table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccccc     $                 =0.25*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll)
cccccc     $                       +table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cccc                  end select
cccc                else
cccc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cccc                endif
cc              enddo
cc            enddo
cc
cc          enddo

cc                if (ii /= ll .and. ii /= kk) then
cc                  gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
cccc     $                    = table(ii,ll,kk)
cc     $                    = hlf*( table(ii,ll,kk)+table(ii,kk,ll) )
cc                else
cc                  if (ii == ll) then
cc                    select case(ii)
cc                    case(1)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk))
cc     $                 =0.25*(table(1,1,kk)-table(2,2,kk)-table(3,3,kk)
cc     $                       +table(1,kk,1)-table(2,kk,2)-table(3,kk,3))
cc                    case(2)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk))
cc     $                 =0.25*(table(2,2,kk)-table(1,1,kk)-table(3,3,kk)
cc     $                       +table(2,kk,2)-table(1,kk,1)-table(3,kk,3))
cc                    case(3)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk))
cc     $                 =0.25*(table(3,3,kk)-table(1,1,kk)-table(2,2,kk)
cc     $                       +table(3,kk,3)-table(1,kk,1)-table(2,kk,2))
cc                    end select
cc                  else
cc                    select case(ii)
cc                    case(1)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll))
cc     $                 =0.25*(table(1,1,ll)-table(2,2,ll)-table(3,3,ll)
cc     $                       +table(1,ll,1)-table(2,ll,2)-table(3,ll,3))
cc                    case(2)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
ccc$$$     $                 = hlf*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll))
cc     $                 =0.25*(table(2,2,ll)-table(1,1,ll)-table(3,3,ll)
cc     $                       +table(2,ll,2)-table(1,ll,1)-table(3,ll,3))
cc                    case(3)
cc                      gmetric%grid(igrid)%Gamma(i,j,k,ii,ll,kk)
ccc$$$     $                 = hlf*(table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
ccc$$$     $                 = hlf*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll))
cc     $                 =0.25*(table(3,3,ll)-table(1,1,ll)-table(2,2,ll)
cc     $                       +table(3,ll,3)-table(1,ll,1)-table(2,ll,2))
cc                    end select
cc
cc                  endif
cc                endif
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif
        
c      START CHECKS

cc        call zeroForceCheck(i,j,k,igrid,gamma
cc     $                     ,gmetric%grid(igrid)%Gamma(i,j,k,:,:,:))

c     End program

      end subroutine gammaZeroForce

c     zeroForceCheck
c     #################################################################
      subroutine zeroForceCheck(i,j,k,igrid,gold,gnew)

c     -----------------------------------------------------------------
c     Postprocess Christoffel symbol to satisfy zero-force condition
c     to machine accuracy
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer(4) :: i,j,k,igrid
        real(8)    :: gold(3,3,3),gnew(3,3,3)
        
c     Local variables

        integer(4) :: ii,ll,kk,mm,jj,ig,jg,kg
        real(8)    :: summ,summ2,dh,jac,jacp,jacm,const,hlf
        real(8)    :: hess(3,3,3),gsup(3,3),gsub(3,3)
     $               ,gsupp(3,3),gsupm(3,3)
        logical    :: alt_eom

c     Begin program

        write (*,*)
        write (*,*) 'Grid:',igrid,'  Grid node:',i,j,k

        call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

        gsup  = gmetric%grid(igrid)%gsup(i,j,k,:,:)
        gsub  = gmetric%grid(igrid)%gsub(i,j,k,:,:)

        hlf = 0.5

c      Check difference between correction and prediction

        hess = gnew - gold

        summ2= sqrt(sum(gold*gold))
        summ = sqrt(sum(hess*hess))/summ2

        write (*,*) 'Correction size:',summ

c      Check cancellation property g^lk(d_il G^j_kj-g_lm g^ij G^m_jk) = 0

        do ii=1,3 !cycle through Chistoffel symbols
          do ll=1,3
            do kk=1,3
              hess(ii,ll,kk) =
     .               delta(ii,ll)*(gnew(1,kk,1)
     .                            +gnew(2,kk,2)
     .                            +gnew(3,kk,3))

              summ=0d0
              do mm=1,3
                do jj=1,3
                  summ = summ + gsub(ll,mm)*gsup(ii,jj)*gnew(mm,jj,kk)
                enddo
              enddo

              hess(ii,ll,kk) = hess(ii,ll,kk) - summ
            enddo
          enddo

          summ = 0d0
          do ll=1,3
            do kk=1,3
              summ = summ + gsup(ll,kk)*hess(ii,ll,kk)
            enddo
          enddo

          write (*,*) 'Cancellation property comp.#',ii,':',summ
        enddo

c      Check zero-force property: g^lk G^i_lk + nabla_m(g^mi) = 0

        do ii=1,3 !cycle through Chistoffel symbols

          alt_eom = (ii == 2 .and. alt__eom())

          if (alt_eom) then
            jac  = gmetric%grid(igrid)%jac(i,j,k)
            const= 1d0
          else
            jac  = 1d0
            const= 0d0
          endif

          summ = 0d0

          dh = grid_params%dxh(ig)
          if (alt_eom) then
            jacp = hlf*(jac+gmetric%grid(igrid)%jac(i+1,j,k))
            jacm = hlf*(jac+gmetric%grid(igrid)%jac(i-1,j,k))
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i+1,j,k,:,:))
          gsupm = jacm
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i-1,j,k,:,:))

          summ = summ +(gsupp(1,ii)-gsupm(1,ii))/dh/jac

          dh = grid_params%dyh(jg)
          if (alt_eom) then
            jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j+1,k))
            jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j-1,k))
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j+1,k,:,:))
          gsupm = jacm
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j-1,k,:,:))

          summ = summ +(gsupp(2,ii)-gsupm(2,ii))/dh/jac

          dh = grid_params%dzh(kg)
          if (alt_eom) then
            jacp = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k+1))
            jacm = hlf*(jac+gmetric%grid(igrid)%jac(i,j,k-1))
          else
            jacp=1d0
            jacm=1d0
          endif
          gsupp = jacp
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k+1,:,:))
          gsupm = jacm
     $         *hlf*(gsup+gmetric%grid(igrid)%gsup(i,j,k-1,:,:))

          summ = summ +(gsupp(3,ii)-gsupm(3,ii))/dh/jac

          summ2 = 0d0
          do ll=1,3
            do kk=1,3
              summ2 = summ2 + gsup(ll,kk)*gnew(ii,ll,kk)
     $                - const*gsup(ll,ii)*gnew(kk,ll,kk)
            enddo
          enddo

          write (*,*) 'Zero force comp. #',ii,':',summ,summ2,summ+summ2
cc          write (*,*) 'Zero force   property comp.#',ii,':',summ+summ2
        enddo

c     End program

      end subroutine zeroForceCheck

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
      function delta(i,j)

        integer(4) :: i,j
        real(8)    :: delta

        if (i == j) then
          delta = 1d0
        else
          delta = 0d0
        endif

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

        vec = (/ c1,c2,c3 /)

        if (covariant) then
          cx = dot_product(vec,gmetric%grid(igx)%cov(i,j,k,:,1))
          cy = dot_product(vec,gmetric%grid(igx)%cov(i,j,k,:,2))
          cz = dot_product(vec,gmetric%grid(igx)%cov(i,j,k,:,3))
        else
          cx = dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,:,1))
          cy = dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,:,2))
          cz = dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,:,3))
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

        vec = (/ cx,cy,cz /)

        if (covariant) then

          c1 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,1,:))
          c2 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,2,:))
          c3 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,3,:))
            
        else

          c1 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,1,:))
          c2 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,2,:))
          c3 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,3,:))

        endif

      end subroutine transformVectorToCurvilinear

c     transformFromCurvToCurv
c     #################################################################
      subroutine transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv,half_elem)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at grid coordinates (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        integer(4),optional,intent(IN) :: half_elem

        real(8)    :: cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical    :: tocnv

c     Local variables

        integer(4) :: he,ip
        real(8)    :: cov(3),cnv(3),tensor(3,3)

c     Begin program

        if (PRESENT(half_elem)) then
          he = half_elem
        else
          he = 0
        endif

        select case(he)
        case(1)
          ip = i+1
          if (isSP(ip,j,k,igx,igy,igz)) then !We are right @ SP ghost cell
            if (tocnv) then
              tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
            else
              tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
            endif
          elseif (isSP(i,j,k,igx,igy,igz)) then !We are right SP
            if (tocnv) then
              tensor = (gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .                 *gmetric%grid(igx)%jac (ip,j,k)
     .                 +gmetric%grid(igx)%gsup(i ,j,k,:,:)
     .                 *gmetric%grid(igx)%jac (i ,j,k))
     .                 /(gmetric%grid(igx)%jac(ip,j,k)
     .                  +gmetric%grid(igx)%jac(i ,j,k))
            else
              tensor = (gmetric%grid(igx)%gsub(ip,j,k,:,:)
     .                 *gmetric%grid(igx)%jac (ip,j,k)
     .                 +gmetric%grid(igx)%gsub(i ,j,k,:,:)
     .                 *gmetric%grid(igx)%jac (i ,j,k))
     .                 /(gmetric%grid(igx)%jac(ip,j,k)
     .                  +gmetric%grid(igx)%jac(i ,j,k))
            endif
          else
            if (tocnv) then
              tensor = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .                     +gmetric%grid(igx)%gsup(i ,j,k,:,:))
            else
              tensor = 0.5*(gmetric%grid(igx)%gsub(ip,j,k,:,:)
     .                     +gmetric%grid(igx)%gsub(i ,j,k,:,:))
            endif
          endif
        case(2)
          if (tocnv) then
            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j+1,k,:,:)
     .                   +gmetric%grid(igx)%gsup(i,j  ,k,:,:))
          else                                        
            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j+1,k,:,:)
     .                   +gmetric%grid(igx)%gsub(i,j  ,k,:,:))
          endif
        case(3)
          if (tocnv) then
            tensor = 0.5*(gmetric%grid(igx)%gsup(i,j,k+1,:,:)
     .                   +gmetric%grid(igx)%gsup(i,j,k  ,:,:))
          else                                          
            tensor = 0.5*(gmetric%grid(igx)%gsub(i,j,k+1,:,:)
     .                   +gmetric%grid(igx)%gsub(i,j,k  ,:,:))
          endif
        case(0)
          if (tocnv) then
            tensor = gmetric%grid(igx)%gsup(i,j,k,:,:)
          else
            tensor = gmetric%grid(igx)%gsub(i,j,k,:,:)
          endif
        end select

        if (tocnv) then
          cnv  = matmul(tensor,(/ cov1,cov2,cov3 /))
          cnv1 = cnv(1)
          cnv2 = cnv(2)
          cnv3 = cnv(3)
        else
          cov  = matmul(tensor,(/ cnv1,cnv2,cnv3 /))
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

cc        if (isSP2(i,j,k,igx,igy,igz)) vol = 0.25*(2+sqrt(2d0))*vol

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

        use error

        use grid_operations

        implicit none

        integer(4),private :: nxx,nyy,nzz

        integer(4) :: mg_ratio

        real(8),private :: pi

        type :: dim_pack
          logical    :: pack
cc          real(8)    :: xi
          real(8)    :: xp
          real(8)    :: dx0
cc          integer(4) :: n1
        end type dim_pack

        type :: grid_pack
          type(dim_pack),dimension(3) :: dim
        end type grid_pack

        type(grid_pack) :: g_pack

      contains

c     createGrid
c     #################################################################
      subroutine createGrid(imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
     .                     ,gpack)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities from:
c       * imin,imax: global grid limits in X
c       * jmin,jmax: global grid limits in Y
c       * kmin,kmax: global grid limits in Z
c       * nx,ny,nz: global grid size
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
        type(grid_pack),optional :: gpack

c     Local variables

        integer(4) :: i,n1,n2,n3
        integer(4) :: ngrdx,ngrdy,ngrdz
        type(grid_pack) :: g_pack

c     Begin program

        if (PRESENT(gpack)) then
          g_pack = gpack
        else
          g_pack%dim(:)%pack = .false.
        endif

        nxx = imax-imin+1
        nyy = jmax-jmin+1
        nzz = kmax-kmin+1

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
          grid_params%ilo(i) = grid_params%ilo(i-1)/mg_ratio + 1
cc          grid_params%ilo (i) = imin
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
          grid_params%jlo(i) = grid_params%jlo(i-1)/mg_ratio + 1
cc          grid_params%jlo (i) = jmin
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
          grid_params%klo(i) = grid_params%klo(i-1)/mg_ratio + 1
cc          grid_params%klo (i) = kmin
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

        call createLogicalGrid(grid_params%nxgl,grid_params%mg_ratio_x
     .                        ,grid_params%xx,grid_params%dx
     .                        ,grid_params%dxh,grid_params%nxv
     .                        ,grid_params%ngrid,grid_params%istartx
     .                        ,grid_params%ilo,grid_params%ihi
     .                        ,xmin,xmax,bcond(1),bcond(2)
     .                        ,g_pack%dim(1))

        call createLogicalGrid(grid_params%nygl,grid_params%mg_ratio_y
     .                        ,grid_params%yy,grid_params%dy
     .                        ,grid_params%dyh,grid_params%nyv
     .                        ,grid_params%ngrid,grid_params%istarty
     .                        ,grid_params%jlo,grid_params%jhi
     .                        ,ymin,ymax,bcond(3),bcond(4)
     .                        ,g_pack%dim(2))

        call createLogicalGrid(grid_params%nzgl,grid_params%mg_ratio_z
     .                        ,grid_params%zz,grid_params%dz
     .                        ,grid_params%dzh,grid_params%nzv
     .                        ,grid_params%ngrid,grid_params%istartz
     .                        ,grid_params%klo,grid_params%khi
     .                        ,zmin,zmax,bcond(5),bcond(6)
     .                        ,g_pack%dim(3))

c     Initialize MPI for BCs

#if defined(petsc)
        call initMPI
#endif

c     Store grid metric parameters in grid metric structure

        if (checkGridDatabase()) call defineGridMetric(gmetric)

      end subroutine createGrid

c     createLogicalGrid
c     #################################################################
      subroutine createLogicalGrid (nglobal,mgr,xx,dx,dxh,nx,ngrid
     .                             ,istart,ilo,ihi,lmin,lmax,bcs1,bcs2
     .                             ,dpack)

        implicit none

c     Call variables

        integer(4) :: ngrid,nglobal(ngrid),nx(ngrid),istart(ngrid)
     $               ,mgr(ngrid),ilo(ngrid),ihi(ngrid),bcs1,bcs2
        real(8)    :: xx(*),dx(*),dxh(*),lmin,lmax

        type(dim_pack) :: dpack

c     Local variables
        
        integer(4) :: ig,i,isig,ng
        real(8)    :: dh,length,lmax_loc,lmin_loc,lstart
     .               ,dx0,dxl,dxr,alpha1=1d0,alpha2=1d0
        logical    :: sing_point,spt1,spt2

c     Begin program

        spt1 = (bcs1 == SP)
        spt2 = (bcs1 == SP2)
        sing_point = spt1 .or. spt2

c     Find local domain limits lmin_loc,lmax_loc

        ig = 1  !Finest grid level

c     Find local length

        dx0 = dpack%dx0
        dxl = dpack%dx0
        dxr = dpack%dx0

c     Find grid quantities

        do ig = 1,ngrid

          isig = istart(ig)

          ng   = nglobal(ig)

          !Find grid spacing
          call findLocalGrid(ig,xx(isig),ng)

          !Find integer mesh spacings
          do i = 0,nx(ig)
            dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
          enddo

          !Find half mesh spacings
          if (sing_point .and. ilo(1) == 1) then
            dxh(1 + isig) = 0.5*(xx(2+isig)+xx(1+isig))
          else
            dxh(1 + isig) = 0.5*(dx(1 + isig) + dx(0   + isig))
          endif

          do i = 2,nx(ig)
            dxh(i + isig) = 0.5*(dx(i + isig) + dx(i-1 + isig))
          enddo

          if (bcs1 == PER) then
            dxh(0        + isig) = dxh(nx(ig) + isig)
            dxh(nx(ig)+1 + isig) = dxh(1      + isig)
          elseif (sing_point .and. ilo(1) == 1) then
            dxh(0        + isig) = dx(0      + isig)
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)*0.5
          else
            dxh(0        + isig) = dx(0      + isig)*0.5
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)*0.5
          endif

        enddo

      contains

c     findLocalGrid
c     #################################################################
      subroutine findLocalGrid(ig,xx,ng)

c     -----------------------------------------------------------------
c     Finds grid on local processor
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: ig,ng
        real(8)    :: xx(*)

c     Local variables

        real(8)    :: xg(0:ng+1)

c     Begin program

c     Find global grid

        call findGlobalGrid(ig,xg,ng)

c     Assign local grid

        xx(1:nx(ig)+2) = xg(ilo(ig)-1:ihi(ig)+1)

c     End program

      end subroutine findLocalGrid

c     findGlobalGrid
c     #################################################################
      subroutine findGlobalGrid(ig,xg,ng)

c     -----------------------------------------------------------------
c     Finds global grid
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: ig,ng
        real(8)    :: xg(0:ng+1)

c     Local variables

        integer(4) :: n1,n2
        real(8)    :: dh,length,lstart,lend,xp,xi,xil,xir
        logical    :: pack_l,pack_r

c     Begin program

        length = lmax-lmin

        if (dpack%pack) then

          dh = length/ng
          lstart = lmin
          lend   = lmax

          xp = dpack%xp

          pack_l = .false.
          pack_r = .false.

          if (xp < lend .and. xp > lstart) then

            n1 = (xp - lmin)/length*ng

            n2 = ng - n1

            dxl = alpha1*dxl  !Modify initial interval for coarser grids

cc            call packGrid(ig,n1,lstart,xp,-dxl,xg(1:n1),alpha1,xi)
            call packGrid(ig,n1,lstart,xp-dxl/2,-dxl,xg(1:n1),alpha1
     .                   ,xil)

            dxr = alpha2*dxr  !Modify initial interval for coarser grids

cc            call packGrid(ig,n2,xp+dxr,lend,dxr,xg(n1+1:ng),alpha2,xi)
            call packGrid(ig,n2,xp+dxr/2,lend,dxr,xg(n1+1:ng),alpha2
     .                   ,xir)

          else

            dx0 = alpha1*dx0  !Modify initial interval for coarser grids

            if (xp <= lstart) then  !Pack to the left
              pack_l = .true.
              call packGrid(ig,ng,lstart,lend, dx0,xg(1:ng),alpha1,xi)
            else                    !Pack to the right
              pack_r = .true.
              call packGrid(ig,ng,lstart,lend,-dx0,xg(1:ng),alpha1,xi)
            endif

            xil = xi
            xir = xi

          endif

          !X0 ghost cell
          if (sing_point .and. ilo(1) == 1) then
            xg(0) = 1d-8*dh
          elseif (bcs1 == PER) then
            xg(0) = xg(1) - (xg(ng)-xg(ng-1))
          elseif (abs(bcs1) == SYM) then
            xg(0) = xg(1) - (xg(2 )-xg(1   ))
          else
            if (pack_l) then
              xg(0) = xg(1) - (xg(2 )-xg(1   ))
            else
              xg(0) = xg(1) - (xg(2 )-xg(1   ))*xil
            endif
          endif

          !Xn ghost cell
          if (bcs2 == PER) then
            xg(ng+1) = xg(ng) + (xg(2) -xg(1)   )
          elseif (abs(bcs2) == SYM) then
            xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))
          else
            if (pack_r) then
              xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))
            else
              xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))*xir
            endif
          endif

cc          write (*,*) xp,lmin,lmax,lstart,lend,dh,n1,n2,ng
cc          write (*,*) xg
cc          stop

        else  !Uniform grid

          dh = length/ng

          if (bcs1 == PER .and. bcs2 == PER) then
            dh = length/ng
            lstart = lmin
            lend   = lmax
          elseif (bcs1 == SYM .and. bcs2 == SYM) then
            write (*,*) 'here'
            dh = length/(ng-1)
            lstart = lmin
            lend   = lmax
          elseif (bcs1 == SYM .and. bcs2 /= SYM) then
            dh = length/(ng - 0.5)
            lstart = lmin
            lend   = lmax-dh/2.
          elseif (bcs1 /= SYM .and. bcs2 == SYM) then
            dh = length/(ng - 0.5)
            lstart = lmin+dh/2.
            lend   = lmax
          elseif (spt2) then
            dh = sqrt(2d0)*length/(1+sqrt(2d0)*(ng - 0.5))
            lstart = lmin+dh/sqrt(2d0)
            lend   = lmax-dh/2.
          else
            dh = length/ng
            lstart = lmin+dh/2.
            lend   = lmax-dh/2.
          endif

          xg(1) = lstart
          do i = 2,ng+1
            xg(i) = xg(i-1) + dh
          enddo
          if (sing_point) then
            xg(0) = 1d-8*dh
          else
            xg(0) = xg(1) - dh
          endif

        endif

c     End program

      end subroutine findGlobalGrid

c     packGrid
c     #################################################################
      subroutine packGrid(ig,nn,lmin,lmax,dx0,xx,alpha,xi)

c     -----------------------------------------------------------------
c     Finds packed grid, packing to the right if dx0<0 and to the 
c     left if dx0>0.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nn,ig
        real(8)    :: lmin,lmax,dx0,xx(nn),alpha,xi

c     Local variables

        integer(4) :: i
        real(8)    :: rhs

c     Begin program

        rhs = (lmax-lmin)/abs(dx0)

        call find_xi(nn-1,rhs,xi)

        if (dx0 > 0d0) then !Pack to the left

          xx(1) = lmin
          do i = 2,nn
            xx(i) = xx(i-1) + dx0*xi**(i-2)
          enddo
            
        else              !Pack to the right

          xx(nn) = lmax
          do i = nn,2,-1
            xx(i-1) = xx(i) - abs(dx0)*xi**(nn-i)
          enddo

        endif

        if (mgr(ig) > 1 ) then
          alpha = (xi**mgr(ig)-1)/(xi-1)
        else
          alpha = 1
        endif

c     End program

      end subroutine packGrid

c     find_xi
c     #################################################################
      subroutine find_xi(nn,rhs,xi)

c     -----------------------------------------------------------------
c     Finds geometric factor xi using Newton's method according to the
c     formula:
c              (xi^nn-1)/(xi-1) = rhs
c     The Newton iteration is initialized by using xi=1+eps, and
c     expanding in eps to second order in the formula above to find:
c              eps_0 ~ 2*(rhs-nn)/nn^2
c     This result is used to initialize a single Picard iteration,
c              xi_0 ~ exp(log(1+eps_0*rhs)/nn)
c     which is used as the first guess for Newton.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nn
        real(8)    :: rhs,xi

c     Local variables

        integer(4) :: i,nit=20
        real(8)    :: ijac,dxi,res

c     Begin program

cc        write (*,*) nn,rhs

        !Initial guess
        xi = 2*(rhs-nn)/nn**2  !Initial guess
        xi = exp(log(1+xi*rhs)/nn)

cc        write (*,*) 'xi0',xi

        do i=1,nit
cc          res = rhs-(xi**nn-1)/(xi-1)
cc          ijac = (xi-1)**2/((nn-1)*xi**nn-nn*xi**(nn-1)+1)
          res = rhs-0.5*(xi**nn+xi**(nn+1)-2)/(xi-1)
          ijac = 2*(xi-1)**2/(nn*(xi**(nn+1)-xi**(nn-1))-2*xi**nn+2)

cc          write (*,*) 'Newton it.',i,' res. in find_xi',abs(res)

          dxi = res*ijac
          xi = xi + dxi
cc          write (*,*) 'xi',xi
          if (abs(res) < 1d-10) exit
        enddo

        if (i >= nit) then
          write (*,*) 'Newton it. in find_xi did not converge'
          write (*,*) 'Aborting...'
          stop
        endif

c     End program

      end subroutine find_xi

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
        case ('cyl','hel','hl2')
          if (xmin /= 0d0 .and. bcSP()) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif
        case ('tor')
          if (xmin /= 0d0 .and. bcSP()) then
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

c     Local variables

        integer(4) :: i,j,n1,n2,n3
        integer(4) :: igx,igy,igz,isig

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

        call metricTensorCheck(igx,igy,igz
     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)

cc        call hessianCheck(igx,igy,igz
cc     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)

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

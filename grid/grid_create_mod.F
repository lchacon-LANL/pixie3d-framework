
c module grid_create
c #####################################################################
      module grid_create

        use grid_metric_init

        implicit none

        integer(4),private :: nxx,nyy,nzz

        integer(4) :: mg_ratio

        real(8),private :: pi

        type :: dim_pack
          logical    :: pack = .false.
          real(8)    :: xp = 0d0
          real(8)    :: dx0 = 0d0
        end type dim_pack

        type :: grid_pack
          type(dim_pack),dimension(3) :: dim
        end type grid_pack

        type(grid_pack),save :: g_pack

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

        !Check for packing
        if (PRESENT(gpack)) then
          g_pack = gpack
        else
          g_pack%dim(:)%pack = .false.
        endif

        !Initialize global quantities
        nxgl = nx
        nygl = ny
        nzgl = nz

        !Initialize local quantities
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
     .                     ,grid_params%xg,grid_params%xx,grid_params%dx
     .                     ,grid_params%dxh,grid_params%nxv
     .                     ,grid_params%ngrid,grid_params%istartx
     .                     ,grid_params%ilo,grid_params%ihi
     .                     ,xmin,xmax,bcond(1),bcond(2)
     .                     ,g_pack%dim(1))

        call createLogicalGrid(grid_params%nygl,grid_params%mg_ratio_y
     .                     ,grid_params%yg,grid_params%yy,grid_params%dy
     .                     ,grid_params%dyh,grid_params%nyv
     .                     ,grid_params%ngrid,grid_params%istarty
     .                     ,grid_params%jlo,grid_params%jhi
     .                     ,ymin,ymax,bcond(3),bcond(4)
     .                     ,g_pack%dim(2))

        call createLogicalGrid(grid_params%nzgl,grid_params%mg_ratio_z
     .                     ,grid_params%zg,grid_params%zz,grid_params%dz
     .                     ,grid_params%dzh,grid_params%nzv
     .                     ,grid_params%ngrid,grid_params%istartz
     .                     ,grid_params%klo,grid_params%khi
     .                     ,zmin,zmax,bcond(5),bcond(6)
     .                     ,g_pack%dim(3))

c     Initialize check grid output

        if (check_grid) then
          write (*,*)
          write (*,*) 'Coordinate system: ',coords
          write (*,*)
          write (*,*) 'Number of grid levels:'
          write (*,*) 'nx: ',grid_params%ngrid
     .            ,'   ny: ',grid_params%ngrid
     .            ,'   nz: ',grid_params%ngrid
        endif

c     Create fortran PETSc parallel grid hierarchy

#if defined(petsc)
      call createPETScGrid
#endif

c     Store grid metric parameters in grid metric structure

        if (checkMapDatabase()) call defineGridMetric(gmetric)

c     Perform grid checks

        if (check_grid) call checkGrid

      end subroutine createGrid

c     createLogicalGrid
c     #################################################################
      subroutine createLogicalGrid (nglobal,mgr,xg,xx,dx,dxh,nx,ngrid
     .                             ,istart,ilo,ihi,lmin,lmax,bcs1,bcs2
     .                             ,dpack)

        implicit none

c     Call variables

        integer(4) :: ngrid,nglobal(ngrid),nx(ngrid),istart(ngrid)
     $               ,mgr(ngrid),ilo(ngrid),ihi(ngrid),bcs1,bcs2
        real(8)    :: xg(*),xx(*),dx(*),dxh(*),lmin,lmax

        type(dim_pack) :: dpack

c     Local variables
        
        integer(4) :: ig,i,isig,ng
        real(8)    :: dh,length,lmax_loc,lmin_loc,lstart
     .               ,dx0,dxl,dxr,alpha1=1d0,alpha2=1d0

c     Begin program

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
          dxh(1 + isig) = 0.5*(dx(1 + isig) + dx(0   + isig))

          do i = 2,nx(ig)
            dxh(i + isig) = 0.5*(dx(i + isig) + dx(i-1 + isig))
          enddo

          if (bcs1 == PER) then
            dxh(0        + isig) = dxh(nx(ig) + isig)
            dxh(nx(ig)+1 + isig) = dxh(1      + isig)
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

        real(8)    :: xxg(0:ng+1)

c     Begin program

c     Find global grid

        call findGlobalGrid(ig,xxg,ng)

c     Assign global grid

        if (ig == 1) xg(1:ng+2) = xxg(0:ng+1)

c     Assign local grid

        xx(1:nx(ig)+2) = xxg(ilo(ig)-1:ihi(ig)+1)

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
        logical    :: pack_l,pack_r,fv_l,fv_r

c     Begin program

        length = lmax-lmin

        fv_l = .not.(abs(bcs1) == SYM)  !FV discr. at opposite boundary except for SYM BC
        fv_r = .not.(abs(bcs2) == SYM)  !FV discr. at opposite boundary except for SYM BC

        if (dpack%pack) then

          dh = length/ng
          lstart = lmin
          lend   = lmax

          xp = dpack%xp

          pack_l = .false.
          pack_r = .false.

          if (xp < lend .and. xp > lstart) then

            n1 = max(int((xp - lmin)/length*ng),2)

            n2 = ng - n1

            dxl = alpha1*dxl  !Modify initial interval for coarser grids

cc            call packGrid(ig,n1,lstart,xp,-dxl,xg(1:n1),alpha1,xi)
            call packGrid(ig,n1,lstart,xp-dxl/2,-dxl,xg(1:n1),alpha1
     .                   ,xil,fv_l)

            dxr = alpha2*dxr  !Modify initial interval for coarser grids

cc            call packGrid(ig,n2,xp+dxr,lend,dxr,xg(n1+1:ng),alpha2,xi)
            call packGrid(ig,n2,xp+dxr/2,lend,dxr,xg(n1+1:ng),alpha2
     .                   ,xir,fv_r)

          else

            dx0 = alpha1*dx0  !Modify initial interval for coarser grids

            if (xp <= lstart) then  !Pack to the left
              pack_l = .true.
              call packGrid(ig,ng,lstart,lend, dx0,xg(1:ng),alpha1,xi
     .                     ,fv_r)
            else                    !Pack to the right
              pack_r = .true.
              call packGrid(ig,ng,lstart,lend,-dx0,xg(1:ng),alpha1,xi
     .                     ,fv_l)
            endif

            xil = xi
            xir = xi

          endif

          !X0 ghost cell
          if (bcs1 == PER) then
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

        else  !Uniform grid

          dh = length/ng

          if (bcs1 == PER .and. bcs2 == PER) then
            dh = length/ng
            lstart = lmin
            lend   = lmax
          elseif (bcs1 == SYM .and. bcs2 == SYM) then
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
          else
            dh = length/ng
            lstart = lmin+dh/2.
            lend   = lmax-dh/2.
          endif

          xg(1) = lstart
          do i = 2,ng+1
            xg(i) = xg(i-1) + dh
          enddo
cc          if (bcSP()) then
cc            xg(0) = xg(1) - dh + 1d-10
cc          else
            xg(0) = xg(1) - dh
cc          endif

        endif

c     End program

      end subroutine findGlobalGrid

c     packGrid
c     #################################################################
      subroutine packGrid(ig,nn,lmin,lmax,dx0,xx,alpha,xi,fv)

c     -----------------------------------------------------------------
c     Finds packed grid, packing to the right if dx0<0 and to the 
c     left if dx0>0. If fv=.true., we pack so that the grid in the
c     unpacked boundary meets the boundary at a face. Otherwise, we
c     pack so that a grid point meets that boundary.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nn,ig
        real(8)    :: lmin,lmax,dx0,xx(nn),alpha,xi
        logical    :: fv

c     Local variables

        integer(4) :: i
        real(8)    :: rhs

c     Begin program

        rhs = (lmax-lmin)/abs(dx0)

        call find_xi(nn-1,rhs,xi,fv)

        if (dx0 > 0d0) then !Pack to the left

          xx(1) = lmin
          do i = 2,nn
            xx(i) = xx(i-1) + dx0*xi**(i-2)
          enddo
            
        else                !Pack to the right

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
      subroutine find_xi(nn,rhs,xi,fv)

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
        logical    :: fv

c     Local variables

        integer(4),parameter :: nit=20
        integer(4) :: i
        real(8)    :: ijac,dxi(nit),res(nit)

c     Begin program

c     Initial guess

        xi = 2*(rhs-nn)/nn**2
        xi = exp(log(1+xi*rhs)/nn)

c     Newton iteration

        do i=1,nit
          if (fv) then  !Do packing in FV fashion at unpacked boundary
            res(i) = rhs-0.5*(xi**nn+xi**(nn+1)-2)/(xi-1)
            ijac = 2*(xi-1)**2/(nn*(xi**(nn+1)-xi**(nn-1))-2*xi**nn+2)
          else          !Position both grid points at both boundaries
            res(i) = rhs-(xi**nn-1)/(xi-1)
            ijac = (xi-1)**2/((nn-1)*xi**nn-nn*xi**(nn-1)+1)
          endif

          dxi(i) = res(i)*ijac
          xi = xi + dxi(i)

          if (abs(res(i)) < 1d-10 .or. (abs(dxi(i)) < 1d-15) ) exit
        enddo

c     Convergence check

        if (i >= nit) then
          write (*,*)
          write (*,*) 'Newton it. in grid-packing routine find_xi'
     .               ,' did not converge'
          write (*,*)
          write (*,*) 'Convergence history:'
          do i=1,nit
            write (*,*) 'Newton it.',i,' residual',abs(res(i))
     .                 ,'update',abs(dxi(i))
          enddo
          write (*,*)
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
          numerical_grid = .false.
        case ('scl')
          numerical_grid = .true.
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
c     Dumps grid info and checks metric tensors.
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

        call gridInfo('X',igx,grid_params%istartx,grid_params%nxv
     .               ,grid_params%xx,grid_params%dx,grid_params%dxh)

        call gridInfo('Y',igy,grid_params%istarty,grid_params%nyv
     .               ,grid_params%yy,grid_params%dy,grid_params%dyh)

        call gridInfo('Z',igz,grid_params%istartz,grid_params%nzv
     .               ,grid_params%zz,grid_params%dz,grid_params%dzh)

        if (checkMapDatabase()) then
          call metricTensorCheck(1,1,1
cc          call metricTensorCheck(igx,igy,igz
     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)
        endif

cc        if (checkAnalMapDatabase()) then
cc          call hessianCheck(igx,igy,igz
cc     .               ,grid_params%nxv,grid_params%nyv,grid_params%nzv)
cc        endif

        call pstop('checkGrid','End of grid check')

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
        write (*,*) 'Grid levels:',igx,igy,igz
        write (*,*)

        do i = 1,nxv(igx)   !Start at 2 to avoid singular points
          do j = 1,nyv(igy)
            do k = 1,nzv(igz)

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = grid_params%xx(ig)
              y1 = grid_params%yy(jg)
              z1 = grid_params%zz(kg)

              gup   = gmetric%grid(igx)%gsup(i,j,k,:,:)
              gdown = gmetric%grid(igx)%gsub(i,j,k,:,:)

              tensor = matmul(gup,gdown)

cc              write (*,5) 'Grid point: (',x1,',',y1,',',z1,')'
cc 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)
              write (*,5) 'Grid point: (',i,',',j,',',k,')'
 5            format (/,a,i3,a,i3,a,i3,a)

              write (*,*)
              write (*,*) 'Gup'
              write (*,10) gup(1,1:3)
              write (*,10) gup(2,1:3)
              write (*,10) gup(3,1:3)

              write (*,*)
              write (*,*) 'Gdown'
              write (*,10) gdown(1,1:3)
              write (*,10) gdown(2,1:3)
              write (*,10) gdown(3,1:3)

              write (*,*)
              write (*,*) 'Cov: cov1      cov2      cov3'
              write (*,10) gmetric%grid(igx)%cov(i,j,k,1:3,1)
              write (*,10) gmetric%grid(igx)%cov(i,j,k,1:3,2)
              write (*,10) gmetric%grid(igx)%cov(i,j,k,1:3,3)

              write (*,*)
              write (*,*) 'Cnv: cnv1      cnv2      cnv3'
              write (*,10) gmetric%grid(igx)%cnv(i,j,k,1:3,1)
              write (*,10) gmetric%grid(igx)%cnv(i,j,k,1:3,2)
              write (*,10) gmetric%grid(igx)%cnv(i,j,k,1:3,3)

              write (*,*)
              write (*,*) 'Jac:',gmetric%grid(igx)%jac(i,j,k)

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

      end module grid_create

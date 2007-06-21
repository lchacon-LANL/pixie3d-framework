
c module grid_create
c #####################################################################
      module grid_create

        use grid_metric_init

        implicit none

        integer    :: mg_ratio

        type :: dim_pack
          logical    :: pack
          real(8)    :: xp
          real(8)    :: dx0
        end type dim_pack

        type :: grid_pack
          type(dim_pack),dimension(3) :: dim
        end type grid_pack

        type(grid_pack),save :: g_pack

        private :: initMG,createLogicalGrid
     $            ,consistencyCheck,checkGrid

      contains

c     createGrid
c     #################################################################
      subroutine createGrid(nx,ny,nz,gpack)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities from:
c       * nx,ny,nz: global grid size
c       * gpack: defines grid packing
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nx,ny,nz
        type(grid_pack) :: gpack

c     Local variables

        integer    :: imin,imax,jmin,jmax,kmin,kmax

c     Begin program

        pi = acos(-1d0)

c     Initialize MPI and PETSc (parallel version)

#if defined(petsc)
      call initFortranMPI(nx,ny,nz)
cc#if !defined(petsc_c)
      call setASMflag(.false.)
cc#endif
#endif

c     Initialize MG

        call initMG(nx,ny,nz,grid_params)

c     Create MG grid hierarchy

        call createMGgrid(gpack,grid_params)

c     Create fortran PETSc parallel grid hierarchy (parallel version)

#if defined(petsc)
        call createPETScGrid
#endif

c     Store grid metric parameters in grid metric structure

        if (checkMapDatabase()) call defineGridMetric(gmetric)

c     Perform grid checks

        if (check_grid) call checkGrid

      end subroutine createGrid

c     initMG
c     #################################################################
      subroutine initMG(nx,ny,nz,gdef)

c     -----------------------------------------------------------------
c     Initializes MG quantities and allocates MG strucutures. On input:
c       * nx,ny,nz: global grid size
c       * gdef (type grid_mg_def): variable that stores grid definition.
c     In parallel, local grid sizes are read from PETSc's DA context.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nx,ny,nz
        type(grid_mg_def) :: gdef

c     Local variables

        integer    :: i,n1,n2,n3,nxx,nyy,nzz
        integer    :: ngrdx,ngrdy,ngrdz

c     Begin program

        !Initialize local quantities
#if defined(petsc)
        nxx = dactx(1)%xm
        nyy = dactx(1)%ym
        nzz = dactx(1)%zm
#else
        nxx = nx
        nyy = ny
        nzz = nz
#endif

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
     .                            ,nx,ny,nz,gdef)

      end subroutine initMG

c     createMGgrid
c     #################################################################
      subroutine createMGgrid(gpack,gdef)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities from:
c       * gpack (type grid_pack): stores grid packing info
c       * gdef (type grid_mg_def): stores grid definition.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
        type(grid_pack) :: gpack
        type(grid_mg_def) :: gdef

c     Local variables

        integer    :: i,n1,n2,n3,nxx,nyy,nzz
        integer    :: ngrdx,ngrdy,ngrdz
        type(grid_pack) :: g_pack

c     Begin program

        !Save packing configuration
        g_pack = gpack

        !Initialize local quantities
#if defined(petsc)
        nxx = dactx(1)%xm
        nyy = dactx(1)%ym
        nzz = dactx(1)%zm
#else
        nxx = gdef%nglx
        nyy = gdef%ngly
        nzz = gdef%nglz
#endif

        !Grid levels
        ngrdx = gdef%ngrdx
        ngrdy = gdef%ngrdy
        ngrdz = gdef%ngrdz

c     Initialize MG grid-size arrays

        gdef%mg_ratio_x = 1
        gdef%mg_ratio_y = 1
        gdef%mg_ratio_z = 1

        gdef%nxgl(1) = gdef%nglx
#if defined(petsc)
        gdef%nxv (1) = dactx(1)%xm
        gdef%ilo (1) = dactx(1)%xs
        gdef%ihi (1) = dactx(1)%xe
#else
        gdef%nxv (1) = nxx
        gdef%ilo (1) = 1
        gdef%ihi (1) = nxx
#endif
        do i = 2,ngrdx
          gdef%nxgl(i) = gdef%nxgl(i-1) / mg_ratio
          gdef%nxv (i) = gdef%nxv (i-1) / mg_ratio
          gdef%ilo (i) = gdef%ilo(i-1)/mg_ratio + 1
          gdef%ihi (i) = gdef%ilo(i)+gdef%nxv(i)-1
          gdef%mg_ratio_x(i-1) = mg_ratio
        enddo
        do i = ngrdx+1,gdef%ngrid
          gdef%nxgl(i) = gdef%nxgl(i-1)
          gdef%nxv (i) = gdef%nxv (i-1)
          gdef%ilo (i) = gdef%ilo (i-1)
          gdef%ihi (i) = gdef%ihi (i-1)
        enddo

        gdef%nygl(1) = gdef%ngly
#if defined(petsc)
        gdef%nyv (1) = dactx(1)%ym
        gdef%jlo (1) = dactx(1)%ys
        gdef%jhi (1) = dactx(1)%ye
#else
        gdef%nyv (1) = nyy
        gdef%jlo (1) = 1
        gdef%jhi (1) = nyy
#endif
        do i = 2,ngrdy
          gdef%nygl(i) = gdef%nygl(i-1) / mg_ratio
          gdef%nyv (i) = gdef%nyv(i-1) / mg_ratio
          gdef%jlo (i) = gdef%jlo(i-1)/mg_ratio + 1
          gdef%jhi (i) = gdef%jlo(i)+gdef%nyv(i)-1
          gdef%mg_ratio_y(i-1) = mg_ratio
        enddo
        do i = ngrdy+1,gdef%ngrid
          gdef%nygl(i) = gdef%nygl(i-1)
          gdef%nyv (i) = gdef%nyv (i-1)
          gdef%jlo (i) = gdef%jlo (i-1)
          gdef%jhi (i) = gdef%jhi (i-1)
        enddo

        gdef%nzgl(1) = gdef%nglz
#if defined(petsc)
        gdef%nzv (1) = dactx(1)%zm
        gdef%klo (1) = dactx(1)%zs
        gdef%khi (1) = dactx(1)%ze
#else
        gdef%nzv (1) = nzz
        gdef%klo (1) = 1
        gdef%khi (1) = nzz
#endif
        do i = 2,ngrdz
          gdef%nzgl(i) = gdef%nzgl(i-1) / mg_ratio
          gdef%nzv (i) = gdef%nzv (i-1) / mg_ratio
          gdef%klo (i) = gdef%klo(i-1)/mg_ratio + 1
          gdef%khi (i) = gdef%klo(i)+gdef%nzv(i)-1
          gdef%mg_ratio_z(i-1) = mg_ratio
        enddo
        do i = ngrdz+1,gdef%ngrid
          gdef%nzgl(i) = gdef%nzgl(i-1)
          gdef%nzv (i) = gdef%nzv (i-1)
          gdef%klo (i) = gdef%klo (i-1)
          gdef%khi (i) = gdef%khi (i-1)
        enddo

c     Initialize MG pointer arrays

        gdef%istartx(1) = 1
        do i = 2,gdef%ngrid
          gdef%istartx(i) = gdef%istartx(i-1)
     .                            +(gdef%nxv(i-1)+2)
        enddo

        gdef%istarty(1) = 1
        do i = 2,gdef%ngrid
          gdef%istarty(i) = gdef%istarty(i-1)
     .                            +(gdef%nyv(i-1)+2)
        enddo

        gdef%istartz(1) = 1
        do i = 2,gdef%ngrid
          gdef%istartz(i) = gdef%istartz(i-1)
     .                            +(gdef%nzv(i-1)+2)
        enddo

        gdef%istartp(1) = 1
        do i = 2,gdef%ngrid
          gdef%istartp(i) = gdef%istartp(i-1)
     .                            +gdef%nxv(i-1)
     .                            *gdef%nyv(i-1)
     .                            *gdef%nzv(i-1)
        enddo

c     Set grid parameters

        gdef%params = gparams

c     Consistency checks

        call consistencyCheck(nxx,nyy,nzz)

c     Define uniform logical LOCAL grid on ALL grid levels

        call createLogicalGrid(gdef%nxgl,gdef%mg_ratio_x
     .                     ,gdef%xg,gdef%xx,gdef%dx
     .                     ,gdef%dxh,gdef%nxv
     .                     ,gdef%ngrid,gdef%istartx
     .                     ,gdef%ilo,gdef%ihi
     .                     ,xmin,xmax,bcond(1),bcond(2)
     .                     ,g_pack%dim(1))

        call createLogicalGrid(gdef%nygl,gdef%mg_ratio_y
     .                     ,gdef%yg,gdef%yy,gdef%dy
     .                     ,gdef%dyh,gdef%nyv
     .                     ,gdef%ngrid,gdef%istarty
     .                     ,gdef%jlo,gdef%jhi
     .                     ,ymin,ymax,bcond(3),bcond(4)
     .                     ,g_pack%dim(2))

        call createLogicalGrid(gdef%nzgl,gdef%mg_ratio_z
     .                     ,gdef%zg,gdef%zz,gdef%dz
     .                     ,gdef%dzh,gdef%nzv
     .                     ,gdef%ngrid,gdef%istartz
     .                     ,gdef%klo,gdef%khi
     .                     ,zmin,zmax,bcond(5),bcond(6)
     .                     ,g_pack%dim(3))

c     Store physical dimensions of local domain

        gdef%lxmin = gdef%xg(gdef%ilo(1))
        gdef%lxmax = gdef%xg(gdef%ihi(1))

        gdef%lymin = gdef%yg(gdef%jlo(1))
        gdef%lymax = gdef%yg(gdef%jhi(1))

        gdef%lzmin = gdef%zg(gdef%klo(1))
        gdef%lzmax = gdef%zg(gdef%khi(1))

c     Initialize grid check output

        if (check_grid) then
          write (*,*)
          write (*,*) 'Coordinate system: ',coords
          write (*,*)
          write (*,*) 'Number of grid levels:'
          write (*,*) 'nx: ',gdef%ngrid
     .            ,'   ny: ',gdef%ngrid
     .            ,'   nz: ',gdef%ngrid
        endif

      end subroutine createMGgrid

c     createLogicalGrid
c     #################################################################
      subroutine createLogicalGrid (nglobal,mgr,xg,xx,dx,dxh,nx,ngrid
     .                             ,istart,ilo,ihi,lmin,lmax,bcs1,bcs2
     .                             ,dpack)

c     -----------------------------------------------------------------
c     Creates logical grid parameters (positions, spacing, etc.) using
c     BC info. 
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: ngrid,nglobal(ngrid),nx(ngrid),istart(ngrid)
     $               ,mgr(ngrid),ilo(ngrid),ihi(ngrid),bcs1,bcs2
        real(8)    :: xg(*),xx(*),dx(*),dxh(*),lmin,lmax

        type(dim_pack) :: dpack

c     Local variables
        
        integer    :: ig,i,isig,ng
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
            dxh(0        + isig) = dx(0      + isig)*0.5 !Needed for extrapolation
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)*0.5 !Needed for extrapolation
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

        integer    :: ig,ng
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

        integer    :: ig,ng
        real(8)    :: xg(0:ng+1)

c     Local variables

        integer    :: n1,n2
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

        integer    :: nn,ig
        real(8)    :: lmin,lmax,dx0,xx(nn),alpha,xi
        logical    :: fv

c     Local variables

        integer    :: i
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

        integer    :: nn
        real(8)    :: rhs,xi
        logical    :: fv

c     Local variables

        integer   ,parameter :: nit=20
        integer    :: i
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
      subroutine consistencyCheck(nxx,nyy,nzz)

c     -----------------------------------------------------------------
c     Checks consistency of grid parameters
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: nxx,nyy,nzz

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
          numerical_grid = .true.   !To deal properly with ghost cells in atanh profile
cc          numerical_grid = .false.
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

          major_r = gparams(1)

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
          if (nxx == 1 .and. ((xmax-xmin)>1d-3)) then
            write (*,*)
            write (*,*) 'Warning: large grid spacing in',
     .                  ' ignorable direction'
            write (*,*)
          endif
          if (nyy == 1 .and. ((ymax-ymin)>1d-3)) then
            write (*,*)
            write (*,*) 'Warning: large grid spacing in',
     .                  ' ignorable direction'
            write (*,*)
          endif
          if (nzz == 1 .and. ((zmax-zmin)>1d-3)) then
            write (*,*)
            write (*,*) 'Warning: large grid spacing in',
     .                  ' ignorable direction'
            write (*,*)
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

        integer    :: i,j,n1,n2,n3
        integer    :: igx,igy,igz,isig

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
        integer    :: ig,istart(*),nxv(*)
        real(8)    :: xx(*),dx(*),dxh(*)

c     Local variables

        integer    :: i,j,isig

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

        integer    :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

c     Local variables

        integer    :: i,j,k,i1,j1,k1,ig,jg,kg,ih
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

        integer    :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

c     Local variables

        integer    :: i,j,k,i1,j1,k1,ig,jg,kg,ih
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

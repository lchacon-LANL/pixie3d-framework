
c module grid_create
c #####################################################################
      module grid_create

        use grid_metric_init

        implicit none

        integer    :: mg_ratio

        private :: initGridSt,createLogicalGrid
     $            ,consistencyCheck,checkMGGrid

      contains

c     createGrid
c     #################################################################
      subroutine createGrid(nx,ny,nz,gpack,gdef)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities from:
c       * nx,ny,nz: global grid size
c       * gpack: defines grid packing
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_pack),dimension(3) :: gpack

        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer :: imin,imax,jmin,jmax,kmin,kmax
        integer :: nx,ny,nz,mpi_comm

c     Begin program

        pi = acos(-1d0)

c     Initialize global quantities (from gdef, if existing)

        if (associated(gdef)) then
          nx = gdef%nglx
          ny = gdef%ngly
          nz = gdef%nglz
#if defined(petsc)
          mpi_comm = gdef%mpi_comm
        else
          mpi_comm = MPI_COMM_WORLD
#endif
        endif

c     Initialize MPI and PETSc (parallel version)

#if defined(petsc)
        call initFortranMPI(mpi_comm,nx,ny,nz)
        call setup_petsc_BC
#endif

c     Initialize grid structure

        call initGridSt(nx,ny,nz,gdef
#if defined(petsc) || defined(samrai)
     .                 ,mpi_comm=mpi_comm
#endif
#if defined(samrai)
     .                 ,coarsen_in_x=.false.
     .                 ,coarsen_in_y=.false.
     .                 ,coarsen_in_z=.false.
#endif
     .                 )

c     Create MG grid hierarchy

        gdef%g_pack = gpack
        gdef%params = gparams

        call createMGgrid(gdef)

        grid_params => gdef    !For backward compatibility

#if defined(petsc)
c     Create fortran PETSc parallel grid hierarchy (parallel version)

        call createPETScGrid(gdef,0)

c     Create MPI communicators

        call create_MPI_comms
#endif

c     Initialize grid check

        if (check_grid) then

          if (my_rank == 0) then
            write (*,*)
            write (*,*) 'Checking grid...'
          endif

          open(unit=funit,file='grid_p'//trim(int2char(my_rank))//'.txt'
     .        ,status='unknown')

          call checkMGGrid(gdef)
        endif

c     Store grid metric parameters in grid metric structure

        if (checkAnalMapDatabase() .or. coords == 'pck') then
          call defineGridMetric(gdef)
          gmetric => gdef%gmetric !For backward compatibility
        else
          nullify(gmetric)
        endif

      end subroutine createGrid

c     destroyGrid
c     #################################################################
      subroutine destroyGrid(gdef)

c     -----------------------------------------------------------------
c     Deallocates grid memory space
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer :: imin,imax,jmin,jmax,kmin,kmax
        integer :: nx,ny,nz

c     Begin program

c     Destroy metrics storage

        call deallocateGridMetric(gdef%gmetric)

c     Destroy MG grid hierarchy

        call deallocateGridStructure(gdef)

      end subroutine destroyGrid

c     initGridSt
c     #################################################################
      subroutine initGridSt(nxg,nyg,nzg,gdef,mpi_comm
     .                     ,coarsen_in_x,coarsen_in_y,coarsen_in_z)

c     -----------------------------------------------------------------
c     Initializes MG quantities and allocates MG strucutures. On input:
c       * nx,ny,nz: global grid size
c       * gdef (type grid_mg_def): variable that stores grid definition.
c     In parallel, local grid sizes are read from PETSc's DA context.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nxg,nyg,nzg
        type(grid_mg_def),pointer :: gdef
        integer,optional :: mpi_comm
        logical,optional :: coarsen_in_x,coarsen_in_y,coarsen_in_z

c     Local variables

        integer    :: i,n1,n2,n3,nxx,nyy,nzz
        integer    :: ngrdx,ngrdy,ngrdz
        logical    :: cx,cy,cz

c     Begin program

        if (PRESENT(coarsen_in_x)) then
          cx = coarsen_in_x
        else
          cx = .true.
cc          cx = .not.bcSP()
        endif

        if (PRESENT(coarsen_in_y)) then
          cy = coarsen_in_y
        else
          cy = .true.
        endif

        if (PRESENT(coarsen_in_z)) then
          cz = coarsen_in_z
        else
          cz = .true.
        endif

c     Find adequate number of grid levels (for MG)

        if (cx) then
          if (npx == 0)
     .         call pstop('initGridSt','Undefined proc. alloc. in X')

          nxx = nxg/npx
          ngrdx = MG_levels(nxx,mg_ratio)
        else
          ngrdx = 1
        endif

        if (cy) then
          if (npy == 0)
     .         call pstop('initGridSt','Undefined proc. alloc. in Y')

          nyy = nyg/npy
          ngrdy = MG_levels(nyy,mg_ratio)
        else
          ngrdy = 1
        endif

        if (cz) then
          if (npz == 0)
     .         call pstop('initGridSt','Undefined proc. alloc. in Z')

          nzz = nzg/npz
          ngrdz = MG_levels(nzz,mg_ratio)
        else
          ngrdz = 1
        endif

ccc     Estimate local domain sizes
cc
cc        if (npx == 0 .and. npy == 0 .and. npz == 0) then
cc          call pstop('initGridSt','Undefined processor alloc.')
cc        endif
cc
cc        nxx = nxg/npx
cc        nyy = nyg/npy
cc        nzz = nzg/npz
cc
ccc     Find adequate number of grid levels (for MG)
cc
cc        n1 = int(dlog(1d0*nxx)/dlog(1d0*mg_ratio)+0.001)
cc        n2 = int(dlog(1d0*nyy)/dlog(1d0*mg_ratio)+0.001)
cc        n3 = int(dlog(1d0*nzz)/dlog(1d0*mg_ratio)+0.001)
cc
cc        if (cx) then
cc          ngrdx = max(n1-1,1)
cc          do i = ngrdx,1,-1
cc            n1 = nxx/mg_ratio**(i-1)
cc            if (n1*mg_ratio**(i-1).eq.nxx) exit
cc          enddo
cc          ngrdx = i
cc        else
cc          ngrdx = 1
cc        endif
cc
cc        if (cy) then
cc          ngrdy = max(n2-1,1)
cc          do i = ngrdy,1,-1
cc            n2 = nyy/mg_ratio**(i-1)
cc            if (n2*mg_ratio**(i-1).eq.nyy) exit
cc          enddo
cc          ngrdy = i
cc        else
cc          ngrdy = 1
cc        endif
cc
cc        if (cz) then
cc          ngrdz = max(n3-1,1)
cc          do i = ngrdz,1,-1
cc            n3 = nzz/mg_ratio**(i-1)
cc            if (n3*mg_ratio**(i-1).eq.nzz) exit
cc          enddo
cc          ngrdz = i
cc        else
cc          ngrdz = 1
cc        endif

c     Allocate grid storage structure

        call allocateGridStructure(ngrdx,ngrdy,ngrdz
     .                            ,nxg,nyg,nzg,gdef,mpi_comm=mpi_comm)

      contains

c     MG_levels
c     #################################################################
      function MG_levels(nx,mgr) result(ngrd)

c     -----------------------------------------------------------------
c     Determines available MG levels based on fine grid resolution nx
c     and MG ratio mgr
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,mgr,ngrd

c     Local variables

      integer :: nl,i

c     Begin program

      nl = int(dlog(1d0*nx)/dlog(1d0*mgr)+0.001)

      ngrd = max(nl-1,1)

      do i = ngrd,1,-1
        nl = nx/mgr**(i-1)
        if (nl*mgr**(i-1) == nx) exit
      enddo

      ngrd = i

      end function MG_levels

      end subroutine initGridSt

c     createMGgrid
c     #################################################################
      subroutine createMGgrid(gdef)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities from:
c       * gpack (type grid_pack): stores grid packing info
c       * gdef (type grid_mg_def): stores grid definition.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer :: i

c     Begin program

c     Initialize MG grid-size arrays

        !Select refinement ratio (mg_ratio is input from module)
        gdef%mg_ratio_x = 1
        gdef%mg_ratio_x(1:gdef%ngrdx-1) = mg_ratio

        gdef%mg_ratio_y = 1
        gdef%mg_ratio_y(1:gdef%ngrdy-1) = mg_ratio

        gdef%mg_ratio_z = 1
        gdef%mg_ratio_z(1:gdef%ngrdz-1) = mg_ratio

        !Create MG coarse grids
        gdef%nxgl(1) = gdef%nglx
        gdef%nxv (1) = gdef%nlx
        gdef%ilo (1) = gdef%ilog
        gdef%ihi (1) = gdef%ihig
        do i = 2,gdef%ngrid
          gdef%nxgl(i) = gdef%nxgl(i-1)/gdef%mg_ratio_x(i-1)
          gdef%nxv (i) = gdef%nxv (i-1)/gdef%mg_ratio_x(i-1)
          if (gdef%mg_ratio_x(i-1) /= 1) then
            gdef%ilo (i) = gdef%ilo (i-1)/gdef%mg_ratio_x(i-1)+1
            gdef%ihi (i) = gdef%ilo(i)+gdef%nxv(i)-1
          else
            gdef%ilo (i) = gdef%ilo (i-1)
            gdef%ihi (i) = gdef%ihi (i-1)
          endif
        enddo

        gdef%nygl(1) = gdef%ngly
        gdef%nyv (1) = gdef%nly
        gdef%jlo (1) = gdef%jlog
        gdef%jhi (1) = gdef%jhig
        do i = 2,gdef%ngrid
          gdef%nygl(i) = gdef%nygl(i-1)/gdef%mg_ratio_y(i-1)
          gdef%nyv (i) = gdef%nyv (i-1)/gdef%mg_ratio_y(i-1)
          if (gdef%mg_ratio_y(i-1) /= 1) then
            gdef%jlo (i) = gdef%jlo (i-1)/gdef%mg_ratio_y(i-1)+1
            gdef%jhi (i) = gdef%jlo(i)+gdef%nyv(i)-1
          else
            gdef%jlo (i) = gdef%jlo (i-1)
            gdef%jhi (i) = gdef%jhi (i-1)
          endif
        enddo

        gdef%nzgl(1) = gdef%nglz
        gdef%nzv (1) = gdef%nlz
        gdef%klo (1) = gdef%klog
        gdef%khi (1) = gdef%khig
        do i = 2,gdef%ngrid
          gdef%nzgl(i) = gdef%nzgl(i-1)/gdef%mg_ratio_z(i-1)
          gdef%nzv (i) = gdef%nzv (i-1)/gdef%mg_ratio_z(i-1)
          if (gdef%mg_ratio_z(i-1) /= 1) then
            gdef%klo (i) = gdef%klo (i-1)/gdef%mg_ratio_z(i-1)+1
            gdef%khi (i) = gdef%klo(i)+gdef%nzv(i)-1
          else
            gdef%klo (i) = gdef%klo (i-1)
            gdef%khi (i) = gdef%khi (i-1)
          endif
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

        if (gdef%gxmin == UNDEF_DP) gdef%gxmin = xmin
        if (gdef%gxmax == UNDEF_DP) gdef%gxmax = xmax
        if (gdef%gymin == UNDEF_DP) gdef%gymin = ymin
        if (gdef%gymax == UNDEF_DP) gdef%gymax = ymax
        if (gdef%gzmin == UNDEF_DP) gdef%gzmin = zmin
        if (gdef%gzmax == UNDEF_DP) gdef%gzmax = zmax

c     Consistency checks

        call consistencyCheck(gdef)

c     Define uniform logical LOCAL grid on ALL grid levels

        call createLogicalGrid(gdef%nxgl,gdef%mg_ratio_x
     .                     ,gdef%xg,gdef%xx,gdef%dx
     .                     ,gdef%dxh,gdef%nxv
     .                     ,gdef%ngrid,gdef%istartx
     .                     ,gdef%ilo,gdef%ihi
     .                     ,gdef%gxmin,gdef%gxmax,bcond(1),bcond(2)
     .                     ,gdef%g_pack(1))

        call createLogicalGrid(gdef%nygl,gdef%mg_ratio_y
     .                     ,gdef%yg,gdef%yy,gdef%dy
     .                     ,gdef%dyh,gdef%nyv
     .                     ,gdef%ngrid,gdef%istarty
     .                     ,gdef%jlo,gdef%jhi
     .                     ,gdef%gymin,gdef%gymax,bcond(3),bcond(4)
     .                     ,gdef%g_pack(2))

        call createLogicalGrid(gdef%nzgl,gdef%mg_ratio_z
     .                     ,gdef%zg,gdef%zz,gdef%dz
     .                     ,gdef%dzh,gdef%nzv
     .                     ,gdef%ngrid,gdef%istartz
     .                     ,gdef%klo,gdef%khi
     .                     ,gdef%gzmin,gdef%gzmax,bcond(5),bcond(6)
     .                     ,gdef%g_pack(3))

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

        integer :: ngrid,nglobal(ngrid),nx(ngrid),istart(ngrid)
     $            ,mgr(ngrid),ilo(ngrid),ihi(ngrid),bcs1,bcs2
        real(8) :: xg(*),xx(*),dx(*),dxh(*),lmin,lmax

        type(grid_pack) :: dpack

c     Local variables
    
        type(grid_pack) :: lpack
    
        integer :: ig,i,isig,ng,coarse_rat

c     Begin program

        ig = 1  !Finest grid level

        lpack = dpack

        lpack%pack = lpack%pack.and.(coords == 'car')  !Logical packing only for Cartesian

c     Find grid quantities

        do ig = 1,ngrid

          isig = istart(ig)

          ng   = nglobal(ig)

          !Find grid spacing
          call findLocalGrid(ig,xx(isig),ng)

          !Find integer mesh spacings
          dx(isig:isig+nx(ig)+1) = 0d0
          do i = 0,nx(ig)
            dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
          enddo

          !Find half mesh spacings
          dxh(isig:isig+nx(ig)+1) = 0d0
          do i = 1,nx(ig)
            dxh(i + isig) = 0.5*(dx(i + isig) + dx(i-1 + isig))
          enddo

          if (bcs1 == PER) then
            dxh(0        + isig) = dxh(nx(ig) + isig)
            dxh(nx(ig)+1 + isig) = dxh(1      + isig)
          else
cc            dxh(0        + isig) = dx(0      + isig)*0.5 !Needed for extrapolation
cc            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)*0.5 !Needed for extrapolation
            dxh(0        + isig) = dx(0      + isig)
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)
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

        integer :: ig,ng
        real(8) :: xx(*)

c     Local variables

        real(8) :: xxg(0:ng+1)

c     Begin program

        if (ig == 1) then
          coarse_rat = 1
        else
          coarse_rat = coarse_rat*mgr(ig-1)
        endif

c     Find global grid

        call findGlobalGrid(lmin,lmax,bcs1,bcs2,xxg,lpack,coarse_rat)

cc        write (*,*) 'Coarsening ratios=',mgr
cc        write (*,*) 'Grid level=',ig
cc        write (*,*) 'Domain limits=',lmin,lmax
cc        write (*,*) xxg
cc        write (*,*)
cc        pause

c     Assign global grid

        if (ig == 1) xg(1:ng+2) = xxg(0:ng+1)

c     Assign local grid

        xx(1:nx(ig)+2) = xxg(ilo(ig)-1:ihi(ig)+1)

c     End program

      end subroutine findLocalGrid

      end subroutine createLogicalGrid

c     consistencyCheck
c     #################################################################
      subroutine consistencyCheck(gdef)

c     -----------------------------------------------------------------
c     Checks consistency of grid parameters. On input, MG grid structure.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer :: nxx,nyy,nzz
        real(8) :: major_r,xmin,xmax,ymin,ymax,zmin,zmax

c     Begin program

        nxx = gdef%nglx
        nyy = gdef%ngly
        nzz = gdef%nglz

        xmin = gdef%gxmin
        xmax = gdef%gxmax
                         
        ymin = gdef%gymin
        ymax = gdef%gymax
                         
        zmin = gdef%gzmin
        zmax = gdef%gzmax

        !Consistency
        numerical_grid = .true.  !Needed because of vertex metrics

        select case (coords)
        case ('car')
cc          numerical_grid = .false.
        case ('scl')
          numerical_grid = .true.   !To deal properly with ghost cells

          gdef%params = 0d0

          !Save packing configuration
          if (gdef%g_pack(1)%pack) then
            gdef%params(1) = gdef%g_pack(1)%xp
            gdef%params(2) = gdef%g_pack(1)%dx0/(xmax-xmin)*nxx
          endif

          if (gdef%g_pack(2)%pack) then
            gdef%params(3) = gdef%g_pack(2)%xp
            gdef%params(4) = gdef%g_pack(2)%dx0/(ymax-ymin)*nyy
          endif

          if (gdef%g_pack(3)%pack) then
            gdef%params(5) = gdef%g_pack(3)%xp
            gdef%params(6) = gdef%g_pack(3)%dx0/(zmax-zmin)*nzz
          endif
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

          major_r = gdef%params(1)

          if (major_r < xmax) then
            write (*,*) 'Ill-defined toroidal coordinate system'
            write (*,*) 'Major radius < minor radius'
            write (*,*) 'Aborting'
            stop
          endif
        end select

        !Ensure ignorable directions are small for numerical computation
        !  of grid parameters
        if (numerical_grid .and. my_rank == 0) then
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

c     checkMGGrid
c     #################################################################
      subroutine checkMGGrid(gdef)

c     -----------------------------------------------------------------
c     Dumps grid info and checks metric tensors.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer    :: i,j,n1,n2,n3
        integer    :: igx,igy,igz,isig

c     Begin program

c     Grid check output

        write (funit,*)
        write (funit,*) 'Coordinate system: ',coords
        write (funit,*)
        write (funit,*) 'Number of grid levels:'
        write (funit,*) 'nx: ',gdef%ngrid
     .              ,'   ny: ',gdef%ngrid
     .              ,'   nz: ',gdef%ngrid

c     Multigrid parameters

        igx = gdef%ngrid
        igy = gdef%ngrid
        igz = gdef%ngrid

        call gridInfo('X',igx,gdef%istartx,gdef%nxv
     .               ,gdef%xx,gdef%dx,gdef%dxh)

        call gridInfo('Y',igy,gdef%istarty,gdef%nyv
     .               ,gdef%yy,gdef%dy,gdef%dyh)

        call gridInfo('Z',igz,gdef%istartz,gdef%nzv
     .               ,gdef%zz,gdef%dz,gdef%dzh)

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

        write (funit,*)
        write (funit,*) '***************************'
        write (funit,*) 'MG grid in ',char,'-axis'
        write (funit,*) '***************************'
        do i = 1,ig
          isig = istart(i)
          write (funit,*)
          write (funit,*) '************* Grid ',i,' **************'
          write (funit,*)
          write (funit,*) 'Size ',nxv(i)
          write (funit,*) 'MG pointer: ',isig
          write (funit,*) 'Grid nodes'
          do j = isig,isig+nxv(i)+1
            write (funit,10) 'node :',j-isig
     .                  ,'   position: ',xx(j)
     .                  ,'   int dh: ',dx(j)
     .                  ,'   half dh: ',dxh(j)
 10         format (a,i3,a,f10.3,a,f10.3,a,f10.3)
          enddo
        enddo

c     End program

      end subroutine gridInfo

      end subroutine checkMGGrid

      end module grid_create

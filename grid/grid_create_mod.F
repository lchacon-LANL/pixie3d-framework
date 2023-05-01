c module grid_create
c #####################################################################
      module grid_create

        use grid_metric_init

        implicit none

        integer :: mg_ratio=2  !Default MG coarsening ratio

        logical :: coarsen_X=.true.
     .            ,coarsen_Y=.true.
     .            ,coarsen_Z=.true.

        private :: createLogicalGrid
     .            ,createMGgrid
     $            ,consistencyCheck!,initGridSt,checkGrid

#if !defined(petsc)
        integer :: npx=1,npy=1,npz=1  !Serial processor allocation
#else
        integer :: npx=0,npy=0,npz=0  !Parallel processor allocation
#endif
        logical,private :: map_grid_pack=.false.

        real(8),parameter :: LMIN_NEQ_ONE=1d-3
      contains

c     createGrid
c     #################################################################
      subroutine createGrid(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,gdef)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities for global mesh
c     of size nx x ny x nz. Results are stored in grid structure gdef.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: nx,ny,nz
        type(grid_mg_def),pointer :: gdef
        real(8) :: xmin,xmax,ymin,ymax,zmin,zmax

c     Local variables

        integer :: imin,imax,jmin,jmax,kmin,kmax
        integer :: mpi_comm

c     Begin program

c     Initialize global quantities (from gdef, if existing) and MPI

        if (associated(gdef)) then
          nx = gdef%nglx
          ny = gdef%ngly
          nz = gdef%nglz
#if defined(petsc)
          mpi_comm = gdef%mpi_comm
          call processorAlloc(mpi_comm,nx,ny,nz,gdef%mx,gdef%my,gdef%mz)
#endif
        else
#if defined(petsc)
          mpi_comm = MPI_COMM_WORLD
          call processorAlloc(mpi_comm,nx,ny,nz,npx,npy,npz)
#endif
        endif

c     Create grid structures

        call create_gdef(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,gdef
     .                  ,mpi_comm)

c     Legacy assignments (backward compatibility)

#if defined(petsc)
        call initMPI(mpi_comm,np,my_rank)
#endif

      end subroutine createGrid

c     destroyGrid
c     #################################################################
      recursive subroutine destroyGrid(gdef)

c     -----------------------------------------------------------------
c     Deallocates grid memory space
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gdef

c     Local variables

c     Begin program        

        if (.not.associated(gdef)) return

c     Deallocate coarser levels

        call destroyGrid(gdef%g_crse_def)

c     Deallocate MPI COMMs
        
#if defined(petsc)
        call destroy_MPI_comms(gdef)
#endif

c     Deallocate grid structure
        
        call deallocateGridStructure(gdef)

      end subroutine destroyGrid

c     create_gdef
c     ######################################################################
      subroutine create_gdef(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax
     .                    ,g_def,mpi_comm,coarsen_x,coarsen_y,coarsen_z
     .                    ,np_x,np_y,np_z,pack_x,pack_y,pack_z)
 
      implicit none

c     ----------------------------------------------------------------------
c     Initialize PARALLEL grid structure
c     ----------------------------------------------------------------------
     
c     Call variables

      integer :: nx,ny,nz,mpi_comm
      type(grid_mg_def),pointer,INTENT(IN)  :: g_def
      real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
      logical,optional :: coarsen_x,coarsen_y,coarsen_z
     .                   ,pack_x,pack_y,pack_z
      integer,optional :: np_x,np_y,np_z

c     Local variables

      integer :: il

      logical :: exit_loop

      real(8),pointer,dimension(:,:,:,:) :: carg

      type(grid_mg_def),pointer :: g_def_sav => null()

      type(grid_pack),dimension(3) :: g_pack

      logical :: crsn_x,crsn_y,crsn_z

c     Begin program
      
c     Initialize grid structure

      call get_grid_pack_config(g_pack)

      if (PRESENT(pack_x)) g_pack(1)%pack = pack_x
      if (PRESENT(pack_y)) g_pack(2)%pack = pack_y
      if (PRESENT(pack_z)) g_pack(3)%pack = pack_z

      if (nx == 1) g_pack(1)%pack = .false.
      if (ny == 1) g_pack(2)%pack = .false.
      if (nz == 1) g_pack(3)%pack = .false.

      g_pack(1)%dxl = g_pack(1)%dx0
      g_pack(1)%dxr = g_pack(1)%dx0
      g_pack(2)%dxl = g_pack(2)%dx0
      g_pack(2)%dxr = g_pack(2)%dx0
      g_pack(3)%dxl = g_pack(3)%dx0
      g_pack(3)%dxr = g_pack(3)%dx0
      
      call init_gdef(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax
     .               ,g_pack,gparams,g_def
#if defined(petsc) || defined(samrai)
     .               ,mpi_comm=mpi_comm
     .               ,np_x=np_x
     .               ,np_y=np_y
     .               ,np_z=np_z
#endif
#if defined(samrai)
     .               ,coarsen_in_x=.false.
     .               ,coarsen_in_y=.false.
     .               ,coarsen_in_z=.false.
#else
     .               ,coarsen_in_x=coarsen_x
     .               ,coarsen_in_y=coarsen_y
     .               ,coarsen_in_z=coarsen_z
     .               )
#endif

c     Store grid metric parameters in grid metric structure

      if (checkAnalMapDatabase()) then
        call defineGridMetric(g_def)
      else
        if (my_rank == 0 .and. coords /= 'ext') then
          write (*,*) 'WARNING: Map "',coords,'" not in database'
        endif
      endif

c     Create beyond-processor sub-grids for coarse MG solve

#if defined(petsc) && defined(coarse_MG) && !defined(samrai)
      g_def%mpi_comm_crse = MPI_COMM_WORLD

      il = 0
      g_def_sav => g_def
      do
        il = il + 1

        exit_loop= create_coarse_gdef(il,g_def_sav,g_def_sav%g_crse_def
     .                               ,cx=coarsen_x
     .                               ,cy=coarsen_y
     .                               ,cz=coarsen_z
     .                               )
        if (exit_loop) exit

        g_def_sav => g_def_sav%g_crse_def
      enddo

      if (coarse_mg_debug) call coarse_mg_debug_io
#endif

c     Perform grid check

      if (check_grid.and.checkAnalMapDatabase()) call checkGrid(g_def)

c     End program

#if defined(petsc) && defined(coarse_MG)
      contains

c     coarse_mg_debug_io
c     ######################################################################
      subroutine coarse_mg_debug_io

      implicit none

      integer :: il

      if (my_rank == 0) then
        write (*,*) '-------------'
        write (*,*) 'Final summary'
        write (*,*) '-------------'
      endif

      il = 0

      g_def_sav => g_def

      if (my_rank == 0) write (*,*) 'Fine levels',g_def_sav%ngrid
      if (my_rank == 0) write (*,*) "Fine X dims",g_def_sav%nxv
      if (my_rank == 0) write (*,*) "Fine Y dims",g_def_sav%nyv
      if (my_rank == 0) write (*,*) "Fine Z dims",g_def_sav%nzv

      do while(associated(g_def_sav%g_crse_def))
        il = il + 1

        if (my_rank == 0) then
          write (*,*) "************************"
          write (*,*) "Level "//trim(int2char(il))//" associated"
          write (*,*) "************************"
        endif

        g_def_sav => g_def_sav%g_crse_def

        if (my_rank == 0) write (*,*) 'Coarse levels',g_def_sav%ngrid
        if (my_rank == 0) write (*,*) "Coarse X dims",g_def_sav%nxv
        if (my_rank == 0) write (*,*) "Coarse Y dims",g_def_sav%nyv
        if (my_rank == 0) write (*,*) "Coarse Z dims",g_def_sav%nzv

        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      enddo

      end subroutine coarse_mg_debug_io
#endif
        
      end subroutine create_gdef

#if defined(petsc) && defined(coarse_MG)
c     create_coarse_gdef
c     ######################################################################
      function create_coarse_gdef(il,g_def,g_crse_def,cx,cy,cz)
     .         result(exit_loop)
 
      implicit none

c     ----------------------------------------------------------------------
c     Initialize grid structure over COARSE mesh
c     ----------------------------------------------------------------------
     
c     Call variables

      integer :: il
      logical :: exit_loop
      type(grid_mg_def),pointer,INTENT(INOUT)  :: g_def
      type(grid_mg_def),pointer,INTENT(OUT) :: g_crse_def
      logical,optional :: cx,cy,cz

c     Local variables

      integer :: mpi_comm_splt,mpi_comm_crse,mgcf
      integer :: igr,nxg,nyg,nzg,npx,npy,npz,ierr

c     Begin program

      exit_loop = create_MG_comms(g_def,mpi_comm_splt,mpi_comm_crse)

      igr = g_def%ngrid  !Last grid level

c     Compute # procs involved

      mgcf = MG_get_proc_coarsen_factor()

      npx = max(g_def%mx/mgcf,1)
      npy = max(g_def%my/mgcf,1)
      npz = max(g_def%mz/mgcf,1)

c     Global mesh

      nxg = g_def%nxgl(igr)
      nyg = g_def%nygl(igr)
      nzg = g_def%nzgl(igr)

c     Init coarse mesh structure

      call init_gdef(nxg,nyg,nzg
     .              ,g_def%gxmin,g_def%gxmax
     .              ,g_def%gymin,g_def%gymax
     .              ,g_def%gzmin,g_def%gzmax
     .              ,g_def%g_pack,g_def%params
     .              ,g_crse_def
     .              ,np_x=npx,np_y=npy,np_z=npz
     .              ,mpi_comm=mpi_comm_splt
     .              ,coarsen_in_x=cx
     .              ,coarsen_in_y=cy
     .              ,coarsen_in_z=cz
     .              )

      g_crse_def%mpi_comm_crse = mpi_comm_crse

c     Gather map and define grid metrics

      if (checkAnalMapDatabase()) then
        call defineGridMetric(g_crse_def,ierr=ierr)
        if (g_crse_def%ngrid < 1.and.ierr /= 0) then !Could not evaluate metrics at top level
          exit_loop = .true.
          call destroyGrid(g_crse_def)
        endif
      endif

c     Diagnostics

      if (coarse_mg_debug) then
cc        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
        if (my_rank == 0) then
          write (*,*) "************************"
          write (*,*) "Level "//trim(int2char(il))//" created"
          write (*,*) "************************"
        endif
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
      endif

c     End program
        
      end function create_coarse_gdef
#endif

c     createMGMetricHierarchy
c     #################################################################
      subroutine createMGMetricHierarchy(g_def,xf)

        implicit none

c     ----------------------------------------------------------------------
c     Creates grid metric hierarchy in g_def when providing external map xcar
c     ----------------------------------------------------------------------

c     Call variables

        type(grid_mg_def),pointer,INTENT(INOUT)  :: g_def

        real(8),dimension(:,:,:,:),optional :: xf

c     Local variables

        integer :: nxl,nyl,nzl,ierr,dim
        real(8),pointer,dimension(:,:,:,:) :: xc
        type(grid_mg_def),pointer :: gc_def => null()
     .                              ,gf_def => null()

c     Begin program

        if (.not.associated(g_def)) then
           call pstop("createGridMetricHierarchy"
     $               ,"Target Metric not allocated")
        endif
        
        if (PRESENT(xf)) then
          call defineGridMetric(g_def,xcar=xf,ierr=ierr)

          if (ierr /= 0 .and. my_rank == 0) then
            write (*,*) ' >>> Discarding last mesh due to jac < 0'
            return
          endif
        endif

c     Cycle through grid hierarchies

        gf_def => g_def
        gc_def => g_def%g_crse_def

        do while (associated(gc_def))

          if (my_rank == 0)
     .         write (*,*) '   >>> New coarse grid hierarchy'

          nxl = gc_def%nxv(1)
          nyl = gc_def%nyv(1)
          nzl = gc_def%nzv(1)
          allocate(xc(0:nxl+1,0:nyl+1,0:nzl+1,3))

          !Gather coarsest mesh over sub-procs
          call find_global(gf_def%gmetric%grid(gf_def%ngrid)%car,xc
#if defined(petsc) 
     .                    ,mpi_comm=gc_def%mpi_comm_crse
#endif
     .                    )

          call defineGridMetric(gc_def,xcar=xc,ierr=ierr)

          deallocate(xc)

          !Check for negative Jac, and if so deallocate subsequent grid levels and exit
          if (ierr /= 0) then
            if (my_rank == 0) then
              write (*,*) ' >>> defineGridMetric error',ierr
              write (*,*) ' >>> Discarding leftover meshes due to jac<0'
            endif
            call destroyGrid(gf_def%g_crse_def)
            exit                
          endif

          gf_def => gc_def
          gc_def => gf_def%g_crse_def

        enddo

        nullify(gf_def,gc_def)

c     End program

      end subroutine createMGMetricHierarchy

c     init_gdef
c     #################################################################
      subroutine init_gdef(nxg,nyg,nzg,gxmin,gxmax,gymin,gymax,gzmin
     .     ,gzmax,g_pack,gparams,gdef,mpi_comm,np_x,np_y,np_z
     .     ,coarsen_in_x,coarsen_in_y,coarsen_in_z)

c     -----------------------------------------------------------------
c     Initializes grid structure, including MG pointers and PETSc
c     hierarchical mesh. On input:
c       * nxg,nyg,nzg: global grid size
c       * gxmin,gxmax,gymin,gymax,gzmin,gzmax: global grid size
c       * gdef (type grid_mg_def): variable that stores grid definition.
c       * g_pack: grid packing info
c       * gparams: grid setup info
c     Optional parameters:
c       * mpi_comm: MPI communicator
c       * np_x,np_y,np_z: processor distribution
c       * coarsen_in_x,coarsen_in_y,coarsen_in_z: whether to coarsen
c         mesh in a given axis
c     In parallel, local grid sizes are read from PETSc's DA context.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: nxg,nyg,nzg

        type(grid_mg_def),pointer :: gdef

        real(8) :: gxmin,gxmax,gymin,gymax,gzmin,gzmax

        type(grid_pack),dimension(3) :: g_pack

        real(8) :: gparams(ngparams)

        integer,optional :: mpi_comm,np_x,np_y,np_z
        logical,optional :: coarsen_in_x,coarsen_in_y,coarsen_in_z

c     Local variables

        integer :: i,n1,n2,n3,nxx,nyy,nzz,mx,my,mz
        integer :: ngrdx,ngrdy,ngrdz
        logical :: cx,cy,cz

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

        if (PRESENT(np_x)) then
          mx = np_x
        else
          mx = npx
        endif

        if (PRESENT(np_y)) then
          my = np_y
        else
          my = npy
        endif

        if (PRESENT(np_z)) then
          mz = np_z
        else
          mz = npz
        endif

c     Find adequate number of grid levels (for MG)

        if (cx) then
          if (mx == 0)
     .         call pstop('init_gdef','Undefined proc. alloc. in X')

          if (mx > nxg) call pstop('init_gdef','Too many procs in X')

          nxx = nxg/mx
          ngrdx = MG_levels(nxx,mg_ratio,mx)
        else
          ngrdx = 1
        endif

        if (cy) then
          if (my == 0)
     .         call pstop('init_gdef','Undefined proc. alloc. in Y')

          if (my > nyg) call pstop('init_gdef','Too many procs in Y')

          nyy = nyg/my
          ngrdy = MG_levels(nyy,mg_ratio,my)
        else
          ngrdy = 1
        endif

        if (cz) then
          if (mz == 0)
     .         call pstop('init_gdef','Undefined proc. alloc. in Z')

          if (mz > nzg) call pstop('init_gdef','Too many procs in Z')

          nzz = nzg/mz
          ngrdz = MG_levels(nzz,mg_ratio,mz)
        else
          ngrdz = 1
        endif

c     Allocate grid storage structure

        call allocateGridStructure(ngrdx,ngrdy,ngrdz,mx,my,mz
     .                            ,nxg,nyg,nzg,gdef,mpi_comm=mpi_comm)

c     Adjust domain size

        !Default sizes
        if (gxmin == 0d0 .and. gxmax == 0d0) then
          gxmax = 2*pi
          gxmin = 0d0
        endif

        if (gzmin == 0d0 .and. gymax == 0d0) then
          gymax = 2*pi
          gymin = 0d0
        endif

        if (gzmin == 0d0 .and. gzmax == 0d0) then
          gzmax = 2*pi
          gzmin = 0d0
        endif

c$$$        if (mx == 1 .and. (.not.cx)) then
c$$$          gdef%gxmin = xmin
c$$$          gdef%gxmax = xmax/npx
c$$$        endif
c$$$
c$$$        if (my == 1 .and. (.not.cy)) then
c$$$          gdef%gymin = ymin
c$$$          gdef%gymax = ymax/npy
c$$$        endif
c$$$
c$$$        if (mz == 1 .and. (.not.cz)) then
c$$$          gdef%gzmin = zmin
c$$$          gdef%gzmax = zmax/npz
c$$$        endif

        if (gdef%gxmin == UNDEF_DP) gdef%gxmin = gxmin
        if (gdef%gxmax == UNDEF_DP) gdef%gxmax = gxmax
        if (gdef%gymin == UNDEF_DP) gdef%gymin = gymin
        if (gdef%gymax == UNDEF_DP) gdef%gymax = gymax
        if (gdef%gzmin == UNDEF_DP) gdef%gzmin = gzmin
        if (gdef%gzmax == UNDEF_DP) gdef%gzmax = gzmax

c     Create MG grid hierarchy

        gdef%g_pack = g_pack
        gdef%params = gparams

        map_grid_pack = setup_grid_packing(coords,gdef)

        call createMGgrid(gdef)

c     Create fortran PETSc parallel grid hierarchy (parallel version)

#if defined(petsc)
        call createPETScGrid(gdef,0)
        call setup_petsc_BC(gdef)   !Set PETSc BC flags
        call create_MPI_comms(gdef) !Set MPI comms
#endif

      end subroutine init_gdef

c     createMGgrid
c     #################################################################
      subroutine createMGgrid(gdef)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities in MG grid hierarchy.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

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

c     Consistency checks

        call consistencyCheck(gdef)

c     Define logical LOCAL grid on ALL grid levels

        call createLogicalGrid(gdef%nxgl,gdef%mg_ratio_x
     .                        ,gdef%xg,gdef%xx,gdef%dx
     .                        ,gdef%dxh,gdef%nxv
     .                        ,gdef%ngrid,gdef%istartx
     .                        ,gdef%ilo,gdef%ihi
     .                        ,gdef%gxmin,gdef%gxmax
     .                        ,bcond(1),bcond(2)
     .                        ,gdef%g_pack(1))

        gdef%lxmin = gdef%xx(1)
        gdef%lxmax = gdef%xx(gdef%nxv(1)+2)

        call createLogicalGrid(gdef%nygl,gdef%mg_ratio_y
     .                     ,gdef%yg,gdef%yy,gdef%dy
     .                     ,gdef%dyh,gdef%nyv
     .                     ,gdef%ngrid,gdef%istarty
     .                     ,gdef%jlo,gdef%jhi
     .                     ,gdef%gymin,gdef%gymax,bcond(3),bcond(4)
     .                     ,gdef%g_pack(2))

        gdef%lymin = gdef%yy(1)
        gdef%lymax = gdef%yy(gdef%nyv(1)+2)

        call createLogicalGrid(gdef%nzgl,gdef%mg_ratio_z
     .                     ,gdef%zg,gdef%zz,gdef%dz
     .                     ,gdef%dzh,gdef%nzv
     .                     ,gdef%ngrid,gdef%istartz
     .                     ,gdef%klo,gdef%khi
     .                     ,gdef%gzmin,gdef%gzmax,bcond(5),bcond(6)
     .                     ,gdef%g_pack(3))

        gdef%lzmin = gdef%zz(1)
        gdef%lzmax = gdef%zz(gdef%nzv(1)+2)

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
    
        integer :: ig,i,isig,ng

        real(8),allocatable,dimension(:) :: xxg
        
        real(8) :: xi_l,xi_r,xi_0

        real(8) :: nxi_l,nxi_r,nxi_0

c     Begin program

c     Find grid quantities

        do ig = 1,ngrid

          isig = istart(ig)

          ng   = nglobal(ig)

          !Find global mesh
          allocate(xxg(0:ng+1))
          call findGlobalGrid(lmin,lmax,bcs1,bcs2,dpack,xxg
     .                       ,xi_l=xi_l,xi_r=xi_r,xi_0=xi_0)
          
c$$$          if (my_rank == 0) then
c$$$            write (*,*) 'Grid level=',ig
c$$$            write (*,*) 'Domain limits=',lmin,lmax
c$$$            write (*,*) 'Coarsening ratio  =',mgr(ig)
c$$$            write (*,*) 'Packing spacings  =',dpack%dxl
c$$$     .                                       ,dpack%dx0
c$$$     .                                       ,dpack%dxr
c$$$            write (*,*) 'Geometric  factors=',xi_l
c$$$     .                                       ,xi_0
c$$$     .                                       ,xi_r
c$$$c$$$            if (dpack%unif_nx == 0) then
c$$$c$$$              write (*,*) 'Coarsening factors=',nxi_l
c$$$c$$$     .                                         ,nxi_0
c$$$c$$$     .                                         ,nxi_r
c$$$c$$$            else
c$$$c$$$              write (*,*) 'Coarser uniform nx=',dpack%unif_nx
c$$$c$$$            endif
c$$$            do i=0,ng+1
c$$$              write (*,*) 'Mesh',i,xxg(i)
c$$$            enddo
c$$$            write (*,*)
c$$$          endif
c$$$          pause

          !Compute new packing spacings for current grid and store in
          !dpack for MG coarsening
          if (ig /= ngrid) then

            if (dpack%unif_nx == 0) then
              nxi_l = 1d0 ; nxi_r = 1d0 ; nxi_0 = 1d0
              do i=2,mgr(ig)
                nxi_l = nxi_l + 2*(mgr(ig)-i+1)/mgr(ig)*xi_l**(i-1)
                nxi_r = nxi_r + 2*(mgr(ig)-i+1)/mgr(ig)*xi_r**(i-1)
                nxi_0 = nxi_0 + 2*(mgr(ig)-i+1)/mgr(ig)*xi_0**(i-1)
              enddo

              dpack%dx0 = dpack%dx0*nxi_0
              dpack%dxl = dpack%dxl*nxi_l
              dpack%dxr = dpack%dxr*nxi_r
            else
              dpack%unif_nx = dpack%unif_nx/mgr(ig)
              dpack%left_nx = dpack%left_nx/mgr(ig)
              dpack%rite_nx = dpack%rite_nx/mgr(ig)
            endif
          endif

          !Assign top-level global grid
          if (ig == 1) xg(1:ng+2) = xxg(0:ng+1)

          !Assign local grid at every level
          xx(isig:isig+nx(ig)+1) = xxg(ilo(ig)-1:ihi(ig)+1)

          !Find integer mesh spacings
c$$$          dx(isig:isig+nx(ig)+1) = 0d0
c$$$          do i = 0,nx(ig)
c$$$            dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
c$$$          enddo
c$$$          dx(nx(ig)+1+isig) = dx(nx(ig)+isig)
c$$$     .                      *(dx(nx(ig)+isig)/dx(nx(ig)-1+isig))

          dx(isig:isig+nx(ig)+1) = 0d0
          do i = 0,nx(ig)+1
            if (ilo(ig)+i <= ng+1) then
              dx(i+isig) = xxg(ilo(ig)+i) - xxg(ilo(ig)+i-1)
            else
              dx(i+isig) = dx(i-1+isig)
            endif
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
            dxh(0        + isig) = dx (0      + isig)
            dxh(nx(ig)+1 + isig) = dx (nx(ig) + isig)
          endif

          deallocate(xxg)
          
        enddo

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

c     Coordinate-specific checks

cc        if (map_grid_pack) numerical_grid = .true.

        select case (coords)
        case ('s3d','hsn')
          numerical_grid = .true.
        case ('car')
          call set_map(map_grid_pack)
          if (.not.have_map()) numerical_grid = .false.
        case ('cyl','hel','hl2')
          if (xmin /= 0d0 .and. bcSP()) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif
        case ('tor','sha','sh2','sh3','stl','tsq')
          if (xmin /= 0d0 .and. bcSP()) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif

          if (    coords=='sh2'
     .        .or.coords=='sh3'
     .        .or.coords=='stl'
     .        .or.coords=='tsq') numerical_grid = .true.

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
        if (numerical_grid .and. coords /= 'car') then
          if    ((nxx == 1 .and. ((xmax-xmin)>LMIN_NEQ_ONE))
     .       .or.(nyy == 1 .and. ((ymax-ymin)>LMIN_NEQ_ONE))
     .       .or.(nzz == 1 .and. ((zmax-zmin)>LMIN_NEQ_ONE))) then
            write (*,*) nxx,nyy,nzz
            write (*,*) xmax,xmin,ymax,ymin,zmax,zmin
            call pstop('consistencyCheck'
     .                ,'Error: large grid spacing in'//
     .                 ' ignorable direction')
          endif
        endif

      end subroutine consistencyCheck

c     checkGrid
c     #################################################################
      subroutine checkGrid(g_def)

c     -----------------------------------------------------------------
c     Dumps grid info and checks metric tensors.
c     -----------------------------------------------------------------

        use io

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: funit,il

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin

        character(100) :: file
        character(4) :: level

        type(grid_mg_def),pointer :: gdef => null()

c     Begin program

        xmin = g_def%gxmin
        xmax = g_def%gxmax

        ymin = g_def%gymin
        ymax = g_def%gymax

        zmin = g_def%gzmin
        zmax = g_def%gzmax

        if (my_rank == 0) then
          write (*,*)
          write (*,*) 'Checking grid...'
        endif

        file = 'grid_p'//trim(int2char(my_rank))//'.txt'

        funit = find_unit(11111)
        open(unit=funit,file=trim(file),status='unknown')

        il = 0
        gdef => g_def
        do while(associated(gdef))
          level = 'cl'//trim(int2char(il))

          write (funit,*)
          write (funit,*)'*************************************'
          write (funit,*)'Grid check in COARSE grid level '
     .                  //trim(int2char(il))
          write (funit,*)'*************************************'
          write (funit,*)
          write (funit,*)

          if (my_rank == 0) then
            write (*,*)
            write (*,*)'*************************************'
            write (*,*)'Grid check in COARSE grid level '
     .             //trim(int2char(il))
            write (*,*)'*************************************'
            write (*,*)
            write (*,*)
          endif
          
          call MGcheck

          call numericalGridCoeffsCheck

          call graphicsCheck_cc

          call graphicsCheck_vrtx

          call hessianZeroForceCheck(1)

cc          call symmCheck(1)

          il = il + 1

          gdef => gdef%g_crse_def
        end do

        close(funit)

        call pstop("checkGrid","Finished grid check")

c     End program

      contains

c     MGcheck
c     #################################################################
      subroutine MGcheck

c     -----------------------------------------------------------------
c     Dumps grid info and checks metric tensors.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

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

        call MG_grid_info('X',igx,gdef%istartx,gdef%nxv
     .               ,gdef%xx,gdef%dx,gdef%dxh)

        call MG_grid_info('Y',igy,gdef%istarty,gdef%nyv
     .               ,gdef%yy,gdef%dy,gdef%dyh)

        call MG_grid_info('Z',igz,gdef%istartz,gdef%nzv
     .               ,gdef%zz,gdef%dz,gdef%dzh)

c     End program

      end subroutine MGcheck

c     MG_grid_info
c     #################################################################  
      subroutine MG_grid_info(char,ig,istart,nxv,xx,dx,dxh)

        implicit none

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

      end subroutine MG_grid_info

c$$$c     hessianCheck
c$$$c     #################################################################
c$$$      subroutine hessianCheck(igx,igy,igz,nxv,nyv,nzv,gdef)
c$$$
c$$$        implicit  none
c$$$
c$$$c     Call variables
c$$$
c$$$        integer    :: igx,igy,igz,nxv(*),nyv(*),nzv(*)
c$$$
c$$$        type(grid_mg_def),pointer :: gdef
c$$$
c$$$c     Local variables
c$$$
c$$$        integer    :: i,j,k,i1,j1,k1,ig,jg,kg,ih
c$$$        real(8)    :: x1,y1,z1
c$$$        real(8)    :: hess(3,3,3),hess_cnv(3,3,3),table(3,3,3)
c$$$
c$$$c     Begin program
c$$$
c$$$        write (funit,*) 
c$$$        write (funit,*) '*************************'
c$$$        write (funit,*) '     Hessian check'
c$$$        write (funit,*) '*************************'
c$$$        write (funit,*)
c$$$
c$$$        do i = 2,nxv(igx)   !Start at 2 to avoid singular points
c$$$          do j = 1,nyv(igy)
c$$$            do k = 1,nzv(igz)
c$$$
c$$$              call getMGmap(gdef,i,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$              x1 = gdef%xx(ig)
c$$$              y1 = gdef%yy(jg)
c$$$              z1 = gdef%zz(kg)
c$$$
c$$$              hess = -christ_2knd(gdef,i,j,k,igx,igy,igz)
c$$$
c$$$              hess_cnv(1,:,:) = hessian_cnv(gdef,1,x1,y1,z1)
c$$$              hess_cnv(2,:,:) = hessian_cnv(gdef,2,x1,y1,z1)
c$$$              hess_cnv(3,:,:) = hessian_cnv(gdef,3,x1,y1,z1)
c$$$
c$$$              do i1=1,3
c$$$                do j1=1,3
c$$$                  do k1=1,3
c$$$                    table(i1,j1,k1) = hess_cnv(i1,j1,k1)
c$$$     .                               -delta(i1,k1)*(hess(1,j1,1)
c$$$     .                                             +hess(2,j1,2)
c$$$     .                                             +hess(3,j1,3))
c$$$     .                               +hess(k1,i1,j1)
c$$$                  enddo
c$$$                enddo
c$$$              enddo
c$$$
c$$$              write (funit,5) 'Grid point: (',x1,',',y1,',',z1,')'
c$$$ 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)
c$$$
c$$$              write (funit,*) 
c$$$              write (funit,*) 'Hessian relation'
c$$$
c$$$              do ih=1,3
c$$$                write (funit,*) 
c$$$                write (funit,10) table(ih,1,1:3)
c$$$                write (funit,10) table(ih,2,1:3)
c$$$                write (funit,10) table(ih,3,1:3)
c$$$ 10             format (3f10.3)
c$$$              enddo
c$$$
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine hessianCheck

c     metricTensorCheck
c     #################################################################
      subroutine metricTensorCheck(igx,igy,igz,nxv,nyv,nzv,gdef)

        implicit  none

c     Call variables

        integer    :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

        type(grid_mg_def),pointer :: gdef

c     Local variables

        integer    :: i,j,k,i1,j1,k1,ig,jg,kg,ih
        real(8)    :: x1,y1,z1,check
        real(8)    :: gup(3,3),gdown(3,3),tensor(3,3)
        logical    :: cartesian

c     Begin program

        check = 0d0

        write (funit,*) 
        write (funit,*) '*************************'
        write (funit,*) '   Metric Tensor check'
        write (funit,*) '*************************'
        write (funit,*)
        write (funit,*) 'Grid levels:',igx,igy,igz
        write (funit,*)

        do k = 1,nzv(igz)
          do j = 1,nyv(igy)
            do i = 1,nxv(igx)   !Start at 2 to avoid singular points

              call getMGmap(gdef,i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = gdef%xx(ig)
              y1 = gdef%yy(jg)
              z1 = gdef%zz(kg)

              gup   = gdef%gmetric%grid(igx)%gsup(i,j,k,:,:)
              gdown = gdef%gmetric%grid(igx)%gsub(i,j,k,:,:)

              tensor = matmul(gup,gdown)

cc              write (funit,5) 'Grid point: (',x1,',',y1,',',z1,')'
cc 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)
              write (funit,5) 'Grid point: (',i,',',j,',',k,')'
 5            format (/,a,i3,a,i3,a,i3,a)

              write(funit,*)
              write(funit,*) 'Gup'
              write(funit,10) gup(1,1:3)
              write(funit,10) gup(2,1:3)
              write(funit,10) gup(3,1:3)

              write(funit,*)
              write(funit,*) 'Gdown'
              write(funit,10) gdown(1,1:3)
              write(funit,10) gdown(2,1:3)
              write(funit,10) gdown(3,1:3)

              write(funit,*)
              write(funit,*) 'Cov: cov1      cov2      cov3'
              write(funit,10) gdef%gmetric%grid(igx)%cov(i,j,k,1:3,1)
              write(funit,10) gdef%gmetric%grid(igx)%cov(i,j,k,1:3,2)
              write(funit,10) gdef%gmetric%grid(igx)%cov(i,j,k,1:3,3)

              write(funit,*)
              write(funit,*) 'Cnv: cnv1      cnv2      cnv3'
              write(funit,10) gdef%gmetric%grid(igx)%cnv(i,j,k,1:3,1)
              write(funit,10) gdef%gmetric%grid(igx)%cnv(i,j,k,1:3,2)
              write(funit,10) gdef%gmetric%grid(igx)%cnv(i,j,k,1:3,3)

              write(funit,*)
              write(funit,*) 'Jac:',gdef%gmetric%grid(igx)%jac(i,j,k)

              write(funit,*) 
              write(funit,*) 'Metric tensor product'
              write(funit,10) tensor(1,1:3)
              write(funit,10) tensor(2,1:3)
              write(funit,10) tensor(3,1:3)
 10           format (3f10.3)

            enddo
          enddo
        enddo

c     End program

      end subroutine metricTensorCheck

c     numericalGridCoeffsCheck
c     #################################################################
      subroutine numericalGridCoeffsCheck

c     -----------------------------------------------------------------
c     Test of numerical calculation vs. analytical calculation of
c     grid quantities.
c     -----------------------------------------------------------------

      implicit none

      integer :: i,j,k,igrid,igx,igy,igz
     .          ,ilo,ihi,jlo,jhi,klo,khi
     .          ,ilom,ihip,jlom,jhip,klom,khip

      real(8) :: mag,mag1,mag2(3,3),mag3(3,3),mag4(3,3,3)

c     Begin program

      if (anal_map.and.(.not.numerical_grid)) then
        write(funit,*)
        write(funit,*)'***********************************************'
        write(funit,*)'Checking accuracy of numerical grid quantities.'
        write(funit,*)'***********************************************'
        write(funit,*)
      else
        write(funit,*)
        write(funit,*)'***********************************************'
        write(funit,*)'Cannot check numerical grid quantities:'
        write(funit,*)'   analytical mapping unavailable'
        write(funit,*)'***********************************************'
        write(funit,*)
        return
      endif

      do igrid=1,gdef%ngrid

        igx = igrid
        igy = igrid
        igz = igrid

        ilo = 1
        jlo = 1
        klo = 1
        ihi = gdef%nxv(igrid)
        jhi = gdef%nyv(igrid)
        khi = gdef%nzv(igrid)

        ilom = ilo-1
        jlom = jlo-1
        klom = klo-1
        ihip = ihi+1
        jhip = jhi+1
        khip = khi+1

        mag1 = 0d0
        mag2 = 0d0
        mag3 = 0d0
        mag4 = 0d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              mag1 = mag1
     .             + (gdef%gmetric%grid(igrid)%jac  (i,j,k)
     .               -jacobian(gdef,i,j,k,igx,igy,igz))**2
              mag2 = mag2
     .             + (gdef%gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .               -g_sub   (gdef,i,j,k,igx,igy,igz))**2
              mag3 = mag3
     .             + (gdef%gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .               -g_sup   (gdef,i,j,k,igx,igy,igz))**2
              mag4 = mag4
     .             + (gdef%gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
     .               -christ_2knd(gdef,i,j,k,igx,igy,igz))**2
            enddo
          enddo
        enddo

        write(funit,*) 'Grid level:',igrid
          
        mag = sqrt(    mag1/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1))
        write(funit,*)'   Jacobian tst=',mag
        mag = sum(sqrt(mag2/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1)))/9.
        write(funit,*)'   Gsub tst    =',mag
        mag = sum(sqrt(mag3/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1)))/9.
        write(funit,*)'   Gsup tst    =',mag
        mag = sum(sqrt(mag4/(ihi-ilo+1)/(jhi-jlo+1)/(khi-klo+1)))/27.
        write(funit,*)'   Gamma tst   =',mag
        write(funit,*)

      enddo

      end subroutine numericalGridCoeffsCheck

c     graphicsCheck_cc
c     #################################################################
      subroutine graphicsCheck_cc

c     -----------------------------------------------------------------
c     Test of numerical calculation vs. analytical calculation of
c     grid quantities.
c     -----------------------------------------------------------------

      use xdraw_io

      implicit none

      integer    :: nx,ny,nz,i,j,k,l,m,n,igrid,ivar
     .             ,nxg,nyg,nzg,ig,jg,kg
      real(8)    :: tensor(3,3,3),vec(3)

      integer   ,parameter :: nvar=100
      character(30) :: descr(nvar),drawfile,debugfile

      logical :: dump_diff

      real(8),allocatable,dimension(:,:,:) :: dbl
      real(8),pointer,dimension(:,:,:) :: jac
      real(8),pointer,dimension(:,:,:,:,:) :: cov,cnv

c     Begin program

c     Output

      write(funit,*)
      write(funit,*) 'Dumping CELL-CENTERED grid debugging graphics...'
      write(funit,*)

      if (my_rank == 0) then
        write(*  ,*)
        write(*  ,*) 'Dumping CELL-CENTERED grid debugging graphics...'
        write(*  ,*)
      endif

cc      dump_diff = anal_map.and.numerical_grid.and.(.not.map_grid_pack)
      dump_diff = .false.

      if (dump_diff) then
        if (my_rank == 0)
     .    write(*    ,*)
     .       '>> Analytical mapping available; dumping diffs...'
        write(funit,*)
     .       '>> Analytical mapping available; dumping diffs...'
      else
        if (.not.numerical_grid) then
          if (my_rank == 0)
     .      write(*    ,*)'>> Dumping analytical qtys...'
          write(funit,*)'>> Dumping analytical qtys...'
        else
          if (my_rank == 0)
     .      write(*    ,*)'>> Dumping numerical qtys...'
          write(funit,*)'>> Dumping numerical qtys...'
        endif
      endif
      write(funit,*)
      if (my_rank == 0) write(*,*)

c     Dump quantities

      do igrid=1,gdef%ngrid

        jac => gdef%gmetric%grid(igrid)%jac
        cov => gdef%gmetric%grid(igrid)%cov
        cnv => gdef%gmetric%grid(igrid)%cov

        nx = gdef%nxv(igrid)
        ny = gdef%nyv(igrid)
        nz = gdef%nzv(igrid)

        nxg = gdef%nxgl(igrid)
        nyg = gdef%nygl(igrid)
        nzg = gdef%nzgl(igrid)

        do k=0,nz+1

        write(funit,*)
        write(funit,*) '   ******************************'
        write(funit,*) '   Grid level',igrid,'on slice k=',k
        write(funit,*) '   ******************************'
        write(funit,*)

        debugfile = trim(level)//'_grid_gl'//trim(int2char(igrid))
     .                         //'_k='//trim(int2char(k))//'.bin'

        !Allocate dummy array
        allocate(dbl(0:nx+1,0:ny+1,0:nz+1))

        !Jacobian
        write(funit,*) '   Dumping Jacobian'
        ivar = 1
        descr(ivar) = 'jac'

        do j=0,ny+1
          do i=0,nx+1
            dbl(i,j,1) = gdef%gmetric%grid(igrid)%jac(i,j,k)
            if (dump_diff) dbl(i,j,1) = dbl(i,j,1)
     .                    -jacobian(gdef,i,j,k,igrid,igrid,igrid)
          enddo
        enddo

        call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl,0)

cc        !Volume
cc        dbg = 0d0
cc        if (igrid > 1) then
cc          do j=1,ny
cc            do i=1,nx
cc              dbg(i,j,1,1) =
cc     .       sum(gdef%gmetric%grid(igrid-1)%dvol(2*i-1:2*i,2*j-1:2*j,k))
cc            enddo
cc          enddo
cc          dbg(1:nx,1:ny,1,1)=dbg(1:nx,1:ny,1,1)
cc     .            -gdef%gmetric%grid(igrid)%dvol(1:nx,1:ny,k)
cc          write (*,*) 'Grid level=',igrid,'  Sum=',sum(dbg)
cc        else
cc          dbg(:,:,1,1) = gdef%gmetric%grid(igrid)%dvol(:,:,k)
cc        endif
cc        write(funit,*) '   Dumping Cell Volume'
cc        ivar = ivar+1
cc        descr(ivar) = 'dvol'
cc        call contour(dbg(:,:,1,1),nx+2,ny+2,0d0,xmax,0d0,ymax,1,110)

        !Metric tensors
        do l=1,3
          do m=l,3
            write(funit,*) '   Dumping Gsub component',l,m
            ivar = ivar+1
            descr(ivar)='gsub('//trim(int2char(l))//','
     .                         //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%gsub(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac(i,j,k)
cc     .               gdef%gmetric%grid(igrid)%gsub_v(i,j,k,l,m)
cc     .              *gdef%gmetric%grid(igrid)%jac_v(i,j,k)
                if (dump_diff) then
                  tensor(:,:,1)  = g_sub(gdef,i,j,k,igrid,igrid,igrid)
                  dbl(i,j,1) = dbl(i,j,1)
     .                    -tensor(l,m,1)
     .                    *jacobian(gdef,i,j,k,igrid,igrid,igrid)
                endif
              enddo
            enddo

            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1)
          enddo
        enddo

        do l=1,3
          do m=l,3
            write(funit,*) '   Dumping Gsup component',l,m
            ivar = ivar+1
            descr(ivar)='gsup('//trim(int2char(l))//','
     .                         //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
cc     .               gdef%gmetric%grid(igrid)%gsup_v(i,j,k,l,m)
cc     .              *gdef%gmetric%grid(igrid)%jac_v(i,j,k)
     .               gdef%gmetric%grid(igrid)%gsup(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac(i,j,k)
                if (dump_diff) then
                  tensor(:,:,1)  = g_sup(gdef,i,j,k,igrid,igrid,igrid)
                  dbl(i,j,1) = dbl(i,j,1)
     .                        -tensor(l,m,1)
     .                    *jacobian(gdef,i,j,k,igrid,igrid,igrid)
                endif
              enddo
            enddo
            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1)
          enddo
        enddo

        !Fundamental property: div(J.cov)=0
        dbl = 0d0
        do l=1,3
          write(funit,*) '   Dumping div(J.cov)',l
          ivar = ivar+1
          descr(ivar)='div(J.cov('//trim(int2char(l))//'))'

          if (k /= 0 .and. k /= nz+1) then
            do j=1,ny
              do i=1,nx
                call getMGmap(gdef,i,j,k,igrid,igrid,igrid,ig,jg,kg)

                dbl(i,j,1) =
     .                (jac(i+1,j,k)*cov(i+1,j,k,1,l)
     .                -jac(i-1,j,k)*cov(i-1,j,k,1,l))/(2*gdef%dxh(ig))
     .               +(jac(i,j+1,k)*cov(i,j+1,k,2,l)
     .                -jac(i,j-1,k)*cov(i,j-1,k,2,l))/(2*gdef%dyh(jg))
     .               +(jac(i,j,k+1)*cov(i,j,k+1,3,l)
     .                -jac(i,j,k-1)*cov(i,j,k-1,3,l))/(2*gdef%dzh(kg))
              enddo
            enddo
          endif

          call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     $                      ,1)
        enddo

        !Cnv vectors
        do l=1,3
          do m=1,3
            write(funit,*) '   Dumping Cnv vector',l,' component',m
            ivar = ivar+1
            descr(ivar)='Cnv_'//trim(int2char(l))//'('
     .                        //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%cnv(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac(i,j,k)
                if (dump_diff) then
                  vec =
     .               contravariantVector(gdef,l,i,j,k,igrid,igrid,igrid)
     .                 *jacobian(gdef,i,j,k,igrid,igrid,igrid)
                  dbl(i,j,1) = dbl(i,j,1)-vec(m)
                endif
              enddo
            enddo
            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1)
          enddo
        enddo

        !Cov vectors
        do l=1,3
          do m=1,3
            write(funit,*) '   Dumping Cov vector',l,' component',m
            ivar = ivar+1
            descr(ivar)='Cov_'//trim(int2char(l))//'('
     .                        //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%cov(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac(i,j,k)
                if (dump_diff) then
                  vec =
     .               covariantVector(gdef,l,i,j,k,igrid,igrid,igrid)
     .                 *jacobian(gdef,i,j,k,igrid,igrid,igrid)
                  dbl(i,j,1) = dbl(i,j,1)-vec(m)
                endif
              enddo
            enddo
            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1)
          enddo
        enddo

        !Christoffle symbols
        do l=1,3
          do m=1,3
            do n=m,3
              write(funit,*) '   Dumping Gamma component',l,m,n
              ivar = ivar+1
              descr(ivar)='Gamma'//trim(int2char(l))//'('
     .                           //trim(int2char(m))//','
     .                           //trim(int2char(n))//')'

              do j=0,ny+1
                do i=0,nx+1
                  dbl(i,j,1) =
     .                 gdef%gmetric%grid(igrid)%Gamma(i,j,k,l,m,n)
     .                *gdef%gmetric%grid(igrid)%jac(i,j,k)
                  if (dump_diff) then
                    tensor  = christ_2knd(gdef,i,j,k,igrid,igrid,igrid)
                    dbl(i,j,1) = dbl(i,j,1)
     .                    -tensor(l,m,n)
     .                    *jacobian(gdef,i,j,k,igrid,igrid,igrid)
                  endif
                enddo
              enddo
              call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg
     .                            ,dbl,1)
            enddo
          enddo
        enddo

        deallocate(dbl)

        if (my_rank == 0) then
          drawfile  = 'draw'//trim(level)//'_gl'//trim(int2char(igrid))
     .               //'_k='//trim(int2char(k))//'.in'
          !Create draw*.in file
          call createDrawInCfile(ivar,debugfile,'Grid debug'
     .                          ,'t','x','y',descr,'-c -X0 -L57'
     .                          ,drawfile)
        endif

        enddo      !k loop
     
      enddo        !Grid level loop

      if (my_rank == 0) then
        write(*,*)
        write(*,*) 'Done!'
        write(*,*) 
      endif
      
      write(funit,*)
      write(funit,*) 'Done!'
      write(funit,*) 

      end subroutine graphicsCheck_cc

c     graphicsCheck_vrtx
c     #################################################################
      subroutine graphicsCheck_vrtx

c     -----------------------------------------------------------------
c     Test of numerical calculation vs. analytical calculation of
c     grid vectex quantities.
c     -----------------------------------------------------------------

      use xdraw_io

      implicit none

      integer    :: nx,ny,nz,i,j,k,l,m,n,igrid,ivar
     .             ,nxg,nyg,nzg
      real(8)    :: tensor(3,3,3),vec(3)

      integer   ,parameter :: nvar=50
      character(30) :: descr(nvar),drawfile,debugfile

      logical :: dump_diff

      real(8),allocatable,dimension(:,:,:) :: dbl

c     Begin program

c     Output

      write(funit,*)
      write(funit,*) 'Dumping VERTEX grid debugging graphics...'
      write(funit,*)

      if (my_rank == 0) then
        write(*  ,*)
        write(*  ,*) 'Dumping VERTEX grid debugging graphics...'
        write(*  ,*)
      endif

cc      dump_diff = anal_map.and.numerical_grid.and.(.not.map_grid_pack)
      dump_diff = .false.

      if (dump_diff) then
        if (my_rank == 0)
     .    write(*    ,*)
     .       '>> Analytical mapping available; dumping diffs...'
        write(funit,*)
     .       '>> Analytical mapping available; dumping diffs...'
      else
        if (.not.numerical_grid) then
          if (my_rank == 0)
     .      write(*    ,*)'>> Dumping analytical qtys...'
          write(funit,*)'>> Dumping analytical qtys...'
        else
          if (my_rank == 0)
     .      write(*    ,*)'>> Dumping numerical qtys...'
          write(funit,*)'>> Dumping numerical qtys...'
        endif
      endif
      write(funit,*)
      if (my_rank == 0) write(*,*)

c     Dump quantities

      do igrid=1,gdef%ngrid

        nx = gdef%nxv(igrid)
        ny = gdef%nyv(igrid)
        nz = gdef%nzv(igrid)

        nxg = gdef%nxgl(igrid)
        nyg = gdef%nygl(igrid)
        nzg = gdef%nzgl(igrid)

        do k=0,nz+1

        write(funit,*)
        write(funit,*) '   ******************************'
        write(funit,*) '   Grid level',igrid,'on slice k=',k
        write(funit,*) '   ******************************'
        write(funit,*)

        debugfile = trim(level)//'_vrtx_gl'//trim(int2char(igrid))
     .                         //'_k='//trim(int2char(k))//'.bin'

        !Allocate dummy array
        allocate(dbl(0:nx+1,0:ny+1,0:nz+1))

        !Jacobian
        write(funit,*) '   Dumping Jacobian'
        ivar = 1
        descr(ivar) = 'jac'

        do j=0,ny+1
          do i=0,nx+1
            dbl(i,j,1) = gdef%gmetric%grid(igrid)%jac_v(i,j,k)
            if (dump_diff) dbl(i,j,1) = dbl(i,j,1)
     .             -jacobian(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
          enddo
        enddo

        call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl,0
     .                      ,no_bc=.true.)

        !Metric tensors
        do l=1,3
          do m=l,3
            write(funit,*) '   Dumping Gsub component',l,m
            ivar = ivar+1
            descr(ivar)='gsub('//trim(int2char(l))//','
     .                         //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%gsub_v(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac_v(i,j,k)
                if (dump_diff) then
                  tensor(:,:,1)  =
     .                 g_sub(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
                  dbl(i,j,1) = dbl(i,j,1)
     .                    -tensor(l,m,1)
     .             *jacobian(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
                endif
              enddo
            enddo

            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1,no_bc=.true.)
          enddo
        enddo

        do l=1,3
          do m=l,3
            write(funit,*) '   Dumping Gsup component',l,m
            ivar = ivar+1
            descr(ivar)='gsup('//trim(int2char(l))//','
     .                         //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%gsup_v(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac_v(i,j,k)
                if (dump_diff) then
                  tensor(:,:,1)  =
     .                 g_sup(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
                  dbl(i,j,1) = dbl(i,j,1)
     .                        -tensor(l,m,1)
     .             *jacobian(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
                endif
              enddo
            enddo
            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1,no_bc=.true.)
          enddo
        enddo

        !Cnv vectors
        do l=1,3
          do m=1,3
            write(funit,*) '   Dumping Cnv vector',l,' component',m
            ivar = ivar+1
            descr(ivar)='Cnv_'//trim(int2char(l))//'('
     .                        //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%cnv_v(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac_v(i,j,k)
                if (dump_diff) then
                  vec =
     .               contravariantVector(gdef,l,i,j,k,igrid,igrid,igrid
     .                                  ,vertex=.true.)
     .             *jacobian(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
                  dbl(i,j,1) = dbl(i,j,1)-vec(m)
                endif
              enddo
            enddo
            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1,no_bc=.true.)
          enddo
        enddo

        !Cov vectors
        do l=1,3
          do m=1,3
            write(funit,*) '   Dumping Cov vector',l,' component',m
            ivar = ivar+1
            descr(ivar)='Cov_'//trim(int2char(l))//'('
     .                        //trim(int2char(m))//')'

            do j=0,ny+1
              do i=0,nx+1
                dbl(i,j,1) =
     .               gdef%gmetric%grid(igrid)%cov_v(i,j,k,l,m)
     .              *gdef%gmetric%grid(igrid)%jac_v(i,j,k)
                if (dump_diff) then
                  vec =
     .               covariantVector(gdef,l,i,j,k,igrid,igrid,igrid
     .                                  ,vertex=.true.)
     .             *jacobian(gdef,i,j,k,igrid,igrid,igrid,vertex=.true.)
                  dbl(i,j,1) = dbl(i,j,1)-vec(m)
                endif
              enddo
            enddo
            call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                          ,1,no_bc=.true.)
          enddo
        enddo

        deallocate(dbl)

        if (my_rank == 0) then
          drawfile  = 'draw'//trim(level)//'_vrtx_gl'
     .               //trim(int2char(igrid))
     .               //'_k='//trim(int2char(k))//'.in'
          !Create draw*.in file
          call createDrawInCfile(ivar,debugfile,'Grid debug'
     .                          ,'t','x','y',descr,'-c -X0 -L57'
     .                          ,drawfile)
        endif

        enddo      !k loop
     
      enddo        !Grid level loop

      if (my_rank == 0) then
        write(*,*)
        write(*,*) 'Done!'
        write(*,*) 
      endif
      
      write(funit,*)
      write(funit,*) 'Done!'
      write(funit,*) 

      end subroutine graphicsCheck_vrtx

c     hessianZeroForceCheck
c     #################################################################
      subroutine hessianZeroForceCheck(igrid)

c     -----------------------------------------------------------------
c     Check zero-force condition of Christoffle symbols
c     -----------------------------------------------------------------

        implicit  none

c     Call variables

        integer :: igrid
        
c     Local variables

        integer :: ii,ll,kk,mm,jj,ig,jg,kg,nx,ny,nz,i,j,k
        real(8) :: summ,summ2,dh,jac,jacp,jacm,const
        real(8) :: hess(3,3,3),gsup(3,3),gsub(3,3)
     $            ,gsupp(3,3),gsupm(3,3),gamm(3,3,3)
        logical :: alt_eom

c     Begin program

        write (funit,*)
        write (funit,*) '********************************'
        write (funit,*) 'Checking zero-force condition...'
        write (funit,*) '********************************'
        write (funit,*)

        nx = gdef%nxv(igrid)
        ny = gdef%nyv(igrid)
        nz = gdef%nzv(igrid)

        do k=1,nz
          do j=1,ny
            do i=1,nx

              write(funit,*)
              write(funit,*) 'Grid:',igrid,'  Grid node:',i,j,k

              call getMGmap(gdef,i,j,k,igrid,igrid,igrid,ig,jg,kg)

              gsup  = gdef%gmetric%grid(igrid)%gsup(i,j,k,:,:)
              gsub  = gdef%gmetric%grid(igrid)%gsub(i,j,k,:,:)
              gamm  = gdef%gmetric%grid(igrid)%gamma(i,j,k,:,:,:)

c             Check difference between correction and prediction

cc              hess = gamm - gold
cc
cc              summ2= sqrt(sum(gold*gold))
cc              if (summ2 /= 0d0) then
cc                summ = sqrt(sum(hess*hess))/summ2
cc              else
cc                summ = 0d0
cc              endif
cc
cc              write(funit,*) '  Correction size:',summ

c             Check cancellation property g^lk(d_il G^j_kj-g_lm g^ij G^m_jk) = 0

              do ii=1,3 !cycle through Chistoffel symbols
                do ll=1,3
                  do kk=1,3
                    hess(ii,ll,kk) =
     .                     delta(ii,ll)*(gamm(1,kk,1)
     .                                  +gamm(2,kk,2)
     .                                  +gamm(3,kk,3))

                    summ=0d0
                    do mm=1,3
                      do jj=1,3
                        summ = summ
     .                       + gsub(ll,mm)*gsup(ii,jj)*gamm(mm,jj,kk)
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

                write(funit,'(a,i1,a,e10.3)')
     .             '  Cancel. of corrector term for Gamma^',ii,':',summ
              enddo

c             Check zero-force property: g^lk G^i_lk + nabla_m(g^mi) = 0

              do ii=1,3 !cycle through Chistoffel symbols

                alt_eom = (ii == 2 .and. alt__eom())

                if (alt_eom) then
                  jac  = gdef%gmetric%grid(igrid)%jac(i,j,k)
                  const= 1d0
                else
                  jac  = 1d0
                  const= 0d0
                endif

                summ = 0d0

                dh = 2*gdef%dxh(ig)
                if (alt_eom) then
                  jacp = gdef%gmetric%grid(igrid)%jac(i+1,j,k)
                  jacm = gdef%gmetric%grid(igrid)%jac(i-1,j,k)
                else
                  jacp=1d0
                  jacm=1d0
                endif
                gsupp=jacp*gdef%gmetric%grid(igrid)%gsup(i+1,j,k,:,:)
                gsupm=jacm*gdef%gmetric%grid(igrid)%gsup(i-1,j,k,:,:)

                summ = summ +(gsupp(1,ii)-gsupm(1,ii))/dh/jac

                dh = 2*gdef%dyh(jg)
                if (alt_eom) then
                  jacp = gdef%gmetric%grid(igrid)%jac(i,j+1,k)
                  jacm = gdef%gmetric%grid(igrid)%jac(i,j-1,k)
                else
                  jacp=1d0
                  jacm=1d0
                endif
                gsupp=jacp*gdef%gmetric%grid(igrid)%gsup(i,j+1,k,:,:)
                gsupm=jacm*gdef%gmetric%grid(igrid)%gsup(i,j-1,k,:,:)

                summ = summ +(gsupp(2,ii)-gsupm(2,ii))/dh/jac

                dh = 2*gdef%dzh(kg)
                if (alt_eom) then
                  jacp = gdef%gmetric%grid(igrid)%jac(i,j,k+1)
                  jacm = gdef%gmetric%grid(igrid)%jac(i,j,k-1)
                else
                  jacp=1d0
                  jacm=1d0
                endif
                gsupp=jacp*gdef%gmetric%grid(igrid)%gsup(i,j,k+1,:,:)
                gsupm=jacm*gdef%gmetric%grid(igrid)%gsup(i,j,k-1,:,:)

                summ = summ +(gsupp(3,ii)-gsupm(3,ii))/dh/jac

                summ2 = 0d0
                do ll=1,3
                  do kk=1,3
                    summ2 = summ2 + gsup(ll,kk)*gamm(ii,ll,kk)
     $                      - const*gsup(ll,ii)*gamm(kk,ll,kk)
                  enddo
                enddo

                write(funit,'(a,i1,a,e10.3)')
     .              '  Zero force cancellation for i=',ii,':',summ+summ2
              enddo

            enddo
          enddo
        enddo

c     End program

      end subroutine hessianZeroForceCheck

c     symmCheck
c     #################################################################
      subroutine symmCheck(dim)

c     -----------------------------------------------------------------
c     Symmetry test on map (only for applicable cases).
c     -----------------------------------------------------------------

      use xdraw_io

      implicit none

c     Call variables

      integer :: dim

c     Local variables

      integer    :: nx,ny,nz,i,j,k,l,m,n,igrid,ivar
     .             ,nxg,nyg,nzg
      real(8)    :: tensor(3,3,3),vec(3)
      character(30) :: debugfile,drawfile

      integer   ,parameter :: nvar=50
      character(30) :: descr(nvar)

      logical :: dump_diff

      real(8),allocatable,dimension(:,:,:) :: dbl

c     Begin program

c     Output

      write(funit,*)
      write(funit,*) '**********************************'
      write(funit,*) 'Dumping symmetry test graphics...'
      write(funit,*) '**********************************'
      write(funit,*)

      write(*    ,*)
      write(*    ,*) '**********************************'
      write(*    ,*) 'Dumping symmetry test graphics...'
      write(*    ,*) '**********************************'
      write(*    ,*)

c     Dump quantities

      do igrid=1,gdef%ngrid

        nx = gdef%nxv(igrid)
        ny = gdef%nyv(igrid)
        nz = gdef%nzv(igrid)

        nxg = gdef%nxgl(igrid)
        nyg = gdef%nygl(igrid)
        nzg = gdef%nzgl(igrid)

        do k=0,nz+1

        write(funit,*)
        write(funit,*) '   ******************************'
        write(funit,*) '   Grid level',igrid,'on slice k=',k
        write(funit,*) '   ******************************'
        write(funit,*)

        !Open .bin file
        debugfile = 'symm_grid_gl'//trim(int2char(igrid))
     .                 //'k='//trim(int2char(k))//'.bin'
        drawfile = 'drawsymm_gl'//trim(int2char(igrid))
     .              //'_k='//trim(int2char(k))//'.in'

        !Allocate dummy array
        allocate(dbl(0:nx+1,0:ny+1,0:nz+1))

        !Compare map along dimension dim
        select case(dim)
        case(1)
          write(funit,*) '   Dumping SYMM test along X'
          ivar = 1
          descr(ivar) = 'X symm test'

          do j=0,ny+1
            do i=0,nx+1
              dbl(i,j,1)=
     .             (gdef%gmetric%grid(igrid)%car(i     ,j,k,dim)
     .             +gdef%gmetric%grid(igrid)%car(nx+1-i,j,k,dim))
     .             -(xmax-xmin)
            enddo
          enddo

          call dumpGlobalField(gdef,debugfile,nx,ny,nz,nxg,nyg,nzg,dbl
     .                        ,0)

        case default
          call pstop('symmCheck','Dimension not implemented')
        end select

        deallocate(dbl)

        !Create draw*.in file
        if (my_rank == 0) then
          call createDrawInCfile(ivar,debugfile,'Grid debug'
     .                          ,'t','x','y',descr,'-c -X0 -L57'
     .                          ,drawfile)
        endif

        enddo
      enddo

      write(funit,*)
      write(funit,*) 'Done!'
      write(funit,*) 

      end subroutine symmCheck

      end subroutine checkGrid

      end module grid_create


c module grid
c #####################################################################
      module grid

        use grid_create

        use grid_operations

        use grid_mg

        use grid_diff_ops

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin !3D global domain dimensions

        integer :: nxd,nyd,nzd,nxdp,nydp,nzdp
     .            ,nxl,nyl,nzl,nxlp,nylp,nzlp

        integer :: ilo ,ihi ,jlo ,jhi ,klo ,khi
     .            ,ilom,ihip,jlom,jhip,klom,khip
     .            ,ilog ,ihig ,jlog ,jhig ,klog ,khig
     .            ,ilomg,ihipg,jlomg,jhipg,klomg,khipg

        integer :: gcw=1  !Ghost cell width

      contains

c     readGridConfig
c     #################################################################
      subroutine readGridConfig(file,dump)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: file
      logical :: dump

c     Local variables

      integer        :: dim(1),loc(1),ioerr,uinput=1234
      character*(3)  :: bcs(6)
      type(grid_pack):: gp1,gp2,gp3

c     Namelists

      !Problem setup
      namelist /setup/ nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .               ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .               ,npx,npy,npz,gp1,gp2,gp3

c     Begin program

c     Set defaults

      nxd      = 64            ! Mesh points in x-direction
      nyd      = 64            ! Mesh points in y-direction
      nzd      = 64            ! Mesh points in z-direction

      coords   = 'car'         ! Coordinate system (see grid_anal_mod.F)

      xmax     = 1d0           ! Length in x-direction
      ymax     = 1d0           ! Length in y-direction
      zmax     = 1d0           ! Length in z-direction

      xmin     = 0d0           ! Length in x-direction
      ymin     = 0d0           ! Length in y-direction
      zmin     = 0d0           ! Length in z-direction

      gparams  = 0d0           ! Array with additional grid parameters (grid-dependent)

      numerical_grid = .false. ! Whether grid metrics are calculated numerically (.true.)
                               !   or analytically.

      bcs      = (/ 'def','def','per','per','per','per' /) 
                               ! Defines topological boundary conditions:
                               !   + 'def' = default
                               !     (see applyBoundaryCondtions.F)
                               !   + 'per' = periodic
                               !   + 'spt' = singular point
                               !   + 'sym' = symmetry 
                               !     (homogeneous Neumann/Dirichlet
                               !      for even/odd quantities)

      gp1%pack = .false.       ! Do not pack in X-direction
      gp1%xp   = 0d0
      gp1%dx0  = 0d0
      gp2%pack = .false.       ! Do not pack in Y-direction
      gp2%xp   = 0d0
      gp2%dx0  = 0d0
      gp3%pack = .false.       ! Do not pack in Z-direction
      gp3%xp   = 0d0
      gp3%dx0  = 0d0
                               ! To select packing, one needs to set the fields
                               ! of gp1, gp2, gp3 as follows:
                               !    gp* = pack,xp,dx0
                               ! where:
                               !   + pack (logical): whether to pack
                               !   + xp (real): where to pack
                               !   + dx0 (real): initial grid spacing (at xp)

      mg_ratio = 2             ! MG coarsening ratio

c     Read initialization parameters

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,setup,iostat=ioerr)

      if (ioerr/=0)call pstop('readGridInput','Problem reading "setup"')

      close(unit=uinput)

c     Initialize grid packing structure

      call init_grid_pack_config(gp1,gp2,gp3)

c     Consistency check

      !Adequate BCs for collapsed dimension
      if (nxd == 1) bcs(1:2) = 'per'
      if (nyd == 1) bcs(3:4) = 'per'
      if (nzd == 1) bcs(5:6) = 'per'

c     Translate boundary conditions

      ioerr = initBCS(bcs)

      if (ioerr /= 0.and.my_rank == 0) then
        loc = 1 - mod(minloc(bcond),2)
        dim = 1+(minloc(bcond) - (1+loc))/2
        write (*,*) 'Error in defining boundary conditions'
        write (*,*) 'Undefined boundary condition in axis',dim,
     .              ', location',loc
        write (*,*) 'Aborting'
        write (*,*) bcond
        stop
      endif

c     Write initialization parameters to standard output

      if (my_rank == 0.and.dump) then
        write (*,*) '------------------'
        write (*,*) 'GRID configuration'
        write (*,*) '------------------'
        write (*,setup)
      endif

c     End program

      end subroutine readGridConfig

      end module grid

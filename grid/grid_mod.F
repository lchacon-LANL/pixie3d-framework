
c module grid
c #####################################################################
      module grid

        use grid_create

        use grid_operations

        use grid_mg

        use grid_diff_ops

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin !3D global domain dimensions

        integer :: nxd,nyd,nzd

      contains

c     readGridConfig
c     #################################################################
      subroutine readGridConfig(file,dump,no_error)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: file
      logical :: dump
      logical,optional :: no_error

c     Local variables

      integer        :: dim(1),loc(1),ioerr,uinput=1234
      character*(3)  :: bcs(6)
      type(grid_pack):: gp1,gp2,gp3
      logical        :: no_err

c     Namelists

      !Problem setup
      namelist /setup/ nxd,nyd,nzd,coords,bcs,xmax,ymax,zmax
     .               ,xmin,ymin,zmin,gparams,mg_ratio,numerical_grid
     .               ,npx,npy,npz,gp1,gp2,gp3

c     Begin program

      if(PRESENT(no_error)) then
        no_err = no_error
      else
        no_err = .false.
      endif
      
c     Set defaults

      nxd      = 1             ! Mesh points in x-direction
      nyd      = 1             ! Mesh points in y-direction
      nzd      = 1             ! Mesh points in z-direction

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

      call init_grid_pack_config(gp1)  ! To select packing, one needs to set the fields
      call init_grid_pack_config(gp2)  ! of gp1, gp2, gp3 as follows:                  
      call init_grid_pack_config(gp3)  !    gp* = pack,xp,dx0                          
                                       ! where:                                        
                                       !   + pack (logical): whether to pack           
                                       !   + xp (real): where to pack                  
                                       !   + dx0 (real): initial grid spacing (at xp)  

      mg_ratio = 2             ! MG coarsening ratio

c     Read initialization parameters

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,setup,iostat=ioerr)

      no_err = no_err.or.(ioerr==0)
      
      if (.not.no_err)
     $     call pstop('readGridConfig','Problem reading "setup"')

      close(unit=uinput)

c     Initialize grid packing structure

      call set_grid_pack_config(gp1,gp2,gp3)

c     Consistency check

      !Adequate BCs for collapsed dimension
      if (nxd == 1) bcs(1:2) = 'per'
      if (nyd == 1) bcs(3:4) = 'per'
      if (nzd == 1) bcs(5:6) = 'per'

c     Translate boundary conditions

      ioerr = initBCs(bcs)

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

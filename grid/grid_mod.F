
c module grid
c #####################################################################
      module grid

        use grid_create

        use grid_operations

        use grid_mg

        use grid_diff_ops

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
      subroutine readGridConfig(file)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: file

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

      g_pack(1) = gp1
      g_pack(2) = gp2
      g_pack(3) = gp3

c     Consistency check

      !Adequate BCs for collapsed dimension
      if (nxd == 1) bcs(1:2) = 'per'
      if (nyd == 1) bcs(3:4) = 'per'
      if (nzd == 1) bcs(5:6) = 'per'

c     Translate boundary conditions

      bcond = -1

      where (bcs == 'def')
        bcond = DEF
      elsewhere (bcs == 'per')
        bcond = PER
      elsewhere (bcs == 'spt')
        bcond = SP
      elsewhere (bcs == 'sym')
        bcond = FSYM  !Symmetry at faces (for conservation)
      elsewhere (bcs == 'equ')
        bcond = EQU
      end where

      if (minval(bcond) < 0) then
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

      if (my_rank == 0) then
        write (*,*) '------------------'
        write (*,*) 'GRID configuration'
        write (*,*) '------------------'
        write (*,setup)
      endif

c     End program

      end subroutine readGridConfig

c     crossProduct_stg_mesh
c     ###############################################################
      function crossProduct_stg_mesh(igrid,vec1,vertex1,vec2,vertex2
     .                              ,cnv_in,apply_bc)
     .         result(vec3)
c     ---------------------------------------------------------------
c     Calculates cross product of staggered vectors, and returns a
c     staggered vector at cell faces. Staggered input vectors can be
c     cell-based (vertex=.false.) or vertex-based (vertex=.true.). Also,
c     they can be contravariant (cnv_in=.true.) or covariant
c     (cnv_in=.false.). If the former, they return a covariant vector,
c     and contravariant otherwise.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: vec1(0:,0:,0:,:)
     .          ,vec2(0:,0:,0:,:)
     .          ,vec3(0:size(vec1,1)-1
     .               ,0:size(vec1,2)-1
     .               ,0:size(vec1,3)-1
     .               ,  size(vec1,4))

      logical :: vertex1,vertex2,cnv_in,apply_bc

c     Local variables

      integer :: nx,ny,nz,bcnd(6,3)

c     Begin program

      nx = size(vec1,1)-2
      ny = size(vec1,2)-2
      nz = size(vec1,3)-2

      vec3 = 0d0

      if (.not.vertex1.and.vertex2) then
        !At (i+1/2,j,k)
        vec3(1:nx,1:ny,1:nz,1) =(0.25*(vec1(1:nx  ,1:ny  ,1:nz,2)
     .                                +vec1(2:nx+1,1:ny  ,1:nz,2)
     .                                +vec1(1:nx  ,0:ny-1,1:nz,2)
     .                                +vec1(2:nx+1,0:ny-1,1:nz,2))
     .                          * 0.5*(vec2(1:nx,1:ny  ,1:nz,3)
     .                                +vec2(1:nx,0:ny-1,1:nz,3))
     .                         -  0.5*(vec1(1:nx  ,1:ny,1:nz,3)
     .                                +vec1(2:nx+1,1:ny,1:nz,3))
     .                               * vec2(1:nx  ,1:ny,1:nz,2)   )

        !At (i,j+1/2,k)
        vec3(1:nx,1:ny,1:nz,2) = (0.5*(vec1(1:nx,1:ny  ,1:nz,3)
     .                                +vec1(1:nx,2:ny+1,1:nz,3))
     .                               * vec2(1:nx,1:ny  ,1:nz,1)
     .                           -.25*(vec1(1:nx  ,1:ny  ,1:nz,1)
     .                                +vec1(1:nx  ,2:ny+1,1:nz,1)
     .                                +vec1(0:nx-1,1:ny  ,1:nz,1)
     .                                +vec1(0:nx-1,2:ny+1,1:nz,1))
     .                           *0.5*(vec2(1:nx  ,1:ny,1:nz,3)
     .                                +vec2(0:nx-1,1:ny,1:nz,3)))

        !At (i,j,k)
        vec3(1:nx,1:ny,1:nz,3) =(0.5*(vec1(1:nx  ,1:ny,1:nz,1)
     .                               +vec1(0:nx-1,1:ny,1:nz,1))
     .                          *0.5*(vec2(1:nx  ,1:ny,1:nz,2)
     .                               +vec2(0:nx-1,1:ny,1:nz,2))
     .                          -0.5*(vec1(1:nx,1:ny  ,1:nz,2)
     .                               +vec1(1:nx,0:ny-1,1:nz,2))
     .                          *0.5*(vec2(1:nx,1:ny  ,1:nz,1)
     .                               +vec2(1:nx,0:ny-1,1:nz,1)))
      else
        call pstop("crossProduct_stg_mesh"
     .            ,"Case not implemented yet")
      endif

      !Extrapolate to ghost cells
      if (apply_bc) then
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond
        where (bcnd == DEF) bcnd = EXT
        call setMGBC_stg(0,3,nx,ny,nz,igrid,vec3,bcnd,iorder=1)
      endif

      !Apply Jacobian factor
      if (cnv_in) then
        vec3(:,:,:,1) = vec3(:,:,:,1)/gmetric%grid(igrid)%jac
        vec3(:,:,:,2) = vec3(:,:,:,2)/gmetric%grid(igrid)%jac
        vec3(:,:,:,3) = vec3(:,:,:,3)/gmetric%grid(igrid)%jac
      endif

c     End program

      end function crossProduct_stg_mesh

c     scalarProduct_stg_mesh
c     ###############################################################
      function scalarProduct_stg_mesh(igrid,vec1,vertex1,vec2,vertex2
     .                               ,apply_bc)
     .         result(scl)
c     ---------------------------------------------------------------
c     Calculates dot product of staggered vectors, and returns a
c     three scalar values at cell faces. Staggered input vectors can be
c     cell-based (vertex=.false.) or vertex-based (vertex=.true.).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer :: igrid

      real(8) :: vec1(0:,0:,0:,:)
     .          ,vec2(0:,0:,0:,:)
     .          ,scl (0:size(vec1,1)-1
     .               ,0:size(vec1,2)-1
     .               ,0:size(vec1,3)-1
     .               ,3)

      logical :: vertex1,vertex2
      logical,optional :: apply_bc

c     Local variables

      integer :: nx,ny,nz,bcnd(6,3)
      logical :: app_bc

c     Begin program

      if (PRESENT(apply_bc)) then
        app_bc = apply_bc
      else
        app_bc = .false.
      endif

      nx = size(vec1,1)-2
      ny = size(vec1,2)-2
      nz = size(vec1,3)-2

      scl = 0d0

      if (.not.vertex1.and.vertex2) then
!     At (i+1/2,j,k)
         scl(1:nx,1:ny,1:nz,1) = vec1(1:nx  ,1:ny  ,1:nz,1)*
     .                   0.25d0*(vec2(1:nx  ,1:ny  ,1:nz,1)+
     .                           vec2(1:nx  ,0:ny-1,1:nz,1)+         
     .                           vec2(2:nx+1,0:ny-1,1:nz,1)+         
     .                           vec2(2:nx+1,1:ny  ,1:nz,1)) +
     .                           vec2(1:nx  ,1:ny  ,1:nz,2)*
     .                   0.25d0*(vec1(1:nx  ,1:ny  ,1:nz,2)+
     .                           vec1(1:nx  ,0:ny-1,1:nz,2)+         
     .                           vec1(2:nx+1,0:ny-1,1:nz,2)+         
     .                           vec1(2:nx+1,1:ny  ,1:nz,2)) + 
     .                    0.5d0*(vec1(1:nx  ,1:ny  ,1:nz,3)+
     .                           vec1(2:nx+1,1:ny  ,1:nz,3))*
     .                    0.5d0*(vec2(1:nx  ,1:ny  ,1:nz,3)+
     .                           vec2(1:nx  ,0:ny-1,1:nz,3))

!     At (i,j+1/2,k)
         scl(1:nx,1:ny,1:nz,2) = vec2(1:nx  ,1:ny  ,1:nz,1)*
     .                   0.25d0*(vec1(1:nx  ,1:ny  ,1:nz,1)+
     .                           vec1(1:nx  ,2:ny+1,1:nz,1)+
     .                           vec1(0:nx-1,2:ny+1,1:nz,1)+
     .                           vec1(0:nx-1,1:ny  ,1:nz,1)) +
     .                           vec1(1:nx  ,1:ny  ,1:nz,2)*
     .                   0.25d0*(vec2(1:nx  ,1:ny  ,1:nz,2)+
     .                           vec2(1:nx  ,2:ny+1,1:nz,2)+ 
     .                           vec2(0:nx-1,2:ny+1,1:nz,2)+ 
     .                           vec2(0:nx-1,1:ny  ,1:nz,2)) +
     .                    0.5d0*(vec2(1:nx  ,1:ny  ,1:nz,3)+
     .                           vec2(0:nx-1,1:ny  ,1:nz,3))*
     .                    0.5d0*(vec1(1:nx  ,1:ny  ,1:nz,3)+
     .                           vec1(1:nx  ,2:ny+1,1:nz,3))


!     At (i,j,k)
         scl(1:nx,1:ny,1:nz,3) = 
     .                    0.5d0*(vec2(1:nx  ,1:ny  ,1:nz,1)+
     .                           vec2(1:nx  ,0:ny-1,1:nz,1))*
     .                    0.5d0*(vec1(1:nx  ,1:ny  ,1:nz,1)+
     .                           vec1(0:nx-1,1:ny  ,1:nz,1)) +
     .                    0.5d0*(vec1(1:nx  ,1:ny  ,1:nz,2)+
     .                           vec1(1:nx  ,0:ny-1,1:nz,2))*
     .                    0.5d0*(vec2(1:nx  ,1:ny  ,1:nz,2)+
     .                           vec2(0:nx-1,1:ny  ,1:nz,2)) +
     .                           vec1(1:nx  ,1:ny  ,1:nz,3)*
     .                   0.25d0*(vec2(1:nx  ,1:ny  ,1:nz,3)+
     .                           vec2(0:nx-1,1:ny  ,1:nz,3)+ 
     .                           vec2(0:nx-1,0:ny-1,1:nz,3)+ 
     .                           vec2(1:nx  ,0:ny-1,1:nz,3)) 

      else
        call pstop("scalarProduct_stg_mesh"
     .            ,"Case not implemented yet")
      endif

      !Apply Jacobian factor
      scl(:,:,:,1) = scl(:,:,:,1)/gmetric%grid(igrid)%jac
      scl(:,:,:,2) = scl(:,:,:,2)/gmetric%grid(igrid)%jac
      scl(:,:,:,3) = scl(:,:,:,3)/gmetric%grid(igrid)%jac

      !Extrapolate to ghost cells
      if (app_bc) then
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond
        where (bcnd == DEF) bcnd = EXT
        call setMGBC_stg(0,3,nx,ny,nz,igrid,scl,bcnd,iorder=1)
      endif

c     End program

      end function scalarProduct_stg_mesh

      end module grid

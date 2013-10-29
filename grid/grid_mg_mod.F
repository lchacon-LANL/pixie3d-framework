c module grid_mg
c######################################################################
      module grid_mg

        use grid_mpi

        use grid_operations

        integer    :: ngrdx,ngrdy,ngrdz,ngrid,igmax

cc        integer   ,dimension(:),allocatable ::
cc     .                     istart,ntotv,istartp,ntotvp
cc     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
cc     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z
        integer   ,dimension(:),pointer ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z

        type(grid_mg_def),pointer :: MGgrid,MGgrid_sv

        type :: mg_ctx
          integer :: ngrdx,ngrdy,ngrdz,ngrid,igmax
     .              ,imin,imax,jmin,jmax,kmin,kmax
          integer,pointer,dimension(:) ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z
          integer   ,pointer,dimension(:)  :: iline       !Restrict ops. to i=iline in MG
          integer   ,pointer,dimension(:)  :: jline       !Restrict ops. to j=jline in MG
          integer   ,pointer,dimension(:)  :: kline       !Restrict ops. to k=kline in MG
        end type mg_ctx

        type(mg_ctx),pointer :: mg_ctx_temp

      contains

ccc     allocPointers
ccc     #################################################################
cc      subroutine allocPointers(neq,fpointers,mg_grid)
cc
ccc     -----------------------------------------------------------------
ccc     Initializes pointers for block MG. In call sequence:
ccc       * neq (in), integer: number of coupled equations to be solved
ccc       * mg_grid (in/out), type(grid_mg_def),pointer: contains definition of MG grid
ccc            level structure.
ccc       * fpointers (out), logical: indicates if pointers were allocated
ccc            or not.
ccc     -----------------------------------------------------------------
cc
cc      implicit none        !For safe fortran
cc
ccc     Call variables
cc
cc      integer             :: neq
cc      logical,intent(OUT) :: fpointers
cc      type(grid_mg_def),pointer,optional :: mg_grid
cc
ccc     Local variables
cc
cc      integer    :: i
cc
ccc     Begin program
cc
cc      fpointers = .false.
cc
cc      if (PRESENT(mg_grid)) then
cc        MGgrid=>mg_grid
cc      else
cc        if (.not.associated(MGgrid)) MGgrid => grid_params
cc      endif
cc
cc      ngrdx = MGgrid%ngrdx
cc      ngrdy = MGgrid%ngrdy
cc      ngrdz = MGgrid%ngrdz
cc      ngrid = MGgrid%ngrid
cc
ccc     Check if MG pointers are allocated
cc
cc      if (.not.allocated(istart)) then
cc
cc        allocate(istart (ngrid+1),ntotv (ngrid)
cc     .          ,istartp(ngrid+1),ntotvp(ngrid)
cc     .          ,istartb(ngrid+1),ntotb (ngrid)
cc     .          ,nxv(ngrid),nyv(ngrid),nzv(ngrid),nblock(ngrid)
cc     .          ,mg_ratio_x(ngrid),mg_ratio_y(ngrid),mg_ratio_z(ngrid))
cc
cc        mg_ratio_x = MGgrid%mg_ratio_x
cc        mg_ratio_y = MGgrid%mg_ratio_y
cc        mg_ratio_z = MGgrid%mg_ratio_z
cc
cc        nxv = MGgrid%nxv
cc        nyv = MGgrid%nyv
cc        nzv = MGgrid%nzv
cc
cc        istartp(1) = 1
cc        istart (1) = 1
cc        istartb(1) = 1
cc        nblock (1) = mg_ratio_x(1)*mg_ratio_y(1)*mg_ratio_z(1)
cc        ntotvp (1) = nxv(1)*nyv(1)*nzv(1)
cc        ntotv  (1) = neq*ntotvp(1)
cc        ntotb  (1) = nblock(1)*ntotv(1)
cc        do i = 2,ngrid
cc          nblock (i) = mg_ratio_x(i)*mg_ratio_y(i)*mg_ratio_z(i)
cc          ntotvp (i) = nxv(i)*nyv(i)*nzv(i)
cc          ntotv  (i) = neq*ntotvp(i)
cc          ntotb  (i) = nblock(i)*ntotv(i)
cc        enddo
cc        do i = 2,ngrid+1
cc          istartp(i) = istartp(i-1) + ntotvp(i-1)
cc          istart (i) = istart (i-1) + ntotv (i-1)
cc          istartb(i) = istartb(i-1) + ntotb (i-1)
cc        enddo
cc
cc        fpointers = .true.
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine allocPointers
cc
ccc     deallocPointers
ccc     #################################################################
cc      subroutine deallocPointers(fpointers)
cc
ccc     -----------------------------------------------------------------
ccc     Deallocates pointers for block MG if fpointers (in, logical) is
ccc     true.
ccc     -----------------------------------------------------------------
cc
cc      implicit none        !For safe fortran
cc
cc      logical :: fpointers
cc
cc      if (fpointers) then
cc        deallocate(istart,ntotv,istartp,ntotvp
cc     .            ,istartb,ntotb,nxv,nyv,nzv,nblock
cc     .            ,mg_ratio_x,mg_ratio_y,mg_ratio_z)
cc      endif
cc
cc      end subroutine deallocPointers
cc
cc
ccc     limits
ccc     ###############################################################
cc      subroutine limits(elem,nx,ny,nz,igr,imin,imax,jmin,jmax,kmin,kmax)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Finds limits on loops for matvec routines. In call sequence:
ccc       * elem: if positive, defines lexicographic position of a
ccc          grid node.
ccc       * nx,ny,nz (in): grid dimensions
ccc       * imin,imax,jmin,jmax,kmin,kmax: loop limits.
ccc          If elem <=0, these sample the whole grid.
ccc          If elem > 0, these sample stencil around grid node defined
ccc             by elem.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax,igr
cc
ccc     Local variables
cc
cc      integer    :: el1
cc
ccc     Begin program
cc
cc      if (elem > 0) then      !Give node coordinates on grid
cc        el1  = mod(elem,nx*ny)
cc        if (el1 == 0) el1 = nx*ny
cc
cc        imin = mod(el1 ,nx)
cc        if (imin == 0) imin = nx
cc        imax = imin
cc
cc        jmin = 1 + (el1 - imin)/nx
cc        jmax = jmin
cc
cc        kmin = 1 + (elem - imin - nx*(jmin-1))/(nx*ny)
cc        kmax = kmin
cc      elseif (elem <= 0) then !Give grid coordinates
cc        imin = 1
cc        imax = nx
cc        jmin = 1
cc        jmax = ny
cc        kmin = 1
cc        kmax = nz
cc
ccc       Check for plane/line definitions
cc
cc        if ( MGgrid%iline(igr) /= 0 ) then
cc          imin = MGgrid%iline(igr)
cc          imax = imin
cc        endif
cc        if ( MGgrid%jline(igr) /= 0 ) then
cc          jmin = MGgrid%jline(igr)
cc          jmax = jmin
cc        endif
cc        if ( MGgrid%kline(igr) /= 0 ) then
cc          kmin = MGgrid%kline(igr)
cc          kmax = kmin
cc        endif
cc      endif
cc
ccc     End program
cc
cc      end subroutine limits

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,fpointers,mgctx)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in): number of coupled equations to be solved
c       * fpointers (out): indicates if pointers were allocated or not.
c       * mgctx (in/out): contains definition of MG grid level structure.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer             :: neq
      logical,intent(OUT) :: fpointers
      type(mg_ctx),pointer,optional :: mgctx

c     Local variables

c     Begin program

      if (.not.associated(MGgrid)) MGgrid => grid_params

      if (PRESENT(mgctx)) then

cc        nullify(mg_ctx_sv)

        if (associated(mgctx)) then
cc          if (associated(mg_ctx_temp)) mg_ctx_sv => mg_ctx_temp  !Save old mg_ctx_temp
          mg_ctx_temp => mgctx
          fpointers = .false.
        else
          call allocMGPointers(neq,fpointers,mg_ctx_temp)
        endif

      else
        call allocMGPointers(neq,fpointers,mg_ctx_temp)
      endif

      istart  => mg_ctx_temp%istart
      istartp => mg_ctx_temp%istartp
      istartb => mg_ctx_temp%istartb
      ntotv   => mg_ctx_temp%ntotv
      ntotvp  => mg_ctx_temp%ntotvp
      ntotb   => mg_ctx_temp%ntotb
      nxv     => mg_ctx_temp%nxv
      nyv     => mg_ctx_temp%nyv
      nzv     => mg_ctx_temp%nzv
      nblock  => mg_ctx_temp%nblock
      mg_ratio_x => mg_ctx_temp%mg_ratio_x
      mg_ratio_y => mg_ctx_temp%mg_ratio_y
      mg_ratio_z => mg_ctx_temp%mg_ratio_z

      ngrid = mg_ctx_temp%ngrid

c     End program

      end subroutine allocPointers

c     deallocPointers
c     #################################################################
      subroutine deallocPointers(fpointers)

c     -----------------------------------------------------------------
c     Deallocates pointers for block MG if fpointers (in, logical) is
c     true.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers

      if (fpointers) then
        call deallocMGPointers(fpointers,mg_ctx_temp)

        nullify(istart,ntotv,istartp,ntotvp
     .         ,istartb,ntotb,nxv,nyv,nzv,nblock
     .         ,mg_ratio_x,mg_ratio_y,mg_ratio_z,MGgrid)
      endif

      !Recover parent MG config
cc      if (associated(mg_ctx_sv)) then
cc        mg_ctx_temp => mg_ctx_sv
cc        nullify(mg_ctx_sv)
cc      endif

      end subroutine deallocPointers

c     limits
c     ###############################################################
      subroutine limits(elem,nx,ny,nz,igr,imin,imax,jmin,jmax,kmin,kmax)
      implicit none
c     ---------------------------------------------------------------
c     Finds limits on loops for matvec routines. In call sequence:
c       * elem: if positive, defines lexicographic position of a
c          grid node.
c       * nx,ny,nz (in): grid dimensions
c       * imin,imax,jmin,jmax,kmin,kmax: loop limits.
c          If elem <=0, these sample the whole grid.
c          If elem > 0, these sample stencil around grid node defined
c             by elem.
c     ---------------------------------------------------------------

c     Call variables

      integer :: elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax,igr

c     Local variables

      integer :: el1

c     Begin program

      call MGlimits(elem,igr,mg_ctx_temp)

      imin = mg_ctx_temp%imin
      jmin = mg_ctx_temp%jmin
      kmin = mg_ctx_temp%kmin

      imax = mg_ctx_temp%imax
      jmax = mg_ctx_temp%jmax
      kmax = mg_ctx_temp%kmax

c     End program

      end subroutine limits

c     allocMGPointers
c     #################################################################
      subroutine allocMGPointers(neq,fpointers,mgctx)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in), integer: number of coupled equations to be solved
c       * ngrid: number of grid levels
c       * mgx,mgy,mgz: mg coarsening ratios
c       * nxv,nyv,nzv: dimensions for grids in all levels
c       * fpointers (out), logical: indicates if pointers were allocated
c            or not.
c       * mgctx: MG context
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer             :: neq
      logical,intent(OUT) :: fpointers

      type(mg_ctx),pointer:: mgctx

c     Local variables

      integer    :: i
      logical    :: default

c     Begin program

      fpointers = .not.associated(mgctx)

c     Check if MG pointers are allocated

      if (fpointers) then  !Defaults
        allocate(mgctx)

        mgctx%ngrid = grid_params%ngrid

        allocate(mgctx%istart    (mgctx%ngrid+1)
     .          ,mgctx%istartb   (mgctx%ngrid+1)
     .          ,mgctx%istartp   (mgctx%ngrid+1)
     .          ,mgctx%ntotv     (mgctx%ngrid)
     .          ,mgctx%ntotvp    (mgctx%ngrid)
     .          ,mgctx%ntotb     (mgctx%ngrid)
     .          ,mgctx%nxv       (mgctx%ngrid)
     .          ,mgctx%nyv       (mgctx%ngrid)
     .          ,mgctx%nzv       (mgctx%ngrid)
     .          ,mgctx%iline     (mgctx%ngrid)
     .          ,mgctx%jline     (mgctx%ngrid)
     .          ,mgctx%kline     (mgctx%ngrid)
     .          ,mgctx%nblock    (mgctx%ngrid)
     .          ,mgctx%mg_ratio_x(mgctx%ngrid)
     .          ,mgctx%mg_ratio_y(mgctx%ngrid)
     .          ,mgctx%mg_ratio_z(mgctx%ngrid))

        mgctx%iline = grid_params%iline
        mgctx%jline = grid_params%jline
        mgctx%kline = grid_params%kline

        mgctx%mg_ratio_x = grid_params%mg_ratio_x
        mgctx%mg_ratio_y = grid_params%mg_ratio_y
        mgctx%mg_ratio_z = grid_params%mg_ratio_z

        mgctx%nxv = grid_params%nxv
        mgctx%nyv = grid_params%nyv
        mgctx%nzv = grid_params%nzv

        mgctx%istartp(1) = 1
        mgctx%istart (1) = 1
        mgctx%istartb(1) = 1
        mgctx%nblock (1) = mgctx%mg_ratio_x(1)
     .                    *mgctx%mg_ratio_y(1)
     .                    *mgctx%mg_ratio_z(1)
        mgctx%ntotvp (1) = mgctx%nxv(1)*mgctx%nyv(1)*mgctx%nzv(1)
        mgctx%ntotv  (1) = neq*mgctx%ntotvp(1)
        mgctx%ntotb  (1) = mgctx%nblock(1)*mgctx%ntotv(1)
        do i = 2,mgctx%ngrid
          mgctx%nblock (i) = mgctx%mg_ratio_x(i)
     .                      *mgctx%mg_ratio_y(i)
     .                      *mgctx%mg_ratio_z(i)
          mgctx%ntotvp (i) = mgctx%nxv(i)*mgctx%nyv(i)*mgctx%nzv(i)
          mgctx%ntotv  (i) = neq*mgctx%ntotvp(i)
          mgctx%ntotb  (i) = mgctx%nblock(i)*mgctx%ntotv(i)
        enddo
        do i = 2,mgctx%ngrid+1
          mgctx%istartp(i) = mgctx%istartp(i-1) + mgctx%ntotvp(i-1)
          mgctx%istart (i) = mgctx%istart (i-1) + mgctx%ntotv (i-1)
          mgctx%istartb(i) = mgctx%istartb(i-1) + mgctx%ntotb (i-1)
        enddo

      endif

c     End program

      end subroutine allocMGPointers

c     deallocMGPointers
c     #################################################################
      subroutine deallocMGPointers(fpointers,mgctx)

c     -----------------------------------------------------------------
c     Deallocates pointers for block MG if fpointers (in, logical) is
c     true.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers
      type(mg_ctx),pointer :: mgctx

      if (fpointers .and. associated(mgctx)) then
        deallocate(mgctx%istart ,mgctx%ntotv
     .            ,mgctx%istartp,mgctx%ntotvp
     .            ,mgctx%istartb,mgctx%ntotb ,mgctx%nblock
     .            ,mgctx%nxv,mgctx%nyv,mgctx%nzv
     .            ,mgctx%iline,mgctx%jline,mgctx%kline
     .            ,mgctx%mg_ratio_x,mgctx%mg_ratio_y,mgctx%mg_ratio_z)
        deallocate(mgctx)
      endif

      end subroutine deallocMGPointers

c     MGlimits
c     ###############################################################
      subroutine MGlimits(elem,igr,mgctx)
      implicit none
c     ---------------------------------------------------------------
c     Finds limits on loops for matvec routines. In call sequence:
c       * elem: if positive, defines lexicographic position of a
c          grid node.
c       * nx,ny,nz (in): grid dimensions
c       * imin,imax,jmin,jmax,kmin,kmax: loop limits.
c          If elem <=0, these sample the whole grid.
c          If elem > 0, these sample stencil around grid node defined
c             by elem.
c     ---------------------------------------------------------------

c     Call variables

      integer    :: elem,igr
      type(mg_ctx),pointer :: mgctx

c     Local variables

      integer    :: el1,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax

c     Begin program

      nx = mgctx%nxv(igr)
      ny = mgctx%nyv(igr)
      nz = mgctx%nzv(igr)

      if (elem > 0) then      !Give node coordinates on grid
        el1  = mod(elem,nx*ny)
        if (el1 == 0) el1 = nx*ny

        imin = mod(el1 ,nx)
        if (imin == 0) imin = nx
        imax = imin

        jmin = 1 + (el1 - imin)/nx
        jmax = jmin

        kmin = 1 + (elem - imin - nx*(jmin-1))/(nx*ny)
        kmax = kmin
      elseif (elem <= 0) then !Give grid coordinates
        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
        kmin = 1
        kmax = nz

c       Check for plane/line definitions

        if ( mgctx%iline(igr) /= 0 ) then
          imin = mgctx%iline(igr)
          imax = imin
        endif
        if ( mgctx%jline(igr) /= 0 ) then
          jmin = mgctx%jline(igr)
          jmax = jmin
        endif
        if ( mgctx%kline(igr) /= 0 ) then
          kmin = mgctx%kline(igr)
          kmax = kmin
        endif
      endif

      mgctx%imin = imin
      mgctx%jmin = jmin
      mgctx%kmin = kmin

      mgctx%imax = imax
      mgctx%jmax = jmax
      mgctx%kmax = kmax

c     End program

      end subroutine MGlimits

      end module grid_mg

c module grid_mg
c######################################################################
      module grid_mg

        use grid_mpi

        use grid_operations

        integer :: ngrdx,ngrdy,ngrdz,ngrid

        integer,dimension(:),pointer ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z

        type(grid_mg_def),pointer :: MGgrid    => null()
     $                              ,MGgrid_sv => null()

        type :: mg_ctx
          type(grid_mg_def),pointer :: g_def
          integer :: mpi_comm
          integer :: neq,ngrdx,ngrdy,ngrdz,ngrid,igmax
     .              ,imin,imax,jmin,jmax,kmin,kmax
          integer,pointer,dimension(:) ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z
          integer   ,pointer,dimension(:)  :: iline       !Restrict ops. to i=iline in MG
          integer   ,pointer,dimension(:)  :: jline       !Restrict ops. to j=jline in MG
          integer   ,pointer,dimension(:)  :: kline       !Restrict ops. to k=kline in MG
          real(8) :: lxmin,lxmax,lymin,lymax,lzmin,lzmax
        end type mg_ctx

        type(mg_ctx),pointer :: mg_ctx_temp => null()

      contains

c     MG_levels
c     #################################################################
      function MG_levels(nx,mgr) result(ngrd)

c     -----------------------------------------------------------------
c     Determines available MG levels based on grid resolution nx
c     and MG ratio mgr
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,mgr,ngrd

c     Local variables

      integer :: nl,i

c     Begin program

      nl = int(dlog(1d0*nx)/dlog(1d0*mgr)+0.001)

      ngrd = max(nl-1,1)  !Back off one level

      do i = ngrd,1,-1
        nl = nx/mgr**(i-1)
        if (nl*mgr**(i-1) == nx) exit
      enddo

      ngrd = i

      end function MG_levels

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,fpointers,mgctx,g_def)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in): number of coupled equations to be solved
c       * fpointers (out): indicates if pointers were allocated or not.
c       * mgctx (in/out): contains definition of MG grid level structure.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer :: neq
      logical,intent(OUT) :: fpointers
      type(mg_ctx),pointer,optional :: mgctx
      type(grid_mg_def),pointer,optional :: g_def

c     Local variables

c     Begin program

      if (.not.associated(MGgrid)) then
        if (PRESENT(g_def)) then
          MGgrid => g_def
        else
          MGgrid => grid_params
        endif
      endif

      if (PRESENT(mgctx)) then

        if (associated(mgctx)) then
          fpointers = .false.

          MGgrid => mgctx%g_def
        else
          call allocMGPointers(neq,fpointers,mgctx)
        endif

        mg_ctx_temp => mgctx

      else
        call allocMGPointers(neq,fpointers,mg_ctx_temp)
        if (PRESENT(mgctx)) mgctx => mg_ctx_temp
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
      subroutine deallocPointers(fpointers,mgctx)

c     -----------------------------------------------------------------
c     Deallocates pointers for block MG if fpointers (in, logical) is
c     true.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers
      type(mg_ctx),pointer,optional :: mgctx

      if (fpointers) then
        if (PRESENT(mgctx)) then
          call deallocMGPointers(fpointers,mgctx)
          nullify(mg_ctx_temp)
        else
          call deallocMGPointers(fpointers,mg_ctx_temp)
        endif

        nullify(istart,ntotv,istartp,ntotvp
     .         ,istartb,ntotb,nxv,nyv,nzv,nblock
     .         ,mg_ratio_x,mg_ratio_y,mg_ratio_z,MGgrid)
      endif

      end subroutine deallocPointers

c     allocMGPointers
c     #################################################################
      subroutine allocMGPointers(neq,fpointers,mgctx)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in), integer: number of coupled equations to be solved
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

        mgctx%g_def => MGgrid

        mgctx%neq = neq

        mgctx%ngrid = MGgrid%ngrid

        mgctx%igmax = mgctx%ngrid

        mgctx%mpi_comm = MGgrid%mpi_comm

        mgctx%lxmin = MGgrid%lxmin
        mgctx%lxmax = MGgrid%lxmax
        mgctx%lymin = MGgrid%lymin
        mgctx%lymax = MGgrid%lymax
        mgctx%lzmin = MGgrid%lzmin
        mgctx%lzmax = MGgrid%lzmax

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

        mgctx%iline = MGgrid%iline
        mgctx%jline = MGgrid%jline
        mgctx%kline = MGgrid%kline

        mgctx%mg_ratio_x = MGgrid%mg_ratio_x
        mgctx%mg_ratio_y = MGgrid%mg_ratio_y
        mgctx%mg_ratio_z = MGgrid%mg_ratio_z

        mgctx%nxv = MGgrid%nxv
        mgctx%nyv = MGgrid%nyv
        mgctx%nzv = MGgrid%nzv

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
        nullify(mgctx%g_def)
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
c          If elem <=0, these sample the whole grid.
c          If elem > 0, these sample stencil around grid node defined
c             by elem.
c       * igr: grid level
c       * mgctx: MG context
c     ---------------------------------------------------------------

c     Call variables

      integer :: elem,igr
      type(mg_ctx),pointer :: mgctx

c     Local variables

      integer :: el1,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax

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

      end module grid_mg

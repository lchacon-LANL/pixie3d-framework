c module grid_mg
c######################################################################
      module grid_mg

        use grid_mpi

        use grid_operations

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

        type(grid_mg_def),pointer,private :: g_mg_def

#if defined(coarse_MG)
        integer,private :: mg_proc_coarsen_factor=4  !Coarsening ratio in # procs/dim
                                                     !in parallel MG solve
        integer,private :: MPI_COMM_SAV

        logical :: coarse_mg_debug=.false.
#endif

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

c     allocMGPointers
c     #################################################################
      subroutine allocMGPointers(neq,g_def,mgctx,fpointers)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in), integer: number of coupled equations to be solved
c       * fpointers (out), logical: indicates if pointers were allocated
c            or not.
c       * mgctx: MG context
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer,intent(IN)  :: neq
      logical,intent(OUT) :: fpointers

      type(mg_ctx),pointer,intent(OUT):: mgctx
      type(grid_mg_def),pointer,intent(IN) :: g_def

c     Local variables

      integer    :: i
      logical    :: default

c     Begin program

      fpointers = .not.associated(mgctx)

c     Check if MG pointers are allocated

      if (fpointers) then  !Defaults
        allocate(mgctx)

        mgctx%g_def => g_def

        mgctx%neq = neq

        mgctx%ngrid = g_def%ngrid

        mgctx%igmax = mgctx%ngrid

#if defined(petsc)
        mgctx%mpi_comm = g_def%mpi_comm
#endif

        mgctx%lxmin = g_def%lxmin
        mgctx%lxmax = g_def%lxmax
        mgctx%lymin = g_def%lymin
        mgctx%lymax = g_def%lymax
        mgctx%lzmin = g_def%lzmin
        mgctx%lzmax = g_def%lzmax

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

        mgctx%iline = g_def%iline
        mgctx%jline = g_def%jline
        mgctx%kline = g_def%kline

        mgctx%mg_ratio_x = g_def%mg_ratio_x
        mgctx%mg_ratio_y = g_def%mg_ratio_y
        mgctx%mg_ratio_z = g_def%mg_ratio_z

        mgctx%nxv = g_def%nxv
        mgctx%nyv = g_def%nyv
        mgctx%nzv = g_def%nzv

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

cc      if (fpointers .and. associated(mgctx)) then
      if (fpointers) then
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
      subroutine limits(mgctx,elem,igr,imin,imax,jmin,jmax,kmin,kmax)
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

      integer :: elem,imin,imax,jmin,jmax,kmin,kmax,igr

      type(mg_ctx),pointer :: mgctx

c     Local variables

      integer :: el1

c     Begin program

      call MGlimits(elem,igr,mgctx)

      imin = mgctx%imin
      jmin = mgctx%jmin
      kmin = mgctx%kmin

      imax = mgctx%imax
      jmax = mgctx%jmax
      kmax = mgctx%kmax

c     End program

      end subroutine limits

#if defined(petsc) && defined(coarse_MG)
c     MG_set_proc_coarsen_factor
c     #################################################################
      subroutine MG_set_proc_coarsen_factor(mgcf)

c     -----------------------------------------------------------------
c     Determines available MG levels based on grid resolution nx
c     and MG ratio mgr
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: mgcf

c     Local variables

c     Begin program

      mg_proc_coarsen_factor = mgcf

      end subroutine MG_set_proc_coarsen_factor

c     MG_get_proc_coarsen_factor
c     #################################################################
      function MG_get_proc_coarsen_factor() result(mgcf)

c     -----------------------------------------------------------------
c     Determines available MG levels based on grid resolution nx
c     and MG ratio mgr
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: mgcf

c     Local variables

c     Begin program

      mgcf = mg_proc_coarsen_factor

      end function MG_get_proc_coarsen_factor

c     selectRank_MG_coarse
c     #####################################################################
      subroutine selectRank_MG_coarse(split_key)

c     ---------------------------------------------------------------------
c     Groups processors in slices COARSE MG domains.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(INOUT) :: split_key

c     Local variables

        integer :: npx_new,npy_new,npz_new,ip,jp,kp,npx,npy,npz

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Set global limits of local domain l_lim(dim,loc)

        l_lim(1,0) = g_mg_def%ilo(1)
        l_lim(2,0) = g_mg_def%jlo(1)
        l_lim(3,0) = g_mg_def%klo(1)
        l_lim(1,1) = g_mg_def%ihi(1)
        l_lim(2,1) = g_mg_def%jhi(1)
        l_lim(3,1) = g_mg_def%khi(1)

c     Set global limits of global domain g_lim(dim,loc)

        g_lim(1,0) = 1
        g_lim(2,0) = 1
        g_lim(3,0) = 1
        g_lim(1,1) = g_mg_def%nxgl(1)
        g_lim(2,1) = g_mg_def%nygl(1)
        g_lim(3,1) = g_mg_def%nzgl(1)

c     Find coarse domain config

        npx = g_mg_def%mx
        npy = g_mg_def%my
        npz = g_mg_def%mz
        
        !New processor count
        npx_new = max(npx/mg_proc_coarsen_factor,1)
        npy_new = max(npy/mg_proc_coarsen_factor,1)
        npz_new = max(npz/mg_proc_coarsen_factor,1)

        !Coordinates in new processor mesh
        ip = (npx_new*l_lim(1,0)/g_lim(1,1))
        jp = (npy_new*l_lim(2,0)/g_lim(2,1))
        kp = (npz_new*l_lim(3,0)/g_lim(3,1))

c$$$        !At bottom of hierarchy separate all procs (SELF)
c$$$        if (npx_new == 1.and.npx /= 1) ip = my_rank
c$$$        if (npy_new == 1.and.npy /= 1) jp = my_rank
c$$$        if (npz_new == 1.and.npz /= 1) kp = my_rank

c     Select rank

        !Group processors based on location in x, y and z directions
        !Organize in 3D COARSE processor mesh: i-1 + ny*(j-1) + ny*nz*(k-1)
        split_key = split_key + ip  + npy_new*jp
     .                              + npy_new*npz_new*kp

        !At bottom of hierarchy consider SELF
c$$$        if (npx_new == 1.and.npy_new == 1.and.npz_new == 1) then
c$$$          split_key = my_rank
c$$$        endif

c$$$        write(*,*) 'Proc=',my_rank,'mpi_comm_crse X'
c$$$     .             ,npx_new,npx,l_lim(1,0),g_lim(1,1),ip,split_key

c     End program

      end subroutine selectRank_MG_coarse

c     selectRank_MG_split
c     #####################################################################
      subroutine selectRank_MG_split(split_key)

c     ---------------------------------------------------------------------
c     Groups processors in checkerboard pattern over coarse MG domains.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer,INTENT(INOUT) :: split_key

c     Local variables

        integer :: my_rank_l,np_l

c     Begin program

        if (split_key == MPI_UNDEFINED) then
          split_key = 1
        endif

c     Select rank

        !Group processors based on RANK in COARSE MG communicator
        call MPI_Comm_size(MPI_COMM_SAV,np_l,mpierr)
        call MPI_Comm_rank(MPI_COMM_SAV,my_rank_l,mpierr)

        split_key = split_key + my_rank_l

        !At bottom of hierarchy consider SELF
        if (np_l == 1) split_key = my_rank

c$$$        write(*,*) 'Proc=',my_rank,'mpi_comm_split X',split_key

c     End program

      end subroutine selectRank_MG_split

c     create_MG_comms
c     #################################################################
      function create_MG_comms(g_def,mpi_comm_splt,mpi_comm_crse)
     .         result(end_of_hierarchy)

c     -----------------------------------------------------------------
c     Creates recursive coarse-proc MG MPI communicators
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        logical :: end_of_hierarchy

        integer,INTENT(OUT) :: mpi_comm_splt,mpi_comm_crse

        type(grid_mg_def),pointer,INTENT(IN)  :: g_def

c     Local variables

        integer :: npx,npy,npz,np_l

c     Begin program

        g_mg_def => g_def

        end_of_hierarchy=.false.

        !Create CRSE comm (for gather within domains)
        call splitMPIComm  (mpi_comm_crse,selectRank_MG_coarse
     .                     ,MPI_COMM_ORIG=g_mg_def%mpi_comm)
        if (coarse_mg_debug) then
          call check_MPI_COMM(mpi_comm_crse,'MPI_COMM_MG_CRSE')
        endif

        MPI_COMM_SAV = mpi_comm_crse

        !Create SPLIT comm (for concurrent solve across domains)
        call splitMPIComm  (mpi_comm_splt,selectRank_MG_split
     .                     ,MPI_COMM_ORIG=g_mg_def%mpi_comm)
        if (coarse_mg_debug) then
          call check_MPI_COMM(mpi_comm_splt,'MPI_COMM_MG_SPLT')
        endif

        !Check for termination (when no more splitting is possible)
        call MPI_Comm_size(MPI_COMM_SPLT,np_l,mpierr)

        end_of_hierarchy = (np_l == 1)

c     End program

      end function create_MG_comms
#endif

      end module grid_mg

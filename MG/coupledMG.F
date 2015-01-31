c TODO list:
c
c 1) Implement possibility of forming matrix from matvec routine
c    for very expensive matvec's. This requires implementing a
c    proxy matvec routine that selects the matrix-free or regular
c    matrix-vector multiplies. As the diagonal, there should be
c    the possibility of forming this matrix within the suite 
c    or taking it from the outside.
c
c 2) In weighed Jacobi, replace omega factors by a matrix-vector
c    operation.
c
c 3) Document i/o variables in all subroutines consistently.

c module mg_ops
c######################################################################
      module mg_ops

        use grid_mg

        use bc_ops

        use grid_create, ONLY: grid_pack,g_pack,createMGgrid

        use mlsolverSetup

        integer, parameter :: BOX=1,STAR=2

      contains

c     getMGvcomp
c     #################################################################
      function getMGvcomp(mgctx,igr,igptr,i,j,k,ieq) result(ijkg)

c     -----------------------------------------------------------------
c     Gives MG vector component corresponing to equation ieq on
c     coordinates (i,j,k) and grid level igr.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(mg_ctx),pointer :: mgctx
        integer :: i,j,k,ieq,igr,ijkg,igptr

c     Local variables

        integer :: nx,ny,nz

c     Begin program

        nx = mgctx%g_def%nxv(igr)
        ny = mgctx%g_def%nyv(igr)
        nz = mgctx%g_def%nzv(igr)

        ijkg = mgctx%neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
     .       + mgctx%istart(igptr) - 1

      end function getMGvcomp

c     MGvecadd
c     ###################################################################
      subroutine MGvecadd(mgctx,igr,coef1,vec1,coef2,vec2)

c     -------------------------------------------------------------------
c     Performs the vector add operation vec1 <- coef1*vec1 + coef2*vec2,
c     but restricted on a grid patch at grid level igr. The grid patch is
c     determined by a call to the "limits" routine. In call:
c       * igr: grid level
c       * neq: number of variables contained in vectors.
c       * ntot: vector dimension
c       * coef1,vec1: on input, first term of sum. On output, vec1 contains
c           sum result.
c       * coef2,vec2: second term of sum.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(mg_ctx),pointer :: mgctx
      integer :: igr
      real(8) :: coef1,coef2,vec1(:),vec2(:)

c     Local variables

      integer :: i,j,k,iii,ieq,nx,ny,nz

c     Begin program

cc      nnx = MGgrid%nxv(igr)
cc      nny = MGgrid%nyv(igr)
cc      nnz = MGgrid%nzv(igr)
cc
cc      call limits(0,nnx,nny,nnz,igr,imin,imax,jmin,jmax,kmin,kmax)

      nx = mgctx%g_def%nxv(igr)
      ny = mgctx%g_def%nyv(igr)
      nz = mgctx%g_def%nzv(igr)

      call MGlimits(0,igr,mgctx)

      do k=mgctx%kmin,mgctx%kmax
        do j=mgctx%jmin,mgctx%jmax
          do i=mgctx%imin,mgctx%imax
            do ieq=1,mgctx%neq
              iii=mgctx%neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
              vec1(iii) = coef1*vec1(iii) + coef2*vec2(iii)
            enddo
          enddo
        enddo
      enddo

      end subroutine MGvecadd

c     MGdot
c     ###################################################################
      function MGdot(mgctx,igr,vec1,vec2) result(dot)

c     -------------------------------------------------------------------
c     Performs scalar product (vec1,vec2),but restricted on a grid patch
c     at grid level igr. The grid patch is determined by a call to the
c     "limits" routine. In call:
c       * igr: grid level
c       * neq: number of variables contained in vectors.
c       * ntot: vector dimension
c       * vec1: vector, first term of scalar product.
c       * vec2: vector, second term of scalar product.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(mg_ctx),pointer :: mgctx
      integer    :: igr
      real(8)    :: vec1(:),vec2(:),dot

c     Local variables

      integer    :: i,j,k,iii,ieq,nx,ny,nz
      real(8)    :: ldot

c     Begin program

cc      nnx = MGgrid%nxv(igr)
cc      nny = MGgrid%nyv(igr)
cc      nnz = MGgrid%nzv(igr)
cc
cc      call limits(0,nnx,nny,nnz,igr,imin,imax,jmin,jmax,kmin,kmax)

cc      write (*,*) 'DIAG -- MGdot',nnx,nny,nnz
cc     .         ,imin,imax,jmin,jmax,kmin,kmax

      nx = mgctx%g_def%nxv(igr)
      ny = mgctx%g_def%nyv(igr)
      nz = mgctx%g_def%nzv(igr)

      call MGlimits(0,igr,mgctx)

      dot = 0d0

      do k=mgctx%kmin,mgctx%kmax
        do j=mgctx%jmin,mgctx%jmax
          do i=mgctx%imin,mgctx%imax
            do ieq=1,mgctx%neq
              iii=mgctx%neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
              dot = dot + vec1(iii)*vec2(iii)
            enddo
          enddo
        enddo
      enddo

#if defined(petsc)
      if (.not.asm) then
        ldot = dot
        call MPI_Allreduce(ldot,dot,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,mgctx%mpi_comm,mpierr)
      endif
#endif

      end function MGdot

c     scheduleColors
c     ###################################################################
      subroutine scheduleColors(mgctx,igrid,ncolors,nc1,nc2,nc3
     .                         ,i1min,i1max,j1min,j1max,k1min,k1max
     .                         ,stencil_width,stencil_type
     .                         ,zebra_dir)

c     -------------------------------------------------------------------
c     This routine schedules colors according to grid dimensions. In
c     call sequence:
c         * igrid(input): grid level
c         * ncolors (input/output): number of colors
c         * nc1,nc2,nc3 (output): number of colors in each grid dimension
c         * i1min,i1max,j1min,j1max,k1min,k1max(output): local grid limits
c              on which to perform colored relaxation.
c         * stencil_width (input, integer, optional): stencil width in each
c              dimension (default = 3)
c         * stencil_type (input, char, optional): stencil type (BOX=1, STAR=2;
c              default = box)
c         * zebra_dir (input, integer, optional): indicates direction of
c              zebra (line) relaxation.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: igrid,ncolors,nc1,nc2,nc3
     .          ,i1min,i1max,j1min,j1max,k1min,k1max

      integer,optional :: stencil_width,stencil_type,zebra_dir

c     Local variables
      
      integer :: nxf,nyf,nzf,nxl,nyl,nzl,nclrs
     .          ,imin,imax,jmin,jmax,kmin,kmax
     .          ,stncl_wdt,stncl_tpe,dd,ddd,zb_dir

c     Begin program

      if (PRESENT(stencil_width)) then
        stncl_wdt = stencil_width
      else
        stncl_wdt = 3
      endif

      if (PRESENT(stencil_type)) then
        stncl_tpe = stencil_type
      else
        stncl_tpe = BOX
      endif

      if (PRESENT(zebra_dir)) then
        zb_dir = zebra_dir
      else
        zb_dir = 0
      endif

      nclrs = stncl_wdt/2+1

cc      !Global grid size
cc      nxf = nxv(igrid)
cc      nyf = nyv(igrid)
cc      nzf = nzv(igrid)

      !Local grid size (different from global in line smoothing)
      nxl = mgctx%nxv(igrid)
      nyl = mgctx%nyv(igrid)
      nzl = mgctx%nzv(igrid)

      !Find grid limits
cc      call limits(0,nxf,nyf,nzf,igrid,imin,imax,jmin,jmax,kmin,kmax)
      call MGlimits(0,igrid,mgctx)

      imin = mgctx%imin
      jmin = mgctx%jmin
      kmin = mgctx%kmin

      imax = mgctx%imax
      jmax = mgctx%jmax
      kmax = mgctx%kmax

      if (nxl > 1 .and. nyl > 1 .and. nzl > 1) then

        dd = 3
        if (zb_dir /= 0) dd = dd-1

        select case(zb_dir)
        case(2)
          i1min = jmin
          i1max = jmax
          j1min = kmin
          j1max = kmax
          k1min = imin
          k1max = imax
        case(3)
          i1min = kmin
          i1max = kmax
          j1min = imin
          j1max = imax
          k1min = jmin
          k1max = jmax
        case default
          i1min = imin
          i1max = imax
          j1min = jmin
          j1max = jmax
          k1min = kmin
          k1max = kmax
        end select

      elseif (nxl > 1 .and. nyl > 1 .and. nzl == 1) then

        dd = 2
        if (zb_dir /= 0) dd = dd-1

        select case(zb_dir)
        case(2)
          i1min = jmin
          i1max = jmax
          j1min = imin
          j1max = imax
          k1min = kmin
          k1max = kmax
        case(3)
          nc1 = 0
          nc2 = 0
          nc3 = 0
          return
        case default
          i1min = imin
          i1max = imax
          j1min = jmin
          j1max = jmax
          k1min = kmin
          k1max = kmax
        end select

      elseif (nxl > 1 .and. nyl == 1 .and. nzl > 1) then

        dd = 2
        if (zb_dir /= 0) dd = dd-1

        select case(zb_dir)
        case(1)
          i1min = imin
          i1max = imax
          j1min = kmin
          j1max = kmax
          k1min = jmin
          k1max = jmax
        case(2)
          nc1 = 0
          nc2 = 0
          nc3 = 0
          return
        case default
          i1min = kmin
          i1max = kmax
          j1min = imin
          j1max = imax
          k1min = jmin
          k1max = jmax
        end select

      elseif (nxl == 1 .and. nyl > 1 .and. nzl > 1) then

        dd = 2
        if (zb_dir /= 0) dd = dd-1

        select case(zb_dir)
        case(1)
          nc1 = 0
          nc2 = 0
          nc3 = 0
          return
        case(3)
          i1min = kmin
          i1max = kmax
          j1min = jmin
          j1max = jmax
          k1min = imin
          k1max = imax
        case default
          i1min = jmin
          i1max = jmax
          j1min = kmin
          j1max = kmax
          k1min = imin
          k1max = imax
        end select

      elseif (nxl > 1 .and. nyl == 1 .and. nzl == 1) then

        dd = 1

        select case(zb_dir)
        case(2,3)
          nc1 = 0
          nc2 = 0
          nc3 = 0
          return
        case default
          i1min = imin
          i1max = imax
          j1min = jmin
          j1max = jmax
          k1min = kmin
          k1max = kmax
        end select

      elseif (nxl == 1 .and. nyl > 1 .and. nzl == 1) then

        dd = 1

        select case(zb_dir)
        case(1,3)
          nc1 = 0
          nc2 = 0
          nc3 = 0
          return
        case default
          i1min = jmin
          i1max = jmax
          j1min = kmin
          j1max = kmax
          k1min = imin
          k1max = imax
        end select

      elseif (nxl == 1 .and. nyl == 1 .and. nzl > 1) then

        dd = 1

        select case(zb_dir)
        case(1,2)
          nc1 = 0
          nc2 = 0
          nc3 = 0
          return
        case default
          i1min = kmin
          i1max = kmax
          j1min = imin
          j1max = imax
          k1min = jmin
          k1max = jmax
        end select

      endif

      nclrs = nclrs**dd

      !Account for special case of 2 colors in 3-pt star stencil
      if (  (ncolors == 2 .or. stncl_tpe == STAR)
     .     .and. stncl_wdt == 3 .and. zb_dir == 0) nclrs = 2
cc      nclrs = min(ncolors,nclrs)
cc      nclrs = ncolors  !Trust input

cc      write (*,*) dd
cc      pause

      !Allocate colors per dimension (in this order: nc1, nc2, nc3)

cc      select case(zb_dir)
cc      case(1)
cc        nc1 = 1
cc        ddd = max(dd,1)
cc        nc2 = max(1,nint((1d0*nclrs/nc1)**(1d0/ddd)))
cc        nc3 = nclrs/nc1/nc2
cc      case(2)
cc        nc2 = 1
cc        ddd = max(dd,1)
cc        nc1 = max(1,nint((1d0*nclrs/nc1)**(1d0/ddd)))
cc        nc3 = nclrs/nc1/nc2
cc      case(3)
cc        nc3 = 1
cc        ddd = max(dd,1)
cc        nc1 = max(1,nint((1d0*nclrs/nc1)**(1d0/ddd)))
cc        nc2 = nclrs/nc1/nc2
cc      case default
cc        nc1 = max(2,nint((1d0*nclrs)**(1d0/dd)))
cc        ddd = max(dd - 1,1)
cc        nc2 = max(1,nint((1d0*nclrs/nc1)**(1d0/ddd)))
cc        nc3 = nclrs/nc1/nc2
cc      end select

      if (zb_dir /= 0) then
        nc1 = 1
        ddd = max(dd,1)
        nc2 = max(1,nint((1d0*nclrs/nc1)**(1d0/ddd)))
        nc3 = nclrs/nc1/nc2
      else
        nc1 = max(2,nint((1d0*nclrs)**(1d0/dd)))
        ddd = max(dd - 1,1)
        nc2 = max(1,nint((1d0*nclrs/nc1)**(1d0/ddd)))
        nc3 = nclrs/nc1/nc2
      endif

      if (nc1*nc2*nc3 /= nclrs) then
        if (my_rank == 0) then
          write (*,*) 'Input ncolors=',ncolors
          write (*,*) 'Grid=',igrid
          write (*,*) 'Dimensions =',dd
          write (*,*) 'Colors in direction 1=',nc1
          write (*,*) 'Colors in direction 2=',nc2
          write (*,*) 'Colors in direction 3=',nc3
          write (*,*) 'Total number of colors=',nclrs
        endif
        call pstop('scheduleColors',
     .             'Unsupported number of colors')
      endif

      !Save number of colors
cc      ncolors = nclrs

cc      select case(nclrs)
cc      case(2)
cc        nc1 = 2
cc        nc2 = 1
cc        nc3 = 1
cc      case(3)
cc        nc1 = 3
cc        nc2 = 1
cc        nc3 = 1
cc      case(4)
cc        nc1 = 2
cc        nc2 = 2
cc        nc3 = 1
cc      case(8)
cc        nc1 = 2
cc        nc2 = 2
cc        nc3 = 2
cc      case(9)
cc        nc1 = 3
cc        nc2 = 3
cc        nc3 = 1
cc      case(27)
cc      case default
cc        call pstop('scheduleColors',
cc     .             'Unsupported number of colors')
cc      end select

c     End program

      end subroutine scheduleColors

c     scheduleIndices
c     ###################################################################
      subroutine scheduleIndices(mgctx,igrid,i1,j1,k1,i,j,k,zebra_dir)

c     -------------------------------------------------------------------
c     This routine schedules colored indices according to grid dimensions
c     and to the orderings in scheduleColors. In call:
c        * igrid: grid level
c        * i1,j1,k1: logical colored indices
c        * i,j,k: actual grid indices
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: igrid,i1,j1,k1,i,j,k
      integer,optional :: zebra_dir

c     Local variables

      integer :: nxl,nyl,nzl,zb_dir

c     Begin program

      if (PRESENT(zebra_dir)) then
        zb_dir = zebra_dir
      else
        zb_dir = 0
      endif

      !Local grid size (different than global in line smoothing)
      nxl = mgctx%nxv(igrid)
      nyl = mgctx%nyv(igrid)
      nzl = mgctx%nzv(igrid)

c     Swap indices

      if (nxl > 1 .and. nyl > 1 .and. nzl > 1) then

        select case(zb_dir)
        case(2)
          i = k1
          j = i1
          k = j1
        case(3)
          i = j1
          j = k1
          k = i1
        case default
          i = i1
          j = j1
          k = k1
        end select

      elseif (nxl > 1 .and. nyl > 1 .and. nzl == 1) then

        select case(zb_dir)
        case(2)
          i = j1
          j = i1
          k = k1
        case(3)
          call pstop('scheduleIndices','Wrong zebra direction')
        case default
          i = i1
          j = j1
          k = k1
        end select

      elseif (nxl > 1 .and. nyl == 1 .and. nzl > 1) then

        select case(zb_dir)
        case(1)
          i = i1
          j = k1
          k = j1
        case(2)
          call pstop('scheduleIndices','Wrong zebra direction')
        case default
          i = j1
          j = k1
          k = i1
        end select

      elseif (nxl == 1 .and. nyl > 1 .and. nzl > 1) then

        select case(zb_dir)
        case(1)
          call pstop('scheduleIndices','Wrong zebra direction')
        case(3)
          i = k1
          j = j1
          k = i1
        case default
          i = k1
          j = i1
          k = j1
        end select

      elseif (nxl > 1 .and. nyl == 1 .and. nzl == 1) then

        select case(zb_dir)
        case(2,3)
          call pstop('scheduleIndices','Wrong zebra direction')
        case default
          i = i1
          j = j1
          k = k1
        end select

      elseif (nxl == 1 .and. nyl > 1 .and. nzl == 1) then

        select case(zb_dir)
        case(1,3)
          call pstop('scheduleIndices','Wrong zebra direction')
        case default
          i = k1
          j = i1
          k = j1
        end select

      elseif (nxl == 1 .and. nyl == 1 .and. nzl > 1) then

        select case(zb_dir)
        case(1,2)
          call pstop('scheduleIndices','Wrong zebra direction')
        case default
          i = j1
          j = k1
          k = i1
        end select

      endif

      end subroutine scheduleIndices

      end module mg_ops

c module mg_Xfer
c ######################################################################
      module mg_Xfer

        use mg_ops

        use variables

        real(8),private :: eps=0d0

        type :: garray
          real(8),pointer,dimension(:,:,:,:) :: array
        end type garray

        type :: mg_array
          type(grid_mg_def),pointer :: g_def
          type(garray),pointer,dimension(:) :: grid
        end type mg_array

        logical :: is__cnv,is__vec,res_is_vec

        INTERFACE ASSIGNMENT (=)
          module procedure equateMGArray
        END INTERFACE

        INTERFACE fillArray
          module procedure fillArray1,fillArray2
        END INTERFACE

        INTERFACE fillArray_nobc
          module procedure fillArray1_nobc,fillArray2_nobc
     .                    ,fillArray3_nobc
        END INTERFACE

        INTERFACE fillVec
          module procedure fillVec1,fillVec2,fillVec3
        END INTERFACE

        private :: restrictArray,restrictArray_w_BCs!,mapMGVectorToArray

      contains

c     allocateMGArray
c     #################################################################
      subroutine allocateMGarray(neq,mgarray,g_def)

        implicit none

c     Call variables

        integer :: neq
        type(mg_array)  :: mgarray
        type(grid_mg_def),pointer,optional :: g_def

c     Local variables

        integer :: igrid,nxp,nyp,nzp
        type(grid_mg_def),pointer :: gl_def

c     Begin program

        if (PRESENT(g_def)) then
          gl_def => g_def
        else
          gl_def => grid_params
        endif
cc        if (.not.associated(MGgrid)) MGgrid => grid_params

        if (.not.associated(mgarray%grid)) then
          allocate(mgarray%grid(gl_def%ngrid))
          mgarray%g_def => gl_def
          do igrid=1,gl_def%ngrid
            nxp = gl_def%nxv(igrid)+1
            nyp = gl_def%nyv(igrid)+1
            nzp = gl_def%nzv(igrid)+1
            allocate(mgarray%grid(igrid)%array(0:nxp,0:nyp,0:nzp,neq))
            mgarray%grid(igrid)%array = 0d0
          enddo
        endif

c     End program

      end subroutine allocateMGarray

c     deallocateMGArray
c     #################################################################
      subroutine deallocateMGArray(mgarray)

        implicit none

c     Call variables

        type(mg_array)  :: mgarray

c     Local variables

        integer          :: igrid

c     Begin program

        if (associated(mgarray%grid)) then
          nullify(mgarray%g_def)
          do igrid=1,size(mgarray%grid)
            if (associated(mgarray%grid(igrid)%array)) then
              deallocate(mgarray%grid(igrid)%array)
            endif
          enddo
          deallocate(mgarray%grid)
        endif

c     End program

      end subroutine deallocateMGArray

c     equateMGArray
c     #################################################################
      subroutine equateMGArray(mg2,mg1)

        implicit none

c     Call variables

        type(mg_array),intent(IN)  :: mg1
        type(mg_array),intent(OUT) :: mg2

c     Local variables

        integer :: igrid,neq

c     Begin program

        neq = size(mg1%grid(1)%array,4)

        call allocateMGArray(neq,mg2)

        mg2%g_def => mg1%g_def

        do igrid=1,mg1%g_def%ngrid
          mg2%grid(igrid)%array = mg1%grid(igrid)%array
        enddo

c     End program

      end subroutine equateMGArray

c     mapArrayToMGVector
c     ##################################################################
      subroutine mapArrayToMGVector(mgctx,array,mgvector,igr,ismgv)
c     ------------------------------------------------------------------
c     Maps array into a MG vector excluding ghost nodes. In call sequence:
c       * mgctx: MG context
c       * array (in): contains variables in all grid nodes at grid level
c           igr.
c       * mgvector (out): array values stored in MG vector format, w/o
c           ghost cell values.
c       * igr (in): grid level where array is defined.
c     ------------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx
      integer :: igr
      real(8) :: mgvector(:)
     .          ,array(0:,0:,0:,:)
      logical :: ismgv

c     Local variables

      integer :: i,j,k,ii,ieq,neq,nx,ny,nz,igptr
cc      logical :: fpointers

c     Begin program

cc      call allocPointers(neq,fpointers)

      nx = size(array,1)-2
      ny = size(array,2)-2
      nz = size(array,3)-2
      neq = size(array,4)

      !Choose correct mapping for non-MG vectors
      if (ismgv) then
        igptr = igr
      else
        igptr = 1
      endif

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            do ieq=1,neq
              ii = getMGvcomp(mgctx,igr,igptr,i,j,k,ieq)
              mgvector(ii) = array(i,j,k,ieq)
            enddo
          enddo
        enddo
      enddo

cc      call deallocPointers(fpointers)

c     End program

      end subroutine mapArrayToMGVector

c     mapMGVectorToArray
c     ##################################################################
      subroutine mapMGVectorToArray(mgctx,gpos,mgvector,array,igr,ismgv)
c     ------------------------------------------------------------------
c     Maps values of MG vector mgvector corresponding to grid level igr
c     into array, w/o ghost cells. If vector is not a MG vector, set
c     igrid=1. In call sequence:
c       * mgctx: MG context
c       * gpos (in): if positive, map only locally in stencil of grid node
c           determined by gpos. Otherwise, map whole grid.
c       * mgvector (in): vector to be mapped
c       * array (out): mapped array
c       * igr (in): grid level
c       * ismgv (in): logical variable that determines whether mgvector
c           is indeed a MG vector or a regular vector.
c     ------------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: igr,gpos
      real(8) :: mgvector(:),array(0:,0:,0:,:)
      logical :: ismgv

c     Local variables

      integer :: neq,nx,ny,nz,ieq,i,j,k,ii,offset,igptr
      integer :: imin ,imax ,jmin ,jmax ,kmin ,kmax
     .          ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
      logical :: fpointers

c     Begin program

      call MGlimits(gpos,igr,mgctx)

      imin = mgctx%imin
      jmin = mgctx%jmin
      kmin = mgctx%kmin

      imax = mgctx%imax
      jmax = mgctx%jmax
      kmax = mgctx%kmax

      neq  = mgctx%neq

c     Define stencil

      nx = mgctx%g_def%nxv(igr)
      ny = mgctx%g_def%nyv(igr)
      nz = mgctx%g_def%nzv(igr)

      offset = 1

      iimin = max(imin-offset,1)
      iimax = min(imax+offset,nx)
      jjmin = max(jmin-offset,1)
      jjmax = min(jmax+offset,ny)
      kkmin = max(kmin-offset,1)
      kkmax = min(kmax+offset,nz)

c     Map vector to array

      !Choose correct mapping for non-MG vectors
      if (ismgv) then
        igptr = igr
      else
        igptr = 1
      endif

      array = 0d0

      do k = kkmin,kkmax
        do j = jjmin,jjmax
          do i = iimin,iimax
            do ieq=1,neq
              ii = getMGvcomp(mgctx,igr,igptr,i,j,k,ieq)
              array(i,j,k,ieq) = mgvector(ii)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapMGVectorToArray

c     restrictMGArray
c     #################################################################
      subroutine restrictMGArray(icmp,neq,mgarray,bcnd,igrid,order
     .                          ,order_bc,iscnv,isvec,result_is_vec
     .                          ,equilb,mgvcov)
c     -----------------------------------------------------------------
c     Restricts MG array in all grids with ghost nodes. If equilibrium
c     is not provided, then extrapolation is done at EQU boundaries.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer :: neq,icmp,bcnd(6,neq),order,igrid
      type(mg_array) :: mgarray

      integer,optional :: order_bc
      logical,optional,intent(IN) :: iscnv,isvec,result_is_vec
      type(mg_array),optional,intent(IN)  :: equilb
      type(mg_array),optional,intent(OUT) :: mgvcov

c     Local variables

      integer :: igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
      logical :: fpointers

c     Begin program

      if (PRESENT(isvec)) then
        is__vec = isvec
      else
        if (PRESENT(mgvcov)) then
          is__vec = .true.
        else
          is__vec = (neq == 3 .or. neq == 6)
        endif
      endif

      if (PRESENT(iscnv)) then
        is__cnv = iscnv
      else
        is__cnv = is__vec   !Contravariant representation by default
      endif

      if (PRESENT(result_is_vec)) then
        res_is_vec = result_is_vec
      else
        res_is_vec = is__vec   !Contravariant representation by default
      endif

c     Consistency check

      if (size(mgarray%grid(igrid)%array,4) /= neq) then
        write (*,*) 'Cannot restrict MG array: ',
     .              'inconsistent number of components'
        write (*,*) neq,size(mgarray%grid)
        write (*,*) 'Aborting...'
        stop
      endif

c     Restrict array

      do igc=igrid+1,mgarray%g_def%ngrid
        igf = igc-1

        nxf = mgarray%g_def%nxv(igf)
        nyf = mgarray%g_def%nyv(igf)
        nzf = mgarray%g_def%nzv(igf)
        nxc = mgarray%g_def%nxv(igc)
        nyc = mgarray%g_def%nyv(igc)
        nzc = mgarray%g_def%nzv(igc)

        call restrictArray_w_BCs(icmp,neq
     .       ,igf,nxf,nyf,nzf,mgarray%grid(igf)%array
     .       ,igc,nxc,nyc,nzc,mgarray%grid(igc)%array
     .       ,order,.false.,bcnd
     .       ,order_bc=order_bc
     .       ,equilb=equilb
     .       ,mgvcov=mgvcov
     .       ,g_def=mgarray%g_def)

      enddo

c     End program

      end subroutine restrictMGArray

c     interpArray_setup
c     #################################################################
      subroutine interpArray_setup(x1,y1,z1,array1
     .                            ,tx,ty,tz,work,bcoef)
c     -----------------------------------------------------------------
c     This is a restriction routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * order (input): order of interpolation
c       * x1,y1,z1 (input): source interpolation nodes (1D)
c       * array1 (input): source data (3D)
c       * tx,ty,tz,work,bcoef (out): spline arrays
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      real(8) :: x1(:),y1(:),z1(:),array1(:,:,:)
      real(8) :: tx(:),ty(:),tz(:),work(:),bcoef(:,:,:)

c     Local variables

      integer :: nx,ny,nz,kx,ky,kz,flg
      real(8),allocatable,dimension(:) :: q

c     Begin program

      flg = 0

      nx = size(bcoef,1)
      ny = size(bcoef,2)
      nz = size(bcoef,3)

      kx = size(tx)-nx
      ky = size(ty)-ny
      kz = size(tz)-nz

      if (nx == 1 .and. ny == 1) then
        allocate(q((2*kz-1)*nz))
        call dbknot(z1,nz,kz,tz)
        call dbintk(z1,array1(1,1,:),tz,nz,kz,bcoef(1,1,:),q,work)
        deallocate(q)
      elseif (nx == 1 .and. nz == 1) then
        allocate(q((2*ky-1)*ny))
        call dbknot(y1,ny,ky,ty)
        call dbintk(y1,array1(1,:,1),ty,ny,ky,bcoef(1,:,1),q,work)
        deallocate(q)
      elseif (ny == 1 .and. nz == 1) then
        allocate(q((2*kx-1)*nx))
        call dbknot(x1,nx,kx,tx)
        call dbintk(x1,array1(:,1,1),tx,nx,kx,bcoef(:,1,1),q,work)
        deallocate(q)
      elseif (nx == 1) then
        call db2ink(y1,ny,z1,nz,array1(1,:,:)
     .             ,ny,ky,kz,ty,tz,bcoef(1,:,:),work,flg)
      elseif (ny == 1) then
        call db2ink(x1,nx,z1,nz,array1(:,1,:)
     .             ,nx,kx,kz,tx,tz,bcoef(:,1,:),work,flg)
      elseif (nz == 1) then
        call db2ink(x1,nx,y1,ny,array1(:,:,1)
     .             ,nx,kx,ky,tx,ty,bcoef(:,:,1),work,flg)
      else
        call db3ink(x1,nx,y1,ny,z1,nz,array1
     .             ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)
      endif

c     End program

      end subroutine interpArray_setup

c     interpArray_spline
c     #################################################################
      subroutine interpArray_spline(x1,y1,z1,tx,ty,tz,work,bcoef
     .                             ,x2,y2,z2,array2,d_x,d_y,d_z)
c     -----------------------------------------------------------------
c     This is a restriction routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * tx,ty,tz,work,bcoef (input): spline arrays
c       * x1,y1,z1 (input): source interpolation nodes (1D)
c       * x2,y2,z2 (output): target interpolation nodes (1D)
c       * array2 (output): target interpolated data (3D)
c       * d_x,d_y,d_z (optional input): order of derivatives in X,Y,Z
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      real(8) :: tx(:),ty(:),tz(:),work(:)
      real(8) :: bcoef(:,:,:)

      real(8) :: x1(:),y1(:),z1(:),x2(:),y2(:),z2(:),array2(:,:,:)

      integer, optional :: d_x,d_y,d_z

c     Local variables

      integer :: i,j,k,nx,ny,nz,nx2,ny2,nz2

      real(8) :: xx,yy,zz,val

      integer :: kx,ky,kz,inbv,derx,dery,derz

      real(8) :: dbvalu,db2val,db3val
      external   dbvalu,db2val,db3val

c     Begin program

      if (PRESENT(d_x)) then
        derx = d_x
      else
        derx = 0
      endif

      if (PRESENT(d_y)) then
        dery = d_y
      else
        dery = 0
      endif

      if (PRESENT(d_z)) then
        derz = d_z
      else
        derz = 0
      endif

c     Compute storage based on source data

      nx = size(bcoef,1)
      ny = size(bcoef,2)
      nz = size(bcoef,3)

      kx = size(tx)-nx
      ky = size(ty)-ny
      kz = size(tz)-nz

c     Perform spline interpolation

      nx2 = size(array2,1)
      ny2 = size(array2,2)
      nz2 = size(array2,3)

      if (eps == 0d0) eps = 1d4*findRoundOff()

      do k = 1,nz2
        do j = 1,ny2
          do i = 1,nx2

            xx = x2(i)
            yy = y2(j)
            zz = z2(k)

            !Do nothing if outside of original domain (otherwise --> error in spliners)
            if (  xx < x1(1) - eps
     .       .or. xx > x1(nx)+ eps
     .       .or. yy < y1(1) - eps
     .       .or. yy > y1(ny)+ eps
     .       .or. zz < z1(1) - eps
     .       .or. zz > z1(nz)+ eps ) then
              cycle
            endif

            if (nx == 1 .and. ny == 1) then
              inbv=1
              val = dbvalu(tz,bcoef(1,1,:),nz,kz,derz,zz,inbv,work)
            elseif (nx == 1 .and. nz == 1) then
              inbv=1
              val = dbvalu(ty,bcoef(1,:,1),ny,ky,dery,yy,inbv,work)
            elseif (ny == 1 .and. nz == 1) then
              inbv=1
              val = dbvalu(tx,bcoef(:,1,1),nx,kx,derx,xx,inbv,work)
            elseif (nx == 1) then
              val=db2val(yy,zz,dery,derz,ty,tz,ny,nz,ky,kz,bcoef(1,:,:)
     .                  ,work)
            elseif (ny == 1) then
              val=db2val(xx,zz,derx,derz,tx,tz,nx,nz,kx,kz,bcoef(:,1,:)
     .                  ,work)
            elseif (nz == 1) then
              val=db2val(xx,yy,derx,dery,tx,ty,nx,ny,kx,ky,bcoef(:,:,1)
     .                  ,work)
            else
              val=db3val(xx,yy,zz,derx,dery,derz,tx,ty,tz,nx,ny,nz
     .                      ,kx,ky,kz,bcoef,work)
            endif

            array2(i,j,k) = val

          enddo
        enddo
      enddo

c     End program

      end subroutine interpArray_spline

c     interpArray
c     #################################################################
      subroutine interpArray(order,x1,y1,z1,array1,x2,y2,z2,array2)
c     -----------------------------------------------------------------
c     This is a restriction routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * order (input): order of interpolation
c       * x1,y1,z1 (input): source interpolation nodes (1D)
c       * array1 (input): source data (3D)
c       * x2,y2,z2 (input): target interpolation nodes (1D)
c       * array2 (output): target interpolated data (3D)
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer :: order

      real(8) :: x1(:),y1(:),z1(:),x2(:),y2(:),z2(:)

      real(8) :: array1(:,:,:),array2(:,:,:)

c     Local variables

      integer :: nx,ny,nz

      integer :: kx,ky,kz,dim
      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

c     Begin program

c     Compute storage based on source data

      nx = size(array1,1)
      ny = size(array1,2)
      nz = size(array1,3)

      kx = min(order+1,nx-1)
      ky = min(order+1,ny-1)
      kz = min(order+1,nz-1)

      dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

      allocate(tx(nx+kx))
      allocate(ty(ny+ky))
      allocate(tz(nz+kz))
      allocate(work(dim))
      allocate(bcoef(nx,ny,nz))

c     Setup up spline coeffs working array

      call interpArray_setup(x1,y1,z1,array1,tx,ty,tz,work,bcoef)

c     Perform spline interpolation

      call interpArray_spline(x1,y1,z1,tx,ty,tz,work,bcoef
     .                       ,x2,y2,z2,array2)

      deallocate(tx,ty,tz,work,bcoef)

c     End program

      end subroutine interpArray

c     restrictArray_w_BCs
c     #################################################################
      subroutine restrictArray_w_BCs(icmp,neq,igf,nxf,nyf,nzf,arrayf
     .                                       ,igc,nxc,nyc,nzc,arrayc
     .                              ,order,volf,bcnd,order_bc,equilb
     $                              ,mgvcov,g_def)
c     -----------------------------------------------------------------
c     Restricts array to array in all grids (with ghost nodes),
c     starting at grid igf.
c     -----------------------------------------------------------------

      use setMGBC_interface

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
     .             ,order,bcnd(6,neq),icmp
      real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
     .             ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)
      logical    :: volf

      integer,optional :: order_bc
      type(mg_array),optional,intent(IN)  :: equilb
      type(mg_array),optional,intent(OUT) :: mgvcov
      type(grid_mg_def),pointer,optional :: g_def

c     Local variables

      integer    :: igridf,igridc,ieq,ierr,ord_bc
     .             ,nxxf,nyyf,nzzf,nxxc,nyyc,nzzc
     .             ,bcmod(6,neq),icmp2(max(neq/3,1))

      real(8),allocatable,dimension(:,:,:,:) :: arrf,arrc
      logical    :: fpointers,have_equl,have_cov

c     Begin program

      call allocPointers(neq,fpointers,g_def=g_def)

      have_equl = PRESENT(equilb)
      have_cov = PRESENT(mgvcov)

      if (PRESENT(order_bc)) then
        ord_bc = order_bc
      else
        ord_bc = order
      endif

c     Restrict MG vectors

cc      allocate(arrf(0:nxf+1,0:nyf+1,0:nzf+1,neq))
      allocate(arrc(0:nxf+1,0:nyf+1,0:nzf+1,neq))

      arrc = arrayf

      do igridc = igf+1,igc

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxxf = MGgrid%nxv(igridf)
        nyyf = MGgrid%nyv(igridf)
        nzzf = MGgrid%nzv(igridf)

        nxxc = MGgrid%nxv(igridc)
        nyyc = MGgrid%nyv(igridc)
        nzzc = MGgrid%nzv(igridc)

c       Allocate coarse mesh vector

cc        deallocate(arrf)
        allocate(arrf(0:nxxf+1,0:nyyf+1,0:nzzf+1,neq))

        arrf = arrc

        deallocate(arrc)
        allocate(arrc(0:nxxc+1,0:nyyc+1,0:nzzc+1,neq))

c       Restrict vector

        do ieq = 1,neq
          call restrictArray(mg_ctx_temp
     .                      ,arrf(:,:,:,ieq)
     .                      ,arrc(:,:,:,ieq)
     .                      ,order,igridf,volf)
        enddo

        deallocate(arrf)

      enddo

      arrayc = arrc

      deallocate(arrc)
cc      deallocate(arrf,arrc)

c     Impose BCs

      if (icmp /= 0) then
        icmp2 = icmp
        if (have_equl) then
          if (have_cov) then
            call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
     .                  ,arr0=equilb%grid(igc)%array
     .                  ,arr_cov=mgvcov%grid(igc)%array
     .                  ,icomp=icmp2,is_cnv=is__cnv,is_vec=is__vec
     .                  ,result_is_curv=res_is_vec,iorder=ord_bc
     .                  ,g_def=g_def)
          else
            call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
     .                  ,arr0=equilb%grid(igc)%array
     .                  ,icomp=icmp2,is_cnv=is__cnv,is_vec=is__vec
     .                  ,result_is_curv=res_is_vec,iorder=ord_bc
     .                  ,g_def=g_def)
          endif
        else
          bcmod = bcnd
          where (abs(bcmod) == EQU) bcmod = EXT
          if (have_cov) then
            call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcmod
     .                  ,arr_cov=mgvcov%grid(igc)%array
     .                  ,icomp=icmp2,is_cnv=is__cnv,is_vec=is__vec
     .                  ,result_is_curv=res_is_vec,iorder=ord_bc
     .                  ,g_def=g_def)
          else
            call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcmod
     .                  ,icomp=icmp2,is_cnv=is__cnv,is_vec=is__vec
     .                  ,result_is_curv=res_is_vec,iorder=ord_bc
     .                  ,g_def=g_def)
          endif
        endif
      endif

c     Deallocate vectors

      call deallocPointers(fpointers)

c     End program

      end subroutine restrictArray_w_BCs

c     restrictArray
c     #################################################################
      subroutine restrictArray(mgctx,arrayf,arrayc,order,igf,volw)
c     -----------------------------------------------------------------
c     This is a restriction routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * arrayc (real): array in coarse grid
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * arrayf (real): array in fine grid
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igf (int): fine grid level identifier
c       * volw (logical): whether values are volume-weighed.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: order,igf

      real(8) :: arrayf(0:,0:,0:)
     .          ,arrayc(0:,0:,0:)

      logical :: volw

c     Local variables
 
      integer :: nxc,nyc,nzc,nxf,nyf,nzf,ic,if,jc,jf,kc,kf,igc
     .          ,icg,jcg,kcg,ifg,jfg,kfg
      integer :: iminc,imaxc,jminc,jmaxc,kminc,kmaxc
     .          ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf

      real(8) :: volt,vol,xc

      logical :: fpointers

c     Interpolation

      real(8) :: ff

      real(8),allocatable,dimension(:) :: xx,yy,zz,xxc,yyc,zzc

c     Begin program

      nxc = size(arrayc,1)-2
      nyc = size(arrayc,2)-2
      nzc = size(arrayc,3)-2

      nxf = size(arrayf,1)-2
      nyf = size(arrayf,2)-2
      nzf = size(arrayf,3)-2

      call allocPointers(1,fpointers,mgctx=mgctx)

      igc = igf + 1

      arrayc = 0d0

c     Agglomeration

      if (order.eq.0) then

c$$$        call limits(0,nxc,nyc,nzc,igc
c$$$     .             ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)
        call MGlimits(0,igc,mgctx)

        iminc = mgctx%imin
        jminc = mgctx%jmin
        kminc = mgctx%kmin

        imaxc = mgctx%imax
        jmaxc = mgctx%jmax
        kmaxc = mgctx%kmax

        do kc = kminc,kmaxc
          do jc = jminc,jmaxc
            do ic = iminc,imaxc

              volt   = 0d0

              do kf = mg_ratio_z(igf)*(kc-1)+1,mg_ratio_z(igf)*kc
                do jf = mg_ratio_y(igf)*(jc-1)+1,mg_ratio_y(igf)*jc
                  do if = mg_ratio_x(igf)*(ic-1)+1,mg_ratio_x(igf)*ic

                    if (.not.volw) then
                      vol =mgctx%g_def%gmetric%grid(igf)%dlvol(if,jf,kf)  !Logical volume
                      arrayc(ic,jc,kc) = arrayc(ic,jc,kc)
     .                                 + arrayf(if,jf,kf)*vol
                      volt = volt + vol
                    else
                      arrayc(ic,jc,kc) = arrayc(ic,jc,kc)
     .                                 + arrayf(if,jf,kf)
                    endif

                  enddo
                enddo
              enddo

              if (.not.volw) arrayc(ic,jc,kc) = arrayc(ic,jc,kc)
     .                    /mgctx%g_def%gmetric%grid(igc)%dlvol(ic,jc,kc)

            enddo
          enddo
        enddo

c       Interpolation

      else

c       Define interpolation nodes

c$$$        call limits(0,nxf,nyf,nzf,igf
c$$$     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf)
        call MGlimits(0,igf,mgctx)

        iminf = mgctx%imin
        jminf = mgctx%jmin
        kminf = mgctx%kmin

        imaxf = mgctx%imax
        jmaxf = mgctx%jmax
        kmaxf = mgctx%kmax

        allocate(xx(1:imaxf-iminf+1)
     .          ,yy(1:jmaxf-jminf+1)
     .          ,zz(1:kmaxf-kminf+1))

        call getMGmap(iminf,jminf,kminf,igf,igf,igf,ifg,jfg,kfg
     .               ,g_def=mgctx%g_def)

        xx = mgctx%g_def%xx(ifg:ifg+imaxf-iminf)
        yy = mgctx%g_def%yy(jfg:jfg+jmaxf-jminf)
        zz = mgctx%g_def%zz(kfg:kfg+kmaxf-kminf)

        call limits(0,nxc,nyc,nzc,igc
     .             ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

        allocate(xxc(1:imaxc-iminc+1)
     .          ,yyc(1:jmaxc-jminc+1)
     .          ,zzc(1:kmaxc-kminc+1))

        call getMGmap(iminc,jminc,kminc,igc,igc,igc,icg,jcg,kcg
     .               ,g_def=mgctx%g_def)

        xxc = mgctx%g_def%xx(icg:icg+imaxc-iminc)
        yyc = mgctx%g_def%yy(jcg:jcg+jmaxc-jminc)
        zzc = mgctx%g_def%zz(kcg:kcg+kmaxc-kminc)

c       Take out volume factors

        if (volw) then
          do kf = kminf,kmaxf
            do jf = jminf,jmaxf
              do if = iminf,imaxf
                call getMGmap(if,jf,kf,igf,igf,igf,ifg,jfg,kfg
     .                       ,g_def=mgctx%g_def)

                arrayf(if,jf,kf) = arrayf(if,jf,kf)
     .                    /mgctx%g_def%gmetric%grid(igf)%dlvol(if,jf,kf)
              enddo
            enddo
          enddo
        endif

c       Calculate interpolation

        call interpArray(order,xx, yy ,zz
     .       ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .                        ,xxc,yyc,zzc
     .       ,arrayc(iminc:imaxc,jminc:jmaxc,kminc:kmaxc))

c       Apply volume factors

        if (volw) then
          do kc = kminc,kmaxc
            do jc = jminc,jmaxc
              do ic = iminc,imaxc
                arrayc(ic,jc,kc) = arrayc(ic,jc,kc)
     .                    *mgctx%g_def%gmetric%grid(igc)%dlvol(ic,jc,kc)
              enddo
            enddo
          enddo
        endif

        deallocate(xx,yy,zz,xxc,yyc,zzc)

      endif

c     End program

      call deallocPointers(fpointers,mgctx=mgctx)

      end subroutine restrictArray

c     restrictArrayToMGVector
c     #################################################################
      subroutine restrictArrayToMGVector(array,mgvector,igr0,order,volf)
c     -----------------------------------------------------------------
c     Restricts array to mgvector in all grids (without ghost nodes),
c     starting at grid igr0.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: order,igr0
      real(8)    :: array(0:,0:,0:,:)
      real(8)    :: mgvector(:)
      logical    :: volf

c     Local variables

      integer    :: neq,nx,ny,nz,ieq,nxf,nyf,nzf,nxc,nyc,nzc
     .             ,igridc,igridf,isigf,isigc,ntotc,ntotf
      logical    :: fpointers

c     Begin program

      nx  = size(array,1)-2
      ny  = size(array,2)-2
      nz  = size(array,3)-2
      neq = size(array,4)

      call allocPointers(neq,fpointers)

c     Consistency check

      nxf = mg_ctx_temp%g_def%nxv(igr0)
      nyf = mg_ctx_temp%g_def%nyv(igr0)
      nzf = mg_ctx_temp%g_def%nzv(igr0)

      if (nxf /= nx .or. nyf /= ny .or. nzf /= nz) then
        write (*,*) 'Grid mismatch in restrictArray:'
        write (*,*) 'Aborting...'
        stop
      endif

c     Map array in initial grid onto MG vector

      call mapArrayToMGVector(mg_ctx_temp,array,mgvector,igr0,.true.)

c     Restrict array to coarser grids

      do igridc = igr0+1,mg_ctx_temp%ngrid

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxf = mg_ctx_temp%g_def%nxv(igridf)
        nyf = mg_ctx_temp%g_def%nyv(igridf)
        nzf = mg_ctx_temp%g_def%nzv(igridf)

        nxc = mg_ctx_temp%g_def%nxv(igridc)
        nyc = mg_ctx_temp%g_def%nyv(igridc)
        nzc = mg_ctx_temp%g_def%nzv(igridc)

        isigc = mg_ctx_temp%istart(igridc)
        isigf = mg_ctx_temp%istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Restrict MG vector

        call crestrict(mg_ctx_temp,mgvector(isigf:isigf+ntotf-1)
     .                      ,mgvector(isigc:isigc+ntotc-1)
     .                ,order,igridf,volf)
      enddo

      call deallocPointers(fpointers)

      end subroutine restrictArrayToMGVector

c     restrictMGVector
c     #################################################################
      subroutine restrictMGVector(mgctx,mgvector,igr0,igr1,order,volf)
c     -----------------------------------------------------------------
c     Restricts vector to mgvector from grid igr0 to grid igr1.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer    :: order,igr0,igr1
      real(8)    :: mgvector(:)
      logical    :: volf

c     Local variables

      integer    :: neq,ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf
     .             ,isigf,isigc,ntotc,ntotf
cc      logical    :: fpointers

c     Begin program

cc      call allocPointers(neq,fpointers)

c     Consistency check

      nxf = mgctx%g_def%nxv(igr0)
      nyf = mgctx%g_def%nyv(igr0)
      nzf = mgctx%g_def%nzv(igr0)

      neq = mgctx%neq

      if (igr0 > igr1) then
        call pstop('restrictVectorToMGVector'
     .            ,'Starting grid is coarser than final grid')
      endif

c     Restrict vector to coarser grids

      do igridc = igr0+1,igr1

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxf = mgctx%g_def%nxv(igridf)
        nyf = mgctx%g_def%nyv(igridf)
        nzf = mgctx%g_def%nzv(igridf)

        nxc = mgctx%g_def%nxv(igridc)
        nyc = mgctx%g_def%nyv(igridc)
        nzc = mgctx%g_def%nzv(igridc)

        isigc = mgctx%istart(igridc)
        isigf = mgctx%istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Restrict MG vector

        call crestrict(mgctx,mgvector(isigf:isigf+ntotf-1)
     .                      ,mgvector(isigc:isigc+ntotc-1)
     .                ,order,igridf,volf)
      enddo

cc      call deallocPointers(fpointers)

      end subroutine restrictMGVector

c     prolongMGVector
c     #################################################################
      subroutine prolongMGVector(mgctx,mgvector,igr0,igr1,order,bcnd)
c     -----------------------------------------------------------------
c     Prolongs vector to mgvector from grid igr0 to grid igr1.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer    :: order,igr0,igr1,bcnd(:,:)
      real(8)    :: mgvector(:)

c     Local variables

      integer    :: neq,ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf
     .             ,isigf,isigc,ntotc,ntotf
cc      logical    :: fpointers

c     Begin program

cc      call allocPointers(neq,fpointers)

c     Consistency check

      nxc = mgctx%g_def%nxv(igr0)
      nyc = mgctx%g_def%nyv(igr0)
      nzc = mgctx%g_def%nzv(igr0)

      neq = mgctx%neq

      if (igr0 < igr1) then
        call pstop('prolongMGVector'
     .            ,'Starting grid is coarser than final grid')
      endif

c     Prolong vector to coarser grids

      do igridf = igr0-1,igr1,-1

        igridc = igridf+1

c       Characterize coarse and fine grids

        nxf = mgctx%g_def%nxv(igridf)
        nyf = mgctx%g_def%nyv(igridf)
        nzf = mgctx%g_def%nzv(igridf)

        nxc = mgctx%g_def%nxv(igridc)
        nyc = mgctx%g_def%nyv(igridc)
        nzc = mgctx%g_def%nzv(igridc)

        isigc = mgctx%istart(igridc)
        isigf = mgctx%istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Prolong MG vector

        call cprolong(mgctx,mgvector(isigc:isigc+ntotc-1)
     .                     ,mgvector(isigf:isigf+ntotf-1)
     .               ,order,igridc,bcnd)
      enddo

cc      call deallocPointers(fpointers,mgctx=mgctx)

      end subroutine prolongMGVector

c     restrictDerivedType
c     #################################################################
      subroutine restrictDerivedType(igf,vf,igc,vc,order)
c     -----------------------------------------------------------------
c     Restricts array to array in all grids (with ghost nodes),
c     starting at grid igf.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer :: order,igf,igc

      type(var_array),pointer :: vf
      type(var_array),pointer :: vc

c     Local variables

      integer :: nxf,nyf,nzf,nxc,nyc,nzc,ieq

      real(8),pointer,dimension(:,:,:) :: arrayf,arrayc
      integer,pointer,dimension(:) :: bcnd

      logical :: fpointers,have_equl

c     Begin program

c     Find grid sizes

      nxf = size(vf%array_var(1)%array,1)-2
      nyf = size(vf%array_var(1)%array,2)-2
      nzf = size(vf%array_var(1)%array,3)-2

      nxc = size(vc%array_var(1)%array,1)-2
      nyc = size(vc%array_var(1)%array,2)-2
      nzc = size(vc%array_var(1)%array,3)-2
      
c     Restrict

      do ieq=1,vf%nvar
        arrayf => vf%array_var(ieq)%array
        arrayc => vc%array_var(ieq)%array
        bcnd   => vf%array_var(ieq)%bconds
        call restrictArray_w_BCs(ieq,1,igf,nxf,nyf,nzf,arrayf
     .                                ,igc,nxc,nyc,nzc,arrayc
     .                          ,order,.false.,bcnd)
      enddo

      nullify(arrayf,arrayc,bcnd)

c     End program

      end subroutine restrictDerivedType

c     prolongDerivedType
c     #################################################################
      subroutine prolongDerivedType(mgctx,igc,vc,igf,vf,order)
c     -----------------------------------------------------------------
c     Prolongs array to array in all grids (with ghost nodes),
c     starting at grid igf.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: order,igf,igc

      type(var_array),pointer :: vf
      type(var_array),pointer :: vc

c     Local variables

      integer :: ieq

      real(8),pointer,dimension(:,:,:) :: arrayf,arrayc
      integer,pointer,dimension(:) :: bcnd

      logical :: fpointers,have_equl

c     Begin program

c     Prolong

      do ieq=1,vf%nvar
        arrayf => vf%array_var(ieq)%array
        arrayc => vc%array_var(ieq)%array
        call prolongArray(mgctx,igc,arrayc,igf,arrayf,order)
      enddo

      nullify(arrayf,arrayc,bcnd)

c     End program

      end subroutine prolongDerivedType

c     crestrict
c     #################################################################
      subroutine crestrict(mgctx,xf,xc,order,igf,volw)
c     -----------------------------------------------------------------
c     This is a restriction routine for system MG, with arbitrary order
c     of interpolation.
c
c     In call sequence, we have:
c       * mgctx: MG context
c       * xf (real): vector in fine grid
c       * xc (real): vector in coarse grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igf (int): fine grid level identifier
c       * volw (logical): whether vectors contain volume fractions.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: order,igf

      real(8) :: xc(:),xf(:)

      logical :: volw

c     Local variables
 
      real(8),allocatable,dimension(:,:,:,:) :: arrayf,arrayc

      integer :: nxc,nyc,nzc,nxf,nyf,nzf,ieq,igc,neq

c     Begin program

      igc = igf + 1

      nxf = mgctx%g_def%nxv(igf)
      nyf = mgctx%g_def%nyv(igf)
      nzf = mgctx%g_def%nzv(igf)

      nxc = mgctx%g_def%nxv(igc)
      nyc = mgctx%g_def%nyv(igc)
      nzc = mgctx%g_def%nzv(igc)

      neq = mgctx%neq

      allocate(arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
     .        ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq))

c     Map vector to array

      call mapMGVectorToArray(mgctx,0,xf,arrayf,igf,.false.)

c     Restrict array(s)

      arrayc = 0d0

      do ieq = 1,mgctx%neq
        call restrictArray(mgctx
     .                    ,arrayf(:,:,:,ieq)
     .                    ,arrayc(:,:,:,ieq)
     .                    ,order,igf,volw)
      enddo

c     Map arrays to vector

      call mapArrayToMGVector(mgctx,arrayc,xc,igc,.false.)

c     End program

      deallocate(arrayc,arrayf)

      end subroutine crestrict

c     cprolong
c     #################################################################
      subroutine cprolong(mgctx,xc,xf,order,igc,bcnd)
c     -----------------------------------------------------------------
c     This is a prolongation routine for system MG, with arbitrary order
c     of interpolation.
c
c     In call sequence, we have:
c       * mgctx: MG context
c       * xc (real): vector in coarse grid (not a MG vector)
c       * xf (real): vector in fine grid (not a MG vector)
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igc (int): coarse grid level identifier
c       * bcnd (int array): boundary condition info.
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: order,igc,bcnd(:,:)
      real(8) :: xc(:),xf(:)

c     Local variables
 
      integer :: neq,nxc,nyc,nzc,nxf,nyf,nzf,igf

      real(8),allocatable,dimension(:,:,:,:) :: arrayf,arrayc

c     Begin program

      igf = igc - 1

      nxf = mgctx%g_def%nxv(igf)
      nyf = mgctx%g_def%nyv(igf)
      nzf = mgctx%g_def%nzv(igf)

      nxc = mgctx%g_def%nxv(igc)
      nyc = mgctx%g_def%nyv(igc)
      nzc = mgctx%g_def%nzv(igc)

      neq = mgctx%neq

      allocate(arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
     .        ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq))

c     Unpack vector into array

      call mapMGVectorToArray(mgctx,0,xc,arrayc,igc,.false.)

c     Prolong array

      call cprolongArray

c     Repack prolonged vector

      call mapArrayToMGVector(mgctx,arrayf,xf,igf,.false.)

c     End program

      deallocate(arrayc,arrayf)

      contains

c     cprolongArray
c     #################################################################
      subroutine cprolongArray

      use setMGBC_interface

      implicit none            ! For safe Fortran

c     Call variables

c     Local variables
 
      real(8),allocatable :: sp_sum(:,:),dv_sum(:,:),l_sum(:,:)

      integer :: ic,jc,if,jf,iic,iif,i,j,k,ig,jg,kg,ieq

      logical :: cartesian,isvec,iscnv

c     Begin program

c     Impose boundary conditions (external)

      isvec = (neq == 3)                !Initialize isvec (may be modified by setMGBC)

      cartesian = (bcSP() .and. isvec)  !Return Cartesian variables for prolongation

      if (order > 0 .or. cartesian) then
        call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
cc     .              ,iorder=min(order,3)             !now inherited from previous calls
     .              ,result_is_curv=(.not.cartesian)
     .              ,out_iscurv=isvec,out_iscnv=iscnv  !<= These define whether quantities are 
     .              ,g_def=mgctx%g_def)                !   vectors or scalars. 
      endif

c     Perform average along azimuthal direction for SP configurations

      if (isSP2(1,igc) .and. isvec) then

        allocate(sp_sum(0:nzc+1,neq),dv_sum(0:nzc+1,neq))

        do ieq=1,neq
          do k=0,nzc+1
            sp_sum(k,ieq)= sum(arrayc(1,1:nyc,k,ieq)
     .                 *mgctx%g_def%gmetric%grid(igc)%dlvol(1,1:nyc,k))
            dv_sum(k,ieq)=
     .              sum(mgctx%g_def%gmetric%grid(igc)%dlvol(1,1:nyc,k))
          enddo
        enddo

#if defined(petsc)
        allocate(l_sum(0:nzc+1,neq))
        l_sum = sp_sum
        call MPI_Allreduce(l_sum,sp_sum,neq*(nzc+2)
     .                    ,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)

        l_sum = dv_sum
        call MPI_Allreduce(l_sum,dv_sum,neq*(nzc+2)
     .                    ,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)

        deallocate(l_sum)
#endif

        do ieq=1,neq
          do k=0,nzc+1
            arrayc(1,:,k,ieq) = sp_sum(k,ieq)/dv_sum(k,ieq)
            arrayc(0,:,k,ieq) = sp_sum(k,ieq)/dv_sum(k,ieq)
          enddo
        enddo

        deallocate(sp_sum,dv_sum)
      endif

c     Prolong arrays (use scalar prolongation)

      arrayf = 0d0

      do ieq=1,neq
        call prolongArray(mgctx,igc,arrayc(:,:,:,ieq)
     .                         ,igf,arrayf(:,:,:,ieq)
     .                   ,order)
      enddo

c     Transform back to curvilinear (if isvec=.true. and cartesian=.true.)

      if (isvec .and. cartesian) then
        call XformToCurv_mesh2(igf,arrayf,iscnv)
      endif

c     End program

      end subroutine cprolongArray

      end subroutine cprolong

c     prolongArray
c     #################################################################
      subroutine prolongArray(mgctx,igc,arrayc,igf,arrayf,order)
c     -----------------------------------------------------------------
c     This is a prolongation routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * arrayf (real): array in fine grid
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * arrayc (real): array of values in coarse grid
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igc (int): coarse grid level identifier
c       * igf (int): fine grid level identifier
c     -----------------------------------------------------------------

      implicit none            ! For safe Fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: igc,igf,order
      real(8) :: arrayc(0:,0:,0:),arrayf(0:,0:,0:)

c     Local variables

      integer :: ic,if,jc,jf,kc,kf 
     .          ,icg,jcg,kcg,ifg,jfg,kfg
      integer :: iminc,imaxc,jminc,jmaxc,kminc,kmaxc
     .          ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf
      real(8) :: mag,vol

      logical :: fpointers

      real(8),allocatable,dimension(:) :: xx,yy,zz,xxf,yyf,zzf

c     Begin program

      call allocPointers(1,fpointers,mgctx=mgctx)

cc      call limits(0,nxc,nyc,nzc,igc,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

      call MGlimits(0,igc,mgctx) 

      iminc = mgctx%imin
      jminc = mgctx%jmin
      kminc = mgctx%kmin

      imaxc = mgctx%imax
      jmaxc = mgctx%jmax
      kmaxc = mgctx%kmax

c     Injection

      if (order.eq.0) then

        if (igf /= igc - 1) then
           call pstop('prolongArray'
     $               ,'Cannot perform injection: igf /= igc -1')
        endif

        do kc = kminc,kmaxc
          do jc = jminc,jmaxc
            do ic = iminc,imaxc

              do kf = mg_ratio_z(igf)*(kc-1)+1,mg_ratio_z(igf)*kc
                do jf = mg_ratio_y(igf)*(jc-1)+1,mg_ratio_y(igf)*jc
                  do if = mg_ratio_x(igf)*(ic-1)+1,mg_ratio_x(igf)*ic
                    arrayf(if,jf,kf) = arrayc(ic,jc,kc)
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo

      else

c Interpolation

c     Setup dimension vectors

        allocate(xx(1:imaxc-iminc+3)
     .          ,yy(1:jmaxc-jminc+3)
     .          ,zz(1:kmaxc-kminc+3))

        call getMGmap(iminc,jminc,kminc,igc,igc,igc,icg,jcg,kcg
     .               ,g_def=mgctx%g_def)

        xx = mgctx%g_def%xx(icg-1:icg+imaxc-iminc+1)
        yy = mgctx%g_def%yy(jcg-1:jcg+jmaxc-jminc+1)
        zz = mgctx%g_def%zz(kcg-1:kcg+kmaxc-kminc+1)

cc        call limits(0,nxf,nyf,nzf,igf
cc     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf)
        call MGlimits(0,igf,mgctx) 

        iminf = mgctx%imin
        jminf = mgctx%jmin
        kminf = mgctx%kmin

        imaxf = mgctx%imax
        jmaxf = mgctx%jmax
        kmaxf = mgctx%kmax

        allocate(xxf(1:imaxf-iminf+3)
     .          ,yyf(1:jmaxf-jminf+3)
     .          ,zzf(1:kmaxf-kminf+3))

        call getMGmap(iminf,jminf,kminf,igf,igf,igf,ifg,jfg,kfg
     .               ,g_def=mgctx%g_def)

        xxf = mgctx%g_def%xx(ifg-1:ifg+imaxf-iminf+1)
        yyf = mgctx%g_def%yy(jfg-1:jfg+jmaxf-jminf+1)
        zzf = mgctx%g_def%zz(kfg-1:kfg+kmaxf-kminf+1)

c     Calculate interpolation

        call interpArray(order,xx, yy ,zz ,arrayc
     .                        ,xxf,yyf,zzf,arrayf)

c     Free work arrays

        deallocate(xx,yy,zz,xxf,yyf,zzf)

      endif

c End program

      call deallocPointers(fpointers,mgctx=mgctx)

      end subroutine prolongArray

c     fillArray1
c     ###################################################################
      subroutine fillArray1(igrid,neq,x,arr,bcs,iorder,arr_cov,bc
     .                    ,arr0,is_vec,is_cnv,gpos,icomp,result_is_curv)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        use setMGBC_interface

        implicit none

c     Call variables

        integer :: igrid,neq,iorder,bcs(:,:)
        real(8),dimension(:,:) :: x
        real(8),dimension(:,:,:,:) :: arr

        integer,optional :: gpos,icomp(:)
        real(8),optional,dimension(:,:,:,:) :: arr_cov,arr0
        logical,optional :: is_vec,is_cnv,result_is_curv,bc

c     Local variables

        integer :: ieq,nx,ny,nz,gps
        logical :: isvec,iscnv,apply_bc
        logical :: fpointers

c     Begin program

        call allocPointers(1,fpointers)

        nx = MGgrid%nxv(igrid)
        ny = MGgrid%nyv(igrid)
        nz = MGgrid%nzv(igrid)

        if (PRESENT(is_vec)) then
          isvec = is_vec
        else
          isvec = (neq == 3 .or. neq == 6)
        endif

        if (PRESENT(is_cnv)) then
          iscnv = is_cnv
        else
          iscnv = .true.
        endif

        if (PRESENT(gpos)) then
          gps = max(0,gpos)
        else
          gps = 0
        endif

        if (PRESENT(bc)) then
          apply_bc = bc
        else
          apply_bc = .true.
        endif

        if (PRESENT(arr_cov).and.(.not.iscnv)) then
          do ieq=1,neq
            call mapMGVectorToArray(mg_ctx_temp,gps,x(:,ieq)
     .                           ,arr_cov(:,:,:,ieq:ieq),igrid,.false.)
          enddo
        else
          do ieq=1,neq
            call mapMGVectorToArray(mg_ctx_temp,gps,x(:,ieq)
     .                           ,arr(:,:,:,ieq:ieq),igrid,.false.)
          enddo
        endif

        if (PRESENT(arr0).and.PRESENT(arr_cov)) then
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr_cov=arr_cov,arr0=arr0
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        elseif (PRESENT(arr_cov)) then
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr_cov=arr_cov
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        elseif (PRESENT(arr0)) then
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr0=arr0
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        else
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        endif

      call deallocPointers(fpointers)

c     End program

      end subroutine fillArray1

c     fillArray2
c     ###################################################################
      subroutine fillArray2(igrid,neq,x,arr,bcs,iorder,arr_cov,arr0
     .                     ,is_vec,is_cnv,gpos,icomp,result_is_curv)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        use setMGBC_interface

        implicit none

c     Call variables

        integer :: igrid,neq,iorder,bcs(:,:)
        real(8),dimension(:) :: x
        real(8),dimension(:,:,:,:) :: arr

        integer,optional :: gpos,icomp(:)
        real(8),optional,dimension(:,:,:,:) :: arr_cov,arr0
        logical,optional :: is_vec,is_cnv,result_is_curv

c     Local variables

        integer :: ieq,nx,ny,nz,gps
        logical :: isvec,iscnv
        logical :: fpointers

c     Begin program

        call allocPointers(neq,fpointers)

        nx = MGgrid%nxv(igrid)
        ny = MGgrid%nyv(igrid)
        nz = MGgrid%nzv(igrid)

        if (PRESENT(is_vec)) then
          isvec = is_vec
        else
          isvec = (neq == 3 .or. neq == 6)
        endif

        if (PRESENT(is_cnv)) then
          iscnv = is_cnv
        else
          iscnv = .true.
        endif

        if (PRESENT(gpos)) then
          gps = max(0,gpos)
        else
          gps = 0
        endif

        if (PRESENT(arr_cov).and.(.not.iscnv)) then
          call mapMGVectorToArray(mg_ctx_temp,gps,x,arr_cov,igrid
     .                           ,.false.)
        else
          call mapMGVectorToArray(mg_ctx_temp,gps,x,arr,igrid,.false.)
        endif

        if (PRESENT(arr0).and.PRESENT(arr_cov)) then
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr_cov=arr_cov,arr0=arr0
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        elseif (PRESENT(arr_cov)) then
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr_cov=arr_cov
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        elseif (PRESENT(arr0)) then
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder,arr0=arr0
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        else
          call setMGBC(gps,neq,nx,ny,nz,igrid,arr,bcs
     .              ,icomp=icomp,is_vec=isvec,is_cnv=iscnv
     .              ,iorder=iorder
     .              ,result_is_curv=result_is_curv
     $              ,g_def=MGgrid)
        endif

        call deallocPointers(fpointers)

c     End program

      end subroutine fillArray2

c     fillArray1_nobc
c     ###################################################################
      subroutine fillArray1_nobc(igrid,neq,x,arr,gpos)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,gpos,neq
        real(8),dimension(:,:) :: x
        real(8),dimension(:,:,:,:) :: arr

c     Local variables

        integer :: ieq
        logical :: fpointers

c     Begin program

        call allocPointers(1,fpointers)

        do ieq=1,neq
          call mapMGVectorToArray(mg_ctx_temp,gpos,x(:,ieq)
     .                           ,arr(:,:,:,ieq:ieq),igrid,.false.)
        enddo

        call deallocPointers(fpointers)

c     End program

      end subroutine fillArray1_nobc

c     fillArray2_nobc
c     ###################################################################
      subroutine fillArray2_nobc(igrid,neq,x,arr,gpos)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,gpos,neq
        real(8),dimension(:) :: x
        real(8),dimension(:,:,:,:) :: arr

c     Local variables

        integer :: ieq

        logical :: fpointers

c     Begin program

        call allocPointers(neq,fpointers)

        call mapMGVectorToArray(mg_ctx_temp,gpos,x,arr,igrid,.false.)

        call deallocPointers(fpointers)

c     End program

      end subroutine fillArray2_nobc

c     fillArray3_nobc
c     ###################################################################
      subroutine fillArray3_nobc(igrid,x,arr,gpos)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,gpos,neq
        real(8),dimension(:) :: x
        real(8),dimension(:,:,:) :: arr

c     Local variables

        logical :: fpointers
        real(8),dimension(size(arr,1),size(arr,2),size(arr,3),1) :: arr1

c     Begin program

        call allocPointers(1,fpointers)

        call mapMGVectorToArray(mg_ctx_temp,gpos,x,arr1,igrid,.false.)
        arr = arr1(:,:,:,1)

        call deallocPointers(fpointers)

c     End program

      end subroutine fillArray3_nobc

c     fillVec1
c     ###################################################################
      subroutine fillVec1(igrid,neq,arr,x)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,neq
        real(8),dimension(:) :: x
        real(8),dimension(:,:,:,:) :: arr

c     Local variables

        logical :: fpointers

c     Begin program

        call allocPointers(neq,fpointers)

        call mapArrayToMGVector(mg_ctx_temp,arr,x,igrid,.false.)

        call deallocPointers(fpointers)

c     End program

      end subroutine fillVec1

c     fillVec2
c     ###################################################################
      subroutine fillVec2(igrid,neq,arr,x)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,neq
        real(8),dimension(:,:) :: x
        real(8),dimension(:,:,:,:) :: arr

c     Local variables

        integer :: ieq
        logical :: fpointers

c     Begin program

        call allocPointers(1,fpointers)

        do ieq=1,neq
          call mapArrayToMGVector(mg_ctx_temp,arr(:,:,:,ieq:ieq)
     .                           ,x(:,ieq),igrid,.false.)
        enddo

        call deallocPointers(fpointers)

c     End program

      end subroutine fillVec2

c     fillVec3
c     ###################################################################
      subroutine fillVec3(igrid,arr,x)

c     -------------------------------------------------------------------
c     Fills array arr (with BCs) using vector x.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,neq
        real(8),dimension(:) :: x
        real(8),dimension(:,:,:) :: arr

c     Local variables

        logical :: fpointers
        real(8),dimension(size(arr,1),size(arr,2),size(arr,3),1) :: arr1

c     Begin program

        call allocPointers(1,fpointers)

        arr1(:,:,:,1) = arr
        call mapArrayToMGVector(mg_ctx_temp,arr1,x,igrid,.false.)

        call deallocPointers(fpointers)

c     End program

      end subroutine fillVec3

      end module mg_Xfer

c module mg_io
c######################################################################
      module mg_io

        use mg_Xfer

        use grid_debug

        use math

      contains

c     MGplot
c     #################################################################
      subroutine MGplot(mgctx,igr,mgv,iflag,ofile)

c     -----------------------------------------------------------------
c     Plots MG vector using xdraw. In call:
c       * neq: number of variables contained in MG vector
c       * mgv: MG vector
c       * igr: grid level to be plotted.
c       * iflag: whether to initialize output file (iflag=0) or not.
c       * ofile: output file.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(mg_ctx),pointer :: mgctx
        integer    :: iflag,igr
        real(8)    :: mgv(:)
        character*(*) :: ofile

c     Local variables

        integer    :: neq,nx,ny,nz,ieq,nunit
        real(8),allocatable,dimension(:,:,:,:) :: debug

c     Begin program

        nunit = find_unit(110)

        neq= mgctx%neq

        nx = mgctx%g_def%nxv(igr)
        ny = mgctx%g_def%nyv(igr)
        nz = mgctx%g_def%nzv(igr)

c     Dump plots

        allocate(debug(0:nx+1,0:ny+1,0:nz+1,neq))

        if (iflag == 0) then
          open(nunit,file=trim(ofile),form='unformatted'
     .        ,status='replace')
        else
          open(nunit,file=trim(ofile),form='unformatted'
     .        ,status='old',position='append')
        endif

        call mapMGVectorToArray(mgctx,0,mgv,debug,igr,.true.)

        do ieq=1,neq
          call C_contour(debug(1:nx,1:ny,1:nz,ieq)
     .                ,mgctx%g_def%lxmin,mgctx%g_def%lxmax
     .                ,mgctx%g_def%lymin,mgctx%g_def%lymax
     .                ,mgctx%g_def%lzmin,mgctx%g_def%lzmax
     .                ,iflag+ieq-1,nunit)
        enddo

c     End program

        deallocate(debug)
        close(nunit)

      end subroutine MGplot

c     vec_plot
c     #####################################################################
      subroutine vec_plot(nvar,vec,igrid,file,title)

c     ---------------------------------------------------------------------
c     Plots vector with nvar variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,nvar
        real(8) :: vec(:)
        character(*) :: file,title

c     Local variables

        character(4)  :: dbg_desc(nvar)
        character(50) :: binfile,drawfile

        integer :: i,ieq
        logical :: fpointers

        type(mg_ctx),pointer:: mgctx => null()

c     Begin program

        binfile =trim(file)//'_p'//trim(int2char(my_rank))//'.bin'
        drawfile='draw'
     .         //trim(file)//'_p'//trim(int2char(my_rank))//'.in'

        do i=1,nvar
          dbg_desc(i) = trim(file)//trim(int2char(i))
        enddo

        call createDrawInCfile(nvar,trim(binfile),trim(title)
     .          ,'t','x','y',dbg_desc,'-c -X0 -L57'
     .          ,trim(drawfile))

        call allocPointers(nvar,fpointers,mgctx=mgctx)

        call MGplot(mgctx,igrid,vec,0,binfile)

        call deallocPointers(fpointers,mgctx=mgctx)

c     End program

      end subroutine vec_plot

c     MGplot_gl
c     #################################################################
      subroutine MGplot_gl(mgctx,igr,mgv,iflag,ofile)

c     -----------------------------------------------------------------
c     Plots MG vector using xdraw. In call:
c       * neq: number of variables contained in MG vector
c       * mgv: MG vector
c       * igr: grid level to be plotted.
c       * iflag: whether to initialize output file (iflag=0) or not.
c       * ofile: output file.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(mg_ctx),pointer :: mgctx
        integer :: iflag,igr
        real(8) :: mgv(:)
        character*(*) :: ofile

c     Local variables

        integer :: nx,ny,nz,nxg,nyg,nzg,ieq,nunit,neq
        real(8),allocatable,dimension(:,:,:,:) :: debugl
        real(8),allocatable,dimension(:,:,:)   :: debugg

c     Begin program

        nunit = find_unit(110)

        neq= mgctx%neq

        nx = mgctx%g_def%nxv(igr)
        ny = mgctx%g_def%nyv(igr)
        nz = mgctx%g_def%nzv(igr)

        nxg = mgctx%g_def%nxgl(igr)
        nyg = mgctx%g_def%nygl(igr)
        nzg = mgctx%g_def%nzgl(igr)

c     Dump plots

        allocate(debugg(0:nxg+1,0:nyg+1,0:nzg+1)
     .          ,debugl(0:nx +1,0:ny +1,0:nz +1,neq))

        call mapMGVectorToArray(mgctx,0,mgv,debugl,igr,.true.)

        if (my_rank == 0) then

          if (iflag == 0) then
            open(nunit,file=trim(ofile),form='unformatted'
     .          ,status='replace')
          else
            open(nunit,file=trim(ofile),form='unformatted'
     .          ,status='old',position='append')
          endif

        endif

        do ieq=1,neq
          call find_global_nobc(debugl(1:nx ,1:ny ,1:nz ,ieq)
     .                         ,debugg(1:nxg,1:nyg,1:nzg))

          if (my_rank == 0) then  
            call C_contour(debugg(1:nxg,1:nyg,1:nzg)
     .                  ,mgctx%g_def%gxmin,mgctx%g_def%gxmax
     .                  ,mgctx%g_def%gymin,mgctx%g_def%gymax
     .                  ,mgctx%g_def%gzmin,mgctx%g_def%gzmax
     .                  ,iflag+ieq-1,nunit)
          endif
        enddo

c     End program

        deallocate(debugg,debugl)
        close(nunit)

      end subroutine MGplot_gl

c     vec_plot_gl
c     #####################################################################
      subroutine vec_plot_gl(nvar,vec,igrid,file,title)

c     ---------------------------------------------------------------------
c     Plots vector with nvar variables.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igrid,nvar
        real(8) :: vec(:)
        character(*) :: file,title

c     Local variables

        character(4)  :: dbg_desc(nvar)
        character(50) :: binfile,drawfile

        integer :: i,ieq
        logical :: fpointers

        type(mg_ctx),pointer:: mgctx => null()

c     Begin program

        binfile =trim(file)//'.bin'
        drawfile='draw'//trim(file)//'.in'

        do i=1,nvar
          dbg_desc(i) = trim(file)//trim(int2char(i))
        enddo

        if (my_rank == 0) then
          call createDrawInCfile(nvar,trim(binfile),trim(title)
     .          ,'t','x','y',dbg_desc,'-c -X0 -L57'
     .          ,trim(drawfile))
        endif

        call allocPointers(nvar,fpointers,mgctx=mgctx)

        call MGplot_gl(mgctx,igrid,vec,0,binfile)

        call deallocPointers(fpointers,mgctx=mgctx)

cc        call MGplot_gl(nvar,vec,igrid,0,binfile)

c     End program

      end subroutine vec_plot_gl

      end module mg_io

c     module mg_solver
c     ##################################################################
      module mg_solver

        use mg_io, ONLY: MGplot
        use mg_Xfer

        logical :: vbr_mg
     .            ,galerkin=.false.

        private :: find_mf_diag_std,find_mf_diag_colored

      contains

c     blockSolve
c     #################################################################
      subroutine blockSolve(size,mat,icol,rhs,x)

c     -----------------------------------------------------------------
c     Solves block systems using a direct solve approach. Requires
c     linking with LAPACK.
c
c     In the call sequence, we have:
c       * size (integer): block size
c       * mat  (real array): block matrix
c       * icol (integer) : number of columns of rhs
c       * rhs  (real array): rhs of equation
c       * x    ( "     "   ): solution (output)
c       * smoother (character): identifies calling subroutine (JB, GS)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: size,icol

        real(8)    :: mat(size,size),rhs(size,icol),x(size,icol)

c     Local variables

        integer    :: ipiv(size),info

        real(8)    :: mat2(size,size)

        external dgesv

c     Begin program

        mat2 = mat !Avoid overwritting mat
        x = rhs

        call dgesv(size,icol,mat2,size,ipiv,x,size,info) !LAPACK routine

        call LAPACK_error(info)

      end subroutine blockSolve

c     blockInv
c     #################################################################
      subroutine blockInv(size,mat)

c     -----------------------------------------------------------------
c     Inverts block systems using a direct solve approach. Requires
c     linking with LAPACK.
c
c     In the call sequence, we have:
c       * size (integer): block size
c       * mat  (real array): block matrix on input, inverse on output
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: size

        real(8)    :: mat(size,size)

c     Local variables

        integer    :: ipiv(size),info,lwork

        real(8),allocatable,dimension(:) :: work

        external   dgetri,dgetrf

c     Begin program

c     Find LU decomposition

        call dgetrf(size,size,mat,size,ipiv,info) !LAPACK routine

        call LAPACK_error(info)

c     Invert matrix

        !Workspace query
        lwork=-1
cc        lwork=size
        allocate(work(1))
        call dgetri(size,mat,size,ipiv,work,lwork,info) !LAPACK routine
        lwork=floor(work(1))
        deallocate(work)

        !Matrix inversion
        allocate(work(lwork))
        call dgetri(size,mat,size,ipiv,work,lwork,info) !LAPACK routine
        deallocate(work)

        call LAPACK_error(info)

      end subroutine blockInv

c     LAPACK_error
c     #################################################################
      subroutine LAPACK_error(info)

c     -----------------------------------------------------------------
c     Error routine for LAPACK calls
c     -----------------------------------------------------------------

        use math

        implicit none

c     Call variables

        integer    :: info

c     Local variables

c     Begin program

        if (info /= 0) then
          if (info < 0) then
            messg = 'Problem in factorization in argument '
     .             //int2char(-info)
            call pstop('blockSolve',messg)
          else
            call pstop('blockSolve','Matrix is singular')
          endif
        endif

c     End program

      end subroutine LAPACK_error

c     cSolver
c     ##################################################################
      subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess,matvec
     .                  ,vol_wgt
     .                  ,dg,cvrg_tst,tol,iters
     .                  ,gm_driver
     .                  ,cg_driver
     .                  ,ks_it
     .                  ,sing_matrix
     .                  ,mg_grid_def
     .                  ,mg_order_res
     .                  ,mg_order_prol
     .                  ,mg_vcyc
     .                  ,mg_cycle
     .                  ,mg_coarse_grid_size
     .                  ,mg_gm_coarse_solve
     .                  ,mg_smooth
     .                  ,mg_galerkin
     .                  ,mg_debug
     .                  ,sm_ncolors
     .                  ,sm_omega
     .                  ,sm_it
     .                  ,sm_line_relax
     .                  ,sm_zebra_relax
     .                  ,sm_zebra_omega
     .                  ,sm_zebra_it
     .                  ,sm_zebra_dir
     .                  )
c     ---------------------------------------------------------------
c     This subroutine solves a coupled system of neq equations. 
c     In call sequence:
c       * neq: number of coupled equations
c       * ntotp: number of mesh points
c       * b: rhs
c       * x: solution
c       * bcnd: boundary condition defs.
c       * igrid: MG grid level (igrid=1 is finest level)
c       * out: level of output information
c       * guess: whether a non-trivial initial guess is provided
c               (iguess=1) or not (iguess=0)
c       * matvec (external): matrix-vector product definition.
c       * vol_wgt: whether residuals are volume-weighed or not.
c
c     Optional variables:
c       * dg: matrix neq*neq diagonal block (for stationary its).
c       * cvrg_tst: whether we want to perform a convergence test
c       * tol: relative convergence tolerance
c       * iters (output): number of iterations
c       * gm_driver: whether we want GMRES as a driver
c       * cg_driver: whether we want CG    as a driver
c       * ks_it: # of krylov subspace iterations
c       * sing_mat: whether GMRES is dealing with a singular matrix or not
c
c       * mg_grid_def: MG grid definition structure
c       * mg_order_res: restriction order
c       * mg_order_prol: restriction order
c       * mg_smooth: determines smoothing (jb,gs,gm; default jb)
c       * mg_vcyc: maximum number of MG V-cycles (default 1)
c       * mg_cycle: number of mu_cycles in MG (default 1)
c       * mg_galerkin: whether to do Galerkin coarsening of operator
c                  (true) or rediscretization (default: false)
c       * mg_coarse_grid_size: minimum coarse mesh (roughly mg_ratio^mg_cgs per
c                 direction).
c       * mg_gm_coarse_solve: whether we do GMRES at the coarsest grid
c       * mg_debug: turn on debugging in MG
c
c       * sm_it: smoother iterations
c       * sm_omega: damping for JB smoother
c       * sm_ncolors: number of colors in grid (for GS).
c       * sm_line_relax: whether we want line relaxation
c       * sm_zebra_relax: whether we want zebra GS relaxation
c       * sm_zebra_omega: inner relaxation factor in ZEBRA smoothing (between directions)
c       * sm_zebra_it: inner ZEBRA iterations (per direction).
c       * sm_zebra_dir: ZEBRA preferred direction (1,2,3; do all if zero)
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        integer,intent(IN) :: neq,ntotp,igrid,bcnd(6,neq),out,guess
        real(8) :: x(ntotp,neq),b(ntotp,neq)
        logical :: vol_wgt

        real(8),optional,intent(IN) :: dg(neq,2*neq*ntotp)
     .                                ,sm_omega
     .                                ,sm_zebra_omega
     .                                ,tol
        integer,optional,intent(IN) :: mg_order_res
     .                                ,mg_order_prol
     .                                ,mg_vcyc,mg_cycle
     .                                ,mg_coarse_grid_size
     .                                ,sm_it
     .                                ,sm_zebra_it
     .                                ,sm_zebra_dir
     .                                ,sm_ncolors
     .                                ,ks_it
        logical,optional,intent(IN) :: sm_line_relax
     .                                ,sm_zebra_relax
     .                                ,mg_gm_coarse_solve
     .                                ,mg_galerkin
     .                                ,mg_debug
     .                                ,gm_driver
     .                                ,sing_matrix
     .                                ,cg_driver
     .                                ,cvrg_tst
        integer,optional,intent(OUT) :: iters

        character(2),optional :: mg_smooth

        type(grid_mg_def),optional,pointer :: mg_grid_def

        external   matvec

c     Local variables

        integer    :: ntot,iorder_res,iorder_prol,vcyc,mu_cyc,nclrs,smit
     .               ,mgmin,i,ieq,gmit,zit,zdir
        real(8)    :: xi(ntotp*neq),bi(ntotp*neq),omg,rtol,zomg
        real(8), target :: ddg(neq,2*neq*ntotp)
        logical    :: gm_crse_solve,cvrg_test,lrelax,zrelax
     .               ,gm_drvr,smat,cg_drvr,glkin,bc_chk

        character(2) :: smth

c     Begin program

        res_vol_wgt = vol_wgt          !Initialize volume-weighing variable

c     Check input integrity

        bc_chk = .false.

        do i=1,6
          do ieq=1,neq
            bc_chk = (bcnd(i,ieq) < 0 .or. bcnd(i,ieq) > BCTOT)
          enddo
        enddo

        if (bc_chk) call pstop('cSolver','Wrong BC input')

c     Process optional arguments

        if (PRESENT(dg)) ddg = dg

        !Solver options <<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if (PRESENT(gm_driver)) then
          gm_drvr = gm_driver
        else
          gm_drvr = .false.
        endif

        if (PRESENT(sing_matrix)) then
          smat = sing_matrix
        else
          smat = .false.
        endif

        if (PRESENT(cg_driver)) then
          cg_drvr = cg_driver
        else
          cg_drvr = .false.
        endif

        if(PRESENT(ks_it)) then
          gmit = ks_it
        else
          gmit = 10
        endif

        if (PRESENT(cvrg_tst)) then
          cvrg_test = cvrg_tst
        else
          cvrg_test = .false.
        endif

        if(PRESENT(tol)) then
          rtol = tol
        else
          rtol = 1d-1
        endif

        !Multigrid options <<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if(PRESENT(mg_order_res)) then
          iorder_res = mg_order_res
        else
          iorder_res = 2
        endif

        if(PRESENT(mg_order_prol)) then
          iorder_prol = mg_order_prol
        else
          iorder_prol = 2
        endif

        if(PRESENT(mg_vcyc)) then
          vcyc = mg_vcyc
        else
          vcyc = 1
        endif

        if(PRESENT(mg_cycle)) then
          mu_cyc = mg_cycle
        else
          mu_cyc = 1
        endif

        if(PRESENT(mg_smooth)) then
          smth = mg_smooth
        else
          smth = 'jb'
        endif

        if(PRESENT(mg_galerkin)) then
          glkin = mg_galerkin
        else
          glkin = .false.
        endif

        if (PRESENT(mg_coarse_grid_size)) then
          mgmin = mg_coarse_grid_size
        else
          mgmin = 2
        endif

        if (PRESENT(mg_gm_coarse_solve)) then
          gm_crse_solve = mg_gm_coarse_solve
        else
          gm_crse_solve = .false.
        endif

        !Smoother options <<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if (PRESENT(sm_line_relax)) then
          lrelax = sm_line_relax
        else
          lrelax = .false.
        endif

        if (PRESENT(sm_zebra_relax)) then
          zrelax = sm_zebra_relax
        else
          zrelax = .false.
        endif

        if (PRESENT(sm_zebra_omega)) then
          zomg = sm_zebra_omega
        else
          zomg = 0.8
        endif

        if (PRESENT(sm_zebra_it)) then
          zit = sm_zebra_it
        else
          zit = 3
        endif

        if (PRESENT(sm_zebra_dir)) then
          zdir = sm_zebra_dir
        else
          zdir = 0  !Do all directions
        endif

        if(PRESENT(sm_it)) then
          smit = sm_it
        else
          smit = 2
        endif

        if(PRESENT(sm_omega)) then
          omg = sm_omega
        else
          if (smth == 'gs') then
            omg = 1d0
          else
            omg = 0.7
          endif
        endif

        if (PRESENT(sm_ncolors)) then
          nclrs = sm_ncolors
        else
          nclrs = 8
        endif

        if (zrelax) then
          smth = 'gs'
          omg  = 1d0
        endif

c     Gather variables for coupled solve

        do i=1,ntotp
          do ieq=1,neq
            xi(neq*(i-1)+ieq) = x(i,ieq)
            bi(neq*(i-1)+ieq) = b(i,ieq)
          enddo
        enddo

        if (cvrg_test) then

          bi = 0d0

          call random_number(xi)
cc          xi = 1d0

          if (guess == 0) then
            call pstop('cSolver','guess must be =1 for cvrgnce test')
          endif

cc          write (*,*) 'DIAG -- cSolver: JB/GS solve'
cc
ccc         Initialize solver
cc
cc          call solverInit
cc
ccc         Setup solver
cc
cc          call solverOptionsInit
cc
cc          solverOptions%iter    = smit
cc          solverOptions%tol     = rtol
cc
cc          if (PRESENT(dg)) solverOptions%diag => ddg
cc
cc          select case (smth)
cc          case('gs')
cc            solverOptions%omega   = 1d0
cc            solverOptions%ncolors = nclrs
cc            solverOptions%mg_zebra_relax = zrelax
cc            solverOptions%mg_zebra_prefd = 0
cccc            if(bcSP()) solverOptions%mg_zebra_prefd = 1
cc            solverOptions%mg_zebra_it    = 2
cc            solverOptions%mg_zebra_omega = 0.8
cc            call assembleSolverHierarchy('gs')
cc          case('gm')
cc            solverOptions%stp_test = 1 
cc            solverOptions%krylov_subspace = solverOptions%iter 
cc            solverOptions%singular_matrix = smat
cc            call assembleSolverHierarchy('gm')
cc            call assembleSolverHierarchy('id')
cc          case default !JB
cc            solverOptions%omega   = omg
cc            call assembleSolverHierarchy('jb')
cc          end select
cc
ccc         Invoke solver
cc
cc          ntot=neq*ntotp
cc          call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,1,out,1)
cc
ccc         Kill solver
cc
cc          call solverKill
cc
ccc         Scatter solution
cc
cc          do i = 1,ntotp
cc            do ieq=1,neq
cc              x(i,ieq) = xi(neq*(i-1)+ieq)
cc            enddo
cc          enddo
cc
cc          return
        endif

c     Solve coupled MG

c       Initialize solver

        call solverInit

c       Check whether GMRES is driving solve

        if (gm_drvr .or. cg_drvr) then

          call solverOptionsInit

          solverOptions%iter    = gmit
          solverOptions%tol     = rtol

          solverOptions%stp_test = 1 
          solverOptions%krylov_subspace = solverOptions%iter 

          if (gm_drvr) then
            solverOptions%singular_matrix = smat
            call assembleSolverHierarchy('gm')
          else
            call assembleSolverHierarchy('cg')
          endif

          if (vcyc == 0) then !Identity preconditioner
            call assembleSolverHierarchy('id')
          endif
        endif

c       MG solver options

        call solverOptionsInit

        solverOptions%tol      = rtol

        if (gm_drvr.or.cg_drvr) then
          solverOptions%vcyc   = 1
        else
          solverOptions%vcyc   = vcyc
        endif
        solverOptions%mg_gmin  = mgmin
        solverOptions%orderres = iorder_res
        solverOptions%orderprol= iorder_prol
        solverOptions%mg_mu    = mu_cyc
        solverOptions%mg_galerkin = glkin
cc        solverOptions%vol_res  = vol_wgt

        if (PRESENT(mg_debug)) solverOptions%mg_debug = mg_debug

        if (PRESENT(dg)) solverOptions%diag => ddg

        if (PRESENT(mg_grid_def)) solverOptions%mg_grid_def=>mg_grid_def

        solverOptions%ncolors  = nclrs

        !Coarse solve
        if (gm_crse_solve) then
        !With GMRES as coarse solver (defined below)
          if (         (gm_drvr  .and. smth/='gm')
     .       .or.((.not.gm_drvr) .and. smth=='gm') ) then
            solverOptions%mg_coarse_solver_depth = 4
          elseif (gm_drvr .and. smth=='gm') then
            solverOptions%mg_coarse_solver_depth = 5
          else
            solverOptions%mg_coarse_solver_depth = 3
          endif
        else
          solverOptions%mg_coarse_solver_depth = 0 !Use smoother as coarse solver
cc          write (*,*) 'DIAG -- cSolver'
cc          solverOptions%mg_coarse_solver_depth = 5
        endif

        !Vertex relaxation
cc        solverOptions%vertex_based_relax = .true.

        !Plane/line relaxation
        if (lrelax) then
cc          solverOptions%mg_gmin         = mgmin
cc          solverOptions%orderres        = ores
cc          solverOptions%orderprol       = oprol
cc          solverOptions%mg_mu           = mu_cyc
          solverOptions%mg_line_relax   = .true.
          solverOptions%mg_line_nsweep  = 2
          solverOptions%mg_line_vcyc    = 100
          solverOptions%mg_line_tol     = 1d-1
          solverOptions%mg_line_omega   = 1d0
          solverOptions%mg_line_x       = .true.
          solverOptions%mg_line_y       = .true.
          solverOptions%mg_line_z       = .false.
          solverOptions%mg_line_solve   = "jb"
          solverOptions%mg_line_coarse_solver_depth = 0
        endif

        call assembleSolverHierarchy('mg')

c       MG smoother

        call solverOptionsInit

        solverOptions%iter            = smit
        solverOptions%tol             = rtol
        solverOptions%singular_matrix = smat

        if (lrelax) then
          solverOptions%iter  = 2
        endif

        select case (smth)
        case('gs')
          solverOptions%omega   = omg
          solverOptions%ncolors = nclrs
          solverOptions%mg_zebra_relax = zrelax
          solverOptions%mg_zebra_prefd = zdir
cc          if(bcSP()) solverOptions%mg_zebra_prefd = 1
          solverOptions%mg_zebra_it    = zit
          solverOptions%mg_zebra_omega = zomg

cc          solverOptions%iter=ceiling(1.*smit/solverOptions%mg_zebra_it)
          solverOptions%iter=smit

          call assembleSolverHierarchy('gs')
        case('gm')
          solverOptions%stp_test        = 1 
          solverOptions%krylov_subspace = solverOptions%iter 

          call assembleSolverHierarchy('gm')

          call assembleSolverHierarchy('id')
        case default !JB
          solverOptions%omega   = omg
          call assembleSolverHierarchy('jb')
        end select

c       MG coarsest grid solve

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = min(100,2*neq*ntotp)
        solverOptions%iter            = 100
        solverOptions%stp_test        = 1 
        solverOptions%omega           = 1d0
        solverOptions%singular_matrix = smat

        call assembleSolverHierarchy('gm')
        call assembleSolverHierarchy('id') !GMRES preconditioner

c       Coarsest grid solve for inner line/plane MG

        call solverOptionsInit

cc        if (PRESENT(dg)) solverOptions%diag => ddg

        solverOptions%omega   = 1d0
        solverOptions%iter    = 100
        solverOptions%tol     = 1d-4
        solverOptions%ncolors = nclrs

cc        call assembleSolverHierarchy('jb')
        call assembleSolverHierarchy('gs')

c       Invoke solver

        ntot=neq*ntotp
        call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,guess,out,1)

c       Get output data

        call getSolverOptions(1)

        if(PRESENT(iters)) iters = solverOptions%iter_out

c       Kill solver

        call solverKill

c       Scatter solution for output

        do i = 1,ntotp
          do ieq=1,neq
            x(i,ieq) = xi(neq*(i-1)+ieq)
          enddo
        enddo

c     End program

      end subroutine cSolver

c     MGmatvec
c     ###############################################################
      subroutine MGmatvec(mgctx,gposc,x,y,igc,bcnd,matvec)
c     ---------------------------------------------------------------
c     This subroutine calculates the matrix vector product defined
c     symbolically in matvec, in two ways:
c       * Rediscretizing (i.e., using matvec straight)
c       * Galerkin (preprocessing w/ prolongation, and postprocessing
c                   w/ restriction).
c
c     In call:
c      * mgctx: MG context
c      * gposc: vector index of position on the numerical grid:
c            + If gposc = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gposc = 0, all the grid is considered.
c            + If gposc < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gposc) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntotc: total number of unknowns: neq*nx*ny*nz
c      * x(ntotc): input vector
c      * y(ntotc): output vector
c      * igc: grid level
c      * bcnd: boundary conditions on x vector.
c      * matvec: external with symbolic matrix-vector code.
c      * volw: whether residuals are volume weighed
c     ---------------------------------------------------------------

      implicit none

c Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: igc,gposc,bcnd(:,:)
      real(8) :: x(:),y(:)

      external   matvec

c Local variables

      integer :: neq,nxf,nyf,nzf,ntotf,igf,gposf,if,jf,kf
     .          ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf
      integer :: ntotc,nxc,nyc,nzc,isigc,ic,jc,kc
     .          ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc
      integer :: mgrx,mgry,mgrz,ig
      real(8),allocatable,dimension(:) :: MGx,MGy

      logical :: fpointers

c Begin program

      neq = mgctx%neq

      nxc = mgctx%g_def%nxv(igc)
      nyc = mgctx%g_def%nyv(igc)
      nzc = mgctx%g_def%nzv(igc)

      ntotc = size(x)

      call allocPointers(neq,fpointers,mgctx=mgctx)

      if (igc > 1 .and. galerkin) then

        isigc  = mgctx%istart(igc)

c     Allocate MG vector

        igf = 1

        nxf = mgctx%g_def%nxv(igf)
        nyf = mgctx%g_def%nyv(igf)
        nzf = mgctx%g_def%nzv(igf)

        ntotf = neq*nxf*nyf*nzf

        allocate(MGx(2*ntotf),MGy(2*ntotf))

c     Prolong vector to finest grid

        MGx(isigc:isigc+ntotc-1) = x

        call prolongMGVector(mgctx,MGx,igc,igf,0,bcnd)

c     Perform matvec in finest grid

        if (gposc == 0) then

          call matvec(0,neq,ntotf,MGx(1:ntotf),MGy(1:ntotf),igf,bcnd)

        else

          !Find coarse cell position (ic,jc,jc)
          call limits(gposc,nxc,nyc,nzc,igc
     .               ,ic,imaxc,jc,jmaxc,kc,kmaxc)

          !Define finest grid stencil
          mgrx = product(mgctx%mg_ratio_x(igf:igc-1))
          mgry = product(mgctx%mg_ratio_y(igf:igc-1))
          mgrz = product(mgctx%mg_ratio_z(igf:igc-1))

          iminf = mgrx*(ic-1)+1
          imaxf = mgrx* ic

          jminf = mgry*(jc-1)+1
          jmaxf = mgry* jc

          kminf = mgrz*(kc-1)+1
          kmaxf = mgrz* kc

          !Perform matvec on fine grid stencil
          do kf = kminf,kmaxf
            do jf = jminf,jmaxf
              do if = iminf,imaxf
                gposf = if + nxf*(jf-1) + nxf*nyf*(kf-1)
                call matvec(gposf,neq,ntotf,MGx(1:ntotf),MGy(1:ntotf)
     .                     ,igf,bcnd)
              enddo
            enddo
          enddo
          
        endif

c     Restrict residual

        call restrictMGVector(mgctx,MGy,igf,igc,0,res_vol_wgt)

c     Map coarse grid residual

        y = MGy(isigc:isigc+ntotc-1)

c     Deallocate MG vectors

        deallocate(MGx,MGy)

      else

        call matvec(gposc,neq,ntotc,x,y,igc,bcnd)

      endif

      call deallocPointers(fpointers,mgctx=mgctx)

c     End program

      end subroutine MGmatvec

c     find_mf_diag_std
c     ###############################################################
      subroutine find_mf_diag_std(mgctx,matvec,igrid,bbcnd,diag1)
c     ---------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec
c     for all grids, commencing with grid "igrid".
c     ---------------------------------------------------------------

      implicit none      !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer    :: bbcnd(:,:),igrid

      real(8)    :: diag1(:,:)

      external      matvec

c     Local variables

      real(8)    :: x1( size(diag1,2)/2),dummy(size(diag1,2)/2)
     .             ,mat(size(diag1,1),size(diag1,1))
     .             ,mat2(size(diag1,1),size(diag1,1))
     $             ,idet
      integer    :: neq,ntot,ii,jj,nn,ig,iig,isig,igmax
      integer    :: igr,ieq,alloc_stat
      logical    :: fpointers

c     Begin program

      neq  = size(diag1,1)
      ntot = size(diag1,2)/2

c     Allocate MG pointers

      call allocPointers(neq,fpointers,mgctx=mgctx)

cc      if (fpointers) igmax=MGgrid%ngrid  !Routine not called from JB,GS,MG
      igmax = mgctx%igmax

c     Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        call pstop('find_mf_diag_std'
     .            ,'Error in input of find_mf_diag_stc')
      endif

c     Form diagonal

      do igr = igrid,igmax

c     Form diagonal terms for smoother

        nn  = ntotv(igr)

        isig = istart(igr)

        x1(1:nn) = 0d0
        dummy(1:nn) = 0d0

c     Finds block diagonals for neq equations.

        do ii = 1,ntotvp(igr)

          jj  = (ii-1)*neq + isig - 1

          do ieq = 1,neq

c         Find column vector corresponding to grid node ii and equation ieq

            x1(neq*(ii-1) + ieq) = 1d0

            call MGmatvec(mgctx,ii,x1,dummy,igr,bbcnd,matvec)

            x1(neq*(ii-1) + ieq) = 0d0

c         Fill diagonal

            diag1(ieq,jj+1:jj+neq) = dummy((ii-1)*neq+1:ii*neq)

          enddo

        enddo

c     Invert diagonal and store in diag1

        select case (neq)
        case (1)

          do ii = 1,ntotvp(igr)
            ig = ii + isig - 1
            diag1(neq,ig) = 1d0/diag1(neq,ig)
          enddo

cc        case (2)
cc
cc          do ii = 1,ntotvp(igr)
cc            iig = neq*(ii - 1) + isig - 1
cc            idet = 1d0/(diag1(1,iig+1)*diag1(2,iig+2)
cc     .                 -diag1(2,iig+1)*diag1(1,iig+2))
cc            mat(1,1) = diag1(2,iig+2)*idet
cc            mat(1,2) =-diag1(2,iig+1)*idet
cc            mat(2,1) =-diag1(1,iig+2)*idet
cc            mat(2,2) = diag1(1,iig+1)*idet
cc
cc            diag1(:,iig+1:iig+neq) = mat
cc          enddo

        case default

          !Invert
          do ii = 1,ntotvp(igr)
            iig = neq*(ii - 1) + isig - 1
            call blockInv(neq,diag1(:,iig+1:iig+neq))
          enddo

        end select

      enddo

c     Deallocate pointers

      call deallocPointers(fpointers,mgctx=mgctx)

c End program

      end subroutine find_mf_diag_std

c     find_mf_diag_colored
c     ###############################################################
      subroutine find_mf_diag_colored(mgctx,matvec,igrid,bbcnd,dg,nclrs)
c     ---------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec
c     for all grids, commencing with grid "igrid". This routine
c     employs a colored arrangement of the grid that allows one
c     to find the diagonal using full-vector matvecs, instead of
c     point-wise matvecs.
c     ---------------------------------------------------------------

      implicit none      !For safe fortran

c     Call variables

      type(mg_ctx),pointer :: mgctx

      integer :: bbcnd(:,:),igrid,nclrs

      real(8) :: dg(:,:)

      external   matvec

c     Local variables

      real(8)    :: x1( size(dg,2)/2),dummy(size(dg,2)/2)
     .             ,mat(size(dg,1),size(dg,1))
     .             ,mat2(size(dg,1),size(dg,1))
     $             ,idet,coeff
      integer :: neq,ntot,ii,jj,nn,ig,iig,isig,bbcnd2(6,size(dg,1))
      integer :: igr,ieq,alloc_stat
      logical :: fpointers

      integer :: i,j,k,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3,igmax
      integer :: i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max

c     Begin program

      neq  = size(dg,1)
      ntot = size(dg,2)/2

c     Allocate MG pointers

      call allocPointers(neq,fpointers,mgctx=mgctx)

cc      if (fpointers) igmax=MGgrid%ngrid  !Routine not called from JB,GS,MG
      igmax = mgctx%igmax

c     Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        call pstop('find_mf_diag_colored'
     .            ,'Error in input of find_mf_diag_colored')
      endif

c     Form diagonal

      dg = 0d0

      if (nclrs == 1) then

        call find_mf_diag_std(mgctx,matvec,igrid,bbcnd,dg)

      else

        do igr = igrid,igmax

c       Form diagonal terms for smoother

          nn  = ntotv(igr)

          isig = istart(igr)

          x1   (1:nn) = 0d0
          dummy(1:nn) = 0d0

          bbcnd2 = bbcnd
          where (bbcnd == SP) bbcnd2 = DEF !Avoid SP bc's for forming diagonal (they're nonlocal)

c       Finds block diagonals for neq equations.

          call scheduleColors(mgctx,igr,nclrs,nrbg1,nrbg2,nrbg3
     .                       ,i1min,i1max,j1min,j1max,k1min,k1max)

          do ieq = 1,neq

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                  do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                    do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                        call scheduleIndices(mgctx,igr,i1,j1,k1,i,j,k)

                        ii = i +nxv(igr)*(j-1) +nxv(igr)*nyv(igr)*(k-1)

                        x1((ii-1)*neq+ieq) = 1d0
                      enddo
                    enddo
                  enddo

                  call MGmatvec(mgctx,0,x1,dummy,igr,bbcnd2,matvec)

                  do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                    do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                        call scheduleIndices(mgctx,igr,i1,j1,k1,i,j,k)

                        ii = i +nxv(igr)*(j-1) +nxv(igr)*nyv(igr)*(k-1)

                        jj = (ii-1)*neq + isig - 1

                        x1((ii-1)*neq+ieq) = 0d0

                        dg(ieq,jj+1:jj+neq) = dummy((ii-1)*neq+1:ii*neq)
                      enddo
                    enddo
                  enddo

                enddo
              enddo
            enddo

          enddo

c       Take care of points around singular point

cccc          if (isSP(1,1,1,igr,igr,igr)) then
cc          if (bcSP()) then
cc
cc#if defined(petsc)
cccc            call pstop('find_mf_diag_colored'
cccc     .                ,'Does not work in parallel')
cc#endif
cc
cc            !Concentrate around singular point
cc            i = 1
cc            coeff = 0d0
cc            if (isSP2(i,igr)) coeff = 1d0
cc
cccc            MGgrid%iline(igr) = i !Set line smoothing around SP
cccc
cccc            call scheduleColors(igr,nclrs,nrbg1,nrbg2,nrbg3
cccc     .                         ,i1min,i1max,j1min,j1max,k1min,k1max)
cc
cccc            call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
cccc     .                 ,i1min,i1max,j1min,j1max,k1min,k1max)
cc
cc            do ieq = 1,neq
cc
cc              do k=k1min,k1max  !Need to do 2 colors in k; not yet implemented
cc                do j=j1min,j1max
cc
cccc              do irbg3 = 1,nrbg3
cccc                do irbg2 = 1,nrbg2
cccc                  do irbg1 = 1,nrbg1
cccc
cccc                  do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
cccc                    do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
cccc                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1
cc
cc                  ii = i +nxv(igr)*(j-1) +nxv(igr)*nyv(igr)*(k-1)
cc
cc                  !Perform matvec
cc
cc                  x1(neq*(ii-1) + ieq) = coeff
cc
cc                  call MGmatvec(ii,neq,nn,x1,dummy,igr,bbcnd,matvec)
cc
cc                  x1(neq*(ii-1) + ieq) = 0d0
cc
cc                  !Fill diagonal
cc
cc                  jj  = (ii-1)*neq + isig - 1
cc
cc                  if (coeff == 1d0)
cc     .               dg(ieq,jj+1:jj+neq) = dummy((ii-1)*neq+1:ii*neq)
cc
cc                enddo
cc              enddo
cc
cc            enddo
cc
cc          endif

c       Invert diagonal and store in dg

          select case (neq)
          case (1)

            do ii = 1,ntotvp(igr)
              ig = ii + isig - 1
              dg(neq,ig) = 1d0/dg(neq,ig)
            enddo

cc        case (2)
cc
cc          do ii = 1,ntotvp(igr)
cc            iig = neq*(ii - 1) + isig - 1
cc            idet = 1d0/(dg(1,iig+1)*dg(2,iig+2)
cc     .                 -dg(2,iig+1)*dg(1,iig+2))
cc            mat(1,1) = dg(2,iig+2)*idet
cc            mat(1,2) =-dg(2,iig+1)*idet
cc            mat(2,1) =-dg(1,iig+2)*idet
cc            mat(2,2) = dg(1,iig+1)*idet
cc
cc            dg(:,iig+1:iig+neq) = mat
cc          enddo

          case default

            !Invert
            do ii = 1,ntotvp(igr)
              iig = neq*(ii - 1) + isig - 1
              call blockInv(neq,dg(:,iig+1:iig+neq))
            enddo

          end select

        enddo

      endif

c     Deallocate pointers

      call deallocPointers(fpointers,mgctx=mgctx)

c     End program

      end subroutine find_mf_diag_colored

c     find_mf_diag
c     ###############################################################
      subroutine find_mf_diag(neq,nblk,ntot,matvec,igrid,bbcnd
     .                       ,dg,nclrs,volw,glrkin,mgctx,g_def)
c     ---------------------------------------------------------------
c     Finds diagonal elements for relaxation. This is done matrix-free
c     using subroutine matvec for all grids, commencing with grid "igrid".
c     
c     Careful: for vertex-based relax., diagonal not formed for coarsest
c     grid since no additional grid levels are available.
c     MG needs additional coarsest level solver (such as GMRES).
c
c     In call sequence:
c       * neq: number of coupled equations
c       * nblk: number of fine grid cells per coarse cell.
c           If nblk=1, we use standard single-node relaxation.
c           Otherwise, we use vertex-based relaxation.
c       * ntot: dimension of vectors (neq*nx*ny*nz)
c       * matvec: external indicating matrix-vector routine.
c       * igrid: current grid level
c       * bbcnd: boundary condition information (passed on to matvec)
c       * dg: returns inverse of diagonal.
c       * nclrs: number of colors for colored algorithm
c       * volw: whether residuals are volume-weighed
c       * glrkin: whether we use galerkin coarsening
c       * mgctx (optional): MG context
c       * g_def (optional): grid definition context
c     ---------------------------------------------------------------

      implicit none      !For safe fortran

c     Call variables

      integer  :: neq,nblk,ntot,bbcnd(6,neq),igrid,nclrs

      real(8)  :: dg(neq*nblk,2*ntot*nblk)

      logical  :: volw,glrkin

      external    matvec

      type(mg_ctx),pointer,optional :: mgctx

      type(grid_mg_def),pointer,optional :: g_def

c     Local variables

      real(8) :: x1(ntot),dummy(ntot),mag(3)
      integer :: ii,jj,nn,ig,isig,size,igmax
      integer :: iv,jv,kv,if,jf,kf,nxf,nyf,nzf
      integer :: iii,icolb,irowb,icol,iiv,iig,ivg,irowmin,irowmax
     .          ,iblock ,jblock ,kblock
     .          ,iblock2,jblock2,kblock2,nblkg

      integer :: igr,igc,ieq,alloc_stat
      logical :: fpointers

c     Begin program

c     Set internal flags

      res_vol_wgt = volw
      galerkin = glrkin

c     Allocate MG pointers

      call allocPointers(neq,fpointers,mgctx=mgctx,g_def=g_def)

cc      if (fpointers) igmax=MGgrid%ngrid  !Routine not called from JB,GS,MG
      igmax = mg_ctx_temp%igmax

c     Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        call pstop('find_mf_diag'
     .            ,'Error in input of find_mf_diag')
      endif

c     Form diagonal

c     STANDARD RELAXATION
      if (nblk == 1) then

        call find_mf_diag_colored(mg_ctx_temp,matvec,igrid,bbcnd,dg
     .                           ,nclrs)

C     VERTEX-BASED RELAXATION
      else

        dg = 0d0

        do igr = igrid,igmax

          nxf = mgctx%g_def%nxv(igr)
          nyf = mgctx%g_def%nyv(igr)
          nzf = mgctx%g_def%nzv(igr)

          nblkg = nblock(igr)

          nn    = ntotv(igr)

          isig  = istartb(igr)

          size  = neq*nblkg

          x1(1:nn) = 0d0

c       Sample vertices

          do kv=1,max(nzf-mg_ratio_z(igr)+1,1)
            do jv=1,max(nyf-mg_ratio_y(igr)+1,1)
              do iv=1,max(nxf-mg_ratio_x(igr)+1,1)

              !This defines vertex numbers (lexicographic)
              iiv = iv + nxf*(jv-1) + nxf*nyf*(kv-1)
              ivg  = (iiv-1)*neq*nblkg + isig -1

              !Sample columns in diagonal block
              do kblock=1,mg_ratio_z(igr)
                do jblock=1,mg_ratio_y(igr)
                  do iblock=1,mg_ratio_x(igr)

                    !This defines block column index (lexicographic)
                    icolb = iblock + mg_ratio_x(igr)*(jblock-1)
     .                             + mg_ratio_x(igr)
     .                              *mg_ratio_y(igr)*(kblock-1)

                    !This defines cell positions around vertex (iv,jv,kv) (lexicographic)
                    if = iv-1 + iblock
                    jf = jv-1 + jblock
                    kf = kv-1 + kblock

                    ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                    !Sample rows in diagonal block
                    do kblock2=1,mg_ratio_z(igr)
                      do jblock2=1,mg_ratio_y(igr)
                        do iblock2=1,mg_ratio_x(igr)
 
                          !This defines block row index
                          irowb =iblock2 + mg_ratio_x(igr)*(jblock2-1)
     .                                   + mg_ratio_x(igr)
     .                                    *mg_ratio_y(igr)*(kblock2-1)

                          !This defines cell position around vertex (iv,jv,kv) (lexicographic)
                          if = iv-1 + iblock2
                          jf = jv-1 + jblock2
                          kf = kv-1 + kblock2

                          jj  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                          !Sample equations
                          do ieq = 1,neq

                            !This defines column index of system matrix
                            iii  = neq*(ii   -1) + ieq
                            !This defines column index in diagonal block
                            icol = neq*(icolb-1) + ieq
                            !This defines min and max row indices in diagonal block
                            irowmin = ivg + neq*(irowb-1)+1
                            irowmax = ivg + neq* irowb

                            !Set column vector corresponding to grid node ii and equation ieq
                            x1(iii) = 1d0

                            !Find matrix components corresponding to column iii and rows
                            !((jj-1)*neq+1) to (jj*neq)
                            call MGmatvec(mg_ctx_temp,jj,x1,dummy,igr
     .                                   ,bbcnd,matvec)

                            !Reset vector x1 to zero
                            x1(iii) = 0d0

                            !Store results in diagonal
                            dg(icol,irowmin:irowmax)
     .                                = dummy((jj-1)*neq+1:jj*neq)
                          enddo

                        enddo
                      enddo
                    enddo

                  enddo
                enddo
              enddo

c             Invert diagonal and store in dg

              call blockInv(size,dg(1:size,ivg+1:ivg+size))

              enddo
            enddo
          enddo

        enddo

      endif

c     Deallocate pointers

      call deallocPointers(fpointers,mgctx=mgctx)

c     End program

      end subroutine find_mf_diag

c     find_mf_mat
c     ###############################################################
      subroutine find_mf_mat(neq,ntot,matvec,igrid,bbcnd,mat,mgctx)
c     ---------------------------------------------------------------
c     Finds all matrix elements matrix-free using subroutine matvec
c     for all grids, commencing with grid "igrid".
c     ---------------------------------------------------------------

      implicit none      !For safe fortran

c     Call variables

      integer    :: neq,ntot,bbcnd(6,neq),igrid

      real(8)    :: mat(ntot,ntot)

      external      matvec

      type(mg_ctx),pointer,optional :: mgctx

c     Local variables

      real(8)    :: x1(ntot),dummy(ntot)
      integer    :: ii,jj,nn
      integer    :: igr,ieq,alloc_stat
      logical    :: fpointers

c     Begin program

      igr = igrid

c     Allocate MG pointers

      call allocPointers(neq,fpointers,mgctx=mgctx)

c     Consistency check

      nn  = ntotv(igr)

      if (nn /= ntot) then
        call pstop('find_mf_mat'
     .            ,'Error in input of find_mf_mat')
      endif

c     Form diagonal 

      x1 = 0d0

      do ii = 1,nn

c       Find column vector corresponding to grid node ii and equation ieq

          x1(ii) = 1d0

          call MGmatvec(mg_ctx_temp,0,x1,dummy,igr,bbcnd,matvec)

          x1(ii) = 0d0

c       Fill matrix

          mat(:,ii) = dummy

      enddo

c     Deallocate pointers

      call deallocPointers(fpointers,mgctx=mgctx)

c     End program

      end subroutine find_mf_mat

c     findBaseVector
c     ###############################################################
      subroutine findBaseVector(ii,ieq,neq,ntot,x1,coef)
c     ---------------------------------------------------------------
c     Finds base vector corresponding to grid node ii and equation ieq.
c     Assumes that x1 has been initialized to zero elsewhere.
c     ---------------------------------------------------------------

      implicit none      !For safe fortran

c     Call variables

      integer    :: neq,ii,ieq,ntot

      real(8)    :: x1(ntot),coef

c     Local variables

c     Begin program

      x1(neq*(ii-1) + ieq) = coef

c     End program

      end subroutine findBaseVector

      end module mg_solver

c mg
c#######################################################################
      recursive subroutine mg(neq,ntot,y,x,matvec,options,igrid,bcnd
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free coupled MG routine to solve
c     Ax = y. Call variables:
c       * neq: number of equations
c       * ntot: total vector size (grid dimension*neq)
c       * y,x: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level to start from in MG applications
c       * bcnd: integer array containing BC information for unknown.
c       * guess: 0->no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c
c    Grid convention: finest grid has igrid=1, coarsest has igrid=ngrid
c
c--------------------------------------------------------------------

      use mg_solver

      implicit none       !For safe fortran

c Call variables

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer :: iter,igridmin,vcyc,crsedpth,isig,ncolors
      integer :: orderres,orderprol,alloc_stat
     .          ,line_nsweep,line_crse_depth,line_vcyc

      integer :: guess2,outc,mu
      integer :: izero,i,j,k,ii,ivcyc,nblk,igmax

      real(8) :: xx(2*ntot),yy(2*ntot),wrk(2*ntot),rr(ntot),x0(ntot)
      real(8) :: rr0,rr1,mag,mag1,mag2,mgtol,line_tol,line_omega

      real(8) :: epsmac=0d0,abstol

      real(8), pointer, dimension(:,:) :: diag

      logical :: fdiag,line_relax,fpointers,line_x,line_y,line_z
     .          ,mg_debug,prnt

      character(2) :: line_solve

      type(mg_ctx),pointer :: mgctx

c Begin program

      if (epsmac == 0d0) epsmac = find_mach_eps()

      abstol = epsmac*sqrt(1.*ntot)

      prnt = (my_rank == 0) .or. asm

c Consistency check

      call optionsConsistencyCheck

c Assign options

      igridmin        = options%mg_gmin
      vcyc            = options%vcyc
      mgtol           = options%tol
      orderres        = options%orderres
      orderprol       = options%orderprol
      fdiag           = options%fdiag
cc      volf            = options%vol_res
      crsedpth        = options%mg_coarse_solver_depth
      mu              = options%mg_mu
      vbr_mg          = options%mg_vertex_relax
      galerkin        = options%mg_galerkin
      diag            => options%diag
      mgctx           => options%mgctx
      MGgrid          => options%mg_grid_def

      line_relax      = options%mg_line_relax
      line_vcyc       = options%mg_line_vcyc
      line_nsweep     = options%mg_line_nsweep
      line_tol        = options%mg_line_tol
      line_omega      = options%mg_line_omega
      line_crse_depth = options%mg_line_coarse_solver_depth 
      line_solve      = options%mg_line_solve
      line_x          = options%mg_line_x
      line_y          = options%mg_line_y
      line_z          = options%mg_line_z

      ncolors         = options%ncolors

      mg_debug        = options%mg_debug

c Set pointers and find ngrid

      call allocPointers(neq,fpointers,mgctx=mgctx)

      if (fpointers.and.out.ge.2.and.my_rank==0)
     .     write (*,*) 'Allocating pointers...'

c Check limits

      mgctx%igmax = ngrid - igridmin + 2  !Determines minimum resolution

      igmax = mgctx%igmax

      if (igmax.lt.igrid) then
        call pstop('mg',' Grid chosen is below minimum grid resolution')
      endif

      if (out.ge.2.and.igmax.gt.igrid.and.prnt) write (*,5)

      outc = out
      if (igmax.gt.igrid) outc = out - 2

c Find diagonal for smoothers

      if (vbr_mg) then
        nblk = nblock(igrid)
      else
        nblk = 1
      endif

cc      if (fdiag) then
cc
cc        if (allocated(diag)) then !Diagonal is known
cc          if(out.ge.2.and.prnt) write(*,*)'Diagonal already allocated'
cc          fdiag = .false.
cc        elseif (associated(options%diag)) then !Diagonal provided externally
cc          if(out.ge.2.and.prnt) write(*,*)'Diagonal externally provided'
cc          allocate(diag(neq*nblk,2*ntot*nblk))
cc          diag = options%diag
cc        else                      !Form diagonal
cc          if (out.ge.2.and.prnt) write (*,*) 'Forming diagonal...'
cc          allocate(diag(neq*nblk,2*ntot*nblk))
cc          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
cc     .                     ,ncolors,volf,galerkin)
cc          if (out.ge.2.and.prnt) write (*,*) 'Finished!'
cc        endif
cc
cc      endif

      if (associated(diag)) then !Diagonal provided externally
        if(out.ge.2.and.prnt) write(*,*)'Diagonal externally provided'
        fdiag = .false.
      else                      !Form diagonal
        if (out.ge.2.and.prnt) write (*,*) 'Forming diagonal...'
        allocate(diag(neq*nblk,2*ntot*nblk))
        call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
     .                   ,ncolors,res_vol_wgt,galerkin,mgctx=mgctx)
        if (out.ge.2.and.prnt) write (*,*) 'Finished!'
        fdiag = .true.
      endif

c Initialize local solution vector (xx) and local r.h.s. (yy)

      xx  = 0d0
      yy  = 0d0
      wrk = 0d0

      if (guess.eq.0) then
        x = 0d0
      else
        xx(istart(igrid):istart(igrid+1)-1) = x(:)
      endif

      yy(istart(igrid):istart(igrid+1)-1) = y(:)

c Compute initial residual and check convergence

      if (guess.eq.0) then
        rr0 = sqrt(MGdot(mgctx,igrid,y,y))
      else
        call MGmatvec(mgctx,0,x,rr,igrid,bcnd,matvec)
        call MGvecadd(mgctx,igrid,-1d0,rr,1d0,y)
        rr0 = sqrt(MGdot(mgctx,igrid,rr,rr))
      endif

      if (rr0.lt.abstol) then
        if (out.ge.1.and.prnt) then
          write (*,*) 'Initial solution seems exact in MG'
cc          write (*,'(a,1pe10.2,a,e10.2)') '    Residual=',rr0
cc     .          ,' < limit =',1d-16*ntot
        endif
        call killmg
        return
      endif

      rr1 = rr0

c Start mu-cycle

      guess2 = 1

      do ivcyc = 1,vcyc

        if (outc.ge.1.and.igmax.gt.igrid.and.prnt) then
          if (mu == 1) write (*,7) ivcyc
          if (mu == 2) write (*,8) ivcyc
        endif

        x0 = xx(istart(igrid):istart(igrid)+ntot-1)  !Save initial solution

c     Perform mu-cycle recursively

        if (igrid.eq.igmax) then
          call smooth (igrid)
          exit
        else
          call mcycle(igrid)
        endif

c     Check MG convergence

        call MGmatvec(mgctx,0,xx(istart(igrid):istart(igrid)+ntot-1)
     .               ,rr,igrid,bcnd,matvec)

        call MGvecadd(mgctx,igrid,-1d0,rr,1d0,y)

        mag = sqrt(MGdot(mgctx,igrid,rr,rr))

        mag1 = mag/rr1

        x0 = xx(istart(igrid):istart(igrid)+ntot-1)-x0
        mag2 = sqrt(MGdot(mgctx,igrid,x0,x0)/ntot)

        if (out.ge.3.and.prnt) then
          write (*,10) mag,mag1
        elseif (out.ge.2.and.prnt) then
          if (mu == 1) write (*,20) mag,mag2,mag1,ivcyc
          if (mu == 2) write (*,21) mag,mag2,mag1,ivcyc
        endif

        rr1 = mag

        if (mag/rr0 < mgtol .or. mag < 1d-20*ntot) exit

      enddo

c Map solution from local vector xx to external vector x

      x(:) = xx(istart(igrid):istart(igrid+1)-1)

c MG convergence info

      mag1 = mag/rr0

      if (igmax.gt.igrid) then
        if (out.eq.1.and.prnt) then
          if (mu == 1) write (*,20) mag,mag2,mag1,min(ivcyc,vcyc)
          if (mu == 2) write (*,21) mag,mag2,mag1,min(ivcyc,vcyc)
        elseif (out.ge.2.and.vcyc.gt.1.and.prnt) then
          write (*,*) 
          write (*,*) 'Final MG convergence info:'
          if (mu == 1) write (*,20) mag,mag2,mag1,min(ivcyc,vcyc)
          if (mu == 2) write (*,21) mag,mag2,mag1,min(ivcyc,vcyc)
          write (*,*)
        endif
      endif

c End program

      options%tol_out = mag1

      call killmg

      return

 5    format (/,' MG method output:')
 7    format (/,' MG V-cycle #:',i3)
 8    format (/,' MG W-cycle #:',i3)
 10   format (  ' MG residual:',1p,1e12.4,'; Ratio:',1e12.4)
 20   format (  ' MG residual:',1p,1e12.4,'; MG update:',1e12.4,
     .          '; MG ratio:',1e12.4,'; V-cycle #:',i3)
 21   format (  ' MG residual:',1p,1e12.4,'; MG update:',1e12.4,
     .          '; MG ratio:',1e12.4,'; W-cycle #:',i3)

      contains

c     optionsConsistencyCheck
c     ###################################################################
      subroutine optionsConsistencyCheck

c     -------------------------------------------------------------------
c       Checks consistency of MG input options.
c     -------------------------------------------------------------------

        if (     options%mg_vertex_relax
     .      .and.options%mg_line_relax) then
          if (prnt)
     .       write (*,*) 'Invalid setting for MG relaxation: cannot do'
          if (prnt)
     .       write (*,*) 'vertex-based and line-based relaxation',
     .                ' simultaneously!'
          call pstop('optionsConsistencyCheck in mg','Aborting...')
        endif

        if (     options%mg_coarse_solver_depth == 0
     .      .and.options%mg_vertex_relax) then
          if (prnt)
     .         write (*,*) 'Coarsest level solver not defined'
          if (prnt)
     .         write (*,*) 'Cannot do vertex-based relaxation in MG'
          call pstop('optionsConsistencyCheck','Aborting...')
        endif
          
      end subroutine optionsConsistencyCheck

c     mcycle
c     ###################################################################
      recursive subroutine mcycle(igr)

c     -------------------------------------------------------------------
c       Performs MG M-cycle (V-cycle, m=1, W-cycle, m=2) recursively
c       starting at grid igr.
c     -------------------------------------------------------------------

        use math

        use grid_debug, ONLY:createDrawInCfile

        implicit none

c     Call variables

        integer    :: igr

c     Local variables

        integer    :: igc,isigc,isig,nn,nnc,imu,iplot
        character(50) :: dbg_desc(4*neq)

cc        external     matvec

c     Begin program

        igc   = igr+1

        nn    = ntotv(igr)
        nnc   = ntotv(igc)

        isig  = istart(igr)
        isigc = istart(igc)

c     Setup debugging plots

        if (mg_debug) then

          do iplot=1,neq
            dbg_desc(iplot)       = 'F-sol_'//trim(int2char(iplot))
     .                              //' after pre-smoothing'
            dbg_desc(iplot+neq)   = 'F-sol_'//trim(int2char(iplot))
     .                              //' correction'
            dbg_desc(iplot+2*neq) = 'F-sol_'//trim(int2char(iplot))
     .                              //' after V-cycle'
          enddo

          call createDrawInCfile(3*neq,'finesol.bin','Fine grid Sol'
     .          ,'t','x','y',dbg_desc,'-c -X0 -L57'
     .          ,'drawfsol.in')

          do iplot=1,neq
            dbg_desc(iplot)       = 'F-res_'//trim(int2char(iplot))
     .                              //' before pre-smoothing'
            dbg_desc(iplot+neq)   = 'F-res_'//trim(int2char(iplot))
     .                              //' after pre-smoothing'
            dbg_desc(iplot+2*neq) = 'F-res_'//trim(int2char(iplot))
     .                              //' after coarse updt'
            dbg_desc(iplot+3*neq) = 'F-res_'//trim(int2char(iplot))
     .                              //' after post-smoothing'
          enddo

          call createDrawInCfile(4*neq,'fineres.bin','Fine grid Res'
     .          ,'t','x','y',dbg_desc,'-c -X0 -L57'
     .          ,'drawfres.in')

          do iplot=1,neq
            dbg_desc(iplot)       = 'C-res_'//trim(int2char(iplot))
     .                              //' after restriction'
            dbg_desc(iplot+neq)   = 'C-sol_'//trim(int2char(iplot))
     .                              //' after coarse solve'
            dbg_desc(iplot+2*neq) = 'C-res_'//trim(int2char(iplot))
     .                              //' after coarse solve'
          enddo

          call createDrawInCfile(3*neq,'coarse.bin','Coarse grid'
     .          ,'t','x','y',dbg_desc,'-c -X0 -L57'
     .          ,'drawcrse.in')

        endif

c     Relax error/solution on grid number igr/igrid (find new xx)

        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGmatvec(mgctx,0,xx (isig:isig+nn-1)
     .                         ,wrk(isig:isig+nn-1),igr,bcnd,matvec)
          call MGvecadd(mgctx,igr,-1d0,wrk(isig:isig+nn-1)
     .                           , 1d0,yy (isig:isig+nn-1))
          call MGplot(mgctx,igr,wrk,0,'fineres.bin')
        endif

        call smooth(igr)

        !Debugging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGplot(mgctx,igr,xx,0,'finesol.bin')
        endif

c     Evaluate residual (ie wrk = yy - A xx = yy - wrk )

        wrk(isig:isig+nn-1) = 0d0

        call MGmatvec(mgctx,0,xx (isig:isig+nn-1)
     .                       ,wrk(isig:isig+nn-1),igr,bcnd,matvec)

        call MGvecadd(mgctx,igr,-1d0,wrk(isig:isig+nn-1)
     .                         , 1d0,yy (isig:isig+nn-1))

        !Debuging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGplot(mgctx,igr,wrk,1,'fineres.bin')
        endif

c     Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) to a coarser grid

        call crestrict(mgctx,wrk(isig :isig +nn -1)
     .                      ,yy (isigc:isigc+nnc-1)
     .                ,orderres,igr,res_vol_wgt)

        !Debugging plots
        if (igc == igmax.and.mg_debug) then
          write (*,*) 'Plotting here at level',igc
          call MGplot(mgctx,igc,yy,0,'coarse.bin')
        endif

c     Initialize solution on coarse grid

        xx(isigc:isigc+nnc-1) = 0d0

c     If on coarsest grid, solve for error, else descend a grid level

cc#if defined(petsc)
cc        if ((np > 1)
cc     .      .and. (    MGgrid%ngrdx == igc    !Grid-level based criterion
cc     .             .or.MGgrid%ngrdy == igc
cc     .             .or.MGgrid%ngrdz == igc)
cccc     .      .and. (    MGgrid%nxgl(igc) >= 32 !Coarse grid size criterion
cccc     .             .or.MGgrid%nygl(igc) >= 32
cccc     .             .or.MGgrid%nzgl(igc) >= 32
cccc     .             .or.MGgrid%nxgl(igc)
cccc     .                *MGgrid%nygl(igc)
cccc     .                *MGgrid%nzgl(igc) >= 1024)
cc     .  ) then
cc
cc          write (*,*) 'DIAG -- mcycle at coarseGridSolve: proc=',my_rank
cc
cc          call coarseGridSolve(igc)     !Perform parallel coarse grid solve
cc
cc        elseif (igc.eq.igmax) then
cc#else
        if (igc.eq.igmax) then
cc#endif
          call coarseSolve(igc)         !Perform coarse grid solve
        else
          do imu = 1,mu
            call mcycle(igc)
          enddo
        endif

        !Debugging plots
        if (igc == igmax.and.mg_debug) then
          write (*,*) 'Plotting here at level',igc
          call MGplot(mgctx,igc,xx,1,'coarse.bin')

          call MGmatvec(mgctx,0,xx (isigc:isigc+nnc-1)
     .                         ,wrk(isigc:isigc+nnc-1),igc,bcnd,matvec)
          call MGvecadd(mgctx,igc,-1d0,wrk(isigc:isigc+nnc-1)
     .                           , 1d0,yy (isigc:isigc+nnc-1))
          call MGplot(mgctx,igc,wrk,1,'coarse.bin')
        endif

c     Cycle back up to grid igr updating errors (xx)
c     with fixed R.H.S. (yy)

c     Prolong error (wrk = P * xx_1) to a finer grid

        call cprolong(mgctx,xx (isigc:isigc+nnc-1)
     .                     ,wrk(isig :isig +nn -1)
     .               ,orderprol,igc,bcnd)

c     Update existing solution on grid igr (i.e. xx_igr): xx_igr = xx_igr + wrk 

        call MGvecadd(mgctx,igr,1d0,xx (isig:isig+nn-1)
     .                         ,1d0,wrk(isig:isig+nn-1))

        !Debugging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGplot(mgctx,igr,wrk,1,'finesol.bin')
          call MGplot(mgctx,igr,xx ,1,'finesol.bin')

          call MGmatvec(mgctx,0,xx (isig:isig+nn-1)
     .                         ,wrk(isig:isig+nn-1),igr,bcnd,matvec)
          call MGvecadd(mgctx,igr,-1d0,wrk(isig:isig+nn-1)
     .                           , 1d0,yy (isig:isig+nn-1))
          call MGplot(mgctx,igr,wrk,1,'fineres.bin')
        endif

c     Relax updated solution on igr (i.e. xx_igr)

        call smooth(igr)

        !Debugging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr

          call MGmatvec(mgctx,0,xx (isig:isig+nn-1)
     .                         ,wrk(isig:isig+nn-1),igr,bcnd,matvec)
          call MGvecadd(mgctx,igr,-1d0,wrk(isig:isig+nn-1)
     .                           , 1d0,yy (isig:isig+nn-1))
          call MGplot(mgctx,igr,wrk,1,'fineres.bin')
        endif

c     End program

      end subroutine mcycle

c     smooth
c     #####################################################################
      recursive subroutine smooth(igr)

c     ---------------------------------------------------------------------
c     Performs MG smoothing at grid igr. Smoothing may be point or
c     plane/line smoothing, depending on external variable line_relax.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igr

c     Local variables

        integer    :: nn,isig,depth1
        type (solver_unit) :: solver

c     Begin program

        if (outc.ge.1.and.prnt) write (*,*) 'Grid Level',igr

        if (.not.line_relax) then
          nn   = ntotv (igr)
          isig = istart(igr)

          depth1 = depth + 1

          !Inherit relevant previous solver info
          call readSolverHierarchy (solver,depth1)
          solver%options%diag => diag
          solver%options%mgctx=> mgctx
          call writeSolverHierarchy(solver,depth1)

          !Call solver
          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                  ,guess2,outc,depth1)
        else
          call linesmooth(igr)
        endif

c     End program

      end subroutine smooth

c     linesmooth
c     #####################################################################
      recursive subroutine linesmooth(igr)

c     ---------------------------------------------------------------------
c     This routine performs plane/line smoothing instead of point
c     smoothing.  This is done recursively, by calling 2D MG for planes
c     and 1D MG for lines. When reducing the dimensionality of MG,
c     vectors remain of the same length as in the original grid, and
c     only matvec and smoothing operations use a decreased
c     dimensionality. This allows one to employ the SAME matvec and
c     smoothing routines as for the original grid, with minor
c     modifications that restrict operations to the proper plane/line.
c     
c     As vector dimensions are the same when going from 3d MG to 2D MG
c     and 1D MG, MG pointers are not reallocated (and are consistent
c     with the original grid, defined in MGgrid, and stored in
c     mg_ctx_temp).  However, the grid definition
c     (nxv,nyv,nzv,iline,jline,kline) contained in mg_ctx_temp is
c     "collapsed", and this is the way MG knows it should operate on a
c     restricted subspace of mesh points. Therefore, node positioning on
c     the grid follows MGgrid, whereas matvec and smoothing
c     operations follow mg_ctx_temp.
c     -----------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igr

c     Local variables

        integer    :: nn,isig,depth1,i,j,k,iii,ig,iter,nsweep
     .               ,imin,imax,jmin,jmax,kmin,kmax
        integer    :: nxl,nyl,nzl

        real(8)    :: mag,mag0,omega
        real(8),allocatable,dimension(:) :: dx,rr

        type(grid_mg_def),pointer :: line_mg_grid

        type(mg_ctx),pointer :: line_mg_ctx,mg_ctx_sv

        logical :: fpointers

c     Begin program

        nn   = ntotv (igr)
        isig = istart(igr)

        !Line smoothing dimensions
        nxl = mgctx%nxv(igr)
        nyl = mgctx%nyv(igr)
        nzl = mgctx%nzv(igr)

        nsweep = line_nsweep

        omega  = line_omega

        nullify(line_mg_grid)
        nullify(line_mg_ctx)

c     Schedule planes/lines

c       Check if at least 2D problem (based on existing configuration)
        if (    (nxl > 1 .and. nyl > 1)
     .      .or.(nxl > 1 .and. nzl > 1)
     .      .or.(nyl > 1 .and. nzl > 1) ) then

c         Initialize quantities

          allocate(rr(nn),dx(nn))

          i = 1
          j = 1
          k = 1

c         Find initial residual

          call MGmatvec(mgctx,0,xx(isig:isig+nn-1),rr,igr,bcnd,matvec)
          call MGvecadd(mgctx,igr,-1d0,rr,1d0,yy(isig:isig+nn-1)) !rr = yy(isig:isig+nn-1)-rr

          if (outc > 0) then
            mag0 = sqrt(MGdot(mgctx,igr,rr,rr))
            if (outc > 1.and.prnt) write (*,*)
            if (prnt) write (*,7) mag0
          endif

c         Planes/Lines in X-direction

          call allocMGPointers(neq,fpointers,line_mg_ctx)

          if (nxl > 1 .and. line_x) then

            !Collapse MG grid levels
            line_mg_ctx%nxv = 1         !Collapse grid in x-direction
            line_mg_ctx%mg_ratio_x = 1  !Collapse grid in x-direction

            do iter=1,nsweep

              if (outc > 1.and.prnt) then
                write (*,*)
                write (*,*) '******** X-relax sweep #',iter,'********'
              endif

              do i=1,nxl

                !Define plane/line for subsequent MG solvers at all grids
                line_mg_ctx%iline(igr) = i
                do ig = igr+1,line_mg_ctx%ngrid
                  line_mg_ctx%iline(ig) = ( line_mg_ctx%iline(ig-1)
     .                                 +mg_ctx_temp%mg_ratio_x(ig-1)-1 )
     .                                 /mg_ctx_temp%mg_ratio_x(ig-1)
                enddo

cc                write (*,*) 'DIAG -- linesmooth'
cc     .                     ,line_mg_ctx%iline(igr:line_mg_ctx%ngrid)

                !Output info
                if (outc > 1.and.prnt) then
                  write (*,*)
                  write (*,'(a,i4,a,i4,a,i4)')
     .                 'PLANE SOLVE: i=',i,'; ny x nz:',nyl,'x',nzl
                endif

                dx = 0d0

                !Save parent MG pointers
                mg_ctx_sv => mg_ctx_temp

                !Solve plane/line (find solution update dx)
                call lineSolver(line_mg_ctx,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call MGvecadd(line_mg_ctx,igr,1d0,xx(isig:isig+nn-1)
     .                       ,omega,dx)

                !Recover grid pointers
                mg_ctx_temp => mg_ctx_sv
                nullify(mg_ctx_sv)

              enddo

              !Update residual
              call MGmatvec(mgctx,0,xx(isig:isig+nn-1),rr,igr,bcnd
     .                     ,matvec)
              call MGvecadd(mgctx,igr,-1d0,rr,1d0,yy(isig:isig+nn-1)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(MGdot(mgctx,igr,rr,rr))
              if (outc > 1.and.prnt) write (*,5)
              if (prnt) write (*,10) mag,mag/mag0
              if (outc > 1.and.prnt) write (*,6)
              mag0 = mag
            endif

            !Recover configuration
            line_mg_ctx%iline = 0
            line_mg_ctx%nxv = mg_ctx_temp%nxv
            line_mg_ctx%mg_ratio_x = mg_ctx_temp%mg_ratio_x

          endif

c         Planes/Lines in y-direction

          if (nyl > 1 .and. line_y) then

cc            !Update residual
cc            call MGmatvec(0,neq,nn,xx(isig),rr,igr,bcnd,matvec)
cc            call MGvecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            !Collapse MG grid levels
            line_mg_ctx%nyv = 1         !Collapse grid in y-direction
            line_mg_ctx%mg_ratio_y = 1  !Collapse grid in y-direction

            do iter=1,nsweep

              if (outc > 1.and.prnt) then
                write (*,*)
                write (*,*) '******** Y-relax sweep #',iter,'********'
              endif

              do j=1,nyl

                !Define plane/line for subsequent MG solvers at all grids
                line_mg_ctx%jline(igr) = j
                do ig = igr+1,line_mg_ctx%ngrid
                  line_mg_ctx%jline(ig) = ( line_mg_ctx%jline(ig-1)
     .                                 +mg_ctx_temp%mg_ratio_y(ig-1)-1 )
     .                                 /mg_ctx_temp%mg_ratio_y(ig-1)
                enddo

cc                write (*,*) 'DIAG -- linesmooth'
cc     .                     ,line_mg_ctx%jline(igr:line_mg_ctx%ngrid)

                !Output info
                if (outc > 1.and.prnt) then
                  write (*,*)
                  write (*,'(a,i4,a,i4,a,i4)')
     .                 'PLANE SOLVE: j=',j,'; nx x nz:',nxl,'x',nzl
                endif

                dx = 0d0

                !Save parent MG pointers
                mg_ctx_sv => mg_ctx_temp

                !Solve plane/line (find solution update dx)
                call lineSolver(line_mg_ctx,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call MGvecadd(line_mg_ctx,igr,1d0,xx(isig:isig+nn-1)
     .                       ,omega,dx)

                !Recover grid pointers
                mg_ctx_temp => mg_ctx_sv
                nullify(mg_ctx_sv)

              enddo

              !Update residual
              call MGmatvec(mgctx,0,xx(isig:isig+nn-1),rr,igr,bcnd
     .                     ,matvec)
              call MGvecadd(mgctx,igr,-1d0,rr,1d0,yy(isig:isig+nn-1)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(MGdot(mgctx,igr,rr,rr))
              if (outc > 1.and.prnt) write (*,5)
              if (prnt) write (*,20) mag,mag/mag0
              if (outc > 1.and.prnt) write (*,6)
              mag0 = mag
            endif

            !Recover configuration
            line_mg_ctx%jline = 0
            line_mg_ctx%nyv = mg_ctx_temp%nyv
            line_mg_ctx%mg_ratio_y = mg_ctx_temp%mg_ratio_y

          endif

c         Planes/Lines in z-direction

          if (nzl > 1 .and. line_z) then

            !Collapse MG grid levels
            line_mg_ctx%nzv = 1         !Collapse grid in z-direction
            line_mg_ctx%mg_ratio_z = 1  !Collapse grid in z-direction

            do iter=1,nsweep

              if (outc > 1.and.prnt) then
                write (*,*)
                write (*,*) '******** Z-relax sweep #',iter,'********'
              endif

              do k=1,nzl

                !Define plane/line for subsequent MG solvers at all grids
                line_mg_ctx%kline(igr) = k
                do ig = igr+1,line_mg_ctx%ngrid
                  line_mg_ctx%kline(ig) = ( line_mg_ctx%kline(ig-1)
     .                                 +mg_ctx_temp%mg_ratio_z(ig-1)-1 )
     .                                 /mg_ctx_temp%mg_ratio_z(ig-1)
                enddo

cc                write (*,*) 'DIAG -- linesmooth'
cc     .                     ,line_mg_ctx%kline(igr:line_mg_ctx%ngrid)

                !Output info
                if (outc > 1.and.prnt) then
                  write (*,*)
                  write (*,'(a,i4,a,i4,a,i4)')
     .                 'PLANE SOLVE: k=',k,'; nx x ny:',nxl,'x',nyl
                endif

                dx = 0d0

                !Save parent MG pointers
                mg_ctx_sv => mg_ctx_temp

                !Solve plane/line (find solution update dx)
                call lineSolver(line_mg_ctx,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call MGvecadd(line_mg_ctx,igr,1d0,xx(isig:isig+nn-1)
     .                                       ,omega,dx)

                !Recover grid pointers
                mg_ctx_temp => mg_ctx_sv
                nullify(mg_ctx_sv)

              enddo

              !Update residual
              call MGmatvec(mgctx,0,xx(isig:isig+nn-1),rr,igr,bcnd
     .                     ,matvec)
              call MGvecadd(mgctx,igr,-1d0,rr,1d0,yy(isig:isig+nn-1)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(MGdot(mgctx,igr,rr,rr))
              if (outc > 1.and.prnt) write (*,5)
              if (prnt) write (*,30) mag,mag/mag0
              if (outc > 1.and.prnt) write (*,6)
              mag0 = mag
            endif

            !Recover configuration
            line_mg_ctx%kline = 0
            line_mg_ctx%nzv = mg_ctx_temp%nzv
            line_mg_ctx%mg_ratio_z = mg_ctx_temp%mg_ratio_z

          endif

c         Deallocate variables

          call deallocMGPointers(fpointers,line_mg_ctx)

          deallocate(rr,dx)

c       Else smooth 1D problem (for MG line smoother)
        else

          if (outc > 1.and.prnt) then
            write (*,*)
            if (nxl > 1) then
              write (*,*) '1D LINE SOLVE on j =',mg_ctx_temp%jline(igr)
     .               ,', k =',mg_ctx_temp%kline(igr),' **** i = 1 :',nxl
            elseif (nyl > 1) then
              write (*,*) '1D LINE SOLVE on i =',mg_ctx_temp%iline(igr)
     .               ,', k =',mg_ctx_temp%kline(igr),' **** j = 1 :',nyl
            else
              write (*,*) '1D LINE SOLVE on i =',mg_ctx_temp%iline(igr)
     .               ,', j =',mg_ctx_temp%jline(igr),' **** k = 1 :',nzl
            endif
          endif

          depth1 = depth + 1

          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                  ,guess2,outc,depth1)

        endif

c     End program

 5    format(/,' *****************************************************')
 6    format(  ' *****************************************************')
 7    format('  Initial residual:',1p,1e10.2)

 10   format('  X-plane relax. Residual:',1p,1e10.2,'; Ratio:',1e10.2)

 20   format('  Y-plane relax. Residual:',1p,1e10.2,'; Ratio:',1e10.2)

 30   format('  Z-plane relax. Residual:',1p,1e10.2,'; Ratio:',1e10.2)

      end subroutine linesmooth
        
c     lineSolver
c     #####################################################################
      recursive subroutine lineSolver(mgctx,b,x,igr,guess)

c     ---------------------------------------------------------------------
c     This routine performs a recursive MG solve on selected planes/lines.
c     In call:
c       * nn: total grid size at grid level igr
c       * mgctx: MG grid definition structure that defines local grid patch
c       * b: rhs vector
c       * x: solution vector
c       * igr: grid level
c       * guess: whether initial guess is provided (guess=1) or not (guess=0).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer        :: igr,guess
        real(8)        :: b(:),x(:)
cc        type(grid_mg_def),pointer :: mgctx
        type(mg_ctx),pointer :: mgctx
        logical        :: prelax

c     Local variables

        integer        :: nxl,nyl,nzl,nn

c     Begin program

        nn = size(b)

        !Line smoothing dimensions
        nxl = mgctx%nxv(igr)
        nyl = mgctx%nyv(igr)
        nzl = mgctx%nzv(igr)

c     Schedule solvers

ccc       Recursive MG (at least 2D or if 1D MG is wanted)
cc        if  (    (nxl > 1 .and. nyl > 1)
cc     .       .or.(nxl > 1 .and. nzl > 1)
cc     .       .or.(nyl > 1 .and. nzl > 1)) then
cc
cc          call lineMGsolver(nn,mgctx,b,x,igr,guess)
cc
ccc       GMRES 1D solver (matvec HARDWIRED to v_mtvc for now)
cccc        elseif (line_solve == "gm") then
cccc
cccc          call lineGMsolver(nn,mgctx,b,x,igr,guess)
cc
ccc       MG 1D solver
cc        elseif (line_solve == "mg") then
cc
cc          call lineMGsolver(nn,mgctx,b,x,igr,guess,oned_solve=.true.)
cc
ccc       GS 1D solver
cc        elseif (line_solve == "gs") then
cc
cc          call lineGSsolver(nn,mgctx,b,x,igr,guess)
cc
ccc       JB 1D solver
cc        elseif (line_solve == "jb") then

          call lineJBsolver(mgctx,b,x,igr,guess)

cc        else
cc
cc          call pstop('lineSolver','Line solver undefined')
cc
cc        endif

c     End program

      end subroutine lineSolver

ccc     lineGMsolver
ccc     #####################################################################
cc      recursive subroutine lineGMsolver(nn,mgctx,b,x,igr,guess)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine performs a GMRES solve on selected lines.
ccc     In call:
ccc       * nn: total grid size at grid level igr
ccc       * mgctx: MG grid definition structure that defines local grid patch
ccc       * b: rhs vector
ccc       * x: solution vector
ccc       * igr: grid level
ccc       * guess: whether initial guess is provided (guess=1) or not (guess=0).
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer        :: igr,nn,guess
cc        real(8)        :: b(nn),x(nn)
cc        type(grid_mg_def),pointer :: mgctx
cc        logical        :: prelax
cc
ccc     Local variables
cc
cc        integer    :: nxl,nyl,nzl,nnl,depth1
cc        integer    :: imin,imax,jmin,jmax,kmin,kmax
cc     .               ,i,j,k,iii,il,jl,kl,iil,ieq
cc
cc        real(8),allocatable,dimension(:) :: bb,xx
cc
cc        type (solver_options) :: options
cc
ccc     Begin program
cc
ccc     Find local problem size
cc
cc        nxl = mgctx%nxv(igr)
cc        nyl = mgctx%nyv(igr)
cc        nzl = mgctx%nzv(igr)
cc
cc        nnl = nxl*nyl*nzl*neq
cc
cc        allocate(bb(nnl),xx(nnl))
cc
cc        xx = 0d0
cc
ccc     Map global to local vectors
cc
cc        !Use global grid limits for this
cc        call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
cc     .             ,imin,imax,jmin,jmax,kmin,kmax)
cc
cc        do k=kmin,kmax
cc          do j=jmin,jmax
cc            do i=imin,imax
cc              do ieq=1,neq
cc                !Global vector location
cc                iii=getMGvcomp(i,j,k,nxv(igr),nyv(igr),nzv(igr),1
cc     .                        ,ieq,neq)
cc                !Local indices
cc                il = i - imin + 1
cc                jl = j - jmin + 1
cc                kl = k - kmin + 1
cc                iil=getMGvcomp(il,jl,kl,nxl,nyl,nzl,1,ieq,neq)
cc
cc                bb(iil) = b(iii)
cc                if (guess == 1) xx(iil) = x(iii)
cc              enddo
cc            enddo
cc          enddo
cc        enddo
cc
ccc     Configure GMRES solve (need to initialize ALL relevant options)
cc
cc        options%sym_test        = .false.
cc
cc        options%tol             = line_tol
cc
cc        options%krylov_subspace = nnl
cc        options%iter            = nnl
cc        options%stp_test        = 1 
cc
ccc     Call GMRES (proxy routine lmtvc is defined in mg_ops module)
cc
cc        depth1 =depth + 1
cc        call gm(neq,nnl,bb,xx,lmtvc,options,igr,bcnd,guess,outc-1
cc     .         ,depth1)
cc
ccc     Map solution to global grid
cc
cc        do k=kmin,kmax
cc          do j=jmin,jmax
cc            do i=imin,imax
cc              do ieq=1,neq
cc                !Global vector location
cc                iii=getMGvcomp(i,j,k,nxv(igr),nyv(igr),nzv(igr),1
cc     .                        ,ieq,neq)
cc                !Local indices
cc                il = i - imin + 1
cc                jl = j - jmin + 1
cc                kl = k - kmin + 1
cc                iil=getMGvcomp(il,jl,kl,nxl,nyl,nzl,1,ieq,neq)
cc
cc                x(iii) = xx(iil)
cc              enddo
cc            enddo
cc          enddo
cc        enddo
cc
ccc     Deallocate memory
cc
cc        deallocate(bb,xx)
cc
cc      end subroutine lineGMsolver

ccc     lineGSsolver
ccc     #####################################################################
cc      recursive subroutine lineGSsolver(nn,mgctx,b,x,igr,guess)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine performs a recursive MG solve on selected planes/lines.
ccc     In call:
ccc       * nn: total grid size at grid level igr
ccc       * mgctx: MG grid definition structure that defines local grid patch
ccc       * b: rhs vector
ccc       * x: solution vector
ccc       * igr: grid level
ccc       * guess: whether initial guess is provided (guess=1) or not (guess=0).
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer        :: igr,nn,guess
cc        real(8)        :: b(nn),x(nn)
cccc        type(grid_mg_def),pointer :: mgctx
cc        type(mg_ctx),pointer :: mgctx
cc        logical        :: prelax
cc
ccc     Local variables
cc
cc        integer        :: istart_sv (size(istart) )
cc     .                   ,istartp_sv(size(istartp))
cc     .                   ,istartb_sv(size(istartb))
cc
cc        real(8),allocatable,dimension(:,:) :: old_diag
cc
cc        type (solver_options) :: options
cc
cc        logical               :: fpointers
cc
ccc     Begin program
cc
ccc     Save parent MG grid configuration (shared with filial MG via mg_ops module)
cc
cccc        MGgrid_sv%ngrdx = 0
cccc        MGgrid_sv%ngrdy = 0
cccc        MGgrid_sv%ngrdz = 0
cccc        MGgrid_sv%ngrid = 0
cccc
cccc        !Save pointers
cccc        MGgrid_sv  = MGgrid
cc        istart_sv  = istart
cc        istartp_sv = istartp
cc        istartb_sv = istartb
cc
cc        !Save diagonal
cc        allocate(old_diag(size(diag,1),size(diag,2)))
cc        old_diag = diag
cc        deallocate(diag)
cc
ccc     'Prime' grid structure variable to pass to filial MG call
cc
cc        options%mg_grid_def%ngrdx = 0
cc        options%mg_grid_def%ngrdy = 0
cc        options%mg_grid_def%ngrdz = 0
cc        options%mg_grid_def%ngrid = 0
cc
cccc        options%mg_grid_def = MGgrid_sv
cc        call equateGridStructure(options%mg_grid_def,MGgrid_sv)
cc
ccc     Shift MG pointers to current grid level igr
cc
cc        call transferMGPointers(igr)
cc
ccc     Allocate new diagonal and transfer elements
cc
cc        allocate(diag(size(old_diag,1)
cc     .               ,size(old_diag,2)-istart_sv(igr)+1))
cc
cc        diag = old_diag(:,istart_sv(igr):size(old_diag,2))
cc
ccc     Configure recursive plane/line GS solve
cc
cc        options%iter                        = 100
cc        options%ncolors                     = 4
cc        options%omega                       = 1.0           
cc
cc        options%tol                         = line_tol
cc        options%vcyc                        = line_vcyc
cc        options%mg_gmin                     = igridmin
cc        options%orderres                    = orderres
cc        options%orderprol                   = orderprol
cc        options%mg_mu                       = mu
cc        options%vol_res                     = volf        
cc
cc        options%mg_coarse_solver_depth      = line_crse_depth
cc        options%mg_line_relax               = .false.
cc        options%mg_line_nsweep              = line_nsweep
cc        options%mg_line_vcyc                = line_vcyc
cc        options%mg_line_tol                 = line_tol
cc        options%mg_line_omega               = line_omega
cc        options%mg_line_coarse_solver_depth = line_crse_depth
cccc        options%mg_grid_def                 = mg_grid
cc        call equateGridStructure(options%mg_grid_def,mg_grid)
cc
cc        options%vertex_based_relax          = vbr_mg
cc        options%fdiag                       = fdiag
cc
ccc     Call plane/line GS
cc
cc        call gs(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc-1,depth)
cc
ccc     Recover configuration for parent level MG
cc
cc        !Recover pointers
cc        istart  = istart_sv
cc        istartp = istartp_sv
cc        istartb = istartb_sv
cccc        MGgrid  = mg_grid
cc        call equateGridStructure(MGgrid,mg_grid)
cc
cc        !Recover diagonal
cc        deallocate(diag)
cc        allocate(diag(size(old_diag,1),size(old_diag,2)))
cc        diag = old_diag
cc        deallocate(old_diag)
cc
ccc     Deallocate memory
cc
cc        call deallocateGridStructure(options%mg_grid_def)
cc
cc      end subroutine lineGSsolver

c     lineJBsolver
c     #####################################################################
      recursive subroutine lineJBsolver(mgctx,b,x,igr,guess)

c     ---------------------------------------------------------------------
c     This routine performs a recursive JB solve on selected planes/lines.
c     In call:
c       * nn: total grid size at grid level igr
c       * mgctx: MG grid definition structure that defines local grid patch
c       * b: rhs vector
c       * x: solution vector
c       * igr: grid level
c       * guess: whether initial guess is provided (guess=1) or not (guess=0).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer        :: igr,guess
        real(8)        :: b(:),x(:)
        type(mg_ctx),pointer :: mgctx
        logical        :: prelax

c     Local variables

        integer :: nn
        type (solver_options) :: options

c     Begin program

        nn = size(b)

        nullify(options%mg_grid_def)

c     Configure recursive plane/line JB solve

        options%iter                        = 2
        options%omega                       = 0.7

        options%tol                         = line_tol
        options%vcyc                        = line_vcyc
        options%mg_gmin                     = igridmin
        options%orderres                    = orderres
        options%orderprol                   = orderprol
        options%mg_mu                       = mu
cc        options%vol_res                     = volf        

        options%mg_coarse_solver_depth      = line_crse_depth
        options%mg_line_relax               = .false.
        options%mg_line_nsweep              = line_nsweep
        options%mg_line_vcyc                = line_vcyc
        options%mg_line_tol                 = line_tol
        options%mg_line_omega               = line_omega
        options%mg_line_coarse_solver_depth = line_crse_depth

        options%mg_vertex_relax             = vbr_mg

        options%mgctx                       => mgctx
        options%diag                        => diag

c     Call plane/line JB

        call jb(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc-1,depth)

c     Deallocate memory

        nullify(options%mg_grid_def,options%mgctx)
        nullify(options%diag)

      end subroutine lineJBsolver

ccc     lineMGsolver
ccc     #####################################################################
cc      recursive subroutine lineMGsolver(nn,mg_grid,b,x,igr,guess
cc     .                                 ,oned_solve)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine performs a recursive MG solve on selected planes/lines.
ccc     In call:
ccc       * nn: total grid size at grid level igr
ccc       * mg_grid: MG grid definition structure that defines local grid patch
ccc       * b: rhs vector
ccc       * x: solution vector
ccc       * igr: grid level
ccc       * guess: whether initial guess is provided (guess=1) or not (guess=0).
ccc       * 1d_solve (optional):whether we are doing 1D MG.
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer        :: igr,nn,guess
cc        real(8)        :: b(nn),x(nn)
cc        type(grid_mg_def),pointer :: mg_grid
cc        logical,optional :: oned_solve
cc
ccc     Local variables
cc
cc        integer        :: istart_sv (size(istart) )
cc     .                   ,istartp_sv(size(istartp))
cc     .                   ,istartb_sv(size(istartb))
cc
cc        real(8),allocatable,dimension(:,:) :: old_diag
cc
cc        type (solver_options) :: options
cc
cc        logical               :: fpointers,oned_slv
cc
ccc     Begin program
cc
cc        if (PRESENT(oned_solve)) then
cc          oned_slv = oned_solve
cc        else
cc          oned_slv = .false.
cc        endif
cc
ccc     Save parent MG grid configuration (shared with filial MG via mg_ops module)
cc
cc        !Save MG pointers
cc        istart_sv  = istart
cc        istartp_sv = istartp
cc        istartb_sv = istartb
cc
cc        !Save diagonal
cc        allocate(old_diag(size(diag,1),size(diag,2)))
cc        old_diag = diag
cc        deallocate(diag)
cc
ccc     'Prime' grid structure variable to pass to filial MG call
cc
cccc        options%mg_grid_def%ngrdx = 0
cccc        options%mg_grid_def%ngrdy = 0
cccc        options%mg_grid_def%ngrdz = 0
cccc        options%mg_grid_def%ngrid = 0
cc        nullify(options%mg_grid_def)
cc
cccc        call equateGridStructure(options%mg_grid_def,MGgrid_sv)
cc
ccc     Shift MG pointers to current grid level igr
cc
cc        call transferMGPointers(igr)
cc
ccc     Allocate new diagonal and transfer elements
cc
cc        allocate(diag(size(old_diag,1)
cc     .               ,size(old_diag,2)-istart_sv(igr)+1))
cc
cc        diag = old_diag(:,istart_sv(igr):size(old_diag,2))
cc
ccc     Configure recursive plane/line MG solve
cc
cc        options%tol                         = line_tol
cc        options%vcyc                        = line_vcyc
cc        options%mg_gmin                     = igridmin
cc        if (oned_slv) then
cc          options%orderres                  = 0
cc          options%orderprol                 = 0
cc        else
cc          options%orderres                  = orderres
cc          options%orderprol                 = orderprol
cc        endif
cc        options%mg_mu                       = mu
cc        options%vol_res                     = volf        
cc
cc        options%mg_coarse_solver_depth      = line_crse_depth
cc        options%mg_line_relax               = .true.
cc        options%mg_line_nsweep              = line_nsweep
cc        options%mg_line_vcyc                = line_vcyc
cc        options%mg_line_tol                 = line_tol
cc        options%mg_line_omega               = line_omega
cc        options%mg_line_coarse_solver_depth = line_crse_depth
cc        options%mg_line_x                   = line_x
cc        options%mg_line_y                   = line_y
cc        options%mg_line_z                   = line_z
cc        options%mg_grid_def                 => mg_grid
cccc        call equateGridStructure(options%mg_grid_def,mg_grid)
cc
cc        options%vertex_based_relax          = vbr_mg
cc        options%fdiag                       = fdiag
cccc        options%diag                        => new_diag
cc
ccc     Call plane/line MG
cc
cc        call mg(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc-1,depth)
cc
ccc     Recover configuration for parent level MG
cc
cc        !Recover pointers
cc        istart  = istart_sv
cc        istartp = istartp_sv
cc        istartb = istartb_sv
cccc        MGgrid  = mg_grid
cc        call equateGridStructure(MGgrid,mg_grid)
cccc        MGgrid  = MGgrid_sv
cccc        call deallocateGridStructure(MGgrid_sv)
cc
cc        !Recover diagonal
cc        deallocate(diag)
cc        allocate(diag(size(old_diag,1),size(old_diag,2)))
cc        diag = old_diag
cc        deallocate(old_diag)
cc
ccc     Deallocate memory
cc
cc        call deallocateGridStructure(options%mg_grid_def)
cc
cc      end subroutine lineMGsolver

c     transferMGPointers
c     #################################################################
      subroutine transferMGPointers(igr)

c     -----------------------------------------------------------------
c     Transfers MG pointers from upper grid level to subsequent grid
c     levels. This is done by shifting pointer arrays by one grid,
c     removing the finest grid level each time. The MG grid level
c     definition mg_grid_sv is used as a reference.
c     -----------------------------------------------------------------

        implicit none             !For safe fortran

c     Call variables

        integer        :: igr

c     Local variables

        integer        :: ig

c     Begin program

c     If not at finest grid level, shift MG pointers

        if (igr > 1) then

          do ig=ngrid,igr+1,-1
            istartp(ig) = istartp(ig)-istartp(igr)+1
            istart (ig) = istart (ig)-istart (igr)+1
            istartb(ig) = istartb(ig)-istartb(igr)+1
          enddo
          istartp(igr) = 1
          istart (igr) = 1
          istartb(igr) = 1

        endif

c     End program

      end subroutine transferMGPointers

c     saveMGgrid
c     #####################################################################
      recursive subroutine saveMGgrid(mg_grid)

c     ---------------------------------------------------------------------
c     This routine saves current MG grid configuration for recursive MG
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: mg_grid

c     Local variables

c     Begin program

        MGgrid_sv  => MGgrid
        MGgrid     => mg_grid

c     End program

      end subroutine saveMGgrid

c     restoreMGgrid
c     #####################################################################
      subroutine restoreMGgrid

c     ---------------------------------------------------------------------
c     This routine saves current MG grid configuration for recursive MG
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        MGgrid  => MGgrid_sv
        nullify(MGgrid_sv)

c     End program

      end subroutine restoreMGgrid

c     coarseSolve
c     ###################################################################
      recursive subroutine coarseSolve(igr)

c     -------------------------------------------------------------------
c     Performs solve at coarsest grid level igr.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: igr

c     Local variables

        integer :: nn,isig,depth1

c     Begin program

        if (outc.ge.1.and.prnt) write (*,*) 'Grid Level',igr

cc#if defined(petsc)
cccc        if ((np > 1) .and.
cccc     .      (    MGgrid%ngrdx == igr
cccc     .       .or.MGgrid%ngrdy == igr
cccc     .       .or.MGgrid%ngrdz == igr)) then   !Perform coarse grid solve
cccc
cccc        else
cc
cc          nn   = ntotv(igr)
cc          isig = istart(igr)
cc
cc          if (crsedpth /= 0) then
cc            depth1 = crsedpth
cc          else
cc            depth1 = depth + 1
cc          endif
cc
cc          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
cc     .                ,guess2,outc,depth1)
cc
cccc        endif
cc#else
        nn   = ntotv(igr)
        isig = istart(igr)

        if (crsedpth /= 0) then
          depth1 = crsedpth
        else
          depth1 = depth + 1
        endif

        call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                ,guess2,outc,depth1)
cc#endif

      end subroutine coarseSolve

ccc     coarseGridSolve
ccc     #####################################################################
cc      recursive subroutine coarseGridSolve(igr)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine performs a recursive MG solve on parallel coarse grid.
ccc     In call:
ccc       * igr: grid level
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer        :: igr
cc
ccc     Local variables
cc
cc        integer        :: nn,isig,nx,ny,nz,nc,i,j,k,il,jl,kl
cc     .                   ,ii,iig,guess,istart_sv(size(istart))
cc     .                   ,ig_cs,ig_sv,isig_cs,isig_sv,nsv,ncs
cc
cc        real(8),allocatable,dimension(:,:) :: diag_sv
cc
cc        real(8),allocatable,dimension(:)   :: b_global,x_global
cc
cc        type(grid_mg_def),pointer :: gparams_sv
cc
cc        type (solver_options) :: options
cc
cc        logical               :: fpointers
cc
ccc     Begin program
cc
cc        nn   = ntotv(igr)
cc        isig = istart(igr)
cc
ccc     Save parent MG grid configuration (shared with filial MG via mg_ops module)
cc
cc        !Save pointers
cc        MGgrid_sv%ngrdx = 0
cc        MGgrid_sv%ngrdy = 0
cc        MGgrid_sv%ngrdz = 0
cc        MGgrid_sv%ngrid = 0
cccc        MGgrid_sv  = MGgrid
cc        call equateGridStructure(MGgrid_sv,MGgrid)
cc
cc        call deallocateGridStructure(MGgrid)
cc
cc        gparams_sv%ngrdx = 0
cc        gparams_sv%ngrdy = 0
cc        gparams_sv%ngrdz = 0
cc        gparams_sv%ngrid = 0
cccc        gparams_sv  = grid_params
cc        call equateGridStructure(gparams_sv,grid_params)
cc
cc        call deallocateGridStructure(grid_params)
cc
cc        istart_sv  = istart
cc
cc        !Save diagonal
cc        allocate(diag_sv(size(diag,1),size(diag,2)))
cc        diag_sv = diag
cc        deallocate(diag)
cc
ccc     Find coarse global grid dimensions
cc
cc        nx = MGgrid_sv%nxgl(igr)
cc        ny = MGgrid_sv%nygl(igr)
cc        nz = MGgrid_sv%nzgl(igr)
cc
cc        nc = nx*ny*nz*neq
cc
ccc     Ensure serial treatment of boundary conditions
cc
cc        !MISSING STUFF HERE
cc
ccc     Create coarse global grid structure variable to pass to filial MG call
cc
cc        call createMGgrid(grid_params)  !May be wrong interface...
cc          
ccc     Reallocate MG pointers to coarse grid level igr
cc
cc        call deallocPointers(.true.)
cc        call allocPointers(neq,fpointers)
cc
ccc     Allocate coarse global quantities
cc
cc        allocate(b_global(nc),x_global(nc))
cc
cc        x_global = 0d0
cc
ccc     Gather rhs b in all processors
cc
cc#if defined(petsc)
cc        call MPI_Allgather(yy(isig),nn,MPI_DOUBLE_PRECISION
cc     .                    ,b_global,nn,MPI_DOUBLE_PRECISION
cc     .                    ,MPI_COMM_WORLD,mpierr)
ccc diag ****
cccc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cccc     .             ,'local rhs =',yy(isig:isig+nn-1)
cc
cccc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cccc     .             ,'global rhs=',b_global
cccc        call pstop('coarseGridSolve','')
ccc diag ****
cc#endif
cc
ccccc     Allocate new diagonal and transfer elements
cccc
cccc        allocate(diag(neq,2*nc))
cccc
ccccc     Gather diag in all processors and all grids
cccc
cccc#if defined(petsc)
cccc        do ig_sv=igr,gparams_sv%ngrid
cccc          ig_cs   = ig_sv - igr + 1
cccc          isig_sv = istart_sv(ig_sv)
cccc          isig_cs = istart(ig_cs)
cccc          nsv     = istart_sv(ig_sv+1)-istart_sv(ig_sv)
cccc          ncs     = istart   (ig_cs+1)-istart   (ig_cs)
cccc
ccccc diag ****
cccc          write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cccc     .               ,'global grid=',ig_sv,'local grid=',ig_cs
cccc     .               ,'size vector transferred=',nsv
cccc     .               ,'size vector received=',ncs
cccc     .               ,'global pointers=',istart_sv
cccc     .               ,'local pointers=',istart
ccccc diag ****
cccc          call MPI_Allgather(diag_sv(:,isig_sv:isig_sv+nsv-1)
cccc     .                              ,neq*nsv,MPI_DOUBLE_PRECISION
cccc     .                      ,diag   (:,isig_cs:isig_cs+ncs-1)
cccc     .                              ,neq*nsv,MPI_DOUBLE_PRECISION
cccc     .                      ,MPI_COMM_WORLD,mpierr)
cccc        enddo
cccc
ccccc diag ****
cccc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cccc     .             ,'local diag =',diag_sv(1,isig:isig+nn-1)
cccc
cccccc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cccccc     .             ,'global diag=',diag(1,1:nc)
cccc        call pstop('coarseGridSolve','')
ccccc diag ****
cccc#endif
cc
ccc     Configure recursive plane/line MG solve
cc
cc        options%tol                         = mgtol
cc        options%vcyc                        = vcyc
cc        options%mg_gmin                     = igridmin
cc        options%orderres                    = orderres
cc        options%orderprol                   = orderprol
cc        options%mg_mu                       = mu
cccc        options%vol_res                     = volf        
cc        options%mg_line_relax               = .false.
cc        options%mg_coarse_solver_depth      = crsedpth 
cc        options%mg_vertex_relax             = vbr_mg
cc        options%mg_galerkin                 = galerkin
cc        options%ncolors                     = ncolors
cccc        options%mg_grid_def                 = grid_params
cc        call equateGridStructure(options%mg_grid_def,grid_params)
cc
cc        options%fdiag                       = .true.
cc
ccc     Call plane/line MG
cc
ccc diag ****
cc        if (my_rank == 0) call writeGridStructure(grid_params)
cc        write (*,*) size(b_global),size(x_global),nc,neq
cc        call pstop('coarseGridSolve','')
ccc diag ****
cc
cc        guess = 0
cc        call mg(neq,nc,b_global,x_global,matvec,options,1,bcnd
cc     .         ,guess,outc-1,depth)
cc
ccc     Scatter x_global
cc
cc        do k = MGgrid_sv%klo(igr),MGgrid_sv%khi(igr)
cc          do j = MGgrid_sv%jlo(igr),MGgrid_sv%jhi(igr)
cc            do i = MGgrid_sv%ilo(igr),MGgrid_sv%ihi(igr)
cc              call fromGlobalToLocalLimits(i,j,k,il,jl,kl,igr,igr,igr)
cc              ii  = vecPos(neq,il,jl,kl,igr,igr,igr) + isig - 1
cc              iig = neq*(i-1 + nx*(j-1) + nx*ny*(k-1))
cc
cc              xx(ii+1:ii+neq) = x_global(iig+1:iig+neq)
cc            enddo
cc          enddo
cc        enddo
cc
ccc diag ****
cc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cc     .             ,xx(isig:isig+nn)
cc        call pstop('coarseGridSolve','')
ccc diag ****
cc
ccc     Recover configuration for parent level MG
cc
cc        !Recover MG grid def
cc        call restoreMGgrid
cc
cc        !Recover grid_params
cc        call deallocateGridStructure(grid_params)
cccc        grid_params = gparams_sv
cc        call equateGridStructure(grid_params,gparams_sv)
cc        call deallocateGridStructure(gparams_sv)
cc
cc        !Recover pointers
cc        call deallocPointers(fpointers)
cc        call allocPointers(neq,fpointers)
cc
cc        !Recover diagonal
cc        if (.not.options%fdiag) deallocate(diag)
cc        allocate(diag(size(diag_sv,1),size(diag_sv,2)))
cc        diag = diag_sv
cc        deallocate(diag_sv)
cc
ccc     Deallocate memory
cc
cc        call deallocateGridStructure(options%mg_grid_def)
cc
cc      end subroutine coarseGridSolve

c     killmg
c     ###################################################################
      recursive subroutine killmg

        implicit none

        if (fdiag) then
          deallocate(diag)
        else
          nullify(diag)
        endif

        call deallocPointers(fpointers,mgctx=mgctx)

      end subroutine killmg

      end subroutine mg

c jb
c#######################################################################
      recursive subroutine jb(neq,ntot,rr,zz,matvec,options,igrid,bcnd
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free Jacobi routine to solve Azz = rr. Call variables:
c       * neq: number of equations (system JB)
c       * ntot: vectors dimension (neq*grid dimension)
c       * rr,zz: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0 -> no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mg_solver

      implicit none       !For safe fortran

c Call variables

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer    :: iter,alloc_stat,isig,itr,nn,ieq,irbg1,nrbg1,igridc
      integer    :: imin,imax,jmin,jmax,kmin,kmax,iming,itr2,ng
      real(8)    :: omega0,omega10,omega01,tol
      logical    :: fdiag,fpointers,vbr,prnt,smat

      integer    :: i,j,k,iv,jv,kv,if,jf,kf,nxf,nyf,nzf,iiv,ivg,ncolors
      integer    :: ii,iii,iib,iiib,iic,iig,iblock,jblock,kblock,nblk,i1
      real(8)    :: mag0=1d0,mag1,mag,magsp,magsp0
     $             ,lmag(1+neq),yy(ntot),ones(ntot)

      real(8),allocatable, dimension(:)  :: dummy,rhs

      real(8), pointer, dimension(:,:) :: diag

      type(mg_ctx),pointer :: mgctx

c Begin program

      prnt = (my_rank == 0) .or. asm

      if (out.ge.2.and.prnt) write (*,*)

c Read solver configuration

      call optionsConsistencyCheck

      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag
      vbr    = options%mg_vertex_relax
      ncolors= options%ncolors
      diag   => options%diag
      mgctx  => options%mgctx
      smat   = options%singular_matrix

c Allocate pointers

      call allocPointers(neq,fpointers,mgctx=mgctx)

      if (fpointers) then          !JB is NOT called from MG
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2.and.prnt)
     .     write (*,*) 'Allocating pointers...'

c Initialize auxiliary variables

      !Grid dimensions
      nxf = mgctx%g_def%nxv(igrid)
      nyf = mgctx%g_def%nyv(igrid)
      nzf = mgctx%g_def%nzv(igrid)

      if (vbr.or.vbr_mg) then
        nblk = nblock(igrid)
        isig = istartb(igrid)
        vbr  = .true.
      else
        nblk = 1
        isig = istart(igrid)
      endif

c Find diagonal for smoothers

      if (associated(diag)) then !Diagonal provided externally
        if(out.ge.2.and.prnt) write(*,*)'Diagonal externally provided'
        fdiag = .false.
      else                      !Form diagonal
        if (out.ge.2.and.prnt) write (*,*) 'Forming diagonal...'
        allocate(diag(neq*nblk,2*ntot*nblk))
        call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
     .                   ,ncolors,res_vol_wgt,galerkin,mgctx=mgctx)
        if (out.ge.2.and.prnt) write (*,*) 'Finished!'
        fdiag = .true.
      endif

c Preparation for iteration

      nn = ntot

      ones = 1d0

      if (guess.eq.0) then
        zz = 0d0
        mag0 = sqrt(MGdot(mgctx,igrid,rr,rr)) !sqrt(yy*yy)
      else
        call MGmatvec(mgctx,0,zz,yy,igrid,bcnd,matvec)
        call MGvecadd(mgctx,igrid,-1d0,yy,1d0,rr) !yy=rr-yy
        mag0 = sqrt(MGdot(mgctx,igrid,yy,yy)) !sqrt(yy*yy)
       endif

      mag = mag0

      if (out.ge.2.and.prnt) write(*,10) 0,mag,1d0

      if (mag0.lt.1d-20*nn) then
        mag0 = 1d0
        mag1 = 1d0
        iter = 0
      endif

c Jacobi iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag1 = mag

        call MGmatvec(mgctx,0,zz,yy,igrid,bcnd,matvec)

c       VERTEX-BASED RELAXATION
        if (vbr) then

          call vertex_JB

c       STANDARD RELAXATION
        else

          call std_JB

        endif

c     Eliminate null space

cc        if (smat) zz = zz - dot2(zz,ones)/dot2(ones,ones)*ones

c     Check convergence

        if (out.ge.2.and.prnt) write(*,20) itr,mag,mag/mag0,mag/mag1

        if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.2.and.prnt) then
        write (*,*)
      elseif (out.eq.1.and.prnt) then
        write (*,10) itr,mag,mag/mag0
      endif

      options%iter_out = itr
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag) then
        deallocate(diag)
      else
        nullify(diag)
      endif

      call deallocPointers(fpointers,mgctx=mgctx)

 10   format (' JB Iteration:',i4,'; Residual:',1p,1e10.2,
     .        '; Ratio:',1e10.2)
 20   format (' JB Iteration:',i4,'; Residual:',1p,1e10.2,
     .        '; Ratio:',1e10.2,'; Damping:',1e10.2)

      contains

c     optionsConsistencyCheck
c     ###################################################################
      subroutine optionsConsistencyCheck

        if (options%mg_vertex_relax .and. options%mg_line_relax) then
          if (prnt)
     .       write (*,*) 'Invalid setting for JB relaxation: cannot do'
          if (prnt)
     .       write (*,*) 'vertex-based and line-based relaxation',
     .                  ' simultaneously!'
          call pstop('optionsConsistencyCheck in jb','Aborting...')
        endif
          
      end subroutine optionsConsistencyCheck

c     std_JB
c     ###################################################################
      subroutine std_JB

        mag  = 0d0

        call limits(0,nxf,nyf,nzf,igrid,imin,imax,jmin,jmax,kmin,kmax)

        do k=kmin,kmax
          do j=jmin,jmax
            do i=imin,imax

              ii = i + nxf*(j-1) + nxf*nyf*(k-1)

              iii = neq*(ii-1)

              !Find new residual
              yy(iii+1:iii+neq) = rr(iii+1:iii+neq) - yy(iii+1:iii+neq)

              !Multiply by D^-1 (stored in diag)
              iig = iii + isig - 1
              dummy = matmul(transpose(diag(:,iig+1:iig+neq))
     $                                ,yy(iii+1:iii+neq))

              !Update zz
              zz(iii+1:iii+neq) = zz(iii+1:iii+neq) + omega0*dummy

            enddo
          enddo
        enddo

        mag = sqrt(MGdot(mgctx,igrid,yy,yy))       !sqrt(yy*yy)

cc          !Take care of points around singular point
cc          if (isSP(1,1,1,igrid,igrid,igrid)) then
cc#if defined(petsc)
cc            if (np > 1) call pstop('jb','Does not work in parallel')
cc#endif
cc            do i1=3,1,-1
cc
cc              !plane JB
cc              mg_ctx_temp%iline(igrid) = i1 !Set i for line smoothing around SP
cc
cc              call limits(0,nxf,nyf,nzf,igrid
cc     .             ,imin,imax,jmin,jmax,kmin,kmax)
cc
cc              do itr2=1,5 !JB iteration around SP
cc
cc                magsp = 0d0
cc
cc                call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)
cc
cc                do k=kmin,kmax
cc                  do j=jmin,jmax
cc                    do i=imin,imax
cc                      !nxf,nyf,nzf are the original grid's even in line relaxation.
cc                      ii = i + nxf*(j-1) + nxf*nyf*(k-1)
cc
cc                      iii = neq*(ii-1)
cc
cc                      !Find new residual
cc                      do ieq = 1,neq
cc                        rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
cc                      enddo
cc
cc                      !Multiply by D^-1 (stored in diag)
cc                      iig = iii + isig - 1
cc                      dummy =
cc     .                   matmul(transpose(diag(:,iig+1:iig+neq))
cc     .                         ,rhs)
cc
cc                      !Update solution zz
cc                      do ieq=1,neq
cc                        zz(iii+ieq) = zz(iii+ieq) + 0.8*dummy(ieq)
cc                      enddo
cc
cc                      magsp = magsp + sum(rhs*rhs)
cc                    enddo
cc                  enddo
cc                enddo
cc
cc                if (itr2 == 1) then
cc                  magsp0 = magsp
cc                else
cc                  if (out.ge.3.and.prnt)
cc     .              write (*,*) '  SP smoothing in JB, iter=',itr2
cc     .                             ,' Ratio=',sqrt(magsp/magsp0)
cc                  if (sqrt(magsp/magsp0) < 0.5) exit
cc                endif
cc
cc              enddo
cc
cc              if (out.eq.2.and.prnt)
cc     .           write (*,*) '  SP smoothing in GS: line = ',i1
cc     .                       ,'iter=',itr2,'Ratio=',sqrt(magsp/magsp0)
cc
cc            enddo
cc
cc            mg_ctx_temp%iline(igrid) = 0 !Unset line smoothing
cc
cc          endif

      end subroutine std_JB

c     vertex_JB
c     ###################################################################
      subroutine vertex_JB

        mag  = 0d0

        !Vertex sampling
        do kv=1,max(nzf-mg_ratio_z(igrid)+1,1)
          do jv=1,max(nyf-mg_ratio_y(igrid)+1,1)
            do iv=1,max(nxf-mg_ratio_x(igrid)+1,1)

            !Find new residual
            do kblock=1,mg_ratio_z(igrid)
              do jblock=1,mg_ratio_y(igrid)
                do iblock=1,mg_ratio_x(igrid)

                  iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                         + mg_ratio_x(igrid)
     .                          *mg_ratio_y(igrid)*(kblock-1)

                  if = iv-1 + iblock
                  jf = jv-1 + jblock
                  kf = kv-1 + kblock

                  ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                  do ieq = 1,neq
                    iii  = neq*(ii -1) + ieq
                    iiib = neq*(iib-1) + ieq
                    rhs(iiib) = rr(iii) - yy(iii)
                  enddo

                enddo
              enddo
            enddo

           !Multiply by D^-1 (stored in diag)
            iiv = iv + nxf*(jv-1) + nxf*nyf*(kv-1)
            ivg  = (iiv-1)*neq*nblk + isig -1
            dummy = matmul(transpose(diag(:,ivg+1:ivg+neq*nblk)),rhs)

           !Update zz
            do kblock=1,mg_ratio_z(igrid)
              do jblock=1,mg_ratio_y(igrid)
                do iblock=1,mg_ratio_x(igrid)

                  iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                         + mg_ratio_x(igrid)
     .                          *mg_ratio_y(igrid)*(kblock-1)

                  if = iv-1 + iblock
                  jf = jv-1 + jblock
                  kf = kv-1 + kblock

                  ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                  do ieq = 1,neq
                    iii  = neq*(ii -1) + ieq
                    iiib = neq*(iib-1) + ieq
                    zz(iii) = zz(iii) + omega0*dummy(iiib)
                  enddo

                enddo
              enddo
            enddo

            mag = mag + sum(rhs*rhs)

            enddo
          enddo
        enddo

#if defined(petsc)
        if (.not.asm) then
          lmag = mag
          call MPI_Allreduce(lmag,mag,1,MPI_DOUBLE_PRECISION
     .                       ,MPI_SUM,mgctx%mpi_comm,mpierr)
        endif
#endif

        mag = sqrt(mag)

      end subroutine vertex_JB

      end subroutine jb

c gs
c#######################################################################
      recursive subroutine gs(neq,ntot,rr,zz,matvec,options,igrid,bcnd
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free Gauss-Seidel routine to solve A*zz=rr. Call variables:
c       * neq: number of equations (system GS)
c       * ntot: vectors dimension (neq*grid dimension)
c       * rr,zz: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0 -> no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mg_solver

      implicit none       !For safe fortran

c Call variables

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer    :: iter,alloc_stat,isig,itr,itr2,nn,ieq,iig
     .             ,zebra_dir,minzd,maxzd,zebra_it,zd,zit,ng
      integer    :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
      integer    :: i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max

      real(8)    :: omega0,tol,zebra_omg
      logical    :: fdiag,fpointers,vbr,zebra,prnt,smat

      integer    :: ic,jc,kc,iv,jv,kv,if,jf,kf,nxf,nyf,nzf
      integer    :: ii,iii,iib,iiib,iiv,ivg,iblock,jblock,kblock,nblk
      real(8)    :: mag0,mag1,magsp,magsp0,mag,lmag,magz,magz0
     $             ,yy(ntot),ones(ntot)

      integer    :: ncolors,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3

      real(8), pointer, dimension(:,:) :: diag

      real(8),allocatable, dimension(:) :: dummy,rhs

      type(mg_ctx),pointer :: mgctx => null()

      character(7) :: gs_type

c Begin program

      prnt = (my_rank == 0) .or. asm

      if (out.ge.2.and.prnt) write (*,*)

c Read solver configuration

      call optionsConsistencyCheck

      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag
      vbr    = options%mg_vertex_relax

      zebra     = options%mg_zebra_relax
      zebra_dir = options%mg_zebra_prefd
      zebra_it  = options%mg_zebra_it
      zebra_omg = options%mg_zebra_omega

      ncolors= options%ncolors
      diag   => options%diag
      mgctx  => options%mgctx

c Set pointers

      call allocPointers(neq,fpointers,mgctx=mgctx)

      if (fpointers) then          !GS is NOT called from MG
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2.and.prnt)
     .     write (*,*) 'Allocated pointers.'

c Initialize auxiliary variables

      !Grid size
      nxf = mgctx%g_def%nxv(igrid)
      nyf = mgctx%g_def%nyv(igrid)
      nzf = mgctx%g_def%nzv(igrid)

      if (vbr.or.vbr_mg) then
        nblk = nblock(igrid)
        isig = istartb(igrid)
        vbr  = .true.
      else
        nblk = 1
        isig = istart(igrid)
      endif

c Find diagonal for smoothers

      if (associated(diag)) then !Diagonal provided externally
        if(out.ge.2.and.prnt) write(*,*)'Diagonal externally provided'
        fdiag = .false.
      else                      !Form diagonal
        if (out.ge.2.and.prnt) write (*,*) 'Forming diagonal...'
        allocate(diag(neq*nblk,2*ntot*nblk))
        call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
     .                   ,ncolors,res_vol_wgt,galerkin,mgctx=mgctx)
        if (out.ge.2.and.prnt) write (*,*) 'Finished!'
        fdiag = .true.
      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) then
        zz = 0d0
        mag0 = sqrt(MGdot(mgctx,igrid,rr,rr)) !sqrt(yy*yy)
      else
        call MGmatvec(mgctx,0,zz,yy,igrid,bcnd,matvec)
        call MGvecadd(mgctx,igrid,-1d0,yy,1d0,rr) !yy=rr-yy
        mag0 = sqrt(MGdot(mgctx,igrid,yy,yy)) !sqrt(yy*yy)
       endif

      mag = mag0

      if (out.ge.2.and.prnt) write(*,10) '',0,mag,1d0

      if (mag0.lt.1d-20*nn) then
        mag0 = 1d0
        mag1 = 1d0
        iter = 0
      endif

      ones = 1d0

c GS iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag1 = mag

c       ----------
c       Colored GS
c       ---------- 
        if (ncolors > 1) then

c         VERTEX-BASED RELAXATION
          if (vbr) then

            gs_type = 'C-VERTX'

            call colored_vrtx_GS

c         ZEBRA LINE RELAXATION
          elseif (zebra) then

            gs_type = 'C-ZEBRA'

            call zebra_GS

c         STANDARD RELAXATION
          else

            gs_type = 'Colored'

            call colored_GS

          endif

c       -----------
c       Standard GS
c       ----------- 
        else

c         VERTEX-BASED RELAXATION
          if (vbr) then
            gs_type = 'L-VERTX'

            call vertex_GS

c         STANDARD RELAXATION
          else
 
           gs_type = 'Lexicog'

           call lexic_GS

          endif

        endif

c     Eliminate null space

cc        if (smat) zz = zz - dot2(zz,ones)/dot2(ones,ones)*ones

c     Check convergence

        if (out.ge.2.and.prnt)
     $       write(*,20) gs_type,itr,mag,mag/mag0,mag/mag1

        if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.2.and.prnt) then
        write (*,*)
      elseif (out.eq.1.and.prnt) then
        write (*,10) gs_type,itr,mag,mag/mag0
      endif

      options%iter_out = itr
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag) then
        deallocate(diag)
      else
        nullify(diag)
      endif

      call deallocPointers(fpointers,mgctx=mgctx)

      return

 10   format (a,' GS Iteration:',i4,'; Residual:',1p,1e10.2,
     .        '; Ratio:',1e10.2)
 20   format (a,' GS Iteration:',i4,'; Residual:',1p,1e10.2,
     .        '; Ratio:',1e10.2,'; Damping:',1e10.2)

      contains

c     optionsConsistencyCheck
c     ###################################################################
      subroutine optionsConsistencyCheck

        if ( (options%mg_vertex_relax.and.options%mg_zebra_relax))then
          if (prnt) write (*,*) 'Invalid relaxation choice for GS'
          call pstop('optionsConsistencyCheck in gs','Aborting...')
        elseif (options%mg_vertex_relax .or. vbr_mg) then
          options%ncolors = 4
        elseif (options%mg_line_relax) then
          options%ncolors = 1
cc        elseif (options%mg_zebra_relax) then
cc          options%omega = 0.7   !Do damped JB along lines
cccc          options%omega = 1.
        endif
          
      end subroutine optionsConsistencyCheck

c     zebra_GS
c     ###################################################################
      subroutine zebra_GS

        mag = 0d0

        if (zebra_dir == 0) then
          minzd = 1
          maxzd = 3
        else
          minzd = zebra_dir
          maxzd = zebra_dir
        endif

        do zd=minzd,maxzd  !Cycle in zebra directions

          call scheduleColors(mgctx,igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                     ,i1min,i1max,j1min,j1max,k1min,k1max
     .                     ,zebra_dir=zd)

          if (out.ge.2.and.prnt) then
            if (nrbg1/=0 .and. nrbg2/=0 .and. nrbg3/=0) then
              write (*,'(a,i1)') 'Zebra GS, line = ',zd
            else
              write (*,'(a,i1)') 'Skipping GS zebra line = ',zd
            endif
          endif

          if (out.ge.3.and.prnt) then
            if (nrbg1/=0 .and. nrbg2/=0 .and. nrbg3/=0) then
              write (*,*) '  colors =',nrbg1,nrbg2,nrbg3
              write (*,*) '  ilimits=',i1min,i1max
              write (*,*) '  jlimits=',j1min,j1max
              write (*,*) '  klimits=',k1min,k1max
            endif
          endif

          if (     (minzd==maxzd)
     .        .and.(nrbg1==0.and.nrbg2==0.and.nrbg3==0)) then
            if (prnt) then
              write (*,*) 'zebra dir=',zd,' range=',minzd,maxzd
              write (*,*) 'colors =',nrbg1,nrbg2,nrbg3
              write (*,*) 'ilimits=',i1min,i1max
              write (*,*) 'jlimits=',j1min,j1max
              write (*,*) 'klimits=',k1min,k1max
            endif
            call pstop('gs','No zebra GS relaxation is possible!')
          endif

          !Cycle in colors
          do irbg3 = 1,nrbg3
           do irbg2 = 1,nrbg2
            do irbg1 = 1,nrbg1

              !Inner ZEBRA iterations
              do zit=1,zebra_it 

                magz = 0d0

                call MGmatvec(mgctx,0,zz,yy,igrid,bcnd,matvec)

                do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                 do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                  do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                    call scheduleIndices(mgctx,igrid,i1,j1,k1,i,j,k
     .                                  ,zebra_dir=zd)

                    if (out.ge.4.and.prnt)
cc                    if (prnt)
     .                 write (*,'(a,i4,a,i2,a,3i4)')
     .                         'DIAG -- zebra GS: it= ',zit
     .                        ,', dir=',zd
     .                        ,', grid pos=',i,j,k

c diag ****
cc                    ii = i + nxf*(j-1) + nxf*nyf*(k-1)
cc                    call MGmatvec(-ii,neq,ntot,zz,yy,igrid,bcnd,matvec)
c diag ****

                    !nxf,nyf,nzf are the original grid's even in line relaxation.
                    iii = neq*(i-1 + nxf*(j-1) + nxf*nyf*(k-1))

                    !Find new residual
                    rhs = rr(iii+1:iii+neq) - yy(iii+1:iii+neq)

                    !Multiply by D^-1 (stored in diag)
                    iig = iii + isig - 1
                    dummy=matmul(transpose(diag(:,iig+1:iig+neq)),rhs)

                    !Update zz
                    zz(iii+1:iii+neq)=zz(iii+1:iii+neq)+zebra_omg*dummy

                    magz = magz + sum(rhs*rhs)  !Zebra color accumulator (sum over same color)

cc                    summ = summ + zz(iii+1:iii+neq)

                  enddo
                 enddo
                enddo

                if (zit == 1) magz0 = magz

                if (out.ge.3.and.prnt) then
                  write (*,'(a,i2,a,i4,a,i2,a,1pe10.2)')
     .                            '   --> Inner zebra GS: dir=',zd
     .                           ,'; it= ',zit
     .                           ,'; color=',irbg1
     .                              + nrbg1*(irbg2-1)
     .                              + nrbg2*(irbg3-1)
     .                           ,'; residual=',sqrt(magz/magz0)
                endif

              enddo

              mag = mag + magz !Update outer accumulator (sum over colors)

              if (out == 2.and.prnt) then
                write (*,'(a,i2,a,i4,a,i2,a,1pe10.4)')
     .                            '   --> Inner zebra GS: dir=',zd
     .                           ,'; it= ',zebra_it
     .                           ,'; color=',irbg1
     .                              + nrbg1*(irbg2-1)
     .                              + nrbg2*(irbg3-1)
     .                           ,'; residual=',sqrt(magz/magz0)
              endif

             enddo
           enddo
          enddo

        enddo

#if defined(petsc)
        if (.not.asm) then
          lmag = mag
          call MPI_Allreduce(lmag,mag,1,MPI_DOUBLE_PRECISION
     .                       ,MPI_SUM,mgctx%mpi_comm,mpierr)
        endif
#endif

        mag = sqrt(mag)

      end subroutine zebra_GS

c     colored_vrtx_GS
c     ###################################################################
      subroutine colored_vrtx_GS

        mag = 0d0

        call scheduleColors(mgctx,igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                     ,i1min,i1max,j1min,j1max,k1min,k1max)

        do irbg3 = 1,nrbg3
          do irbg2 = 1,nrbg2
            do irbg1 = 1,nrbg1

            call MGmatvec(mgctx,0,zz,yy,igrid,bcnd,matvec)

            !Vertex sampling
            do kv=1+mod((    irbg3-1),nrbg3)
     .                   ,max(nzf-mg_ratio_z(igrid)+1,1),nrbg3
              do jv=1+mod((kv   +irbg2+irbg1-1),nrbg2)
     .                     ,max(nyf-mg_ratio_y(igrid)+1,1),nrbg2
                do iv=1+mod((jv+kv+irbg1-1),nrbg1)
     .                       ,max(nxf-mg_ratio_x(igrid)+1,1),nrbg1

                 !Find new residual
                  do kblock=1,mg_ratio_z(igrid)
                    do jblock=1,mg_ratio_y(igrid)
                      do iblock=1,mg_ratio_x(igrid)

                        iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                               + mg_ratio_x(igrid)
     .                                *mg_ratio_y(igrid)*(kblock-1)

                        if = iv-1 + iblock
                        jf = jv-1 + jblock
                        kf = kv-1 + kblock

                        ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                        do ieq = 1,neq
                          iii  = neq*(ii -1) + ieq
                          iiib = neq*(iib-1) + ieq
                          rhs(iiib) = rr(iii) - yy(iii)
                        enddo

                      enddo
                    enddo
                  enddo

                  !Multiply by D^-1 (stored in diag)
                  iiv = iv + nxf*(jv-1) + nxf*nyf*(kv-1)
                  ivg  = (iiv-1)*neq*nblk + isig -1
                  dummy=matmul(transpose(diag(:,ivg+1:ivg+neq*nblk))
     .                        ,rhs)

                  !Update zz
                  do kblock=1,mg_ratio_z(igrid)
                    do jblock=1,mg_ratio_y(igrid)
                      do iblock=1,mg_ratio_x(igrid)

                        iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                               + mg_ratio_x(igrid)
     .                                *mg_ratio_y(igrid)*(kblock-1)

                        if = iv-1 + iblock
                        jf = jv-1 + jblock
                        kf = kv-1 + kblock

                        ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                        do ieq = 1,neq
                          iii  = neq*(ii -1) + ieq
                          iiib = neq*(iib-1) + ieq
                          zz(iii) = zz(iii) + omega0*dummy(iiib)
                        enddo

                      enddo
                    enddo
                  enddo

                  mag = mag + sum(rhs*rhs)

                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo

#if defined(petsc)
        if (.not.asm) then
          lmag = mag
          call MPI_Allreduce(lmag,mag,1,MPI_DOUBLE_PRECISION
     .                       ,MPI_SUM,mgctx%mpi_comm,mpierr)
        endif
#endif

        mag = sqrt(mag)

      end subroutine colored_vrtx_GS

c     colored_GS
c     ###################################################################
      subroutine colored_GS

        mag = 0d0

        call scheduleColors(mgctx,igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                     ,i1min,i1max,j1min,j1max,k1min,k1max)

        if (out.ge.3.and.prnt) then
          write (*,*) '  colors=',nrbg1,nrbg2,nrbg3
          write (*,*) '  ilimits=',i1min,i1max
          write (*,*) '  jlimits=',j1min,j1max
          write (*,*) '  klimits=',k1min,k1max
        endif

        do irbg3 = 1,nrbg3
          do irbg2 = 1,nrbg2
            do irbg1 = 1,nrbg1

            call MGmatvec(mgctx,0,zz,yy,igrid,bcnd,matvec)

            do k1=k1min+mod((    irbg3-1),nrbg3),k1max,nrbg3
              do j1=j1min+mod((k1+irbg2+irbg1-1),nrbg2),j1max,nrbg2
                do i1=i1min+mod((j1+k1+irbg1-1),nrbg1),i1max,nrbg1

                  call scheduleIndices(mgctx,igrid,i1,j1,k1,i,j,k)

                  !nxf,nyf,nzf are the original grid's even in line relaxation.
                  iii = neq*(i-1 + nxf*(j-1) + nxf*nyf*(k-1))

                  !Find new residual
                  yy(iii+1:iii+neq) = rr(iii+1:iii+neq)
     $                              - yy(iii+1:iii+neq)

                  !Multiply by D^-1 (stored in diag)
                  iig = iii + isig - 1
                  dummy=matmul(transpose(diag(:,iig+1:iig+neq))
     .                        ,yy(iii+1:iii+neq))

                  !Update zz (except around SP)
                  zz(iii+1:iii+neq) = zz(iii+1:iii+neq) + omega0*dummy

                enddo
              enddo

            enddo

            enddo
          enddo
        enddo

        mag = sqrt(MGdot(mgctx,igrid,yy,yy))

      end subroutine colored_GS

c     vertex_GS
c     ###################################################################
      subroutine vertex_GS

        mag = 0d0

        !Vertex sampling
        do kv=1,max(nzf-mg_ratio_z(igrid)+1,1)
          do jv=1,max(nyf-mg_ratio_y(igrid)+1,1)
            do iv=1,max(nxf-mg_ratio_x(igrid)+1,1)

            !Find new residual
            do kblock=1,mg_ratio_z(igrid)
              do jblock=1,mg_ratio_y(igrid)
                do iblock=1,mg_ratio_x(igrid)

                iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                       + mg_ratio_x(igrid)
     .                        *mg_ratio_y(igrid)*(kblock-1)

                if = iv-1 + iblock
                jf = jv-1 + jblock
                kf = kv-1 + kblock

                ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                call MGmatvec(mgctx,-ii,zz,yy,igrid,bcnd,matvec)

                do ieq = 1,neq
                  iii  = neq*(ii -1) + ieq
                  iiib = neq*(iib-1) + ieq
                  rhs(iiib) = rr(iii) - yy(iii)
                enddo

                enddo
              enddo
            enddo

            !Multiply by D^-1 (stored in diag)
            iiv = iv + nxf*(jv-1) + nxf*nyf*(kv-1)
            ivg  = (iiv-1)*neq*nblk + isig -1
            dummy=matmul(transpose(diag(:,ivg+1:ivg+neq*nblk)),rhs)

            !Update zz
            do kblock=1,mg_ratio_z(igrid)
              do jblock=1,mg_ratio_y(igrid)
                do iblock=1,mg_ratio_x(igrid)

                iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                       + mg_ratio_x(igrid)
     .                        *mg_ratio_y(igrid)*(kblock-1)

                if = iv-1 + iblock
                jf = jv-1 + jblock
                kf = kv-1 + kblock

                ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                do ieq = 1,neq
                  iii  = neq*(ii -1) + ieq
                  iiib = neq*(iib-1) + ieq
                  zz(iii) = zz(iii) + omega0*dummy(iiib)
                enddo

                enddo
              enddo
            enddo

            mag = mag + sum(rhs*rhs)

            enddo
          enddo
        enddo

#if defined(petsc)
        if (.not.asm) then
          lmag = mag
          call MPI_Allreduce(lmag,mag,1,MPI_DOUBLE_PRECISION
     .                       ,MPI_SUM,mgctx%mpi_comm,mpierr)
        endif
#endif

        mag = sqrt(mag)

      end subroutine vertex_GS

c     lexic_GS
c     ###################################################################
      subroutine lexic_GS

        mag = 0d0

        call limits(0,nxf,nyf,nzf,igrid
     .             ,imin,imax,jmin,jmax,kmin,kmax)

c       Forward pass

        do k=kmin,kmax
          do j=jmin,jmax
            do i=imin,imax

              !nxf,nyf,nzf are the original grid's even in line relaxation.
              ii = i + nxf*(j-1) + nxf*nyf*(k-1)

              call MGmatvec(mgctx,-ii,zz,yy,igrid,bcnd,matvec)

              iii = neq*(ii-1)

              !Find new residual
              yy(iii+1:iii+neq) = rr(iii+1:iii+neq) - yy(iii+1:iii+neq)

              !Multiply by D^-1 (stored in diag)
              iig = iii + isig - 1
              dummy = matmul(transpose(diag(:,iig+1:iig+neq))
     $                                ,yy(iii+1:iii+neq))

              !Update solution zz
              zz(iii+1:iii+neq) = zz(iii+1:iii+neq) + omega0*dummy

            enddo
          enddo
        enddo
        
c       Backward pass

        do k=kmax,kmin,-1
          do j=jmax,jmin,-1
            do i=imax,imin,-1

              ii = i + nxf*(j-1) + nxf*nyf*(k-1)

              call MGmatvec(mgctx,-ii,zz,yy,igrid,bcnd,matvec)

              iii = neq*(ii-1)

              !Find new residual
              yy(iii+1:iii+neq) = rr(iii+1:iii+neq) - yy(iii+1:iii+neq)

              !Multiply by D^-1
              iig = iii + isig - 1
              dummy = matmul(transpose(diag(:,iig+1:iig+neq))
     $                                ,yy(iii+1:iii+neq))

              !Update solution zz
              zz(iii+1:iii+neq) = zz(iii+1:iii+neq) + omega0*dummy

            enddo
          enddo
        enddo

        mag = sqrt(MGdot(mgctx,igrid,yy,yy))

      end subroutine lexic_GS

      end subroutine gs

c     symm_test
c     ###############################################################
      subroutine symm_test(neq,igrid,matvec,bcnd)
c     ---------------------------------------------------------------
c     Performs symmetry test of matvec on grid "igrid".
c     ---------------------------------------------------------------

      use mg_solver

      implicit none     !For safe fortran

c     Call variables

      integer    :: neq,igrid,bcnd(6,neq)

      external      matvec

c     Local variables

      real(8),allocatable,dimension(:)::x1,dummy,dummy2

      real(8)    :: dd1,dd2,err
      integer    :: nx,ny,nz,nn,ii,jj,i1,j1,i2,j2,ix1,iy1,ix2,iy2,ieq

      logical :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Initialize variables

      nx = MGgrid%nxv(igrid)
      ny = MGgrid%nyv(igrid)
      nz = MGgrid%nzv(igrid)

      nn = neq*nx*ny*nz

      write (*,*) 'Performing symmetry test of system matrix ',
     .            'on grid:',nx,'x',ny,'x',nz,'...'

      allocate(x1(nn),dummy(nn),dummy2(nn))

      err = 0d0

      do ii = 1,nn
        x1    (ii) = 0d0
        dummy (ii) = 0d0
        dummy2(ii) = 0d0
      enddo

c     Check symmetry

      do ii = 1,nn/neq

        do ieq =1,neq

c       Find column vector ii

          call findBaseVector(ii,ieq,neq,nn,x1,1d0)

          call MGmatvec(mg_ctx_temp,0,x1,dummy,igrid,bcnd,matvec)

          call findBaseVector(ii,ieq,neq,nn,x1,0d0)

c       Compare column vector ii with corresponding row vector (intersect in
c       diagonal)

          do jj = ii,nn/neq

            call findBaseVector(jj,ieq,neq,nn,x1,1d0)

c THIS DOESN'T WORK ANYMORE. NEED TO FIX.
cc            call matvec(ii,neq,nn,x1,dummy2,igrid,bcnd)
            call MGmatvec(mg_ctx_temp,0,x1,dummy2,igrid,bcnd,matvec)

            call findBaseVector(jj,ieq,neq,nn,x1,0d0)

            dd1 = abs(dummy(jj) - dummy2(ii))
            if(abs(dummy(jj)).gt.1d-15.or.abs(dummy2(ii)).gt.1d-15) then
              write(*,15) jj,ii,dummy(jj),ii,jj,dummy2(ii),dd1
     .             ,100*dd1/max(abs(dummy(jj)),abs(dummy2(ii)))
              err = err + dd1
            endif

          enddo

        enddo
      enddo

      write (*,20) err

      deallocate (x1,dummy,dummy2)

      call deallocPointers(fpointers)

      stop

c     End program

 15   format ('(',i3,',',i3,'):',1pe10.2,'; (',i3,',',i3,'):',e10.2,
     .        '  Error:',e10.2,'  %error:',0pf7.2)
 20   format (/,'Total error:',1pe10.3)

      end subroutine symm_test

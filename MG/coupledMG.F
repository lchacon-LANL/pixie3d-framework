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

c module mg_internal
c######################################################################
      module mg_internal

        use grid_mg

        use bc_ops

        use grid_create, ONLY: grid_pack,g_pack,createMGgrid

        use mlsolverSetup

        integer, parameter :: BOX=1,STAR=2

        logical :: vbr_mg,galerkin=.false.,volf,prnt=.true.

      contains

c     getMGvcomp
c     #################################################################
      function getMGvcomp(i,j,k,nx,ny,nz,igr,ieq,neq) result(ijkg)

c     -----------------------------------------------------------------
c     Gives MG vector component corresponing to equation ieq on
c     coordinates (i,j,k) and grid level igr.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igr,ieq,neq,ijkg,nx,ny,nz

c     Local variables

        logical    :: fpointers

c     Begin program

        ijkg = neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
     .       + istart(igr) - 1

      end function getMGvcomp

c     mapArrayToMGVector
c     ##################################################################
      subroutine mapArrayToMGVector(neq,nx,ny,nz,array,mgvector,igr)
c     ------------------------------------------------------------------
c     Maps array into a MG vector excluding ghost nodes. In call sequence:
c       * neq (in): number of variables contained in array
c       * nx,ny,nz (in): grid dimensions
c       * array (in): contains variables in all grid nodes at grid level
c           igr.
c       * mgvector (out): array values stored in MG vector format, w/o
c           ghost cell values.
c       * igr (in): grid level where array is defined.
c     ------------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,igr,nx,ny,nz
      real(8)    :: mgvector(*),array(0:nx+1,0:ny+1,0:nz+1,neq)

c     Local variables

      integer    :: i,j,k,ii,ieq
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            do ieq=1,neq
cc              ii = neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
cc     .           + istart(igr) - 1
              ii = getMGvcomp(i,j,k,nx,ny,nz,igr,ieq,neq)
              mgvector(ii) = array(i,j,k,ieq)
            enddo
          enddo
        enddo
      enddo

      call deallocPointers(fpointers)

c     End program

      end subroutine mapArrayToMGVector

c     mapMGVectorToArray
c     ##################################################################
      subroutine mapMGVectorToArray(gpos,neq,mgvector,nx,ny,nz
     .                             ,array,igr,ismgvec)
c     ------------------------------------------------------------------
c     Maps values of MG vector mgvector corresponding to grid level igr
c     into array, w/o ghost cells. If vector is not a MG vector, set
c     igrid=1. In call sequence:
c       * gpos (in): if positive, map only locally in stencil of grid node
c           determined by gpos. Otherwise, map whole grid.
c       * neq (in): number of variables contained in mgvector
c       * mgvector (in): vector to be mapped
c       * nx,ny,nz( in): grid dimensions
c       * array (out): mapped array
c       * igr (in): grid level
c       * ismgvec (in): logical variable that determines whether mgvector
c           is indeed a MG vector or a regular vector.
c     ------------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,igr,nx,ny,nz,gpos
      real(8)    :: mgvector(*),array(0:nx+1,0:ny+1,0:nz+1,neq)
      logical    :: ismgvec

c     Local variables

      integer    :: ieq,i,j,k,ii,offset,igptr
      integer    :: imin ,imax ,jmin ,jmax ,kmin ,kmax
     .             ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

      call limits(gpos,nx,ny,nz,igr,imin,imax,jmin,jmax,kmin,kmax)

c     Define stencil

      offset = 1

      iimin = max(imin-offset,1)
      iimax = min(imax+offset,nx)
      jjmin = max(jmin-offset,1)
      jjmax = min(jmax+offset,ny)
      kkmin = max(kmin-offset,1)
      kkmax = min(kmax+offset,nz)

c     Map vector to array

      !Choose correct mapping for non-MG vectors
      if (ismgvec) then
        igptr = igr
      else
        igptr = 1
      endif

      array = 0d0

      do k = kkmin,kkmax
        do j = jjmin,jjmax
          do i = iimin,iimax
            do ieq=1,neq
cc              ii = neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
cc     .           + istart(igptr) - 1
              ii = getMGvcomp(i,j,k,nx,ny,nz,igptr,ieq,neq)
              array(i,j,k,ieq) = mgvector(ii)
            enddo
          enddo
        enddo
      enddo

      call deallocPointers(fpointers)

c     End program

      end subroutine mapMGVectorToArray

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

c     MGplot
c     #################################################################
      subroutine MGplot(neq,mgv,igr,iflag,ofile)

c     -----------------------------------------------------------------
c     Plots MG vector using xdraw. In call:
c       * neq: number of variables contained in MG vector
c       * mgv: MG vector
c       * igr: grid level to be plotted.
c       * iflag: whether to initialize output file (iflag=0) or not.
c       * ofile: output file.
c     -----------------------------------------------------------------

        use grid_debug

        implicit none

c     Call variables

        integer    :: iflag,igr,neq
        real(8)    :: mgv(*)
        character*(*) :: ofile

c     Local variables

        integer    :: nx,ny,nz,ieq,nunit
        real(8),allocatable,dimension(:,:,:,:) :: debug
        logical    :: fpointers

c     Begin program

        call allocPointers(neq,fpointers)

        nunit = 110

        nx = grid_params%nxv(igr)
        ny = grid_params%nyv(igr)
        nz = grid_params%nzv(igr)

c     Dump plots

        allocate(debug(0:nx+1,0:ny+1,0:nz+1,neq))

        if (iflag == 0) then
          open(nunit,file=trim(ofile),form='unformatted'
     .        ,status='replace')
        else
          open(nunit,file=trim(ofile),form='unformatted'
     .        ,status='old',position='append')
        endif

        call mapMGVectorToArray(0,neq,mgv,nx,ny,nz,debug,igr,.true.)
        do ieq=1,neq
cc          call contour(debug(1:nx,1:ny,1,ieq),nx,ny
cc     .                ,grid_params%gxmin,grid_params%gxmax
cc     .                ,grid_params%gymin,grid_params%gymax
cc     .                ,iflag+ieq-1,nunit)
          call C_contour(debug(1:nx,1:ny,1:nz,ieq)
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,grid_params%gzmin,grid_params%gzmax
     .                ,iflag+ieq-1,nunit)
        enddo

c     End program

        deallocate(debug)
        call deallocPointers(fpointers)
        close(nunit)

      end subroutine MGplot

c     MGvecadd
c     ###################################################################
      subroutine MGvecadd(igr,neq,ntot,coef1,vec1,coef2,vec2)

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

      integer    :: igr,neq,ntot
      real(8)    :: coef1,coef2,vec1(ntot),vec2(ntot)

c     Local variables

      integer    :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq,nnx,nny,nnz

c     Begin program

      nnx = grid_params%nxv(igr)
      nny = grid_params%nyv(igr)
      nnz = grid_params%nzv(igr)

      call limits(0,nnx,nny,nnz,igr,imin,imax,jmin,jmax,kmin,kmax)

cc      write (*,*) 'DIAG -- MGvecadd',nnx,nny,nnz
cc     .      ,imin,imax,jmin,jmax,kmin,kmax

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            do ieq=1,neq
              iii=getMGvcomp(i,j,k,nnx,nny,nnz,1,ieq,neq)
              vec1(iii) = coef1*vec1(iii) + coef2*vec2(iii)
            enddo
          enddo
        enddo
      enddo

      end subroutine MGvecadd

c     MGdot
c     ###################################################################
      function MGdot(igr,neq,ntot,vec1,vec2) result(dot)

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

      integer    :: igr,neq,ntot
      real(8)    :: vec1(ntot),vec2(ntot),dot

c     Local variables

      integer    :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq,nnx,nny,nnz
      real(8)    :: ldot

c     Begin program

      nnx = grid_params%nxv(igr)
      nny = grid_params%nyv(igr)
      nnz = grid_params%nzv(igr)

      call limits(0,nnx,nny,nnz,igr,imin,imax,jmin,jmax,kmin,kmax)

cc      write (*,*) 'DIAG -- MGdot',nnx,nny,nnz
cc     .         ,imin,imax,jmin,jmax,kmin,kmax

      dot = 0d0

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            do ieq=1,neq
              iii=getMGvcomp(i,j,k,nnx,nny,nnz,1,ieq,neq)
              dot = dot + vec1(iii)*vec2(iii)
            enddo
          enddo
        enddo
      enddo

#if defined(petsc)
      if (.not.asm) then
        ldot = dot
        call MPI_Allreduce(ldot,dot,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      endif
#endif

      end function MGdot

c     scheduleColors
c     ###################################################################
      subroutine scheduleColors(igrid,ncolors,nc1,nc2,nc3
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

      integer    :: igrid,ncolors,nc1,nc2,nc3
     .             ,i1min,i1max,j1min,j1max,k1min,k1max

      integer,optional :: stencil_width,stencil_type,zebra_dir

c     Local variables
      
      integer    :: nxf,nyf,nzf,nxl,nyl,nzl,nclrs
     .             ,imin,imax,jmin,jmax,kmin,kmax
     .             ,stncl_wdt,stncl_tpe,dd,ddd,zb_dir

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

      !Global grid size
      nxf = nxv(igrid)
      nyf = nyv(igrid)
      nzf = nzv(igrid)

      !Local grid size (different from global in line smoothing)
cc      nxl = MGgrid%nxv(igrid)
cc      nyl = MGgrid%nyv(igrid)
cc      nzl = MGgrid%nzv(igrid)
      nxl = mg_ctx_temp%nxv(igrid)
      nyl = mg_ctx_temp%nyv(igrid)
      nzl = mg_ctx_temp%nzv(igrid)

      !Find grid limits
      call limits(0,nxf,nyf,nzf,igrid,imin,imax,jmin,jmax,kmin,kmax)

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
      subroutine scheduleIndices(igrid,i1,j1,k1,i,j,k,zebra_dir)

c     -------------------------------------------------------------------
c     This routine schedules colored indices according to grid dimensions
c     and to the orderings in scheduleColors. In call:
c        * igrid: grid level
c        * i1,j1,k1: logical colored indices
c        * i,j,k: actual grid indices
c     -------------------------------------------------------------------

      implicit none

      integer    :: igrid,i1,j1,k1,i,j,k
      integer,optional :: zebra_dir

      integer    :: nxl,nyl,nzl,zb_dir

c     Begin program

      if (PRESENT(zebra_dir)) then
        zb_dir = zebra_dir
      else
        zb_dir = 0
      endif

      !Local grid size
cc      nxl = MGgrid%nxv(igrid)
cc      nyl = MGgrid%nyv(igrid)
cc      nzl = MGgrid%nzv(igrid)
      nxl = mg_ctx_temp%nxv(igrid)
      nyl = mg_ctx_temp%nyv(igrid)
      nzl = mg_ctx_temp%nzv(igrid)

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

      end module mg_internal

c module mgarraySetup
c ######################################################################
      module mgarraySetup

        use mg_internal

        type :: garray
          real(8),pointer,dimension(:,:,:,:) :: array
        end type garray

        type :: mg_array
          type(garray),pointer,dimension(:) :: grid
        end type mg_array

        logical :: is__cnv,is__vec,have_equl,res_is_vec

        type(mg_array) :: equl

        INTERFACE ASSIGNMENT (=)
          module procedure equateMGArray
        END INTERFACE

      contains

c     allocateMGArray
c     #################################################################
      subroutine allocateMGarray(neq,mgarray)

        implicit none

c     Call variables

        integer         :: neq
        type(mg_array)  :: mgarray

c     Local variables

        integer         :: igrid,nxp,nyp,nzp

c     Begin program

        if (.not.associated(mgarray%grid)) then
          allocate(mgarray%grid(grid_params%ngrid))
          do igrid=1,grid_params%ngrid
            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1
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
          do igrid=1,grid_params%ngrid
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

        integer          :: igrid,neq

c     Begin program

        neq = size(mg1%grid(1)%array,4)

        call allocateMGArray(neq,mg2)

        do igrid=1,grid_params%ngrid
          mg2%grid(igrid)%array = mg1%grid(igrid)%array
        enddo

c     End program

      end subroutine equateMGArray

c     restrictMGArray
c     #################################################################
      subroutine restrictMGArray(icmp,neq,mgarray,bcnd,igrid,order
     .                          ,iscnv,isvec,result_is_vec,equilb)
c     -----------------------------------------------------------------
c     Restricts MG array in all grids with ghost nodes.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer        :: neq,icmp,bcnd(6,neq),order,igrid
      type(mg_array) :: mgarray
      logical,optional,intent(IN) :: iscnv,isvec,result_is_vec
      type(mg_array),optional,intent(IN) :: equilb

c     Local variables

      integer        :: igf,nxf,nyf,nzf,igc,nxc,nyc,nzc

c     Begin program

      if (PRESENT(iscnv)) then
        is__cnv = iscnv
      else
        is__cnv = .true.   !Contravariant representation by default
      endif

      if (PRESENT(isvec)) then
        is__vec = isvec
      else
        is__vec = .true.   !Contravariant representation by default
      endif

      if (PRESENT(result_is_vec)) then
        res_is_vec = result_is_vec
      else
        res_is_vec = is__vec   !Contravariant representation by default
      endif

      if (PRESENT(equilb)) then
        have_equl = .true.
        equl = equilb
      else
        have_equl = .false.
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

      do igc=igrid+1,grid_params%ngrid
        igf = igc-1

        nxf = grid_params%nxv(igf)
        nyf = grid_params%nyv(igf)
        nzf = grid_params%nzv(igf)
        nxc = grid_params%nxv(igc)
        nyc = grid_params%nyv(igc)
        nzc = grid_params%nzv(igc)

        call restrictArrayToArray(icmp,neq
     .       ,igf,nxf,nyf,nzf,mgarray%grid(igf)%array
     .       ,igc,nxc,nyc,nzc,mgarray%grid(igc)%array
     .       ,order,.false.,bcnd)

      enddo

      if (have_equl) call deallocateMGArray(equl)

c     End program

      end subroutine restrictMGArray

c     restrictArrayToArray
c     #################################################################
      subroutine restrictArrayToArray(icmp,neq,igf,nxf,nyf,nzf,arrayf
     .                                        ,igc,nxc,nyc,nzc,arrayc
     .                               ,order,volf,bcnd)
c     -----------------------------------------------------------------
c     Restricts array to array in all grids (with ghost nodes),
c     starting at grid igf.
c     -----------------------------------------------------------------

      use setMGBC_interface
cc      use imposeBCinterface

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
     .             ,order,bcnd(6,neq),icmp
      real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
     .             ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)
      logical    :: volf

c     Local variables

      integer    :: igridf,igridc,isigf,isigc,i,j,k,ii,ieq
     .             ,nxxf,nyyf,nzzf,nxxc,nyyc,nzzc,ntotc,ntotf
     .             ,bcmod(6,neq),icmp2(max(neq/3,1))

      real(8),allocatable,dimension(:) :: vecf,vecc
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Consistency check

      nxxf = nxv(igf)
      nyyf = nyv(igf)
      nzzf = nzv(igf)

      if (nxf /= nxxf .or. nyf /= nyyf .or. nzf /= nzzf) then
        write (*,*) 'Grid mismatch in restrictArrayToArray:'
        write (*,*) 'Aborting...'
        stop
      endif

c     Allocate vectors

      nxxc = nxxf
      nyyc = nyyf
      nzzc = nzzf

      ntotf = neq*nxxf*nyyf*nzzf
      ntotc = ntotf

      allocate(vecf(ntotf))
      allocate(vecc(ntotc))

c     Map arrays onto MG vector

      !Set igrid=1 since vecf is NOT a MG vector
      call mapArrayToMGVector(neq,nxxf,nyyf,nzzf,arrayf,vecc,1)

c     Restrict MG vectors

      do igridc = igf+1,igc

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxxf = grid_params%nxv(igridf)
        nyyf = grid_params%nyv(igridf)
        nzzf = grid_params%nzv(igridf)

        nxxc = grid_params%nxv(igridc)
        nyyc = grid_params%nyv(igridc)
        nzzc = grid_params%nzv(igridc)

        ntotf = neq*nxxf*nyyf*nzzf
        ntotc = neq*nxxc*nyyc*nzzc

c       Allocate coarse mesh vector

        deallocate(vecf)
        allocate(vecf(ntotf))

        vecf = vecc

        deallocate(vecc)
        allocate(vecc(ntotc))

c       Restrict vector

        call crestrict(neq,vecc,ntotc,nxxc,nyyc,nzzc
     .                    ,vecf,ntotf,nxxf,nyyf,nzzf
cc     .                ,order,igridf,volf)
     .                ,0,igridf,volf)

      enddo

c     Map vector to array

      call mapMGVectorToArray(0,neq,vecc,nxc,nyc,nzc,arrayc,igc,.false.)

      if (icmp /= 0) then
        icmp2 = icmp
        if (have_equl) then
          call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
     .                ,arr0=equl%grid(igc)%array
     .                ,icomp=icmp2,is_cnv=is__cnv,is_vec=is__vec
     .                ,result_is_vec=res_is_vec,iorder=order)
        else
          bcmod = bcnd
          where (bcnd == EQU)
            bcmod = EXT
          end where
          call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcmod
     .              ,icomp=icmp2,is_cnv=is__cnv,is_vec=is__vec
     .              ,result_is_vec=res_is_vec,iorder=order)
cc     .              ,is_cnv=is__cnv)
        endif
      endif

c     Deallocate vectors

      deallocate(vecf,vecc)

      call deallocPointers(fpointers)

c     End program

      end subroutine restrictArrayToArray

c     restrictArrayToMGVector
c     #################################################################
      subroutine restrictArrayToMGVector(neq,nx,ny,nz,array,mgvector
     .                                  ,igr0,order,volf)
c     -----------------------------------------------------------------
c     Restricts array to mgvector in all grids (without ghost nodes),
c     starting at grid igr0.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,nx,ny,nz,order,igr0
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1,neq)
      real(8)    :: mgvector(*)
      logical    :: volf

c     Local variables

      integer    :: ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf
     .             ,isigf,isigc,ntotc,ntotf
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Consistency check

      nxf = nxv(igr0)
      nyf = nyv(igr0)
      nzf = nzv(igr0)

      if (nxf /= nx .or. nyf /= ny .or. nzf /= nz) then
        write (*,*) 'Grid mismatch in restrictArray:'
        write (*,*) 'Aborting...'
        stop
      endif

c     Map array in initial grid onto MG vector

      call mapArrayToMGVector(neq,nx,ny,nz,array,mgvector,igr0)

c     Restrict array to coarser grids

      do igridc = igr0+1,grid_params%ngrid

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxf = nxv(igridf)
        nyf = nyv(igridf)
        nzf = nzv(igridf)

        nxc = nxv(igridc)
        nyc = nyv(igridc)
        nzc = nzv(igridc)

        isigc = istart(igridc)
        isigf = istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Restrict MG vector

        call crestrict(neq,mgvector(isigc),ntotc,nxc,nyc,nzc
     .                    ,mgvector(isigf),ntotf,nxf,nyf,nzf
     .                ,order,igridf,volf)
      enddo

      call deallocPointers(fpointers)

      end subroutine restrictArrayToMGVector

c     restrictMGVector
c     #################################################################
      subroutine restrictMGVector(neq,mgvector,igr0,igr1,order,volf)
c     -----------------------------------------------------------------
c     Restricts vector to mgvector from grid igr0 to grid igr1.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,order,igr0,igr1
      real(8)    :: mgvector(*)
      logical    :: volf

c     Local variables

      integer    :: ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf
     .             ,isigf,isigc,ntotc,ntotf
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Consistency check

      nxf = nxv(igr0)
      nyf = nyv(igr0)
      nzf = nzv(igr0)

      if (igr0 > igr1) then
        call pstop('restrictVectorToMGVector'
     .            ,'Starting grid is coarser than final grid')
      endif

c     Restrict vector to coarser grids

      do igridc = igr0+1,igr1

        igridf = igridc-1

c       Characterize coarse and fine grids

        nxf = nxv(igridf)
        nyf = nyv(igridf)
        nzf = nzv(igridf)

        nxc = nxv(igridc)
        nyc = nyv(igridc)
        nzc = nzv(igridc)

        isigc = istart(igridc)
        isigf = istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Restrict MG vector

        call crestrict(neq,mgvector(isigc),ntotc,nxc,nyc,nzc
     .                    ,mgvector(isigf),ntotf,nxf,nyf,nzf
     .                ,order,igridf,volf)
      enddo

      call deallocPointers(fpointers)

      end subroutine restrictMGVector

c     prolongMGVector
c     #################################################################
      subroutine prolongMGVector(neq,mgvector,igr0,igr1,order,bcnd)
c     -----------------------------------------------------------------
c     Prolongs vector to mgvector from grid igr0 to grid igr1.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer    :: neq,order,igr0,igr1,bcnd(6,neq)
      real(8)    :: mgvector(*)

c     Local variables

      integer    :: ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf
     .             ,isigf,isigc,ntotc,ntotf
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

c     Consistency check

      nxc = nxv(igr0)
      nyc = nyv(igr0)
      nzc = nzv(igr0)

      if (igr0 < igr1) then
        call pstop('prolongMGVector'
     .            ,'Starting grid is coarser than final grid')
      endif

c     Prolong vector to coarser grids

      do igridf = igr0-1,igr1,-1

        igridc = igridf+1

c       Characterize coarse and fine grids

        nxf = nxv(igridf)
        nyf = nyv(igridf)
        nzf = nzv(igridf)

        nxc = nxv(igridc)
        nyc = nyv(igridc)
        nzc = nzv(igridc)

        isigc = istart(igridc)
        isigf = istart(igridf)

        ntotf = neq*nxf*nyf*nzf
        ntotc = neq*nxc*nyc*nzc

c       Prolong MG vector

        call cprolong(neq,mgvector(isigf),ntotf,nxf,nyf,nzf
     .                   ,mgvector(isigc),ntotc,nxc,nyc,nzc
     .               ,order,igridc,bcnd)
      enddo

      call deallocPointers(fpointers)

      end subroutine prolongMGVector

      end module mgarraySetup

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

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer    :: iter,igridmin,vcyc,crsedpth,isig,ncolors
      integer    :: orderres,orderprol,alloc_stat
     .             ,line_nsweep,line_crse_depth,line_vcyc

      integer    :: guess2,outc,mu
      integer    :: izero,i,j,k,ii,ivcyc,nblk

      real(8)    :: xx(2*ntot),yy(2*ntot),wrk(2*ntot),rr(ntot),x0(ntot)
      real(8)    :: rr0,rr1,mag,mag1,mag2,mgtol,line_tol,line_omega

      real(8)    :: epsmac=0d0,abstol

      real(8), pointer, dimension(:,:) :: diag

      logical    :: fdiag,line_relax,fpointers,line_x,line_y,line_z
     .             ,mg_debug

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
      volf            = options%vol_res
      crsedpth        = options%mg_coarse_solver_depth
      mu              = options%mg_mu
      vbr_mg          = options%mg_vertex_relax
      galerkin        = options%mg_galerkin
      MGgrid          => options%mg_grid_def
      diag            => options%diag
      mgctx           => options%mgctx

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

      igmax = ngrid - igridmin + 2  !Determines minimum resolution

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
     .                   ,ncolors,volf,galerkin)
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
        rr0 = sqrt(MGdot(igrid,neq,ntot,y,y))
      else
        call MGmatvec(0,neq,ntot,x,rr,igrid,bcnd,matvec)
        call MGvecadd(igrid,neq,ntot,-1d0,rr,1d0,y)
        rr0 = sqrt(MGdot(igrid,neq,ntot,rr,rr))
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

        call MGmatvec(0,neq,ntot,xx(istart(igrid)),rr,igrid,bcnd
     .               ,matvec)

        call MGvecadd(igrid,neq,ntot,-1d0,rr,1d0,y)

        mag = sqrt(MGdot(igrid,neq,ntot,rr,rr))

        mag1 = mag/rr1

        x0 = xx(istart(igrid):istart(igrid)+ntot-1)-x0
        mag2 = sqrt(MGdot(igrid,neq,ntot,x0,x0)/ntot)

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
          call MGmatvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd,matvec)
          call MGvecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))
          call MGplot(neq,wrk,igr,0,'fineres.bin')
        endif

        call smooth(igr)

        !Debugging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGplot(neq,xx ,igr,0,'finesol.bin')
        endif

c     Evaluate residual (ie wrk = yy - A xx = yy - wrk )

        wrk(isig:isig+nn-1) = 0d0

        call MGmatvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd,matvec)

        call MGvecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))

        !Debuging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGplot(neq,wrk,igr,1,'fineres.bin')
        endif

c     Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) to a coarser grid

        call crestrict(neq
     .                ,yy(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .                ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .                ,orderres,igr,volf)

        !Debugging plots
        if (igc == igmax.and.mg_debug) then
          write (*,*) 'Plotting here at level',igc
          call MGplot(neq,yy,igc,0,'coarse.bin')
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
          call MGplot(neq,xx,igc,1,'coarse.bin')

          call MGmatvec(0,neq,nnc,xx(isigc),wrk(isigc),igc,bcnd,matvec)
          call MGvecadd(igc,neq,nnc,-1d0,wrk(isigc),1d0,yy(isigc))
          call MGplot(neq,wrk,igc,1,'coarse.bin')
        endif

c     Cycle back up to grid igr updating errors (xx)
c     with fixed R.H.S. (yy)

c     Prolong error (wrk = P * xx_1) to a finer grid

        call cprolong(neq
     .               ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .               ,xx(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .               ,orderprol,igc,bcnd)

c     Update existing solution on grid igr (i.e. xx_igr): xx_igr = xx_igr + wrk 

        call MGvecadd(igr,neq,nn,1d0,xx(isig),1d0,wrk(isig))

        !Debugging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr
          call MGplot(neq,wrk,igr,1,'finesol.bin')
          call MGplot(neq,xx ,igr,1,'finesol.bin')

          call MGmatvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd,matvec)
          call MGvecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))
          call MGplot(neq,wrk,igr,1,'fineres.bin')
        endif

c     Relax updated solution on igr (i.e. xx_igr)

        call smooth(igr)

        !Debugging plots
        if (igr == igmax-1.and.mg_debug) then
          write (*,*) 'Plotting here at level',igr

          call MGmatvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd,matvec)
          call MGvecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))
          call MGplot(neq,wrk,igr,1,'fineres.bin')
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
     .                    ,guess2,outc,depth1)
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
c     with the original grid, defined in grid_params, and stored in
c     mg_ctx_temp).  However, the grid definition
c     (nxv,nyv,nzv,iline,jline,kline) contained in mg_ctx_temp is
c     "collapsed", and this is the way MG knows it should operate on a
c     restricted subspace of mesh points. Therefore, node positioning on
c     the grid follows grid_params, whereas matvec and smoothing
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

        logical        :: fpointers

c     Begin program

        nn   = ntotv (igr)
        isig = istart(igr)

cc        nxl = MGgrid%nxv
cc        nyl = MGgrid%nyv
cc        nzl = MGgrid%nzv
        nxl = mg_ctx_temp%nxv(igr)
        nyl = mg_ctx_temp%nyv(igr)
        nzl = mg_ctx_temp%nzv(igr)

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

          call MGmatvec(0,neq,nn,xx(isig),rr,igr,bcnd,matvec)
          call MGvecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

          if (outc > 0) then
            mag0 = sqrt(MGdot(igr,neq,nn,rr,rr))
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
                call lineSolver(nn,line_mg_ctx,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call MGvecadd(igr,neq,nn,1d0,xx(isig),omega,dx)

                !Recover grid pointers
                mg_ctx_temp => mg_ctx_sv
                nullify(mg_ctx_sv)

              enddo

              !Update residual
              call MGmatvec(0,neq,nn,xx(isig),rr,igr,bcnd,matvec)
              call MGvecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(MGdot(igr,neq,nn,rr,rr))
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
                call lineSolver(nn,line_mg_ctx,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call MGvecadd(igr,neq,nn,1d0,xx(isig),omega,dx)

                !Recover grid pointers
                mg_ctx_temp => mg_ctx_sv
                nullify(mg_ctx_sv)

              enddo

              !Update residual
              call MGmatvec(0,neq,nn,xx(isig),rr,igr,bcnd,matvec)
              call MGvecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(MGdot(igr,neq,nn,rr,rr))
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
                call lineSolver(nn,line_mg_ctx,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call MGvecadd(igr,neq,nn,1d0,xx(isig),omega,dx)

                !Recover grid pointers
                mg_ctx_temp => mg_ctx_sv
                nullify(mg_ctx_sv)

              enddo

              !Update residual
              call MGmatvec(0,neq,nn,xx(isig),rr,igr,bcnd,matvec)
              call MGvecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(MGdot(igr,neq,nn,rr,rr))
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
     .         ,guess2,outc,depth1)

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
      recursive subroutine lineSolver(nn,mg_grid,b,x,igr,guess)

c     ---------------------------------------------------------------------
c     This routine performs a recursive MG solve on selected planes/lines.
c     In call:
c       * nn: total grid size at grid level igr
c       * mg_grid: MG grid definition structure that defines local grid patch
c       * b: rhs vector
c       * x: solution vector
c       * igr: grid level
c       * guess: whether initial guess is provided (guess=1) or not (guess=0).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer        :: igr,nn,guess
        real(8)        :: b(nn),x(nn)
cc        type(grid_mg_def),pointer :: mg_grid
        type(mg_ctx),pointer :: mg_grid
        logical        :: prelax

c     Local variables

        integer        :: nxl,nyl,nzl

c     Begin program

        nxl = mg_grid%nxv(igr)
        nyl = mg_grid%nyv(igr)
        nzl = mg_grid%nzv(igr)

c     Schedule solvers

ccc       Recursive MG (at least 2D or if 1D MG is wanted)
cc        if  (    (nxl > 1 .and. nyl > 1)
cc     .       .or.(nxl > 1 .and. nzl > 1)
cc     .       .or.(nyl > 1 .and. nzl > 1)) then
cc
cc          call lineMGsolver(nn,mg_grid,b,x,igr,guess)
cc
ccc       GMRES 1D solver (matvec HARDWIRED to v_mtvc for now)
cccc        elseif (line_solve == "gm") then
cccc
cccc          call lineGMsolver(nn,mg_grid,b,x,igr,guess)
cc
ccc       MG 1D solver
cc        elseif (line_solve == "mg") then
cc
cc          call lineMGsolver(nn,mg_grid,b,x,igr,guess,oned_solve=.true.)
cc
ccc       GS 1D solver
cc        elseif (line_solve == "gs") then
cc
cc          call lineGSsolver(nn,mg_grid,b,x,igr,guess)
cc
ccc       JB 1D solver
cc        elseif (line_solve == "jb") then

          call lineJBsolver(nn,mg_grid,b,x,igr,guess)

cc        else
cc
cc          call pstop('lineSolver','Line solver undefined')
cc
cc        endif

c     End program

      end subroutine lineSolver

ccc     lineGMsolver
ccc     #####################################################################
cc      recursive subroutine lineGMsolver(nn,mg_grid,b,x,igr,guess)
cc
ccc     ---------------------------------------------------------------------
ccc     This routine performs a GMRES solve on selected lines.
ccc     In call:
ccc       * nn: total grid size at grid level igr
ccc       * mg_grid: MG grid definition structure that defines local grid patch
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
cc        type(grid_mg_def),pointer :: mg_grid
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
cc        nxl = mg_grid%nxv(igr)
cc        nyl = mg_grid%nyv(igr)
cc        nzl = mg_grid%nzv(igr)
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
ccc     Call GMRES (proxy routine lmtvc is defined in mg_internal module)
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
cc      recursive subroutine lineGSsolver(nn,mg_grid,b,x,igr,guess)
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
ccc     ---------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer        :: igr,nn,guess
cc        real(8)        :: b(nn),x(nn)
cccc        type(grid_mg_def),pointer :: mg_grid
cc        type(mg_ctx),pointer :: mg_grid
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
ccc     Save parent MG grid configuration (shared with filial MG via mg_internal module)
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
      recursive subroutine lineJBsolver(nn,mg_grid,b,x,igr,guess)

c     ---------------------------------------------------------------------
c     This routine performs a recursive JB solve on selected planes/lines.
c     In call:
c       * nn: total grid size at grid level igr
c       * mg_grid: MG grid definition structure that defines local grid patch
c       * b: rhs vector
c       * x: solution vector
c       * igr: grid level
c       * guess: whether initial guess is provided (guess=1) or not (guess=0).
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer        :: igr,nn,guess
        real(8)        :: b(nn),x(nn)
        type(mg_ctx),pointer :: mg_grid
        logical        :: prelax

c     Local variables

        type (solver_options) :: options

c     Begin program

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
        options%vol_res                     = volf        

        options%mg_coarse_solver_depth      = line_crse_depth
        options%mg_line_relax               = .false.
        options%mg_line_nsweep              = line_nsweep
        options%mg_line_vcyc                = line_vcyc
        options%mg_line_tol                 = line_tol
        options%mg_line_omega               = line_omega
        options%mg_line_coarse_solver_depth = line_crse_depth

        options%mg_vertex_relax             = vbr_mg

        options%mgctx                       => mg_grid
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
ccc     Save parent MG grid configuration (shared with filial MG via mg_internal module)
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

        integer    :: igr

c     Local variables

        integer    :: nn,isig,depth1

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

c     coarseGridSolve
c     #####################################################################
      recursive subroutine coarseGridSolve(igr)

c     ---------------------------------------------------------------------
c     This routine performs a recursive MG solve on parallel coarse grid.
c     In call:
c       * igr: grid level
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer        :: igr

c     Local variables

        integer        :: nn,isig,nx,ny,nz,nc,i,j,k,il,jl,kl
     .                   ,ii,iig,guess,istart_sv(size(istart))
     .                   ,ig_cs,ig_sv,isig_cs,isig_sv,nsv,ncs

        real(8),allocatable,dimension(:,:) :: diag_sv

        real(8),allocatable,dimension(:)   :: b_global,x_global

        type(grid_mg_def),pointer :: gparams_sv

        type (solver_options) :: options

        logical               :: fpointers

c     Begin program

        nn   = ntotv(igr)
        isig = istart(igr)

c     Save parent MG grid configuration (shared with filial MG via mg_internal module)

        !Save pointers
        MGgrid_sv%ngrdx = 0
        MGgrid_sv%ngrdy = 0
        MGgrid_sv%ngrdz = 0
        MGgrid_sv%ngrid = 0
cc        MGgrid_sv  = MGgrid
        call equateGridStructure(MGgrid_sv,MGgrid)

        call deallocateGridStructure(MGgrid)

        gparams_sv%ngrdx = 0
        gparams_sv%ngrdy = 0
        gparams_sv%ngrdz = 0
        gparams_sv%ngrid = 0
cc        gparams_sv  = grid_params
        call equateGridStructure(gparams_sv,grid_params)

        call deallocateGridStructure(grid_params)

        istart_sv  = istart

        !Save diagonal
        allocate(diag_sv(size(diag,1),size(diag,2)))
        diag_sv = diag
        deallocate(diag)

c     Find coarse global grid dimensions

        nx = MGgrid_sv%nxgl(igr)
        ny = MGgrid_sv%nygl(igr)
        nz = MGgrid_sv%nzgl(igr)

        nc = nx*ny*nz*neq

c     Ensure serial treatment of boundary conditions

        !MISSING STUFF HERE

c     Create coarse global grid structure variable to pass to filial MG call

        call createMGgrid(grid_params)  !May be wrong interface...
          
c     Reallocate MG pointers to coarse grid level igr

        call deallocPointers(.true.)
        call allocPointers(neq,fpointers)

c     Allocate coarse global quantities

        allocate(b_global(nc),x_global(nc))

        x_global = 0d0

c     Gather rhs b in all processors

#if defined(petsc)
        call MPI_Allgather(yy(isig),nn,MPI_DOUBLE_PRECISION
     .                    ,b_global,nn,MPI_DOUBLE_PRECISION
     .                    ,MPI_COMM_WORLD,mpierr)
c diag ****
cc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cc     .             ,'local rhs =',yy(isig:isig+nn-1)

cc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cc     .             ,'global rhs=',b_global
cc        call pstop('coarseGridSolve','')
c diag ****
#endif

ccc     Allocate new diagonal and transfer elements
cc
cc        allocate(diag(neq,2*nc))
cc
ccc     Gather diag in all processors and all grids
cc
cc#if defined(petsc)
cc        do ig_sv=igr,gparams_sv%ngrid
cc          ig_cs   = ig_sv - igr + 1
cc          isig_sv = istart_sv(ig_sv)
cc          isig_cs = istart(ig_cs)
cc          nsv     = istart_sv(ig_sv+1)-istart_sv(ig_sv)
cc          ncs     = istart   (ig_cs+1)-istart   (ig_cs)
cc
ccc diag ****
cc          write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cc     .               ,'global grid=',ig_sv,'local grid=',ig_cs
cc     .               ,'size vector transferred=',nsv
cc     .               ,'size vector received=',ncs
cc     .               ,'global pointers=',istart_sv
cc     .               ,'local pointers=',istart
ccc diag ****
cc          call MPI_Allgather(diag_sv(:,isig_sv:isig_sv+nsv-1)
cc     .                              ,neq*nsv,MPI_DOUBLE_PRECISION
cc     .                      ,diag   (:,isig_cs:isig_cs+ncs-1)
cc     .                              ,neq*nsv,MPI_DOUBLE_PRECISION
cc     .                      ,MPI_COMM_WORLD,mpierr)
cc        enddo
cc
ccc diag ****
cc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cc     .             ,'local diag =',diag_sv(1,isig:isig+nn-1)
cc
cccc        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
cccc     .             ,'global diag=',diag(1,1:nc)
cc        call pstop('coarseGridSolve','')
ccc diag ****
cc#endif

c     Configure recursive plane/line MG solve

        options%tol                         = mgtol
        options%vcyc                        = vcyc
        options%mg_gmin                     = igridmin
        options%orderres                    = orderres
        options%orderprol                   = orderprol
        options%mg_mu                       = mu
        options%vol_res                     = volf        
        options%mg_line_relax               = .false.
        options%mg_coarse_solver_depth      = crsedpth 
        options%mg_vertex_relax             = vbr_mg
        options%mg_galerkin                 = galerkin
        options%ncolors                     = ncolors
cc        options%mg_grid_def                 = grid_params
        call equateGridStructure(options%mg_grid_def,grid_params)

        options%fdiag                       = .true.

c     Call plane/line MG

c diag ****
        if (my_rank == 0) call writeGridStructure(grid_params)
        write (*,*) size(b_global),size(x_global),nc,neq
        call pstop('coarseGridSolve','')
c diag ****

        guess = 0
        call mg(neq,nc,b_global,x_global,matvec,options,1,bcnd
     .         ,guess,outc-1,depth)

c     Scatter x_global

        do k = MGgrid_sv%klo(igr),MGgrid_sv%khi(igr)
          do j = MGgrid_sv%jlo(igr),MGgrid_sv%jhi(igr)
            do i = MGgrid_sv%ilo(igr),MGgrid_sv%ihi(igr)
              call fromGlobalToLocalLimits(i,j,k,il,jl,kl,igr,igr,igr)
              ii  = vecPos(neq,il,jl,kl,igr,igr,igr) + isig - 1
              iig = neq*(i-1 + nx*(j-1) + nx*ny*(k-1))

              xx(ii+1:ii+neq) = x_global(iig+1:iig+neq)
            enddo
          enddo
        enddo

c diag ****
        write (*,*) 'DIAG -- coarseGridSolve: proc=',my_rank
     .             ,xx(isig:isig+nn)
        call pstop('coarseGridSolve','')
c diag ****

c     Recover configuration for parent level MG

        !Recover MG grid def
        call restoreMGgrid

        !Recover grid_params
        call deallocateGridStructure(grid_params)
cc        grid_params = gparams_sv
        call equateGridStructure(grid_params,gparams_sv)
        call deallocateGridStructure(gparams_sv)

        !Recover pointers
        call deallocPointers(fpointers)
        call allocPointers(neq,fpointers)

        !Recover diagonal
        if (.not.options%fdiag) deallocate(diag)
        allocate(diag(size(diag_sv,1),size(diag_sv,2)))
        diag = diag_sv
        deallocate(diag_sv)

c     Deallocate memory

        call deallocateGridStructure(options%mg_grid_def)

      end subroutine coarseGridSolve

c     killmg
c     ###################################################################
      recursive subroutine killmg

        implicit none

        if (fdiag) then
          deallocate(diag)
        else
          nullify(diag)
        endif

        call deallocPointers(fpointers)
cc        call deallocateGridStructure(MGgrid)
        nullify(MGgrid)

      end subroutine killmg

      end subroutine mg

*deck cprolong
c######################################################################
      subroutine cprolong(neq,xf,ntotf,nxf,nyf,nzf
     .                       ,xc,ntotc,nxc,nyc,nzc
     .                   ,order,igc,bcnd)
c----------------------------------------------------------------------
c     This is a prolongation routine for system MG, with arbitrary order
c     of interpolation.
c
c     In call sequence, we have:
c       * neq (int): number of equations
c       * xf (real): vector in fine grid (not a MG vector)
c       * ntotf (int): dimension of xf.
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * xc (real): vector in coarse grid (not a MG vector)
c       * ntotc (int): dimension of xc.
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igc (int): coarse grid level identifier
c       * bcnd (int array): boundary condition info.
c----------------------------------------------------------------------

      use mg_internal

      use setMGBC_interface

      implicit none            ! For safe Fortran

c Call variables

      integer    :: neq,ntotc,nxc,nyc,nzc,ntotf,nxf,nyf,nzf,order,igc
     .             ,bcnd(6,neq)
      real(8)    :: xc(ntotc),xf(ntotf)

c Local variables
 
      real(8)    :: xxf(ntotf/neq,neq),car(3)
     .             ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)

      real(8),allocatable :: arrayf(:,:,:,:),sp_sum(:,:)
     .                      ,dv_sum(:,:),l_sum(:,:)

      integer    :: ic,jc,if,jf,iic,iif,i,j,k,ig,jg,kg
     .             ,ieq,nntotc,nntotf,igf

      logical    :: cartesian,isvec,iscnv

c Begin program

      nntotc = ntotc/neq
      nntotf = ntotf/neq

c Unpack vector into array

      call mapMGVectorToArray(0,neq,xc,nxc,nyc,nzc,arrayc,igc,.false.)

c Impose boundary conditions (external)

      isvec = (neq > 1)                   !Initialize isvec (may be modified by setMGBC)

      cartesian = (bcSP() .and. neq > 1)  !Return Cartesian variables for prolongation

      if (order > 0 .or. cartesian) then
        call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
     .              ,iorder=min(order,3)
     .              ,result_is_vec=(.not.cartesian)
     .              ,out_isvec=isvec,out_iscnv=iscnv)  !<= These define whether quantities are
                                                       !   vectors or scalars.
      endif

c Perform average along azimuthal direction for SP configurations

      if (isSP2(1,igc) .and. isvec) then

        allocate(sp_sum(0:nzc+1,neq),dv_sum(0:nzc+1,neq))

        do ieq=1,neq
          do k=0,nzc+1
            sp_sum(k,ieq)= sum(arrayc(1,1:nyc,k,ieq)
     .                        *gmetric%grid(igc)%dvol(1,1:nyc,k))
            dv_sum(k,ieq)= sum(gmetric%grid(igc)%dvol(1,1:nyc,k))
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

c Prolong arrays (use scalar prolongation)

      xxf = 0d0

      do ieq=1,neq

        call prolong(xxf   (:    ,ieq),nxf,nyf,nzf
     .              ,arrayc(:,:,:,ieq),nxc,nyc,nzc
     .              ,order,igc)

cc        !Repack prolonged vector
cc        do i=1,nntotf
cc          xf(neq*(i-1)+ieq)=xxf(i)
cc        enddo

      enddo

c Transform back to curvilinear (if isvec=.true. and cartesian=.true.)

      if (isvec .and. cartesian) then

        igf = igc - 1

        do k = 1,nzf
          do j = 1,nyf
            do i = 1,nxf
              do ieq=1,neq/3
                iif = vecPos(1,i,j,k,igf,igf,igf)+1

                car = xxf(iif,(ieq-1)*3+1:ieq*3)

                call transformVectorToCurvilinear
     .             (i,j,k,igf,igf,igf
     .             ,car(1),car(2),car(3)
     .             ,iscnv
     .             ,xxf(iif,(ieq-1)*3+1)
     .             ,xxf(iif,(ieq-1)*3+2)
     .             ,xxf(iif,(ieq-1)*3+3))
              enddo
            enddo
          enddo
        enddo

      endif

c Repack prolonged vector

      do ieq=1,neq
        do i=1,nntotf
          xf(neq*(i-1)+ieq)=xxf(i,ieq)
        enddo
      enddo

c End program

      end subroutine cprolong

*deck prolong
c######################################################################
      subroutine prolong(xf,nxf,nyf,nzf,arrayc,nxc,nyc,nzc,order,igc)
c----------------------------------------------------------------------
c     This is a prolongation routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * xf (real): vector in fine grid
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * arrayc (real): array of values in coarse grid
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igc (int): coarse grid level identifier
c----------------------------------------------------------------------

      use mg_internal

      implicit none            ! For safe Fortran

c Call variables

      integer    :: nxc,nyc,nzc,nxf,nyf,nzf,order,igc
      real(8)    :: arrayc(0:nxc+1,0:nyc+1,0:nzc+1),xf(nxf*nyf*nzf)

c Local variables

      integer    :: ic,if,jc,jf,kc,kf,igf,iic,iif
     .             ,icg,jcg,kcg,ifg,jfg,kfg
      integer    :: iminc,imaxc,jminc,jmaxc,kminc,kmaxc
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf
      real(8)    :: mag,vol

      logical    :: fpointers

c Extrapolation

      real(8)    :: xxf,yyf,zzf,ff

      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer    :: kx,ky,kz,nx,ny,nz,dim,flg
      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db3val
      external      db3val

c Begin program

      call allocPointers(1,fpointers)

c Define fine grid

      igf = igc - 1

      call limits(0,nxc,nyc,nzc,igc,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

c Injection

      if (order.eq.0) then

        do kc = kminc,kmaxc
          do jc = jminc,jmaxc
            do ic = iminc,imaxc

              do kf = mg_ratio_z(igf)*(kc-1)+1,mg_ratio_z(igf)*kc
                do jf = mg_ratio_y(igf)*(jc-1)+1,mg_ratio_y(igf)*jc
                  do if = mg_ratio_x(igf)*(ic-1)+1,mg_ratio_x(igf)*ic
                    iif = if + nxf*(jf-1) + nxf*nyf*(kf-1)
                    xf(iif) = arrayc(ic,jc,kc)
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

        call getMGmap(iminc,jminc,kminc,igc,igc,igc,icg,jcg,kcg)
        xx = grid_params%xx(icg-1:icg+imaxc-iminc+1)
        yy = grid_params%yy(jcg-1:jcg+jmaxc-jminc+1)
        zz = grid_params%zz(kcg-1:kcg+kmaxc-kminc+1)

c     Prepare 3d spline interpolation

        flg = 0
        nx = imaxc-iminc+3
        ny = jmaxc-jminc+3
        nz = kmaxc-kminc+3
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(tx(nx+kx))
        allocate(ty(ny+ky))
        allocate(tz(nz+kz))
        allocate(work(dim))
        allocate(bcoef(nx,ny,nz))

        call db3ink(xx,nx,yy,ny,zz,nz
     .        ,arrayc(iminc-1:imaxc+1,jminc-1:jmaxc+1,kminc-1:kmaxc+1)
     .        ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)

c     Interpolate

        call limits(0,nxf,nyf,nzf,igf
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf)

        do kf = kminf,kmaxf
          do jf = jminf,jmaxf
            do if = iminf,imaxf
              iif = if + nxf*(jf-1) + nxf*nyf*(kf-1)

              call getMGmap(if,jf,kf,igf,igf,igf,ifg,jfg,kfg)
              xxf = grid_params%xx(ifg)
              yyf = grid_params%yy(jfg)
              zzf = grid_params%zz(kfg)

              xf(iif) = db3val(xxf,yyf,zzf,0,0,0,tx,ty,tz,nx,ny,nz
     .                        ,kx,ky,kz,bcoef,work)
            enddo
          enddo
        enddo

        deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

      endif

ccc Postprocessing for SP systems
cc
cc      call limits(0,nxf,nyf,nzf,igf
cc     .           ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf)
cc
cc      if (isSP(iminf,jminf,kminf,igf,igf,igf)) then
cc
cc        if = iminf
cc
cc        do kf = kminf,kmaxf
cc          mag = 0d0
cc          vol = 0d0
cc          do jf = jminf,jmaxf
cc            iif = if + nxf*(jf-1) + nxf*nyf*(kf-1)
cc
cc            mag = mag + gmetric%grid(igf)%dvol(if,jf,kf)*xf(iif)
cc            vol = vol + gmetric%grid(igf)%dvol(if,jf,kf)
cc          enddo
cc          xf(iif) = mag/vol
cc        enddo
cc
cc      endif

c End program

      call deallocPointers(fpointers)

      end subroutine prolong

*deck crestrict
c######################################################################
      subroutine crestrict(neq,xc,ntotc,nxc,nyc,nzc
     .                        ,xf,ntotf,nxf,nyf,nzf
     .                    ,order,igf,volw)
c----------------------------------------------------------------------
c     This is a restriction routine for system MG, with arbitrary order
c     of interpolation.
c
c     In call sequence, we have:
c       * neq (int): number of equations
c       * xc (real): vector in coarse grid
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * xf (real): vector in fine grid
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igf (int): fine grid level identifier
c       * bcnd (int array): boundary condition info.
c       * volw (logical): whether vectors contain volume fractions.
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer    :: ntotc,nxc,nyc,nzc,ntotf,nxf,nyf,nzf,order,igf,neq

      real(8)    :: xc(ntotc),xf(ntotf)

      logical    :: volw

c Local variables
 
      real(8)    :: xxc(ntotc/neq),xxf(ntotf/neq)

      integer    :: ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

      nntotc = ntotc/neq
      nntotf = ntotf/neq

c Restrict arrays

      do ieq = 1,neq

c     Unpack vector

        do i = 1,nntotf
          xxf(i) = xf(neq*(i-1)+ieq)
        enddo

c     Use scalar restriction (xxf -> xxc)

        xxc = 0d0
        call restrict(xxc,nxc,nyc,nzc
     .               ,xxf,nxf,nyf,nzf
     .               ,order,igf,volw)

c     Repack restricted vector

        do i=1,nntotc
          xc(neq*(i-1)+ieq)=xxc(i)
        enddo

      enddo

c End program

      end subroutine crestrict

*deck restrict
c######################################################################
      subroutine restrict(xc,nxc,nyc,nzc,xf,nxf,nyf,nzf,order,igf,volw)
c----------------------------------------------------------------------
c     This is a restriction routine for a single quantity, with
c     arbitrary order of interpolation.
c
c     In call sequence, we have:
c       * xc (real): vector in coarse grid
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * xf (real): vector in fine grid
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igf (int): fine grid level identifier
c       * volw (logical): whether values are volume-weighed.
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer    :: nxc,nyc,nzc,nxf,nyf,nzf,order,igf

      real(8)    :: xc(nxc*nyc*nzc),xf(nxf*nyf*nzf)

      logical    :: volw

c Local variables
 
      integer    :: ic,if,jc,jf,kc,kf,iif,iic,igc,inbv
     .             ,icg,jcg,kcg,ifg,jfg,kfg
      integer    :: iminc,imaxc,jminc,jmaxc,kminc,kmaxc
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf

      real(8)    :: volt,vol,arrayf(0:nxf+1,0:nyf+1,0:nzf+1)

      logical    :: fpointers

c Interpolation

      real(8)    :: xxc,yyc,zzc,ff

      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer    ::  kx,ky,kz,nx,ny,nz,dim,flg
      real(8), dimension(:),allocatable:: tx,ty,tz,work,q
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: dbvalu,db2val,db3val
      external      dbvalu,db2val,db3val

c Begin program

      call allocPointers(1,fpointers)

      igc = igf + 1

c Agglomeration

      if (order.eq.0) then

        call limits(0,nxc,nyc,nzc,igc
     .             ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

        do kc = kminc,kmaxc
          do jc = jminc,jmaxc
            do ic = iminc,imaxc

              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)

              volt   = 0d0
              xc(iic)= 0d0
              do kf = mg_ratio_z(igf)*(kc-1)+1,mg_ratio_z(igf)*kc
                do jf = mg_ratio_y(igf)*(jc-1)+1,mg_ratio_y(igf)*jc
                  do if = mg_ratio_x(igf)*(ic-1)+1,mg_ratio_x(igf)*ic
                    iif = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                    if (.not.volw) then
                      vol = gmetric%grid(igf)%dvol(if,jf,kf)
                      xc(iic) = xc(iic) + xf(iif)*vol
                      volt = volt + vol
                    else
                      xc(iic) = xc(iic) + xf(iif)
                    endif

                  enddo
                enddo
              enddo

              if (.not.volw) xc(iic) = xc(iic)/volt

            enddo
          enddo
        enddo

c Interpolation

      else

        call limits(0,nxf,nyf,nzf,igf
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf)

c     Map vectors into arrays

        allocate(xx(1:imaxf-iminf+1)
     .          ,yy(1:jmaxf-jminf+1)
     .          ,zz(1:kmaxf-kminf+1))

        call getMGmap(iminf,jminf,kminf,igf,igf,igf,ifg,jfg,kfg)

        xx = grid_params%xx(ifg:ifg+imaxf-iminf)
        yy = grid_params%yy(jfg:jfg+jmaxf-jminf)
        zz = grid_params%zz(kfg:kfg+kmaxf-kminf)

        call mapMGVectorToArray(0,1,xf,nxf,nyf,nzf,arrayf,igf,.false.)

c     Renormalize without volume fractions

        if (volw) then
          do kf = kminf,kmaxf
            do jf = jminf,jmaxf
              do if = iminf,imaxf

                call getMGmap(if,jf,kf,igf,igf,igf,ifg,jfg,kfg)

                arrayf(if,jf,kf) = arrayf(if,jf,kf)
     .                            /gmetric%grid(igf)%dvol(if,jf,kf)
cc     .                            /volume(if,jf,kf,igf,igf,igf)
cc     .                            /(grid_params%dxh(ifg)
cc     .                             *grid_params%dyh(jfg)
cc     .                             *grid_params%dzh(kfg))

              enddo
            enddo
          enddo
        endif

c     Calculate interpolation

        flg = 0
        nx = imaxf-iminf+1
        ny = jmaxf-jminf+1
        nz = kmaxf-kminf+1
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(tx(nx+kx))
        allocate(ty(ny+ky))
        allocate(tz(nz+kz))
        allocate(work(dim))
        allocate(bcoef(nx,ny,nz))

        if (nx == 1 .and. ny == 1) then
          allocate(q((2*kz-1)*nz))
          inbv=1
          call dbknot(zz,nz,kz,tz)
          call dbintk(zz,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,tz,nz,kz,bcoef(1,1,:),q,work)
          deallocate(q)
        elseif (nx == 1 .and. nz == 1) then
          allocate(q((2*ky-1)*ny))
          inbv=1
          call dbknot(yy,ny,ky,ty)
          call dbintk(yy,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,ty,ny,ky,bcoef(1,:,1),q,work)
          deallocate(q)
        elseif (ny == 1 .and. nz == 1) then
          allocate(q((2*kx-1)*nx))
          inbv=1
          call dbknot(xx,nx,kx,tx)
          call dbintk(xx,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,tx,nx,kx,bcoef(:,1,1),q,work)
          deallocate(q)
        elseif (nx == 1) then
          call db2ink(yy,ny,zz,nz
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,ny,ky,kz,ty,tz,bcoef(1,:,:),work,flg)
        elseif (ny == 1) then
          call db2ink(xx,nx,zz,nz
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,nx,kx,kz,tx,tz,bcoef(:,1,:),work,flg)
        elseif (nz == 1) then
          call db2ink(xx,nx,yy,ny
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,nx,kx,ky,tx,ty,bcoef(:,:,1),work,flg)
        else
          call db3ink(xx,nx,yy,ny,zz,nz
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)
        endif

        call limits(0,nxc,nyc,nzc,igc
     .             ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

        do kc = kminc,kmaxc
          do jc = jminc,jmaxc
            do ic = iminc,imaxc

              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)

              call getMGmap(ic,jc,kc,igc,igc,igc,icg,jcg,kcg)
              xxc = grid_params%xx(icg)
              yyc = grid_params%yy(jcg)
              zzc = grid_params%zz(kcg)

              if (nx == 1 .and. ny == 1) then
                xc(iic) = dbvalu(tz,bcoef(1,1,:),nz,kz,0,zzc,inbv,work)
              elseif (nx == 1 .and. nz == 1) then
                xc(iic) = dbvalu(ty,bcoef(1,:,1),ny,ky,0,yyc,inbv,work)
              elseif (ny == 1 .and. nz == 1) then
                xc(iic) = dbvalu(tx,bcoef(:,1,1),nx,kx,0,xxc,inbv,work)
              elseif (nx == 1) then
                xc(iic)=db2val(yyc,zzc,0,0,ty,tz,ny,nz,ky,kz
     .                        ,bcoef(1,:,:),work)
              elseif (ny == 1) then
                xc(iic)=db2val(xxc,zzc,0,0,tx,tz,nx,nz,kx,kz
     .                        ,bcoef(:,1,:),work)
              elseif (nz == 1) then
                xc(iic)=db2val(xxc,yyc,0,0,tx,ty,nx,ny,kx,ky
     .                        ,bcoef(:,:,1),work)
              else
                xc(iic)=db3val(xxc,yyc,zzc,0,0,0,tx,ty,tz,nx,ny,nz
     .                        ,kx,ky,kz,bcoef,work)
              endif

              if (volw) xc(iic)=xc(iic)*gmetric%grid(igc)%dvol(ic,jc,kc)

            enddo
          enddo
        enddo

        deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

      endif

c End program

      call deallocPointers(fpointers)

      end subroutine restrict

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

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer    :: iter,alloc_stat,isig,itr,nn,ieq,irbg1,nrbg1,igridc
      integer    :: imin,imax,jmin,jmax,kmin,kmax,iming,itr2
      real(8)    :: omega0,omega10,omega01,tol
      logical    :: fdiag,fpointers,vbr

      integer    :: i,j,k,iv,jv,kv,if,jf,kf,nxf,nyf,nzf,iiv,ivg,ncolors
      integer    :: ii,iii,iib,iiib,iic,iig,iblock,jblock,kblock,nblk,i1
      real(8)    :: mag0=1d0,mag1,mag,lmag,yy(ntot),magsp,magsp0

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

c Allocate pointers

      call allocPointers(neq,fpointers,mgctx=mgctx)

      if (fpointers) then          !JB is NOT called from MG
        igmax = igrid
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2.and.prnt)
     .     write (*,*) 'Allocating pointers...'

c Initialize auxiliary variables

      !Global dimensions
cc      nxf = nxv(igrid)
cc      nyf = nyv(igrid)
cc      nzf = nzv(igrid)

      nxf = grid_params%nxv(igrid)
      nyf = grid_params%nyv(igrid)
      nzf = grid_params%nzv(igrid)

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
     .                   ,ncolors,volf,galerkin)
        if (out.ge.2.and.prnt) write (*,*) 'Finished!'
        fdiag = .true.
      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c Jacobi iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag  = 0d0

        call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)

c       VERTEX-BASED RELAXATION
        if (vbr) then

          !Vertex sampling
          do kv=1,max(nzf-mg_ratio_z(igrid)+1,1)
            do jv=1,max(nyf-mg_ratio_y(igrid)+1,1)
              do iv=1,max(nxf-mg_ratio_x(igrid)+1,1)

              !Find new residual
              do kblock=1,mg_ratio_z(igrid)
                do jblock=1,mg_ratio_y(igrid)
                  do iblock=1,mg_ratio_x(igrid)

                    iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                           + mg_ratio_x(igrid)
     .                            *mg_ratio_y(igrid)*(kblock-1)

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
     .                           + mg_ratio_x(igrid)
     .                            *mg_ratio_y(igrid)*(kblock-1)

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

c       STANDARD RELAXATION
        else

          call limits(0,nxf,nyf,nzf,igrid,imin,imax,jmin,jmax,kmin,kmax)

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                ii = i + nxf*(j-1) + nxf*nyf*(k-1)

                iii = neq*(ii-1)

                !Find new residual
                do ieq = 1,neq
                  rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
                enddo

                !Multiply by D^-1 (stored in diag)
                iig = iii + isig - 1
                dummy = matmul(transpose(diag(:,iig+1:iig+neq)),rhs)

                !Update zz
                do ieq=1,neq
                  zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
                enddo

                mag = mag + sum(rhs*rhs)

              enddo
            enddo
          enddo

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

        endif

c     Check convergence

#if defined(petsc)
cc        if (.not.asm .and. out >= 1) then  !Do reduction only for debugging
        if (.not.asm) then
          lmag = mag
          call MPI_Allreduce(lmag,mag,1,MPI_DOUBLE_PRECISION
     .                       ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc        else
cc          mag = mag0                 !Do NOT do reduction in parallel (convergence is set by #its)
        endif
#endif

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          mag1 = mag
          if (mag.lt.1d-20*nn) then
            mag0 = 1d0
            mag1 = 1d0
            exit
          endif
          if (out.ge.2.and.prnt) write (*,10)itr-1,mag,mag/mag0
        else
          if (out.ge.2.and.prnt) write (*,20)itr-1,mag,mag/mag0,mag/mag1
          mag1 = mag
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.1) then
        call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)
        call MGvecadd(igrid,neq,ntot,-1d0,yy,1d0,rr)  !yy=rr-yy
        mag = sqrt(MGdot(igrid,neq,ntot,yy,yy))       !sqrt(yy*yy)

        if (itr.eq.0) then
          mag0 = mag
          mag1 = mag
        endif
      endif

      if (out.ge.2.and.prnt) then
        write (*,20) itr,mag,mag/mag0,mag/mag1
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

      call deallocPointers(fpointers)

      return

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

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer    :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer    :: iter,alloc_stat,isig,itr,itr2,nn,ieq,iig
     .             ,zebra_dir,minzd,maxzd,zebra_it,zd,zit
      integer    :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
      integer    :: i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max

      real(8)    :: omega0,tol,zebra_omg
      logical    :: fdiag,fpointers,vbr,zebra

      integer    :: ic,jc,kc,iv,jv,kv,if,jf,kf,nxf,nyf,nzf
      integer    :: ii,iii,iib,iiib,iiv,ivg,iblock,jblock,kblock,nblk
      real(8)    :: mag0=1d0,mag1,magsp,magsp0,mag,lmag,yy(ntot)

      integer    :: ncolors,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3

      real(8), pointer, dimension(:,:) :: diag

      real(8),allocatable, dimension(:) :: dummy,rhs

cc      type(grid_mg_def),pointer :: line_mg_grid

      type(mg_ctx),pointer :: mgctx

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
        igmax = igrid
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2.and.prnt)
     .     write (*,*) 'Allocated pointers.'

c Initialize auxiliary variables

      !Global grid size
cc      nxf = nxv(igrid)
cc      nyf = nyv(igrid)
cc      nzf = nzv(igrid)
      nxf = grid_params%nxv(igrid)
      nyf = grid_params%nyv(igrid)
      nzf = grid_params%nzv(igrid)

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
     .                   ,ncolors,volf,galerkin)
        if (out.ge.2.and.prnt) write (*,*) 'Finished!'
        fdiag = .true.
      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c GS iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag  = 0d0

c       ----------
c       Colored GS
c       ---------- 
        if (ncolors > 1) then

c         VERTEX-BASED RELAXATION
          if (vbr) then

            gs_type = 'C-VERTX'

            call scheduleColors(igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                         ,i1min,i1max,j1min,j1max,k1min,k1max)

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)

                !Vertex sampling
                do kv=1+mod((    irbg3-1),nrbg3)
     .                       ,max(nzf-mg_ratio_z(igrid)+1,1),nrbg3
                  do jv=1+mod((kv   +irbg2+irbg1-1),nrbg2)
     .                         ,max(nyf-mg_ratio_y(igrid)+1,1),nrbg2
                    do iv=1+mod((jv+kv+irbg1-1),nrbg1)
     .                           ,max(nxf-mg_ratio_x(igrid)+1,1),nrbg1

                     !Find new residual
                      do kblock=1,mg_ratio_z(igrid)
                        do jblock=1,mg_ratio_y(igrid)
                          do iblock=1,mg_ratio_x(igrid)

                            iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                                   + mg_ratio_x(igrid)
     .                                    *mg_ratio_y(igrid)*(kblock-1)

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
     .                            ,rhs)

                      !Update zz
                      do kblock=1,mg_ratio_z(igrid)
                        do jblock=1,mg_ratio_y(igrid)
                          do iblock=1,mg_ratio_x(igrid)

                            iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                                   + mg_ratio_x(igrid)
     .                                    *mg_ratio_y(igrid)*(kblock-1)

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

c         ZEBRA LINE RELAXATION
          elseif (zebra) then

            gs_type = 'C-ZEBRA'

            if (zebra_dir == 0) then
              minzd = 1
              maxzd = 3
            else
              minzd = zebra_dir
              maxzd = zebra_dir
            endif

            do zd=minzd,maxzd  !Cycle in zebra directions

              call scheduleColors(igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                         ,i1min,i1max,j1min,j1max,k1min,k1max
     .                         ,zebra_dir=zd)

              if (out.ge.2.and.prnt) then
                if (nrbg1/=0 .and. nrbg2/=0 .and. nrbg3/=0) then
                  write (*,*) 'Zebra GS, line = ',zd
                else
                  write (*,*) 'Skipping GS zebra line = ',zd
                endif
              endif

              if (out.ge.3.and.prnt) then
                if (nrbg1/=0 .and. nrbg2/=0 .and. nrbg3/=0) then
                  write (*,*) '  colors=',nrbg1,nrbg2,nrbg3
                  write (*,*) '  ilimits=',i1min,i1max
                  write (*,*) '  jlimits=',j1min,j1max
                  write (*,*) '  klimits=',k1min,k1max
                endif
              endif

              if (     (minzd==maxzd)
     .            .and.(nrbg1==0.and.nrbg2==0.and.nrbg3==0)) then
                write (*,*) 'zebra dir=',zd,' range=',minzd,maxzd
                write (*,*) 'colors=',nrbg1,nrbg2,nrbg3
                write (*,*) 'ilimits=',i1min,i1max
                write (*,*) 'jlimits=',j1min,j1max
                write (*,*) 'klimits=',k1min,k1max
                call pstop('gs','No zebra GS relaxation is possible!')
              endif

              do irbg3 = 1,nrbg3
               do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                  do zit=1,zebra_it

                    call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)

                    if (out.ge.3.and.prnt) then
                      write (*,'(a,i2,a,i2,a,i2,a,3i4)')
     .                                'Inner zebra GS: dir=',zd
     .                               ,'; it= ',zit
     .                               ,'; color=',irbg1
     .                                  + nrbg1*(irbg2-1)
     .                                  + nrbg2*(irbg3-1)
cc     .                               ,'; grid pos=',i,j,k
                    endif

                    do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                     do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                        call scheduleIndices(igrid,i1,j1,k1,i,j,k
     .                                      ,zebra_dir=zd)

cc                        write (*,*) 'DIAG -- zebra GS: it= ',zebra_it
cc     .                             ,', dir=',zd
cc     .                             ,', grid pos=',i,j,k
cc
cc                        ii = i + nxf*(j-1) + nxf*nyf*(k-1)
cc                        
cc                        call MGmatvec(-ii,neq,ntot,zz,yy,igrid,bcnd
cc     .                               ,matvec)

cc                        if (zd==2 .and. isSP2(i,igrid,dim=1)) cycle

                        !nxf,nyf,nzf are the original grid's even in line relaxation.
                        iii = neq*(i-1 + nxf*(j-1) + nxf*nyf*(k-1))

                        !Find new residual
                        do ieq = 1,neq
                          rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
                        enddo

                        !Multiply by D^-1 (stored in diag)
                        iig = iii + isig - 1
                        dummy=matmul(transpose(diag(:,iig+1:iig+neq))
     .                              ,rhs)

                        !Update zz
                        do ieq=1,neq
                          zz(iii+ieq) = zz(iii+ieq)+zebra_omg*dummy(ieq)
                        enddo

                        if (zit == zebra_it) mag = mag + sum(rhs*rhs)

                      enddo
                     enddo
                    enddo

                  enddo

                 enddo
               enddo
              enddo

            enddo

c         STANDARD RELAXATION
          else

            gs_type = 'Colored'

            call scheduleColors(igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                         ,i1min,i1max,j1min,j1max,k1min,k1max)


            if (out.ge.3.and.prnt) then
              write (*,*) '  colors=',nrbg1,nrbg2,nrbg3
              write (*,*) '  ilimits=',i1min,i1max
              write (*,*) '  jlimits=',j1min,j1max
              write (*,*) '  klimits=',k1min,k1max
            endif

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)

                do k1=k1min+mod((    irbg3-1),nrbg3),k1max,nrbg3
                  do j1=j1min+mod((k1+irbg2+irbg1-1),nrbg2),j1max,nrbg2
                    do i1=i1min+mod((j1+k1+irbg1-1),nrbg1),i1max,nrbg1

                      call scheduleIndices(igrid,i1,j1,k1,i,j,k)

cc                      if (.not.isSP(i,j,k,igrid,igrid,igrid)) then
                        !nxf,nyf,nzf are the original grid's even in line relaxation.
                        iii = neq*(i-1 + nxf*(j-1) + nxf*nyf*(k-1))

                        !Find new residual
                        do ieq = 1,neq
                          rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
                        enddo

                        !Multiply by D^-1 (stored in diag)
                        iig = iii + isig - 1
                        dummy=matmul(transpose(diag(:,iig+1:iig+neq))
     .                              ,rhs)

                        !Update zz (except around SP)
                        do ieq=1,neq
                          zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
                        enddo

                        mag = mag + sum(rhs*rhs)
cc                      endif

                    enddo
                  enddo

                enddo

                enddo
              enddo
            enddo

          endif

c       -----------
c       Standard GS
c       ----------- 
        else

c         VERTEX-BASED RELAXATION
          if (vbr) then
            gs_type = 'L-VERTX'

            !Vertex sampling
            do kv=1,max(nzf-mg_ratio_z(igrid)+1,1)
              do jv=1,max(nyf-mg_ratio_y(igrid)+1,1)
                do iv=1,max(nxf-mg_ratio_x(igrid)+1,1)

                !Find new residual
                do kblock=1,mg_ratio_z(igrid)
                  do jblock=1,mg_ratio_y(igrid)
                    do iblock=1,mg_ratio_x(igrid)

                    iib = iblock + mg_ratio_x(igrid)*(jblock-1)
     .                           + mg_ratio_x(igrid)
     .                            *mg_ratio_y(igrid)*(kblock-1)

                    if = iv-1 + iblock
                    jf = jv-1 + jblock
                    kf = kv-1 + kblock

                    ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                    call MGmatvec(-ii,neq,ntot,zz,yy,igrid,bcnd,matvec)

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
     .                           + mg_ratio_x(igrid)
     .                            *mg_ratio_y(igrid)*(kblock-1)

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

c         STANDARD RELAXATION
          else
 
           gs_type = 'Lexicog'

c           Forward pass
            call limits(0,nxf,nyf,nzf,igrid
     .                 ,imin,imax,jmin,jmax,kmin,kmax)

            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  !nxf,nyf,nzf are the original grid's even in line relaxation.
                  ii = i + nxf*(j-1) + nxf*nyf*(k-1)

                  call MGmatvec(-ii,neq,ntot,zz,yy,igrid,bcnd,matvec)

                  iii = neq*(ii-1)

                  !Find new residual
                  do ieq = 1,neq
                    rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
                  enddo

                  !Multiply by D^-1 (stored in diag)
                  iig = iii + isig - 1
                  dummy = matmul(transpose(diag(:,iig+1:iig+neq)),rhs)

                  !Update solution zz
                  do ieq=1,neq
                    zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
                  enddo

                  mag = mag + sum(rhs*rhs)

                enddo
              enddo
            enddo

ccc           Backward pass
cc            do k=kmax,kmin,-1
cc              do j=jmax,jmin,-1
cc                do i=imax,imin,-1
cc
cc                  ii = i + nxf*(j-1) + nxf*nyf*(k-1)
cc
cc                  call MGmatvec(-ii,neq,ntot,zz,yy,igrid,bcnd,matvec)
cc
cc                  iii = neq*(ii-1)
cc
cc                  !Find new residual
cc                  do ieq = 1,neq
cc                    rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
cc                  enddo
cc
cc                  !Multiply by D^-1
cc                  iig = iii + isig - 1
cc                  dummy = matmul(transpose(diag(:,iig+1:iig+neq)),rhs)
cc
cc                  !Update solution zz
cc                  do ieq=1,neq
cc                    zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
cc                  enddo
cc
cc                  mag = mag + sum(rhs*rhs)
cc
cc                enddo
cc              enddo
cc            enddo

          endif

        endif

cc#if 0
ccc       Take care of points around singular point (using JB)
cc
cc        if (isSP(1,1,1,igrid,igrid,igrid)) then
cc#if defined(petsc)
cc          if (np > 1) call pstop('gs','Does not work in parallel')
cc#endif
cc          do i1=3,1,-1
cc
cc            mg_ctx_temp%iline(igrid) = i1   !Set line smoothing around SP
cc
cc            call limits(0,nxf,nyf,nzf,igrid
cc     .                 ,imin,imax,jmin,jmax,kmin,kmax)
cc
cccc            write (*,*) imin,imax,jmin,jmax,kmin,kmax
cc
cc            do itr2=1,5    !JB iteration around SP
cc
cc              magsp = 0d0
cc
cc              call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)
cc
cc              do k=kmin,kmax
cc                do j=jmin,jmax
cc                  do i=imin,imax
cc                    !nxf,nyf,nzf are the original grid's even in line relaxation.
cc                    ii = i + nxf*(j-1) + nxf*nyf*(k-1)
cc
cc                    iii = neq*(ii-1)
cc
cc                    !Find new residual
cc                    do ieq = 1,neq
cc                      rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
cc                    enddo
cc
cc                    !Multiply by D^-1 (stored in diag)
cc                    iig = iii + isig - 1
cc                    dummy = matmul(transpose(diag(:,iig+1:iig+neq))
cc     .                            ,rhs)
cc
cc                    !Update solution zz
cc                    do ieq=1,neq
cc                      zz(iii+ieq) = zz(iii+ieq) + 0.7*dummy(ieq)
cc                    enddo
cc
cc                    magsp = magsp + sum(rhs*rhs)
cc                  enddo
cc                enddo
cc              enddo
cc
cc              if (itr2 == 1) then
cc                magsp0 = magsp
cc              else
cc                if (out.ge.3.and.prnt)
cc     .               write (*,*) '  SP smoothing in GS: line = ',i1
cc     .                   ,'iter=',itr2,'Ratio=',sqrt(magsp/magsp0)
cc                if (   sqrt(magsp/magsp0) < 0.5
cc     .             .or.sqrt(magsp/magsp0) > 1d0) exit
cc              endif
cc
cc            enddo
cc
cc            if (out.ge.3.and.prnt) write (*,*)
cc
cc            if (out.eq.2.and.prnt)
cc     .        write (*,*) '  SP smoothing in GS: line = ',i1
cc     .                   ,'iter=',itr2,'Ratio=',sqrt(magsp/magsp0)
cc
cccc            if (magsp/magsp0 > 1d0) stop
cc          enddo
cc
cc          mg_ctx_temp%iline(igrid) = 0 !Unset line smoothing
cc
cccc        !GS
cccc        call limits(0,nxf,nyf,nzf,igrid
cccc     .             ,imin,imax,jmin,jmax,kmin,kmax)
cccc
cccc        i = 1
cccc
cccc        do k=kmin,kmax
cccc         do itr2=1,10
cccc          do j=jmin,jmax
cccc
cccc            !nxf,nyf,nzf are the original grid's even in line relaxation.
cccc            ii = i + nxf*(j-1) + nxf*nyf*(k-1)
cccc
cccc            call MGmatvec(-ii,neq,ntot,zz,yy,igrid,bcnd,matvec)
cccc
cccc            iii = neq*(ii-1)
cccc
cccc            !Find new residual
cccc            do ieq = 1,neq
cccc              rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
cccc            enddo
cccc
cccc            !Multiply by D^-1 (stored in diag)
cccc            iig = iii + isig - 1
cccc            dummy = matmul(transpose(diag(:,iig+1:iig+neq)),rhs)
cccc
cccc            !Update solution zz
cccc            do ieq=1,neq
cccc              zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
cccc            enddo
cccc
cccc            mag = mag + sum(rhs*rhs)
cccc
cccc          enddo
cccc         enddo
cccc        enddo
cc
cc        endif
cc#endif

c     Check convergence

#if defined(petsc)
cc        if (.not.asm .and. out >= 1) then   !Do reduction only for debugging
        if (.not.asm) then
           lmag = mag
           call MPI_Allreduce(lmag,mag,1,MPI_DOUBLE_PRECISION
     .                       ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc        else
cc          mag = mag0                  !Do NOT do reduction in parallel (convergence is set by #its)
        endif
#endif

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          mag1 = mag
          if (mag.lt.1d-20*nn) then
            mag0 = 1d0
            mag1 = 1d0
            exit
          endif
          if (out.ge.2.and.prnt) write(*,10)gs_type,itr-1,mag,mag/mag0
        else
          if (out.ge.2.and.prnt) write(*,20)gs_type,itr-1,mag,mag/mag0
     .                                     ,mag/mag1
          mag1 = mag
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.1) then
        call MGmatvec(0,neq,ntot,zz,yy,igrid,bcnd,matvec)
        call MGvecadd(igrid,neq,ntot,-1d0,yy,1d0,rr)  !yy=rr-yy
        mag = sqrt(MGdot(igrid,neq,ntot,yy,yy))       !sqrt(yy*yy)
      endif

      if (out.ge.2.and.prnt) then
        write (*,20) gs_type,itr,mag,mag/mag0,mag/mag1
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

      call deallocPointers(fpointers)

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

      end subroutine gs

c find_mf_diag_std
c####################################################################
      subroutine find_mf_diag_std(neq,ntot,matvec,igrid,bbcnd,diag1)
c--------------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec
c     for all grids, commencing with grid "igrid".
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer    :: neq,ntot,bbcnd(6,neq),igrid

      real(8)    :: diag1(neq,2*ntot)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot),mat(neq,neq),mat2(neq,neq)
     $             ,idet
      integer    :: ii,jj,nn,ig,iig,isig
      integer    :: igr,ieq,alloc_stat
      logical    :: fpointers

c Begin program

c Allocate MG pointers

      call allocPointers(neq,fpointers)

      if (fpointers) igmax=grid_params%ngrid  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        call pstop('find_mf_diag_std'
     .            ,'Error in input of find_mf_diag_stc')
      endif

c Form diagonal

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

            call MGmatvec(ii,neq,nn,x1,dummy,igr,bbcnd,matvec)

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

c Deallocate pointers

      call deallocPointers(fpointers)

c End program

      end subroutine find_mf_diag_std

c find_mf_diag_colored
c####################################################################
      subroutine find_mf_diag_colored(neq,ntot,matvec,igrid,bbcnd,diag1
     .                               ,ncolors)
c--------------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec
c     for all grids, commencing with grid "igrid". This routine
c     employs a colored arrangement of the grid that allows one
c     to find the diagonal using full-vector matvecs, instead of
c     point-wise matvecs.
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer    :: neq,ntot,bbcnd(6,neq),igrid,ncolors

      real(8)    :: diag1(neq,2*ntot)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot),mat(neq,neq),mat2(neq,neq)
     .             ,idet,coeff
      integer    :: ii,jj,nn,ig,iig,isig,bbcnd2(6,neq)
      integer    :: igr,ieq,alloc_stat
      logical    :: fpointers

      integer    :: i,j,k,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3
      integer    :: i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max

c Begin program

c Allocate MG pointers

      call allocPointers(neq,fpointers)

      if (fpointers) igmax=grid_params%ngrid  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        call pstop('find_mf_diag_colored'
     .            ,'Error in input of find_mf_diag_colored')
      endif

c Form diagonal

      diag1 = 0d0

      if (ncolors == 1) then

        call find_mf_diag_std(neq,ntot,matvec,igrid,bbcnd,diag1)

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

          call scheduleColors(igr,ncolors,nrbg1,nrbg2,nrbg3
     .                       ,i1min,i1max,j1min,j1max,k1min,k1max)

          do ieq = 1,neq

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                  do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                    do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                        call scheduleIndices(igr,i1,j1,k1,i,j,k)

                        ii = i +nxv(igr)*(j-1) +nxv(igr)*nyv(igr)*(k-1)

                        x1((ii-1)*neq+ieq) = 1d0
                      enddo
                    enddo
                  enddo

                  call MGmatvec(0,neq,nn,x1,dummy,igr,bbcnd2,matvec)

                  do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                    do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                        call scheduleIndices(igr,i1,j1,k1,i,j,k)

                        ii = i +nxv(igr)*(j-1) +nxv(igr)*nyv(igr)*(k-1)

                        jj = (ii-1)*neq + isig - 1

                        x1((ii-1)*neq+ieq) = 0d0

                        diag1(ieq,jj+1:jj+neq)
     .                                     = dummy((ii-1)*neq+1:ii*neq)
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
cccc            call scheduleColors(igr,ncolors,nrbg1,nrbg2,nrbg3
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
cc     .               diag1(ieq,jj+1:jj+neq) = dummy((ii-1)*neq+1:ii*neq)
cc
cc                enddo
cc              enddo
cc
cc            enddo
cc
cc          endif

c       Invert diagonal and store in diag1

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

      endif

c Deallocate pointers

      call deallocPointers(fpointers)

c End program

      end subroutine find_mf_diag_colored

c find_mf_diag
c####################################################################
      subroutine find_mf_diag(neq,nblk,ntot,matvec,igrid,bbcnd
     .                       ,diag1,ncolors,volw,glrkin)
c--------------------------------------------------------------------
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
c       * diag1: returns inverse of diagonal.
c       * ncolors: number of colors for colored algorithm
c       * volw: whether residuals are volume-weighed
c       * glrkin: whether we use galerkin coarsening
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer    :: neq,nblk,ntot,bbcnd(6,neq),igrid,ncolors

      real(8)    :: diag1(neq*nblk,2*ntot*nblk)

      logical    :: volw,glrkin

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot),mag(3)
      integer    :: ii,jj,nn,ig,isig,size
      integer    :: iv,jv,kv,if,jf,kf,nxf,nyf,nzf
      integer    :: iii,icolb,irowb,icol,iiv,iig,ivg,irowmin,irowmax
     .             ,iblock ,jblock ,kblock
     .             ,iblock2,jblock2,kblock2,nblkg

      integer    :: igr,igc,ieq,alloc_stat
      logical    :: fpointers

c Begin program

      volf = volw
      galerkin = glrkin

c Allocate MG pointers

      call allocPointers(neq,fpointers)

      if (fpointers) igmax=grid_params%ngrid  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        call pstop('find_mf_diag'
     .            ,'Error in input of find_mf_diag')
      endif

c Form diagonal

c     STANDARD RELAXATION
      if (nblk == 1) then

        call find_mf_diag_colored(neq,ntot,matvec,igrid,bbcnd,diag1
     .                           ,ncolors)

C     VERTEX-BASED RELAXATION
      else

        diag1 = 0d0

        do igr = igrid,igmax

          nxf = nxv(igr)
          nyf = nyv(igr)
          nzf = nzv(igr)

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
                            call MGmatvec(jj,neq,nn,x1,dummy,igr,bbcnd
     .                                   ,matvec)

                            !Reset vector x1 to zero
                            x1(iii) = 0d0

                            !Store results in diagonal
                            diag1(icol,irowmin:irowmax)
     .                                = dummy((jj-1)*neq+1:jj*neq)
                          enddo

                        enddo
                      enddo
                    enddo

                  enddo
                enddo
              enddo

c             Invert diagonal and store in diag1

              call blockInv(size,diag1(1:size,ivg+1:ivg+size))

              enddo
            enddo
          enddo

        enddo

      endif

c Deallocate pointers

      call deallocPointers(fpointers)

c End program

      end subroutine find_mf_diag

c find_mf_mat
c####################################################################
      subroutine find_mf_mat(neq,ntot,matvec,igrid,bbcnd,mat)
c--------------------------------------------------------------------
c     Finds all matrix elements matrix-free using subroutine matvec
c     for all grids, commencing with grid "igrid".
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer    :: neq,ntot,bbcnd(6,neq),igrid

      real(8)    :: mat(ntot,ntot)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot)
      integer    :: ii,jj,nn
      integer    :: igr,ieq,alloc_stat
      logical    :: fpointers

c Begin program

      igr = igrid

c Allocate MG pointers

      call allocPointers(neq,fpointers)

c Consistency check

      nn  = ntotv(igr)

      if (nn /= ntot) then
        call pstop('find_mf_mat'
     .            ,'Error in input of find_mf_mat')
      endif

c Form diagonal 

      x1 = 0d0

      do ii = 1,nn

c       Find column vector corresponding to grid node ii and equation ieq

          x1(ii) = 1d0

          call MGmatvec(0,neq,nn,x1,dummy,igr,bbcnd,matvec)

          x1(ii) = 0d0

c       Fill matrix

          mat(:,ii) = dummy

      enddo

c Deallocate pointers

      call deallocPointers(fpointers)

c End program

      end subroutine find_mf_mat

c     symm_test
c     ###############################################################
      subroutine symm_test(neq,igrid,matvec,bcnd)
c     ---------------------------------------------------------------
c     Performs symmetry test of matvec on grid "igrid".
c     ---------------------------------------------------------------

      use mg_internal

      implicit none     !For safe fortran

c     Call variables

      integer    :: neq,igrid,bcnd(6,neq)

      external      matvec

c     Local variables

      real(8),allocatable,dimension(:)::x1,dummy,dummy2

      real(8)    :: dd1,dd2,err
      integer    :: nx,ny,nz,nn,ii,jj,i1,j1,i2,j2,ix1,iy1,ix2,iy2,ieq

c     Begin program

c     Initialize variables

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

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

          call MGmatvec(0,neq,nn,x1,dummy,igrid,bcnd,matvec)

          call findBaseVector(ii,ieq,neq,nn,x1,0d0)

c       Compare column vector ii with corresponding row vector (intersect in
c       diagonal)

          do jj = ii,nn/neq

            call findBaseVector(jj,ieq,neq,nn,x1,1d0)

c THIS DOESN'T WORK ANYMORE. NEED TO FIX.
cc            call matvec(ii,neq,nn,x1,dummy2,igrid,bcnd)
            call MGmatvec(0,neq,nn,x1,dummy2,igrid,bcnd,matvec)

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

      stop

c     End program

      deallocate (x1,dummy,dummy2)

 15   format ('(',i3,',',i3,'):',1pe10.2,'; (',i3,',',i3,'):',e10.2,
     .        '  Error:',e10.2,'  %error:',0pf7.2)
 20   format (/,'Total error:',1pe10.3)

      end subroutine symm_test

c findBaseVector
c####################################################################
      subroutine findBaseVector(ii,ieq,neq,ntot,x1,coef)
c--------------------------------------------------------------------
c     Finds base vector corresponding to grid node ii and equation ieq.
c     Assumes that x1 has been initialized to zero elsewhere.
c--------------------------------------------------------------------

      implicit none      !For safe fortran

c Call variables

      integer    :: neq,ii,ieq,ntot

      real(8)    :: x1(ntot),coef

c Local variables

c Begin program

      x1(neq*(ii-1) + ieq) = coef

c End program

      end subroutine findBaseVector

c MGmatvec
c####################################################################
      subroutine MGmatvec(gposc,neq,ntotc,x,y,igc,bcnd,matvec)
c--------------------------------------------------------------------
c     This subroutine calculates the matrix vector product defined
c     symbolically in matvec, in two ways:
c       * Rediscretizing (i.e., using matvec straight)
c       * Galerkin (preprocessing w/ prolongation, and postprocessing
c                   w/ restriction).
c
c     In call:
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
c--------------------------------------------------------------------

      use mg_internal

      use mgarraySetup

      implicit none

c Call variables

      integer    :: neq,ntotc,igc,gposc,bcnd(6,neq)
      real(8)    :: x(ntotc),y(ntotc),minv,idx,idy,idz,
     .              idxcnv,idycnv,idzcnv,vxx,vyy,vzz

      external      matvec

c Local variables

      integer    :: nxf,nyf,nzf,ntotf,igf,gposf,if,jf,kf
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf
      integer    :: nxc,nyc,nzc,isigc,ic,jc,kc
     .             ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc
      integer    :: mgrx,mgry,mgrz,ig
      real(8),allocatable,dimension(:) :: MGx,MGy

      logical    :: fpointers

c Begin program

      call allocPointers(neq,fpointers)

      if (igc > 1 .and. galerkin) then

        isigc  = istart(igc)

cc        nxc = mg_ctx_temp%nxv(igc)
cc        nyc = mg_ctx_temp%nyv(igc)
cc        nzc = mg_ctx_temp%nzv(igc)
        nxc = grid_params%nxv(igc)
        nyc = grid_params%nyv(igc)
        nzc = grid_params%nzv(igc)

c     Allocate MG vector

        igf = 1

cc        nxf = mg_ctx_temp%nxv(igf)
cc        nyf = mg_ctx_temp%nyv(igf)
cc        nzf = mg_ctx_temp%nzv(igf)
        nxf = grid_params%nxv(igf)
        nyf = grid_params%nyv(igf)
        nzf = grid_params%nzv(igf)

        ntotf = neq*nxf*nyf*nzf

        allocate(MGx(2*ntotf),MGy(2*ntotf))

c     Prolong vector to finest grid

        MGx(isigc:isigc+ntotc-1) = x

        call prolongMGVector(neq,MGx,igc,igf,0,bcnd)

c     Perform matvec in finest grid

        if (gposc == 0) then

          call matvec(0,neq,ntotf,MGx(1:ntotf),MGy(1:ntotf),igf,bcnd)

        else

          !Find coarse cell position (ic,jc,jc)
          call limits(gposc,nxc,nyc,nzc,igc
     .               ,ic,imaxc,jc,jmaxc,kc,kmaxc)

          !Define finest grid stencil
          mgrx = product(mg_ratio_x(igf:igc-1))
          mgry = product(mg_ratio_y(igf:igc-1))
          mgrz = product(mg_ratio_z(igf:igc-1))

          iminf = mgrx*(ic-1)+1
          imaxf = mgrx* ic

          jminf = mgry*(jc-1)+1
          jmaxf = mgry* jc

          kminf = mgrz*(kc-1)+1
          kmaxf = mgrz* kc

cc          write (*,*) 'DIAG -- MGmatvec',igf,igc,mgrx,mgry,mgrz
cc          write (*,*) 'DIAG -- MGmatvec',ic,jc,kc
cc          write (*,*) 'DIAG -- MGmatvec',iminf,imaxf,jminf,jmaxf
cc     .                                  ,kminf,kmaxf
cc          write (*,*)

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

        call restrictMGVector(neq,MGy,igf,igc,0,volf)

c     Map coarse grid residual

        y = MGy(isigc:isigc+ntotc-1)

c     Deallocate MG vectors

        deallocate(MGx,MGy)

      else

        call matvec(gposc,neq,ntotc,x,y,igc,bcnd)

      endif

      call deallocPointers(fpointers)

c     End program

      end subroutine MGmatvec

c     module mg_solver
c     ##################################################################
      module mg_solver

        use mg_internal

      contains

c     cSolver
c     ##################################################################
      subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess,matvec
     .                  ,vol_wgt
     .                  ,dg,cvrg_tst,tol,iters
     .                  ,gm_driver
     .                  ,gm_sing_matrix
     .                  ,cg_driver
     .                  ,ks_it
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
c       * gm_sing_mat: whether GMRES is dealing with a singular matrix or not
c       * cg_driver: whether we want CG    as a driver
c       * ks_it: # of krylov subspace iterations
c
c       * mg_order_res: restriction order
c       * mg_order_prol: restriction order
c       * mg_smooth: determines smoothing (jb,gs,gm; default jb)
c       * mg_vcyc: maximum number of MG V-cycles (default 1)
c       * mg_cycle: number of mu_cycles in MG (default 1)
c       * mg_galerkin: whether to do Galerkin coarsening of operator
c                  (true) or rediscretization (default: false)
c       * mg_coarse_grid_size: minimum coarse mesh (roughly mg_ratio^mg_cgs per
c                 direction).
c       * mg_debug: turn on debugging in MG
c
c       * sm_it: smoother iterations
c       * sm_omega: damping for JB smoother
c       * sm_ncolors: number of colors in grid (for GS).
c       * sm_line_relax: whether we want line relaxation
c       * sm_zebra_relax: whether we want zebra GS relaxation
c       * sm_zebra_omega: inner relaxation factor in ZEBRA smoothing (between directions)
c       * sm_zebra_it: inner ZEBRA iterations (per direction).
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
     .                                ,sm_ncolors
     .                                ,ks_it
        logical,optional,intent(IN) :: sm_line_relax
     .                                ,sm_zebra_relax
     .                                ,mg_gm_coarse_solve
     .                                ,mg_galerkin
     .                                ,mg_debug
     .                                ,gm_driver
     .                                ,gm_sing_matrix
     .                                ,cg_driver
     .                                ,cvrg_tst
        integer,optional,intent(OUT) :: iters

        character(2),optional :: mg_smooth

        external   matvec

c     Local variables

        integer    :: ntot,iorder_res,iorder_prol,vcyc,mu_cyc,nclrs,smit
     .               ,mgmin,i,ieq,gmit,zit
        real(8)    :: xi(ntotp*neq),bi(ntotp*neq),omg,rtol,zomg
        real(8), target :: ddg(neq,2*neq*ntotp)
        logical    :: gm_crse_solve,cvrg_test,lrelax,zrelax
     .               ,gm_drvr,gm_smat,cg_drvr,glkin

        character(2) :: smth

c     Begin program

c     Process optional arguments

        if (PRESENT(dg)) ddg = dg

        !Solver options <<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if (PRESENT(gm_driver)) then
          gm_drvr = gm_driver
        else
          gm_drvr = .false.
        endif

        if (PRESENT(gm_sing_matrix)) then
          gm_smat = gm_sing_matrix
        else
          gm_smat = .false.
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
          zit = 1
        endif

        if(PRESENT(sm_it)) then
          smit = sm_it
        else
          smit = 3
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
cc            solverOptions%singular_matrix = gm_smat
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
            solverOptions%singular_matrix = gm_smat
            call assembleSolverHierarchy('gm')
          else
            call assembleSolverHierarchy('cg')
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
        solverOptions%vol_res  = vol_wgt
        solverOptions%mg_galerkin = glkin

        if (PRESENT(mg_debug)) solverOptions%mg_debug = mg_debug

        if (PRESENT(dg)) solverOptions%diag => ddg

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

        solverOptions%iter    = smit
        solverOptions%tol     = rtol

        if (lrelax) then
          solverOptions%iter  = 2
        endif

        select case (smth)
        case('gs')
          solverOptions%omega   = omg
          solverOptions%ncolors = nclrs
          solverOptions%mg_zebra_relax = zrelax
          solverOptions%mg_zebra_prefd = 0
cc          if(bcSP()) solverOptions%mg_zebra_prefd = 1
          solverOptions%mg_zebra_it    = zit
          solverOptions%mg_zebra_omega = zomg

cc          solverOptions%iter=ceiling(1.*smit/solverOptions%mg_zebra_it)
          solverOptions%iter=smit

          call assembleSolverHierarchy('gs')
        case('gm')
          solverOptions%stp_test        = 1 
          solverOptions%krylov_subspace = solverOptions%iter 
          solverOptions%singular_matrix = gm_smat

          call assembleSolverHierarchy('gm')

          call assembleSolverHierarchy('id')
        case default !JB
          solverOptions%omega   = omg
          call assembleSolverHierarchy('jb')
        end select

c       MG coarsest grid solve

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = 100
        solverOptions%iter            = 100
        solverOptions%stp_test        = 1 
        solverOptions%omega           = 1d0
        solverOptions%singular_matrix = gm_smat

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

      end module mg_solver

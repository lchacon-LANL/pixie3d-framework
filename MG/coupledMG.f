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

c module setMGBC_interface
c######################################################################
      module setMGBC_interface

        INTERFACE
          subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd,arr_cov
     .                      ,arr0,icomp,is_cnv,is_vec,result_is_vec
     .                      ,iorder)
            integer(4) :: nnx,nny,nnz,neq,bcnd(6,neq),iig,gpos
            real(8)    :: array(0:nnx+1,0:nny+1,0:nnz+1,neq)
            real(8),optional,intent(INOUT) ::
     .                    arr_cov(0:nnx+1,0:nny+1,0:nnz+1,neq)
            real(8),optional,intent(IN) ::
     .                    arr0   (0:nnx+1,0:nny+1,0:nnz+1,neq)
            integer(4),optional,intent(IN) :: icomp,iorder
            logical   ,optional,intent(IN) :: is_cnv,is_vec
     .                                       ,result_is_vec
          end subroutine setMGBC
        END INTERFACE

      end module setMGBC_interface

c module mg_internal
c######################################################################
      module mg_internal

        use grid

        use mlsolverSetup

        integer(4) :: ngrdx,ngrdy,ngrdz,ngrid,igmax

        real(8), allocatable, dimension(:,:) :: diag

        integer(4),dimension(:),allocatable ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z

        logical :: vbr_mg

        type (grid_def) :: MGgrid,MGgrid_sv

      contains

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,fpointers,mg_grid)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in), integer: number of coupled equations to be solved
c       * mg_grid (in/out), type(grid_def): contains definition of MG grid
c            level structure.
c       * fpointers (out), logical: indicates if pointers were allocated
c            or not.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer(4)          :: neq
      logical,intent(OUT) :: fpointers
      type(grid_def),optional,intent(IN) :: mg_grid

c     Local variables

      integer(4) :: i

c     Begin program

      fpointers = .false.

      if (PRESENT(mg_grid)) then
        MGgrid=mg_grid
      else
        if (.not.associated(MGgrid%xx)) MGgrid=grid_params     !Failsafe grid definition
      endif

      ngrdx = MGgrid%ngrdx
      ngrdy = MGgrid%ngrdy
      ngrdz = MGgrid%ngrdz
      ngrid = MGgrid%ngrid

c     Check if MG pointers are allocated

      if (.not.allocated(istart)) then

        allocate(istart (ngrid),ntotv (ngrid)
     .          ,istartp(ngrid),ntotvp(ngrid)
     .          ,istartb(ngrid),ntotb (ngrid)
     .          ,nxv(ngrid),nyv(ngrid),nzv(ngrid),nblock(ngrid)
     .          ,mg_ratio_x(ngrid),mg_ratio_y(ngrid),mg_ratio_z(ngrid))

        mg_ratio_x = MGgrid%mg_ratio_x
        mg_ratio_y = MGgrid%mg_ratio_y
        mg_ratio_z = MGgrid%mg_ratio_z

        nxv = MGgrid%nxv
        nyv = MGgrid%nyv
        nzv = MGgrid%nzv

        istartp(1) = 1
        istart (1) = 1
        istartb(1) = 1
        nblock (1) = mg_ratio_x(1)*mg_ratio_y(1)*mg_ratio_z(1)
        ntotvp (1) = nxv(1)*nyv(1)*nzv(1)
        ntotv  (1) = neq*ntotvp(1)
        ntotb  (1) = nblock(1)*ntotv(1)
        do i = 2,ngrid
          istartp(i) = istartp(i-1) + ntotvp(i-1)
          istart (i) = istart (i-1) + ntotv (i-1)
          istartb(i) = istartb(i-1) + ntotb (i-1)
          nblock (i) = mg_ratio_x(i)*mg_ratio_y(i)*mg_ratio_z(i)
          ntotvp (i) = nxv(i)*nyv(i)*nzv(i)
          ntotv  (i) = neq*ntotvp(i)
          ntotb  (i) = nblock(i)*ntotv(i)
        enddo

        fpointers = .true.

      endif

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
        deallocate(istart,ntotv,istartp,ntotvp
     .            ,istartb,ntotb,nxv,nyv,nzv,nblock
     .            ,mg_ratio_x,mg_ratio_y,mg_ratio_z)
      endif

      end subroutine deallocPointers

c     getMGvcomp
c     #################################################################
      function getMGvcomp(i,j,k,nx,ny,nz,igr,ieq,neq) result(ijkg)

c     -----------------------------------------------------------------
c     Gives MG vector component corresponing to equation ieq on
c     coordinates (i,j,k) and grid level igr.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igr,ieq,neq,ijkg,nx,ny,nz

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

      integer(4) :: neq,igr,nx,ny,nz
      real(8)    :: mgvector(*),array(0:nx+1,0:ny+1,0:nz+1,neq)

c     Local variables

      integer(4) :: i,j,k,ii,ieq
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

      integer(4) :: neq,igr,nx,ny,nz,gpos
      real(8)    :: mgvector(*),array(0:nx+1,0:ny+1,0:nz+1,neq)
      logical    :: ismgvec

c     Local variables

      integer(4) :: ieq,i,j,k,ii,offset,igptr
      integer(4) :: imin ,imax ,jmin ,jmax ,kmin ,kmax
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

      integer(4) :: elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax,igr

c     Local variables

      integer(4) :: el1

c     Begin program

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

        if ( MGgrid%iline(igr) /= 0 ) then
          imin = MGgrid%iline(igr)
          imax = imin
        endif
        if ( MGgrid%jline(igr) /= 0 ) then
          jmin = MGgrid%jline(igr)
          jmax = jmin
        endif
        if ( MGgrid%kline(igr) /= 0 ) then
          kmin = MGgrid%kline(igr)
          kmax = kmin
        endif
      endif

c     End program

      end subroutine limits

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

        integer(4) :: size,icol

        real(8)    :: mat(size,size),rhs(size,icol),x(size,icol)

c     Local variables

        integer(4) :: ipiv(size),info

        real(8)    :: mat2(size,size)

        external dgesv

c     Begin program

        mat2 = mat !Avoid overwritting mat
        x = rhs

        call dgesv(size,icol,mat2,size,ipiv,x,size,info) !LAPACK routine

        if (info /= 0) then
          if (info < 0) then
            write (*,*) 'Problem in factorization in argument',-info
            write (*,*) 'Aborting'
            stop
          else
            write (*,*) 'Matrix is singular'
            write (*,*) 'Aborting'
            stop
          endif
        endif

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

        integer(4) :: size

        real(8)    :: mat(size,size),matinv(size,size)

c     Local variables

        integer(4) :: ipiv(size),info,lwork

        real(8),allocatable,dimension(:) :: work

        external   dgetri,dgetrf

c     Begin program

c     Find LU decomposition

        call dgetrf(size,size,mat,size,ipiv,info) !LAPACK routine

        call error(info)

c     Invert matrix

        !Workspace query
        lwork=-1
        allocate(work(1))
        call dgetri(size,mat,size,ipiv,work,lwork,info) !LAPACK routine
        lwork=floor(work(1))
        deallocate(work)

        !Matrix inversion
        allocate(work(lwork))
        call dgetri(size,mat,size,ipiv,work,lwork,info) !LAPACK routine
        deallocate(work)

        call error(info)

      contains

      subroutine error(info)

        integer(4) :: info

        if (info /= 0) then
          if (info < 0) then
            write (*,*) 'Problem in factorization in argument',-info
            write (*,*) 'Aborting'
            stop
          else
            write (*,*) 'Matrix is singular'
            write (*,*) 'Aborting'
            stop
          endif
        endif

      end subroutine error

      end subroutine blockInv

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

        implicit none

c     Call variables

        integer(4) :: iflag,igr,neq
        real(8)    :: mgv(*)
        character*(*) :: ofile

c     Local variables

        integer(4) :: nx,ny,nz,ieq,nunit
        real(8),allocatable,dimension(:,:,:,:) :: debug
        logical    :: fpointers

c     Begin program

        call allocPointers(neq,fpointers)

        nunit = 110

        nx = nxv(igr)
        ny = nyv(igr)
        nz = nzv(igr)

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
          call contour(debug(1:nx,1:ny,1,ieq),nx,ny,xmin,xmax,ymin,ymax
     .                ,iflag+ieq-1,nunit)
        enddo

cc        if (igr == igmax) then
cc          write (*,*) debug(nx/2,1:ny,1,1)
cc          write (*,*) debug(nx/2,1:ny,1,2)
cc          write (*,*) debug(nx/2,1:ny,1,3)
cc        endif

c     End program

        deallocate(debug)
        call deallocPointers(fpointers)
        close(nunit)

      end subroutine MGplot

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c     ---------------------------------------------------------------------
c     Contours 2D array in xdraw format. In call:
c       * arr: 2D array to be plotted
c       * nx,ny: dimensions of array
c       * xmin,xmax,ymin,ymax: 2D domain limits
c       * iopt: whether to initialize xdraw plot (iopt=0) or not.
c       * nunit: integer file identifier.
c     ---------------------------------------------------------------------

c     Call variables

      integer(4) :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer(4) :: i,j

c     Begin program

      if(iopt == 0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4)
     .              ,real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*)

c     End program

      end subroutine contour

c     vecadd
c     ###################################################################
      subroutine vecadd(igr,neq,ntot,coef1,vec1,coef2,vec2)

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

      integer(4) :: igr,neq,ntot
      real(8)    :: coef1,coef2,vec1(ntot),vec2(ntot)

c     Local variables

      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq

c     Begin program

      call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
     .               ,imin,imax,jmin,jmax,kmin,kmax)

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            do ieq=1,neq
              iii=getMGvcomp(i,j,k,nxv(igr),nyv(igr),nzv(igr),1,ieq,neq)
cc              iii = neq*(i-1 + nxv(igr)*(j-1)
cc     .                       + nxv(igr)*nyv(igr)*(k-1)) + ieq
              vec1(iii) = coef1*vec1(iii) + coef2*vec2(iii)
            enddo
          enddo
        enddo
      enddo

      end subroutine vecadd

c     dot
c     ###################################################################
      function dot(igr,neq,ntot,vec1,vec2)

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

      integer(4) :: igr,neq,ntot
      real(8)    :: vec1(ntot),vec2(ntot),dot

c     Local variables

      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq

c     Begin program

      call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
     .               ,imin,imax,jmin,jmax,kmin,kmax)

      dot = 0d0

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            do ieq=1,neq
              iii=getMGvcomp(i,j,k,nxv(igr),nyv(igr),nzv(igr),1,ieq,neq)
              dot = dot + vec1(iii)*vec2(iii)
            enddo
          enddo
        enddo
      enddo

      end function dot

ccc     lmtvc
ccc     ###############################################################
cc      subroutine lmtvc(gpos,neq,ntot,xl,yl,igr,bcnd)
ccc     ---------------------------------------------------------------
ccc     This subroutine is a proxy matvec routine to perform local
ccc     matvec on global grid for GMRES routine
ccc     ---------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: neq,ntot,igr,gpos,bcnd(6,neq)
cc        real(8)    :: xl(ntot),yl(ntot)
cc
ccc     Local variables
cc
cc        integer(4) :: nn,nxl,nyl,nzl
cc        integer(4) :: imin,imax,jmin,jmax,kmin,kmax
cc        integer(4) :: i,j,k,il,jl,kl,iii,iil,gloc,ieq
cc
cc        real(8),allocatable,dimension(:) :: xg,yg
cc
cc        external v_mtvc
cc
ccc     Begin program
cc
ccc     Global parameters
cc
cc        nn = ntotv(igr)
cc
cc        allocate(xg(nn),yg(nn))
cc
cc        xg = 0d0
cc        yg = 0d0
cc
ccc     Local parameters
cc
cc        nxl = MGgrid%nxv(igr)
cc        nyl = MGgrid%nyv(igr)
cc        nzl = MGgrid%nzv(igr)
cc
ccc     Map local vector to global vector
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
cc                xg(iii) = xl(iil)
cc              enddo
cc            enddo
cc          enddo
cc        enddo
cc
ccc     Perform mat-vec and map directly to local grid
cc
cc        do k=kmin,kmax
cc          do j=jmin,jmax
cc            do i=imin,imax
cc              do ieq=1,neq
cc                !Global vector and node locations
cc                gloc=getMGvcomp(i,j,k,nxv(igr),nyv(igr),nzv(igr),1,1,1)
cc                iii =getMGvcomp(i,j,k,nxv(igr),nyv(igr),nzv(igr),1
cc     .                         ,ieq,neq)
cc                !Local indices
cc                il = i - imin + 1
cc                jl = j - jmin + 1
cc                kl = k - kmin + 1
cc                iil=getMGvcomp(il,jl,kl,nxl,nyl,nzl,1,ieq,neq)
cc
cc                !Perform stencil-wise matvec
cccc                call matvec(gloc,neq,nn,xg,yg,igr,bcnd)
ccc THIS IS HARDWIRED FOR NOW!!!!!!!!!!
cc                call v_mtvc(gloc,neq,nn,xg,yg,igr,bcnd)
cc
cc                yl(iil) = yg(iii)
cc              enddo
cc            enddo
cc          enddo
cc        enddo
cc
ccc     End program
cc
cc        deallocate(xg)
cc
cc      end subroutine lmtvc 

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

        logical :: is__cnv,is__vec,have_equl

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

        integer(4)      :: neq
        type(mg_array)  :: mgarray

c     Local variables

        integer(4)      :: igrid,nxp,nyp,nzp

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
     .                          ,iscnv,isvec,equil)
c     -----------------------------------------------------------------
c     Restricts MG array in all grids with ghost nodes.
c     -----------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4)     :: neq,icmp,bcnd(6,neq),order,igrid
      type(mg_array) :: mgarray
      logical,optional,intent(IN) :: iscnv,isvec
      type(mg_array),optional,intent(IN) :: equil

c     Local variables

      integer(4)     :: igf,nxf,nyf,nzf,igc,nxc,nyc,nzc

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

      if (PRESENT(equil)) then
        have_equl = .true.
        equl = equil
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

      implicit none    !For safe fortran

c     Call variables

      integer(4) :: neq,igf,nxf,nyf,nzf,igc,nxc,nyc,nzc
     .             ,order,bcnd(6,neq),icmp
      real(8)    :: arrayf(0:nxf+1,0:nyf+1,0:nzf+1,neq)
     .             ,arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)
      logical    :: volf

c     Local variables

      integer(4) :: igridf,igridc,isigf,isigc,i,j,k,ii,ieq
     .             ,nxxf,nyyf,nzzf,nxxc,nyyc,nzzc,ntotc,ntotf
     .             ,bcmod(6,neq)

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
     .                ,order,igridf,volf)

      enddo

c     Map vector to array

      call mapMGVectorToArray(0,neq,vecc,nxc,nyc,nzc,arrayc,igc,.false.)

      if (icmp /= 0) then
        if (have_equl) then
          call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
     .                ,arr0=equl%grid(igc)%array
     .                ,icomp=icmp,is_cnv=is__cnv,is_vec=is__vec
     .                ,iorder=order)
        else
          bcmod = bcnd
          where (bcnd == EQU)
            bcmod = EXT
          end where
          call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcmod,icomp=icmp
     .              ,is_cnv=is__cnv,is_vec=is__vec,iorder=order)
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

      integer(4) :: neq,nx,ny,nz,order
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1,neq)
      real(8)    :: mgvector(*)
      logical    :: volf

c     Local variables

      integer(4) :: ieq,nxf,nyf,nzf,nxc,nyc,nzc,igridc,igridf,igr0
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

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: iter,igridmin,vcyc,crsedpth,isig,ncolors
      integer(4) :: orderres,orderprol,alloc_stat
     .             ,line_nsweep,line_crse_depth,line_vcyc

      integer(4) :: guess2,outc,mu
      integer(4) :: izero,i,j,k,ii,ivcyc,nblk

      real(8)    :: xx(2*ntot),yy(2*ntot),wrk(2*ntot),rr(ntot)
      real(8)    :: rr0,rr1,mag,mag1,mgtol,line_tol,line_omega

      logical    :: fdiag,line_relax,fpointers,volf
     .             ,line_x,line_y,line_z

      character(2) :: line_solve

c Begin program

c Consistency check

      call optionsConsistencyCheck

c Assign options

      igridmin        = options%igridmin
      vcyc            = options%vcyc
      mgtol           = options%tol
      orderres        = options%orderres
      orderprol       = options%orderprol
      fdiag           = options%fdiag
      volf            = options%vol_res
      crsedpth        = options%mg_coarse_solver_depth
      mu              = options%mg_mu
      vbr_mg          = options%vertex_based_relax
      MGgrid          = options%mg_grid_def

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

c Consistency check

      if (crsedpth == 0 .and. vbr_mg) then
        write (*,*) 'Coarsest level solver not defined'
        write (*,*) 'Cannot do vertex-based relaxation in MG'
        write (*,*) 'Aborting...'
        stop
      endif

c Set pointers and find ngrid

      call allocPointers(neq,fpointers)

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

c Check limits

      igmax = ngrid - igridmin + 2  !Determines minimum resolution

      if (igmax.lt.igrid) then
        write (*,*) ' Grid chosen is below minimum grid resolution'
        write (*,*) ' Aborting...'
        stop
      endif

      if (out.ge.2.and.igmax.gt.igrid) write (*,5)

      outc = out
      if (igmax.gt.igrid) outc = out - 2

c Find diagonal for smoothers

      if (vbr_mg) then
        nblk = nblock(igrid)
      else
        nblk = 1
      endif

      if (fdiag) then

        if (allocated(diag)) then !Diagonal is known
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          allocate(diag(neq*nblk,2*ntot*nblk))
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          allocate(diag(neq*nblk,2*ntot*nblk))
          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
     .                     ,ncolors)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

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
        rr0 = sqrt(dot(igrid,neq,ntot,y,y))
      else
        call matvec(0,neq,ntot,x,rr,igrid,bcnd)
        call vecadd(igrid,neq,ntot,-1d0,rr,1d0,y)
        rr0 = sqrt(dot(igrid,neq,ntot,rr,rr))
      endif

      if (rr0.lt.1d-16*ntot) then
        if (out.ge.1) then
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

        if (outc.ge.1.and.igmax.gt.igrid) then
          if (mu == 1) write (*,7) ivcyc
          if (mu == 2) write (*,8) ivcyc
        endif

c     Perform mu-cycle recursively

        if (igrid.eq.igmax) then
          call smooth (igrid)
          exit
        else
          call mcycle(igrid)
        endif

c     Check MG convergence

        call matvec(0,neq,ntot,xx(istart(igrid)),rr,igrid,bcnd)

        call vecadd(igrid,neq,ntot,-1d0,rr,1d0,y)

        mag = sqrt(dot(igrid,neq,ntot,rr,rr))

        mag1 = mag/rr1

        if (out.ge.3) then
          write (*,10) mag,mag1
        elseif (out.ge.2) then
          if (mu == 1) write (*,20) mag,mag1,ivcyc
          if (mu == 2) write (*,21) mag,mag1,ivcyc
        endif

        rr1 = mag

        if (mag/rr0 < mgtol .or. mag < 1d-20*ntot) exit

      enddo

c Map solution from local vector xx to external vector x

      x(:) = xx(istart(igrid):istart(igrid+1)-1)

c MG convergence info

      mag1 = mag/rr0

      if (igmax.gt.igrid) then
        if (out.eq.1) then
          if (mu == 1) write (*,20) mag,mag1,min(ivcyc,vcyc)
          if (mu == 2) write (*,21) mag,mag1,min(ivcyc,vcyc)
        elseif (out.ge.2.and.vcyc.gt.1) then
          write (*,*) 
          write (*,*) 'Final MG convergence info:'
          if (mu == 1) write (*,20) mag,mag1,min(ivcyc,vcyc)
          if (mu == 2) write (*,21) mag,mag1,min(ivcyc,vcyc)
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
 20   format (  ' MG residual:',1p,1e12.4,'; Ratio:',1e12.4,
     .          '; V-cycle #:',i3)
 21   format (  ' MG residual:',1p,1e12.4,'; Ratio:',1e12.4,
     .          '; W-cycle #:',i3)

      contains

c     optionsConsistencyCheck
c     ###################################################################
      subroutine optionsConsistencyCheck

c     -------------------------------------------------------------------
c       Checks consistency of MG input options.
c     -------------------------------------------------------------------

        if (options%vertex_based_relax .and. options%mg_line_relax) then
          write (*,*) 'Invalid setting for MG relaxation: cannot do'
          write (*,*) 'vertex-based and line-based relaxation',
     .                ' simultaneously!'
          write (*,*) 'Aborting...'
          stop
        endif

        if (      options%mg_coarse_solver_depth == 0
     .      .and. options%vertex_based_relax          ) then
          write (*,*) 'Coarsest level solver not defined'
          write (*,*) 'Cannot do vertex-based relaxation in MG'
          write (*,*) 'Aborting...'
          stop
        endif
          
      end subroutine optionsConsistencyCheck

c     mcycle
c     ###################################################################
      recursive subroutine mcycle(igr)

c     -------------------------------------------------------------------
c       Performs MG M-cycle (V-cycle, m=1, W-cycle, m=2) recursively
c       starting at grid igr.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igr

c     Local variables

        integer(4) :: igc,isigc,isig,nn,nnc,imu

c     Begin program

        igc   = igr+1

        nn    = ntotv(igr)
        nnc   = ntotv(igc)

        isig  = istart(igr)
        isigc = istart(igc)

c     Relax error/solution on grid number igr/igrid (find new xx)

c diag ****
cccc        if (igr == igmax-1) then
cccc          write (*,*) 'Plotting here at level',igr
cccc          call matvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd)
cccc          call vecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))
cccc          call MGplot(neq,wrk,igr,0,'fine2.bin')
cccc        endif
c diag ****

        call smooth(igr)

c diag ****
cc        if (igr == igmax-1) then
cc          call MGplot(neq,xx ,igr,0,'fine.bin')
cc        endif
c diag ****

c     Evaluate residual (ie wrk = yy - A xx = yy - wrk )

        wrk(isig:isig+nn-1) = 0d0

        call matvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd)

        call vecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))

c diag ****
cc        if (igr == igmax-1) then
cc          write (*,*) 'Plotting here at level',igr
cc          call MGplot(neq,wrk,igr,0,'fine2.bin')
cc        endif
c diag ****

c     Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) to a coarser grid

        call crestrict(neq
     .                ,yy(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .                ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .                ,orderres,igr,volf)

c diag ****
cc        if (igc == igmax) then
cc          write (*,*) 'Plotting here at level',igc
cc          call MGplot(neq,yy,igc,0,'coarse.bin')
cccc          write (*,*) yy(isigc:isigc+ntotv(igc)-1)
cc        endif
c diag ****

c     Initialize solution on coarse grid

        xx(isigc:isigc+nnc-1) = 0d0

c     If on coarsest grid, solve for error, else descend a grid level

        if (igc.eq.igmax) then
          call coarseSolve(igc)
        else
          do imu = 1,mu
            call mcycle(igc)
          enddo
        endif

c diag ****
cc        if (igc == igmax) then
cc          write (*,*) 'Plotting here at level',igc
cc          call MGplot(neq,xx,igc,1,'coarse.bin')
cc
cc          call matvec(0,neq,nnc,xx(isigc),wrk(isigc),igc,bcnd)
cc          call vecadd(igc,neq,nnc,-1d0,wrk(isigc),1d0,yy(isigc))
cc          call MGplot(neq,wrk,igc,1,'coarse.bin')
cc        endif
c diag ****

c     Cycle back up to grid igr updating errors (xx)
c     with fixed R.H.S. (yy)

c     Prolong error (wrk = P * xx_1) to a finer grid

        call cprolong(neq
     .               ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .               ,xx(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .               ,orderprol,igc,bcnd)

c     Update existing solution on grid igr (i.e. xx_igr): xx_igr = xx_igr + wrk 

c diag ****
cc        xx(isig:isig+nn-1) = 0d0
cc        wrk(isig:isig+nn-1) = 0d0
c diag ****

        call vecadd(igr,neq,nn,1d0,xx(isig),1d0,wrk(isig))

c diag ****
cc        if (igr == igmax-1) then
cc          write (*,*) 'Plotting here at level',igr
cc          call MGplot(neq,wrk,igr,1,'fine.bin')
cc          call MGplot(neq,xx ,igr,1,'fine.bin')
cc
cc          call matvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd)
cc          call vecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))
cc          call MGplot(neq,wrk,igr,1,'fine2.bin')
cc        endif
c diag ****

c     Relax updated solution on igr (i.e. xx_igr)

        call smooth(igr)

c diag ****
cc        if (igr == igmax-1) then
cc          write (*,*) 'Plotting here at level',igr
cc
cc          call matvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd)
cc          call vecadd(igr,neq,nn,-1d0,wrk(isig),1d0,yy(isig))
cc          call MGplot(neq,wrk,igr,1,'fine2.bin')
cc        endif
c diag ****


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

        integer(4) :: igr

c     Local variables

        integer(4) :: nn,isig,depth1

c     Begin program

        if (outc.ge.1) write (*,*) 'Grid Level',igr

        if (.not.line_relax) then
          nn   = ntotv (igr)
          isig = istart(igr)

          depth1 = depth + 1

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
c     This routine performs plane/line smoothing instead of point smoothing.
c     This is done recursively, by calling 2D MG for planes and 1D MG for
c     lines. When reducing the dimensionality of MG, vectors remain of the
c     same length as in the original grid, and only matvec and smoothing 
c     operations use a decreased dimensionality. This allows one to employ
c     the SAME matvec and smoothing routines as for the original grid, with
c     minor modifications that restrict operations to the proper plane/line.
c
c     As vector dimensions are the same when going from 3d MG to 2D MG and
c     1D MG, MG pointers are not reallocated (and are consistent with the
c     original grid, defined in grid_params). However, the grid definition
c     contained in MGgrid is "collapsed", and this is the way MG knows it
c     should operate on a restricted subspace of mesh points. Therefore,
c     node positioning on the grid follows grid_params (i.e., nxv,nyv,nzv,etc)
c     whereas matvec and smoothing operations follow MGgrid.
c     -----------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igr

c     Local variables

        integer(4) :: nn,isig,depth1,i,j,k,iii,ig,iter,nsweep
     .               ,imin,imax,jmin,jmax,kmin,kmax
        integer(4) :: nxl(ngrid),nyl(ngrid),nzl(ngrid)

        real(8)    :: mag,mag0,omega
        real(8),allocatable,dimension(:) :: dx,rr

        type(grid_def) :: line_mg_grid

        logical        :: fpointers

c     Begin program

        nn   = ntotv (igr)
        isig = istart(igr)

        nxl = MGgrid%nxv
        nyl = MGgrid%nyv
        nzl = MGgrid%nzv

        nsweep = line_nsweep

        omega  = line_omega

        line_mg_grid%ngrdx = 0
        line_mg_grid%ngrdy = 0
        line_mg_grid%ngrdz = 0
        line_mg_grid%ngrid = 0

c     Schedule planes/lines

c       Check if at least 2D problem
        if (    (nxl(igr) > 1 .and. nyl(igr) > 1)
     .      .or.(nxl(igr) > 1 .and. nzl(igr) > 1)
     .      .or.(nyl(igr) > 1 .and. nzl(igr) > 1) ) then

c         Initialize quantities

          allocate(rr(nn),dx(nn))

          i = 1
          j = 1
          k = 1

c         Find initial residual

          call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
          call vecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

          if (outc > 0) then
            mag0 = sqrt(dot(igr,neq,nn,rr,rr))
            if (outc > 1) write (*,*)
            write (*,7) mag0
          endif

c         Planes/Lines in X-direction

cc            do iter=1,nsweep

          if (nxl(igr) > 1 .and. line_x) then

            !Collapse MG grid levels
            line_mg_grid            = MGgrid !Copy existing grid info
            line_mg_grid%nxv        = 1 !Collapse grid in x-direction
            line_mg_grid%mg_ratio_x = 1 !Collapse grid in x-direction

            do iter=1,nsweep

              if (outc > 1) then
                write (*,*)
                write (*,*) '******** X-relax sweep #',iter,'********'
              endif

              do i=1,nxl(igr)

                !Define plane/line for subsequent MG solvers at all grids
                line_mg_grid%iline(igr) = i
                do ig = igr+1,ngrid
                  line_mg_grid%iline(ig) = ( line_mg_grid%iline(ig-1)
     .                                      +MGgrid%mg_ratio_x(ig-1)-1 )
     .                                     /MGgrid%mg_ratio_x(ig-1)
                enddo

                !Output info
                if (outc > 1) then
                  write (*,*)
                  write (*,*) 'PLANE SOLVE: i=',i
     .                   ,'; ny x nz:',nyl(igr),'x',nzl(igr)
                endif

                dx = 0d0

                !Restrict grid pointers
                call saveMGgrid(line_mg_grid)

                !Solve plane/line (find solution update dx)
                call lineSolver(nn,line_mg_grid,rr,dx,igr,0)

                !Update solution (x = x + dx)
                call vecadd(igr,neq,nn,1d0,xx(isig),omega,dx)

                !Recover grid pointers
                call restoreMGgrid

              enddo

              !Update residual
              call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
              call vecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(dot(igr,neq,nn,rr,rr))
              if (outc > 1) write (*,5)
              write (*,10) mag,mag/mag0
              if (outc > 1) write (*,6)
              mag0 = mag
            endif

            call deallocateGridStructure(line_mg_grid)

          endif

c         Planes/Lines in Y-direction

          if (nyl(igr) > 1 .and. line_y) then

            line_mg_grid            = MGgrid !Copy existing grid info
            line_mg_grid%nyv        = 1 !Collapse grid in y-direction
            line_mg_grid%mg_ratio_y = 1 !Collapse grid in y-direction

            do iter=1,nsweep

              if (outc > 1) then
                write (*,*)
                write (*,*) '******** Y-relax sweep #',iter,'********'
              endif

              do j=1,nyl(igr)

                !Defines plane/line for subsequent MG solvers at all grids
                line_mg_grid%jline(igr) = j
                do ig = igr+1,ngrid
                  line_mg_grid%jline(ig) = ( line_mg_grid%jline(ig-1)
     .                                      +MGgrid%mg_ratio_y(ig-1)-1 )
     .                                     /MGgrid%mg_ratio_y(ig-1)
                enddo

                !Output info
                if (outc > 1) then
                  write (*,*)
                  write (*,*) 'PLANE SOLVE: j=',j
     .                   ,'; nx x nz:',nxl(igr),'x',nzl(igr)
                endif

                dx = 0d0

                !Restrict grid pointers
                call saveMGgrid(line_mg_grid)

                !Solve plane/line (find solution update)
                call lineSolver(nn,line_mg_grid,rr,dx,igr,0)

                !Update solution
                call vecadd(igr,neq,nn,1d0,xx(isig),omega,dx)

                !Recover grid pointers
                call restoreMGgrid

              enddo

              !Find residual
              call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
              call vecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(dot(igr,neq,nn,rr,rr))
              if (outc > 1) write (*,5)
              write (*,20) mag,mag/mag0
              if (outc > 1) write (*,6)
              mag0 = mag
            endif

            call deallocateGridStructure(line_mg_grid)

          endif

c         Planes/Lines in Z-direction

          if (nzl(igr) > 1 .and. line_z) then

            line_mg_grid            = MGgrid !Copy existing grid info
            line_mg_grid%nzv        = 1 !Collapse grid in z-direction
            line_mg_grid%mg_ratio_z = 1 !Collapse grid in z-direction

            do iter=1,nsweep

              if (outc > 1) then
                write (*,*)
                write (*,*) '******** Z-relax sweep #',iter,'********'
              endif

              do k=1,nzl(igr)

                !Defines plane/line for subsequent MG solvers at all grids
                line_mg_grid%kline(igr) = k
                do ig = igr+1,ngrid
                  line_mg_grid%kline(ig) = ( line_mg_grid%kline(ig-1)
     .                                      +MGgrid%mg_ratio_z(ig-1)-1 )
     .                                     /MGgrid%mg_ratio_z(ig-1)
                enddo

                !Output info
                if (outc > 1) then
                  write (*,*)
                  write (*,*) 'PLANE SOLVE: k=',k
     .                   ,'; nx x ny:',nxl(igr),'x',nyl(igr)
                endif

                dx = 0d0

                !Restrict grid pointers
                call saveMGgrid(line_mg_grid)

                !Solve plane/line (find solution update)
                call lineSolver(nn,line_mg_grid,rr,dx,igr,0)

                !Update solution
                call vecadd(igr,neq,nn,1d0,xx(isig),omega,dx)

                !Recover grid pointers
                call restoreMGgrid

              enddo

              !Find residual
              call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
              call vecadd(igr,neq,nn,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

            enddo

            if (outc > 0) then
              mag = sqrt(dot(igr,neq,nn,rr,rr))
              if (outc > 1) write (*,5)
              write (*,30) mag,mag/mag0
              if (outc > 1) write (*,6)
            endif

            call deallocateGridStructure(line_mg_grid)

          endif

cc            enddo

c         Deallocate variables

          deallocate(rr,dx)

c       Else smooth 1D problem (for MG line smoother)
        else

          if (outc > 1) then
            write (*,*)
            if (nxl(igr) > 1) then
              write (*,*) '1D LINE SOLVE on j =',mggrid%jline(igr)
     .               ,', k =',mggrid%kline(igr),' **** i = 1 :',nxl(igr)
            elseif (nyl(igr) > 1) then
              write (*,*) '1D LINE SOLVE on i =',mggrid%iline(igr)
     .               ,', k =',mggrid%kline(igr),' **** j = 1 :',nyl(igr)
            else
              write (*,*) '1D LINE SOLVE on i =',mggrid%iline(igr)
     .               ,', j =',mggrid%jline(igr),' **** k = 1 :',nzl(igr)
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

        integer(4)     :: igr,nn,guess
        real(8)        :: b(nn),x(nn)
        type(grid_def) :: mg_grid
        logical        :: prelax

c     Local variables

        integer(4)     :: nxl,nyl,nzl

c     Begin program

        nxl = mg_grid%nxv(igr)
        nyl = mg_grid%nyv(igr)
        nzl = mg_grid%nzv(igr)

c     Schedule solvers

c       Recursive MG (at least 2D or if 1D MG is wanted)\
        if  (    (nxl > 1 .and. nyl > 1)
     .       .or.(nxl > 1 .and. nzl > 1)
     .       .or.(nyl > 1 .and. nzl > 1)) then

          call lineMGsolver(nn,mg_grid,b,x,igr,guess)

c       GMRES 1D solver (matvec HARDWIRED to v_mtvc for now)
cc        elseif (line_solve == "gm") then
cc
cc          call lineGMsolver(nn,mg_grid,b,x,igr,guess)

c       MG 1D solver
        elseif (line_solve == "mg") then

          call lineMGsolver(nn,mg_grid,b,x,igr,guess,oned_solve=.true.)

c       GS 1D solver
        elseif (line_solve == "gs") then

          call lineGSsolver(nn,mg_grid,b,x,igr,guess)

c       JB 1D solver
        elseif (line_solve == "jb") then

          call lineJBsolver(nn,mg_grid,b,x,igr,guess)

        else

          write (*,*) 'Line solver undefined'
          write (*,*) 'Aborting...'
          stop

        endif

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
cc        integer(4)     :: igr,nn,guess
cc        real(8)        :: b(nn),x(nn)
cc        type(grid_def) :: mg_grid
cc        logical        :: prelax
cc
ccc     Local variables
cc
cc        integer(4) :: nxl,nyl,nzl,nnl,depth1
cc        integer(4) :: imin,imax,jmin,jmax,kmin,kmax
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

c     lineGSsolver
c     #####################################################################
      recursive subroutine lineGSsolver(nn,mg_grid,b,x,igr,guess)

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

        integer(4)     :: igr,nn,guess
        real(8)        :: b(nn),x(nn)
        type(grid_def) :: mg_grid
        logical        :: prelax

c     Local variables

        integer(4)     :: istart_sv(ngrid)
     .                   ,istartp_sv(ngrid)
     .                   ,istartb_sv(ngrid)

        real(8),allocatable,dimension(:,:) :: old_diag

        type (solver_options) :: options

        logical               :: fpointers

c     Begin program

c     Save parent MG grid configuration (shared with filial MG via mg_internal module)

cc        MGgrid_sv%ngrdx = 0
cc        MGgrid_sv%ngrdy = 0
cc        MGgrid_sv%ngrdz = 0
cc        MGgrid_sv%ngrid = 0
cc
cc        !Save pointers
cc        MGgrid_sv  = MGgrid
        istart_sv  = istart
        istartp_sv = istartp
        istartb_sv = istartb

        !Save diagonal
        allocate(old_diag(size(diag,1),size(diag,2)))
        old_diag = diag
        deallocate(diag)

c     'Prime' grid structure variable to pass to filial MG call

        options%mg_grid_def%ngrdx = 0
        options%mg_grid_def%ngrdy = 0
        options%mg_grid_def%ngrdz = 0
        options%mg_grid_def%ngrid = 0

        options%mg_grid_def = MGgrid_sv
          
c     Shift MG pointers to current grid level igr

        call transferMGPointers(igr)

c     Allocate new diagonal and transfer elements

        allocate(diag(size(old_diag,1)
     .               ,size(old_diag,2)-istart_sv(igr)+1))

        diag = old_diag(:,istart_sv(igr):size(old_diag,2))

c     Configure recursive plane/line GS solve

        options%iter                        = 100
        options%ncolors                     = 4
        options%omega                       = 1.0           

        options%tol                         = line_tol
        options%vcyc                        = line_vcyc
        options%igridmin                    = igridmin
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
        options%mg_grid_def                 = mg_grid

        options%vertex_based_relax          = vbr_mg
        options%fdiag                       = fdiag

c     Call plane/line GS

        call gs(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc-1,depth)

c     Recover configuration for parent level MG

        !Recover pointers
        istart  = istart_sv
        istartp = istartp_sv
        istartb = istartb_sv
        MGgrid  = mg_grid

        !Recover diagonal
        deallocate(diag)
        allocate(diag(size(old_diag,1),size(old_diag,2)))
        diag = old_diag
        deallocate(old_diag)

c     Deallocate memory

        call deallocateGridStructure(options%mg_grid_def)

      end subroutine lineGSsolver

c     lineJBsolver
c     #####################################################################
      recursive subroutine lineJBsolver(nn,mg_grid,b,x,igr,guess)

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

        integer(4)     :: igr,nn,guess
        real(8)        :: b(nn),x(nn)
        type(grid_def) :: mg_grid
        logical        :: prelax

c     Local variables

        integer(4)     :: istart_sv(ngrid)
     .                   ,istartp_sv(ngrid)
     .                   ,istartb_sv(ngrid)

        real(8),allocatable,dimension(:,:) :: old_diag

        type (solver_options) :: options

        logical               :: fpointers

c     Begin program

c     Save parent MG grid configuration (shared with filial MG via mg_internal module)

        !Save pointers
        istart_sv  = istart
        istartp_sv = istartp
        istartb_sv = istartb

        !Save diagonal
        allocate(old_diag(size(diag,1),size(diag,2)))
        old_diag = diag
        deallocate(diag)

c     'Prime' grid structure variable to pass to filial MG call

        options%mg_grid_def%ngrdx = 0
        options%mg_grid_def%ngrdy = 0
        options%mg_grid_def%ngrdz = 0
        options%mg_grid_def%ngrid = 0

        options%mg_grid_def = MGgrid_sv
          
c     Shift MG pointers to current grid level igr

        call transferMGPointers(igr)

c     Allocate new diagonal and transfer elements

        allocate(diag(size(old_diag,1)
     .               ,size(old_diag,2)-istart_sv(igr)+1))

        diag = old_diag(:,istart_sv(igr):size(old_diag,2))

c     Configure recursive plane/line JB solve

        options%iter                        = 100
        options%omega                       = 1d0           

        options%tol                         = line_tol
        options%vcyc                        = line_vcyc
        options%igridmin                    = igridmin
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
        options%mg_grid_def                 = mg_grid

        options%vertex_based_relax          = vbr_mg
        options%fdiag                       = fdiag

c     Call plane/line JB

        call jb(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc-1,depth)

c     Recover configuration for parent level MG

        !Recover pointers
        istart  = istart_sv
        istartp = istartp_sv
        istartb = istartb_sv
        MGgrid  = mg_grid

        !Recover diagonal
        deallocate(diag)
        allocate(diag(size(old_diag,1),size(old_diag,2)))
        diag = old_diag
        deallocate(old_diag)

c     Deallocate memory

        call deallocateGridStructure(options%mg_grid_def)

      end subroutine lineJBsolver

c     lineMGsolver
c     #####################################################################
      recursive subroutine lineMGsolver(nn,mg_grid,b,x,igr,guess
     .                                 ,oned_solve)

c     ---------------------------------------------------------------------
c     This routine performs a recursive MG solve on selected planes/lines.
c     In call:
c       * nn: total grid size at grid level igr
c       * mg_grid: MG grid definition structure that defines local grid patch
c       * b: rhs vector
c       * x: solution vector
c       * igr: grid level
c       * guess: whether initial guess is provided (guess=1) or not (guess=0).
c       * 1d_solve (optional):whether we are doing 1D MG.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4)     :: igr,nn,guess
        real(8)        :: b(nn),x(nn)
        type(grid_def) :: mg_grid
        logical,optional :: oned_solve

c     Local variables

        integer(4)     :: istart_sv(ngrid)
     .                   ,istartp_sv(ngrid)
     .                   ,istartb_sv(ngrid)

        real(8),allocatable,dimension(:,:) :: old_diag

        type (solver_options) :: options

        logical               :: fpointers,oned_slv

c     Begin program

        if (PRESENT(oned_solve)) then
          oned_slv = oned_solve
        else
          oned_slv = .false.
        endif

c     Save parent MG grid configuration (shared with filial MG via mg_internal module)

cc        MGgrid_sv%ngrdx = 0
cc        MGgrid_sv%ngrdy = 0
cc        MGgrid_sv%ngrdz = 0
cc        MGgrid_sv%ngrid = 0
cc
cc        !Save pointers
cc        MGgrid_sv  = MGgrid
        istart_sv  = istart
        istartp_sv = istartp
        istartb_sv = istartb

        !Save diagonal
        allocate(old_diag(size(diag,1),size(diag,2)))
        old_diag = diag
        deallocate(diag)

c     'Prime' grid structure variable to pass to filial MG call

        options%mg_grid_def%ngrdx = 0
        options%mg_grid_def%ngrdy = 0
        options%mg_grid_def%ngrdz = 0
        options%mg_grid_def%ngrid = 0

        options%mg_grid_def = MGgrid_sv
          
c     Shift MG pointers to current grid level igr

        call transferMGPointers(igr)

c     Allocate new diagonal and transfer elements

        allocate(diag(size(old_diag,1)
     .               ,size(old_diag,2)-istart_sv(igr)+1))

        diag = old_diag(:,istart_sv(igr):size(old_diag,2))

c     Configure recursive plane/line MG solve

        options%tol                         = line_tol
        options%vcyc                        = line_vcyc
        options%igridmin                    = igridmin
        if (oned_slv) then
          options%orderres                  = 0
          options%orderprol                 = 0
        else
          options%orderres                  = orderres
          options%orderprol                 = orderprol
        endif
        options%mg_mu                       = mu
        options%vol_res                     = volf        

        options%mg_coarse_solver_depth      = line_crse_depth
        options%mg_line_relax               = .true.
        options%mg_line_nsweep              = line_nsweep
        options%mg_line_vcyc                = line_vcyc
        options%mg_line_tol                 = line_tol
        options%mg_line_omega               = line_omega
        options%mg_line_coarse_solver_depth = line_crse_depth
        options%mg_line_x                   = line_x
        options%mg_line_y                   = line_y
        options%mg_line_z                   = line_z
        options%mg_grid_def                 = mg_grid

        options%vertex_based_relax          = vbr_mg
        options%fdiag                       = fdiag

c     Call plane/line MG

        call mg(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc-1,depth)

c     Recover configuration for parent level MG

        !Recover pointers
        istart  = istart_sv
        istartp = istartp_sv
        istartb = istartb_sv
        MGgrid  = mg_grid
cc        MGgrid  = MGgrid_sv
cc        call deallocateGridStructure(MGgrid_sv)

        !Recover diagonal
        deallocate(diag)
        allocate(diag(size(old_diag,1),size(old_diag,2)))
        diag = old_diag
        deallocate(old_diag)

c     Deallocate memory

        call deallocateGridStructure(options%mg_grid_def)

      end subroutine lineMGsolver

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

        integer(4)     :: igr

c     Local variables

        integer(4)     :: ig

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

        type(grid_def) :: mg_grid

c     Local variables

c     Begin program

        MGgrid_sv%ngrdx = 0
        MGgrid_sv%ngrdy = 0
        MGgrid_sv%ngrdz = 0
        MGgrid_sv%ngrid = 0

        MGgrid_sv  = MGgrid
        MGgrid     = mg_grid

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

        call deallocateGridStructure(MGgrid)

        MGgrid  = MGgrid_sv

        call deallocateGridStructure(MGgrid_sv)

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

        integer(4) :: igr

c     Local variables

        integer(4) :: nn,isig,depth1

c     Begin program

        if (outc.ge.1) write (*,*) 'Grid Level',igr

        nn   = ntotv(igr)
        isig = istart(igr)

        if (crsedpth /= 0) then
          depth1 = crsedpth
        else
          depth1 = depth + 1
        endif

        call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                ,guess2,outc,depth1)

      end subroutine coarseSolve

c     killmg
c     ###################################################################
      recursive subroutine killmg

        implicit none

        if (fdiag)  deallocate(diag)
        call deallocPointers(fpointers)
        call deallocateGridStructure(MGgrid)

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

      integer(4) :: neq,ntotc,nxc,nyc,nzc,ntotf,nxf,nyf,nzf,order,igc
     .             ,bcnd(6,neq)
      real(8)    :: xc(ntotc),xf(ntotf)

c Local variables
 
      real(8)    :: xxf(ntotf/neq),arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)

      integer(4) :: ic,jc,if,jf,iic,iif,i,j,k,ig,jg,kg
     .             ,ieq,nntotc,nntotf,igf

c Diag

c Begin program

      nntotc = ntotc/neq
      nntotf = ntotf/neq

c Unpack vector into array

      call mapMGVectorToArray(0,neq,xc,nxc,nyc,nzc,arrayc,igc,.false.)

c Impose boundary conditions (external)

      call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd
     .            ,iorder=min(order,3))

c Prolong arrays

      do ieq=1,neq

c     Use scalar prolongation (arrayc -> xxf)

        xxf = 0d0

        call prolong(xxf              ,nxf,nyf,nzf
     .              ,arrayc(:,:,:,ieq),nxc,nyc,nzc
     .              ,order,igc)

c     Repack prolonged vector

        do i=1,nntotf
          xf(neq*(i-1)+ieq)=xxf(i)
        enddo

      enddo

c End program

      end subroutine cprolong

*deck prolong
c######################################################################
      subroutine prolong(xf,nxf,nyf,nzf,arrayc,nxc,nyc,nzc
     .                  ,order,igc)
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

      integer(4) :: nxc,nyc,nzc,nxf,nyf,nzf,order,igc
      real(8)    :: arrayc(0:nxc+1,0:nyc+1,0:nzc+1),xf(nxf*nyf*nzf)

c Local variables

      integer(4) :: ic,if,jc,jf,kc,kf,igf,iic,iif
     .             ,icg,jcg,kcg,ifg,jfg,kfg
      integer(4) :: iminc,imaxc,jminc,jmaxc,kminc,kmaxc
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf

      logical    :: fpointers

c Extrapolation

      real(8)    :: xxf,yyf,zzf,ff

      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer(4) :: kx,ky,kz,nx,ny,nz,dim,flg
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
        xx = MGgrid%xx(icg-1:icg+imaxc-iminc+1)
        yy = MGgrid%yy(jcg-1:jcg+jmaxc-jminc+1)
        zz = MGgrid%zz(kcg-1:kcg+kmaxc-kminc+1)

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
              xxf = MGgrid%xx(ifg)
              yyf = MGgrid%yy(jfg)
              zzf = MGgrid%zz(kfg)

              xf(iif) = db3val(xxf,yyf,zzf,0,0,0,tx,ty,tz,nx,ny,nz
     .                        ,kx,ky,kz,bcoef,work)

            enddo
          enddo
        enddo

        deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

      endif

c End program

      call deallocPointers(fpointers)

      end subroutine prolong

*deck crestrict
c######################################################################
      subroutine crestrict(neq,xc,ntotc,nxc,nyc,nzc
     .                        ,xf,ntotf,nxf,nyf,nzf
     .                    ,order,igf,volf)
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
c       * volf (logical): whether vectors contain volume fractions.
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer(4) :: ntotc,nxc,nyc,nzc,ntotf,nxf,nyf,nzf,order,igf,neq

      real(8)    :: xc(ntotc),xf(ntotf)

      logical    :: volf

c Local variables
 
      real(8)    :: xxc(ntotc/neq),xxf(ntotf/neq)

      integer(4) :: ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

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
     .               ,order,igf,volf)

c     Repack restricted vector

        do i=1,nntotc
          xc(neq*(i-1)+ieq)=xxc(i)
        enddo

      enddo

c End program

      end subroutine crestrict

*deck restrict
c######################################################################
      subroutine restrict(xc,nxc,nyc,nzc,xf,nxf,nyf,nzf
     .                   ,order,igf,volf)
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
c       * volf (logical): whether values are volume-weighed.
c----------------------------------------------------------------------

      use mg_internal

      implicit none        !For safe fortran

c Call variables

      integer(4) :: nxc,nyc,nzc,nxf,nyf,nzf,order,igf

      real(8)    :: xc(nxc*nyc*nzc),xf(nxf*nyf*nzf)

      logical    :: volf

c Local variables
 
      integer(4) :: ic,if,jc,jf,kc,kf,iif,iic,igc,inbv
     .             ,icg,jcg,kcg,ifg,jfg,kfg
      integer(4) :: iminc,imaxc,jminc,jmaxc,kminc,kmaxc
     .             ,iminf,imaxf,jminf,jmaxf,kminf,kmaxf

      real(8)    :: volt,vol,arrayf(0:nxf+1,0:nyf+1,0:nzf+1)

      logical    :: fpointers

c Interpolation

      real(8)    :: xxc,yyc,zzc,ff

      real(8),allocatable,dimension(:) :: xx,yy,zz

      integer(4) ::  kx,ky,kz,nx,ny,nz,dim,flg
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

                    if (.not.volf) then
                      vol = volume(if,jf,kf,igf,igf,igf)
                      xc(iic) = xc(iic) + xf(iif)*vol
                      volt = volt + vol
                    else
                      xc(iic) = xc(iic) + xf(iif)
                    endif

                  enddo
                enddo
              enddo

              if (.not.volf) xc(iic) = xc(iic)/volt

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

        xx = MGgrid%xx(ifg:ifg+imaxf-iminf)
        yy = MGgrid%yy(jfg:jfg+jmaxf-jminf)
        zz = MGgrid%zz(kfg:kfg+kmaxf-kminf)

        call mapMGVectorToArray(0,1,xf,nxf,nyf,nzf,arrayf,igf,.false.)

c     Renormalize without volume fractions

        if (volf) then
          do kf = kminf,kmaxf
            do jf = jminf,jmaxf
              do if = iminf,imaxf
                arrayf(if,jf,kf) = arrayf(if,jf,kf)
     .                            /volume(if,jf,kf,igf,igf,igf)
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
              xxc = MGgrid%xx(icg)
              yyc = MGgrid%yy(jcg)
              zzc = MGgrid%zz(kcg)

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

              if (volf) xc(iic) = xc(iic)*volume(ic,jc,kc,igc,igc,igc)

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

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: iter,alloc_stat,isig,itr,nn,ieq,irbg1,nrbg1,igridc
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax,iming
      real(8)    :: omega0,omega10,omega01,tol
      logical    :: fdiag,fpointers,vbr

      integer(4) :: i,j,k,iv,jv,kv,if,jf,kf,nxf,nyf,nzf,iiv,ivg,ncolors
      integer(4) :: ii,iii,iib,iiib,iic,iig,iblock,jblock,kblock,nblk
      real(8)    :: mag0,mag1,mag,yy(ntot)

      real(8),allocatable, dimension(:)  :: dummy,rhs

c Begin program

      if (out.ge.2) write (*,*)

c Allocate pointers

      call allocPointers(neq,fpointers)

      if (fpointers) then          !JB is NOT called from MG
        igmax = igrid
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

c Read solver configuration

      call optionsConsistencyCheck

      iter   = options%iter
      omega0 = options%omega
c$$$      omega10= options%omega10
c$$$      omega01= options%omega01
      tol    = options%tol
      fdiag  = options%fdiag
      vbr    = options%vertex_based_relax
      ncolors= options%ncolors

c Initialize auxiliary variables

      nxf = nxv(igrid)
      nyf = nyv(igrid)
      nzf = nzv(igrid)

      if (vbr.or.vbr_mg) then
        nblk = nblock(igrid)
        isig = istartb(igrid)
        vbr  = .true.
      else
        nblk = 1
        isig = istart(igrid)
      endif

c Find diagonal for smoothers

      if (fdiag) then

        if (allocated(diag)) then !Diagonal is known
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          allocate(diag(neq*nblk,2*ntot*nblk))
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          allocate(diag(neq*nblk,2*ntot*nblk))
          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
     .                     ,ncolors)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c Jacobi iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag  = 0d0

        call matvec(0,neq,ntot,zz,yy,igrid,bcnd)

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
              dummy = matmul(diag(:,ivg+1:ivg+neq*nblk),rhs)

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
                dummy = matmul(diag(:,iig+1:iig+neq),rhs)

                !Update zz
                do ieq=1,neq
                  zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
                enddo

                mag = mag + sum(rhs*rhs)

              enddo
            enddo
          enddo

        endif

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          mag1 = mag
          if (out.ge.2) write (*,10) itr-1,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr-1,mag,mag/mag0,mag/mag1
          mag1 = mag
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.1) then
        call matvec(0,neq,ntot,zz,yy,igrid,bcnd)
        call vecadd(igrid,neq,ntot,-1d0,yy,1d0,rr)  !yy=rr-yy
        mag = sqrt(dot(igrid,neq,ntot,yy,yy))       !sqrt(yy*yy)

        if (itr.eq.0) then
          mag0 = mag
          mag1 = mag
        endif
      endif

      if (out.ge.2) then
        write (*,20) itr,mag,mag/mag0,mag/mag1
        write (*,*)
      elseif (out.eq.1) then
        write (*,10) itr,mag,mag/mag0
      endif

      options%iter_out = itr
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag)  deallocate(diag)
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

        if (options%vertex_based_relax .and. options%mg_line_relax) then
          write (*,*) 'Invalid setting for JB relaxation: cannot do'
          write (*,*) 'vertex-based and line-based relaxation',
     .                ' simultaneously!'
          write (*,*) 'Aborting...'
          stop
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

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: iter,alloc_stat,isig,itr,nn,ieq,iig
      integer(4) :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
      integer(4) :: i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max

      real(8)    :: omega0,tol
      logical    :: fdiag,fpointers,vbr

      integer(4) :: ic,jc,kc,iv,jv,kv,if,jf,kf,nxf,nyf,nzf,nxl,nyl,nzl
      integer(4) :: ii,iii,iib,iiib,iiv,ivg,iblock,jblock,kblock,nblk
      real(8)    :: mag0,mag1,mag,yy(ntot)

      integer(4) :: ncolors,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3

      real(8),allocatable, dimension(:) :: dummy,rhs

c Begin program

      if (out.ge.2) write (*,*)

c Set pointers

      call allocPointers(neq,fpointers)

      if (fpointers) then          !GS is NOT called from MG
        igmax = igrid
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2) write (*,*) 'Allocated pointers.'

c Read solver configuration

      call optionsConsistencyCheck

      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag
      ncolors= options%ncolors
      vbr    = options%vertex_based_relax

c Initialize auxiliary variables

      !Global grid size
      nxf = nxv(igrid)
      nyf = nyv(igrid)
      nzf = nzv(igrid)

      !Local grid size (different from global in line smoothing)
      nxl = MGgrid%nxv(igrid)
      nyl = MGgrid%nyv(igrid)
      nzl = MGgrid%nzv(igrid)

      if (vbr.or.vbr_mg) then
        nblk = nblock(igrid)
        isig = istartb(igrid)
        vbr  = .true.
      else
        nblk = 1
        isig = istart(igrid)
      endif

c Find diagonal for smoothers

      if (fdiag) then

        if (allocated(diag)) then !Diagonal is known
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          allocate(diag(neq*nblk,2*ntot*nblk))
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          allocate(diag(neq*nblk,2*ntot*nblk))
          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag
     .                     ,ncolors)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

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

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                call matvec(0,neq,ntot,zz,yy,igrid,bcnd)

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
                      dummy = matmul(diag(:,ivg+1:ivg+neq*nblk),rhs)

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

c         STANDARD RELAXATION
          else

            call scheduleColors(igrid,ncolors,nrbg1,nrbg2,nrbg3
     .                         ,i1min,i1max,j1min,j1max,k1min,k1max)

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                call matvec(0,neq,ntot,zz,yy,igrid,bcnd)

                do k1=k1min+mod((    irbg3-1),nrbg3),k1max,nrbg3
                  do j1=j1min+mod((k1+irbg2+irbg1-1),nrbg2),j1max,nrbg2
                    do i1=i1min+mod((j1+k1+irbg1-1),nrbg1),i1max,nrbg1

                      call scheduleIndices(igrid,i1,j1,k1,i,j,k)

                      !nxf,nyf,nzf are the original grid's even in line relaxation.

                      iii = neq*(i-1 + nxf*(j-1) + nxf*nyf*(k-1))

                      !Find new residual
                      do ieq = 1,neq
                        rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
                      enddo

                      !Multiply by D^-1 (stored in diag)
                      iig = iii + isig - 1
                      dummy = matmul(diag(:,iig+1:iig+neq),rhs)

                      !Update zz
                      do ieq=1,neq
                        zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
                      enddo

                      mag = mag + sum(rhs*rhs)

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

                    call matvec(-ii,neq,ntot,zz,yy,igrid,bcnd)

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
                dummy = matmul(diag(:,ivg+1:ivg+neq*nblk),rhs)

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

c           Forward pass
            call limits(0,nxf,nyf,nzf,igrid
     .                 ,imin,imax,jmin,jmax,kmin,kmax)

            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  !nxf,nyf,nzf are the original grid's even in line relaxation.
                  ii = i + nxf*(j-1) + nxf*nyf*(k-1)

                  call matvec(-ii,neq,ntot,zz,yy,igrid,bcnd)

                  iii = neq*(ii-1)

                  !Find new residual
                  do ieq = 1,neq
                    rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
                  enddo

                  !Multiply by D^-1 (stored in diag)
                  iig = iii + isig - 1
                  dummy = matmul(diag(:,iig+1:iig+neq),rhs)

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
cc                  call matvec(-ii,neq,ntot,zz,yy,igrid,bcnd)
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
cc                  dummy = matmul(diag(:,iig+1:iig+neq),rhs)
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

c     Check convergence

cc        call matvec(0,neq,ntot,zz,yy,igrid,bcnd)
cc        call vecadd(igrid,neq,ntot,-1d0,yy,1d0,rr)  !yy=rr-yy
cc        mag = sqrt(dot(igrid,neq,ntot,yy,yy))       !sqrt(yy*yy)

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          mag1 = mag
          if (out.ge.2) write (*,10) itr-1,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr-1,mag,mag/mag0,mag/mag1
          mag1 = mag
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.1) then
        call matvec(0,neq,ntot,zz,yy,igrid,bcnd)
        call vecadd(igrid,neq,ntot,-1d0,yy,1d0,rr)  !yy=rr-yy
        mag = sqrt(dot(igrid,neq,ntot,yy,yy))       !sqrt(yy*yy)
      endif

      if (out.ge.2) then
        write (*,20) itr,mag,mag/mag0,mag/mag1
        write (*,*)
      elseif (out.eq.1) then
        write (*,10) itr,mag,mag/mag0
      endif

      options%iter_out = itr
      options%tol_out  = mag/mag0

c End program

      deallocate (dummy,rhs)

      if (fdiag)  deallocate(diag)
      call deallocPointers(fpointers)

      return

 10   format (' GS Iteration:',i4,'; Residual:',1p,1e10.2,
     .        '; Ratio:',1e10.2)
 20   format (' GS Iteration:',i4,'; Residual:',1p,1e10.2,
     .        '; Ratio:',1e10.2,'; Damping:',1e10.2)

      contains

c     optionsConsistencyCheck
c     ###################################################################
        subroutine optionsConsistencyCheck

        if (options%vertex_based_relax .and. options%mg_line_relax) then
          write (*,*) 'Invalid setting for GS relaxation: cannot do'
          write (*,*) 'vertex-based and line-based relaxation',
     .                ' simultaneously!'
          write (*,*) 'Aborting...'
          stop
        elseif (options%vertex_based_relax .or. vbr_mg) then
          options%ncolors = 4
        elseif (options%mg_line_relax) then
          options%ncolors = 1
        endif
          
        end subroutine optionsConsistencyCheck

      end subroutine gs

c scheduleColors
c #######################################################################
      subroutine scheduleColors(igrid,ncolors,nc1,nc2,nc3
     .                         ,i1min,i1max,j1min,j1max,k1min,k1max)

c -----------------------------------------------------------------------
c     This routine schedules colors according to grid dimensions. In
c     call sequence:
c         * igrid: grid level (input)
c         * ncolors: maximum number of colors (input)
c         * nc1,nc2,nc3: number of colors in each grid dimension (output)
c         * i1min,i1max,j1min,j1max,k1min,k1max: local grid limits
c              to perform relaxation on (output)
c -----------------------------------------------------------------------

      use mg_internal

      implicit none

      integer(4) :: igrid,ncolors,nc1,nc2,nc3
     .             ,i1min,i1max,j1min,j1max,k1min,k1max
      
      integer(4) :: nxf,nyf,nzf,nxl,nyl,nzl,nclrs
     .             ,imin,imax,jmin,jmax,kmin,kmax

c Begin program

      nclrs = ncolors

      !Global grid size
      nxf = nxv(igrid)
      nyf = nyv(igrid)
      nzf = nzv(igrid)

      !Local grid size (different from global in line smoothing)
      nxl = MGgrid%nxv(igrid)
      nyl = MGgrid%nyv(igrid)
      nzl = MGgrid%nzv(igrid)

      call limits(0,nxf,nyf,nzf,igrid,imin,imax,jmin,jmax,kmin,kmax)

      if (nxl > 1 .and. nyl > 1 .and. nzl > 1) then

        i1min = imin
        i1max = imax
        j1min = jmin
        j1max = jmax
        k1min = kmin
        k1max = kmax

      elseif (nxl > 1 .and. nyl > 1 .and. nzl == 1) then

        if (nclrs > 4) nclrs = 4

        i1min = imin
        i1max = imax
        j1min = jmin
        j1max = jmax
        k1min = kmin
        k1max = kmax

      elseif (nxl > 1 .and. nyl == 1 .and. nzl > 1) then

        if (nclrs > 4) nclrs = 4

        i1min = kmin
        i1max = kmax
        j1min = imin
        j1max = imax
        k1min = jmin
        k1max = jmax

      elseif (nxl == 1 .and. nyl > 1 .and. nzl > 1) then

        if (nclrs > 4) nclrs = 4

        i1min = jmin
        i1max = jmax
        j1min = kmin
        j1max = kmax
        k1min = imin
        k1max = imax

      elseif (nxl > 1 .and. nyl == 1 .and. nzl == 1) then

        if (nclrs > 2) nclrs = 2

        i1min = imin
        i1max = imax
        j1min = jmin
        j1max = jmax
        k1min = kmin
        k1max = kmax

      elseif (nxl == 1 .and. nyl > 1 .and. nzl == 1) then

        if (nclrs > 2) nclrs = 2

        i1min = jmin
        i1max = jmax
        j1min = kmin
        j1max = kmax
        k1min = imin
        k1max = imax

      elseif (nxl == 1 .and. nyl == 1 .and. nzl > 1) then

        if (nclrs > 2) nclrs = 2

        i1min = kmin
        i1max = kmax
        j1min = imin
        j1max = imax
        k1min = jmin
        k1max = jmax

      endif

      select case(nclrs)
      case(2)
        nc1 = 2
        nc2 = 1
        nc3 = 1
      case(4)
        nc1 = 2
        nc2 = 2
        nc3 = 1
      case(8)
        nc1 = 2
        nc2 = 2
        nc3 = 2
      case default
        write (*,*) 'Unsupported number of colors in GS'
        write (*,*) 'Aborting...'
        stop
      end select

c End program

      end subroutine scheduleColors

c scheduleIndices
c #######################################################################
      subroutine scheduleIndices(igrid,i1,j1,k1,i,j,k)

c -----------------------------------------------------------------------
c     This routine schedules colored indices according to grid dimensions
c     and to the orderings in scheduleColors. In call:
c        * igrid: grid level
c        * i1,j1,k1: logical colored indices
c        * i,j,k: actual grid indices
c -----------------------------------------------------------------------

      use mg_internal

      implicit none

      integer(4) :: igrid,i1,j1,k1,i,j,k

      integer(4) :: nxl,nyl,nzl

c Begin program

      !Local grid size
      nxl = MGgrid%nxv(igrid)
      nyl = MGgrid%nyv(igrid)
      nzl = MGgrid%nzv(igrid)

      if (nxl > 1 .and. nyl > 1 .and. nzl > 1) then

        i = i1
        j = j1
        k = k1

      elseif (nxl > 1 .and. nyl > 1 .and. nzl == 1) then

        i = i1
        j = j1
        k = k1

      elseif (nxl > 1 .and. nyl == 1 .and. nzl > 1) then

        i = j1
        j = k1
        k = i1

      elseif (nxl == 1 .and. nyl > 1 .and. nzl > 1) then

        i = k1
        j = i1
        k = j1

      elseif (nxl > 1 .and. nyl == 1 .and. nzl == 1) then

        i = i1
        j = j1
        k = k1

      elseif (nxl == 1 .and. nyl > 1 .and. nzl == 1) then

        i = k1
        j = i1
        k = j1

      elseif (nxl == 1 .and. nyl == 1 .and. nzl > 1) then

        i = j1
        j = k1
        k = i1

      endif

      end subroutine scheduleIndices

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

      integer(4) :: neq,ntot,bbcnd(6,neq),igrid

      real(8)    :: diag1(neq,2*ntot)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot),mat(neq,neq),mat2(neq,neq)
     $             ,det
      integer(4) :: ii,jj,nn,ig,iig,isig
      integer(4) :: igr,ieq,alloc_stat
      logical    :: fpointers

c Begin program

c Allocate MG pointers

      call allocPointers(neq,fpointers)

      if (fpointers) igmax=grid_params%ngrid  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        write (*,*) 'Error in input of find_mf_diag_stc'
        write (*,*) 'Aborting...'
        stop
      endif

c Form diagonal

      do igr = igrid,igmax

c     Form diagonal terms for smoother

        nn  = ntotv(igr)

        isig = istart(igr)

        x1(1:nn) = 0d0

c     Finds block diagonals for neq equations.

        do ii = 1,ntotvp(igr)

          jj  = (ii-1)*neq + isig - 1

          do ieq = 1,neq

c         Find column vector corresponding to grid node ii and equation ieq

            x1(neq*(ii-1) + ieq) = 1d0

            call matvec(ii,neq,nn,x1,dummy,igr,bbcnd)

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

        case (2)

          do ii = 1,ntotvp(igr)
            iig = neq*(ii - 1) + isig - 1
            det = diag1(1,iig+1)*diag1(2,iig+2)
     .           -diag1(2,iig+1)*diag1(1,iig+2)
            mat(1,1) = diag1(2,iig+2)/det
            mat(1,2) =-diag1(2,iig+1)/det
            mat(2,1) =-diag1(1,iig+2)/det
            mat(2,2) = diag1(1,iig+1)/det

            diag1(:,iig+1:iig+neq) = mat
          enddo

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

      integer(4) :: neq,ntot,bbcnd(6,neq),igrid,ncolors

      real(8)    :: diag1(neq,2*ntot)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot),mat(neq,neq),mat2(neq,neq)
     .             ,det
      integer(4) :: ii,jj,nn,ig,iig,isig
      integer(4) :: igr,ieq,alloc_stat
      logical    :: fpointers

      integer(4) :: i,j,k,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3
      integer(4) :: i1,j1,k1,i1min,i1max,j1min,j1max,k1min,k1max

c Begin program

c Allocate MG pointers

      call allocPointers(neq,fpointers)

      if (fpointers) igmax=grid_params%ngrid  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        write (*,*) 'Error in input of find_mf_diag_colored'
        write (*,*) 'Aborting...'
        stop
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

          x1(1:nn) = 0d0

c       Finds block diagonals for neq equations.

          call scheduleColors(igr,ncolors,nrbg1,nrbg2,nrbg3
     .                       ,i1min,i1max,j1min,j1max,k1min,k1max)

          do ieq = 1,neq

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

cc                  do k=1+mod((irbg3-1),nrbg3),nzv(igr),nrbg3
cc                    do j=1+mod((k+irbg2+irbg1-1),nrbg2),nyv(igr),nrbg2
cc                      do i=1+mod((j+k+irbg1-1),nrbg1),nxv(igr),nrbg1
                  do k1=k1min+mod(irbg3-1,nrbg3),k1max,nrbg3
                    do j1=j1min+mod(k1+irbg2+irbg1-1,nrbg2),j1max,nrbg2
                      do i1=i1min+mod(j1+k1+irbg1-1,nrbg1),i1max,nrbg1

                        call scheduleIndices(igr,i1,j1,k1,i,j,k)

                        ii = i +nxv(igr)*(j-1) +nxv(igr)*nyv(igr)*(k-1)

                        jj = (ii-1)*neq + isig - 1

                        x1((ii-1)*neq+ieq) = 1d0
                      enddo
                    enddo
                  enddo

                  call matvec(0,neq,nn,x1,dummy,igr,bbcnd)

cc                  do k=1+mod((irbg3-1),nrbg3),nzv(igr),nrbg3
cc                    do j=1+mod((k+irbg2+irbg1-1),nrbg2),nyv(igr),nrbg2
cc                      do i=1+mod((j+k+irbg1-1),nrbg1),nxv(igr),nrbg1
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

c       Invert diagonal and store in diag1

          select case (neq)
          case (1)

            do ii = 1,ntotvp(igr)
              ig = ii + isig - 1
              diag1(neq,ig) = 1d0/diag1(neq,ig)
            enddo

          case (2)

            do ii = 1,ntotvp(igr)
              iig = neq*(ii - 1) + isig - 1
              det = diag1(1,iig+1)*diag1(2,iig+2)
     .             -diag1(2,iig+1)*diag1(1,iig+2)
              mat(1,1) = diag1(2,iig+2)/det
              mat(1,2) =-diag1(2,iig+1)/det
              mat(2,1) =-diag1(1,iig+2)/det
              mat(2,2) = diag1(1,iig+1)/det

              diag1(:,iig+1:iig+neq) = mat
            enddo

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
     .                       ,diag1,ncolors)
c--------------------------------------------------------------------
c     Finds diagonal elements for relaxation. Thisis done matrix-free
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
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer(4) :: neq,nblk,ntot,bbcnd(6,neq),igrid,ncolors

      real(8)    :: diag1(neq*nblk,2*ntot*nblk)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot),mag(3)
      integer(4) :: ii,jj,nn,ig,isig,size
      integer(4) :: iv,jv,kv,if,jf,kf,nxf,nyf,nzf
      integer(4) :: iii,icolb,irowb,icol,iiv,iig,ivg,irowmin,irowmax
     .             ,iblock ,jblock ,kblock
     .             ,iblock2,jblock2,kblock2,nblkg

      integer(4) :: igr,igc,ieq,alloc_stat
      logical    :: fpointers

c Begin program

c Allocate MG pointers

      call allocPointers(neq,fpointers)

      if (fpointers) igmax=grid_params%ngrid  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        write (*,*) 'Error in input of find_mf_diag'
        write (*,*) 'Aborting...'
        stop
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
                            call matvec(jj,neq,nn,x1,dummy,igr,bbcnd)

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

      integer(4) :: neq,ntot,bbcnd(6,neq),igrid

      real(8)    :: mat(ntot,ntot)

      external      matvec

c Local variables

      real(8)    :: x1(ntot),dummy(ntot)
      integer(4) :: ii,jj,nn
      integer(4) :: igr,ieq,alloc_stat
      logical    :: fpointers

c Begin program

      igr = igrid

c Allocate MG pointers

      call allocPointers(neq,fpointers)

c Consistency check

      nn  = ntotv(igr)

      if (nn /= ntot) then
        write (*,*) 'Error in input of find_mf_mat'
        write (*,*) 'Aborting...'
        stop
      endif

c Form diagonal 

      x1 = 0d0

      do ii = 1,nn

c       Find column vector corresponding to grid node ii and equation ieq

          x1(ii) = 1d0

          call matvec(0,neq,nn,x1,dummy,igr,bbcnd)

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

      integer(4) :: neq,igrid,bcnd(6,neq)

      external      matvec

c     Local variables

      real(8),allocatable,dimension(:)::x1,dummy,dummy2

      real(8)    :: dd1,dd2,err
      integer(4) :: nx,ny,nz,nn,ii,jj,i1,j1,i2,j2,ix1,iy1,ix2,iy2,ieq

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

          call matvec(0,neq,nn,x1,dummy,igrid,bcnd)

          call findBaseVector(ii,ieq,neq,nn,x1,0d0)

c       Compare column vector ii with corresponding row vector (intersect in
c       diagonal)

          do jj = ii,nn/neq

            call findBaseVector(jj,ieq,neq,nn,x1,1d0)

            call matvec(ii,neq,nn,x1,dummy2,igrid,bcnd)

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
 20   format (/,'Total relative error:',1pe10.3)

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

      integer(4) :: neq,ii,ieq,ntot

      real(8)    :: x1(ntot),coef

c Local variables

      integer(4) :: iii

c Begin program

      iii = neq*(ii-1) + ieq
      x1(iii) = coef

c End program

      end subroutine findBaseVector



c TODO list:
c
c 2) In weighed Jacobi, replace omega factors by a matrix-vector
c    operation.
c
c 3) Document i/o variables in all subroutines consistently.

c module mg_internal
c######################################################################
      module mg_internal

        use mlsolverSetup

        integer(4) :: ngrdx,ngrdy,ngrdz,ngrid,igmax

        real(8), allocatable, dimension(:,:) :: diag

        integer(4),dimension(:),allocatable ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z

        logical :: vbr_mg

        type (grid_def) :: MGgrid

      contains

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,mg_grid,fpointers)

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
      type(grid_def)      :: mg_grid

c     Local variables

      integer(4) :: i

c     Begin program

      fpointers = .false.

      if (.not.associated(mg_grid%xx)) mg_grid=grid_params     !Failsafe grid definition

      ngrdx = mg_grid%ngrdx
      ngrdy = mg_grid%ngrdy
      ngrdz = mg_grid%ngrdz
      ngrid = mg_grid%ngrid

c     Check if MG pointers are allocated

      if (.not.allocated(istart)) then

        allocate(istart (ngrid),ntotv (ngrid)
     .          ,istartp(ngrid),ntotvp(ngrid)
     .          ,istartb(ngrid),ntotb (ngrid)
     .          ,nxv(ngrid),nyv(ngrid),nzv(ngrid),nblock(ngrid)
     .          ,mg_ratio_x(ngrid),mg_ratio_y(ngrid),mg_ratio_z(ngrid))

        mg_ratio_x = mg_grid%mg_ratio_x
        mg_ratio_y = mg_grid%mg_ratio_y
        mg_ratio_z = mg_grid%mg_ratio_z

        nxv = mg_grid%nxv
        nyv = mg_grid%nyv
        nzv = mg_grid%nzv

cc        istartp = mg_grid%istartp

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

cc      else
cc
cc        !This is the only pointer redefined if pointers are allocated,
cc        !for the case of recursive plane/line smoothing
cc        mg_ratio_x = mg_grid%mg_ratio_x
cc        mg_ratio_y = mg_grid%mg_ratio_y
cc        mg_ratio_z = mg_grid%mg_ratio_z

      endif

c     End program

      end subroutine allocPointers

c     deallocPointers
c     #################################################################
      subroutine deallocPointers(fpointers)

c     -----------------------------------------------------------------
c     Deallocates pointers for block MG in fpointers (in, logical) is
c     true.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers

      if (fpointers) deallocate(istart,ntotv,istartp,ntotvp
     .                         ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                         ,mg_ratio_x,mg_ratio_y,mg_ratio_z)

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

      call allocPointers(neq,MGgrid,fpointers)

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
c       * mgvector (in): MG vector
c       * nx,ny,nz( in): grid dimensions
c       * array (out): mapped array
c       * igr (in): grid level
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

      call allocPointers(neq,MGgrid,fpointers)

      call limits(gpos,nx,ny,nz,igr,imin,imax,jmin,jmax,kmin,kmax)
cc      write (*,*) 'mapMGvector Loop limits:'
cc     .     ,imin,imax,jmin,jmax,kmin,kmax

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
      elseif (elem <= 0) then !Give 
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

        external   :: dgetri,dgetrf

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
c     Plots MG vector "mgv" at grid level "igr" using xdraw.
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

        call allocPointers(neq,MGgrid,fpointers)

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
          call contour(debug(1:nx,1:ny,1,ieq),nx,ny,0d0,xmax,0d0,ymax
     .                ,iflag+ieq-1,nunit)
        enddo

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
c     Contours arr in xdraw format
c     Notes:
c      put the next 2 lines in main
c      open(unit=nunit,file='contour.bin',form='unformatted') before
c      close(unit=nunit) after
c     ---------------------------------------------------------------------

c     Call variables

      integer(4) :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer(4) :: i,j

c     Begin program

      if(iopt.eq.0) then
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
      subroutine vecadd(igr,neq,coef1,vec1,coef2,vec2)

c     -------------------------------------------------------------------
c     Performs vec1 <- coef1*vec1 + coef2*vec2
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: igr,neq
      real(8)    :: coef1,coef2,vec1(*),vec2(*)

c     Local variables

      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
     .             ,i,j,k,iii,ieq

c     Begin program

      call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
     .               ,imin,imax,jmin,jmax,kmin,kmax)

cc          write (*,*) imin,imax,jmin,jmax,kmin,kmax

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            do ieq=1,neq
              iii = neq*(i-1 + nxv(igr)*(j-1)
     .                       + nxv(igr)*nyv(igr)*(k-1)) + ieq
              vec1(iii) = coef1*vec1(iii) + coef2*vec2(iii)
            enddo
          enddo
        enddo
      enddo

      end subroutine vecadd

c     dot
c     ###################################################################
      real(8) function dot(igr,neq,vec1,vec2)

c     -------------------------------------------------------------------
c     Performs scalar product (vec1,vec2)
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: igr,neq
      real(8)    :: vec1(*),vec2(*)

c     Local variables

      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
     .                 ,i,j,k,iii,ieq

c     Begin program

      call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
     .               ,imin,imax,jmin,jmax,kmin,kmax)

cc          write (*,*) imin,imax,jmin,jmax,kmin,kmax

      dot = 0d0

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            do ieq=1,neq
              iii = neq*(i-1 + nxv(igr)*(j-1)
     .                       + nxv(igr)*nyv(igr)*(k-1)) + ieq
              dot = dot + vec1(iii)*vec2(iii)
            enddo
          enddo
        enddo
      enddo

      end function dot

      end module mg_internal

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

      integer(4) :: iter,igridmin,vcyc,crsedpth,isig
      integer(4) :: orderres,orderprol,alloc_stat

      integer(4) :: guess2,outc,mu
      integer(4) :: izero,i,j,ii,ivcyc,nblk

      real(8)    :: xx(2*ntot),yy(2*ntot),wrk(2*ntot),rr(ntot)
      real(8)    :: rr0,rr1,mag,mag1,mgtol

      logical    :: fdiag,line_relax,fpointers,volf

c Begin program

      igridmin    = options%igridmin
      vcyc        = options%vcyc
      mgtol       = options%tol
      orderres    = options%orderres
      orderprol   = options%orderprol
      fdiag       = options%fdiag
      volf        = options%vol_res
      crsedpth    = options%mg_coarse_solver_depth
      mu          = options%mg_mu
      vbr_mg      = options%vertex_based_relax
      MGgrid      = options%mg_grid_def
      line_relax = options%mg_line_relax

c Consistency check

      if (crsedpth == 0 .and. vbr_mg) then
        write (*,*) 'Coarsest level solver not defined'
        write (*,*) 'Cannot do vertex-based relaxation in MG'
        write (*,*) 'Aborting...'
        stop
      endif

c Set pointers and find ngrid

      call allocPointers(neq,MGgrid,fpointers)

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
          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag)
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
        rr0 = sqrt(dot(igrid,neq,y,y))
      else
        call matvec(0,neq,ntot,x,rr,igrid,bcnd)
        call vecadd(igrid,neq,-1d0,rr,1d0,y)
        rr0 = sqrt(dot(igrid,neq,rr,rr))
      endif

      if (rr0.lt.1d-16*ntot) then
        if (out.ge.1) write (*,*) 'Initial solution seems exact in MG'
        call killmg
        return
      endif

cc      write (*,*) 'initial residual', rr0

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

        call vecadd(igrid,neq,-1d0,rr,1d0,y)
        mag = sqrt(dot(igrid,neq,rr,rr))

        mag1 = mag/rr1

        if (out.ge.3) then
          write (*,10) mag,mag1
        elseif (out.ge.2) then
          if (mu == 1) write (*,20) mag,mag1,ivcyc
          if (mu == 2) write (*,21) mag,mag1,ivcyc
        endif

c diag ****
cc      plot = 'n'
cc      write (*,*) 'Plot now?'
cc      read (*,'(a)') plot
cc      if (plot == 'y') then
cc        call MGplot(neq,xx,igrid,0,'debug.bin')
cc        call MGplot(neq,rr,1    ,1,'debug.bin')
cc        stop
cc      endif
cc      if (ivcyc == 4) then
cc        call MGplot(neq,xx,igrid,0,'debug.bin')
cc        call MGplot(neq,rr,1    ,1,'debug.bin')
cc        stop
cc      endif
c diag ****

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
 10   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4)
 20   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4,
     .          '; V-cycle #:'i3)
 21   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4,
     .          '; W-cycle #:'i3)

      contains

c     mcycle
c     ###################################################################
        recursive subroutine mcycle(igr)

          implicit none
          integer(4) :: igr,igc,isigc,isig,nn,nnc,imu

c diag ********        
          character*(1) :: plot
c diag ********

c       Begin program

          igc   = igr+1

          nn    = ntotv(igr)
          nnc   = ntotv(igc)

          isig  = istart(igr)
          isigc = istart(igc)

c       Relax error/solution on grid number igr/igrid (find new xx)

          call smooth(igr)

c       Evaluate residual (ie wrk = yy - A xx = yy - wrk )

          wrk(isig:isig+nn-1) = 0d0

          call matvec(0,neq,nn,xx(isig),wrk(isig),igr,bcnd)

          call vecadd(igr,neq,-1d0,wrk(isig),1d0,yy(isig))

cc          write (*,*) 'pause'
cc          pause

c diag ****
cc      plot = 'n'
cc      write (*,*) 'Plot now?'
cc      read (*,'(a)') plot
cc      if (plot == 'y') then
cc        call MGplot(neq,xx ,igr,0,'fine.bin')
cc        call MGplot(neq,wrk,igr,1,'fine.bin')
cc      endif
c diag ****

c       Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) to a coarser grid

          call crestrict(neq
     .                  ,yy(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .                  ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .                  ,orderres,igr,volf)

c diag ****
cc      if (plot == 'y') call MGplot(neq,yy,igc,0,'coarse.bin')
c diag ****

c       Initialize solution on coarse grid

          xx(isigc:isigc+nnc-1) = 0d0

c       If on coarsest grid, solve for error, else descend a grid level

          if (igc.eq.igmax) then
            call coarseSolve(igc)
          else
            do imu = 1,mu
              call mcycle(igc)
            enddo
          endif

c       Cycle back up to grid igr updating errors (xx)
c       with fixed R.H.S. (yy)

c diag ****
cc      if (plot == 'y') call MGplot(neq,xx,igc,1,'coarse.bin')
c diag ****

c       Prolong error (wrk = P * xx_1) to a finer grid

          call cprolong(neq
     .                 ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .                 ,xx(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .                 ,orderprol,igc,bcnd)

ccc diag ****
cc      if (plot == 'y') call MGplot(neq,wrk,igr,1,'fine.bin')
cccc      if (plot == 'y') stop
ccc diag ****

c       Update existing solution on grid igr (i.e. xx_igr): xx_igr = xx_igr + wrk 

          call vecadd(igr,neq,1d0,xx(isig),1d0,wrk(isig))

c       Relax updated solution on igr (i.e. xx_igr)

          call smooth(igr)

ccc diag ****
cc      if (plot == 'y') call MGplot(neq,xx,igr,1,'fine.bin')
cc      if (plot == 'y') stop
ccc diag ****

        end subroutine mcycle

c       smooth
c       ###################################################################
        recursive subroutine smooth(igr)

          implicit none
          integer :: igr,nn,isig,depth1

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

        end subroutine smooth

c       linesmooth
c       ###################################################################
        recursive subroutine linesmooth(igr)

c       -------------------------------------------------------------------
c       This routine performs plane/line smoothing instead of point smoothing.
c       This is done recursively, by calling 2D MG for planes and 1D MG for
c       lines. When reducing the dimensionality of MG, vectors remain of the
c       same length as in the original grid, and only matvec and smoothing 
c       operations use a decreased dimensionality. This allows one to employ
c       the SAME matvec and smoothing routines as for the original grid, with
c       minor modifications that restrict operations to the proper plane/line.
c
c       As vector dimensions are the same when going from 3d MG to 2D MG and
c       1D MG, MG pointers are not reallocated (and are consistent with the
c       original grid, defined in grid_params). However, the grid definition
c       contained in MGgrid is "collapsed", and this is the way MG knows it
c       should operate on a restricted subspace of mesh points. Therefore,
c       node positioning on the grid follows grid_params (i.e., nxv,nyv,nzv,etc)
c       whereas matvec and smoothing operations follow MGgrid.
c       ---------------------------------------------------------------------

          implicit none

c       Call variables

          integer(4) :: igr

c       Local variables

          integer(4) :: nn,isig,depth1,i,j,k,iii,ig,iter,nsweep
     .                 ,imin,imax,jmin,jmax,kmin,kmax
          integer(4),allocatable,dimension(:) :: nxl,nyl,nzl

          real(8)    :: mag,mag0,omega
          real(8),allocatable,dimension(:) :: dx,rr

          type(grid_def) :: line_mg_grid

          logical        :: fpointers

c       Begin program

          nn   = ntotv (igr)
          isig = istart(igr)

          allocate(nxl(ngrid),nyl(ngrid),nzl(ngrid))

          nxl = MGgrid%nxv
          nyl = MGgrid%nyv
          nzl = MGgrid%nzv

          nsweep = 1

          omega = 1d0

c       Schedule planes/lines

c         Check if at least 2D problem
          if (    (nxl(igr) > 1 .and. nyl(igr) > 1)
     .        .or.(nxl(igr) > 1 .and. nzl(igr) > 1)
     .        .or.(nyl(igr) > 1 .and. nzl(igr) > 1) ) then

c           Initialize quantities

            allocate(rr(nn),dx(nn))

cc            rr = 0d0
cc            dx = 0d0

cc            mag = 0d0

            i = 1
            j = 1
            k = 1

cc            do iter=1,nsweep

c           Planes/Lines in X-direction

            if (nxl(igr) > 1) then

              !Collapse MG grid levels
              line_mg_grid            = MGgrid !Copy existing grid info
              line_mg_grid%nxv        = 1 !Collapse grid in x-direction
              line_mg_grid%mg_ratio_x = 1 !Collapse grid in x-direction

              do iter=1,nsweep

               do i=1,nxl(igr)

                !Define plane/line for subsequent MG solvers at all grids
                line_mg_grid%iline(igr) = i
                do ig = igr+1,ngrid
                  line_mg_grid%iline(ig) = ( line_mg_grid%iline(ig-1)
     .                                      +MGgrid%mg_ratio_x(ig-1)-1 )
     .                                     /MGgrid%mg_ratio_x(ig-1)
                enddo

                !Find residual
                call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
                call vecadd(igr,neq,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

                if (outc > 0) then

                  !Find magnitude of residual and initialize reference
cc                  mag = sqrt(dot(igr,neq,rr,rr))
                  if (i == 1 .and. iter == 1)
     .                 mag0 = sqrt(dot(igr,neq,rr,rr))

                  !Output info
                  write (*,*)
                  write (*,*) 'PLANE SOLVE: i=',i
     .                   ,'; ny x nz:',nyl(igr),nzl(igr)
                endif

cc                call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
cc     .                     ,imin,imax,jmin,jmax,kmin,kmax)
cc
cc                write (*,*)
cc                write (*,*) '***************************************'
cc                write (*,*) 'i=',i
cc     .                   ,'; Solve plane ny x nz:',nyl(igr),nzl(igr)
cc     .                   ,'; Problem size:',neq*nyl(igr)*nzl(igr)
cc     .                   ,'; Loop limits:',imin,imax,jmin,jmax,kmin,kmax
cc                write (*,*) '***************************************'
cc                pause

                !Solve plane/line (find solution update dx)
                call lineMGsolver(nn,line_mg_grid,rr,dx,igr,0)
cc                call lineMGsolver(nn,line_mg_grid,yy(isig),xx(isig),igr
cc     .                           ,1)

                !Update solution (x = x + dx)
                call vecadd(igr,neq,1d0,xx(isig),omega,dx)

               enddo

              enddo

              if (outc > 0) then
                mag = sqrt(dot(igr,neq,rr,rr))
                write (*,10) mag,mag/mag0
              endif

            endif

c           Planes/Lines in Y-direction

            if (nyl(igr) > 1) then

              line_mg_grid            = MGgrid !Copy existing grid info
              line_mg_grid%nyv        = 1      !Collapse grid in y-direction
              line_mg_grid%mg_ratio_y = 1      !Collapse grid in y-direction

              do iter=1,nsweep

               do j=1,nyl(igr)

                !Defines plane/line for subsequent MG solvers at all grids
                line_mg_grid%jline(igr) = j
                do ig = igr+1,ngrid
                  line_mg_grid%jline(ig) = ( line_mg_grid%jline(ig-1)
     .                                      +MGgrid%mg_ratio_y(ig-1)-1 )
     .                                     /MGgrid%mg_ratio_y(ig-1)
                enddo

                !Find residual
                call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
                call vecadd(igr,neq,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

                if (outc > 0) then

                  !Find magnitude of residual and initialize reference
cc                  mag = sqrt(dot(igr,neq,rr,rr))
                  if (j == 1 .and. iter == 1)
     .                 mag0 = sqrt(dot(igr,neq,rr,rr))

                  !Output info
                  write (*,*)
                  write (*,*) 'PLANE SOLVE: j=',j
     .                   ,'; nx x nz:',nxl(igr),nzl(igr)
                endif

cc                call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
cc     .                     ,imin,imax,jmin,jmax,kmin,kmax)
cc
cc                write (*,*)
cc                write (*,*) '***************************************'
cc                write (*,*) 'j=',j
cc     .                   ,'; Solve plane nx x nz:',nxl(igr),nzl(igr)
cc     .                   ,'; Problem size:',neq*nxl(igr)*nzl(igr)
cc     .                   ,'; Loop limits:',imin,imax,jmin,jmax,kmin,kmax
cc                write (*,*) '***************************************'
cc                pause

                !Solve plane/line (find solution update)
                call lineMGsolver(nn,line_mg_grid,rr,dx,igr,0)
cc                call lineMGsolver(nn,line_mg_grid,yy(isig),xx(isig),igr
cc     .                           ,1)

                !Update solution
                call vecadd(igr,neq,1d0,xx(isig),omega,dx)

               enddo

              enddo

              if (outc > 0) then
                mag = sqrt(dot(igr,neq,rr,rr))
                write (*,20) mag,mag/mag0
              endif
            endif

c           Planes/Lines in Z-direction

            if (nzl(igr) > 1) then

              line_mg_grid            = MGgrid !Copy existing grid info
              line_mg_grid%nzv        = 1      !Collapse grid in z-direction
              line_mg_grid%mg_ratio_z = 1      !Collapse grid in z-direction

              do iter=1,nsweep

               do k=1,nzl(igr)

                !Defines plane/line for subsequent MG solvers at all grids
                line_mg_grid%kline(igr) = k
                do ig = igr+1,ngrid
                  line_mg_grid%kline(ig) = ( line_mg_grid%kline(ig-1)
     .                                      +MGgrid%mg_ratio_z(ig-1)-1 )
     .                                     /MGgrid%mg_ratio_z(ig-1)
                enddo

                !Find residual
                call matvec(0,neq,nn,xx(isig),rr,igr,bcnd)
                call vecadd(igr,neq,-1d0,rr,1d0,yy(isig)) !rr = yy(isig:isig+nn-1)-rr

                if (outc > 0) then

                  !Find magnitude of residual and initialize reference
cc                  mag = sqrt(dot(igr,neq,rr,rr))
                  if (k == 1 .and. iter == 1)
     .                 mag0 = sqrt(dot(igr,neq,rr,rr))

                  !Output info
                  write (*,*)
                  write (*,*) 'PLANE SOLVE: k=',k
     .                   ,'; nx x ny:',nxl(igr),nyl(igr)
                endif

cc                call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
cc     .                     ,imin,imax,jmin,jmax,kmin,kmax)
cc                write (*,*) 'k=',k
cc     .                   ,'; Solve plane nx x ny:',nxl(igr),nyl(igr)
cc     .                   ,'; Problem size:',neq*nxl(igr)*nyl(igr)
cc     .                   ,'; Loop limits:',imin,imax,jmin,jmax,kmin,kmax
cc                pause

                !Solve plane/line (find solution update)
                call lineMGsolver(nn,line_mg_grid,rr,dx,igr,0)
cc                call lineMGsolver(nn,line_mg_grid,yy(isig),xx(isig),igr,1)

                !Update solution
                call vecadd(igr,neq,1d0,xx(isig),omega,dx)
               enddo

              enddo

              if (outc > 0) then
                mag = sqrt(dot(igr,neq,rr,rr))
                write (*,30) mag,mag/mag0
              endif
            endif

cc            enddo

c           Deallocate variables

            deallocate(rr,dx)

c         Else smooth 1D problem
          else

cc            call limits(0,nxv(igr),nyv(igr),nzv(igr),igr
cc     .                     ,imin,imax,jmin,jmax,kmin,kmax)
cc            write (*,*)   'Smooth line:',nxl(igr),nyl(igr),nzl(igr)
cc     .                 ,'; Problem size:',neq*nxl(igr)*nyl(igr)*nzl(igr)
cc     .                 ,'; Loop limits:',imin,imax,jmin,jmax,kmin,kmax
cc            pause

            if (outc > 1) then
              write (*,*)
              write (*,*) '1D LINE SMOOTH on',nxl(igr),nyl(igr),nzl(igr)
            endif

            depth1 = depth + 1

            call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                    ,guess2,outc,depth1)

          endif

 10   format(/,'******************************************************',
     .       /,' X-plane relax. Residual:',1p1e10.2,
     .         '; Ratio:',1p1e10.2,/
     .        ,'******************************************************')
 20   format(/,'******************************************************',
     .       /,' Y-plane relax. Residual:',1p1e10.2,
     .         '; Ratio:',1p1e10.2/
     .        ,'******************************************************')
 30   format(/,'******************************************************',
     .       /,' Z-plane relax. Residual:',1p1e10.2,
     .         '; Ratio:',1p1e10.2/
     .        ,'******************************************************')

        end subroutine linesmooth

c       lineMGsolver
c       ###################################################################
        recursive subroutine lineMGsolver(nn,mg_grid,b,x,igr,guess)

c       -------------------------------------------------------------------
c       This routine performs a recursive MG solve on selected planes/lines.
c       ---------------------------------------------------------------------

          implicit none

c       Call variables

          integer(4)     :: igr,nn,guess
          real(8)        :: b(nn),x(nn)
          type(grid_def) :: mg_grid
          logical        :: prelax

c       Local variables

          integer(4)     :: istart_sv(ngrid)
     .                     ,istartp_sv(ngrid)
     .                     ,istartb_sv(ngrid)

          type (solver_options) :: options

          type(grid_def)        :: MGgrid_sv

          logical               :: fpointers

c       Begin program

c       Save parent MG grid configuration

          MGgrid_sv  = MGgrid
          istart_sv  = istart
          istartp_sv = istartp
          istartb_sv = istartb

c       Shift MG pointers to current grid level igr

          call transferPointers(igr)

c       Deallocate pointers of parent MG call

cc          call deallocPointers(fpointers) 

c       Configure recursive plane/line MG solve

          options%tol                    = mgtol
          options%vcyc                   = 1
          options%igridmin               = igridmin
cc          options%igridmin               = 2
          options%orderres               = orderres
          options%orderprol              = orderprol
          options%mg_mu                  = mu
          options%vol_res                = volf        
cc          options%diag                   => dg
          options%mg_coarse_solver_depth = crsedpth
          options%mg_coarse_solver_depth = 0
          options%mg_line_relax         = .true.
          options%mg_grid_def            = mg_grid
          options%vertex_based_relax     = vbr_mg

c       Call plane/line MG

          call mg(neq,nn,b,x,matvec,options,igr,bcnd,guess,outc,depth)

c       Recover grid definition for parent level MG

          MGgrid  = MGgrid_sv
          istart  = istart_sv
          istartp = istartp_sv
          istartb = istartb_sv

c       Allocate pointers of parent MG call

cc          call allocPointers(neq,MGgrid,fpointers)

        end subroutine lineMGsolver

c       transferPointers
c       ###############################################################
        subroutine transferPointers(igr)

c       ---------------------------------------------------------------
c       Transfers MG pointers from upper grid level to subsequent grid
c       levels. This is done by shifting pointer arrays by one grid,
c       removing the finest grid level each time. The MG grid level
c       definition mg_grid_sv is used as a reference.
c       ---------------------------------------------------------------

          implicit none        !For safe fortran

c       Call variables

          integer(4)     :: igr

c       Local variables

          integer(4)     :: ig

c       Begin program

c       If not at finest grid level, shift MG pointers

          if (igr > 1) then

            do ig=igr+1,ngrid
              istartp(ig) = istartp(ig)-istartp(ig-1)
              istart (ig) = istart (ig)-istart (ig-1)
              istartb(ig) = istartb(ig)-istartb(ig-1)
            enddo
            istartp(igr) = 1
            istart (igr) = 1
            istartb(igr) = 1

          endif

c       End program

        end subroutine transferPointers

c       coarseSolve
c       ###################################################################
        recursive subroutine coarseSolve(igr)

          implicit none
          integer :: igr,nn,isig,depth1

          if (outc.ge.1) write (*,*) 'Grid Level',igr

          nn   = ntotv(igr)
          isig = istart(igr)

          if (crsedpth /= 0) then
            depth1 = crsedpth
          else
            depth1 = depth + 1
          endif

          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                  ,guess2,outc,depth1)

        end subroutine coarseSolve

c       killmg
c       ###################################################################
        recursive subroutine killmg

          implicit none

          if (fdiag)  deallocate(diag)
          call deallocPointers(fpointers)

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
c       * nxf,nyf,nzf(int): dimensions of fine grid
c       * xc (real): vector in coarse grid (not a MG vector)
c       * nxc,nyc,nzc(int): dimensions of coarse grid
c       * order (int): order of interpolation (0-arbitrary)
c           If order = 0, it employs simple injection.
c           If order > 0, it employs spline interpolation.
c       * igc (int): coarse grid level identifier
c       * bcnd (int array): boundary condition info.
c----------------------------------------------------------------------

      use mg_internal

      implicit none            ! For safe Fortran

c Call variables

      integer(4) :: neq,ntotc,nxc,nyc,nzc,ntotf,nxf,nyf,nzf,order,igc
     .             ,bcnd(6,neq)
      real(8)    :: xc(ntotc),xf(ntotf)

c Local variables
 
      real(8)    :: xxf(ntotf/neq),arrayc(0:nxc+1,0:nyc+1,0:nzc+1,neq)

      integer(4) :: ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

      nntotc = ntotc/neq
      nntotf = ntotf/neq

c Unpack vector into array

      call mapMGVectorToArray(0,neq,xc,nxc,nyc,nzc,arrayc,igc,.false.)

c Impose boundary conditions (external)

      call setMGBC(0,neq,nxc,nyc,nzc,igc,arrayc,bcnd)

c Restric arrays

      do ieq=1,neq

c     Use scalar prolongation (arrayc -> xxf)

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

cc      real(8)    :: xx(nxc+2),yy(nyc+2),zz(nzc+2)

      integer(4) ::  kx,ky,kz,nx,ny,nz,dim,flg
      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db3val
      external      db3val

c Begin program

      call allocPointers(1,MGgrid,fpointers)

c Define fine grid

      igf = igc - 1

      call limits(0,nxc,nyc,nzc,igc,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

c Injection

      if (order.eq.0) then

cc        do kc = 1,nzc
cc          do jc = 1,nyc
cc            do ic = 1,nxc
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

cc        call getMGmap(1,1,1,igc,igc,igc,icg,jcg,kcg)
cc        xx(1:nxc+2) = MGgrid%xx(icg-1:icg+nxc)
cc        yy(1:nyc+2) = MGgrid%yy(jcg-1:jcg+nyc)
cc        zz(1:nzc+2) = MGgrid%zz(kcg-1:kcg+nzc)

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

cc        do kf = 1,nzf
cc          do jf = 1,nyf
cc            do if = 1,nxf
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

      call allocPointers(1,MGgrid,fpointers)

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
cc          do kf = 1,nzf
cc            do jf = 1,nyf
cc              do if = 1,nxf
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
cc        nx = nxc
cc        ny = nyc
cc        nz = nzc
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
     .               ,tz,nz,kz,bcoef,q,work)
          deallocate(q)
        elseif (nx == 1 .and. nz == 1) then
          allocate(q((2*ky-1)*ny))
          inbv=1
          call dbknot(yy,ny,ky,ty)
          call dbintk(yy,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,ty,ny,ky,bcoef,q,work)
          deallocate(q)
        elseif (ny == 1 .and. nz == 1) then
          allocate(q((2*kx-1)*nx))
          inbv=1
          call dbknot(xx,nx,kx,tx)
          call dbintk(xx,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,tx,nx,kx,bcoef,q,work)
          deallocate(q)
        elseif (nx == 1) then
cc          call db2ink(yy,ny,zz,nz,arrayf(1:nx,1:ny,1:nz)
          call db2ink(yy,ny,zz,nz
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,ny,ky,kz,ty,tz,bcoef,work,flg)
        elseif (ny == 1) then
cc          call db2ink(xx,nx,zz,nz,arrayf(1:nx,1:ny,1:nz)
          call db2ink(xx,nx,zz,nz
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,nx,kx,kz,tx,tz,bcoef,work,flg)
        elseif (nz == 1) then
cc          call db2ink(xx,nx,yy,ny,arrayf(1:nx,1:ny,1:nz)
          call db2ink(xx,nx,yy,ny
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,nx,kx,ky,tx,ty,bcoef,work,flg)
        else
cc          call db3ink(xx,nx,yy,ny,zz,nz,arrayf(1:nx,1:ny,1:nz)
          call db3ink(xx,nx,yy,ny,zz,nz
     .                     ,arrayf(iminf:imaxf,jminf:jmaxf,kminf:kmaxf)
     .               ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)
        endif

        call limits(0,nxc,nyc,nzc,igc
     .             ,iminc,imaxc,jminc,jmaxc,kminc,kmaxc)

        do kc = kminc,kmaxc
          do jc = jminc,jmaxc
            do ic = iminc,imaxc
cc        do kc = 1,nzc
cc          do jc = 1,nyc
cc            do ic = 1,nxc

              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)

              call getMGmap(ic,jc,kc,igc,igc,igc,icg,jcg,kcg)
              xxc = MGgrid%xx(icg)
              yyc = MGgrid%yy(jcg)
              zzc = MGgrid%zz(kcg)

              if (nx == 1 .and. ny == 1) then
                xc(iic) = dbvalu(tz,bcoef,nz,kz,0,zzc,inbv,work)
              elseif (nx == 1 .and. nz == 1) then
                xc(iic) = dbvalu(ty,bcoef,ny,ky,0,yyc,inbv,work)
              elseif (ny == 1 .and. nz == 1) then
                xc(iic) = dbvalu(tx,bcoef,nx,kx,0,xxc,inbv,work)
              elseif (nx == 1) then
                xc(iic)=db2val(yyc,zzc,0,0,ty,tz,ny,nz,ky,kz,bcoef,work)
              elseif (ny == 1) then
                xc(iic)=db2val(xxc,zzc,0,0,tx,tz,nx,nz,kx,kz,bcoef,work)
              elseif (nz == 1) then
                xc(iic)=db2val(xxc,yyc,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)
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

ccc jb
ccc#######################################################################
cc      recursive subroutine jb(neq,ntot,rr,zz,matvec,options,igrid,bcnd
cc     .                       ,guess,out,depth)
ccc--------------------------------------------------------------------
ccc     Matrix-free Jacobi routine to solve Azz = rr. Call variables:
ccc       * neq: number of equations (system JB)
ccc       * ntot: grid dimension
ccc       * rr,zz: rhs, solution vectors
ccc       * matvec: matrix-free matvec product (external)
ccc       * options: structure containing solver defs.
ccc       * igrid: grid level in MG applications
ccc       * guess: 0 -> no initial guess; 1 -> initial guess provided.
ccc       * out: convergence info output on screen if out > 1. 
ccc       * depth: integer specifying solver depth in solver_queue
ccc                definitions
ccc--------------------------------------------------------------------
cc
cc      use mlsolverSetup
cc
cc      use mg_internal
cc
cc      implicit none       !For safe fortran
cc
ccc Call variables
cc
cc      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
cc      real(8)    :: rr(ntot),zz(ntot)
cc
cc      type (solver_options) :: options
cc
cc      external     matvec
cc
ccc Local variables
cc
cc      integer(4) :: iter,alloc_stat,isig
cc      real(8)    :: omega0,omega10,omega01,tol
cc      logical    :: fdiag,fpointers
cc
cc      integer(4) :: i,j,itr,nn,ieq,nx,ny
cc      integer(4) :: ii,iii,iip,iim,jj,jjp,jjm
cc      integer(4) :: ig,ipg,img,jg,jpg,jmg,iig
cc      real(8)    :: mag0,mag1,mag,yy(ntot),delta
cc      real(8)    :: aw,ae,an,as
cc
cc      real(8),allocatable, dimension(:)  :: dummy,rhs
cc
ccc Begin program
cc
ccc Read solver configuration
cc
cc      iter   = options%iter
cc      omega0 = options%omega
cc      omega10= options%omega10
cc      omega01= options%omega01
cc      tol    = options%tol
cc      fdiag  = options%fdiag
cc
cc      if (out.ge.2) write (*,*)
cc
ccc Allocate pointers
cc
cc      call allocPointers(neq,fpointers)
cc
cc      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'
cc
cc      isig = istart(igrid)
cc
ccc Find diagonal for smoothers
cc
cc      if (fdiag) then
cc
cc        if (allocated(diag)) then !Diagonal is known
cc          if (out.ge.2) write (*,*) 'Diagonal already allocated'
cc          fdiag = .false.
cc        elseif (associated(options%diag)) then !Diagonal provided externally
cc          if (out.ge.2) write (*,*) 'Diagonal externally provided'
cc          allocate(diag(neq,2*ntot))
cc          diag = options%diag
cc        else                      !Form diagonal
cc          if (out.ge.2) write (*,*) 'Forming diagonal...'
cc          allocate(diag(neq,2*ntot))
cc          call find_mf_diag(neq,ntot,matvec,igrid,bcnd,diag)
cc          if (out.ge.2) write (*,*) 'Finished!'
cc        endif
cc
cc      endif
cc
ccc Preparation for iteration
cc
cc      nn = ntot
cc
cc      if (guess.eq.0) zz = 0d0
cc
ccc Jacobi iteration
cc
cc      allocate(dummy(neq),rhs(neq))
cc
cc      do itr=1,iter
cc
cc        mag1 = mag
cc        mag  = 0d0
cc
ccc$$$            nx = nxv(igrid)
ccc$$$            ny = nyv(igrid)
ccc$$$
ccc$$$            do i = 1,nx
ccc$$$              do j = 1,ny
ccc$$$                call shiftcoeffs  (i,j,bcnd(1,1))
ccc$$$
ccc$$$                call shiftIndices(i,j,ii,iim,iip,jj,jjm,jjp,nx,ny,1
ccc$$$     .                            ,bcnd(1,1))
ccc$$$                call shiftIndices(i,j,ig,img,ipg,jg,jmg,jpg,nx,ny,isig
ccc$$$     .                            ,bcnd(1,1))
ccc$$$
ccc$$$                zz(ii) =zz(ii)
ccc$$$     .                 +omega0*     (rr(ii) -yy(ii) )/diag(neq,ig )
ccc$$$     .                 +omega10*(ae*(rr(iip)-yy(iip))/diag(neq,ipg)
ccc$$$     .                          +aw*(rr(iim)-yy(iim))/diag(neq,img))
ccc$$$     .                 +omega01*(an*(rr(jjp)-yy(jjp))/diag(neq,jpg)
ccc$$$     .                          +as*(rr(jjm)-yy(jjm))/diag(neq,jmg))
ccc$$$
ccc$$$                mag = mag + (rr(ii)-yy(ii))**2
ccc$$$              enddo
ccc$$$            enddo
ccc$$$
ccc$$$          endif
cc
cccc        if (omega10.ne.0d0.or.omega01.ne.0d0) then
cccc          write (*,*)'Weighed jacobi not available in coupled Jacobi'
cccc          write (*,*)'Aborting...'
cccc          stop
cccc        endif
cc
cc        call matvec(0,ntot,zz,yy,igrid,bcnd)
cc
cc        do ii = 1,nn/neq
cc
cc          iii = neq*(ii-1)
cc
cc          !Find new residual
cc          do ieq = 1,neq
cc            rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
cc          enddo
cc
cc          !Multiply by D^-1 (stored in diag)
cc          iig = iii + isig - 1
cc          dummy = matmul(diag(:,iig+1:iig+neq),rhs)
cc
cc          !Update zz
cc          do ieq=1,neq
cc            zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
cc          enddo
cc
cc          mag = mag + sum(rhs*rhs)
cc        enddo
cc
ccc     Check convergence
cc
cc        mag = sqrt(mag)
cc
cc        if (itr.eq.1) then
cc          mag0 = mag
cc          if (out.ge.2) write (*,10) itr,mag,mag/mag0
cc        else
cc          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
cc          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
cc        endif
cc
cc      enddo
cc
cc      itr = min(itr,iter)
cc
cc      if (out.ge.2) write (*,*)
cc      if (out.eq.1) write (*,10) itr,mag,mag/mag0
cc
cc      options%iter_out = itr
cc      options%tol_out  = mag/mag0
cc
ccc End program
cc
cc      deallocate (dummy,rhs)
cc
cc      if (fdiag)  deallocate(diag)
cc      call deallocPointers(fpointers)
cc
cc      return
cc
cc 10   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
cc     .        '; Ratio:',1p1e10.2)
cc 20   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
cc     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)
cc
cc      end subroutine jb

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
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax
      real(8)    :: omega0,omega10,omega01,tol
      logical    :: fdiag,fpointers,vbr

      integer(4) :: i,j,k,iv,jv,kv,if,jf,kf,nxf,nyf,nzf,iiv,ivg
      integer(4) :: ii,iii,iib,iiib,iic,iig,iblock,jblock,kblock,nblk
      real(8)    :: mag0,mag1,mag,yy(ntot)

      real(8),allocatable, dimension(:)  :: dummy,rhs

c Begin program

c Read solver configuration

      iter   = options%iter
      omega0 = options%omega
c$$$      omega10= options%omega10
c$$$      omega01= options%omega01
      tol    = options%tol
      fdiag  = options%fdiag
      vbr    = options%vertex_based_relax

      if (out.ge.2) write (*,*)

c Allocate pointers

      call allocPointers(neq,MGgrid,fpointers)

      if (fpointers) then          !JB is NOT called from MG
        igmax = igrid
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

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
          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c Jacobi iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag1 = mag
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
cc          write (*,*) 'JB Loop limits:',imin,imax,jmin,jmax,kmin,kmax

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

cc                write (*,*) 'rhs ',rhs
cc                write (*,*) 'dx  ',dummy
cc                write (*,*) 'zz  ',zz(iii+1:iii+neq)
cc                write (*,*) 'yy  ',yy(iii+1:iii+neq)
cc                write (*,*) 'omega',omega0
cc                write (*,*) 'mag ',mag

              enddo
            enddo
          enddo

        endif

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr-1,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr-1,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.1) then
        call matvec(0,neq,ntot,zz,yy,igrid,bcnd)
        call vecadd(igrid,neq,-1d0,yy,1d0,rr)  !yy=rr-yy
        mag = sqrt(dot(igrid,neq,yy,yy))       !sqrt(yy*yy)
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

 10   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

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

      integer(4) :: iter,alloc_stat,isig,itr,nn,ieq,i,j,k
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

      real(8)    :: omega0,tol
      logical    :: fdiag,fpointers,vbr

      integer(4) :: ic,jc,kc,iv,jv,kv,if,jf,kf,nxf,nyf,nzf,iig
      integer(4) :: ii,iii,iib,iiib,iiv,ivg,iblock,jblock,kblock,nblk
      real(8)    :: mag0,mag1,mag,yy(ntot)

      integer(4) :: ncolors,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3

      real(8),allocatable, dimension(:) :: dummy,rhs

c Begin program

c Read solver configuration

      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag
      ncolors= options%ncolors
      vbr    = options%vertex_based_relax

      if (out.ge.2) write (*,*)

c Set pointers

      call allocPointers(neq,MGgrid,fpointers)

      if (fpointers) then          !JB is NOT called from MG
        igmax = igrid
        vbr_mg = .false.
      endif

      if (fpointers.and.out.ge.2) write (*,*) 'Allocated pointers.'

      isig  = istart(igrid)

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
          call find_mf_diag(neq,nblk,ntot,matvec,igrid,bcnd,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c GS iteration

      allocate(dummy(neq*nblk),rhs(neq*nblk))

      do itr=1,iter

        mag1 = mag
        mag  = 0d0

c       Colored GS 
        if (ncolors > 1) then

          select case(ncolors)
          case(2)
            nrbg1 = 2
            nrbg2 = 1
            nrbg3 = 1
          case(4)
            nrbg1 = 2
            nrbg2 = 2
            nrbg3 = 1
          case(8)
            nrbg1 = 2
            nrbg2 = 2
            nrbg3 = 2
          case default
            write (*,*) 'Unsupported number of colors in GS'
            write (*,*) 'Aborting...'
            stop
          end select

c         VERTEX-BASED RELAXATION
          if (vbr) then

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                call matvec(0,neq,ntot,zz,yy,igrid,bcnd)

cc                !Vertex sampling
cc                do kv=1+mod((    irbg3-1),nrbg3)
cc     .                       ,max(nzf-mg_ratio_z(igrid)+1,1),nrbg3
cc                  do jv=1+mod((kv   +irbg2-2),nrbg2)
cc     .                         ,max(nyf-mg_ratio_y(igrid)+1,1),nrbg2
cc                    do iv=1+mod((jv+kv+irbg1-1),nrbg1)
cc     .                           ,max(nxf-mg_ratio_x(igrid)+1,1),nrbg1
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

            do irbg3 = 1,nrbg3
              do irbg2 = 1,nrbg2
                do irbg1 = 1,nrbg1

                call matvec(0,neq,ntot,zz,yy,igrid,bcnd)

cc                do k=1+mod((    irbg3-1),nrbg3),nzv(igrid),nrbg3
cc                  do j=1+mod((k  +irbg2-2),nrbg2),nyv(igrid),nrbg2
cc                    do i=1+mod((j+k+irbg1-1),nrbg1),nxv(igrid),nrbg1
                do k=1+mod((    irbg3-1),nrbg3),nzv(igrid),nrbg3
                  do j=1+mod((k  +irbg2+irbg1-1),nrbg2),nyv(igrid),nrbg2
                    do i=1+mod((j+k+irbg1-1),nrbg1),nxv(igrid),nrbg1

                      iii = neq*(i-1 + nxv(igrid)*(j-1)
     .                               + nxv(igrid)*nyv(igrid)*(k-1))

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

c       Regular GS
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
cc            write (*,*) 'GS Loop limits:',imin,imax,jmin,jmax,kmin,kmax

            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

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

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr-1,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr-1,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

c Calculate final residual and output info

      itr = min(itr,iter)

      if (out.ge.1) then
        call matvec(0,neq,ntot,zz,yy,igrid,bcnd)
        call vecadd(igrid,neq,-1d0,yy,1d0,rr)  !yy=rr-yy
        mag = sqrt(dot(igrid,neq,yy,yy))       !sqrt(yy*yy)
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

 10   format (' GS Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' GS Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

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

      call allocPointers(neq,MGgrid,fpointers)

      if (fpointers) igmax=grid_params%ngrid-1  !Routine not called from JB,GS,MG

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        write (*,*) 'Error in input of find_mf_mat'
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

c find_mf_diag
c####################################################################
      subroutine find_mf_diag(neq,nblk,ntot,matvec,igrid,bbcnd
     .                       ,diag1)
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
c--------------------------------------------------------------------

      use mg_internal

      implicit none      !For safe fortran

c Call variables

      integer(4) :: neq,nblk,ntot,bbcnd(6,neq),igrid

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

      call allocPointers(neq,MGgrid,fpointers)

      if (fpointers) igmax=grid_params%ngrid-1  !Routine not called from JB,GS,MG

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

        call find_mf_diag_std(neq,ntot,matvec,igrid,bbcnd,diag1)

C     VERTEX-BASED RELAXATION
      else

        do igr = igrid,igmax

          nxf = nxv(igr)
          nyf = nyv(igr)
          nzf = nzv(igr)

          nblkg = nblock(igr)

          nn    = ntotv(igr)

          isig  = istartb(igr)

          size  = neq*nblkg

          x1(1:nn) = 0d0

c       Coarse grid guides diagonal formation process for vertex-based relaxation

          do kv=1,max(nzf-mg_ratio_z(igr)+1,1)
            do jv=1,max(nyf-mg_ratio_y(igr)+1,1)
              do iv=1,max(nxf-mg_ratio_x(igr)+1,1)

              !This defines coarse cell numbers (lexicographic)
              iiv = iv + nxf*(jv-1) + nxf*nyf*(kv-1)
              ivg  = (iiv-1)*neq*nblkg + isig -1

              !Sample block diagonal columns
              do kblock=1,mg_ratio_z(igr)
                do jblock=1,mg_ratio_y(igr)
                  do iblock=1,mg_ratio_x(igr)

                    !This defines block column index
                    icolb = iblock + mg_ratio_x(igr)*(jblock-1)
     .                             + mg_ratio_x(igr)
     .                              *mg_ratio_y(igr)*(kblock-1)

                    !This defines fine cell position (lexicographic)
                    if = iv-1 + iblock
                    jf = jv-1 + jblock
                    kf = kv-1 + kblock

                    ii  = if + nxf*(jf-1) + nxf*nyf*(kf-1)

                    !Sample block diagonal rows
                    do kblock2=1,mg_ratio_z(igr)
                      do jblock2=1,mg_ratio_y(igr)
                        do iblock2=1,mg_ratio_x(igr)
 
                          !This defines block row index
                          irowb =iblock2 + mg_ratio_x(igr)*(jblock2-1)
     .                                   + mg_ratio_x(igr)
     .                                    *mg_ratio_y(igr)*(kblock2-1)

                          !This defines fine cell position (lexicographic)
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

                            !Set column vector corresponding to fine grid node ii and equation ieq
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

      call allocPointers(neq,MGgrid,fpointers)

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

      real(8)    :: dd1,dd2,error
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

      error = 0d0

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
              error = error + dd1
            endif

          enddo

        enddo
      enddo

      write (*,20) error

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



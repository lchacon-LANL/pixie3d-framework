c TODO list:
c
c 2) In weighed Jacobi, replace omega factors by a matrix-vector
c    operation.
c
c 3) Document i/o variables in all subroutines consistently.

c module mg_internal
c######################################################################
      module mg_internal

        use grid

        integer(4) :: ngrdx,ngrdy,ngrdz,ngrid,igmax

        real(8), allocatable, dimension(:,:) :: diag

        integer(4),dimension(:),allocatable ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,nxv,nyv,nzv
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z

      contains

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,fpointers)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer(4)          :: neq
      logical,intent(OUT) :: fpointers

c     Local variables

      integer(4) :: i,alloc_stat

c     Begin program

      fpointers = .false.

      ngrdx = grid_params%ngrdx
      ngrdy = grid_params%ngrdy
      ngrdz = grid_params%ngrdz
      ngrid = grid_params%ngrid

      allocate(istart (ngrid),ntotv (ngrid)
     .        ,istartp(ngrid),ntotvp(ngrid)
     .        ,nxv(ngrid),nyv(ngrid),nzv(ngrid)
     .        ,mg_ratio_x(ngrid),mg_ratio_y(ngrid),mg_ratio_z(ngrid)
     .        ,STAT = alloc_stat)

      !Successful memory allocation
      if (alloc_stat == 0) then

        mg_ratio_x = grid_params%mg_ratio_x
        mg_ratio_y = grid_params%mg_ratio_y
        mg_ratio_z = grid_params%mg_ratio_z

        nxv = grid_params%nxv
        nyv = grid_params%nyv
        nzv = grid_params%nzv

        istartp = grid_params%istartp

        istart (1) = 1
        ntotvp (1) = nxv(1)*nyv(1)*nzv(1)
        ntotv  (1) = neq*ntotvp(1)
        do i = 2,ngrid
          istart (i) = istart (i-1) + ntotv (i-1)
          ntotvp (i) = nxv(i)*nyv(i)*nzv(i)
          ntotv  (i) = neq*ntotvp(i)
        enddo

        fpointers = .true.

      endif

c     End program

      end subroutine allocPointers

c     deallocPointers
c     #################################################################
      subroutine deallocPointers(fpointers)

c     -----------------------------------------------------------------
c     Initializes MG and creates 2D uniform grid
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers

      if (fpointers) deallocate(istart,ntotv,istartp,ntotvp,nxv,nyv,nzv
     .                         ,mg_ratio_x,mg_ratio_y,mg_ratio_z)

      end subroutine deallocPointers

c     getMGvcomp
c     #################################################################
      function getMGvcomp(i,j,k,igr,ieq,neq) result(ijkg)

c     -----------------------------------------------------------------
c     Gives MG vector component corresponing to coordinates (i,j,k) on
c     grid igx,igy,igz.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igr,ieq,neq,ijkg

c     Local variables

        integer(4) :: nx,ny,nz
        logical    :: fpointers

c     Begin program

        call allocPointers(neq,fpointers)

        nx = nxv(igr)
        ny = nyv(igr)
        nz = nzv(igr)

        ijkg = neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
     .     + istart(igr) - 1

        call deallocPointers(fpointers)

      end function getMGvcomp

c     mapArrayToMGVector
c     ##################################################################
      subroutine mapArrayToMGVector(neq,nx,ny,nz,array,mgvector,igr)
c     ------------------------------------------------------------------
c     Maps array into a MG vector excluding ghost nodes.
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
              ii = neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
     .           + istart(igr) - 1
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
     .                             ,array,igr)
c     ------------------------------------------------------------------
c     Maps a vector into an array, filling ghost cells.
c     If vector is a MG vector, set igrid to grid level. Otherwise,
c     set igrid=1.
c     ------------------------------------------------------------------

      implicit none    !For safe fortran

c     Call variables

      integer(4) :: neq,igr,nx,ny,nz,gpos
      real(8)    :: mgvector(*),array(0:nx+1,0:ny+1,0:nz+1,neq)

c     Local variables

      integer(4) :: ieq,i,j,k,ii,offset
      integer(4) :: imin ,imax ,jmin ,jmax ,kmin ,kmax
     .             ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

      call limits(gpos,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)

      offset = 1

      iimin = max(imin-offset,1)
      iimax = min(imax+offset,nx)
      jjmin = max(jmin-offset,1)
      jjmax = min(jmax+offset,ny)
      kkmin = max(kmin-offset,1)
      kkmax = min(kmax+offset,nz)

      do k = kkmin,kkmax
        do j = jjmin,jjmax
          do i = iimin,iimax
            do ieq=1,neq
              ii = neq*(i-1 + nx*(j-1) + nx*ny*(k-1)) + ieq
     .           + istart(igr) - 1
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
      subroutine limits(elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax)
      implicit none
c     ---------------------------------------------------------------
c     Finds limits on loops for matvec routines. Used in finding 
c     diagonal from matvec.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax

c     Local variables

      integer(4) :: el1

c     Begin program

      if (elem.eq.0) then
        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
        kmin = 1
        kmax = nz
      else
        el1  = mod(elem,nx*ny)
        if (el1 == 0) el1 = nx*ny

        imin = mod(el1 ,nx)
        if (imin == 0) imin = nx
        imax = imin

        jmin = 1 + (el1 - imin)/nx
        jmax = jmin

        kmin = 1 + (elem - imin - nx*(jmin-1))/nx*ny
        kmax = kmin
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
c     Solves block systems using a direct solve approach. Requires
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

      use mlsolverSetup

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: iter,igridmin,vcyc,crsedpth
      integer(4) :: orderres,orderprol,alloc_stat

      integer(4) :: guess2,outc,mu
      integer(4) :: izero,i,j,ii,ivcyc

      real(8)    :: xx(2*ntot),yy(2*ntot),wrk(2*ntot)
      real(8)    :: rr0,rr1,mag,mag1,mgtol
      real(8)    :: dummy(ntot),rr(ntot)

      logical    :: fdiag,fpointers,volf

c Begin program

      igridmin = options%igridmin
      vcyc     = options%vcyc
      mgtol    = options%tol
      orderres = options%orderres
      orderprol= options%orderprol
      fdiag    = options%fdiag
      volf     = options%vol_res
      crsedpth = options%mg_coarse_solver_depth
      mu       = options%mg_mu

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

      if (fdiag) then

        if (allocated(diag)) then !Diagonal is known
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          allocate(diag(neq,2*ntot))
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          allocate(diag(neq,2*ntot))
          call find_mf_diag(neq,ntot,matvec,igrid,bcnd,diag)
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
        rr0 = sqrt(sum(y*y))
      else
        call matvec(0,ntot,x,dummy,igrid,bcnd)
        rr = y - dummy
        rr0 = sqrt(sum(rr*rr))
      endif

      if (rr0.lt.1d-16*ntot) then
        if (out.ge.1) write (*,*) 'Initial solution seems exact in MG'
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

c     Perform Mcycle recursively

        if (igrid.eq.igmax) then
          call smooth (igrid)
          exit
        else
          call mcycle(igrid)
        endif

c     Check MG convergence

        call matvec(0,ntot,xx(istart(igrid)),dummy,igrid,bcnd)

        rr   = y - dummy
        mag  = sqrt(sum(rr*rr))
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

c       Begin program

          igc   = igr+1

          nn    = ntotv(igr)
          nnc   = ntotv(igc)

          isig  = istart(igr)
          isigc = istart(igc)

c       Relax error/solution on grid number igr/igrid (find new xx)

          call smooth(igr)

c       Evaluate residual (ie wrk = yy - A xx = yy - wrk )

          call matvec(0,nn,xx(isig),wrk(isig),igr,bcnd)

          wrk(isig:isig+nn-1) = yy(isig:isig+nn-1) - wrk(isig:isig+nn-1)

c       Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) to a coarser grid

          call crestrict(neq
     .                  ,yy(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .                  ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .                  ,orderres,igr,volf)

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

c       Prolong error (wrk = P * xx_1) to a finer grid

          call cprolong(neq
     .                 ,wrk(isig),ntotv(igr),nxv(igr),nyv(igr),nzv(igr)
     .                 ,xx(isigc),ntotv(igc),nxv(igc),nyv(igc),nzv(igc)
     .                 ,orderprol,igc,bcnd)

c       Update existing error on grid 2 (i.e. xx_2): xx_2 = xx_2 + wrk 

          xx(isig:isig+nn-1) = xx(isig:isig+nn-1) + wrk(isig:isig+nn-1)

c       Relax updated error on igr (i.e. xx_igr)

          call smooth(igr)

        end subroutine mcycle

c       smooth
c       ###################################################################
        recursive subroutine smooth(igr)

          implicit none
          integer :: igr,nn,isig,depth1

          if (outc.ge.1) write (*,*) 'Grid Level',igr

          nn   = ntotv (igr)
          isig = istart(igr)

          depth1 = depth + 1

          call getSolver(neq,nn,yy(isig),xx(isig),matvec,igr,bcnd
     .                  ,guess2,outc,depth1)

        end subroutine smooth

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
        subroutine killmg

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

      !Set grid=1 because xc is NOT a MG vector.
      call mapMGVectorToArray(0,neq,xc,nxc,nyc,nzc,arrayc,1)

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
      logical    :: fpointers

c Extrapolation

      real(8)    :: xxf,yyf,zzf,ff

      real(8)    :: xx(nxc+2),yy(nyc+2),zz(nzc+2)

      integer(4) ::  kx,ky,kz,nx,ny,nz,dim,flg
      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db3val
      external      db3val

c Begin program

      call allocPointers(1,fpointers)

c Define fine grid

      igf = igc - 1

c Injection

      if (order.eq.0) then

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

cc              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)

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

cc        call getMGmap(0,0,0,min(igc,ngrdx),min(igc,ngrdy),min(igc,ngrdz),icg,jcg,kcg)
        call getMGmap(1,1,1,igc,igc,igc,icg,jcg,kcg)
        xx(1:nxc+2) = grid_params%xx(icg-1:icg+nxc)
        yy(1:nyc+2) = grid_params%yy(jcg-1:jcg+nyc)
        zz(1:nzc+2) = grid_params%zz(kcg-1:kcg+nzc)

c     Prepare 3d spline interpolation

        flg = 0
        nx = nxc + 2
        ny = nyc + 2
        nz = nzc + 2
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(tx(nx+kx))
        allocate(ty(ny+ky))
        allocate(tz(nz+kz))
        allocate(work(dim))
        allocate(bcoef(nx,ny,nz))

        call db3ink(xx,nx,yy,ny,zz,nz,arrayc,nx,ny,kx,ky,kz,tx,ty,tz
     .             ,bcoef,work,flg)

c     Interpolate

        do kf = 1,nzf
          do jf = 1,nyf
            do if = 1,nxf

              iif = if + nxf*(jf-1) + nxf*nyf*(kf-1)

cc              call getMGmap(if,jf,kf,min(igf,ngrdx),min(igf,ngrdy)
cc     .                     ,min(igf,ngrdz),ifg,jfg,kfg)
              call getMGmap(if,jf,kf,igf,igf,igf,ifg,jfg,kfg)
              xxf = grid_params%xx(ifg)
              yyf = grid_params%yy(jfg)
              zzf = grid_params%zz(kfg)

              xf(iif) = db3val(xxf,yyf,zzf,0,0,0,tx,ty,tz,nx,ny,nz
     .                        ,kx,ky,kz,bcoef,work)

            enddo
          enddo
        enddo

        deallocate(tx,ty,tz,work,bcoef)

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
 
      integer(4) :: ic,if,jc,jf,kc,kf,iif,iic,igc
     .             ,icg,jcg,kcg,ifg,jfg,kfg

      real(8)    :: volt,vol,arrayf(0:nxf+1,0:nyf+1,0:nzf+1)

      logical    :: fpointers

c Interpolation

      real(8)    :: xxc,yyc,zzc,ff

      real(8)    :: xx(nxf),yy(nyf),zz(nzf)

      integer(4) ::  kx,ky,kz,nx,ny,nz,dim,flg
      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: db2val,db3val
      external      db2val,db3val

c Begin program

      call allocPointers(1,fpointers)

c Agglomeration

      if (order.eq.0) then

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

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

        igc = igf + 1

c     Map vectors into arrays

        call getMGmap(1,1,1,igf,igf,igf,ifg,jfg,kfg)

        xx = grid_params%xx(ifg:ifg+nxf-1)
        yy = grid_params%yy(jfg:jfg+nyf-1)
        zz = grid_params%zz(kfg:kfg+nzf-1)

        call mapMGVectorToArray(0,1,xf,nxf,nyf,nzf,arrayf,1)

c     Renormalize without volume fractions

        if (volf) then
          do kf = 1,nzf
            do jf = 1,nyf
              do if = 1,nxf
                arrayf(if,jf,kf) = arrayf(if,jf,kf)
     .                            /volume(if,jf,kf,igf,igf,igf)
              enddo
            enddo
          enddo
        endif

c     Calculate interpolation

        flg = 0
        nx = nxf
        ny = nyf
        nz = nzf
        kx = min(order+1,nx-1)
        ky = min(order+1,ny-1)
        kz = min(order+1,nz-1)

        dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))

        allocate(tx(nx+kx))
        allocate(ty(ny+ky))
        allocate(tz(nz+kz))
        allocate(work(dim))
        allocate(bcoef(nx,ny,nz))

        if (nx == 1) then
          call db2ink(yy,ny,zz,nz,arrayf(1:nx,1:ny,1:nz),ny,ky,kz,ty,tz
     .               ,bcoef,work,flg)
        elseif (ny == 1) then
          call db2ink(xx,nx,zz,nz,arrayf(1:nx,1:ny,1:nz),nx,kx,kz,tx,tz
     .               ,bcoef,work,flg)
        elseif (nz == 1) then
          call db2ink(xx,nx,yy,ny,arrayf(1:nx,1:ny,1:nz),nx,kx,ky,tx,ty
     .               ,bcoef,work,flg)
        else
          call db3ink(xx,nx,yy,ny,zz,nz,arrayf(1:nx,1:ny,1:nz)
     .               ,nx,ny,kx,ky,kz,tx,ty,tz,bcoef,work,flg)
        endif

        do kc = 1,nzc
          do jc = 1,nyc
            do ic = 1,nxc

              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)

              call getMGmap(ic,jc,kc,igc,igc,igc,icg,jcg,kcg)
              xxc = grid_params%xx(icg)
              yyc = grid_params%yy(jcg)
              zzc = grid_params%zz(kcg)

              if (nx == 1) then
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

        deallocate(tx,ty,tz,work,bcoef)

      endif

c End program

      call deallocPointers(fpointers)

      end subroutine restrict

cc*deck restrict
ccc######################################################################
cc      subroutine restrict(xc,nxc,nyc,nzc,arrayf,nxf,nyf,nzf
cc     .                   ,order,igf,volf)
ccc----------------------------------------------------------------------
ccc     This is a restriction routine for a single quantity, with
ccc     arbitrary order of interpolation.
ccc
ccc     In call sequence, we have:
ccc       * xc (real): vector in coarse grid
ccc       * nxc,nyc,nzc(int): dimensions of coarse grid
ccc       * arrayf (real): array of values in fine grid
ccc       * nxf,nyf,nzf(int): dimensions of fine grid
ccc       * order (int): order of interpolation (0-arbitrary)
ccc           If order = 0, it employs simple injection.
ccc           If order > 0, it employs spline interpolation.
ccc       * igf (int): fine grid level identifier
ccc       * bbcnd (int array): boundary condition info.
ccc       * volf (logical): whether vectors contain volume fractions.
ccc----------------------------------------------------------------------
cc
cc      use mg_internal
cc
cc      implicit none        !For safe fortran
cc
ccc Call variables
cc
cc      integer(4) :: nxc,nyc,nzc,nxf,nyf,nzf,order,igf
cc
cc      real(8)    :: xc(nxc*nyc*nzc),arrayf(0:nxf+1,0:nyf+1,0:nzf+1)
cc
cc      logical    :: volf
cc
ccc Local variables
cc 
cc      integer(4) :: ic,if,jc,jf,kc,kf,iif,iic,igc
cc     .             ,icg,jcg,kcg,ifg,jfg,kfg
cc
cc      real(8)    :: volt,vol
cc
cc      logical    :: fpointers
cc
ccc Interpolation
cc
cc      real(8)    :: xxc,yyc,zzc,ff
cc
cc      real(8)    :: xx(nxf+2),yy(nyf+2),zz(nzf+2)
cc
cc      integer(4) ::  kx,ky,kz,nx,ny,nz,dim,flg
cc      real(8), dimension(:),allocatable:: tx,ty,tz,work
cc      real(8), dimension(:,:,:),allocatable:: bcoef
cc
cc      real(8)    :: db3val
cc      external      db3val
cc
ccc Begin program
cc
cc      call allocPointers(1,fpointers)
cc
ccc Agglomeration
cc
cc      if (order.eq.0) then
cc
cc        do kc = 1,nzc
cc          do jc = 1,nyc
cc            do ic = 1,nxc
cc
cc              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)
cc
cc              xc(iic)=0d0
cc              do kf = mg_ratio_z(igf)*(kc-1)+1,mg_ratio_z(igf)*kc
cc                do jf = mg_ratio_y(igf)*(jc-1)+1,mg_ratio_y(igf)*jc
cc                  do if = mg_ratio_x(igf)*(ic-1)+1,mg_ratio_x(igf)*ic
cccc                    iif = if + nxf*(jf-1) + nxf*nyf*(kf-1)
cc
cc                    if (.not.volf) then
cc                      vol = volume(if,jf,kf,min(igf,ngrdx)
cc     .                            ,min(igf,ngrdy),min(igf,ngrdz))
cc                      xc(iic) = xc(iic) + arrayf(if,jf,kf)*vol
cc                      volt = volt + vol
cc                    else
cc                      xc(iic) = xc(iic) + arrayf(if,jf,kf)
cc                    endif
cc
cc                  enddo
cc                enddo
cc              enddo
cc
cc              if (.not.volf) xc(iic) = xc(iic)/volt
cc            enddo
cc          enddo
cc        enddo
cc
ccc Interpolation
cc
cc      else
cc
cc        igc = igf + 1
cc
ccc     Map vectors into arrays
cc
cc        call getMGmap(0,0,0,min(igf,ngrdx),min(igf,ngrdy),min(igf,ngrdz)
cc     .               ,ifg,jfg,kfg)
cc        xx(1:nxf+2) = grid_params%xx(ifg:ifg+nxf+1)
cc        yy(1:nyf+2) = grid_params%yy(jfg:jfg+nyf+1)
cc        zz(1:nzf+2) = grid_params%zz(kfg:kfg+nzf+1)
cc
ccc     Renormalize without volume fractions
cc
cc        if (volf) then
cc          do kf = 0,nzf+1
cc            do jf = 0,nyf+1
cc              do if = 0,nxf+1
cc                arrayf(if,jf,kf) = arrayf(if,jf,kf)
cc     .                        /volume(if,jf,kf,min(igf,ngrdx)
cc     .                               ,min(igf,ngrdy),min(igf,ngrdz))
cc              enddo
cc            enddo
cc          enddo
cc        endif
cc
ccc     Calculate interpolation
cc
cc        flg = 0
cc        kx = order+1
cc        ky = order+1
cc        kz = order+1
cc        nx = nxf + 2
cc        ny = nyf + 2
cc        nz = nzf + 2
cc        dim = nx*ny*nz + max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))
cc
cc        allocate(tx(nx+kx))
cc        allocate(ty(ny+ky))
cc        allocate(tz(nz+kz))
cc        allocate(work(dim))
cc        allocate(bcoef(nx,ny,nz))
cc
cc        call db3ink(xx,nx,yy,ny,zz,nz,arrayf,nx,ny,kx,ky,kz,tx,ty,tz
cc     .             ,bcoef,work,flg)
cc
cc        do kc = 1,nzc
cc          do jc = 1,nyc
cc            do ic = 1,nxc
cc
cc              iic = ic + nxc*(jc-1) + nxc*nyc*(kc-1)
cc
cc              call getMGmap(ic,jc,kc,min(igc,ngrdx),min(igc,ngrdy)
cc     .                     ,min(igc,ngrdz),icg,jcg,kcg)
cc              xxc = grid_params%xx(icg)
cc              yyc = grid_params%yy(jcg)
cc              zzc = grid_params%zz(kcg)
cc
cc              xc(iic) = db3val(xxc,yyc,zzc,0,0,0,tx,ty,tz,nx,ny,nz
cc     .                        ,kx,ky,kz,bcoef,work)
cc
cc              if (volf) then
cc                xc(iic) = xc(iic)
cc     .                   *volume(ic,jc,kc,min(igc,ngrdx)
cc     .                          ,min(igc,ngrdy),min(igc,ngrdz))
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        deallocate(tx,ty,tz,work,bcoef)
cc
cc      endif
cc
ccc End program
cc
cc      call deallocPointers(fpointers)
cc
cc      end subroutine restrict

c jb
c#######################################################################
      recursive subroutine jb(neq,ntot,rr,zz,matvec,options,igrid,bcnd
     .                       ,guess,out,depth)
c--------------------------------------------------------------------
c     Matrix-free Jacobi routine to solve Azz = rr. Call variables:
c       * neq: number of equations (system JB)
c       * ntot: grid dimension
c       * rr,zz: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0 -> no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: iter,alloc_stat,isig
      real(8)    :: omega0,omega10,omega01,tol
      logical    :: fdiag,fpointers

      integer(4) :: i,j,itr,nn,ieq,nx,ny
      integer(4) :: ii,iii,iip,iim,jj,jjp,jjm
      integer(4) :: ig,ipg,img,jg,jpg,jmg,iig
      real(8)    :: mag0,mag1,mag,yy(ntot),delta
      real(8)    :: aw,ae,an,as

      real(8),allocatable, dimension(:)  :: dummy,rhs

c Begin program

c Read solver configuration

      iter   = options%iter
      omega0 = options%omega
      omega10= options%omega10
      omega01= options%omega01
      tol    = options%tol
      fdiag  = options%fdiag

      if (out.ge.2) write (*,*)

c Allocate pointers

      call allocPointers(neq,fpointers)

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

      isig = istart(igrid)

c Find diagonal for smoothers

      if (fdiag) then

cc        allocate(diag(neq,2*ntot),STAT = alloc_stat)

cc        if (alloc_stat.ne.0) then !Diagonal is known (failed alloc)
        if (allocated(diag)) then !Diagonal is known
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          allocate(diag(neq,2*ntot))
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          allocate(diag(neq,2*ntot))
          call find_mf_diag(neq,ntot,matvec,igrid,bcnd,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c Jacobi iteration

      allocate(dummy(neq),rhs(neq))

      do itr=1,iter

        mag1 = mag
        mag  = 0d0

c$$$            nx = nxv(igrid)
c$$$            ny = nyv(igrid)
c$$$
c$$$            do i = 1,nx
c$$$              do j = 1,ny
c$$$                call shiftcoeffs  (i,j,bcnd(1,1))
c$$$
c$$$                call shiftIndices(i,j,ii,iim,iip,jj,jjm,jjp,nx,ny,1
c$$$     .                            ,bcnd(1,1))
c$$$                call shiftIndices(i,j,ig,img,ipg,jg,jmg,jpg,nx,ny,isig
c$$$     .                            ,bcnd(1,1))
c$$$
c$$$                zz(ii) =zz(ii)
c$$$     .                 +omega0*     (rr(ii) -yy(ii) )/diag(neq,ig )
c$$$     .                 +omega10*(ae*(rr(iip)-yy(iip))/diag(neq,ipg)
c$$$     .                          +aw*(rr(iim)-yy(iim))/diag(neq,img))
c$$$     .                 +omega01*(an*(rr(jjp)-yy(jjp))/diag(neq,jpg)
c$$$     .                          +as*(rr(jjm)-yy(jjm))/diag(neq,jmg))
c$$$
c$$$                mag = mag + (rr(ii)-yy(ii))**2
c$$$              enddo
c$$$            enddo
c$$$
c$$$          endif

cc        if (omega10.ne.0d0.or.omega01.ne.0d0) then
cc          write (*,*)'Weighed jacobi not available in coupled Jacobi'
cc          write (*,*)'Aborting...'
cc          stop
cc        endif

        call matvec(0,ntot,zz,yy,igrid,bcnd)

        do ii = 1,nn/neq

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

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

      itr = min(itr,iter)

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) itr,mag,mag/mag0

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
c       * ntot: grid dimension
c       * rr,zz: rhs, solution vectors
c       * matvec: matrix-free matvec product (external)
c       * options: structure containing solver defs.
c       * igrid: grid level in MG applications
c       * guess: 0 -> no initial guess; 1 -> initial guess provided.
c       * out: convergence info output on screen if out > 1. 
c       * depth: integer specifying solver depth in solver_queue
c                definitions
c--------------------------------------------------------------------

      use mlsolverSetup

      use mg_internal

      implicit none       !For safe fortran

c Call variables

      integer(4) :: neq,ntot,igrid,guess,out,depth,bcnd(6,neq)
      real(8)    :: rr(ntot),zz(ntot)

      type (solver_options) :: options

      external     matvec

c Local variables

      integer(4) :: iter,alloc_stat,isig
      real(8)    :: omega0,tol
      logical    :: fdiag,fpointers,fbcnd

      integer(4) :: i,j,k,itr,nn,ieq
      integer(4) :: ii,iii,iig,dum
      integer(4) :: ncolors,irbg1,irbg2,irbg3,nrbg1,nrbg2,nrbg3
      real(8)    :: mag0,mag1,mag,yy(ntot)

      real(8),allocatable, dimension(:) :: dummy,rhs

c Begin program

c Read solver configuration

      iter   = options%iter
      omega0 = options%omega
      tol    = options%tol
      fdiag  = options%fdiag
      ncolors= options%ncolors

      if (out.ge.2) write (*,*)

c Set pointers

      call allocPointers(neq,fpointers)

      if (fpointers.and.out.ge.2) write (*,*) 'Allocating pointers...'

      isig = istart(igrid)

c Find diagonal for smoothers

      if (fdiag) then

cc        allocate(diag(neq,2*ntot),STAT = alloc_stat)

cc        if (alloc_stat.ne.0) then !Diagonal is known (failed alloc)
        if (allocated(diag)) then !Diagonal is known
          if (out.ge.2) write (*,*) 'Diagonal already allocated'
          fdiag = .false.
        elseif (associated(options%diag)) then !Diagonal provided externally
          if (out.ge.2) write (*,*) 'Diagonal externally provided'
          allocate(diag(neq,2*ntot))
          diag = options%diag
        else                      !Form diagonal
          if (out.ge.2) write (*,*) 'Forming diagonal...'
          allocate(diag(neq,2*ntot))
          call find_mf_diag(neq,ntot,matvec,igrid,bcnd,diag)
          if (out.ge.2) write (*,*) 'Finished!'
        endif

      endif

c Preparation for iteration

      nn = ntot

      if (guess.eq.0) zz = 0d0

c GS iteration

      allocate(dummy(neq),rhs(neq))

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

c         Colored grid loops

          do irbg1 = 1,nrbg1
            do irbg2 = 1,nrbg2
              do irbg3 = 1,nrbg3

                call matvec(0,ntot,zz,yy,igrid,bcnd)

                do k=1+mod((    irbg3-1),nrbg3),nzv(igrid),nrbg3
                  do i=1+mod((k  +irbg2-2),nrbg2),nxv(igrid),nrbg2
                    do j=1+mod((i+k+irbg1-1),nrbg1),nyv(igrid),nrbg1

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

c       Regular GS
        else

c         Forward pass
          do ii = 1,nn/neq

            call matvec(-ii,ntot,zz,yy,igrid,bcnd)

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

ccc         Backward pass
cc          do ii = nn/neq,1,-1
cc
cc            call matvec(-ii,ntot,zz,yy,igrid,bcnd)
cc
cc            iii = neq*(ii-1)
cc
cc            !Find new residual
cc            do ieq = 1,neq
cc              rhs(ieq) = rr(iii+ieq) - yy(iii+ieq)
cc            enddo
cc
cc            !Multiply by D^-1
cc            iig = iii + isig - 1
cc            dummy = matmul(diag(:,iig+1:iig+neq),rhs)
cc
cc            !Update solution zz
cc            do ieq=1,neq
cc              zz(iii+ieq) = zz(iii+ieq) + omega0*dummy(ieq)
cc            enddo
cc
cc            mag = mag + sum(rhs*rhs)
cc
cc          enddo

        endif

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol.or.mag.lt.1d-20*nn) exit
        endif

      enddo

      itr = min(itr,iter)

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) itr,mag,mag/mag0

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

c find_mf_diag
c####################################################################
      subroutine find_mf_diag(neq,ntot,matvec,igrid,bbcnd,diag1)
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
     $             ,delta
      integer(4) :: ii,jj,nn,ig,iig,isig
      integer(4) :: igr,ieq,alloc_stat
      logical    :: fpointers

c Begin program

c Allocate MG pointers

      call allocPointers(neq,fpointers)

c Consistency check

      nn  = ntotv(igrid)

      if (nn /= ntot) then
        write (*,*) 'Error in input of find_mf_mat'
        write (*,*) 'Aborting...'
        stop
      endif

c Form diagonal

      do igr = igrid,grid_params%ngrid

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

            call matvec(ii,nn,x1,dummy,igr,bbcnd)

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
            delta = diag1(1,iig+1)*diag1(2,iig+2)
     .             -diag1(2,iig+1)*diag1(1,iig+2)
            mat(1,1) = diag1(2,iig+2)/delta
            mat(1,2) =-diag1(2,iig+1)/delta
            mat(2,1) =-diag1(1,iig+2)/delta
            mat(2,2) = diag1(1,iig+1)/delta

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

          call matvec(0,nn,x1,dummy,igr,bbcnd)

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

          call matvec(0,nn,x1,dummy,igrid,bcnd)

          call findBaseVector(ii,ieq,neq,nn,x1,0d0)

c       Compare column vector ii with corresponding row vector (intersect in
c       diagonal)

          do jj = ii,nn/neq

            call findBaseVector(jj,ieq,neq,nn,x1,1d0)

            call matvec(ii,nn,x1,dummy2,igrid,bcnd)

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




c module mk_pc
c####################################################################
      module mk_pc

        use mk

        use nk, ONLY: matrixFreeMatVec

        !PC setup
        real(8), allocatable,dimension(:,:)   :: mkdiag

        type(mg_array),target :: gmk_phi,gx0,gxv0,gtjac,gmk_vel

        integer,private :: mk_nsweep,mk_maxvcyc,mk_precpass,mk_debug_it
     .                    ,mk_si_it,mk_ores,mk_oprol,mk_mg_coarse_size
     .                    ,mk_ncolors

        real(8),private :: mk_mgtol,mk_omeg_smth,mk_omega

        character(2),private :: mk_smooth

        logical,private :: mk_gm_smooth
     .                    ,mk_mg_glrkin
     .                    ,mk_gm_coarse
     .                    ,mk_line_relax
     .                    ,mk_zebra_relax
     .                    ,mk_mg_debug
     .                    ,mk_cvrg_test
     .                    ,mk_asm_PC
     .                    ,mk_mf_pc

      contains

c     MK_setupPC
c     ###############################################################
      subroutine MK_setupPC(igrid,phi)

c     ---------------------------------------------------------------
c     Sets up MK preconditioner variables and operators.
c     ---------------------------------------------------------------

cc      use xdraw_io

      implicit none

c     Call variables

      integer :: igrid
      real(8) :: phi(:,:,:)

c     Local variables

      integer :: i,j,k,ig,jg,kg,ieq,igr,ijk,ijkg
     .          ,igx,igy,igz,nx,ny,nz

      integer :: order,alloc_stat,bcnd(6)

      integer :: nn,isig

      real(8), allocatable, dimension(:) :: xpc

      real(8) :: dummy(3),x1,y1

c     Debug

cc      real(8),allocatable,dimension(:,:,:) :: debug2
cc      character(13) :: debugfile

c     Begin program

      call MK_PC_params

      igx = igrid
      igy = igrid
      igz = igrid

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

      ngrid = mk_gr_params%ngrid

      order = 0  !Order of restriction for MG arrays

c     Allocate variables

      call MK_allocPCvar

c     Extract arrays and BC's

      call MK_setup_grid(igrid,phi)

      bcnd = mk_bcs

      where(mk_bcs == NEU) bcnd = EXT
cc      call setASM_PC_BC(igx,igy,igz,bcnd,.true.,IIBC=EXT) !THIS STATES NOT TO COMMUNICATE COEFF.

c     Store phi in all grids

      gmk_phi%grid(igrid)%array(:,:,:,1) = phi

      call restrictMGArray(IMK,1,gmk_phi,mk_bcs,igrid,0
     .                    ,isvec=.false.)

c     Store xold in all grids (This should only be done once in the beginning!)

      gx0%grid(igrid)%array = x0

      do igr=2,mk_gr_params%ngrid
        call coarsen_map(igr,mk_gr_params
     .                  ,gx0%grid(igr-1)%array
     .                  ,gx0%grid(igr  )%array)
      enddo

c     Store xv_k (current vertex positions) in all grids

      do igr=1,mk_gr_params%ngrid
        do k=0,mk_gr_params%nzv(igr)
          do j=0,mk_gr_params%nyv(igr)
            do i=0,mk_gr_params%nxv(igr)
              gxv0%grid(igr)%array(i,j,k,:)
     .             =findPos(i,j,k,igr,igr,igr
     .                     ,mk_gr_params%nxv(igr)
     .                     ,mk_gr_params%nyv(igr)
     .                     ,mk_gr_params%nzv(igr)
     .                     ,gx0    %grid(igr)%array
     .                     ,gmk_phi%grid(igr)%array
     .                     ,1)
            enddo
          enddo
        enddo
      enddo

c     Store grid velocity in all grids

      gmk_vel%grid(igrid)%array = XformToCar(1,gvel,.false.)

      do igr=2,mk_gr_params%ngrid
        call coarsen_map(igr,mk_gr_params
     .                  ,gmk_vel%grid(igr-1)%array
     .                  ,gmk_vel%grid(igr  )%array,bcs=mk_map_bcs)
      enddo

      do igr=1,mk_gr_params%ngrid
        gmk_vel%grid(igr)%array
     .       = XformToCurv(igr,gmk_vel%grid(igr)%array,.true.)
      enddo

c     Store current target jacobian in all grids

      gtjac%grid(igrid)%array(:,:,:,1) = tjac

      call restrictMGArray(IMK,1,gtjac,bcnd,igrid,order,isvec=.false.)

c     Find diagonal

      if (mk_3d) then
        call find_mf_diag_colored(1,ntotdp,mk_3d_mtvc,1,mk_bcs
     .                           ,mkdiag,mk_ncolors)
      else
        call find_mf_diag_colored(1,ntotdp,mk_2d_mtvc,1,mk_bcs
     .                           ,mkdiag,mk_ncolors)
      endif

c diag ****
cc      debugfile = 'debug_p'//trim(int2char(my_rank))//'.bin'
cc      open(unit=110,file=debugfile,form='unformatted'
cc     .    ,status='replace')
cc
cc      k = 1
cc
cc      igrid = 3
cc
cc      nx = mk_gr_params%nxv(igrid)
cc      ny = mk_gr_params%nyv(igrid)
cc      nz = mk_gr_params%nzv(igrid)
cc
cc      allocate(debug2(0:nx+1,0:ny+1,0:nz+1))
cc
cc      debug2 = gmk_phi%grid(igrid)%array(:,:,:,1)
cc      call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,0d0,1d0,0d0,1d0,0,110)
cc
cc      do ieq=1,3
cc        debug2 = gx0%grid(igrid)%array(:,:,:,ieq)
cc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cc     .              ,0d0,1d0,0d0,1d0,1,110)
cc      enddo
cc
cccc      do ieq=1,3
cccc        debug2 = 
cccc        call contour(debug2(0:nx+1,0:ny+1,k),nx+2,ny+2
cccc     .              ,0d0,1d0,0d0,1d0,1,110)
cccc      enddo
cc
cc      deallocate(debug2)
cc      close(110)
cc
cc#if defined(petsc)
cc      call PetscEnd(mpierr)
cc#endif
cc      stop
c diag ****

c     End program

      end subroutine MK_setupPC

c     MK_killPC
c     ###################################################################
      subroutine MK_killPC

c     Deallocate variables

      call MK_deallocPCvar

      end subroutine MK_killPC

c     MK_allocPCvar
c     ###################################################################
      subroutine MK_allocPCvar

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat,ntotd2p

c     Begin program

        ntotd2p = 2*ntotdp

        allocate (mkdiag(1,ntotd2p))
      
!      allocate (pc_res(ntotd2p,1))

        call allocateMGArray(1,gmk_phi)

        call allocateMGArray(1,gtjac)

        call allocateMGArray(3,gx0)

        call allocateMGArray(3,gxv0)

        call allocateMGArray(3,gmk_vel)

c     End program

      end subroutine MK_allocPCvar

c     MK_deallocPCvar
c     ###################################################################
      subroutine MK_deallocPCvar

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

      deallocate (mkdiag)

      call deallocateMGArray(gmk_phi)

      call deallocateMGArray(gx0)

      call deallocateMGArray(gxv0)

      call deallocateMGArray(gtjac)

      call deallocateMGArray(gmk_vel)

c     End program

      end subroutine MK_deallocPCvar

c     MK_PC_params
c     ##################################################################
      subroutine MK_PC_params

c     ------------------------------------------------------------------
c     Define setup variables for PC
c     ------------------------------------------------------------------

c     Begin program

      mk_ncolors     = 8      !Number of colors for diagonal calc., smoothing

      mk_asm_PC      =.false. !Whether to do additive Schwartz PC

      mk_gm_smooth   =.false. !Whether to use GMRES smoothing in MG
      mk_gm_coarse   =.false.  !Whether to use GMRES as the coarse solve
      mk_mg_glrkin   =.false. !Whether to use Galerking MG
      mk_line_relax  =.false. !Whether to perform line relaxation in MG
      mk_zebra_relax =.false. !Whether to perform zebra (line) relaxation in MG. If true,
                              !  smoother defaults to GS

      mk_mg_coarse_size = 2   !Minimum coarse grid size (mg_ratio^mg_coarse_size)
!!      mk_mg_coarse_size = 3   !Minimum coarse grid size (mg_ratio^mg_coarse_size)

      mk_mg_debug    = .false.!Whether to go into debugging mode in MG
!!      mk_mg_debug    = .true.!Whether to go into debugging mode in MG
      mk_cvrg_test   = .false.!Whether to perform convergence test of smoother

!!      mk_maxvcyc     = 10      !Maximum number of MG vcycles
      mk_maxvcyc     = 2      !Maximum number of MG vcycles
      mk_nsweep      = 2      !Maximum number of smoothing passes
      mk_ores        = 0      !Restriction order for MG
      mk_oprol       = 1      !Prolongation order for MG
!!      mk_oprol       = 0      !Prolongation order for MG
      mk_mgtol       = 1d-2  !MG tolerance

      mk_smooth      = 'jb'   !Type of smoother ('jg','gs','gm')
      mk_omeg_smth   = 0.7    !Damping for MG smoother (only when 'jb' is chosen)

      !Consistency checks
      if (np > 1) mk_gm_coarse = .false.       !Do not perform GMRES coarse solve in parallel

      if (mk_gm_smooth) mk_smooth = 'gm'

      if (bcSP()) then
        mk_smooth = 'gs'
        mk_ores  = 1
        mk_oprol = 1
        mk_zebra_relax = .true.
      endif

      if (mk_smooth /= 'jb') mk_omeg_smth = 1d0

      end subroutine MK_PC_params

c     MK_applyPC
c     ###################################################################
      subroutine MK_applyPC(igrid,y,x,iout)

c     -------------------------------------------------------------------
c     This subroutine solves P x = y_mk for the vector x.
c     Vector y contains residuals for MK and application.
c     The parameter iout controls output level:
c       * iout <= 0 --> No output
c       * iout >  0 --> Level of output design by user.
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,iout
      real(8) :: x(:),y(:,:)

c     Local variables

      integer :: i,j,k,ieq,igx,igy,igz,nx,ny,nz,nn,neq

      real(8) :: vol,mag,mag0,mag_n,ratio

      integer :: ii,iii,ntotp,guess

      logical :: nrmtst=.true.,pc_diag

c     Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nx = mk_gr_params%nxv(igx)
      ny = mk_gr_params%nyv(igy)
      nz = mk_gr_params%nzv(igz)

      nn  = size(x)
      neq = size(y,2)

      if (nn /= nx*ny*nz) then
        write (*,*)'Grid sizes do not agree in applyPreconditioner'
        write (*,*)'Aborting...'
        stop
      endif

      x = 0d0

      y = y/(1d0+tau/dt)  !Renormalize MK residual to account for time constant (introduced in evaluateNLF)

c     Solve linear shellcmd with MG

      guess = 0

      if (mk_3d) then
        call cSolver(1,nn,y(:,neq),x
     .            ,mk_bcs,igrid,iout,guess,mk_3d_mtvc,vol_wgt
     .            ,tol           = mk_mgtol
     .            ,dg            = mkdiag
     .            ,mg_debug      = mk_mg_debug
     .            ,mg_galerkin   = mk_mg_glrkin
     .            ,mg_order_res  = mk_ores
     .            ,mg_order_prol = mk_oprol
     .            ,mg_coarse_grid_size = mk_mg_coarse_size
     .            ,mg_gm_coarse_solve  = mk_gm_coarse
     .            ,mg_smooth     = mk_smooth
     .            ,mg_vcyc       = mk_maxvcyc
     .            ,sm_it         = mk_nsweep
     .            ,sm_omega      = mk_omeg_smth
     .            ,sm_ncolors    = mk_ncolors
     .            ,sm_line_relax = mk_line_relax
     .            ,sm_zebra_relax= mk_zebra_relax
     .            )
      else
        call cSolver(1,nn,y(:,neq),x
     .            ,mk_bcs,igrid,iout,guess,mk_2d_mtvc,vol_wgt
     .            ,tol           = mk_mgtol
     .            ,dg            = mkdiag
     .            ,mg_debug      = mk_mg_debug
     .            ,mg_galerkin   = mk_mg_glrkin
     .            ,mg_order_res  = mk_ores
     .            ,mg_order_prol = mk_oprol
     .            ,mg_coarse_grid_size = mk_mg_coarse_size
     .            ,mg_gm_coarse_solve  = mk_gm_coarse
     .            ,mg_smooth     = mk_smooth
     .            ,mg_vcyc       = mk_maxvcyc
     .            ,sm_it         = mk_nsweep
     .            ,sm_omega      = mk_omeg_smth
     .            ,sm_ncolors    = mk_ncolors
     .            ,sm_line_relax = mk_line_relax
     .            ,sm_zebra_relax= mk_zebra_relax
     .            )
      endif

c     This redefines grid quantities, including grid velocity

      call MK_correct_RHS(x,y)

c End program

      end subroutine MK_applyPC

c     MK_correct_RHS
c     ##############################################################
      subroutine MK_correct_RHS(dx_mk,rhs)

c     --------------------------------------------------------------
c     Find grid source to physics equation, rhs_mk=U_G.dx_G
c     --------------------------------------------------------------

        implicit none

c     Call variables

        real(8)    :: dx_mk(:),rhs(:,:)

c     Local variables

        integer(4) :: nn,neq,ntot
        real(8)    :: z (size(rhs))
     .               ,rr(size(rhs))
     .               ,scat_vec(size(rhs,1),size(rhs,2))

        type(grid_mg_def),pointer :: grid_sav => null()

c     Begin program

        nn  = size(rhs,1)
        neq = size(rhs,2)
        ntot = nn*neq

        rr = 0d0

        !Gather MK solution update
        scat_vec = 0d0
        scat_vec(:,neq) = dx_mk(:)   !MK solution update

        call gatherVector(z,scat_vec)

        !Save current map
        call equateGridStructure(grid_sav,gv%gparams)
        grid_params => grid_sav

        !Do linear matvec for application (NOTE: this redefines grid quantities!)
        mk_res_eval = .false.
        call matrixFreeMatVec(ntot,z,rr)
        mk_res_eval = .true.

        !Recover map and associated metrics
        grid_params => gv%gparams
        gmetric => grid_params%gmetric !For backward compatibility
        xnew = gmetric%grid(1)%car

        call deallocateGridStructure(grid_sav)

        !Add contribution to application RHS
        call scatterVector(scat_vec,rr)
        rhs(:,1:neq-1) = rhs(:,1:neq-1) - scat_vec(:,1:neq-1) !Update application residual

      end subroutine MK_correct_RHS

c     mk_3d_mtvc
c     ####################################################################
      subroutine mk_3d_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c     --------------------------------------------------------------------
c     This subroutine is a matvec for the linearized Beltrami's equation.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D,
c              27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,ieq,bcmod(6,neq)
     .             ,igx,igy,igz,nxx,nyy,nzz

      integer    :: isig,ijk,ijkg,ip,im,jp,jm,kp,km,ic,jc
      integer    :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:)   :: dphi
      real(8),pointer    ,dimension(:,:,:,:)   :: phi_k,xv_k,tj

      real(8),allocatable,dimension(:,:,:)     :: pjac
      real(8),allocatable,dimension(:,:,:,:)   :: dxcar

      real(8)    :: upwind,nu2,dd0(3,3),dpert(3,3),srce(3),src
     .             ,xx_1,xx_2,eps,flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,fp(neq),fm(neq),modz,modx,xdotz,lap,pos(2)
     .             ,cov(3),cnv(3),dvol

      logical    :: fpointers,is_cnv

c     Begin program

      is_cnv = .true.

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nxx = mk_gr_params%nxv(igx)
      nyy = mk_gr_params%nyv(igy)
      nzz = mk_gr_params%nzv(igz)

c     Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(dphi(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      dphi = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dphi,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dphi,bcnd
     .            ,icomp=(/IMK/),is_vec=.false.,iorder=0)

c     Aliases

      phi_k => gmk_phi%grid(igrid)%array
      xv_k  => gxv0%grid(igrid)%array
      tj    => gtjac%grid(igrid)%array

c     Find auxiliary quantities

      allocate(dxcar(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,pjac (0:nxx+1,0:nyy+1,0:nzz+1))

      do k = kmin-1,kmax
        do j = jmin-1,jmax
          do i = imin-1,imax
            dxcar(i,j,k,:)=findPos(i,j,k,igx,igy,igz,nxx,nyy,nzz
     .                            ,xv_k,dphi,1,abs=.false.)
          enddo
        enddo
      enddo

c     Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            if (vol_wgt) then
              dvol = gmetric%grid(igx)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            ijk  = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg = ijk + isig - 1

            !Linearized MK elliptic operator piece
            lap = mk_det_3d(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                     ,dxcar(:,:,:,1)
     .                     ,xv_k (:,:,:,2)
     .                     ,xv_k (:,:,:,3))
     .           +mk_det_3d(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                     ,xv_k (:,:,:,1)
     .                     ,dxcar(:,:,:,2)
     .                     ,xv_k (:,:,:,3))
     .           +mk_det_3d(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                     ,xv_k (:,:,:,1)
     .                     ,xv_k (:,:,:,2)
     .                     ,dxcar(:,:,:,3))

            lap = lap
     .           +mk_lap     (i,j,k,nxx,nyy,nzz,igx,igy,igz,dphi
     .                       ,enf_bc=.true.)
     .           -mk_lap_vrtx(i,j,k,nxx,nyy,nzz,igx,igy,igz,dphi
     .                       ,enf_bc=.true.)

            !Advective piece due to nonlinear RHS term
            cov = grad(i,j,k,nxx,nyy,nzz,igx,igy,igz,tj(:,:,:,1))
            cnv =-matmul(mk_gmetric%grid(igx)%gsup(i,j,k,:,:),cov)
     .           /mk_gmetric%grid(igx)%jac(i,j,k)

            if (cnv(1) > 0d0) then
              cov(1) = (dphi(i,j,k,1)-dphi(i-1,j,k,1))
     .                 /mk_gr_params%dx(ig-1)
            else
              cov(1) = (dphi(i+1,j,k,1)-dphi(i,j,k,1))
     .                 /mk_gr_params%dx(ig)
            endif
            if (cnv(2) > 0d0) then
              cov(2) = (dphi(i,j,k,1)-dphi(i,j-1,k,1))
     .                 /mk_gr_params%dy(jg-1)
            else
              cov(2) = (dphi(i,j+1,k,1)-dphi(i,j,k,1))
     .                 /mk_gr_params%dy(jg)
            endif
            if (cnv(3) > 0d0) then
              cov(3) = (dphi(i,j,k,1)-dphi(i,j,k-1,1))
     .                 /mk_gr_params%dz(kg-1)
            else
              cov(3) = (dphi(i,j,k+1,1)-dphi(i,j,k,1))
     .                 /mk_gr_params%dz(kg)
            endif

            lap = lap + dot_product(cnv,cov)

            y(ijk) =-dvol*lap
          enddo
        enddo
      enddo

c     End program

      deallocate(dphi,dxcar,pjac)

      nullify(phi_k,xv_k)

      call deallocPointers(fpointers)

      end subroutine mk_3d_mtvc

c     mk_2d_mtvc
c     ####################################################################
      subroutine mk_2d_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c     --------------------------------------------------------------------
c     This subroutine is a matvec for the linearized Beltrami's equation.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D,
c              27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,ieq,bcmod(6,neq)
     .             ,igx,igy,igz,nxx,nyy,nzz

      integer    :: isig,ijk,ijkg,ip,im,jp,jm,kp,km,ic,jc
      integer    :: imin,imax,jmin,jmax,kmin,kmax

      real(8),allocatable,dimension(:,:,:,:)   :: dphi
      real(8),pointer    ,dimension(:,:,:,:)   :: phi_k,xv_k,tj

      real(8),allocatable,dimension(:,:,:)     :: pjac
      real(8),allocatable,dimension(:,:,:,:)   :: dxcar

      real(8)    :: upwind,nu2,dd0(3,3),dpert(3,3),srce(3),src
     .             ,xx_1,xx_2,eps,flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,fp(neq),fm(neq),modz,modx,xdotz,lap,pos(2)
     .             ,cov(3),cnv(3),dvol

      logical    :: fpointers,is_cnv

c     Begin program

      is_cnv = .true.

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

      nxx = mk_gr_params%nxv(igx)
      nyy = mk_gr_params%nyv(igy)
      nzz = mk_gr_params%nzv(igz)

c     Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz,igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(dphi(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      dphi = 0d0

      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,dphi,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,dphi,bcnd
     .            ,icomp=(/IMK/),is_vec=.false.,iorder=0)

c     Aliases

      phi_k => gmk_phi%grid(igrid)%array
      xv_k  => gxv0%grid(igrid)%array
      tj    => gtjac%grid(igrid)%array

c     Find auxiliary quantities

      allocate(dxcar(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .        ,pjac (0:nxx+1,0:nyy+1,0:nzz+1))

      do k = kmin-1,kmax
        do j = jmin-1,jmax
          do i = imin-1,imax
            dxcar(i,j,k,:)=findPos(i,j,k,igx,igy,igz,nxx,nyy,nzz
     .                            ,xv_k,dphi,1,abs=.false.)
          enddo
        enddo
      enddo

c     Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            if (vol_wgt) then
              dvol = gmetric%grid(igx)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            ijk  = i + nxx*(j-1) + nxx*nyy*(k-1)

            ijkg = ijk + isig - 1

            !Linearized MK elliptic operator piece
            lap =mk_det_2d(i,j,k,nxx,nyy,nzz,igx,igy,igz,dxcar(:,:,:,1)
     .                                                  ,xv_k (:,:,:,2))
     .          +mk_det_2d(i,j,k,nxx,nyy,nzz,igx,igy,igz,xv_k (:,:,:,1)
     .                                                  ,dxcar(:,:,:,2))

            lap = lap + mk_lap(i,j,k,nxx,nyy,nzz,igx,igy,igz,dphi
     .                        ,enf_bc=.true.)

            !Advective piece due to nonlinear RHS term
            cov = grad(i,j,k,nxx,nyy,nzz,igx,igy,igz,tj(:,:,:,1))
            cnv =-matmul(mk_gmetric%grid(igx)%gsup(i,j,k,:,:),cov)
     .           /mk_gmetric%grid(igx)%jac(i,j,k)

            if (cnv(1) > 0d0) then
              cov(1) = (dphi(i,j,k,1)-dphi(i-1,j,k,1))
     .                 /mk_gr_params%dx(ig-1)
            else
              cov(1) = (dphi(i+1,j,k,1)-dphi(i,j,k,1))
     .                 /mk_gr_params%dx(ig)
            endif
            if (cnv(2) > 0d0) then
              cov(2) = (dphi(i,j,k,1)-dphi(i,j-1,k,1))
     .                 /mk_gr_params%dy(jg-1)
            else
              cov(2) = (dphi(i,j+1,k,1)-dphi(i,j,k,1))
     .                 /mk_gr_params%dy(jg)
            endif
            if (cnv(3) > 0d0) then
              cov(3) = (dphi(i,j,k,1)-dphi(i,j,k-1,1))
     .                 /mk_gr_params%dz(kg-1)
            else
              cov(3) = (dphi(i,j,k+1,1)-dphi(i,j,k,1))
     .                 /mk_gr_params%dz(kg)
            endif

            lap = lap + dot_product(cnv,cov)

            y(ijk) =-dvol*lap
          enddo
        enddo
      enddo

c     End program

      deallocate(dphi,dxcar,pjac)

      nullify(phi_k,xv_k)

      call deallocPointers(fpointers)

      end subroutine mk_2d_mtvc

      end module mk_pc

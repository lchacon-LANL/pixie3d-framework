#if !defined(PER_SYNC_ONLY)
#define PER_SYNC_ONLY .false.
#endif

c module BCS_variables
c####################################################################
      module BCS_variables

        use oned_int

        use grid_mpi

        use grid_operations

        use math

        use bc_ops

        integer :: nxbc,nybc,nzbc,igxbc,igybc,igzbc

        integer :: bc_debug_unit
        logical :: bc_debug=.false.

        integer,parameter :: VECTOR=1,SCALAR=0,AUX=-1

        real(8),pointer,private,dimension(:,:,:)   :: zeros
        real(8),pointer,private,dimension(:,:,:,:) :: vzeros

      contains

c     set_BC_debug
c     ##############################################################
      subroutine set_BC_debug(bcdebug)

      use io

      implicit none

c     --------------------------------------------------------------
c     Whether to skip corners in BC treatment
c     --------------------------------------------------------------

c     Call variables

      logical :: bcdebug

c     Begin program

      bc_debug = bcdebug

      if (bc_debug) then
        bc_debug_unit=find_unit(2222+my_rank)
        open(unit=bc_debug_unit
     .      ,file='bc_debug_p'//trim(int2char(my_rank))//'.txt'
     .      ,status='unknown')
      endif

c     End program

      end subroutine set_BC_debug

      end module BCS_variables

c module singularBCinterface
c####################################################################
      module singularBCinterface

        use BCS_variables

      contains

c     scalarSingularBC
c     #################################################################
      subroutine scalarSingularBC(gl_def,array)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * array: contains variable on which singular BC is imposed
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gl_def

      real(8) :: array(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer :: i,j,k,ig,jg,kg,jj,jlog,nyg,jglobal,mm

#if defined(petsc)
      real(8),allocatable,dimension(:)   :: local_x,global_x
      real(8),allocatable,dimension(:,:) :: garray
#endif

c     Begin program

c     Check that we are at domain adjacent to singular point

      i = 1

      if (.not.isSP2(gl_def,i,igxbc)) return

c     Check that domain size is a power of 2

      nyg = gl_def%nygl(igybc)

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('scalarSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

c     Parallel Allgather in SP communicator

#if defined(petsc)
      allocate(garray(nyg,0:nzbc+1))

      if (gl_def%my < 1) then
        call pstop('vectorSingularBC','npy < 1')
      elseif (gl_def%my == 1) then
        garray = array(i,1:nyg,:)
      else
        allocate(global_x(nyg*(nzbc+2)))
        call MPI_Allgather(array(i,1:nybc,0:nzbc+1)
     .                             ,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                    ,global_x,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                    ,gl_def%MPI_COMM_SP,mpierr)
        do k=1,gl_def%my
          garray(1+nybc*(k-1):nybc*k,:)
     .       = reshape(global_x(1+nybc*(nzbc+2)*(k-1):nybc*(nzbc+2)*k)
     .                , (/ nybc,nzbc+2 /))
        enddo
        deallocate(global_x)
      endif
#endif

c     Element exchange in poloidal direction

      jlog = gl_def%jlo(igybc)

      mm = 1
      if (coords=='hel'.or.coords=='hl2') mm = gl_def%params(1)

      do k=0,nzbc+1
        do j=0,nybc+1
          jglobal = j + jlog - 1
          jj = mod(jglobal+mm*nyg/2,nyg)
          if (jj == 0) jj = nyg
#if defined(petsc)
          array(i-1,j,k) = garray(jj,k)
#else
          array(i-1,j,k) = array(i,jj,k)
#endif
        enddo
      enddo

c     End program

#if defined(petsc)
      deallocate(garray)
#endif

      end subroutine scalarSingularBC

c     vectorSingularBC
c     #################################################################
      subroutine vectorSingularBC(gl_def,ivar,vec,cov)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gl_def

      integer :: ivar
      real(8) :: vec(0:nxbc+1,0:nybc+1,0:nzbc+1)
      logical :: cov

c     Local variables

      integer :: i,j,k,jj,jlog,nyg,jglobal,mm

#if defined(petsc)
      real(8),allocatable,dimension(:)   :: local_x,global_x
      real(8),allocatable,dimension(:,:) :: garray
#endif

c     Begin program

c     Check that we are at domain adjacent to singular point

      i = 1

      if (.not.isSP2(gl_def,i,igxbc)) return

c     Check that domain size is a power of 2

      nyg = gl_def%nygl(igybc)

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('vectorSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

c     Parallel Allgather in SP communicator

#if defined(petsc)
      allocate(garray(nyg,0:nzbc+1))

      if (gl_def%my < 1) then
        call pstop('vectorSingularBC','npy < 1')
      elseif (gl_def%my == 1) then
        garray = vec(i,1:nyg,:)
      else
        allocate(global_x(nyg*(nzbc+2)))
        call MPI_Allgather(vec(i,1:nybc,0:nzbc+1)
     .                             ,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                    ,global_x,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                    ,gl_def%MPI_COMM_SP,mpierr)
        do k=1,gl_def%my
          garray(1+nybc*(k-1):nybc*k,:)
     .       = reshape(global_x(1+nybc*(nzbc+2)*(k-1):nybc*(nzbc+2)*k)
     .                , (/ nybc,nzbc+2 /))
        enddo
        deallocate(global_x)
      endif
#endif

c     Element exchange in poloidal direction

      mm = 1
      if (coords=='hel'.or.coords=='hl2') mm = gl_def%params(1)

      jlog = gl_def%jlo(igybc)
      if (coords == 'hl2') then
        do k=0,nzbc+1
          do j=0,nybc+1
            jglobal = j + jlog - 1
            jj = mod(jglobal+mm*nyg/2,nyg)
            if (jj == 0) jj = nyg

#if defined(petsc)
            vec(i-1,j,k) =  garray(jj,k)
#else
            vec(i-1,j,k) =  vec(i,jj,k)
#endif
          enddo
        enddo
      else
        do k=0,nzbc+1
          do j=0,nybc+1
            jglobal = j + jlog - 1
            jj = mod(jglobal+mm*nyg/2,nyg)
            if (jj == 0) jj = nyg

            if (cov) then
              select case(ivar)
              case(1)
#if defined(petsc)
                vec(i-1,j,k) = -garray(jj,k)
#else
                vec(i-1,j,k) = -vec(i,jj,k)
#endif
              case(2,3)
#if defined(petsc)
                vec(i-1,j,k) =  garray(jj,k)
#else
                vec(i-1,j,k) =  vec(i,jj,k)
#endif
              end select
            else
              select case(ivar)
              case(1)
#if defined(petsc)
                vec(i-1,j,k) =  garray(jj,k)
#else
                vec(i-1,j,k) =  vec(i,jj,k)
#endif
              case(2,3)
#if defined(petsc)
                vec(i-1,j,k) = -garray(jj,k)
#else
                vec(i-1,j,k) = -vec(i,jj,k)
#endif
              end select
            endif
          enddo
        enddo
      endif

c     End program

#if defined(petsc)
      deallocate(garray)
#endif

      end subroutine vectorSingularBC

c     vectorSingularBC_stg
c     #################################################################
      subroutine vectorSingularBC_stg(gl_def,ivar,vec,cov)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gl_def

      integer :: ivar
      real(8) :: vec(0:nxbc+1,0:nybc+1,0:nzbc+1)
      logical :: cov

c     Local variables

      integer :: i,j,k,jj,jlog,nyg,jglobal,mm

#if defined(petsc)
      real(8),allocatable,dimension(:)   :: local_x,global_x
      real(8),allocatable,dimension(:,:) :: garray
#endif

c     Begin program

c     Check that we are at domain adjacent to singular point

      i = 1

      if (.not.isSP2(gl_def,i,igxbc)) return

c     Do staggered radial component case first (@origin v^1= r*vr = 0)

      if (ivar == 1.and.(.not.cov)) then
        vec(i-1,:,:) = 0d0
        return
      endif

c     Check that domain size is a power of 2

      nyg = gl_def%nygl(igybc)

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('vectorSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

c     Parallel Allgather in SP communicator

#if defined(petsc)
      allocate(garray(nyg,0:nzbc+1))

      if (gl_def%my < 1) then
        call pstop('vectorSingularBC','npy < 1')
      elseif (gl_def%my == 1) then
        garray = vec(i,1:nyg,:)
      else
        allocate(global_x(nyg*(nzbc+2)))
        call MPI_Allgather(vec(i,1:nybc,0:nzbc+1)
     .                             ,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                    ,global_x,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                    ,gl_def%MPI_COMM_SP,mpierr)
        do k=1,gl_def%my
          garray(1+nybc*(k-1):nybc*k,:)
     .       = reshape(global_x(1+nybc*(nzbc+2)*(k-1):nybc*(nzbc+2)*k)
     .                , (/ nybc,nzbc+2 /))
        enddo
        deallocate(global_x)
      endif
#endif

c     Element exchange in poloidal direction

      mm = 1
      if (coords=='hel'.or.coords=='hl2') mm = gl_def%params(1)

      jlog = gl_def%jlo(igybc)
      if (coords == 'hl2') then
        stop 'hl2 staggered option not implemented for hl2 geom'
c$$$        do k=0,nzbc+1
c$$$          do j=0,nybc+1
c$$$            jglobal = j + jlog - 1
c$$$            jj = mod(jglobal+nyg/2,nyg)
c$$$            if (jj == 0) jj = nyg
c$$$
c$$$#if defined(petsc)
c$$$            vec(i-1,j,k) =  garray(jj,k)
c$$$#else
c$$$            vec(i-1,j,k) =  vec(i,jj,k)
c$$$#endif
c$$$          enddo
c$$$        enddo
      else
        do k=0,nzbc+1
          do j=0,nybc+1
            jglobal = j + jlog - 1
            jj = mod(jglobal+mm*nyg/2,nyg)
            if (jj == 0) jj = nyg

            if (cov) then
              select case(ivar)
              case(1)
#if defined(petsc)
                vec(i-1,j,k) = 0.5*(vec(i,j,k)-garray(jj,k)) !Average r+ and r- to origin
#else
                vec(i-1,j,k) = 0.5*(vec(i,j,k) -vec(i,jj,k)) !Average r+ and r- to origin
#endif
              case(2,3)
#if defined(petsc)
                vec(i-1,j,k) =  garray(jj,k)
#else
                vec(i-1,j,k) =  vec(i,jj,k)
#endif
              end select
            else
              select case(ivar)
              case(2,3)
#if defined(petsc)
                vec(i-1,j,k) = -garray(jj,k)
#else
                vec(i-1,j,k) = -vec(i,jj,k)
#endif
              end select
            endif
          enddo
        enddo
      endif

c     End program

#if defined(petsc)
      deallocate(garray)
#endif

      end subroutine vectorSingularBC_stg

c     averageSingularBC
c     ################################################################
      subroutine averageSingularBC(gl_def,igr,vec,order_avg,order_ext
     .                            ,iscnv,iscar)

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gl_def

      real(8) :: vec(0:,0:,0:,:)

      integer :: igr,order_avg,order_ext

      logical :: iscnv,iscar

c     Local variables
      
      integer :: nx,ny,nz,neq,i,j,k,ig,jg,kg,ieq,one,two,thr
     .          ,oavg,oext,icomp,ivar,nr

      real(8) :: curv(3),car(3)

      real(8),allocatable :: sp_sum(:,:),dv_sum(:,:),l_sum(:,:)

      real(8),pointer,dimension(:) :: xx

c     Begin program

      neq = size(vec,4)

      nx = gl_def%nxv(igr)
      ny = gl_def%nyv(igr)
      nz = gl_def%nzv(igr)

c     Prepare interpolations

      oext = min(order_ext,nx)
      oavg = min(order_avg,nx)
      nr = max(oext,oavg)

      one = 1 ; two = 2; thr = 3

      if (nr <= 1) then
        two = one
        thr = one
      elseif (nr == 2) then
        thr = two
      endif

      call getMGmap(gl_def,1,1,1,igr,igr,igr,ig,jg,kg)

      xx => gl_def%xx

c     Regularize innermost cell layer (closest to SP)

      allocate(sp_sum(0:nz+1,3),dv_sum(0:nz+1,3))

      do ivar=1,neq/3

        !Transform to Cartesian
        if (.not.iscar) then
          icomp = 3*(ivar-1)

          do k = 0,nz+1
            do j = 0,ny+1
              do i = 1,nr+1
                curv = vec(i,j,k,icomp+1:icomp+3)
                vec(i,j,k,icomp+1:icomp+3)
     .             = XformToCar_ijk(gl_def,i,j,k,igr,curv,.not.iscnv)
              enddo
            enddo
          enddo
        endif

        !Compute weighed average
        do ieq=1,3

          icomp = ieq+3*(ivar-1)

cc          vec(0,:,:,icomp) = vec(1,:,:,icomp)
          vec(0,:,:,icomp) =
     .        quad_int(xx(ig),xx(ig+one),xx(ig+two),xx(ig+thr)
     .                ,vec(1    ,:,:,icomp),vec(1+one,:,:,icomp)
     .                ,vec(1+two,:,:,icomp),vec(1+thr,:,:,icomp)
     .                ,0d0,oavg)

          do k=0,nz+1
            sp_sum(k,ieq)= sum(vec(0,1:ny,k,icomp)
     .                *gl_def%gmetric%grid(igr)%dvol(1,1:ny,k))
            dv_sum(k,ieq)=
     .             sum(gl_def%gmetric%grid(igr)%dvol(1,1:ny,k))
          enddo
        enddo

#if defined(petsc)
        allocate(l_sum(0:nz+1,3))
        l_sum = sp_sum
        call MPI_Allreduce(l_sum,sp_sum,3*(nz+2)
     .                    ,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gl_def%MPI_COMM_SP,mpierr)

        l_sum = dv_sum
        call MPI_Allreduce(l_sum,dv_sum,3*(nz+2)
     .                    ,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gl_def%MPI_COMM_SP,mpierr)

        deallocate(l_sum)
#endif

        !Assign back to innermost cell center
        do ieq=1,3
          icomp = ieq+3*(ivar-1)

          do k=0,nz+1
            vec(0,:,k,icomp) = sp_sum(k,ieq)/dv_sum(k,ieq)
          enddo

          vec(1,:,:,icomp) =
     .        quad_int(0d0,xx(ig+one),xx(ig+two),xx(ig+thr)
     .                ,vec(0    ,:,:,icomp),vec(1+one,:,:,icomp)
     .                ,vec(1+two,:,:,icomp),vec(1+thr,:,:,icomp)
     .                ,xx(ig),oext)
        enddo

        !Transform back to Curvilinear
        if (.not.iscar) then
          icomp = 3*(ivar-1)

          do k = 0,nz+1
            do j = 0,ny+1
              do i = 1,nr+1
                car = vec(i,j,k,icomp+1:icomp+3)
                vec(i,j,k,icomp+1:icomp+3)
     .             = XformToCurv_ijk(gl_def,i,j,k,igr,car,iscnv)
              enddo
            enddo
          enddo
        endif

      enddo

      deallocate(sp_sum,dv_sum)
        
c     End program

      end subroutine averageSingularBC

      end module singularBCinterface

c module imposeBCinterface
c #####################################################################
      module imposeBCinterface

        use singularBCinterface

        use oned_int

        INTERFACE setBC
          module procedure imposeBConScalar,imposeBConVector
        end INTERFACE

        INTERFACE setBC_ez
          module procedure imposeBConScalar_ez,imposeBConVector_ez
        end INTERFACE

        INTERFACE XferBCs
          module procedure XferBCs_scl, XferBCs_vec
        END INTERFACE 

        integer  :: iimax,iimin,jjmax,jjmin,kkmax,kkmin
        integer  :: iimxp,iimnm,jjmxp,jjmnm,kkmxp,kkmnm

        integer,private :: bctype

        logical,private :: skip_corners=.false.
cc     .                    ,force_corners=.false.

        !Define order of interpolation for EDGE and CORNER treatment
        !(linear to preserve nearest-neighbors stencil; this is 
        ! important for colored diagonal formation)
        integer,private :: order=1

        logical,private :: stg_vector=.false.

        private :: setASMflags,interpolateToFace
     .            ,fillEdges,fillCorners,do__nothing

      contains

c     set_stg_vector
c     ##############################################################
      subroutine set_stg_vector(stg)

      implicit none

c     --------------------------------------------------------------
c     Whether to skip corners in BC treatment
c     --------------------------------------------------------------

c     Call variables

      logical :: stg

c     Begin program

      stg_vector = stg

c     End program

      end subroutine set_stg_vector

c     get_stg_vector
c     ##############################################################
      function get_stg_vector() result(stg)

      implicit none

c     --------------------------------------------------------------
c     Whether to skip corners in BC treatment
c     --------------------------------------------------------------

c     Call variables

      logical :: stg

c     Begin program

      stg = stg_vector

c     End program

      end function get_stg_vector

c     set_skip_corners
c     ##############################################################
      subroutine set_skip_corners(sc)

      implicit none

c     --------------------------------------------------------------
c     Whether to skip corners in BC treatment
c     --------------------------------------------------------------

c     Call variables

      logical :: sc

c     Begin program

      skip_corners = sc

c     End program

      end subroutine set_skip_corners

c     set_force_corners
c     ##############################################################
      subroutine set_force_corners(fc)

      implicit none

c     --------------------------------------------------------------
c     Whether to force corners in BC treatment
c     --------------------------------------------------------------

c     Call variables

      logical :: fc

c     Begin program

cc      force_corners = fc

c     End program

      end subroutine set_force_corners

c     findBCLoopLimits
c     ##############################################################
      subroutine findBCLoopLimits(gbc_def,dim,loc,i1,i2,j1,j2,k1,k2
     .                                   ,imn,imx,jmn,jmx,kmn,kmx)
c     --------------------------------------------------------------
c     Finds local limits for BC loops. On input, it takes global 
c     limits (i1,i2), (j1,j2), (k1,k2) and outputs local loop limits
c     (imn,imx), (jmn,jmx), (kmn,kmx).
c     --------------------------------------------------------------

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer,intent(IN)  :: dim,loc,i1,i2,j1,j2,k1,k2
      integer,intent(OUT) :: imn,imx,jmn,jmx,kmn,kmx

c     Begin program

      imn = i1 - gbc_def%ilo(igxbc) + 1
      imx = i2 - gbc_def%ilo(igxbc) + 1

      jmn = j1 - gbc_def%jlo(igybc) + 1
      jmx = j2 - gbc_def%jlo(igybc) + 1

      kmn = k1 - gbc_def%klo(igzbc) + 1
      kmx = k2 - gbc_def%klo(igzbc) + 1

      select case(dim)
      case (1)

        call selectLimits(loc,nxbc,imn,imx)

      case(2)

        call selectLimits(loc,nybc,jmn,jmx)

      case(3)

        call selectLimits(loc,nzbc,kmn,kmx)

      end select

c     End program

      contains

c     selectLimits
c     #####################################################################
      subroutine selectLimits(loc,nl,ill,ilr)

c     ---------------------------------------------------------------------
c     Selects LOCAL limits for BC treatment at boundary location determined
c     by loc.
c     ---------------------------------------------------------------------

        implicit none

        integer :: loc,nl,ill,ilr

c     Begin program

        select case(loc)
        case(0)
          ill = min(ill,1)
          ilr = ill
        case(1)
          ilr = max(ilr,nl)
          ill = ilr
        end select

      end subroutine selectLimits

      end subroutine findBCLoopLimits

c     do__nothing
c     #####################################################################
      function do__nothing(gl_def,dim,loc,i1,i2,j1,j2,k1,k2)

c     ---------------------------------------------------------------------
c     Selects whether to fill or not ghost cells of a given local domain
c     based on its position in the global grid.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        logical :: do__nothing

        type(grid_mg_def),pointer :: gl_def

        integer,intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2

c     Local variables

c     Begin program

        select case(dim)
        case(1)
          call checkLimits(i1,i2,gl_def%nxgl(igxbc))
        case(2)
          call checkLimits(j1,j2,gl_def%nygl(igybc))
        case(3)
          call checkLimits(k1,k2,gl_def%nzgl(igzbc))
        end select

c     End program

      contains

c     checkLimits
c     #####################################################################
      subroutine checkLimits(igl,igr,ng)

c     ---------------------------------------------------------------------
c     Selects LOCAL limits for BC treatment at boundary location determined
c     by (dim,loc).
c     ---------------------------------------------------------------------

        implicit none

        integer :: igl,igr,ng

c     Begin program

        if (loc == 0) then
          do__nothing = (igl > 1 )
        else
          do__nothing = (igr < ng)
        endif

      end subroutine checkLimits

      end function do__nothing

c     selectFace
c     #####################################################################
      function selectFace(g_def,ibc)

c     ---------------------------------------------------------------------
c     Determines whether boundary location determined by ibc(=1-6) is
c     selected for boundary treatment.
c     ---------------------------------------------------------------------

      implicit none

      type(grid_mg_def),pointer :: g_def

      logical    :: selectFace

      integer    :: dim,loc,ibc

c     Begin program

      dim = (ibc+1)/2
      loc = ibc + 1 - 2*dim

      selectFace= (ibc >= 1 .and. ibc <= 6).and.(.not.
     .   do__nothing(g_def,dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax))

      end function selectFace

c     setASMflags
c     ######################################################################
      subroutine setASMflags(g_def,asm_flag)

c     -----------------------------------------------------------------
c     Sets flag to indicate Additive Schwartz method.
c     -----------------------------------------------------------------

      implicit none

      type(grid_mg_def),pointer :: g_def

      logical :: asm_flag

      asm = asm_flag

      asm_dir(1) = asm.and.(g_def%mx > 1)
      asm_dir(2) = asm.and.(g_def%my > 1)
      asm_dir(3) = asm.and.(g_def%mz > 1)

      if (bc_debug) then
        write (bc_debug_unit,*) 'Proc ',my_rank,': asm    =',asm
        write (bc_debug_unit,*) 'Proc ',my_rank,': asm_dir=',asm_dir
      endif

      end subroutine setASMflags

c     setup_petsc_BC
c     ######################################################################
      subroutine setup_petsc_BC(g_def)

c     -----------------------------------------------------------------
c     Sets flags to indicate that BCs are fully provided by PETSc.
c     -----------------------------------------------------------------

      type(grid_mg_def),pointer :: g_def

      call setASMflags(g_def,.false.)

      end subroutine setup_petsc_BC

c     setASM_PC_BC
c     ######################################################################
      subroutine setASM_PC_BC(gl_def,igx,igy,igz,bcs,pc_asm,IIBC)

c     -----------------------------------------------------------------
c     Sets boundary conditions for Additive Schwartz method PC.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gl_def

      integer :: igx,igy,igz,bcs(:,:)

      integer,optional :: IIBC

      logical :: pc_asm

c     Local variables

      integer :: dim,loc,ibc,igmin,igmax,jgmin,jgmax,kgmin,kgmax

      logical :: dn

c     Begin program

c     Reset Additive-Schwarz-Method flags

      if ((np == 1) .or. (.not.pc_asm)) then
        call setASMflags(gl_def,.false.)
        return
      else
        call setASMflags(gl_def,.true.)
      endif

c     Reset BC grid level

      igxbc = igx
      igybc = igy
      igzbc = igz

c     Find global limits of local domain

      igmin = gl_def%ilo(igx)
      igmax = gl_def%ihi(igx)
      jgmin = gl_def%jlo(igy)
      jgmax = gl_def%jhi(igy)
      kgmin = gl_def%klo(igz)
      kgmax = gl_def%khi(igz)

c     Determine ASM BCs 

      if (PRESENT(IIBC)) then
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            !Select logical and periodic boundaries
            dn=do__nothing(gl_def,dim,loc,igmin,igmax,jgmin,jgmax
     $                    ,kgmin,kgmax)
cc     .         .or.((bcond(ibc)==PER).and.asm_dir(dim))  !Deal with interrupted periodic boundaries
     .         .or.(bcond(ibc)==PER)  !Deal with interrupted periodic boundaries

            if (dn) then
              where (bcs(ibc,:) /= 0) 
                bcs(ibc,:) = IIBC*bcs(ibc,:)/abs(bcs(ibc,:))
              elsewhere
                bcs(ibc,:) = 0
              end where
            endif
          enddo
        enddo
      else
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
cc            if (bcond(ibc) == PER.and.asm_dir(dim)) then !Take care of PER boundaries in ASM
            if (bcond(ibc) == PER) then !Take care of PER boundaries in ASM
              bcs(ibc,:) = EXT
            endif
          enddo
        enddo
      endif

      end subroutine setASM_PC_BC

c     set_vec_symm_bc
c     ###################################################################
      subroutine set_vec_symm_bc(vbcs,normal_parity)

      implicit none

c     -------------------------------------------------------------------
c     Processes symmetry BCs for vectors based on parity of normal
c     component.
c     -------------------------------------------------------------------

c     Call variables

      integer :: vbcs(6,3),normal_parity

c     Local variables

      integer :: diml,locl,ibcl

c     Begin program

      do ibcl = 1,6

        locl = abs(mod(ibcl,2)-2)-1
        diml = 1 + (ibcl-1-locl)/2

        if (vbcs(ibcl,diml) == FSYM) then
          if (normal_parity == 1) then
            vbcs(ibcl,:)    =-DIR !Tangential components
            vbcs(ibcl,diml) = NEU !Normal component
          else
            vbcs(ibcl,:)    =-NEU !Tangential components
            vbcs(ibcl,diml) = DIR !Normal component
          endif
        endif

      enddo

c     End program

      end subroutine set_vec_symm_bc

c     imposeBConScalar_ez
c     #################################################################
      subroutine imposeBConScalar_ez(gbc_def,igr,array,bcnd,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: igr,bcnd(6)
      real(8)    :: array (0:,0:,0:)
      integer   ,optional,intent(IN) :: iorder

c     Local variables

      integer :: nx,ny,nz
      real(8) :: array0(0:size(array,1)-1
     .                 ,0:size(array,2)-1
     .                 ,0:size(array,3)-1)

c     Begin program

      nx = gbc_def%nxv(igr)
      ny = gbc_def%nyv(igr)
      nz = gbc_def%nzv(igr)

      array0 = 0d0

      call imposeBConScalar(gbc_def,0,nx,ny,nz,array,array0,bcnd
     .                     ,igr,igr,igr,iorder=iorder)

c     End program

      end subroutine imposeBConScalar_ez

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(gbc_def,ieq,nx,ny,nz,array,array0,bcnd
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: ieq,bcnd(6),nx,ny,nz,iigx,iigy,iigz
      real(8)    :: array (0:,0:,0:)
     .             ,array0(0:,0:,0:)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

c     Local variables

      integer :: dim,loc,ibc,xtrap_o
      logical :: petsc__bc

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

c     Initialize interpolation order (else use last set)

      if (PRESENT(iorder)) then
        order = iorder
cc      else
cc        order = 1
      endif

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize GLOBAL limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = gbc_def%ilo(iigx)
        iimax = gbc_def%ihi(iigx)
        jjmin = gbc_def%jlo(iigy)
        jjmax = gbc_def%jhi(iigy)
        kkmin = gbc_def%klo(iigz)
        kkmax = gbc_def%khi(iigz)
      endif

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,0)  
      iimxp =  min(iimax + 1,gbc_def%nxgl(iigx)+1)
      jjmnm =  max(jjmin - 1,0)  
      jjmxp =  min(jjmax + 1,gbc_def%nygl(iigy)+1)
      kkmnm =  max(kkmin - 1,0)  
      kkmxp =  min(kkmax + 1,gbc_def%nzgl(iigz)+1)
#endif

c     Fill PETSc ghost cells and PER BCs

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc =  (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => gbc_def%dactx(igxbc)
        call fillPetscGhostCells(da_ctx,array)
      endif
#elif !defined(samrai)
      call sync_bc(.true.)
#endif

      if (bcPER(1).and.bcPER(2).and.bcPER(3)) return !PETSc deals with all periodic domains

c     Impose BCs

      do bctype=2,BCLIM             !Do not do PER=1; already done
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (abs(bcnd(ibc)) == bctype) then
              call FillGhostNodes_scl(gbc_def,ieq,1,dim,loc,bctype
     .                               ,array,array0)
            endif
          enddo
        enddo
      enddo

c     Fill edges

      xtrap_o = max(2,order)

      call fillEdges(gbc_def,array,array0,bcnd,xtrap_o)

c     Fill corners

      call fillCorners(gbc_def,array,bcnd,xtrap_o)

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_bc(PER_SYNC_ONLY)
#endif

c     End program

      contains

c     sync_bc
c     #################################################################
      subroutine sync_bc(only_per)

c     -----------------------------------------------------------------
c     Synchronizes boundaries in parallel.
c     -----------------------------------------------------------------

      implicit none

c     Call variables
      
      logical :: only_per
      
c     Local variables

c     Begin program

      !Periodic: order by dimension
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          if (abs(bcnd(ibc)) == PER) then
            call FillGhostNodes_scl(gbc_def,ieq,1,dim,loc,PER,array
     .                             ,array0)
          endif
        enddo
      enddo

      if (only_per) return

      !Other: order by location
      do loc=0,1
        do dim=1,3
          ibc = (1+loc)+2*(dim-1)
          if (abs(bcnd(ibc)) /= PER) then
            call FillGhostNodes_scl(gbc_def,ieq,1,dim,loc,SYNC,array
     .                             ,array0)
          endif
        enddo
      enddo

      end subroutine sync_bc

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(gbc_def,fcomp,ncomp,nx,ny,nz
     .                           ,v_cnv,v_cov,var0
     .                           ,bcnd,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder,is_cnv,is_vec,sp_reg)
c     -----------------------------------------------------------------
c     Imposes BC on vector field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: ncomp,bcnd(6,ncomp),fcomp,nx,ny,nz,iigx,iigy,iigz
      real(8)    :: v_cnv(0:,0:,0:,:)
     .             ,var0 (0:,0:,0:,:)
     .             ,v_cov(0:,0:,0:,:)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder
      logical   ,optional,intent(IN) :: is_cnv,is_vec,sp_reg

c     Local variables

      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin
      integer    :: ivar,ieq,ibc,loc,dim,bctype,xtrap_o
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: cov_to_cnv(6,3),to_cnv,is__vec,petsc__bc
     .             ,regularize_sp

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

      !By default, we assume input are vector components
      if (PRESENT(is_vec)) then
        is__vec = is_vec
      else
        is__vec = .true.
      endif

      if (is__vec .and. ncomp < 3) then
        call pstop('imposeBConVector'
     .            ,'Cannot impose BC on vector; # components < 3')
      endif

      !By default, we assume input is in cnv representation
      if (PRESENT(is_cnv)) then
        to_cnv = .not.is_cnv
      else
        to_cnv = .false.
      endif

      !By default, we assume input is in cnv representation
      if (PRESENT(sp_reg)) then
        regularize_sp = sp_reg
      else
        regularize_sp = .false.
      endif

c     Initialize interpolation order (else use last set)

      if (PRESENT(iorder)) then
        order = iorder
cc      else
cc        order = 1
      endif

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize GLOBAL limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = gbc_def%ilo(iigx)
        iimax = gbc_def%ihi(iigx)
        jjmin = gbc_def%jlo(iigy)
        jjmax = gbc_def%jhi(iigy)
        kkmin = gbc_def%klo(iigz)
        kkmax = gbc_def%khi(iigz)
      endif

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,0)  
      iimxp =  min(iimax + 1,gbc_def%nxgl(iigx)+1)
      jjmnm =  max(jjmin - 1,0)  
      jjmxp =  min(jjmax + 1,gbc_def%nygl(iigy)+1)
      kkmnm =  max(kkmin - 1,0)  
      kkmxp =  min(kkmax + 1,gbc_def%nzgl(iigz)+1)
#endif

c     Prepare SP BCs

      if (isSP2(gbc_def,1,igxbc).and.regularize_sp) then
        call averageSingularBC(gbc_def,igxbc,v_cnv,0,2,.true. ,.false.)
        call averageSingularBC(gbc_def,igxbc,v_cov,0,2,.false.,.false.)
      endif

c     Fill PETSc ghost cells and PER BCs

      call set_LOG_bc(v_cnv,.true.)
      call set_LOG_bc(v_cov,.false.)

      if (bcPER(1).and.bcPER(2).and.bcPER(3)) return

c     Set SP BCs (if appropriate)
      
      call set_TOPOL_bc(SP,v_cnv,.true.)
      call set_TOPOL_bc(SP,v_cov,.false.)

c     Impose physical BCs

      cov_to_cnv = .false.

      do bctype=1,BCLIM
        if (bctype == SP.or.bctype == PER) cycle !Already done
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              ieq = ivar + fcomp - 1

              if (abs(bcnd(ibc,ivar)) == bctype) then

                if (bcnd(ibc,ivar) < 0 .and. is__vec) then
                  call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                                   ,bctype,v_cov,var0,cnv=.false.)
                  cov_to_cnv(ibc,ivar) = .true.
                else
                  call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                                   ,bctype,v_cnv,var0)
                endif

              endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components

      call synchronize(ncomp,v_cnv,v_cov,bcnd,cov_to_cnv)

c     Fill edges (for both cnv and cov components)

      xtrap_o = max(2,order)

      do ivar=1,ncomp
        call fillEdges(gbc_def,v_cov(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar),xtrap_o)
        call fillEdges(gbc_def,v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar),xtrap_o)
      enddo

c     Fill corners (for both cnv and cov components)

      do ivar=1,ncomp
        call fillCorners(gbc_def,v_cov(:,:,:,ivar),bcnd(:,ivar),xtrap_o)
        call fillCorners(gbc_def,v_cnv(:,:,:,ivar),bcnd(:,ivar),xtrap_o)
      enddo

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_bc(PER_SYNC_ONLY)
#endif

c     End program

      contains

c     set_LOG_bc
c     #################################################################
      subroutine set_LOG_bc(vec,cnv)

      implicit none
      
c     -----------------------------------------------------------------
c     Set BC "bctype" independently for cov and cnv vectors
c     -----------------------------------------------------------------

      real(8) :: vec(:,:,:,:)
      logical :: cnv,petsc__bc
      integer :: ieq

#if defined(petsc) && !defined(samrai)
        !Determine whether we need to fill ghost cells using PETSc
        petsc__bc = (.not.asm)

        !Fill ghost cells
        if (petsc__bc) then
          da_ctx => gbc_def%dactx(igxbc)
          do ieq=1,ncomp
            call fillPetscGhostCells(da_ctx,vec(:,:,:,ieq))
          enddo
        endif
#elif !defined(samrai)
        call set_TOPOL_bc(PER,vec,cnv)
#endif

      end subroutine set_LOG_bc

c     set_TOPOL_bc
c     #################################################################
      subroutine set_TOPOL_bc(bctype,vec,cnv)

      implicit none
      
c     -----------------------------------------------------------------
c     Set topological BCs independently for cov and cnv vectors
c     -----------------------------------------------------------------

      integer :: bctype
      real(8) :: vec(:,:,:,:)
      logical :: cnv

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            ieq = ivar + fcomp - 1

            if (abs(bcnd(ibc,ivar)) == bctype) then
              call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                               ,bctype,vec,var0,cnv=cnv)
            endif
          enddo

        enddo
      enddo

      end subroutine set_TOPOL_bc

c     sync_bc
c     #################################################################
      subroutine sync_bc(only_per)

c     Call variables

      logical :: only_per

c     Local variables
      
c     Begin program

      !Periodic: order by dimension first
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            ieq = ivar + fcomp - 1

            if (abs(bcnd(ibc,ivar)) == PER) then
              call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                               ,PER,v_cnv,var0)
              call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                               ,PER,v_cov,var0,cnv=.false.)
            endif
          enddo

        enddo
      enddo

      if (only_per) return

      !Other: order by location first
      do loc=0,1
        do dim=1,3
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            ieq = ivar + fcomp - 1

            if (abs(bcnd(ibc,ivar)) /= PER) then
              call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                               ,SYNC,v_cnv,var0)
              call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                               ,SYNC,v_cov,var0,cnv=.false.)
            endif
          enddo

        enddo
      enddo

      end subroutine sync_bc

c     synchronize
c     #################################################################
      subroutine synchronize(ncomp,v_cnv,v_cov,bcond,cov_to_cnv,dim,loc)
c     -----------------------------------------------------------------
c     Finds all vector components at physical boundaries.
c
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravariant
c     and covariant components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ncomp,bcond(6,ncomp)
      real(8)    :: v_cnv(0:nxbc+1,0:nybc+1,0:nzbc+1,ncomp)
     .             ,v_cov(0:nxbc+1,0:nybc+1,0:nzbc+1,ncomp)
      logical    :: cov_to_cnv(6,3)
      integer,optional :: dim,loc

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,ig,jg,kg,ivar,ibc,ic1,ic2
     .             ,ld,hd,ll,hl,id,il
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => gbc_def%gmetric

      if (PRESENT(dim)) then
        ld = dim
        hd = dim
      else
        ld = 1
        hd = 3
      endif

      if (PRESENT(loc)) then
        ll = loc
        hl = loc
      else
        ll = 0
        hl = 1
      endif

      do id = ld,hd
        do il = ll,hl
          ibc = (1+il)+2*(id-1)

          if (.not.selectFace(gbc_def,ibc)
     .       .or.bcPER(id)
     .       .or.bcSP(ibc)) cycle !Do nothing for logical, periodic, or SP boundaries

          !Do NOT change: SAMRAI depends on it
          call findBCLoopLimits(gbc_def
     .                         ,id,il,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                               ,imin,imax
     .                               ,jmin,jmax
     .                               ,kmin,kmax)

          if (     (.not.cov_to_cnv(ibc,1))
     .        .and.(.not.cov_to_cnv(ibc,2))
     .        .and.(.not.cov_to_cnv(ibc,3))) then  !All contravariant: find covariant

            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax
                  v_cov(i,j,k,:) =
     .                XformToCov_ijk(gbc_def,i,j,k,igxbc,v_cnv(i,j,k,:))
                enddo
              enddo
            enddo

          else   !Synchronize (assumes normal component is contravariant)

cc            if (.not.selectFace(gbc_def,ibc)) cycle !Do nothing for logical boundaries

            do ivar = 1,ncomp

              if (ivar == id) then !Do nothing for normal component

                cycle

              else

cc                if (.not.selectFace(gbc_def,ibc)) cycle

                ic1 = mod(id+1,3)
                if (ic1 == 0) ic1 = 3

                ic2 = mod(id+2,3)
                if (ic2 == 0) ic2 = 3

                if ( cov_to_cnv(ibc,ic1).and.cov_to_cnv(ibc,ic2) ) then  !Both are covariant

                  do k=kmin,kmax
                    do j=jmin,jmax
                      do i=imin,imax

                        gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)

                        v_cov(i,j,k,id) =
     .                       -(gsuper(id,ic1)*v_cov(i,j,k,ic1)
     .                        +gsuper(id,ic2)*v_cov(i,j,k,ic2)
     .                        -v_cnv(i,j,k,id))/gsuper(id,id)

                        v_cnv(i,j,k,:) = XformToCnv_ijk(gbc_def,i,j,k
     .                                           ,igxbc,v_cov(i,j,k,:))
                      enddo
                    enddo
                  enddo

                elseif (.not.cov_to_cnv(ibc,ic1)) then   !ic1 is contravariant, ic2 is covariant

                  do k=kmin,kmax
                    do j=jmin,jmax
                      do i=imin,imax

                        gsub = gmetric%grid(igxbc)%gsub(i,j,k,:,:)

                        v_cnv(i,j,k,ic2) =
     .                       -(gsub(ic2,ic1)*v_cnv(i,j,k,ic1)
     .                        +gsub(ic2,id )*v_cnv(i,j,k,id)
     .                        -v_cov(i,j,k,ic2))/gsub(ic2,ic2)

                        v_cov(i,j,k,:) = XformToCov_ijk(gbc_def,i,j,k
     .                                           ,igxbc,v_cnv(i,j,k,:))
                      enddo
                    enddo
                  enddo

                elseif (.not.cov_to_cnv(ibc,ic2)) then   !ic2 is contravariant, ic1 is covariant

                  do k=kmin,kmax
                    do j=jmin,jmax
                      do i=imin,imax

                        gsub = gmetric%grid(igxbc)%gsub(i,j,k,:,:)

                        v_cnv(i,j,k,ic1) =
     .                       -(gsub(ic1,ic2)*v_cnv(i,j,k,ic2)
     .                        +gsub(ic1,id )*v_cnv(i,j,k,id)
     .                        -v_cov(i,j,k,ic1))/gsub(ic1,ic1)

                        v_cov(i,j,k,:) = XformToCov_ijk(gbc_def,i,j,k
     .                                           ,igxbc,v_cnv(i,j,k,:))
                      enddo
                    enddo
                  enddo

                endif

             endif
            enddo

          endif

        enddo
      enddo

c     End program

      end subroutine synchronize

      end subroutine imposeBConVector

c     imposeBConVector_ez
c     #################################################################
      subroutine imposeBConVector_ez(gbc_def,igrid,vec,is_cnv,bcnd
     .                              ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on vector field (easy interface: no cov-cnv coupling)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer :: bcnd(:,:),igrid,iorder
      real(8) :: vec (0:,0:,0:,:)

      logical :: is_cnv

c     Local variables

      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin,nx,ny,nz
      integer    :: ncomp,ivar,ibc,loc,dim,bctype,xtrap_o
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: petsc__bc

      real(8),pointer :: vec0(:,:,:,:)

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

      nx = size(vec,1)-2
      ny = size(vec,2)-2
      nz = size(vec,3)-2
      ncomp = size(vec,4)

      allocate(vec0(0:nx+1,0:ny+1,0:nz+1,3))
      vec0 = 0d0

c     Initialize BC grid quantities

      igxbc = igrid
      igybc = igrid
      igzbc = igrid

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize GLOBAL limits

      iimin = gbc_def%ilo(igxbc)
      iimax = gbc_def%ihi(igxbc)
      jjmin = gbc_def%jlo(igybc)
      jjmax = gbc_def%jhi(igybc)
      kkmin = gbc_def%klo(igzbc)
      kkmax = gbc_def%khi(igzbc)

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,0)  
      iimxp =  min(iimax + 1,gbc_def%nxgl(igxbc)+1)
      jjmnm =  max(jjmin - 1,0)  
      jjmxp =  min(jjmax + 1,gbc_def%nygl(igybc)+1)
      kkmnm =  max(kkmin - 1,0)  
      kkmxp =  min(kkmax + 1,gbc_def%nzgl(igzbc)+1)
#endif

c     Fill PETSc ghost cells and PER BCs

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc = (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => gbc_def%dactx(igxbc)

        do ivar=1,ncomp
          call fillPetscGhostCells(da_ctx,vec(:,:,:,ivar))
        enddo
      endif
#elif !defined(samrai)
      call sync_bc(.true.)
#endif

      if (bcPER(1).and.bcPER(2).and.bcPER(3)) return

c     Impose physical BCs

      do bctype=2,BCLIM             !Do not do PER=1; already done 
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              if (abs(bcnd(ibc,ivar)) == bctype) then
                call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                                 ,bctype,vec,vec0,cnv=is_cnv)
              endif
            enddo

          enddo
        enddo
      enddo

c     Fill edges 

      xtrap_o = max(2,order)

      do ivar=1,ncomp
        call fillEdges(gbc_def,vec(:,:,:,ivar),vec0(:,:,:,ivar)
     .                ,bcnd(:,ivar),xtrap_o)
      enddo

c     Fill corners 

      do ivar=1,ncomp
        call fillCorners(gbc_def,vec(:,:,:,ivar),bcnd(:,ivar),xtrap_o)
      enddo

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_bc(PER_SYNC_ONLY)
#endif

      deallocate(vec0)

c     End program

      contains

c     sync_bc
c     #################################################################
      subroutine sync_bc(only_per)

c     Call variables

      logical :: only_per

c     Local variables
      
c     Begin program

      !Periodic: order by dimension first
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            if (abs(bcnd(ibc,ivar)) == PER) then
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,PER,vec,vec0,cnv=is_cnv)
            endif
          enddo

        enddo
      enddo

      if (only_per) return

      !Other: order by location first
      do loc=0,1
        do dim=1,3
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            if (abs(bcnd(ibc,ivar)) /= PER) then
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,SYNC,vec,vec0,cnv=is_cnv)
            endif
          enddo

        enddo
      enddo

      end subroutine sync_bc

      end subroutine imposeBConVector_ez

c     imposeBConVector_stg
c     #################################################################
      subroutine imposeBConVector_stg(gbc_def,fcomp,ncomp,nx,ny,nz
     .          ,v_cnv,v_cov,var0,bcnd,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .          ,iorder,is_cnv)
c     -----------------------------------------------------------------
c     Imposes BC on staggered vector field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ncomp,bcnd(6,ncomp),fcomp,nx,ny,nz,iigx,iigy,iigz
      real(8) :: v_cnv(0:nx+1,0:ny+1,0:nz+1,ncomp)
     .          ,var0 (0:nx+1,0:ny+1,0:nz+1,ncomp)
     .          ,v_cov(0:nx+1,0:ny+1,0:nz+1,ncomp)
      integer,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder
      logical,optional,intent(IN) :: is_cnv

      type(grid_mg_def),pointer :: gbc_def

c     Local variables

      integer :: i,j,k,imax,imin,jmax,jmin,kmax,kmin,xtrap_o
      integer :: ivar,ieq,ibc,loc,dim,bctype
     .          ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical :: cov_to_cnv(6,3),to_cnv,flg_sc

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

      !Indicate staggered vector treatment
      call set_stg_vector(.true.)

      !By default, we assume input is in cnv representation
      if (PRESENT(is_cnv)) then
        to_cnv = .not.is_cnv
      else
        to_cnv = .false.
      endif

      flg_sc = .false.
      if (.not.skip_corners) then
        skip_corners=.true.     !Skip BC corner treatment
        flg_sc = .true.
      endif

c     Initialize interpolation order (else use last set)

      if (PRESENT(iorder)) then
        order = iorder
cc      else
cc        order = 1
      endif

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize GLOBAL limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = gbc_def%ilo(iigx)
        iimax = gbc_def%ihi(iigx)
        jjmin = gbc_def%jlo(iigy)
        jjmax = gbc_def%jhi(iigy)
        kkmin = gbc_def%klo(iigz)
        kkmax = gbc_def%khi(iigz)
      endif

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,0)  
      iimxp =  min(iimax + 1,gbc_def%nxgl(iigx)+1)
      jjmnm =  max(jjmin - 1,0)  
      jjmxp =  min(jjmax + 1,gbc_def%nygl(iigy)+1)
      kkmnm =  max(kkmin - 1,0)  
      kkmxp =  min(kkmax + 1,gbc_def%nzgl(iigz)+1)
#endif

c     Find cnv components

      if (to_cnv) then

        !Cov ghost cells
        call set_LOG_bc(v_cov,.false.)
        call set_TOPOL_bc(SP,v_cov,.false.)

        cov_to_cnv = .false.

        do bctype=2,BCLIM       !Do not do PER=1; already done by PETSc or SAMRAI
          if (bctype == SP) cycle
          call set_bc(bctype,v_cov,.false.)
        enddo

        !Xform to cnv
        v_cnv =  XformToCnv_stg_mesh(gbc_def,igxbc,v_cov)

        !Cnv ghost cells
        call set_LOG_bc(v_cnv,.true.)
        call set_TOPOL_bc(SP,v_cnv,.true.)

        do bctype=2,BCLIM       !Do not do PER=1; already done by PETSc or SAMRAI
          if (bctype == SP) cycle
          call set_bc(bctype,v_cnv,.true.)
        enddo

      else

        !Cnv ghost cells
        call set_LOG_bc(v_cnv,.true.)
        call set_TOPOL_bc(SP,v_cnv,.true.)

        do bctype=2,BCLIM       !Do not do PER=1; already done by PETSc or SAMRAI
          if (bctype == SP) cycle
          call set_bc(bctype,v_cnv,.true.)
        enddo

        !Xform to Cov
        v_cov =  XformToCov_stg_mesh(gbc_def,igxbc,v_cnv)

        !Cov ghost cells
        call set_LOG_bc(v_cov,.false.)
        call set_TOPOL_bc(SP,v_cov,.false.)

        cov_to_cnv = .false.

        do bctype=2,BCLIM       !Do not do PER=1; already done by PETSc or SAMRAI
          if (bctype == SP) cycle
          call set_bc(bctype,v_cov,.false.)
        enddo

      endif

c     Synchronize covariant and contravariant components (only bottom bdrys)

      call synchronize_stg(ncomp,v_cnv,v_cov,bcnd,cov_to_cnv)

c$$$c     Fill edges (for both cnv and cov components)
c$$$
c$$$      xtrap_o = max(2,order)
c$$$
c$$$      do ivar=1,ncomp
c$$$        call fillEdges(gbc_def,v_cov(:,:,:,ivar),var0(:,:,:,ivar)
c$$$     .                ,bcnd(:,ivar),xtrap_o)
c$$$        call fillEdges(gbc_def,v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
c$$$     .                ,bcnd(:,ivar),xtrap_o)
c$$$      enddo
c$$$
c$$$c     Fill corners (for both cnv and cov components)
c$$$
c$$$      do ivar=1,ncomp
c$$$        call fillCorners(gbc_def,v_cov(:,:,:,ivar),bcnd(:,ivar),xtrap_o)
c$$$        call fillCorners(gbc_def,v_cnv(:,:,:,ivar),bcnd(:,ivar),xtrap_o)
c$$$      enddo

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_bc(PER_SYNC_ONLY)
#endif

c     End program

      call set_stg_vector(.false.)
      if (flg_sc) skip_corners=.false. !Reset skip_corners

      contains

c     set_LOG_bc
c     #################################################################
      subroutine set_LOG_bc(vec,cnv)

      implicit none
      
c     -----------------------------------------------------------------
c     Set BC "bctype" independently for cov and cnv vectors
c     -----------------------------------------------------------------

      real(8) :: vec(:,:,:,:)
      logical :: cnv,petsc__bc
      integer :: ieq

#if defined(petsc) && !defined(samrai)
        !Determine whether we need to fill ghost cells using PETSc
        petsc__bc = (.not.asm)

        !Fill ghost cells
        if (petsc__bc) then
          da_ctx => gbc_def%dactx(igxbc)
          do ieq=1,ncomp
            call fillPetscGhostCells(da_ctx,vec(:,:,:,ieq))
          enddo
        endif
#elif !defined(samrai)
        call set_TOPOL_bc(PER,vec,cnv)
#endif

      end subroutine set_LOG_bc

c     set_TOPOL_bc
c     #################################################################
      subroutine set_TOPOL_bc(bctype,vec,cnv)

      implicit none
      
c     -----------------------------------------------------------------
c     Set topological BCs independently for cov and cnv vectors
c     -----------------------------------------------------------------

      integer :: bctype
      real(8) :: vec(:,:,:,:)
      logical :: cnv

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            ieq = ivar + fcomp - 1

            if (abs(bcnd(ibc,ivar)) == bctype) then
              call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                               ,bctype,vec,var0,cnv=cnv)
            endif
          enddo

        enddo
      enddo

      end subroutine set_TOPOL_bc

c     set_bc
c     #################################################################
      subroutine set_bc(bctype,vec,cnv)

      implicit none
      
c     -----------------------------------------------------------------
c     Set BC "bctype" independently for cov and cnv vectors
c     -----------------------------------------------------------------

c     Call variables

      integer :: bctype
      real(8) :: vec(:,:,:,:)
      logical :: cnv

c     Local variables

      integer :: dim,loc,ibc,ivar,ieq

c     Begin program

      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            ieq = ivar + fcomp - 1

            if (abs(bcnd(ibc,ivar)) == bctype) then
              if ((.not.cnv).and.bcnd(ibc,ivar) < 0) then
                call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                                 ,bctype,vec,var0,cnv=.false.)
                cov_to_cnv(ibc,ivar) = .true.
              elseif (cnv.and.bcnd(ibc,ivar) > 0) then
                call FillGhostNodes_vec(gbc_def,ieq,ivar,ncomp,dim,loc
     .                                 ,bctype,vec,var0,cnv=.true.)
              endif
            endif
          enddo

        enddo
      enddo

      end subroutine set_bc

c     sync_bc
c     #################################################################
      subroutine sync_bc(only_per)

c     Call variables

      logical :: only_per

c     Local variables
      
c     Begin program

      !Periodic: order by dimension first
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            if (abs(bcnd(ibc,ivar)) == PER) then
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,PER,v_cnv,var0,cnv=.true.)
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,PER,v_cov,var0,cnv=.false.)
            endif
          enddo

        enddo
      enddo

      if (only_per) return

      !Other: order by location first
      do loc=0,1
        do dim=1,3
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            if (abs(bcnd(ibc,ivar)) /= PER) then
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,SYNC,v_cnv,var0,cnv=.true.)
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,SYNC,v_cov,var0,cnv=.false.)
            endif
          enddo

        enddo
      enddo

      end subroutine sync_bc

c     synchronize_stg
c     #################################################################
      subroutine synchronize_stg(ncomp,v_cnv,v_cov,bcond,cov_to_cnv,dim
     $                          ,loc)
c     -----------------------------------------------------------------
c     Finds all vector components at physical boundaries.
c
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravariant
c     and covariant components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ncomp,bcond(6,ncomp)
      real(8) :: v_cnv(0:,0:,0:,:)
     .          ,v_cov(0:,0:,0:,:)
      logical :: cov_to_cnv(6,3)
      integer,optional :: dim,loc

c     Local variables

      integer :: imax,imin,jmax,jmin,kmax,kmin
      integer :: i,j,k,ig,jg,kg,ivar,ibc,ic1,ic2
     .          ,ld,hd,ll,hl,id,il,im,jm,km
      real(8) :: x1,x2,x3,gsup,gsub,cnv(3),cov(3)

      real(8),dimension(:,:,:,:,:),pointer :: gtensor

c     Begin program

      gtensor => gbc_def%gmetric%grid(igxbc)%gsup_v

      if (PRESENT(dim)) then
        ld = dim
        hd = dim
      else
        ld = 1
        hd = 3
      endif

      if (PRESENT(loc)) then
        ll = loc
        hl = loc
      else
        ll = 0
        hl = 1
      endif

      do id = ld,hd
        do il = ll,hl
          ibc = (1+il)+2*(id-1)

          if ((.not.selectFace(gbc_def,ibc))
     .       .or.bcPER(id)
     .       .or.bcSP(ibc)) cycle !Do nothing for logical, periodic, or SP boundaries

          !Do NOT change: SAMRAI depends on it
          call findBCLoopLimits(gbc_def
     .                         ,id,il,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                               ,imin   ,imax
     .                               ,jmin   ,jmax
     .                               ,kmin   ,kmax)

          if (     (.not.cov_to_cnv(ibc,1))
     .        .and.(.not.cov_to_cnv(ibc,2))
     .        .and.(.not.cov_to_cnv(ibc,3))) then  !All contravariant: find covariant

            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax
                   v_cov(i,j,k,:)=
     .                 XformToCov_stg_ijk(gbc_def,i,j,k,igxbc,v_cnv)
                enddo
              enddo
            enddo

          else     !Synchronize (assumes normal component is contravariant)

            ic1 = mod(id+1,3)
            if (ic1 == 0) ic1 = 3

            ic2 = mod(id+2,3)
            if (ic2 == 0) ic2 = 3

            if ( cov_to_cnv(ibc,ic1).and.cov_to_cnv(ibc,ic2) ) then  !Both tangentials are covariant

              do k=kmin,kmax
                do j=jmin,jmax
                  do i=imin,imax

                    im = i-1
                    jm = j-1
                    km = k-1

                    select case(id)
                    case(1)
                      if (jm < 0 .or. km < 0) cycle
                      gsup = 0.25*(gtensor(i,j ,k ,1,1)
     .                            +gtensor(i,jm,k ,1,1)
     .                            +gtensor(i,j ,km,1,1)
     .                            +gtensor(i,jm,km,1,1))
                    case(2)
                      if (im < 0 .or. km < 0) cycle
                      gsup = 0.25*(gtensor(i ,j,k ,2,2)
     .                            +gtensor(im,j,k ,2,2)
     .                            +gtensor(i ,j,km,2,2)
     .                            +gtensor(im,j,km,2,2))
                    case(3)
                      if (im < 0 .or. jm < 0) cycle
                      gsup = 0.25*(gtensor(i ,j ,k,3,3)
     .                            +gtensor(im,j ,k,3,3)
     .                            +gtensor(i ,jm,k,3,3)
     .                            +gtensor(im,jm,k,3,3))
                    end select

                    v_cov(i,j,k,id) = 0d0

                    cnv =XformToCnv_stg_ijk(gbc_def,i,j,k,igxbc,v_cov)
                    v_cov(i,j,k,id) =(v_cnv(i,j,k,id) - cnv(id))/gsup
                  enddo
                enddo
              enddo

              !Update contravariant tangential components
              do k=kmin,kmax
                do j=jmin,jmax
                  do i=imin,imax
                    cnv =XformToCnv_stg_ijk(gbc_def,i,j,k,igxbc,v_cov)
                    v_cnv(i,j,k,ic1) = cnv(ic1)
                    v_cnv(i,j,k,ic2) = cnv(ic2)
                  enddo
                enddo
              enddo

            elseif (.not.cov_to_cnv(ibc,ic1)) then   !ic1 is contravariant, ic2 is covariant

              call pstop('sync_stg'
     .                  ,'Code untested: ic1 cnv, ic2 cov')

cc              do k=kmin,kmax
cc                do j=jmin,jmax
cc                  do i=imin,imax
cc
cc                    im = max(i-1,0)
cc                    jm = max(j-1,0)
cc                    km = max(k-1,0)
cc
cccc                    if(ic2/=3.and.(i>nxbc.or.j>nybc.or.k>nzbc))cycle  !2D fix for now
cc
cc                    select case(ic2)
cc                    case(1)
cc                      gsub =
cc     .               .5*(gmetric%grid(igxbc)%gsub_v(i,j ,k,ic2,ic2)
cc     .                  +gmetric%grid(igxbc)%gsub_v(i,jm,k,ic2,ic2))
cc                    case(2)
cc                      gsub =
cc     .               .5*(gmetric%grid(igxbc)%gsub_v(im,j,k,ic2,ic2)
cc     .                  +gmetric%grid(igxbc)%gsub_v(i ,j,k,ic2,ic2))
cc                    case(3)
cc                      gsub = 1d0
cc                    end select
cc
cc                    v_cnv(i,j,k,ic2) = 0d0
cc                    cov =
cc     .                 XformToCov_stg_ijk(gbc_def,i,j,k,igxbc,v_cnv)
cc                    v_cnv(i,j,k,ic2)=-cov(ic2)-v_cov(i,j,k,ic2)/gsub
cc
cc                    v_cov(i,j,k,:) =
cc     .                 XformToCov_stg_ijk(gbc_def,i,j,k,igxbc,v_cnv)
cc
cc                  enddo
cc                enddo
cc              enddo

            elseif (.not.cov_to_cnv(ibc,ic2)) then   !ic2 is contravariant, ic1 is covariant

              call pstop('sync_stg'
     .                  ,'Code untested: ic1 cov, ic2 cnv')

cc              do k=kmin,kmax
cc                do j=jmin,jmax
cc                  do i=imin,imax
cc
cc                    im = max(i-1,0)
cc                    jm = max(j-1,0)
cc                    km = max(k-1,0)
cc
cccc                    if(ic1/=3.and.(i>nxbc.or.j>nybc.or.k>nzbc))cycle
cc
cc                    select case(ic1)
cc                    case(1)
cc                      gsub =
cc     .               .5*(gmetric%grid(igxbc)%gsub_v(i,j ,k,ic1,ic1)
cc     .                  +gmetric%grid(igxbc)%gsub_v(i,jm,k,ic1,ic1))
cc                    case(2)
cc                      gsub =
cc     .               .5*(gmetric%grid(igxbc)%gsub_v(im,j,k,ic1,ic1)
cc     .                  +gmetric%grid(igxbc)%gsub_v(i ,j,k,ic1,ic1))
cc                    case(3)
cc                      gsub = 1d0
cc                    end select
cc
cc                    v_cnv(i,j,k,ic1) = 0d0
cc                    cov =
cc     .                 XformToCov_stg_ijk(gbc_def,i,j,k,igxbc,v_cnv)
cc                    v_cnv(i,j,k,ic1)=-cov(ic1)-v_cov(i,j,k,ic1)/gsub
cc
cc                    v_cov(i,j,k,:) =
cc     .                 XformToCov_stg_ijk(gbc_def,i,j,k,igxbc,v_cnv)
cc
cc                  enddo
cc                enddo
cc              enddo

            endif

          endif
        enddo
      enddo

c     End program

      end subroutine synchronize_stg

      end subroutine imposeBConVector_stg

c     imposeBConVector_stg_ez
c     #################################################################
      subroutine imposeBConVector_stg_ez(gbc_def,igrid,vec,is_cnv,bcnd
     .                                  ,iorder,fcomp,vec0)
c     -----------------------------------------------------------------
c     Imposes BC on vector field (easy interface: no cov-cnv coupling)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer :: bcnd(:,:),igrid,iorder
      real(8) :: vec(0:,0:,0:,:)

      logical :: is_cnv

      integer,optional :: fcomp
      real(8),optional,target :: vec0(0:,0:,0:,:)

c     Local variables

      integer :: i,j,k,imax,imin,jmax,jmin,kmax,kmin,nx,ny,nz
      integer :: ncomp,ivar,ibc,loc,dim,bctype,xtrap_o
     .          ,ilmin,jlmin,klmin,ilmax,jlmax,klmax,fcmp
      logical :: petsc__bc,flg_sc

      real(8),pointer :: vec_0(:,:,:,:)

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

      !Indicate staggered vector treatment
      call set_stg_vector(.true.)

      !Skip BC corner treatment
      flg_sc = .false.
      if (.not.skip_corners) then
        skip_corners=.true.
        flg_sc = .true.
      endif

      if (PRESENT(fcomp)) then
        fcmp = fcomp
      else
        fcmp = 0
      endif

c     Initialize BC grid quantities

      nx = size(vec,1)-2
      ny = size(vec,2)-2
      nz = size(vec,3)-2
      ncomp = size(vec,4)

      if (PRESENT(vec0)) then
        vec_0 => vec0
      else
        allocate(vec_0(0:nx+1,0:ny+1,0:nz+1,3))
        vec_0 = 0d0
      endif

c     Initialize GLOBAL limits

      igxbc = igrid
      igybc = igrid
      igzbc = igrid

      nxbc = nx
      nybc = ny
      nzbc = nz

      iimin = gbc_def%ilo(igxbc)
      iimax = gbc_def%ihi(igxbc)
      jjmin = gbc_def%jlo(igybc)
      jjmax = gbc_def%jhi(igybc)
      kkmin = gbc_def%klo(igzbc)
      kkmax = gbc_def%khi(igzbc)

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,0)  
      iimxp =  min(iimax + 1,gbc_def%nxgl(igxbc)+1)
      jjmnm =  max(jjmin - 1,0)  
      jjmxp =  min(jjmax + 1,gbc_def%nygl(igybc)+1)
      kkmnm =  max(kkmin - 1,0)  
      kkmxp =  min(kkmax + 1,gbc_def%nzgl(igzbc)+1)
#endif

c     Fill PETSc ghost cells and PER BCs

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc = (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => gbc_def%dactx(igxbc)

        do ivar=1,ncomp
          call fillPetscGhostCells(da_ctx,vec(:,:,:,ivar))
        enddo
      endif
#elif !defined(samrai)
      call sync_bc(.true.)
#endif

      !Exit if all periodic domains 
      if (bcPER(1).and.bcPER(2).and.bcPER(3)) then
        call set_stg_vector(.false.)
        if (flg_sc) skip_corners=.false. !Reset skip_corners
        return
      endif

c     Impose physical BCs

      do bctype=2,BCLIM             !Do not do PER=1; already done 
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              if (abs(bcnd(ibc,ivar)) == bctype) then
                call FillGhostNodes_vec(gbc_def,fcmp,ivar,ncomp,dim,loc
     .                                 ,bctype,vec,vec_0,cnv=is_cnv)
              endif
            enddo

          enddo
        enddo
      enddo

c$$$c     Fill edges 
c$$$
c$$$      xtrap_o = max(2,order)
c$$$
c$$$      do ivar=1,ncomp
c$$$        call fillEdges(gbc_def,vec(:,:,:,ivar),vec_0(:,:,:,ivar)
c$$$     .                ,bcnd(:,ivar),xtrap_o)
c$$$      enddo
c$$$
c$$$c     Fill corners 
c$$$
c$$$      do ivar=1,ncomp
c$$$        call fillCorners(gbc_def,vec(:,:,:,ivar),bcnd(:,ivar),xtrap_o)
c$$$      enddo

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_bc(PER_SYNC_ONLY)
#endif

      if (PRESENT(vec0)) then
        nullify(vec_0)
      else
        deallocate(vec_0)
      endif

c     End program

      call set_stg_vector(.false.)
      if (flg_sc) skip_corners=.false. !Reset skip_corners

      contains

c     sync_bc
c     #################################################################
      subroutine sync_bc(only_per)

c     Call variables

      logical :: only_per

c     Local variables
      
c     Begin program

      !Periodic: order by dimension first
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            if (abs(bcnd(ibc,ivar)) == PER) then
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,PER,vec,vec_0,cnv=is_cnv)
            endif
          enddo

        enddo
      enddo

      if (only_per) return

      !Other: order by location first
      do loc=0,1
        do dim=1,3
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            if (abs(bcnd(ibc,ivar)) /= PER) then
              call FillGhostNodes_vec(gbc_def,0,ivar,ncomp,dim,loc
     .                               ,SYNC,vec,vec_0,cnv=is_cnv)
            endif
          enddo

        enddo
      enddo

      end subroutine sync_bc

      end subroutine imposeBConVector_stg_ez

c     fillGhostNodes_scl
c     ###############################################################
      subroutine fillGhostNodes_scl(gbc_def,ieq,ivar,dim,loc,bctype
     .                             ,array,array0)

c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * ieq    -> equation identifier
c       * ivar   -> vector component
c       * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c       * loc    -> location in dimension (0 -> right, 1 -> left)
c       * bctype -> type of BC (dirichlet, neumann, periodic, etc.)
c       * array  -> real array with ghost-nodes
c       * array0 -> auxiliary real array
c     ---------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: ieq,dim,loc,bctype,ivar
      real(8)    :: array (0:,0:,0:)
     .             ,array0(0:,0:,0:)

c     Local variables

      integer    :: neq,ibc,i,j,k,ig,jg,kg
      logical    :: cov

      real(8),allocatable,dimension(:,:,:) :: zeros

      INTERFACE
        subroutine scalarDirichletBC(g_def,arr,arr0,ieq,dim,loc,order)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: ieq,dim,loc,order
          real(8) :: arr (0:,0:,0:)
     .              ,arr0(0:,0:,0:)
        end subroutine scalarDirichletBC

        subroutine scalarNeumannBC(g_def,arr,arr0,ieq,dim,loc)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: ieq,dim,loc
          real(8) :: arr (0:,0:,0:)
     .              ,arr0(0:,0:,0:)
        end subroutine scalarNeumannBC
      END INTERFACE

c     Begin program

      ibc = (1+loc)+2*(dim-1)
      
c     Determine boundary limits

      if (.not.selectFace(gbc_def,ibc)) return

c     Find BC update

      select case(bctype)
      case(PER)
        call periodicBC(gbc_def,array,dim,loc)
      case(EQU)
        call scalarDirichletBC(gbc_def,array,array0
     .                        ,ieq,dim,loc,order)
      case(DIR)
        allocate(zeros(0:nxbc+1,0:nybc+1,0:nzbc+1))
        zeros = 0d0
        call scalarDirichletBC(gbc_def,array,zeros
     .                        ,ieq,dim,loc,order)
        deallocate(zeros)
      case(NEU)
        call scalarNeumannBC(gbc_def,array,array0
     .                      ,ieq,dim,loc)
      case(EXT)
        call fillFaces(gbc_def,array,dim,loc,order)
      case(SP)
        call scalarSingularBC(gbc_def,array)
      case(IFC)
        call interpolateToFace(gbc_def,array,dim,loc
     .                        ,order)
      case(SYNC) 
        call syncBC(gbc_def,array,dim,loc)
      case default
        call pstop('fillGhostNodes'
     .           ,'BC '//trim(int2char(bctype))//' not implemented')
      end select

c     End

      end subroutine fillGhostNodes_scl

c     fillGhostNodes_vec
c     ###############################################################
      subroutine fillGhostNodes_vec(gbc_def,ieq,ivar,nvar,dim,loc,bctype
     .                             ,array,array0,cnv)

c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * ieq    -> equation identifier
c       * ivar   -> vector component
c       * nvar   -> vector dimension
c       * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c       * loc    -> location in dimension (0 -> right, 1 -> left)
c       * bctype -> type of BC (dirichlet, neumann, periodic, etc.)
c       * array  -> real array with ghost-nodes
c       * array0 -> auxiliary real array
c       * cnv    -> whether vector components are contravariant (optional)
c     ---------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: ieq,dim,loc,bctype,nvar,ivar
      real(8)    :: array (0:,0:,0:,:)
     .             ,array0(0:,0:,0:,:)
      logical,optional :: cnv

c     Local variables

      integer    :: neq,ibc,i,j,k,ig,jg,kg
      logical    :: cov

      real(8),allocatable,dimension(:,:,:,:) :: zeros

      INTERFACE
        subroutine scalarDirichletBC(g_def,arr,arr0,ieq,dim,loc,order)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: ieq,dim,loc,order
          real(8) :: arr (0:,0:,0:)
     .              ,arr0(0:,0:,0:)
        end subroutine scalarDirichletBC

        subroutine vectorDirichletBC(g_def,i,arr,arr0,ieq,dim,loc,order)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: ieq,dim,loc,order,i
          real(8) :: arr (0:,0:,0:,:)
     .              ,arr0(0:,0:,0:,:)
        end subroutine vectorDirichletBC

        subroutine scalarNeumannBC(g_def,arr,arr0,ieq,dim,loc)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: ieq,dim,loc
          real(8) :: arr (0:,0:,0:)
     .              ,arr0(0:,0:,0:)
        end subroutine scalarNeumannBC

        subroutine vectorNeumannBC(g_def,i,arr,arr0,ieq,dim,loc)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: ieq,dim,loc,i
          real(8) :: arr (0:,0:,0:,:)
     .              ,arr0(0:,0:,0:,:)
        end subroutine vectorNeumannBC
      END INTERFACE

c     Begin program

      if (PRESENT(cnv)) then
        cov = .not.cnv
      else
        cov = .false.
      endif

      ibc = (1+loc)+2*(dim-1)
      
c     Determine boundary limits

      if (.not.selectFace(gbc_def,ibc)) return

c     Find BC update

      if (nvar == 1) then
        select case(bctype)
        case(PER)
          call periodicBC(gbc_def,array(:,:,:,nvar),dim,loc)
        case(EQU)
          call scalarDirichletBC(gbc_def,array(:,:,:,nvar)
     .                          ,array0(:,:,:,nvar)
     .                          ,ieq,dim,loc,order)
        case(DIR)
          allocate(zeros(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar))
          zeros = 0d0
          call scalarDirichletBC(gbc_def,array(:,:,:,nvar)
     .                          ,zeros(:,:,:,nvar)
     .                          ,ieq,dim,loc,order)
          deallocate(zeros)
        case(NEU)
          call scalarNeumannBC(gbc_def,array (:,:,:,nvar)
     .                        ,array0(:,:,:,nvar)
     .                        ,ieq,dim,loc)
        case(EXT)
          call fillFaces(gbc_def,array(:,:,:,nvar),dim,loc,order)
        case(SP)
          call scalarSingularBC(gbc_def,array(:,:,:,nvar))
        case(IFC)
          call interpolateToFace(gbc_def,array(:,:,:,nvar),dim,loc
     .                          ,order)
        case(SYNC) 
          call syncBC(gbc_def,array(:,:,:,nvar),dim,loc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC '//trim(int2char(bctype))//' not implemented')
        end select
      else
        select case(bctype)
        case(PER)
          call periodicBC(gbc_def,array(:,:,:,ivar),dim,loc)
        case(EQU)
          call vectorDirichletBC(gbc_def,ivar,array,array0,ieq,dim,loc
     .                          ,order)
        case(DIR)
          allocate(zeros(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar))
          zeros = 0d0
          call vectorDirichletBC(gbc_def,ivar,array,zeros,ieq,dim,loc
     .                          ,order)
          deallocate(zeros)
        case(NEU)
          call vectorNeumannBC(gbc_def,ivar,array,array0,ieq,dim,loc)
        case(EXT)
          call fillFaces(gbc_def,array(:,:,:,ivar),dim,loc,order)
        case(SP)
          if (get_stg_vector()) then
           call vectorSingularBC_stg(gbc_def,ivar,array(:,:,:,ivar),cov)
          else
           call vectorSingularBC    (gbc_def,ivar,array(:,:,:,ivar),cov)
          endif
        case(IFC)
          call interpolateToFace(gbc_def,array(:,:,:,ivar),dim,loc
     .                          ,order)
        case(SYNC) 
          call syncBC(gbc_def,array(:,:,:,ivar),dim,loc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC '//trim(int2char(bctype))//' not implemented')
        end select
      endif

c     End

      end subroutine fillGhostNodes_vec

c     syncBC
c     #################################################################
      subroutine syncBC(gbc_def,array,dim,loc)
c     -----------------------------------------------------------------
c     Synchronizes ghostcells in parallel. On input:
c        * array -> 3D array to impose BC's on.
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1 -> right)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer :: dim,loc
      real(8),dimension(0:,0:,0:) :: array

c     Local variables

      integer :: ibc,flg,nullp,currp
      integer,dimension(-1:1,-1:1,-1:1) :: pn

#if defined(petsc)
      integer :: dest,tag=0,ireq,request(10),status(MPI_STATUS_SIZE,10)

      real(8),volatile :: dumzin1 (nzbc+2)
     .                   ,dumzin2 (nzbc+2)
     .                   ,dumyin1 (nybc+2)
     .                   ,dumyin2 (nybc+2)
     .                   ,dumxin1 (nxbc+2)
     .                   ,dumxin2 (nxbc+2)

      real(8),volatile :: dumzout1(nzbc+2)
     .                   ,dumzout2(nzbc+2)
     .                   ,dumyout1(nybc+2)
     .                   ,dumyout2(nybc+2)
     .                   ,dumxout1(nxbc+2)
     .                   ,dumxout2(nxbc+2)
#endif

c     Begin program

      ibc = (1+loc)+2*(dim-1)

#if !defined(petsc)
      !Synchronize PER corners
      select case(dim)
      case(1)
        if (bcPER(2)) then
          array(0   ,0     ,:) = array(0   ,nybc,:)
          array(0   ,nybc+1,:) = array(0   ,1   ,:)
          array(nxbc,0     ,:) = array(nxbc,nybc,:)
          array(nxbc,nybc+1,:) = array(nxbc,1   ,:)
        endif
        if (bcPER(3)) then
          array(0   ,:,0     ) = array(0   ,:,nzbc)
          array(0   ,:,nzbc+1) = array(0   ,:,1   )
          array(nxbc,:,0     ) = array(nxbc,:,nzbc)
          array(nxbc,:,nzbc+1) = array(nxbc,:,1   )
        endif
      case(2)
        if (bcPER(1)) then
          array(0     ,0   ,:) = array(nxbc,0   ,:)
          array(nxbc+1,0   ,:) = array(1   ,0   ,:)
          array(0     ,nybc,:) = array(nxbc,nybc,:)
          array(nxbc+1,nybc,:) = array(1   ,nybc,:)
        endif
        if (bcPER(3)) then
          array(:,0   ,0     ) = array(:,0   ,nzbc)
          array(:,0   ,nzbc+1) = array(:,0   ,1   )
          array(:,nybc,0     ) = array(:,nybc,nzbc)
          array(:,nybc,nzbc+1) = array(:,nybc,1   )
        endif
      case(3)
        if (bcPER(1)) then
          array(0     ,:,0   ) = array(nxbc,:,0   )
          array(nxbc+1,:,0   ) = array(1   ,:,0   )
          array(0     ,:,nzbc) = array(nxbc,:,nzbc)
          array(nxbc+1,:,nzbc) = array(1   ,:,nzbc)
        endif
        if (bcPER(2)) then
          array(:,0     ,0   ) = array(:,nybc,0   )
          array(:,nybc+1,0   ) = array(:,1   ,0   )
          array(:,0     ,nzbc) = array(:,nybc,nzbc)
          array(:,nybc+1,nzbc) = array(:,1   ,nzbc)
        endif
      end select

#else

      dumxin1 = 0d0 ; dumxout1 = 0d0
      dumyin1 = 0d0 ; dumyout1 = 0d0
      dumzin1 = 0d0 ; dumzout1 = 0d0
      dumxin2 = 0d0 ; dumxout2 = 0d0
      dumyin2 = 0d0 ; dumyout2 = 0d0
      dumzin2 = 0d0 ; dumzout2 = 0d0
      
      pn = gbc_def%dactx(igxbc)%proc_neighbors

      nullp = MPI_PROC_NULL
      currp = pn(0,0,0)

      ireq = 0

      !Schedule receives first (flg=0) and sends later (flg=1)
      do flg=0,1
        select case(loc)
        case(0)

          if     (dim==1.and.isBdry(gbc_def,1,igxbc,ibc)) then
            call d1(0,flg)
          elseif (dim==2.and.isBdry(gbc_def,1,igybc,ibc)) then
            call d2(0,flg)
          elseif (dim==3.and.isBdry(gbc_def,1,igzbc,ibc)) then
            call d3(0,flg)
          endif
          
        case(1)

          if     (dim==1.and.isBdry(gbc_def,nxbc,igxbc,ibc)) then
            call d1(nxbc+1,flg)
          elseif (dim==2.and.isBdry(gbc_def,nybc,igybc,ibc)) then
            call d2(nybc+1,flg)
          elseif (dim==3.and.isBdry(gbc_def,nzbc,igzbc,ibc)) then
            call d3(nzbc+1,flg)
          endif

        end select

c$$$        if (flg == 0) then
c$$$          write (*,*) my_rank,dim,loc,'recv',ireq
c$$$        else
c$$$          write (*,*) my_rank,dim,loc,'send',ireq
c$$$        endif

      enddo

c$$$      write (*,*) my_rank,'Number of MPI requests',ireq

c     End program

      contains

c     d1
c     ################################################################
      subroutine d1(i,flg)

      implicit none

c     Call variables
      
      integer :: flg,i

c     Local variables

      logical :: down,up,top,bottom

c     Begin program

      down = .false. ; up = .false. ; top = .false. ; bottom = .false.
      
      !Down
      dest = pn(0,-1,0)
      down = (dest /= currp).and.(dest /= nullp)
      if (down) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumzout1,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumzin1 = array(i,1,0:nzbc+1)
          call MPI_Isend(dumzin1,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Up
      dest = pn(0,1,0)
      up = (dest /= currp).and.(dest /= nullp)
      if (up) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumzout2,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumzin2 = array(i,nybc,0:nzbc+1)
          call MPI_Isend(dumzin2,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Bottom
      dest = pn(0,0,-1)
      bottom = (dest /= currp).and.(dest /= nullp)
      if (bottom) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumyout1,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumyin1 = array(i,0:nybc+1,1)
          call MPI_Isend(dumyin1,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Top
      dest = pn(0,0,1)
      top = (dest /= currp).and.(dest /= nullp)
      if (top) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumyout2,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumyin2 = array(i,0:nybc+1,nzbc)
          call MPI_Isend(dumyin2,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      if(ireq>0.and.flg==1) then
        call MPI_Waitall(ireq,request,status,mpierr)

c$$$        if (pn( 0,-1, 0)/=nullp.and.pn( 0,-1, 0)/=currp)
c$$$     .       array(i,0     ,0:nzbc+1) = dumzout1
c$$$        if (pn( 0, 1, 0)/=nullp.and.pn( 0, 1, 0)/=currp)
c$$$     .       array(i,nybc+1,0:nzbc+1) = dumzout2
c$$$        if (pn( 0, 0,-1)/=nullp.and.pn( 0, 0,-1)/=currp)
c$$$     .       array(i,0:nybc+1,0     ) = dumyout1
c$$$        if (pn( 0, 0, 1)/=nullp.and.pn( 0, 0, 1)/=currp)
c$$$     .       array(i,0:nybc+1,nzbc+1) = dumyout2

        if (down) then
          if (pn(0,-1,0)==pn(0,1,0)) then  !Special case where both neighbors are same proc (PER)
            array(i,0     ,0:nzbc+1) = dumzout2
          else
            array(i,0     ,0:nzbc+1) = dumzout1
          endif
        endif

        if (up) then
          if (pn(0,-1,0)==pn(0,1,0)) then
            array(i,nybc+1,0:nzbc+1) = dumzout1
          else
            array(i,nybc+1,0:nzbc+1) = dumzout2
          endif
        endif

        if (bottom) then
          if (pn(0,0,-1)==pn(0,0,1)) then
            array(i,0:nybc+1,0     ) = dumyout2
          else
            array(i,0:nybc+1,0     ) = dumyout1
          endif
        endif

        if (top) then
          if (pn(0,0,-1)==pn(0,0,1)) then
            array(i,0:nybc+1,nzbc+1) = dumyout1
          else
            array(i,0:nybc+1,nzbc+1) = dumyout2
          endif
        endif

      endif

      end subroutine d1

c     d2
c     ################################################################
      subroutine d2(j,flg)

      implicit none

c     Call variables
      
      integer :: flg,j

c     Local variables
      
      logical :: left,right,top,bottom

c     Begin program

      left = .false. ; right = .false. ; top = .false. ; bottom =.false.

      !Left
      dest = pn(-1,0,0)
      left = (dest /= currp).and.(dest /= nullp)
      if (left) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumzout1,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumzin1 = array(1,j,0:nzbc+1)
          call MPI_Isend(dumzin1,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Right
      dest = pn(1,0,0)
      right = (dest /= currp).and.(dest /= nullp)
      if (right) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumzout2,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumzin2 = array(nxbc,j,0:nzbc+1)
          call MPI_Isend(dumzin2,(nzbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Bottom
      dest = pn(0,0,-1)
      bottom = (dest /= currp).and.(dest /= nullp)
      if (bottom) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumxout1,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumxin1 = array(0:nxbc+1,j,1)
          call MPI_Isend(dumxin1,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Top
      dest = pn(0,0,1)
      top = (dest /= currp).and.(dest /= nullp)
      if (top) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumxout2,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumxin2 = array(0:nxbc+1,j,nzbc)
          call MPI_Isend(dumxin2,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      if(ireq>0.and.flg==1) then
        call MPI_Waitall(ireq,request,status,mpierr)

c$$$        if (pn(-1, 0, 0)/=nullp.and.pn(-1, 0, 0)/=currp)
c$$$     .       array(0     ,j,0:nzbc+1) = dumzout1
c$$$        if (pn( 1, 0, 0)/=nullp.and.pn( 1, 0, 0)/=currp)
c$$$     .       array(nxbc+1,j,0:nzbc+1) = dumzout2
c$$$        if (pn( 0, 0,-1)/=nullp.and.pn( 0, 0,-1)/=currp)
c$$$     .       array(0:nxbc+1,j,0     ) = dumxout1
c$$$        if (pn( 0, 0, 1)/=nullp.and.pn( 0, 0, 1)/=currp)
c$$$     .       array(0:nxbc+1,j,nzbc+1) = dumxout2

        if (left) then
          if (pn(-1,0,0)==pn(1,0,0)) then
            array(0     ,j,0:nzbc+1) = dumzout2
          else
            array(0     ,j,0:nzbc+1) = dumzout1
          endif
        endif

        if (right) then
          if (pn(-1,0,0)==pn(1,0,0)) then
            array(nxbc+1,j,0:nzbc+1) = dumzout1
          else
            array(nxbc+1,j,0:nzbc+1) = dumzout2
          endif
        endif

        if (bottom) then
          if (pn(0,0,-1)==pn(0,0,1)) then
            array(0:nxbc+1,j,0     ) = dumxout2
          else
            array(0:nxbc+1,j,0     ) = dumxout1
          endif
        endif

        if (top) then
          if (pn(0,0,-1)==pn(0,0,1)) then
            array(0:nxbc+1,j,nzbc+1) = dumxout1
          else
            array(0:nxbc+1,j,nzbc+1) = dumxout2
          endif
        endif
      endif

      end subroutine d2

c     d3
c     ################################################################
      subroutine d3(k,flg)

      implicit none

c     Call variables

      integer :: flg,k

c     Local variables
      
      logical :: left,right,down,up

c     Begin program

      left = .false. ; right = .false. ; down = .false. ; up = .false.

      !Left
      dest = pn(-1,0,0)
      left = (dest /= currp).and.(dest /= nullp)
      if (left) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumyout1,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumyin1 = array(1,0:nybc+1,k)
          call MPI_Isend(dumyin1,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Right
      dest = pn(1,0,0)
      right = (dest /= currp).and.(dest /= nullp)
      if (right) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumyout2,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumyin2 = array(nxbc,0:nybc+1,k)
          call MPI_Isend(dumyin2,(nybc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Down
      dest = pn(0,-1,0)
      down = (dest /= currp).and.(dest /= nullp)
      if (down) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumxout1,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumxin1 = array(0:nxbc+1,1,k)
          call MPI_Isend(dumxin1,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      !Up
      dest = pn(0,1,0)
      up = (dest /= currp).and.(dest /= nullp)
      if (up) then
        if (flg == 0) then
          ireq = ireq+1
          call MPI_Irecv(dumxout2,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        else
          ireq = ireq+1
          dumxin2 = array(0:nxbc+1,nybc,k)
          call MPI_Isend(dumxin2,(nxbc+2)
     .         ,MPI_DOUBLE_PRECISION,dest,ibc
     .         ,gbc_def%MPI_COMM,request(ireq)
     .         ,mpierr)
        endif
      endif

      if(ireq>0.and.flg==1) then
        call MPI_Waitall(ireq,request,status,mpierr)
 
c$$$        if (pn(-1, 0, 0)/=nullp.and.pn(-1, 0, 0)/=currp)
c$$$     .       array(0     ,0:nybc+1,k) = dumyout1
c$$$        if (pn( 1, 0, 0)/=nullp.and.pn( 1, 0, 0)/=currp)
c$$$     .       array(nxbc+1,0:nybc+1,k) = dumyout2
c$$$        if (pn( 0,-1, 0)/=nullp.and.pn( 0,-1, 0)/=currp)
c$$$     .       array(0:nxbc+1,0     ,k) = dumxout1
c$$$        if (pn( 0, 1, 0)/=nullp.and.pn( 0, 1, 0)/=currp)
c$$$     .       array(0:nxbc+1,nybc+1,k) = dumxout2

        if (left) then
          if (pn(-1,0,0)==pn(1,0,0)) then
            array(0     ,0:nybc+1,k) = dumyout2
          else
            array(0     ,0:nybc+1,k) = dumyout1
          endif
        endif

        if (right) then
          if (pn(-1,0,0)==pn(1,0,0)) then
            array(nxbc+1,0:nybc+1,k) = dumyout1
          else
            array(nxbc+1,0:nybc+1,k) = dumyout2
          endif
        endif

        if (down) then
          if (pn(0,-1,0)==pn(0,1,0)) then
            array(0:nxbc+1,0     ,k) = dumxout2
          else
            array(0:nxbc+1,0     ,k) = dumxout1
          endif
        endif

        if (up) then
          if (pn(0,-1,0)==pn(0,1,0)) then
            array(0:nxbc+1,nybc+1,k) = dumxout1
          else
            array(0:nxbc+1,nybc+1,k) = dumxout2
          endif
        endif

      endif
      
      end subroutine d3

#endif

      end subroutine syncBC

c     periodicBC
c     #################################################################
      subroutine periodicBC(gbc_def,array,dim,loc)
c     -----------------------------------------------------------------
c     Imposes periodic BC. On input:
c        * array -> 3D array to impose BC's on.
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer :: dim,loc
      real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer :: ibc,i,j,k

      integer :: np_l,my_rank_l,dest,tag=0

      real(8),pointer,dimension(:,:) :: rhs

c     Begin program

      ibc = (1+loc)+2*(dim-1)

#if !defined(petsc)

      call per_BC(ibc,array)

#else

      select case(dim)
      case(1)

        if (  gbc_def%dactx(igxbc)%proc_neighbors( 0,0,0)
     .      ==gbc_def%dactx(igxbc)%proc_neighbors( 1,0,0)) then

          call per_BC(ibc,array)

        elseif (loc==0.and.isBdry(gbc_def,1,igxbc,ibc)) then

          allocate(rhs(2,nzbc+2))
          !Left
          dest = gbc_def%dactx(igxbc)%proc_neighbors(-1,0,0)
          rhs(1,:) = array(1,0     ,:)
          rhs(2,:) = array(1,nybc+1,:)
          call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
     .                ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                ,gbc_def%MPI_COMM,MPI_STATUS_IGNORE
     .                ,mpierr)
          array(0,0     ,:) = rhs(1,:)
          array(0,nybc+1,:) = rhs(2,:)
          deallocate(rhs)
          
        elseif (loc==1.and.isBdry(gbc_def,nxbc,igxbc,ibc)) then

          allocate(rhs(2,nzbc+2))
          !Right
          dest = gbc_def%dactx(igxbc)%proc_neighbors(1,0,0)
          rhs(1,:) = array(nxbc,0     ,:)
          rhs(2,:) = array(nxbc,nybc+1,:)
          call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
     .                ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                ,gbc_def%MPI_COMM,MPI_STATUS_IGNORE
     .                ,mpierr)
          array(nxbc+1,0     ,:) = rhs(1,:)
          array(nxbc+1,nybc+1,:) = rhs(2,:)
          deallocate(rhs)
          
        endif
        
      case(2)

        if (  gbc_def%dactx(igybc)%proc_neighbors(0, 0,0)
     .      ==gbc_def%dactx(igybc)%proc_neighbors(0, 1,0)) then

          call per_BC(ibc,array)

        elseif (loc==0.and.isBdry(gbc_def,1,igybc,ibc)) then
            
          allocate(rhs(2,nzbc+2))
          !Down
          dest = gbc_def%dactx(igybc)%proc_neighbors(0,-1,0)
          rhs(1,:) = array(0     ,1,:)
          rhs(2,:) = array(nxbc+1,1,:)
          call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
     .               ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .               ,gbc_def%MPI_COMM,MPI_STATUS_IGNORE
     .               ,mpierr)
          array(0     ,0,:) = rhs(1,:)
          array(nxbc+1,0,:) = rhs(2,:)
          deallocate(rhs)

        elseif (loc==1.and.isBdry(gbc_def,nybc,igybc,ibc)) then

          allocate(rhs(2,nzbc+2))
          !Up
          dest = gbc_def%dactx(igybc)%proc_neighbors(0,1,0)
          rhs(1,:) = array(0     ,nybc,:)
          rhs(2,:) = array(nxbc+1,nybc,:)
          call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
     .               ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .               ,gbc_def%MPI_COMM,MPI_STATUS_IGNORE
     .               ,mpierr)
          array(0     ,nybc+1,:) = rhs(1,:)
          array(nxbc+1,nybc+1,:) = rhs(2,:)

          deallocate(rhs)

        endif

      case(3)

        if (  gbc_def%dactx(igzbc)%proc_neighbors(0,0, 0)
     .      ==gbc_def%dactx(igzbc)%proc_neighbors(0,0, 1)) then

          call per_BC(ibc,array)

        elseif (loc==0.and.isBdry(gbc_def,1,igzbc,ibc)) then

          allocate(rhs(2,nxbc+2))
          !Bottom
          dest = gbc_def%dactx(igzbc)%proc_neighbors(0,0,-1)
          rhs(1,:) = array(:,0     ,1)
          rhs(2,:) = array(:,nybc+1,1)
          call MPI_Sendrecv_replace(rhs,2*(nxbc+2)
     .                 ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                 ,gbc_def%MPI_COMM,MPI_STATUS_IGNORE
     .                 ,mpierr)
          array(:,0     ,0) = rhs(1,:) 
          array(:,nybc+1,0) = rhs(2,:)
          deallocate(rhs)
  
        elseif (loc==1.and.isBdry(gbc_def,nzbc,igzbc,ibc)) then

          allocate(rhs(2,nxbc+2))
          !Top
          dest = gbc_def%dactx(igzbc)%proc_neighbors(0,0,1)
          rhs(1,:) = array(:,0     ,nzbc)
          rhs(2,:) = array(:,nybc+1,nzbc)
          call MPI_Sendrecv_replace(rhs,2*(nxbc+2)
     .                 ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                 ,gbc_def%MPI_COMM,MPI_STATUS_IGNORE
     .                 ,mpierr)
          array(:,0     ,nzbc+1) = rhs(1,:) 
          array(:,nybc+1,nzbc+1) = rhs(2,:)
          deallocate(rhs)

        endif

      end select

c$$$      if (.not.associated(gbc_def)) then
c$$$        call pstop("periodicBC"
c$$$     .            ,"This routine should not be called standalone")
c$$$      endif
c$$$
c$$$      if (gbc_def%MPI_COMM_PER(dim) /= MPI_COMM_NULL) then
c$$$        call MPI_Comm_size(gbc_def%MPI_COMM_PER(dim),np_l     ,mpierr)
c$$$        call MPI_Comm_rank(gbc_def%MPI_COMM_PER(dim),my_rank_l,mpierr)
c$$$
c$$$        if (np_l > 2.or.np_l < 1) then
c$$$          call pstop('periodicBC'
c$$$     .         ,'Periodic comm in dimension '//trim(int2char(dim))
c$$$     .       //' has wrong rank='//trim(int2char(np_l)))
c$$$        elseif (np_l == 2) then
c$$$          dest = 1 - my_rank_l
c$$$        endif
c$$$
c$$$        if (np_l == 1) then
c$$$
c$$$          call per_BC(ibc,array)
c$$$
c$$$        else
c$$$
c$$$          select case(dim)
c$$$          case(1)
c$$$
c$$$            allocate(rhs(2,nzbc+2))
c$$$
c$$$            if (loc==0.and.isBdry(gbc_def,1,igxbc,ibc)) then
c$$$cc              rhs = array(1,:,:)
c$$$              rhs(1,:) = array(1,0     ,:)
c$$$              rhs(2,:) = array(1,nybc+1,:)
c$$$              call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
c$$$     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
c$$$     .                   ,gbc_def%MPI_COMM_PER(dim),MPI_STATUS_IGNORE
c$$$     .                   ,mpierr)
c$$$cc              rhs = array(0,:,:)
c$$$              array(0,0     ,:) = rhs(1,:)
c$$$              array(0,nybc+1,:) = rhs(2,:)
c$$$
c$$$            elseif (loc==1.and.isBdry(gbc_def,nxbc,igxbc,ibc)) then
c$$$cc              rhs = array(nxbc,:,:)
c$$$              rhs(1,:) = array(nxbc,0     ,:)
c$$$              rhs(2,:) = array(nxbc,nybc+1,:)
c$$$              call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
c$$$     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
c$$$     .                   ,gbc_def%MPI_COMM_PER(dim),MPI_STATUS_IGNORE
c$$$     .                   ,mpierr)
c$$$cc              array(nxbc+1,:,:) = rhs
c$$$              array(nxbc+1,0     ,:) = rhs(1,:)
c$$$              array(nxbc+1,nybc+1,:) = rhs(2,:)
c$$$            endif
c$$$
c$$$            deallocate(rhs)
c$$$
c$$$          case(2)
c$$$
c$$$            allocate(rhs(2,nzbc+2))
c$$$
c$$$            if (loc==0.and.isBdry(gbc_def,1,igybc,ibc)) then
c$$$cc              rhs = array(:,1,:)
c$$$              rhs(1,:) = array(0     ,1,:)
c$$$              rhs(2,:) = array(nxbc+1,1,:)
c$$$              call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
c$$$     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
c$$$     .                   ,gbc_def%MPI_COMM_PER(dim),MPI_STATUS_IGNORE
c$$$     .                   ,mpierr)
c$$$cc              array(:,0,:) = rhs
c$$$              array(0     ,0,:) = rhs(1,:)
c$$$              array(nxbc+1,0,:) = rhs(2,:)
c$$$            elseif (loc==1.and.isBdry(gbc_def,nybc,igybc,ibc)) then
c$$$cc              rhs = array(:,nybc,:)
c$$$              rhs(1,:) = array(0     ,nybc,:)
c$$$              rhs(2,:) = array(nxbc+1,nybc,:)
c$$$              call MPI_Sendrecv_replace(rhs,2*(nzbc+2)
c$$$     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
c$$$     .                   ,gbc_def%MPI_COMM_PER(dim),MPI_STATUS_IGNORE
c$$$     .                   ,mpierr)
c$$$cc              array(:,nybc+1,:) = rhs
c$$$              array(0     ,nybc+1,:) = rhs(1,:)
c$$$              array(nxbc+1,nybc+1,:) = rhs(2,:)
c$$$            endif
c$$$
c$$$            deallocate(rhs)
c$$$
c$$$          case(3)
c$$$
c$$$            allocate(rhs(2,nxbc+2))
c$$$
c$$$            if (loc==0.and.isBdry(gbc_def,1,igzbc,ibc)) then
c$$$cc              rhs = array(:,:,1)
c$$$              rhs(1,:) = array(:,0     ,1)
c$$$              rhs(2,:) = array(:,nybc+1,1)
c$$$              call MPI_Sendrecv_replace(rhs,2*(nxbc+2)
c$$$     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
c$$$     .                   ,gbc_def%MPI_COMM_PER(dim),MPI_STATUS_IGNORE
c$$$     .                   ,mpierr)
c$$$cc              array(:,:,0) = rhs
c$$$              array(:,0     ,0) = rhs(1,:) 
c$$$              array(:,nybc+1,0) = rhs(2,:)
c$$$            elseif (loc==1.and.isBdry(gbc_def,nzbc,igzbc,ibc)) then
c$$$cc              rhs = array(:,:,nzbc)
c$$$              rhs(1,:) = array(:,0     ,nzbc)
c$$$              rhs(2,:) = array(:,nybc+1,nzbc)
c$$$              call MPI_Sendrecv_replace(rhs,2*(nxbc+2)
c$$$     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
c$$$     .                   ,gbc_def%MPI_COMM_PER(dim),MPI_STATUS_IGNORE
c$$$     .                   ,mpierr)
c$$$cc              array(:,:,nzbc+1) = rhs
c$$$              array(:,0     ,nzbc+1) = rhs(1,:) 
c$$$              array(:,nybc+1,nzbc+1) = rhs(2,:)
c$$$            endif
c$$$
c$$$            deallocate(rhs)
c$$$          end select
c$$$        endif
c$$$
c$$$      endif
#endif

c     End program

      end subroutine periodicBC

c     extrapolateBC
c     #################################################################
      subroutine extrapolateBC(gbc_def,array,array0,dim,loc,order,atfce)
c     -----------------------------------------------------------------
c     Fills ghost nodes by extrapolation across relevant boundary.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gbc_def

        integer :: order,dim,loc
        real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .            ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

        logical, optional :: atfce

c     Local variables

        integer :: imax,imin,jmax,jmin,kmax,kmin,ibc
     .            ,order_int(3)

        integer :: i,j,k,ig,jg,kg,one,two,thr

        real(8) :: x0

        real(8),pointer,dimension(:) :: xx,yy,zz

        logical :: atface

c     Begin program

        ibc = (1+loc)+2*(dim-1)

        order_int(1) = min(order,nxbc)
        order_int(2) = min(order,nybc)
        order_int(3) = min(order,nzbc)

        one = 1 ; two = 2; thr = 3

        if (PRESENT(atfce)) then
          atface = atfce
        else
          atface = .true.
        endif

c     Setup grid aliases

        xx => gbc_def%xx
        yy => gbc_def%yy
        zz => gbc_def%zz

c     Extrapolate

        !Do not change: SAMRAI depends on it
        call findBCLoopLimits(gbc_def,dim,loc,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                       ,imin ,imax ,jmin ,jmax ,kmin ,kmax )

        if (.not.bcSYM(ibc)) then

          !Safeguards to prevent stepping out of domain
          if (order_int(dim) <= 1) then
            two = one
            thr = one
          elseif (order_int(dim) == 2) then
            thr = two
          endif

          select case (ibc)
          case (1)

            i = imin

            call getMGmap(gbc_def,i,1,1,igxbc,igybc,igzbc,ig,jg,kg)

            if (atface) then
              x0 = 0.5*(xx(ig+one)+xx(ig))
              array(i,:,:)= quad_int(x0,xx(ig+one),xx(ig+two),xx(ig+thr)
     .                      ,array0(i    ,:,:),array(i+one,:,:)
     .                      ,array (i+two,:,:),array(i+thr,:,:)
     .                      ,xx(ig),order_int(dim) )
            else
              array(i,:,:)= array0(i,:,:)
            endif

          case (2)

            i = imax

            call getMGmap(gbc_def,i,1,1,igxbc,igybc,igzbc,ig,jg,kg)

            if (atface) then
              x0 = 0.5*(xx(ig-one)+xx(ig))
              array(i,:,:)= quad_int(x0,xx(ig-one),xx(ig-two),xx(ig-thr)
     .                      ,array0(i    ,:,:),array(i-one,:,:)
     .                      ,array (i-two,:,:),array(i-thr,:,:)
     .                      ,xx(ig),order_int(dim) )
            else
              array(i,:,:)= array0(i,:,:)
            endif

          case (3)

            j = jmin

            call getMGmap(gbc_def,1,j,1,igxbc,igybc,igzbc,ig,jg,kg)

            if (atface) then
              x0 = 0.5*(yy(jg+one)+yy(jg))
              array(:,j,:)= quad_int(x0,yy(jg+one),yy(jg+two),yy(jg+thr)
     .                      ,array0(:,j    ,:),array(:,j+one,:)
     .                      ,array (:,j+two,:),array(:,j+thr,:)
     .                      ,yy(jg),order_int(dim) )
            else
              array(:,j,:)=array0(:,j,:)
            endif

          case (4)

            j = jmax

            call getMGmap(gbc_def,1,j,1,igxbc,igybc,igzbc,ig,jg,kg)

            if (atface) then
              x0 = 0.5*(yy(jg-one)+yy(jg))

              array(:,j,:)= quad_int(x0,yy(jg-one),yy(jg-two),yy(jg-thr)
     .                      ,array0(:,j    ,:),array(:,j-one,:)
     .                      ,array (:,j-two,:),array(:,j-thr,:)
     .                      ,yy(jg),order_int(dim) )
            else
              array(:,j,:)=array0(:,j,:)
            endif

          case (5)

            k = kmin

            call getMGmap(gbc_def,1,1,k,igxbc,igybc,igzbc,ig,jg,kg)

            if (atface) then
              x0 = 0.5*(zz(kg+one)+zz(kg))

              array(:,:,k)= quad_int(x0,zz(kg+one),zz(kg+two),zz(kg+thr)
     .                      ,array0(:,:,k    ),array(:,:,k+one)
     .                      ,array (:,:,k+two),array(:,:,k+thr)
     .                      ,zz(kg),order_int(dim) )
            else
              array(:,:,k)= array0(:,:,k)
            endif

          case (6)

            k = kmax

            call getMGmap(gbc_def,1,1,k,igxbc,igybc,igzbc,ig,jg,kg)

            if (atface) then
              x0 = 0.5*(zz(kg-one)+zz(kg))
 
              array(:,:,k)= quad_int(x0,zz(kg-one),zz(kg-two),zz(kg-thr)
     .                      ,array0(:,:,k    ),array(:,:,k-one)
     .                      ,array (:,:,k-two),array(:,:,k-thr)
     .                      ,zz(kg),order_int(dim) )
            else
             array(:,:,k)= array0(:,:,k)
            endif

          case default

            call pstop('extrapolateBC'
     .            ,'Boundary '//int2char(ibc)//' non existent')

          end select

        else  !Extrapolate for ODD BCs

          select case (ibc)
          case (1)
            i = imin
            array(i,:,:) = -array(i+two,:,:)
          case (2)
            i = imax
            array(i,:,:) = -array(i-two,:,:)
          case (3)
            j = jmin
            array(:,j,:) = -array(:,j+two,:)
          case (4)
            j = jmax
            array(:,j,:) = -array(:,j-two,:)
          case (5)
            k = kmin
            array(:,:,k) = -array(:,:,k+two)
          case (6)
            k = kmax
            array(:,:,k) = -array(:,:,k-two)
          case default
            call pstop('extrapolateBC'
     .            ,'Boundary '//int2char(ibc)//' non existent')
          end select

        endif

      end subroutine extrapolateBC

c     interpolateToFace
c     #################################################################
      subroutine interpolateToFace(gbc_def,array,dim,loc,order)
c     -----------------------------------------------------------------
c     Interpolates values to face for EQU BC. On input:
c        * array -> 3D array to impose BC's on.
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gbc_def

        integer :: dim,loc,order
        real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

        integer :: imax,imin,jmax,jmin,kmax,kmin,ibc
     .            ,order_int(3)

        integer :: i,j,k,ig,jg,kg,one,two,thr
        real(8),pointer,dimension(:) :: xx,yy,zz

cc        type(grid_mg_def),pointer :: gl_def

c     Begin program

        ibc = (1+loc)+2*(dim-1)

        order_int(1) = min(order,nxbc)
        order_int(2) = min(order,nybc)
        order_int(3) = min(order,nzbc)

        one = 1; two = 2 ; thr = 3

        !Safeguards to prevent stepping out of domain
        if (order_int(dim) <= 1) then
          two = one
          thr = one
        elseif (order_int(dim) == 2) then
          thr = two
        endif

c     Setup grid aliases

        xx => gbc_def%xx
        yy => gbc_def%yy
        zz => gbc_def%zz

c     Extrapolate

        !Do NOT change: SAMRAI depends on it
        call findBCLoopLimits(gbc_def,dim,loc,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                       ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

        select case (ibc)
        case (1)

          i = imin

          call getMGmap(gbc_def,i,1,1,igxbc,igybc,igzbc,ig,jg,kg)

          array(i,:,:) = quad_int(xx(ig    ),xx(ig+one)
     .                      ,xx(ig+two),xx(ig+thr)
     .                      ,array(i    ,:,:),array(i+one,:,:)
     .                      ,array(i+two,:,:),array(i+thr,:,:)
cc     .                      ,xx(ig)+dxh(ig),order )
     .                      ,0.5*(xx(ig)+xx(ig+one)),order )

        case (2)

          i = imax

          call getMGmap(gbc_def,i,1,1,igxbc,igybc,igzbc,ig,jg,kg)

          array(i,:,:) =
     .              quad_int(xx(ig    ),xx(ig-one)
     .                      ,xx(ig-two),xx(ig-thr)
     .                      ,array(i    ,:,:),array(i-one,:,:)
     .                      ,array(i-two,:,:),array(i-thr,:,:)
cc     .                      ,xx(ig)-dxh(ig),order )
     .                      ,0.5*(xx(ig)+xx(ig-one)),order )

        case (3)

          j = jmin

          call getMGmap(gbc_def,1,j,1,igxbc,igybc,igzbc,ig,jg,kg)

          array(:,j,:) = quad_int(yy(jg    ),yy(jg+one)
     .                      ,yy(jg+two),yy(jg+thr)
     .                      ,array(:,j    ,:),array(:,j+one,:)
     .                      ,array(:,j+two,:),array(:,j+thr,:)
cc     .                      ,yy(jg)+dyh(jg),order )
     .                      ,0.5*(yy(jg)+yy(jg+one)),order )

        case (4)

          j = jmax

          call getMGmap(gbc_def,1,j,1,igxbc,igybc,igzbc,ig,jg,kg)

          array(:,j,:) = quad_int(yy(jg    ),yy(jg-one)
     .                      ,yy(jg-two),yy(jg-thr)
     .                      ,array(:,j    ,:),array(:,j-one,:)
     .                      ,array(:,j-two,:),array(:,j-thr,:)
cc     .                      ,yy(jg)-dyh(jg),order )
     .                      ,0.5*(yy(jg)+yy(jg-one)),order )

        case (5)

          k = kmin

          call getMGmap(gbc_def,1,1,k,igxbc,igybc,igzbc,ig,jg,kg)

          array(:,:,k) = quad_int(zz(kg    ),zz(kg+one)
     .                      ,zz(kg+two),zz(kg+thr)
     .                      ,array(:,:,k    ),array(:,:,k+one)
     .                      ,array(:,:,k+two),array(:,:,k+thr)
cc     .                      ,zz(kg)+dzh(kg),order )
     .                      ,0.5*(zz(kg)+zz(kg+one)),order )

        case (6)

          k = kmax

          call getMGmap(gbc_def,1,1,k,igxbc,igybc,igzbc,ig,jg,kg)

          array(:,:,k) = quad_int(zz(kg    ),zz(kg-one)
     .                      ,zz(kg-two),zz(kg-thr)
     .                      ,array(:,:,k    ),array(:,:,k-one)
     .                      ,array(:,:,k-two),array(:,:,k-thr)
cc     .                      ,zz(kg)-dzh(kg),order )
     .                      ,0.5*(zz(kg)+zz(kg-one)),order )

        end select

c     End program

      end subroutine interpolateToFace

c     fillFaces
c     #################################################################
      subroutine fillFaces(gbc_def,array,dim,loc,order)
c     -----------------------------------------------------------------
c     Fills ghost cells in faces by extrapolation (6 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     This routine is only used when no good information is available
c     to fill the ghost cells of a particular quantity (for instance,
c     during MG restriction of arrays that contain EQU BCs).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: dim,loc,order
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: ibcl(3),type
      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,bcond(6)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      bcond = 0
      bcond(ibc) = EXT

      ibcl = 0

      type = 1
      
c     Determine boundary limits

      if (.not.selectFace(gbc_def,ibc)) return

      !Do NOT change: SAMRAI depends on it
      call findBCLoopLimits(gbc_def,dim,loc,iimin-1,iimax+1
     .                             ,jjmin-1,jjmax+1
     .                             ,kkmin-1,kkmax+1
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (bc_debug) then
        write(bc_debug_unit,111) '   fillFaces: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
 111    format (a,i4,a,i4,a,i4,a,6i4)
      endif

c     Select face

      select case (ibc)
      case (1) !On face #1, fill face i=0

        ito   = imin
        ifrom = min(ito+1+order,nxbc)

        imn   = min(ito,ifrom)
        imx   = max(ito,ifrom)

        ibcl(1) = 1

        do j=jmin,jmax
          do k=kmin,kmax
            call extrapolateToCorner(gbc_def,ibcl,ito,j,k,imn,j,k
     .                ,imx,j,k,array(imn:imx,j:j,k:k),bcond,order,type)
          enddo
        enddo

      case (2) !On face #2, fill face i=nx+1

        ito   = imax
        ifrom = max(ito-1-order,1)

        imn   = min(ito,ifrom)
        imx   = max(ito,ifrom)

        ibcl(1) = 2

        do j=jmin,jmax
          do k=kmin,kmax
            call extrapolateToCorner(gbc_def,ibcl,ito,j,k,imn,j,k
     .                ,imx,j,k,array(imn:imx,j:j,k:k),bcond,order,type)
          enddo
        enddo

      case (3) !On face #3, fill face j=0

        jto   = jmin
        jfrom = min(jto+1+order,nybc)

        jmn   = min(jto,jfrom)
        jmx   = max(jto,jfrom)

        ibcl(2) = 3

        do k=kmin,kmax
          do i=imin,imax
              call extrapolateToCorner(gbc_def,ibcl,i,jto,k,i,jmn,k
     .                ,i,jmx,k,array(i:i,jmn:jmx,k:k),bcond,order,type)
          enddo
        enddo

      case (4) !On face #4, fill face j=ny+1

        jto   = jmax
        jfrom = max(jto-1-order,1)

        jmn   = min(jto,jfrom)
        jmx   = max(jto,jfrom)

        ibcl(2) = 4

        do k=kmin,kmax
          do i=imin,imax
            call extrapolateToCorner(gbc_def,ibcl,i,jto,k,i,jmn,k
     .                ,i,jmx,k,array(i:i,jmn:jmx,k:k),bcond,order,type)
          enddo
        enddo

      case (5) !On face #5, fill face k=0

        kto   = kmin
        kfrom = min(kto+1+order,nzbc)

        kmn   = min(kto,kfrom)
        kmx   = max(kto,kfrom)

        ibcl(3) = 5

        do j=jmin,jmax
          do i=imin,imax
            call extrapolateToCorner(gbc_def,ibcl,i,j,kto,i,j,kmn
     .                ,i,j,kmx,array(i:i,j:j,kmn:kmx),bcond,order,type)
          enddo
        enddo

      case (6) !On face #6, fill face k=nz+1

        kto   = kmax
        kfrom = max(kto-1-order,1)

        kmn   = min(kto,kfrom)
        kmx   = max(kto,kfrom)

        ibcl(3) = 6

        do j=jmin,jmax
          do i=imin,imax
            call extrapolateToCorner(gbc_def,ibcl,i,j,kto,i,j,kmn
     .                ,i,j,kmx,array(i:i,j:j,kmn:kmx),bcond,order,type)
          enddo
        enddo

      end select

c     End program

      end subroutine fillFaces

c     fillEdges
c     #################################################################
      subroutine fillEdges(gbc_def,array,array0,bcond,order)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * bcond  -> boundary condition array
c        * order  -> order of extrapolation
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: bcond(6),order
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,ord
     .             ,dim,loc,ibcl(3),type

c     Begin program

      type = 2
      
      do dim=1,3
      do loc=0,1

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      call findBCLoopLimits(gbc_def,dim,loc
     .                             ,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      ord = max(order,1)+1  !Does not affect order of interpolation in extrapolateToCorner

      if (bc_debug) then
        write(bc_debug_unit,111) '>>> fillEdges: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
     .            ,'; Order =',order
 111    format (a,i4,a,i4,a,i4,a,6i4,a,i4)
      endif

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill edges y=jmin-1,y=jmax+1, vary z

        if (imin == 1) then
          ito   = imin-1
          ifrom = min(imin+ord,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1
     .       .and.(selectFace(gbc_def,1).and.selectFace(gbc_def,3)))then !Check faces #1,#3
            jto   = jmin-1
            jfrom = min(jmin+ord,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(3)==EQU) .and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/1,3,0/)
                call extrapolateToCorner(gbc_def,ibcl,ito,jto,k,imn,jmn
     .              ,k,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order
     .              ,type)
              endif
            enddo
          endif

          if (jmax == nybc
     .       .and.(selectFace(gbc_def,1).and.selectFace(gbc_def,4)))then !Check faces #1,#4
            jto   = jmax+1
            jfrom = max(jmax-ord,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/1,4,0/)
                call extrapolateToCorner(gbc_def,ibcl,ito,jto,k,imn,jmn
     .              ,k,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order
     .              ,type)
              endif
            enddo
          endif

        endif

      case (2) !On face #2, fill edges y=jmin-1,y=jmax+1, vary z

        if (imax == nxbc) then
          ito   = imax+1
          ifrom = max(imax-ord,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1
     .       .and.(selectFace(gbc_def,2).and.selectFace(gbc_def,3)))then !Check faces #2,#3
            jto   = jmin-1
            jfrom = min(jmin+ord,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(3)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/2,3,0/)
                call extrapolateToCorner(gbc_def,ibcl,ito,jto,k,imn,jmn
     .              ,k,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order
     .              ,type)
              endif
            enddo
          endif

          if (jmax == nybc
     .       .and.(selectFace(gbc_def,2).and.selectFace(gbc_def,4)))then !Check faces #2,#4
            jto   = jmax+1
            jfrom = max(jmax-ord,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/2,4,0/)
                call extrapolateToCorner(gbc_def,ibcl,ito,jto,k,imn,jmn
     .              ,k,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order
     .              ,type)
              endif
            enddo
          endif

        endif

      case (3) !On face #3, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+ord,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1
     .       .and.(selectFace(gbc_def,3).and.selectFace(gbc_def,5)))then !Check faces #3,#5

            kto   = kmin-1
            kfrom = min(kmin+ord,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,3,5/)
                call extrapolateToCorner(gbc_def,ibcl,i,jto,kto,i,jmn
     .            ,kmn,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

          if (kmax == nzbc
     .       .and.(selectFace(gbc_def,3).and.selectFace(gbc_def,6)))then !Check faces #3,#6

            kto   = kmax+1
            kfrom = max(kmax-ord,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,3,6/)
                call extrapolateToCorner(gbc_def,ibcl,i,jto,kto,i,jmn
     .            ,kmn,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

        endif

      case (4) !On face #4, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-ord,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1
     .       .and.(selectFace(gbc_def,4).and.selectFace(gbc_def,5)))then !Check faces #4,#5

            kto   = kmin-1
            kfrom = min(kmin+ord,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,4,5/)
                call extrapolateToCorner(gbc_def,ibcl,i,jto,kto,i,jmn
     .            ,kmn,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

          if (kmax == nzbc
     .       .and.(selectFace(gbc_def,4).and.selectFace(gbc_def,6)))then !Check faces #4,#6

            kto   = kmax+1
            kfrom = max(kmax-ord,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,4,6/)
                call extrapolateToCorner(gbc_def,ibcl,i,jto,kto,i,jmn
     .            ,kmn,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

        endif

      case (5) !On face #5, fill edges x=imin-1,x=imax+1, vary y

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+ord,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1
     .       .and.(selectFace(gbc_def,5).and.selectFace(gbc_def,1)))then !Check faces #5,#1

            ito   = imin-1
            ifrom = min(imin+ord,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/1,0,5/)
                call extrapolateToCorner(gbc_def,ibcl,ito,j,kto,imn,j
     .            ,kmn,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

          if (imax == nxbc
     .       .and.(selectFace(gbc_def,5).and.selectFace(gbc_def,2)))then !Check faces #5,#2

            ito   = imax+1
            ifrom = max(imax-ord,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/2,0,5/)
                call extrapolateToCorner(gbc_def,ibcl,ito,j,kto,imn,j
     .            ,kmn,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

        endif

      case (6) !On face #6, fill edges x=imin-1,x=imax+1, vary y

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-ord,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1
     .       .and.(selectFace(gbc_def,6).and.selectFace(gbc_def,1)))then !Check faces #6,#1

            ito   = imin-1
            ifrom = min(imin+ord,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/1,0,6/)
                call extrapolateToCorner(gbc_def,ibcl,ito,j,kto,imn,j
     .            ,kmn,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

          if (imax == nxbc
     .       .and.(selectFace(gbc_def,6).and.selectFace(gbc_def,2)))then !Check faces #6,#2

            ito   = imax+1
            ifrom = max(imax-ord,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/2,0,6/)
                call extrapolateToCorner(gbc_def,ibcl,ito,j,kto,imn,j
     .            ,kmn,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order
     .            ,type)
              endif
            enddo
          endif

        endif

      end select

      enddo
      enddo

c     End program

      end subroutine fillEdges

c     fillCorners
c     #################################################################
      subroutine fillCorners(gbc_def,array,bcond,order)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer    :: bcond(6),order
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

c     Begin program

      call fillCorner(0     ,0     ,0     ,1,3,5) !Corner 1

      call fillCorner(0     ,nybc+1,0     ,1,4,5) !Corner 2

      call fillCorner(nxbc+1,0     ,0     ,2,3,5) !Corner 3

      call fillCorner(nxbc+1,nybc+1,0     ,2,4,5) !Corner 4

      call fillCorner(0     ,0     ,nzbc+1,1,3,6) !Corner 5

      call fillCorner(0     ,nybc+1,nzbc+1,1,4,6) !Corner 6

      call fillCorner(nxbc+1,0     ,nzbc+1,2,3,6) !Corner 7

      call fillCorner(nxbc+1,nybc+1,nzbc+1,2,4,6) !Corner 8

c     End program

      contains

c     fillCorner
c     #####################################################################
      subroutine fillCorner(icor,jcor,kcor,f1,f2,f3)

      implicit none

c     Call variables

      integer, INTENT(IN) :: icor,jcor,kcor,f1,f2,f3

c     Local variables

      integer :: ito,jto,kto,ifrom,jfrom,kfrom
     .          ,imn,imx,jmn,jmx,kmn,kmx,ord,ibcl(3),type

      logical :: sel_f(3)

c     Begin program

      type = 3
      
      ord = max(order,1)+1   !Does not affect order of interpolation in extrapolateToCorner

      ito = icor
      if (icor == 0) then
        ifrom = min(1+ord,nxbc)
      elseif (icor == nxbc+1) then
        ifrom = max(nxbc-ord,1)
      endif

      jto = jcor
      if (jcor == 0) then
        jfrom = min(1+ord,nybc)
      elseif (jcor == nybc+1) then
        jfrom = max(nybc-ord,1)
      endif

      kto = kcor
      if (kcor == 0) then
        kfrom = min(1+ord,nzbc)
      elseif (kcor == nzbc+1) then
        kfrom = max(nzbc-ord,1)
      endif

      if (( selectFace(gbc_def,f1)
     . .and.selectFace(gbc_def,f2)
     . .and.selectFace(gbc_def,f3))) then

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)

        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if ( bc_debug ) then
          write (bc_debug_unit,111) '>>> fillCorners: Proc',my_rank
     .                 ,'; Corner =',ito,jto,kto
     .                 ,'; From =',ifrom,jfrom,kfrom
     .                 ,'; Order =',order

 111      format (a,i4,a,3i4,a,3i4,a,i4)
        endif

        ibcl  = (/ f1, f2, f3 /)

        call extrapolateToCorner(gbc_def,ibcl,ito,jto,kto
     .                   ,imn,jmn,kmn,imx,jmx,kmx
     .             ,array(imn:imx,jmn:jmx,kmn:kmx),bcond,order,type)

      endif

      end subroutine fillCorner

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(gbc_def,ibcl,ito,jto,kto
     .            ,imn,jmn,kmn,imx,jmx,kmx,array,bcnd,order,type,array0)
c     -----------------------------------------------------------------
c     Routine that extrapolates information to edge or vertex corners.
c     It is the workhorse of the fillEdge and fillCorner routines.
c     In the call sequence:
c        * ibcl(3): boundary identifiers of converging faces (integer)
c        * ito,jto,kto: coordinates of ghost cell to be filled (integer)
c        * imn,jmn,kmn: node coordinates identifying start of
c                       interpolation region (integer)
c        * imx,jmx,kmx: node coordinates identifying end of interpolation
c                       region (integer)
c        * array: local array defined on interpolation region, including
c                 ghost cell (ito,jto,kto) (real)
c        * bcnd: array specifying BCs
c        * order: order of extrapolation
c        * type: 1 -> face, 2-> edge, 3-> corner
c        * array0, optional: array specifying EQU BCs
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: gbc_def

      integer :: ibcl(3),ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx,bcnd(6)
     .          ,order,type
      real(8) :: array(imn:imx,jmn:jmx,kmn:kmx)
      real(8),optional,intent(IN) :: array0(imn:imx,jmn:jmx,kmn:kmx)
      
c     Local variables

      integer :: i,j,k,ig,jg,kg,inx,iny,inz,order_int(3),offset,count_bc
     .          ,imni,imxi,jmni,jmxi,kmni,kmxi,i_bc,i1(3),i2(3)

      real(8) :: valx(1),valy(1),valz(1),xto(1),yto(1),zto(1),curv(3)
     $          ,idx,idy,idz,num

      logical :: flg(3),select_bc(3),face,edge,corner!,get_out

c     Extrapolation

      real(8),dimension(10) :: xx,yy,zz

c     Begin program
      
cc      xx = 0d0 ; yy = 0d0 ; zz = 0d0

c     Initialize variables

c$$$      count_bc = count(ibcl /= 0)
c$$$      face   = (count_bc == 1)
c$$$      edge   = (count_bc == 2)
c$$$      corner = (count_bc == 3)
      face   = (type == 1)
      edge   = (type == 2)
      corner = (type == 3)

c     Check for exceptions

cc      get_out = (.not.force_corners).and.checkExceptions()
      if (checkExceptions()) return

c     Identify directions for extrapolation.

      if (bc_debug) then
        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' Identifying directions extrap'
      endif

      inx = imx-imn
      iny = jmx-jmn
      inz = kmx-kmn

      !Order in each direction
      !Directions 1,2,3 are defined as normal to faces.
      order_int(1) = min(order,max(inx-1,0))
      order_int(2) = min(order,max(iny-1,0))
      order_int(3) = min(order,max(inz-1,0))

      !Smoothness (curvature) indicator
      curv = 1d10

      if (order >= 1) then
        call getMGmap(gbc_def,ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)

        if (inx > 2) then
          if (ito == imn) then
            idx = 1d0/gbc_def%dx(ig+2)
            curv(1) = (  array(ito+1,jto,kto)
     $                -2*array(ito+2,jto,kto)
     $                  +array(ito+3,jto,kto))*idx*idx
          else
            idx = 1d0/gbc_def%dx(ig-2)
            curv(1) = (  array(ito-1,jto,kto)
     $                -2*array(ito-2,jto,kto)
     $                  +array(ito-3,jto,kto))*idx*idx
          endif
        endif

        if (iny > 2) then
          if (jto == jmn) then
            idy = 1d0/gbc_def%dy(jg+2)
            curv(2) = (  array(ito,jto+1,kto)
     $                -2*array(ito,jto+2,kto)
     $                  +array(ito,jto+3,kto))*idy*idy
          else                                              
            idy = 1d0/gbc_def%dy(jg-2)
            curv(2) = (  array(ito,jto-1,kto)
     $                -2*array(ito,jto-2,kto)
     $                  +array(ito,jto-3,kto))*idy*idy
          endif
        endif

        if (inz > 2) then
          if (kto == kmn) then
            idz = 1d0/gbc_def%dz(kg+2)
            curv(3) = (  array(ito,jto,kto+1)
     $                -2*array(ito,jto,kto+2)
     $                  +array(ito,jto,kto+3))*idz*idz
          else                                                  
            idz = 1d0/gbc_def%dz(kg-2)
            curv(3) = (  array(ito,jto,kto-1)
     $                -2*array(ito,jto,kto-2)
     $                  +array(ito,jto,kto-3))*idz*idz
          endif
        endif

        curv(1) = abs(curv(1))
        curv(2) = abs(curv(2))
        curv(3) = abs(curv(3))
      endif

      !Do not interpolate **along** direction if not enough points,
      !Directions 1,2,3 are defined as normal to faces.
      flg(1) = (inx > 0)
      flg(2) = (iny > 0)
      flg(3) = (inz > 0)

c     Treat special boundaries

      if (bc_debug) then
        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' Treating special boundaries'
      endif

      !Prepare cyclic directions
      i1(1) = 2 ; i2(1) = 3
      i1(2) = 3 ; i2(2) = 1
      i1(3) = 1 ; i2(3) = 2

      call select_along

      call select_across

c     Check smoothness if various interpolation directions are selected

      call select_smooth

c     Enforce symmetry extrapolation at corners

      do i=1,3
        i_bc = ibcl(i)
        if (i_bc <= 0) cycle
        select_bc(i)=(abs(bcond(i_bc))==FSYM.and.abs(bcnd(i_bc))==NEU)
     .               .and.selectFace(gbc_def,i_bc)
        if (select_bc(i)) order_int(i) = 0
      enddo

c     Check some direction for interpolation exists

      if (    (.not.flg(1)).and.(.not.flg(2)).and.(.not.flg(3))
     .    .or. bc_debug) then
        write (bc_debug_unit,*)
        write (bc_debug_unit,*) '>>> ExtrapolateToCorner, Proc=',my_rank
        write (bc_debug_unit,*) '    Extrapolation directions',flg
        write (bc_debug_unit,*) '    Boundary type: Face=',face
     .                         ,'      Edge=',edge,' Corner=',corner
        write (bc_debug_unit,*) '    Boundaries involved',ibcl
        write (bc_debug_unit,*) '    Order of extrapolation',order
        write (bc_debug_unit,*) '    Smoothness indicators',curv
        if (.not.bc_debug)
     .     call pstop('extrapolateToCorner'
     .            ,'No extrapolation directions left')
      endif

c     Initialize accumulators

      valx = 0d0
      valy = 0d0
      valz = 0d0

      num = 0d0

c     X extrapolation

      if (flg(1)) then

        if (ito == imn) then
         call getMGmap(gbc_def,imn+1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
         xx(1:inx) = gbc_def%xx(ig:ig+inx-1)
         imni = imn+1
         imxi = imx
        else
         call getMGmap(gbc_def,imx-1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
         xx(1:inx) = gbc_def%xx(ig-(inx-1):ig)
         imni = imn
         imxi = imx-1
        endif

        call getMGmap(gbc_def,ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        xto = gbc_def%xx(ig)

        call IntDriver1d(inx,xx(1:inx),array(imni:imxi,jto,kto)
     .                  ,1,xto,valx,order_int(1))

        num = num + 1d0

      endif

c     Y extrapolation

      if (flg(2)) then

        if (jto == jmn) then
         call getMGmap(gbc_def,imn,jmn+1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
         yy(1:iny) = gbc_def%yy(jg:jg+iny-1)
         jmni = jmn+1
         jmxi = jmx
        else
         call getMGmap(gbc_def,imn,jmx-1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
         yy(1:iny) = gbc_def%yy(jg-(iny-1):jg)
         jmni = jmn
         jmxi = jmx-1
        endif

        call getMGmap(gbc_def,ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        yto = gbc_def%yy(jg)

        call IntDriver1d(iny,yy(1:iny),array(ito,jmni:jmxi,kto)
     .                  ,1,yto,valy,order_int(2))

        num = num + 1d0

      endif

c     Z extrapolation

      if (flg(3)) then

        if (kto == kmn) then
         call getMGmap(gbc_def,imn,jmn,kmn+1,igxbc,igybc,igzbc,ig,jg,kg)
         zz(1:inz) = gbc_def%zz(kg:kg+inz-1)
         kmni = kmn+1
         kmxi = kmx
        else
         call getMGmap(gbc_def,imn,jmn,kmx-1,igxbc,igybc,igzbc,ig,jg,kg)
         zz(1:inz) = gbc_def%zz(kg-(inz-1):kg)
         kmni = kmn
         kmxi = kmx-1
        endif

        call getMGmap(gbc_def,ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        zto = gbc_def%zz(kg)

        call IntDriver1d(inz,zz(1:inz),array(ito,jto,kmni:kmxi)
     .                  ,1,zto,valz,order_int(3))

        num = num + 1d0

      endif

c     Average extrapolation results and fill ghost cell

      array(ito,jto,kto) = (valx(1) + valy(1) + valz(1))/num

      if (bc_debug) then
       if(face) write(bc_debug_unit,222)
     .                      ' >>> extrapolateToCorner - face, Proc='
     .                      ,my_rank,num,flg,ito,jto,kto,valx,valy,valz
       if(edge) write(bc_debug_unit,222)
     .                      ' >>> extrapolateToCorner - edge, Proc='
     .                      ,my_rank,num,flg,ito,jto,kto,valx,valy,valz
       if(corner)write(bc_debug_unit,222)
     .                      ' >>> extrapolateToCorner - corner,Proc='
     .                      ,my_rank,num,flg,ito,jto,kto,valx,valy,valz
 222   format(a,2i3,3l2,3i3,3e14.5)
      endif

c     End program

      contains

c     checkExceptions
c     ######################################################################
      function checkExceptions() result(get_out)

      implicit none

      logical :: get_out

      get_out = .false.

c     Check for "do nothing" conditions

      if (bc_debug) then
        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' Check for "do nothing" conditions'
        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' Skip corners',(skip_corners.and.(.not.face))
      endif

      get_out = (skip_corners.and.(.not.face))

      if (get_out) return
      
c     Do nothing for SP,DEF,PER, SYM (not FSYM) boundary faces, edges, and corners

      if (bc_debug) then
        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' SP boundaries',
     .       (bcChk(ibcl(1),SP,bcnd)
     .    .or.bcChk(ibcl(2),SP,bcnd)
     .    .or.bcChk(ibcl(3),SP,bcnd))

        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' DEF boundaries',
     .       (bcChk(ibcl(1),DEF,bcnd)
     .    .or.bcChk(ibcl(2),DEF,bcnd)
     .    .or.bcChk(ibcl(3),DEF,bcnd))

        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' SYM boundaries',
     .       (bcChk(ibcl(1),SYM,bcnd)
     .    .or.bcChk(ibcl(2),SYM,bcnd)
     .    .or.bcChk(ibcl(3),SYM,bcnd))

        write(bc_debug_unit,*)
     .       '>>> extrapolateToCorner, Proc='
     .      ,my_rank,' PER boundaries',
     .       (bcChk(ibcl(1),PER,bcnd)
     .    .or.bcChk(ibcl(2),PER,bcnd)
     .    .or.bcChk(ibcl(3),PER,bcnd))
      endif

      get_out = (bcChk(ibcl(1),PER,bcnd)
     .       .or.bcChk(ibcl(2),PER,bcnd)
     .       .or.bcChk(ibcl(3),PER,bcnd))

      if (get_out) return

      get_out = (bcChk(ibcl(1),SP,bcnd)
     .       .or.bcChk(ibcl(2),SP,bcnd)
     .       .or.bcChk(ibcl(3),SP,bcnd))

      if (get_out) return

      get_out = (bcChk(ibcl(1),SYM,bcnd)
     .       .or.bcChk(ibcl(2),SYM,bcnd)
     .       .or.bcChk(ibcl(3),SYM,bcnd))

      if (get_out) return
      
      get_out = (bcChk(ibcl(1),DEF,bcnd)
     .       .or.bcChk(ibcl(2),DEF,bcnd)
     .       .or.bcChk(ibcl(3),DEF,bcnd))

      if (get_out) return

c$$$c     Do nothing for doubly-periodic edges (PETSc does it)
c$$$
c$$$      if (edge) then
c$$$        if ((ibcl(1) == 0)
c$$$     .      .and.(bcChk(ibcl(2),PER).and.bcChk(ibcl(3),PER))) return
c$$$        if ((ibcl(2) == 0)
c$$$     .      .and.(bcChk(ibcl(3),PER).and.bcChk(ibcl(1),PER))) return
c$$$        if ((ibcl(3) == 0)
c$$$     .      .and.(bcChk(ibcl(1),PER).and.bcChk(ibcl(2),PER))) return
c$$$      endif

c     Do nothing for Dirichlet edges (extrapolateBC does it), except when BOTH are

      if (edge) then
        get_out = ((ibcl(1) == 0)
     .  .and.((selectFace(gbc_def,ibcl(3))
     .       .and.(bcChk(ibcl(3),DIR,bcnd)
     .         .or.bcChk(ibcl(3),EQU,bcnd)))
     .    .or.(selectFace(gbc_def,ibcl(2))
     .       .and.(bcChk(ibcl(2),DIR,bcnd)
     .         .or.bcChk(ibcl(2),EQU,bcnd))))
     .       )

        if (get_out) return

        get_out = ((ibcl(2) == 0)
     .  .and.((selectFace(gbc_def,ibcl(3))
     .       .and.(bcChk(ibcl(3),DIR,bcnd)
     .         .or.bcChk(ibcl(3),EQU,bcnd)))
     .    .or.(selectFace(gbc_def,ibcl(1))
     .       .and.(bcChk(ibcl(1),DIR,bcnd)
     .         .or.bcChk(ibcl(1),EQU,bcnd))))
     .       )

        if (get_out) return

        get_out = ((ibcl(3) == 0)
     .  .and.((selectFace(gbc_def,ibcl(1))
     .       .and.(bcChk(ibcl(1),DIR,bcnd)
     .         .or.bcChk(ibcl(1),EQU,bcnd)))
     .    .or.(selectFace(gbc_def,ibcl(2))
     .       .and.(bcChk(ibcl(2),DIR,bcnd)
     .         .or.bcChk(ibcl(2),EQU,bcnd))))
     .       )

        if (get_out) return

      endif

      end function checkExceptions

c     checkFaces
c     ######################################################################
      subroutine checkFaces(ibcl)

      implicit none

      integer :: ibcl(3)
      logical :: sel_f(3)

c     -----------------------------------------------------------------------
c     Checks whether faces in ibcl is at external boundary or not, and sets
c     ibcl < 0 for those that are not.
c     -----------------------------------------------------------------------

      sel_f = (/selectFace(gbc_def,ibcl(1))
     .         ,selectFace(gbc_def,ibcl(2))
     .         ,selectFace(gbc_def,ibcl(3))/)

      where (.not.sel_f) ibcl = -ibcl

      end subroutine checkFaces

c     select_along
c     ######################################################################
      subroutine select_along

      implicit none

c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs ALONG
c     selected boundary only.
c         1) Interpolate-to-face boundary
c         2) Equilibrium boundary (removed 12/13/06, added back 11/10/08)
c         3) Discard internal processor boundaries (not real boundaries)
c     -----------------------------------------------------------------------

      integer :: bcl
      
      select_bc = .false.
      do i=1,3
c$$$        i_bc = ibcl(i)
c$$$        if (i_bc <= 0) cycle
c$$$        select_bc(i) =(    (abs(bcnd(i_bc)) == IFC)
c$$$     .                 .or.(abs(bcnd(i_bc)) == EQU)
c$$$     .                ).and.selectFace(gbc_def,i_bc)
        if (ibcl(i) <= 0) cycle
        bcl = abs(bcnd(ibcl(i)))
        select_bc(i) =((bcl == IFC).or.(bcl == EQU))
     .                .and.selectFace(gbc_def,ibcl(i))
      enddo

      do i=1,3
        
        !Identify cases in which an edge has two or more directions
        !flagged for elimination:
        !  can't eliminate both interpolation directions simultaneously;
        !  resolve indetermination by choosing smoothest direction
        if (select_bc(i).and.(select_bc(i1(i)).or.select_bc(i2(i))))then

          if (select_bc(i1(i)).and.select_bc(i2(i))) then
            if (curv(i1(i)) <= curv(i2(i))) then
              select_bc(i1(i)) = .false.
            else
              select_bc(i2(i)) = .false.
            endif
          endif

          if (select_bc(i1(i))) then
            if (curv(i) <= curv(i1(i))) then
              cycle
            else
              flg(i) = .not.select_bc(i)
            endif
          else
            if (curv(i) <= curv(i2(i))) then
              cycle
            else
              flg(i) = .not.select_bc(i)
            endif
          endif

        elseif (select_bc(i)) then
          flg(i) = .not.select_bc(i)
        endif

      enddo

      end subroutine select_along

c     select_across
c     ######################################################################
      subroutine select_across
      implicit none
c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs ACROSS
c     selected boundary only.
c         1) Symmetry boundary (NEU+order=0; SYM removed on 6/1/12)
c         2) Discard internal processor boundaries (not real boundaries)
c     -----------------------------------------------------------------------

      integer :: bcl
      
      select_bc = .false.
      do i=1,3
c$$$        i_bc = ibcl(i)
c$$$        if (i_bc <= 0) cycle
c$$$        select_bc(i) =(     (abs(bcnd(i_bc)) == NEU .and. order == 0)
c$$$cc     .                 .or.(abs(bcond(i_bc))== FSYM .and. abs(bcnd(i_bc)) == NEU)
c$$$     .                ).and.selectFace(gbc_def,i_bc)
        if (ibcl(i) <= 0) cycle
        bcl = abs(bcnd(ibcl(i)))
        select_bc(i) =(   (bcl==NEU .and.order==0)
     .                .or.(bcl==FSYM.and.bcl==NEU))
     .                .and.selectFace(gbc_def,ibcl(i))
      enddo

      do i=1,3
        if (select_bc(i) .and. flg(i)) then
          flg(i1(i)) = .false.
          flg(i2(i)) = .false.
        endif
      enddo

      end subroutine select_across

c     select_smooth
c     ######################################################################
      subroutine select_smooth

      implicit none

c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs along
c     boundaries with comparable smoothness
c     -----------------------------------------------------------------------

      do i=1,3

        !Identify cases in which an edge has two or more directions
        !flagged for inclusion, and check for smoothness
        if (flg(i1(i)).and.flg(i2(i))) then
          if (curv(i1(i)) > 5*curv(i2(i))) then
            flg(i1(i)) = .false.
          elseif (curv(i2(i)) > 5*curv(i1(i))) then
            flg(i2(i)) = .false.
          endif
        endif

      enddo

      end subroutine select_smooth

      end subroutine extrapolateToCorner

c     XferBCs_scl
c     ###################################################################
      subroutine XferBCs_scl(arr1,arr2,average)

      implicit none

c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form arr1 to arr2
c     -------------------------------------------------------------------

c     Input variables

      real(8) :: arr1(0:,0:,0:)
      real(8) :: arr2(0:,0:,0:)
      logical :: average

c     Local variables

      integer :: nxp,nyp,nzp

c     Begin program

      nxp = size(arr1,1)-1
      nyp = size(arr1,2)-1
      nzp = size(arr1,3)-1

      if (average) then
        arr2 = arr1
        arr2(nxp,:,:) = 0.5*(arr2(nxp  ,:,:)
     .                      +arr2(nxp-1,:,:))
        arr2(0  ,:,:) = 0.5*(arr2(1    ,:,:)
     .                      +arr2(0    ,:,:))

        arr2(:,nyp,:) = 0.5*(arr2(:,nyp  ,:)
     .                      +arr2(:,nyp-1,:))
        arr2(:,0  ,:) = 0.5*(arr2(:,1    ,:)
     .                      +arr2(:,0    ,:))

        arr2(:,:,nzp) = 0.5*(arr2(:,:,nzp  )
     .                      +arr2(:,:,nzp-1))
        arr2(:,:,0  ) = 0.5*(arr2(:,:,1    )
     .                      +arr2(:,:,0    ))
      else
        !X boundaries
        arr2(0  ,:,:) = arr1(0  ,:,:)
        arr2(nxp,:,:) = arr1(nxp,:,:)

        !Y boundaries
        arr2(:,0  ,:) = arr1(:,0  ,:)
        arr2(:,nyp,:) = arr1(:,nyp,:)

        !Z boundaries
        arr2(:,:,0  ) = arr1(:,:,0  )
        arr2(:,:,nzp) = arr1(:,:,nzp)
      endif

c     End program

      end subroutine XferBCs_scl

c     XferBCs_vec
c     ###################################################################
      subroutine XferBCs_vec(vec1,vec2,average)

      implicit none
c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form vec1 to vec2
c     -------------------------------------------------------------------

c     Input variables

      real(8) :: vec1(0:,0:,0:,:)
      real(8) :: vec2(0:,0:,0:,:)
      logical :: average

c     Local variables

      integer :: nxp,nyp,nzp

c     Begin program

      nxp = size(vec1,1)-1
      nyp = size(vec1,2)-1
      nzp = size(vec1,3)-1

      if (average) then
        vec2 = vec1
        vec2(nxp,:,:,:) = 0.5*(vec2(nxp  ,:,:,:)
     .                        +vec2(nxp-1,:,:,:))
        vec2(0  ,:,:,:) = 0.5*(vec2(1    ,:,:,:)
     .                        +vec2(0    ,:,:,:))

        vec2(:,nyp,:,:) = 0.5*(vec2(:,nyp  ,:,:)
     .                        +vec2(:,nyp-1,:,:))
        vec2(:,0  ,:,:) = 0.5*(vec2(:,1    ,:,:)
     .                        +vec2(:,0    ,:,:))

        vec2(:,:,nzp,:) = 0.5*(vec2(:,:,nzp  ,:)
     .                        +vec2(:,:,nzp-1,:))
        vec2(:,:,0  ,:) = 0.5*(vec2(:,:,1    ,:)
     .                        +vec2(:,:,0    ,:))
      else
        !X boundaries
        vec2(0  ,:,:,:) = vec1(0  ,:,:,:)
        vec2(nxp,:,:,:) = vec1(nxp,:,:,:)

        !Y boundaries
        vec2(:,0  ,:,:) = vec1(:,0  ,:,:)
        vec2(:,nyp,:,:) = vec1(:,nyp,:,:)

        !Z boundaries
        vec2(:,:,0  ,:) = vec1(:,:,0  ,:)
        vec2(:,:,nzp,:) = vec1(:,:,nzp,:)
      endif

c     End program

      end subroutine XferBCs_vec

      end module imposeBCinterface

c module setMGBC_interface
c######################################################################
      module setMGBC_interface

        use grid_mg

        use imposeBCinterface

      contains

c     setMGBC
c     ###############################################################
      subroutine setMGBC(gbc_def,gpos,neq,nnx,nny,nnz,iig,array,bcnd
     .                  ,arr_cov,arr0,icomp,is_cnv,is_vec,is_curv
     .                  ,result_is_curv,out_iscnv,out_iscurv,iorder
     .                  ,regularize_sp)
c     ---------------------------------------------------------------
c     Interfaces BC routines with MG code, for preconditioning.
c     On input:
c       *gbc_def: grid definitions for MG
c       *gpos: node number (whole grid if gpos=0)
c       *neq,nnx,nny,nnz: number of degrees of freedom.
c       *iig: grid level
c       *array: array containing values of magnitudes for all grid
c               points. The magnitudes may be scalars or vector
c               components (the latter assumed contravariant by default)
c       *bcnd: array containing BC information
c
c     Optional parameters:
c       *arr_cov: covariant components of vector components.
c       *arr0: array containing dirichlet BCs (zero by default)
c       *icomp: integer array specifying magnitude of interest (if not
c               provided, use last used value).
c       *is_vec: logical variable specifying whether array
c          represents a curvilinear vector or a set of scalars
c          (if not provided, use last used value):
c             -If is_vec=.true., array is taken as curvilinear vector
c              components (covariant or contravariant, according
c              to variable is_cnv).
c             -If is_vec=.false., array is taken as scalars.
c       *is_curv: whether vector is provided in Cartesian (.false.) or
c                 curvilinear representation
c       *is_cnv: logical variable indicating whether array
c                contains covariant or contravariant components
c                (if not provided, use last used value).
c       *result_is_curv: logical variable specifying whether output
c                        in array should be curvilinear (.true.) or
c                        not (default=is_curv)
c       *out_iscurv (intent OUT): specifies whether output is curvilinear
c                                vector (needed when is_curv is inherited
c                                from a previous call).
c       *out_iscnv (intent OUT): specifies representation of output
c                                vector (needed when is_cnv is inherited
c                                from a previous call).
c       *iorder: order of inter/extrapolation in BC routines (quadratic
c                by default).
c       *regularize_sp: whether to regularize around SP
c
c     On output, array is returned either in the same curvilinear
c     representation or in Cartesian, if specified by "result_is_vec".
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nnx,nny,nnz,neq,bcnd(6,neq),iig,gpos

      real(8),target :: array(0:nnx+1,0:nny+1,0:nnz+1,neq)

      real(8),optional,intent(INOUT),target ::
     .                         arr_cov(0:nnx+1,0:nny+1,0:nnz+1,neq)
      real(8),optional,intent(IN),target ::
     .                         arr0   (0:nnx+1,0:nny+1,0:nnz+1,neq)

      integer,optional,intent(IN) :: icomp(max(neq,3)/3),iorder

      logical,optional,intent(IN) ::is_cnv,is_vec,is_curv,result_is_curv
     .                             ,regularize_sp
      logical,optional,intent(OUT)::out_iscnv,out_iscurv

      type(grid_mg_def),pointer :: gbc_def

c     Local variables

      integer :: i,j,k,stencil_width,ieq,neql
     .          ,imng,imxg,jmng,jmxg,kmng,kmxg
     .          ,imnl,imxl,jmnl,jmxl,kmnl,kmxl

      integer,save :: order=-1 !Needed to "remember" order of extrapolation at
                               !boundaries from previous calls when not explicitly
                               !set (exits with error when unset)

      integer,save :: ivar(2)  !Needed to "remember" ivar from previous calls
                               !when icomp is not provided.

      logical,save :: iscnv,iscurv !Needed to "remember' these from previous calls
                                   !when they are not provided.

      logical      :: res_is_curv,fpointers,vector

      real(8)      :: car(3),curv(3)

      real(8),pointer,dimension(:,:,:,:) :: v_cov,v_cnv,vv0

      character(3) :: cvec

      type(mg_ctx),pointer :: mgctx => null()

c     Begin program

      neql = neq !To avoid problems if neq is hardwired

c     Optional arguments

      if (PRESENT(icomp))  ivar(1:size(icomp))  = icomp

      if (PRESENT(iorder)) order = iorder

c     Consistency check

      if (     gbc_def%nxv(iig) /= nnx
     .    .or. gbc_def%nyv(iig) /= nny
     .    .or. gbc_def%nzv(iig) /= nnz) then
        write (*,*) gbc_def%nxv(iig), nnx
        write (*,*) gbc_def%nyv(iig), nny
        write (*,*) gbc_def%nzv(iig), nnz
        call pstop('setMGBC','Grid sizes do not agree')
      endif

      if (order < 0) then
        call pstop('setMGBC','Order of extrapolation not set')
      endif

c     Allocate MG pointers (if not allocated earlier)

      call allocMGPointers(neql,gbc_def,mgctx,fpointers)
      
c     Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(mgctx,gpos,iig,imnl,imxl,jmnl,jmxl,kmnl,kmxl)

c     Check LOCAL limits (return if not close to local domain boundaries)

      if (     (imnl > 1 .and. imxl < nnx)
     .    .and.(jmnl > 1 .and. jmxl < nny)
     .    .and.(kmnl > 1 .and. kmxl < nnz)) return

c     Find GLOBAL limits (required for BC routine setBC)

      imng = imnl + gbc_def%ilo(iig)-1
      imxg = imxl + gbc_def%ilo(iig)-1
      jmng = jmnl + gbc_def%jlo(iig)-1
      jmxg = jmxl + gbc_def%jlo(iig)-1
      kmng = kmnl + gbc_def%klo(iig)-1
      kmxg = kmxl + gbc_def%klo(iig)-1

c     Select operation

      if (PRESENT(is_vec)) then
        vector = is_vec
      else
        vector = (neql == 3 .or. neql == 6) !Backward compatibility with previous
                                           ! versions of this routine
      endif

      if (vector) then

        if (PRESENT(is_cnv )) iscnv = is_cnv

        if (PRESENT(is_curv)) then
          iscurv = is_curv
        else
          iscurv = .true.
        endif

        if (PRESENT(result_is_curv)) then
          res_is_curv = result_is_curv
        else
          res_is_curv = iscurv
        endif

        if (.not.iscurv) iscnv=.true.

        cvec = 'cov'
        if (iscnv) cvec = 'cnv'

        if (PRESENT(arr0)) then
          vv0 => arr0
        else
          allocate(vv0(0:nnx+1,0:nny+1,0:nnz+1,neql))
          vv0 = 0d0
        endif

        if (.not.iscurv) then  !Transform to curvilinear

          do ieq=1,neql/3
            call transformVector(gbc_def,iig,iig,iig
     .           ,imnl-1,imxl+1,jmnl-1,jmxl+1,kmnl-1,kmxl+1
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,1+3*(ieq-1))
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,2+3*(ieq-1))
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,3+3*(ieq-1))
     .                          ,'car',cvec)
          enddo

        endif

        if (PRESENT(arr_cov)) then
          v_cnv => array
          v_cov => arr_cov
        else
          if (iscnv) then
            v_cnv => array
            allocate(v_cov(0:nnx+1,0:nny+1,0:nnz+1,neql))
            v_cov = 0d0
          else
            v_cov => array
            allocate(v_cnv(0:nnx+1,0:nny+1,0:nnz+1,neql))
            v_cnv = 0d0
          endif
        endif

        if (iscnv) then
          do ieq=1,neql/3
            v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .          = XformVector(gbc_def,iig,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .                        ,.not.iscnv)
          enddo
        else
          do ieq=1,neql/3
            v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .          = XformVector(gbc_def,iig,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .                        ,.not.iscnv)
          enddo
        endif

        do ieq=1,neql/3
          call setBC(gbc_def,ivar(ieq),3,nnx,nny,nnz
     .            ,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .            ,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .            ,vv0  (:,:,:,1+3*(ieq-1):3*ieq)
     .            ,bcnd (:    ,1+3*(ieq-1):3*ieq)
     .            ,iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,is_cnv=iscnv,is_vec=.true.,iorder=order
     .            ,sp_reg=regularize_sp)
        enddo

        if (.not.PRESENT(arr0)) deallocate(vv0)

        if (.not.PRESENT(arr_cov)) then
          if (iscnv) then
            deallocate(v_cov)
          else
            deallocate(v_cnv)
          endif
        endif

        if (.not.res_is_curv) then  !Transform to cartesian

          do ieq=1,neql/3
            call transformVector(gbc_def,iig,iig,iig
     .           ,imnl-1,imxl+1,jmnl-1,jmxl+1,kmnl-1,kmxl+1
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,1+3*(ieq-1))
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,2+3*(ieq-1))
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,3+3*(ieq-1))
     .                          ,cvec,'car')
          enddo

        endif

        !For output
        if (PRESENT(out_iscnv )) out_iscnv =iscnv
        if (PRESENT(out_iscurv)) out_iscurv=iscurv

      else

        if (PRESENT(arr0)) then
          vv0 => arr0
        else
          allocate(vv0(0:nnx+1,0:nny+1,0:nnz+1,neql))
          vv0 = 0d0
        endif

        do ieq=1,neql
          call setBC(gbc_def,ivar(1),nnx,nny,nnz
     $              ,array(:,:,:,ieq),vv0(:,:,:,ieq)
     .              ,bcnd(:,ieq),iig,iig,iig
     .              ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .              ,iorder=order)
        enddo

        if (.not.PRESENT(arr0)) deallocate(vv0)

      endif

c     Deallocate MG pointers (if allocated in this routine)

      call deallocMGPointers(fpointers,mgctx)

c     End

      end subroutine setMGBC

c     setMGBC_stg
c     ###############################################################
      subroutine setMGBC_stg(gbc_def,gpos,neq,nnx,nny,nnz,iig,array,bcnd
     .                  ,arr_cov,arr0,icomp,is_cnv,is_vec,is_curv
     .                  ,result_is_curv,out_iscnv,out_iscurv,iorder)
c     ---------------------------------------------------------------
c     Interfaces staggered vector BC routines with MG code, for 
c     preconditioning.
c
c     On input:
c       *gpos: node number (whole grid if gpos=0)
c       *neq,nnx,nny,nnz: number of degrees of freedom.
c       *iig: grid level
c       *array: array containing values of magnitudes for all grid
c               points. The magnitudes may be scalars or vector
c               components (the latter assumed contravariant by default)
c       *bcnd: array containing BC information
c
c     Optional parameters:
c       *arr_cov: covariant components of vector components.
c       *arr0: array containing dirichlet BCs (zero by default)
c       *icomp: integer array specifying magnitude of interest (if not
c               provided, use last used value).
c       *is_vec: logical variable specifying whether array
c          represents a curvilinear vector or a set of scalars
c          (if not provided, use last used value):
c             -If is_vec=.true., array is taken as curvilinear vector
c              components (covariant or contravariant, according
c              to variable is_cnv).
c             -If is_vec=.false., array is taken as scalars.
c       *is_curv: whether vector is provided in Cartesian (.false.) or
c                 curvilinear representation
c       *is_cnv: logical variable indicating whether array
c                contains covariant or contravariant components
c                (if not provided, use last used value).
c       *result_is_curv: logical variable specifying whether output
c                        in array should be curvilinear (.true.) or
c                        not (default=is_curv)
c       *out_iscurv (intent OUT): specifies whether output is curvilinear
c                                vector (needed when is_curv is inherited
c                                from a previous call).
c       *out_iscnv (intent OUT): specifies representation of output
c                                vector (needed when is_cnv is inherited
c                                from a previous call).
c       *iorder: order of inter/extrapolation in BC routines (quadratic
c                by default).
c       *g_def: grid definitions for MG
c
c     On output, array is returned either in the same curvilinear
c     representation or in Cartesian, if specified by "result_is_vec".
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nnx,nny,nnz,neq,bcnd(6,neq),iig,gpos

      real(8),target :: array(0:nnx+1,0:nny+1,0:nnz+1,neq)

      real(8),optional,intent(INOUT),target ::
     .                         arr_cov(0:nnx+1,0:nny+1,0:nnz+1,neq)
      real(8),optional,intent(IN),target ::
     .                         arr0   (0:nnx+1,0:nny+1,0:nnz+1,neq)

      integer,optional,intent(IN) :: icomp(max(neq,3)/3),iorder

      logical,optional,intent(IN) ::is_cnv,is_vec,is_curv,result_is_curv
      logical,optional,intent(OUT)::out_iscnv,out_iscurv

      type(grid_mg_def),pointer :: gbc_def

c     Local variables

      integer :: i,j,k,stencil_width,ieq,neql
     .          ,imng,imxg,jmng,jmxg,kmng,kmxg
     .          ,imnl,imxl,jmnl,jmxl,kmnl,kmxl

      integer,save :: order=-1 !Needed to "remember" order of extrapolation at
                               !boundaries from previous calls when not explicitly
                               !set (exits with error when unset)

      integer,save :: ivar(2)  !Needed to "remember" ivar from previous calls
                               !when icomp is not provided.

      logical,save :: iscnv,iscurv !Needed to "remember' these from previous calls
                                   !when they are not provided.

      logical      :: res_is_curv,fpointers,vector

      real(8)      :: car(3),curv(3)

      real(8),pointer,dimension(:,:,:,:) :: v_cov,v_cnv,vv0

      character(3) :: cvec

      type(mg_ctx),pointer :: mgctx => null()

c     Begin program

      neql = neq !To avoid problems if neq is hardwired

      if (PRESENT(is_vec)) then
        vector = is_vec
      else
        vector = (neql == 3 .or. neql == 6) !Backward compatibility with previous
                                           ! versions of this routine
      endif

      if (.not.vector) then
        call pstop('setMGBC_stg','Input is not a vector')
      endif

c     Optional arguments

      if (PRESENT(icomp))  ivar(1:size(icomp))  = icomp

      if (PRESENT(iorder)) order = iorder

c     Consistency check

      if (     gbc_def%nxv(iig) /= nnx
     .    .or. gbc_def%nyv(iig) /= nny
     .    .or. gbc_def%nzv(iig) /= nnz) then
        write (*,*) gbc_def%nxv(iig), nnx
        write (*,*) gbc_def%nyv(iig), nny
        write (*,*) gbc_def%nzv(iig), nnz
        call pstop('setMGBC_stg','Grid sizes do not agree')
      endif

      if (order < 0) then
        call pstop('setMGBC_stg','Order of extrapolation not set')
      endif

c     Allocate MG pointers (if not allocated earlier)

      call allocMGPointers(neql,gbc_def,mgctx,fpointers)

c     Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(mgctx,gpos,iig,imnl,imxl,jmnl,jmxl,kmnl,kmxl)

c     Check LOCAL limits (return if not close to local domain boundaries)

      if (     (imnl > 1 .and. imxl < nnx)
     .    .and.(jmnl > 1 .and. jmxl < nny)
     .    .and.(kmnl > 1 .and. kmxl < nnz)) return

c     Find GLOBAL limits (required for BC routine setBC)

      imng = imnl + gbc_def%ilo(iig)-1
      imxg = imxl + gbc_def%ilo(iig)-1
      jmng = jmnl + gbc_def%jlo(iig)-1
      jmxg = jmxl + gbc_def%jlo(iig)-1
      kmng = kmnl + gbc_def%klo(iig)-1
      kmxg = kmxl + gbc_def%klo(iig)-1

c     Select operation

      if (PRESENT(is_cnv )) iscnv = is_cnv

      if (PRESENT(is_curv)) then
        iscurv = is_curv
      else
        iscurv = .true.
      endif

      if (PRESENT(result_is_curv)) then
        res_is_curv = result_is_curv
      else
        res_is_curv = iscurv
      endif

      if (.not.iscurv) iscnv=.true.

      cvec = 'cov'
      if (iscnv) cvec = 'cnv'

      if (PRESENT(arr0)) then
        vv0 => arr0
      else
        allocate(vv0(0:nnx+1,0:nny+1,0:nnz+1,neql))
        vv0 = 0d0
      endif

      if (.not.iscurv) then  !Transform to curvilinear

        do ieq=1,neql/3
          call  XformVector_stg_mesh2(gbc_def,iig
     .              ,array(:,:,:,1+3*(ieq-1):3*ieq),'car',cvec)
        enddo

      endif

      if (PRESENT(arr_cov)) then
        v_cnv => array
        v_cov => arr_cov
      else
        if (iscnv) then
          v_cnv => array
          allocate(v_cov(0:nnx+1,0:nny+1,0:nnz+1,neql))
          v_cov = 0d0
        else
          v_cov => array
          allocate(v_cnv(0:nnx+1,0:nny+1,0:nnz+1,neql))
          v_cnv = 0d0
        endif
      endif

      if (iscnv) then
        do ieq=1,neql/3
          v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .          = XformToCov_stg_mesh(gbc_def,iig
     .                               ,v_cnv(:,:,:,1+3*(ieq-1):3*ieq))
        enddo
      else
        do ieq=1,neql/3
          v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .         = XformToCnv_stg_mesh(gbc_def,iig
     .                              ,v_cov(:,:,:,1+3*(ieq-1):3*ieq))
        enddo
      endif

      do ieq=1,neql/3
        call imposeBConVector_stg(gbc_def,ivar(ieq),3,nnx,nny,nnz
     .          ,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .          ,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .          ,vv0  (:,:,:,1+3*(ieq-1):3*ieq)
     .          ,bcnd (:    ,1+3*(ieq-1):3*ieq)
     .          ,iig,iig,iig
     .          ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .          ,is_cnv=iscnv,iorder=order)
      enddo

      if (.not.PRESENT(arr0)) deallocate(vv0)

      if (.not.PRESENT(arr_cov)) then
        if (iscnv) then
          deallocate(v_cov)
        else
          deallocate(v_cnv)
        endif
      endif

      if (.not.res_is_curv) then

        do ieq=1,neql/3
          call  XformVector_stg_mesh2(gbc_def,iig
     .              ,array(:,:,:,1+3*(ieq-1):3*ieq),cvec,'car')
        enddo

      endif

      !For output
      if (PRESENT(out_iscnv )) out_iscnv =iscnv
      if (PRESENT(out_iscurv)) out_iscurv=iscurv

c     Deallocate MG pointers (if allocated in this routine)

      call deallocMGPointers(fpointers,mgctx)

c     End

      end subroutine setMGBC_stg

      end module setMGBC_interface

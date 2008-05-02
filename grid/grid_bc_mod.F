c module constants
c ######################################################################
      module constants

        real(8) :: pi
        real(8), allocatable, dimension(:,:,:)   :: zeros,ones
        real(8), allocatable, dimension(:,:,:,:) :: vzeros

      end module constants

c module BCS_variables
c####################################################################
      module BCS_variables

        use oned_int

        use grid_mpi

        use grid_operations

        use math

        use constants

cc        integer    :: imax,imin,jmax,jmin,kmax,kmin

        integer    :: nnvar,nxbc,nybc,nzbc,igxbc,igybc,igzbc

        real(8),allocatable,dimension(:,:) :: rhs

        logical :: symm=.false.,bc_debug=.false.

      end module BCS_variables

c module singularBCinterface
c####################################################################
      module singularBCinterface

        use BCS_variables

      contains

c     scalarSingularBC2
c     #################################################################
      subroutine scalarSingularBC2(array,order)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * array: contains variable on which singular BC is imposed
c        * order: order of interpolation towards singular point
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: order
      real(8)    :: array(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,order1,ilog,ihig
      real(8)    :: lavg_q,lavg_vol,avg_q,avg_vol,rho0(1),vol,x0(1)

c     Begin program

      ilog = grid_params%ilo(1)
      ihig = grid_params%ihi(1)

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Set interpolation order

      order1 = min(order,ihig-ilog)

c     Find scalar value at singular point by averaging neighboring cells

      do k=1,nzbc

c       Perform local sum
        lavg_q   = 0d0
        lavg_vol = 0d0

        do j=1,nybc
          call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
          vol = gmetric%grid(igxbc)%dvol(i,j,k)
          x0  = grid_params%xx(ig-1)

          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,array(i:i+order1+1,j,k)
     .                    ,1,x0,rho0,order1)

          lavg_q   = lavg_q   + vol*rho0(1)
          lavg_vol = lavg_vol + vol
        enddo

c       Perform global sum (and broadcast it to all processors)

#if defined(petsc)
        call MPI_Allreduce(lavg_q  ,avg_q  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_vol,avg_vol,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
#else
        avg_q   = lavg_q
        avg_vol = lavg_vol
#endif

c       Perform average

        array(0,:,k) = avg_q/avg_vol

      enddo

c     End program

      end subroutine scalarSingularBC2

c     scalarSingularBC
c     #################################################################
      subroutine scalarSingularBC(array)
c     -----------------------------------------------------------------
c     Imposes singular point BC. On input:
c        * array: contains variable on which singular BC is imposed
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: array(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,jj,ilog,jlog,nyg,jglobal

#if defined(petsc)
      real(8),allocatable,dimension(:)   :: local_x,global_x
      real(8),allocatable,dimension(:,:) :: garray
#endif

c     Begin program

      ilog = grid_params%ilo(igxbc)
      jlog = grid_params%jlo(igybc)

      if (npy < 1) call pstop('scalarSingularBC','npy < 1')

      nyg = nybc*npy

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Parallel Allgather in SP communicator

#if defined(petsc)
cc      local_x = reshape(array(i,1:nybc,1:nzbc), (/ nybc*nzbc /))
cc      write (*,*) 'DIAG -- scalarSPBC; proc=',my_rank
cc     .           ,array(i,1:nybc,1:nzbc)
cc      call MPI_Allgather(local_x ,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,global_x,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_COMM_SP,mpierr)
cc      garray = reshape(global_x, (/ nyg,nzbc /))
cc      write (*,*) 'DIAG -- scalarSPBC; proc=',my_rank
cc     .           ,garray
cc      stop

      allocate(global_x(nyg*nzbc),garray(nyg,nzbc))
      call MPI_Allgather(array(i,1:nybc,1:nzbc)
     .                           ,nybc*nzbc,MPI_DOUBLE_PRECISION
     .                  ,global_x,nybc*nzbc,MPI_DOUBLE_PRECISION
     .                  ,MPI_COMM_SP,mpierr)
      do k=1,npy
        garray(1+nybc*(k-1):nybc*k,:)
     .     = reshape(global_x(1+nybc*nzbc*(k-1):nybc*nzbc*k)
     .              , (/ nybc,nzbc /))
      enddo
      deallocate(global_x)

cc      call MPI_Allgather(array(i,1:nybc,1:nzbc)
cc     .                         ,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,garray,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_COMM_SP,mpierr)
#endif

c     Element exchange in poloidal direction

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('scalarSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

      do k=1,nzbc
        do j=1,nybc
          jglobal = j + jlog - 1
          jj = mod(jglobal+nyg/2,nyg)
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

c     vectorSingularBC2
c     #################################################################
      subroutine vectorSingularBC2(vec,cov,order)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: order
      real(8)    :: vec(0:nxbc+1,0:nybc+1,0:nzbc+1,3)
      logical    :: cov

c     Local variables

      integer    :: i,j,k

      integer    :: ic,ig,jg,kg,order1,ilog,ihig
      real(8)    :: vol,lavg_x,lavg_y,lavg_z,lavg_vol
     .                 , avg_x, avg_y, avg_z, avg_vol
      real(8)    :: daxdr(1),daydr(1),dummy1,dummy2
      real(8),allocatable,dimension(:) :: ax0,ay0,az0,cx,cy,cz,xx

c     External

c     Begin program

      ilog = grid_params%ilo(1)
      ihig = grid_params%ihi(1)

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Set interpolation order

      order1 = min(order,ihig-ilog)

c     Allocate quantities for interpolation

      allocate(ax0(nzbc),ay0(nzbc),az0(nzbc))
      allocate(cx(0:order1+1),cy(0:order1+1),cz(0:order1+1)
     .        ,xx(0:order1+1))

      call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

      xx(0:order1+1) = grid_params%xx(ig-1:ig+order1)

c     Find average cartesian coordinates

      do k=1,nzbc

c       Perform local sum
        lavg_x   = 0d0
        lavg_y   = 0d0
        lavg_z   = 0d0
        lavg_vol = 0d0

        do j=1,nybc
          call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)
          vol = gmetric%grid(igxbc)%dvol(i,j,k)

          do ic=1,order1+1
            call transformVectorToCartesian(ic,j,k,igxbc,igybc,igzbc
     .           ,vec(ic,j,k,1),vec(ic,j,k,2),vec(ic,j,k,3),cov
     .           ,cx(ic),cy(ic),cz(ic))
          enddo

          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,cx(1:order1+1),1,grid_params%xx(ig-1:ig-1)
     .                    ,cx(0),order1)
          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,cy(1:order1+1),1,grid_params%xx(ig-1:ig-1)
     .                    ,cy(0),order1)
          call IntDriver1d(order1+1,grid_params%xx(ig:ig+order1)
     .                    ,cz(1:order1+1),1,grid_params%xx(ig-1:ig-1)
     .                    ,cz(0),order1)

          lavg_x   = lavg_x   + vol*cx(0)
          lavg_y   = lavg_y   + vol*cy(0)
          lavg_z   = lavg_z   + vol*cz(0)
          lavg_vol = lavg_vol + vol

        enddo

c       Perform global sum (and broadcast it to all processors)

#if defined(petsc)
        call MPI_Allreduce(lavg_x  ,avg_x  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_y  ,avg_y  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_z  ,avg_z  ,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
        call MPI_Allreduce(lavg_vol,avg_vol,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_SP,mpierr)
#else
        avg_x   = lavg_x
        avg_y   = lavg_y
        avg_z   = lavg_z
        avg_vol = lavg_vol
#endif

c       Perform average

        ax0(k) = avg_x/avg_vol
        ay0(k) = avg_y/avg_vol
        az0(k) = avg_z/avg_vol

      enddo

      i = ilog-1
      do k=1,nzbc
        do j=0,nybc+1           !This also fills ghost cells
          call transformVectorToCurvilinear(i,j,k,igxbc,igybc,igzbc
     .                 ,ax0(k),ay0(k),az0(k),cov
     .                 ,vec(i,j,k,1),vec(i,j,k,2),vec(i,j,k,3))
        enddo
      enddo

      deallocate(ax0,ay0,az0)
      deallocate(cx,cy,cz,xx)

c     End program

      end subroutine vectorSingularBC2

c     vectorSingularBC
c     #################################################################
      subroutine vectorSingularBC(ivar,vec,cov)
c     -----------------------------------------------------------------
c     Averages vector components around singular point and calculates
c     curvilinear components at singular point.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ivar
      real(8)    :: vec(0:nxbc+1,0:nybc+1,0:nzbc+1)
      logical    :: cov

c     Local variables

      integer    :: i,j,k,jj,ilog,jlog,nyg,jglobal

#if defined(petsc)
      real(8),allocatable,dimension(:)   :: local_x,global_x
      real(8),allocatable,dimension(:,:) :: garray
#endif

c     External

c     Begin program

      ilog = grid_params%ilo(igxbc)
      jlog = grid_params%jlo(igybc)

      if (npy < 1) call pstop('vectorSingularBC','npy < 1')

      nyg = grid_params%nygl(igybc)

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Parallel Allgather in SP communicator

#if defined(petsc)
cc      allocate(local_x(nybc*nzbc),global_x(nyg*nzbc),garray(nyg,nzbc))
cc      local_x = reshape(vec(i,1:nybc,1:nzbc), (/ nybc*nzbc /))
      allocate(global_x(nyg*nzbc),garray(nyg,nzbc))
      call MPI_Allgather(vec(i,1:nybc,1:nzbc)
     .                           ,nybc*nzbc,MPI_DOUBLE_PRECISION
     .                  ,global_x,nybc*nzbc,MPI_DOUBLE_PRECISION
     .                  ,MPI_COMM_SP,mpierr)
      do k=1,npy
        garray(1+nybc*(k-1):nybc*k,:)
     .     = reshape(global_x(1+nybc*nzbc*(k-1):nybc*nzbc*k)
     .              , (/ nybc,nzbc /))
      enddo
      deallocate(global_x)
cc      deallocate(local_x,global_x)

cc      allocate(local_x(nybc),global_x(nyg),garray(nyg,nzbc))
cc      do k=1,nzbc
cc        local_x = vec(i,1:nybc,k)
cc        call MPI_Allgather(local_x ,nybc,MPI_DOUBLE_PRECISION
cc     .                    ,global_x,nybc,MPI_DOUBLE_PRECISION
cc     .                    ,MPI_COMM_SP,mpierr)
cc        garray(:,k) = global_x
cc      enddo
cc      deallocate(local_x,global_x)

cc      allocate(garray(nyg,nzbc))
cc      call MPI_Allgather(vec(i,1:nybc,1:nzbc)
cc     .                         ,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,garray,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_COMM_SP,mpierr)
#endif

c     Element exchange in poloidal direction

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('vectorSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

      if (coords == 'hl2') then
        do k=1,nzbc
          do j=1,nybc
            jglobal = j + jlog - 1
            jj = mod(jglobal+nyg/2,nyg)
            if (jj == 0) jj = nyg

#if defined(petsc)
            vec(i-1,j,k) =  garray(jj,k)
#else
            vec(i-1,j,k) =  vec(i,jj,k)
#endif
          enddo
        enddo
      else
        do k=1,nzbc
          do j=1,nybc
            jglobal = j + jlog - 1
            jj = mod(jglobal+nyg/2,nyg)
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

      end module singularBCinterface

c module imposeBCinterface
c #####################################################################
      module imposeBCinterface

        use singularBCinterface

        INTERFACE setBC
          module procedure imposeBConScalar,imposeBConVector
        end INTERFACE

cc        INTERFACE
cc          subroutine scalarDirichletBC(array,array0,ieq,dim,loc,order)
cc            use BCS_variables
cc            integer    :: ieq,dim,loc,order
cc            real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
cc     .                   ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)
cc          end subroutine scalarDirichletBC
cc        end INTERFACE
cc
cc        INTERFACE
cc          subroutine vectorDirichletBC(ivar,array,array0,ieq,dim,loc
cc     .                                ,order)
cc            use BCS_variables
cc            integer    :: ieq,dim,loc,order,ivar
cc            real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)
cc     .                   ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,*)
cc          end subroutine vectorDirichletBC
cc        end INTERFACE
cc
cc        INTERFACE
cc          subroutine scalarNeumannBC(array,ieq,dim,loc)
cc            use BCS_variables
cc            integer    :: ieq,dim,loc
cc            real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
cc          end subroutine scalarNeumannBC
cc        end INTERFACE
cc
cc        INTERFACE
cc          subroutine vectorNeumannBC(ivar,array,ieq,dim,loc)
cc            use BCS_variables
cc            integer    :: ieq,dim,loc,ivar
cc            real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,*)
cc          end subroutine vectorNeumannBC
cc        end INTERFACE

        integer    :: iimax,iimin,jjmax,jjmin,kkmax,kkmin
        integer    :: iimxp,iimnm,jjmxp,jjmnm,kkmxp,kkmnm

        integer   ,private :: bctype,dim,loc,ibc,i,j,k

        !Define order of interpolation for EDGE and CORNER treatment
        !(linear to preserve nearest-neighbors stencil; this is 
        ! important for colored diagonal formation)
        integer   ,private :: order=1

        real(8),allocatable,dimension(:,:,:,:) :: v_cov,v_cnv,vv0

      contains

c     findBCLoopLimits
c     ##############################################################
      subroutine findBCLoopLimits(dim,loc,i1 ,i2 ,j1 ,j2 ,k1 ,k2
     .                                   ,imn,imx,jmn,jmx,kmn,kmx)
c     --------------------------------------------------------------
c     Finds local limits for BC loops. On input, it takes global 
c     limits (i1,i2), (j1,j2), (k1,k2) and outputs local loop limits
c     (imn,imx), (jmn,jmx), (kmn,kmx).
c     --------------------------------------------------------------

c     Call variables

      integer   ,intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2
      integer   ,intent(OUT):: imn,imx,jmn,jmx,kmn,kmx

c     Begin program

      select case(dim)
      case (1)

        call selectLimits(i1,i2,grid_params%ilo(igxbc)
     .                   ,nxbc,grid_params%nxgl(igxbc),imn,imx)

        jmn = j1 - grid_params%jlo(igybc) + 1
        jmx = j2 - grid_params%jlo(igybc) + 1

        kmn = k1 - grid_params%klo(igzbc) + 1
        kmx = k2 - grid_params%klo(igzbc) + 1

      case(2)
        imn = i1 - grid_params%ilo(igxbc) + 1
        imx = i2 - grid_params%ilo(igxbc) + 1

        call selectLimits(j1,j2,grid_params%jlo(igybc)
     .                   ,nybc,grid_params%nygl(igybc),jmn,jmx)

        kmn = k1 - grid_params%klo(igzbc) + 1
        kmx = k2 - grid_params%klo(igzbc) + 1

      case(3)
        imn = i1 - grid_params%ilo(igxbc) + 1
        imx = i2 - grid_params%ilo(igxbc) + 1

        jmn = j1 - grid_params%jlo(igybc) + 1
        jmx = j2 - grid_params%jlo(igybc) + 1

        call selectLimits(k1,k2,grid_params%klo(igzbc)
     .                   ,nzbc,grid_params%nzgl(igzbc),kmn,kmx)

      end select

c     End program

      contains

c     selectLimits
c     #####################################################################
      subroutine selectLimits(igl,igr,igmin,nxl,nxg,ill,ilr)

c     ---------------------------------------------------------------------
c     Selects LOCAL limits for BC treatment at boundary location determined
c     by (dim,loc).
c     ---------------------------------------------------------------------

        implicit none

        integer    :: igl,igr,igmin,nxl,nxg,ill,ilr

c     Begin program

cc        do_nothing = .false.

        select case(loc)
        case(0)
cc          if (igl > 1) do_nothing = .true.
          ill = min(igl - igmin + 1,1)
          ilr = ill
        case(1)
cc          if (igr < nxg) do_nothing = .true.
          ill = max(igr - igmin + 1,nxl)
          ilr = ill
        end select

      end subroutine selectLimits

      end subroutine findBCLoopLimits

c     do__nothing
c     #####################################################################
      function do__nothing(dim,loc,i1,i2,j1,j2,k1,k2)

c     ---------------------------------------------------------------------
c     Selects whether to fill or not ghost cells of a given local domain
c     based on its position in the global grid.
c     ---------------------------------------------------------------------

        implicit none

        logical    :: do__nothing
        integer   ,intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2

        integer    :: ibc

c     Begin program

        do__nothing = .false.
        if (asm) return      !Return if additive Schwartz method is used (PC)

        select case(dim)
        case(1)
          call checkLimits(i1,i2,grid_params%nxgl(igxbc))
        case(2)
          call checkLimits(j1,j2,grid_params%nygl(igybc))
        case(3)
          call checkLimits(k1,k2,grid_params%nzgl(igzbc))
        end select

        !Do nothing for periodic BCs in parallel (PETSc does it)
cc        ibc = (1+loc)+2*(dim-1)
cccc        if ((.not.do__nothing).and.(bcond(ibc) == PER.and.asm_dir(dim)))
cc        if((.not.do__nothing).and.(bcond(ibc)==PER).and.(np>1))
cc     .       do__nothing =.true. 

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

        integer    :: igl,igr,ng

c     Begin program

        select case(loc)
        case(0)
          if (igl > 1 ) do__nothing = .true.
        case(1)
          if (igr < ng) do__nothing = .true.
        end select

      end subroutine checkLimits

      end function do__nothing

c     selectFace
c     #####################################################################
      function selectFace(ibc)

c     ---------------------------------------------------------------------
c     Determines whether boundary location determined by ibc(=1-6) is
c     selected for boundary treatment.
c     ---------------------------------------------------------------------

      implicit none

      logical    :: selectFace

      integer    :: dim,loc,ibc

c     Begin program

      loc = abs(mod(ibc,2)-2)-1
      dim = 1 + (ibc-1-loc)/2

      selectFace= 
     .     .not.do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax)
cc     .     .not.do__nothing(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp)

cc      selectFace = .true.

cc      if (.not.selectFace) write (*,*) 'DIAG -- selectFace'
cc     .                          ,my_rank,dim,loc,ibc, selectFace

      end function selectFace

c     setup_petsc_BC
c     ######################################################################
      subroutine setup_petsc_BC

c     -----------------------------------------------------------------
c     Sets flags to indicate that parallel BCs are NOT to be imposed
c     on the fortran side (are fully provided by PETSc).
c     -----------------------------------------------------------------

      implicit none

      par_bc = .false.  !No parallel BC (provided by PETSc)
      asm    = .false.  !No ASM (no extrapolation to corner, edges)

      end subroutine setup_petsc_BC

c     setASMflag
c     ######################################################################
      subroutine setASMflag(asm_flag)

c     -----------------------------------------------------------------
c     Sets flag to indicate Additive Schwartz method.
c     -----------------------------------------------------------------

      implicit none

      logical :: asm_flag

      asm = asm_flag

      par_bc = .not.asm

      end subroutine setASMflag

c     setASM_PC_BC
c     ######################################################################
      subroutine setASM_PC_BC(igx,igy,igz,bcs,IIBC)

c     -----------------------------------------------------------------
c     Sets boundary conditions for Additive Schwartz method PC.
c     -----------------------------------------------------------------

      implicit none

      integer    :: igx,igy,igz,bcs(:,:)

      integer   ,optional :: IIBC

      integer    :: dim,loc,ibc,igmin,igmax,jgmin,jgmax,kgmin,kgmax

      logical    :: dn

      if (np == 1) then
        call setASMflag(.false.)
        return
      endif

c     Reset BC grid level

      igxbc = igx
      igybc = igy
      igzbc = igx

c     Find grid global limits

      igmin = grid_params%ilo(igx)
      igmax = grid_params%ihi(igx)
      jgmin = grid_params%jlo(igy)
      jgmax = grid_params%jhi(igy)
      kgmin = grid_params%klo(igz)
      kgmax = grid_params%khi(igz)

c     Determine ASM BCs (requires ASM flag set to false)

      call setASMflag(.false.)

      if (PRESENT(IIBC)) then
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            dn=do__nothing(dim,loc,igmin,igmax,jgmin,jgmax,kgmin,kgmax)
     .         .or.((bcond(ibc)==PER).and.asm_dir(dim))  !Deal with interrupted periodic boundaries

cc            if (dn) bcs(ibc,:) = IIBC
            if (dn) bcs(ibc,:) = IIBC*bcs(ibc,:)/abs(bcs(ibc,:))
          enddo
        enddo
      else
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (bcond(ibc) == PER.and.asm_dir(dim)) then !Take care of PER boundaries in ASM
              bcs(ibc,:) = EXT
            endif
          enddo
        enddo
      endif

c     Set ASM flag to true to perform ASM in PC

      call setASMflag(.true.)

      end subroutine setASM_PC_BC

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(ieq,nx,ny,nz,array,array0,bcond
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ieq,bcond(6),nx,ny,nz,iigx,iigy,iigz
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

c     Local variables

      logical :: petsc__bc

c     Begin program

c     Initialize interpolation order

      if (PRESENT(iorder)) then
        order = iorder
      else
        order = 2
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
        iimin = grid_params%ilo(iigx)
        iimax = grid_params%ihi(iigx)
        jjmin = grid_params%jlo(iigy)
        jjmax = grid_params%jhi(iigy)
        kkmin = grid_params%klo(iigz)
        kkmax = grid_params%khi(iigz)
      endif

      iimnm =  max(iimin - 1,1)  
      iimxp =  min(iimax + 1,grid_params%nxgl(iigx))
      jjmnm =  max(jjmin - 1,1)  
      jjmxp =  min(jjmax + 1,grid_params%nygl(iigy))
      kkmnm =  max(kkmin - 1,1)  
      kkmxp =  min(kkmax + 1,grid_params%nzgl(iigz))

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc =      (np > 1)
     .           .and. (    ((.not.asm) .and. par_bc)    !For PC (not ASM)
     .                  .or.(ieq > NVAR             ))   !For auxiliary variables

      !Fill ghost cells
      if (petsc__bc) call fillPetscGhostCells(array,igxbc)
#endif

c     Impose BCs

      do bctype=1,BCLIM             !Enforces a particular order in the BCs
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (abs(bcond(ibc)) == bctype) then
              if (bcond(ibc) < 0 .and. bctype == SYM) then
                symm = .true.
                call FillGhostNodes(ieq,1,1,dim,loc,DIR,array,array0)
                symm = .false.
              elseif (bcond(ibc) > 0 .and. bctype == SYM) then
                symm = .true.
                call FillGhostNodes(ieq,1,1,dim,loc,NEU,array,array0)
                symm = .false.
              else
                call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
              endif
            endif
          enddo
        enddo
      enddo

c     Fill edges

      do dim=1,3
        do loc=0,1
          call fillEdges(array,array0,dim,loc,bcond)
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          call fillCorners(array,array0,dim,loc,bcond)
        enddo
      enddo

c     Periodic boundary synchronizing

cc      bctype=PER
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc          if (bcond(ibc) == bctype) then
cc            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
cc          endif
cc        enddo
cc      enddo

c     End program

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(fcomp,ncomp,nx,ny,nz,v_cnv,v_cov,var0
     .                          ,bcond,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                          ,iorder,is_cnv,is_vec)
c     -----------------------------------------------------------------
c     Imposes BC on vector field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ncomp,bcond(6,ncomp),fcomp,nx,ny,nz,iigx,iigy,iigz
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,ncomp)
     .             ,var0 (0:nx+1,0:ny+1,0:nz+1,ncomp)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,ncomp)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder
      logical   ,optional,intent(IN) :: is_cnv,is_vec

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: ivar,ieq,ibc,loc,dim,bctype
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: cov_to_cnv,to_cnv,is__vec,petsc__bc

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

c     Initialize interpolation order

      if (PRESENT(iorder)) then
        order = iorder
      else
        order = 2
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
        iimin = grid_params%ilo(iigx)
        iimax = grid_params%ihi(iigx)
        jjmin = grid_params%jlo(iigy)
        jjmax = grid_params%jhi(iigy)
        kkmin = grid_params%klo(iigz)
        kkmax = grid_params%khi(iigz)
      endif

      iimnm =  max(iimin - 1,1)  
      iimxp =  min(iimax + 1,grid_params%nxgl(iigx))
      jjmnm =  max(jjmin - 1,1)  
      jjmxp =  min(jjmax + 1,grid_params%nygl(iigy))
      kkmnm =  max(kkmin - 1,1)  
      kkmxp =  min(kkmax + 1,grid_params%nzgl(iigz))

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc =      (np > 1)
     .           .and. (    ((.not.asm) .and. par_bc)    !For PC (not ASM)
     .                  .or.(fcomp > NVAR           ))   !For auxiliary variables

      !Fill ghost cells
      if (petsc__bc) then
        if (is__vec .and. to_cnv) then
          do ieq=1,ncomp
            call fillPetscGhostCells(v_cov(:,:,:,ieq),igxbc)
          enddo
        else
          do ieq=1,ncomp
            call fillPetscGhostCells(v_cnv(:,:,:,ieq),igxbc)
          enddo
        endif
      endif
#endif

c     Find complete curvilinear representation

      if (is__vec) then

#if defined(petsc)
       call fromGlobalToLocalLimits(iimin,jjmin,kkmin,ilmin,jlmin,klmin
     $                             ,igxbc,igybc,igzbc)
       call fromGlobalToLocalLimits(iimax,jjmax,kkmax,ilmax,jlmax,klmax
     $                             ,igxbc,igybc,igzbc)

       do k=klmin-1,klmax+1
         do j=jlmin-1,jlmax+1
           do i=ilmin-1,ilmax+1
#else
       do k=kkmin-1,kkmax+1
         do j=jjmin-1,jjmax+1
           do i=iimin-1,iimax+1
#endif
              call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),to_cnv)
           enddo
         enddo
       enddo

      endif

c     Impose BCs

      cov_to_cnv = .false.

      do bctype=1,BCLIM            !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              ieq = ivar + fcomp - 1

              if (abs(bcond(ibc,ivar)) == bctype) then

                if (bcond(ibc,ivar) < 0 .and. is__vec) then
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cov,var0,cnv=.false.)
                  cov_to_cnv = .true.
                else
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cnv,var0)
                endif

             endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components => contravariant

      if (cov_to_cnv) call synchronize_cnv(ncomp,v_cnv,v_cov,bcond)

c     Fill edges

      do dim=1,3
        do loc=0,1
          do ivar=1,ncomp
            call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                    ,dim,loc,bcond(:,ivar))
          enddo
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          do ivar=1,ncomp
            call fillCorners(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                      ,dim,loc,bcond(:,ivar))
          enddo
        enddo
      enddo

c     Periodic boundary synchronizing (need to change periodicBC)

cc      write (*,*) 'DIAG -- imposeBConVector'
cc      bctype=PER
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc
cc          do ivar = 1,ncomp
cc            ieq = ivar + fcomp - 1
cc
cc            if (abs(bcond(ibc,ivar)) == bctype) then
cc              call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype,v_cnv
cc     .                               ,var0)
cc            endif
cc          enddo
cc
cc        enddo
cc      enddo

c     Find covariant components at ALL boundary ghost cells

      if (is__vec) then

        do dim = 1,3
          do loc = 0,1

          call findBCLoopLimits(dim,loc,iimin-1,iimax+1
     .                                 ,jjmin-1,jjmax+1
     .                                 ,kkmin-1,kkmax+1
     .                                 ,imin,imax,jmin,jmax,kmin,kmax)

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax
                call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
              enddo
            enddo
          enddo

          enddo
        enddo

c       Impose BCs on covariant components (if not Cartesian)

cc        if (coords /= 'car') then
cc
ccc         Impose SP BCs
cc
cc          bctype=SP
cc          do dim=1,3
cc            do loc=0,1
cc              ibc = (1+loc)+2*(dim-1)
cc
cc              do ivar = 1,ncomp
cc                ieq = ivar + fcomp - 1
cc
cc                if (abs(bcond(ibc,ivar)) == bctype) then
cc                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
cc     .                               ,v_cov,var0,cnv=.false.)
cc                endif
cc              enddo
cc
cc            enddo
cc          enddo
cc
ccc         Fill edges
cc
cc          do dim=1,3
cc            do loc=0,1
cc              do ivar=1,ncomp
cc                call fillEdges(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
cc     .                        ,dim,loc,bcond(:,ivar))
cc              enddo
cc            enddo
cc          enddo
cc
ccc         Fill corners
cc
cc          do dim=1,3
cc            do loc=0,1
cc              do ivar=1,ncomp
cc                call fillCorners(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
cc     .                          ,dim,loc,bcond(:,ivar))
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif

      endif

c     End program

      end subroutine imposeBConVector

c     synchronize_cnv
c     #################################################################
      subroutine synchronize_cnv(ncomp,v_cnv,v_cov,bcond)
c     -----------------------------------------------------------------
c     Finds all contravariant components at boundaries with BC in terms
c     of covariant components.
c
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravariant
c     components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ncomp,bcond(6,ncomp)
      real(8)    :: v_cnv(0:nxbc+1,0:nybc+1,0:nzbc+1,ncomp)
     .             ,v_cov(0:nxbc+1,0:nybc+1,0:nzbc+1,ncomp)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,dim,loc,ig,jg,kg,ivar,ibc
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)

c     Begin program

      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp

            if (ivar == dim) then  !Select tangential components
              cycle
            elseif (bcond(ibc,ivar) < 0) then

              if (.not.selectFace(ibc)) cycle

              call findBCLoopLimits(dim,loc,iimin,iimax
     .                                     ,jjmin,jjmax
     .                                     ,kkmin,kkmax
     .                                     ,imin,imax,jmin
     .                                     ,jmax,kmin,kmax)

              select case (ibc)
              case (1)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i-1,j,k,:,:)

                      v_cov(i-1,j,k,1) = -(gsuper(1,2)*v_cov(i-1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i-1,j,k,3)
     .                                    -v_cnv(i-1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i-1,j,k
     .               ,igxbc,igybc,igzbc
     .               ,v_cov(i-1,j,k,1),v_cov(i-1,j,k,2),v_cov(i-1,j,k,3)
     .               ,v_cnv(i-1,j,k,1),v_cnv(i-1,j,k,2),v_cnv(i-1,j,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (2)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i+1,j,k,:,:)

                      v_cov(i+1,j,k,1) = -(gsuper(1,2)*v_cov(i+1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i+1,j,k,3)
     .                                    -v_cnv(i+1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i+1,j,k
     .               ,igxbc,igybc,igzbc
     .               ,v_cov(i+1,j,k,1),v_cov(i+1,j,k,2),v_cov(i+1,j,k,3)
     .               ,v_cnv(i+1,j,k,1),v_cnv(i+1,j,k,2),v_cnv(i+1,j,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (3)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j-1,k,:,:)

                      v_cov(i,j-1,k,2) = -(gsuper(2,1)*v_cov(i,j-1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j-1,k,3)
     .                                    -v_cnv(i,j-1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j-1,k
     .               ,igxbc,igybc,igzbc
     .               ,v_cov(i,j-1,k,1),v_cov(i,j-1,k,2),v_cov(i,j-1,k,3)
     .               ,v_cnv(i,j-1,k,1),v_cnv(i,j-1,k,2),v_cnv(i,j-1,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (4)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j+1,k,:,:)

                      v_cov(i,j+1,k,2) = -(gsuper(2,1)*v_cov(i,j+1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j+1,k,3)
     .                                    -v_cnv(i,j+1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j+1,k
     .               ,igxbc,igybc,igzbc
     .               ,v_cov(i,j+1,k,1),v_cov(i,j+1,k,2),v_cov(i,j+1,k,3)
     .               ,v_cnv(i,j+1,k,1),v_cnv(i,j+1,k,2),v_cnv(i,j+1,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (5)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j,k-1,:,:)

                      v_cov(i,j,k-1,3) = -(gsuper(3,1)*v_cov(i,j,k-1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k-1,2)
     .                                    -v_cnv(i,j,k-1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k-1
     .               ,igxbc,igybc,igzbc
     .               ,v_cov(i,j,k-1,1),v_cov(i,j,k-1,2),v_cov(i,j,k-1,3)
     .               ,v_cnv(i,j,k-1,1),v_cnv(i,j,k-1,2),v_cnv(i,j,k-1,3)
     .               ,.true.)
                    enddo
                  enddo
                enddo

              case (6)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j,k+1,:,:)

                      v_cov(i,j,k+1,3) = -(gsuper(3,1)*v_cov(i,j,k+1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k+1,2)
     .                                    -v_cnv(i,j,k+1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k+1
     .               ,igxbc,igybc,igzbc
     .               ,v_cov(i,j,k+1,1),v_cov(i,j,k+1,2),v_cov(i,j,k+1,3)
     .               ,v_cnv(i,j,k+1,1),v_cnv(i,j,k+1,2),v_cnv(i,j,k+1,3)
     .               ,.true.)
                    enddo
                  enddo
                enddo

              end select

            endif
          enddo

        enddo
      enddo

c     End program

      end subroutine synchronize_cnv

c     fillGhostNodes
c     ###############################################################
      subroutine fillGhostNodes(ieq,ivar,nvar,dim,loc,bctype
     .                         ,array,array0,cnv)

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

      integer    :: ieq,dim,loc,bctype,nvar,ivar
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,nvar)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar)
      logical,optional :: cnv

c     Local variables

cc      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: neq,ibc,i,j,k,ig,jg,kg
      logical    :: cov

c     Begin program

      nnvar = nvar   !This transfers info to vectorDirichlet (needed? Check)

      if (PRESENT(cnv)) then
        cov = .not.cnv
      else
        cov = .false.
      endif

      ibc = (1+loc)+2*(dim-1)
      
c     Determine boundary limits

      if (.not.selectFace(ibc)) return

cccc      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
cc      call findBCLoopLimits(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
cc     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

c     Allocate rhs

      select case(dim)
      case(1)
        allocate(rhs(0:nybc+1,0:nzbc+1))
      case(2)
        allocate(rhs(0:nxbc+1,0:nzbc+1))
      case(3)
        allocate(rhs(0:nxbc+1,0:nybc+1))
      end select

c     Find BC update

      if (nvar == 1) then
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,nvar),dim,loc)
        case(EQU)
          call extrapolateBC(array(:,:,:,nvar),array0(:,:,:,nvar)
cc     .                      ,dim,loc,1)
     .                      ,dim,loc,order)
        case(DIR)
          call scalarDirichletBC(array(:,:,:,nvar),zeros,ieq,dim,loc,1)
        case(NEU)
          call scalarNeumannBC(array(:,:,:,nvar),ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,nvar),dim,loc)
        case(SP)
          call scalarSingularBC(array(:,:,:,nvar))
        case(IFC)
          call interpolateToFace(array(:,:,:,nvar),dim,loc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC'//int2char(bctype)//' not implemented')
        end select
      else
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,ivar),dim,loc)
        case(EQU)
          call vectorDirichletBC(ivar,array,array0,ieq,dim,loc,order)
cc          call vectorDirichletBC(ivar,array,array0,ieq,dim,loc,1)
        case(DIR)
          call vectorDirichletBC(ivar,array,vzeros,ieq,dim,loc,1)
        case(NEU)
          call vectorNeumannBC(ivar,array,ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,ivar),dim,loc)
        case(SP)
          call vectorSingularBC(ivar,array(:,:,:,ivar),cov)
        case(IFC)
          call interpolateToFace(array(:,:,:,ivar),dim,loc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC'//int2char(bctype)//' not implemented')
        end select
      endif

      deallocate(rhs)

c     End

      end subroutine fillGhostNodes

c     periodicBC
c     #################################################################
      subroutine periodicBC(array,dim,loc)
c     -----------------------------------------------------------------
c     Imposes periodic BC. On input:
c        * array -> 3D array to impose BC's on.
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      if (np == 1 .or. asm) then

#if 1
      call findBCLoopLimits(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ibc)
      case (1)                  !x0
       if (.not.asm_dir(1))
     . array(0     ,jmin:jmax,kmin:kmax)=array(nxbc,jmin:jmax,kmin:kmax)
      case (2)                  !x1                                       
       if (.not.asm_dir(1))
     . array(nxbc+1,jmin:jmax,kmin:kmax)=array(1   ,jmin:jmax,kmin:kmax)
      case (3)                  !y0                                       
       if (.not.asm_dir(2))
     . array(imin:imax,0     ,kmin:kmax)=array(imin:imax,nybc,kmin:kmax)
      case (4)                  !y1                                       
       if (.not.asm_dir(2))
     . array(imin:imax,nybc+1,kmin:kmax)=array(imin:imax,1   ,kmin:kmax)
      case (5)                  !z0                                       
       if (.not.asm_dir(3))
     . array(imin:imax,jmin:jmax,0     )=array(imin:imax,jmin:jmax,nzbc)
      case (6)                  !z1                                       
       if (.not.asm_dir(3))
     . array(imin:imax,jmin:jmax,nzbc+1)=array(imin:imax,jmin:jmax,1   )
      case default
       call pstop('periodicBC'
     .           ,'Boundary '//int2char(ibc)//' non existent')
      end select

#else

      write (*,*) 'DIAG -- periodicBC'   !Need to change imposeBConVector

      select case (ibc)
      case (1)                  !x0
       array(0     ,:,:)=array(nxbc,:,:)
      case (2)                  !x1
       array(nxbc+1,:,:)=array(1   ,:,:)
      case (3)                  !y0
       array(:,0     ,:)=array(:,nybc,:)
      case (4)                  !y1
       array(:,nybc+1,:)=array(:,1   ,:)
      case (5)                  !z0
       array(:,:,0     )=array(:,:,nzbc)
      case (6)                  !z1
       array(:,:,nzbc+1)=array(:,:,1   )
      case default
       write (*,*) 'Boundary',ibc,' non existent'
       stop
      end select

#endif

      endif

c     End program

      end subroutine periodicBC

c     extrapolateBC
c     #######################################################################
      subroutine extrapolateBC(array,array0,dim,loc,order)
c     -----------------------------------------------------------------
c     Fills ghost nodes by extrapolation across relevant boundary.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: order,dim,loc
        real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .               ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

        integer    :: imax,imin,jmax,jmin,kmax,kmin,ibc

        integer    :: i,j,k,ig,jg,kg
        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

c     Begin program

        ibc = (1+loc)+2*(dim-1)

        rhs = 0d0

c     Setup grid aliases

        xx => grid_params%xx
        yy => grid_params%yy
        zz => grid_params%zz

        dxh => grid_params%dxh
        dyh => grid_params%dyh
        dzh => grid_params%dzh

        dx => grid_params%dx
        dy => grid_params%dy
        dz => grid_params%dz

c     Extrapolate

        call findBCLoopLimits(dim,loc
     .                       ,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                       ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

        if (.not.symm) then
          
          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

                select case (ibc)
                case (1)

                  rhs(j,k) = 
     .              quad_int(xx(ig-1)+dxh(ig-1),xx(ig),xx(ig+1),xx(ig+2)
     .                       ,array0(i-1,j,k),array(i,j,k)
     .                       ,array (i+1,j,k),array(i+2,j,k)
     .                       ,xx(ig-1),order )
                case (2)

                  rhs(j,k) =
     .              quad_int(xx(ig+1)-dxh(ig+1),xx(ig),xx(ig-1),xx(ig-2)
     .                        ,array0(i+1,j,k),array(i,j,k)
     .                        ,array (i-1,j,k),array(i-2,j,k)
     .                        ,xx(ig+1),order )
                case (3)

                  rhs(i,k) =
     .              quad_int(yy(jg-1)+dyh(jg-1),yy(jg),yy(jg+1),yy(jg+2)
     .                        ,array0(i,j-1,k),array(i,j,k)
     .                        ,array (i,j+1,k),array(i,j+2,k)
     .                        ,yy(jg-1),order )
                case (4)

                  rhs(i,k) =
     .              quad_int(yy(jg+1)-dyh(jg+1),yy(jg),yy(jg-1),yy(jg-2)
     .                        ,array0(i,j+1,k),array(i,j,k)
     .                        ,array (i,j-1,k),array(i,j-2,k)
     .                        ,yy(jg+1),order )
                case (5)

                  rhs(i,j) =
     .              quad_int(zz(kg-1)+dzh(kg-1),zz(kg),zz(kg+1),zz(kg+2)
     .                        ,array0(i,j,k-1),array(i,j,k)
     .                        ,array (i,j,k+1),array(i,j,k+2)
     .                        ,zz(kg-1),order )
                case (6)

                  rhs(i,j) =
     .              quad_int(zz(kg+1)-dzh(kg+1),zz(kg),zz(kg-1),zz(kg-2)
     .                        ,array0(i,j,k+1),array(i,j,k)
     .                        ,array (i,j,k-1),array(i,j,k-2)
     .                        ,zz(kg+1),order )
                end select
              enddo
            enddo
          enddo

        else

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

              select case (ibc)
              case (1)

                rhs(j,k) = -array(i+1,j,k)

              case (2)

                rhs(j,k) = -array(i-1,j,k)

              case (3)

                rhs(i,k) = -array(i,j+1,k)

              case (4)

                rhs(i,k) = -array(i,j-1,k)

              case (5)

                rhs(i,j) = -array(i,j,k+1)

              case (6)

                rhs(i,j) = -array(i,j,k-1)

              end select

              enddo
            enddo
          enddo

        endif

c       Update BC ghost nodes

        select case (ibc)
        case (1)                !x0
         array(0     ,jmin:jmax,kmin:kmax) = rhs(jmin:jmax,kmin:kmax)
        case (2)               !x1
         array(nxbc+1,jmin:jmax,kmin:kmax) = rhs(jmin:jmax,kmin:kmax)
        case (3)                !y0
         array(imin:imax,0     ,kmin:kmax) = rhs(imin:imax,kmin:kmax)
        case (4)                !y1
         array(imin:imax,nybc+1,kmin:kmax) = rhs(imin:imax,kmin:kmax)
        case (5)                !z0
         array(imin:imax,jmin:jmax,0     ) = rhs(imin:imax,jmin:jmax)
        case (6)                !z1
         array(imin:imax,jmin:jmax,nzbc+1) = rhs(imin:imax,jmin:jmax)
        case default
         messg = 'Boundary '//int2char(ibc)//' non existent'
         call pstop('extrapolate',messg)
        end select

      end subroutine extrapolateBC

c     interpolateToFace
c     #################################################################
      subroutine interpolateToFace(array,dim,loc)
c     -----------------------------------------------------------------
c     Interpolates values to face for EQU BC. On input:
c        * array -> 3D array to impose BC's on.
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

cc      if (np == 1 .or. asm) then

      call findBCLoopLimits(dim,loc,iimnm,iimxp,jjmnm,jjmxp,kkmnm,kkmxp
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      select case (ibc)
      case (1)                  !x0
cc       if (.not.asm_dir(1))
cc     . array(0     ,jmin:jmax,kmin:kmax)=
       array(0     ,jmin:jmax,kmin:kmax)=
     .      0.5*(array(0,jmin:jmax,kmin:kmax)
     .          +array(1,jmin:jmax,kmin:kmax))
      case (2)                  !x1                                       
cc       if (.not.asm_dir(1))
cc     . array(nxbc+1,jmin:jmax,kmin:kmax)=
       array(nxbc+1,jmin:jmax,kmin:kmax)=
     .      0.5*(array(nxbc+1,jmin:jmax,kmin:kmax)
     .          +array(nxbc  ,jmin:jmax,kmin:kmax))
      case (3)                  !y0                                       
cc       if (.not.asm_dir(2))
cc     . array(imin:imax,0     ,kmin:kmax)=
       array(imin:imax,0     ,kmin:kmax)=
     .      0.5*(array(imin:imax,0,kmin:kmax)
     .          +array(imin:imax,1,kmin:kmax))
      case (4)                  !y1                                       
cc       if (.not.asm_dir(2))
cc     . array(imin:imax,nybc+1,kmin:kmax)=
       array(imin:imax,nybc+1,kmin:kmax)=
     .      0.5*(array(imin:imax,nybc+1,kmin:kmax)
     .          +array(imin:imax,nybc  ,kmin:kmax))
      case (5)                  !z0                                       
cc       if (.not.asm_dir(3))
cc     . array(imin:imax,jmin:jmax,0     )=
       array(imin:imax,jmin:jmax,0     )=
     .      0.5*(array(imin:imax,jmin:jmax,0)
     .          +array(imin:imax,jmin:jmax,1))
      case (6)                  !z1                                       
cc       if (.not.asm_dir(3))
cc     . array(imin:imax,jmin:jmax,nzbc+1)=
       array(imin:imax,jmin:jmax,nzbc+1)=
     .      0.5*(array(imin:imax,jmin:jmax,nzbc+1)
     .          +array(imin:imax,jmin:jmax,nzbc  ))
      case default
       call pstop('interpolateToFace'
     .           ,'Boundary '//int2char(ibc)//' non existent')
      end select

cc      endif

c     End program

      end subroutine interpolateToFace

c     fillFaces
c     #################################################################
      subroutine fillFaces(array,dim,loc)
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

      integer    :: dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,bcond(6)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      bcond = 0
      bcond(ibc) = EXT

c     Determine boundary limits

      if (.not.selectFace(ibc)) return

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (bc_debug) then
cc        write (*,*) 'Proc',my_rank,iimin,iimax,jjmin,jjmax,kkmin,kkmax
        write(*,*) '  fillFaces: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
      endif

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill face i=0

        if (imin == 1) then

          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibc,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j,k),bcond)
            enddo
          enddo

        endif

      case (2) !On face #2, fill face i=nx+1

        if (imax == nxbc) then

          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibc,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j,k),bcond)
            enddo
          enddo

        endif

      case (3) !On face #3, fill face j=0

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          do k=kmin,kmax
            do i=imin,imax
                call extrapolateToCorner(ibc,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i,jmn:jmx,k),bcond)
            enddo
          enddo

        endif

      case (4) !On face #4, fill face j=ny+1

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          do k=kmin,kmax
            do i=imin,imax
              call extrapolateToCorner(ibc,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i,jmn:jmx,k),bcond)
            enddo
          enddo

        endif

      case (5) !On face #5, fill face k=0

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibc,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i,j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      case (6) !On face #6, fill face k=nz+1

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibc,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i,j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      end select

c     End program

      end subroutine fillFaces

c     fillEdges
c     #################################################################
      subroutine fillEdges(array,array0,dim,loc,bcond)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: dim,loc,bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

cc      if (.not.selectFace(ibc)) return

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (bc_debug) then
cc        write (*,*) 'Proc',my_rank,iimin,iimax,jjmin,jjmax,kkmin,kkmax
        write(*,*) '  fillEdges: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
      endif

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill edges y=jmin-1,y=jmax+1, vary z

        if (imin == 1) then
          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1.and.(selectFace(ibc).or.selectFace(3))) then  !Check faces #1,#3
            jto   = jmin-1
            jfrom = min(jmin+order,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(3)==EQU) .and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc.and.(selectFace(ibc).or.selectFace(4))) then  !Check faces #1,#4
            jto   = jmax+1
            jfrom = max(jmax-order,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

        endif

      case (2) !On face #2, fill edges y=jmin-1,y=jmax+1, vary z

        if (imax == nxbc) then
          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1.and.(selectFace(ibc).or.selectFace(3))) then  !Check faces #2,#3
            jto   = jmin-1
            jfrom = min(jmin+order,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(3)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc.and.(selectFace(ibc).or.selectFace(4))) then  !Check faces #2,#4
            jto   = jmax+1
            jfrom = max(jmax-order,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

        endif

      case (3) !On face #3, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1.and.(selectFace(ibc).or.selectFace(5))) then  !Check faces #3,#5

            kto   = kmin-1
            kfrom = min(kmin+order,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc.and.(selectFace(ibc).or.selectFace(6))) then  !Check faces #3,#6

            kto   = kmax+1
            kfrom = max(kmax-order,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (4) !On face #4, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1.and.(selectFace(ibc).or.selectFace(5))) then  !Check faces #4,#5

            kto   = kmin-1
            kfrom = min(kmin+order,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc.and.(selectFace(ibc).or.selectFace(6))) then  !Check faces #4,#6

            kto   = kmax+1
            kfrom = max(kmax-order,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (5) !On face #5, fill edges x=imin-1,x=imax+1, vary y

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1.and.(selectFace(ibc).or.selectFace(1))) then  !Check faces #5,#1

            ito   = imin-1
            ifrom = min(imin+order,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc.and.(selectFace(ibc).or.selectFace(2))) then  !Check faces #5,#2

            ito   = imax+1
            ifrom = max(imax-order,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (6) !On face #6, fill edges x=imin-1,x=imax+1, vary y

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1.and.(selectFace(ibc).or.selectFace(1))) then  !Check faces #6,#1

            ito   = imin-1
            ifrom = min(imin+order,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc.and.(selectFace(ibc).or.selectFace(2))) then  !Check faces #6,#2

            ito   = imax+1
            ifrom = max(imax-order,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      end select

c     End program

      end subroutine fillEdges

c     fillCorners
c     #################################################################
      subroutine fillCorners(array,array0,dim,loc,bcond)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: dim,loc,bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

cc      if (.not.selectFace(ibc)) return

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

c     Select corner

      if (      imin == 1  .and. jmin == 1  .and. kmin == 1
     .    .and.(selectFace(1).or.selectFace(3).or.selectFace(5)) ) then   !Check concurrent faces
        ifrom = min(imin+order,nxbc)
        jfrom = min(jmin+order,nybc)
        kfrom = min(kmin+order,nzbc)
        ito   = imin-1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write (*,*) '  fillCorners: Proc',my_rank
     .               ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .          ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imin == 1  .and. jmax == nybc .and. kmin == 1
     .    .and.(selectFace(1).or.selectFace(4).or.selectFace(5))) then

        ifrom = min(imin+order,nxbc)
        jfrom = max(jmax-order,1   )
        kfrom = min(kmin+order,nzbc)
        ito   = imin-1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .           ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmin == 1  .and. kmin == 1
     .    .and.(selectFace(2).or.selectFace(3).or.selectFace(5))) then

        ifrom = max(imax-order,1   )
        jfrom = min(jmin+order,nybc)
        kfrom = min(kmin+order,nzbc)
        ito   = imax+1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmax == nybc .and. kmin == 1
     .    .and.(selectFace(2).or.selectFace(4).or.selectFace(5))) then

        ifrom = max(imax-order,1   )
        jfrom = max(jmax-order,1   )
        kfrom = min(kmin+order,nzbc)
        ito   = imax+1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imin == 1  .and. jmin == 1  .and. kmax == nzbc
     .    .and.(selectFace(1).or.selectFace(3).or.selectFace(6))) then

        ifrom = min(imin+order,nxbc)
        jfrom = min(jmin+order,nybc)
        kfrom = max(kmax-order,1   )
        ito   = imin-1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imin == 1  .and. jmax == nybc .and. kmax == nzbc
     .    .and.(selectFace(1).or.selectFace(4).or.selectFace(6))) then

        ifrom = min(imin+order,nxbc)
        jfrom = max(jmax-order,1   )
        kfrom = max(kmax-order,1   )
        ito   = imin-1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmin == 1  .and. kmax == nzbc
     .    .and.(selectFace(2).or.selectFace(3).or.selectFace(6))) then

        ifrom = max(imax-order,1   )
        jfrom = min(jmin+order,nybc)
        kfrom = max(kmax-order,1   )
        ito   = imax+1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (     imax == nxbc .and. jmax == nybc .and. kmax == nzbc
     .    .and.(selectFace(2).or.selectFace(4).or.selectFace(6))) then

        ifrom = max(imax-order,1 )
        jfrom = max(jmax-order,1 )
        kfrom = max(kmax-order,1 )
        ito   = imax+1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if (bc_debug) then
          write(*,*) '  fillCorners: Proc',my_rank
     .              ,'; Corner =',ito,jto,kto
        endif

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)

      endif

c     End program

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .                              ,imx,jmx,kmx,array,bcond,array0)
c     -----------------------------------------------------------------
c     Routine that extrapolates information to edge or vertex corners.
c     It is the workhorse of the fillEdge and fillCorner routines.
c     In the call sequence:
c        * ibc: boundary identifier (integer)
c        * ito,jto,kto: coordinates of ghost cell to be filled (integer)
c        * imn,jmn,kmn: node coordinates identifying start of
c                       interpolation region (integer)
c        * imx,jmx,kmx: node coordinates identifying end of interpolation
c                       region (integer)
c        * array: local array defined on interpolation region, including
c                 ghost cell (ito,jto,kto) (real)
c     The order of extrapolation 'order' is defined in the module header.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ibc,ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx,bcond(6)
      real(8)    :: array(imn:imx,jmn:jmx,kmn:kmx)
      real(8),optional,intent(IN) :: array0(imn:imx,jmn:jmx,kmn:kmx)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,inx,iny,inz
     .             ,order_int(3),offset
     .             ,imni,imxi,jmni,jmxi,kmni,kmxi,num
     .             ,ibcl(3),i_bc,i1,i2

      real(8)    :: valx(1),valy(1),valz(1),xto(1),yto(1),zto(1)
     .             ,smooth(3)

      logical    :: flg(3),select_bc(3),have_equ

c     Extrapolation

      real(8),allocatable,dimension(:) :: xx,yy,zz

c     Begin program

      have_equ = PRESENT(array0)

c     Identify directions for extrapolation.

      inx = imx-imn
      iny = jmx-jmn
      inz = kmx-kmn

      flg = .true.

      !Order in each direction
      order_int(1) = min(order,max(inx-1,0))
      order_int(2) = min(order,max(iny-1,0))
      order_int(3) = min(order,max(inz-1,0))

      !Smoothness indicator
      smooth = 0d0
      if (order >= 2) then
        if (order_int(1) == order) then
          if (ito == imn) then
            smooth(1) = array(ito+1,jto,kto)-array(ito+2,jto,kto)
          else
            smooth(1) = array(ito-1,jto,kto)-array(ito-2,jto,kto)
          endif
        endif

        if (order_int(2) == order) then
          if (jto == jmn) then
            smooth(2) = array(ito,jto+1,kto)-array(ito,jto+2,kto)
          else                                              
            smooth(2) = array(ito,jto-1,kto)-array(ito,jto-2,kto)
          endif
        endif

        if (order_int(3) == order) then
          if (kto == kmn) then
            smooth(3) = array(ito,jto,kto+1)-array(ito,jto,kto+2)
          else                                                  
            smooth(3) = array(ito,jto,kto-1)-array(ito,jto,kto-2)
          endif
        endif

        smooth = abs(smooth)
      endif

      !Do not interpolate **along** direction if
      !   1) not enough points,
      if (inx == 0) flg(1) = .false.
      if (iny == 0) flg(2) = .false.
      if (inz == 0) flg(3) = .false.

      !Select boundaries involved 
      ![ibc = (1+loc)+2*(dim-1), with dim=1,2,3, and loc=0(left),1(right)]

      ibcl = 0

      if (flg(1) .and. flg(2) .and. flg(3)) then !Corner

        if (ito == 0 .and. jto == 0 .and. kto == 0) then
          ibcl(1) = 1
          ibcl(2) = 3
          ibcl(3) = 5
        elseif (ito == 0 .and. jto == 0 .and. kto /= 0) then
          ibcl(1) = 1
          ibcl(2) = 3
          ibcl(3) = 6
        elseif (ito == 0 .and. jto /= 0 .and. kto == 0) then
          ibcl(1) = 1
          ibcl(2) = 4
          ibcl(3) = 5
        elseif (ito /= 0 .and. jto == 0 .and. kto == 0) then
          ibcl(1) = 2
          ibcl(2) = 3
          ibcl(3) = 5
        elseif (ito /= 0 .and. jto /= 0 .and. kto == 0) then
          ibcl(1) = 2
          ibcl(2) = 4
          ibcl(3) = 5
        elseif (ito /= 0 .and. jto == 0 .and. kto /= 0) then
          ibcl(1) = 2
          ibcl(2) = 3
          ibcl(3) = 6
        elseif (ito == 0 .and. jto /= 0 .and. kto /= 0) then
          ibcl(1) = 1
          ibcl(2) = 4
          ibcl(3) = 6
        elseif (ito /= 0 .and. jto /= 0 .and. kto /= 0) then
          ibcl(1) = 2
          ibcl(2) = 4
          ibcl(3) = 6
        endif

      elseif (   (flg(1) .and. flg(2))
     $       .or.(flg(1) .and. flg(3))
     $       .or.(flg(2) .and. flg(3))) then !Edge

        select case(ibc)
        case (1,2)
          ibcl(1) = ibc
          if (jto == 0) then
            ibcl(2) = 3
          else
            ibcl(2) = 4
          endif
        case (3,4)
          ibcl(2) = ibc
          if (kto == 0) then
            ibcl(3) = 5
          else
            ibcl(3) = 6
          endif
        case (5,6)
          ibcl(3) = ibc
          if (ito == 0) then
            ibcl(1) = 1
          else
            ibcl(1) = 2
          endif
        end select

      else !Face

        select case(ibc)
        case (1,2)
          ibcl(1) = ibc
        case (3,4)
          ibcl(2) = ibc
        case (5,6)
          ibcl(3) = ibc
        end select

      endif
        
      !Treat special boundaries

      !Identify special boundaries to interpolate **along**
      !   1) SP boundary
      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc == 0) cycle
        select_bc(i) =    (bcSP(ibc=i_bc) .and. bcond(i_bc)==SP)
cc     .                .or.(bcond(i_bc) == EQU)  !Removed 12/13/06
      enddo

      call select_along

      !Identify special boundaries to interpolate **across**
      !   1) periodic boundary
      !   2) Symmetry boundary
      !   3) equilibrium boundary
      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc == 0) cycle
        select_bc(i) =    (abs(bcond(i_bc)) == SYM)
     .                .or.(    bcond(i_bc)  == PER)
     .                .or.(    bcond(i_bc)  == EQU) !Added 12/13/06
      enddo

      call select_across

      !Check some direction for interpolation exists

      if ((.not.flg(1)).and.(.not.flg(2)).and.(.not.flg(3))) then
        write (*,*) 'Extrapolation directions',flg
        write (*,*) 'Boundaries involved',ibcl
        write (*,*) 'No extrapolation directions left'
        call pstop('extrapolateToCorner','')
      endif

c     Initialize accumulators

      valx = 0d0
      valy = 0d0
      valz = 0d0

      num = 0

c     X extrapolation

      if (flg(1)) then

        allocate(xx(inx))

        if (ito == imn) then
          call getMGmap(imn+1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          xx = grid_params%xx(ig:ig+inx-1)
          imni = imn+1
          imxi = imx
        else
          call getMGmap(imx-1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          xx = grid_params%xx(ig-(inx-1):ig)
          imni = imn
          imxi = imx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        xto = grid_params%xx(ig)

        if (abs(bcond(ibcl(1))) == SYM) then  !Exception for symmetry BC
          if (ito == imn) then
            valx(1) = sign(1,bcond(1))*array(ito+2,jto,kto)
          else
            valx(1) = sign(1,bcond(1))*array(ito-2,jto,kto)
          endif
        else
          call IntDriver1d(inx,xx,array(imni:imxi,jto,kto)
     .                  ,1,xto,valx,order_int(1))
        endif

        num = num + 1

        deallocate(xx)
      endif

c     Y extrapolation

      if (flg(2)) then

        allocate(yy(iny))

        if (jto == jmn) then
          call getMGmap(imn,jmn+1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          yy = grid_params%yy(jg:jg+iny-1)
          jmni = jmn+1
          jmxi = jmx
        else
          call getMGmap(imn,jmx-1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          yy = grid_params%yy(jg-(iny-1):jg)
          jmni = jmn
          jmxi = jmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        yto = grid_params%yy(jg)

        if (abs(bcond(ibcl(2))) == SYM) then  !Exception for symmetry BC
          if (jto == jmn) then
            valy(1) = sign(1,bcond(1))*array(ito,jto+2,kto)
          else                                                   
            valy(1) = sign(1,bcond(1))*array(ito,jto-2,kto)
          endif
        else
          call IntDriver1d(iny,yy,array(ito,jmni:jmxi,kto)
     .                  ,1,yto,valy,order_int(2))
        endif

        num = num + 1

        deallocate(yy)

      endif

c     Z extrapolation

      if (flg(3)) then

        allocate(zz(inz))

        if (kto == kmn) then
          call getMGmap(imn,jmn,kmn+1,igxbc,igybc,igzbc,ig,jg,kg)
          zz = grid_params%zz(kg:kg+inz-1)
          kmni = kmn+1
          kmxi = kmx
        else
          call getMGmap(imn,jmn,kmx-1,igxbc,igybc,igzbc,ig,jg,kg)
          zz = grid_params%zz(kg-(inz-1):kg)
          kmni = kmn
          kmxi = kmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        zto = grid_params%zz(kg)

        if (abs(bcond(ibcl(3))) == SYM) then  !Exception for symmetry BC
          if (kto == kmn) then
            valz(1) = sign(1,bcond(1))*array(ito,jto,kto+2)
          else                                                   
            valz(1) = sign(1,bcond(1))*array(ito,jto,kto-2)
          endif
        else
          call IntDriver1d(inz,zz,array(ito,jto,kmni:kmxi)
     .                  ,1,zto,valz,order_int(3))
        endif

        num = num + 1

        deallocate(zz)
      endif

c     Average extrapolation results and fill ghost cell

      array(ito,jto,kto) = (valx(1) + valy(1) + valz(1))/num

c     End program

      contains

c     select_along
c     ######################################################################
      subroutine select_along
c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs ALONG
c     selected boundary only.
c     -----------------------------------------------------------------------

      do i=1,3

        !Identify other directions
        i1 = mod(i+1,3)
        i2 = mod(i+2,3)
        if (i1 == 0) i1 = 3
        if (i2 == 0) i2 = 3

        !Identify cases in which an edge has two or more selected directions:
        !  can't eliminate both interpolation directions simultaneously;
        !  resolve indetermination by choosing smoothest direction
        if (select_bc(i) .and. (select_bc(i1) .or. select_bc(i2))) then

cc          if (select_bc(i1).and.select_bc(i2)) then
cc            if (order_int(i1) <= order_int(i2)) then
cc              select_bc(i1) = .false.
cc            else
cc              select_bc(i2) = .false.
cc            endif
cc          endif
cc
cc          if (select_bc(i1)) then
cc            if (order_int(i) <= order_int(i1)) then
cc              cycle
cc            else
cc              flg(i) = .not.select_bc(i)
cc            endif
cc          else
cc            if (order_int(i) <= order_int(i2)) then
cc              cycle
cc            else
cc              flg(i) = .not.select_bc(i)
cc            endif
cc          endif

          if (select_bc(i1).and.select_bc(i2)) then
            if (smooth(i1) <= smooth(i2)) then
              select_bc(i1) = .false.
            else
              select_bc(i2) = .false.
            endif
          endif

          if (select_bc(i1)) then
            if (smooth(i) <= smooth(i1)) then
              cycle
            else
              flg(i) = .not.select_bc(i)
            endif
          else
            if (smooth(i) <= smooth(i2)) then
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
c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs ACROSS
c     selected boundary only.
c     -----------------------------------------------------------------------

      do i=1,3

        !Identify other directions
        i1 = mod(i+1,3)
        i2 = mod(i+2,3)
        if (i1 == 0) i1 = 3
        if (i2 == 0) i2 = 3

        if (select_bc(i) .and. flg(i)) then
          flg(i1) = .false.
          flg(i2) = .false.
        endif

      enddo

      end subroutine select_across

      end subroutine extrapolateToCorner

c     quad_int
c     #################################################################
      real(8) function quad_int(x0,x1,x2,x3,y0,y1,y2,y3,x,order)
     .        result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: order
      real(8)    :: x0,x1,x2,x3,y0,y1,y2,y3,x

c     Local variables

c     Begin program

      select case (order)
      case (3)
        y = y0*(x-x1)*(x-x2)*(x-x3)/(x0-x1)/(x0-x2)/(x0-x3)
     .     +y1*(x-x0)*(x-x2)*(x-x3)/(x1-x0)/(x1-x2)/(x1-x3)
     .     +y2*(x-x0)*(x-x1)*(x-x3)/(x2-x0)/(x2-x1)/(x2-x3)
     .     +y3*(x-x0)*(x-x1)*(x-x2)/(x3-x0)/(x3-x1)/(x3-x2)
      case (2)
        y = y0*(x-x1)*(x-x2)/(x0-x1)/(x0-x2)
     .     +y1*(x-x0)*(x-x2)/(x1-x0)/(x1-x2)
     .     +y2*(x-x0)*(x-x1)/(x2-x0)/(x2-x1)
      case (1)
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      case (0)
        y = y0
      case default
        messg = 'Order of interpolation not implemented'
        call pstop('quad_int',messg)
      end select

c     End program

      end function quad_int

      end module imposeBCinterface

c module setMGBC_interface
c######################################################################
      module setMGBC_interface

        use grid_mg

        use imposeBCinterface

      contains

c     setMGBC
c     ###############################################################
      subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd,arr_cov
     .                  ,arr0,icomp,is_cnv,is_vec,result_is_vec,iorder)
c     ---------------------------------------------------------------
c     Interfaces BC routines with MG code, for preconditioning.
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
c       *icomp: integer specifying magnitude of interest (if not
c               provided, use last used value).
c       *is_cnv: logical variable indicating whether array
c                contains covariant or contravariant components
c                (if not provided, use last used value).
c       *is_vec: logical variable specifying whether array
c          represents a curvilinear vector or a set of scalars
c          (if not provided, use last used value):
c             -If is_vec=.true., array is taken as curvilinear vector
c              components (covariant or contravariant, according
c              to variable is_cnv).
c             -If is_vec=.false., array is taken as coupled scalars
c              (also Cartesian vector components).
c       *result_is_vec: logical variable specifying whether output
c                       in array should be curvilinear (.true.) or
c                       not (default=is_vec)
c       *iorder: order of inter/extrapolation in BC routines (linear
c                by default).
c
c     On output, array is returned in the same representation.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: nnx,nny,nnz,neq,bcnd(6,neq),iig,gpos

      real(8)    :: array(0:nnx+1,0:nny+1,0:nnz+1,neq)

      real(8),optional,intent(INOUT) ::
     .                         arr_cov(0:nnx+1,0:nny+1,0:nnz+1,neq)
      real(8),optional,intent(IN) ::
     .                         arr0   (0:nnx+1,0:nny+1,0:nnz+1,neq)

      integer   ,optional,intent(IN) :: icomp,iorder

      logical   ,optional,intent(IN) :: is_cnv,is_vec,result_is_vec

c     Local variables

      integer    :: i,j,k,stencil_width,order,ieq
     .             ,imng,imxg,jmng,jmxg,kmng,kmxg
     .             ,imnl,imxl,jmnl,jmxl,kmnl,kmxl

      integer   ,save :: ivar  !Needed to "remember' ivar from previous calls
                               !when icomp is not provided.

      logical,save    :: iscnv,isvec !Needed to "remember' these from previous calls
                                     !when they are not provided.

      logical         :: res_is_vec,fpointers

      real(8)    :: car(3),curv(3)

c     Begin program

c     Optional arguments

      if (PRESENT(icomp))  ivar  = icomp

      if (PRESENT(is_cnv)) iscnv = is_cnv

      if (PRESENT(is_vec)) isvec = is_vec

      if (PRESENT(result_is_vec)) then
        res_is_vec = result_is_vec
      else
        res_is_vec = isvec
      endif

      if (PRESENT(iorder)) then
        order = iorder
      else
        order = 2
      endif

c     Consistency check

      if (     grid_params%nxv(iig) /= nnx
     .    .or. grid_params%nyv(iig) /= nny
     .    .or. grid_params%nzv(iig) /= nnz) then
        write (*,*) 'Grid sizes do not agree in setMGBC'
        write (*,*) 'Aborting...'
        stop
      endif

      if (.not.isvec) iscnv=.true.

c     Allocate MG pointers (if not allocated earlier)

      call allocPointers(neq,fpointers)

c     Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(gpos,nnx,nny,nnz,iig,imnl,imxl,jmnl,jmxl,kmnl,kmxl)

c     Check LOCAL limits (return if not close to local domain boundaries)

      if (     (imnl > 1 .and. imxl < nnx)
     .    .and.(jmnl > 1 .and. jmxl < nny)
     .    .and.(kmnl > 1 .and. kmxl < nnz)) return

c     Find GLOBAL limits (required for BC routine setBC)

cc      stencil_width = 1
cc
cc      imng = max(imnl-stencil_width,1)   + grid_params%ilo(iig)-1
cc      imxg = min(imxl+stencil_width,nnx) + grid_params%ilo(iig)-1
cc      jmng = max(jmnl-stencil_width,1)   + grid_params%jlo(iig)-1
cc      jmxg = min(jmxl+stencil_width,nny) + grid_params%jlo(iig)-1
cc      kmng = max(kmnl-stencil_width,1)   + grid_params%klo(iig)-1
cc      kmxg = min(kmxl+stencil_width,nnz) + grid_params%klo(iig)-1

      imng = imnl + grid_params%ilo(iig)-1
      imxg = imxl + grid_params%ilo(iig)-1
      jmng = jmnl + grid_params%jlo(iig)-1
      jmxg = jmxl + grid_params%jlo(iig)-1
      kmng = kmnl + grid_params%klo(iig)-1
      kmxg = kmxl + grid_params%klo(iig)-1

c     Select operation

      select case (neq)
      case(1)

        allocate(vv0(0:nnx+1,0:nny+1,0:nnz+1,1))

        if (PRESENT(arr0)) then
cc          vv0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc     .         = arr0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
          vv0 = arr0
        else
cc          vv0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
          vv0 = 0d0
        endif

        call setBC(ivar,nnx,nny,nnz,array(:,:,:,1),vv0(:,:,:,1)
     .            ,bcnd(:,neq),iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,iorder=order)

        deallocate(vv0)

      case(2)

        allocate(vv0(0:nnx+1,0:nny+1,0:nnz+1,2))

        if (PRESENT(arr0)) then
          vv0 = arr0
        else
          vv0 = 0d0
        endif

        do ieq=1,neq
          call setBC(ivar,nnx,nny,nnz,array(:,:,:,ieq),vv0(:,:,:,ieq)
     .            ,bcnd(:,ieq),iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,iorder=order)
        enddo

        deallocate(vv0)

      case(3)

        allocate(v_cnv (0:nnx+1,0:nny+1,0:nnz+1,neq)
     .          ,v_cov (0:nnx+1,0:nny+1,0:nnz+1,neq)
     .          ,vv0   (0:nnx+1,0:nny+1,0:nnz+1,neq))

        v_cnv(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
        v_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
cc        v_cnv = 0d0
cc        v_cov = 0d0

        if (PRESENT(arr0)) then
          vv0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
     .       = arr0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc          vv0 = arr0
        else
          vv0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
cc          vv0 = 0d0
        endif

        if (.not.isvec) then

          do k = kmnl-1,kmxl+1
            do j = jmnl-1,jmxl+1
              do i = imnl-1,imxl+1
                car = array(i,j,k,:)
                call transformVectorToCurvilinear
     .               (i,j,k,iig,iig,iig
     .               ,car(1),car(2),car(3)
     .               ,.not.iscnv
     .               ,array(i,j,k,1),array(i,j,k,2),array(i,j,k,3))
              enddo
            enddo
          enddo

        endif

        if (PRESENT(arr_cov)) then
cc          v_cnv(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc     .         = array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc          v_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc     .         = arr_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
          v_cnv = array
          v_cov = arr_cov
        else
          if (iscnv) then
cc            v_cnv(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc     .           = array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
            v_cnv = array
          else
cc            v_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc     .           = array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
            v_cov = array
          endif
        endif

        call setBC(ivar,neq,nnx,nny,nnz,v_cnv,v_cov,vv0,bcnd
     .            ,iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,is_cnv=iscnv,is_vec=.true.,iorder=order)

        if (isvec .or. res_is_vec) then

          if (PRESENT(arr_cov)) then
            array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
     .           = v_cnv(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
            arr_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
     .           = v_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc            array = v_cnv
cc            arr_cov = v_cov
          else
            if (iscnv) then
              array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
     .             = v_cnv(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc              array = v_cnv
            else
              array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
     .             = v_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc              array = v_cov
            endif
          endif

        else

          do k = kmnl-1,kmxl+1
            do j = jmnl-1,jmxl+1
              do i = imnl-1,imxl+1
                curv = v_cnv(i,j,k,:)
                call transformVectorToCartesian
     .               (i,j,k,iig,iig,iig
     .               ,curv(1),curv(2),curv(3)
     .               ,.false.
     .               ,array(i,j,k,1),array(i,j,k,2),array(i,j,k,3))
              enddo
            enddo
          enddo

        endif

        deallocate(v_cov,v_cnv,vv0)

      case default
        write (*,*) 'Number of equations not implemented in setMGBC'
        write (*,*) 'Aborting...'
        stop
      end select

c     Deallocate MG pointers (if allocated in this routine)

      call deallocPointers(fpointers)

c     End

      end subroutine setMGBC

      end module setMGBC_interface

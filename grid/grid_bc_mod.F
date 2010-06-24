c module bc_ops
c #####################################################################
      module bc_ops

c     -----------------------------------------------------------------
c     This module defines various BC-related operators
c     -----------------------------------------------------------------

        use grid_def

      contains

c     isBdry
c     #################################################################
      function isBdry(i,igr,ibc)

c     -----------------------------------------------------------------
c     Detects whether we are at a boundary or not. Boundary is identified
c     by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

        integer    :: i,igr,ibc
        logical    :: isBdry

        integer    :: ilog,itst,diml,locl

        locl = abs(mod(ibc,2)-2)-1
        diml = 1 + (ibc-1-locl)/2

        select case(diml)
        case(1)
          ilog = grid_params%ilo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nxgl(igr)
          end select
        case(2)
          ilog = grid_params%jlo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nygl(igr)
          end select
        case(3)
          ilog = grid_params%klo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nzgl(igr)
          end select
        end select

        isBdry = (i+ilog-1 == itst)

      end function isBdry

c     isBC
c     #################################################################
      function isBC(i,igr,ibc,BC)

        implicit none

        integer    :: i,igr,ibc,BC
        logical    :: isBC

        isBC = (isBdry(i,igr,ibc) .and. abs(bcond(ibc)) == BC)

      end function isBC

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        implicit none

        integer    :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point = (i+grid_params%ilo(igx)-1==1).and.bcSP()

      end function isSP

c     isSP2
c     #################################################################
      function isSP2(i,igr,dim,loc,ibc)

        implicit none

        integer    :: i,igr
        logical    :: isSP2
        integer,optional :: dim,loc,ibc

        integer    :: ibcl,ilog,itst,diml,locl

        if (PRESENT(dim) .and.PRESENT(loc)) then
          diml = dim
          locl = loc
          ibcl = (1+locl)+2*(diml-1)
        elseif (PRESENT(ibc)) then
          ibcl = ibc
        else
          ibcl = 1
        endif

        if (ibcl < 1 .or. ibcl > 6) then
          isSP2 = .false.
          return
        else
          isSP2 = isBC(i,igr,ibcl,SP)
        endif

      end function isSP2

c     isSYM
c     #################################################################
      function isSYM(i,igr,dim,loc,ibc)

        implicit none

        integer    :: i,igr
        logical    :: isSYM
        integer,optional :: dim,loc,ibc

        integer    :: ibcl,ilog,itst,diml,locl

        if (PRESENT(dim) .and.PRESENT(loc)) then
          diml = dim
          locl = loc
          ibcl = (1+locl)+2*(diml-1)
        elseif (PRESENT(ibc)) then
          ibcl = ibc
          locl = abs(mod(ibcl,2)-2)-1
          diml = 1 + (ibcl-1-locl)/2
        else
          write (*,*) 'Error in isSYM'
          write (*,*) 'Aborting...'
          stop
        endif

        isSYM = isBC(i,igr,ibcl,SYM)

      end function isSYM

      end module bc_ops

c module BCS_variables
c####################################################################
      module BCS_variables

        use oned_int

        use grid_mpi

        use grid_operations

        use math

        use bc_ops

        integer    :: nnvar,nxbc,nybc,nzbc,igxbc,igybc,igzbc

        real(8),allocatable,dimension(:,:) :: rhs

        logical :: symm=.false.,bc_debug=.false.

        integer,parameter :: VECTOR=1,SCALAR=0,AUX=-1

        real(8), pointer, dimension(:,:,:)   :: zeros,ones
        real(8), pointer, dimension(:,:,:,:) :: vzeros

      end module BCS_variables

c module singularBCinterface
c####################################################################
      module singularBCinterface

        use BCS_variables

      contains

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
cc      allocate(global_x((nyg+2*npy)*nzbc),garray(0:nyg+1,nzbc))
cc      call MPI_Allgather(array(i,0:nybc+1,nzbc)
cc     .                           ,(nybc+2)*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,global_x,(nybc+2)*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_COMM_SP,mpierr)
cc      do k=1,npy
cc        garray(nybc*(k-1):nybc*k+1,:)
cc     .     = reshape(global_x(1+(nybc+2)*nzbc*(k-1)
cc     .                         :(nybc+2)*nzbc* k)
cc     .              , (/ nybc+2,nzbc /))
cc      enddo
cc      deallocate(global_x)

cc      allocate(global_x(nyg*nzbc),garray(nyg,nzbc))
cc      call MPI_Allgather(array(i,1:nybc,1:nzbc)
cc     .                           ,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,global_x,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_COMM_SP,mpierr)
cc      do k=1,npy
cc        garray(1+nybc*(k-1):nybc*k,:)
cc     .     = reshape(global_x(1+nybc*nzbc*(k-1):nybc*nzbc*k)
cc     .              , (/ nybc,nzbc /))
cc      enddo
cc      deallocate(global_x)
      allocate(global_x(nyg*(nzbc+2)),garray(nyg,0:nzbc+1))
      call MPI_Allgather(array(i,1:nybc,0:nzbc+1)
     .                           ,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                  ,global_x,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                  ,MPI_COMM_SP,mpierr)
      do k=1,npy
        garray(1+nybc*(k-1):nybc*k,:)
     .     = reshape(global_x(1+nybc*(nzbc+2)*(k-1):nybc*(nzbc+2)*k)
     .              , (/ nybc,nzbc+2 /))
      enddo
      deallocate(global_x)
#endif

c     Element exchange in poloidal direction

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('scalarSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

      do k=0,nzbc+1
        do j=0,nybc+1
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
cc      allocate(global_x(nyg*nzbc),garray(nyg,nzbc))
cc      call MPI_Allgather(vec(i,1:nybc,1:nzbc)
cc     .                           ,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,global_x,nybc*nzbc,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_COMM_SP,mpierr)
cc      do k=1,npy
cc        garray(1+nybc*(k-1):nybc*k,:)
cc     .     = reshape(global_x(1+nybc*nzbc*(k-1):nybc*nzbc*k)
cc     .              , (/ nybc,nzbc /))
cc      enddo
cc      deallocate(global_x)

      allocate(global_x(nyg*(nzbc+2)),garray(nyg,0:nzbc+1))
      call MPI_Allgather(vec(i,1:nybc,0:nzbc+1)
     .                           ,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                  ,global_x,nybc*(nzbc+2),MPI_DOUBLE_PRECISION
     .                  ,MPI_COMM_SP,mpierr)
      do k=1,npy
        garray(1+nybc*(k-1):nybc*k,:)
     .     = reshape(global_x(1+nybc*(nzbc+2)*(k-1):nybc*(nzbc+2)*k)
     .              , (/ nybc,nzbc+2 /))
      enddo
      deallocate(global_x)
#endif

c     Element exchange in poloidal direction

      if (mod(nyg,2) /= 0 .and. nyg /= 1) then
        call pstop('vectorSingularBC'
     .            ,'Poloidal dimension not a multiple of 2')
      endif

      if (coords == 'hl2') then
        do k=0,nzbc+1
          do j=0,nybc+1
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
        do k=0,nzbc+1
          do j=0,nybc+1
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

        integer   ,private :: bctype
cc        integer   ,private :: dim,loc,ibc,i,j,k

        !Define order of interpolation for EDGE and CORNER treatment
        !(linear to preserve nearest-neighbors stencil; this is 
        ! important for colored diagonal formation)
        integer   ,private :: order=1

        real(8),allocatable,dimension(:,:,:,:) :: v_cov,v_cnv,vv0

        private :: setASMflags

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

        select case(loc)
        case(0)
          ill = min(igl - igmin + 1,1)
          ilr = ill
        case(1)
          ill = max(igr - igmin + 1,nxl)
          ilr = ill
        end select

      end subroutine selectLimits

      end subroutine findBCLoopLimits

c     do__nothing
c     #####################################################################
      function do__nothing(dim,loc,i1,i2,j1,j2,k1,k2,asm_bc)

c     ---------------------------------------------------------------------
c     Selects whether to fill or not ghost cells of a given local domain
c     based on its position in the global grid.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        logical :: do__nothing
        integer,intent(IN) :: dim,loc,i1,i2,j1,j2,k1,k2
        logical,optional :: asm_bc

c     Local variables

        integer :: ibc
        logical :: sing_face,asmbc

c     Begin program

        do__nothing = .false.

        if (PRESENT(asm_bc)) then
          asmbc = asm_bc
        else
          asmbc = asm   !Setup externally
        endif

        if (asmbc) return      !Return if additive Schwartz method is used (PC)

        select case(dim)
        case(1)
          call checkLimits(i1,i2,grid_params%nxgl(igxbc))
        case(2)
          call checkLimits(j1,j2,grid_params%nygl(igybc))
        case(3)
          call checkLimits(k1,k2,grid_params%nzgl(igzbc))
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

      if (ibc <= 0 .or. ibc > 6) then
        selectFace = .false.
        return
      endif

      loc = abs(mod(ibc,2)-2)-1
      dim = 1 + (ibc-1-loc)/2

      selectFace= 
     .     .not.do__nothing(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax)

      end function selectFace

c     setASMflags
c     ######################################################################
      subroutine setASMflags(asm_flag)

c     -----------------------------------------------------------------
c     Sets flag to indicate Additive Schwartz method.
c     -----------------------------------------------------------------

      implicit none

      logical :: asm_flag

      asm = asm_flag

      asm_dir(1) = asm.and.(npx > 1)
      asm_dir(2) = asm.and.(npy > 1)
      asm_dir(3) = asm.and.(npz > 1)

      end subroutine setASMflags

c     setup_petsc_BC
c     ######################################################################
      subroutine setup_petsc_BC

c     -----------------------------------------------------------------
c     Sets flags to indicate that BCs are fully provided by PETSc.
c     -----------------------------------------------------------------

      call setASMflags(.false.)

      end subroutine setup_petsc_BC

c     setASM_PC_BC
c     ######################################################################
      subroutine setASM_PC_BC(igx,igy,igz,bcs,pc_asm,IIBC)

c     -----------------------------------------------------------------
c     Sets boundary conditions for Additive Schwartz method PC.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igx,igy,igz,bcs(:,:)

      integer,optional :: IIBC

      logical :: pc_asm

c     Local variables

      integer :: dim,loc,ibc,igmin,igmax,jgmin,jgmax,kgmin,kgmax

      logical :: dn

c     Begin program

c     Reset Additive-Schwarz-Method flags

      if ((np == 1) .or. (.not.pc_asm)) then
        call setASMflags(.false.)
        return
      else
        call setASMflags(.true.)
      endif

c     Reset BC grid level

      igxbc = igx
      igybc = igy
      igzbc = igz

c     Find global limits of local domain

      igmin = grid_params%ilo(igx)
      igmax = grid_params%ihi(igx)
      jgmin = grid_params%jlo(igy)
      jgmax = grid_params%jhi(igy)
      kgmin = grid_params%klo(igz)
      kgmax = grid_params%khi(igz)

c     Determine ASM BCs 

      if (PRESENT(IIBC)) then
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            dn=do__nothing(dim,loc,igmin,igmax,jgmin,jgmax,kgmin,kgmax
     .                    ,asm_bc=.false.)
cc     .         .or.((bcond(ibc)==PER).and.asm_dir(dim))  !Deal with interrupted periodic boundaries
     .         .or.(bcond(ibc)==PER)  !Deal with interrupted periodic boundaries

            if (dn) bcs(ibc,:) = IIBC*bcs(ibc,:)/abs(bcs(ibc,:))
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

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(ieq,nx,ny,nz,array,array0,bcnd
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ieq,bcnd(6),nx,ny,nz,iigx,iigy,iigz
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

c     Local variables

      integer :: dim,loc,ibc
      logical :: petsc__bc

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

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

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,1)  
      iimxp =  min(iimax + 1,grid_params%nxgl(iigx))
      jjmnm =  max(jjmin - 1,1)  
      jjmxp =  min(jjmax + 1,grid_params%nygl(iigy))
      kkmnm =  max(kkmin - 1,1)  
      kkmxp =  min(kkmax + 1,grid_params%nzgl(iigz))
#endif

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc =  (np > 1)  .and. (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => grid_params%dactx(igxbc)
        call fillPetscGhostCells(da_ctx,array)
      endif
#endif

c     Impose BCs

      do bctype=1,BCLIM             !Enforces a particular order in the BCs
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (abs(bcnd(ibc)) == bctype) then
cc              if (bcnd(ibc) < 0 .and. bctype == SYM) then
cc                symm = .true.
cc                call FillGhostNodes(ieq,1,1,dim,loc,DIR,array,array0)
cc                symm = .false.
cc              elseif (bcnd(ibc) > 0 .and. bctype == SYM) then
cc                symm = .true.
cc                call FillGhostNodes(ieq,1,1,dim,loc,NEU,array,array0)
cc                symm = .false.
cc              else
cc                call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
cc              endif

              symm = (bcond(ibc) == SYM)
              call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
              symm = .false.
            endif
          enddo
        enddo
      enddo

c     Fill edges

cc      do dim=1,3
cc        do loc=0,1
cc          call fillEdges(array,array0,dim,loc,bcnd)
cc        enddo
cc      enddo

      call fillEdges(array,array0,bcnd)

c     Fill corners

cc      do dim=1,3
cc        do loc=0,1
cc          call fillCorners(array,array0,dim,loc,bcnd)
cc        enddo
cc      enddo
      call fillCorners(array,array0,bcnd)

#if defined(PER_BC_SYNC)
c     Periodic boundary synchronizing

      bctype=PER
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          if (bcnd(ibc) == bctype) then
            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
          endif
        enddo
      enddo
#endif

c     End program

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(fcomp,ncomp,nx,ny,nz,v_cnv,v_cov,var0
     .                          ,bcnd,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                          ,iorder,is_cnv,is_vec)
c     -----------------------------------------------------------------
c     Imposes BC on vector field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ncomp,bcnd(6,ncomp),fcomp,nx,ny,nz,iigx,iigy,iigz
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,ncomp)
     .             ,var0 (0:nx+1,0:ny+1,0:nz+1,ncomp)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,ncomp)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder
      logical   ,optional,intent(IN) :: is_cnv,is_vec

c     Local variables

      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin
      integer    :: ivar,ieq,ibc,loc,dim,bctype
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: cov_to_cnv(6,3),to_cnv,is__vec,petsc__bc

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

#if defined(samrai)
      iimnm =  iimin 
      iimxp =  iimax
      jjmnm =  jjmin
      jjmxp =  jjmax
      kkmnm =  kkmin
      kkmxp =  kkmax
#else
      iimnm =  max(iimin - 1,1)  
      iimxp =  min(iimax + 1,grid_params%nxgl(iigx))
      jjmnm =  max(jjmin - 1,1)  
      jjmxp =  min(jjmax + 1,grid_params%nygl(iigy))
      kkmnm =  max(kkmin - 1,1)  
      kkmxp =  min(kkmax + 1,grid_params%nzgl(iigz))
#endif

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc = (np > 1).and.(.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => grid_params%dactx(igxbc)
        if (is__vec .and. to_cnv) then
          do ieq=1,ncomp
            call fillPetscGhostCells(da_ctx,v_cov(:,:,:,ieq))
          enddo
        else
          do ieq=1,ncomp
            call fillPetscGhostCells(da_ctx,v_cnv(:,:,:,ieq))
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

      do bctype=1,BCLIM            !Enforces a particular order in the BCs (see grid_def_mod.F)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              ieq = ivar + fcomp - 1

              if (abs(bcnd(ibc,ivar)) == bctype) then

                symm = (bcond(ibc) == SYM)

                if (bcnd(ibc,ivar) < 0 .and. is__vec) then
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cov,var0,cnv=.false.)
                  cov_to_cnv(ibc,ivar) = .true.
                else
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cnv,var0)
                endif

                symm = .false.

             endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components

      call synchronize(ncomp,v_cnv,v_cov,bcnd,cov_to_cnv)

c     Fill edges (for both cnv and cov components)

cc      do dim=1,3
cc        do loc=0,1
cc          do ivar=1,ncomp
cc            call fillEdges(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
cc     .                    ,dim,loc,bcnd(:,ivar))
cc            call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
cc     .                    ,dim,loc,bcnd(:,ivar))
cc          enddo
cc        enddo
cc      enddo

      do ivar=1,ncomp
        call fillEdges(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar))
        call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar))
      enddo

c     Fill corners (for both cnv and cov components)

cc      do dim=1,3
cc        do loc=0,1
cc          do ivar=1,ncomp
cc            call fillCorners(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
cc     .                      ,dim,loc,bcnd(:,ivar))
cc            call fillCorners(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
cc     .                      ,dim,loc,bcnd(:,ivar))
cc          enddo
cc        enddo
cc      enddo

      do ivar=1,ncomp
        call fillCorners(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
     .                  ,bcnd(:,ivar))
        call fillCorners(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                  ,bcnd(:,ivar))
      enddo

#if defined(PER_BC_SYNC)    
c     Periodic boundary synchronizing

      bctype=PER
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,ncomp
            ieq = ivar + fcomp - 1

            if (abs(bcnd(ibc,ivar)) == bctype) then
              call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype,v_cnv
     .                         ,var0)
              call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype,v_cov
     .                         ,var0,cnv=.false.)
            endif
          enddo

        enddo
      enddo
#endif

c     Find covariant components at ALL boundary ghost cells

cc      if (is__vec) then
cc
cc        do dim = 1,3
cc          do loc = 0,1
cc
cc          call findBCLoopLimits(dim,loc,iimin-1,iimax+1
cc     .                                 ,jjmin-1,jjmax+1
cc     .                                 ,kkmin-1,kkmax+1
cc     .                                 ,imin,imax,jmin,jmax,kmin,kmax)
cc
cc          do i=imin,imax
cc            do j=jmin,jmax
cc              do k=kmin,kmax
cc                call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
cc     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
cc     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
cc              enddo
cc            enddo
cc          enddo
cc
cc          enddo
cc        enddo
cc
ccc       Impose BCs on covariant components (if not Cartesian)
cc
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
cc                if (abs(bcnd(ibc,ivar)) == bctype) then
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
cc     .                        ,dim,loc,bcnd(:,ivar))
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
cc     .                          ,dim,loc,bcnd(:,ivar))
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif
cc
cc      endif

c     End program

      end subroutine imposeBConVector

c     synchronize
c     #################################################################
      subroutine synchronize(ncomp,v_cnv,v_cov,bcond,cov_to_cnv)
c     -----------------------------------------------------------------
c     Finds all contravariant components at boundaries with BC in terms
c     of covariant components.
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

c     Local variables

      integer    :: imax,imin,jmax,jmin,kmax,kmin
      integer    :: i,j,k,dim,loc,ig,jg,kg,ivar,ibc,ic1,ic2
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)

c     Begin program

      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)

          call findBCLoopLimits(dim,loc,iimin-1,iimax+1
     .                                 ,jjmin-1,jjmax+1
     .                                 ,kkmin-1,kkmax+1
     .                                 ,imin,imax
     .                                 ,jmin,jmax
     .                                 ,kmin,kmax)

          if (     (.not.cov_to_cnv(ibc,1))
     .        .and.(.not.cov_to_cnv(ibc,2))
     .        .and.(.not.cov_to_cnv(ibc,3))) then  !All contravariant: find covariant

            do i=imin,imax
              do j=jmin,jmax
                do k=kmin,kmax
                  call transformFromCurvToCurv(i,j,k
     .                  ,igxbc,igybc,igzbc
     .                  ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                  ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                  ,.false.)
                enddo
              enddo
            enddo

          else   !Synchronize (assumes normal component is contravariant)

            do ivar = 1,ncomp

              if (ivar == dim) then !Do nothing for normal component

                cycle

              else

                if (.not.selectFace(ibc)) cycle

                ic1 = mod(dim+1,3)
                if (ic1 == 0) ic1 = 3

                ic2 = mod(dim+2,3)
                if (ic2 == 0) ic2 = 3

                if ( cov_to_cnv(ibc,ic1).and.cov_to_cnv(ibc,ic2) ) then  !Both are covariant

                  do i=imin,imax
                    do j=jmin,jmax
                      do k=kmin,kmax

                        gsuper = gmetric%grid(igxbc)%gsup(i,j,k,:,:)

                        v_cov(i,j,k,dim) =
     .                       -(gsuper(dim,ic1)*v_cov(i,j,k,ic1)
     .                        +gsuper(dim,ic2)*v_cov(i,j,k,ic2)
     .                        -v_cnv(i,j,k,dim))/gsuper(dim,dim)

                        call transformFromCurvToCurv(i,j,k
     .                    ,igxbc,igybc,igzbc
     .                    ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                    ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                    ,.true.)

                      enddo
                    enddo
                  enddo

                elseif (.not.cov_to_cnv(ibc,ic1)) then   !ic1 is contravariant, ic2 is covariant

                  do i=imin,imax
                    do j=jmin,jmax
                      do k=kmin,kmax

                        gsub = gmetric%grid(igxbc)%gsub(i,j,k,:,:)

                        v_cnv(i,j,k,ic2) =
     .                       -(gsub(ic2,ic1)*v_cnv(i,j,k,ic1)
     .                        +gsub(ic2,dim)*v_cnv(i,j,k,dim)
     .                        -v_cov(i,j,k,ic2))/gsub(ic2,ic2)

                        call transformFromCurvToCurv(i,j,k
     .                    ,igxbc,igybc,igzbc
     .                    ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                    ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                    ,.false.)

                      enddo
                    enddo
                  enddo

                elseif (.not.cov_to_cnv(ibc,ic2)) then   !ic2 is contravariant, ic1 is covariant

                  do i=imin,imax
                    do j=jmin,jmax
                      do k=kmin,kmax

                        gsub = gmetric%grid(igxbc)%gsub(i,j,k,:,:)

                        v_cnv(i,j,k,ic1) =
     .                       -(gsub(ic1,ic2)*v_cnv(i,j,k,ic2)
     .                        +gsub(ic1,dim)*v_cnv(i,j,k,dim)
     .                        -v_cov(i,j,k,ic1))/gsub(ic1,ic1)

                        call transformFromCurvToCurv(i,j,k
     .                    ,igxbc,igybc,igzbc
     .                    ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                    ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                    ,.false.)

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
     .                      ,dim,loc,order)
        case(DIR)
          call scalarDirichletBC(array(:,:,:,nvar),zeros,ieq,dim,loc
     .                          ,order)
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
     .             ,'BC '//trim(int2char(bctype))//' not implemented')
        end select
      else
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,ivar),dim,loc)
        case(EQU)
          call vectorDirichletBC(ivar,array,array0,ieq,dim,loc,order)
        case(DIR)
          call vectorDirichletBC(ivar,array,vzeros,ieq,dim,loc,order)
        case(NEU)
          call vectorNeumannBC(ivar,array,array0,ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,ivar),dim,loc)
        case(SP)
          call vectorSingularBC(ivar,array(:,:,:,ivar),cov)
        case(IFC)
          call interpolateToFace(array(:,:,:,ivar),dim,loc)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC '//trim(int2char(bctype))//' not implemented')
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

#if !defined(samrai)
      ibc = (1+loc)+2*(dim-1)

      if (np == 1) then
cc      if (np == 1 .or. asm) then  !Do periodic BCs for ASM when possible

#if !defined(PER_BC_SYNC)

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

#endif

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

        integer    :: i,j,k,ig,jg,kg,one=1,two=2
        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

c     Begin program

        ibc = (1+loc)+2*(dim-1)

        rhs = 0d0

        if (order == 1) two = one  !To prevent stepping over domain

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

     .              quad_int(0.5*(xx(ig-one)+xx(ig)),xx(ig)
     .                       ,xx(ig+one),xx(ig+two)
     .                       ,array0(i-one,j,k),array(i    ,j,k)
     .                       ,array (i+one,j,k),array(i+two,j,k)
     .                       ,xx(ig-one),order )
                case (2)

                  rhs(j,k) =

     .              quad_int(0.5*(xx(ig+one)+xx(ig)),xx(ig)
     .                        ,xx(ig-one),xx(ig-two)
     .                        ,array0(i+one,j,k),array(i,j    ,k)
     .                        ,array (i-one,j,k),array(i-two,j,k)
     .                        ,xx(ig+one),order )
                case (3)

                  rhs(i,k) =

     .              quad_int(0.5*(yy(jg-one)+yy(jg)),yy(jg)
     .                        ,yy(jg+one),yy(jg+two)
     .                        ,array0(i,j-one,k),array(i,j    ,k)
     .                        ,array (i,j+one,k),array(i,j+two,k)
     .                        ,yy(jg-one),order )
                case (4)

                  rhs(i,k) =

     .              quad_int(0.5*(yy(jg+one)+yy(jg)),yy(jg)
     .                        ,yy(jg-one),yy(jg-two)
     .                        ,array0(i,j+one,k),array(i,j    ,k)
     .                        ,array (i,j-one,k),array(i,j-two,k)
     .                        ,yy(jg+one),order )
                case (5)

                  rhs(i,j) =

     .              quad_int(0.5*(zz(kg-one)+zz(kg)),zz(kg)
     .                        ,zz(kg+one),zz(kg+one)
     .                        ,array0(i,j,k-one),array(i,j,k    )
     .                        ,array (i,j,k+one),array(i,j,k+two)
     .                        ,zz(kg-one),order )
                case (6)

                  rhs(i,j) =

     .              quad_int(0.5*(zz(kg+one)+zz(kg)),zz(kg)
     .                        ,zz(kg-one),zz(kg-two)
     .                        ,array0(i,j,k+one),array(i,j,k    )
     .                        ,array (i,j,k-one),array(i,j,k-two)
     .                        ,zz(kg+one),order )
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

        do k=kmin,kmax
          do j=jmin,jmax
            do i=imin,imax

              call getMGmap(i,j,k,igxbc,igybc,igzbc,ig,jg,kg)

              select case (ibc)
              case (1)
                rhs(j,k) = 
     .            quad_int(xx(ig-1),xx(ig),xx(ig+1),xx(ig+2)
     .                      ,array(i-1,j,k),array(i,j,k)
     .                      ,array(i+1,j,k),array(i+2,j,k)
     .                      ,xx(ig-1)+dxh(ig-1),order )
              case (2)
                rhs(j,k) =
     .            quad_int(xx(ig+1),xx(ig),xx(ig-1),xx(ig-2)
     .                      ,array(i+1,j,k),array(i,j,k)
     .                      ,array(i-1,j,k),array(i-2,j,k)
     .                      ,xx(ig+1)-dxh(ig+1),order )
              case (3)

                rhs(i,k) =
     .            quad_int(yy(jg-1),yy(jg),yy(jg+1),yy(jg+2)
     .                      ,array(i,j-1,k),array(i,j,k)
     .                      ,array(i,j+1,k),array(i,j+2,k)
     .                      ,yy(jg-1)+dyh(jg-1),order )
              case (4)

                rhs(i,k) =
     .            quad_int(yy(jg+1),yy(jg),yy(jg-1),yy(jg-2)
     .                      ,array(i,j+1,k),array(i,j,k)
     .                      ,array(i,j-1,k),array(i,j-2,k)
     .                      ,yy(jg+1)-dyh(jg+1),order )
              case (5)

                rhs(i,j) =
     .            quad_int(zz(kg-1),zz(kg),zz(kg+1),zz(kg+2)
     .                      ,array(i,j,k-1),array(i,j,k)
     .                      ,array(i,j,k+1),array(i,j,k+2)
     .                      ,zz(kg-1)+dzh(kg-1),order )
              case (6)

                rhs(i,j) =
     .            quad_int(zz(kg+1),zz(kg),zz(kg-1),zz(kg-2)
     .                      ,array(i,j,k+1),array(i,j,k)
     .                      ,array(i,j,k-1),array(i,j,k-2)
     .                      ,zz(kg+1)-dzh(kg+1),order )
              end select
            enddo
          enddo
        enddo

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

      integer    :: ibcl(3)
      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,bcond(6)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      bcond = 0
      bcond(ibc) = EXT

      ibcl = 0

c     Determine boundary limits

      if (.not.selectFace(ibc)) return

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (bc_debug) then
        write(*,111) '   fillFaces: Proc',my_rank
     .            ,'; Boundary: dim =',dim,' loc =',loc
     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
 111    format (a,i4,a,i4,a,i4,a,6i4)
      endif

c     Select face

      select case (ibc)
      case (1) !On face #1, fill face i=0

        if (imin == 1) then

          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          ibcl(1) = 1

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibcl,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j:j,k:k),bcond)
            enddo
          enddo

        endif

      case (2) !On face #2, fill face i=nx+1

        if (imax == nxbc) then

          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          ibcl(1) = 2

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibcl,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j:j,k:k),bcond)
            enddo
          enddo

        endif

      case (3) !On face #3, fill face j=0

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          ibcl(2) = 3

          do k=kmin,kmax
            do i=imin,imax
                call extrapolateToCorner(ibcl,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i:i,jmn:jmx,k:k),bcond)
            enddo
          enddo

        endif

      case (4) !On face #4, fill face j=ny+1

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          ibcl(2) = 4

          do k=kmin,kmax
            do i=imin,imax
              call extrapolateToCorner(ibcl,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i:i,jmn:jmx,k:k),bcond)
            enddo
          enddo

        endif

      case (5) !On face #5, fill face k=0

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          ibcl(3) = 5

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibcl,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i:i,j:j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      case (6) !On face #6, fill face k=nz+1

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          ibcl(3) = 6

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibcl,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i:i,j:j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      end select

c     End program

      end subroutine fillFaces

ccc     fillEdges
ccc     #################################################################
cc      subroutine fillEdges(array,array0,dim,loc,bcond)
ccc     -----------------------------------------------------------------
ccc     Fills ghost cells in edges (12 in 3D). On input:
ccc        * array  -> 3D array to impose BC's on.
ccc        * array0 -> 3D array containing equilibrium BCs
ccc        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
ccc        * loc    -> location in dimension (0 -> right, 1 -> left)
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: dim,loc,bcond(6)
cc      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
cc     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)
cc
ccc     Local variables
cc
cc      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin
cc
cc      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
cc     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,ord
cc
ccc     Begin program
cc
cc      ibc = (1+loc)+2*(dim-1)
cc
ccc     Determine boundary limits
cc
cccc      if (.not.selectFace(ibc)) return
cc
cc      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
cc     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)
cc
cc      ord = max(order,1)  !Does not affect order of interpolation in extrapolateToCorner
cc
cc      if (bc_debug) then
cc        write(*,111) '   fillEdges: Proc',my_rank
cc     .            ,'; Boundary: dim =',dim,' loc =',loc
cc     .            ,'; Local limits:',imin,imax,jmin,jmax,kmin,kmax
cc     .            ,'; Order =',order
cc 111    format (a,i4,a,i4,a,i4,a,6i4,a,i4)
cc      endif
cc
ccc     Select edge
cc
cc      select case (ibc)
cc      case (1) !On face #1, fill edges y=jmin-1,y=jmax+1, vary z
cc
cc        if (imin == 1) then
cc          ito   = imin-1
cc          ifrom = min(imin+ord,nxbc)
cc
cc          imn   = min(ito,ifrom)
cc          imx   = max(ito,ifrom)
cc
cc          if (jmin == 1.and.(selectFace(ibc).or.selectFace(3))) then  !Check faces #1,#3
cc            jto   = jmin-1
cc            jfrom = min(jmin+ord,nybc)
cc
cc            jmn   = min(jto,jfrom)
cc            jmx   = max(jto,jfrom)
cc
cc            do k=kmin,kmax
cc              if ((bcond(1)==EQU.or.bcond(3)==EQU) .and.order == 0) then
cc                array(ito,jto,k) = array0(ito,jto,k)
cc              else
cccc                write (*,*)
cccc                write (*,*) order,ito,ifrom,jto,jfrom,kto,kfrom
cccc                write (*,*) nxbc,nybc
cc                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
cc     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc          if (jmax == nybc.and.(selectFace(ibc).or.selectFace(4))) then  !Check faces #1,#4
cc            jto   = jmax+1
cc            jfrom = max(jmax-ord,1)
cc
cc            jmn   = min(jto,jfrom)
cc            jmx   = max(jto,jfrom)
cc
cc            do k=kmin,kmax
cc              if ((bcond(1)==EQU.or.bcond(4)==EQU).and.order == 0) then
cc                array(ito,jto,k) = array0(ito,jto,k)
cc              else
cc                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
cc     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc        endif
cc
cc      case (2) !On face #2, fill edges y=jmin-1,y=jmax+1, vary z
cc
cc        if (imax == nxbc) then
cc          ito   = imax+1
cc          ifrom = max(imax-ord,1)
cc
cc          imn   = min(ito,ifrom)
cc          imx   = max(ito,ifrom)
cc
cc          if (jmin == 1.and.(selectFace(ibc).or.selectFace(3))) then  !Check faces #2,#3
cc            jto   = jmin-1
cc            jfrom = min(jmin+ord,nybc)
cc
cc            jmn   = min(jto,jfrom)
cc            jmx   = max(jto,jfrom)
cc
cc            do k=kmin,kmax
cc              if ((bcond(2)==EQU.or.bcond(3)==EQU).and.order == 0) then
cc                array(ito,jto,k) = array0(ito,jto,k)
cc              else
cc                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
cc     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc          if (jmax == nybc.and.(selectFace(ibc).or.selectFace(4))) then  !Check faces #2,#4
cc            jto   = jmax+1
cc            jfrom = max(jmax-ord,1)
cc
cc            jmn   = min(jto,jfrom)
cc            jmx   = max(jto,jfrom)
cc
cc            do k=kmin,kmax
cc              if ((bcond(2)==EQU.or.bcond(4)==EQU).and.order == 0) then
cc                array(ito,jto,k) = array0(ito,jto,k)
cc              else
cc                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
cc     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc        endif
cc
cc      case (3) !On face #3, fill edges z=kmin-1,z=kmax+1, vary x
cc
cc        if (jmin == 1) then
cc
cc          jto   = jmin-1
cc          jfrom = min(jmin+ord,nybc)
cc
cc          jmn   = min(jto,jfrom)
cc          jmx   = max(jto,jfrom)
cc
cc          if (kmin == 1.and.(selectFace(ibc).or.selectFace(5))) then  !Check faces #3,#5
cc
cc            kto   = kmin-1
cc            kfrom = min(kmin+ord,nzbc)
cc
cc            kmn   = min(kto,kfrom)
cc            kmx   = max(kto,kfrom)
cc
cc            do i=imin,imax
cc              if ((bcond(3)==EQU.or.bcond(5)==EQU).and.order == 0) then
cc                array(i,jto,kto) = array0(i,jto,kto)
cc              else
cc                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
cc     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc          if (kmax == nzbc.and.(selectFace(ibc).or.selectFace(6))) then  !Check faces #3,#6
cc
cc            kto   = kmax+1
cc            kfrom = max(kmax-ord,1)
cc
cc            kmn   = min(kto,kfrom)
cc            kmx   = max(kto,kfrom)
cc
cc            do i=imin,imax
cc              if ((bcond(3)==EQU.or.bcond(6)==EQU).and.order == 0) then
cc                array(i,jto,kto) = array0(i,jto,kto)
cc              else
cc                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
cc     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc        endif
cc
cc      case (4) !On face #4, fill edges z=kmin-1,z=kmax+1, vary x
cc
cc        if (jmax == nybc) then
cc
cc          jto   = jmax+1
cc          jfrom = max(jmax-ord,1)
cc
cc          jmn   = min(jto,jfrom)
cc          jmx   = max(jto,jfrom)
cc
cc          if (kmin == 1.and.(selectFace(ibc).or.selectFace(5))) then  !Check faces #4,#5
cc
cc            kto   = kmin-1
cc            kfrom = min(kmin+ord,nzbc)
cc
cc            kmn   = min(kto,kfrom)
cc            kmx   = max(kto,kfrom)
cc
cc            do i=imin,imax
cc              if ((bcond(4)==EQU.or.bcond(5)==EQU).and.order == 0) then
cc                array(i,jto,kto) = array0(i,jto,kto)
cc              else
cc                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
cc     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc          if (kmax == nzbc.and.(selectFace(ibc).or.selectFace(6))) then  !Check faces #4,#6
cc
cc            kto   = kmax+1
cc            kfrom = max(kmax-ord,1)
cc
cc            kmn   = min(kto,kfrom)
cc            kmx   = max(kto,kfrom)
cc
cc            do i=imin,imax
cc              if ((bcond(4)==EQU.or.bcond(6)==EQU).and.order == 0) then
cc                array(i,jto,kto) = array0(i,jto,kto)
cc              else
cc                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
cc     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc        endif
cc
cc      case (5) !On face #5, fill edges x=imin-1,x=imax+1, vary y
cc
cc        if (kmin == 1) then
cc
cc          kto   = kmin-1
cc          kfrom = min(kmin+ord,nzbc)
cc
cc          kmn   = min(kto,kfrom)
cc          kmx   = max(kto,kfrom)
cc
cc          if (imin == 1.and.(selectFace(ibc).or.selectFace(1))) then  !Check faces #5,#1
cc
cc            ito   = imin-1
cc            ifrom = min(imin+ord,nxbc)
cc
cc            imn   = min(ito,ifrom)
cc            imx   = max(ito,ifrom)
cc
cc            do j=jmin,jmax
cc              if ((bcond(5)==EQU.or.bcond(1)==EQU).and.order == 0) then
cc                array(ito,j,kto) = array0(ito,j,kto)
cc              else
cc                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
cc     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc          if (imax == nxbc.and.(selectFace(ibc).or.selectFace(2))) then  !Check faces #5,#2
cc
cc            ito   = imax+1
cc            ifrom = max(imax-ord,1)
cc
cc            imn   = min(ito,ifrom)
cc            imx   = max(ito,ifrom)
cc
cc            do j=jmin,jmax
cc              if ((bcond(5)==EQU.or.bcond(2)==EQU).and.order == 0) then
cc                array(ito,j,kto) = array0(ito,j,kto)
cc              else
cc                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
cc     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc        endif
cc
cc      case (6) !On face #6, fill edges x=imin-1,x=imax+1, vary y
cc
cc        if (kmax == nzbc) then
cc
cc          kto   = kmax+1
cc          kfrom = max(kmax-ord,1)
cc
cc          kmn   = min(kto,kfrom)
cc          kmx   = max(kto,kfrom)
cc
cc          if (imin == 1.and.(selectFace(ibc).or.selectFace(1))) then  !Check faces #6,#1
cc
cc            ito   = imin-1
cc            ifrom = min(imin+ord,nxbc)
cc
cc            imn   = min(ito,ifrom)
cc            imx   = max(ito,ifrom)
cc
cc            do j=jmin,jmax
cc              if ((bcond(6)==EQU.or.bcond(1)==EQU).and.order == 0) then
cc                array(ito,j,kto) = array0(ito,j,kto)
cc              else
cc                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
cc     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc          if (imax == nxbc.and.(selectFace(ibc).or.selectFace(2))) then  !Check faces #6,#2
cc
cc            ito   = imax+1
cc            ifrom = max(imax-ord,1)
cc
cc            imn   = min(ito,ifrom)
cc            imx   = max(ito,ifrom)
cc
cc            do j=jmin,jmax
cc              if ((bcond(6)==EQU.or.bcond(2)==EQU).and.order == 0) then
cc                array(ito,j,kto) = array0(ito,j,kto)
cc              else
cc                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
cc     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
cc              endif
cc            enddo
cc          endif
cc
cc        endif
cc
cc      end select
cc
ccc     End program
cc
cc      end subroutine fillEdges

c     fillEdges
c     #################################################################
      subroutine fillEdges(array,array0,bcond)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin

      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,ord
     .             ,dim,loc,ibcl(3)

c     Begin program

      do dim=1,3
      do loc=0,1

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      ord = max(order,1)  !Does not affect order of interpolation in extrapolateToCorner

      if (bc_debug) then
        write(*,111) '   fillEdges: Proc',my_rank
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

          if (jmin == 1.and.(selectFace(1).or.selectFace(3))) then  !Check faces #1,#3
            jto   = jmin-1
            jfrom = min(jmin+ord,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(3)==EQU) .and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/1,3,0/)
                call extrapolateToCorner(ibcl,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc.and.(selectFace(1).or.selectFace(4))) then  !Check faces #1,#4
            jto   = jmax+1
            jfrom = max(jmax-ord,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/1,4,0/)
                call extrapolateToCorner(ibcl,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond)
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

          if (jmin == 1.and.(selectFace(2).or.selectFace(3))) then  !Check faces #2,#3
            jto   = jmin-1
            jfrom = min(jmin+ord,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(3)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/2,3,0/)
                call extrapolateToCorner(ibcl,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc.and.(selectFace(2).or.selectFace(4))) then  !Check faces #2,#4
            jto   = jmax+1
            jfrom = max(jmax-ord,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(4)==EQU).and.order == 0) then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                ibcl = (/2,4,0/)
                call extrapolateToCorner(ibcl,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond)
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

          if (kmin == 1.and.(selectFace(3).or.selectFace(5))) then  !Check faces #3,#5

            kto   = kmin-1
            kfrom = min(kmin+ord,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,3,5/)
                call extrapolateToCorner(ibcl,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc.and.(selectFace(3).or.selectFace(6))) then  !Check faces #3,#6

            kto   = kmax+1
            kfrom = max(kmax-ord,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,3,6/)
                call extrapolateToCorner(ibcl,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond)
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

          if (kmin == 1.and.(selectFace(4).or.selectFace(5))) then  !Check faces #4,#5

            kto   = kmin-1
            kfrom = min(kmin+ord,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(5)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,4,5/)
                call extrapolateToCorner(ibcl,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc.and.(selectFace(4).or.selectFace(6))) then  !Check faces #4,#6

            kto   = kmax+1
            kfrom = max(kmax-ord,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(6)==EQU).and.order == 0) then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                ibcl = (/0,4,6/)
                call extrapolateToCorner(ibcl,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond)
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

          if (imin == 1.and.(selectFace(5).or.selectFace(1))) then  !Check faces #5,#1

            ito   = imin-1
            ifrom = min(imin+ord,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/1,0,5/)
                call extrapolateToCorner(ibcl,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc.and.(selectFace(5).or.selectFace(2))) then  !Check faces #5,#2

            ito   = imax+1
            ifrom = max(imax-ord,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/2,0,5/)
                call extrapolateToCorner(ibcl,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond)
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

          if (imin == 1.and.(selectFace(6).or.selectFace(1))) then  !Check faces #6,#1

            ito   = imin-1
            ifrom = min(imin+ord,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(1)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/1,0,6/)
                call extrapolateToCorner(ibcl,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc.and.(selectFace(6).or.selectFace(2))) then  !Check faces #6,#2

            ito   = imax+1
            ifrom = max(imax-ord,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(2)==EQU).and.order == 0) then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                ibcl = (/2,0,6/)
                call extrapolateToCorner(ibcl,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      end select

      enddo
      enddo

c     End program

      end subroutine fillEdges

ccc     fillCorners
ccc     #################################################################
cc      subroutine fillCorners(array,array0,dim,loc,bcond)
ccc     -----------------------------------------------------------------
ccc     Imposes ghost cells in corners (8 in 3D). On input:
ccc        * array -> 3D array to impose BC's on.
ccc        * array0 -> 3D array containing equilibrium BCs
ccc        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
ccc        * loc    -> location in dimension (0 -> right, 1 -> left)
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: dim,loc,bcond(6)
cc      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
cc     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)
cc
ccc     Local variables
cc
cc      integer    :: imax,imin,jmax,jmin,kmax,kmin
cc
cc      integer    :: ito,jto,kto,ifrom,jfrom,kfrom
cc     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,ord
cc
ccc     Begin program
cc
cc      ibc = (1+loc)+2*(dim-1)
cc
cc      ito= UNDEF_INT
cc      jto= UNDEF_INT
cc      kto= UNDEF_INT
cc
ccc     Determine boundary limits
cc
cccc      if (.not.selectFace(ibc)) return
cc
cc      call findBCLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
cc     .                             ,imin ,imax ,jmin ,jmax ,kmin ,kmax)
cc
cc      ord = max(order,1)   !Does not affect order of interpolation in extrapolateToCorner
cc
ccc     Select corner
cc
cc      if (      imin == 1  .and. jmin == 1  .and. kmin == 1
cc     .    .and.(selectFace(1).or.selectFace(3).or.selectFace(5)) ) then   !Check concurrent faces
cc        ifrom = min(imin+ord,nxbc)
cc        jfrom = min(jmin+ord,nybc)
cc        kfrom = min(kmin+ord,nzbc)
cc        ito   = imin-1
cc        jto   = jmin-1
cc        kto   = kmin-1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc         endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .          ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imin == 1  .and. jmax == nybc .and. kmin == 1
cc     .    .and.(selectFace(1).or.selectFace(4).or.selectFace(5))) then
cc
cc        ifrom = min(imin+ord,nxbc)
cc        jfrom = max(jmax-ord,1   )
cc        kfrom = min(kmin+ord,nzbc)
cc        ito   = imin-1
cc        jto   = jmax+1
cc        kto   = kmin-1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .           ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imax == nxbc .and. jmin == 1  .and. kmin == 1
cc     .    .and.(selectFace(2).or.selectFace(3).or.selectFace(5))) then
cc
cc        ifrom = max(imax-ord,1   )
cc        jfrom = min(jmin+ord,nybc)
cc        kfrom = min(kmin+ord,nzbc)
cc        ito   = imax+1
cc        jto   = jmin-1
cc        kto   = kmin-1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imax == nxbc .and. jmax == nybc .and. kmin == 1
cc     .    .and.(selectFace(2).or.selectFace(4).or.selectFace(5))) then
cc
cc        ifrom = max(imax-ord,1   )
cc        jfrom = max(jmax-ord,1   )
cc        kfrom = min(kmin+ord,nzbc)
cc        ito   = imax+1
cc        jto   = jmax+1
cc        kto   = kmin-1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imin == 1  .and. jmin == 1  .and. kmax == nzbc
cc     .    .and.(selectFace(1).or.selectFace(3).or.selectFace(6))) then
cc
cc        ifrom = min(imin+ord,nxbc)
cc        jfrom = min(jmin+ord,nybc)
cc        kfrom = max(kmax-ord,1   )
cc        ito   = imin-1
cc        jto   = jmin-1
cc        kto   = kmax+1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imin == 1  .and. jmax == nybc .and. kmax == nzbc
cc     .    .and.(selectFace(1).or.selectFace(4).or.selectFace(6))) then
cc
cc        ifrom = min(imin+ord,nxbc)
cc        jfrom = max(jmax-ord,1   )
cc        kfrom = max(kmax-ord,1   )
cc        ito   = imin-1
cc        jto   = jmax+1
cc        kto   = kmax+1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imax == nxbc .and. jmin == 1  .and. kmax == nzbc
cc     .    .and.(selectFace(2).or.selectFace(3).or.selectFace(6))) then
cc
cc        ifrom = max(imax-ord,1   )
cc        jfrom = min(jmin+ord,nybc)
cc        kfrom = max(kmax-ord,1   )
cc        ito   = imax+1
cc        jto   = jmin-1
cc        kto   = kmax+1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc      endif
cc
cc      if (     imax == nxbc .and. jmax == nybc .and. kmax == nzbc
cc     .    .and.(selectFace(2).or.selectFace(4).or.selectFace(6))) then
cc
cc        ifrom = max(imax-ord,1 )
cc        jfrom = max(jmax-ord,1 )
cc        kfrom = max(kmax-ord,1 )
cc        ito   = imax+1
cc        jto   = jmax+1
cc        kto   = kmax+1
cc
cc        imn   = min(ito,ifrom)
cc        jmn   = min(jto,jfrom)
cc        kmn   = min(kto,kfrom)
cc        imx   = max(ito,ifrom)
cc        jmx   = max(jto,jfrom)
cc        kmx   = max(kto,kfrom)
cc
cc        if ( bc_debug
cc     .      .and. (ito /= UNDEF_INT
cc     .        .and.jto /= UNDEF_INT
cc     .        .and.kto /= UNDEF_INT) ) then
cc          write (*,111) '   fillCorners: Proc',my_rank
cc     .                 ,'; Corner =',ito,jto,kto
cc     .                 ,'; Order =',order
cc        endif
cc
cc        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
cc     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
cc
cc      endif
cc
cc      if (bc_debug) then
cc        write (*,*)
cc#if defined(petsc)
cc        call MPI_barrier(MPI_COMM_WORLD,mpierr)
cc#endif
cccc        stop
cc      endif
cc
ccc     End program
cc
cc 111  format (a,i4,a,3i4,a,i4)
cc
cc      end subroutine fillCorners

c     fillCorners
c     #################################################################
      subroutine fillCorners(array,array0,bcond)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

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
     .          ,imn,imx,jmn,jmx,kmn,kmx,ord,ibcl(3)

      logical :: sel_f(3)

c     Begin program

      ord = max(order,1)   !Does not affect order of interpolation in extrapolateToCorner

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

      if ((selectFace(f1).or.selectFace(f2).or.selectFace(f3))) then

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)

        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        if ( bc_debug ) then
          write (*,111) '   fillCorners: Proc',my_rank
     .                 ,'; Corner =',ito,jto,kto
     .                 ,'; From =',ifrom,jfrom,kfrom
     .                 ,'; Order =',order

 111      format (a,i4,a,3i4,a,3i4,a,i4)
        endif

        ibcl  = (/ f1, f2, f3 /)

        call extrapolateToCorner(ibcl,ito,jto,kto
     .                   ,imn,jmn,kmn,imx,jmx,kmx
     .             ,array(imn:imx,jmn:jmx,kmn:kmx),bcond)

      endif

      end subroutine fillCorner

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(ibcl,ito,jto,kto,imn,jmn,kmn
     .                              ,imx,jmx,kmx,array,bcond,array0)
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
c     The order of extrapolation 'order' is defined in the module header.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ibcl(3),ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx,bcond(6)
      real(8)    :: array(imn:imx,jmn:jmx,kmn:kmx)
      real(8),optional,intent(IN) :: array0(imn:imx,jmn:jmx,kmn:kmx)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,inx,iny,inz
     .             ,order_int(3),offset
     .             ,imni,imxi,jmni,jmxi,kmni,kmxi,num
     .             ,i_bc,i1,i2

      real(8)    :: valx(1),valy(1),valz(1),xto(1),yto(1),zto(1)
     .             ,smooth(3)

      logical    :: flg(3),select_bc(3),have_equ
     .             ,face,edge,corner

c     Extrapolation

      real(8),allocatable,dimension(:) :: xx,yy,zz

c     Begin program

c     Do nothing for SP boundary face, edge, and corners

      if (     isSP2(ito+1,1,ibc=ibcl(1))
     .    .or. isSP2(ito+1,1,ibc=ibcl(2))
     .    .or. isSP2(ito+1,1,ibc=ibcl(3))) return

c     Initialize variables

      have_equ = PRESENT(array0)

      face   = (count(ibcl /= 0) == 1)
      edge   = (count(ibcl /= 0) == 2)
      corner = (count(ibcl /= 0) == 3)

c     Identify directions for extrapolation.

      inx = imx-imn
      iny = jmx-jmn
      inz = kmx-kmn

      !Order in each direction
      !Directions 1,2,3 are defined as normal to faces.
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

      !Do not interpolate **along** direction if not enough points,
      !Directions 1,2,3 are defined as normal to faces.
      flg = .true.

      if (inx == 0) flg(1) = .false.
      if (iny == 0) flg(2) = .false.
      if (inz == 0) flg(3) = .false.

      !Select boundaries involved for extrapolation
      ![ibc = (1+loc)+2*(dim-1), with dim=1,2,3, and loc=0(left),1(right)]
      !Directions 1,2,3 are defined as normal to faces.

cc      ibcl = 0
cc
cc      if (flg(1) .and. flg(2) .and. flg(3)) then !Corner
cc
cc        corner = .true.
cc
cc        if (ito == 0 .and. jto == 0 .and. kto == 0) then
cc          ibcl(1) = 1
cc          ibcl(2) = 3
cc          ibcl(3) = 5
cc        elseif (ito == 0 .and. jto == 0 .and. kto /= 0) then
cc          ibcl(1) = 1
cc          ibcl(2) = 3
cc          ibcl(3) = 6
cc        elseif (ito == 0 .and. jto /= 0 .and. kto == 0) then
cc          ibcl(1) = 1
cc          ibcl(2) = 4
cc          ibcl(3) = 5
cc        elseif (ito /= 0 .and. jto == 0 .and. kto == 0) then
cc          ibcl(1) = 2
cc          ibcl(2) = 3
cc          ibcl(3) = 5
cc        elseif (ito /= 0 .and. jto /= 0 .and. kto == 0) then
cc          ibcl(1) = 2
cc          ibcl(2) = 4
cc          ibcl(3) = 5
cc        elseif (ito /= 0 .and. jto == 0 .and. kto /= 0) then
cc          ibcl(1) = 2
cc          ibcl(2) = 3
cc          ibcl(3) = 6
cc        elseif (ito == 0 .and. jto /= 0 .and. kto /= 0) then
cc          ibcl(1) = 1
cc          ibcl(2) = 4
cc          ibcl(3) = 6
cc        elseif (ito /= 0 .and. jto /= 0 .and. kto /= 0) then
cc          ibcl(1) = 2
cc          ibcl(2) = 4
cc          ibcl(3) = 6
cc        endif
cc
cc      elseif (   (flg(1) .and. flg(2))
cc     $       .or.(flg(1) .and. flg(3))
cc     $       .or.(flg(2) .and. flg(3))) then !Edge
cc
cc        edge = .true.
cc
cc        select case(ibc)
cc        case (1,2)
cc          ibcl(1) = ibc
cc          if (jto == 0) then
cc            ibcl(2) = 3
cc          else
cc            ibcl(2) = 4
cc          endif
cc        case (3,4)
cc          ibcl(2) = ibc
cc          if (kto == 0) then
cc            ibcl(3) = 5
cc          else
cc            ibcl(3) = 6
cc          endif
cc        case (5,6)
cc          ibcl(3) = ibc
cc          if (ito == 0) then
cc            ibcl(1) = 1
cc          else
cc            ibcl(1) = 2
cc          endif
cc        end select
cc
cc      else !Face
cc
cc        face = .true.
cc
cc        select case(ibc)
cc        case (1,2)
cc          ibcl(1) = ibc
cc        case (3,4)
cc          ibcl(2) = ibc
cc        case (5,6)
cc          ibcl(3) = ibc
cc        end select
cc
cc      endif

c     Treat special boundaries

      !Identify special boundaries to interpolate **along**
      !   1) SP boundary
      !   2) Interpolate-to-face boundary
      !   3) Equilibrium boundary (Removed 12/13/06, added back 11/10/08)
      !   4) Discard internal processor boundaries (not real boundaries)
      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc <= 0) cycle
        select_bc(i) =(    (abs(bcond(i_bc)) == SP )
     .                 .or.(abs(bcond(i_bc)) == IFC)
     .                 .or.(abs(bcond(i_bc)) == EQU)
     .                ).and.selectFace(i_bc)
      enddo

      call select_along

      !Identify special boundaries to interpolate **across**
      !   1) Periodic boundary
      !   2) Symmetry boundary (either SYM or NEU+order=0)
      !   3) Equilibrium boundary (Removed 11/10/08; was interfering with parallelization)
      !   4) Discard internal processor boundaries (not real boundaries)
      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc <= 0) cycle
        select_bc(i) =(    (    bcond(i_bc)  == PER)
     .                .or. (abs(bcond(i_bc)) == SYM)
     .                .or. (abs(bcond(i_bc)) == NEU .and. order == 0)
     .                ).and.selectFace(i_bc)
      enddo

cc      if (edge) write (*,*)  'DIAG -- extrapolateToCorner, select_bc'
cc     .                       ,select_bc
      call select_across

c     Check some direction for interpolation exists

      if ((.not.flg(1)).and.(.not.flg(2)).and.(.not.flg(3))) then
        write (*,*) 'Extrapolation directions',flg
        write (*,*) 'Boundary type: Face=',face,' Edge=',edge
     .             ,' Corner=',corner
        write (*,*) 'Boundaries involved',ibcl
        call pstop('extrapolateToCorner'
     .            ,'No extrapolation directions left')
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

        if (isSYM(ito+1,igxbc,ibc=ibcl(1))) then
          valx(1) = sign(1,bcond(ibcl(1)))*array(ito+2,jto,kto)
        elseif (isSYM(ito-1,igxbc,ibc=ibcl(1))) then
          valx(1) = sign(1,bcond(ibcl(1)))*array(ito-2,jto,kto)
        else
          call IntDriver1d(inx,xx,array(imni:imxi,jto,kto)
     .                    ,1,xto,valx,order_int(1))
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

        if (isSYM(jto+1,igybc,ibc=ibcl(2))) then
          valy(1) = sign(1,bcond(ibcl(2)))*array(ito,jto+2,kto)
        elseif (isSYM(jto-1,igybc,ibc=ibcl(2))) then
          valy(1) = sign(1,bcond(ibcl(2)))*array(ito,jto-2,kto)
        else
          call IntDriver1d(iny,yy,array(ito,jmni:jmxi,kto)
     .                    ,1,yto,valy,order_int(2))
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

        if (isSYM(kto+1,igzbc,ibc=ibcl(3))) then
          valz(1) = sign(1,bcond(ibcl(3)))*array(ito,jto,kto+2)
        elseif (isSYM(kto-1,igzbc,ibc=ibcl(3))) then
          valz(1) = sign(1,bcond(ibcl(3)))*array(ito,jto,kto-2)
        else
          call IntDriver1d(inz,zz,array(ito,jto,kmni:kmxi)
     .                    ,1,zto,valz,order_int(3))
        endif

        num = num + 1

        deallocate(zz)
      endif

c     Average extrapolation results and fill ghost cell

      array(ito,jto,kto) = (valx(1) + valy(1) + valz(1))/num

      if (bc_debug) then
       if(face) write(*,222)'DIAG -- extrapolateToCorner - face, Proc='
     .                      ,my_rank,num,flg,ito,jto,kto,valx,valy,valz
       if(edge) write(*,222)'DIAG -- extrapolateToCorner - edge, Proc='
     .                      ,my_rank,num,flg,ito,jto,kto,valx,valy,valz
       if(corner)write(*,222)'DIAG -- extrapolateToCorner -corner,Proc='
     .                      ,my_rank,num,flg,ito,jto,kto,valx,valy,valz
 222   format(a,2i3,3l2,3i3,3e14.5)
      endif

c     End program

      contains

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

      sel_f = (/selectFace(ibcl(1))
     .         ,selectFace(ibcl(2))
     .         ,selectFace(ibcl(3))/)

      where (.not.sel_f) ibcl = -ibcl

      end subroutine checkFaces

c     select_along
c     ######################################################################
      subroutine select_along

      implicit none

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

        !Identify cases in which an edge has two or more directions
        !flagged for elimination:
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
      function quad_int(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: order
      real(8)    :: x0,x1,x2,x3,y0,y1,y2,y3,x,y

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
     .                  ,arr0,icomp,is_cnv,is_vec,result_is_vec
     .                  ,out_iscnv,out_isvec,iorder)
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
c       *out_isvec (intent OUT): specifies whether input is curvilinear
c                                vector (needed when is_vec is inherited
c                                from a previous call).
c       *out_iscnv (intent OUT): specifies representation of output
c                                vector (needed when is_cnv is inherited
c                                from a previous call).
c       *iorder: order of inter/extrapolation in BC routines (quadratic
c                by default).
c
c     On output, array is returned either in the same curvilinear
c     representation or in Cartesian, if specified by "result_is_vec".
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
      logical   ,optional,intent(OUT):: out_iscnv,out_isvec

c     Local variables

      integer    :: i,j,k,stencil_width,order,ieq,neql
     .             ,imng,imxg,jmng,jmxg,kmng,kmxg
     .             ,imnl,imxl,jmnl,jmxl,kmnl,kmxl

      integer,save :: ivar  !Needed to "remember' ivar from previous calls
                               !when icomp is not provided.

      logical,save :: iscnv,isvec !Needed to "remember' these from previous calls
                                     !when they are not provided.

      logical      :: res_is_vec,fpointers

      real(8)      :: car(3),curv(3)

c     Begin program

      neql = neq !To avoid problems if neq is hardwired

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

      if (PRESENT(out_iscnv)) out_iscnv=iscnv
      if (PRESENT(out_isvec)) out_isvec=isvec

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

      call allocPointers(neql,fpointers)

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

      select case (neql)
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
     .            ,bcnd(:,neql),iig,iig,iig
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

        do ieq=1,neql
          call setBC(ivar,nnx,nny,nnz,array(:,:,:,ieq),vv0(:,:,:,ieq)
     .            ,bcnd(:,ieq),iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,iorder=order)
        enddo

        deallocate(vv0)

      case(3,6)

        allocate(v_cnv (0:nnx+1,0:nny+1,0:nnz+1,neql)
     .          ,v_cov (0:nnx+1,0:nny+1,0:nnz+1,neql)
     .          ,vv0   (0:nnx+1,0:nny+1,0:nnz+1,neql))

cc        v_cnv(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
cc        v_cov(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
        v_cnv = 0d0
        v_cov = 0d0

        if (PRESENT(arr0)) then
cc          vv0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
cc     .       = arr0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:)
          vv0 = arr0
        else
cc          vv0(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,:) = 0d0
          vv0 = 0d0
        endif

        if (.not.isvec) then

          do k = kmnl-1,kmxl+1
            do j = jmnl-1,jmxl+1
              do i = imnl-1,imxl+1
                do ieq=1,neql/3
                  car = array(i,j,k,1+3*(ieq-1):3*ieq)
                  call transformVectorToCurvilinear
     .               (i,j,k,iig,iig,iig
     .               ,car(1),car(2),car(3)
     .               ,.not.iscnv
     .               ,array(i,j,k,1+3*(ieq-1))
     .               ,array(i,j,k,2+3*(ieq-1))
     .               ,array(i,j,k,3+3*(ieq-1)))
                enddo
              enddo
            enddo
          enddo

        endif

        if (PRESENT(arr_cov)) then
          v_cnv = array
          v_cov = arr_cov
        else
          if (iscnv) then
            v_cnv = array
          else
            v_cov = array
          endif
        endif

        do ieq=1,neql/3
          call setBC(ivar,3,nnx,nny,nnz
     .            ,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .            ,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .            ,vv0  (:,:,:,1+3*(ieq-1):3*ieq)
     .            ,bcnd (:    ,1+3*(ieq-1):3*ieq)
     .            ,iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,is_cnv=iscnv,is_vec=.true.,iorder=order)
        enddo

cc        if (isvec .or. res_is_vec) then
        if (res_is_vec) then  !Overrule isvec when res_is_vec is specified

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
                do ieq=1,neql/3
                  curv = v_cnv(i,j,k,1+3*(ieq-1):3*ieq)
                  call transformVectorToCartesian
     .                 (i,j,k,iig,iig,iig
     .                 ,curv(1),curv(2),curv(3)
     .                 ,.false.
     .                 ,array(i,j,k,1+3*(ieq-1))
     .                 ,array(i,j,k,2+3*(ieq-1))
     .                 ,array(i,j,k,3+3*(ieq-1)))
                enddo
              enddo
            enddo
          enddo

        endif

        deallocate(v_cov,v_cnv,vv0)

      case default
        write (*,*) 'Number of equations',neql
     .             ,'not implemented in setMGBC'
        write (*,*) 'Aborting...'
        stop
      end select

c     Deallocate MG pointers (if allocated in this routine)

      call deallocPointers(fpointers)

c     End

      end subroutine setMGBC

      end module setMGBC_interface

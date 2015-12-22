c module bc_ops
c #####################################################################
      module bc_ops

c     -----------------------------------------------------------------
c     This module defines various BC-related operators
c     -----------------------------------------------------------------

        use grid_def

        type(grid_mg_def),pointer :: gbc_def => null()

      contains

c     isBdry
c     #################################################################
      function isBdry(i,igr,ibc)

c     -----------------------------------------------------------------
c     Detects whether we are at a PHYSICAL boundary or not. Boundary is
c     identified by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,igr,ibc
        logical :: isBdry

c     Local variables

        integer :: ilog,itst,diml,locl

        type(grid_mg_def),pointer :: gl_def

c     Begin program

        locl = abs(mod(ibc,2)-2)-1
        diml = 1 + (ibc-1-locl)/2

        if (associated(gbc_def)) then
          gl_def => gbc_def
        else
          gl_def => grid_params
        endif

        select case(diml)
        case(1)
          ilog = gl_def%ilo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = gl_def%nxgl(igr)
          end select
        case(2)
          ilog = gl_def%jlo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = gl_def%nygl(igr)
          end select
        case(3)
          ilog = gl_def%klo(igr)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = gl_def%nzgl(igr)
          end select
        end select

        isBdry = (i+ilog-1 == itst)

      end function isBdry

c     isBC
c     #################################################################
      function isBC(i,igr,ibc,BC,bcnd)

c     -----------------------------------------------------------------
c     Checks whether point i at grid level igr, boundary location
c     ibc=(1,6) (see isBdry for definition) is of boundary type BC.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,igr,ibc,BC
        logical :: isBC
        integer,optional :: bcnd(6)

c     Local variables

        integer :: bcc(6)

c     Begin program

        if (ibc < 1 .or. ibc > 6) then
          isBC = .false.
          return
        endif

        if (PRESENT(bcnd)) then
          bcc = bcnd
        else
          bcc = bcond
        endif

        isBC = (isBdry(i,igr,ibc) .and. abs(bcc(ibc)) == BC)

      end function isBC

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        implicit none

        integer    :: i,j,k,igx,igy,igz
        logical    :: sing_point

c     Local variables

        type(grid_mg_def),pointer :: gl_def

c     Begin program

        if (associated(gbc_def)) then
          gl_def => gbc_def
        else
          gl_def => grid_params
        endif

        sing_point = (i+gl_def%ilo(igx)-1==1).and.bcSP()

      end function isSP

c     isSP2
c     #################################################################
      function isSP2(i,igr,dim,loc,ibc)

c     -----------------------------------------------------------------
c     Checks whether point at grid level igr, dimension dim(=1,2,3)
c     and location loc(=0,1) is at SP boundary.
c     -----------------------------------------------------------------

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

c     -----------------------------------------------------------------
c     Checks whether point at grid level igr, dimension dim(=1,2,3)
c     and location loc(=0,1) is at SYM boundary
c     -----------------------------------------------------------------

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

c     isSYMBdry
c     #################################################################
      function isSYMBdry(dim,loc,ibc)

c     -----------------------------------------------------------------
c     Checks whether boundary at dimension dim(=1,2,3)
c     and location loc(=0,1) is a SYM boundary
c     -----------------------------------------------------------------

        implicit none

        logical    :: isSYMBdry
        integer,optional :: dim,loc,ibc

        integer    :: ibcl,ilog,itst,diml,locl,igr

        type(grid_mg_def),pointer :: gl_def

        if (PRESENT(dim) .and.PRESENT(loc)) then
          diml = dim
          locl = loc
          ibcl = (1+locl)+2*(diml-1)
        elseif (PRESENT(ibc)) then
          ibcl = ibc
          locl = abs(mod(ibcl,2)-2)-1
          diml = 1 + (ibcl-1-locl)/2
        else
          write (*,*) 'Error in isSYMBdry'
          write (*,*) 'Aborting...'
          stop
        endif

        if (associated(gbc_def)) then
          gl_def => gbc_def
        else
          gl_def => grid_params
        endif

        igr = 1 !Use fine grid for test

        select case(diml)
        case(1)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = gl_def%nxv(igr)
          end select
        case(2)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = gl_def%nyv(igr)
          end select
        case(3)
          select case(locl)
          case (0)
            itst = 1
          case (1)
            itst = gl_def%nzv(igr)
          end select
        end select

        isSYMBdry = isBC(itst,igr,ibcl,SYM)

      end function isSYMBdry

      end module bc_ops

c module BCS_variables
c####################################################################
      module BCS_variables

        use oned_int

        use grid_mpi

        use grid_operations

        use math

        use bc_ops

        integer :: nxbc,nybc,nzbc,igxbc,igybc,igzbc

        logical :: bc_debug=.false.

        integer,parameter :: VECTOR=1,SCALAR=0,AUX=-1

        real(8),pointer,private,dimension(:,:,:)   :: zeros
        real(8),pointer,private,dimension(:,:,:,:) :: vzeros

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

      type(grid_mg_def),pointer :: gl_def

c     Begin program

      if (associated(gbc_def)) then
        gl_def => gbc_def
      else
        gl_def => grid_params
      endif

      ilog = gl_def%ilo(igxbc)
      jlog = gl_def%jlo(igybc)

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

      nullify(gl_def)

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

      type(grid_mg_def),pointer :: gl_def

c     External

c     Begin program

      if (associated(gbc_def)) then
        gl_def => gbc_def
      else
        gl_def => grid_params
      endif

      ilog = gl_def%ilo(igxbc)
      jlog = gl_def%jlo(igybc)

      if (npy < 1) call pstop('vectorSingularBC','npy < 1')

      nyg = gl_def%nygl(igybc)

c     Check that we are domain adjacent to singular point

      if (ilog /= 1) then
        return
      else
        i = ilog
      endif

c     Parallel Allgather in SP communicator

#if defined(petsc)
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

      nullify(gl_def)

      end subroutine vectorSingularBC

      end module singularBCinterface

c module imposeBCinterface
c #####################################################################
      module imposeBCinterface

        use singularBCinterface

        INTERFACE setBC
          module procedure imposeBConScalar,imposeBConVector
        end INTERFACE

        INTERFACE quad_int
          module procedure quad_int_0,quad_int_1,quad_int_2
        end INTERFACE

        INTERFACE XferBCs
          module procedure XferBCs_scl, XferBCs_vec
        END INTERFACE 

        integer  :: iimax,iimin,jjmax,jjmin,kkmax,kkmin
        integer  :: iimxp,iimnm,jjmxp,jjmnm,kkmxp,kkmnm

        integer,private :: bctype

        logical,private :: skip_corners=.false.
     .                    ,force_corners=.false.

        !Define order of interpolation for EDGE and CORNER treatment
        !(linear to preserve nearest-neighbors stencil; this is 
        ! important for colored diagonal formation)
        integer,private :: order=1

        private :: setASMflags,interpolateToFace
     .            ,fillEdges,fillCorners

      contains

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

      force_corners = fc

c     End program

      end subroutine set_force_corners

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

        type(grid_mg_def),pointer :: gl_def

c     Begin program

        if (associated(gbc_def)) then
          gl_def => gbc_def
        else
          gl_def => grid_params
        endif

        do__nothing = .false.

        if (PRESENT(asm_bc)) then
          asmbc = asm_bc
        else
          asmbc = asm   !Setup externally
        endif

        if (asmbc) return      !Return if additive Schwartz method is used (PC)

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

      if (bc_debug) then
        write (*,*) 'Proc ',my_rank,': asm    =',asm
        write (*,*) 'Proc ',my_rank,': asm_dir=',asm_dir
      endif

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

      type(grid_mg_def),pointer :: gl_def

c     Begin program

      if (associated(gbc_def)) then
        gl_def => gbc_def
      else
        gl_def => grid_params
      endif

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
            dn=do__nothing(dim,loc,igmin,igmax,jgmin,jgmax,kgmin,kgmax
     .                    ,asm_bc=.false.)
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

      nullify(gl_def)

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

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(ieq,nx,ny,nz,array,array0,bcnd
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder,g_def)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ieq,bcnd(6),nx,ny,nz,iigx,iigy,iigz
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)
      integer   ,optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

      type(grid_mg_def),optional,pointer :: g_def

c     Local variables

      integer :: dim,loc,ibc
      logical :: petsc__bc

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

      call set_skip_corners(.false.) !Skip BC corner treatment

c     Initialize interpolation order (else use last set)

      if (PRESENT(iorder)) then
        order = iorder
cc      else
cc        order = 1
      endif

c     Initialize BC grid quantities

      if (PRESENT(g_def)) then
        gbc_def => g_def
      else
        gbc_def => grid_params
      endif

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

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc =  (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => gbc_def%dactx(igxbc)
        call fillPetscGhostCells(da_ctx,array)

        if (bcPER(1).and.bcPER(2).and.bcPER(3)) return !PETSc deals with all periodic domains
      endif
#endif

c     Impose BCs

#if defined(petsc) || defined(samrai)
      do bctype=2,BCLIM             !Do not do PER=1; already done by PETSc or SAMRAI
#else
      do bctype=1,BCLIM             !Enforces a particular order in the BCs
#endif
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (abs(bcnd(ibc)) == bctype) then
              call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
            endif
          enddo
        enddo
      enddo

c     Fill edges

      call fillEdges(array,array0,bcnd,order)

c     Fill corners

      call fillCorners(array,bcnd,order)

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_PER
#endif

c     End program

      nullify(gbc_def)

      contains

c     sync_PER
c     #################################################################
      subroutine sync_PER

c     -----------------------------------------------------------------
c     Synchronizes periodic boundaries.
c     -----------------------------------------------------------------

      implicit none

c     Local variables

c     Begin program

      bctype=PER
      do dim=1,3
        do loc=0,1
          ibc = (1+loc)+2*(dim-1)
          if (abs(bcnd(ibc)) == bctype) then
            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
          endif
        enddo
      enddo

      end subroutine sync_PER

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(fcomp,ncomp,nx,ny,nz,v_cnv,v_cov,var0
     .                          ,bcnd,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                          ,iorder,is_cnv,is_vec,g_def)
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

      type(grid_mg_def),optional,pointer :: g_def

c     Local variables

      integer    :: i,j,k,imax,imin,jmax,jmin,kmax,kmin
      integer    :: ivar,ieq,ibc,loc,dim,bctype
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: cov_to_cnv(6,3),to_cnv,is__vec,petsc__bc

#if defined(petsc)
      type(petsc_da_ctx),pointer :: da_ctx
#endif

c     Begin program

      call set_skip_corners(.false.) !Skip BC corner treatment

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

c     Initialize interpolation order (else use last set)

      if (PRESENT(iorder)) then
        order = iorder
cc      else
cc        order = 1
      endif

c     Initialize BC grid quantities

      if (PRESENT(g_def)) then
        gbc_def => g_def
      else
        gbc_def => grid_params
      endif

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

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc = (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => gbc_def%dactx(igxbc)
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

c     Find complete curvilinear representation at boundaries

      if (is__vec) then
        if (to_cnv) then
          call XformVector_BC(igxbc,v_cov,v_cnv,to_cnv)
        else
          call XformVector_BC(igxbc,v_cnv,v_cov,to_cnv)
        endif
      endif

c     Exit if all periodic domains (PETSc deals with all periodic domains)

#if defined(petsc)
      if (petsc__bc.and.bcPER(1).and.bcPER(2).and.bcPER(3)) return
#endif

c     Impose BCs

      cov_to_cnv = .false.

#if defined(petsc) || defined(samrai)
      do bctype=2,BCLIM             !Do not do PER=1; already done by PETSc or SAMRAI
#else
      do bctype=1,BCLIM             !Enforces a particular order in the BCs
#endif
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              ieq = ivar + fcomp - 1

              if (abs(bcnd(ibc,ivar)) == bctype) then

                if (bcnd(ibc,ivar) < 0 .and. is__vec) then
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cov,var0,cnv=.false.)
                  cov_to_cnv(ibc,ivar) = .true.
                else
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cnv,var0)
                endif

              endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components

      call synchronize(ncomp,v_cnv,v_cov,bcnd,cov_to_cnv)

c     Fill edges (for both cnv and cov components)

      do ivar=1,ncomp
        call fillEdges(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar),order)
        call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar),order)
      enddo

c     Fill corners (for both cnv and cov components)

      do ivar=1,ncomp
        call fillCorners(v_cov(:,:,:,ivar),bcnd(:,ivar),order)
        call fillCorners(v_cnv(:,:,:,ivar),bcnd(:,ivar),order)
      enddo

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_PER
#endif

c     End program

      nullify(gbc_def)

      contains

c     sync_PER
c     #################################################################
      subroutine sync_PER

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

      end subroutine sync_PER

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

c     Begin program

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

          if (.not.selectFace(ibc)) cycle !Do nothing for logical boundaries

          !Do NOT change: SAMRAI depends on it
          call findBCLoopLimits(id,il,iimin-1,iimax+1
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
                  call transformFromCurvToCurv(i,j,k
     .                  ,igxbc,igybc,igzbc
     .                  ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                  ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                  ,.false.)
                enddo
              enddo
            enddo

          else   !Synchronize (assumes normal component is contravariant)

cc            if (.not.selectFace(ibc)) cycle !Do nothing for logical boundaries

            do ivar = 1,ncomp

              if (ivar == id) then !Do nothing for normal component

                cycle

              else

cc                if (.not.selectFace(ibc)) cycle

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

                        call transformFromCurvToCurv(i,j,k
     .                    ,igxbc,igybc,igzbc
     .                    ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                    ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                    ,.true.)

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

                        call transformFromCurvToCurv(i,j,k
     .                    ,igxbc,igybc,igzbc
     .                    ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .                    ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3)
     .                    ,.false.)

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

      end subroutine imposeBConVector

c     imposeBConVector_stg
c     #################################################################
      subroutine imposeBConVector_stg(fcomp,ncomp,nx,ny,nz,v_cnv,v_cov
     .                     ,var0,bcnd,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                     ,iorder,is_cnv,is_vec,g_def)
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
      logical,optional,intent(IN) :: is_cnv,is_vec

      type(grid_mg_def),optional,pointer :: g_def

c     Local variables

      integer :: i,j,k,imax,imin,jmax,jmin,kmax,kmin!,order
      integer :: ivar,ieq,ibc,loc,dim,bctype
     .          ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical :: cov_to_cnv(6,3),to_cnv,is__vec,petsc__bc

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

      call set_skip_corners(.true.) !Skip BC corner treatment

c     Initialize interpolation order (else use last set)

      if (PRESENT(iorder)) then
        order = iorder
cc      else
cc        order = 1
      endif

c     Initialize BC grid quantities

      if (PRESENT(g_def)) then
        gbc_def => g_def
      else
        gbc_def => grid_params
      endif

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

c     Fill PETSc ghost cells

#if defined(petsc) && !defined(samrai)
      !Determine whether we need to fill ghost cells using PETSc
      petsc__bc = (.not.asm)

      !Fill ghost cells
      if (petsc__bc) then
        da_ctx => gbc_def%dactx(igxbc)
        do ieq=1,ncomp
          call fillPetscGhostCells(da_ctx,v_cov(:,:,:,ieq))
          call fillPetscGhostCells(da_ctx,v_cnv(:,:,:,ieq))
        enddo
      endif
#endif

c     Exit if all periodic domains (PETSc deals with all periodic domains)

#if defined(petsc)
      if (petsc__bc.and.bcPER(1).and.bcPER(2).and.bcPER(3)) return
#endif

c     Impose BCs

      cov_to_cnv = .false.

#if defined(petsc) || defined(samrai)
      do bctype=2,BCLIM             !Do not do PER=1; already done by PETSc or SAMRAI
#else
      do bctype=1,BCLIM             !Enforces a particular order in the BCs
#endif
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,ncomp
              ieq = ivar + fcomp - 1

              if (abs(bcnd(ibc,ivar)) == bctype) then

                if (bcnd(ibc,ivar) < 0 .and. is__vec) then
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cov,var0,cnv=.false.)
                  cov_to_cnv(ibc,ivar) = .true.
                else
                  call FillGhostNodes(ieq,ivar,ncomp,dim,loc,bctype
     .                               ,v_cnv,var0)
                endif

              endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components

      call synchronize_stg(ncomp,v_cnv,v_cov,bcnd,cov_to_cnv)

c     Fill edges (for both cnv and cov components)

      do ivar=1,ncomp
        call fillEdges(v_cov(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar),order)
        call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                ,bcnd(:,ivar),order)
      enddo

c     Fill corners (for both cnv and cov components)

      do ivar=1,ncomp
        call fillCorners(v_cov(:,:,:,ivar),bcnd(:,ivar),order)
        call fillCorners(v_cnv(:,:,:,ivar),bcnd(:,ivar),order)
      enddo

c     Periodic boundary synchronizing

#if defined(PER_BC_SYNC) && !defined(samrai)
      call sync_PER
#endif

c     End program

      nullify(gbc_def)

      contains

c     sync_PER
c     #################################################################
      subroutine sync_PER

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

      end subroutine sync_PER

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
      real(8) :: v_cnv(0:nxbc+1,0:nybc+1,0:nzbc+1,ncomp)
     .          ,v_cov(0:nxbc+1,0:nybc+1,0:nzbc+1,ncomp)
      logical :: cov_to_cnv(6,3)
      integer,optional :: dim,loc

c     Local variables

      integer :: imax,imin,jmax,jmin,kmax,kmin
      integer :: i,j,k,ig,jg,kg,ivar,ibc,ic1,ic2
     .          ,ld,hd,ll,hl,id,il,im,jm,km
      real(8) :: x1,x2,x3,gsup,gsub,cnv(3),cov(3)

c     Begin program

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

          if (.not.selectFace(ibc).or.bcPER(id)) cycle !Do nothing for logical or periodic boundaries

          !Do NOT change: SAMRAI depends on it
          call findBCLoopLimits(id,il,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                               ,imin   ,imax
     .                               ,jmin   ,jmax
     .                               ,kmin   ,kmax)

          if (     (.not.cov_to_cnv(ibc,1))
     .        .and.(.not.cov_to_cnv(ibc,2))
     .        .and.(.not.cov_to_cnv(ibc,3))) then           !All contravariant: find covariant

            !This code is suspect
cc            call pstop('sync_stg','Code is suspect')
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax
                   v_cov(i,j,k,:)= XformToCov_stg_ijk(i,j,k,igxbc,v_cnv)
                enddo
              enddo
            enddo

          else   !Synchronize (assumes normal component is contravariant)

cc            if (.not.selectFace(ibc)) cycle !Do nothing for logical boundaries

            do ivar = 1,ncomp

              if (ivar == id) then !Do nothing for normal component

                cycle

              else

cc                if (.not.selectFace(ibc)) cycle

                ic1 = mod(id+1,3)
                if (ic1 == 0) ic1 = 3

                ic2 = mod(id+2,3)
                if (ic2 == 0) ic2 = 3

                if ( cov_to_cnv(ibc,ic1).and.cov_to_cnv(ibc,ic2) ) then  !Both are covariant

                  do k=kmin,kmax
                    do j=jmin,jmax
                      do i=imin,imax

                        im = max(i-1,0)
                        jm = max(j-1,0)
                        km = max(k-1,0)

cc                        if(id/=3.and.(i>nxbc.or.j>nybc.or.k>nzbc)) cycle  !2D fix for now

                        select case(id)
                        case(1)
                          gsup =
     .                    0.5*(gmetric%grid(igxbc)%gsup_v(i,j ,k,id,id)
     .                        +gmetric%grid(igxbc)%gsup_v(i,jm,k,id,id))
                        case(2)
                          gsup =
     .                    0.5*(gmetric%grid(igxbc)%gsup_v(im,j,k,id,id)
     .                        +gmetric%grid(igxbc)%gsup_v(i ,j,k,id,id))
                        case(3)
                          gsup = 1d0
                        end select

                        v_cov(i,j,k,id) = 0d0

                        cnv = XformToCnv_stg_ijk(i,j,k,igxbc,v_cov)
                        v_cov(i,j,k,id) =-cnv(id) - v_cnv(i,j,k,id)/gsup

                        v_cnv(i,j,k,:) = XformToCnv_stg_ijk(i,j,k,igxbc
     .                                                     ,v_cov)

                      enddo
                    enddo
                  enddo

                elseif (.not.cov_to_cnv(ibc,ic1)) then   !ic1 is contravariant, ic2 is covariant

                  call pstop('sync_stg','Code untested')

                  do k=kmin,kmax
                    do j=jmin,jmax
                      do i=imin,imax

                        im = max(i-1,0)
                        jm = max(j-1,0)
                        km = max(k-1,0)

cc                        if(ic2/=3.and.(i>nxbc.or.j>nybc.or.k>nzbc))cycle  !2D fix for now

                        select case(ic2)
                        case(1)
                          gsub =
     .                   .5*(gmetric%grid(igxbc)%gsub_v(i,j ,k,ic2,ic2)
     .                      +gmetric%grid(igxbc)%gsub_v(i,jm,k,ic2,ic2))
                        case(2)
                          gsub =
     .                   .5*(gmetric%grid(igxbc)%gsub_v(im,j,k,ic2,ic2)
     .                      +gmetric%grid(igxbc)%gsub_v(i ,j,k,ic2,ic2))
                        case(3)
                          gsub = 1d0
                        end select

                        v_cnv(i,j,k,ic2) = 0d0
                        cov = XformToCov_stg_ijk(i,j,k,igxbc,v_cnv)
                        v_cnv(i,j,k,ic2)=-cov(ic2)-v_cov(i,j,k,ic2)/gsub

                        v_cov(i,j,k,:) = XformToCov_stg_ijk(i,j,k,igxbc
     .                                                     ,v_cnv)

                      enddo
                    enddo
                  enddo

                elseif (.not.cov_to_cnv(ibc,ic2)) then   !ic2 is contravariant, ic1 is covariant

                  call pstop('sync_stg','Code untested')

                  do k=kmin,kmax
                    do j=jmin,jmax
                      do i=imin,imax

                        im = max(i-1,0)
                        jm = max(j-1,0)
                        km = max(k-1,0)

cc                        if(ic1/=3.and.(i>nxbc.or.j>nybc.or.k>nzbc))cycle

                        select case(ic1)
                        case(1)
                          gsub =
     .                   .5*(gmetric%grid(igxbc)%gsub_v(i,j ,k,ic1,ic1)
     .                      +gmetric%grid(igxbc)%gsub_v(i,jm,k,ic1,ic1))
                        case(2)
                          gsub =
     .                   .5*(gmetric%grid(igxbc)%gsub_v(im,j,k,ic1,ic1)
     .                      +gmetric%grid(igxbc)%gsub_v(i ,j,k,ic1,ic1))
                        case(3)
                          gsub = 1d0
                        end select

                        v_cnv(i,j,k,ic1) = 0d0
                        cov = XformToCov_stg_ijk(i,j,k,igxbc,v_cnv)
                        v_cnv(i,j,k,ic1)=-cov(ic1)-v_cov(i,j,k,ic1)/gsub

                        v_cov(i,j,k,:) = XformToCov_stg_ijk(i,j,k,igxbc
     .                                                     ,v_cnv)

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

      end subroutine synchronize_stg

      end subroutine imposeBConVector_stg

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

      integer    :: neq,ibc,i,j,k,ig,jg,kg
      logical    :: cov

      real(8),allocatable,dimension(:,:,:,:) :: zeros

c     Begin program

      if (PRESENT(cnv)) then
        cov = .not.cnv
      else
        cov = .false.
      endif

      ibc = (1+loc)+2*(dim-1)
      
c     Determine boundary limits

      if (.not.selectFace(ibc)) return

c     Find BC update

      if (nvar == 1) then
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,nvar),dim,loc)
        case(EQU)
          call extrapolateBC(array (:,:,:,nvar)
     $                      ,array0(:,:,:,nvar)
     .                      ,dim,loc,order)
        case(DIR)
          allocate(zeros(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar))
          zeros = 0d0
          call scalarDirichletBC(array(:,:,:,nvar)
     .                          ,zeros(:,:,:,nvar)
cc     .                          ,array0(:,:,:,nvar)
     .                          ,ieq,dim,loc,order)
          deallocate(zeros)
        case(NEU)
          call scalarNeumannBC(array (:,:,:,nvar)
     .                        ,array0(:,:,:,nvar)
     .                        ,ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,nvar),dim,loc,order)
        case(SP)
          call scalarSingularBC(array(:,:,:,nvar))
        case(IFC)
          call interpolateToFace(array(:,:,:,nvar),dim,loc,order)
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
          allocate(zeros(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar))
          zeros = 0d0
          call vectorDirichletBC(ivar,array,zeros,ieq,dim,loc,order)
cc          call vectorDirichletBC(ivar,array,array0,ieq,dim,loc,order)
          deallocate(zeros)
        case(NEU)
          call vectorNeumannBC(ivar,array,array0,ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,ivar),dim,loc,order)
        case(SP)
          call vectorSingularBC(ivar,array(:,:,:,ivar),cov)
        case(IFC)
          call interpolateToFace(array(:,:,:,ivar),dim,loc,order)
        case default
          call pstop('fillGhostNodes'
     .             ,'BC '//trim(int2char(bctype))//' not implemented')
        end select
      endif

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

      integer :: dim,loc
      real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer :: ibc

      integer :: np_l,my_rank_l,dest,tag=0

      real(8),pointer,dimension(:,:) :: rhs

c     Begin program

      ibc = (1+loc)+2*(dim-1)

#if !defined(petsc)

      call per_BC(ibc,array)

#else

      if (MPI_COMM_PER(dim) /= MPI_COMM_NULL) then
        call MPI_Comm_size(MPI_COMM_PER(dim),np_l     ,mpierr)
        call MPI_Comm_rank(MPI_COMM_PER(dim),my_rank_l,mpierr)

        if (np_l > 2.or.np_l < 1) then
          call pstop('sync_PER'
     .         ,'Periodic comm has wrong rank='//int2char(np_l))
        elseif (np_l == 2) then
          dest = 1 - my_rank_l
        endif

        if (np_l == 1) then

          call per_BC(ibc,array)

        else
          
          select case(dim)
          case(1)

            allocate(rhs(size(array,2),size(array,3)))

            if (loc==0.and.isBdry(1,igxbc,ibc)) then
cc              write (*,*) 'DIAG -- sync_PER   ;   '
cc     .         ,'dim=',dim,'my_rank=',my_rank,'Ranks=',np_l
cc     .         ,'Source=',my_rank_l,'Dest=',dest,'left'

              rhs = array(1,:,:)
              call MPI_Sendrecv_replace(rhs,size(array,2)*size(array,3)
     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                   ,MPI_COMM_PER(dim),MPI_STATUS_IGNORE,mpierr)
              array(0,:,:) = rhs
            elseif (loc==1.and.isBdry(nxbc,igxbc,ibc)) then
cc              write (*,*) 'DIAG -- sync_PER   ;   '
cc     .         ,'dim=',dim,'my_rank=',my_rank,'Ranks=',np_l
cc     .         ,'Source=',my_rank_l,'Dest=',dest,'right'

              rhs = array(nxbc,:,:)
              call MPI_Sendrecv_replace(rhs,size(array,2)*size(array,3)
     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                   ,MPI_COMM_PER(dim),MPI_STATUS_IGNORE,mpierr)
              array(nxbc+1,:,:) = rhs
            endif

            deallocate(rhs)

          case(2)

            allocate(rhs(size(array,1),size(array,3)))

            if (loc==0.and.isBdry(1,igybc,ibc)) then
              rhs = array(:,1,:)
              call MPI_Sendrecv_replace(rhs,size(array,1)*size(array,3)
     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                   ,MPI_COMM_PER(dim),MPI_STATUS_IGNORE,mpierr)
              array(:,0,:) = rhs
            elseif (loc==1.and.isBdry(nybc,igybc,ibc)) then
              rhs = array(:,nybc,:)
              call MPI_Sendrecv_replace(rhs,size(array,1)*size(array,3)
     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                   ,MPI_COMM_PER(dim),MPI_STATUS_IGNORE,mpierr)
              array(:,nybc+1,:) = rhs
            endif

            deallocate(rhs)

          case(3)

            allocate(rhs(size(array,1),size(array,2)))

            if (loc==0.and.isBdry(1,igzbc,ibc)) then
              rhs = array(:,:,1)
              call MPI_Sendrecv_replace(rhs,size(array,1)*size(array,2)
     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                   ,MPI_COMM_PER(dim),MPI_STATUS_IGNORE,mpierr)
              array(:,:,0) = rhs
            elseif (loc==1.and.isBdry(nzbc,igzbc,ibc)) then
              rhs = array(:,:,nzbc)
              call MPI_Sendrecv_replace(rhs,size(array,1)*size(array,2)
     .                   ,MPI_DOUBLE_PRECISION,dest,dim,dest,dim
     .                   ,MPI_COMM_PER(dim),MPI_STATUS_IGNORE,mpierr)
              array(:,:,nzbc+1) = rhs
            endif

            deallocate(rhs)
          end select
        endif

      endif
#endif

c     End program

      end subroutine periodicBC

c     extrapolateBC
c     #################################################################
      subroutine extrapolateBC(array,array0,dim,loc,order)
c     -----------------------------------------------------------------
c     Fills ghost nodes by extrapolation across relevant boundary.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: order,dim,loc
        real(8) :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .            ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

        integer :: imax,imin,jmax,jmin,kmax,kmin,ibc
     .            ,order_int(3)

        integer :: i,j,k,ig,jg,kg,one,two,thr
        real(8),pointer,dimension(:) :: xx,yy,zz

c     Begin program

        ibc = (1+loc)+2*(dim-1)

        order_int(1) = min(order,nxbc)
        order_int(2) = min(order,nybc)
        order_int(3) = min(order,nzbc)

        one = 1 ; two = 2; thr = 3

c     Setup grid aliases

        xx => gbc_def%xx
        yy => gbc_def%yy
        zz => gbc_def%zz

c     Extrapolate

        !Do not change: SAMRAI depends on it
        call findBCLoopLimits(dim,loc,iimin-1,iimax+1
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

            call getMGmap(i,1,1,igxbc,igybc,igzbc,ig,jg,kg
     .                   ,g_def=gbc_def)

            array(i,:,:) = 
     .              quad_int(0.5*(xx(ig+one)+xx(ig)),xx(ig+one)
     .                      ,xx(ig+two),xx(ig+thr)
     .                      ,array0(i    ,:,:),array(i+one,:,:)
     .                      ,array (i+two,:,:),array(i+thr,:,:)
     .                      ,xx(ig),order_int(dim) )

          case (2)

            i = imax

            call getMGmap(i,1,1,igxbc,igybc,igzbc,ig,jg,kg
     .                   ,g_def=gbc_def)

            array(i,:,:) =
     .              quad_int(0.5*(xx(ig-one)+xx(ig)),xx(ig-one)
     .                      ,xx(ig-two),xx(ig-thr)
     .                      ,array0(i    ,:,:),array(i-one,:,:)
     .                      ,array (i-two,:,:),array(i-thr,:,:)
     .                      ,xx(ig),order_int(dim) )

          case (3)

            j = jmin

            call getMGmap(1,j,1,igxbc,igybc,igzbc,ig,jg,kg
     .                   ,g_def=gbc_def)

            array(:,j,:) = quad_int(0.5*(yy(jg+one)+yy(jg)),yy(jg+one)
     .                      ,yy(jg+two),yy(jg+thr)
     .                      ,array0(:,j    ,:),array(:,j+one,:)
     .                      ,array (:,j+two,:),array(:,j+thr,:)
     .                      ,yy(jg),order_int(dim) )

          case (4)

            j = jmax

            call getMGmap(1,j,1,igxbc,igybc,igzbc,ig,jg,kg
     .                   ,g_def=gbc_def)

            array(:,j,:) = quad_int(0.5*(yy(jg-one)+yy(jg)),yy(jg-one)
     .                      ,yy(jg-two),yy(jg-thr)
     .                      ,array0(:,j    ,:),array(:,j-one,:)
     .                      ,array (:,j-two,:),array(:,j-thr,:)
     .                      ,yy(jg),order_int(dim) )

          case (5)

            k = kmin

            call getMGmap(1,1,k,igxbc,igybc,igzbc,ig,jg,kg
     .                   ,g_def=gbc_def)

            array(:,:,k) = quad_int(0.5*(zz(kg+one)+zz(kg)),zz(kg+one)
     .                      ,zz(kg+two),zz(kg+thr)
     .                      ,array0(:,:,k    ),array(:,:,k+one)
     .                      ,array (:,:,k+two),array(:,:,k+thr)
     .                      ,zz(kg),order_int(dim) )

          case (6)

            k = kmax

            call getMGmap(1,1,k,igxbc,igybc,igzbc,ig,jg,kg
     .                   ,g_def=gbc_def)

            array(:,:,k) =
     .              quad_int(0.5*(zz(kg-one)+zz(kg)),zz(kg-one)
     .                      ,zz(kg-two),zz(kg-thr)
     .                      ,array0(:,:,k    ),array(:,:,k-one)
     .                      ,array (:,:,k-two),array(:,:,k-thr)
     .                      ,zz(kg),order_int(dim) )

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
      subroutine interpolateToFace(array,dim,loc,order)
c     -----------------------------------------------------------------
c     Interpolates values to face for EQU BC. On input:
c        * array -> 3D array to impose BC's on.
c        * dim -> dimension we are imposing BC on (X,Y,Z)
c        * loc -> boundary location (0 -> left, 1->right)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

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

cc        if (associated(gbc_def)) then
cc          gl_def => gbc_def
cc        else
cc          gl_def => grid_params
cc        endif

        xx => gbc_def%xx
        yy => gbc_def%yy
        zz => gbc_def%zz

c     Extrapolate

        !Do NOT change: SAMRAI depends on it
        call findBCLoopLimits(dim,loc,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                       ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

        select case (ibc)
        case (1)

          i = imin

          call getMGmap(i,1,1,igxbc,igybc,igzbc,ig,jg,kg,g_def=gbc_def)

          array(i,:,:) = quad_int(xx(ig    ),xx(ig+one)
     .                      ,xx(ig+two),xx(ig+thr)
     .                      ,array(i    ,:,:),array(i+one,:,:)
     .                      ,array(i+two,:,:),array(i+thr,:,:)
cc     .                      ,xx(ig)+dxh(ig),order )
     .                      ,0.5*(xx(ig)+xx(ig+one)),order )

        case (2)

          i = imax

          call getMGmap(i,1,1,igxbc,igybc,igzbc,ig,jg,kg,g_def=gbc_def)

          array(i,:,:) =
     .              quad_int(xx(ig    ),xx(ig-one)
     .                      ,xx(ig-two),xx(ig-thr)
     .                      ,array(i    ,:,:),array(i-one,:,:)
     .                      ,array(i-two,:,:),array(i-thr,:,:)
cc     .                      ,xx(ig)-dxh(ig),order )
     .                      ,0.5*(xx(ig)+xx(ig-one)),order )

        case (3)

          j = jmin

          call getMGmap(1,j,1,igxbc,igybc,igzbc,ig,jg,kg,g_def=gbc_def)

          array(:,j,:) = quad_int(yy(jg    ),yy(jg+one)
     .                      ,yy(jg+two),yy(jg+thr)
     .                      ,array(:,j    ,:),array(:,j+one,:)
     .                      ,array(:,j+two,:),array(:,j+thr,:)
cc     .                      ,yy(jg)+dyh(jg),order )
     .                      ,0.5*(yy(jg)+yy(jg+one)),order )

        case (4)

          j = jmax

          call getMGmap(1,j,1,igxbc,igybc,igzbc,ig,jg,kg,g_def=gbc_def)

          array(:,j,:) = quad_int(yy(jg    ),yy(jg-one)
     .                      ,yy(jg-two),yy(jg-thr)
     .                      ,array(:,j    ,:),array(:,j-one,:)
     .                      ,array(:,j-two,:),array(:,j-thr,:)
cc     .                      ,yy(jg)-dyh(jg),order )
     .                      ,0.5*(yy(jg)+yy(jg-one)),order )

        case (5)

          k = kmin

          call getMGmap(1,1,k,igxbc,igybc,igzbc,ig,jg,kg,g_def=gbc_def)

          array(:,:,k) = quad_int(zz(kg    ),zz(kg+one)
     .                      ,zz(kg+two),zz(kg+thr)
     .                      ,array(:,:,k    ),array(:,:,k+one)
     .                      ,array(:,:,k+two),array(:,:,k+thr)
cc     .                      ,zz(kg)+dzh(kg),order )
     .                      ,0.5*(zz(kg)+zz(kg+one)),order )

        case (6)

          k = kmax

          call getMGmap(1,1,k,igxbc,igybc,igzbc,ig,jg,kg,g_def=gbc_def)

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
      subroutine fillFaces(array,dim,loc,order)
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

      integer    :: dim,loc,order
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

      !Do NOT change: SAMRAI depends on it
      call findBCLoopLimits(dim,loc,iimin-1,iimax+1
     .                             ,jjmin-1,jjmax+1
     .                             ,kkmin-1,kkmax+1
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

        ito   = imin
        ifrom = min(ito+1+order,nxbc)

        imn   = min(ito,ifrom)
        imx   = max(ito,ifrom)

        ibcl(1) = 1

        do j=jmin,jmax
          do k=kmin,kmax
            call extrapolateToCorner(ibcl,ito,j,k,imn,j,k
     .                ,imx,j,k,array(imn:imx,j:j,k:k),bcond,order)
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
            call extrapolateToCorner(ibcl,ito,j,k,imn,j,k
     .                ,imx,j,k,array(imn:imx,j:j,k:k),bcond,order)
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
              call extrapolateToCorner(ibcl,i,jto,k,i,jmn,k
     .                ,i,jmx,k,array(i:i,jmn:jmx,k:k),bcond,order)
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
            call extrapolateToCorner(ibcl,i,jto,k,i,jmn,k
     .                ,i,jmx,k,array(i:i,jmn:jmx,k:k),bcond,order)
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
            call extrapolateToCorner(ibcl,i,j,kto,i,j,kmn
     .                ,i,j,kmx,array(i:i,j:j,kmn:kmx),bcond,order)
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
            call extrapolateToCorner(ibcl,i,j,kto,i,j,kmn
     .                ,i,j,kmx,array(i:i,j:j,kmn:kmx),bcond,order)
          enddo
        enddo

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
      subroutine fillEdges(array,array0,bcond,order)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * bcond  -> boundary condition array
c        * order  -> order of extrapolation
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: bcond(6),order
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

      ord = max(order,1)+1  !Does not affect order of interpolation in extrapolateToCorner

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
     .                ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order)
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
     .                ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order)
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
     .                ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order)
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
     .                ,imx,jmx,k,array(imn:imx,jmn:jmx,k:k),bcond,order)
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
     .                ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order)
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
     .                ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order)
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
     .                ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order)
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
     .               ,i,jmx,kmx,array(i:i,jmn:jmx,kmn:kmx),bcond,order)
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
     .                ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order)
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
     .               ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order)
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
     .               ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order)
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
     .               ,imx,j,kmx,array(imn:imx,j:j,kmn:kmx),bcond,order)
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
      subroutine fillCorners(array,bcond,order)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

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
     .          ,imn,imx,jmn,jmx,kmn,kmx,ord,ibcl(3)

      logical :: sel_f(3)

c     Begin program

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
     .             ,array(imn:imx,jmn:jmx,kmn:kmx),bcond,order)

      endif

      end subroutine fillCorner

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(ibcl,ito,jto,kto,imn,jmn,kmn
     .                             ,imx,jmx,kmx,array,bcnd,order,array0)
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
c        * array0: array specifying EQU BCs
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ibcl(3),ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx,bcnd(6)
     .          ,order
      real(8) :: array(imn:imx,jmn:jmx,kmn:kmx)
      real(8),optional,intent(IN) :: array0(imn:imx,jmn:jmx,kmn:kmx)

c     Local variables

      integer :: i,j,k,ig,jg,kg,inx,iny,inz,order_int(3),offset
     .          ,imni,imxi,jmni,jmxi,kmni,kmxi,num,i_bc,i1,i2

      real(8) :: valx(1),valy(1),valz(1),xto(1),yto(1),zto(1),curv(3)
     $          ,dx,dy,dz

      logical :: flg(3),select_bc(3),face,edge,corner

cc      type(grid_mg_def),pointer :: gl_def

c     Extrapolation

      real(8),dimension(10) :: xx,yy,zz

c     Begin program

cc      if (associated(gbc_def)) then
cc        gl_def => gbc_def
cc      else
cc        gl_def => grid_params
cc      endif

c     Initialize variables

      face   = (count(ibcl /= 0) == 1)
      edge   = (count(ibcl /= 0) == 2)
      corner = (count(ibcl /= 0) == 3)

      if (skip_corners.and.(.not.face)) return

c     Do nothing for SP boundary faces, edges, and corners

      if (   (bcSP(ibc=ibcl(1),bcnd=bcnd)
     .    .or.bcSP(ibc=ibcl(2),bcnd=bcnd)
     .    .or.bcSP(ibc=ibcl(3),bcnd=bcnd))
     .    .and.(.not.force_corners)) return

c     Do nothing for DEF boundary faces, edges, and corners

      if (   (bcChk(ibcl(1),DEF,bcnd=bcnd)
     .    .or.bcChk(ibcl(2),DEF,bcnd=bcnd)
     .    .or.bcChk(ibcl(3),DEF,bcnd=bcnd))
     .    .and.(.not.force_corners)) return

c     Do nothing for SYM (not FSYM) boundary faces, edges, and corners
c     (taken care of by scalarNeumann and scalarDirichlet routines)

      if (    isSYMBdry(ibc=ibcl(1))
     .    .or.isSYMBdry(ibc=ibcl(2))
     .    .or.isSYMBdry(ibc=ibcl(3))
     .    .and.(.not.force_corners)) return

c     Do nothing for doubly-periodic edges (PETSc does it)

      if (edge) then
        if ((ibcl(1) == 0)
     .      .and.(bcChk(ibcl(2),PER).and.bcChk(ibcl(3),PER))) return
        if ((ibcl(2) == 0)
     .      .and.(bcChk(ibcl(3),PER).and.bcChk(ibcl(1),PER))) return
        if ((ibcl(3) == 0)
     .      .and.(bcChk(ibcl(1),PER).and.bcChk(ibcl(2),PER))) return
      endif

c     Do nothing for Dirichlet edges (extrapolateBC does it)

      if (edge) then
        if ((ibcl(1) == 0)
     .  .and.((selectFace(ibcl(3)).and.(bcChk(ibcl(3),DIR,bcnd=bcnd)
     .                              .or.bcChk(ibcl(3),EQU,bcnd=bcnd)))
     .    .or.(selectFace(ibcl(2)).and.(bcChk(ibcl(2),DIR,bcnd=bcnd)
     .                              .or.bcChk(ibcl(2),EQU,bcnd=bcnd)))))
     .       return
        if ((ibcl(2) == 0)
     .  .and.((selectFace(ibcl(3)).and.(bcChk(ibcl(3),DIR,bcnd=bcnd)
     .                              .or.bcChk(ibcl(3),EQU,bcnd=bcnd)))
     .    .or.(selectFace(ibcl(1)).and.(bcChk(ibcl(1),DIR,bcnd=bcnd)
     .                              .or.bcChk(ibcl(2),EQU,bcnd=bcnd)))))
     .       return
        if ((ibcl(3) == 0)
     .  .and.((selectFace(ibcl(1)).and.(bcChk(ibcl(1),DIR,bcnd=bcnd)
     .                              .or.bcChk(ibcl(1),EQU,bcnd=bcnd)))
     .    .or.(selectFace(ibcl(2)).and.(bcChk(ibcl(2),DIR,bcnd=bcnd)
     .                              .or.bcChk(ibcl(2),EQU,bcnd=bcnd)))))
     .       return
      endif

c     Identify directions for extrapolation.

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
        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg
     .               ,g_def=gbc_def)

        if (inx > 2) then
          if (ito == imn) then
            dx = gbc_def%dx(ig+2)
            curv(1) = (  array(ito+1,jto,kto)
     $                -2*array(ito+2,jto,kto)
     $                  +array(ito+3,jto,kto))*dx**(-2)
          else
            dx = gbc_def%dx(ig-2)
            curv(1) = (  array(ito-1,jto,kto)
     $                -2*array(ito-2,jto,kto)
     $                  +array(ito-3,jto,kto))*dx**(-2)
          endif
        endif

        if (iny > 2) then
          if (jto == jmn) then
            dy = gbc_def%dy(jg+2)
            curv(2) = (  array(ito,jto+1,kto)
     $                -2*array(ito,jto+2,kto)
     $                  +array(ito,jto+3,kto))*dy**(-2)
          else                                              
            dy = gbc_def%dy(jg-2)
            curv(2) = (  array(ito,jto-1,kto)
     $                -2*array(ito,jto-2,kto)
     $                  +array(ito,jto-3,kto))*dy**(-2)
          endif
        endif

        if (inz > 2) then
          if (kto == kmn) then
            dz = gbc_def%dz(kg+2)
            curv(3) = (  array(ito,jto,kto+1)
     $                -2*array(ito,jto,kto+2)
     $                  +array(ito,jto,kto+3))*dz**(-2)
          else                                                  
            dz = gbc_def%dz(kg-2)
            curv(3) = (  array(ito,jto,kto-1)
     $                -2*array(ito,jto,kto-2)
     $                  +array(ito,jto,kto-3))*dz**(-2)
          endif
        endif

        curv = abs(curv)
      endif

      !Do not interpolate **along** direction if not enough points,
      !Directions 1,2,3 are defined as normal to faces.
      flg = .true.

      if (inx == 0) flg(1) = .false.
      if (iny == 0) flg(2) = .false.
      if (inz == 0) flg(3) = .false.

c     Treat special boundaries

      call select_along

      call select_across

c     Check smoothness if various interpolation directions are selected

      call select_smooth

c     Enforce symmetry extrapolation at corners

      do i=1,3
        i_bc = ibcl(i)
        if (i_bc <= 0) cycle
        select_bc(i)=(abs(bcond(i_bc))==FSYM.and.abs(bcnd(i_bc))==NEU)
     .               .and.selectFace(i_bc)
        if (select_bc(i)) order_int(i) = 0
      enddo

c     Check some direction for interpolation exists

      if (    (.not.flg(1)).and.(.not.flg(2)).and.(.not.flg(3))
     .    .or. bc_debug) then
        write (*,*) 'Extrapolation directions',flg
        write (*,*) 'Boundary type: Face=',face,' Edge=',edge
     .             ,' Corner=',corner
        write (*,*) 'Boundaries involved',ibcl
        write (*,*) 'Order of extrapolation',order
        write (*,*) 'Smoothness indicators',curv
        if (.not.bc_debug)
     .     call pstop('extrapolateToCorner'
     .            ,'No extrapolation directions left')
      endif

c     Initialize accumulators

      valx = 0d0
      valy = 0d0
      valz = 0d0

      num = 0

c     X extrapolation

      if (flg(1)) then

        if (ito == imn) then
          call getMGmap(imn+1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg
     .                 ,g_def=gbc_def)
          xx(1:inx) = gbc_def%xx(ig:ig+inx-1)
          imni = imn+1
          imxi = imx
        else
          call getMGmap(imx-1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg
     .                 ,g_def=gbc_def)
          xx(1:inx) = gbc_def%xx(ig-(inx-1):ig)
          imni = imn
          imxi = imx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg
     .               ,g_def=gbc_def)
        xto = gbc_def%xx(ig)

        call IntDriver1d(inx,xx(1:inx),array(imni:imxi,jto,kto)
     .                  ,1,xto,valx,order_int(1))

        num = num + 1

      endif

c     Y extrapolation

      if (flg(2)) then

        if (jto == jmn) then
          call getMGmap(imn,jmn+1,kmn,igxbc,igybc,igzbc,ig,jg,kg
     .                 ,g_def=gbc_def)
          yy(1:iny) = gbc_def%yy(jg:jg+iny-1)
          jmni = jmn+1
          jmxi = jmx
        else
          call getMGmap(imn,jmx-1,kmn,igxbc,igybc,igzbc,ig,jg,kg
     .                 ,g_def=gbc_def)
          yy(1:iny) = gbc_def%yy(jg-(iny-1):jg)
          jmni = jmn
          jmxi = jmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg
     .               ,g_def=gbc_def)
        yto = gbc_def%yy(jg)

        call IntDriver1d(iny,yy(1:iny),array(ito,jmni:jmxi,kto)
     .                  ,1,yto,valy,order_int(2))

        num = num + 1

      endif

c     Z extrapolation

      if (flg(3)) then

        if (kto == kmn) then
          call getMGmap(imn,jmn,kmn+1,igxbc,igybc,igzbc,ig,jg,kg
     .                 ,g_def=gbc_def)
          zz(1:inz) = gbc_def%zz(kg:kg+inz-1)
          kmni = kmn+1
          kmxi = kmx
        else
          call getMGmap(imn,jmn,kmx-1,igxbc,igybc,igzbc,ig,jg,kg
     .                 ,g_def=gbc_def)
          zz(1:inz) = gbc_def%zz(kg-(inz-1):kg)
          kmni = kmn
          kmxi = kmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg
     .               ,g_def=gbc_def)
        zto = gbc_def%zz(kg)

        call IntDriver1d(inz,zz(1:inz),array(ito,jto,kmni:kmxi)
     .                  ,1,zto,valz,order_int(3))

        num = num + 1

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
c         1) SP boundary
c         2) Interpolate-to-face boundary
c         3) Equilibrium boundary (removed 12/13/06, added back 11/10/08)
c         4) Discard internal processor boundaries (not real boundaries)
c     -----------------------------------------------------------------------

      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc <= 0) cycle
        select_bc(i) =(    (abs(bcnd(i_bc)) == SP )
     .                 .or.(abs(bcnd(i_bc)) == IFC)
     .                 .or.(abs(bcnd(i_bc)) == EQU)
     .                ).and.selectFace(i_bc)
      enddo

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
            if (curv(i1) <= curv(i2)) then
              select_bc(i1) = .false.
            else
              select_bc(i2) = .false.
            endif
          endif

          if (select_bc(i1)) then
            if (curv(i) <= curv(i1)) then
              cycle
            else
              flg(i) = .not.select_bc(i)
            endif
          else
            if (curv(i) <= curv(i2)) then
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
c         1) Periodic boundary
c         2) Symmetry boundary (NEU+order=0; SYM removed on 6/1/12)
c         3) Equilibrium boundary (removed 11/10/08; was interfering with
c            parallelization)
c         4) Discard internal processor boundaries (not real boundaries)
c     -----------------------------------------------------------------------

      select_bc = .false.
      do i=1,3
        i_bc = ibcl(i)
        if (i_bc <= 0) cycle
        select_bc(i) =(     (abs(bcnd(i_bc)) == NEU .and. order == 0)
cc     .                 .or.(abs(bcond(i_bc))== FSYM .and. abs(bcnd(i_bc)) == NEU)
     .                 .or. (    bcnd(i_bc)  == PER)
     .                ).and.selectFace(i_bc)
      enddo

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

c     select_smooth
c     ######################################################################
      subroutine select_smooth

      implicit none

c     -----------------------------------------------------------------------
c     Selects interpolation directions so that interpolation occurs along
c     boundaries with comparable smoothness
c     -----------------------------------------------------------------------

      do i=1,3

        !Identify other directions
        i1 = mod(i+1,3)
        i2 = mod(i+2,3)
        if (i1 == 0) i1 = 3
        if (i2 == 0) i2 = 3

        !Identify cases in which an edge has two or more directions
        !flagged for inclusion, and check for smoothness
        if (flg(i1).and.flg(i2)) then
          if (curv(i1) > 5*curv(i2)) then
            flg(i1) = .false.
          elseif (curv(i2) > 5*curv(i1)) then
            flg(i2) = .false.
          endif
        endif

      enddo

      end subroutine select_smooth

      end subroutine extrapolateToCorner

c     quad_int_0
c     #################################################################
      function quad_int_0(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: order
      real(8) :: x0,x1,x2,x3,y0,y1,y2,y3,x,y

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
        call pstop('quad_int_0',messg)
      end select

c     End program

      end function quad_int_0

c     quad_int_1
c     #################################################################
      function quad_int_1(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: order
      real(8) :: x0,x1,x2,x3,x
     .          ,y0(:),y1(:),y2(:),y3(:),y(size(y0))

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
        call pstop('quad_int_1',messg)
      end select

c     End program

      end function quad_int_1

c     quad_int_2
c     #################################################################
      function quad_int_2(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: order
      real(8) :: x0,x1,x2,x3,x
     .          ,y0(:,:),y1(:,:),y2(:,:),y3(:,:)
     .          ,y(size(y0,1),size(y0,2))

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
        call pstop('quad_int_2',messg)
      end select

c     End program

      end function quad_int_2

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
        arr2(nxp,:,:) = 0.5*(arr1(nxp  ,:,:)
     .                      +arr1(nxp-1,:,:))
        arr2(0  ,:,:) = 0.5*(arr1(1    ,:,:)
     .                      +arr1(0    ,:,:))

        arr2(:,nyp,:) = 0.5*(arr1(:,nyp  ,:)
     .                      +arr1(:,nyp-1,:))
        arr2(:,0  ,:) = 0.5*(arr1(:,1    ,:)
     .                      +arr1(:,0    ,:))

        arr2(:,:,nzp) = 0.5*(arr1(:,:,nzp  )
     .                      +arr1(:,:,nzp-1))
        arr2(:,:,0  ) = 0.5*(arr1(:,:,1    )
     .                      +arr1(:,:,0    ))
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
        vec2(nxp,:,:,:) = 0.5*(vec1(nxp  ,:,:,:)
     .                        +vec1(nxp-1,:,:,:))
        vec2(0  ,:,:,:) = 0.5*(vec1(1    ,:,:,:)
     .                        +vec1(0    ,:,:,:))

        vec2(:,nyp,:,:) = 0.5*(vec1(:,nyp  ,:,:)
     .                        +vec1(:,nyp-1,:,:))
        vec2(:,0  ,:,:) = 0.5*(vec1(:,1    ,:,:)
     .                        +vec1(:,0    ,:,:))

        vec2(:,:,nzp,:) = 0.5*(vec1(:,:,nzp  ,:)
     .                        +vec1(:,:,nzp-1,:))
        vec2(:,:,0  ,:) = 0.5*(vec1(:,:,1    ,:)
     .                        +vec1(:,:,0    ,:))
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
      subroutine setMGBC(gpos,neq,nnx,nny,nnz,iig,array,bcnd,arr_cov
     .                  ,arr0,icomp,is_cnv,is_vec,is_curv,result_is_curv
     .                  ,out_iscnv,out_iscurv,iorder,g_def)
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

      type(grid_mg_def),optional,pointer :: g_def

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

      logical      :: res_is_curv,fpointers,vector,isvec

      real(8)      :: car(3),curv(3)

      real(8),pointer,dimension(:,:,:,:) :: v_cov,v_cnv,vv0

      character(3) :: cvec

c     Begin program

      neql = neq !To avoid problems if neq is hardwired

c     Optional arguments

      if (PRESENT(icomp))  ivar(1:size(icomp))  = icomp

      if (PRESENT(iorder)) order = iorder

      if (PRESENT(g_def)) then
        gbc_def => g_def
      else
        gbc_def => grid_params
      endif

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

      call allocPointers(neql,fpointers,g_def=gbc_def)

c     Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(gpos,nnx,nny,nnz,iig,imnl,imxl,jmnl,jmxl,kmnl,kmxl)

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
        isvec = is_vec
      else
        isvec = (neql == 3 .or. neql == 6) !Backward compatibility with previous
                                           ! versions of this routine
      endif

      vector = (neql == 3 .or. neql == 6).and.isvec

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
            call transformVector(iig,iig,iig
     .           ,imnl-1,imxl+1,jmnl-1,jmxl+1,kmnl-1,kmxl+1
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,1+3*(ieq-1))
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,2+3*(ieq-1))
     .     ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,3+3*(ieq-1))
     .                            ,'car',cvec)
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
     .           = XformVector(iig,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .                        ,.not.iscnv)
          enddo
        else
          do ieq=1,neql/3
            v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .           = XformVector(iig,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .                        ,.not.iscnv)
          enddo
        endif

        do ieq=1,neql/3
          call setBC(ivar(ieq),3,nnx,nny,nnz
     .            ,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .            ,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .            ,vv0  (:,:,:,1+3*(ieq-1):3*ieq)
     .            ,bcnd (:    ,1+3*(ieq-1):3*ieq)
     .            ,iig,iig,iig
     .            ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .            ,is_cnv=iscnv,is_vec=.true.,iorder=order
     $            ,g_def=g_def)
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
            call transformVector(iig,iig,iig
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
          call setBC(ivar(1),nnx,nny,nnz
     $              ,array(:,:,:,ieq),vv0(:,:,:,ieq)
     .              ,bcnd(:,ieq),iig,iig,iig
     .              ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .              ,iorder=order,g_def=g_def)
        enddo

        if (.not.PRESENT(arr0)) deallocate(vv0)

      endif

c     Deallocate MG pointers (if allocated in this routine)

      call deallocPointers(fpointers)

      nullify(gbc_def)

c     End

      end subroutine setMGBC

c     setMGBC_stg
c     ###############################################################
      subroutine setMGBC_stg(gpos,neq,nnx,nny,nnz,iig,array,bcnd,arr_cov
     .                  ,arr0,icomp,is_cnv,is_vec,is_curv,result_is_curv
     .                  ,out_iscnv,out_iscurv,iorder,g_def)
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

      type(grid_mg_def),optional,pointer :: g_def

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

      logical      :: res_is_curv,fpointers,vector,isvec

      real(8)      :: car(3),curv(3)

      real(8),pointer,dimension(:,:,:,:) :: v_cov,v_cnv,vv0

      character(3) :: cvec

c     Begin program

      neql = neq !To avoid problems if neq is hardwired

      if (PRESENT(is_vec)) then
        isvec = is_vec
      else
        isvec = (neql == 3 .or. neql == 6) !Backward compatibility with previous
                                           ! versions of this routine
      endif

      vector = (neql == 3 .or. neql == 6).and.isvec

      if (.not.vector) then
        call pstop('setMGBC_stg','Input is not a vector')
      endif

c     Optional arguments

      if (PRESENT(icomp))  ivar(1:size(icomp))  = icomp

      if (PRESENT(iorder)) order = iorder

      if (PRESENT(g_def)) then
        gbc_def => g_def
      else
        gbc_def => grid_params
      endif

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

      call allocPointers(neql,fpointers,g_def=gbc_def)

c     Find grid node LOCAL position (if gpos > 0, else return local domain limits)

      call limits(gpos,nnx,nny,nnz,iig,imnl,imxl,jmnl,jmxl,kmnl,kmxl)

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

        call pstop('setMGBC_stg'
     $       ,'Staggered transformation from Cartesian not implemented')

c$$$        do ieq=1,neql/3
c$$$          call transformVector(iig,iig,iig
c$$$     .         ,imnl-1,imxl+1,jmnl-1,jmxl+1,kmnl-1,kmxl+1
c$$$     .   ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,1+3*(ieq-1))
c$$$     .   ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,2+3*(ieq-1))
c$$$     .   ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,3+3*(ieq-1))
c$$$     .                          ,'car',cvec)
c$$$        enddo

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
     .         = XformToCov_stg_mesh(iig,v_cnv(:,:,:,1+3*(ieq-1):3*ieq))
        enddo
      else
        do ieq=1,neql/3
          v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .         = XformToCnv_stg_mesh(iig,v_cov(:,:,:,1+3*(ieq-1):3*ieq))
        enddo
      endif

      do ieq=1,neql/3
        call imposeBConVector_stg(ivar(ieq),3,nnx,nny,nnz
     .          ,v_cnv(:,:,:,1+3*(ieq-1):3*ieq)
     .          ,v_cov(:,:,:,1+3*(ieq-1):3*ieq)
     .          ,vv0  (:,:,:,1+3*(ieq-1):3*ieq)
     .          ,bcnd (:    ,1+3*(ieq-1):3*ieq)
     .          ,iig,iig,iig
     .          ,i1=imng,i2=imxg,j1=jmng,j2=jmxg,k1=kmng,k2=kmxg
     .          ,is_cnv=iscnv,is_vec=.true.,iorder=order
     $          ,g_def=g_def)
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

        call pstop('setMGBC_stg'
     $       ,'Staggered transformation to Cartesian not implemented')

c$$$        do ieq=1,neql/3
c$$$          call transformVector(iig,iig,iig
c$$$     .         ,imnl-1,imxl+1,jmnl-1,jmxl+1,kmnl-1,kmxl+1
c$$$     .   ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,1+3*(ieq-1))
c$$$     .   ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,2+3*(ieq-1))
c$$$     .   ,array(imnl-1:imxl+1,jmnl-1:jmxl+1,kmnl-1:kmxl+1,3+3*(ieq-1))
c$$$     .                        ,cvec,'car')
c$$$        enddo

      endif

      !For output
      if (PRESENT(out_iscnv )) out_iscnv =iscnv
      if (PRESENT(out_iscurv)) out_iscurv=iscurv

c     Deallocate MG pointers (if allocated in this routine)

      call deallocPointers(fpointers)

      nullify(gbc_def)

c     End

      end subroutine setMGBC_stg

      end module setMGBC_interface

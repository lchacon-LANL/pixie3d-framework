c module parameters
c #####################################################################
      module parameters

      integer    :: neqd
     .             ,nxd,nyd,nzd,nxdp,nydp,nzdp
     .             ,nxl,nyl,nzl,nxlp,nylp,nzlp

      integer    :: ilo ,ihi ,jlo ,jhi ,klo ,khi
     .             ,ilom,ihip,jlom,jhip,klom,khip
     $             ,ilog ,ihig ,jlog ,jhig ,klog ,khig
     .             ,ilomg,ihipg,jlomg,jhipg,klomg,khipg
     .             ,gcw

      integer    :: ntotd,ntotdp,ntimemax

      end module parameters

c module variable_setup
c #####################################################################
      module variable_setup

        use parameters

        use grid

        implicit none
        
        integer :: NAUXS=0,NAUXV=0  !Number of auxiliary scalars, vectors
                                    !(defined by application)

        type :: var_def
          integer :: bconds(6)          !Boundary conditions
          real(8),pointer,dimension(:,:,:) :: array ! contains values at grid points
          character(20) :: descr
        end type var_def

        type :: vec_def
          integer :: bconds(6,3)        !Boundary conditions, 3 refers to the 3 components
          real(8),pointer,dimension(:,:,:,:) :: vec ! vector (i,j,k,d), d being the component
          character(20) :: descr        !Description
          logical       :: cnv          !Whether it is contravariant or not
        end type vec_def

        type :: var_array
          integer :: nvar               !Number of variables
          type (var_def),pointer,dimension(:) :: array_var ! array of variables
        end type var_array

        type :: aux_array
          integer :: nvar,nvec          !Number of variables, vectors
          type (vec_def),pointer,dimension(:) :: vec_list   ! vector auxillary components
          type (var_def),pointer,dimension(:) :: var_list   ! scalar auxillary components
        end type aux_array

        type :: bcgrp
          integer  :: nbc_seq                         ! number of variables in bc sequence
          integer,pointer,dimension(:,:):: bc_seq     ! array showing sequence in which to set bc's 
        end type bcgrp

        type :: patch
          integer  :: nbc_grp                         ! number of BC bundled groups
          type(bcgrp),pointer,dimension(:):: bc_grp   ! array of bundled variables for BC treatment

          type(var_array)  ,pointer     :: u_0        ! equilibrium data
          type(var_array)  ,pointer     :: u_n        ! data at time n
          type(var_array)  ,pointer     :: u_np       ! data at n+1
          type(var_array)  ,pointer     :: u_nm       ! data at n-1
          type(var_array)  ,pointer     :: u_graph    ! auxiliary for plotting, all u variables same number of components
          type(grid_mg_def),pointer     :: gparams    ! contains all grid info
          type(aux_array)  ,pointer     :: aux        ! auxillary variables
        end type patch

        type (patch),pointer :: gv => null()

        type(var_array),pointer :: u_0 => null()

      contains

c     setVectorDimensions
c     #################################################################
      subroutine setVectorDimensions

        !Ghost cell width
        gcw = gv%gparams%gcw

        !Global problem info
        nxdp = gv%gparams%nglx + gcw
        nydp = gv%gparams%ngly + gcw
        nzdp = gv%gparams%nglz + gcw

        !Define offset quantities (global grid)
        ihig = gv%gparams%ihig
        ilog = gv%gparams%ilog
        jhig = gv%gparams%jhig
        jlog = gv%gparams%jlog
        khig = gv%gparams%khig
        klog = gv%gparams%klog

        ihipg = gv%gparams%ihig + gcw
        ilomg = gv%gparams%ilog - gcw
        jhipg = gv%gparams%jhig + gcw
        jlomg = gv%gparams%jlog - gcw
        khipg = gv%gparams%khig + gcw
        klomg = gv%gparams%klog - gcw

        !Local problem info
        nxl  = gv%gparams%nlx
        nyl  = gv%gparams%nly
        nzl  = gv%gparams%nlz

        !Define local limits
        ilo  = 1
        ihi  = nxl
        jlo  = 1
        jhi  = nyl
        klo  = 1
        khi  = nzl

        !Define offset quantities (local grid)
        ihip = ihi + gcw
        ilom = ilo - gcw
        jhip = jhi + gcw
        jlom = jlo - gcw
        khip = khi + gcw
        klom = klo - gcw

        nxlp = nxl + gcw
        nylp = nyl + gcw
        nzlp = nzl + gcw

        !Define local vector dimensions
        ntotdp = nxl*nyl*nzl
        ntotd  = neqd*ntotdp

      end subroutine setVectorDimensions

c     allocateGlobalVar
c     #################################################################
      subroutine allocateGlobalVar(gv)

        implicit none

c     Call variables

        type(patch),pointer  :: gv

c     Local variables

        integer          :: ieq

c     Begin program

        if (.not.associated(gv)) then
          allocate(gv)
          nullify(gv%bc_grp)
          nullify(gv%u_np,gv%u_n,gv%u_nm,gv%u_0,gv%u_graph)
          nullify(gv%gparams)
          nullify(gv%aux)
        endif

c     End program

      end subroutine allocateGlobalVar

c     allocateAuxStruct
c     #################################################################
      subroutine allocateAuxStruct(nvar,nvec,aux)

        implicit none

c     Call variables

        integer :: nvar,nvec
        type(aux_array),pointer  :: aux

c     Local variables

        integer          :: ieq

c     Begin program

        if (.not.associated(aux)) then
          allocate(aux)

          aux%nvar = nvar
          aux%nvec = nvec

          if (nvar > 0) then
            allocate(aux%var_list(nvar))

            do ieq=1,nvar
              allocate(aux%var_list(ieq)
     $                     %array(ilom:ihip,jlom:jhip,klom:khip))
              aux%var_list(ieq)%array = 0d0
              aux%var_list(ieq)%bconds = bcond
            enddo
          endif

          if (nvec > 0) then
            allocate(aux%vec_list(nvec))

            do ieq=1,nvec
              allocate(aux%vec_list(ieq)
     $                    %vec(ilom:ihip,jlom:jhip,klom:khip,3))
              aux%vec_list(ieq)%vec = 0d0
              aux%vec_list(ieq)%bconds(:,1) = bcond
              aux%vec_list(ieq)%bconds(:,2) = bcond
              aux%vec_list(ieq)%bconds(:,3) = bcond
            enddo
          endif
        endif

c     End program

      end subroutine allocateAuxStruct

c     deallocateAuxStruct
c     #################################################################
      subroutine deallocateAuxStruct(aux)

        implicit none

c     Call variables

        type(aux_array),pointer  :: aux

c     Local variables

        integer :: ieq,nvar,nvec

c     Begin program

        if (associated(aux)) then

          nvar = aux%nvar
          nvec = aux%nvec

          if (nvar > 0) then
            do ieq=1,nvar
              deallocate(aux%var_list(ieq)%array)
            enddo

            deallocate(aux%var_list)
          endif

          if (nvec > 0) then
            do ieq=1,nvec
              deallocate(aux%vec_list(ieq)%vec)
            enddo

            deallocate(aux%vec_list)
          endif

          deallocate(aux)

        endif

c     End program

      end subroutine deallocateAuxStruct

c     allocateDerivedType
c     #################################################################
      subroutine allocateDerivedType(varray,igx,igy,igz)

        implicit none

c     Call variables

        type(var_array),pointer  :: varray
        integer, optional :: igx,igy,igz

c     Local variables

        integer          :: ieq,imn,imx,jmn,jmx,kmn,kmx

c     Begin program

        if (.not.associated(varray)) then
          allocate(varray)

          varray%nvar = neqd

          if (PRESENT(igx)) then
            imn = 1-gcw
            imx = gv%gparams%nxv(igx)+gcw
          else
            imn = ilom
            imx = ihip
          endif

          if (PRESENT(igy)) then
            jmn = 1-gcw
            jmx = gv%gparams%nyv(igy)+gcw
          else
            jmn = jlom
            jmx = jhip
          endif

          if (PRESENT(igz)) then
            kmn = 1-gcw
            kmx = gv%gparams%nzv(igz)+gcw
          else
            kmn = klom
            kmx = khip
          endif

          allocate(varray%array_var(varray%nvar))
          do ieq=1,varray%nvar
            allocate(varray%array_var(ieq)
     $                     %array(imn:imx,jmn:jmx,kmn:kmx))
cc     $                     %array(ilom:ihip,jlom:jhip,klom:khip))
            varray%array_var(ieq)%array = 0d0
          enddo
        endif

c     End program

      end subroutine allocateDerivedType

c     deallocateDerivedType
c     #################################################################
      subroutine deallocateDerivedType(varray)

        implicit none

c     Call variables

        type(var_array),pointer  :: varray

c     Local variables

        integer          :: ieq,ierr

c     Begin program

cc        deallocate(varray)

        if (associated(varray)) then
          if (associated(varray%array_var)) then
            do ieq=1,neqd
              if (associated(varray%array_var(ieq)%array)) then
                deallocate(varray%array_var(ieq)%array)
              endif
            enddo
            deallocate(varray%array_var)
          endif
          deallocate(varray)
cc        nullify(varray)
        endif

c     End program

      end subroutine deallocateDerivedType

c     equateDerivedType
c     #################################################################
      subroutine equateDerivedType(varray2,varray1)
      implicit none
c     -----------------------------------------------------------------
c     Performs varray2 = varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),pointer :: varray1
        type(var_array),pointer :: varray2

c     Local variables

        integer    :: ieq

c     Begin program

        if (associated(varray2)) then
           if (varray2%nvar /= neqd) nullify(varray2)
        endif

        call allocateDerivedType(varray2)

        varray2%nvar = varray1%nvar
        do ieq=1,varray2%nvar
          varray2%array_var(ieq)%bconds = varray1%array_var(ieq)%bconds
          varray2%array_var(ieq)%array  = varray1%array_var(ieq)%array
          varray2%array_var(ieq)%descr  = varray1%array_var(ieq)%descr
        enddo

c     End program

      end subroutine equateDerivedType

c     initializeDerivedType
c     #################################################################
      subroutine initializeDerivedType(varray,igx,igy,igz)

c     -----------------------------------------------------------------
c     Initializes fields of varray2 using a template (u_0) except for
c     numeric arrays.
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),pointer  :: varray
        integer, optional :: igx,igy,igz

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray,igx=igx,igy=igy,igz=igz)

        varray%nvar = u_0%nvar
        do ieq=1,varray%nvar
          varray%array_var(ieq)%bconds = u_0%array_var(ieq)%bconds
          varray%array_var(ieq)%descr  = u_0%array_var(ieq)%descr
          varray%array_var(ieq)%array  = 0d0
        enddo

c     End program

      end subroutine initializeDerivedType

c     AXPYDerivedType
c     #################################################################
      subroutine AXPYDerivedType(alpha2,varray2,alpha1,varray1,varray3)

c     -----------------------------------------------------------------
c     Performs varray3 = alpha2*varray2 + alpha1*varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),pointer :: varray3
        type(var_array),pointer :: varray1,varray2
        real(8) :: alpha1,alpha2

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray3)

        varray3%nvar = varray2%nvar
        do ieq=1,varray3%nvar
          varray3%array_var(ieq)%bconds = varray2%array_var(ieq)%bconds
          varray3%array_var(ieq)%array  =
     .          alpha2*varray2%array_var(ieq)%array
     .         +alpha1*varray1%array_var(ieq)%array
          varray3%array_var(ieq)%descr  = varray2%array_var(ieq)%descr
        enddo

c     End program

      end subroutine AXPYDerivedType

c     varPack
c     #################################################################
      subroutine varPack(array,bcs,desc,veq,varray)

c     -----------------------------------------------------------------
c     Fills structure fields of varray with pertinent data corresponding
c     to variable veq: boundary conditions, bcs; variable description,
c     desc; numerical data, array).
c     -----------------------------------------------------------------

c     Call variables

        integer      :: veq,bcs(6)
        real(8)      :: array(ilom:ihip,jlom:jhip,klom:khip)
        character(*) :: desc

        type(var_array),pointer :: varray

c     Local variables

        call allocateDerivedType(varray)

c     Begin program

        varray%array_var(veq)%bconds = bcs
        varray%array_var(veq)%array  = array
        varray%array_var(veq)%descr  = desc

c     End program

      end subroutine varPack

c     varUnPack
c     #################################################################
      subroutine varUnPack(varray,veq,array,bcs,desc)

c     Call variables

        integer       :: nx,ny,nz,veq,bcs(6)
        real(8)       :: array(ilom:ihip,jlom:jhip,klom:khip)
        character(20) :: desc

        type(var_array),pointer :: varray

c     Local variables

c     Begin program

        bcs   = varray%array_var(veq)%bconds
        array = varray%array_var(veq)%array
        desc  = trim(varray%array_var(veq)%descr)

c     End program

      end subroutine varUnPack

c     mapStructureToVector
c     #################################################################
      subroutine mapStructureToVector(x,varray,igx,igy,igz)

c     -----------------------------------------------------------------
c     Maps structure into vector (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(out)      :: x(:)
      integer, optional :: igx,igy,igz

      type(var_array),pointer  :: varray

c     Local variables

      integer :: i,j,k,ii,ieq,neq,iigx=1,iigy=1,iigz=1

c     Begin program

      if (PRESENT(igx)) iigx=igx
      if (PRESENT(igy)) iigy=igy
      if (PRESENT(igz)) iigz=igz

c     Unpack vector x

      neq = varray%nvar

      do k = 1,gv%gparams%nzv(iigz)
        do j = 1,gv%gparams%nyv(iigy)
          do i = 1,gv%gparams%nxv(iigx)
            ii = vecPos(neq,i,j,k,iigx,iigy,iigz)
            do ieq=1,neq
              x(ii+ieq) = varray%array_var(ieq)%array(i,j,k)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapStructureToVector

c     mapVectorToStructure
c     #################################################################
      subroutine mapVectorToStructure(varray,x,igx,igy,igz)

c     -----------------------------------------------------------------
c     Maps vector solution into structure.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

cc      real(8),intent(in)      :: x(ntotd)
      real(8),dimension(:),intent(in) :: x(:)
      integer, optional :: igx,igy,igz

      type(var_array),pointer :: varray

c     Local variables

      integer    :: i,j,k,ii,ieq,neq,iigx=1,iigy=1,iigz=1

c     Begin program

      if (PRESENT(igx)) iigx=igx
      if (PRESENT(igy)) iigy=igy
      if (PRESENT(igz)) iigz=igz

c     Initialize varray

      call initializeDerivedType(varray)

c     Unpack vector x

      neq = varray%nvar

      do k = 1,gv%gparams%nzv(iigz)
        do j = 1,gv%gparams%nyv(iigy)
          do i = 1,gv%gparams%nxv(iigx)
            ii = vecPos(neq,i,j,k,iigx,iigy,iigz)
            do ieq=1,neq
              varray%array_var(ieq)%array(i,j,k) = x(ii+ieq)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapVectorToStructure

      end module variable_setup

c module variables
c #####################################################################
      module variables

        use variable_setup

        type(var_array),pointer :: u_np => null()
     .                            ,u_n  => null()
     .                            ,u_nm => null()
     .                            ,u_ic => null()
     .                            ,u_graph => null()
  
      contains

c     allocateStructures
c     #################################################################
      subroutine allocateStructures

c     Call variables

c     Local variables

c     Begin program

        call allocateDerivedType(gv%u_np)
        call allocateDerivedType(gv%u_nm)
        call allocateDerivedType(gv%u_n )
        call allocateDerivedType(gv%u_0)
        call allocateDerivedType(gv%u_graph)

        u_0     => gv%u_0
        u_np    => gv%u_np
        u_n     => gv%u_n
        u_nm    => gv%u_nm
        u_graph => gv%u_graph

c     End program

      end subroutine allocateStructures

c     deallocateStructures
c     #################################################################
      subroutine deallocateStructures

c     Call variables

c     Local variables

c     Begin program

        call deallocateDerivedType(gv%u_np)
        call deallocateDerivedType(gv%u_nm)
        call deallocateDerivedType(gv%u_n )
        call deallocateDerivedType(gv%u_0)
        call deallocateDerivedType(gv%u_graph)

c     End program

      end subroutine deallocateStructures

      end module variables

ccc module generalPurposeFunctions
ccc ######################################################################
cc      module generalPurposeFunctions
cc
cc        use grid
cc
cc      contains
cc
ccc     integral
ccc     ################################################################
cc      function integral(nx,ny,nz,array,igx,igy,igz,avg)
cc
ccc     ---------------------------------------------------------------
ccc     Integrates array(i,j,k) on domain (nx)x(ny)x(nz).
ccc     ---------------------------------------------------------------
cc
cc      implicit none
cc
ccc    Call variables
cc
cc      integer :: igx,igy,igz,nx,ny,nz
cc      real(8) :: array(0:nx+1,0:ny+1,0:nz+1),integral
cc      logical :: avg
cc
ccc     Local variables
cc
cc      integer :: i,j,k
cc
cc      real(8) :: tvolume,vol,lvolume,lintegral
cc
ccc     Begin program
cc
ccc     Integrate
cc
cc      lintegral = 0d0
cc      lvolume   = 0d0
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            vol = gmetric%grid(igx)%dvol(i,j,k)
cc
cc            if (isSYM(i,igx,1,0)) vol = 0.5*vol
cc            if (isSYM(i,igx,1,1)) vol = 0.5*vol
cc            if (isSYM(j,igy,2,0)) vol = 0.5*vol
cc            if (isSYM(j,igy,2,1)) vol = 0.5*vol
cc            if (isSYM(k,igz,3,0)) vol = 0.5*vol
cc            if (isSYM(k,igz,3,1)) vol = 0.5*vol
cc
cc            lintegral = lintegral + array(i,j,k)*vol
cc            lvolume = lvolume + vol
cc          enddo
cc        enddo
cc      enddo
cc
cc#if defined(petsc)
cc      call MPI_Allreduce(lintegral,integral,1,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc      if (avg)
cc     .call MPI_Allreduce(lvolume  ,tvolume ,1,MPI_DOUBLE_PRECISION
cc     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc#else
cc      integral = lintegral
cc      tvolume  = lvolume
cc#endif
cc
cc      if (avg) integral = integral/tvolume
cc
ccc     End 
cc
cc      end function integral
cc
cc      end module generalPurposeFunctions

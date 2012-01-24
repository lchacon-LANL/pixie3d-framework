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
          integer,pointer,dimension(:,:):: bc_dep_list ! Array indicating inter-variable dependencies;
                                                       !   dimensioned (NAUXS+NAUXV+NEQD,3).
                                                       ! First dimension cycles over number of deps.
                                                       ! Second dimension provides information about coupling:
                                                       !  + First component identifies variable; its sign indicates
                                                       !      auxiliary (<0) or dependent (>0) variable.
                                                       !  + Second component indicates scalar (0) or vector (1)
                                                       !  + Third component indicates trivial (0) or non-trivial (1)
          integer,pointer,dimension(:,:):: dom_dep_list! Array indicating inter-variable dependencies;
                                                       !   dimensioned (NAUXS+NAUXV+NEQD,3).
                                                       ! First dimension cycles over number of deps.
                                                       ! Second dimension provides information about coupling:
                                                       !  + First component identifies variable; its sign indicates
                                                       !      auxiliary (<0) or dependent (>0) variable.
                                                       !  + Second component indicates scalar (0) or vector (1)
                                                       !  + Third component indicates trivial (0) or non-trivial (1)
          logical :: time_0             !Whether this is a time=0 variable
          logical :: time_n             !Whether this is a previous time step variable
        end type var_def

        type :: vec_def
          integer :: bconds(6,3)        !Boundary conditions, 3 refers to the 3 components
          real(8),pointer,dimension(:,:,:,:) :: vec ! vector (i,j,k,d), d being the component
          character(20) :: descr        !Description
          logical       :: cnv          !Whether it is contravariant or not
          integer,pointer,dimension(:,:):: bc_dep_list ! Array indicating inter-variable dependencies;
                                                       !   dimensioned (NAUXS+NAUXV+NEQD,3).
                                                       ! First dimension cycles over number of deps.
                                                       ! Second dimension provides information about coupling:
                                                       !  + First component identifies variable; its sign indicates
                                                       !      auxiliary (<0) or dependent (>0) variable.
                                                       !  + Second component indicates scalar (0) or vector (1)
                                                       !  + Third component indicates trivial (0) or non-trivial (1)
          integer,pointer,dimension(:,:):: dom_dep_list! Array indicating inter-variable dependencies;
                                                       !   dimensioned (NAUXS+NAUXV+NEQD,3).
                                                       ! First dimension cycles over number of deps.
                                                       ! Second dimension provides information about coupling:
                                                       !  + First component identifies variable; its sign indicates
                                                       !      auxiliary (<0) or dependent (>0) variable.
                                                       !  + Second component indicates scalar (0) or vector (1)
                                                       !  + Third component indicates trivial (0) or non-trivial (1)
          logical :: time_0             !Whether this is a time=0 variable
          logical :: time_n             !Whether this is a previous time step variable
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

          type(var_array)  ,pointer :: u_0        ! equilibrium data
          type(var_array)  ,pointer :: u_n        ! data at time n
          type(var_array)  ,pointer :: u_np       ! data at n+1
          type(var_array)  ,pointer :: u_nm       ! data at n-1

          type(aux_array)  ,pointer :: aux        ! auxiliary variables
          type(grid_mg_def),pointer :: gparams    ! contains all grid info
        end type patch

        type (patch),pointer :: gv => null()

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
          nullify(gv%u_np,gv%u_n,gv%u_nm,gv%u_0)
          nullify(gv%gparams)
          nullify(gv%aux)
        endif

c     End program

      end subroutine allocateGlobalVar

c     deallocateGlobalVar
c     #################################################################
      subroutine deallocateGlobalVar(gv)

        implicit none

c     Call variables

        type(patch),pointer  :: gv

c     Local variables

        integer :: ierr,igr

c     Begin program

        call deallocateAuxStruct(gv%aux)

        call destroyGrid(gv%gparams)

        if (associated(gv%bc_grp)) then
          do igr=1,size(gv%bc_grp)
            deallocate(gv%bc_grp(igr)%bc_seq,STAT=ierr)
          enddo
          deallocate(gv%bc_grp,STAT=ierr)
        endif

        deallocate(gv,STAT=ierr)

c     End program

      end subroutine deallocateGlobalVar

c     allocateAuxStruct
c     #################################################################
      subroutine allocateAuxStruct(aux)

        implicit none

c     Call variables

        type(aux_array),pointer  :: aux

c     Local variables

        integer :: ieq,nvar,nvec

c     Begin program

        nvar = NAUXS
        nvec = NAUXV

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
              aux%var_list(ieq)%bconds= bcond
              aux%var_list(ieq)%descr = 'unset'

              allocate(aux%var_list(ieq)%bc_dep_list(nvar+nvec+neqd,3))
              aux%var_list(ieq)%bc_dep_list = 0
              allocate(aux%var_list(ieq)%dom_dep_list(nvar+nvec+neqd,3))
              aux%var_list(ieq)%dom_dep_list = 0
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
              aux%vec_list(ieq)%descr       = 'unset'

              allocate(aux%vec_list(ieq)%bc_dep_list(nvar+nvec+neqd,3))
              aux%vec_list(ieq)%bc_dep_list = 0
              allocate(aux%vec_list(ieq)%dom_dep_list(nvar+nvec+neqd,3))
              aux%vec_list(ieq)%dom_dep_list = 0
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
              deallocate(aux%var_list(ieq)%bc_dep_list)
              deallocate(aux%var_list(ieq)%dom_dep_list)
            enddo

            deallocate(aux%var_list)
          endif

          if (nvec > 0) then
            do ieq=1,nvec
              deallocate(aux%vec_list(ieq)%vec)
              deallocate(aux%vec_list(ieq)%bc_dep_list)
              deallocate(aux%vec_list(ieq)%dom_dep_list)
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

        integer :: ieq,imn,imx,jmn,jmx,kmn,kmx

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
            varray%array_var(ieq)%array = 0d0
            varray%array_var(ieq)%bconds= bcond
            varray%array_var(ieq)%descr = 'unset'

            allocate(varray%array_var(ieq)
     .                     %bc_dep_list(NAUXS+NAUXV+neqd,3))
            varray%array_var(ieq)%bc_dep_list = 0
            allocate(varray%array_var(ieq)
     .                     %dom_dep_list(NAUXS+NAUXV+neqd,3))
            varray%array_var(ieq)%dom_dep_list = 0
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

        integer          :: ieq,ierr,istat

c     Begin program

        if (associated(varray)) then
          do ieq=1,size(varray%array_var)
            deallocate(varray%array_var(ieq)%array,stat=istat)
            deallocate(varray%array_var(ieq)%bc_dep_list,stat=istat)
            deallocate(varray%array_var(ieq)%dom_dep_list,stat=istat)
          enddo

          deallocate(varray%array_var,stat=istat)

          deallocate(varray,stat=istat)
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

cc        if (associated(varray2)) then
cccc           if (varray2%nvar /= neqd) nullify(varray2)
cc           if (varray2%nvar /= neqd) call deallocateDerivedType(varray2)
cc        endif

        call allocateDerivedType(varray2)

        varray2%nvar = varray1%nvar

        do ieq=1,varray2%nvar
          varray2%array_var(ieq)%bconds =varray1%array_var(ieq)%bconds
          varray2%array_var(ieq)%array  =varray1%array_var(ieq)%array
          varray2%array_var(ieq)%descr  =varray1%array_var(ieq)%descr

          varray2%array_var(ieq)%bc_dep_list
     .                              =varray1%array_var(ieq)%bc_dep_list
          varray2%array_var(ieq)%dom_dep_list
     .                              =varray1%array_var(ieq)%dom_dep_list
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

        varray%nvar = gv%u_0%nvar

        do ieq=1,varray%nvar
          varray%array_var(ieq)%bconds = gv%u_0%array_var(ieq)%bconds
          varray%array_var(ieq)%descr  = gv%u_0%array_var(ieq)%descr
          varray%array_var(ieq)%bc_dep_list
     .                              = gv%u_0%array_var(ieq)%bc_dep_list
          varray%array_var(ieq)%dom_dep_list
     .                              = gv%u_0%array_var(ieq)%dom_dep_list
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
          varray3%array_var(ieq)%bconds =varray2%array_var(ieq)%bconds
          varray3%array_var(ieq)%array  =
     .          alpha2*varray2%array_var(ieq)%array
     .         +alpha1*varray1%array_var(ieq)%array
          varray3%array_var(ieq)%descr  =varray2%array_var(ieq)%descr

          varray3%array_var(ieq)%bc_dep_list
     .                              =varray2%array_var(ieq)%bc_dep_list
          varray3%array_var(ieq)%dom_dep_list
     .                              =varray2%array_var(ieq)%dom_dep_list
        enddo

c     End program

      end subroutine AXPYDerivedType

ccc     varPack
ccc     #################################################################
cc      subroutine varPack(array,bcs,desc,veq,varray)
cc
ccc     -----------------------------------------------------------------
ccc     Fills structure fields of varray with pertinent data corresponding
ccc     to variable veq: boundary conditions, bcs; variable description,
ccc     desc; numerical data, array).
ccc     -----------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer      :: veq,bcs(6)
cc        real(8)      :: array(ilom:ihip,jlom:jhip,klom:khip)
cc        character(*) :: desc
cc
cc        type(var_array),pointer :: varray
cc
ccc     Local variables
cc
cc        call allocateDerivedType(varray)
cc
ccc     Begin program
cc
cc        varray%array_var(veq)%bconds = bcs
cc        varray%array_var(veq)%array  = array
cc        varray%array_var(veq)%descr  = desc
cc
ccc     End program
cc
cc      end subroutine varPack
cc
ccc     varUnPack
ccc     #################################################################
cc      subroutine varUnPack(varray,veq,array,bcs,desc)
cc
ccc     Call variables
cc
cc        integer       :: nx,ny,nz,veq,bcs(6)
cc        real(8)       :: array(ilom:ihip,jlom:jhip,klom:khip)
cc        character(20) :: desc
cc
cc        type(var_array),pointer :: varray
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        bcs   = varray%array_var(veq)%bconds
cc        array = varray%array_var(veq)%array
cc        desc  = trim(varray%array_var(veq)%descr)
cc
ccc     End program
cc
cc      end subroutine varUnPack

c     mapStructureToVector
c     #################################################################
      subroutine mapStructureToVector(x,varray,igr)

c     -----------------------------------------------------------------
c     Maps structure into vector (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(out)      :: x(:)
      integer, optional :: igr

      type(var_array),pointer  :: varray

c     Local variables

      integer :: i,j,k,ii,ieq,neq,iigr=1

c     Begin program

      if (PRESENT(igr)) iigr=igr

c     Unpack vector x

      neq = varray%nvar

      do k = 1,gv%gparams%nzv(iigr)
        do j = 1,gv%gparams%nyv(iigr)
          do i = 1,gv%gparams%nxv(iigr)
            ii = vecPos(neq,i,j,k,iigr,iigr,iigr)
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
      subroutine mapVectorToStructure(varray,x,igr)

c     -----------------------------------------------------------------
c     Maps vector solution into structure.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),dimension(:),intent(in) :: x(:)
      integer, optional :: igr

      type(var_array),pointer :: varray

c     Local variables

      integer    :: i,j,k,ii,ieq,neq,iigr=1

c     Begin program

      if (PRESENT(igr)) iigr=igr

c     Initialize varray

      call initializeDerivedType(varray)

c     Unpack vector x

      neq = varray%nvar

      do k = 1,gv%gparams%nzv(iigr)
        do j = 1,gv%gparams%nyv(iigr)
          do i = 1,gv%gparams%nxv(iigr)
            ii = vecPos(neq,i,j,k,iigr,iigr,iigr)
            do ieq=1,neq
              varray%array_var(ieq)%array(i,j,k) = x(ii+ieq)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapVectorToStructure

c     gatherVector
c     ###################################################################
      subroutine gatherVector(xout,xin)

c     -------------------------------------------------------------------
c     Collects PC solution vectors into a single vector.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN)  :: xin(:,:)
        real(8),intent(OUT) :: xout(:)

c     Local variables

        integer    :: i,j,k,ii,iii,ieq

c     Begin program

        do ieq=1,size(xin,2)
          do ii=1,size(xin,1)
            iii = ieq + size(xin,2)*(ii-1)
            xout(iii) = xin(ii,ieq)
          enddo
        enddo

c     End program

      end subroutine gatherVector

c     scatterVector
c     ###################################################################
      subroutine scatterVector(xout,xin)

c     -------------------------------------------------------------------
c     Scatters PC residual into separate residuals, one per equation
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN)  :: xin(:)
        real(8),intent(OUT) :: xout(:,:)

c     Local variables

        integer :: i,j,k,ii,iii,ieq
        real(8) :: dvol

c     Begin program

        do ieq=1,size(xout,2)
          do ii=1,size(xout,1)
            iii = ieq + size(xout,2)*(ii-1)
            xout(ii,ieq) = xin(iii)
          enddo
        enddo

c     End program

      end subroutine scatterVector

      end module variable_setup

c module variables
c #####################################################################
      module variables

        use variable_setup

        type(var_array),pointer :: u_np => null()
     .                            ,u_n  => null()
     .                            ,u_nm => null()
     .                            ,u_ic => null()
     .                            ,u_0  => null()

      contains

c     allocateStructures
c     #################################################################
      subroutine allocateStructures

c     Call variables

c     Local variables

c     Begin program

        call allocateDerivedType(gv%u_n )
        call allocateDerivedType(gv%u_0)

        u_0     => gv%u_0
        u_n     => gv%u_n

#if !defined(samrai)
        call allocateDerivedType(gv%u_np)
        call allocateDerivedType(gv%u_nm)

        u_np    => gv%u_np
        u_nm    => gv%u_nm
#endif

c     End program

      end subroutine allocateStructures

c     deallocateStructures
c     #################################################################
      subroutine deallocateStructures

c     Call variables

c     Local variables

c     Begin program

        call deallocateDerivedType(gv%u_n )
        call deallocateDerivedType(gv%u_0)

#if !defined(samrai)
        call deallocateDerivedType(gv%u_np)
        call deallocateDerivedType(gv%u_nm)
#endif

c     End program

      end subroutine deallocateStructures

      end module variables

c module constants
c ######################################################################
      module constants

        real(8) :: pi
        real(8), allocatable, dimension(:,:,:)   :: zeros,ones
        real(8), allocatable, dimension(:,:,:,:) :: vzeros

      end module constants

c module parameters
c #####################################################################
      module parameters

      integer(4) :: neqd
     .             ,nxd,nyd,nzd,nxdp,nydp,nzdp
     .             ,nxl,nyl,nzl,nxlp,nylp,nzlp

      integer(4) :: ilo ,ihi ,jlo ,jhi ,klo ,khi
     .             ,ilom,ihip,jlom,jhip,klom,khip
     $             ,ilog ,ihig ,jlog ,jhig ,klog ,khig
     .             ,ilomg,ihipg,jlomg,jhipg,klomg,khipg

      integer(4) :: ntotd,ntotdp,ntimemax

      contains

c     setVectorDimensions
c     #################################################################
      subroutine setVectorDimensions

      !Global problem info
        nxdp = nxd+1
        nydp = nyd+1
        nzdp = nzd+1

cc        ntotdp  = nxd*nyd*nzd
cc        ntotd   = neqd*ntotdp

        !Define offset quantities (global grid)
        ihipg = ihig+1
        ilomg = ilog-1
        jhipg = jhig+1
        jlomg = jlog-1
        khipg = khig+1
        klomg = klog-1

      !Local problem info (ihig,ilog,etc. initialized elsewhere)

        nxl  = ihig-ilog+1
        nyl  = jhig-jlog+1
        nzl  = khig-klog+1

        !Define local limits
        ilo  = 1
        ihi  = nxl
        jlo  = 1
        jhi  = nyl
        klo  = 1
        khi  = nzl

        !Define offset quantities (local grid)
        ihip = ihi+1
        ilom = ilo-1
        jhip = jhi+1
        jlom = jlo-1
        khip = khi+1
        klom = klo-1

        nxlp = nxl+1
        nylp = nyl+1
        nzlp = nzl+1

        !Define local vector dimensions
        ntotdp = nxl*nyl*nzl
        ntotd  = neqd*ntotdp

      end subroutine setVectorDimensions

c$$$c     vecPos
c$$$c     #################################################################
c$$$      function vecPos(neq,i,j,k,igx,igy,igz) result (ijk)
c$$$
c$$$        implicit none
c$$$
c$$$        integer(4) :: ijk,i,j,k,neq,igx,igy,igz
c$$$
c$$$        ijk = neq*(i-ilo + nxl*(j-jlo) + nxl*nyl*(k-klo))
c$$$      
c$$$      end function vecPos
c$$$
c$$$c     fromGlobalToLocalLimits
c$$$c     #################################################################
c$$$      subroutine fromGlobalToLocalLimits(ig,jg,kg,il,jl,kl,igx,igy,igz)
c$$$
c$$$        implicit none
c$$$
c$$$        integer(4) :: ig,jg,kg,il,jl,kl,igx,igy,igz
c$$$
c$$$        il = ig - ilog + 1
c$$$        jl = jg - jlog + 1
c$$$        kl = kg - klog + 1
c$$$
c$$$      end subroutine fromGlobalToLocalLimits

      end module parameters

#if defined(petsc)
c module petsc_var_def
c #####################################################################
      module petsc_var_def

        use parameters

        implicit none

        !NVAR is defined via C preprocessor flag
        type :: petsc_var
          real(8) :: var(NVAR)
        end type petsc_var

c$$$        type :: petsc_array
c$$$          type(petsc_var),pointer,dimension(:,:,:) :: array
c$$$        end type petsc_array
c$$$
c$$$      contains
c$$$
c$$$c     allocatePetcsType
c$$$c     #################################################################
c$$$      subroutine allocatePetscType(petscarray)
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(petsc_array)  :: petscarray
c$$$
c$$$c     Local variables
c$$$
c$$$        integer(4)         :: ieq
c$$$
c$$$c     Begin program
c$$$
c$$$        !Error check
c$$$        if (neqd /= NVAR) then
c$$$          write (*,*) 'Different number of equations in Petsc and F90'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        endif
c$$$
c$$$        if (.not.associated(petscarray%array)) then
c$$$          allocate(petscarray%array(ilom:ihip,jlom:jhip,klom:khip))
c$$$          do ieq=1,neqd
c$$$            petscarray%array(:,:,:)%var(ieq) = 0d0
c$$$          enddo
c$$$        endif
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine allocatePetscType
c$$$
c$$$c     deallocatePetscType
c$$$c     #################################################################
c$$$      subroutine deallocatePetscType(petscarray)
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(petsc_array)  :: petscarray
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$        if (associated(petscarray%array)) deallocate(petscarray%array)
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine deallocatePetscType
c$$$
c$$$c     equatePetscType
c$$$c     #################################################################
c$$$      subroutine equatePetscType(petscarray2,petscarray1)
c$$$      implicit none
c$$$c     -----------------------------------------------------------------
c$$$c     Performs petscarray2 = petscarray1
c$$$c     -----------------------------------------------------------------
c$$$
c$$$c     Call variables
c$$$
c$$$        type(petsc_array),intent(in)  :: petscarray1
c$$$        type(petsc_array),intent(out) :: petscarray2
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$        call allocatePetscType(petscarray2)
c$$$
c$$$        petscarray2%array = petscarray1%array
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine equatePetscType
c$$$
c$$$c     initializePetscType
c$$$c     #################################################################
c$$$      subroutine initializePetscType(petscarray)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Initializes fields of petscarray.
c$$$c     -----------------------------------------------------------------
c$$$
c$$$c     Call variables
c$$$
c$$$        type(petsc_array)  :: petscarray
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$        call allocatePetscType(petscarray)
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine initializePetscType
c$$$
c$$$c     substractPetscType
c$$$c     #################################################################
c$$$      function substractPetscType(petscarray2,petscarray1)
c$$$     .         result (petscarray3)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Performs petscarray3 = petscarray2 - petscarray1
c$$$c     -----------------------------------------------------------------
c$$$
c$$$c     Call variables
c$$$
c$$$        type(petsc_array),intent(in)  :: petscarray1,petscarray2
c$$$        type(petsc_array)             :: petscarray3
c$$$
c$$$c     Local variables
c$$$
c$$$        integer(4)                    :: ieq
c$$$
c$$$c     Begin program
c$$$
c$$$        call allocatePetscType(petscarray3)
c$$$
c$$$        do ieq=1,neqd
c$$$          petscarray3%array(:,:,:)%var(ieq) =
c$$$     .           petscarray2%array(:,:,:)%var(ieq)
c$$$     .         - petscarray1%array(:,:,:)%var(ieq)
c$$$        enddo
c$$$
c$$$c     End program
c$$$
c$$$      end function substractPetscType
c$$$
c$$$c     mapPetscToVector
c$$$c     #################################################################
c$$$      subroutine mapPetscToVector(x,petscarray)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Maps petsc array petscarray into vector x (without ghost cells)
c$$$c     -----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8),intent(out) :: x(ntotd)
c$$$
c$$$      type (petsc_var),intent(in)  :: petscarray(ilo:ihi,jlo:jhi,klo:khi)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer(4) :: i,j,k,il,jl,kl,ii,ieq
c$$$
c$$$c     Begin program
c$$$
c$$$c     Form vector x
c$$$
c$$$      do k = klo,khi
c$$$        do j = jlo,jhi
c$$$          do i = ilo,ihi
c$$$            call fromGlobalToLocalLimits(i,j,k,il,jl,kl,1,1,1)
c$$$            ii = vecPos(neqd,il,jl,kl,1,1,1)
c$$$            do ieq = 1,neqd
c$$$              x(ii+ieq) = petscarray(i,j,k)%var(ieq)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine mapPetscToVector
c$$$
c$$$c     mapVectorToPetsc
c$$$c     #################################################################
c$$$      subroutine mapVectorToPetsc(petscarray,x)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Maps vector solution into Petsc arrays.
c$$$c     -----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8),intent(in)             :: x(ntotd)
c$$$
c$$$      type (petsc_array),intent(out) :: petscarray
c$$$
c$$$c     Local variables
c$$$
c$$$      integer(4) :: i,j,k,ii,ieq
c$$$
c$$$c     Begin program
c$$$
c$$$c     Initialize petscarray
c$$$
c$$$      call initializePetscType(petscarray)
c$$$
c$$$c     Unpack vector x
c$$$
c$$$      do k = klo,khi
c$$$        do j = jlo,jhi
c$$$          do i = ilo,ihi
c$$$            ii = vecPos(neqd,i,j,k)
c$$$            do ieq = 1,neqd
c$$$              petscarray%array(i,j,k)%var(ieq) = x(ii+ieq)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine mapVectorToPetsc

      end module petsc_var_def
#endif

c module variable_setup
c #####################################################################
      module variable_setup

        use parameters

        use grid

#if defined(petsc)
        use petsc_var_def
#endif

        implicit none

        type :: var_def
          integer :: bconds(6)          !Boundary conditions
          real(8),pointer,dimension(:,:,:) :: array
          character(20):: descr
        end type var_def

        type :: var_array
          integer :: nvar               !Number of variables
          type (var_def),pointer,dimension(:) :: array_var
        end type var_array

        type (var_array) :: u_0

        INTERFACE ASSIGNMENT (=)
          module procedure equateDerivedType
     .                    ,mapStructureToVector
     .                    ,mapVectorToStructure
        END INTERFACE

cc        INTERFACE OPERATOR (-)
cc          module procedure substractDerivedType
cc        END INTERFACE

      contains

c     allocateDerivedType
c     #################################################################
      subroutine allocateDerivedType(varray)

        implicit none

c     Call variables

        type(var_array)  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        varray%nvar = neqd

        if (.not.associated(varray%array_var)) then
          allocate(varray%array_var(neqd))
          do ieq=1,neqd
            allocate(varray%array_var(ieq)
     $                     %array(ilom:ihip,jlom:jhip,klom:khip))
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

        type(var_array)  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        if (associated(varray%array_var)) then
          do ieq=1,neqd
            if (associated(varray%array_var(ieq)%array)) then
              deallocate(varray%array_var(ieq)%array)
            endif
          enddo
          deallocate(varray%array_var)
        endif

c     End program

      end subroutine deallocateDerivedType

c     writeDerivedType
c     #################################################################
      subroutine writeDerivedType(varray,unit,frmt)

c     Call variables

        type(var_array)  :: varray
        integer(4)       :: unit
        logical          :: frmt

c     Local variables

        integer          :: ieq

c     Begin program

        if (frmt) then
#if defined(petsc)
          write (unit,*) ilomg,ihipg,jlomg,jhipg,klomg,khipg
#else
          write (unit,*) ilom,ihip,jlom,jhip,klom,khip
#endif
          write (unit,*) varray%nvar
          do ieq=1,varray%nvar
            write(unit,*) varray%array_var(ieq)%array
            write(unit,*) varray%array_var(ieq)%bconds
            write(unit,*) varray%array_var(ieq)%descr
          enddo
        else
#if defined(petsc)
          write (unit) ilomg,ihipg,jlomg,jhipg,klomg,khipg
#else
          write (unit) ilom,ihip,jlom,jhip,klom,khip
#endif
          write (unit) varray%nvar
          do ieq=1,varray%nvar
            write(unit) varray%array_var(ieq)%array
            write(unit) varray%array_var(ieq)%bconds
            write(unit) varray%array_var(ieq)%descr
          enddo
        endif

c     End program

      end subroutine writeDerivedType

c     readDerivedType
c     #################################################################
      subroutine readDerivedType(varray,unit,format,ierr)

c     Call variables

        type(var_array)  :: varray
        integer(4)       :: unit,ierr
        logical          :: format

c     Local variables

        integer          :: ieq,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        ierr = 0

        if (format) then
#if defined(petsc)
          read (unit,*,iostat=ierr,end=100)
     .               ilomg,ihipg,jlomg,jhipg,klomg,khipg
          call fromGlobalToLocalLimits(ilomg,jlomg,klomg,ilom,jlom,klom
     .                                ,1,1,1)
          call fromGlobalToLocalLimits(ihipg,jhipg,khipg,ihip,jhip,khip
     .                                ,1,1,1)
#else
          read (unit,*,iostat=ierr,end=100)ilom,ihip,jlom,jhip,klom,khip
#endif
          read (unit,*,iostat=ierr,end=100)varray%nvar
          if (ierr /= 0) goto 200
          do ieq=1,varray%nvar
            read(unit,*,iostat=ierr,end=100)varray%array_var(ieq)
     .                 %array(ilom:ihip,jlom:jhip,klom:khip)
            if (ierr /= 0) goto 200
            read(unit,*,iostat=ierr,end=100)varray%array_var(ieq)%bconds
            if (ierr /= 0) goto 200
            read(unit,*,iostat=ierr,end=100)varray%array_var(ieq)%descr
            if (ierr /= 0) goto 200
          enddo
        else
#if defined(petsc)
          read (unit,iostat=ierr,end=100)
     .               ilomg,ihipg,jlomg,jhipg,klomg,khipg
          call fromGlobalToLocalLimits(ilomg,jlomg,klomg,ilom,jlom,klom
     .                                ,1,1,1)
          call fromGlobalToLocalLimits(ihipg,jhipg,khipg,ihip,jhip,khip
     .                                ,1,1,1)
cc          write (*,*) ilom,ihip,jlom,jhip,klom,khip
#else
          read (unit,iostat=ierr,end=100)ilom,ihip,jlom,jhip,klom,khip
#endif
          read (unit,iostat=ierr,end=100) varray%nvar
          if (ierr /= 0) goto 200
          do ieq=1,varray%nvar
            read(unit,iostat=ierr,end=100) varray%array_var(ieq)
     .                 %array(ilom:ihip,jlom:jhip,klom:khip)
            if (ierr /= 0) goto 200
            read(unit,iostat=ierr,end=100) varray%array_var(ieq)%bconds
            if (ierr /= 0) goto 200
            read(unit,iostat=ierr,end=100) varray%array_var(ieq)%descr
            if (ierr /= 0) goto 200
          enddo
cc            write (*,*) varray%array_var(1)
cc     .                 %array(ilom:ihip,jlom:jhip,klom:khip)
cc            stop
        endif

        return

 100    ierr=-2  !EOF
        return

 200    ierr=-1  !Error, but not EOF

c     End program

      end subroutine readDerivedType

c     equateDerivedType
c     #################################################################
      subroutine equateDerivedType(varray2,varray1)
      implicit none
c     -----------------------------------------------------------------
c     Performs varray2 = varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),intent(in)  :: varray1
        type(var_array),intent(out) :: varray2

c     Local variables

        integer(4) :: ieq

c     Begin program

        if (varray2%nvar /= neqd) nullify(varray2%array_var)

        call allocateDerivedType(varray2)

        varray2%nvar = varray1%nvar
        do ieq=1,varray2%nvar
          varray2%array_var(ieq)%bconds = varray1%array_var(ieq)%bconds
cc          write (*,*) ieq,shape(varray1%array_var(ieq)%array)
cc          write (*,*) ieq,shape(varray2%array_var(ieq)%array)
          varray2%array_var(ieq)%array  = varray1%array_var(ieq)%array
          varray2%array_var(ieq)%descr  = varray1%array_var(ieq)%descr
        enddo

c     End program

      end subroutine equateDerivedType

c     initializeDerivedType
c     #################################################################
      subroutine initializeDerivedType(varray)

c     -----------------------------------------------------------------
c     Initializes fields of varray2 using a template (u_0) except for
c     numeric arrays.
c     -----------------------------------------------------------------

c     Call variables

        type(var_array)  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray)

        varray%nvar = u_0%nvar
        do ieq=1,varray%nvar
          varray%array_var(ieq)%bconds = u_0%array_var(ieq)%bconds
          varray%array_var(ieq)%descr  = u_0%array_var(ieq)%descr
cc          varray%array_var(ieq)%array  = 0d0
        enddo

c     End program

      end subroutine initializeDerivedType

ccc     substractDerivedType
ccc     #################################################################
cc      function substractDerivedType(varray2,varray1) result (varray3)
cc
ccc     -----------------------------------------------------------------
ccc     Performs varray3 = varray2 - varray1
ccc     -----------------------------------------------------------------
cc
ccc     Call variables
cc
cc        type(var_array)             :: varray3
cc        type(var_array),intent(in)  :: varray1,varray2
cc
ccc     Local variables
cc
cc        integer          :: ieq
cc
ccc     Begin program
cc
cc        call allocateDerivedType(varray3)
cc
cc        varray3%nvar = varray2%nvar
cc        do ieq=1,varray3%nvar
cc          varray3%array_var(ieq)%bconds = varray2%array_var(ieq)%bconds
cc          varray3%array_var(ieq)%array  = varray2%array_var(ieq)%array
cc     .                                   -varray1%array_var(ieq)%array
cc          varray3%array_var(ieq)%descr  = varray2%array_var(ieq)%descr
cc        enddo
cc
ccc     End program
cc
cc      end function substractDerivedType

c     substractDerivedType
c     #################################################################
      subroutine substractDerivedType(varray2,varray1,varray3)

c     -----------------------------------------------------------------
c     Performs varray3 = varray2 - varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),intent(OUT) :: varray3
        type(var_array),intent(IN)  :: varray1,varray2

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray3)

        varray3%nvar = varray2%nvar
        do ieq=1,varray3%nvar
          varray3%array_var(ieq)%bconds = varray2%array_var(ieq)%bconds
          varray3%array_var(ieq)%array  = varray2%array_var(ieq)%array
     .                                   -varray1%array_var(ieq)%array
          varray3%array_var(ieq)%descr  = varray2%array_var(ieq)%descr
        enddo

c     End program

      end subroutine substractDerivedType

c     varPack
c     #################################################################
      subroutine varPack(array,bcs,desc,veq,varray)

c     -----------------------------------------------------------------
c     Fills structure fields of varray with pertinent data corresponding
c     to variable veq: boundary conditions, bcs; variable description,
c     desc; numerical data, array).
c     -----------------------------------------------------------------

c     Call variables

        integer(4)   :: veq,bcs(6)
        real(8)      :: array(ilom:ihip,jlom:jhip,klom:khip)
        character(*) :: desc

        type (var_array) :: varray

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

        integer(4)    :: nx,ny,nz,veq,bcs(6)
        real(8)       :: array(ilom:ihip,jlom:jhip,klom:khip)
        character(20) :: desc

        type (var_array) :: varray

c     Local variables

c     Begin program

        bcs   = varray%array_var(veq)%bconds
        array = varray%array_var(veq)%array
        desc  = trim(varray%array_var(veq)%descr)

c     End program

      end subroutine varUnPack

c     mapStructureToVector
c     #################################################################
      subroutine mapStructureToVector(x,varray)

c     -----------------------------------------------------------------
c     Maps structure into vector (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(out) :: x(ntotd)

      type (var_array),intent(in)  :: varray

c     Local variables

      integer(4) :: i,j,k,ii,ieq,neq

c     Begin program

c     Unpack vector x

      neq = varray%nvar

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neq,i,j,k,1,1,1)
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
      subroutine mapVectorToStructure(varray,x)

c     -----------------------------------------------------------------
c     Maps vector solution into structure.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(in)           :: x(ntotd)

      type (var_array),intent(out) :: varray

c     Local variables

      integer(4) :: i,j,k,ii,ieq,neq

c     Begin program

c     Initialize varray

      call initializeDerivedType(varray)

c     Unpack vector x

      neq = varray%nvar

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neq,i,j,k,1,1,1)
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

        type (var_array),target :: u_np,u_n,u_graph,u_ic

      contains

c     allocateStructures
c     #################################################################
      subroutine allocateStructures

c     Call variables

c     Local variables

c     Begin program

cc        call allocateDerivedType(u_np)
cc        call allocateDerivedType(u_n )
cc        call allocateDerivedType(u_0)
cccc        call allocateDerivedType(utmp)
cc        call allocateDerivedType(u_graph)

c     End program

      end subroutine allocateStructures

      end module variables

c module generalPurposeFunctions
c ######################################################################
      module generalPurposeFunctions

        use grid

      contains

c     integral
c     ################################################################
      function integral(nx,ny,nz,array,igx,igy,igz,avg)

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on domain (nx)x(ny)x(nz).
c     ---------------------------------------------------------------

      implicit none

c    Call variables

      integer(4) :: igx,igy,igz,nx,ny,nz
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1),integral
      logical    :: avg

c     Local variables

      integer(4) :: i,j,k

      real(8)    :: tvolume,vol,lvolume,lintegral

c     Begin program

c     Integrate

      lintegral = 0d0
      lvolume   = 0d0

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            vol = volume(i,j,k,igx,igy,igz)
            lintegral = lintegral + array(i,j,k)*vol
            lvolume = lvolume + vol
          enddo
        enddo
      enddo

#if defined(petsc)
      call MPI_Allreduce(lintegral,integral,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      call MPI_Allreduce(lvolume  ,tvolume ,1,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
      integral = lintegral
      tvolume  = lvolume
#endif

      if (avg) integral = integral/tvolume

c     End 

      end function integral

c     int2char
c     ################################################################
      function int2char(n) result (chr)

      implicit none

      integer(4)   :: n
      character(10):: chr

      integer(4)   :: i,exp,k,j
      character(3) :: c

      if (abs(n) > 0) then
         exp = int(log(float(n))/log(1d1))
      else
         exp = 0
      endif

      if (n >= 0) then
        chr=''
      else
        chr='-'
      endif

      k = abs(n)
      do i=exp,0,-1
         j = k/10**i
         c = achar(48+j)
         chr = trim(chr)//trim(c)
         if (i > 0 .and. j /= 0) k = mod(k,j*10**i)
      enddo

      end function int2char

      end module generalPurposeFunctions

c module icond
c ######################################################################
      module icond

        implicit none

        integer :: nh1,nh2,nh3

        logical :: odd,random

        real(8),dimension(10) :: pert

      end module icond

c module timeStepping
c ######################################################################
      module timeStepping

        use parameters

        integer(4) :: numtime,ndstep,nrstep,sm_pass,sm_flag,inewtime

        real(8)    :: dt,dtbase,cnfactor,cnf_d,alpha,gammat
        real(8)    :: tmax,dstep,rstep,time,dfreq,dtexp,damp,dt0

cc#if !defined(petsc)
        real(8),allocatable,dimension(:) :: fold,fsrc
cc#endif

        real(8)    :: vx_max,vy_max,vz_max,bx_max,by_max,bz_max

        logical    :: timecorr,restart,source

        real(8)    :: tmrst

        integer(4) :: itime,nrst

        real(8),allocatable,dimension(:) :: cnf,one_over_dt

      end module timeStepping

c module newtongm
c ######################################################################
      module newtongm

        integer(4) ::  maxitnwt,maxitgm,iguess,maxksp,method,global

        real(8)    ::  tolgm,rtol,atol

      end module newtongm

c module counters
c ######################################################################
      module counters

        integer(4) :: itnewt,itgmres,itwhistler
        integer(4) :: gmres_tot,newt_tot,wh_tot

      end module counters

c module iosetup
c ######################################################################
      module iosetup

        integer(4)     :: ilevel

        integer(4)     :: urecord,uinput,uhdf5
        character*(20) :: recordfile,inputfile,hdf5file

        logical        :: debug

      contains

c     readRecordFile
c     #################################################################
      subroutine readRecordFile(unit,itime,time,dt,varray,ierr)

c     -----------------------------------------------------------------
c     Reads record file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer(4) :: ierr,itime,unit
      real(8)    :: time,dt
      type (var_array):: varray

c     Begin program

      ierr = 0

      read (unit,iostat=ierr,end=100) time
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) itime
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) dt
      if (ierr /= 0) goto 200

      call readDerivedType(varray,unit,.false.,ierr)

      return

 100  ierr=-2  !EOF
      return

 200  ierr=-1  !Error, but not EOF

c     End

      end subroutine readRecordFile

c     writeRecordFile
c     #################################################################
      subroutine writeRecordFile(unit,itime,time,dt,varray)

c     -----------------------------------------------------------------
c     Writes record file
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer(4) :: itime,unit
      real(8)    :: time,dt
      type (var_array):: varray

c     Begin program

      write (unit) time
      write (unit) itime
      write (unit) dt

      call writeDerivedType(varray,unit,.false.)

c     End

      end subroutine writeRecordFile

c     closeRecordFile
c     #################################################################
      subroutine closeRecordFile

c     -----------------------------------------------------------------
c     Closes record file
c     -----------------------------------------------------------------

      implicit none

c     Begin program

      close(urecord)

c     End

      end subroutine closeRecordFile

      end module iosetup

c module precond_setup
c ######################################################################
      module precond_setup

        integer(4)    :: nsweep,maxvcyc,precpass
        real(8)       :: mgtol

        character*(10):: precon

      end module precond_setup

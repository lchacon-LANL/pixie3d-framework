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

c     vecPos
c     #################################################################
      function vecPos(neq,i,j,k) result (ijk)

        implicit none

        integer(4) :: ijk,i,j,k,neq

        ijk = neq*(i-ilo + nxl*(j-jlo) + nxl*nyl*(k-klo))
      
      end function vecPos

c     fromGlobalToLocalLimits
c     #################################################################
      subroutine fromGlobalToLocalLimits(ig,jg,kg,il,jl,kl)

        implicit none

        integer(4) :: ig,jg,kg,il,jl,kl

        il = ig - ilog + 1
        jl = jg - jlog + 1
        kl = kg - klog + 1

      end subroutine fromGlobalToLocalLimits

      end module parameters

#if defined(petsc)
c module petsc_var_setup
c #####################################################################
      module petsc_var_setup

        use parameters

        implicit none

        !NVAR is defined via C preprocessor flag
        type :: petsc_var
          real(8) :: var(NVAR)
        end type petsc_var

        type :: petsc_array
          type(petsc_var),pointer,dimension(:,:,:) :: array
        end type petsc_array

      contains

c     allocatePetcsType
c     #################################################################
      subroutine allocatePetscType(petscarray)

        implicit none

c     Call variables

        type(petsc_array)  :: petscarray

c     Local variables

        integer(4)         :: ieq

c     Begin program

        !Error check
        if (neqd /= NVAR) then
          write (*,*) 'Different number of equations in Petsc and F90'
          write (*,*) 'Aborting...'
          stop
        endif

        if (.not.associated(petscarray%array)) then
          allocate(petscarray%array(ilom:ihip,jlom:jhip,klom:khip))
          do ieq=1,neqd
            petscarray%array(:,:,:)%var(ieq) = 0d0
          enddo
        endif

c     End program

      end subroutine allocatePetscType

c     deallocatePetscType
c     #################################################################
      subroutine deallocatePetscType(petscarray)

        implicit none

c     Call variables

        type(petsc_array)  :: petscarray

c     Local variables

c     Begin program

        if (associated(petscarray%array)) deallocate(petscarray%array)

c     End program

      end subroutine deallocatePetscType

c     equatePetscType
c     #################################################################
      subroutine equatePetscType(petscarray2,petscarray1)
      implicit none
c     -----------------------------------------------------------------
c     Performs petscarray2 = petscarray1
c     -----------------------------------------------------------------

c     Call variables

        type(petsc_array),intent(in)  :: petscarray1
        type(petsc_array),intent(out) :: petscarray2

c     Local variables

c     Begin program

        call allocatePetscType(petscarray2)

        petscarray2%array = petscarray1%array

c     End program

      end subroutine equatePetscType

c     initializePetscType
c     #################################################################
      subroutine initializePetscType(petscarray)

c     -----------------------------------------------------------------
c     Initializes fields of petscarray.
c     -----------------------------------------------------------------

c     Call variables

        type(petsc_array)  :: petscarray

c     Local variables

c     Begin program

        call allocatePetscType(petscarray)

c     End program

      end subroutine initializePetscType

c     substractPetscType
c     #################################################################
      function substractPetscType(petscarray2,petscarray1)
     .         result (petscarray3)

c     -----------------------------------------------------------------
c     Performs petscarray3 = petscarray2 - petscarray1
c     -----------------------------------------------------------------

c     Call variables

        type(petsc_array),intent(in)  :: petscarray1,petscarray2
        type(petsc_array)             :: petscarray3

c     Local variables

        integer(4)                    :: ieq

c     Begin program

        call allocatePetscType(petscarray3)

        do ieq=1,neqd
          petscarray3%array(:,:,:)%var(ieq) =
     .           petscarray2%array(:,:,:)%var(ieq)
     .         - petscarray1%array(:,:,:)%var(ieq)
        enddo

c     End program

      end function substractPetscType

c     mapPetscToVector
c     #################################################################
      subroutine mapPetscToVector(x,petscarray)

c     -----------------------------------------------------------------
c     Maps petsc array petscarray into vector x (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(out) :: x(ntotd)

      type (petsc_array),intent(in)  :: petscarray

c     Local variables

      integer(4) :: i,j,k,ii,ieq

c     Begin program

c     Form vector x

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k)
            do ieq = 1,neqd
              x(ii+ieq) = petscarray%array(i,j,k)%var(ieq)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapPetscToVector

c     mapVectorToPetsc
c     #################################################################
      subroutine mapVectorToPetsc(petscarray,x)

c     -----------------------------------------------------------------
c     Maps vector solution into Petsc arrays.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(in)             :: x(ntotd)

      type (petsc_array),intent(out) :: petscarray

c     Local variables

      integer(4) :: i,j,k,ii,ieq

c     Begin program

c     Initialize petscarray

      call initializePetscType(petscarray)

c     Unpack vector x

      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            ii = vecPos(neqd,i,j,k)
            do ieq = 1,neqd
              petscarray%array(i,j,k)%var(ieq) = x(ii+ieq)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapVectorToPetsc

      end module petsc_var_setup
#endif

c module variable_setup
c #####################################################################
      module variable_setup

        use parameters

#if defined(petsc)
        use petsc_var_setup
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
#if defined(petsc)
     .                    ,mapPetscToStructure
     .                    ,mapStructureToPetsc
#endif
        END INTERFACE

        INTERFACE OPERATOR (-)
          module procedure substractDerivedType
#if defined(petsc)
     .                    ,substractPetscType
#endif
        END INTERFACE

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
cc            write (*,*) ieq,shape(varray%array_var(ieq)%array)
          enddo
        endif

cc        do ieq=1,neqd
cc          if (.not.associated(varray%array_var(ieq)%array)) then
cc            allocate(varray%array_var(ieq)
cc     $                     %array(ilom:ihip,jlom:jhip,klom:khip))
cc            varray%array_var(ieq)%array = 0d0
cccc          write (*,*) ieq,shape(varray%array_var(ieq)%array)
cc          endif
cc        enddo

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

        integer          :: ieq

c     Begin program

        ierr = 0

        if (format) then
#if defined(petsc)
          read (unit,*,iostat=ierr,end=100)
     .               ilomg,ihipg,jlomg,jhipg,klomg,khipg
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

c     substractDerivedType
c     #################################################################
      function substractDerivedType(varray2,varray1) result (varray3)

c     -----------------------------------------------------------------
c     Performs varray3 = varray2 - varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),intent(in)  :: varray1,varray2
        type(var_array)             :: varray3

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

      end function substractDerivedType

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
            ii = vecPos(neq,i,j,k)
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
            ii = vecPos(neq,i,j,k)
            do ieq=1,neq
              varray%array_var(ieq)%array(i,j,k) = x(ii+ieq)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapVectorToStructure

#if defined(petsc)
c     mapPetscToStructure
c     #################################################################
      subroutine mapPetscToStructure(varray,petscarray)

c     -----------------------------------------------------------------
c     Maps structure into Petsc format.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type (var_array  ),intent(out)  :: varray
      type (petsc_array),intent(in)   :: petscarray

c     Local variables

      integer(4) :: i,j,k,ieq

c     Begin program

c     Initialize structure

      call initializeDerivedType(varray)

c     Perform mapping

      do k = klom,khip
        do j = jlom,jhip
          do i = ilom,ihip
            do ieq = 1,neqd
              varray%array_var(ieq)%array(i,j,k)
     .             = petscarray%array(i,j,k)%var(ieq)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapPetscToStructure

c     mapStructureToPetsc
c     #################################################################
      subroutine mapStructureToPetsc(petscarray,varray)

c     -----------------------------------------------------------------
c     Maps structure into petscarray.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type (var_array  ),intent(in ) :: varray
      type (petsc_array),intent(out) :: petscarray

c     Local variables

      integer(4) :: i,j,k,ieq

c     Begin program

c     Initialize petsc array

      call initializePetscType(petscarray)

c     Perform mapping

      do k = klom,khip
        do j = jlom,jhip
          do i = ilom,ihip
            do ieq = 1,neqd
              petscarray%array(i,j,k)%var(ieq)
     .             = varray%array_var(ieq)%array(i,j,k)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapStructureToPetsc
#endif

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

#if defined(petsc)
        use mpi
#endif

      contains

c     integral
c     ################################################################
      real(8) function integral(nx,ny,nz,array,igx,igy,igz,avg)

c     ---------------------------------------------------------------
c     Integrates array(i,j,k) on domain (nx)x(ny)x(nz).
c     ---------------------------------------------------------------

      implicit none

c    Call variables

      integer(4) :: igx,igy,igz,nx,ny,nz
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1)
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
         if (i > 0 ) k = mod(k,j*10**i)
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

#if !defined(petsc)
        real(8),allocatable,dimension(:) :: fold,fsrc
#endif

        real(8)    :: vx_max,vy_max,vz_max,bx_max,by_max,bz_max

        logical    :: timecorr,restart,source

        real(8)    :: tmrst

        integer(4) :: itime,nrst

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

#if defined(petsc)
        use local_mpi
#endif

        integer(4)     :: ilevel

        integer(4)     :: urecord,uinput,uhdf5
        character*(20) :: recordfile,inputfile,hdf5file

        logical        :: debug

      contains

c     readRecord
c     #################################################################
      subroutine readRecord(unit,itime,time,dt,varray,ierr)

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

      end subroutine readRecord

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

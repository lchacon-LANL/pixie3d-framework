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

#if defined(petsc)
        !NVAR is defined via C preprocessor flag
        type :: petsc_var
          real(8) :: var(NVAR)
        end type petsc_var
#endif

        type :: var_def
          integer :: bconds(6)          !Boundary conditions
          real(8),pointer,dimension(:,:,:) :: array
          character(20) :: descr
        end type var_def

        type :: vec_def
          integer :: bconds(6,3)        !Boundary conditions
          real(8),pointer,dimension(:,:,:,:) :: vec
          character(20) :: descr        !Description
          logical       :: cnv          !Whether it is contravariant or not
        end type vec_def

        type :: var_array
          integer :: nvar               !Number of variables
          type (var_def),pointer,dimension(:) :: array_var
        end type var_array

        type :: aux_array
          integer :: nvar,nvec          !Number of variables, vectors
          type (vec_def),pointer,dimension(:) :: vec_list
          type (var_def),pointer,dimension(:) :: var_list
        end type aux_array

        type :: patch
          integer  :: nbc_seq
          integer,pointer,dimension(:,:):: bc_seq
          type(var_array)  ,pointer     :: u_0
          type(var_array)  ,pointer     :: u_n
          type(var_array)  ,pointer     :: u_np
          type(var_array)  ,pointer     :: u_nm
          type(var_array)  ,pointer     :: u_graph
          type(grid_mg_def),pointer     :: gparams
          type(aux_array)  ,pointer     :: aux
        end type patch

        type (patch),pointer :: gv

        type(var_array),pointer :: u_0

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
            enddo
          endif

          if (nvec > 0) then
            allocate(aux%vec_list(nvec))

            do ieq=1,nvec
              allocate(aux%vec_list(ieq)
     $                    %vec(ilom:ihip,jlom:jhip,klom:khip,3))
              aux%vec_list(ieq)%vec = 0d0
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

        integer          :: ieq,nvar,nvec

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
      subroutine allocateDerivedType(varray)

        implicit none

c     Call variables

        type(var_array),pointer  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        if (.not.associated(varray)) then
          allocate(varray)

          varray%nvar = neqd

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

c     writeDerivedType
c     #################################################################
      subroutine writeDerivedType(varray,unit,frmt)

c     Call variables

        type(var_array),pointer  :: varray
        integer          :: unit
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
            write(unit,*) varray%array_var(ieq)
     .                          %array(ilom:ihip,jlom:jhip,klom:khip)
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
            write(unit) varray%array_var(ieq)
     .                        %array(ilom:ihip,jlom:jhip,klom:khip)
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

        type(var_array),pointer  :: varray
        integer          :: unit,ierr
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
        write (*,*) 'readDerivedType: EOF in file unit',unit
        return

 200    ierr=-1  !Error, but not EOF
        write (*,*) 'readDerivedType: Error reading file unit',unit

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
      subroutine initializeDerivedType(varray)

c     -----------------------------------------------------------------
c     Initializes fields of varray2 using a template (u_0) except for
c     numeric arrays.
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),pointer  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray)

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
      subroutine mapStructureToVector(x,varray)

c     -----------------------------------------------------------------
c     Maps structure into vector (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(out)      :: x(ntotd)

      type(var_array),pointer  :: varray

c     Local variables

      integer    :: i,j,k,ii,ieq,neq

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

      real(8),intent(in)      :: x(ntotd)

      type(var_array),pointer :: varray

c     Local variables

      integer    :: i,j,k,ii,ieq,neq

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

        type(var_array),pointer :: u_np,u_n,u_nm,u_graph,u_ic
  
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

      integer    :: igx,igy,igz,nx,ny,nz
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1),integral
      logical    :: avg

c     Local variables

      integer    :: i,j,k

      real(8)    :: tvolume,vol,lvolume,lintegral

c     Begin program

c     Integrate

      lintegral = 0d0
      lvolume   = 0d0

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            vol = gmetric%grid(igx)%dvol(i,j,k)

            if (isSYM(i,igx,1,0)) vol = 0.5*vol
            if (isSYM(i,igx,1,1)) vol = 0.5*vol
            if (isSYM(j,igy,2,0)) vol = 0.5*vol
            if (isSYM(j,igy,2,1)) vol = 0.5*vol
            if (isSYM(k,igz,3,0)) vol = 0.5*vol
            if (isSYM(k,igz,3,1)) vol = 0.5*vol

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

      end module generalPurposeFunctions

ccc module iosetup
ccc ######################################################################
cc      module iosetup
cc
cc        integer        :: ilevel
cc
cc        integer        :: urecord=10000,usave=20000,uinput=25
cc        character*(50) :: recordfile='record.bin',recordsave,inputfile
cc
cc        integer        :: nfiles,ifile
cc        integer       ,allocatable,dimension(:) :: murecord,musave
cc        character*(50),allocatable,dimension(:) :: mrecfile,mrecsave
cc
cc        logical        :: debug,clean=.false.,tst_flg=.false.
cc
cc      contains
cc
ccc     readRecordFile
ccc     #################################################################
cc      subroutine readRecordFile(unit,itime,time,dt,varray,ierr)
cc
ccc     -----------------------------------------------------------------
ccc     Reads record file
ccc     -----------------------------------------------------------------
cc
cc      use variables
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: ierr,itime,unit
cc      real(8)    :: time,dt
cc      type(var_array),pointer:: varray
cc
ccc     Begin program
cc
cc      ierr = 0
cc
cc      read (unit,iostat=ierr,end=100) time
cc      if (ierr /= 0) goto 200
cc
cc      read (unit,iostat=ierr,end=100) itime
cc      if (ierr /= 0) goto 200
cc
cc      read (unit,iostat=ierr,end=100) dt
cc      if (ierr /= 0) goto 200
cc
cc      call readDerivedType(varray,unit,.false.,ierr)
cc
cc      return
cc
cc 100  ierr=-2  !EOF
cccc      write (*,*) 'readRecordFile: EOF in file unit',unit
cc      return
cc
cc 200  ierr=-1  !Error, but not EOF
cc      write (*,*) 'readRecordFile: Error reading file unit',unit
cc
ccc     End
cc
cc      end subroutine readRecordFile
cc
ccc     writeRecordFile
ccc     #################################################################
cc      subroutine writeRecordFile(unit,itime,time,dt,varray)
cc
ccc     -----------------------------------------------------------------
ccc     Writes record file
ccc     -----------------------------------------------------------------
cc
cc      use variables
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer    :: itime,unit
cc      real(8)    :: time,dt
cc      type(var_array),pointer:: varray
cc
ccc     Begin program
cc
cc      write (unit) time
cc      write (unit) itime
cc      write (unit) dt
cc
cc      call writeDerivedType(varray,unit,.false.)
cc
ccc     End
cc
cc      end subroutine writeRecordFile
cc
ccc     closeRecordFile
ccc     #################################################################
cc      subroutine closeRecordFile
cc
ccc     -----------------------------------------------------------------
ccc     Closes record file
ccc     -----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Begin program
cc
cc      close(urecord)
cc
ccc     End
cc
cc      end subroutine closeRecordFile
cc
cc      end module iosetup

ccc module graphics_variables
ccc ######################################################################
cc      module graphics_variables
cc
cc        use xdraw_io
cc
cc        integer        :: ngroups
cc
cc        real(8), allocatable, dimension(:) :: xl,yl,zl
cc
cc        type :: graph_var_def
cc          real(8),pointer,dimension(:,:,:) :: array
cc          character(20)                    :: descr
cc          character(20)                    :: vector_name=''
cc        end type graph_var_def
cc
cc        type :: graph_group
cc          type (graph_var_def),dimension(ngraph) :: array_graph
cc          logical                                :: cartesian
cc          character(20)                          :: descr
cc        end type graph_group
cc
cc        type (graph_group),pointer,dimension(:) :: graph
cc
cc        integer   ,allocatable,dimension(:) :: nqty
cc
cc        logical    :: plot
cc
cc#if defined(petsc)
cc        logical    :: hdf_plot=.true.
cc#else
cc        logical    :: hdf_plot=.false.
cc#endif
cc
cc        integer    :: iming,imaxg,jming,jmaxg,kming,kmaxg,igroup
cc     .               ,iggx,iggy,iggz,iig,jjg,kkg
cc
cc        integer    :: sel_diag (xdraw_cont_lim)
cc     .               ,sel_graph(xdraw_cont_lim)
cc
cc        integer    :: ndplot
cc        real(8)    :: dplot,tmplot
cc
cc        integer    :: iplot,jplot,kplot
cc
cc        !Profile conf
cc        type :: prof_info
cc          integer        :: line = 1
cc          character(3)   :: label = 'x'
cc          integer        :: coords(3) = (/1,1,1/)
cc        end type prof_info
cc
cc        type(prof_info),save :: prof_conf
cc
cc        !Contour conf
cc        type :: cont_info
cc          integer        :: plane = 3
cc          character(3)   :: label(2) = (/'x','y'/)
cc          integer        :: coords(3)= (/1,1,1/)
cc        end type cont_info
cc
cc        type(cont_info),save :: cont_conf
cc
cc      end module graphics_variables
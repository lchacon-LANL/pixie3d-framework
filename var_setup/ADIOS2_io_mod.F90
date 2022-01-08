! module ADIOS2_io
! ######################################################################
      module ADIOS2_mod

#if defined(ADIOS2)
        use io
        
        use adios2

        use variables

        implicit none

        type(adios2_adios) :: adios2obj
        type(adios2_engine) :: engine, rengine
        type(adios2_io) :: aio
        !! only true for the first "write" step at each run
        logical :: isfirst = .true.

        ! ADIOS2 needs a varia        !!! ADIOS2 variables
        !!! Used to dump main record file (to be postprocessed)
        character(1024) :: recordfile=ADIOS2_FILE

        ! ADIOS needs a variableble containing the MPI_WORLD_COMM communicator
        integer :: adios2_world_comm

        !!! variables for reading 
        integer*8 :: adios2_fh, adios2_gh   ! file handler and group handler
        integer*8 :: adios2_tstart, adios2_tend, adios2_ntsteps ! timesteps available
        integer*8,dimension(0:3) :: adios2_vardims
        integer   :: adios2_varcount       ! number of variables in file (~32)
        integer   :: adios2_attrcount      ! number of attributes in file (~16)
        integer   :: adios2_nvar           ! number of /var/ variables in file (~8)
        integer   :: adios2_tread          ! 0 ... ntsteps-1
        integer   :: adios2_err            ! error handler
        integer   :: adios2_attrtype,adios2_attrsize,adios2_vartype,adios2_ndim,adios2_tdim
        character(30), dimension(50) :: adios2_varlist   ! list of all vars 
        character(30), dimension(50) :: adios2_attrlist  ! list of all attributes 
        character(30), dimension(50) :: adios2_varnames  ! list of /var/ vars 

        logical,private :: adios2_debug=.true.,adios2_append_recordfile=.false.

#if !defined(ADIOS2_BUFFER_MB)
        !ADIOS2 IO buffer 
        integer,private :: ADIOS2_BUFFER_MB=100  !In MB
#endif

#if !defined(ADIOS2_METHOD)
        !ADIOS2 method 
        character(20) :: ADIOS2_METHOD="MPI"
#endif     
      contains

!     setADIOS2AppendIOMode
!     ##############################################################
      subroutine setADIOS2AppendIOMode

        implicit none

!       Call variables

!       Local variables

        adios2_append_recordfile = .true. !Set ADIOS2 recordfile to append

      end subroutine setADIOS2AppendIOMode

!     init_ADIOS2_IO
!     #################################################################
      function init_ADIOS2_IO() result(ierr)

!     -----------------------------------------------------------------
!     Initializes ADIOS2 IO according to adios_config.xml
!     -----------------------------------------------------------------

      implicit none

!     Global variables

      integer :: ierr
      
!     Local variables

!     Begin program

!     Initialize ADIOS2 for componentized I/O

      if (.not.is_file('adios_config.xml')) then
         call pstop('init_ADIOS2_IO','ADIOS2 config file not present')
      else
         call MPI_Comm_dup (MPI_COMM_WORLD,adios2_world_comm,ierr)
         if (.not.adios2obj%valid) then
            call adios2_logging('init')
            call adios2_init(adios2obj,'adios_config.xml', adios2_world_comm,.true.,ierr)
            call adios2_check_err(ierr, 'Problem in init')
         endif
      endif

      end function init_ADIOS2_IO

!     destroy_ADIOS2_IO
!     #################################################################
      function destroy_ADIOS2_IO() result(ierr)

!     -----------------------------------------------------------------
!     Initializes ADIOS IO according to adios_config.xml
!     -----------------------------------------------------------------

      implicit none

!     Global variables

      integer :: ierr

      call adios2_logging('Terminating IO')

      if (.not.isfirst) then
         call adios2_logging('write close')
         call adios2_close(engine, ierr)
         call adios2_check_err(ierr, 'Problem in write close')
         call adios2_logging('finalize')
         call adios2_finalize(adios2obj, ierr)
      endif
      call MPI_Comm_free(adios2_world_comm,ierr)

      end function destroy_ADIOS2_IO

!     writeADIOS2RecordFile
!     #################################################################
      subroutine writeADIOS2RecordFile(file,itime,time,dt,gammat,varray,init)

!     -----------------------------------------------------------------
!     Writes record file
!     -----------------------------------------------------------------

      implicit none

!     Call variables

      integer    :: itime
      real(8)    :: time,dt,gammat
      type(var_array),pointer:: varray
      character(*) :: file
      logical,optional :: init

!     Local variables

      logical :: rinit

      integer*8    :: handle, totalsize, groupsize
      integer      :: err
      character(2) :: mode='a'//char(0)
      integer      :: adios2_mode
      logical      :: addl_write  ! true: write names and bconds

      integer      :: xsize, ysize, zsize
      integer      :: xloghost, xhighost
      integer      :: yloghost, yhighost
      integer      :: zloghost, zhighost
      integer      :: xoffset, yoffset, zoffset

      type(adios2_variable) :: var
      integer       :: istatus

!     Begin program

      if (PRESENT(init)) then
         rinit = init
      else
         rinit = .not.adios2_append_recordfile
      endif

!     Offsets to include ghosts
!     1. internal processes: no ghost cells written out

      xloghost = 1
      yloghost = 1
      zloghost = 1
      xhighost = -1
      yhighost = -1
      zhighost = -1
      
!     2. processes that has some face 

      !x=0 face
      if (ilomg == 0) xloghost = 0

      !x=1 face
      if (ihipg == nxd+1) xhighost = 0

      !y=0 face
      if (jlomg == 0) yloghost = 0

      !y=1 face
      if (jhipg == nyd+1) yhighost = 0

      !z=0 face
      if (klomg == 0) zloghost = 0

      !z=1 face
      if (khipg == nzd+1) zhighost = 0

!     calculate size actually written for each dimension
      xsize = ihip+xhighost-ilom-xloghost+1  
      ysize = jhip+yhighost-jlom-yloghost+1  
      zsize = khip+zhighost-klom-zloghost+1  

!     calculate offsets in the global array
      xoffset = ilomg + xloghost  
      yoffset = jlomg + yloghost  
      zoffset = klomg + zloghost  
        
!     Create/Append adios file

      !!jyc: rinit is true only with no restart. With restart, rinit is false.
      if (rinit) then
         mode = 'w'//char(0)
         adios2_mode = adios2_mode_write
      else 
         mode = 'a'//char(0)
         adios2_mode = adios2_mode_append
      endif

      !!jyc: isfirst will be true only once at each run.
      addl_write=.false.
      if (isfirst) then
         isfirst = .false.

         if (my_rank == 0.and.rinit) then
            ! Process 0 in the first timestep writes 
            ! the boundary conditions and variable names
            addl_write = .true.
         endif

         call adios2_declare_io (aio, adios2obj, "record", err)
         call adios2_set_engine (aio, "BP4", err)
         call adios2_define_variable (var, aio, "time"  , adios2_type_dp, err)
         call adios2_define_variable (var, aio, "itime" , adios2_type_integer4, err)
         call adios2_define_variable (var, aio, "dt"    , adios2_type_dp, err)
         call adios2_define_variable (var, aio, "gammat", adios2_type_dp, err)
         call defineDerivedTypeADIOS2 (aio, varray)

         call adios2_logging('write access mode: '//trim(mode))
         call adios2_logging('write open')
         call adios2_open(engine,aio,file,adios2_mode,adios2_world_comm,err)
         call adios2_check_err(err,'Could not open file for writing')
      else
         call adios2_logging('not first write call: no attrs written')
      endif

      call adios2_logging('begin write step')
      call adios2_begin_step(engine, adios2_step_mode_append, err)
      call adios2_check_err(err, 'Problem in begin write step')

      if (my_rank == 0) call adios2_put(engine,"time"  ,time  ,err)
      if (my_rank == 0) call adios2_put(engine,"itime" ,itime ,err)
      if (my_rank == 0) call adios2_put(engine,"dt"    ,dt    ,err)
      if (my_rank == 0) call adios2_put(engine,"gammat",gammat,err)

      call writeDerivedTypeADIOS2(engine, varray, addl_write)
      
      !End adios step
      call adios2_logging('end write step')
      call adios2_end_step(engine, err)
      call adios2_check_err(err, 'Problem end write step')

      end subroutine writeADIOS2RecordFile

!     defineDerivedTypeADIOS2
!     #################################################################
      subroutine defineDerivedTypeADIOS2(aio, varray)

        implicit none

!     Call variables
        type(adios2_io)          :: aio
        integer*8                :: handle
        type(var_array),pointer  :: varray

!     Local variables

        integer          :: ieq, err, &
                            xoffset, yoffset, zoffset, &
                            xsize, ysize, zsize

        integer          :: xloghost, xhighost
        integer          :: yloghost, yhighost
        integer          :: zloghost, zhighost
        integer,dimension(6*varray%nvar) :: bconds  ! all boundary conditions as one array
        character(20)    :: vvar      ! /var/v<idx> 
        character(20)    :: vname     ! /name/v<idx> 
        type(adios2_variable) :: var
        type(adios2_attribute) :: attr

!     Begin program

!   global size of the output array = nxd+2 x nyd+2 x nzd+2 with ghost cells on the faces

!   ihipg = nxd+1    -> x=1 face
!   ihipg = nxd+1    -> x=1 face
!   j... -> y face
!   k... -> z face
!
!   global position of an internal processor in the 3D mesh
!     from  ilog = ilomg+1  to ihig = ihipg-1
!     from  jlog = jlomg+1  to jhig = jhipg-1
!     from  klog = klomg+1  to khig = khipg-1
!   global position for boundary proc
!     if (ilomg == 0) ilog = ilomg
!     if (ihipg == nxd+1) ihig = ihipg
!
!   local sizes: ilom:ihip   but ilom+1:ihip-1 for internal processor
!   local indices start from 0 for ghost cell

!     Offsets to include ghosts
!     1. internal processes: no ghost cells written out
        xloghost = 1
        yloghost = 1
        zloghost = 1
        xhighost = -1
        yhighost = -1
        zhighost = -1
      
!     2. processes that has some face 
!       x=0 face
        if (ilomg == 0) xloghost = 0

!       x=1 face
        if (ihipg == nxd+1) xhighost = 0

!       y=0 face
        if (jlomg == 0) yloghost = 0

!       y=1 face
        if (jhipg == nyd+1) yhighost = 0

!       z=0 face
        if (klomg == 0) zloghost = 0

!       z=1 face
        if (khipg == nzd+1) zhighost = 0

!       calculate offsets in the global array
        xoffset = ilomg + xloghost  
        yoffset = jlomg + yloghost  
        zoffset = klomg + zloghost  

!       calculate size actually written for each dimension
        xsize = ihip+xhighost-ilom-xloghost+1  
        ysize = jhip+yhighost-jlom-yloghost+1  
        zsize = khip+zhighost-klom-zloghost+1 

        call adios2_define_variable (var, aio, "nvar", adios2_type_integer4, err)

        do ieq=1,varray%nvar
            write (vvar,  '("/var/v",i0)')  ieq
            ! write /var/v<ieq>
            call adios2_define_variable (var, aio, vvar, adios2_type_dp, 3, &
                  int((/ nxd+2,nyd+2,nzd+2 /), kind=8), &
                  int((/ xoffset,yoffset,zoffset /), kind=8), &
                  int((/ xsize,ysize,zsize /), kind=8), .true., err)

            write (vname, '("/name/v",i0,"/description")') ieq
            call adios2_define_attribute (attr, aio, vname, "Name of "//trim(vvar), err)

            ! write /name/v<ieq>
            write (vname, '("/name/v",i0)') ieq
            call adios2_define_variable (var, aio, vname, adios2_type_string, err)
        enddo

        call adios2_define_variable (var, aio, "bconds", adios2_type_integer4, 1, &
             int((/ varray%nvar*6 /), kind=8), &
             int((/ 0 /), kind=8), &
             int((/ varray%nvar*6 /), kind=8), .true., err)                 

      end subroutine defineDerivedTypeADIOS2

!     writeDerivedTypeADIOS2
!     #################################################################
      subroutine writeDerivedTypeADIOS2(engine, varray, addl_write)

        implicit none

!     Call variables
        type(adios2_engine)      :: engine
        integer*8                :: handle
        type(var_array),pointer  :: varray
        logical                  :: addl_write

!     Local variables

        integer          :: ieq, err, &
                            xoffset, yoffset, zoffset, &
                            xsize, ysize, zsize

        integer          :: xloghost, xhighost
        integer          :: yloghost, yhighost
        integer          :: zloghost, zhighost
        integer,dimension(6*varray%nvar) :: bconds  ! all boundary conditions as one array
        character(20)    :: vvar      ! /var/v<idx> 
        character(20)    :: vname     ! /name/v<idx> 

!     Begin program

!        write (*,*) "# addl_write=",addl_write
        
!   global size of the output array = nxd+2 x nyd+2 x nzd+2 with ghost cells on the faces

!   ihipg = nxd+1    -> x=1 face
!   ihipg = nxd+1    -> x=1 face
!   j... -> y face
!   k... -> z face
!
!   global position of an internal processor in the 3D mesh
!     from  ilog = ilomg+1  to ihig = ihipg-1
!     from  jlog = jlomg+1  to jhig = jhipg-1
!     from  klog = klomg+1  to khig = khipg-1
!   global position for boundary proc
!     if (ilomg == 0) ilog = ilomg
!     if (ihipg == nxd+1) ihig = ihipg
!
!   local sizes: ilom:ihip   but ilom+1:ihip-1 for internal processor
!   local indices start from 0 for ghost cell

!     Offsets to include ghosts
!     1. internal processes: no ghost cells written out
        xloghost = 1
        yloghost = 1
        zloghost = 1
        xhighost = -1
        yhighost = -1
        zhighost = -1
      
!     2. processes that has some face 
!       x=0 face
        if (ilomg == 0) xloghost = 0

!       x=1 face
        if (ihipg == nxd+1) xhighost = 0

!       y=0 face
        if (jlomg == 0) yloghost = 0

!       y=1 face
        if (jhipg == nyd+1) yhighost = 0

!       z=0 face
        if (klomg == 0) zloghost = 0

!       z=1 face
        if (khipg == nzd+1) zhighost = 0

!       calculate offsets in the global array
        xoffset = ilomg + xloghost  
        yoffset = jlomg + yloghost  
        zoffset = klomg + zloghost  

        xsize = ihip+xhighost-ilom-xloghost+1
        ysize = jhip+yhighost-jlom-yloghost+1
        zsize = khip+zhighost-klom-zloghost+1

        if (my_rank == 0) call adios2_put (engine, "nvar", varray%nvar, err)

        do ieq=1,varray%nvar
            write (vvar,  '("/var/v",i0)')  ieq
            ! write /var/v<ieq>
            call adios2_put (engine, vvar, &
                  varray%array_var(ieq)%array( &
                        ilom+xloghost:ihip+xhighost, &
                        jlom+yloghost:jhip+yhighost, &
                        klom+zloghost:khip+zhighost), &
                        adios2_mode_sync, err)

            if (addl_write) then
                ! write /name/v<ieq>
                write (vname, '("/name/v",i0)') ieq
                call adios2_put (engine, vname, trim(varray%array_var(ieq)%descr)//char(0), adios2_mode_sync, err)
                if (adios2_debug) then
                   write (*, '(" ADIOS2 INFO: attr write ",a,": ",a)') trim(vname), trim(varray%array_var(ieq)%descr)
                endif
            endif
        enddo

        if (addl_write) then 
            do ieq=1,varray%nvar
                bconds( (ieq-1)*6+1:ieq*6 ) = varray%array_var(ieq)%bconds(:)
                if (adios2_debug) then
                  write (*,'(" ADIOS2 INFO: attr write bconds=(",6i3")")') bconds( (ieq-1)*6+1:ieq*6 )
                endif
            enddo
            call adios2_put (engine, "bconds", bconds, adios2_mode_sync, err)
        endif

      end subroutine writeDerivedTypeADIOS2

!     openADIOS2RecordFileForRead
!     #################################################################
      subroutine openADIOS2RecordFileForRead(ierr,file)!,nvar,nx,ny,nz)

        implicit none

!     Call variables

        integer,intent(out) :: ierr
        character(*),optional :: file
!        integer,optional,intent(out) :: nvar,nx,ny,nz

!     Local variables

        integer :: gcnt, vcnt, acnt, nullidx, ts, lastts,i
        character (len=20), dimension(1)  :: gnamelist  ! list of groups (=1)
        character(len=len(recordfile)) :: rfile
        logical :: inq_adios_var

        type(adios2_variable) :: var
        type(adios2_attribute) :: attr
        real(8) :: tmp

!     Begin program

        if (PRESENT(file)) then
          rfile = file
        else
          rfile = recordfile
        endif

        call MPI_Comm_dup (MPI_COMM_WORLD, adios2_world_comm, ierr)

        !Open ADIOS file
        call adios2_declare_io(aio, adios2obj, 'record.read', ierr)

        call adios2_logging('read-only open')
        call adios2_open(rengine, aio,trim(rfile), adios2_mode_read, adios2_world_comm, ierr)
        call adios2_check_err(ierr, 'Problem in read-only open')
        
        call adios2_get(rengine, 'time', tmp, ierr)

!         call adios2_steps(adios2_ntsteps,rengine,ierr)

!         !Inquire ADIOS file
!         call adios_inq_file(adios_fh,vcnt,acnt,adios_tstart,adios2_ntsteps,gnamelist,ierr)

!         adios_tend = adios_tstart + adios2_ntsteps - 1

!         !Open ADIOS group
!         call adios_gopen (adios_fh,adios_gh,gnamelist(1),adios_varcount,adios_attrcount,ierr)
         

!         !Inquire ADIOS group
!         call adios_inq_group(adios_gh,adios_varlist,adios_attrlist,ts,lastts,ierr)

!         ! count /var/ variables
!         inq_adios_var = (PRESENT(nx).or.PRESENT(ny).or.PRESENT(nz))
!         adios2_nvar=0
!         do i=1,adios_varcount
!           if (index(adios_varlist(i),"/var/") == 1) then
!             adios2_nvar = adios2_nvar+1
!             !Check problem size
!             call adios_inq_var(adios_gh,trim(adios_varlist(i)),adios_vartype,adios_ndim,adios_vardims,adios_tdim,ierr)
! c$$$                  write (*,*) "Var="//adios_varlist(i) //" Dimensions: "
! c$$$     .                 , adios_vardims,adios_ndim,adios_tdim 
!           endif
!         enddo

        adios2_tread = 0  ! next time read first timestep

        !!jyc: With Adios2, the number of variables can vary step by step.
        !!jyc: We cannot have this info at this stage
        
        ! if (adios2_debug) then
        !   write (*,"(a,i0,a,a)") " Proc=",my_rank,": opened ",trim(rfile)

        !   if (my_rank == 0) then
        !     write (*,"(a, i5)") " Number of vars: ", adios_varcount
        !     write (*,*) "Variables:"
        !     do i=1,adios_varcount
        !       write (*,"(i5, a, a)") i,")  ", trim(adios_varlist(i))
        !     enddo
        !     write (*,"(a, 3i5)") " Dimensions: ", adios_vardims(1:3)
        !   endif
        ! endif

        !if (adios2_debug) then
        !   print *, 'PRESENT(nvar)', PRESENT(nvar)
        !   print *, 'PRESENT(nx)', PRESENT(nx), nx
        !   print *, 'PRESENT(ny)', PRESENT(ny), ny
        !   print *, 'PRESENT(nz)', PRESENT(nz), nz
        !endif
        !!jyc: is this critical?
        !if (PRESENT(nvar)) nvar = adios2_nvar
        !if (PRESENT(nx)) nx = adios_vardims(1)-2
        !if (PRESENT(ny)) ny = adios_vardims(2)-2
        !if (PRESENT(nz)) nz = adios_vardims(3)-2

      end subroutine openADIOS2RecordFileForRead

!     closeADIOS2RecordFileForRead
!     #################################################################
      function closeADIOS2RecordFileForRead() result(ierr)

        implicit none

        integer :: ierr

        call adios2_logging('read close')
        call adios2_close(rengine, ierr)
        call adios2_check_err(ierr, 'Problem in read close')

      end function closeADIOS2RecordFileForRead

!     openADIOS2RestartFileForRead
!     #################################################################
      function openADIOS2RestartFileForRead(file) result(ufile)

!     -----------------------------------------------------------------
!     Opens ADIOS restart file. Does NOT use "adios_config.xml" file.
!     -----------------------------------------------------------------

      implicit none

!     Global variables

      integer :: ufile
      character(*),optional :: file
      integer :: ierr

!     Local variables

!     Begin program

      ufile = find_unit(12345)

      call MPI_Comm_dup (MPI_COMM_WORLD,adios2_world_comm,ierr)
      call adios2_logging('read init')
      call adios2_init(adios2obj,'adios_config.xml',adios2_world_comm,.true.,ierr)
      call adios2_check_err(ierr, 'Problem in read init')

      call openADIOS2RecordFileForRead(adios2_err,file=file)
      if (adios2_err /= 0) then
        call pstop('openRestartFileForRead','Error reading ADIOS restart file')
      endif

      end function openADIOS2RestartFileForRead

!     readADIOS2RecordFile
!     #################################################################
      function readADIOS2RecordFile(unit,itime,time,dt,gammat,varray) result(ierr)

!     -----------------------------------------------------------------
!     Reads record file
!     -----------------------------------------------------------------

      implicit none

!     Call variables

      integer    :: ierr,itime,unit
      real(8)    :: time,dt,gammat

      type(var_array),pointer :: varray

!     Local variables

      integer   :: tmp(1),lerr(1),gerr(1)

      !ADIOS variables
      integer*8 :: start(1), readcount(1)  ! dimensions are 64bit in adios
      integer   :: istatus

!     Begin program

      ierr = 0

      if (adios2_debug.and.(my_rank.eq.0)) then
         write (*,"(a,i5)") " ADIOS2 INFO: read time=", adios2_tread
      endif

      call adios2_logging('begin read step')
      call adios2_begin_step(rengine, adios2_step_mode_read, 0.0, istatus, ierr)
      call adios2_check_err(ierr, 'Problem in begin read step')
      
      if ((istatus.eq.0).and.(ierr.eq.0)) then
         call adios2_get(rengine,"time",time,ierr)
         call adios2_get(rengine,"itime",itime,ierr)
         call adios2_get(rengine,"dt",dt,ierr)
         call adios2_get(rengine,"gammat",gammat,ierr)
         call readDerivedTypeADIOS2(rengine,varray,adios2_tread,ierr)
         call adios2_end_step(rengine, ierr)
         call adios2_logging('end read step')

         adios2_tread = adios2_tread + 1
      else
         ierr = -2
      endif

      end function readADIOS2RecordFile

!     readDerivedTypeADIOS2
!     #################################################################
      subroutine readDerivedTypeADIOS2(engine,varray,step,ierr)

        implicit none

!     Call variables
        type(adios2_engine) :: engine
        !integer*8                :: gh  ! adios group handler (to read data)
        type(var_array),pointer  :: varray
        integer                  :: step
        integer, intent(out)     :: ierr

!     Local variables

        integer      :: ieq,ilom,ihip,jlom,jhip,klom,khip
        integer*8, dimension(4) :: start, readcount ! dimensions are 64bit in adios
        integer      :: vrank, vtype, vtimed, err  ! adios_inq_var() outputs
        integer*8, dimension(10) :: dims ! adios_inq_var() output
        integer,allocatable :: bconds(:)  ! all boundary conditions as one array
        character(20):: vvar      ! /var/v<idx> 
        character(20):: vname     ! /name/v<idx> 
        character(len(varray%array_var(1)%descr))    :: desc      ! name of variable read from file
        integer*8    :: n ! bytes to read

        integer      :: iloghost, ihighost
        integer      :: jloghost, jhighost
        integer      :: kloghost, khighost

        integer      :: xsize, ysize, zsize
        integer      :: xoffset, yoffset, zoffset
        type(adios2_variable) :: var

        logical      :: firstread

        character(len=100), save :: desc_sv(20)
        integer, save :: bconds_sv(120)
        
!     Begin program

        ierr = 0

        firstread = (step == 0)

!       Set offsets to exclude ghosts

        iloghost = 0          ! gcw
        jloghost = 0          ! gcw
        kloghost = 0          ! gcw
        ihighost = 0          ! gcw
        jhighost = 0          ! gcw
        khighost = 0          ! gcw

        call fromGlobalToLocalLimits(gv%gparams,1,ilomg,jlomg,klomg,ilom,jlom,klom)
        call fromGlobalToLocalLimits(gv%gparams,1,ihipg,jhipg,khipg,ihip,jhip,khip)

        if (adios2_debug) &
            write (*,"(a,i0,a,i4,i4,i4,i4,i4,i4)") " ADIOS2 INFO: rank=",my_rank, &
            " read: local indices: ",ilom,ihip,jlom,jhip,klom,khip

        ! calculate size actually written for each dimension
        xsize = ihip+ihighost-ilom-iloghost+1  
        ysize = jhip+jhighost-jlom-jloghost+1  
        zsize = khip+khighost-klom-kloghost+1  

        ! calculate offsets in the global array
        xoffset = ilomg + iloghost  
        yoffset = jlomg + jloghost  
        zoffset = klomg + kloghost  

        start = (/ilomg+iloghost,jlomg+jloghost,klomg+kloghost,step/)
        readcount = (/ ihip-ihighost - ilom-iloghost + 1, &
                       jhip-jhighost - jlom-jloghost + 1, &
                       khip-khighost - klom-kloghost + 1, &
                       1 /)
        n = readcount(1)*readcount(2)*readcount(3)*8
        if (adios2_debug) &
            write (*,"(a,i0,a,i0,a,4i4,a,4i4)") " ADIOS2 INFO: rank=",my_rank, &
            " read vars n=",n," start= ", start," count= ",readcount

        call adios2_get(engine,"nvar",adios2_nvar, adios2_mode_sync,ierr)
        varray%nvar = adios2_nvar

        do ieq=1,adios2_nvar
         if (firstread) then
            desc = "" ; desc_sv(ieq) = ""
            ! read in name of Nth variable
            write (vname, '("/name/v",i0)') ieq
            call adios2_get (engine, vname, desc, adios2_mode_sync, ierr)
            varray%array_var(ieq)%descr = desc
            desc_sv(ieq) = trim(desc)
!!$            if (my_rank == 0) write (*,*) "fr desc",ieq,trim(desc)," ",desc_sv(ieq)
            if (adios2_debug.and.my_rank==0) &
                 write (*,"(a,a,a,a,a,i0)") " ADIOS2 INFO: first read ",trim(vname), &
                 ": [",trim(desc),"]"
         else
!!$            if (my_rank == 0) write (*,*) "aft desc",ieq,trim(desc_sv(ieq))
            varray%array_var(ieq)%descr = trim(desc_sv(ieq))
         endif

         ! read in data of Nth variable
         write (vvar,  '("/var/v",i0)')  ieq
         call adios2_inquire_variable(var, aio, vvar, ierr)
         call adios2_set_selection(var, 3, &
                 int((/ xoffset,yoffset,zoffset /), kind=8), &
                 int((/ xsize,ysize,zsize /), kind=8), ierr)
         call adios2_get(engine, var, &
                            varray%array_var(ieq) &
                               %array(ilom+iloghost:ihip+ihighost, &
                                      jlom+jloghost:jhip+jhighost, &
                                      klom+kloghost:khip+khighost), &
                            ierr)
         if (.not. ierr.eq.0) then
            write (*,"(a,a,i0)") "ERROR: could not read ", trim(vvar), ierr
            ierr = -1
         endif

         if (adios2_debug.and.my_rank==0) write (*,"(a,a,a,i0)") " ADIOS2 INFO: read ",trim(vvar)
       enddo

       if (firstread) then
          allocate(bconds(6*adios2_nvar))
          call adios2_get (engine, "bconds", bconds, adios2_mode_sync, ierr)
          if (.not.ierr.eq.0) then
             write (*,"(a,i0,a,i0)") "ERROR: could not read bconds err=", ierr
             ierr = -1
          endif
          bconds_sv = 0
          bconds_sv(1:size(bconds)) = bconds
          do ieq=1,adios2_nvar
             varray%array_var(ieq)%bconds = bconds( (ieq-1)*6+1:ieq*6 ) 
             if (adios2_debug.and.my_rank==0) &
                  write (*,"(a,6i3)") " ADIOS2 INFO: read bcond=",bconds( (ieq-1)*6+1:ieq*6 )
          enddo
          deallocate(bconds)
       else
          do ieq=1,adios2_nvar
!!$             if (my_rank == 0) write (*,*) "aft BCs",ieq,bconds_sv( (ieq-1)*6+1:ieq*6 )
             varray%array_var(ieq)%bconds = bconds_sv( (ieq-1)*6+1:ieq*6 )
          enddo
       endif

!     End program

      end subroutine readDerivedTypeADIOS2

      subroutine adios2_check_err(errno, msg)
        implicit none
        integer :: errno
        character(*) :: msg

        if (errno.ne.0) then
           write (*,'(a,a,a,i3,a,i4,a)') ' ADIOS2 ERROR: ',trim(msg),' (errno=',errno,', rank=',my_rank,')'
        endif

        if (pmax(dble(errno))/=0d0)  call pstop("","")

      end subroutine adios2_check_err

      subroutine adios2_logging(msg)
        implicit none
        character(*) :: msg

        if (adios2_debug.and.(my_rank.eq.0)) then
           write (*,*) 'ADIOS2 INFO: ',trim(msg)
        endif
      end subroutine adios2_logging
#endif
      end module ADIOS2_mod

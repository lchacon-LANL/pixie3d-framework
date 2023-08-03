! module ADIOS2_io
! ######################################################################
      module ADIOS2_mod

#if defined(ADIOS2)
        use io
        
        use adios2

        use variables

        implicit none

        type(adios2_adios)  :: adios2obj
        type(adios2_engine) :: engine, rengine
        type(adios2_io)     :: aio

        !!! Used to dump main record file (to be postprocessed)
        character(1024) :: recordfile=ADIOS2_FILE

        integer :: adios2_world_comm     ! MPI communicator

        integer :: adios2_err            ! error handler

        logical,private :: adios2_debug=.false.

      contains

!     init_ADIOS2_IO
!     #################################################################
      function init_ADIOS2_IO(obj,io,varray) result(ierr)

!     -----------------------------------------------------------------
!     Initializes ADIOS2 IO according to adios_config.xml
!     -----------------------------------------------------------------

      implicit none

!     Global variables

      integer :: ierr
      type(adios2_adios) :: obj
      type(adios2_io) :: io
      type(var_array),pointer :: varray

!     Local variables

      type(adios2_variable) :: var
      integer      :: err

!     Begin program

!     Initialize ADIOS2 for componentized I/O

      if (.not.is_file('adios_config.xml')) then
         call pstop('init_ADIOS2_IO','ADIOS2 config file not present')
      else
         call MPI_Comm_dup (MPI_COMM_WORLD,adios2_world_comm,ierr)
         if (.not.obj%valid) then
            call adios2_logging('init')
            call adios2_init(obj,'adios_config.xml', adios2_world_comm,.true.,ierr)
            call adios2_check_err(ierr, 'Problem in init')

            call adios2_declare_io (io, obj, "record", err)
!!$            call adios2_set_engine (io, "BP4", err)
            call adios2_define_variable (var, io, "time"  , adios2_type_dp, err)
            call adios2_define_variable (var, io, "itime" , adios2_type_integer4, err)
            call adios2_define_variable (var, io, "dt"    , adios2_type_dp, err)
            call adios2_define_variable (var, io, "gammat", adios2_type_dp, err)
            call defineDerivedTypeADIOS2 (io, varray)
         endif
      endif

      end function init_ADIOS2_IO

!     destroy_ADIOS2_IO
!     #################################################################
      function destroy_ADIOS2_IO(obj,eng) result(ierr)

!     -----------------------------------------------------------------
!     Initializes ADIOS IO according to adios_config.xml
!     -----------------------------------------------------------------

      implicit none

!     Global variables

      integer :: ierr
      type(adios2_adios) :: obj
      type(adios2_engine) :: eng

      call adios2_logging('Terminating IO')

      call adios2_logging('write close')
      call adios2_close(eng, ierr)
      call adios2_check_err(ierr, 'Problem in write close')
      call adios2_logging('Finalize ADIOS2')
      call adios2_finalize(obj, ierr)

      call MPI_Comm_free(adios2_world_comm,ierr)

      end function destroy_ADIOS2_IO

!     defineDerivedTypeADIOS2
!     #################################################################
      subroutine defineDerivedTypeADIOS2(aio, varray)

        implicit none

!     Call variables
        type(adios2_io)          :: aio
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
      subroutine writeDerivedTypeADIOS2(eng, varray)

        implicit none

!     Call variables
        type(adios2_engine)      :: eng
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

        character(200)   :: msg
        
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

        if (my_rank == 0) call adios2_put (eng, "nvar", varray%nvar, err)

        do ieq=1,varray%nvar
            write (vvar,  '("/var/v",i0)')  ieq
            ! write /var/v<ieq>
            call adios2_put (eng, vvar, &
                  varray%array_var(ieq)%array( &
                        ilom+xloghost:ihip+xhighost, &
                        jlom+yloghost:jhip+yhighost, &
                        klom+zloghost:khip+zhighost), &
                        adios2_mode_sync, err)

            if (my_rank == 0) then
                ! write /name/v<ieq>
                write (vname, '("/name/v",i0)') ieq
                call adios2_put (eng, vname, trim(varray%array_var(ieq)%descr)//char(0), &
                                 adios2_mode_sync, err)
                write (msg, '(" attr write ",a,": ",a)') trim(vname), &
                         trim(varray%array_var(ieq)%descr)
                call adios2_logging(msg)
            endif
        enddo

        if (my_rank == 0) then
           do ieq=1,varray%nvar
              bconds( (ieq-1)*6+1:ieq*6 ) = varray%array_var(ieq)%bconds(:)
              write (msg,'(" attr write bconds=(",6i3")")') &
                      bconds( (ieq-1)*6+1:ieq*6 )
              call adios2_logging(msg)
           enddo
           call adios2_put (eng, "bconds", bconds, adios2_mode_sync, err)
        endif
        
      end subroutine writeDerivedTypeADIOS2

!     readDerivedTypeADIOS2
!     #################################################################
      subroutine readDerivedTypeADIOS2(io,eng,varray,step,ierr)

        implicit none

!     Call variables
        type(adios2_io)          :: io
        type(adios2_engine)      :: eng
        type(var_array),pointer  :: varray
        integer                  :: step
        integer, intent(out)     :: ierr

!     Local variables

        integer      :: ieq,ilom,ihip,jlom,jhip,klom,khip
        integer*8, dimension(4) :: start, readcount ! dimensions are 64bit in adios
        integer,allocatable :: bconds(:)  ! all boundary conditions as one array
        character(20):: vvar      ! /var/v<idx> 
        character(20):: vname     ! /name/v<idx> 
        character(len(varray%array_var(1)%descr)) :: desc      ! name of variable read from file
        integer*8    :: n ! bytes to read

        integer      :: iloghost, ihighost
        integer      :: jloghost, jhighost
        integer      :: kloghost, khighost

        integer      :: xsize, ysize, zsize
        integer      :: xoffset, yoffset, zoffset
        type(adios2_variable) :: var

        logical      :: firstread

        integer      :: adios2_nvar

        character(200) :: msg

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

        if (firstread) then
           call adios2_get(eng,"nvar",adios2_nvar, adios2_mode_sync,ierr)
           if (adios2_nvar /= varray%nvar) then
              write (*,*) adios2_nvar,varray%nvar
              call pstop("readADIOS2Derivedtype" &
                        ,"Incorrect # of variables in ADIOS2 file")
           endif
        endif

        if (adios2_debug) write (*,"(a,i0,a,i0)") " ADIOS2 INFO: rank=",my_rank, &
            " read # vars=",varray%nvar

        do ieq=1,varray%nvar
         if (firstread) then
            desc = ""
            ! read in name of Nth variable
            write (vname, '("/name/v",i0)') ieq
            call adios2_get (eng, vname, desc, adios2_mode_sync, ierr)
            varray%array_var(ieq)%descr = desc
            write(msg,"(a,a,a,a,a,i0)") "first read ",trim(vname), ": [",trim(desc),"]"
            call adios2_logging(msg)
         endif

         ! read in data of Nth variable
         write (vvar,  '("/var/v",i0)')  ieq
         call adios2_inquire_variable(var, io, vvar, ierr)
         call adios2_set_selection(var, 3, &
                 int((/ xoffset,yoffset,zoffset /), kind=8), &
                 int((/ xsize,ysize,zsize /), kind=8), ierr)
         call adios2_get(eng, var, &
                            varray%array_var(ieq) &
                               %array(ilom+iloghost:ihip+ihighost, &
                                      jlom+jloghost:jhip+jhighost, &
                                      klom+kloghost:khip+khighost), &
                            ierr)

         write (msg,"(a,a,i0)") "could not read ", trim(vvar), ierr
         call adios2_check_err(ierr,msg)

         write (msg,"(a,a)") "read ",trim(vvar)
         call adios2_logging(msg)
       enddo

       allocate(bconds(6*varray%nvar))
       call adios2_get (eng, "bconds", bconds, adios2_mode_sync, ierr)
       call adios2_check_err(ierr,"could not read bconds")

       do ieq=1,varray%nvar
          varray%array_var(ieq)%bconds = bconds( (ieq-1)*6+1:ieq*6 )
          write (msg,"(a,6i3)") "read bcond=",bconds( (ieq-1)*6+1:ieq*6 )
          call adios2_logging(msg)
       enddo
       deallocate(bconds)

!     End program

      end subroutine readDerivedTypeADIOS2

!     openADIOS2RestartFileForRead
!     #################################################################
      function openADIOS2RestartFileForRead(file) result(ufile)

!     -----------------------------------------------------------------
!     Opens ADIOS restart file.
!     -----------------------------------------------------------------

      implicit none

!     Global variables

      integer :: ufile
      character(*),optional :: file
      integer :: ierr

!     Local variables

      character(len=len(recordfile)) :: rfile

!     Begin program

      ufile = find_unit(12345)

      if (PRESENT(file)) then
         rfile = file
      else
         rfile = recordfile
      endif

      ierr = openADIOS2FileForRead(adios2obj,aio,rengine,rfile)

      if (ierr /= 0) then
         call pstop('openADIOS2RestartFileForRead' &
                   ,'Error reading ADIOS restart file')
      endif

      end function openADIOS2RestartFileForRead

!     closeADIOS2RestartFileForRead
!     #################################################################
      function closeADIOS2RestartFileForRead() result(ierr)

        implicit none

        integer :: ierr

        ierr = closeADIOS2FileForRead(adios2obj,rengine)

        if (ierr /= 0) then
           call pstop('closeADIOS2RestartFileForRead' &
                     ,'Error closing ADIOS restart file')
        endif

      end function closeADIOS2RestartFileForRead

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

      integer :: ierr
      logical :: rinit

!     Begin program

      if (PRESENT(init)) then
         rinit = init
      else
         rinit = .false.
      endif

      ierr = writeADIOS2File(aio,engine,file,itime,time,varray,rinit &
                            ,dt=dt,gammat=gammat)

      if (ierr /= 0) then
         call pstop('writeADIOS2RecordFile','Error writing ADIOS file')
      endif

      end subroutine writeADIOS2RecordFile
    
!     readADIOS2RecordFile
!     #################################################################
      function readADIOS2RecordFile(itime,time,dt,gammat,varray) &
               result(ierr)

!     -----------------------------------------------------------------
!     Reads record file
!     -----------------------------------------------------------------

      implicit none

!     Call variables

      integer    :: ierr,itime
      real(8)    :: time,dt,gammat

      type(var_array),pointer :: varray

!     Local variables

!     Begin program

      ierr = readADIOS2File(aio,rengine,itime,time,varray &
                           ,dt=dt,gammat=gammat)

      end function readADIOS2RecordFile

!     openADIOS2FileForRead
!     #################################################################
      function openADIOS2FileForRead(obj,io,eng,file) result(ierr)

        implicit none

!     Call variables

        type(adios2_adios) :: obj
        type(adios2_engine) :: eng
        type(adios2_io) :: io

        integer :: ierr
        character(*) :: file

!     Local variables

        real(8) :: tmp

!     Begin program

        call adios2_logging('read-only open')

        call MPI_Comm_dup (MPI_COMM_WORLD, adios2_world_comm, ierr)
        call adios2_logging('read init')
        call adios2_init(obj,'adios_config.xml',adios2_world_comm,.true.,ierr)
        call adios2_check_err(ierr, 'Problem in read init')
        
        !Open ADIOS file
        call adios2_logging('declare io')
        call adios2_declare_io(io, obj, 'record.read', ierr)

        call adios2_logging('open file')
        call adios2_open(eng,io,trim(file),adios2_mode_read,adios2_world_comm,ierr)
        call adios2_check_err(ierr, 'Problem in read-only open')
        
        call adios2_get(eng, 'time', tmp, ierr)

      end function openADIOS2FileForRead

!     closeADIOS2FileForRead
!     #################################################################
      function closeADIOS2FileForRead(obj,eng) result(ierr)

        implicit none
        
        type(adios2_adios)  :: obj
        type(adios2_engine) :: eng

        integer :: ierr

        call adios2_logging('Terminating IO')

        call adios2_logging('read close')
        call adios2_close(eng, ierr)
        call adios2_check_err(ierr, 'Problem in read close')

        call adios2_logging('Finalize ADIOS2')
        call adios2_finalize(obj, ierr)

        call MPI_Comm_free(adios2_world_comm,ierr)

      end function closeADIOS2FileForRead

!     writeADIOS2File
!     #################################################################
      function writeADIOS2File(io,eng,file,itime,time,varray,init,dt,gammat) &
               result(err)

!     -----------------------------------------------------------------
!     Writes record file
!     -----------------------------------------------------------------

      implicit none

!     Call variables

      type(adios2_io)     :: io
      type(adios2_engine) :: eng
      integer             :: itime,err
      real(8)             :: time
      real(8),optional    :: dt,gammat
      character(*)        :: file
      logical             :: init
      type(var_array),pointer :: varray

!     Local variables

      character(2) :: mode='a'//char(0)
      integer      :: adios2_mode
      logical      :: addl_write  ! true: write names and bconds

!     Begin program
        
!     Create/Append adios file

      !!jyc: init is true only with no restart. With restart, init is false.
      if (init) then
         mode = 'w'//char(0)
         adios2_mode = adios2_mode_write
      else 
         mode = 'a'//char(0)
         adios2_mode = adios2_mode_append
      endif

      call adios2_logging('write access mode: '//trim(mode))
      call adios2_logging('write open')

      if (.not.eng%valid) then
         call adios2_open(eng,io,file,adios2_mode,adios2_world_comm,err)
         call adios2_check_err(err,'Could not open file for writing')
      else
         call adios2_logging('not first write call: no attrs written')
      endif

      call adios2_logging('begin write step')
      call adios2_begin_step(eng, adios2_step_mode_append, err)
      call adios2_check_err(err, 'Problem in begin write step')

      if (my_rank == 0) call adios2_put(eng,"time"  ,time  ,err)
      if (my_rank == 0) call adios2_put(eng,"itime" ,itime ,err)
      if (my_rank == 0.and.PRESENT(dt    )) &
           call adios2_put(eng,"dt"    ,dt    ,err)
      if (my_rank == 0.and.PRESENT(gammat)) &
           call adios2_put(eng,"gammat",gammat,err)

      call writeDerivedTypeADIOS2(eng, varray)
      
      !End adios step
      call adios2_logging('end write step')
      call adios2_end_step(eng, err)
      call adios2_check_err(err, 'Problem end write step')

      end function writeADIOS2File

!     readADIOS2File
!     #################################################################
      function readADIOS2File(io,eng,itime,time,varray,dt,gammat) &
               result(ierr)

!     -----------------------------------------------------------------
!     Reads ADIOS2 file
!     -----------------------------------------------------------------

      implicit none

!     Call variables

      type(adios2_io)     :: io
      type(adios2_engine) :: eng
      integer    :: ierr,itime
      real(8)    :: time
      real(8),optional :: dt,gammat
      
      type(var_array),pointer :: varray

!     Local variables

      integer   :: istatus

!     Begin program

      ierr = 0

      call adios2_logging('begin read step')
      call adios2_begin_step(eng, adios2_step_mode_read, 0.0, istatus, ierr)
      call adios2_check_err(ierr, 'Problem in begin read step')
      
      if ((istatus.eq.0).and.(ierr.eq.0)) then
         call adios2_get(eng,"time",time,ierr)
         call adios2_get(eng,"itime",itime,ierr)
         if (PRESENT(dt    )) call adios2_get(eng,"dt",dt,ierr)
         if (PRESENT(gammat)) call adios2_get(eng,"gammat",gammat,ierr)
         if (adios2_debug.and.(my_rank.eq.0)) then
            write (*,"(a,i5)") " ADIOS2 INFO: read time=", itime
         endif
         call readDerivedTypeADIOS2(io,eng,varray,itime,ierr)
         call adios2_end_step(eng, ierr)
         call adios2_logging('end read step')
      else
         ierr = -2
      endif

      end function readADIOS2File

!     adios2_check_err
!     #################################################################
      subroutine adios2_check_err(errno, msg)
        implicit none
        integer :: errno
        character(*) :: msg

        character(1000) :: msg2

        if (errno.ne.0) then
           write (msg2,'(a,a,a,i3,a,i4,a)') ' ADIOS2 ERROR: ',trim(msg),&
                ' (errno=',errno,', rank=',my_rank,')'
        endif

        if (ipmax(errno)/=0)  call pstop("",msg2)

      end subroutine adios2_check_err

!     adios2_logging
!     #################################################################
      subroutine adios2_logging(msg)
        implicit none
        character(*) :: msg

        if (adios2_debug.and.(my_rank.eq.0)) then
           write (*,*) 'ADIOS2 INFO: ',trim(msg)
        endif
      end subroutine adios2_logging
#endif
      end module ADIOS2_mod

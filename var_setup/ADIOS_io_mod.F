c module ADIOS_io
c ######################################################################
      module ADIOS_io

#if defined(adios)
        use io

        use variables

        !!! ADIOS variables
        !!! Used to dump main record file (to be postprocessed)
        character*(1024) :: recordfile='record.bp'

        ! ADIOS needs a variable containing the MPI_WORLD_COMM communicator
        integer :: adios_world_comm

        !!! variables for reading 
        integer*8 :: adios_fh, adios_gh   ! file handler and group handler
        integer   :: adios_varcount       ! number of variables in file (~32)
        integer   :: adios_attrcount      ! number of attributes in file (~16)
        character(30), dimension(50) :: adios_varlist  ! list of all vars 
        character(30), dimension(50) :: adios_attrlist  ! list of all attributes 
        integer   :: adios_nvar           ! number of /var/ variables in file (~8)
        character(30), dimension(50) :: adios_varnames  ! list of /var/ vars 
        integer   :: adios_tstart, adios_tend, adios_ntsteps ! timesteps available
        integer   :: adios_tread          ! 0 ... ntsteps-1
        integer   :: adios_err            ! error handler
        integer   :: adios_attrtype,adios_attrsize,adios_vartype
     .              ,adios_ndim,adios_tdim
        integer*8,dimension(0:3) :: adios_vardims

        logical,private :: adios_debug=.false.
     $                    ,adios_append_recordfile=.false.

#if !defined(ADIOS_BUFFER_MB)
        !ADIOS IO buffer (for restart, postprocessing; main app set in adios_config.xml)
        integer,private :: ADIOS_BUFFER_MB=100
#endif
      contains

c     setADIOSAppendIOMode
c     ##############################################################
      subroutine setADIOSAppendIOMode

        implicit none

c       Call variables

c       Local variables

        adios_append_recordfile = .true. !Set ADIOS recordfile to append

      end subroutine setADIOSAppendIOMode

c     init_ADIOS_IO
c     #################################################################
      function init_ADIOS_IO() result(ierr)

c     -----------------------------------------------------------------
c     Initializes ADIOS IO according to adios_config.xml
c     -----------------------------------------------------------------

      implicit none

c     Global variables

      integer :: ierr
      
c     Local variables

c     Begin program

c     Initialize ADIOS for componentized I/O

      if (.not.is_file('adios_config.xml')) then
        call pstop('setupRecordFileForWrite'
     .            ,'ADIOS config file not present')
      else
        call MPI_Comm_dup (MPI_COMM_WORLD,adios_world_comm,adios_err)

        if (set_adios_buffer(ADIOS_BUFFER_MB) /= 0) then
          call pstop("init_ADIOS_IO"
     .              ,"Could not setup ADIOS buffer")
        endif

        call ADIOS_init('adios_config.xml',adios_world_comm,adios_err);
      endif

      contains

c     set_adios_buffer
c     #################################################################
      function set_adios_buffer(buffer) result(ierr)

        implicit none

        integer :: buffer,ierr

        character(60) :: text,line

        if (my_rank == 0) then
          text = int2char(buffer)
          line = '    <buffer size-MB=\"'//trim(text)//'\"/>'
          ierr = replace_line('adios_config.xml',"buffer size-MB",line)
        endif

#if defined(petsc) || defined(petsc_c)
        call MPI_Bcast(ierr,1,MPI_INTEGER,0,ADIOS_WORLD_COMM,mpierr)
#endif

      end function set_adios_buffer

      end function init_ADIOS_IO

c     openADIOSRecordFileForRead
c     #################################################################
      subroutine openADIOSRecordFileForRead(ierr,file,nvar,nx,ny,nz)

        implicit none

c     Call variables

        integer,intent(out) :: ierr
        character(*),optional :: file
        integer,optional,intent(out) :: nvar,nx,ny,nz

c     Local variables

        integer :: gcnt, vcnt, acnt, nullidx, ts, lastts,i
        character (len=20), dimension(1)  :: gnamelist  ! list of groups (=1)
        character(len=1024) :: rfile
        logical :: inq_adios_var

c     Begin program

        if (PRESENT(file)) then
          rfile = file
        else
          rfile = recordfile
        endif

        call MPI_Comm_dup (MPI_COMM_WORLD, adios_world_comm, ierr)

        !Open ADIOS file
        call adios_fopen(adios_fh,rfile,adios_world_comm,gcnt,ierr)

        !Inquire ADIOS file
        call adios_inq_file(adios_fh,vcnt,acnt,adios_tstart
     .                     ,adios_ntsteps,gnamelist,ierr)

        adios_tend = adios_tstart + adios_ntsteps - 1

        !Open ADIOS group
        call adios_gopen (adios_fh,adios_gh,gnamelist(1)
     .                   ,adios_varcount,adios_attrcount,ierr)

        !Inquire ADIOS group
        call adios_inq_group(adios_gh,adios_varlist,adios_attrlist
     .                      ,ts,lastts,ierr)

        ! count /var/ variables
        inq_adios_var = (PRESENT(nx).or.PRESENT(ny).or.PRESENT(nz))
        adios_nvar=0
        do i=1,adios_varcount
          if (index(adios_varlist(i),"/var/") == 1) then
            adios_nvar = adios_nvar+1
            !Check problem size
            call adios_inq_var(adios_gh,trim(adios_varlist(i))
     .               ,adios_vartype,adios_ndim,adios_vardims,adios_tdim
     .               ,ierr)
c$$$                  write (*,*) "Var="//adios_varlist(i) //" Dimensions: "
c$$$     .                 , adios_vardims,adios_ndim,adios_tdim 
          endif
        enddo
        
cc        do i=1,adios_attrcount
cc          call adios_inq_attr(adios_gh,adios_attrlist(i),adios_attrtype
cc     .                       ,adios_attrsize,ierr)
cc          write (*,*) adios_attrlist(i),adios_attrtype,adios_attrsize
cc        enddo

        adios_tread = 0  ! next time read first timestep

        if (adios_debug) then
          write (*,"(a,i0,a,a)")
     .       " Proc=",my_rank,": opened ",trim(rfile)

          if (my_rank == 0) then
            write (*,"(a, i5)") " Number of vars: ", adios_varcount
            write (*,*), "Variables:"
            do i=1,adios_varcount
              write (*,"(i5, a, a)") i,")  ", trim(adios_varlist(i))
            enddo
            write (*,"(a, 3i5)") " Dimensions: ", adios_vardims(1:3)
          endif
        endif

        if (PRESENT(nvar)) nvar = adios_nvar
        !C ordering for array dimensions
        if (PRESENT(nx)) nx = adios_vardims(1)-2
        if (PRESENT(ny)) ny = adios_vardims(2)-2
        if (PRESENT(nz)) nz = adios_vardims(3)-2

      end subroutine openADIOSRecordFileForRead

c     closeADIOSRecordFileForRead
c     #################################################################
      subroutine closeADIOSRecordFileForRead(ierr)

        implicit none

        integer, intent(out) :: ierr

        call adios_gclose(adios_gh, ierr)
        call adios_fclose(adios_fh, ierr)

        if (adios_debug) write (*,*), "Closed ", trim(recordfile)

      end subroutine closeADIOSRecordFileForRead

c     openADIOSRestartFileForRead
c     #################################################################
      function openADIOSRestartFileForRead(file) result(ufile)

c     -----------------------------------------------------------------
c     Opens ADIOS restart file. Does NOT use "adios_config.xml" file.
c     -----------------------------------------------------------------

      implicit none

c     Global variables

      integer :: ufile
      character(*),optional :: file

c     Local variables

c     Begin program

cc      ufile = urecord
      ufile = find_unit(12345)

      call MPI_Comm_dup (MPI_COMM_WORLD, adios_world_comm,adios_err)
      call adios_init_noxml (adios_world_comm, adios_err)

      !allocate buffer for ADIOS
      call adios_set_max_buffer_size(ADIOS_BUFFER_MB,adios_err)

      call openADIOSRecordFileForRead(adios_err,file=file
     .                               ,nx=nxd,ny=nyd,nz=nzd)

      if (adios_err /= 0) then
        call pstop('openRestartFileForRead'
     $            ,'Error reading ADIOS restart file')
      endif

      end function openADIOSRestartFileForRead

c     readADIOSRecordFile
c     #################################################################
      function readADIOSRecordFile(unit,itime,time,dt,gammat,varray)
     .         result(ierr)

c     -----------------------------------------------------------------
c     Reads record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ierr,itime,unit
      real(8)    :: time,dt,gammat

      type(var_array),pointer :: varray

c     Local variables

      integer    :: tmp(1),lerr(1),gerr(1)

      !ADIOS variables
      integer*8       :: start(1), readcount(1)  ! dimensions are 64bit in adios
      integer*8       :: m ! actually read bytes
      logical, save   :: firstread = .true.

c     Begin program

      ierr = 0

      !firstread = (adios_tread == 0)

      if (adios_debug)
     .  write (*,"(a,i0,a,i5)"), "rank=",my_rank,
     .       " adios read time=", adios_tread

      if (adios_tread < adios_ntsteps) then
          ! read 1 item from the 1D arrays, from position adios_tread
          start = adios_tread 
          readcount = 1

          if (adios_debug)
     .      write (*,"(a,i0,a,i5,a,i5)"), "rank=",my_rank,
     .           ": read in time.  start=",start,"count=",readcount

          call adios_read_var(adios_gh,"time",start,readcount,time,m)
          start = adios_tread 

          if (adios_debug)
     .      write (*,"(a,i0,a,i5,a,i5)"), "rank=",my_rank,
     .           ": read in itime. start=",start,"count=",readcount

          call adios_read_var(adios_gh,"itime",start,readcount,itime,m)
          start = adios_tread 

          if (adios_debug)
     .      write (*,"(a,i0,a,i5,a,i5)"), "rank=",my_rank,
     .           ": read in dt.   start=",start,"count=",readcount

          call adios_read_var(adios_gh,"dt",start,readcount,dt,m)
          start = adios_tread 

          if (adios_debug)
     .      write (*,"(a,i0,a,i5,a,i5)"), "rank=",my_rank,
     .           ": read in gammat.   start=",start,"count=",readcount

          call adios_read_var(adios_gh,"gammat",start,readcount,gammat
     .                       ,m)

          if (firstread) then
            start = adios_tread 
            call adios_read_var(adios_gh,"nvar",start,readcount,tmp,m)

            varray%nvar = tmp(1)
            ! consistency check
            if (adios_nvar.ne.varray%nvar) then
                write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .                 " ERROR: adios_nvar=",adios_nvar,
     .                 " != varray%nvar=",varray%nvar
            endif
          endif

          if (adios_debug)
     .      write (*,"(a,i0,a,i5,a,f8.2,a,i2)")
     .          "rank=",my_rank
     .         ," read itime=",itime
     .         ," time=",time
     .         ," nvar=",varray%nvar

          call readDerivedTypeADIOS(adios_gh,varray,adios_tread
     .                             ,firstread,ierr)

          if (adios_debug) then
            write (*,*), 'rank=',my_rank,' ADIOS read returns ierr',ierr
            write (*,*)
          endif

          adios_tread = adios_tread + 1
      else  ! EOF
         ierr = -2
      endif

      firstread = .false.

      end function readADIOSRecordFile

c     writeADIOSRecordFile
c     #################################################################
      subroutine writeADIOSRecordFile(file,itime,time,dt,gammat,varray
     .                               ,init)

c     -----------------------------------------------------------------
c     Writes record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: itime
      real(8)    :: time,dt,gammat
      type(var_array),pointer:: varray
      character(*) :: file
      logical,optional :: init

c     Local variables

      logical :: rinit

      integer*8    :: handle, totalsize, groupsize
      integer      :: err
      character(2) :: mode='a'//char(0)
      logical      :: addl_write  ! true: write names and bconds

      integer      :: xsize, ysize, zsize
      integer      :: xloghost, xhighost
      integer      :: yloghost, yhighost
      integer      :: zloghost, zhighost

c     Begin program

      if (PRESENT(init)) then
         rinit = init
      else
         rinit = .not.adios_append_recordfile
      endif

c     Offsets to include ghosts
c     1. internal processes: no ghost cells written out

      xloghost = 1
      yloghost = 1
      zloghost = 1
      xhighost = -1
      yhighost = -1
      zhighost = -1
      
c     2. processes that has some face 

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

      !calculate size actually written for each dimension
      xsize = ihip+xhighost-ilom-xloghost+1  
      ysize = jhip+yhighost-jlom-yloghost+1  
      zsize = khip+zhighost-klom-zloghost+1  

c     Create/Append adios file

      addl_write=.false.
cc      if (record_open_num.eq.0) then
      if (rinit) then
          mode = 'w'//char(0)
          if (my_rank == 0) then
c             Process 0 (arbitrary choice) 
c             in the first timestep writes 
c             the boundary conditions and variable names
              addl_write = .true.
          endif
          adios_append_recordfile = .true.
      else 
          mode = 'a'//char(0)
      endif

cc      record_open_num = record_open_num+1
      groupsize = 10 * 4 +     ! 12 integers
     .  2 * 8 +                ! 2 doubles + nvar*double 3D arrays
     .  varray%nvar * 8 * xsize*ysize*zsize

      if (addl_write) then
         groupsize = groupsize + 4 + 4 + 4 + ! namelen
     .     varray%nvar * 6 * 4 +  ! bconds
     .     varray%nvar * 20       ! descr (var names)
      endif

      if (adios_debug.and.my_rank==0) then
        write (*,*) ' ADIOS file open...'
      endif

      call adios_open(handle,'record',file,mode
     .               ,adios_world_comm,err)

      if (err /= 0) then
        write (*,*) 'Problem in writeRecordFile'
        write (*,*) 'rank=',my_rank,'  ERROR in "adios_open"'
        stop
      endif

      if (adios_debug.and.my_rank==0) then
        write (*,*) 'ADIOS file access mode=',mode
        write (*,*) 'ADIOS open: adios_file=',trim(file)
      endif

      if (adios_debug) then
        write (*,111) my_rank, ihip, ilom, jhip, jlom, khip, klom
 111    format ('rank= ',i0,' ihip=',i0,' ilom=',i0,' jhip='
     .         ,i0,' jlom=',i0,' khip=',i0,' klom=',i0)
        write (*,*) 'rank=',my_rank,' ADIOS open: handle=',handle
      endif

c     Inquire IO sizes

      if (adios_debug.and.my_rank==0) then
        write (*,*) ' ADIOS file inquire...'
      endif

      call adios_group_size(handle,groupsize,totalsize,err)

      if (err /= 0) then
        write (*,*) 'Problem in writeRecordFile'
        write (*,*) 'rank=',my_rank,'  ERROR in "adios_group_size"'
        stop
      endif

      if (adios_debug) then
        write (*,*) 'rank=',my_rank,' ADIOS group_size=',groupsize
        write (*,*) 'rank=',my_rank,' ADIOS total size=',totalsize
      endif

c     Write all variables to buffer

      if (adios_debug.and.my_rank==0) then
        write (*,*) ' ADIOS file write...'
      endif

      if (my_rank == 0) then
        call adios_write(handle,"time"  ,time  ,err)
        call adios_write(handle,"itime" ,itime ,err)
        call adios_write(handle,"dt"    ,dt    ,err)
        call adios_write(handle,"gammat",gammat,err)
      endif

      call writeDerivedTypeADIOS(handle,varray,addl_write)

c     Close adios file which starts I/O writing from the buffer

      if (adios_debug.and.my_rank==0) then
        write (*,*) ' ADIOS file close...'
      endif

      call adios_close (handle, err)

      if (err /= 0) then
        write (*,*) 'Problem in writeRecordFile'
        write (*,*) 'rank=',my_rank,'  ERROR in "adios_close"'
        stop
      endif

c     End

      end subroutine writeADIOSRecordFile

c     writeDerivedTypeADIOS
c     #################################################################
      subroutine writeDerivedTypeADIOS(handle, varray, addl_write)

        implicit none

c     Call variables

        integer*8                :: handle
        type(var_array),pointer  :: varray
        logical                  :: addl_write

c     Local variables

        integer          :: ieq, err, 
     .                      xoffset, yoffset, zoffset,
     .                      xsize, ysize, zsize

        integer          :: xloghost, xhighost
        integer          :: yloghost, yhighost
        integer          :: zloghost, zhighost
        integer,dimension(6*varray%nvar) :: bconds  ! all boundary conditions as one array
        character(20)    :: vvar      ! /var/v<idx> 
        character(20)    :: vname     ! /name/v<idx> 

c     Begin program

c   global size of the output array = nxd+2 x nyd+2 x nzd+2 with ghost cells on the faces

c   ihipg = nxd+1    -> x=1 face
c   ihipg = nxd+1    -> x=1 face
c   j... -> y face
c   k... -> z face
c
c   global position of an internal processor in the 3D mesh
c     from  ilog = ilomg+1  to ihig = ihipg-1
c     from  jlog = jlomg+1  to jhig = jhipg-1
c     from  klog = klomg+1  to khig = khipg-1
c   global position for boundary proc
c     if (ilomg == 0) ilog = ilomg
c     if (ihipg == nxd+1) ihig = ihipg
c
c   local sizes: ilom:ihip   but ilom+1:ihip-1 for internal processor
c   local indices start from 0 for ghost cell

c     Offsets to include ghosts
c     1. internal processes: no ghost cells written out
        xloghost = 1
        yloghost = 1
        zloghost = 1
        xhighost = -1
        yhighost = -1
        zhighost = -1
      
c     2. processes that has some face 
c       x=0 face
        if (ilomg == 0) xloghost = 0

c       x=1 face
        if (ihipg == nxd+1) xhighost = 0

c       y=0 face
        if (jlomg == 0) yloghost = 0

c       y=1 face
        if (jhipg == nyd+1) yhighost = 0

c       z=0 face
        if (klomg == 0) zloghost = 0

c       z=1 face
        if (khipg == nzd+1) zhighost = 0

c       calculate offsets in the global array
        xoffset = ilomg + xloghost  
        yoffset = jlomg + yloghost  
        zoffset = klomg + zloghost  


c       write auxiliary variables (adios writer needs it),
c          will not be accessible from file (everybody writes these vars)
        call adios_write (handle, "nvar", varray%nvar, err)
        call adios_write (handle, "nxd+2", nxd+2, err)
        call adios_write (handle, "nyd+2", nyd+2, err)
        call adios_write (handle, "nzd+2", nzd+2, err)
        call adios_write (handle, "xoffset", ilomg+xloghost, err)
        call adios_write (handle, "yoffset", jlomg+yloghost, err)
        call adios_write (handle, "zoffset", klomg+zloghost, err)
        call adios_write (handle, "xsize", 
     .                    ihip+xhighost-ilom-xloghost+1, err)
        call adios_write (handle, "ysize", 
     .                    jhip+yhighost-jlom-yloghost+1, err)
        call adios_write (handle, "zsize", 
     .                    khip+zhighost-klom-zloghost+1, err)
        if (addl_write) then
            call adios_write (handle, "namelen", 
     .                        len(varray%array_var(1)%descr), err)
        endif

        if (adios_debug) then
          print"('r',i1,' offset=(',i3,',',i3,',',i3,')')", my_rank,
     .           ilomg+xloghost, jlomg+yloghost, klomg+zloghost
          print"('r',i1,' size=(',i3,',',i3,',',i3,')')", my_rank,
     .           ihip+xhighost-ilom-xloghost+1, 
     .           jhip+yhighost-jlom-yloghost+1,
     .           khip+zhighost-klom-zloghost+1
          print"('r',i1,' i=',i3,':',i3,',',i3
     .          ,':',i3,',',i3,':',i3,'')",
     .           my_rank,
     .           ilom+xloghost, ihip+xhighost,
     .           jlom+yloghost, jhip+yhighost,
     .           klom+zloghost, khip+zhighost
        endif

        do ieq=1,varray%nvar
            write (vvar,  '("/var/v",i0)')  ieq
            ! write /var/v<ieq>
            call adios_write (handle, vvar,
     .          varray%array_var(ieq)%array(
     .               ilom+xloghost:ihip+xhighost,
     .               jlom+yloghost:jhip+yhighost,
     .               klom+zloghost:khip+zhighost), 
     .          err)

            if (addl_write) then
                ! write /name/v<ieq>
                write (vname, '("/name/v",i0)') ieq
                call adios_write (handle, vname,
     .              varray%array_var(ieq)%descr, err)
            endif
        enddo

c     write all boundary conditions into one array buffer
        if (addl_write) then
            do ieq=1,varray%nvar
                bconds( (ieq-1)*6+1:ieq*6 ) = 
     .                  varray%array_var(ieq)%bconds(:)
            enddo
            call adios_write (handle, "nvar*6", varray%nvar*6, err)
            call adios_write (handle, "bconds", bconds, err)
        endif

c     End program

      end subroutine writeDerivedTypeADIOS

c     readDerivedTypeADIOS
c     #################################################################
      subroutine readDerivedTypeADIOS(gh,varray,step,firstread,ierr)

        implicit none

c     Call variables
        integer*8                :: gh  ! adios group handler (to read data)
        type(var_array),pointer  :: varray
        integer                  :: step
        logical                  :: firstread
        integer, intent(out)     :: ierr

c     Local variables

        integer          :: ieq,ilom,ihip,jlom,jhip,klom,khip
        integer*8, dimension(4) :: start, readcount ! dimensions are 64bit in adios
        integer*8        :: zero, nlen ! 64bit start/count for reading the description
        integer          :: vrank, vtype, vtimed, err  ! adios_inq_var() outputs
        integer*8, dimension(10) :: dims ! adios_inq_var() output
        integer,dimension(6*adios_nvar) :: bconds  ! all boundary conditions as one array
        character(20)    :: vvar      ! /var/v<idx> 
        character(20)    :: vname     ! /name/v<idx> 
        character(len(varray%array_var(1)%descr))    :: desc      ! name of variable read from file
        integer*8        :: n,m ! bytes to read and actually read

        integer          :: iloghost, ihighost
        integer          :: jloghost, jhighost
        integer          :: kloghost, khighost

c     Begin program

          ierr = 0

c     Set offsets to include ghosts
c     1. internal processes: no ghost cells will be read in
          iloghost = 0          ! gcw
          jloghost = 0          ! gcw
          kloghost = 0          ! gcw
          ihighost = 0          ! gcw
          jhighost = 0          ! gcw
          khighost = 0          ! gcw
      
c     2. processes that has some face have ghost cells 
c       x=0 face
          if (ilomg == 0) iloghost = 0

c       x=1 face
          if (ihipg == nxd+1) ihighost = 0

c       y=0 face
          if (jlomg == 0) jloghost = 0

c       y=1 face
          if (jhipg == nyd+1) jhighost = 0

c       z=0 face
          if (klomg == 0) kloghost = 0

c       z=1 face
          if (khipg == nzd+1) khighost = 0


c        write (*,"(a,i1,a,i4,i4,i4,i4,i4,i4)") "rank=",my_rank,
c     .      " global indices again: ",
c     .      ilomg,ihipg,jlomg,jhipg,klomg,khipg

          call fromGlobalToLocalLimits(gv%gparams,1
     .         ,ilomg,jlomg,klomg,ilom,jlom,klom)
          call fromGlobalToLocalLimits(gv%gparams,1
     .         ,ihipg,jhipg,khipg,ihip,jhip,khip)

          if (adios_debug)
     .      write (*,"(a,i0,a,i4,i4,i4,i4,i4,i4)") "rank=",my_rank,
     .      " read: local indices: ",ilom,ihip,jlom,jhip,klom,khip

          start = (/ilomg+iloghost,jlomg+jloghost,klomg+kloghost,step/)
          readcount = (/ ihip-ihighost - ilom-iloghost + 1,
     .                   jhip-jhighost - jlom-jloghost + 1,
     .                   khip-khighost - klom-kloghost + 1,
     .                   1 /)
          n = readcount(1)*readcount(2)*readcount(3)*8
          if (adios_debug)
     .      write (*,"(a,i0,a,i0,a,4i4,a,4i4)") "rank=",my_rank,
     .      " read vars n=",n," start= ", start," count= ",readcount

          do ieq=1,adios_nvar
            if (firstread) then
                zero = 0
                nlen = len(varray%array_var(1)%descr)
                ! read in name of Nth variable
                write (vname, '("/name/v",i0)') ieq
                !write (*,"(a,a,a,i0)") "read ",trim(vname),": len=", nlen
                call adios_read_var (gh, vname, zero, nlen, desc, m)
                varray%array_var(ieq)%descr = desc
                if (adios_debug)
     .            write (*,"(a,a,a,a,a,i0)") "read ",trim(vname),
     .             ": [",trim(desc),"] m=",m
            endif

            ! read in data of Nth variable
            write (vvar,  '("/var/v",i0)')  ieq
            call adios_read_var (gh, vvar, 
     .                       start, readcount,
     .                       varray%array_var(ieq)
     .                          %array(ilom+iloghost:ihip+ihighost,
     .                                 jlom+jloghost:jhip+jhighost,
     .                                 klom+kloghost:khip+khighost),
     .                       m)
            n = readcount(1)*readcount(2)*readcount(3)*readcount(4)*8
            if (.not. m .eq. n) then
                write (*,"(a,a,a,i0,a,i0)") 
     .               "ERROR: could not read ", trim(vvar),
     .               ". Read ",m," bytes instead of ", n
                ierr = -1
            endif

            if (adios_debug)
     .        write (*,"(a,a,a,i0)"), "read ",trim(vvar),"  bytes=",m
          enddo

          if (firstread) then
            nlen = adios_nvar*6
            n = adios_nvar*6*4
            call adios_read_var (gh, "bconds", zero, nlen, 
     .                            bconds, m)
            if (m < n) then
                write (*,"(a,i0,a,i0)") 
     .               "ERROR: could not read bconds. Read ",
     .               m," bytes instead of ",n
                ierr = -1
            endif
            do ieq=1,adios_nvar
              varray%array_var(ieq)%bconds(:) = 
     .             bconds( (ieq-1)*6+1:ieq*6 ) 
            enddo
          endif

c     End program

      end subroutine readDerivedTypeADIOS

#endif
      end module ADIOS_io

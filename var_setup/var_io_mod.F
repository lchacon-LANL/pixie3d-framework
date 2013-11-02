
c module var_io
c ######################################################################
      module var_io

        use io

        use file_io

        use variables

        integer :: ilevel

        integer :: urecord=10000,usave=20000,uinput=25
        character*(1024) :: recordsave,inputfile

        integer,pointer,dimension(:) :: murecord,musave

        logical :: clean=.false.
     $            ,tst_flg=.false.
     $            ,append_recordfile=.false.

cc        ! number of opens -> first create, then append file
cc        integer :: record_open_num=0

#if defined(adios)
        !!! ADIOS variables
        character*(1024) :: recordfile='record.bp'
        ! ADIOS needs a variable containing the MPI_WORLD_COMM communicator
        integer         :: adios_world_comm

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

        logical,private :: adios_debug=.false.
#else
        character*(1024) :: recordfile='record.bin'
#endif

      contains

#if defined(adios)
c     openRecordFileAdiosForRead
c     #################################################################
      subroutine openRecordFileAdiosForRead(ierr,file,nvar,nx,ny,nz)

        implicit none

c     Call variables

        integer,intent(out) :: ierr
        character(*),optional :: file
        integer,optional,intent(out) :: nvar,nx,ny,nz

c     Local variables

        integer :: gcnt, vcnt, acnt, nullidx, ts, lastts,i
        character (len=20), dimension(1)  :: gnamelist  ! list of groups (=1)
        character(len=1024) :: rfile

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
        if (ADIOS_VER > 1.21) then
          call adios_inq_group(adios_gh,adios_varlist,adios_attrlist
     .                        ,ts,lastts,ierr)
        else
          call adios_inq_group(adios_gh,adios_varlist,adios_attrlist
     .                        ,ierr)
        endif

        ! count /var/ variables
        adios_nvar=0
        do i=1,adios_varcount
            if (index(adios_varlist(i),"/var/") == 1) then
                adios_nvar = adios_nvar+1
cc                ! get rid of trailing char(0) in names
cc                nullidx=index(adios_varlist(i),char(0))
cc                if (nullidx == 0) nullidx=len(trim(adios_varlist(i)))+1
cc                adios_varnames(adios_nvar) = 
cc     .                    adios_varlist(i)(1:nullidx-1)
            endif
        enddo
        
        adios_tread = 0  ! next time read first timestep

        write (*,"(a,i0,a,a)")
     .       " Proc=",my_rank,": opened ",trim(rfile)

        if (adios_debug) then
          if (my_rank == 0) then
            write (*,"(a, i5)") " Number of vars: ", adios_varcount
            write (*,*), "Variables:"
            do i=1,adios_varcount
              write (*,"(i5, a, a)") i,")  ", trim(adios_varlist(i))
            enddo
          endif
        endif

        if (PRESENT(nvar)) nvar = adios_nvar

      end subroutine openRecordFileAdiosForRead

c     closeRecordFileAdiosForRead
c     #################################################################
      subroutine closeRecordFileAdiosForRead(ierr)

        implicit none

        integer, intent(out) :: ierr

        call adios_gclose(adios_gh, ierr)
        call adios_fclose(adios_fh, ierr)

        if (adios_debug) write (*,*), "Closed ", trim(recordfile)

      end subroutine closeRecordFileAdiosForRead

c     writeDerivedTypeAdios
c     #################################################################
      subroutine writeDerivedTypeAdios(handle, varray, addl_write)

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

      end subroutine writeDerivedTypeAdios

c     readDerivedTypeAdios
c     #################################################################
      subroutine readDerivedTypeAdios(gh,varray,step,firstread,ierr)

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

          call fromGlobalToLocalLimits(ilomg,jlomg,klomg,ilom,jlom,klom
     .                                ,1,1,1)
          call fromGlobalToLocalLimits(ihipg,jhipg,khipg,ihip,jhip,khip
     .                                ,1,1,1)

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

      end subroutine readDerivedTypeAdios

#endif 

c     openRestartFileForRead
c     #################################################################
      subroutine openRestartFileForRead(file,lneq)

c     -----------------------------------------------------------------
c     Opens restart file (either POSIX, or ADIOS).
c     -----------------------------------------------------------------

      implicit none

c     Global variables

      character(*),optional :: file
      integer,optional,intent(out) :: lneq

c     Local variables

      integer :: nx,ny,nz,ierr,neq
      character(len=1024) :: rfile

c     Begin program

#if !defined(adios)

      if (PRESENT(file)) then
        rfile = file
      else
        rfile = recordfile
      endif

      if (PRESENT(lneq)) then
         call probe_recordfile_nvar(trim(rfile),lneq)
cc         lneq = neq
      endif

      call openPOSIXBinFile(urecord,trim(rfile),status='old')

      read (urecord) nx,ilog,ihig
      read (urecord) ny,jlog,jhig
      read (urecord) nz,klog,khig

      if (associated(gv)) then
#if defined(petsc) 
         if (nx /= nxl .or. ny /= nyl .or. nz /= nzl) then
#else
         if (nx /= nxd .or .ny /= nyd .or. nz /= nzd) then
#endif
           call pstop('openRestartFileForRead'
     $            ,'Grid meshes do not agree; cannot restart')

         endif
      else
        nxl = nx
        nyl = ny
        nzl = nz
      endif

#else

      call adios_init_noxml (adios_err)

      ! allocate 100MB buffer for ADIOS
      call adios_allocate_buffer (100,adios_err)

      call openRecordFileAdiosForRead(adios_err,file=file,nvar=lneq)

      if (adios_err /= 0) then
        call pstop('openRestartFileForRead'
     $            ,'Error reading ADIOS restart file')
      endif

#endif

      end subroutine openRestartFileForRead

c     closeRestartFileForRead
c     #################################################################
      subroutine closeRestartFileForRead

c     -----------------------------------------------------------------
c     Close restart file (either POSIX, or ADIOS).
c     -----------------------------------------------------------------

      implicit none

c     Global variables

c     Local variables

c     Begin program

#if !defined(adios)
      close(urecord)
#else
      call closeRecordFileAdiosForRead(adios_err)

      if (adios_err /= 0) then
        call pstop('closeRestartFileForRead'
     $            ,'Error closing ADIOS restart file')
      endif
#endif

c     End program

      end subroutine closeRestartFileForRead

c     openPOSIXBinFile
c     #################################################################
      subroutine openPOSIXBinFile(unit,file,status,position,check
     .                           ,suffix)

c     -----------------------------------------------------------------
c     Opens POSIX binary file.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer,INTENT(IN) :: unit
      character(*),INTENT(IN) :: file
      character(*),optional,INTENT(IN) :: status,position
      logical,optional :: check,suffix

c     Local variables

cc      integer :: unit1
      character*(len(file)+10) :: file1
      character*(10) :: stat
      logical :: chk,sffx

c     Begin program

cc      unit1 = unit + my_rank

      if (PRESENT(status)) then
        stat=status
      else
        stat='unknown'
      endif

      if (PRESENT(check)) then
        chk = check
      else
        chk = .false.
      endif

      if (PRESENT(suffix)) then
        sffx = suffix
      else
        sffx = .true.
      endif

      if (sffx) then
        file1=trim(file)//'_proc'//trim(int2char(my_rank))
      else
        file1=trim(file)
      endif

      if (chk.and.(.not.is_file(file1))) then
        call pstop('openBinFile','Cannot find file '//trim(file1))
      endif

      if (PRESENT(position)) then
        open(unit=unit,file=trim(file1),form='unformatted'
     .       ,status=trim(stat),position=trim(position))
      else
        open(unit=unit,file=trim(file1),form='unformatted'
     .       ,status=trim(stat))
      endif

c     End

      end subroutine openPOSIXBinFile

c     initPOSIXBinFile
c     ##############################################################
      subroutine initPOSIXBinFile(unit,file)

        implicit none

c       Call variables

        integer :: unit
        character(*) :: file

c       Local variables

        integer :: ierr

        if (my_rank == 0) ierr = rm_files(file)

#if defined(petsc) 
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

        !Initialize record file
        call openPOSIXBinFile(unit,file,status='unknown')

        write (unit) nxl,ilog,ihig
        write (unit) nyl,jlog,jhig
        write (unit) nzl,klog,khig

      end subroutine initPOSIXBinFile

c     setRecordfileToAppend
c     ##############################################################
      subroutine setRecordFileToAppend

        implicit none

c       Call variables

c       Local variables

        append_recordfile = .true. !Set recordfile to append

      end subroutine setRecordfileToAppend

c     readRecordFile
c     #################################################################
      function readRecordFile(unit,itime,time,dt,varray) result(ierr)

c     -----------------------------------------------------------------
c     Reads record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ierr,itime,unit
      real(8)    :: time,dt

      type(var_array),pointer:: varray

c     Local variables

      integer    :: tmp(1)

#if defined(adios)
      !ADIOS variables
      integer*8       :: start(1), readcount(1)  ! dimensions are 64bit in adios
      integer*8       :: m ! actually read bytes
      logical, save   :: firstread = .true.
#endif

c     Begin program

      ierr = 0

#if defined(adios)
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

          call readDerivedTypeAdios(adios_gh,varray,adios_tread
     .                             ,firstread,ierr)

          if (adios_debug) then
            write (*,*), 'rank=',my_rank,' Adios read returns ierr',ierr
            write (*,*)
          endif

          adios_tread = adios_tread + 1
      else  ! EOF
         ierr = -2
      endif

      firstread = .false.

#else
      read (unit,iostat=ierr,end=100) time
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) itime
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) dt
      if (ierr /= 0) goto 200

      call readDerivedType(varray,unit,.false.,ierr)

      return

 100  ierr=-2  !EOF
cc      write (*,*) 'readRecordFile: EOF in file unit',unit
      return

 200  ierr=-1  !Error, but not EOF
      write (*,*) 'readRecordFile: Error reading file unit',unit

c     End
#endif

      end function readRecordFile

c     writeRecordFile
c     #################################################################
      subroutine writeRecordFile(file,itime,time,dt,varray,init)

c     -----------------------------------------------------------------
c     Writes record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: itime
      real(8)    :: time,dt
      type(var_array),pointer:: varray
      character(*) :: file
      logical,optional :: init

c     Local variables

      logical :: rinit

#if defined(adios)
c     ADIOS variables
      integer*8    :: handle, totalsize, groupsize
      integer      :: err
      character(2) :: mode='a'//char(0)
      logical      :: addl_write  ! true: write names and bconds

      integer      :: xsize, ysize, zsize
      integer      :: xloghost, xhighost
      integer      :: yloghost, yhighost
      integer      :: zloghost, zhighost
#else
      integer      :: unit=1111
#endif

c     Begin program

      if (PRESENT(init)) then
         rinit = init
      else
         rinit = .not.append_recordfile
      endif

#if defined(adios)

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
          append_recordfile = .true.
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
        write (*,*) 'ADIOS open: adios_file=',file
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
        call adios_write(handle,"time" ,time ,err)
        call adios_write(handle,"itime",itime,err)
        call adios_write(handle,"dt"   ,dt   ,err)
      endif

      call writeDerivedTypeAdios(handle,varray,addl_write)

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

#else

      !Initialize POSIX files
      if (rinit) then
        call initPOSIXBinFile(unit,file)
        append_recordfile = .true.
      else
        call openPOSIXBinFile(unit,file,status='old'
     .                       ,position='append'
     .                       ,check=.true.)
      endif

      write (unit) time
      write (unit) itime
      write (unit) dt

      call writeDerivedType(varray,unit,.false.)

      close(unit)
cc#if !defined(absoft_) && !defined(absoft) && !defined(lf95)
cc      call flush(unit)
cc#endif

#endif

c     End

      end subroutine writeRecordFile

c     openRecordFiles
c     #################################################################
      subroutine openRecordFiles(base_unit,base_fname,munit
     .                          ,ilg,ihg,jlg,jhg,klg,khg)

c     -----------------------------------------------------------------
c     Opens input binary files, one per processor.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: base_unit,ilg,ihg,jlg,jhg,klg,khg
      character(*),INTENT(IN) :: base_fname
      integer,pointer,dimension(:) :: munit

c     Local variables

      integer :: ierr,ifl,nf,send_buf(3),rec_buf(3)
      integer,allocatable,dimension(:) :: nx,ny,nz,ilg_,ihg_
     .                                   ,jlg_,jhg_,klg_,khg_

      logical :: check,lchck

c     Begin program

c     Open graphics files

      ierr = openFiles(my_rank,np,base_unit,base_fname,munit)

      if (ierr /= 0) then
          call pstop('openRecordFiles'
     .              ,'Cannot find record files')
      else
#if defined(petsc)
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif
        if (my_rank == 0) then
          write (*,*) 'Finished!'
          write (*,*)
        endif
      endif

c     Read init info

      nf = size(munit)

      allocate(nx(nf),ny(nf),nz(nf)
     .        ,ilg_(nf),ihg_(nf)
     .        ,jlg_(nf),jhg_(nf)
     .        ,klg_(nf),khg_(nf))

      do ifl=1,nf
        read (munit(ifl)) nx(ifl),ilg_(ifl),ihg_(ifl)
        read (munit(ifl)) ny(ifl),jlg_(ifl),jhg_(ifl)
        read (munit(ifl)) nz(ifl),klg_(ifl),khg_(ifl)
      enddo

c     Consistency check

#if defined(petsc)
      if (nf /= 1) then
          call pstop('openRecordfiles'
     .          ,'Cannot deal with multiple record files in parallel')
      endif

      ihg = ihg_(1)
      ilg = ilg_(1)
      jhg = jhg_(1)
      jlg = jlg_(1)
      khg = khg_(1)
      klg = klg_(1)

      nxl = ihg-ilg+1
      nyl = jhg-jlg+1
      nzl = khg-klg+1

      if (associated(gv%gparams)) then
        check =(     nxl /= gv%gparams%nlx
     .          .or. nyl /= gv%gparams%nly
     .          .or. nzl /= gv%gparams%nlz )

        lchck = check

        call MPI_Allreduce(lchck,check,1
     .                ,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpierr)

        if (check) then
          write (*,*) my_rank,nxl,gv%gparams%nlx
     .               ,nyl,gv%gparams%nly
     .               ,nzl,gv%gparams%nlz
          call pstop('openRecordFiles','Grid meshes do not agree')
        endif
      else
        send_buf = (/ ihg,jhg,khg /)
        call MPI_Allreduce(send_buf,rec_buf,3
     .                ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)

        nxd = rec_buf(1)
        nyd = rec_buf(2)
        nzd = rec_buf(3)
      endif
#else
      !Set global limits (NOT a parallel run)
      if (associated(gv%gparams)) then

        if (    (nxd /= maxval(ihg_)-minval(ilg_)+1)
     .      .or.(nyd /= maxval(jhg_)-minval(jlg_)+1)
     .      .or.(nzd /= maxval(khg_)-minval(klg_)+1) ) then
          call pstop('openRecordFiles','Grid size does not agree')
        endif

        ihg = nxd
        ilg = 1
        jhg = nyd
        jlg = 1
        khg = nzd
        klg = 1

      else

        nxd = maxval(ihg_)-minval(ilg_)+1
        nyd = maxval(jhg_)-minval(jlg_)+1
        nzd = maxval(khg_)-minval(klg_)+1

        ihg = nxd
        ilg = 1
        jhg = nyd
        jlg = 1
        khg = nzd
        klg = 1

      endif

#endif

      deallocate(nx,ny,nz
     .          ,ilg_,ihg_
     .          ,jlg_,jhg_
     .          ,klg_,khg_)
 
c     End program

      end subroutine openRecordFiles

c     closeRecordFiles
c     #################################################################
      subroutine closeRecordFiles

c     -----------------------------------------------------------------
c     Close record files
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      call closeFiles(murecord)

c     End programs

      end subroutine closeRecordFiles

c     writeDerivedType
c     #################################################################
      subroutine writeDerivedType(varray,unit,frmt)

        implicit none

c     Call variables

        type(var_array),pointer  :: varray
        integer :: unit
        logical :: frmt

c     Local variables

        integer :: ieq

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
            write(unit,*) varray%array_var(ieq)%bc_dep_list
            write(unit,*) varray%array_var(ieq)%dom_dep_list
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
            write(unit) varray%array_var(ieq)%bc_dep_list
            write(unit) varray%array_var(ieq)%dom_dep_list
          enddo
        endif

c     End program

      end subroutine writeDerivedType

c     readDerivedType
c     #################################################################
      subroutine readDerivedType(varray,unit,format,ierr)

        implicit none

c     Call variables

        type(var_array),pointer  :: varray
        integer :: unit,ierr
        logical :: format

c     Local variables

        integer :: ieq,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        ierr = 0

        if (format) then
#if defined(petsc)
          read(unit,*,iostat=ierr,end=100)
     .               ilomg,ihipg,jlomg,jhipg,klomg,khipg

          ihip = ihipg - ilog + 1 !ilog known from recordfile header
          jhip = jhipg - jlog + 1 !jlog known from recordfile header
          khip = khipg - klog + 1 !klog known from recordfile header
          ilom = ilomg - ilog + 1 !ilog known from recordfile header
          jlom = jlomg - jlog + 1 !jlog known from recordfile header
          klom = klomg - klog + 1 !klog known from recordfile header
#else
          read(unit,*,iostat=ierr,end=100) ilom,ihip,jlom,jhip,klom,khip
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
            read(unit,*,iostat=ierr,end=100)
     .           varray%array_var(ieq)%bc_dep_list
            if (ierr /= 0) goto 200
            read(unit,*,iostat=ierr,end=100)
     .           varray%array_var(ieq)%dom_dep_list
            if (ierr /= 0) goto 200
          enddo
        else
#if defined(petsc)
          read (unit,iostat=ierr,end=100)
     .               ilomg,ihipg,jlomg,jhipg,klomg,khipg

          ihip = ihipg - ilog + 1 !ilog known from recordfile header
          jhip = jhipg - jlog + 1 !jlog known from recordfile header
          khip = khipg - klog + 1 !klog known from recordfile header
          ilom = ilomg - ilog + 1 !ilog known from recordfile header
          jlom = jlomg - jlog + 1 !jlog known from recordfile header
          klom = klomg - klog + 1 !klog known from recordfile header
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
            read(unit,iostat=ierr,end=100)
     .           varray%array_var(ieq)%bc_dep_list
            if (ierr /= 0) goto 200
            read(unit,iostat=ierr,end=100)
     .           varray%array_var(ieq)%dom_dep_list
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

      contains

      end subroutine readDerivedType

c     probe_recordfile_nvar
c     #################################################################
      subroutine probe_recordfile_nvar(file,nvar)

        implicit none

c     Call variables

        character(*) :: file
        integer :: nvar

c     Local variables

        integer :: unit,nx,ny,nz,i1,i2,j1,j2,k1,k2,ierr,it
        real(8) :: t,ddt

c     Begin program

c     Open files

        unit = 1111

        call openPOSIXBinFile(unit,trim(file),status='old')

c     Read preamble

        read (unit,iostat=ierr,end=100) nx,i1,i2
        if (ierr /= 0) goto 200
        read (unit,iostat=ierr,end=100) ny,j1,j2
        if (ierr /= 0) goto 200
        read (unit,iostat=ierr,end=100) nz,k1,k2
        if (ierr /= 0) goto 200
        read (unit,iostat=ierr,end=100) t
        if (ierr /= 0) goto 200
        read (unit,iostat=ierr,end=100) it
        if (ierr /= 0) goto 200
        read (unit,iostat=ierr,end=100) ddt
        if (ierr /= 0) goto 200
        read (unit,iostat=ierr,end=100) i1,i2,j1,j2,k1,k2
        if (ierr /= 0) goto 200

c     Read number of variables

        read (unit) nvar

c     Close file

        close(unit)

        return

c     End program

 100    write (*,*) 'probe_recordfile_nvar: EOF in file unit ',unit
        stop

 200    write (*,*) 'probe_recordfile_nvar: Error reading file unit '
     $              ,unit
        stop

      end subroutine probe_recordfile_nvar

c     dumpVectorToFile
c     #################################################################
      subroutine dumpVectorToFile(unit,vec,frmt)

        implicit none

c     Call variables

        real(8),dimension(:,:,:,:) :: vec
        integer :: unit
        logical :: frmt

c     Local variables

c     Begin program

        if (frmt) then
          write (unit,*) ilomg,ihipg,jlomg,jhipg,klomg,khipg

          write (unit,*) vec
        else
          write (unit) ilomg,ihipg,jlomg,jhipg,klomg,khipg

          write (unit) vec
        endif

c     End program

      end subroutine dumpVectorToFile

c     readVectorFromFiles
c     #################################################################
      function readVectorFromFiles(unit,vec,format) result(ierr)

        implicit none

c     Call variables

        real(8),dimension(0:,0:,0:,:) :: vec
        integer :: unit,ierr
        logical :: format

c     Local variables

        integer :: ieq,ilom,ihip,jlom,jhip,klom,khip

c     Begin program

        ierr = 0

        if (format) then
          read(unit,*,iostat=ierr,end=100) ilom,ihip,jlom,jhip,klom,khip

          read(unit,*,iostat=ierr,end=100)
     .         vec(ilom:ihip,jlom:jhip,klom:khip,:)
          if (ierr /= 0) goto 200
        else
          read(unit,iostat=ierr,end=100) ilom,ihip,jlom,jhip,klom,khip

          read(unit,iostat=ierr,end=100)
     .         vec(ilom:ihip,jlom:jhip,klom:khip,:)
          if (ierr /= 0) goto 200
        endif

        return

 100    ierr=-2  !EOF
        write (messg,*) 'readVector: EOF in file unit',unit
        call pstop('readVector',messg)

 200    ierr=-1  !Error, but not EOF
        write (messg,*) 'Error reading file unit',unit
        call pstop('readVector',messg)

c     End program

      end function readVectorFromFiles

c     readTimeStep
c     #################################################################
      subroutine readTimeStep(itime,time,dt,varray,ierr)

c     -----------------------------------------------------------------
c     Reads time step info from single or multiple files.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray
      integer :: ierr,itime
      real(8) :: time,dt

c     Local variables

      integer :: ifile

c     Begin program

      if (associated(murecord)) then
        do ifile=1,size(murecord)
          ierr = readRecordFile(murecord(ifile),itime,time,dt,varray)
          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) cycle !EOF
        enddo
      else
        ierr = readRecordFile(urecord,itime,time,dt,varray)
      endif

c     End program

      end subroutine readTimeStep

      end module var_io

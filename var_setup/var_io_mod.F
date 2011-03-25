
c module iosetup
c ######################################################################
      module iosetup

        use io

        use file_io

        use variables

        integer :: ilevel

        integer        :: urecord=10000,usave=20000,uinput=25
        character*(50) :: recordfile='record.bin',recordsave,inputfile

        integer,pointer,dimension(:) :: murecord,musave

        logical :: debug,clean=.false.,tst_flg=.false.

#if defined(adios)
        !!! ADIOS variables
        character*(10) :: adiosfile='record.bp'
        ! ADIOS needs a variable containing the MPI_WORLD_COMM communicator
        integer        :: adios_world_comm
        ! number of opens -> first create, then append file
        integer        :: adios_open_num=0

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

        logical,private   :: adios_debug=.false.
#endif

      contains

#if defined(adios)
c     openRecordFileAdiosForRead(ierr)
c     #################################################################
      subroutine openRecordFileAdiosForRead(ierr)
        integer, intent(out) :: ierr

        integer :: gcnt, vcnt, acnt, nullidx 
        character (len=20), dimension(1)  :: gnamelist  ! list of groups (=1)

        call MPI_Comm_dup (MPI_COMM_WORLD, adios_world_comm, ierr)
        call adios_fopen (adios_fh, adiosfile, adios_world_comm, 
     .                     gcnt, ierr)
        call adios_inq_file (adios_fh,vcnt,acnt,
     .                    adios_tstart,adios_ntsteps,gnamelist,ierr)

        adios_tstep = adios_tstart + adios_ntsteps - 1
        call adios_gopen (adios_fh, adios_gh, gnamelist(1), 
     .                     adios_varcount, adios_attrcount, ierr)
        call adios_inq_group (adios_gh, adios_varlist, 
     .                       adios_attrlist, ierr)

        ! count /var/ variables
        adios_nvar=0
        do i=1,adios_varcount
            if (index(adios_varlist(i),"/var/") == 1) then
                adios_nvar = adios_nvar+1
        !        ! get rid of trailing char(0) in names
        !        nullidx=index(adios_varlist(i),char(0))
        !        if (nullidx == 0) nullidx=len(trim(adios_varlist(i)))+1
        !        adios_varnames(adios_nvar) = 
     .  !                  adios_varlist(i)(1:nullidx-1)
            endif
        enddo
        
        adios_tread = 0  ! next time read first timestep

        write (*,"(a,i0,a,a)"), " Proc=",my_rank,": opened ", adiosfile

        if (adios_debug) then
          if (my_rank == 0) then
            write (*,"(a, i5)") " Number of vars: ", adios_varcount
            write (*,*), "Variables:"
            do i=1,adios_varcount
              write (*,"(i5, a, a)") i,")  ", trim(adios_varlist(i))
            enddo
          endif
        endif

      end subroutine openRecordFileAdiosForRead

c     closeRecordFileAdiosForRead(ierr)
c     #################################################################
      subroutine closeRecordFileAdiosForRead(ierr)
        integer, intent(out) :: ierr

        call adios_gclose(adios_gh, ierr)
        call adios_fclose(adios_fh, ierr)

        if (adios_debug) write (*,*), "Closed ", adiosfile

      end subroutine closeRecordFileAdiosForRead

c     writeDerivedTypeAdios
c     #################################################################
      subroutine writeDerivedTypeAdios(handle, varray, addl_write)

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
        iloghost = 0 ! gcw
        jloghost = 0 ! gcw
        kloghost = 0 ! gcw
        ihighost = 0 ! gcw
        jhighost = 0 ! gcw
        khighost = 0 ! gcw
      
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

c     openBinFile
c     #################################################################
      subroutine openBinFile(unit,file,status,position,check)

c     -----------------------------------------------------------------
c     Opens binary file 
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer,INTENT(IN) :: unit
      character(*),INTENT(IN) :: file
      character(*),optional,INTENT(IN) :: status,position
      logical,optional :: check

c     Local variables

cc      integer :: unit1
      character*(len(file)+10) :: file1
      character*(10) :: stat
      logical :: chk

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

      file1=trim(file)//'_proc'//trim(int2char(my_rank))

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

      end subroutine openBinFile

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
     .  write (*,"(a,i0,a,i2)"), "rank=",my_rank,
     .       " adios read time=", adios_tread
      if (adios_tread < adios_ntsteps) then
          ! read 1 item from the 1D arrays, from position adios_tread
          start = adios_tread 
          readcount = 1
          if (adios_debug)
     .      write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .           ": read in time. start=",start,"count=",readcount
          call adios_read_var (adios_gh, "time", 
     .                     start, readcount, time, m)
          start = adios_tread 
          if (adios_debug)
     .      write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .           ": read in itime. start=",start,"count=",readcount
          call adios_read_var (adios_gh, "itime", 
     .                     start, readcount, itime, m)
          start = adios_tread 
          if (adios_debug)
     .      write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .           ": read in dt. start=",start,"count=",readcount

          call adios_read_var (adios_gh, "dt",  
     .                     start, readcount, dt, m)

          if (firstread) then
                start = adios_tread 
                call adios_read_var (adios_gh, "nvar", 
     .                              start, readcount, tmp, m)
                varray%nvar = tmp(1)
                ! consistency check
                if (adios_nvar.ne.varray%nvar) then
                    write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .                     " ERROR: adios_nvar=",adios_nvar,
     .                     " != varray%nvar=",varray%nvar
                endif
          endif

          if (adios_debug)
     .      write (*,"(a,i0,a,i2,a,f6.2,a,i2)"), "rank=",my_rank,
     .         ": read itime=",itime," time=",time," nvar=",varray%nvar

          call readDerivedTypeAdios(adios_gh,varray,adios_tread,
     .           firstread,ierr)

          adios_tread = adios_tread + 1
      else  ! EOF
         ierr = -2
      endif

      firstread = .false.
      if (adios_debug) write (*,*), "Adios read returns ierr",ierr

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
      subroutine writeRecordFile(unit,itime,time,dt,varray)

c     -----------------------------------------------------------------
c     Writes record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: itime,unit
      real(8)    :: time,dt
      type(var_array),pointer:: varray

c     Local variables

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
#endif

c     Begin program

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

c       calculate size actually written for each dimension
        xsize = ihip+xhighost-ilom-xloghost+1  
        ysize = jhip+yhighost-jlom-yloghost+1  
        zsize = khip+zhighost-klom-zloghost+1  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Create/Append adios file
      addl_write=.false.
      if (adios_open_num.eq.0) then
          mode = 'w'//char(0)
          if (my_rank == 0) then
c             Process 0 (arbitrary choice) 
c             in the first timestep writes 
c             the boundary conditions and variable names
              addl_write = .true.
          endif
      else 
          mode = 'a'//char(0)
      endif

      adios_open_num = adios_open_num+1
      groupsize = 10 * 4 +     ! 12 integers
     .  2 * 8 +                ! 2 doubles + nvar*double 3D arrays
     .  varray%nvar * 8 * xsize*ysize*zsize

      if (addl_write) then
         groupsize = groupsize + 4 + 4 + 4 + ! namelen
     .     varray%nvar * 6 * 4 +  ! bconds
     .     varray%nvar * 20       ! descr (var names)
      endif

c      write (*,"('rank= ',i0,' ihip=',i0,' ilom=',i0,' jhip=',
c     . i0,' jlom=',i0,' khip=',i0,' klom=',i0)") my_rank, ihip, ilom, 
c     . jhip, jlom, khip, klom

      call adios_open (handle, 'record', adiosfile, mode,
     .                 adios_world_comm, err)

c      print *, "ADIOS: adios_group_size =", groupsize
      call adios_group_size(handle, groupsize, totalsize, err)
c      print *, "ADIOS: write time", time

c     Write all variables to buffer
      if (my_rank == 0) then
          call adios_write (handle, "time", time, err)
          call adios_write (handle, "itime", itime, err)
          call adios_write (handle, "dt", dt, err)
      endif
      call writeDerivedTypeAdios(handle, varray, addl_write)

c     Close adios file which starts I/O writing from the buffer
c      print *, "ADIOS: close", time

      call adios_close (handle, err)

#else
      write (unit) time
      write (unit) itime
      write (unit) dt

      call writeDerivedType(varray,unit,.false.)
#endif

#if !defined(absoft_) && !defined(absoft) && !defined(lf95)
      call flush(unit)
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
          call pstop('openRecordfiles'
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
          call pstop('openRecordFiles','Grid meshes do not agree')
        endif
      else
        send_buf = (/ nxl,nyl,nzl /)
        call MPI_Allreduce(send_buf,rec_buf,3
     .                ,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)

        nxd = rec_buf(1)
        nyd = rec_buf(2)
        nzd = rec_buf(3)
      endif
#else
      !Set global limits (NOT a parallel run)
      if (associated(gv%gparams)) then
        ihg = nxd
        ilg = 1
        jhg = nyd
        jlg = 1
        khg = nzd
        klg = 1

        if (nf == 1) then
          if (nx(1) /= nxd .or. ny(1) /= nyd .or. nz(1) /= nzd) then
            call pstop('openRecordFiles','Grid size does not agree')
          endif
        endif
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

      end subroutine readDerivedType

c     dumpVectorToFile
c     #################################################################
      subroutine dumpVectorToFile(unit,vec,frmt)

        implicit none

c     Call variables

        real(8),dimension(:,:,:,:) :: vec
        integer          :: unit
        logical          :: frmt

c     Local variables

c     Begin program

        if (frmt) then
cc#if defined(petsc)
          write (unit,*) ilomg,ihipg,jlomg,jhipg,klomg,khipg
cc#else
cc          write (unit,*) ilom,ihip,jlom,jhip,klom,khip
cc#endif
          write (unit,*) vec
        else
cc#if defined(petsc)
          write (unit) ilomg,ihipg,jlomg,jhipg,klomg,khipg
cc          write (*,*) ilomg,ihipg,jlomg,jhipg,klomg,khipg
cc#else
cc          write (unit) ilom,ihip,jlom,jhip,klom,khip
cc#endif
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
cc     .            ,ilomg,ihipg,jlomg,jhipg,klomg,khipg

c     Begin program

        ierr = 0

        if (format) then
cc#if defined(petsc)
cc          read (unit,*,iostat=ierr,end=100)
cc     .           ilomg,ihipg,jlomg,jhipg,klomg,khipg
cc          call fromGlobalToLocalLimits(ilomg,jlomg,klomg,ilom,jlom,klom
cc     .                                ,1,1,1)
cc          call fromGlobalToLocalLimits(ihipg,jhipg,khipg,ihip,jhip,khip
cc     .                                ,1,1,1)
cc#else
          read (unit,*,iostat=ierr,end=100)ilom,ihip,jlom,jhip,klom,khip
cc#endif
          read (unit,*,iostat=ierr,end=100)
     .         vec(ilom:ihip,jlom:jhip,klom:khip,:)
          if (ierr /= 0) goto 200
        else
cc#if defined(petsc)
cc          read (unit,iostat=ierr,end=100)
cc     .         ilomg,ihipg,jlomg,jhipg,klomg,khipg
cc          call fromGlobalToLocalLimits(ilomg,jlomg,klomg,ilom,jlom,klom
cc     .                                ,1,1,1)
cc          call fromGlobalToLocalLimits(ihipg,jhipg,khipg,ihip,jhip,khip
cc     .                                ,1,1,1)
cc#else
          read (unit,iostat=ierr,end=100)ilom,ihip,jlom,jhip,klom,khip
cc#endif
          read (unit,iostat=ierr,end=100)
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

      end module iosetup

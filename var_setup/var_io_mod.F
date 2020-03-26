
c module var_io
c ######################################################################
      module var_io

        use io

        use file_io

        use variables

#if defined(adios)
        use ADIOS_io
#endif
#if defined(ADIOS2)
        use ADIOS2_mod
#endif
#if !defined(adios) && !defined(ADIOS2)
        !character*(1024) :: recordfile='record.bin'
#endif

        integer :: ilevel

        integer :: urecord=10000,usave=20000,uinput=25,u_mapf=3000
        character(1024) :: recordsave,inputfile
        character(30) :: map_file="mapfile"

        integer,pointer,dimension(:) :: murecord,musave

        logical :: clean=.false.
     $            ,tst_flg=.false.
     $            ,posix_append_recordfile=.false.

      contains

c     finalize_IO
c     #################################################################
      subroutine finalize_IO

c     -----------------------------------------------------------------
c     Finalizes IO for write
c     -----------------------------------------------------------------

      implicit none

c     Global variables

c     Local variables

      integer :: ierr
      
#if defined(adios)
      ierr = destroy_ADIOS_IO()
#endif

#if defined(ADIOS2)
      ierr = destroy_ADIOS2_IO()
#endif

c     Begin program

      end subroutine finalize_IO

c     openRestartFileForRead
c     #################################################################
      function openRestartFileForRead(file) result(ufile)

c     -----------------------------------------------------------------
c     Opens restart file (either POSIX, or ADIOS).
c     -----------------------------------------------------------------

      implicit none

c     Global variables

      integer :: ufile
      character(*),optional :: file

c     Local variables

c     Begin program

#if defined(adios)
      ufile = openADIOSRestartFileForRead(file=file)
#endif
#if defined(ADIOS2)
      ufile = openADIOS2RestartFileForRead(file="record.bp")
#endif
#if !defined(adios) && !defined(ADIOS2)
      ufile = openPOSIXRestartFileForRead(file=file)
#endif

      end function openRestartFileForRead

c     closeRestartFileForRead
c     #################################################################
      subroutine closeRestartFileForRead(ufile)

c     -----------------------------------------------------------------
c     Close restart file (either POSIX, or ADIOS).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ufile

c     Local variables

c     Begin program

#if defined(adios)
      adios_err = closeADIOSRecordFileForRead()
      if (adios_err /= 0) then
        call pstop('closeRestartFileForRead'
     $            ,'Error closing ADIOS restart file')
      endif
      adios_err = destroy_ADIOS_IO()
#endif
#if defined(ADIOS2)
      adios2_err = closeADIOS2RecordFileForRead()
      if (adios2_err /= 0) then
        call pstop('closeRestartFileForRead'
     $            ,'Error closing ADIOS restart file')
      endif
      adios2_err = destroy_ADIOS2_IO()
#endif
#if !defined(adios) && !defined(ADIOS2)
      close(ufile)
#endif

c     End program

      end subroutine closeRestartFileForRead

c     openPOSIXRestartFileForRead
c     #################################################################
      function openPOSIXRestartFileForRead(file) result(ufile)

c     -----------------------------------------------------------------
c     Opens POSIX restart file.
c     -----------------------------------------------------------------

      implicit none

c     Global variables

      integer :: ufile
      character(*),optional :: file
cc      integer,optional,intent(out) :: lneq

c     Local variables

      integer :: nx,ny,nz,ierr,neq,ndepvl,ndepsl,nauxvl,nauxsl
     .          ,ilg,ihg,jlg,jhg,klg,khg
      character(len=1024) :: rfile

c     Begin program

      if (PRESENT(file)) then
        rfile = file
      else
         rfile = recordfile
      endif

      ufile = find_unit(urecord)
      
c$$$      call openRecordFiles(ufile,trim(rfile)
c$$$     .                    ,ilg,ihg,jlg,jhg,klg,khg)
c$$$
c$$$      nx = ihg-ilg+1
c$$$      ny = jhg-jlg+1
c$$$      nz = khg-klg+1

      call openPOSIXBinFile(ufile,trim(rfile),status='old')

      read (ufile,iostat=ierr) ndepvl,ndepsl,nauxvl,nauxsl
      if (ierr == 0) then

        if (NDEPV == 0) then
           NDEPV=ndepvl
        elseif (NDEPV /= ndepvl) then
           call pstop('openPOSIXRestartFileForRead'
     $            ,'# Dep. vec. variables do not agree; cannot restart')
        endif

        if (NDEPS == 0) then
           NDEPS=ndepsl
        elseif (NDEPS /= ndepsl) then
           call pstop('openPOSIXRestartFileForRead'
     $            ,'# Dep. scl. variables do not agree; cannot restart')
        endif

        if (NAUXV == 0) then
           NAUXV=nauxvl
        elseif (NAUXV /= nauxvl) then
           call pstop('openPOSIXRestartFileForRead'
     $            ,'# Aux. vec. variables do not agree; cannot restart')
        endif

        if (NAUXS == 0) then
           NAUXS=nauxsl
        elseif (NAUXS /= nauxsl) then
           call pstop('openPOSIXRestartFileForRead'
     $            ,'# Aux. scl. variables do not agree; cannot restart')
        endif

      else  !COMPATIBILITY MODE

        rewind(ufile)

      endif

      read (ufile) nx,ilg,ihg
      read (ufile) ny,jlg,jhg
      read (ufile) nz,klg,khg

cc      if (associated(gv%gparams)) then
      if (associated(gv)) then
#if defined(petsc) 
         if (nx /= nxl .or. ny /= nyl .or. nz /= nzl) then
#else
         if (nx /= nxd .or .ny /= nyd .or. nz /= nzd) then
#endif
           write (*,*) nx,nxd,ny,nyd,nz,nzd
           call pstop('openPOSIXRestartFileForRead'
     $            ,'Grid meshes do not agree; cannot restart')

         endif
      else
        nxl = nx
        nyl = ny
        nzl = nz

#if defined(petsc)
        nxd = nx*npx
        nyd = ny*npy
        nzd = nz*npz
#else
        nxd = nx
        nyd = ny
        nzd = nz
#endif
      endif

      end function openPOSIXRestartFileForRead

c     openPOSIXBinFile
c     #################################################################
      subroutine openPOSIXBinFile(unit,file,status,position,suffix)

c     -----------------------------------------------------------------
c     Opens POSIX binary file.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer,INTENT(IN) :: unit
      character(*),INTENT(IN) :: file
      character(*),optional,INTENT(IN) :: status,position
      logical,optional :: suffix

c     Local variables

      character*(len(file)+10) :: file1
      character*(10) :: stat
      logical :: chk,sffx

c     Begin program

      if (PRESENT(status)) then
        stat=status
      else
        stat='unknown'
      endif

      chk = (trim(stat) == 'old')

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

#if defined(petsc) || defined(use_pic_mpi)
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

        !Initialize record file
        call openPOSIXBinFile(unit,file)!,status='unknown',check=.false.)

        write (unit) NDEPV,NDEPS,NAUXV,NAUXS

        write (unit) nxl,ilog,ihig
        write (unit) nyl,jlog,jhig
        write (unit) nzl,klog,khig

      end subroutine initPOSIXBinFile

c     setAppendIOMode
c     ##############################################################
      subroutine setAppendIOMode

        implicit none

c       Call variables

c       Local variables

#if defined(adios)
        call setADIOSAppendIOMode()      !Set ADIOS recordfile to append
#endif
#if defined(ADIOS2)
        call setADIOS2AppendIOMode()     !Set ADIOS2 recordfile to append
#endif
#if !defined(adios) && !defined(ADIOS2)
        posix_append_recordfile = .true. !Set POSIX recordfile to append
#endif
      end subroutine setAppendIOMode

c     readPOSIXRecordFile
c     #################################################################
      function readPOSIXRecordFile(unit,itime,time,dt,gammat,varray)
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

c     Begin program

      ierr = 0

      read (unit,iostat=ierr,end=100) time
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) itime
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) dt
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) gammat
      if (ierr /= 0) goto 200

      call readDerivedType(varray,unit,.false.,ierr)

      goto 300

 100  ierr=-2  !EOF

cc      write (*,*) 'readRecordFile: EOF in rank',my_rank
      goto 300

 200  ierr=-1  !Error, but not EOF
cc      write (*,*) 'readRecordFile: Error reading file unit',unit

 300  continue

#if defined(petsc) || defined(petsc_c)
      lerr = ierr
      call MPI_Allreduce(lerr,gerr,1
     .                ,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpierr)
      ierr = gerr(1)
#endif

c     End

      end function readPOSIXRecordFile

c     readRecordFile
c     #################################################################
      function readRecordFile(unit,itime,time,dt,gammat,varray)
     .         result(ierr)

c     -----------------------------------------------------------------
c     Reads record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ierr,itime,unit
      real(8)    :: time,dt,gammat

      type(var_array),pointer:: varray

c     Local variables

c     Begin program

      ierr = 0

#if defined(adios)
      ierr = readADIOSRecordFile(unit,itime,time,dt,gammat,varray)
#endif
#if defined(ADIOS2)
      ierr = readADIOS2RecordFile(unit,itime,time,dt,gammat,varray)
#endif
#if !defined(adios) && !defined(ADIOS2)
      ierr = readPOSIXRecordFile(unit,itime,time,dt,gammat,varray)
#endif

c     End

      end function readRecordFile

c     writePOSIXRecordFile
c     #################################################################
      subroutine writePOSIXRecordFile(file,itime,time,dt,gammat,varray
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

      integer      :: unit=1111

c     Begin program

      if (PRESENT(init)) then
         rinit = init
      else
         rinit = .not.posix_append_recordfile
      endif

      !Initialize POSIX files
      if (rinit) then
        call initPOSIXBinFile(unit,file)
        posix_append_recordfile = .true.
      else
        call openPOSIXBinFile(unit,file
     .                       ,status='old'
     .                       ,position='append')
      endif

      write (unit) time
      write (unit) itime
      write (unit) dt
      write (unit) gammat

      call writeDerivedType(varray,unit,.false.)

      close(unit)

c     End

      end subroutine writePOSIXRecordFile

c     writeRecordFile
c     #################################################################
      subroutine writeRecordFile(file,itime,time,dt,gammat,varray,init)

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

c     Begin program

#if defined(adios)
      call writeADIOSRecordFile(file,itime,time,dt,gammat,varray
     .                         ,init=init)
#endif
#if defined(ADIOS2)
      call writeADIOS2RecordFile('out.bp',itime,time,dt,gammat,varray
     .                         ,init=init)
#endif
#if !defined(adios) && !defined(ADIOS2)
      call writePOSIXRecordFile(file,itime,time,dt,gammat,varray
     .                         ,init=init)
#endif

c     End

      end subroutine writeRecordFile

c     openRecordFiles
c     #################################################################
      subroutine openRecordFiles(base_unit,base_fname
     .                          ,ilg,ihg,jlg,jhg,klg,khg)

c     -----------------------------------------------------------------
c     Opens multiple input binary files for serial postprocessing of
c     parallel run.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: base_unit,ilg,ihg,jlg,jhg,klg,khg
     .          ,ndepvl,ndepsl,nauxvl,nauxsl
      character(*),INTENT(IN) :: base_fname

c     Local variables

      integer :: ierr,ifl,nf,send_buf(3),rec_buf(3)
      integer,allocatable,dimension(:) :: nx,ny,nz,ilg_,ihg_
     .                                   ,jlg_,jhg_,klg_,khg_

      logical :: check,lchck

c     Begin program

c     Open graphics files

      ierr = openFiles(my_rank,np,base_unit,base_fname,murecord)

      if (ierr /= 0) then
          call pstop('openRecordFiles'
     .              ,'Cannot find record files')
      else
#if defined(petsc) || defined(use_pic_mpi)
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif
        if (my_rank == 0) then
          write (*,*) 'Finished!'
          write (*,*)
        endif
      endif

c     Read init info

      nf = size(murecord)

      allocate(nx(nf),ny(nf),nz(nf)
     .        ,ilg_(nf),ihg_(nf)
     .        ,jlg_(nf),jhg_(nf)
     .        ,klg_(nf),khg_(nf))

      do ifl=1,nf
        read (murecord(ifl),iostat=ierr) ndepvl,ndepsl,nauxvl,nauxsl

        if (ierr == 0) then

          if (NDEPV == 0) then
             NDEPV=ndepvl
          elseif (NDEPV /= ndepvl) then
             call pstop('openRecordFiles'
     $            ,'# Dep. vec. variables do not agree; cannot proceed')
          endif

          if (NDEPS == 0) then
             NDEPS=ndepsl
          elseif (NDEPS /= ndepsl) then
             call pstop('openRecordFiles'
     $            ,'# Dep. scl. variables do not agree; cannot proceed')
          endif

          if (NAUXV == 0) then
             NAUXV=nauxvl
          elseif (NAUXV /= nauxvl) then
             call pstop('openRecordFiles'
     $            ,'# Aux. vec. variables do not agree; cannot proceed')
          endif

          if (NAUXS == 0) then
             NAUXS=nauxsl
          elseif (NAUXS /= nauxsl) then
             call pstop('openRecordFiles'
     $            ,'# Aux. scl. variables do not agree; cannot proceed')
          endif

        else  !COMPATIBILITY MODE

          rewind(murecord(ifl))

        endif

        read (murecord(ifl)) nx(ifl),ilg_(ifl),ihg_(ifl)
        read (murecord(ifl)) ny(ifl),jlg_(ifl),jhg_(ifl)
        read (murecord(ifl)) nz(ifl),klg_(ifl),khg_(ifl)
      enddo

c     Consistency check

#if defined(petsc)
      if (np > 1) then
        call setParallelLimits
      else
        call setSerialLimits
      endif
#else
      call setSerialLimits
#endif

      deallocate(nx,ny,nz
     .          ,ilg_,ihg_
     .          ,jlg_,jhg_
     .          ,klg_,khg_)
 
c     End program

      contains

#if defined(petsc)
c     setParallelLimits
c     ################################################################
      subroutine setParallelLimits

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

cc      if (associated(gv%gparams)) then
      if (associated(gv)) then
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
          call pstop('setParallelLimits','Grid meshes do not agree')
        endif
      else
        send_buf = (/ ihg,jhg,khg /)
        call MPI_Allreduce(send_buf,rec_buf,3
     .                ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)

        nxd = rec_buf(1)
        nyd = rec_buf(2)
        nzd = rec_buf(3)
      endif

      end subroutine setParallelLimits
#endif

c     setSerialLimits
c     ################################################################
      subroutine setSerialLimits

      !Set global limits (NOT a parallel run)
cc      if (associated(gv%gparams)) then
      if (associated(gv)) then

        if (    (nxd /= maxval(ihg_)-minval(ilg_)+1)
     .      .or.(nyd /= maxval(jhg_)-minval(jlg_)+1)
     .      .or.(nzd /= maxval(khg_)-minval(klg_)+1) ) then
cc          write (*,*) ihg_
cc          write (*,*) jhg_
cc          write (*,*) khg_
cc
cc          write (*,*) maxval(ihg_)-minval(ilg_)+1
cc     .               ,maxval(jhg_)-minval(jlg_)+1
cc     .               ,maxval(khg_)-minval(klg_)+1
          call pstop('openRecordFiles','Grid sizes do not agree')
        endif

      else

        nxd = maxval(ihg_)-minval(ilg_)+1
        nyd = maxval(jhg_)-minval(jlg_)+1
        nzd = maxval(khg_)-minval(klg_)+1

      endif
      
      ihg = nxd
      ilg = 1
      jhg = nyd
      jlg = 1
      khg = nzd
      klg = 1

      end subroutine setSerialLimits

      end subroutine openRecordFiles

c     openGraphicsFiles
c     #################################################################
      subroutine openGraphicsFiles

c     -----------------------------------------------------------------
c     Open record files
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

#if defined(adios)
      urecord = openADIOSRestartFileForRead()
#endif
#if defined(ADIOS2)
      urecord = openADIOS2RestartFileForRead()
#endif
#if !defined(adios) && !defined(ADIOS2)
      if (np > 1) then
        urecord = openPOSIXRestartFileForRead()
      else
        call openRecordFiles(urecord,recordfile
     .                      ,ilog,ihig,jlog,jhig,klog,khig)
      endif
#endif
c$$$      if (np > 1) then
c$$$        urecord = openRestartFileForRead()
c$$$      else
c$$$#if defined(adios)
c$$$        urecord = openRestartFileForRead()
c$$$#else
c$$$        call openRecordFiles(urecord,recordfile
c$$$     .                      ,ilog,ihig,jlog,jhig,klog,khig)
c$$$#endif
c$$$      endif

      end subroutine openGraphicsFiles

c     closeGraphicsFiles
c     #################################################################
      subroutine closeGraphicsFiles

c     -----------------------------------------------------------------
c     Close record files
c     -----------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

#if defined(adios)
      adios_err = closeADIOSRecordFileForRead()
      if (adios_err /= 0) then
        call pstop('closeRestartFileForRead'
     $            ,'Error closing ADIOS restart file')
      endif
      adios_err = destroy_ADIOS_IO()
#endif
#if defined(ADIOS2)
      adios2_err = closeADIOS2RecordFileForRead()
      if (adios2_err /= 0) then
        call pstop('closeRestartFileForRead'
     $            ,'Error closing ADIOS restart file')
      endif
      adios2_err = destroy_ADIOS2_IO()
#endif
#if !defined(adios) && !defined(ADIOS2)
      if (associated(murecord)) then
        call closeFiles(murecord)
      else
        close(urecord)
      endif
#endif

c     End programs

      end subroutine closeGraphicsFiles

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

        !Read number of variables
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
        write (messg,*) 'EOF in file unit',unit
        call pstop('readVectorFromFiles',messg)

 200    ierr=-1  !Error, but not EOF
        write (messg,*) 'Error reading file unit',unit
        call pstop('readVectorFromFiles',messg)

c     End program

      end function readVectorFromFiles

c     readTimeStep
c     #################################################################
      subroutine readTimeStep(unit,itime,time,dt,gammat,varray,ierr)

c     -----------------------------------------------------------------
c     Reads time step info from single or multiple files.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: varray
      integer :: ierr,itime,unit
      real(8) :: time,dt,gammat

c     Local variables

      integer :: ifile

c     Begin program

      if (associated(murecord)) then
        do ifile=1,size(murecord)
          ierr = readRecordFile(murecord(ifile),itime,time,dt,gammat
     .                         ,varray)
          if (ierr == -1) cycle !Error, but not EOF
          if (ierr == -2) cycle !EOF
        enddo
      else
        ierr = readRecordFile(unit,itime,time,dt,gammat,varray)
      endif

c     End program

      end subroutine readTimeStep

c     readExtMap
c     ###############################################################
      subroutine readExtMap(g_def)

c     -----------------------------------------------------------------
c     Reads external map using POSIX files
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        real(8),pointer,dimension(:,:,:,:) :: ext_map

        integer,pointer,dimension(:) :: mumap

        integer :: nxg,nyg,nzg,ierr,ifile,igrid

c     Begin program

        igrid = 1

        nxg = g_def%nxgl(igrid)
        nyg = g_def%nygl(igrid)
        nzg = g_def%nzgl(igrid)

        if (len(trim(map_file)) == 0) then
          call pstop('readExtMap'
     .              ,'Map file not specified')
        endif

        ierr = openFiles(my_rank,np,u_mapf,map_file,mumap)

        if (ierr /= 0) then
          call pstop('readExtMap'
     .              ,'Problem reading map files')
        endif

c     Fill GLOBAL map

        allocate(ext_map(0:nxg+1,0:nyg+1,0:nzg+1,3))

        do ifile=1,size(mumap)
          ierr=readVectorFromFiles(mumap(ifile),ext_map,.false.)
        enddo

        call closeFiles(mumap)

c     Fill LOCAL map and related metrics at ALL grid levels

        call defineGridMetric(g_def
     .           ,xcar=ext_map(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
     .           ,ierr=ierr)

        deallocate(ext_map)

        if (ierr == 0) call createMGMetricHierarchy(g_def)
        
      end subroutine readExtMap

c     writeExtMap
c     ###############################################################
      subroutine writeExtMap(g_def)

c     -----------------------------------------------------------------
c     Writes external map using POSIX files
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: ierr

c     Begin program

        if (len(trim(map_file)) == 0) then
          call pstop('writeExtMap'
     .              ,'Map file not specified')
        endif

        if (my_rank == 0) ierr=rm_files(map_file)

#if defined(petsc) || defined(use_pic_mpi)
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

        call openPOSIXBinFile(u_mapf,map_file)!,status='unknown')
        call dumpVectorToFile(u_mapf,g_def%gmetric%grid(1)%car,.false.)
        close(u_mapf)

      end subroutine writeExtMap

      end module var_io

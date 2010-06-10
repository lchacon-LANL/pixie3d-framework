
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

      contains

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
cc      write (*,*) 'readRecordFile: EOF in file unit',unit
      return

 200  ierr=-1  !Error, but not EOF
      write (*,*) 'readRecordFile: Error reading file unit',unit

c     End

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

c     Begin program

      write (unit) time
      write (unit) itime
      write (unit) dt

      call writeDerivedType(varray,unit,.false.)

#if !defined(absoft_) && !defined(absoft)
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
      integer,pointer,dimension(:),INTENT(OUT) :: munit

c     Local variables

      integer :: nx,ny,nz,ierr,ifile

      logical :: check,lchck

c     Begin program

c     Open graphics files

      ierr = openFiles(my_rank,base_unit,base_fname,munit)

      if (ierr /= 0) then
          call pstop('openRecordfiles'
     .              ,'Cannot find record files')
      endif

c     Read init info

      do ifile=1,size(munit)
        read (munit(ifile)) nx,ilg,ihg
        read (munit(ifile)) ny,jlg,jhg
        read (munit(ifile)) nz,klg,khg
      enddo

c     Consistency check

#if defined(petsc)
      nxl = ihg-ilg+1
      nyl = jhg-jlg+1
      nzl = khg-klg+1

      check =(     nxl /= gv%gparams%nlx
     .        .or. nyl /= gv%gparams%nly
     .        .or. nzl /= gv%gparams%nlz )
      lchck = check

      call MPI_Allreduce(lchck,check,1
     .                ,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpierr)

      if (check) then
        call pstop('openRecordFiles','Grid meshes do not agree')
      endif
#else
      !Set global limits (NOT a parallel run)
      ihg = nxd
      ilg = 1
      jhg = nyd
      jlg = 1
      khg = nzd
      klg = 1

      if (size(munit) == 1) then
        if (nx /= nxd .or. ny /= nyd .or. nz /= nzd) then
          call pstop('openRecordFiles','Grid size does not agree')
        endif
      endif
#endif

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

cc      if (my_rank == 0) then
cc        write (*,*)
cc        write (*,*) 'Closing clean-up files...'
cc        write (*,*)
cc      endif

      call closeFiles(murecord)

c     End programs

      end subroutine closeRecordFiles

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

c     dumpVector
c     #################################################################
      subroutine dumpVector(unit,vec,frmt)

c     Call variables

        real(8),dimension(:,:,:,:) :: vec
        integer          :: unit
        logical          :: frmt

c     Local variables

c     Begin program

        if (frmt) then
#if defined(petsc)
          write (unit,*) ilomg,ihipg,jlomg,jhipg,klomg,khipg
#else
          write (unit,*) ilom,ihip,jlom,jhip,klom,khip
#endif
          write (unit,*) vec
        else
#if defined(petsc)
          write (unit) ilomg,ihipg,jlomg,jhipg,klomg,khipg
#else
          write (unit) ilom,ihip,jlom,jhip,klom,khip
#endif
          write (unit) vec
        endif

c     End program

      end subroutine dumpVector

c     readVector
c     #################################################################
      function readVector(unit,vec,format) result(ierr)

c     Call variables

        real(8),dimension(0:,0:,0:,:) :: vec
        integer :: unit,ierr
        logical :: format

c     Local variables

        integer :: ieq,ilom,ihip,jlom,jhip,klom,khip

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
          read (unit,*,iostat=ierr,end=100)
     .         vec(ilom:ihip,jlom:jhip,klom:khip,:)
          if (ierr /= 0) goto 200
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
          read (unit,iostat=ierr,end=100)
     .         vec(ilom:ihip,jlom:jhip,klom:khip,:)
          if (ierr /= 0) goto 200
        endif

        return

 100    ierr=-2  !EOF
        write (messg,*) 'readVector: EOF in file unit',unit
        call pstop('readVector',messg)
cc        return

 200    ierr=-1  !Error, but not EOF
        write (messg,*) 'Error reading file unit',unit
        call pstop('readVector',messg)
cc        return

c     End program

      end function readVector

      end module iosetup

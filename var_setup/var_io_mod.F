
c module iosetup
c ######################################################################
      module iosetup

        use io

        use variables

        integer        :: ilevel

        integer        :: urecord=10000,usave=20000,uinput=25
        character*(50) :: recordfile='record.bin',recordsave,inputfile

        integer        :: nfiles,ifile
        integer       ,allocatable,dimension(:) :: murecord,musave
        character*(50),allocatable,dimension(:) :: mrecfile,mrecsave

        logical        :: debug,clean=.false.,tst_flg=.false.

      contains

c     readRecordFile
c     #################################################################
      subroutine readRecordFile(unit,itime,time,dt,varray,ierr)

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

      end subroutine readRecordFile

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

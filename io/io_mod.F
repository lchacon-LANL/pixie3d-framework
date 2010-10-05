c module io
c ########################################################################
      module io

        integer, parameter, private :: sp = kind(1.0)
        integer, parameter, public ::
     .       dp = selected_real_kind(2*precision(1.0_sp))
        integer, parameter, public ::
     .       qp_preferred = selected_real_kind(2*precision(1.0_dp))
        integer, parameter, public ::
     .       qp = (1+sign(1,qp_preferred))/2*qp_preferred
     .           +(1-sign(1,qp_preferred))/2*dp

      contains 

c     shellcmd
c     ####################################################################
      function shellcmd(command) result (ierr)

c     --------------------------------------------------------------------
c     Runs shell command from fortran.
c     --------------------------------------------------------------------

#if defined(samrai)
      use iso_c_binding
      interface
         integer(c_int) function system(sysin) bind(c)
         use iso_c_binding
         character(kind=c_char), dimension(*) :: sysin
         end function system
         end interface
#endif
      
c     Call variables
         
         integer :: ierr
         character(*) :: command
       
c     Local variables

#if !defined(samrai) && !defined(gfortran)
        integer    :: system
        external      system
#endif

c     Begin program

        ierr = system(trim(command))

      end function shellcmd

c     is_file
c     ####################################################################
      function is_file(file)

c     --------------------------------------------------------------------
c     Checks whether "file" exists.
c     --------------------------------------------------------------------

c     Call variables

      logical :: is_file
      character(*) :: file

c     Local variables

      integer :: ierr

c     Begin program

      ierr = shellcmd('test -f '//trim(file))

      is_file = (ierr == 0)

c     End program

      end function is_file

c     rm_files
c     ####################################################################
      function rm_files(file) result(ierr)

c     --------------------------------------------------------------------
c     Removes all files starting with "file".
c     --------------------------------------------------------------------

c     Call variables

      character(*) :: file

      integer :: ierr

c     Local variables

      character(200) :: command

c     Begin program

      command = 'rm -f '//trim(file)//'* > /dev/null'
      ierr = shellcmd(trim(command))

      if (ierr /= 0) then
        write (*,*) 'IO error in rm_files=',ierr
        write (*,*) 'Command: ',trim(command)
      endif

c     End program

      end function rm_files

      end module io

c module file_io
c ########################################################################
      module file_io

      use io

      use math

      contains

c     openFiles
c     #################################################################
      function openFiles(myrank,np,base_unit,base_fname,munit)
     .         result(ierr)

c     -----------------------------------------------------------------
c     Opens files with basename "base_fname" and base unit "base_unit",
c     with suffix "_proc#", with "#" the processor number.
c     Stores file units in "munits" and file names in "mfile"
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: myrank,np,base_unit,ierr
      character*(*),INTENT(IN) :: base_fname
      integer,pointer,dimension(:) :: munit

c     Local variables

      integer :: ifile,nfiles
      logical :: check,lchck
      character(len(base_fname)+10) :: file

c     Begin program

c     Find graphics files

      if (is_file(base_fname)) then  !Serial run, single file (for backward compatibility)

        nfiles = 1

        ierr = shellcmd('mv '//trim(base_fname)//' '
     .                       //trim(base_fname)//'_proc0')

      else

        if (myrank == 0) then
          write (*,*)
          write (*,*) 'Opening '//trim(base_fname)//' files...'
        endif

        !Find out number of graphics files
cc        nfiles = -1
cc        ierr   = 0
cc        do while (ierr == 0)
cc          nfiles=nfiles+1
cc          ierr = shellcmd('test -f '//trim(base_fname)//'_proc'
cc     .                 //trim(int2char(nfiles)))
cc        enddo

        nfiles = 0
        do while (is_file(trim(base_fname)//'_proc'
     .                  //trim(int2char(nfiles))))
          nfiles=nfiles+1
        enddo

      endif

c     Check for errors

      ierr = 0

      if (nfiles == 0) then
        ierr = 1
        return
      elseif (nfiles == np) then
        nfiles = 1  !One file per processor
      endif

c     Allocate arrays and open files

      allocate(munit(nfiles))

      if (nfiles == 1) then
        munit(1) = base_unit+myrank
        file = trim(base_fname)//'_proc'//trim(int2char(myrank))
        write (*,*) '  => ',trim(file),', unit=',munit(1)

        open(munit(1),file=trim(file)
     $      ,form='unformatted',status='old')
      else
        do ifile=1,nfiles
          munit(ifile) = base_unit+ifile
          file = trim(base_fname)//'_proc'//trim(int2char(ifile-1))
          write (*,*) '  => ',trim(file),', unit=',munit(ifile)

          open(munit(ifile),file=trim(file)
     $        ,form='unformatted',status='old')
        enddo
      endif

c End programs

      end function openFiles

c     closeFiles
c     #################################################################
      subroutine closeFiles(munit)

c     -----------------------------------------------------------------
c     Close graphics files
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer,pointer,dimension(:) :: munit

c     Local variables

      integer :: ifile

c     Begin program

      if (associated(munit)) then
        do ifile=1,size(munit)
          close(munit(ifile))
        enddo
        deallocate(munit)
      endif
        
c     End programs

      end subroutine closeFiles

      end module file_io

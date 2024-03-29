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

        integer, parameter :: SUCCESS=0, NON_FATAL=1, FATAL=2

        INTERFACE is_open
          module procedure is_open_file,is_open_unit
        END INTERFACE

#if defined(Darwin)
        character(4) :: sed='gsed'
#else
        character(4) :: sed='sed'
#endif

#if defined(sed_nosymlinks)
        character(30) :: sed_options=' -i'
#else
        character(30) :: sed_options=' -i --follow-symlinks'
#endif
      contains 

c     command_files
c     ##################################################################
      subroutine command_files(narg,command)

c     -----------------------------------------------------------------
c     Gets files from command input.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer      :: narg
      character(*) :: command(narg)

c     Local variables

      integer       :: iarg

c     Begin program

      do iarg = 1,narg
        call getarg(iarg,command(iarg))
      enddo

c     End program

      end subroutine command_files

c     shellcmd
c     ####################################################################
      function shellcmd(command) result (ierr)

c     --------------------------------------------------------------------
c     Runs shell command from fortran.
c     --------------------------------------------------------------------

      implicit none
      
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

      implicit none

c     Call variables

      logical :: is_file
      character(*) :: file

c     Local variables

      integer :: ierr

c     Begin program

      INQUIRE(FILE=trim(file),exist=is_file)

c     End program

      end function is_file

#if defined(ifort)
c     is_dir
c     ####################################################################
      function is_dir(dir)

c     --------------------------------------------------------------------
c     Checks whether "dir" exists.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      logical :: is_dir
      character(*) :: dir

c     Local variables

      integer :: ierr

c     Begin program

      INQUIRE(DIRECTORY=trim(dir),exist=is_dir)

c     End program

      end function is_dir
#else
c     is_dir
c     ####################################################################
      function is_dir(dir)

c     --------------------------------------------------------------------
c     Checks whether "file" exists.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      logical :: is_dir
      character(*) :: dir

c     Local variables

      integer :: ierr

c     Begin program

      INQUIRE(FILE=trim(dir),exist=is_dir)

c     End program

      end function is_dir
#endif
      
c     is_open_file
c     ####################################################################
      function is_open_file(file)

c     --------------------------------------------------------------------
c     Checks whether "file" has been opened.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      logical :: is_open_file
      character(*) :: file

c     Local variables

      integer :: ierr

c     Begin program

      INQUIRE(FILE=trim(file),opened=is_open_file)

c     End program

      end function is_open_file

c     is_open_unit
c     ####################################################################
      function is_open_unit(ufile)

c     --------------------------------------------------------------------
c     Checks whether unit "ufile" has been opened.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      logical :: is_open_unit
      integer :: ufile

c     Local variables

c     Begin program

      INQUIRE(unit=ufile,opened=is_open_unit)

c     End program

      end function is_open_unit

c     find_unit
c     ####################################################################
      function find_unit(seed) result(ufile)

c     --------------------------------------------------------------------
c     Finds available file unit to open.
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: seed
      integer :: ufile

c     Local variables

c     Begin program

      ufile = seed
      do while(is_open_unit(ufile))
        ufile = ufile+1
      enddo

c     End program

      end function find_unit

c     rm_files
c     ####################################################################
      function rm_files(file) result(ierr)

c     --------------------------------------------------------------------
c     Removes all files starting with "file".
c     --------------------------------------------------------------------

      implicit none

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

c     add_line_after_txt
c     ####################################################################
      function add_line_after_txt(file,text,line) result(ierr)

c     --------------------------------------------------------------------
c     Adds line "line" after text "text" in ascii file "file"
c     --------------------------------------------------------------------

      implicit none

      integer :: ierr

      character(*) :: file,text,line

      character(200):: command

c     Begin program

      command=
     .       sed //' '//trim(sed_options)//' -e"/'
     .                //trim(text)//'/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'a\\' // achar(39) //' -e"'
#else
     .       // 'a\'  // achar(39) //' -e"'
#endif
     .       // trim(line)   // '" -e"}" '    !'

     .       // trim(file)

        ierr = shellcmd(trim(command))

       end function add_line_after_txt

c     replace_line
c     ####################################################################
      function replace_line(file,text,line) result(ierr)

c     --------------------------------------------------------------------
c     Replaces line containing "text" with line "line" in ascii file "file"
c     --------------------------------------------------------------------

      implicit none

      integer :: ierr

      character(*) :: file,text,line

      character(200):: command

c     Begin program

      command=sed //' '//trim(sed_options)//' -e"/'
     .                 //trim(text)//'/{" -e' // achar(39)
#if defined(pgf90) || defined(xlf)
     .       // 'c\\' // achar(39) //' -e"'
#else
     .       // 'c\'  // achar(39) //' -e"'
#endif
     .       // trim(line)   // '" -e"}" '   !'
     .       // trim(file)

      ierr = shellcmd(trim(command))

      end function replace_line

c     delete_line
c     ####################################################################
      function delete_line(file,text) result(ierr)

c     --------------------------------------------------------------------
c     Deletes line containing "text" in ascii file "file"
c     --------------------------------------------------------------------

      implicit none

      integer :: ierr

      character(*) :: file,text

      character(200):: command

c     Begin program

      command=sed //' '// trim(sed_options)//' /'//trim(text)//'/d '
     $                 // trim(file)

      ierr = shellcmd(trim(command))

      end function delete_line

c     lines_in_file
c     ####################################################################
      function lines_in_file(file) result(nlines)

c     --------------------------------------------------------------------
c     Computes number of lines in "file"
c     --------------------------------------------------------------------

      implicit none

      integer :: nlines

      character(*) :: file

      integer :: ifile,i,reason
      
c     Begin program

      ifile = find_unit(1234)
      open(unit=ifile,FILE=trim(file),action="read")

      !how many
      i = 0
      do
       read(ifile,FMT='(a)',iostat=reason)
       if (reason/=0) EXIT
       i = i+1
      end do

      nlines = i

      close(ifile)
c$$$      write(*,'(a,I0)') "Number of lines in file: " , nlines
      
      end function lines_in_file
      
c     sstop
c     ################################################################
      subroutine sstop(rank,routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        integer :: rank
        character(*)  :: routine, message

        integer :: ierr
        character(200) :: command

c     Begin program

        if (rank == 0) then
          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)
        endif

#if defined(petsc)
        call MPI_Finalize(ierr)
#endif

        call app_exit(rank,FATAL)  !Fatal termination

c     End program

      end subroutine sstop

c     push_exit_code
c     ################################################################
      subroutine push_exit_code(code)

c     ---------------------------------------------------------------
c     Pushes application exit code to file for script processing
c     ---------------------------------------------------------------

        implicit none

        integer :: code,funit

c     Begin program

        funit = find_unit(123)
        open(unit=funit,file="exit_code.txt",status="unknown")
        write (funit,10) code
 10     format(i0)
        close(funit)

c     End program

      end subroutine push_exit_code

c     app_exit
c     ################################################################
      subroutine app_exit(rank,code)

c     ---------------------------------------------------------------
c     Ends application
c     ---------------------------------------------------------------

        implicit none

        integer :: rank,code

c     Begin program

      if (rank == 0) then
        call push_exit_code(code)
      endif

      call exit(code)

c     End program

      end subroutine app_exit
      
c     int2char
c     ################################################################
      function int2char(n) result (chr)

      implicit none

      integer      :: n
      character(10):: chr

      integer      :: i,exp,k,j
      character(3) :: c

      if (abs(n) > 0) then
         exp = int(log(float(n))/log(1d1))
      else
         exp = 0
      endif

      if (n >= 0) then
        chr=''
      else
        chr='-'
      endif

      k = abs(n)
      do i=exp,0,-1
         j = k/10**i
         c = achar(48+j)
         chr = trim(chr)//trim(c)
         if (i > 0 .and. j /= 0) k = mod(k,j*10**i)
      enddo

      end function int2char

      end module io

c module file_io
c ########################################################################
      module file_io

      use io

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
      logical :: check,lchck,single_file,posix_read
      character(len(base_fname)+10) :: file

c     Begin program

c     Find graphics files

      if (is_file(base_fname)) then  !Single file (for backward compatibility)

        nfiles = 1

cc        ierr = shellcmd('ln -s '//trim(base_fname)//' '
cc     .                       //trim(base_fname)//'_proc0')
        single_file = .true.
      else

        if (myrank == 0) then
          write (*,*)
          write (*,*) 'Opening '//trim(base_fname)//' files...'
        endif

        !Find out number of graphics files
        nfiles = 0
        do while (is_file(trim(base_fname)//'_proc'
     .                  //trim(int2char(nfiles))))
          nfiles=nfiles+1
        enddo

        single_file = .false.
      endif

c     Check for errors

      ierr = 0

      if (nfiles == 0) then  !No files!
        ierr = 1
        return
      elseif (nfiles == np) then
        posix_read = .true.  !One file per processor
        nfiles = 1
      else
        posix_read = .false. !Async file reading
      endif

c     Allocate arrays and open files

      allocate(munit(nfiles))

      if (posix_read) then
        munit(1) = find_unit(base_unit+myrank)
        file = trim(base_fname)//'_proc'//trim(int2char(myrank))
        write (*,*) '  => ',trim(file),', unit=',munit(1)

        open(munit(1),file=trim(file)
     $      ,form='unformatted',status='old')
      else
        if (single_file) then
          ifile = 1
          munit(ifile) = find_unit(base_unit + myrank)
          file = base_fname
          write (*,*) '  => ',trim(file),', unit=',munit(ifile)

          open(munit(ifile),file=trim(file)
     $        ,form='unformatted',status='old')
        else
          do ifile=1,nfiles
            munit(ifile) = find_unit(base_unit+ifile)
            file = trim(base_fname)//'_proc'//trim(int2char(ifile-1))
            write (*,*) '  => ',trim(file),', unit=',munit(ifile)

            open(munit(ifile),file=trim(file)
     $          ,form='unformatted',status='old')
          enddo
        endif
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

c     getFiles
c     #################################################################
      function getFiles(bname,munit,info) result(ierr)

c     -----------------------------------------------------------------
c     Get files with basename bname
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr,info
      integer,pointer,dimension(:) :: munit
      character(*) :: bname

c     Local variables

      integer :: ifile
      real(8) :: r
      integer :: i,reason,NFiles
      character(LEN=100) :: dumstr

c     Begin program

      !get the files
      ierr = shellcmd('ls ./'//trim(bname)//'* > fileContents.txt')

      if (ierr /= 0) return
      
      NFiles = lines_in_file('fileContents.txt')

      if (NFiles == 0) then
        ierr = 1
        return
      endif

      !Open files
      ifile = find_unit(1234)
      open(ifile,file='fileContents.txt',status='old')
      allocate(munit(NFiles))
      do i = 1,NFiles
        read(ifile,'(a)') dumstr
        if (info > 0) write(*,'(a,a)') " Opening file ",trim(dumstr)
        munit(i) = find_unit(12345)
        open(unit=munit(i),file=trim(dumstr),status='old')
      end do 

      close(ifile)

      ierr = shellcmd('rm fileContents.txt')
      
c     End programs

      end function getFiles

      end module file_io

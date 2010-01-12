c module io_mod
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
#if !defined(samrai)
#if !defined(gfortran)
        integer    :: system
        external      system
#endif
#endif
c     Begin program

        ierr = system(trim(command))

      end function shellcmd

c     is_file
c     ####################################################################
      function is_file(file)

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

c     rm_file
c     ####################################################################
      subroutine rm_file(file)

c     Call variables

      character(*) :: file

c     Local variables

      integer :: ierr

c     Begin program

      ierr = shellcmd('rm -f '//trim(file)//'* > /dev/null')

c     End program

      end subroutine rm_file

      end module io

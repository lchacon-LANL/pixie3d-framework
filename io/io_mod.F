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

c     Call variables

        integer :: ierr
        character(*) :: command

c     Local variables

#if !defined(gfortran)
        integer    :: system
        external      system
#endif

c     Begin program

        ierr = system(trim(command))

      end function shellcmd

      end module io

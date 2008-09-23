c module io_mod
c ########################################################################
      module io_mod

#if defined(absoft_) || defined(absoft)
        integer    :: system
        external      system
#endif

cc      contains
cc
cc#if defined(gfortran)
cc
cc      function mysystem(command)
cc
cc      integer :: mysystem
cc      character(*) :: command
cc
cc      call system (command,mysystem)
cc
cc      end function mysystem
cc
cc#else
cc
cc#if !defined(g95)
cc
cc      function mysystem(command)
cc
cc      integer :: mysystem,system
cc      character(*) :: command
cc
cc      external system
cc
cc      mysystem = system(command)
cc
cc      end function mysystem
cc
cc#else
cc
cc      function mysystem(command)
cc
cc      integer :: mysystem,system
cc      character(*) :: command
cc
cc      mysystem = system(command)
cc
cc      end function mysystem
cc
cc#endif
cc
cc#endif
      end module io_mod

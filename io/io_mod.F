c module io_mod
c ########################################################################
      module io_mod

#if !defined(gfortran)
        integer    :: system
        external      system
#endif

      end module io_mod

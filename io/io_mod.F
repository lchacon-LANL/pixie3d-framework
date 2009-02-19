c module io_mod
c ########################################################################
      module io

#if !defined(gfortran)
        integer    :: system
        external      system
#endif

      end module io

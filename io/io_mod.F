c module io
c ########################################################################
      module io

        integer    :: system
#if !defined(g95)
        external      system
#endif

      end module io

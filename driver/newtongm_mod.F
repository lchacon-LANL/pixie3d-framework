c module newtongm
c ######################################################################
      module newtongm

        use nk_setup

#if !defined(petsc_c)
        use nk_mod
#else
        integer    :: jit=1
#endif

      end module newtongm

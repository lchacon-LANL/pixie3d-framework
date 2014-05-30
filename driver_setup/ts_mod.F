
c module timeStepping
c ######################################################################
      module timeStepping

#if defined(samrai)
        use ts_base
#else
        use prj_int
#endif

      end module timeStepping

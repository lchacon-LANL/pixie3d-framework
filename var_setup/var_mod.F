
c module variables
c #####################################################################
      module variables

        use var_bc

#if defined(samrai)
        integer :: bc_order = 1
#endif

      end module variables

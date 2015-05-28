
c module timeStepping
c ######################################################################
      module timeStepping

#if defined(samrai)
        use ts_io
#else
        use prj_int

      contains

c     TS_update
c     ######################################################################
      function TS_update(vnp,iout) result(ierr)

c     ----------------------------------------------------------------------
c     Performs TS update
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr,iout

      type(var_array),pointer :: vnp

c     Local variables
      
      INTERFACE
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use var_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c     Begin program
      
#if defined(petsc)
        !Begin profiling NL solver
        call switch_profile_stage(iNL)
#endif

c     Evaluate old-time NLF

      old_time_eval = .true.    !This informs nlfunction that we are processing
                                !  n time level info

      call evaluate_NLF(1,u_n,fold)

      old_time_eval = .false.

c     Store BDF2 old-time solution (in u_np)

      if (bdf2) call AXPYDerivedType(-cn/cnp,u_n,-cnm/cnp,u_nm,vnp)

c     Perform TS update

      if (eigen) then
        !Eigenvalue solver
        ierr = eigenSolution (vnp,iout)
      elseif (prj_int_ts) then
        ierr = prj_int_update(vnp,iout)
      else
        !Time step
        ierr = updateSolution(vnp,iout)
      endif

c     End program

      end function TS_update
#endif

      end module timeStepping

c module ts_io
c ######################################################################
      module ts_io

        use ts_base

        use graphics_io, ONLY: no_graphics_io

      contains

c     TS_chk_dump
c     ###############################################################
      function TS_chk_dump(etime,reset_counters) result(dump)

c     ---------------------------------------------------------------
c     Check for end of time loop
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      logical :: dump
      real(8),optional,INTENT(OUT) :: etime
      logical,optional :: reset_counters

c     Local variables

      real(8) :: e_time
      logical :: reset

c     Begin program

      if (PRESENT(reset_counters)) then
        reset = reset_counters
      else
        reset = .false.
      endif

#if defined(petsc)
      call PetscGetTime   (prof_time(3),prof_err)
      e_time = prof_time(3)-prof_time(1)
cc        if (my_rank == 0) write (*,*) e_time, quota
#else
      e_time = 0d0
#endif

#if defined(RFX)
      dump = (nrst.eq.ndstep.or.tmrst.ge.         dstep.or.e_time>quota)
#else
      dump = (nrst.eq.ndstep.or.tmrst.ge.0.999999*dstep.or.e_time>quota)
#endif

      if (reset.and.dump) nrst = 0

      dump = dump.or.relax_init_grid
      if (reset.and.dump.and.itime>init_grid_offset) tmrst=tmrst-dstep

      if (PRESENT(etime)) etime = e_time

c     End program

      end function TS_chk_dump

c     TS_dump
c     ##################################################################
      subroutine TS_dump(varray)

      implicit none

c     Call variables

      type(var_array),pointer :: varray

c     Local variables

      INTERFACE
        subroutine perturbEquilibrium(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

c     Begin program

        if (TS_chk_dump(etime=exec_time,reset_counters=.true.)
     .     .and.(.not.no_graphics_io)) then
cc          if (.not.relax_init_grid) nrst  = 0
cc          if(itime-init_grid_offset > 0) tmrst=tmrst-dstep
          call writeRecordFile(recordfile,itime,time,dt,gammat,varray)
        endif

        if (relax_init_grid) then
          relax_init_grid = .false.
cc          if (source) call get_src(varray,.true.) !Update source with interpolated equilibrium
          call equateDerivedType (u_n,varray)     !Save interpolated (unperturbed) equilibrium
          call perturbEquilibrium(varray)         !Perturb interpolated equilibrium
        endif

c     End program
      
      end subroutine TS_dump

c     TS_final_io
c     ##################################################################
      subroutine TS_final_io
      
      implicit none

c     Call variables

c     Local variables

      integer :: omp_threads

c     Begin program

c     Average explicit timestep

      if (cnfactor > 5d-1.and.my_rank == 0.and.(itime > inewtime-1))then
        dtexp = dtexp/(itime - inewtime+1)
        write (*,'(a,1pe11.4)') 'Average explicit time step: ',dtexp
      endif

c     PETSc CPU diagnostics

#if defined(petsc)
      call PetscGetTime   (prof_time(3),prof_err)
      call PetscGetCPUTime(prof_time(4),prof_err)

      omp_threads = 1
#ifdef _OPENMP
!$OMP PARALLEL
!$      omp_threads = omp_get_num_threads()
!$OMP END PARALLEL
#endif
      if (my_rank == 0 .and. .not.(tst_flg)) then
        write (*,*) 'CPU statistics'
        write (*,*)
        write (*,*) 'Elapsed time =', prof_time(3)-prof_time(1)
        write (*,*) 'CPU time     =',(prof_time(4)-prof_time(2))*1d4
     .                              /omp_threads
        write (*,*) 'COMM % time  ='
     .      ,(prof_time(3)-prof_time(1)
     .      -(prof_time(4)-prof_time(2))*1d4/omp_threads)
     .      /(prof_time(3)-prof_time(1))*100
      endif
#endif

c     Final statistics

      if (my_rank == 0) then
cc        prec_tot = gmres_tot + iguess*newt_tot

        write(*,300) 

        if (itime-inewtime+1 > 0 .and. newt_tot > 0) then
          write(*,310) itime,float(newt_tot )/(itime-inewtime+1)
     .                      ,float(gmres_tot)/(itime-inewtime+1)
     .                      ,float(gmres_tot)/newt_tot
        elseif (itime-inewtime+1 == 0.and.newt_tot > 0) then
          write(*,310) itime,0d0
     .                      ,0d0
     .                      ,float(gmres_tot)/newt_tot
        elseif (itime-inewtime+1 > 0.and.newt_tot == 0) then
          write(*,310) itime,float(newt_tot )/(itime-inewtime+1)
     .                      ,float(gmres_tot)/(itime-inewtime+1)
     .                      ,0d0
        else
          write(*,310) itime,0d0,0d0,0d0
        endif
      endif

c     End program

 300  format (/,'Final statistics',/,/,
     .          ' itime  Newt/ndt  GMRES/ndt  GMRES/Newt')
 310  format (i5,3x,f7.1,4x,f7.1,4x,f7.1)

      end subroutine TS_final_io
      
      end module ts_io

c module ts_io
c ######################################################################
      module ts_io

        use ts_base

        use graphics_io, ONLY: no_graphics_io,nplot,ndplot,tplot,dplot

      contains

c     TS_init_timers
c     ####################################################################
      subroutine TS_init_timers

c     --------------------------------------------------------------------
c     Initialize time stepping counters
c     --------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

c     PETSc CPU counters and profiling
      
#if defined(petsc)
      call PetscGetTime   (prof_time(1),prof_err)
      call PetscGetCPUTime(prof_time(2),prof_err)
#endif

c     End program

      end subroutine TS_init_timers

c     TS_end_timers
c     ####################################################################
      subroutine TS_end_timers

c     --------------------------------------------------------------------
c     Initialize time stepping counters
c     --------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

c     PETSc CPU counters and profiling
      
#if defined(petsc)
      call PetscGetTime   (prof_time(3),prof_err)
      call PetscGetCPUTime(prof_time(4),prof_err)
#endif

c     End program

      end subroutine TS_end_timers

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
cc      if (my_rank == 0) write (*,*) "Time counters",e_time, quota
#else
      e_time = 0d0
#endif

#if defined(RFX)
      dump = (nrst.eq.ndstep.or.tmrst.ge.dstep.or.e_time>quota)
#else
      dump = (nrst.eq.ndstep)
     .   .or.(sign(1d0,dt)*(tmrst-dstep)>=1d-10)
     .   .or.(e_time>quota)
#endif

      if (reset.and.dump) nrst = 0

      dump = dump.or.relax_init_grid

      if (reset.and.dump.and.itime>init_grid_offset) then
        tmrst=tmrst-sign(1d0,dt)*dstep
      endif

      if (PRESENT(etime)) etime = e_time

c     End program

      end function TS_chk_dump

c     TS_dump
c     ##################################################################
      subroutine TS_dump(varray,dump_now)

      implicit none

c     Call variables

      type(var_array),pointer :: varray
      logical, optional :: dump_now
      
c     Local variables

      logical :: dump
      
      INTERFACE
        subroutine perturbEquilibrium(varray)
        use var_setup
        type(var_array),pointer :: varray
        end subroutine perturbEquilibrium
      END INTERFACE

c     Begin program

      if (PRESENT(dump_now)) then
        dump = dump_now
      else
        dump = .false.
      endif
      
      if ((TS_chk_dump(etime=exec_time,reset_counters=.true.).or.dump)
     .     .and.(.not.no_graphics_io)) then
        call writeRecordFile(recordfile,itime,time,dt,gammat,varray)
      endif

#if defined(MK)
      if (relax_init_grid) then
        relax_init_grid = .false.
cc          if (source) call get_src(varray,.true.) !Update source with interpolated equilibrium
        call equateDerivedType (u_n,varray) !Save interpolated (unperturbed) equilibrium
        call perturbEquilibrium(varray) !Perturb interpolated equilibrium
      endif
#endif

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

c     CPU diagnostics

#if defined(petsc)
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

c     TS_chk_plot_dump
c     ###############################################################
      function TS_chk_plot_dump(reset) result(dump)

c     ---------------------------------------------------------------
c     Check for end of time loop
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      logical :: dump,reset

c     Local variables

c     Begin program

      dump = (nplot == ndplot.and.nplot > 0)
cc     .   .or.((time-tmplot) >= dplot-0.9*dt .and. dplot > 0d0)
     .   .or.(sign(1d0,dt)*(time-tmplot)>=dplot-.9*abs(dt).and.dplot>0.)
     .   .or.(sign(1d0,dt)*(tmax-time  )<=1d-10           .and.dplot>0.)
     .   .or.(time == 0d0)

      if (reset.and.dump) then
        if (itime.gt.0) then
          tplot = time
          nplot = itime
        endif
      endif

c     End program

      end function TS_chk_plot_dump

      end module ts_io

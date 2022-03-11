c OMP
c ######################################################################    
      module OMP

#ifdef _OPENMP
      use omp_lib
#endif

      integer :: thr_tot=-1,thr_id=-1,ilup=0
!$OMP THREADPRIVATE(thr_tot,thr_id,ilup)
      integer :: itot
!$OMP THREADPRIVATE(itot)
  
      integer,parameter :: max_thr = 300

      contains

c     set_omp_thread_id
c     ###############################################################
      subroutine set_omp_thread_id(thrnum,thrtot)

        implicit none

c     Call variables

        integer,optional :: thrnum,thrtot

c     Local variables

        logical :: io_l=.false.

c     Begin program

        if (PRESENT(thrtot)) then
          thr_tot = thrtot
        else
#ifdef _OPENMP
          thr_tot = omp_get_max_threads()
#else
          thr_tot = 1
#endif
        endif

        if (PRESENT(thrnum)) then
          thr_id = thrnum
        else
#ifdef _OPENMP
          thr_id = omp_get_thread_num()
#endif
          thr_id = max(thr_id,0)
        endif

      end subroutine set_omp_thread_id

c     get_omp_thread_id
c     ###############################################################
      subroutine get_omp_thread_id(thrnum,thrtot)

        implicit none

c     Call variables

        integer :: thrnum,thrtot

c     Local variables

c     Begin program

        if (.not.is_omp_thread_id()) call set_omp_thread_id

        thrtot = thr_tot
        thrnum = thr_id

      end subroutine get_omp_thread_id

c     is_omp_thread_id
c     ###############################################################
      function is_omp_thread_id() result(is_omp)

        implicit none

        logical :: is_omp

        is_omp = (thr_tot > 0 .and. thr_id > 0)
        
      end function is_omp_thread_id

c     check_omp_numthreads
c     ###############################################################
      subroutine check_omp_numthreads(io,thrtot)
    
        implicit none

c     Global variables

        logical :: io
        integer,INTENT(OUT),optional :: thrtot

c     Local variables

        integer :: thr_tot
        integer :: omp_get_num_threads

c     Begin program

#ifdef _OPENMP
        thr_tot = omp_get_max_threads()
c$$$!$omp parallel default(shared)
c$$$!$omp critical
c$$$        thr_tot = omp_get_num_threads()
c$$$!$omp end critical
c$$$!$omp end parallel
#else
        thr_tot = 1
#endif

        if (thr_tot == 0) then
          write (*,*) "Define OMP_NUM_THREADS environment variable"
          write (*,*) "Aborting in routine check_omp_numthreads"
          stop
        else
          if (io) then
            write (*,*)          '******************'
            write (*,'(a,i3,a)') ' Using ',thr_tot,' threads'
            write (*,*)          '******************'
          endif
        endif

        if (PRESENT(thrtot)) thrtot = thr_tot

      end subroutine check_omp_numthreads

c$$$#ifndef _OPENMP
c$$$      function omp_get_num_threads() result(i)
c$$$      integer :: i
c$$$      i = 1
c$$$      end function omp_get_num_threads
c$$$#endif

      end module OMP

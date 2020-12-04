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
          thr_tot = omp_get_num_threads()
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

        call check_omp_thread_id()

        thrtot = thr_tot
        thrnum = thr_id

      end subroutine get_omp_thread_id

c     check_omp_thread_id
c     ###############################################################
      subroutine check_omp_thread_id()

        implicit none

        if (thr_tot < 0 .or. thr_id < 0) then
          write (*,*) "Thread ID not initialized"
          write (*,*) "Aborting"
          stop
        endif
        
      end subroutine check_omp_thread_id

c     check_omp_numthreads
c     ###############################################################
      subroutine check_omp_numthreads(io)
    
        implicit none

c     Global variables

        logical :: io

c     Local variables

        integer,save :: thr_tot
!$OMP THREADPRIVATE(thr_tot)
        integer :: omp_get_num_threads

c     Begin program

#ifdef _OPENMP
!$omp parallel
        thr_tot = omp_get_num_threads()
!$omp end parallel
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

      end subroutine check_omp_numthreads

c$$$#ifndef _OPENMP
c$$$      function omp_get_num_threads() result(i)
c$$$      integer :: i
c$$$      i = 1
c$$$      end function omp_get_num_threads
c$$$#endif

      end module OMP

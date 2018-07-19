c OMP
c ######################################################################    
      module OMP

#ifdef _OPENMP
      use omp_lib
#endif

      use grid, ONLY: pstop

      integer :: thr_tot=-1,thr_num=-1
!$OMP THREADPRIVATE(thr_tot,thr_num)
  
      contains

c     set_omp_thread_id
c     ###############################################################
      subroutine set_omp_thread_id(thrnum,thrtot)

        implicit none

c     Call variables

        integer,optional,INTENT(IN) :: thrnum,thrtot

c     Local variables

        logical :: io_l=.false.

c     Begin program

        if (PRESENT(thrtot)) then
          thr_tot = thrtot
        else
          thr_tot = 1
#ifdef _OPENMP
          thr_tot = omp_get_num_threads()
#endif
        endif

        if (PRESENT(thrnum)) then
          thr_num = thrnum
        else
#ifdef _OPENMP
          thr_num = omp_get_thread_num()
#endif
          thr_num = max(thr_num,0)
        endif

      end subroutine set_omp_thread_id

c     get_omp_thread_id
c     ###############################################################
      subroutine get_omp_thread_id(thrnum,thrtot)

        implicit none

c     Call variables

        integer,INTENT(OUT) :: thrnum,thrtot

c     Local variables

c     Begin program

#ifdef _OPENMP
        thrtot = omp_get_num_threads()
        thrnum = omp_get_thread_num()
#else
        thrtot = 1
        thrnum = 0
#endif

        call set_omp_thread_id()

      end subroutine get_omp_thread_id

c     check_omp_thread_id
c     ###############################################################
      subroutine check_omp_thread_id()

        implicit none

        if (thr_tot < 0 .or. thr_num < 0) then
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

        integer :: thr_tot
cc        character*10 :: ev
        integer :: omp_get_num_threads

c     Begin program

        thr_tot = 1
        
c$$$!$      call getenv('OMP_NUM_THREADS',ev)
c$$$!$      read(ev,'(i2)') thr_tot
#ifdef _OPENMP
!$omp parallel
        thr_tot = omp_get_num_threads()
!$omp end parallel                                                             
#endif

        if (thr_tot == 0) then
          call pstop("check_omp_numthreads"
     .         ,"Define OMP_NUM_THREADS environment variable")
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

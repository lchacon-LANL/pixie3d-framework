c OMP
c ######################################################################    
      module OMP

!$    use omp_lib

      use grid, ONLY: pstop

      integer :: thr_tot=-1,thr_num=-1,ilup=0
!$OMP THREADPRIVATE(thr_tot,thr_num,ilup)
      integer :: itot
!$OMP THREADPRIVATE(itot)
  
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
          thr_tot = 1
!$        thr_tot = omp_get_num_threads()
        endif

        if (PRESENT(thrnum)) then
          thr_num = thrnum
        else
!$        thr_num = omp_get_thread_num()
          thr_num = max(thr_num,0)
        endif

      end subroutine set_omp_thread_id

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
!$omp parallel
!$      thr_tot = omp_get_num_threads()
!$omp end parallel                                                             

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

      end module OMP

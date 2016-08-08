
! module lyap
! ######################################################################
      module lyapn

        use xdraw_io

        use grid, ONLY: np,my_rank,pstop,messg

        use file_io

        use OMP

!lyap private variables
      logical :: lyap,connect
!$OMP THREADPRIVATE(lyap,connect)
       logical :: auto_corr
!$OMP THREADPRIVATE(auto_corr)
       integer :: intersect_x, intersect_y, intersect_z
!$OMP THREADPRIVATE(intersect_x,intersect_y,intersect_z)
       integer :: intersect_tor
!$OMP THREADPRIVATE(intersect_tor)
       integer :: max_inter_y,max_inter_x,max_inter_z
!$OMP THREADPRIVATE(max_inter_y,max_inter_x,max_inter_z)
       real(8) :: atc_tol,ss_min,deltass
!$OMP THREADPRIVATE(atc_tol,ss_min,deltass)
       real(8) :: cnctth
!$OMP THREADPRIVATE(cnctth)
       real(8) :: slyp, slyp_old
!$OMP THREADPRIVATE(slyp,slyp_old)
       real(8) :: big_n=100000000.d0
!$OMP THREADPRIVATE(big_n)

!#! Marco, attention
       character(200) :: init_cond_file

       integer, dimension(:), allocatable :: cnctx
!connect private variables

!#!auto_corr
!auto_corr private variables
       logical :: ex_chi
!$OMP THREADPRIVATE(ex_chi)
       real(8) :: ss_check
!$OMP THREADPRIVATE(ss_check)

!#! local variables

      contains

!#! IO routines

!     lyap_io_init
!     #################################################################
      subroutine lyap_io_init

!     -----------------------------------------------------------------
!     Sets up Lyapunov output files
!     -----------------------------------------------------------------

        implicit none

!     Call variables


!     Local variables

        integer :: ierr,ic
        integer :: unitlp
        character(50) :: drawin,dummy

        character(50) :: ssx,ssy,ssz,ssrz

!     Begin program

        ssx = 'split/lyap_'//int2char(ilup)
  
        unitlp= 56027+ilup
        open(unit=unitlp,file=trim(ssx),action='write'
     .        ,form='unformatted')

        !#! Marco attenzione
        slyp = 0.d0
        slyp_old = 0.d0

      end subroutine lyap_io_init

      subroutine lyap_write_ci(ci,ninit,str)

!     -----------------------------------------------------------------
!     Write the initial conditions on a text file
!     -----------------------------------------------------------------

       implicit none

!     Call variables
       character :: str
       integer :: ninit
       real(8) :: ci(1:ninit)

!     Local variables
       character(20) :: filename
       integer :: i

!     Begin program

       filename='init'//trim(str)//'.dat'

       open(file=trim(filename),unit=657,action='write')
        do i=1,ninit
         write(657,'(f13.8)') ci(i)
        enddo
       close(unit=657)

      end subroutine


!     lyap_io_end
!     #################################################################
      subroutine lyap_io_end

!     -----------------------------------------------------------------
!     Closes lyap output files
!     -----------------------------------------------------------------

        implicit none

!     Call variables

        integer :: nunitx

!     Begin program

          nunitx= 56027+ilup
          close(nunitx)

      end subroutine lyap_io_end

!     connect_io_init
!     #################################################################
      subroutine connect_io_init(file)

!     -----------------------------------------------------------------
!     Sets up connect IO
!     -----------------------------------------------------------------

        implicit none

!     Call variables

        character(*) :: file

!     Local variables

        integer :: cpos,nfiles,ifile,ith,iproc,ierr,ic

        character(50),dimension(:),allocatable :: cnct_file

        integer :: omp_get_num_threads

!     Begin program

        thr_tot = 1

!$omp parallel
!$      thr_tot = omp_get_num_threads()
!$omp end parallel                                                             

        cpos = index(trim(file),'.bin') - 1

        nfiles = np*thr_tot
        allocate(cnct_file(0:nfiles-1))

        do ith=0,thr_tot-1
          do iproc=0,np-1
            ifile = ith+thr_tot*iproc
            cnct_file(ifile) = 'connect_proc'//int2char(ifile)
          enddo
        enddo

!!# Marco: *.bin files opening
        allocate(cnctx(0:thr_tot-1))
        do ic = 0, thr_tot - 1
         cnctx(ic) = 91027+ic
        enddo

        do ic = 0, thr_tot - 1
         ifile = ic + thr_tot*my_rank 
         open(unit=cnctx(ic),file=trim(cnct_file(ifile))
     .        ,form='unformatted')
	enddo

        deallocate(cnct_file)

      end subroutine connect_io_init

!     connect_io_end
!     #################################################################
      subroutine connect_io_end

!     -----------------------------------------------------------------
!     Closes XDRAW surface of section output files
!     -----------------------------------------------------------------

        implicit none

!     Call variables

!     Local variables

        integer :: ith

!     Begin program

        do ith=0,thr_tot-1
          close(cnctx(ith))
        enddo

        deallocate(cnctx)

      end subroutine connect_io_end

c     auto_corr_io_init
c     ###############################################################
      subroutine auto_corr_io_init(deltass, ss_min, x00,y00,z00
     .                             ,ss_check)

!     -----------------------------------------------------------------
!     Initialize some files needed for auto_correlation computation
!     -----------------------------------------------------------------

        implicit none

!     Call variables
        real(8) :: deltass, ss_min
        real(8), intent(out) :: ss_check
        real(8) :: x00, y00, z00

!     Local variables
	character(100) :: outfile1,outfile2
        character(100) :: outfile3,outfile4,outfile5
        character(100) :: str_fin_len

        write(str_fin_len,'(i0)') ilup	

!#! REMIND: depending on which is the selected 'time-like' variable
!choose wheter to open outfile{1,3,5} (respectively {L,th,z})
        outfile1='split/auto_cors'//trim(str_fin_len)
        outfile2='split/auto_corr'//trim(str_fin_len)
        outfile3='split/auto_corth'//trim(str_fin_len)
        outfile4='split/deltass'
        outfile5='split/auto_corz'//trim(str_fin_len)
!#! to be uncommented if time is arclength
!         open(unit=173+ilup,file=outfile1
!     .        ,form='unformatted'
!     .        ,action='write',status='unknown')
!#! to be uncommented if time is theta
!         open(unit=200173+ilup,file=outfile3
!     .        ,form='unformatted'
!     .        ,action='write',status='unknown')
!#! to be uncommented if time is zeta
         open(unit=400173+ilup,file=outfile5
     .        ,form='unformatted'
     .        ,action='write',status='unknown')

!#! never comment!
         open(unit=100173+ilup,file=outfile2
     .        ,form='unformatted'
     .        ,action='write',status='unknown')
!#! never comment!
         inquire(file = outfile4, exist = ex_chi)
         if (.not.ex_chi) then
         open(unit=17,file=outfile4
     .        ,form='unformatted'
     .        ,action='write',status='unknown')

          write(*,*) 'inside', deltass
          write(17) deltass
          write(17) ss_min
         close(17)
         else
          continue
         endif
!#! REMIND: depending on which is the selected 'time-like' variable
!choose wheter to change sscheck
!#! to be uncommented if time is arclength
!#! this is valid when the time is the orbit arclength
!         ss_check = ss_min
!#! to be uncommented if time is theta
!#! this is valid when the time is the poloidal angle
!         ss_check = y00
!#! to be uncommented if time is zeta
!#! this is valid when the time is the z
         ss_check = z00
         write(100173+ilup) x00,y00,z00
         write(*,*) 'ilup',ilup,x00,y00,z00

      end subroutine auto_corr_io_init

c     auto_corr_read_ext
c     ###############################################################
      subroutine auto_corr_read_ext(init_cond_file,temp4,itot)

!     -----------------------------------------------------------------
!     Read initial conditions from external file
!     -----------------------------------------------------------------

        implicit none

!     Call variables
        character*200 :: init_cond_file
        real(8), intent(out) :: temp4(100000,0:2)
        integer, intent(out) :: itot

!     Local variables
        character*200 :: tfile
        real(8) :: tot, temp0
        real(8) :: temp3(0:2)
        integer :: iost, part

!#! initial conditions from file, to have initial conditions on the same
!#! flux surface
       part=1
       tfile = trim(init_cond_file)
       tot=0.
       write(*,*) '^^^^^^^auto_corr^^^^^^^'
       open(unit=657,file=trim(tfile),
     .      form='unformatted',status='unknown')
       read(657,iostat=iost) temp0
       write(*,*) 'psirange= ',temp0
       read(657,iostat=iost) temp0
       write(*,*) 'max chi= ',temp0
       read(657,iostat=iost) temp0
       write(*,*) 'level= ',temp0
       write(*,*) '^^^^^^^auto_corr^^^^^^^'

       do
        read(657,iostat=iost) temp3(0:2)
        tot = tot + 1.
        temp4(part,:)=temp3
        part = part + 1
        if (iost .ne. 0) then
         exit
        endif
       enddo
       close(657)

       itot=int(tot)
!       allocate(x0(itot))
!       allocate(y0(itot))
!       allocate(z0(itot))
!       x0(1:itot)=temp4(1:itot,0)
!       y0(1:itot)=temp4(1:itot,1)
!       z0(1:itot)=temp4(1:itot,2)

      end subroutine auto_corr_read_ext

c     auto_corr_io_end
c     ###############################################################
      subroutine auto_corr_io_end(lssx,lssy,lssz,rzss,x0,y0,z0,ss)

!     -----------------------------------------------------------------
!     Closes auto_corr output files
!     -----------------------------------------------------------------

        implicit none

!     Call variables

        real(8) :: x0,y0,z0,ss
        logical :: lssx,lssy,lssz,rzss

!       Write final length of the orbit in a vector
!#! REMIND: depending on which is the selected 'time-like' variable
!choose wheter to close {173,200173,400173} (respectively {L,th,z})
!#! to be uncommented if time is arclength
!          close(173+ilup)
!          fin_len(nfin_len) = ss
!          nfin_len = nfin_len + 1
!#! never comment
        close(100173+ilup)
!#! to be uncommented if time is theta
!          close(200173+ilup)
!#! to be uncommented if time is zeta
        close(400173+ilup)

        if (lssx) then
          write(*,'(a,3f12.4,3i6,f12.4)') ' # of intersect with the X,y,z
     .    surf',x0,y0,z0,intersect_x,intersect_y,intersect_z,ss
        endif
      
        if (lssy) then
          write(*,'(a,3f12.4,3i6,f12.4)') ' # of intersect with the x,Y,z
     .    surf',x0,y0,z0,intersect_x,intersect_y,intersect_z,ss
        endif
      
        if (lssz) then
          write(*,'(a,3f12.4,3i6,f12.4)') ' # of intersect with the x,y,Z
     .    surf',x0,y0,z0,intersect_x,intersect_y,intersect_z,ss
        endif
      
        if (rzss) then
          write(*,'(a,3f12.4,4i6,f12.4)') ' # of intersect with the x,y,z
     .    and TOR surf',x0,y0,z0,intersect_x,intersect_y,intersect_z,
     .       intersect_tor,ss
        endif

      end subroutine auto_corr_io_end

!#! Lyapunov routines

c     lyap_write
c     ###############################################################
      subroutine lyap_write(x,ss)

!     ------------------------------------------------------------------------
!     Write orbit radial position when needed
!#! attention: manually change resolution
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables

        real(8) , intent(in) :: x,ss

!     Local variables

        if (mod(ss,0.003d0) .lt. 1.e-3 ) then
          if (ss-slyp_old .gt. 0.0015) then 
            write(56027+ilup) x,ss
            slyp = ss
          endif
        endif
        
        slyp_old = slyp

      end subroutine lyap_write

c     connect_write
c     ###############################################################
      subroutine connect_write(x_0,y_0,z_0,ss_0,thr_num)

!     ------------------------------------------------------------------------
!     Write orbit position when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
        real(8) :: x_0, y_0, z_0, ss_0
        integer :: thr_num

!     Local variables

        write(*,'(a,i2,i2,4f15.8)') 'a1, thr_num ',my_rank, thr_num
     .           ,x_0,y_0,z_0,ss_0
        write(cnctx(thr_num)) x_0,y_0,z_0,ss_0

      end subroutine connect_write

c     auto_corr_write
c     ###############################################################
      subroutine auto_corr_write(x,y,z,atc_tol,deltass,ss_chk)

!     ------------------------------------------------------------------------
!     Write orbit position when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
        real(8) :: x, y, z     
        real(8) :: deltass, atc_tol
        real(8), intent(inout) :: ss_chk

!     Local variables

!#! old, wrong, version with time=arc-lenght
!#! remember to uncomment the proper definition of ss_chk
!#! remember to open the proper output files
!!            write(*,*) 'ss',ss,mod(ss,deltass)
!	   if (ss .gt. ss_min .and. abs(mod(ss,deltass)) .lt. atc_tol) then
!	    !check if step in ss is greater than deltass
!            if ((ss-ss_chk) .gt. 3./2.*deltass) then
!	     write(*,*) 'Reduce time step or raise atc_tol'
!	     exit
!            endif
!            if ((ss - ss_chk) .gt. deltass/2.) then
!!             write(*,*) 'condition',ss,ss-ss_chk
!	     write(173+ilup) ss
!	     write(100173+ilup) x,y,z
!	     ss_chk = ss
!            endif
!	   endif
!	  endif

!#! to be uncommented if time is theta
!#! new version with time=theta
!#! remember to uncomment the proper definition of ss_chk
!#! remember to open the proper output files
!!            write(*,'(a,2f11.6,i7)') 'th',y,abs(mod(y,deltass)),itime
!!	   if (y .gt. y0 .and. abs(mod(y,deltass)) .lt. atc_tol) then
!	   if (abs(mod(y,deltass)) .lt. atc_tol) then
!
!!#! part of code if theta is growing
!	    !check if step in ss is greater than deltass
!            if (abs((y-ss_chk)) .gt. 3./2.*deltass) then
!	     write(*,*) 'Reduce time step or raise atc_tol'
!	     exit
!            endif
!            if (abs((y - ss_chk)) .gt. deltass/2.) then
!!             if ((y - ss_chk .gt. .0)) then
!!!              write(*,*) 'forward',y,y-ss_chk
!!             else
!!!              write(*,*) 'backward',y,y-ss_chk
!!             endif
!	     write(200173+ilup) y
!	     write(100173+ilup) x,y,z
!	     ss_chk = y
!            endif
!	   endif
!	  endif

!#! to be uncommented if time is zeta
!#! new version with time=zeta
!#! remember to uncomment the proper definition of ss_chk
!#! remember to open the proper output files
!	   if (z .gt. z0 .and. abs(mod(z,deltass)) .lt. atc_tol) then
	   if (abs(mod(z,deltass)) .lt. atc_tol) then
	    !check if step in ss is greater than deltass
            if (abs((z-ss_chk)) .gt. 3./2.*deltass) then
	     write(*,*) 'Reduce time step or raise atc_tol'
	     !exit
             return
            endif
            if (abs(z - ss_chk) .gt. deltass/2.) then
	     write(400173+ilup) z
	     write(100173+ilup) x,y,z
	     ss_chk = z
            endif
	   endif
      end subroutine auto_corr_write

c     auto_corr_counters_init
c     ###############################################################
      subroutine auto_corr_counters_init()
      
!     ------------------------------------------------------------------------
!      Init autocorrelation counters
!     ------------------------------------------------------------------------

      implicit none

      intersect_x = 0
      intersect_y = 0
      intersect_z = 0
      intersect_tor = 0

      end subroutine auto_corr_counters_init

c     auto_corr_intersx_update
c     ###############################################################
      subroutine auto_corr_intersx_update(intersect_x)

!     ------------------------------------------------------------------------
!      Stop the orbit integration when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
      integer :: intersect_x

!     Local variables
 
      intersect_x = intersect_x + 1

      end subroutine auto_corr_intersx_update

c     auto_corr_intersy_update
c     ###############################################################
      subroutine auto_corr_intersy_update(intersect_y)

!     ------------------------------------------------------------------------
!      Stop the orbit integration when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
      integer :: intersect_y

!     Local variables
 
      intersect_y = intersect_y + 1

      end subroutine auto_corr_intersy_update

c     auto_corr_intersz_update
c     ###############################################################
      subroutine auto_corr_intersz_update(intersect_z)

!     ------------------------------------------------------------------------
!      Stop the orbit integration when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
      integer :: intersect_z

!     Local variables
 
      intersect_z = intersect_z + 1

      end subroutine auto_corr_intersz_update

c     auto_corr_interstor_update
c     ###############################################################
      subroutine auto_corr_interstor_update(intersect_tor)

!     ------------------------------------------------------------------------
!      Stop the orbit integration when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
      integer :: intersect_tor

!     Local variables
 
      intersect_tor = intersect_tor + 1

      end subroutine auto_corr_interstor_update

c     auto_corr_inters_stop
c     ###############################################################
      function auto_corr_inters_stop(lssx,lssy,lssz) result(orbit_exit)

!     ------------------------------------------------------------------------
!      Stop the orbit integration when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables

        logical :: lssx,lssy,lssz,orbit_exit

!     Local variables

        orbit_exit = (lssx .and. (intersect_x .ge. max_inter_x))
     $           .or.(lssy .and. (intersect_y .ge. max_inter_y))
     $           .or.(lssz .and. (intersect_z .ge. max_inter_z))
 
      end function auto_corr_inters_stop

      end module lyapn

c     set_lyapn_counter
c     ###############################################################
      subroutine set_lyapn_counter(ilop)

       	use lyapn

        implicit none

	integer :: ilop

	ilup = ilop

      end subroutine set_lyapn_counter

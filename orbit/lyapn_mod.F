
! module lyap
! ######################################################################
      module lyapn

        use xdraw_io

        use grid, ONLY: np,my_rank,pstop,messg

        use file_io

        use OMP

        logical :: lssx,lssy,lssz,rzss,ss_io=.false.

!private variables
       logical :: lyap=.false.,connect=.false.,ridge=.false.
!$OMP THREADPRIVATE(lyap,connect,ridge)
       logical :: auto_corr=.false., ac_rtrn
!$OMP THREADPRIVATE(auto_corr,ac_rtrn)
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
       integer :: delta_y,delta_z
!$OMP THREADPRIVATE(delta_y,delta_z)
       integer :: ilup = 0

!#! Marco, attention
       character(200) :: init_cond_file

       integer, dimension(:), allocatable :: cnctx
       integer, dimension(:), allocatable :: ridgex

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
!     read_ext_file
!     #################################################################
      subroutine read_ext_file(init_cond_file,nextci,temp4,itot)

!     -----------------------------------------------------------------
!     Read initial conditions from external file
!     -----------------------------------------------------------------

        implicit none

!     Call variables
        character*200 :: init_cond_file
        integer(8), intent(in) :: nextci
        real(8), intent(out) :: temp4(0:nextci,0:2)
        integer, intent(out) :: itot

!     Local variables
        character*200 :: tfile
        real(8) :: tot, temp0
        real(8) :: temp3(0:2)
        integer :: iost,part

!#! initial conditions from file, to have initial conditions on the same
!#! flux surface
       tfile = trim(init_cond_file)
       tot=0.
       part = 1
       open(unit=657,file=trim(tfile),
     .      form='unformatted',status='unknown')

       do
        read(657,iostat=iost) temp3(0:2)
!        write(*,*) 'DIAG_read',temp3
        tot = tot + 1.
        temp4(part,:)=temp3
        part = part + 1
        if (iost .ne. 0) then
         exit
        endif
       enddo
       close(657)

       itot=int(tot)

      end subroutine read_ext_file

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

!     lyap_write_ci
!     #################################################################
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

!     ridge_io_init
!     #################################################################
      subroutine ridge_io_init(file,y_surf,z_surf,rz_ss)

!     -----------------------------------------------------------------
!     Sets up ridge IO
!     -----------------------------------------------------------------

        implicit none

!     Call variables

        character(*) :: file
        logical :: y_surf, z_surf
        logical, optional :: rz_ss

!     Local variables

        integer :: cpos,nfiles,ifile,ith,iproc,ierr,ic

        character(50),dimension(:),allocatable :: ridge_file

!     Begin program
        thr_tot = 1

#ifdef _OPENMP
!$omp parallel
        thr_tot = omp_get_num_threads()
!$omp end parallel                                                             
#endif

        if (PRESENT(rz_ss)) then
          rzss = rz_ss
        else
          rzss = .false.
        endif

        ss_io = .true.   !Turn on SS XDRAW IO

        cpos = index(trim(file),'.bin') - 1

        nfiles = np*thr_tot
        write(*,*) "DIAG_ridge_io", np, thr_tot

        allocate(ridge_file(0:nfiles-1))

        if (y_surf) then
         do ith=0,thr_tot-1
           do iproc=0,np-1
             ifile = ith+thr_tot*iproc
             ridge_file(ifile) = 'ridgey_proc'//int2char(ifile)
           enddo
         enddo

!!# Marco: *.bin files opening
         allocate(ridgex(0:thr_tot-1))
         do ic = 0, thr_tot - 1
          ridgex(ic) = 51127+ic
         enddo
 
         do ic = 0, thr_tot - 1
          ifile = ic + thr_tot*my_rank 
          open(unit=ridgex(ic),file=trim(ridge_file(ifile))
     .        ,form='unformatted')
         enddo
        endif

        if (z_surf) then
         do ith=0,thr_tot-1
           do iproc=0,np-1
             ifile = ith+thr_tot*iproc
             ridge_file(ifile) = 'ridgez_proc'//int2char(ifile)
           enddo
         enddo

!!# Marco: *.bin files opening
         allocate(ridgex(0:thr_tot-1))
         do ic = 0, thr_tot - 1
          ridgex(ic) = 51127+ic
         enddo
!         write(*,*) 'DIAGa',ridgex 
         do ic = 0, thr_tot - 1
          ifile = ic + thr_tot*my_rank 
          open(unit=ridgex(ic),file=trim(ridge_file(ifile))
     .        ,form='unformatted')
         enddo
        endif

        deallocate(ridge_file)

      end subroutine ridge_io_init

!     ridge_io_end
!     #################################################################
      subroutine ridge_io_end

!     -----------------------------------------------------------------
!     Closes XDRAW surface of section output files
!     -----------------------------------------------------------------

        implicit none

!     Call variables

!     Local variables

        integer :: ith

!     Begin program

        do ith=0,thr_tot-1
          close(ridgex(ith))
        enddo

        deallocate(ridgex)

      end subroutine ridge_io_end

!     ridge_write
!     #################################################################
      subroutine ridge_write(x_1,x_2,thrnum)

!     ------------------------------------------------------------------------
!     Write intersections when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
        real(8) :: x_1, x_2
        integer :: thrnum

!     Local variables

!        write(*,'(a,i4,i4,i7,2f18.8)') 'DIAGw1, thr_num ',my_rank
!     .           ,thrnum
!     .           ,ridgex(thrnum),x_1,x_2
        write(ridgex(thrnum)) x_1,x_2

      end subroutine ridge_write

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

!#! Marco, post-processing routines to determine the Lyapunov exponent.

!     orbit_couples
!     #################################################################
      subroutine lyap_orbit_couples(idxx,idyy,idzz)

       implicit none
       
       integer, intent(in) :: idxx,idyy,idzz
       integer :: idxaa, idxbb, k
       character(100) :: tfile

       if (idxx .gt. 2 .or. idxx .lt. 0) then 
        write(*,*) 'Lyapunov exponent can be computed only between a
     .              couple of orbits! ciao!',idxx
        stop
       endif
       if (idxx .eq. 2) then

!#! Marco, write the coupled orbits, in case a new run of lypnv2 is
!needed (mainly to tune the parameters alfa, delta, nmax etc)
        tfile='couples'//trim(int2char(my_rank))//'.bin'
        open(unit=my_rank+16674,file=trim(tfile),
     .     form='unformatted',status='unknown',
     .     action='write')

        do k=0, (idxx*idyy*idzz) / 2 - 1
         write(my_rank+16674) my_rank*idxx*idyy*idzz + 2 * k 
     .            , my_rank*idxx*idyy*idzz + 2 * k + 1
         write(*,*) 'cp  ', my_rank*idxx*idyy*idzz + 2 * k 
     .            , my_rank*idxx*idyy*idzz + 2 * k + 1
         idxaa = my_rank*idxx*idyy*idzz + 2 * k
         idxbb = my_rank*idxx*idyy*idzz + 2 * k + 1
         call lypnv2(idxaa,idxbb)
        enddo

        close(my_rank+16674)

       else !(idxx eq 1)
        do k=0, (idxx*idyy*idzz) - 1
!#! to be determined
         write(*,*) k, k + idxx*idyy*idzz , 'attenzione, write
     .             code, though this should be right'
         idxaa = k
         idxbb = k + idxx*idyy*idzz
         call lypnv2(idxaa,idxbb)
         exit
        enddo 
       endif

!#! final computation, must be done serially, at the end of the parallel
!part
      
 
      end subroutine lyap_orbit_couples

!     lypnv2
!     #################################################################
      subroutine lypnv2(idxa,idxb)

       implicit none

       integer, intent(in) :: idxa, idxb

!#! local variables 
       character(100) :: tfile, string, orb0, orb1, kk
       integer :: ninx, maxq
       integer :: io, i, j, k, p, q, unitt, a
 
       real(8) :: dumr, dums, deltad
       real(8) :: dumx(0:1)
!#! dum matrix: 
!#!             dimension1: internal index of the couple (0,1)
!#!             dimension2: radial coordinate, arc-length
!#!             dimension3: steps
       real(8),allocatable :: dum(:,:,:)
       real(8),allocatable :: couples(:), doubling(:)
       real(8) :: mut_dist, catch_len, arc_len
       integer :: qmax(0:1), n_int

!#! parameters for Lyapunov calculation
!#! deltan=alfa^n * delta
       real(8) :: delta, alfa, tol
       integer :: nmax
       real(8), allocatable :: deltan(:)

       alfa = 1.10d0
       delta = 5.e-5
       nmax = 97

       allocate(dum(0:1,0:1,0:1000000))
       allocate(deltan(0:nmax))
       allocate(couples(0:nmax))
       allocate(doubling(0:nmax))
      
       do i=0,nmax
        deltan(i) = alfa**(i) * delta
       enddo
 
!!#! read the files
       q=0
       io=0
       open(unit=my_rank+235,file=
     .      'split/lyap_'//trim(int2char(idxa)),
     .      form='unformatted',action='read',
     .      status='unknown')

       do while (io .eq. 0)
        read(my_rank+235,IOSTAT=io) dumx
        dum(0,0,q) = dumx(0)
        dum(0,1,q) = dumx(1)
       
        q = q + 1 
        if (q .gt. size(dum,3)) then
         write(*,*) 'lyap: enlarge temporary array'
         exit
        endif
       
        if (io .lt. 0) then 
         qmax(0) = q
!#! test to see if the particles arrive at the fixed max-arclenght in
!the same number of integration steps
!          write(*,*) i,p,qmax
         continue
        endif

       enddo
       close(my_rank+235)

       q=0
       p=1
       open(unit=my_rank+236,file=
     .      'split/lyap_'//trim(int2char(idxb)),
     .      form='unformatted',action='read',
     .      status='unknown')

       io = 0
       do while (io .eq. 0)
        read(my_rank+236,IOSTAT=io) dumx
        dum(1,0,q) = dumx(0)
        dum(1,1,q) = dumx(1)
       
        q = q + 1 
       
        if (io .lt. 0) then 
         qmax(1) = q
!#! test to see if the particles arrive at the fixed max-arclenght in
!the same number of integration steps
!          write(*,*) i,p,qmax(i,p)
         continue
        endif

       enddo
       close(my_rank+236)

!#! define the tollerance used to look for the positions where the
!distance between the two particles forming a couple are equal to
!deltan(i)
       tol = 1.e-6
       couples = 0.d0
   
       if (my_rank .eq. 0) then
        tfile='split/par_lypnv'
        open(unit=26574,file=trim(tfile),
     .      form='unformatted',status='unknown',
     .      action='write')
         write(26574) alfa
         write(26574) delta
         write(26574) nmax
         write(26574) tol
        close(26574)
       endif

       write(*,*) 'lyapunov computation'
!#! calculate mutual distance
       maxq = max(qmax(0),qmax(1))
       catch_len = deltan(1)
       n_int = 1

       lp2: do q=0, maxq-1
         mut_dist = abs(dum(1,0,q) - dum(0,0,q))
         arc_len = (dum(0,1,q) + dum(1,1,q)) / 2.d0
         deltad = abs(dum(0,1,q)-dum(1,1,q)) / arc_len * 100.d0
         if (deltad .gt. 50.d0) then
          write(*,*) 'attention reduce time step',deltad,n_int
          exit
         endif
!#! in couples(i,n_int) it is written the arclength needed to arrive to
!deltan(n_int)
!#! once the mutual distance deltan(n_int) is reached, look for the next
!one, deltan(n_int+1)
         do while ( abs(mut_dist) .gt. catch_len ) 
          couples(n_int) = arc_len
!          write(*,'(a,i4,2f13.7)') '*+*',n_int,catch_len,couples(n_int)
          n_int = n_int + 1
          catch_len = deltan(n_int)
          if (n_int .gt. int(nmax)) then
!           write(*,*) 'out from the cycle',n_int,int(nmax)
           exit lp2
          endif
         enddo
       end do lp2
        
       write(*,*) 'N_sep reached',n_int


!!#! in the new version of couples the element couples(i,n_int)
!!represents the arclength needed to pass from deltan(n_int-1) to
!!deltan(n_int)
       doubling(0) = 0.d0
       do q=1,nmax
        doubling(q) = couples(q)-couples(q-1)
       enddo

       tfile='split/coupl'//trim(int2char(idxa))
       open(unit=my_rank+6574,file=tfile,
     .     form='unformatted',status='unknown',
     .     action='write')
       write(unit=my_rank+6574) nmax
       write(unit=my_rank+6574) n_int
       do q=0,nmax-1
        write(unit=my_rank+6574) couples(q),doubling(q)
!        write(*,'(i3,2f13.6)') q,couples(q),doubling(q)
       enddo
       close(unit=my_rank+6574)
!       write(kk,'(i0)') idxa 
!       tfile='split/blyap'//trim(kk)
!       open(unit=6575,file=tfile,
!     .     form='unformatted',status='unknown',
!     .     action='write')
!!       write(unit=6575) nmax
!!       write(unit=6575) n_int
!       do q=0,nmax
!        write(unit=6575) couples(q)
!       enddo
!       close(unit=6575)
       deallocate(deltan,couples,doubling,dum)

      end subroutine lypnv2

!     lypnv3
!     #################################################################
      subroutine lypnv3(ncoupl,npp)
  
       implicit none

       integer, intent(in) :: ncoupl, npp

       character(100) :: tfile
       integer :: i,maxn,n_int0,n_int,iunit,q,k,l,j
       integer :: i0, i1, nmax
       real(8) :: delta, alfa, tol
       real(8) :: dum0, dum1
       integer, allocatable :: couples(:)
       real(8), allocatable :: deltan(:)
       real(8), allocatable :: arc(:,:,:)
       real(8), allocatable :: lyap(:), coupl(:)



       tfile='split/par_lypnv'
       open(unit=26574,file=trim(tfile),
     .     form='unformatted',status='unknown',
     .     action='read')
        read(26574) alfa
        read(26574) delta
        read(26574) nmax
        read(26574) tol
       close(26574)
!       write(*,*) alfa,delta,nmax,tol

       allocate(deltan(0:nmax))
       do i=0,nmax
        deltan(i) = alfa**(i) * delta
       enddo
    
       allocate(couples(npp*(ncoupl+1)))
       do k=0, npp-1
        tfile='couples'//trim(int2char(k))//'.bin'
        open(unit=k+16574,file=trim(tfile),
     .       form='unformatted',status='unknown',
     .       action='read')
         do l=0,ncoupl
          read(k+16574) i0,i1
          couples(k*(ncoupl+1)+l) = i0
         end do
        close(k+16574)
       end do

!#! Marco: compute the maximum number of separation distances reached by
!all the couples of orbits
       n_int = 1
       allocate(arc(0:1,0:size(couples,1),0:nmax))
       do q=0, size(couples,1) - 1

        tfile='split/coupl'//trim(int2char(2*q))
!        write(*,*) 'tfile',tfile
        iunit = 26575+q
        open(unit=iunit,file=trim(tfile),
     .     form='unformatted',status='unknown',
     .     action='read')
!#! Marco, read the split/coupl* files. Their names are coupl(2*k) with
!k=[0:ncoupl]
        read(unit=iunit) maxn
        read(unit=iunit) n_int0
        if (n_int0 .gt. n_int) then
         n_int = n_int0
        endif

!#! Marco. arc(0,*,*) contains 'couples'
!#! Marco. arc(1,*,*) contains 'doubling'
        do k=0,nmax-1
         read(unit=iunit) dum0, dum1
         arc(0,q,k) = dum0
         arc(1,q,k) = dum1
        enddo

        close(iunit)
       enddo

       allocate(lyap(0:n_int),coupl(0:n_int))
      
       do j=0, n_int 
        lyap(j) = sum(arc(1,:,j)) / (max(1,size(arc(1,:,j))))
        coupl(j) = sum(arc(0,:,j)) / (max(1,size(arc(0,:,j))))
       enddo

       open(1,file='delta.dat',form='unformatted')
         write(1) n_int
         write(1) deltan(0:n_int)
       close(1)
       write(*,*) '****************'
       open(2,file='lyap.dat',form='unformatted')
         write(2) lyap(0:n_int)
       close(2)
       write(*,*) '****************'
       open(3,file='coupl.dat',form='unformatted')
         write(3) coupl(0:n_int)
       close(3)
       write(*,*) '****************'

       deallocate(deltan, coupl, arc, lyap)
      end subroutine lypnv3


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

!     Begin program

        thr_tot = 1

#ifdef _OPENMP
!$omp parallel
        thr_tot = omp_get_num_threads()
!$omp end parallel                                                             
#endif

        cpos = index(trim(file),'.bin') - 1

        nfiles = np*thr_tot
        allocate(cnct_file(0:nfiles-1))

        do ith=0,thr_tot-1
          do iproc=0,np-1
            ifile = ith+thr_tot*iproc
            cnct_file(ifile) =
     .           'connect_proc'//trim(int2char(ifile))//'.bin'
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
        outfile1='split/auto_cors'//int2char(ilup)
        outfile2='split/auto_corr'//int2char(ilup)
        outfile3='split/auto_corth'//int2char(ilup)
        outfile4='split/deltass'
        outfile5='split/auto_corz'//int2char(ilup)
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

        write(*,*)
        write(*,'(a,i2,i2,4e15.7)') 'a1, thr_num ',my_rank, thr_num
     .           ,x_0,y_0,z_0,ss_0
        write(cnctx(thr_num)) x_0,y_0,z_0,ss_0

      end subroutine connect_write

c     auto_corr_write
c     ###############################################################
      subroutine auto_corr_write(x,y,z,atc_tol,deltass,ss_chk,rtrn)

!     ------------------------------------------------------------------------
!     Write orbit position when needed
!     ------------------------------------------------------------------------
       
        implicit none

!     Call variables
        real(8) :: x, y, z     
        real(8) :: deltass, atc_tol
        real(8), intent(inout) :: ss_chk
        logical, intent(out) :: rtrn

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
             rtrn = .true.
!             return
            endif
!            if (abs(z - ss_chk) .gt. deltass/2.) then
            if (abs(z - ss_chk) .gt. deltass/2. 
     .          .and. (.not. rtrn)) then
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

!#! routine per ridges
c     write_final_cond
c     ###############################################################
      subroutine write_final_cond(zsurf,numberic,NstationFiles
     .    ,ddeltar,ddeltath)

      implicit none 

      ! Local variables
      real(8) :: r
      real(8), intent(in) :: zsurf,ddeltar,ddeltath ! input
      real(8), allocatable,dimension(:,:) :: real_points_start
      real(8), allocatable,dimension(:,:) :: real_points_end
      integer :: number_ci_for_thread!,thrtt
      integer, parameter :: n_intersection = 1
      integer :: n_data_tot,dumm,counter,qq,tenmod
      real(8), dimension(2) :: dati
      character(50) l_ll
      real(8), allocatable,dimension(:,:) :: punti
      integer :: i,ll,ii,aa,reason,iStation
      integer, intent(in) :: NstationFiles
      character(LEN=100), dimension(:), allocatable :: stationFileNames
      character(50) :: filename 
      integer, intent(in):: numberic !input 

!      call system('ls ridgez* > fileContents.txt')
!      open(3331,FILE='fileContents.txt',action="read")
!
!      !how many
!      i = 0
!
!      do
!       read(3331,FMT='(a)',iostat=reason) r
!       write(*,*) r
!       if (reason/=0) EXIT
!       i = i+1
!      end do
!      close(3331)
!
!      NstationFiles = i
      number_ci_for_thread = int(numberic/NstationFiles)
      write(*,*) 'NstationFiles',NstationFiles,number_ci_for_thread
!     .            ,thrtt

!      open(unit=290,file='data.dat',status='old',
!     .       form='formatted',action='write')
      open(unit=190,file='final_conditions.bin',status='unknown',
     .       form='unformatted',action='write')


!------------- inizio a leggere e scrivere------------

      allocate (punti(0:1,0:number_ci_for_thread
     .           *(1+n_intersection) - 1))


      allocate (real_points_start(0:1,
     .          0:(number_ci_for_thread
     .           *(1+n_intersection) - 1)/10))

      allocate (real_points_end(0:1,
     .          0:(number_ci_for_thread
     .           *(1+n_intersection) - 1)/10))
!#! Marco, set to zero the vectors
      punti(:,:) = 0.d0
      real_points_start(:,:) = 0.d0
      real_points_end(:,:) = 0.d0
      counter = 0

      do ll = 0,NstationFiles-1
     
       filename = 'ridgez_proc'//int2char(ll)
       write(*,*) 'DIAGfilename',filename,ll

       open(unit=180,file=filename,status='old',
     .       form='unformatted',action='read')

!        write(*,*) 'DIAG2',number_ci_for_thread,n_intersection
!        do ii=0,number_ci_for_thread*(1+n_intersection) - 1
!         read(180,iostat=dumm), dati
!!         write(*,*) 'dati  ', ii,dati     
!         if (dumm .eq. 0) then
!          punti(0:1,ii) = dati
!         else
!          exit
!         endif
!        enddo

        do ii=0,number_ci_for_thread*(1+n_intersection) - 1
         read(180,iostat=dumm), dati
!         write(*,*) 'dati  ', ii,dati     
         if (dumm .eq. 0) then
          punti(0:1,ii) = dati
         else 
          tenmod = mod(ii,10)
!#! Marco, nel file ridgez_proc* devono esserci number_ci_for_thread*2
!valori scritti
          counter = counter + (number_ci_for_thread*2- (ii-tenmod))
          write(*,*) 'error reading, skipping at i=',ii
          write(*,*) 'number of ci not written '
     .       , number_ci_for_thread*2 - (ii-tenmod)
          !#! remove broken couple c.i.+advanced
           punti(0,ii-tenmod:ii) = 0.5d0 
           punti(1,ii-tenmod:ii) = 0.5d0
           do qq=ii-10-tenmod,ii 
            write(*,*) 'check',qq,punti(0,qq)
           enddo
          exit
         endif
        enddo



!        do ii=0,number_ci_for_thread*(1+n_intersection) - 1
!         read(180,iostat=dumm), dati
!!         write(*,*) 'dati  ', ii,dati     
!         if (dumm .eq. 0) then
!          punti(0:1,ii) = dati
!         else 
!          counter = counter + (number_ci_for_thread*2- ii)
!          write(*,*) 'error reading, skipping at i=',ii,counter
!          exit
!         endif
!        enddo

        n_data_tot = number_ci_for_thread*(1+n_intersection)
        ! end of ll-th thread data reading
        !get real data

!        deltar = abs(punti(0,1) - punti(0,0))
!        deltath = abs(punti(1,3) - punti(1,0))
        real_points_start = punti(0:1,0:n_data_tot:10)
        real_points_end = punti(0:1,1:n_data_tot:10)

        do aa = 0,size(real_points_end(0,:)) -1
               
         write(190), real_points_end(0,aa),
     .               real_points_end(1,aa), zsurf

         write(190), real_points_end(0,aa) + ddeltar,
     .               real_points_end(1,aa), zsurf

         write(190), real_points_end(0,aa) - ddeltar,
     .               real_points_end(1,aa), zsurf

         write(190), real_points_end(0,aa),
     .               real_points_end(1,aa) + ddeltath, 
     .               zsurf

         write(190), real_points_end(0,aa),
     .               real_points_end(1,aa) - ddeltath,
     .               zsurf


!         write(290,'(2f13.7)'),
!     .    real_points_end(0,aa),real_points_end(1,aa)
!
!         write(290,'(2f13.7)'), real_points_end(0,aa) + deltar,
!     .               real_points_end(1,aa)
!
!         write(290,'(2f13.7)'), real_points_end(0,aa) - deltar,
!     .               real_points_end(1,aa)
!
!         write(290,'(2f13.7)'), real_points_end(0,aa),
!     .               real_points_end(1,aa) + deltath
!
!         write(290,'(2f13.7)'), real_points_end(0,aa),
!     .               real_points_end(1,aa) - deltath
!

        enddo
        close(180) ! close reading file
       enddo ! close reading loop
       close(190) !close writing file
       write(*,*) 'Ho mancato ',counter,' condizioni iniziali su '
     .        ,numberic

      end subroutine write_final_cond

      end module lyapn

c     set_lyapn_counter
c     ###############################################################
      subroutine set_lyapn_counter(ilop)

       	use lyapn

        implicit none

	integer :: ilop

	ilup = ilop

      end subroutine set_lyapn_counter

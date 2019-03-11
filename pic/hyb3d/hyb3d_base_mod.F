#include "pic_defs.h"

c pic_io
c ######################################################################    
      module pic_io

      use math

      use io, ONLY: find_unit,is_file
      
#if defined(petsc)
      use grid, ONLY: MPI_COMM_WORLD,MPI_DOUBLE_PRECISION,MPI_SUM
     .               ,MPI_INTEGER,mpierr
#endif

      integer, parameter :: nsp_max=5

      type :: particle
        SEQUENCE
        real(_prec) :: x_np(_Npg,3)
        real(_prec) :: v_np(_Npg,3)
        real(_prec) :: x_n(_Npg,3)
        real(_prec) :: v_n(_Npg,3)
        integer :: ijk_np(_Npg), ijk_n(_Npg)
#if _Npg == 4
        integer :: inpdum(_Npg),indum(_Npg) !Padding in double prec
#endif
        real(_prec) :: w_np(_Npg),w_n(_Npg)
      end type particle

      type :: species
        real(8) :: q,m,q_m
        type(particle),allocatable,dimension(:) :: pcles
      end type species

      type(species),allocatable,dimension(:) :: spcs

      integer,parameter :: Imass=1,Imom=2,Ien=3,Icmx=4,Icmy=5,Icmz=6
      real(8) :: m_np(_Dim+5),m_n(_Dim+5),PoynFlux,cm_ex(3)

      integer :: nip=0                   !Number of inactive particles

      integer :: gnip,pgroup
 
      integer :: nproc=1,my_rank_pic=0

      contains

c     alloc_pcles
c     ##################################################################
      subroutine alloc_pcles(n_sp,n_p)

      implicit none

c     Call variables

      integer :: n_sp
      integer(8) :: n_p(n_sp)

c     Local variables

      integer :: is,np_s,np_l,np_r

c     Begin program

c     Allocate pcle quantities

      allocate(spcs(n_sp))

      do is=1,n_sp
        np_s = n_p(is)
        if(    mod(np_s,pgroup).ne.0) then
          print *, "adjust Np (", np_s ,") of species ", is,
     $           " to be a multiple of pgroup (", pgroup ,")"
          stop
        end if

        !distribute particles to mpi-processes
        np_l = np_s/(nproc*pgroup)
        np_r = mod(np_s/pgroup,nproc)
        if(my_rank_pic.lt.np_r) then
          np_l = np_l+1
        end if
        allocate(spcs(is)%pcles((np_l*pgroup/_Npg)))
      enddo

      !Change pcle buffer pointer to local storage
      if(mod(gnip,pgroup).ne.0 ) then
        print *, "adjust nip (", gnip ,") to be a multiple of pgroup ("
     $        , pgroup ,")"
        stop
      end if

      nip = gnip/(nproc*pgroup)
      np_r= mod(gnip/pgroup,nproc)

      if(my_rank_pic.lt.np_r) then
        nip = (nip + 1)*pgroup  !Group pcles
      else
        nip = nip*pgroup        !Group pcles
      end if

      nip = nip/_Npg

c     End program

      end subroutine alloc_pcles

c     dealloc_pcles
c     ##################################################################
      subroutine dealloc_pcles

      implicit none

c     Call variables

c     Local variables

      integer :: is,n_sp

c     Begin program

      n_sp = size(spcs)

      do is=1,n_sp
        deallocate(spcs(is)%pcles)
      enddo

      deallocate(spcs)

c     End program

      end subroutine dealloc_pcles

c     pic_read_pcles
c     ######################################################################
      function pic_read_pcles(time,no_time_chk,file) result(ierror)

c     ----------------------------------------------------------------------
c     Reads particle info for restart
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierror
      real(8) :: time
      character(*),optional :: file

      logical,optional :: no_time_chk

c     Local variables

      integer :: isp,ip,np,ufile,nsp
      integer(8) :: npcl(nsp_max)
      character(30) :: pcle_rstrt_file
      real(8) :: pcle_time
      logical :: tchk

c     Begin program

      if (PRESENT(no_time_chk)) then
        tchk = .not.no_time_chk
      else
        tchk = .true.
      endif

      ierror = 0

      if (PRESENT(file)) then
        pcle_rstrt_file=trim(file)//trim(int2char(my_rank_pic))//'.bin'
      else
        pcle_rstrt_file='pcle_rstrt_proc'//trim(int2char(my_rank_pic))
     .                //'.bin'
      endif
      
      if (.not.is_file(pcle_rstrt_file)) then
        call p_stop("pic_read_pcles","PCLE restart file does not exist")
      endif

      ufile = find_unit(1234)

      open(ufile,file=trim(pcle_rstrt_file),form='unformatted'
     .    ,status='unknown')

      if (my_rank_pic == 0)
     .     write (*,*) ' Reading PCLES restart file(s)...'

      read(ufile,iostat=ierror,end=100) pcle_time
      if (ierror /= 0) goto 200

      if (tchk) then
        if (pcle_time /= time)
     .       call p_stop("perturbEq","Wrong restart time")
      else
        time = pcle_time
      endif

      read(ufile,iostat=ierror,end=100) nip,nsp,npcl,pgroup,gnip
      if (ierror /= 0) goto 200

      read(ufile,iostat=ierror,end=100) m_np,m_n,cm_ex,PoynFlux
      if (ierror /= 0) goto 200

      if (.not.allocated(spcs)) call alloc_pcles(nsp,npcl(1:nsp))

      do isp=1,nsp
         read (ufile,iostat=ierror,end=100) spcs(isp)%q
         if (ierror /= 0) goto 200
         read (ufile,iostat=ierror,end=100) spcs(isp)%m
         if (ierror /= 0) goto 200
         read (ufile,iostat=ierror,end=100) spcs(isp)%q_m
         if (ierror /= 0) goto 200

         np = npcl(isp)
         do ip=1,np
           read (ufile,iostat=ierror,end=100) spcs(isp)%pcles(ip)
           if (ierror /= 0) goto 200
         end do
      enddo

      close(ufile)

      if (my_rank_pic == 0) write (*,*) ' Done!'

      return

 100  ierror=-2                   !EOF
      write (*,*) 'pic_read_pcles: EOF in file unit',ufile
      return

 200  ierror=-1                   !Error, but not EOF
      write (*,*) 'pic_read_pcles: Error reading file unit',ufile

c     End program

      end function pic_read_pcles

c     pic_dump_pcles
c     ######################################################################
      function pic_dump_pcles(time,file) result(ierror)

c     ----------------------------------------------------------------------
c     Reads particle info for restart
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierror
      real(8) :: time
      character(*),optional :: file

c     Local variables

      integer :: isp,ip,np,ufile,nsp
      integer(8) :: npcl(nsp_max)
      character(30) :: pcle_rstrt_file

c     Begin program

      ierror = 0

      if (PRESENT(file)) then
        pcle_rstrt_file=trim(file)//trim(int2char(my_rank_pic))
     .       //'.bin'
      else
        pcle_rstrt_file='pcle_rstrt_proc'//trim(int2char(my_rank_pic))
     .       //'.bin'
      endif

      ufile = find_unit(1234)

      open(ufile,file=trim(pcle_rstrt_file),form='unformatted'
     .    ,status='unknown')

      nsp = size(spcs)
      do isp=1,nsp
        npcl(isp) = size(spcs(isp)%pcles)
      enddo

      write (ufile) time
      write (ufile) nip,nsp,npcl,pgroup,gnip
      write (ufile) m_np,m_n,cm_ex,PoynFlux

      do isp=1,nsp
         write (ufile) spcs(isp)%q
         write (ufile) spcs(isp)%m
         write (ufile) spcs(isp)%q_m

         np = npcl(isp)
         do ip=1,np
           write (ufile) spcs(isp)%pcles(ip)
         end do
      enddo

      close(ufile)

c     End program

      end function pic_dump_pcles

c     p_stop
c     ################################################################
      subroutine p_stop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        character(*)  :: routine, message

c     Begin program

        if (my_rank_pic == 0) then
          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)
        endif
#if defined(petsc)
        call PetscEnd(mpierr)
#endif

      end subroutine p_stop

      end module pic_io

c pic_BCs
c ######################################################################    
      module pic_BCs

      use grid, ONLY: bcond,PER,NEU,DIR,DEF,SP,OPN

      use pic_io, ONLY: p_stop

      !B.C.
      integer, parameter :: pPER=PER,pSP=SP,pOPN=OPN,pREF=-2

      integer :: pbcond(6)
      real(8) :: vinject(3,2*_Dimx)=0d0  !flow speed at boundaries for open b.c.

      contains

c     defineBCs_pcles
c     ###############################################################
      subroutine defineBCs_pcles
c     ---------------------------------------------------------------
c     Defines boundary conditions of particles
c     ---------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

c     Begin program

      pbcond = bcond  !Topological BCs

      where (pbcond == DEF) pbcond = pREF

c     End program

      end subroutine defineBCs_pcles

      end module pic_BCs

c pic_shape_f
c ######################################################################    
      module pic_shape_f

      real(8),private :: ihx=0d0,ihy=0d0,ihz=0d0
      real(8),private :: ihxsq,ihysq,ihzsq

      contains

c     xform_pcle_idx_scl
c     ##################################################################
      subroutine xform_pcle_idx_scl(ijk,nx,ny,i,j,k) 

      implicit none

      integer :: ijk,nx,ny,i,j,k

!      j = (ijk-1)/nx + 1
!      i = ijk - nx*(j-1)
!cc      i  = mod(ijk, nx) + 1

!ADAM 3D:
      k = (ijk-1)/(nx*ny) + 1
      j = ((ijk-nx*ny*(k-1))-1)/nx + 1
      i = ijk - nx*ny*(k-1) - nx*(j-1)

      end subroutine xform_pcle_idx_scl

c     xform_pcle_idx_sort
c     ##################################################################
      subroutine xform_pcle_idx_sort(ijk,nx,i,j) ! Adam. Need to modify this too - single particle.

      implicit none

      integer :: ijk,nx,ny,i,j,k

      j = (ijk-1)/nx + 1
      i = ijk - nx*(j-1)
cc      i  = mod(ijk, nx) + 1

      end subroutine xform_pcle_idx_sort

c     xform_pcle_idx
c     ##################################################################
!      subroutine xform_pcle_idx(ijk,nx,i,j)
      subroutine xform_pcle_idx(ijk,nx,ny,i,j,k)

      implicit none

!      integer :: ijk(:),nx,i(:),j(:)
      integer :: ijk(:),nx,ny,i(:),j(:),k(:)

!      j = (ijk-1)/nx + 1
!      i = ijk - nx*(j-1)

!Adam 3D:
      k = (ijk-1)/(nx*ny) + 1
      j = ((ijk-nx*ny*(k-1))-1)/nx + 1
      i = ijk - nx*ny*(k-1) - nx*(j-1)

cc      i  = mod(ijk, nx) + 1

      end subroutine xform_pcle_idx

c     setup_pic_shape_f
c     ##################################################################
      subroutine setup_pic_shape_f(ih_x,ih_y,ih_z)

      implicit none

      real(8) :: ih_x,ih_y,ih_z

      ihx = ih_x
      ihy = ih_y
      ihz = ih_z

      ihxsq = ihx*ihx!1d0/(hx*hx)
      ihysq = ihy*ihy!1d0/(hy*hy)
      ihzsq = ihz*ihz!1d0/(hz*hz)

      end subroutine setup_pic_shape_f

c     first order shape function (scalar)
c     ##################################################################
      function shape1_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x
      real(_prec) :: ss(2)

c     Local variables

c     Begin program

      ss(1) = 1d0 - x
      ss(2) = x

c     End program

      end function shape1_offset_scalar

c     2nd order shape function
c     ##################################################################
      function shape2_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x
      real(_prec) :: ss(3)

c     Local variables

c     Begin program

      ss(1) = 0.5d0*(1d0-x)*(1d0-x)
      ss(2) = 0.75d0-(x-0.5d0)*(x-0.5d0)
      ss(3) = 0.5d0*x*x    

c     End program

      end function shape2_offset_scalar

c     Dx_shape2_offset
c     ##################################################################
      function Dx_shape2_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x,ss(3)

c     Local variables

c     Begin program

      ss(1) = (x-1d0)*ihx
      ss(2) = (-2d0*x+1d0)*ihx
      ss(3) = x*ihx

c     End program

      end function Dx_shape2_offset_scalar

c     Dy_shape2_offset
c     ##################################################################
      function Dy_shape2_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x,ss(3)

c     Local variables

c     Begin program

      ss(1) = (x-1d0)*ihy
      ss(2) = (-2d0*x+1d0)*ihy
      ss(3) = x*ihy

c     End program

      end function Dy_shape2_offset_scalar

c     Dz_shape2_offset
c     ##################################################################
      function Dz_shape2_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x,ss(3)

c     Local variables

c     Begin program

      ss(1) = (x-1d0)*ihz
      ss(2) = (-2d0*x+1d0)*ihz
      ss(3) = x*ihz

c     End program                                                                                                                                           
      end function Dz_shape2_offset_scalar

c     Dxx_shape2_offset
c     ##################################################################
      function Dxx_shape2_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x,ss(3)

c     Local variables

c     Begin program

      ss(1) = ihxsq
      ss(2) =-2d0*ihxsq
      ss(3) = ihxsq

c     End program

      end function Dxx_shape2_offset_scalar

c     Dyy_shape2_offset
c     ##################################################################
      function Dyy_shape2_offset_scalar(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x,ss(3)

c     Local variables

c     Begin program

      ss(1) = ihysq
      ss(2) =-2d0*ihysq
      ss(3) = ihysq

c     End program

      end function Dyy_shape2_offset_scalar

c     first order shape function
c     ##################################################################
      subroutine shape1_offset(x,ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,2)

c     Local variables

c     Begin program

      ss(:,1) = 1d0 - x
      ss(:,2) = x

c     End program

      end subroutine shape1_offset

c     2nd order shape function
c     ##################################################################
      subroutine shape2_offset(x,ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,3)

c     Local variables

c     Begin program

      ss(:,1) = 0.5d0*(1d0-x)*(1d0-x)
      ss(:,2) = 0.75d0-(x-0.5d0)*(x-0.5d0)
      ss(:,3) = 0.5d0*x*x    

c     End program

      end subroutine shape2_offset

c     Dx_shape2_offset
c     ##################################################################
      function Dx_shape2_offset(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,3)

c     Local variables

c     Begin program

      ss(:,1) = (x-1d0)*ihx
      ss(:,2) = (-2d0*x+1d0)*ihx
      ss(:,3) = x*ihx

c     End program

      end function Dx_shape2_offset

c     Dy_shape2_offset
c     ##################################################################
      function Dy_shape2_offset(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,3)

c     Local variables

c     Begin program

      ss(:,1) = (x-1d0)*ihy
      ss(:,2) = (-2d0*x+1d0)*ihy
      ss(:,3) = x*ihy

c     End program

      end function Dy_shape2_offset

c     Dz_shape2_offset
c     ##################################################################
      function Dz_shape2_offset(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,3)

c     Local variables

c     Begin program

      ss(:,1) = (x-1d0)*ihz
      ss(:,2) = (-2d0*x+1d0)*ihz
      ss(:,3) = x*ihz

c     End program

      end function Dz_shape2_offset

c     Dxx_shape2_offset
c     ##################################################################
      function Dxx_shape2_offset(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,3)

c     Local variables

c     Begin program

      ss(:,1) = ihxsq
      ss(:,2) =-2d0*ihxsq
      ss(:,3) = ihxsq

c     End program

      end function Dxx_shape2_offset

c     Dyy_shape2_offset
c     ##################################################################
      function Dyy_shape2_offset(x) result(ss)

      implicit none

c     Call variables

      real(_prec) :: x(_Npg),ss(_Npg,3)

c     Local variables

c     Begin program

      ss(:,1) = ihysq
      ss(:,2) =-2d0*ihysq
      ss(:,3) = ihysq

c     End program

      end function Dyy_shape2_offset

! ADAM: MAY NEED Dzz for curvilinear computation of gamma - however gamma not currently used.

c     Optimized matrix dot product (2) (C. Mikida)
c     ##################################################################
cc      function dot_product_2_pcle(A,B) result(C)
      subroutine dot_product_2_pcle(A,B,C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,2), B(_Npg,2) 
      real(_prec) :: C(_Npg)

c     Local variables
      integer :: i,j

      C(:) = 0.0d0

      do i=1, 2
         !VECTOR SIMD
          do j=1,_Npg
            C(j) = C(j) + A(j,i)*B(j,i)
          end do
      end do

      end subroutine dot_product_2_pcle
cc      end function dot_product_2_pcle

c     Optimized matrix dot product (3) (C. Mikida)
c     ##################################################################
cc      function dot_product_3_pcle(A,B) result(C)
      subroutine dot_product_3_pcle(A,B,C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,3), B(_Npg,3) 
      real(_prec) :: C(_Npg)

c     Local variables
      integer :: i,j

      !VECTOR SIMD
      do j=1,_Npg
        C(j) = 0.0d0
      enddo
      
      do i=1, 3
        !VECTOR SIMD
        do j=1,_Npg
            C(j) = C(j) + A(j,i)*B(j,i)
        end do
      end do

      end subroutine dot_product_3_pcle
cc      end function dot_product_3_pcle

c     Optimized matrix-vector multiplication (2) (C. Mikida)
c     ##################################################################
cc      function opt_matvec_2_pcle(A,B) result(C)
      subroutine opt_matvec_2_pcle(A,B,C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,3), B(_Npg,3,2) 
      real(_prec) :: C(_Npg,2)

c     Local variables
      integer :: i,j,k

      !VECTOR SIMD
      do k=1,_Npg
        C(k,:) = 0.0d0
      enddo

      do i=1, 2
         do j=1, 3
            !VECTOR SIMD
            do k=1,_Npg
               C(k,i) = C(k,i) + A(k,j)*B(k,j,i)
            end do
         end do
      end do
 
      end subroutine opt_matvec_2_pcle
cc      end function opt_matvec_2_pcle

c     Optimized matrix-vector multiplication (2) (C. Mikida)
c     ##################################################################
      subroutine opt_matvec_23_pcle(A,B,C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,2,3), B(_Npg,3) 
      real(_prec) :: C(_Npg,2)

c     Local variables
      integer :: i, j, k

      !VECTOR SIMD
      do k=1,_Npg
        C(k,:) = 0.0d0
      enddo

      do i=1, 2
         do j=1, 3
           !VECTOR SIMD
            do k=1,_Npg
               C(k,i) = C(k,i) + B(k,j)*A(k,i,j)
            end do
         end do
      end do
 
      end subroutine opt_matvec_23_pcle

c     Optimized matrix-vector multiplication (2x2) (C. Mikida)
c     ##################################################################
      subroutine opt_matvec_22_pcle(A,B,C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,2), B(_Npg,2,3) 
      real(_prec) :: C(_Npg,3)

c     Local variables
      integer :: i, j, k

      !VECTOR SIMD
      do k=1,_Npg
        C(k,:) = 0.0d0
      enddo
      
      do i=1, 3
         do j=1, 2
           !VECTOR SIMD
            do k=1,_Npg
               C(k,i) = C(k,i) + A(k,j)*B(k,j,i)
            end do
         end do
      end do
      
      end subroutine opt_matvec_22_pcle

c     Optimized matrix-vector multiplication (3) (C. Mikida)
c     ##################################################################
      subroutine opt_matvec_tr_pcle(A,B,C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,3), B(_Npg,3,3) 
      real(_prec) :: C(_Npg,3)

c     Local variables
      integer :: i, j, k

      !VECTOR SIMD
      do k=1,_Npg
        C(k,:) = 0.0d0
      enddo

#if defined(__GNU)
      do i=1, 3
         do j=1, 3
            !VECTOR SIMD
            do k=1,_Npg
               C(k,i) = C(k,i) + A(k,j)*B(k,j,i)
            end do
         end do
      end do
#else
      !VECTOR SIMD
      do k=1,_Npg
         C(k,:) = C(k,:) + matmul(A(k,:),B(k,:,:))
      end do
#endif
      end subroutine opt_matvec_tr_pcle

c     Optimized vector-vector multiplication (3x3) (C. Mikida)
c     ##################################################################
      function opt_vecvec_3_pcle(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,3,1), B(_Npg,1,3) 
      real(_prec) :: C(_Npg,3,3)

c     Local variables
      integer :: i, j, k

      do i=1, 3
         do j=1, 3
            !VECTOR SIMD
            do k=1,_Npg
               C(k,i,j) = A(k,i,1)*B(k,1,j)
            end do
         end do
      end do

      end function opt_vecvec_3_pcle

c     Optimized vector-vector-vector outer product (3x3x3) (Adam)
c     ##################################################################
      function opt_vecvecvec_3_pcle(A,B,C) result(D)
      implicit none

c     Call variables                                                                       
      real(_prec) :: A(_Npg,3,1,1), B(_Npg,1,3,1), C(_Npg,1,1,3)
      real(_prec) :: D(_Npg,3,3,3)

c     Local variables                                                                                                                          
      integer :: i, j, k, ii

      do i=1, 3
         do j=1, 3
            do k=1, 3
               !VECTOR SIMD
               do ii=1,_Npg
                  D(ii,i,j,k) = A(ii,i,1,1)*B(ii,1,j,1)*C(ii,1,1,k)
               end do
            end do
         end do
      end do
               
      end function opt_vecvecvec_3_pcle

c     Optimized vector-vector multiplication (2x3) (C. Mikida)
c     ##################################################################
      function opt_vecvec_23_pcle(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,2,1), B(_Npg,1,3) 
      real(_prec) :: C(_Npg,2,3)

c     Local variables
      integer :: i, j, k

      do i=1, 2
         do j=1, 3
            !VECTOR SIMD
            do k=1,_Npg
               C(k,i,j) = A(k,i,1)*B(k,1,j)
            end do
         end do
      end do

      end function opt_vecvec_23_pcle

c     Optimized vector-vector multiplication (3x2) (C. Mikida)
c     ##################################################################
      function opt_vecvec_32_pcle(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(_Npg,3,1), B(_Npg,1,2) 
      real(_prec) :: C(_Npg,3,2)

c     Local variables
      integer :: i, j, k

      do i=1, 3
         do j=1, 2
            !VECTOR SIMD
            do k=1,_Npg
               C(k,i,j) = A(k,i,1)*B(k,1,j)
            end do
         end do
      end do

      end function opt_vecvec_32_pcle

c     Matrix-matrix (yz) contraction to 2-vector (2x3x3*3x3->2) (Adam)
c     ################################################################
      function opt_matmat_233_yz_pcle(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(_Npg,2,3,3), B(_Npg,3,3)
      real(_prec) :: C(_Npg,2)
c     Local variables
      integer :: i, j, k, ii

      C(:,:) = 0.0d0

      do k=1, 3
         do j=1, 3
            do i=1, 2
               !VECTOR SIMD
               do ii=1,_Npg
                  C(ii,i) = C(ii,i)+A(ii,i,j,k)*B(ii,j,k) ! Sum over y&z.
               end do
            end do
         end do
      end do

      end function opt_matmat_233_yz_pcle

c     Matrix-matrix (zx) contraction to 2-vector (3x2x3*3x3->2) (Adam)
c     ################################################################
      function opt_matmat_323_zx_pcle(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(_Npg,3,2,3), B(_Npg,3,3)
      real(_prec) :: C(_Npg,2)
c     Local variables
      integer :: i, j, k, ii

      C(:,:) = 0.0d0

      do k=1, 3
         do j=1, 2
            do i=1, 3
               !VECTOR SIMD
               do ii=1,_Npg
                  C(ii,j) = C(ii,j)+A(ii,i,j,k)*B(ii,k,i) ! Sum over z&x.
               end do
            end do
         end do
      end do

      end function opt_matmat_323_zx_pcle

c     Matrix-matrix (xy) contraction to 2-vector (3x3x2*3x3->2) (Adam)
c     ################################################################
      function opt_matmat_332_xy_pcle(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(_Npg,3,3,2), B(_Npg,3,3)
      real(_prec) :: C(_Npg,2)
c     Local variables
      integer :: i, j, k, ii

      C(:,:) = 0.0d0

      do k=1, 2
         do j=1, 3
            do i=1, 3
               !VECTOR SIMD
               do ii=1,_Npg
                  C(ii,k) = C(ii,k)+A(ii,i,j,k)*B(ii,i,j) ! Sum over x&y.
               end do
            end do
         end do
      end do

      end function opt_matmat_332_xy_pcle

c     Matrix-matrix (yz) contraction to 3-vector (3x3x3*3x3->3x1x1) (Adam)
c     ################################################################
      function opt_matmat_333_yz_pcle(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(_Npg,3,3,3), B(_Npg,3,3)
      real(_prec) :: C(_Npg,3)
c     Local variables
      integer :: i, j, k, ii

      C(:,:) = 0.0d0

      do k=1, 3
         do j=1, 3
            do i=1, 3
               !VECTOR SIMD
               do ii=1,_Npg
                  C(ii,i) = C(ii,i)+A(ii,i,j,k)*B(ii,j,k) ! Sum over y&z.
               end do
            end do
         end do
      end do

      end function opt_matmat_333_yz_pcle

c     Matrix-matrix (zx) contraction to 3-vector (3x3x3*3x3->3) (Adam)
c     ################################################################
      function opt_matmat_333_zx_pcle(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(_Npg,3,3,3), B(_Npg,3,3)
      real(_prec) :: C(_Npg,3)
c     Local variables
      integer :: i, j, k, ii

      C(:,:) = 0.0d0

      do k=1, 3
         do j=1, 3
            do i=1, 3
               !VECTOR SIMD
               do ii=1,_Npg
                  C(ii,j) = C(ii,j)+A(ii,i,j,k)*B(ii,k,i) ! Sum over z&x.
               end do
            end do
         end do
      end do

      end function opt_matmat_333_zx_pcle

c     Matrix-matrix (xy) contraction to 3-vector (3x3x3*3x3->3) (Adam)
c     ################################################################
      function opt_matmat_333_xy_pcle(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(_Npg,3,3,3), B(_Npg,3,3)
      real(_prec) :: C(_Npg,3)
c     Local variables
      integer :: i, j, k, ii

      C(:,:) = 0.0d0

      do k=1, 3
         do j=1, 3
            do i=1, 3
               !VECTOR SIMD
               do ii=1,_Npg
                  C(ii,k) = C(ii,k)+A(ii,i,j,k)*B(ii,i,j) ! Sum over x&y.
               end do
            end do
         end do
      end do

      end function opt_matmat_333_xy_pcle

c     Matrix-matrix contraction to scalar (3x3x3*3x3x3->1) (Adam)
c     ################################################################
      function opt_matmat_333_xyz_pcle(A,B) result(C)
      implicit none
c     Call variables                                                                                                                                                                                     
      real(_prec) :: A(_Npg,3,3,3), B(_Npg,3,3,3)
      real(_prec) :: C(_Npg)
c     Local variables
                                                                                                                                                        
      integer :: i, j, k, ii

      C(:) = 0.0d0

      do k=1, 3
         do j=1, 3
            do i=1, 3
               !VECTOR SIMD                
               do ii=1,_Npg
                  C(ii) = C(ii)+A(ii,i,j,k)*B(ii,i,j,k) 
               end do
            end do
         end do
      end do

      end function opt_matmat_333_xyz_pcle

c     Matrix-matrix contraction to scalar (3x3x3*3x3x3->1) (Adam)
c     ################################################################
      function opt_matmat_333_xyz(A,B) result(C)
      implicit none
c     Call variables
      real(_prec) :: A(3,3,3), B(3,3,3)
      real(_prec) :: C

c     Local variables
      integer :: i, j, k

      C = 0.0d0

      do k=1, 3
         do j=1, 3
            do i=1, 3
                  C = C+A(i,j,k)*B(i,j,k)
            end do
         end do
      end do

      end function opt_matmat_333_xyz

c     Optimized matrix-vector multiplication (2) (C. Mikida)
c     ##################################################################
      function opt_matvec_2(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(3), B(3,2) 
      real(_prec) :: C(2)

c     Local variables
      integer :: i, j

      C = 0.0d0

      do i=1, 2
         !VECTOR SIMD REDUCTION(+:C)
         do j=1, 3
            C(i) = C(i) + A(j)*B(j,i)
         end do
      end do
 
      end function opt_matvec_2

c     Optimized matrix-vector multiplication (2) (C. Mikida)
c     ##################################################################
      function opt_matvec_23(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(2,3), B(3) 
      real(_prec) :: C(2)

c     Local variables
      integer :: i, j

      C = 0.0d0

      do i=1, 2
         !VECTOR SIMD REDUCTION(+:C)
         do j=1, 3
            C(i) = C(i) + B(j)*A(i,j)
         end do
      end do
 
      end function opt_matvec_23

c     Optimized matrix-vector multiplication (2x2) (C. Mikida)
c     ##################################################################
      function opt_matvec_22(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(2), B(2,3) 
      real(_prec) :: C(3)

c     Local variables
      integer :: i, j

      C = 0.0d0

      do i=1, 3
         !VECTOR SIMD REDUCTION(+:C)
         do j=1, 2
            C(i) = C(i) + A(j)*B(j,i)
         end do
      end do
      
      end function opt_matvec_22

c     Optimized matrix-vector multiplication (3) (C. Mikida)
c     ##################################################################
      function opt_matvec_tr(A,B) result(C)

      implicit none

c     Call variables

      real(_prec) :: A(3), B(3,3) 
      real(_prec) :: C(3)

c     Local variables

      integer :: i, j

      C(:) = 0.0d0

      do i=1, 3
         !VECTOR SIMD REDUCTION(+:C)
         do j=1, 3
            C(i) = C(i) + A(j)*B(j,i)
         end do
      end do
 
      end function opt_matvec_tr

c     Optimized matrix-vector multiplication
c     ##################################################################
      function opt_matvec(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: B(3), A(3,3) 
      real(_prec) :: C(3)

c     Local variables
      integer :: i, j

      C(:) = 0.0d0

      do i=1, 3
         !VECTOR SIMD REDUCTION(+:C)
         do j=1, 3
            C(i) = C(i) + A(i,j)*B(j)
         end do
      end do
 
      end function opt_matvec

c     Optimized vector-vector multiplication (3x3) (C. Mikida)
c     ##################################################################
      function opt_vecvec_3(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(3,1), B(1,3) 
      real(_prec) :: C(3,3)

c     Local variables
      integer :: i, j

      C(:,:) = 0.0d0

      do i=1, 3
         !VECTOR SIMD
         do j=1, 3
            C(i,j) = A(i,1)*B(1,j)
         end do
      end do

      end function opt_vecvec_3

c     Outer product of 3 vectors (3x3x3) (Adam)
c     ##################################################################
      function opt_vecvecvec_3(A,B,C) result(D)

      implicit none
c     Call variables

      real(_prec) :: A(3,1,1), B(1,3,1), C(1,1,3)
      real(_prec) :: D(3,3,3)

c     Local variables 
      integer :: i, j, k

      D(:,:,:) = 0.0d0

      do i=1, 3
         do j=1, 3
            do k=1, 3
               D(i,j,k) = A(i,1,1)*B(1,j,1)*C(1,1,k)
            end do
         end do
      end do

      end function opt_vecvecvec_3

c     Optimized vector-vector multiplication (2x3) (C. Mikida)
c     ##################################################################
      function opt_vecvec_23(A,B) result(C)

      implicit none

c     Call variables
      real(_prec) :: A(2,1), B(1,3) 
      real(_prec) :: C(2,3)

c     Local variables
      integer :: i, j

      C(:,:) = 0.0d0

      do i=1, 2
         !VECTOR SIMD
         do j=1, 3
            C(i,j) = A(i,1)*B(1,j)
         end do
      end do

      end function opt_vecvec_23

c     Outer product of three vectors (2x3x3) (Adam)
c     ##################################################################
      function opt_vecvecvec_233(A,B,C) result(D)

      implicit none
c     Call variables
      real(_prec) :: A(2,1,1), B(1,3,1), C(1,1,3)
      real(_prec) :: D(2,3,3)

c     Local variables
      integer :: i, j, k

      D(:,:,:) = 0.0d0

      do i=1, 2
         do j=1, 3
         !VECTOR SIMD
            do k=1, 3
               D(i,j,k) = A(i,1,1)*B(1,j,1)*C(1,1,k)
            end do
         end do
      end do
      
      end function opt_vecvecvec_233

c     Optimized vector-vector multiplication (3x2) (C. Mikida)
c     ##################################################################
      function opt_vecvec_32(A,B) result(C)

      implicit none

c     Call variables

      real(_prec) :: A(3,1), B(1,2) 
      real(_prec) :: C(3,2)

c     Local variables
      integer :: i, j

      C(:,:) = 0.0d0

      do i=1, 3
         !VECTOR SIMD
         do j=1, 2
            C(i,j) = A(i,1)*B(1,j)
         end do
      end do

      end function opt_vecvec_32

c     Outer product of three vectors (3x2x3) (Adam)
c     ##################################################################
      function opt_vecvecvec_323(A,B,C) result(D)

      implicit none
c     Call variables                                                                                                                           
      real(_prec) :: A(3,1,1), B(1,2,1), C(1,1,3)
      real(_prec) :: D(3,2,3)

c     Local variables
      integer :: i, j, k

      D(:,:,:) = 0.0d0

      do i=1, 3
         do j=1, 2
         !VECTOR SIMD
            do k=1, 3
               D(i,j,k) = A(i,1,1)*B(1,j,1)*C(1,1,k)
            end do
         end do
      end do


      end function opt_vecvecvec_323

c     Outer product of 3 vectors (3x3x2) (Adam)
c     ##################################################################
      function opt_vecvecvec_332(A,B,C) result(D)

      implicit none
c     Call variables                                                                                                                           
      real(_prec) :: A(3,1,1), B(1,3,1), C(1,1,2)
      real(_prec) :: D(3,3,2)

c     Local variables
      integer :: i, j, k

      D(:,:,:) = 0.0d0

      do i=1, 3
         do j=1, 3
         !VECTOR SIMD
            do k=1, 2
               D(i,j,k) = A(i,1,1)*B(1,j,1)*C(1,1,k)
            end do
         end do
      end do

      end function opt_vecvecvec_332

c     Matrix-vector outer product (2x3x3) (Adam)
c     ################################################################
      function opt_omatvec_233(A,B) result(C)

      implicit none
c     Call variables                                                                                                                           
      real(_prec) :: A(3,3), B(2,1,1)
      real(_prec) :: C(2,3,3)

c     Local variables
      integer :: i, j, k

      C(:,:,:) = 0.0d0

      do k=1, 3
         do j=1, 3
         !VECTOR SIMD
            do i=1, 2
               C(i,j,k) = B(i,1,1)*A(j,k)
            end do
         end do
      end do

      end function opt_omatvec_233

c     Matrix-vector outer product (3x2x3) (Adam)
c     ################################################################
      function opt_omatvec_323(A,B) result(C)

      implicit none
c     Call variables
      real(_prec) :: A(3,3), B(1,2,1)
      real(_prec) :: C(3,2,3)

c     Local variables
      integer :: i, j, k

      C(:,:,:) = 0.0d0

      do k=1, 3
         do j=1, 2
         !VECTOR SIMD
            do i=1, 3
               C(i,j,k) = B(1,j,1)*A(k,i) ! Need in xyz, from A_zx.
            end do
         end do
      end do

      end function opt_omatvec_323

c     Matrix-vector outer product (3x3x2) (Adam)
c     ################################################################
      function opt_omatvec_332(A,B) result(C)

      implicit none
c     Call variables
      real(_prec) :: A(3,3), B(1,1,2)
      real(_prec) :: C(3,3,2)

c     Local variables
      integer :: i, j, k

      C(:,:,:) = 0.0d0

      do k=1, 2
         do j=1, 3
         !VECTOR SIMD
            do i=1, 3
               C(i,j,k) = B(1,1,k)*A(i,j) ! Need in xyz, from A_xy.
            end do
         end do
      end do

      end function opt_omatvec_332

      end module pic_shape_f

c module pic_map
c ######################################################################
      module pic_map

      use math

      use grid_operations

      use grid_splines

      use pic_shape_f

      use pic_BCs

      real(8),pointer,dimension(:,:,:)   :: jac_g
      real(8),pointer,dimension(:,:,:,:) :: map_g
      integer,private :: nxl,nyl,nzl,nnl,nxg,nyg,nzg,nng
      real(8),private :: ihx,ihy,ihz,h_x,h_y,h_z
      logical,private :: global_splines=.false.
      logical,protected :: no_map=.true.

      real(8),private :: x_min,x_max,y_min,y_max,z_min,z_max

      contains

c     setup_map
c     ##################################################################
      subroutine setup_map(nxd,nyd,nzd,jac,no__map,xmap,
     .     xin,yin,zin,xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz)

      implicit none

c     Call variables
      integer :: nxd,nyd,nzd
      real(8), intent(in) :: jac   (0:,0:,0:)
     .                      ,xmap  (0:,0:,0:,:)
     .                      ,xin(0:),yin(0:),zin(0:)
      real(8), intent(out) :: hx,hy,hz
      real(8), intent(inout):: xmin,xmax,ymin,ymax,zmin,zmax
      logical, intent(in) :: no__map

      no_map  = no__map

      nxg = nxd
      nyg = nyd
      nzg = nzd

      call alloc_map(nxg,nyg,nzg,nxl,nyl,nzl)

      jac_g  = jac (:,:,:)
      map_g  = xmap(:,:,:,:)

      if (.not.no_map.and.global_splines) then
        call setupSplines(nxg+2,nyg+2,nzg+2,xin,yin,zin,3,bcond)
        call splineJ(jac )
        call splineX(xmap)
      endif

      !Set limits to faces
      xmax = 0.5*(xin(nxg)+xin(nxg+1))
      xmin = 0.5*(xin(0)  +xin(1)    )
      ymax = 0.5*(yin(nyg)+yin(nyg+1))
      ymin = 0.5*(yin(0)  +yin(1)    )
      zmax = 0.5*(zin(nzg)+zin(nzg+1))
      zmin = 0.5*(zin(0)  +zin(1)    )

      hx  = (xmax-xmin)/dble(nxg)
      hy  = (ymax-ymin)/dble(nyg)
      hz  = (zmax-zmin)/dble(nzg)

      !Assign internal variables

      h_x = hx
      h_y = hy
      h_z = hz

      ihx = 1d0/hx
      ihy = 1d0/hy
      ihz = 1d0/hz

      x_min = xmin
      x_max = xmax
      y_min = ymin
      y_max = ymax
      z_min = zmin
      z_max = zmax

      end subroutine setup_map

c     get_inverse_h
c     ##################################################################
      subroutine get_inverse_h(ihx0,ihy0,ihz0)

      implicit none

      real(8) :: ihx0,ihy0,ihz0

      if (ihx == 0d0 .or. ihy == 0d0 .or. ihz == 0d0) then
         write (*,*) "Inverse spacing undefined"
         stop
      endif

      ihx0 = ihx
      ihy0 = ihy
      ihz0 = ihz

      end subroutine get_inverse_h

c     get_h
c     ##################################################################
      subroutine get_h(hx0,hy0,hz0)

      implicit none

      real(8) :: hx0,hy0,hz0

      if (h_x == 0d0 .or. h_y == 0d0 .or. h_z == 0d0) then
         write (*,*) "Inverse spacing undefined"
         stop
      endif

      hx0 = h_x
      hy0 = h_y
      hz0 = h_z

      end subroutine get_h

c     get_global_grid_size
c     ##################################################################
      subroutine get_global_grid_size(nx_g,ny_g,nz_g)

      implicit none

      integer :: nx_g,ny_g,nz_g

      nx_g = nxg
      ny_g = nyg
      nz_g = nzg

      end subroutine get_global_grid_size

c     alloc_map
c     ##################################################################
      subroutine alloc_map(nxd,nyd,nzd,nx,ny,nz)

      implicit none

c     Call variables

      integer :: nxd,nyd,nzd,nx,ny,nz

c     Local variables

c     Begin program

      allocate(jac_g(0:nxg+1,0:nyg+1,0:nzg+1))
      allocate(map_g(0:nxg+1,0:nyg+1,0:nzg+1,3))

c     End program

      end subroutine alloc_map

c     dealloc_map
c     ##################################################################
      subroutine dealloc_map

      deallocate(jac_g,map_g)

      end subroutine dealloc_map

c     find_vcurv
c     ##################################################################
      function find_vcurv(i,j,k,xp,vp,gamma) result(vcurv) 

      implicit none

c     Call variables

      integer :: i(_Npg),j(_Npg),k(_Npg)
      real(_prec) :: xp(_Npg,3)
      real(_prec) :: vp(_Npg,3)
      real(_prec) :: vcurv(_Npg,3)
      real(_prec),optional :: gamma(_Npg,3,3,3)

c     Local variables

      real(_prec) :: xlog(_Npg,3),cov(_Npg,3,3)
      integer :: ii

c     Begin program

      xlog(:,1) = xp(:,1)*ihx
      xlog(:,2) = xp(:,2)*ihy
      xlog(:,3) = xp(:,3)*ihz

      call transf_vector_Bspl(i,j,k,xlog,cov,gamma=gamma)

      !VECTOR SIMD
      do ii=1,_Npg  
#ifdef __GNU
        vcurv(ii,:) = opt_matvec(cov(ii,:,:),vp(ii,:))
#else               
        vcurv(ii,:) = matmul    (cov(ii,:,:),vp(ii,:))
#endif
      end do

c     End program

      end function find_vcurv

c     transf_vector_Bspl
c     #################################################################
      subroutine transf_vector_Bspl(i0,j0,k0,xlog,cov,cnv,gamma)

c     -----------------------------------------------------------------
c     Transform Cartesian vectors to curvilinear using B-splines.
c     Employ Luis's convention.
c     -----------------------------------------------------------------

      implicit none

c     Input variables
     
      integer, intent(in) :: i0(_Npg),j0(_Npg),k0(_Npg) !particle cell-indices
      real(_prec) :: xlog(_Npg,3) !normalized particle offset position
      real(_prec),INTENT(OUT) :: cov(_Npg,3,3)
      real(_prec),optional :: cnv(_Npg,3,3),gamma(_Npg,3,3,3) !Christoffel symbol of 2nd kind

c     Local variables

      integer :: ierr,ii

      real(8) :: vec(_Npg,3,3),iJac(_Npg) !dum3(_Npg,3)
      real(8) :: Sx2(_Npg,3,1,1),Sy2(_Npg,1,3,1),Sz2(_Npg,1,1,3)
      real(8) :: DSx(_Npg,3,1,1),DSy(_Npg,1,3,1),DSz(_Npg,1,1,3)
      real(8) :: DSxx(_Npg,3),DSyy(_Npg,3),DSzz(_Npg,3)
      real(8) :: map_g_1(_Npg,3,3,3),map_g_2(_Npg,3,3,3),
     .     map_g_3(_Npg,3,3,3)
      real(8) :: dum333(_Npg,3,3,3)

c     Begin program

c     Tangential vectors

      call shape2_offset(xlog(:,1),Sx2(:,:,1,1))
      call shape2_offset(xlog(:,2),Sy2(:,1,:,1))
      call shape2_offset(xlog(:,3),Sz2(:,1,1,:))

      DSx(:,:,1,1) = Dx_shape2_offset(xlog(:,1))
      DSy(:,1,:,1) = Dy_shape2_offset(xlog(:,2))
      DSz(:,1,1,:) = Dz_shape2_offset(xlog(:,3))

      !VECTOR SIMD
      do ii=1,_Npg
         map_g_1(ii,:,:,:) = map_g(i0(ii)-1:i0(ii)+1,j0(ii)-1:j0(ii)+1,
     .                             k0(ii)-1:k0(ii)+1,1)
         map_g_2(ii,:,:,:) = map_g(i0(ii)-1:i0(ii)+1,j0(ii)-1:j0(ii)+1,
     .                             k0(ii)-1:k0(ii)+1,2)
         map_g_3(ii,:,:,:) = map_g(i0(ii)-1:i0(ii)+1,j0(ii)-1:j0(ii)+1,
     .                             k0(ii)-1:k0(ii)+1,3)
      enddo

      ! Adam: vec(1:Npg,:,:) = (dx/dxi, dy/dxi, dz/dxi; dx/deta, dy/deta, dz/deta; dx/dzeta, dy/dzeta, dz/dzeta) 
      dum333(:,:,:,:) = opt_vecvecvec_3_pcle(DSx,Sy2,Sz2)
      vec(:,1,1) = opt_matmat_333_xyz_pcle(dum333,map_g_1)
      vec(:,1,2) = opt_matmat_333_xyz_pcle(dum333,map_g_2)
      vec(:,1,3) = opt_matmat_333_xyz_pcle(dum333,map_g_3)

      dum333(:,:,:,:) = opt_vecvecvec_3_pcle(Sx2,DSy,Sz2)
      vec(:,2,1) = opt_matmat_333_xyz_pcle(dum333,map_g_1)
      vec(:,2,2) = opt_matmat_333_xyz_pcle(dum333,map_g_2)
      vec(:,2,3) = opt_matmat_333_xyz_pcle(dum333,map_g_3)

      dum333(:,:,:,:) = opt_vecvecvec_3_pcle(Sx2,Sy2,DSz)
      vec(:,3,1) = opt_matmat_333_xyz_pcle(dum333,map_g_1)
      vec(:,3,2) = opt_matmat_333_xyz_pcle(dum333,map_g_2)
      vec(:,3,3) = opt_matmat_333_xyz_pcle(dum333,map_g_3)

!      call opt_matvec_tr_pcle(DSx,map_g_1,dum3)
!      call dot_product_3_pcle(dum3,Sy2,vec(:,1,1))
!      call opt_matvec_tr_pcle(DSx,map_g_2,dum3)
!      call dot_product_3_pcle(dum3,Sy2,vec(:,1,2))
!      vec(:,1,3)=0d0

!      call opt_matvec_tr_pcle(Sx2,map_g_1,dum3)
!      call dot_product_3_pcle(dum3,Dsy,vec(:,2,1))
!      call opt_matvec_tr_pcle(Sx2,map_g_2,dum3)
!      call dot_product_3_pcle(dum3,DSy,vec(:,2,2))
!      vec(:,2,3)=0d0

!      !(careful) z is orthogonal to the (x,y) plane
!      vec(:,3,1) = 0d0
!      vec(:,3,2) = 0d0
!      vec(:,3,3) = 1d0

c     Find the inverse Jacobian

      call idet(vec,iJac)

c     Divide tangential vectors by Jacobian (Luis's convention)

      if(present(cnv)) then
        !VECTOR SIMD
        do ii=1,_Npg
          cnv(ii,:,:) = vec(ii,:,:)*iJac(ii)
        enddo
      end if

c     Normal vectors (J*grad(xi)=(dx/deta x dx/dzeta))
      call cross_product_vec(vec(:,2,:),vec(:,3,:),cov(:,1,:))
      call cross_product_vec(vec(:,3,:),vec(:,1,:),cov(:,2,:))
      call cross_product_vec(vec(:,1,:),vec(:,2,:),cov(:,3,:))

c     Christoffel symbols
      if(present(gamma)) call compute_gamma

c     Find normal vectors

      !VECTOR SIMD
      do ii=1,_Npg
        cov(ii,:,:) = cov(ii,:,:)*iJac(ii)
      enddo

      contains

c     compute_gamma
c     #################################################################
      subroutine compute_gamma ! ADAM: NOT UPDATED TO 3D. NOT CURRENTLY USED.

      implicit none

c     -----------------------------------------------------------------
c     Compute Christoffel symbols
c     -----------------------------------------------------------------

      gamma = 0d0

c$$$
c$$$      !Gamma(:,1,2)
c$$$      DSxx = Dxx_shape2_offset(xlog(:,1))
c$$$      DSyy = Dyy_shape2_offset(xlog(:,2))
c$$$c$$$         DSxx(1,:)=Dxx_shape2_offset_scalar(xlog(1))
c$$$c$$$         DSyy(1,:)=Dyy_shape2_offset_scalar(xlog(2))
c$$$c$$$
c$$$c$$$cc         call evalX(xlog(1),xlog(2),xlog(3),x,y,z
c$$$c$$$cc     .             ,ierr,xder=1,yder=1,zder=0)  
c$$$
c$$$      call opt_matvec_tr_pcle(DSx,map_g_1,dum3)
c$$$      call dot_product_3_pcle(dum3,Dsy,vec(:,1,1))
c$$$      call opt_matvec_tr_pcle(DSx,map_g_2,dum3)
c$$$      call dot_product_3_pcle(dum3,Dsy,vec(:,1,2))
c$$$      vec(:,1,3)=0d0
c$$$
c$$$      !VECTOR SIMD
c$$$      do k=1,_Npg
c$$$#ifdef __GNU
c$$$        gamma(k,:,1,2) = opt_matvec(cov(k,:,:),vec(k,1,:))
c$$$#else
c$$$        gamma(k,:,1,2) = matmul    (cov(k,:,:),vec(k,1,:))
c$$$#endif
c$$$        gamma(k,:,2,1) = gamma(k,:,1,2)
c$$$      enddo
c$$$c$$$#ifdef __GNU
c$$$c$$$         x = dot_product(opt_matvec_tr(DSx(1,:),
c$$$c$$$     .                  map_g_1),DSy(1,:))
c$$$c$$$         y = dot_product(opt_matvec_tr(DSx(1,:),
c$$$c$$$     .                  map_g_2),DSy(1,:))
c$$$c$$$#else
c$$$c$$$         x = dot_product(MATMUL(DSx(1,:),map_g_1),DSy(1,:))
c$$$c$$$         y = dot_product(MATMUL(DSx(1,:),map_g_2),DSy(1,:))
c$$$c$$$#endif
c$$$c$$$         z = 0d0
c$$$c$$$         gamma(:,1,2) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
c$$$c$$$         gamma(:,2,1) = gamma(:,1,2)
c$$$
c$$$         !Gamma(:,1,3)
c$$$c$$$c$$$         call evalX(xlog(1),xlog(2),xlog(3),x,y,z
c$$$c$$$c$$$     .             ,ierr,xder=1,yder=0,zder=1)  
c$$$c$$$         x = 0d0 !dot_product(MATMUL(DSx,map_g(i0-1:i0+1,j0-1:j0+1,1)),Sy2)
c$$$c$$$         y = 0d0 !dot_product(MATMUL(DSx,map_g(i0-1:i0+1,j0-1:j0+1,2)),Sy2)
c$$$c$$$         z = 0d0
c$$$c$$$         
c$$$c$$$         gamma(:,1,3) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
c$$$c$$$         gamma(:,3,1) = gamma(:,1,3)
c$$$c$$$
c$$$         !Gamma(:,2,3)
c$$$c$$$c$$$         call evalX(xlog(1),xlog(2),xlog(3),x,y,z
c$$$c$$$c$$$     .             ,ierr,xder=0,yder=1,zder=1)  
c$$$c$$$         x = 0d0 !dot_product(MATMUL(Sx2,map_g(i0-1:i0+1,j0-1:j0+1,1)),DSy)
c$$$c$$$         y = 0d0 !dot_product(MATMUL(Sx2,map_g(i0-1:i0+1,j0-1:j0+1,2)),DSy)
c$$$c$$$         z = 0d0
c$$$c$$$
c$$$c$$$         gamma(:,2,3) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
c$$$c$$$         gamma(:,3,2) = gamma(:,2,3)
c$$$
c$$$      !Gamma(:,1,1)
c$$$      call opt_matvec_tr_pcle(DSxx,map_g_1,dum3)
c$$$      call dot_product_3_pcle(dum3,Sy2,vec(:,1,1))
c$$$      call opt_matvec_tr_pcle(DSxx,map_g_2,dum3)
c$$$      call dot_product_3_pcle(dum3,Sy2,vec(:,1,2))
c$$$      vec(:,1,3)=0d0
c$$$      !VECTOR SIMD
c$$$      do k=1,_Npg
c$$$#ifdef __GNU
c$$$        gamma(k,:,1,1) = opt_matvec(cov(k,:,:),vec(k,1,:))
c$$$#else
c$$$        gamma(k,:,1,1) = matmul    (cov(k,:,:),vec(k,1,:))
c$$$#endif
c$$$      enddo
c$$$c$$$#ifdef __GNU
c$$$c$$$         x = dot_product(opt_matvec_tr(DSxx(1,:),
c$$$c$$$     .                   map_g_1),Sy2(1,:))
c$$$c$$$         y = dot_product(opt_matvec_tr(DSxx(1,:),
c$$$c$$$     .                   map_g_2),Sy2(1,:))
c$$$c$$$#else
c$$$c$$$         x = dot_product(MATMUL(DSxx(1,:),map_g_1),Sy2(1,:))
c$$$c$$$         y = dot_product(MATMUL(DSxx(1,:),map_g_2),Sy2(1,:))
c$$$c$$$#endif
c$$$c$$$         z = 0d0
c$$$c$$$         gamma(:,1,1) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
c$$$
c$$$      !Gamma(:,2,2)
c$$$      if(nyg>1) then
c$$$        call opt_matvec_tr_pcle(Sx2,map_g_1,dum3)
c$$$        call dot_product_3_pcle(dum3,DSyy,vec(:,2,1))
c$$$        call opt_matvec_tr_pcle(Sx2,map_g_2,dum3)
c$$$        call dot_product_3_pcle(dum3,DSyy,vec(:,2,2))
c$$$        vec(:,2,3)=0d0
c$$$        !VECTOR SIMD
c$$$        do k=1,_Npg
c$$$#ifdef __GNU
c$$$          gamma(k,:,2,2) = opt_matvec(cov(k,:,:),vec(k,2,:))
c$$$#else
c$$$          gamma(k,:,2,2) = matmul    (cov(k,:,:),vec(k,2,:))
c$$$#endif
c$$$        enddo
c$$$c$$$c$$$            call evalX(xlog(1),xlog(2),xlog(3),x,y,z
c$$$c$$$c$$$     .           ,ierr,xder=0,yder=2,zder=0)  
c$$$c$$$#ifdef __GNU
c$$$c$$$          x = dot_product(opt_matvec_tr(Sx2(1,:),
c$$$c$$$     .                    map_g_1),DSyy(1,:))
c$$$c$$$          y = dot_product(opt_matvec_tr(Sx2(1,:),
c$$$c$$$     .                    map_g_2),DSyy(1,:))
c$$$c$$$#else
c$$$c$$$          x = dot_product(MATMUL(Sx2(1,:),map_g_1),DSyy(1,:))
c$$$c$$$          y = dot_product(MATMUL(Sx2(1,:),map_g_2),DSyy(1,:))
c$$$c$$$#endif
c$$$c$$$          z = 0d0            
c$$$c$$$          gamma(:,2,2) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
c$$$      end if

      end subroutine compute_gamma

      end subroutine transf_vector_Bspl

c     transf_vector_Bspl_scalar
c     #################################################################
      subroutine transf_vector_Bspl_scalar(i0,j0,k0,xlog,cov,cnv,jac) 

c     -----------------------------------------------------------------
c     Transform Cartesian vectors to curvilinear using B-splines.
c     Employ Luis's convention.
c     -----------------------------------------------------------------

      implicit none

c     Input variables
     
      integer, intent(in) :: i0,j0,k0 !particle cell-indices
      real(_prec) :: xlog(3) !normalized particle offset position
      real(_prec),INTENT(OUT) :: cov(3,3)
      real(_prec),optional :: cnv(3,3),jac

c     Local variables

      integer :: ierr,k

      real(_prec) :: vec(3,3),dum333(3,3,3),iJac,jacob
      real(_prec) :: Sx2(3,1,1),Sy2(1,3,1),Sz2(1,1,3)
      real(_prec) :: DSx(3,1,1),DSy(1,3,1),DSz(1,1,3)
      real(_prec) :: DSxx(3),DSyy(3),DSzz(3)
      real(_prec) :: map_g_1(3,3,3),map_g_2(3,3,3),map_g_3(3,3,3)

c     Begin program

c     Tangential vectors (2D for now)

      Sx2(:,1,1) = shape2_offset_scalar(xlog(1))
      Sy2(1,:,1) = shape2_offset_scalar(xlog(2))
      Sz2(1,1,:) = shape2_offset_scalar(xlog(3))

      DSx(:,1,1) = Dx_shape2_offset_scalar(xlog(1))
      DSy(1,:,1) = Dy_shape2_offset_scalar(xlog(2))
      DSz(1,1,:) = Dz_shape2_offset_scalar(xlog(3))

      map_g_1 = map_g(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,1)
      map_g_2 = map_g(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,2)
      map_g_3 = map_g(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,3)

      dum333(:,:,:) = opt_vecvecvec_3(DSx,Sy2,Sz2)
      vec(1,1) = opt_matmat_333_xyz(dum333,map_g_1)
      vec(1,2) = opt_matmat_333_xyz(dum333,map_g_2)
      vec(1,3) = opt_matmat_333_xyz(dum333,map_g_3)      

      dum333(:,:,:) = opt_vecvecvec_3(Sx2,DSy,Sz2)
      vec(2,1) = opt_matmat_333_xyz(dum333,map_g_1)
      vec(2,2) = opt_matmat_333_xyz(dum333,map_g_2)
      vec(2,3) = opt_matmat_333_xyz(dum333,map_g_3)

      dum333(:,:,:) = opt_vecvecvec_3(Sx2,Sy2,DSz)
      vec(3,1) = opt_matmat_333_xyz(dum333,map_g_1)
      vec(3,2) = opt_matmat_333_xyz(dum333,map_g_2)
      vec(3,3) = opt_matmat_333_xyz(dum333,map_g_3)

!      dum3 = opt_matvec_tr(DSx,map_g_1)
!      vec(1,1) = dot_product(dum3,Sy2)
!      dum3 = opt_matvec_tr(DSx,map_g_2)
!      vec(1,2) = dot_product(dum3,Sy2)
!      vec(1,3)=0.0

!      dum3 = opt_matvec_tr(Sx2,map_g_1)
!      vec(2,1) = dot_product(dum3,Dsy)
!      dum3 = opt_matvec_tr(Sx2,map_g_2)
!      vec(2,2) = dot_product(dum3,DSy)
!      vec(2,3)=0.0

!      !(careful) z is orthogonal to the (x,y) plane
!      vec(3,1) = 0.0
!      vec(3,2) = 0.0
!      vec(3,3) = 1.0

c     Find the inverse Jacobian

      jacob = det(vec)

      if (PRESENT(jac)) jac = jacob

c     Divide tangential vectors by Jacobian (Luis's convention)

      iJac = 1.0/jacob

      if(present(cnv)) cnv = vec*iJac

c     Normal vectors (J*grad(xi)=(dx/deta x dx/dzeta))

      cov(1,:) = cross_product(vec(2,:),vec(3,:))
      cov(2,:) = cross_product(vec(3,:),vec(1,:))
      cov(3,:) = cross_product(vec(1,:),vec(2,:))

      cov = cov*iJac

      end subroutine transf_vector_Bspl_scalar

c     cross_product_vec
c     #################################################################
      subroutine cross_product_vec(vec1,vec2,vec3)

c     -----------------------------------------------------------------
c     Perform cross product of vectors vec1 and vec2: vec3 = vec1 x vec2,
c     where the vectors are in Cartesian coordinates. Vectorized.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: vec1(_Npg,3),vec2(_Npg,3),vec3(_Npg,3)

c     Local variables

        integer :: k

c     Begin program

        !VECTOR SIMD
        do k=1,_Npg
          vec3(k,1) = vec1(k,2)*vec2(k,3)-vec1(k,3)*vec2(k,2)
          vec3(k,2) = vec1(k,3)*vec2(k,1)-vec1(k,1)*vec2(k,3)
          vec3(k,3) = vec1(k,1)*vec2(k,2)-vec1(k,2)*vec2(k,1)
        enddo

      end subroutine cross_product_vec

c     idet
c     #################################################################
      subroutine idet(tensor,ideter)

c     -----------------------------------------------------------------
c     Compute determinant (vectorized)
c     -----------------------------------------------------------------

        implicit none

        real(8) :: tensor(_Npg,3,3)

        real(8) :: ideter(_Npg)

c     Local variables

        integer :: k
        real(8) :: det

c     Begin program

        !VECTOR SIMD
        do k=1,_Npg
          det = tensor(k,1,1)*tensor(k,2,2)*tensor(k,3,3)
     .        + tensor(k,3,2)*tensor(k,2,1)*tensor(k,1,3)
     .        + tensor(k,1,2)*tensor(k,2,3)*tensor(k,3,1)
     .        - tensor(k,1,3)*tensor(k,2,2)*tensor(k,3,1)
     .        - tensor(k,1,1)*tensor(k,2,3)*tensor(k,3,2)
     .        - tensor(k,3,3)*tensor(k,1,2)*tensor(k,2,1)
          ideter(k) = 1d0/det
        enddo

      end subroutine idet

c     cov_vector_spl
c     #################################################################
      function cov_vector_spl(curv) result(covv)

      implicit none

c     Input variables
      
      real(8) :: curv(3),covv(3,3)

c     Local variables

      real(8) :: vec(3,3)
      
      real(8) :: x,y,z, Jac
      integer :: ierr

c     Begin program

!     standard covriant base vector
      call evalX(curv(1),curv(2),curv(3),x,y,z
     .          ,ierr,xder=1,yder=0,zder=0) !x=dx/dxi,y=dy/dxi,z=dz/dxi
      vec(1,:) = (/ x, y, z /)
      call evalX(curv(1),curv(2),curv(3),x,y,z
     .          ,ierr,xder=0,yder=1,zder=0) 
      vec(2,:) = (/ x, y, z /)
      call evalX(curv(1),curv(2),curv(3),x,y,z
     .          ,ierr,xder=0,yder=0,zder=1) 
      vec(3,:) = (/ x, y, z /)

!     find the Jacobian
      Jac = det(vec)

      covv(1,:) = cross_product(vec(2,:),vec(3,:))
      covv(2,:) = cross_product(vec(3,:),vec(1,:))
      covv(3,:) = cross_product(vec(1,:),vec(2,:))

      covv = covv/Jac             !Luis's convention

!      covv(3,1) = 0d0 ! ADAM: COMMENTED OUT - THINK THIS WAS TO HARDCODE e_z in 2D.
!      covv(3,2) = 0d0 ! ACTUALLY THIS FUNCTION DOESN'T SEEM TO BE USED.
!      covv(3,3) = 1d0

      end function cov_vector_spl

c     transf_vector_spl
c     #################################################################
      subroutine transf_vector_spl(curv,cov,cnv,gamma)
!     Luis's convention

      implicit none

c     Input variables
      
      real(8) :: curv(3),cov(3,3)
      real(8), optional :: cnv(3,3)            !Contravariant vector basis
     .                    ,gamma(3,3,3)        !Christoffel symbol of 2nd kind

c     Local variables

      real(8) :: vec(3,3)
      
      real(8) :: x,y,z, Jac
      integer :: ierr

c     Begin program

      if(curv(1)>x_max) curv(1) = x_max
      if(curv(2)>y_max) curv(2) = y_max
      if(curv(3)>z_max) curv(3) = z_max

      if(curv(1)<x_min) curv(1) = x_min
      if(curv(2)<y_min) curv(2) = y_min
      if(curv(3)<z_min) curv(3) = z_min

!     standard covriant base vector
      call evalX(curv(1),curv(2),curv(3),x,y,z
     .          ,ierr,xder=1,yder=0,zder=0) !x=dx/dxi,y=dy/dxi,z=dz/dxi
      vec(1,:) = (/ x, y, z /)
      call evalX(curv(1),curv(2),curv(3),x,y,z
     .          ,ierr,xder=0,yder=1,zder=0) 
      vec(2,:) = (/ x, y, z /)
      call evalX(curv(1),curv(2),curv(3),x,y,z
     .          ,ierr,xder=0,yder=0,zder=1) 
      vec(3,:) = (/ x, y, z /)

!     find the Jacobian
      Jac = det(vec)

      cov(1,:) = cross_product(vec(2,:),vec(3,:))
      cov(2,:) = cross_product(vec(3,:),vec(1,:))
      cov(3,:) = cross_product(vec(1,:),vec(2,:))

      if(present(cnv)) then
         cnv = vec/Jac          !Luis's convention
      end if

      if(present(gamma)) then
         call evalX(curv(1),curv(2),curv(3),x,y,z
     .             ,ierr,xder=1,yder=1,zder=0)  
         gamma(:,1,2) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
         call evalX(curv(1),curv(2),curv(3),x,y,z
     .             ,ierr,xder=1,yder=0,zder=1)  
         gamma(:,1,3) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
         call evalX(curv(1),curv(2),curv(3),x,y,z
     .             ,ierr,xder=0,yder=1,zder=1)  
         gamma(:,2,3) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)

         gamma(:,2,1) = gamma(:,1,2)
         gamma(:,3,1) = gamma(:,1,3)
         gamma(:,3,2) = gamma(:,2,3)

         call evalX(curv(1),curv(2),curv(3),x,y,z
     .             ,ierr,xder=2,yder=0,zder=0)  
         gamma(:,1,1) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
         if(nyg>1) then
            call evalX(curv(1),curv(2),curv(3),x,y,z
     .           ,ierr,xder=0,yder=2,zder=0)  
            gamma(:,2,2) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
         else
            gamma(:,2,2) = 0d0
         end if

         if(nzg>1) then
            call evalX(curv(1),curv(2),curv(3),x,y,z
     .           ,ierr,xder=0,yder=0,zder=2)  
            gamma(:,3,3) = x*cov(:,1) + y*cov(:,2) + z*cov(:,3)
         else
            gamma(:,3,3) = 0d0  
         end if
      end if

      cov = cov/Jac             !Luis's convention

!      cov(3,1) = 0d0 ! ADAM: Commented out.
!      cov(3,2) = 0d0
!      cov(3,3) = 1d0

      end subroutine transf_vector_spl

      end module pic_map

c$$$c pic_moments
c$$$c ######################################################################    
c$$$      module pic_moments
c$$$      
c$$$      use pic_io
c$$$
c$$$      use pic_BCs
c$$$
c$$$      use pic_shape_f
c$$$
c$$$      use pic_map
c$$$
c$$$      use grid, ONLY: find_local,setbc,setmgbc_stg
c$$$
c$$$      use variables, ONLY: gv
c$$$
c$$$      !Moments
c$$$      integer,parameter:: Ins=0,Ipx=1,Ipy=2,Ipz=3,IS11=4,IS22=5,IS33=6
c$$$      integer,private :: order_moms=1
c$$$
c$$$      integer,pointer,dimension(:,:) :: mom_bcs
c$$$      real(8),pointer,dimension(:,:,:,:,:) :: moms,momsn,moms_g
c$$$
c$$$      integer,private :: n_sp,nxl,nyl,nzl,nnl,nxg,nyg,nzg,nng
c$$$
c$$$      real(8) :: idlvol
c$$$
c$$$      contains
c$$$
c$$$c     alloc_moms
c$$$c     ##################################################################
c$$$      subroutine alloc_moms(nxd,nyd,nzd,nx,ny,nz,nsp)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: nxd,nyd,nzd,nx,ny,nz,nsp
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$      n_sp = nsp  !Particles not allocated yet
c$$$
c$$$      nxl = nx; nxg = nxd;
c$$$      nyl = ny; nyg = nyd;
c$$$      nzl = nz; nzg = nzd;
c$$$
c$$$      allocate(moms  (0:nx +1,0:ny +1,0:nz +1,n_sp,0:6)) !0=>n,1=>px,2=>py,3=>pz,4=>vx*vx,5=>vy*vy,6=>vz*vz
c$$$      allocate(momsn (0:nx +1,0:ny +1,0:nz +1,n_sp,0:6))
c$$$      allocate(moms_g(0:nxd+1,0:nyd+1,0:nzd+1,n_sp,0:6))
c$$$
c$$$      allocate(mom_bcs(6,0:6))
c$$$
c$$$      moms = 0d0
c$$$      momsn= 0d0
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine alloc_moms
c$$$
c$$$c     dealloc_moms
c$$$c     ##################################################################
c$$$      subroutine dealloc_moms
c$$$
c$$$      deallocate(moms,momsn,mom_bcs,moms_g)
c$$$
c$$$      end subroutine dealloc_moms
c$$$
c$$$c     defineBCs_moms
c$$$c     ###############################################################
c$$$      subroutine defineBCs_moms
c$$$c     ---------------------------------------------------------------
c$$$c     Defines boundary conditions of physical quantities.
c$$$c     On input:
c$$$c       * varray
c$$$c     On output:
c$$$c       * bbcs -> real array of size (6,neq) containing BC setup:
c$$$c           + bbcs(1) ---> at x0
c$$$c           + bbcs(2) ---> at x1
c$$$c           + bbcs(3) ---> at y0
c$$$c           + bbcs(4) ---> at y1
c$$$c           + bbcs(5) ---> at z0
c$$$c           + bbcs(6) ---> at z1
c$$$c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c$$$c     a negative BC identifier means that BCs are to be imposed on
c$$$c     covariant components instead of on (default) contravariant comps.
c$$$c     ---------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$      order_moms = 1
c$$$
c$$$      mom_bcs = DEF
c$$$
c$$$      !for cell-centered moments
c$$$      mom_bcs(:,Ins) = pbcond
c$$$      where (mom_bcs(:,Ins)==pOPN) mom_bcs(:,Ins ) = NEU !Cell centered
c$$$      where (mom_bcs(:,Ins)==pREF) mom_bcs(:,Ins ) = NEU !Cell centered
c$$$
c$$$      mom_bcs(:,IS11) = pbcond
c$$$      where (mom_bcs(:,IS11)==pOPN
c$$$     .   .or.mom_bcs(:,IS11)==pREF) mom_bcs(:,IS11) = NEU !Cell centered
c$$$
c$$$      mom_bcs(:,IS22) = pbcond
c$$$      where (mom_bcs(:,IS22)==pOPN
c$$$     .   .or.mom_bcs(:,IS22)==pREF) mom_bcs(:,IS22) = NEU !Cell centered
c$$$
c$$$      mom_bcs(:,IS33) = pbcond
c$$$      where (mom_bcs(:,IS33)==pOPN
c$$$     .   .or.mom_bcs(:,IS33)==pREF) mom_bcs(:,IS33) = NEU !Cell centered
c$$$
c$$$      !for staggered moments
c$$$      mom_bcs(:,Ipx) = pbcond
c$$$      where (mom_bcs(:,Ipx) == pOPN) mom_bcs(:,Ipx) = DIR !All
c$$$
c$$$c$$$      where (mom_bcs(3:6,Ipx) == pREF) mom_bcs(3:6,Ipx) = NEU !Y,Z
c$$$c$$$      where (mom_bcs(1:2,Ipx) == pREF) mom_bcs(1:2,Ipx) = DIR !X
c$$$      where (mom_bcs(:,Ipx) == pREF) mom_bcs(:,Ipx) = DIR !All
c$$$
c$$$      mom_bcs(:,Ipy) = pbcond
c$$$      where (mom_bcs(:,Ipy) == pOPN) mom_bcs(:,Ipy) = DIR !All
c$$$
c$$$c$$$      where (mom_bcs(1:2,Ipy) == pREF) mom_bcs(1:2,Ipy) = NEU !X
c$$$c$$$      where (mom_bcs(3:4,Ipy) == pREF) mom_bcs(3:4,Ipy) = DIR !Y
c$$$c$$$      where (mom_bcs(5:6,Ipy) == pREF) mom_bcs(5:6,Ipy) = NEU !Z
c$$$      where (mom_bcs(:,Ipy) == pREF) mom_bcs(:,Ipy) = DIR !All
c$$$
c$$$      mom_bcs(:,Ipz) = pbcond
c$$$      where (mom_bcs(:,Ipz) == pOPN) mom_bcs(:,Ipz) = DIR !All
c$$$
c$$$c$$$      where (mom_bcs(1:4,Ipz) == pREF) mom_bcs(1:4,Ipz) = NEU !X,Y
c$$$c$$$      where (mom_bcs(5:6,Ipz) == pREF) mom_bcs(5:6,Ipz) = DIR !Z
c$$$      where (mom_bcs(:,Ipz) == pREF) mom_bcs(:,Ipz) = DIR !All
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine defineBCs_moms
c$$$
c$$$c     gather_moments
c$$$c     ##################################################################
c$$$      subroutine gather_moments
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: im
c$$$
c$$$c     Begin program
c$$$
c$$$
c$$$!      call gather_2d_moments(moms_g(:,:,1,1:n_sp,:))
c$$$!      moms_g(:,:,0,:,:) = moms_g(:,:,1,:,:)
c$$$!      moms_g(:,:,2,:,:) = moms_g(:,:,1,:,:)
c$$$
c$$$      call gather_3d_moments(moms_g(:,:,:,:,:))
c$$$
c$$$      do im=0,size(moms,5)-1
c$$$       call find_local(gv%gparams,1,moms_g(:,:,:,:,im),moms(:,:,:,:,im))
c$$$      enddo
c$$$
c$$$      call applyBC_moms
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine gather_moments
c$$$
c$$$c     gather_2d_moments
c$$$c     ##################################################################
c$$$      subroutine gather_2d_moments(mments)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8), intent(out) :: mments(0:,0:,:,0:)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: isp,ip,ii,jj,i,j,numspcs,numpcles
c$$$      integer :: i_np(_Npg), j_np(_Npg), k_np(_Npg)
c$$$
c$$$      real(8) :: garr_n  (0:nxg+1,0:nyg+1)
c$$$     $          ,garr_vx (0:nxg+1,0:nyg+1)
c$$$     $          ,garr_vy (0:nxg+1,0:nyg+1)
c$$$     $          ,garr_vz (0:nxg+1,0:nyg+1)
c$$$     $          ,garr_vx2(0:nxg+1,0:nyg+1)
c$$$     $          ,garr_vy2(0:nxg+1,0:nyg+1)
c$$$     $          ,garr_vz2(0:nxg+1,0:nyg+1)
c$$$
c$$$      real(8) :: My2(1,3),Mx2(3,1),My1(1,2)
c$$$      real(8) :: Mx1(2,1),Mxy2(3,3),My1_int(1,2)
c$$$      real(_prec) :: xp(_Npg),yp(_Npg),v2,w_pcle(_Npg)
c$$$      real(8) :: vcurv(_Npg,3)
c$$$      real(8) :: ihx,ihy,ihz
c$$$
c$$$c     Begin program
c$$$
c$$$      call get_inverse_h(ihx,ihy,ihz)
c$$$
c$$$      idlvol = ihx*ihy*ihz
c$$$
c$$$      numspcs = size(spcs)
c$$$
c$$$      do isp=1,numspcs
c$$$        numpcles = size(spcs(isp)%pcles)
c$$$
c$$$        garr_n = 0d0
c$$$        garr_vx = 0d0
c$$$        garr_vy = 0d0
c$$$        garr_vz = 0d0
c$$$
c$$$        garr_vx2= 0d0
c$$$        garr_vy2= 0d0
c$$$        garr_vz2= 0d0
c$$$        
c$$$!     note that there is some issue with valgrind and libgomp.
c$$$!$omp parallel DEFAULT(SHARED) private(ip,xp,yp,i_np,j_np,vcurv,i
c$$$!$omp.                                ,Mx1,Mx2,My1,My2,Mxy2,v2,w_pcle)
c$$$!$      call set_omp_thread_id()
c$$$!$omp DO 
c$$$!$omp.reduction(+:garr_n,garr_vx,garr_vy,garr_vz
c$$$!$omp.           ,garr_vx2,garr_vy2,garr_vz2
c$$$!$omp.         ) 
c$$$        do ip=nip+1,numpcles
c$$$
c$$$          xp = spcs(isp)%pcles(ip)%x_n(:,1)*ihx
c$$$          yp = spcs(isp)%pcles(ip)%x_n(:,2)*ihy
c$$$
c$$$          w_pcle = spcs(isp)%pcles(ip)%w_n
c$$$
c$$$          call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
c$$$     .                                          ,i_np,j_np,k_np)   
c$$$          
c$$$          if (no_map) then
c$$$             vcurv = spcs(isp)%pcles(ip)%v_np
c$$$          else
c$$$             vcurv = find_vcurv(i_np,j_np,k_np,spcs(isp)%pcles(ip)%x_n
c$$$     .                                        ,spcs(isp)%pcles(ip)%v_n)
c$$$          endif
c$$$
c$$$          do i=1,_Npg
c$$$
c$$$             Mx1(:,1) = shape1_offset_scalar(xp(i))
c$$$             My1(1,:) = shape1_offset_scalar(yp(i))
c$$$             Mx2(:,1) = shape2_offset_scalar(xp(i))
c$$$             My2(1,:) = shape2_offset_scalar(yp(i))
c$$$             
c$$$#ifdef __GNU
c$$$             Mxy2(:,:)= opt_vecvec_3(Mx2,My2)
c$$$#else
c$$$             Mxy2(:,:)= matmul(Mx2,My2)
c$$$#endif
c$$$
c$$$             !accumulate n (quadratic in x,y)
c$$$             garr_n(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) = 
c$$$     .       garr_n(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1)
c$$$     .       + Mxy2*w_pcle(i)
c$$$
c$$$#ifdef __GNU
c$$$             !accumulate vx (linear in x, 2nd in y)
c$$$             garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) = 
c$$$     .       garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) + 
c$$$     .         vcurv(i,1)*opt_vecvec_23(Mx1,My2)*w_pcle(i)
c$$$
c$$$             !accumulate vy (2nd in x, linear in y)
c$$$             garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) = 
c$$$     .       garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) + 
c$$$     .         vcurv(i,2)*opt_vecvec_32(Mx2,My1)*w_pcle(i)
c$$$#else
c$$$             !accumulate vx (linear in x, 2nd in y)
c$$$             garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) = 
c$$$     .       garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) + 
c$$$     .         vcurv(i,1)*matmul(Mx1,My2)*w_pcle(i)
c$$$
c$$$             !accumulate vy (2nd in x, linear in y)
c$$$             garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) = 
c$$$     .       garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) + 
c$$$     .         vcurv(i,2)*matmul(Mx2,My1)*w_pcle(i)
c$$$#endif
c$$$
c$$$             !accumulate vz (2nd in x,y)
c$$$             garr_vz(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) =
c$$$     .       garr_vz(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) + 
c$$$     .         vcurv(i,3)*Mxy2*w_pcle(i)
c$$$          
c$$$             !accumulate v^2 (2nd in x,y, cell-centered)
c$$$             v2 = spcs(isp)%pcles(ip)%v_n(i,1)
c$$$     .           *spcs(isp)%pcles(ip)%v_n(i,1)
c$$$             garr_vx2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) = 
c$$$     .       garr_vx2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) + 
c$$$     .         v2*(Mxy2*w_pcle(i))
c$$$
c$$$             v2 = spcs(isp)%pcles(ip)%v_n(i,2)
c$$$     .           *spcs(isp)%pcles(ip)%v_n(i,2)
c$$$             garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) = 
c$$$     .       garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) + 
c$$$     .         v2*(Mxy2*w_pcle(i))     
c$$$
c$$$             v2 = spcs(isp)%pcles(ip)%v_n(i,3)
c$$$     .           *spcs(isp)%pcles(ip)%v_n(i,3)
c$$$             garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) = 
c$$$     .       garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) + 
c$$$     .         v2*(Mxy2*w_pcle(i))    
c$$$          
c$$$          end do
c$$$
c$$$        enddo
c$$$!$omp end do 
c$$$!$omp end parallel 
c$$$
c$$$!       just the moments without b.c., the b.c. should be done where needed.
c$$$        mments(:,:,isp,IS11) = garr_vx2
c$$$        mments(:,:,isp,IS22) = garr_vy2
c$$$        mments(:,:,isp,IS33) = garr_vz2
c$$$        mments(:,:,isp,Ipx ) = garr_vx
c$$$        mments(:,:,isp,Ipy ) = garr_vy
c$$$        mments(:,:,isp,Ipz ) = garr_vz
c$$$        mments(:,:,isp,Ins ) = garr_n
c$$$
c$$$      enddo
c$$$
c$$$c     Wrap moments around boundaries
c$$$
c$$$      call wrap_moments_nodes(mments,wrap_currents=.true.)
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine gather_2d_moments
c$$$
c$$$c     gather_3d_moments
c$$$c     ##################################################################
c$$$      subroutine gather_3d_moments(mments)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8), intent(out) :: mments(0:,0:,0:,:,0:)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: isp,ip,ii,jj,i,j,numspcs,numpcles
c$$$      integer :: i_np(_Npg), j_np(_Npg), k_np(_Npg)
c$$$
c$$$      real(8) :: garr_n  (0:nxg+1,0:nyg+1,0:nzg+1)
c$$$     $          ,garr_vx (0:nxg+1,0:nyg+1,0:nzg+1)
c$$$     $          ,garr_vy (0:nxg+1,0:nyg+1,0:nzg+1)
c$$$     $          ,garr_vz (0:nxg+1,0:nyg+1,0:nzg+1)
c$$$     $          ,garr_vx2(0:nxg+1,0:nyg+1,0:nzg+1)
c$$$     $          ,garr_vy2(0:nxg+1,0:nyg+1,0:nzg+1)
c$$$     $          ,garr_vz2(0:nxg+1,0:nyg+1,0:nzg+1)
c$$$
c$$$      real(_prec) :: Mx1(2,1,1),My1(1,2,1),Mz1(1,1,2)
c$$$      real(_prec) :: Mx2(3,1,1),My2(1,3,1),Mz2(1,1,3)
c$$$      real(_prec) :: Myz2(3,3),Mzx2(3,3),Mxy2(3,3)
c$$$      real(_prec) :: M_w_prod(3,3,3)
c$$$      real(_prec) :: dumA(3,1), dumB(1,3), dumC(1,3)
c$$$
c$$$
c$$$!      real(8) :: My2(1,3),Mx2(3,1),My1(1,2)
c$$$!      real(8) :: Mx1(2,1),Mxy2(3,3),My1_int(1,2)
c$$$      real(_prec) :: xp(_Npg),yp(_Npg),zp(_Npg),v2,w_pcle(_Npg)
c$$$      real(8) :: vcurv(_Npg,3)
c$$$      real(8) :: ihx,ihy,ihz
c$$$
c$$$c     Begin program
c$$$
c$$$      call get_inverse_h(ihx,ihy,ihz)
c$$$
c$$$      idlvol = ihx*ihy*ihz
c$$$
c$$$      numspcs = size(spcs)
c$$$
c$$$      do isp=1,numspcs
c$$$        numpcles = size(spcs(isp)%pcles)
c$$$
c$$$        garr_n = 0d0
c$$$        garr_vx = 0d0
c$$$        garr_vy = 0d0
c$$$        garr_vz = 0d0
c$$$
c$$$        garr_vx2= 0d0
c$$$        garr_vy2= 0d0
c$$$        garr_vz2= 0d0
c$$$        
c$$$!     note that there is some issue with valgrind and libgomp.
c$$$!$omp parallel DEFAULT(SHARED) private(ip,xp,yp,zp,i_np,j_np,vcurv,i
c$$$!$omp.                                ,Mx1,Mx2,My1,My2,Mz1,Mz2,Myz2
c$$$!$omp.                                ,Mzx2,Mxy2,M_w_prod,v2,w_pcle)
c$$$!$      call set_omp_thread_id()
c$$$!$omp DO 
c$$$!$omp.reduction(+:garr_n,garr_vx,garr_vy,garr_vz
c$$$!$omp.           ,garr_vx2,garr_vy2,garr_vz2
c$$$!$omp.         ) 
c$$$        do ip=nip+1,numpcles
c$$$
c$$$          xp = spcs(isp)%pcles(ip)%x_n(:,1)*ihx
c$$$          yp = spcs(isp)%pcles(ip)%x_n(:,2)*ihy
c$$$          zp = spcs(isp)%pcles(ip)%x_n(:,3)*ihz
c$$$
c$$$          w_pcle = spcs(isp)%pcles(ip)%w_n
c$$$
c$$$          call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
c$$$     .                                          ,i_np,j_np,k_np)   
c$$$          
c$$$          if (no_map) then
c$$$             vcurv = spcs(isp)%pcles(ip)%v_np
c$$$          else
c$$$             vcurv = find_vcurv(i_np,j_np,k_np,spcs(isp)%pcles(ip)%x_n
c$$$     .                                   ,spcs(isp)%pcles(ip)%v_n)
c$$$          endif
c$$$
c$$$          do i=1,_Npg
c$$$
c$$$             Mx1(:,1,1) = shape1_offset_scalar(xp(i))
c$$$             My1(1,:,1) = shape1_offset_scalar(yp(i))
c$$$             Mz1(1,1,:) = shape1_offset_scalar(zp(i))
c$$$
c$$$
c$$$             Mx2(:,1,1) = shape2_offset_scalar(xp(i))
c$$$             My2(1,:,1) = shape2_offset_scalar(yp(i))
c$$$             Mz2(1,1,:) = shape2_offset_scalar(zp(i))
c$$$
c$$$ccc!#ifdef __GNU
c$$$ccc             Mxy2(:,:)= opt_vecvec_3(Mx2,My2)
c$$$ccc!#else
c$$$ccc             Mxy2(:,:)= matmul(Mx2,My2)
c$$$ccc!#endif
c$$$             
c$$$             Mxy2(:,:) = opt_vecvec_3(Mx2(:,:,1),My2(:,:,1))
c$$$
c$$$             dumA(:,1) = My2(1,:,1)
c$$$             dumB(1,:) = Mz2(1,1,:)
c$$$             Myz2(:,:) = opt_vecvec_3(dumA,dumB)
c$$$
c$$$             dumA(:,1) = Mz2(1,1,:)
c$$$             dumB(1,:) = Mx2(:,1,1)
c$$$             Mzx2(:,:) = opt_vecvec_3(dumA,dumB)
c$$$
c$$$             M_w_prod(:,:,:) = opt_vecvecvec_3(Mx2,My2,Mz2)
c$$$
c$$$
c$$$
c$$$             !accumulate n (quadratic in x,y)
c$$$             garr_n(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .             ,k_np(i)-1:k_np(i)+1) = 
c$$$     .       garr_n(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .             ,k_np(i)-1:k_np(i)+1)
c$$$     .       + M_w_prod*w_pcle(i)
c$$$
c$$$cc#ifdef __GNU
c$$$cc             !accumulate vx (linear in x, 2nd in y)
c$$$cc             garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) = 
c$$$cc     .       garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) + 
c$$$cc     .         vcurv(i,1)*opt_vecvec_23(Mx1,My2)*w_pcle(i)
c$$$cc
c$$$cc             !accumulate vy (2nd in x, linear in y)
c$$$cc             garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) = 
c$$$cc     .       garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) + 
c$$$cc     .         vcurv(i,2)*opt_vecvec_32(Mx2,My1)*w_pcle(i)
c$$$cc#else
c$$$cc             !accumulate vx (linear in x, 2nd in y)
c$$$cc             garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) = 
c$$$cc     .       garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1) + 
c$$$cc     .         vcurv(i,1)*matmul(Mx1,My2)*w_pcle(i)
c$$$cc
c$$$cc             !accumulate vy (2nd in x, linear in y)
c$$$cc             garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) = 
c$$$cc     .       garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)) + 
c$$$cc     .         vcurv(i,2)*matmul(Mx2,My1)*w_pcle(i)
c$$$cc#endif
c$$$
c$$$             !accumulate vx (linear in x, 2nd in y, 2nd in z)
c$$$             garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1
c$$$     .              ,k_np(i)-1:k_np(i)+1) =
c$$$     .       garr_vx(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1
c$$$     .              ,k_np(i)-1:k_np(i)+1) +
c$$$     .            vcurv(i,1)*opt_omatvec_233(Myz2,Mx1)*w_pcle(i)
c$$$
c$$$             !accumulate vy (2nd in x, linear in y, 2nd in z)
c$$$             garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)
c$$$     .            ,k_np(i)-1:k_np(i)+1) =                                                                            
c$$$     .       garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)
c$$$     .            ,k_np(i)-1:k_np(i)+1) +                                                                            
c$$$     .            vcurv(i,2)*opt_omatvec_323(Mzx2,My1)*w_pcle(i)
c$$$
c$$$             !accumulate vz (2nd in x, 2nd in y, linear in z)
c$$$             garr_vz(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)) =
c$$$     .       garr_vz(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)) + 
c$$$     .            vcurv(i,3)*opt_omatvec_332(Mxy2,Mz1)*w_pcle(i)
c$$$          
c$$$             !accumulate v^2 (2nd in x,y,z cell-centered)
c$$$             v2 = spcs(isp)%pcles(ip)%v_n(i,1)
c$$$     .           *spcs(isp)%pcles(ip)%v_n(i,1)
c$$$             garr_vx2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)+1) = 
c$$$     .       garr_vx2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)+1) + v2*(M_w_prod*w_pcle(i))
c$$$
c$$$             v2 = spcs(isp)%pcles(ip)%v_n(i,2)
c$$$     .           *spcs(isp)%pcles(ip)%v_n(i,2)
c$$$             garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)+1) =
c$$$     .       garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)+1) + v2*(M_w_prod*w_pcle(i))
c$$$
c$$$             v2 = spcs(isp)%pcles(ip)%v_n(i,3)
c$$$     .           *spcs(isp)%pcles(ip)%v_n(i,3)
c$$$             garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)+1) =
c$$$     .       garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
c$$$     .            ,k_np(i)-1:k_np(i)+1) + v2*(M_w_prod*w_pcle(i))
c$$$
c$$$
c$$$c             v2 = spcs(isp)%pcles(ip)%v_n(i,2)
c$$$c     .           *spcs(isp)%pcles(ip)%v_n(i,2)
c$$$c             garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) = 
c$$$c     .       garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) + 
c$$$c     .         v2*(Mxy2*w_pcle(i))     
c$$$
c$$$c             v2 = spcs(isp)%pcles(ip)%v_n(i,3)
c$$$c     .           *spcs(isp)%pcles(ip)%v_n(i,3)
c$$$c             garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) = 
c$$$c     .       garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1) + 
c$$$c     .         v2*(Mxy2*w_pcle(i))    
c$$$          
c$$$          end do
c$$$
c$$$        enddo
c$$$!$omp end do 
c$$$!$omp end parallel 
c$$$
c$$$!       just the moments without b.c., the b.c. should be done where needed.
c$$$        mments(:,:,:,isp,IS11) = garr_vx2
c$$$        mments(:,:,:,isp,IS22) = garr_vy2
c$$$        mments(:,:,:,isp,IS33) = garr_vz2
c$$$        mments(:,:,:,isp,Ipx ) = garr_vx
c$$$        mments(:,:,:,isp,Ipy ) = garr_vy
c$$$        mments(:,:,:,isp,Ipz ) = garr_vz
c$$$        mments(:,:,:,isp,Ins ) = garr_n
c$$$
c$$$      enddo
c$$$
c$$$c     Wrap moments around boundaries
c$$$
c$$$      call wrap_moments_nodes3d(mments,wrap_currents=.true.)
c$$$
c$$$c     End program
c$$$
c$$$      end subroutine gather_3d_moments
c$$$
c$$$
c$$$c     wrap_moments_nodes
c$$$c     ##################################################################
c$$$      subroutine wrap_moments_nodes(mments,wrap_currents)
c$$$      
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8), intent(inout) :: mments(0:,0:,:,0:)
c$$$      logical, optional :: wrap_currents
c$$$
c$$$c     local variables
c$$$
c$$$      real(8),allocatable,dimension(:) :: send_buf,rec_buf
c$$$      logical :: wcurr
c$$$
c$$$c     Begin program
c$$$
c$$$      if (PRESENT(wrap_currents)) then
c$$$        wcurr = wrap_currents
c$$$      else
c$$$        wcurr = .true.
c$$$      endif
c$$$
c$$$      if (pbcond(1) == pPER) then
c$$$         mments(1  ,:,:,:) = mments(1,:,:,:)+mments(nxg+1,:,:,:)
c$$$         mments(nxg,:,:,:) = mments(0,:,:,:)+mments(nxg  ,:,:,:)
c$$$         mments(nxg+1,:,:,:) = mments(1  ,:,:,:)
c$$$         mments(0    ,:,:,:) = mments(nxg,:,:,:)
c$$$      else
c$$$        if (wcurr) then
c$$$          !Jx
c$$$          if (     pbcond(1) == pREF
c$$$     .        .or. pbcond(1) == pOPN) mments(0  ,:,:,Ipx) = 0d0
c$$$          if (     pbcond(2) == pREF
c$$$     .        .or. pbcond(2) == pOPN) mments(nxg,:,:,Ipx) = 0d0
c$$$
c$$$          !Recover left-over deposition for second-order splines at ghost cell
c$$$          ! for staggered quantities
c$$$          !Jy
c$$$          if (pbcond(1) == pOPN)
c$$$     .      mments(1  ,:,:,Ipy) = mments(1    ,:,:,Ipy)
c$$$     .                           -mments(0    ,:,:,Ipy)
c$$$          if (pbcond(1) == pREF)
c$$$     .      mments(1  ,:,:,Ipy) = mments(1    ,:,:,Ipy)
c$$$     .                           +mments(0    ,:,:,Ipy)
c$$$          if (pbcond(2) == pOPN)
c$$$     .      mments(nxg,:,:,Ipy) = mments(nxg  ,:,:,Ipy)
c$$$     .                           -mments(nxg+1,:,:,Ipy)
c$$$          if (pbcond(2) == pREF)
c$$$     .      mments(nxg,:,:,Ipy) = mments(nxg  ,:,:,Ipy)
c$$$     .                           +mments(nxg+1,:,:,Ipy)
c$$$
c$$$          !Jz
c$$$          if (pbcond(1) == pOPN)
c$$$     .      mments(1  ,:,:,Ipz) = mments(1    ,:,:,Ipz)
c$$$     .                           -mments(0    ,:,:,Ipz)
c$$$          if (pbcond(1) == pREF)
c$$$     .      mments(1  ,:,:,Ipz) = mments(1    ,:,:,Ipz)
c$$$     .                           +mments(0    ,:,:,Ipz)
c$$$
c$$$          if (pbcond(2) == pOPN)
c$$$     .      mments(nxg,:,:,Ipz) = mments(nxg  ,:,:,Ipz)
c$$$     .                           -mments(nxg+1,:,:,Ipz)
c$$$          if (pbcond(2) == pREF)
c$$$     .      mments(nxg,:,:,Ipz) = mments(nxg  ,:,:,Ipz)
c$$$     .                           +mments(nxg+1,:,:,Ipz)
c$$$
c$$$          !Rho
c$$$          mments(1  ,:,:,Ins) = mments(1    ,:,:,Ins)
c$$$     .                         +mments(0    ,:,:,Ins)
c$$$          mments(nxg,:,:,Ins) = mments(nxg  ,:,:,Ins)
c$$$     .                         +mments(nxg+1,:,:,Ins)
c$$$        endif
c$$$
c$$$        !S11
c$$$        mments(1  ,:,:,IS11) = mments(1    ,:,:,IS11)
c$$$     .                        +mments(0    ,:,:,IS11)
c$$$        mments(nxg,:,:,IS11) = mments(nxg  ,:,:,IS11)
c$$$     .                        +mments(nxg+1,:,:,IS11)
c$$$        !S22
c$$$        mments(1  ,:,:,IS22) = mments(1    ,:,:,IS22)
c$$$     .                        +mments(0    ,:,:,IS22)
c$$$        mments(nxg,:,:,IS22) = mments(nxg  ,:,:,IS22)
c$$$     .                        +mments(nxg+1,:,:,IS22)
c$$$        !S33
c$$$        mments(1  ,:,:,IS33) = mments(1    ,:,:,IS33)
c$$$     .                        +mments(0    ,:,:,IS33)
c$$$        mments(nxg,:,:,IS33) = mments(nxg  ,:,:,IS33)
c$$$     .                        +mments(nxg+1,:,:,IS33)
c$$$
c$$$      endif
c$$$
c$$$      if (pbcond(3) == pPER) then
c$$$         mments(:  ,1,:,:) = mments(:,1,:,:)+mments(:,nyg+1,:,:)
c$$$         mments(:,nyg,:,:) = mments(:,0,:,:)+mments(:,nyg  ,:,:)
c$$$         mments(:,nyg+1,:,:) = mments(:,1  ,:,:)
c$$$         mments(:,0    ,:,:) = mments(:,nyg,:,:)
c$$$      else
c$$$        if (wcurr) then
c$$$          !Jy
c$$$          if (     pbcond(3) == pREF
c$$$     .        .or. pbcond(3) == pOPN) mments(:,0  ,:,Ipy) = 0d0
c$$$          if (     pbcond(4) == pREF               
c$$$     .        .or. pbcond(4) == pOPN) mments(:,nyg,:,Ipy) = 0d0
c$$$
c$$$          !Recover left-over deposition for second-order splines at ghost cell
c$$$          ! for staggered quantities
c$$$          !Jx
c$$$          if (pbcond(3) == pOPN)
c$$$     .      mments(:,1  ,:,Ipx) = mments(:,1    ,:,Ipx)
c$$$     .                           -mments(:,0    ,:,Ipx)
c$$$          if (pbcond(3) == pREF)
c$$$     .      mments(:,1  ,:,Ipx) = mments(:,1    ,:,Ipx)
c$$$     .                           +mments(:,0    ,:,Ipx)
c$$$
c$$$          if (pbcond(4) == pOPN)
c$$$     .      mments(:,nyg,:,Ipx) = mments(:,nyg  ,:,Ipx)
c$$$     .                           -mments(:,nyg+1,:,Ipx)
c$$$          if (pbcond(4) == pREF)
c$$$     .      mments(:,nyg,:,Ipx) = mments(:,nyg  ,:,Ipx)
c$$$     .                           +mments(:,nyg+1,:,Ipx)
c$$$
c$$$          !Jz
c$$$          if (pbcond(3) == pOPN)
c$$$     .      mments(:,1  ,:,Ipz) = mments(:,1    ,:,Ipz)
c$$$     .                           -mments(:,0    ,:,Ipz)
c$$$          if (pbcond(3) == pREF)
c$$$     .      mments(:,1  ,:,Ipz) = mments(:,1    ,:,Ipz)
c$$$     .                           +mments(:,0    ,:,Ipz)
c$$$
c$$$          if (pbcond(4) == pOPN)
c$$$     .      mments(:,nyg,:,Ipz) = mments(:,nyg  ,:,Ipz)
c$$$     .                           -mments(:,nyg+1,:,Ipz)
c$$$          if (pbcond(4) == pREF)
c$$$     .      mments(:,nyg,:,Ipz) = mments(:,nyg  ,:,Ipz)
c$$$     .                           +mments(:,nyg+1,:,Ipz)
c$$$
c$$$          !Rho
c$$$          mments(:,1  ,:,Ins) = mments(:,1    ,:,Ins)
c$$$     .                         +mments(:,0    ,:,Ins)
c$$$          mments(:,nyg,:,Ins) = mments(:,nyg  ,:,Ins)
c$$$     .                         +mments(:,nyg+1,:,Ins)
c$$$#if defined(PIC_CM_DIAG)
c$$$#endif
c$$$        endif
c$$$
c$$$        !S11
c$$$        mments(:,1  ,:,IS11) = mments(:,1    ,:,IS11)
c$$$     .                        +mments(:,0    ,:,IS11)
c$$$        mments(:,nyg,:,IS11) = mments(:,nyg  ,:,IS11)
c$$$     .                        +mments(:,nyg+1,:,IS11)
c$$$        !S22                            
c$$$        mments(:,1  ,:,IS22) = mments(:,1    ,:,IS22)
c$$$     .                        +mments(:,0    ,:,IS22)
c$$$        mments(:,nyg,:,IS22) = mments(:,nyg  ,:,IS22)
c$$$     .                        +mments(:,nyg+1,:,IS22)
c$$$        !S33                            
c$$$        mments(:,1  ,:,IS33) = mments(:,1    ,:,IS33)
c$$$     .                        +mments(:,0    ,:,IS33)
c$$$        mments(:,nyg,:,IS33) = mments(:,nyg  ,:,IS33)
c$$$     .                        +mments(:,nyg+1,:,IS33)
c$$$
c$$$      endif
c$$$
c$$$#if defined(petsc)
c$$$      allocate(send_buf(size(mments))
c$$$     $        ,rec_buf (size(mments)))
c$$$      send_buf = reshape(mments, (/size(mments)/))
c$$$      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
c$$$     $                  ,MPI_DOUBLE_PRECISION
c$$$     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
c$$$      mments(0:,0:,:,0:) = reshape(rec_buf,shape(mments))
c$$$      deallocate(send_buf,rec_buf)
c$$$#endif
c$$$
c$$$c     end program
c$$$
c$$$      end subroutine wrap_moments_nodes
c$$$
c$$$c     wrap_moments_nodes3d
c$$$c     ##################################################################
c$$$      subroutine wrap_moments_nodes3d(mments,wrap_currents)
c$$$      
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8), intent(inout) :: mments(0:,0:,0:,:,0:)
c$$$      logical, optional :: wrap_currents
c$$$
c$$$c     local variables
c$$$
c$$$      real(8),allocatable,dimension(:) :: send_buf,rec_buf
c$$$      logical :: wcurr
c$$$
c$$$c     Begin program
c$$$
c$$$      if (PRESENT(wrap_currents)) then
c$$$        wcurr = wrap_currents
c$$$      else
c$$$        wcurr = .true.
c$$$      endif
c$$$
c$$$      if (pbcond(1) == pPER) then
c$$$         mments(1  ,:,:,:,:) = mments(1,:,:,:,:)+mments(nxg+1,:,:,:,:)
c$$$         mments(nxg,:,:,:,:) = mments(0,:,:,:,:)+mments(nxg  ,:,:,:,:)
c$$$         mments(nxg+1,:,:,:,:) = mments(1  ,:,:,:,:)
c$$$         mments(0    ,:,:,:,:) = mments(nxg,:,:,:,:)
c$$$      else
c$$$        if (wcurr) then
c$$$          !Jx
c$$$          if (     pbcond(1) == pREF
c$$$     .        .or. pbcond(1) == pOPN) mments(0  ,:,:,:,Ipx) = 0d0
c$$$          if (     pbcond(2) == pREF
c$$$     .        .or. pbcond(2) == pOPN) mments(nxg,:,:,:,Ipx) = 0d0
c$$$
c$$$          !Recover left-over deposition for second-order splines at ghost cell
c$$$          ! for staggered quantities
c$$$          !Jy
c$$$          if (pbcond(1) == pOPN)
c$$$     .      mments(1  ,:,:,:,Ipy) = mments(1    ,:,:,:,Ipy)
c$$$     .                           -mments(0    ,:,:,:,Ipy)
c$$$          if (pbcond(1) == pREF)
c$$$     .      mments(1  ,:,:,:,Ipy) = mments(1    ,:,:,:,Ipy)
c$$$     .                           +mments(0    ,:,:,:,Ipy)
c$$$          if (pbcond(2) == pOPN)
c$$$     .      mments(nxg,:,:,:,Ipy) = mments(nxg  ,:,:,:,Ipy)
c$$$     .                           -mments(nxg+1,:,:,:,Ipy)
c$$$          if (pbcond(2) == pREF)
c$$$     .      mments(nxg,:,:,:,Ipy) = mments(nxg  ,:,:,:,Ipy)
c$$$     .                           +mments(nxg+1,:,:,:,Ipy)
c$$$
c$$$          !Jz
c$$$          if (pbcond(1) == pOPN)
c$$$     .      mments(1  ,:,:,:,Ipz) = mments(1    ,:,:,:,Ipz)
c$$$     .                           -mments(0    ,:,:,:,Ipz)
c$$$          if (pbcond(1) == pREF)
c$$$     .      mments(1  ,:,:,:,Ipz) = mments(1    ,:,:,:,Ipz)
c$$$     .                           +mments(0    ,:,:,:,Ipz)
c$$$
c$$$          if (pbcond(2) == pOPN)
c$$$     .      mments(nxg,:,:,:,Ipz) = mments(nxg  ,:,:,:,Ipz)
c$$$     .                           -mments(nxg+1,:,:,:,Ipz)
c$$$          if (pbcond(2) == pREF)
c$$$     .      mments(nxg,:,:,:,Ipz) = mments(nxg  ,:,:,:,Ipz)
c$$$     .                           +mments(nxg+1,:,:,:,Ipz)
c$$$
c$$$          !Rho
c$$$          mments(1  ,:,:,:,Ins) = mments(1    ,:,:,:,Ins)
c$$$     .                         +mments(0    ,:,:,:,Ins)
c$$$          mments(nxg,:,:,:,Ins) = mments(nxg  ,:,:,:,Ins)
c$$$     .                         +mments(nxg+1,:,:,:,Ins)
c$$$        endif
c$$$
c$$$        !S11
c$$$        mments(1  ,:,:,:,IS11) = mments(1    ,:,:,:,IS11)
c$$$     .                        +mments(0    ,:,:,:,IS11)
c$$$        mments(nxg,:,:,:,IS11) = mments(nxg  ,:,:,:,IS11)
c$$$     .                        +mments(nxg+1,:,:,:,IS11)
c$$$        !S22
c$$$        mments(1  ,:,:,:,IS22) = mments(1    ,:,:,:,IS22)
c$$$     .                        +mments(0    ,:,:,:,IS22)
c$$$        mments(nxg,:,:,:,IS22) = mments(nxg  ,:,:,:,IS22)
c$$$     .                        +mments(nxg+1,:,:,:,IS22)
c$$$        !S33
c$$$        mments(1  ,:,:,:,IS33) = mments(1    ,:,:,:,IS33)
c$$$     .                        +mments(0    ,:,:,:,IS33)
c$$$        mments(nxg,:,:,:,IS33) = mments(nxg  ,:,:,:,IS33)
c$$$     .                        +mments(nxg+1,:,:,:,IS33)
c$$$
c$$$      endif
c$$$
c$$$      if (pbcond(3) == pPER) then
c$$$         mments(:  ,1,:,:,:) = mments(:,1,:,:,:)+mments(:,nyg+1,:,:,:)
c$$$         mments(:,nyg,:,:,:) = mments(:,0,:,:,:)+mments(:,nyg  ,:,:,:)
c$$$         mments(:,nyg+1,:,:,:) = mments(:,1  ,:,:,:)
c$$$         mments(:,0    ,:,:,:) = mments(:,nyg,:,:,:)
c$$$      else
c$$$        if (wcurr) then
c$$$          !Jy
c$$$          if (     pbcond(3) == pREF
c$$$     .        .or. pbcond(3) == pOPN) mments(:,0  ,:,:,Ipy) = 0d0
c$$$          if (     pbcond(4) == pREF               
c$$$     .        .or. pbcond(4) == pOPN) mments(:,nyg,:,:,Ipy) = 0d0
c$$$
c$$$          !Recover left-over deposition for second-order splines at ghost cell
c$$$          ! for staggered quantities
c$$$          !Jx
c$$$          if (pbcond(3) == pOPN)
c$$$     .      mments(:,1  ,:,:,Ipx) = mments(:,1    ,:,:,Ipx)
c$$$     .                           -mments(:,0    ,:,:,Ipx)
c$$$          if (pbcond(3) == pREF)
c$$$     .      mments(:,1  ,:,:,Ipx) = mments(:,1    ,:,:,Ipx)
c$$$     .                           +mments(:,0    ,:,:,Ipx)
c$$$
c$$$          if (pbcond(4) == pOPN)
c$$$     .      mments(:,nyg,:,:,Ipx) = mments(:,nyg  ,:,:,Ipx)
c$$$     .                           -mments(:,nyg+1,:,:,Ipx)
c$$$          if (pbcond(4) == pREF)
c$$$     .      mments(:,nyg,:,:,Ipx) = mments(:,nyg  ,:,:,Ipx)
c$$$     .                           +mments(:,nyg+1,:,:,Ipx)
c$$$
c$$$          !Jz
c$$$          if (pbcond(3) == pOPN)
c$$$     .      mments(:,1  ,:,:,Ipz) = mments(:,1    ,:,:,Ipz)
c$$$     .                           -mments(:,0    ,:,:,Ipz)
c$$$          if (pbcond(3) == pREF)
c$$$     .      mments(:,1  ,:,:,Ipz) = mments(:,1    ,:,:,Ipz)
c$$$     .                           +mments(:,0    ,:,:,Ipz)
c$$$
c$$$          if (pbcond(4) == pOPN)
c$$$     .      mments(:,nyg,:,:,Ipz) = mments(:,nyg  ,:,:,Ipz)
c$$$     .                           -mments(:,nyg+1,:,:,Ipz)
c$$$          if (pbcond(4) == pREF)
c$$$     .      mments(:,nyg,:,:,Ipz) = mments(:,nyg  ,:,:,Ipz)
c$$$     .                           +mments(:,nyg+1,:,:,Ipz)
c$$$
c$$$          !Rho
c$$$          mments(:,1  ,:,:,Ins) = mments(:,1    ,:,:,Ins)
c$$$     .                         +mments(:,0    ,:,:,Ins)
c$$$          mments(:,nyg,:,:,Ins) = mments(:,nyg  ,:,:,Ins)
c$$$     .                         +mments(:,nyg+1,:,:,Ins)
c$$$#if defined(PIC_CM_DIAG)
c$$$#endif
c$$$        endif
c$$$
c$$$        !S11
c$$$        mments(:,1  ,:,:,IS11) = mments(:,1    ,:,:,IS11)
c$$$     .                        +mments(:,0    ,:,:,IS11)
c$$$        mments(:,nyg,:,:,IS11) = mments(:,nyg  ,:,:,IS11)
c$$$     .                        +mments(:,nyg+1,:,:,IS11)
c$$$        !S22                            
c$$$        mments(:,1  ,:,:,IS22) = mments(:,1    ,:,:,IS22)
c$$$     .                        +mments(:,0    ,:,:,IS22)
c$$$        mments(:,nyg,:,:,IS22) = mments(:,nyg  ,:,:,IS22)
c$$$     .                        +mments(:,nyg+1,:,:,IS22)
c$$$        !S33                            
c$$$        mments(:,1  ,:,:,IS33) = mments(:,1    ,:,:,IS33)
c$$$     .                        +mments(:,0    ,:,:,IS33)
c$$$        mments(:,nyg,:,:,IS33) = mments(:,nyg  ,:,:,IS33)
c$$$     .                        +mments(:,nyg+1,:,:,IS33)
c$$$
c$$$      endif
c$$$
c$$$      if (pbcond(5) == pPER) then
c$$$         mments(:,:,1  ,:,:) = mments(:,:,1,:,:)+mments(:,:,nzg+1,:,:)
c$$$         mments(:,:,nzg,:,:) = mments(:,:,0,:,:)+mments(:,:,nzg  ,:,:)
c$$$         mments(:,:,nzg+1,:,:) = mments(:,:,1  ,:,:)
c$$$         mments(:,:,0    ,:,:) = mments(:,:,nzg,:,:)
c$$$      else
c$$$         call pstop('readinput','Only per boundary implemented in z')
c$$$      endif
c$$$
c$$$#if defined(petsc)
c$$$      allocate(send_buf(size(mments))
c$$$     $        ,rec_buf (size(mments)))
c$$$      send_buf = reshape(mments, (/size(mments)/))
c$$$      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
c$$$     $                  ,MPI_DOUBLE_PRECISION
c$$$     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
c$$$      mments(0:,0:,0:,:,0:) = reshape(rec_buf,shape(mments))
c$$$      deallocate(send_buf,rec_buf)
c$$$#endif
c$$$
c$$$c     end program
c$$$
c$$$      end subroutine wrap_moments_nodes3d
c$$$
c$$$c     applyBC_moms
c$$$c     ##################################################################
c$$$      subroutine applyBC_moms
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: im,isp,nx,ny,nz,igrid,order,icomp
c$$$      real(8),pointer,dimension(:,:,:) :: dummy
c$$$
c$$$c     Begin program
c$$$
c$$$      igrid = 1
c$$$
c$$$      icomp = 8 !Make sure we do not interfere with dependent vars
c$$$
c$$$      nx = size(moms,1)-2
c$$$      ny = size(moms,2)-2
c$$$      nz = size(moms,3)-2
c$$$
c$$$      allocate(dummy(0:nx+1,0:ny+1,0:nz+1))
c$$$
c$$$      dummy = 0d0
c$$$
c$$$      !Densities
c$$$      do isp = 1,n_sp
c$$$        call setBC(gv%gparams,icomp,nx,ny,nz
c$$$     .          ,moms(:,:,:,isp,Ins)
c$$$     .          ,dummy
c$$$     .          ,mom_bcs(:,Ins)
c$$$     .          ,igrid,igrid,igrid,iorder=order_moms)
c$$$      enddo
c$$$
c$$$      !Momenta
c$$$      do isp = 1,n_sp
c$$$        call setMGBC_stg(gv%gparams,0,3,nx,ny,nz,igrid
c$$$     $              ,moms(:,:,:,isp,IPx:IPz)
c$$$     $              ,mom_bcs(:,IPx:IPz)
c$$$     .              ,is_cnv=.true.,icomp=(/icomp+1/)
c$$$     .              ,iorder=order_moms)
c$$$      enddo
c$$$
c$$$      !Temperatures
c$$$      do im=IS11,IS33
c$$$        do isp = 1,n_sp
c$$$          call setBC(gv%gparams,icomp,nx,ny,nz
c$$$     .              ,moms(:,:,:,isp,im)
c$$$     .              ,dummy
c$$$     .              ,mom_bcs(:,im)
c$$$     .              ,igrid,igrid,igrid,iorder=order_moms)
c$$$        enddo
c$$$      enddo
c$$$
c$$$      deallocate(dummy)
c$$$
c$$$      end subroutine applyBC_moms
c$$$
c$$$c     find_rho_from_pcles
c$$$c     ##################################################################
c$$$      subroutine find_rho_from_pcles(rho)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8) :: rho(0:,0:,0:)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: isp
c$$$      real(8) :: rho0,ihx,ihy,ihz
c$$$
c$$$c     Define current
c$$$
c$$$      call get_inverse_h(ihx,ihy,ihz)
c$$$
c$$$      idlvol = ihx*ihy*ihz
c$$$
c$$$      rho(:,:,:) = 0d0
c$$$      do isp = 1,n_sp
c$$$        rho = rho + spcs(isp)%q*moms(:,:,:,isp,Ins)*idlvol
c$$$      enddo
c$$$
c$$$      end subroutine find_rho_from_pcles
c$$$
c$$$c     find_J_from_pcles
c$$$c     ############################################################
c$$$      subroutine find_J_from_pcles(jcnv)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8) :: jcnv(0:,0:,0:,:)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: isp
c$$$      real(8) :: ihx,ihy,ihz
c$$$
c$$$c     Define current
c$$$
c$$$      call get_inverse_h(ihx,ihy,ihz)
c$$$
c$$$      idlvol = ihx*ihy*ihz
c$$$
c$$$      jcnv = 0d0
c$$$      do isp = 1,n_sp
c$$$        jcnv(:,:,:,1) = jcnv(:,:,:,1)
c$$$     .                + spcs(isp)%q*moms(:,:,:,isp,Ipx)*idlvol
c$$$        jcnv(:,:,:,2) = jcnv(:,:,:,2)
c$$$     .                + spcs(isp)%q*moms(:,:,:,isp,Ipy)*idlvol
c$$$        jcnv(:,:,:,3) = jcnv(:,:,:,3)
c$$$     .                + spcs(isp)%q*moms(:,:,:,isp,Ipz)*idlvol
c$$$      enddo
c$$$
c$$$      end subroutine find_J_from_pcles
c$$$
c$$$      end module pic_moments

c hybrid moments
c ######################################################################
      module hybrid_moments

      use pic_io

      use pic_BCs

      use pic_shape_f

      use pic_map

      use grid, ONLY: find_local,setbc,setmgbc,setmgbc_stg ! stg for DPIC mass consv. diag.

      use variables, ONLY: gv

      !Moments
      integer,parameter:: Ins=0,Ipx=1,Ipy=2,Ipz=3,IS11=4,IS22=5,IS33=6
#if defined(HYB_MASS_DIAG)
      integer,parameter:: InDPIC=7,IpxDPIC=8,IpyDPIC=9,IpzDPIC=10
#endif

      integer,private :: order_moms=1

      integer,pointer,dimension(:,:) :: mom_bcs
      real(8),pointer,dimension(:,:,:,:,:) :: moms,momsn,moms_g

      integer,private :: n_sp,nxl,nyl,nzl,nnl,nxg,nyg,nzg,nng

      real(8) :: idlvol

      contains

c     alloc_moms
c     ##################################################################
      subroutine alloc_moms(nxd,nyd,nzd,nx,ny,nz,nsp)

      implicit none

c     Call variables

      integer :: nxd,nyd,nzd,nx,ny,nz,nsp

c     Local variables

c     Begin program

      n_sp = nsp  !Particles not allocated yet

      nxl = nx; nxg = nxd;
      nyl = ny; nyg = nyd;
      nzl = nz; nzg = nzd;

#if defined(HYB_MASS_DIAG)
      allocate(moms  (0:nx +1,0:ny +1,0:nz +1,n_sp,0:10)) !0=>n,1=>px,2=>py,3=>pz,4=>vx*vx,5=>vy*vy,6=>vz*vz
      allocate(momsn (0:nx +1,0:ny +1,0:nz +1,n_sp,0:10)) !7=>n,8=>px,9=>py,10=>pz using DPIC definitions
      allocate(moms_g(0:nxd+1,0:nyd+1,0:nzd+1,n_sp,0:10))
      allocate(mom_bcs(6,0:10))
#else
      allocate(moms  (0:nx +1,0:ny +1,0:nz +1,n_sp,0:6)) !0=>n,1=>px,2=>py,3=>pz,4=>vx*vx,5=>vy*vy,6=>vz*vz
      allocate(momsn (0:nx +1,0:ny +1,0:nz +1,n_sp,0:6))
      allocate(moms_g(0:nxd+1,0:nyd+1,0:nzd+1,n_sp,0:6))
      allocate(mom_bcs(6,0:6))
#endif

      moms = 0d0
      momsn= 0d0

c     End program

      end subroutine alloc_moms

c     dealloc_moms
c     ##################################################################
      subroutine dealloc_moms

      deallocate(moms,momsn,mom_bcs,moms_g)

      end subroutine dealloc_moms

c     defineBCs_moms
c     ###############################################################
      subroutine defineBCs_moms(bc_rho,bc_vx,bc_vy,bc_vz,bc_T)
c     ---------------------------------------------------------------
c     Defines boundary conditions of physical quantities.
c     On input:
c       * varray
c     On output:
c       * bbcs -> real array of size (6,neq) containing BC setup:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: bc_rho(6),bc_T(6),bc_vx(6),bc_vy(6),bc_vz(6)

c     Local variables

c     Begin program

      order_moms = 1

      mom_bcs = DEF

      mom_bcs(:,Ins) = bc_rho

      mom_bcs(:,IS11)= bc_T
      mom_bcs(:,IS22)= bc_T
      mom_bcs(:,IS33)= bc_T

      mom_bcs(:,Ipx) = bc_vx
      mom_bcs(:,Ipy) = bc_vy
      mom_bcs(:,Ipz) = bc_vz

#if defined(HYB_MASS_DIAG)
      mom_bcs(:,InDPIC)  = bc_rho
      mom_bcs(:,IpxDPIC) = bc_vx
      mom_bcs(:,IpyDPIC) = bc_vy
      mom_bcs(:,IpzDPIC) = bc_vz
#endif

c     End program

      end subroutine defineBCs_moms

c     gather_moments
c     ##################################################################
      subroutine gather_moments

      implicit none

c     Call variables

c     Local variables

      integer :: im

c     Begin program

      call gather_3d_moments(moms_g(:,:,:,:,:))

      do im=0,size(moms,5)-1
       call find_local(gv%gparams,1,moms_g(:,:,:,:,im),moms(:,:,:,:,im))
      enddo

      call applyBC_moms

c     End program

      end subroutine gather_moments

c     gather_3d_moments
c     ##################################################################
      subroutine gather_3d_moments(mments)

      implicit none

c     Call variables

      real(8), intent(out) :: mments(0:,0:,0:,:,0:)

c     Local variables

      integer :: isp,ip,ii,jj,i,j,numspcs,numpcles
      integer :: i_np(_Npg), j_np(_Npg), k_np(_Npg)

      real(8) :: garr_n  (0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vx (0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vy (0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vz (0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vx2(0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vy2(0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vz2(0:nxg+1,0:nyg+1,0:nzg+1)
#if defined(HYB_MASS_DIAG)
     $          ,garr_nDPIC (0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vxDPIC(0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vyDPIC(0:nxg+1,0:nyg+1,0:nzg+1)
     $          ,garr_vzDPIC(0:nxg+1,0:nyg+1,0:nzg+1)

      real(_prec) :: Mx1(2,1,1),My1(1,2,1),Mz1(1,1,2)
      real(_prec) :: Myz2(3,3),Mzx2(3,3),Mxy2(3,3)
      real(_prec) :: dumA(3,1), dumB(1,3), dumC(1,3)

#endif
      real(_prec) :: Mx2(3,1,1),My2(1,3,1),Mz2(1,1,3)
      real(_prec) :: M_w_prod(3,3,3)

      real(_prec) :: xp(_Npg),yp(_Npg),zp(_Npg),v2,w_pcle(_Npg)
      real(8) :: vcurv(_Npg,3)
      real(8) :: ihx,ihy,ihz

c     Begin program

      call get_inverse_h(ihx,ihy,ihz)

      idlvol = ihx*ihy*ihz

      numspcs = size(spcs)

      do isp=1,numspcs
        numpcles = size(spcs(isp)%pcles)

        garr_n = 0d0
        garr_vx = 0d0
        garr_vy = 0d0
        garr_vz = 0d0

        garr_vx2= 0d0
        garr_vy2= 0d0
        garr_vz2= 0d0

#if defined(HYB_MASS_DIAG)
        garr_nDPIC  = 0d0
        garr_vxDPIC = 0d0
        garr_vyDPIC = 0d0
        garr_vzDPIC = 0d0
#endif

!     note that there is some issue with valgrind and libgomp.
!$omp parallel DEFAULT(SHARED) private(ip,xp,yp,zp,i_np,j_np,k_np,vcurv,
#if defined(HYB_MASS_DIAG)
!$omp.                                 Mx1,My1,Mz1,Myz2,Mzx2,Mxy2,
#endif
!$omp.                                 i,Mx2,My2,Mz2,M_w_prod,v2,w_pcle)
!$      call set_omp_thread_id()
!$omp DO 
!$omp.reduction(+:garr_n,garr_vx,garr_vy,garr_vz
!$omp.           ,garr_vx2,garr_vy2,garr_vz2
#if defined(HYB_MASS_DIAG)
!$omp.           ,garr_nDPIC,garr_vxDPIC,garr_vyDPIC,garr_vzDPIC
#endif
!$omp.         ) 
        do ip=nip+1,numpcles

          xp = spcs(isp)%pcles(ip)%x_n(:,1)*ihx
          yp = spcs(isp)%pcles(ip)%x_n(:,2)*ihy
          zp = spcs(isp)%pcles(ip)%x_n(:,3)*ihz

          w_pcle = spcs(isp)%pcles(ip)%w_n

          call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
     .                                          ,i_np,j_np,k_np)
      
          if (no_map) then
             vcurv = spcs(isp)%pcles(ip)%v_np
          else
!             call pstop('gather_3d_moments', 'No curvilinear yet')
             vcurv = find_vcurv(i_np,j_np,k_np,spcs(isp)%pcles(ip)%x_n
     .                                   ,spcs(isp)%pcles(ip)%v_n)
          endif

          do i=1,_Npg
             Mx2(:,1,1) = shape2_offset_scalar(xp(i))
             My2(1,:,1) = shape2_offset_scalar(yp(i))
             Mz2(1,1,:) = shape2_offset_scalar(zp(i))
             M_w_prod(:,:,:) = opt_vecvecvec_3(Mx2,My2,Mz2)

#if defined(HYB_MASS_DIAG)
            Mx1(:,1,1) = shape1_offset_scalar(xp(i))
            My1(1,:,1) = shape1_offset_scalar(yp(i))
            Mz1(1,1,:) = shape1_offset_scalar(zp(i))
            Mxy2(:,:) = opt_vecvec_3(Mx2(:,:,1),My2(:,:,1))

            dumA(:,1) = My2(1,:,1)
            dumB(1,:) = Mz2(1,1,:)
            Myz2(:,:) = opt_vecvec_3(dumA,dumB)

            dumA(:,1) = Mz2(1,1,:)
            dumB(1,:) = Mx2(:,1,1)
            Mzx2(:,:) = opt_vecvec_3(dumA,dumB)
#endif


             garr_n(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) = 
     .       garr_n(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1)
     .       + M_w_prod*w_pcle(i)

             garr_vx(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) = 
     .       garr_vx(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) + vcurv(i,1)*M_w_prod*w_pcle(i)

             garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) = 
     .       garr_vy(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) + vcurv(i,2)*M_w_prod*w_pcle(i)

             garr_vz(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) = 
     .       garr_vz(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .             ,k_np(i)-1:k_np(i)+1) + vcurv(i,3)*M_w_prod*w_pcle(i)

             v2 = spcs(isp)%pcles(ip)%v_n(i,1)
     .           *spcs(isp)%pcles(ip)%v_n(i,1)
             garr_vx2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) = 
     .       garr_vx2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) + v2*(M_w_prod*w_pcle(i))

             v2 = spcs(isp)%pcles(ip)%v_n(i,2)
     .           *spcs(isp)%pcles(ip)%v_n(i,2)
             garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) =
     .       garr_vy2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) + v2*(M_w_prod*w_pcle(i))

             v2 = spcs(isp)%pcles(ip)%v_n(i,3)
     .           *spcs(isp)%pcles(ip)%v_n(i,3)
             garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) =
     .       garr_vz2(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) + v2*(M_w_prod*w_pcle(i))


#if defined(HYB_MASS_DIAG)
             garr_nDPIC(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1) =
     .       garr_nDPIC(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)+1)
     .            + M_w_prod*w_pcle(i)

             !accumulate vx (linear in x, 2nd in y, 2nd in z)
             garr_vxDPIC(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1
     .              ,k_np(i)-1:k_np(i)+1) =
     .       garr_vxDPIC(i_np(i)-1:i_np(i),j_np(i)-1:j_np(i)+1
     .              ,k_np(i)-1:k_np(i)+1) +
     .            vcurv(i,1)*opt_omatvec_233(Myz2,Mx1)*w_pcle(i)

             !accumulate vy (2nd in x, linear in y, 2nd in z)
             garr_vyDPIC(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)
     .            ,k_np(i)-1:k_np(i)+1) =
     .       garr_vyDPIC(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)
     .            ,k_np(i)-1:k_np(i)+1) +
     .            vcurv(i,2)*opt_omatvec_323(Mzx2,My1)*w_pcle(i)

             !accumulate vz (2nd in x, 2nd in y, linear in z)
             garr_vzDPIC(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)) =
     .       garr_vzDPIC(i_np(i)-1:i_np(i)+1,j_np(i)-1:j_np(i)+1
     .            ,k_np(i)-1:k_np(i)) +
     .            vcurv(i,3)*opt_omatvec_332(Mxy2,Mz1)*w_pcle(i)
#endif

          end do
        enddo
!$omp end do 
!$omp end parallel 

!       just the moments without b.c., the b.c. should be done where needed.
        mments(:,:,:,isp,IS11) = garr_vx2
        mments(:,:,:,isp,IS22) = garr_vy2
        mments(:,:,:,isp,IS33) = garr_vz2
        mments(:,:,:,isp,Ipx ) = garr_vx
        mments(:,:,:,isp,Ipy ) = garr_vy
        mments(:,:,:,isp,Ipz ) = garr_vz
        mments(:,:,:,isp,Ins ) = garr_n

#if defined(HYB_MASS_DIAG)
        mments(:,:,:,isp,InDPIC ) = garr_n
        mments(:,:,:,isp,IpxDPIC) = garr_vxDPIC
        mments(:,:,:,isp,IpyDPIC) = garr_vyDPIC
        mments(:,:,:,isp,IpzDPIC) = garr_vzDPIC
#endif

      enddo

c     Wrap moments around boundaries

      call wrap_moments_nodes3d(mments,wrap_currents=.true.)

c     End program

      end subroutine gather_3d_moments

c     wrap_moments_nodes3d
c     ##################################################################
      subroutine wrap_moments_nodes3d(mments,wrap_currents)
      
      implicit none

c     Call variables

      real(8), intent(inout) :: mments(0:,0:,0:,:,0:)
      logical, optional :: wrap_currents

c     local variables

      real(8),allocatable,dimension(:) :: send_buf,rec_buf
!      logical :: wcurr

c     Begin program

!      if (PRESENT(wrap_currents)) then
!        wcurr = wrap_currents
!      else
!        wcurr = .true.
!      endif

      if (pbcond(1) == pPER) then
         mments(1  ,:,:,:,:) = mments(1,:,:,:,:)+mments(nxg+1,:,:,:,:)
         mments(nxg,:,:,:,:) = mments(0,:,:,:,:)+mments(nxg  ,:,:,:,:)
         mments(nxg+1,:,:,:,:) = mments(1  ,:,:,:,:)
         mments(0    ,:,:,:,:) = mments(nxg,:,:,:,:)

      else
         call pstop('wrap_moments_nodes3d'
     .             ,'Currently only pPER for hybrid bounds')
      endif

      if (pbcond(3) == pPER) then
         mments(:  ,1,:,:,:) = mments(:,1,:,:,:)+mments(:,nyg+1,:,:,:)
         mments(:,nyg,:,:,:) = mments(:,0,:,:,:)+mments(:,nyg  ,:,:,:)
         mments(:,nyg+1,:,:,:) = mments(:,1  ,:,:,:)
         mments(:,0    ,:,:,:) = mments(:,nyg,:,:,:)
      else
         call pstop('wrap_moments_nodes3d'
     .             ,'Currently only pPER for hybrid bounds')
      endif

      if (pbcond(5) == pPER) then
         mments(:,:,1  ,:,:) = mments(:,:,1,:,:)+mments(:,:,nzg+1,:,:)
         mments(:,:,nzg,:,:) = mments(:,:,0,:,:)+mments(:,:,nzg  ,:,:)
         mments(:,:,nzg+1,:,:) = mments(:,:,1  ,:,:)
         mments(:,:,0    ,:,:) = mments(:,:,nzg,:,:)
      else
         call pstop('wrap_moments_nodes3d'
     .             ,'Currently only pPER for hybrid bounds')
      endif

#if defined(petsc)
      allocate(send_buf(size(mments))
     $        ,rec_buf (size(mments)))
      send_buf = reshape(mments, (/size(mments)/))
      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      mments(0:,0:,0:,:,0:) = reshape(rec_buf,shape(mments))
      deallocate(send_buf,rec_buf)
#endif

c     end program

      end subroutine wrap_moments_nodes3d
 
c     applyBC_moms
c     ##################################################################
      subroutine applyBC_moms

      implicit none

c     Call variables

c     Local variables

      integer :: im,isp,nx,ny,nz,igrid,order,icomp
      real(8),pointer,dimension(:,:,:) :: dummy

c     Begin program

      igrid = 1

      icomp = 30 !Make sure we do not interfere with dependent vars ! ADAM: THERE MAY BE MORE IN PIX3D.

      nx = size(moms,1)-2
      ny = size(moms,2)-2
      nz = size(moms,3)-2

      allocate(dummy(0:nx+1,0:ny+1,0:nz+1))

      dummy = 0d0

      !Densities
      do isp = 1,n_sp
        call setBC(gv%gparams,icomp,nx,ny,nz
     .          ,moms(:,:,:,isp,Ins)
     .          ,dummy
     .          ,mom_bcs(:,Ins)
     .          ,igrid,igrid,igrid,iorder=order_moms)
      enddo

      !Momenta
      do isp = 1,n_sp
        call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid
     $              ,moms(:,:,:,isp,IPx:IPz)
     $              ,mom_bcs(:,IPx:IPz)
     .              ,is_cnv=.true.,icomp=(/icomp+1/)
     .              ,iorder=order_moms)
      enddo

      !Temperatures
      do im=IS11,IS33
        do isp = 1,n_sp
          call setBC(gv%gparams,icomp,nx,ny,nz
     .              ,moms(:,:,:,isp,im)
     .              ,dummy
     .              ,mom_bcs(:,im)
     .              ,igrid,igrid,igrid,iorder=order_moms)
        enddo
      enddo


#if defined(HYB_MASS_DIAG)
      icomp=35 ! Adam: Check this.
      !Densities
      do isp = 1,n_sp
        call setBC(gv%gparams,icomp,nx,ny,nz
     .          ,moms(:,:,:,isp,InDPIC)
     .          ,dummy
     .          ,mom_bcs(:,InDPIC)
     .          ,igrid,igrid,igrid,iorder=order_moms)
      enddo

      !Momenta
      do isp = 1,n_sp
        call setMGBC_stg(gv%gparams,0,3,nx,ny,nz,igrid
     $              ,moms(:,:,:,isp,IpxDPIC:IpzDPIC)
     $              ,mom_bcs(:,IpxDPIC:IpzDPIC)
     .              ,is_cnv=.true.,icomp=(/icomp+1/)
     .              ,iorder=order_moms)
      enddo
#endif


      deallocate(dummy)

      end subroutine applyBC_moms

c     find_rho_from_pcles
c     ##################################################################
      subroutine find_rho_from_pcles(rho)

      implicit none

c     Call variables

      real(8) :: rho(0:,0:,0:)

c     Local variables

      integer :: isp
      real(8) :: rho0,ihx,ihy,ihz

c     Define current

      call get_inverse_h(ihx,ihy,ihz)

      idlvol = ihx*ihy*ihz

      rho(:,:,:) = 0d0
      do isp = 1,n_sp
        rho = rho + spcs(isp)%q*moms(:,:,:,isp,Ins)*idlvol
      enddo

      end subroutine find_rho_from_pcles

c     find_P_from_pcles
c     ############################################################
      subroutine find_P_from_pcles(pcnv)

      implicit none

c     Call variables

      real(8) :: pcnv(0:,0:,0:,:)

c     Local variables

      integer :: isp
      real(8) :: ihx,ihy,ihz

c     Define current

      call get_inverse_h(ihx,ihy,ihz)

      idlvol = ihx*ihy*ihz

      pcnv = 0d0
      do isp = 1,n_sp
        pcnv(:,:,:,1) = pcnv(:,:,:,1)
     .                + spcs(isp)%q*moms(:,:,:,isp,Ipx)*idlvol
        pcnv(:,:,:,2) = pcnv(:,:,:,2)
     .                + spcs(isp)%q*moms(:,:,:,isp,Ipy)*idlvol
        pcnv(:,:,:,3) = pcnv(:,:,:,3)
     .                + spcs(isp)%q*moms(:,:,:,isp,Ipz)*idlvol
      enddo

      end subroutine find_P_from_pcles

c     find_rho_from_pcles_DPIC
c     ##################################################################
      subroutine find_rho_from_pcles_DPIC(rho)

      implicit none

c     Call variables

      real(8) :: rho(0:,0:,0:)

c     Local variables

      integer :: isp
      real(8) :: rho0,ihx,ihy,ihz

c     Define current

      call get_inverse_h(ihx,ihy,ihz)

      idlvol = ihx*ihy*ihz

      rho(:,:,:) = 0d0
#if defined(HYB_MASS_DIAG)
      do isp = 1,n_sp
        rho = rho + spcs(isp)%m*moms(:,:,:,isp,InDPIC)*idlvol
      enddo
#endif

      end subroutine find_rho_from_pcles_DPIC

c     find_P_from_pcles_DPIC
c     ############################################################
      subroutine find_P_from_pcles_DPIC(pcnv)

      implicit none

c     Call variables

      real(8) :: pcnv(0:,0:,0:,:)

c     Local variables

      integer :: isp
      real(8) :: ihx,ihy,ihz

c     Define current

      call get_inverse_h(ihx,ihy,ihz)

      idlvol = ihx*ihy*ihz

      pcnv = 0d0
#if defined(HYB_MASS_DIAG)
      do isp = 1,n_sp
        pcnv(:,:,:,1) = pcnv(:,:,:,1)
     .                + spcs(isp)%m*moms(:,:,:,isp,IpxDPIC)*idlvol
        pcnv(:,:,:,2) = pcnv(:,:,:,2)
     .                + spcs(isp)%m*moms(:,:,:,isp,IpyDPIC)*idlvol
        pcnv(:,:,:,3) = pcnv(:,:,:,3)
     .                + spcs(isp)%m*moms(:,:,:,isp,IpzDPIC)*idlvol
      enddo
#endif

      end subroutine find_P_from_pcles_DPIC

      end module hybrid_moments

c pic_mass_matrix
c ######################################################################    
      module pic_mass_matrix ! ADAM: HAVENT MODIFIED THIS FOR 3D OR CURV YET.

      use math

      use OMP

      use pic_io

      use pic_map

      use pic_BCs

      use pic_shape_f

      use grid, ONLY: fromLocalToGlobalLimits
     .               ,fromGlobalToLocalLimits
     .               ,find_global_nobc

      use variables, ONLY: gv
      
      use mg_solver, ONLY: cSolver

      logical :: mass_matrix_solve =.false. !Whether we solve mass matrix for Equ matching

      real(8),private :: My2(1,3),Mx2(3,1)
     .                  ,My1(1,2),Mx1(2,1)
     .                  ,Mxy2(3,3)
     $                  ,Mx1y2(2,3)
      real(8),pointer,dimension(:,:),private :: sol,rr
      real(8),pointer,dimension(:,:,:,:),private :: solg

      !Mass matrix for exact density sampling
      real(8),pointer,dimension(:,:,:) :: MM2l
      integer :: it_count,iiout=0

      real(8),pointer,dimension(:,:,:,:) :: Jns,Jpz,Jpx

      real(8) :: Jnpc_avg(nsp_max)        !average Jac*np/cell

      integer(8) :: npcles(nsp_max)=0        !Number of particles
      integer    :: npcles_cell(nsp_max)=0   !Number of particles/cell

      integer,private :: n_sp,nxl,nyl,nzl,nnl,nxg,nyg,nzg,nng
      real(8),private :: ihx,ihy,ihz

      integer,allocatable,dimension(:,:) :: npc_int

      type(particle),pointer,dimension(:) :: pcles

      contains

c     setup_mass_matrix
c     ######################################################################
      subroutine setup_mass_matrix(nx_l,ny_l,nz_l,nx_g,ny_g,nz_g)

      implicit none

c     Call variables
      
      integer :: nx_g,ny_g,nz_g,nx_l,ny_l,nz_l,nsp

c     Begin program


      call get_inverse_h(ihx,ihy,ihz)

      n_sp = size(spcs)

      nxl = nx_l
      nyl = ny_l
      nzl = nz_l
      nnl = nxl*nyl*nzl

      nxg = nx_g
      nyg = ny_g
      nzg = nz_g
      nng = nxg*nyg*nzg

      allocate(MM2l(nnl,3**_Dimx,n_sp) !Global
     $        ,rr  (nnl,n_sp)          !Local
     .        ,sol (nnl,n_sp)          !Local
     .        ,solg(nxg,nyg,nzg,n_sp) )!Global

      end subroutine setup_mass_matrix

c     destroy_mass_matrix
c     ######################################################################
      subroutine destroy_mass_matrix

      deallocate(rr,sol,MM2l,solg)

      end subroutine destroy_mass_matrix

c     find_pcle_weights
c     ######################################################################
      subroutine find_pcle_weights(mments)

      implicit none

c     Call variables

      real(8), intent(in) :: mments(0:,0:,:,0:)

      integer :: bcs(6,n_sp),igl,jgl,kgl,i,j,k,ii,ip,isp,numspcs
     .          ,numpcles

      real(8) :: tol_out,sump(n_sp),lsump(n_sp)
      integer :: i_n(_Npg), j_n(_Npg), k_n(_Npg)

      external massmat_mtvc

c     Begin program

c     Build GLOBAL 9-point stencil Mass Matrix

      call build_mass_matrix(.false.)

c     Build LOCAL rhs taking Jacobian into account

      rr  = 0d0

      do k=1,nzl
        do j=1,nyl
          do i=1,nxl
            ii = i + nxl*(j-1) + nxl*nyl*(k-1)
            call fromLocalToGlobalLimits(gv%gparams,1,i,j,k,igl,jgl,kgl)
            !solve for delta weight
            rr(ii,:) = Jnpc_avg(1:n_sp)*Jns(igl,jgl,1,:)-mments(i,j,:,0)
          enddo
        enddo
      enddo

c     Solve 2D mass matrix problem

      numspcs = size(spcs)

      do isp=1,numspcs
        bcs(:,isp) = pbcond
      enddo
      where (bcs == pREF .or. bcs == pOPN) bcs = NEU

      iiout = 0 !No output
      call cSolver(gv%gparams,numspcs,nnl,rr,sol,bcs
     .            ,1,iiout,0,massmat_mtvc,.false.
     .            ,tol           = 1d-12
     .            ,gm_driver     = .true.
     .            ,ks_it         = 100
     .            !Identity preconditioner
cc     .            ,mg_vcyc       = 0
     .            !GS preconditioner
     .            ,mg_vcyc       = 1
     .            ,mg_smooth_only= .true. !Only at fine level
     .            ,mg_smooth     = 'gs'
     .
     .            ,sm_it         = 1
     .            ,sm_omega      = 1d0
     .            !Output
     .            ,iters         = it_count
     .            ,conv_tol      = tol_out
     .            )  

c$$$      do isp=1,n_sp
c$$$        call vec_plot_gl(1,rr(:,isp),1,"rhs"//trim(int2char(isp))
c$$$     .       ,"Density rhs")
c$$$        call vec_plot_gl(1,sol(:,isp),1,"dens"//trim(int2char(isp))
c$$$     .       ,"Density match")
c$$$      enddo
c$$$      stop 'after weight solve'

      !Diagnostics
      if (my_rank_pic == 0) then
        if (iiout == 0) then
          write (*,*)
          write (*,'(a,i3,a,1pe10.3)')
     .     ' Density mass matrix solve converged in ',it_count,' its'
     .    //' at tol=',tol_out
        endif
      endif

c     Gather global solution, sol => solg
      
      call find_global_nobc(reshape(sol,(/nxl,nyl,nzl,n_sp/)),solg)

c     Particle accounting

      sump = 0d0
      do k = 1,nzg
        do j = 1,nyg
          do i = 1,nxg
            ii  = i + nxg*(j-1) + nxg*nyg*(k-1)
            sump = sump + npc_int(ii,:)*solg(i,j,k,:)
          enddo
       enddo
      enddo

#if defined(petsc)
      lsump = sump
      call MPI_Allreduce(lsump,sump,size(sump)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      do isp=1,n_sp
        sump(isp) = pgroup*sump(isp)/npcles(isp)
      enddo

      if (my_rank_pic == 0) then
        write (*,*) 
        write (*,*) 'Particle weight accounting:',sump+1d0
      endif

      !Adjust particle accounting
cc      do isp=1,n_sp
cc        solg(:,:,:,isp) = solg(:,:,:,isp)/sump(isp)
cccc        Jnpc_avg(isp) = Jnpc_avg(isp)/sump(isp)
cccc        charge(isp) = sign(1d0,q_m(isp))/Jnpc_avg(isp)*(hx*hy*hz)
cccc        mass(isp) = charge(isp)/q_m(isp)
cc      enddo

c     Assign particle weights

      do isp=1, numspcs
        numpcles = size(spcs(isp)%pcles)
!$OMP PARALLEL DEFAULT(SHARED) private(ip,i,i_n,j_n)
!$OMP DO
        do ip=1,numpcles
          if(ip.le.nip) then    !reserved inflow particles
            spcs(isp)%pcles(ip)%w_n = 1d0
          else
            call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
     .                                               ,i_n,j_n,k_n)
            do i=1,_Npg
              spcs(isp)%pcles(ip)%w_n(i) = spcs(isp)%pcles(ip)%w_n(i)
     .                                    +solg(i_n(i),j_n(i),1,isp)
            enddo
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL
      end do

c     Free work space

      end subroutine find_pcle_weights

c     find_pcle_drift_velocity (vz only, which is cell-centered)
c     ######################################################################
      subroutine find_pcle_drift_velocity(mments)

      implicit none

      real(8), intent(in) :: mments(0:,0:,:,0:)

      integer :: bcs(6,n_sp)
      integer :: ci(9),ipcl,numpcles,numspcs

      integer :: i,j,k,im,ip,jm,jp,ijk,ijkim,ijkip,ijkjm,ijkjp
     .          ,ijkimjm,ijkimjp,ijkipjm,ijkipjp,igl,jgl,kgl
     .          ,ijkgl,ii,isp
      integer :: i_n(_Npg),j_n(_Npg),k_n(_Npg)

      real(8) :: send_buf(n_sp),rec_buf(n_sp),sump(n_sp),tol_out

      external massmat_mtvc

c     Begin program

      if (maxval(abs(Jpz(:,:,1,:)))==0d0) return

c     Build GLOBAL 9-point stencil matrix

      call build_mass_matrix(.true.)

c     Build LOCAL rhs taking Jacobian into account

      rr  = 0d0

      do k=1,nzl
        do j=1,nyl
          do i=1,nxl
            ii = i + nxl*(j-1) + nxl*nyl*(k-1)
            call fromLocalToGlobalLimits(gv%gparams,1,i,j,k,igl,jgl,kgl)
            rr(ii,:) = Jnpc_avg(1:n_sp)*Jpz(igl,jgl,1,:)-mments(i,j,:,3) !No charge included!
          enddo
        enddo
      enddo

c     Solve 2D mass matrix problem

      numspcs = size(spcs)

      do isp=1,numspcs
        bcs(:,isp) = pbcond
      enddo
      where (bcs == pREF) bcs = NEU
      where (bcs == pOPN) bcs = DIR

      iiout = 0 !No output
      call cSolver(gv%gparams,numspcs,nnl,rr,sol,bcs
     .            ,1,iiout,0,massmat_mtvc,.false.
     .            ,tol           = 1d-12
     .            ,gm_driver     = .true.
     .            ,ks_it         = 100
     .            !Identity preconditioner
cc     .            ,mg_vcyc       = 0
     .            !GS preconditioner
     .            ,mg_vcyc       = 1
     .            ,mg_smooth_only= .true. !Only at fine level
     .            ,mg_smooth     = 'gs'
     .
     .            ,sm_it         = 1
     .            ,sm_omega      = 1d0
     .            !Output
     .            ,iters         = it_count
     .            ,conv_tol      = tol_out)  

c$$$      do isp=1,n_sp
c$$$        call vec_plot_gl(1,rr(:,isp),1,"rhs"//trim(int2char(isp))
c$$$     .       ,"Pz rhs")
c$$$        call vec_plot_gl(1,sol(:,isp),1,"pz"//trim(int2char(isp))
c$$$     .       ,"Pz match")
c$$$      enddo
c$$$      stop 'after drift vel solve'

      if (my_rank_pic == 0) then
        if (iiout == 0) then
          write (*,*)
          write (*,'(a,i3,a,1pe10.3)')
     .        ' Current mass matrix solve converged in ',it_count,' its'
     .       //' at tol=',tol_out
        endif
      endif

c     Gather global solution, sol => solg

      call find_global_nobc(reshape(sol,(/nxl,nyl,nzl,n_sp/)),solg)

c     Assign particle drift velocity vz (Stardard curvilinear convention, i.e., no Jac)

      do isp=1, numspcs 
        numpcles = size(spcs(isp)%pcles)
!$OMP PARALLEL DEFAULT(SHARED) private(ip,i,i_n,j_n)
!$OMP DO
        do ip=1,numpcles
          if(ip.le.nip) then    !reserved inflow particles
            spcs(isp)%pcles(ip)%v_n (:,3) = 0d0
            spcs(isp)%pcles(ip)%v_np(:,3) = spcs(isp)%pcles(ip)%v_n(:,3)
          else
            call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
     .                                               ,i_n,j_n,k_n)
            do i=1, _Npg
              spcs(isp)%pcles(ip)%v_n(i,3) =spcs(isp)%pcles(ip)%v_n(i,3)
     .                                     + solg(i_n(i),j_n(i),1,isp)
cc     .                                     -v0_z(isp)
              spcs(isp)%pcles(ip)%v_np(i,3)=spcs(isp)%pcles(ip)%v_n(i,3)
            end do
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL 
      end do

      end subroutine find_pcle_drift_velocity

c     find_pcle_drift_velocity_x (vx only)
c     ######################################################################
      subroutine find_pcle_drift_velocity_x(mments)

      implicit none

      real(8), intent(in) :: mments(0:,0:,:,0:)

      integer :: bcs(6,n_sp)
      integer :: ci(9),ipcl,numpcles,numspcs

      integer :: i,j,k,im,ip,jm,jp,ijk,ijkim,ijkip,ijkjm,ijkjp
     .          ,ijkimjm,ijkimjp,ijkipjm,ijkipjp,igl,jgl,kgl
     .          ,ijkgl,ii,isp
      integer :: i_n(_Npg),j_n(_Npg),k_n(_Npg)

      real(8) :: send_buf(n_sp),rec_buf(n_sp),sump(n_sp),tol_out
     .           ,xp
      integer :: inp,jnp,imp,jmp

      external massmat_mtvc

c     Begin program

      if (maxval(abs(Jpx(:,:,1,:)))==0d0) return

c     Build GLOBAL 9-point stencil matrix

      call build_mass_matrix_x(.true.)

c     Build LOCAL rhs taking Jacobian into account

      rr  = 0d0

      do k=1,nzl
        do j=1,nyl
          do i=1,nxl
            ii = i + nxl*(j-1) + nxl*nyl*(k-1)
            call fromLocalToGlobalLimits(gv%gparams,1,i,j,k,igl,jgl,kgl)
            rr(ii,:) = Jnpc_avg(1:n_sp)*Jpx(igl,jgl,1,:)-mments(i,j,:,1) 
          enddo
        enddo
      enddo

c     Solve 2D mass matrix problem

      numspcs = size(spcs)

      do isp=1,numspcs
        bcs(:,isp) = pbcond
      enddo
      where (bcs == pREF) bcs = NEU
      where (bcs == pOPN) bcs = DIR

      iiout = 0 !No output
      call cSolver(gv%gparams,numspcs,nnl,rr,sol,bcs
     .            ,1,iiout,0,massmat_mtvc,.false.
     .            ,tol           = 1d-12
     .            ,gm_driver     = .true.
     .            ,ks_it         = 100
     .            !Identity preconditioner
cc     .            ,mg_vcyc       = 0
     .            !GS preconditioner
     .            ,mg_vcyc       = 1
     .            ,mg_smooth_only= .true. !Only at fine level
     .            ,mg_smooth     = 'gs'
     .
     .            ,sm_it         = 1
     .            ,sm_omega      = 1d0
     .            !Output
     .            ,iters         = it_count
     .            ,conv_tol      = tol_out)  

c$$$      do isp=1,n_sp
c$$$        call vec_plot_gl(1,rr(:,isp),1,"rhs"//trim(int2char(isp))
c$$$     .       ,"Pz rhs")
c$$$        call vec_plot_gl(1,sol(:,isp),1,"pz"//trim(int2char(isp))
c$$$     .       ,"Pz match")
c$$$      enddo
c$$$      stop 'after drift vel solve'

      if (my_rank_pic == 0) then
        if (iiout == 0) then
          write (*,*)
          write (*,'(a,i3,a,1pe10.3)')
     .        ' Current mass matrix solve converged in ',it_count,' its'
     .       //' at tol=',tol_out
        endif
      endif

c     Gather global solution, sol => solg

      call find_global_nobc(reshape(sol,(/nxl,nyl,nzl,n_sp/)),solg)

c     Assign particle drift velocity vz (Stardard curvilinear convention, i.e., no Jac)

      do isp=1, numspcs 
        numpcles = size(spcs(isp)%pcles)
!$OMP PARALLEL DEFAULT(SHARED) private(ip,i,i_n,j_n,imp)
!$OMP DO
        do ip=1,numpcles
          if(ip > nip) then    !reserved inflow particles
            call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
     .                                               ,i_n,j_n,k_n)
            do i=1, _Npg
              imp = i_n(i)-1
              !periodic b.c.
              if(imp==0) imp = nxg
              xp = spcs(isp)%pcles(ip)%x_n(i,1)*ihx
              if(xp>0.5d0) then
                 spcs(isp)%pcles(ip)%v_n(i,1) =
     .                    spcs(isp)%pcles(ip)%v_n(i,1)
     .                   +solg(i_n(i),j_n(i),1,isp)
              else
                 spcs(isp)%pcles(ip)%v_n(i,1) =
     .                    spcs(isp)%pcles(ip)%v_n(i,1)
     .                   +solg(imp   ,j_n(i),1,isp)
              end if
              spcs(isp)%pcles(ip)%v_np(i,1)
     .                  = spcs(isp)%pcles(ip)%v_n(i,1)
	    end do
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL
      enddo

      end subroutine find_pcle_drift_velocity_x

c     build_mass_matrix
c     ######################################################################
      subroutine build_mass_matrix(include_pcle_weight)

      implicit none

c     ----------------------------------------------------------------------
c     Build GLOBAL 9-point stencil Mass Matrix
c     ----------------------------------------------------------------------

c     Call variables

      logical :: include_pcle_weight

c     Local variables

      integer :: ci(9),ipcl,isp,numspcs,numpcles

      integer :: i,j,k,im,ip,jm,jp,ijk,ijkim,ijkip,ijkjm,ijkjp
     .          ,ijkimjm,ijkimjp,ijkipjm,ijkipjp,igl,jgl,kgl
     .          ,ijkgl,i_npg
      real(_prec) :: w_pcle,xp,yp

      real(8) :: My2(3,1),Mx2(1,3),My1(1,2),Mx1(1,2),Mxy2(3,3)

      real(8),allocatable,dimension(:) :: send_buf,rec_buf

      real(8),allocatable,dimension(:,:,:) :: MM2

c     Begin program

      allocate(MM2(nng,3**_Dimx,n_sp)) !Global

      MM2 = 0d0

      numspcs = size(spcs)

      do isp=1,numspcs
        numpcles = size(spcs(isp)%pcles)
!$omp parallel DEFAULT(SHARED) private(ipcl,i,j,ip,im,jp,jm,ci,ijk
!$omp.         ,ijkim,ijkip,ijkjm,ijkjp,ijkimjm,ijkimjp,ijkipjm
!$omp.         ,ijkipjp,Mx2,My2,Mxy2,xp,yp,i_npg,w_pcle)
!$         call set_omp_thread_id()
!$omp DO 
!$omp.reduction(+:MM2) 
         do ipcl=nip+1,numpcles
           do i_npg=1,_Npg
              w_pcle = spcs(isp)%pcles(ipcl)%w_n(i_npg)

              xp = spcs(isp)%pcles(ipcl)%x_n(i_npg,1)*ihx
              yp = spcs(isp)%pcles(ipcl)%x_n(i_npg,2)*ihy
              Mx2(1,:) = shape2_offset_scalar(xp)
              My2(:,1) = shape2_offset_scalar(yp)
#ifdef __GNU
              Mxy2(:,:)= opt_vecvec_3(My2,Mx2)
#else
              Mxy2(:,:)= matmul(My2,Mx2)
#endif

              if (include_pcle_weight) Mxy2 = Mxy2*w_pcle

              call xform_pcle_idx_scl(spcs(isp)%pcles(ipcl)%ijk_n(i_npg)
     .                               ,nxg,nyg,i,j,k)

              ip = i + 1
              im = i - 1
              jp = j + 1
              jm = j - 1

              !Set up the contributing matrix index
              ci(7)=3; ci(8)=2; ci(9)=1
              ci(4)=6; ci(5)=5; ci(6)=4
              ci(1)=9; ci(2)=8; ci(3)=7

              if (i==1) then
                if (pbcond(1)==pPER) then
                  im = nxg
                else
                  im = i
                  ci(1) = ci(1)-2
                  ci(4) = ci(4)-2
                  ci(7) = ci(7)-2
                endif              
              endif

              if (i==nxg) then
                if (pbcond(2)==pPER) then
                  ip = 1
                else
                  ip = i
                  ci(3) = ci(3)+2
                  ci(6) = ci(6)+2
                  ci(9) = ci(9)+2                
                endif
              endif

              if (j==1) then
                if (pbcond(3)==pPER) then
                  jm = nyg
                else
                  jm = j
                  ci(1) = ci(1)-6 
                  ci(2) = ci(2)-6 
                  ci(3) = ci(3)-6                 
                endif
              endif

              if (j==nyg) then
                if (pbcond(4)==pPER) then
                  jp = 1
                else
                  jp = j
                  ci(7) = ci(7)+6 
                  ci(8) = ci(8)+6 
                  ci(9) = ci(9)+6 
                endif
              endif
            
              ijk = i + nxg*(j-1)

              ijkip   = ijk + ip-i
              ijkim   = ijk + im-i
              ijkjp   = ijk + (jp-j)*nxg
              ijkjm   = ijk + (jm-j)*nxg
              ijkimjm = ijk + im-i + (jm-j)*nxg
              ijkipjm = ijk + ip-i + (jm-j)*nxg
              ijkimjp = ijk + im-i + (jp-j)*nxg
              ijkipjp = ijk + ip-i + (jp-j)*nxg

              MM2(ijkipjp,ci(9),isp)=MM2(ijkipjp,ci(9),isp)+Mxy2(3,3)
              MM2(ijkjp  ,ci(8),isp)=MM2(ijkjp  ,ci(8),isp)+Mxy2(2,3)
              MM2(ijkimjp,ci(7),isp)=MM2(ijkimjp,ci(7),isp)+Mxy2(1,3)
                                                     
              MM2(ijkip  ,ci(6),isp)=MM2(ijkip  ,ci(6),isp)+Mxy2(3,2)
              MM2(ijk    ,ci(5),isp)=MM2(ijk    ,ci(5),isp)+Mxy2(2,2)
              MM2(ijkim  ,ci(4),isp)=MM2(ijkim  ,ci(4),isp)+Mxy2(1,2)
                                                     
              MM2(ijkipjm,ci(3),isp)=MM2(ijkipjm,ci(3),isp)+Mxy2(3,1)
              MM2(ijkjm  ,ci(2),isp)=MM2(ijkjm  ,ci(2),isp)+Mxy2(2,1)
              MM2(ijkimjm,ci(1),isp)=MM2(ijkimjm,ci(1),isp)+Mxy2(1,1)
            end do
         end do
!$omp end do 
!$omp end parallel 
      end do

      !Wrap GLOBAL mass matrix
#if defined(petsc)
      allocate(send_buf(size(mm2))
     $        ,rec_buf (size(mm2)))
      send_buf = reshape(mm2, (/size(mm2)/))
      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      mm2 = reshape(rec_buf,shape(mm2))
      deallocate(send_buf,rec_buf)
#endif

      !Assign LOCAL mass matrix
      do j=1,nyl
        do i=1,nxl
          call fromLocalToGlobalLimits(gv%gparams,1,i,j,1,igl,jgl,k)

          ijkgl = igl + nxg*(jgl-1)
          ijk   = i   + nxl*(j  -1)
             
          MM2l(ijk,:,:) = MM2(ijkgl,:,:)
        enddo
      enddo

      deallocate(MM2)

      end subroutine build_mass_matrix

c     build_mass_matrix_x
c     ######################################################################
      subroutine build_mass_matrix_x(include_pcle_weight)

      implicit none

c     ----------------------------------------------------------------------
c     Build GLOBAL 9-point stencil Mass Matrix for staggered grid
c     component in x (e.g. jx)
c     ----------------------------------------------------------------------

c     Call variables

      logical :: include_pcle_weight

c     Local variables

      integer :: ci(9),ipcl,isp,numspcs,numpcles

      integer :: i,j,k,im,ip,jm,jp,ijk,ijkim,ijkip,ijkjm,ijkjp
     .          ,ijkimjm,ijkimjp,ijkipjm,ijkipjp,igl,jgl,kgl
     .          ,ijkgl,i_npg
      real(_prec) :: w_pcle,xp,yp

      real(8) :: My2(3,1),Mx2(1,3),My1(1,2),Mx1(1,2),Mxy2(3,3)

      real(8),allocatable,dimension(:) :: send_buf,rec_buf

      real(8),allocatable,dimension(:,:,:) :: MM2

c     Begin program

      allocate(MM2(nng,3**_Dimx,n_sp)) !Global

      MM2 = 0d0

      numspcs = size(spcs)

      do isp=1,numspcs
        numpcles = size(spcs(isp)%pcles)
!$omp parallel DEFAULT(SHARED) private(ipcl,i,j,ip,im,jp,jm,ci,ijk
!$omp.         ,ijkim,ijkip,ijkjm,ijkjp,ijkimjm,ijkimjp,ijkipjm
!$omp.         ,ijkipjp,Mx2,My2,Mxy2,xp,yp,i_npg,w_pcle)
!$         call set_omp_thread_id()
!$omp DO 
!$omp.reduction(+:MM2) 
         do ipcl=nip+1,numpcles
           do i_npg=1,_Npg
              w_pcle = spcs(isp)%pcles(ipcl)%w_n(i_npg)

              xp = spcs(isp)%pcles(ipcl)%x_n(i_npg,1)*ihx
              yp = spcs(isp)%pcles(ipcl)%x_n(i_npg,2)*ihy
              Mx2(1,:) = shape2_offset_scalar(xp)
              My2(:,1) = shape2_offset_scalar(yp)
#ifdef __GNU
              Mxy2(:,:)= opt_vecvec_3(My2,Mx2)
#else
              Mxy2(:,:)= matmul(My2,Mx2)
#endif

              if (include_pcle_weight) Mxy2 = Mxy2*w_pcle

              call xform_pcle_idx_scl(spcs(isp)%pcles(ipcl)%ijk_n(i_npg)
     .                               ,nxg,nyg,i,j,k)

              ip = i + 1
              im = i - 1
              jp = j + 1
              jm = j - 1

              !Set up the contributing matrix index
              ci(7)=3; ci(8)=2; ci(9)=1
              ci(4)=6; ci(5)=5; ci(6)=4
              ci(1)=9; ci(2)=8; ci(3)=7

              if (i==1) then
                if (pbcond(1)==pPER) then
                  im = nxg
                else
                  print *,"May not be correct building mass-matrix 
     .                     with relfectiv b.c.. Stop."
                  stop
                  im = i
                  ci(1) = ci(1)-2
                  ci(4) = ci(4)-2
                  ci(7) = ci(7)-2
                endif              
              endif

              if (i==nxg) then
                if (pbcond(2)==pPER) then
                  ip = 1
                else
                  ip = i
                  ci(3) = ci(3)+2
                  ci(6) = ci(6)+2
                  ci(9) = ci(9)+2                
                endif
              endif

              if (j==1) then
                if (pbcond(3)==pPER) then
                  jm = nyg
                else
                  jm = j
                  ci(1) = ci(1)-6 
                  ci(2) = ci(2)-6 
                  ci(3) = ci(3)-6                 
                endif
              endif

              if (j==nyg) then
                if (pbcond(4)==pPER) then
                  jp = 1
                else
                  jp = j
                  ci(7) = ci(7)+6 
                  ci(8) = ci(8)+6 
                  ci(9) = ci(9)+6 
                endif
              endif
            
              ijk = i + nxg*(j-1)

              ijkip   = ijk + ip-i
              ijkim   = ijk + im-i
              ijkjp   = ijk + (jp-j)*nxg
              ijkjm   = ijk + (jm-j)*nxg
              ijkimjm = ijk + im-i + (jm-j)*nxg
              ijkipjm = ijk + ip-i + (jm-j)*nxg
              ijkimjp = ijk + im-i + (jp-j)*nxg
              ijkipjp = ijk + ip-i + (jp-j)*nxg

              if(xp>0.5d0) then   !    J=j
                MM2(ijkjp  ,ci(8),isp) = MM2(ijkjp  ,ci(8),isp)
     .                                 + Mx1y2(2,3)
                MM2(ijkimjp,ci(7),isp) = MM2(ijkimjp,ci(7),isp)
     .                                 + Mx1y2(1,3)
                                                           
                MM2(ijk    ,ci(5),isp) = MM2(ijk    ,ci(5),isp)
     .                                 + Mx1y2(2,2)
                MM2(ijkim  ,ci(4),isp) = MM2(ijkim  ,ci(4),isp)
     .                                 + Mx1y2(1,2)
                                                           
                MM2(ijkjm  ,ci(2),isp) = MM2(ijkjm  ,ci(2),isp)
     .                                 + Mx1y2(2,1)
                MM2(ijkimjm,ci(1),isp) = MM2(ijkimjm,ci(1),isp)
     .                                 + Mx1y2(1,1)

              else               !     J=j-1

                MM2(ijkjp  ,ci(8)-1,isp) = MM2(ijkjp  ,ci(8)-1,isp)
     .                                   + Mx1y2(2,3)
                MM2(ijkimjp,ci(7)-1,isp) = MM2(ijkimjp,ci(7)-1,isp)
     .                                   + Mx1y2(1,3)
                                                           
                MM2(ijk    ,ci(5)-1,isp) = MM2(ijk    ,ci(5)-1,isp)
     .                                   + Mx1y2(2,2)
                MM2(ijkim  ,ci(4)-1,isp) = MM2(ijkim  ,ci(4)-1,isp)
     .                                   + Mx1y2(1,2)
                                                           
                MM2(ijkjm  ,ci(2)-1,isp) = MM2(ijkjm  ,ci(2)-1,isp)
     .                                   + Mx1y2(2,1)
                MM2(ijkimjm,ci(1)-1,isp) = MM2(ijkimjm,ci(1)-1,isp)
     .                                   + Mx1y2(1,1)
                          
              end if

            end do
         end do
!$omp end do 
!$omp end parallel 
      end do

      !Wrap GLOBAL mass matrix
#if defined(petsc)
      allocate(send_buf(size(mm2))
     $        ,rec_buf (size(mm2)))
      send_buf = reshape(mm2, (/size(mm2)/))
      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      mm2 = reshape(rec_buf,shape(mm2))
      deallocate(send_buf,rec_buf)
#endif

      !Assign LOCAL mass matrix
      do j=1,nyl
        do i=1,nxl
          call fromLocalToGlobalLimits(gv%gparams,1,i,j,1,igl,jgl,k)

          ijkgl = igl + nxg*(jgl-1)
          ijk   = i   + nxl*(j  -1)
             
          MM2l(ijk,:,:) = MM2(ijkgl,:,:)
        enddo
      enddo

      deallocate(MM2)

      end subroutine build_mass_matrix_x

      end module pic_mass_matrix

c pic_setup
c ######################################################################    
      module pic_setup ! Adam: May want to create hybrid_setup.

      use hybrid_moments

      use pic_mass_matrix

      use grid, ONLY: npx,npy,npz,gparams,np
#if defined(petsc)
     .               ,add_profile_stage
#endif

      use ts_io, ONLY:time,dt,ts_chk_dump

      use var_io, ONLY: tst_flg

      real(8) :: omega_p (nsp_max)=1d0    !Electron plasma freq
     .          ,mass    (nsp_max)=1d0    !Pcle mass
     .          ,charge  (nsp_max)=1d0    !Pcle charge
     .          ,q_m     (nsp_max)=-1d0   !Charge-mass ratio

      real(8),private :: w_bp (nsp_max)  !weight of background particles

      real(8),allocatable,dimension(:) :: pxx,pyy,pzz
      real(8),allocatable,dimension(:) :: npc_acm
      integer,allocatable,dimension(:,:) :: mnpc_int

      integer,private :: nps(nsp_max)
      integer,private :: nxg,nyg,nzg,nng,nxl,nyl,nzl,nnl
      real(8),private :: hx,hy,hz,dh(3),Lx,Ly,Lz
      real(8),private :: ihx,ihy,ihz,ih(3)
      real(8),private :: mach_eps,two_mach_eps
      real(8),private :: xmax,xmin,ymax,ymin,zmax,zmin

      character*(20),private :: test = ''

      integer :: pcle_push_log  !Handle for profiling

      logical :: exp_pic

      !Diagnostics
      integer :: npa_int
     .          ,npc1_0
     .          ,npl                  !Lost particles

      !Input parameters defaults
      integer :: n_sp=nsp_max              !Number of species
     .          ,min_pcles=0               !Minimun number of particles

      integer :: n_b=0                     !Number and index of beam species
      integer :: i_b=0
      real(8) :: r_b=0d0                   !Ratio of the beam

      real(8) :: C_bkgrnd(nsp_max)=0d0     !Density percentage of background particles

      real(8) :: v_thx(nsp_max)=0d0        !Thermal velocity
     .          ,v_thy(nsp_max)=0d0        !Thermal velocity
     .          ,v_thz(nsp_max)=0d0        !Thermal velocity
     .          ,v0_x (nsp_max)=0d0        !Drift velocity
     .          ,v0_y (nsp_max)=0d0        !Drift velocity
     .          ,v0_z (nsp_max)=0d0        !Drift velocity
     .          ,theta=0d0                 !Angle B forms with x-axis, for rotated anisotropic distributions

      real(8) :: nsubstep(nsp_max)=1       !for fixed number of substeps
     .          ,eps_omg=0d0               !Perturbation frequency 
     .          ,eps_pic=0d0               !Perturbation of particle positions
     .          ,eps_v  =0d0               !Perturbation in velocity

      logical :: phase_space_plot =.false. !Whether to plot phase space or not
     .          ,quiet_start      =.false. !For using quiet start 
     .          ,eps_prop         =.false. !Perturbation propagate
     .          ,particle_sort    =.true.  !Whether to sort particles

      contains

c     readPICInput
c     #################################################################
      function readPICInput(file,dump) result(ioerr)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ioerr
      character(*) :: file
      logical :: dump

c     Local variables

      integer :: uinput=1234,isp
      real(8) :: lam

c     Namelists

      namelist /picdef/ npcles,v0_x,v0_y,v0_z,v_thx,v_thy,v_thz,theta
     .     ,n_sp,omega_p,q_m,min_pcles,npcles_cell
     .     ,quiet_start,nsubstep,eps_prop,eps_omg,eps_pic,eps_v
     .     ,n_b,i_b,r_b,nip,vinject,C_bkgrnd
     .     ,phase_space_plot,mass_matrix_solve,particle_sort

c     Begin program

      ioerr = 0

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank_pic, mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, mpierr)
#endif

c     Read initialization parameters

      do while(is_open_unit(uinput))
        uinput = uinput + 1
      enddo

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,picdef,iostat=ioerr)

      close(unit=uinput)

      if (ioerr/=0) return

c     Consistency check

      !>>Numbers
      if (n_sp == 0) call pstop('readInput','Zero number of species')

      select case(test)
      case('tm')
        lam = gparams(2)
!        write (*,*) v0_z(1),v0_z(2)
c$$$        if (v0_z(2) /= v_thz(2)/lam .or. v0_z(1) /= -v0_z(2)
c$$$     .      .or. v0_z(1) == 0d0 .or. v0_z(2) == 0d0) then
c$$$          if (my_rank_pic == 0)
c$$$     .         write (*,*) "WARNING: problem with v0_z setup"
c$$$        endif
      end select

c$$$      if (min_pcles > 0) then
c$$$        call pstop('readInput','min_pcles > 0; '//
c$$$     .          'need to fix issue with min_pcles and drift velocity')
c$$$      endif

c     Write initialization parameters to standard output

      if (my_rank_pic == 0.and.dump) then
        write (*,*) '-----------------'
        write (*,*) 'PIC configuration'
        write (*,*) '-----------------'
        write (*,picdef)
      endif

c     PETSc profiling

#if defined(petsc)
      pcle_push_log = add_profile_stage("Pcle push")
#endif

c     End program

      end function readPICInput

c$$$c     get_spline_order
c$$$c     ##################################################################
c$$$      function get_spline_order() result(ss)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: ss
c$$$
c$$$c     Begin program
c$$$
c$$$      ss = sp_order
c$$$
c$$$c     End program
c$$$
c$$$      end function get_spline_order

c     setup_pic_3d
c     ##################################################################
      subroutine setup_pic_3d(inputfile,problem,nxin,xin,nyin
     $                       ,yin,nzin,zin,explicit_pic,no__map,dt
     .                       ,jac,xmap,pcle_rstrt)

      implicit none

c     Call variables

      character(*) :: problem
      integer :: nxin,nyin,nzin
      real(8) :: xin(0:nxin+1),yin(0:nyin+1),zin(0:nzin+1),dt
      logical :: explicit_pic,pcle_rstrt,no__map

      real(8), intent(in) :: jac   (0:,0:,0:)
     .                      ,xmap  (0:,0:,0:,:)

      character(*) :: inputfile

c     Local variables

      integer :: groupn,ip,ierr,numspcs,isp
      real(8) :: n0,tot_bp

      type(particle),pointer,dimension(:) :: pcles

c     Begin program

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank_pic, mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, mpierr)
#endif

      test = problem

c     Initialize local variables

      exp_pic = explicit_pic

      !find machine epsilon
#if _prec == 4
      mach_eps=findRoundOff_sp()
#else
      mach_eps=findRoundOff()
#endif
      two_mach_eps = 2d0*mach_eps 

c     Get number of OMP threads and provide output

      if (my_rank_pic == 0) call check_omp_numthreads(.not.tst_flg)

c     Import the grid info

c           [i=1]
c     x---|===x===|---x ... |---x---|---x
c     0   0   1   1   2         nx  nx  nx+1
c        xmin                      xmax           
      nxg = nxin
      nxl = nxg/npx

      nyg = nyin
      nyl = nyg/npy

      nzg = nzin
      nzl = nzg/npz
      
      nng = nxg*nyg*nzg
      nnl = nxl*nyl*nzl

      !Allocate cell center position arrays (logical space)
      allocate(pxx(0:nxg+1))
      pxx = xin
cc      allocate(pxx(0:nxg))
cc      pxx = 0.5*(xin(0:nxg)+xin(1:nxg+1))

      allocate(pyy(0:nyg+1))
      pyy = yin
cc      allocate(pyy(0:nyg))
cc      pyy = 0.5*(yin(0:nyg)+yin(1:nyg+1))

      allocate(pzz(0:nzg+1))
      pzz = zin
cc      allocate(pzz(0:nzg+1))
cc      pzz = 0.5*(zin(0:nzg)+zin(0:nzg+1))

c     Initialize mesh quantities

      call setup_map(nxg,nyg,nzg,jac,no__map,xmap,xin,yin,zin
     .              ,xmin,xmax,ymin,ymax,zmin,zmax,hx,hy,hz) !Limits, and hx/y/z are output

c     Initialize mesh lengths and volumes

      Lx = xmax-xmin
      Ly = ymax-ymin
      Lz = zmax-zmin

      ihx = 1d0/hx
      ihy = 1d0/hy
      ihz = 1d0/hz

      dh = (/hx,hy,hz/)
      ih = (/ihx,ihy,ihz/)

c     Initialize pcle shape f

      call setup_pic_shape_f(ihx,ihy,ihz)

c     Initialize pcle accounting qtys

      if (sum(npcles) == 0) npcles = npcles_cell*nxg*nyg*nzg

      gnip = nip

      do isp = 1,n_sp
        if (npcles(isp) <= gnip) then
          call pstop('readInput'
     .              ,'Pcle buffer exceeds total pcle budget')
        endif
      enddo

      if(quiet_start) then
!     group of Antithetic variables
         groupn = 8
      else
         groupn = 1
      end if

      pgroup = max(groupn,_Npg)

      if (_Npg > 8.and.quiet_start) then
        call pstop('setup_pic_2d',' _Npg > 8')
      endif

c     Find PIC physical quantities

      allocate(Jns(0:nxg+1,0:nyg+1,0:nzg+1,n_sp)
     $        ,Jpz(0:nxg+1,0:nyg+1,0:nzg+1,n_sp)
     $        ,Jpx(0:nxg+1,0:nyg+1,0:nzg+1,n_sp))

      do isp=1,n_sp
        Jns(:,:,:,isp)  = moms_g(:,:,:,isp,Ins)*jac(:,:,:)
      enddo
      Jpz = moms_g(:,:,:,:,IPz)
      Jpx = moms_g(:,:,:,:,IPx)

      do isp=1,n_sp
         omega_p(isp) = sqrt(abs(q_m(isp))) !careful
         tot_bp = min_pcles*nng*pgroup
         nps(isp) = npcles(isp) - gnip - tot_bp

         Jnpc_avg(isp) = (nps(isp))/sum(Jns(1:nxg,1:nyg,1:nzg,isp))
         n0 = Jnpc_avg(isp)/(hx*hy*hz)     !Add volume factor

c$$$         charge(isp) = sign(1d0,q_m(isp))/n0
c$$$         mass(isp) = charge(isp)/q_m(isp)
         mass(isp) = 1d0/n0
         charge(isp) = mass(isp)*q_m(isp)

         if (tot_bp > 0) then
           w_bp(isp) = nps(isp)*C_bkgrnd(isp)/tot_bp !This stems from defining n_b = C_b*int(Jns)/int(J)
         else
           w_bp(isp) = 0d0
         endif
      enddo

c     Allocate particles

      call alloc_pcles(n_sp,npcles)

c     Initialize particle quantities

      !Setup shape functions
      call setup_pic_shape_f(ihx,ihy,ihz)

      !Allocate particle accumulators
      allocate(npc_acm(0:nng+1),npc_int(nng,n_sp),mnpc_int(nng,n_sp))

      npc_acm = 0 ; npc_int = 0; mnpc_int = 0 ; npl  = 0

      m_np = 0d0
      m_n  = 0d0

      PoynFlux = 0d0

c     Initialize particles

      if (pcle_rstrt) then
        ierr = pic_read_pcles(time)

        if (ierr /= 0) then
          write (*,*) 'Error reading restart file'
          stop
        endif

!        call gather_moments
        call gather_moments ! Adam: Flag here to switch.
        momsn = moms
      else
        numspcs = size(spcs)
        do isp=1,numspcs
          spcs(isp)%q   = charge(isp)
          spcs(isp)%m   = mass(isp)
          spcs(isp)%q_m = q_m(isp)
        enddo

        call load_pcles_offset
      endif

cc      call setup_moms(nxg,nyg,nzg,nxl,nyl,nzl,n_sp,hx,hy,hz)

c     End program

      end subroutine setup_pic_3d

c     kill_pic
c     ##################################################################
      subroutine kill_pic

      call dealloc_pcles

      call dealloc_moms

      deallocate(npc_acm,npc_int,mnpc_int)
      
      deallocate(pxx,pyy,pzz)

      deallocate(Jns,Jpz,Jpx)

      call dealloc_map

      if (.not.no_map) call killSplines

      end subroutine kill_pic

c$$$c     defineBCs_PIC
c$$$c     ##################################################################
c$$$      subroutine defineBCs_PIC
c$$$
c$$$      implicit none
c$$$
c$$$      call defineBCs_pcles
c$$$
c$$$      call defineBCs_moms
c$$$
c$$$      end subroutine defineBCs_PIC

c     load_pcles_offset  
c     ##################################################################
      subroutine load_pcles_offset

      use iso_c_binding, only: c_int, c_double
      use F_call_C
      
      use hammersley_seq
      
      implicit none
c     ------------------------------------------------------------------
c     Load particles in phase space
c     ------------------------------------------------------------------

c     Call variables

c     Local variables

      real(8) :: v_0x (1:n_sp+n_b)
     .          ,v_0y (1:n_sp+n_b)
     .          ,v_0z (1:n_sp+n_b)

      integer :: ip,ip_ng,ipl,tot_np,ipc,is,ipb
     .           ,dimt,ipc0,ncall,i,j,k,ii,cii,offset
     .           ,ipx,ipy,ipz,ixyz,np3,npc_max(n_sp),igl,jgl,kgl
     .           ,tot_np_scan,npc_scan(nng,n_sp),numspcs,numpcles
      real(8) :: rx1(_Npg),rx2(_Npg),rx3(_Npg)
     $          ,rx (_Npg),ry (_Npg),rz (_Npg)
     .          ,xp(_Npg),yp(_Npg),zp(_Npg),dxp,vx2,vt2
     .          ,signx,signy,signz,pp,pr,betaj
      real(8) :: vpar(_Npg),vperp(_Npg)
      real(8) :: v_tot(1:n_sp+n_b,_Dim)
      integer, save :: seed
      real(8) :: vt(3), v0(3)
      real(8),allocatable,dimension(:,:) :: r4

c     Begin program

      v_0x = v0_x(1:n_sp+n_b)
      v_0y = v0_y(1:n_sp+n_b)
      v_0z = v0_z(1:n_sp+n_b)

      call init_random_seed(nproc,my_rank_pic,.not.tst_flg)

      dimt = 6!_Dim+_Dimx ! 2D3V !ADAM - used in problem_init files.
      numspcs = size(spcs)
      do is=1,numspcs
         call compute_pcle_spat_dist(is,npc_int(:,is),npc_scan(:,is)
     $                              ,npc_max(is),mnpc_int(:,is))
      end do

      if(quiet_start) then

        select case(test)
        case('ic','tm')
#include "problem_init/icLoadPartls_qs.f"
        case('slow_shock')
          call pstop("load_pcles_offset"
     .              ,"Slow shock quiet start not implemented")
        case('kaw_ii')
#include "problem_init/kawiiLoadPartls_qs.f"
        case('landau','2strm','langdon') ! took out 'iaw'
#include "problem_init/UniformLoadPartls_qs.f"
        case('aniso')
#include "problem_init/AnisoLoadPartls_qs.f"
        case default  
#include "problem_init/DefaultLoadPartls_qs.f"


        end select

      else                      !quiet_start

        select case(test)
        case('ic','tm')
#include "problem_init/icLoadPartls.f"
        case('slow_shock')
#include "problem_init/slowShockLoadPartls.f"
        case('kaw_ii')   
#include "problem_init/kawiiLoadPartls.f"
        case('landau','2strm','iaw','langdon')! took out 'iaw'.
#include "problem_init/UniformLoadPartls.f"
        case default
#include "problem_init/DefaultLoadPartls.f"
        end select

      end if                    !without quiet_start

      call free_random_seed()

      contains

cc      include 'problem_init/qsRandn.f'

c     compute_pcle_spat_dist
c     ##################################################################
      subroutine compute_pcle_spat_dist(isp,npc_int_s,npc_scan_s
     $                                 ,npc_max_s,mnpc_int_s)

      implicit none

c     ------------------------------------------------------------------
c     compute particle spatial distribution (LOCAL particles, GLOBAL mesh)
c     ------------------------------------------------------------------

      integer :: npc_max_s,npc_int_s(:),npc_scan_s(:),isp,gtot_np,npc
     .          ,mnpc_int_s(:)
      logical :: low_pcles = .false.
      integer :: bnpc_avg,gtot_bnp, gnpc,npc_r,npc_r0,mnpc_r
      real(8) :: gnpc_avg
      
c     Begin program

      tot_np = size(spcs(isp)%pcles)
      tot_np = tot_np*_Npg
      gtot_np = npcles(isp)    !all particles         

!     Always add some uniform base particles

      bnpc_avg = min_pcles      !per pgroup*cell
      gtot_bnp = bnpc_avg*nng*pgroup
      betaj = 1d0
      gnpc_avg = (gtot_np-gnip-gtot_bnp)/
     $            sum(Jns(1:nxg,1:nyg,1:nzg,isp)) ! ADAM: MADE 3D.
!     .           *jac(1:nxg,1:nyg,1:nzg)**(betaj-1d0)) 

      if (gnpc_avg.le.0) then
        write (*,'(a,i4)')
     .       'Not enough particles to satisfy minimum per cell of'
     .       ,min_pcles*pgroup
        write (*,'(a,i8,a)')
     .       'Need at least', min_pcles*nxg*nyg*nzg*pgroup,' pcles'
        stop
      endif

!     first sweep accumulates (per pgroup) globally

      npc_acm(0) = gnip/pgroup
      do k = 1,nzg
         do j = 1,nyg
            do i = 1,nxg
               ii  = i + nxg*(j-1) + nxg*nyg*(k-1)
               npc_acm(ii) = npc_acm(ii-1)
     .                     +(gnpc_avg/pgroup)*Jns(i,j,k,isp)
     .                     + bnpc_avg! *jac(i,j,k)**(betaj-1d0) + bnpc_avg
            end do
         end do
      end do      
      
!     Check the total accumulated particles

      if(nint(npc_acm(ii)).ne.gtot_np/pgroup) then
         write (*,*) ii,nng,gnip,nip,npc_acm(ii),gtot_np/pgroup
         print *
     $        ,"The accumulated number of particles does not add up"
         print *," to the prescribed one, stop."
         stop
      end if

!     second sweep allocate particles to each cell

      npc_r0 = 0
      npc_max_s = 0 
      do k = 1,nzg
         do j = 1,nyg
            do i = 1,nxg
               ii  = i + nxg*(j-1) + nxg*nyg*(k-1)
               gnpc = nint(npc_acm(ii)) - nint(npc_acm(ii-1))
               npc_int_s(ii) = gnpc / nproc
               mnpc_int_s(ii)= min_pcles / nproc
               npc_r = mod(gnpc,nproc)
               mnpc_r= mod(min_pcles,nproc)
               if((     (my_rank_pic-npc_r0)>=0
     .             .and.(my_rank_pic-npc_r0)<npc_r)
     $        .or.(     (my_rank_pic      -npc_r0)<0
     .             .and.(my_rank_pic+nproc-npc_r0)<npc_r)) then
                  npc_int_s(ii) = npc_int_s(ii) + 1
               end if
               if((     (my_rank_pic-npc_r0)>=0
     .             .and.(my_rank_pic-npc_r0)<mnpc_r)
     $        .or.(     (my_rank_pic      -npc_r0)<0
     .             .and.(my_rank_pic+nproc-npc_r0)<mnpc_r)) then
                  mnpc_int_s(ii) = mnpc_int_s(ii) + 1
               end if

               npc_r0 = mod(npc_r0+npc_r,nproc)
               if (npc_max_s<npc_int_s(ii)) npc_max_s = npc_int_s(ii)
cc               if(ii==8) print *,my_rank_pic,ii,mnpc_int_s(ii),mnpc_r
            end do
         end do
      end do


cc      print *,"after arranging particles:",my_rank_pic,sum(npc_int_s),tot_np

#if defined(petsc)
!     find the index of counting particles per mpi process
      call MPI_SCAN(npc_int_s,npc_scan_s,nng,MPI_INTEGER,MPI_SUM
     $     ,MPI_COMM_WORLD,mpierr)         
      do ii = 1,nng
         npc_scan_s(ii) = npc_scan_s(ii) - npc_int_s(ii) +
     $        nint(npc_acm(ii-1)) 
      end do
#else
      do ii = 1,nng
         npc_scan_s(ii) = nint(npc_acm(ii-1))
      end do
#endif

      npc1_0 = npc_int_s(1) +npc_int_s(2)
      
      end subroutine compute_pcle_spat_dist

      end subroutine load_pcles_offset

c     perturb_pcles_offset
c     ##################################################################
      subroutine perturb_pcles_offset(nh1,nh2)

      implicit none

c     Call variables

      integer :: nh1,nh2

c     Local variables

      integer :: isp,ip,ii,jj,i,j,k,im,numspcs,numpcles
      integer :: i_np(_Npg), j_np(_Npg), k_np(_Npg)
      real(_prec) :: kx,kxp,xn,xp,v2,xp2,yp,kyp,ky,dxp
     .              ,x_np(_Npg)

      real(8) :: vcurv(3),cov(3,3),xpp(3),xlog(3)
      real(8) :: const_pi

c     Begin program

      const_pi = pi

c     Initialize particle new info

      numspcs = size(spcs)
      do isp=1, numspcs
         numpcles = size(spcs(isp)%pcles)
!     note that there is some issue with valgrind and libgomp.
!$omp parallel DEFAULT(SHARED) private(ip)
!$       call set_omp_thread_id()
!$OMP DO
         do ip=1,numpcles
           spcs(isp)%pcles(ip)%ijk_np = spcs(isp)%pcles(ip)%ijk_n
           spcs(isp)%pcles(ip)%x_np   = spcs(isp)%pcles(ip)%x_n
           spcs(isp)%pcles(ip)%v_np   = spcs(isp)%pcles(ip)%v_n
           spcs(isp)%pcles(ip)%w_np   = spcs(isp)%pcles(ip)%w_n
        end do
!$OMP END DO
!$omp end parallel 
      end do

c     Find particle weights and drift velocities

      if (mass_matrix_solve) then
        call gather_moments

        call setup_mass_matrix(nxl,nyl,nzl,nxg,nyg,nzg)

        call find_pcle_weights(moms(:,:,1,1:n_sp,:))

        call gather_moments

        call find_pcle_drift_velocity(moms(:,:,1,1:n_sp,:))
c        call find_pcle_drift_velocity_x(moms(:,:,1,1:n_sp,:))

        call destroy_mass_matrix
      end if

C    perturb particles after matching given equilibrium

      select case(trim(test))
      case('ic','tm')
!#include "problem_init/icPerturb.f"
      case('dtm','dgem')
      case('slow_shock')
!#include "problem_init/slowShockPerturb.f"
      case('kaw_ii')
        if(quiet_start) then
!#include "problem_init/kawiiPerturb.f"
        end if
c$$$      case('bpnch')
c$$$        if(quiet_start) then
c$$$            call perturb_pcle_xv(1
c$$$     .                          ,(/nh1,nh2,0/)
c$$$     .                          ,(/0d0,0d0,0d0/)
c$$$     .                          ,(/0d0,0d0,eps_pic/))
c$$$        end if
      case('weibel','wbl_bc','aniso','bpnch')
        if(quiet_start) then
           print*,'eps_v=',eps_v
            call perturb_pcle_xv(1
     .                          ,(/nh1,nh2,0/)
     .                          ,(/0d0,0d0,eps_v/)
     .                          ,(/0d0,0d0,eps_pic/),delta_fun=.true.)
        end if
      case('2strm','iaw','iash','msw','landau','langdon','whsl2')
        do isp=1,n_sp
          call perturb_pcle_xv(isp
     .                        ,(/nh1,nh2,0/)
     .                        ,(/eps_pic,0d0,0d0/)
     .                        ,(/eps_v,0d0,0d0/),delta_fun=.false.)!(/0d0,0d0,eps_v/))
        enddo
      case default
        print *, "test is not set correctly in the input, "
     .         //"quit perturb_pcles_offset."
        stop
      end select

c     Gather moments

c     call gather_moments
      call gather_moments ! Adam: Flag here.

c     Diagnostics and memory freeing

      if (my_rank_pic == 0) then
         write (*,*)
         print *,'Total number of particles (npcles-nip)',
     .           ' calculated from zeroth moment:'
         do isp=1,n_sp
cc            print *,size(spcs(isp)%pcles),charge(isp),mass(isp)            
           write (*,'(a,i3,f20.10,a,f20.10)')
     .          ' >>Species=',isp,sum(moms_g(1:nxg,1:nyg,1:nzg,isp,Ins))
     .         ,' vs. expected',(1+c_bkgrnd(isp))*nps(isp)
         enddo
      end if

c     End program

      contains

c     perturb_pcle_xv
c     ######################################################################
      subroutine perturb_pcle_xv(isp,mode,eps_pos,eps_vel,delta_fun)

      implicit none

c     Call variables

      integer :: isp,mode(3),dim
      real(8) :: eps_pos(3),eps_vel(3)
      logical,optional :: delta_fun

c     Local variables

      real(_prec) :: k0(3),xg(3),xp(3),kxp(3),dx(_Npg,3),dv(_Npg,3)
      integer :: ip,numpcles,j
      logical :: delta

c     Begin program

      if (PRESENT(delta_fun)) then
        delta = delta_fun
      else
        delta = .false.
      endif

      numpcles = size(spcs(isp)%pcles)

      k0 = 2d0*pi*dble(mode)/(/Lx,Ly,Lz/)

      if (my_rank == 0) then
        print*,'delta=',delta
        print*,'eps_vel=',eps_vel
      endif

      do ip=1,numpcles

        call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg
     .                                        ,i_np,j_np,k_np)

        !VECTOR SIMD
        do j=1,_Npg

           xg = (/pxx(i_np(j)),pyy(j_np(j)),pzz(k_np(j))/)

           xp = spcs(isp)%pcles(ip)%x_n(j,:) + xg  - 0.5d0*dh

           kxp = k0*xp      

           where (kxp /= 0d0) 
             dx(j,:) = eps_pos*cos(kxp)
           elsewhere
             dx(j,:) = 0d0
           end where
           spcs(isp)%pcles(ip)%x_n(j,:) = spcs(isp)%pcles(ip)%x_n(j,:)
     .                                  + dx(j,:)

           if(delta) then
             if (i_np(j)==nxg/2+1) then!.and.j_np(j)==nyg/2+1) then
               spcs(isp)%pcles(ip)%v_n(j,:)
     .              = spcs(isp)%pcles(ip)%v_n(j,:) + eps_vel(:)
               spcs(isp)%pcles(ip)%v_np(j,:)
     .              = spcs(isp)%pcles(ip)%v_n(j,:)
             endif
           else
!             where (kxp(1) /= 0d0)
              select case(trim(test))
              case('whsl2')
               dv(j,2) = -eps_vel(1)*cos(kxp(1))
               dv(j,3) = eps_vel(1)*sin(kxp(1))

              case default
               dv(j,:) = eps_vel*sin(kxp(1))!eps_vel*cos(kxp)
              end select

 !            elsewhere
 !              dv(j,:) = 0d0
 !            end where

 !            print*,'kxp=',kxp
 !            print*,'dv=',dv

             spcs(isp)%pcles(ip)%v_n(j,:)
     .            = spcs(isp)%pcles(ip)%v_n(j,:) + dv(j,:)
             spcs(isp)%pcles(ip)%v_np(j,:)
     .            = spcs(isp)%pcles(ip)%v_n(j,:)
           end if

        end do

        call cell_crossing(dx,spcs(isp)%pcles(ip)%x_n
     .                       ,spcs(isp)%pcles(ip)%v_n,i_np,j_np,k_np)

        spcs(isp)%pcles(ip)%ijk_n = i_np+nxg*(j_np-1)+nxg*nyg*(k_np-1)        

        spcs(isp)%pcles(ip)%x_np =  spcs(isp)%pcles(ip)%x_n

      end do

      end subroutine perturb_pcle_xv

      end subroutine perturb_pcles_offset

c     cell_crossing
c     ##################################################################
      subroutine cell_crossing(dx,xnp,vnp,inp,jnp,knp)

      implicit none

c     Call variables

      real(_prec),intent(in)    :: dx(:,:)
      real(_prec),intent(inout) :: xnp(:,:)
      real(_prec),intent(inout) :: vnp(:,:)
      integer,intent(inout) :: inp(:),jnp(:),knp(:)


c     Local variables

      integer :: nsize(3),in0(size(inp),3),i,dim,nxg,nyg,nzg
      real(8) :: hx,hy,hz,dh(3)

c     Begin program

      call get_global_grid_size(nxg,nyg,nzg)
      call get_h(hx,hy,hz)

      dh = (/hx,hy,hz/)

      nsize = (/nxg,nyg,nzg/)
      in0(:,1)   = inp
      in0(:,2)   = jnp
      in0(:,3)   = knp

c     Impose BCs

      do dim=1,3!_Dimx
        do i=1,_Npg
          if(dx(i,dim)>0 .and. xnp(i,dim).ge.dh(dim)) then
            do while(xnp(i,dim).ge.dh(dim))
              in0(i,dim) = in0(i,dim)+1
              xnp(i,dim) = xnp(i,dim) - dh(dim)
            enddo
            if(in0(i,dim)>nsize(dim)) !BC @ location=1
     .         call pcles_bc_offset(dim,1,pbcond,nsize(dim),xnp(i,:)
     .                             ,vnp(i,:),in0(i,:))
          else if(dx(i,dim)<0.and. xnp(i,dim).le.0) then
            do while(xnp(i,dim).le.0d0)
              in0(i,dim) = in0(i,dim)-1
              xnp(i,dim) = xnp(i,dim) + dh(dim)
            enddo
            if(in0(i,dim)<1)         !BC @ location=0
     .         call pcles_bc_offset(dim,0,pbcond,nsize(dim),xnp(i,:)
     .                             ,vnp(i,:),in0(i,:))
          end if
        enddo
      end do

c     Transfer location info

      inp = in0(:,1)
      jnp = in0(:,2)
      knp = in0(:,3)

      end subroutine cell_crossing

c     pcles_bc_offset
c     ##################################################################
      subroutine pcles_bc_offset(dim,loc,pbc,nn,xp,vp,ip)

      implicit none

c     ------------------------------------------------------------------
c     Impose particle boundary conditions along dimension dim=(1,2,3)
c     at location loc (0, left; 1, right). Same convention as framework.
c     ------------------------------------------------------------------

c     Call variables

      integer,intent(in)   :: pbc(6),dim,loc,nn
      real(8),intent(inout):: xp(3),vp(3)
      integer,intent(inout):: ip(3)

c     Local variables

      integer :: ibc,dim1,dim2
      real(8) :: hx,hy,hz,ihx,ihy,ihz,dh(3)
      real(_prec) :: xlog(3),cov(3,3),vdn,inmag

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      call get_h(hx,hy,hz)
      call get_inverse_h(ihx,ihy,ihz)
      dh = (/hx,hy,hz/)

      select case(dim)
      case(1)
        dim1 = 2 ; dim2 = 3
      case(2)
        dim1 = 1 ; dim2 = 3
      case(3)
        dim1 = 1 ; dim2 = 2
      end select

      select case(pbc(ibc))
      case(pPER)
!     periodic
         if (loc == 0) ip(dim) = nn
         if (loc == 1) ip(dim) = 1  
      case(pREF,pSP)
!     reflective, force it to be on the boundary, before which xp=(0:eps).
         if (loc == 0) ip(dim) = ip(dim) + 1
         if (loc == 1) ip(dim) = ip(dim) - 1 
         xp(dim) = loc*dh(dim)
         if (no_map) then
           vp(dim) =-vp(dim)
         else
           xlog(1) = xp(1)*ihx
           xlog(2) = xp(2)*ihy
           xlog(3) = xp(3)*ihz

           call transf_vector_Bspl_scalar(ip(1),ip(2),ip(3),xlog,cov) 
           inmag = 1d0/sum(cov(dim,:)*cov(dim,:))  ! |n|
           vdn = dot_product(cov(dim,:),vp)        ! v.n
cc           if (vdn*(-1)**loc < 0d0) vp = vp - 2d0*vdn*cov(dim,:)*inmag ! vp = vp - 2*(v.n)n/|n|^2
           vp = vp - 2d0*vdn*cov(dim,:)*inmag ! vp = vp - 2*(v.n)n/|n|^2
         endif
      case(pOPN)
         if (loc == 0) ip(dim) = ip(dim) + 1
         if (loc == 1) ip(dim) = ip(dim) - 1 
         xp(dim) = loc*dh(dim)
         vp(dim) = 2d0*vinject(dim,ibc) - vp(dim) !assuming a symmetric distribution
         vp(dim1) = -vp(dim1)
         vp(dim2) = -vp(dim2)
      case default
         print *, "Specify the particle BC, quit"
         stop
      end select

c     End program

      end subroutine pcles_bc_offset

      end module pic_setup

c pic_diag
c ######################################################################    
      module pic_diag

      use pic_setup

      contains

c     init_pcle_diag
c     ###################################################################
      subroutine init_pcle_diag(Acovg,B0,cm,lm)

      implicit none

c     Call variables

      real(8) :: Acovg(0:,0:,0:,:),B0(3)
      real(8) :: cm(3),lm(3)

c     Local variables

      real(8) :: lm_np(9),gm_np(9)

c     Begin program

      !Canonical momenta
#if defined(PIC_CM_DIAG)
      cm = canon_moment(Acovg(:,:,1,:))
#else
      cm = 0d0
      cm_ex = 0d0
#endif

      !Linear momenta
      lm = 0d0
      lm(1) = moments(Ipx)

      !Parallel transfer
#if defined(petsc)
      lm_np(1:3) = lm
      lm_np(4:6) = cm
      lm_np(7:9) = cm_ex
      call MPI_AllReduce(lm_np,gm_np,9,MPI_DOUBLE_PRECISION
     $                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      lm = gm_np(1:3)
      cm = gm_np(4:6)
      cm_ex = gm_np(7:9)
#endif

c     End program

      end subroutine init_pcle_diag

c     moment
c     ##################################################################
      function moment(pcles,order,iv) result(m)

      implicit none

c     Call variables

      integer :: order,iv
      type(particle),dimension(:) :: pcles
      real(_prec) :: m

c     Local variables

      integer :: ip, numpcles, i

c     Begin program

      m = 0d0

      numpcles = size(pcles)

      do ip=nip+1,numpcles
        do i=1, _Npg
           m = m + pcles(ip)%v_np(i,iv)**order * pcles(ip)%w_np(i)
        end do
      enddo

c     End program

      end function moment

c     moments
c     ##################################################################
      function moments(id) result(m)

      implicit none

c     Call variables

      integer :: id
      real(_prec) :: m

c     Local variables

      integer :: is,numspcs

c     Begin program

      m = 0d0
      
      numspcs = size(spcs)

      select case(id)
      case(Ins)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,0,1)*spcs(is)%q
        enddo
      case(Ipx)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,1,1)*spcs(is)%m
        enddo
      case(Ipy)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,1,2)*spcs(is)%m
        enddo
      case(Ipz)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,1,3)*spcs(is)%m
        enddo
      case(IS11)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,2,1)*spcs(is)%m
        enddo
      case(IS22)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,2,2)*spcs(is)%m
        enddo
      case(IS33)
        do is=1,numspcs
          m = m + moment(spcs(is)%pcles,2,3)*spcs(is)%m
        enddo
      end select

c     End program

      end function moments

#if defined(PIC_CM_DIAG)
c     canon moment
c     ##################################################################
      function canon_moment(A) result(cm)

      implicit none

c     Call variables 

      real(8), intent(in) :: A(0:,0:,1:)
      real(8) :: cm(3)

c     Local variables

      integer :: ip,is,i,nxg,nyg,nzg
      integer :: numspcs,numpcles
      integer :: i0(_Npg),j0(_Npg),k0(_Npg)
      real(_prec) :: xp(_Npg),yp(_Npg),cml(3),q_sp,m_sp
      real(8) :: ihx,ihy,ihz

c     Begin program

      call get_global_grid_size(nxg,nyg,nzg)
      call get_inverse_h(ihx,ihy,ihz)

      if (exp_pic) then
        cm = 0d0
        return
      else
        cml = 0d0
        numspcs = size(spcs)
        do is=1, numspcs
          numpcles = size(spcs(is)%pcles)    
          q_sp = spcs(is)%q
          m_sp = spcs(is)%m
!$OMP PARALLEL DEFAULT(SHARED) private(ip,xp,yp,i0,j0,i)
!$OMP DO 
!$OMP.REDUCTION (+:cml)
          do ip=nip+1,numpcles
            xp = spcs(is)%pcles(ip)%x_np(:,1)*ihx
            yp = spcs(is)%pcles(ip)%x_np(:,2)*ihy
            call xform_pcle_idx(spcs(is)%pcles(ip)%ijk_np,nxg,nyg
     .                                                  ,i0,j0,k0)
            do i=1,_Npg
               cml = cml + cmp(q_sp,i0(i),j0(i),xp(i),yp(i),A) 
               cml = cml + m_sp*spcs(is)%pcles(ip)%v_np(i,:)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL 
          cm = cml
        end do
      end if

c     End program

      end function canon_moment

c     cmp
c     ##################################################################
      function cmp(q,i0,j0,xp,yp,A) result(cm)
c     ------------------------------------------------------------------
c     canonical momentum per particle
c     ------------------------------------------------------------------

      implicit none

c     Call variables 

      real(8), intent(in) :: A(0:,0:,1:)
      real(8) :: cm(3)
      integer :: i0, j0
      real(_prec) :: xp, yp, q

c     Local variables

      real(8) :: A_1(2,3),A_2(3,2),A_3(3,3),Ap
      real(8) :: Sx1(1,2),Sx2(1,3),Sy1(1,2),Sy2(1,3)

c     Begin program

      Sx1(1,:) = shape1_offset_scalar(xp)
      Sy1(1,:) = shape1_offset_scalar(yp)
      Sy2(1,:) = shape2_offset_scalar(yp)
      Sx2(1,:) = shape2_offset_scalar(xp)

      A_1(:,:) = A(i0-1:i0  , j0-1:j0+1, 1)
      A_2(:,:) = A(i0-1:i0+1, j0-1:j0  , 2)
      A_3(:,:) = A(i0-1:i0+1, j0-1:j0+1, 3)

#ifdef __GNU
      Ap = dot_product(opt_matvec_23(A_1,Sy2(1,:)),Sx1(1,:))
cc      cm(1) =  pcle%v_np(1)*m + Ap*q
      cm(1) =  Ap*q

      Ap = dot_product(opt_matvec_2(Sx2(1,:),A_2),Sy1(1,:))
cc      cm(2) =  pcle%v_np(2)*m + Ap*q               
      cm(2) =  Ap*q

      Ap = dot_product(opt_matvec_tr(Sx2(1,:),A_3),Sy2(1,:))
cc      cm(3) =  pcle%v_np(3)*m + Ap*q 
      cm(3) =  Ap*q
#else               
      Ap = dot_product(MATMUL(A_1,Sy2(1,:)),Sx1(1,:))
cc      cm(1) =  pcle%v_np(1)*m + Ap*q
      cm(1) =  Ap*q

      Ap = dot_product(MATMUL(Sx2(1,:),A_2),Sy1(1,:))
cc      cm(2) =  pcle%v_np(2)*m + Ap*q               
      cm(2) =  Ap*q

      Ap = dot_product(MATMUL(Sx2(1,:),A_3),Sy2(1,:))
cc      cm(3) =  pcle%v_np(3)*m + Ap*q 
      cm(3) =  Ap*q
#endif

      end function cmp

c$$$c     canonical momentum change due to external B0
c$$$c     ##################################################################
c$$$      function canon_moment_exchg(B0) result(cmc)
c$$$
c$$$c     ------------------------------------------------------------------
c$$$c     Adjusts canonical momenta diagnostic for external magnetic field
c$$$c     contributions. DOES NOT WORK IN PERIODIC MAPPED MESHES.
c$$$C     ------------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables 
c$$$
c$$$      real(8) :: cmc(3),B0(3)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: ip,is,numspcs,numpcles,i
c$$$      real(_prec) :: cmcl(3),q_sp
c$$$
c$$$c     Begin program
c$$$
c$$$      if (exp_pic) then
c$$$         cmc = 0d0
c$$$      else if(.not.no_map
c$$$     .        .and.bcPER(1).and.bcPER(2).and.bcPER(3)) then
c$$$         cmc = 0d0
c$$$      else
c$$$         cmcl = 0d0
c$$$
c$$$         numspcs = size(spcs)
c$$$         do is=1,numspcs
c$$$            numpcles = size(spcs(is)%pcles)
c$$$            q_sp = spcs(is)%q
c$$$!$OMP PARALLEL DEFAULT(SHARED) private(ip)
c$$$!$OMP DO 
c$$$!$OMP.REDUCTION (+:cmcl)
c$$$            do ip=nip+1,numpcles
c$$$              write (*,*) spcs(is)%pcles(ip)%d_x
c$$$              !VECTOR SIMD
c$$$              do i=1,_Npg
c$$$                 cmcl(2) = cmcl(2)
c$$$     .                  +q_sp*B0(3)*spcs(is)%pcles(ip)%d_x(i,1)
c$$$     .                  -q_sp*B0(1)*spcs(is)%pcles(ip)%d_x(i,3)
c$$$                 cmcl(3) = cmcl(3)
c$$$     .                  +q_sp*B0(1)*spcs(is)%pcles(ip)%d_x(i,2)
c$$$     .                  -q_sp*B0(2)*spcs(is)%pcles(ip)%d_x(i,1)
c$$$              end do
c$$$            end do
c$$$!$OMP END DO
c$$$!$OMP END PARALLEL 
c$$$         end do
c$$$
c$$$         cmc = cmcl
c$$$      end if
c$$$
c$$$c     End program                                                  
c$$$
c$$$      end function canon_moment_exchg
#endif

c     find vd T parallel and perpendicular to B0 as a function of x
c     ###################################################################
      subroutine find_vdTx(nxg,NDim,nS,theta,vd,T)

      implicit none

c     Call variables

      integer, intent(in) :: NDim,nS,nxg
      real(8), intent(in) :: theta
      real(8), intent(out) :: vd(nxg,NDim,nS),T(nxg,NDim,nS)
      
c     Local variables

      integer :: i,j,k,isp,ip,Np(nS),ci(_Npg),Npc(nxg,nS),n_sp
      real(8) :: vx(_Npg),vy(_Npg),vz(_Npg)
      real(8) :: v(_Npg,3)
      real(8) :: ms(nS)

      real(8),allocatable,dimension(:) :: send_buf,rec_buf

c     Begin program

      n_sp = size(spcs)

      vd = 0d0
      T  = 0d0
      Npc= 0
      
      do isp = 1,n_sp
        ms(isp) = spcs(isp)%m
        Np(isp) = size(spcs(isp)%pcles)
        do ip = nip+1, size(spcs(isp)%pcles)
          ci = spcs(isp)%pcles(ip)%ijk_np
          vx = spcs(isp)%pcles(ip)%v_np(:,1)
          vy = spcs(isp)%pcles(ip)%v_np(:,2)
          vz = spcs(isp)%pcles(ip)%v_np(:,3)
          v(:,1) = vx*cos(theta) + vz*sin(theta) !parallel
          v(:,2) = vy
          v(:,3) = vz*cos(theta) - vx*sin(theta) !perpendicular
          !VECTOR SIMD
          do k=1,_Npg
            Npc(ci(k),isp) = Npc(ci(k),isp)+1
            vd(ci(k),:,isp) = vd(ci(k),:,isp) + v(k,:)
            T (ci(k),:,isp) = T (ci(k),:,isp) + v(k,:)*v(k,:)
          enddo
        end do                  !ip
      end do                    !isp

#if defined(petsc)
      allocate(send_buf(size(T))
     $        ,rec_buf (size(T)))
      send_buf = reshape(T, (/size(T)/))
      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      T = reshape(rec_buf,shape(T))

      send_buf = reshape(vd, (/size(vd)/))
      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      vd = reshape(rec_buf,shape(vd))

      send_buf = reshape(Npc, (/size(Npc)/))
      call MPI_Allreduce(send_buf,rec_buf,size(send_buf)
     $                  ,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      Npc = reshape(rec_buf,shape(Npc))

      deallocate(send_buf,rec_buf)
#endif

      do i=1,nxg
        do j=1,NDim        
          vd(i,j,1:n_sp) = vd(i,j,1:n_sp)/Npc(i,1:n_sp)
          T (i,j,1:n_sp) =(T (i,j,1:n_sp)/Npc(i,1:n_sp)
     .                    -vd(i,j,1:n_sp)*vd(i,j,1:n_sp))*ms(1:n_sp)
        end do                  !j
      end do                    !i

c     end program
      
      end subroutine find_vdTx

c$$$c     find total vd T parallel and perpendicular to B0
c$$$c     ###################################################################
c$$$      subroutine find_vdT(spcs,vd,T,NDim,nS,theta)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      type(species) :: spcs(:)
c$$$      integer, intent(in) :: NDim,nS
c$$$      real(8), intent(in) :: theta
c$$$      real(8), intent(out) :: vd(NDim,nS),T(NDim,nS)
c$$$      
c$$$c     Local variables
c$$$
c$$$      integer :: i,isp,ip,Np(nS)
c$$$      type(particle),pointer,dimension(:) :: pcles
c$$$      real(8) :: vx,vy,vz
c$$$      real(8) :: v(3) ! parall, y, perp
c$$$      real(8) :: ms(nS)
c$$$
c$$$c     start program
c$$$
c$$$      vd = 0d0
c$$$      T  = 0d0
c$$$
c$$$      do isp = 1,n_sp
c$$$         pcles => spcs(isp)%pcles
c$$$         ms(isp  ) = spcs(isp)%m
c$$$         if(isp==i_b) then
c$$$            Np(isp  ) = ipb*2
c$$$            Np(isp+1) = size(pcles)-Np(isp)
c$$$            ms(isp+1) = spcs(isp)%m
c$$$         else
c$$$            Np(isp) = size(pcles)
c$$$         end if
c$$$         do ip = 1, size(pcles) !/2
c$$$            vx = pcles(ip)%v_np(1)
c$$$            vy = pcles(ip)%v_np(2)
c$$$            vz = pcles(ip)%v_np(3)
c$$$            v(1) = vx*cos(theta) + vz*sin(theta)
c$$$            v(2) = vy
c$$$            v(3) = -vx*sin(theta) + vz*cos(theta) !careful
c$$$            if(isp==i_b.and.ip>ipb) then
c$$$               vd(:,isp+1) = vd(:,isp+1) + v
c$$$               T (:,isp+1) = T (:,isp+1) + v*v
c$$$            else
c$$$               vd(:,isp) = vd(:,isp) + v
c$$$               T (:,isp) = T (:,isp) + v*v
c$$$            end if
c$$$
c$$$c$$$            ip2 = size(pcles)-ip+1
c$$$c$$$            vx = pcles(ip2)%v_np(1)
c$$$c$$$            vy = pcles(ip2)%v_np(2)
c$$$c$$$            vz = pcles(ip2)%v_np(3)
c$$$c$$$            v(1) = vx*cos(theta) + vz*sin(theta)
c$$$c$$$            v(2) = vy
c$$$c$$$            v(3) = -vx*sin(theta) + vz*cos(theta)
c$$$c$$$            if(isp==i_b.and.ip>ipb) then
c$$$c$$$               vd(:,isp+1) = vd(:,isp+1) + v
c$$$c$$$               T (:,isp+1) = T (:,isp+1) + v*v
c$$$c$$$            else
c$$$c$$$               vd(:,isp) = vd(:,isp) + v
c$$$c$$$               T (:,isp) = T (:,isp) + v*v
c$$$c$$$            end if
c$$$
c$$$         end do
c$$$      end do
c$$$
c$$$      do i = 1,NDim
c$$$         vd(i,:) = vd(i,:)/Np
c$$$         T (i,1:nS)=(T (i,1:nS)/Np(1:nS)-vd(i,1:nS)*vd(i,1:nS))*ms(1:nS)
c$$$      end do
c$$$
c$$$c     end program
c$$$c$$$      print *, T
c$$$c$$$      print *, sqrt(T(1:NDim,1)/ms(1)), sqrt(T(1:NDim,2)/ms(2)), 
c$$$c$$$     .     sqrt(T(1:NDim,3)/ms(3))
c$$$c$$$      stop
c$$$      end subroutine find_vdT

      end module pic_diag

c pic_sort
c ######################################################################    
      module pic_sort

      use pic_diag

      contains

c     sort_particles
c     ##################################################################
      subroutine sort_particles(nxg,nng,n_sp,Npart_cell)

      implicit none

c     ------------------------------------------------------------------
c     Sort particles by cell index
c     ------------------------------------------------------------------

c     Call variables

      integer, intent(in) :: nxg,nng,n_sp
      integer,optional,intent(out) :: Npart_cell(nng,n_sp)

c     Local variables

      integer :: isp,ipc,ipl,ip,ijk(_Npg),i,k,kl,dum
      integer :: part_cell_index(nng)
     .          ,locl_cell_index(nng) !assume that particles see the whole domain
      integer :: ti_0,ti_1,Np     
      real(8) :: vnp
      integer :: nthreads,tid,ct

c     Begin program

      if (my_rank_pic == 0.and.(.not.tst_flg)) then
        write (*,'(a)',ADVANCE = "NO") "Sorting particles..."
      endif

      do isp=1,size(spcs)

        part_cell_index = 0
        locl_cell_index = 0

        Np = size(spcs(isp)%pcles)

        !Scan first

!$OMP PARALLEL DEFAULT(SHARED) firstprivate(locl_cell_index) 
!$OMP. private(tid,nthreads,ct,i)
cc!$      nthreads = omp_get_num_threads()
cc!$      tid      = omp_get_thread_num()
        call set_omp_thread_id()
        call get_omp_thread_id(tid,nthreads)
!$OMP DO private(ijk,k) schedule(static,2) 
        do ipc=nip+1,Np
          ijk = spcs(isp)%pcles(ipc)%ijk_np
          do k=1,_Npg
            locl_cell_index(ijk(k)) = locl_cell_index(ijk(k))+1
          enddo
        end do
!$OMP END DO nowait 

!$OMP critical !one thread at a time
        part_cell_index = part_cell_index  + locl_cell_index
!$OMP end critical

!$OMP barrier

        !Prefix sum
!$OMP single
        do i = 2, nng
          part_cell_index(i) = part_cell_index(i)
     $                        +part_cell_index(i-1)
        end do
!$OMP end single

        do ct = nthreads-1, 0, -1
          if(ct == tid) then
            do i = 1, nng
              part_cell_index(i) = part_cell_index(i)
     $                            -locl_cell_index(i) 
              locl_cell_index(i) = part_cell_index(i)+1
            end do
          else
!$OMP barrier
          end if
        end do

        !Copy next
!$OMP DO PRIVATE(ijk,ip,ipl,k,kl) schedule(static,2)
        do ipc=nip+1,Np
          ijk = spcs(isp)%pcles(ipc)%ijk_np
            
          do k=1,_Npg
            ip = locl_cell_index(ijk(k))

            call xform_pcle_idx_sort(ip,_Npg,kl,ipl)
            ipl = ipl + nip

            !Save old positions and velocities only
            spcs(isp)%pcles(ipl)%x_n(kl,:)
     .           = spcs(isp)%pcles(ipc)%x_np(k,:)
            spcs(isp)%pcles(ipl)%v_n(kl,:)
     .           = spcs(isp)%pcles(ipc)%v_np(k,:)
            spcs(isp)%pcles(ipl)%ijk_n(kl)
     .           = spcs(isp)%pcles(ipc)%ijk_np(k)
            spcs(isp)%pcles(ipl)%w_n(kl)
     .           = spcs(isp)%pcles(ipc)%w_np(k)

            locl_cell_index(ijk(k))=locl_cell_index(ijk(k))+1
          enddo
        end do
!$OMP END DO
!$OMP END PARALLEL 

        if(present(Npart_cell)) Npart_cell(:,isp) = part_cell_index
         
      end do !species

      if (my_rank_pic == 0.and.(.not.tst_flg)) write (*,*) "Done!"

      end subroutine sort_particles

      end module pic_sort

c massmat_mtvc
c##################################################################
      subroutine massmat_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
c------------------------------------------------------------------
c     This subroutine is a mass matrix matvec, y = A.x.
c     In call:
c      * g_def: grid structure defs
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnd: boundary conditions on x vector.
c------------------------------------------------------------------

      use setMGBC_interface

      use mg_Xfer, ONLY: fillArray

      use pic_mass_matrix, ONLY: MM2l

      implicit none

c Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8) :: x(ntot),y(ntot)

c Local variables

      integer :: i,j,k,ig,jg,kg,isig,ijkg,nxx,nyy,nzz,ibc
      integer :: imin,imax,jmin,jmax,kmin,kmax,order
     .          ,im,ip,jm,jp,ijk,ijkim,ijkip,ijkjm,ijkjp
     .          ,ijkimjm,ijkimjp,ijkipjm,ijkipjp,ijkgl
     .          ,nxg,nyg,nzg,igl,jgl,kgl
      real(8),allocatable,dimension(:,:,:,:) :: xarr

      logical :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c Begin program
      
      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

      nxg = g_def%nxgl(igrid)
      nyg = g_def%nygl(igrid)
      nzg = g_def%nzgl(igrid)

c Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c Map vector x to array for processing

      allocate(xarr(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      xarr = 0d0

      order = 0 !BCs for mass matrix are Neumann or Periodic

      call fillArray(g_def,igrid,neq,x,xarr,bcnd,order,gpos=gpos
     .              ,is_vec=.false.)

c Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk = i + nxx*(j-1) + nxx*nyy*(k-1)

c$$$            call fromLocalToGlobalLimits(i,j,k,igl,jgl,kgl
c$$$     .                                  ,igrid,igrid,igrid)
c$$$
c$$$            ijkgl = igl + nxg*(jgl-1) + nxg*nyg*(kgl-1)
c$$$
c$$$            y(neq*(ijk-1)+1:neq*ijk) = xarr(i  ,j  ,k,:)*MM2(ijkgl,5,:)
c$$$     .                               + xarr(i+1,j  ,k,:)*MM2(ijkgl,6,:)
c$$$     .                               + xarr(i-1,j  ,k,:)*MM2(ijkgl,4,:)
c$$$     .                               + xarr(i  ,j+1,k,:)*MM2(ijkgl,8,:)
c$$$     .                               + xarr(i  ,j-1,k,:)*MM2(ijkgl,2,:)
c$$$     .                               + xarr(i-1,j-1,k,:)*MM2(ijkgl,1,:)
c$$$     .                               + xarr(i+1,j-1,k,:)*MM2(ijkgl,3,:)
c$$$     .                               + xarr(i-1,j+1,k,:)*MM2(ijkgl,7,:)
c$$$     .                               + xarr(i+1,j+1,k,:)*MM2(ijkgl,9,:)

            y(neq*(ijk-1)+1:neq*ijk) = xarr(i  ,j  ,k,:)*MM2l(ijk,5,:)
     .                               + xarr(i+1,j  ,k,:)*MM2l(ijk,6,:)
     .                               + xarr(i-1,j  ,k,:)*MM2l(ijk,4,:)
     .                               + xarr(i  ,j+1,k,:)*MM2l(ijk,8,:)
     .                               + xarr(i  ,j-1,k,:)*MM2l(ijk,2,:)
     .                               + xarr(i-1,j-1,k,:)*MM2l(ijk,1,:)
     .                               + xarr(i+1,j-1,k,:)*MM2l(ijk,3,:)
     .                               + xarr(i-1,j+1,k,:)*MM2l(ijk,7,:)
     .                               + xarr(i+1,j+1,k,:)*MM2l(ijk,9,:)
          enddo
        enddo
      enddo

c End program

      deallocate(xarr)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine massmat_mtvc



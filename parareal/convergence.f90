
      program parareal_convergence_test

      use oned_int

      use var_io

      implicit none

      !energy_old : energy from previous iteration, energy is in third column
      !energy_new : energy from this iteration, energy is in third column
      !conv_result: file for result of convergence test

      character(len = 1024) :: file_diag_old, file_diag_new, file_conv_result

      integer :: uold,unew,uout,ierr,nr=1000 !nr=number of rows in file
      integer :: i,j,ic,nnew,nold,ndiag,nc,nn

      real(8),dimension(:,:),allocatable::diag_new,diag_old
      real(8),dimension(:),allocatable :: dum
      real(8) :: err_tol,err_diag,L1,Linf

      logical :: debug=.true.

      !namelist
      namelist /conv_tst/ file_diag_old,file_diag_new,file_conv_result,err_tol

      !Read configuration
      open(unit=uinput,file='conv.in',status='old')
      read(uinput,conv_tst,iostat=ierr)
      close(uinput)

      if (debug) then
         print*, 'file from prev iteration = ', trim(file_diag_old)
         print*, 'file from this iteration = ', trim(file_diag_new)
         print*, 'file for convergence result = ', trim(file_conv_result)
         write (*,*)
      endif

      !Open files
      uold = 111
      unew = 222
      uout = 333

      open(uold,file= trim(file_diag_old),status='unknown')
      open(unew,file= trim(file_diag_new),status='unknown')
      open(uout,file= trim(file_conv_result),status='unknown')

      !Read number of diagnostics
      read(uold,*) ndiag

      if (debug) write (*,*) 'Number of diagnostics=',ndiag

      !Allocate arrays
      nc = ndiag+1
      allocate(diag_new(nr,nc),diag_old(nr,nc),dum(nr))

      diag_new=0d0
      diag_old=0d0
      dum = 0d0

      !Read OLD data
      if (debug) write (*,*)
      do  i=1,nr
         read(uold,*,err=300,end=300) diag_old(i,1:nc)
         if (debug) write (*,*) diag_old(i,1:nc)
      end do
 300  nold = i-1
      if (debug) write(*,*) "number of lines in diag_old", nold

      !Read NEW data
      if (debug) write (*,*)
      read(unew,*) ndiag
      do  i=1,nr
         read(unew,*,err=301,end=301) diag_new(i,1:nc)
         if (debug) write (*,*) diag_new(i,1:nc)
      end do
 301  nnew = i-1
      if (debug) write(*,*) "number of lines in diag_new", nnew

      !Compute diagnostic (interpolate in time in case old and new data are on different time meshes)
      if (nnew == 1) then
        err_diag = 0d0
        dum = 0d0
        nn = 1

        do ic=2,nc
           dum(1:nn) = abs((diag_old(1:nn,ic)-dum(1:nn))/diag_old(1:nn,ic))

           L1   = sum   (dum(1:nn))/dble(nn)
           Linf = maxval(dum(1:nn))
           if (debug) then
              write (*,*)
              write (*,'(a,i4)') 'Error diagnostic ',ic-1
              write (*,*) 'L-inf=',Linf
              write (*,*) 'L-1  =',L1
              write (*,*)
           endif

           err_diag = max(err_diag,L1)  !Use L-1 norm for now
        enddo

      elseif (nnew > 1) then

        err_diag = 0d0
        dum = 0d0
        nn = min(nnew,nold)
        do ic=2,nc
           if (nnew > nold) then
              call IntDriver1d (nnew,diag_new(1:nnew,1),diag_new(1:nnew,ic) &
                               ,nold,diag_old(1:nold,1),dum     (1:nold)  ,2)

              dum(1:nold) = abs((diag_old(1:nold,ic)-dum(1:nold))/diag_old(1:nold,ic))
           else
              call IntDriver1d (nold,diag_old(1:nold,1),diag_old(1:nold,ic) &
                               ,nnew,diag_new(1:nnew,1),dum     (1:nnew)  ,2)            

              dum(1:nnew) = abs((diag_new(1:nnew,ic)-dum(1:nnew))/diag_new(1:nnew,ic))
           endif

           L1   = sum   (dum(1:nn))/dble(nn)
           Linf = maxval(dum(1:nn))
           if (debug) then
              write (*,*)
              write (*,'(a,i4)') 'Error diagnostic ',ic-1
              write (*,*) 'L-inf=',Linf
              write (*,*) 'L-1  =',L1
              write (*,*)
           endif

           err_diag = max(err_diag,L1)  !Use L-1 norm for now
        enddo
      endif

      if (err_diag < err_tol) then
         write(uout,"(i4)") 0
      else
         write(uout,"(i4)") 1
      endif
      write(uout,"(ES18.8)") err_diag

      !Close files
      close(unew)
      close(uold)
      close(uout)

      !Deallocate memory
      deallocate(diag_old,diag_new,dum)

      end program parareal_convergence_test


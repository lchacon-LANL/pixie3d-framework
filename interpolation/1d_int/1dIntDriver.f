c IntDriver1d
c#######################################################################
      subroutine IntDriver1d (nx,x,vec,nv,x1,vec1,order)
      implicit none       !For safe Fortran
c***********************************************************************
c     Interpolates vector vec along coordinate x to vector vec1 along x1,
c     with interpolation order "order". Currently, we consider first
c     order (order=1), second order (2) and cubic splines (3). 
c
c     This routine requires linking with cubic spline routines dpchsp
c     (spline preparation) and dpchfe (interpolant evaluation).
c***********************************************************************

c Call variables

      integer*4      nx,nv,order
      real*8         x(nx),vec(nx),x1(nv),vec1(nv)

c Local variables

      real*8         dxx
      integer*4      i,j

c Cubic splines

      integer ::     incfd
      parameter     (incfd = 1)  !Stride to take for input vectors in spline routine

      integer ::     ic(2),nwk,ierr
      real(8) ::     vc(2),f(incfd,nx),d(incfd,nx)
      real(8), allocatable, dimension(:) :: wk
      logical ::     skip

c Externals

      real*8         q_int,l_int
      external       q_int,l_int

c Begin program

      select case (order)
      case (0) !Constant interpolation

        do i = 1,nv
          call locatep (x,nx,x1(i),j)
          if (j == 0) then
            j = 1 !Extrapolation on the left
          elseif (j == nx) then
            j = nx !Extrapolation on the right
          endif
          vec1(i) = vec(j)
        enddo

      case (1) !Linear interpolation

        do i = 1,nv
          call locatep (x,nx,x1(i),j)
          if (j == 0) then
            j = 1 !Extrapolation on the left
          elseif (j == nx) then
            j = nx !Extrapolation on the right
          endif
          vec1(i) = l_int (nx,vec,x,x1(i),j)
        enddo

      case (2) !Quadratic interpolation

        do i = 1,nv
          call locatep (x,nx,x1(i),j)
          if (j == 0) then
            j = 1 !Extrapolation on the left
          elseif (j == nx) then
            j = nx !Extrapolation on the right
          endif
          vec1(i) = q_int (nx,vec,x,x1(i),j)
        enddo
        
      case (3) !Cubic splines

        ic(1) = 0  !Default boundary condition on left
        ic(2) = 0  !Default boundary condition of right
        f(incfd,:) = vec(:)
        nwk = 2*nx

        allocate(wk(nwk))
        call dpchsp(ic,vc,nx,x,f,d,incfd,wk,nwk,ierr)

        if (ierr.ne.0) then
          write (*,*) 'Errors in cubic spline preparation'
          write (*,*) 'Aborting..'
          stop
        endif

        skip = .false.
        call dpchfe(nx,x,f,d,incfd,skip,nv,x1,vec1,ierr)

        if (ierr.lt.0) then
          write (*,*) 'Error',ierr,' in cubic spline interpolation'
          write (*,*) 'Aborting..'
          stop
c$$$        elseif (ierr.gt.0) then
c$$$          write (*,*) 'Warning: extrapolation in cubic spline in'
c$$$     .                ,ierr,' points'
        endif

        deallocate(wk)

      case default

        write (*,*) 'Order of interpolation not available'
        write (*,*) 'Aborting...'
        stop

      end select

c End program

      return
      end

c locatep
c#######################################################################
      subroutine locatep (xx,n,x,j)
      implicit      none                          ! for safe FORTRAN

c-----------------------------------------------------------------------
c     Given an array xx(1:n) and given a value x, returns a value j
c     such that x is between xx(j) and xx(j+1). xx(1:n) must be 
c     monotonic, either decreasing of increasing. j = 0 or j = n is
c     returned to indicate that x is out of range.
c-----------------------------------------------------------------------

c Variables in call

      integer*4     j,n
      real*8        x,xx(n)
 
c Local variables

      integer*4     jl,jm,ju

c Begin program

      jl = 0
      ju = n + 1

      do j = 1,n
        if (ju-jl.le.1) exit
        jm = (ju+jl)/2
        if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
          jl = jm
        else
          ju = jm
        endif
      enddo

      if (x.eq.xx(1)) then 
        j = 1
      elseif (x.eq.xx(n)) then
        j = n - 1
      else
        j = jl
      endif

c End

      return
      end

c q_int
c#######################################################################
      real*8  function q_int (nx,vec,x,x1,j)
      implicit      none                          ! for safe FORTRAN
c-----------------------------------------------------------------------
c     Interpolate vec(x) at x1 using the nodes (j-1),j,(j+1).
c-----------------------------------------------------------------------

c Variables in call

      integer*4     nx,j
      real*8        x1,x(nx),vec(nx)
 
c Local variables

      integer*4     jm,jp,jj

c Begin program

      if (j == nx) j = nx -1
      if (j == 1 ) j = 2

      jm = j-1
      jj = j
      jp = j+1

      q_int = 
     .   vec(j)* (x1 -x(jm))*(x1 - x(jp))/(x(jj)-x(jm))/(x(jj)-x(jp))
     . + vec(jp)*(x1 -x(jm))*(x1 - x(jj))/(x(jp)-x(jm))/(x(jp)-x(jj))
     . + vec(jm)*(x1 -x(jj))*(x1 - x(jp))/(x(jm)-x(jj))/(x(jm)-x(jp))

c End

      return
      end

c l_int
c#######################################################################
      real*8  function l_int (nx,vec,x,x1,j)
      implicit      none                          ! for safe FORTRAN
c-----------------------------------------------------------------------
c     Interpolate vec(x) at x1 using the nodes j,(j+1).
c-----------------------------------------------------------------------

c Variables in call

      integer*4     nx,j
      real*8        x1,x(nx),vec(nx)
 
c Local variables

      integer*4     jp,jj

c Begin program

      if (j == nx) j = nx - 1  !Extrapolation on the right

      jj = j
      jp = j+1

      l_int = vec(j) *(x1 - x(jp))/(x(jj)-x(jp))
     .      + vec(jp)*(x1 - x(jj))/(x(jp)-x(jj))

c End

      return
      end


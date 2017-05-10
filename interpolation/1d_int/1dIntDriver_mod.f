c module oned_int
c ######################################################################
      module oned_int

        use io

        use math

cc      private :: l_int,q_int

        INTERFACE quad_int
          module procedure quad_int_0,quad_int_1,quad_int_2
        end INTERFACE

      contains

c     IntDriver1d
c     ##################################################################
      subroutine IntDriver1d (nx,x,vec,nv,x1,vec1,order,deriv)
      implicit none       !For safe Fortran
c     ******************************************************************
c     Interpolates vector vec along coordinate x to vector vec1 along x1,
c     with interpolation order "order". Currently, we consider zeroth
c     order (order=0), first order (order=1), second order (2) and cubic
c     splines (3). The derivative of order "deriv" is returned.
c
c     This routine requires linking with cubic spline routines dpchsp
c     (spline preparation) and dpchfe (interpolant evaluation), and
c     other SLATEC routines.
c     ******************************************************************

c     Call variables

      integer :: nx,nv,order
      real(8) :: x(nx),vec(nx),x1(nv),vec1(nv)
      integer,optional :: deriv

c     Local variables

      real(8) :: dxx
      integer :: i,j,derv

      character(80) :: messg

c     Begin program

      if (PRESENT(deriv)) then
        derv = deriv
      else
        derv = 0
      endif

      select case (order)
      case (0) !Constant interpolation

        if (derv == 0) then
          do i = 1,nv
            call locatep (x,nx,x1(i),j)
            if (j == 0) then
              j = 1 !Extrapolation on the left
            endif
            vec1(i) = vec(j)
          enddo
        else
          vec1 = 0d0
        endif

      case (1) !Linear interpolation

        do i = 1,nv
          vec1(i) = lin_int(nx,vec,x,x1(i),derv)
        enddo

      case (2) !Quadratic interpolation

        do i = 1,nv
          vec1(i) = q_int(nx,vec,x,x1(i),derv)
        enddo
        
      case (3) !Cubic splines

        if (   maxval(x1) > maxval(x)
     .     .or.minval(x1) < minval(x)) then  !Extrapolation needed
          if (derv > 1) then
            messg = 'Cannot compute this order derivative'
            call sstop(0,'IntDriver1d',messg)
          else
            call cubic_int
          endif
        else
          call slatec
        endif

      case default

        if (   maxval(x1) > maxval(x)
     .     .or.minval(x1) < minval(x)) then  !Extrapolation needed
          messg = 'Cannot extrapolate at this order'
          call sstop(0,'IntDriver1d',messg)
        else
          call slatec
        endif

      end select

c     End program

      contains

c     slatec
c     #################################################################
      subroutine slatec

c     -----------------------------------------------------------------
c     Slatec arbitrary-order spline interpolation
c     -----------------------------------------------------------------

      implicit none

c     Local variables

      integer ::  kx,dim,flg,inbv
      real(8), dimension(:),allocatable :: tx,work,q
      real(8), dimension(:),allocatable :: bcoef

      real(8) :: dbvalu
      external   dbvalu

c     Begin program

      flg = 0

      kx = min(order+1,nx-1)
      dim = nx + 2*kx*(nx+1)

      allocate(tx(nx+kx))
      allocate(q((2*kx-1)*nx))
      allocate(work(dim))
      allocate(bcoef(nx))

      inbv=1
      call dbknot(x,nx,kx,tx)
      call dbintk(x,vec,tx,nx,kx,bcoef,q,work)
      deallocate(q)

      do i=1,nv
        vec1(i) = dbvalu(tx,bcoef,nx,kx,derv,x1(i),inbv,work)
      enddo

      deallocate(tx,work,bcoef)

      end subroutine slatec

c     cubic_int
c     #################################################################
      subroutine cubic_int

c     -----------------------------------------------------------------
c     Cubic spline interpolation (can extrapolate)
c     -----------------------------------------------------------------

      implicit none

c     Local variables

      integer,parameter :: incfd = 1  !Stride to take for input vectors in spline routine

      integer    :: ic(2),nwk,ierr
      real(8)    :: vc(2)!,f(incfd,nx),d(incfd,nx),dummy(nv)
      real(8), allocatable, dimension(:) :: wk,dummy
      real(8), allocatable, dimension(:,:) :: f,d
      logical    :: skip

c     Begin program

      nwk = 2*nx
      allocate(wk(nwk),f(incfd,nx),d(incfd,nx),dummy(nv))

      ic(1) = 0  !Default boundary condition on left
      ic(2) = 0  !Default boundary condition of right
      f(incfd,:) = vec(:)

      call dpchsp(ic,vc,nx,x,f,d,incfd,wk,nwk,ierr)

      if (ierr.ne.0) then
        messg = 'Errors in cubic spline preparation'
        call sstop(0,'cubic_int',messg)
      endif

      skip = .false.
      if (derv == 0) then
        call dpchfe(nx,x,f,d,incfd,skip,nv,x1,vec1,ierr)
      elseif (derv == 1) then
        call dpchfd(nx,x,f,d,incfd,skip,nv,x1,dummy,vec1,ierr)
      else
        messg = 'Cannot provide derivatives of this order'
        call sstop(0,'cubic_int',messg)
      endif

      if (ierr.lt.0) then
        messg = 'Error'//int2char(ierr)
     .         //' in cubic spline interpolation'
        call sstop(0,'cubic_int',messg)
      endif

      deallocate(wk,f,d,dummy)

      end subroutine cubic_int

      end subroutine IntDriver1d

c     locatep
c     ##################################################################
      subroutine locatep (xx,n,x,j)
      implicit      none                          ! for safe FORTRAN

c     ------------------------------------------------------------------
c     Given an array xx(1:n) and given a value x, returns a value j
c     such that x is between xx(j) and xx(j+1). xx(1:n) must be 
c     monotonic, either decreasing of increasing. j = 0 or j = n is
c     returned to indicate that x is out of range.
c     ------------------------------------------------------------------

c     Variables in call

      integer    :: j,n
      real(8)    :: x,xx(n)
 
c     Local variables

      integer    :: jl,jm,ju

c     Begin program

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

c     End

      end subroutine locatep

c     q_int
c     ##################################################################
      function q_int (nx,vec,x,x1,derv)
      implicit      none                          ! for safe FORTRAN
c     ------------------------------------------------------------------
c     Quadratically interpolate vec(x) or its derivatives at x1 using
c     the nodes (j-1),j,(j+1).
c     ------------------------------------------------------------------

c     Variables in call

      integer    :: nx,derv
      real(8)    :: x1,x(nx),vec(nx),q_int
 
c     Local variables

      integer    :: jm,jp,jj,j
      real(8)    :: den1,den2,den3

c     Begin program

      call locatep (x,nx,x1,j)

      j = max(2,min(j,nx-1))  !Extrapolation on both sides

      jm = j-1
      jj = j
      jp = j+1

      den1 = 1./((x(jj)-x(jm))*(x(jj)-x(jp)))
      den2 = 1./((x(jp)-x(jm))*(x(jp)-x(jj)))
      den3 = 1./((x(jm)-x(jj))*(x(jm)-x(jp)))

      select case(derv)
      case(-1)                  !Integral
        q_int =   vec(j )*den1*(x1**3/3
     .                        -(x(jm)+x(jp))*x1**2/2
     .                         +x(jm)*x(jp) *x1)
     .          + vec(jp)*den2*(x1**3/3
     .                        -(x(jm)+x(jj))*x1**2/2
     .                         +x(jm)*x(jj) *x1)
     .          + vec(jm)*den3*(x1**3/3
     .                        -(x(jp)+x(jj))*x1**2/2
     .                         +x(jp)*x(jj) *x1)
      case(0)
        q_int =   vec(j )*den1*(x1 -x(jm))*(x1 - x(jp))
     .          + vec(jp)*den2*(x1 -x(jm))*(x1 - x(jj))
     .          + vec(jm)*den3*(x1 -x(jj))*(x1 - x(jp))
      case(1)
        q_int =   vec(j )*den1*((x1 -x(jm)) + (x1 -x(jp)))
     .          + vec(jp)*den2*((x1 -x(jm)) + (x1 -x(jj)))
     .          + vec(jm)*den3*((x1 -x(jj)) + (x1 -x(jp)))
      case(2)
        q_int = 2*( vec(j )*den1
     .            + vec(jp)*den2
     .            + vec(jm)*den3)
      case default
        q_int = 0d0
      end select

c     End

      end function q_int

c     lin_int
c     ##################################################################
      function lin_int (nx,vec,x,x1,derv)
      implicit      none                          ! for safe FORTRAN
c     ------------------------------------------------------------------
c     Linearly interpolate vec(x) or its derivatives at x1 using the
c     nodes j,(j+1).
c     ------------------------------------------------------------------

c     Variables in call

      integer    :: nx,derv
      real(8)    :: x1,x(nx),vec(nx),lin_int
 
c     Local variables

      integer    :: jp,jj,j

c     Begin program

      call locatep (x,nx,x1,j)
      
      if (j == 0) then
        j = 1                   !Extrapolation on the left
      elseif (j == nx) then
        j = nx -1               !Extrapolation on the right
      endif

cc      if (j == nx) j = nx - 1  !Extrapolation on the right

      jj = j
      jp = j+1

      select case(derv)
      case(-1)                  !Integral
        lin_int = vec(j) *(x1 - x(jp))**2/(x(jj)-x(jp))*0.5
     .          + vec(jp)*(x1 - x(jj))**2/(x(jp)-x(jj))*0.5
      case(0)
        lin_int = vec(j) *(x1 - x(jp))/(x(jj)-x(jp))
     .          + vec(jp)*(x1 - x(jj))/(x(jp)-x(jj))
      case(1)
        lin_int = vec(j) /(x(jj)-x(jp))
     .          + vec(jp)/(x(jp)-x(jj))
      case(2)
        lin_int = 0d0
      end select

c     End

      end function lin_int

c     quad_int_0
c     #################################################################
      function quad_int_0(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: order
      real(8) :: x0,x1,x2,x3,y0,y1,y2,y3,x,y

c     Local variables

      character(80) :: messg

c     Begin program

      select case (order)
      case (3)
        y = y0*(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))
     .     +y1*(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))
     .     +y2*(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))
     .     +y3*(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2))
      case (2)
        y = y0*(x-x1)*(x-x2)/((x0-x1)*(x0-x2))
     .     +y1*(x-x0)*(x-x2)/((x1-x0)*(x1-x2))
     .     +y2*(x-x0)*(x-x1)/((x2-x0)*(x2-x1))
      case (1)
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      case (0)
        y = y0
      case default
        messg = 'Order of interpolation not implemented'
        call sstop(0,'quad_int_0',messg)
      end select

c     End program

      end function quad_int_0

c     quad_int_1
c     #################################################################
      function quad_int_1(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: order
      real(8) :: x0,x1,x2,x3,x
     .          ,y0(:),y1(:),y2(:),y3(:),y(size(y0))

c     Local variables

      character(80) :: messg

c     Begin program

      select case (order)
      case (3)
        y = y0*(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))
     .     +y1*(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))
     .     +y2*(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))
     .     +y3*(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2))
      case (2)
        y = y0*(x-x1)*(x-x2)/((x0-x1)*(x0-x2))
     .     +y1*(x-x0)*(x-x2)/((x1-x0)*(x1-x2))
     .     +y2*(x-x0)*(x-x1)/((x2-x0)*(x2-x1))
      case (1)
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      case (0)
        y = y0
      case default
        messg = 'Order of interpolation not implemented'
        call sstop(0,'quad_int_1',messg)
      end select

c     End program

      end function quad_int_1

c     quad_int_2
c     #################################################################
      function quad_int_2(x0,x1,x2,x3,y0,y1,y2,y3,x,order) result(y)
c     -----------------------------------------------------------------
c     Interpolation (extrapolation) routine, up to cubic order.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: order
      real(8) :: x0,x1,x2,x3,x
     .          ,y0(:,:),y1(:,:),y2(:,:),y3(:,:)
     .          ,y(size(y0,1),size(y0,2))

c     Local variables

      character(80) :: messg

c     Begin program

      select case (order)
      case (3)
        y = y0*(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))
     .     +y1*(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))
     .     +y2*(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))
     .     +y3*(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2))
      case (2)
        y = y0*(x-x1)*(x-x2)/((x0-x1)*(x0-x2))
     .     +y1*(x-x0)*(x-x2)/((x1-x0)*(x1-x2))
     .     +y2*(x-x0)*(x-x1)/((x2-x0)*(x2-x1))
      case (1)
        y = y0*(x-x1)/(x0-x1)
     .     +y1*(x-x0)/(x1-x0)
      case (0)
        y = y0
      case default
        messg = 'Order of interpolation not implemented'
        call sstop(0,'quad_int_2',messg)
      end select

c     End program

      end function quad_int_2

      end module oned_int


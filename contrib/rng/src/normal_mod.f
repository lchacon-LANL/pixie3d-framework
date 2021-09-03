      MODULE MV_NORMAL

      USE hammersley_seq

      PRIVATE    ! Make everything in module private by default 
      PUBLIC :: SETUP_MV_NORMAL,dinvnorm,multinormal_sample_rnd
     .         ,r8po_fa
cc     .         ,multinormal_sample_rnd_gr8

      !Input parameters
      INTEGER :: m                 ! Number of dimensions
     .          ,n                 ! Number of points
     .          ,seed              ! Seed for the rng
      REAL(8), pointer, DIMENSION(:)   :: mu=>NULL() ! Mean vector
      REAL(8), pointer, DIMENSION(:,:) :: a =>NULL() ! variance/covariance matrix
                           
      
      CONTAINS

      SUBROUTINE SETUP_MV_NORMAL(M_0,N_0,SEED_0,MU_0,A_0,INFO)

      IMPLICIT NONE
      !CALL VARIABLES
      INTEGER, INTENT(IN) :: M_0,N_0,SEED_0,INFO
      REAL(8), target,INTENT(IN) :: MU_0(:),A_0(:,:)

      M = M_0
      N = N_0
      SEED = SEED_0
      MU => MU_0
      A  => A_0

      if (INFO /= 0) then
        write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
        write ( *, '(a,i8)' ) '  Number of points N = ', n
        write ( *, '(a,i12)' ) '  Seed SEED = ', seed
        call r8vec_print ( m, mu, '  Mean vector MU:' )
        call r8mat_print ( m, m, a, '  Variance/covariance matrix A:' )
      endif
      
      END SUBROUTINE SETUP_MV_NORMAL


      subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ) a(n)
      integer ( kind = 4 ) i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
         write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
      end do

      return
      end subroutine r8vec_print

      subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      real ( kind = 8 ) a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end subroutine r8mat_print

      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
      implicit none

      integer ( kind = 4 ), parameter :: incx = 5
      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      real ( kind = 8 ) a(m,n)
      character ( len = 14 ) ctemp(incx)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) i2hi
      integer ( kind = 4 ) i2lo
      integer ( kind = 4 ) ihi
      integer ( kind = 4 ) ilo
      integer ( kind = 4 ) inc
      integer ( kind = 4 ) j
      integer ( kind = 4 ) j2
      integer ( kind = 4 ) j2hi
      integer ( kind = 4 ) j2lo
      integer ( kind = 4 ) jhi
      integer ( kind = 4 ) jlo
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

         j2hi = j2lo + incx - 1
         j2hi = min ( j2hi, n )
         j2hi = min ( j2hi, jhi )

         inc = j2hi + 1 - j2lo

         write ( *, '(a)' ) ' '

         do j = j2lo, j2hi
            j2 = j + 1 - j2lo
            write ( ctemp(j2), '(i8,6x)' ) j
         end do

         write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
         write ( *, '(a)' ) '  Row'
         write ( *, '(a)' ) ' '

         i2lo = max ( ilo, 1 )
         i2hi = min ( ihi, m )

         do i = i2lo, i2hi

            do j2 = 1, inc

               j = j2lo - 1 + j2

               if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
                  write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
               else
                  write ( ctemp(j2), '(g14.6)' ) a(i,j)
               end if

            end do

            write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

         end do

      end do

      return
      end subroutine r8mat_print_some

      subroutine multinormal_sample_rnd (x)

!*****************************************************************************80
!
!! MULTINORMAL_SAMPLE_RND samples a multivariate normal distribution.
!
!  Discussion:
!
!    The multivariate normal distribution for the M dimensional vector X
!    has the form:
!
!      pdf(X) = (2*pi*det(A))**(-M/2) * exp(-0.5*(X-MU)'*inverse(A)*(X-MU))
!
!    where MU is the mean vector, and A is a positive definite symmetric
!    matrix called the variance-covariance matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the space.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) A(M,M), the variance-covariance 
!    matrix.  A must be positive definite symmetric.
!
!    Input, real ( kind = 8 ) MU(M), the mean vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the random number seed.
!
!    Output, real ( kind = 8 ) X(M), the points.
!
      implicit none

      real ( kind = 8 ) x(m,n)
      integer ( kind = 4 ) info
      integer ( kind = 4 ) i,j,k,jj,jjj,imax,jmax,kmax
      real ( kind = 8 ) r(m,m),rn(m),sign(3)
      integer nr,jr
      logical qs,group
!     
!     Compute the upper triangular Cholesky factor R of the variance-covariance
!     matrix.
!     
      r(1:m,1:m) = a(1:m,1:m)

      call r8po_fa ( m, r, info )

      if ( info /= 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'MULTINORMAL_SAMPLE - Fatal error!'
         write ( *, '(a)' ) 
     .    '  The variance-covariance matrix is not positive definite 
     .       symmetric.'
         stop
      end if
!     
!     Get an MxN matrix of samples of the 1D normal distribution with mean 0
!     and variance 1.  
!     
      call r8vec_normal_01 ( m*n, seed, x(1:m,1:n) )

!     Compute R' * X.
!     We actually carry out this computation in the equivalent form X' * R.
      
      do j = 1, n
         x(1:m,j) = mu(1:m) + matmul ( x(1:m,j), r(1:m,1:m) )
      end do

      end subroutine multinormal_sample_rnd

      subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is 
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
      implicit none

      integer ( kind = 4 ) n

      integer ( kind = 4 ) m
      integer ( kind = 4 ), save :: made = 0
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) r(n+1)
!     real ( kind = 8 ) r8_uniform_01
      integer ( kind = 4 ), save :: saved = 0
      integer ( kind = 4 ) seed
      real ( kind = 8 ) x(n)
      integer ( kind = 4 ) x_hi_index
      integer ( kind = 4 ) x_lo_index
      real ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
      if ( n < 0 ) then
         n = made
         made = 0
         saved = 0
         y = 0.0D+00
         return
      else if ( n == 0 ) then
         return
      end if
!     
!     Record the range of X we need to fill in.
!     
      x_lo_index = 1
      x_hi_index = n
!     
!     Use up the old value, if we have it.
!     
      if ( saved == 1 ) then
         x(1) = y
         saved = 0
         x_lo_index = 2
      end if
!     
!     Maybe we don't need any more values.
!     
      if ( x_hi_index - x_lo_index + 1 == 0 ) then
!     
!     If we need just one new value, do that here to avoid null arrays.
!     
      else if ( x_hi_index - x_lo_index + 1 == 1 ) then

         r(1) = r8_uniform_01 ( seed )

         if ( r(1) == 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
            write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
            stop
         end if

         r(2) = r8_uniform_01 ( seed )

         x(x_hi_index) = 
     .    sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
         y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi
     $        * r(2) )

         saved = 1

         made = made + 2
!     
!     If we require an even number of values, that's easy.
!     
      else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

         m = ( x_hi_index - x_lo_index + 1 ) / 2

         call r8vec_uniform_01 ( 2 * m, seed, r )

         x(x_lo_index:x_hi_index-1:2) = 
     .    sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) 
     .    * cos ( 2.0D+00 * pi * r(2:2*m:2) )

         x(x_lo_index+1:x_hi_index:2) = 
     .    sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) 
     .    * sin ( 2.0D+00 * pi * r(2:2*m:2) )

         made = made + x_hi_index - x_lo_index + 1
!     
!     If we require an odd number of values, we generate an even number,
!     and handle the last pair specially, storing one in X(N), and
!     saving the other for later.
!     
      else

         x_hi_index = x_hi_index - 1

         m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

         call r8vec_uniform_01 ( 2 * m, seed, r )

         x(x_lo_index:x_hi_index-1:2) = 
     .    sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) 
     .    * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

         x(x_lo_index+1:x_hi_index:2) = 
     .    sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) 
     .    * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

         x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) 
     .    * cos ( 2.0D+00 * pi * r(2*m) )

         y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) 
     .    * sin ( 2.0D+00 * pi * r(2*m) )

         saved = 1

         made = made + x_hi_index - x_lo_index + 2

      end if

      return
      end subroutine r8vec_normal_01

      subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
      implicit none

      integer ( kind = 4 ) n,nh

      integer ( kind = 4 ) i
      integer ( kind = 4 ) k
      integer ( kind = 4 ) seed
      real ( kind = 8 ) r(n),rn

      if ( seed == 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
         write ( *, '(a)' ) '  Input value of SEED = 0.'
         stop
      end if

      do i = 1, n
c$$$      nh = n/2
c$$$      do i = 1,nh

         k = seed / 127773

         seed = 16807 * ( seed - k * 127773 ) - k * 2836

         if ( seed < 0 ) then
            seed = seed + 2147483647
         end if

cc         r(i) = real ( seed, kind = 8 ) * 4.656612875D-10
         rn = real ( seed, kind = 8 ) * 4.656612875D-10
         r(i) = (ceiling(rn*dble(n))-0.5d0)/dble(n)
cc         r(i) = (ceiling(rn*dble(nh))-0.5d0)/dble(n)
c$$$         r(n-i) = 1d0 - r(i)
cc         call random_number(rn)

c$$$         r(i) = (dble(i)-rn)/dble(n)


      end do

      return
      end subroutine r8vec_uniform_01

      function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
      implicit none

      integer ( kind = 4 ) k
      real ( kind = 8 ) r8_uniform_01
      integer ( kind = 4 ) seed

      if ( seed == 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
         write ( *, '(a)' ) '  Input value of SEED = 0.'
         stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
         seed = seed + 2147483647
      end if
!     
!     Although SEED can be represented exactly as a 32 bit integer,
!     it generally cannot be represented exactly as a 32 bit real number!
!     
      r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

      return
      end function r8_uniform_01

      subroutine r8po_fa ( n, a, info )

!*****************************************************************************80
!
!! R8PO_FA factors an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite 
!    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
!    upper triangular matrix, so it will be in R8GE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The positive definite symmetric matrix A has a Cholesky factorization
!    of the form:
!
!      A = R' * R
!
!    where R is an upper triangular matrix with positive elements on
!    its diagonal.  This routine overwrites the matrix A with its
!    factor R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the matrix in R8PO storage.
!    On output, the Cholesky factor R in R8GE storage.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!
      implicit none

      integer ( kind = 4 ) n

      real ( kind = 8 ) a(n,n)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) info
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      real ( kind = 8 ) s

      do j = 1, n

         do k = 1, j - 1
            a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k
     $           ,k)
         end do

         s = a(j,j) - sum ( a(1:j-1,j)**2 )

         if ( s <= 0.0D+00 ) then
            info = j
            return
         end if

         a(j,j) = sqrt ( s )

      end do

      info = 0
!     
!     Since the Cholesky factor is stored in R8GE format, be sure to
!     zero out the lower triangle.
!     
      do i = 1, n
         do j = 1, i-1
            a(i,j) = 0.0D+00
         end do
      end do

      return
      end subroutine r8po_fa

! normal inverse
!     ##################################################################
      real*8 function dinvnorm(p)
      real*8 p,p_low,p_high
      real*8 a1,a2,a3,a4,a5,a6
      real*8 b1,b2,b3,b4,b5
      real*8 c1,c2,c3,c4,c5,c6
      real*8 d1,d2,d3,d4
      real*8 z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
 201     q=dsqrt(-2*dlog(p))
         z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q
     .        +d4)*q+1)
      goto 204
 301     if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
 202     q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r
     .     +b4)*r+b5)*r+1)
      goto 204
 302     if((p.gt.p_high).and.(p.lt.1)) goto 203
 203        q=dsqrt(-2*dlog(1-p))
            z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)
     .           *q+d4)*q+1)
204   dinvnorm=z
      return
      end function dinvnorm

      END MODULE MV_NORMAL

c module eigen_mod
c ###################################################################
      module eigen_mod

      type :: eigen_options
        integer :: my_rank=0
        integer :: ksmax=10
        integer :: armax=10
        real(8) :: rtol=1d-4
        real(8) :: atol=1d-4
        real(8) :: step=1d0
        real(8) :: eigv
        integer :: iter_out
      end type eigen_options

      type(eigen_options),save :: eigen_conf

      contains

#if defined(arpack)
c     arnoldi_eig
c     ###############################################################
      subroutine arnoldi_eig(neq,ntot,x,out,matvec)
c     ---------------------------------------------------------------
c     Performs Arnoldi iteration to find largest real eigenvalue of
c     Jacobian. Uses "matvec" for matrix-vector products
c     (provided by user). 
c
c     Call parameters:
c       * neq: number of equations
c       * ntot: dimension of vector of unknowns.
c       * x: on input, initial guess; on output, eigenvector.
c       * out: level of output:
c              out = 0 --> No output
c              out = 1 --> Newton convergence output
c              out = 2 --> Previous plus GMRES convergence output
c              out = 3 --> Previous plus residual information output
c              out > 4 --> Previous plus preconditioner output (controlled
c                          in preconditioner routine).
c       * ierr: error flag:
c              ierr =-1 --> Initial guess is exact solution
c              ierr = 0 --> No error
c              ierr = 1 --> No Newton convergence to prescribed tolerance in
c                           prescribed number of iterations (ntit_max_rej)
c              ierr = 2 --> Newton converged, but took too many iterations
c                           (>  ntit_max_acc)
c
c     --------------------------------------------------------------
c     ARPACK documentation (from dnsimp.f example)
c     --------------------------------------------------------------
c
c     This example program is intended to illustrate the
c     simplest case of using ARPACK in considerable detail.
c     This code may be used to understand basic usage of ARPACK
c     and as a template for creating an interface to ARPACK.
c
c     This code shows how to use ARPACK to find a few eigenvalues
c     (lambda) and corresponding eigenvectors (x) for the standard
c     eigenvalue problem:
c
c                        A*x = lambda*x
c
c     where A is a n by n real nonsymmetric matrix.
c
c     The main points illustrated here are
c
c        1) How to declare sufficient memory to find NEV
c           eigenvalues of largest magnitude.  Other options
c           are available.
c
c        2) Illustration of the reverse communication interface
c           needed to utilize the top level ARPACK routine DNAUPD
c           that computes the quantities needed to construct
c           the desired eigenvalues and eigenvectors(if requested).
c
c        3) How to extract the desired eigenvalues and eigenvectors
c           using the ARPACK routine DNEUPD.
c
c     The only thing that must be supplied in order to use this
c     routine on your problem is to change the array dimensions
c     appropriately, to specify WHICH eigenvalues you want to compute
c     and to supply a matrix-vector product
c
c                         w <-  Av
c
c     in place of the call to AV( )  below.
c
c     Once usage of this routine is understood, you may wish to explore
c     the other available options to improve convergence, to solve generalized
c     problems, etc.  Look at the file ex-nonsym.doc in DOCUMENTS directory.
c     This codes implements
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is obtained from the standard central difference
c         discretization of the convection-diffusion operator 
c                 (Laplacian u) + rho*(du / dx)
c         on the unit square, with zero Dirichlet boundary condition.
c
c     ... OP = A  and  B = I.
c     ... Assume "call av (nx,x,y)" computes y = A*x
c     ... Use mode 1 of DNAUPD.
c
c\BeginLib
c
c\Routines called:
c     dnaupd  ARPACK reverse communication interface routine.
c     dneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     av      Matrix vector multiplication routine that computes A*x.
c     tv      Matrix vector multiplication routine that computes T*x, 
c             where T is a tridiagonal matrix.  It is used in routine
c             av.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: nsimp.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c     ---------------------------------------------------------------

      implicit none   !For safe fortran

c     Call variables

      integer    :: ntot,neq,out

      real(8)    :: x(ntot)

      external   :: matvec

      INTERFACE
        subroutine matvec(n,x,y)
        integer :: n
        real(8),dimension(n) :: x,y
        end subroutine
      END INTERFACE

c---------------------------------------------------------------------------
c$$$c
c$$$c     %------------------------------------------------------%
c$$$c     | Storage Declarations:                                |
c$$$c     |                                                      |
c$$$c     | The maximum dimensions for all arrays are            |
c$$$c     | set here to accommodate a problem size of            |
c$$$c     | N .le. MAXN                                          |
c$$$c     |                                                      |
c$$$c     | NEV is the number of eigenvalues requested.          |
c$$$c     |     See specifications for ARPACK usage below.       |
c$$$c     |                                                      |
c$$$c     | NCV is the largest number of basis vectors that will |
c$$$c     |     be used in the Implicitly Restarted Arnoldi      |
c$$$c     |     Process.  Work per major iteration is            |
c$$$c     |     proportional to N*NCV*NCV.                       |
c$$$c     |                                                      |
c$$$c     | You must set:                                        |
c$$$c     |                                                      |
c$$$c     | MAXN:   Maximum dimension of the A allowed.          |
c$$$c     | MAXNEV: Maximum NEV allowed.                         |
c$$$c     | MAXNCV: Maximum NCV allowed.                         |
c$$$c     %------------------------------------------------------%
c$$$c
c$$$      integer           maxn, maxnev, maxncv, ldv
c$$$      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
      integer :: ldv
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer :: iparam(11), ipntr(14)
      logical,allocatable,dimension(:)   :: select
      real(8),allocatable,dimension(:)   :: ax,resid,workd,workev,workl
      real(8),allocatable,dimension(:,:) :: d,v
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr,
     &                  j, ishfts, maxitr, mode1, nconv, ar_it
      Double precision
     &                  tol, sigmar, sigmai
      logical           first, rvec
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision
     &                  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy 
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic         abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting mnaupd = 1.                    |
c     %-------------------------------------------------%
c
      include 'debug.h'

      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1
      mnaup2 = 0
      mneigh = 0
      mneupd = 0
c
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      n = ntot
c
c     %-----------------------------------------------%
c     |                                               |
c     | Specifications for ARPACK usage are set       |
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |
c     |       computed.                               |
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization.                          |
c     |                                               |
c     |    3) This is a standard problem.             |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude.                      |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in DNAUPD for the     |
c     |       other options SM, LR, SR, LI, SI.       |
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 2 <= NCV <= MAXNCV             |
c     |                                               |
c     %-----------------------------------------------%
c
      nev   = 1     !Find only one eigenvalue
      ncv   = eigen_conf%ksmax
      bmat  = 'I'
      which = 'SR'  !EV with smallest real part

c     Allocate arrays

      ldv = n

      allocate(select(ncv),ax(n),d(ncv,3),resid(n)
     &        ,v(ldv,ncv),workd(3*n) 
     &        ,workev(3*ncv),workl(3*ncv*ncv+6*ncv))

c$$$      if ( n .gt. maxn ) then
c$$$         print *, ' ERROR with _NSIMP: N is greater than MAXN '
c$$$         go to 9000
c$$$      else if ( nev .gt. maxnev ) then
c$$$         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
c$$$         go to 9000
c$$$      else if ( ncv .gt. maxncv ) then
c$$$         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
c$$$         go to 9000
c$$$      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DNAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DNAUPD. (see usage below)                 |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DNAUPD.                                |
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     |
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID).  |
c     |                                                     |
c     | The work array WORKL is used in DNAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = 3*ncv**2+6*ncv 
      tol    = eigen_conf%rtol 
      ido    = 0
      if (sum(x*x) == 0d0) then
         info = 0
      else
         info = 1 ; resid = x !Provide initial guess
      endif
c     
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1).             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = eigen_conf%armax
      mode1 = 1
c
      iparam(1) = ishfts
c
      iparam(3) = maxitr
c
      iparam(7) = mode1
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) | 
c     %-------------------------------------------%
c
      if (my_rank == 0) then
         write (*,*) 'Computing Eigenvalues...'
      endif

      ar_it = 0

 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         if (out > 0 .and. my_rank == 0) then
           ar_it = ar_it + 1
           write (*,*) 'Arnoldi iteration...',ar_it
         endif

         call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, 
     &                 v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &                 info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- Op*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               | 
c           %-------------------------------------------%
c
cc            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
            call matvec(n,workd(ipntr(1):ipntr(1)+n-1)
     $                   ,workd(ipntr(2):ipntr(2)+n-1))
c$$$            workd(ipntr(2):ipntr(2)+n-1) =
c$$$     $           -workd(ipntr(2):ipntr(2)+n-1)
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         endif
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in DNAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation of _naupd'
         print *, ' '
         call pstop('arnoldi_eig','')
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        |                                           |
c        | The routine DNEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1,)                                   |
c        |                                           |
c        %-------------------------------------------%
c
         if (my_rank == 0) then
            write (*,*) 'Computing Eigenvectors...'
         endif

         rvec = .true.
c
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, 
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
c
c        %------------------------------------------------%
c        | The real parts of the eigenvalues are returned |
c        | in the first column of the two dimensional     |
c        | array D, and the IMAGINARY part are returned   |
c        | in the second column of D.  The corresponding  |
c        | eigenvectors are returned in the first         |
c        | NCONV (= IPARAM(5)) columns of the two         |
c        | dimensional array V if requested.  Otherwise,  |
c        | an orthogonal basis for the invariant subspace |
c        | corresponding to the eigenvalues in D is       |
c        | returned in V.                                 |
c        %------------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of DNEUPD. |
c           %------------------------------------%
c
            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
c
         else
c
            first = .true.
            nconv =  iparam(5)
            do 20 j=1, nconv
c
c              %---------------------------%
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*x - lambda*x ||   |
c              |                           |
c              | for the NCONV accurately  |
c              | computed eigenvalues and  |
c              | eigenvectors.  (IPARAM(5) |
c              | indicates how many are    |
c              | accurate to the requested |
c              | tolerance)                |
c              %---------------------------%
c
               if (d(j,2) .eq. zero)  then
c
c                 %--------------------%
c                 | Ritz value is real |
c                 %--------------------%
c
cc                  call av(nx, v(1,j), ax)
                  call matvec(n,v(1:n,j),ax)
cc                  v(1:n,j) = -v(1:n,j)

                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
c
               else if (first) then
c
c                 %------------------------%
c                 | Ritz value is complex. |
c                 | Residual of one Ritz   |
c                 | value of the conjugate |
c                 | pair is computed.      |
c                 %------------------------%
c
cc                  call av(nx, v(1,j), ax)
                  call matvec(n,v(1:n,j),ax)
cc                  v(1:n,j) = -v(1:n,j)

                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
cc                  call av(nx, v(1,j+1), ax)
                  call matvec(n,v(1:n,j+1),ax)
cc                  v(1:n,j) = -v(1:n,j+1)

                  call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
c
 20         continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call dmout(6, nconv, 3, d, ncv, -6,
     &           'Ritz values (Real, Imag) and residual residuals')
         end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit',
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
c
         print *, ' '
         print *, ' _NSIMP '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', 
     &              nconv 
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
c
      end if

c     Set output eigenvalue

      x = v(:,1)
c
c     %---------------------------%
c     | Done with program dnsimp. |
c     %---------------------------%
c
 9000 continue

      deallocate(select,ax,d,resid,v,workd,workev,workl)

      end subroutine arnoldi_eig
#endif

c     lehoucq_eig
c     ###############################################################
      subroutine lehoucq_eig(ntot,x,iguess,out,ierr,matvec,applyPC,ldot)
c     ---------------------------------------------------------------
c     Performs Lehoucq's iteration to find largest real eigenvalue of
c     Jacobian
c
c     Call parameters:
c       * ntot: dimension of vector of unknowns.
c       * x: on input, initial guess; on output, eigenvector.
c       * iguess: whether x contains initial guess (iguess = 1).
c       * out: level of output:
c              out = 0 --> No output
c              out > 0 --> Preconditioner output (controlled
c                          in preconditioner routine).
c       * ierr: on output, indicates error level:
c              ierr = 0 --> no error
c              ierr = 1 --> matrix not positive definite
c       * matvec: matrix-vector multiply routine (provided by user,
c                with interface matvec(n,x,y) to give y=A*x).
c       * applyPC: preconditioner application routine, with
c                interface applyPC(n,x,y,out-1) to give y=PC(x), and
c                out indicating debugging output level.
c       * ldot: optional dot-product routine, with interface ldot(x,y).

c     ---------------------------------------------------------------

      implicit none   !For safe fortran

c     Call variables

      integer    :: ntot,out,iguess,ierr

      real(8)    :: x(ntot)

      INTERFACE
        subroutine matvec(n,x,y)
        integer :: n
        real(8),dimension(n) :: x,y
        end subroutine
      END INTERFACE

      INTERFACE
        subroutine applyPC(n,x,y,out)
        integer :: n,out
        real(8),dimension(n) :: x,y
        end subroutine
      END INTERFACE

      optional :: ldot

      INTERFACE
        function ldot(x,y)
        real(8) :: ldot,x(:),y(:)
        end function
      END INTERFACE

c     Local variables

      real(8)    :: ax(ntot),dum(ntot),xn(ntot),mageig
     $             ,thetan,norm2,theta,h,mag0,magp,magr,magdp
     $             ,half=5d-1,magpAp
      integer    :: itr,itrsi,eigsub
      parameter (eigsub = 1)

      logical :: ext_dot

cc      real(8) :: pr(ntot,eigsub),pi(ntot,eigsub)

c     Begin program

      ierr = 0

      ext_dot = PRESENT(ldot)

      h = eigen_conf%step

c     Initialize iteration vector

      if (iguess == 0) call random_number(x)

c     Normalize iteration vector

      if (ext_dot) then
        norm2 = ldot(x,x)
      else
        norm2 = dot_product(x,x)
      endif

      x = x/sqrt(norm2)

c     Start iteration

      itr = 0
      magdp = 1d0

      do
        !Find matvec
        call matvec(ntot,x,ax)

        !Find |x|
        if (ext_dot) then
          norm2 = ldot(x,x)
        else
          norm2 = dot_product(x,x)
        endif

        !Form Rayleigh quotient
        thetan = theta
        if (ext_dot) then
          magpAp = ldot(x,ax)
        else
          magpAp = dot_product(x,ax)
        endif
        theta  = magpAp/norm2
        mageig = abs(theta-thetan)

        if (magpAp < 0d0) then
          ierr = 1
          return
        endif

        !Form eigenvalue residual (normalized to keep in manifold)
        dum = (theta*x - ax)/sqrt(magpAp)

        !Renormalize p=x
        magp  = sqrt(norm2)

        !Find norm of residual
        if (ext_dot) then
          magr = sqrt(ldot(dum,dum)/ntot)
        else
          magr = sqrt(dot_product(dum,dum)/ntot)
        endif

        !Diagnostics
        if (itr == 0) mag0 = magr

        if (eigen_conf%my_rank == 0) then
           write(*,'(a,i5,1p,a,e10.3,a,e10.3,a,e10.3,a,e10.3)')
     $                 'It=',itr,', eig=',theta
     $                ,', |p|=',magp,', |d_eig|=',mageig
     $                ,', |res|/|res0|',magr/mag0
        endif

        !Check convergence
        if (     magr/mag0 < eigen_conf%rtol
     $      .or. mageig    < eigen_conf%atol) exit

        itr = itr + 1

c$$$        !---> SEMI-IMPLICIT eigenvector iteration
c$$$        xn = x
c$$$        itrsi = 0
c$$$
c$$$        do
c$$$          itrsi = itrsi + 1
c$$$
c$$$          !Save current solution
c$$$          dum = x
c$$$
c$$$          !Form eigenvector new guess at (n+1/2)
c$$$          x = (xn-half*h*ax)/(1d0-half*h*theta) 
c$$$
c$$$          !Check convergence
c$$$          dum = x - dum
c$$$          magr = sqrt(dot(ntot,dum,dum)/ntot)
c$$$
c$$$          if (out > 2) write(*,*) 'SI it:',itrsi,'; residual =',magr
c$$$
c$$$          if (magr < 1d-10) exit
c$$$
c$$$          !Form Raleigh quotient
c$$$          magr  = dot(ntot,x,x)
c$$$          theta = dot(ntot,x,ax)/magr
c$$$
c$$$          !Find matvec
c$$$cc          call matvec(ntot,x,ax)
c$$$cc          write (*,*) 'h =',h,' theta=',theta
c$$$        enddo
c$$$
c$$$cc        if (out == 1) write (*,*) 'SI its:',itrsi,'; residual =',magr
c$$$        write (*,*) 'SI its:',itrsi,'; residual =',magr
c$$$
c$$$        x = 2*x-xn !Find x(n+1)

        !---> Preconditioner
        call applyPC(ntot,dum,ax,out-1)

        x = x + h*ax

      enddo

c     Store eigenvalues and other info

      eigen_conf%eigv = theta

      eigen_conf%iter_out = itr

c     End program

      end subroutine lehoucq_eig

      end module eigen_mod

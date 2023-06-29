c module ts_if
c ######################################################################
      module ts_if

      use ts_io

#if !defined(samrai)
      use nk, ONLY: nk_conf,matrixFreeMatVec,jfnk,evaluateNewtonResidual
     .             ,NOERR,CUTDT,anderson,kgm

      contains

c     implicit_ts
c     #################################################################
      function implicit_ts(x,nlit,lit,iout) result(ierr)

c     -----------------------------------------------------------------
c     Performs implicit advance of solution x (on input, x=xold)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr,nlit,lit,iout

      real(8) :: x(:)

c     Local variables

      integer :: nn,igrid

c     Begin program

      igrid = 1

      ierr = NOERR

      nn = size(x)

c     Newton-Krylov iteration

      nk_conf%etak_meth      = method
      nk_conf%ksmax          = maxksp
      nk_conf%gmmax          = maxitgm
      nk_conf%nwt_max_it_acc = ceiling(0.7*maxitnwt)
      nk_conf%nwt_max_it_rej = maxitnwt
      nk_conf%global_meth    = global

      nk_conf%eta0   = tolgm
      nk_conf%damp   = damp
      if (global == 2) nk_conf%pdt0   = dt
      if (atol  > 0d0) nk_conf%atol   = atol
      if (rtol  > 0d0) nk_conf%rtol   = rtol
      if (stol  > 0d0) nk_conf%stol   = stol
      if (mf_eps> 0d0) nk_conf%mf_eps = mf_eps
      nk_conf%mf_first_order = mf_first_order

      nk_conf%deferred_correction = nk_def_corr

      nk_conf%implicit_mid_point = implicit_mid_point

      nk_conf%krylov_method='fg'
cc      nk_conf%krylov_method='gm'
      nk_conf%fpa=fpa
      nk_conf%fpa_it=fpa_it
      nk_conf%fpa_left_pc=fpa_left_pc

      nk_conf%vol_wgt = vol_wgt

      if (fpa) then
        call anderson(neqd,nn,x,iguess,iout,ierr)
      elseif(.not.nk_is_gmres) then
        call jfnk    (neqd,nn,x,iguess,iout,ierr)
      else
        call kgm     (neqd,nn,x,iguess,iout,ierr)
      endif

c     Store convergence info

      lit  = nk_conf%gm_it_out
      nlit = nk_conf%nwt_it_out

c     End program

      end function implicit_ts

c     explicit_ts
c     #################################################################
      function explicit_ts(x,itgm,iout) result(ierr)

c     -----------------------------------------------------------------
c     Performs explicit advance of solution x (on input, x=xold)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr,iout,itgm

      real(8) :: x(:)

c     Local variables

      integer :: nn,igrid

c     Begin program

      igrid = 1

      vol_wgt = .false.

      ierr = NOERR

      nn = size(x)

c     Predictor step

c$$$      old_time_eval = .true. !This informs nlfunction that we are processing
c$$$                             !  n time level info
c$$$
c$$$      call evaluateNonlinearFunction(igrid,nn,x_old,fold)
c$$$
c$$$      old_time_eval = .false.

      x = x_old - dt*(fold-fsrc)*cnfactor

c     Corrector step

      call evaluateNonlinearFunction(igrid,nn,x,fold)

      x = x_old - dt*(fold-fsrc)

      itgm = 0

c     End program

      end function explicit_ts

c     updateSolution
c     ######################################################################
      function updateSolution(vnp,iout,ddt) result(ierr)

c     ----------------------------------------------------------------------
c     Performs time advance of solution from vn (input in vnp) to vnp
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr,iout
      real(8),optional :: ddt

      type(var_array),pointer :: vnp

c     Local variables

      integer :: itnw,itgm
      real(8) :: x(ntotd),dt_save

      INTERFACE
#if defined(EXT_EXPLICIT)
        function app_explicit(x,itgm,iout)
        use parameters
        real(8) :: x(ntotd)
        integer :: itgm,iout,app_explicit
        end function app_explicit
#endif
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use var_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c     Begin program

      if (PRESENT(ddt)) then
        dt_save = dt
        dt = ddt
      endif

      if (dt == 0d0) call pstop('updateSolution','dt=0!!!')

      ierr = NOERR

      itgmres = 0
      itnewt = 0

c     Define TS parameters

      if (bdf2) then
        bdfp  = cnp
        bdfn  = cn
        bdfnm = cnm
        cnf   = 0d0
      else
        bdfp  = 1d0
        bdfn  =-1d0
        bdfnm = 0d0
        cnf   = cnfactor
      endif
      one_over_dt = 1d0/dt

      call app_TS_Parameters(dt,one_over_dt,cnf)  !External

#if defined(MK)
      if (mk_grid) then
        if (relax_init_grid) then
          cnf = 0d0
          bdfp = 0d0
          bdfn = 0d0
          bdfnm = 0d0
          one_over_dt = 0d0
        else
          cnf  (neqd) = -tau/dt
          bdfp (neqd) = 0d0
          bdfn (neqd) = 0d0
          bdfnm(neqd) = 0d0
          one_over_dt(neqd) = 0d0
        endif
      endif
#endif

c     Time step update

      if (cnfactor <= 5d-1) then
        old_time_eval = .true. !Processing n time level info
        call evaluate_NLF(1,vnp,fold)
        old_time_eval = .false.

#if !defined(NO_BDF2_XOLD)
        call AXPYDerivedType(-cn/cnp,u_n,-cnm/cnp,u_nm,vnp)
#endif
        call mapStructureToVector(gv%gparams,x,vnp) !BDF2 ref

        x_old = x

        ierr = implicit_ts(x,itnw,itgm,iout)
      else
        old_time_eval = .true. !Processing n time level info
        call evaluate_NLF(1,vnp,fold)
        old_time_eval = .false.

        call mapStructureToVector(gv%gparams,x,vnp) !Old time info

        x_old = x

#if defined(EXT_EXPLICIT)
        ierr = app_explicit(x,itgm,iout)
#else
        ierr = explicit_ts (x,itgm,iout)
#endif
        itnw = 0
      endif

c     Map solution to vnp (if no error)

      if (ierr == NOERR .or. ierr == CUTDT) then
        call mapVectorToStructure(gv%gparams,vnp,x)

        call applyBC(1,vnp,gv%aux)

#if defined(MK)
        if (mk_grid) call MK_applyBC(1,vnp%array_var(neqd)%array)
#endif
      endif

c     Update counters

      itgmres = itgmres + itgm
      itnewt  = itnewt  + itnw

c     End program

      if (PRESENT(ddt)) dt = dt_save

      end function updateSolution

c     eigenSolution
c     ######################################################################
      function eigenSolution(vnp,iout) result(ierr)

c     ----------------------------------------------------------------------
c     Performs Eigen solution of Jacobian system. In call:
c       - v0: equilibrium state
c       - vnp: on input, initial guess; on output, eigenfunction
c     ----------------------------------------------------------------------

      use eigen_mod

      implicit none

c     Call variables

      integer :: ierr,iout

      type(var_array),pointer :: vnp

c     Local variables

      real(8) :: x(ntotd)

      INTERFACE
        subroutine applyPC(n,x,y,out)
        integer :: n,out
        real(8),dimension(n) :: x,y
        end subroutine
      END INTERFACE

c     Begin program

      ierr = NOERR

      call mapStructureToVector(gv%gparams,x ,vnp)

c     Setup vectors

      !Setup matrix-free vectors (xk, rk)
      allocate(nk_conf%xk(ntotd),nk_conf%rk(ntotd))

      call mapStructureToVector(gv%gparams,nk_conf%xk,u_0)

      call evaluateNewtonResidual(ntotd,nk_conf%xk,nk_conf%rk)

      !Setup preconditioner
      call setupPC(ntotd,nk_conf%xk)

      !Setup initial guess for eigenfunction
      x = x - nk_conf%xk

c     Eigensolve

cc#if defined(arpack)
cc        call arnoldi_eig(neqd,ntotd,x,iout)
cc#endif
      eigen_conf%my_rank = my_rank

      eigen_conf%ksmax   = maxksp
      eigen_conf%armax   = maxitgm
      if (rtol > 0d0) eigen_conf%rtol = rtol
      if (atol > 0d0) eigen_conf%atol = atol
      eigen_conf%step    = eig_dt

      iguess = 1
      call lehoucq_eig(ntotd,x,iguess,iout,ierr
     .                ,matrixFreeMatVec
     .                ,applyPC
     .                ,ldot=dot2)

      select case(ierr)
      case(1)
        call pstop('lehoucq_eig','Matrix not positive definite')
      end select

c     Output

      if (my_rank == 0) then
         write (*,*)
         write (*,*) '**************************************'
         write (*,'(a,1p,1e10.3,a)') ' ****** Eigenvalue is='
     $           ,1d0/dt-eigen_conf%eigv,' ******' 
         write (*,*) '**************************************'
      endif

      itgmres = eigen_conf%iter_out

c     Map eigenvector

      x = x + nk_conf%xk

c     Free memory

      deallocate(nk_conf%xk,nk_conf%rk)
 
      call mapVectorToStructure(gv%gparams,vnp,x)

c     End program

      end function eigenSolution

#endif

      end module ts_if

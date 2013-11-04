
c module ts_if
c ######################################################################
      module ts_if

      use ts_setup

      use nk, ONLY: nk_conf,matrixFreeMatVec,jfnk,evaluateNewtonResidual
     .             ,NOERR,CUTDT

      integer :: gmres_tot,newt_tot,itnewt,itgmres

      contains

c     implicit_ts
c     #################################################################
      function implicit_ts(neq,x,nlit,lit,iout) result(ierr)

c     -----------------------------------------------------------------
c     Performs implicit advance of solution x (on input, x=xold)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: neq,ierr,nlit,lit,iout

      real(8) :: x(:)

c     Local variables

      integer :: nn

c     Begin program

      ierr = NOERR

      nn = size(x)

      x_old = x

c     Evaluate nonlinear function at old time for theta scheme

      old_time_eval = .true.     !This informs nlfunction that we are processing
                                 !  n time level info

      if (.not.bdf2) call evaluateNonlinearFunction(nn,x,fold)

      old_time_eval = .false.

c     Newton-Krylov iteration

      nk_conf%etak_meth      = method
      nk_conf%ksmax          = maxksp
      nk_conf%gmmax          = maxitgm
      nk_conf%nwt_max_it_acc = 0.7*maxitnwt
      nk_conf%nwt_max_it_rej = maxitnwt
      nk_conf%global_meth    = global

      nk_conf%eta0   = tolgm
      nk_conf%damp   = damp
      nk_conf%pdt0   = dt0
      if (atol  > 0d0) nk_conf%atol   = atol
      if (rtol  > 0d0) nk_conf%rtol   = rtol
      if (stol  > 0d0) nk_conf%stol   = stol
      if (mf_eps> 0d0) nk_conf%mf_eps = mf_eps

      nk_conf%krylov_method='fg'
cc      nk_conf%krylov_method='gm'
      nk_conf%fpa=fpa
      nk_conf%fpa_it=fpa_it

      nk_conf%vol_wgt = vol_wgt

      call jfnk(neq,nn,x,iguess,iout,ierr)

c     Store convergence info

      lit  = nk_conf%gm_it_out
      nlit = nk_conf%nwt_it_out

c     End program

      end function implicit_ts

c     explicit_ts
c     #################################################################
      function explicit_ts(x) result(ierr)

c     -----------------------------------------------------------------
c     Performs explicit advance of solution x (on input, x=xold)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ierr

      real(8) :: x(:)

c     Local variables

      integer :: nn

c     Begin program

      ierr = NOERR

      nn = size(x)

      x_old = x

c     Predictor step

      old_time_eval = .true. !This informs nlfunction that we are processing
                             !  n time level info

      call evaluateNonlinearFunction(nn,x,fold)

      old_time_eval = .false.

      x = x_old - dt*(fold-fsrc)*cnfactor

c     Corrector step

      call evaluateNonlinearFunction(nn,x,fold)

      x = x_old - dt*(fold-fsrc)

c     End program

      end function explicit_ts

      end module ts_if

c module ts_base
c ######################################################################
      module ts_base

        use ts_if

        use variables 

#if defined(MK)
        use mk
#endif

#if !defined(petsc_c)
        type(var_array),pointer :: dn,dnm,dnp
        real(8) :: l_np,l_n,et1,et2

      contains

c     findGuess
c     ###############################################################
      subroutine findGuess(vn,vnp)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: vn,vnp

c     Local variables

c     Begin program

      if (itime <= 2) then
        !Previous time step
        vnp = vn
cc      elseif (itime == 2) then
cc        !Linear extrapolation
cc        l_np = dt/dt_n
cc        call AXPYDerivedType(1d0,vn,l_np,dn,vnp)
cc        dnp = dn   !Initialize recursion relation
      else
        !Quadratic extrapolation
        l_np = dt/dt_n
        et1 = l_np*(1d0+l_n*(1+l_np)/(1+l_n))
        et2 =-l_np*l_n**2*(1+l_np)/(1+l_n)
        call AXPYDerivedType(et1,dn,et2,dnm,dnp) !Recursion relation
        call AXPYDerivedType(1d0,vn,1d0,dnp,vnp) !vnp = vn + dnp
      endif

cc      write (*,*) itime,l_np,l_n,et1,et2

cc      if (itime <= 2) then
cc        !Previous time step
cc        vnp = vn
cccc      elseif (itime == 2) then
cccc        !FE predictor
cccc        dn = dt*old_f(:,1)
cccc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp)
cc      else
cc        !AB@ predictor
cc        old_f(:,2) = 0.5*dt*(3*old_f(:,1)-old_f(:,2))
cc        dn = old_f(:,2)
cc        call AXPYDerivedType(1d0,vn,-1d0,dn,vnp) !vnp = vn - dt/2*(3*f_n - f_n-1)
cc      endif

c     End program

      end subroutine findGuess

c     storeTSinfo
c     ###############################################################
      subroutine storeTSinfo(vn,vnp)

c     ---------------------------------------------------------------
c     Advance recursion relation
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(var_array),pointer :: vn,vnp

c     Local variables

c     Begin program

      l_n = l_np

      call AXPYDerivedType(1d0,vnp,-1d0,vn,dnp)

      dnm = dn
      dn  = dnp

cc      old_f(:,2) = old_f(:,1)
cc      call evaluateNonlinearFunction(vnp,old_f(:,1))

c     End program

      end subroutine storeTSinfo

c     correctTimeStep
c     ###############################################################
      subroutine correctTimeStep(itm,ierr,dt_exp)

c     ---------------------------------------------------------------
c     Correct time step
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer          :: ierr,itm
      real(8),optional :: dt_exp

c     Local variables

c     Begin program

c     Estimate local growth rate

      if (sm_flag == 1) call calculate_gammat

c     Find new time step

      call findNewTimeStep(itm,ierr,dt_exp=dt_exp)

c     End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

c     Calculation of local growth rate for CN

        integer     :: ieq

        real(8)     :: dmag1,dmag2,dpert,mag(neqd)

        real(8)     :: array(ilom:ihip,jlom:jhip,klom:khip)

c     Begin program

        do ieq=1,neqd

          array = (u_np%array_var(ieq)%array
     .            -u_0 %array_var(ieq)%array )
          array = array*array

          dpert = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)

          dpert = sqrt(dpert)

          array = (u_n%array_var(ieq)%array
     .            -u_0%array_var(ieq)%array )
          array = array*array

          dmag1 = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)

          array = (u_np%array_var(ieq)%array
     .            +u_n %array_var(ieq)%array
     .         -2.*u_0 %array_var(ieq)%array )
          array = array*array

          dmag2 = integral(1,1,1,nxl,nyl,nzl,array,average=.true.)

          if (dpert/=sqrt(dmag1).and.dmag2 /=0d0) then
            mag(ieq) = .5*dt*sqrt(dmag2)/(dpert-sqrt(dmag1))
          else
            mag(ieq) = 1e30
          endif

        enddo

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

      end subroutine correctTimeStep

#else

      contains

c     correctTimeStepPetsc
c     ######################################################################
      subroutine correctTimeStepPetsc(dn,dnh,dnp,ierr,dt_to_c)

      use parameters

      implicit none

c     Call variables

      integer     :: ierr
      real(8)     :: dn(neqd),dnh(neqd),dnp(neqd),dt_to_c

c     Local variables

c     Begin program

c     Estimate local growth rate

      if (sm_flag == 1) call calculate_gammat

c     Find new time step

      call findNewTimeStep(itime+1,ierr)

c     Update counters (only if timeStep is successful)

      itime  = itime + 1

      time   = time  + dt
      tmrst  = tmrst + dt
      nrst   = nrst  + 1

      dt_to_c = dt

c     End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

c     -----------------------------------------------------------------------
c     Calculation of local growth rate for CN
c    ------------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer     :: ieq

        real(8)     :: mag(neqd)

c     Begin program

        where (dnp /= dn .and. dnh /= 0d0)
          mag = dt*dnh/(dnp-dn)
        elsewhere
          mag = 1e30
        end where

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

      end subroutine correctTimeStepPetsc

#endif

c     findNewTimeStep
c     ####################################################################
      subroutine findNewTimeStep(itm,ierr,dt_exp)

c     --------------------------------------------------------------------
c     Find new time step
c     --------------------------------------------------------------------

      implicit none

c     Call variables

      integer    ::  ierr,itm
      real(8),optional :: dt_exp

c     Local variables

c     Begin program

      if (tn == tnm) then
         dtold = dt
         dt_n  = dt
      else
         dtold = tn-tnm
         dt_n  = tn-tnm
      endif

      if (PRESENT(dt_exp)) then
        dt = dt_exp
      else
        call adapt_dt

        call calculate_cnfactor
      endif

      alpha = 1d0 - cnfactor

c     End program

      contains

c     adapt_dt
c     #######################################################################
      subroutine adapt_dt

        real(8) ::    coef1=0.8   !Time subcycling coefficient
     .               ,coef2=1.05  !Time recovery   coefficient

        if (ierr > 0) then   !Error in time stepping

          dt = dt*coef1
          if (dt < 1d-3*dtbase) then
            call pstop('adapt_dt','Time step too small')
          endif
          if (my_rank == 0 .and. ierr == CUTDT) write (*,240)
          if (my_rank == 0) write (*,400) dt

        else                 !No error in time stepping

          select case (sm_flag)
          case(1)  !Rannacher time stepping

            if (   (itm <= sm_pass  .and.(.not.timecorr))
     .         .or.(itm <= sm_pass+1.and.      timecorr )) then
              dt = 0.5*dtbase
            elseif ((itm == sm_pass+1.and.(.not.timecorr))
     .          .or.(itm == sm_pass+2.and.      timecorr )) then
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
            endif

          case(2)  !BDF time stepping

            if (timecorr) then
              if  (itm == 2) then
                dt = 0.1*dtbase
              elseif (itm == 3) then
                dt = dtbase/coef2**15  !Quiet BDF2 start to ramp up in a few time steps
              else
                dt = min(dtbase,dt*coef2)
              endif
            else
              if  (itm == 1) then
                dt = 0.1*dtbase
                dtold = dt
              elseif (itm == 2) then
                dt = dtbase - dt
              elseif (itm == 3) then
                dt = dtbase
              else
                dt = min(dtbase,dt*coef2)
              endif
            endif

          case default   !Default time stepping
            if (    (itm == 1.and.(.not.timecorr))
     .          .or.(itm == 2.and.      timecorr )) then
              dt = dtbase
            else
              dt = min(dtbase,dt*coef2)
            endif
          end select

        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step: ',1p,e10.4)

      end subroutine adapt_dt

c     calculate_cnfactor
c     #######################################################################
      subroutine calculate_cnfactor

      if (sm_flag == 1) then      !Rannacher

        if (itm == 1.and.timecorr) then
          cnfactor = .5
        elseif (    (itm <= sm_pass+1 .and.      timecorr )
     .          .or.(itm <= sm_pass   .and.(.not.timecorr)))then
          cnfactor = .3
        elseif (    (itm >= sm_pass+1.and.(.not.timecorr))
     .          .or.(itm >= sm_pass+2.and.      timecorr )) then
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      elseif (sm_flag == 2) then  !BDF2

        cnfactor = 0d0

        bdf2 = (itm > 1 .and. (.not.timecorr) )
     .    .or. (itm > 2 .and.       timecorr  .and. ierr == NOERR)

        if (bdf2) then
          !All these are affected by factor 1/dt (done in evaluateNLF)
          cnp = (2*dt+dtold)/(dt+dtold)
          cn  = -(dt+dtold)/dtold
          cnm = -cnp-cn    !Sum should be zero identically
        else
          cnp = 1d0
          cn  =-1d0
          cnm = 0d0
        endif

      endif

      end subroutine calculate_cnfactor

      end subroutine findNewTimeStep

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

c     Begin program

      if (PRESENT(ddt)) then
        dt_save = dt
        dt = ddt
      endif

      ierr = NOERR

      call mapStructureToVector(x,vnp)

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

      if (dt == 0d0) call pstop('updateSolution','dt=0!!!')

c     Time step update

      if (cnfactor <= 5d-1) then
        ierr = implicit_ts(neqd,x,itnw,itgm,iout)
      else
        ierr = explicit_ts(x)
        itnewt  = 1
        itgmres = 0
      endif

c     Map solution to vnp (if no error)

      if (ierr == NOERR .or. ierr == CUTDT) then
        call mapVectorToStructure(vnp,x)

        call applyBC(1,vnp,it=itime)

#if defined(MK)
        if (mk_grid) call MK_applyBC(1,vnp%array_var(neqd)%array)
#endif
      endif

c     Update counters

      itgmres = itgmres + itgm
      itnewt  = itnewt  + itnw

      gmres_tot = gmres_tot + itgmres
      newt_tot  = newt_tot  + itnewt

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

      call mapStructureToVector(x ,vnp)

c     Setup vectors

      !Setup matrix-free vectors (xk, rk)
      allocate(nk_conf%xk(ntotd),nk_conf%rk(ntotd))

      call mapStructureToVector(nk_conf%xk,u_0)

      call evaluateNewtonResidual(ntotd,nk_conf%xk,nk_conf%rk)

      !Setup preconditioner
      call setupPC(ntotd,nk_conf%xk,1)

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
 
      call mapVectorToStructure(vnp,x)

c     End program

      end function eigenSolution

      end module ts_base

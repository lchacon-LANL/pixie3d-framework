c setupPC
c####################################################################
      subroutine setupPC(ntot,x)

c--------------------------------------------------------------------
c     Finds velocity and magnetic field components in all grids
c     for matrix-free preconditioning. Here, x is current Newton
c     state vector.
c--------------------------------------------------------------------

      use grid

#if defined(MK)
      use mk_pc
#else
      use variables
#endif

      use ts_setup

      implicit none

c Call variables

      integer :: ntot
      real(8) :: x(ntot)

c Local variables

      integer :: igrid,ieq

      type (var_array),pointer :: varray => null()

      INTERFACE
        subroutine setup_app_PC(igrid,varray)
        use var_setup
        integer    :: igrid
        type(var_array),pointer :: varray
        end subroutine setup_app_PC
      END INTERFACE

c Begin program

      if (pc_id) return

      igrid = 1

c Change profiling stage

#if defined(petsc)
      call switch_profile_stage(iPCSTP)
#endif

c Copy solution state info

      call mapVectorToStructure(gv%gparams,varray,x)

c Setup MK preconditioner

#if defined(MK)
      if (mk_grid) then
        call MK_setupPC(igrid,varray%array_var(neqd)%array)
        varray%nvar = varray%nvar-1  !Substract MK equation
      endif
#endif

c Setup application PC

      call setup_app_PC(igrid,varray)

c Free memory

      call deallocateDerivedType(varray)

c End PC setup profiling

#if defined(petsc)
      call switch_profile_stage(iNL)
#endif

c End program

      end subroutine setupPC

c applyPC
c###################################################################
      subroutine applyPC(ntot,y,x,iout)

c-------------------------------------------------------------------
c     This subroutine solves P x = y for the vector x. 
c     The parameter iout controls output level:
c       * iout <= 0 --> No output
c       * iout >  0 --> Level of output design by user.
c-------------------------------------------------------------------

      use var_io

      use timeStepping

#if defined(MK)
      use mk_pc
#else
      use nk
#endif

      use mg_io

      use profile_petsc

      implicit none

c Call variables

      integer :: ntot,iout
      real(8) :: x(ntot),y(ntot)

c Local variables

      real(8) :: xxx(ntot/neqd,neqd),yyy(ntot/neqd,neqd),dvol

      integer :: ii,iii,igrid,nn,neq,ieq,nx,ny,nz,i,j,k,ig,jg,kg

      logical :: pc_diag

c Begin program

      igrid = 1

      nx = gv%gparams%nxv(igrid)
      ny = gv%gparams%nyv(igrid)
      nz = gv%gparams%nzv(igrid)

      nn = nx*ny*nz

#if defined(petsc)
      !Profiling stage PC
      call switch_profile_stage(iPCSLV)
#endif

c *******************************************************************
c     Identity preconditioner
c *******************************************************************

      if (pc_id) then

        if (iout > 1 .and. my_rank == 0) then
          write (*,*) '--> Performing IDENTITY PC'
        endif

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                if (one_over_dt(ieq) == 0d0) then
                  x(iii) = y(iii)
                else
                  if (vol_wgt) then
                    x(iii) = y(iii)*dt
     .                   /gv%gparams%gmetric%grid(igrid)%dlvol(i,j,k)
                  else
                    x(iii) = y(iii)*dt
                  endif
                endif
              enddo

            enddo
          enddo
        enddo

      else

#if defined(MK)
        if (mk_grid) then
          neq = neqd-1   !Substract MK equation
        else
          neq = neqd
        endif
#else
       neq = neqd
#endif

        xxx = 0d0

c     Scatter residuals

        call scatterVector(yyy,y)

c     MK PC

#if defined(MK)
        if (mk_grid) call MK_applyPC(igrid,yyy,xxx(:,neqd),iout)  !Last equation is MK
#endif

c     Application PC

#if defined(MK)
        if (relax_init_grid) then !Trivial PC for spline EQU
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid
     .                       ,ig,jg,kg)

                if (vol_wgt) then
                  if (mk_nc) then
                    dvol = gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
                  else
                    dvol = gv%gparams%gmetric%grid(igrid)%dlvol(i,j,k)
                  endif
                else
                  dvol = 1d0
                endif

                xxx(ii,1:neq) = yyy(ii,1:neq)/dvol
              enddo
            enddo
          enddo
        else
          call apply_app_PC(igrid,nn,neq,yyy(:,1:neq),xxx(:,1:neq),iout)
        endif
#else
        call apply_app_PC(igrid,nn,neq,yyy(:,1:neq),xxx(:,1:neq),iout)
#endif

c     Gather solution in "x"

        call gatherVector(x,xxx)

      endif

c     PC diagnostics

      if (pc_debug) then
        if (my_rank == 0) then
          write (*,'(a)',ADVANCE='no')
     .             '  ==> Run PC diagnostics (t or f)? '
          read(*,*) pc_diag
        endif

#if defined(petsc) || defined(petsc_c)
        call MPI_Bcast(pc_diag,1,MPI_LOGICAL,0,MPI_COMM_WORLD
     .                ,mpierr)
#endif
        if (pc_diag) call PCdiagnostics
      endif

c     End PC profiling

#if defined(petsc)
      !Begin profiling NL solver
      call switch_profile_stage(iNL)
#endif

c End program

      contains

c     PCdiagnostics
c     ##############################################################
      subroutine PCdiagnostics

c     --------------------------------------------------------------
c     Performs diagnostics on PC
c     --------------------------------------------------------------

        implicit none

        real(8) :: mag,mag0,rr(ntot)

c     Begin program

c     Plot PC solution

        call vec_plot_gl(gv%gparams,neqd,-x,igrid,'sol','Solution')

c     Plot PC RHS

        call vec_plot_gl(gv%gparams,neqd,-y,igrid,'rhs','RHS')

c     Plot corresponding linear residual to PC solution

#if !defined(petsc_c)
        if (fpa) then
          mag0 = rms(x)

          if (my_rank == 0) then
            write (*,*)
            write (*,*) 'FPA RESIDUAL'
            write (*,*) 'RMS dsol=',mag0
          endif

          call vec_plot_gl(gv%gparams,neqd,x,igrid,'res','Residual')

          !Compute magnitude of individual residual eqs
          call scatterVector(yyy,x)
          do ieq=1,neqd
            mag = rms(yyy(:,ieq))
            if (my_rank == 0) write (*,*) 'RMS sol update',ieq,':',mag
          enddo

        else
          mag0 = rms(y)

          call matrixFreeMatVec(ntot,x,rr)
          y = y - rr
          mag = rms(y)

          call vec_plot_gl(gv%gparams,neqd,y,igrid,'res','Residual')

          if (my_rank == 0) then
            write (*,*)
            write (*,*) 'LINEAR RESIDUAL'
            write (*,*) 'RMS res=',mag
            write (*,*) 'Rel res=',mag/mag0
          endif

          !Compute magnitude of individual residual eqs
          call scatterVector(yyy,y)
          do ieq=1,neqd
            mag = rms(yyy(:,ieq))
            if (my_rank == 0) write (*,*) 'RMS residual eqn',ieq,':',mag
          enddo
        endif
#endif

        call pstop('PCdiag','Finished diagnostics...')

      end subroutine PCdiagnostics

      end subroutine applyPC

c killPC
c###########################################################################
      subroutine killPC

#if defined(MK)
      use mk_pc
#endif

      use ts_setup

      implicit none
      
c     Deallocate variables

      if (pc_id) return

#if defined(MK)
      if (mk_grid) call MK_killPC
#endif

      call kill_app_PC

      end subroutine killPC

      program plot3d

c##########################################################################
c     Postprocess output information from driver and creates graphics
c     in various formats.
c##########################################################################

      use parameters

      use variables

      use graphics

      use ts_setup

      use grid

      use var_io

      use acc_tests

#if defined(MK)
      use mk
#endif

      implicit none

c Local variables

      integer        :: i,j,k,ig,jg,kg
      integer        :: ierr,gierr,itime_old
      character*(40) :: command
      logical        :: exitloop,petsc_flg

c Begin program

c Initialize PETSc (parallel version)

#if defined(petsc)
      call PetscInitialize(PETSC_NULL_CHARACTER,mpierr)
#endif

c Initialize computation

      call initializeCalculation

c Initialize counters

      nplot  = 0
      tmplot = 0d0

c Time loop

      do

c     Update counters

        call TS_update_preCounters(time)

c     Read next record 

        itime_old = itime
        
        call readTimeStep(urecord,itime,time,dt,gammat,u_np,ierr)

c     Check for alternate restart files and read again if present
        
        if(switch_dumpfile(itime_old,itime,urecord)) then
          call readTimeStep(urecord,itime,time,dt,gammat,u_np,ierr)
        endif
        
c     Gather error info

#if defined(petsc)
        call MPI_Allreduce(ierr,gierr,1,MPI_INTEGER,MPI_SUM
     .                    ,MPI_COMM_WORLD,mpierr)
        exitloop = (gierr /= 0)
#else
        exitloop = (ierr  /= 0)
#endif

c     Check for error

        if (exitloop) exit

c     Update counters

        nplot  = nplot + 1

c     MK setup

#if defined(MK)
        !Setup MK grid
        if(mk_grid) then

          old_time_eval = .false.
          xnm  = xold
          xold = xnew

          call MK_setup_grid(1,u_np%array_var(neqd)%array)

cc          !Grid velocity in Cartesian geometry
cc          if (bdf2) then
cc            gvel =-(cnp*xnew+cn*xold+cnm*xnm)/dt
cc          else
cc            gvel =-(xnew-xold)/dt
cc          endif

          call MK_xfer_equ(1,u_ic)
        endif
#endif

c     Output per time step

        call output

c     Store previous record (for time derivatives in diagnostics)

        call equateDerivedType(u_n,u_np)

c     Check for tmax

        if ((tmax > 0d0 .and. time.ge.tmax))  exit

      enddo       !End of time loop

c Numerical accuracy tests

      if (temp_acc) call acc_test_time
      if (spat_acc) call acc_test_space

c Finalize graphics

      call finalizeGraphics

      call finalizeDiagnostics

      call closeGraphicsFiles

      call fortranDestroy

#if defined(petsc)
      call PetscFinalize(mpierr)
#endif

c End program

      end program

c output
c######################################################################
      subroutine output

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use ts_io

      use variables

      use graphics

      implicit none

c Call variables

c Local variables

      logical :: plot

c Interfaces

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use var_setup
        integer    :: iigx,iigy,iigz
        logical    :: init
        type(var_array),pointer :: varray
        end subroutine evaluateDiagnostics
      END INTERFACE

c Begin program

c Post-process solution for diagnostics and dump plots

      if (.not.skip_diag)
     .       call evaluateDiagnostics(u_np,1,1,1,.false.)

      plot = TS_chk_plot_dump(.true.)

      if (plot) then

        if (skip_diag) 
     .       call evaluateDiagnostics(u_np,1,1,1,.false.)

        !Dump time step plots
        call dumpTimeStepPlots

      endif

c Information output

      if (my_rank == 0) then

        if (itime.eq.0) then

          write (*,10) nxd,nyd,nzd

          write (*,*) 
          write (*,100)

        endif

        if (plot) then
          write (*,120) itime,time,' --> dumped plots'
        else
          write (*,110) itime,time
        endif

      endif

c End program

 10   format (/,'  Grid mesh:................',i4,'x',i4,'x',i4)
 100  format ('   itime      time')
 110  format (i7,2x,1p,e10.2)
 120  format (i7,2x,1p,e10.2,a)

      end subroutine output

c initializeCalculation
c######################################################################
      subroutine initializeCalculation

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use ts_setup

      use var_io

      use graphics

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr,igx,igy,igz
      logical :: petsc_flg

c Interfaces

      INTERFACE
        subroutine initializeDiagnostics(varray,iigx,iigy,iigz)
        use var_setup
        integer :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine initializeDiagnostics
      END INTERFACE

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use var_setup
        integer :: iigx,iigy,iigz
        logical :: init
        type(var_array),pointer :: varray
        end subroutine evaluateDiagnostics
      END INTERFACE

c Begin program

      igx = 1 ; igy = 1 ; igz = 1

      call TS_init_counters

c Read user initializations

      call readGraphicsInput

#if defined(MK)
      if (mk_grid) neqd = neqd+1  !Add MK equation
#endif

c Runtime initializations

#if defined(petsc)
      !Read PETSc runtime options
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-npx",npx,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-npy",npy,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-npz",npz,petsc_flg,ierr)

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-tmax",tmax,petsc_flg,ierr)

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-dplot",dplot,petsc_flg,ierr)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-ndplot",ndplot,petsc_flg,ierr)

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-numtime",numtime,petsc_flg,ierr)
      if (petsc_flg) tmax = 0d0

      call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER
     .   ,"-ndstep",ndstep,petsc_flg,ierr)
      if (petsc_flg) dstep = 0d0
#endif

c Allocate global variable

      call allocateGlobalVar(gv)

c Open graphics files

      urecord = openGraphicsFiles(recordfile)

c Read geometry (if appropriate)

cc      if (coords == 'ext') call readExtMap(gv%gparams)

c Read equilibrium u_0 (Do NOT postprocess for BC)

      call readTimeStep(urecord,itime,time,dt,gammat,u_0,ierr)

      if (ierr /= 0) call pstop('initializeCalculation',
     .                          'Unable to read equilibrium')

c Define application arrays (external)

      call allocateApplicationVariables     !External

c Initialize MK grid

#if defined(MK)
      if (mk_grid) then
        call MK_alloc
        call MK_initGrid(igx,bcs=bcond)
cc        call MK_initGrid(1)
      endif
#endif

c Set plotting range

      call setGraphicsRange

c Read perturbed solution u_ic

      call equateDerivedType(u_ic,u_0)

c$$$      call readTimeStep(urecord,itime,time,dt,gammat,u_ic,ierr)

c Initialize diagnostics and related variables

#if defined(MK)
      if (mk_grid) then
        call MK_spline_equ(igx,u_ic) !Spline equilibrium
        call MK_setup_grid(igx,u_ic%array_var(neqd)%array)
      endif
#endif

      call equateDerivedType(u_n,u_ic)

#if defined(MK)
      if (mk_grid) u_ic%nvar = u_ic%nvar-1
#endif
      call initializeDiagnostics(u_ic,igx,igy,igz)    !External

#if defined(MK)
      if (mk_grid) u_ic%nvar = u_ic%nvar+1
#endif

      call equateDerivedType(u_np,u_ic)

c Initialize graphics

      call initializeGraphics

c End program

      end subroutine initializeCalculation

c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use var_io

      use grid

      use graphics

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

c Local variables

      integer :: ierr,igr

c Begin program

      call deallocateStructures

      call deallocateDerivedType(u_ic)

      call deallocateApplicationVariables   !External

#if defined(MK)
      if (mk_grid) call MK_dealloc
#endif

      call deallocateGlobalVar(gv)

c End program

      end subroutine fortranDestroy

      program plot3d

c##########################################################################
c     Postprocess output information from 3DMHD and creates graphics
c     in various formats.
c##########################################################################

      use parameters

      use variables

      use graphics

      use ts_setup

      use grid

      use iosetup

      use acc_tests

#if defined(MK)
      use mk
#endif

      implicit none

c Local variables

      integer        :: i,j,k,ig,jg,kg
      integer        :: ierr,gierr,nplot,igx,igy,igz
      character*(40) :: command
      logical        :: exitloop

c Interfaces

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use variable_setup
        integer    :: iigx,iigy,iigz
        logical    :: init
        type(var_array),pointer :: varray
        end subroutine evaluateDiagnostics
      END INTERFACE

c Begin program

      igx = 1
      igy = 1
      igz = 1

c Initialize 

      call initializeCalculation(igx,igy,igz)

c Initialize counters

      nplot  = 0
      tmplot = 0d0

c Time loop

      do

c diag ***** Redefine perturbations for preconditioner debugging
cc        write (*,*) 'Redefined perts for PC debugging in 3dplot'
cc        call equateDerivedType(u_ic,u_n)
cc        call applyBC(u_ic,igx,igy,igz)
c diag *****

c     Read next record

        call readTimeStep(u_np,ierr)

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

c     Time level plots (xdraw)

        plot = .false.

        if (    nplot        == ndplot
     .     .or.(time-tmplot) >= dplot-0.9*dt
     .     .or. time         == 0d0) then

#if defined(MK)
          !Setup MK grid
          if(mk_grid) then

            old_time_eval = .false.
            xnm  = xold
            xold = xnew

            call MK_setup_grid(igx,u_np%array_var(neqd)%array)

cc            !Grid velocity in Cartesian geometry
cc            if (bdf2) then
cc              gvel =-(cnp*xnew+cn*xold+cnm*xnm)/dt
cc            else
cc              gvel =-(xnew-xold)/dt
cc            endif

            call MK_xfer_equ(igx,u_ic)
          endif
#endif
          !Post-process solution for diagnostics
          call evaluateDiagnostics(u_np,igx,igy,igz,.false.)

          !Dump time step plots
          plot  = .true.
          nplot = 0
          if (itime.gt.0) tmplot = time
          call dumpTimeStepPlots

          !Store previous record (for time derivatives in diagnostics)
          call equateDerivedType(u_n,u_np)

        endif

c     Output per time step

        call output(plot)

c     Check for tmax

cc        if ((tmax > 0d0 .and. time.ge.tmax*0.99999)
cc     .      .or. (tmax < 0d0 .and. itime == numtime))  exit
        if ((tmax > 0d0 .and. time.ge.tmax*0.99999))  exit

      enddo       !End of time loop

c Numerical accuracy tests

      if (temp_acc) call acc_test_time
      if (spat_acc) call acc_test_space

c Finalize graphics

      call finalizeGraphics

      call finalizeDiagnostics

#if defined(adios)
      call closeRecordFileAdiosForRead(adios_err)
      call adios_finalize(my_rank,adios_err)
#else
      call closeRecordFiles
#endif

      call fortranDestroy

#if defined(petsc)
      call PetscFinalize(mpierr)
#endif

c End program

      end

c output
c######################################################################
      subroutine output(plot)

c----------------------------------------------------------------------
c     Writes program output to standard output
c----------------------------------------------------------------------

      use ts_setup

      use parameters

      use grid

      implicit none

c Call variables

      logical :: plot

c Local variables

      integer     :: ngrd(3)

c Begin program

      if (my_rank == 0) then

        ngrd =(/grid_params%ngrdx,grid_params%ngrdy,grid_params%ngrdz/)

        if (itime.eq.0) then

          write (*,10) nxd,nyd,nzd,ngrd

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

 10   format (/,'  Grid mesh:................',i4,'x',i4,'x',i4
     .        /,'  Number of grid levels:....',i4,',',i2,',',i2)
 100  format ('   itime      time')
 110  format (i7,2x,1p,e10.2)
 120  format (i7,2x,1p,e10.2,a)

      end subroutine output

c initializeCalculation
c######################################################################
      subroutine initializeCalculation(igx,igy,igz)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use ts_setup

      use iosetup

      use graphics

#if defined(MK)
      use mk
#endif

      implicit none

c Call variables

      integer    :: igx,igy,igz

c Local variables

      integer    :: ierr
      logical    :: petsc_flg

c Interfaces

      INTERFACE
        subroutine initializeDiagnostics(varray,iigx,iigy,iigz)
        use variable_setup
        integer          :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine initializeDiagnostics
      END INTERFACE

      INTERFACE
        subroutine evaluateDiagnostics(varray,iigx,iigy,iigz,init)
        use variable_setup
        integer    :: iigx,iigy,iigz
        logical    :: init
        type(var_array),pointer :: varray
        end subroutine evaluateDiagnostics
      END INTERFACE

c Begin program

c Initialize PETSc (parallel version)

#if defined(petsc)
      call PetscInitialize(PETSC_NULL_CHARACTER,mpierr)

      !Read PETSc runtime options
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npx',npx
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npy',npy
     .                       ,petsc_flg,ierr)
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-npz',npz
     .                       ,petsc_flg,ierr)
#endif

c Initialize ADIOS

#if defined(adios) && defined(petsc)
      call adios_init_noxml (adios_err)
      ! allocate 100MB buffer for ADIOS
      call adios_allocate_buffer (100,adios_err)
#endif

c Read user initializations

      call readGraphicsInput

#if defined(MK)
      if (mk_grid) neqd = neqd+1  !Add MK equation
#endif

c Allocate global variable

      call allocateGlobalVar(gv)
cc      allocate(gv)
cc      nullify(gv%bc_grp)
cc      nullify(gv%u_np,gv%u_n,gv%u_nm,gv%u_0)
cc      nullify(gv%gparams)
cc      nullify(gv%aux)
cc
ccc Initialize MG and create grid
cc
cc      call createGrid(nxd,nyd,nzd,gv%gparams)

c Open graphics files (this initializes local limits, among other things)

#if defined(adios)
      call openRecordFileAdiosForRead(adios_err)
#else
      call openRecordFiles(urecord,recordfile,murecord
     .                    ,ilog,ihig,jlog,jhig,klog,khig)
#endif

ccc Re-initialize vector dimensions
cc
cc      call setVectorDimensions
cc
ccc Allocate records
cc
cc      call allocateStructures

c Read equilibrium u_0 (Do NOT postprocess for BC)

      call readTimeStep(u_0,ierr)

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

      call setGraphicsRange(igx,igy,igz)

c Initialize diagnostics and related variables

      call equateDerivedType(u_ic,u_0)

#if defined(MK)
      if (mk_grid) then
        call MK_spline_equ(igx,u_ic) !Spline equilibrium
        call MK_setup_grid(igx,u_ic%array_var(neqd)%array)
      endif
#endif

      call equateDerivedType(u_n ,u_0)

#if defined(MK)
      if (mk_grid) u_ic%nvar = u_ic%nvar-1
#endif
      call initializeDiagnostics(u_ic,igx,igy,igz)    !External

#if defined(MK)
      if (mk_grid) u_ic%nvar = u_ic%nvar+1
#endif

      call equateDerivedType(u_np,u_ic)
      call equateDerivedType(u_n ,u_ic)

c Initialize graphics

      call initializeGraphics(igx,igy,igz,bcond)

c End program

      end subroutine initializeCalculation

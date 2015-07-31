c module TS_solver_setup
c ######################################################################
      module TS_solver_setup

        use grid, ONLY: np,my_rank,pstop

        !NKS setup
        integer :: maxitnwt = 20          !Maximum number of Newton its. (if zero, maxitnwt
                                          !  is determined in code)
     .            ,maxitgm  = 15          !Maximum GMRES iterations
     .            ,iguess   = 0           !Whether preconditioner is used to give initial
                                          !  guess to GMRES (when =1)
     .
     .            ,maxksp   = 15          !Maximum krylov subspace dimension
     .            ,method   = 1           !Inexact Newton method:
                                          !  + 0: constant forcing
                                          !  + other: adaptive (Eisenstat-Walker)
     .
     .            ,global   = 0           !Newton's globalization method:
                                          !  + 0 --> no globalization
                                          !  + 1 --> linesearch backtracking 
                                          !  + 2 --> pseudo-transient

        real(8) :: tolgm=0d0              !Inexact Newton parameter                 
     .            ,rtol=0d0               !Newton relative convergence tolerance    
     .            ,atol=0d0               !Newton absolute convergence tolerance    
     .            ,stol=0d0               !Newton update convergence tolerance      
     .            ,mf_eps=1d-8            !Newtom matrix-free differencing parameter
     .            ,damp=1d0               !Damping parameter in Newton
     .            ,dt0=1d30               !Initial pseudo-transient time step

        logical :: mf_first_order=.true.  !Whether to use first or second-order MF formula

        logical :: nk_def_corr=.false. !Whether we use deferred correction or not

        !FPA configuration
        logical :: fpa=.false.            !Whether to use Carlson's Fixed-Point-Accelerator instead of JFNK
        integer :: fpa_it=5               !Number of FPA iterations

        !PC setup
        integer :: pc_iter        = 1     !Number of SI iterations in preconditioner
     .            ,mg_ores        = 0     !Restriction order for MG
     .            ,mg_oprol       = 2     !Prolongation order for MG
     .            ,mg_coarse_size = 2     !Minimum coarse grid size (mg_ratio^mg_coarse_size)
     .            ,mg_vcyc        = 1     !Maximum number of MG V-cycles
     .            ,sm_zebra_it    = 1     !Number of relaxations in ZEBRA smoothing
     .            ,sm_zebra_dir   = 0     !Preferred smoothing direction in ZEBRA smoothing (all if zero)
     .            ,sm_ncolors     = 8     !Number of colors for diagonal calc., smoothing
     .            ,sm_iter        = 5     !Number of MG smoothing passes

        real(8) :: pc_tol        = 1d-3   !PC solvers convergence tolerance
     .            ,sm_omega      = 0.7    !Damping for MG smoother (only when 'jb' is chosen)
     .            ,sm_zebra_omg  = 1d0    !Damping for ZEBRA smoothing

        character(2) :: pc_type  = 'id'   !Type of preconditioner (application dependent)
     .                 ,sm_type  = 'jb'   !Type of smoother ('jg','gs','gm')

        logical :: mg_glrkin     =.false. !Whether to use Galerking MG
     .            ,mg_gm_coarse  =.false. !Whether to use GMRES as the coarse solve
     .            ,sm_line_relax =.false. !Whether to perform line relaxation in MG
     .            ,sm_zebra_relax=.false. !Whether to perform ZEBRA line relaxation in MG
     .            ,mg_debug      =.false. !Whether to go into debugging mode in MG
     .            ,mg_cvrg_test  =.false. !Whether to perform convergence test of smoother
     .            ,pc_asm        =.false. !Whether to do additive Schwartz PC
     .            ,pc_debug      =.false. !PC debugging flag
     .            ,pc_id         =.false. !Whether we are doing identity PC
     .            ,pc_gm_driver  =.false. !Whether GMRES drives PC solves or not
     .            ,pc_as_solver  =.false. !Whether to use the PC as solver (bypassing GMRES)
     .            ,vol_wgt       =.true.  !Whether residuals are volume-weighed or not

      contains

c     TS_solver_consistency_chk
c     ##################################################################
      subroutine TS_solver_consistency_chk

      implicit none

      pc_id = (pc_type == 'id')

cc      if (pc_id) iguess = 0

cc      if (np > 1) mg_gm_coarse = .false.       !Do not perform GMRES coarse solve in parallel

      if (sm_type /= 'jb') sm_omega = 1d0

      if (maxitnwt<0) maxitnwt = max(floor(1.5*log(rtol)/log(tolgm)),10)

      end subroutine TS_solver_consistency_chk

c     readNKS_PCInput
c     #################################################################
      function readNKS_PCInput(file) result(ioerr)

c     -----------------------------------------------------------------
c     Initializes NK & PC parameters. On output,
c       ioerr = 0 if successful
c       ioerr = -1 if problems with NKS
c       ioerr = -2 if problems with PC
c       ioerr = -3 if problems with both
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ioerr
      character(*) :: file

c     Local variables

      integer :: uinput=1234

c     Namelists

      !Newton-Krylov
      namelist /nks/   rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,maxitgm,global,method,damp,dt0,vol_wgt,fpa
     .                   ,iguess,pc_as_solver,nk_def_corr

      !Preconditioner
      namelist /pc/    pc_type,mg_vcyc,sm_iter,pc_tol,pc_iter
     .                   ,pc_asm,sm_type,sm_omega,sm_zebra_relax
     .                   ,pc_debug,mg_ores,mg_oprol

c     Begin program

ccc     Set defaults (set in module now)
cc
cc      !NK parameters
cc      tolgm    = 8d-1          ! Inexact Newton parameter
cc      rtol     = 1d-4          ! Newton relative convergence tolerance
cc      atol     = 0d0           ! Newton absolute convergence tolerance
cc      stol     = 0d0           ! Newton update convergence tolerance
cc      mf_eps   = 1d-6          ! Newtom matrix-free differencing parameter
cc      maxitnwt = 20            ! Maximum number of Newton its. (if zero, maxitnwt
cc                               !        is determined in code)
cc      maxksp   = 15            ! Maximum krylov subspace dimension
cc      maxitgm  = maxksp        ! Maximum GMRES iterations
cc      method   = 1             ! Inexact Newton method:
cc                               !   + 0: constant forcing
cc                               !   + other: adaptive (Eisenstat-Walker)
cc
cc      global   = 0             ! Newton's globalization method:
cc                               !   + 0 --> no globalization
cc                               !   + 1 --> linesearch backtracking 
cc                               !   + 2 --> pseudo-transient
cc
cc      damp     = 1d0           ! Damping parameter in Newton (see nk_mod.f doc))
cc      dt0      = 1d30          ! Initial pseudo-transient time step (" " ")     
cc      iguess   = 0             ! Whether preconditioner is used to give initial
cc                               !   guess to GMRES (when =1) 
cc
cc      pc_as_solver= .false.    ! Whether we use PC as solver, instead of preconditioned GMRES
cc
cc      vol_wgt  = .true.        ! Whether residual is volume weighed or not
cc
cc      fpa      = .false.       ! Whether to use Carlson's Fixed-Point-Accelerator instead of JFNK
cc
cc      !Preconditioner parameters
cc      pc_type  = 'id'          ! Type of preconditioner. Currently:
cc                               !   - 'id': identity (default)
cc                               !   - 's1': SI without flow
cc                               !   - 's2': SI with flow
cc      pc_iter  = 1             ! Number of SI iterations in preconditioner
cc      pc_tol   = 1d-3          ! PC solvers convergence tolerance
cc
cc      pc_asm   = .false.       ! Whether we are doing additive Schwartz PC (in parallel)
cc      pc_debug = .false.       ! PC debugging flag
cc
cc      mg_vcyc  = 1             ! Maximum number of MG V-cycles
cc      mg_ores  = 0             ! Restriction order for MG
cc      mg_oprol = 2             ! Prolongation order for MG
cc
cc      sm_iter  = 5             ! Number of MG smoothing passes
cc      sm_zebra_relax = .false. ! Whether to use ZEBRA relaxation

c     Read initialization parameters

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,nks,iostat=ioerr)
      if (ioerr/=0) then
        ioerr = -1
cc        call pstop('readInput','Problem reading "nks"')
      endif

      read(uinput,pc,iostat=ioerr)
      if (ioerr/=0) then
        ioerr = ioerr-2
cc        call pstop('readInput','Problem reading "pc"')
      endif

      close(unit=uinput)

c     Consistency check

      if (pc_as_solver) then
         iguess = 1
         maxitgm = 0
cc         maxitnwt = 10
         rtol = 1d-1
      endif

      call TS_solver_consistency_chk

ccc     Write initialization parameters to standard output
cc
cc      if (my_rank == 0) then
cc        write (*,*) '----------------- '
cc        write (*,*) 'NKS configuration'
cc        write (*,*) '----------------- '
cc        write (*,nks)
cc
cc        write (*,*) '----------------- '
cc        write (*,*) 'PC  configuration'
cc        write (*,*) '----------------- '
cc        write (*,pc)
cc      endif

c     End program

      end function readNKS_PCInput

      end module TS_solver_setup

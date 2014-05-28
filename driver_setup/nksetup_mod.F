c module nk_setup
c ######################################################################
      module nk_setup

        use grid, ONLY: np,my_rank,pstop

        integer :: maxitnwt
     .            ,maxitgm
     .            ,iguess
     .            ,maxksp
     .            ,method
     .            ,global

        real(8) :: tolgm=0d0
     .            ,rtol=0d0
     .            ,atol=0d0
     .            ,stol=0d0
     .            ,mf_eps=1d-8
     .            ,damp=1d0
     .            ,dt0=1d30

        !FPA configuration
        logical :: fpa=.false.
        integer :: fpa_it=5

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

c     nk_setup_consistency
c     ##################################################################
      subroutine nk_setup_consistency

      implicit none

      pc_id = (pc_type == 'id')

cc      if (pc_id) iguess = 0

      if (np > 1) mg_gm_coarse = .false.       !Do not perform GMRES coarse solve in parallel

      if (sm_type /= 'jb') sm_omega = 1d0

      end subroutine nk_setup_consistency

c     readNKS_PCInput
c     #################################################################
      subroutine readNKS_PCInput(file)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: file

c     Local variables

      integer :: ioerr,uinput=1234

c     Namelists

      !Newton-Krylov
      namelist /nks/   rtol,atol,stol,mf_eps,maxitnwt,tolgm,maxksp
     .                   ,maxitgm,global,method,damp,dt0,vol_wgt,fpa
     .                   ,iguess,pc_as_solver

      !Preconditioner
      namelist /pc/    pc_type,mg_vcyc,sm_iter,pc_tol,pc_iter
     .                   ,pc_asm,sm_type,sm_omega,sm_zebra_relax
     .                   ,pc_debug,mg_ores,mg_oprol

c     Begin program

c     Set defaults

      !NK parameters
      tolgm    = 8d-1          ! Inexact Newton parameter
      rtol     = 1d-4          ! Newton relative convergence tolerance
      atol     = 0d0           ! Newton absolute convergence tolerance
      stol     = 0d0           ! Newton update convergence tolerance
      mf_eps   = 1d-6          ! Newtom matrix-free differencing parameter
      maxitnwt = 20            ! Maximum number of Newton its. (if zero, maxitnwt
                               !        is determined in code)
      maxksp   = 15            ! Maximum krylov subspace dimension
      maxitgm  = maxksp        ! Maximum GMRES iterations
      method   = 1             ! Inexact Newton method:
                               !   + 0: constant forcing
                               !   + other: adaptive (Eisenstat-Walker)

      global   = 0             ! Newton's globalization method:
                               !   + 0 --> no globalization
                               !   + 1 --> linesearch backtracking 
                               !   + 2 --> pseudo-transient

      damp     = 1d0           ! Damping parameter in Newton (see nk_mod.f doc))
      dt0      = 1d30          ! Initial pseudo-transient time step (" " ")
      iguess   = 0             ! Whether preconditioner is used to give initial
                               !   guess to GMRES (when =1) 

      pc_as_solver= .false.    ! Whether we use PC as solver, instead of preconditioned GMRES

      vol_wgt  = .true.        ! Whether residual is volume weighed or not

      fpa      = .false.       ! Whether to use Carlson's Fixed-Point-Accelerator instead of JFNK

      !Preconditioner parameters
      pc_type  = 'id'          ! Type of preconditioner. Currently:
                               !   - 'id': identity (default)
                               !   - 's1': SI without flow
                               !   - 's2': SI with flow
      pc_iter  = 1             ! Number of SI iterations in preconditioner
      pc_tol   = 1d-3          ! PC solvers convergence tolerance

      pc_asm   = .false.       ! Whether we are doing additive Schwartz PC (in parallel)
      pc_debug = .false.       ! PC debugging flag

      mg_vcyc  = 1             ! Maximum number of MG V-cycles
      mg_ores  = 0             ! Restriction order for MG
      mg_oprol = 2             ! Prolongation order for MG

      sm_iter  = 5             ! Number of MG smoothing passes
      sm_zebra_relax = .false. ! Whether to use ZEBRA relaxation

c     Read initialization parameters

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,nks,iostat=ioerr)
      if (ioerr/=0) call pstop('readInput','Problem reading "nks"')

      read(uinput,pc,iostat=ioerr)
      if (ioerr/=0) call pstop('readInput','Problem reading "pc"')

      close(unit=uinput)

c     Consistency check

      if (pc_as_solver) then
         iguess = 1
         maxitgm = 0
cc         maxitnwt = 10
         rtol = 1d-1
      endif

      call nk_setup_consistency

c     Write initialization parameters to standard output

      if (my_rank == 0) then
        write (*,*) '----------------- '
        write (*,*) 'NKS configuration'
        write (*,*) '----------------- '
        write (*,nks)

        write (*,*) '----------------- '
        write (*,*) 'PC  configuration'
        write (*,*) '----------------- '
        write (*,pc)
      endif

c     End program

      end subroutine readNKS_PCInput

      end module nk_setup

c module mlsolverSetup
c######################################################################
      module mlsolverSetup

        use grid_mg

        implicit none

        type :: solver_options
          !Global quantities
          integer    :: iter
          real(8)    :: tol,atol,stol
cc          logical    :: vol_res

          !Stationary iterative methods quantities
          real(8)    :: omega
          real(8)    :: omega10
          real(8)    :: omega01
          integer    :: ncolors
          logical    :: fdiag
          real(8), pointer, dimension(:,:) :: diag

          !Krylov methods quantities
          integer    :: stp_test
          logical    :: sym_test
          integer    :: ngrd_tst
          integer    :: krylov_subspace
          logical    :: singular_matrix

          !MG quantities
          integer      :: vcyc
          integer      :: orderres
          integer      :: orderprol
          integer      :: mg_coarse_solver_depth
          integer      :: mg_coarse_grid_level
cc          integer      :: mg_coarse_grid_res
          integer      :: mg_mu

          logical      :: mg_line_relax
          integer      :: mg_line_nsweep
          integer      :: mg_line_vcyc
          integer      :: mg_line_coarse_solver_depth
          logical      :: mg_line_x
          logical      :: mg_line_y
          logical      :: mg_line_z
          character(2) :: mg_line_solve
          real(8)      :: mg_line_tol
          real(8)      :: mg_line_omega

          logical      :: mg_vertex_relax

          logical      :: mg_zebra_relax
          integer      :: mg_zebra_prefd
          integer      :: mg_zebra_it
          real(8)      :: mg_zebra_omega

          logical      :: mg_galerkin
          logical      :: mg_debug
          logical      :: mg_smooth_only

          logical      :: is_coarse_proc_MG_solve

          type(grid_mg_def),pointer  :: mg_grid_def

          type(mg_ctx),pointer :: mgctx

          !Output quantities
          integer    :: iter_out
          real(8)    :: tol_out
        end type solver_options

        type:: solver_unit
          character*2            :: solver
          type (solver_options)  :: options
        end type solver_unit

        type:: node_type
          type (solver_unit)         :: solver_def
          type (node_type), pointer  :: next_solver
        end type node_type

        type:: queue_type
          type (node_type), pointer:: front
          type (node_type), pointer:: rear
        end type queue_type

        type (queue_type) :: solver_queue

        type (solver_options) :: solverOptions

      contains

c     find_mach_eps
c     ###############################################################
      function find_mach_eps() result(epsmac)

      implicit none

c     ---------------------------------------------------------------
c     Finds machine round-off constant epsmac
c     ---------------------------------------------------------------

c     Call variables

      real(8) :: epsmac

c     Local variables

      real(8) :: mag,mag2

      mag = 1d0
      do
        epsmac = mag
        mag = mag/2
        mag2 = 1d0 + mag
        if (.not.(mag2 > 1d0)) exit
      enddo

      end function find_mach_eps

c     solverOptionsInit
c     ###################################################################
        subroutine solverOptionsInit

c       Initializes solver options

          !Test options
          solverOptions%sym_test = .false.         !Whether to perform symmetry 
                                                   !  test on matvec operator
          solverOptions%ngrd_tst = 1               !Grid level to perform symmetry test

          !Generic options
          solverOptions%iter  = 10                 !Number of iterations
          solverOptions%tol   = 0d-5               !Relative convergence tolerance
          solverOptions%atol  = 0d-5               !Absolute convergence tolerance
          solverOptions%stol  = 0d-5               !Update convergence tolerance

          !MG and smoother options
          solverOptions%vcyc     = 1               !Number of V-cycles (MG)
cc          solverOptions%mg_coarse_grid_res = 2     !Minimum grid level considered (mg_ratio^()) (MG)
          solverOptions%mg_coarse_grid_level = 0   !Grid level considered to be coarsest (MG)
          solverOptions%orderres = 0               !Interpolation order in restriction (MG)
          solverOptions%orderprol= 0               !Interpolation order in prolongation (MG)
          solverOptions%fdiag    = .true.          !Whether to form matrix diagonal
                                                   !  for smoothing
cc          solverOptions%vol_res  = .true.          !Whether residual contains volume information
          nullify(solverOptions%diag)              !Diagonal not provided externally
          solverOptions%omega = 1d0                !Relaxation parameter
          solverOptions%omega10= 0d0               !Weighed Jacobi relaxation parameter
          solverOptions%omega01= 0d0               !Weighed Jacobi relaxation parameter
          solverOptions%ncolors= 1                 !Number of colors in colored GS
          solverOptions%mg_coarse_solver_depth= 0  !Identifies coarse solver for MG 
                                                   !  (if =0, use smoother)
          solverOptions%mg_mu  = 1                 !Identifies MG cycle: V-cycle (mu=1),
                                                   !                     W-cycle (mu=2)
          solverOptions%mg_line_relax=.false.      !Whether to do point-wise 
                                                   !  relaxation (false) or line-wise 
                                                   !  relaxation (true).
          solverOptions%mg_line_nsweep = 1         !Number of line relaxation sweeps
          solverOptions%mg_line_vcyc   = 10        !Number of V-cycles in MG line solver
          solverOptions%mg_line_tol    = 1d-1      !Tolerance for line solves
          solverOptions%mg_line_omega  = 1d0       !Damping parameter for line relax
          solverOptions%mg_line_coarse_solver_depth= 0  !Identifies coarse solver for line solves
                                                        !  (if =0, use smoother)
          solverOptions%mg_line_solve  = "mg"      !Which line solver to use (mg,gm,gs,jb)
          solverOptions%mg_line_x      = .true.    !Whether to do lines in x-direction
          solverOptions%mg_line_y      = .true.    !Whether to do lines in y-direction
          solverOptions%mg_line_z      = .true.    !Whether to do lines in z-direction

          solverOptions%mg_vertex_relax= .false.   !Whether to do vertex-based relax.

          solverOptions%mg_zebra_relax = .false.   !Whether to do zebra relax.
          solverOptions%mg_zebra_prefd = 0         !Whether there is a preferred zebra direction
          solverOptions%mg_zebra_it    = 1         !Number of smoothing passes per zebra line
          solverOptions%mg_zebra_omega = 1d0       !Relaxation constant per zebra line

          solverOptions%mg_galerkin    =.false.    !Whether to do Galerkin coarsening (true)
                                                   !  or rediscretization (false)

          solverOptions%mg_smooth_only =.false.    !Whether to perform only smoothing

          solverOptions%mg_debug       =.false.    !Whether to turn on MG graphic debugging

          solverOptions%mg_grid_def => null()      !Defines default MG grid levels def.

          solverOptions%mgctx => null()            !Nullify MG context

          solverOptions%is_coarse_proc_MG_solve = .false. !Whether we are doing coarse-proc MG solve

          !Krylov methods options
          solverOptions%stp_test = 0               !Stopping criterion (CG, GMRES)
                                                   !If one, use initial residual; 
                                                   !  else, use rhs.
          solverOptions%krylov_subspace = 15       !Krylov subspace dimension (GMRES)

          solverOptions%singular_matrix = .false.  !Matrices are regular by default

          !Output
          solverOptions%iter_out = 0               !Number of iterations (output)
          solverOptions%tol_out  = 0d0             !Convergence achieved (output)

        end subroutine solverOptionsInit

c     xferSolverOptions
c     ###################################################################
        subroutine xferSolverOptions(sopts_bckup,sopts)

        type (solver_options) :: sopts_bckup,sopts

c       Initializes solver options

          !Test options
          sopts_bckup%sym_test = sopts%sym_test
                                                   
          sopts_bckup%ngrd_tst = sopts%ngrd_tst

          !Generic options
          sopts_bckup%iter  = sopts%iter
          sopts_bckup%tol   = sopts%tol 
          sopts_bckup%atol  = sopts%atol
          sopts_bckup%stol  = sopts%stol

          !MG and smoother options
          sopts_bckup%vcyc                 = sopts%vcyc                
cc          sopts_bckup%mg_coarse_grid_res   = sopts%mg_coarse_grid_res  
          sopts_bckup%mg_coarse_grid_level = sopts%mg_coarse_grid_level
          sopts_bckup%orderres             = sopts%orderres            
          sopts_bckup%orderprol            = sopts%orderprol           
          sopts_bckup%fdiag                = sopts%fdiag               
                                                   
cc          sopts_bckup%vol_res  = .true.        
          sopts_bckup%diag => sopts%diag
          sopts_bckup%omega   = sopts%omega  
          sopts_bckup%omega10 = sopts%omega10
          sopts_bckup%omega01 = sopts%omega01
          sopts_bckup%ncolors = sopts%ncolors
          sopts_bckup%mg_coarse_solver_depth
     .        = sopts%mg_coarse_solver_depth
                                                   
          sopts_bckup%mg_mu  = sopts%mg_mu
                                                   
          sopts_bckup%mg_line_relax = sopts%mg_line_relax 
                                                   
          sopts_bckup%mg_line_nsweep = sopts%mg_line_nsweep
          sopts_bckup%mg_line_vcyc   = sopts%mg_line_vcyc  
          sopts_bckup%mg_line_tol    = sopts%mg_line_tol   
          sopts_bckup%mg_line_omega  = sopts%mg_line_omega 
          sopts_bckup%mg_line_coarse_solver_depth
     .        = sopts%mg_line_coarse_solver_depth
                                                        
          sopts_bckup%mg_line_solve  = sopts%mg_line_solve
          sopts_bckup%mg_line_x      = sopts%mg_line_x    
          sopts_bckup%mg_line_y      = sopts%mg_line_y    
          sopts_bckup%mg_line_z      = sopts%mg_line_z    

          sopts_bckup%mg_vertex_relax= sopts%mg_vertex_relax

          sopts_bckup%mg_zebra_relax = sopts%mg_zebra_relax 
          sopts_bckup%mg_zebra_prefd = sopts%mg_zebra_prefd 
          sopts_bckup%mg_zebra_it    = sopts%mg_zebra_it    
          sopts_bckup%mg_zebra_omega = sopts%mg_zebra_omega 

          sopts_bckup%mg_galerkin    = sopts%mg_galerkin
          sopts_bckup%mg_smooth_only = sopts%mg_smooth_only
                                                   
          sopts_bckup%mg_debug       = sopts%mg_debug 

          sopts_bckup%mg_grid_def    => sopts%mg_grid_def

          sopts_bckup%mgctx          => sopts%mgctx   

          sopts_bckup%is_coarse_proc_MG_solve
     .                               = sopts%is_coarse_proc_MG_solve

          !Krylov methods options
          sopts_bckup%stp_test = sopts%stp_test    
                                                   
          sopts_bckup%krylov_subspace = sopts%krylov_subspace

          sopts_bckup%singular_matrix = sopts%singular_matrix  

          !Output
          sopts_bckup%iter_out = sopts%iter_out            
          sopts_bckup%tol_out  = sopts%tol_out            

        end subroutine xferSolverOptions

c     solverInit
c     ###################################################################
        subroutine solverInit

c       Initializes solver

          call init_queue (solver_queue)

        end subroutine solverInit

c     solverKill
c     ###################################################################
        subroutine solverKill

c       Deallocates solver queue

          implicit none
    
          type (node_type), pointer :: temp

        !Begin

          do while (.not.is_empty(solver_queue))
            temp => solver_queue%front%next_solver

            deallocate (solver_queue%front)

            solver_queue%front => temp
          enddo

        end subroutine solverKill
     
c     assembleSolverHierarchy
c     ###################################################################
        subroutine assembleSolverHierarchy (solver)

c       Assembles solver queue

          character*2 :: solver

          type (solver_unit) :: solver_def

        !Begin

          solver_def%solver = solver
          solver_def%options = solverOptions

          call put_in_rear(solver_queue,solver_def)

        end subroutine assembleSolverHierarchy

c     readSolverHierarchy
c     ###################################################################
        subroutine readSolverHierarchy (solver_def,depth)

c       Reads TOP solver definition from solver hierarchy

          type (solver_unit) :: solver_def
          integer            :: depth

        !Begin

          call read_node (solver_queue,solver_def,depth)

        end subroutine readSolverHierarchy

c     writeSolverHierarchy
c     ###################################################################
        subroutine writeSolverHierarchy (solver_def,depth)

c       Modifies solver definition at depth 'depth' in solver hierarchy

          type (solver_unit) :: solver_def
          integer :: depth

        !Begin

          call write_node (solver_queue,solver_def,depth)

        end subroutine writeSolverHierarchy

c     getSolverOptions
c     ###################################################################
        subroutine getSolverOptions(depth)

c       Gets solver options at depth 'depth'

          integer :: depth
          type (solver_unit) :: solver_def

        !Begin

          call readSolverHierarchy (solver_def,depth)
          solverOptions = solver_def%options

        end subroutine getSolverOptions

c     printSolverHierarchy
c     ###################################################################
cc        subroutine printSolverHierarchy

c       Reads TOP solver definition from solver hierarchy

cc          type (solver_unit) :: solver_def
cc          integer               :: depth,idepth

        !Begin

cccc          idepth = 0

cc          depth = count_elements(solver_queue)

cc          do idepth = 1,depth

cccc            if (is_empty(solver_queue)) then
cccc              write (*,*) 'No more solver definitions present'
cccc              write (*,*) 'Aborting...'
cccc              stop
cccc              exit
cccc            endif

cccc            idepth = idepth + 1

cccc            call take_from_front (solver_queue,solver_def)

cc            call read_node (solver_queue,solver_def,idepth)

cc            write (*,10) idepth
cc            write (*,*) solver_def
cc            write (*,*)
cc          enddo

cc 10       format ('Solver level:',i3,/,'Options:')

cc        end subroutine printSolverHierarchy

c     init_queue
c     ###################################################################
        subroutine init_queue(q)

c       Initializes pointer queue

          type (queue_type), intent (out) :: q

        !Begin

          nullify (q%front,q%rear)

        end subroutine init_queue

c     put_in_rear
c     ###################################################################
        subroutine put_in_rear(q,buffer)

c       Appends a node to a queue

          type (queue_type), intent(in out):: q
          type (solver_unit), intent(in)   :: buffer

        !Begin

          if (.not.associated(q%front)) then

            allocate (q%front)
            q%rear => q%front

          else

            allocate(q%rear%next_solver)
            q%rear => q%rear%next_solver

          endif

          q%rear%solver_def = buffer

          nullify(q%rear%next_solver)

        end subroutine put_in_rear

c     count_elements
c     ###################################################################
        function count_elements(q) result (count)

c       Counts number of elements in a queue

          integer    :: count
          type (queue_type), intent (in) :: q
          type (node_type), pointer      :: node_ptr

        !Begin

          count = 0

          node_ptr => q%front

        !Traverse the list

          do while (associated(node_ptr))
            count = count + 1
            node_ptr => node_ptr%next_solver
          enddo

        end function count_elements

c     countSolverLevels
c     ###################################################################
        subroutine countSolverLevels(count)

c       Counts number of elements in a queue

          integer    :: count

        !Begin

          count = count_elements(solver_queue)

        end subroutine countSolverLevels

c     take_from_front
c     ###################################################################
        subroutine take_from_front (q, buffer)

c       Deletes front node from queue

          type (queue_type), intent (in out) ::q
          type (solver_unit), intent (out)   ::buffer

          type (node_type), pointer :: temp

        !Begin

          buffer = q%front%solver_def

          temp => q%front%next_solver

          deallocate (q%front)

          q%front => temp

        end subroutine take_from_front

c     read_node
c     ###################################################################
        subroutine read_node (q, buffer,depth)

c       Deletes front node from queue

          integer    :: depth
          type (queue_type), intent (in out) ::q
          type (solver_unit), intent (out)   ::buffer

          integer    :: count
          type (node_type), pointer :: node_ptr

        !Begin

          count = 0

          node_ptr => q%front

        !Traverse the list

          do while (associated(node_ptr))
            count = count + 1
            if (count == depth) exit
            node_ptr => node_ptr%next_solver
          enddo

        !Read solver definition

          if (count == depth) then
            buffer = node_ptr%solver_def
          else
            write (*,*) 'Solver depth =',depth,' unreachable'
            write (*,*) 'Aborting...'
            stop
          endif

        end subroutine read_node

c     write_node
c     ###################################################################
        subroutine write_node (q,buffer,depth)

c       Deletes front node from queue

          integer    :: depth
          type (queue_type), intent (in out) ::q
          type (solver_unit), intent (in)   ::buffer

          integer    :: count
          type (node_type), pointer :: node_ptr

        !Begin

          count = 0

          node_ptr => q%front

        !Traverse the list

          do while (associated(node_ptr))
            count = count + 1
            if (count == depth) exit
            node_ptr => node_ptr%next_solver
          enddo

        !Write solver definition

          if (count == depth) then
cc            node_ptr%solver_def = buffer
            node_ptr%solver_def%solver = buffer%solver
            call xferSolverOptions(node_ptr%solver_def%options
     .                            ,buffer%options)
          else
            write (*,*) 'Solver depth =',depth,' unreachable'
            write (*,*) 'Aborting...'
            stop
          endif

        end subroutine write_node

c     is_empty
c     ###################################################################
        function is_empty (q) result (empty)

c       Determines if queue is empty

          logical empty

          type (queue_type), intent (in) :: q

          empty = .not.associated(q%front)

        end function is_empty

      end module mlsolverSetup

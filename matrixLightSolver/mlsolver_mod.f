c module mlsolverSetup
c######################################################################
      module mlsolverSetup

        use grid

        implicit none

        type:: solver_options
          integer         :: iter
          integer         :: iter_out
          double precision:: tol
          double precision:: tol_out
          double precision:: omega
          double precision:: omega10
          double precision:: omega01
          integer         :: ncolors
          integer         :: vcyc
          integer         :: igridmin
          integer         :: stp_test
          logical         :: sym_test
          logical         :: fdiag
          logical         :: vol_res
          logical         :: mg_line_relax
          logical         :: vertex_based_relax
          integer         :: krylov_subspace
          integer         :: orderres
          integer         :: orderprol
          integer         :: mg_coarse_solver_depth
          integer         :: mg_mu
          type(grid_def)  :: mg_grid_def
          double precision, pointer, dimension(:,:) :: diag
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

c     solverOptionsInit
c     ###################################################################
        subroutine solverOptionsInit

c       Initializes solver options

          !Test options
          solverOptions%sym_test = .false.         !Whether to perform symmetry 
                                                   !  test on matvec operator

          !Generic options
          solverOptions%iter  = 10                 !Number of iterations
          solverOptions%tol   = 1d-5               !Convergence tolerance

          !MG and smoother options
          solverOptions%vcyc     = 1               !Number of V-cycles (MG)
          solverOptions%igridmin = 2               !Minimum grid level considered (MG)
          solverOptions%orderres = 0               !Interpolation order in restriction (MG)
          solverOptions%orderprol= 0               !Interpolation order in prolongation (MG)
          solverOptions%fdiag    = .true.          !Whether to form matrix diagonal
                                                   !  for smoothing
          solverOptions%vol_res  = .false.         !Whether residual contains volume information
          nullify(solverOptions%diag)              !Diagonal not provided externally
          solverOptions%omega = 1d0                !Relaxation parameter
          solverOptions%omega10= 0d0               !Weighed Jacobi relaxation parameter
          solverOptions%omega01= 0d0               !Weighed Jacobi relaxation parameter
          solverOptions%ncolors= 1                 !Number of colors in colored GS
          solverOptions%mg_coarse_solver_depth= 0  !Identifies coarse solver for MG 
                                                   !  (if =0, use smoother)
          solverOptions%mg_mu  = 1                 !Identifies MG cycle: V-cycle (mu=1),
                                                   !                     W-cycle (mu=2)
          solverOptions%mg_line_relax=.false.      !Specifies whether to do point-wise 
                                                   !  relaxation (false) or line-wise 
                                                   !  relaxation (true).
          solverOptions%vertex_based_relax=.false. !Specifies whether to do vertex-based relax.
          solverOptions%mg_grid_def = grid_params  !Defines default MG grid levels def.

          !Krylov methods options
          solverOptions%stp_test = 0               !Stopping criterion (CG, GMRES)
                                                   !If one, use initial residual; 
                                                   !  else, use rhs.
          solverOptions%krylov_subspace = 15       !Krylov subspace dimension (GMRES)

          !Output
          solverOptions%iter_out = 0               !Number of iterations (output)
          solverOptions%tol_out  = 0d0             !Convergence achieved (output)

        end subroutine solverOptionsInit

c     solverInit
c     ###################################################################
        subroutine solverInit

c       Initializes solver

          call init_queue (solver_queue)

        end subroutine solverInit

c     solverKill
c     ###################################################################
        subroutine solverKill

c       Kills solver

          if (associated(solver_queue%front) )
     .             deallocate(solver_queue%front)
cc          if (associated(solver_queue%rear ) )
cc     .             deallocate(solver_queue%rear)

cc          deallocate (solver_queue%front,solver_queue%rear)

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
          integer            :: depth

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
cc          integer            :: depth,idepth

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

          integer :: count
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

          integer :: count

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

          integer :: depth
          type (queue_type), intent (in out) ::q
          type (solver_unit), intent (out)   ::buffer

          integer :: count
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

          integer :: depth
          type (queue_type), intent (in out) ::q
          type (solver_unit), intent (in)   ::buffer

          integer :: count
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
            node_ptr%solver_def = buffer
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

c module variables_bc
c #########################################################################
      module variables_bc

        use BCS_variables
        use variable_setup

      contains

#if !defined(samrai) && !defined(OLD_BC)

c     applyBC
c     ###############################################################
      subroutine applyBC(igr,varray,it,bc_order)
c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igr
      integer,optional :: it,bc_order

      type(var_array),pointer :: varray

c     Local variables

      integer :: i,order,igrp,itm

c     Interfaces

      INTERFACE
        subroutine applyDepVarBC(igr,bc_inf,varray,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          integer,optional :: order
        end subroutine applyDepVarBC
      END INTERFACE

      INTERFACE
        subroutine applyAuxVarBC(igr,bc_inf,varray,vaux,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyAuxVarBC
      END INTERFACE

      INTERFACE
        subroutine preproc_AuxVar(igr,bc_inf,varray,vaux,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine preproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine postproc_AuxVar(igr,bc_inf,varray,vaux)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine postproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use variable_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

      if (PRESENT(it)) then
        itm = it
      else
        itm = 0
      endif

      if (PRESENT(bc_order))  then
        order = bc_order
      else
        order = 1
      endif

c     Set application BCs (external)

      call setAppBCs(varray,gv%aux)

c     Set BC schedule

      call setBCsched(varray,gv%aux,itm)

c     Impose BCs

      do igrp = 1,gv%nbc_grp
        !Fill auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call preproc_AuxVar(igr,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                         ,gv%aux,order=order)
          endif
        enddo

        !Impose BCs on group
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX < 0) then     !Process dep. variables
            call applyDepVarBC(igr,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                        ,order=order)
          elseif (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then !Process aux. variables
            call applyAuxVarBC(igr,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                        ,gv%aux,order=order)
          endif
        enddo

        !Postprocess auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call postproc_AuxVar(igr,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                          ,gv%aux)
          endif
        enddo

      enddo

c     End program

      end subroutine applyBC

c     applyBCplot
c     ##############################################################
      subroutine applyBCplot(igr,varray,it,bc_order)
c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igr
      type(var_array),pointer :: varray
      integer, optional :: it,bc_order

c     Local variables

      integer :: i,order,igrp,itm

c     Interfaces

      INTERFACE
        subroutine preproc_AuxVar(igr,bc_inf,varray,vaux,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine preproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine postproc_AuxVar(igr,bc_inf,varray,vaux)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine postproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use variable_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

      if (PRESENT(it)) then
        itm = it
      else
        itm = 0
      endif

      if (PRESENT(bc_order))  then
        order = bc_order
      else
        order = 1
      endif

c     Set application BCs (external)

      call setAppBCs(varray,gv%aux)

c     Set BC schedule

      call setBCsched(varray,gv%aux,itm)

c     Impose BCs

      do igrp = 1,gv%nbc_grp
        !Fill auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call preproc_AuxVar(igr,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                         ,gv%aux,order=order)
          endif
        enddo

        !Postprocess auxiliary quantities
        do i = 1,gv%bc_grp(igrp)%nbc_seq
          if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
            call postproc_AuxVar(igr,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                          ,gv%aux)
          endif
        enddo

      enddo

c     End program

      end subroutine applyBCplot

#else

#if defined(OLD_BC)

c     applyBC
c     ###############################################################
      subroutine applyBC(igr,varray,it)
c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igr
      integer,optional :: it

      type(var_array),pointer :: varray

c     Local variables

c     Interfaces

      INTERFACE
        subroutine apply_App_BC(igr,varray)
          use variable_setup
          integer    :: igr
          type(var_array),pointer :: varray
        end subroutine apply_App_BC
      END INTERFACE

c     Begin program

      call apply_App_BC(igr,varray)

c     End program

      end subroutine applyBC
#endif

#endif

c     setBCsched
c     ###################################################################
      subroutine setBCsched(varray,vaux,it)

c     -------------------------------------------------------------------
c     Sets BC schedule for dependent and auxiliary variables, which
c     are identified by a 2-dim integer array, bc_seq. The first
c     dimension is of size nbc_seq, and indicates order of treatment
c     lexicographically (i.e., first component of first dimension is
c     to be treated first, second component of first dimension is
c     to be treated second, and so on).
c
c     The second dimension is of size two, and has the following convention:
c       * First component of the second dimension identifies variable,
c         and its sign indicates whether the variable is dependent (>0) or
c         auxiliary (<0).
c       * Second component of the second dimension identifies whether
c         variable is a scalar (=0) or a vector (=1).
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: it

      type(var_array),pointer :: varray

      type(aux_array),pointer :: vaux

c     Local variables

      logical :: dump=.false.

c     Begin program

#if defined(samrai)
      call create_bc_schedule(varray,vaux,gv%nbc_grp,gv%bc_grp
     .                       ,.true.,.false.,dump,.false.)
#else
cc      dump = (it==0).and.(my_rank==0)
      call create_bc_schedule(varray,vaux,gv%nbc_grp,gv%bc_grp
     .                       ,.false.,.false.,dump,.false.)
#endif

c     End program

      end subroutine setBCsched

c     create_bc_schedule
c     ######################################################################
      subroutine create_bc_schedule(u,aux,nbc_grp,bc_grp,inherit
     .                             ,interior_only,print_used,reverse)

c     ----------------------------------------------------------------------
c     Creates BC schedule based on dependencies within domain and at
c     boundaries. Input/output explained below.  Contributed by Mark
c     Berrill, 2011
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      !Input variables
      type(var_array) :: u    ! We need nvar, bc_dep_list, dom_dep_list, descr for each variable
      type(aux_array) :: aux  ! We need nvar, nvec, bc_dep_list, dom_dep_list, descr for each variable

      logical :: inherit      ! Do we want to inherit the dependencies for the groups (set to false)
      logical :: interior_only! Do we want to allow some variables to fill interiors only
                              !     (potential performance increase if true)
      logical :: print_used   ! Do we want to print how the variables are used
      logical :: reverse      ! Do we want to reverse the default order within a group
                              !     (useful for debugging, dependencies will always be preserved)

      !Output variables
      integer :: nbc_grp                    ! Number of boundary condition groups
      type(bcgrp), pointer :: bc_grp(:)     ! Number of boundary condition groups

c     Local variables

      integer :: N, i, j, k, m, nn, N_groups
      integer :: error, n_var(4)
      logical :: test
      INTEGER, POINTER :: dep_list(:,:)
      integer,pointer :: bc_dep(:,:) ,dom_dep(:,:), tot_dep(:,:)    ! Pointers used to store the dependency tables
      integer,pointer :: n_seq(:), bc_groups(:,:), used(:)          ! Pointers used to store the internal group information
      integer,pointer :: tmp1(:), tmp2(:,:)                         ! Pointers used for temporary memory

c     Begin program

      n_var(1) = u%nvar     ! Number of scalar dependent variables
      n_var(2) = 0          ! Number of vector dependent variables
      n_var(3) = aux%nvar   ! Number of scalar auxiliary variables
      n_var(4) = aux%nvec   ! Number of vector auxiliary variables
      N = n_var(1) + n_var(2) + n_var(3) + n_var(4)

c     Search for invalid or duplicate entries (this is bad practice and may indicate an error)
      error = 0;
      do i = 1,n_var(1)
         if ( size(u%array_var(i)%bc_dep_list,1)/=N .OR. 
     &        size(u%array_var(i)%dom_dep_list,1)/=N ) then
            call pstop('create_bc_schedule',
     &         'Unexpected size for dependency list')
         endif
         error = check_dep_list(N,u%array_var(i)%bc_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in dependent variable (bc_dep_list)')
         end if
         error = check_dep_list(N,u%array_var(i)%dom_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in dependent variable (dom_dep_list)')
         end if
      end do
      do i = 1,n_var(2)
         !!!! Place holder for dependent vectors
      end do
      do i = 1,n_var(3)
         if ( size(aux%var_list(i)%bc_dep_list,1)/=N .OR. 
     &        size(aux%var_list(i)%dom_dep_list,1)/=N ) then
            call pstop('create_bc_schedule',
     &         'Unexpected size for dependency list')
         endif
         error = check_dep_list(N,aux%var_list(i)%bc_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in scalar auxiliary '//
     &      'variable (bc_dep_list)')
         endif
         error = check_dep_list(N,aux%var_list(i)%dom_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in scalar auxiliary '//
     &      'variable (dom_dep_list)')
         endif
      enddo
      do i = 1,n_var(4)
         if ( size(aux%vec_list(i)%bc_dep_list,1)/=N .OR. 
     &        size(aux%vec_list(i)%dom_dep_list,1)/=N ) then
            call pstop('create_bc_schedule',
     &         'Unexpected size for dependency list')
         endif
         error = check_dep_list(N,aux%vec_list(i)%bc_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in vector auxiliary '//
     &      'variable (bc_dep_list)')
         endif
         error = check_dep_list(N,aux%vec_list(i)%dom_dep_list)
         if ( error .ne. 0 ) then
            call pstop('create_bc_schedule',
     &      'Invalid entry found in vector auxiliary '//
     &      'variable (dom_dep_list)')
         endif
      enddo

c     Create complete dependency tables
      allocate ( bc_dep(N,N) )
      allocate ( dom_dep(N,N) )
      allocate ( tot_dep(N,N) )
      bc_dep = 0
      dom_dep = 0
      tot_dep = 0
c     Add the dependencies for the scalar dependent variables
      do i = 1,n_var(1)
         j = 1
         dep_list => u%array_var(i)%bc_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            bc_dep(i,k) = dep_list(j,3)+1
            j = j+1
            if (j > N) exit
         enddo
         j = 1
         dep_list => u%array_var(i)%dom_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            dom_dep(i,k) = dep_list(j,3)+1
            j = j+1
            if (j > N) exit
         enddo
      enddo
c     Add the dependencies for the vector dependent variables
      do i = 1,n_var(2)
         !!!! Place holder for dependent vectors
      enddo
c     Add the dependencies for the scalar auxiliary variables
      do i = 1,n_var(3)
         j = 1
         dep_list => aux%var_list(i)%bc_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            bc_dep(i+n_var(1)+n_var(2),k) = dep_list(j,3)+1
            j = j+1
         enddo
         j = 1
         dep_list => aux%var_list(i)%dom_dep_list
        do  while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            dom_dep(i+n_var(1)+n_var(2),k) = dep_list(j,3)+1
            j = j+1
         enddo
      enddo
c     Add the dependencies for the vector auxiliary variables
      do i = 1,n_var(4)
         j = 1
         dep_list => aux%vec_list(i)%bc_dep_list
         do while ( dep_list(j,1)/=0 ) 
            k = get_dep_index( dep_list(j,:), n_var )
            bc_dep(i+n_var(1)+n_var(2)+n_var(3),k) = dep_list(j,3)+1
            j = j+1
         enddo
         j = 1
         dep_list => aux%vec_list(i)%dom_dep_list
         do while ( dep_list(j,1)/=0 )
            k = get_dep_index( dep_list(j,:), n_var )
            dom_dep(i+n_var(1)+n_var(2)+n_var(3),k) = dep_list(j,3)+1
            j = j+1
         enddo
      enddo

c     Check which variables are used from the dependency list and fill the total dependency table
c     (0: Variable is not needed, 1: Only interior needs to be filled, 2: Interior and BC are needed)
      allocate ( used(N) )
      used = 0
      tot_dep = 0
c     All dependent variables need to have BC applied 
c     Note: we do not need to respect the interior dependencies for determining the groups and group order)
      nn = n_var(1)+n_var(2)
      used(1:nn) = 2
      tot_dep(1:nn,:) = bc_dep(1:nn,:)
c     Set the dependencies needed for the dependent variables (note: we need to respect the domain and bc dependencies)
      do i = nn+1,N
         used(i) = MAX(MAXVAL(dom_dep(1:nn,i)),MAXVAL(bc_dep(1:nn,i)))
         if ( (.NOT.interior_only) .AND. used(i)==1 ) then
c           We are requiring that we will fill bc on all used variables
            used(i) = 2
         endif
         if ( used(i)==1 ) then
c           We only need to fill the interior domain
            tot_dep(i,:) = dom_dep(i,:)
         else if ( used(i)==2 ) then
c           We only need to fill the interior domain
            do j = 1,N
               tot_dep(i,j) = MAX(dom_dep(i,j),bc_dep(i,j))
            enddo
         endif
      enddo
c     Set the dependencies needed for all the used variables
      test = .TRUE.
      do while ( test )
         test = .FALSE.
         do i = 1,N
            if ( used(i)>0 ) then
               do j = 1,N
                  if ( tot_dep(i,j) > used(j) ) then
                     used(j) = tot_dep(i,j)
                     if ( (.NOT.interior_only) .AND. used(j)==1 ) then
c                       We are requiring that we will fill bc on all used variables
                        used(j) = 2
                     endif
                     if ( used(j)==1 ) then
c                       We only need to fill the interior domain
                        tot_dep(j,:) = dom_dep(j,:)
                     else if ( used(j)==2 ) then
c                       We only need to fill the interior domain
                        do k = 1,N
                           tot_dep(j,k) = MAX(dom_dep(j,k),bc_dep(j,k))
                        enddo
                     endif
                     test = .TRUE.
                  endif
               enddo
            endif
         enddo
      enddo
c     Remove self-dependencies (this can create issues when we do the sorting)
      do i = 1,N
         tot_dep(i,i) = 0
      enddo

c     List the variables that are used
      if ( print_used ) then
         allocate ( tmp1(4) )
c        Determine the number of variables with bc, interior, unused, and dep
         tmp1 = 0
         tmp1(1) = n_var(1)+n_var(2)
         do i = n_var(1)+n_var(2)+1,N
            if ( used(i)==2 ) then
               tmp1(2) = tmp1(2)+1
            else if ( used(i)==1 ) then
               tmp1(3) = tmp1(3)+1
            else if ( used(i)==0 ) then
               tmp1(4) = tmp1(4)+1
            endif
         enddo
         write (*,100) N, tmp1(1), tmp1(2), tmp1(3), tmp1(4)
100      format ( 1x, 1I2, ' total variables: ',
     &                1I2, ' dependent variables, ',
     &                1I2, ' bc applied, ',
     &                1I2, ' interiors only, ',
     &                1I2, ' unused' )
c        Write the auxiliary variables and how they are used
         write (*,*) 'Dependent variables: '
         do i = 1,SUM(n_var(1:2))
             call print_var_name(u,aux,i)
         enddo
         write (*,*) 'Variables that need bc applied: '
         do i = SUM(n_var(1:2)),N
            if ( used(i)==2 ) then
               call print_var_name(u,aux,i)
            endif
         enddo
         write (*,*) 'Variables that only need interiors filled: '
         do i = SUM(n_var(1:2)),N
c LC 1/25/11, small bug:            if ( used(i)==2 ) then
            if ( used(i)==1 ) then
               call print_var_name(u,aux,i)
            endif
         enddo
         write (*,*) 'Variables that are unused: '
         do i = SUM(n_var(1:2)),N
            if ( used(i)==0 ) then
               call print_var_name(u,aux,i)
            endif
         enddo
         deallocate ( tmp1 )
      endif

c     Check if any variables have no dependencies and are used (this is an error)
c     Note that if a variable only depends on its self, it is OK (it might depend on the equilibrium data)
      do i = n_var(1)+n_var(2)+1,N      ! Only auxiliary variables need dependencies
         if ( MAXVAL(tot_dep(i,:))==0 .AND. used(i)>0 .AND. 
     &        bc_dep(i,i)==0 .AND. tot_dep(i,i)==0 ) then
            call pstop('create_bc_schedule',
     &        'Variables are used but do not contain any dependencies')
         endif
      enddo

c     Create a first group containing all used variables
      allocate ( n_seq(N+1) )
      allocate ( bc_groups(N+1,N) )
      n_seq = 0
      bc_groups = 0
      N_groups = 1
      do i = 1,N
         if ( used(i)>0 ) then
            bc_groups(1,n_seq(1)+1) = i
            n_seq(1) = n_seq(1)+1
         endif
      enddo

c     While there are any dependencies, create new groups consisting of those dependencies
      i = 1
      do
c        Identify variables that can be in the same group and variables that must be in an earlier group
         do j = n_seq(i),1,-1
            test = .FALSE.
            do k = 1,n_seq(i)
               if ( tot_dep(bc_groups(i,k),bc_groups(i,j))==2 ) then
                  test = .TRUE.
               endif
            end do
            if ( test ) then
               n_seq(i+1) = n_seq(i+1)+1
               bc_groups(i+1,n_seq(i+1)) = bc_groups(i,j)
               N_groups = i+1
               do k = j+1,n_seq(i)
                  bc_groups(i,k-1) = bc_groups(i,k)
               enddo
               bc_groups(i,n_seq(i)) = 0
               n_seq(i) = n_seq(i)-1
            endif
         enddo
c        Move any variables that are needed by the variables in the earlier group to the earlier group
         test = .TRUE.
         do while ( test ) 
            test = .FALSE.
            j = 1
            do while (j <= n_seq(i+1))
               k = 1
               do while (k <= n_seq(i))
                  if ( tot_dep(bc_groups(i+1,j),bc_groups(i,k))>0 ) then
                     n_seq(i+1) = n_seq(i+1)+1
                     bc_groups(i+1,n_seq(i+1)) = bc_groups(i,k)                     
                     do m = k+1,n_seq(i)
                        bc_groups(i,m-1) = bc_groups(i,m)
                     enddo
                     bc_groups(i,n_seq(i)) = 0
                     n_seq(i) = n_seq(i)-1
                     test = .TRUE.
                  endif
                  k = k + 1
               enddo
               j = j + 1
            enddo 
         enddo
         if ( n_seq(i)==0 ) then
c           All variables are in the earlier group, there must be a circular dependence
            write (*,*) 'Problem separating groups, ', 
     &         'circular dependence detected between: '
            do j = 1,n_seq(i+1)
               call print_var_name(u,aux,bc_groups(i+1,j))
            enddo
            call pstop('create_bc_schedule',
     &        'error creating schedule, circular dependence detected')
         endif 
         if ( n_seq(i+1)==0 ) then
c           We are done creating new groups
            exit
         endif
         if ( i==N ) then
c           We have had too many iterations (this should never occur)
            call pstop('create_bc_schedule','Infinite loop')
         endif
         i = i+1
      enddo

c     Change the group order so it reflects the order in which the variables need to be processed
      allocate ( tmp1(N+1) )
      allocate ( tmp2(N+1,N) )
      tmp1 = n_seq
      tmp2 = bc_groups
      do i = 1,N_groups
         bc_groups(i,:) = tmp2(N_groups-i+1,:)
         n_seq(i) = tmp1(N_groups-i+1)
      end do
      deallocate ( tmp1 )
      deallocate ( tmp2 )

c     Sort the variables in each group in ascending order (to create a consistent default order)
      allocate ( tmp1(N) )
      do i = 1,N_groups
         tmp1 = 0
         do j = 1,n_seq(i)
            tmp1(bc_groups(i,j)) = 1
         enddo
         j = 1
         do k = 1,N
            if ( tmp1(k)==1 ) then
               bc_groups(i,j) = k
               j = j+1
            endif
         enddo
      end do
      deallocate ( tmp1 )

c     Reverse the default order within a group (useful for debugging, dependencies will always be preserved)
      if ( reverse ) then
         do i = 1,N_groups
            allocate ( tmp1(n_seq(i)) )
            tmp1 = bc_groups(i,1:n_seq(i))
            do j = 1,n_seq(i)
               bc_groups(i,j) = tmp1(n_seq(i)-j+1)
            enddo
            deallocate ( tmp1 )
         end do
      endif

c     Sort the schedules within groups to reflect the internal dependencies
      do i = 1,N_groups
         do j = 1,n_seq(i)
c           Get the dependencies between the remaining members of the group
            nn = n_seq(i)-j+1
            allocate ( tmp1(nn) )
            allocate ( tmp2(nn,nn) )
            tmp1 = bc_groups(i,j:n_seq(i))
            do k = 1,nn
               do m = 1,nn
                  tmp2(k,m) = tot_dep(tmp1(k),tmp1(m))
               end do
            end do
c           Check that no members of the group depend critically on any other members
            if ( MAXVAL(tmp2)==2 ) then
               call pstop('create_bc_schedule',
     &            'Error in dependency detected')
            end if
c           Find the first variable which does not depend on the others
            k = -1
            do m = nn,1,-1
               if ( MAXVAL(tmp2(m,:))==0 ) then
                  k = m
               endif
            enddo
            if ( k==-1 ) then 
c              The variables all depend weakly on each other, we only need to
c              obey the bc dependencies within a group and re-search
               do k = 1,nn
                  do m = 1,nn
                     tmp2(k,m) = bc_dep(tmp1(k),tmp1(m))
                  end do
                  tmp2(k,k) = 0
               end do
               do m = nn,1,-1
                  if ( MAXVAL(tmp2(m,:))==0 ) then
                     k = m
                  end if
               end do
            end if
            if ( k==-1) then
c              Unable to find a variable that does not depend on the others, a circular dependence is present
               write (*,*) 'Problem ordering groups, ', 
     &            'circular dependence detected between: '
               do m = 1,nn
                  call print_var_name(u,aux,tmp1(m))
               enddo
               call pstop('create_bc_schedule',
     &            'Error sorting group, circular dependence detected')
            end if
            bc_groups(i,j) = tmp1(k)
            bc_groups(i,j+1:j+k-1) = tmp1(1:k-1)
            bc_groups(i,j+k:n_seq(i)) = tmp1(k+1:nn)            
            deallocate ( tmp1 )
            deallocate ( tmp2 )
         end do
      end do
      
c     If we need inherite the dependencies from previous groups, add them
c     This is needed for temporary patches when the src and dst ids are not the same.  
c     In this case we assume all of the interior data is compied sucessfully, 
c     and we only need to worry about bc values on the temporary patch.
      if ( inherit ) then
         allocate ( tmp1(N) )
         do i = 2,N_groups
c           Get an array that indicates if each variable is needed by any member of the current group
            tmp1(:) = 0
            do j = 1,n_seq(i)
               k = bc_groups(i,j)
               if ( used(k) == 2 ) then
                  do m = 1,N
                     if ( bc_dep(k,m)>0 ) then
                        tmp1(m) = 1
                     endif
                  end do
               endif
            enddo
c           We need any dependencies that are necessary to fill the bc
c           for each dependency (interior dependencies should not be needed)
            test = .TRUE.
            do while ( test ) 
               test = .FALSE.
               do j = 1,N
                  if ( tmp1(j)>0 ) then
                     do k = 1,N
                        if ( bc_dep(j,k)>0 .AND. tmp1(k)==0 
     &                        .AND. used(k)==2 ) then
                           tmp1(k) = 1
                           test = .TRUE.
                        endif
                     enddo
                  endif
               enddo
            enddo
c           Remove any dependency variables from the current group (these are met by being members of the group)
            do k = 1,n_seq(i)
               m = bc_groups(i,k)
               tmp1(m) = 0
            enddo
c           Make space for the new variables we will add to the group
            nn = SUM(tmp1)       ! This is the number of variables we will be adding to the group
            do k = n_seq(i),1,-1
               bc_groups(i,k+nn) = bc_groups(i,k)
            enddo
            bc_groups(i,1:nn) = 0
            n_seq(i) = n_seq(i)+nn
c           Add the variables to the group preserving their orders
            nn = 1
            do j = 1,i-1
               do k = 1,n_seq(j)
                  m = bc_groups(j,k)
                  if ( tmp1(m)>0 ) then
                     bc_groups(i,nn) = m
                     tmp1(m) = 0
                     nn = nn+1
                  endif
               enddo
            enddo
c           Check the group
            if ( MINVAL(bc_groups(i,1:n_seq(i)))<1 ) then
               call pstop('create_bc_schedule',
     &            'Error inheriting dependencies')
            endif
         enddo
         deallocate ( tmp1 )
      endif

c     Deallocate existing 
      if ( associated(bc_grp) ) then
         do j = 1,size(bc_grp)
            if ( associated(bc_grp(j)%bc_seq) ) then
               deallocate ( bc_grp(j)%bc_seq )
            endif
         enddo
         deallocate ( bc_grp )
      endif

c     Create the output groups
      nbc_grp = N_groups
      allocate ( bc_grp(nbc_grp) )
      do i = 1,nbc_grp
         bc_grp(i)%nbc_seq = n_seq(i)
         allocate ( bc_grp(i)%bc_seq(n_seq(i),3) )
         do j = 1,n_seq(i)
            if ( bc_groups(i,j)<=n_var(1) ) then
               bc_grp(i)%bc_seq(j,1) = bc_groups(i,j)
               bc_grp(i)%bc_seq(j,2) = 0
            else if ( bc_groups(i,j)<=SUM(n_var(1:2)) ) then
               bc_grp(i)%bc_seq(j,1) = bc_groups(i,j)-n_var(1)
               bc_grp(i)%bc_seq(j,2) = 1
            else if ( bc_groups(i,j)<=SUM(n_var(1:3)) ) then
               bc_grp(i)%bc_seq(j,1) = -(bc_groups(i,j)-SUM(n_var(1:2)))
               bc_grp(i)%bc_seq(j,2) = 0
            else if ( bc_groups(i,j)<=SUM(n_var(1:4)) ) then
               bc_grp(i)%bc_seq(j,1) = -(bc_groups(i,j)-SUM(n_var(1:3)))
               bc_grp(i)%bc_seq(j,2) = 1
            else
               call pstop('create_bc_schedule','This should not occur')
            endif
            if ( used(bc_groups(i,j))==1 ) then
               bc_grp(i)%bc_seq(j,3) = 0
            else if ( used(bc_groups(i,j))==2 ) then
               bc_grp(i)%bc_seq(j,3) = 1
            else
               call pstop('create_bc_schedule','This should not occur')
            endif
         enddo
      enddo

      if (print_used) then
        write (*,*) 'BC sequence'
        do i = 1,nbc_grp
          write (*,*) 'Group=',i
          do j = 1,n_seq(i)
            call print_var_name(u,aux,bc_groups(i,j))
            if (bc_grp(i)%bc_seq(j,2) == 0) then
              write (*,*) '    Scalar=',bc_grp(i)%bc_seq(j,1)
            else
              write (*,*) '    Vector=',bc_grp(i)%bc_seq(j,1)
            endif
          enddo
          write (*,*)
        enddo
      endif

c     Free temporary memory
      deallocate ( used )
      deallocate ( bc_dep )
      deallocate ( dom_dep )
      deallocate ( tot_dep )
      deallocate ( n_seq )
      deallocate ( bc_groups )

      contains

c     check_dep_list
c     ################################################################
      integer function check_dep_list(N,dep_list)

      implicit none

      INTEGER :: dep_list(N,3)
      INTEGER i, j, N

      check_dep_list = 0
      i = 1
      do while ( dep_list(i,1).ne.0) 
         if ( dep_list(i,1)>0 ) then
            if ( dep_list(i,2).ne.0 ) then
               check_dep_list = 1       ! Invalid entry found
            endif
         else if ( dep_list(i,1)<0 ) then
            if ( dep_list(i,2)<0 .or. dep_list(i,2)<0 ) then
               check_dep_list = 1       ! Invalid entry found
            endif
         endif
         do j = 1,i-1
            if ( (dep_list(i,1)==dep_list(j,1)) .AND. 
     &           (dep_list(i,2)==dep_list(j,2)) ) then
               check_dep_list = 2       ! Duplicate entry found
            endif
         enddo
         i = i+1
         if (i > N) exit
      enddo

      end function check_dep_list

c     get_dep_index
c     ################################################################
      function get_dep_index(dep_var,n_var)

      implicit none

      INTEGER :: get_dep_index
      INTEGER, INTENT(IN) :: dep_var(3)
      INTEGER, INTENT(IN) :: n_var(4)
      get_dep_index = -1
      if ( dep_var(1)>0 .AND. dep_var(2)==0 ) then
c        Dependency is a scalar dependent variable
         get_dep_index = dep_var(1)
      else if ( dep_var(1)>0 .AND. dep_var(2)==1 ) then
c        Dependency is a vector dependent variable
         get_dep_index = n_var(1) + dep_var(1)
      else if ( dep_var(1)<0 .AND. dep_var(2)==0 ) then
c        Dependency is a vector dependent variable
         get_dep_index = n_var(1) + n_var(2) - dep_var(1)
      else if ( dep_var(1)<0 .AND. dep_var(2)==1 ) then
c        Dependency is a vector dependent variable
         get_dep_index = n_var(1) + n_var(2) + n_var(3) - dep_var(1)
      endif
      if ( get_dep_index<=0 ) then
c        Bad case
         call pstop('get_dep_index','Bad index')
      endif
      end function get_dep_index

c     print_var_name
c     ################################################################
      subroutine print_var_name(u,aux,k)

      use variable_setup, only: var_array, aux_array

      implicit none
      type(var_array) :: u
      type(aux_array) :: aux
      INTEGER :: k, n_var(4), N
      n_var(1) = u%nvar     ! Number of scalar dependent variables
      n_var(2) = 0          ! Number of vector dependent variables
      n_var(3) = aux%nvar   ! Number of scalar auxiliary variables
      n_var(4) = aux%nvec   ! Number of vector auxiliary variables
      N = n_var(1) + n_var(2) + n_var(3) + n_var(4)
      if ( k <= n_var(1) ) then
         write (*,*) '   ', u%array_var(k)%descr
      else if ( k <= SUM(n_var(1:2)) ) then
         !!!! Place holder for dependent vector
      else if ( k <= SUM(n_var(1:3)) ) then
         write (*,*) '   ', aux%var_list(k-SUM(n_var(1:2)))%descr
      else if ( k <= SUM(n_var(1:4)) ) then
         write (*,*) '   ', aux%vec_list(k-SUM(n_var(1:3)))%descr
      endif
      end subroutine print_var_name

#if defined(samrai)
c     pstop
c     ################################################################
      subroutine pstop(routine,message)

c     ---------------------------------------------------------------
c     Stops program at "routine" with "message"
c     ---------------------------------------------------------------

        implicit none

        character(*)  :: routine, message

c     Begin program

          write (*,*)
          write (*,*) trim(message)
          write (*,*) 'Program stopped at routine ',trim(routine)

        call abort

        stop

      end subroutine pstop
#endif

      end subroutine create_bc_schedule

      end module variables_bc

c module variables
c #####################################################################
      module variables

        use variables_bc

#if defined(samrai)
        integer :: bc_order = 1
#endif

      end module variables

#if defined(samrai)
c     applybc
c     #################################################################
      subroutine applybc(patch_var,nbc_seq,bc_seq)

c     -----------------------------------------------------------------
c     Applies the boundary conditions for SAMRAI
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer :: nbc_seq, bc_seq(nbc_seq,3)
      type(patch), TARGET :: patch_var

c     Local variables

      integer :: i
      type(var_array),pointer :: varray => null()

c     Interfaces

      INTERFACE
        subroutine applyDepVarBC(igr,bc_inf,varray,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          integer,optional :: order
        end subroutine applyDepVarBC
      END INTERFACE

      INTERFACE
        subroutine applyAuxVarBC(igr,bc_inf,varray,vaux,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine applyAuxVarBC
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use variable_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

c     Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

      call setAppBCs(varray,gv%aux)
cc      call allocAuxVariables(gv%aux)

c     Set boundary conditions

      !Impose BCs on group
      do i = 1,nbc_seq
        if (bc_seq(i,1)*AUX < 0) then !Process dep. variable
          call applyDepVarBC(1,bc_seq(i,:),varray,order=bc_order)
        elseif (bc_seq(i,1)*AUX > 0) then !Process aux. variable
          call applyAuxVarBC(1,bc_seq(i,:),varray,gv%aux
     .                      ,order=bc_order)
        endif
      enddo

c     Nullify aux pointers

cc      call deallocLocalAuxVars
      call destroyAppBCs

c     End program

      end subroutine applybc

c     postproc
c     #################################################################
      subroutine postproc(patch_var,nbc_seq,bc_seq)

c     -----------------------------------------------------------------
c     Applies the postprocessor for SAMRAI
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer :: nbc_seq, bc_seq(nbc_seq,3)
      type(patch), TARGET :: patch_var

c     Local variables

      integer :: i
      type(var_array),pointer :: varray => null()

c     Interfaces

      INTERFACE
        subroutine postproc_AuxVar(igr,bc_inf,varray,vaux)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine postproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use variable_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

c     Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions
      varray => gv%u_n
      u_0    => gv%u_0

      call setAppBCs(varray,gv%aux)
cc      call allocAuxVariables(gv%aux)

c     Postprocess auxiliary quantities

      do i = 1,nbc_seq
        if (bc_seq(i,1)*AUX > 0) then
          call postproc_AuxVar(1,bc_seq(i,:),varray,gv%aux)
        endif
      enddo

c     Nullify aux pointers

cc      call deallocLocalAuxVars
      call destroyAppBCs

c     End program

      end subroutine postproc

c     initializeAuxVar
c     #################################################################
      subroutine initializeAuxVar(patch_var,igrp)

c     -----------------------------------------------------------------
c     Fills auxiliary variable arrays
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      type(patch), TARGET :: patch_var
      integer:: igrp,i

c     Local variables

      type(var_array),pointer :: varray => null()

      INTERFACE
        subroutine preproc_AuxVar(igr,bc_inf,varray,vaux,order)
          use variable_setup
          integer    :: igr,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
          integer,optional :: order
        end subroutine preproc_AuxVar
      END INTERFACE

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use variable_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

c     Setup pointers

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      varray => gv%u_n
      u_0    => gv%u_0

      call setAppBCs(varray,gv%aux)
cc      call allocAuxVariables(gv%aux)

c     Fill auxiliary variable arrays

      do i = 1,gv%bc_grp(igrp)%nbc_seq
        if (gv%bc_grp(igrp)%bc_seq(i,1)*AUX > 0) then
          call preproc_AuxVar(1,gv%bc_grp(igrp)%bc_seq(i,:),varray
     .                       ,gv%aux,order=bc_order)
        endif
      enddo

c     Nullify pointers

cc      call deallocAuxVariables(gv%aux)
      call destroyAppBCs

      end subroutine initializeAuxVar

c     setupVarInitSeq
c     #################################################################
      subroutine setupVarInitSeq(patch_var,it)
c     -----------------------------------------------------------------
c     Sets up application BC for SAMRAI
c     -----------------------------------------------------------------

      use variables

      implicit none

c     Call variables

      integer :: it

      type(patch), TARGET :: patch_var

c     Local variables

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use variable_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c     Begin program

      gv => patch_var

c     Set application BCs (external)

      call setAppBCs(gv%u_n,gv%aux)

c     Set BC schedule

      call setBCsched(gv%u_n,gv%aux,it)

c     End program

      end subroutine setupVarInitSeq

#endif

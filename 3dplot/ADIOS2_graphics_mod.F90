
! module ADIOS2_graphics
! ##################################################################
      module ADIOS2_graphics

#if defined(ADIOS2)
        use graphics_io

        use var_io

        use ts_setup

        use adios2

#if !defined(ADIOS2)
        ! almost empty mod definition to be able to compile and link
        character(20) :: adios2_fname  ! set in ./plot/defineGraphics.F

#else
!       generic constants
        integer, parameter :: ARRAY_RANK = 3
        integer, parameter :: ATTR_RANK = 1
        integer, parameter :: VEXP_RANK = 1
        integer, parameter :: strlength = 20

        integer, parameter :: adios2_integer = 2 
        integer, parameter :: adios2_double = 6 
        integer, parameter :: adios2_string = 9 

        logical :: moving_mesh = .true.

!       ADIOS2 function call variables
        integer        :: num_open    ! first create, then append file
        character(20)  :: adios2_fname="adios2.bp" ! set in ./plot/defineGraphics.F
        integer        :: adios2_rank  ! ?=? my_rank ???

        logical,private :: adios2_debug=.false.

!       ADIOS2 no-xml definition variables
        integer*8      :: grp
        character(strlength)  :: grp_name = 'plot'

        integer, dimension(ATTR_RANK ) :: attr_dim = (/ARRAY_RANK/)
        integer, dimension(VEXP_RANK ) :: vexp_dim = (/VEXP_RANK/)
        integer, dimension(7)       :: array_dim_i 

        character(strlength), dimension(ARRAY_RANK) :: coords_descr = (/'X','Y','Z'/)
        character(strlength) :: cell_coords_gname = 'cells'
        character(strlength) :: node_coords_gname = 'nodes'
        character(strlength) :: attr_name = 'coords'

        character(strlength), dimension(ARRAY_RANK) :: attr_data

!       VisIt expression
        character(300) :: vexp

      contains 

!     initializeADIOS2Plot
!     ##############################################################
      subroutine initializeADIOS2Plot

      implicit none

!     Call variables
      
!     Local variables

      integer :: iattr, ierr

!     Begin program

      if (adios2_debug) write (*,'("%%%% initializeADIOS2Plot")') 

!     Setup plot dump variables

!     Attach location of node positions as attribute (relative path to timestep)

      do iattr=1,attr_dim(1)
        attr_data(iattr) = trim('/'//trim(node_coords_gname)//'/'//trim(coords_descr(iattr)))
      enddo

!     Check whether we are dealing with moving mesh
!c      moving_mesh = .not.checkMapDatabase()
      ! The mesh definitions with time below work only with moving_mesh=true 
      moving_mesh = .true.

      ! FIXME: isProc(0) is always equal to adios2_rank==0 ??
      !        then we do not need adios2_rank at all.
      call MPI_Comm_rank(adios2_world_comm,adios2_rank,mpierr)

!     Define adios group and transport method

      if (adios2_debug) &
         write (*,'("%%%%     declare adios group ",a)') grp_name 

      !call adios2_declare_group (grp,grp_name,"iter",1,ierr);

      if (adios2_debug) write (*,'("%%%%     select method")')

      !call adios2_select_method (grp, ADIOS2_METHOD, "","",ierr);

      num_open = 0

!     Exit if not in rank=0

      if (adios2_rank.ne.0) return

!     Create VISIT expressions (to identify vectors)

      if (adios2_debug) write (*,'("%%%%     selectExpressions")')  
      call selectExpressions

      if (adios2_debug) write (*,'("%%%% Done initializeADIOS2Plot")') 

!     End program

      end subroutine initializeADIOS2Plot

!     selectExpressions
!     #################################################################
      subroutine selectExpressions

      implicit none

!     -----------------------------------------------------------------
!     This routine identifies variables that are vector components,
!     and creates VISIT expressions that allows visit to recognize
!     them as vectors.
!     -----------------------------------------------------------------

      integer   :: igrp,ivar,ivp,ivm

      character(strlength) :: str1,str2,str3
      logical      :: first

      first = .true.
      vexp = ''

      do igrp=1,ngroups
        if (graph(igrp)%cartesian) then
          do ivar=1,ngraph
            ivm = max(ivar-1,1)
            ivp = min(ivar+1,ngraph)
            str1 = trim(graph(igrp)%array_graph(ivm )%vector_name)
            str2 = trim(graph(igrp)%array_graph(ivar)%vector_name)
            str3 = trim(graph(igrp)%array_graph(ivp )%vector_name)
            if (str2 /='' .and. (first.or.str1 /= str2)) then
              if (.not.first) then
                vexp = trim(vexp)//' ;'
              else
                first = .false.
              endif
              vexp=trim(vexp)//' '//trim(str2)//':{<'//trim(graph(igrp)%descr)//'/'//trim(graph(igrp)%array_graph(ivar)%descr)
            elseif (str2 /='' .and. str3 /= str2) then
              vexp= trim(vexp)//'>,<'//trim(graph(igrp)%descr)//'/'//trim(graph(igrp)%array_graph(ivar)%descr)//'>}'
            elseif (str2 /='') then
              vexp= trim(vexp)//'>,<'//trim(graph(igrp)%descr)//'/'//trim(graph(igrp)%array_graph(ivar)%descr)
            endif
          enddo
        endif
      enddo

!c      write (*,*) vexp

      end subroutine selectExpressions

!     writeADIOS2Plotfile
!     ##################################################################
      subroutine writeADIOS2Plotfile
!     ------------------------------------------------------------------
!     Write the output date using HDF5 HyperSlab Column
!     ------------------------------------------------------------------

      implicit none

!     Call variables

!     Local variables

      integer   :: i,j,k,ieq,icoords
      integer   :: slth
      integer   :: toffset
      integer*8 :: handle, totalsize
      integer*8 :: groupsize0  ! extra groupsize on rank 0
      integer   :: ierr
      integer*8 :: varid

      character*(strlength) :: attname, attdimstr
      character*128         :: basepath,path
      character*1           :: openmode

      integer :: iattr,igrp
      integer, dimension(ARRAY_RANK) :: count,node_count
      integer, dimension(ARRAY_RANK) :: offset
      character*128  :: node_ldim_str,cell_ldim_str,offset_str

      real(8),allocatable,dimension(:,:,:,:) :: cell,node

!     Variables for temporary fix of ADIOS2 bug
      integer :: nx_temp,ny_temp,nz_temp
      
!     Saved variables (calculated in first call)
      integer*8,save :: groupsize

      logical, save :: isfirst = .true.
      type(adios2_engine), save :: engine
      type(adios2_io), save :: io
      type(adios2_variable) :: var
      type(adios2_attribute) :: attr
      integer :: err, istatus

!     Begin program

      if (adios2_debug) write (*,'("%%%% writeADIOS2Plot")') 

!     Create the output file and data groups

      toffset = 0
      if (itime == -1) toffset = 1 

      if (num_open == 0) then
          openmode = 'w'
      else
          openmode = 'a'
      endif

      slth = len(trim(vexp))

      if (num_open == 0) then
         if (my_rank.eq.0) print *, 'ADIOS2: Writing restart_3db from ', trim(adios2_fname)
         call adios2_declare_io(io, adios2obj, "pixplot", ierr)
         
         ! Define global dimension variables a the first time
         call adios2_define_variable (var, io, "nxd", adios2_type_integer4, ierr)
         call adios2_define_variable (var, io, "nyd", adios2_type_integer4, ierr)
         call adios2_define_variable (var, io, "nzd", adios2_type_integer4, ierr)
         call adios2_define_variable (var, io, "nxd+1", adios2_type_integer4, ierr)
         call adios2_define_variable (var, io, "nyd+1", adios2_type_integer4, ierr)
         call adios2_define_variable (var, io, "nzd+1", adios2_type_integer4, ierr)
         
         if (adios2_rank.eq.0) then
            call adios2_define_attribute (attr, io, "/schema/name", "Pixie", ierr)
            if (adios2_debug) &
                 write (*,'("rank",i3,"%%%%    wrote attribute schema")') my_rank
            
            ! time array
            call adios2_define_variable (var, io, "time", adios2_type_dp, ierr) 
            call adios2_define_variable (var, io, "itime", adios2_type_integer4, ierr) 
            
            if (slth > 0) then
               call adios2_define_variable (var, io, "visit_expressions", adios2_type_string, ierr)
            endif
         endif
      endif

      basepath = ""

      write (offset_str   , '(i0,",",i0,",",i0)') ilomg,jlomg,klomg
      write (cell_ldim_str, '(i0,",",i0,",",i0,",iter")') nxl,nyl,nzl
      !c      write (node_ldim_str, '(i0,",",i0,",",i0)') nxl+1,nyl+1,nzl+1
      ! Temporary fix to ADIOS bug (L. Chacon, 10/6/10)
      ! Don't know if it is still a problem in Adios2 (J. Choi, 3/6/20)
      nx_temp = nxl
      if (isBdry(gv%gparams,nxl,1,2)) nx_temp = nxl + 1
      ny_temp = nyl
      if (isBdry(gv%gparams,nyl,1,4)) ny_temp = nyl + 1
      nz_temp = nzl
      if (isBdry(gv%gparams,nzl,1,6)) nz_temp = nzl + 1
      
      write (node_ldim_str, '(i0,",",i0,",",i0,",iter")') nx_temp,ny_temp,nz_temp

!$$$      write (*,*) isBdry(gv%gparams,nxl,1,2),nx_temp,ihip-ilo+1
!$$$      write (*,*) isBdry(gv%gparams,nyl,1,4),ny_temp,jhip-jlo+1
!$$$      write (*,*) isBdry(gv%gparams,nzl,1,6),nz_temp,khip-klo+1
! End of temporary fix

      if (adios2_debug) &
           write (*,'("rank",i3,"%% offset=",a," cell=",a," node=",a)') &
           my_rank, trim(offset_str), trim(cell_ldim_str), &
           trim(node_ldim_str)
      
      ! Define coordinate variables and calculate size 
      if (num_open == 0) then
          ! define cells/X, cells/Y, cells/Z,
          path="/"//trim(cell_coords_gname) 
          do icoords=1,ARRAY_RANK
             ! cell_ldim_str, "nxd,nyd,nzd", offset_str
             print *, my_rank, '#10:def:', trim(path)//'/'//trim(coords_descr(icoords))
             call adios2_define_variable (var, io, trim(path)//'/'//trim(coords_descr(icoords)), &
                  adios2_type_dp, 3, &
                  int((/ nxd, nyd, nzd /), kind=8), &
                  int((/ ilomg, jlomg, klomg /), kind=8), &
                  int((/ nxl, nyl, nzl /), kind=8), .true., ierr)
             
          enddo

          ! define nodes/X, nodes/Y, nodes/Z,
          path="/"//trim(node_coords_gname) 
          do icoords=1,ARRAY_RANK
             ! node_ldim_str, "nxd+1,nyd+1,nzd+1", offset_str
             print *, my_rank, '#20:def:', trim(path)//'/'//trim(coords_descr(icoords))
             call adios2_define_variable (var, io, trim(path)//'/'//trim(coords_descr(icoords)), &
                  adios2_type_dp, 3, &
                  int((/ nxd+1, nyd+1, nzd+1 /), kind=8), &
                  int((/ ilomg, jlomg, klomg /), kind=8), &
                  int((/ nx_temp,ny_temp,nz_temp /), kind=8), .true., ierr)
          enddo

          ! Define data variables and calculate size
          do igrp=1,ngroups
              ! set path in file for all variables
              basepath = "/"//graph(igrp)%descr
              if (graph(igrp)%cartesian) then
                 ! Define Coordinate Attribute data 
                 !write (attdimstr, '(i0)') attr_dim(1)
                 call adios2_define_attribute (attr, io, &
                      trim(basepath)//"/"//trim(attr_name)//'/'//"ncoord", &
                      attr_dim(1), ierr)

                 do iattr=1,attr_dim(1)
                    write (attname,'("coord",i0)') iattr
                    call adios2_define_attribute (attr, io, &
                         trim(basepath)//"/"//trim(attr_name)//'/'//trim(attname), &
                         attr_data(iattr), ierr)
                 enddo

                 if (adios2_debug) write (*,'("%%%%     defined coordinate attributes")') 

                 ! Define nodal data variables    
                 do ieq=1,nqty(igrp)
                    ! node_ldim_str, "nxd+1,nyd+1,nzd+1", offset_str
                    print *, my_rank, '#30:def:', trim(basepath)//'/'//trim(graph(igrp)%array_graph(ieq)%descr)
                    call adios2_define_variable (var, io, &
                         trim(basepath)//'/'//trim(graph(igrp)%array_graph(ieq)%descr), &
                         adios2_type_dp, 3, &
                         int((/ nxd+1, nyd+1, nzd+1 /), kind=8), &
                         int((/ ilomg, jlomg, klomg /), kind=8), &
                         int((/ nx_temp,ny_temp,nz_temp /), kind=8), .true., ierr)
                 enddo
              else
                 ! Define cell data variables
                 do ieq=1,nqty(igrp)
                    ! cell_ldim_str, "nxd,nyd,nzd", offset_str
                    print *, my_rank, '#40:def:', igrp, nqty(igrp), ieq, &
                         trim(basepath)//'/'//trim(graph(igrp)%array_graph(ieq)%descr)
                    call adios2_define_variable (var, io, &
                         trim(basepath)//'/'//trim(graph(igrp)%array_graph(ieq)%descr), &
                         adios2_type_dp, 3, &
                         int((/ nxd, nyd, nzd /), kind=8), &
                         int((/ ilomg, jlomg, klomg /), kind=8), &
                         int((/ nxl, nyl, nzl /), kind=8), .true., ierr)
                 enddo
              endif
           enddo

           ! open only once
           print *, 'adios2_fname', adios2_fname
           call adios2_open(engine, io, adios2_fname, adios2_mode_write, adios2_world_comm, ierr)
      endif  ! first open

      call adios2_begin_step(engine, adios2_step_mode_append, 0.0, istatus, ierr)

      
!     Write the variables 
!       (they are buffered and I/O starts at adios2_close())

      if (num_open == 0) then
          ! First writing
          !if (isProc(0)) then
          if (adios2_rank.eq.0) then
              if (slth > 0) then
                  ! Write VisIt expressions from proc 0
                  call adios2_put(engine, "visit_expressions", vexp, ierr)
                  if (adios2_debug) &
                       write (*,'("rank",i3,"%%%%    wrote vexp l=",i0)') my_rank, slth
              endif
          endif
      endif  ! first writing

      ! Write dimension variables (must be each time in adios)
      call adios2_put(engine,"nxd",nxd,ierr)
      call adios2_put(engine,"nyd",nyd,ierr)
      call adios2_put(engine,"nzd",nzd,ierr)
      call adios2_put(engine,"nxd+1",nxd+1,ierr)
      call adios2_put(engine,"nyd+1",nyd+1,ierr)
      call adios2_put(engine,"nzd+1",nzd+1,ierr)
      !if (isProc(0)) then
      if (adios2_rank.eq.0) then
          call adios2_put(engine,"time",time,ierr)
          call adios2_put(engine,"itime",itime+toffset,ierr)
      endif
      if (adios2_debug) &
           write (*,'("rank",i3,"%%%%     wrote dimensions")') my_rank

      ! First writing or always if the mesh is changing
      if (moving_mesh .or. num_open == 0) then
         ! Write the coordinate data

          allocate (cell(ilom:ihip,jlom:jhip,klom:khip,3) &
                   ,node(ilom:ihip,jlom:jhip,klom:khip,3))

          ! Cell positions
          do k=klom,khip
            do j=jlom,jhip
              do i=ilom,ihip
                 call getCartesianCoordinates(gv%gparams,i,j,k,iggx,iggy &
                      ,iggz,iig,jjg,kkg,cell(i,j,k,1),cell(i,j,k,2) &
                      ,cell(i,j,k,3))
              enddo
            enddo
          enddo

          ! Define coordinate variables and calculate size

          path="/"//trim(cell_coords_gname) 
          do icoords=1,ARRAY_RANK
             call adios2_put(engine, &
                  trim(path)//'/'//trim(coords_descr(icoords)), &
                  cell(ilo:ihi,jlo:jhi,klo:khi,icoords),adios2_mode_sync,ierr)
          enddo
          if (adios2_debug) &
               write (*,'("rank",i3,"%%%%     wrote cell positions")') my_rank

          ! Node positions
          call interpolateNodalData

          path="/"//trim(node_coords_gname) 
          do icoords=1,ARRAY_RANK
             call adios2_put(engine, &
                  trim(path)//'/'//trim(coords_descr(icoords)), &
                  node(ilo:nx_temp,jlo:ny_temp,klo:nz_temp,icoords),adios2_mode_sync,ierr)
          enddo
          if (adios2_debug) &
               write (*,'("rank",i3,"%%%%     wrote node positions")') my_rank

          deallocate(cell,node)

      endif ! first or all writes (coordinate data)

!     Write the actual data

      do igrp=1,ngroups

!         set path in file for all variables
          basepath = "/"//graph(igrp)%descr

          if (graph(igrp)%cartesian) then

!             Write nodal data
              allocate (cell(ilom:ihip,jlom:jhip,klom:khip,1) &
                       ,node(ilom:ihip,jlom:jhip,klom:khip,1))

              do ieq=1,nqty(igrp)
                  ! Interpolate nodal values
                  cell(:,:,:,1) = graph(igrp)%array_graph(ieq)%array
                  call interpolateNodalData   ! calculates node array
                
                  ! write nodal data
                  call adios2_put(engine, &
                       trim(basepath)//"/"//trim(graph(igrp)%array_graph(ieq)%descr), &
                       node(ilo:nx_temp,jlo:ny_temp,klo:nz_temp,1),adios2_mode_sync,ierr)

              enddo
              if (adios2_debug) &
                   write (*,'("%%%%     wrote nodal data")') 

              deallocate(cell,node)

          else ! not Cartesian

!             Write cell-centered data
              do ieq=1,nqty(igrp)
                  call adios2_put(engine, &
                       trim(basepath)//"/"//trim(graph(igrp)%array_graph(ieq)%descr), &
                       graph(igrp)%array_graph(ieq)%array(ilo:ihi,jlo:jhi,klo:khi),adios2_mode_sync,ierr)
              enddo
              if (adios2_debug) &
                   write (*,'("%%%%     wrote cell data")') 

          endif

      enddo ! 1..ngroups

!     Release resources

      call adios2_end_step(engine, ierr)
      !call adios2_close(engine, ierr)
      num_open = num_open+1

      if (adios2_debug) write (*,'("%%%% Done writeADIOS2Plot")') 

!     End program

      contains

!     interpolateNodalData
!     ####################################################################
      subroutine interpolateNodalData

!     ---------------------------------------------------------------------
!     Interpolates cell-centered data to nodal data
!     ---------------------------------------------------------------------

        implicit none

        do k=klo,khip
          do j=jlo,jhip
            do i=ilo,ihip
!c              if (isSP2(i,1,ibc=1)) then
!c                node(i,j,k,:)=
!c                   0.5*(sum(cell(i,jlo:jhi,k-1,:),dim=1)/(jhi-jlo+1)
!c                       +sum(cell(i,jlo:jhi,k  ,:),dim=1)/(jhi-jlo+1))
!c              else
                node(i,j,k,:) = 0.125 &
                              *(cell(i  ,j  ,k  ,:)+cell(i-1,j  ,k  ,:) &
                              +cell(i  ,j-1,k  ,:)+cell(i-1,j-1,k  ,:) &
                              +cell(i  ,j  ,k-1,:)+cell(i-1,j  ,k-1,:) &
                              +cell(i  ,j-1,k-1,:)+cell(i-1,j-1,k-1,:))
!c              endif
            enddo
          enddo
        enddo

      end subroutine interpolateNodalData

      end subroutine writeADIOS2Plotfile

#endif
#endif
      end module ADIOS2_graphics

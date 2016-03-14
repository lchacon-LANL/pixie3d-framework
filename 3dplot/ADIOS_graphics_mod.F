
c module ADIOS_graphics
c ##################################################################
      module ADIOS_graphics

        use graphics_io

        use var_io

        use ts_setup

#if !defined(adios)
        ! almost empty mod definition to be able to compile and link
        character(20) :: adios_fname  ! set in ./plot/defineGraphics.F

#else
c       generic constants
        integer, parameter :: ARRAY_RANK = 3
        integer, parameter :: ATTR_RANK = 1
        integer, parameter :: VEXP_RANK = 1
        integer, parameter :: strlength = 20

        integer, parameter :: adios_integer = 2 
        integer, parameter :: adios_double = 6 
        integer, parameter :: adios_string = 9 

        logical :: moving_mesh = .true.

c       ADIOS function call variables
        integer        :: num_open    ! first create, then append file
        character(20)  :: adios_fname="adios.bp" ! set in ./plot/defineGraphics.F
        integer        :: adios_rank  ! ?=? my_rank ???

        logical,private :: adios_debug=.false.

c       ADIOS no-xml definition variables
        integer*8      :: grp
        character(strlength)  :: grp_name = 'plot'

        integer, dimension(ARRAY_RANK) ::
     .                         array_dim,node_array_dim
        integer, dimension(ATTR_RANK ) ::
     .                         attr_dim = (/ARRAY_RANK/)
        integer, dimension(VEXP_RANK ) ::
     .                         vexp_dim = (/VEXP_RANK/)
        integer, dimension(7)       :: array_dim_i 

        character(strlength), dimension(ARRAY_RANK) :: coords_descr
     .                                            = (/'X','Y','Z'/)
        character(strlength) :: cell_coords_gname = 'cells'
        character(strlength) :: node_coords_gname = 'nodes'
        character(strlength) :: attr_name = 'coords'

        character(strlength), dimension(ARRAY_RANK) :: attr_data

c       VisIt expression
        character(300) :: vexp

      contains 

c     initializeADIOSPlot
c     ##############################################################
      subroutine initializeADIOSPlot

      implicit none

c     Call variables
      
c     Local variables

      integer :: iattr, ierr

c     Begin program

      if (adios_debug) write (*,'("%%%% initializeADIOSPlot")') 

c     Setup plot dump variables

      array_dim      = (/nxd,nyd,nzd/)
      node_array_dim = (/nxd+1,nyd+1,nzd+1/)

c     Attach location of node positions as attribute (relative path to timestep)

      do iattr=1,attr_dim(1)
        attr_data(iattr) = trim('/'//trim(node_coords_gname)//
     $                     '/'//trim(coords_descr(iattr)))
      enddo

c     Check whether we are dealing with moving mesh
cc      moving_mesh = .not.checkMapDatabase()
      ! The mesh definitions with time below work only with moving_mesh=true 
      moving_mesh = .true.

      ! FIXME: isProc(0) is always equal to adios_rank==0 ??
      !        then we do not need adios_rank at all.
      call MPI_Comm_rank(adios_world_comm,adios_rank,mpierr)


c     Define adios group and transport method
      if (adios_debug)
     .   write (*,'("%%%%     declare adios group ",a)') grp_name 
      !! FIXME: after adios 1.2.1, adios_declare_group has one more
       ! argument: (ADIOS_FLAG stats), 1 for true, 2 for false
      if (ADIOS_VER <= 1.21) then
        call adios_declare_group (grp,grp_name,"iter",ierr);
      else
        call adios_declare_group (grp,grp_name,"iter",1,ierr);
      endif

      if (adios_debug) write (*,'("%%%%     select method MPI")')

      call adios_select_method (grp, "MPI", "","",ierr);

      num_open = 0

c     Exit if not in rank=0
      !if (.not.isProc(0)) return
      if (adios_rank.ne.0) return

c     Create VISIT expressions (to identify vectors)
      if (adios_debug) write (*,'("%%%%     selectExpressions")')  
      call selectExpressions

      if (adios_debug) write (*,'("%%%% Done initializeADIOSPlot")') 

c     End program
      end subroutine initializeADIOSPlot

c     selectExpressions
c     #################################################################
      subroutine selectExpressions

      implicit none

c     -----------------------------------------------------------------
c     This routine identifies variables that are vector components,
c     and creates VISIT expressions that allows visit to recognize
c     them as vectors.
c     -----------------------------------------------------------------

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
              vexp=trim(vexp)//' '//trim(str2)//':{<'
     .           //trim(graph(igrp)%descr)
     .           //'/'//trim(graph(igrp)%array_graph(ivar)%descr)
            elseif (str2 /='' .and. str3 /= str2) then
              vexp= trim(vexp)//'>,<'//trim(graph(igrp)%descr)//
     .           '/'//trim(graph(igrp)%array_graph(ivar)%descr)//'>}'
            elseif (str2 /='') then
              vexp= trim(vexp)//'>,<'//trim(graph(igrp)%descr)//
     .           '/'//trim(graph(igrp)%array_graph(ivar)%descr)
            endif
          enddo
        endif
      enddo

cc      write (*,*) vexp

      end subroutine selectExpressions

c     writeADIOSPlotfile
c     ##################################################################
      subroutine writeADIOSPlotfile
c     ------------------------------------------------------------------
c     Write the output date using HDF5 HyperSlab Column
c     ------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer        :: i,j,k,ieq,icoords
      integer        :: slth
      integer        :: toffset
      integer*8      :: handle, totalsize
      integer*8      :: groupsize0  ! extra groupsize on rank 0
      integer        :: ierr

      character*(strlength) :: attname, attdimstr
      character*128         :: basepath,path
      character*1           :: openmode

      integer :: iattr,igrp
      integer, dimension(ARRAY_RANK) :: count,node_count
      integer, dimension(ARRAY_RANK) :: offset
      character*128  :: node_ldim_str, cell_ldim_str, offset_str

      real(8),allocatable,dimension(:,:,:,:) :: cell,node

c     Variables for temporary fix of ADIOS bug
      integer :: nx_temp,ny_temp,nz_temp
      
c     Saved variables (calculated in first call)
      integer*8,save :: groupsize

c     Begin program

      if (adios_debug) write (*,'("%%%% writeADIOSPlot")') 

c     Create the output file and data groups

      toffset = 0
      if (itime == -1) toffset = 1 

      if (num_open == 0) then
          openmode = 'w'
      else
          openmode = 'a'
      endif

      slth = len(trim(vexp))

      call adios_open (handle, grp_name, adios_fname, openmode,
     .                 adios_world_comm, ierr)

      groupsize0 = 0
c     Define variables and
c     Calculate/approx. the size of data we will going to write out
      if (num_open == 0) then
          groupsize = 0
          ! Define global dimension variables a the first time
          call adios_define_var (grp,"nxd","/dimension",
     .                           adios_integer,"","","",ierr)
          call adios_define_var (grp,"nyd","/dimension",
     .                           adios_integer,"","","",ierr)
          call adios_define_var (grp,"nzd","/dimension",
     .                           adios_integer,"","","",ierr)
          call adios_define_var (grp,"nxd+1","/dimension",
     .                           adios_integer,"","","",ierr)
          call adios_define_var (grp,"nyd+1","/dimension",
     .                           adios_integer,"","","",ierr)
          call adios_define_var (grp,"nzd+1","/dimension",
     .                           adios_integer,"","","",ierr)
          groupsize = groupsize + 12*4  ! aux/ and dimension/ variables
 
          !if (isProc(0)) then
          if (adios_rank.eq.0) then
              call adios_define_attribute (grp, "name","/schema",
     .                                 adios_string,"Pixie","",ierr)
              if (adios_debug)
     .          write (*,'("rank",i3,"%%%%    wrote attribute schema")') 
     .             my_rank

              ! time array
              call adios_define_var (grp, "time", "",
     .               adios_double,"iter","","", ierr) 
              call adios_define_var (grp, "itime", "",
     .               adios_integer,"iter","","", ierr) 

              if (slth > 0) then
                  call adios_define_var (grp, "visit_expressions","",
     .                                   adios_string,"","","",ierr)
                  groupsize0 = groupsize0 + len(trim(vexp))+1 ! vexp
              endif
          endif
      endif

      basepath = ""

      write (offset_str   , '(i0,",",i0,",",i0)') ilomg,jlomg,klomg
      write (cell_ldim_str, '(i0,",",i0,",",i0,",iter")') nxl,nyl,nzl
cc      write (node_ldim_str, '(i0,",",i0,",",i0)') nxl+1,nyl+1,nzl+1
c Temporary fix to ADIOS bug (L. Chacon, 10/6/10)
      nx_temp = nxl
      if (isBdry(gv%gparams,nxl,1,2)) nx_temp = nxl + 1
      ny_temp = nyl
      if (isBdry(gv%gparams,nyl,1,4)) ny_temp = nyl + 1
      nz_temp = nzl
      if (isBdry(gv%gparams,nzl,1,6)) nz_temp = nzl + 1

      write (node_ldim_str, '(i0,",",i0,",",i0,",iter")')
     .        nx_temp,ny_temp,nz_temp

cc      write (*,*) isBdry(nxl,1,2),nx_temp,ihip-ilo+1
cc      write (*,*) isBdry(nyl,1,4),ny_temp,jhip-jlo+1
cc      write (*,*) isBdry(nzl,1,6),nz_temp,khip-klo+1
c End of temporary fix

      if (adios_debug)
     .  write (*,'("rank",i3,"%% offset=",a," cell=",a," node=",a)') 
     .       my_rank, trim(offset_str), trim(cell_ldim_str), 
     .       trim(node_ldim_str)

      ! Define coordinate variables and calculate size 
      if (num_open == 0) then
          ! define cells/X, cells/Y, cells/Z,
          path="/"//trim(cell_coords_gname) 
          do icoords=1,ARRAY_RANK
              call adios_define_var (grp, coords_descr(icoords), 
     .                 path,
     .                 adios_double, 
     .                 cell_ldim_str, "nxd,nyd,nzd",
     .                 offset_str, ierr)
          enddo
          groupsize = groupsize 
     .      + (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1) *ARRAY_RANK*8 ! cell coords

          ! define nodes/X, nodes/Y, nodes/Z,
          path="/"//trim(node_coords_gname) 
          do icoords=1,ARRAY_RANK
              call adios_define_var (grp, coords_descr(icoords), 
     .               path,
     .               adios_double, 
     .               node_ldim_str, "nxd+1,nyd+1,nzd+1", 
     .               offset_str, ierr)
          enddo
          groupsize = groupsize 
cc     .      + (ihip-ilo+1) * (jhip-jlo+1) * (khip-klo+1) *ARRAY_RANK*8 ! node coords
c Temporary fix to ADIOS bug (L. Chacon 10/6/10)
     .      + nx_temp * ny_temp * nz_temp *ARRAY_RANK*8 ! node coords
c End of temporary fix

          ! Define data variables and calculate size
          do igrp=1,ngroups
              ! set path in file for all variables
              basepath = "/"//graph(igrp)%descr
              if (graph(igrp)%cartesian) then
                  ! Define Coordinate Attribute data 
                  write (attdimstr, '(i0)') attr_dim(1)
                  call adios_define_attribute (grp, "ncoord",
     .                  trim(basepath)//"/"//trim(attr_name),
     .                  adios_integer,attdimstr,"",ierr)
                  do iattr=1,attr_dim(1)
                      write (attname,'("coord",i0)') iattr
                      call adios_define_attribute (grp, attname,
     .                     trim(basepath)//"/"//trim(attr_name),
     .                     adios_string,attr_data(iattr),"",ierr)
                  enddo

                  if (adios_debug) write (*,
     .              '("%%%%     defined coordinate attributes")') 

                  ! Define nodal data variables    
                  do ieq=1,nqty(igrp)
                     call adios_define_var (grp, 
     .                        graph(igrp)%array_graph(ieq)%descr, 
     .                        basepath, adios_double, 
     .                        node_ldim_str, "nxd+1,nyd+1,nzd+1",
     .                        offset_str, ierr)
                  enddo

                  groupsize = groupsize +
cc     .             nqty(igrp)*(ihip-ilo+1)*(jhip-jlo+1)*(khip-klo+1)*8 ! node data
c Temporary fix to ADIOS bug (L. Chacon 10/6/10)
     .              nqty(igrp)*nx_temp*ny_temp*nz_temp*8 ! node data
c End of temporary fix
            else
                  ! Define cell data variables
                  do ieq=1,nqty(igrp)
                     call adios_define_var (grp, 
     .                        graph(igrp)%array_graph(ieq)%descr, 
     .                        basepath, adios_double, 
     .                        cell_ldim_str, "nxd,nyd,nzd", 
     .                        offset_str, ierr)
                  enddo

                  groupsize = groupsize +
     .              nqty(igrp)*(ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*8 ! cell data
              endif
          enddo

      endif  ! first open

      ! Tell ADIOS about the size of all data in this process
      call adios_group_size (handle, groupsize+groupsize0, 
     .                       totalsize, ierr)
      if (adios_debug)
     .  write (*,'("rank",i3," %%%%  groupsize=",i0,"  totalsize=",i0)') 
     .    my_rank, groupsize+groupsize0, totalsize


c     Write the variables 
c       (they are buffered and I/O starts at adios_close())

      if (num_open == 0) then
          ! First writing
          !if (isProc(0)) then
          if (adios_rank.eq.0) then
              if (slth > 0) then
                  ! Write VisIt expressions from proc 0
                  call adios_write(handle,"visit_expressions",vexp,ierr)
                  if (adios_debug)
     .              write (*,'("rank",i3,"%%%%    wrote vexp l=",i0)') 
     .                    my_rank, slth
              endif
          endif
      endif  ! first writing

      ! Write dimension variables (must be each time in adios)
      call adios_write(handle,"nxd",nxd,ierr)
      call adios_write(handle,"nyd",nyd,ierr)
      call adios_write(handle,"nzd",nzd,ierr)
      call adios_write(handle,"nxd+1",nxd+1,ierr)
      call adios_write(handle,"nyd+1",nyd+1,ierr)
      call adios_write(handle,"nzd+1",nzd+1,ierr)
      !if (isProc(0)) then
      if (adios_rank.eq.0) then
          call adios_write(handle,"time",time,ierr)
          call adios_write(handle,"itime",itime+toffset,ierr)
      endif
      if (adios_debug)
     .  write (*,'("rank",i3,"%%%%     wrote dimensions")') my_rank


      ! First writing or always if the mesh is changing
      if (moving_mesh .or. num_open == 0) then

c         Write the coordinate data

          allocate (cell(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,node(ilom:ihip,jlom:jhip,klom:khip,3))

          ! Cell positions
          do k=klom,khip
            do j=jlom,jhip
              do i=ilom,ihip
                call getCartesianCoordinates(gv%gparams,i,j,k,iggx,iggy
     .               ,iggz,iig,jjg,kkg,cell(i,j,k,1),cell(i,j,k,2)
     .               ,cell(i,j,k,3))
              enddo
            enddo
          enddo

          ! Define coordinate variables and calculate size

          path="/"//trim(cell_coords_gname) 
          do icoords=1,ARRAY_RANK
              call adios_write(handle, 
     .                trim(path)//"/"//coords_descr(icoords),
     .                cell(ilo:ihi,jlo:jhi,klo:khi,icoords),ierr)
          enddo
          if (adios_debug)
     .      write (*,'("rank",i3,"%%%%     wrote cell positions")') 
     .         my_rank

          ! Node positions
          call interpolateNodalData

          path="/"//trim(node_coords_gname) 
          do icoords=1,ARRAY_RANK
              call adios_write(handle, 
     .            trim(path)//"/"//coords_descr(icoords),
cc     .            node(ilo:ihip,jlo:jhip,klo:khip,icoords), ierr)
c Temporary fix to ADIOS bug (L. Chacon 10/6/10)
     .           node(ilo:nx_temp,jlo:ny_temp,klo:nz_temp,icoords),ierr)
c End of temporary fix
          enddo
          if (adios_debug)
     .      write (*,'("rank",i3,"%%%%     wrote node positions")') 
     .          my_rank

          deallocate(cell,node)

      endif ! first or all writes (coordinate data)

c     Write the actual data

      do igrp=1,ngroups

c         set path in file for all variables
          basepath = "/"//graph(igrp)%descr

          if (graph(igrp)%cartesian) then

c             Write nodal data
              allocate (cell(ilom:ihip,jlom:jhip,klom:khip,1)
     .                 ,node(ilom:ihip,jlom:jhip,klom:khip,1))

              do ieq=1,nqty(igrp)
                  ! Interpolate nodal values
                  cell(:,:,:,1) = graph(igrp)%array_graph(ieq)%array
                  call interpolateNodalData   ! calculates node array
                
                  ! write nodal data
                  call adios_write(handle, 
     .                  trim(basepath)//"/"//
     .                      graph(igrp)%array_graph(ieq)%descr,
cc     .                  node(ilo:ihip,jlo:jhip,klo:khip,1), ierr)
c Temporary fix to ADIOS bug (L. Chacon 10/6/10)
     .           node(ilo:nx_temp,jlo:ny_temp,klo:nz_temp,1),ierr)
c End of temporary fix

              enddo
              if (adios_debug)
     .            write (*,'("%%%%     wrote nodal data")') 

              deallocate(cell,node)

          else ! not Cartesian

c             Write cell-centered data
              do ieq=1,nqty(igrp)
                  call adios_write(handle, 
     .                  trim(basepath)//"/"//
     .                      graph(igrp)%array_graph(ieq)%descr,
     .                  graph(igrp)%array_graph(ieq)
     .                     %array(ilo:ihi,jlo:jhi,klo:khi), ierr)
              enddo
              if (adios_debug)
     .            write (*,'("%%%%     wrote cell data")') 

          endif

      enddo ! 1..ngroups

c     Release resources

      call adios_close(handle, ierr)
      num_open = num_open+1

      if (adios_debug) write (*,'("%%%% Done writeADIOSPlot")') 

c     End program

      contains

c     interpolateNodalData
c     ####################################################################
      subroutine interpolateNodalData

c     ---------------------------------------------------------------------
c     Interpolates cell-centered data to nodal data
c     ---------------------------------------------------------------------

        implicit none

        do k=klo,khip
          do j=jlo,jhip
            do i=ilo,ihip
cc              if (isSP2(i,1,ibc=1)) then
cc                node(i,j,k,:)=
cc     .             0.5*(sum(cell(i,jlo:jhi,k-1,:),dim=1)/(jhi-jlo+1)
cc     .                 +sum(cell(i,jlo:jhi,k  ,:),dim=1)/(jhi-jlo+1))
cc              else
                node(i,j,k,:) = 0.125
     .                       *(cell(i  ,j  ,k  ,:)+cell(i-1,j  ,k  ,:)
     .                        +cell(i  ,j-1,k  ,:)+cell(i-1,j-1,k  ,:)
     .                        +cell(i  ,j  ,k-1,:)+cell(i-1,j  ,k-1,:)
     .                        +cell(i  ,j-1,k-1,:)+cell(i-1,j-1,k-1,:))
cc              endif
            enddo
          enddo
        enddo

      end subroutine interpolateNodalData

      end subroutine writeADIOSPlotfile

#endif

      end module ADIOS_graphics

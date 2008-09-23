
c module HDF5_graphics
c ##################################################################
      module HDF5_graphics

        use graphics_io

        use variables

        use ts_setup

        use generalPurposeFunctions

#if !defined(hdf5)

        character(20) :: hdf5_fname

#else
        use HDF5

        integer        :: toffset=0

        integer(hsize_t) :: izero=0

c       HDF5 ID variables
        integer(hid_t) :: group_id,subgroup_id
        integer(hid_t) :: file_id        
        integer(hid_t) :: dataset_id     
        integer(hid_t) :: datatype_id     
        integer(hid_t) :: dataspace_id,node_dataspace_id
        integer(hid_t) :: attribute_id
        integer(hid_t) :: dataspace_attr_id
        integer(hid_t) :: attr_type_id

c       variables for parallel HDF5
        integer(hid_t) :: plist_id
        integer(hid_t) :: hs_dataspace_id
        integer(hid_t) :: p_dataspace_id,p_node_dataspace_id

        integer        :: hdferr,igrp

c       generic HDF5 variables
        integer   , parameter :: ARRAY_RANK = 3
        integer   , parameter :: ATTR_RANK = 1
        integer   , parameter :: VEXP_RANK = 1
        integer(hsize_t), parameter :: strlength = 20

        integer(hsize_t), dimension(ARRAY_RANK) ::
     .                         array_dim,node_array_dim
        integer(hsize_t), dimension(ATTR_RANK ) ::
     .                         attr_dim = (/ARRAY_RANK/)
        integer(hsize_t), dimension(VEXP_RANK ) ::
     .                         vexp_dim = (/VEXP_RANK/)
        integer(hsize_t), dimension(7)       :: array_dim_i 

        character(20), dimension(ARRAY_RANK) :: coords_descr
     .                                            = (/'X','Y','Z'/)
        character(20)                        :: cell_coords_gname
     .                                            = 'cell_coords'
        character(20)                        :: node_coords_gname
     .                                            = 'node_coords'
        character(20)                        :: attr_name = 'coords'

        character(strlength), dimension(ARRAY_RANK) :: attr_data

        character(20) :: hdf5_fname

      contains 

c     initializeHDF5
c     ##############################################################
      subroutine initializeHDF5

      implicit none

c     Call variables
      
c     Local variables

      integer(hsize_t)    :: iattr,slth
      character(300) :: vexp

c     Begin program

c     Setup plot dump variables

      array_dim      = (/nxd,nyd,nzd/)
      node_array_dim = (/nxd+1,nyd+1,nzd+1/)

      !Attach location of node positions as attribute
      do iattr=1,attr_dim(1)
        attr_data(iattr) = trim('/'//trim(node_coords_gname)//
     $                     '/'//trim(coords_descr(iattr)))
      enddo

c     Exit if not in rank=0

      if (.not.isProc(0)) return

c     Initialize FORTRAN interface
      call h5open_f(hdferr)

c     Create HDF5 file

      !Create a new file using default properties.
      call h5fcreate_f(hdf5_fname, H5F_ACC_TRUNC_F, file_id, hdferr)
      
      !Terminate access to the file.
      call h5fclose_f(file_id, hdferr)

      !Open .h5 file using default properties
      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,hdferr
     $              ,H5P_DEFAULT_F)

c     Create VISIT expressions (to identify vectors)

      call selectExpressions

      slth = len(trim(vexp))

      if (slth > 0) then

        !Create dataspace
        call h5screate_simple_f(VEXP_RANK,vexp_dim,dataspace_id,hdferr)

        !Create datatype
        call h5tcopy_f(H5T_NATIVE_CHARACTER,datatype_id,hdferr)
        call h5tset_size_f(datatype_id,slth,hdferr)
        call h5tset_strpad_f(datatype_id,H5T_STR_NULLTERM_F,hdferr)

        !Create dataset
        call h5dcreate_f(file_id,'/visit_expressions',datatype_id
     .                  ,dataspace_id,dataset_id,hdferr,H5P_DEFAULT_F)

        !Write to dataset
        call h5dwrite_f(dataset_id,datatype_id,trim(vexp),vexp_dim
     .                 ,hdferr)

        !Close dataset
        call h5dclose_f(dataset_id,hdferr)

        !Close datatype
        call h5tclose_f(datatype_id,hdferr)

        !Close dataspace
        call h5sclose_f(dataspace_id,hdferr)

      endif

c     Create groups for coordinate data

      !Cell center positions
      call h5gcreate_f(file_id,cell_coords_gname,group_id,hdferr,izero)

      call h5gclose_f(group_id,hdferr)
      !Node positions
      call h5gcreate_f(file_id,node_coords_gname,group_id,hdferr,izero)
      call h5gclose_f(group_id,hdferr)

c     Close fortran interface

      call h5fclose_f(file_id, hdferr)   !Close .h5 file

      call h5close_f(hdferr)

c     End program

      contains

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

      character(20) :: str1,str2,str3
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

      end subroutine initializeHDF5

c     createGroups
c     #################################################################
      subroutine createGroups(groupname)

      implicit none

c     Call variables
      
cc      integer            :: step
      character*(20)     :: groupname

c     Local variables

      integer            :: ieq

c     Begin program

      if (.not.isProc(0)) return

c     Initialize FORTRAN interface

      call h5open_f(hdferr)

c     Open .h5 file using default properties.

      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,hdferr,
     $               H5P_DEFAULT_F)

ccc     Create a group for coordinates data
cc
cc      if (step == 0) then
cc        !Cell center positions
cc        call h5gcreate_f(file_id,cell_coords_gname,group_id,hdferr,0)
cc        call h5gclose_f(group_id,hdferr)
cc        !Node positions
cc        call h5gcreate_f(file_id,node_coords_gname,group_id,hdferr,0)
cc        call h5gclose_f(group_id,hdferr)
cc      endif

c     Create groups for computation data

      call h5gcreate_f(file_id,groupname,group_id,hdferr, izero)

      do igrp=1,ngroups
        call h5gcreate_f(group_id,graph(igrp)%descr,subgroup_id,
     $                  hdferr,izero)
        call h5gclose_f(subgroup_id, hdferr)
      enddo

      call h5gclose_f(group_id, hdferr)

c     Terminate access to the file.
      call h5fclose_f(file_id, hdferr)

      call h5close_f(hdferr)

c     End program

      end subroutine createGroups

c     writeHDF5file
c     ##################################################################
      subroutine writeHDF5file
c     ------------------------------------------------------------------
c     Write the output date using HDF5 HyperSlab Column
c     ------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer        :: i,j,k,ieq,icoords

      character*(20) :: timestep

      integer(hsize_t), dimension(ARRAY_RANK) :: count,node_count
      integer(hsize_t), dimension(ARRAY_RANK) :: offset

      real(8),allocatable,dimension(:,:,:,:) :: cell,node

c     Begin program

c     Create the output file and data groups

      if (itime == -1) toffset = 1

      timestep = trim('Timestep_'//int2char(itime+toffset))

      call createGroups(timestep)

c     Initialize FORTRAN interface
      call h5open_f(hdferr)

#if defined(petsc)
c PARALLEL
c     Set up file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdferr)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,
     $                        hdferr)

c     Open an existing file.
      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,hdferr,
     $               access_prp = plist_id)
      call h5pclose_f(plist_id,hdferr)
#else
c SERIAL
c     Open an existing file.
      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,hdferr,
     $               H5P_DEFAULT_F)
#endif

c     Create the dataspace of actual data and its attributes

      !Cell-centered data
      call h5screate_simple_f(ARRAY_RANK,array_dim,dataspace_id,hdferr)

      !Node-centered data
      call h5screate_simple_f(ARRAY_RANK,node_array_dim
     .                       ,node_dataspace_id,hdferr)

      !Attribute data
      call h5screate_simple_f(ATTR_RANK,attr_dim,dataspace_attr_id
     $                       ,hdferr)

#if defined(petsc)
c PARALLEL
c     Each process defines dataset in memory and writes it to the hyperslab
c     in the file.

      offset(1) = ilomg
      offset(2) = jlomg
      offset(3) = klomg

      !Cell-centered data
      count(1) = nxl
      count(2) = nyl
      count(3) = nzl

      call h5screate_simple_f(ARRAY_RANK,count,p_dataspace_id,hdferr)

      !Node-centered data
      node_count(1) = nxl+1
      node_count(2) = nyl+1
      node_count(3) = nzl+1

      call h5screate_simple_f(ARRAY_RANK,node_count,p_node_dataspace_id
     .                       ,hdferr)

c     Create property list for collective dataset write.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdferr)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdferr)
#endif

c     Write the coordinate data

      if ((itime+toffset) == 0) then

        allocate (cell(ilom:ihip,jlom:jhip,klom:khip,3)
     .           ,node(ilom:ihip,jlom:jhip,klom:khip,3))

        !Cell center positions
        call h5gopen_f(file_id,cell_coords_gname,group_id,hdferr) 

        do k=klom,khip
          do j=jlom,jhip
            do i=ilom,ihip
              call getCartesianCoordinates(i,j,k,iggx,iggy,iggz
     .          ,iig,jjg,kkg,cell(i,j,k,1),cell(i,j,k,2),cell(i,j,k,3))
            enddo
          enddo
        enddo

        do icoords=1,ARRAY_RANK

#if defined(petsc)
c PARALLEL
          call h5dcreate_f(group_id,coords_descr(icoords),
     $                     H5T_IEEE_F64BE,
     $                     dataspace_id,dataset_id,hdferr)

c         Select hyperslab in the file.
          call h5dget_space_f(dataset_id,hs_dataspace_id,hdferr)
          call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F,
     $                               offset,count,hdferr)    

c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                 cell(ilo:ihi,jlo:jhi,klo:khi,icoords),
     $                 array_dim_i,
     $                 hdferr,file_space_id = hs_dataspace_id,
     $                 mem_space_id = p_dataspace_id,
     $                 xfer_prp = plist_id)

          call h5sclose_f(hs_dataspace_id,hdferr)
#else
c SERIAL
          call h5dcreate_f(group_id,coords_descr(icoords),
     $                     H5T_IEEE_F64BE,
     $                     dataspace_id,dataset_id,hdferr)

c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                    cell(ilo:ihi,jlo:jhi,klo:khi,icoords), !must specify index ranges
     $                    array_dim,hdferr)
          
#endif
          call h5dclose_f(dataset_id,hdferr)
        enddo

        call h5gclose_f(group_id, hdferr)

        !Node positions
        call h5gopen_f(file_id,node_coords_gname,group_id,hdferr) 

        call interpolateNodalData

        do icoords=1,ARRAY_RANK
          call h5dcreate_f(group_id,coords_descr(icoords),
     $                     H5T_IEEE_F64BE,
     $                     node_dataspace_id,dataset_id,hdferr)

#if defined(petsc)
c PARALLEL
c         Select hyperslab in the file.
          call h5dget_space_f(dataset_id,hs_dataspace_id,hdferr)
          call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F,
     $                               offset,node_count,hdferr)    

c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                 node(ilo:ihip,jlo:jhip,klo:khip,icoords),
     $                 array_dim_i,
     $                 hdferr,
     .                 file_space_id = hs_dataspace_id,
     $                 mem_space_id = p_node_dataspace_id,
     $                 xfer_prp = plist_id)

          call h5sclose_f(hs_dataspace_id,hdferr)
#else
c SERIAL
c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                    node(ilo:ihip,jlo:jhip,klo:khip,icoords),
     $                    node_array_dim,hdferr)
#endif

          call h5dclose_f(dataset_id,hdferr)
        enddo

        call h5gclose_f(group_id, hdferr)

        deallocate(cell,node)

      endif

c     Write the actual data

      call h5gopen_f(file_id,timestep,group_id,hdferr)

      do igrp=1,ngroups

c       Open an existing group of the specified file.
        call h5gopen_f(group_id,graph(igrp)%descr,subgroup_id,hdferr)

        if (graph(igrp)%cartesian) then

c         Write "Attribute data"

          call h5tcopy_f(H5T_NATIVE_CHARACTER,attr_type_id,hdferr)
          call h5tset_size_f(attr_type_id,strlength,hdferr)

          call h5acreate_f(subgroup_id,attr_name,attr_type_id,
     $                     dataspace_attr_id,attribute_id,
     $                     hdferr)

          call h5awrite_f(attribute_id,attr_type_id,
     $                    attr_data,
     $                    attr_dim,hdferr)
          call h5aclose_f(attribute_id,hdferr)

          call h5tclose_f(attr_type_id,hdferr)

c         Write nodal data

          allocate (cell(ilom:ihip,jlom:jhip,klom:khip,1)
     .             ,node(ilom:ihip,jlom:jhip,klom:khip,1))

          do ieq=1,nqty(igrp)

c           Create the dataset
            call h5dcreate_f(subgroup_id
     $                      ,graph(igrp)%array_graph(ieq)%descr
     $                      ,H5T_IEEE_F64BE
     $                      ,node_dataspace_id,dataset_id,hdferr)

c           Interpolate nodal values

            cell(:,:,:,1) = graph(igrp)%array_graph(ieq)%array

cc            write (*,*) cell(:,:,:,1)
cc            pause
            call interpolateNodalData

cc            write (*,*) node(:,:,:,1)
cc            pause

#if defined(petsc)
c PARALLEL
c           Select hyperslab in the file.
            call h5dget_space_f(dataset_id,hs_dataspace_id,hdferr)
            call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F
     $                                ,offset,node_count,hdferr)    

c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,node(ilo:ihip,jlo:jhip,klo:khip,1)
     $                     ,array_dim_i
     $                     ,hdferr
     $                     ,file_space_id = hs_dataspace_id
     $                     ,mem_space_id = p_node_dataspace_id
     $                     ,xfer_prp = plist_id)

c           Release resources
            call h5sclose_f(hs_dataspace_id,hdferr)
#else
c SERIAL
c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,node(ilo:ihip,jlo:jhip,klo:khip,1)
     $                     ,node_array_dim,hdferr)
#endif

c         Release resources

            call h5dclose_f(dataset_id,hdferr)

          enddo

          deallocate(cell,node)

        else

c         Write cell-centered data

          do ieq=1,nqty(igrp)

c           Create the dataset
            call h5dcreate_f(subgroup_id
     $                      ,graph(igrp)%array_graph(ieq)%descr
     $                      ,H5T_IEEE_F64BE
     $                      ,dataspace_id,dataset_id,hdferr)

#if defined(petsc)
c PARALLEL
c           Select hyperslab in the file.
            call h5dget_space_f(dataset_id,hs_dataspace_id,hdferr)
            call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F
     $                                ,offset,count,hdferr)    

c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,graph(igrp)%array_graph(ieq)
     $                         %array(ilo:ihi,jlo:jhi,klo:khi)
     $                     ,array_dim_i
     $                     ,hdferr,file_space_id = hs_dataspace_id
     $                     ,mem_space_id = p_dataspace_id
     $                     ,xfer_prp = plist_id)

c           Release resources
            call h5sclose_f(hs_dataspace_id,hdferr)
#else
c SERIAL
c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,graph(igrp)%array_graph(ieq)
     $                                 %array(ilo:ihi,jlo:jhi,klo:khi) !must specify index ranges
     $                     ,array_dim,hdferr)

#endif

c         Release resources

            call h5dclose_f(dataset_id,hdferr)

          enddo

        endif

        call h5gclose_f(subgroup_id, hdferr)

      enddo

c     Release resources

      call h5gclose_f(group_id, hdferr)

#if defined(petsc)
c PARALLEL
      call h5sclose_f(p_dataspace_id,hdferr)
      call h5sclose_f(p_node_dataspace_id,hdferr)
      call h5pclose_f(plist_id,hdferr)
#endif

      call h5sclose_f(dataspace_id,hdferr)
      call h5sclose_f(node_dataspace_id,hdferr)
      call h5fclose_f(file_id,hdferr)
      call h5close_f (hdferr)

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
              if (i == ilo .and. ilog == 1 .and. bcSP()) then
                node(i,j,k,:) = 0.25*(cell(i-1,j  ,k  ,:)
     .                               +cell(i-1,j-1,k  ,:)
     .                               +cell(i-1,j  ,k-1,:)
     .                               +cell(i-1,j-1,k-1,:))
              else
                node(i,j,k,:) = 0.125
     .                       *(cell(i  ,j  ,k  ,:)+cell(i-1,j  ,k  ,:)
     .                        +cell(i  ,j-1,k  ,:)+cell(i-1,j-1,k  ,:)
     .                        +cell(i  ,j  ,k-1,:)+cell(i-1,j  ,k-1,:)
     .                        +cell(i  ,j-1,k-1,:)+cell(i-1,j-1,k-1,:))
              endif
            enddo
          enddo
        enddo

      end subroutine interpolateNodalData

      end subroutine writeHDF5file

#endif

      end module HDF5_graphics

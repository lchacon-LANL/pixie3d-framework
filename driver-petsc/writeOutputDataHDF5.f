c     module HDF5_graphics
c     ##############################################################
      module HDF5_graphics

        use HDF5

        include 'mpif.h'
      
        integer(hid_t) :: group_id       
        integer(hid_t) :: attribute_id   
        integer(hid_t) :: dataspace_att_id
        integer(hid_t) :: file_id        
        integer(hid_t) :: dataset_id     
        integer(hid_t) :: dataspace_id   
        integer(hid_t) :: plist_id       
        integer(hid_t) :: p_dataspace_id 
        integer(hid_t) :: hs_dataspace_id

        integer(4)     :: error

        integer(hsize_t), dimension(3) :: array_dim
        integer(hsize_t), dimension(1) :: attr_dim = (/1/)
        integer(4), dimension(7) :: array_dim_i 
        
        integer(4), parameter  :: ARRAY_RANK = 3
        integer(4), parameter  :: ATTR_RANK = 1

      contains 

c     createGroup
c     #################################################################
      subroutine createGroup(groupname, filename)

      implicit none

c     call variables
      
      integer(4)         :: step
      character*(*)      :: filename,groupname

c     local variables

      integer(4)     :: ieq,my_rank,mpierr

      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)

      if (my_rank /= 0) return

c     Initialize FORTRAN interface

      call h5open_f(error)

c     Create a new file using default properties.

      call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error,
     $               H5P_DEFAULT_F)

c     Create a new group.

      call h5gcreate_f(file_id,groupname,group_id,error, 0)
      call h5gclose_f(group_id, error)

c     Terminate access to the file.
      call h5fclose_f(file_id, error)

      call h5close_f(error)

c     End program

      end subroutine createGroup

c     int2char
c     ##################################################################
      function int2char(n) result (chr)

      integer(4)   :: n
      character(10):: chr

      integer(4)   :: i,j,exp,k
      character(3) :: c

      if (n > 0) then
         exp = int(log(float(n))/log(1d1))
      else
         exp = 0
      endif

      chr=''
      k = n
      do i=exp,0,-1
         j = k/10**i
         c = achar(48+j)
         chr = trim(chr)//trim(c)
         if (i > 0 ) k = mod(k,j*10**i)
      enddo

      end function int2char

      end module HDF5_graphics

c     initializeHDF5file
c     ##############################################################
      subroutine initializeHDF5file(filename)

      use HDF5_graphics

      implicit none

c     call variables
      
      character(len=256) :: filename

c     local variables

      integer(4)     :: ieq

c     Initialize FORTRAN interface
      call h5open_f(error)

c     Create a new file using default properties.
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

c     Terminate access to the file.
      call h5fclose_f(file_id, error)

      call h5close_f(error)

c     End program

      end subroutine initializeHDF5file

c     writeHDF5file
c     ##################################################################
      subroutine writeHDF5file(imin,imax,jmin,jmax,kmin,kmax
     $                        ,xm,ym,zm,step,dt,filename,xvec)

c     ------------------------------------------------------------------
c     Write the output date using HDF5 HyperSlab Column
c     ------------------------------------------------------------------

      use HDF5_graphics

      use variables

      implicit none

c     Call variables

      integer(4)     :: imin,imax,jmin,jmax,kmin,kmax,xm,ym,zm,step
      integer(4)     :: nx,ny,nz
      real(8)        :: dt
      character(len=256) :: filename

      type(petsc_var) :: xvec (imin:imax,jmin:jmax,kmin:kmax)

c     Local variables

      integer(4)     :: i,j,k,ieq,my_rank,mpierr

      character*(20) :: groupname
      character*(20) :: attr_name = "Actual time"

      real(8)        :: attr_data

      integer(hsize_t), dimension(ARRAY_RANK) :: count(3)
      integer(hsize_t), dimension(ARRAY_RANK) :: offset(3)

c     Begin program

      array_dim = (/nxd,nyd,nzd/)
      array_dim_i = (/nxd,nyd,nzd,0,0,0,0/)
      attr_name = "Actual Time"

c     Create the output file and data groups

      groupname = int2char(step)
      groupname = trim('Timestep '//int2char(step))

      call createGroup(groupname, filename)

c     Initialize FORTRAN interface
      call h5open_f(error)

c     Set up file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,
     $                        error)

c     Open an existing file.
      call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error,
     $               access_prp = plist_id)
      call h5pclose_f(plist_id,error)

c     Write the actual data in each group

c     Create the dataspace of actual data and its attributes
c     ARRAY_RANK = rank of array,array_dim = dimension of array
      call h5screate_simple_f(ARRAY_RANK,array_dim,dataspace_id,error)

      call h5screate_simple_f(ATTR_RANK,attr_dim,dataspace_att_id,error) 

      attr_data = step * dt

c     Open an existing group of the specified file.
      call h5gopen_f(file_id,groupname,group_id,error)

c     Each process defines dataset in memory and writes it to the hyperslab
c     in the file.
      
      count(1) = xm
      count(2) = ym
      count(3) = zm
      offset(1) = imin-1
      offset(2) = jmin-1
      offset(3) = kmin-1

      call h5screate_simple_f(ARRAY_RANK,count,p_dataspace_id,error)

c     Create property list for collective dataset write.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error) !omit two optional arguments

      do ieq=1,neqd

c     Create the dataset
         call h5dcreate_f(group_id,u_0%array_var(ieq)%descr,
     $                    H5T_IEEE_F64BE,
     $                    dataspace_id,dataset_id,error)
    
c     Select hyperslab in the file.
         call h5dget_space_f(dataset_id,hs_dataspace_id,error)
         call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F,
     $                              offset,count,error)

c     Write data to the dataset
         call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                   xvec(imin:imax,jmin:jmax,kmin:kmax)%var(ieq),
     $                   array_dim_i,
     $                   error,file_space_id = hs_dataspace_id,
     $                   mem_space_id = p_dataspace_id,
     $                   xfer_prp = plist_id)

c     Release resources
         call h5dclose_f(dataset_id,error)
         call h5sclose_f(hs_dataspace_id,error)

      enddo
  
c     Release resources
      call h5sclose_f(p_dataspace_id,error)
      call h5pclose_f(plist_id,error)
      call h5sclose_f(dataspace_id,error)
      call h5sclose_f(dataspace_att_id,error)
      call h5gclose_f(group_id, error)
      call h5fclose_f(file_id,error)

      call h5close_f(error)

c     End program

      end subroutine writeHDF5file

c     readHDF5file
c     ##################################################################
      subroutine readHDF5file(imin,imax,jmin,jmax,kmin,kmax
     $                        ,xm,ym,zm,step,dt,filename,xvec)

c     ------------------------------------------------------------------
c     Write the output date using HDF5 HyperSlab Column
c     ------------------------------------------------------------------

      use HDF5_graphics

      use variables

      implicit none

c     Call variables

      integer(4)     :: imin,imax,jmin,jmax,kmin,kmax,xm,ym,zm,step
      integer(4)     :: nx,ny,nz
      real(8)        :: dt
      character(len=256) :: filename

      type(petsc_var) :: xvec (imin:imax,jmin:jmax,kmin:kmax)

c     Local variables

      integer(4)     :: i,j,k,ieq,my_rank,mpierr

      character*(20) :: groupname
      character*(20) :: attr_name = "Actual time"

      real(8)        :: attr_data

      integer(hsize_t), dimension(ARRAY_RANK) :: count(3)
      integer(hsize_t), dimension(ARRAY_RANK) :: offset(3)

c     Begin program

      array_dim_i = (/nxd,nyd,nzd,0,0,0,0/)

c     Initialize FORTRAN interface
      call h5open_f(error)

c     Set up file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,
     $                        error)

c     Open an existing file.
      call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error,
     $               access_prp = plist_id)
      call h5pclose_f(plist_id,error)

c     Read the actual data in each group

      attr_data = step * dt

      groupname = int2char(step)
      groupname = 'Timestep '//groupname

c     Open an existing group of the specified file.
      call h5gopen_f(file_id,groupname,group_id,error)

c     Each process defines dataset in memory and writes it to the hyperslab
c     in the file.
      
      count(1) = xm
      count(2) = ym
      count(3) = zm
      offset(1) = imin-1
      offset(2) = jmin-1
      offset(3) = kmin-1

      call h5screate_simple_f(ARRAY_RANK,count,p_dataspace_id,error)

c     Create property list for collective dataset write.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error) !omit two optional arguments

      do ieq=1,neqd

c     Open the existing dataset
         call h5dopen_f(group_id,u_0%array_var(ieq)%descr,
     $                  dataset_id,error)
    
c     Select hyperslab in the file.
         call h5dget_space_f(dataset_id,hs_dataspace_id,error)
         call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F,
     $                              offset,count,error)

c     Write data to the dataset
         call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                  xvec(imin:imax,jmin:jmax,kmin:kmax)%var(ieq),
     $                  array_dim_i,
     $                  error,file_space_id = hs_dataspace_id,
     $                  mem_space_id = p_dataspace_id,
     $                  xfer_prp = plist_id)

c     Release resources
         call h5dclose_f(dataset_id,error)
         call h5sclose_f(hs_dataspace_id,error)

      enddo
  
c     Release resources
      call h5sclose_f(p_dataspace_id,error)
      call h5pclose_f(plist_id,error)
      call h5gclose_f(group_id, error)
      call h5fclose_f(file_id,error)

      call h5close_f(error)

c     End program

      end subroutine readHDF5file

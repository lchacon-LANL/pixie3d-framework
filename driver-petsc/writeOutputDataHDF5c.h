/* Write the output data using HDF5 library */

#undef __FUNCT__
#define __FUNCT__ "initializeOutputFile"
int initializeOutputFile(char* datasetname[], char* filename)
{
  int      ivar;

  hid_t    file_id;     /* file identifier  */
  herr_t   status;


  /* Create a new file using default properties. */
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Terminate access to the file. */
  status = H5Fclose(file_id); assert (status >= 0);

  return 0;
}
#undef __FUNCT__
#define __FUNCT__ "createGroupDataset"
int createGroup(int step, char* filename)
{
  hid_t      group_id;             /* group identifier                  */
  hid_t      file_id;              /* file identifier                   */
  herr_t     status;

  char       groupname[256];

  file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);assert(file_id >= 0);

  sprintf(groupname, "TimeStep%d", step);
  /* Create the group */
  group_id = H5Gcreate(file_id, groupname, 0); assert(group_id >= 0);

  /* Release resources */
  status = H5Gclose(group_id); assert (status >= 0);
  status = H5Fclose(file_id); assert (status >= 0);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeOutputDataRoutineHSR"
int writeOutputDataRoutineHSR(int nx, int ny, int nz, int xs, int ys, int zs,
			      int xm, int ym, int zm,
			      int step, double dt,
			      char* datasetname[], char* filename,
			      Field*** xvec)
{
  hid_t      group_id;             /* group identifier                  */
  hid_t      attribute_id;         /* attribute identifier              */
  hid_t      dataspace_att_id;     /* attribute dataspace id            */
  hid_t      file_id;              /* file identifier                   */
  hid_t      dataset_id;           /* dataset identifier of array       */
  hid_t      dataspace_id;         /* dataspace id of array             */
  hid_t      plist_id;             /* property list identifier          */
  hid_t      p_dataspace_id;       /* dataspace id for a parallel       */
  hid_t      hs_dataspace_id;      /* dataspace id in a hyperslab       */

  herr_t     status;

  hsize_t    array_dim[] = {nz,ny,nx};       /* array dimensions         */
  hsize_t    attr_dim[] = {1};               /* attribute data dimension */

  int        ARRAY_RANK = 3;       /* rank of array                     */
  int        ATT_RANK = 1;         /* rank of attribute                 */
  
  int        i,j,k,ivar,my_rank;

  char       groupname[256];       /* dataset name of each timestep     */
  char       attr_name[256];       /* dataset name of attribute         */

  double     attr_data;

  hsize_t    count[ARRAY_RANK];              /* hyperslab selection parameters */
  hsize_t    offset[ARRAY_RANK];

  typedef struct var_def {
    double   array[zm][ym][xm];
    /* int      bcs[bcsdim];
       char     descr[descrdim];*/
  } var_def;

  var_def*    varray;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  sprintf(attr_name, "Actual Time");

  varray = (var_def*)malloc(NVAR * sizeof(var_def));

  for (ivar = 0; ivar < NVAR; ivar++) {
    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
	for (i = xs; i < xs+xm; i++) {
	  varray[ivar].array[k-zs][j-ys][i-xs] = xvec[k][j][i].var[ivar];
	}
      }
    }
  }

  /*************************************************************************/
  /* Create the output file and data groups                                */
  /*************************************************************************/

  /* Set up file access property list with parallel I/O access */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);assert(plist_id >= 0);
  status = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(status >= 0);

  /* Open an existing file. */
  file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);assert(file_id >= 0);
  status = H5Pclose(plist_id); assert(plist_id >= 0);

  /*************************************************************************/
  /* Write the actual data in each group                                   */
  /*************************************************************************/

  /* Create the dataspace of actual data and its attributes */
  /* ARRAY_RANK = rank of array, array_dim = dimension of array */
  dataspace_id = H5Screate_simple(ARRAY_RANK, array_dim, NULL); 
  assert(dataspace_id >= 0);

  sprintf(groupname, "TimeStep%d", step);

  /* Open the group */ 
  group_id = H5Gopen(file_id, groupname);assert(group_id >= 0);

  for (ivar = 0; ivar < NVAR; ivar++) {
    /* Create the dataset */
    dataset_id = H5Dcreate(group_id, datasetname[ivar], H5T_IEEE_F64BE,
			   dataspace_id, H5P_DEFAULT);assert(dataset_id >= 0);

    /* Each process defines dataset in memory and writes it to the hyperslab
       in the file.*/
      
    count[0] = zm;
    count[1] = ym;
    count[2] = xm;
    offset[0] = zs;
    offset[1] = ys;
    offset[2] = xs;
    p_dataspace_id = H5Screate_simple(ARRAY_RANK, count, NULL);
    assert(p_dataspace_id >= 0);
    
    /* Select hyperslab in the file. */
    hs_dataspace_id = H5Dget_space(dataset_id);
    status = H5Sselect_hyperslab(hs_dataspace_id, H5S_SELECT_SET, offset, 
				 NULL, count, NULL); assert(status >= 0);

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);assert(dataset_id >= 0);
    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    assert(status >= 0);
    
    /* To write dataset independently use
     * H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); */
    
    /* Write data to the dataset */
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, p_dataspace_id, 
		      hs_dataspace_id, plist_id,
		      &varray[ivar].array[0][0][0]);
    assert(status >= 0 );

    /* Release resources */
    status = H5Dclose(dataset_id); assert(status >= 0); 
    status = H5Sclose(p_dataspace_id); assert(status >= 0);
    status = H5Sclose(hs_dataspace_id); assert(status >= 0);
    status = H5Pclose(plist_id); assert(status >= 0);
  }
  
  /* Release resources */
  free(varray);
  status = H5Sclose(dataspace_id); assert (status >= 0);
  status = H5Gclose(group_id); assert (status >= 0);
  status = H5Fclose(file_id); assert (status >= 0);

  /* Return value */
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readOutputDataRoutineHSR"
int readOutputDataRoutineHSR(int nx, int ny, int nz, int xs, int ys, int zs,
			     int xm, int ym, int zm,
			     int step, double dt,
			     char* datasetname[], char* filename,
			     Field*** xvec)
{
  hid_t      group_id;             /* group identifier                  */
  hid_t      attribute_id;         /* attribute identifier              */
  hid_t      dataspace_att_id;     /* attribute dataspace id            */
  hid_t      file_id;              /* file identifier                   */
  hid_t      dataset_id;           /* dataset identifier of array       */
  hid_t      dataspace_id;         /* dataspace id of array             */
  hid_t      plist_id;             /* property list identifier          */
  hid_t      p_dataspace_id;       /* dataspace id for a parallel       */
  hid_t      hs_dataspace_id;      /* dataspace id in a hyperslab       */

  herr_t     status;

  hsize_t    array_dim[] = {nz,ny,nx};       /* array dimensions         */
  hsize_t    attr_dim[] = {1};               /* attribute data dimension */

  int        ARRAY_RANK = 3;       /* rank of array                     */
  int        ATT_RANK = 1;         /* rank of attribute                 */
  
  int        i,j,k,ivar,my_rank;

  char       groupname[256];       /* dataset name of each timestep     */
  char       attr_name[256];       /* dataset name of attribute         */

  double     attr_data;

  hsize_t    count[ARRAY_RANK];              /* hyperslab selection parameters */
  hsize_t    offset[ARRAY_RANK];

  double     dummy_array[zm][ym][xm];


  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  sprintf(attr_name, "Actual Time");


  /*************************************************************************/
  /* Create the output file and data groups                                */
  /*************************************************************************/

  /* Set up file access property list with parallel I/O access */
  plist_id = H5Pcreate(H5P_FILE_ACCESS);assert(plist_id >= 0);
  status = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(status >= 0);

  /* Open an existing file. */
  file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);assert(file_id >= 0);
  status = H5Pclose(plist_id); assert(plist_id >= 0);

  /*************************************************************************/
  /* Write the actual data in each group                                   */
  /*************************************************************************/

  /* Create the dataspace of actual data and its attributes */
  /* ARRAY_RANK = rank of array, array_dim = dimension of array */
  /*dataspace_id = H5Screate_simple(ARRAY_RANK, array_dim, NULL); 
    assert(dataspace_id >= 0);*/

  sprintf(groupname, "TimeStep%d", step);

  /* Open the existing group */ 
  group_id = H5Gopen(file_id, groupname);assert(group_id >= 0);

  for (ivar = 0; ivar < NVAR; ivar++) {
    /* Open the existing dataset */
    dataset_id = H5Dopen(group_id, datasetname[ivar]);assert(dataset_id >= 0);

    /* Each process defines dataset in memory and writes it to the hyperslab
       in the file.*/
      
    count[0] = zm;
    count[1] = ym;
    count[2] = xm;
    offset[0] = zs;
    offset[1] = ys;
    offset[2] = xs;
    p_dataspace_id = H5Screate_simple(ARRAY_RANK, count, NULL);
    assert(p_dataspace_id >= 0);
    
    /* Select hyperslab in the file. */
    hs_dataspace_id = H5Dget_space(dataset_id);
    status = H5Sselect_hyperslab(hs_dataspace_id, H5S_SELECT_SET, offset, 
				 NULL, count, NULL); assert(status >= 0);

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);assert(dataset_id >= 0);
    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    assert(status >= 0);
    
    /* To write dataset independently use
     * H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); */
    
    /* Write data to the dataset */
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, p_dataspace_id, 
		     hs_dataspace_id, plist_id,
		     dummy_array);
    assert(status >= 0 );

    for (k = zs; k < zs+zm; k++) {
      for (j = ys; j < ys+ym; j++) {
	for (i = xs; i < xs+xm; i++) {
	  xvec[k][j][i].var[ivar] = dummy_array[k-zs][j-ys][i-xs];
	}
      }
    }

    /* Release resources */
    status = H5Dclose(dataset_id); assert(status >= 0); 
    status = H5Sclose(p_dataspace_id); assert(status >= 0);
    status = H5Sclose(hs_dataspace_id); assert(status >= 0);
    status = H5Pclose(plist_id); assert(status >= 0);
  }
  
  /* Release resources */
  /*  status = H5Sclose(dataspace_id); assert (status >= 0);*/
  status = H5Gclose(group_id); assert (status >= 0);
  status = H5Fclose(file_id); assert (status >= 0);

  /* Return value */
  return 0;
}

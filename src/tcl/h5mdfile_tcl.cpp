#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include "parser.hpp"
#include "communication.hpp"
#include <stdio.h>
#include <iostream>

#ifdef H5MD

#include "hdf5.h"

typedef char h5string[1000];
class H5mdfile
{
	public:
	hid_t    		file_id; // File id returned from H5Fopen
	hid_t    		dataset_id; // Dataset id returned from H5Dopen2 or H5Dcreate2
	hid_t    		dataspace_id; // Dataspace id returned from H5Dget_space
	hid_t    		dataspace_simple_id; // Dataspace id returned from H5_Screate_simple
	hid_t    		dataset_type_id; // Dataset data type id
	hid_t    		group_id; // Group id
	hid_t  			prop_id; // Property id returned from H5Pcreate and H5Dget_create_plist

	hsize_t*  	dims; // Dimensions of complete respectively extended dataset (extended dataset is part of complete dataset)
	hsize_t*  	dimstotal; // Dimensions of complete dataset
	hsize_t*  	maxdims; // Maximal dimensions of complete respectively extended dataset (always set to UNLIMITED)
	hsize_t*  	offset; // Offset for extended dataset in accordance to complete dataset
	hsize_t*  	chunk_dims; // Dimensions of  chunked dataset (complete dataset split in smaller pieces)

	int 				dataset_rank; // Rank of dataset (number of the single dimensions)
	int 				chunk_rank; // Rank of chunked dataset (number of the single dimensions)
	int					datatype_size; // Size of data type
	herr_t   		status; // Return error value of several functions

	int 				dset_data_size; // Length of complete respectively extended dataset array (length of all dimension serialized to one array)
	int					string_maxlength; // Maximal string length

	void*				dset_data; // Array of complete respectively extended dataset (all dimension serialized to one array)
	int*				dset_data_int; // Dataset array (type: integer)
	float*			dset_data_float; // Dataset array (type: float)
	double*			dset_data_double; // Dataset array (type: double)
	h5string*		dset_data_string; // Dataset array (type: string)
	int 				dset_data_singlevalue_int; // Single value of complete respectively extended dataset array (type: int)
	float 			dset_data_singlevalue_float; // Single value of complete respectively extended dataset array (type: float)
	double 			dset_data_singlevalue_double; // Single value of complete respectively extended dataset array (type: double)
	char*				dset_data_singlevalue_string; // Single value of complete respectively extended dataset array (type: string)

	// Constructor
	H5mdfile();

	// Close h5-dataset
	int H5_Dclose(int argc, char **argv, Tcl_Interp *interp);

	// Create h5-dataset and return dataset_id
	int H5_Dcreate2(int argc, char **argv, Tcl_Interp *interp);

	// Extend h5-dataset
	int H5_Dextend(int argc, char **argv, Tcl_Interp *interp);

	// Open h5-dataset and return dataset_id
	int H5_Dopen2(int argc, char **argv, Tcl_Interp *interp);

	// Read h5-dataset and write to dset_data
	int H5_Dread(int argc, char **argv, Tcl_Interp *interp);

	// Write dataset array to h5-file
	int H5_Dwrite(int argc, char **argv, Tcl_Interp *interp);

	// Close h5-file
	int H5_Fclose(int argc, char **argv, Tcl_Interp *interp);

	// Create h5-file and return file_id
	int H5_Fcreate(int argc, char **argv, Tcl_Interp *interp);

	// Open h-5 file and return file_id
	int H5_Fopen(int argc, char **argv, Tcl_Interp *interp);

	// Close h5-group
	int H5_Gclose(int argc, char **argv, Tcl_Interp *interp);

	// Create h5-group and return group_id
	int H5_Gcreate2(int argc, char **argv, Tcl_Interp *interp);

	// Open h5-group and return group_id
	int H5_Gopen2(int argc, char **argv, Tcl_Interp *interp);

	// Close property
	int H5_Pclose(int argc, char **argv, Tcl_Interp *interp);

	// Set properties for chunk
	int H5_Pset_chunk(int argc, char **argv, Tcl_Interp *interp);

	// Read single value from dataset and return to Tcl
	int H5_read_value(int argc, char **argv,Tcl_Interp *interp);

	// Close dataspace
	int H5_Sclose(int argc, char **argv, Tcl_Interp *interp);

	// Create simple dataspace
	int H5_Screate_simple(int argc, char **argv, Tcl_Interp *interp);

	// Selects a hyperslab region to add to the current selected region
	int H5_Sselect_hyperslab(int argc, char **argv, Tcl_Interp *interp);

	// Write single value from Tcl to dataset
	int H5_write_value(int argc, char **argv, Tcl_Interp *interp);

	// Free dataset memory
	int H5_free_memory(int argc, char **argv, Tcl_Interp *interp);

	// Get dataset dimensions
	int get_dataset_dims(int argc, char **argv, Tcl_Interp *interp);

	// Flush data to disc
	int H5_Fflush(int argc, char **argv, Tcl_Interp *interp);
};

H5mdfile h5mdfile;

int tclcommand_h5mdfile(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
	
	/* Parse Tcl commands */
	if (!strncmp(argv[1], "H5Dclose", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dclose\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Dclose(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Dcreate2", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dcreate2 \"path_to_dset\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Dcreate2(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Dextend", strlen(argv[1])))
	{
		if(strncmp(argv[2], "dims", strlen(argv[2])))
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dextend dims {<dims>}\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Dextend(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Dopen2", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dopen2 \"path_to_dset\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Dopen2(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Dread", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dread\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Dread(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Dwrite", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dwrite\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Dwrite(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Fclose", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Fclose\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Fclose(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Fcreate", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Fcreate \"path_to_file\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Fcreate(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Fopen", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Fopen \"path_to_file\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Fopen(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Gclose", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Gclose\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Gclose(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Gcreate2", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Gcreate2 \"path_to_group\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Gcreate2(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Gopen2", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Gopen2 \"path_to_group\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Gopen2(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Pclose", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Pclose\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Pclose(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Pset_chunk", strlen(argv[1])))
	{
		if(strncmp(argv[2], "dims", strlen(argv[2])))
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Pset_chunk dims {<dims>}\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Pset_chunk(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5_read_value", strlen(argv[1])))
	{
		if(strncmp(argv[2], "index", strlen(argv[2])))
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5_read_value index {<index>}\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_read_value(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Sclose", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Sclose\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Sclose(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Screate_simple", strlen(argv[1])))
	{
		if(strncmp(argv[2], "type", strlen(argv[2])) or strncmp(argv[4], "dims", strlen(argv[4])))
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Screate_simple type <type> dims {<dims>}\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Screate_simple(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Sselect_hyperslab", strlen(argv[1])))
	{
		if(strncmp(argv[2], "offset", strlen(argv[2])))
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Sselect_hyperslab offset {<offset>}\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Sselect_hyperslab(argc, argv, interp);
	}

	if (!strncmp(argv[1], "H5_write_value", strlen(argv[1])))
	{
		if(strncmp(argv[2], "value", strlen(argv[2])) or strncmp(argv[4], "index", strlen(argv[4])))
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5_write_value value <value> index {<index>}\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_write_value(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5_free_memory", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5_free_memory\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_free_memory(argc, argv, interp);
	}

	if (!strncmp(argv[1], "get_dataset_dims", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile get_dataset_dims\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.get_dataset_dims(argc, argv, interp);
	}
	
	if (!strncmp(argv[1], "H5_Fflush", strlen(argv[1])))
	{
		if(argc!=2)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5_Fflush\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Fflush(argc, argv, interp);
	}
	return TCL_ERROR;
}

H5mdfile::H5mdfile()
{
	/* Constructor */
	dims = new hsize_t[32];
	maxdims = new hsize_t[32];
  	dimstotal = new hsize_t[32];
	chunk_dims = new hsize_t[32];
	offset = new hsize_t[32];
	dset_data=NULL;
}
int H5mdfile::H5_Dclose(int argc, char **argv, Tcl_Interp *interp)
{
	/* End access to the dataset and release resources used by it */
	status = H5Dclose(dataset_id);
	return TCL_OK;
}
int H5mdfile::H5_Dcreate2(int argc, char **argv, Tcl_Interp *interp)
{
  /* Create the dataset */
	dataset_id = H5Dcreate2(file_id, argv[2], dataset_type_id, dataspace_simple_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
	dataspace_id = H5Dget_space (dataset_id);
	return TCL_OK;
}
int H5mdfile::H5_Dextend(int argc, char **argv, Tcl_Interp *interp)
{
	/* Extend dataset to higher dimensions */
	for(int i=0;i<dataset_rank;i++)
	{

	 if(atoi(argv[3+i])>(int)dims[i])
	 {
		 dims[i]=atoi(argv[3+i])-dimstotal[i];
	 }
	 dimstotal[i] = atoi(argv[3+i]);
	}

	status = H5Dextend(dataset_id, dimstotal);
	return TCL_OK;
}
int H5mdfile::H5_Dopen2(int argc, char **argv, Tcl_Interp *interp)
{
	/* Open an existing dataset */
	dataset_id = H5Dopen2(file_id, argv[2], H5P_DEFAULT);

	// Dataset properties
	dataspace_id = H5Dget_space(dataset_id);
	dataset_type_id = H5Dget_type(dataset_id);
	dataspace_simple_id = H5S_ALL;
	datatype_size = H5Tget_size(dataset_type_id);
	dataset_rank = H5Sget_simple_extent_dims(dataspace_id,dims,maxdims);
	dataset_rank = H5Sget_simple_extent_dims(dataspace_id,dimstotal,maxdims);
	prop_id = H5Dget_create_plist (dataset_id);
	if (H5D_CHUNKED == H5Pget_layout (prop_id))
	  chunk_rank = H5Pget_chunk (prop_id, dataset_rank, chunk_dims);
	// Dataset size
	dset_data_size=1;
	for(int i=0;i<dataset_rank;i++)
	{
	   dset_data_size*=dims[i];
	}
  return TCL_OK;
}
int H5mdfile::H5_Dread(int argc, char **argv, Tcl_Interp *interp)
{
	/* Read h5-dataset and write to dataset array values */
	// Allocate memeory
	if(dset_data!=NULL) free(dset_data);

	if(H5Tequal(dataset_type_id, H5T_NATIVE_FLOAT))
	{
	   dset_data=(float*) malloc(dset_data_size*sizeof(float));
	   memset(dset_data,0,dset_data_size*sizeof(float));
	}
	else if(H5Tequal(dataset_type_id, H5T_NATIVE_DOUBLE))
	{
	   dset_data=(double*) malloc(dset_data_size*sizeof(double));
	   memset(dset_data,0,dset_data_size*sizeof(double));
	}
	else if(H5Tequal(dataset_type_id, H5T_NATIVE_INT))
	{
	   dset_data=(int*) malloc(dset_data_size*sizeof(int));
	   memset(dset_data,0,dset_data_size*sizeof(int));
	}
	else if(H5Tequal(dataset_type_id, H5T_C_S1))
	{
	   dset_data = (h5string*) malloc(dset_data_size * sizeof(h5string));
	}
	else
	{
		Tcl_AppendResult(interp, "\nh5mdfile: No data type in H5_Dread given\n",(char *) NULL);
		return TCL_ERROR;
	}
	// Read h5-dataset
	status = H5Dread(dataset_id, dataset_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,dset_data);
	return TCL_OK;
}
int H5mdfile::H5_Dwrite(int argc, char **argv, Tcl_Interp *interp)
{
	/* Write dataset array to h5-file dataset */
	status = H5Dwrite(dataset_id, dataset_type_id, dataspace_simple_id, dataspace_id, H5P_DEFAULT, dset_data);
	return TCL_OK;
}
int H5mdfile::H5_Fclose(int argc, char **argv, Tcl_Interp *interp)
{
	/* Close h5-file */
	status = H5Fclose(file_id);
  return TCL_OK;
}
int H5mdfile::H5_Fcreate(int argc, char **argv, Tcl_Interp *interp)
{
	/* Create a new h5-file using default properties */
	file_id = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	return TCL_OK;
}
int H5mdfile::H5_Fopen(int argc, char **argv, Tcl_Interp *interp)
{
	/* Open an existing h5-file */
	file_id = H5Fopen(argv[2], H5F_ACC_RDWR, H5P_DEFAULT);
	return TCL_OK;
}
int H5mdfile::H5_Gclose(int argc, char **argv, Tcl_Interp *interp)
{
	/* Close h5-group */
	status = H5Gclose(group_id);
    return TCL_OK;
}
int H5mdfile::H5_Gcreate2(int argc, char **argv, Tcl_Interp *interp)
{
	/* Create a new h5-group using default properties */
	group_id = H5Gcreate2(file_id, argv[2], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	return TCL_OK;
}
int H5mdfile::H5_Gopen2(int argc, char **argv, Tcl_Interp *interp)
{
	/* Open an existing h5-group */
	group_id = H5Gopen2(file_id, argv[2], H5P_DEFAULT);
	return TCL_OK;
}
int H5mdfile::H5_Pclose(int argc, char **argv, Tcl_Interp *interp)
{
	/* Terminate access to the data space */
	status = H5Pclose(prop_id);
	return TCL_OK;
}
int H5mdfile::H5_Pset_chunk(int argc, char **argv, Tcl_Interp *interp)
{
	/* Set properties for chunk */
	for(int i=0;i<dataset_rank;i++)
	{
	   chunk_dims[i] = atoi(argv[3+i]);
	}
	prop_id = H5Pcreate (H5P_DATASET_CREATE);
	status = H5Pset_chunk (prop_id, dataset_rank, chunk_dims);
	return TCL_OK;
}
int H5mdfile::H5_read_value(int argc, char **argv,Tcl_Interp *interp)
{
	/* Read value from dataset array and print it to Tcl */
	// Get array index
	int index=0;
	if(dataset_rank>=1) index+=atoi(argv[2+dataset_rank-0]);
	if(dataset_rank>=2) index+=atoi(argv[2+dataset_rank-1])*dims[dataset_rank-1];else goto label_index;
	if(dataset_rank>=3) index+=atoi(argv[2+dataset_rank-2])*dims[dataset_rank-1]*dims[dataset_rank-2];else goto label_index;
	if(dataset_rank>=4) index+=atoi(argv[2+dataset_rank-3])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3];else goto label_index;
	if(dataset_rank>=5) index+=atoi(argv[2+dataset_rank-4])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4];else goto label_index;
	if(dataset_rank>=6) index+=atoi(argv[2+dataset_rank-5])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4]*dims[dataset_rank-5];else goto label_index;
	if(dataset_rank>=7) index+=atoi(argv[2+dataset_rank-6])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4]*dims[dataset_rank-5]*dims[dataset_rank-6];else goto label_index;
	if(dataset_rank>=8) index+=atoi(argv[2+dataset_rank-7])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4]*dims[dataset_rank-5]*dims[dataset_rank-6]*dims[dataset_rank-7];else goto label_index;
	label_index:
	// Read single value from dataset array and print to Tcl
	if(H5Tequal(dataset_type_id, H5T_NATIVE_FLOAT))
	{
		dset_data_float = static_cast<float*>(dset_data);
		dset_data_singlevalue_float = dset_data_float[index];

		char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
		Tcl_PrintDouble(interp, (double)dset_data_singlevalue_float, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type_id, H5T_NATIVE_DOUBLE))
	{
		dset_data_double = static_cast<double*>(dset_data);
		dset_data_singlevalue_double = dset_data_double[index];

		char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
		Tcl_PrintDouble(interp, dset_data_singlevalue_double, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type_id, H5T_NATIVE_INT))
	{
		dset_data_int = static_cast<int*>(dset_data);
		dset_data_singlevalue_int = dset_data_int[index];

		char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
		Tcl_PrintDouble(interp, (double)dset_data_singlevalue_int, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type_id, H5T_C_S1))
	{
		dset_data_string = static_cast<h5string*>(dset_data);
		dset_data_singlevalue_string = dset_data_string[index];
		dataspace_simple_id = H5S_ALL;
		dataspace_id = H5S_ALL;
		Tcl_AppendResult(interp, dset_data_singlevalue_string, (char *)NULL);
		return TCL_OK;
	}
	return TCL_ERROR;
}
int H5mdfile::H5_Sclose(int argc, char **argv, Tcl_Interp *interp)
{
	/* Terminate access to the data space */
	status = H5Sclose(dataspace_simple_id);
  return TCL_OK;
}
int H5mdfile::H5_Screate_simple(int argc, char **argv, Tcl_Interp *interp)
{
  /*  Create simple dataspace for dataset */

	// Get rank and dimension
	dataset_rank=0;
	dset_data_size=1;
	while(!(argv[5+dataset_rank]==NULL))
	{
	   dset_data_size*=atoi(argv[5+dataset_rank]);
	   dataset_rank++;
	}

	for(int i=0;i<dataset_rank;i++)
	{
	   dims[i] = atoi(argv[5+i]);
	   maxdims[i] = H5S_UNLIMITED;
	   if(dims[i]>dimstotal[i]) dimstotal[i] = dims[i];
	}
	// Cretae simple dataspace
	dataspace_simple_id = H5Screate_simple(dataset_rank, dims, maxdims);

	// Allocate memory for dataset array
	if(dset_data!=NULL) free(dset_data);

	if(!strncmp(argv[3], "float", strlen(argv[3])))
	{
	   dataset_type_id = H5T_NATIVE_FLOAT;
	   dset_data=(float*) malloc(dset_data_size*sizeof(float));
	   memset(dset_data,0,dset_data_size*sizeof(float));
	}
	else if(!strncmp(argv[3], "double", strlen(argv[3])))
	{
	   dataset_type_id = H5T_NATIVE_DOUBLE;
	   dset_data=(double*) malloc(dset_data_size*sizeof(double));
	   memset(dset_data,0,dset_data_size*sizeof(double));
	}
	else if(!strncmp(argv[3], "int", strlen(argv[3])))
	{
	   dataset_type_id = H5T_NATIVE_INT;
	   dset_data=(int*) malloc(dset_data_size*sizeof(int));
	   memset(dset_data,0,dset_data_size*sizeof(int));
	}
	else if(!strncmp(argv[3], "str", strlen(argv[3])))
	{
	   dataset_type_id = H5Tcopy(H5T_C_S1);
	   dset_data = (h5string*) malloc(dset_data_size * sizeof(h5string));
	}
	else
	{
		Tcl_AppendResult(interp, "\nh5mdfile: Wrong type in H5_Screate_simple chosen\n",(char *) NULL);
		return TCL_ERROR;
	}
  return TCL_OK;
}
int H5mdfile::H5_Sselect_hyperslab(int argc, char **argv, Tcl_Interp *interp)
{
  /* Select a hyperslab region to extend to the current selected region */
	for(int i=0;i<dataset_rank;i++)
	{
		offset[i] = atoi(argv[3+i]);
	}
	dataspace_id = H5Dget_space (dataset_id);
	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL,dims, NULL);
	return TCL_OK;
}
int H5mdfile::H5_write_value(int argc, char **argv, Tcl_Interp *interp)
{
	/* Read value from Tcl and write it to dataset array */
	// Get array index
	int index=0;
	if(dataset_rank>=1) index+=atoi(argv[4+dataset_rank-0]);
	if(dataset_rank>=2) index+=atoi(argv[4+dataset_rank-1])*dims[dataset_rank-1];else goto label_index;
	if(dataset_rank>=3) index+=atoi(argv[4+dataset_rank-2])*dims[dataset_rank-1]*dims[dataset_rank-2];else goto label_index;
	if(dataset_rank>=4) index+=atoi(argv[4+dataset_rank-3])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3];else goto label_index;
	if(dataset_rank>=5) index+=atoi(argv[4+dataset_rank-4])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4];else goto label_index;
	if(dataset_rank>=6) index+=atoi(argv[4+dataset_rank-5])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4]*dims[dataset_rank-5];else goto label_index;
	if(dataset_rank>=7) index+=atoi(argv[4+dataset_rank-6])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4]*dims[dataset_rank-5]*dims[dataset_rank-6];else goto label_index;
	if(dataset_rank>=8) index+=atoi(argv[4+dataset_rank-7])*dims[dataset_rank-1]*dims[dataset_rank-2]*dims[dataset_rank-3]*dims[dataset_rank-4]*dims[dataset_rank-5]*dims[dataset_rank-6]*dims[dataset_rank-7];else goto label_index;
	label_index:
	// Write single value from Tcl to dataset array
	if(H5Tequal(dataset_type_id, H5T_NATIVE_FLOAT))
	{
		dset_data_float = static_cast<float*>(dset_data);
		dset_data_float[index]=(float)atof(argv[3]);
	}
	if(H5Tequal(dataset_type_id, H5T_NATIVE_DOUBLE))
	{
		dset_data_double = static_cast<double*>(dset_data);
		dset_data_double[index]=atof(argv[3]);
	}
	if(H5Tequal(dataset_type_id, H5T_NATIVE_INT))
	{
		dset_data_int = static_cast<int*>(dset_data);
		dset_data_int[index]=atoi(argv[3]);
	}
	if(H5Tequal(dataset_type_id, H5T_C_S1))
	{
		dset_data_string = static_cast<h5string*>(dset_data);
		strcpy(dset_data_string[index], argv[3]);
	}
	return TCL_OK;
}
int H5mdfile::H5_free_memory(int argc, char **argv, Tcl_Interp *interp)
{
	/* Free allocated memory from dataset array */
	free(dset_data);
	dset_data=NULL;
	return TCL_OK;
}
int H5mdfile::get_dataset_dims(int argc, char **argv, Tcl_Interp *interp)
{
	hid_t new_dataspace = H5Dget_space(dataset_id);	//dataspace handle
	int new_rank_dataset = H5Sget_simple_extent_ndims(new_dataspace);
	hsize_t dims_out[new_rank_dataset];
	H5Sget_simple_extent_dims(new_dataspace, dims_out, NULL);
	char buffer[32+TCL_INTEGER_SPACE];
	for(int i =0;i<new_rank_dataset;i++){
		sprintf(buffer, "%d ", (int)dims_out[i]);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
	}
	return TCL_OK;
}
int H5mdfile::H5_Fflush(int argc, char **argv, Tcl_Interp *interp)
{
	H5Fflush(dataset_id, H5F_SCOPE_LOCAL);
	return TCL_OK;
}
#endif

#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>
#include "parser.hpp"
#include "communication.hpp"
#include "hdf5.h"
#include <stdio.h>

#ifdef H5MD

class H5mdfile
{
public:
   hid_t    file_id; // File id returned from H5Fopen
   hid_t    dataset_id; // Dataset id returned from H5Dopen2 or H5Dcreate2
   hid_t    dataspace_id; // Dataspace id returned from H5Dget_space
   hid_t    dataspace_simple_id; // Dataspace id returned from H5_Screate_simple
   hid_t    dataspace_simple_ext_id; // Dataspace id returned from H5_Screate_simple for extended dataset
   hid_t    dataset_type_id; // Dataset data type id
   hid_t    group_id; // Group id
   hid_t  	prop_id; // Property id returned from H5Pcreate and H5Dget_create_plist

   hsize_t*  dims, *maxdims; // Array of dimensions and maximum dimensions for basic dataset
   hsize_t*  dimsext; // Array of dimensions for extended dataset
   hsize_t*  dimstotal; // Array of dimensions for basic plus extended dataset
   hsize_t*  offset; // Offset for extended dataset in accordance to basic dataset
   hsize_t*  chunk_dims; //  Array of dimensions for chunked dataset

   int 		dataset_rank; // Rank of dataset
   int 		chunk_rank; // Rank of chunked dataset
   int		datatype_size; // Size of data type
   herr_t   status; // Return error value of several functions

   int 		dset_data_size; // Size serial array for the basic dataset data
   int 		dset_data_size_ext; // Size serial array for the extended dataset data
   int		string_maxlength; // Maximal string length


   void*	dset_data; // Dataset array
   int*		dset_data_int; // Dataset array (type: integer)
   float*	dset_data_float; // Dataset array (type: float)
   double*	dset_data_double; // Dataset array (type: double)
   char**	dset_data_string; // Dataset array (type: string)
   int 		dset_data_singlevalue_int; // Single value from basic dataset array (type: integer)
   float 	dset_data_singlevalue_float; // Single value from basic dataset array (type: float)
   double 	dset_data_singlevalue_double; // Single value from basic dataset array (type: double)
   char*	dset_data_singlevalue_string; // Single value from basic dataset array (type: string)


   int H5_Dclose(int argc, char **argv, Tcl_Interp *interp);
   int H5_Dcreate2(int argc, char **argv, Tcl_Interp *interp);
   int H5_Dextend(int argc, char **argv, Tcl_Interp *interp);
   int H5_Dopen2(int argc, char **argv, Tcl_Interp *interp);
   int H5_Dread(int argc, char **argv, Tcl_Interp *interp);
   int H5_Dwrite(int argc, char **argv, Tcl_Interp *interp);
   int H5_Fclose(int argc, char **argv, Tcl_Interp *interp);
   int H5_Fcreate(int argc, char **argv, Tcl_Interp *interp);
   int H5_Fopen(int argc, char **argv, Tcl_Interp *interp);
   int H5_Gclose(int argc, char **argv, Tcl_Interp *interp);
   int H5_Gcreate2(int argc, char **argv, Tcl_Interp *interp);
   int H5_Gopen2(int argc, char **argv, Tcl_Interp *interp);
   int H5_Pclose(int argc, char **argv, Tcl_Interp *interp);
   int H5_Pset_chunk(int argc, char **argv, Tcl_Interp *interp);
   int H5_read_value(int argc, char **argv,Tcl_Interp *interp);
   int H5_Sclose(int argc, char **argv, Tcl_Interp *interp);
   int H5_Screate_simple(int argc, char **argv, Tcl_Interp *interp);
   int H5_Sselect_hyperslab(int argc, char **argv, Tcl_Interp *interp);
   int H5_write_value(int argc, char **argv, Tcl_Interp *interp);
};

H5mdfile h5mdfile;

int tclcommand_h5mdfile(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{

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
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dcreate2 \"\path_to_dset\"\n",(char *) NULL);
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
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Dopen2 \"\path_to_dset\"\n",(char *) NULL);
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
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Fcreate \"\path_to_file\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Fcreate(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Fopen", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Fopen \"\path_to_file\"\n",(char *) NULL);
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
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Gcreate2 \"\path_to_group\"\n",(char *) NULL);
			return TCL_ERROR;
		}
		return h5mdfile.H5_Gcreate2(argc, argv, interp);
	}
	if (!strncmp(argv[1], "H5Gopen2", strlen(argv[1])))
	{
		if(argc!=3)
		{
			Tcl_AppendResult(interp, "\nExpected: h5mdfile H5Gopen2 \"\path_to_group\"\n",(char *) NULL);
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
}

int H5mdfile::H5_Dclose(int argc, char **argv, Tcl_Interp *interp)
{
	// End access to the dataset and release resources used by it
	status = H5Dclose(dataset_id);
	return TCL_OK;
}
int H5mdfile::H5_Dcreate2(int argc, char **argv, Tcl_Interp *interp)
{
   // Create the dataset
	dataset_id = H5Dcreate2(file_id, argv[2], dataset_type_id, dataspace_simple_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
	dataspace_id = H5Dget_space (dataset_id);
	return TCL_OK;
}
int H5mdfile::H5_Dextend(int argc, char **argv, Tcl_Interp *interp)
{
	//Extend dataset to higher dimensions
   dataset_rank=0;
   dset_data_size_ext=1;

   while(!(argv[3+dataset_rank]==NULL))
   {
	   dset_data_size_ext*=atoi(argv[3+dataset_rank]);
	   dataset_rank++;
   }

   dimstotal = new hsize_t[atoi(argv[dataset_rank])];
   dimsext = new hsize_t[atoi(argv[dataset_rank])];

   for(int i=0;i<dataset_rank;i++)
   {
	   dimstotal[i] = atoi(argv[3+i]);
	   if(dimstotal[i]>dims[i])
	   {
		   dimsext[i]=dimstotal[i]-dims[i];
		   dims[i]=dimstotal[i];
	   }
	   else
	   {
		   dimsext[i]=dims[i];
		   dims[i]=dimstotal[i];
	   }
   }

   status = H5Dextend(dataset_id, dimstotal);
   return TCL_OK;
}
int H5mdfile::H5_Dopen2(int argc, char **argv, Tcl_Interp *interp)
{
	// Open an existing dataset
	dataset_id = H5Dopen2(file_id, argv[2], H5P_DEFAULT);

	// Dataset properties
	dataspace_id = H5Dget_space(dataset_id);
	dataset_type_id = H5Dget_type(dataset_id);
	datatype_size = H5Tget_size(dataset_type_id);
	dims = new hsize_t[32];
	maxdims = new hsize_t[32];
	dataset_rank = H5Sget_simple_extent_dims(dataspace_id,dims,maxdims);
	prop_id = H5Dget_create_plist (dataset_id);
	if (H5D_CHUNKED == H5Pget_layout (prop_id))
	  chunk_rank = H5Pget_chunk (prop_id, dataset_rank, chunk_dims);

	// Total length of dataset array as the product of all dimensions
	dset_data_size=1;
	for(int i=0;i<dataset_rank;i++)
	{
	   dset_data_size*=dims[i];
	}

    //Reserve memeory
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
	   dset_data = (char**) malloc(dset_data_size * sizeof(char*));
	   //for (int i = 0; i < dset_data_size; i++)
		   //dset_data[i] = malloc((string_maxlength) * sizeof(char));
	}
	else
	{
		Tcl_AppendResult(interp, "\nh5mdfile: No data type in H5_Dopen2 given\n",(char *) NULL);
		return TCL_ERROR;
	}
   return TCL_OK;
}
int H5mdfile::H5_Dread(int argc, char **argv, Tcl_Interp *interp)
{
	// Read dataset from h5-file and write it to dataset array values
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
	   dset_data = (char**) malloc(dset_data_size * sizeof(char*));
	   //for (int i = 0; i < dset_data_size; i++)
		   //dset_data[i] = malloc((string_maxlength) * sizeof(char));
	}
	else
	{
		Tcl_AppendResult(interp, "\nh5mdfile: No data type in H5_Dread given\n",(char *) NULL);
		return TCL_ERROR;
	}
   /* Read the dataset. */
   status = H5Dread(dataset_id, dataset_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,dset_data);
   return TCL_OK;
}
int H5mdfile::H5_Dwrite(int argc, char **argv, Tcl_Interp *interp)
{
	// Write dataset array to h5-file dataset
	hsize_t* testdims = new hsize_t[32];
	hsize_t* testmaxdims = new hsize_t[32];
	int testdataset_rank = H5Sget_simple_extent_dims(dataspace_id,testdims,testmaxdims);

	if(dimsext==NULL) // No dimension extension used
	{
		status = H5Dwrite(dataset_id, dataset_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
	}
	else // Dimension extension used
	{
		status = H5Dwrite(dataset_id, dataset_type_id, dataspace_simple_ext_id, dataspace_id, H5P_DEFAULT, dset_data);
	}
	return TCL_OK;
}
int H5mdfile::H5_Fclose(int argc, char **argv, Tcl_Interp *interp)
{
	// Close the file
	status = H5Fclose(file_id);
    return TCL_OK;
}
int H5mdfile::H5_Fcreate(int argc, char **argv, Tcl_Interp *interp)
{
   // Create a new file using default properties
   file_id = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Fopen(int argc, char **argv, Tcl_Interp *interp)
{
   // Open an existing file
   file_id = H5Fopen(argv[2], H5F_ACC_RDWR, H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Gclose(int argc, char **argv, Tcl_Interp *interp)
{
	// Close the group
	status = H5Gclose(group_id);
    return TCL_OK;
}
int H5mdfile::H5_Gcreate2(int argc, char **argv, Tcl_Interp *interp)
{
   // Create a new group using default properties
   group_id = H5Gcreate2(file_id, argv[2], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Gopen2(int argc, char **argv, Tcl_Interp *interp)
{
   // Open an existing group
   group_id = H5Gopen2(file_id, argv[2], H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Pclose(int argc, char **argv, Tcl_Interp *interp)
{
	// Terminate access to the data space
	status = H5Pclose(prop_id);
    return TCL_OK;
}
int H5mdfile::H5_Pset_chunk(int argc, char **argv, Tcl_Interp *interp)
{
	// Set properties for chunk
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
	// Read value from dataset array and print it to Tcl
	int index;
	if(dataset_rank==1) index=atoi(argv[3]);
	else if(dataset_rank==2) index=atoi(argv[3])*dims[1] + atoi(argv[4]);
	else if(dataset_rank==3) index=atoi(argv[3])*dims[1]*dims[2] + atoi(argv[4])*dims[2] + atoi(argv[5]);
	else if(dataset_rank==4) index=atoi(argv[3])*dims[1]*dims[2]*dims[3] + atoi(argv[4])*dims[2]*dims[3] + atoi(argv[5])*dims[3] + atoi(argv[6]);
	else if(dataset_rank==5) index=atoi(argv[3])*dims[1]*dims[2]*dims[3]*dims[4] + atoi(argv[4])*dims[2]*dims[3]*dims[4] + atoi(argv[5])*dims[3]*dims[4] + atoi(argv[6])*dims[4] + atoi(argv[7]);


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

		char buffer[string_maxlength];
		Tcl_PrintDouble(interp, (double)dset_data_singlevalue_int, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type_id, H5T_C_S1))
	{
		dset_data_string = static_cast<char**>(dset_data);
		dset_data_singlevalue_string = dset_data_string[index];
		//		char* buffer[string_maxlength];
		//		buffer=dset_data_singlevalue_string;
		Tcl_AppendResult(interp, dset_data_singlevalue_string, (char *)NULL);
		return TCL_OK;
	}
}
int H5mdfile::H5_Sclose(int argc, char **argv, Tcl_Interp *interp)
{
	// Terminate access to the data space
	status = H5Sclose(dataspace_simple_id);
    return TCL_OK;
}
int H5mdfile::H5_Screate_simple(int argc, char **argv, Tcl_Interp *interp)
{
    // Create the data space for the dataset
	if(dimsext==NULL) // No dimension extension used
	{
		dataset_rank=0;
		dset_data_size=1;

		while(!(argv[5+dataset_rank]==NULL))
		{
		   dset_data_size*=atoi(argv[5+dataset_rank]);
		   dataset_rank++;
		}
		dims = new hsize_t[atoi(argv[dataset_rank])];
		maxdims = new hsize_t[atoi(argv[dataset_rank])];
		chunk_dims = new hsize_t[atoi(argv[dataset_rank])];
		for(int i=0;i<dataset_rank;i++)
		{
		   dims[i] = atoi(argv[5+i]);
		   maxdims[i] = H5S_UNLIMITED;
		   chunk_dims[i] = 1;
		}
		dataspace_simple_id = H5Screate_simple(dataset_rank, dims, maxdims);

		/* Data set array*/
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
		   dataset_type_id = H5T_C_S1;
		   dset_data = (char**) malloc(dset_data_size * sizeof(char*));
			//		   for (int i = 0; i < dset_data_size; i++)
			//			   dset_data[i] = malloc((string_maxlength) * sizeof(char));
		}
		else
		{
			Tcl_AppendResult(interp, "\nh5mdfile: Wrong type in H5_Screate_simple chosen\n",(char *) NULL);
			return TCL_ERROR;
		}
	}
	else // Dimension extension used
	{
		dataspace_simple_ext_id = H5Screate_simple(dataset_rank, dimsext, NULL);
		/* Data set array*/
		if(!strncmp(argv[3], "float", strlen(argv[3])))
		{
		   dataset_type_id = H5T_NATIVE_FLOAT;
		   dset_data=(float*) malloc(dset_data_size_ext*sizeof(float));
		   memset(dset_data,0,dset_data_size_ext*sizeof(float));
		}
		else if(!strncmp(argv[3], "double", strlen(argv[3])))
		{
		   dataset_type_id = H5T_NATIVE_DOUBLE;
		   dset_data=(double*) malloc(dset_data_size_ext*sizeof(double));
		   memset(dset_data,0,dset_data_size_ext*sizeof(double));
		}
		else if(!strncmp(argv[3], "int", strlen(argv[3])))
		{
		   dataset_type_id = H5T_NATIVE_INT;
		   dset_data=(int*) malloc(dset_data_size_ext*sizeof(int));
		   memset(dset_data,0,dset_data_size_ext*sizeof(int));
		}
		else if(!strncmp(argv[3], "str", strlen(argv[3])))
		{
			   dataset_type_id = H5T_C_S1;
			   dset_data = (char**) malloc(dset_data_size_ext * sizeof(char*));
				//			   for (int i = 0; i < dset_data_size_ext; i++)
				//				   dset_data[i] = malloc((string_maxlength) * sizeof(char));
		}
		else
		{
			Tcl_AppendResult(interp, "\nh5mdfile: Wrong type in H5_Screate_simple chosen\n",(char *) NULL);
			return TCL_ERROR;
		}
	}
   return TCL_OK;
}
int H5mdfile::H5_Sselect_hyperslab(int argc, char **argv, Tcl_Interp *interp)
{
    // Selects a hyperslab region to extend to the current selected region
	int rank=0;
	while(!(argv[3+rank]==NULL))
	{
		rank++;
	}
	offset = new hsize_t[atoi(argv[rank])];
	for(int i=0;i<rank;i++)
	{
		offset[i] = atoi(argv[3+i]);
	}
	dataspace_id = H5Dget_space (dataset_id);
	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL,dimsext, NULL);
	return TCL_OK;
}
int H5mdfile::H5_write_value(int argc, char **argv, Tcl_Interp *interp)
{
	// Read value from Tcl and write it to dataset array
	int index;

	if(dimsext==NULL) // No dimension extension used
	{
		if(dataset_rank==1) index=atoi(argv[5]);
		else if(dataset_rank==2) index=atoi(argv[5])*dims[1] + atoi(argv[6]);
		else if(dataset_rank==3) index=atoi(argv[5])*dims[1]*dims[2] + atoi(argv[6])*dims[2] + atoi(argv[7]);
		else if(dataset_rank==4) index=atoi(argv[5])*dims[1]*dims[2]*dims[3] + atoi(argv[6])*dims[2]*dims[3] + atoi(argv[7])*dims[3] + atoi(argv[8]);
		else if(dataset_rank==5) index=atoi(argv[5])*dims[1]*dims[2]*dims[3]*dims[4] + atoi(argv[6])*dims[2]*dims[3]*dims[4] + atoi(argv[7])*dims[3]*dims[4] + atoi(argv[8])*dims[4] + atoi(argv[9]);

	}
	else // Dimension extension used
	{
		if(dataset_rank==1) index=atoi(argv[5]);
		else if(dataset_rank==2) index=atoi(argv[5])*dimsext[1] + atoi(argv[6]);
		else if(dataset_rank==3) index=atoi(argv[5])*dimsext[1]*dimsext[2] + atoi(argv[6])*dimsext[2] + atoi(argv[7]);
		else if(dataset_rank==4) index=atoi(argv[5])*dimsext[1]*dimsext[2]*dimsext[3] + atoi(argv[6])*dimsext[2]*dimsext[3] + atoi(argv[7])*dimsext[3] + atoi(argv[8]);
		else if(dataset_rank==5) index=atoi(argv[5])*dimsext[1]*dimsext[2]*dimsext[3]*dimsext[4] + atoi(argv[6])*dimsext[2]*dimsext[3]*dimsext[4] + atoi(argv[7])*dimsext[3]*dimsext[4] + atoi(argv[8])*dimsext[4] + atoi(argv[9]);
	}


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
		dset_data_string = static_cast<char**>(dset_data);
		dset_data_string[index]=argv[3];
	}

	return TCL_OK;
}


#endif


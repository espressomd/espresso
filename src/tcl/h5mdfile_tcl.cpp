#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include "parser.hpp"
#include "communication.hpp"
#include "/usr/include/hdf5.h"

#ifdef H5MD

#include <stdio.h>


class H5mdfile
{
public:
   hid_t       file_id, dataset_id, dataspace_id,dataset_type;
   hsize_t     *dims;
   herr_t      status;
   int *dset_data_int;
   float *dset_data_float;
   double *dset_data_double;
   int dset_data_int_size;
   int H5_Fopen(char **argv);
   int H5_Dopen2(char **argv);
   int H5_Fcreate(char **argv);
   int H5_Screate_simple(char **argv);
   int H5_Dcreate2(char **argv);
   int H5_Dwrite(char **argv);
   int H5_Dread(char **argv);
   int H5_Dclose(char **argv);
   int H5_Sclose(char **argv);
   int H5_Fclose(char **argv);
   int Value_to_Dset_int(char **argv);
   void *Tcl_to_Dset(char* f);
   void *Dset_to_Tcl(char* f);
};

H5mdfile h5mdfile;

int tclcommand_h5mdfile(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
	//Open
	if (!strncmp(argv[1], "H5Fopen", strlen(argv[1])))
	{
	  return h5mdfile.H5_Fopen(argv);
	}
	if (!strncmp(argv[1], "H5Dopen2", strlen(argv[1])))
	{
	  return h5mdfile.H5_Dopen2(argv);
	}
	//Create
	if (!strncmp(argv[1], "H5Fcreate", strlen(argv[1])))
	{
	  return h5mdfile.H5_Fcreate(argv);
	}
	if (!strncmp(argv[1], "H5Screate_simple", strlen(argv[1])))
	{
	  return h5mdfile.H5_Screate_simple(argv);
	}
	if (!strncmp(argv[1], "H5Dcreate2", strlen(argv[1])))
	{
	  return h5mdfile.H5_Dcreate2(argv);
	}
	//Read/Write
	if (!strncmp(argv[1], "H5Dwrite", strlen(argv[1])))
	{
	  return h5mdfile.H5_Dwrite(argv);
	}
	if (!strncmp(argv[1], "H5Dread", strlen(argv[1])))
	{
	  return h5mdfile.H5_Dread(argv);
	}
	if (!strncmp(argv[1], "Value_to_Dset_int", strlen(argv[1])))
	{
	  return h5mdfile.Value_to_Dset_int(argv);
	}
	//Close
	if (!strncmp(argv[1], "H5Dclose", strlen(argv[1])))
	{
	  return h5mdfile.H5_Dclose(argv);
	}
	if (!strncmp(argv[1], "H5Sclose", strlen(argv[1])))
	{
	  return h5mdfile.H5_Sclose(argv);
	}
	if (!strncmp(argv[1], "H5Fclose", strlen(argv[1])))
	{
	  return h5mdfile.H5_Fclose(argv);
	}
}


int H5mdfile::H5_Fopen(char **argv)
{
   /* Open an existing file. */
   file_id = H5Fopen(argv[2], H5F_ACC_RDWR, H5P_DEFAULT);
   return TCL_OK;
};
int H5mdfile::H5_Dopen2(char **argv)
{
   /* Open an existing dataset. */
   dataset_id = H5Dopen2(file_id, argv[2], H5P_DEFAULT);
   return TCL_OK;
};
int H5mdfile::H5_Fcreate(char **argv)
{
   /* Create a new file using default properties. */
   file_id = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   return TCL_OK;
};
int H5mdfile::H5_Screate_simple(char **argv)
{
   /* Create the data space for the dataset. */
   dims = new hsize_t[atoi(argv[3])];
   int number_dims=0;
   dset_data_int_size=1;
   while(!(argv[3+number_dims]==NULL))
   {
	   dset_data_int_size*=atoi(argv[3+number_dims]);
	   dims[number_dims] = atoi(argv[3+number_dims]);
	   number_dims++;
   }
   dataspace_id = H5Screate_simple(number_dims, dims, NULL);
   /* Data set array*/
   dset_data_int=new int[dset_data_int_size];
   memset(dset_data_int,0,dset_data_int_size*sizeof(int));
   return TCL_OK;
};
int H5mdfile::H5_Dcreate2(char **argv)
{
   /* Create the dataset. */
	if(strncmp(argv[3], "float", strlen(argv[2]))) dataset_type = H5T_NATIVE_FLOAT;
	if(strncmp(argv[3], "double", strlen(argv[2]))) dataset_type = H5T_NATIVE_DOUBLE;
	if(strncmp(argv[3], "int", strlen(argv[2]))) dataset_type = H5T_NATIVE_INT;

	dataset_id = H5Dcreate2(file_id, argv[2], dataset_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return TCL_OK;
};
int H5mdfile::H5_Dwrite(char **argv)
{
	/* Write the dataset. */
	if(strncmp(argv[2], "float", strlen(argv[2]))) dataset_type = H5T_NATIVE_FLOAT;
	if(strncmp(argv[2], "double", strlen(argv[2]))) dataset_type = H5T_NATIVE_DOUBLE;
	if(strncmp(argv[2], "int", strlen(argv[2]))) dataset_type = H5T_NATIVE_INT;
	printf("XXXXXXXXX1\n");
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data_int);
	printf("XXXXXXXXX2\n");
	return TCL_OK;
};
int H5mdfile::H5_Dread(char **argv)
{
	/* Read the dataset. */
	if(strncmp(argv[3], "float", strlen(argv[2]))) dataset_type = H5T_NATIVE_FLOAT;
	if(strncmp(argv[3], "double", strlen(argv[2]))) dataset_type = H5T_NATIVE_DOUBLE;
	if(strncmp(argv[3], "int", strlen(argv[2]))) dataset_type = H5T_NATIVE_INT;

	status = H5Dread(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,dset_data_int);
	return TCL_OK;
};
int H5mdfile::H5_Dclose(char **argv)
{
	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);
	return TCL_OK;
};
int H5mdfile::H5_Sclose(char **argv)
{
	/* Terminate access to the data space. */
	status = H5Sclose(dataspace_id);
    return TCL_OK;
};
int H5mdfile::H5_Fclose(char **argv)
{
	/* Close the file. */
	status = H5Fclose(file_id);
    return TCL_OK;
};
int H5mdfile::Value_to_Dset_int(char **argv)
{
	/* Write value to dataset array */
	int index;
	if(sizeof(dims)/sizeof(dims[0])==1) index=atoi(argv[5]);
//	if(sizeof(dims)/sizeof(dims[0])==2) index=atoi(argv[5])*dims[0] + atoi(argv[6]);
//	if(sizeof(dims)/sizeof(dims[0])==3) index=atoi(argv[5])*dims[0]*dims[1] + atoi(argv[6])*dims[0] + atoi(argv[7]);
//	if(sizeof(dims)/sizeof(dims[0])==4) index=atoi(argv[5])*dims[0]*dims[1]*dims[2] + atoi(argv[6])*dims[0]*dims[1] + atoi(argv[7])*dims[0] + atoi(argv[8]);
	dset_data_int[index]=atoi(argv[3]);
	return TCL_OK;
};
void *H5mdfile::Tcl_to_Dset(char* f)
{
	return 0;
};
void *H5mdfile::Dset_to_Tcl(char* f)
{
	return 0;
};
#endif

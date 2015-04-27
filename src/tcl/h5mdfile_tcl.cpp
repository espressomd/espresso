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
   int number_dims;
   herr_t      status;
   //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   void *dset_data_void;
   int *dset_data_int;
   float *dset_data_float;
   double *dset_data_double;
   char **dset_data_string;
   int dset_data_void_size;
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
}
int H5mdfile::H5_Dopen2(char **argv)
{
   /* Open an existing dataset. */
   dataset_id = H5Dopen2(file_id, argv[2], H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Fcreate(char **argv)
{
   /* Create a new file using default properties. */
   file_id = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Screate_simple(char **argv)
{
   /* Create the data space for the dataset. */
   number_dims=0;
   dset_data_void_size=1;

   while(!(argv[4+number_dims]==NULL))
   {
	   dset_data_void_size*=atoi(argv[4+number_dims]);

	   number_dims++;
   }
   dims = new hsize_t[atoi(argv[number_dims])];
   for(int i=0;i<number_dims;i++)
   {
	   dims[i] = atoi(argv[4+i]);
   }

   dataspace_id = H5Screate_simple(number_dims, dims, NULL);

   /* Data set array*/
   if(!strncmp(argv[2], "float", strlen(argv[3])))
   {
	   dataset_type = H5T_NATIVE_FLOAT;
	   dset_data_void=new float[dset_data_void_size];
	   memset(dset_data_void,0,dset_data_void_size*sizeof(float));
   }
   if(!strncmp(argv[2], "double", strlen(argv[3])))
   {
	   dataset_type = H5T_NATIVE_DOUBLE;
	   dset_data_void=new double[dset_data_void_size];
	   memset(dset_data_void,0,dset_data_void_size*sizeof(double));
   }
   if(!strncmp(argv[2], "int", strlen(argv[3])))
   {
	   dataset_type = H5T_NATIVE_INT;
	   dset_data_void=new int[dset_data_void_size];
	   memset(dset_data_void,0,dset_data_void_size*sizeof(int));
   }
   if(!strncmp(argv[2], "str", strlen(argv[3])))
   {

//	   dataset_type = H5T_NATIVE_SCHAR;
//	   dset_data_void=new char*[dset_data_void_size];
	   //TODO memset(dset_data_void,0,dset_data_void_size*sizeof(int));
   }


   return TCL_OK;
}
int H5mdfile::H5_Dcreate2(char **argv)
{
   /* Create the dataset. */
	dataset_id = H5Dcreate2(file_id, argv[2], dataset_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return TCL_OK;
}
int H5mdfile::H5_Dwrite(char **argv)
{
	/* Write the dataset. */
	status = H5Dwrite(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data_void);

	return TCL_OK;
}
int H5mdfile::H5_Dread(char **argv)
{
	/* Read the dataset. */
	status = H5Dread(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,dset_data_void);
	return TCL_OK;
}
int H5mdfile::H5_Dclose(char **argv)
{
	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);
	return TCL_OK;
}
int H5mdfile::H5_Sclose(char **argv)
{
	/* Terminate access to the data space. */
	status = H5Sclose(dataspace_id);
    return TCL_OK;
}
int H5mdfile::H5_Fclose(char **argv)
{
	/* Close the file. */
	status = H5Fclose(file_id);
    return TCL_OK;
}
int H5mdfile::Value_to_Dset_int(char **argv)
{
	/* Write value to dataset array */
	int index;
	if(number_dims==1) index=atoi(argv[5]);
	if(number_dims==2) index=atoi(argv[5])*dims[1] + atoi(argv[6]);
	if(number_dims==3) index=atoi(argv[5])*dims[1]*dims[2] + atoi(argv[6])*dims[2] + atoi(argv[7]);
	if(number_dims==4) index=atoi(argv[5])*dims[1]*dims[2]*dims[3] + atoi(argv[6])*dims[2]*dims[3] + atoi(argv[7])*dims[3] + atoi(argv[8]);
	if(number_dims==5) index=atoi(argv[5])*dims[1]*dims[2]*dims[3]*dims[4] + atoi(argv[6])*dims[2]*dims[3]*dims[4] + atoi(argv[7])*dims[3]*dims[4] + atoi(argv[8])*dims[4] + atoi(argv[9]);


	if(dataset_type==H5T_NATIVE_FLOAT)
	{
		dset_data_float = static_cast<float*>(dset_data_void);
		dset_data_float[index]=(float)atof(argv[3]);
	}
	if(dataset_type==H5T_NATIVE_DOUBLE)
	{
		dset_data_double = static_cast<double*>(dset_data_void);
		dset_data_double[index]=atof(argv[3]);
	}
	if(dataset_type==H5T_NATIVE_INT)
	{
		dset_data_int = static_cast<int*>(dset_data_void);
		dset_data_int[index]=atoi(argv[3]);
	}
	if(dataset_type==H5T_NATIVE_INT)
	{
//		dset_data_int = static_cast<int*>(dset_data_void);
//		dset_data_int[index]=atoi(argv[3]);
	}
	printf("XXX1111\n");

	return TCL_OK;
}

#endif

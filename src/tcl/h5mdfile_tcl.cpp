#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include "parser.hpp"
#include "communication.hpp"
#include "/usr/include/hdf5.h"//XXXXXXXXXXXX

#ifdef H5MD

#include <stdio.h>
#include <iostream>
using namespace std;

class H5mdfile
{
public:
   hid_t    file_id, dataset_id, dataspace_id,dataset_type;
   hsize_t  *dims, *maxdims;
   int 		dataset_rank;
   int		datatype_size;
   herr_t   status;
   int 		dset_data_size;


   void 	*dset_data;
   int 		*dset_data_int;
   float 	*dset_data_float;
   double 	*dset_data_double;
   char 	**dset_data_string;
   int 		dset_data_singlevalue_int;
   float 	dset_data_singlevalue_float;
   double 	dset_data_singlevalue_double;

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
   int H5_write_value(char **argv);
   int H5_read_value(char **argv,Tcl_Interp *interp);
   void *Tcl_to_Dset(char* f);
   void *Dset_to_Tcl(char* f);
};

H5mdfile h5mdfile;

int tclcommand_h5mdfile(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
	//Open/Close
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
	if (!strncmp(argv[1], "H5Fclose", strlen(argv[1])))
	{
		return h5mdfile.H5_Fclose(argv);
	}
	if (!strncmp(argv[1], "H5Dclose", strlen(argv[1])))
	{
		return h5mdfile.H5_Dclose(argv);
	}
	if (!strncmp(argv[1], "H5Sclose", strlen(argv[1])))
	{
		return h5mdfile.H5_Sclose(argv);
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
	if (!strncmp(argv[1], "H5_write_value", strlen(argv[1])))
	{
		return h5mdfile.H5_write_value(argv);
	}
	if (!strncmp(argv[1], "H5_read_value", strlen(argv[1])))
	{
		return h5mdfile.H5_read_value(argv,interp);
	}
}

//Open/Close
int H5mdfile::H5_Fopen(char **argv)
{
   /* Open an existing file. */
   file_id = H5Fopen(argv[2], H5F_ACC_RDWR, H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Dopen2(char **argv)
{printf("dataset_type:%i\n",dataset_type);
   /* Open an existing dataset. */
   dataset_id = H5Dopen2(file_id, argv[2], H5P_DEFAULT);
   /* Dataset properties */
   dataspace_id = H5Dget_space(dataset_id);
   dataset_type = H5Dget_type(dataset_id);
   datatype_size = H5Tget_size(dataset_type);
   dims = new hsize_t[5];//TODO nicht manuell zuornden
   dataset_rank = H5Sget_simple_extent_dims(dataspace_id,dims,maxdims);
   /* Total size of dataset */
   dset_data_size=1;
   for(int i=0;i<dataset_rank;i++)
   {
	   dset_data_size*=dims[i];
   }
   return TCL_OK;
}
int H5mdfile::H5_Fclose(char **argv)
{
	/* Close the file. */
	status = H5Fclose(file_id);
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
//Create
int H5mdfile::H5_Fcreate(char **argv)
{
   /* Create a new file using default properties. */
   file_id = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   return TCL_OK;
}
int H5mdfile::H5_Screate_simple(char **argv)
{
   /* Create the data space for the dataset. */
   dataset_rank=0;
   dset_data_size=1;

   while(!(argv[4+dataset_rank]==NULL))
   {
	   dset_data_size*=atoi(argv[4+dataset_rank]);
	   dataset_rank++;
   }
   dims = new hsize_t[atoi(argv[dataset_rank])];
   for(int i=0;i<dataset_rank;i++)
   {
	   dims[i] = atoi(argv[4+i]);
   }

   dataspace_id = H5Screate_simple(dataset_rank, dims, NULL);

   /* Data set array*/
   if(!strncmp(argv[2], "float", strlen(argv[3])))
   {
	   dataset_type = H5T_NATIVE_FLOAT;
	   dset_data=new float[dset_data_size];
	   memset(dset_data,0,dset_data_size*sizeof(float));
   }
   if(!strncmp(argv[2], "double", strlen(argv[3])))
   {
	   dataset_type = H5T_NATIVE_DOUBLE;
	   dset_data=new double[dset_data_size];
	   memset(dset_data,0,dset_data_size*sizeof(double));
   }
   if(!strncmp(argv[2], "int", strlen(argv[3])))
   {
	   dataset_type = H5T_NATIVE_INT;
	   dset_data=new int[dset_data_size];
	   memset(dset_data,0,dset_data_size*sizeof(int));
   }
   if(!strncmp(argv[2], "str", strlen(argv[3])))
   {
  	   dataset_type = H5T_NATIVE_SCHAR;
	   //dset_data=new char*[dset_data_size];
	   //TODO memset(dset_data,0,dset_data_size*sizeof(int));
   }

   return TCL_OK;
}
int H5mdfile::H5_Dcreate2(char **argv)
{
   /* Create the dataset. */
	dataset_id = H5Dcreate2(file_id, argv[2], dataset_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	return TCL_OK;
}
//Read/Write
int H5mdfile::H5_Dwrite(char **argv)
{
	/* Write the dataset. */
	status = H5Dwrite(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

	return TCL_OK;
}
int H5mdfile::H5_Dread(char **argv)
{
   /* Data set array*/
   if(H5Tequal(dataset_type, H5T_NATIVE_FLOAT))
   {
	   dataset_type = H5T_NATIVE_FLOAT;
	   dset_data=new float[dset_data_size];
	   memset(dset_data,0,dset_data_size*sizeof(float));
   }
   if(H5Tequal(dataset_type, H5T_NATIVE_DOUBLE))
   {
	   dataset_type = H5T_NATIVE_DOUBLE;
	   dset_data=new double[dset_data_size];
	   memset(dset_data,0,dset_data_size*sizeof(double));
   }
   if(H5Tequal(dataset_type, H5T_NATIVE_INT))
   {
	   dataset_type = H5T_NATIVE_INT;
	   dset_data=new int[dset_data_size];
	   memset(dset_data,0,dset_data_size*sizeof(int));
   }
   if(H5Tequal(dataset_type, H5T_NATIVE_SCHAR))
   {
	   //dataset_type = H5T_NATIVE_SCHAR;
	   //dset_data=new char*[dset_data_size];
	   //TODO memset(dset_data,0,dset_data_size*sizeof(int));
   }
   /* Read the dataset. */
   status = H5Dread(dataset_id, dataset_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,dset_data);
   return TCL_OK;
}
int H5mdfile::H5_write_value(char **argv)
{
	/* Write value to dataset array */
	int index;
	if(dataset_rank==1) index=atoi(argv[5]);
	if(dataset_rank==2) index=atoi(argv[5])*dims[1] + atoi(argv[6]);
	if(dataset_rank==3) index=atoi(argv[5])*dims[1]*dims[2] + atoi(argv[6])*dims[2] + atoi(argv[7]);
	if(dataset_rank==4) index=atoi(argv[5])*dims[1]*dims[2]*dims[3] + atoi(argv[6])*dims[2]*dims[3] + atoi(argv[7])*dims[3] + atoi(argv[8]);
	if(dataset_rank==5) index=atoi(argv[5])*dims[1]*dims[2]*dims[3]*dims[4] + atoi(argv[6])*dims[2]*dims[3]*dims[4] + atoi(argv[7])*dims[3]*dims[4] + atoi(argv[8])*dims[4] + atoi(argv[9]);


	if(H5Tequal(dataset_type, H5T_NATIVE_FLOAT))
	{
		dset_data_float = static_cast<float*>(dset_data);
		dset_data_float[index]=(float)atof(argv[3]);
	}
	if(H5Tequal(dataset_type, H5T_NATIVE_DOUBLE))
	{
		dset_data_double = static_cast<double*>(dset_data);
		dset_data_double[index]=atof(argv[3]);
	}
	if(H5Tequal(dataset_type, H5T_NATIVE_INT))
	{
		dset_data_int = static_cast<int*>(dset_data);
		dset_data_int[index]=atoi(argv[3]);
	}
	if(H5Tequal(dataset_type, H5T_NATIVE_SCHAR))
	{
//		dset_data_int = static_cast<int*>(dset_data);
//		dset_data_int[index]=atoi(argv[3]);
	}

	return TCL_OK;
}
int H5mdfile::H5_read_value(char **argv,Tcl_Interp *interp)
{
	/* Write value to dataset array */
	int index;
	if(dataset_rank==1) index=atoi(argv[5]);
	if(dataset_rank==2) index=atoi(argv[5])*dims[1] + atoi(argv[6]);
	if(dataset_rank==3) index=atoi(argv[5])*dims[1]*dims[2] + atoi(argv[6])*dims[2] + atoi(argv[7]);
	if(dataset_rank==4) index=atoi(argv[5])*dims[1]*dims[2]*dims[3] + atoi(argv[6])*dims[2]*dims[3] + atoi(argv[7])*dims[3] + atoi(argv[8]);
	if(dataset_rank==5) index=atoi(argv[5])*dims[1]*dims[2]*dims[3]*dims[4] + atoi(argv[6])*dims[2]*dims[3]*dims[4] + atoi(argv[7])*dims[3]*dims[4] + atoi(argv[8])*dims[4] + atoi(argv[9]);


	if(H5Tequal(dataset_type, H5T_NATIVE_FLOAT))
	{
		dset_data_float = static_cast<float*>(dset_data);
		dset_data_singlevalue_float = dset_data_float[index];

		char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
		Tcl_PrintDouble(interp, (double)dset_data_singlevalue_float, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type, H5T_NATIVE_DOUBLE))
	{
		dset_data_double = static_cast<double*>(dset_data);
		dset_data_singlevalue_double = dset_data_double[index];

		char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
		Tcl_PrintDouble(interp, dset_data_singlevalue_double, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type, H5T_NATIVE_INT))
	{
		dset_data_int = static_cast<int*>(dset_data);
		dset_data_singlevalue_int = dset_data_int[index];

		char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
		Tcl_PrintDouble(interp, (double)dset_data_singlevalue_int, buffer);
		Tcl_AppendResult(interp, buffer, (char *)NULL);
		return TCL_OK;
	}
	if(H5Tequal(dataset_type, H5T_NATIVE_SCHAR))
	{
//		dset_data_int = static_cast<int*>(dset_data);
//		dset_data_int[index]=atoi(argv[3]);
	}


}
#endif

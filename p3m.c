/** \file p3m.c
 *
 *  Calculation of long range part of the coulomb interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "grid.h"
#include "particle_data.h"
#include "utils.h"
#include "communication.h"
#include "p3m.h"


/************************************************
 * data types
 ************************************************/

/** Structure for local mesh parameters. */
typedef struct {
  /* local mesh characterization. */
  int dim[3];       /** dimension (size) of local mesh. */
  int l_ind[3];     /** coordinate index of lower left corner of the 
			local mesh in the global mesh. */
  double left[3];   /** coordinate of lower left corner of local mesh */
  double off[3];    /** offset of the first local mesh point (lower left
			corner from my_left[3]. */
  int inner[3];        /** dimension of mesh inside node domain. */
  int left_margin[3];  /** number of left margin mesh points. */
  int rigth_margin[3]; /** number of right margin mesh points. */


} local_mesh;

/** Structure for send/recv meshs. */
typedef struct {
  int s_size[6];     /** sizes for send buffers. */
  int s_ld[6][3];    /** left down corners of sub meshs to send. */
  int s_ur[6][3];    /** up right corners of sub meshs to send. */
  int r_size[6];     /** sizes for recv buffers. */
  int r_ld[6][3];    /** left down corners of sub meshs to recv. */
  int r_ur[6][3];    /** up right corners of sub meshs to recv. */
} send_mesh;

/************************************************
 * variables
 ************************************************/

/** p3m parameters. */
p3m_struct p3m;

/** local mesh. */
local_mesh lm;
/** send/recv mesh sizes */
send_mesh srm;

/************************************************
 * privat functions
 ************************************************/

/************************************************
 * public functions
 ************************************************/

void   P3M_init()
{
  MPI_Barrier(MPI_COMM_WORLD);   
  fprintf(stderr,"%d: P3M Parameters: Bjerrum = %f, alpha = %f, r_cut = %f\n"
	  ,this_node,p3m.bjerrum,p3m.alpha,p3m.r_cut);
  fprintf(stderr,"        mesh (%d, %d, %d) CAO = %d\n"
	  ,p3m.mesh[0],p3m.mesh[1],p3m.mesh[2],p3m.cao);
  fprintf(stderr,"        meshoff (%f, %f, %f) epsiolon = %f\n"
	  ,p3m.mesh_off[0],p3m.mesh_off[1],p3m.mesh_off[2],p3m.epsilon);
  
  MPI_Barrier(MPI_COMM_WORLD);   
}

void   P3M_perform()
{

}

void   P3M_exit()
{

}

int bjerrum_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "Bjerrum length must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.bjerrum = data;
  mpi_bcast_parameter(FIELD_BJERRUM);
  return (TCL_OK);
}

int p3malpha_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0.0 || data > 1.0) {
    Tcl_AppendResult(interp, "P3M alpha must be in interval [0.0,1.0]", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.alpha = data;
  mpi_bcast_parameter(FIELD_P3M_ALPHA);
  return (TCL_OK);
}

int p3mrcut_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "P3M r_cut must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.r_cut = data;
  mpi_bcast_parameter(FIELD_P3M_RCUT);
  return (TCL_OK);
}

int p3mmesh_callback(Tcl_Interp *interp, void *_data)
{
  int i, *data = (int *)_data;
  if (data[0] < 0 || data[1] < 0 || data[2] < 0) {
    Tcl_AppendResult(interp, "mesh sizes must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  for(i=0;i<3;i++) p3m.mesh[i] = data[i];
  mpi_bcast_parameter(FIELD_P3M_MESH);
  return (TCL_OK);
}

int p3mcao_callback(Tcl_Interp *interp, void *_data)
{
  int data = *(int *)_data;
  if (data < 0 || data > 7) {
    Tcl_AppendResult(interp, "P3M CAO must be in interval [0,7]", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.cao = data;
  mpi_bcast_parameter(FIELD_P3M_CAO);
  return (TCL_OK);
}

int p3mepsilon_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  p3m.epsilon = data;
  mpi_bcast_parameter(FIELD_P3M_EPSILON);
  return (TCL_OK);
}

int p3mmeshoff_callback(Tcl_Interp *interp, void *_data)
{
  int i;
  double *data = (double *)_data;
  if (data[0] < 0.0 || data[0] > 1.0) {
    Tcl_AppendResult(interp, "P3M mesh offsets must be in interval [0.0,1.0]", (char *) NULL);
    return (TCL_ERROR);
  }  for(i=0;i<3;i++) p3m.mesh_off[i] = data[i];
  mpi_bcast_parameter(FIELD_P3M_MESH_OFF);
  return (TCL_OK);
}


// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file p3m.c  P3M algorithm for long range coulomb interaction.
 *
 *  For more information about the p3m algorithm,
 *  see \ref p3m.h "p3m.h"
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "integrate.h"

#include "global.h"
#include "debug.h"
#include "grid.h"
#include "integrate.h"
#include "particle_data.h"
#include "utils.h"
#include "communication.h"
#include "fft.h"
#ifdef NPT
#include "pressure.h"
#endif
#include "p3m.h"
#include "thermostat.h"
#include "cells.h"
#include "tuning.h"

#ifdef ELECTROSTATICS

/************************************************
 * DEFINES
 ************************************************/
   
/** increment size of charge assignment fields. */
#define CA_INCREMENT 50       

/* MPI tags for the p3m communications: */
/** Tag for communication in P3M_init() -> send_calc_mesh(). */
#define REQ_P3M_INIT   200
/** Tag for communication in gather_fft_grid(). */
#define REQ_P3M_GATHER 201
/** Tag for communication in spread_force_grid(). */
#define REQ_P3M_SPREAD 202

/************************************************
 * data types
 ************************************************/

/** Structure for local mesh parameters. */
typedef struct {
  /* local mesh characterization. */
  int dim[3];       /** dimension (size) of local mesh. */
  int size;         /** number of local mesh points. */
  int ld_ind[3];    /** index of lower left corner of the 
			local mesh in the global mesh. */
  double ld_pos[3]; /** position of the first local mesh point. */
  int inner[3];     /** dimension of mesh inside node domain. */
  int in_ld[3];     /** inner left down grid point */
  int in_ur[3];     /** inner up right grid point + (1,1,1)*/
  int margin[6];    /** number of margin mesh points. */
  int r_margin[6];  /** number of margin mesh points from neighbour nodes */
} local_mesh;

/** Structure for send/recv meshs. */
typedef struct {
  int s_dim[6][3];   /** dimension of sub meshs to send. */
  int s_ld[6][3];    /** left down corners of sub meshs to send. */
  int s_ur[6][3];    /** up right corners of sub meshs to send. */
  int s_size[6];     /** sizes for send buffers. */
  int r_dim[6][3];   /** dimensionof sub meshs to recv. */
  int r_ld[6][3];    /** left down corners of sub meshs to recv. */
  int r_ur[6][3];    /** up right corners of sub meshs to recv. */
  int r_size[6];     /** sizes for recv buffers. */
  int max;           /** maximal size for send/recv buffers. */
} send_mesh;

/************************************************
 * variables
 ************************************************/

p3m_struct p3m = { 0.0, 0.0, 
		   {0,0,0}, {P3M_MESHOFF, P3M_MESHOFF, P3M_MESHOFF}, 
		   0, P3M_N_INTERPOL, 0.0, P3M_EPSILON, 
		   {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}, 0.0, 0.0 };

/** number of charged particles (only on master node). */
int p3m_sum_qpart=0;
/** Sum of square of charges (only on master node). */
double p3m_sum_q2 = 0.0;
/** square of sum of charges (only on master node). */
double p3m_square_sum_q = 0.0;


/** local mesh. */
static local_mesh lm;
/** send/recv mesh sizes */
static send_mesh  sm;

/** size of linear array for local CA/FFT mesh . */
int    ca_mesh_size;
/** real space mesh (local) for CA/FFT.*/
double *rs_mesh = NULL;
/** k space mesh (local) for k space calculation and FFT.*/
double *ks_mesh = NULL;

/** Field to store grid points to send. */
double *send_grid = NULL; 
/** Field to store grid points to recv */
double *recv_grid = NULL;
/** Allocation size of send_grid and recv_grid. */
int send_recv_grid_size=0;

/** interpolation of the charge assignment function. */
double *int_caf[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
/** position shift for calc. of first assignment mesh point. */
double pos_shift;

/** help variable for calculation of aliasing sums */
double *meshift = NULL;
/** Spatial differential operator in k-space. We use an i*k differentiation. */
double *d_op = NULL;
/** Optimal influence function (k-space) */
double *g = NULL;

/** number of charged particles on the node. */
int ca_num=0;
/** Charge fractions for mesh assignment. */
double *ca_frac = NULL;
/** first mesh point for charge assignment. */
int *ca_fmp = NULL;

/** number of permutations in k_space */
int ks_pnum;

/** \name Private Functions */
/************************************************************/
/*@{*/

/** Initializes the (inverse) mesh constant \ref p3m_struct::a (\ref p3m_struct::ai) 
    and the cutoff for charge assignment \ref p3m_struct::cao_cut, which has to be
    done by \ref P3M_init once and by \ref P3M_scaleby_box_l whenever the \ref box_l changed.
*/
void P3M_init_a_ai_cao_cut();

/** checks for correctness of the cao_cut, necessary when the box length changes */
int P3M_sanity_checks_boxl();

/** Calculate the spacial position of the left down mesh point of the local mesh, to be
    stored in \ref local_mesh::ld_pos; function called by \ref calc_local_ca_mesh once
    and by \ref P3M_scaleby_box_l whenever the \ref box_l changed. */
void calc_lm_ld_pos();

/** Calculates the dipole term */
double calc_dipole_term(int force_flag, int energy_flag);

/** Calculates properties of the local FFT mesh for the 
    charge assignment process. */
void calc_local_ca_mesh();

/** print local mesh content. 
    \param l local mesh structure.
*/
void p3m_print_local_mesh(local_mesh l);

/** Calculates the properties of the send/recv sub-meshes of the local FFT mesh. 
 *  In order to calculate the recv sub-meshes there is a communication of 
 *  the margins between neighbouring nodes. */ 
void calc_send_mesh();

/** print send mesh content. 
 *  \param sm send mesh structure.
*/
void p3m_print_send_mesh(send_mesh sm);

/** Gather FFT grid.
 *  After the charge assignment Each node needs to gather the
 *  information for the FFT grid in his spatial domain.
 */
void gather_fft_grid();

/** Spread force grid.
 *  After the k-space calculations each node needs to get all force
 *  information to reassigne the forces from the grid to the
 *  particles.
 */
void spread_force_grid();

/** realloc charge assignment fields. */
void realloc_ca_fields(int newsize);

/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  ouput array with dimension dim[3] at start position start[3].  
 *
 *  \param in          Pointer to first element of input block data.
 *  \param out         Pointer to first element of output grid.
 *  \param start       Start position of block in output grid.
 *  \param size        Dimensions of the block
 *  \param dim         Dimensions of the output grid.
*/
void add_block(double *in, double *out, int start[3], int size[3], int dim[3]);

/** Interpolates the P-th order charge assignment function from
 * Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
 * are also tabulated in Deserno/Holm. */
void interpolate_charge_assignment_function();

/** shifts the mesh points by mesh/2 */
void calc_meshift();

/** Calculates the Fourier transformed differential operator.  
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *           i.e. the prefactor i*2*PI/L is missing! */
void calc_differential_operator();
/** Calculates the optimal influence function of Hockney and Eastwood. 
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft_plan[3].mesh and fft_plan[3].start).

 *  See also: Hockney/Eastwood 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  Deserno/Holm. */
void calc_influence_function();

/** Calculates the aliasing sums for the optimal influence function.
 *
 * Calculates the aliasing sums in the nominator and denominator of
 * the expression for the optimal influence function (see
 * Hockney/Eastwood: 8-22, p. 275).  
 *
 * \param  n           n-vector for which the aliasing sum is to be performed.
 * \param  nominator   aliasing sums in the nominator.
 * \return denominator aliasing sum in the denominator
 */
MDINLINE double perform_aliasing_sums(int n[3], double nominator[3]);
/*@}*/


/** \name P3M Tuning Functions (private)*/
/************************************************************/
/*@{*/

/** Calculates the real space contribution to the rms error in the force (as described 
   by Kolafa and Perram). 
   \param box_size size of cubic simulation box. 
   \param prefac   Prefactor of coulomb interaction.
   \param r_cut_iL rescaled real space cutoff for p3m method.
   \param n_c_part number of charged particles in the system.
   \param sum_q2   sum of square of charges in the system
   \param alpha_L  rescaled ewald splitting parameter.
   \return real space error
*/
double P3M_real_space_error(double box_size, double prefac, double r_cut_iL, 
			    int n_c_part, double sum_q2, double alpha_L);

/** Calculate the analytic expression of the error estimate for the
    P3M method in the book of Hockney and Eastwood (Eqn. 8.23) in
    order to obtain the rms error in the force for a system of N
    randomly distributed particles in a cubic box (k space part).
    \param box_size size of cubic simulation box. 
    \param prefac   Prefactor of coulomb interaction.
    \param mesh     number of mesh points in one direction.
    \param cao      charge assignment order.
    \param n_c_part number of charged particles in the system.
    \param sum_q2   sum of square of charges in the system
    \param alpha_L  rescaled ewald splitting parameter.
    \return reciprocal (k) space error
*/

double P3M_k_space_error(double box_size, double prefac, int mesh, 
			 int cao, int n_c_part, double sum_q2, double alpha_L);

/** One of the aliasing sums used by \ref P3M_k_space_error. 
    (fortunately the one which is most important (because it converges
    most slowly, since it is not damped exponentially)) can be
    calculated analytically. The result (which depends on the order of
    the spline interpolation) can be written as an even trigonometric
    polynomial. The results are tabulated here (The employed formula
    is Eqn. 7.66 in the book of Hockney and Eastwood). */
double analytic_cotangent_sum(int n, int mesh, int cao);

/** aliasing sum used by \ref P3M_k_space_error. */
void P3M_tune_aliasing_sums(int nx, int ny, int nz, 
			    int mesh, int cao, double alpha_L , 
			    double *alias1, double *alias2);
/*@}*/


/************************************************************/

int printP3MToResult(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_PrintDouble(interp, p3m.r_cut, buffer);
  Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.mesh[0]);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",p3m.cao);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.accuracy, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  Tcl_AppendResult(interp, "} {coulomb epsilon ", (char *) NULL);
  if (p3m.epsilon == P3M_EPSILON_METALLIC)
    Tcl_AppendResult(interp, " metallic ", (char *) NULL);
  else {
    Tcl_PrintDouble(interp, p3m.epsilon, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  }
  sprintf(buffer,"%d",p3m.inter);
  Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.mesh_off[0], buffer);
  Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.mesh_off[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, p3m.mesh_off[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

void p3m_set_tune_params(double r_cut, int mesh, int cao,
			 double alpha, double accuracy)
{
  if (r_cut >= 0) {
    p3m.r_cut    = r_cut;
    p3m.r_cut_iL = r_cut*box_l_i[0];
  }

  if (mesh >= 0)
    p3m.mesh[2] = p3m.mesh[1] = p3m.mesh[0] = mesh;

  if (cao >= 0)
    p3m.cao = cao;

  if (alpha >= 0) {
    p3m.alpha   = alpha;
    p3m.alpha_L = alpha*box_l[0];
  }

  if (accuracy >= 0)
    p3m.accuracy = accuracy;

  mpi_bcast_coulomb_params();
}

int p3m_set_params(double r_cut, int mesh, int cao,
		   double alpha, double accuracy)
{
  if(r_cut < 0)
    return -1;

  if(mesh < 0)
    return -2;

  if(cao < 1 || cao > 7 || cao > mesh)
    return -3;

  p3m.r_cut    = r_cut;
  p3m.r_cut_iL = r_cut*box_l_i[0];
  p3m.mesh[2]  = p3m.mesh[1] = p3m.mesh[0] = mesh;
  p3m.cao      = cao;

  if (alpha > 0) {
    p3m.alpha   = alpha;
    p3m.alpha_L = alpha*box_l[0];
  }
  else
    if (alpha != -1.0)
      return -4;

  if (accuracy >= 0)
    p3m.accuracy = accuracy;
  else
    if (accuracy != -1.0)
      return -5;

  mpi_bcast_coulomb_params();

  return 0;
}

int p3m_set_mesh_offset(double x, double y, double z)
{
  if(x < 0.0 || x > 1.0 ||
     y < 0.0 || y > 1.0 ||
     z < 0.0 || z > 1.0 )
    return TCL_ERROR;

  p3m.mesh_off[0] = x;
  p3m.mesh_off[1] = y;
  p3m.mesh_off[2] = z;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

int p3m_set_eps(double eps)
{
  p3m.epsilon = eps;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

int p3m_set_ninterpol(int n)
{
  if (n < 0)
    return TCL_ERROR;

  p3m.inter = n;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

int inter_parse_p3m_tune_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int mesh, cao;
  double r_cut, accuracy;

  mesh = cao = -1;
  r_cut = accuracy = -1.0;

  while(argc > 0) {

    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(r_cut) && r_cut > 0)) {
	Tcl_AppendResult(interp, "r_cut expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("mesh")) {
      if(! (argc > 1 && ARG1_IS_I(mesh) && mesh >= -1)) {
	Tcl_AppendResult(interp, "mesh expects an integer >= -1",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("cao")) {
      if(! (argc > 1 && ARG1_IS_I(cao) && cao >= -1 && cao < 7)) {
	Tcl_AppendResult(interp, "cao expects an integer between -1 and 7",
			 (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("accuracy")) {
      if(! (argc > 1 && ARG1_IS_D(accuracy) && accuracy > 0)) {
	Tcl_AppendResult(interp, "accuracy expects an positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } else {
      Tcl_AppendResult(interp, "Unknown p3m tune parameter \"",argv[0],"\"",
                       (char *) NULL);
      return TCL_ERROR;
    }
    
    argc -= 2;
    argv += 2;
  }

  p3m_set_tune_params(r_cut, mesh, cao, -1.0, accuracy);

  if(P3M_tune_parameters(interp) == TCL_ERROR) 
    return TCL_ERROR;
  
  return TCL_OK;
}

int inter_parse_p3m(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha, accuracy = -1.0;
  int mesh, cao, i;

  coulomb.method = COULOMB_P3M;
    
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 ||
     PERIODIC(1) == 0 ||
     PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with Coulomb P3M",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if (argc < 1) {
    Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> p3m tune | <r_cut> <mesh> <cao> [<alpha> [accuracy]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    Tcl_AppendResult(interp, "Node grid not suited for Coulomb P3M. Node grid must be sorted, largest first.", (char *) NULL);
    return TCL_ERROR;  
  }

  if (ARG0_IS_S("tune"))
    return inter_parse_p3m_tune_params(interp, argc-1, argv+1);
      
  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc < 3 || argc > 5) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb <bjerrum> p3m <r_cut> <mesh> <cao> [<alpha> [accuracy]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if((! ARG_IS_I(1, mesh)) || (! ARG_IS_I(2, cao))) {
    Tcl_AppendResult(interp, "integer expected", (char *) NULL);
    return TCL_ERROR;
  }
	
  if(argc > 3) {
    if(! ARG_IS_D(3, alpha))
      return TCL_ERROR;
  }
  else {
    Tcl_AppendResult(interp, "Automatic p3m tuning not implemented.",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(argc > 4) {
    if(! ARG_IS_D(4, accuracy)) {
      Tcl_AppendResult(interp, "double expected", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if ((i = p3m_set_params(r_cut, mesh, cao, alpha, accuracy)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
      break;
    case -3:
      Tcl_AppendResult(interp, "cao must be between 1 and 7 and less than mesh",
		       (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "accuracy must be positive", (char *) NULL);
      break;
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;

  }

  return TCL_OK;
}

int inter_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int i; double d1, d2, d3;

  Tcl_ResetResult(interp);

  while (argc > 0) {
    /* p3m parameter: inter */
    if (ARG0_IS_S("n_interpol")) {
      
      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG1_IS_I(i)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 INTEGER parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (p3m_set_ninterpol(i) == TCL_ERROR) {
	Tcl_AppendResult(interp, argv[0], " argument must be positive",
			 (char *) NULL);
	return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;
    }
    
    /* p3m parameter: mesh_off */
    else if (ARG0_IS_S("mesh_off")) {
      
      if(argc < 4) {
	Tcl_AppendResult(interp, argv[0], " needs 3 parameters",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
      if ((! ARG_IS_D(1, d1)) ||
	  (! ARG_IS_D(2, d2)) ||
	  (! ARG_IS_D(3, d3)))
	{
	  Tcl_AppendResult(interp, argv[0], " needs 3 DOUBLE parameters",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      if (p3m_set_mesh_offset(d1, d2 ,d3) == TCL_ERROR)
	{
	  Tcl_AppendResult(interp, argv[0], " parameters have to be between 0.0 an 1.0",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      argc -= 4;
      argv += 4;
    }
    
    /* p3m parameter: epsilon */
    else if(ARG0_IS_S( "epsilon")) {

      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }

      if (ARG1_IS_S("metallic")) {
	d1 = P3M_EPSILON_METALLIC;
      }
      else {
	if (! ARG1_IS_D(d1)) {
	  Tcl_AppendResult(interp, argv[0], " needs 1 DOUBLE parameter or \"metallic\"",
			   (char *) NULL);
	  return TCL_ERROR;
	}
	
	if (p3m_set_eps(d1) == TCL_ERROR) {
	  Tcl_AppendResult(interp, argv[0], " There is no error msg yet!",
			   (char *) NULL);
	  return TCL_ERROR;
	}
      }

      argc -= 2;
      argv += 2;	    
    }
    else {
      Tcl_AppendResult(interp, "Unknown coulomb p3m parameter: \"",argv[0],"\"",(char *) NULL);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}

void P3M_scaleby_box_l() {
  p3m.r_cut = p3m.r_cut_iL* box_l[0];
  p3m.alpha = p3m.alpha_L * box_l_i[0];
  P3M_init_a_ai_cao_cut();
  calc_lm_ld_pos();
  P3M_sanity_checks_boxl();
}

void calc_lm_ld_pos() {
  int i; 
  /* spacial position of left down mesh point */
  for(i=0;i<3;i++) {
    lm.ld_pos[i] = (lm.ld_ind[i]+ p3m.mesh_off[i])*p3m.a[i];
  }
}

void P3M_init_a_ai_cao_cut() {
  int i;
  for(i=0;i<3;i++) {
    p3m.ai[i]      = (double)p3m.mesh[i]/box_l[i]; 
    p3m.a[i]       = 1.0/p3m.ai[i];
    p3m.cao_cut[i] = 0.5*p3m.a[i]*p3m.cao;
  }
}

int P3M_sanity_checks_boxl() {
  char *errtxt;
  int i;
  for(i=0;i<3;i++) {
    /* check k-space cutoff */
    if(p3m.cao_cut[i] >= 0.5*box_l[i]) {
      errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      sprintf(errtxt,"{P3M_init: k-space cutoff %f is larger than half of box dimension %f} ",p3m.cao_cut[i],box_l[i]);
      return 1;
    }
    if(p3m.cao_cut[i] >= local_box_l[i]) {
      errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      sprintf(errtxt,"{P3M_init: k-space cutoff %f is larger than local box dimension %f} ",p3m.cao_cut[i],local_box_l[i]);
      return 1;
    }
  }
  return 0;
}

int P3M_sanity_checks()
{
  char *errtxt;

  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    sprintf(errtxt, "{P3M requires periodicity 1 1 1} ");
    return 1;
  }
  
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    errtxt = runtime_error(128);
    sprintf(errtxt, "{P3M at present requires the domain decomposition cell system} ");
    return 1;
  }
  
  if( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
    errtxt = runtime_error(128);
    sprintf(errtxt,"{P3M requires a cubic box} ");
    return 1;
  }

  if( (p3m.mesh[0] != p3m.mesh[1]) || (p3m.mesh[1] != p3m.mesh[2]) ) {
    errtxt = runtime_error(128);
    sprintf(errtxt, "{P3M requires a cubic mesh} ");
    return 1;
  }

  if (P3M_sanity_checks_boxl()) return 1;

  return 0;
}

void   P3M_init()
{
  int n;

  if(coulomb.bjerrum == 0.0) {       
    p3m.r_cut    = 0.0;
    p3m.r_cut_iL = 0.0;
    if(this_node==0) 
      P3M_TRACE(fprintf(stderr,"0: P3M_init: Bjerrum length is zero.\n");
		fprintf(stderr,"   Electrostatics switched off!\n"));
  }
  else {  
    if (P3M_sanity_checks()) return;

    if( p3m.mesh[0] == 0 || p3m.cao == 0 ) {
      P3M_TRACE(fprintf(stderr,"%d: P3M_init: Not enough data: return\n",this_node));
      return;
    }

    P3M_TRACE(fprintf(stderr,"%d: P3M_init: \n",this_node));

    P3M_TRACE(fprintf(stderr,"%d: mesh=%d, cao=%d, mesh_off=(%f,%f,%f)\n",this_node,p3m.mesh[0],p3m.cao,p3m.mesh_off[0],p3m.mesh_off[1],p3m.mesh_off[2]));

    /* initializes the (inverse) mesh constant p3m.a (p3m.ai) and the cutoff for charge assignment p3m.cao_cut */
    P3M_init_a_ai_cao_cut();

    /* initialize ca fields to size CA_INCREMENT: ca_frac and ca_fmp */
    ca_num = 0;
    if(ca_num < CA_INCREMENT) {
      ca_num = 0;
      realloc_ca_fields(CA_INCREMENT);
    }
 
    calc_local_ca_mesh();
    calc_send_mesh();
    P3M_TRACE(p3m_print_local_mesh(lm));
    /* DEBUG */
    for(n=0;n<n_nodes;n++) {
      /* MPI_Barrier(MPI_COMM_WORLD); */
      if(n==this_node) P3M_TRACE(p3m_print_send_mesh(sm));
    }
    if(sm.max != send_recv_grid_size) {
      send_recv_grid_size=sm.max;
      send_grid = (double *) realloc(send_grid, sizeof(double)*sm.max);
      recv_grid = (double *) realloc(recv_grid, sizeof(double)*sm.max);
    }

    interpolate_charge_assignment_function();
    /* position offset for calc. of first meshpoint */
    pos_shift = (double)((p3m.cao-1)/2) - (p3m.cao%2)/2.0;
    P3M_TRACE(fprintf(stderr,"%d: pos_shift = %f\n",this_node,pos_shift)); 
 
    /* FFT */
    P3M_TRACE(fprintf(stderr,"%d: rs_mesh ADR=%p\n",this_node,rs_mesh));
    ca_mesh_size = fft_init(&rs_mesh,lm.dim,lm.margin,&ks_pnum);
    /* rs_mesh = (double *) realloc(rs_mesh, ca_mesh_size*sizeof(double)); */
    ks_mesh = (double *) realloc(ks_mesh, ca_mesh_size*sizeof(double));

    P3M_TRACE(fprintf(stderr,"%d: rs_mesh ADR=%p\n",this_node,rs_mesh));
 
    /* k-space part: */
    calc_differential_operator();
    calc_influence_function();


    P3M_count_charged_particles();

    P3M_TRACE(fprintf(stderr,"%d: p3m initialized\n",this_node));
  }
}

double P3M_calc_kspace_forces(int force_flag, int energy_flag)
{
  Cell *cell;
  Particle *p;
  int i,c,n, np,d,d_rs,i0,i1,i2,ind,j[3];
  double q;
  /* position of a particle in local mesh units */
  double pos[3];
  /* index of first assignment mesh point; argument for interpolated caf */
  int first[3],arg[3];
  /* index, index jumps for rs_mesh array */
  int q_ind, q_m_off, q_s_off;
  /* (2*p3m.inter) + 1 */
  int inter2;
  /* tmp variables */
  double tmp0,tmp1;
  /* charged particle counter, charge fraction counter */
  int cp_cnt=0, cf_cnt=0;
  /* Prefactor for force */
  double force_prefac;
  /* k space energy */
  double k_space_energy=0.0, node_k_space_energy=0.0;

  P3M_TRACE(fprintf(stderr,"%d: p3m_perform: \n",this_node));

  /* prepare local FFT mesh */
  for(n=0; n<lm.size; n++) rs_mesh[n] = 0.0;
  q_m_off = (lm.dim[2] - p3m.cao);
  q_s_off = lm.dim[2] * (lm.dim[1] - p3m.cao);
  
  /* === charge assignment === */
  force_prefac = coulomb.prefactor / (double)(p3m.mesh[0]*p3m.mesh[1]*p3m.mesh[2]);
  inter2 = (p3m.inter*2)+1;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if( (q=p[i].p.q) != 0.0 ) {

	/* particle position in mesh coordinates */
	for(d=0;d<3;d++) {
	  pos[d]   = ((p[i].r.p[d]-lm.ld_pos[d])*p3m.ai[d]) - pos_shift;
	  first[d] = (int) pos[d];
	  ca_fmp[(3*cp_cnt)+d] = first[d];
	  arg[d]   = (int) ((pos[d]-first[d])*inter2);

#ifdef ADDITIONAL_CHECKS
	  if( pos[d]<0.0 ) {
	    fprintf(stderr,"%d: rs_mesh underflow! (P%d at %f)\n",
		    this_node,p[i].p.identity,p[i].r.p[d]);
	    fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		    this_node,my_left[d],my_right[d]);	    
	  }
	  if( (first[d]+p3m.cao) > lm.dim[d] ) {
	    fprintf(stderr,"%d: rs_mesh overflow! dir=%d (P_id=%d at %f) first=%d\n",
		    this_node,d,p[i].p.identity,p[i].r.p[d],first[d]);
	    fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		    this_node,my_left[d],my_right[d]);
	  }    
#endif

	}

	/* charge assignment */
	q_ind = first[2] + lm.dim[2]*(first[1] + (lm.dim[1]*first[0]));
	for(i0=0; i0<p3m.cao; i0++) {
	  tmp0 = q * int_caf[i0][arg[0]];
	  for(i1=0; i1<p3m.cao; i1++) {
	    tmp1 = tmp0 * int_caf[i1][arg[1]];
	    for(i2=0; i2<p3m.cao; i2++) {
	      ca_frac[cf_cnt] = tmp1 * int_caf[i2][arg[2]];
	      /*
		P3M_TRACE(fprintf(stderr,"%d Frac %d = %f add tp [%d,%d,%d] (%d)\n",
		this_node,cf_cnt,ca_frac[cf_cnt],
		first[0]+i0,first[1]+i1,first[2]+i2,q_ind));
	      */
	      rs_mesh[q_ind++] += ca_frac[cf_cnt++];
	    }
	    q_ind += q_m_off;
	  }
	  q_ind += q_s_off;
	}
	cp_cnt++;
	if( (cp_cnt+1)>ca_num ) realloc_ca_fields(cp_cnt+1);
      }
    }
  }
  if( (cp_cnt+CA_INCREMENT)<ca_num ) realloc_ca_fields(cp_cnt+CA_INCREMENT);

  /* Gather information for FFT grid inside the nodes domain (inner local mesh) */
  gather_fft_grid();


  /* === Perform forward 3D FFT (Charge Assignment Mesh) === */
  fft_perform_forw(rs_mesh);

  /* === K Space Calculations === */
  P3M_TRACE(fprintf(stderr,"%d: p3m_perform: k-Space\n",this_node));

  /* === K Space Energy Calculation  === */
  if(energy_flag) {
    ind = 0;
    for(i=0; i<fft_plan[3].new_size; i++) {
      node_k_space_energy += g[i] * ( SQR(rs_mesh[ind]) + SQR(rs_mesh[ind+1]) );
      ind += 2;
    }
    node_k_space_energy *= force_prefac * box_l[0] / (4.0*PI) ;
 
    MPI_Reduce(&node_k_space_energy, &k_space_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(this_node==0) {
      /* self energy correction */
      k_space_energy -= coulomb.prefactor*(p3m_sum_q2*p3m.alpha_L*box_l_i[0] * wupii);
      /* net charge correction */
      k_space_energy -= coulomb.prefactor*(p3m_square_sum_q*PI / (2.0*box_l[0]*SQR(p3m.alpha_L)));
    }
  }

  /* === K Space Force Calculation  === */
  if(force_flag) {
    /* Force preparation */
    ind = 0;
    for(i=0; i<fft_plan[3].new_size; i++) {
      ks_mesh[ind] = g[i] * rs_mesh[ind]; ind++;
      ks_mesh[ind] = g[i] * rs_mesh[ind]; ind++;
    }
    
    /* === 3 Fold backward 3D FFT (Force Component Meshs) === */
    
    /* Force component loop */
    for(d=0;d<3;d++) {  
      /* direction in k space: */
      d_rs = (d+ks_pnum)%3;
      /* srqt(-1)*k differentiation */
      ind=0;
      for(j[0]=0; j[0]<fft_plan[3].new_mesh[0]; j[0]++) {
	for(j[1]=0; j[1]<fft_plan[3].new_mesh[1]; j[1]++) {
	  for(j[2]=0; j[2]<fft_plan[3].new_mesh[2]; j[2]++) {
	    /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */ 
	    rs_mesh[ind] = -(ks_mesh[ind+1]*d_op[ j[d]+fft_plan[3].start[d] ]); ind++;
	    rs_mesh[ind] =   ks_mesh[ind-1]*d_op[ j[d]+fft_plan[3].start[d] ];  ind++;
	  }
	}
      }
      /* Back FFT force componenet mesh */
      fft_perform_back(rs_mesh);
      /* redistribute force componenet mesh */
      spread_force_grid();
      /* Assign force component from mesh to particle */
      cp_cnt=0; cf_cnt=0;
      for (c = 0; c < local_cells.n; c++) {
	cell = local_cells.cell[c];
	p  = cell->part;
	np = cell->n;
	for(i=0; i<np; i++) { 
	  if( (q=p[i].p.q) != 0.0 ) {
#ifdef ADDITIONAL_CHECKS
	    double db_fsum=0.0;
#endif
	    
	    q_ind = ca_fmp[(3*cp_cnt)+2] + lm.dim[2]*(ca_fmp[(3*cp_cnt)+1] + (lm.dim[1]*ca_fmp[(3*cp_cnt)+0]));
	    for(i0=0; i0<p3m.cao; i0++) {
	      for(i1=0; i1<p3m.cao; i1++) {
		for(i2=0; i2<p3m.cao; i2++) {
#ifdef ADDITIONAL_CHECKS
		  db_fsum += force_prefac*ca_frac[cf_cnt]*rs_mesh[q_ind];
#endif
#ifdef NPT
		  if(integ_switch == INTEG_METHOD_NPT_ISO) {
		    // nptiso.p_vir[0] -= force_prefac*ca_frac[cf_cnt]*rs_mesh[q_ind++]; 
		  }
#endif
		  p[i].f.f[d_rs] -= force_prefac*ca_frac[cf_cnt]*rs_mesh[q_ind++]; 
		  cf_cnt++;
		}
		q_ind += q_m_off;
	      }
	      q_ind += q_s_off;
	    }
	    cp_cnt++;
#ifdef ADDITIONAL_CHECKS
	    /*	    if(fabs(db_fsum)> 1.0) fprintf(stderr,"%d: Part %d: k-space-force = %e\n",this_node,p[i].p.identity,db_fsum); */
#endif
	    ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: P3M  f = (%.3e,%.3e,%.3e) in dir %d add %.5f\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],d_rs,-db_fsum));
	  }
	}
      }
    }
  }

  if (p3m.epsilon != P3M_EPSILON_METALLIC)
    k_space_energy -= calc_dipole_term(force_flag, energy_flag);

  return k_space_energy;
}

void   P3M_exit()
{
  int i;
  /* free memory */
  free(ca_frac);
  free(ca_fmp);
  free(send_grid);
  free(recv_grid);
  free(rs_mesh);
  free(ks_mesh); 
  free(rs_mesh); 
  for(i=0; i<p3m.cao; i++) free(int_caf[i]);
}

/************************************************************/

double calc_dipole_term(int force_flag, int energy_flag)
{
  int np, c, i, j;
  Particle *part;
  double pref = coulomb.prefactor*4*M_PI*box_l_i[0]*box_l_i[1]*box_l_i[2]/(2*p3m.epsilon + 1);
  double lcl_dm[3], gbl_dm[3];
  double en;

  for (j = 0; j < 3; j++)
    lcl_dm[j] = 0;

  for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;
    for (i = 0; i < np; i++) {
      for (j = 0; j < 3; j++)
	/* dipole moment with unfolded coordinates */
	lcl_dm[j] += part[i].p.q*(part[i].r.p[j] + part[i].l.i[j]*box_l[j]);
    }
  }
  
  MPI_Allreduce(lcl_dm, gbl_dm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (energy_flag)
    en = 0.5*pref*(SQR(gbl_dm[0]) + SQR(gbl_dm[1]) + SQR(gbl_dm[2]));
  else
    en = 0;
  if (force_flag) {
    for (j = 0; j < 3; j++)
      gbl_dm[j] *= pref;
    for (c = 0; c < local_cells.n; c++) {
      np   = local_cells.cell[c]->n;
      part = local_cells.cell[c]->part;
      for (i = 0; i < np; i++)
	for (j = 0; j < 3; j++)
	  part[i].f.f[j] -= gbl_dm[j]*part[i].p.q;
    }
  }
  return en;
}

/************************************************************/

void calc_local_ca_mesh() {
  int i;
  int ind[3];

  /* inner left down grid point (global index) */
  for(i=0;i<3;i++) lm.in_ld[i] = (int)ceil(my_left[i]*p3m.ai[i]-p3m.mesh_off[i]);
  /* inner up right grid point (global index) */
  for(i=0;i<3;i++) lm.in_ur[i] = (int)floor(my_right[i]*p3m.ai[i]-p3m.mesh_off[i]);
  
  /* correct roundof errors at boundary */
  for(i=0;i<3;i++) {
    if((my_right[i]*p3m.ai[i]-p3m.mesh_off[i])-lm.in_ur[i]<ROUND_ERROR_PREC) lm.in_ur[i]--;
    if(1.0+(my_left[i]*p3m.ai[i]-p3m.mesh_off[i])-lm.in_ld[i]<ROUND_ERROR_PREC) lm.in_ld[i]--;
  }
  /* inner grid dimensions */
  for(i=0;i<3;i++) lm.inner[i] = lm.in_ur[i] - lm.in_ld[i] + 1;
  /* index of left down grid point in global mesh */
  for(i=0;i<3;i++) 
    lm.ld_ind[i]=(int)ceil((my_left[i]-p3m.cao_cut[i]-skin)*p3m.ai[i]-p3m.mesh_off[i]);
  /* spacial position of left down mesh point */
  calc_lm_ld_pos();
  /* left down margin */
  for(i=0;i<3;i++) lm.margin[i*2] = lm.in_ld[i]-lm.ld_ind[i];
  /* up right grid point */
  for(i=0;i<3;i++) ind[i]=(int)floor((my_right[i]+p3m.cao_cut[i]+skin)*p3m.ai[i]-p3m.mesh_off[i]);
  /* correct roundof errors at up right boundary */
  for(i=0;i<3;i++)
    if(((my_right[i]+p3m.cao_cut[i]+skin)*p3m.ai[i]-p3m.mesh_off[i])-ind[i]==0) ind[i]--;
  /* up right margin */
  for(i=0;i<3;i++) lm.margin[(i*2)+1] = ind[i] - lm.in_ur[i];

  /* grid dimension */
  lm.size=1; 
  for(i=0;i<3;i++) {lm.dim[i] = ind[i] - lm.ld_ind[i] + 1; lm.size*=lm.dim[i];}
  /* reduce inner grid indices from global to local */
  for(i=0;i<3;i++) lm.in_ld[i] = lm.margin[i*2];
  for(i=0;i<3;i++) lm.in_ur[i] = lm.margin[i*2]+lm.inner[i];
}

void p3m_print_local_mesh(local_mesh l) 
{
  fprintf(stderr,"%d: local_mesh: dim=(%d,%d,%d), size=%d\n",this_node,
	  l.dim[0],l.dim[1],l.dim[2],l.size);
  fprintf(stderr,"%d:    ld_ind=(%d,%d,%d), ld_pos=(%f,%f,%f)\n",this_node,
	  l.ld_ind[0],l.ld_ind[1],l.ld_ind[2],
	  l.ld_pos[0],l.ld_pos[1],l.ld_pos[2]);
  fprintf(stderr,"%d:    inner=(%d,%d,%d) [(%d,%d,%d)-(%d,%d,%d)]\n",this_node,
	  l.inner[0],l.inner[1],l.inner[2],
	  l.in_ld[0],l.in_ld[1],l.in_ld[2],
	  l.in_ur[0],l.in_ur[1],l.in_ur[2]);
  fprintf(stderr,"%d:    margin = (%d,%d, %d,%d, %d,%d)\n",this_node,
	  l.margin[0],l.margin[1],l.margin[2],l.margin[3],l.margin[4],l.margin[5]);
  fprintf(stderr,"%d:    r_margin=(%d,%d, %d,%d, %d,%d)\n",this_node,
	  l.r_margin[0],l.r_margin[1],l.r_margin[2],l.r_margin[3],l.r_margin[4],l.r_margin[5]);
}

void calc_send_mesh()
{
  int i,j,evenodd;
  int done[3]={0,0,0};
  MPI_Status status;
  /* send grids */
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      /* left */
      sm.s_ld[i*2][j] = 0 + done[j]*lm.margin[j*2];
      if(j==i) sm.s_ur[i*2][j] = lm.margin[j*2]; 
      else     sm.s_ur[i*2][j] = lm.dim[j]-done[j]*lm.margin[(j*2)+1];
      /* right */
      if(j==i) sm.s_ld[(i*2)+1][j] = lm.in_ur[j];
      else     sm.s_ld[(i*2)+1][j] = 0 + done[j]*lm.margin[j*2];
      sm.s_ur[(i*2)+1][j] = lm.dim[j] - done[j]*lm.margin[(j*2)+1];
    }   
    done[i]=1;
  }
  sm.max=0;
  for(i=0;i<6;i++) {
    sm.s_size[i] = 1;
    for(j=0;j<3;j++) {
      sm.s_dim[i][j] = sm.s_ur[i][j]-sm.s_ld[i][j];
      sm.s_size[i] *= sm.s_dim[i][j];
    }
    if(sm.s_size[i]>sm.max) sm.max=sm.s_size[i];
  }
  /* communication */
  for(i=0;i<6;i++) {
    if(i%2==0) j = i+1;
    else       j = i-1;
    if(node_neighbors[i] != this_node) {
      /* two step communication: first all even positions than all odd */
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[i/2]+evenodd)%2==0)
	  MPI_Send(&(lm.margin[i]), 1, MPI_INT, 
		   node_neighbors[i],REQ_P3M_INIT,MPI_COMM_WORLD);
	else
	  MPI_Recv(&(lm.r_margin[j]), 1, MPI_INT,
		   node_neighbors[j],REQ_P3M_INIT,MPI_COMM_WORLD,&status);    
      }
    }
    else {
      lm.r_margin[j] = lm.margin[i];
    }
  }
  /* recv grids */
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) {
      if(j==i) {
	sm.r_ld[ i*2   ][j] = sm.s_ld[ i*2   ][j] + lm.margin[2*j];
	sm.r_ur[ i*2   ][j] = sm.s_ur[ i*2   ][j] + lm.r_margin[2*j];
	sm.r_ld[(i*2)+1][j] = sm.s_ld[(i*2)+1][j] - lm.r_margin[(2*j)+1];
	sm.r_ur[(i*2)+1][j] = sm.s_ur[(i*2)+1][j] - lm.margin[(2*j)+1];
      }
      else {
	sm.r_ld[ i*2   ][j] = sm.s_ld[ i*2   ][j];
	sm.r_ur[ i*2   ][j] = sm.s_ur[ i*2   ][j];
	sm.r_ld[(i*2)+1][j] = sm.s_ld[(i*2)+1][j];
	sm.r_ur[(i*2)+1][j] = sm.s_ur[(i*2)+1][j];
      }
    }
  for(i=0;i<6;i++) {
    sm.r_size[i] = 1;
    for(j=0;j<3;j++) {
      sm.r_dim[i][j] = sm.r_ur[i][j]-sm.r_ld[i][j];
      sm.r_size[i] *= sm.r_dim[i][j];
    }
    if(sm.r_size[i]>sm.max) sm.max=sm.r_size[i];
  }
}

void p3m_print_send_mesh(send_mesh sm) 
{
  int i;
  fprintf(stderr,"%d: send_mesh: max=%d\n",this_node,sm.max);
  for(i=0;i<6;i++) {
    fprintf(stderr,"%d:  dir=%d: s_dim (%d,%d,%d)  s_ld (%d,%d,%d) s_ur (%d,%d,%d) s_size=%d\n",this_node,i,sm.s_dim[i][0],sm.s_dim[i][1],sm.s_dim[i][2],sm.s_ld[i][0],sm.s_ld[i][1],sm.s_ld[i][2],sm.s_ur[i][0],sm.s_ur[i][1],sm.s_ur[i][2],sm.s_size[i]);
    fprintf(stderr,"%d:         r_dim (%d,%d,%d)  r_ld (%d,%d,%d) r_ur (%d,%d,%d) r_size=%d\n",this_node,sm.r_dim[i][0],sm.r_dim[i][1],sm.r_dim[i][2],sm.r_ld[i][0],sm.r_ld[i][1],sm.r_ld[i][2],sm.r_ur[i][0],sm.r_ur[i][1],sm.r_ur[i][2],sm.r_size[i]);
  }
}

void gather_fft_grid()
{
  int s_dir,r_dir,evenodd;
  MPI_Status status;
  double *tmp_ptr;

  P3M_TRACE(fprintf(stderr,"%d: gather_fft_grid:\n",this_node));

  /* direction loop */
  for(s_dir=0; s_dir<6; s_dir++) {
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */ 
    if(sm.s_size[s_dir]>0) 
      pack_block(rs_mesh, send_grid, sm.s_ld[s_dir], sm.s_dim[s_dir], lm.dim, 1);
    /* communication */
    if(node_neighbors[s_dir] != this_node) {
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[s_dir/2]+evenodd)%2==0) {
	  if(sm.s_size[s_dir]>0) 
	    MPI_Send(send_grid, sm.s_size[s_dir], MPI_DOUBLE, 
		     node_neighbors[s_dir], REQ_P3M_GATHER, MPI_COMM_WORLD);
	}
	else {
	  if(sm.r_size[r_dir]>0) 
	    MPI_Recv(recv_grid, sm.r_size[r_dir], MPI_DOUBLE, 
		     node_neighbors[r_dir], REQ_P3M_GATHER, MPI_COMM_WORLD, &status); 	    
	}
      }
    }
    else {
      tmp_ptr = recv_grid;
      recv_grid = send_grid;
      send_grid = tmp_ptr;
    }
    /* add recv block */
    if(sm.r_size[r_dir]>0) {
      add_block(recv_grid, rs_mesh, sm.r_ld[r_dir], sm.r_dim[r_dir], lm.dim); 
    }
  }
}

void spread_force_grid()
{
  int s_dir,r_dir,evenodd;
  MPI_Status status;
  double *tmp_ptr;
  P3M_TRACE(fprintf(stderr,"%d: spread_force_grid:\n",this_node));

  /* direction loop */
  for(s_dir=5; s_dir>=0; s_dir--) {
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */ 
    if(sm.s_size[s_dir]>0) 
      pack_block(rs_mesh, send_grid, sm.r_ld[r_dir], sm.r_dim[r_dir], lm.dim, 1);
    /* communication */
    if(node_neighbors[r_dir] != this_node) {
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[r_dir/2]+evenodd)%2==0) {
	  if(sm.r_size[r_dir]>0) 
	    MPI_Send(send_grid, sm.r_size[r_dir], MPI_DOUBLE, 
		     node_neighbors[r_dir], REQ_P3M_SPREAD, MPI_COMM_WORLD);
   	}
	else {
	  if(sm.s_size[s_dir]>0) 
	    MPI_Recv(recv_grid, sm.s_size[s_dir], MPI_DOUBLE, 
		     node_neighbors[s_dir], REQ_P3M_SPREAD, MPI_COMM_WORLD, &status); 	    
	}
      }
    }
    else {
      tmp_ptr = recv_grid;
      recv_grid = send_grid;
      send_grid = tmp_ptr;
    }
    /* un pack recv block */
    if(sm.s_size[s_dir]>0) {
      unpack_block(recv_grid, rs_mesh, sm.s_ld[s_dir], sm.s_dim[s_dir], lm.dim, 1); 
    }
  }
}

void realloc_ca_fields(int newsize)
{
  int incr = 0;
  if( newsize > ca_num ) incr = (newsize - ca_num)/CA_INCREMENT +1;
  else if( newsize < ca_num ) incr = (newsize - ca_num)/CA_INCREMENT +1;
  incr *= CA_INCREMENT;
  if(incr != 0) {
    P3M_TRACE(fprintf(stderr,"%d: realloc_ca_fields: old_size=%d -> new_size=%d\n",this_node,ca_num,ca_num+incr));
    ca_num += incr;
    if(ca_num<CA_INCREMENT) ca_num = CA_INCREMENT;
    ca_frac = (double *)realloc(ca_frac, p3m.cao*p3m.cao*p3m.cao*ca_num*sizeof(double));
    ca_fmp  = (int *)realloc(ca_fmp, 3*ca_num*sizeof(int));
  }  
}

void add_block(double *in, double *out, int start[3], int size[3], int dim[3])
{
  /* fast,mid and slow changing indices */
  int f,m,s;
  /* linear index of in grid, linear index of out grid */
  int li_in=0,li_out=0;
  /* offsets for indizes in output grid */
  int m_out_offset,s_out_offset;

  li_out = start[2] + ( dim[2]*( start[1] + (dim[1]*start[0]) ) );
  m_out_offset  = dim[2] - size[2];
  s_out_offset  = (dim[2] * (dim[1] - size[1]));

  for(s=0 ;s<size[0]; s++) {
    for(m=0; m<size[1]; m++) {
      for(f=0; f<size[2]; f++) {
	out[li_out++] += in[li_in++];
      }
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

void interpolate_charge_assignment_function()
{
  /* REMARK: This function is taken unchanged from polysim9 by M. Deserno. */
  double dInterpol=(double)p3m.inter, x;
  long   i;

  P3M_TRACE(fprintf(stderr,"%d - interpolating (%d) the order-%d charge assignment function\n",
		    this_node,p3m.inter,p3m.cao));

  for(i=0;i<p3m.cao;i++) 
    int_caf[i] = (double *) realloc(int_caf[i], sizeof(double)*(2*p3m.inter+1));

  switch (p3m.cao) {
  case 1 : { 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = 1.0;
    }
  } break;
  case 2 : { 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = 0.5-x;
      int_caf[1][i+p3m.inter] = 0.5+x;
    }
  } break;
  case 3 : { 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = 0.5*SQR(0.5 - x);
      int_caf[1][i+p3m.inter] = 0.75 - SQR(x);
      int_caf[2][i+p3m.inter] = 0.5*SQR(0.5 + x);
    }
  } break;
  case 4 :{ 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
      int_caf[1][i+p3m.inter] = (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
      int_caf[2][i+p3m.inter] = (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
      int_caf[3][i+p3m.inter] = ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    }
  } break;
  case 5 : {
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
      int_caf[1][i+p3m.inter] = ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
      int_caf[2][i+p3m.inter] = (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
      int_caf[3][i+p3m.inter] = ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
      int_caf[4][i+p3m.inter] = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    }
  } break;
  case 6 : {
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
      int_caf[1][i+p3m.inter] = (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
      int_caf[2][i+p3m.inter] = (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
      int_caf[3][i+p3m.inter] = (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
      int_caf[4][i+p3m.inter] = (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
      int_caf[5][i+p3m.inter] = (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    }
  } break;
  case 7 : {
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
      int_caf[1][i+p3m.inter] = (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
      int_caf[2][i+p3m.inter] = (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
      int_caf[3][i+p3m.inter] = ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
      int_caf[4][i+p3m.inter] = (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
      int_caf[5][i+p3m.inter] = (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
      int_caf[6][i+p3m.inter] = (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    }
  } break;
  default :{
    fprintf(stderr,"%d: Charge assignment order %d unknown.\n",this_node,p3m.cao);
  }
  }
}

void calc_meshift(void)
{
  int i;
  double dmesh;

  dmesh = (double)p3m.mesh[0];
  meshift = (double *) realloc(meshift, p3m.mesh[0]*sizeof(double));

  for (i=0; i<p3m.mesh[0]; i++) meshift[i] = i - dround(i/dmesh)*dmesh; 
}

void calc_differential_operator()
{
  int i;
  double dmesh;

  dmesh = (double)p3m.mesh[0];
  d_op = (double *) realloc(d_op, p3m.mesh[0]*sizeof(double));

  for (i=0; i<p3m.mesh[0]; i++) 
    d_op[i] = (double)i - dround((double)i/dmesh)*dmesh;

  d_op[p3m.mesh[0]/2] = 0;
}

void calc_influence_function()
{
  int i,n[3],ind;
  int end[3];
  int size=1;
  double fak1,fak2;
  double nominator[3]={0.0,0.0,0.0},denominator=0.0;

  calc_meshift();

  for(i=0;i<3;i++) {
    size *= fft_plan[3].new_mesh[i];
    end[i] = fft_plan[3].start[i] + fft_plan[3].new_mesh[i];
  }
  g = (double *) realloc(g, size*sizeof(double));

  fak1  = p3m.mesh[0]*p3m.mesh[0]*p3m.mesh[0]*2.0/(box_l[0]*box_l[0]);

  for(n[0]=fft_plan[3].start[0]; n[0]<end[0]; n[0]++) 
    for(n[1]=fft_plan[3].start[1]; n[1]<end[1]; n[1]++) 
      for(n[2]=fft_plan[3].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-fft_plan[3].start[2]) + fft_plan[3].new_mesh[2] * ((n[1]-fft_plan[3].start[1]) + (fft_plan[3].new_mesh[1]*(n[0]-fft_plan[3].start[0])));
	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	  g[ind] = 0.0;
	else if( (n[0]%(p3m.mesh[0]/2)==0) && 
		 (n[1]%(p3m.mesh[0]/2)==0) && 
		 (n[2]%(p3m.mesh[0]/2)==0) )
	  g[ind] = 0.0;
	else {
	  denominator = perform_aliasing_sums(n,nominator);
	  fak2 =  d_op[n[0]]*nominator[0] + d_op[n[1]]*nominator[1] + d_op[n[2]]*nominator[2];  
	  fak2 /= ( ( SQR(d_op[n[0]])+SQR(d_op[n[1]])+SQR(d_op[n[2]]) ) * SQR(denominator) );
	  g[ind] = fak1*fak2;
	}
      }
}

MDINLINE double perform_aliasing_sums(int n[3], double nominator[3])
{
  int i;
  double denominator=0.0;
  /* lots of temporary variables... */
  double sx,sy,sz,f1,f2,f3,mx,my,mz,nmx,nmy,nmz,nm2,expo;
  double limit = 30;

  for(i=0;i<3;i++) nominator[i]=0.0;
  f1 = 1.0/(double)p3m.mesh[0];
  f2 = SQR(PI/(p3m.alpha_L));

  for(mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = meshift[n[0]] + p3m.mesh[0]*mx;
    sx  = pow(sinc(f1*nmx),2.0*p3m.cao);
    for(my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = meshift[n[1]] + p3m.mesh[0]*my;
      sy  = sx*pow(sinc(f1*nmy),2.0*p3m.cao);
      for(mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
	nmz = meshift[n[2]] + p3m.mesh[0]*mz;
	sz  = sy*pow(sinc(f1*nmz),2.0*p3m.cao);
	
	nm2          =  SQR(nmx)+SQR(nmy)+SQR(nmz);
	expo         =  f2*nm2;
	f3           =  (expo<limit) ? sz*exp(-expo)/nm2 : 0.0;

	nominator[0] += f3*nmx; 
	nominator[1] += f3*nmy; 
	nominator[2] += f3*nmz; 
	denominator  += sz;
      }
    }
  }
  return denominator;
}

/************************************************
 * Functions for P3M Parameter tuning
 * This tuning is based on P3M_tune by M. Deserno
 ************************************************/

#define P3M_TUNE_MAX_CUTS 50

int P3M_tune_parameters(Tcl_Interp *interp)
{
  int i,ind, try=0, best_try=0, n_cuts;
  double r_cut_iL, r_cut_iL_min  , r_cut_iL_max, r_cut_iL_best=0, cuts[P3M_TUNE_MAX_CUTS], cut_start;
  int    mesh    , mesh_min      , mesh_max    , mesh_best=0;
  int    cao     , cao_min       , cao_max     , cao_best=0;
  double alpha_L , alpha_L_best=0, accuracy    , accuracy_best=0;
  double mesh_size, k_cut;
  double rs_err, rs_err_best=0, ks_err, ks_err_best=0;
  double int_time=0, min_time=1e20, int_num;
  char b1[TCL_DOUBLE_SPACE + 12],b2[TCL_DOUBLE_SPACE + 12],b3[TCL_DOUBLE_SPACE + 12];
 
  P3M_TRACE(fprintf(stderr,"%d: P3M_tune_parameters\n",this_node));
  
  /* preparation */
  mpi_bcast_event(P3M_COUNT_CHARGES);

  /* calculate r_cut_iL tune range */
  if(p3m.r_cut_iL == 0.0) { 
    n_cuts = P3M_TUNE_MAX_CUTS;
    for(i=0;i<n_cuts;i++) {
      if(min_local_box_l == min_box_l) 
	cuts[i] = min_local_box_l/(i+2.0)-(skin);
      else 
	cuts[i] = min_local_box_l/(i+1.0)-(skin);
      cuts[i]*=box_l_i[0];
      if( cuts[i] <= 0.0 ) {
	n_cuts = i; 
	break;
      } 
    }
    r_cut_iL_max = cuts[0];
    r_cut_iL_min = cuts[n_cuts-1];
  }
  else { 
    n_cuts = 1;
    r_cut_iL_min = r_cut_iL_max = p3m.r_cut_iL; 
    cuts[0] = p3m.r_cut_iL;
  }
  /* calculate mesh tune range */
  if(p3m.mesh[0] == 0 ) {
    double expo;
    expo = log(pow((double)p3m_sum_qpart,(1.0/3.0)))/log(2.0);
    mesh_min = (int)(pow(2.0,(double)((int)expo))+0.1);
    mesh_max = mesh_min*4;
    if(mesh_min < 8) { mesh_min = 8; mesh_max = 16; }
  }
  else { mesh_min = mesh_max = p3m.mesh[0]; }
  /* calculate cao tune range */
  if(p3m.cao == 0) { cao_min = 1; cao_max = 7; }
  else             { cao_min = cao_max = p3m.cao; }

  /* Print Status */
  sprintf(b1,"%.5e",p3m.accuracy);
  Tcl_AppendResult(interp, "P3M tune parameters: Accuracy goal = ",b1,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, box_l[0], b1);
  sprintf(b2,"%d",p3m_sum_qpart);
  Tcl_PrintDouble(interp, p3m_sum_q2, b3);
  Tcl_AppendResult(interp, "System: box_l = ",b1,", # charged part = ",b2," Sum[q_i^2] = ",b3,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, r_cut_iL_min, b1);  Tcl_PrintDouble(interp, r_cut_iL_max, b2);
  Tcl_AppendResult(interp, "Range for p3m.r_cut_iL: [",b1,"-",b2,"]","\n", (char *) NULL);
  sprintf(b1,"%d",mesh_min);  sprintf(b2,"%d",mesh_max);
  Tcl_AppendResult(interp, "Range for p3m.mesh:     [",b1,"-",b2,"]","\n", (char *) NULL);
  sprintf(b1,"%d",cao_min);  sprintf(b2,"%d",cao_max);
  Tcl_AppendResult(interp, "Range for p3m.cao:      [",b1,"-",b2,"]","\n\n", (char *) NULL);
  Tcl_AppendResult(interp, "set mesh cao r_cut_iL     alpha_L      err          ks_err     rs_err     time [ms]\n", (char *) NULL);

  /* Tuning Loops */
  for(mesh = mesh_min; mesh <= mesh_max; mesh*=2) { /* mesh loop */
    cut_start = box_l[0] * box_l_i[0];
    if(mesh <= 32 || p3m_sum_qpart > 2000) int_num=5; else int_num=1;
    for(cao = cao_min; cao <= cao_max; cao++) {     /* cao loop */
      mesh_size = box_l[0]/(double)mesh;
      k_cut =  mesh_size*cao/2.0;
      if(cao < mesh && k_cut < dmin(min_box_l,min_local_box_l)-skin) {
	ind=0;
	for(i=0;i< n_cuts -1;i++) {
	  if(cut_start <= cuts[i]) ind=i+1;
	}
	while (ind < n_cuts) {           /* r_cut_iL loop */
	  r_cut_iL = cuts[ind];
	  /* calc maximal real space error for setting */
	  rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,0);
	  if(sqrt(2.0)*rs_err > p3m.accuracy) {
	    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
	    alpha_L = sqrt(log(sqrt(2.0)*rs_err/p3m.accuracy)) / r_cut_iL;
	    /* calculate real space and k space error for this alpha_L */
	    rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,alpha_L);
	    ks_err = P3M_k_space_error(box_l[0],coulomb.prefactor,mesh,cao,p3m_sum_qpart,p3m_sum_q2,alpha_L);
	    accuracy = sqrt(SQR(rs_err)+SQR(ks_err));
	    /* check if this matches the accuracy goal */
	    if(accuracy <= p3m.accuracy) {
	      cut_start = cuts[ind];
	      /* broadcast p3m parameters for test run */
	      p3m.r_cut_iL = r_cut_iL;
	      p3m.mesh[0]  = p3m.mesh[1] = p3m.mesh[2] = mesh;
	      p3m.cao      = cao;
	      p3m.alpha_L  = alpha_L;
	      P3M_scaleby_box_l();
	      /* initialize p3m structures */
	      mpi_bcast_coulomb_params();
	      /* perform force calculation test */
	      int_time = time_force_calc(int_num);
	      try++;
	      P3M_TRACE(fprintf(stderr,"%d ",try));
	      /* print result */
	      sprintf(b1,"%-3d",try); sprintf(b2,"%-4d",mesh); sprintf(b3,"%-3d",cao);
	      Tcl_AppendResult(interp, b1," ", b2," ", b3," ", (char *) NULL);
	      sprintf(b1,"%.5e",r_cut_iL); sprintf(b2,"%.5e",alpha_L); sprintf(b3,"%.5e",accuracy);
	      Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
	      sprintf(b1,"%.3e",rs_err); sprintf(b2,"%.3e",ks_err); sprintf(b3,"%-8d",(int)int_time);
	      Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);
	      if(int_time <= min_time  && r_cut_iL > 0) {
		min_time      = int_time;
		r_cut_iL_best = r_cut_iL;
		mesh_best     = mesh;
		cao_best      = cao;
		alpha_L_best  = alpha_L;
		accuracy_best = sqrt(SQR(rs_err)+SQR(ks_err));
		rs_err_best   = rs_err;
		ks_err_best   = ks_err;
		best_try      = try;
	      }
	    }
	  }
	  ind++;
	}
      }
    }
  }
  P3M_TRACE(fprintf(stderr,"\n"));
  if(try==0) {
    Tcl_AppendResult(interp, "\nFailed to tune P3M parameters to required accuracy ! \n", (char *) NULL);
    return (TCL_ERROR);
  }

  /* set tuned p3m parameters */
  p3m.r_cut_iL = r_cut_iL_best;
  p3m.mesh[0]  = p3m.mesh[1] = p3m.mesh[2] = mesh_best;
  p3m.cao      = cao_best;
  p3m.alpha_L  = alpha_L_best;
  p3m.accuracy = accuracy_best;
  P3M_scaleby_box_l();
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  /* Tell the user about the outcome */
  sprintf(b1,"%d",try);
  Tcl_AppendResult(interp, "\nTune results of ",b1," trials:\n", (char *) NULL);
  sprintf(b1,"%-3d",best_try); sprintf(b2,"%-4d",mesh_best); sprintf(b3,"%-3d",cao_best);
  Tcl_AppendResult(interp, b1," ", b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",r_cut_iL_best); sprintf(b2,"%.5e",alpha_L_best); sprintf(b3,"%.5e",accuracy_best);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);
  sprintf(b1,"%.3e",rs_err_best); sprintf(b2,"%.3e",ks_err_best); sprintf(b3,"%-8d",(int)min_time);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"\n", (char *) NULL);
  sprintf(b1,"%g",coulomb.bjerrum); sprintf(b2,"%g",p3m.r_cut); sprintf(b3,"%d",mesh_best); 
  Tcl_AppendResult(interp, "=> inter coulomb ", b1, " p3m ", b2, " ", b3, (char *) NULL);
  sprintf(b1,"%d",cao_best); sprintf(b2,"%g",p3m.alpha); sprintf(b3,"%g",accuracy_best);
  Tcl_AppendResult(interp, " ", b1," ", b2," ", b3," \n", (char *) NULL);

  return (TCL_OK);
}

void P3M_count_charged_particles()
{  
  Cell *cell;
  Particle *part;
  int i,c,np;
  double node_sums[3], tot_sums[3];

  for(i=0;i<3;i++) { node_sums[i]=0.0; tot_sums[i]=0.0;}

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
      if( part[i].p.q != 0.0 ) {
	node_sums[0] += 1.0;
	node_sums[1] += SQR(part[i].p.q);
	node_sums[2] += part[i].p.q;
      }
    }
  }
  
  MPI_Reduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(this_node==0) {
    p3m_sum_qpart    = (int)(tot_sums[0]+0.1);
    p3m_sum_q2       = tot_sums[1];
    p3m_square_sum_q = tot_sums[2];
  }

}


double P3M_real_space_error(double box_size, double prefac, double r_cut_iL, 
			    int n_c_part, double sum_q2, double alpha_L)
{
  return (2.0*prefac*sum_q2*exp(-SQR(r_cut_iL*alpha_L))) / (sqrt(n_c_part*r_cut_iL)*box_size*box_size);
}

double P3M_k_space_error(double box_size, double prefac, int mesh, 
			 int cao, int n_c_part, double sum_q2, double alpha_L)
{
  int  nx, ny, nz;
  double he_q = 0.0;
  double alias1, alias2, n2, cs;

  /*   fprintf(stderr,"KS ERR: box %.8f, prefac %.8f, mesh %d, cao %d\n",
	  box_size,prefac,mesh,cao); 
  fprintf(stderr,"        n_c_part %d, sum_q2 %.8f,  alpha_L %.8f\n", 
  n_c_part,sum_q2,alpha_L); */

  for (nx=-mesh/2; nx<mesh/2; nx++)
    for (ny=-mesh/2; ny<mesh/2; ny++)
      for (nz=-mesh/2; nz<mesh/2; nz++)
	if((nx!=0) || (ny!=0) || (nz!=0)) {
	  n2 = SQR(nx) + SQR(ny) + SQR(nz);
	  cs = analytic_cotangent_sum(nx,mesh,cao)*
 	       analytic_cotangent_sum(ny,mesh,cao)*
	       analytic_cotangent_sum(nz,mesh,cao);
	  P3M_tune_aliasing_sums(nx,ny,nz,mesh,cao,alpha_L,&alias1,&alias2);
	  he_q += (alias1  -  SQR(alias2/cs) / n2);
	  /* fprintf(stderr,"%d %d %d he_q = %.20f %.20f %.20f %.20f\n",nx,ny,nz,he_q,cs,alias1,alias2); */
	}
  return 2.0*prefac*sum_q2*sqrt(he_q/(double)n_c_part) / SQR(box_size);
}

double analytic_cotangent_sum(int n, int mesh, int cao)
{
  double c, res=0.0;
  c = SQR(cos(PI*(double)n/(double)mesh));

  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  default : {
    fprintf(stderr,"%d: INTERNAL_ERROR: The value %d for the interpolation order should not occur!\n",this_node, cao);
    errexit();
  }
  }
  
  return res;
}

void P3M_tune_aliasing_sums(int nx, int ny, int nz, 
			    int mesh, int cao, double alpha_L , 
			    double *alias1, double *alias2)
{

  int    mx,my,mz;
  double nmx,nmy,nmz;
  double fnmx,fnmy,fnmz;

  double ex,ex2,nm2,U2,factor1,factor2;

  factor1 = SQR(PI/(alpha_L));
  factor2 = 1.0/(double)mesh;

  *alias1 = *alias2 = 0.0;
  for (mx=-P3M_BRILLOUIN; mx<=P3M_BRILLOUIN; mx++) {
    fnmx = factor2 * (nmx = nx + mx*mesh);
    for (my=-P3M_BRILLOUIN; my<=P3M_BRILLOUIN; my++) {
      fnmy = factor2 * (nmy = ny + my*mesh);
      for (mz=-P3M_BRILLOUIN; mz<=P3M_BRILLOUIN; mz++) {
	fnmz = factor2 * (nmz = nz + mz*mesh);
	
	nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	ex2 = SQR( ex = exp(-factor1*nm2) );
	
	U2 = pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex2 / nm2;
	*alias2 += U2 * ex * (nx*nmx + ny*nmy + nz*nmz) / nm2;
	/* fprintf(stderr,"%d %d %d : %d %d %d alias %.20f %.20f\n",nx,ny,nz,mx,my,mz,*alias1,*alias2); */
      }
    }
  }
}


/************************************************
 * Debug functions printing p3m structures 
 ************************************************/

void p3m_print_p3m_struct(p3m_struct ps) {
  fprintf(stderr,"%d: p3m_struct: \n",this_node);
  fprintf(stderr,"   alpha_L=%f, r_cut_iL=%f \n",
	  ps.alpha_L,ps.r_cut_iL);
  fprintf(stderr,"   mesh=(%d,%d,%d), mesh_off=(%.4f,%.4f,%.4f)\n",
	  ps.mesh[0],ps.mesh[1],ps.mesh[2],
	  ps.mesh_off[0],ps.mesh_off[1],ps.mesh_off[2]);
  fprintf(stderr,"   cao=%d, inter=%d, epsilon=%f\n",
	  ps.cao,ps.inter,ps.epsilon);
  fprintf(stderr,"   cao_cut=(%f,%f,%f)\n",
	  ps.cao_cut[0],ps.cao_cut[1],ps.cao_cut[2]);
  fprintf(stderr,"   a=(%f,%f,%f), ai=(%f,%f,%f)\n",
	  ps.a[0],ps.a[1],ps.a[2],ps.ai[0],ps.ai[1],ps.ai[2]);
}

#endif

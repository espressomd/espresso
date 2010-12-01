/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file p3m-charges.c  P3M algorithm for long range coulomb interaction.
 *
 *  For more information about the p3m algorithm,
 *  see \ref p3m.h "p3m.h"
 *  see \ref p3m-charges.c  "p3m-charges.c"
 *  see \ref p3m-charges.h  "p3m-charges.h"
 *  see \ref p3m-dipoles.c  "p3m-dipoles.c"
 *  see \ref p3m-dipoles.h  "p3m-dipoles.h"
 *  see \ref p3m-assignment.c  "p3m-assignment.c"
*/

#ifdef ELECTROSTATICS

/* only include from within p3m.c */
#ifndef P3M_C_CURRENT
#error never compile this file file directly, it is part of p3m.c
#endif

/************************************************
 * DEFINES
 ************************************************/

/* MPI tags for the charge-charge p3m communications: */
/** Tag for communication in P3M_init() -> send_calc_mesh(). */
#define REQ_P3M_INIT   200
/** Tag for communication in gather_fft_grid(). */
#define REQ_P3M_GATHER 201
/** Tag for communication in spread_force_grid(). */
#define REQ_P3M_SPREAD 202


/***********************/

/** interpolation of the charge assignment function. */
  double *int_caf[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
/** position shift for calc. of first assignment mesh point. */
  double pos_shift;
/** help variable for calculation of aliasing sums */
  double *meshift = NULL;
/** Spatial differential operator in k-space. We use an i*k differentiation. */
  double *d_op = NULL;
/** Force optimised influence function (k-space) */
  double *g_force = NULL;
/** Energy optimised influence function (k-space) */
double *g_energy = NULL;
/** number of charged particles on the node. */
int ca_num=0;
/** Charge fractions for mesh assignment. */
double *ca_frac = NULL;
/** index of first mesh point for charge assignment. */
int *ca_fmp = NULL;
/** number of permutations in k_space */
int ks_pnum;


/** number of charged particles (only on master node). */
int p3m_sum_qpart=0;
/** Sum of square of charges (only on master node). */
double p3m_sum_q2 = 0.0;
/** square of sum of charges (only on master node). */
double p3m_square_sum_q = 0.0;

/** local mesh. */
local_mesh lm;

/** send/recv mesh sizes */
send_mesh  sm;


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


/** \name Private Functions */
/************************************************************/
/*@{*/


/** Calculates for charges the properties of the send/recv sub-meshes of the local FFT mesh. 
 *  In order to calculate the recv sub-meshes there is a communication of 
 *  the margins between neighbouring nodes. */ 
void calc_send_mesh();


/** Initializes the (inverse) mesh constant \ref p3m_struct::a (\ref p3m_struct::ai) 
    and the cutoff for charge assignment \ref p3m_struct::cao_cut, which has to be
    done by \ref P3M_init once and by \ref P3M_scaleby_box_l_charges whenever the \ref box_l changed.
*/
void P3M_init_a_ai_cao_cut(void);


/** Calculate the spacial position of the left down mesh point of the local mesh, to be
    stored in \ref local_mesh::ld_pos; function called by \ref calc_local_ca_mesh once
    and by \ref P3M_scaleby_box_l_charges whenever the \ref box_l changed. */
void calc_lm_ld_pos(void);


/** Calculates the dipole term */
double calc_dipole_term(int force_flag, int energy_flag);

/** Gather FFT grid.
 *  After the charge assignment Each node needs to gather the
 *  information for the FFT grid in his spatial domain.
 */
void gather_fft_grid(double* mesh);

/** Spread force grid.
 *  After the k-space calculations each node needs to get all force
 *  information to reassigne the forces from the grid to the
 *  particles.
 */
void spread_force_grid(double* mesh);

/** realloc charge assignment fields. */
void realloc_ca_fields(int newsize);

/** Initializes the (inverse) mesh constant \ref p3m_struct::a (\ref p3m_struct::ai) 
    and the cutoff for charge assignment \ref p3m_struct::cao_cut, which has to be
    done by \ref P3M_init once and by \ref P3M_scaleby_box_l_charges whenever the \ref box_l changed.
*/
void P3M_init_a_ai_cao_cut(void);

/** checks for correctness for charges in P3M of the cao_cut, necessary when the box length changes */
int P3M_sanity_checks_boxl(void);

/** Calculate the spacial position of the left down mesh point of the local mesh, to be
    stored in \ref local_mesh::ld_pos; function called by \ref calc_local_ca_mesh once
    and by \ref P3M_scaleby_box_l_charges whenever the \ref box_l changed. */
void calc_lm_ld_pos(void);

/** Calculates properties of the local FFT mesh for the 
    charge assignment process. */
void calc_local_ca_mesh(void);


/** Interpolates the P-th order charge assignment function from
 * Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
 * are also tabulated in Deserno/Holm. */
void interpolate_charge_assignment_function(void);

/** shifts the mesh points by mesh/2 */
void calc_meshift(void);

/** Calculates the Fourier transformed differential operator.  
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *           i.e. the prefactor i*2*PI/L is missing! */
void calc_differential_operator(void);

/** Calculates the optimal influence function of Hockney and Eastwood. 
 * (optimised for force calculations)
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft_plan[3].mesh and fft_plan[3].start).
 *
 *  See also: Hockney/Eastwood 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  Deserno/Holm. */
void calc_influence_function_force(void);

/** Calculates the influence function optimized for the energy and the
    self energy correction.  */
void calc_influence_function_energy(void);


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
MDINLINE double perform_aliasing_sums_force(int n[3], double nominator[3]);
MDINLINE double perform_aliasing_sums_energy(int n[3]);

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



/** aliasing sum used by \ref P3M_k_space_error. */
void P3M_tune_aliasing_sums(int nx, int ny, int nz, 
			    int mesh, double mesh_i, int cao, double alpha_L_i, 
			    double *alias1, double *alias2);

void p3m_set_tune_params(double r_cut, int mesh, int cao,
			 double alpha, double accuracy, int n_interpol)
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

  if (n_interpol != -1)
    p3m.inter = n_interpol;

  coulomb.prefactor = (temperature > 0) ? temperature*coulomb.bjerrum : coulomb.bjerrum;
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




int tclcommand_inter_coulomb_parse_p3m_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive)
{
  int mesh = -1, cao = -1, n_interpol = -1;
  double r_cut = -1, accuracy = -1;

  while(argc > 0) {
    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(r_cut) && r_cut >= -1)) {
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
	Tcl_AppendResult(interp, "accuracy expects a positive double",
			 (char *) NULL);
	return TCL_ERROR;
      }

    } else if (ARG0_IS_S("n_interpol")) {
      if (! (argc > 1 && ARG1_IS_I(n_interpol) && n_interpol >= 0)) {
	Tcl_AppendResult(interp, "n_interpol expects an nonnegative integer",
			 (char *) NULL);
	return TCL_ERROR;
      }
    }
    /* unknown parameter. Probably one of the optionals */
    else break;
    
    argc -= 2;
    argv += 2;
  }
  p3m_set_tune_params(r_cut, mesh, cao, -1.0, accuracy, n_interpol);

  /* check for optional parameters */
  if (argc > 0) {
    if (tclcommand_inter_coulomb_parse_p3m_opt_params(interp, argc, argv) == TCL_ERROR)
      return TCL_ERROR;
  }

  if (adaptive) {
    if(tclcommand_inter_coulomb_print_p3m_adaptive_tune_parameteres(interp) == TCL_ERROR) 
      return TCL_ERROR;
  }
  else {
    if(tclcommand_inter_coulomb_print_p3m_tune_parameteres(interp) == TCL_ERROR) 
      return TCL_ERROR;
  }

  return TCL_OK;
}




int tclcommand_inter_coulomb_parse_p3m(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha, accuracy = -1.0;
  int mesh, cao, i;

  if (coulomb.method != COULOMB_P3M && coulomb.method != COULOMB_ELC_P3M)
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
    Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> p3m tune | <r_cut> <mesh> <cao> [<alpha> [<accuracy>]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    Tcl_AppendResult(interp, "Node grid not suited for Coulomb P3M. Node grid must be sorted, largest first.", (char *) NULL);
    return TCL_ERROR;  
  }

  if (ARG0_IS_S("tune"))
    return tclcommand_inter_coulomb_parse_p3m_tune(interp, argc-1, argv+1, 0);

  if (ARG0_IS_S("tunev2"))
    return tclcommand_inter_coulomb_parse_p3m_tune(interp, argc-1, argv+1, 1);
      
  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc < 3 || argc > 5) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb <bjerrum> p3m <r_cut> <mesh> <cao> [<alpha> [<accuracy>]]",
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




int tclcommand_inter_coulomb_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv)
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
      else if (! ARG1_IS_D(d1)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 DOUBLE parameter or \"metallic\"",
	                 (char *) NULL);
	return TCL_ERROR;
      }
	
      if (p3m_set_eps(d1) == TCL_ERROR) {
        Tcl_AppendResult(interp, argv[0], " There is no error msg yet!",
                         (char *) NULL);
        return TCL_ERROR;
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

/************************************* method ********************************/
/*****************************************************************************/

void interpolate_charge_assignment_function()
{
  double dInterpol = 0.5 / (double)p3m.inter;
  int i;
  long j;

  if (p3m.inter == 0) return;

  P3M_TRACE(fprintf(stderr,"%d - interpolating (%d) the order-%d charge assignment function\n",
		    this_node,p3m.inter,p3m.cao));

  p3m.inter2 = 2*p3m.inter + 1;

  for (i=0; i < p3m.cao; i++) {
    /* allocate memory for interpolation array */
    int_caf[i] = (double *) realloc(int_caf[i], sizeof(double)*(2*p3m.inter+1));

    /* loop over all interpolation points */
    for (j=-p3m.inter; j<=p3m.inter; j++)
      int_caf[i][j+p3m.inter] = P3M_caf(i, j*dInterpol,p3m.cao);
  }
  
}

/* assign the charges */
void P3M_charge_assign()
{
  Cell *cell;
  Particle *p;
  int i,c,np;
  /* charged particle counter, charge fraction counter */
  int cp_cnt=0;
  /* prepare local FFT mesh */
  for(i=0; i<lm.size; i++) rs_mesh[i] = 0.0;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if( p[i].p.q != 0.0
	  ) {
	P3M_assign_charge(p[i].p.q, p[i].r.p, cp_cnt);
	cp_cnt++;
      }
    }
  }
  P3M_shrink_wrap_charge_grid(cp_cnt);
  
}

/* assign the forces obtained from k-space */
static void P3M_assign_forces(double force_prefac, int d_rs) 
{
  Cell *cell;
  Particle *p;
  int i,c,np,i0,i1,i2;
  double q;
#ifdef ONEPART_DEBUG
  double db_fsum=0 ; /* TODO: db_fsum was missing and code couldn't compile. Now it has the arbitrary value of 0, fix it. */ 
#endif
  /* charged particle counter, charge fraction counter */
  int cp_cnt=0, cf_cnt=0;
  /* index, index jumps for rs_mesh array */
  int q_ind;
  int q_m_off = (lm.dim[2] - p3m.cao);
  int q_s_off = lm.dim[2] * (lm.dim[1] - p3m.cao);

  cp_cnt=0; cf_cnt=0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i=0; i<np; i++) { 
      if( (q=p[i].p.q) != 0.0 ) {
	q_ind = ca_fmp[cp_cnt];
	for(i0=0; i0<p3m.cao; i0++) {
	  for(i1=0; i1<p3m.cao; i1++) {
	    for(i2=0; i2<p3m.cao; i2++) {
	      p[i].f.f[d_rs] -= force_prefac*ca_frac[cf_cnt]*rs_mesh[q_ind++]; 
	      cf_cnt++;
	    }
	    q_ind += q_m_off;
	  }
	  q_ind += q_s_off;
	}
	cp_cnt++;

	ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: P3M  f = (%.3e,%.3e,%.3e) in dir %d add %.5f\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],d_rs,-db_fsum));
      }
    }
  }
}








double P3M_calc_kspace_forces_for_charges(int force_flag, int energy_flag) 
{
  int i,d,d_rs,ind,j[3];
  /**************************************************************/
  /* Prefactor for force */
  double force_prefac;
  /* k space energy */
  double k_space_energy=0.0, node_k_space_energy=0.0;

  P3M_TRACE(fprintf(stderr,"%d: p3m_perform: \n",this_node));

  force_prefac = coulomb.prefactor / (double)(p3m.mesh[0]*p3m.mesh[1]*p3m.mesh[2]);

  /* Gather information for FFT grid inside the nodes domain (inner local mesh) */
  /* and Perform forward 3D FFT (Charge Assignment Mesh). */
  if (p3m_sum_q2 > 0) {
    gather_fft_grid(rs_mesh);
    fft_perform_forw(rs_mesh);
    }
//Note: after these calls, the grids are in the order yzx and not xyz anymore!!!

  /* === K Space Calculations === */
  P3M_TRACE(fprintf(stderr,"%d: p3m_perform: k-Space\n",this_node));

  /* === K Space Energy Calculation  === */
  if(energy_flag) {
/*********************
   Coulomb energy
**********************/
   if (p3m_sum_q2 > 0) {
    ind = 0;
    for(i=0; i<fft_plan[3].new_size; i++) {
      // Use the energy optimized influence function for energy!
      node_k_space_energy += g_energy[i] * ( SQR(rs_mesh[ind]) + SQR(rs_mesh[ind+1]) );
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

  } //if (energy_flag)

  /* === K Space Force Calculation  === */
  if(force_flag) {
 /***************************
    COULOMB FORCES (k-space)
****************************/ 
    if (p3m_sum_q2 > 0) {
    /* Force preparation */
    ind = 0;

    /* apply the influence function */
    for(i=0; i<fft_plan[3].new_size; i++) {
      ks_mesh[ind] = g_force[i] * rs_mesh[ind]; ind++;
      ks_mesh[ind] = g_force[i] * rs_mesh[ind]; ind++;
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
      /* Back FFT force component mesh */
      fft_perform_back(rs_mesh);
      /* redistribute force component mesh */
      spread_force_grid(rs_mesh);
      /* Assign force component from mesh to particle */
      P3M_assign_forces(force_prefac, d_rs);
    }
   }  // if(p3m_sum_q2>0)

  } // if(force_flag)

  if (p3m.epsilon != P3M_EPSILON_METALLIC) {
    k_space_energy += calc_dipole_term(force_flag, energy_flag);
  }

  return k_space_energy;
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

void gather_fft_grid(double* themesh)
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
      pack_block(themesh, send_grid, sm.s_ld[s_dir], sm.s_dim[s_dir], lm.dim, 1);
      
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
      add_block(recv_grid, themesh, sm.r_ld[r_dir], sm.r_dim[r_dir], lm.dim); 
    }
  }
}


void spread_force_grid(double* themesh)
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
      pack_block(themesh, send_grid, sm.r_ld[r_dir], sm.r_dim[r_dir], lm.dim, 1);
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
      unpack_block(recv_grid, themesh, sm.s_ld[s_dir], sm.s_dim[s_dir], lm.dim, 1); 
    }
  }
}


void realloc_ca_fields(int newsize)
{
  newsize = ((newsize + CA_INCREMENT - 1)/CA_INCREMENT)*CA_INCREMENT;
  if (newsize == ca_num) return;
  if (newsize < CA_INCREMENT) newsize = CA_INCREMENT;

  P3M_TRACE(fprintf(stderr,"%d: realloc_ca_fields: old_size=%d -> new_size=%d\n",this_node,ca_num,newsize));
  ca_num = newsize;
  ca_frac = (double *)realloc(ca_frac, p3m.cao3*ca_num*sizeof(double));
  ca_fmp  = (int *)realloc(ca_fmp, ca_num*sizeof(int));
    
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

void calc_influence_function_force()
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
  g_force = (double *) realloc(g_force, size*sizeof(double));

  fak1  = p3m.mesh[0]*p3m.mesh[0]*p3m.mesh[0]*2.0/(box_l[0]*box_l[0]);

  for(n[0]=fft_plan[3].start[0]; n[0]<end[0]; n[0]++) 
    for(n[1]=fft_plan[3].start[1]; n[1]<end[1]; n[1]++) 
      for(n[2]=fft_plan[3].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-fft_plan[3].start[2]) + fft_plan[3].new_mesh[2] * ((n[1]-fft_plan[3].start[1]) + (fft_plan[3].new_mesh[1]*(n[0]-fft_plan[3].start[0])));

	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	  g_force[ind] = 0.0;
	else if( (n[0]%(p3m.mesh[0]/2)==0) && 
		 (n[1]%(p3m.mesh[0]/2)==0) && 
		 (n[2]%(p3m.mesh[0]/2)==0) )
	  g_force[ind] = 0.0;
	else {
	  denominator = perform_aliasing_sums_force(n,nominator);
	  fak2 =  d_op[n[0]]*nominator[0] + d_op[n[1]]*nominator[1] + d_op[n[2]]*nominator[2];  
	  fak2 /= ( ( SQR(d_op[n[0]])+SQR(d_op[n[1]])+SQR(d_op[n[2]]) ) * SQR(denominator) );
	  g_force[ind] = fak1*fak2;
	}
      }
}

MDINLINE double perform_aliasing_sums_force(int n[3], double nominator[3])
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

void calc_influence_function_energy()
{
  int i,n[3],ind;
  int end[3];
  int size=1;
  double fak1;

  calc_meshift();

  for(i=0;i<3;i++) {
    size *= fft_plan[3].new_mesh[i];
    end[i] = fft_plan[3].start[i] + fft_plan[3].new_mesh[i];
  }
  g_energy = (double *) realloc(g_energy, size*sizeof(double));

  fak1  = p3m.mesh[0]*p3m.mesh[0]*p3m.mesh[0]*2.0/(box_l[0]*box_l[0]);

  for(n[0]=fft_plan[3].start[0]; n[0]<end[0]; n[0]++) 
    for(n[1]=fft_plan[3].start[1]; n[1]<end[1]; n[1]++) 
      for(n[2]=fft_plan[3].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-fft_plan[3].start[2]) + fft_plan[3].new_mesh[2] * ((n[1]-fft_plan[3].start[1]) + (fft_plan[3].new_mesh[1]*(n[0]-fft_plan[3].start[0])));
	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	  g_energy[ind] = 0.0;
	else if( (n[0]%(p3m.mesh[0]/2)==0) && 
		 (n[1]%(p3m.mesh[0]/2)==0) && 
		 (n[2]%(p3m.mesh[0]/2)==0) )
	  g_energy[ind] = 0.0;
	else {
	  g_energy[ind] = fak1*perform_aliasing_sums_energy(n);
	}
      }
}

MDINLINE double perform_aliasing_sums_energy(int n[3])
{
  double nominator=0.0,denominator=0.0;
  /* lots of temporary variables... */
  double sx,sy,sz,f1,f2,f3,mx,my,mz,nmx,nmy,nmz,nm2,expo;
  double limit = 30;

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

	nominator += f3; 
	denominator  += sz;
      }
    }
  }
  return nominator/SQR(denominator);
}


/************************************************
 * Functions for P3M Parameter tuning
 * This tuning is based on P3M_tune by M. Deserno
 ************************************************/

#define P3M_TUNE_MAX_CUTS 50
/* TODO: separate tcl from nontcl thingies */
int tclcommand_inter_coulomb_print_p3m_tune_parameteres(Tcl_Interp *interp)
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
 
  P3M_TRACE(fprintf(stderr,"%d: tclcommand_inter_coulomb_print_p3m_tune_parameteres\n",this_node));
  
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
	  //Beginning of JJCP modification 15/5/2006 
	    rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,0);
	  //Original line -> rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,0);

	  if(sqrt(2.0)*rs_err > p3m.accuracy) {
	    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
	    alpha_L = sqrt(log(sqrt(2.0)*rs_err/p3m.accuracy)) / r_cut_iL;
	    /* calculate real space and k space error for this alpha_L */
	    alpha_L = sqrt(log(sqrt(2.0)*rs_err/p3m.accuracy)) / r_cut_iL;
	    rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,alpha_L);
	    ks_err = P3M_k_space_error(box_l[0],coulomb.prefactor,mesh,cao,p3m_sum_qpart,p3m_sum_q2,alpha_L);
	  //End of JJCP modification 15/5/2006 
	  // Original line -> rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,alpha_L);
	  // Original line ->  ks_err = P3M_k_space_error(box_l[0],coulomb.prefactor,mesh,cao,p3m_sum_qpart,p3m_sum_q2,alpha_L);
	  //End of JJCP modification 15/5/2006 
	    accuracy = sqrt(SQR(rs_err)+SQR(ks_err));
	    /* check if this matches the accuracy goal */
	    if(accuracy <= p3m.accuracy) {
	      cut_start = cuts[ind];
	      /* broadcast p3m parameters for test run */
	      p3m.r_cut_iL = r_cut_iL;
	      p3m.mesh[0]  = p3m.mesh[1] = p3m.mesh[2] = mesh;
	      p3m.cao      = cao;
	      p3m.alpha_L  = alpha_L;
	      P3M_scaleby_box_l_charges();
	      /* initialize p3m structures */
	      mpi_bcast_coulomb_params();
	      /* perform force calculation test */
	      int_time = time_force_calc(int_num);
	      if (int_time == -1)
		return TCL_ERROR;
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
  P3M_scaleby_box_l_charges();
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

/** get the minimal error for this combination of parameters. In fact, the real space error is tuned such that it
    contributes half of the total error, and then the Fourier space error is calculated. Returns the error and the
    optimal alpha, or 0 if this combination does not work at all */
static double get_accuracy(int mesh, int cao, double r_cut_iL, double *_alpha_L, double *_rs_err, double *_ks_err)
{
  double rs_err, ks_err;
  double alpha_L;
  P3M_TRACE(fprintf(stderr, "get_accuracy: mesh %d, cao %d, r_cut %f ", mesh, cao, r_cut_iL));

  /* calc maximal real space error for setting */
  /* calc maximal real space error for setting */
  //Beginning of JJCP modification 15/5/2006 
  rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,0);
  
  //Original line ->  rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,0);
  
    if(M_SQRT2*rs_err > p3m.accuracy) {
     /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
        alpha_L = sqrt(log(M_SQRT2*rs_err/p3m.accuracy)) / r_cut_iL;
    
    }

  else
    /* even alpha=0 is ok, however, we cannot choose it since it kills the k-space error formula.
       Anyways, this very likely NOT the optimal solution */
    alpha_L = 0.1;

  *_alpha_L = alpha_L;
  /* calculate real space and k space error for this alpha_L */
    rs_err = P3M_real_space_error(box_l[0],coulomb.prefactor,r_cut_iL,p3m_sum_qpart,p3m_sum_q2,alpha_L);
    ks_err = P3M_k_space_error(box_l[0],coulomb.prefactor,mesh,cao,p3m_sum_qpart,p3m_sum_q2,alpha_L);
  
 //End of JJCP modification 15/5/2006 
  *_rs_err = rs_err;
  *_ks_err = ks_err;
  P3M_TRACE(fprintf(stderr, "resulting: %f -> %f %f\n", alpha_L, rs_err, ks_err));
  return sqrt(SQR(rs_err)+SQR(ks_err));
}

/** get the optimal alpha and the corresponding computation time for fixed mesh, cao, r_cut and alpha */
static double p3m_mcr_time(int mesh, int cao, double r_cut_iL, double alpha_L)
{
  /* rounded up 2000/n_charges timing force evaluations */
  int int_num = (1999 + p3m_sum_qpart)/p3m_sum_qpart;
  

  /* broadcast p3m parameters for test run */
  p3m.r_cut_iL = r_cut_iL;
  p3m.mesh[0]  = p3m.mesh[1] = p3m.mesh[2] = mesh;
  p3m.cao      = cao;
  p3m.alpha_L  = alpha_L;
  P3M_scaleby_box_l_charges();
  /* initialize p3m structures */
  mpi_bcast_coulomb_params();
  /* perform force calculation test */
  return time_force_calc(int_num);    
}

/** get the optimal alpha and the corresponding computation time for fixed mesh, cao. The r_cut is determined via
    a simple bisection. Returns -1 if the force evaluation does not work, -2 if there is no valid r_cut, and -3 if
    the charge assigment order is to large for this grid */
/* TODO: separate tcl from nontcl thingies */
static double p3m_mc_time(int mesh, int cao,
			  double r_cut_iL_min, double r_cut_iL_max, double *_r_cut_iL,
			  double *_alpha_L, double *_accuracy)
{
  double int_time;
  double r_cut_iL;
  double rs_err, ks_err, mesh_size, k_cut;
  int i, n_cells;

  /* initial checks. */
  mesh_size = box_l[0]/(double)mesh;
  k_cut =  mesh_size*cao/2.0;
  P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_mc_time: mesh=%d, cao=%d, rmin=%f, rmax=%f\n",
		    mesh, cao, r_cut_iL_min, r_cut_iL_max));
  if(cao >= mesh || k_cut >= dmin(min_box_l,min_local_box_l) - skin) {
    return -P3M_TUNE_CAOTOLARGE;
  }

  /* Either low and high boundary are equal (for fixed cut), or the low border is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary fails, there is no possible r_cut */
  if ((*_accuracy = get_accuracy(mesh, cao, r_cut_iL_max, _alpha_L, &rs_err, &ks_err)) > p3m.accuracy) {
    /* print result */
    return -P3M_TUNE_NOCUTOFF;
  }

  for (;;) {
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_mc_time: interval [%f,%f]\n", r_cut_iL_min, r_cut_iL_max));
    r_cut_iL = 0.5*(r_cut_iL_min + r_cut_iL_max);

    if (r_cut_iL_max - r_cut_iL_min < P3M_RCUT_PREC)
      break;

    /* bisection */
    if (get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err) > p3m.accuracy)
      r_cut_iL_min = r_cut_iL;
    else
      r_cut_iL_max = r_cut_iL;
  }
  /* final result is always the upper interval boundary, since only there
     we know that the desired minimal accuracy is obtained */
  *_r_cut_iL = r_cut_iL = r_cut_iL_max;

  /* check whether we are running P3M+ELC, and whether we leave a reasonable gap space */
  if (coulomb.method == COULOMB_ELC_P3M && elc_params.gap_size <= 1.1*r_cut_iL*box_l[0]) {
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_mc_time: mesh %d cao %d r_cut %f reject r_cut %f > gap %f\n", mesh, cao, r_cut_iL,
		      2*r_cut_iL*box_l[0], elc_params.gap_size));
    return -P3M_TUNE_ELCTEST;
  }

  /* check whether this radius is too large, so that we would use less cells than allowed */
  n_cells = 1;
  for (i = 0; i < 3; i++)
    n_cells *= (int)(floor(local_box_l[i]/(r_cut_iL*box_l[0] + skin)));
  if (n_cells < min_num_cells) {
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_mc_time: mesh %d cao %d r_cut %f reject n_cells %d\n", mesh, cao, r_cut_iL, n_cells));

    return -P3M_TUNE_CUTOFF_TOO_LARGE;
  }

  int_time = p3m_mcr_time(mesh, cao, r_cut_iL, *_alpha_L);
  if (int_time == -1) {
    return -P3M_TUNE_FAIL;
  }

  *_accuracy = get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_mc_time: mesh %d cao %d r_cut %f time %f\n", mesh, cao, r_cut_iL, int_time));

  return int_time;
}

/** get the optimal alpha and the corresponding computation time for fixed mesh. *cao
    should contain an initial guess, which is then adapted by stepping up and down. Returns the time
    upon completion, -1 if the force evaluation does not work, and -2 if the accuracy cannot be met */
static double tclcommand_inter_coulomb_print_p3m_m_time(Tcl_Interp *interp, int mesh,
			 int cao_min, int cao_max, int *_cao,
			 double r_cut_iL_min, double r_cut_iL_max, double *_r_cut_iL,
			 double *_alpha_L, double *_accuracy)
{
  double best_time = -1, tmp_time, tmp_r_cut_iL, tmp_alpha_L=0.0, tmp_accuracy=0.0;
  /* in which direction improvement is possible. Initially, we dont know it yet. */
  int final_dir = 0;
  int cao = *_cao;

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: mesh=%d, cao_min=%d, cao_max=%d, rmin=%f, rmax=%f\n",
		    mesh, cao_min, cao_max, r_cut_iL_min, r_cut_iL_max));
  /* the initial step sets a timing mark. If there is no valid r_cut, we can only try
     to increase cao to increase the obtainable precision of the far formula. */
  do {
    tmp_time = p3m_mc_time(mesh, cao,  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out if the force evaluation is not working */
    if (tmp_time == -1) return -1;
    /* cao is too large for this grid, but still the accuracy cannot be achieved, give up */
    if (tmp_time == -3) {
      P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: no possible cao found\n"));
      return -2;
    }
    /* we have a valid time, start optimising from there */
    if (tmp_time >= 0) {
      best_time  = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L  = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao      = cao;
      break;
    }
    /* the required accuracy could not be obtained, try higher caos. Therefore optimisation can only be
       obtained with even higher caos, but not lower ones */
    P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: doesn't give precision, step up\n"));
    cao++;
    final_dir = 1;
  }
  while (cao <= cao_max);
  /* with this mesh, the required accuracy cannot be obtained. */
  if (cao > cao_max) return -2;

  /* at the boundaries, only the opposite direction can be used for optimisation */
  if (cao == cao_min)      final_dir = 1;
  else if (cao == cao_max) final_dir = -1;

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: final constraints dir %d\n", final_dir));

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time =
        p3m_mc_time(mesh, cao + final_dir,  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
      /* bail out on errors, as usual */
      if (tmp_time == -1) return -P3M_TUNE_FAIL;
      /* in this direction, we cannot optimise, since we get into precision trouble */
      if (tmp_time < 0) continue;

      if (tmp_time < best_time) {
	best_time  = tmp_time;
	*_r_cut_iL = tmp_r_cut_iL;
	*_alpha_L  = tmp_alpha_L;
	*_accuracy = tmp_accuracy;
	*_cao      = cao + final_dir;
      }
    }
    /* choose the direction which was optimal, if any of the two */
    if      (dir_times[0] == best_time) { final_dir = -1; }
    else if (dir_times[2] == best_time) { final_dir = 1; }
    else {
      /* no improvement in either direction, however if one is only marginally worse, we can still try*/
      /* down is possible and not much worse, while up is either illegal or even worse */
      if ((dir_times[0] >= 0 && dir_times[0] < best_time + P3M_TIME_GRAN) &&
	  (dir_times[2] < 0 || dir_times[2] > dir_times[0]))
	final_dir = -1;
      /* same for up */
      else if ((dir_times[2] >= 0 && dir_times[2] < best_time + P3M_TIME_GRAN) &&
	       (dir_times[0] < 0 || dir_times[0] > dir_times[2]))
	final_dir = 1;
      else {
	/* really no chance for optimisation */
	P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: mesh=%d final cao=%d time=%f\n",mesh, cao, best_time));
	return best_time;
      }
    }
    /* we already checked the initial cao and its neighbor */
    cao += 2*final_dir;
  }
  else {
    /* here some constraint is active, and we only checked the initial cao itself */
    cao += final_dir;
  }

  P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: optimise in direction %d\n", final_dir));

  /* move cao into the optimisation direction until we do not gain anymore. */
  for (; cao >= cao_min && cao <= cao_max; cao += final_dir) {
    tmp_time = tclcommand_inter_coulomb_print_p3m_mc_time(interp, mesh, cao,  r_cut_iL_min, r_cut_iL_max,
			   &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out on errors, as usual */
    if (tmp_time == -1) return -1;
    /* if we cannot meet the precision anymore, give up */
    if (tmp_time < 0) break;

    if (tmp_time < best_time) {
      best_time  = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L  = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao      = cao;
    }
    /* no hope of further optimisation */
    else if (tmp_time > best_time + P3M_TIME_GRAN)
      break;
  }
  P3M_TRACE(fprintf(stderr, "tclcommand_inter_coulomb_print_p3m_m_time: mesh=%d final cao=%d r_cut=%f time=%f\n",mesh, *_cao, *_r_cut_iL, best_time));
  return best_time;
}

int p3m_adaptive_tune() {
    int    mesh_max,                   mesh     = -1, tmp_mesh;
  double r_cut_iL_min, r_cut_iL_max, r_cut_iL = -1, tmp_r_cut_iL=0.0;
  int    cao_min, cao_max,           cao      = -1, tmp_cao;

  double                             alpha_L  = -1, tmp_alpha_L=0.0;
  double                             accuracy = -1, tmp_accuracy=0.0;
  double                            time_best=1e20, tmp_time;

  P3M_TRACE(fprintf(stderr,"%d: p3m_adaptive_tune\n",this_node));
  
  /* preparation */
  mpi_bcast_event(P3M_COUNT_CHARGES);

  /* parameter ranges */
  if (p3m.mesh[0] == 0 ) {
    double expo;
    expo = log(pow((double)p3m_sum_qpart,(1.0/3.0)))/log(2.0);

    tmp_mesh = (int)(pow(2.0,(double)((int)expo))+0.1);
    /* this limits the tried meshes if the accuracy cannot
       be obtained with smaller meshes, but normally not all these
       meshes have to be tested */
    mesh_max = tmp_mesh * 256;
    /* avoid using more than 1 GB of FFT arrays (per default, see config.h) */
    if (mesh_max > P3M_MAX_MESH)
      mesh_max = P3M_MAX_MESH;
  }
  else {
    tmp_mesh = mesh_max = p3m.mesh[0];
  }

  if(p3m.r_cut_iL == 0.0) {
    r_cut_iL_min = 0;
    r_cut_iL_max = min_local_box_l/2 - skin;
    r_cut_iL_min *= box_l_i[0];
    r_cut_iL_max *= box_l_i[0];
  }
  else {
    r_cut_iL_min = r_cut_iL_max = p3m.r_cut_iL;
  }

  if(p3m.cao == 0) {
    cao_min = 1;
    cao_max = 7;
    cao = 3;
  }
  else {
    cao_min = cao_max = cao = p3m.cao;
  }

  /* mesh loop */
  for (;tmp_mesh <= mesh_max; tmp_mesh *= 2) {
    tmp_cao = cao;
    
#warning FIX ME!
    
/*    tmp_time = tclcommand_inter_coulomb_print_p3m_m_time(interp, tmp_mesh,
			  cao_min, cao_max, &tmp_cao,
			  r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL,
			  &tmp_alpha_L, &tmp_accuracy); */
    /* some error occured during the tuning force evaluation */
    if (tmp_time == -1) return P3M_TUNE_FAIL;
    /* this mesh does not work at all */
    if (tmp_time < 0) continue;

    /* the optimum r_cut for this mesh is the upper limit for higher meshes,
       everything else is slower */
    r_cut_iL_max = tmp_r_cut_iL;

    /* new optimum */
    if (tmp_time < time_best) {
      time_best = tmp_time;
      mesh      = tmp_mesh;
      cao       = tmp_cao;
      r_cut_iL  = tmp_r_cut_iL;
      alpha_L   = tmp_alpha_L;
      accuracy  = tmp_accuracy;
    }
    /* no hope of further optimisation */
    else if (tmp_time > time_best + P3M_TIME_GRAN)
      break;
  }
  
  P3M_TRACE(fprintf(stderr,"finshed tuning\n"));
  if(time_best == 1e20) {
    return -1;
  }

  /* set tuned p3m parameters */
  p3m.r_cut_iL = r_cut_iL;
  p3m.mesh[0]  = p3m.mesh[1] = p3m.mesh[2] = mesh;
  p3m.cao      = cao;
  p3m.alpha_L  = alpha_L;
  p3m.accuracy = accuracy;
  P3M_scaleby_box_l_charges();
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  return 0;
}
  


int tclcommand_inter_coulomb_print_p3m_adaptive_tune_parameteres(Tcl_Interp *interp)
{
  char
    b1[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 12],
    b2[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 12],
    b3[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE + 17];
 
  P3M_TRACE(fprintf(stderr,"%d: tclcommand_inter_coulomb_print_p3m_adaptive_tune_parameteres\n",this_node));

  if (skin == -1) {
    Tcl_AppendResult(interp, "p3m cannot be tuned, since the skin is not yet set", (char *) NULL);
    return TCL_ERROR;
  }

  /* Print Status */
  sprintf(b1,"%.5e",p3m.accuracy);
  Tcl_AppendResult(interp, "P3M tune parameters: Accuracy goal = ",b1,"\n", (char *) NULL);
  Tcl_PrintDouble(interp, box_l[0], b1);

  sprintf(b2,"%d",p3m_sum_qpart);
  Tcl_PrintDouble(interp, p3m_sum_q2, b3);
  Tcl_AppendResult(interp, "System: box_l = ",b1,", # charged part = ",b2," Sum[q_i^2] = ",b3,"\n", (char *) NULL);

  if (p3m_sum_qpart == 0) {
    Tcl_AppendResult(interp, "no charged particles in the system, cannot tune P3M", (char *) NULL);
    return (TCL_ERROR);
  }
  
  if(!p3m_adaptive_tune()) {  
    Tcl_AppendResult(interp, "failed to tune P3M parameters to required accuracy", (char *) NULL);
    return (TCL_ERROR);
  }
  
  /* Tell the user about the outcome */
  Tcl_AppendResult(interp, "\nresulting parameters:\n", (char *) NULL);
  sprintf(b2,"%-4d",p3m.mesh[0]); sprintf(b3,"%-3d",p3m.cao);
  Tcl_AppendResult(interp, b2," ", b3," ", (char *) NULL);
  sprintf(b1,"%.5e",p3m.r_cut_iL); sprintf(b2,"%.5e",p3m.alpha_L); sprintf(b3,"%.5e",p3m.accuracy);
  Tcl_AppendResult(interp, b1,"  ", b2,"  ", b3,"  ", (char *) NULL);

  return (TCL_OK);  
}

void P3M_count_charged_particles()
{  
  Cell *cell;
  Particle *part;
  int i,c,np;
  double node_sums[3], tot_sums[3];
  for(i=0;i<3;i++)
    { node_sums[i]=0.0; tot_sums[i]=0.0;}

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
  
  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  p3m_sum_qpart    = (int)(tot_sums[0]+0.1);
  p3m_sum_q2       = tot_sums[1];
  p3m_square_sum_q = SQR(tot_sums[2]);
}


double P3M_real_space_error(double box_size, double prefac, double r_cut_iL, 
			    int n_c_part, double sum_q2, double alpha_L)
{
// Modification done just to be sure about n_c_part,  JJCP 15/5/06
//  return (2.0*prefac*sum_q2*exp(-SQR(r_cut_iL*alpha_L))) / (sqrt(n_c_part*r_cut_iL)*box_size*box_size);
  return (2.0*prefac*sum_q2*exp(-SQR(r_cut_iL*alpha_L))) / (sqrt((double)n_c_part*r_cut_iL)*box_size*box_size);
}

double P3M_k_space_error(double box_size, double prefac, int mesh, 
			 int cao, int n_c_part, double sum_q2, double alpha_L)
{
  int  nx, ny, nz;
  double he_q = 0.0, mesh_i = 1./mesh, alpha_L_i = 1./alpha_L;
  double alias1, alias2, n2, cs;

  for (nx=-mesh/2; nx<mesh/2; nx++)
    for (ny=-mesh/2; ny<mesh/2; ny++)
      for (nz=-mesh/2; nz<mesh/2; nz++)
	if((nx!=0) || (ny!=0) || (nz!=0)) {
	  n2 = SQR(nx) + SQR(ny) + SQR(nz);
	  cs = analytic_cotangent_sum(nx,mesh_i,cao)*
 	       analytic_cotangent_sum(ny,mesh_i,cao)*
	       analytic_cotangent_sum(nz,mesh_i,cao);
	  P3M_tune_aliasing_sums(nx,ny,nz,mesh,mesh_i,cao,alpha_L_i,&alias1,&alias2);
	  he_q += (alias1  -  SQR(alias2/cs) / n2);
	  /* fprintf(stderr,"%d %d %d he_q = %.20f %.20f %.20f %.20f\n",nx,ny,nz,he_q,cs,alias1,alias2); */
	}
  return 2.0*prefac*sum_q2*sqrt(he_q/(double)n_c_part) / SQR(box_size);
}


void P3M_tune_aliasing_sums(int nx, int ny, int nz, 
			    int mesh, double mesh_i, int cao, double alpha_L_i, 
			    double *alias1, double *alias2)
{

  int    mx,my,mz;
  double nmx,nmy,nmz;
  double fnmx,fnmy,fnmz;

  double ex,ex2,nm2,U2,factor1;

  factor1 = SQR(PI*alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (mx=-P3M_BRILLOUIN; mx<=P3M_BRILLOUIN; mx++) {
    fnmx = mesh_i * (nmx = nx + mx*mesh);
    for (my=-P3M_BRILLOUIN; my<=P3M_BRILLOUIN; my++) {
      fnmy = mesh_i * (nmy = ny + my*mesh);
      for (mz=-P3M_BRILLOUIN; mz<=P3M_BRILLOUIN; mz++) {
	fnmz = mesh_i * (nmz = nz + mz*mesh);
	
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




/************************************************************/

void calc_local_ca_mesh() {
  int i;
  int ind[3];
  /* total skin size */
  double full_skin[3];
  
  for(i=0;i<3;i++)
    full_skin[i]= p3m.cao_cut[i]+skin+p3m.additional_mesh[i];

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
    lm.ld_ind[i]=(int)ceil((my_left[i]-full_skin[i])*p3m.ai[i]-p3m.mesh_off[i]);
  /* spacial position of left down mesh point */
  calc_lm_ld_pos();
  /* left down margin */
  for(i=0;i<3;i++) lm.margin[i*2] = lm.in_ld[i]-lm.ld_ind[i];
  /* up right grid point */
  for(i=0;i<3;i++) ind[i]=(int)floor((my_right[i]+full_skin[i])*p3m.ai[i]-p3m.mesh_off[i]);
  /* correct roundof errors at up right boundary */
  for(i=0;i<3;i++)
    if(((my_right[i]+full_skin[i])*p3m.ai[i]-p3m.mesh_off[i])-ind[i]==0) ind[i]--;
  /* up right margin */
  for(i=0;i<3;i++) lm.margin[(i*2)+1] = ind[i] - lm.in_ur[i];

  /* grid dimension */
  lm.size=1; 
  for(i=0;i<3;i++) {lm.dim[i] = ind[i] - lm.ld_ind[i] + 1; lm.size*=lm.dim[i];}
  /* reduce inner grid indices from global to local */
  for(i=0;i<3;i++) lm.in_ld[i] = lm.margin[i*2];
  for(i=0;i<3;i++) lm.in_ur[i] = lm.margin[i*2]+lm.inner[i];

  lm.q_2_off  = lm.dim[2] - p3m.cao;
  lm.q_21_off = lm.dim[2] * (lm.dim[1] - p3m.cao);
 
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
  int i, ret = 0;
  for(i=0;i<3;i++) {
    /* check k-space cutoff */
    if(p3m.cao_cut[i] >= 0.5*box_l[i]) {
      errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      ERROR_SPRINTF(errtxt,"{039 P3M_init: k-space cutoff %g is larger than half of box dimension %g} ",p3m.cao_cut[i],box_l[i]);
      ret = 1;
    }
    if(p3m.cao_cut[i] >= local_box_l[i]) {
      errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      ERROR_SPRINTF(errtxt,"{040 P3M_init: k-space cutoff %g is larger than local box dimension %g} ",p3m.cao_cut[i],local_box_l[i]);
      ret = 1;
    }
  }
  
  
  return ret;
}




int P3M_sanity_checks()
{
  char *errtxt;
  int ret = 0;

  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{041 P3M requires periodicity 1 1 1} ");
    ret = 1;
  }
  
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{042 P3M at present requires the domain decomposition cell system} ");
    ret = 1;
  }
  
  if( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{043 P3M requires a cubic box} ");
    ret = 1;
  }

  if( (p3m.mesh[0] != p3m.mesh[1]) || (p3m.mesh[1] != p3m.mesh[2]) ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{044 P3M requires a cubic mesh} ");
    ret = 1;
  }
   

  if (P3M_sanity_checks_boxl()) ret = 1;

  if( p3m.mesh[0] == 0) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{045 P3M_init: mesh size is not yet set} ");
    ret = 1;
  }
  if( p3m.cao == 0) {
    errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
    ERROR_SPRINTF(errtxt,"{046 P3M_init: cao is not yet set} ");
    ret = 1;
  }
  if (skin == -1) {
    errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
    ERROR_SPRINTF(errtxt,"{047 P3M_init: skin is not yet set} ");
    ret = 1;
  }
    
  return ret;
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



/************************************************/

void P3M_scaleby_box_l_charges() {

  p3m.r_cut = p3m.r_cut_iL* box_l[0];
  p3m.alpha = p3m.alpha_L * box_l_i[0];
   P3M_init_a_ai_cao_cut();
  calc_lm_ld_pos();
  P3M_sanity_checks_boxl(); 
}

/************************************************/



void   P3M_init_charges() {
  int n;

  if(coulomb.bjerrum == 0.0) {       
    
      if(coulomb.bjerrum == 0.0) {
           p3m.r_cut    = 0.0;
           p3m.r_cut_iL = 0.0;
          if(this_node==0) 
             P3M_TRACE(fprintf(stderr,"0: P3M_init: Bjerrum length is zero.\n");
	   fprintf(stderr,"   Electrostatics switched off!\n"));
      }
      
  } else {  
    P3M_TRACE(fprintf(stderr,"%d: P3M_init_charges: \n",this_node));

    if (P3M_sanity_checks()) return;

    P3M_TRACE(fprintf(stderr,"%d: P3M_init_charges: starting\n",this_node));

    P3M_TRACE(fprintf(stderr,"%d: mesh=%d, cao=%d, mesh_off=(%f,%f,%f)\n",this_node,p3m.mesh[0],p3m.cao,p3m.mesh_off[0],p3m.mesh_off[1],p3m.mesh_off[2]));
    p3m.cao3 = p3m.cao*p3m.cao*p3m.cao;


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

    if (p3m.inter > 0)
      interpolate_charge_assignment_function();
  
    /* position offset for calc. of first meshpoint */
    pos_shift = (double)((p3m.cao-1)/2) - (p3m.cao%2)/2.0;
    P3M_TRACE(fprintf(stderr,"%d: pos_shift = %f\n",this_node,pos_shift)); 
 
    /* FFT */
    P3M_TRACE(fprintf(stderr,"%d: rs_mesh ADR=%p\n",this_node,rs_mesh));
 
    ca_mesh_size = fft_init(&rs_mesh,lm.dim,lm.margin,&ks_pnum);
    ks_mesh = (double *) realloc(ks_mesh, ca_mesh_size*sizeof(double));
    

    P3M_TRACE(fprintf(stderr,"%d: rs_mesh ADR=%p\n",this_node,rs_mesh));
  
    /* k-space part: */
    calc_differential_operator();
    
    
    calc_influence_function_force();
    calc_influence_function_energy();

    P3M_count_charged_particles();

    P3M_TRACE(fprintf(stderr,"%d: p3m-charges  initialized\n",this_node));
  }


}

#endif




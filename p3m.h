// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef P3M_H 
#define P3M_H
/** \file p3m.h   P3M algorithm for long range coulomb interaction.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  Further reading: 
 *  <ul>
 *  <li> P.P. Ewald,
 *       <i>Die Berechnung optischer und elektrostatischer Gitterpotentiale</i>,
 *       Ann. Phys. (64) 253-287, 1921
 *  <li> R. W. Hockney and J. W. Eastwood, 
 *       <i>Computer Simulation Using Particles</i>,
 *       IOP, London, 1988
 *  <li> M. Deserno and C. Holm,
 *       <i>How to mesh up {E}wald sums. I. + II.</i>,
 *       J. Chem. Phys. (109) 7678, 1998; (109) 7694, 1998
 *  <li> M. Deserno, C. Holm and H. J. Limbach,
 *       <i>How to mesh up {E}wald sums. </i>,
 *       in Molecular Dynamics on Parallel Computers,
 *       Ed. R. Esser et al., World Scientific, Singapore, 2000
 *  <li> M. Deserno,
 *       <i>Counterion condensation for rigid linear polyelectrolytes</i>,
 *       PhdThesis, Universit{\"a}t Mainz, 2000
 *  </ul>
 *
 *  For more information about the p3m algorithm,
 *  see \ref p3m.c "p3m.c"
 */

#include "utils.h"
#include "interaction_data.h"
#include "integrate.h"

#ifdef ELECTROSTATICS

/** This value for p3m.epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0

/************************************************
 * data types
 ************************************************/

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
  /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double alpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double r_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    mesh[3];
  /** offset of the first mesh point (lower left 
      corner) from the coordinate origin ([0,1[). */
  double mesh_off[3];
  /** charge assignment order ([0,7]). */
  int    cao;
  /** number of interpolation points for charge assignment function */
  int    inter;
  /** Accuracy of the actual parameter set. */
  double accuracy;

  /** epsilon of the "surrounding dielectric". */
  double epsilon;
  /** Cutoff for charge assignment. */
  double cao_cut[3];
  /** mesh constant. */
  double a[3];
  /** inverse mesh constant. */
  double ai[3];
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double alpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double r_cut;
  /** full size of the interpolated assignment function */
  int inter2;
  /** number of points unto which a single charge is interpolated, i.e. p3m.cao^3 */
  int cao3;
  /** additional points around the charge assignment mesh, for method like dielectric ELC
      creating virtual charges. */
  double additional_mesh[3];
} p3m_struct;

/** Structure for local mesh parameters. */
typedef struct {
  /* local mesh characterization. */
  /** dimension (size) of local mesh. */
  int dim[3];
  /** number of local mesh points. */
  int size;
  /** index of lower left corner of the 
      local mesh in the global mesh. */
  int ld_ind[3];
  /** position of the first local mesh point. */
  double ld_pos[3];
  /** dimension of mesh inside node domain. */
  int inner[3];
  /** inner left down grid point */
  int in_ld[3];
  /** inner up right grid point + (1,1,1)*/
  int in_ur[3];
  /** number of margin mesh points. */
  int margin[6];
  /** number of margin mesh points from neighbour nodes */
  int r_margin[6];
  /** offset between mesh lines of the last dimension */
  int q_2_off;
  /** offset between mesh lines of the two last dimensions */
  int q_21_off;
} local_mesh;

/** Structure for send/recv meshs. */
typedef struct {
  /** dimension of sub meshs to send. */
  int s_dim[6][3];
  /** left down corners of sub meshs to send. */
  int s_ld[6][3];
  /** up right corners of sub meshs to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6];
  /** dimensionof sub meshs to recv. */
  int r_dim[6][3];
  /** left down corners of sub meshs to recv. */
  int r_ld[6][3];
  /** up right corners of sub meshs to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6];
  /** maximal size for send/recv buffers. */
  int max;
} send_mesh;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** P3M parameters. */
extern p3m_struct p3m;
/** local mesh. */
extern local_mesh lm;
/** send/recv mesh sizes */
extern send_mesh  sm;

//The following are defined in p3m.c :

#define CA_INCREMENT 32  
extern int p3m_sum_qpart;
extern double p3m_sum_q2;
extern double p3m_square_sum_q;
extern int ca_num;
extern double  *ca_frac;
extern int *ca_fmp;
extern double *rs_mesh;
extern void realloc_ca_fields(int newsize);

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/// print the p3m parameters to the interpreters result
int printP3MToResult(Tcl_Interp *interp);

/// parse the optimization parameters of p3m and the tuner
int inter_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/// parse the basic p3m parameters
int inter_parse_p3m(Tcl_Interp * interp, int argc, char ** argv);

/// sanity checks
int P3M_sanity_checks();

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm.
 */
void P3M_init();

/** Calculate number of charged particles, the sum of the squared
    charges and the squared sum of the charges. */
void P3M_count_charged_particles();

/** Updates \ref p3m_struct::alpha and \ref p3m_struct::r_cut if \ref box_l changed. */
void P3M_scaleby_box_l();

/** Tune P3M parameters to desired accuracy.

    Usage:
    \verbatim inter coulomb <bjerrum> p3m tune accuracy <value> [r_cut <value> mesh <value> cao <value>] \endverbatim

    The parameters are tuned to obtain the desired accuracy in best
    time, by running mpi_integrate(0) for several parameter sets.

    The function utilizes the analytic expression of the error estimate 
    for the P3M method in the book of Hockney and Eastwood (Eqn. 8.23) in 
    order to obtain the rms error in the force for a system of N randomly 
    distributed particles in a cubic box.
    For the real space error the estimate of Kolafa/Perram is used. 

    Parameter range if not given explicit values: For \ref p3m_struct::r_cut_iL
    the function uses the values (\ref min_local_box_l -\ref #skin) /
    (n * \ref box_l), n being an integer (this implies the assumption that \ref
    p3m_struct::r_cut_iL is the largest cutoff in the system!). For \ref
    p3m_struct::mesh the function uses the two values which matches best the
    equation: number of mesh point = number of charged particles. For
    \ref p3m_struct::cao the function considers all possible values.

    For each setting \ref p3m_struct::alpha_L is calculated assuming that the
    error contributions of real and reciprocal space should be equal.

    After checking if the total error fulfils the accuracy goal the
    time needed for one force calculation (including verlet list
    update) is measured via \ref mpi_integrate(0).

    The function returns a log of the performed tuning.

    The function is based on routines of the program HE_Q.c written by M. Deserno.
 */
int P3M_tune_parameters(Tcl_Interp *interp);

/** a probably faster adaptive tuning method. Uses the same error estimates and parameters as
    \ref P3M_adaptive_tune_parameters, but a different strategy for finding the optimum. The algorithm
    basically determines the mesh, cao and then the real space cutoff, in this nested order.

    For each mesh, the cao optimal for the mesh tested previously is used as an initial guess,
    and the algorithm tries whether increasing or decreasing it leads to a better solution. This
    is efficient, since the optimal cao only changes little with the meshes in general.

    The real space cutoff for a given mesh and cao is determined via a bisection on the error estimate,
    which determines where the error estimate equals the required accuracy. Therefore the smallest 
    possible, i.e. fastest real space cutoff is determined.

    Both the search over mesh and cao stop to search in a specific direction once the computation time is
    significantly higher than the currently known optimum.

    Compared to \ref P3M_tune_parameters, this function will test more parameters sets for efficiency, but
    the error estimate is calculated less often. In general this should be faster and give better results.
 */
int P3M_adaptive_tune_parameters(Tcl_Interp *interp);

/** assign the physical charges using the tabulated charge assignment function.
    If store_ca_frac is true, then the charge fractions are buffered in cur_ca_fmp and
    cur_ca_frac. */
void P3M_charge_assign();

/** assign a single charge into the current charge grid. cp_cnt gives the a running index,
    which may be smaller than 0, in which case the charge is assumed to be virtual and is not
    stored in the ca_frac arrays. */
MDINLINE void P3M_assign_charge(double q,
				double real_pos[3],
#ifdef DIPOLES
				double mu,
				double dip[3],
#endif
				int cp_cnt)
{
  /* we do not really want to export these, but this function should be inlined */
  double P3M_caf(int i, double x);
  void realloc_ca_fields(int size);

  extern int    *ca_fmp;
  extern double *ca_frac;
  extern double *int_caf[7];
  extern double pos_shift;
  extern double *rs_mesh;
#ifdef DIPOLES
  extern double *rs_mesh_dip[3];
#endif

  int d, i0, i1, i2;
  double tmp0, tmp1;
  /* position of a particle in local mesh units */
  double pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* distance to nearest mesh point */
  double dist[3];
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  int q_ind = 0;
  double cur_ca_frac_val, *cur_ca_frac;

  // make sure we have enough space
  if (cp_cnt >= ca_num) realloc_ca_fields(cp_cnt + 1);
  // do it here, since realloc_ca_fields may change the address of ca_frac
  cur_ca_frac = ca_frac + p3m.cao3*cp_cnt;

  if (p3m.inter == 0) {
    for(d=0;d<3;d++) {
      /* particle position in mesh coordinates */
      pos    = ((real_pos[d]-lm.ld_pos[d])*p3m.ai[d]) - pos_shift;
      /* nearest mesh point */
      nmp  = (int)pos;
      /* distance to nearest mesh point */
      dist[d] = (pos-nmp)-0.5;
      /* 3d-array index of nearest mesh point */
      q_ind = (d == 0) ? nmp : nmp + lm.dim[d]*q_ind;

#ifdef ADDITIONAL_CHECKS
      if( pos < -skin*p3m.ai[d] ) {
	fprintf(stderr,"%d: rs_mesh underflow! (pos %f)\n", this_node, real_pos[d]);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node,my_left[d] - skin, my_right[d] + skin);	    
      }
      if( (nmp + p3m.cao) > lm.dim[d] ) {
	fprintf(stderr,"%d: rs_mesh overflow! (pos %f, nmp=%d)\n", this_node, real_pos[d],nmp);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node, my_left[d] - skin, my_right[d] + skin);
      }
#endif
    }
    if (cp_cnt >= 0) ca_fmp[cp_cnt] = q_ind;
    
    for(i0=0; i0<p3m.cao; i0++) {
      tmp0 = P3M_caf(i0, dist[0]);
      for(i1=0; i1<p3m.cao; i1++) {
	tmp1 = tmp0 * P3M_caf(i1, dist[1]);
	for(i2=0; i2<p3m.cao; i2++) {
#ifdef DIPOLES
	  cur_ca_frac_val = tmp1 * P3M_caf(i2, dist[2]);
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  if (q  != 0.0) rs_mesh[q_ind] += q * cur_ca_frac_val;
	  if (mu != 0.0) {
	    rs_mesh_dip[0][q_ind] += dip[0] * cur_ca_frac_val;
	    rs_mesh_dip[1][q_ind] += dip[1] * cur_ca_frac_val;
	    rs_mesh_dip[2][q_ind] += dip[2] * cur_ca_frac_val;
	  }
#else
	  // In the case without dipoles, ca_frac[] contains an additional factor q!
	  cur_ca_frac_val = q * tmp1 * P3M_caf(i2, dist[2]);
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  rs_mesh[q_ind] += cur_ca_frac_val;
#endif
	  q_ind++;
	}
	q_ind += lm.q_2_off;
      }
      q_ind += lm.q_21_off;
    }
  }
  else {
    /* particle position in mesh coordinates */
    for(d=0;d<3;d++) {
      pos    = ((real_pos[d]-lm.ld_pos[d])*p3m.ai[d]) - pos_shift;
      nmp    = (int) pos;
      arg[d] = (int) ((pos - nmp)*p3m.inter2);
      /* for the first dimension, q_ind is always zero, so this shifts correctly */
      q_ind = nmp + lm.dim[d]*q_ind;

#ifdef ADDITIONAL_CHECKS
      if( pos < -skin*p3m.ai[d] ) {
	fprintf(stderr,"%d: rs_mesh underflow! (pos %f)\n", this_node, real_pos[d]);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node,my_left[d] - skin, my_right[d] + skin);	    
      }
      if( (nmp + p3m.cao) > lm.dim[d] ) {
	fprintf(stderr,"%d: rs_mesh overflow! (pos %f, nmp=%d)\n", this_node, real_pos[d],nmp);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node, my_left[d] - skin, my_right[d] + skin);
      }
#endif
    }
    if (cp_cnt >= 0) ca_fmp[cp_cnt] = q_ind;

    for(i0=0; i0<p3m.cao; i0++) {
      tmp0 = int_caf[i0][arg[0]];
      for(i1=0; i1<p3m.cao; i1++) {
	tmp1 = tmp0 * int_caf[i1][arg[1]];
	for(i2=0; i2<p3m.cao; i2++) {
#ifdef DIPOLES
	  cur_ca_frac_val = tmp1 * int_caf[i2][arg[2]];
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  if (q  != 0.0) rs_mesh[q_ind] += q * cur_ca_frac_val;
	  if (mu != 0.0) {
	    rs_mesh_dip[0][q_ind] += dip[0] * cur_ca_frac_val;
	    rs_mesh_dip[1][q_ind] += dip[1] * cur_ca_frac_val;
	    rs_mesh_dip[2][q_ind] += dip[2] * cur_ca_frac_val;
	  }
#else
	  //In the case without dipoles, ca_frac[] contains an additional factor q!
	  cur_ca_frac_val = q * tmp1 * int_caf[i2][arg[2]];
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  rs_mesh[q_ind] += cur_ca_frac_val;
#endif
	  q_ind++;
	}
	q_ind += lm.q_2_off;
      }
      q_ind += lm.q_21_off;
    }
  }
}

/** shrink wrap the charge grid */
MDINLINE void P3M_shrink_wrap_charge_grid(int n_charges) {
  /* we do not really want to export these */
  void realloc_ca_fields(int size);
  
  if( n_charges < ca_num ) realloc_ca_fields(n_charges);
}

/** Calculate the k-space contribution to the coulomb interaction
    forces. */ 
double P3M_calc_kspace_forces(int force_flag, int energy_flag);

/** Calculate real space contribution of coulomb pair forces.
    If NPT is compiled in, it returns the energy, which is needed for NPT. */
MDINLINE double add_p3m_coulomb_pair_force(double chgfac, double *d,double dist2,double dist,double force[3])
{
  int j;
  double fac1,fac2, adist, erfc_part_ri;
  if(dist < p3m.r_cut) {
    if (dist > 0.0){		//Vincent
      adist = p3m.alpha * dist;
#if USE_ERFC_APPROXIMATION
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac1 = coulomb.prefactor * chgfac  * exp(-adist*adist);
      fac2 = fac1 * (erfc_part_ri + 2.0*p3m.alpha*wupii) / dist2;
#else
      erfc_part_ri = erfc(adist) / dist;
      fac1 = coulomb.prefactor * chgfac;
      fac2 = fac1 * (erfc_part_ri + 2.0*p3m.alpha*wupii*exp(-adist*adist)) / dist2;
#endif
      for(j=0;j<3;j++)
	force[j] += fac2 * d[j];
      ESR_TRACE(fprintf(stderr,"%d: RSE: Pair dist=%.3f: force (%.3e,%.3e,%.3e)\n",this_node,
			dist,fac2*d[0],fac2*d[1],fac2*d[2]));
#ifdef NPT
      return fac1 * erfc_part_ri;
#endif
    }
  }
  return 0.0;
}

#ifdef DIPOLES

/** Calculate real space contribution of p3m dipolar pair forces and torques.
    If NPT is compiled in, it returns the energy, which is needed for NPT. */
MDINLINE double add_p3m_dipolar_pair_force(Particle *p1, Particle *p2,
					   double *d,double dist2,double dist,double force[3])
{
  int j;
  double fac1,fac2, adist, erfc_part_ri, coeff, exp_adist2, dist2i;
  double mimj, mir, mjr;
  double B_r, C_r, D_r;
  double alpsq = p3m.alpha * p3m.alpha;
  double mixmj[3], mixr[3], mjxr[3];


  if(dist < p3m.r_cut) {
  if(dist > 0){
    adist = p3m.alpha * dist;
#if USE_ERFC_APPROXIMATION
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac1 = coulomb.prefactor * p1->p.dipm*p2->p.dipm * exp(-adist*adist);
      fac2 = fac1 * (erfc_part_ri + 2.0*p3m.alpha*wupii) / dist2;
#else
      erfc_part_ri = erfc(adist) / dist;
      fac1 = coulomb.prefactor * p1->p.dipm*p2->p.dipm;
      fac2 = fac1 * (erfc_part_ri + 2.0*p3m.alpha*wupii*exp(-adist*adist)) / dist2;
#endif

  //Calculate scalar multiplications for vectors mi, mj, rij
  mimj = p1->r.dip[0]*p2->r.dip[0] + p1->r.dip[1]*p2->r.dip[1] + p1->r.dip[2]*p2->r.dip[2];
  mir = p1->r.dip[0]*d[0] + p1->r.dip[1]*d[1] + p1->r.dip[2]*d[2];
  mjr = p2->r.dip[0]*d[0] + p2->r.dip[1]*d[1] + p2->r.dip[2]*d[2];

  coeff = 2.0*p3m.alpha*wupii;
  dist2i = 1 / dist2;
  exp_adist2 = exp(-adist*adist);

  if(p3m.accuracy > 5e-06)
    B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
  else
    B_r = (erfc(adist)/dist + coeff * exp_adist2) * dist2i;
  C_r = (3*B_r + 2*alpsq*coeff*exp_adist2) * dist2i;
  D_r = (5*C_r + 4*coeff*alpsq*alpsq*exp_adist2) * dist2i;

  // Calculate real-space forces
  for(j=0;j<3;j++)
    force[j] += (mimj*d[j] + p1->r.dip[j]*mjr + p2->r.dip[j]*mir) * C_r - mir*mjr*D_r*d[j] ;

  //Calculate vector multiplications for vectors mi, mj, rij

  mixmj[0] = p1->r.dip[1]*p2->r.dip[2] - p1->r.dip[2]*p2->r.dip[1];
  mixmj[1] = p1->r.dip[2]*p2->r.dip[0] - p1->r.dip[0]*p2->r.dip[2];
  mixmj[2] = p1->r.dip[0]*p2->r.dip[1] - p1->r.dip[1]*p2->r.dip[0];

  mixr[0] = p1->r.dip[1]*d[2] - p1->r.dip[2]*d[1];
  mixr[1] = p1->r.dip[2]*d[0] - p1->r.dip[0]*d[2];
  mixr[2] = p1->r.dip[0]*d[1] - p1->r.dip[1]*d[0];

  mjxr[0] = p2->r.dip[1]*d[2] - p2->r.dip[2]*d[1];
  mjxr[1] = p2->r.dip[2]*d[0] - p2->r.dip[0]*d[2];
  mjxr[2] = p2->r.dip[0]*d[1] - p2->r.dip[1]*d[0];

  // Calculate real-space torques

  for(j=0;j<3;j++){
    p1->f.torque[j] += -mixmj[j]*B_r + mixr[j]*mjr*C_r;
    p2->f.torque[j] +=  mixmj[j]*B_r + mjxr[j]*mir*C_r;
  }

#ifdef NPT
  return fac1 * ( mimj*B_r - mir*mjr * C_r );
#endif
  }}
  return 0.0;
}
#endif  /* ifdef DIPOLES */

/** Calculate real space contribution of coulomb pair energy. */
MDINLINE double p3m_coulomb_pair_energy(double chgfac, double *d,double dist2,double dist)
{
  double adist, erfc_part_ri;

  if(dist < p3m.r_cut) {
    adist = p3m.alpha * dist;
#if USE_ERFC_APPROXIMATION
    erfc_part_ri = AS_erfc_part(adist) / dist;
    return coulomb.prefactor*chgfac*erfc_part_ri*exp(-adist*adist);
#else
    erfc_part_ri = erfc(adist) / dist;
    return coulomb.prefactor*chgfac*erfc_part_ri;
#endif
  }
  return 0.0;
}

#ifdef DIPOLES
/** Calculate real space contribution of dipolar pair energy. */
MDINLINE double p3m_dipolar_pair_energy(Particle *p1, Particle *p2,
				      double *d,double dist2,double dist)
{
  double fac1, adist, erfc_part_ri, coeff, exp_adist2, dist2i;
  double mimj, mir, mjr;
  double B_r, C_r;
  double alpsq = p3m.alpha * p3m.alpha;
 
  if(dist < p3m.r_cut) {
  if(dist > 0){
    adist = p3m.alpha * dist;
#if USE_ERFC_APPROXIMATION
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac1 = coulomb.prefactor * p1->p.dipm*p2->p.dipm; /* *exp(-adist*adist); */
#else
      erfc_part_ri = erfc(adist) / dist;
      fac1 = coulomb.prefactor * p1->p.dipm*p2->p.dipm;
#endif

  //Calculate scalar multiplications for vectors mi, mj, rij
  mimj = p1->r.dip[0]*p2->r.dip[0] + p1->r.dip[1]*p2->r.dip[1] + p1->r.dip[2]*p2->r.dip[2];
  mir = p1->r.dip[0]*d[0] + p1->r.dip[1]*d[1] + p1->r.dip[2]*d[2];
  mjr = p2->r.dip[0]*d[0] + p2->r.dip[1]*d[1] + p2->r.dip[2]*d[2];

  coeff = 2.0*p3m.alpha*wupii;
  dist2i = 1 / dist2;
  exp_adist2 = exp(-adist*adist);

  if(p3m.accuracy > 5e-06)
    B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
  else
    B_r = (erfc(adist)/dist + coeff * exp_adist2) * dist2i;
  
  C_r = (3*B_r + 2*alpsq*coeff*exp_adist2) * dist2i;

  /*
  printf("(%4i %4i) pair energy = %f (B_r=%15.12f C_r=%15.12f)\n",p1->p.identity,p2->p.identity,fac1*(mimj*B_r-mir*mjr*C_r),B_r,C_r);
  */
  
  return fac1 * ( mimj*B_r - mir*mjr * C_r );

  }}
  return 0.0;
}
#endif /* ifdef DIPOLES */

/** Clean up P3M memory allocations. */
void P3M_exit();

/*@}*/
#endif

#endif

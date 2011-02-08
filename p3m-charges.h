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

/** \file p3m-charges.h   P3M algorithm for long range coulomb interaction.
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
 *  see \ref p3m.h "p3m.h"
 *  see \ref p3m-charges.c  "p3m-charges.c"
 *  see \ref p3m-charges.h  "p3m-charges.h"
 *  see \ref p3m-dipoles.c  "p3m-dipoles.c"
 *  see \ref p3m-dipoles.h  "p3m-dipoles.h"
 *  see \ref p3m-assignement.c  "p3m-assignement.c"
 */

#ifdef ELECTROSTATICS

/* only include from within p3m.h */
#ifndef P3M_H_CURRENT
#error never include this file file directly, include p3m.h
#endif

/** local mesh. */
extern local_mesh lm;
/** send/recv mesh sizes */
extern send_mesh  sm;

//The following are defined in p3m-charges.c :

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

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm for charge-charge interactions.
 */
void  P3M_init_charges(void);

/** Updates \ref p3m_struct::alpha and \ref p3m_struct::r_cut if \ref box_l changed. */
void P3M_scaleby_box_l_charges();

/// parse the optimization parameters of p3m and the tuner
int tclcommand_inter_coulomb_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/// parse the basic p3m parameters
int tclcommand_inter_coulomb_parse_p3m(Tcl_Interp * interp, int argc, char ** argv);

/** compute the k-space part of forces and energies for the charge-charge interaction  **/
double P3M_calc_kspace_forces_for_charges(int force_flag, int energy_flag);

/** computer the k-space part of the stress tensor **/
void P3M_calc_kspace_stress (double* stress);

/// sanity checks
int P3M_sanity_checks();

/** checks for correctness for charges in P3M of the cao_cut, necessary when the box length changes */
int P3M_sanity_checks_boxl();

/** Calculate number of charged particles, the sum of the squared
    charges and the squared sum of the charges. */
void P3M_count_charged_particles();

/** Error Codes for p3m tuning (version 2) :
    P3M_TUNE_FAIL: force evaluation failes,
    P3M_TUNE_NO_CUTOFF: could not finde a valid realspace cutoff radius,
    P3M_TUNE_CAOTOLARGE: Charge asignment order to large for mesh size,
    P3M_TUNE_ELCTEST: conflict with ELC gap size.
*/

enum P3M_TUNE_ERROR { P3M_TUNE_FAIL = 1, P3M_TUNE_NOCUTOFF = 2, P3M_TUNE_CAOTOLARGE = 4, P3M_TUNE_ELCTEST = 8, P3M_TUNE_CUTOFF_TOO_LARGE = 16 };

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
int tclcommand_inter_coulomb_print_p3m_tune_parameteres(Tcl_Interp *interp);

/** a probably faster adaptive tuning method. Uses the same error estimates and parameters as
    \ref tclcommand_inter_coulomb_print_p3m_adaptive_tune_parameteres, but a different strategy for finding the optimum. The algorithm
    basically determines the mesh, cao and then the real space cutoff, in this nested order.

    For each mesh, the cao optimal for the mesh tested previously is used as an initial guess,
    and the algorithm tries whether increasing or decreasing it leads to a better solution. This
    is efficient, since the optimal cao only changes little with the meshes in general.

    The real space cutoff for a given mesh and cao is determined via a bisection on the error estimate,
    which determines where the error estimate equals the required accuracy. Therefore the smallest 
    possible, i.e. fastest real space cutoff is determined.

    Both the search over mesh and cao stop to search in a specific direction once the computation time is
    significantly higher than the currently known optimum.

    Compared to \ref tclcommand_inter_coulomb_print_p3m_tune_parameteres, this function will test more parameters sets for efficiency, but
    the error estimate is calculated less often. In general this should be faster and give better results.
 */
int tclcommand_inter_coulomb_print_p3m_adaptive_tune_parameteres(Tcl_Interp *interp);

/** assign the physical charges using the tabulated charge assignment function.
    If store_ca_frac is true, then the charge fractions are buffered in cur_ca_fmp and
    cur_ca_frac. */
void P3M_charge_assign();

/** assign a single charge into the current charge grid. cp_cnt gives the a running index,
    which may be smaller than 0, in which case the charge is assumed to be virtual and is not
    stored in the ca_frac arrays. */
    
MDINLINE void P3M_assign_charge(double q,
				double real_pos[3],
				int cp_cnt)
{
  /* we do not really want to export these, but this function should be inlined */
  double P3M_caf(int i, double x,int cao_value);
  void realloc_ca_fields(int size);

  extern int    *ca_fmp;
  extern double *ca_frac;
  extern double *int_caf[7];
  extern double pos_shift;
  extern double *rs_mesh;

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
      tmp0 = P3M_caf(i0, dist[0],p3m.cao);
      for(i1=0; i1<p3m.cao; i1++) {
	tmp1 = tmp0 * P3M_caf(i1, dist[1],p3m.cao);
	for(i2=0; i2<p3m.cao; i2++) {
	  cur_ca_frac_val = q * tmp1 * P3M_caf(i2, dist[2], p3m.cao);
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  rs_mesh[q_ind] += cur_ca_frac_val;
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
	  cur_ca_frac_val = q * tmp1 * int_caf[i2][arg[2]];
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  rs_mesh[q_ind] += cur_ca_frac_val;
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

#endif

/*@}*/

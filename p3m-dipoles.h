// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

/** \file p3m-dipoles.h   P3M algorithm for long range magnetic dipole-dipole interaction.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the dipolar Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  Further reading: 
 *  <ul>
 *  <li> J.J. Cerda, P3M for dipolar interactions. J. Chem. Phys, 129, xxx ,(2008).
 *  </ul>
 *
 *  For more information about the p3m algorithm,
 *  see \ref p3m.h "p3m.h"
 *  see \ref p3m-charges.c  "p3m-charges.c"
 *  see \ref p3m-charges.h  "p3m-charges.h"
 *  see \ref p3m-dipoles.c  "p3m-dipoles.c"
 *  see \ref p3m-dipoles.h  "p3m-dipoles.h"
 *  see \ref p3m-assignment.c  "p3m-assignment.c"
 */

#ifdef MAGNETOSTATICS

/* only include from within p3m.h */
#ifndef P3M_H_CURRENT
#error never include this file file directly, include p3m.h
#endif

/** local mesh. */
extern local_mesh Dlm;
/** send/recv mesh sizes */
extern send_mesh  Dsm;

//The following are defined in p3m-dipoles.c :

extern int Dca_num;
extern double  *Dca_frac;
extern int *Dca_fmp;
extern double *Drs_mesh;
extern void Drealloc_ca_fields(int newsize);

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm for dipole-dipole interactions.
 */
void  P3M_init_dipoles(void);

/** Updates \ref p3m_struct::alpha and \ref p3m_struct::r_cut if \ref box_l changed. */
void P3M_scaleby_box_l_dipoles();


/// parse the optimization parameters of p3m and the tuner
int Dinter_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/// parse the basic p3m parameters
int Dinter_parse_p3m(Tcl_Interp * interp, int argc, char ** argv);


/// sanity checks
int DP3M_sanity_checks();

/** checks for correctness for magnetic dipoles in P3M of the cao_cut, necessary when the box length changes */
int DP3M_sanity_checks_boxl();


/** Tune dipolar P3M parameters to desired accuracy.

    Usage:
    \verbatim inter dipolar <bjerrum> p3m tune accuracy <value> [r_cut <value> mesh <value> cao <value>] \endverbatim

    The parameters are tuned to obtain the desired accuracy in best
    time, by running mpi_integrate(0) for several parameter sets.

    The function utilizes the analytic expression of the error estimate 
    for the dipolar P3M method see JCP,2008 paper by J.J.Cerda et al in 
    order to obtain the rms error in the force for a system of N randomly 
    distributed particles in a cubic box.
    For the real space error the estimate of Kolafa/Perram is used. 

    Parameter range if not given explicit values: For \ref p3m_struct::Dr_cut_iL
    the function uses the values (\ref min_local_box_l -\ref #skin) /
    (n * \ref box_l), n being an integer (this implies the assumption that \ref
    p3m_struct::Dr_cut_iL is the largest cutoff in the system!). For \ref
    p3m_struct::Dmesh the function uses the two values which matches best the
    equation: number of mesh point = number of magnetic dipolar particles. For
    \ref p3m_struct::cao the function considers all possible values.

    For each setting \ref p3m_struct::Dalpha_L is calculated assuming that the
    error contributions of real and reciprocal space should be equal.

    After checking if the total error fulfils the accuracy goal the
    time needed for one force calculation (including verlet list
    update) is measured via \ref mpi_integrate(0).

    The function returns a log of the performed tuning.

    The function is based on routines for charges.
 */
int DP3M_tune_parameters(Tcl_Interp *interp);

/** a probably faster adaptive tuning method. Uses the same error estimates and parameters as
    \ref DP3M_adaptive_tune_parameters, but a different strategy for finding the optimum. The algorithm
    basically determines the mesh, cao and then the real space cutoff, in this nested order.

    For each mesh, the cao optimal for the mesh tested previously is used as an initial guess,
    and the algorithm tries whether increasing or decreasing it leads to a better solution. This
    is efficient, since the optimal cao only changes little with the meshes in general.

    The real space cutoff for a given mesh and cao is determined via a bisection on the error estimate,
    which determines where the error estimate equals the required accuracy. Therefore the smallest 
    possible, i.e. fastest real space cutoff is determined.

    Both the search over mesh and cao stop to search in a specific direction once the computation time is
    significantly higher than the currently known optimum.

    Compared to \ref DP3M_tune_parameters, this function will test more parameters sets for efficiency, but
    the error estimate is calculated less often. In general this should be faster and give better results.
 */
int DP3M_adaptive_tune_parameters(Tcl_Interp *interp);

/** assign the physical dipoles using the tabulated assignment function.
    If Dstore_ca_frac is true, then the charge fractions are buffered in Dcur_ca_fmp and
    Dcur_ca_frac. */
void P3M_dipole_assign();


/** compute the k-space part of forces and energies for the magnetic dipole-dipole interaction  */
double P3M_calc_kspace_forces_for_dipoles(int force_flag, int energy_flag);


/** Calculate number of magnetic  particles, the sum of the squared
    charges and the squared sum of the charges. */

void P3M_count_magnetic_particles();


/** assign a single dipole into the current charge grid. cp_cnt gives the a running index,
    which may be smaller than 0, in which case the charge is assumed to be virtual and is not
    stored in the Dca_frac arrays. */
    
MDINLINE void P3M_assign_dipole(double real_pos[3],double mu, double dip[3],int cp_cnt)
{
  /* we do not really want to export these, but this function should be inlined */
  double P3M_caf(int i, double x, int cao_value);
  void Drealloc_ca_fields(int size);

  extern int    *Dca_fmp;
  extern double *Dca_frac;
  extern double *Dint_caf[7];
  extern double Dpos_shift;
  extern double *Drs_mesh_dip[3];

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
  /* index, index jumps for Drs_mesh array */
  int q_ind = 0;
  double cur_ca_frac_val, *cur_ca_frac;

  // make sure we have enough space
  if (cp_cnt >= Dca_num) Drealloc_ca_fields(cp_cnt + 1);
  // do it here, since realloc_ca_fields may change the address of Dca_frac
  cur_ca_frac = Dca_frac + p3m.Dcao3*cp_cnt;

  if (p3m.Dinter == 0) {
    for(d=0;d<3;d++) {
      /* particle position in mesh coordinates */
      pos    = ((real_pos[d]-Dlm.ld_pos[d])*p3m.Dai[d]) - Dpos_shift;
      /* nearest mesh point */
      nmp  = (int)pos;
      /* distance to nearest mesh point */
      dist[d] = (pos-nmp)-0.5;
      /* 3d-array index of nearest mesh point */
      q_ind = (d == 0) ? nmp : nmp + Dlm.dim[d]*q_ind;

#ifdef ADDITIONAL_CHECKS
      if( pos < -skin*p3m.Dai[d] ) {
	fprintf(stderr,"%d: dipolar Drs_mesh underflow! (pos %f)\n", this_node, real_pos[d]);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node,my_left[d] - skin, my_right[d] + skin);	    
      }
      if( (nmp + p3m.Dcao) > Dlm.dim[d] ) {
	fprintf(stderr,"%d: dipolar Drs_mesh overflow! (pos %f, nmp=%d)\n", this_node, real_pos[d],nmp);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node, my_left[d] - skin, my_right[d] + skin);
      }
#endif
    }
    if (cp_cnt >= 0) Dca_fmp[cp_cnt] = q_ind;
    
    for(i0=0; i0<p3m.Dcao; i0++) {
      tmp0 = P3M_caf(i0, dist[0], p3m.Dcao);
      for(i1=0; i1<p3m.Dcao; i1++) {
	tmp1 = tmp0 * P3M_caf(i1, dist[1],p3m.Dcao);
	for(i2=0; i2<p3m.Dcao; i2++) {
	  cur_ca_frac_val = tmp1 * P3M_caf(i2, dist[2],p3m.Dcao);
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  if (mu != 0.0) {
	    Drs_mesh_dip[0][q_ind] += dip[0] * cur_ca_frac_val;
	    Drs_mesh_dip[1][q_ind] += dip[1] * cur_ca_frac_val;
	    Drs_mesh_dip[2][q_ind] += dip[2] * cur_ca_frac_val;
	  }
	  q_ind++;
	}
	q_ind += Dlm.q_2_off;
      }
      q_ind += Dlm.q_21_off;
    }
  }
  else {
    /* particle position in mesh coordinates */
    for(d=0;d<3;d++) {
      pos    = ((real_pos[d]-Dlm.ld_pos[d])*p3m.Dai[d]) - Dpos_shift;
      nmp    = (int) pos;
      arg[d] = (int) ((pos - nmp)*p3m.Dinter2);
      /* for the first dimension, q_ind is always zero, so this shifts correctly */
      q_ind = nmp + Dlm.dim[d]*q_ind;

#ifdef ADDITIONAL_CHECKS
      if( pos < -skin*p3m.Dai[d] ) {
	fprintf(stderr,"%d: dipolar Drs_mesh underflow! (pos %f)\n", this_node, real_pos[d]);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node,my_left[d] - skin, my_right[d] + skin);	    
      }
      if( (nmp + p3m.Dcao) > Dlm.dim[d] ) {
	fprintf(stderr,"%d: dipolar Drs_mesh overflow! (pos %f, nmp=%d)\n", this_node, real_pos[d],nmp);
	fprintf(stderr,"%d: allowed coordinates: %f - %f\n",
		this_node, my_left[d] - skin, my_right[d] + skin);
      }
#endif
    }
    if (cp_cnt >= 0) Dca_fmp[cp_cnt] = q_ind;

    for(i0=0; i0<p3m.Dcao; i0++) {
      tmp0 = Dint_caf[i0][arg[0]];
      for(i1=0; i1<p3m.Dcao; i1++) {
	tmp1 = tmp0 * Dint_caf[i1][arg[1]];
	for(i2=0; i2<p3m.Dcao; i2++) {
	  cur_ca_frac_val = tmp1 * Dint_caf[i2][arg[2]];
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  if (mu != 0.0) {
	    Drs_mesh_dip[0][q_ind] += dip[0] * cur_ca_frac_val;
	    Drs_mesh_dip[1][q_ind] += dip[1] * cur_ca_frac_val;
	    Drs_mesh_dip[2][q_ind] += dip[2] * cur_ca_frac_val;
	  }
	  q_ind++;
	}
	q_ind += Dlm.q_2_off;
      }
      q_ind += Dlm.q_21_off;
    }
  }
}

/** shrink wrap the dipoles grid */
MDINLINE void DP3M_shrink_wrap_dipole_grid(int n_dipoles) {
  /* we do not really want to export these */
  void Drealloc_ca_fields(int size);
  
  if( n_dipoles < Dca_num ) Drealloc_ca_fields(n_dipoles);
}

/** Calculate real space contribution of p3m dipolar pair forces and torques.
    If NPT is compiled in, it returns the energy, which is needed for NPT. */
MDINLINE double add_p3m_dipolar_pair_force(Particle *p1, Particle *p2,
					   double *d,double dist2,double dist,double force[3])
{
  int j;
  double fac1, adist, erfc_part_ri, coeff, exp_adist2, dist2i;
  double mimj, mir, mjr;
  double B_r, C_r, D_r;
  double alpsq = p3m.Dalpha * p3m.Dalpha;
  double mixmj[3], mixr[3], mjxr[3];

  if(dist < p3m.Dr_cut && dist > 0) {
    adist = p3m.Dalpha * dist;
    #if USE_ERFC_APPROXIMATION
       erfc_part_ri = AS_erfc_part(adist) / dist;
    #else
       erfc_part_ri = erfc(adist) / dist;
    #endif

  //Calculate scalar multiplications for vectors mi, mj, rij
  mimj = p1->r.dip[0]*p2->r.dip[0] + p1->r.dip[1]*p2->r.dip[1] + p1->r.dip[2]*p2->r.dip[2];
  mir = p1->r.dip[0]*d[0] + p1->r.dip[1]*d[1] + p1->r.dip[2]*d[2];
  mjr = p2->r.dip[0]*d[0] + p2->r.dip[1]*d[1] + p2->r.dip[2]*d[2];

  coeff = 2.0*p3m.Dalpha*wupii;
  dist2i = 1 / dist2;
  exp_adist2 = exp(-adist*adist);

  if(p3m.Daccuracy > 5e-06)
    B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
  else
    B_r = (erfc(adist)/dist + coeff * exp_adist2) * dist2i;
    
  C_r = (3*B_r + 2*alpsq*coeff*exp_adist2) * dist2i;
  D_r = (5*C_r + 4*coeff*alpsq*alpsq*exp_adist2) * dist2i;

  // Calculate real-space forces
  for(j=0;j<3;j++)
    force[j] += coulomb.Dprefactor *((mimj*d[j] + p1->r.dip[j]*mjr + p2->r.dip[j]*mir) * C_r - mir*mjr*D_r*d[j]) ;

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
    p1->f.torque[j] += coulomb.Dprefactor *(-mixmj[j]*B_r + mixr[j]*mjr*C_r);
    p2->f.torque[j] += coulomb.Dprefactor *( mixmj[j]*B_r + mjxr[j]*mir*C_r);
  }

#ifdef NPT
#if USE_ERFC_APPROXIMATION
  fac1 = coulomb.Dprefactor * p1->p.dipm*p2->p.dipm * exp(-adist*adist);
#else
  fac1 = coulomb.Dprefactor * p1->p.dipm*p2->p.dipm;
#endif
  return fac1 * ( mimj*B_r - mir*mjr * C_r );
#endif
  }
  return 0.0;
}

/** Calculate real space contribution of dipolar pair energy. */
MDINLINE double p3m_dipolar_pair_energy(Particle *p1, Particle *p2,
				      double *d,double dist2,double dist)
{
  double /* fac1,*/ adist, erfc_part_ri, coeff, exp_adist2, dist2i;
  double mimj, mir, mjr;
  double B_r, C_r;
  double alpsq = p3m.Dalpha * p3m.Dalpha;
 
  if(dist < p3m.Dr_cut && dist > 0) {
    adist = p3m.Dalpha * dist;
    /*fac1 = coulomb.Dprefactor;*/

#if USE_ERFC_APPROXIMATION
    erfc_part_ri = AS_erfc_part(adist) / dist;
    /*  fac1 = coulomb.Dprefactor * p1->p.dipm*p2->p.dipm; IT WAS WRONG */ /* *exp(-adist*adist); */ 
#else
    erfc_part_ri = erfc(adist) / dist;
    /* fac1 = coulomb.Dprefactor * p1->p.dipm*p2->p.dipm;  IT WAS WRONG*/
#endif

    //Calculate scalar multiplications for vectors mi, mj, rij
    mimj = p1->r.dip[0]*p2->r.dip[0] + p1->r.dip[1]*p2->r.dip[1] + p1->r.dip[2]*p2->r.dip[2];
    mir = p1->r.dip[0]*d[0] + p1->r.dip[1]*d[1] + p1->r.dip[2]*d[2];
    mjr = p2->r.dip[0]*d[0] + p2->r.dip[1]*d[1] + p2->r.dip[2]*d[2];

    coeff = 2.0*p3m.Dalpha*wupii;
    dist2i = 1 / dist2;
    exp_adist2 = exp(-adist*adist);

    if(p3m.Daccuracy > 5e-06)
      B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
    else
      B_r = (erfc(adist)/dist + coeff * exp_adist2) * dist2i;
  
    C_r = (3*B_r + 2*alpsq*coeff*exp_adist2) * dist2i;

    /*
      printf("(%4i %4i) pair energy = %f (B_r=%15.12f C_r=%15.12f)\n",p1->p.identity,p2->p.identity,fac1*(mimj*B_r-mir*mjr*C_r),B_r,C_r);
    */
  
    /* old line return fac1 * ( mimj*B_r - mir*mjr * C_r );*/
    return coulomb.Dprefactor * ( mimj*B_r - mir*mjr * C_r );
  }
  return 0.0;
}

/*@}*/
#endif /* of defined(MAGNETOSTATICS) */

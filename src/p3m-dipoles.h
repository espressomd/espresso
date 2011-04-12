/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

/** \file p3m-dipoles.h   P3M algorithm for long range magnetic dipole-dipole interaction.
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
/** dipolar p3m parser */
int tclcommand_inter_magnetic_parse_p3m(Tcl_Interp * interp, int argc, char ** argv);

/** dipolar p3m parser, optional parameters */
int tclcommand_inter_magnetic_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm for dipole-dipole interactions.
 */
void  P3M_init_dipoles(void);

/** Updates \ref p3m_struct::alpha and \ref p3m_struct::r_cut if \ref box_l changed. */
void P3M_scaleby_box_l_dipoles();

/// sanity checks
int DP3M_sanity_checks();

/** assign the physical dipoles using the tabulated assignment function.
    If Dstore_ca_frac is true, then the charge fractions are buffered in Dcur_ca_fmp and
    Dcur_ca_frac. */
void P3M_dipole_assign(void);


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
#ifdef NPT
  double fac1;
#endif
  double adist, erfc_part_ri, coeff, exp_adist2, dist2i;
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
#ifdef ROTATION
  for(j=0;j<3;j++){
    p1->f.torque[j] += coulomb.Dprefactor *(-mixmj[j]*B_r + mixr[j]*mjr*C_r);
    p2->f.torque[j] += coulomb.Dprefactor *( mixmj[j]*B_r + mjxr[j]*mir*C_r);
  }
#endif
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

/** \file forces.c Force calculation.
 *
 *  For more information see \ref forces.h "forces.h".
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "forces.h"
#include "debug.h"
#include "thermostat.h"
#include "communication.h"
#include "ghosts.h" 
#include "verlet.h"
#include "utils.h"
#include "interaction_data.h"
#include "grid.h"
#include "p3m.h"
/* include the force files */
#include "fene.h"

double minimum_part_dist = -1;

/************************************************************/

void force_init()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d interaction types\n",
		      this_node,n_interaction_types));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_particle_types));
}

/************************************************************/

void force_calc()
{
  int i,j;
  int ind1, ind2, ind3;
  double d[3], dist2, dist;
  IA_parameters *ia_params;
  extern Particle *particles;
  extern int *verletList;
  double frac2,frac6,r_off;
  double fac, adist;
  /* bonded interactions */
  int type_num;
  /* angle bend interaction */
  double vec1[3],vec2[3],d1i,d2i,cosine,f1,f2;
  /* electrostatic */
  double  erfc_part_ri;

  FORCE_TRACE(fprintf(stderr,"%d: force_calc: for %d (P %d,G %d)\n",this_node,n_particles+n_ghosts,n_particles,n_ghosts));

  minimum_part_dist = box_l[0] + box_l[1] + box_l[2];

  /* initialize forces with thermostat forces */
  friction_thermo();

  /* initialize ghost forces with zero */
  for(i=n_particles;i<n_particles+n_ghosts;i++)
    for(j=0;j<3;j++)
      particles[i].f[j] = 0.0;

  /* calculate bonded interactions (loop local particles) */
  for(ind1=0; ind1<n_particles; ind1++) {
    i=0;
    while(i<particles[ind1].n_bonds) {
      type_num = particles[ind1].bonds[i];
      switch(bonded_ia_params[type_num].type) {
      case BONDED_IA_FENE:
	ind2 = local_index[particles[ind1].bonds[i+1]];
	add_fene_pair_force(ind1, ind2, type_num);
	i+=2; break;
      case BONDED_IA_ANGLE:
	ind2 = local_index[particles[ind1].bonds[i+1]];
	ind3 = local_index[particles[ind1].bonds[i+2]];
	if(ind2 == -1 || ind3 == -1) {
	  fprintf(stderr,"WARNING: Bonded atoms %d, %d and %d not on the same node\n",
		  particles[ind1].identity,particles[ind1].bonds[i+1],particles[ind1].bonds[i+2]); 
	  i+=3; break;
	}
	/* angle bend force */
	cosine=0.0;
	/* vector from particles[ind2] to particles[ind1] */
	for(j=0;j<3;j++) vec1[j] = particles[ind2].p[j] - particles[ind1].p[j];
	dist2 = SQR(vec1[0]) + SQR(vec1[1]) + SQR(vec1[2]);
	d1i = 1.0 / sqrt(dist2);
	for(j=0;j<3;j++) vec1[j] *= d1i;
	/* vector from particles[ind3] to particles[ind1] */
 	for(j=0;j<3;j++) vec2[j] = particles[ind3].p[j] - particles[ind1].p[j];
	dist2 = SQR(vec2[0]) + SQR(vec2[1]) + SQR(vec2[2]);
	d2i = 1.0 / sqrt(dist2);
	for(j=0;j<3;j++) vec2[j] *= d2i;
	/* scalar produvt of vec1 and vec2 */
	for(j=0;j<3;j++) cosine += vec1[j] * vec2[j];
	/* apply bend forces */
	for(j=0;j<3;j++) {
	  f1 = bonded_ia_params[type_num].p.angle.bend * (vec2[j] - cosine * vec1[j]) * d1i;
	  f2 = bonded_ia_params[type_num].p.angle.bend * (vec1[j] - cosine * vec2[j]) * d2i;
	  particles[ind2].f[j] -= f1;
	  particles[ind1].f[j] += (f1+f2);
	  particles[ind3].f[j] -= f2;
	}
 	i+=3; break;
      default :
	fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",particles[ind1].identity);
	i = particles[ind1].n_bonds; 
	break;
     }
    }
  }

  /* calculate non bonded interactions (loop verlet list) */
  for(i=0;i<2*n_verletList;i=i+2) {
    ind1 = verletList[i];
    ind2 = verletList[i+1];
    ia_params = get_ia_param(particles[ind1].type,particles[ind2].type);
    for(j=0;j<3;j++)
      d[j] = particles[ind1].p[j] - particles[ind2].p[j];
    dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
    dist = sqrt(dist2);

    /* lennnard jones */

    if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
      r_off = dist - ia_params->LJ_offset;
      if(r_off>0.0) {
#ifdef LJ_WARN_WHEN_CLOSE
	if (r_off < 0.9*ia_params->LJ_sig) {
	  fprintf(stderr, "Lennard-Jones warning: particles getting close\n");
	}
#endif
	frac2 = SQR(ia_params->LJ_sig/r_off);
	frac6 = frac2*frac2*frac2;
	fac = 48.* ia_params->LJ_eps * frac6*(frac6 - 0.5)*frac2;
	for(j=0;j<3;j++) {
	  particles[ind1].f[j] += fac * d[j];
	  particles[ind2].f[j] -= fac * d[j];
	}
      }
    }
    
    /* real space coulomb */
    if(dist < p3m.r_cut) {
      adist = p3m.alpha * dist;
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac = p3m.bjerrum * particles[ind1].q * particles[ind2].q  * 
	exp(-adist*adist) * (erfc_part_ri + 2.0*p3m.alpha/1.772453851) / dist2;
      particles[ind1].f[0] += fac * d[0];
      particles[ind1].f[1] += fac * d[1];
      particles[ind1].f[2] += fac * d[2];
      particles[ind2].f[0] -= fac * d[0];
      particles[ind2].f[1] -= fac * d[1];
      particles[ind2].f[2] -= fac * d[2];
    }

    /* ramp */
    if(dist < ia_params->ramp_cut) {
      if (dist < 1e-4) {
	particles[ind1].f[0] += ia_params->ramp_force;
	particles[ind1].f[1] += 0;
	particles[ind1].f[2] += 0;
	particles[ind2].f[0] -= ia_params->ramp_force;
	particles[ind2].f[1] -= 0;
	particles[ind2].f[2] -= 0;
      }
      else {
	fac = ia_params->ramp_force/dist;
	particles[ind1].f[0] += fac * d[0];
	particles[ind1].f[1] += fac * d[1];
	particles[ind1].f[2] += fac * d[2];
	particles[ind2].f[0] -= fac * d[0];
	particles[ind2].f[1] -= fac * d[1];
	particles[ind2].f[2] -= fac * d[2];
      }
    }
    /* minimal particle distance calculation */
    if (dist < minimum_part_dist)
      minimum_part_dist = dist;
  }

  /* calculate k-space part of electrostatic interaction. */
  if(p3m.bjerrum != 0.0) P3M_calc_kspace_forces();
}

/************************************************************/

void force_exit()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_exit:\n",this_node));
}



/**************************************************/
/*******************  FORCES.C  *******************/
/**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "debug.h"
#include "thermostat.h"
/* #include "p3m_parallel.h" */
#include "communication.h"
#include "ghosts.h" 
#include "verlet.h"
#include "utils.h"
#include "interaction_data.h"
#include "grid.h"
#include "p3m.h"

double minimum_part_dist = -1;

void force_init()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d interaction types\n",
		      this_node,n_interaction_types));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_particle_types));
}

void force_calc()
{
  int i,j;
  int id1, id2, id3;
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

  ia_params = get_ia_param(0,0);

  minimum_part_dist = box_l[0] + box_l[1] + box_l[2];

  FORCE_TRACE(fprintf(stderr,"%d: interactions type:(0,0):\n",this_node));
  FORCE_TRACE(fprintf(stderr,"    LJ: cut=%f, eps=%f, off=%f\n",
		      ia_params->LJ_cut,ia_params->LJ_eps,ia_params->LJ_offset));
  FORCE_TRACE(fprintf(stderr,"    ES: cut=%f,\n",p3m.r_cut));
  FORCE_TRACE(fprintf(stderr,"    RAMP: cut=%f\n",ia_params->ramp_cut));


  /* initialize forces with thermostat forces */
  friction_thermo();
  /* initialize ghost forces with zero */
  for(i=n_particles;i<n_particles+n_ghosts;i++)
    for(j=0;j<3;j++)
      particles[i].f[j] = 0.0;

  /* calculate bonded interactions (loop local particles) */
  for(id1=0; id1<n_particles; id1++) {
    i=0;
    while(i<particles[id1].n_bonds) {
      type_num = particles[id1].bonds[i];
      switch(bonded_ia_params[type_num].type) {
      case BONDED_IA_FENE:
	id2 = local_index[particles[id1].bonds[i+1]];
	/* check particle existence */
	if(id2 == -1) {
	  fprintf(stderr,"WARNING: Bonded atoms %d and %d not on the same node\n",
		  particles[id1].identity,particles[id1].bonds[i+1] ); 
	  i+=2; break;
	}
	/* fene force */
	for(j=0;j<3;j++) d[j] = particles[id1].p[j] - particles[id2].p[j];
	dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
	fac = bonded_ia_params[type_num].p.fene.k_fene;
	fac /= (1.0 - dist2/SQR(bonded_ia_params[type_num].p.fene.r_fene));
	for(j=0;j<3;j++) {
	  particles[id1].f[j] -= fac*d[j];
	  particles[id2].f[j] += fac*d[j];
	}
	i+=2; break;
      case BONDED_IA_ANGLE:
	id2 = local_index[particles[id1].bonds[i+1]];
	id3 = local_index[particles[id1].bonds[i+2]];
	if(id2 == -1 || id3 == -1) {
	  fprintf(stderr,"WARNING: Bonded atoms %d, %d and %d not on the same node\n",
		  particles[id1].identity,particles[id1].bonds[i+1],particles[id1].bonds[i+2]); 
	  i+=3; break;
	}
	/* angle bend force */
	cosine=0.0;
	/* vector from particles[id2] to particles[id1] */
	for(j=0;j<3;j++) vec1[j] = particles[id2].p[j] - particles[id1].p[j];
	dist2 = SQR(vec1[0]) + SQR(vec1[1]) + SQR(vec1[2]);
	d1i = 1.0 / sqrt(dist2);
	for(j=0;j<3;j++) vec1[j] *= d1i;
	/* vector from particles[id3] to particles[id1] */
 	for(j=0;j<3;j++) vec2[j] = particles[id3].p[j] - particles[id1].p[j];
	dist2 = SQR(vec2[0]) + SQR(vec2[1]) + SQR(vec2[2]);
	d2i = 1.0 / sqrt(dist2);
	for(j=0;j<3;j++) vec2[j] *= d2i;
	/* scalar produvt of vec1 and vec2 */
	for(j=0;j<3;j++) cosine += vec1[j] * vec2[j];
	/* apply bend forces */
	for(j=0;j<3;j++) {
	  f1 = bonded_ia_params[type_num].p.angle.bend * (vec2[j] - cosine * vec1[j]) * d1i;
	  f2 = bonded_ia_params[type_num].p.angle.bend * (vec1[j] - cosine * vec2[j]) * d2i;
	  particles[id2].f[j] -= f1;
	  particles[id1].f[j] += (f1+f2);
	  particles[id3].f[j] -= f2;
	}
 	i+=3; break;
      default :
	fprintf(stderr,"WARNING: Bonds of atom %d unkwown\n",particles[id1].identity);
	i = particles[id1].n_bonds; 
	break;
     }
    }
  }

  /* calculate non bonded interactions (loop verlet list) */
  for(i=0;i<2*n_verletList;i=i+2) {
    id1 = verletList[i];
    id2 = verletList[i+1];
    ia_params = get_ia_param(particles[id1].type,particles[id2].type);
    for(j=0;j<3;j++)
      d[j] = particles[id1].p[j] - particles[id2].p[j];
    dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
    dist = sqrt(dist2);
    FORCE_TRACE(if(dist>0.01) fprintf(stderr,"P%d id(%d,%d) loc(%d, %d): dist %f\n",
			i,particles[id1].identity,particles[id2].identity,id1,id2,dist));
    /* lennnard jones */

    if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
      FORCE_TRACE(fprintf(stderr, "%d: LJ %d(%d) %d(%d) %f @ %f %f %f - %f %f %f\n", this_node,
	      id1, particles[id1].identity, id2, particles[id2].identity, dist,
	      particles[id1].p[0], particles[id1].p[1], particles[id1].p[2],
	      particles[id2].p[0], particles[id2].p[1], particles[id2].p[2]));
      
      r_off = dist - ia_params->LJ_offset;
      if(r_off>0.0) {
	frac2 = SQR(ia_params->LJ_sig/r_off);
	frac6 = frac2*frac2*frac2;
	fac = 48.* ia_params->LJ_eps * frac6*(frac6 - 0.5)*frac2;
	for(j=0;j<3;j++) {
	  particles[id1].f[j] += fac * d[j];
	  particles[id2].f[j] -= fac * d[j];
	}
      }
    }
    
    /* real space coulomb */
    if(dist < p3m.r_cut2) {
      adist = p3m.alpha * dist;
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac = p3m.bjerrum * particles[id1].q * particles[id2].q  * 
	exp(-adist*adist) * (erfc_part_ri + 2.0*p3m.alpha/1.772453851) / dist2;
      particles[id1].f[0] += fac * d[0];
      particles[id1].f[1] += fac * d[1];
      particles[id1].f[2] += fac * d[2];
      particles[id2].f[0] -= fac * d[0];
      particles[id2].f[1] -= fac * d[1];
      particles[id2].f[2] -= fac * d[2];
    }

    /* ramp */
    if(dist < ia_params->ramp_cut) {
      if (dist < minimum_part_dist)
	minimum_part_dist = dist;

      if (dist < 1e-4) {
	particles[id1].f[0] += ia_params->ramp_force;
	particles[id1].f[1] += 0;
	particles[id1].f[2] += 0;
	particles[id2].f[0] -= ia_params->ramp_force;
	particles[id2].f[1] -= 0;
	particles[id2].f[2] -= 0;
      }
      else {
	fac = ia_params->ramp_force/dist;
	particles[id1].f[0] += fac * d[0];
	particles[id1].f[1] += fac * d[1];
	particles[id1].f[2] += fac * d[2];
	particles[id2].f[0] -= fac * d[0];
	particles[id2].f[1] -= fac * d[1];
	particles[id2].f[2] -= fac * d[2];
      }
    }
  }

  /* calculate k-space part of electrostatic interaction. */
  if(p3m.bjerrum != 0.0) P3M_perform();
}

void force_exit()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_exit:\n",this_node));
}

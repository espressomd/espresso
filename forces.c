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
#include "cells.h"
#include "ghosts.h" 
#include "verlet.h"
#include "utils.h"
#include "interaction_data.h"

double Bjerrum;

void force_init()
{

  Bjerrum = 1.68;

  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d interaction types\n",
		      this_node,n_interaction_types));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_particle_types));
}

void force_calc()
{
  int i,j;
  int id1, id2;
  double d[3], dist2, dist;
  IA_parameters *ia_params;
  extern Particle *particles;
  extern int *verletList;
  double frac6,r_off;
  double fac, adist;
  double coul_r_cut = 20.;
  double Bjerrum = 1.68, alpha = 1.;
  double  erfc_part_ri;

  FORCE_TRACE(fprintf(stderr,"%d: force_calc: for %d (P %d,G %d)\n",this_node,n_particles+n_ghosts,n_particles,n_ghosts));

  ia_params = get_ia_param(0,0);
  FORCE_TRACE(fprintf(stderr,"%d: interactions type:(0,0):\n",this_node));
  FORCE_TRACE(fprintf(stderr,"    LJ: cut=%f, eps=%f\n",
		      ia_params->LJ_cut,ia_params->LJ_eps));
  FORCE_TRACE(fprintf(stderr,"    ES: cut=%f,\n",coul_r_cut));
  FORCE_TRACE(fprintf(stderr,"    RAMP: cut=%f\n",ia_params->ramp_cut));
  /* define 'new old' coordinates*/
  for(i=0;i<n_particles+n_ghosts;i++)
    for(j=0;j<3;j++)
	particles[i].p_old[j] = particles[i].p[j];

  /* initialize forces with thermostat forces */
  friction_thermo();

  for(i=n_particles;i<n_particles+n_ghosts;i++)
    for(j=0;j<3;j++)
      particles[i].f[j] = 0;

  /* loop verlet list */
  for(i=0;i<2*n_verletList;i=i+2) {
    id1 = verletList[i];
    id2 = verletList[i+1];
    ia_params = get_ia_param(particles[id1].type,particles[id2].type);
    for(j=0;j<3;j++) d[j] = particles[id1].p[j] - particles[id2].p[j];
    dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
    dist = sqrt(dist2);
    FORCE_TRACE(if(dist<0.01) fprintf(stderr,"P%d (%d,%d): dist %f\n",
			i,particles[id1].identity,particles[id2].identity,dist));
    /* lennnard jones */

    if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
      r_off = dist - ia_params->LJ_offset;
      if(r_off>0.0) {
	frac6 = SQR(ia_params->LJ_sig/r_off);
	frac6 *= SQR(frac6);
	fac = 24.* ia_params->LJ_eps * (2.*SQR(frac6) - frac6 + ia_params->LJ_shift)/r_off;
	for(j=0;j<3;j++) {
	  particles[id1].f[j] += fac * d[j]/dist;
	  particles[id2].f[j] -= fac * d[j]/dist;
	}
      }
    }
    
    /* real space coulomb */
    /*
    if(dist < coul_r_cut) {
      adist = alpha * dist;
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac = Bjerrum * particles[id1].q * particles[id2].q  * 
	exp(-adist*adist) * (erfc_part_ri + 2.0*alpha/1.772453851) / dist2;
      particles[id1].f[0] += fac * d[0];
      particles[id1].f[1] += fac * d[1];
      particles[id1].f[2] += fac * d[2];
      particles[id2].f[0] -= fac * d[0];
      particles[id2].f[1] -= fac * d[1];
      particles[id2].f[2] -= fac * d[2];
    }
    */

    /* ramp */
    if(dist < ia_params->ramp_cut) {
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
}

void force_exit()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_exit:\n",this_node));
}

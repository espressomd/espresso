/**************************************************/
/*******************  FORCES.C  *******************/
/**************************************************/

#include "forces.h"
#include "debug.h"
#include "global.h"

void force_init()
{

  Bjerrum = 1.68;

  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d interaction types\n",
		      this_node,n_particle_types));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_interaction_types));
}

void force_calc()
{
  int i,j;
  int id1, id2;
  double d[3], dist2, dist;
  IA_parameters *ia_params;
  extern int n_total_particles;
  extern Particle *particles;
  extern int *verletList;
  double frac6,r_off;
  double fac, adist;
  double coul_r_cut = 20.;
  double Bjerrum = 1.68, alpha = 1.;
  double  erfc_part_ri;


  for(i=0;i<n_total_particles;i++)
    for(j=0;j<3;j++)
      particles[i].f[j] = 0;
      
  for(i=0;i<max_particles;i=i+2)
    {
      id1 = verletList[i];
      id2 = verletList[i+1];
      ia_params = get_ia_param(id1,id2);
      
      d[0] = particles[id1].p[0] - particles[id2].p[0]; //d[0] -= dround(d[0]/box_l[0]) * box_l[0];
      d[1] = particles[id1].p[1] - particles[id2].p[1]; //d[1] -= dround(d[1]/box_l[1]) * box_l[1];
      d[2] = particles[id1].p[2] - particles[id2].p[2]; //d[2] -= dround(d[2]/box_l[2]) * box_l[2];
      dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
      dist = sqrt(dist2);

      if(dist < ia_params->LJ_cut + ia_params->LJ_offset)
	{
	  r_off = dist - ia_params->LJ_offset;
	  frac6 = SQR(ia_params->LJ_sig/r_off);
	  frac6 *= SQR(frac6);
	  fac = 24.* ia_params->LJ_eps * (2.*SQR(frac6) - frac6 + ia_params->LJ_shift)/r_off;
	  particles[id1].f[0] += fac * d[0]/dist;
	  particles[id1].f[1] += fac * d[1]/dist;
	  particles[id1].f[2] += fac * d[2]/dist;
	  particles[id2].f[0] -= fac * d[0]/dist;
	  particles[id2].f[1] -= fac * d[1]/dist;
	  particles[id2].f[2] -= fac * d[2]/dist;
	  
	}
      if(dist < coul_r_cut)
	{
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
    }
   

      
      
  FORCE_TRACE(fprintf(stderr,"%d: force_calc:\n",this_node));
}

void force_exit()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_exit:\n",this_node));
}

#include "SubtLj.hpp"
#include "interaction_data.hpp" // for iaparams LJ

int Bond::SubtLj::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const {
  int i;
  double fac_lj=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if(dist >= m_r)
    return 1;

  ia_params = get_ia_param(p1->p.type,p2->p.type);
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) { 
    r_off = dist - ia_params->LJ_offset;

    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac_lj = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);			  
    }
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac_lj = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
    }
  } 

  for(i=0;i<3;i++)
    force[i] = -fac_lj*dx[i];



  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac_lj %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(dist2),fac_lj));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac_lj %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(dist2),fac_lj));

  return 0;
}

int Bond::SubtLj::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const {
  double energy=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  
  if(dist >= m_r)
    return 1;
  
  ia_params = get_ia_param(p1->p.type,p2->p.type);
  
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
    r_off = dist - ia_params->LJ_offset;
    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      energy = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }   
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      energy = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
  }
  *_energy = -energy;
  return 0;
}

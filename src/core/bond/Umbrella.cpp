#include "Umbrella.hpp"
#include "debug.hpp"

int Bond::Umbrella::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
					  double force[3]) const {
  double distn;
  double fac=0.0;
  distn = dx[m_dir];
  fac = -m_k * (distn - m_r);
  force[m_dir] += fac;
      
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,
    "%d: OPT: umbrella f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,distn,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,
    "%d: OPT: umbrella f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,distn,fac));
  return 0;

}

int Bond::Umbrella::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], 
					   double *_energy) const {
  double distn;
  distn = dx[m_dir];  
  *_energy = 0.5 * m_k * Utils::sqr(distn - m_r);
  return 0;

}

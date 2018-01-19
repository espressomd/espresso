#include "Quartic.hpp"
#include "debug.hpp"

//---QUARTIC BOND---
int Bond::Quartic::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const {
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  // printf("Quartic dist2 %e, dist %e\n", dist2, dist);

  if ((m_r_cut > 0.0) &&
      (dist > m_r_cut)) 
    return 1;

  dr = dist - m_r;

  fac = (m_k0 * dr + m_k1 * dr * dr * dr)/dist;
  
  for(i=0;i<3;i++)
    force[i] = -fac*dx[i];

  //  printf("Quartic (%d-%d), dist %e, dx %e %e %e, dr %e, f %e %e %e\n", p1->p.identity, p2->p.identity, dist, dx[0], dx[1], dx[2], dr, force[0], force[1], force[2]);

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: QUARTIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: QUARTIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));

  return 0;
}

int Bond::Quartic::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((m_r_cut > 0.0) && 
      (dist > m_r_cut)) 
    return 1;
 
  double dr2 = SQR(dist -m_r);

  *_energy = 0.5*m_k0*dr2 + 0.25 * m_k1 * SQR(dr2);
  return 0;
}

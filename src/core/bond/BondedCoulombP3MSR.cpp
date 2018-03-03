#include "BondedCoulombP3MSR.hpp"

#include "config.hpp"
#include "utils.hpp"
#include "debug.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "p3m.hpp"

/** Computes the BONDED_COULOMB_P3M_SR pair force and adds this
    force to the particle forces (see \ref interaction_data.cpp). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond parameters.
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
*/
int Bond::BondedCoulombP3MSR::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3],
						     double force[3]) const
{
#ifdef P3M
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  if (dist < p3m.params.r_cut) {
    //Set to zero because p3m adds forces
    force[0] = force[1] = force[2] = 0.;

    p3m_add_pair_force(m_q1q2, dx, dist2, dist, force);
      
    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: BONDED_COULOMB_P3M_SR f = (%.3e,%.3e,%.3e) with part id=%d at dist %f\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: BONDED_COULOMB_P3M_SR f = (%.3e,%.3e,%.3e) with part id=%d at dist %f\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2));
  }   
#endif
  return 0;
  
}
int Bond::BondedCoulombP3MSR::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3],
						      double *_energy) const
{

#ifdef P3M
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  *_energy = p3m_pair_energy(m_q1q2, dist);
#endif
  return 0;
  
}

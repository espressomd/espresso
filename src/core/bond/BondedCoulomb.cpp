#include "BondedCoulomb.hpp"
#include "debug.hpp"

//---BONDED_COULOMB---
int Bond::BondedCoulomb::calc_bonded_pair_force(Particle *p1, Particle *p2,
                                                double dx[3],
                                                double force[3]) const {
#ifdef ELECTROSTATICS
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  fac = m_prefactor * p1->p.q * p2->p.q / (dist * dist2);

  for (i = 0; i < 3; i++)
    force[i] = fac * dx[i];
  ONEPART_TRACE(if (p1->p.identity == check_id) fprintf(
      stderr, "%d: OPT: BONDED_COULOMB f = (%.3e,%.3e,%.3e) with part id=%d at "
              "dist %f fac %.3e\n",
      this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2], p2->p.identity, dist2,
      fac));
  ONEPART_TRACE(if (p2->p.identity == check_id) fprintf(
      stderr, "%d: OPT: BONDED_COULOMB f = (%.3e,%.3e,%.3e) with part id=%d at "
              "dist %f fac %.3e\n",
      this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2], p1->p.identity, dist2,
      fac));
#endif
  return 0;
}

int Bond::BondedCoulomb::calc_bonded_pair_energy(Particle *p1, Particle *p2,
                                                 double dx[3],
                                                 double *_energy) const {
#ifdef ELECTROSTATICS
  double dist = sqrt(sqrlen(dx));
  *_energy = m_prefactor * p1->p.q * p2->p.q / dist;
#endif
  return 0;
}

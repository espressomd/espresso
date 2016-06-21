#include "ChargedRod.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"

namespace Constraints {

void ChargedRod::add_energy(Particle const *p, double const *folded_pos,
                            Observable_stat &energy) const {
#ifdef ELECTROSTATICS
  int i;
  double vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for (i = 0; i < 2; i++) {
    vec[i] = folded_pos[i] - center[i];
    c_dist_2 += SQR(vec[i]);
  }

  if ((coulomb.prefactor != 0.0) && (p->p.q != 0.0) && (lambda != 0.0)) {
    energy.coulomb[0] +=
        coulomb.prefactor * p->p.q * lambda *
        (-log(c_dist_2 * SQR(box_l_i[2])) + 2 * (M_LN2 - c_gamma));
  }
#endif
}

void ChargedRod::add_force(Particle *p, double const *folded_pos) const {
#ifdef ELECTROSTATICS
  double vec[2], c_dist_2;

  c_dist_2 = 0.0;
  for (int i = 0; i < 2; i++) {
    vec[i] = folded_pos[i] - center[i];
    c_dist_2 += SQR(vec[i]);
  }

  if ((coulomb.prefactor != 0.0) && (p->p.q != 0.0) && (lambda != 0.0)) {
    const double fac = 2 * coulomb.prefactor * lambda * p->p.q / c_dist_2;
    p->f.f[0] += fac * vec[0];
    p->f.f[1] += fac * vec[1];
    total_force[0] -= fac * vec[0];
    total_force[1] -= fac * vec[1];
  }
#endif
};
}

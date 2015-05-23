#include "ChargedPlate.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

namespace Constraints {

  void ChargedPlate::add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy) {
#ifdef ELECTROSTATICS
    if (coulomb.prefactor != 0.0 && p->p.q != 0.0 && sigma != 0.0)
      energy.coulomb[0] += -2*M_PI*coulomb.prefactor*sigma*p->p.q*fabs(folded_pos[2] - pos);
#endif
  }

  void ChargedPlate::add_force(Particle *p, const double *folded_pos) {
#ifdef ELECTROSTATICS
    double f;

    if (coulomb.prefactor != 0.0 && p->p.q != 0.0 && sigma != 0.0) {
      f = 2*M_PI*coulomb.prefactor*sigma*p->p.q;
      if (folded_pos[2] < pos)
        f = -f;
      p->f.f[2]  += f;
      total_force[2] -= f;
    }
#endif

  };
}

#ifndef SUBT_LJ_BOND_CLASS_H
#define SUBT_LJ_BOND_CLASS_H
#include "PairBond.hpp"

namespace Bond {
  class SubtLj : public PairBond {
  public:
    SubtLj() {m_bondtype = BondType::BONDED_IA_SUBT_LJ;}
    // Member function
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const override;
  };

}
#endif

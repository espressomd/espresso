#ifndef VIRTUAL_BOND
#define VIRTUAL_BOND
#include "Bond.hpp"

namespace Bond{

  class VirtualBond : public Bond {

  public:
    //constructor
    VirtualBond() : Bond(1) {m_bondtype = BondType::BONDED_IA_VIRTUAL_BOND;}
    //functions from bond
    int add_bonded_force(Particle *p1, int bl_id) override;
    int add_bonded_energy(Particle *p1, int bl_id, double *_energy) override;

  };

}

#endif

#ifndef MEMBRANE_COLLISION_BOND_CLASS_H
#define MEMBRANE_COLLISION_BOND_CLASS_H
#include "FourParticleBond.hpp"

namespace Bond {

  class MembraneCollision : public FourParticleBond {
    
  public:
    //constructor
    MembraneCollision(){m_bondtype = BondType::BONDED_IA_OIF_OUT_DIRECTION;}
    
    //force calculation
    int calc_bonded_four_particle_force(Particle *p1, Particle *p2, Particle *p3, 
					       Particle *p4, double force[3], double force2[3], 
					       double force3[3], double force4[3]) const override;
    //energy calculation
    int calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						Particle *p4, double *_energy) const override;

  };

}

#endif

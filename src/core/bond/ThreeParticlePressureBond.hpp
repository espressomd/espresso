#ifndef THREE_PARTICLE_PRESSURE_BOND
#define THREE_PARTICLE_PRESSURE_BOND
#include "ThreeParticleBond.hpp"

namespace Bond {

  class ThreeParticlePressureBond : public ThreeParticleBond {
  public:
    ThreeParticlePressureBond() : ThreeParticleBond() {}
    virtual ~ThreeParticlePressureBond() = default;

    //for add_three_body_bonded_stress
    int add_three_body_pressure(Particle *p1, int bl_id);

    // three body force calculation for pressure
    virtual int calc_3body_forces(Particle *p_mid, Particle *p_left,
				  Particle *p_right, double force1[3],
				  double force2[3], double force3[3]) const=0;

  };

}

#endif

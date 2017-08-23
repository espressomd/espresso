#ifndef THREE_PARTICLE_BOND_H
#define THREE_PARTICLE_BOND_H
#include "Bond.hpp"

namespace Bond {
  class ThreeParticleBond : public Bond {
  public:
    virtual ~ThreeParticleBond() = default;

    // new virtual methods for three particle bond calculation implemented
    // in concrete classes
    // force calculation
    virtual int add_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const=0;
    //energy calculation
    virtual int add_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						 double *_energy) const=0;

    // write forces to particles
    // there is a default implementation of this function
    // but it can be overwritten, because there are forces, which
    // sum up the forces in a different way
    virtual void write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const;

    // general bond calculation functions of abstract class
    // p1: particle, bl_id: id number of bond in bl.e
    // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
    int add_bonded_force(Particle *p1, int bl_id) const override;
    int add_bonded_energy(Particle *p1, int bl_id, double* _energy) const override;

    // get the number of bond partners
    int get_number_of_bond_partners() const override;
  };

}

#endif

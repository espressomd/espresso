#ifndef PAIR_BOND_H
#define PAIR_BOND_H
#include "Bond.hpp"
#include <iostream>

namespace Bond{
  class PairBond : public Bond {
  public:
    // constructor
    PairBond() : Bond(1) {}
    // virtual destructor
    virtual ~PairBond() = default;
    
    // new virtual methods for pair bond calculation implemented
    // in concrete classes
    // force calculation
    virtual int add_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			       double force[3]) const=0;
    //energy calculation
    virtual int add_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], 
				double *_energy) const=0;
    
    // function, which writes forces to particles
    // there is a default implementation of this function
    // but it can be overwritten, because there are forces, which
    // sum up the forces in a different way
    virtual void write_force_to_particle(Particle *p1, Particle *p2, double force[3]) const;


    // general bond calculation functions of abstract class
    // p1: particle, bl_id: id number of bond in bl.e
    // return value: 0: ok, 1: bond broken, 2: return from "add_bonded_force" in forces_inline.cpp
    int add_bonded_force(Particle *p1, int bl_id) override;
    int add_bonded_energy(Particle *p1, int bl_id, double *_energy) override;

    //pressure calculation
    int add_virial(Particle *p1, int bl_id) override;
    //for int local_stress_tensor_calc in pressure.hpp
    int calc_pair_force(Particle *p1, Particle *p2, int bl_id, double force[3]) override;

  };
  
}
#endif

#ifndef OVERLAP_BOND_ANGLE_BOND_CLASS_H
#define OVERLAP_BOND_ANGLE_BOND_CLASS_H
#include "ThreeParticleBond.hpp"
#include "Overlap.hpp"

namespace Bond {
  class OverlapBondAngle : public ThreeParticleBond, public Overlap {
  public:
    // constructor
    OverlapBondAngle(char* filename, OverlappedBondedInteraction type, double maxval, 
		      int noverlaps, double* para_a, double* para_b, double* para_c) : 
      Overlap{filename, type, maxval, noverlaps, para_a, para_b, para_c} 
    {m_bondtype = BondType::BONDED_IA_OVERLAPPED;}

    //force *
    int calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const override;
    //energy *
    int calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						 double *_energy) const override;

    boost::any get_bond_parameters_from_bond() const override;
  };
}
#endif

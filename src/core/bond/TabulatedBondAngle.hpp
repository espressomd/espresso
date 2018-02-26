#ifndef TABULATED_BOND_ANGLE_BOND_CLASS_H
#define TABULATED_BOND_ANGLE_BOND_CLASS_H
#include "ThreeParticlePressureBond.hpp"
#include "Tabulated.hpp"

namespace Bond {
  class TabulatedBondAngle : public ThreeParticlePressureBond, public Tabulated {
  public:

    //constructor
    TabulatedBondAngle(TabulatedPotential tab_pot) : 
      Tabulated{std::move(tab_pot), TabulatedBondedInteraction::TAB_BOND_ANGLE}
    {m_bondtype = BondType::BONDED_IA_TABULATED;}

    //force *
    int calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const override;
    //energy *
    int calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						 double *_energy) const override;
    // for pressure.hpp
    int calc_3body_forces(Particle *p_mid, Particle *p_left,
			   Particle *p_right, double force1[3], 
			   double force2[3], double force3[3]) const override;
  };
}
#endif

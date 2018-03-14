#ifndef ANGLE_COSSQUARE_BOND_CLASS_H
#define ANGLE_COSSQUARE_BOND_CLASS_H
#include "ThreeParticlePressureBond.hpp"
#include "Angle.hpp"

namespace Bond {

  class AngleCosSquare : public ThreeParticlePressureBond, public Angle {
  public:
    //using Angle::Angle;//inherit constructor from Angle class
    AngleCosSquare(double bend, double phi0) : Angle{bend, phi0} 
    {m_bondtype = BondType::BONDED_IA_ANGLE_COSSQUARE;}

    //force *
    int calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3],
					 double force2[3]) const override;
    //energy *
    int calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
					  double *_energy) const override;
    //virtual function from Angle
    int calc_3body_forces(Particle *p_mid, Particle *p_left,
			  Particle *p_right, double force1[3], 
			  double force2[3], double force3[3]) const override;
    
    boost::any get_bond_parameters_from_bond() const override;

    //variables
  private:
    const double m_cos_phi0 = cos(m_phi0);
  };
}
#endif

#ifndef ANGLE_COSINE_BOND_CLASS_H
#define ANGLE_COSINE_BOND_CLASS_H
#include "ThreeParticlePressureBond.hpp"
#include "Angle.hpp"

namespace Bond {

  class AngleCosine : public ThreeParticlePressureBond, public Angle {
  public:
    //using Angle::Angle;//inherit constructor from Angle class
    AngleCosine(double bend, double phi0) : Angle{bend, phi0} 
    {m_bondtype = BondType::BONDED_IA_ANGLE_COSINE;}

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
    double m_cos_phi0 = cos(m_phi0);
    double m_sin_phi0 = sin(m_phi0);
  };
}

#endif

#ifndef IBM_TRIEL_BOND_CLASS_H
#define IBM_TRIEL_BOND_CLASS_H
#include "ThreeParticleBond.hpp"
#include "CutoffBond.hpp"

namespace Bond {
  class IbmTriel : public ThreeParticleBond, public CutoffBond {
  public:

    //constructor
    IbmTriel(double a1_i, double a2_i, double b1_i, double b2_i, double l0_i, double lp0_i, 
	     double sinPhi0_i, double cosPhi0_i, double area0_i, double maxdist_i, 
	     double elasticLaw_i, double k1_i, double k2_i) :
      CutoffBond(maxdist_i), m_a1{a1_i}, m_a2{a2_i}, m_b1{b1_i}, m_b2{b2_i}, m_l0{l0_i},
      m_lp0{lp0_i}, m_sinPhi0{sinPhi0_i}, m_cosPhi0{cosPhi0_i},
      m_area0{area0_i}, m_maxdist{maxdist_i}, m_elasticLaw{elasticLaw_i}, m_k1{k1_i}, m_k2{k2_i}
    {m_bondtype = BondType::BONDED_IA_IBM_TRIEL;}

    //force *
    int calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3],
					 double force2[3]) const override;
    //energy *
    int calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
					  double *_energy) const override;

    //reset params *
    int ResetParams( const double k1, const double l0);

  private:

    // Internal function
    void RotateForces(const double f1_rot[2], const double f2_rot[2], double f1[3], double f2[3],
		      double v12[3], double v13[3]) const;

    // Specific stuff
    double m_a1;
    double m_a2;
    double m_b1;
    double m_b2;
    double m_l0;
    double m_lp0;
    double m_sinPhi0;
    double m_cosPhi0;
    double m_area0;
    double m_maxdist;
    double m_elasticLaw;
    // Always store two constants, for NeoHookean only k1 is used
    double m_k1;
    double m_k2;

  };
}
#endif

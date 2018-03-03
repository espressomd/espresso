#ifndef BONDED_COULOMB_P3M_SR
#define BONDED_COULOMB_P3M_SR
#include "PairBond.hpp"

namespace Bond{

  class BondedCoulombP3MSR : public PairBond{

  public:
    BondedCoulombP3MSR(double q1q2) : PairBond(), m_q1q2{q1q2}
    {m_bondtype = BondType::BONDED_IA_BONDED_COULOMB_P3M_SR;}

    // Member function
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3],
			       double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3],
				double *_energy) const override;

  private:
    double m_q1q2;
    
  };
  
}

#endif

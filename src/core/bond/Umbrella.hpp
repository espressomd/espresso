#ifndef UMBRELLA_BOND_CLASS_H
#define UMBRELLA_BOND_CLASS_H
#include "PairBond.hpp"

namespace Bond {
  class Umbrella : public PairBond {
  public:
    Umbrella(double k_i, int dir_i, double r_i) : m_k{k_i}, m_dir{dir_i}, m_r{r_i} {m_bondtype = BondType::BONDED_IA_UMBRELLA;}
    //force calculation
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			       double force[3]) const override;
    //energy calculation
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], 
				double *_energy) const override;

    double &k(){return m_k;}
    int &dir(){return m_dir;}
    double &r(){return m_r;}
    
  private:
    double m_k;
    int m_dir;
    double m_r;

  };
}
#endif

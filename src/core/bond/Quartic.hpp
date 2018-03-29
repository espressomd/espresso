#ifndef QUARTIC_BOND_CLASS_H
#define QUARTIC_BOND_CLASS_H
#include"PairBond.hpp"
/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

// quartic bond
namespace Bond {
class Quartic : public PairBond {
public: 
  Quartic(double k0_i, double k_1_i, double r_i, double r_cut_i) : m_k0{k0_i}, m_k1{k_1_i}, m_r{r_i}, m_r_cut{r_cut_i} {m_bondtype = BondType::BONDED_IA_QUARTIC;}
    // Member function
  int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3])
    const override;
  int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy)
    const override;

  double &k0(){return m_k0;}
  double &k1(){return m_k1;}
  double &r(){return m_r;}
  double &r_cut(){return m_r_cut;}
  
  //bond parameters
private:
  double m_k0;
  double m_k1;
  double m_r;
  double m_r_cut;
};
}
#endif

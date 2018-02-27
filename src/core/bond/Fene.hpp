#ifndef FENE_BOND_CLASS_H
#define FENE_BOND_CLASS_H
#include"PairCutoffBond.hpp"

/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

//FENE Bond
namespace Bond {
class Fene : public PairCutoffBond {
public:
  //constructor
  // initializer list -> faster
  Fene(double r0_input, double drmax_input, double drmax2_input, double drmax2i_input, double k_input) : PairCutoffBond(&m_cutoff), m_r0{r0_input}, m_drmax{drmax_input}, m_drmax2{drmax2_input}, m_drmax2i{drmax2i_input}, m_k{k_input} {m_bondtype = BondType::BONDED_IA_FENE; m_cutoff=m_r0 + m_drmax;}
  // Member function
  int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const override;
  int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const override;

  double m_cutoff;
  
private:
  // specific bond parameters
  double m_r0;
  double m_drmax;
  double m_drmax2;
  double m_drmax2i;
  double m_k;

};
}
#endif

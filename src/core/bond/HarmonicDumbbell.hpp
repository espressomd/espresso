#ifndef HARMONIC_DUMBBELL_BOND_CLASS_H
#define HARMONIC_DUMBBELL_BOND_CLASS_H
#include"PairBond.hpp"

/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

// Harmonic dumbbell bond
namespace Bond {
class HarmonicDumbbell : public PairBond {
public: 
  HarmonicDumbbell(double k1_i, double k_2_i, double r_i, double r_cut_i)
    : m_k1{k1_i}, m_k2{k_2_i}, m_r{r_i}, m_r_cut{r_cut_i}
  {m_bondtype = BondType::BONDED_IA_HARMONIC_DUMBBELL;}
  // Member function
  int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3])
    const override;
  int calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy)
    const override;

  boost::any get_bond_parameters_from_bond() const override;

  double &k1(){return m_k1;}
  double &k2(){return m_k2;}
  double &r(){return m_r;}
  double &r_cut(){return m_r_cut;}

  //bond parameters
private:
  double m_k1;
  double m_k2;
  double m_r;
  double m_r_cut;

};
}
#endif


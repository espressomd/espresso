#ifndef HARMONIC_DUMBBELL_BOND_CLASS_H
#define HARMONIC_DUMBBELL_BOND_CLASS_H
#include"Bond.hpp"
/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/

// Harmonic dumbbell bond
namespace Bond {
class HarmonicDumbbell : public Bond {
public: 
  HarmonicDumbbell(double k1_i, double k_2_i, double r_i, double r_cut_i) : m_k1{k1_i}, m_k2{k_2_i}, m_r{r_i}, m_r_cut{r_cut_i} {}
    // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const override;
  //bond parameters
private:
  double m_k1;
  double m_k2;
  double m_r;
  double m_r_cut;
};
}
#endif

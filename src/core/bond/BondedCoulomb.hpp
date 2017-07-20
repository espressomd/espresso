#ifndef BONDED_COULOMB_BOND_CLASS_H
#define BONDED_COULOMB_BOND_CLASS_H
#include"Bond.hpp"
/*
The definition of the concrete classes.
only possible to inherit public from abstact class!
*/
namespace Bond {
class BondedCoulomb : public Bond {
public:
  BondedCoulomb(double prefactor_i) : m_prefactor{prefactor_i} {}
    // Member function
  int add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const override;
  int add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const override;
  //bond parameters
private:
  double m_prefactor;
};
}
#endif

#ifndef TABULATED_BOND_LENGTH_BOND_CLASS_H
#define TABULATED_BOND_LENGTH_BOND_CLASS_H
#include "PairBond.hpp"
#include "Tabulated.hpp"

namespace Bond {
  class TabulatedBondLength : public PairBond, public Tabulated {
  public: 
    //constructor
    //using Tabulated::Tabulated; //inherit constructor
    TabulatedBondLength(TabulatedBondedInteraction tab_type, char* filename, double minval, 
			double maxval, int npoints, double invstepsize, double* f, double* e) : 
      Tabulated{tab_type, filename, minval, maxval, npoints, invstepsize, f, e} 
    {m_bondtype = BondType::BONDED_IA_TABULATED;}

    // Member function
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			      double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, 
			       double dx[3], double *_energy) const override;

  };
}

#endif

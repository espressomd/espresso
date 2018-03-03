#ifndef TABULATED_BOND_LENGTH_BOND_CLASS_H
#define TABULATED_BOND_LENGTH_BOND_CLASS_H
#include "PairBond.hpp"
#include "Tabulated.hpp"
#include "CutoffBond.hpp"

namespace Bond {
  class TabulatedBondLength : public PairBond, public Tabulated, public CutoffBond {
  public: 
    //constructor
    //using Tabulated::Tabulated; //inherit constructor
    TabulatedBondLength(TabulatedPotential tab_pot) : 
      Tabulated{std::move(tab_pot), TabulatedBondedInteraction::TAB_BOND_LENGTH},
      CutoffBond(tab_pot.cutoff())
    {m_bondtype = BondType::BONDED_IA_TABULATED;}
    // Member function
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			      double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, 
			       double dx[3], double *_energy) const override;

  };
}

#endif

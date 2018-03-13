#ifndef OVERLAP_BOND_LENGTH_BOND_CLASS_H
#define OVERLAP_BOND_LENGTH_BOND_CLASS_H
#include "PairBond.hpp"
#include "Overlap.hpp"
#include "CutoffBond.hpp"

namespace Bond {

  class OverlapBondLength : public PairBond, public Overlap, public CutoffBond {
  public:
    // constructor
    OverlapBondLength(char* filename, OverlappedBondedInteraction type, double maxval, 
		      int noverlaps, double* para_a, double* para_b, double* para_c) :
      Overlap{filename, type, maxval, noverlaps, para_a, para_b, para_c},
      CutoffBond(maxval)
    {m_bondtype = BondType::BONDED_IA_OVERLAPPED;}

    // Member function
    int calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], 
			      double force[3]) const override;
    int calc_bonded_pair_energy(Particle *p1, Particle *p2, 
			       double dx[3], double *_energy) const override;

    boost::any get_bond_parameters_from_bond() const override;

  };

}

#endif

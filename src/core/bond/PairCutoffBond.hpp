#ifndef PAIR_CUTOFF_BOND_H
#define PAIR_CUTOFF_BOND_H
#include"PairBond.hpp"

namespace Bond{

  class PairCutoffBond : public PairBond {

    PairCutoffBond(double *cutoff) : m_cutoff{cutoff} {}
    ~PairCutoffBond() = default;

    double get_cutoff(){return *m_cutoff;}

  private:
    double* m_cutoff;
    
  };
  
}

#endif

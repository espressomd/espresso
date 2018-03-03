#ifndef PAIR_CUTOFF_BOND_H
#define PAIR_CUTOFF_BOND_H

namespace Bond{

  class CutoffBond{

  public:
    CutoffBond(double cutoff) : m_cutoff{cutoff} {}
    virtual ~CutoffBond() = default;

    double get_cutoff() const {return m_cutoff;}

  private:
    double m_cutoff;
    
  };
  
}

#endif

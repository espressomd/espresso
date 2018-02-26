#ifndef TABULATED_BOND_CLASS_H
#define TABULATED_BOND_CLASS_H
#include "TabulatedBondedInteraction.hpp" //for TabulatedBondedInteraction
#include "TabulatedPotential.hpp" // For TabulatedPotential Class

namespace Bond {

  // template class
  class Tabulated {
  public:
    //constructor
    Tabulated(TabulatedPotential tab_pot, TabulatedBondedInteraction tab_type) : 
      m_tab_pot{std::move(tab_pot)}, m_tab_type{tab_type} {}

    virtual ~Tabulated() = default;

    //functions
    TabulatedBondedInteraction get_tab_type(){return m_tab_type;};
    
    //variables
    TabulatedPotential m_tab_pot;
    TabulatedBondedInteraction m_tab_type;
    
    
  };

}

#endif

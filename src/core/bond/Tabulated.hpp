#ifndef TABULATED_BOND_CLASS_H
#define TABULATED_BOND_CLASS_H
#include "TabulatedBondedInteraction.hpp" //for TabulatedBondedInteraction
#include "TabulatedPotential.hpp" // For TabulatedPotential Class

namespace Bond {

  // template class
  class Tabulated {
  public:
    //constructor
    Tabulated(double min, double max, std::vector<double> const &energy,
	      std::vector<double> const &force, TabulatedBondedInteraction tab_type) :
      m_tab_type{tab_type}
    {

      //initialize tab potential
      m_tab_pot.invstepsize = static_cast<double>(force.size() -1)/(max - min);
      m_tab_pot.force_tab = force;
      m_tab_pot.energy_tab = energy;

      switch(tab_type){
      case TabulatedBondedInteraction::TAB_BOND_LENGTH:
	m_tab_pot.minval = min;
	m_tab_pot.maxval = max;
	break;
      case TabulatedBondedInteraction::TAB_BOND_ANGLE:
	m_tab_pot.minval = 0.0;
	m_tab_pot.maxval = PI + ROUND_ERROR_PREC;
	break;
      case TabulatedBondedInteraction::TAB_BOND_DIHEDRAL:
	m_tab_pot.minval = 0.0;
	m_tab_pot.maxval = 2.0 * PI + ROUND_ERROR_PREC;
	break;
      case TabulatedBondedInteraction::TAB_UNKNOWN:
	m_tab_pot.minval = min;
	m_tab_pot.maxval = max;
	runtimeError("Unsupported tabulated type: unknown!");
	break;
      }
      
    }
    
    virtual ~Tabulated() = default;

    //functions
    TabulatedBondedInteraction get_tab_type(){return m_tab_type;};
    
    //variables
    TabulatedPotential m_tab_pot;
    TabulatedBondedInteraction m_tab_type;

    double &min(){return m_tab_pot.minval;}
    double &max(){return m_tab_pot.maxval;}
    std::vector<double> &energy(){return m_tab_pot.energy_tab;}
    std::vector<double> &force(){return m_tab_pot.force_tab;}
    
    
  };

}

#endif

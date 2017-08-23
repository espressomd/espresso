#ifndef TABULATED_BOND_CLASS_H
#define TABULATED_BOND_CLASS_H
#include "TabulatedBondedInteraction.hpp" //for TabulatedBondedInteraction
#include <stdlib.h> //free

namespace Bond {

  // template class
  class Tabulated {
  public:
    //constructor
    Tabulated(TabulatedBondedInteraction tab_type, char* filename, double minval, double maxval,
	      int npoints, double invstepsize, double* f, double* e) :
      m_tab_type{tab_type}, m_filename{filename}, m_minval{minval}, m_maxval{maxval}, 
      m_npoints{npoints}, m_invstepsize{invstepsize}, m_f{f}, m_e{e} {}

    virtual ~Tabulated(){
      free(m_filename);
      free(m_f);
      free(m_e);
    }

    // functions for class that inherit from Tabulated
    double bonded_tab_force_lookup(double val) const;
    double bonded_tab_energy_lookup(double val) const;

    //variables
    const TabulatedBondedInteraction m_tab_type;
    char* m_filename;
    const double m_minval;
    const double m_maxval;
    const int m_npoints;
    const double m_invstepsize;
    double* m_f;
    double* m_e;
    
    
  };

}

#endif

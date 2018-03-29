#ifndef OVERLAP_BOND_CLASS_H
#define OVERLAP_BOND_CLASS_H
#include "OverlappedBondedInteraction.hpp"
#include <stdlib.h> //free

namespace Bond {
  class Overlap {
  public:
    Overlap(char* filename, OverlappedBondedInteraction type, double maxval, int noverlaps, 
	    double* para_a, double* para_b, double* para_c) : 
      m_filename{filename}, m_type{type}, m_maxval{maxval}, m_noverlaps{noverlaps},
      m_para_a{para_a}, m_para_b{para_b}, m_para_c{para_c} {}

    virtual ~Overlap(){
      free(m_filename);
      free(m_para_a);
      free(m_para_b);
      free(m_para_c);
    }

    //variables
    char* m_filename;
    const OverlappedBondedInteraction m_type;
    double m_maxval;
    int m_noverlaps;
    double* m_para_a;
    double* m_para_b;
    double* m_para_c;

    char* filename(){return m_filename;}
    double &maxval(){return m_maxval;}
    int &noverlaps(){return m_noverlaps;}
    double &para_a(){return *m_para_a;}
    double &para_b(){return *m_para_b;}
    double &para_c(){return *m_para_c;}

  };
}

#endif

#ifndef OVERLAP_BOND_CLASS_H
#define OVERLAP_BOND_CLASS_H
#include "OverlappedBondedInteraction.hpp"
#include <stdlib.h> //free

namespace Bond {
  class Overlap {
  public:
    Overlap(std::string filename, OverlappedBondedInteraction type, double maxval, int noverlaps, 
	    std::vector<double> para_a, std::vector<double> para_b, std::vector<double> para_c) : 
      m_filename{filename}, m_type{type}, m_maxval{maxval}, m_noverlaps{noverlaps},
      m_para_a{para_a}, m_para_b{para_b}, m_para_c{para_c} {}

    virtual ~Overlap()=default;

    //variables
    std::string m_filename;
    const OverlappedBondedInteraction m_type;
    double m_maxval;
    int m_noverlaps;
    std::vector<double> m_para_a;
    std::vector<double> m_para_b;
    std::vector<double> m_para_c;

    std::string &filename(){return m_filename;}
    double &maxval(){return m_maxval;}
    int &noverlaps(){return m_noverlaps;}
    std::vector<double> &para_a(){return m_para_a;}
    std::vector<double> &para_b(){return m_para_b;}
    std::vector<double> &para_c(){return m_para_c;}

  };
}

#endif

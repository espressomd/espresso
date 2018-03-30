#ifndef ANGLE_BOND_CLASS_H
#define ANGLE_BOND_CLASS_H
#include "particle_data.hpp"

// Angle class from which AngleHarmonic, AngleCosine and AngleCossquare inherit from
// This is just a template class for them
namespace Bond {

  class Angle {
  public:
    // constructor
    Angle(double bend, double phi0) : m_bend{bend}, m_phi0{phi0} {}

    virtual ~Angle() = default;

    // variables
    double m_bend;
    double m_phi0;

    double &bend(){return m_bend;}
    double &phi0(){return m_phi0;}
    
  };

}


#endif

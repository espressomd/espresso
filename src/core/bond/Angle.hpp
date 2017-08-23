#ifndef ANGLE_BOND_CLASS_H
#define ANGLE_BOND_CLASS_H
#include "particle_data.hpp"

// Angle class from which AngleHarmonic, AngleCosine and AngleCossquare inherit from
// This is just a template class for them
// virtual function calc_angle_3body_forces can be called in pressure.hpp
// -> this has to be generalized for the other 3body force calculations!
// all three classes have the same variables and arguments in calc_angle_..._3body_forces(...)
namespace Bond {

  class Angle {
  public:
    // constructor
    Angle(double bend, double phi0) : m_bend{bend}, m_phi0{phi0} {}

    virtual ~Angle() = default;

    // pure virtual for the 3 different classes
    virtual void calc_angle_3body_forces(Particle *p_mid, Particle *p_left,
					 Particle *p_right, double force1[3], double force2[3], double force3[3]) const=0;

    // variables
    const double m_bend;
    const double m_phi0;
  };

}


#endif

#ifndef __EXT_MAG_FIELD_HPP
#define __EXT_MAG_FIELD_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ExternalMagneticField : public Constraint {
    ExternalMagneticField(double _ext_magn_field[3]) {
      ext_magn_field[0] = _ext_magn_field[0];
      ext_magn_field[1] = _ext_magn_field[1];
      ext_magn_field[2] = _ext_magn_field[2];
    }
    void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual std::string name() { return Constraint::name() + std::string("ExternalMagneticField"); }
    double ext_magn_field[3];
  };
}

#endif

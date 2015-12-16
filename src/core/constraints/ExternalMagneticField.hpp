#ifndef __EXT_MAG_FIELD_HPP
#define __EXT_MAG_FIELD_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ExternalMagneticField : public Constraint {
    ExternalMagneticField() = default;
    ExternalMagneticField(double _ext_magn_field[3]) : ext_magn_field(Vector3d(_ext_magn_field)) {}
    
    void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual const std::string name() { return Constraint::name() + std::string("ExternalMagneticField"); }
    Vector3d ext_magn_field;

    /** Parsing stuff */
    Parameters get_parameters();
    Parameters all_parameters() const;
    void set_parameter(const std::string &name, const Variant &value);    
  };
}

#endif

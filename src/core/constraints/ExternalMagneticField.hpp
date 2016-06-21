#ifndef __EXT_MAG_FIELD_HPP
#define __EXT_MAG_FIELD_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ExternalMagneticField : public Constraint {
    ExternalMagneticField() = default;
    ExternalMagneticField(Utils::Vector3d _ext_magn_field) : ext_magn_field(_ext_magn_field) {}
    
    void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual const std::string name() { return Constraint::name() + std::string("ExternalMagneticField"); }
    Utils::Vector3d ext_magn_field;

    /** Parsing stuff */
    std::map<std::string, Variant> get_parameters();
    Parameters all_parameters() const;
    void set_parameter(const std::string &name, const Variant &value);    
  };
}

#endif

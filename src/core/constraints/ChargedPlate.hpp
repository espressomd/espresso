#ifndef __CHARGED_PLATE_HPP
#define __CHARGED_PLATE_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ChargedPlate : public Constraint {
    /** Constraint stuff */
    ChargedPlate() : pos(0.0), sigma(0.0) {}
    ChargedPlate(double _pos, double _sigma) : pos(_pos), sigma(_sigma) { }
    void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual const std::string name() const { return Constraint::name() + std::string("ChargedPlate"); }
    double pos;
    double sigma;

    /** Parsing stuff */
    std::map<std::string, Variant> get_parameters();
    Parameters all_parameters() const;
    void set_parameter(const std::string &name, const Variant &value);
  };
}

#endif

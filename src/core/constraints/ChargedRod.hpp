#ifndef __CHARGED_ROD_HPP
#define __CHARGED_ROD_HPP

#include "Constraint.hpp"
#include "Vector.hpp"

namespace Constraints {
  struct ChargedRod : public Constraint {
    ChargedRod() : ChargedRod(0., 0., 0.) {}
    ChargedRod(double center_x, double center_y, double _lambda) : lambda(_lambda) {
      center[0] = center_x;
      center[1] = center_y;
    }
    /** Constraint stuff */
    void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual const std::string name() const { return Constraint::name() + std::string("ChargedRod"); }
    Utils::Vector3d center;
    double lambda;

   /** Parsing stuff */
    std::map<std::string, Variant> get_parameters();
    Parameters all_parameters() const;
    void set_parameter(const std::string &name, const Variant &value);
  };
}

#endif

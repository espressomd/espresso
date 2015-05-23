#ifndef __CHARGED_ROD_HPP
#define __CHARGED_ROD_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ChargedPlate : public Constraint {
    ChargedPlate(double _pos, double _sigma) : pos(_pos), sigma(_sigma) { }
    void add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual std::string name() { return Constraint::name() + std::string("ChargedPlate"); }
    double pos;
    double sigma;
  };
}

#endif

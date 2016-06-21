#ifndef CONSTRAINTS_CHARGED_ROD_HPP
#define CONSTRAINTS_CHARGED_ROD_HPP

#include <vector>

#include "Constraint.hpp"

namespace Constraints {
struct ChargedRod : public Constraint {
  void add_energy(Particle const *p, double const *folded_pos,
                  Observable_stat &energy) const override;
  void add_force(Particle *p, double const *folded_pos) const override;

  Vector2d center;
  double lambda;

private:
  static constexpr double c_gamma{0.57721566490153286060651209008};
};
}

#endif

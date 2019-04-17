/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CONSTRAINTS_EXTERNAL_FIELD_HPP
#define CONSTRAINTS_EXTERNAL_FIELD_HPP

#include "Constraint.hpp"
#include "field_coupling/ForceField.hpp"

namespace Constraints {
/**
 * @brief Constraint interface for ExternalField::ForceField.
 */
template <typename Coupling, typename Field>
class ExternalField : public Constraint {
  FieldCoupling::ForceField<Coupling, Field> impl;

public:
  template <typename... Args>
  ExternalField(Args &&... args) : impl(std::forward<Args>(args)...) {}

  const Coupling &coupling() const { return impl.coupling(); }
  const Field &field() const { return impl.field(); }

  void add_energy(const Particle &, const Utils::Vector3d &, double t,
                  Observable_stat &) const override {}

  ParticleForce force(const Particle &p, const Utils::Vector3d &folded_pos,
                      double t) override {
    return impl.force(p, folded_pos, t);
  }

  bool fits_in_box(Utils::Vector3d const &box) const override {
    return impl.field().fits_in_box(box);
  }
};
} /* namespace Constraints */

#endif

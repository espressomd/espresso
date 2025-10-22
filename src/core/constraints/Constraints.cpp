/*
 * Copyright (C) 2010-2024 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Constraints.hpp"
#include "BoxGeometry.hpp"
#include "Constraint.hpp"
#include "Observable_stat.hpp"
#include "system/System.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Constraints {

void Constraints::add(std::shared_ptr<Constraint> const &constraint) {
  auto &system = get_system();
  if (not constraint->fits_in_box(system.box_geo->length())) {
    throw std::runtime_error("Constraint not compatible with box size.");
  }
  assert(not contains(constraint));
  m_constraints.emplace_back(constraint);
  system.on_constraint_change();
}
void Constraints::remove(std::shared_ptr<Constraint> const &constraint) {
  auto &system = get_system();
  assert(contains(constraint));
  std::erase(m_constraints, constraint);
  system.on_constraint_change();
}
void Constraints::add_forces(ParticleRange &particles, double time) const {
  if (m_constraints.empty())
    return;

  reset_forces();
  auto const &box_geo = *get_system().box_geo;

  for (auto &p : particles) {
    auto const pos = box_geo.folded_position(p.pos());
    ParticleForce force{};
    for (auto const &constraint : *this) {
      force += constraint->force(p, pos, time);
    }

    p.force_and_torque() += force;
  }
}

void Constraints::add_energy(ParticleRange const &particles, double time,
                             Observable_stat &obs_energy) const {
  if (m_constraints.empty())
    return;

  auto const &box_geo = *get_system().box_geo;

  for (auto const &p : particles) {
    auto const pos = box_geo.folded_position(p.pos());

    for (auto const &constraint : *this) {
      constraint->add_energy(p, pos, time, obs_energy);
    }
  }
}

} // namespace Constraints

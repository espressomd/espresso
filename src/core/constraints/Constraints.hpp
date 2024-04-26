/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef CORE_CONSTRAINTS_CONSTRAINTS_HPP
#define CORE_CONSTRAINTS_CONSTRAINTS_HPP

#include "BoxGeometry.hpp"
#include "Observable_stat.hpp"
#include "system/System.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Constraints {
template <class ParticleRange, class Constraint> class Constraints {
  using container_type = std::vector<std::shared_ptr<Constraint>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  void reset_forces() const {
    for (auto const &constraint : *this) {
      constraint->reset_force();
    }
  }

  container_type m_constraints;

public:
  bool contains(std::shared_ptr<Constraint> const &constraint) const noexcept {
    return std::find(begin(), end(), constraint) != end();
  }
  void add(std::shared_ptr<Constraint> const &constraint) {
    auto &system = System::get_system();
    auto const &box_geo = *system.box_geo;
    if (not constraint->fits_in_box(box_geo.length())) {
      throw std::runtime_error("Constraint not compatible with box size.");
    }
    assert(not contains(constraint));
    m_constraints.emplace_back(constraint);
    system.on_constraint_change();
  }
  void remove(std::shared_ptr<Constraint> const &constraint) {
    auto &system = System::get_system();
    assert(contains(constraint));
    m_constraints.erase(std::remove(begin(), end(), constraint), end());
    system.on_constraint_change();
  }

  iterator begin() { return m_constraints.begin(); }
  iterator end() { return m_constraints.end(); }
  const_iterator begin() const { return m_constraints.begin(); }
  const_iterator end() const { return m_constraints.end(); }

  void add_forces(BoxGeometry const &box_geo, ParticleRange &particles,
                  double time) const {
    if (m_constraints.empty())
      return;

    reset_forces();

    for (auto &p : particles) {
      auto const pos = box_geo.folded_position(p.pos());
      ParticleForce force{};
      for (auto const &constraint : *this) {
        force += constraint->force(p, pos, time);
      }

      p.force_and_torque() += force;
    }
  }

  void add_energy(BoxGeometry const &box_geo, ParticleRange const &particles,
                  double time, Observable_stat &obs_energy) const {
    if (m_constraints.empty())
      return;

    for (auto const &p : particles) {
      auto const pos = box_geo.folded_position(p.pos());

      for (auto const &constraint : *this) {
        constraint->add_energy(p, pos, time, obs_energy);
      }
    }
  }

  void veto_boxl_change() const {
    if (not m_constraints.empty()) {
      throw std::runtime_error("The box size can not be changed because there "
                               "are active constraints.");
    }
  }

  void on_boxl_change() const { veto_boxl_change(); }
};
} // namespace Constraints

#endif

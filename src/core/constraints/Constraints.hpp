/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "grid.hpp"
#include "statistics.hpp"

#include <memory>
#include <vector>

void on_constraint_change();

namespace Constraints {
template <class ParticleRange, class Constraint> class Constraints {
  using container_type = std::vector<std::shared_ptr<Constraint>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  void reset_forces() const {
    for (auto const &c : *this) {
      c->reset_force();
    }
  }

  container_type m_constraints;

public:
  void add(std::shared_ptr<Constraint> const &c) {
    if (not c->fits_in_box(box_geo.length())) {
      throw std::runtime_error("Constraint not compatible with box size.");
    }

    m_constraints.emplace_back(c);
    on_constraint_change();
  }
  void remove(std::shared_ptr<Constraint> const &c) {
    m_constraints.erase(
        std::remove(m_constraints.begin(), m_constraints.end(), c),
        m_constraints.end());
    on_constraint_change();
  }

  iterator begin() { return m_constraints.begin(); }
  iterator end() { return m_constraints.end(); }
  const_iterator begin() const { return m_constraints.begin(); }
  const_iterator end() const { return m_constraints.end(); }

  void add_forces(ParticleRange &particles, double t) const {
    if (m_constraints.empty())
      return;

    reset_forces();

    for (auto &p : particles) {
      auto const pos = folded_position(p.r.p, box_geo);
      ParticleForce force{};
      for (auto const &c : *this) {
        force += c->force(p, pos, t);
      }

      p.f += force;
    }
  }

  void add_energy(const ParticleRange &particles, double t,
                  Observable_stat &energy) const {
    for (auto &p : particles) {
      auto const pos = folded_position(p.r.p, box_geo);

      for (auto const &c : *this) {
        c->add_energy(p, pos, t, energy);
      }
    }
  }

  void on_boxl_change() const {
    if (not this->empty()) {
      throw std::runtime_error("The box size can not be changed because there "
                               "are active constraints.");
    }
  }
};
} // namespace Constraints

#endif

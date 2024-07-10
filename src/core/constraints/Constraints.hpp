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

#pragma once

#include "Constraint.hpp"
#include "system/Leaf.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

class ParticleRange;
class Observable_stat;

namespace Constraints {
class Constraints : public System::Leaf<Constraints> {
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
  void add(std::shared_ptr<Constraint> const &constraint);
  void remove(std::shared_ptr<Constraint> const &constraint);

  iterator begin() { return m_constraints.begin(); }
  iterator end() { return m_constraints.end(); }
  const_iterator begin() const { return m_constraints.begin(); }
  const_iterator end() const { return m_constraints.end(); }

  void add_forces(ParticleRange &particles, double time) const;

  void add_energy(ParticleRange const &particles, double time,
                  Observable_stat &obs_energy) const;

  void veto_boxl_change() const {
    if (not m_constraints.empty()) {
      throw std::runtime_error("The box size can not be changed because there "
                               "are active constraints.");
    }
  }

  void on_boxl_change() const { veto_boxl_change(); }
};
} // namespace Constraints

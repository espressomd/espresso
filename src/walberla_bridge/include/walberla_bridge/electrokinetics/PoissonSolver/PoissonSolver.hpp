/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include <walberla_bridge/LatticeWalberla.hpp>

#include <cstddef>
#include <memory>
#include <utility>

namespace walberla {

class PoissonSolver {
private:
  std::shared_ptr<LatticeWalberla> m_lattice;
  double m_permittivity;

public:
  PoissonSolver(std::shared_ptr<LatticeWalberla> lattice, double permittivity)
      : m_lattice(std::move(lattice)), m_permittivity(permittivity) {}
  virtual ~PoissonSolver() = default;

  virtual void reset_charge_field() = 0;

  virtual void add_charge_to_field(std::size_t id, double valency,
                                   bool is_double_precision) = 0;

  [[nodiscard]] virtual std::size_t get_potential_field_id() const noexcept = 0;

  void set_permittivity(double permittivity) noexcept {
    m_permittivity = permittivity;
  }

  [[nodiscard]] double get_permittivity() const noexcept {
    return m_permittivity;
  }

  [[nodiscard]] auto const &get_lattice() const noexcept { return *m_lattice; }

  virtual void solve() = 0;
};

} // namespace walberla

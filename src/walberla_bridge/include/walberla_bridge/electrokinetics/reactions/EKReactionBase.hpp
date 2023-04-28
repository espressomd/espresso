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

#include "EKReactant.hpp"
#include <walberla_bridge/LatticeWalberla.hpp>

#include <memory>
#include <utility>
#include <vector>

namespace walberla {

class EKReactionBase {
private:
  std::vector<std::shared_ptr<EKReactant>> m_reactants;
  double m_coefficient;

  std::shared_ptr<LatticeWalberla> m_lattice;

public:
  EKReactionBase(std::shared_ptr<LatticeWalberla> lattice,
                 std::vector<std::shared_ptr<EKReactant>> reactants,
                 double coefficient)
      : m_reactants(std::move(reactants)), m_coefficient(coefficient),
        m_lattice(std::move(lattice)) {}

  virtual ~EKReactionBase() = default;

  void set_coefficient(double coefficient) noexcept {
    m_coefficient = coefficient;
  }
  [[nodiscard]] double get_coefficient() const noexcept {
    return m_coefficient;
  }
  [[nodiscard]] auto get_lattice() const noexcept { return m_lattice; }
  [[nodiscard]] auto get_reactants() const noexcept { return m_reactants; }

  virtual void perform_reaction() = 0;
};

} // namespace walberla

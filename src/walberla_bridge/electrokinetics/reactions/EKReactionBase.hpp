/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_WALBERLA_BRIDGE_EKREACTIONBASE_HPP
#define ESPRESSO_WALBERLA_BRIDGE_EKREACTIONBASE_HPP

#include "EKReactant.hpp"
#include "LatticeWalberla.hpp"

#include <memory>

namespace walberla {
template <typename FloatType> class EKReactionBase {
private:
  std::vector<std::shared_ptr<EKReactant<FloatType>>> m_reactants;
  FloatType m_coefficient;

  std::shared_ptr<LatticeWalberla> m_lattice;

public:
  EKReactionBase(std::shared_ptr<LatticeWalberla> lattice,
                 std::vector<std::shared_ptr<EKReactant<FloatType>>> reactants,
                 FloatType coefficient)
      : m_reactants(std::move(reactants)), m_coefficient(coefficient),
        m_lattice(std::move(lattice)) {}

  void set_coefficient(FloatType coefficient) noexcept {
    m_coefficient = coefficient;
  }
  [[nodiscard]] FloatType get_coefficient() const noexcept {
    return m_coefficient;
  }
  [[nodiscard]] auto get_lattice() const noexcept { return m_lattice; }
  [[nodiscard]] auto get_reactants() const noexcept { return m_reactants; }

  virtual void perform_reaction() = 0;
};
} // namespace walberla

#endif // ESPRESSO_WALBERLA_BRIDGE_EKREACTIONBASE_HPP

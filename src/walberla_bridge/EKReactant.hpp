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

#ifndef ESPRESSO_SRC_WALBERLA_BRIDGE_EKREACTANT_HPP
#define ESPRESSO_SRC_WALBERLA_BRIDGE_EKREACTANT_HPP

#include "EKinWalberlaBase.hpp"

namespace walberla {

template <typename FloatType> class EKReactant {
private:
  std::shared_ptr<EKinWalberlaBase<FloatType>> m_ekspecies;
  FloatType m_stoech_coeff;
  FloatType m_order;

public:
  EKReactant(std::shared_ptr<EKinWalberlaBase<FloatType>> ekspecies,
             FloatType stoech_coeff, FloatType order)
      : m_ekspecies(std::move(ekspecies)), m_stoech_coeff(stoech_coeff),
        m_order(order) {}

  void set_stoech_coefficient(FloatType stoech_coeff) {
    m_stoech_coeff = stoech_coeff;
  }

  [[nodiscard]] FloatType get_stoech_coeff() const { return m_stoech_coeff; }

  void set_order(FloatType order) { m_order = order; }

  [[nodiscard]] FloatType get_order() const { return m_order; }

  void set_species(std::shared_ptr<EKinWalberlaBase<FloatType>> ekspecies) {
    m_ekspecies = std::move(ekspecies);
  }

  [[nodiscard]] auto get_species() const { return m_ekspecies; }
};
} // namespace walberla

#endif // ESPRESSO_SRC_WALBERLA_BRIDGE_EKREACTANT_HPP

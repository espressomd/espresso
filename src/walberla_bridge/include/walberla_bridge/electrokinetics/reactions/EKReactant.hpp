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

#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include <memory>
#include <utility>

namespace walberla {

class EKReactant {
private:
  std::shared_ptr<EKinWalberlaBase> m_ekspecies;
  double m_stoech_coeff;
  double m_order;

public:
  EKReactant(std::shared_ptr<EKinWalberlaBase> ekspecies, double stoech_coeff,
             double order)
      : m_ekspecies(std::move(ekspecies)), m_stoech_coeff(stoech_coeff),
        m_order(order) {}

  void set_stoech_coefficient(double stoech_coeff) noexcept {
    m_stoech_coeff = stoech_coeff;
  }

  [[nodiscard]] double get_stoech_coeff() const noexcept {
    return m_stoech_coeff;
  }

  void set_order(double order) noexcept { m_order = order; }

  [[nodiscard]] double get_order() const noexcept { return m_order; }

  void set_species(std::shared_ptr<EKinWalberlaBase> ekspecies) noexcept {
    m_ekspecies = std::move(ekspecies);
  }

  [[nodiscard]] auto get_species() const noexcept { return m_ekspecies; }
};

} // namespace walberla

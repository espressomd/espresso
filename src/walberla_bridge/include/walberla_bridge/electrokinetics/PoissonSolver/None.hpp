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

#include "PoissonSolver.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>

#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>

#include <cstddef>
#include <memory>
#include <utility>

namespace walberla {

template <typename FloatType> class None : public PoissonSolver {
private:
  BlockDataID m_potential_field_id;

  using PotentialField = GhostLayerField<FloatType, 1>;

public:
  explicit None(std::shared_ptr<LatticeWalberla> lattice)
      : PoissonSolver(std::move(lattice), 0.0) {
    m_potential_field_id = field::addToStorage<PotentialField>(
        get_lattice().get_blocks(), "potential field", 0.0, field::fzyx,
        get_lattice().get_ghost_layers());
  }
  ~None() override = default;

  void reset_charge_field() override {}
  void add_charge_to_field(std::size_t, double, bool) override {}

  [[nodiscard]] std::size_t get_potential_field_id() const noexcept override {
    return m_potential_field_id;
  }

  void solve() override {}
};

} // namespace walberla

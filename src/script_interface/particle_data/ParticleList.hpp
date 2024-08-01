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

#pragma once

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"
#include "script_interface/system/Leaf.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Particles {

class ParticleList : public System::Leaf {
  std::weak_ptr<CellSystem::CellSystem> m_cell_structure;
  std::weak_ptr<Interactions::BondedInteractions> m_bonded_ias;

  auto get_cell_structure() {
    auto ptr = m_cell_structure.lock();
    assert(ptr != nullptr);
    return ptr;
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  void do_construct(VariantMap const &) override {}

  void attach(std::weak_ptr<CellSystem::CellSystem> cell_structure,
              std::weak_ptr<Interactions::BondedInteractions> bonded_ias) {
    m_cell_structure = cell_structure;
    m_bonded_ias = bonded_ias;
  }
};

} // namespace Particles
} // namespace ScriptInterface

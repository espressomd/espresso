/*
 * Copyright (C) 2021-2026 The ESPResSo project
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
#include "script_interface/particle_data/ParticleHandle.hpp"
#include "script_interface/system/System.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace ReactionMethods {

class ReactionAlgorithm : public AutoParameters<ReactionAlgorithm> {
  std::shared_ptr<System::System> m_system;
  std::shared_ptr<CellSystem::CellSystem> m_cell_system;

public:
  void do_construct(VariantMap const &params) override;
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

protected:
  std::shared_ptr<Particles::ParticleModifier> m_particle_modifier;

private:
  std::string get_internal_state() const override {
    throw std::runtime_error("Reaction methods do not support checkpointing");
  }
};

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

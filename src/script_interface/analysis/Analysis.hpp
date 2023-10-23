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

#include "ObservableStat.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/system/Leaf.hpp"

#include <string>

namespace ScriptInterface {
namespace Analysis {

class Analysis : public System::Leaf {
  std::shared_ptr<ObservableStat> m_obs_stat;

  /** @brief Check if a particle type exists. */
  void check_particle_type(int p_type) const;

  void do_construct(VariantMap const &params) override {
    m_obs_stat = std::make_shared<ObservableStat>();
  }

  void on_bind_system(::System::System &) override {
    m_obs_stat->bind_system(m_system.lock());
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &parameters) override;
};

} // namespace Analysis
} // namespace ScriptInterface

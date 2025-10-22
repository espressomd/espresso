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
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/system/Leaf.hpp"

#include "Integrator.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

class IntegratorHandle : public AutoParameters<IntegratorHandle, System::Leaf> {
  std::shared_ptr<Integrator> m_instance;
  std::shared_ptr<VariantMap> m_params;

  void on_bind_system(::System::System &system) override;

public:
  IntegratorHandle();

  void do_construct(VariantMap const &params) override {
    m_params = std::make_shared<VariantMap>(params);
  }
};

} // namespace Integrators
} // namespace ScriptInterface

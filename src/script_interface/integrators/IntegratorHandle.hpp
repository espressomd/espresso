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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_INTEGRATORS_INTEGRATOR_HANDLE_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_INTEGRATORS_INTEGRATOR_HANDLE_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "Integrator.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

class IntegratorHandle : public AutoParameters<IntegratorHandle> {
  std::shared_ptr<Integrator> m_instance;

public:
  IntegratorHandle();

  void do_construct(VariantMap const &params) override;
};

} // namespace Integrators
} // namespace ScriptInterface

#endif

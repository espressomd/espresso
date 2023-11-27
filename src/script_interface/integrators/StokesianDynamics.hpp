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

#include "config/config.hpp"

#ifdef STOKESIAN_DYNAMICS

#include "Integrator.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "core/stokesian_dynamics/sd_interface.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

class StokesianDynamics : public AutoParameters<StokesianDynamics, Integrator> {
  std::shared_ptr<::StokesianDynamicsParameters> m_instance;

public:
  StokesianDynamics();

  void do_construct(VariantMap const &params) override;
  void activate() override;

  ::StokesianDynamicsParameters const &get_instance() const {
    return *m_instance;
  }
};

} // namespace Integrators
} // namespace ScriptInterface

#endif // STOKESIAN_DYNAMICS

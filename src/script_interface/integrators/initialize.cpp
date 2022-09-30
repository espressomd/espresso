/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

#include "initialize.hpp"

#include "BrownianDynamics.hpp"
#include "IntegratorHandle.hpp"
#include "SteepestDescent.hpp"
#include "StokesianDynamics.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletIsoNPT.hpp"
#include "config/config.hpp"

namespace ScriptInterface {
namespace Integrators {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<IntegratorHandle>("Integrators::IntegratorHandle");
  om->register_new<BrownianDynamics>("Integrators::BrownianDynamics");
  om->register_new<SteepestDescent>("Integrators::SteepestDescent");
#ifdef STOKESIAN_DYNAMICS
  om->register_new<StokesianDynamics>("Integrators::StokesianDynamics");
#endif // STOKESIAN_DYNAMICS
  om->register_new<VelocityVerlet>("Integrators::VelocityVerlet");
#ifdef NPT
  om->register_new<VelocityVerletIsoNPT>("Integrators::VelocityVerletIsoNPT");
#endif // NPT
}

} // namespace Integrators
} // namespace ScriptInterface

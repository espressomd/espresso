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

#include "IntegratorHandle.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "BrownianDynamics.hpp"
#include "SteepestDescent.hpp"
#include "StokesianDynamics.hpp"
#include "VelocityVerlet.hpp"
#include "VelocityVerletIsoNPT.hpp"

#include "core/PropagationMode.hpp"
#include "core/integrators/Propagation.hpp"
#include "core/system/System.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

IntegratorHandle::IntegratorHandle() {
  add_parameters({
      {"time_step",
       [this](Variant const &v) {
         context()->parallel_try_catch(
             [&]() { get_system().set_time_step(get_value<double>(v)); });
       },
       [this]() { return get_system().get_time_step(); }},
      {"time",
       [&](Variant const &v) {
         get_system().set_sim_time(get_value<double>(v));
       },
       [this]() { return get_system().get_sim_time(); }},
      {"force_cap",
       [this](Variant const &v) {
         get_system().set_force_cap(get_value<double>(v));
       },
       [this]() { return get_system().get_force_cap(); }},
      {"integrator",
       [this](Variant const &v) {
         auto const old_instance = m_instance;
         m_instance = get_value<std::shared_ptr<Integrator>>(v);
         if (old_instance) {
           old_instance->deactivate();
         }
         m_instance->bind_system(m_system.lock());
         m_instance->activate();
       },
       [this]() {
         switch (get_system().propagation->integ_switch) {
         case INTEG_METHOD_STEEPEST_DESCENT:
           return Variant{
               std::dynamic_pointer_cast<SteepestDescent>(m_instance)};
#ifdef NPT
         case INTEG_METHOD_NPT_ISO:
           return Variant{
               std::dynamic_pointer_cast<VelocityVerletIsoNPT>(m_instance)};
#endif
         case INTEG_METHOD_BD:
           return Variant{
               std::dynamic_pointer_cast<BrownianDynamics>(m_instance)};
#ifdef STOKESIAN_DYNAMICS
         case INTEG_METHOD_SD:
           return Variant{
               std::dynamic_pointer_cast<StokesianDynamics>(m_instance)};
#endif // STOKESIAN_DYNAMICS
         default: {
           auto ptr = std::dynamic_pointer_cast<VelocityVerlet>(m_instance);
           assert(ptr.get());
           return Variant{ptr};
         }
         }
       }},
  });
}

void IntegratorHandle::on_bind_system(::System::System &system) {
  auto const &params = *m_params;
  for (auto const &key : get_parameter_insertion_order()) {
    if (params.count(key) != 0ul) {
      if (not(key == "time_step" and
              system.propagation->integ_switch == INTEG_METHOD_NVT and
              system.get_time_step() == -1. and
              is_type<double>(params.at(key)) and
              get_value<double>(is_type<double>(params.at(key))) == -1.)) {
        do_set_parameter(key, params.at(key));
      }
    }
  }
  m_params.reset();
}

} // namespace Integrators
} // namespace ScriptInterface

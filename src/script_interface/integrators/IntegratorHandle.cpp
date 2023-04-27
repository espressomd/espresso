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

#include "core/forcecap.hpp"
#include "core/integrate.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Integrators {

IntegratorHandle::IntegratorHandle() {
  add_parameters({
      {"integrator",
       [this](Variant const &v) {
         m_instance = get_value<std::shared_ptr<Integrator>>(v);
         m_instance->activate();
       },
       [this]() {
         switch (::integ_switch) {
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
      {"time_step",
       [this](Variant const &v) {
         context()->parallel_try_catch(
             [&]() { set_time_step(get_value<double>(v)); });
       },
       []() { return get_time_step(); }},
      {"time", [](Variant const &v) { set_time(get_value<double>(v)); },
       []() { return get_sim_time(); }},
      {"force_cap",
       [](Variant const &v) { set_force_cap(get_value<double>(v)); },
       []() { return get_force_cap(); }},
  });
}

static bool checkpoint_filter(typename VariantMap::value_type const &kv) {
  /* When loading from a checkpoint file, defer the integrator object to last,
   * and skip the time_step if it is -1 to avoid triggering sanity checks.
   */
  return kv.first == "integrator" or
         (kv.first == "time_step" and ::integ_switch == INTEG_METHOD_NVT and
          get_time_step() == -1. and is_type<double>(kv.second) and
          get_value<double>(kv.second) == -1.);
}

void IntegratorHandle::do_construct(VariantMap const &params) {
  for (auto const &kv : params) {
    if (not checkpoint_filter(kv)) {
      do_set_parameter(kv.first, kv.second);
    }
  }
  do_set_parameter("integrator", params.at("integrator"));
}

} // namespace Integrators
} // namespace ScriptInterface

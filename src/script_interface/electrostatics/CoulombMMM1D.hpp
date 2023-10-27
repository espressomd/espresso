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

#ifdef ELECTROSTATICS

#include "Actor.hpp"

#include "core/electrostatics/mmm1d.hpp"

#include "script_interface/get_value.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Coulomb {

class CoulombMMM1D : public Actor<CoulombMMM1D, ::CoulombMMM1D> {

public:
  CoulombMMM1D() {
    add_parameters({
        {"is_tuned", AutoParameter::read_only,
         [this]() { return actor()->is_tuned(); }},
        {"far_switch_radius", AutoParameter::read_only,
         [this]() { return actor()->far_switch_radius; }},
        {"maxPWerror", AutoParameter::read_only,
         [this]() { return actor()->maxPWerror; }},
        {"timings", AutoParameter::read_only,
         [this]() { return actor()->tune_timings; }},
        {"verbose", AutoParameter::read_only,
         [this]() { return actor()->tune_verbose; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([this, &params]() {
      m_actor = std::make_shared<CoreActorClass>(
          get_value<double>(params, "prefactor"),
          get_value<double>(params, "maxPWerror"),
          get_value<double>(params, "far_switch_radius"),
          get_value<int>(params, "timings"),
          get_value<bool>(params, "verbose"));
    });
    set_charge_neutrality_tolerance(params);
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ELECTROSTATICS

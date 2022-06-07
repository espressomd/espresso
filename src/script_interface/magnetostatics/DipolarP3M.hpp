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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_MAGNETOSTATICS_DIPOLAR_P3M_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_MAGNETOSTATICS_DIPOLAR_P3M_HPP

#include "config.hpp"

#ifdef DP3M

#include "Actor.hpp"

#include "core/magnetostatics/dp3m.hpp"

#include "script_interface/get_value.hpp"

#include <memory>

namespace ScriptInterface {
namespace Dipoles {

class DipolarP3M : public Actor<DipolarP3M, ::DipolarP3M> {
  bool m_tune;

public:
  DipolarP3M() {
    add_parameters({
        {"alpha_L", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.alpha_L; }},
        {"r_cut_iL", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.r_cut_iL; }},
        {"mesh", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.mesh; }},
        {"mesh_off", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.mesh_off; }},
        {"cao", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.cao; }},
        {"accuracy", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.accuracy; }},
        {"epsilon", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.epsilon; }},
        {"a", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.a; }},
        {"alpha", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.alpha; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return actor()->dp3m.params.r_cut; }},
        {"is_tuned", AutoParameter::read_only,
         [this]() { return actor()->is_tuned(); }},
        {"verbose", AutoParameter::read_only,
         [this]() { return actor()->tune_verbose; }},
        {"timings", AutoParameter::read_only,
         [this]() { return actor()->tune_timings; }},
        {"tune", AutoParameter::read_only, [this]() { return m_tune; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_tune = get_value<bool>(params, "tune");
    context()->parallel_try_catch([&]() {
      auto p3m = P3MParameters{get_value<bool>(params, "tune"),
                               get_value<double>(params, "epsilon"),
                               get_value<double>(params, "r_cut"),
                               get_value<Utils::Vector3i>(params, "mesh"),
                               get_value<Utils::Vector3d>(params, "mesh_off"),
                               get_value<int>(params, "cao"),
                               get_value<double>(params, "alpha"),
                               get_value<double>(params, "accuracy")};
      m_actor = std::make_shared<CoreActorClass>(
          std::move(p3m), get_value<double>(params, "prefactor"),
          get_value<int>(params, "timings"),
          get_value<bool>(params, "verbose"));
    });
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DP3M
#endif

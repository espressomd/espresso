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

#ifdef DP3M

#include "Actor.hpp"

#include "core/magnetostatics/dp3m.hpp"
#include "core/magnetostatics/dp3m.impl.hpp"
#include "core/p3m/FFTBackendLegacy.hpp"

#include "script_interface/get_value.hpp"

#include <memory>

namespace ScriptInterface {
namespace Dipoles {

template <Arch Architecture>
class DipolarP3M : public Actor<DipolarP3M<Architecture>, ::DipolarP3M> {
  bool m_tune;
  bool m_single_precision;

public:
  using Base = Actor<DipolarP3M<Architecture>, ::DipolarP3M>;
  using Base::actor;
  using Base::add_parameters;
  using Base::context;

protected:
  using Base::m_actor;

public:
  DipolarP3M() {
    add_parameters({
        {"single_precision", AutoParameter::read_only,
         [this]() { return m_single_precision; }},
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
    m_single_precision = get_value<bool>(params, "single_precision");
    static_assert(Architecture == Arch::CPU, "GPU not implemented");
    context()->parallel_try_catch([&]() {
      auto p3m = P3MParameters{!get_value_or<bool>(params, "is_tuned", !m_tune),
                               get_value<double>(params, "epsilon"),
                               get_value<double>(params, "r_cut"),
                               get_value<Utils::Vector3i>(params, "mesh"),
                               get_value<Utils::Vector3d>(params, "mesh_off"),
                               get_value<int>(params, "cao"),
                               get_value<double>(params, "alpha"),
                               get_value<double>(params, "accuracy")};
      make_handle(m_single_precision, std::move(p3m),
                  get_value<double>(params, "prefactor"),
                  get_value<int>(params, "timings"),
                  get_value<bool>(params, "verbose"));
    });
  }

private:
  template <class... Args>
  void make_handle(bool single_precision, Args &&...args) {
    if (single_precision) {
      m_actor = new_dp3m_handle<float, Architecture, FFTBackendLegacy>(
          std::forward<Args>(args)...);
    } else {
      m_actor = new_dp3m_handle<double, Architecture, FFTBackendLegacy>(
          std::forward<Args>(args)...);
    }
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DP3M

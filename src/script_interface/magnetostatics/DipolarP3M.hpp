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
#include "core/p3m/FFTBuffersLegacy.hpp"

#include "script_interface/get_value.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef FFTW3_H
#error "The FFTW3 library shouldn't be visible in this translation unit"
#endif

namespace ScriptInterface {
namespace Dipoles {

template <Arch Architecture>
class DipolarP3M : public Actor<DipolarP3M<Architecture>, ::DipolarP3M> {
  int m_tune_timings;
  bool m_tune;
  bool m_tune_verbose;

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
         [this]() { return not actor()->is_double_precision(); }},
        {"alpha_L", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.alpha_L; }},
        {"r_cut_iL", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.r_cut_iL; }},
        {"mesh", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.mesh; }},
        {"mesh_off", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.mesh_off; }},
        {"cao", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.cao; }},
        {"accuracy", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.accuracy; }},
        {"epsilon", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.epsilon; }},
        {"a", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.a; }},
        {"alpha", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.alpha; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return actor()->dp3m_params.r_cut; }},
        {"is_tuned", AutoParameter::read_only,
         [this]() { return actor()->is_tuned(); }},
        {"verbose", AutoParameter::read_only,
         [this]() { return m_tune_verbose; }},
        {"timings", AutoParameter::read_only,
         [this]() { return m_tune_timings; }},
        {"tune", AutoParameter::read_only, [this]() { return m_tune; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_tune = get_value<bool>(params, "tune");
    m_tune_timings = get_value<int>(params, "timings");
    m_tune_verbose = get_value<bool>(params, "verbose");
    auto const single_precision = get_value<bool>(params, "single_precision");
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
      make_handle(single_precision, std::move(p3m),
                  get_value<double>(params, "prefactor"), m_tune_timings,
                  m_tune_verbose);
    });
  }

private:
  template <typename FloatType, class... Args>
  void make_handle_impl(Args &&...args) {
    m_actor = new_dp3m_handle<FloatType, Architecture, FFTBackendLegacy,
                              FFTBuffersLegacy>(std::forward<Args>(args)...);
  }
  template <class... Args>
  void make_handle(bool single_precision, Args &&...args) {
    if (single_precision) {
      make_handle_impl<float, Args...>(std::forward<Args>(args)...);
    } else {
      make_handle_impl<double, Args...>(std::forward<Args>(args)...);
    }
  }
};

} // namespace Dipoles
} // namespace ScriptInterface

#endif // DP3M

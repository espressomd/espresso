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

#ifdef ESPRESSO_P3M

#include "Actor.hpp"

#include "core/electrostatics/p3m.hpp"
#include "core/electrostatics/p3m.impl.hpp"

#include "script_interface/get_value.hpp"

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace ScriptInterface {
namespace Coulomb {

template <Arch Architecture>
class CoulombP3M : public Actor<CoulombP3M<Architecture>, ::CoulombP3M> {
  int m_tune_timings;
  bool m_tune;
  bool m_tune_verbose;
  bool m_check_complex_residuals;
  bool m_single_precision;
  std::pair<std::optional<int>, std::optional<int>> m_tune_limits;

public:
  using Base = Actor<CoulombP3M<Architecture>, ::CoulombP3M>;
  using Base::actor;
  using Base::add_parameters;
  using Base::context;

protected:
  using Base::m_actor;
  using Base::set_charge_neutrality_tolerance;

public:
  CoulombP3M() {
    add_parameters({
        {"single_precision", AutoParameter::read_only,
         [this]() { return not actor()->is_double_precision(); }},
        {"alpha_L", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.alpha_L; }},
        {"r_cut_iL", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.r_cut_iL; }},
        {"mesh", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.mesh; }},
        {"mesh_off", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.mesh_off; }},
        {"cao", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.cao; }},
        {"accuracy", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.accuracy; }},
        {"epsilon", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.epsilon; }},
        {"a", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.a; }},
        {"alpha", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.alpha; }},
        {"r_cut", AutoParameter::read_only,
         [this]() { return actor()->p3m_params.r_cut; }},
        {"is_tuned", AutoParameter::read_only,
         [this]() { return actor()->is_tuned(); }},
        {"verbose", AutoParameter::read_only,
         [this]() { return m_tune_verbose; }},
        {"timings", AutoParameter::read_only,
         [this]() { return m_tune_timings; }},
        {"tune_limits", AutoParameter::read_only,
         [this]() {
           auto const &[range_min, range_max] = m_tune_limits;
           std::vector<Variant> retval = {
               range_min ? Variant{*range_min} : Variant{None{}},
               range_max ? Variant{*range_max} : Variant{None{}},
           };
           return retval;
         }},
        {"tune", AutoParameter::read_only, [this]() { return m_tune; }},
        {"check_complex_residuals", AutoParameter::read_only,
         [this]() { return m_check_complex_residuals; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_tune = get_value<bool>(params, "tune");
    m_tune_timings = get_value<int>(params, "timings");
    m_tune_verbose = get_value<bool>(params, "verbose");
    m_tune_limits = {std::nullopt, std::nullopt};
    if (params.contains("tune_limits")) {
      auto const &variant = params.at("tune_limits");
      std::size_t range_length = 0u;
      if (is_type<std::vector<int>>(variant)) {
        auto const range = get_value<std::vector<int>>(variant);
        range_length = range.size();
        if (range_length == 2u) {
          m_tune_limits = {range[0u], range[1u]};
        }
      } else {
        auto const range = get_value<std::vector<Variant>>(variant);
        range_length = range.size();
        if (range_length == 2u) {
          if (not is_none(range[0u])) {
            m_tune_limits.first = get_value<int>(range[0u]);
          }
          if (not is_none(range[1u])) {
            m_tune_limits.second = get_value<int>(range[1u]);
          }
        }
      }
      context()->parallel_try_catch([&]() {
        if (range_length != 2u) {
          throw std::invalid_argument("Parameter 'tune_limits' needs 2 values");
        }
        if (m_tune_limits.first and *m_tune_limits.first <= 0) {
          throw std::domain_error("Parameter 'tune_limits' must be > 0");
        }
        if (m_tune_limits.second and *m_tune_limits.second <= 0) {
          throw std::domain_error("Parameter 'tune_limits' must be > 0");
        }
      });
    }
    m_check_complex_residuals =
        get_value<bool>(params, "check_complex_residuals");
    auto const single_precision = get_value<bool>(params, "single_precision");
    context()->parallel_try_catch([&]() {
      if (Architecture == Arch::GPU and not single_precision) {
        throw std::invalid_argument(
            "P3M GPU only implemented in single-precision mode");
      }
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
                  m_tune_verbose, m_tune_limits, m_check_complex_residuals);
    });
    set_charge_neutrality_tolerance(params);
  }

private:
  template <typename FloatType, class... Args>
  void make_handle_impl(Args &&...args) {
    m_actor =
        new_coulomb_p3m<FloatType, Architecture>(std::forward<Args>(args)...);
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

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ESPRESSO_P3M

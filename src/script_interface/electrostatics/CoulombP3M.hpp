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

#include <config/config.hpp>

#ifdef ESPRESSO_P3M

#include "Actor.hpp"

#include "core/electrostatics/p3m.hpp"

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
  TuningParameters m_tuning;
  bool m_tune;
  bool m_single_precision;

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
         [this]() { return m_tuning.verbose; }},
        {"timings", AutoParameter::read_only,
         [this]() { return m_tuning.timings; }},
        {"tune_limits", AutoParameter::read_only,
         [this]() {
#if defined(__clang__) and defined(__cray__)
           auto const &range_min = m_tune_limits.first;
           auto const &range_max = m_tune_limits.second;
#else
           auto const &[range_min, range_max] = m_tuning.limits;
#endif
           std::vector<Variant> retval = {
               range_min ? Variant{*range_min} : Variant{None{}},
               range_max ? Variant{*range_max} : Variant{None{}},
           };
           return retval;
         }},
        {"tune", AutoParameter::read_only, [this]() { return m_tune; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_tune = get_value<bool>(params, "tune");
    m_tuning.timings = get_value<int>(params, "timings");
    m_tuning.verbose = get_value<bool>(params, "verbose");
    m_tuning.limits = {std::nullopt, std::nullopt};
    if (params.contains("tune_limits")) {
      auto const &variant = params.at("tune_limits");
      std::size_t range_length = 0u;
      if (is_type<std::vector<int>>(variant)) {
        auto const range = get_value<std::vector<int>>(variant);
        range_length = range.size();
        if (range_length == 2u) {
          m_tuning.limits = {range[0u], range[1u]};
        }
      } else {
        auto const range = get_value<std::vector<Variant>>(variant);
        range_length = range.size();
        if (range_length == 2u) {
          if (not is_none(range[0u])) {
            m_tuning.limits.first = get_value<int>(range[0u]);
          }
          if (not is_none(range[1u])) {
            m_tuning.limits.second = get_value<int>(range[1u]);
          }
        }
      }
      context()->parallel_try_catch([&]() {
        if (range_length != 2u) {
          throw std::invalid_argument("Parameter 'tune_limits' needs 2 values");
        }
        if (m_tuning.limits.first and *m_tuning.limits.first <= 0) {
          throw std::domain_error("Parameter 'tune_limits' must be > 0");
        }
        if (m_tuning.limits.second and *m_tuning.limits.second <= 0) {
          throw std::domain_error("Parameter 'tune_limits' must be > 0");
        }
      });
    }
    auto const single_precision = get_value<bool>(params, "single_precision");
    context()->parallel_try_catch([&]() {
      auto p3m = P3MParameters{!get_value_or<bool>(params, "is_tuned", !m_tune),
                               get_value<double>(params, "epsilon"),
                               get_value<double>(params, "r_cut"),
                               get_value<Utils::Vector3i>(params, "mesh"),
                               get_value<Utils::Vector3d>(params, "mesh_off"),
                               get_value<int>(params, "cao"),
                               get_value<double>(params, "alpha"),
                               get_value<double>(params, "accuracy")};
      m_actor = new_coulomb_p3m_heffte(std::move(p3m), m_tuning,
                                       get_value<double>(params, "prefactor"),
                                       single_precision, Architecture);
    });
    set_charge_neutrality_tolerance(params);
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ESPRESSO_P3M

/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_LBPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_LBPROFILEOBSERVABLE_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "Observable.hpp"
#include "core/observables/LBProfileObservable.hpp"
#include "script_interface/get_value.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace Observables {

template <typename CoreLBObs>
class LBProfileObservable
    : public AutoParameters<LBProfileObservable<CoreLBObs>, Observable> {
  using Base = AutoParameters<LBProfileObservable<CoreLBObs>, Observable>;

public:
  static_assert(
      std::is_base_of_v<::Observables::LBProfileObservable, CoreLBObs>);
  using Base::Base;
  LBProfileObservable() {
    this->add_parameters(
        {{"n_x_bins", AutoParameter::read_only,
          [this]() {
            return static_cast<int>(profile_observable()->n_bins()[0]);
          }},
         {"n_y_bins", AutoParameter::read_only,
          [this]() {
            return static_cast<int>(profile_observable()->n_bins()[1]);
          }},
         {"n_z_bins", AutoParameter::read_only,
          [this]() {
            return static_cast<int>(profile_observable()->n_bins()[2]);
          }},
         {"min_x", AutoParameter::read_only,
          [this]() { return profile_observable()->limits()[0].first; }},
         {"min_y", AutoParameter::read_only,
          [this]() { return profile_observable()->limits()[1].first; }},
         {"min_z", AutoParameter::read_only,
          [this]() { return profile_observable()->limits()[2].first; }},
         {"max_x", AutoParameter::read_only,
          [this]() { return profile_observable()->limits()[0].second; }},
         {"max_y", AutoParameter::read_only,
          [this]() { return profile_observable()->limits()[1].second; }},
         {"max_z", AutoParameter::read_only,
          [this]() { return profile_observable()->limits()[2].second; }},
         {"sampling_delta_x", AutoParameter::read_only,
          [this]() { return profile_observable()->sampling_delta[0]; }},
         {"sampling_delta_y", AutoParameter::read_only,
          [this]() { return profile_observable()->sampling_delta[1]; }},
         {"sampling_delta_z", AutoParameter::read_only,
          [this]() { return profile_observable()->sampling_delta[2]; }},
         {"sampling_offset_x", AutoParameter::read_only,
          [this]() { return profile_observable()->sampling_offset[0]; }},
         {"sampling_offset_y", AutoParameter::read_only,
          [this]() { return profile_observable()->sampling_offset[1]; }},
         {"sampling_offset_z", AutoParameter::read_only,
          [this]() { return profile_observable()->sampling_offset[2]; }},
         {"allow_empty_bins", AutoParameter::read_only,
          [this]() { return profile_observable()->allow_empty_bins; }}});
  }

  void do_construct(VariantMap const &params) override {
    ObjectHandle::context()->parallel_try_catch([&]() {
      m_observable = std::make_shared<CoreLBObs>(
          get_value_or<double>(params, "sampling_delta_x", 1.),
          get_value_or<double>(params, "sampling_delta_y", 1.),
          get_value_or<double>(params, "sampling_delta_z", 1.),
          get_value_or<double>(params, "sampling_offset_x", 0.),
          get_value_or<double>(params, "sampling_offset_y", 0.),
          get_value_or<double>(params, "sampling_offset_z", 0.),
          get_value<int>(params, "n_x_bins"),
          get_value<int>(params, "n_y_bins"),
          get_value<int>(params, "n_z_bins"),
          get_value<double>(params, "min_x"),
          get_value<double>(params, "max_x"),
          get_value<double>(params, "min_y"),
          get_value<double>(params, "max_y"),
          get_value<double>(params, "min_z"),
          get_value<double>(params, "max_z"),
          get_value_or<bool>(params, "allow_empty_bins", false));
    });
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "edges") {
      std::vector<Variant> variant_edges;
      std::ranges::copy(profile_observable()->edges(),
                        std::back_inserter(variant_edges));
      return variant_edges;
    }
    return Base::do_call_method(method, parameters);
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

  virtual std::shared_ptr<::Observables::LBProfileObservable>
  profile_observable() const {
    return m_observable;
  }

private:
  std::shared_ptr<CoreLBObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <boost/range/algorithm.hpp>

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
    m_observable =
        make_shared_from_args<CoreLBObs, double, double, double, double, double,
                              double, int, int, int, double, double, double,
                              double, double, double, bool>(
            params, "sampling_delta_x", "sampling_delta_y", "sampling_delta_z",
            "sampling_offset_x", "sampling_offset_y", "sampling_offset_z",
            "n_x_bins", "n_y_bins", "n_z_bins", "min_x", "max_x", "min_y",
            "max_y", "min_z", "max_z", "allow_empty_bins");
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "edges") {
      std::vector<Variant> variant_edges;
      boost::copy(profile_observable()->edges(),
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

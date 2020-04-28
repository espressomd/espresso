/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include <boost/range/algorithm.hpp>
#include <iterator>
#include <memory>

#include "Observable.hpp"
#include "core/observables/LBProfileObservable.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename CoreLBObs>
class LBProfileObservable
    : public AutoParameters<LBProfileObservable<CoreLBObs>, Observable> {
  using Base = AutoParameters<LBProfileObservable<CoreLBObs>, Observable>;

public:
  static_assert(
      std::is_base_of<::Observables::LBProfileObservable, CoreLBObs>::value,
      "");
  using Base::Base;
  LBProfileObservable() {
    this->add_parameters(
        {{"n_x_bins",
          [this](const Variant &v) {
            profile_observable()->n_x_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(profile_observable()->n_x_bins);
          }},
         {"n_y_bins",
          [this](const Variant &v) {
            profile_observable()->n_y_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(profile_observable()->n_y_bins);
          }},
         {"n_z_bins",
          [this](const Variant &v) {
            profile_observable()->n_z_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(profile_observable()->n_z_bins);
          }},
         {"min_x",
          [this](const Variant &v) {
            profile_observable()->min_x = get_value<double>(v);
          },
          [this]() { return profile_observable()->min_x; }},
         {"min_y",
          [this](const Variant &v) {
            profile_observable()->min_y = get_value<double>(v);
          },
          [this]() { return profile_observable()->min_y; }},
         {"min_z",
          [this](const Variant &v) {
            profile_observable()->min_z = get_value<double>(v);
          },
          [this]() { return profile_observable()->min_z; }},
         {"max_x",
          [this](const Variant &v) {
            profile_observable()->max_x = get_value<double>(v);
          },
          [this]() { return profile_observable()->max_x; }},
         {"max_y",
          [this](const Variant &v) {
            profile_observable()->max_y = get_value<double>(v);
          },
          [this]() { return profile_observable()->max_y; }},
         {"max_z",
          [this](const Variant &v) {
            profile_observable()->max_z = get_value<double>(v);
          },
          [this]() { return profile_observable()->max_z; }},
         {"sampling_delta_x",
          [this](const Variant &v) {
            profile_observable()->sampling_delta_x = get_value<double>(v);
          },
          [this]() { return profile_observable()->sampling_delta_x; }},
         {"sampling_delta_y",
          [this](const Variant &v) {
            profile_observable()->sampling_delta_y = get_value<double>(v);
          },
          [this]() { return profile_observable()->sampling_delta_y; }},
         {"sampling_delta_z",
          [this](const Variant &v) {
            profile_observable()->sampling_delta_z = get_value<double>(v);
          },
          [this]() { return profile_observable()->sampling_delta_z; }},
         {"sampling_offset_x",
          [this](const Variant &v) {
            profile_observable()->sampling_offset_x = get_value<double>(v);
          },
          [this]() { return profile_observable()->sampling_offset_x; }},
         {"sampling_offset_y",
          [this](const Variant &v) {
            profile_observable()->sampling_offset_y = get_value<double>(v);
          },
          [this]() { return profile_observable()->sampling_offset_y; }},
         {"sampling_offset_z",
          [this](const Variant &v) {
            profile_observable()->sampling_offset_z = get_value<double>(v);
          },
          [this]() { return profile_observable()->sampling_offset_z; }},
         {"allow_empty_bins",
          [this](const Variant &v) {
            profile_observable()->allow_empty_bins = get_value<bool>(v);
          },
          [this]() { return profile_observable()->allow_empty_bins; }}});
  }

  void construct(VariantMap const &params) override {
    m_observable =
        make_shared_from_args<CoreLBObs, double, double, double, double, double,
                              double, int, int, int, double, double, double,
                              double, double, double, bool>(
            params, "sampling_delta_x", "sampling_delta_y", "sampling_delta_z",
            "sampling_offset_x", "sampling_offset_y", "sampling_offset_z",
            "n_x_bins", "n_y_bins", "n_z_bins", "min_x", "min_y", "min_z",
            "max_x", "max_y", "max_z", "allow_empty_bins");
  }

  Variant call_method(std::string const &method,
                      VariantMap const &parameters) override {
    if (method == "edges") {
      std::vector<Variant> variant_edges;
      boost::copy(profile_observable()->edges(),
                  std::back_inserter(variant_edges));
      return variant_edges;
    }
    return Base::call_method(method, parameters);
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

/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SCRIPT_INTERFACE_OBSERVABLES_LBOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_LBOBSERVABLE_HPP

#include "auto_parameters/AutoParameters.hpp"
#include <memory>

#include "core/observables/LBObservable.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class LBObservable
    : virtual public AutoParameters<LBObservable<CoreObs>, Observable> {
public:
  static_assert(std::is_base_of<::Observables::LBObservable, CoreObs>::value,
                "");
  LBObservable() : m_observable(std::make_shared<CoreObs>()) {
    this->add_parameters(
        {{"sampling_delta_x",
          [this](const Variant &v) {
            lb_observable()->sampling_delta_x = get_value<double>(v);
          },
          [this]() { return lb_observable()->sampling_delta_x; }},
         {"sampling_delta_y",
          [this](const Variant &v) {
            lb_observable()->sampling_delta_y = get_value<double>(v);
          },
          [this]() { return lb_observable()->sampling_delta_y; }},
         {"sampling_delta_z",
          [this](const Variant &v) {
            lb_observable()->sampling_delta_z = get_value<double>(v);
          },
          [this]() { return lb_observable()->sampling_delta_z; }},
         {"sampling_offset_x",
          [this](const Variant &v) {
            lb_observable()->sampling_offset_x = get_value<double>(v);
          },
          [this]() { return lb_observable()->sampling_offset_x; }},
         {"sampling_offset_y",
          [this](const Variant &v) {
            lb_observable()->sampling_offset_y = get_value<double>(v);
          },
          [this]() { return lb_observable()->sampling_offset_y; }},
         {"sampling_offset_z",
          [this](const Variant &v) {
            lb_observable()->sampling_offset_z = get_value<double>(v);
          },
          [this]() { return lb_observable()->sampling_offset_z; }},
         {"allow_empty_bins",
          [this](const Variant &v) {
            lb_observable()->allow_empty_bins = get_value<bool>(v);
          },
          [this]() { return lb_observable()->allow_empty_bins; }}});
  }

  std::shared_ptr<::Observables::LBObservable> lb_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

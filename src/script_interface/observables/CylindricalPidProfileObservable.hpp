/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_CYLINDRICALPIDPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_CYLINDRICALPIDPROFILEOBSERVABLE_HPP

#include <memory>

#include "ScriptInterface.hpp"
#include "auto_parameters/AutoParameters.hpp"

#include "Observable.hpp"
#include "core/observables/CylindricalDensityProfile.hpp"
#include "core/observables/CylindricalFluxDensityProfile.hpp"
#include "core/observables/CylindricalLBFluxDensityProfileAtParticlePositions.hpp"
#include "core/observables/CylindricalLBVelocityProfileAtParticlePositions.hpp"
#include "core/observables/CylindricalPidProfileObservable.hpp"
#include "core/observables/CylindricalVelocityProfile.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class CylindricalPidProfileObservable
    : public AutoParameters<CylindricalPidProfileObservable<CoreObs>,
                            Observable> {
public:
  static_assert(std::is_base_of<::Observables::CylindricalPidProfileObservable,
                                CoreObs>::value,
                "");
  CylindricalPidProfileObservable()
      : m_observable(std::make_shared<CoreObs>()) {
    this->add_parameters({
        {"ids",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->ids() =
               get_value<std::vector<int>>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->ids(); }},
        {"center",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->center =
               get_value<::Utils::Vector3d>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->center; }},
        {"axis",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->axis =
               get_value<std::string>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->axis; }},
        {"n_r_bins",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->n_r_bins = get_value<int>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->n_r_bins; }},
        {"n_phi_bins",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->n_phi_bins = get_value<int>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->n_phi_bins; }},
        {"n_z_bins",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->n_z_bins = get_value<int>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->n_z_bins; }},
        {"min_r",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->min_r = get_value<double>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->min_r; }},
        {"min_phi",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->min_phi = get_value<double>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->min_phi; }},
        {"min_z",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->min_z = get_value<double>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->min_z; }},
        {"max_r",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->max_r = get_value<double>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->max_r; }},
        {"max_phi",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->max_phi = get_value<double>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->max_phi; }},
        {"max_z",
         [this](const Variant &v) {
           cylindrical_pid_profile_observable()->max_z = get_value<double>(v);
         },
         [this]() { return cylindrical_pid_profile_observable()->max_z; }},
    });
  };

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

  virtual std::shared_ptr<::Observables::CylindricalPidProfileObservable>
  cylindrical_pid_profile_observable() const {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

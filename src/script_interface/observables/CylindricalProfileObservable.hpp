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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP

#include "ScriptInterface.hpp"

#include <memory>

#include "Observable.hpp"
#include "core/observables/CylindricalDensityProfile.hpp"
#include "core/observables/CylindricalVelocityProfile.hpp"
#include "core/observables/CylindricalFluxDensityProfile.hpp"
#include "core/observables/CylindricalLBFluxDensityProfileAtParticlePositions.hpp"
#include "core/observables/CylindricalLBVelocityProfileAtParticlePositions.hpp"
#include "core/observables/CylindricalProfileObservable.hpp"

namespace ScriptInterface {
namespace Observables {

class CylindricalProfileObservable : public Observable {
public:
  VariantMap get_parameters() const override {
    return {{"ids", cylindrical_profile_observable()->ids()},
            {"center", cylindrical_profile_observable()->center},
            {"axis", cylindrical_profile_observable()->axis},
            {"n_r_bins", cylindrical_profile_observable()->n_r_bins},
            {"n_phi_bins", cylindrical_profile_observable()->n_phi_bins},
            {"n_z_bins", cylindrical_profile_observable()->n_z_bins},
            {"min_r", cylindrical_profile_observable()->min_r},
            {"min_phi", cylindrical_profile_observable()->min_phi},
            {"min_z", cylindrical_profile_observable()->min_z},
            {"max_r", cylindrical_profile_observable()->max_r},
            {"max_phi", cylindrical_profile_observable()->max_phi},
            {"max_z", cylindrical_profile_observable()->max_z}};
  }

  ParameterMap valid_parameters() const override {
    return {{"ids", {ParameterType::INT_VECTOR, true}},
            {"center", {ParameterType::DOUBLE_VECTOR, true}},
            {"axis", {ParameterType::STRING, true}},
            {"n_r_bins", {ParameterType::INT, true}},
            {"n_phi_bins", {ParameterType::INT, true}},
            {"n_z_bins", {ParameterType::INT, true}},
            {"min_r", {ParameterType::DOUBLE, true}},
            {"min_phi", {ParameterType::DOUBLE, true}},
            {"min_z", {ParameterType::DOUBLE, true}},
            {"max_r", {ParameterType::DOUBLE, true}},
            {"max_phi", {ParameterType::DOUBLE, true}},
            {"max_z", {ParameterType::DOUBLE, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    SET_PARAMETER_HELPER("ids", cylindrical_profile_observable()->ids());
    SET_PARAMETER_HELPER("center", cylindrical_profile_observable()->center);
    SET_PARAMETER_HELPER("axis", cylindrical_profile_observable()->axis);
    SET_PARAMETER_HELPER("n_r_bins",
                         cylindrical_profile_observable()->n_r_bins);
    SET_PARAMETER_HELPER("n_phi_bins",
                         cylindrical_profile_observable()->n_phi_bins);
    SET_PARAMETER_HELPER("n_z_bins",
                         cylindrical_profile_observable()->n_z_bins);
    SET_PARAMETER_HELPER("min_r", cylindrical_profile_observable()->min_r);
    SET_PARAMETER_HELPER("min_phi", cylindrical_profile_observable()->min_phi);
    SET_PARAMETER_HELPER("min_z", cylindrical_profile_observable()->min_z);
    SET_PARAMETER_HELPER("max_r", cylindrical_profile_observable()->max_r);
    SET_PARAMETER_HELPER("max_phi", cylindrical_profile_observable()->max_phi);
    SET_PARAMETER_HELPER("max_z", cylindrical_profile_observable()->max_z);
  }

  virtual std::shared_ptr<::Observables::CylindricalProfileObservable>
  cylindrical_profile_observable() const = 0;
};

#define NEW_CYLINDRICAL_PROFILE_OBSERVABLE(obs_name)                           \
  class obs_name : public CylindricalProfileObservable {                       \
  public:                                                                      \
    obs_name() : m_observable(new ::Observables::obs_name()){};                \
                                                                               \
    std::shared_ptr<::Observables::Observable> observable() const override {   \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
    std::shared_ptr<::Observables::CylindricalProfileObservable>               \
    cylindrical_profile_observable() const override {                          \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
  private:                                                                     \
    std::shared_ptr<::Observables::obs_name> m_observable;                     \
  };

NEW_CYLINDRICAL_PROFILE_OBSERVABLE(CylindricalDensityProfile)
NEW_CYLINDRICAL_PROFILE_OBSERVABLE(CylindricalVelocityProfile)
NEW_CYLINDRICAL_PROFILE_OBSERVABLE(CylindricalFluxDensityProfile)
NEW_CYLINDRICAL_PROFILE_OBSERVABLE(
    CylindricalLBFluxDensityProfileAtParticlePositions)
NEW_CYLINDRICAL_PROFILE_OBSERVABLE(
    CylindricalLBVelocityProfileAtParticlePositions)

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

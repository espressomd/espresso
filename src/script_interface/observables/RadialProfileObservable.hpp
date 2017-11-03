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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_RADIALPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_RADIALPROFILEOBSERVABLE_HPP

#include "ScriptInterface.hpp"

#include <memory>

#include "Observable.hpp"
#include "core/observables/RadialProfileObservable.hpp"
#include "core/observables/RadialFluxDensityProfile.hpp"

namespace ScriptInterface {
namespace Observables {

class RadialProfileObservable : public Observable {
public:
  const std::string name() const override {
    return "Observables::RadialProfileObservable";
  };

  VariantMap get_parameters() const override {
    return {{"ids", radial_profile_observable()->ids()},
            {"center", radial_profile_observable()->center},
            {"n_r_bins", radial_profile_observable()->n_r_bins},
            {"n_phi_bins", radial_profile_observable()->n_phi_bins},
            {"n_z_bins", radial_profile_observable()->n_z_bins},
            {"min_r", radial_profile_observable()->min_r},
            {"min_phi", radial_profile_observable()->min_phi},
            {"min_z", radial_profile_observable()->min_z},
            {"max_r", radial_profile_observable()->max_r},
            {"max_phi", radial_profile_observable()->max_phi},
            {"max_z", radial_profile_observable()->max_z}};
  };

  ParameterMap valid_parameters() const override {
    return {{"ids", {ParameterType::INT_VECTOR, true}},
            {"center", {ParameterType::DOUBLE_VECTOR, true}},
            {"n_r_bins", {ParameterType::INT, true}},
            {"n_phi_bins", {ParameterType::INT, true}},
            {"n_z_bins", {ParameterType::INT, true}},
            {"min_r", {ParameterType::DOUBLE, true}},
            {"min_phi", {ParameterType::DOUBLE, true}},
            {"min_z", {ParameterType::DOUBLE, true}},
            {"max_r", {ParameterType::DOUBLE, true}},
            {"max_phi", {ParameterType::DOUBLE, true}},
            {"max_z", {ParameterType::DOUBLE, true}}};
  };

  void set_parameter(std::string const &name, Variant const &value) override {
    SET_PARAMETER_HELPER("ids", radial_profile_observable()->ids());
    SET_PARAMETER_HELPER("n_r_bins", radial_profile_observable()->n_r_bins);
    SET_PARAMETER_HELPER("n_phi_bins", radial_profile_observable()->n_phi_bins);
    SET_PARAMETER_HELPER("n_z_bins", radial_profile_observable()->n_z_bins);
    SET_PARAMETER_HELPER("min_r", radial_profile_observable()->min_r);
    SET_PARAMETER_HELPER("min_phi", radial_profile_observable()->min_phi);
    SET_PARAMETER_HELPER("min_z", radial_profile_observable()->min_z);
    SET_PARAMETER_HELPER("max_r", radial_profile_observable()->max_r);
    SET_PARAMETER_HELPER("max_phi", radial_profile_observable()->max_phi);
    SET_PARAMETER_HELPER("max_z", radial_profile_observable()->max_z);
  };

  virtual std::shared_ptr<::Observables::RadialProfileObservable>
  radial_profile_observable() const = 0;
};

#define NEW_RADIAL_PROFILE_OBSERVABLE(obs_name)                                \
  class obs_name : public RadialProfileObservable {                            \
  public:                                                                      \
    obs_name() : m_observable(new ::Observables::obs_name()){};                \
                                                                               \
    const std::string name() const override {                                  \
      return "Observables::" #obs_name;                                        \
    }                                                                          \
                                                                               \
    std::shared_ptr<::Observables::Observable> observable() const override {   \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
    std::shared_ptr<::Observables::RadialProfileObservable>                    \
    radial_profile_observable() const override {                               \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
  private:                                                                     \
    std::shared_ptr<::Observables::obs_name> m_observable;                     \
  };

NEW_RADIAL_PROFILE_OBSERVABLE(RadialFluxDensityProfile);

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

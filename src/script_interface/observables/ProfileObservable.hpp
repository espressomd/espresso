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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "ScriptInterface.hpp"

#include <memory>

#include "Observable.hpp"
#include "core/observables/DensityProfile.hpp"
#include "core/observables/FluxDensityProfile.hpp"
#include "core/observables/ForceDensityProfile.hpp"
#include "core/observables/LBVelocityProfile.hpp"
#include "core/observables/ProfileObservable.hpp"

namespace ScriptInterface {
namespace Observables {

class ProfileObservable : public Observable {
public:
  VariantMap get_parameters() const override {
    return {{"ids", profile_observable()->ids},

            {"n_x_bins", profile_observable()->n_x_bins},
            {"n_y_bins", profile_observable()->n_y_bins},
            {"n_z_bins", profile_observable()->n_z_bins},

            {"min_x", profile_observable()->min_x},
            {"min_y", profile_observable()->min_y},
            {"min_z", profile_observable()->min_z},

            {"max_x", profile_observable()->max_x},
            {"max_y", profile_observable()->max_y},
            {"max_z", profile_observable()->max_z}};
  };

  ParameterMap valid_parameters() const override {
    return {{"ids", {ParameterType::INT_VECTOR, true}},
            {"n_x_bins", {ParameterType::INT, true}},
            {"n_y_bins", {ParameterType::INT, true}},
            {"n_z_bins", {ParameterType::INT, true}},
            {"min_x", {ParameterType::DOUBLE, true}},
            {"min_y", {ParameterType::DOUBLE, true}},
            {"min_z", {ParameterType::DOUBLE, true}},
            {"max_x", {ParameterType::DOUBLE, true}},
            {"max_y", {ParameterType::DOUBLE, true}},
            {"max_z", {ParameterType::DOUBLE, true}}};
  };

  void set_parameter(std::string const &name, Variant const &value) override {
    SET_PARAMETER_HELPER("ids", profile_observable()->ids);

    SET_PARAMETER_HELPER("n_x_bins", profile_observable()->n_x_bins);
    SET_PARAMETER_HELPER("n_y_bins", profile_observable()->n_y_bins);
    SET_PARAMETER_HELPER("n_z_bins", profile_observable()->n_z_bins);

    SET_PARAMETER_HELPER("min_x", profile_observable()->min_x);
    SET_PARAMETER_HELPER("min_y", profile_observable()->min_y);
    SET_PARAMETER_HELPER("min_z", profile_observable()->min_z);

    SET_PARAMETER_HELPER("max_x", profile_observable()->max_x);
    SET_PARAMETER_HELPER("max_y", profile_observable()->max_y);
    SET_PARAMETER_HELPER("max_z", profile_observable()->max_z);
  };

  virtual std::shared_ptr<::Observables::ProfileObservable>
  profile_observable() const = 0;
};

#define NEW_PROFILE_OBSERVABLE(obs_name)                                       \
  class obs_name : public ProfileObservable {                                  \
  public:                                                                      \
    obs_name() : m_observable(new ::Observables::obs_name()){};                \
                                                                               \
    std::shared_ptr<::Observables::Observable> observable() const override {   \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
    std::shared_ptr<::Observables::ProfileObservable>                          \
    profile_observable() const override {                                      \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
  private:                                                                     \
    std::shared_ptr<::Observables::obs_name> m_observable;                     \
  };

NEW_PROFILE_OBSERVABLE(DensityProfile)
NEW_PROFILE_OBSERVABLE(ForceDensityProfile)
NEW_PROFILE_OBSERVABLE(FluxDensityProfile)
NEW_PROFILE_OBSERVABLE(LBVelocityProfile)

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

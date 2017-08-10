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
#include "core/observables/ProfileObservable.hpp"
#include "core/observables/DensityProfile.hpp"
#include "core/observables/ForceDensityProfile.hpp"
#include "core/observables/FluxDensityProfile.hpp"
#include "core/observables/LbVelocityProfile.hpp"

namespace ScriptInterface {
namespace Observables {

class ProfileObservable : public Observable {
public:
  const std::string name() const override {
    return "Observables::ProfileObservable";
  };

  VariantMap get_parameters() const override {
    return {{"ids", profile_observable()->ids},

            {"xbins", profile_observable()->xbins},
            {"ybins", profile_observable()->ybins},
            {"zbins", profile_observable()->zbins},

            {"minx", profile_observable()->minx},
            {"miny", profile_observable()->miny},
            {"minz", profile_observable()->minz},

            {"maxx", profile_observable()->maxx},
            {"maxy", profile_observable()->maxy},
            {"maxz", profile_observable()->maxz}};
  };

  ParameterMap valid_parameters() const override {
    return {{"ids", {ParameterType::INT_VECTOR, true}},
            {"xbins", {ParameterType::INT, true}},
            {"ybins", {ParameterType::INT, true}},
            {"zbins", {ParameterType::INT, true}},
            {"minx", {ParameterType::DOUBLE, true}},
            {"miny", {ParameterType::DOUBLE, true}},
            {"minz", {ParameterType::DOUBLE, true}},
            {"maxx", {ParameterType::DOUBLE, true}},
            {"maxy", {ParameterType::DOUBLE, true}},
            {"maxz", {ParameterType::DOUBLE, true}}};
  };

  void set_parameter(std::string const &name, Variant const &value) override {
    SET_PARAMETER_HELPER("ids", profile_observable()->ids);

    SET_PARAMETER_HELPER("xbins", profile_observable()->xbins);
    SET_PARAMETER_HELPER("ybins", profile_observable()->ybins);
    SET_PARAMETER_HELPER("zbins", profile_observable()->zbins);

    SET_PARAMETER_HELPER("minx", profile_observable()->minx);
    SET_PARAMETER_HELPER("miny", profile_observable()->miny);
    SET_PARAMETER_HELPER("minz", profile_observable()->minz);

    SET_PARAMETER_HELPER("maxx", profile_observable()->maxx);
    SET_PARAMETER_HELPER("maxy", profile_observable()->maxy);
    SET_PARAMETER_HELPER("maxz", profile_observable()->maxz);
  };

  virtual std::shared_ptr<::Observables::ProfileObservable>
  profile_observable() const = 0;
};

#define NEW_PROFILE_OBSERVABLE(obs_name)                                       \
  class obs_name : public ProfileObservable {                                  \
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
    std::shared_ptr<::Observables::ProfileObservable>                          \
    profile_observable() const override {                                      \
      return m_observable;                                                     \
    }                                                                          \
                                                                               \
  private:                                                                     \
    std::shared_ptr<::Observables::obs_name> m_observable;                     \
  };

NEW_PROFILE_OBSERVABLE(DensityProfile);
NEW_PROFILE_OBSERVABLE(ForceDensityProfile);
NEW_PROFILE_OBSERVABLE(FluxDensityProfile);
NEW_PROFILE_OBSERVABLE(LBVelocityProfile);

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

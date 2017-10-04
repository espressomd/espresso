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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP

#include "ScriptInterface.hpp"

#include <memory>

#include "Observable.hpp"
#include "core/observables/ComForce.hpp"
#include "core/observables/ComPosition.hpp"
#include "core/observables/ComVelocity.hpp"
#include "core/observables/Current.hpp"
#include "core/observables/DipoleMoment.hpp"
#include "core/observables/MagneticDipoleMoment.hpp"
#include "core/observables/ParticleAngularMomentum.hpp"
#include "core/observables/ParticleBodyAngularMomentum.hpp"
#include "core/observables/ParticleBodyVelocities.hpp"
#include "core/observables/ParticleCurrents.hpp"
#include "core/observables/ParticleForces.hpp"
#include "core/observables/ParticlePositions.hpp"
#include "core/observables/ParticleVelocities.hpp"
#include "core/observables/PidObservable.hpp"

namespace ScriptInterface {
namespace Observables {

class PidObservable : public Observable {
public:
  const std::string name() const override {
    return "Observables::PidObservable";
  };

  VariantMap get_parameters() const override {
    return {{"ids", pid_observable()->ids()}};
  };

  ParameterMap valid_parameters() const override {
    return {{"ids", {ParameterType::INT_VECTOR, true}}};
  };

  void set_parameter(std::string const &name, Variant const &value) override {
    SET_PARAMETER_HELPER("ids", pid_observable()->ids());
  };

  virtual std::shared_ptr<::Observables::Observable> observable() const = 0;
  virtual std::shared_ptr<::Observables::PidObservable>
  pid_observable() const = 0;
};

#define NEW_PID_OBSERVABLE(obs_name)                                           \
  class obs_name : public PidObservable {                                      \
  public:                                                                      \
    obs_name() : m_observable(new ::Observables::obs_name()){};                \
                                                                               \
    const std::string name() const override {                                  \
      return "Observables::" #obs_name;                                        \
    }                                                                          \
                                                                               \
    std::shared_ptr<::Observables::Observable> observable() const override {   \
      return m_observable;                                                     \
    };                                                                         \
                                                                               \
    std::shared_ptr<::Observables::PidObservable>                              \
    pid_observable() const override {                                          \
      return m_observable;                                                     \
    };                                                                         \
                                                                               \
  private:                                                                     \
    std::shared_ptr<::Observables::obs_name> m_observable;                     \
  };

NEW_PID_OBSERVABLE(ParticlePositions);
NEW_PID_OBSERVABLE(ParticleVelocities);
NEW_PID_OBSERVABLE(ParticleForces);
NEW_PID_OBSERVABLE(ParticleBodyVelocities);
NEW_PID_OBSERVABLE(ParticleAngularMomentum);
NEW_PID_OBSERVABLE(ParticleBodyAngularMomentum);
NEW_PID_OBSERVABLE(ParticleCurrent);
NEW_PID_OBSERVABLE(Current);
NEW_PID_OBSERVABLE(DipoleMoment);
NEW_PID_OBSERVABLE(MagneticDipoleMoment);
NEW_PID_OBSERVABLE(ComPosition);
NEW_PID_OBSERVABLE(ComVelocity);
NEW_PID_OBSERVABLE(ComForce);

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

/*
  Copyright (C) 2015,2016 The ESPResSo project

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

#include "initialize.hpp"
#include "ScriptInterface.hpp"

#include "CylindricalPidProfileObservable.hpp"
#include "CylindricalLBProfileObservable.hpp"
#include "LBProfileObservable.hpp"
#include "ParamlessObservable.hpp"
#include "PidObservable.hpp"
#include "PidProfileObservable.hpp"
#include "ProfileObservable.hpp"

#include "core/observables/ComForce.hpp"
#include "core/observables/ComPosition.hpp"
#include "core/observables/ComVelocity.hpp"
#include "core/observables/Current.hpp"
#include "core/observables/DipoleMoment.hpp"
#include "core/observables/MagneticDipoleMoment.hpp"
#include "core/observables/ParticleAngularVelocities.hpp"
#include "core/observables/ParticleBodyAngularVelocities.hpp"
#include "core/observables/ParticleBodyVelocities.hpp"
#include "core/observables/ParticleForces.hpp"
#include "core/observables/ParticlePositions.hpp"
#include "core/observables/ParticleVelocities.hpp"
#include "core/observables/CylindricalLBVelocityProfile.hpp"
#include "core/observables/LBVelocityProfile.hpp"


namespace ScriptInterface {
namespace Observables {

#define REGISTER(name)                                                         \
  ScriptInterface::register_new<name>("Observables::" #name "");

#define REGISTER_PROFILE_OBS(name)                                                 \
  ScriptInterface::register_new<ProfileObservable<::Observables::name>>(           \
      "Observables::" #name "");

#define REGISTER_PID_OBS(name)                                                 \
  ScriptInterface::register_new<PidObservable<::Observables::name>>(           \
      "Observables::" #name "");

#define REGISTER_PID_PROFILE_OBS(name)                                                 \
  ScriptInterface::register_new<PidProfileObservable<::Observables::name>>(           \
      "Observables::" #name "");

#define REGISTER_CYLPID_PROFILE_OBS(name)                                                 \
  ScriptInterface::register_new<CylindricalPidProfileObservable<::Observables::name>>(           \
      "Observables::" #name "");

#define REGISTER_CYLLB_OBS(name)                                                 \
  ScriptInterface::register_new<CylindricalLBProfileObservable<::Observables::name>>(           \
      "Observables::" #name "");

#define REGISTER_LB_OBS(name)                                                 \
  ScriptInterface::register_new<LBProfileObservable<::Observables::name>>(           \
      "Observables::" #name "");

void initialize() {
  // Manual registration:
  //  ScriptInterface::register_new<ScriptInterface::Observables::ParticleVelocities>::
  //    register_new("Observables::ParticleVelocities");

  REGISTER(StressTensor);
  REGISTER_PID_OBS(ParticlePositions);
  REGISTER_PID_OBS(ParticleVelocities);
  REGISTER_PID_OBS(ParticleForces);
  REGISTER_PID_OBS(ParticleBodyVelocities);
  REGISTER_PID_OBS(ParticleAngularVelocities);
  REGISTER_PID_OBS(ParticleBodyAngularVelocities);
  REGISTER_PID_OBS(Current);
  REGISTER_PID_OBS(DipoleMoment);
  REGISTER_PID_OBS(MagneticDipoleMoment);
  REGISTER_PID_OBS(ComPosition);
  REGISTER_PID_OBS(ComVelocity);
  REGISTER_PID_OBS(ComForce);
  REGISTER_PID_PROFILE_OBS(DensityProfile);
  REGISTER_PID_PROFILE_OBS(ForceDensityProfile);
  REGISTER_PID_PROFILE_OBS(FluxDensityProfile);
  REGISTER_LB_OBS(LBVelocityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalDensityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalVelocityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalFluxDensityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalLBFluxDensityProfileAtParticlePositions);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalLBVelocityProfileAtParticlePositions);
  REGISTER_CYLLB_OBS(CylindricalLBVelocityProfile);

#undef REGISTER
#undef REGISTER_PID_OBS
}
} /* namespace Observables */
} /* namespace ScriptInterface */

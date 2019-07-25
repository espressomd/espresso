/*
  Copyright (C) 2015-2018 The ESPResSo project

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
#include "CylindricalLBProfileObservable.hpp"
#include "CylindricalPidProfileObservable.hpp"
#include "LBProfileObservable.hpp"
#include "ParamlessObservable.hpp"
#include "PidObservable.hpp"
#include "PidProfileObservable.hpp"
#include "ProfileObservable.hpp"
#include "ScriptInterface.hpp"
#include "config.hpp"

#include "core/observables/ComForce.hpp"
#include "core/observables/ComPosition.hpp"
#include "core/observables/ComVelocity.hpp"
#include "core/observables/Current.hpp"
#include "core/observables/CylindricalDensityProfile.hpp"
#include "core/observables/CylindricalFluxDensityProfile.hpp"
#include "core/observables/CylindricalLBFluxDensityProfileAtParticlePositions.hpp"
#include "core/observables/CylindricalLBVelocityProfile.hpp"
#include "core/observables/CylindricalLBVelocityProfileAtParticlePositions.hpp"
#include "core/observables/CylindricalVelocityProfile.hpp"
#include "core/observables/DipoleMoment.hpp"
#include "core/observables/LBVelocityProfile.hpp"
#include "core/observables/MagneticDipoleMoment.hpp"
#include "core/observables/ParticleAngles.hpp"
#include "core/observables/ParticleAngularVelocities.hpp"
#include "core/observables/ParticleBodyAngularVelocities.hpp"
#include "core/observables/ParticleBodyVelocities.hpp"
#include "core/observables/ParticleDihedrals.hpp"
#include "core/observables/ParticleDistances.hpp"
#include "core/observables/ParticleForces.hpp"
#include "core/observables/ParticlePositions.hpp"
#include "core/observables/ParticleVelocities.hpp"

namespace ScriptInterface {
namespace Observables {

/** @name %Observables registration
 *  Convenience macro functions to automatize the registration of observable
 *  interfaces via @ref ScriptInterface::register_new<T>()
 */
/*@{*/

/** Register a @ref ScriptInterface::Observables::ParamlessObservableInterface
 *  "ParamlessObservableInterface"
 */
#define REGISTER(name)                                                         \
  ScriptInterface::register_new<name>("Observables::" #name "");

/** Register a @ref ScriptInterface::Observables::ProfileObservable
 *  "ProfileObservable"
 */
#define REGISTER_PROFILE_OBS(name)                                             \
  ScriptInterface::register_new<ProfileObservable<::Observables::name>>(       \
      "Observables::" #name "");

/** Register a @ref ScriptInterface::Observables::PidObservable
 *  "PidObservable"
 */
#define REGISTER_PID_OBS(name)                                                 \
  ScriptInterface::register_new<PidObservable<::Observables::name>>(           \
      "Observables::" #name "");

/** Register a @ref ScriptInterface::Observables::PidProfileObservable
 *  "PidProfileObservable"
 */
#define REGISTER_PID_PROFILE_OBS(name)                                         \
  ScriptInterface::register_new<PidProfileObservable<::Observables::name>>(    \
      "Observables::" #name "");

/** Register a @ref
 *  ScriptInterface::Observables::CylindricalPidProfileObservable
 *  "CylindricalPidProfileObservable"
 */
#define REGISTER_CYLPID_PROFILE_OBS(name)                                      \
  ScriptInterface::register_new<                                               \
      CylindricalPidProfileObservable<::Observables::name>>(                   \
      "Observables::" #name "");

/** Register a @ref ScriptInterface::Observables::CylindricalLBProfileObservable
 *  "CylindricalLBProfileObservable"
 */
#define REGISTER_CYLLB_OBS(name)                                               \
  ScriptInterface::register_new<                                               \
      CylindricalLBProfileObservable<::Observables::name>>(                    \
      "Observables::" #name "");

/** Register an @ref ScriptInterface::Observables::LBProfileObservable
 *  "LBProfileObservable"
 */
#define REGISTER_LB_OBS(name)                                                  \
  ScriptInterface::register_new<LBProfileObservable<::Observables::name>>(     \
      "Observables::" #name "");
/*@}*/

void initialize() {
  // Manual registration:
  //  ScriptInterface::register_new<ScriptInterface::Observables::ParticleVelocities>::
  //    register_new("Observables::ParticleVelocities");

  REGISTER(StressTensor);
  REGISTER_PID_OBS(ParticlePositions);
  REGISTER_PID_OBS(ParticleVelocities);
  REGISTER_PID_OBS(ParticleForces);
  REGISTER_PID_OBS(ParticleBodyVelocities);
#ifdef ROTATION
  REGISTER_PID_OBS(ParticleAngularVelocities);
  REGISTER_PID_OBS(ParticleBodyAngularVelocities);
#endif
  REGISTER_PID_OBS(Current);
#ifdef ELECTROSTATICS
  REGISTER_PID_OBS(DipoleMoment);
#endif
#ifdef DIPOLES
  REGISTER_PID_OBS(MagneticDipoleMoment);
#endif
  REGISTER_PID_OBS(ComPosition);
  REGISTER_PID_OBS(ComVelocity);
  REGISTER_PID_OBS(ComForce);
  REGISTER_PID_OBS(ParticleDistances);
  REGISTER_PID_OBS(ParticleAngles);
  REGISTER_PID_OBS(ParticleDihedrals);
  REGISTER_PID_PROFILE_OBS(DensityProfile);
  REGISTER_PID_PROFILE_OBS(ForceDensityProfile);
  REGISTER_PID_PROFILE_OBS(FluxDensityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalDensityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalVelocityProfile);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalFluxDensityProfile);
#ifdef DPD
  REGISTER(DPDStress)
#endif

  REGISTER(LBFluidStress);
  REGISTER_CYLPID_PROFILE_OBS(
      CylindricalLBFluxDensityProfileAtParticlePositions);
  REGISTER_CYLPID_PROFILE_OBS(CylindricalLBVelocityProfileAtParticlePositions);
  REGISTER_CYLLB_OBS(CylindricalLBVelocityProfile);
  REGISTER_LB_OBS(LBVelocityProfile);

#undef REGISTER
#undef REGISTER_PID_OBS
}
} /* namespace Observables */
} /* namespace ScriptInterface */

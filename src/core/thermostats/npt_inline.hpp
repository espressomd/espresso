/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef THERMOSTATS_NPT_INLINE_HPP
#define THERMOSTATS_NPT_INLINE_HPP

#include "config/config.hpp"

#ifdef NPT

#include "random.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <cstddef>

/** Add velocity-dependent noise and friction for NpT-sims to the particle's
 *  velocity
 *  @tparam step       Which half time step to integrate (1 or 2)
 *  @param npt_iso     Parameters
 *  @param vel         particle velocity
 *  @param p_identity  particle identity
 *  @return noise added to the velocity, already rescaled by
 *          dt/2 (contained in prefactors)
 */
template <std::size_t step>
inline Utils::Vector3d
friction_therm0_nptiso(IsotropicNptThermostat const &npt_iso,
                       Utils::Vector3d const &vel, int p_identity) {
  static_assert(step == 1 or step == 2, "NPT only has 2 integration steps");
  constexpr auto const salt =
      (step == 1) ? RNGSalt::NPTISO0_HALF_STEP1 : RNGSalt::NPTISO0_HALF_STEP2;
  if (npt_iso.pref_noise_0 > 0.0) {
    return npt_iso.pref_rescale_0 * vel +
           npt_iso.pref_noise_0 *
               Random::noise_uniform<salt>(npt_iso.rng_counter(),
                                           npt_iso.rng_seed(), p_identity);
  }
  return npt_iso.pref_rescale_0 * vel;
}

/** Add p_diff-dependent noise and friction for NpT-sims to \ref
 *  NptIsoParameters::p_diff
 */
inline double friction_thermV_nptiso(IsotropicNptThermostat const &npt_iso,
                                     double p_diff) {
  if (npt_iso.pref_noise_V > 0.0) {
    return npt_iso.pref_rescale_V * p_diff +
           npt_iso.pref_noise_V *
               Random::noise_uniform<RNGSalt::NPTISOV, 1>(
                   npt_iso.rng_counter(), npt_iso.rng_seed(), 0);
  }
  return npt_iso.pref_rescale_V * p_diff;
}

#endif // NPT
#endif

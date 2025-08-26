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

#ifdef ESPRESSO_NPT

#include "random.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <cstddef>

/** Add velocity-dependent noise and friction for NpT-sims to the particle's
 *  velocity;
 *  @f$ p(t) = p(t) \exp(- \gamma_0 dt / m)
 *                      + \sqrt{k_B T (1 - \exp(-2 \gamma_0 dt)}N(0,1) @f$
 *
 *  @param npt_iso     Parameters
 *  @param vel         particle velocity
 *  @param mass        particle mass
 *  @param p_identity  particle identity
 *  @return velocity added noise
 */
inline Utils::Vector3d
propagate_therm0_nptiso(IsotropicNptThermostat const &npt_iso,
                        Utils::Vector3d const &vel, double mass,
                        int p_identity) {
  if (npt_iso.gamma0 > 0.0) {
    return npt_iso.pref_rescale_0.at(mass) * vel +
           npt_iso.pref_noise_0.at(mass) *
               Random::noise_gaussian<RNGSalt::NPTISO_PARTICLE>(
                   npt_iso.rng_counter(), npt_iso.rng_seed(), p_identity);
  }
  return npt_iso.pref_rescale_0.at(mass) * vel;
}

/** Added noise and friction for NpT-sims to \ref NptIsoParameters::p_epsilon;
 *  @f$ p_{\epsilon} = p_{epsilon}(t) \exp(- \gamma_V dt / W)
 *                     + \sqrt{k_B T (1 - \exp(-2 \gamma_V dt / W)}N(0,1) @f$
 *
 *  @param npt_iso     Parameters
 *  @param p_epsilon    conjugate momentum of volume
 *  @param piston      piston mass
 *  @return conjugate momentum of volume added noise
 */
inline double propagate_thermV_nptiso(IsotropicNptThermostat const &npt_iso,
                                      double p_epsilon, double piston) {
  if (npt_iso.gammav > 0.0) {
    return npt_iso.pref_rescale_V * p_epsilon +
           npt_iso.pref_noise_V *
               Random::noise_gaussian<RNGSalt::NPTISO_VOLUME, 1>(
                   npt_iso.rng_counter(), npt_iso.rng_seed(), 0)[0];
  }
  return npt_iso.pref_rescale_V * p_epsilon;
}

#endif // ESPRESSO_NPT
#endif

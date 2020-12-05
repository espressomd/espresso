/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef ESPRESSO_DIPOLE_INLINE_HPP
#define ESPRESSO_DIPOLE_INLINE_HPP

#include "config.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"

namespace Dipole {
// forces_inline
inline ParticleForce pair_force(Particle const &p1, Particle const &p2,
                                Utils::Vector3d const &d, double dist,
                                double dist2) {
  ParticleForce pf{};
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall trough
  case DIPOLAR_P3M: {
    pf = dp3m_pair_force(p1, p2, d, dist2, dist);
    break;
  }
#endif /*ifdef DP3M */
  default:
    break;
  }
  return pf;
}

// energy_inline
inline double pair_energy(Particle const &p1, Particle const &p2,
                          Utils::Vector3d const &d, double dist, double dist2) {
  double energy = 0;
  if (dipole.method != DIPOLAR_NONE) {
    switch (dipole.method) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
      // fall trough
    case DIPOLAR_P3M:
      energy = dp3m_pair_energy(p1, p2, d, dist2, dist);
      break;
#endif
    default:
      energy = 0;
    }
  }
  return energy;
}

} // namespace Dipole

#endif // ESPRESSO_DIPOLE_INLINE_HPP

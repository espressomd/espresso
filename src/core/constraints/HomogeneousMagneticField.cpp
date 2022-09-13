/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include "HomogeneousMagneticField.hpp"

#include "Observable_stat.hpp"
#include "Particle.hpp"
#include "config.hpp"

#include <utils/Vector.hpp>

namespace Constraints {

ParticleForce HomogeneousMagneticField::force(const Particle &p,
                                              const Utils::Vector3d &, double) {
#ifdef DIPOLES
  return {{}, vector_product(p.calc_dip(), m_field)};
#else
  return {};
#endif
}

void HomogeneousMagneticField::add_energy(const Particle &p,
                                          const Utils::Vector3d &, double,
                                          Observable_stat &obs_energy) const {
#ifdef DIPOLES
  obs_energy.dipolar[0] += -1.0 * m_field * p.calc_dip();
#endif
}

} // namespace Constraints

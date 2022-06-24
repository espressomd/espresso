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
#ifndef OBSERVABLES_PARTICLEVELOCITIES_HPP
#define OBSERVABLES_PARTICLEVELOCITIES_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

/** Extract particle velocities.
 *  For \f$n\f$ particles, return \f$3 n\f$ velocities ordered as
 *  \f$(v_x^1, v_y^1, v_z^1, \dots, v_x^n, v_y^n, v_z^n)\f$.
 */
using ParticleVelocities = ParticleObservable<ParticleObservables::Velocities>;

} // namespace Observables
#endif

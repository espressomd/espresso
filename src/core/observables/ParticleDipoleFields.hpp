/*
 * Copyright (C) 2023 The ESPResSo project
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
#ifndef OBSERVABLES_PARTICLEDIPOLEFIELDS_HPP
#define OBSERVABLES_PARTICLEDIPOLEFIELDS_HPP

#include "config/config.hpp"

#include "PidObservable.hpp"
#include "energy.hpp"

#include <vector>

namespace Observables {

/** Extract particle dipole fields.
 *  For \f$n\f$ particles, return \f$3 n\f$ dipole fields ordered as
 *  \f$(h_d1_x, h_d1_y, h_d1_z, \dots, h_dn_x, h_dn_y, h_dn_z)\f$.
 */
class ParticleDipoleFields
    : public ParticleObservable<ParticleObservables::DipoleFields> {
public:
  using ParticleObservable<
      ParticleObservables::DipoleFields>::ParticleObservable;
  std::vector<double>
  evaluate(ParticleReferenceRange particles,
           const ParticleObservables::traits<Particle> &traits) const override {
#ifdef DIPOLE_FIELD_TRACKING
    mpi_calc_long_range_fields();
#endif
    return ParticleObservable<ParticleObservables::DipoleFields>::evaluate(
        particles, traits);
  }
};

} // namespace Observables
#endif

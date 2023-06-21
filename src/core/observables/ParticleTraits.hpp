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
#ifndef OBSERVABLES_PARTICLE_TRAITS
#define OBSERVABLES_PARTICLE_TRAITS

#include "Particle.hpp"
#include "config/config.hpp"
#include "rotation.hpp"

namespace ParticleObservables {
/**
 * Template specialization for `Particle`. The traits mechanism is used to get
 * indirect access to particle properties. This helps making the implementation
 * of observables independent of the particle type.
 */
template <> struct traits<Particle> {
  auto position(Particle const &p) const { return p.pos(); }
  auto velocity(Particle const &p) const { return p.v(); }
  auto force(Particle const &p) const { return p.force(); }
  auto mass(Particle const &p) const {
#ifdef VIRTUAL_SITES
    // we exclude virtual particles since their mass does not have a meaning
    if (p.is_virtual())
      return decltype(p.mass()){};
#endif
    return p.mass();
  }
  auto charge(Particle const &p) const { return p.q(); }
  auto dipole_moment(Particle const &p) const {
#ifdef DIPOLES
    return p.calc_dip();
#else
    return Utils::Vector3d{};
#endif
  }
  auto dipole_field(Particle const &p) const {
#ifdef DIPOLE_FIELD_TRACKING
    return p.dip_fld();
#else
    return Utils::Vector3d{};
#endif
  }
  auto velocity_body(Particle const &p) const {
#ifdef ROTATION
    return convert_vector_space_to_body(p, p.v());
#else
    return Utils::Vector3d{};
#endif
  }
  auto angular_velocity(Particle const &p) const {
#ifdef ROTATION
    return convert_vector_body_to_space(p, p.omega());
#else
    return Utils::Vector3d{};
#endif
  }
  auto angular_velocity_body(Particle const &p) const {
#ifdef ROTATION
    return p.omega();
#else
    return Utils::Vector3d{};
#endif
  }
  auto director(Particle const &p) const {
#ifdef ROTATION
    return p.calc_director();
#else
    return Utils::Vector3d{{0., 0., 1.}};
#endif
  }
};

} // namespace ParticleObservables

#endif

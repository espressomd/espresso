/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#include "config.hpp"

namespace GenObs {
template <> struct traits<Particle> {
  auto position(Particle const &p) const { return p.r.p; }
  auto velocity(Particle const &p) const { return p.m.v; }
  auto mass(Particle const &p) const {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      return decltype(p.p.mass){};
#endif
    return p.p.mass;
  }
  auto charge(Particle const &p) const { return p.p.q; }
  auto force(Particle const &p) const {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      return decltype(p.f.f){};
#endif
    return p.f.f;
  }
  auto dipole_moment(Particle const &p) const {
#if defined(ROTATION) && defined(DIPOLES)
    return p.calc_dip();
#else
    return Utils::Vector3d{};
#endif
  }
};

} // namespace GenObs

#endif

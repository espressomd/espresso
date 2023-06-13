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
#ifndef OBSERVABLES_PROPERTIES_HPP
#define OBSERVABLES_PROPERTIES_HPP

/** @file properties.hpp
 * This file contains interface functionality for extracting particle properties
 * via a particle traits class.
 */
namespace ParticleObservables {
template <class DoF> struct traits;

namespace detail {
template <class T> struct decay { using type = typename std::decay_t<T>; };

template <class U> struct decay<std::reference_wrapper<U>> {
  using type = std::decay_t<U>;
};

template <class T> using decay_t = typename decay<T>::type;
} // namespace detail

template <class Particle>
using default_traits = traits<detail::decay_t<Particle>>;

struct Force {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.force(p);
  }
};

struct Position {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.position(p);
  }
};

struct Velocity {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.velocity(p);
  }
};

struct Director {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.director(p);
  }
};

struct BodyVelocity {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.velocity_body(p);
  }
};

struct AngularVelocity {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.angular_velocity(p);
  }
};

struct BodyAngularVelocity {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.angular_velocity_body(p);
  }
};

struct Mass {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.mass(p);
  }
};

struct Charge {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.charge(p);
  }
};

struct DipoleMoment {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.dipole_moment(p);
  }
};

struct DipoleField {
  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits particle_traits = {}) const {
    return particle_traits.dipole_field(p);
  }
};
} // namespace ParticleObservables

#endif // OBSERVABLES_PROPERTIES_HPP

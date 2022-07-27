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
#ifndef TESTS_MOCK_HPP
#define TESTS_MOCK_HPP

namespace Testing {

struct Particle {
  double position = 1.;
  double m_vel = 2.;
  double m_force = 2.1;
  double m_dipole_moment = 1.2;
  double mass = 3.;
  double charge = 0.0;
  auto const &velocity() const { return m_vel; }
  auto const &dipole_moment() const { return m_dipole_moment; }
};
} // namespace Testing

namespace ParticleObservables {

template <> struct traits<Testing::Particle> {
  using Particle = Testing::Particle;

  double position(Particle const &p) const { return p.position; }
  double velocity(Particle const &p) const { return p.velocity(); }
  double mass(Particle const &p) const { return p.mass; }
  double charge(Particle const &p) const { return p.charge; }
  double dipole_moment(Particle const &p) const { return p.dipole_moment(); }
};
} // namespace ParticleObservables

#endif

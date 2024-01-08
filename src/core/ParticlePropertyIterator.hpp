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

#pragma once

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <utils/Vector.hpp>

namespace ParticlePropertyRange {

namespace detail {
template <class Kernel>
auto create_transform_range(ParticleRange const &particles, Kernel kernel) {
  auto transform_iterator_begin =
      boost::make_transform_iterator(particles.begin(), kernel);
  auto transform_iterator_end =
      boost::make_transform_iterator(particles.end(), kernel);
  return boost::make_iterator_range<decltype(transform_iterator_begin)>(
      transform_iterator_begin, transform_iterator_end);
}
} // namespace detail

inline auto unfolded_pos_range(ParticleRange const &particles,
                               BoxGeometry const &box) {
  auto return_unfolded_pos = [&box](Particle &p) {
    return ::detail::unfolded_position(p.pos(), p.image_box(), box.length());
  };
  return detail::create_transform_range(particles, return_unfolded_pos);
}

inline auto pos_range(ParticleRange const &particles) {
  auto return_pos = [](Particle &p) -> Utils::Vector3d & { return p.pos(); };
  return detail::create_transform_range(particles, return_pos);
}

inline auto charge_range(ParticleRange const &particles) {
  auto return_charge = [](Particle &p) -> double & { return p.q(); };
  return detail::create_transform_range(particles, return_charge);
}

inline auto force_range(ParticleRange const &particles) {
  auto return_force = [](Particle &p) -> Utils::Vector3d & {
    return p.force();
  };
  return detail::create_transform_range(particles, return_force);
}

} // namespace ParticlePropertyRange

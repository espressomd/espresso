/*
 * Copyright (C) 2016-2019 The ESPResSo project
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
#include "CylindricalLBVelocityProfileAtParticlePositions.hpp"

#include "BoxGeometry.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"

#include <utils/Histogram.hpp>
#include <utils/Span.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <cstddef>
#include <vector>

namespace Observables {
std::vector<double> CylindricalLBVelocityProfileAtParticlePositions::evaluate(
    Utils::Span<std::reference_wrapper<const Particle>> particles,
    const ParticleObservables::traits<Particle> &traits) const {
  Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());

  for (auto p : particles) {
    auto const pos = folded_position(traits.position(p), box_geo);
    auto const v = lb_lbfluid_get_interpolated_velocity(pos) *
                   lb_lbfluid_get_lattice_speed();

    histogram.update(
        Utils::transform_coordinate_cartesian_to_cylinder(
            pos - transform_params->center(), transform_params->axis(),
            transform_params->orientation()),
        Utils::transform_vector_cartesian_to_cylinder(
            v, transform_params->axis(), pos - transform_params->center()));
  }

  // normalize by number of hits per bin
  auto hist_tmp = histogram.get_histogram();
  auto tot_count = histogram.get_tot_count();
  std::transform(hist_tmp.begin(), hist_tmp.end(), tot_count.begin(),
                 hist_tmp.begin(), [](auto hi, auto ci) {
                   return ci > 0 ? hi / static_cast<double>(ci) : 0.;
                 });
  return hist_tmp;
}

} // namespace Observables

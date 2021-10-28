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
#ifndef OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP

#include "BoxGeometry.hpp"
#include "CylindricalPidProfileObservable.hpp"
#include "grid.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace Observables {
class CylindricalFluxDensityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;

  std::vector<double>
  evaluate(Utils::Span<std::reference_wrapper<const Particle>> particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());

    // Write data to the histogram
    for (auto p : particles) {
      auto const pos = folded_position(traits.position(p), box_geo) -
                       transform_params->center();
      histogram.update(
          Utils::transform_coordinate_cartesian_to_cylinder(
              pos, transform_params->axis(), transform_params->orientation()),
          Utils::transform_vector_cartesian_to_cylinder(
              traits.velocity(p), transform_params->axis(), pos));
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
  std::vector<std::size_t> shape() const override {
    auto const b = n_bins();
    return {b[0], b[1], b[2], 3};
  }
};

} // Namespace Observables

#endif

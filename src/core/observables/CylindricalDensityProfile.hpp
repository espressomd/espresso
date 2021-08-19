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
#ifndef OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP

#include "CylindricalPidProfileObservable.hpp"

#include "BoxGeometry.hpp"
#include "grid.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <cstddef>
#include <vector>

namespace Observables {
class CylindricalDensityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;
  std::vector<double>
  evaluate(Utils::Span<std::reference_wrapper<const Particle>> particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    Utils::CylindricalHistogram<double, 1> histogram(n_bins(), limits());

    for (auto p : particles) {
      histogram.update(Utils::transform_coordinate_cartesian_to_cylinder(
          folded_position(traits.position(p), box_geo) -
              transform_params->center(),
          transform_params->axis(), transform_params->orientation()));
    }

    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif

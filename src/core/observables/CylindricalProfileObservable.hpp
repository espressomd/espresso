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
#ifndef OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP

#include "ProfileObservable.hpp"

#include <utils/Vector.hpp>
#include <utils/math/abs.hpp>
#include <utils/math/cylindrical_transformation_parameters.hpp>
#include <utils/math/make_lin_space.hpp>

#include <boost/range/algorithm.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace Observables {

/** Cylindrical profile observable */
class CylindricalProfileObservable : public ProfileObservable {
public:
  CylindricalProfileObservable(
      std::shared_ptr<Utils::CylindricalTransformationParameters>
          transform_params,
      int n_r_bins, int n_phi_bins, int n_z_bins, double min_r, double max_r,
      double min_phi, double max_phi, double min_z, double max_z)
      : ProfileObservable(n_r_bins, n_phi_bins, n_z_bins, min_r, max_r, min_phi,
                          max_phi, min_z, max_z),
        transform_params(std::move(transform_params)) {}

  std::shared_ptr<Utils::CylindricalTransformationParameters> transform_params;
};

} // Namespace Observables
#endif

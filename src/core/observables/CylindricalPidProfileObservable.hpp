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
#ifndef OBSERVABLES_CYLINDRICALPIDPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPIDPROFILEOBSERVABLE_HPP

#include <utility>

#include "CylindricalProfileObservable.hpp"
#include "PidObservable.hpp"

namespace Observables {

class CylindricalPidProfileObservable : public PidObservable,
                                        public CylindricalProfileObservable {
public:
  CylindricalPidProfileObservable(
      std::vector<int> const &ids,
      std::shared_ptr<Utils::CylindricalTransformationParameters>
          transform_params,
      int n_r_bins, int n_phi_bins, int n_z_bins, double min_r, double max_r,
      double min_phi, double max_phi, double min_z, double max_z)
      : PidObservable(ids),
        CylindricalProfileObservable(std::move(transform_params), n_r_bins,
                                     n_phi_bins, n_z_bins, min_r, max_r,
                                     min_phi, max_phi, min_z, max_z) {}
};

} // Namespace Observables
#endif

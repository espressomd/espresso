/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP

#include <utils/Vector.hpp>

#include <string>

namespace Observables {
class CylindricalProfileObservable {
public:
  ::Utils::Vector3d center;
  std::string axis;
  double min_r, max_r;
  double min_phi, max_phi;
  double min_z, max_z;
  // Number of bins for each coordinate.
  int n_r_bins, n_phi_bins, n_z_bins;
  double r_bin_size() const { return (max_r - min_r) / n_r_bins; }
  double phi_bin_size() const { return (max_phi - min_phi) / n_phi_bins; }
  double z_bin_size() const { return (max_z - min_z) / n_z_bins; }
};

} // Namespace Observables
#endif

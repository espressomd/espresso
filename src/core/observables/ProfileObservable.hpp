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
#ifndef OBSERVABLES_PROFILEOBSERVABLE_HPP
#define OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "Observable.hpp"

namespace Observables {

// Observable which acts on a given list of particle ids
class ProfileObservable : virtual public Observable {
public:
  ProfileObservable(double min_x, double max_x, double min_y, double max_y,
                    double min_z, double max_z, int n_x_bins, int n_y_bins,
                    int n_z_bins)
      : min_x(min_x), max_x(max_x), min_y(min_y), max_y(max_y), min_z(min_z),
        max_z(max_z), n_x_bins(n_x_bins), n_y_bins(n_y_bins),
        n_z_bins(n_z_bins) {}
  double min_x, max_x;
  double min_y, max_y;
  double min_z, max_z;
  int n_x_bins, n_y_bins, n_z_bins;
  int n_values() const override { return n_x_bins * n_y_bins * n_z_bins; };
};

} // Namespace Observables
#endif

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
#ifndef CORE_UNIT_TESTS_COMMON_GAUSSIAN_HPP
#define CORE_UNIT_TESTS_COMMON_GAUSSIAN_HPP

#include "utils/Vector.hpp"

#include <boost/multi_array.hpp>

#include <cmath>

inline double gaussian(Utils::Vector3d x, Utils::Vector3d x0, double sigma) {
  return std::exp(-((x - x0).norm2() / (2. * sigma * sigma)));
}

inline Utils::Vector3d del_gaussian(Utils::Vector3d x, Utils::Vector3d x0,
                                    double sigma) {
  return -(x - x0) * gaussian(x, x0, sigma) / (sigma * sigma);
}

inline boost::multi_array<double, 3> gaussian_field(int size, Utils::Vector3d h,
                                                    Utils::Vector3d origin,
                                                    Utils::Vector3d x0,
                                                    double sigma) {
  boost::multi_array<double, 3> data(Utils::Vector3i{10, 10, 10});
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++) {
        auto const x = origin + Utils::Vector3d{i * h[0], j * h[1], k * h[2]};
        data[i][j][k] = gaussian(x, x0, sigma);
      }

  return data;
}

#endif

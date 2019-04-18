/*
  Copyright (C) 2019 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::raster test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/raster.hpp>

/* The basic functionality is tested by just
 * recording the x values for every grid point. */
BOOST_AUTO_TEST_CASE(basic_test) {
  using Utils::raster;
  using Utils::Vector3d;
  using Utils::Vector3i;

  Vector3d const origin{1., 2., 3.};
  Vector3d const h{2., 3., 4.};
  Vector3i size{10, 11, 12};

  auto const res = raster<Vector3d>(origin, h, size, [](auto i) { return i; });

  for (int i = 0; i < size[0]; i++)
    for (int j = 0; j < size[1]; j++)
      for (int k = 0; k < size[2]; k++) {
        auto const expected = origin + Vector3d{i * h[0], j * h[1], k * h[2]};
        BOOST_CHECK(res[i][j][k] == expected);
      }
}
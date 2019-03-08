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
#define BOOST_TEST_MODULE d3q19 test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "grid_based_algorithms/lb-d3q19.hpp"

BOOST_AUTO_TEST_CASE(d3q19) {
  for (int i = 0; i < 19; ++i) {
    for (int j = 0; j < 19; ++j) {
      BOOST_CHECK(D3Q19::e_ki[i][j] == D3Q19::e_ki_transposed[j][i]);
    }
  }
}

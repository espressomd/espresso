/*
   Copyright (C) 2018 The ESPResSo project

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

/** \file NumeratedContainer_test.cpp Unit tests for the
 * Utils::NumeratedContainer class.
 *
 */

#define BOOST_TEST_MODULE Utils::permute test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/math/permute_ifield.hpp"

#include <array>
#include <numeric>

BOOST_AUTO_TEST_CASE(check_permu) {
  std::array<int, 11> f;
  std::iota(f.begin(), f.end(), 0);

  for (int i = -11; i <= 11; i++) {
    std::array<int, 11> h(f);
    Utils::permute_ifield(h.data(), h.size(), i);

    for (int j = 0; j < 11; j++)
      BOOST_CHECK(h[j] == ((i + j + 11) % 11));
  }
}

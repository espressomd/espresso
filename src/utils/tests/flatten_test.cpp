/*
 * Copyright (C) 2015-2020 The ESPResSo project
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

#define BOOST_TEST_MODULE flatten test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/flatten.hpp>

#include <algorithm>
#include <array>
#include <iterator>
#include <list>
#include <vector>

BOOST_AUTO_TEST_CASE(flatten_) {
  using Utils::flatten;

  /* not nested */
  {
    const std::array<int, 4> in = {{1, 2, 3, 4}};
    std::array<int, 4> out{};
    flatten(in, out.begin());
    BOOST_CHECK_EQUAL_COLLECTIONS(in.begin(), in.end(), out.begin(), out.end());
  }

  /* nested */
  {
    const std::array<std::array<int, 2>, 2> in{{{{1, 2}}, {{3, 4}}}};
    std::array<int, 4> out{};
    flatten(in, out.begin());

    const std::array<int, 4> expected = {{1, 2, 3, 4}};
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(), out.begin(),
                                  out.end());
  }
  {
    const std::vector<int> in = {{1, 2, 3, 4}};
    std::vector<int> out;
    flatten(in, std::back_inserter(out));
    BOOST_CHECK_EQUAL_COLLECTIONS(in.begin(), in.end(), out.begin(), out.end());
  }
  {
    const std::vector<int> in = {{1, 2, 3, 4}};
    std::list<int> out;
    flatten(in, std::front_inserter(out));
    BOOST_CHECK_EQUAL_COLLECTIONS(in.rbegin(), in.rend(), out.begin(),
                                  out.end());
  }
}

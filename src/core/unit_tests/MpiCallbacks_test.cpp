/*
  Copyright (C) 2016-2018 The ESPResSo project
    Max-Planck-Institute for Polymer Research, Theory Group

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

/** \file
 * Unit tests for the MpiCallbacks class.
 *
 */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE MpiCallbacks test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "../MpiCallbacks.hpp"

#include <boost/mpi.hpp>

struct Archive {
    std::vector<int> values;
    int cnt = 0;

    void operator&(int &i) {
      values.push_back(i);
      i = cnt++;
    }
};

BOOST_AUTO_TEST_CASE(tuple_helpers) {
  auto t = Communication::detail::tuple<int, int, int>{{1, 3, 5}};

  Archive ar;
  t.serialize(ar, 0);

  BOOST_CHECK_EQUAL(ar.values[0], 1);
  BOOST_CHECK_EQUAL(ar.values[1], 3);
  BOOST_CHECK_EQUAL(ar.values[2], 5);

  BOOST_CHECK_EQUAL(std::get<0>(t.t), 0);
  BOOST_CHECK_EQUAL(std::get<1>(t.t), 1);
  BOOST_CHECK_EQUAL(std::get<2>(t.t), 2);
}

BOOST_AUTO_TEST_CASE(callback_handle) {
  ;
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

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

#define BOOST_TEST_MODULE Utils::type_id test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/type_id.hpp"

BOOST_AUTO_TEST_CASE(type_id_test) {
  using Utils::type_id;

  auto const t_i = type_id<int>();
  auto const t_f = type_id<float>();

  /* Different value for different types */
  BOOST_CHECK_NE(t_i, t_f);
  /* Stable value for one type */
  BOOST_CHECK_EQUAL(t_i, type_id<int>());
}
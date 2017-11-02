/*
  Copyright (C) 2017 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::keys test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/utils/keys.hpp"

#include <algorithm>
#include <map>

BOOST_AUTO_TEST_CASE(values) {
  using pair_type = std::map<int,int>::value_type;

  auto pairs = std::map<int, int>{{31, 2}, {2, 63}, {23, 9}, {4, 9}};
  auto keys = Utils::keys(pairs);

  BOOST_CHECK(
      std::equal(keys.begin(), keys.end(), pairs.begin(),
                 [](int k, pair_type const &kv) { return k == kv.first; }));
}


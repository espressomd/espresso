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

#define BOOST_TEST_MODULE Utils::cartesian_product test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/cartesian_product.hpp"

#include <array>
#include <string>

BOOST_AUTO_TEST_CASE(cart_prod) {
    std::array<std::string, 2> strings{"Apple", "Banana"};
    std::array<int, 4> numbers{5,6,7};

    std::vector<std::pair<std::string, int>> prod;
    Utils::cartesian_product([&prod](const std::string&s, int i) {
        prod.emplace_back(s, i);
    }, strings, numbers);

    auto it = prod.begin();
    for(auto const&s: strings) {
        for(auto const&i: numbers) {
            auto const expected = std::pair<std::string, int>{s, i};
            BOOST_CHECK(expected == *it++);
        }
    }

    /* empty */
    {
        auto called = false;
        Utils::cartesian_product([&](){ called = true; });
        BOOST_CHECK(not called);
    }
}

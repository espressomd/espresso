/*
  Copyright (C) 2018 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::uniform test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/uniform.hpp>

BOOST_AUTO_TEST_CASE(limits) {
  BOOST_CHECK_SMALL(Utils::uniform(0ul),
                    std::numeric_limits<double>::epsilon());
  BOOST_CHECK_EQUAL(1., Utils::uniform(std::numeric_limits<uint64_t>::max()));
  BOOST_CHECK_EQUAL(Utils::uniform(0ul) - Utils::uniform(5ul),
                    Utils::uniform(10000ul) - Utils::uniform(10005ul));
}

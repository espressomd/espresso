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

#define BOOST_TEST_MODULE Utils::Batch test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Span.hpp"
using Utils::Span;

#include <vector>
#include <numeric>

BOOST_AUTO_TEST_CASE(ctor) {
  /* from ptr + size */
  {
  std::vector<int> v(23);

  auto s = Span<int>(v.data(), v.size());

  BOOST_CHECK(v.size() == s.size());
  BOOST_CHECK(v.data() == s.data());
  }

  /* from ptr + size */
  {
    const std::vector<int> v(23);

    auto s = Span<const int>(v.data(), v.size());

    BOOST_CHECK(v.size() == s.size());
    BOOST_CHECK(v.data() == s.data());
  }
}

BOOST_AUTO_TEST_CASE(iterators) {
  int dummy;
  auto const p = &dummy;
  auto const size = 11u;

  auto s = Span<int>(p, size);

  BOOST_CHECK(s.begin() == p);
  BOOST_CHECK(s.cbegin() == p);
  BOOST_CHECK(s.end() == (p + size));
  BOOST_CHECK(s.cend() == (p + size));
}

BOOST_AUTO_TEST_CASE(element_access) {
  std::vector<int> v(11);
  std::iota(v.begin(), v.end(), 5);

  auto s = Span<int>(v.data(), v.size());

  for (Span<int>::size_type i = 0; i < s.size(); ++i) {
    BOOST_CHECK(v.at(i) == s[i]);
    BOOST_CHECK(v.at(i) == s.at(i));
  }

  BOOST_CHECK_THROW(s.at(s.size()), std::out_of_range);
}


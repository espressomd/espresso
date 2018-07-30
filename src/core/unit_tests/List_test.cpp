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

/** \file Vector_test.cpp Unit tests for the Utils::Vector class.
 *
*/

#define BOOST_TEST_MODULE List test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iterator>
#include <numeric>
#include <sstream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "core/utils/List.hpp"
#include "core/utils/serialization/List.hpp"

using List = Utils::List<int>;

BOOST_AUTO_TEST_CASE(constructors) {
  /* List::List() */
  {
    auto l = List();
    BOOST_CHECK(l.size() == 0);
    BOOST_CHECK(l.empty());
    BOOST_CHECK(l.begin() == l.end());
    BOOST_CHECK(l.capacity() == 0);
  }

  /* List::List(size_type) */
  {
    auto l = List(35);
    BOOST_CHECK(l.size() == 35);
    BOOST_CHECK(std::distance(l.begin(), l.end()) == 35);
    BOOST_CHECK(l.capacity() == 35);
  }

  /* List::List(size_type, value_type const&) */
  {
    auto l = List(35, 42);
    BOOST_CHECK(l.size() == 35);
    BOOST_CHECK(std::distance(l.begin(), l.end()) == 35);
    BOOST_CHECK(l.capacity() == 35);

    for (auto const &e : l) {
      BOOST_CHECK(e == 42);
    }
  }

  /* List::List(List const&) */
  {
    auto l = List(23);
    std::iota(l.begin(), l.end(), 0);

    auto u(l);
    BOOST_CHECK(std::equal(l.begin(), l.end(), u.begin()));
  }

  /* List::operator=(List const&) */
  {
    auto l = List(23);
    std::iota(l.begin(), l.end(), 0);

    List u;
    u = l;
    BOOST_CHECK(std::equal(l.begin(), l.end(), u.begin()));
  }

  /* List::List(List &&) */
  {
    auto l = List(23);
    std::iota(l.begin(), l.end(), 0);

    auto v(l);
    auto const d_ptr = v.data();

    auto u(std::move(v));

    /* Check that there was no realloc */
    BOOST_CHECK(d_ptr == u.data());
    BOOST_CHECK(std::equal(l.begin(), l.end(), u.begin()));
  }

  /* List::operator=(List &&) */
  {
    auto l = List(23);
    std::iota(l.begin(), l.end(), 0);

    auto v(l);
    auto const d_ptr = v.data();

    List u;
    u = std::move(v);

    /* Check that there was no realloc */
    BOOST_CHECK(d_ptr == u.data());
    BOOST_CHECK(std::equal(l.begin(), l.end(), u.begin()));
  }

  /* List::List(std::initializer_list) */
  {
    auto il = {1, 2, 3, 4};
    auto l = List(il);

    BOOST_CHECK(l.size() == il.size());
    BOOST_CHECK(std::equal(l.begin(), l.end(), il.begin()));
  }
}

BOOST_AUTO_TEST_CASE(iterators) {
  /* List::begin() */
  {
    auto l = List(23);

    BOOST_CHECK(l.begin() == l.data());
  }

  /* List::begin() const */
  {
    auto const l = List(23);

    BOOST_CHECK(l.begin() == l.data());
  }

  /* List::end() */
  {
    auto l = List(23);

    BOOST_CHECK(l.end() == l.data() + l.size());
  }

  /* List::end() const */
  {
    auto const l = List(23);

    BOOST_CHECK(l.end() == l.data() + l.size());
  }
}

BOOST_AUTO_TEST_CASE(front_back) {
  /* List::front() */
  {
    auto l = List(22);
    std::iota(l.begin(), l.end(), 11);

    BOOST_CHECK(l.front() == 11);
  }

  /* List::front() const */
  {
    auto i = List(22);
    std::iota(i.begin(), i.end(), 11);
    auto const l = i;

    BOOST_CHECK(l.front() == *l.begin());
  }

  /* List::back() */
  {
    auto l = List(22);
    std::iota(l.begin(), l.end(), 11);

    BOOST_CHECK(l.back() == 32);
  }

  /* List::back() const */
  {
    auto i = List(22);
    std::iota(i.begin(), i.end(), 11);

    auto const l = i;
    BOOST_CHECK(l.back() == 32);
  }
}

BOOST_AUTO_TEST_CASE(clear) {
  /* List::clear() with empty list */
  {
    auto l = List();

    l.clear();
    BOOST_CHECK(l.empty());
    BOOST_CHECK(l.begin() == l.end());
  }

  /* List::clear() with full list */
  {
    auto l = List(14);

    l.clear();
    BOOST_CHECK(l.empty());
    BOOST_CHECK(l.begin() == l.end());
  }
}

BOOST_AUTO_TEST_CASE(reserve) {
  /* List::reserve() with size < capacity */
  {
    auto l = List(128);

    l.reserve(127);
    BOOST_CHECK(l.capacity() == 128);
    BOOST_CHECK(l.size() == 128);
  }

  /* List::reserve() with size > capacity */
  {
    auto l = List(32);

    l.reserve(127);
    BOOST_CHECK(l.capacity() == 127);
  }
}

BOOST_AUTO_TEST_CASE(resize) {
  /* List::resize() with size < capacity */
  {
    auto l = List(128);

    l.resize(127);
    BOOST_CHECK(l.size() == 127);
    BOOST_CHECK(l.capacity() == 127);
    l.resize(42);
    BOOST_CHECK(l.size() == 42);
    BOOST_CHECK(l.capacity() == 42);
    l.resize(0);
    BOOST_CHECK(l.size() == 0);
    BOOST_CHECK(l.capacity() == 0);
  }

  /* List::resize() with size > capacity */
  {
    auto l = List(32);

    l.resize(127);
    BOOST_CHECK(l.size() == 127);
    BOOST_CHECK(l.capacity() == 127);
  }

  /* List::resize() with size == capacity */
  {
    auto l = List(32);

    auto const d_ptr = l.data();

    l.resize(32);
    BOOST_CHECK(l.size() == 32);
    BOOST_CHECK(l.capacity() == 32);
    BOOST_CHECK(d_ptr == l.data());
  }
}

BOOST_AUTO_TEST_CASE(operator_brackets) {
  /* List::operator[](size_type) */
  {
    auto l = List(32);
    std::iota(l.begin(), l.end(), 0);

    for (List::size_type i = 0; i < l.size(); i++) {
      BOOST_CHECK(l[i] == i);
    }
  }

  /* List::operator[](size_type) const */
  {
    auto i = List(32);
    std::iota(i.begin(), i.end(), 0);
    auto const l = i;

    for (List::size_type i = 0; i < l.size(); i++) {
      BOOST_CHECK(l[i] == i);
    }
  }
}

BOOST_AUTO_TEST_CASE(comparison) {
  /* List::operator{!,=}=(List const&), true */
  {
    auto l = List(31);
    std::iota(l.begin(), l.end(), 0);
    auto m = List(31);
    std::iota(m.begin(), m.end(), 0);

    BOOST_CHECK(l == m);
    BOOST_CHECK(not(l != m));
  }

  /* List::operator{!,=}=(List const&), wrong size */
  {
    auto l = List(31);
    std::iota(l.begin(), l.end(), 0);
    auto m = List(12);
    std::iota(m.begin(), m.end(), 0);

    BOOST_CHECK(not(l == m));
    BOOST_CHECK(l != m);
  }

  /* List::operator{!,=}=(List const&), wrong values */
  {
    auto l = List(31);
    std::iota(l.begin(), l.end(), 0);
    auto m = List(31);
    std::iota(m.begin(), m.end(), 0);
    m[11] = 0;

    BOOST_CHECK(not(l == m));
    BOOST_CHECK(l != m);
  }
}

BOOST_AUTO_TEST_CASE(push_back) {
  /* List::push_back(value_type const&) */
  {
    auto l = List(32);
    std::iota(l.begin(), l.end(), 0);

    for (int i = 32; i < 35; i++) {
      l.push_back(i);
    }

    BOOST_CHECK(l.size() == 35);

    for (int i = 0; i < l.size(); i++) {
      BOOST_CHECK(l[i] == i);
    }
  }
}

BOOST_AUTO_TEST_CASE(emplace_back) {
  /* List::emplace_back(...) */
  {
    auto l = List(31);

    l.emplace_back(31);

    BOOST_CHECK(l.size() == 32);
    BOOST_CHECK(l.back() == 31);
  }
}

BOOST_AUTO_TEST_CASE(erase) {
  /* List::erase(iterator first, iterator last) front */
  {
    auto l = List(32);
    std::iota(l.begin(), l.end(), 0);

    auto r = l.erase(l.begin(), l.begin() + 5);
    BOOST_CHECK(l.front() == 5);
    BOOST_CHECK(l.size() == (32 - 5));
    BOOST_CHECK(*r == 31);
    BOOST_CHECK(l.back() == 31);
  }

  /*  List::erase(iterator first, iterator last) back */
  {
    auto l = List(32);
    std::iota(l.begin(), l.end(), 0);

    auto r = l.erase(l.end() - 4, l.end());
    BOOST_CHECK(*r == (31 - 4));
    BOOST_CHECK(l.size() == (32 - 4));
    BOOST_CHECK(l.front() == 0);
    BOOST_CHECK(l.back() == (31 - 4));
  }

  /*  List::erase(iterator first, iterator last) middle */
  {
    auto l = List(32);
    std::iota(l.begin(), l.end(), 0);

    auto r = l.erase(l.begin() + 2, l.end() - 2);
    BOOST_CHECK(*r == 31);
    BOOST_CHECK(l.size() == 4);
    BOOST_CHECK(l.front() == 0);
    BOOST_CHECK(l.back() == 31);
  }

  /* List::erase(iterator first, iterator last) whole list */
  {
    auto l = List(32);

    auto r = l.erase(l.begin(), l.end());
    BOOST_CHECK(l.empty());
    BOOST_CHECK(r = l.begin());
  }
}

BOOST_AUTO_TEST_CASE(shrink_to_fit) {
  /* List::shrink_to_fit() size == capacity */
  {
    auto l = List(32);
    auto const d_ptr = l.data();

    l.shrink_to_fit();
    BOOST_CHECK(l.capacity() == 32);
    BOOST_CHECK(l.size() == 32);
    BOOST_CHECK(l.data() == d_ptr);
  }

  /* List::shrink_to_fit() size < capacity */
  {
    auto l = List(32);
    l.reserve(128);

    l.shrink_to_fit();
    BOOST_CHECK(l.capacity() == 32);
    BOOST_CHECK(l.size() == 32);
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  auto l = List(32);
  std::iota(l.begin(), l.end(), 0);

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << l;

  boost::archive::text_iarchive in_ar(stream);
  auto m = List();
  in_ar >> m;

  BOOST_CHECK(l.size() == m.size());
  BOOST_CHECK(std::equal(l.begin(), l.end(), m.begin()));
}

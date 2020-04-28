/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE Variant test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/Variant.hpp"
using namespace ScriptInterface;

#include "script_interface/get_value.hpp"

BOOST_AUTO_TEST_CASE(is_a) {
  BOOST_CHECK(is_type<None>(Variant(None{})));
  BOOST_CHECK(is_type<bool>(Variant(bool{})));
  BOOST_CHECK(is_type<int>(Variant(int{})));
}

BOOST_AUTO_TEST_CASE(none_is_default) { BOOST_CHECK(is_type<None>(Variant{})); }

BOOST_AUTO_TEST_CASE(make_from_args_test) {
  struct C {
    int i;

    C() = default;
    C(int i, double, std::string s) : i{i} { s.clear(); }
  };

  {
    VariantMap vals;

    auto c = make_from_args<C>(vals);
    c.i = 5;

    BOOST_CHECK(5 == c.i);
  }

  {
    VariantMap vals{{"a", 1.3}, {"b", 5}, {"c", std::string("c")}};

    auto c = make_from_args<C, int, double, std::string>(vals, "b", "a", "c");

    BOOST_CHECK(5 == c.i);
  }

  /* Missing argument */
  {
    VariantMap vals{{"a", 1.3}, {"b", 5}, {"c", std::string("c")}};

    BOOST_CHECK_THROW(
        (make_from_args<C, int, double, std::string>(vals, "b", "a", "d")),
        std::out_of_range);
  }

  /* Wrong type */
  {
    VariantMap vals{{"a", 1.3}, {"b", 5.0}, {"c", std::string("c")}};

    BOOST_CHECK_THROW(
        (make_from_args<C, int, double, std::string>(vals, "b", "a", "c")),
        std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(make_shared_from_args_test) {
  struct C {
    int i;

    C() = default;
    C(int i, double, std::string s) : i{i} { s.clear(); }
  };

  {
    VariantMap vals;

    auto c = make_shared_from_args<C>(vals);
    c->i = 5;
  }

  {
    VariantMap vals{{"a", 1.3}, {"b", 5}, {"c", std::string("c")}};

    auto c =
        make_shared_from_args<C, int, double, std::string>(vals, "b", "a", "c");

    BOOST_CHECK(5 == c->i);
  }

  /* Missing argument */
  {
    VariantMap vals{{"a", 1.3}, {"b", 5}, {"c", std::string()}};

    BOOST_CHECK_THROW((make_shared_from_args<C, int, double, std::string>(
                          vals, "b", "a", "d")),
                      std::out_of_range);
  }

  /* Wrong type */
  {
    VariantMap vals{{"a", 1.3}, {"b", 5.0}, {"c", std::string("c")}};

    BOOST_CHECK_THROW((make_shared_from_args<C, int, double, std::string>(
                          vals, "b", "a", "c")),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(pack_pair_test) {
  using ScriptInterface::pack_pair;
  using ScriptInterface::unpack_pair;

  const std::pair<int, double> p{4, 3.14};

  BOOST_CHECK((p == unpack_pair<int, double>(pack_pair(p))));
}

BOOST_AUTO_TEST_CASE(pack_map_test) {
  using ScriptInterface::pack_map;
  using ScriptInterface::unpack_map;

  const std::unordered_map<double, int> map{{3.3, 2}, {4.2, 9}};

  BOOST_CHECK((map == unpack_map<double, int>(pack_map(map))));
}

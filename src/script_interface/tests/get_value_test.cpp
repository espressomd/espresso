/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE ScriptInterface::get_value test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/get_value.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>

BOOST_AUTO_TEST_CASE(default_case) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  auto const s = std::string{"Abc"};
  auto const v = Variant(s);

  BOOST_CHECK(s == get_value<std::string>(v));
}

BOOST_AUTO_TEST_CASE(conversions) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;
  using ScriptInterface::detail::allow_conversion;

  static_assert(allow_conversion<int, int>::value, "");
  static_assert(allow_conversion<double, int>::value, "");
  static_assert(not allow_conversion<int, double>::value, "");

  BOOST_CHECK_EQUAL(3.1415, get_value<double>(3.1415));
  BOOST_CHECK_EQUAL(double(3), get_value<double>(3));
}

BOOST_AUTO_TEST_CASE(static_vector) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  /* From same type */
  {
    Variant v = std::vector<Variant>({1., 2., 3.});
    auto const expected = Utils::Vector3d{1., 2., 3.};
    BOOST_CHECK(get_value<Utils::Vector3d>(v) == expected);
  }

  /* Conversion applied */
  {
    Variant v = std::vector<Variant>({1, 2, 3});
    auto const expected = Utils::Vector3d{1, 2, 3};
    BOOST_CHECK(get_value<Utils::Vector3d>(v) == expected);
  }
}

BOOST_AUTO_TEST_CASE(heap_vector) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  /* From same type */
  {
    Variant v = std::vector<Variant>({1., 2., 3.});
    auto const expected = std::vector<double>{1., 2., 3.};
    BOOST_CHECK(get_value<std::vector<double>>(v) == expected);
  }

  /* Conversion applied */
  {
    Variant v = std::vector<Variant>({1, 2, 3});
    auto const expected = std::vector<double>{1, 2, 3};
    BOOST_CHECK(get_value<std::vector<double>>(v) == expected);
  }
}

BOOST_AUTO_TEST_CASE(get_value_from_map) {
  using ScriptInterface::get_value;
  using ScriptInterface::get_value_or;
  using ScriptInterface::Variant;
  using ScriptInterface::VariantMap;

  VariantMap map{{"a", 13}, {"e", 3.1}, {"f", "s"}};

  BOOST_CHECK(3.1 == get_value<double>(map, "e"));
  BOOST_CHECK(13 == get_value_or(map, "a", -1));
  BOOST_CHECK(-1 == get_value_or(map, "nope", -1));
}

BOOST_AUTO_TEST_CASE(get_map_value) {
  using ScriptInterface::get_map;
  using ScriptInterface::Variant;

  std::unordered_map<int, Variant> const map_variant{{1, 1.5}, {2, 2.5}};
  std::unordered_map<int, double> const map = get_map<int, double>(map_variant);
  BOOST_CHECK_EQUAL(map.at(1), 1.5);
  BOOST_CHECK_EQUAL(map.at(2), 2.5);

  std::unordered_map<int, Variant> const mixed{{1, 1}, {2, std::string("2")}};
  BOOST_CHECK_THROW((get_map<int, double>(mixed)), std::exception);
}

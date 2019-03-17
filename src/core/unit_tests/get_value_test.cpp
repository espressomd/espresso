/*
  Copyright (C) 2017-2018 The ESPResSo project

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

#define BOOST_TEST_MODULE ScriptInterface::get_value test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "get_value.hpp"

BOOST_AUTO_TEST_CASE(default_case) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  auto const s = std::string{"Abc"};
  auto const v = Variant(s);

  BOOST_CHECK(s == get_value<std::string>(v));
}

BOOST_AUTO_TEST_CASE(static_vector) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  Variant v = std::vector<double>({1., 2., 3.});
  Vector3d r = get_value<Vector3d>(v);
  BOOST_CHECK((Vector3d{1., 2., 3.} == r));
}

BOOST_AUTO_TEST_CASE(empty_vector) { // NOLINT
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  Variant v = std::vector<Variant>{};        // NOLINT
  auto vec = get_value<std::vector<int>>(v); // NOLINT

  BOOST_CHECK(std::vector<int>{} == vec);
}

// BOOST_AUTO_TEST_CASE(script_object) {
//   using ScriptInterface::get_value;
//   using ScriptInterface::Variant;
//   using ScriptInterface::ObjectId;
//   using so_ptr = std::shared_ptr<ScriptInterface::ScriptInterfaceBase>;

//   {
//     auto v = Variant{ObjectId{}};
//     BOOST_CHECK(nullptr == get_value<so_ptr>(v));
//   }
// }

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

/*
 * Copyright (C) 2018-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE Utils::strcat_alloc test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/strcat_alloc.hpp>
using Utils::strcat_alloc;

BOOST_AUTO_TEST_CASE(empty_left) {
  const char *right = "right";
  auto res = strcat_alloc(nullptr, right);

  BOOST_REQUIRE(res);
  BOOST_REQUIRE(strlen(res) == strlen(right));
  BOOST_CHECK(0 == strcmp(res, right));
  free(res);
}

BOOST_AUTO_TEST_CASE(empty_right) {
  char left[] = "left";
  auto res = strcat_alloc(left, nullptr);

  BOOST_REQUIRE(res == left);
  BOOST_REQUIRE(strlen(res) == strlen(left));
  BOOST_CHECK(0 == strcmp(res, left));
}

BOOST_AUTO_TEST_CASE(nonempty) {
  const char *left = "left";
  const char *right = "right";
  auto *str = static_cast<char *>(Utils::malloc(5));
  strncpy(str, left, 5);

  auto res = strcat_alloc(str, right);

  BOOST_REQUIRE(res);
  BOOST_REQUIRE(strlen(res) == (strlen(left) + strlen(right)));
  BOOST_CHECK(0 == strcmp(res, "leftright"));
  free(res);
}

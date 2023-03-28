/*
 * Copyright (C) 2020-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE ScriptInterface::Exception test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <script_interface/Exception.hpp>

#include <type_traits>

BOOST_AUTO_TEST_CASE(ctor) {
  /* Exception can be formed from a string */
  static_assert(
      std::is_constructible_v<ScriptInterface::Exception, std::string>);
  /* Exception can be formed from a char constant */
  static_assert(
      std::is_constructible_v<ScriptInterface::Exception, const char *>);
  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(what_) {
  /* The what method returns the message from construction */
  auto const msg = std::string("error message");
  BOOST_CHECK_EQUAL(msg, ScriptInterface::Exception(msg).what());
}

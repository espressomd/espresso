/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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

/** \file RuntimeError_test.cpp Unit tests for the ErrorHandling::RuntimeError
 * class.
 *
*/

#include <sstream>
#include <string>

#define BOOST_TEST_MODULE RuntimeError test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "../RuntimeError.hpp"

using std::string;
using ErrorHandling::RuntimeError;

/** Check constructor and getters */
BOOST_AUTO_TEST_CASE(values) {
  RuntimeError::ErrorLevel level = RuntimeError::ErrorLevel::WARNING;
  string what("Test error");
  int who(5);
  string function("Test_function");
  string file("Test_file.cpp");
  int line(42);
  RuntimeError err(level, who, what, function, file, line);

  BOOST_CHECK(level == err.level());
  BOOST_CHECK(what == err.what());
  BOOST_CHECK(who == err.who());
  BOOST_CHECK(function == err.function());
  BOOST_CHECK(file == err.file());
  BOOST_CHECK(line == err.line());
}

/** Check copy ctor */
BOOST_AUTO_TEST_CASE(def_ctor_and_assignment) {
  RuntimeError::ErrorLevel level = RuntimeError::ErrorLevel::WARNING;
  string what("Test error");
  int who(5);
  string function("Test_function");
  string file("Test_file.cpp");
  int line(42);
  RuntimeError err(level, who, what, function, file, line);

  /** Copy ctor */
  RuntimeError err2(err); // NOLINT

  BOOST_CHECK(level == err2.level());
  BOOST_CHECK(what == err2.what());
  BOOST_CHECK(who == err2.who());
  BOOST_CHECK(function == err2.function());
  BOOST_CHECK(file == err2.file());
  BOOST_CHECK(line == err2.line());
}

/** Check the serialization */
BOOST_AUTO_TEST_CASE(serialization) {
  std::stringstream ss;
  boost::archive::text_oarchive oa(ss);

  RuntimeError::ErrorLevel level = RuntimeError::ErrorLevel::WARNING;
  string what("Test error");
  int who(7);
  string function("Test_function");
  string file("Test_file.cpp");
  int line(21);
  RuntimeError err(level, who, what, function, file, line);

  /** Serialize to string stream */
  oa << err;

  boost::archive::text_iarchive ia(ss);
  RuntimeError err2;

  /** Deserialize into empty instance */
  ia >> err2;

  /** Check that the result is equal to the original instance */
  BOOST_CHECK((err.level() == err2.level()) && (err.who() == err2.who()) &&
              (err.what() == err2.what()) &&
              (err.function() == err2.function()) &&
              (err.file() == err2.file()) && (err.line() == err2.line()));
}

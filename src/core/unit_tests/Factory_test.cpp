/*
  Copyright (C) 2015,2016 The ESPResSo project
  
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

/** \file Factory_test.cpp Unit tests for the Utils::Factory class.
 *
*/

#define BOOST_TEST_MODULE Factory test
#include <boost/test/included/unit_test.hpp>

#include "utils/Factory.hpp"
#include "testing_utils/HelloWorld.hpp"

using namespace Testing;
using namespace std;

using BaseType = HelloWorld;
using Factory = Utils::Factory<BaseType>;

BOOST_AUTO_TEST_CASE(test_factory) {
  Factory::register_new("HelloWorld", Factory::builder<HelloWorld>);
  Factory::register_new("foo", Factory::builder<HelloWorld>);
  Factory::register_new("bar", Factory::builder<HelloWorld>);

  HelloWorld *p = Factory::make("HelloWorld");

  BOOST_CHECK(Factory::has_builder("HelloWorld"));
  BOOST_CHECK(p != nullptr);
  BOOST_CHECK((HelloWorld::constructed_once && p->instance) == true);
  delete p;

  BOOST_CHECK(Factory::has_builder("foo") and Factory::has_builder("bar"));
  p = Factory::make("bar");
  BOOST_CHECK(p != nullptr);
  BOOST_CHECK((HelloWorld::constructed_once && p->instance) == true);
  delete p;  
}

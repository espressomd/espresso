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

#include "../utils/Factory.hpp"

using namespace std;

class HelloWorld {
 public:
  HelloWorld() {
    constructed = true;
    instance = true;
  }
  static bool constructed;
  bool instance;
};

bool HelloWorld::constructed = false;

using BaseType = HelloWorld;
using Factory = Utils::Factory<BaseType>;

BOOST_AUTO_TEST_CASE(test_factory) {
  Factory::Instance().register_new("HelloWorld", Factory::builder<HelloWorld>);

  HelloWorld *p = Factory::Instance().make("HelloWorld");

  BOOST_CHECK(Factory::Instance().has_builder("HelloWorld"));
  BOOST_CHECK(p != nullptr);
  BOOST_CHECK((HelloWorld::constructed && p->instance) == true);
  
  delete p;
}

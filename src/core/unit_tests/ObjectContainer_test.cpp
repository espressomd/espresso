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

/** \file ObjectContainer_test.cpp Unit tests for the Utils::ObjectContainer class.
 *
*/

#define BOOST_TEST_MODULE Factory test
#include <boost/test/included/unit_test.hpp>

#include "ObjectContainer.hpp"

BOOST_AUTO_TEST_CASE(test_ObjectContainer) {
  ObjectContainer<int> container;

  int a = 1, b = 2, c = 3, d = 4, e = 5;

  BOOST_CHECK(container.add(a) == 0);
  BOOST_CHECK(container.add(b) == 1);
  BOOST_CHECK(container.add(c) == 2);
  BOOST_CHECK(container.add(d) == 3);
  BOOST_CHECK(container.add(e) == 4);

  container.remove(1);
  BOOST_CHECK(container.add(b) == 1);

  container.remove(4);
  BOOST_CHECK(container.add(e) == 4);

  container.remove(2);
  container.remove(3);
  BOOST_CHECK(container.add(c) == 2);
  BOOST_CHECK(container.add(d) == 3);

  container.remove(3);
  container.remove(4);
  container.remove(2);
  container.remove(0);
  container.remove(1);
  
  BOOST_CHECK(container.add(a) == 0);
  BOOST_CHECK(container.add(b) == 1);
  BOOST_CHECK(container.add(c) == 2);
  BOOST_CHECK(container.add(d) == 3);
  BOOST_CHECK(container.add(e) == 4);

  BOOST_CHECK(container[0] == a);
  BOOST_CHECK(container[4] == e);
}

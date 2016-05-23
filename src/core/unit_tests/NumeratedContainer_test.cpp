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

/** \file NumeratedContainer_test.cpp Unit tests for the Utils::NumeratedContainer class.
 *
*/

#define BOOST_TEST_MODULE NumeratedContainerTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <exception>

#include "utils/NumeratedContainer.hpp"

using namespace Utils;

BOOST_AUTO_TEST_CASE(test_NumeratedContainer) {
  NumeratedContainer<int> container;

  int a = 18, b = 7, c = 3, d = 35325, e = 588;

  /** Check adding elements */
  BOOST_CHECK(container.add(a) == 0);
  BOOST_CHECK(container.add(b) == 1);
  BOOST_CHECK(container.add(c) == 2);
  BOOST_CHECK(container.add(d) == 3);
  BOOST_CHECK(container.add(e) == 4);

  /** Check removing and id reuse */
  container.remove(1);
  BOOST_CHECK(container.add(b) == 1);

  container.remove(4);
  BOOST_CHECK(container.add(e) == 4);

  container.remove(2);
  container.remove(3);
  BOOST_CHECK(container.add(c) == 2);
  BOOST_CHECK(container.add(d) == 3);

  /** Check out-of-order remove */
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

  /** Check values */
  BOOST_CHECK(container[0] == a);
  BOOST_CHECK(container[4] == e);

  /** Check that out-of-bounds check works */
  BOOST_CHECK_THROW(container[5], std::out_of_range);
}

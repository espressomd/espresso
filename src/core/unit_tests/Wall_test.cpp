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

#include <limits>
#include <cmath>

#define BOOST_TEST_MODULE Wall test
#include <boost/test/included/unit_test.hpp>

#include "shapes/Wall.hpp"

BOOST_AUTO_TEST_CASE(name) {
  Shapes::Wall w;

  BOOST_CHECK(w.name() == "Shape::Wall");
}

/** Tests the generic setters and getters. */
BOOST_AUTO_TEST_CASE(script_interface) {
  Shapes::Wall w;    
  std::map<std::string, Variant> test_parameters;

  /* Wall has parameters normal and dist */
  std::vector<double> normal{0., 1., 0.};
  test_parameters["normal"] = normal;
  test_parameters["dist"] = 5.0;
  
  w.set_parameters(test_parameters);

  BOOST_CHECK(boost::get<double>(w.get_parameter("dist")) == 5.0);

  /* Get value as vector */
  auto ret_normal = boost::get<std::vector<double> >(w.get_parameter("normal"));
  /* Check if equal up to numerical precision (wall has normalized the vector,
     so it need not be bitwise equal). */
  for(int i = 0; i < 3; i++) {
    BOOST_CHECK(std::abs(ret_normal[i] - normal[i]) <= std::numeric_limits<double>::epsilon());
  }
}

/* @TODO: Functional unit test of the distance function */

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

/** \file
 * Unit tests for the Particle struct.
 *
 */

#define BOOST_TEST_MODULE Particle test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "serialization/Particle.hpp"

BOOST_AUTO_TEST_CASE(comparison) {
  {
    Particle p, q;

    p.identity() = 1;
    q.identity() = 2;

    BOOST_CHECK(p != q);
    BOOST_CHECK(not(p == q));
  }

  {
    Particle p, q;

    p.identity() = 2;
    q.identity() = 2;

    BOOST_CHECK(not(p != q));
    BOOST_CHECK(p == q);
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  auto p = Particle();

  Utils::List<int> bl = {1, 2, 3, 4};
  Utils::List<int> el = {5, 6, 7, 8};

  p.p.identity = 15;
  p.bl = bl;
#ifdef EXCLUSIONS
  p.el = el;
#endif

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << p;

  boost::archive::text_iarchive in_ar(stream);
  auto q = Particle();
  in_ar >> q;

  BOOST_CHECK(q.p.identity == p.p.identity);
  BOOST_CHECK(q.bl == bl);

#ifdef EXCLUSIONS
  BOOST_CHECK(q.el == el);
#endif
}

/*
   Copyright (C) 2019 The ESPResSo project

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

#define BOOST_TEST_MODULE Algorithm::periodic_fold test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "algorithm/periodic_fold.hpp"

BOOST_AUTO_TEST_CASE(with_image_count) {
  using Algorithm::periodic_fold;

  /* Inside */
  {
    auto const x = 5.;
    auto const i = 9;
    auto const res = periodic_fold(x, i, 10.);
    BOOST_CHECK_EQUAL(res.first, x);
    BOOST_CHECK_EQUAL(res.second, i);
  }

  /* Outside left */
  {
    auto const x = -5.;
    auto const i = 9;
    auto const box = 10.;
    auto const res = periodic_fold(x, i, box);
    BOOST_CHECK_EQUAL(res.first - box, x);
    BOOST_CHECK_EQUAL(res.second, i - 1);
  }

  /* Left boundary */
  {
    auto const x = 0.;
    auto const i = 9;
    auto const box = 10.;
    auto const res = periodic_fold(x, i, box);
    BOOST_CHECK_EQUAL(res.first, x);
    BOOST_CHECK_EQUAL(res.second, i);
  }

  /* Right boundary */
  {
    auto const x = 10.;
    auto const i = 9;
    auto const box = x;
    auto const res = periodic_fold(x, i, box);
    BOOST_CHECK_EQUAL(res.first, x - box);
    BOOST_CHECK_EQUAL(res.second, i + 1);
  }

  /* Pathological (NaN) */
  {
    auto const x = std::nan("");
    auto const i = 9;
    auto const box = 10.;
    auto const res = periodic_fold(x, i, box);
    BOOST_CHECK(std::isnan(res.first));
  }

  /* Overflow right */
  {
    auto const x =
        (100. * static_cast<double>(std::numeric_limits<int>::max()));
    auto const box = 10.;
    int const i = std::numeric_limits<int>::max() - 10;
    auto const res = periodic_fold(x, i, box);
    BOOST_CHECK_EQUAL(res.second, std::numeric_limits<int>::max());
  }

  /* Overflow left */
  {
    auto const x =
        (100. * static_cast<double>(std::numeric_limits<int>::min()));
    auto const box = 10.;
    int const i = std::numeric_limits<int>::min() + 10;
    auto const res = periodic_fold(x, i, box);
    BOOST_CHECK_EQUAL(res.second, std::numeric_limits<int>::min());
  }

  /* Corner left */
  {
    auto const x = std::nextafter(0., -1.);
    auto const box = 10.;
    auto const res = periodic_fold(x, 0, box);
    BOOST_CHECK(res.first >= 0.);
    BOOST_CHECK(res.first < box);
    BOOST_CHECK(std::abs(res.first - x + res.second * box) <=
                std::numeric_limits<double>::epsilon());
  }

  /* Corner right */
  {
    auto const box = 10.;
    auto const x = std::nextafter(box, 2 * box);
    auto const res = periodic_fold(x, 0, box);
    BOOST_CHECK(res.first >= 0.);
    BOOST_CHECK(res.first < box);
    BOOST_CHECK(std::abs(res.first - x + res.second * box) <=
                std::numeric_limits<double>::epsilon());
  }
}

BOOST_AUTO_TEST_CASE(without_image_count) {
  using Algorithm::periodic_fold;

  /* Inside */
  {
    auto const x = 5.;
    auto const res = periodic_fold(x, 10.);
    BOOST_CHECK_EQUAL(res, x);
  }

  /* Outside left */
  {
    auto const x = -5.;
    auto const box = 10.;
    auto const res = periodic_fold(x, box);
    BOOST_CHECK_EQUAL(res - box, x);
  }

  /* Left boundary */
  {
    auto const x = 0.;
    auto const box = 10.;
    auto const res = periodic_fold(x, box);
    BOOST_CHECK_EQUAL(res, x);
  }

  /* Right boundary */
  {
    auto const x = 10.;
    auto const box = x;
    auto const res = periodic_fold(x, box);
    BOOST_CHECK_EQUAL(res, x - box);
  }

  /* Pathological (NaN value) */
  {
    auto const x = std::nan("");
    auto const box = 10.;
    auto const res = periodic_fold(x, box);
    BOOST_CHECK(std::isnan(res));
  }

  /* Pathological (NaN box) */
  {
    auto const x = 5.;
    auto const box = std::nan("");
    auto const res = periodic_fold(x, box);
    BOOST_CHECK(std::isnan(res));
  }

  /* Corner left */
  {
    auto const x = std::nextafter(0., -1.);
    auto const box = 10.;
    auto const res = periodic_fold(x, box);
    BOOST_CHECK(res >= 0.);
    BOOST_CHECK(res < box);
  }

  /* Corner right */
  {
    auto const box = 10.;
    auto const x = std::nextafter(box, 2 * box);
    auto const res = periodic_fold(x, box);
    BOOST_CHECK(res >= 0.);
    BOOST_CHECK(res < box);
  }
}

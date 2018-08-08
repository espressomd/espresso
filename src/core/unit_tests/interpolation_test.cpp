/*
  Copyright (C) 2016 The ESPResSo project

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

#define BOOST_TEST_MODULE Utils::interpolate test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/interpolation/bspline_3d.hpp"
using Utils::Interpolation::bspline_3d;
using Utils::Interpolation::bspline_3d_accumulate;
using Utils::Interpolation::detail::pos_shift;
using Utils::Interpolation::detail::ll_and_dist;

#include "common/gaussian.hpp"

#include <limits>

BOOST_AUTO_TEST_CASE(pos_shift_test) {
  BOOST_CHECK(pos_shift<2>() == 0.0);  
  BOOST_CHECK(pos_shift<3>() == 0.5);
}

BOOST_AUTO_TEST_CASE(ll_and_dist_test) {
  std::array<double, 3> dist;
  std::array<int, 3> ll;

  std::tie(ll, dist) = ll_and_dist<3>(
      Vector3d{.1, .2, .3}, Vector3d{0.5, 0.5, 0.5}, Vector3d{-1., -3., -2.});

  /* Pos with offset is {1.1, 3.2, 2.3} */
  /* nmp is {2, 6, 5} @ {1.0, 3.0, 2.5} */
  /* distance in lattice units is +.1/.5, +.2/.5, -.2/.5 = {.2, .4, -.4} */
  BOOST_CHECK_CLOSE(dist[0], .2, 1e-13);
  BOOST_CHECK_CLOSE(dist[1], .4, 1e-13);
  BOOST_CHECK_CLOSE(dist[2], -.4, 1e-13);

  BOOST_CHECK(ll[0] == 1);
  BOOST_CHECK(ll[1] == 5);
  BOOST_CHECK(ll[2] == 4);
}

BOOST_AUTO_TEST_CASE(number_of_points) {
  int count = 0;
  auto counter = [&count](const std::array<int, 3> &, double) { count++; };

  bspline_3d<5>({}, counter, {1., 1., 1.}, {});

  BOOST_CHECK(5 * 5 * 5 == count);
}

BOOST_AUTO_TEST_CASE(sum_of_weights) {
  double sum = 0.0;
  auto summer = [&sum](const std::array<int, 3> &, double w) { sum += w; };

  bspline_3d<5>({.1, .2, .3}, summer, {4, 5, 6}, {7, 8, 9});

  BOOST_CHECK_CLOSE(sum, 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(nearest_point) {
  std::array<int, 3> nmp;
  double weight;
  auto save_ind = [&nmp, &weight](const std::array<int, 3> &ind, double w) {
    nmp = ind;
    weight = w;
  };

  bspline_3d<1>({.1, .2, .3}, save_ind, {0.5, 0.5, 0.5}, {});

  BOOST_CHECK((std::array<int, 3>{{0, 0, 1}} == nmp));
  BOOST_CHECK_CLOSE(weight, 1., std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(interpolation_points_3) {
  std::vector<std::array<int, 3>> int_points;

  auto save_ind = [&int_points](const std::array<int, 3> &ind, double w) {
    int_points.push_back(ind);
  };

  bspline_3d<3>({5., 6., 7.}, save_ind, /* grid spacing */ {1., 2., 3.},
                /* offset */ {10., 0., 15.});

  /* pos - offset = {-5., 6., -8} */
  /* nmp = {-5, 3, -3 } @ pos {-5., 6., -9.} */
  /* minus order / 2 (= 1) = {-6, 2, -4} */
  std::array<int, 3> lower_left = {{-5, 2, -3}};
  /* upper limit = {-6, 2, -4} + order = {-3, 5, -1} */

  auto it = int_points.begin();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        auto const expected = std::array<int, 3>{
		{lower_left[0] + i, lower_left[1] + j, lower_left[2] + k}};
        BOOST_CHECK(*it == expected);
        ++it;
      }
}

BOOST_AUTO_TEST_CASE(interpolation_points_2) {
  std::vector<std::array<int, 3>> int_points;

  auto save_ind = [&int_points](const std::array<int, 3> &ind, double w) {
    int_points.push_back(ind);
  };

  bspline_3d<2>({5., 6., 7.}, save_ind, /* grid spacing */ {1., 2., 3.},
                /* offset */ {10., 0., 15.});

  /* pos - offset = {-5., 6., -8} */
  /* shited pos = {-4.5, 6.5, -7.5} */
  /* nmp = {-4, 4, -2 } @ pos {-4., 8., -6.} */
  /* ll = nmp - order / 2 (= 1) = {-5, 3, -3} */
  std::array<int, 3> lower_left = {{-5, 3, -2}};
  /* upper limit = {-5, 2, -3} + order = {-3, 4, -1} */

  auto it = int_points.begin();
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        auto const expected = std::array<int, 3>{
		{lower_left[0] + i, lower_left[1] + j, lower_left[2] + k}};
        BOOST_CHECK(*it == expected);
        ++it;
      }
}

BOOST_AUTO_TEST_CASE(interpolation_weights_2) {
  auto check_weight = [](const std::array<int, 3> &ind, double w) {
    auto const expected = Utils::bspline<2>(ind[0], 0.1 - 0.5) *
                          Utils::bspline<2>(ind[1], 0.2 - 0.5) *
                          Utils::bspline<2>(ind[2], 0.3 - 0.5);
    BOOST_CHECK_CLOSE(w, expected, 1e-12);
  };

  bspline_3d<2>({.1, .2, .3}, check_weight,
                /* grid spacing */ {1., 1., 1.},
                /* offset */ {0., 0., 0.});
}

BOOST_AUTO_TEST_CASE(interpolation_weights_3) {
  auto check_weight = [](const std::array<int, 3> &ind, double w) {
    auto const expected = Utils::bspline<3>(ind[0], 0.1) *
                          Utils::bspline<3>(ind[1], 0.2) *
                          Utils::bspline<3>(ind[2], 0.3);
    BOOST_CHECK_CLOSE(w, expected, 1e-12);
  };

  bspline_3d<3>({1.1, 1.2, 1.3}, check_weight,
                /* grid spacing */ {1., 1., 1.},
                /* offset */ {0., 0., 0.});
}

BOOST_AUTO_TEST_CASE(interpolation_accumulate) {
  auto const w_acc = bspline_3d_accumulate(
      {1.1, 1.2, 1.3}, [](const std::array<int, 3> &) { return 2.; },
      /* grid spacing */ {1., 1., 1.},
      /* offset */ {0., 0., 0.}, 7, 0.0);

  BOOST_CHECK_CLOSE(w_acc, 2., 1e-12);
}

BOOST_AUTO_TEST_CASE(interpolation_integration_test_odd) {
  const int order = 5;
  const Vector3d grid_spacing = {.1, .2, .3};
  const Vector3d origin = {-1., 2., -1.4};
  const int n_nodes = 10;

  auto const x0 = origin + 0.57 * 10. * grid_spacing;
  auto const sigma = 4.;

  auto const data = gaussian_field(n_nodes, grid_spacing, origin, x0, sigma);

  auto const p = Vector3d{-.4, 3.14, 0.1};
  auto const interpolated_value = bspline_3d_accumulate(
      p, [&data](const std::array<int, 3> &ind) { return data(ind); },
      grid_spacing, origin, order, 0.0);

  auto const exact_value = gaussian(p, x0, sigma);

  BOOST_CHECK_CLOSE(interpolated_value, exact_value, 1.);
}

BOOST_AUTO_TEST_CASE(interpolation_integration_test_even) {
  const int order = 6;
  const Vector3d grid_spacing = {.1, .2, .3};
  const Vector3d origin = {-1., 2., -1.4};
  const int n_nodes = 10;

  auto const x0 = origin + 0.57 * 10. * grid_spacing;
  auto const sigma = 4.;

  auto const data = gaussian_field(n_nodes, grid_spacing, origin, x0, sigma);

  auto const p = Vector3d{-.4, 3.14, 0.1};
  auto const interpolated_value = bspline_3d_accumulate(
      p, [&data](const std::array<int, 3> &ind) { return data(ind); },
      grid_spacing, origin, order, 0.0);

  auto const exact_value = gaussian(p, x0, sigma);

  BOOST_CHECK_CLOSE(interpolated_value, exact_value, 1.);
}

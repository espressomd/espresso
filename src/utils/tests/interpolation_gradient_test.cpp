/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE Utils::interpolate_gradient test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/interpolation/bspline_3d_gradient.hpp"
using Utils::Interpolation::bspline_3d_gradient;
using Utils::Interpolation::bspline_3d_gradient_accumulate;

#include "utils/math/gaussian.hpp"
#include "utils/raster.hpp"

BOOST_AUTO_TEST_CASE(number_of_points) {
  int count = 0;
  auto counter = [&count](const std::array<int, 3> &, const Utils::Vector3d &) {
    count++;
  };

  bspline_3d_gradient<5>({}, counter, {1., 1., 1.}, {});

  BOOST_CHECK(5 * 5 * 5 == count);
}

BOOST_AUTO_TEST_CASE(sum_of_weights) {
  const Utils::Vector3d grid_spacing{.4, .5, .6};

  Utils::Vector3d sum_const{};
  Utils::Vector3d sum_lin{};
  auto summer = [&sum_const, &sum_lin, &grid_spacing](
                    const std::array<int, 3> &ind, const Utils::Vector3d &w) {
    sum_const += w;
    sum_lin +=
        {w[0] * ind[0] * grid_spacing[0], w[1] * ind[1] * grid_spacing[1],
         w[2] * ind[2] * grid_spacing[2]};
  };

  bspline_3d_gradient<5>({.1, .2, .3}, summer, grid_spacing, {7., 8., 9.});

  for (int i = 0; i < 3; i++) {
    BOOST_CHECK_SMALL(sum_const[i], 1.e-12);
    BOOST_CHECK_CLOSE(sum_lin[i], 1., 1e-12);
  }
}

BOOST_AUTO_TEST_CASE(interpolation_weights_2) {
  auto check_weight = [](const std::array<int, 3> &ind,
                         const Utils::Vector3d &w) {
    auto const expected =
        Utils::Vector3d{Utils::bspline_d<2>(ind[0], 0.1 - 0.5) *
                            Utils::bspline<2>(ind[1], 0.2 - 0.5) *
                            Utils::bspline<2>(ind[2], 0.3 - 0.5),
                        Utils::bspline<2>(ind[0], 0.1 - 0.5) *
                            Utils::bspline_d<2>(ind[1], 0.2 - 0.5) *
                            Utils::bspline<2>(ind[2], 0.3 - 0.5),
                        Utils::bspline<2>(ind[0], 0.1 - 0.5) *
                            Utils::bspline<2>(ind[1], 0.2 - 0.5) *
                            Utils::bspline_d<2>(ind[2], 0.3 - 0.5)};

    for (int i = 0; i < 3; i++) {
      BOOST_CHECK_CLOSE(w[i], expected[i], 1e-12);
    }
  };

  bspline_3d_gradient<2>({.1, .2, .3}, check_weight,
                         /* grid spacing */ {1., 1., 1.},
                         /* offset */ {0., 0., 0.});
}

BOOST_AUTO_TEST_CASE(interpolation_weights_3) {
  auto check_weight = [](const std::array<int, 3> &ind,
                         const Utils::Vector3d &w) {
    auto const expected = Utils::Vector3d{
        Utils::bspline_d<3>(ind[0], 0.1) * Utils::bspline<3>(ind[1], 0.2) *
            Utils::bspline<3>(ind[2], 0.3),
        Utils::bspline<3>(ind[0], 0.1) * Utils::bspline_d<3>(ind[1], 0.2) *
            Utils::bspline<3>(ind[2], 0.3),
        Utils::bspline<3>(ind[0], 0.1) * Utils::bspline<3>(ind[1], 0.2) *
            Utils::bspline_d<3>(ind[2], 0.3)};

    for (int i = 0; i < 3; i++) {
      BOOST_CHECK_CLOSE(w[i], expected[i], 1e-12);
    }
  };

  bspline_3d_gradient<3>({1.1, 1.2, 1.3}, check_weight,
                         /* grid spacing */ {1., 1., 1.},
                         /* offset */ {0., 0., 0.});
}

BOOST_AUTO_TEST_CASE(interpolation_gradient_integration_test_odd) {
  constexpr int order = 5;
  const Utils::Vector3d grid_spacing = {.1, .2, .3};
  const Utils::Vector3d origin = {-1., 2., -1.4};
  const int n_nodes = 10;

  auto const x0 = origin + 0.57 * n_nodes * grid_spacing;
  auto const sigma = 4.;

  auto const data = Utils::raster<double>(
      origin, grid_spacing, Utils::Vector3i::broadcast(n_nodes),
      [&](auto x) { return gaussian(x, x0, sigma); });

  auto const p = Utils::Vector3d{-.4, 3.14, 0.1};
  auto const interpolated_value = bspline_3d_gradient_accumulate<order>(
      p, [&data](const std::array<int, 3> &ind) { return data(ind); },
      grid_spacing, origin, Utils::Vector3d{});

  auto const exact_value = del_gaussian(p, x0, sigma);

  BOOST_CHECK_SMALL((interpolated_value - exact_value).norm(), 1.e-4);
}

BOOST_AUTO_TEST_CASE(interpolation_gradient_vec_integration_test_odd) {
  constexpr int order = 5;
  const Utils::Vector3d grid_spacing = {.1, .2, .3};
  const Utils::Vector3d origin = {-1., 2., -1.4};
  const int n_nodes = 10;

  auto const a = origin + 0.37 * n_nodes * grid_spacing;
  Utils::Vector3d x0[2] = {0.12 * a, -3. * a};
  auto const sigma = Utils::Vector2d{2., 3.};

  boost::multi_array<Utils::Vector2d, 3> data(Utils::Vector3i{10, 10, 10});
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++) {
        auto const &h = grid_spacing;
        auto const x = origin + Utils::Vector3d{i * h[0], j * h[1], k * h[2]};
        data[i][j][k] = {gaussian(x, x0[0], sigma[0]),
                         gaussian(x, x0[1], sigma[1])};
      }

  auto const p = Utils::Vector3d{-.4, 3.14, 0.1};
  auto const interpolated_value = bspline_3d_gradient_accumulate<order>(
      p, [&data](const std::array<int, 3> &ind) { return data(ind); },
      grid_spacing, origin, Utils::Vector<Utils::Vector3d, 2>{});

  const Utils::Vector<Utils::Vector3d, 2> exact_value = {
      del_gaussian(p, x0[0], sigma[0]), del_gaussian(p, x0[1], sigma[1])};

  BOOST_CHECK_SMALL((interpolated_value[0] - exact_value[0]).norm(), 1.e-2);
  BOOST_CHECK_SMALL((interpolated_value[1] - exact_value[1]).norm(), 1.e-4);
}

BOOST_AUTO_TEST_CASE(interpolation_gradient_integration_test_even) {
  constexpr int order = 6;
  const Utils::Vector3d grid_spacing = {.1, .2, .3};
  const Utils::Vector3d origin = {-1., 2., -1.4};
  const int n_nodes = 10;

  auto const x0 = origin + 0.57 * n_nodes * grid_spacing;
  auto const sigma = 4.;

  auto const data = Utils::raster<double>(
      origin, grid_spacing, Utils::Vector3i::broadcast(n_nodes),
      [&](auto x) { return gaussian(x, x0, sigma); });

  auto const p = Utils::Vector3d{-.4, 3.14, 0.1};
  auto const interpolated_value = bspline_3d_gradient_accumulate<order>(
      p, [&data](const std::array<int, 3> &ind) { return data(ind); },
      grid_spacing, origin, Utils::Vector3d{});

  auto const exact_value = del_gaussian(p, x0, sigma);

  BOOST_CHECK_SMALL((interpolated_value - exact_value).norm(), 1.e-4);
}

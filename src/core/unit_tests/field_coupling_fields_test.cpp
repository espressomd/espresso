/*
Copyright (C) 2010-2018 The ESPResSo project

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

#define BOOST_TEST_MODULE AutoParameter test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "field_coupling/fields/AffineMap.hpp"
#include "field_coupling/fields/Constant.hpp"
#include "field_coupling/fields/Interpolated.hpp"
#include "field_coupling/fields/PlaneWave.hpp"
#include "field_coupling/fields/jacobian_type.hpp"

#include "utils/interpolation/bspline_3d.hpp"
#include "utils/interpolation/bspline_3d_gradient.hpp"
#include <utils/math/gaussian.hpp>
#include <utils/raster.hpp>

#include <limits>

using namespace FieldCoupling::Fields;

BOOST_AUTO_TEST_CASE(jacobian_type_test) {
  using FieldCoupling::Fields::detail::jacobian_type;
  using std::is_same;

  static_assert(is_same<jacobian_type<double, 1>, Utils::Vector3d>::value, "");
  static_assert(is_same<jacobian_type<double, 2>,
                        Utils::Vector<Utils::Vector3d, 2>>::value,
                "");
}

BOOST_AUTO_TEST_CASE(constant_scalar_field) {
  using Field = Constant<double, 1>;

  /* Types */
  {
    static_assert(std::is_same<Field::value_type, double>::value, "");
    static_assert(std::is_same<Field::jacobian_type, Utils::Vector3d>::value,
                  "");
  }

  /* ctor */
  {
    const double val = 1.23;
    Field field(val);

    BOOST_CHECK(val == field.value());
  }

  /* setter */
  {
    Field field(0.);

    const double val = 1.23;
    field.value() = val;

    BOOST_CHECK(val == field.value());
  }

  /* Field value */
  {
    Field field(5.);

    BOOST_CHECK(5. == field({1., 2., 3.}));
  }

  /* Gradient */
  {
    Field field(5.);

    BOOST_CHECK(
        (Utils::Vector3d{0.0, 0.0, 0.0} == field.jacobian({1., 2., 3.})));
  }
}

BOOST_AUTO_TEST_CASE(constant_vector_field) {
  using Field = Constant<double, 2>;

  /* Types */
  {
    static_assert(std::is_same<Field::value_type, Utils::Vector2d>::value, "");
    static_assert(std::is_same<Field::jacobian_type,
                               Utils::Vector<Utils::Vector3d, 2>>::value,
                  "");
  }

  /* ctor */
  {
    const Utils::Vector2d val = {1.23, 4.56};
    Field field(val);

    BOOST_CHECK(val == field.value());
  }

  /* setter */
  {
    Field field({});

    const Utils::Vector2d val = {1.23, 4.56};
    field.value() = val;

    BOOST_CHECK(val == field.value());
  }

  /* Field value */
  {
    const Utils::Vector2d field_val = {5., 6.};

    Field field(field_val);

    BOOST_CHECK(field_val == field({1., 2., 3.}));
  }

  /* Gradient */
  {
    const auto zero = Field::jacobian_type{};
    Field field({5., 6.});

    BOOST_CHECK(zero == field.jacobian({1., 2., 3.}));
  }
}

BOOST_AUTO_TEST_CASE(affine_scalar_field) {
  using Field = AffineMap<double, 1>;

  /* Types */
  {
    static_assert(std::is_same<Field::value_type, double>::value, "");
    static_assert(std::is_same<Field::jacobian_type,
                               Utils::Vector<Field::value_type, 3>>::value,
                  "");
  }

  /* ctor */
  {
    const Utils::Vector3d A = {1., 2., 3.};
    const double b = 4.;
    Field field(A, b);

    BOOST_CHECK(A == field.A());
    BOOST_CHECK(b == field.b());
  }

  /* setter */
  {
    Field field({}, {});

    const Utils::Vector3d A = {1., 2., 3.};
    const double b = 4.;
    field.A() = A;
    field.b() = b;

    BOOST_CHECK(A == field.A());
    BOOST_CHECK(b == field.b());
  }

  /* Field value */
  {
    const Utils::Vector3d A = {1., 2., 3.};
    const double b = 4.;
    Field field(A, b);

    const Utils::Vector3d x = {1., 2., 3.};

    BOOST_CHECK((A * x + b) == field(x));
  }

  /* Gradient */
  {
    const Utils::Vector3d A = {1., 2., 3.};
    const double b = 4.;
    Field field(A, b);

    BOOST_CHECK(A == field.jacobian({1., 2., 3.}));
  }
}

template <size_t N, size_t M, typename T>
using Matrix = Utils::Vector<Utils::Vector<T, M>, N>;

BOOST_AUTO_TEST_CASE(affine_vector_field) {
  using Field = AffineMap<double, 2>;

  /* Types */
  {
    static_assert(std::is_same<Field::value_type, Utils::Vector2d>::value, "");
    static_assert(
        std::is_same<Field::jacobian_type, Matrix<2, 3, double>>::value, "");
  }

  /* Field value unshifted */
  {
    const Utils::Vector<Utils::Vector3d, 2> A = {{1., 2., 3}, {4., 5., 6.}};
    const Utils::Vector2d b = {7., 8.};
    Field field(A, b);

    const Utils::Vector3d x = {1., 1., 1.};

    auto const res = field(x);

    BOOST_CHECK((A[0][0] + A[0][1] + A[0][2] + b[0]) == res[0]);
    BOOST_CHECK((A[1][0] + A[1][1] + A[1][2] + b[1]) == res[1]);
  }

  /* Gradient */
  {
    const Utils::Vector<Utils::Vector3d, 2> A = {{1., 2., 3}, {4., 5., 6.}};
    const Utils::Vector2d b = {7., 8.};
    Field field(A, b);

    BOOST_CHECK(A == field.jacobian({1., 2., 3.}));
  }
}

BOOST_AUTO_TEST_CASE(interpolated_scalar_field) {
  using Field = Interpolated<double, 1>;

  /* Types */
  {
    static_assert(std::is_same<Field::value_type, double>::value, "");
    static_assert(std::is_same<Field::jacobian_type, Utils::Vector3d>::value,
                  "");
  }

  /* Ctor */
  {
    boost::multi_array<double, 3> data(Utils::Vector3i{10, 11, 12});
    data[5][5][5] = -1. / 12.;

    const Utils::Vector3d grid_spacing = {.1, .2, .3};
    const Utils::Vector3d origin = {-1., 2., -3.};

    Field field(data, grid_spacing, origin);

    BOOST_CHECK(data == field.field_data());
    BOOST_CHECK(grid_spacing == field.grid_spacing());
    BOOST_CHECK(origin == field.origin());
  }

  /* field value */
  {
    using Utils::Interpolation::bspline_3d_accumulate;

    const Utils::Vector3d grid_spacing = {.1, .2, .3};
    const Utils::Vector3d origin = {-1., 2., -1.4};
    auto const n_nodes = Utils::Vector3i{10, 10, 10};

    auto const x0 = origin + 0.5 * 10 * grid_spacing;
    auto const sigma = 2.;

    auto const data =
        Utils::raster<double>(origin, grid_spacing, n_nodes,
                              [&](auto x) { return gaussian(x, x0, sigma); });

    Field field(data, grid_spacing, origin);

    auto const p = Utils::Vector3d{-.4, 3.14, 0.1};

    auto const interpolated_value = bspline_3d_accumulate<2>(
        p, [&data](const std::array<int, 3> &ind) { return data(ind); },
        grid_spacing, origin, 0.0);

    auto const field_value = field(p);

    BOOST_CHECK(std::abs(interpolated_value - field_value) <
                std::numeric_limits<double>::epsilon());
  }

  /* jacobian value */
  {
    using Utils::Interpolation::bspline_3d_gradient_accumulate;

    const Utils::Vector3d grid_spacing = {.1, .2, .3};
    const Utils::Vector3d origin = {-1., 2., -1.4};
    auto const n_nodes = 10;

    auto const x0 = origin + 0.57 * n_nodes * grid_spacing;
    auto const sigma = 2.;

    auto const data = Utils::raster<double>(
        origin, grid_spacing, Utils::Vector3i::broadcast(n_nodes),
        [&](auto x) { return gaussian(x, x0, sigma); });

    Field field(data, grid_spacing, origin);

    auto const p = Utils::Vector3d{-.4, 3.14, 0.1};

    auto const field_value = field.jacobian(p);

    auto const interpolated_value = bspline_3d_gradient_accumulate<2>(
        p, [&data](const std::array<int, 3> &ind) { return data(ind); },
        grid_spacing, origin, Utils::Vector3d{});

    BOOST_CHECK((interpolated_value - field_value).norm() <
                2 * std::numeric_limits<double>::epsilon());
  }
}

BOOST_AUTO_TEST_CASE(interpolated_vector_field) {
  using Field = Interpolated<double, 2>;

  /* Types */
  {
    static_assert(std::is_same<Field::value_type, Utils::Vector2d>::value, "");
    static_assert(std::is_same<Field::jacobian_type,
                               Utils::Vector<Utils::Vector3d, 2>>::value,
                  "");
  }

  /* field value */
  {
    using Utils::Interpolation::bspline_3d_accumulate;

    const Utils::Vector3d grid_spacing = {.1, .2, .3};
    const Utils::Vector3d origin = {-1., 2., -1.4};
    const int n_nodes = 10;

    auto const a = origin + 0.37 * n_nodes * grid_spacing;
    Utils::Vector3d x0[2] = {0.12 * a, -3. * a};
    auto const sigma = Utils::Vector2d{2., 3.};

    boost::multi_array<Utils::Vector2d, 3> data(
        Utils::Vector3i{n_nodes, n_nodes, n_nodes});
    for (int i = 0; i < n_nodes; i++)
      for (int j = 0; j < n_nodes; j++)
        for (int k = 0; k < n_nodes; k++) {
          auto const &h = grid_spacing;
          auto const x = origin + Utils::Vector3d{i * h[0], j * h[1], k * h[2]};
          data[i][j][k] = {gaussian(x, x0[0], sigma[0]),
                           gaussian(x, x0[1], sigma[1])};
        }

    Field field(data, grid_spacing, origin);

    auto const p = Utils::Vector3d{-.4, 3.14, 0.1};

    auto const field_value = field(p);

    auto const interpolated_value = bspline_3d_accumulate<2>(
        p, [&data](const std::array<int, 3> &ind) { return data(ind); },
        grid_spacing, origin, Utils::Vector2d{});

    BOOST_CHECK_SMALL((interpolated_value - field_value).norm(),
                      std::numeric_limits<double>::epsilon());
  }

  /* jacobian value */
  {
    using Utils::Interpolation::bspline_3d_gradient_accumulate;

    const Utils::Vector3d grid_spacing = {.1, .2, .3};
    const Utils::Vector3d origin = {-1., 2., -1.4};
    const int n_nodes = 10;

    auto const a = origin + 0.37 * n_nodes * grid_spacing;
    Utils::Vector3d x0[2] = {0.12 * a, -3. * a};
    auto const sigma = Utils::Vector2d{2., 3.};

    boost::multi_array<Utils::Vector2d, 3> data(
        Utils::Vector3i{n_nodes, n_nodes, n_nodes});
    for (int i = 0; i < n_nodes; i++)
      for (int j = 0; j < n_nodes; j++)
        for (int k = 0; k < n_nodes; k++) {
          auto const &h = grid_spacing;
          auto const x = origin + Utils::Vector3d{i * h[0], j * h[1], k * h[2]};
          data[i][j][k] = {gaussian(x, x0[0], sigma[0]),
                           gaussian(x, x0[1], sigma[1])};
        }

    Field field(data, grid_spacing, origin);

    auto const p = Utils::Vector3d{-.4, 3.14, 0.1};

    auto const field_value = field.jacobian(p);

    auto const interpolated_value = bspline_3d_gradient_accumulate<2>(
        p, [&data](const std::array<int, 3> &ind) { return data(ind); },
        grid_spacing, origin, Field::jacobian_type{});

    BOOST_CHECK_SMALL((interpolated_value[0] - field_value[0]).norm(),
                      std::numeric_limits<double>::epsilon());
    BOOST_CHECK_SMALL((interpolated_value[1] - field_value[1]).norm(),
                      std::numeric_limits<double>::epsilon());
  }
}

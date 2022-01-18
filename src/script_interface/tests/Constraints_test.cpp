/*
 * Copyright (C) 2022 The ESPResSo project
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

#define BOOST_TEST_MODULE Constraints test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/constraints/ExternalField.hpp"
#include "script_interface/constraints/ExternalPotential.hpp"
#include "script_interface/constraints/HomogeneousMagneticField.hpp"
#include "script_interface/constraints/ShapeBasedConstraint.hpp"
#include "script_interface/constraints/couplings.hpp"
#include "script_interface/constraints/fields.hpp"

#include "script_interface/get_value.hpp"
#include "script_interface/shapes/NoWhere.hpp"
#include "script_interface/shapes/Shape.hpp"
#include "script_interface/shapes/Sphere.hpp"

#include <utils/Vector.hpp>
#include <utils/as_const.hpp>

#include <limits>
#include <memory>

using namespace ScriptInterface;

BOOST_AUTO_TEST_CASE(shape_based_constraint) {
  ScriptInterface::Constraints::ShapeBasedConstraint constraint;
  {
    // check const and non-const access
    BOOST_TEST(constraint.constraint()->fits_in_box({}));
    BOOST_TEST(Utils::as_const(constraint).constraint()->fits_in_box({}));
    BOOST_TEST(constraint.shape_based_constraint()->fits_in_box({}));
  }
  {
    // check shape setter and getter
    auto const shape = std::make_shared<ScriptInterface::Shapes::NoWhere>();
    constraint.set_parameter("shape", shape);
    auto const shape_si =
        get_value<std::shared_ptr<ScriptInterface::Shapes::Shape>>(
            constraint.get_parameter("shape"));
    auto const shape_core = shape_si->shape();
    BOOST_REQUIRE_EQUAL(shape_si, shape);
    // check distance calculation
    auto constexpr inf = std::numeric_limits<double>::infinity();
    auto const inf_vec = Utils::Vector3d::broadcast(inf);
    double dist;
    Utils::Vector3d vec;
    shape_core->calculate_dist({}, dist, vec);
    BOOST_CHECK_EQUAL(dist, inf);
    BOOST_TEST(vec == inf_vec, boost::test_tools::per_element());
  }
  {
    // check shape setter and getter
    auto shape = std::make_shared<ScriptInterface::Shapes::Sphere>();
    shape->set_parameter("radius", 2.);
    constraint.set_parameter("shape", shape);
    auto const shape_si =
        get_value<std::shared_ptr<ScriptInterface::Shapes::Shape>>(
            constraint.get_parameter("shape"));
    auto const shape_core = shape_si->shape();
    BOOST_REQUIRE_EQUAL(shape_si, shape);
    // check distance calculation
    auto const vec_ref = Utils::Vector3d{0.5, 0., 0.};
    auto const pos = Utils::Vector3d{2.5, 0., 0.};
    auto const variant =
        get_value<std::vector<Variant>>(shape_si->do_call_method(
            "calc_distance", VariantMap{{"position", pos}}));
    auto const dist1 = get_value<double>(variant[0]);
    auto const vec1 = get_value<Utils::Vector3d>(variant[1]);
    double dist2;
    Utils::Vector3d vec2;
    shape_core->calculate_dist(pos, dist2, vec2);
    BOOST_CHECK_EQUAL(dist1, 0.5);
    BOOST_CHECK_EQUAL(dist2, 0.5);
    BOOST_TEST(vec1 == vec_ref, boost::test_tools::per_element());
    BOOST_TEST(vec2 == vec_ref, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_CASE(field_constraints) {
  using namespace FieldCoupling::Coupling;
  using namespace FieldCoupling::Fields;
  using namespace ScriptInterface::Constraints;
  using Gravity = ExternalField<Mass, Constant<double, 3>>;
  {
    // check getters and setters
    HomogeneousMagneticField field;
    auto const h_ref = Utils::Vector3d{1., 2., 3.};
    field.set_parameter("H", Variant{h_ref});
    auto const h = get_value<Utils::Vector3d>(field.get_parameter("H"));
    BOOST_TEST(h == h_ref, boost::test_tools::per_element());
    // check const and non-const access
    BOOST_TEST(field.constraint()->fits_in_box({}));
    BOOST_TEST(Utils::as_const(field).constraint()->fits_in_box({}));
  }
  {
    // check constructor and getters
    Gravity field;
    auto const gravity_constant = Utils::Vector3d{1., 2., 3.};
    field.do_construct({{"value", Variant{gravity_constant}}});
    auto const g = get_value<Utils::Vector3d>(field.get_parameter("value"));
    BOOST_TEST(g == gravity_constant, boost::test_tools::per_element());
    // check const and non-const access
    BOOST_TEST(field.constraint()->fits_in_box({}));
    BOOST_TEST(Utils::as_const(field).constraint()->fits_in_box({}));
  }
}

BOOST_AUTO_TEST_CASE(potential_constraints) {
  using namespace FieldCoupling::Coupling;
  using namespace FieldCoupling::Fields;
  using namespace ScriptInterface::Constraints;
  using ElectricPotential = ExternalPotential<Charge, Interpolated<double, 1>>;
  {
    // check constructor and getters
    ElectricPotential potential;
    auto const grid_spacing_ref = Utils::Vector3d{{1., 1., 1.}};
    auto const field_shape_ref = Utils::Vector3i{{1, 2, 3}};
    auto const field_codim_ref = 1;
    auto const field_data_ref = std::vector<double>{1., 2., 3., 4., 5., 6.};
    potential.do_construct(
        {{"_field_shape", Variant{field_shape_ref}},
         {"_field_codim", Variant{field_codim_ref}},
         {"_field_data", Variant{std::vector<Variant>(field_data_ref.begin(),
                                                      field_data_ref.end())}},
         {"grid_spacing", Variant{grid_spacing_ref}}});
    auto const grid_spacing =
        get_value<Utils::Vector3d>(potential.get_parameter("grid_spacing"));
    auto const field_shape =
        get_value<Utils::Vector3i>(potential.get_parameter("_field_shape"));
    auto const field_codim =
        get_value<int>(potential.get_parameter("_field_codim"));
    auto const field_data =
        get_value<std::vector<double>>(potential.get_parameter("_field_data"));
    BOOST_TEST(grid_spacing == grid_spacing_ref,
               boost::test_tools::per_element());
    BOOST_TEST(field_shape == field_shape_ref,
               boost::test_tools::per_element());
    BOOST_TEST(field_data == field_data_ref, boost::test_tools::per_element());
    BOOST_CHECK_EQUAL(field_codim, field_codim_ref);
    // check const and non-const access
    BOOST_TEST(potential.constraint()->fits_in_box({}));
    BOOST_TEST(Utils::as_const(potential).constraint()->fits_in_box({}));
  }
  {
    // check exception mechanism
    ElectricPotential potential;
    auto const grid_spacing = Utils::Vector3d{{1., 1., 1.}};
    auto const field_data = std::vector<Variant>{1., 2., 3., 4., 5., 6.};
    auto const wrong_params1 =
        VariantMap{{"_field_shape", Variant{Utils::Vector3i{{0, 2, 3}}}},
                   {"_field_codim", Variant{1}},
                   {"_field_data", field_data},
                   {"grid_spacing", Variant{grid_spacing}}};
    auto const wrong_params2 =
        VariantMap{{"_field_shape", Variant{Utils::Vector3i{{1, 2, 3}}}},
                   {"_field_codim", Variant{5}},
                   {"_field_data", field_data},
                   {"grid_spacing", Variant{grid_spacing}}};
    BOOST_CHECK_THROW(potential.do_construct(wrong_params1),
                      std::runtime_error);
    BOOST_CHECK_THROW(potential.do_construct(wrong_params2),
                      std::runtime_error);
  }
}

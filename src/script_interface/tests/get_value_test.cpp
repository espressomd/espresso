/*
 * Copyright (C) 2017-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE ScriptInterface::get_value test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/ObjectHandle.hpp"
#include "script_interface/get_value.hpp"

#include <cassert>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

BOOST_AUTO_TEST_CASE(default_case) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  {
    auto const s = std::string{"Abc"};
    auto const v = Variant(s);

    BOOST_CHECK_EQUAL(get_value<std::string>(v), s);
  }
  {
    auto const vec = Utils::Vector<double, 3>{1., 2., 3.};
    auto const var = Variant{vec};

    BOOST_CHECK_EQUAL((get_value<Utils::Vector<double, 3>>(var)), vec);
    BOOST_CHECK_EQUAL((get_value<Utils::Vector3<double>>(var)), vec);
    BOOST_CHECK_EQUAL((get_value<Utils::Vector3d>(var)), vec);
  }
}

BOOST_AUTO_TEST_CASE(conversions) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;
  using ScriptInterface::detail::allow_conversion;

  static_assert(allow_conversion<int, int>::value);
  static_assert(allow_conversion<double, int>::value);
  static_assert(not allow_conversion<int, double>::value);

  BOOST_CHECK_EQUAL(get_value<double>(3.1415), 3.1415);
  BOOST_CHECK_EQUAL(get_value<double>(3), 3.);
  BOOST_CHECK_EQUAL(get_value<double>(Variant{3}), 3.);
}

BOOST_AUTO_TEST_CASE(static_vector) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  /* From same type */
  {
    Variant v = std::vector<Variant>({1., 2., 3.});
    auto const expected = Utils::Vector3d{1., 2., 3.};
    BOOST_CHECK(get_value<Utils::Vector3d>(v) == expected);
  }
  {
    Variant v = std::vector<double>({1., 2., 3.});
    auto const expected = Utils::Vector3d{1., 2., 3.};
    BOOST_CHECK(get_value<Utils::Vector3d>(v) == expected);
  }

  /* Conversion applied */
  {
    Variant v = std::vector<Variant>({1, 2, 3});
    auto const expected = Utils::Vector3d{1, 2, 3};
    BOOST_CHECK(get_value<Utils::Vector3d>(v) == expected);
  }
}

BOOST_AUTO_TEST_CASE(heap_vector) {
  using ScriptInterface::get_value;
  using ScriptInterface::make_vector_of_variants;
  using ScriptInterface::Variant;

  /* From same type */
  {
    Variant v = std::vector<Variant>({1., 2., 3.});
    auto const expected = std::vector<double>{1., 2., 3.};
    BOOST_CHECK(get_value<std::vector<double>>(v) == expected);
  }
  {
    Variant v = std::vector<double>({1., 2., 3.});
    auto const expected = std::vector<double>{1., 2., 3.};
    BOOST_CHECK(get_value<std::vector<double>>(v) == expected);
  }
  {
    auto const vec_dbl = std::vector<double>({1., 2., 3.});
    auto const vec_var = make_vector_of_variants(vec_dbl);
    auto const expected = std::vector<double>{1., 2., 3.};
    BOOST_CHECK(get_value<std::vector<double>>(Variant{vec_var}) == expected);
  }

  /* Conversion applied */
  {
    Variant v = std::vector<Variant>({1, 2, 3});
    auto const expected = std::vector<double>{1, 2, 3};
    BOOST_CHECK(get_value<std::vector<double>>(v) == expected);
  }
}

BOOST_AUTO_TEST_CASE(get_value_from_map) {
  using ScriptInterface::get_value;
  using ScriptInterface::get_value_or;
  using ScriptInterface::Variant;
  using ScriptInterface::VariantMap;

  VariantMap map{{"a", 13}, {"e", 3.1}, {"f", "s"}};

  BOOST_CHECK(3.1 == get_value<double>(map, "e"));
  BOOST_CHECK(13 == get_value_or(map, "a", -1));
  BOOST_CHECK(-1 == get_value_or(map, "nope", -1));
  BOOST_CHECK_THROW((get_value<int>(map, "unknown")), std::exception);
}

BOOST_AUTO_TEST_CASE(unordered_map) {
  using ScriptInterface::get_value;
  using ScriptInterface::make_unordered_map_of_variants;
  using ScriptInterface::Variant;

  auto const var = Variant{std::unordered_map<int, Variant>{{1, 1}, {2, 2.5}}};
  {
    auto const map_var = get_value<std::unordered_map<int, Variant>>(var);
    BOOST_CHECK_EQUAL(get_value<int>(map_var.at(1)), 1);
    BOOST_CHECK_EQUAL(get_value<double>(map_var.at(2)), 2.5);
  }
  {
    auto const map_dbl = get_value<std::unordered_map<int, double>>(var);
    BOOST_CHECK_EQUAL(map_dbl.at(1), 1.0);
    BOOST_CHECK_EQUAL(map_dbl.at(2), 2.5);
  }
  {
    auto const map_dbl_input = get_value<std::unordered_map<int, double>>(var);
    auto const map_var = make_unordered_map_of_variants(map_dbl_input);
    auto const map_dbl =
        get_value<std::unordered_map<int, double>>(Variant{map_var});
    BOOST_CHECK_EQUAL(map_dbl.at(1), 1.0);
    BOOST_CHECK_EQUAL(map_dbl.at(2), 2.5);
  }
}

auto exception_message_predicate(std::string const &pattern) {
  return [=](std::exception const &ex) {
    boost::test_tools::predicate_result result = true;
    std::string const what = ex.what();
    std::smatch match;
    if (!std::regex_search(what, match, std::regex(pattern))) {
      result = false;
      result.message() << "Error message \"" << what << "\" "
                       << "doesn't match pattern \"" << pattern << "\"";
    }
    return result;
  };
}

BOOST_AUTO_TEST_CASE(check_exceptions) {
  using ScriptInterface::get_value;
  using ScriptInterface::Variant;

  assert(!!exception_message_predicate("A")(std::runtime_error("A")));
  assert(!exception_message_predicate("A")(std::runtime_error("B")));

  using so_ptr_t = std::shared_ptr<ScriptInterface::ObjectHandle>;

  auto const so_obj = so_ptr_t();
  auto const msg_prefix = std::string("Provided argument of type ");
  auto const variant_sip_name =
      "ScriptInterface::Variant\\{" + Utils::demangle<so_ptr_t>() + "\\}";

  {
    // basic types
    auto const obj_variant = Variant{so_obj};
    auto const obj_variant_pattern = Utils::demangle<so_ptr_t>();
    auto const what = msg_prefix + "'" + obj_variant_pattern + "'";
    auto const predicate_nullptr =
        exception_message_predicate(what + " is a null pointer");
    auto const predicate_conversion =
        exception_message_predicate(what + " is not convertible to 'int'");
    BOOST_CHECK_EXCEPTION(get_value<so_ptr_t>(obj_variant), std::exception,
                          predicate_nullptr);
    BOOST_CHECK_EXCEPTION(get_value<int>(obj_variant), std::exception,
                          predicate_conversion);
  }
  {
    // vectors
    auto const int_variant = Variant{1.5};
    auto const vec_variant = Variant{std::vector<Variant>{{so_obj}}};
    auto const vec_variant_pattern = "std::vector<" + variant_sip_name + ">";
    auto const what = msg_prefix + "'" + vec_variant_pattern + "'";
    auto const predicate_nullptr = exception_message_predicate(
        what + " contains a value that is a null pointer");
    auto const predicate_conversion_containee = exception_message_predicate(
        what + " is not convertible to 'std::vector<int>' because"
               " it contains a value that is not convertible to 'int'");
    auto const predicate_conversion = exception_message_predicate(
        msg_prefix + "'double' is not convertible to 'std::vector<int>'");
    BOOST_CHECK_EXCEPTION(get_value<std::vector<so_ptr_t>>(vec_variant),
                          std::exception, predicate_nullptr);
    BOOST_CHECK_EXCEPTION(get_value<std::vector<int>>(vec_variant),
                          std::exception, predicate_conversion_containee);
    BOOST_CHECK_EXCEPTION(get_value<std::vector<int>>(int_variant),
                          std::exception, predicate_conversion);
  }
  {
    // unordered maps with integral key
    auto const map_variant =
        Variant{std::unordered_map<int, Variant>{{1, so_obj}}};
    auto const map_variant_pattern =
        "std::unordered_map<int, " + variant_sip_name + ">";
    auto const what = msg_prefix + "'" + map_variant_pattern + "'";
    auto const predicate_nullptr = exception_message_predicate(
        what + " contains a value that is a null pointer");
    auto const predicate_conversion = exception_message_predicate(
        what +
        " is not convertible to 'std::unordered_map<int, double>' because"
        " it contains a value that is not convertible to 'int' or 'double'");
    BOOST_CHECK_EXCEPTION(
        (get_value<std::unordered_map<int, so_ptr_t>>(map_variant)),
        std::exception, predicate_nullptr);
    BOOST_CHECK_EXCEPTION(
        (get_value<std::unordered_map<int, double>>(map_variant)),
        std::exception, predicate_conversion);
  }
  {
    // unordered maps with string key
    auto const map_variant =
        Variant{std::unordered_map<std::string, Variant>{{"key", so_obj}}};
    auto const map_variant_pattern =
        "std::unordered_map<std::string, " + variant_sip_name + ">";
    auto const what = msg_prefix + "'" + map_variant_pattern + "'";
    auto const predicate_nullptr = exception_message_predicate(
        what + " contains a value that is a null pointer");
    auto const predicate_conversion = exception_message_predicate(
        what +
        " is not convertible to 'std::unordered_map<std::string, int>' because"
        " it contains a value that is not convertible to 'std::string' or "
        "'int'");
    BOOST_CHECK_EXCEPTION(
        (get_value<std::unordered_map<std::string, so_ptr_t>>(map_variant)),
        std::exception, predicate_nullptr);
    BOOST_CHECK_EXCEPTION(
        (get_value<std::unordered_map<std::string, int>>(map_variant)),
        std::exception, predicate_conversion);
  }
  {
    using Utils::Vector3d;
    std::unordered_map<int, Variant> const mixed{{1, 1}, {2, std::string("2")}};
    std::vector<Variant> const v_var = {1., 2.};
    std::vector<double> const v_dbl = {1., 2.};
    BOOST_CHECK_THROW((get_value<Vector3d>(Variant{v_var})), std::exception);
    BOOST_CHECK_THROW((get_value<Vector3d>(Variant{v_dbl})), std::exception);
    Utils::Vector4d const quat{};
    BOOST_CHECK_THROW((get_value<Vector3d>(Variant{quat})), std::exception);
    BOOST_CHECK_THROW((get_value<Vector3d>(Variant{1.})), std::exception);
    BOOST_CHECK_THROW((get_value<std::vector<int>>(Variant{})), std::exception);
    BOOST_CHECK_THROW((get_value<int>(Variant{v_dbl})), std::exception);
    BOOST_CHECK_THROW((get_value<std::unordered_map<int, Variant>>(Variant{})),
                      std::exception);
    BOOST_CHECK_THROW((get_value<std::unordered_map<int, int>>(Variant{mixed})),
                      std::exception);
    BOOST_CHECK_THROW(
        (get_value<std::unordered_map<std::string, int>>(Variant{mixed})),
        std::exception);
  }
}

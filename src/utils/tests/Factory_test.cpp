/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE Factory test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Factory.hpp"

#include <stdexcept>
#include <string>

struct TestClass {
  virtual void method() = 0;
  virtual ~TestClass() = default;
};

struct DerivedTestClass : public TestClass {
  void method() override {}
};

BOOST_AUTO_TEST_CASE(make) {
  Utils::Factory<TestClass> factory;
  auto const derived_class_name = std::string{"derived_test_class"};

  // Register construction function
  factory.register_new<DerivedTestClass>(derived_class_name);

  // Check registration of construction function
  BOOST_REQUIRE(factory.has_builder(derived_class_name));

  // Make a derived object
  auto o = factory.make(derived_class_name);
  BOOST_REQUIRE(o);
  o->method();

  // Check for correct type name
  BOOST_CHECK_EQUAL(factory.type_name(*o.get()), derived_class_name);

  // Check for correct (derived) type
  BOOST_CHECK(dynamic_cast<DerivedTestClass *>(o.get()) != nullptr);

  // Make an unknown object
  BOOST_CHECK(not factory.has_builder("unknown"));
  BOOST_CHECK_THROW(factory.make("unknown"), std::domain_error);
}

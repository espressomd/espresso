/*
 * Copyright (C) 2015-2019 The ESPResSo project
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

/* Unit tests for the Utils::Factory class.
 * The factory is tested by registering different types of classes
 * with it (first test), and then checking if instances of those classes can be
 * made via the Factory (second test).
 */

#define BOOST_TEST_MODULE Factory test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Factory.hpp"

namespace Testing {
struct TestClass {
  virtual void method() = 0;
  virtual ~TestClass() = default;
};

struct DerivedTestClass : public TestClass {
  void method() override {}
};

struct OtherDerivedTestClass : public TestClass {
  void method() override {}
};
} /* namespace Testing */

/* Check registration of construction functions */
BOOST_AUTO_TEST_CASE(regiser_class) {
  using namespace Testing;
  Utils::Factory<TestClass> factory;

  /* Overload with explicit builder */
  factory.register_new("derived_test_class",
                       []() { return new DerivedTestClass(); });
  /* Overload with default builder */
  factory.register_new<OtherDerivedTestClass>(
      "other_derived_class");

  /* Both builders should be present. */
  BOOST_CHECK(factory.has_builder("derived_test_class"));
  BOOST_CHECK(factory.has_builder("other_derived_class"));
}

/* Check object construction. */
BOOST_AUTO_TEST_CASE(make) {
  using namespace Testing;

  using namespace Testing;
  Utils::Factory<TestClass> factory;

  /* Overload with explicit builder */
  factory.register_new("derived_test_class",
                       []() { return new DerivedTestClass(); });

  /* Make a derived object */
  auto o = factory.make("derived_test_class");
  BOOST_CHECK(o);

  /* Check for correct (derived) type */
  BOOST_CHECK(dynamic_cast<DerivedTestClass *>(o.get()) != nullptr);
}

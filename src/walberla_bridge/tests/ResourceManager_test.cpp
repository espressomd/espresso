/*
 * Copyright (C) 2023 The ESPResSo project
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
#define BOOST_TEST_MODULE resources manager
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <walberla_bridge/utils/ResourceManager.hpp>

#include <memory>
#include <string>
#include <vector>

namespace Testing {

static std::vector<std::string> logger;

template <char Name> class LogWriter {

public:
  LogWriter() { logger.emplace_back(std::string(1, Name) + "()"); }
  ~LogWriter() { logger.emplace_back("~" + std::string(1, Name) + "()"); }
};

} // namespace Testing

BOOST_AUTO_TEST_CASE(destruction_order) {
  // instantiate three resources in a specific order
  auto obj_a = std::make_shared<Testing::LogWriter<'A'>>();
  auto obj_b = std::make_shared<Testing::LogWriter<'B'>>();
  auto obj_c = std::make_shared<Testing::LogWriter<'C'>>();
  auto obj_d = std::make_shared<Testing::LogWriter<'D'>>();
  BOOST_REQUIRE_EQUAL(Testing::logger.size(), 4ul);
  BOOST_CHECK_EQUAL(Testing::logger[0], "A()");
  BOOST_CHECK_EQUAL(Testing::logger[1], "B()");
  BOOST_CHECK_EQUAL(Testing::logger[2], "C()");
  BOOST_CHECK_EQUAL(Testing::logger[3], "D()");

  // lock resources in some order (but *not* the reverse order of construction!)
  auto manager = std::make_unique<ResourceManager>();
  manager->acquire_lock(obj_c);
  manager->acquire_lock(obj_a);
  manager->acquire_lock(obj_d);
  manager->acquire_lock(obj_b);
  BOOST_REQUIRE_EQUAL(Testing::logger.size(), 4ul);

  // resetting the local shared pointers should not make the resources expire
  obj_a.reset();
  obj_b.reset();
  obj_c.reset();
  obj_d.reset();
  BOOST_REQUIRE_EQUAL(Testing::logger.size(), 4ul);

  // the manager should free the resources in the reverse order of their locking
  manager.reset();
  BOOST_REQUIRE_EQUAL(Testing::logger.size(), 8ul);
  BOOST_CHECK_EQUAL(Testing::logger[4], "~B()");
  BOOST_CHECK_EQUAL(Testing::logger[5], "~D()");
  BOOST_CHECK_EQUAL(Testing::logger[6], "~A()");
  BOOST_CHECK_EQUAL(Testing::logger[7], "~C()");
}

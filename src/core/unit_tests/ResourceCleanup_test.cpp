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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ResourceCleanup test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "communication.hpp"
#include "cuda/utils.hpp"
#include "system/GpuParticleData.hpp"
#include "system/ResourceCleanup.hpp"
#include "system/System.hpp"

#include <boost/mpi.hpp>

#include <memory>
#include <vector>

class MyClass {
  std::vector<int> m_data;
  void deallocate() { m_data.clear(); }
  using Cleanup = ResourceCleanup::Attorney<&MyClass::deallocate>;
  friend Cleanup;

public:
  MyClass() { m_data = std::vector<int>(5); }
  ~MyClass() { deallocate(); }
  auto size() const { return m_data.size(); }
  template <class... Args> static auto make_shared(Args... args) {
    auto obj = std::make_shared<MyClass>(args...);
    System::get_system().cleanup_queue.push<Cleanup>(obj);
    return obj;
  }
};

BOOST_AUTO_TEST_CASE(checks) {
  auto system = std::make_shared<::System::System>();
  System::set_system(system);

  BOOST_REQUIRE_EQUAL(system->cleanup_queue.size(), 0);
  BOOST_REQUIRE_EQUAL(system->cleanup_queue.empty(), true);
  system->init();

#ifdef CUDA
  if (system->gpu.has_compatible_device()) {
    // allocate device memory to populate the cleanup queue
    system->gpu.enable_property(GpuParticleData::prop::pos);
    system->gpu.update();
    BOOST_REQUIRE_EQUAL(system->cleanup_queue.size(), 1);
    BOOST_REQUIRE_EQUAL(system->cleanup_queue.empty(), false);
  }
#endif

  auto const obj = MyClass::make_shared();
  BOOST_REQUIRE_EQUAL(system->cleanup_queue.empty(), false);
  BOOST_REQUIRE_EQUAL(obj->size(), 5);
  system.reset();
  System::reset_system();
  BOOST_REQUIRE_EQUAL(obj->size(), 0);
}

int main(int argc, char **argv) {
  auto mpi_env = mpi_init(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

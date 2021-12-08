/*
 * Copyright (C) 2021 The ESPResSo project
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
#define BOOST_TEST_MODULE EspressoSystemInterface test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "particle_data.hpp"
#include "virtual_sites.hpp"
#include "virtual_sites/VirtualSitesOff.hpp"

#include <boost/mpi.hpp>

#include <cassert>
#include <memory>

namespace espresso {
auto &system = EspressoSystemInterface::Instance();
}

inline void check_uninitialized_device_pointers() {
  BOOST_CHECK_EQUAL(espresso::system.eGpu(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.rGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.dipGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.fGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.torqueGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.qGpuBegin(), nullptr);
}

#ifdef CUDA

#include "cuda_init.hpp"
#include "cuda_utils.hpp"

/* Decorator to run a unit test depending on GPU availability. */
boost::test_tools::assertion_result has_gpu(boost::unit_test::test_unit_id) {
  bool has_compatible_gpu = false;
  try {
    cuda_check_device();
    has_compatible_gpu = true;
  } catch (cuda_runtime_error const &err) {
  }
  return has_compatible_gpu;
}

BOOST_AUTO_TEST_CASE(check_with_gpu, *boost::unit_test::precondition(has_gpu)) {
  check_uninitialized_device_pointers();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);

  auto const pid = 1;
  place_particle(pid, {0., 0., 0.});
  set_particle_type(pid, 0);

  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  espresso::system.update();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  espresso::system.requestParticleStructGpu();
  espresso::system.update();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 1);

  // check position split
  BOOST_TEST(espresso::system.hasRGpu());
  espresso::system.requestRGpu();
  espresso::system.update();
  BOOST_TEST(espresso::system.rGpuBegin() != nullptr);

  // check force split
  BOOST_TEST(espresso::system.hasFGpu());
  espresso::system.requestFGpu();
  espresso::system.update();
  BOOST_TEST(espresso::system.fGpuBegin() != nullptr);

  // check torque split
#ifdef ROTATION
  BOOST_CHECK(espresso::system.hasTorqueGpu());
  espresso::system.requestTorqueGpu();
  espresso::system.update();
  BOOST_TEST(espresso::system.torqueGpuBegin() != nullptr);
#else
  BOOST_CHECK(!espresso::system.hasTorqueGpu());
  BOOST_CHECK_THROW(espresso::system.requestTorqueGpu(), std::runtime_error);
#endif

  // check charge split
#ifdef ELECTROSTATICS
  BOOST_CHECK(espresso::system.hasQGpu());
  espresso::system.requestQGpu();
  espresso::system.update();
  BOOST_TEST(espresso::system.qGpuBegin() != nullptr);
#else
  BOOST_CHECK(!espresso::system.hasQGpu());
  BOOST_CHECK_THROW(espresso::system.requestQGpu(), std::runtime_error);
#endif

  // check dipole split
#ifdef DIPOLES
  BOOST_CHECK(espresso::system.hasDipGpu());
  espresso::system.requestDipGpu();
  espresso::system.update();
  BOOST_TEST(espresso::system.dipGpuBegin() != nullptr);
#else
  BOOST_CHECK(!espresso::system.hasDipGpu());
  BOOST_CHECK_THROW(espresso::system.requestDipGpu(), std::runtime_error);
#endif

  // clear device memory
  remove_particle(pid);
  espresso::system.update();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
}

#else // CUDA

BOOST_AUTO_TEST_CASE(check_without_cuda) {
  check_uninitialized_device_pointers();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  BOOST_TEST(!espresso::system.hasRGpu());
  BOOST_TEST(!espresso::system.hasDipGpu());
  BOOST_TEST(!espresso::system.hasFGpu());
  BOOST_TEST(!espresso::system.hasTorqueGpu());
  BOOST_TEST(!espresso::system.hasQGpu());
  BOOST_CHECK_THROW(espresso::system.requestRGpu(), std::runtime_error);
  BOOST_CHECK_THROW(espresso::system.requestDipGpu(), std::runtime_error);
  BOOST_CHECK_THROW(espresso::system.requestFGpu(), std::runtime_error);
  BOOST_CHECK_THROW(espresso::system.requestTorqueGpu(), std::runtime_error);
  BOOST_CHECK_THROW(espresso::system.requestQGpu(), std::runtime_error);
}

#endif // CUDA

int main(int argc, char **argv) {
  auto mpi_env = mpi_init(argc, argv);

  // this unit test only runs on 1 core
  assert(boost::mpi::communicator().size() == 1);

  // initialize the MpiCallbacks framework
  Communication::init(mpi_env);
#ifdef VIRTUAL_SITES
  set_virtual_sites(std::make_shared<VirtualSitesOff>());
#endif
  mpi_loop();

  espresso::system.init();

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

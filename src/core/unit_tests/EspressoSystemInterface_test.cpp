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
  BOOST_CHECK_EQUAL(espresso::system.rGpuEnd(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.dipGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.dipGpuEnd(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.fGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.fGpuEnd(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.torqueGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.torqueGpuEnd(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.qGpuBegin(), nullptr);
  BOOST_CHECK_EQUAL(espresso::system.qGpuEnd(), nullptr);
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

  // features compiled in
  auto has_feature_rotation = false;
  auto has_feature_electrostatics = false;
  auto has_feature_dipoles = false;
#ifdef ROTATION
  has_feature_rotation = true;
#endif
#ifdef ELECTROSTATICS
  has_feature_electrostatics = true;
#endif
#ifdef DIPOLES
  has_feature_dipoles = true;
#endif

  auto const pid = 1;
  place_particle(pid, {0., 0., 0.});
  set_particle_type(pid, 0);

  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  espresso::system.update();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  BOOST_TEST(espresso::system.requestParticleStructGpu());
  espresso::system.update();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 1);

  // check position split
  BOOST_TEST(espresso::system.hasRGpu());
  BOOST_TEST(espresso::system.requestRGpu());
  espresso::system.update();
  BOOST_TEST(espresso::system.rGpuBegin() != nullptr);
  BOOST_TEST(espresso::system.rGpuEnd() != nullptr);
  BOOST_CHECK_EQUAL(espresso::system.rGpuEnd() - espresso::system.rGpuBegin(),
                    3);

  // check force split
  BOOST_TEST(espresso::system.hasFGpu());
  BOOST_TEST(espresso::system.requestFGpu());
  espresso::system.update();
  BOOST_TEST(espresso::system.fGpuBegin() != nullptr);
  BOOST_TEST(espresso::system.fGpuEnd() != nullptr);
  BOOST_CHECK_EQUAL(espresso::system.fGpuEnd() - espresso::system.fGpuBegin(),
                    3);

  // check torque split
  BOOST_CHECK_EQUAL(espresso::system.hasTorqueGpu(), has_feature_rotation);
  BOOST_CHECK_EQUAL(espresso::system.requestTorqueGpu(), has_feature_rotation);
#ifdef ROTATION
  espresso::system.update();
  BOOST_TEST(espresso::system.torqueGpuBegin() != nullptr);
  BOOST_TEST(espresso::system.torqueGpuEnd() != nullptr);
  BOOST_CHECK_EQUAL(
      espresso::system.torqueGpuEnd() - espresso::system.torqueGpuBegin(), 3);
#endif

  // check charge split
  BOOST_CHECK_EQUAL(espresso::system.hasQGpu(), has_feature_electrostatics);
  BOOST_CHECK_EQUAL(espresso::system.requestQGpu(), has_feature_electrostatics);
#ifdef ELECTROSTATICS
  espresso::system.update();
  BOOST_TEST(espresso::system.qGpuBegin() != nullptr);
  BOOST_TEST(espresso::system.qGpuEnd() != nullptr);
  BOOST_CHECK_EQUAL(espresso::system.qGpuEnd() - espresso::system.qGpuBegin(),
                    3);
#endif

  // check dipole split
  BOOST_CHECK_EQUAL(espresso::system.hasDipGpu(), has_feature_dipoles);
  BOOST_CHECK_EQUAL(espresso::system.requestDipGpu(), has_feature_dipoles);
#ifdef DIPOLES
  espresso::system.update();
  BOOST_TEST(espresso::system.dipGpuBegin() != nullptr);
  BOOST_TEST(espresso::system.dipGpuEnd() != nullptr);
  BOOST_CHECK_EQUAL(
      espresso::system.dipGpuEnd() - espresso::system.dipGpuBegin(), 3);
#endif

  // clear device memory
  remove_particle(pid);
  espresso::system.update();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  BOOST_CHECK_EQUAL(espresso::system.rGpuBegin(), espresso::system.rGpuEnd());
  BOOST_CHECK_EQUAL(espresso::system.fGpuBegin(), espresso::system.fGpuEnd());
#ifdef ELECTROSTATICS
  BOOST_CHECK_EQUAL(espresso::system.qGpuBegin(), espresso::system.qGpuEnd());
#endif
#ifdef DIPOLES
  BOOST_CHECK_EQUAL(espresso::system.dipGpuBegin(),
                    espresso::system.dipGpuEnd());
#endif
#ifdef ROTATION
  BOOST_CHECK_EQUAL(espresso::system.torqueGpuBegin(),
                    espresso::system.torqueGpuEnd());
#endif
}

#else // CUDA

BOOST_AUTO_TEST_CASE(check_without_cuda) {
  check_uninitialized_device_pointers();
  BOOST_CHECK_EQUAL(espresso::system.npart_gpu(), 0);
  BOOST_TEST(!espresso::system.hasRGpu());
  BOOST_TEST(!espresso::system.requestRGpu());
  BOOST_TEST(!espresso::system.hasDipGpu());
  BOOST_TEST(!espresso::system.requestDipGpu());
  BOOST_TEST(!espresso::system.hasFGpu());
  BOOST_TEST(!espresso::system.requestFGpu());
  BOOST_TEST(!espresso::system.hasTorqueGpu());
  BOOST_TEST(!espresso::system.requestTorqueGpu());
  BOOST_TEST(!espresso::system.hasQGpu());
  BOOST_TEST(!espresso::system.requestQGpu());
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

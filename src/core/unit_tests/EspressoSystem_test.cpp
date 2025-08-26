/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE System test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "EspressoCoreGlobalConfig.hpp"
#include "ParticleFactory.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cells.hpp"
#include "particle_node.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"
#include "system/System.impl.hpp"

#include "cuda/init.hpp"
#include "cuda/utils.hpp"

#include <boost/mpi.hpp>

#include <cassert>
#include <memory>

struct GlobalConfig : public EspressoCoreGlobalConfig {
  ~GlobalConfig() { ::System::reset_system(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

/* Decorator to run a unit test depending on GPU availability. */
boost::test_tools::assertion_result has_gpu(boost::unit_test::test_unit_id) {
  bool has_compatible_gpu = false;
#ifdef ESPRESSO_CUDA
  invoke_skip_cuda_exceptions([&]() {
    cuda_check_device();
    has_compatible_gpu = true;
  });
#endif // ESPRESSO_CUDA
  return has_compatible_gpu;
}

BOOST_FIXTURE_TEST_CASE(check_with_gpu, ParticleFactory,
                        *boost::unit_test::precondition(has_gpu)) {
#ifdef ESPRESSO_CUDA
  auto const rank = boost::mpi::communicator().rank();

  auto system = ::System::System::create();
  System::set_system(system);
  system->set_cell_structure_topology(CellStructureType::REGULAR);
  auto &gpu = system->gpu;

  // check uninitialized device pointers
  BOOST_CHECK_EQUAL(gpu.get_energy_device(), nullptr);
  BOOST_CHECK_EQUAL(gpu.get_particle_positions_device(), nullptr);
#ifdef ESPRESSO_DIPOLES
  BOOST_CHECK_EQUAL(gpu.get_particle_dipoles_device(), nullptr);
#endif
  BOOST_CHECK_EQUAL(gpu.get_particle_forces_device(), nullptr);
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(gpu.get_particle_torques_device(), nullptr);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  BOOST_CHECK_EQUAL(gpu.get_particle_charges_device(), nullptr);
#endif
  BOOST_CHECK_EQUAL(gpu.n_particles(), 0);

  auto const p_id = 1;
  create_particle({0., 0., 0.}, p_id, 0);

  BOOST_CHECK_EQUAL(gpu.n_particles(), 0);
  gpu.update();
  BOOST_CHECK_EQUAL(gpu.n_particles(), 0);

  // check position split
  gpu.enable_property(GpuParticleData::prop::pos);
  gpu.update();
  BOOST_CHECK_EQUAL(gpu.n_particles(), (rank == 0) ? 1 : 0);
  if (rank == 0) {
    BOOST_TEST(gpu.get_particle_positions_device() != nullptr);
  } else {
    BOOST_TEST(gpu.get_particle_positions_device() == nullptr);
  }

  // check force split
  gpu.enable_property(GpuParticleData::prop::force);
  gpu.update();
  if (rank == 0) {
    BOOST_TEST(gpu.get_particle_forces_device() != nullptr);
  } else {
    BOOST_TEST(gpu.get_particle_forces_device() == nullptr);
  }

  // check torque split
#ifdef ESPRESSO_ROTATION
  gpu.enable_property(GpuParticleData::prop::torque);
  gpu.update();
  if (rank == 0) {
    BOOST_TEST(gpu.get_particle_torques_device() != nullptr);
  } else {
    BOOST_TEST(gpu.get_particle_torques_device() == nullptr);
  }
#endif

  // check charge split
#ifdef ESPRESSO_ELECTROSTATICS
  gpu.enable_property(GpuParticleData::prop::q);
  gpu.update();
  if (rank == 0) {
    BOOST_TEST(gpu.get_particle_charges_device() != nullptr);
  } else {
    BOOST_TEST(gpu.get_particle_charges_device() == nullptr);
  }
#endif

  // check dipole split
#ifdef ESPRESSO_DIPOLES
  gpu.enable_property(GpuParticleData::prop::dip);
  gpu.update();
  if (rank == 0) {
    BOOST_TEST(gpu.get_particle_dipoles_device() != nullptr);
  } else {
    BOOST_TEST(gpu.get_particle_dipoles_device() == nullptr);
  }
#endif

  // clear device memory
  remove_particle(p_id);
  gpu.update();
  BOOST_CHECK_EQUAL(gpu.n_particles(), 0);

  clear_particles();
  System::reset_system();
#endif // ESPRESSO_CUDA
}

BOOST_AUTO_TEST_SUITE_END()

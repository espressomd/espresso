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

/* Unit tests for the particle tracking mechanism. */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Particle tracking test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ParticleFactory.hpp"

#include "communication.hpp"
#include "particle_data.hpp"

#include <boost/mpi.hpp>

#include <limits>
#include <memory>
#include <stdexcept>

// Check the mechanism that tracks particles of a certain type and the
// function that selects a random particle in the pool of tracked particles.
BOOST_FIXTURE_TEST_CASE(particle_type_map_test, ParticleFactory) {
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // particle properties
  int const type = 10;
  int const pid = 1;

  // exception for untracked particle ids
  BOOST_CHECK_THROW(number_of_particles_with_type(type), std::runtime_error);

  // exception for negative particle ids
  BOOST_CHECK_THROW(init_type_map(-10), std::runtime_error);

  // check particle counting
  init_type_map(type);
  BOOST_CHECK_EQUAL(number_of_particles_with_type(type), 0);
  create_particle(pid, type);
  BOOST_CHECK_EQUAL(number_of_particles_with_type(type), 1);

  // exception for random index that exceeds the number of particles
  BOOST_CHECK_THROW(get_random_p_id(type, 10), std::runtime_error);

  // check particle selection
  BOOST_CHECK_EQUAL(get_random_p_id(type, 0), pid);
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

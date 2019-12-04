/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef RANDOM_H
#define RANDOM_H

/** \file

    A random generator
*/

#include "errorhandling.hpp"

#include <Random123/philox.h>
#include <utils/Vector.hpp>
#include <utils/u32_to_u64.hpp>
#include <utils/uniform.hpp>

#include <random>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * @brief Salt for the RNGs
 *
 * This is to avoid correlations between the
 * noise on the particle coupling and the fluid
 * thermalization.
 */
enum class RNGSalt : uint64_t {
  FLUID = 0,
  PARTICLES,
  LANGEVIN,
  LANGEVIN_ROT,
  SALT_DPD,
  THERMALIZED_BOND
};

namespace Random {
/**
 * @brief 3d uniform vector noise.
 *
 * This uses the Philox PRNG, the state is controlled
 * by the counter, the salt and two keys.
 * If any of the keys and salt differ, the noise is
 * not correlated between two calls along the same counter
 * sequence.
 *
 */
template <RNGSalt salt>
Utils::Vector3d v_noise(uint64_t counter, int key1, int key2 = 0) {

  using rng_type = r123::Philox4x64;
  using ctr_type = rng_type::ctr_type;
  using key_type = rng_type::key_type;

  const ctr_type c{{counter, static_cast<uint64_t>(salt)}};

  auto const id1 = static_cast<uint32_t>(key1);
  auto const id2 = static_cast<uint32_t>(key2);
  const key_type k{id1, id2};

  auto const noise = rng_type{}(c, k);

  using Utils::uniform;
  return Utils::Vector3d{uniform(noise[0]), uniform(noise[1]),
                         uniform(noise[2])} -
         Utils::Vector3d::broadcast(0.5);
}

extern std::mt19937 generator;
extern std::uniform_real_distribution<double> uniform_real_distribution;
extern bool user_has_seeded;
inline void unseeded_error() {
  runtimeErrorMsg() << "Please seed the random number generator.\nESPResSo "
                       "can choose one for you with set_random_state_PRNG().";
}

/**
 * @brief checks the seeded state and throws error if unseeded
 *
 **/
inline void check_user_has_seeded() {
  static bool unseeded_error_thrown = false;
  if (!user_has_seeded && !unseeded_error_thrown) {
    unseeded_error_thrown = true;
    unseeded_error();
  }
}

/**
 * @brief Set seed of random number generators on each node.
 *
 * @param cnt   Unused.
 * @param seeds A vector of seeds, must be at least n_nodes long.
 **/
void mpi_random_seed(int cnt, std::vector<int> &seeds);

/**
 * @brief Gets a string representation of the state of all
 *        the nodes.
 */
std::string mpi_random_get_stat();

/**
 * @brief Set the seeds on all the node to the state represented
 *        by the string.
 * The string representation must be one that was returned by
 * mpi_random_get_stat.
 */
void mpi_random_set_stat(const std::vector<std::string> &stat);

/**
 * @brief
 * Get the state size of the random number generator
 */
int get_state_size_of_generator();

/**
 * @brief Initialize PRNG with MPI rank as seed.
 */
void init_random();

/**
 * @brief Initialize PRNG with user-provided seed.
 *
 * @param seed seed
 */
void init_random_seed(int seed);

} // namespace Random

/**
 * @brief Draws a random real number from the uniform distribution in the range
 * [0,1)
 */

inline double d_random() {
  using namespace Random;
  check_user_has_seeded();
  return uniform_real_distribution(generator);
}

#endif

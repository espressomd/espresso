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
 *  Random number generation using Philox.
 */

#include "errorhandling.hpp"

#include <Random123/philox.h>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>
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
  BROWNIAN_WALK,
  BROWNIAN_INC,
  BROWNIAN_ROT_INC,
  BROWNIAN_ROT_WALK,
  NPTISO0_HALF_STEP1,
  NPTISO0_HALF_STEP2,
  NPTISOV,
  SALT_DPD,
  THERMALIZED_BOND
};

namespace Random {
/**
 * @brief get 4 random uint 64 from the Philox RNG
 *
 * This uses the Philox PRNG, the state is controlled
 * by the counter, the salt and two keys.
 * If any of the keys and salt differ, the noise is
 * not correlated between two calls along the same counter
 * sequence.
 *
 */
template <RNGSalt salt>
Utils::Vector<uint64_t, 4> philox_4_uint64s(uint64_t counter, int key1,
                                            int key2 = 0) {

  using rng_type = r123::Philox4x64;
  using ctr_type = rng_type::ctr_type;
  using key_type = rng_type::key_type;

  const ctr_type c{{counter, static_cast<uint64_t>(salt)}};

  auto const id1 = static_cast<uint32_t>(key1);
  auto const id2 = static_cast<uint32_t>(key2);
  const key_type k{id1, id2};

  auto const res = rng_type{}(c, k);
  return {res[0], res[1], res[2], res[3]};
}

/**
 * @brief Uniform noise.
 *
 * Mean = 0, variance = 1 / 12.
 * This uses the Philox PRNG, the state is controlled
 * by the counter, the salt and two keys.
 * If any of the keys and salt differ, the noise is
 * not correlated between two calls along the same counter
 * sequence.
 *
 */
template <RNGSalt salt> double noise(uint64_t counter, int key1, int key2 = 0) {

  auto const noise = philox_4_uint64s<salt>(counter, key1, key2);

  using Utils::uniform;
  return uniform(noise[0]) - 0.5;
}

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

  auto const noise = philox_4_uint64s<salt>(counter, key1, key2);

  using Utils::uniform;
  return Utils::Vector3d{uniform(noise[0]), uniform(noise[1]),
                         uniform(noise[2])} -
         Utils::Vector3d::broadcast(0.5);
}

/** @brief Generator for Gaussian random 3d vector.
 *
 * Mean = 0, standard deviation = 1.0
 * Based on the Philox RNG using 4x64 bits.
 * The Box-Muller transform is used to convert from uniform to normal
 * distribution. The transform is only valid, if the uniformly distributed
 * random numbers are not zero (approx one in 2^64). To avoid this case,
 * such numbers are replaced by std::numeric_limits<double>::min()
 * This breaks statistics in rare cases but allows for consistent RNG
 * counters across MPI ranks.
 *
 * @param counter counter for the random number generation
 * @param key1 key for random number generation
 * @param key2 key for random number generation
 * @tparam salt (decorrelates different thermostat types)
 *
 * @return 3D vector of Gaussian random numbers.
 *
 */
template <RNGSalt salt>
inline Utils::Vector3d v_noise_g(uint64_t counter, int key1, int key2 = 0) {

  auto const noise = philox_4_uint64s<salt>(counter, key1, key2);
  using Utils::uniform;
  Utils::Vector4d u{uniform(noise[0]), uniform(noise[1]), uniform(noise[2]),
                    uniform(noise[3])};

  // Replace zeros from uniformly distributed numbers (see doc string)
  static const double epsilon = std::numeric_limits<double>::min();
  for (int i = 0; i < 4; i++)
    if (u[i] < epsilon)
      u[i] = epsilon;

  // Box muller transform code adapted from
  // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform

  static const double two_pi = 2.0 * Utils::pi();
  return {sqrt(-2.0 * log(u[0])) * cos(two_pi * u[1]),
          sqrt(-2.0 * log(u[0])) * sin(two_pi * u[1]),
          sqrt(-2.0 * log(u[2])) * cos(two_pi * u[3])};
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
 */
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
 */
void mpi_random_seed(int cnt, std::vector<int> &seeds);

/**
 * @brief Gets a string representation of the state of all the nodes.
 */
std::string mpi_random_get_stat();

/**
 * @brief Set the seeds on all the node to the state represented
 *        by the string.
 * The string representation must be one that was returned by
 * @ref mpi_random_get_stat.
 */
void mpi_random_set_stat(const std::vector<std::string> &stat);

/**
 * @brief Get the state size of the random number generator
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

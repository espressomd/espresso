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

#include <Random123/philox.h>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/u32_to_u64.hpp>
#include <utils/uniform.hpp>

#include <random>
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
 * @brief Generator for random uniform noise.
 *
 * Mean = 0, variance = 1 / 12.
 * This uses the Philox PRNG, the state is controlled
 * by the counter, the salt and two keys.
 * If any of the keys and salt differ, the noise is
 * not correlated between two calls along the same counter
 * sequence.
 *
 * @tparam salt RNG salt
 * @tparam N    Size of the noise vector
 *
 * @return Vector of uniform random numbers.
 */
template <RNGSalt salt, size_t N = 3,
          std::enable_if_t<(N > 1) and (N <= 4), int> = 0>
auto noise_uniform(uint64_t counter, int key1, int key2 = 0) {

  auto const integers = philox_4_uint64s<salt>(counter, key1, key2);
  Utils::VectorXd<N> noise{};
  std::transform(integers.begin(), integers.begin() + N, noise.begin(),
                 [](size_t value) { return Utils::uniform(value) - 0.5; });
  return noise;
}

template <RNGSalt salt, size_t N, std::enable_if_t<N == 1, int> = 0>
auto noise_uniform(uint64_t counter, int key1, int key2 = 0) {

  auto const integers = philox_4_uint64s<salt>(counter, key1, key2);
  return Utils::uniform(integers[0]) - 0.5;
}

/** @brief Generator for Gaussian noise.
 *
 * Mean = 0, standard deviation = 1.0.
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
 * @tparam salt decorrelates different thermostat types
 *
 * @return Vector of Gaussian random numbers.
 *
 */
template <RNGSalt salt, size_t N = 3,
          class = std::enable_if_t<(N >= 1) and (N <= 4)>>
auto noise_gaussian(uint64_t counter, int key1, int key2 = 0) {

  auto const integers = philox_4_uint64s<salt>(counter, key1, key2);
  static const double epsilon = std::numeric_limits<double>::min();

  constexpr size_t M = (N <= 2) ? 2 : 4;
  Utils::VectorXd<M> u{};
  std::transform(integers.begin(), integers.begin() + M, u.begin(),
                 [](size_t value) {
                   auto u = Utils::uniform(value);
                   return (u < epsilon) ? epsilon : u;
                 });

  // Box-Muller transform code adapted from
  // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
  // optimizations: the modulo is cached (logarithms are expensive), the
  // sin/cos are evaluated simultaneously by gcc or separately by Clang
  Utils::VectorXd<N> noise{};
  constexpr double two_pi = 2.0 * Utils::pi();
  auto const modulo = sqrt(-2.0 * log(u[0]));
  auto const angle = two_pi * u[1];
  noise[0] = modulo * cos(angle);
  if (N > 1) {
    noise[1] = modulo * sin(angle);
  }
  if (N > 2) {
    auto const modulo = sqrt(-2.0 * log(u[2]));
    auto const angle = two_pi * u[3];
    noise[2] = modulo * cos(angle);
    if (N > 3) {
      noise[3] = modulo * sin(angle);
    }
  }
  return noise;
}

/** Mersenne Twister with warmup.
 *  The first 100'000 values of Mersenne Twister generators are often heavily
 *  correlated @cite panneton06a. This utility function discards the first
 *  1'000'000 values.
 *
 *  @param seed RNG seed
 */
template <typename T> std::mt19937 mt19937(T &&seed) {
  std::mt19937 generator(seed);
  generator.discard(1'000'000);
  return generator;
}

} // namespace Random

#endif

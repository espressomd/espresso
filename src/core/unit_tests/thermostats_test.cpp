/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

/* Unit tests for thermostats. */

#define BOOST_TEST_MODULE Thermostats test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "config/config.hpp"
#include "random.hpp"
#include "random_test.hpp"
#include "thermostat.hpp"
#include "thermostats/brownian_inline.hpp"
#include "thermostats/langevin_inline.hpp"
#include "thermostats/npt_inline.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <tuple>

// multiply by 100 because BOOST_CHECK_CLOSE takes a percentage tolerance,
// and by 8 to account for error accumulation in thermostat functions
auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();

Particle particle_factory() {
  Particle p{};
  p.id() = 0;
  p.force() = {1.0, 2.0, 3.0};
#ifdef ROTATION
  p.torque() = 4.0 * p.force();
#endif
  return p;
}

template <typename T, typename... Args> T thermostat_factory(Args... args) {
  T thermostat = T{};
#ifdef PARTICLE_ANISOTROPY
  thermostat.gamma = {3.0, 5.0, 7.0};
#else
  thermostat.gamma = 2.0;
#endif
  thermostat.gamma_rotation = 3.0 * thermostat.gamma;
  thermostat.rng_initialize(0);
  thermostat.recalc_prefactors(args...);
  return thermostat;
}

BOOST_AUTO_TEST_CASE(test_brownian_dynamics) {
  constexpr double time_step = 0.1;
  constexpr double kT = 3.0;
  auto const brownian = thermostat_factory<BrownianThermostat>(kT);
  auto const dispersion =
      hadamard_division(particle_factory().force(), brownian.gamma);

  /* check translation */
  {
    auto const p = particle_factory();
    auto const ref = time_step * dispersion;
    auto const out = bd_drag(brownian.gamma, p, time_step);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

  /* check translational velocity */
  {
    auto const p = particle_factory();
    auto const ref = dispersion;
    auto const out = bd_drag_vel(brownian.gamma, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

  /* check walk translation */
  {
    auto const p = particle_factory();
    auto const sigma = sqrt(brownian.gamma / (2.0 * kT));
    auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_WALK>(0, 0, 0);
    auto const ref = hadamard_division(noise, sigma) * sqrt(time_step);
    auto const out = bd_random_walk(brownian, p, time_step, kT);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

  /* check walk translational velocity */
  {
    auto const p = particle_factory();
    auto const sigma = sqrt(kT);
    auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_INC>(0, 0, 0);
    auto const ref = sigma * noise / sqrt(p.mass());
    auto const out = bd_random_walk_vel(brownian, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

#ifdef ROTATION
  auto const dispersion_rotation =
      hadamard_division(particle_factory().torque(), brownian.gamma_rotation);

  /* check rotation */
  {
    auto p = particle_factory();
    p.set_cannot_rotate_all_axes();
    p.set_can_rotate_around(0, true);
    auto const phi = time_step * dispersion_rotation[0];
    auto const out = bd_drag_rot(brownian.gamma_rotation, p, time_step);
    BOOST_CHECK_CLOSE(out[0], std::cos(phi / 2), tol);
    BOOST_CHECK_CLOSE(out[1], std::sin(phi / 2), tol);
    BOOST_CHECK_CLOSE(out[2], 0.0, tol);
    BOOST_CHECK_CLOSE(out[3], 0.0, tol);
  }

  /* check rotational velocity */
  {
    auto p = particle_factory();
    p.set_can_rotate_all_axes();
    auto const ref = dispersion_rotation;
    auto const out = bd_drag_vel_rot(brownian.gamma_rotation, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

  /* check walk rotation */
  {
    auto p = particle_factory();
    p.set_cannot_rotate_all_axes();
    p.set_can_rotate_around(0, true);
    auto const sigma = sqrt(brownian.gamma_rotation / (2.0 * kT));
    auto const noise =
        Random::noise_gaussian<RNGSalt::BROWNIAN_ROT_INC>(0, 0, 0);
    auto const phi = hadamard_division(noise, sigma)[0] * sqrt(time_step);
    auto const out = bd_random_walk_rot(brownian, p, time_step, kT);
    BOOST_CHECK_CLOSE(out[0], std::cos(phi / 2), tol);
    BOOST_CHECK_CLOSE(out[1], std::sin(phi / 2), tol);
    BOOST_CHECK_CLOSE(out[2], 0.0, tol);
    BOOST_CHECK_CLOSE(out[3], 0.0, tol);
  }

  /* check walk rotational velocity */
  {
    auto p = particle_factory();
    p.set_can_rotate_all_axes();
    auto const sigma = sqrt(kT);
    auto const noise =
        Random::noise_gaussian<RNGSalt::BROWNIAN_ROT_WALK>(0, 0, 0);
    auto const ref = sigma * noise;
    auto const out = bd_random_walk_vel_rot(brownian, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }
#endif // ROTATION
}

BOOST_AUTO_TEST_CASE(test_langevin_dynamics) {
  constexpr double time_step = 0.1;
  constexpr double kT = 3.0;
  auto const langevin = thermostat_factory<LangevinThermostat>(kT, time_step);
  auto const prefactor_squared = 24.0 * kT / time_step;

  /* check translation */
  {
    auto p = particle_factory();
    p.v() = {1.0, 2.0, 3.0};
    auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN>(0, 0, 0);
    auto const pref = sqrt(prefactor_squared * langevin.gamma);
    auto const ref = hadamard_product(-langevin.gamma, p.v()) +
                     hadamard_product(pref, noise);
    auto const out = friction_thermo_langevin(langevin, p, time_step, kT);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

#ifdef ROTATION
  /* check rotation */
  {
    auto p = particle_factory();
    p.omega() = {1.0, 2.0, 3.0};
    auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(0, 0, 0);
    auto const pref = sqrt(prefactor_squared * langevin.gamma_rotation);
    auto const ref = hadamard_product(-langevin.gamma_rotation, p.omega()) +
                     hadamard_product(pref, noise);
    auto const out =
        friction_thermo_langevin_rotation(langevin, p, time_step, kT);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }
#endif // ROTATION
}

BOOST_AUTO_TEST_CASE(test_noise_statistics) {
  constexpr double time_step = 1.0;
  constexpr double kT = 2.0;
  constexpr std::size_t const sample_size = 10'000;
  auto thermostat = thermostat_factory<LangevinThermostat>(kT, time_step);
  auto p1 = particle_factory();
  auto p2 = particle_factory();
  p1.id() = 0;
  p2.id() = 1;

  auto const correlation = std::get<3>(noise_statistics(
      [&p1, &p2, &thermostat]() -> std::array<VariantVectorXd, 3> {
        thermostat.rng_increment();
        return {{friction_thermo_langevin(thermostat, p1, time_step, kT),
                 -friction_thermo_langevin(thermostat, p1, time_step, kT),
                 friction_thermo_langevin(thermostat, p2, time_step, kT)}};
      },
      sample_size));
  for (std::size_t i = 0; i < correlation.size(); ++i) {
    for (std::size_t j = i; j < correlation.size(); ++j) {
      double expected;
      if (i == j) {
        expected = 1.0;
      } else if (i < 3 and j == i + 3) {
        expected = -1.0;
      } else {
        expected = 0.0;
      }
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, expected, 3e-2));
    }
  }
}

BOOST_AUTO_TEST_CASE(test_brownian_randomness) {
  constexpr double time_step = 1.0;
  constexpr double kT = 2.0;
  constexpr std::size_t const sample_size = 10'000;
  auto thermostat = thermostat_factory<BrownianThermostat>(kT);
  auto p = particle_factory();
#ifdef ROTATION
  p.set_can_rotate_all_axes();
  constexpr std::size_t N = 4;
#else
  constexpr std::size_t N = 2;
#endif

  auto const correlation = std::get<3>(noise_statistics(
      [&p, &thermostat]() -> std::array<VariantVectorXd, N> {
        thermostat.rng_increment();
        return {{
            bd_random_walk(thermostat, p, time_step, kT),
            bd_random_walk_vel(thermostat, p),
#ifdef ROTATION
            bd_random_walk_rot(thermostat, p, time_step, kT),
            bd_random_walk_vel_rot(thermostat, p),
#endif
        }};
      },
      sample_size));
  for (std::size_t i = 0; i < correlation.size(); ++i) {
    for (std::size_t j = i + 1; j < correlation.size(); ++j) {
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, 0.0, 4e-2));
    }
  }
}

BOOST_AUTO_TEST_CASE(test_langevin_randomness) {
  constexpr double time_step = 1.0;
  constexpr double kT = 2.0;
  constexpr std::size_t const sample_size = 10'000;
  auto thermostat = thermostat_factory<LangevinThermostat>(kT, time_step);
  auto p = particle_factory();
#ifdef ROTATION
  constexpr std::size_t N = 2;
#else
  constexpr std::size_t N = 1;
#endif

  auto const correlation = std::get<3>(noise_statistics(
      [&p, &thermostat]() -> std::array<VariantVectorXd, N> {
        thermostat.rng_increment();
        return {{
            friction_thermo_langevin(thermostat, p, time_step, kT),
#ifdef ROTATION
            friction_thermo_langevin_rotation(thermostat, p, time_step, kT),
#endif
        }};
      },
      sample_size));
  for (std::size_t i = 0; i < correlation.size(); ++i) {
    for (std::size_t j = i + 1; j < correlation.size(); ++j) {
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, 0.0, 3e-2));
    }
  }
}

#ifdef NPT
BOOST_AUTO_TEST_CASE(test_npt_iso_randomness) {
  extern int thermo_switch;
  thermo_switch |= THERMO_NPT_ISO;
  constexpr double time_step = 1.0;
  constexpr double kT = 2.0;
  constexpr std::size_t const sample_size = 10'000;
  IsotropicNptThermostat thermostat{};
  thermostat.rng_initialize(0);
  thermostat.gamma0 = 2.0;
  thermostat.gammav = 0.1;
  thermostat.recalc_prefactors(kT, 1.0, time_step);
  auto p = particle_factory();

  auto const correlation = std::get<3>(noise_statistics(
      [&p, &thermostat]() -> std::array<VariantVectorXd, 3> {
        thermostat.rng_increment();
        return {{
            friction_therm0_nptiso<1>(thermostat, p.v(), 0),
            friction_therm0_nptiso<2>(thermostat, p.v(), 0),
            friction_thermV_nptiso(thermostat, 1.5),
        }};
      },
      sample_size));
  for (std::size_t i = 0; i < correlation.size(); ++i) {
    for (std::size_t j = i + 1; j < correlation.size(); ++j) {
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, 0.0, 2e-2));
    }
  }
}
#endif // NPT

BOOST_AUTO_TEST_CASE(test_predicate) {
  std::vector<std::vector<double>> const correlation = {{1., 0.}, {0., 1.}};
  auto const is_true = correlation_almost_equal(correlation, 0, 0, 1., 1e-10);
  auto const is_false = correlation_almost_equal(correlation, 0, 1, 1., 1e-10);
  BOOST_REQUIRE(is_true);
  BOOST_REQUIRE(!is_false);
  BOOST_CHECK_EQUAL(is_true.message(), "");
  BOOST_CHECK_EQUAL(is_false.message(),
                    "The correlation coefficient M[0][1]{0} "
                    "differs from 1 by 1 (> 1e-10)");
}

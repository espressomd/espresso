/*
 * Copyright (C) 2019 The ESPResSo project
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
#include <cmath>
#include <limits>

#include <utils/Vector.hpp>

#include "Particle.hpp"
#include "integrators/brownian_inline.hpp"
#include "integrators/langevin_inline.hpp"
#include "random_test.hpp"
#include "thermostat.hpp"

extern double time_step;
extern double temperature;

// multiply by 100 because BOOST_CHECK_CLOSE takes a percentage tolerance,
// and by 6 to account for error accumulation in thermostat functions
auto const tol = 6 * 100 * std::numeric_limits<double>::epsilon();

Particle particle_factory() {
  Particle p{};
  p.p.identity = 0;
  p.f.f = {1.0, 2.0, 3.0};
#ifdef ROTATION
  p.f.torque = 4.0 * p.f.f;
#endif
  return p;
}

template <typename T> T thermostat_factory() {
  T thermostat = T{};
#ifdef PARTICLE_ANISOTROPY
  thermostat.gamma = {3.0, 5.0, 7.0};
#else
  thermostat.gamma = 2.0;
#endif
  thermostat.gamma_rotation = 3.0 * thermostat.gamma;
  thermostat.recalc_prefactors();
  thermostat.rng_initialize(0);
  return thermostat;
}

BOOST_AUTO_TEST_CASE(test_brownian_dynamics) {
  time_step = 0.1;
  temperature = 3.0;
  auto const brownian = thermostat_factory<BrownianThermostat>();
  auto const dispersion =
      hadamard_division(particle_factory().f.f, brownian.gamma);

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
    auto const sigma = sqrt(brownian.gamma / (2.0 * temperature));
    auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_WALK>(0, 0);
    auto const ref = hadamard_division(noise, sigma) * sqrt(time_step);
    auto const out = bd_random_walk(brownian, p, time_step);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

  /* check walk translational velocity */
  {
    auto const p = particle_factory();
    auto const sigma = sqrt(temperature);
    auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_INC>(0, 0);
    auto const ref = sigma * noise / sqrt(p.p.mass);
    auto const out = bd_random_walk_vel(brownian, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

#ifdef ROTATION
  auto const dispersion_rotation =
      hadamard_division(particle_factory().f.torque, brownian.gamma_rotation);

  /* check rotation */
  {
    auto p = particle_factory();
    p.p.rotation = ROTATION_X;
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
    p.p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;
    auto const ref = dispersion_rotation;
    auto const out = bd_drag_vel_rot(brownian.gamma_rotation, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

  /* check walk rotation */
  {
    auto p = particle_factory();
    p.p.rotation = ROTATION_X;
    auto const sigma = sqrt(brownian.gamma_rotation / (2.0 * temperature));
    auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_ROT_INC>(0, 0);
    auto const phi = hadamard_division(noise, sigma)[0] * sqrt(time_step);
    auto const out = bd_random_walk_rot(brownian, p, time_step);
    BOOST_CHECK_CLOSE(out[0], std::cos(phi / 2), tol);
    BOOST_CHECK_CLOSE(out[1], std::sin(phi / 2), tol);
    BOOST_CHECK_CLOSE(out[2], 0.0, tol);
    BOOST_CHECK_CLOSE(out[3], 0.0, tol);
  }

  /* check walk rotational velocity */
  {
    auto p = particle_factory();
    p.p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;
    auto const sigma = sqrt(temperature);
    auto const noise = Random::noise_gaussian<RNGSalt::BROWNIAN_ROT_WALK>(0, 0);
    auto const ref = sigma * noise;
    auto const out = bd_random_walk_vel_rot(brownian, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }
#endif // ROTATION
}

BOOST_AUTO_TEST_CASE(test_langevin_dynamics) {
  time_step = 0.1;
  temperature = 3.0;
  auto const langevin = thermostat_factory<LangevinThermostat>();
  auto const prefactor_squared = 24.0 * temperature / time_step;

  /* check translation */
  {
    auto p = particle_factory();
    p.m.v = {1.0, 2.0, 3.0};
    auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN>(0, 0);
    auto const pref = sqrt(prefactor_squared * langevin.gamma);
    auto const ref = hadamard_product(-langevin.gamma, p.m.v) +
                     hadamard_product(pref, noise);
    auto const out = friction_thermo_langevin(langevin, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }

#ifdef ROTATION
  /* check rotation */
  {
    auto p = particle_factory();
    p.m.omega = {1.0, 2.0, 3.0};
    auto const noise = Random::noise_uniform<RNGSalt::LANGEVIN_ROT>(0, 0);
    auto const pref = sqrt(prefactor_squared * langevin.gamma_rotation);
    auto const ref = hadamard_product(-langevin.gamma_rotation, p.m.omega) +
                     hadamard_product(pref, noise);
    auto const out = friction_thermo_langevin_rotation(langevin, p);
    BOOST_CHECK_CLOSE(out[0], ref[0], tol);
    BOOST_CHECK_CLOSE(out[1], ref[1], tol);
    BOOST_CHECK_CLOSE(out[2], ref[2], tol);
  }
#endif // ROTATION
}

BOOST_AUTO_TEST_CASE(test_noise_statistics) {
  time_step = 1.0;
  temperature = 2.0;
  constexpr size_t const sample_size = 1'000'000;
  auto thermostat = thermostat_factory<LangevinThermostat>();
  auto p1 = particle_factory();
  auto p2 = particle_factory();
  p1.p.identity = 0;
  p1.p.identity = 1;

  auto const correlation = std::get<3>(noise_statistics(
      std::function<std::vector<VariantVectorXd>()>(
          [&p1, &p2, &thermostat]() -> std::vector<VariantVectorXd> {
            thermostat.rng_increment();
            return {{friction_thermo_langevin(thermostat, p1),
                     -friction_thermo_langevin(thermostat, p1),
                     friction_thermo_langevin(thermostat, p2)}};
          }),
      sample_size));
  for (size_t i = 0; i < correlation.size(); ++i) {
    for (size_t j = i; j < correlation.size(); ++j) {
      double expected;
      if (i == j) {
        expected = 1.0;
      } else if (i < 3 and j == i + 3) {
        expected = -1.0;
      } else {
        expected = 0.0;
      }
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, expected, 2e-3));
    }
  }
}

BOOST_AUTO_TEST_CASE(test_brownian_randomness) {
  time_step = 1.0;
  temperature = 2.0;
  constexpr size_t const sample_size = 500'000;
  auto thermostat = thermostat_factory<BrownianThermostat>();
  auto p = particle_factory();
#ifdef ROTATION
  p.p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;
#endif

  auto const correlation = std::get<3>(
      noise_statistics(std::function<std::vector<VariantVectorXd>()>(
                           [&p, &thermostat]() -> std::vector<VariantVectorXd> {
                             thermostat.rng_increment();
                             return {{
                                 bd_random_walk(thermostat, p, time_step),
                                 bd_random_walk_vel(thermostat, p),
#ifdef ROTATION
                                 bd_random_walk_rot(thermostat, p, time_step),
                                 bd_random_walk_vel_rot(thermostat, p),
#endif
                             }};
                           }),
                       sample_size));
  for (size_t i = 0; i < correlation.size(); ++i) {
    for (size_t j = i + 1; j < correlation.size(); ++j) {
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, 0.0, 5e-2));
    }
  }
}

BOOST_AUTO_TEST_CASE(test_langevin_randomness) {
  time_step = 1.0;
  temperature = 2.0;
  constexpr size_t const sample_size = 500'000;
  auto thermostat = thermostat_factory<LangevinThermostat>();
  auto p = particle_factory();

  auto const correlation = std::get<3>(noise_statistics(
      std::function<std::vector<VariantVectorXd>()>(
          [&p, &thermostat]() -> std::vector<VariantVectorXd> {
            thermostat.rng_increment();
            return {{
                friction_thermo_langevin(thermostat, p),
#ifdef ROTATION
                friction_thermo_langevin_rotation(thermostat, p),
#endif
            }};
          }),
      sample_size));
  for (size_t i = 0; i < correlation.size(); ++i) {
    for (size_t j = i + 1; j < correlation.size(); ++j) {
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, 0.0, 5e-2));
    }
  }
}

#ifdef NPT
BOOST_AUTO_TEST_CASE(test_npt_iso_randomness) {
  extern int thermo_switch;
  thermo_switch |= THERMO_NPT_ISO;
  time_step = 1.0;
  temperature = 2.0;
  constexpr size_t const sample_size = 500'000;
  IsotropicNptThermostat thermostat{};
  thermostat.rng_initialize(0);
  thermostat.gamma0 = 2.0;
  thermostat.gammav = 0.1;
  thermostat.recalc_prefactors(1.0);
  auto p = particle_factory();

  auto const correlation = std::get<3>(noise_statistics(
      std::function<std::vector<VariantVectorXd>()>(
          [&p, &thermostat]() -> std::vector<VariantVectorXd> {
            thermostat.rng_increment();
            return {{
                friction_therm0_nptiso<1>(thermostat, p.m.v, 0),
                friction_therm0_nptiso<2>(thermostat, p.m.v, 0),
                friction_thermV_nptiso(thermostat, 1.5),
            }};
          }),
      sample_size));
  for (size_t i = 0; i < correlation.size(); ++i) {
    for (size_t j = i + 1; j < correlation.size(); ++j) {
      BOOST_CHECK(correlation_almost_equal(correlation, i, j, 0.0, 5e-2));
    }
  }
}
#endif // NPT

/*
 * Copyright (C) 2021-2023 The ESPResSo project
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
#define BOOST_TEST_MODULE Walberla kernels
#define BOOST_TEST_DYN_LINK

#include "config/config.hpp"

#ifdef WALBERLA

#include <boost/test/unit_test.hpp>

#include "../src/lattice_boltzmann/generated_kernels/Dynamic_UBB_double_precision.h"
#include "../src/lattice_boltzmann/generated_kernels/Dynamic_UBB_single_precision.h"
#include "../src/lattice_boltzmann/generated_kernels/FieldAccessorsDoublePrecision.h"
#include "../src/lattice_boltzmann/generated_kernels/FieldAccessorsSinglePrecision.h"

#include <walberla_bridge/utils/walberla_utils.hpp>

#include <utils/Vector.hpp>

#include <cmath>
#include <limits>

bool operator!=(
    const walberla::lbm::Dynamic_UBB_single_precision::IndexInfo &lhs,
    const walberla::lbm::Dynamic_UBB_single_precision::IndexInfo &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    const walberla::lbm::Dynamic_UBB_double_precision::IndexInfo &lhs,
    const walberla::lbm::Dynamic_UBB_double_precision::IndexInfo &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    const walberla::lbm::Dynamic_UBB_single_precision::IndexVectors &lhs,
    const walberla::lbm::Dynamic_UBB_single_precision::IndexVectors &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    const walberla::lbm::Dynamic_UBB_double_precision::IndexVectors &lhs,
    const walberla::lbm::Dynamic_UBB_double_precision::IndexVectors &rhs) {
  return not(lhs == rhs);
}

BOOST_AUTO_TEST_CASE(dynamic_ubb) {
  using Dynamic_UBB_f = walberla::lbm::Dynamic_UBB_single_precision;
  using Dynamic_UBB_d = walberla::lbm::Dynamic_UBB_double_precision;

  // check IndexInfo
  auto vel1_f = Dynamic_UBB_f::IndexInfo(1, 2, 3, 0);
  auto vel2_f = Dynamic_UBB_f::IndexInfo(1, 2, 3, 0);
  auto vel1_d = Dynamic_UBB_d::IndexInfo(1, 2, 3, 0);
  auto vel2_d = Dynamic_UBB_d::IndexInfo(1, 2, 3, 0);
  vel1_f.vel_0 = vel2_f.vel_0 = 1.0f;
  vel1_f.vel_1 = vel2_f.vel_1 = 2.0f;
  vel1_f.vel_2 = vel2_f.vel_2 = 3.0f;
  vel1_d.vel_0 = vel2_d.vel_0 = 1.0;
  vel1_d.vel_1 = vel2_d.vel_1 = 2.0;
  vel1_d.vel_2 = vel2_d.vel_2 = 3.0;
  BOOST_TEST((vel1_f == vel2_f));
  BOOST_TEST((vel1_d == vel2_d));
  vel2_f.vel_2 += 1.0f;
  vel2_d.vel_2 += 1.0;
  BOOST_TEST((vel1_f != vel2_f));
  BOOST_TEST((vel1_d != vel2_d));
  vel2_f.vel_2 += 1.0f;
  vel2_d.vel_2 += 1.0;
  BOOST_TEST((vel1_f != vel2_f));
  BOOST_TEST((vel1_d != vel2_d));

  // check IndexVector
  auto vec1_f = Dynamic_UBB_f::IndexVectors();
  auto vec2_f = Dynamic_UBB_f::IndexVectors();
  auto vec1_d = Dynamic_UBB_d::IndexVectors();
  auto vec2_d = Dynamic_UBB_d::IndexVectors();
  vec1_f.indexVector(Dynamic_UBB_f::IndexVectors::Type::ALL).push_back(vel1_f);
  vec2_f.indexVector(Dynamic_UBB_f::IndexVectors::Type::ALL).push_back(vel1_f);
  vec1_d.indexVector(Dynamic_UBB_d::IndexVectors::Type::ALL).push_back(vel1_d);
  vec2_d.indexVector(Dynamic_UBB_d::IndexVectors::Type::ALL).push_back(vel1_d);
  BOOST_TEST((vec1_f == vec2_f));
  BOOST_TEST((vec1_d == vec2_d));
  vec1_f.indexVector(Dynamic_UBB_f::IndexVectors::Type::ALL).push_back(vel1_f);
  vec2_f.indexVector(Dynamic_UBB_f::IndexVectors::Type::ALL).push_back(vel2_f);
  vec1_d.indexVector(Dynamic_UBB_d::IndexVectors::Type::ALL).push_back(vel1_d);
  vec2_d.indexVector(Dynamic_UBB_d::IndexVectors::Type::ALL).push_back(vel2_d);
  BOOST_TEST((vec1_f != vec2_f));
  BOOST_TEST((vec1_d != vec2_d));
}

static auto clamp_zero(double value) {
  auto constexpr epsilon = std::numeric_limits<float>::epsilon();
  return (std::abs(value) < 5.f * epsilon) ? 0. : value;
}

static auto clamp_zero(float value) {
  auto constexpr epsilon = std::numeric_limits<float>::epsilon();
  return (std::abs(value) < 5.f * epsilon) ? 0.f : value;
}

BOOST_AUTO_TEST_CASE(macroscopic_accessor_equilibrium_distribution) {
  using namespace walberla::stencil;
  using namespace walberla::lbm::accessor;

  auto const x = std::sqrt(1. / 3.);
  auto const u = Utils::Vector3d::broadcast(x);
  auto const u_f = walberla::to_vector3<float>(u);
  auto const u_d = walberla::to_vector3<double>(u);
  auto const rho_f = 0.2f;
  auto const rho_d = 0.2;
  auto const tol_f = 100.f * 5e-7f;
  auto const tol_d = 100. * 5e-9;

  {
    auto const direction = Direction::C;
    auto const ref_d = rho_d * 1. / 3. * (1. - u.norm2());
    auto const ref_f = static_cast<float>(ref_d);
    auto const pop_f = EquilibriumDistribution::get(direction, u_f, rho_f);
    auto const pop_d = EquilibriumDistribution::get(direction, u_d, rho_d);
    BOOST_CHECK_CLOSE(clamp_zero(pop_f), clamp_zero(ref_f), tol_f);
    BOOST_CHECK_CLOSE(clamp_zero(pop_d), clamp_zero(ref_d), tol_d);
  }
  {
    auto const ref_d = rho_d * (1. / 18. - 1. / 6. * (x * x - x));
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::N, Direction::E, Direction::T}) {
      auto const pop_f = EquilibriumDistribution::get(direction, u_f, rho_f);
      auto const pop_d = EquilibriumDistribution::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(clamp_zero(pop_f), clamp_zero(ref_f), tol_f);
      BOOST_CHECK_CLOSE(clamp_zero(pop_d), clamp_zero(ref_d), tol_d);
    }
  }
  {
    auto const ref_d = rho_d * (1. / 18. - 1. / 6. * (x + x * x));
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::S, Direction::W, Direction::B}) {
      auto const pop_f = EquilibriumDistribution::get(direction, u_f, rho_f);
      auto const pop_d = EquilibriumDistribution::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(clamp_zero(pop_f), clamp_zero(ref_f), tol_f);
      BOOST_CHECK_CLOSE(clamp_zero(pop_d), clamp_zero(ref_d), tol_d);
    }
  }
  {
    auto const ref_d = rho_d * (1. / 36. - 1. / 12. * x * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::NW, Direction::SE, Direction::TS,
                                 Direction::TW, Direction::BN, Direction::BE}) {
      auto const pop_f = EquilibriumDistribution::get(direction, u_f, rho_f);
      auto const pop_d = EquilibriumDistribution::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(clamp_zero(pop_f), clamp_zero(ref_f), tol_f);
      BOOST_CHECK_CLOSE(clamp_zero(pop_d), clamp_zero(ref_d), tol_d);
    }
  }
  {
    auto const ref_d = rho_d * (1. / 36. + 5. / 12. * x * x + 2. / 12. * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::NE, Direction::TN, Direction::TE}) {
      auto const pop_f = EquilibriumDistribution::get(direction, u_f, rho_f);
      auto const pop_d = EquilibriumDistribution::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(clamp_zero(pop_f), clamp_zero(ref_f), tol_f);
      BOOST_CHECK_CLOSE(clamp_zero(pop_d), clamp_zero(ref_d), tol_d);
    }
  }
  {
    auto const ref_d = rho_d * (1. / 36. + 5. / 12. * x * x - 2. / 12. * x);
    auto const ref_f = static_cast<float>(ref_d);
    for (auto const direction : {Direction::SW, Direction::BS, Direction::BW}) {
      auto const pop_f = EquilibriumDistribution::get(direction, u_f, rho_f);
      auto const pop_d = EquilibriumDistribution::get(direction, u_d, rho_d);
      BOOST_CHECK_CLOSE(clamp_zero(pop_f), clamp_zero(ref_f), tol_f);
      BOOST_CHECK_CLOSE(clamp_zero(pop_d), clamp_zero(ref_d), tol_d);
    }
  }
}

#endif

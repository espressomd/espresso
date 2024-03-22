/*
 * Copyright (C) 2023 The ESPResSo project
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
#define BOOST_TEST_MODULE waLBerla EK kernels
#define BOOST_TEST_DYN_LINK

#include "config/config.hpp"

#ifdef WALBERLA

#include <boost/test/unit_test.hpp>

#include "../src/electrokinetics/generated_kernels/Dirichlet_double_precision.h"
#include "../src/electrokinetics/generated_kernels/Dirichlet_single_precision.h"
#include "../src/electrokinetics/generated_kernels/FixedFlux_double_precision.h"
#include "../src/electrokinetics/generated_kernels/FixedFlux_single_precision.h"

bool operator!=(
    walberla::pystencils::Dirichlet_single_precision::IndexInfo const &lhs,
    walberla::pystencils::Dirichlet_single_precision::IndexInfo const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::Dirichlet_double_precision::IndexInfo const &lhs,
    walberla::pystencils::Dirichlet_double_precision::IndexInfo const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::Dirichlet_single_precision::IndexVectors const &lhs,
    walberla::pystencils::Dirichlet_single_precision::IndexVectors const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::Dirichlet_double_precision::IndexVectors const &lhs,
    walberla::pystencils::Dirichlet_double_precision::IndexVectors const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::FixedFlux_single_precision::IndexInfo const &lhs,
    walberla::pystencils::FixedFlux_single_precision::IndexInfo const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::FixedFlux_double_precision::IndexInfo const &lhs,
    walberla::pystencils::FixedFlux_double_precision::IndexInfo const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::FixedFlux_single_precision::IndexVectors const &lhs,
    walberla::pystencils::FixedFlux_single_precision::IndexVectors const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::pystencils::FixedFlux_double_precision::IndexVectors const &lhs,
    walberla::pystencils::FixedFlux_double_precision::IndexVectors const &rhs) {
  return not(lhs == rhs);
}

BOOST_AUTO_TEST_CASE(dirichlet_density) {
  using Dirichlet_f = walberla::pystencils::Dirichlet_single_precision;
  using Dirichlet_d = walberla::pystencils::Dirichlet_double_precision;

  // check IndexInfo
  auto dens1_f = Dirichlet_f::IndexInfo(1, 2, 3, 0);
  auto dens2_f = Dirichlet_f::IndexInfo(1, 2, 3, 0);
  auto dens1_d = Dirichlet_d::IndexInfo(1, 2, 3, 0);
  auto dens2_d = Dirichlet_d::IndexInfo(1, 2, 3, 0);
  dens1_f.value = dens2_f.value = 2.0f;
  dens1_d.value = dens2_d.value = 2.0;
  BOOST_TEST((dens1_f == dens2_f));
  BOOST_TEST((dens1_d == dens2_d));
  dens2_f.value += 1.0f;
  dens2_d.value += 1.0;
  BOOST_TEST((dens1_f != dens2_f));
  BOOST_TEST((dens1_d != dens2_d));
  dens2_f.value += 1.0f;
  dens2_d.value += 1.0;
  BOOST_TEST((dens1_f != dens2_f));
  BOOST_TEST((dens1_d != dens2_d));

  // check IndexVector
  auto vec1_f = Dirichlet_f::IndexVectors();
  auto vec2_f = Dirichlet_f::IndexVectors();
  auto vec1_d = Dirichlet_d::IndexVectors();
  auto vec2_d = Dirichlet_d::IndexVectors();
  vec1_f.indexVector(Dirichlet_f::IndexVectors::Type::ALL).push_back(dens1_f);
  vec2_f.indexVector(Dirichlet_f::IndexVectors::Type::ALL).push_back(dens1_f);
  vec1_d.indexVector(Dirichlet_d::IndexVectors::Type::ALL).push_back(dens1_d);
  vec2_d.indexVector(Dirichlet_d::IndexVectors::Type::ALL).push_back(dens1_d);
  BOOST_TEST((vec1_f == vec2_f));
  BOOST_TEST((vec1_d == vec2_d));
  vec1_f.indexVector(Dirichlet_f::IndexVectors::Type::ALL).push_back(dens1_f);
  vec2_f.indexVector(Dirichlet_f::IndexVectors::Type::ALL).push_back(dens2_f);
  vec1_d.indexVector(Dirichlet_d::IndexVectors::Type::ALL).push_back(dens1_d);
  vec2_d.indexVector(Dirichlet_d::IndexVectors::Type::ALL).push_back(dens2_d);
  BOOST_TEST((vec1_f != vec2_f));
  BOOST_TEST((vec1_d != vec2_d));
}

BOOST_AUTO_TEST_CASE(dirichlet_flux) {
  using FixedFlux_f = walberla::pystencils::FixedFlux_single_precision;
  using FixedFlux_d = walberla::pystencils::FixedFlux_double_precision;

  // check IndexInfo
  auto flux1_f = FixedFlux_f::IndexInfo(1, 2, 3, 0);
  auto flux2_f = FixedFlux_f::IndexInfo(1, 2, 3, 0);
  auto flux1_d = FixedFlux_d::IndexInfo(1, 2, 3, 0);
  auto flux2_d = FixedFlux_d::IndexInfo(1, 2, 3, 0);
  flux1_f.flux_0 = flux2_f.flux_0 = 1.0f;
  flux1_f.flux_1 = flux2_f.flux_1 = 2.0f;
  flux1_f.flux_2 = flux2_f.flux_2 = 3.0f;
  flux1_d.flux_0 = flux2_d.flux_0 = 1.0;
  flux1_d.flux_1 = flux2_d.flux_1 = 2.0;
  flux1_d.flux_2 = flux2_d.flux_2 = 3.0;
  BOOST_TEST((flux1_f == flux2_f));
  BOOST_TEST((flux1_d == flux2_d));
  flux2_f.flux_2 += 1.0f;
  flux2_d.flux_2 += 1.0;
  BOOST_TEST((flux1_f != flux2_f));
  BOOST_TEST((flux1_d != flux2_d));
  flux2_f.flux_2 += 1.0f;
  flux2_d.flux_2 += 1.0;
  BOOST_TEST((flux1_f != flux2_f));
  BOOST_TEST((flux1_d != flux2_d));

  // check IndexVector
  auto vec1_f = FixedFlux_f::IndexVectors();
  auto vec2_f = FixedFlux_f::IndexVectors();
  auto vec1_d = FixedFlux_d::IndexVectors();
  auto vec2_d = FixedFlux_d::IndexVectors();
  vec1_f.indexVector(FixedFlux_f::IndexVectors::Type::ALL).push_back(flux1_f);
  vec2_f.indexVector(FixedFlux_f::IndexVectors::Type::ALL).push_back(flux1_f);
  vec1_d.indexVector(FixedFlux_d::IndexVectors::Type::ALL).push_back(flux1_d);
  vec2_d.indexVector(FixedFlux_d::IndexVectors::Type::ALL).push_back(flux1_d);
  BOOST_TEST((vec1_f == vec2_f));
  BOOST_TEST((vec1_d == vec2_d));
  vec1_f.indexVector(FixedFlux_f::IndexVectors::Type::ALL).push_back(flux1_f);
  vec2_f.indexVector(FixedFlux_f::IndexVectors::Type::ALL).push_back(flux2_f);
  vec1_d.indexVector(FixedFlux_d::IndexVectors::Type::ALL).push_back(flux1_d);
  vec2_d.indexVector(FixedFlux_d::IndexVectors::Type::ALL).push_back(flux2_d);
  BOOST_TEST((vec1_f != vec2_f));
  BOOST_TEST((vec1_d != vec2_d));
}

#endif

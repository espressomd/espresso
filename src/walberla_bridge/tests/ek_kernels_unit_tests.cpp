/*
 * Copyright (C) 2023-2026 The ESPResSo project
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
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include "../src/electrokinetics/generated_kernels/DirichletDoublePrecision.h"
#include "../src/electrokinetics/generated_kernels/DirichletSinglePrecision.h"
#include "../src/electrokinetics/generated_kernels/FixedFlux_double_precision.h"
#include "../src/electrokinetics/generated_kernels/FixedFlux_single_precision.h"

#include <walberla_bridge/walberla_init.hpp>

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/FlagUID.h>
#include <field/GhostLayerField.h>

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>

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
  using Dirichlet_f = walberla::ek::DirichletSinglePrecision;
  using Dirichlet_d = walberla::ek::DirichletDoublePrecision;

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

/**
 * Run the generated Dirichlet boundary on a field with one ghost layer:
 * the value of the callback is stored in the boundary cells (those of the
 * interior and of the ghost layer), in the order of the cells in memory, and
 * the other cells are left alone.
 */
template <typename Dirichlet, typename FloatType> void check_dirichlet_sweep() {
  using namespace walberla;
  using Field = field::GhostLayerField<FloatType, 1>;
  using Flags = field::FlagField<uint8_t>;
  auto const blocks = blockforest::createUniformBlockGrid(
      uint_t(1), uint_t(1), uint_t(1), uint_t(4), uint_t(3), uint_t(5), 1.,
      uint_t(1), uint_t(1), uint_t(1), false, false, false);
  auto const flags_id = field::addFlagFieldToStorage<Flags>(blocks, "flags", 1);
  auto const field_id = field::addToStorage<Field>(
      blocks, "density", FloatType(-1), field::fzyx, 1);
  auto &block = *blocks->begin();
  auto *flags = block.getData<Flags>(flags_id);
  auto *density = block.getData<Field>(field_id);
  auto const domain = FlagUID("domain");
  auto const boundary = FlagUID("boundary");
  flags->registerFlag(domain);
  flags->registerFlag(boundary);
  for (auto it = flags->beginWithGhostLayerXYZ(); it != flags->end(); ++it) {
    flags->addFlag(it.x(), it.y(), it.z(), flags->getFlag(domain));
  }
  // cells in the interior, in the ghost layer and in a corner of the field
  std::vector<Cell> const cells = {
      Cell(0, 0, 0), Cell(1, 1, 3), Cell(3, 2, 4),   Cell(-1, 1, 2),
      Cell(2, 3, 1), Cell(4, 3, 5), Cell(-1, -1, -1)};
  for (auto const &cell : cells) {
    flags->addFlag(cell.x(), cell.y(), cell.z(), flags->getFlag(boundary));
  }
  auto const value = [](Cell const &cell) {
    return FloatType(1 + 0.5 * cell.x() + 0.25 * cell.y() + 0.125 * cell.z());
  };
  std::function<FloatType(Cell const &,
                          shared_ptr<StructuredBlockForest> const &, IBlock &)>
      callback = [&value](Cell const &cell,
                          shared_ptr<StructuredBlockForest> const &,
                          IBlock &) { return value(cell); };
  Dirichlet dirichlet(blocks, field_id, callback);
  dirichlet.template fillFromFlagField<Flags>(blocks, flags_id, boundary,
                                              domain);
  // one link for each boundary cell, in the order of the cells in memory
  auto const &links = dirichlet.getIndexVector(&block);
  BOOST_REQUIRE_EQUAL(links.size(), cells.size());
  for (std::size_t i = 1; i < links.size(); ++i) {
    auto const key = [&](auto const &l) {
      return std::make_tuple(l.z, l.y, l.x);
    };
    BOOST_CHECK(key(links[i - 1]) < key(links[i]));
  }
  for (auto const &link : links) {
    BOOST_CHECK_EQUAL(link.dir, 0);
    BOOST_CHECK(std::ranges::find(cells, Cell(link.x, link.y, link.z)) !=
                cells.end());
    BOOST_CHECK_EQUAL(link.value, value(Cell(link.x, link.y, link.z)));
  }
  // the interior links are those of the cells of the interior
  BOOST_CHECK_EQUAL(
      dirichlet.getIndexVector(&block, Dirichlet::Type::INNER).size(), 1u);
  BOOST_CHECK_EQUAL(
      dirichlet.getIndexVector(&block, Dirichlet::Type::OUTER).size(),
      cells.size() - 1u);
  dirichlet(&block);
  for (auto it = density->beginWithGhostLayerXYZ(); it != density->end();
       ++it) {
    auto const cell = Cell(it.x(), it.y(), it.z());
    auto const is_boundary = std::ranges::find(cells, cell) != cells.end();
    BOOST_CHECK_EQUAL(*it, is_boundary ? value(cell) : FloatType(-1));
  }
}

BOOST_AUTO_TEST_CASE(dirichlet_sweep) {
  check_dirichlet_sweep<walberla::ek::DirichletSinglePrecision, float>();
  check_dirichlet_sweep<walberla::ek::DirichletDoublePrecision, double>();
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  walberla::mpi_init();
  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

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
#include "../src/electrokinetics/generated_kernels/FixedFluxDoublePrecision.h"
#include "../src/electrokinetics/generated_kernels/FixedFluxSinglePrecision.h"
#include "../src/electrokinetics/generated_kernels/PackInfoScalarDoublePrecision.h"
#include "../src/electrokinetics/generated_kernels/PackInfoScalarSinglePrecision.h"
#include "../src/electrokinetics/reactions/ReactionKernelSelector.hpp"

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <blockforest/communication/UniformBufferedScheme.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/FlagUID.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>
#include <stencil/D3Q27.h>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

bool operator!=(walberla::ek::FixedFluxSinglePrecision::IndexInfo const &lhs,
                walberla::ek::FixedFluxSinglePrecision::IndexInfo const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(walberla::ek::FixedFluxDoublePrecision::IndexInfo const &lhs,
                walberla::ek::FixedFluxDoublePrecision::IndexInfo const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::ek::FixedFluxSinglePrecision::IndexVectors const &lhs,
    walberla::ek::FixedFluxSinglePrecision::IndexVectors const &rhs) {
  return not(lhs == rhs);
}

bool operator!=(
    walberla::ek::FixedFluxDoublePrecision::IndexVectors const &lhs,
    walberla::ek::FixedFluxDoublePrecision::IndexVectors const &rhs) {
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
  using FixedFlux_f = walberla::ek::FixedFluxSinglePrecision;
  using FixedFlux_d = walberla::ek::FixedFluxDoublePrecision;

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

namespace {
/** Species of a reaction: the density field and its precision. */
struct ReactionSpecies {
  std::size_t density_id;
  bool double_precision;
  [[nodiscard]] std::size_t get_density_id() const noexcept {
    return density_id;
  }
  [[nodiscard]] bool is_double_precision() const noexcept {
    return double_precision;
  }
};

/** Reactant with the interface that the reaction kernel selector needs. */
struct Reactant {
  std::shared_ptr<ReactionSpecies> species;
  double stoech_coeff;
  double order;
  [[nodiscard]] auto get_species() const noexcept { return species; }
  [[nodiscard]] double get_stoech_coeff() const noexcept {
    return stoech_coeff;
  }
  [[nodiscard]] double get_order() const noexcept { return order; }
};
} // namespace

/**
 * Run the selected reaction kernel of @c num_reactants reactants on a block
 * with one ghost layer: the densities of the reactants change by
 * <tt>stoech_i * coefficient * prod_j rho_j^order_j</tt> in the cells of the
 * interior (bulk) or in the cells of the index list (indexed, which are in
 * the interior and in the ghost layer), and the other cells are left alone.
 */
template <typename FloatType, bool Indexed>
void check_reaction_kernel(std::size_t const num_reactants) {
  using namespace walberla;
  using Field = field::GhostLayerField<FloatType, 1>;
  namespace Selector = walberla::detail::ReactionKernelSelector;
  constexpr auto arch = lbmpy::Arch::CPU;
  auto const blocks = blockforest::createUniformBlockGrid(
      uint_t(1), uint_t(1), uint_t(1), uint_t(4), uint_t(3), uint_t(5), 1.,
      uint_t(1), uint_t(1), uint_t(1), false, false, false);
  auto &block = *blocks->begin();
  auto const initial = [](std::size_t i, Cell const &cell) {
    return FloatType(1.5 + 0.1 * static_cast<double>(i) +
                     0.02 * (cell.x() + 1) + 0.03 * (cell.y() + 1) +
                     0.01 * (cell.z() + 1));
  };
  auto const coefficient = 0.05;
  std::vector<std::shared_ptr<Reactant>> reactants;
  std::vector<Field *> densities;
  for (std::size_t i = 0; i < num_reactants; ++i) {
    auto const id = field::addToStorage<Field>(
        blocks, "density " + std::to_string(i), FloatType(0), field::fzyx, 1);
    auto *density = block.getData<Field>(id);
    for (auto it = density->beginWithGhostLayerXYZ(); it != density->end();
         ++it) {
      *it = initial(i, Cell(it.x(), it.y(), it.z()));
    }
    densities.push_back(density);
    auto const sign = (i % 2 == 0) ? -1. : 1.;
    reactants.push_back(std::make_shared<Reactant>(
        Reactant{std::make_shared<ReactionSpecies>(ReactionSpecies{
                     std::size_t(id), std::is_same_v<FloatType, double>}),
                 sign * (1. + 0.5 * static_cast<double>(i)),
                 0.5 + 0.5 * static_cast<double>(i)}));
  }
  // cells in the interior, in the ghost layer and in the corners of the field
  std::vector<Cell> const cells = {
      Cell(0, 0, 0), Cell(1, 1, 3), Cell(3, 2, 4),   Cell(-1, 1, 2),
      Cell(2, 3, 1), Cell(4, 3, 5), Cell(-1, -1, -1)};
  auto const reacts = [&](Cell const &cell) {
    if constexpr (Indexed) {
      return std::ranges::find(cells, cell) != cells.end();
    } else {
      return densities[0]->xyzSize().contains(cell);
    }
  };
  std::function<void(IBlock *)> kernel;
  if constexpr (Indexed) {
    Selector::IndexList<arch> index_list(*blocks);
    for (auto const &cell : cells) {
      index_list.getVector(block).emplace_back(cell);
    }
    kernel =
        Selector::get_indexed_kernel<arch>(reactants, coefficient, index_list);
  } else {
    kernel = Selector::get_bulk_kernel<arch>(reactants, coefficient);
  }
  kernel(&block);
  // the tolerance is relative to the terms of the sum, because the densities
  // can cancel to values close to zero
  auto const epsilon = std::is_same_v<FloatType, double> ? 1e-12 : 1e-6;
  for (auto it = densities[0]->beginWithGhostLayerXYZ();
       it != densities[0]->end(); ++it) {
    auto const cell = Cell(it.x(), it.y(), it.z());
    auto rate = coefficient;
    for (std::size_t j = 0; j < num_reactants; ++j) {
      rate *= std::pow(static_cast<double>(initial(j, cell)),
                       reactants[j]->get_order());
    }
    for (std::size_t i = 0; i < num_reactants; ++i) {
      auto const before = static_cast<double>(initial(i, cell));
      auto const change =
          reacts(cell) ? reactants[i]->get_stoech_coeff() * rate : 0.;
      BOOST_CHECK_SMALL(static_cast<double>(densities[i]->get(cell)) -
                            (before + change),
                        epsilon * (std::abs(before) + std::abs(change)));
    }
  }
}

BOOST_AUTO_TEST_CASE(reaction_kernel_selector) {
  for (std::size_t num_reactants = 1; num_reactants <= 5; ++num_reactants) {
    BOOST_TEST_CONTEXT("reactants: " << num_reactants) {
      check_reaction_kernel<double, false>(num_reactants);
      check_reaction_kernel<float, false>(num_reactants);
      check_reaction_kernel<double, true>(num_reactants);
      check_reaction_kernel<float, true>(num_reactants);
    }
  }
}

// The dispatch of the GPU kernels is not unit-tested here, because a .cpp file
// cannot instantiate the GPU traits; it is covered by the CUDA build and by
// the Python tests.
BOOST_AUTO_TEST_CASE(reaction_kernel_selector_size_limit) {
  namespace Selector = walberla::detail::ReactionKernelSelector;
  auto const check_message = [](std::runtime_error const &error) {
    return std::string(error.what()) ==
           "reactions of this size are not implemented!";
  };
  auto const blocks = walberla::blockforest::createUniformBlockGrid(
      walberla::uint_t(1), walberla::uint_t(1), walberla::uint_t(1),
      walberla::uint_t(4), walberla::uint_t(3), walberla::uint_t(5), 1.,
      walberla::uint_t(1), walberla::uint_t(1), walberla::uint_t(1), false,
      false, false);
  Selector::IndexList<lbmpy::Arch::CPU> const index_list(*blocks);
  auto const species =
      std::make_shared<ReactionSpecies>(ReactionSpecies{std::size_t{0}, true});
  for (auto const num_reactants : {std::size_t{0}, std::size_t{6}}) {
    std::vector<std::shared_ptr<Reactant>> const reactants(
        num_reactants, std::make_shared<Reactant>(Reactant{species, 1., 1.}));
    BOOST_CHECK_EXCEPTION(
        Selector::get_bulk_kernel<lbmpy::Arch::CPU>(reactants, 1.),
        std::runtime_error, check_message);
    BOOST_CHECK_EXCEPTION(Selector::get_indexed_kernel<lbmpy::Arch::CPU>(
                              reactants, 1., index_list),
                          std::runtime_error, check_message);
  }
}

/**
 * @brief Communicate the ghost layers of a scalar field with the generic
 * waLBerla pack info and with the generated one, and count the cells of the
 * two fields (with their ghost layers) that differ bitwise.
 */
template <typename PackInfo>
static std::size_t compare_scalar_pack_infos(
    walberla::uint_t n_ghost_layers,
    walberla::blockforest::LocalCommunicationMode local_mode) {
  using namespace walberla;
  using Field = typename PackInfo::Field_T;
  using FT = typename Field::value_type;
  static_assert(std::is_same_v<Field, field::GhostLayerField<FT, 1u>>);
  // 2x2x1 blocks of one process, periodic like the lattice of the bridge
  auto const blocks = blockforest::createUniformBlockGrid(
      uint_t(2), uint_t(2), uint_t(1), uint_t(5), uint_t(4), uint_t(6), 1.,
      uint_t(1), uint_t(1), uint_t(1), true, true, true);
  auto const generic_id = field::addToStorage<Field>(
      blocks, "generic", FT{0}, field::fzyx, n_ghost_layers);
  auto const generated_id = field::addToStorage<Field>(
      blocks, "generated", FT{0}, field::fzyx, n_ghost_layers);
  for (auto &block : *blocks) {
    auto const origin = blocks->getBlockCellBB(block);
    for (auto const id : {generic_id, generated_id}) {
      auto &field = *block.getData<Field>(id);
      for (auto it = field.beginWithGhostLayer(); it != field.end(); ++it) {
        auto const is_ghost = it.x() < 0 or it.y() < 0 or it.z() < 0 or
                              it.x() >= cell_idx_c(field.xSize()) or
                              it.y() >= cell_idx_c(field.ySize()) or
                              it.z() >= cell_idx_c(field.zSize());
        // a different value in every cell of the domain
        *it = is_ghost ? FT(-1)
                       : FT(1) + FT(0.25) * FT(origin.xMin() + it.x()) +
                             FT(0.0625) * FT(origin.yMin() + it.y()) +
                             FT(0.015625) * FT(origin.zMin() + it.z()) *
                                 FT(origin.zMin() + it.z() + 1);
      }
    }
  }
  using Scheme =
      blockforest::communication::UniformBufferedScheme<stencil::D3Q27>;
  Scheme generic_scheme(blocks);
  Scheme generated_scheme(blocks);
  generic_scheme.setLocalMode(local_mode);
  generated_scheme.setLocalMode(local_mode);
  generic_scheme.addPackInfo(
      std::make_shared<field::communication::PackInfo<Field>>(generic_id));
  // all ghost layers, like the generic pack info
  generated_scheme.addPackInfo(
      std::make_shared<PackInfo>(generated_id, n_ghost_layers));
  generic_scheme();
  generated_scheme();
  std::size_t n_differences = 0;
  for (auto &block : *blocks) {
    auto &generic = *block.getData<Field>(generic_id);
    auto &generated = *block.getData<Field>(generated_id);
    auto it_generated = generated.beginWithGhostLayer();
    for (auto it = generic.beginWithGhostLayer(); it != generic.end();
         ++it, ++it_generated) {
      // no ghost cell is left at its initial value
      BOOST_REQUIRE(*it != FT(-1));
      n_differences += std::memcmp(&*it, &*it_generated, sizeof(FT)) != 0;
    }
  }
  return n_differences;
}

BOOST_AUTO_TEST_CASE(pack_info_scalar) {
  using walberla::blockforest::BUFFER;
  using walberla::blockforest::START;
  using PackInfo_d = walberla::pystencils::PackInfoScalarDoublePrecision;
  using PackInfo_f = walberla::pystencils::PackInfoScalarSinglePrecision;
  for (auto const mode : {START, BUFFER}) {
    for (auto const n_ghost_layers : {1u, 2u}) {
      BOOST_CHECK_EQUAL(compare_scalar_pack_infos<PackInfo_d>(
                            walberla::uint_c(n_ghost_layers), mode),
                        0u);
      BOOST_CHECK_EQUAL(compare_scalar_pack_infos<PackInfo_f>(
                            walberla::uint_c(n_ghost_layers), mode),
                        0u);
    }
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  walberla::mpi_init();
  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

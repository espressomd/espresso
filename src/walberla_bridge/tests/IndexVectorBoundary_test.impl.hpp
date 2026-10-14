/*
 * Copyright (C) 2026 The ESPResSo project
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

/**
 * @file
 * Tests of @ref walberla::IndexVectorBoundary with a boundary class that
 * has a trivial kernel, on the host and (with
 * @c ESPRESSO_TEST_BUILD_FOR_GPU) on the device. The kernel doubles the PDF
 * of each link, adds a value to it, and stores the old PDF, the direction
 * and the coordinates of the link as forces: the results of two runs
 * can be compared bitwise, and each force identifies its link.
 */

#define BOOST_TEST_MODULE "index vector boundary test"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include "../src/utils/IndexVectorBoundary.hpp"

#include <walberla_bridge/walberla_init.hpp>

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/math/Vector3.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/GhostLayerField.h>
#include <stencil/D3Q19.h>

#include <waLBerlaDefinitions.h>

#if defined(WALBERLA_BUILD_WITH_CUDA) and defined(ESPRESSO_TEST_BUILD_FOR_GPU)
#define TEST_ON_DEVICE
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/ErrorChecking.h>
#include <gpu/FieldCopy.h>
#include <gpu/GPUField.h>
#include <gpu/GPUWrapper.h>
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <vector>

using namespace walberla;

namespace Testing {

constexpr uint_t gl = 2u;
using HostField = field::GhostLayerField<double, 19>;
using Flags = field::FlagField<uint8_t>;
using Stencil = stencil::D3Q19;

struct IndexInfo {
  int32_t x;
  int32_t y;
  int32_t z;
  int32_t dir;
  double value;
  IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_)
      : x(x_), y(y_), z(z_), dir(dir_), value(0.) {}
  bool operator==(IndexInfo const &o) const {
    return x == o.x and y == o.y and z == o.z and dir == o.dir and
           value == o.value;
  }
};

struct ForceStruct {
  double F_0 = 0.;
  double F_1 = 0.;
  double F_2 = 0.;
  bool operator==(ForceStruct const &o) const {
    return F_0 == o.F_0 and F_1 == o.F_1 and F_2 == o.F_2;
  }
};

#if defined(TEST_ON_DEVICE)

using PdfField = gpu::GPUField<double>;
constexpr auto memory = IndexVectorMemory::Device;

__global__ void link_kernel(double *pdfs, IndexInfo const *links,
                            ForceStruct *forces, std::size_t size,
                            std::ptrdiff_t stride_x, std::ptrdiff_t stride_y,
                            std::ptrdiff_t stride_z, std::ptrdiff_t stride_f) {
  auto const i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < size) {
    auto const &link = links[i];
    auto &pdf = pdfs[link.x * stride_x + link.y * stride_y + link.z * stride_z +
                     link.dir * stride_f];
    forces[i].F_0 = pdf;
    forces[i].F_1 = link.dir;
    forces[i].F_2 = link.x + link.y + link.z;
    pdf = 2. * pdf + link.value;
  }
}

#else

using PdfField = HostField;
constexpr auto memory = IndexVectorMemory::Host;

#endif

class Boundary : public IndexVectorBoundary<Boundary, IndexInfo, ForceStruct,
                                            Stencil, Vector3<double>, memory> {
  using Base = IndexVectorBoundary<Boundary, IndexInfo, ForceStruct, Stencil,
                                   Vector3<double>, memory>;

public:
  using PdfField_T = PdfField;
  using Base::Base;

  static void assignValue(IndexInfo &element, Vector3<double> const &value) {
    element.value = value[0] + 0.1 * value[1] + 0.01 * value[2];
  }

#if defined(TEST_ON_DEVICE)
  void linkKernel(PdfField *pdfs, IndexInfo *links, ForceStruct *forces,
                  std::size_t size, gpuStream_t stream) {
    auto const threads = 64u;
    auto const blocks = static_cast<unsigned>((size + threads - 1) / threads);
    link_kernel<<<blocks, threads, 0, stream>>>(
        pdfs->dataAt(0, 0, 0, 0), links, forces, size, pdfs->xStride(),
        pdfs->yStride(), pdfs->zStride(), pdfs->fStride());
    WALBERLA_GPU_CHECK(gpuGetLastError());
  }
#else
  void linkKernel(PdfField *pdfs, IndexInfo *links, ForceStruct *forces,
                  std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
      auto const &link = links[i];
      auto &pdf = pdfs->get(link.x, link.y, link.z, link.dir);
      forces[i].F_0 = pdf;
      forces[i].F_1 = link.dir;
      forces[i].F_2 = link.x + link.y + link.z;
      pdf = 2. * pdf + link.value;
    }
  }
#endif
};

using Forces = Boundary::ForceVector;
using Links = Boundary::IndexVectors;

/** Value callback that records the boundary cells it is called for. */
struct Recorder {
  std::vector<Cell> *cells;
  Vector3<double> operator()(Cell const &cell,
                             shared_ptr<StructuredBlockForest> const &,
                             IBlock &) const {
    cells->push_back(cell);
    return Vector3<double>(double(cell.x()), double(cell.y()),
                           double(cell.z()));
  }
};

/** Two blocks with a flag field and two copies of a PDF field. */
struct Setup {
  shared_ptr<StructuredBlockForest> blocks;
  BlockDataID flags_id;
  BlockDataID host_pdfs_id;
  BlockDataID pdfs_id[2];
  std::mt19937 rng{42u};

  Setup() {
    blocks = blockforest::createUniformBlockGrid(
        uint_t(2), uint_t(1), uint_t(1), uint_t(6), uint_t(5), uint_t(4), 1.,
        uint_t(1), uint_t(1), uint_t(1), false, false, false);
    flags_id = field::addFlagFieldToStorage<Flags>(blocks, "flags", gl);
    host_pdfs_id =
        field::addToStorage<HostField>(blocks, "pdfs", 0., field::fzyx, gl);
    std::uniform_real_distribution<double> distribution(0.02, 0.12);
    for (auto &block : *blocks) {
      auto *pdfs = block.getData<HostField>(host_pdfs_id);
      for (auto it = pdfs->beginWithGhostLayerXYZ(cell_idx_c(gl));
           it != pdfs->end(); ++it) {
        for (uint_t f = 0; f < 19; ++f) {
          it.getF(f) = distribution(rng);
        }
      }
      auto *flags = block.getData<Flags>(flags_id);
      flags->registerFlag(FlagUID("domain"));
      flags->registerFlag(FlagUID("boundary"));
    }
    for (auto &id : pdfs_id) {
      id = make_pdfs();
    }
  }

  BlockDataID make_pdfs() {
#if defined(TEST_ON_DEVICE)
    return gpu::addGPUFieldToStorage<HostField>(blocks, host_pdfs_id, "device");
#else
    return field::addCloneToStorage<HostField>(blocks, host_pdfs_id, "clone");
#endif
  }

  /** Set the domain flag everywhere and the boundary flag at random. */
  void set_flags(double boundary_fraction) {
    std::uniform_real_distribution<double> distribution(0., 1.);
    for (auto &block : *blocks) {
      auto *flags = block.getData<Flags>(flags_id);
      auto const domain = flags->getFlag(FlagUID("domain"));
      auto const boundary = flags->getFlag(FlagUID("boundary"));
      for (auto it = flags->beginWithGhostLayerXYZ(cell_idx_c(gl));
           it != flags->end(); ++it) {
        *it = domain;
        if (distribution(rng) < boundary_fraction) {
          addFlag(it, boundary);
        }
      }
    }
  }

  /** Set the boundary flag of one cell of the first block. */
  void set_flags(std::vector<Cell> const &boundary_cells) {
    for (auto &block : *blocks) {
      auto *flags = block.getData<Flags>(flags_id);
      auto const domain = flags->getFlag(FlagUID("domain"));
      auto const boundary = flags->getFlag(FlagUID("boundary"));
      for (auto it = flags->beginWithGhostLayerXYZ(cell_idx_c(gl));
           it != flags->end(); ++it) {
        *it = domain;
      }
      if (&block == &*blocks->begin()) {
        for (auto const &cell : boundary_cells) {
          addFlag(flags->get(cell), boundary);
        }
      }
    }
  }

  void fill(Boundary &boundary) {
    boundary.fillFromFlagField<Flags>(blocks, flags_id, FlagUID("boundary"),
                                      FlagUID("domain"));
  }

  /** Values of the PDF field @p i, including the ghost layers. */
  std::vector<double> pdfs(int i) {
#if defined(TEST_ON_DEVICE)
    return values(make_host_copy(pdfs_id[i]));
#else
    return values(pdfs_id[i]);
#endif
  }

  /** Values of the PDF field before any sweep. */
  std::vector<double> initial_pdfs() { return values(host_pdfs_id); }

  std::vector<double> values(BlockDataID host_id) {
    std::vector<double> result;
    for (auto &block : *blocks) {
      auto *pdfs = block.getData<HostField>(host_id);
      for (auto it = pdfs->beginWithGhostLayerXYZ(cell_idx_c(gl));
           it != pdfs->end(); ++it) {
        for (uint_t f = 0; f < 19; ++f) {
          result.push_back(it.getF(f));
        }
      }
    }
    return result;
  }

#if defined(TEST_ON_DEVICE)
  BlockDataID make_host_copy(BlockDataID device_id) {
    auto const host_id =
        field::addToStorage<HostField>(blocks, "copy", 0., field::fzyx, gl);
    gpu::fieldCpy<HostField, PdfField>(blocks, host_id, device_id);
    return host_id;
  }
#endif
};

/**
 * Links of type @p type of a block that the flag field implies, in the order
 * of the stencil and then of the cells, computed without the class under test.
 */
std::vector<IndexInfo> expected_links(Setup &setup, IBlock &block,
                                      Links::Type type) {
  std::vector<IndexInfo> links;
  auto *flags = block.getData<Flags>(setup.flags_id);
  auto const domain = flags->getFlag(FlagUID("domain"));
  auto const boundary = flags->getFlag(FlagUID("boundary"));
  auto const is_inner = [&](auto const &it) {
    return it.x() >= 1 and it.x() < cell_idx_c(flags->xSize()) - 1 and
           it.y() >= 1 and it.y() < cell_idx_c(flags->ySize()) - 1 and
           it.z() >= 1 and it.z() < cell_idx_c(flags->zSize()) - 1;
  };
  for (auto d = Stencil::begin(); d != Stencil::end(); ++d) {
    for (auto it = flags->beginWithGhostLayerXYZ(cell_idx_c(gl - 1));
         it != flags->end(); ++it) {
      auto const neighbor =
          Cell(it.x() + stencil::cx[*d], it.y() + stencil::cy[*d],
               it.z() + stencil::cz[*d]);
      auto const selected =
          type == Links::ALL or (type == Links::INNER) == is_inner(it);
      if (selected and isFlagSet(it, domain) and not isFlagSet(it, boundary) and
          isFlagSet(flags->get(neighbor), boundary)) {
        IndexInfo link(it.x(), it.y(), it.z(), int32_t(Stencil::idx[*d]));
        Boundary::assignValue(link, Vector3<double>(double(neighbor.x()),
                                                    double(neighbor.y()),
                                                    double(neighbor.z())));
        links.push_back(link);
      }
    }
  }
  return links;
}

template <typename T> bool is_zero(std::vector<T> const &forces) {
  return std::ranges::all_of(forces,
                             [](auto const &f) { return f == ForceStruct{}; });
}

} // namespace Testing

using namespace Testing;

BOOST_AUTO_TEST_CASE(links_of_one_boundary_cell) {
  Setup setup;
  // one boundary cell at the center of the first block, which has 6x5x4
  // cells: all its 18 neighbors are INNER cells except for the neighbors
  // in z, because z = 2 is next to the last cell layer
  setup.set_flags({Cell(2, 2, 1)});
  std::vector<Cell> cells;
  Boundary boundary(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "test");
  setup.fill(boundary);
  auto &block = *setup.blocks->begin();
  auto const &all = boundary.getIndexVector(&block);
  BOOST_REQUIRE_EQUAL(all.size(), 18u);
  // one link per direction of D3Q19 except for the rest direction, ordered
  // by the index of the direction, that points to the boundary cell
  std::size_t i = 0;
  for (auto d = Stencil::begin(); d != Stencil::end(); ++d) {
    if (stencil::cx[*d] == 0 and stencil::cy[*d] == 0 and
        stencil::cz[*d] == 0) {
      continue;
    }
    BOOST_CHECK_EQUAL(all[i].dir, int32_t(Stencil::idx[*d]));
    BOOST_CHECK_EQUAL(all[i].x + stencil::cx[*d], 2);
    BOOST_CHECK_EQUAL(all[i].y + stencil::cy[*d], 2);
    BOOST_CHECK_EQUAL(all[i].z + stencil::cz[*d], 1);
    BOOST_CHECK_EQUAL(all[i].value, 2. + 0.1 * 2. + 0.01 * 1.);
    // the value callback is called for the boundary cell, once per link
    BOOST_CHECK(cells[i] == Cell(2, 2, 1));
    ++i;
  }
  BOOST_CHECK_EQUAL(cells.size(), 18u);
  // the fluid cells with z = 0 are not INNER
  auto const outer =
      std::ranges::count_if(all, [](auto const &link) { return link.z == 0; });
  BOOST_CHECK_EQUAL(outer, 5);
  BOOST_CHECK_EQUAL(boundary.getIndexVector(&block, Links::OUTER).size(), 5u);
  BOOST_CHECK_EQUAL(boundary.getIndexVector(&block, Links::INNER).size(), 13u);
  // the other block has no boundary cell
  auto const &other = *(++setup.blocks->begin());
  BOOST_CHECK(boundary.getIndexVector(&other).empty());
}

BOOST_AUTO_TEST_CASE(index_vectors_from_flag_field) {
  for (auto const fraction : {0.02, 0.2, 0.5}) {
    Setup setup;
    setup.set_flags(fraction);
    std::vector<Cell> cells;
    Boundary boundary(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "test");
    setup.fill(boundary);
    std::size_t n_links = 0;
    for (auto &block : *setup.blocks) {
      for (auto type : {Links::ALL, Links::INNER, Links::OUTER}) {
        auto const expected = expected_links(setup, block, type);
        BOOST_CHECK(boundary.getIndexVector(&block, type) == expected);
        BOOST_CHECK_EQUAL(boundary.getForceVector(&block, type).size(),
                          expected.size());
      }
      n_links += boundary.getIndexVector(&block).size();
    }
    BOOST_CHECK_GT(n_links, 0u);
    BOOST_CHECK_EQUAL(cells.size(), n_links);
  }
}

BOOST_AUTO_TEST_CASE(fill_resets_the_forces) {
  Setup setup;
  setup.set_flags(0.3);
  std::vector<Cell> cells;
  Boundary boundary(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "test");
  setup.fill(boundary);
  for (auto &block : *setup.blocks) {
    boundary.run(&block);
  }
  auto &block = *setup.blocks->begin();
  auto *force_vector = block.getData<Forces>(boundary.getForceVectorID());
  force_vector->syncCPU();
  BOOST_REQUIRE(not is_zero(boundary.getForceVector(&block)));
  // fewer links: the forces of the old links are gone
  setup.set_flags(0.1);
  setup.fill(boundary);
  for (auto type : {Links::ALL, Links::INNER, Links::OUTER}) {
    auto const &forces = boundary.getForceVector(&block, type);
    BOOST_CHECK_EQUAL(forces.size(),
                      boundary.getIndexVector(&block, type).size());
    BOOST_CHECK(is_zero(forces));
  }
  // the device copies match the new vectors: the kernel doesn't fail
  boundary.run(&block);
  force_vector->syncCPU();
  BOOST_CHECK(not is_zero(boundary.getForceVector(&block)));
}

BOOST_AUTO_TEST_CASE(inner_and_outer_are_equivalent_to_run) {
  for (auto const fraction : {0.05, 0.3, 0.7}) {
    Setup setup;
    setup.set_flags(fraction);
    std::vector<Cell> cells;
    Boundary run_all(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "run");
    Boundary run_split(setup.blocks, setup.pdfs_id[1], Recorder{&cells},
                       "split");
    setup.fill(run_all);
    setup.fill(run_split);
    for (auto &block : *setup.blocks) {
      run_all.run(&block);
      run_split.inner(&block);
      run_split.outer(&block);
    }
    BOOST_CHECK(setup.pdfs(0) == setup.pdfs(1));
    BOOST_CHECK(setup.pdfs(0) != setup.initial_pdfs());
    for (auto &block : *setup.blocks) {
      auto *all = block.getData<Forces>(run_all.getForceVectorID());
      auto *split = block.getData<Forces>(run_split.getForceVectorID());
      all->syncCPU(Forces::ALL);
      split->syncCPU(Forces::INNER);
      split->syncCPU(Forces::OUTER);
      auto const &links = run_all.getIndexVector(&block);
      auto const &inner_links = run_split.getIndexVector(&block, Links::INNER);
      auto const &forces = run_all.getForceVector(&block);
      auto const &inner_forces = run_split.getForceVector(&block, Links::INNER);
      auto const &outer_forces = run_split.getForceVector(&block, Links::OUTER);
      BOOST_REQUIRE_EQUAL(inner_forces.size() + outer_forces.size(),
                          forces.size());
      // the force on each link is the same, whatever the sweep
      std::size_t i_inner = 0;
      std::size_t i_outer = 0;
      for (std::size_t i = 0; i < links.size(); ++i) {
        auto const is_inner =
            i_inner < inner_links.size() and links[i] == inner_links[i_inner];
        auto const &force =
            is_inner ? inner_forces[i_inner++] : outer_forces[i_outer++];
        BOOST_CHECK(force == forces[i]);
      }
      // the sweeps of a type don't change the forces of the other types
      BOOST_CHECK(is_zero(run_all.getForceVector(&block, Links::INNER)));
      BOOST_CHECK(is_zero(run_all.getForceVector(&block, Links::OUTER)));
      BOOST_CHECK(is_zero(run_split.getForceVector(&block)));
    }
  }
}

BOOST_AUTO_TEST_CASE(total_force) {
  Setup setup;
  setup.set_flags(0.3);
  std::vector<Cell> cells;
  Boundary boundary(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "test");
  setup.fill(boundary);
  for (auto &block : *setup.blocks) {
    // zero before the first sweep
    BOOST_CHECK(boundary.getForce(&block) == Vector3<double>(0.));
    boundary.run(&block);
    boundary.inner(&block);
    boundary.outer(&block);
    for (auto type : {Links::ALL, Links::INNER, Links::OUTER}) {
      auto *force_vector = block.getData<Forces>(boundary.getForceVectorID());
      force_vector->syncCPU(type);
      Vector3<double> expected(0.);
      for (auto const &force : boundary.getForceVector(&block, type)) {
        expected += Vector3<double>(force.F_0, force.F_1, force.F_2);
      }
      auto const force = boundary.getForce(&block, type);
      BOOST_CHECK_GT(force.length(), 0.);
      BOOST_CHECK(force == expected);
    }
  }
}

BOOST_AUTO_TEST_CASE(sweeps) {
  Setup setup;
  setup.set_flags(0.3);
  std::vector<Cell> cells;
  Boundary direct(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "direct");
  Boundary wrapped(setup.blocks, setup.pdfs_id[1], Recorder{&cells}, "wrapped");
  setup.fill(direct);
  setup.fill(wrapped);
  auto const sweep = wrapped.getSweep();
  auto const inner_sweep = wrapped.getInnerSweep();
  auto const outer_sweep = wrapped.getOuterSweep();
  for (auto &block : *setup.blocks) {
    direct(&block);
    direct.inner(&block);
    direct.outer(&block);
    sweep(&block);
    inner_sweep(&block);
    outer_sweep(&block);
  }
  BOOST_CHECK(setup.pdfs(0) == setup.pdfs(1));
}

BOOST_AUTO_TEST_CASE(block_without_links) {
  Setup setup;
  setup.set_flags(std::vector<Cell>{});
  std::vector<Cell> cells;
  Boundary boundary(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "test");
  setup.fill(boundary);
  BOOST_CHECK(cells.empty());
  for (auto &block : *setup.blocks) {
    auto *force_vector = block.getData<Forces>(boundary.getForceVectorID());
    // all vectors are empty, which doesn't dereference anything on the device
    force_vector->syncCPU();
    force_vector->syncCPU(Forces::INNER);
    force_vector->syncCPU(Forces::OUTER);
    force_vector->syncGPU();
    for (auto type : {Links::ALL, Links::INNER, Links::OUTER}) {
      BOOST_CHECK(boundary.getIndexVector(&block, type).empty());
      BOOST_CHECK(boundary.getForceVector(&block, type).empty());
      BOOST_CHECK(force_vector->empty(type));
      BOOST_CHECK(boundary.getForce(&block, type) == Vector3<double>(0.));
    }
    // there is nothing to do
    boundary.run(&block);
    boundary.inner(&block);
    boundary.outer(&block);
  }
  BOOST_CHECK(setup.pdfs(0) == setup.pdfs(1));
}

BOOST_AUTO_TEST_CASE(flags_that_do_not_exist) {
  Setup setup;
  setup.set_flags(0.3);
  std::vector<Cell> cells;
  Boundary boundary(setup.blocks, setup.pdfs_id[0], Recorder{&cells}, "test");
  setup.fill(boundary);
  auto &block = *setup.blocks->begin();
  auto const links = boundary.getIndexVector(&block);
  BOOST_REQUIRE(not links.empty());
  // the vectors stay as they are if a flag is not registered
  boundary.fillFromFlagField<Flags>(setup.blocks, setup.flags_id,
                                    FlagUID("unknown"), FlagUID("domain"));
  BOOST_CHECK(boundary.getIndexVector(&block) == links);
}

BOOST_AUTO_TEST_CASE(index_vectors_compare) {
  Links lhs;
  Links rhs;
  BOOST_CHECK(lhs == rhs);
  lhs.indexVector(Links::INNER).emplace_back(1, 2, 3, 4);
  BOOST_CHECK(lhs != rhs);
  rhs.indexVector(Links::INNER).emplace_back(1, 2, 3, 4);
  BOOST_CHECK(lhs == rhs);
  rhs.indexVector(Links::INNER).back().value = 1.;
  BOOST_CHECK(lhs != rhs);
  BOOST_CHECK(lhs.pointerCpu(Links::ALL) == nullptr);
  BOOST_CHECK(lhs.pointerCpu(Links::INNER) != nullptr);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  walberla::mpi_init();
  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

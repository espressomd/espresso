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
 * Tests of the generated flux kernels of the electrokinetics (the diffusive,
 * advective and continuity kernels, the friction coupling, the fixed flux
 * boundary and the accessors of the flux field) on the CPU. They check the
 * physics (analytic fluxes of linear fields, conservation of mass) and the
 * layout of the flux field, which has 13 components for each cell: the flux
 * through the face between the cell @c c and @c c + e along @c e, for the 13
 * directions @c e of @ref stored_directions. The fluxes of the other faces
 * of a cell are stored in its neighbors, with the opposite sign.
 */

#define BOOST_TEST_MODULE waLBerla EK flux kernels
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include "../src/electrokinetics/ek_kernels.hpp"

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Vector3.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/FlagUID.h>
#include <field/GhostLayerField.h>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace {

using namespace walberla;
using Dir = std::array<int, 3>;
using Size = Vector3<uint_t>;

/** Directions of the components of the flux field. */
constexpr std::array<Dir, 13> stored_directions = {{{0, 1, 0},
                                                    {1, 0, 0},
                                                    {0, 0, 1},
                                                    {1, 1, 0},
                                                    {1, -1, 0},
                                                    {0, 1, 1},
                                                    {1, 0, 1},
                                                    {0, 1, -1},
                                                    {1, 0, -1},
                                                    {1, 1, 1},
                                                    {-1, 1, 1},
                                                    {1, -1, 1},
                                                    {1, 1, -1}}};

/** Directions to the 26 neighbors of a cell. */
std::vector<Dir> all_directions() {
  std::vector<Dir> directions;
  for (int x = -1; x <= 1; ++x) {
    for (int y = -1; y <= 1; ++y) {
      for (int z = -1; z <= 1; ++z) {
        if (x != 0 or y != 0 or z != 0) {
          directions.push_back({x, y, z});
        }
      }
    }
  }
  return directions;
}

double length(Dir const &d) {
  return std::sqrt(double(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
}

double dot(std::array<double, 3> const &v, Dir const &d) {
  return v[0] * d[0] + v[1] * d[1] + v[2] * d[2];
}

/** Normalization of the fluxes of the flux kernels. */
double flux_normalization() {
  double sum = 0.;
  for (auto const &d : stored_directions) {
    sum += length(d);
  }
  return sum / 3.;
}

/** Flux through the face between @c c and @c c + d along @c d. */
template <typename Field>
double face_flux(Field const &flux, Cell const &c, Dir const &d) {
  for (std::size_t k = 0; k < stored_directions.size(); ++k) {
    if (stored_directions[k] == d) {
      return double(flux.get(c.x(), c.y(), c.z(), cell_idx_t(k)));
    }
  }
  Dir const inverse = {-d[0], -d[1], -d[2]};
  for (std::size_t k = 0; k < stored_directions.size(); ++k) {
    if (stored_directions[k] == inverse) {
      return -double(
          flux.get(c.x() + d[0], c.y() + d[1], c.z() + d[2], cell_idx_t(k)));
    }
  }
  throw std::logic_error("not a direction to a neighbor");
}

/** Set the flux through the face between @c c and @c c + e, with @c e stored.
 */
template <typename Field>
void set_face_flux(Field &flux, Cell const &c, std::size_t k, double value) {
  flux.get(c.x(), c.y(), c.z(), cell_idx_t(k)) =
      static_cast<typename Field::value_type>(value);
}

/** Blocks and fields of the electrokinetics kernels. */
template <typename T> struct Fixture {
  using Scalar = field::GhostLayerField<T, 1>;
  using Vector = field::GhostLayerField<T, 3>;
  using Flux = field::GhostLayerField<T, 13>;
  using Kernels = detail::KernelTrait<T, lbmpy::Arch::CPU>;

  std::shared_ptr<StructuredBlockForest> blocks;
  BlockDataID density_id, potential_id, velocity_id, force_id, flux_id;
  uint_t ghost_layers;
  Vector3<int64_t> grid_size;

  Fixture(Size const &n_blocks, Size const &n_cells, uint_t n_ghost_layers,
          T flux_value = T(0))
      : ghost_layers(n_ghost_layers),
        grid_size(int64_t(n_blocks[0] * n_cells[0]),
                  int64_t(n_blocks[1] * n_cells[1]),
                  int64_t(n_blocks[2] * n_cells[2])) {
    blocks = blockforest::createUniformBlockGrid(
        n_blocks[0], n_blocks[1], n_blocks[2], n_cells[0], n_cells[1],
        n_cells[2], 1., false, true, true, true);
    auto const add_scalar = [&](char const *name) {
      return field::addToStorage<Scalar>(blocks, name, T(0), field::fzyx,
                                         n_ghost_layers);
    };
    auto const add_vector = [&](char const *name) {
      return field::addToStorage<Vector>(blocks, name, T(0), field::fzyx,
                                         n_ghost_layers);
    };
    density_id = add_scalar("density");
    potential_id = add_scalar("potential");
    velocity_id = add_vector("velocity");
    force_id = add_vector("force");
    flux_id = field::addToStorage<Flux>(blocks, "flux", flux_value, field::fzyx,
                                        n_ghost_layers);
  }

  IBlock &block() { return *blocks->begin(); }

  /** Set a field in all cells, including the ghost layers, to a function
   *  of the global coordinates of the cells and of the component. */
  template <typename F, typename Function>
  void fill(BlockDataID id, Function &&function) {
    for (auto &block : *blocks) {
      auto const lower = blocks->getBlockCellBB(block).min();
      auto *field = block.template getData<F>(id);
      for (auto it = field->beginWithGhostLayerXYZ(); it != field->end();
           ++it) {
        for (cell_idx_t f = 0; f < cell_idx_t(field->fSize()); ++f) {
          auto const value =
              function(int64_t(lower[0] + it.x()), int64_t(lower[1] + it.y()),
                       int64_t(lower[2] + it.z()), int(f));
          field->get(it.x(), it.y(), it.z(), f) = static_cast<T>(value);
        }
      }
    }
  }

  /** Cells of a block that the flux kernels run on. */
  CellInterval flux_interval(IBlock &block) {
    auto interval = block.template getData<Scalar>(density_id)->xyzSize();
    interval.expand(1);
    return interval;
  }

  /** Global coordinates of the lower corner of a block. */
  Cell lower_corner(IBlock &block) {
    return blocks->getBlockCellBB(block).min();
  }

  /** Global coordinates of the first cell of the interval of a flux kernel
   *  (the RNG of a thermalized kernel depends on them). */
  Vector3<int64_t> block_offset(IBlock &block) {
    auto const lower = lower_corner(block) + flux_interval(block).min();
    return Vector3<int64_t>(lower[0], lower[1], lower[2]);
  }

  /** Cells of the interior of a block. */
  std::vector<Cell> interior(IBlock &block) {
    std::vector<Cell> cells;
    auto const interval = block.template getData<Scalar>(density_id)->xyzSize();
    for (auto const &cell : interval) {
      cells.push_back(cell);
    }
    return cells;
  }
};

/** Smallest tolerance for the comparison of single precision values. */
template <typename T> constexpr double tolerance() {
  return std::is_same_v<T, double> ? 1e-12 : 2e-5;
}

using Vec3 = std::array<double, 3>;

/** Check the fluxes through all the faces of the interior cells. */
template <typename T, typename Expected>
void check_all_faces(Fixture<T> &fixture, Expected &&expected) {
  auto &block = fixture.block();
  auto const *flux =
      block.template getData<typename Fixture<T>::Flux>(fixture.flux_id);
  double max_error = 0., max_value = 0.;
  for (auto const &cell : fixture.interior(block)) {
    for (auto const &d : all_directions()) {
      auto const reference = expected(cell, d);
      auto const value = face_flux(*flux, cell, d);
      max_error = std::max(max_error, std::abs(value - reference));
      max_value = std::max(max_value, std::abs(reference));
    }
  }
  BOOST_CHECK_GT(max_value, 0.);
  BOOST_CHECK_LE(max_error, tolerance<T>() * max_value);
}

template <typename T> void check_diffusive_flux_of_linear_density() {
  using F = Fixture<T>;
  F fixture(Size(1, 1, 1), Size(4, 5, 6), 2);
  Vec3 const slope = {0.3, -0.2, 0.5};
  T const diffusion = T(0.37);
  fixture.template fill<typename F::Scalar>(
      fixture.density_id, [&](auto x, auto y, auto z, auto) {
        return 2. + slope[0] * x + slope[1] * y + slope[2] * z;
      });
  typename F::Kernels::DiffusiveFluxKernel kernel(fixture.density_id,
                                                  fixture.flux_id, diffusion);
  auto &block = fixture.block();
  fixture.template fill<typename F::Flux>(
      fixture.flux_id, [](auto, auto, auto, auto) { return -777.; });
  kernel.runOnCellInterval(&block, fixture.flux_interval(block));
  auto const a0 = flux_normalization();
  check_all_faces(fixture, [&](Cell const &, Dir const &d) {
    // the gradient of a linear field is exact on all the faces
    return -double(diffusion) * dot(slope, d) / length(d) / a0;
  });
  // only the interior and the first ghost layer are computed
  auto const *flux = block.template getData<typename F::Flux>(fixture.flux_id);
  for (auto it = flux->beginWithGhostLayerXYZ(); it != flux->end(); ++it) {
    auto const outer = it.x() < -1 or it.y() < -1 or it.z() < -1 or
                       it.x() > 4 or it.y() > 5 or it.z() > 6;
    if (outer) {
      for (cell_idx_t k = 0; k < 13; ++k) {
        BOOST_REQUIRE_EQUAL(flux->get(it.x(), it.y(), it.z(), k), T(-777));
      }
    }
  }
}

template <typename T> void check_diffusive_flux_with_electrostatics() {
  using F = Fixture<T>;
  F fixture(Size(1, 1, 1), Size(4, 5, 6), 1);
  Vec3 const slope = {0.3, -0.2, 0.5};
  Vec3 const gradient = {0.1, 0.25, -0.15};
  Vec3 const ext_efield = {0.05, -0.1, 0.2};
  T const diffusion = T(0.37), kT = T(1.3), valency = T(2.);
  auto const density = [&](double x, double y, double z) {
    return 2. + slope[0] * x + slope[1] * y + slope[2] * z;
  };
  fixture.template fill<typename F::Scalar>(
      fixture.density_id,
      [&](auto x, auto y, auto z, auto) { return density(x, y, z); });
  fixture.template fill<typename F::Scalar>(
      fixture.potential_id, [&](auto x, auto y, auto z, auto) {
        return gradient[0] * x + gradient[1] * y + gradient[2] * z;
      });
  typename F::Kernels::DiffusiveFluxKernelElectrostatic kernel(
      fixture.density_id, fixture.flux_id, fixture.potential_id, diffusion,
      Vector3<T>(T(ext_efield[0]), T(ext_efield[1]), T(ext_efield[2])), kT,
      valency);
  auto &block = fixture.block();
  kernel.runOnCellInterval(&block, fixture.flux_interval(block));
  auto const a0 = flux_normalization();
  check_all_faces(fixture, [&](Cell const &c, Dir const &d) {
    // the density at the center of the face is the mean of the two cells
    auto const rho_face =
        density(c.x() + 0.5 * d[0], c.y() + 0.5 * d[1], c.z() + 0.5 * d[2]);
    Vec3 field;
    for (int i = 0; i < 3; ++i) {
      field[i] = gradient[i] - ext_efield[i];
    }
    auto const D = double(diffusion);
    return (-D * dot(slope, d) -
            D / double(kT) * rho_face * double(valency) * dot(field, d)) /
           length(d) / a0;
  });
}

BOOST_AUTO_TEST_CASE(diffusive_flux_of_a_linear_density) {
  check_diffusive_flux_of_linear_density<double>();
  check_diffusive_flux_of_linear_density<float>();
}

BOOST_AUTO_TEST_CASE(diffusive_flux_with_electrostatics) {
  check_diffusive_flux_with_electrostatics<double>();
  check_diffusive_flux_with_electrostatics<float>();
}

/** Random fluxes of the faces of a periodic block (and of its periodic
 *  images in the ghost layer). */
template <typename T> void randomize_periodic_flux(Fixture<T> &fixture) {
  std::mt19937_64 generator(42);
  std::uniform_real_distribution<double> distribution(-0.1, 0.1);
  auto &block = fixture.block();
  auto *flux =
      block.template getData<typename Fixture<T>::Flux>(fixture.flux_id);
  auto const size = fixture.grid_size;
  for (auto const &cell : fixture.interior(block)) {
    for (cell_idx_t k = 0; k < 13; ++k) {
      flux->get(cell.x(), cell.y(), cell.z(), k) = T(distribution(generator));
    }
  }
  auto const wrap = [](int64_t i, int64_t n) {
    return cell_idx_t(((i % n) + n) % n);
  };
  for (auto it = flux->beginWithGhostLayerXYZ(); it != flux->end(); ++it) {
    for (cell_idx_t k = 0; k < 13; ++k) {
      flux->get(it.x(), it.y(), it.z(), k) =
          flux->get(wrap(it.x(), size[0]), wrap(it.y(), size[1]),
                    wrap(it.z(), size[2]), k);
    }
  }
}

template <typename T> void check_continuity() {
  using F = Fixture<T>;
  F fixture(Size(1, 1, 1), Size(4, 5, 6), 1);
  fixture.template fill<typename F::Scalar>(
      fixture.density_id, [](auto x, auto y, auto z, auto) {
        return 1. + 0.01 * x + 0.02 * y + 0.03 * z;
      });
  randomize_periodic_flux(fixture);
  auto &block = fixture.block();
  auto *density =
      block.template getData<typename F::Scalar>(fixture.density_id);
  auto const *flux = block.template getData<typename F::Flux>(fixture.flux_id);
  // rho -= outflow - inflow
  std::map<std::tuple<int, int, int>, double> expected;
  double total_before = 0.;
  for (auto const &c : fixture.interior(block)) {
    double outflow = 0.;
    for (cell_idx_t k = 0; k < 13; ++k) {
      auto const &e = stored_directions[k];
      outflow += double(flux->get(c.x(), c.y(), c.z(), k)) -
                 double(flux->get(c.x() - e[0], c.y() - e[1], c.z() - e[2], k));
    }
    expected[{c.x(), c.y(), c.z()}] = double(density->get(c)) - outflow;
    total_before += double(density->get(c));
  }
  typename F::Kernels::ContinuityKernel kernel(fixture.density_id,
                                               fixture.flux_id);
  kernel(&block);
  double total_after = 0., max_error = 0.;
  for (auto const &c : fixture.interior(block)) {
    auto const value = double(density->get(c));
    total_after += value;
    max_error =
        std::max(max_error, std::abs(value - expected[{c.x(), c.y(), c.z()}]));
  }
  BOOST_CHECK_LE(max_error, 10. * tolerance<T>());
  // the fluxes of a periodic block conserve mass
  BOOST_CHECK_SMALL(total_after - total_before, 1e3 * tolerance<T>());
}

BOOST_AUTO_TEST_CASE(continuity_of_the_fluxes) {
  check_continuity<double>();
  check_continuity<float>();
}

template <typename T> void check_friction_coupling() {
  using F = Fixture<T>;
  F fixture(Size(1, 1, 1), Size(4, 5, 6), 1);
  randomize_periodic_flux(fixture);
  fixture.template fill<typename F::Vector>(
      fixture.force_id, [](auto x, auto y, auto z, auto f) {
        return 0.1 * x - 0.2 * y + 0.3 * z + f;
      });
  auto &block = fixture.block();
  auto *force = block.template getData<typename F::Vector>(fixture.force_id);
  auto const *flux = block.template getData<typename F::Flux>(fixture.flux_id);
  T const diffusion = T(0.37), kT = T(1.3), lb_density = T(0.8);
  std::vector<std::array<double, 3>> before;
  for (auto const &c : fixture.interior(block)) {
    before.push_back({double(force->get(c, 0)), double(force->get(c, 1)),
                      double(force->get(c, 2))});
  }
  typename F::Kernels::FrictionCouplingKernel kernel(
      fixture.flux_id, fixture.force_id, diffusion, kT, lb_density);
  kernel(&block);
  // f += kT / (2 D rho_lb) * sum of the fluxes through the 26 faces times
  // the direction, which is kT / (D rho_lb) times the flux vector
  auto const factor =
      double(kT) / (2. * double(diffusion) * double(lb_density));
  std::size_t index = 0;
  double max_error = 0., max_value = 0.;
  for (auto const &c : fixture.interior(block)) {
    for (cell_idx_t i = 0; i < 3; ++i) {
      double sum = 0.;
      for (auto const &d : all_directions()) {
        sum += face_flux(*flux, c, d) * d[i];
      }
      auto const reference = before[index][i] + factor * sum;
      max_error =
          std::max(max_error, std::abs(double(force->get(c, i)) - reference));
      max_value = std::max(max_value, std::abs(factor * sum));
    }
    ++index;
  }
  BOOST_CHECK_GT(max_value, 0.);
  BOOST_CHECK_LE(max_error, tolerance<T>() * 10. * max_value);
}

/** A uniform flux has no divergence: the density doesn't change at all. */
template <typename T> void check_continuity_of_a_uniform_flux() {
  using F = Fixture<T>;
  F fixture(Size(1, 1, 1), Size(4, 5, 6), 1);
  fixture.template fill<typename F::Scalar>(
      fixture.density_id, [](auto x, auto y, auto z, auto) {
        return 3e-4 + 1e-6 * x + 2e-6 * y + 3e-6 * z;
      });
  fixture.template fill<typename F::Flux>(
      fixture.flux_id,
      [](auto, auto, auto, auto k) { return 1.6524072e-7 * (1. + 0.37 * k); });
  auto &block = fixture.block();
  auto const *density =
      block.template getData<typename F::Scalar>(fixture.density_id);
  std::vector<T> before;
  for (auto const &c : fixture.interior(block)) {
    before.push_back(density->get(c));
  }
  typename F::Kernels::ContinuityKernel kernel(fixture.density_id,
                                               fixture.flux_id);
  for (int step = 0; step < 100; ++step) {
    kernel(&block);
  }
  std::size_t index = 0;
  for (auto const &c : fixture.interior(block)) {
    BOOST_REQUIRE_EQUAL(density->get(c), before[index++]);
  }
}

BOOST_AUTO_TEST_CASE(continuity_of_a_uniform_flux) {
  check_continuity_of_a_uniform_flux<double>();
  check_continuity_of_a_uniform_flux<float>();
}

BOOST_AUTO_TEST_CASE(friction_coupling_of_the_fluxes) {
  check_friction_coupling<double>();
  check_friction_coupling<float>();
}

template <typename T> void check_advective_flux() {
  using F = Fixture<T>;
  F fixture(Size(1, 1, 1), Size(5, 5, 5), 1);
  auto &block = fixture.block();
  auto const run = [&](Vec3 const &velocity, auto const &density) {
    fixture.template fill<typename F::Scalar>(
        fixture.density_id,
        [&](auto x, auto y, auto z, auto) { return density(x, y, z); });
    fixture.template fill<typename F::Vector>(
        fixture.velocity_id,
        [&](auto, auto, auto, auto f) { return velocity[f]; });
    fixture.template fill<typename F::Flux>(
        fixture.flux_id, [](auto, auto, auto, auto) { return 0.; });
    typename F::Kernels::AdvectiveFluxKernel advection(
        fixture.density_id, fixture.flux_id, fixture.velocity_id);
    advection.runOnCellInterval(&block, fixture.flux_interval(block));
    typename F::Kernels::ContinuityKernel continuity(fixture.density_id,
                                                     fixture.flux_id);
    continuity(&block);
  };
  auto const *density =
      block.template getData<typename F::Scalar>(fixture.density_id);
  auto const *flux = block.template getData<typename F::Flux>(fixture.flux_id);
  // a uniform density in a uniform flow is stationary, and the flux is not zero
  run({0.05, -0.03, 0.02}, [](auto, auto, auto) { return 1.5; });
  double max_flux = 0.;
  for (auto it = flux->beginWithGhostLayerXYZ(1); it != flux->end(); ++it) {
    for (cell_idx_t k = 0; k < 13; ++k) {
      max_flux = std::max(
          max_flux, std::abs(double(flux->get(it.x(), it.y(), it.z(), k))));
    }
  }
  BOOST_CHECK_GT(max_flux, 1e-3);
  for (auto const &c : fixture.interior(block)) {
    BOOST_CHECK_CLOSE(double(density->get(c)), 1.5, 1e4 * tolerance<T>());
  }
  // no flow, no change
  run({0., 0., 0.},
      [](auto x, auto y, auto z) { return 1. + 0.1 * x + 0.2 * y + 0.3 * z; });
  for (auto const &c : fixture.interior(block)) {
    BOOST_CHECK_CLOSE(double(density->get(c)),
                      1. + 0.1 * c.x() + 0.2 * c.y() + 0.3 * c.z(),
                      1e4 * tolerance<T>());
  }
  // a linear density profile along the flow moves with the velocity: the
  // flux through the faces along x is the density of the cell upstream times
  // the velocity
  run({0.1, 0., 0.}, [](auto x, auto, auto) { return 1. + 0.2 * x; });
  for (auto const &c : fixture.interior(block)) {
    BOOST_CHECK_CLOSE(double(density->get(c)), 1. + 0.2 * c.x() - 0.1 * 0.2,
                      1e4 * tolerance<T>());
  }
}

BOOST_AUTO_TEST_CASE(advective_flux_of_a_uniform_flow) {
  check_advective_flux<double>();
  check_advective_flux<float>();
}

/** Run a flux kernel with the RNG on all the blocks of a fixture. */
template <typename F, typename Kernel>
void run_thermalized_kernel(F &fixture, Kernel &kernel, int64_t time_step) {
  kernel.time_step() = time_step;
  for (auto &block : *fixture.blocks) {
    kernel.block_offset() = fixture.block_offset(block);
    kernel.runOnCellInterval(&block, fixture.flux_interval(block));
  }
}

/** Coordinate in the periodic grid. */
int64_t wrap(int64_t i, int64_t n) { return ((i % n) + n) % n; }

/**
 * Check that the fluxes of the thermalized kernel don't depend on the
 * decomposition of the domain: the faces of the interior cells of two blocks
 * have the same flux as in one block that has the same cells, in particular
 * the faces between a ghost cell and an interior cell. The density is a
 * function of the coordinates in the periodic grid, so that the ghost layers
 * don't need to be communicated.
 */
template <typename T, bool electrostatic> void check_decomposition() {
  using F = Fixture<T>;
  auto const grid = Vector3<int64_t>(8, 4, 4);
  auto const density = [&grid](int64_t x, int64_t y, int64_t z, auto) {
    return 1. + 0.01 * double(wrap(x, grid[0]) + 7 * wrap(y, grid[1]) +
                              3 * wrap(z, grid[2]));
  };
  auto const potential = [&grid](int64_t x, int64_t y, int64_t z, auto) {
    return 0.1 * double(wrap(x, grid[0]) * wrap(y, grid[1]) - wrap(z, grid[2]));
  };
  auto const run = [&](F &fixture, int64_t time_step) {
    fixture.template fill<typename F::Scalar>(fixture.density_id, density);
    fixture.template fill<typename F::Scalar>(fixture.potential_id, potential);
    auto const grid_size = fixture.grid_size;
    if constexpr (electrostatic) {
      typename F::Kernels::DiffusiveFluxKernelElectrostaticThermalized kernel(
          fixture.density_id, fixture.flux_id, fixture.potential_id,
          Vector3<int64_t>(0), T(0.3), Vector3<T>(T(0.1), T(-0.2), T(0.3)),
          grid_size, T(1.5), 42, time_step, T(2));
      run_thermalized_kernel(fixture, kernel, time_step);
    } else {
      typename F::Kernels::DiffusiveFluxKernelThermalized kernel(
          fixture.density_id, fixture.flux_id, Vector3<int64_t>(0), T(0.3),
          grid_size, 42, time_step);
      run_thermalized_kernel(fixture, kernel, time_step);
    }
  };
  // the fluxes of the faces in one block, by global cell and direction
  F one(Size(1, 1, 1), Size(8, 4, 4), 1);
  run(one, 5);
  std::map<std::tuple<int64_t, int64_t, int64_t, int, int, int>, double> fluxes;
  auto const &one_block = one.block();
  auto const *one_flux =
      one_block.template getData<typename F::Flux>(one.flux_id);
  for (auto const &cell : one.interior(one.block())) {
    for (auto const &d : all_directions()) {
      fluxes[{cell.x(), cell.y(), cell.z(), d[0], d[1], d[2]}] =
          face_flux(*one_flux, cell, d);
    }
  }
  F two(Size(2, 1, 1), Size(4, 4, 4), 1);
  run(two, 5);
  std::size_t n_compared = 0;
  for (auto &block : *two.blocks) {
    auto const offset = two.lower_corner(block);
    auto const *flux = block.template getData<typename F::Flux>(two.flux_id);
    for (auto const &cell : two.interior(block)) {
      for (auto const &d : all_directions()) {
        auto const global = cell + offset;
        auto const reference =
            fluxes.at({global.x(), global.y(), global.z(), d[0], d[1], d[2]});
        BOOST_REQUIRE_EQUAL(face_flux(*flux, cell, d), reference);
        ++n_compared;
      }
    }
  }
  BOOST_CHECK_EQUAL(n_compared, 8u * 4u * 4u * 26u);
  // the noise depends on the time step
  F other(Size(1, 1, 1), Size(8, 4, 4), 1);
  run(other, 6);
  auto const *other_flux =
      other.block().template getData<typename F::Flux>(other.flux_id);
  std::size_t n_equal = 0;
  for (auto const &cell : other.interior(other.block())) {
    for (auto const &d : all_directions()) {
      n_equal += face_flux(*other_flux, cell, d) ==
                 fluxes.at({cell.x(), cell.y(), cell.z(), d[0], d[1], d[2]});
    }
  }
  BOOST_CHECK_EQUAL(n_equal, 0u);
}

BOOST_AUTO_TEST_CASE(thermalized_flux_is_independent_of_the_decomposition) {
  check_decomposition<double, false>();
  check_decomposition<float, false>();
  check_decomposition<double, true>();
  check_decomposition<float, true>();
}

/**
 * Check the statistics of the noise: in a uniform density without
 * electrostatics, the flux through a face along the direction @c e has zero
 * mean and the variance <tt>2 D rho / |e| / A0</tt>, and the fluxes of the
 * different faces are not correlated.
 */
template <typename T> void check_noise_statistics() {
  using F = Fixture<T>;
  auto constexpr n = 16u;
  F fixture(Size(1, 1, 1), Size(n, n, n), 1);
  auto const density = T(1.5), diffusion = T(0.2);
  fixture.template fill<typename F::Scalar>(
      fixture.density_id, [&](auto, auto, auto, auto) { return density; });
  typename F::Kernels::DiffusiveFluxKernelThermalized kernel(
      fixture.density_id, fixture.flux_id, Vector3<int64_t>(0), diffusion,
      fixture.grid_size, 42, 0);
  run_thermalized_kernel(fixture, kernel, 3);
  auto const *flux =
      fixture.block().template getData<typename F::Flux>(fixture.flux_id);
  auto const cells = fixture.interior(fixture.block());
  auto const a0 = flux_normalization();
  std::vector<std::vector<double>> samples(13);
  for (std::size_t k = 0; k < 13; ++k) {
    for (auto const &c : cells) {
      samples[k].push_back(double(flux->get(c.x(), c.y(), c.z(), k)));
    }
  }
  auto const count = double(cells.size());
  for (std::size_t k = 0; k < 13; ++k) {
    auto const mean =
        std::accumulate(samples[k].begin(), samples[k].end(), 0.) / count;
    double variance = 0.;
    for (auto x : samples[k]) {
      variance += (x - mean) * (x - mean);
    }
    variance /= count;
    auto const expected = 2. * double(diffusion) * double(density) /
                          length(stored_directions[k]) / a0;
    // 5 standard deviations of the mean and 8% for the variance
    BOOST_CHECK_SMALL(mean, 5. * std::sqrt(expected / count));
    BOOST_CHECK_CLOSE(variance, expected, 8.);
  }
  for (std::size_t k = 1; k < 13; ++k) {
    double covariance = 0.;
    for (std::size_t i = 0; i < samples[0].size(); ++i) {
      covariance += samples[0][i] * samples[k][i];
    }
    auto const correlation =
        covariance / count /
        std::sqrt(2. * double(diffusion) * double(density) / a0 * 2. *
                  double(diffusion) * double(density) / a0 /
                  length(stored_directions[0]) / length(stored_directions[k]));
    BOOST_CHECK_SMALL(correlation, 0.08);
  }
}

BOOST_AUTO_TEST_CASE(statistics_of_the_thermal_noise) {
  check_noise_statistics<double>();
  check_noise_statistics<float>();
}

/**
 * Run the fixed flux boundary on a block with the boundary cells @c cells
 * and the values @c value of the callback, and check the flux field: every
 * link of a domain cell @c c in the interior with a boundary cell @c b =
 * c + d sets the flux through the face between them to
 * <tt>v(b) . d / 9</tt> along @c d, and nothing else is written.
 */
template <typename FixedFlux, typename T, typename Value>
void check_fixed_flux_sweep(std::vector<Cell> const &cells, Value &&value,
                            std::size_t expected_n_links) {
  using F = Fixture<T>;
  using Flags = field::FlagField<uint8_t>;
  F fixture(Size(1, 1, 1), Size(6, 5, 7), 1, T(-777));
  auto &block = fixture.block();
  auto const flags_id =
      field::addFlagFieldToStorage<Flags>(fixture.blocks, "flags", 1);
  auto *flags = block.template getData<Flags>(flags_id);
  auto const domain = FlagUID("domain");
  auto const boundary = FlagUID("boundary");
  flags->registerFlag(domain);
  flags->registerFlag(boundary);
  for (auto it = flags->beginWithGhostLayerXYZ(); it != flags->end(); ++it) {
    flags->addFlag(it.x(), it.y(), it.z(), flags->getFlag(domain));
  }
  for (auto const &cell : cells) {
    flags->addFlag(cell.x(), cell.y(), cell.z(), flags->getFlag(boundary));
  }
  std::function<Vector3<T>(Cell const &,
                           shared_ptr<StructuredBlockForest> const &, IBlock &)>
      callback = [&value](Cell const &cell,
                          shared_ptr<StructuredBlockForest> const &, IBlock &) {
        auto const v = value(cell);
        return Vector3<T>(T(v[0]), T(v[1]), T(v[2]));
      };
  FixedFlux fixed_flux(fixture.blocks, fixture.flux_id, callback);
  fixed_flux.template fillFromFlagField<Flags>(fixture.blocks, flags_id,
                                               boundary, domain);
  BOOST_CHECK_EQUAL(fixed_flux.getIndexVector(&block).size(), expected_n_links);
  fixed_flux(&block);
  // the expected flux field: the faces of the links
  std::map<std::tuple<int, int, int, std::size_t>, double> expected;
  auto const interval =
      block.template getData<typename F::Scalar>(fixture.density_id)->xyzSize();
  std::size_t n_links = 0;
  for (auto const &b : cells) {
    for (auto const &d : all_directions()) {
      Cell const c(b.x() - d[0], b.y() - d[1], b.z() - d[2]);
      if (not interval.contains(c)) {
        continue;
      }
      ++n_links;
      auto const flux = dot(value(b), d) / 9.;
      std::size_t k = 0;
      Cell owner = c;
      double sign = 1.;
      auto const stored =
          std::find(stored_directions.begin(), stored_directions.end(), d);
      if (stored != stored_directions.end()) {
        k = std::size_t(stored - stored_directions.begin());
      } else {
        Dir const inverse = {-d[0], -d[1], -d[2]};
        k = std::size_t(std::find(stored_directions.begin(),
                                  stored_directions.end(), inverse) -
                        stored_directions.begin());
        owner = b;
        sign = -1.;
      }
      auto const key = std::make_tuple(owner.x(), owner.y(), owner.z(), k);
      auto const found = expected.find(key);
      if (found != expected.end()) {
        // two links of boundary cells write the same face
        BOOST_REQUIRE_EQUAL(found->second, sign * flux);
      }
      expected[key] = sign * flux;
    }
  }
  BOOST_CHECK_EQUAL(n_links, expected_n_links);
  auto const *flux = block.template getData<typename F::Flux>(fixture.flux_id);
  for (auto it = flux->beginWithGhostLayerXYZ(); it != flux->end(); ++it) {
    for (std::size_t k = 0; k < 13; ++k) {
      auto const found = expected.find({it.x(), it.y(), it.z(), k});
      auto const value_in_field =
          double(flux->get(it.x(), it.y(), it.z(), cell_idx_t(k)));
      if (found == expected.end()) {
        BOOST_REQUIRE_EQUAL(value_in_field, -777.);
      } else {
        BOOST_REQUIRE_CLOSE(value_in_field, found->second, 1e-4);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(fixed_flux_sweep) {
  // isolated boundary cells, in the interior and in the ghost layer, with
  // different values; 26 links each, but the ones of the ghost cell
  // from the interior cells only (a link starts in the interior)
  std::vector<Cell> const cells = {Cell(1, 1, 1), Cell(4, 3, 5),
                                   Cell(-1, 2, 3)};
  auto const value = [](Cell const &cell) {
    return Vec3{0.5 + 0.25 * cell.x(), -1. + 0.125 * cell.y(),
                0.75 - 0.5 * cell.z()};
  };
  // the ghost cell (-1, 2, 3) has the link of the interior cell (0, 2, 3)
  // in the direction (-1, 0, 0) and the links of its 8 neighbors of x = 0
  auto const n_links = 26u + 26u + 9u;
  check_fixed_flux_sweep<typename Fixture<double>::Kernels::FixedFlux, double>(
      cells, value, n_links);
  check_fixed_flux_sweep<typename Fixture<float>::Kernels::FixedFlux, float>(
      cells, value, n_links);
}

BOOST_AUTO_TEST_CASE(fixed_flux_of_adjacent_boundary_cells) {
  // the face between two boundary cells is written by two links, which
  // store the same value if the value of the two cells is the same
  std::vector<Cell> const cells = {Cell(2, 2, 2), Cell(3, 2, 2), Cell(3, 3, 2)};
  auto const value = [](Cell const &) { return Vec3{0.25, -0.5, 0.75}; };
  // 3 cells with 26 links each, 2 of the links point to a cell that
  // is a boundary cell itself for the first and 4 for the second
  auto const n_links = 3u * 26u;
  check_fixed_flux_sweep<typename Fixture<double>::Kernels::FixedFlux, double>(
      cells, value, n_links);
}

BOOST_AUTO_TEST_CASE(flux_accessors) {
  using F = Fixture<double>;
  namespace accessor = ek::accessor;
  F fixture(Size(1, 1, 1), Size(3, 4, 5), 1);
  auto &block = fixture.block();
  auto *flux = block.template getData<F::Flux>(fixture.flux_id);
  // initialize and get
  std::array<double, 13> values;
  std::iota(values.begin(), values.end(), 0.5);
  accessor::Flux::initialize(flux, values);
  auto const interval = flux->xyzSize();
  for (auto it = flux->beginWithGhostLayerXYZ(); it != flux->end(); ++it) {
    for (cell_idx_t k = 0; k < 13; ++k) {
      BOOST_REQUIRE_EQUAL(flux->get(it.x(), it.y(), it.z(), k),
                          values[std::size_t(k)]);
    }
  }
  BOOST_CHECK(accessor::Flux::get(flux, Cell(1, 2, 3)) == values);
  // all the 13 components of every cell of an interval, in the order of the
  // cells
  auto const all = accessor::Flux::get(flux, interval);
  BOOST_REQUIRE_EQUAL(all.size(), 13u * interval.numCells());
  for (std::size_t i = 0; i < interval.numCells(); ++i) {
    for (std::size_t k = 0; k < 13; ++k) {
      BOOST_REQUIRE_EQUAL(all[13u * i + k], values[k]);
    }
  }
  // the flux vector of a uniform flux density is 9 times the density,
  // since the sum of d d^T over the 26 neighbors is 18 times the identity
  Vec3 const uniform = {0.1, -0.2, 0.3};
  for (auto it = flux->beginWithGhostLayerXYZ(); it != flux->end(); ++it) {
    for (cell_idx_t k = 0; k < 13; ++k) {
      flux->get(it.x(), it.y(), it.z(), k) =
          dot(uniform, stored_directions[std::size_t(k)]);
    }
  }
  auto const vector = accessor::Flux::get_vector(flux, Cell(1, 2, 3));
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(vector[i], 9. * uniform[i], 1e-9);
  }
  auto const vectors = accessor::Flux::get_vector(flux, interval);
  BOOST_REQUIRE_EQUAL(vectors.size(), 3u * interval.numCells());
  for (std::size_t i = 0; i < interval.numCells(); ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      BOOST_CHECK_CLOSE(vectors[3u * i + j], 9. * uniform[j], 1e-9);
    }
  }
  // a flux in the face between two cells is part of the flux vector of both
  accessor::Flux::initialize(flux, std::array<double, 13>{});
  set_face_flux(*flux, Cell(1, 2, 3), 1, 2.);
  BOOST_CHECK_CLOSE(accessor::Flux::get_vector(flux, Cell(1, 2, 3))[0], 1.,
                    1e-9);
  BOOST_CHECK_CLOSE(accessor::Flux::get_vector(flux, Cell(2, 2, 3))[0], 1.,
                    1e-9);
  BOOST_CHECK_SMALL(accessor::Flux::get_vector(flux, Cell(0, 2, 3))[0], 1e-12);
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  walberla::mpi_init();
  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

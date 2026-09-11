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

#define BOOST_TEST_MODULE HaloPlanValidator test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/AtomDecomposition.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/HybridDecomposition.hpp"
#include "cell_system/RegularDecomposition.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include <boost/mpi.hpp>

#include <functional>
#include <optional>
#include <set>
#include <span>
#include <vector>

BOOST_TEST_GLOBAL_FIXTURE(EspressoCoreGlobalConfig);

namespace {
BoxGeometry make_periodic_box(double box_l) {
  BoxGeometry box;
  box.set_length(Utils::Vector3d::broadcast(box_l));
  box.set_periodic(0u, true);
  box.set_periodic(1u, true);
  box.set_periodic(2u, true);
  return box;
}

RegularDecomposition make_dd(Utils::Vector3i const &node_grid, double box_l,
                             double range) {
  ::communicator.set_node_grid(node_grid);
  auto const box = make_periodic_box(box_l);
  auto const local_box = LocalBox::make_regular_decomposition(
      box.length(), ::communicator.calc_node_index(), ::communicator.node_grid);
  return {::communicator.comm, range, box, local_box, std::nullopt};
}

// The AtomDecomposition (n-square) covers every ghost via its collective
// broadcast/reduce section, so the validator must recognise collective-covered
// ghosts as covered. The box lifetime must outlive the returned decomposition;
// AtomDecomposition holds a BoxGeometry const& (see AtomDecomposition.hpp).
AtomDecomposition make_atom_dd(BoxGeometry const &box) {
  return {::communicator.comm, box};
}

// The HybridDecomposition combines the regular child's point-to-point section
// with the n-square child's collective section, so its plan exercises BOTH
// coverage paths at once.
HybridDecomposition make_hybrid_dd(Utils::Vector3i const &node_grid,
                                   BoxGeometry const &box,
                                   double cutoff_regular, double skin,
                                   std::set<int> const &n_square_types) {
  ::communicator.set_node_grid(node_grid);
  auto const local_box = LocalBox::make_regular_decomposition(
      box.length(), ::communicator.calc_node_index(), ::communicator.node_grid);
  return {
      ::communicator.comm, cutoff_regular, skin, []() { return false; }, box,
      local_box,           n_square_types};
}
} // namespace

namespace utf = boost::unit_test;
static bool has_4_mpi_ranks() { return boost::mpi::communicator{}.size() == 4; }
static bool has_1_mpi_rank() { return boost::mpi::communicator{}.size() == 1; }

// ========================
// single-rank local checks
// ========================
// These use hand-made cells and are rank-independent, so they pass at NUM_PROC
// 4 without modification.

BOOST_AUTO_TEST_CASE(detects_defects) {
  using namespace GhostComm;
  // 1 local cell whose only ghost neighbor is `ghost`.
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});
  std::vector<Cell *> locals{&local}, ghosts{&ghost};

  // Mark interior/boundary so the consistency check in validate_halo_plan
  // can correctly classify `local` as boundary (has a ghost neighbor).
  mark_boundary_cells(locals, ghosts);

  HaloPlan good;
  good.neighbors.push_back(NeighborComm{1, {SendRegion{&local, {}}}, {&ghost}});
  BOOST_CHECK(validate_halo_plan(good, locals, ghosts).empty());

  // uncovered ghost (neighborship-match + coverage failure)
  HaloPlan uncovered; // empty plan, ghost never filled
  BOOST_CHECK(not validate_halo_plan(uncovered, locals, ghosts).empty());

  // duplicate peer
  HaloPlan duppeer = good;
  duppeer.neighbors.push_back(
      NeighborComm{1, {SendRegion{&local, {}}}, {&ghost}});
  BOOST_CHECK(not validate_halo_plan(duppeer, locals, ghosts).empty());

  // double-filled ghost
  HaloPlan dbl = good;
  dbl.neighbors.push_back(NeighborComm{2, {SendRegion{&local, {}}}, {&ghost}});
  BOOST_CHECK(not validate_halo_plan(dbl, locals, ghosts).empty());

  // send/recv size mismatch
  HaloPlan badshape = good;
  badshape.neighbors[0].recv.clear();
  BOOST_CHECK(not validate_halo_plan(badshape, locals, ghosts).empty());

  // recv target outside ghost set
  Cell alien;
  HaloPlan alien_recv = good;
  alien_recv.neighbors[0].recv[0] = &alien;
  auto violations = validate_halo_plan(alien_recv, locals, ghosts);
  BOOST_CHECK_MESSAGE(not violations.empty(),
                      "Expected recv-outside-ghost violation");

  // isolated shape mismatch (send.size() != recv.size())
  HaloPlan shape_only = good;
  shape_only.neighbors[0].send.push_back(SendRegion{&local, {}});
  // send.size()==2, recv.size()==1 -> shape violation; ghost still covered once
  auto shape_violations = validate_halo_plan(shape_only, locals, ghosts);
  BOOST_CHECK_MESSAGE(not shape_violations.empty(),
                      "Expected shape-mismatch violation");
}

// ======================
// real-plan local checks
// ======================
// Build a real RegularDecomposition on 4 ranks and assert that the local
// validator (coverage + neighborship-match + peer-uniqueness + shape) passes.
// This is the guard that the wiring in RegularDecomposition.cpp won't
// false-trip in CI (which builds with ESPRESSO_ADDITIONAL_CHECKS on).
BOOST_AUTO_TEST_CASE(local_checks_pass_real_decomposition,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  auto const dd = make_dd({2, 2, 1}, 8., 1.2);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(
      violations.empty(),
      "Expected no local-check violations for real decomposition");
}

// Single-rank real RegularDecomposition. make_halo_plan() returns an empty
// plan on one rank (periodic neighbourships are wired directly between local
// cells in init_cell_interactions), yet mark_cells() still classifies the
// halo-layer cells as ghosts. Those ghosts are never referenced by a local
// cell's neighbor stencil, so the validator must NOT flag them as uncovered.
// This is the regression that aborted particle_reduction_test,
// p3m_test_1_mpi_ranks, EspressoSystemStandAlone_test_1_mpi_ranks etc. under a
// checks-on (RelWithAssert) build. Runs only when launched on a single rank.
BOOST_AUTO_TEST_CASE(local_checks_pass_single_rank_decomposition,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_1_mpi_rank();
                     })) {
  using namespace GhostComm;
  auto const dd = make_dd({1, 1, 1}, 8., 1.2);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  // Sanity: the single-rank plan is intentionally empty (no p2p, no local),
  // but there are ghost cells present.
  BOOST_CHECK(plan->neighbors.empty());
  BOOST_CHECK(plan->local.empty());
  BOOST_CHECK(not dd.ghost_cells().empty());
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no local-check violations for a single-rank "
                      "RegularDecomposition (unreferenced halo-layer ghosts)");
}

// =========================
// collective-covered ghosts
// =========================
// A real AtomDecomposition (n-square) covers ALL of its ghosts via the
// collective broadcast/reduce section, not via point-to-point recv/dst
// targets. Before the fix the validator false-positived on every such ghost
// ("never filled"). These cases run the validator directly in the checks-OFF
// build, so they catch the false-positives without needing a checks-on build.
BOOST_AUTO_TEST_CASE(atom_decomposition_collective_coverage,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  ::communicator.set_node_grid({4, 1, 1});
  auto const box = make_periodic_box(8.);
  auto const dd = make_atom_dd(box);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no local-check violations for real "
                      "AtomDecomposition (collective-covered ghosts)");
  BOOST_CHECK(validate_halo_plan_symmetry(*plan).empty());
}

// A real HybridDecomposition covers the regular ghosts via point-to-point and
// the n-square ghosts via the collective section — both paths at once.
BOOST_AUTO_TEST_CASE(hybrid_decomposition_mixed_coverage,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  auto const box = make_periodic_box(8.);
  auto const dd = make_hybrid_dd({2, 2, 1}, box, 1.2, 0.1, std::set<int>{0});
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations =
      validate_halo_plan(*plan, dd.local_cells(), dd.ghost_cells());
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no local-check violations for real "
                      "HybridDecomposition (mixed p2p + collective coverage)");
  BOOST_CHECK(validate_halo_plan_symmetry(*plan).empty());
}

// ============================================
// wrap-aware boundary classification on 1 rank
// ============================================
// On a single MPI rank all periodic axes are "wrap axes": no ghost cells are
// generated for them (init_cell_interactions wires local cells directly across
// the periodic boundary), so the base ghost-neighbour rule would leave those
// cells interior.  The wrap predicate must catch this and mark the
// first and last local layers along each periodic axis as boundary.
//
// Correctness argument: every cell in the first or last local layer (ghost-grid
// coord 1 or cell_grid[i] along axis i) has at least one neighbour in the
// opposite boundary layer.  After mark_boundary_cells with the wrap predicate
// those cells must be boundary.  Cells strictly in the interior (all ghost-grid
// coords in [2, cell_grid[i]-1] along every axis) have no wrap neighbour and
// must remain interior — PROVIDED the cell grid is large enough to have such
// interior cells (cell_grid[i] >= 3 on every axis).
BOOST_AUTO_TEST_CASE(single_rank_wrap_axis_cells_are_boundary,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_1_mpi_rank();
                     })) {
  using namespace GhostComm;
  // Use a box large enough for cell_grid >= 3 along every axis so that
  // genuinely interior cells exist (they would be missed by the wrap rule).
  // range=1.0 with box_l=6.0 gives cell_grid={6,6,6}.
  auto const dd = make_dd({1, 1, 1}, 6., 1.);

  // Every cell in the first or last local layer along any axis must be
  // boundary due to the periodic wrap.
  auto const &gcg = dd.ghost_cell_grid;
  auto const &cg = dd.cell_grid;

  bool any_wrap_layer_interior = false;

  for (Cell const *c : dd.local_cells()) {
    // Compute ghost-grid coordinate.
    int const idx = static_cast<int>(c - dd.cells.data());
    int const gx = idx % gcg[0];
    int const gy = (idx / gcg[0]) % gcg[1];
    int const gz = idx / (gcg[0] * gcg[1]);

    bool const in_wrap_layer = (gx == 1 || gx == cg[0] || gy == 1 ||
                                gy == cg[1] || gz == 1 || gz == cg[2]);

    if (in_wrap_layer and not c->is_boundary()) {
      any_wrap_layer_interior = true;
    }
  }

  BOOST_CHECK_MESSAGE(
      not any_wrap_layer_interior,
      "single-rank: all cells in first/last local layer must be boundary "
      "(periodic-wrap rule)");

  // Sanity: cells strictly in the interior must remain interior.
  if (cg[0] >= 3 && cg[1] >= 3 && cg[2] >= 3) {
    bool found_interior = false;
    for (Cell const *c : dd.local_cells()) {
      int const idx = static_cast<int>(c - dd.cells.data());
      int const gx = idx % gcg[0];
      int const gy = (idx / gcg[0]) % gcg[1];
      int const gz = idx / (gcg[0] * gcg[1]);
      bool const strictly_interior =
          (gx >= 2 && gx <= cg[0] - 1 && gy >= 2 && gy <= cg[1] - 1 &&
           gz >= 2 && gz <= cg[2] - 1);
      if (strictly_interior) {
        BOOST_CHECK_MESSAGE(not c->is_boundary(),
                            "strictly interior cell must not be boundary");
        found_interior = true;
      }
    }
    BOOST_CHECK_MESSAGE(found_interior,
                        "expected at least one strictly-interior cell "
                        "(cell_grid >= 3 on all axes)");
  }
}

// Degenerate wrap case: cell_grid==1 along at least one periodic axis.
//
// When node_grid[i]==1 && periodic[i] && cell_grid[i]==1 the single cell
// layer on that axis is its own periodic neighbour.  init_cell_interactions()
// drops the self-pair (ind1==ind2), so neighbors().all() has no entry along
// that axis and the wrap predicate in mark_boundary_cells() is never invoked.
// Without the post-predicate fixup every local cell would be left INTERIOR —
// wrong, because each cell interacts with itself across the periodic boundary.
//
// box_l=1.0, range=2.0 -> cell_grid[i]=floor(1/2)=0 -> clamped to 1 on all
// three axes.  Every local cell must be boundary.
BOOST_AUTO_TEST_CASE(single_rank_cell_grid_one_all_cells_are_boundary,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_1_mpi_rank();
                     })) {
  using namespace GhostComm;
  // range > box_l -> cell_grid clamps to {1,1,1}
  auto const dd = make_dd({1, 1, 1}, 1., 2.);
  // Verify the precondition: cell_grid should be {1,1,1}.
  BOOST_REQUIRE_EQUAL(dd.cell_grid[0], 1);
  BOOST_REQUIRE_EQUAL(dd.cell_grid[1], 1);
  BOOST_REQUIRE_EQUAL(dd.cell_grid[2], 1);
  // Every local cell must be boundary (self-interaction across periodic
  // boundary).
  for (Cell const *c : dd.local_cells()) {
    BOOST_CHECK_MESSAGE(c->is_boundary(),
                        "cell_grid==1 wrap axis: every local cell must be "
                        "boundary (self-interaction across periodic boundary)");
  }
  flush_runtime_errors_local(); // clear interaction range error messages
}

// ==========================
// cross-rank symmetry checks
// ==========================

// Good case: real RegularDecomposition on 4 ranks must be symmetric.
BOOST_AUTO_TEST_CASE(symmetry_good_real_decomposition,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  auto const dd = make_dd({2, 2, 1}, 8., 1.2);
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);
  auto violations = validate_halo_plan_symmetry(*plan);
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "Expected no symmetry violations for real decomposition");
}

// Asymmetric case: every rank posts a NeighborComm to peer=(me+1)%n with
// send=[1 cell] and recv=[] (expects 0 from peer), but peer sends 1 to me.
// Each rank should detect a mismatch for peer (me-1+n)%n.
BOOST_AUTO_TEST_CASE(symmetry_detects_mismatch,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator comm;
  int const n = comm.size();
  int const me = comm.rank();

  // Dummy cell for the send region.
  Cell dummy;

  HaloPlan asym;
  asym.comm = comm;
  // Every rank sends 1 item to (me+1)%n but declares recv=0 from that peer.
  // So rank (me+1)%n will receive my send count=1 via all_to_all but see that
  // my declared recv count from (me+1)%n is 0, which doesn't match what
  // (me+1)%n actually sends.  Each rank therefore detects a mismatch.
  int const send_peer = (me + 1) % n;
  asym.neighbors.push_back(
      NeighborComm{send_peer,
                   {SendRegion{&dummy, {}}}, // send.size() == 1
                   {}});                     // recv.size() == 0

  auto violations = validate_halo_plan_symmetry(asym);
  BOOST_CHECK_MESSAGE(not violations.empty(),
                      "Expected symmetry violation: recv=0 but peer sends 1");
}

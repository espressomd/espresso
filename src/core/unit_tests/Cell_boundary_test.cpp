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

#define BOOST_TEST_MODULE Cell_boundary test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "cell_system/Cell.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include <algorithm>
#include <string>
#include <vector>

// ===================
// mark_boundary_cells
// ===================

// A local cell whose only neighbors are other local cells must be interior.
BOOST_AUTO_TEST_CASE(local_only_neighbors_is_interior) {
  Cell local_a, local_b, local_c;
  // local_a has two neighbors, both local (not ghost)
  local_a.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_b, &local_c}, {});

  std::vector<Cell *> locals{&local_a, &local_b, &local_c};
  std::vector<Cell *> ghosts{};

  GhostComm::mark_boundary_cells(locals, ghosts);

  BOOST_CHECK_MESSAGE(!local_a.is_boundary(),
                      "cell with only local neighbors must be interior");
  BOOST_CHECK_MESSAGE(!local_b.is_boundary(),
                      "local_b (no neighbors listed) must be interior");
  BOOST_CHECK_MESSAGE(!local_c.is_boundary(),
                      "local_c (no neighbors listed) must be interior");
}

// A local cell that has at least one ghost neighbor must be boundary.
BOOST_AUTO_TEST_CASE(ghost_neighbor_makes_boundary) {
  Cell local_a, local_b, ghost;
  // local_a has one local neighbor and one ghost neighbor
  local_a.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_b, &ghost}, {});
  // local_b has only a local neighbor
  local_b.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_a}, {});

  std::vector<Cell *> locals{&local_a, &local_b};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);

  BOOST_CHECK_MESSAGE(local_a.is_boundary(),
                      "cell with a ghost neighbor must be boundary");
  BOOST_CHECK_MESSAGE(!local_b.is_boundary(),
                      "cell without ghost neighbors must be interior");
}

// Calling mark_boundary_cells a second time (rebuild) must yield the same
// result (reset + re-mark = idempotent).
BOOST_AUTO_TEST_CASE(mark_boundary_cells_is_idempotent) {
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});
  std::vector<Cell *> locals{&local};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);
  bool const first_result = local.is_boundary();

  GhostComm::mark_boundary_cells(locals, ghosts);
  bool const second_result = local.is_boundary();

  BOOST_CHECK_EQUAL(first_result, second_result);
  BOOST_CHECK(first_result); // should be boundary
}

// ==============
// wrap-predicate
// ==============

// A pair of local cells connected by a wrap-predicate that always fires must
// both be marked boundary, even without any ghost neighbors.
BOOST_AUTO_TEST_CASE(wrap_predicate_makes_boundary) {
  Cell local_a, local_b;
  // local_a and local_b are mutual local neighbors (no ghosts).
  local_a.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_b}, {});
  local_b.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_a}, {});

  std::vector<Cell *> locals{&local_a, &local_b};
  std::vector<Cell *> ghosts{};

  // Predicate that fires for every pair: simulates the wrap axis case.
  auto always_wrap = [](Cell const *, Cell const *) { return true; };

  GhostComm::mark_boundary_cells(locals, ghosts, always_wrap);

  BOOST_CHECK_MESSAGE(
      local_a.is_boundary(),
      "cell with a wrap-predicate neighbor must be boundary (no ghost needed)");
  BOOST_CHECK_MESSAGE(
      local_b.is_boundary(),
      "cell with a wrap-predicate neighbor must be boundary (no ghost needed)");
}

// A pair of fully-local cells with NO ghost neighbors and a wrap predicate
// that never fires must remain interior.
BOOST_AUTO_TEST_CASE(no_wrap_no_ghost_stays_interior) {
  Cell local_a, local_b;
  local_a.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_b}, {});
  local_b.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_a}, {});

  std::vector<Cell *> locals{&local_a, &local_b};
  std::vector<Cell *> ghosts{};

  // Predicate that never fires: no wrap axis.
  auto never_wrap = [](Cell const *, Cell const *) { return false; };

  GhostComm::mark_boundary_cells(locals, ghosts, never_wrap);

  BOOST_CHECK_MESSAGE(
      !local_a.is_boundary(),
      "cell with only local non-wrap neighbors must be interior");
  BOOST_CHECK_MESSAGE(
      !local_b.is_boundary(),
      "cell with only local non-wrap neighbors must be interior");
}

// ==================================
// validator overlap-safety invariant
// ==================================

// An interior cell that appears as a NeighborComm send source must trigger
// the overlap-safety invariant violation.
BOOST_AUTO_TEST_CASE(validator_fires_for_interior_cell_in_send_region) {
  using namespace GhostComm;
  // Two local cells: local_a is interior (no ghost neighbors),
  // local_b is boundary (ghost neighbor).
  Cell local_a, local_b, ghost;
  local_a.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&local_b}, {});
  local_b.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_a, &ghost}, {});

  std::vector<Cell *> locals{&local_a, &local_b};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);
  BOOST_REQUIRE(!local_a.is_boundary()); // pre-condition: local_a is interior
  BOOST_REQUIRE(local_b.is_boundary());  // pre-condition: local_b is boundary

  // Build a plan that covers the ghost (so coverage checks pass), but
  // incorrectly lists local_a (interior) as a send source.
  HaloPlan plan;
  plan.neighbors.push_back(
      NeighborComm{1,
                   {SendRegion{&local_a, {}}, // ← interior cell!
                    SendRegion{&local_b, {}}},
                   {&ghost, &ghost}});
  // Fix double-fill by using a separate ghost for the second recv slot.
  Cell ghost2;
  ghosts.push_back(&ghost2);
  plan.neighbors[0].recv = {&ghost, &ghost2};

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_overlap =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("overlap-safety invariant") != std::string::npos;
      });
  BOOST_CHECK_MESSAGE(found_overlap,
                      "expected overlap-safety invariant violation; got: " +
                          (violations.empty() ? "(none)" : violations.front()));
}

// An interior cell that appears as a LocalComm src must trigger the
// overlap-safety invariant violation.
BOOST_AUTO_TEST_CASE(validator_fires_for_interior_cell_in_local_comm_src) {
  using namespace GhostComm;
  Cell local_interior, local_boundary, ghost;
  local_interior.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_boundary}, {});
  local_boundary.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&local_interior, &ghost}, {});

  std::vector<Cell *> locals{&local_interior, &local_boundary};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);
  BOOST_REQUIRE(!local_interior.is_boundary());
  BOOST_REQUIRE(local_boundary.is_boundary());

  // Plan with a LocalComm that uses local_interior as src -> violation.
  HaloPlan plan;
  plan.local.push_back(LocalComm{&local_interior, &ghost, {}});

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_overlap =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("overlap-safety invariant") != std::string::npos;
      });
  BOOST_CHECK_MESSAGE(found_overlap,
                      "expected overlap-safety invariant violation; got: " +
                          (violations.empty() ? "(none)" : violations.front()));
}

// A correctly-marked boundary cell as send source must NOT trigger the
// overlap-safety invariant (sanity check that we don't over-fire).
BOOST_AUTO_TEST_CASE(validator_silent_for_boundary_cell_in_send_region) {
  using namespace GhostComm;
  Cell local_boundary, ghost;
  local_boundary.m_neighbors =
      Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});

  std::vector<Cell *> locals{&local_boundary};
  std::vector<Cell *> ghosts{&ghost};

  GhostComm::mark_boundary_cells(locals, ghosts);
  BOOST_REQUIRE(local_boundary.is_boundary());

  HaloPlan plan;
  plan.neighbors.push_back(
      NeighborComm{1, {SendRegion{&local_boundary, {}}}, {&ghost}});

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_overlap =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("overlap-safety invariant") != std::string::npos;
      });
  BOOST_CHECK_MESSAGE(!found_overlap,
                      "boundary cell as send source must not trigger overlap "
                      "invariant; violations: " +
                          (violations.empty() ? "(none)" : violations.front()));
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "expected no violations for a correct plan");
}

// ===========================
// validator consistency check
// ===========================

// A local cell whose is_boundary()==false but that has a ghost neighbor must
// trigger the "interior cell has a ghost neighbor" consistency violation.
BOOST_AUTO_TEST_CASE(validator_fires_for_mismarked_interior_cell) {
  using namespace GhostComm;
  Cell local, ghost;
  // Wire ghost as a neighbor of local, but do NOT call mark_boundary_cells,
  // so local.is_boundary() stays false (the default).
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});

  std::vector<Cell *> locals{&local};
  std::vector<Cell *> ghosts{&ghost};

  // Build a minimal valid plan that covers the ghost so that the coverage
  // checks pass — we only want to test the consistency check in isolation.
  HaloPlan plan;
  plan.neighbors.push_back(NeighborComm{1, {SendRegion{&local, {}}}, {&ghost}});

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_consistency =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("interior cell has a ghost neighbor") !=
               std::string::npos;
      });
  BOOST_CHECK_MESSAGE(
      found_consistency,
      "expected 'interior cell has a ghost neighbor' violation; "
      "got: " +
          (violations.empty() ? "(none)" : violations.front()));
}

// A correctly-marked boundary cell must NOT trigger the consistency violation.
BOOST_AUTO_TEST_CASE(validator_silent_for_correctly_marked_boundary_cell) {
  using namespace GhostComm;
  Cell local, ghost;
  local.m_neighbors = Neighbors<Cell *>(std::vector<Cell *>{&ghost}, {});

  std::vector<Cell *> locals{&local};
  std::vector<Cell *> ghosts{&ghost};

  // Properly mark the cell as boundary first.
  GhostComm::mark_boundary_cells(locals, ghosts);
  BOOST_REQUIRE(local.is_boundary()); // pre-condition for this test

  HaloPlan plan;
  plan.neighbors.push_back(NeighborComm{1, {SendRegion{&local, {}}}, {&ghost}});

  auto violations = validate_halo_plan(plan, locals, ghosts);
  bool found_consistency =
      std::ranges::any_of(violations, [](std::string const &v) {
        return v.find("interior cell has a ghost neighbor") !=
               std::string::npos;
      });
  BOOST_CHECK_MESSAGE(
      !found_consistency,
      "expected no consistency violation for a correctly-marked "
      "boundary cell; got: " +
          (violations.empty() ? "(none)" : violations.front()));
  BOOST_CHECK_MESSAGE(violations.empty(),
                      "expected no violations at all for a correct plan");
}

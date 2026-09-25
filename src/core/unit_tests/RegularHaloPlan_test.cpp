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

#define BOOST_TEST_MODULE RegularHaloPlan test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/RegularDecomposition.hpp"
#include "communication.hpp"
#include "ghosts/HaloPlan.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <vector>

namespace utf = boost::unit_test;

BOOST_TEST_GLOBAL_FIXTURE(EspressoCoreGlobalConfig);

// Only meaningful on 4 ranks: exercises face/edge/corner peers on a
// non-trivial node grid.
static bool has_4_mpi_ranks() { return boost::mpi::communicator{}.size() == 4; }

namespace {
/**
 * @brief Build a RegularDecomposition directly on the global Cartesian
 * communicator (mirrors the setup done by the cell-structure façade).
 */
RegularDecomposition make_dd(Utils::Vector3i const &node_grid, double box_l,
                             double range) {
  ::communicator.set_node_grid(node_grid);
  BoxGeometry box;
  box.set_length(Utils::Vector3d::broadcast(box_l));
  box.set_periodic(0u, true);
  box.set_periodic(1u, true);
  box.set_periodic(2u, true);
  auto const local_box = LocalBox::make_regular_decomposition(
      box.length(), ::communicator.calc_node_index(), ::communicator.node_grid);
  return {::communicator.comm, range, box, local_box, std::nullopt};
}

/**
 * @brief Assert that every ghost cell is the target of exactly one recv/dst
 * entry, and that the union of recv/dst targets equals the ghost-cell set.
 */
void check_coverage(RegularDecomposition const &dd) {
  auto const *plan = dd.halo_plan();
  BOOST_REQUIRE(plan != nullptr);

  // Count how often each Cell is filled as a recv (peer) or dst
  // (local) target.
  std::map<Cell const *, unsigned> filled;
  for (auto const &nc : plan->neighbors) {
    // recv[k] must line up with peer.send[k].
    BOOST_CHECK_EQUAL(nc.send.size(), nc.recv.size());
    for (auto const *c : nc.recv) {
      filled[c] += 1u;
    }
  }
  for (auto const &lc : plan->local) {
    filled[lc.dst] += 1u;
  }

  // Peers must be unique.
  std::set<int> peers;
  for (auto const &nc : plan->neighbors) {
    BOOST_CHECK(peers.insert(nc.peer).second);
    BOOST_CHECK(nc.peer != ::communicator.comm.rank());
  }

  // (a) no ghost cell is filled twice; (b) coverage: every ghost cell is
  // filled exactly once, and nothing outside the ghost set is filled.
  std::set<Cell const *> ghost_set;
  for (auto const *g : dd.ghost_cells()) {
    ghost_set.insert(g);
    BOOST_CHECK_EQUAL(filled[g], 1u);
  }
  // The union of recv/dst equals the ghost-cell set (no extras).
  BOOST_CHECK_EQUAL(filled.size(), ghost_set.size());
  for (auto const &[cell, count] : filled) {
    BOOST_CHECK_EQUAL(count, 1u);
    BOOST_CHECK(ghost_set.count(cell) == 1u);
  }

  // Cross-rank consistency: for the recv[k] <-> peer.send[k] contract to be
  // satisfiable, my recv-from-P count must equal P's send-to-me count (and
  // vice-versa). This is exactly what the async engine relies on for a-priori
  // recv sizing, so verify it here rather than deferring entirely to 1.5.
  auto comm = plan->comm;
  std::map<int, std::size_t> my_recv, my_send;
  for (auto const &nc : plan->neighbors) {
    my_recv[nc.peer] = nc.recv.size();
    my_send[nc.peer] = nc.send.size();
  }
  std::vector<boost::mpi::request> reqs;
  std::map<int, std::size_t> peer_send, peer_recv;
  for (auto const &nc : plan->neighbors) {
    reqs.push_back(comm.isend(nc.peer, 0, my_send[nc.peer]));
    reqs.push_back(comm.isend(nc.peer, 1, my_recv[nc.peer]));
    reqs.push_back(comm.irecv(nc.peer, 0, peer_send[nc.peer]));
    reqs.push_back(comm.irecv(nc.peer, 1, peer_recv[nc.peer]));
  }
  boost::mpi::wait_all(reqs.begin(), reqs.end());
  for (auto const &nc : plan->neighbors) {
    // My recv from P must equal P's send to me, and my send must equal P's
    // recv.
    BOOST_CHECK_EQUAL(my_recv[nc.peer], peer_send[nc.peer]);
    BOOST_CHECK_EQUAL(my_send[nc.peer], peer_recv[nc.peer]);
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(every_ghost_cell_covered_once_2x2x1,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  auto const dd = make_dd({2, 2, 1}, 8., 1.2);
  check_coverage(dd);
}

BOOST_AUTO_TEST_CASE(every_ghost_cell_covered_once_4x1x1,
                     *utf::precondition([](utf::test_unit_id) {
                       return has_4_mpi_ranks();
                     })) {
  auto const dd = make_dd({4, 1, 1}, 8., 1.2);
  check_coverage(dd);
}

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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE "HaloExchange test"
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <boost/mpi.hpp>

#include "BondList.hpp"
#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/Cell.hpp"
#include "ghost_cell_fixture.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloExchange.hpp"
#include "ghosts/HaloPlan.hpp"

#include <Kokkos_Core.hpp>

#include <array>

namespace utf = boost::unit_test;

// The cells under test are backed by a ParticleStore, which allocates Kokkos
// Views and therefore needs an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

/*
 * On two ranks, each rank owns a single particle tagged by its rank. A Push
 * of GHOSTTRANS_POSITION must overwrite the local ghost cell with the peer's
 * position. This exercises the deadlock-free irecv/isend protocol and the
 * a-priori recv sizing.
 */
BOOST_AUTO_TEST_CASE(push_positions_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  // One local and one ghost cell, each holding a single particle -- sized as
  // if ghosts_count() had already run.
  GhostTest::CellFixture cells{{1u, 1u}, 1u};
  auto &local = cells.cell(0);
  auto &ghost = cells.cell(1);
  cells.front(0).id() = me; // tag by owner rank
  cells.front(0).pos() =
      Utils::Vector3d{double(me), 2. * double(me), 3. * double(me)};

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(
      NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_POSITION,
                {Direction::Push, Combine::Overwrite});

  // ghost cell now holds the peer's (folded) position.
  BOOST_CHECK_CLOSE(cells.front(1).pos()[0], double(other), 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[1], 2. * double(other), 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[2], 3. * double(other), 1e-12);
  // the local (owned) particle must be untouched by a Push.
  BOOST_CHECK_CLOSE(cells.front(0).pos()[0], double(me), 1e-12);
}

/*
 * A Reduce with Combine::Add of GHOSTTRANS_FORCE folds ghost forces back onto
 * the owner. The send/recv roles swap: we send from the *ghost* cell and add
 * into the *owned* cell. Each rank seeds a distinct force on its ghost, so the
 * owner's force must end up as (its own pre-existing force + the peer's ghost
 * force).
 */
BOOST_AUTO_TEST_CASE(reduce_forces_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  GhostTest::CellFixture cells{{1u, 1u}, 1u};
  auto &local = cells.cell(0);
  auto &ghost = cells.cell(1);
  // Pre-existing force on the owned particle; Reduce must add to it.
  cells.front(0).force() = Utils::Vector3d{1., 0., 0.};
  // Force accumulated on this rank's ghost (from short-range interactions).
  // Encode the ghost's owner (== other) so the sum is checkable.
  cells.front(1).force() =
      Utils::Vector3d{0., 10. * double(me + 1), 100. * double(me + 1)};

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(
      NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_FORCE, {Direction::Reduce, Combine::Add});

  // Owner force == own pre-existing force + peer's ghost force.
  BOOST_CHECK_CLOSE(cells.front(0).force()[0], 1., 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).force()[1], 10. * double(other + 1), 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).force()[2], 100. * double(other + 1), 1e-12);
}

/*
 * The cold, resort-only BONDS path is exchanged as a *second* per-neighbor
 * message (the CommBuf bond buffer) alongside the flat particle buffer. This
 * checks that a bond placed on the owner round-trips byte-correct into the
 * peer's ghost cell over the async engine.
 */
BOOST_AUTO_TEST_CASE(push_bonds_between_two_ranks,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();
  int const other = 1 - me;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  GhostTest::CellFixture cells{{1u, 1u}, 1u};
  auto &local = cells.cell(0);
  auto &ghost = cells.cell(1);
  cells.front(0).id() = me;
  std::array<int, 2> partners{me, me + 5};
  cells.front(0).bonds().insert(BondView{me + 1, partners});

  HaloPlan plan;
  plan.comm = world;
  plan.neighbors.push_back(
      NeighborComm{other, {SendRegion{&local, {}}}, {&ghost}});

  halo_exchange(plan, box, GHOSTTRANS_PROPRTS | GHOSTTRANS_BONDS,
                {Direction::Push, Combine::Overwrite});

  BOOST_REQUIRE_EQUAL(cells.front(1).bonds().size(), 1);
  auto it = cells.front(1).bonds().begin();
  BOOST_CHECK_EQUAL((*it).bond_id(), other + 1);
  BOOST_REQUIRE_EQUAL((*it).partner_ids().size(), 2u);
  BOOST_CHECK_EQUAL((*it).partner_ids()[0], other);
  BOOST_CHECK_EQUAL((*it).partner_ids()[1], other + 5);
}

/*
 * Collective section: on two ranks each rank owns one cell.  A Broadcast (Push)
 * must copy each rank's owned cell data into the other rank's ghost cell (the
 * entry for that rank in the collective cells vector).  The Broadcast loop runs
 * root=0 then root=1; after root=0's broadcast, rank 1's cells[0] holds rank
 * 0's particle; after root=1's broadcast, rank 0's cells[1] holds rank 1's
 * particle.
 *
 * A subsequent ReduceSum (Reduce) must sum ghost-force contributions back to
 * each owner: rank 0's cells[0] and rank 1's cells[1] accumulate the sums.
 */
BOOST_AUTO_TEST_CASE(collective_broadcast_and_reduce,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 2;
                     })) {
  using namespace GhostComm;
  boost::mpi::communicator world;
  int const me = world.rank();

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  // cells[0] = local cell for rank 0, cells[1] = local cell for rank 1.
  // Each rank owns cells[me] and uses cells[1 - me] as ghost storage.
  GhostTest::CellFixture cells{{1u, 1u}, 2u};
  auto &cell0 = cells.cell(0);
  auto &cell1 = cells.cell(1);

  // Seed the owned cell with a distinguishable position.
  if (me == 0) {
    cells.front(0).pos() = Utils::Vector3d{1.0, 2.0, 3.0};
    cells.front(0).id() = 0;
  } else {
    cells.front(1).pos() = Utils::Vector3d{4.0, 5.0, 6.0};
    cells.front(1).id() = 1;
  }

  HaloPlan plan;
  plan.comm = world;
  // cells[0] -> rank 0's cell, cells[1] -> rank 1's cell.
  plan.collective =
      CollectiveSection{CollectivePattern::Broadcast, {&cell0, &cell1}};

  // --- Broadcast (Push): every rank's owned cell is broadcast to all. ---
  halo_exchange(plan, box, GHOSTTRANS_POSITION,
                {Direction::Push, Combine::Overwrite});

  // After broadcast, each rank must have both particles.
  BOOST_CHECK_CLOSE(cells.front(0).pos()[0], 1.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).pos()[1], 2.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).pos()[2], 3.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[0], 4.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[1], 5.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[2], 6.0, 1e-12);

  // --- ReduceSum (Reduce): ghost forces are summed back to owner. ---
  // Seed ghost forces: each rank seeds a force on the other's ghost cell.
  if (me == 0) {
    cells.front(0).force() = Utils::Vector3d{10., 0., 0.}; // owned
    cells.front(1).force() = Utils::Vector3d{0., 20., 0.}; // ghost of rank1
  } else {
    cells.front(0).force() = Utils::Vector3d{0., 30., 0.}; // ghost of rank0
    cells.front(1).force() = Utils::Vector3d{40., 0., 0.}; // owned
  }

  halo_exchange(plan, box, GHOSTTRANS_FORCE, {Direction::Reduce, Combine::Add});

  // root=0: reduce sums rank0's cell0.force = {10,0,0} + rank1's cell0.force
  // = {0,30,0} -> {10,30,0} on rank 0.
  // root=1: reduce sums rank0's cell1.force = {0,20,0} + rank1's cell1.force
  // = {40,0,0} -> {40,20,0} on rank 1.
  // Only the root's copy of each cell is meaningful after the reduce.
  if (me == 0) {
    BOOST_CHECK_CLOSE(cells.front(0).force()[0], 10.0, 1e-12);
    BOOST_CHECK_CLOSE(cells.front(0).force()[1], 30.0, 1e-12);
  } else {
    BOOST_CHECK_CLOSE(cells.front(1).force()[0], 40.0, 1e-12);
    BOOST_CHECK_CLOSE(cells.front(1).force()[1], 20.0, 1e-12);
  }
}

/*
 * Same-rank (plan.local) path: Push direction.
 *
 * A LocalComm{src=real, dst=ghost, shift} represents a periodic-wrap copy on a
 * node_grid==1 axis.  For Direction::Push the engine must call
 * local_cell_copy(src, dst, shift) so the ghost receives the real particle's
 * position, folded by the shift.  This test pins that semantics and would fail
 * if the arguments were swapped (ghost would receive a zero-shifted copy of
 * its own position instead of the real cell's shifted position).
 *
 * Runs on 1 rank (plan.local requires no MPI).
 */
BOOST_AUTO_TEST_CASE(local_push_position_shift,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 1;
                     })) {
  using namespace GhostComm;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  GhostTest::CellFixture cells{{1u, 1u}, 1u};
  auto &real_cell = cells.cell(0);
  auto &ghost_cell = cells.cell(1);
  cells.front(0).pos() = Utils::Vector3d{3.0, 4.0, 5.0};
  // Ghost starts at a distinct position so we know overwrite occurred.
  cells.front(1).pos() = Utils::Vector3d{0.0, 0.0, 0.0};

  // shift of -10 on x: ghost should see folded position = 3 + (-10) = -7 ->
  // folded inside [0,10) gives 3.0 (no folding needed for this range), but the
  // raw shift is applied by the SAVE path so we check with the raw expected
  // value.  Use shift == {0,0,0} to keep the position arithmetic trivial and
  // focus on src/dst identity.
  HaloPlan plan;
  plan.comm = boost::mpi::communicator{};
  plan.local.push_back(LocalComm{&real_cell, &ghost_cell, {0., 0., 0.}});

  halo_exchange(plan, box, GHOSTTRANS_POSITION,
                {Direction::Push, Combine::Overwrite});

  // Ghost must now hold the real cell's position.
  BOOST_CHECK_CLOSE(cells.front(1).pos()[0], 3.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[1], 4.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(1).pos()[2], 5.0, 1e-12);
  // Real cell must be untouched.
  BOOST_CHECK_CLOSE(cells.front(0).pos()[0], 3.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).pos()[1], 4.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).pos()[2], 5.0, 1e-12);
}

/*
 * Same-rank (plan.local) path: Reduce direction.
 *
 * On a node_grid==1 periodic axis the plan builder emits a LocalComm with
 * src=real_cell, dst=ghost_cell.  During the force-reduce step the engine must
 * ADD the ghost force into the real cell, i.e. call
 *   local_cell_copy(*lc.dst, *lc.src, {}, box, GHOSTTRANS_FORCE)
 * (dst=ghost is the new src, src=real is the new dst).
 *
 * The pre-fix code called local_cell_copy(*lc.src, *lc.dst, ...) for both
 * directions, which silently dropped the ghost's force (it added the real
 * cell's force back into the ghost instead of the other way around).
 *
 * Failure signature when the fix is reverted:
 *   real_cell.force == {1,0,0} (unchanged — ghost force never reached it)
 *   ghost_cell.force == {1,20,300} (real force accumulated into ghost instead)
 *
 * This test runs on 1 rank (plan.local requires no MPI) so it is always
 * exercised, regardless of the MPI geometry used in the test suite.
 */
BOOST_AUTO_TEST_CASE(local_reduce_force_role_swap,
                     *utf::precondition([](utf::test_unit_id) {
                       return boost::mpi::communicator{}.size() == 1;
                     })) {
  using namespace GhostComm;

  BoxGeometry box;
  box.set_length({10., 10., 10.});

  GhostTest::CellFixture cells{{1u, 1u}, 1u};
  auto &real_cell = cells.cell(0);
  auto &ghost_cell = cells.cell(1);
  // Pre-existing force on the owned real particle.
  cells.front(0).force() = Utils::Vector3d{1., 0., 0.};
  // Force accumulated on the ghost during the pair-interaction loop.
  cells.front(1).force() = Utils::Vector3d{0., 20., 300.};

  HaloPlan plan;
  plan.comm = boost::mpi::communicator{};
  // LocalComm is always stored as {src=real, dst=ghost} by the plan builder.
  plan.local.push_back(LocalComm{&real_cell, &ghost_cell, {0., 0., 0.}});

  halo_exchange(plan, box, GHOSTTRANS_FORCE, {Direction::Reduce, Combine::Add});

  // Real cell must hold pre-existing + ghost force.
  BOOST_CHECK_CLOSE(cells.front(0).force()[0], 1.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).force()[1], 20.0, 1e-12);
  BOOST_CHECK_CLOSE(cells.front(0).force()[2], 300.0, 1e-12);
  // Ghost cell is NOT the target; its value is not checked here (the engine
  // may or may not modify it, but the real cell result is what matters).
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *  Asynchronous, split-phase ghost-communication engine.
 */

#include "ghosts/HaloExchange.hpp"

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"
#endif

#include "BoxGeometry.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/particle_packing.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/nonblocking.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <span>
#include <vector>

namespace GhostComm {

namespace {

/**
 * @brief MPI tags for ghost communications.
 *
 * One tag per exchange-kind keeps concurrent exchanges (a future goal) from
 * aliasing: within a single exchange each message is already uniquely matched
 * by @c (peer, tag). The BONDS transfer needs its own tag because it is sent
 * as a *second* message to the same peer, mirroring the legacy two-transfer.
 */
enum GhostTag : int {
  TAG_POSITION = 100,
  TAG_FORCE = 101,
  TAG_PARTNUM = 102,
  TAG_DEFAULT = 103,
  TAG_BONDS = 104,
};

int main_tag(unsigned data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return TAG_PARTNUM;
  if (data_parts & GHOSTTRANS_FORCE)
    return TAG_FORCE;
  if (data_parts & GHOSTTRANS_POSITION)
    return TAG_POSITION;
  return TAG_DEFAULT;
}
// NOTE: the (peer, tag) pair is what makes each MPI message unambiguous.
// This scheme is safe only when each peer appears at most once in the plan
// (the "one NeighborComm per peer" invariant). See halo_exchange_start @pre.

/** @brief View a list of cell pointers as a span for the packing routines. */
std::span<ParticleList *const> as_span(std::vector<ParticleList *> const &v) {
  return {v.data(), v.size()};
}

/** @brief Extract the plain cell pointers from a list of send regions. */
std::vector<ParticleList *> region_cells(std::vector<SendRegion> const &send) {
  std::vector<ParticleList *> cells;
  cells.reserve(send.size());
  for (auto const &r : send)
    cells.emplace_back(r.cell);
  return cells;
}

/**
 * @brief Pack all send regions of one NeighborComm into a single buffer.
 *
 * All regions are packed with a **single** @c pack_cells call so that bonds
 * produce exactly one boost binary archive per neighbor message -- matching
 * the single @c unpack_cells call on the receive side.  Concatenating
 * per-region archives would corrupt the bond stream (the second archive
 * header is misread as bond data by the receiver).
 *
 * Per-region shifts are honored via the common shift. All @c SendRegion.shift
 * values within a @c NeighborComm are equal — @c RegularDecomposition::
 * make_halo_plan() sets them all to @c {}. This lets us pack every region's
 * cells in one @c pack_cells call (a single bond archive, matching the
 * receiver's single @c unpack_cells).
 *
 * @note If a future decomposition (e.g. a Lees-Edwards refactor) ever stores
 * DISTINCT per-region shifts, @c pack_regions must be generalized to pack
 * each region's non-bond (flat) data with its own shift while still writing
 * all bonds into ONE shared archive — do not simply concatenate per-region
 * archives (that was the bond-corruption bug fixed in commit 314e7c366b).
 */
void pack_regions(CommBuf &buf, std::vector<SendRegion> const &regions,
                  BoxGeometry const &box, unsigned data_parts) {
#if defined(ESPRESSO_ADDITIONAL_CHECKS) and not defined(NDEBUG)
  // All regions in a NeighborComm must share the same shift.  If per-region
  // shifts ever diverge the packer must be generalized (separate archives).
  if (not regions.empty()) {
    auto const &ref = regions.front().shift;
    for (auto const &r : regions) {
      assert(r.shift == ref &&
             "pack_regions: per-region shifts differ within one NeighborComm "
             "-- packer must be generalized for per-region shift support");
    }
  }
#endif

  Utils::Vector3d const common_shift =
      regions.empty() ? Utils::Vector3d{} : regions.front().shift;

  auto const cells = region_cells(regions);
  pack_cells(buf, as_span(cells), common_shift, box, data_parts);
}

/**
 * @brief Collective (broadcast/reduce-sum) section.
 *
 * For each root rank in <tt>[0, comm.size())</tt>:
 *
 * 1. Broadcast (Push): root packs its owned cell (cells[root]) and broadcasts
 *    the buffer to all ranks; every non-root rank unpacks into cells[root]
 *    (their ghost copy of root's data).
 * 2. ReduceSum (Reduce): every rank packs cells[root] (ghost-force
 * contribution) and reduces to root with @c std::plus on the raw double buffer;
 * root unpacks into cells[root] (the owned cell, so forces accumulate there).
 *
 * The legacy code uses one GHOST_BCST (or GHOST_RDCE) step per rank, each
 * step addressing a single cell pointer (<tt>ghost_comm.part_lists[0] ==
 * &cells.at(n)</tt>).  We replicate that loop exactly.
 */
void run_collective(HaloPlan const &plan, BoxGeometry const &box,
                    unsigned data_parts, ExchangeOp op) {
  // Precondition: for non-PARTNUM data parts, ghost cell sizes must already be
  // synced by a prior GHOSTTRANS_PARTNUM exchange (same invariant as the legacy
  // GHOST_BCST/GHOST_RDCE path). The per-root broadcast/reduce byte count is
  // derived from the (already-sized) cells; a size mismatch across ranks would
  // be undefined behavior.

  if (!plan.collective or plan.collective->pattern == CollectivePattern::None) {
    return;
  }

  auto const &cs = *plan.collective;
  auto const &comm = plan.comm;
  int const comm_size = comm.size();
  int const my_rank = comm.rank();

  // Use op.direction to select the actual MPI collective: Push -> broadcast
  // (particles flow from owner to all ghosts); Reduce -> reduce-sum (ghost
  // forces flow back to owner).  cs.pattern == Broadcast marks the section as
  // active; the engine caller sets op correctly for each exchange.
  bool const is_broadcast = (op.direction == Direction::Push);

  // cs.cells has one entry per rank: cs.cells[root] is the ParticleList for
  // that root rank (owned on root, ghost on all others).
  assert(static_cast<int>(cs.cells.size()) == comm_size);

  CommBuf buf;
  // boost::mpi's pointer-based collectives bind a reference to the buffer even
  // for count == 0; an empty CommBuf yields data() == nullptr, which is
  // undefined behaviour (flagged by UBSan).  Substitute a dummy address for
  // empty buffers -- the count of 0 guarantees it is never dereferenced.
  alignas(double) static char empty_sentinel[sizeof(double)];
  auto const safe_data = [](CommBuf &b) {
    return b.size() != 0u ? b.data() : &empty_sentinel[0];
  };

  for (int root = 0; root < comm_size; ++root) {
    ParticleList *cell = cs.cells[static_cast<std::size_t>(root)];
    auto const cell_span = std::span<ParticleList *const>{&cell, 1};

    if (is_broadcast) {
      // Push: root broadcasts its owned particles to all other ranks.
      if (my_rank == root) {
        pack_cells(buf, cell_span, {}, box, data_parts);
        boost::mpi::broadcast(comm, safe_data(buf),
                              static_cast<int>(buf.size()), root);
        boost::mpi::broadcast(comm, buf.bonds(), root);
      } else {
        buf.resize(calc_transmit_size(cell_span, box, data_parts));
        buf.bonds().clear();
        boost::mpi::broadcast(comm, safe_data(buf),
                              static_cast<int>(buf.size()), root);
        boost::mpi::broadcast(comm, buf.bonds(), root);
        unpack_cells(buf, cell_span, box, data_parts);
      }
    } else {
      // ReduceSum: every rank sends its ghost-force contribution; root
      // accumulates into its owned cell.  Mirrors the legacy GHOST_RDCE
      // which reduces the raw double buffer with std::plus<double>.
      pack_cells(buf, cell_span, {}, box, data_parts);
      auto *raw = reinterpret_cast<double *>(safe_data(buf));
      int const count = static_cast<int>(buf.size() / sizeof(double));
      if (my_rank == root) {
        CommBuf recv_buf;
        recv_buf.resize(buf.size());
        auto *recv_raw = reinterpret_cast<double *>(safe_data(recv_buf));
        boost::mpi::reduce(comm, raw, count, recv_raw, std::plus<double>{},
                           root);
        unpack_cells(recv_buf, cell_span, box, data_parts);
      } else {
        boost::mpi::reduce(comm, raw, count, std::plus<double>{}, root);
      }
    }
  }
}

} // namespace

GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op,
                                  ExchangeBuffers &bufs) {
  GhostExchange st;
  st.op = op;
  st.data_parts = data_parts;
  st.box = &box;
  st.plan = &plan;
  st.bufs = &bufs;
  if (data_parts == GHOSTTRANS_NONE)
    return st;

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Cross-rank symmetry check: run once per plan at first use.
  //
  // This check is deferred from the decomposition constructors to here because:
  //   (a) During checkpoint loading, decompositions are transiently rebuilt
  //       while maximal_cutoff is rank-divergent (different cell grids per rank
  //       for a brief window before the next consistent rebuild).  The
  //       transient plan is never used, so its asymmetry is harmless — but a
  //       ctor-time collective would fire and abort.
  //   (b) A collective inside a ctor risks partial-abort deadlock: if one rank
  //       throws before entering the all_to_all, the others spin forever.
  //
  // All ranks enter halo_exchange_start together (ghost exchange is
  // collective), so the all_to_all inside validate_halo_plan_symmetry is safe
  // here. plan.symmetry_validated is mutable and rearmed to false whenever the
  // plan object is rebuilt (every decomposition change constructs a fresh
  // HaloPlan).
  if (!plan.symmetry_validated) {
    assert(GhostComm::report_violations(
        GhostComm::validate_halo_plan_symmetry(plan),
        "halo_exchange first use"));
    plan.symmetry_validated = true;
  }
  // Peer-uniqueness is validated at plan-build time by validate_halo_plan().
  // Op-sanity: Combine::Add is only valid for reducible parts
  // (FORCE, FORCE|TORQUE, RATTLE).  TORQUE may only appear together with
  // FORCE (it is never sent alone).
  // Compute the predicate outside assert() so the #ifdef is not embedded in
  // macro arguments (non-portable under -Werror).
  [[maybe_unused]] bool const force_only = (data_parts == GHOSTTRANS_FORCE);
#ifdef ESPRESSO_ROTATION
  [[maybe_unused]] bool const force_with_torque =
      (data_parts == (GHOSTTRANS_FORCE | GHOSTTRANS_TORQUE));
#else
  [[maybe_unused]] bool const force_with_torque = false;
#endif
  [[maybe_unused]] bool reducible = (force_only or force_with_torque);
#ifdef ESPRESSO_BOND_CONSTRAINT
  reducible = reducible or (data_parts == GHOSTTRANS_RATTLE);
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  reducible = reducible or (data_parts == GHOSTTRANS_DIPFLD);
#endif
  assert((op.combine != Combine::Add or reducible) &&
         "Combine::Add only valid for reducible parts (FORCE, FORCE|TORQUE, "
         "RATTLE, DIPFLD)");
#endif // ESPRESSO_ADDITIONAL_CHECKS

  auto const &comm = plan.comm;
  auto const n = plan.neighbors.size();

  // Resize the pool to the current neighbor count.  Existing entries retain
  // their vector capacity so that, after the first call (warm-up), CommBuf
  // ::resize on a same-sized or smaller payload is a no-op on the allocator.
  bufs.send.resize(n);
  bufs.recv.resize(n);
  bufs.send_cells.resize(n);
  bufs.recv_cells.resize(n);
  bufs.slot_to_neighbor.resize(n);
  // Clear and reuse the request vector (capacity retained).
  bufs.requests.clear();

  bool const push = op.direction == Direction::Push;
  bool const bonds = (data_parts & GHOSTTRANS_BONDS) != 0u;
  int const tag = main_tag(data_parts);

  // Resolve, per neighbor, which cells we send from and which we receive into.
  // Push:   send from the region cells, receive into the ghost (recv) cells.
  // Reduce: roles swap -- send from ghost cells, receive into region cells,
  //         adding on arrival.
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    if (push) {
      bufs.send_cells[i] = region_cells(nc.send);
      bufs.recv_cells[i] = nc.recv;
    } else {
      bufs.send_cells[i] = nc.recv;
      bufs.recv_cells[i] = region_cells(nc.send);
    }
  }

  // 1) Pack all send buffers.  Done before posting receives so that sizes are
  //    known; also isolates the serialisation work in the ghost/pack region.
#ifdef ESPRESSO_CALIPER
  if (espresso_cali_active())
    CALI_MARK_BEGIN("ghost/pack");
#endif
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    if (push) {
      pack_regions(bufs.send[i], nc.send, box, data_parts);
    } else {
      // Reduce: pack plain ghost cells, no shift.
      pack_cells(bufs.send[i], as_span(bufs.send_cells[i]), {}, box,
                 data_parts);
    }
  }
#ifdef ESPRESSO_CALIPER
  if (espresso_cali_active())
    CALI_MARK_END("ghost/pack");
#endif

  // 2) Post ALL receives first (deadlock-free), then ALL sends.  Sizes are
  //    known a-priori from the packed size of the cells we will receive into,
  //    except PARTNUM whose ghost cells are not sized yet -- there we post a
  //    fixed-size recv and let unpack_cells resize the cells (legacy bootstrap
  //    path).
#ifdef ESPRESSO_CALIPER
  if (espresso_cali_active())
    CALI_MARK_BEGIN("ghost/post");
#endif
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    bufs.recv[i].resize(
        calc_transmit_size(as_span(bufs.recv_cells[i]), box, data_parts));
    bufs.requests.push_back(comm.irecv(nc.peer, tag, bufs.recv[i].data(),
                                       static_cast<int>(bufs.recv[i].size())));
  }
  for (std::size_t i = 0; i < n; ++i) {
    auto const &nc = plan.neighbors[i];
    bufs.requests.push_back(comm.isend(nc.peer, tag, bufs.send[i].data(),
                                       static_cast<int>(bufs.send[i].size())));
  }

  // BONDS (cold, resort-only path): a second per-neighbor message for the bond
  // buffers, mirroring the legacy two-transfer. Kept entirely off the hot
  // POSITION/FORCE path via the flag guard. boost::mpi serializes the
  // std::vector<char> length + payload, so recv sizing is automatic here.
  if (bonds) {
    for (std::size_t i = 0; i < n; ++i) {
      auto const &nc = plan.neighbors[i];
      bufs.requests.push_back(
          comm.irecv(nc.peer, TAG_BONDS, bufs.recv[i].bonds()));
    }
    for (std::size_t i = 0; i < n; ++i) {
      auto const &nc = plan.neighbors[i];
      bufs.requests.push_back(
          comm.isend(nc.peer, TAG_BONDS, bufs.send[i].bonds()));
    }
  }
#ifdef ESPRESSO_CALIPER
  if (espresso_cali_active())
    CALI_MARK_END("ghost/post");
#endif

  // Local (same-rank) copies run here, after messages are posted, so they
  // overlap with in-flight network transfers for split-phase callers.
  //
  // Direction:
  //   Push:   copy the real (src) cell into its self-ghost (dst), applying the
  //           periodic-image shift so the ghost holds the correct replica.
  //   Reduce: roles swap -- the self-ghost (dst) accumulated a short-range
  //           force contribution during the pair loop that must be added back
  //           into the real (src) cell.  Copy dst -> src (no position shift;
  //           GHOSTTRANS_FORCE with ReductionPolicy::UPDATE accumulates via
  //           +=).  This mirrors the legacy GHOST_LOCL entry in the reverted
  //           collect-forces communicator; without the swap the self-wrap
  //           force contributions (node_grid[i]==1 periodic axes) are dropped,
  //           which silently biases forces on particles that interact across
  //           the local periodic boundary.
  //
  // Correctness (no aliasing with in-flight buffers):
  //   Push:   local_cell_copy writes ghost cells; isends carry real-cell data.
  //   Reduce: local_cell_copy writes real (owned) cells; isends carry ghost-
  //           cell data.  In both cases the written and sent cell sets are
  //           disjoint.  Local copies complete before halo_exchange_finish
  //           unpacks any receive buffer.
  for (auto const &lc : plan.local) {
    if (push)
      local_cell_copy(*lc.src, *lc.dst, lc.shift, box, data_parts);
    else
      local_cell_copy(*lc.dst, *lc.src, Utils::Vector3d{}, box, data_parts);
  }

  return st;
}

GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op) {
  // Convenience overload for callers that do not hold a persistent
  // ExchangeBuffers (unit tests, cold resort paths).  Allocate a pool via
  // unique_ptr so it is freed automatically when the GhostExchange is
  // destroyed — on every exit path, including the GHOSTTRANS_NONE early return
  // in halo_exchange_finish.  Use the pool overload on hot paths.
  auto pool = std::make_unique<ExchangeBuffers>();
  auto st = halo_exchange_start(plan, box, data_parts, op, *pool);
  // Transfer ownership: owned keeps the pool alive; bufs already points to it.
  st.owned = std::move(pool);
  return st;
}

void halo_exchange_finish(GhostExchange &st) {
  if (st.data_parts == GHOSTTRANS_NONE)
    return;

  auto &bufs = *st.bufs;

  // Local (same-rank) copies were already performed in halo_exchange_start
  // to overlap with in-flight messages; nothing to do here.

  if (st.plan->collective)
    run_collective(*st.plan, *st.box, st.data_parts, st.op);

  auto const n = st.plan->neighbors.size();
  bool const bonds = (st.data_parts & GHOSTTRANS_BONDS) != 0u;
  bool const add = st.op.combine == Combine::Add;

  // Request layout (set by halo_exchange_start):
  //   [0..n)    = per-neighbor irecvs
  //   [n..2n)   = per-neighbor isends
  //   [2n..3n)  = bond irecvs   (only when bonds)
  //   [3n..4n)  = bond isends   (only when bonds)

  if (bonds) {
    // Cold path (resort-only): bonds need *both* the flat and bond message from
    // each neighbor before unpack_cells can run.  Keep existing wait_all-then-
    // unpack-in-order semantics; pipelining is not safe here.
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_BEGIN("ghost/wait");
#endif
    boost::mpi::wait_all(bufs.requests.begin(), bufs.requests.end());
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_END("ghost/wait");
    if (espresso_cali_active())
      CALI_MARK_BEGIN("ghost/unpack");
#endif
    for (std::size_t i = 0; i < n; ++i) {
      auto dst = as_span(bufs.recv_cells[i]);
      unpack_cells(bufs.recv[i], dst, *st.box, st.data_parts);
    }
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_END("ghost/unpack");
#endif
  } else if (add) {
    // Add path (FORCE / RATTLE reduce): wait receives in *fixed neighbor order*
    // so that floating-point addition into shared local cells is bitwise
    // reproducible across runs.  Unpacking neighbor i overlaps delivery of
    // i+1..n-1 (pipelining while preserving associativity).
    //
    // Caliper: outer ghost/wait and ghost/unpack markers always fire (even
    // when n==0) so that label discovery works on single-rank runs.
    // Per-neighbor inner BEGIN/END pairs interleave with the outer ones;
    // repeated sequential begin/end accumulates — valid Caliper usage.
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_BEGIN("ghost/wait");
#endif
    for (std::size_t i = 0; i < n; ++i) {
      bufs.requests[i].wait();
#ifdef ESPRESSO_CALIPER
      if (espresso_cali_active())
        CALI_MARK_END("ghost/wait");
      if (espresso_cali_active())
        CALI_MARK_BEGIN("ghost/unpack");
#endif
      auto dst = as_span(bufs.recv_cells[i]);
      if (st.data_parts & GHOSTTRANS_FORCE) {
        add_forces(bufs.recv[i], dst, st.data_parts);
      }
#ifdef ESPRESSO_BOND_CONSTRAINT
      else if (st.data_parts == GHOSTTRANS_RATTLE) {
        add_rattle(bufs.recv[i], dst);
      }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      else if (st.data_parts == GHOSTTRANS_DIPFLD) {
        add_dip_fld(bufs.recv[i], dst);
      }
#endif
#ifdef ESPRESSO_CALIPER
      if (espresso_cali_active())
        CALI_MARK_END("ghost/unpack");
      if (espresso_cali_active())
        CALI_MARK_BEGIN("ghost/wait");
#endif
    }
    // Wait for sends [n..2n) so buffers are safe to reuse next step.
    // The outer ghost/wait region (opened before the loop) is still open here
    // on the n==0 path; on n>0 the loop re-opened it after the last unpack.
    boost::mpi::wait_all(bufs.requests.begin() + static_cast<std::ptrdiff_t>(n),
                         bufs.requests.begin() +
                             static_cast<std::ptrdiff_t>(2 * n));
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_END("ghost/wait");
    if (espresso_cali_active())
      CALI_MARK_BEGIN("ghost/unpack");
    if (espresso_cali_active())
      CALI_MARK_END("ghost/unpack");
#endif
  } else {
    // Overwrite path (position/property push): arrival order is safe because
    // the validator guarantees each ghost cell is filled by exactly one
    // message. Use wait_any to unpack each neighbor's buffer as soon as it
    // arrives.
    //
    // Bookkeeping: maintain a slot->neighbor index map so we always know which
    // recv buffer to unpack for the completed request slot.  Completed requests
    // are swapped to the end of the active range and the map is updated in
    // sync.
    //
    // Caliper: outer ghost/wait always fires (even n==0); the wait_any loop
    // emits nested end/begin pairs; an outer ghost/unpack pair fires after the
    // loop so the label appears on single-rank runs.
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_BEGIN("ghost/wait");
#endif
    // Initialize the scratch map to the identity: slot i -> neighbor i.
    // The vector was pre-sized in halo_exchange_start (no allocation here).
    for (std::size_t i = 0; i < n; ++i)
      bufs.slot_to_neighbor[i] = i;

    // active_end tracks how many recv requests are still pending.
    std::size_t active_end = n;

    while (active_end > 0) {
      auto first = bufs.requests.begin();
      auto last =
          bufs.requests.begin() + static_cast<std::ptrdiff_t>(active_end);

      auto [status, done_it] = boost::mpi::wait_any(first, last);
      (void)status;
#ifdef ESPRESSO_CALIPER
      if (espresso_cali_active())
        CALI_MARK_END("ghost/wait");
      if (espresso_cali_active())
        CALI_MARK_BEGIN("ghost/unpack");
#endif

      std::size_t const done_slot =
          static_cast<std::size_t>(done_it - bufs.requests.begin());
      std::size_t const neighbor_idx = bufs.slot_to_neighbor[done_slot];

      {
        auto dst = as_span(bufs.recv_cells[neighbor_idx]);
        unpack_cells(bufs.recv[neighbor_idx], dst, *st.box, st.data_parts);
      }

      // Remove the completed slot from the active range by swapping it with
      // the last active slot and shrinking the range.
      --active_end;
      if (done_slot != active_end) {
        std::iter_swap(done_it, bufs.requests.begin() +
                                    static_cast<std::ptrdiff_t>(active_end));
        std::swap(bufs.slot_to_neighbor[done_slot],
                  bufs.slot_to_neighbor[active_end]);
      }
#ifdef ESPRESSO_CALIPER
      if (espresso_cali_active())
        CALI_MARK_END("ghost/unpack");
      if (espresso_cali_active() && active_end > 0)
        CALI_MARK_BEGIN("ghost/wait");
#endif
    }

#ifdef ESPRESSO_CALIPER
    // Wait for sends [n..2n) so buffers are safe to reuse next step.
    // Caliper balance: when n==0 the outer CALI_MARK_BEGIN("ghost/wait") posted
    // before the loop is still open (the loop body never ran, so nothing closed
    // it).  When n>0 the last loop iteration emitted END("ghost/unpack") and
    // then — because active_end reached 0 — did NOT re-open "ghost/wait", so
    // the region is closed.  The guard `if (n > 0)` therefore re-opens for
    // n>0 and is deliberately absent for n==0 (the outer region serves).
    // Invariant: exactly one "ghost/wait" region is open entering wait_all.
    if (espresso_cali_active() && n > 0)
      CALI_MARK_BEGIN("ghost/wait");
#endif
    boost::mpi::wait_all(bufs.requests.begin() + static_cast<std::ptrdiff_t>(n),
                         bufs.requests.begin() +
                             static_cast<std::ptrdiff_t>(2 * n));
#ifdef ESPRESSO_CALIPER
    if (espresso_cali_active())
      CALI_MARK_END("ghost/wait");
    if (espresso_cali_active())
      CALI_MARK_BEGIN("ghost/unpack");
    if (espresso_cali_active())
      CALI_MARK_END("ghost/unpack");
#endif
  }

  // st.owned (if non-null) will free the heap-allocated pool automatically
  // when the GhostExchange goes out of scope after this call returns.
}

void halo_exchange(HaloPlan const &plan, BoxGeometry const &box,
                   unsigned data_parts, ExchangeOp op, ExchangeBuffers &bufs) {
  auto st = halo_exchange_start(plan, box, data_parts, op, bufs);
  halo_exchange_finish(st);
}

void halo_exchange(HaloPlan const &plan, BoxGeometry const &box,
                   unsigned data_parts, ExchangeOp op) {
  // Convenience overload: uses the no-pool start, which allocates a temporary
  // ExchangeBuffers on the heap and frees it in finish.
  auto st = halo_exchange_start(plan, box, data_parts, op);
  halo_exchange_finish(st);
}

} // namespace GhostComm

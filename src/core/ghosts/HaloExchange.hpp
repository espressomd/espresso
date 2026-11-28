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
#pragma once

/**
 * @file
 * Asynchronous, split-phase ghost-communication engine.
 *
 * Drive a @ref GhostComm::HaloPlan with non-blocking point-to-point MPI
 * (plus same-rank copies), reusing the byte-identical serialization from
 * @ref particle_packing.hpp. The engine is split into a
 * @c halo_exchange_start phase (post all receives, pack and post all sends)
 * and a @c halo_exchange_finish phase (same-rank copies, wait, unpack/reduce)
 * so that callers can overlap other work with in-flight messages.
 *
 * Deadlock-freedom is guaranteed by posting every @c irecv before any
 * @c isend and matching each message by @c (peer, tag); message sizes are
 * known a-priori from the receiving cells' packed size (see
 * @ref GhostComm::calc_transmit_size), except for the PARTNUM bootstrap.
 */

#include "BoxGeometry.hpp"
#include "ghosts/HaloPlan.hpp"
#include "ghosts/particle_packing.hpp"

#include <boost/mpi/request.hpp>

#include <memory>
#include <vector>

namespace GhostComm {

/**
 * @brief Persistent per-neighbor buffer pool for halo exchanges.
 *
 * Owns the heap-allocated send/recv buffers, MPI-request slots, and
 * cell-pointer scratch arrays for one set of neighbor communications.
 * Holding this object across multiple calls to @c halo_exchange allows
 * the underlying @c std::vector storage to be reused via @c resize —
 * after the first exchange (warm-up) the vectors retain capacity and
 * subsequent calls incur no heap allocations on the hot path.
 *
 * Typical ownership: held as a member of the object that drives ghost
 * communication (e.g. @c CellStructure::m_ghost_buffers) and passed to
 * @c halo_exchange / @c halo_exchange_start as a mutable ref.
 *
 * Thread-safety: a single instance must not be used from concurrent threads.
 * The sequential ghost-exchange loop (each exchange completes via
 * @c halo_exchange_finish before the next starts) is therefore safe.
 */
struct ExchangeBuffers {
  /** Per-neighbor packed send buffers (index-aligned with plan->neighbors). */
  std::vector<CommBuf> send;
  /** Per-neighbor recv buffers (index-aligned with plan->neighbors). */
  std::vector<CommBuf> recv;
  /** Outstanding non-blocking send/recv requests (cleared before each use). */
  std::vector<boost::mpi::request> requests;
  /**
   * Scratch cell-pointer arrays for the packing routines.  For the Reduce
   * direction send/recv roles swap; we need the plain @c ParticleList* from
   * the @c SendRegion list.  These must outlive the pack/unpack calls.
   */
  std::vector<std::vector<ParticleList *>> send_cells;
  std::vector<std::vector<ParticleList *>> recv_cells;
  /**
   * Scratch index map for the Overwrite (wait_any) path in
   * @c halo_exchange_finish: maps active request slot -> original neighbor
   * index.  Sized to the current neighbor count in @c halo_exchange_start
   * alongside the other per-neighbor vectors; capacity is retained across
   * calls so that no heap allocation occurs on the hot position-push path
   * after the first (warm-up) exchange.
   */
  std::vector<std::size_t> slot_to_neighbor;
};

/**
 * @brief Opaque handle for one in-flight halo exchange.
 *
 * Records the per-call metadata (op, data parts, geometry, plan) and holds a
 * pointer to the @ref ExchangeBuffers that back the in-flight messages.  When
 * a persistent @ref ExchangeBuffers is supplied by the caller (via the pool
 * overload), @c owned is null and @c bufs points to the caller's pool (no
 * ownership).  When the no-pool convenience overload is used, @c owned holds
 * the heap-allocated pool and @c bufs == @c owned.get(); the pool is freed
 * automatically when the handle is destroyed — on every exit path, including
 * the @c GHOSTTRANS_NONE early return in @c halo_exchange_finish.
 *
 * Created by @c halo_exchange_start and consumed exactly once by
 * @c halo_exchange_finish.
 */
struct GhostExchange {
  ExchangeOp op{};
  unsigned data_parts = 0u;
  BoxGeometry const *box = nullptr;
  HaloPlan const *plan = nullptr;
  /**
   * Non-null only when the no-pool overload allocated the buffer pool.
   * Destroyed (and the pool freed) when the @c GhostExchange goes out of
   * scope, guaranteeing cleanup on every exit path.
   */
  std::unique_ptr<ExchangeBuffers> owned;
  /** Non-owning pointer to the active buffer pool (caller's or @c owned). */
  ExchangeBuffers *bufs = nullptr;
};

/**
 * @brief Begin a halo exchange using a caller-owned buffer pool.
 *
 * The buffers in @p bufs are resized to the current neighbor count (retaining
 * capacity) so that, after the first call (warm-up), no heap allocation
 * occurs on the POSITION / FORCE hot path.
 *
 * @param plan        Communication plan (peers, regions, local copies).
 * @param box         Box geometry for position folding.
 * @param data_parts  Bitmask of GHOSTTRANS_* flags to transfer.
 * @param op          Direction (Push/Reduce) and combine mode (Overwrite/Add).
 * @param bufs        Persistent buffer pool (must outlive the returned handle).
 *
 * @pre Each peer rank appears at most once in @p plan.neighbors. Multiple
 *      send/receive regions to the same peer **must** be folded into that
 *      peer's single @ref NeighborComm (as the plan builders do). This is
 *      what makes the @c (peer, data-part tag) message matching unambiguous:
 *      without this invariant two @c NeighborComm entries sharing the same
 *      peer and tag would cross-match, silently corrupting data or deadlocking.
 *      The plan builders (e.g. RegularDecomposition) always produce unique
 *      peers; the invariant is asserted at runtime when
 *      @c ESPRESSO_ADDITIONAL_CHECKS is defined.
 */
GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op,
                                  ExchangeBuffers &bufs);

/**
 * @brief Convenience overload that allocates a temporary @ref ExchangeBuffers.
 *
 * Use this form when buffer reuse across calls is not needed (e.g. unit tests,
 * cold resort paths).  Equivalent to the original single-call allocation.
 */
GhostExchange halo_exchange_start(HaloPlan const &plan, BoxGeometry const &box,
                                  unsigned data_parts, ExchangeOp op);

/**
 * @brief Complete a halo exchange: run same-rank copies (overlapping the
 *        in-flight messages), wait for all requests, then unpack or reduce.
 */
void halo_exchange_finish(GhostExchange &state);

/**
 * @brief Blocking wrapper using a caller-owned buffer pool (no per-call alloc
 *        after warm-up).
 */
void halo_exchange(HaloPlan const &plan, BoxGeometry const &box,
                   unsigned data_parts, ExchangeOp op, ExchangeBuffers &bufs);

/**
 * @brief Blocking convenience wrapper: start + finish (allocates a temporary
 *        @ref ExchangeBuffers; use the pool overload on hot paths).
 */
void halo_exchange(HaloPlan const &plan, BoxGeometry const &box,
                   unsigned data_parts, ExchangeOp op);

} // namespace GhostComm

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

#include "particle_node.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_store/MigrationPack.hpp"
#include "particle_store/ParticleStore.hpp"
#include "rotation.hpp"
#include "system/System.hpp"

#include <utils/Cache.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/collectives/scatter.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

constexpr auto some_tag = 42;

/** @brief Mapping particle ids to MPI ranks. */
static std::unordered_map<int, int> particle_node;

static auto &get_cell_structure() {
  return *System::get_system().cell_structure;
}

/**
 * @brief Keep track of the largest particle id.
 * This book-keeping variable is necessary to make particle insertion run
 * in constant time. Traversing the @ref particle_node to find the largest
 * particle id scales with O(N) and traversing the local cells in parallel
 * followed by a reduction scales with O(N^2).
 */
static int max_seen_pid = -1;

static auto rebuild_needed() {
  auto is_rebuild_needed = ::particle_node.empty();
  boost::mpi::broadcast(::comm_cart, is_rebuild_needed, 0);
  return is_rebuild_needed;
}

static void mpi_synchronize_max_seen_pid_local() {
  boost::mpi::broadcast(::comm_cart, ::max_seen_pid, 0);
}

namespace {
/* Limit cache to 100 MiB */
std::size_t const max_cache_size = (100ul * 1048576ul) / sizeof(Particle);
Utils::Cache<int, Particle> particle_fetch_cache(max_cache_size);

/**
 * @brief Snapshot store backing the fetch cache.
 *
 * Every @ref Particle placed into @ref particle_fetch_cache is a view over
 * this head-node-local store. A particle received from a worker rank is
 * unpacked directly into a row here via @ref MigrationPack::unpack_rows; the
 * cached value is then a view over that row (all accessor reads go through
 * the store columns). Mirrors @c
 * Constraints::ShapeBasedConstraint::m_part_rep_store.
 *
 * Pointer-stability design (IMPORTANT):
 * @ref Utils::Cache is backed by a @c std::unordered_map<int,const Particle>,
 * whose nodes are stable across rehash but whose @c put() EVICTS a RANDOM
 * other element when the cache is full (@c drop_random_element). We therefore
 * do NOT keep raw @c Particle* addresses to re-assign on growth (an evicted
 * entry's address would dangle). Instead the store has a FIXED capacity per
 * invalidation epoch: it is lazily built once (on the first attach after an
 * invalidation) sized to the cache's capacity, and rows are handed out
 * monotonically. If the row counter would exceed the capacity, we simply drop
 * the whole cache and store (@ref reset_fetch_cache_store) and rebuild lazily
 * on the next attach — a cache is always safe to drop, and this avoids any
 * pointer bookkeeping or store-growth relocation. The capacity is capped so a
 * single epoch never over-allocates. */
ParticleStore fetch_cache_store;
/** Next free row in @ref fetch_cache_store; -1 while the store is unbuilt. */
int fetch_cache_store_next_row = -1;
/** Capacity @ref fetch_cache_store was built with in the current epoch. */
std::size_t fetch_cache_store_capacity = 0u;
/** Co-ownership of the Kokkos runtime (mirrors @ref CellStructure's
 *  @c m_kokkos_handle and the analogous member of
 *  @ref Constraints::ShapeBasedConstraint). @c fetch_cache_store
 *  holds Kokkos Views that must be released before @c Kokkos::finalize().
 *  Captured lazily when the store first allocates (particle_node code runs
 *  well after Kokkos init, but we capture defensively the same way). */
std::shared_ptr<KokkosHandle> fetch_cache_store_kokkos_handle;

/** Upper bound on rows allocated per invalidation epoch, so a single epoch
 *  never over-allocates the store. The store never needs more rows than the
 *  cache can hold entries, so this is min(cache capacity, a sane cap). */
std::size_t fetch_cache_store_capacity_target() {
  constexpr std::size_t cap = 65536u;
  return std::min(particle_fetch_cache.max_size(), cap);
}

/** Release the snapshot store's columns and reset the epoch. Kept in lock-step
 *  with @ref particle_fetch_cache invalidation: dropping cached particles whose
 *  store rows we release simultaneously is exactly correct. */
void reset_fetch_cache_store() {
  fetch_cache_store.release_columns();
  fetch_cache_store_next_row = -1;
  fetch_cache_store_capacity = 0u;
  // Drop the Kokkos runtime co-ownership: the store now holds no Views, so
  // keeping the handle would defer Kokkos::finalize() to this translation
  // unit's static teardown, where the store static is destroyed AFTER the
  // handle static (reverse construction order) and its (already released)
  // Views would otherwise be touched post-finalize. Releasing both together
  // keeps them in lock-step.
  fetch_cache_store_kokkos_handle.reset();
}

/** Allocate the fixed-capacity columns for a fresh epoch and reset the row
 *  counter to 0. Co-owns the Kokkos runtime so the columns can be released even
 *  if this translation unit outlives the last CellStructure (Kokkos is
 *  initialized by the time any particle fetch happens; capture defensively). */
void build_fetch_cache_store() {
  fetch_cache_store_kokkos_handle = ::kokkos_handle;
  fetch_cache_store_capacity = fetch_cache_store_capacity_target();
  fetch_cache_store.begin_rebuild(fetch_cache_store_capacity, 0u);
  fetch_cache_store.finish_rebuild();
  fetch_cache_store_next_row = 0;
}

/**
 * @brief Cache a fetched particle by unpacking its per-field buffer into
 * @ref fetch_cache_store.
 *
 * Grows nothing: the store is a fixed-capacity snapshot for the current
 * invalidation epoch (see the store's doc comment). Builds the store lazily on
 * first use, then hands out one monotonic row per call. The wire @p buffer
 * (produced by @ref mpi_send_particle_data_local's per-field pack) is unpacked
 * directly into that store row and a VIEW over the row is put into the cache.
 * If the store is exhausted mid-epoch, drops the cache and store and starts a
 * fresh epoch, then unpacks into row 0 of the new store. The returned pointer
 * aliases the cache's view entry and stays valid until the cache/store are next
 * invalidated.
 *
 * @param p_id    the fetched particle's id (the cache key).
 * @param buffer  the per-field packed row (exactly one row).
 * @returns       a pointer to the cached view (accessors read the store row).
 */
Particle const *cache_fetched_particle(int p_id,
                                       std::vector<char> const &buffer) {
  auto const exhausted = fetch_cache_store_next_row >= 0 and
                         static_cast<std::size_t>(fetch_cache_store_next_row) >=
                             fetch_cache_store_capacity;
  if (exhausted) {
    // Dropping the cache is always safe (it is a cache).
    particle_fetch_cache.invalidate();
    reset_fetch_cache_store();
  }
  if (fetch_cache_store_next_row < 0) {
    build_fetch_cache_store();
  }
  auto const row = fetch_cache_store_next_row++;
  MigrationPack::unpack_rows(fetch_cache_store, row, buffer);
  // Cache a VIEW over the freshly unpacked store row: its accessors read the
  // store columns directly.
  auto const cached =
      particle_fetch_cache.put(p_id, fetch_cache_store.make_view(row));
  return cached;
}

/**
 * @brief Cache a MULTI-row per-field buffer into @ref fetch_cache_store.
 *
 * Unpacks a buffer of @c n rows (produced by @ref mpi_get_particles_local's
 * per-field pack) into consecutive fetch-cache-store rows and caches a view per
 * row, keyed by the row's own id (read back from the store). Handles store
 * exhaustion the same way as @ref cache_fetched_particle (drop the cache/store
 * and start a fresh epoch). The caller (prefetch) never requests more ids than
 * the cache capacity, so one epoch always suffices; the exhaustion guard is
 * defensive.
 */
void cache_fetched_rows(std::vector<char> const &buffer) {
  if (buffer.empty()) {
    return;
  }
  std::uint64_t count = 0u;
  std::memcpy(&count, buffer.data(), sizeof(count));
  if (count == 0u) {
    return;
  }
  auto const would_exhaust = [&]() {
    return fetch_cache_store_next_row >= 0 and
           static_cast<std::size_t>(fetch_cache_store_next_row) +
                   static_cast<std::size_t>(count) >
               fetch_cache_store_capacity;
  };
  if (would_exhaust()) {
    particle_fetch_cache.invalidate();
    reset_fetch_cache_store();
  }
  if (fetch_cache_store_next_row < 0) {
    build_fetch_cache_store();
  }
  auto const first_row = fetch_cache_store_next_row;
  MigrationPack::unpack_rows(fetch_cache_store, first_row, buffer);
  fetch_cache_store_next_row += static_cast<int>(count);
  for (int k = 0; k < static_cast<int>(count); ++k) {
    auto const row = first_row + k;
    auto const id = fetch_cache_store.id(row);
    particle_fetch_cache.put(id, fetch_cache_store.make_view(row));
  }
}
} // namespace

void invalidate_fetch_cache() {
  particle_fetch_cache.invalidate();
  // Cache and snapshot store lifecycles stay in lock-step: rows we release
  // here back cached particles that we are dropping at the same time.
  reset_fetch_cache_store();
}
std::size_t fetch_cache_max_size() { return particle_fetch_cache.max_size(); }

static void mpi_send_particle_data_local(int p_id) {
  // Owner-side per-field pack: read the particle's live store row directly
  // into a byte buffer. The owning rank's store must be clean so the row is
  // valid. O(1) when clean.
  get_cell_structure().ensure_particle_store_synchronized();
  auto const p = get_cell_structure().get_local_particle(p_id);
  auto const found = p and not p->is_ghost();
  assert(1 == boost::mpi::all_reduce(::comm_cart, static_cast<int>(found),
                                     std::plus<>()) &&
         "Particle not found");
  if (found) {
    std::array<int, 1> const rows{{p->store_row()}};
    std::vector<char> buffer;
    MigrationPack::pack_rows(*p->store(), rows, buffer);
    ::comm_cart.send(0, 42, buffer);
  }
}

REGISTER_CALLBACK(mpi_send_particle_data_local)

Particle get_particle_data(int p_id) {
  auto const pnode = get_particle_node(p_id);

  if (pnode == this_node) {
    // The local path hands out a live particle whose accessor reads need valid
    // ParticleStore rows. O(1) when the store is clean; rank-local.
    // get_particle_data returns a by-value Particle VIEW (16-byte handle);
    // the view aliases the live store row and is valid for the caller's
    // immediate use.
    get_cell_structure().ensure_particle_store_synchronized();
    auto const p = get_cell_structure().get_local_particle(p_id);
    assert(p.has_value());
    return *p;
  }

  /* Query the cache */
  auto const p_ptr = particle_fetch_cache.get(p_id);
  if (p_ptr) {
    return *p_ptr;
  }

  /* Cache miss, fetch the particle,
   * put it into the cache and return a pointer into the cache. */
  Communication::mpiCallbacks().call_all(mpi_send_particle_data_local, p_id);
  // Receive the owner's per-field packed row and unpack it straight into the
  // head-node fetch-cache store; the cached value is a view over that row.
  std::vector<char> buffer;
  ::comm_cart.recv(boost::mpi::any_source, boost::mpi::any_tag, buffer);
  return *cache_fetched_particle(p_id, buffer);
}

static auto
get_local_particle_property(int p_id,
                            Utils::Vector3d (*getter)(Particle const &)) {
  // Ensure the owning rank's particle has a valid ParticleStore row before the
  // getter reads force/torque. O(1) when the store is clean; rank-local.
  get_cell_structure().ensure_particle_store_synchronized();
  auto const p = get_cell_structure().get_local_particle(p_id);
  auto const found = p and not p->is_ghost();
  assert(1 == boost::mpi::all_reduce(::comm_cart, static_cast<int>(found),
                                     std::plus<>()) &&
         "particle not found exactly once");
  auto const local_value = found ? getter(*p) : Utils::Vector3d{};
  return boost::mpi::all_reduce(::comm_cart, local_value, std::plus<>());
}

static void mpi_get_particle_force_local(int p_id) {
  get_local_particle_property(
      p_id, [](Particle const &p) { return Utils::Vector3d(p.force()); });
}

REGISTER_CALLBACK(mpi_get_particle_force_local)

Utils::Vector3d get_particle_force(int p_id) {
  Communication::mpiCallbacks().call(mpi_get_particle_force_local, p_id);
  return get_local_particle_property(
      p_id, [](Particle const &p) { return Utils::Vector3d(p.force()); });
}

#ifdef ESPRESSO_ROTATION
static void mpi_get_particle_torque_lab_local(int p_id) {
  get_local_particle_property(p_id, [](Particle const &p) {
    return convert_vector_body_to_space(p, Utils::Vector3d(p.torque()));
  });
}

REGISTER_CALLBACK(mpi_get_particle_torque_lab_local)

Utils::Vector3d get_particle_torque_lab(int p_id) {
  Communication::mpiCallbacks().call(mpi_get_particle_torque_lab_local, p_id);
  return get_local_particle_property(p_id, [](Particle const &p) {
    return convert_vector_body_to_space(p, Utils::Vector3d(p.torque()));
  });
}
#endif // ESPRESSO_ROTATION

static void mpi_get_particles_local() {
  std::vector<int> local_ids;
  boost::mpi::scatter(comm_cart, local_ids, 0);

  // Per-field pack of the requested live store rows into one byte buffer;
  // gathered to the head as a per-rank buffer.
  auto &cell_structure = get_cell_structure();
  std::vector<int> rows;
  rows.reserve(local_ids.size());
  for (auto const p_id : local_ids) {
    auto const p = cell_structure.get_local_particle(p_id);
    assert(p.has_value());
    rows.push_back(p->store_row());
  }
  std::vector<char> buffer;
  MigrationPack::pack_rows(cell_structure.particle_store(), rows, buffer);
  boost::mpi::gather(comm_cart, buffer, 0);
}

REGISTER_CALLBACK(mpi_get_particles_local)

/**
 * @brief Fetch multiple particles into the head-node fetch cache at once.
 *
 * Groups the requested ids per owning rank, scatters the id lists, and gathers
 * one per-field packed buffer per rank; each buffer is unpacked straight into
 * @c fetch_cache_store (a view per row is cached, keyed by the row's id).
 * Particles are cached in an arbitrary (per-rank) order; every caller reads the
 * cache by id afterwards, so the order does not matter.
 *
 * @param ids The ids of the particles that should be fetched (none local).
 */
static void mpi_get_particles(std::span<const int> ids) {
  Communication::mpiCallbacks().call(mpi_get_particles_local);

  /* Group ids per node */
  static std::vector<std::vector<int>> node_ids(comm_cart.size());
  for (auto &per_node : node_ids) {
    per_node.clear();
  }

  for (auto const &p_id : ids) {
    auto const p_node = get_particle_node(p_id);
    node_ids[p_node].emplace_back(p_id);
  }
  /* We shouldn't be prefetching particles that are already on the head node */
  assert(node_ids[this_node].empty());

  /* Distribute ids to the worker nodes */
  {
    static std::vector<int> ignore;
    boost::mpi::scatter(comm_cart, node_ids, ignore, 0);
    assert(ignore.empty());
  }

  // Gather one per-field packed buffer per rank and unpack each into the fetch
  // cache. The head's own buffer is empty (it requests nothing of itself).
  std::vector<std::vector<char>> node_buffers(comm_cart.size());
  std::vector<char> const empty_buffer;
  boost::mpi::gather(comm_cart, empty_buffer, node_buffers, 0);
  for (auto const &buffer : node_buffers) {
    cache_fetched_rows(buffer);
  }
}

void prefetch_particle_data(std::span<const int> in_ids) {
  /* Nothing to do on a single node. */
  // NOLINTNEXTLINE(clang-analyzer-core.NonNullParamChecker)
  if (comm_cart.size() == 1)
    return;

  static std::vector<int> ids;
  ids.clear();
  auto out_ids = std::back_inserter(ids);

  /* Don't prefetch particles already on the head node or already cached. */
  // TODO(upstream): inverted predicate -- `has(id)` should be `!has(id)`; this
  // never prefetches uncached particles. Pre-existing upstream bug in a cold
  // multi-rank path; fix in a separate upstream PR.
  std::ranges::copy_if(in_ids, out_ids, [](int id) {
    return (get_particle_node(id) != this_node) && particle_fetch_cache.has(id);
  });

  /* Don't prefetch more particles than fit the cache. */
  if (ids.size() > particle_fetch_cache.max_size())
    ids.resize(particle_fetch_cache.max_size());

  /* Fetch the particles (unpacked directly into the fetch cache). */
  mpi_get_particles(ids);
}

static void mpi_who_has_local() {
  static std::vector<int> sendbuf;

  auto local_particles = get_cell_structure().local_particles();
  auto const n_part = static_cast<int>(local_particles.size());
  boost::mpi::gather(comm_cart, n_part, 0);

  if (n_part == 0) {
    mpi_synchronize_max_seen_pid_local();
    return;
  }

  sendbuf.resize(n_part);

  std::transform(local_particles.begin(), local_particles.end(),
                 sendbuf.begin(), [](Particle const &p) { return p.id(); });

  comm_cart.send(0, some_tag, sendbuf);
  mpi_synchronize_max_seen_pid_local();
}

REGISTER_CALLBACK(mpi_who_has_local)

static void mpi_who_has_head() {
  auto local_particles = get_cell_structure().local_particles();

  static std::vector<int> n_parts;
  boost::mpi::gather(comm_cart, static_cast<int>(local_particles.size()),
                     n_parts, 0);

  static std::vector<int> pdata;
  auto const n_nodes = ::comm_cart.size();
  max_seen_pid = -1;

  /* then fetch particle locations */
  for (int pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      for (auto const &p : local_particles) {
        particle_node[p.id()] = this_node;
        max_seen_pid = std::max(max_seen_pid, p.id());
      }
    } else if (n_parts[pnode] > 0) {
      pdata.resize(n_parts[pnode]);
      comm_cart.recv(pnode, some_tag, pdata);
      for (int i = 0; i < n_parts[pnode]; i++) {
        particle_node[pdata[i]] = pnode;
        max_seen_pid = std::max(max_seen_pid, pdata[i]);
      }
    }
  }
  mpi_synchronize_max_seen_pid_local();
}

/**
 * @brief Rebuild the particle index.
 */
static void build_particle_node() {
  Communication::mpiCallbacks().call(mpi_who_has_local);
  mpi_who_has_head();
}

/**
 * @brief Rebuild the particle index.
 */
static void build_particle_node_parallel() {
  if (this_node == 0) {
    mpi_who_has_head();
  } else {
    mpi_who_has_local();
  }
}

int get_particle_node(int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }

  if (particle_node.empty())
    build_particle_node();

  auto const needle = particle_node.find(p_id);

  // Check if particle has a node, if not, we assume it does not exist.
  if (needle == particle_node.end()) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
  return needle->second;
}

int get_particle_node_parallel(int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }

  if (rebuild_needed()) {
    build_particle_node_parallel();
  }

  if (this_node != 0) {
    return -1;
  }

  auto const needle = particle_node.find(p_id);

  // Check if particle has a node, if not, we assume it does not exist.
  if (needle == particle_node.end()) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
  return needle->second;
}

void clear_particle_node() {
  ::max_seen_pid = -1;
  particle_node.clear();
}

/**
 * @brief Calculate the largest particle id.
 * Traversing the @ref particle_node to find the largest particle id
 * scales with O(N). Consider using the cached value in @ref max_seen_pid
 * if possible. This function is only necessary when the cached value is
 * invalidated, for example when removing the particle which has the
 * largest id.
 */
static int calculate_max_seen_id() {
  return std::accumulate(
      particle_node.begin(), particle_node.end(), -1,
      [](int max, auto const &kv) { return std::max(max, kv.first); });
}

/**
 * @brief Create a new particle and attach it to a cell.
 * @param p_id  The identity of the particle to create.
 * @param pos   The particle position.
 * @return Whether the particle was created on that node.
 */
static bool maybe_insert_particle(int p_id, Utils::Vector3d const &pos) {
  auto const &box_geo = *System::get_system().box_geo;
  auto folded_pos = pos;
  auto image_box = Utils::Vector3i{};
  box_geo.fold_position(folded_pos, image_box);

  // Build the new particle into a staging-store row via a view, then hand it
  // to add_local_particle, which stages the underlying row into the home cell.
  auto &cell_structure = get_cell_structure();
  auto new_part = cell_structure.make_new_particle_view();
  new_part.id() = p_id;
  new_part.pos() = folded_pos;
  new_part.image_box() = image_box;

  return cell_structure.add_local_particle(std::move(new_part)).has_value();
}

/**
 * @brief Move particle to a new position.
 * @param p_id  The identity of the particle to move.
 * @param pos   The new particle position.
 * @return Whether the particle was moved from that node.
 */
static bool maybe_move_particle(int p_id, Utils::Vector3d const &pos) {
  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto p = system.cell_structure->get_local_particle(p_id);
  if (not p) {
    return false;
  }
  auto folded_pos = pos;
  auto image_box = Utils::Vector3i{};
  box_geo.fold_position(folded_pos, image_box);
  p->pos() = folded_pos;
  p->image_box() = image_box;
  return true;
}

void remove_particle(int p_id) {
  if (this_node == 0) {
    particle_node[p_id] = -1;
  }
  get_cell_structure().remove_particle(p_id);
  System::get_system().on_particle_change();
  mpi_synchronize_max_seen_pid_local();
  if (this_node == 0) {
    particle_node.erase(p_id);
    if (p_id == ::max_seen_pid) {
      --::max_seen_pid;
      // if there is a gap (i.e. there is no particle with id max_seen_pid - 1,
      // then the cached value is invalidated and has to be recomputed (slow)
      if (not particle_node.contains(::max_seen_pid) or
          particle_node[::max_seen_pid] == -1) {
        ::max_seen_pid = calculate_max_seen_id();
      }
    }
  }
  mpi_synchronize_max_seen_pid_local();
  get_cell_structure().resort_particles(false);
}

void make_new_particle(int p_id, Utils::Vector3d const &pos) {
  if (rebuild_needed()) {
    build_particle_node_parallel();
  }
  auto const has_created = maybe_insert_particle(p_id, pos);
  System::get_system().on_particle_change();

  auto node = -1;
  auto const node_local = (has_created) ? ::comm_cart.rank() : 0;
  boost::mpi::reduce(::comm_cart, node_local, node, std::plus<int>{}, 0);
  if (::this_node == 0) {
    particle_node[p_id] = node;
    max_seen_pid = std::max(max_seen_pid, p_id);
    assert(not has_created or node == 0);
  }
  mpi_synchronize_max_seen_pid_local();
}

void set_particle_pos(int p_id, Utils::Vector3d const &pos) {
  auto const has_moved = maybe_move_particle(p_id, pos);
  get_cell_structure().set_resort_particles(Cells::RESORT_GLOBAL);
  System::get_system().on_particle_change();

  auto success = false;
  boost::mpi::reduce(::comm_cart, has_moved, success, std::plus<bool>{}, 0);
  if (::this_node == 0 and !success) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
}

bool particle_exists(int p_id) {
  if (particle_node.empty())
    build_particle_node();
  return particle_node.contains(p_id);
}

std::vector<int> get_particle_ids() {
  if (particle_node.empty())
    build_particle_node();

  std::vector<int> pids{};
  std::ranges::copy(std::views::keys(particle_node), std::back_inserter(pids));
  std::ranges::sort(pids);

  return pids;
}

std::vector<int> get_particle_ids_parallel() {
  if (rebuild_needed()) {
    build_particle_node_parallel();
  }
  std::vector<int> pids{};
  std::ranges::copy(std::views::keys(particle_node), std::back_inserter(pids));
  boost::mpi::broadcast(::comm_cart, pids, 0);
  return pids;
}

int get_maximal_particle_id() {
  if (rebuild_needed()) {
    build_particle_node_parallel();
  }

  return max_seen_pid;
}

int get_n_part() {
  if (particle_node.empty())
    build_particle_node();

  return static_cast<int>(particle_node.size());
}

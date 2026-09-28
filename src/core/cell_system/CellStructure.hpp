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

#include "cell_system/MigrationStaging.hpp"
#include "cell_system/ParticleDecomposition.hpp"
#include "cell_system/ParticleListOperations.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"
#include "ParticleRange.hpp"
#include "algorithm/link_cell.hpp"
#include "bond_error.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/CellStructureType.hpp"
#include "config/config.hpp"
#include "custom_verlet_list.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloExchange.hpp"
#include "particle_store/ParticleStore.hpp"
#include "particle_store/StoreGenerationGuard.hpp"
#include "system/Leaf.hpp"

#include <utils/Vector.hpp>

#include <boost/container/static_vector.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/range/algorithm/transform.hpp>

#include <Cabana_Core.hpp>
#include <Cabana_NeighborList.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

// forward declarations
struct KokkosHandle;
struct LocalBondState;

template <typename Callable>
concept ParticleCallback = requires(Callable c, Particle &p) {
  { c(p) } -> std::same_as<void>;
};

namespace Cells {
enum Resort : unsigned {
  RESORT_NONE = 0u,
  RESORT_LOCAL = 1u,
  RESORT_GLOBAL = 2u
};

/**
 * @brief Flags to select particle parts for communication.
 */
enum DataPart : unsigned {
  DATA_PART_NONE = 0u,       /**< Nothing */
  DATA_PART_PROPERTIES = 1u, /**< Particle::p */
  DATA_PART_POSITION = 2u,   /**< Particle::r */
  DATA_PART_MOMENTUM = 8u,   /**< Particle::m */
  DATA_PART_FORCE = 16u,     /**< Particle::f */
#ifdef ESPRESSO_BOND_CONSTRAINT
  DATA_PART_RATTLE = 32u, /**< Particle::rattle */
#endif
  DATA_PART_BONDS = 64u, /**< Particle::bonds */
#ifdef ESPRESSO_ROTATION
  DATA_PART_QUAT = 128u,   /**< orientation quaternion (pushed with position) */
  DATA_PART_TORQUE = 256u, /**< torque (reduced with force) */
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  DATA_PART_DIPFLD = 512u, /**< Particle::dip_fld */
#endif
};
} // namespace Cells

/**
 * @brief Map the data parts flags from cells to those
 *        used internally by the ghost communication.
 *
 * @param data_parts data parts flags
 * @return ghost communication flags
 */
unsigned map_data_parts(unsigned data_parts);

namespace Cells {
inline ParticleRange particles(std::span<Cell *const> cells) {
  /* Find first non-empty cell */
  auto first_non_empty = std::ranges::find_if(
      cells, [](const Cell *c) { return not c->particles().empty(); });

  return {CellParticleIterator(first_non_empty, cells.end()),
          CellParticleIterator(cells.end())};
}
} // namespace Cells

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Utils::Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Utils::Vector3d vec21;
  double dist2;
};

namespace detail {
// NOLINTNEXTLINE(bugprone-exception-escape)
struct MinimalImageDistance {
  BoxGeometry const box;

  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(box.get_mi_vector(Utils::Vector3d(p1.pos()),
                                      Utils::Vector3d(p2.pos())));
  }
};

struct EuclidianDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(Utils::Vector3d(p1.pos()) - Utils::Vector3d(p2.pos()));
  }
};
} // namespace detail

#ifdef ESPRESSO_CUDA
/** @brief Opaque persistent device buffers for the GPU short-range pair loop;
 *  defined in forces_lj_device.cu. Forward-declared so no Kokkos headers reach
 *  this widely-included header (that perturbs host FP codegen). */
struct DeviceShortRangeBuffers;
#endif

/** Describes a cell structure / cell system. Contains information
 *  about the communication of cell contents (particles, ghosts, ...)
 *  between different nodes and the relation between particle
 *  positions and the cell system. All other properties of the cell
 *  system which are not common between different cell systems have to
 *  be stored in separate structures.
 */
class CellStructure : public System::Leaf<CellStructure> {
public:
  static constexpr auto vector_length = 1;
  using memory_space = Kokkos::HostSpace;
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  struct AoSoA_pack;
  using ForceType =
      Kokkos::View<double *[3], Kokkos::LayoutRight, memory_space>;
  using VirialType = Kokkos::View<double[3], Kokkos::LayoutRight, memory_space>;
  using ScatterForce =
      Kokkos::Experimental::ScatterView<double *[3], Kokkos::LayoutRight,
                                        memory_space>;
  using ScatterVirial =
      Kokkos::Experimental::ScatterView<double[3], Kokkos::LayoutRight,
                                        memory_space>;
  using ListAlgorithm = Cabana::HalfNeighborTag;
  using ListType =
      CustomVerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D,
                       Cabana::TeamVectorOpTag>;

private:
  /** The local id -> @ref ParticleStore ROW map. Indexed by particle id; the
   *  entry is the store row that @ref get_local_particle resolves the id to, or
   *  @c no_store_row (-1) when the id is not present on this rank. Locals win
   *  over ghost copies of the same id; among ghost copies the first valid row
   *  (in store-row order) wins. Rebuilt wholesale from the (synchronized) store
   *  by @c rebuild_particle_index / @c index_ghost_particles at store-rebuild
   *  cadence (the single write site). An entry stays valid between store
   *  rebuilds; a rebuild renumbers rows and refreshes the whole map.
   *
   *  @ref get_local_particle resolves the id to a row here and returns a fresh
   *  by-value @ref Particle view over that row (a 16-byte handle). */
  std::vector<int> m_id_to_store_row;
  /** Sentinel stored in @c m_id_to_store_row for an id absent on this rank.
   */
  static constexpr int no_store_row = -1;
  /** Implementation of the primary particle decomposition */
  std::unique_ptr<ParticleDecomposition> m_decomposition;
  /** Active type in m_decomposition */
  CellStructureType m_type = CellStructureType::NSQUARE;
  /** One of @ref Cells::Resort, announces the level of resort needed.
   */
  unsigned m_resort_particles = Cells::RESORT_NONE;
  bool m_verlet_skin_set = false;
  bool m_rebuild_verlet_list = true;
  bool m_rebuild_verlet_list_cabana = true;
  /** Monotonic counter bumped each time the Cabana Verlet list is rebuilt; lets
   *  the GPU short-range path copy the list to the device only on rebuild. */
  std::uint64_t m_verlet_list_cabana_generation = 0u;
  /** Interaction pairs as @ref ParticleStore ROW indices. Held across
   *  integration steps until the next Verlet rebuild. Cells do not own stable
   *  @c Particle addresses, so the pairs record the two particles' store rows;
   *  the pair kernels resolve each row back to a view at loop entry (hoisted
   *  store pointer, one view per row). The rows are only valid for the store
   *  @ref ParticleStore::generation() they were recorded at (a rebuild
   *  renumbers rows); that generation is stamped in @c m_verlet_list_store
   *  / @c m_verlet_list_generation and re-checked (debug) before every use
   *  via @ref ParticleStoreGuard::assert_generation. */
  std::vector<std::pair<int, int>> m_verlet_list;
  /** Store identity + generation the rows in @c m_verlet_list were recorded
   *  at. Used only by the debug generation guard; a rebuild between build and
   *  consume without a Verlet rebuild would make the rows stale. */
  ParticleStore const *m_verlet_list_store = nullptr;
  std::uint64_t m_verlet_list_generation = 0u;
  double m_le_pos_offset_at_last_resort = 0.;
  /** @brief Verlet list skin. */
  double m_verlet_skin = 0.;
  double m_verlet_reuse = 0.;
  int m_cached_max_local_particle_id = 0;
  std::size_t m_num_local_particles_cached = 0;
  int m_max_id = 0;
  std::unique_ptr<Kokkos::View<int *, memory_space>> m_id_to_index;
  std::unique_ptr<ForceType> m_local_force;
  std::optional<ScatterForce> m_scatter_force;
#ifdef ESPRESSO_ROTATION
  std::unique_ptr<ForceType> m_local_torque;
  std::optional<ScatterForce> m_scatter_torque;
  /**
   * @brief True if a kernel may have written to @c m_scatter_torque since
   * the last reset. Cleared by the resets; when false, the torque buffers
   * are known to be all-zero and their O(n_threads * N) zeroing and
   * reduction can be skipped.
   */
  bool m_torque_replicas_dirty = false;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  std::unique_ptr<ForceType> m_local_dip_fld;
  std::optional<ScatterForce> m_scatter_dip_fld;
  /** @brief Same contract as @c m_torque_replicas_dirty, for dipole fields. */
  bool m_dip_fld_replicas_dirty = false;
#endif
#ifdef ESPRESSO_NPT
  std::unique_ptr<VirialType> m_local_virial;
  std::optional<ScatterVirial> m_scatter_virial;
  /** @brief Same contract as @c m_torque_replicas_dirty, for the virial. */
  bool m_virial_replicas_dirty = false;
#endif
  std::unique_ptr<LocalBondState> m_bond_state;
  std::unique_ptr<ListType> m_verlet_list_cabana;
#ifdef ESPRESSO_CUDA
  /** Persistent device buffers for the GPU short-range pair loop (opaque;
   *  defined in forces_lj_device.cu). Reused across steps to avoid per-step
   *  device allocations; destroyed before Kokkos::finalize (this owns a
   *  KokkosHandle via the cell structure lifetime). */
  std::shared_ptr<DeviceShortRangeBuffers> m_device_sr_buffers;
#endif
  /**
   * @brief Persistent per-neighbor buffer pool for ghost exchanges.
   *
   * Reused across calls to @ref ghosts_count / @ref ghosts_update /
   * @ref ghosts_reduce_forces so that, after the first (warm-up) exchange,
   * the underlying @c std::vector storage is retained and no per-step heap
   * allocation occurs on the hot ghost-communication path.
   *
   * Mutable because the ghost methods are logically const with respect to
   * particle data but mutate this scratch storage.
   */
  mutable GhostComm::ExchangeBuffers m_ghost_buffers;
  /**
   * @brief Scratch cell-pointer list for filtered particle iteration.
   *
   * Used by @c for_each_interior_particle and @c for_each_boundary_particle
   * to hold the filtered subset of local cells.  Declared mutable so that the
   * const-qualified iteration helpers can write to it; not thread-safe — the
   * two filtered passes must not run concurrently (they never do: interior pass
   * completes before the ghost-reduce finish, which precedes the boundary
   * pass).
   */
  mutable std::vector<Cell *> m_filtered_cells_scratch;
  /**
   * @brief In-flight ghost force reduction started by
   *        @ref ghosts_reduce_forces_start.
   *
   * When set, @ref ghosts_reduce_forces_finish must be called before any
   * resort or decomposition change.  Stored as optional so that the
   * unfinished state is detectable at runtime.
   */
  mutable std::optional<GhostComm::GhostExchange> m_pending_ghost_reduce;
  /** particle properties using individual Kokkos Views */
  std::unique_ptr<AoSoA_pack> m_aosoa;
  /** Pack-ordered list of the particles that participate in the pack /
   *  force-scatter kernels: the local prefix (store rows [0, n_local) in row
   *  order) followed by the deduped ghost tail (first occurrence of each ghost
   *  id not owned by a local). Consumers (force/energy/pressure/icc reductions,
   *  the AoSoA commit) index it in pack order. Pointers into
   *  @c m_unique_particle_views; rebuilt every @ref set_index_map.
   *  Required for the same reason as @c m_pack_index_to_store_row: the deduped
   *  ghost tail is a non-contiguous subset, so pack index != store row on the
   *  tail and the view list is needed. It can be retired only once the ghost
   *  dedup produces a contiguous pack-order ghost range. */
  std::vector<Particle *> m_unique_particles;
  /** Owning backing for @c m_unique_particles: the storage the pack-order
   * pointer list points into -- one by-value @ref Particle view per pack
   * participant, in pack order, rebuilt wholesale by @ref set_index_map. A
   * @c std::vector reused (grown, never shrunk) across rebuilds; the pointers
   * in @c m_unique_particles are only valid until the next @ref set_index_map.
   */
  std::vector<Particle> m_unique_particle_views;
  std::shared_ptr<KokkosHandle> m_kokkos_handle;
  /** Array-based particle storage. */
  ParticleStore m_particle_store;
  /** Migration staging store. A second, small @ref ParticleStore holding rows
   *  in transit between the live store and a remote rank: `extract` copies a
   *  live row into a staging row (@c stage_row); a migrating row that a round
   *  could not deliver locally waits here between exchange rounds; `receive`
   *  unpacks into staging rows; `rebuild` commits them into the live store.
   *  Lazily sized on first @c stage_row and reused (grown, never shrunk) across
   *  resorts, like the fetch-cache store. */
  ParticleStore m_staging_store;
  /** Next free row in @c m_staging_store; the count of currently-staged rows.
   *  0 while empty; reset by @ref clear_staging_store. */
  int m_staging_store_next_row = 0;
  /** Capacity @c m_staging_store was last (re)built with. */
  std::size_t m_staging_store_capacity = 0u;
  /** Pack-index -> store-row translation. Identity on the local prefix; only
   *  the deduped ghost tail is remapped. Rebuilt in @ref set_index_map.
   *  Identity holds on the local prefix by construction; the ghost tail is a
   *  DEDUPED, non-contiguous subset of [n_local, n_total), so the real
   *  translation is required. It can be removed only once the ghost dedup is
   *  expressed as a contiguous row range. */
  Kokkos::View<int *, Kokkos::HostSpace> m_pack_index_to_store_row;
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  /** Store-side derived director view, sized n_total, recomputed from the
   *  quaternion column by @ref update_director_view. Uses the store's vector
   *  layout (@ref ParticleStore::StateVectorLayout) so it aliases into the
   *  pack's DirectorViewType. */
  Kokkos::View<double *[3], ParticleStore::StateVectorLayout, Kokkos::HostSpace>
      m_director_view;
#endif

public:
  CellStructure(BoxGeometry const &box);
  virtual ~CellStructure();

  bool use_verlet_list = true;

  /**
   * @brief Set the store-row entry for a particle id.
   *
   * Records that @p id resolves to store @p row in @c m_id_to_store_row,
   * growing the map as needed. The single write site is the store-rebuild
   * indexing (@c rebuild_particle_index / @c index_ghost_particles).
   *
   * @param id  Particle id (>= 0).
   * @param row Store row the id resolves to.
   */
  void update_particle_index(int id, int row) {
    assert(id >= 0);
    assert(row >= 0);

    if (static_cast<unsigned int>(id) >= m_id_to_store_row.size())
      m_id_to_store_row.resize(static_cast<unsigned int>(id + 1), no_store_row);

    m_id_to_store_row[static_cast<unsigned int>(id)] = row;
  }

  /**
   * @brief Clear the id -> store-row map.
   */
  void clear_particle_index() { m_id_to_store_row.clear(); }

private:
  /**
   * @brief Rebuild the id -> store-row map from the (synchronized) store.
   *
   * Fills @c m_id_to_store_row with the store row of each indexed particle.
   * Must run with a clean store (rows assigned). Indexes LOCALS only; ghost id
   * columns are not valid until ghosts_update runs, so ghosts are indexed
   * separately by @c index_ghost_particles afterwards. "Locals win" over ghost
   * copies of the same id.
   */
  void rebuild_particle_index();

  /**
   * @brief Append the ghost rows to the id -> store-row map.
   *
   * Run AFTER ghosts_update has filled the ghost id columns (a fresh ghost row
   * carries a default id until then). Records the store row for each ghost
   * whose id is not owned by a local (or by an earlier ghost of the same id --
   * first valid row wins).
   */
  void index_ghost_particles();

  /**
   * @brief Build the resort permutation over surviving store rows.
   *
   * Walks the cells in the rebuild order -- local cells in
   * @ref ParticleDecomposition::local_cells span order, then ghost cells; per
   * cell: surviving rows in raw range order (skipping any pending-removed row),
   * then staged rows in @c cell->staged() push order -- and produces:
   *  - @p permutation : one entry per new row. A surviving row's entry is the
   *    OLD store row whose data the permute rebuild moves into that new row; a
   *    staged / fresh-ghost row's entry is -1 (the permute rebuild seeds the
   *    defaults, the caller overwrites a staged local via copy_row).
   *  - @p cell_ranges : the future @c (offset, count) of each cell, in the same
   *    cell order (locals then ghosts). @c offset is the first new row of the
   *    cell, @c count its surviving + staged row total. This is the (offset,
   *    count) collapse @ref ensure_particle_store_synchronized writes back onto
   *    @ref Cell.
   *
   * Iteration order is cells in span order, surviving rows in range order,
   * staged rows in push order. @ref ensure_particle_store_synchronized feeds
   * @p permutation to @ref ParticleStore::permute_rebuild.
   */
  void build_resort_permutation(
      std::vector<int> &permutation,
      std::vector<std::pair<std::size_t, std::size_t>> &cell_ranges) const;

  /** @brief Grow the migration staging store to hold at least @p needed rows,
   *  preserving already-staged rows. Shared by @ref stage_row and
   *  @ref reserve_staging_rows. */
  void ensure_staging_capacity(std::size_t needed);

  /**
   * @brief Stage a new particle (built into a staging row) into a cell.
   *
   * @p p is a VIEW over a row of the creation staging store (built by
   * @ref make_new_particle_view, then populated by the caller). The referenced
   * staging row is staged into @p cell (@ref
   * CellParticleStorage::insert_staged_row); the next rebuild copies it into a
   * committed row. The staging store must stay valid until the rebuild runs
   * (the caller commits immediately via ensure_particle_store_synchronized).
   * The caller must mark the store dirty.
   */
  void append_staged_particle(Cell &cell, Particle &&p) {
    assert(p.store() != nullptr);
    CellParticleStorage::insert_staged_row(cell, *p.store(), p.store_row());
  }

public:
  /**
   * @brief Get a local particle by id.
   *
   * Resolves @p id to a store row via @c m_id_to_store_row and returns a
   * fresh by-value @ref Particle view over that row. The returned view is a
   * 16-byte handle aliasing the store; it is valid until the next store rebuild
   * (which renumbers rows). An absent id yields an empty optional (the nullptr
   * equivalent).
   *
   * @param id Particle to get.
   * @return A view of the particle if it is local (or an indexed ghost),
   *         @c std::nullopt otherwise.
   */
  std::optional<Particle> get_local_particle(int id) {
    assert(id >= 0);

    if (static_cast<unsigned int>(id) >= m_id_to_store_row.size())
      return std::nullopt;

    auto const row = m_id_to_store_row[static_cast<unsigned int>(id)];
    if (row == no_store_row)
      return std::nullopt;

    return m_particle_store.make_view(row);
  }

  /** @overload */
  std::optional<const Particle> get_local_particle(int id) const {
    assert(id >= 0);

    if (static_cast<unsigned int>(id) >= m_id_to_store_row.size())
      return std::nullopt;

    auto const row = m_id_to_store_row[static_cast<unsigned int>(id)];
    if (row == no_store_row)
      return std::nullopt;

    // const overload: the store is logically const here, but make_view needs a
    // mutable store to bind the view. The returned view is const, so no
    // mutation escapes.
    return const_cast<ParticleStore &>(m_particle_store).make_view(row);
  }

  template <class InputRange, class OutputIterator>
  void get_local_particles(InputRange ids, OutputIterator out) {
    std::ranges::transform(ids, out,
                           [this](int id) { return get_local_particle(id); });
  }

  CellStructureType decomposition_type() const { return m_type; }

  /** Maximal cutoff supported by current cell system. */
  Utils::Vector3d max_cutoff() const { return decomposition().max_cutoff(); }

  /** Maximal pair range supported by current cell system. */
  Utils::Vector3d max_range() const { return decomposition().max_range(); }

  ParticleRange local_particles() const {
    return Cells::particles(decomposition().local_cells());
  }

  ParticleRange ghost_particles() const {
    return Cells::particles(decomposition().ghost_cells());
  }

  std::size_t count_local_particles() const {
    std::size_t count = 0;
    for (auto const &cell : m_decomposition->local_cells()) {
      // Count committed rows directly: avoids requiring the cell's store
      // pointer to be wired, and is exactly the number of committed particles
      // (staged-but-uncommitted particles are not counted).
      count += cell->count();
    }
    return count;
  }

  /** @brief whether to use parallel version of @ref for_each_local_particle */
  bool use_parallel_for_each_local_particle() const { return true; }

  /**
   * @brief Run a kernel on all local particles.
   * The kernel is assumed to be thread-safe.
   */
  void for_each_local_particle(ParticleCallback auto &&f,
                               bool parallel = true) const {
    if (parallel and use_parallel_for_each_local_particle()) {
      parallel_for_each_particle_impl(decomposition().local_cells(), f);
      return;
    }
    for (auto &p : local_particles()) {
      f(p);
    }
  }

  /**
   * @brief Run a column kernel on every local particle by STORE ROW.
   *
   * A direct-column launcher for the hot integrator/thermostat loops. Instead
   * of materialising a @ref Particle per element (which pays @c view_host() +
   * address + stride on every accessor), this launcher iterates exactly the
   * same rows in exactly the same order as @ref for_each_local_particle but
   * hands the kernel the raw STORE ROW @c int. The kernel captures its hoisted
   * @c *_view() column handles ONCE (outside the loop) and indexes them by row.
   *
   * Iteration structure matches @c parallel_for_each_particle_impl (multi-cell:
   * parallel over cells, inner serial over @c [offset, offset+count);
   * single-cell: parallel over the cell's rows). Local cells tile
   * @c [0, n_local) contiguously in cell-traversal order (see
   * @ref ensure_particle_store_synchronized), so the visited row set and order
   * match the view path. Virtual sites are local rows and ARE visited -- the
   * kernel guards them per-row (propagation mask), matching the view-path
   * lambda's @c is_virtual() early return. The kernel is assumed to be
   * thread-safe.
   */
  template <typename RowKernel>
  void for_each_local_particle_row(RowKernel &&kernel) const {
    parallel_for_each_local_row_impl(decomposition().local_cells(), kernel);
  }

  /**
   * @brief Run a kernel on interior (non-boundary) local particles only.
   *
   * A cell is interior iff none of its neighbors is a ghost cell
   * (@ref Cell::is_boundary returns false).  This filtered variant uses the
   * same Kokkos parallelization as @ref for_each_local_particle — it builds
   * a filtered cell list and delegates to @c parallel_for_each_particle_impl
   * so that thread-level behavior is identical to the unfiltered path.
   *
   * The kernel is assumed to be thread-safe.
   */
  void for_each_interior_particle(ParticleCallback auto &&f) const {
    auto const all_cells = decomposition().local_cells();
    // Build a filtered list of interior cells (those that are not boundary).
    m_filtered_cells_scratch.clear();
    for (auto *c : all_cells) {
      if (not c->is_boundary())
        m_filtered_cells_scratch.push_back(c);
    }
    if (m_filtered_cells_scratch.empty())
      return;
    std::span<Cell *const> span{m_filtered_cells_scratch};
    if (use_parallel_for_each_local_particle()) {
      parallel_for_each_particle_impl(span, f);
    } else {
      for (auto *c : span)
        for (auto &p : c->particles())
          f(p);
    }
  }

  /**
   * @brief Run a kernel on boundary local particles only.
   *
   * Complement of @c for_each_interior_particle: visits particles in cells
   * where @ref Cell::is_boundary returns true.  Uses the same Kokkos
   * parallelization structure as @ref for_each_local_particle.
   *
   * The kernel is assumed to be thread-safe.
   */
  void for_each_boundary_particle(ParticleCallback auto &&f) const {
    auto const all_cells = decomposition().local_cells();
    // Build a filtered list of boundary cells.
    m_filtered_cells_scratch.clear();
    for (auto *c : all_cells) {
      if (c->is_boundary())
        m_filtered_cells_scratch.push_back(c);
    }
    if (m_filtered_cells_scratch.empty())
      return;
    std::span<Cell *const> span{m_filtered_cells_scratch};
    if (use_parallel_for_each_local_particle()) {
      parallel_for_each_particle_impl(span, f);
    } else {
      for (auto *c : span)
        for (auto &p : c->particles())
          f(p);
    }
  }

  /**
   * @brief Run a kernel on all ghost particles.
   * The kernel is assumed to be thread-safe.
   */
  void for_each_ghost_particle(ParticleCallback auto &&f) const {
    for (auto &p : ghost_particles()) {
      f(p);
    }
  }

private:
  /** Cell system dependent function to find the right cell for a
   *  particle.
   *  \param  p Particle.
   *  \return pointer to cell where to put the particle, nullptr
   *          if the particle does not belong on this node.
   */
  Cell *particle_to_cell(const Particle &p) {
    return decomposition().particle_to_cell(p);
  }
  Cell const *particle_to_cell(const Particle &p) const {
    return decomposition().particle_to_cell(p);
  }

  inline void parallel_for_each_particle_impl(std::span<Cell *const> cells,
                                              ParticleCallback auto &&f) const;

  template <typename RowKernel>
  inline void parallel_for_each_local_row_impl(std::span<Cell *const> cells,
                                               RowKernel &kernel) const;

public:
  /**
   * @brief Add a particle.
   *
   * Moves a particle into the cell system. This adds
   * a particle to the local node, irrespective of where
   * it belongs.
   *
   * @param p Particle to add.
   * @return A view of the particle in the cell system (a by-value
   *         @ref Particle view, valid until the next store rebuild).
   */
  std::optional<Particle> add_particle(Particle &&p);

  /**
   * @brief Add a particle.
   *
   * Moves a particle into the cell system, if it
   * belongs to this node. Otherwise this does not
   * have an effect and the particle is discarded.
   * This can be used to add a particle without
   * knowledge where it should be placed by calling
   * the function on all nodes, it will then add
   * the particle in exactly one place.
   *
   * @param p Particle to add.
   * @return A view of the particle if it is local (a by-value @ref Particle
   *         view), @c std::nullopt otherwise.
   */
  std::optional<Particle> add_local_particle(Particle &&p);

  /**
   * @brief Remove a particle.
   *
   * Removes a particle and all bonds pointing
   * to it. This is a collective call.
   *
   * @param id Id of particle to remove.
   */
  void remove_particle(int id);

  /**
   * @brief Get the maximal particle id on this node.
   *
   * This returns the highest particle id on
   * this node, or -1 if there are no particles on this node.
   */
  int get_max_local_particle_id() const;
  int get_cached_max_local_particle_id() const {
    return m_cached_max_local_particle_id;
  }
  std::size_t get_num_local_particles_cached() const {
    return m_num_local_particles_cached;
  }
  int get_local_pair_bond_numbers() const;
  int get_local_angle_bond_numbers() const;
  int get_local_dihedral_bond_numbers() const;
  void set_local_bond_numbers(int p, int a, int d);
#ifdef ESPRESSO_COLLISION_DETECTION
  void clear_new_bonds();
  void add_new_bond(int bond_id, std::vector<int> const &particle_ids);
  void rebuild_bond_list();
#endif // ESPRESSO_COLLISION_DETECTION

  /**
   * @brief Remove all particles from the cell system.
   *
   * This allows linear time removal of all particles from
   * the system, removing each particle individually would
   * be quadratic.
   */
  void remove_all_particles();

  /**
   * @brief Get the underlying particle decomposition.
   *
   * Should be used solely for informative purposes.
   *
   * @return The active particle decomposition.
   */
  ParticleDecomposition const &decomposition() const {
    return assert(m_decomposition), *m_decomposition;
  }

private:
  ParticleDecomposition &decomposition() {
    return assert(m_decomposition), *m_decomposition;
  }

public:
  /**
   * @brief Increase the local resort level at least to @p level.
   */
  void set_resort_particles(Cells::Resort level) {
    m_resort_particles |= level;
    assert(m_resort_particles >= level);
  }

  /**
   * @brief Get the currently scheduled resort level.
   */
  unsigned get_resort_particles() const { return m_resort_particles; }

  /**
   * @brief Set the resort level to sorted.
   */
  void clear_resort_particles() { m_resort_particles = Cells::RESORT_NONE; }

  /**
   * @brief Check whether a particle has moved further than half the skin
   * since the last Verlet list update, thus requiring a resort.
   * @param additional_offset   Offset which is added to the distance the
   *                            particle has travelled when comparing to half
   *                            the Verlet skin (e.g., for Lees-Edwards BC).
   * @return Whether a resort is needed.
   */
  bool
  check_resort_required(Utils::Vector3d const &additional_offset = {}) const;

  auto get_le_pos_offset_at_last_resort() const {
    return m_le_pos_offset_at_last_resort;
  }

  /**
   * @brief Synchronize number of ghosts.
   */
  void ghosts_count();

  /**
   * @brief Update ghost particles.
   *
   * Update ghost particles with data from the real particles.
   *
   * @param data_parts Particle parts to update, combination of @ref
   * Cells::DataPart
   */
  void ghosts_update(unsigned data_parts);

  /**
   * @brief Update ghost particles, with particle resort if needed.
   *
   * Update ghost particles with data from the real particles.
   * Resort particles if a resort is due.
   *
   * @param data_parts Particle parts to update, combination of @ref
   * Cells::DataPart
   */
  void update_ghosts_and_resort_particle(unsigned data_parts);

  /**
   * @brief Add forces and torques from ghost particles to real particles.
   */
  void ghosts_reduce_forces();

  /**
   * @brief Begin the split-phase ghost force reduction (non-blocking).
   *
   * Posts all MPI sends/receives for the force reduction and returns
   * immediately.  The caller must later call @ref ghosts_reduce_forces_finish.
   * Asserts that no reduction is already in flight.
   *
   * Only meaningful when @c comm_cart.size() > 1.  At one rank the blocking
   * @ref ghosts_reduce_forces must be used instead (the reduction is then a
   * pure local copy with nothing to hide).
   */
  void ghosts_reduce_forces_start();

  /**
   * @brief Complete the split-phase ghost force reduction.
   *
   * Waits for all outstanding MPI requests and unpacks/reduces force data
   * into real particles.  Clears the pending state afterwards.
   *
   * Must be called exactly once after @ref ghosts_reduce_forces_start.
   */
  void ghosts_reduce_forces_finish();

  /** @brief True when a split-phase force reduction is in flight. */
  bool has_pending_ghost_reduce() const {
    return m_pending_ghost_reduce.has_value();
  }

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  /** Add dipole fields from ghost particles to real particles. */
  void ghosts_reduce_dipole_field();

  /** Set dipole fields on all ghosts to zero. */
  void ghosts_reset_dipole_fields() {
    for_each_ghost_particle([](Particle &p) { p.dip_fld() = {}; });
  }
#endif

  /** Set forces and torques on all ghosts to zero. */
  void ghosts_reset_forces() {
    for_each_ghost_particle([](Particle &p) {
      p.force() = {};
#ifdef ESPRESSO_ROTATION
      p.torque() = {};
#endif
    });
  }

#ifdef ESPRESSO_BOND_CONSTRAINT
  /**
   * @brief Add rattle corrections from ghost particles to real particles.
   */
  void ghosts_reduce_rattle_correction();
#endif

  /**
   * @brief Resort particles.
   *
   * @param global_flag  Whether to do a global (all-to-all) resort.
   * @param commit  When true (the default, and the only value used by a
   *   direct/unit-test caller), the just-migrated local particles are
   *   committed to store rows and the id->view index rebuilt before returning,
   *   so the store is clean and the index consistent on return. When false (the
   *   @ref update_ghosts_and_resort_particle hot path), the commit is DEFERRED:
   *   the migrated locals stay STAGED (a cell's @ref Cell::size still counts
   *   them, which is all @ref ghosts_count needs), and a single post-@ref
   *   ghosts_count @ref ensure_particle_store_synchronized then commits locals
   *   AND ghosts in ONE store rebuild. Deferring yields a single store rebuild
   *   per resort instead of two O(N) column copies; no index is read in the
   *   deferred window, so it is safe.
   */
  void resort_particles(bool global_flag, bool commit = true);

  /** @brief Whether the Verlet skin is set. */
  auto is_verlet_skin_set() const { return m_verlet_skin_set; }

  /** @brief Get the Verlet skin. */
  auto get_verlet_skin() const { return m_verlet_skin; }

  /** @brief Set the Verlet skin. */
  void set_verlet_skin(double value);

  /** @brief Set the Verlet skin using a heuristic. */
  void set_verlet_skin_heuristic();

  void update_verlet_stats(int n_steps, int n_verlet_updates) {
    if (n_verlet_updates > 0) {
      m_verlet_reuse = n_steps / static_cast<double>(n_verlet_updates);
    } else {
      m_verlet_reuse = 0.;
    }
  }

  /** @brief Average number of integration steps the Verlet list was re-used */
  auto get_verlet_reuse() const { return m_verlet_reuse; }

  /**
   * @brief Resolve ids to particles.
   *
   * @throws BondResolutionError if one of the ids
   *         was not found.
   *
   * @param partner_ids Ids to resolve.
   * @return Vector of Particle VIEWS held by value.
   *
   * @ref get_local_particle returns a by-value @ref Particle view (a 16-byte
   * handle), so the resolved partners are collected here as by-value views.
   * Callers that need the @c std::span<Particle*> handler contract build a
   * pointer span into this owned buffer (see @c execute_bond_handler), which
   * stays valid for as long as the returned vector lives.
   */
  boost::container::static_vector<Particle, 4>
  resolve_bond_partners(std::span<const int> partner_ids) {
    boost::container::static_vector<Particle, 4> partners;
    for (auto const id : partner_ids) {
      auto view = get_local_particle(id);
      if (not view) {
        throw BondResolutionError{};
      }
      partners.push_back(*view);
    }
    return partners;
  }

private:
  /**
   * @brief Execute kernel for every bond on particle.
   * @tparam Handler Callable, which can be invoked with
   *                 (Particle, int, std::span<Particle *>),
   *                 returning a bool.
   * @param p Particles for whom the bonds are evaluated.
   * @param handler is called for every bond, and handed
   *                p, the bond id and a span with the bond
   *                partners as arguments. Its return value
   *                should indicate if the bond was broken.
   */
  template <class Handler>
  void execute_bond_handler(Particle &p, Handler const &handler) {
    // Debug guard: resolve_bond_partners hands back by-value Particle VIEWS
    // (16-byte handles over store rows); the handler contract is
    // std::span<Particle*>, so we build a pointer span into the owned
    // `partners` buffer -- valid for the duration of the handler call. The
    // views (and thus the rows they alias) are only valid while the store
    // generation is unchanged: a rebuild renumbers the rows. Record the
    // generation on entry and assert it did not move across each handler call
    // (a handler must not mutate topology) to catch a mis-resolved row.
    auto const bond_epoch_generation = m_particle_store.generation();
    for (const BondView bond : p.bonds()) {
      auto const partner_ids = bond.partner_ids();
      try {
        auto partners = resolve_bond_partners(partner_ids);
        boost::container::static_vector<Particle *, 4> partner_ptrs;
        for (auto &partner : partners) {
          partner_ptrs.push_back(std::addressof(partner));
        }
        auto const partners_span =
            std::span(partner_ptrs.data(), partner_ptrs.size());
        auto const bond_broken = handler(p, bond.bond_id(), partners_span);
        ParticleStoreGuard::assert_generation(m_particle_store,
                                              bond_epoch_generation);
        if (bond_broken) {
          bond_broken_error(p.id(), partner_ids);
        }
      } catch (BondResolutionError const &) {
        bond_resolution_error(partner_ids);
      }
    }
  }

  /** @brief Set the particle decomposition, keeping the particles. */
  void set_particle_decomposition(
      std::unique_ptr<ParticleDecomposition> &&decomposition) {
    assert(not m_pending_ghost_reduce.has_value() &&
           "set_particle_decomposition: ghost force reduction is still in "
           "flight — call ghosts_reduce_forces_finish() first");
    clear_particle_index();

    /* Copy every particle out of the OLD decomposition into staging-store rows
     * BEFORE the store is rebuilt: the old cells' rows index into the current
     * main store, which the swap + rebuild below will invalidate. The staging
     * store is an independent store, so its rows survive the main store's
     * rebuild and can seed the retained particles into the new system.
     */
    clear_staging_store();
    std::vector<int> retained_staging_rows;
    for (auto *cell : m_decomposition->local_cells()) {
      cell->set_store(m_particle_store);
      // Stage every LIVE committed row (the CellRowSpan skips pending-removed
      // rows); the new decomposition below re-homes each retained row.
      for (int const live_row : cell->rows()) {
        retained_staging_rows.push_back(stage_row(live_row));
      }
    }

    /* Swap in new cell system */
    std::swap(m_decomposition, decomposition);

    /* Wire the new cells to the store and stage the retained rows into their
     * home cells. Stage directly (NOT via add_particle, which commits per-add
     * -- that would be O(n^2) here); a single rebuild below commits them all at
     * once. A particle with no home cell on this node goes to the first local
     * cell (a global resort will place it), matching add_particle. The home
     * cell is decided from the staged row's position (read from the staging
     * store). */
    for (auto *cell : m_decomposition->local_cells()) {
      cell->set_store(m_particle_store);
    }
    for (auto *cell : m_decomposition->ghost_cells()) {
      cell->set_store(m_particle_store);
    }
    auto const had_retained = not retained_staging_rows.empty();
    for (auto const staging_row : retained_staging_rows) {
      auto const view = m_staging_store.make_view(staging_row);
      // NB: the `decomposition` PARAMETER (now holding the OLD, swapped-out
      // decomposition) shadows the `decomposition()` accessor here; route
      // through the NEW decomposition via m_decomposition explicitly.
      auto const sort_cell = m_decomposition->particle_to_cell(view);
      auto cell = sort_cell ? sort_cell : m_decomposition->local_cells()[0];
      set_resort_particles(sort_cell ? Cells::RESORT_LOCAL
                                     : Cells::RESORT_GLOBAL);
      CellParticleStorage::insert_staged_row(*cell, m_staging_store,
                                             staging_row);
    }

    mark_particle_store_dirty();
    // Commit the retained particles into store rows NOW so they are immediately
    // live (visible to local_particles(), the particle-node bookkeeping, and
    // get_local_particle) after the decomposition swap. Only when there WERE
    // retained particles: an empty initial setup
    // leaves the store column-free (avoids allocating Kokkos columns that a
    // System torn down at static destruction would release post-finalize).
    if (had_retained) {
      ensure_particle_store_synchronized();
    }
  }

public:
  /**
   * @brief Set the particle decomposition to @ref AtomDecomposition.
   */
  void set_atom_decomposition();

  /**
   * @brief Set the particle decomposition to @ref RegularDecomposition.
   *
   * @param range Interaction range.
   * @param fully_connected_boundary neighbor cell directions for Lees-Edwards.
   */
  void set_regular_decomposition(
      double range,
      std::optional<std::pair<int, int>> fully_connected_boundary);

  /**
   * @brief Set the particle decomposition to @ref HybridDecomposition.
   *
   * @param cutoff_regular Interaction cutoff_regular.
   * @param n_square_types Particle types to put into n_square decomposition.
   */
  void set_hybrid_decomposition(double cutoff_regular,
                                std::set<int> n_square_types);

private:
  /**
   * @brief Run link_cell algorithm for local cells.
   *
   * @tparam Kernel Needs to be callable with (Particle, Particle, Distance).
   * @param kernel Pair kernel functor.
   *
   * Iterates the cells' store-ROW bags directly and REBINDS three cached
   * @ref Particle views (p1 + the two partner roles) via
   * @ref Particle::attach_to_store, rather than driving
   * @ref Algorithm::link_cell over @ref RowParticleRange iterators (each such
   * iterator embeds a Particle by value, so the generic algorithm would
   * construct a fresh Particle for every @c std::next(it) copy and every
   * neighbor-range begin()/end() -- tens of thousands of Particle
   * materialisations per Verlet rebuild). The pair ORDER and the p1/p2 role
   * assignment match @ref Algorithm::link_cell (cell self-pairs [i, j>i] then
   * red-neighbour pairs).
   */
  void link_cell(auto kernel) {
    auto const maybe_box = decomposition().minimum_image_distance();
    if (maybe_box) {
      link_cell_rows(detail::MinimalImageDistance{decomposition().box()},
                     kernel);
    } else {
      if (decomposition().box().type() != BoxType::CUBOID) {
        throw std::runtime_error("Non-cuboid box type is not compatible with a "
                                 "particle decomposition that relies on "
                                 "EuclideanDistance for distance calculation.");
      }
      link_cell_rows(detail::EuclidianDistance{}, kernel);
    }
  }

  /** @brief Row-bag link-cell driver reusing three cached views. See
   *  @c link_cell. */
  template <class DistanceFunction, class Kernel>
  void link_cell_rows(DistanceFunction const &df, Kernel &kernel) {
    auto &store = m_particle_store;
    // Three reused views rebound per element: p1 (outer), and the two partner
    // roles (self-cell partner, neighbour partner).
    Particle p1, p_self, p_nb;
    for (auto *cell : decomposition().local_cells()) {
      // Contiguous store-row range. This runs only on a clean store (the caller
      // ran ensure_particle_store_synchronized before any force loop), so there
      // are no pending-removed rows and the raw range IS the live range: index
      // it directly by offset + i (O(1), no skip).
      auto const offset = cell->offset();
      auto const n = cell->count();
      for (std::size_t i = 0u; i < n; ++i) {
        p1.attach_to_store(store, static_cast<int>(offset + i));
        // Pairs within this cell (j > i), same order as Algorithm::link_cell.
        for (std::size_t j = i + 1u; j < n; ++j) {
          p_self.attach_to_store(store, static_cast<int>(offset + j));
          kernel(p1, p_self, df(p1, p_self));
        }
        // Pairs with the red-partition neighbours, same order.
        for (auto *neighbor : cell->neighbors().red()) {
          auto const nb_offset = neighbor->offset();
          auto const nb_n = neighbor->count();
          for (std::size_t k = 0u; k < nb_n; ++k) {
            p_nb.attach_to_store(store, static_cast<int>(nb_offset + k));
            kernel(p1, p_nb, df(p1, p_nb));
          }
        }
      }
    }
  }

public:
  auto get_max_id() const { return m_max_id; }

  void set_kokkos_handle(std::shared_ptr<KokkosHandle> handle);
  void rebuild_local_properties(double pair_cutoff);
  void reset_local_properties();
  /** @brief Zero the local force view and its scatter replicas. */
  void reset_local_force_buffers();

  auto &get_id_to_index() { return *m_id_to_index; }
  auto &get_local_force() { return *m_local_force; }
  auto get_scatter_force() { return *m_scatter_force; }
#ifdef ESPRESSO_ROTATION
  auto &get_local_torque() { return *m_local_torque; }
  auto get_scatter_torque() { return *m_scatter_torque; }
  /** @brief Declare that a kernel scattering into the torque view is about
   *  to run. Must be called before the torque replicas are reduced via the
   *  force loop dispatch.
   */
  void mark_torque_replicas_dirty() { m_torque_replicas_dirty = true; }
  auto torque_replicas_dirty() const { return m_torque_replicas_dirty; }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  auto &get_local_dip_fld() { return *m_local_dip_fld; }
  auto get_scatter_dip_fld() { return *m_scatter_dip_fld; }
  void mark_dip_fld_replicas_dirty() { m_dip_fld_replicas_dirty = true; }
  auto dip_fld_replicas_dirty() const { return m_dip_fld_replicas_dirty; }
#endif
#ifdef ESPRESSO_NPT
  auto &get_local_virial() { return *m_local_virial; }
  auto get_scatter_virial() { return *m_scatter_virial; }
  void mark_virial_replicas_dirty() { m_virial_replicas_dirty = true; }
  auto virial_replicas_dirty() const { return m_virial_replicas_dirty; }
#endif

  auto &get_aosoa() { return *m_aosoa; }
  auto const &get_aosoa() const { return *m_aosoa; }
  auto const &get_unique_particles() const { return m_unique_particles; }
  auto const &get_verlet_list_cabana() const { return *m_verlet_list_cabana; }
  /** @brief Generation of the Cabana Verlet list (bumps on each rebuild). */
  auto verlet_list_cabana_generation() const {
    return m_verlet_list_cabana_generation;
  }
#ifdef ESPRESSO_CUDA
  /** @brief Persistent GPU short-range device buffers handle (opaque). */
  std::shared_ptr<DeviceShortRangeBuffers> &device_sr_buffers() {
    return m_device_sr_buffers;
  }
#endif
  auto &bond_state() { return *m_bond_state; }
  auto const &bond_state() const { return *m_bond_state; }
  void clear_local_properties();
  void clear_bond_properties();

  auto &particle_store() { return m_particle_store; }
  void mark_particle_store_dirty() { m_particle_store.mark_dirty(); }
  /** @brief Rebuild the store row assignment if topology changed.
   *  Purely rank-local; O(1) when the store is clean. */
  void ensure_particle_store_synchronized();

  // -- migration staging store ----------------------------------------------
  /** @brief Direct access to the migration staging store. */
  auto &staging_store() { return m_staging_store; }
  auto const &staging_store() const { return m_staging_store; }
  /** @brief Number of rows currently held in the staging store. */
  int staged_row_count() const { return m_staging_store_next_row; }
  /**
   * @brief Stage a live store row: copy @p live_row of the main store into the
   * next free staging row and return that staging row.
   *
   * Grows the staging store lazily (doubling, preserving already-staged rows)
   * so an arbitrary number of rows can be staged. This is the row-level
   * `extract` helper (@ref ParticleStore::copy_row live -> staging) used by the
   * migration path.
   */
  int stage_row(int live_row);
  /**
   * @brief Reserve @p count fresh, uninitialized staging rows.
   *
   * Grows the staging store (doubling, preserving already-staged rows) so that
   * @p count consecutive rows starting at the returned index are valid store
   * rows, and advances the row counter past them. Used by the migration
   * `receive` path: the decomposition reserves the rows, then
   * @ref MigrationPack::unpack_rows writes the wire buffer into them. Returns
   * the first reserved row index (== the previous @ref staged_row_count).
   */
  int reserve_staging_rows(int count);
  /**
   * @brief Build a fresh, default-seeded new-particle VIEW.
   *
   * Reserves one staging-store row, seeds it to the new-particle defaults
   * (@ref ParticleStore::seed_default_row), and returns a @ref Particle view
   * over it. This is the particle-creation entry point: the caller writes the
   * fields it wants through the returned view, then hands the view to
   * @ref add_particle / @ref add_local_particle, which stage the underlying
   * staging row into the particle's home cell. The view is valid until the next
   * staging-store growth / clear; @ref add_particle consumes it right away.
   */
  Particle make_new_particle_view();
  /** @brief Drop all staged rows (row counter back to zero). The columns are
   *  retained as reusable capacity; @ref release_staging_store frees them. */
  void clear_staging_store() { m_staging_store_next_row = 0; }
  /** @brief Release the staging store's Kokkos columns (teardown / lock-step
   *  with the main store; must run while Kokkos is alive). */
  void release_staging_store() {
    m_staging_store.release_columns();
    m_staging_store_next_row = 0;
    m_staging_store_capacity = 0u;
  }

  /** @brief Build the migration staging handle a decomposition uses to route
   *  migrating particles through this store's staging store.
   *  The @c store pointer is the address-stable staging-store MEMBER (its
   *  internals may swap out on growth, but the member address is fixed), so a
   *  decomposition holding it always packs from the current columns. */
  MigrationStaging make_migration_staging() {
    MigrationStaging staging;
    staging.store = &m_staging_store;
    staging.stage_row = [this](int live_row) { return stage_row(live_row); };
    staging.reserve_rows = [this](int count) {
      return reserve_staging_rows(count);
    };
    staging.clear = [this]() { clear_staging_store(); };
    return staging;
  }

  [[nodiscard]] auto is_verlet_list_cabana_rebuild_needed() const {
    return m_rebuild_verlet_list_cabana;
  }

  /**
   * @brief Update bond storage(m_*_bond_list_kokkos and m_*_bond_id_kokkos).
   * @param pair_count      Index for pair bond storage.
   * @param angle_count     Index for angle bond storage.
   * @param dihedral_count  Index for dihedral bond storage.
   * @param p               Particle pointer.
   */
  void update_bond_storage(int &pair_count, int &angle_count,
                           int &dihedral_count, Particle const &p);

  /**
   * @brief Reset local properties of the Verlet list.
   * @param cutoff    Pair interaction cutoff.
   * @return True if a rebuild is needed.
   */
  [[nodiscard]] auto prepare_verlet_list_cabana(double cutoff) {
    auto const rebuild = is_verlet_list_cabana_rebuild_needed();
    if (rebuild) {
      // If we have to rebuild, we need to count the particles
      set_index_map(); // parallelized index_map
      // Create essential variables for MD
      rebuild_local_properties(cutoff);
    } else {
      // If we do not rebuild we can use the saved map
      reset_local_properties();
    }
    return rebuild;
  }

  void rebuild_verlet_list_cabana(auto &&kernel, bool rebuild_verlet_list) {
    assert(is_verlet_list_cabana_rebuild_needed());
    if (rebuild_verlet_list) {
      kernel(m_decomposition->local_cells(), m_decomposition->box(),
             *m_verlet_list_cabana);
      ++m_verlet_list_cabana_generation;
    }
    m_rebuild_verlet_list_cabana = false;
  }

  void set_index_map();

  /** @brief Pack-index -> store-row translation view.
   *  @c view(i)==i for all @c i < count_local_particles(). */
  auto pack_index_to_store_row() const { return m_pack_index_to_store_row; }
  /** @brief Host view of the ParticleStore position column. */
  auto store_position_view() { return m_particle_store.position_view(); }
  /** @brief Host view of the ParticleStore image-box column. */
  auto store_image_view() { return m_particle_store.image_box_view(); }
  /** @brief Host view of the ParticleStore velocity column. */
  auto store_velocity_view() { return m_particle_store.velocity_view(); }
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  /** @brief Store-side derived director view. */
  auto director_view() const { return m_director_view; }
  /** @brief Recompute the derived director view from the store quaternion
   *  column. Called every force-calc commit step. */
  void update_director_view();
#endif

  /** @brief Point the pack's store-aliased views at the current store columns
   *  and the translation/director views. Cheap; must run after
   *  @ref set_index_map (or a store rebuild) before the pack is used. */
  void bind_pack_store_views();

  inline void cell_list_loop(auto &&kernel) {
    kernel(m_decomposition->local_cells(), m_decomposition->box());
  }

private:
  /** @brief Zero the torque buffers iff a kernel scattered into them since
   *  the last reset (no-op otherwise, and without the ROTATION feature).
   */
  void reset_torque_replicas_if_dirty();
  /** @brief Same contract for the dipolar fields buffers. */
  void reset_dip_fld_replicas_if_dirty();
  /** @brief Same contract for the virial buffers (NPT feature). */
  void reset_virial_replicas_if_dirty();

  /** Non-bonded pair loop with verlet lists.
   *
   * @param pair_kernel Kernel to apply
   * @param verlet_criterion Filter for verlet lists.
   */
  template <class PairKernel, class VerletCriterion>
  void verlet_list_loop(PairKernel pair_kernel,
                        const VerletCriterion &verlet_criterion) {
    /* In this case the verlet list update is attached to
     * the pair kernel, and the verlet list is rebuilt as
     * we go. */
    if (m_rebuild_verlet_list) {
      m_verlet_list.clear();

      // Record each pair as store ROW indices, not Particle pointers: cells do
      // not own stable Particle addresses. The particles handed out by
      // link_cell are store-attached (the caller ran
      // ensure_particle_store_synchronized before any force loop), so
      // store_row() is valid. Rows are only meaningful for the store
      // generation they were recorded at -- stamp it so the consume branch can
      // detect (debug) a rebuild that renumbered them without a Verlet rebuild.
      link_cell([&](Particle &p1, Particle &p2, Distance const &d) {
        if (verlet_criterion(p1, p2, d.dist2)) {
          m_verlet_list.emplace_back(p1.store_row(), p2.store_row());
          pair_kernel(p1, p2, d);
        }
      });

      m_verlet_list_store = &m_particle_store;
      m_verlet_list_generation = m_particle_store.generation();
      m_rebuild_verlet_list = false;
      m_rebuild_verlet_list_cabana = true;
    } else {
      // Debug guard: the stored rows are only valid while the store generation
      // is unchanged. A rebuild between build and consume (without a Verlet
      // rebuild) would alias the wrong particle -- fire here in debug. The
      // production invariant is that every generation bump sets
      // m_rebuild_verlet_list.
      assert(m_verlet_list_store == &m_particle_store);
      ParticleStoreGuard::assert_generation(
          m_particle_store, m_verlet_list_generation, m_verlet_list_store);
      // Hoist the store pointer once; resolve each row -> view at loop entry,
      // in the SAME pair order (p1 then p2) as the old pointer derefs, so the
      // per-pair arithmetic order is byte-for-byte unchanged (identity gate).
      auto &store = m_particle_store;
      auto const maybe_box = decomposition().minimum_image_distance();
      // Reuse two cached views across the whole pair loop, REBOUND to each
      // pair's rows via attach_to_store (two handle-field writes), instead of
      // constructing a fresh Particle per row per pair. Building views is two
      // full Particle materialisations per pair -- hot, since collision
      // detection walks this branch every step.
      Particle p1, p2;
      /* In this case the pair kernel is just run over the verlet list. */
      if (maybe_box) {
        auto const distance_function =
            detail::MinimalImageDistance{decomposition().box()};
        for (auto const &[row1, row2] : m_verlet_list) {
          p1.attach_to_store(store, row1);
          p2.attach_to_store(store, row2);
          pair_kernel(p1, p2, distance_function(p1, p2));
        }
      } else {
        auto const distance_function = detail::EuclidianDistance{};
        for (auto const &[row1, row2] : m_verlet_list) {
          p1.attach_to_store(store, row1);
          p2.attach_to_store(store, row2);
          pair_kernel(p1, p2, distance_function(p1, p2));
        }
      }
    }
  }

public:
  /** Bonded pair loop.
   * @param bond_kernel Kernel to apply
   */
  template <class BondKernel> void bond_loop(BondKernel const &bond_kernel) {
    for (auto &p : local_particles()) {
      execute_bond_handler(p, bond_kernel);
    }
  }

  /** Non-bonded pair loop.
   * @param pair_kernel Kernel to apply
   */
  template <class PairKernel> void non_bonded_loop(PairKernel pair_kernel) {
    link_cell(pair_kernel);
  }

  /** Non-bonded pair loop with potential use
   * of verlet lists.
   * @param pair_kernel Kernel to apply
   * @param verlet_criterion Filter for verlet lists.
   */
  template <class PairKernel, class VerletCriterion>
  void non_bonded_loop(PairKernel pair_kernel,
                       const VerletCriterion &verlet_criterion) {
    if (use_verlet_list) {
      verlet_list_loop(pair_kernel, verlet_criterion);
    } else {
      /* No verlet lists, just run the kernel with pairs from the cells. */
      link_cell(pair_kernel);
    }
  }

  /**
   * @brief Check that particle index is commensurate with particles.
   *
   * For each local particles is checked that has a correct entry
   * in the particles index, and that there are no excess (non-existing)
   * particles in the index.
   */
  void check_particle_index() const;

  /**
   * @brief Check that particles are in the correct cell.
   *
   * This checks for all local particles that the result
   * of particles_to_cell is the cell the particles is
   * actually in, e.g. that the particles are sorted according
   * to particles_to_cell.
   */
  void check_particle_sorting() const;

public:
  /**
   * @brief Find cell a particle is stored in.
   *
   * For local particles, this returns the cell they
   * are stored in, otherwise nullptr is returned.
   *
   * @param p Particle to find cell for
   * @return Cell for particle or nullptr.
   */
  Cell *find_current_cell(const Particle &p) {
    assert(not get_resort_particles());

    if (p.is_ghost()) {
      return nullptr;
    }

    return particle_to_cell(p);
  }

  /**
   * @brief Run kernel on all particles inside local cell and its neighbors.
   *
   * @param p      Particle to find cell for
   * @param kernel Function with signature <tt>double(Particle const&,
   *               Particle const&, Utils::Vector3d const&)</tt>
   * @return false if cell is not found, otherwise true
   */
  template <class Kernel>
  bool run_on_particle_short_range_neighbors(Particle const &p,
                                             Kernel &kernel) {
    auto const cell = find_current_cell(p);

    if (cell == nullptr) {
      return false;
    }

    auto const maybe_box = decomposition().minimum_image_distance();

    if (maybe_box) {
      auto const distance_function =
          detail::MinimalImageDistance{decomposition().box()};
      short_range_neighbor_loop(p, cell, kernel, distance_function);
    } else {
      auto const distance_function = detail::EuclidianDistance{};
      short_range_neighbor_loop(p, cell, kernel, distance_function);
    }
    return true;
  }

private:
  template <class Kernel, class DistanceFunc>
  void short_range_neighbor_loop(Particle const &p1, Cell *const cell,
                                 Kernel &kernel, DistanceFunc const &df) {
    /* Iterate over particles inside cell */
    for (auto const &p2 : cell->particles()) {
      if (p1.id() != p2.id()) {
        auto const vec = df(p1, p2).vec21;
        kernel(p1, p2, vec);
      }
    }
    /* Iterate over all neighbors */
    for (auto const neighbor : cell->neighbors().all()) {
      /* Iterate over particles in neighbors */
      if (neighbor != cell) {
        for (auto const &p2 : neighbor->particles()) {
          auto const vec = df(p1, p2).vec21;
          kernel(p1, p2, vec);
        }
      }
    }
  }
};

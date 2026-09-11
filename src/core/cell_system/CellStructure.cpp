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

#include "cell_system/CellStructure.hpp"

#include "cell_system/AtomDecomposition.hpp"
#include "cell_system/HybridDecomposition.hpp"
#include "cell_system/ParticleDecomposition.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "cell_system/RegularDecomposition.hpp"

#include "BoxGeometry.hpp"
#include "LocalBondState.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "aosoa_pack.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloExchange.hpp"
#include "integrators/Propagation.hpp"
#include "kokkos_helpers.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "particle_enumeration.hpp"
#include "particle_node.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/contains.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/sqr.hpp>
#include <utils/quaternion.hpp>

#ifdef ESPRESSO_CALIPER
#include "caliper_utils.hpp"
#endif

#include <boost/mpi/collectives/all_reduce.hpp>

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

CellStructure::~CellStructure() {
  assert(not m_pending_ghost_reduce.has_value() &&
         "~CellStructure: ghost force reduction still in flight at destruction "
         "— ghosts_reduce_forces_finish() was not called");
  clear_local_properties();
  // Release the ParticleStore columns while Kokkos is still alive; the handle
  // reset below may finalize Kokkos, and leaving the store's Views to
  // member-destruction would free them after Kokkos::finalize.
  m_particle_store.release_columns();
  // Same Kokkos-lifetime constraint for the migration staging store.
  m_staging_store.release_columns();
  // Kokkos handle can only be freed after all Cabana containers have been freed
  m_kokkos_handle.reset();
}

void CellStructure::clear_local_properties() {
  m_scatter_force.reset();
  m_local_force.reset();
#ifdef ESPRESSO_ROTATION
  m_scatter_torque.reset();
  m_local_torque.reset();
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  m_scatter_dip_fld.reset();
  m_local_dip_fld.reset();
#endif
#ifdef ESPRESSO_NPT
  m_scatter_virial.reset();
  m_local_virial.reset();
#endif
  m_id_to_index.reset();
  m_aosoa.reset();
  m_verlet_list_cabana.reset();
  m_bond_state->clear();
  // Release the Kokkos-backed views owned here (translation view and derived
  // director). The pack's aliasing views were just dropped by
  // m_aosoa.reset() above, so releasing the owners now is safe and must happen
  // while Kokkos is still alive (see the destructor's release ordering).
  m_pack_index_to_store_row = Kokkos::View<int *, Kokkos::HostSpace>{};
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  m_director_view = Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                                 Kokkos::HostSpace>{};
#endif
  m_rebuild_verlet_list_cabana = true;
}
void CellStructure::clear_bond_properties() { m_bond_state->reset(); }

void CellStructure::set_kokkos_handle(std::shared_ptr<KokkosHandle> handle) {
  m_kokkos_handle = std::move(handle);
  m_bond_state = std::make_unique<LocalBondState>();
}

static auto estimate_max_counts(double pair_cutoff,
                                std::size_t number_of_unique_particles,
                                double local_box_volume,
                                std::size_t num_local_particles) {
  if (std::isinf(pair_cutoff)) {
    return number_of_unique_particles;
  }
  if (pair_cutoff < 0.) {
    pair_cutoff = 0.;
  }
  // Estimate number of neighbors based on local density and cutoff sphere:
  // volume n_neighbors = rho * (4/3) * pi * r^3, where rho = n_particles /
  // volume
  auto const local_density =
      (local_box_volume > 0. && num_local_particles > 0)
          ? static_cast<double>(num_local_particles) / local_box_volume
          : 0.;
  auto const cutoff_sphere_volume =
      (4. / 3.) * std::numbers::pi * Utils::int_pow<3>(pair_cutoff);
  // account for local fluctuations. Empirical.
  auto const fluctuation_factor = 2.;
  auto max_counts = static_cast<std::size_t>(
      std::ceil(fluctuation_factor * local_density * cutoff_sphere_volume));
  std::size_t constexpr threshold_num = 16;
  if (max_counts < threshold_num) {
    max_counts = std::min(threshold_num, number_of_unique_particles);
  }
  return max_counts;
}

void CellStructure::rebuild_local_properties(double const pair_cutoff) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  assert(m_kokkos_handle);
  auto const num_part = get_unique_particles().size();
  auto const &system = get_system();
  auto const local_box_volume = system.local_geo->volume();
  auto max_counts = estimate_max_counts(pair_cutoff, num_part, local_box_volume,
                                        get_num_local_particles_cached());
#ifdef ESPRESSO_COLLISION_DETECTION
  if (system.has_collision_detection_enabled()) {
    // TODO: use other types of Verlet list data structures
    max_counts = num_part * 2ul;
  }
#endif
  if (m_local_force) { // local properties are reallocated
    if (get_local_force().extent(0) == num_part) {
      // Extents unchanged (always the case with a single MPI rank): zero the
      // existing buffers in place instead of freeing and reallocating the
      // O(n_threads * N) ScatterView scratch on every Verlet rebuild.
      reset_local_force_buffers();
      reset_torque_replicas_if_dirty();
      reset_dip_fld_replicas_if_dirty();
    } else {
      Kokkos::realloc(get_local_force(), num_part);
      // underlying View extent changed -> scratch buffers must be rebuilt
      m_scatter_force.emplace(
          Kokkos::Experimental::create_scatter_view(get_local_force()));
#ifdef ESPRESSO_ROTATION
      Kokkos::realloc(get_local_torque(), num_part);
      // underlying View extent changed -> scratch buffers must be rebuilt
      m_scatter_torque.emplace(
          Kokkos::Experimental::create_scatter_view(get_local_torque()));
      m_torque_replicas_dirty = false;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
      Kokkos::realloc(get_local_dip_fld(), num_part);
      // underlying View extent changed -> scratch buffers must be rebuilt
      m_scatter_dip_fld.emplace(
          Kokkos::Experimental::create_scatter_view(get_local_dip_fld()));
      m_dip_fld_replicas_dirty = false;
#endif
    }
    auto const required_index_size = get_cached_max_local_particle_id() + 1;
    if (get_id_to_index().extent(0) !=
        static_cast<std::size_t>(required_index_size)) {
      Kokkos::realloc(Kokkos::WithoutInitializing, get_id_to_index(),
                      required_index_size);
    }
    kokkos_deep_copy(execution_space{}, get_id_to_index(), -1);
    // Resize particle views using AoSoA_pack's resize method
    m_aosoa->resize(num_part);
    kokkos_deep_copy(execution_space{}, m_aosoa->flags, uint8_t{0});
    m_verlet_list_cabana->reallocData(num_part, max_counts);
  } else { // local properties are initialized
    m_local_force = std::make_unique<ForceType>("local_force", num_part);
    m_scatter_force.emplace(
        Kokkos::Experimental::create_scatter_view(*m_local_force));
#ifdef ESPRESSO_ROTATION
    m_local_torque = std::make_unique<ForceType>("local_torque", num_part);
    m_scatter_torque.emplace(
        Kokkos::Experimental::create_scatter_view(*m_local_torque));
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    m_local_dip_fld = std::make_unique<ForceType>("local_dip_fld", num_part);
    m_scatter_dip_fld.emplace(
        Kokkos::Experimental::create_scatter_view(*m_local_dip_fld));
#endif
    m_id_to_index = std::make_unique<Kokkos::View<int *, memory_space>>(
        Kokkos::view_alloc(execution_space{}, Kokkos::WithoutInitializing,
                           "id_to_index"),
        get_cached_max_local_particle_id() + 1);
    kokkos_deep_copy(execution_space{}, get_id_to_index(), -1);
    // Create AoSoA_pack and initialize with resize
    m_aosoa = std::make_unique<AoSoA_pack>();
    m_aosoa->resize(num_part);
    kokkos_deep_copy(execution_space{}, m_aosoa->flags, uint8_t{0});

    m_verlet_list_cabana =
        std::make_unique<ListType>(0ul, num_part, max_counts);
  }
#ifdef ESPRESSO_NPT
  if (not m_local_virial) {
    m_local_virial = std::make_unique<VirialType>("local_virial");
    m_scatter_virial.emplace(
        Kokkos::Experimental::create_scatter_view(*m_local_virial));
  } else {
    reset_virial_replicas_if_dirty();
  }
#endif
}

void CellStructure::reset_local_force_buffers() {
  kokkos_deep_copy(execution_space{}, get_local_force(), 0.);
  m_scatter_force->reset();
}

void CellStructure::reset_torque_replicas_if_dirty() {
#ifdef ESPRESSO_ROTATION
  if (m_torque_replicas_dirty) {
    kokkos_deep_copy(execution_space{}, get_local_torque(), 0.);
    m_scatter_torque->reset();
    m_torque_replicas_dirty = false;
  }
#endif
}

void CellStructure::reset_dip_fld_replicas_if_dirty() {
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (m_dip_fld_replicas_dirty) {
    kokkos_deep_copy(execution_space{}, get_local_dip_fld(), 0.);
    m_scatter_dip_fld->reset();
    m_dip_fld_replicas_dirty = false;
  }
#endif
}

void CellStructure::reset_virial_replicas_if_dirty() {
#ifdef ESPRESSO_NPT
  if (m_virial_replicas_dirty) {
    kokkos_deep_copy(execution_space{}, get_local_virial(), 0.);
    m_scatter_virial->reset();
    m_virial_replicas_dirty = false;
  }
#endif
}

void CellStructure::reset_local_properties() {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  reset_local_force_buffers();
  reset_torque_replicas_if_dirty();
  reset_dip_fld_replicas_if_dirty();
  reset_virial_replicas_if_dirty();
  // The has-exclusion `flags` column is NOT zeroed here. This runs on every
  // partial-update (no-rebuild) step; the flag is rebuild-cadence data
  // (commit_particle writes it only on rebuild, and exclusion changes force a
  // rebuild), so it is valid across partial steps and must be preserved. It is
  // (re)zeroed and rewritten on the rebuild path (rebuild_local_properties +
  // the full-commit loop). Zeroing it here would drop live exclusion flags for
  // the rest of the step.
}

void CellStructure::update_bond_storage(int &pair_count, int &angle_count,
                                        int &dihedral_count,
                                        Particle const &p) {
  auto &pair_list = m_bond_state->pair_list;
  auto &pair_ids = m_bond_state->pair_ids;
  auto &angle_list = m_bond_state->angle_list;
  auto &angle_ids = m_bond_state->angle_ids;
  auto &dihedral_list = m_bond_state->dihedral_list;
  auto &dihedral_ids = m_bond_state->dihedral_ids;
  for (auto const bond : p.bonds()) {
    auto const partner_ids = bond.partner_ids();
    try {
      auto const partners = resolve_bond_partners(partner_ids);
      // resolve_bond_partners yields by-value Particle views.
      if (partners.size() == 1u) { // pair bonds
        auto p_index = Kokkos::atomic_fetch_add(&pair_count, 1);
        pair_list(p_index, 0) = p.id();
        pair_list(p_index, 1) = partners[0].id();
        pair_ids(p_index) = bond.bond_id();
      } else if (partners.size() == 2u) { // angle bond
        auto a_index = Kokkos::atomic_fetch_add(&angle_count, 1);
        angle_list(a_index, 0) = p.id();
        angle_list(a_index, 1) = partners[0].id();
        angle_list(a_index, 2) = partners[1].id();
        angle_ids(a_index) = bond.bond_id();
      } else if (partners.size() == 3u) { // dihedral bond
        auto d_index = Kokkos::atomic_fetch_add(&dihedral_count, 1);
        dihedral_list(d_index, 0) = p.id();
        dihedral_list(d_index, 1) = partners[0].id();
        dihedral_list(d_index, 2) = partners[1].id();
        dihedral_list(d_index, 3) = partners[2].id();
        dihedral_ids(d_index) = bond.bond_id();
      }
    } catch (BondResolutionError const &) {
      bond_resolution_error(partner_ids);
    }
  }
}

void CellStructure::set_index_map() {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  // Rebuild the pack-order participant list (m_unique_particles) directly from
  // the synchronized store. The pack order is the local prefix (store rows
  // [0, n_local) in row order == pack order, so pack index i == store row i on
  // the local prefix) followed by the deduped ghost tail (first occurrence of
  // each ghost id not owned by a local, in store-row order). The by-value views
  // are materialized into m_unique_particle_views (the owning backing) and
  // m_unique_particles points into it. The dedup reuses the id->store-row map
  // the store rebuild just built (m_id_to_store_row): "locals win, then first
  // valid ghost row wins" is precisely what that map resolves each id to.
  assert(not m_particle_store.is_dirty());
  auto &unique_particles = m_unique_particles;
  auto &unique_particle_views = m_unique_particle_views;
  unique_particles.clear();
  unique_particle_views.clear();
  auto const n_local = m_particle_store.number_of_local_particles();
  auto const n_total = m_particle_store.number_of_particles();
  unique_particle_views.reserve(n_total);
  int max_id = -1;
  int pair_count = 0;
  int angle_count = 0;
  int dihedral_count = 0;
  m_bond_state->reset_counts();
  // Local prefix: every local row participates, in row order.
  for (std::size_t r = 0u; r < n_local; ++r) {
    auto const &p = unique_particle_views.emplace_back(
        m_particle_store.make_view(static_cast<int>(r)));
    max_id = std::max(p.id(), max_id);
    // Count bonds only on the local prefix (never walk ghosts).
    for (auto const bond : p.bonds()) {
      auto const partner_ids = bond.partner_ids();
      if (partner_ids.empty()) {
        continue;
      }
      if (partner_ids.size() == 1u) {
        ++pair_count;
      } else if (partner_ids.size() == 2u) {
        ++angle_count;
      } else if (partner_ids.size() == 3u) {
        ++dihedral_count;
      }
    }
  }
  // Ghost tail: the deduped ghosts, in store-row order. A ghost row
  // participates iff its id resolves (in m_id_to_store_row) back to THIS row --
  // i.e. it is the winning (first, not-shadowed-by-a-local) copy of that id.
  for (std::size_t r = n_local; r < n_total; ++r) {
    auto const id = m_particle_store.id(static_cast<int>(r));
    if (id < 0) {
      continue;
    }
    auto const resolved =
        (static_cast<unsigned int>(id) < m_id_to_store_row.size())
            ? m_id_to_store_row[static_cast<unsigned int>(id)]
            : no_store_row;
    if (resolved != static_cast<int>(r)) {
      continue; // shadowed by a local or by an earlier ghost of the same id
    }
    auto const &p = unique_particle_views.emplace_back(
        m_particle_store.make_view(static_cast<int>(r)));
    max_id = std::max(p.id(), max_id);
  }
  // The backing buffer is now final-sized: take stable pointers into it.
  unique_particles.reserve(unique_particle_views.size());
  for (auto &p : unique_particle_views) {
    unique_particles.emplace_back(std::addressof(p));
  }
  set_local_bond_numbers(pair_count, angle_count, dihedral_count);
  m_bond_state->allocate();
  m_cached_max_local_particle_id = max_id;
  m_num_local_particles_cached = unique_particles.size();

  // Build the pack-index -> store-row translation view. The store must already
  // be synchronized (rows assigned) by the caller. The local prefix is the
  // identity: the pack order matches ensure_particle_store_synchronized's local
  // loop, so pack index i equals store row i for every local particle. Only the
  // deduped ghost tail needs a real translation.
  assert(not m_particle_store.is_dirty());
  auto const n_unique = unique_particles.size();
  if (m_pack_index_to_store_row.extent(0) != n_unique) {
    m_pack_index_to_store_row = Kokkos::View<int *, Kokkos::HostSpace>(
        Kokkos::ViewAllocateWithoutInitializing("pack_index_to_store_row"),
        n_unique);
  }
  for (std::size_t i = 0u; i < n_unique; ++i) {
    m_pack_index_to_store_row(i) = unique_particles[i]->store_row();
    // local prefix must be the identity (see comment above)
    assert(i >= n_local or m_pack_index_to_store_row(i) == static_cast<int>(i));
  }
  static_cast<void>(n_local);
}

void CellStructure::bind_pack_store_views() {
  auto &aosoa = *m_aosoa;
  aosoa.position = m_particle_store.position_view();
  aosoa.image = m_particle_store.image_box_view();
  aosoa.velocity = m_particle_store.velocity_view();
  // id/mass/charge alias the authoritative ParticleStore scalar columns
  // (indexed by store row), like position/velocity, and are read only on cold
  // paths (bond kernels; charge is needed by BondedCoulomb even with no coulomb
  // solver, so it must always be valid). type is a pack-owned contiguous column
  // written on rebuild in commit_particle, and pair_charge/pair_dipm are
  // pack-owned contiguous hot-path columns refreshed per step by
  // refresh_pack_charges/refresh_pack_dipm when the corresponding solver is
  // active; those pack-owned columns are not bound here.
  aosoa.id = m_particle_store.id_view();
#ifdef ESPRESSO_ELECTROSTATICS
  aosoa.charge = m_particle_store.q_view();
#endif
#ifdef ESPRESSO_MASS
  aosoa.mass = m_particle_store.mass_view();
#endif
  aosoa.row_map = m_pack_index_to_store_row;
#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
  aosoa.director = m_director_view;
#endif
}

#if defined(ESPRESSO_GAY_BERNE) or defined(ESPRESSO_DIPOLES)
void CellStructure::update_director_view() {
#ifdef ESPRESSO_CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  auto const n_total = m_particle_store.number_of_particles();
  if (m_director_view.extent(0) != n_total) {
    m_director_view =
        Kokkos::View<double *[3], ParticleStore::StateVectorLayout,
                     Kokkos::HostSpace>(
            Kokkos::ViewAllocateWithoutInitializing("director"), n_total);
  }
  auto quaternion = m_particle_store.quaternion_view();
  auto director = m_director_view;
  Kokkos::parallel_for("update_director_view",
                       Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(
                           std::size_t{0}, n_total),
                       [quaternion, director](std::size_t row) {
                         Utils::Quaternion<double> q;
                         q[0] = quaternion(row, 0);
                         q[1] = quaternion(row, 1);
                         q[2] = quaternion(row, 2);
                         q[3] = quaternion(row, 3);
                         auto const d =
                             Utils::convert_quaternion_to_director(q);
                         director(row, 0) = d[0];
                         director(row, 1) = d[1];
                         director(row, 2) = d[2];
                       });
  Kokkos::fence();
}
#endif

CellStructure::CellStructure(BoxGeometry const &box)
    : m_decomposition{std::make_unique<AtomDecomposition>(box)} {
  mark_particle_store_dirty();
}

// Permutation builder. Walks cells in the rebuild order and emits perm[] (old
// row per surviving new row, -1 for staged/ghost) + the future (offset, count)
// per cell. Rebuild order: local cells in span order, then ghost cells; per
// cell surviving rows in RAW range order, then staged in push order. A row
// dropped mid-epoch is marked pending-removed on the store; it is NOT carried
// into the permutation (the row is dropped by this rebuild), so surviving rows
// keep their store-row order. Must be called BEFORE the permute rebuild swaps
// the generation out (it reads the current-generation pending-removal mask and
// the cells' current ranges).
void CellStructure::build_resort_permutation(
    std::vector<int> &permutation,
    std::vector<std::pair<std::size_t, std::size_t>> &cell_ranges) const {
  permutation.clear();
  cell_ranges.clear();
  auto const &store = m_particle_store;
  // The permutation entry for a new row is the OLD store row it survives from
  // (>= 0), or -1 for a staged / fresh-ghost row that carries no old-generation
  // data. permutation.size() walks the future contiguous store; each cell's
  // future range is [offset, offset + count).
  auto emit_cells = [&permutation, &cell_ranges,
                     &store](std::span<Cell *const> cells) {
    for (auto const *cell : cells) {
      auto const offset = permutation.size();
      // Surviving committed rows: preserve by old row, in raw range order,
      // skipping any dropped (pending-removed) rows.
      auto const raw_offset = cell->offset();
      auto const raw_count = cell->count();
      for (std::size_t k = 0u; k < raw_count; ++k) {
        auto const old_row = static_cast<int>(raw_offset + k);
        if (store.is_pending_removal(old_row)) {
          continue;
        }
        permutation.push_back(old_row);
      }
      // Staged rows (new / migrated / fresh ghost) in push order: no old-
      // generation row, so a sentinel -1 (the permute rebuild seeds defaults;
      // the caller copies a staged local's data into the row via copy_row).
      for (std::size_t k = 0u; k < cell->staged().size(); ++k) {
        permutation.push_back(-1);
      }
      auto const count = permutation.size() - offset;
      cell_ranges.emplace_back(offset, count);
    }
  };
  emit_cells(decomposition().local_cells());
  emit_cells(decomposition().ghost_cells());
}

void CellStructure::ensure_particle_store_synchronized() {
  if (not m_particle_store.is_dirty()) {
    return;
  }
  // Wire every cell to the store (cheap; idempotent) so that Cell::particles()
  // and Cell::store() work throughout the rebuild and the following step.
  auto wire_stores = [this](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      cell->set_store(m_particle_store);
    }
  };
  wire_stores(decomposition().local_cells());
  wire_stores(decomposition().ghost_cells());

  // A cell's content = its surviving committed rows (the raw range minus
  // pending-removed rows) + its staged particles. Count both (the live count)
  // so the permute rebuild sizes the columns for the new generation.
  // rows().size() is the LIVE count (excludes pending-removed).
  auto cell_count = [](Cell const *cell) {
    return cell->rows().size() + cell->staged().size();
  };
  std::size_t n_local = 0u;
  for (auto const *cell : decomposition().local_cells()) {
    n_local += cell_count(cell);
  }
  std::size_t n_ghost = 0u;
  for (auto const *cell : decomposition().ghost_cells()) {
    n_ghost += cell_count(cell);
  }

  // Resort as a pure column permutation:
  //   1. build the permutation (old row per surviving new row, -1 for staged /
  //      fresh-ghost) + the future (offset, count) per cell, walking cells in
  //      the rebuild order (locals then ghosts; per cell surviving rows in raw
  //      range order skipping dropped rows, then staged in push order);
  //   2. permute_rebuild moves every surviving column/sidecar in one pass per
  //      column (contiguous, vectorizable) and seeds the -1 rows to the
  //      new-particle defaults;
  //   3. copy each staged LOCAL's source-store row into its (already
  //      default-seeded) new row via copy_row (fresh-default ghosts stay at
  //      defaults for the ghost exchange to fill);
  //   4. write the (offset, count) back onto each Cell and clear its staging.
  // The order contract is surviving-then-staged per cell, cells in span order.
  std::vector<int> permutation;
  std::vector<std::pair<std::size_t, std::size_t>> cell_ranges;
  build_resort_permutation(permutation, cell_ranges);
  assert(permutation.size() == n_local + n_ghost);

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Self-check of the permuted result: capture the id each new row is EXPECTED
  // to carry, sourced exactly as the permute path sources it -- a survivor
  // keeps the current-generation id at its old row; a staged row takes its
  // source-store id (or -1 for a fresh, source-less ghost). Captured NOW,
  // before permute_rebuild swaps the generation out; verified row-for-row
  // against the store AFTER the rebuild (below). This is the
  // both-paths-in-debug cross-verification flipped to a self-check of the
  // permuted result.
  std::vector<int> debug_expected_ids;
  debug_expected_ids.reserve(permutation.size());
  {
    auto const &store = m_particle_store;
    auto capture_expected =
        [&store, &debug_expected_ids](std::span<Cell *const> cells) {
          for (auto const *cell : cells) {
            auto const raw_offset = cell->offset();
            auto const raw_count = cell->count();
            for (std::size_t k = 0u; k < raw_count; ++k) {
              auto const old_row = static_cast<int>(raw_offset + k);
              if (store.is_pending_removal(old_row)) {
                continue;
              }
              debug_expected_ids.push_back(store.id(old_row));
            }
            for (auto const &staged : cell->staged()) {
              debug_expected_ids.push_back(
                  staged.source_store != nullptr
                      ? staged.source_store->id(staged.source_row)
                      : -1);
            }
          }
        };
    capture_expected(decomposition().local_cells());
    capture_expected(decomposition().ghost_cells());
  }
#endif

  // Capture each staged entry's source (store + row) in the SAME walk order the
  // builder used, paired with the future row it lands in, so the copy_row step
  // can overwrite the default-seeded staged-local rows after the permute. The
  // walk mirrors build_resort_permutation: per cell, surviving (live) rows fill
  // the leading part of the cell's range, then staged rows fill the tail.
  struct StagedCopy {
    ParticleStore *source_store;
    int source_row;
    int destination_row;
  };
  std::vector<StagedCopy> staged_copies;
  {
    std::size_t new_row = 0u;
    auto collect_staged = [this, &new_row,
                           &staged_copies](std::span<Cell *const> cells) {
      for (auto *cell : cells) {
        // Skip the surviving-row leading part (their perm entries are >= 0).
        auto const raw_offset = cell->offset();
        auto const raw_count = cell->count();
        for (std::size_t k = 0u; k < raw_count; ++k) {
          if (not m_particle_store.is_pending_removal(
                  static_cast<int>(raw_offset + k))) {
            ++new_row;
          }
        }
        for (auto const &staged : cell->staged()) {
          if (staged.source_store != nullptr) {
            staged_copies.push_back(StagedCopy{staged.source_store,
                                               staged.source_row,
                                               static_cast<int>(new_row)});
          }
          ++new_row;
        }
      }
    };
    collect_staged(decomposition().local_cells());
    collect_staged(decomposition().ghost_cells());
    assert(new_row == permutation.size());
  }

  // Step 2: permute every surviving column from the retired generation and seed
  // the -1 (staged / fresh-ghost) rows to defaults, in one pass per column.
  m_particle_store.permute_rebuild(permutation, n_local, n_ghost);

  // Step 3: copy each staged LOCAL's source row into its new (default-seeded)
  // row. The source (staging) store rows stay valid until clear_staging_store
  // below.
  for (auto const &copy : staged_copies) {
    m_particle_store.copy_row(*copy.source_store, copy.source_row,
                              copy.destination_row);
  }

  // Step 4: write the (offset, count) collapse back onto each cell and clear
  // its staging area. The ranges tile [0, n_total) contiguously in cell order.
  {
    std::size_t range_index = 0u;
    auto writeback = [&cell_ranges,
                      &range_index](std::span<Cell *const> cells) {
      for (auto *cell : cells) {
        auto const &[offset, count] = cell_ranges[range_index++];
        cell->set_range(offset, count);
        cell->staged().clear();
      }
    };
    writeback(decomposition().local_cells());
    writeback(decomposition().ghost_cells());
    assert(range_index == cell_ranges.size());
  }

  // Every staged row reference has now been copied into a committed row
  // (copy_row above), so the staging store's rows are consumed. Reset its row
  // counter here -- NOT in the decompositions' resort -- because a staged
  // {staging_store, row} reference must stay valid from the moment it is staged
  // until THIS commit reads it. Clearing at resort end would recycle staging
  // rows that the two hybrid children, or the deferred commit=false hot path,
  // still reference. Clearing at commit makes the staging counter monotonic
  // across a resort cycle and reset exactly when the referenced data is
  // consumed. (Fresh-default staged entries carry a null source store and
  // reference no staging row, so this is harmless for ghosts.)
  clear_staging_store();

  // Refresh the id -> store-row index from the freshly assigned rows (locals
  // win over ghost copies of the same id).
  rebuild_particle_index();

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Verify each cell's committed range matches the store: the store id at each
  // range row equals the id of the view at the same position, and the live
  // range has no pending-removed rows (the rebuild cleared the mask). Covers
  // local + ghost.
  auto check_cell_rows = [this](std::span<Cell *const> cells) {
    for (auto *cell : cells) {
      auto const raw_offset = cell->offset();
      auto const raw_count = cell->count();
      std::size_t k = 0u;
      for (auto const &p : cell->particles()) {
        auto const row = static_cast<int>(raw_offset + k);
        auto const stored_id = m_particle_store.id(row);
        if (stored_id != p.id()) {
          throw std::runtime_error(
              "cell range id mismatch at position " + std::to_string(k) +
              ": store.id(row) = " + std::to_string(stored_id) +
              " but view id = " + std::to_string(p.id()));
        }
        ++k;
      }
      if (k != raw_count) {
        throw std::runtime_error("cell range size mismatch (a fresh rebuild "
                                 "must leave no pending-removed rows)");
      }
    }
  };
  check_cell_rows(decomposition().local_cells());
  check_cell_rows(decomposition().ghost_cells());

  // Self-check of the permuted result: the builder produced one permutation
  // entry per new row in the SAME order the permute rebuild filled them, so the
  // store id at each new row must equal the id captured from that entry's
  // source before the swap. A mismatch means the permutation builder and the
  // permute rebuild disagree on the row order.
  if (permutation.size() != m_particle_store.number_of_particles()) {
    throw std::runtime_error(
        "resort permutation size " + std::to_string(permutation.size()) +
        " != store row count " +
        std::to_string(m_particle_store.number_of_particles()));
  }
  for (std::size_t new_row = 0u; new_row < permutation.size(); ++new_row) {
    auto const stored_id = m_particle_store.id(static_cast<int>(new_row));
    if (stored_id != debug_expected_ids[new_row]) {
      throw std::runtime_error("resort permutation id mismatch at new row " +
                               std::to_string(new_row) +
                               ": builder expected id " +
                               std::to_string(debug_expected_ids[new_row]) +
                               " but store holds " + std::to_string(stored_id));
    }
  }
  // The (offset, count) ranges must tile [0, n_total) contiguously and match
  // each cell's committed row span (the collapse written back above).
  {
    std::size_t expected_offset = 0u;
    auto check_ranges = [&expected_offset](std::span<Cell *const> cells) {
      for (auto const *cell : cells) {
        if (cell->offset() != expected_offset) {
          throw std::runtime_error("resort cell range not contiguous");
        }
        expected_offset += cell->count();
      }
    };
    check_ranges(decomposition().local_cells());
    check_ranges(decomposition().ghost_cells());
    if (expected_offset != m_particle_store.number_of_particles()) {
      throw std::runtime_error("resort cell ranges do not tile the store");
    }
  }
#endif
}

void CellStructure::ensure_staging_capacity(std::size_t const needed) {
  // Grow the staging store when it cannot hold @p needed rows. The staging
  // store's own rows must survive the growth, so a fresh larger store is built
  // and every already-staged row is copied into it via copy_row; the small
  // store is then swapped in. Capacity doubles (min 8) to amortize the copies
  // over a batch of stages.
  if (needed <= m_staging_store_capacity) {
    return;
  }
  auto new_capacity = std::max<std::size_t>(m_staging_store_capacity * 2u, 8u);
  while (new_capacity < needed) {
    new_capacity *= 2u;
  }
  ParticleStore grown{};
  grown.begin_rebuild(new_capacity, 0u);
  grown.finish_rebuild();
  for (int r = 0; r < m_staging_store_next_row; ++r) {
    grown.copy_row(m_staging_store, r, r);
  }
  using std::swap;
  swap(m_staging_store, grown);
  m_staging_store_capacity = new_capacity;
}

int CellStructure::stage_row(int const live_row) {
  ensure_staging_capacity(static_cast<std::size_t>(m_staging_store_next_row) +
                          1u);
  auto const staging_row = m_staging_store_next_row++;
  m_staging_store.copy_row(m_particle_store, live_row, staging_row);
  return staging_row;
}

int CellStructure::reserve_staging_rows(int const count) {
  assert(count >= 0);
  ensure_staging_capacity(static_cast<std::size_t>(m_staging_store_next_row) +
                          static_cast<std::size_t>(count));
  auto const first_row = m_staging_store_next_row;
  m_staging_store_next_row += count;
  return first_row;
}

Particle CellStructure::make_new_particle_view() {
  auto const staging_row = reserve_staging_rows(1);
  m_staging_store.seed_default_row(staging_row);
  return m_staging_store.make_view(staging_row);
}

void CellStructure::rebuild_particle_index() {
  assert(not m_particle_store.is_dirty());
  // Rebuild the id -> store-row map wholesale from the (synchronized) store.
  // This is the single write site of m_id_to_store_row.
  m_id_to_store_row.clear();
  auto const n_local = m_particle_store.number_of_local_particles();
  // Index locals (rows [0, n_local)); their id columns are valid.
  for (std::size_t r = 0u; r < n_local; ++r) {
    update_particle_index(m_particle_store.id(static_cast<int>(r)),
                          static_cast<int>(r));
  }
  // Also index the ALREADY-VALID ghost rows: a bare rebuild (e.g. an
  // add_particle commit, the clean-store forces.cpp sync, or resort_particles'
  // end-of-function sync) must not drop the ghost index entries a prior
  // ghosts_update established, or a ghost-id lookup (LB coupling's
  // is_ghost_for_local_particle, collision detection) would see no row. FRESH
  // ghost rows just created by resize_ghost_storage carry a default id (-1) and
  // are skipped by index_ghost_particles; they are picked up by a later
  // index_ghost_particles once ghosts_update fills their id columns.
  index_ghost_particles();
}

void CellStructure::index_ghost_particles() {
  assert(not m_particle_store.is_dirty());
  auto const n_total = m_particle_store.number_of_particles();
  auto const n_local = m_particle_store.number_of_local_particles();
  // Ghost rows [n_local, n_total): record the store row for each ghost whose id
  // is VALID (>= 0; a freshly resized ghost carries -1 until ghosts_update
  // fills it) and not already mapped by a local (locals win) or by an earlier
  // ghost of the same id (first valid row wins). The local prefix of the map is
  // left intact; this appends the ghost tail. Idempotent: a ghost id already
  // present is skipped, so re-running after ghosts_update only adds the
  // newly-valid ghosts.
  for (std::size_t r = n_local; r < n_total; ++r) {
    auto const id = m_particle_store.id(static_cast<int>(r));
    if (id < 0 or get_local_particle(id).has_value()) {
      continue;
    }
    update_particle_index(id, static_cast<int>(r));
  }
}

void CellStructure::check_particle_index() const {
  auto const max_id = get_max_local_particle_id();

  for (auto const &p : local_particles()) {
    auto const id = p.id();

    if (id < 0 or id > max_id) {
      throw std::runtime_error("Particle id out of bounds.");
    }

    // Particles are by-value views over store rows, so the resolved view and
    // the iterated view have DIFFERENT addresses even for the same particle.
    // Check identity by store row instead of by address: both must resolve to
    // the same row of the same store.
    auto const indexed = get_local_particle(id);
    if (not indexed or indexed->store() != p.store() or
        indexed->store_row() != p.store_row()) {
      throw std::runtime_error("Invalid local particle index entry.");
    }
  }

  /* checks: local particle id. The map also holds ghost entries (for ghost-id
   * lookups), so count only the LOCAL (non-ghost) entries when comparing
   * against local_particles(). */
  std::size_t local_part_cnt = 0u;
  for (int n = 0; n < get_max_local_particle_id() + 1; n++) {
    auto const indexed = get_local_particle(n);
    if (indexed) {
      if (indexed->id() != n) {
        throw std::runtime_error("local_particles part has corrupted id.");
      }
      if (not indexed->is_ghost()) {
        local_part_cnt++;
      }
    }
  }

  if (local_part_cnt != local_particles().size()) {
    throw std::runtime_error(
        std::to_string(local_particles().size()) + " parts in cells but " +
        std::to_string(local_part_cnt) + " parts in local_particles");
  }
}

void CellStructure::check_particle_sorting() const {
  for (auto cell : decomposition().local_cells()) {
    for (auto const &p : cell->particles()) {
      if (particle_to_cell(p) != cell) {
        throw std::runtime_error("misplaced particle with id " +
                                 std::to_string(p.id()));
      }
    }
  }
}

void CellStructure::remove_particle(int id) {
  auto remove_all_bonds_to = [id](BondList &bl) {
    for (auto it = bl.begin(); it != bl.end();) {
      if (Utils::contains(it->partner_ids(), id)) {
        it = bl.erase(it);
      } else {
        std::advance(it, 1);
      }
    }
  };

  for (auto cell : decomposition().local_cells()) {
    cell->set_store(m_particle_store);
    // Iterate the committed rows by RAW position. drop_row marks the matching
    // row pending-removed (the live iteration skips it immediately, so
    // the removed particle is invisible before the next rebuild resolves it); a
    // kept particle has any bond to the removed id stripped (written through
    // the view into the store's bond sidecar, preserved by the next rebuild). A
    // row already marked pending (a prior drop) is skipped.
    auto const offset = cell->offset();
    auto const count = cell->count();
    for (std::size_t index = 0u; index < count; ++index) {
      auto const row = static_cast<int>(offset + index);
      if (m_particle_store.is_pending_removal(row)) {
        continue;
      }
      auto view = m_particle_store.make_view(row);
      if (view.id() == id) {
        CellParticleStorage::drop_row(*cell, index);
      } else {
        remove_all_bonds_to(view.bonds());
      }
    }
  }
  // Eagerly invalidate the id -> store-row entry so no stale row is readable
  // through get_local_particle before the next store rebuild. The full map is
  // still rebuilt wholesale at the next ensure_particle_store_synchronized
  // (which renumbers every row); this just closes the intervening window.
  if (id >= 0 and static_cast<unsigned int>(id) < m_id_to_store_row.size()) {
    m_id_to_store_row[static_cast<unsigned int>(id)] = no_store_row;
  }
  mark_particle_store_dirty();
}

std::optional<Particle> CellStructure::add_local_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  if (sort_cell) {
    sort_cell->set_store(m_particle_store);
    auto const id = p.id();
    append_staged_particle(*sort_cell, std::move(p));
    mark_particle_store_dirty();
    // Commit immediately so the new particle is live (indexed, iterable,
    // readable) right after this call. Return a fresh by-value view resolved
    // via the id -> store-row map.
    ensure_particle_store_synchronized();
    return get_local_particle(id);
  }

  return std::nullopt;
}

std::optional<Particle> CellStructure::add_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  /* There is always at least one cell, so if the particle
   * does not belong to a cell on this node we can put it there. */
  auto cell = sort_cell ? sort_cell : decomposition().local_cells()[0];
  cell->set_store(m_particle_store);

  /* If the particle isn't local a global resort may be
   * needed, otherwise a local resort if sufficient. */
  set_resort_particles(sort_cell ? Cells::RESORT_LOCAL : Cells::RESORT_GLOBAL);

  auto const id = p.id();
  append_staged_particle(*cell, std::move(p));
  mark_particle_store_dirty();
  // Commit immediately so the new particle is live right after this call.
  // Return a fresh by-value view resolved via the id -> store-row map.
  ensure_particle_store_synchronized();
  return get_local_particle(id);
}

int CellStructure::get_max_local_particle_id() const {
  // The id -> store-row map is indexed by id; the highest id present on this
  // rank is the highest index carrying a valid (non-sentinel) row.
  auto const it =
      std::ranges::find_if(std::ranges::views::reverse(m_id_to_store_row),
                           [](int const row) { return row != no_store_row; });

  return (it != m_id_to_store_row.rend())
             ? static_cast<int>(m_id_to_store_row.rend() - it) - 1
             : -1;
}

int CellStructure::get_local_pair_bond_numbers() const {
  return m_bond_state->pair_count;
}
int CellStructure::get_local_angle_bond_numbers() const {
  return m_bond_state->angle_count;
}
int CellStructure::get_local_dihedral_bond_numbers() const {
  return m_bond_state->dihedral_count;
}
void CellStructure::set_local_bond_numbers(int pair_value, int angle_value,
                                           int dihedral_value) {
  m_bond_state->set_counts(pair_value, angle_value, dihedral_value);
}
#ifdef ESPRESSO_COLLISION_DETECTION
void CellStructure::clear_new_bonds() { m_bond_state->clear_new_bonds(); }
void CellStructure::add_new_bond(int bond_id,
                                 std::vector<int> const &particle_ids) {
  m_bond_state->add_new_bond(bond_id, particle_ids, get_id_to_index());
}
void CellStructure::rebuild_bond_list() { m_bond_state->rebuild(); }
#endif // ESPRESSO_COLLISION_DETECTION

void CellStructure::remove_all_particles() {
  for (auto cell : decomposition().local_cells()) {
    CellParticleStorage::clear_particles(*cell);
  }
  // Also clear the GHOST cells: they hold ghost copies of the now-deleted
  // particles. If left in place, the next store rebuild would
  // re-index those stale ghost rows (rebuild_particle_index indexes valid
  // ghosts), so get_local_particle would still report the deleted ids as
  // present -- e.g. a fresh add of the same id would wrongly see "already
  // exists". Emptying them keeps the index consistent with the emptied system.
  for (auto cell : decomposition().ghost_cells()) {
    CellParticleStorage::clear_particles(*cell);
  }

  clear_particle_index();
  clear_bond_properties();
  // clear_particles() above only resets each cell's row range and staging
  // buffer; the store still holds the rows, so it must be marked dirty for the
  // next rebuild to drop them.
  mark_particle_store_dirty();
  clear_particle_node();
  get_system().on_particle_change();
}

/* Map the data parts flags from cells to those used internally
 * by the ghost communication */
unsigned map_data_parts(unsigned data_parts) {
  using namespace Cells;

  /* clang-format off */
  return GHOSTTRANS_NONE
         | ((data_parts & DATA_PART_PROPERTIES) ? GHOSTTRANS_PROPRTS : 0u)
         | ((data_parts & DATA_PART_POSITION) ? GHOSTTRANS_POSITION : 0u)
         | ((data_parts & DATA_PART_MOMENTUM) ? GHOSTTRANS_MOMENTUM : 0u)
         | ((data_parts & DATA_PART_FORCE) ? GHOSTTRANS_FORCE : 0u)
#ifdef ESPRESSO_BOND_CONSTRAINT
         | ((data_parts & DATA_PART_RATTLE) ? GHOSTTRANS_RATTLE : 0u)
#endif
#ifdef ESPRESSO_ROTATION
         | ((data_parts & DATA_PART_QUAT) ? GHOSTTRANS_QUAT : 0u)
         | ((data_parts & DATA_PART_TORQUE) ? GHOSTTRANS_TORQUE : 0u)
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
         | ((data_parts & DATA_PART_DIPFLD) ? GHOSTTRANS_DIPFLD : 0u)
#endif
         | ((data_parts & DATA_PART_BONDS) ? GHOSTTRANS_BONDS : 0u);
  /* clang-format on */
}

void CellStructure::ghosts_count() {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  GhostComm::halo_exchange(
      *decomposition().halo_plan(), *get_system().box_geo, GHOSTTRANS_PARTNUM,
      {GhostComm::Direction::Push, GhostComm::Combine::Overwrite},
      m_ghost_buffers);
  // The PARTNUM step resized the ghost cells, i.e. it staged fresh ghost rows;
  // the store must be rebuilt before the DATA step reads committed views.
  mark_particle_store_dirty();
}

void CellStructure::ghosts_update(unsigned data_parts) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  auto const parts = map_data_parts(data_parts);
  GhostComm::halo_exchange(
      *decomposition().halo_plan(), *get_system().box_geo, parts,
      {GhostComm::Direction::Push, GhostComm::Combine::Overwrite},
      m_ghost_buffers);
}

void CellStructure::ghosts_reduce_forces() {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  GhostComm::halo_exchange(
      *decomposition().halo_plan(), *get_system().box_geo,
      get_system().get_force_reduce_ghost_flags(),
      {GhostComm::Direction::Reduce, GhostComm::Combine::Add}, m_ghost_buffers);
}

#ifdef ESPRESSO_CALIPER
// caliper annotation for split phase ghost forces reduction
static cali_id_t ghost_reduce_async_attr() {
  static const cali_id_t id =
      espresso_cali_active()
          ? cali_create_attribute("ghosts_reduce_forces_async",
                                  CALI_TYPE_STRING,
                                  CALI_ATTR_ASVALUE | CALI_ATTR_SCOPE_THREAD)
          : CALI_INV_ID;
  return id;
}
#endif // ESPRESSO_CALIPER

void CellStructure::ghosts_reduce_forces_start() {
  assert(not m_pending_ghost_reduce.has_value() &&
         "ghosts_reduce_forces_start: a reduction is already in flight");
  // BEGIN fires only after emplace succeeds so a throwing start cannot
  // leave the Caliper region open without a matching END.
  m_pending_ghost_reduce.emplace(GhostComm::halo_exchange_start(
      *decomposition().halo_plan(), *get_system().box_geo,
      get_system().get_force_reduce_ghost_flags(),
      {GhostComm::Direction::Reduce, GhostComm::Combine::Add},
      m_ghost_buffers));
#ifdef ESPRESSO_CALIPER
  if (auto id = ghost_reduce_async_attr(); id != CALI_INV_ID)
    cali_begin_string(id, "in_flight");
#endif
}

void CellStructure::ghosts_reduce_forces_finish() {
  assert(m_pending_ghost_reduce.has_value() &&
         "ghosts_reduce_forces_finish: no reduction is in flight");
  try {
    GhostComm::halo_exchange_finish(*m_pending_ghost_reduce);
  } catch (...) {
    // A failed finish cannot be retried: the exchange state is half-consumed
    // (some requests waited, some buffers unpacked). Drop the pending state so
    // a later finish attempt (e.g. the ReduceGuard in integrate.cpp) does not
    // re-run MPI waits on completed requests.
    m_pending_ghost_reduce.reset();
#ifdef ESPRESSO_CALIPER
    if (auto id = ghost_reduce_async_attr(); id != CALI_INV_ID)
      cali_end(id);
#endif
    throw;
  }
#ifdef ESPRESSO_CALIPER
  // END before reset so the region is closed before the optional is cleared.
  if (auto id = ghost_reduce_async_attr(); id != CALI_INV_ID)
    cali_end(id);
#endif
  m_pending_ghost_reduce.reset();
}
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void CellStructure::ghosts_reduce_dipole_field() {
  GhostComm::halo_exchange(
      *decomposition().halo_plan(), *get_system().box_geo, GHOSTTRANS_DIPFLD,
      {GhostComm::Direction::Reduce, GhostComm::Combine::Add}, m_ghost_buffers);
}
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
void CellStructure::ghosts_reduce_rattle_correction() {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  GhostComm::halo_exchange(
      *decomposition().halo_plan(), *get_system().box_geo, GHOSTTRANS_RATTLE,
      {GhostComm::Direction::Reduce, GhostComm::Combine::Add}, m_ghost_buffers);
}
#endif

void CellStructure::resort_particles(bool global_flag, bool commit) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  assert(not m_pending_ghost_reduce.has_value() &&
         "resort_particles: ghost force reduction is still in flight — "
         "call ghosts_reduce_forces_finish() first");
  // The id -> store-row index is rebuilt wholesale from the store by
  // ensure_particle_store_synchronized after the resort (see
  // rebuild_particle_index), so no incremental index bookkeeping is needed
  // during the resort. The `diff` is still collected (the decomposition
  // contract) but only for its "cells touched" side effect; index correctness
  // comes from the rebuild. Clear the stale index now so nothing reads a
  // dangling entry in the resort window (get_local_particle is not called
  // there).
  clear_particle_index();

  // Drop every ghost row. A resort invalidates the ghost layer wholesale, and
  // the id -> store-row map is rebuilt below from ALL live rows, ghosts
  // included (locals win, but an id owned by nobody local still resolves to its
  // deduped ghost row -- that is the map's documented contract). Carrying stale
  // ghosts into that rebuild would resurrect ids the resort just removed: a
  // caller that treats "resolves" as "exists" (e.g. the ParticleHandle
  // creation guard) would then refuse to re-create a freed id on the rank that
  // still held its ghost. This mirrors invalidate_ghosts() in the AoS layout,
  // which drops ghost entries from the index at the same point. The ghosts are
  // re-established by the very next ghosts_count + ghost update.
  for (auto *cell : m_decomposition->ghost_cells()) {
    cell->set_store(m_particle_store);
    CellParticleStorage::resize_ghost_storage(*cell, 0ul);
  }

  std::vector<ParticleChange> diff;

  // Ensure every cell knows the store before the decomposition mutates cell
  // storage (extract_row / insert / ghost resize all need Cell::store()).
  for (auto *cell : m_decomposition->local_cells()) {
    cell->set_store(m_particle_store);
  }
  for (auto *cell : m_decomposition->ghost_cells()) {
    cell->set_store(m_particle_store);
  }

  // Mark the store dirty BEFORE the decomposition's resort: the dirty flag must
  // be truthful DURING the resort window. HybridDecomposition runs an internal
  // PARTNUM+ghost update inside that window, and the columnar ghost paths must
  // see dirty and fall back to the per-particle path there.
  mark_particle_store_dirty();
  m_decomposition->resort(global_flag, diff);
  static_cast<void>(diff);

  // Commit the staged (migrated/new) particles into store rows and rebuild the
  // id -> store-row index: resort cleared the index and staged the moved
  // particles, so a caller reading the index right after resort (e.g. the
  // ADDITIONAL_CHECKS below, or a direct/unit-test resort_particles call) must
  // see a consistent, populated index. Rank-local (no MPI), so it does not
  // affect the collective ordering the caller relies on.
  //
  // The hot-path caller (update_ghosts_and_resort_particle) passes commit=false
  // to DEFER this to a single rebuild AFTER ghosts_count, so locals are not
  // committed here only to be re-copied by a second whole-store rebuild moments
  // later. ghosts_count sizes ghost cells from the SOURCE cell's Cell::size
  // (committed rows + staged), so it does not need the locals committed first;
  // and nothing reads the index in the deferred window.
  if (commit) {
    ensure_particle_store_synchronized();
  }

  auto const &lebc = get_system().box_geo->lees_edwards_bc();
  m_rebuild_verlet_list = true;
  m_rebuild_verlet_list_cabana = true;
  m_le_pos_offset_at_last_resort = lebc.pos_offset;

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // These checks read the id -> store-row index / sorted store, which only
  // exist after the commit above; when the commit is deferred, the caller runs
  // its own sync and the equivalent checks
  // (ensure_particle_store_synchronized's cell-row validation) after
  // ghosts_count.
  if (commit) {
    check_particle_index();
    check_particle_sorting();
  }
#endif
}

void CellStructure::set_atom_decomposition() {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  auto atom = std::make_unique<AtomDecomposition>(::comm_cart, box_geo);
  // Install the per-field migration staging handle: the decomposition's resort
  // routes migrating particles through this store's staging store.
  atom->set_migration_staging(make_migration_staging());
  set_particle_decomposition(std::move(atom));
  m_type = CellStructureType::NSQUARE;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_regular_decomposition(
    double range, std::optional<std::pair<int, int>> fully_connected_boundary) {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  auto regular = std::make_unique<RegularDecomposition>(
      ::comm_cart, range, box_geo, local_geo, fully_connected_boundary);
  regular->set_migration_staging(make_migration_staging());
  set_particle_decomposition(std::move(regular));
  m_type = CellStructureType::REGULAR;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_hybrid_decomposition(double cutoff_regular,
                                             std::set<int> n_square_types) {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  auto hybrid = std::make_unique<HybridDecomposition>(
      ::comm_cart, cutoff_regular, m_verlet_skin,
      [&system]() { return system.get_global_ghost_flags(); }, box_geo,
      local_geo, n_square_types);
  // Let the hybrid decomposition commit staged particles to store rows between
  // its internal resort and ghost communications (see
  // HybridDecomposition::resort). Mark the store dirty first: the ghost resize
  // (resize_ghost_storage) stages ghosts WITHOUT marking dirty, so a bare
  // ensure_particle_store_synchronized would early-return (clean store) and
  // never commit them; forcing the rebuild guarantees the staged particles
  // (migrated locals AND freshly resized ghosts) become committed rows.
  hybrid->set_commit_store([this]() {
    mark_particle_store_dirty();
    ensure_particle_store_synchronized();
  });
  // Install the per-field migration staging handle; the hybrid propagates it to
  // its two child decompositions, which run the wire exchange.
  hybrid->set_migration_staging(make_migration_staging());
  set_particle_decomposition(std::move(hybrid));
  m_type = CellStructureType::HYBRID;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_verlet_skin(double value) {
  assert(value >= 0.);
  m_verlet_skin = value;
  m_verlet_skin_set = true;
  m_rebuild_verlet_list_cabana = true;
  get_system().on_verlet_skin_change();
}

void CellStructure::set_verlet_skin_heuristic() {
  assert(not is_verlet_skin_set());
  auto const max_cut = get_system().maximal_cutoff();
  if (max_cut <= 0.) {
    throw std::runtime_error(
        "cannot automatically determine skin, please set it manually");
  }
  /* maximal skin that can be used without resorting is the maximal
   * range of the cell system minus what is needed for interactions. */
  auto const max_range = std::ranges::min(max_cutoff());
  auto const new_skin = std::min(0.4 * max_cut, max_range - max_cut);
  set_verlet_skin(new_skin);
}

void CellStructure::update_ghosts_and_resort_particle(unsigned data_parts) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  /* data parts that are only updated on resort */
  auto constexpr resort_only_parts =
      Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS;

  auto const global_resort = boost::mpi::all_reduce(
      ::comm_cart, m_resort_particles, std::bit_or<unsigned>());

  if (global_resort != Cells::RESORT_NONE) {
    auto const do_global_resort = (global_resort & Cells::RESORT_GLOBAL) != 0;

    /* Resort cell system. Defer the local commit (commit=false) so it is folded
     * into the SINGLE post-ghosts_count store rebuild below, instead of
     * committing locals here and re-copying every local column a second time in
     * that rebuild. */
    resort_particles(do_global_resort, /*commit=*/false);
    /* After resort the migrated/new locals are STAGED (not yet committed to
     * store rows). ghosts_count sizes every ghost cell from its source (local)
     * cell's Cell::size, which counts committed rows + staged, so staged locals
     * are counted correctly and do not need committing first. Count ghosts
     * (stages ghost slots + marks the store dirty), then rebuild ONCE to commit
     * locals AND the freshly staged ghosts together. */
    ghosts_count();
    /* Rebuild the store rows NOW: resort staged locals and ghosts_count marked
     * the store dirty (staged ghost slots); a dirty store forces the ghost
     * update below onto the slow per-field accessor path (measured +49% in this
     * slot on lj-4rank). After this single rebuild, locals and freshly created
     * ghosts are all attached and the update writes their state straight into
     * the columns by row. */
    ensure_particle_store_synchronized();
    ghosts_update(data_parts);

    /* Index the ghosts now that ghosts_update has filled their id columns (the
     * rebuild inside ensure_particle_store_synchronized indexed only locals,
     * whose ids were valid; freshly created ghost rows carried a default id
     * until this update). Records the store row for each ghost whose id is not
     * already owned by a local. */
    index_ghost_particles();

    /* Particles are now sorted */
    clear_resort_particles();
  } else {
    /* Communication step: ghost information */
    ghosts_update(data_parts & ~resort_only_parts);
  }
}

bool CellStructure::check_resort_required(
    Utils::Vector3d const &additional_offset) const {
  auto const lim = Utils::sqr(m_verlet_skin / 2.) - additional_offset.norm2();

  auto add_partial = [lim](bool &result, Particle const &p) {
    if ((Utils::Vector3d(p.pos()) -
         Utils::Vector3d(p.pos_at_last_verlet_update()))
            .norm2() > lim) {
      result = true;
    }
  };

  auto reduce_op = [](bool &acc, bool const &val) { acc |= val; };

  return reduce_over_local_particles<bool>(*this, add_partial, reduce_op);
}

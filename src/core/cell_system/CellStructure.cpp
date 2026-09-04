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
#include "cell_system/RegularDecomposition.hpp"

#include "BoxGeometry.hpp"
#include "LocalBondState.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "aosoa_pack.hpp"
#include "bonds.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloExchange.hpp"
#include "integrators/Propagation.hpp"
#include "kokkos_helpers.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "particle_enumeration.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

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
  kokkos_deep_copy(execution_space{}, get_aosoa().flags, uint8_t{0});
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
    // Only primary entries drive force/energy calculation; mirror entries
    // (held by non-owning participants) are a query/removal-only detail
    // and would otherwise double-count every bond.
    if (not bond.is_primary()) {
      continue;
    }
    auto const partner_ids = bond.partner_ids();
    try {
      auto const partners = resolve_bond_partners(partner_ids);
      if (partners.size() == 1u) { // pair bonds
        auto p_index = Kokkos::atomic_fetch_add(&pair_count, 1);
        pair_list(p_index, 0) = p.id();
        pair_list(p_index, 1) = partners[0]->id();
        pair_ids(p_index) = bond.bond_id();
      } else if (partners.size() == 2u) { // angle bond
        auto a_index = Kokkos::atomic_fetch_add(&angle_count, 1);
        angle_list(a_index, 0) = p.id();
        angle_list(a_index, 1) = partners[0]->id();
        angle_list(a_index, 2) = partners[1]->id();
        angle_ids(a_index) = bond.bond_id();
      } else if (partners.size() == 3u) { // dihedral bond
        auto d_index = Kokkos::atomic_fetch_add(&dihedral_count, 1);
        dihedral_list(d_index, 0) = p.id();
        dihedral_list(d_index, 1) = partners[0]->id();
        dihedral_list(d_index, 2) = partners[1]->id();
        dihedral_list(d_index, 3) = partners[2]->id();
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
  auto &unique_particles = m_unique_particles;
  unique_particles.clear();
  unique_particles.resize(count_local_particles());
  std::unordered_set<int> registered_index{};
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  int n_threads = execution_space().concurrency();

  m_bond_state->reset_counts();
  // one cache line per thread: these counters are written on every particle,
  // so packing them into shared cache lines makes the sweep bounce lines
  // between L3 domains
  struct alignas(64) PerThreadCounts {
    int max_id = 0;
    int pair = 0;
    int angle = 0;
    int dihedral = 0;
  };
  std::vector<PerThreadCounts> thread_counts(n_threads);

  enumerate_local_particles(*this, [&unique_particles, &thread_counts](
                                       std::size_t index, Particle &p) {
    unique_particles[index] = &p;
    auto &counts = thread_counts[omp_get_thread_num()];
    counts.max_id = std::max(p.id(), counts.max_id);
    // Read the bond list's incrementally-maintained primary-entry counts
    // instead of walking and decoding every entry (mirror entries
    // included) to recount them on every rebuild; update_bond_storage()
    // below still only fills the flat lists from primary entries, and
    // the two must agree on the count or the Kokkos views below overflow.
    auto const &bond_counts = p.bonds().primary_counts();
    counts.pair += bond_counts.pair;
    counts.angle += bond_counts.angle;
    counts.dihedral += bond_counts.dihedral;
  });
  Kokkos::fence();
  int pair_count = 0;
  int angle_count = 0;
  int dihedral_count = 0;
  int max_id = 0;
  for (auto const &counts : thread_counts) {
    pair_count += counts.pair;
    angle_count += counts.angle;
    dihedral_count += counts.dihedral;
    max_id = std::max(counts.max_id, max_id);
  }
  set_local_bond_numbers(pair_count, angle_count, dihedral_count);
  m_bond_state->allocate();
  for (auto &p : ghost_particles()) {
    auto const *local_particle = get_local_particle(p.id());
    if (not local_particle or not local_particle->is_ghost()) {
      continue;
    }
    if (registered_index.contains(p.id())) {
      continue;
    }
    registered_index.insert(p.id());
    unique_particles.emplace_back(&p);
    max_id = std::max(p.id(), max_id);
  }
  registered_index.clear();
  m_cached_max_local_particle_id = max_id;
  m_num_local_particles_cached = unique_particles.size();
}

CellStructure::CellStructure(BoxGeometry const &box)
    : m_decomposition{std::make_unique<AtomDecomposition>(box)} {}

void CellStructure::check_particle_index() const {
  auto const max_id = get_max_local_particle_id();

  for (auto const &p : local_particles()) {
    auto const id = p.id();

    if (id < 0 or id > max_id) {
      throw std::runtime_error("Particle id out of bounds.");
    }

    if (get_local_particle(id) != &p) {
      throw std::runtime_error("Invalid local particle index entry.");
    }
  }

  /* checks: local particle id */
  std::size_t local_part_cnt = 0u;
  for (int n = 0; n < get_max_local_particle_id() + 1; n++) {
    if (get_local_particle(n) != nullptr) {
      local_part_cnt++;
      if (get_local_particle(n)->id() != n) {
        throw std::runtime_error("local_particles part has corrupted id.");
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
  // The particle's own bond list already names every bond it is involved
  // in, as owner or as a mirror-holding participant (see BondList.hpp), so
  // the other participants needing cleanup can be found directly instead
  // of sweeping every local particle. ::remove_bond() erases the matching
  // entry from each of them; the entry on this particle itself is skipped
  // (via `id`), since its whole bond list is discarded below regardless.
  if (auto const *p = get_local_particle(id)) {
    std::vector<std::pair<int, std::vector<int>>> bonds_to_remove;
    for (auto const bond : p->bonds()) {
      std::vector<int> ids = {id};
      std::ranges::copy(bond.partner_ids(), std::back_inserter(ids));
      bonds_to_remove.emplace_back(bond.bond_id(), std::move(ids));
    }
    auto &system = get_system();
    for (auto const &[bond_id, ids] : bonds_to_remove) {
      ::remove_bond(system, bond_id, ids, id);
    }
  }

  for (auto cell : decomposition().local_cells()) {
    auto &parts = cell->particles();
    for (auto it = parts.begin(); it != parts.end();) {
      if (it->id() == id) {
        it = parts.erase(it);
        update_particle_index(id, nullptr);
        update_particle_index(parts);
      } else {
        it++;
      }
    }
  }
}

Particle *CellStructure::add_local_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  if (sort_cell) {
    return std::addressof(
        append_indexed_particle(sort_cell->particles(), std::move(p)));
  }

  return {};
}

Particle *CellStructure::add_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  /* There is always at least one cell, so if the particle
   * does not belong to a cell on this node we can put it there. */
  auto cell = sort_cell ? sort_cell : decomposition().local_cells()[0];

  /* If the particle isn't local a global resort may be
   * needed, otherwise a local resort if sufficient. */
  set_resort_particles(sort_cell ? Cells::RESORT_LOCAL : Cells::RESORT_GLOBAL);

  return std::addressof(
      append_indexed_particle(cell->particles(), std::move(p)));
}

int CellStructure::get_max_local_particle_id() const {
  auto it = std::ranges::find_if(std::ranges::views::reverse(m_particle_index),
                                 [](auto const *p) { return p != nullptr; });

  return (it != m_particle_index.rend()) ? (*it)->id() : -1;
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
    cell->particles().clear();
  }

  m_particle_index.clear();
  clear_bond_properties();
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

namespace {
/**
 * @brief Apply a @ref ParticleChange to a particle index.
 */
struct UpdateParticleIndexVisitor {
  CellStructure *cs;

  void operator()(RemovedParticle rp) const {
    cs->update_particle_index(rp.id, nullptr);
  }
  void operator()(ModifiedList mp) const { cs->update_particle_index(mp.pl); }
};
} // namespace

void CellStructure::resort_particles(bool global_flag) {
#ifdef ESPRESSO_CALIPER
  ESPRESSO_CALI_MARK_FUNCTION;
#endif
  assert(not m_pending_ghost_reduce.has_value() &&
         "resort_particles: ghost force reduction is still in flight — "
         "call ghosts_reduce_forces_finish() first");
  invalidate_ghosts();

  std::vector<ParticleChange> diff;

  m_decomposition->resort(global_flag, diff);

  for (auto d : diff) {
    std::visit(UpdateParticleIndexVisitor{this}, d);
  }

  auto const &lebc = get_system().box_geo->lees_edwards_bc();
  m_rebuild_verlet_list = true;
  m_rebuild_verlet_list_cabana = true;
  m_le_pos_offset_at_last_resort = lebc.pos_offset;

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  check_particle_index();
  check_particle_sorting();
#endif
}

void CellStructure::set_atom_decomposition() {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  set_particle_decomposition(
      std::make_unique<AtomDecomposition>(::comm_cart, box_geo));
  m_type = CellStructureType::NSQUARE;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_regular_decomposition(
    double range, std::optional<std::pair<int, int>> fully_connected_boundary) {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  set_particle_decomposition(std::make_unique<RegularDecomposition>(
      ::comm_cart, range, box_geo, local_geo, fully_connected_boundary));
  m_type = CellStructureType::REGULAR;
  local_geo.set_cell_structure_type(m_type);
  system.on_cell_structure_change();
}

void CellStructure::set_hybrid_decomposition(double cutoff_regular,
                                             std::set<int> n_square_types) {
  auto &system = get_system();
  auto &local_geo = *system.local_geo;
  auto const &box_geo = *system.box_geo;
  set_particle_decomposition(std::make_unique<HybridDecomposition>(
      ::comm_cart, cutoff_regular, m_verlet_skin,
      [&system]() { return system.get_global_ghost_flags(); }, box_geo,
      local_geo, n_square_types));
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

    /* Resort cell system */
    resort_particles(do_global_resort);
    ghosts_count();
    ghosts_update(data_parts);

    /* Add the ghost particles to the index if we don't already
     * have them. */
    for (auto &p : ghost_particles()) {
      if (get_local_particle(p.id()) == nullptr) {
        update_particle_index(p.id(), &p);
      }
    }

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
    if ((p.pos() - p.pos_at_last_verlet_update()).norm2() > lim) {
      result = true;
    }
  };

  auto reduce_op = [](bool &acc, bool const &val) { acc |= val; };

  return reduce_over_local_particles<bool>(*this, add_partial, reduce_op);
}

/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "BondList.hpp"
#include "particle_store/ParticleParameters.hpp"
#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/quaternion.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

struct Particle; // attach_to_store is defined in Particle.hpp

/**
 * @brief Array-based particle storage (structure of arrays).
 *
 * Owns per-particle quantities in a single index space: local particles
 * first (cell-traversal order), ghosts appended. Vector/quaternion fields are
 * component-major (@ref StateVectorLayout) Kokkos columns by default; the
 * remaining scalar and cold parameters live in scalar columns and host
 * sidecars. The columns are the source of truth for attached particles.
 *
 * Rebuild protocol: mark_dirty() on any topology change; the owner
 * (CellStructure) later renumbers the rows via @ref permute_rebuild (or the
 * standalone @ref begin_rebuild / @ref assign_row / @ref finish_rebuild fills).
 * A rebuild preserves surviving rows by their old row index and seeds
 * genuinely-new rows to the new-particle defaults (@ref seed_default_row).
 * Rebuilds are purely rank-local (no MPI).
 */
class ParticleStore {
public:
  /**
   * @brief Memory layout of the vector/quaternion state columns.
   *
   * COMPONENT-MAJOR (LayoutLeft): each component's values are contiguous
   * across particles. This is the default; it measures slightly faster
   * (order-balanced A/B: p3m 4-rank -1.7%, p3m 1-rank -1.3%, lj 1-rank and
   * 4-rank -0.5%, lj-omp neutral) because per-field ghost exchange and column
   * kernels leave no gather-dominated paths that would favour particle-major,
   * and per-component ghost legs benefit from contiguity. Component-major is
   * also the coalescing-friendly layout for device-resident GPU work.
   *
   * Particle-major (LayoutRight) remains selectable with
   * -DESPRESSO_PARTICLE_STORE_PARTICLE_MAJOR for A/B on identical sources; all
   * column access is layout-agnostic (typed views, (row, component) indexing,
   * stride-carrying contexts/proxies) and trajectories are bitwise-identical
   * under either layout.
   */
#ifdef ESPRESSO_PARTICLE_STORE_PARTICLE_MAJOR
  using StateVectorLayout = Kokkos::LayoutRight;
#else
  using StateVectorLayout = Kokkos::LayoutLeft;
#endif
  using Column = Kokkos::DualView<double *[3], StateVectorLayout>;
  using IntegerColumn = Kokkos::DualView<int *[3], StateVectorLayout>;
  using QuaternionColumn = Kokkos::DualView<double *[4], StateVectorLayout>;
  using ScalarColumn = Kokkos::DualView<double *>;
  using ShortColumn = Kokkos::DualView<short *>;
  // Parameter scalar columns: int-typed (id/mol_id/type/propagation) and
  // uint8-typed bitfields (rotation/ext_flag).
  using IntScalarColumn = Kokkos::DualView<int *>;
  using Uint8Column = Kokkos::DualView<std::uint8_t *>;

  // Per-particle friction (gamma/gamma_rot): scalar when isotropic, a 3-vector
  // column when ESPRESSO_PARTICLE_ANISOTROPY selects per-axis friction. The
  // value type and column follow the same switch as the
  // ParticleProperties::gamma member (see Particle.hpp).
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  using GammaValue = Utils::Vector3d;
  using GammaColumn = Column;
#else
  using GammaValue = double;
  using GammaColumn = ScalarColumn;
#endif

  std::size_t number_of_local_particles() const {
    return m_number_of_local_particles;
  }
  std::size_t number_of_ghost_particles() const {
    return m_number_of_ghost_particles;
  }
  std::size_t number_of_particles() const {
    return m_number_of_local_particles + m_number_of_ghost_particles;
  }

  void mark_dirty() { m_dirty = true; }
  /** @brief Monotonic row-assignment generation; bumped by finish_rebuild.
   *  Consumers may cache row-derived data (e.g. ghost-communication row
   *  ranges) tagged with this value and reuse it while it is unchanged. */
  std::uint64_t generation() const { return m_generation; }
  bool is_dirty() const { return m_dirty; }

  // -- mid-epoch pending removal --------------------------------------------
  // Cells are contiguous (offset, count) ranges, so a row cannot be
  // swap-removed from a cell between rebuilds (a range has no shuffling op). A
  // removal instead marks the store row PENDING-REMOVED here: iteration
  // (@ref RowParticleRange / @ref CellRowSpan) skips such rows immediately, and
  // the next permutation rebuild resolves it by NOT including the row in the
  // permutation (the row is dropped). The mask is indexed by store row and is
  // cleared on every rebuild (begin_rebuild / permute_rebuild size it fresh and
  // finish_rebuild leaves it all-clear for the new generation). Invariant:
  // iteration must not see a removed particle before the rebuild renumbers
  // rows.
  void mark_pending_removal(int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    if (m_pending_removal[static_cast<std::size_t>(row)] == char{0}) {
      m_pending_removal[static_cast<std::size_t>(row)] = char{1};
      ++m_pending_removal_count;
    }
  }
  bool is_pending_removal(int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    return m_pending_removal[static_cast<std::size_t>(row)] != char{0};
  }
  /** @brief Whether ANY row is marked pending-removal. O(1). Removal marks
   *  only exist between a particle removal and the next rebuild; on every
   *  steady-state iteration path this is false, letting iterators skip the
   *  per-element mask load (measured at ~40% of the Verlet-build slot on the
   *  lj benchmark when checked per candidate). */
  bool has_pending_removals() const { return m_pending_removal_count != 0u; }

  void begin_rebuild(std::size_t number_of_local_particles,
                     std::size_t number_of_ghost_particles);
  void assign_row(Particle &particle, int row);
  void finish_rebuild();

  /**
   * @brief Copy every per-particle field of one store row into another.
   *
   * Copies EVERY column, POD sidecar and ragged sidecar (bonds / exclusions,
   * under matching ifdefs) of @p source_row in @p source into @p
   * destination_row of THIS store, row-to-row, without going through a @ref
   * Particle. @p source and @c *this may be the same store (a row can be copied
   * within one store) or two different stores. Both rows must be currently
   * valid indices in their respective stores.
   *
   * This is the machinery shared by the migration `extract` (live store ->
   * staging store) and `receive` (staging store -> live store) paths: a full,
   * value-preserving row transfer.
   *
   * The copied field set is IDENTICAL to @ref assign_row's coverage. All four
   * per-field consumers MUST be kept in sync: @ref assign_row, this
   * @ref copy_row, @ref permute_rebuild, and the MigrationPack per-field wire
   * pack (particle_store/MigrationPack.cpp). Any field added to one must be
   * added to the others. The maximal-population round-trip unit test enforces
   * this.
   */
  void copy_row(ParticleStore const &source, int source_row,
                int destination_row);

  /**
   * @brief Rebuild the store by permuting surviving rows.
   *
   * The resort-as-permutation path: contiguous, vectorizable per-column permute
   * kernels (`col_new[new_row] = col_old[perm[new_row]]`). Like
   * @ref begin_rebuild it swaps the current and spare generations (the retired
   * generation becomes the permute READ source) and grows the write target to
   * @p n_local + @p n_ghost rows; it then, for every column / POD sidecar /
   * ragged sidecar, moves the data from the old row named by @p permutation
   * into the new row in one pass per column.
   *
   * @p permutation has exactly @p n_local + @p n_ghost entries, one per new
   * row. `permutation[new_row] >= 0` names the retired-generation row whose
   * data is moved into @c new_row (a SURVIVING row); it must be a valid index
   * in the retired generation. `permutation[new_row] < 0` marks a NON-survivor
   * row -- a staged / brand-new / fresh-ghost row that carries no
   * retired-generation data: it is seeded to the new-particle defaults here
   * (@ref seed_default_row). The caller overwrites a staged local afterwards
   * via @ref copy_row from the source (staging) store; the ghost tail is left
   * at defaults (ghost exchange re-seeds it). The ghost tail is therefore
   * freshly seeded to defaults by this call whenever its permutation entries
   * are negative, matching the fresh-ghost contract of the assign_row rebuild.
   *
   * The permuted field set is IDENTICAL to @ref assign_row / @ref copy_row (the
   * four-way sync note above @ref copy_row applies); the maximal-population
   * permute unit tests enforce it. On @ref finish_rebuild the generation is
   * bumped exactly as for an assign_row rebuild.
   *
   * This is the live resort rebuild path:
   * @ref CellStructure::ensure_particle_store_synchronized builds the
   * permutation and calls this.
   */
  void permute_rebuild(std::span<int const> permutation, std::size_t n_local,
                       std::size_t n_ghost);

  /**
   * @brief Construct a non-owning view @ref Particle bound to a store row.
   *
   * Returns a @ref Particle whose accessors read/write the columns of THIS
   * store at @p row (via the attach handle).  This is the primary view
   * factory: @ref RowParticleRange / @ref RowParticleIterator call it for
   * every element produced by @ref Cell::particles(), and decomposition
   * resorts, @c for_each, and reduction kernels all consume views built here.
   *
   * The row must be a currently-valid index; the returned view aliases the
   * store and is invalidated by the next rebuild (which may renumber or drop
   * the row).
   */
  Particle make_view(int row);

  /**
   * @brief Seed @p row with the default new-particle values.
   *
   * Writes the default value of every column / POD sidecar (and clears the
   * ragged bond/exclusion sidecars) at @p row -- the values a genuinely-new /
   * fresh-ghost particle starts from (id -1, mol_id 0, type 0, propagation
   * SYSTEM_DEFAULT, position/velocity/force zero, quaternion identity, mass 1,
   * gamma -1, etc.). These must match @ref assign_row's seed branch. The
   * new-particle creation path (@ref CellStructure::add_particle via a staging
   * row) seeds a row, then the caller writes the fields it wants through a
   * view;
   * @ref CellParticleStorage::resize_ghost_storage seeds fresh ghost rows.
   */
  void seed_default_row(int row);

  /** Release all Kokkos-backed columns. Must be called while the Kokkos
   *  runtime is still alive (e.g. before Kokkos::finalize); afterwards the
   *  store is empty and dirty. */
  void release_columns() {
    m_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_torque = Column{};
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_rattle_correction = Column{};
#endif
    m_position = Column{};
    m_image_box = IntegerColumn{};
#ifdef ESPRESSO_ROTATION
    m_quaternion = QuaternionColumn{};
#endif
    m_position_at_last_verlet_update = Column{};
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_position_last_time_step = Column{};
#endif
    m_lees_edwards_offset = ScalarColumn{};
    m_lees_edwards_flag = ShortColumn{};
    m_velocity = Column{};
#ifdef ESPRESSO_ROTATION
    m_angular_velocity = Column{};
#endif
    m_old_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_old_torque = Column{};
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_old_rattle_correction = Column{};
#endif
    m_old_position = Column{};
    m_old_image_box = IntegerColumn{};
#ifdef ESPRESSO_ROTATION
    m_old_quaternion = QuaternionColumn{};
#endif
    m_old_position_at_last_verlet_update = Column{};
#ifdef ESPRESSO_BOND_CONSTRAINT
    m_old_position_last_time_step = Column{};
#endif
    m_old_lees_edwards_offset = ScalarColumn{};
    m_old_lees_edwards_flag = ShortColumn{};
    m_old_velocity = Column{};
#ifdef ESPRESSO_ROTATION
    m_old_angular_velocity = Column{};
#endif
    // Parameter columns.
    m_id = IntScalarColumn{};
    m_mol_id = IntScalarColumn{};
    m_type = IntScalarColumn{};
    m_propagation = IntScalarColumn{};
#ifdef ESPRESSO_ROTATION
    m_rotation = Uint8Column{};
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    m_ext_flag = Uint8Column{};
#endif
#ifdef ESPRESSO_MASS
    m_mass = ScalarColumn{};
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    m_q = ScalarColumn{};
#endif
#ifdef ESPRESSO_DIPOLES
    m_dipm = ScalarColumn{};
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    m_rinertia = Column{};
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    m_mu_E = Column{};
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    m_dip_fld = Column{};
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    m_ext_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_ext_torque = Column{};
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    m_gamma = GammaColumn{};
#ifdef ESPRESSO_ROTATION
    m_gamma_rot = GammaColumn{};
#endif
#endif
    // Parameter spares.
    m_old_id = IntScalarColumn{};
    m_old_mol_id = IntScalarColumn{};
    m_old_type = IntScalarColumn{};
    m_old_propagation = IntScalarColumn{};
#ifdef ESPRESSO_ROTATION
    m_old_rotation = Uint8Column{};
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    m_old_ext_flag = Uint8Column{};
#endif
#ifdef ESPRESSO_MASS
    m_old_mass = ScalarColumn{};
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    m_old_q = ScalarColumn{};
#endif
#ifdef ESPRESSO_DIPOLES
    m_old_dipm = ScalarColumn{};
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    m_old_rinertia = Column{};
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    m_old_mu_E = Column{};
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    m_old_dip_fld = Column{};
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    m_old_ext_force = Column{};
#ifdef ESPRESSO_ROTATION
    m_old_ext_torque = Column{};
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    m_old_gamma = GammaColumn{};
#ifdef ESPRESSO_ROTATION
    m_old_gamma_rot = GammaColumn{};
#endif
#endif
    // Host parameter sidecars.
#ifdef ESPRESSO_ENGINE
    m_swimming.clear();
    m_old_swimming.clear();
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    m_magnetodynamics.clear();
    m_old_magnetodynamics.clear();
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    m_vs_relative.clear();
    m_old_vs_relative.clear();
#endif
    // Ragged host sidecars.
    m_bonds_sidecar.clear();
    m_old_bonds_sidecar.clear();
#ifdef ESPRESSO_EXCLUSIONS
    m_exclusions_sidecar.clear();
    m_old_exclusions_sidecar.clear();
#endif
    m_number_of_local_particles = 0u;
    m_number_of_ghost_particles = 0u;
    m_old_number_of_particles = 0u;
    m_pending_removal.clear();
    m_pending_removal_count = 0u;
    m_dirty = true;
  }

  // -- observable columns ---------------------------------------------------
  VectorReference force_reference(int const row) {
    return column_reference(m_force, row);
  }
  Utils::Vector3d force_value(int const row) const {
    return column_value(m_force, row);
  }
#ifdef ESPRESSO_ROTATION
  VectorReference torque_reference(int const row) {
    return column_reference(m_torque, row);
  }
  Utils::Vector3d torque_value(int const row) const {
    return column_value(m_torque, row);
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  // RATTLE/SHAKE correction accumulator: a per-iteration scratch observable,
  // structurally like force -- zeroed each SHAKE iteration, reduced into
  // locals, then applied. Never persisted nor migrated; a genuinely-new row
  // defaults to zero (preserve-or-default like dip_fld / force).
  VectorReference rattle_correction_reference(int const row) {
    return column_reference(m_rattle_correction, row);
  }
  Utils::Vector3d rattle_correction_value(int const row) const {
    return column_value(m_rattle_correction, row);
  }
  auto rattle_correction_view() { return m_rattle_correction.view_host(); }
#endif

  // -- state columns --------------------------------------------------------
  VectorReference position_reference(int const row) {
    return column_reference(m_position, row);
  }
  Utils::Vector3d position_value(int const row) const {
    return column_value(m_position, row);
  }
  IntegerVectorReference image_box_reference(int const row) {
    return integer_column_reference(m_image_box, row);
  }
  Utils::Vector3i image_box_value(int const row) const {
    return integer_column_value(m_image_box, row);
  }
#ifdef ESPRESSO_ROTATION
  QuaternionReference quaternion_reference(int const row) {
    return quaternion_column_reference(m_quaternion, row);
  }
  Utils::Quaternion<double> quaternion_value(int const row) const {
    return quaternion_column_value(m_quaternion, row);
  }
#endif
  VectorReference position_at_last_verlet_update_reference(int const row) {
    return column_reference(m_position_at_last_verlet_update, row);
  }
  Utils::Vector3d position_at_last_verlet_update_value(int const row) const {
    return column_value(m_position_at_last_verlet_update, row);
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  VectorReference position_last_time_step_reference(int const row) {
    return column_reference(m_position_last_time_step, row);
  }
  Utils::Vector3d position_last_time_step_value(int const row) const {
    return column_value(m_position_last_time_step, row);
  }
#endif
  double &lees_edwards_offset(int const row) {
    return scalar_reference(m_lees_edwards_offset, row);
  }
  double lees_edwards_offset(int const row) const {
    return scalar_value(m_lees_edwards_offset, row);
  }
  short &lees_edwards_flag(int const row) {
    return scalar_reference(m_lees_edwards_flag, row);
  }
  short lees_edwards_flag(int const row) const {
    return scalar_value(m_lees_edwards_flag, row);
  }

  // -- host-view getters for kernel wiring ----------------------------------
  auto force_view() { return m_force.view_host(); }
#ifdef ESPRESSO_ROTATION
  auto torque_view() { return m_torque.view_host(); }
#endif
  auto position_view() { return m_position.view_host(); }
  auto image_box_view() { return m_image_box.view_host(); }
#ifdef ESPRESSO_ROTATION
  auto quaternion_view() { return m_quaternion.view_host(); }
#endif
  auto position_at_last_verlet_update_view() {
    return m_position_at_last_verlet_update.view_host();
  }
#ifdef ESPRESSO_BOND_CONSTRAINT
  auto position_last_time_step_view() {
    return m_position_last_time_step.view_host();
  }
#endif
  auto lees_edwards_offset_view() { return m_lees_edwards_offset.view_host(); }
  auto lees_edwards_flag_view() { return m_lees_edwards_flag.view_host(); }

  // -- device mirrors + host<->device sync (GPU-core migration, phase 1) ------
  // Device-space views of the hot state columns, for device-resident kernels;
  // the sync helpers push/pull them explicitly at integration-loop boundaries.
  // On a host-only Kokkos build the device mirror aliases the host mirror, so
  // the syncs are self-copies; these are only exercised on the (opt-in) device
  // execution path, so the host path is unaffected.
  auto position_view_device() { return m_position.view_device(); }
  auto velocity_view_device() { return m_velocity.view_device(); }
  auto force_view_device() { return m_force.view_device(); }
  auto image_box_view_device() { return m_image_box.view_device(); }
  auto id_view_device() { return m_id.view_device(); }
  auto type_view_device() { return m_type.view_device(); }
#ifdef ESPRESSO_ELECTROSTATICS
  auto q_view_device() { return m_q.view_device(); }
#endif
#ifdef ESPRESSO_MASS
  auto mass_view_device() { return m_mass.view_device(); }
#endif

  /** @brief Push the hot state columns to the device mirror. */
  void sync_state_to_device() {
    Kokkos::deep_copy(m_position.view_device(), m_position.view_host());
    Kokkos::deep_copy(m_velocity.view_device(), m_velocity.view_host());
    Kokkos::deep_copy(m_force.view_device(), m_force.view_host());
    Kokkos::deep_copy(m_image_box.view_device(), m_image_box.view_host());
    Kokkos::deep_copy(m_id.view_device(), m_id.view_host());
    Kokkos::deep_copy(m_type.view_device(), m_type.view_host());
#ifdef ESPRESSO_ELECTROSTATICS
    Kokkos::deep_copy(m_q.view_device(), m_q.view_host());
#endif
#ifdef ESPRESSO_MASS
    Kokkos::deep_copy(m_mass.view_device(), m_mass.view_host());
#endif
  }
  /** @brief Pull the device-written state columns back to the host mirror. */
  void sync_state_to_host() {
    Kokkos::deep_copy(m_position.view_host(), m_position.view_device());
    Kokkos::deep_copy(m_velocity.view_host(), m_velocity.view_device());
    Kokkos::deep_copy(m_force.view_host(), m_force.view_device());
    Kokkos::deep_copy(m_image_box.view_host(), m_image_box.view_device());
  }

  // -- momentum columns -----------------------------------------------------
  VectorReference velocity_reference(int const row) {
    return column_reference(m_velocity, row);
  }
  Utils::Vector3d velocity_value(int const row) const {
    return column_value(m_velocity, row);
  }
  auto velocity_view() { return m_velocity.view_host(); }
#ifdef ESPRESSO_ROTATION
  VectorReference angular_velocity_reference(int const row) {
    return column_reference(m_angular_velocity, row);
  }
  Utils::Vector3d angular_velocity_value(int const row) const {
    return column_value(m_angular_velocity, row);
  }
  auto angular_velocity_view() { return m_angular_velocity.view_host(); }
#endif

  // -- parameter columns ----------------------------------------------------
  // Scalar-typed parameters expose a real element reference (like
  // lees_edwards_offset) plus a by-value getter and a host view. Vector-typed
  // parameters expose a VectorReference/value/view triple (like force).

  // int scalars: id, mol_id, type, propagation.
  int &id(int const row) { return scalar_reference(m_id, row); }
  int id(int const row) const { return scalar_value(m_id, row); }
  auto id_view() { return m_id.view_host(); }
  int &mol_id(int const row) { return scalar_reference(m_mol_id, row); }
  int mol_id(int const row) const { return scalar_value(m_mol_id, row); }
  auto mol_id_view() { return m_mol_id.view_host(); }
  int &type(int const row) { return scalar_reference(m_type, row); }
  int type(int const row) const { return scalar_value(m_type, row); }
  auto type_view() { return m_type.view_host(); }
  int &propagation(int const row) {
    return scalar_reference(m_propagation, row);
  }
  int propagation(int const row) const {
    return scalar_value(m_propagation, row);
  }
  auto propagation_view() { return m_propagation.view_host(); }

#ifdef ESPRESSO_ROTATION
  std::uint8_t &rotation(int const row) {
    return scalar_reference(m_rotation, row);
  }
  std::uint8_t rotation(int const row) const {
    return scalar_value(m_rotation, row);
  }
  auto rotation_view() { return m_rotation.view_host(); }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t &ext_flag(int const row) {
    return scalar_reference(m_ext_flag, row);
  }
  std::uint8_t ext_flag(int const row) const {
    return scalar_value(m_ext_flag, row);
  }
  auto ext_flag_view() { return m_ext_flag.view_host(); }
#endif

#ifdef ESPRESSO_MASS
  double &mass(int const row) { return scalar_reference(m_mass, row); }
  double mass(int const row) const { return scalar_value(m_mass, row); }
  auto mass_view() { return m_mass.view_host(); }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double &q(int const row) { return scalar_reference(m_q, row); }
  double q(int const row) const { return scalar_value(m_q, row); }
  auto q_view() { return m_q.view_host(); }
#endif
#ifdef ESPRESSO_DIPOLES
  double &dipm(int const row) { return scalar_reference(m_dipm, row); }
  double dipm(int const row) const { return scalar_value(m_dipm, row); }
  auto dipm_view() { return m_dipm.view_host(); }
#endif

#ifdef ESPRESSO_ROTATIONAL_INERTIA
  VectorReference rinertia_reference(int const row) {
    return column_reference(m_rinertia, row);
  }
  Utils::Vector3d rinertia_value(int const row) const {
    return column_value(m_rinertia, row);
  }
  auto rinertia_view() { return m_rinertia.view_host(); }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  VectorReference mu_E_reference(int const row) {
    return column_reference(m_mu_E, row);
  }
  Utils::Vector3d mu_E_value(int const row) const {
    return column_value(m_mu_E, row);
  }
  auto mu_E_view() { return m_mu_E.view_host(); }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  VectorReference dip_fld_reference(int const row) {
    return column_reference(m_dip_fld, row);
  }
  Utils::Vector3d dip_fld_value(int const row) const {
    return column_value(m_dip_fld, row);
  }
  auto dip_fld_view() { return m_dip_fld.view_host(); }
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  VectorReference ext_force_reference(int const row) {
    return column_reference(m_ext_force, row);
  }
  Utils::Vector3d ext_force_value(int const row) const {
    return column_value(m_ext_force, row);
  }
  auto ext_force_view() { return m_ext_force.view_host(); }
#ifdef ESPRESSO_ROTATION
  VectorReference ext_torque_reference(int const row) {
    return column_reference(m_ext_torque, row);
  }
  Utils::Vector3d ext_torque_value(int const row) const {
    return column_value(m_ext_torque, row);
  }
  auto ext_torque_view() { return m_ext_torque.view_host(); }
#endif
#endif

#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  // gamma/gamma_rot: scalar element-ref when isotropic, VectorReference when
  // ESPRESSO_PARTICLE_ANISOTROPY selects per-axis friction (GammaColumn).
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  VectorReference gamma_reference(int const row) {
    return column_reference(m_gamma, row);
  }
  GammaValue gamma_value(int const row) const {
    return column_value(m_gamma, row);
  }
#else
  double &gamma_reference(int const row) {
    return scalar_reference(m_gamma, row);
  }
  GammaValue gamma_value(int const row) const {
    return scalar_value(m_gamma, row);
  }
#endif
  auto gamma_view() { return m_gamma.view_host(); }
#ifdef ESPRESSO_ROTATION
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  VectorReference gamma_rot_reference(int const row) {
    return column_reference(m_gamma_rot, row);
  }
  GammaValue gamma_rot_value(int const row) const {
    return column_value(m_gamma_rot, row);
  }
#else
  double &gamma_rot_reference(int const row) {
    return scalar_reference(m_gamma_rot, row);
  }
  GammaValue gamma_rot_value(int const row) const {
    return scalar_value(m_gamma_rot, row);
  }
#endif
  auto gamma_rot_view() { return m_gamma_rot.view_host(); }
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

  // -- typed zero-extent dummy views ----------------------------------------
  // When an optional feature (rotation, mass, per-particle friction, ...) is
  // compiled OUT, the integrator kernels still need to NAME a view handle for
  // that feature's column so the shared lambda signature stays fixed; the
  // handle is only ever passed to code paths that are themselves #ifdef'd out.
  // These helpers hand back a zero-extent host view of the CORRECT column type
  // (rather than aliasing an unrelated live column), so naming is safe but any
  // index into one traps under bounds checking (ESPRESSO_ADDITIONAL_CHECKS /
  // Kokkos debug) instead of silently reading the wrong column type. A
  // default-constructed column has extent 0; its host mirror is the dummy.
  auto dummy_vector_view() { return Column{}.view_host(); }
  auto dummy_quaternion_view() { return QuaternionColumn{}.view_host(); }
  auto dummy_scalar_view() { return ScalarColumn{}.view_host(); }
  auto dummy_uint8_view() { return Uint8Column{}.view_host(); }
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  auto dummy_gamma_view() { return GammaColumn{}.view_host(); }
#endif

  // -- host parameter sidecars ----------------------------------------------
  // Cold PODs live in plain std::vector sidecars indexed by store row (not
  // Kokkos columns). Rebuilt with the store (old vector swapped like a column;
  // preserve-by-old-row / seed-default in assign_row). Accessors return
  // references by row.
#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming &swimming(int const row) {
    return sidecar_reference(m_swimming, row);
  }
  ParticleParametersSwimming const &swimming(int const row) const {
    return sidecar_reference(m_swimming, row);
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  ThermalStonerWohlfarthParameters &magnetodynamics(int const row) {
    return sidecar_reference(m_magnetodynamics, row);
  }
  ThermalStonerWohlfarthParameters const &magnetodynamics(int const row) const {
    return sidecar_reference(m_magnetodynamics, row);
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  VirtualSitesRelativeParameters &vs_relative(int const row) {
    return sidecar_reference(m_vs_relative, row);
  }
  VirtualSitesRelativeParameters const &vs_relative(int const row) const {
    return sidecar_reference(m_vs_relative, row);
  }
#endif

  // -- ragged host sidecars -------------------------------------------------
  // Owned per-particle variable-size data (bond list; exclusion list) lives in
  // plain std::vector sidecars indexed by store row, following the POD sidecar
  // machinery (swap/resize in begin_rebuild; preserve-or-seed in assign_row;
  // cleared in release_columns) with ONE difference: the element owns heap
  // storage, so a surviving row is MOVED from the old vector element rather
  // than copied. Accessors return a real element reference (same as the POD
  // sidecars) via the shared sidecar_reference helper.
  BondList &bonds_sidecar_reference(int const row) {
    return sidecar_reference(m_bonds_sidecar, row);
  }
  BondList const &bonds_sidecar_reference(int const row) const {
    return sidecar_reference(m_bonds_sidecar, row);
  }
#ifdef ESPRESSO_EXCLUSIONS
  Utils::compact_vector<int> &exclusions_sidecar_reference(int const row) {
    return sidecar_reference(m_exclusions_sidecar, row);
  }
  Utils::compact_vector<int> const &
  exclusions_sidecar_reference(int const row) const {
    return sidecar_reference(m_exclusions_sidecar, row);
  }
#endif

private:
  // Swap current <-> spare generations and grow every column / sidecar write
  // target to hold @p n_local + @p n_ghost rows. Shared by @ref begin_rebuild
  // (assign_row loop) and @ref permute_rebuild (per-column permute); both leave
  // the retired generation in the m_old_* buffers as the preserve/permute read
  // source and overwrite every logical row before finish_rebuild.
  void swap_and_grow_generations(std::size_t n_local, std::size_t n_ghost);

  // -- proxy factories for the various column kinds -------------------------
  VectorReference column_reference(Column &column, int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    // Bind by const reference (never by non-const reference): view_host()
    // returns `const t_host&` in the vendored Kokkos but `t_host` by value when
    // KOKKOS_ENABLE_DEPRECATED_CODE_4 is on (Kokkos <=4.5 / fedora / macos),
    // and a non-const lvalue reference cannot bind to that by-value rvalue. A
    // Kokkos View has shallow-const semantics, so `operator()` on a const-ref'd
    // view still yields a writable element; the write-through proxy below
    // remains valid. Binding by const reference (not by value) also avoids a
    // View handle copy -- and its reference-count atomic -- on the default
    // (non-deprecated) build, where this per-particle accessor is a hot path.
    auto const &view = column.view_host();
    // Layout-agnostic proxy over one row: address of component 0 plus the
    // component-to-component stride. stride(1) is 1 for LayoutRight (row
    // contiguous) and number-of-rows for LayoutLeft; &view(row, 0) is the
    // correct base under either layout, unlike view.data() + row (which
    // assumes LayoutLeft). stride(1) is the non-deprecated spelling of
    // stride_1() in the installed Kokkos version.
    return {&view(row, 0), view.stride(1)};
  }
  Utils::Vector3d column_value(Column const &column, int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    return {view(row, 0), view(row, 1), view(row, 2)};
  }
  IntegerVectorReference integer_column_reference(IntegerColumn &column,
                                                  int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    // Const-ref binding: see column_reference for the view_host() by-value vs
    // by-reference compatibility rationale and the shallow-const write-through.
    auto const &view = column.view_host();
    return {&view(row, 0), view.stride(1)};
  }
  Utils::Vector3i integer_column_value(IntegerColumn const &column,
                                       int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    return {view(row, 0), view(row, 1), view(row, 2)};
  }
  QuaternionReference quaternion_column_reference(QuaternionColumn &column,
                                                  int const row) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    // Const-ref binding: see column_reference for the view_host() by-value vs
    // by-reference compatibility rationale and the shallow-const write-through.
    auto const &view = column.view_host();
    return {&view(row, 0), view.stride(1)};
  }
  Utils::Quaternion<double>
  quaternion_column_value(QuaternionColumn const &column, int const row) const {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    Utils::Quaternion<double> q;
    q[0] = view(row, 0);
    q[1] = view(row, 1);
    q[2] = view(row, 2);
    q[3] = view(row, 3);
    return q;
  }
  template <class ScalarColumnType>
  auto scalar_reference(ScalarColumnType &column, int const row)
      -> decltype(column.view_host()(row)) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    return column.view_host()(row);
  }
  template <class ScalarColumnType>
  auto scalar_value(ScalarColumnType const &column, int const row) const
      -> std::remove_reference_t<decltype(column.view_host()(row))> {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    return column.view_host()(row);
  }
  // Host sidecar (plain std::vector) element access by row. The return type is
  // deduced from the operator[] on the (possibly const) SidecarVector argument,
  // so a const sidecar yields `POD const &` and a non-const one yields `POD &`;
  // this single helper therefore backs BOTH the const and non-const public
  // sidecar accessors with the correct element constness (the method itself is
  // const because it only reads the vector handle, never mutates the store).
  template <class SidecarVector>
  auto sidecar_reference(SidecarVector &sidecar, int const row) const
      -> decltype(sidecar[static_cast<std::size_t>(row)]) {
    assert(row >= 0 and static_cast<std::size_t>(row) < number_of_particles());
    return sidecar[static_cast<std::size_t>(row)];
  }

  std::size_t m_number_of_local_particles = 0u;
  std::size_t m_number_of_ghost_particles = 0u;
  bool m_dirty = false;
  std::uint64_t m_generation = 0u;

  // -- mid-epoch pending-removal mask ---------------------------------------
  // One byte per store row; 1 = row was dropped from its cell this epoch and is
  // invisible to iteration until the next rebuild resolves it. Sized to
  // number_of_particles() and zeroed by every rebuild (see the accessors and
  // begin_rebuild / permute_rebuild).
  std::vector<char> m_pending_removal;
  /** Number of rows currently marked pending-removal (backs the O(1)
   *  @ref has_pending_removals fast path). */
  std::size_t m_pending_removal_count = 0u;

  // -- current-generation columns -------------------------------------------
  Column m_force;
#ifdef ESPRESSO_ROTATION
  Column m_torque;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  Column m_rattle_correction;
#endif
  Column m_position;
  IntegerColumn m_image_box;
#ifdef ESPRESSO_ROTATION
  QuaternionColumn m_quaternion;
#endif
  Column m_position_at_last_verlet_update;
#ifdef ESPRESSO_BOND_CONSTRAINT
  Column m_position_last_time_step;
#endif
  ScalarColumn m_lees_edwards_offset;
  ShortColumn m_lees_edwards_flag;
  Column m_velocity;
#ifdef ESPRESSO_ROTATION
  Column m_angular_velocity;
#endif

  // -- current-generation parameter columns ---------------------------------
  IntScalarColumn m_id;
  IntScalarColumn m_mol_id;
  IntScalarColumn m_type;
  IntScalarColumn m_propagation;
#ifdef ESPRESSO_ROTATION
  Uint8Column m_rotation;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Uint8Column m_ext_flag;
#endif
#ifdef ESPRESSO_MASS
  ScalarColumn m_mass;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  ScalarColumn m_q;
#endif
#ifdef ESPRESSO_DIPOLES
  ScalarColumn m_dipm;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Column m_rinertia;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Column m_mu_E;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Column m_dip_fld;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Column m_ext_force;
#ifdef ESPRESSO_ROTATION
  Column m_ext_torque;
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  GammaColumn m_gamma;
#ifdef ESPRESSO_ROTATION
  GammaColumn m_gamma_rot;
#endif
#endif

  // -- host parameter sidecars ----------------------------------------------
  // Cold PODs (indexed by store row). Rebuilt with the store: the current
  // vector is swapped with its spare in begin_rebuild (holding the old-row
  // values as the preserve source), grown to the new count, and each row is
  // preserved-by-old-row / seeded-to-default in assign_row.
#ifdef ESPRESSO_ENGINE
  std::vector<ParticleParametersSwimming> m_swimming;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  std::vector<ThermalStonerWohlfarthParameters> m_magnetodynamics;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  std::vector<VirtualSitesRelativeParameters> m_vs_relative;
#endif

  // -- ragged host sidecars -------------------------------------------------
  // Owned variable-size per-particle data (indexed by store row). Rebuilt with
  // the store like the POD sidecars (swap current <-> spare in begin_rebuild,
  // resize to the new count, preserve-by-old-row / seed-default in assign_row)
  // -- but the preserve step MOVES the element out of the old vector (the
  // element owns heap storage) instead of copying it. Ghost/new rows default
  // to empty. The exclusion sidecar exists only when EXCLUSIONS is enabled (its
  // element type is the exact type of Particle::el).
  std::vector<BondList> m_bonds_sidecar;
#ifdef ESPRESSO_EXCLUSIONS
  std::vector<Utils::compact_vector<int>> m_exclusions_sidecar;
#endif

  // -- spare (previous-generation) columns ----------------------------------
  // Capacity-cached double buffering: these are kept alive across rebuilds as
  // the swap-in write target. During a rebuild they hold the just-current data
  // (the preserve_or_seed READ source); after finish_rebuild they hold the
  // retired generation, reused as the spare on the next swap.
  Column m_old_force;
#ifdef ESPRESSO_ROTATION
  Column m_old_torque;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  Column m_old_rattle_correction;
#endif
  Column m_old_position;
  IntegerColumn m_old_image_box;
#ifdef ESPRESSO_ROTATION
  QuaternionColumn m_old_quaternion;
#endif
  Column m_old_position_at_last_verlet_update;
#ifdef ESPRESSO_BOND_CONSTRAINT
  Column m_old_position_last_time_step;
#endif
  ScalarColumn m_old_lees_edwards_offset;
  ShortColumn m_old_lees_edwards_flag;
  Column m_old_velocity;
#ifdef ESPRESSO_ROTATION
  Column m_old_angular_velocity;
#endif

  // -- spare (previous-generation) parameter columns ------------------------
  IntScalarColumn m_old_id;
  IntScalarColumn m_old_mol_id;
  IntScalarColumn m_old_type;
  IntScalarColumn m_old_propagation;
#ifdef ESPRESSO_ROTATION
  Uint8Column m_old_rotation;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Uint8Column m_old_ext_flag;
#endif
#ifdef ESPRESSO_MASS
  ScalarColumn m_old_mass;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  ScalarColumn m_old_q;
#endif
#ifdef ESPRESSO_DIPOLES
  ScalarColumn m_old_dipm;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Column m_old_rinertia;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Column m_old_mu_E;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Column m_old_dip_fld;
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Column m_old_ext_force;
#ifdef ESPRESSO_ROTATION
  Column m_old_ext_torque;
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  GammaColumn m_old_gamma;
#ifdef ESPRESSO_ROTATION
  GammaColumn m_old_gamma_rot;
#endif
#endif

  // -- spare (previous-generation) host sidecars ----------------------------
#ifdef ESPRESSO_ENGINE
  std::vector<ParticleParametersSwimming> m_old_swimming;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  std::vector<ThermalStonerWohlfarthParameters> m_old_magnetodynamics;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  std::vector<VirtualSitesRelativeParameters> m_old_vs_relative;
#endif

  // -- spare (previous-generation) ragged host sidecars ---------------------
  std::vector<BondList> m_old_bonds_sidecar;
#ifdef ESPRESSO_EXCLUSIONS
  std::vector<Utils::compact_vector<int>> m_old_exclusions_sidecar;
#endif

  std::size_t m_old_number_of_particles = 0u;
};

/*
 * Copyright (C) 2017-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE ParticleStore test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "ParticleStoreMaximalPopulation.hpp"

#include "particle_store/ParticleStore.hpp"

#include "BondList.hpp"
#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/quaternion.hpp>

#include <Kokkos_Core.hpp>

#include <array>
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

BOOST_AUTO_TEST_CASE(default_constructed_store_is_empty) {
  ParticleStore const store{};
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 0ul);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 0ul);
}

BOOST_AUTO_TEST_CASE(rebuild_assigns_rows_and_zero_initializes) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.mark_dirty();
  BOOST_CHECK(store.is_dirty());
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();
  BOOST_CHECK(not store.is_dirty());
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 2ul);
  BOOST_CHECK_EQUAL(p0.store_row(), 0);
  BOOST_CHECK_EQUAL(p1.store_row(), 1);
  Utils::Vector3d const f0 = store.force_value(0);
  BOOST_CHECK_EQUAL(f0.norm2(), 0.);
}

BOOST_AUTO_TEST_CASE(rebuild_preserves_values_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();
  store.force_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.force_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};

  // simulate a resort that swaps the two particles' order and adds one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  Utils::Vector3d const f_p1 = store.force_value(p1.store_row());
  Utils::Vector3d const f_p0 = store.force_value(p0.store_row());
  Utils::Vector3d const f_p2 = store.force_value(p2.store_row());
  BOOST_CHECK_EQUAL(f_p1[0], 4.);      // p1 kept its values at its new row
  BOOST_CHECK_EQUAL(f_p0[2], 3.);      // p0 kept its values at its new row
  BOOST_CHECK_EQUAL(f_p2.norm2(), 0.); // new particle zero-initialized
}

BOOST_AUTO_TEST_CASE(ghost_rows_follow_locals) {
  ParticleStore store{};
  Particle p_local{}, p_ghost{};
  store.begin_rebuild(1u, 1u);
  store.assign_row(p_local, 0);
  store.assign_row(p_ghost, 1);
  store.finish_rebuild();
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 1ul);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 1ul);
}

// Force/torque live in the store columns. A particle that migrates to another
// rank carries its force across via a row-to-row ParticleStore::copy_row from
// the sending rank's store into a fresh row on the receiving rank's store.
// This models that copy and asserts the force (and torque) survive it.
BOOST_AUTO_TEST_CASE(rebuild_seeds_migrated_particle_force_from_carrier) {
  // 1) A row in a source store with a known force (models the sending rank
  //    right before a global resort).
  ParticleStore source{};
  source.begin_rebuild(1u, 0u);
  source.seed_default_row(0);
  source.finish_rebuild();
  auto const f_ref = Utils::Vector3d{-1.5, 2.25, -3.75};
  source.force_reference(0) = f_ref;
#ifdef ESPRESSO_ROTATION
  auto const t_ref = Utils::Vector3d{4.5, -5.5, 6.5};
  source.torque_reference(0) = t_ref;
#endif

  // 2) The receiving rank's store copies the source row into a fresh row, so
  //    the force survives the migration.
  ParticleStore target{};
  target.begin_rebuild(1u, 0u);
  target.copy_row(source, 0, 0);
  target.finish_rebuild();
  auto const f_new = target.force_value(0);
  BOOST_CHECK_EQUAL(f_new[0], f_ref[0]);
  BOOST_CHECK_EQUAL(f_new[1], f_ref[1]);
  BOOST_CHECK_EQUAL(f_new[2], f_ref[2]);
#ifdef ESPRESSO_ROTATION
  auto const t_new = target.torque_value(0);
  BOOST_CHECK_EQUAL(t_new[0], t_ref[0]);
  BOOST_CHECK_EQUAL(t_new[2], t_ref[2]);
#endif
}

// State columns: a rank-local rebuild that shuffles the row order must
// preserve each particle's state (position, image box, quaternion, and the
// Lees-Edwards offset) by old row, exactly as the force column does.
// Velocity columns are verified alongside.
BOOST_AUTO_TEST_CASE(rebuild_preserves_state_columns_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.position_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.position_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};
  store.image_box_reference(p0.store_row()) = Utils::Vector3i{-1, 0, 2};
  store.image_box_reference(p1.store_row()) = Utils::Vector3i{7, -8, 9};
  store.position_at_last_verlet_update_reference(p0.store_row()) =
      Utils::Vector3d{0.1, 0.2, 0.3};
  store.position_at_last_verlet_update_reference(p1.store_row()) =
      Utils::Vector3d{0.4, 0.5, 0.6};
  store.lees_edwards_offset(p0.store_row()) = 1.25;
  store.lees_edwards_offset(p1.store_row()) = -2.5;
  store.lees_edwards_flag(p0.store_row()) = static_cast<short>(3);
  store.lees_edwards_flag(p1.store_row()) = static_cast<short>(-4);
  store.velocity_reference(p0.store_row()) = Utils::Vector3d{1.1, 2.2, 3.3};
  store.velocity_reference(p1.store_row()) = Utils::Vector3d{-4.4, 5.5, -6.6};
#ifdef ESPRESSO_ROTATION
  auto const q0 = Utils::Quaternion<double>{{0.5, 0.5, 0.5, 0.5}};
  auto const q1 = Utils::Quaternion<double>{{0., 1., 0., 0.}};
  store.quaternion_reference(p0.store_row()) = q0;
  store.quaternion_reference(p1.store_row()) = q1;
  store.angular_velocity_reference(p0.store_row()) =
      Utils::Vector3d{0.1, 0.2, 0.3};
  store.angular_velocity_reference(p1.store_row()) =
      Utils::Vector3d{-0.4, 0.5, -0.6};
#endif

  // resort that swaps the two particles' order and appends a new one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  Utils::Vector3d const pos_p1 = store.position_value(p1.store_row());
  Utils::Vector3d const pos_p0 = store.position_value(p0.store_row());
  BOOST_CHECK_EQUAL(pos_p1[0], 4.);
  BOOST_CHECK_EQUAL(pos_p0[2], 3.);
  Utils::Vector3i const img_p1 = store.image_box_value(p1.store_row());
  Utils::Vector3i const img_p0 = store.image_box_value(p0.store_row());
  BOOST_CHECK_EQUAL(img_p1[0], 7);
  BOOST_CHECK_EQUAL(img_p0[2], 2);
  Utils::Vector3d const pold_p1 =
      store.position_at_last_verlet_update_value(p1.store_row());
  BOOST_CHECK_EQUAL(pold_p1[1], 0.5);
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p1.store_row()), -2.5);
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p0.store_row()), 1.25);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p1.store_row()),
                    static_cast<short>(-4));
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p0.store_row()),
                    static_cast<short>(3));
  // Velocity preserved by old row; new row zeroed.
  Utils::Vector3d const vel_p1 = store.velocity_value(p1.store_row());
  Utils::Vector3d const vel_p0 = store.velocity_value(p0.store_row());
  Utils::Vector3d const vel_p2 = store.velocity_value(p2.store_row());
  BOOST_CHECK_EQUAL(vel_p1[0], -4.4);
  BOOST_CHECK_EQUAL(vel_p0[2], 3.3);
  BOOST_CHECK_EQUAL(vel_p2.norm2(), 0.);
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat_p1 =
      store.quaternion_value(p1.store_row());
  Utils::Quaternion<double> const quat_p0 =
      store.quaternion_value(p0.store_row());
  BOOST_CHECK_EQUAL(quat_p1[1], 1.);
  BOOST_CHECK_EQUAL(quat_p0[3], 0.5);
  // Angular velocity preserved by old row; new row zeroed.
  Utils::Vector3d const av_p1 = store.angular_velocity_value(p1.store_row());
  Utils::Vector3d const av_p0 = store.angular_velocity_value(p0.store_row());
  Utils::Vector3d const av_p2 = store.angular_velocity_value(p2.store_row());
  BOOST_CHECK_EQUAL(av_p1[1], 0.5);
  BOOST_CHECK_EQUAL(av_p0[2], 0.3);
  BOOST_CHECK_EQUAL(av_p2.norm2(), 0.);
#endif
}

// A genuinely new row (e.g. a fresh ghost, before its field update overwrites
// the column) must default to sensible values that match Particle's member
// defaults. The quaternion in particular must be the IDENTITY (1,0,0,0), NOT
// the Kokkos zero-init (0,0,0,0) which is an invalid quaternion.
// Velocity and angular velocity must default to zero.
BOOST_AUTO_TEST_CASE(new_row_state_defaults) {
  ParticleStore store{};
  Particle p{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  Utils::Vector3d const pos = store.position_value(p.store_row());
  BOOST_CHECK_EQUAL(pos.norm2(), 0.);
  Utils::Vector3i const img = store.image_box_value(p.store_row());
  BOOST_CHECK_EQUAL(img[0], 0);
  BOOST_CHECK_EQUAL(img[1], 0);
  BOOST_CHECK_EQUAL(img[2], 0);
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p.store_row()), 0.);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p.store_row()),
                    static_cast<short>(0));
  // Velocity defaults to zero.
  Utils::Vector3d const vel = store.velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(vel.norm2(), 0.);
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat = store.quaternion_value(p.store_row());
  BOOST_CHECK_EQUAL(quat[0], 1.); // identity, not zero
  BOOST_CHECK_EQUAL(quat[1], 0.);
  BOOST_CHECK_EQUAL(quat[2], 0.);
  BOOST_CHECK_EQUAL(quat[3], 0.);
  // Angular velocity defaults to zero.
  Utils::Vector3d const av = store.angular_velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(av.norm2(), 0.);
#endif
}

// A genuinely-new row's velocity column is seeded to the default (zero) by
// assign_row's not-preserve branch (via seed_default_row). This guards that
// assign_row seeds the velocity column (not skip it) and that the default is
// exactly zero (not garbage from WithoutInitializing).
BOOST_AUTO_TEST_CASE(rebuild_seeds_velocity_from_carrier) {
  // Fresh (detached) particle -> assign_row seeds the defaults.
  Particle p{};
  BOOST_REQUIRE(p.store() == nullptr);

  ParticleStore store{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  Utils::Vector3d const vel = store.velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(vel[0], 0.);
  BOOST_CHECK_EQUAL(vel[1], 0.);
  BOOST_CHECK_EQUAL(vel[2], 0.);
#ifdef ESPRESSO_ROTATION
  Utils::Vector3d const av = store.angular_velocity_value(p.store_row());
  BOOST_CHECK_EQUAL(av[0], 0.);
  BOOST_CHECK_EQUAL(av[1], 0.);
  BOOST_CHECK_EQUAL(av[2], 0.);
#endif
}

// A row that migrates in from another rank must have all its state columns
// carried across via a row-to-row ParticleStore::copy_row from the sending
// rank's store into a fresh row on the receiving rank's store. Set known state
// on a source row, copy it into a target store row, and confirm every state
// field survives the copy.
BOOST_AUTO_TEST_CASE(rebuild_seeds_migrated_particle_state_from_carrier) {
  ParticleStore source{};
  source.begin_rebuild(1u, 0u);
  source.seed_default_row(0);
  source.finish_rebuild();
  source.id(0) = 11;
  source.position_reference(0) = Utils::Vector3d{1.5, -2.5, 3.5};
  source.image_box_reference(0) = Utils::Vector3i{4, -5, 6};
  source.position_at_last_verlet_update_reference(0) =
      Utils::Vector3d{7.5, 8.5, 9.5};
  source.lees_edwards_offset(0) = 12.75;
  source.lees_edwards_flag(0) = static_cast<short>(2);
#ifdef ESPRESSO_ROTATION
  source.quaternion_reference(0) = Utils::Quaternion<double>{{0., 0., 1., 0.}};
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  source.position_last_time_step_reference(0) = Utils::Vector3d{-1., -2., -3.};
#endif

  ParticleStore target{};
  target.begin_rebuild(1u, 0u);
  target.copy_row(source, 0, 0);
  target.finish_rebuild();

  BOOST_CHECK_EQUAL(target.id(0), 11);
  Utils::Vector3d const pos = target.position_value(0);
  BOOST_CHECK_EQUAL(pos[0], 1.5);
  BOOST_CHECK_EQUAL(pos[1], -2.5);
  BOOST_CHECK_EQUAL(pos[2], 3.5);
  Utils::Vector3i const img = target.image_box_value(0);
  BOOST_CHECK_EQUAL(img[0], 4);
  BOOST_CHECK_EQUAL(img[2], 6);
  Utils::Vector3d const pold = target.position_at_last_verlet_update_value(0);
  BOOST_CHECK_EQUAL(pold[1], 8.5);
  BOOST_CHECK_EQUAL(target.lees_edwards_offset(0), 12.75);
  BOOST_CHECK_EQUAL(target.lees_edwards_flag(0), static_cast<short>(2));
#ifdef ESPRESSO_ROTATION
  Utils::Quaternion<double> const quat = target.quaternion_value(0);
  BOOST_CHECK_EQUAL(quat[2], 1.);
  BOOST_CHECK_EQUAL(quat[0], 0.);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  Utils::Vector3d const plast = target.position_last_time_step_value(0);
  BOOST_CHECK_EQUAL(plast[0], -1.);
  BOOST_CHECK_EQUAL(plast[2], -3.);
#endif
}

// Head-node fetch-cache SNAPSHOT-STORE pattern: a batch of rows fetched from a
// worker rank is materialized into a FIXED-capacity store built once per
// invalidation epoch, with rows handed out monotonically via row-to-row
// ParticleStore::copy_row. Build a source store holding a batch of known rows,
// copy each into a consecutive target-store row, and check each keeps its own
// data (never another row's). The head-node cache itself is exercised
// end-to-end by the multi-rank python gates.
BOOST_AUTO_TEST_CASE(snapshot_store_attaches_batch_of_detached_particles) {
  constexpr std::size_t capacity = 4u;

  // Source store with `capacity` rows of known, per-row data.
  ParticleStore source{};
  source.begin_rebuild(capacity, 0u);
  for (std::size_t i = 0u; i < capacity; ++i) {
    source.seed_default_row(static_cast<int>(i));
  }
  source.finish_rebuild();
  for (std::size_t i = 0u; i < capacity; ++i) {
    auto const row = static_cast<int>(i);
    source.id(row) = row;
    source.position_reference(row) =
        Utils::Vector3d{double(i), double(i) + 0.5, double(i) + 0.25};
    source.image_box_reference(row) =
        Utils::Vector3i{int(i), -int(i), 2 * int(i)};
  }

  // Target store: copy each source row into a consecutive target row.
  ParticleStore store{};
  store.begin_rebuild(capacity, 0u);
  int next_row = 0;
  for (std::size_t i = 0u; i < capacity; ++i) {
    store.copy_row(source, static_cast<int>(i), next_row++);
  }
  store.finish_rebuild();

  for (std::size_t i = 0u; i < capacity; ++i) {
    auto const row = static_cast<int>(i);
    BOOST_CHECK_EQUAL(store.id(row), row);
    auto const pos = store.position_value(row);
    BOOST_CHECK_EQUAL(pos[0], double(i));
    BOOST_CHECK_EQUAL(pos[1], double(i) + 0.5);
    auto const img = store.image_box_value(row);
    BOOST_CHECK_EQUAL(img[0], int(i));
    BOOST_CHECK_EQUAL(img[2], 2 * int(i));
    // Reading through a view over the row matches the column.
    BOOST_CHECK_EQUAL(store.make_view(row).pos()[0], double(i));
  }
}

// Parameter columns/sidecars: a rank-local rebuild that shuffles the row order
// must preserve each particle's parameters (representative subset: id, type,
// mass, gamma, ext_force, dip_fld) by old row, and seed a genuinely new row
// to the ParticleProperties defaults, exactly as the state columns do.
BOOST_AUTO_TEST_CASE(rebuild_preserves_parameter_columns_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.id(p0.store_row()) = 10;
  store.id(p1.store_row()) = 20;
  store.type(p0.store_row()) = 3;
  store.type(p1.store_row()) = 7;
#ifdef ESPRESSO_MASS
  store.mass(p0.store_row()) = 2.5;
  store.mass(p1.store_row()) = 4.25;
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  store.gamma_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.gamma_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};
#else
  store.gamma_reference(p0.store_row()) = 1.5;
  store.gamma_reference(p1.store_row()) = 2.5;
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_force_reference(p0.store_row()) = Utils::Vector3d{-1., -2., -3.};
  store.ext_force_reference(p1.store_row()) = Utils::Vector3d{7., 8., 9.};
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  store.dip_fld_reference(p0.store_row()) = Utils::Vector3d{0.1, 0.2, 0.3};
  store.dip_fld_reference(p1.store_row()) = Utils::Vector3d{0.4, 0.5, 0.6};
#endif

  // resort that swaps the two particles' order and appends a new one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.id(p1.store_row()), 20);
  BOOST_CHECK_EQUAL(store.id(p0.store_row()), 10);
  BOOST_CHECK_EQUAL(store.id(p2.store_row()), -1); // default
  BOOST_CHECK_EQUAL(store.type(p1.store_row()), 7);
  BOOST_CHECK_EQUAL(store.type(p0.store_row()), 3);
  BOOST_CHECK_EQUAL(store.type(p2.store_row()), 0); // default
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(store.mass(p1.store_row()), 4.25);
  BOOST_CHECK_EQUAL(store.mass(p0.store_row()), 2.5);
  BOOST_CHECK_EQUAL(store.mass(p2.store_row()), 1.0); // default
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d const g_p1 = store.gamma_value(p1.store_row());
  Utils::Vector3d const g_p0 = store.gamma_value(p0.store_row());
  Utils::Vector3d const g_p2 = store.gamma_value(p2.store_row());
  BOOST_CHECK_EQUAL(g_p1[0], 4.);
  BOOST_CHECK_EQUAL(g_p0[2], 3.);
  BOOST_CHECK_EQUAL(g_p2[0], -1.); // default {-1,-1,-1}
  BOOST_CHECK_EQUAL(g_p2[2], -1.);
#else
  BOOST_CHECK_EQUAL(store.gamma_value(p1.store_row()), 2.5);
  BOOST_CHECK_EQUAL(store.gamma_value(p0.store_row()), 1.5);
  BOOST_CHECK_EQUAL(store.gamma_value(p2.store_row()), -1.); // default
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d const ef_p1 = store.ext_force_value(p1.store_row());
  Utils::Vector3d const ef_p0 = store.ext_force_value(p0.store_row());
  Utils::Vector3d const ef_p2 = store.ext_force_value(p2.store_row());
  BOOST_CHECK_EQUAL(ef_p1[0], 7.);
  BOOST_CHECK_EQUAL(ef_p0[2], -3.);
  BOOST_CHECK_EQUAL(ef_p2.norm2(), 0.); // default {0,0,0}
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d const df_p1 = store.dip_fld_value(p1.store_row());
  Utils::Vector3d const df_p0 = store.dip_fld_value(p0.store_row());
  Utils::Vector3d const df_p2 = store.dip_fld_value(p2.store_row());
  BOOST_CHECK_EQUAL(df_p1[0], 0.4);
  BOOST_CHECK_EQUAL(df_p0[2], 0.3);
  BOOST_CHECK_EQUAL(df_p2.norm2(), 0.); // default {0,0,0}
#endif
}

// uint8 column (rotation / ext_flag): bitfield values must round-trip through
// the DualView<uint8_t*> column, write through the element reference, and
// preserve by old row across a rebuild.
#if defined(ESPRESSO_ROTATION) || defined(ESPRESSO_EXTERNAL_FORCES)
BOOST_AUTO_TEST_CASE(uint8_parameter_column_write_through_and_preserve) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

#ifdef ESPRESSO_ROTATION
  store.rotation(p0.store_row()) = static_cast<std::uint8_t>(0b101u);
  store.rotation(p1.store_row()) = static_cast<std::uint8_t>(0b010u);
  BOOST_CHECK_EQUAL(store.rotation(p0.store_row()),
                    static_cast<std::uint8_t>(0b101u));
  std::uint8_t &rot_ref = store.rotation(p1.store_row());
  rot_ref = static_cast<std::uint8_t>(0b111u);
  BOOST_CHECK_EQUAL(store.rotation(p1.store_row()),
                    static_cast<std::uint8_t>(0b111u));
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_flag(p0.store_row()) = static_cast<std::uint8_t>(0b011u);
  store.ext_flag(p1.store_row()) = static_cast<std::uint8_t>(0b100u);
#endif

  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(store.rotation(p1.store_row()),
                    static_cast<std::uint8_t>(0b111u));
  BOOST_CHECK_EQUAL(store.rotation(p0.store_row()),
                    static_cast<std::uint8_t>(0b101u));
  BOOST_CHECK_EQUAL(store.rotation(p2.store_row()),
                    static_cast<std::uint8_t>(0u)); // default
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK_EQUAL(store.ext_flag(p1.store_row()),
                    static_cast<std::uint8_t>(0b100u));
  BOOST_CHECK_EQUAL(store.ext_flag(p0.store_row()),
                    static_cast<std::uint8_t>(0b011u));
  BOOST_CHECK_EQUAL(store.ext_flag(p2.store_row()),
                    static_cast<std::uint8_t>(0u)); // default
#endif
}
#endif

// Host sidecar (POD std::vector): a POD written into a sidecar row must
// preserve by old row across a rebuild and a genuinely-new row must default to
// a value-constructed POD.
#ifdef ESPRESSO_ENGINE
BOOST_AUTO_TEST_CASE(swimming_sidecar_preserve_and_default) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.swimming(p0.store_row()).f_swim = 1.5;
  store.swimming(p0.store_row()).swimming = true;
  store.swimming(p1.store_row()).f_swim = -2.5;

  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.swimming(p1.store_row()).f_swim, -2.5);
  BOOST_CHECK_EQUAL(store.swimming(p0.store_row()).f_swim, 1.5);
  BOOST_CHECK(store.swimming(p0.store_row()).swimming);
  // genuinely-new row: value-constructed POD default (matches the
  // ParticleParametersSwimming member defaults)
  BOOST_CHECK_EQUAL(store.swimming(p2.store_row()).f_swim, 0.);
  BOOST_CHECK(not store.swimming(p2.store_row()).swimming);
}
#endif

// RATTLE observable column: a Vector3d written into the column must preserve
// by old row across a rank-local rebuild that shuffles the order, and a
// genuinely-new row must default to zero.
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_AUTO_TEST_CASE(rattle_correction_column_preserve_and_default) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.rattle_correction_reference(p0.store_row()) =
      Utils::Vector3d{1.0, 2.0, 3.0};
  store.rattle_correction_reference(p1.store_row()) =
      Utils::Vector3d{-4.0, -5.0, -6.0};

  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  BOOST_CHECK(store.rattle_correction_value(p1.store_row()) ==
              (Utils::Vector3d{-4.0, -5.0, -6.0}));
  BOOST_CHECK(store.rattle_correction_value(p0.store_row()) ==
              (Utils::Vector3d{1.0, 2.0, 3.0}));
  // genuinely-new row: preserve-or-default seeds a zero vector.
  BOOST_CHECK_EQUAL(store.rattle_correction_value(p2.store_row()).norm2(), 0.);
}
#endif

// Ragged bond sidecar: a non-empty BondList written into a sidecar row must
// survive a rank-local rebuild that shuffles the row order, with its contents
// intact. The preserve path MOVES the element out of the old vector (the
// BondList owns heap storage), so this also exercises the move path explicitly:
// the moved element's logical value (bonds + partner ids) is unchanged in the
// new generation. A genuinely-new row defaults to empty.
BOOST_AUTO_TEST_CASE(bonds_sidecar_preserve_moves_intact_and_default) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  // p0 gets two bonds of different arity (heap-backed ragged run): a pair bond
  // (id 7, partners {11}) and an angle bond (id 3, partners {12, 13}).
  {
    auto &bonds0 = store.bonds_sidecar_reference(p0.store_row());
    std::array<int, 1> const pair_partners{11};
    std::array<int, 2> const angle_partners{12, 13};
    bonds0.insert(BondView{7, pair_partners});
    bonds0.insert(BondView{3, angle_partners});
  }
  // p1 gets a single pair bond (id 5, partners {21}).
  {
    auto &bonds1 = store.bonds_sidecar_reference(p1.store_row());
    std::array<int, 1> const partners{21};
    bonds1.insert(BondView{5, partners});
  }

  // Snapshot p0's logical value (id + partner-ids of every bond) before the
  // rebuild moves it, to assert the move left it unchanged.
  auto flatten = [](BondList const &bonds) {
    std::vector<std::vector<int>> out;
    for (auto const bond : bonds) {
      std::vector<int> entry{bond.bond_id()};
      for (auto const pid : bond.partner_ids()) {
        entry.push_back(pid);
      }
      out.push_back(entry);
    }
    return out;
  };
  auto const p0_before = flatten(store.bonds_sidecar_reference(p0.store_row()));
  auto const p1_before = flatten(store.bonds_sidecar_reference(p1.store_row()));
  BOOST_REQUIRE_EQUAL(p0_before.size(), 2u);
  BOOST_REQUIRE_EQUAL(p1_before.size(), 1u);

  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  // The moved elements carry their exact logical value into the new rows.
  auto const p0_after = flatten(store.bonds_sidecar_reference(p0.store_row()));
  auto const p1_after = flatten(store.bonds_sidecar_reference(p1.store_row()));
  BOOST_CHECK(p0_after == p0_before);
  BOOST_CHECK(p1_after == p1_before);
  // Spell out p0's surviving contents to guard the flatten helper itself.
  BOOST_REQUIRE_EQUAL(p0_after.size(), 2u);
  BOOST_CHECK((p0_after[0] == std::vector<int>{7, 11}));
  BOOST_CHECK((p0_after[1] == std::vector<int>{3, 12, 13}));
  // genuinely-new row: empty bond list.
  BOOST_CHECK_EQUAL(store.bonds_sidecar_reference(p2.store_row()).size(), 0u);
}

#ifdef ESPRESSO_EXCLUSIONS
// Ragged exclusion sidecar: an exclusion id list written into a sidecar row
// must preserve (move) by old row across a shuffling rebuild; a genuinely-new
// / ghost row defaults to empty (exclusions are never ghost-transferred).
BOOST_AUTO_TEST_CASE(exclusions_sidecar_preserve_moves_intact_and_default) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  {
    auto &excl0 = store.exclusions_sidecar_reference(p0.store_row());
    excl0.push_back(101);
    excl0.push_back(102);
    excl0.push_back(103);
  }
  {
    auto &excl1 = store.exclusions_sidecar_reference(p1.store_row());
    excl1.push_back(201);
  }
  Utils::compact_vector<int> const p0_before =
      store.exclusions_sidecar_reference(p0.store_row());
  Utils::compact_vector<int> const p1_before =
      store.exclusions_sidecar_reference(p1.store_row());

  Particle p2{}; // a fresh (ghost-like) row: empty exclusion sidecar
  store.mark_dirty();
  store.begin_rebuild(2u, 1u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  BOOST_CHECK(store.exclusions_sidecar_reference(p0.store_row()) == p0_before);
  BOOST_CHECK(store.exclusions_sidecar_reference(p1.store_row()) == p1_before);
  BOOST_CHECK_EQUAL(store.exclusions_sidecar_reference(p2.store_row()).size(),
                    0u);
}
#endif // ESPRESSO_EXCLUSIONS

// A row that migrates in carries its ragged sidecars (bond list, exclusion
// list) via the row-to-row ParticleStore::copy_row transfer. Set non-empty
// sidecars on a source row, copy it into a target store row, and verify the
// sidecar contents survive.
BOOST_AUTO_TEST_CASE(rebuild_seeds_ragged_sidecars_from_carrier) {
  ParticleStore source{};
  source.begin_rebuild(1u, 0u);
  source.seed_default_row(0);
  source.finish_rebuild();
  {
    std::array<int, 1> const partners{99};
    source.bonds_sidecar_reference(0).insert(BondView{2, partners});
  }
#ifdef ESPRESSO_EXCLUSIONS
  source.exclusions_sidecar_reference(0).push_back(55);
#endif

  ParticleStore store{};
  store.begin_rebuild(1u, 0u);
  store.copy_row(source, 0, 0);
  store.finish_rebuild();

  auto const &seeded_bonds = store.bonds_sidecar_reference(0);
  BOOST_REQUIRE_EQUAL(seeded_bonds.size(), 1u);
  auto const bond = *seeded_bonds.begin();
  BOOST_CHECK_EQUAL(bond.bond_id(), 2);
  BOOST_REQUIRE_EQUAL(bond.partner_ids().size(), 1u);
  BOOST_CHECK_EQUAL(bond.partner_ids()[0], 99);
#ifdef ESPRESSO_EXCLUSIONS
  auto const &seeded_excl = store.exclusions_sidecar_reference(0);
  BOOST_REQUIRE_EQUAL(seeded_excl.size(), 1u);
  BOOST_CHECK_EQUAL(seeded_excl[0], 55);
#endif
}

// A genuinely new row's parameter columns/sidecars must default to
// ParticleProperties' member defaults, seeded from a freshly-constructed
// (detached) particle.
BOOST_AUTO_TEST_CASE(new_row_parameter_defaults) {
  ParticleStore store{};
  Particle p{};
  BOOST_REQUIRE(p.store() == nullptr);
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.id(p.store_row()), -1);
  BOOST_CHECK_EQUAL(store.mol_id(p.store_row()), 0);
  BOOST_CHECK_EQUAL(store.type(p.store_row()), 0);
  BOOST_CHECK_EQUAL(store.propagation(p.store_row()),
                    static_cast<int>(PropagationMode::SYSTEM_DEFAULT));
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(store.mass(p.store_row()), 1.0);
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  BOOST_CHECK_EQUAL(store.q(p.store_row()), 0.0);
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d const ri = store.rinertia_value(p.store_row());
  BOOST_CHECK_EQUAL(ri[0], 1.);
  BOOST_CHECK_EQUAL(ri[1], 1.);
  BOOST_CHECK_EQUAL(ri[2], 1.);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK_EQUAL(store.ext_force_value(p.store_row()).norm2(), 0.);
  BOOST_CHECK_EQUAL(store.ext_flag(p.store_row()),
                    static_cast<std::uint8_t>(0u));
#endif
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(store.rotation(p.store_row()),
                    static_cast<std::uint8_t>(0u));
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d const g = store.gamma_value(p.store_row());
  BOOST_CHECK_EQUAL(g[0], -1.);
  BOOST_CHECK_EQUAL(g[2], -1.);
#else
  BOOST_CHECK_EQUAL(store.gamma_value(p.store_row()), -1.);
#endif
#endif
}

// A row that migrates in carries its parameter columns and POD sidecars via
// the row-to-row ParticleStore::copy_row transfer. Set known parameters on a
// source row, copy it into a target store row, and verify they survive.
BOOST_AUTO_TEST_CASE(rebuild_seeds_parameters_from_carrier) {
  ParticleStore source{};
  source.begin_rebuild(1u, 0u);
  source.seed_default_row(0);
  source.finish_rebuild();
  source.id(0) = 42;
  source.type(0) = 5;
#ifdef ESPRESSO_MASS
  source.mass(0) = 3.5;
#endif
#ifdef ESPRESSO_ENGINE
  source.swimming(0).f_swim = 9.5;
#endif

  ParticleStore store{};
  store.begin_rebuild(1u, 0u);
  store.copy_row(source, 0, 0);
  store.finish_rebuild();

  BOOST_CHECK_EQUAL(store.id(0), 42);
  BOOST_CHECK_EQUAL(store.type(0), 5);
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(store.mass(0), 3.5);
#endif
#ifdef ESPRESSO_ENGINE
  BOOST_CHECK_EQUAL(store.swimming(0).f_swim, 9.5);
#endif
}

// make_view factory: a view built for a row reads that row's columns (identity:
// store.make_view(r).id() == store.id(r)) and writes through it (a write via
// the view lands in the column, readable via the store by row).
BOOST_AUTO_TEST_CASE(make_view_identity_and_write_through) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();

  store.id(0) = 10;
  store.id(1) = 20;
  store.position_reference(0) = Utils::Vector3d{1., 2., 3.};
  store.position_reference(1) = Utils::Vector3d{4., 5., 6.};

  // Identity: the view over row r reports the store's id at row r, and the view
  // is attached to this store at that row.
  for (int r = 0; r < 2; ++r) {
    auto view = store.make_view(r);
    BOOST_REQUIRE(view.store() == &store);
    BOOST_CHECK_EQUAL(view.store_row(), r);
    BOOST_CHECK_EQUAL(view.id(), store.id(r));
    BOOST_CHECK_EQUAL(view.pos()[0], store.position_value(r)[0]);
    BOOST_CHECK_EQUAL(view.pos()[2], store.position_value(r)[2]);
  }

  // Write-through: mutating a column through the view is visible via the store
  // (the view aliases the row, it does not snapshot it).
  auto view = store.make_view(0);
  view.pos() = Utils::Vector3d{-7., -8., -9.};
  Utils::Vector3d const back = store.position_value(0);
  BOOST_CHECK_EQUAL(back[0], -7.);
  BOOST_CHECK_EQUAL(back[1], -8.);
  BOOST_CHECK_EQUAL(back[2], -9.);
  // And the OTHER row is untouched (no cross-row aliasing).
  BOOST_CHECK_EQUAL(store.position_value(1)[0], 4.);
}

// Scalar column references write through to the stored value.
BOOST_AUTO_TEST_CASE(scalar_column_references_write_through) {
  ParticleStore store{};
  Particle p{};
  store.begin_rebuild(1u, 0u);
  store.assign_row(p, 0);
  store.finish_rebuild();

  store.lees_edwards_offset(p.store_row()) = 3.14;
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p.store_row()), 3.14);
  double &offset_ref = store.lees_edwards_offset(p.store_row());
  offset_ref = -1.0;
  BOOST_CHECK_EQUAL(store.lees_edwards_offset(p.store_row()), -1.0);

  store.lees_edwards_flag(p.store_row()) = static_cast<short>(5);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p.store_row()),
                    static_cast<short>(5));
  short &flag_ref = store.lees_edwards_flag(p.store_row());
  flag_ref = static_cast<short>(-7);
  BOOST_CHECK_EQUAL(store.lees_edwards_flag(p.store_row()),
                    static_cast<short>(-7));
}

// -- permute_rebuild ----------------------------------------------------------
// The permutation rebuild is the resort-as-permutation kernel: new-row data is
// moved from old_row = permutation[new_row] (a survivor) or seeded to defaults
// (a negative entry: staged / fresh ghost). These tests use the maximal-
// population helpers so EVERY ifdef field is exercised field-for-field; a field
// missing from permute_rebuild's four-way sync list fails here.

using maximal_population::check_row_equal;
using maximal_population::fill_maximal;
using maximal_population::make_store;

namespace {
// Build a store of `count` maximally-populated rows (each row seeded from a
// distinct sentinel so a mis-permuted row is caught). The data lives in the
// current generation; permute_rebuild swaps it to the retired generation and
// permutes from there.
ParticleStore make_maximal_store(std::size_t const count) {
  auto store = make_store(count);
  for (std::size_t r = 0u; r < count; ++r) {
    fill_maximal(store, static_cast<int>(r), 1000. + 7. * double(r));
  }
  return store;
}
} // namespace

// Identity permutation: every row survives in place. The permuted store must be
// field-for-field identical to the reference (the pre-rebuild data). This is
// the removal-free bitwise-identity contract on the store side.
BOOST_AUTO_TEST_CASE(permute_rebuild_identity) {
  constexpr std::size_t n = 5u;
  auto const reference = make_maximal_store(n);
  auto store = make_maximal_store(n);

  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0); // 0,1,2,3,4
  store.permute_rebuild(std::span<int const>{perm}, n, 0u);

  BOOST_CHECK_EQUAL(store.number_of_local_particles(), n);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 0u);
  for (std::size_t r = 0u; r < n; ++r) {
    check_row_equal(reference, static_cast<int>(r), store, static_cast<int>(r));
  }
}

// Reversal permutation: new row r takes old row n-1-r. Verifies the permute
// kernels read the correct source row (not the identity) for every column and
// ragged sidecar.
BOOST_AUTO_TEST_CASE(permute_rebuild_reversal) {
  constexpr std::size_t n = 6u;
  auto const reference = make_maximal_store(n);
  auto store = make_maximal_store(n);

  std::vector<int> perm(n);
  for (std::size_t r = 0u; r < n; ++r) {
    perm[r] = static_cast<int>(n - 1u - r);
  }
  store.permute_rebuild(std::span<int const>{perm}, n, 0u);

  for (std::size_t r = 0u; r < n; ++r) {
    check_row_equal(reference, perm[r], store, static_cast<int>(r));
  }
}

// Random (arbitrary bijection) permutation: each new row must hold exactly the
// reference row named by the permutation, incl. ragged contents (bonds moved
// intact under the move-by-permutation).
BOOST_AUTO_TEST_CASE(permute_rebuild_random) {
  constexpr std::size_t n = 8u;
  auto const reference = make_maximal_store(n);
  auto store = make_maximal_store(n);

  // A fixed arbitrary bijection of [0, 8).
  std::vector<int> const perm{3, 0, 7, 1, 6, 2, 5, 4};
  BOOST_REQUIRE_EQUAL(perm.size(), n);
  store.permute_rebuild(std::span<int const>{perm}, n, 0u);

  for (std::size_t r = 0u; r < n; ++r) {
    check_row_equal(reference, perm[r], store, static_cast<int>(r));
  }
}

// Ragged contents survive a permutation intact: after a reversal, each new
// row's bond list matches the reference row's bond list exactly (already
// covered by check_row_equal, but spelled out for the ragged run explicitly).
BOOST_AUTO_TEST_CASE(permute_rebuild_ragged_intact) {
  constexpr std::size_t n = 4u;
  auto const reference = make_maximal_store(n);
  auto store = make_maximal_store(n);

  std::vector<int> const perm{2, 3, 0, 1};
  store.permute_rebuild(std::span<int const>{perm}, n, 0u);

  for (std::size_t r = 0u; r < n; ++r) {
    using maximal_population::flatten_bonds;
    BOOST_CHECK(
        flatten_bonds(reference.bonds_sidecar_reference(perm[r])) ==
        flatten_bonds(store.bonds_sidecar_reference(static_cast<int>(r))));
#ifdef ESPRESSO_EXCLUSIONS
    BOOST_CHECK(reference.exclusions_sidecar_reference(perm[r]) ==
                store.exclusions_sidecar_reference(static_cast<int>(r)));
#endif
  }
}

// Ghost tail freshly seeded to defaults: negative permutation entries on the
// ghost suffix must produce the new-particle defaults (id -1, position zero,
// quaternion identity, mass 1, empty bonds), NOT stale old-generation data.
// Locals survive by permutation; ghosts are seeded.
BOOST_AUTO_TEST_CASE(permute_rebuild_ghost_tail_seeded) {
  constexpr std::size_t n_local = 3u;
  constexpr std::size_t n_ghost = 2u;
  auto const reference = make_maximal_store(n_local);
  // Working store holds n_local + n_ghost maximal rows so the ghost tail has
  // NON-default data in the retired generation; permute_rebuild must overwrite
  // it with defaults, not preserve it.
  auto store = make_maximal_store(n_local + n_ghost);

  // Locals survive in place; ghost tail entries are -1 (seed defaults).
  std::vector<int> perm(n_local + n_ghost, -1);
  for (std::size_t r = 0u; r < n_local; ++r) {
    perm[r] = static_cast<int>(r);
  }
  store.permute_rebuild(std::span<int const>{perm}, n_local, n_ghost);

  BOOST_CHECK_EQUAL(store.number_of_local_particles(), n_local);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), n_ghost);

  // Locals: field-for-field identical to the reference.
  for (std::size_t r = 0u; r < n_local; ++r) {
    check_row_equal(reference, static_cast<int>(r), store, static_cast<int>(r));
  }
  // Ghost tail: the new-particle defaults (matches seed_default_row).
  auto const defaults = make_store(1u); // one all-defaults row
  for (std::size_t r = n_local; r < n_local + n_ghost; ++r) {
    check_row_equal(defaults, 0, store, static_cast<int>(r));
  }
}

// A staged local (negative entry in the LOCAL range) is seeded to defaults by
// permute_rebuild; the caller overwrites it via copy_row afterwards. Verify
// the seed-then-copy composition reproduces the source row.
BOOST_AUTO_TEST_CASE(permute_rebuild_staged_local_seed_then_copy) {
  constexpr std::size_t n = 4u;
  auto const reference = make_maximal_store(n);
  auto store = make_maximal_store(n);

  // A separate source store holds the staged particle's data.
  auto source = make_store(1u);
  fill_maximal(source, 0, 9999.);

  // New layout: rows 0,1 survive (old 0,1); row 2 is a staged local (-1);
  // row 3 survives (old 2). n_local = 4, no ghosts.
  std::vector<int> const perm{0, 1, -1, 2};
  store.permute_rebuild(std::span<int const>{perm}, n, 0u);
  // Caller commits the staged local: copy the source row into the seeded slot.
  store.copy_row(source, 0, 2);

  check_row_equal(reference, 0, store, 0);
  check_row_equal(reference, 1, store, 1);
  check_row_equal(source, 0, store, 2);
  check_row_equal(reference, 2, store, 3);
}

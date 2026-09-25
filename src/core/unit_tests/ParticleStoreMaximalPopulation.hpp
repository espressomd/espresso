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

#pragma once

// Maximal-population store helpers shared by the field-completeness tests.
// fill_maximal() writes EVERY ifdef-guarded store field to a distinct sentinel
// derived from a seed; check_row_equal() compares two rows field-for-field
// (incl. ragged contents). A field missed by any row-moving path (MigrationPack
// pack/unpack, copy_row, permute_rebuild) FAILS a test rather than being
// silently dropped. Shared between MigrationPack_test and ParticleStore_test.

#include <config/config.hpp>

#include "particle_store/ParticleStore.hpp"

#include "BondList.hpp"
#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/quaternion.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace maximal_population {

// Build an empty store of exactly `count` rows (all defaults) and return it.
inline ParticleStore make_store(std::size_t const count) {
  ParticleStore store{};
  store.begin_rebuild(count, 0u);
  // assign_row every row so the columns/sidecars are seeded to their defaults
  // (WithoutInitializing storage is otherwise garbage).
  std::vector<Particle> ps(count);
  for (std::size_t r = 0u; r < count; ++r) {
    store.assign_row(ps[r], static_cast<int>(r));
  }
  store.finish_rebuild();
  return store;
}

// Distinct-sentinel fill of EVERY ifdef field of one store row. The `seed`
// disambiguates rows (each field's value is derived from it) so a bug that
// copies the wrong row is caught too. `with_bonds`/`with_exclusions` toggle the
// ragged legs so the edge cases (bond-free particles) can reuse this helper.
inline void fill_maximal(ParticleStore &store, int const row, double const seed,
                         bool const with_bonds = true,
                         bool const with_exclusions = true) {
  auto const s = seed;
  auto v = [s](double a, double b, double c) {
    return Utils::Vector3d{s + a, s + b, s + c};
  };

  // POSITION leg.
  store.position_reference(row) = v(1., 2., 3.);
  store.image_box_reference(row) =
      Utils::Vector3i{static_cast<int>(s) + 1, static_cast<int>(s) - 2,
                      static_cast<int>(s) + 3};
#ifdef ESPRESSO_ROTATION
  store.quaternion_reference(row) =
      Utils::Quaternion<double>{{s + 0.1, s + 0.2, s + 0.3, s + 0.4}};
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  store.position_last_time_step_reference(row) = v(4., 5., 6.);
#endif
  store.position_at_last_verlet_update_reference(row) = v(7., 8., 9.);
  store.lees_edwards_offset(row) = s + 0.5;
  store.lees_edwards_flag(row) = static_cast<short>(static_cast<int>(s) + 11);

  // MOMENTUM leg.
  store.velocity_reference(row) = v(10., 11., 12.);
#ifdef ESPRESSO_ROTATION
  store.angular_velocity_reference(row) = v(13., 14., 15.);
#endif

  // PROPRTS leg.
  store.id(row) = static_cast<int>(s) + 100;
  store.mol_id(row) = static_cast<int>(s) + 200;
  store.type(row) = static_cast<int>(s) + 5;
  store.propagation(row) = static_cast<int>(s) + 6;
#ifdef ESPRESSO_ROTATION
  store.rotation(row) = static_cast<std::uint8_t>(static_cast<int>(s) % 7 + 1);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_flag(row) = static_cast<std::uint8_t>(static_cast<int>(s) % 5 + 1);
#endif
#ifdef ESPRESSO_MASS
  store.mass(row) = s + 20.;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  store.q(row) = s - 30.;
#endif
#ifdef ESPRESSO_DIPOLES
  store.dipm(row) = s + 40.;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  store.rinertia_reference(row) = v(16., 17., 18.);
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  store.mu_E_reference(row) = v(19., 20., 21.);
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  store.dip_fld_reference(row) = v(22., 23., 24.);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_force_reference(row) = v(25., 26., 27.);
#ifdef ESPRESSO_ROTATION
  store.ext_torque_reference(row) = v(28., 29., 30.);
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  store.gamma_reference(row) = v(31., 32., 33.);
#ifdef ESPRESSO_ROTATION
  store.gamma_rot_reference(row) = v(34., 35., 36.);
#endif
#else
  store.gamma_reference(row) = s + 50.;
#ifdef ESPRESSO_ROTATION
  store.gamma_rot_reference(row) = s + 60.;
#endif
#endif
#endif
#ifdef ESPRESSO_ENGINE
  store.swimming(row).f_swim = s + 70.;
  store.swimming(row).swimming = true;
  store.swimming(row).is_engine_force_on_fluid = true;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  store.magnetodynamics(row).is_enabled = true;
  store.magnetodynamics(row).phi0 = s + 80.;
  store.magnetodynamics(row).sat_mag = s + 81.;
  store.magnetodynamics(row).ani_fld_inv = s + 82.;
  store.magnetodynamics(row).ani_energy = s + 83.;
  store.magnetodynamics(row).tau0_inv = s + 84.;
  store.magnetodynamics(row).dt_incr = s + 85.;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  store.vs_relative(row).to_particle_id = static_cast<int>(s) + 300;
  store.vs_relative(row).distance = s + 90.;
  store.vs_relative(row).rel_orientation =
      Utils::Quaternion<double>{{s + 0.5, s + 0.6, s + 0.7, s + 0.8}};
  store.vs_relative(row).quat =
      Utils::Quaternion<double>{{s + 0.9, s + 1.0, s + 1.1, s + 1.2}};
#endif

  // FORCE leg.
  store.force_reference(row) = v(37., 38., 39.);
#ifdef ESPRESSO_ROTATION
  store.torque_reference(row) = v(40., 41., 42.);
#endif

  // Ragged legs.
  auto &bonds = store.bonds_sidecar_reference(row);
  bonds.clear();
  if (with_bonds) {
    std::array<int, 1> const pair_partners{static_cast<int>(s) + 401};
    std::array<int, 2> const angle_partners{static_cast<int>(s) + 402,
                                            static_cast<int>(s) + 403};
    bonds.insert(BondView{static_cast<int>(s) % 9 + 1, pair_partners});
    bonds.insert(BondView{static_cast<int>(s) % 9 + 2, angle_partners});
  }
#ifdef ESPRESSO_EXCLUSIONS
  auto &excl = store.exclusions_sidecar_reference(row);
  excl.clear();
  if (with_exclusions) {
    excl.push_back(static_cast<int>(s) + 501);
    excl.push_back(static_cast<int>(s) + 502);
    excl.push_back(static_cast<int>(s) + 503);
  }
#endif
}

// Flatten a bond list to (id, partner...) tuples for order-preserving compare.
inline std::vector<std::vector<int>> flatten_bonds(BondList const &bonds) {
  std::vector<std::vector<int>> out;
  for (auto const bond : bonds) {
    std::vector<int> entry{bond.bond_id()};
    for (auto const pid : bond.partner_ids()) {
      entry.push_back(pid);
    }
    out.push_back(entry);
  }
  return out;
}

// Field-for-field equality of two store rows (incl. ragged contents). This is
// the completeness check: every field fill_maximal() writes is verified here.
inline void check_row_equal(ParticleStore const &a, int const ra,
                            ParticleStore const &b, int const rb) {
  BOOST_CHECK(a.position_value(ra) == b.position_value(rb));
  BOOST_CHECK(a.image_box_value(ra) == b.image_box_value(rb));
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(a.quaternion_value(ra) == b.quaternion_value(rb));
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  BOOST_CHECK(a.position_last_time_step_value(ra) ==
              b.position_last_time_step_value(rb));
#endif
  BOOST_CHECK(a.position_at_last_verlet_update_value(ra) ==
              b.position_at_last_verlet_update_value(rb));
  BOOST_CHECK_EQUAL(a.lees_edwards_offset(ra), b.lees_edwards_offset(rb));
  BOOST_CHECK_EQUAL(a.lees_edwards_flag(ra), b.lees_edwards_flag(rb));

  BOOST_CHECK(a.velocity_value(ra) == b.velocity_value(rb));
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(a.angular_velocity_value(ra) == b.angular_velocity_value(rb));
#endif

  BOOST_CHECK_EQUAL(a.id(ra), b.id(rb));
  BOOST_CHECK_EQUAL(a.mol_id(ra), b.mol_id(rb));
  BOOST_CHECK_EQUAL(a.type(ra), b.type(rb));
  BOOST_CHECK_EQUAL(a.propagation(ra), b.propagation(rb));
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK_EQUAL(a.rotation(ra), b.rotation(rb));
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK_EQUAL(a.ext_flag(ra), b.ext_flag(rb));
#endif
#ifdef ESPRESSO_MASS
  BOOST_CHECK_EQUAL(a.mass(ra), b.mass(rb));
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  BOOST_CHECK_EQUAL(a.q(ra), b.q(rb));
#endif
#ifdef ESPRESSO_DIPOLES
  BOOST_CHECK_EQUAL(a.dipm(ra), b.dipm(rb));
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  BOOST_CHECK(a.rinertia_value(ra) == b.rinertia_value(rb));
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  BOOST_CHECK(a.mu_E_value(ra) == b.mu_E_value(rb));
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  BOOST_CHECK(a.dip_fld_value(ra) == b.dip_fld_value(rb));
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  BOOST_CHECK(a.ext_force_value(ra) == b.ext_force_value(rb));
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(a.ext_torque_value(ra) == b.ext_torque_value(rb));
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  BOOST_CHECK(a.gamma_value(ra) == b.gamma_value(rb));
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(a.gamma_rot_value(ra) == b.gamma_rot_value(rb));
#endif
#endif
#ifdef ESPRESSO_ENGINE
  BOOST_CHECK_EQUAL(a.swimming(ra).f_swim, b.swimming(rb).f_swim);
  BOOST_CHECK_EQUAL(a.swimming(ra).swimming, b.swimming(rb).swimming);
  BOOST_CHECK_EQUAL(a.swimming(ra).is_engine_force_on_fluid,
                    b.swimming(rb).is_engine_force_on_fluid);
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).is_enabled,
                    b.magnetodynamics(rb).is_enabled);
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).phi0, b.magnetodynamics(rb).phi0);
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).sat_mag,
                    b.magnetodynamics(rb).sat_mag);
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).ani_fld_inv,
                    b.magnetodynamics(rb).ani_fld_inv);
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).ani_energy,
                    b.magnetodynamics(rb).ani_energy);
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).tau0_inv,
                    b.magnetodynamics(rb).tau0_inv);
  BOOST_CHECK_EQUAL(a.magnetodynamics(ra).dt_incr,
                    b.magnetodynamics(rb).dt_incr);
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  BOOST_CHECK_EQUAL(a.vs_relative(ra).to_particle_id,
                    b.vs_relative(rb).to_particle_id);
  BOOST_CHECK_EQUAL(a.vs_relative(ra).distance, b.vs_relative(rb).distance);
  BOOST_CHECK(a.vs_relative(ra).rel_orientation ==
              b.vs_relative(rb).rel_orientation);
  BOOST_CHECK(a.vs_relative(ra).quat == b.vs_relative(rb).quat);
#endif

  BOOST_CHECK(a.force_value(ra) == b.force_value(rb));
#ifdef ESPRESSO_ROTATION
  BOOST_CHECK(a.torque_value(ra) == b.torque_value(rb));
#endif

  BOOST_CHECK(flatten_bonds(a.bonds_sidecar_reference(ra)) ==
              flatten_bonds(b.bonds_sidecar_reference(rb)));
#ifdef ESPRESSO_EXCLUSIONS
  BOOST_CHECK(a.exclusions_sidecar_reference(ra) ==
              b.exclusions_sidecar_reference(rb));
#endif
}

} // namespace maximal_population

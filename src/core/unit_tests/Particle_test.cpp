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

#define BOOST_TEST_MODULE "Particle struct test"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <config/config.hpp>

#include "Particle.hpp"
#include "ParticleStoreTestFixture.hpp"
#include "PropagationMode.hpp"

#include <utils/serialization/memcpy_archive.hpp>

#include <type_traits>
#include <vector>

void check_particle_force(ParticleForce const &out, ParticleForce const &ref) {
  BOOST_TEST(out.f == ref.f, boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  BOOST_TEST(out.torque == ref.torque, boost::test_tools::per_element());
#endif
}

BOOST_AUTO_TEST_CASE(comparison) {
  {
    ParticleStoreTestFixture fx{};
    auto p = fx.make();
    auto q = fx.make();

    p.id() = 1;
    q.id() = 2;

    BOOST_CHECK(p != q);
    BOOST_CHECK(not(p == q));
  }

  {
    ParticleStoreTestFixture fx{};
    auto p = fx.make();
    auto q = fx.make();

    p.id() = 2;
    q.id() = 2;

    BOOST_CHECK(not(p != q));
    BOOST_CHECK(p == q);
  }
}

// A Particle is a two-word non-owning view and cannot be boost-serialized.
// Per-field cross-rank transfer is covered by MigrationPack_test.cpp; the
// row-to-row copy is exercised by ParticleStore::copy_row in
// ParticleStore_test.cpp.

namespace Utils {
template <>
struct is_statically_serializable<ParticleProperties> : std::true_type {};
} // namespace Utils

BOOST_AUTO_TEST_CASE(properties_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ParticleProperties>();

  BOOST_CHECK_LE(expected_size, sizeof(ParticleProperties));

  std::vector<char> buf(expected_size);

  auto prop = ParticleProperties{};
  prop.identity = 1234;

  {
    auto oa = Utils::MemcpyOArchive{buf};

    oa << prop;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{buf};
    ParticleProperties out;

    ia >> out;
    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    BOOST_CHECK_EQUAL(out.identity, prop.identity);
  }
}

namespace Utils {
template <>
struct is_statically_serializable<ParticleForce> : std::true_type {};
} // namespace Utils

BOOST_AUTO_TEST_CASE(force_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ParticleForce>();

  BOOST_CHECK_LE(expected_size, sizeof(ParticleForce));

  std::vector<char> buf(expected_size);

  auto pf = ParticleForce{{1, 2, 3}};
#ifdef ESPRESSO_ROTATION
  pf.torque = {4, 5, 6};
#endif

  {
    auto oa = Utils::MemcpyOArchive{buf};

    oa << pf;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{buf};
    ParticleForce out;

    ia >> out;

    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    check_particle_force(out, pf);
  }
}

BOOST_AUTO_TEST_CASE(force_constructors) {

  auto pf = ParticleForce{{1, 2, 3}};
#ifdef ESPRESSO_ROTATION
  pf.torque = {4, 5, 6};
#endif

  // check copy constructor
  {
    ParticleForce out(pf);
    check_particle_force(out, pf);
  }

  // check copy assignment operator
  {
    ParticleForce out; // avoid copy elision
    out = pf;
    check_particle_force(out, pf);
  }
}

// ParticleRattle is an empty type anchor; the correction Vector3d lives in a
// ParticleStore observable column. The correction round-trip is exercised by
// the store column tests (ParticleStore_test.cpp) and the RATTLE ghost path.

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH

void check_particle_tsw(ThermalStonerWohlfarthParameters const &out,
                        ThermalStonerWohlfarthParameters const &ref) {
  BOOST_TEST(out.is_enabled == ref.is_enabled);
  BOOST_TEST(out.phi0 == ref.phi0);
  BOOST_TEST(out.sat_mag == ref.sat_mag);
  BOOST_TEST(out.ani_fld_inv == ref.ani_fld_inv);
  BOOST_TEST(out.ani_energy == ref.ani_energy);
  BOOST_TEST(out.tau0_inv == ref.tau0_inv);
  BOOST_TEST(out.dt_incr == ref.dt_incr);
}

BOOST_AUTO_TEST_CASE(thermal_stoner_wohlfarth_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ThermalStonerWohlfarthParameters>();

  BOOST_CHECK_LE(expected_size, sizeof(ThermalStonerWohlfarthParameters));

  std::vector<char> buf(expected_size);

  auto pr = ThermalStonerWohlfarthParameters{.is_enabled = false,
                                             .phi0 = 1.,
                                             .sat_mag = 2.,
                                             .ani_fld_inv = 3.,
                                             .ani_energy = 4.,
                                             .tau0_inv = 5.,
                                             .dt_incr = 6.};

  {
    auto oa = Utils::MemcpyOArchive{buf};

    oa << pr;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{buf};
    ThermalStonerWohlfarthParameters out;

    ia >> out;

    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    check_particle_tsw(out, pr);
  }
}

BOOST_AUTO_TEST_CASE(thermal_stoner_wohlfarth_constructors) {
  auto pr = ThermalStonerWohlfarthParameters{.is_enabled = false,
                                             .phi0 = 1.,
                                             .sat_mag = 2.,
                                             .ani_fld_inv = 3.,
                                             .ani_energy = 4.,
                                             .tau0_inv = 5.,
                                             .dt_incr = 6.};

  // check copy constructor
  {
    ThermalStonerWohlfarthParameters out(pr);
    check_particle_tsw(out, pr);
  }

  // check copy assignment operator
  {
    ThermalStonerWohlfarthParameters out; // avoid copy elision
    out = pr;
    check_particle_tsw(out, pr);
  }
}
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH

BOOST_AUTO_TEST_CASE(particle_bitfields) {
  ParticleStoreTestFixture fx{};
  auto p = fx.make();

  // check default values
  BOOST_CHECK(not p.has_fixed_coordinates());
  BOOST_CHECK(not p.can_rotate());
  BOOST_CHECK(not p.is_fixed_along(1));
  BOOST_CHECK(not p.can_rotate_around(1));

  // check setting of one axis
#ifdef ESPRESSO_EXTERNAL_FORCES
  p.set_fixed_along(1, true);
  BOOST_CHECK(p.is_fixed_along(1));
  BOOST_CHECK(p.has_fixed_coordinates());
#endif
#ifdef ESPRESSO_ROTATION
  p.set_can_rotate_around(1, true);
  BOOST_CHECK(p.can_rotate_around(1));
  BOOST_CHECK(p.can_rotate());
#endif

  // check that unsetting is properly registered
#ifdef ESPRESSO_EXTERNAL_FORCES
  p.set_fixed_along(1, false);
  BOOST_CHECK(not p.has_fixed_coordinates());
#endif
#ifdef ESPRESSO_ROTATION
  p.set_can_rotate_around(1, false);
  BOOST_CHECK(not p.can_rotate());
#endif

  // check setting of all flags at once
#ifdef ESPRESSO_ROTATION
  p.set_can_rotate_all_axes();
  BOOST_CHECK(p.can_rotate_around(0));
  BOOST_CHECK(p.can_rotate_around(1));
  BOOST_CHECK(p.can_rotate_around(2));
  p.set_cannot_rotate_all_axes();
  BOOST_CHECK(not p.can_rotate());
#endif

  static_assert(
      std::is_same_v<std::underlying_type_t<PropagationMode::PropagationMode>,
                     decltype(ParticleProperties::propagation)>);
}

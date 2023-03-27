/*
 * Copyright (C) 2017-2022 The ESPResSo project
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

/* Unit tests for the Particle struct. */

#define BOOST_TEST_MODULE Particle test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "config/config.hpp"

#include <utils/Span.hpp>
#include <utils/compact_vector.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <algorithm>
#include <array>
#include <sstream>
#include <utility>
#include <vector>

void check_particle_force(ParticleForce const &out, ParticleForce const &ref) {
  BOOST_TEST(out.f == ref.f, boost::test_tools::per_element());
#ifdef ROTATION
  BOOST_TEST(out.torque == ref.torque, boost::test_tools::per_element());
#endif
}

BOOST_AUTO_TEST_CASE(comparison) {
  {
    Particle p, q;

    p.id() = 1;
    q.id() = 2;

    BOOST_CHECK(p != q);
    BOOST_CHECK(not(p == q));
  }

  {
    Particle p, q;

    p.id() = 2;
    q.id() = 2;

    BOOST_CHECK(not(p != q));
    BOOST_CHECK(p == q);
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  auto p = Particle();

  auto const bond_id = 5;
  auto const bond_partners = std::array<const int, 3>{{12, 13, 14}};

  p.id() = 15;
  p.bonds().insert({bond_id, bond_partners});
  p.force() = {1., -2., 3.};
#ifdef ROTATION
  p.torque() = {-4., 5., -6.};
#endif
#ifdef EXCLUSIONS
  std::vector<int> el = {5, 6, 7, 8};
  p.exclusions() = Utils::compact_vector<int>{el.begin(), el.end()};
#endif

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << p;

  boost::archive::text_iarchive in_ar(stream);
  auto q = Particle();
  in_ar >> q;

  auto const &pf = std::as_const(p).force_and_torque();
  BOOST_CHECK(q.id() == p.id());
  BOOST_CHECK((*q.bonds().begin() == BondView{bond_id, bond_partners}));
  BOOST_TEST(q.force() == pf.f, boost::test_tools::per_element());
#ifdef ROTATION
  BOOST_TEST(q.torque() == pf.torque, boost::test_tools::per_element());
#endif
  check_particle_force(q.force_and_torque(), pf);
}

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
    auto oa = Utils::MemcpyOArchive{Utils::make_span(buf)};

    oa << prop;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{Utils::make_span(buf)};
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
#ifdef ROTATION
  pf.torque = {4, 5, 6};
#endif

  {
    auto oa = Utils::MemcpyOArchive{Utils::make_span(buf)};

    oa << pf;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{Utils::make_span(buf)};
    ParticleForce out;

    ia >> out;

    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    check_particle_force(out, pf);
  }
}

BOOST_AUTO_TEST_CASE(force_constructors) {

  auto pf = ParticleForce{{1, 2, 3}};
#ifdef ROTATION
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

#ifdef BOND_CONSTRAINT

void check_particle_rattle(ParticleRattle const &out,
                           ParticleRattle const &ref) {
  BOOST_TEST(out.correction == ref.correction,
             boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(rattle_serialization) {
  auto const expected_size =
      Utils::MemcpyOArchive::packing_size<ParticleRattle>();

  BOOST_CHECK_LE(expected_size, sizeof(ParticleRattle));

  std::vector<char> buf(expected_size);

  auto pr = ParticleRattle{{1, 2, 3}};

  {
    auto oa = Utils::MemcpyOArchive{Utils::make_span(buf)};

    oa << pr;

    BOOST_CHECK_EQUAL(oa.bytes_written(), expected_size);
  }

  {
    auto ia = Utils::MemcpyIArchive{Utils::make_span(buf)};
    ParticleRattle out;

    ia >> out;

    BOOST_CHECK_EQUAL(ia.bytes_read(), expected_size);
    check_particle_rattle(out, pr);
  }
}

BOOST_AUTO_TEST_CASE(rattle_constructors) {
  auto pr = ParticleRattle{{1, 2, 3}};

  // check copy constructor
  {
    ParticleRattle out(pr);
    check_particle_rattle(out, pr);
  }

  // check copy assignment operator
  {
    ParticleRattle out; // avoid copy elision
    out = pr;
    check_particle_rattle(out, pr);
  }
}
#endif // BOND_CONSTRAINT

BOOST_AUTO_TEST_CASE(particle_bitfields) {
  auto p = Particle();

  // check default values
  BOOST_CHECK(not p.has_fixed_coordinates());
  BOOST_CHECK(not p.can_rotate());
  BOOST_CHECK(not p.is_fixed_along(1));
  BOOST_CHECK(not p.can_rotate_around(1));

  // check setting of one axis
#ifdef EXTERNAL_FORCES
  p.set_fixed_along(1, true);
  BOOST_CHECK(p.is_fixed_along(1));
  BOOST_CHECK(p.has_fixed_coordinates());
#endif
#ifdef ROTATION
  p.set_can_rotate_around(1, true);
  BOOST_CHECK(p.can_rotate_around(1));
  BOOST_CHECK(p.can_rotate());
#endif

  // check that unsetting is properly registered
#ifdef EXTERNAL_FORCES
  p.set_fixed_along(1, false);
  BOOST_CHECK(not p.has_fixed_coordinates());
#endif
#ifdef ROTATION
  p.set_can_rotate_around(1, false);
  BOOST_CHECK(not p.can_rotate());
#endif

  // check setting of all flags at once
#ifdef ROTATION
  p.set_can_rotate_all_axes();
  BOOST_CHECK(p.can_rotate_around(0));
  BOOST_CHECK(p.can_rotate_around(1));
  BOOST_CHECK(p.can_rotate_around(2));
  p.set_cannot_rotate_all_axes();
  BOOST_CHECK(not p.can_rotate());
#endif
}

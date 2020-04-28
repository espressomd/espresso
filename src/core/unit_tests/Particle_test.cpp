/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <utils/serialization/memcpy_archive.hpp>

#include "Particle.hpp"

BOOST_AUTO_TEST_CASE(comparison) {
  {
    Particle p, q;

    p.identity() = 1;
    q.identity() = 2;

    BOOST_CHECK(p != q);
    BOOST_CHECK(not(p == q));
  }

  {
    Particle p, q;

    p.identity() = 2;
    q.identity() = 2;

    BOOST_CHECK(not(p != q));
    BOOST_CHECK(p == q);
  }
}

BOOST_AUTO_TEST_CASE(serialization) {
  auto p = Particle();

  auto const bond_id = 5;
  auto const bond_partners = std::array<const int, 3>{12, 13, 14};

  std::vector<int> el = {5, 6, 7, 8};

  p.p.identity = 15;
  p.bonds().insert({bond_id, bond_partners});
#ifdef EXCLUSIONS
  p.exclusions() = el;
#endif

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << p;

  boost::archive::text_iarchive in_ar(stream);
  auto q = Particle();
  in_ar >> q;

  BOOST_CHECK(q.p.identity == p.p.identity);
  BOOST_CHECK((*q.bonds().begin() == BondView{bond_id, bond_partners}));

#ifdef EXCLUSIONS
  BOOST_CHECK(q.exclusions() == el);
#endif
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

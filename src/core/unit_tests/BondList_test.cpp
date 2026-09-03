/*
 * Copyright (C) 2020-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE BondList
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BondList.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/range/algorithm/equal.hpp>

#include <algorithm>
#include <array>
#include <iterator>
#include <sstream>

BOOST_AUTO_TEST_CASE(BondView_) {
  /* Dummy values */
  auto const id = 5;
  auto const partners = std::array<const int, 3>{{12, 13, 14}};

  /* BondView can be constructed from an id and a partner range */
  auto const view = BondView{id, partners};
  /* Values are stored and returned */
  BOOST_CHECK_EQUAL(id, view.bond_id());
  BOOST_CHECK_EQUAL(partners.data(), view.partner_ids().data());
  BOOST_CHECK_EQUAL(partners.size(), view.partner_ids().size());

  /* Comparison ops */
  {
    auto const partners_same = partners;
    auto const partners_different = std::array<const int, 3>{{15, 16}};

    BOOST_CHECK((BondView{id, partners} == BondView{id, partners_same}));
    BOOST_CHECK(not(BondView{id, partners} != BondView{id, partners_same}));
    BOOST_CHECK(
        not(BondView{id, partners} == BondView{id, partners_different}));
    BOOST_CHECK(not(BondView{id, partners} == BondView{id + 1, partners_same}));
    /* Primary is the default; a mirror entry with otherwise identical
     * id/partners compares unequal to it. */
    BOOST_CHECK(
        not(BondView{id, partners} == BondView{id, partners_same, false}));
  }

  /* is_primary() reflects the role the view was constructed with */
  BOOST_CHECK((BondView{id, partners}.is_primary()));
  BOOST_CHECK((BondView{id, partners, true}.is_primary()));
  BOOST_CHECK(not(BondView{id, partners, false}.is_primary()));
}

BOOST_AUTO_TEST_CASE(default_ctor) {
  /* BondList can be default constructed and a default constructed BondList
   * is empty. */
  BOOST_CHECK(BondList().empty());
}

BOOST_AUTO_TEST_CASE(Iterator_dereference_) {
  /* Delimiter -6 encodes bond id 2 as a primary entry:
   * -(2 * (bond_id + 1) + role) with role 0 for primary. */
  auto const dummy_bonds = BondList::storage_type{1, 2, -6};
  auto it = BondList::Iterator(dummy_bonds.begin());

  auto const result = *it;
  auto const expected = BondView{2, {dummy_bonds.data(), 2u}};

  BOOST_CHECK(result == expected);
  BOOST_CHECK(result.is_primary());
}

BOOST_AUTO_TEST_CASE(Iterator_dereference_mirror_) {
  /* Delimiter -7 encodes bond id 2 as a mirror entry (role 1):
   * -(2 * (bond_id + 1) + role). */
  auto const dummy_bonds = BondList::storage_type{1, 2, -7};
  auto it = BondList::Iterator(dummy_bonds.begin());

  auto const result = *it;

  BOOST_CHECK_EQUAL(result.bond_id(), 2);
  BOOST_CHECK(not result.is_primary());
  BOOST_CHECK(
      (std::ranges::equal(result.partner_ids(), std::array<int, 2>{{1, 2}})));
}

BOOST_AUTO_TEST_CASE(Iterator_incement_) {
  auto const dummy_bonds = BondList::storage_type{1, 2, -3, 4, -5, 6, -7};
  auto it = BondList::Iterator(dummy_bonds.begin());

  {
    auto const result = std::next(it);
    auto const expected = BondList::Iterator(std::next(dummy_bonds.begin(), 3));
    BOOST_CHECK(expected == result);
  }

  {
    auto const result = std::next(it, 2);
    auto const expected = BondList::Iterator(std::next(dummy_bonds.begin(), 5));
    BOOST_CHECK(expected == result);
  }
}

BOOST_AUTO_TEST_CASE(insert_) {
  /* Dummy values */
  auto const partners = std::array<int, 3>{{1, 2, 3}};
  auto const bond1 = BondView{1, partners};
  auto const bond2 = BondView{2, partners};

  BondList bl;
  /* A bond can be inserted */
  bl.insert(bond1);
  /* BondList is not empty */
  BOOST_CHECK_EQUAL(bl.empty(), false);
  /* The size is increased */
  BOOST_CHECK_EQUAL(bl.size(), 1);
  /* The proper bond is inserted */
  BOOST_CHECK(*bl.begin() == bond1);
  /* A bond can be inserted */
  bl.insert(bond2);
  /* The size is increased */
  BOOST_CHECK_EQUAL(bl.size(), 2);
  /* The first bond is unchanged */
  BOOST_CHECK(*bl.begin() == bond1);
  /* The new bond is inserted */
  BOOST_CHECK(*std::next(bl.begin()) == bond2);
}

BOOST_AUTO_TEST_CASE(insert_mirror_) {
  /* A mirror (non-primary) entry round-trips through insert()/iteration
   * alongside a primary one, distinguishable via is_primary(). */
  auto const partners = std::array<int, 3>{{1, 2, 3}};
  auto const primary = BondView{1, partners, true};
  auto const mirror = BondView{2, partners, false};

  BondList bl;
  bl.insert(primary);
  bl.insert(mirror);

  BOOST_CHECK_EQUAL(bl.size(), 2);
  BOOST_CHECK(*bl.begin() == primary);
  BOOST_CHECK(bl.begin()->is_primary());
  BOOST_CHECK(*std::next(bl.begin()) == mirror);
  BOOST_CHECK(not std::next(bl.begin())->is_primary());
}

BOOST_AUTO_TEST_CASE(erase_) {
  auto const partners = std::array<int, 3>{{1, 2, 3}};
  auto const bond1 = BondView{1, partners};
  auto const bond2 = BondView{2, partners};
  auto const bond3 = BondView{3, partners};

  BondList bl;
  bl.insert(bond1);
  bl.insert(bond2);
  bl.insert(bond3);

  /* Erase at the beginning */
  {
    auto bl_test = bl;

    /* The first bond can be erased */
    bl_test.erase(bl_test.begin());
    /* The size is reduced */
    BOOST_CHECK_EQUAL(bl_test.size(), bl.size() - 1);
    /* The remaining bonds moved one to the front */
    BOOST_CHECK(*bl_test.begin() == *std::next(bl.begin()));
    BOOST_CHECK(*std::next(bl_test.begin()) == *std::next(bl.begin(), 2));
  }

  /* Erase in the middle */
  {
    auto bl_test = bl;

    /* The first bond can be erased */
    bl_test.erase(std::next(bl_test.begin()));
    /* The size is reduced */
    BOOST_CHECK_EQUAL(bl_test.size(), bl.size() - 1);
    /* The first element is unchanged */
    BOOST_CHECK(*bl_test.begin() == *bl.begin());
    /* The last element moved one up */
    BOOST_CHECK(*std::next(bl_test.begin()) == *std::next(bl.begin(), 2));
  }

  /* Erase in the end */
  {
    auto bl_test = bl;

    /* The first bond can be erased */
    bl_test.erase(std::next(bl_test.begin(), 2));
    /* The size is reduced */
    BOOST_CHECK_EQUAL(bl_test.size(), bl.size() - 1);
    /* The other elements are unchanged */
    BOOST_CHECK(*bl_test.begin() == *bl.begin());
    BOOST_CHECK(*std::next(bl_test.begin()) == *std::next(bl.begin()));
  }
}

BOOST_AUTO_TEST_CASE(clear_) {
  auto const partners = std::array<int, 3>{{1, 2, 3}};
  auto const bond1 = BondView{1, partners};
  auto const bond2 = BondView{2, partners};

  BondList bl;
  bl.insert(bond1);
  bl.insert(bond2);

  /* The bond list can be cleared */
  bl.clear();
  /* Afterwards it is empty. */
  BOOST_CHECK(bl.empty());
}

BOOST_AUTO_TEST_CASE(primary_counts_) {
  auto const pair_partners = std::array<int, 1>{{1}};
  auto const angle_partners = std::array<int, 2>{{1, 2}};
  auto const dihedral_partners = std::array<int, 3>{{1, 2, 3}};

  BondList bl;
  BOOST_CHECK_EQUAL(bl.primary_counts().pair, 0);
  BOOST_CHECK_EQUAL(bl.primary_counts().angle, 0);
  BOOST_CHECK_EQUAL(bl.primary_counts().dihedral, 0);

  /* primary entries are counted by arity */
  bl.insert(BondView{1, pair_partners, true});
  bl.insert(BondView{2, angle_partners, true});
  bl.insert(BondView{3, dihedral_partners, true});
  bl.insert(BondView{4, pair_partners, true});
  BOOST_CHECK_EQUAL(bl.primary_counts().pair, 2);
  BOOST_CHECK_EQUAL(bl.primary_counts().angle, 1);
  BOOST_CHECK_EQUAL(bl.primary_counts().dihedral, 1);

  /* mirror entries of any arity do not contribute */
  bl.insert(BondView{5, pair_partners, false});
  bl.insert(BondView{6, angle_partners, false});
  bl.insert(BondView{7, dihedral_partners, false});
  BOOST_CHECK_EQUAL(bl.primary_counts().pair, 2);
  BOOST_CHECK_EQUAL(bl.primary_counts().angle, 1);
  BOOST_CHECK_EQUAL(bl.primary_counts().dihedral, 1);

  /* erasing a primary entry decrements the matching counter */
  auto it = std::find_if(bl.begin(), bl.end(), [](BondView const &b) {
    return b.bond_id() == 1 and b.is_primary();
  });
  bl.erase(it);
  BOOST_CHECK_EQUAL(bl.primary_counts().pair, 1);
  BOOST_CHECK_EQUAL(bl.primary_counts().angle, 1);
  BOOST_CHECK_EQUAL(bl.primary_counts().dihedral, 1);

  /* erasing a mirror entry does not change any counter */
  it = std::find_if(bl.begin(), bl.end(), [](BondView const &b) {
    return b.bond_id() == 5 and not b.is_primary();
  });
  bl.erase(it);
  BOOST_CHECK_EQUAL(bl.primary_counts().pair, 1);
  BOOST_CHECK_EQUAL(bl.primary_counts().angle, 1);
  BOOST_CHECK_EQUAL(bl.primary_counts().dihedral, 1);

  /* copy and move preserve the counts */
  auto const bl_copy = bl;
  BOOST_CHECK_EQUAL(bl_copy.primary_counts().pair, 1);
  BOOST_CHECK_EQUAL(bl_copy.primary_counts().angle, 1);
  BOOST_CHECK_EQUAL(bl_copy.primary_counts().dihedral, 1);

  auto bl_move_src = bl;
  BondList bl_moved;
  bl_moved = std::move(bl_move_src);
  BOOST_CHECK_EQUAL(bl_moved.primary_counts().pair, 1);
  BOOST_CHECK_EQUAL(bl_moved.primary_counts().angle, 1);
  BOOST_CHECK_EQUAL(bl_moved.primary_counts().dihedral, 1);

  /* clear resets the counts */
  auto bl_cleared = bl;
  bl_cleared.clear();
  BOOST_CHECK_EQUAL(bl_cleared.primary_counts().pair, 0);
  BOOST_CHECK_EQUAL(bl_cleared.primary_counts().angle, 0);
  BOOST_CHECK_EQUAL(bl_cleared.primary_counts().dihedral, 0);
}

BOOST_AUTO_TEST_CASE(primary_counts_serialization_) {
  auto const pair_partners = std::array<int, 1>{{1}};
  auto const angle_partners = std::array<int, 2>{{1, 2}};

  BondList bl;
  bl.insert(BondView{1, pair_partners, true});
  bl.insert(BondView{2, angle_partners, true});
  bl.insert(BondView{3, pair_partners, false});

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  out_ar << bl;

  boost::archive::text_iarchive in_ar(stream);
  BondList bl_restored;
  in_ar >> bl_restored;

  /* deserialization recomputes primary_counts() to match the original */
  BOOST_CHECK_EQUAL(bl_restored.primary_counts().pair, bl.primary_counts().pair);
  BOOST_CHECK_EQUAL(bl_restored.primary_counts().angle,
                    bl.primary_counts().angle);
  BOOST_CHECK_EQUAL(bl_restored.primary_counts().dihedral,
                    bl.primary_counts().dihedral);
}

namespace {
/**
 * @brief Stand-in for the on-disk layout of a pre-versioning (class
 * version 0) BondList, as found in checkpoints written before mirror
 * entries were introduced. Serializes the same way BondList does (a
 * size followed by the raw storage array), but such archives only ever
 * contain primary entries encoded with the legacy delimiter
 * @c -(bond_id+1) -- one bit narrower than the current
 * @c -(2*(bond_id+1)+role) encoding, since the role bit did not exist
 * yet.
 */
struct LegacyBondList {
  BondList::storage_type m_storage;

  template <class Archive>
  void serialize(Archive &ar, unsigned int const /* version */) {
    if (Archive::is_loading::value) {
      std::size_t size{};
      ar & size;
      m_storage.resize(size);
    }
    if (Archive::is_saving::value) {
      auto size = m_storage.size();
      ar & size;
    }
    ar &boost::serialization::make_array(m_storage.data(), m_storage.size());
  }
};
} // namespace
BOOST_CLASS_VERSION(LegacyBondList, 0)

BOOST_AUTO_TEST_CASE(legacy_archive_migration_) {
  /* Hand-build a version-0 encoded bond list with two primary pair
   * bonds: id 2 with partner 7 (legacy delimiter -(2+1) = -3), and id 5
   * with partner 9 (legacy delimiter -(5+1) = -6). */
  LegacyBondList legacy;
  legacy.m_storage = BondList::storage_type{7, -3, 9, -6};

  std::stringstream stream;
  {
    boost::archive::text_oarchive out_ar(stream);
    out_ar << legacy;
  }

  BondList bl;
  {
    /* Loaded as a current-version BondList: serialize() must detect
     * version < 1 and migrate the legacy delimiters on the fly. */
    boost::archive::text_iarchive in_ar(stream);
    in_ar >> bl;
  }

  BOOST_REQUIRE_EQUAL(bl.size(), 2u);
  auto it = bl.begin();
  BOOST_CHECK_EQUAL(it->bond_id(), 2);
  BOOST_CHECK(it->is_primary());
  BOOST_CHECK(
      (std::ranges::equal(it->partner_ids(), std::array<int, 1>{{7}})));
  ++it;
  BOOST_CHECK_EQUAL(it->bond_id(), 5);
  BOOST_CHECK(it->is_primary());
  BOOST_CHECK(
      (std::ranges::equal(it->partner_ids(), std::array<int, 1>{{9}})));

  /* primary_counts(), populated only on load (recompute_primary_counts()),
   * must reflect the migrated bonds. */
  BOOST_CHECK_EQUAL(bl.primary_counts().pair, 2);
  BOOST_CHECK_EQUAL(bl.primary_counts().angle, 0);
  BOOST_CHECK_EQUAL(bl.primary_counts().dihedral, 0);
}

BOOST_AUTO_TEST_CASE(serialization_) {
  auto const partners = std::array<int, 3>{{4, 5, 6}};
  auto const bond1 = BondView{1, partners};
  auto const bond2 = BondView{2, partners};

  BondList bl;
  bl.insert(bond1);
  bl.insert(bond2);

  std::stringstream stream;
  boost::archive::text_oarchive out_ar(stream);
  /* BondList can be serialized */
  out_ar << bl;

  {
    boost::archive::text_iarchive in_ar(stream);
    BondList bl_restored;
    /* BondList can be deserialized */
    in_ar >> bl_restored;

    /* BondList is correctly restored */
    BOOST_CHECK(boost::equal(bl, bl_restored));
  }
}

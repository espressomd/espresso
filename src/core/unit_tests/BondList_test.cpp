/*
 * Copyright (C) 2020 The ESPResSo project
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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <sstream>

#include "BondList.hpp"

BOOST_AUTO_TEST_CASE(BondView_) {
  /* Dummy values */
  auto const id = 5;
  auto const partners = std::array<const int, 3>{12, 13, 14};

  /* BondView can be constructed from an id and a partner range */
  auto const view = BondView{id, partners};
  /* Values are stored and returned */
  BOOST_CHECK_EQUAL(id, view.bond_id());
  BOOST_CHECK_EQUAL(partners.data(), view.partner_ids().data());
  BOOST_CHECK_EQUAL(partners.size(), view.partner_ids().size());

  /* Comparison ops */
  {
    auto const partners_same = partners;
    auto const partners_different = std::array<const int, 3>{15, 16};

    BOOST_CHECK((BondView{id, partners} == BondView{id, partners_same}));
    BOOST_CHECK(not(BondView{id, partners} != BondView{id, partners_same}));
    BOOST_CHECK(
        not(BondView{id, partners} == BondView{id, partners_different}));
    BOOST_CHECK(not(BondView{id, partners} == BondView{id + 1, partners_same}));
  }
}

BOOST_AUTO_TEST_CASE(default_ctor) {
  /* BondList can be default constructed and a default constructed BondList
   * is empty. */
  BOOST_CHECK(BondList().empty());
}

BOOST_AUTO_TEST_CASE(Iterator_dereference_) {
  auto const dummy_bonds = BondList::storage_type{1, 2, -3};
  auto it = BondList::Iterator(dummy_bonds.begin());

  auto const result = *it;
  auto const expected = BondView{2, {dummy_bonds.data(), 2u}};

  BOOST_CHECK(result == expected);
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
  auto const partners = std::array<int, 3>{1, 2, 3};
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

BOOST_AUTO_TEST_CASE(erase_) {
  auto const partners = std::array<int, 3>{1, 2, 3};
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
  auto const partners = std::array<int, 3>{1, 2, 3};
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

BOOST_AUTO_TEST_CASE(serialization_) {
  auto const partners = std::array<int, 3>{4, 5, 6};
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

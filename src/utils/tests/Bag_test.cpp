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
#define BOOST_TEST_MODULE Utils::Bag
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Bag.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/range/algorithm/find.hpp>

#include <memory>
#include <sstream>

BOOST_AUTO_TEST_CASE(constructor_) {
  /* A Bag can be default constructed */
  auto const bag = Utils::Bag<int>();
  /* A default constructed Bag is empty */
  BOOST_CHECK(bag.empty());
}

BOOST_AUTO_TEST_CASE(insert_) {
  /* Copy insert */
  {
    auto const elements = std::array<int, 3>{1, 2, 3};
    auto bag = Utils::Bag<int>();

    /* Elements can be inserted into the bag */
    for (auto e : elements) {
      bag.insert(e);
    }

    /* The size of the bag has increased */
    BOOST_REQUIRE_EQUAL(bag.size(), elements.size());

    /* The elements are in the bag */
    for (auto e : elements) {
      BOOST_CHECK(boost::find(bag, e) != bag.end());
    }
  }

  /* Move insert */
  {
    using MoveOnly = std::unique_ptr<int>;
    auto bag = Utils::Bag<MoveOnly>();

    /* An element can be moved into the bag */
    bag.insert(std::make_unique<int>(5));
    /* The size of the bag can be increased */
    BOOST_REQUIRE_EQUAL(bag.size(), 1);
    /* The element is in the bag */
    BOOST_CHECK_EQUAL(**bag.begin(), 5);
  }
}

BOOST_AUTO_TEST_CASE(erase_) {
  auto const elements = std::array<int, 3>{1, 2, 3};

  {
    /* Given a bag with elements */
    auto bag = Utils::Bag<int>();
    for (auto e : elements) {
      bag.insert(e);
    }

    /* ... the first element can be erased */
    auto const it = bag.erase(bag.begin());
    /* the begin iterator is returned */
    BOOST_CHECK(it == bag.begin());
    /* and the other elements are still in the bag */
    BOOST_CHECK(boost::find(bag, elements[1]) != bag.end());
    BOOST_CHECK(boost::find(bag, elements[2]) != bag.end());
  }

  {
    /* Given a bag with elements */
    auto bag = Utils::Bag<int>();
    for (auto e : elements) {
      bag.insert(e);
    }

    /* ... a middle element can be erased */
    auto const it = bag.erase(bag.begin() + 1);
    /* the correct iterator is returned */
    BOOST_CHECK(it == bag.begin() + 1);
    /* and the other elements are still in the bag */
    BOOST_CHECK(boost::find(bag, elements[0]) != bag.end());
    BOOST_CHECK(boost::find(bag, elements[2]) != bag.end());
  }

  {
    /* Given a bag with elements */
    auto bag = Utils::Bag<int>();
    for (auto e : elements) {
      bag.insert(e);
    }

    /* ... the last element can be erased */
    auto const it = bag.erase(bag.end() - 1);
    /* the correct iterator is returned */
    BOOST_CHECK(it == bag.end());
    /* and the other elements are still in the bag */
    BOOST_CHECK(boost::find(bag, elements[0]) != bag.end());
    BOOST_CHECK(boost::find(bag, elements[1]) != bag.end());
  }
}

BOOST_AUTO_TEST_CASE(size_) {
  auto const elements = std::array<int, 5>{1, 2, 3, 5, 6};

  /* Given a bag with elements */
  auto bag = Utils::Bag<int>();
  for (auto e : elements) {
    bag.insert(e);
  }

  /* The size is the number of elements in the bag. */
  BOOST_CHECK_EQUAL(bag.size(), elements.size());
}

BOOST_AUTO_TEST_CASE(iterator_range_) {
  auto const elements = std::array<int, 5>{1, 2, 3, 5, 6};

  /* Given a bag with elements */
  auto bag = Utils::Bag<int>();
  for (auto e : elements) {
    bag.insert(e);
  }

  /* The range of the non-const iterators spans all elements */
  for (auto const &e : elements) {
    BOOST_CHECK(boost::find(bag, e) != bag.end());
  }
  /* The range of the const iterators spans all elements */
  for (auto const &e : elements) {
    BOOST_CHECK(boost::find(const_cast<const Utils::Bag<int> &>(bag), e) !=
                bag.end());
  }
}

BOOST_AUTO_TEST_CASE(empty_) {
  /* Given a bag with elements */
  auto bag = Utils::Bag<int>();
  bag.insert(5);
  /* It is not empty */
  BOOST_CHECK(not bag.empty());
}

BOOST_AUTO_TEST_CASE(clear_) {
  /* Given a bag with elements */
  auto bag = Utils::Bag<int>();
  for (auto e : {1, 2, 3}) {
    bag.insert(e);
  }

  /* the bag can be cleared */
  bag.clear();
  /* the cleared bag is empty */
  BOOST_CHECK(bag.empty());
}

BOOST_AUTO_TEST_CASE(reserve_) {
  /* Given a bag */
  auto bag = Utils::Bag<int>();

  /* Memory can be reserved */
  const size_t size = 11;
  bag.reserve(size);
  /* and the capacity of the bag is at least
   * the reserved size */
  BOOST_CHECK_GE(bag.capacity(), size);
}

BOOST_AUTO_TEST_CASE(resize_) {
  auto const elements = std::array<int, 5>{1, 2, 3, 5, 6};

  /* Given a bag with elements */
  auto bag = Utils::Bag<int>();
  for (auto e : elements) {
    bag.insert(e);
  }

  /* It can be resized */
  const size_t size = 11;
  bag.resize(size);
  /* The size matches */
  BOOST_CHECK_EQUAL(bag.size(), size);
  /* All the elements are still in the bag */
  for (auto const &e : bag) {
    BOOST_CHECK(boost::find(bag, e) != bag.end());
  }
}

BOOST_AUTO_TEST_CASE(swap_) {
  auto const elements1 = std::array<int, 5>{1, 2, 3};
  auto const elements2 = std::array<int, 5>{1, 2, 3};

  /* Given two bags with elements */
  auto bag1 = Utils::Bag<int>();
  for (auto e : elements1) {
    bag1.insert(e);
  }

  auto bag2 = Utils::Bag<int>();
  for (auto e : elements2) {
    bag2.insert(e);
  }

  /* they can be swapped. */
  swap(bag1, bag2);

  /* The elements are swapped */
  BOOST_CHECK_EQUAL(bag2.size(), elements1.size());
  for (auto const &e : elements1) {
    BOOST_CHECK(boost::find(bag2, e) != bag2.end());
  }

  BOOST_CHECK_EQUAL(bag1.size(), elements2.size());
  for (auto const &e : elements2) {
    BOOST_CHECK(boost::find(bag1, e) != bag1.end());
  }
}

BOOST_AUTO_TEST_CASE(serialize_) {
  auto const elements = std::array<int, 5>{1, 2, 3, 5, 6};

  /* Given a bag with elements */
  auto bag = Utils::Bag<int>();
  for (auto e : elements) {
    bag.insert(e);
  }

  /* The bag can be serialized */
  std::stringstream stream;
  boost::archive::text_oarchive(stream) << bag;

  /* It can be deserialized */
  Utils::Bag<int> restored_bag;
  boost::archive::text_iarchive(stream) >> restored_bag;

  /* The deserialized object contains the same elements */
  BOOST_CHECK(boost::equal(bag, restored_bag));
}
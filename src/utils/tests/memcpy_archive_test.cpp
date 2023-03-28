/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE memcpy archive test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/* This <boost/serialization/version.hpp> include guards against an issue
 * in boost::serialization from boost 1.74.0 that leads to compiler error
 * "explicit specialization of undeclared template struct 'version'" when
 * including <boost/serialization/optional.hpp>. More details in tickets:
 * https://github.com/boostorg/serialization/issues/210
 * https://github.com/boostorg/serialization/issues/217
 */
#include <boost/serialization/version.hpp>

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/serialization/memcpy_archive.hpp>
#include <utils/type_traits.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/optional.hpp>

#include <array>
#include <cstddef>
#include <type_traits>

struct NonTrivial {
  boost::optional<Utils::Vector3d> ov;

  template <class Archive> void serialize(Archive &ar, long int) { ar &ov; }
};

using OpVec = boost::optional<Utils::Vector3d>;

namespace Utils {
template <> struct is_statically_serializable<NonTrivial> : std::true_type {};

template <class T>
struct is_statically_serializable<boost::optional<T>>
    : is_statically_serializable<T> {};
} // namespace Utils

BOOST_AUTO_TEST_CASE(packing_size_test) {
  BOOST_CHECK_EQUAL(Utils::MemcpyIArchive::packing_size<int>(), sizeof(int));
  BOOST_CHECK_EQUAL(Utils::MemcpyIArchive::packing_size<NonTrivial>(),
                    Utils::MemcpyIArchive::packing_size<OpVec>());
}

BOOST_AUTO_TEST_CASE(type_traits) {
  static_assert(Utils::is_statically_serializable<int>::value);
  static_assert(Utils::detail::use_memcpy<int>::value);
  static_assert(not Utils::detail::use_serialize<int>::value);

  static_assert(Utils::is_statically_serializable<OpVec>::value);
  static_assert(not Utils::detail::use_memcpy<OpVec>::value);
  static_assert(Utils::detail::use_serialize<OpVec>::value);

  BOOST_TEST_PASSPOINT();
}

BOOST_AUTO_TEST_CASE(skipping_and_position) {
  std::array<char, 10> buf;

  auto ar = Utils::MemcpyOArchive(Utils::make_span(buf));

  BOOST_CHECK_EQUAL(0, ar.bytes_processed());
  ar.skip(5);
  BOOST_CHECK_EQUAL(5, ar.bytes_processed());
}

BOOST_AUTO_TEST_CASE(memcpy_processing) {
  std::array<char, 10> buf;

  auto const test_number = 5;

  {
    auto oa = Utils::MemcpyOArchive(Utils::make_span(buf));
    oa << test_number;
    BOOST_CHECK_EQUAL(oa.bytes_written(), sizeof(test_number));
    BOOST_CHECK_EQUAL(oa.get_library_version(), 4);
  }

  {
    auto ia = Utils::MemcpyIArchive(Utils::make_span(buf));
    int out;
    ia >> out;
    BOOST_CHECK_EQUAL(out, test_number);
    BOOST_CHECK_EQUAL(ia.bytes_read(), sizeof(test_number));
    BOOST_CHECK_EQUAL(ia.get_library_version(), 4);
  }
}

BOOST_AUTO_TEST_CASE(serializaton_processing) {
  std::array<char, 2 * sizeof(OpVec)> buf;

  const OpVec active = Utils::Vector3d{1., 2., 3.};
  const OpVec inactive = boost::none;
  {
    auto oa = Utils::MemcpyOArchive{Utils::make_span(buf)};
    auto in1 = active;
    auto in2 = inactive;
    oa << in1;
    oa << in2;

    BOOST_CHECK_EQUAL(oa.bytes_written(), buf.size());
  }

  {
    auto ia = Utils::MemcpyIArchive{Utils::make_span(buf)};
    OpVec out1 = Utils::Vector3d{}, out2;
    ia >> out1;
    ia >> out2;

    BOOST_CHECK_EQUAL(ia.bytes_read(), buf.size());

    BOOST_CHECK(out1 == active);
    BOOST_CHECK(out2 == inactive);
  }
}

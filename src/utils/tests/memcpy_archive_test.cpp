/*
 * Copyright (C) 2019 The ESPResSo project
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

#include "utils/serialization/memcpy_archive.hpp"

#include <boost/optional.hpp>
#include <utils/Vector.hpp>

#include <vector>

struct NonTrivial {
  boost::optional<Utils::Vector3d> ov;

  template <class Archive> void serialize(Archive &ar, long int) { ar &ov; }
};

using OpVec = boost::optional<Utils::Vector3d>;

namespace Utils {
template <> struct is_statically_serializable<NonTrivial> : std::true_type {};
} // namespace Utils

BOOST_AUTO_TEST_CASE(packing_size_test) {
  BOOST_CHECK_EQUAL(Utils::MemcpyIArchive::packing_size<int>(), sizeof(int));

  {
    BOOST_CHECK_EQUAL(Utils::MemcpyIArchive::packing_size<NonTrivial>(),
                      Utils::MemcpyIArchive::packing_size<OpVec>());
  }
}

BOOST_AUTO_TEST_CASE(optional_packing) {
  std::vector<char> buf(2 * Utils::MemcpyOArchive::packing_size<OpVec>());

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
    OpVec out1, out2;
    ia >> out1;
    ia >> out2;

    BOOST_CHECK_EQUAL(ia.bytes_read(), buf.size());

    BOOST_CHECK(out1 == active);
    BOOST_CHECK(out2 == inactive);
  }
}

BOOST_AUTO_TEST_CASE(non_trivial_packing) {
  static_assert(not std::is_trivially_copyable<NonTrivial>::value, "");

  std::vector<char> buf(2 * Utils::MemcpyOArchive::packing_size<NonTrivial>());

  auto const active = NonTrivial{Utils::Vector3d{1., 2., 3.}};
  auto const inactive = NonTrivial{boost::none};

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
    NonTrivial out1, out2;
    ia >> out1;
    ia >> out2;

    BOOST_CHECK_EQUAL(ia.bytes_read(), buf.size());

    BOOST_CHECK(out1.ov == active.ov);
    BOOST_CHECK(out2.ov == inactive.ov);
  }
}
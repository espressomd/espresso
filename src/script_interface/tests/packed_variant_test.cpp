/*
 * Copyright (C) 2020-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE PackedVariant test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <script_interface/packed_variant.hpp>

#include <unordered_map>
#include <vector>

BOOST_AUTO_TEST_CASE(object_id_) {
  using ScriptInterface::object_id;

  auto const p1 = reinterpret_cast<ScriptInterface::ObjectHandle *>(1);
  auto const p2 = reinterpret_cast<ScriptInterface::ObjectHandle *>(2);

  /* Check that different objects have different ids */
  BOOST_CHECK_NE(object_id(p1), object_id(p2));
}

BOOST_AUTO_TEST_CASE(PackVisitor_) {
  using ScriptInterface::ObjectId;
  using ScriptInterface::ObjectRef;
  using ScriptInterface::PackedVariant;
  using ScriptInterface::Variant;

  const Variant v = std::vector<Variant>{{5, std::vector<Variant>{
                                                 6,
                                                 ObjectRef{},
                                             }}};

  auto const visitor = ScriptInterface::PackVisitor();

  auto const result = boost::apply_visitor(visitor, v);

  const PackedVariant expected =
      std::vector<PackedVariant>{{5, std::vector<PackedVariant>{
                                         6,
                                         object_id(ObjectRef{}.get()),
                                     }}};

  BOOST_CHECK(result == expected);

  /* Check that the object hast been captured. */
  BOOST_CHECK(visitor.objects().at(object_id(ObjectRef{}.get())) ==
              ObjectRef{});
}

BOOST_AUTO_TEST_CASE(pack_) {
  using ScriptInterface::ObjectRef;
  using ScriptInterface::Variant;

  const Variant v = std::vector<Variant>{{5, std::vector<Variant>{
                                                 6,
                                                 ObjectRef{},
                                             }}};

  auto const expected = boost::apply_visitor(ScriptInterface::PackVisitor(), v);
  auto const result = ScriptInterface::pack(v);

  BOOST_CHECK(expected == result);
}

BOOST_AUTO_TEST_CASE(UnpackVisitor_) {
  using ScriptInterface::ObjectId;
  using ScriptInterface::ObjectRef;
  using ScriptInterface::PackedVariant;
  using ScriptInterface::Variant;

  const PackedVariant v =
      std::vector<PackedVariant>{{5, std::vector<PackedVariant>{
                                         6,
                                         object_id(ObjectRef{}.get()),
                                     }}};

  std::unordered_map<ObjectId, ObjectRef> const &objects{
      {object_id(ObjectRef{}.get()), ObjectRef{}}};
  auto const visitor = ScriptInterface::UnpackVisitor(objects);

  auto const result = boost::apply_visitor(visitor, v);

  const Variant expected = std::vector<Variant>{{5, std::vector<Variant>{
                                                        6,
                                                        ObjectRef{},
                                                    }}};

  BOOST_CHECK(result == expected);
}

BOOST_AUTO_TEST_CASE(unpack_) {
  using ScriptInterface::ObjectId;
  using ScriptInterface::ObjectRef;
  using ScriptInterface::PackedVariant;

  const PackedVariant v =
      std::vector<PackedVariant>{{5, std::vector<PackedVariant>{
                                         6,
                                         object_id(ObjectRef{}.get()),
                                     }}};

  std::unordered_map<ObjectId, ObjectRef> const &objects{
      {object_id(ObjectRef{}.get()), ObjectRef{}}};

  auto const expected =
      boost::apply_visitor(ScriptInterface::UnpackVisitor(objects), v);
  auto const result = ScriptInterface::unpack(v, objects);

  BOOST_CHECK(expected == result);
}

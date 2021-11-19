/*
 * Copyright (C) 2021 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE ObjectList test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/range/algorithm/find.hpp>

#include "script_interface/LocalContext.hpp"
#include "script_interface/ObjectList.hpp"

#include "core/communication.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <memory>
#include <vector>

using namespace ScriptInterface;

struct ObjectListImpl : ObjectList<ObjectHandle> {
  std::vector<ObjectRef> mock_core;

private:
  void add_in_core(const ObjectRef &obj_ptr) override {
    mock_core.push_back(obj_ptr);
  }
  void remove_in_core(const ObjectRef &obj_ptr) override {
    mock_core.erase(std::remove(mock_core.begin(), mock_core.end(), obj_ptr),
                    mock_core.end());
  }
};

BOOST_AUTO_TEST_CASE(default_construction) {
  // A defaulted ObjectList has no elements.
  BOOST_CHECK(ObjectListImpl{}.elements().empty());
}

BOOST_AUTO_TEST_CASE(adding_elements) {
  // Added elements are on the back of the list of elements.
  auto e = ObjectRef{};
  ObjectListImpl list;
  list.add(e);
  BOOST_CHECK(list.elements().back() == e);
  // And is added to the core
  BOOST_CHECK(boost::find(list.mock_core, e) != list.mock_core.end());
}

BOOST_AUTO_TEST_CASE(removing_elements) {
  // An element that is removed from the list is
  // no longer an element of the list.
  auto e = ObjectRef{};
  ObjectListImpl list;
  list.add(e);
  list.remove(e);
  BOOST_CHECK(boost::find(list.elements(), e) == list.elements().end());
  // And is removed from the core
  BOOST_CHECK(boost::find(list.mock_core, e) == list.mock_core.end());
}

BOOST_AUTO_TEST_CASE(clearing_elements) {
  // A cleared list is empty.
  ObjectListImpl list;
  list.add(ObjectRef{});
  list.add(ObjectRef{});
  list.clear();
  BOOST_CHECK(list.elements().empty());
  BOOST_CHECK(list.mock_core.empty());
}

BOOST_AUTO_TEST_CASE(serialization) {
  // In a context
  Utils::Factory<ObjectHandle> f;
  f.register_new<ObjectHandle>("ObjectHandle");
  f.register_new<ObjectListImpl>("ObjectList");
  auto ctx = std::make_shared<LocalContext>(f, 0);
  // A list of some elements
  auto list = std::dynamic_pointer_cast<ObjectListImpl>(
      ctx->make_shared("ObjectList", {}));
  // with a bunch of elements

  list->add(ctx->make_shared("ObjectHandle", {}));
  list->add(ctx->make_shared("ObjectHandle", {}));
  // can be serialized to a string
  auto const s = list->serialize();
  // and can be restored to a list with the same elements
  auto const list2 = std::dynamic_pointer_cast<ObjectListImpl>(
      ObjectHandle::deserialize(s, *ctx));
  BOOST_CHECK(list2->elements().size() == list->elements().size());
  BOOST_CHECK(list2->elements().front()->name() == "ObjectHandle");
  BOOST_CHECK(list2->elements().back()->name() == "ObjectHandle");
  // and the elements are restored to the core
  BOOST_CHECK(list2->mock_core.size() == 2);
  BOOST_CHECK(list2->mock_core.front()->name() == "ObjectHandle");
  BOOST_CHECK(list2->mock_core.back()->name() == "ObjectHandle");
}

int main(int argc, char **argv) {
  auto mpi_env = std::make_shared<boost::mpi::environment>(argc, argv);
  Communication::init(mpi_env);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

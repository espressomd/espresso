/*
 * Copyright (C) 2022 The ESPResSo project
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
#define BOOST_TEST_MODULE Object container MPI guard test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/ObjectList.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using ScriptInterface::ObjectHandle;

namespace Testing {
struct ObjectContainer : ScriptInterface::ObjectList<ObjectHandle> {
  std::vector<std::shared_ptr<ObjectHandle>> objects;

private:
  void add_in_core(std::shared_ptr<ObjectHandle> const &obj_ptr) override {
    objects.push_back(obj_ptr);
  }
  void remove_in_core(std::shared_ptr<ObjectHandle> const &obj_ptr) override {
    objects.erase(std::remove(objects.begin(), objects.end(), obj_ptr),
                  objects.end());
  }
};
} // namespace Testing

BOOST_AUTO_TEST_CASE(parallel_exception) {
  boost::mpi::communicator world;
  auto const obj_ptr = std::make_shared<ScriptInterface::ObjectHandle>();
  auto const predicate = [](std::exception const &ex) {
    std::string message =
        "Non-empty object containers do not support checkpointing "
        "in MPI environments. Container  contains 1 elements.";
    return ex.what() == message;
  };

  Testing::ObjectContainer list;
  BOOST_CHECK_NO_THROW(list.serialize());

  list.add(obj_ptr);
  if (world.size() > 1) {
    BOOST_CHECK_EXCEPTION(list.serialize(), std::runtime_error, predicate);
  }

  list.remove(obj_ptr);
  BOOST_CHECK_NO_THROW(list.serialize());
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

/*
 * Copyright (C) 2025 The ESPResSo project
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

#define BOOST_TEST_MODULE particle_reduction test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "particle_node.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <cassert>
#include <functional>

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    auto system = System::System::create();
    system->set_box_l(Utils::Vector3d{10., 10., 10.});
    system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(system);
    ::make_new_particle(0, Utils::Vector3d{0., 0., 0.});
    auto *const p = system->cell_structure->get_local_particle(0);
    assert(p);
    p->v() = Utils::Vector3d{3., 0., 0.};
#ifdef ESPRESSO_MASS
    p->mass() = 2.;
#endif
  }
  ~GlobalConfig() { ::System::reset_system(); }
};

// Decorator to skip tests if shared-memory parallelism isn't compiled in
boost::test_tools::assertion_result has_shm(boost::unit_test::test_unit_id) {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  return true;
#else
  return false;
#endif
}

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

auto const reduce_op = []<typename T>(T &a, T const &b) { a = a + b; };

BOOST_TEST_DECORATOR(*boost::unit_test::precondition(has_shm))
BOOST_AUTO_TEST_CASE(test_make_kokkos_reduction) {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  auto &system = System::get_system();
  auto const &cell_structure = *system.cell_structure;
  auto const &cells = cell_structure.decomposition().local_cells();

  auto const kernel = [&cells](int i, Utils::Vector3d &res) {
    for (auto const &p : cells[i]->particles()) {
      res += p.mass() * p.v();
    }
  };
  auto const reducer =
      Reduction::make_kokkos_reducer<Utils::Vector3d>(kernel, reduce_op);
  for (auto i = 0; i < cells.size(); ++i) {
    auto ref = Utils::Vector3d{0., 0., 0.};
    kernel(i, ref);
    auto res = Utils::Vector3d{0., 0., 0.};
    reducer(i, res);
    // The lambda `kernel(i, res)` is executed inside `reducer(i, res)`,
    // so make sure that both results are equal.
    BOOST_CHECK_EQUAL(res, ref);
  }
#endif // ESPRESSO_SHARED_MEMORY_PARALLELISM
}

BOOST_AUTO_TEST_CASE(test_reduce_over_local_particles) {
  auto &system = System::get_system();
  auto const &cell_structure = *system.cell_structure;
  auto const *const p = cell_structure.get_local_particle(0);
  assert(p);

  auto const kernel = [](Utils::Vector3d &acc, Particle const &p) {
    acc += p.mass() * p.v();
  };
  auto const ref = p->mass() * p->v();
  auto const res = reduce_over_local_particles<Utils::Vector3d>(
      cell_structure, kernel, reduce_op);
  BOOST_CHECK_EQUAL(res, ref);
}

BOOST_AUTO_TEST_SUITE_END()

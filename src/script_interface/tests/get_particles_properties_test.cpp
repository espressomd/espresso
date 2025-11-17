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

#define BOOST_TEST_MODULE "particle property parallel getters"
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "script_interface/Context.hpp"
#include "script_interface/LocalContext.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/particle_data/ParticleSlice.hpp"

#include "core/Particle.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/unit_tests/EspressoCoreGlobalConfig.hpp"
#include "core/unit_tests/ParticleFactory.hpp"

#include <utils/Factory.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi/graph_communicator.hpp>
#include <boost/test/tools/old/interface.hpp>

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

using namespace ScriptInterface;

namespace espresso {
// ESPResSo system instance
static std::shared_ptr<::System::System> system;
} // namespace espresso

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    espresso::system = ::System::System::create();
    espresso::system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(espresso::system);
  }
  ~GlobalConfig() {
    espresso::system.reset();
    ::System::reset_system();
  }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

BOOST_FIXTURE_TEST_CASE(test_get_particles_properties, ParticleFactory) {
  // particle ids and types
  std::vector<int> pids{0, 1, 2, 3, 4};
  std::vector<int> expected_types{1, 0, 3, 2, 2};

  Utils::Factory<ObjectHandle> f;
  boost::mpi::communicator comm;
  auto ctx = std::make_shared<LocalContext>(f, comm);

  auto const get_type = [](Particle const &p) { return p.type(); };

  // Create particle core objects
  create_particle(Utils::Vector3d{0., 0., 0.}, 1, 0);
  create_particle(Utils::Vector3d{0., 0., 0.}, 0, 1);
  create_particle(Utils::Vector3d{0., 0., 0.}, 2, 3);
  create_particle(Utils::Vector3d{0., 0., 0.}, 4, 2);
  create_particle(Utils::Vector3d{0., 0., 0.}, 3, 2);

  auto const result = Particles::get_particles_properties<int>(
      pids, get_type, ctx.get(), *espresso::system->cell_structure);

  if (comm.rank() == 0) {
    BOOST_CHECK(std::ranges::equal(result, expected_types));
  }
}

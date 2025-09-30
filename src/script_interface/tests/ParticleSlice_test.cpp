/*
 * Copyright (C) 2016-2025 The ESPResSo project
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

#define BOOST_TEST_MODULE ParticleSlice test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/Context.hpp"
#include "script_interface/Exception.hpp"
#include "script_interface/LocalContext.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/particle_data/ParticleHandle.hpp"
#include "script_interface/particle_data/ParticleSlice.hpp"

#include "core/unit_tests/EspressoCoreGlobalConfig.hpp"
#include "core/unit_tests/ParticleFactory.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <cassert>
#include <memory>
#include <optional>
#include <string>
#include <utility>
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

BOOST_FIXTURE_TEST_CASE(set_slice_positions, ParticleFactory) {
  // New positions for slice
  Variant new_pos1 = Utils::Vector3d{1., 0., 0.};
  Variant new_pos2 = Utils::Vector3d{0., 1., 0.};
  std::vector<Variant> new_positions = {new_pos1, new_pos2};

  Utils::Factory<ObjectHandle> f;
  f.register_new<Particles::ParticleSlice>("ParticleSlice");
  f.register_new<Particles::ParticleHandle>("ParticleHandle");
  f.register_new<CellSystem::CellSystem>("CellSystem");

  boost::mpi::communicator comm;
  auto ctx = std::make_shared<LocalContext>(f, comm);

  // Create TestCellSystem object
  auto &&sp_cell_system{ctx->make_shared("CellSystem", {})};
  std::shared_ptr<CellSystem::CellSystem> cell_system =
      std::dynamic_pointer_cast<CellSystem::CellSystem>(sp_cell_system);

  // Parameters for ScriptInterface objects
  VariantMap p1_init_params{{"id", 0}, {"__cell_structure", cell_system}};
  VariantMap p2_init_params{{"id", 1}, {"__cell_structure", cell_system}};
  VariantMap p_slice_init_params{{"id_selection", std::vector<int>{0, 1}},
                                 {"__cell_structure", cell_system}};

  // Create particle core objects
  create_particle(Utils::Vector3d{0., 0., 0.}, 0, 0);
  create_particle(Utils::Vector3d{0., 0., 0.}, 1, 0);

  // Create 2 ParticleHandle, and 1 ParticleSlice object
  auto &&sp_p1_handle{ctx->make_shared("ParticleHandle", {p1_init_params})};
  auto &&sp_p2_handle{ctx->make_shared("ParticleHandle", {p2_init_params})};
  // Particle handle 1
  auto p1_handle =
      std::dynamic_pointer_cast<Particles::ParticleHandle>(sp_p1_handle);
  // Particle handle 2
  auto p2_handle =
      std::dynamic_pointer_cast<Particles::ParticleHandle>(sp_p2_handle);
  // Attach to system
  p1_handle->attach(espresso::system);
  p2_handle->attach(espresso::system);

  // Particle slice
  auto &&sp_p_slice{ctx->make_shared("ParticleSlice", {p_slice_init_params})};
  auto p_slice =
      std::dynamic_pointer_cast<Particles::ParticleSlice>(sp_p_slice);

  // Check if old positions are actually zero
  auto pos1_before_update =
      get_value<Utils::Vector3d>(p1_handle->get_parameter("pos"));
  auto pos2_before_update =
      get_value<Utils::Vector3d>(p2_handle->get_parameter("pos"));
  Utils::Vector3d zero_vector{};
  BOOST_TEST(pos1_before_update == zero_vector);
  BOOST_TEST(pos2_before_update == zero_vector);

  // Set positions of slice
  VariantMap params{{"name", std::string("pos")}, {"values", new_positions}};
  p_slice->call_method("set_param_parallel", params);

  // Check if positions are set as expected
  auto pos1_after_update =
      get_value<Utils::Vector3d>(p1_handle->get_parameter("pos"));
  auto pos2_after_update =
      get_value<Utils::Vector3d>(p2_handle->get_parameter("pos"));
  BOOST_TEST(pos1_after_update == get_value<Utils::Vector3d>(new_pos1));
  BOOST_TEST(pos2_after_update == get_value<Utils::Vector3d>(new_pos2));
}

BOOST_FIXTURE_TEST_CASE(set_particle_types, ParticleFactory) {
  // New positions for slice
  Variant new_type1 = 1;
  Variant new_type2 = 2;
  std::vector<Variant> new_types = {new_type1, new_type2};

  Utils::Factory<ObjectHandle> f;
  f.register_new<Particles::ParticleSlice>("ParticleSlice");
  f.register_new<Particles::ParticleHandle>("ParticleHandle");
  f.register_new<CellSystem::CellSystem>("CellSystem");

  boost::mpi::communicator comm;
  auto ctx = std::make_shared<LocalContext>(f, comm);

  // Create TestCellSystem object
  auto &&sp_cell_system{ctx->make_shared("CellSystem", {})};
  auto cell_system =
      std::dynamic_pointer_cast<CellSystem::CellSystem>(sp_cell_system);
  cell_system->bind_system(espresso::system);

  // Parameters for ScriptInterface objects
  VariantMap p1_init_params{{"id", 0}, {"__cell_structure", cell_system}};
  VariantMap p2_init_params{{"id", 1}, {"__cell_structure", cell_system}};
  VariantMap p_slice_init_params{{"id_selection", std::vector<int>{0, 1}},
                                 {"__cell_structure", cell_system}};

  // Create particle core objects
  create_particle(Utils::Vector3d{0., 0., 0.}, 0, 0);
  create_particle(Utils::Vector3d{0., 0., 0.}, 1, 0);

  // Create 2 ParticleHandle, and 1 ParticleSlice object
  auto &&sp_p1_handle{ctx->make_shared("ParticleHandle", {p1_init_params})};
  auto &&sp_p2_handle{ctx->make_shared("ParticleHandle", {p2_init_params})};
  // Particle handle 1
  auto p1_handle =
      std::dynamic_pointer_cast<Particles::ParticleHandle>(sp_p1_handle);
  // Particle handle 2
  auto p2_handle =
      std::dynamic_pointer_cast<Particles::ParticleHandle>(sp_p2_handle);

  // Particle slice
  auto &&sp_p_slice{ctx->make_shared("ParticleSlice", {p_slice_init_params})};
  std::shared_ptr<Particles::ParticleSlice> p_slice =
      std::dynamic_pointer_cast<Particles::ParticleSlice>(sp_p_slice);

  // Check if old types are actually zero
  auto type1_before_update = get_value<int>(p1_handle->get_parameter("type"));
  auto type2_before_update = get_value<int>(p2_handle->get_parameter("type"));
  BOOST_TEST(type1_before_update == 0);
  BOOST_TEST(type2_before_update == 0);

  // Set types of slice
  VariantMap params{{"name", std::string("type")}, {"values", new_types}};
  p_slice->call_method("set_param_parallel", params);

  // Check if types are set as expected
  auto type1_after_update = get_value<int>(p1_handle->get_parameter("type"));
  auto type2_after_update = get_value<int>(p2_handle->get_parameter("type"));
  BOOST_TEST(type1_after_update == get_value<int>(new_type1));
  BOOST_TEST(type2_after_update == get_value<int>(new_type2));
}

BOOST_FIXTURE_TEST_CASE(particle_modifier, ParticleFactory) {
  using Particles::ParticleModifier;
  Utils::Factory<ObjectHandle> f;
  f.register_new<ParticleModifier>("ParticleModifier");

  boost::mpi::communicator comm;
  auto ctx = std::make_shared<LocalContext>(f, comm);

  {
    auto &&sp{ctx->make_shared("ParticleModifier", {{"id", -2}})};
    auto so = std::dynamic_pointer_cast<ParticleModifier>(sp);
    auto const pid = get_value<int>(so->get_parameter("id"));
    BOOST_CHECK_EQUAL(pid, -2);
    BOOST_CHECK(is_none(so->do_call_method("unknown", {})));
  }

  {
    auto &&sp{ctx->make_shared("ParticleModifier", {})};
    auto so = std::dynamic_pointer_cast<ParticleModifier>(sp);
    auto const pid = get_value<int>(so->get_parameter("id"));
    BOOST_CHECK_EQUAL(pid, -1);
    BOOST_CHECK_EXCEPTION(
        (so->do_call_method("remove_particle", {})), std::logic_error,
        [](auto const &err) {
          return std::string{err.what()} ==
                 "ParticleModifier has no method 'remove_particle'";
        });
  }

  {
    VariantMap const invalid_params = {{"pos", 0}};
    BOOST_CHECK_EXCEPTION(
        (ctx->make_shared("ParticleModifier", invalid_params)),
        std::logic_error, [](auto const &err) {
          return std::string{err.what()} ==
                 "ParticleModifier is not meant to create new particles";
        });
  }
}

BOOST_FIXTURE_TEST_CASE(test_exceptions, ParticleFactory) {
  Variant non_vector_variant{1};

  Utils::Factory<ObjectHandle> f;
  f.register_new<Particles::ParticleSlice>("ParticleSlice");

  boost::mpi::communicator comm;
  auto ctx = std::make_shared<LocalContext>(f, comm);

  VariantMap p_slice_init_params{
      {"id_selection", std::vector<int>{0, 1}},
  };

  // Create particle core objects
  create_particle(Utils::Vector3d{0., 0., 0.}, 0, 0);
  create_particle(Utils::Vector3d{1., 0., 0.}, 1, 0);

  // Create particle slice
  auto &&sp_p_slice{ctx->make_shared("ParticleSlice", {p_slice_init_params})};
  auto p_slice =
      std::dynamic_pointer_cast<Particles::ParticleSlice>(sp_p_slice);

  VariantMap params{{"name", std::string("v")}, {"values", non_vector_variant}};

  BOOST_CHECK_EXCEPTION(p_slice->call_method("set_param_parallel", params),
                        Exception, [](auto const &err) {
                          return std::string{err.what()} ==
                                 "Values must be of type vector, got int";
                        });
}

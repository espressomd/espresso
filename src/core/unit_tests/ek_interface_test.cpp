/*
 * Copyright (C) 2023 The ESPResSo project
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

#define BOOST_TEST_MODULE EK interface test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include "ParticleFactory.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "config/config.hpp"
#include "ek/EKReactions.hpp"
#include "ek/EKWalberla.hpp"
#include "ek/Implementation.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "system/System.hpp"

#ifdef WALBERLA
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKContainer.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/electrokinetics/ek_poisson_none_init.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactant.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>
#endif

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <memory>
#include <stdexcept>
#include <string>

static struct {
  double kT = 1.3E-4;
  double density = 1.4;
  double diffusion = 3e-3;
  double valency = 1.;
  bool advection = true;
  bool friction_coupling = true;
  double tau = 0.01;
  double time_step = 0.01;
  double agrid = 1.;
  double skin = 0.51;
  Utils::Vector3d ext_efield = Utils::Vector3d{{0.01, 0.02, 0.03}};
  Utils::Vector3d box_dimensions = Utils::Vector3d::broadcast(8.);
  Utils::Vector3i grid_dimensions = Utils::Vector3i::broadcast(8);
} params;

namespace espresso {
// ESPResSo system instance
static std::unique_ptr<EspressoSystemStandAlone> system;
// ESPResSo actors
#ifdef WALBERLA
static std::shared_ptr<EK::EKWalberla::ek_container_type> ek_container;
static std::shared_ptr<EK::EKWalberla::ek_reactions_type> ek_reactions;
static std::shared_ptr<EK::EKWalberla> ek_instance;
static std::shared_ptr<LatticeWalberla> ek_lattice;
#endif

static auto make_ek_actor() {
#ifdef WALBERLA
  auto constexpr n_ghost_layers = 1u;
  auto constexpr single_precision = true;
  ek_lattice = std::make_shared<LatticeWalberla>(params.grid_dimensions,
                                                 ::node_grid, n_ghost_layers);
  ek_container = std::make_shared<EK::EKWalberla::ek_container_type>(
      params.tau, new_ek_poisson_none(ek_lattice, single_precision));
  ek_reactions = std::make_shared<EK::EKWalberla::ek_reactions_type>();
  ek_instance = std::make_shared<EK::EKWalberla>(ek_container, ek_reactions);
#endif
}

static void add_ek_actor() {
#ifdef WALBERLA
  System::get_system().ek.set<::EK::EKWalberla>(ek_instance);
#endif
}

static void remove_ek_actor() { System::get_system().ek.reset(); }
} // namespace espresso

#ifdef WALBERLA
namespace walberla {
class EKReactionImpl : public EKReactionBase {
public:
  using EKReactionBase::EKReactionBase;
  using EKReactionBase::get_coefficient;
  using EKReactionBase::get_lattice;
  using EKReactionBase::get_reactants;

  void perform_reaction() override {}
  ~EKReactionImpl() override = default;
};
} // namespace walberla
#endif // WALBERLA

/** Fixture to manage the lifetime of the EK actor. */
struct CleanupActorEK : public ParticleFactory {
  CleanupActorEK() : ParticleFactory() {
    espresso::make_ek_actor();
    espresso::add_ek_actor();
  }

  // NOLINTNEXTLINE(bugprone-exception-escape)
  ~CleanupActorEK() { espresso::remove_ek_actor(); }
};

BOOST_FIXTURE_TEST_SUITE(suite, CleanupActorEK)

static auto get_n_runtime_errors() { return check_runtime_errors_local(); }

#ifdef WALBERLA
BOOST_AUTO_TEST_CASE(ek_interface_walberla) {
  auto &ek = System::get_system().ek;

  {
    // tau setters and getters
    espresso::ek_container->set_tau(2.);
    BOOST_CHECK_EQUAL(ek.get_tau(), 2.);
    BOOST_CHECK_EQUAL(ek.get_steps_per_md_step(1.), 2);
    BOOST_CHECK_EQUAL(ek.get_steps_per_md_step(2.), 1);
    BOOST_CHECK_EQUAL(ek.get_steps_per_md_step(5.), 0);
  }

  {
    // setup a minimal EK model without coupling to LB
    using walberla::EKReactant;
    auto constexpr single_precision = true;
    auto constexpr stoich = 1.;
    auto constexpr order = 2.;
    auto ek_species = new_ek_walberla(
        espresso::ek_lattice, params.diffusion, params.kT, params.valency,
        params.ext_efield, params.density, false, false, single_precision);
    auto ek_reactant = std::make_shared<EKReactant>(ek_species, stoich, order);
    auto ek_reaction = std::make_shared<walberla::EKReactionImpl>(
        espresso::ek_lattice,
        std::vector<std::shared_ptr<EKReactant>>{
            {ek_reactant, ek_reactant, ek_reactant}},
        1.);
    ek_reaction->perform_reaction();

    BOOST_REQUIRE(espresso::ek_reactions->empty());
    BOOST_REQUIRE(espresso::ek_container->empty());
    ek.propagate(); // no-op
    BOOST_REQUIRE_EQUAL(get_n_runtime_errors(), 0);
    BOOST_REQUIRE(espresso::ek_container->empty());
    BOOST_REQUIRE(not espresso::ek_container->contains(ek_species));
    BOOST_REQUIRE(not espresso::ek_reactions->contains(ek_reaction));
    espresso::ek_container->set_tau(2.);
    espresso::ek_container->add(ek_species);
    BOOST_REQUIRE(not espresso::ek_container->empty());
    BOOST_REQUIRE(espresso::ek_container->contains(ek_species));
    ek.propagate(); // no-op
    BOOST_REQUIRE_EQUAL(get_n_runtime_errors(), 0);
    espresso::ek_reactions->add(ek_reaction);
    BOOST_REQUIRE(espresso::ek_reactions->contains(ek_reaction));
    BOOST_REQUIRE(not espresso::ek_reactions->empty());
    espresso::ek_reactions->remove(ek_reaction);
    BOOST_REQUIRE(not espresso::ek_reactions->contains(ek_reaction));
    BOOST_REQUIRE(espresso::ek_reactions->empty());
    espresso::ek_container->remove(ek_species);
    BOOST_REQUIRE(espresso::ek_container->empty());
    BOOST_REQUIRE(not espresso::ek_container->contains(ek_species));
    ek.propagate(); // no-op
    BOOST_REQUIRE_EQUAL(get_n_runtime_errors(), 0);
  }

  {
    // EK prevents changing most of the system state
    BOOST_CHECK_THROW(ek.on_boxl_change(), std::runtime_error);
    BOOST_CHECK_THROW(ek.on_timestep_change(), std::runtime_error);
    BOOST_CHECK_THROW(ek.on_temperature_change(), std::runtime_error);
    BOOST_CHECK_THROW(ek.on_node_grid_change(), std::runtime_error);
  }
}
#endif // WALBERLA

BOOST_AUTO_TEST_CASE(ek_interface_none) {
  auto &ek = System::get_system().ek;

  {
    using EK::NoEKActive;
    ek.reset();
    BOOST_CHECK_THROW(ek.get_tau(), NoEKActive);
    BOOST_CHECK_THROW(ek.propagate(), NoEKActive);
    auto ek_impl = std::make_shared<EK::EKNone>();
    ek.set<EK::EKNone>(ek_impl);
    BOOST_CHECK_THROW(ek.is_ready_for_propagation(), NoEKActive);
    BOOST_CHECK_THROW(ek.propagate(), NoEKActive);
    BOOST_CHECK_THROW(ek.get_tau(), NoEKActive);
    BOOST_CHECK_THROW(ek.sanity_checks(), NoEKActive);
    BOOST_CHECK_THROW(ek.veto_time_step(0.), NoEKActive);
    BOOST_CHECK_THROW(ek.on_cell_structure_change(), NoEKActive);
    BOOST_CHECK_THROW(ek.on_boxl_change(), NoEKActive);
    BOOST_CHECK_THROW(ek.on_node_grid_change(), NoEKActive);
    BOOST_CHECK_THROW(ek.on_timestep_change(), NoEKActive);
    BOOST_CHECK_THROW(ek.on_temperature_change(), NoEKActive);
    BOOST_CHECK_THROW(ek.get_steps_per_md_step(1.), std::exception);
    auto const err_msg = std::string(NoEKActive().what());
    auto const ref_msg = std::string("EK not activated");
    BOOST_CHECK_EQUAL(err_msg, ref_msg);
    ek.reset();
  }
}

BOOST_AUTO_TEST_SUITE_END()

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  espresso::system->set_box_l(params.box_dimensions);
  espresso::system->set_time_step(params.time_step);
  espresso::system->set_skin(params.skin);

  boost::mpi::communicator world;
  assert(world.size() <= 2);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

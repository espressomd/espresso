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

#define BOOST_TEST_MODULE LB particle coupling test
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include "EspressoSystemStandAlone.hpp"
#include "config/config.hpp"
#include "errorhandling.hpp"
#include "grid_based_algorithms/ek_container.hpp"
#include "grid_based_algorithms/ek_reactions.hpp"

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
} // namespace espresso

static auto get_n_runtime_errors() { return check_runtime_errors_local(); }

#ifdef WALBERLA

#include "grid.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKContainer.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/electrokinetics/ek_poisson_none_init.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>

BOOST_AUTO_TEST_CASE(ek_interface_walberla) {
  {
    // tau setters and getters
    BOOST_CHECK_EQUAL(EK::ek_container.get_tau(), 0.);
    BOOST_CHECK_EQUAL(EK::get_tau(), 0.);
    BOOST_CHECK_EQUAL(EK::get_steps_per_md_step(1.), 0);
    EK::ek_container.set_tau(2.);
    BOOST_CHECK_EQUAL(EK::ek_container.get_tau(), 2.);
    BOOST_CHECK_EQUAL(EK::get_tau(), 2.);
    BOOST_CHECK_EQUAL(EK::get_steps_per_md_step(1.), 2);
    BOOST_CHECK_EQUAL(EK::get_steps_per_md_step(2.), 1);
    BOOST_CHECK_EQUAL(EK::get_steps_per_md_step(5.), 0);
  }

  {
    // setup a minimal EK model without coupling to LB
    auto constexpr n_ghost_layers = 1u;
    auto constexpr single_precision = true;
    auto ek_lattice = std::make_shared<LatticeWalberla>(
        params.grid_dimensions, ::node_grid, n_ghost_layers);
    auto ek_species = new_ek_walberla(
        ek_lattice, params.diffusion, params.kT, params.valency,
        params.ext_efield, params.density, false, false, single_precision);
    auto ek_solver_none = new_ek_poisson_none(ek_lattice, single_precision);

    BOOST_REQUIRE(EK::ek_reactions.empty());
    BOOST_REQUIRE(EK::ek_container.empty());
    BOOST_REQUIRE(not EK::ek_container.is_poisson_solver_set());
    EK::propagate(); // no-op
    BOOST_REQUIRE_EQUAL(get_n_runtime_errors(), 0);
    EK::ek_container.set_poisson_solver(ek_solver_none);
    BOOST_REQUIRE(EK::ek_container.is_poisson_solver_set());
    BOOST_REQUIRE(EK::ek_container.empty());
    EK::ek_container.set_tau(0.);
    BOOST_CHECK_THROW(EK::ek_container.add(ek_species), std::runtime_error);
    EK::ek_container.set_tau(2.);
    EK::ek_container.add(ek_species);
    BOOST_REQUIRE(not EK::ek_container.empty());
    EK::propagate(); // no-op
    BOOST_REQUIRE_EQUAL(get_n_runtime_errors(), 0);
    EK::ek_container.remove(ek_species);
    BOOST_REQUIRE(EK::ek_container.empty());
    EK::propagate(); // no-op
    BOOST_REQUIRE_EQUAL(get_n_runtime_errors(), 0);
  }
}

#else // WALBERLA

BOOST_AUTO_TEST_CASE(ek_interface) {
  {
    EK::propagate(); // no-op
    BOOST_CHECK_THROW(EK::get_tau(), NoEKActive);
    BOOST_CHECK_THROW(EK::get_tau(), std::exception);
    BOOST_CHECK_THROW(EK::get_steps_per_md_step(1.), std::exception);
    auto const err_msg = std::string(NoEKActive().what());
    auto const ref_msg = std::string("EK not activated");
    BOOST_CHECK_EQUAL(err_msg, ref_msg);
  }
}

#endif // WALBERLA

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  espresso::system->set_box_l(params.box_dimensions);
  espresso::system->set_time_step(params.time_step);
  espresso::system->set_skin(params.skin);

  boost::mpi::communicator world;
  assert(world.size() <= 2);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

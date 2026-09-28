/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE "Verlet list update test"

#include <config/config.hpp>

#define BOOST_TEST_DYN_LINK
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;
namespace bdata = boost::unit_test::data;

#include "ParticleFactory.hpp"
#include "particle_management.hpp"

#include "EspressoCoreGlobalConfig.hpp"
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "integrate.hpp"
#include "integrators/Propagation.hpp"
#include "integrators/steepest_descent.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "npt.hpp"
#include "particle_node.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <vector>

namespace espresso {
// ESPResSo system instance
static std::shared_ptr<System::System> system;
} // namespace espresso

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    espresso::system = System::System::create();
    espresso::system->set_cell_structure_topology(CellStructureType::REGULAR);
    ::System::set_system(espresso::system);
  }
  ~GlobalConfig() {
    espresso::system.reset();
    ::System::reset_system();
  }
};

// Decorator to skip tests when the number of MPI ranks isn't 4
boost::test_tools::assertion_result has_4_mpi_ranks(utf::test_unit_id) {
  boost::mpi::communicator world;
  return world.size() == 4;
}

// Decorator to skip tests if Lennard-Jones isn't compiled in
boost::test_tools::assertion_result has_lj(utf::test_unit_id) {
#ifdef ESPRESSO_LENNARD_JONES
  return true;
#else
  return false;
#endif
}

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

namespace Testing {
/**
 * Helper class to setup an integrator and particle properties such that the
 * particle is on a collision course with another particle on the x-axis.
 */
struct IntegratorHelper : public ParticleFactory {
  IntegratorHelper() = default;
  virtual ~IntegratorHelper() = default;
  /** Set integrator parameters. */
  virtual void set_integrator() const = 0;
  /** Set particle to move along the x-axis. */
  virtual void set_particle_properties(int) const = 0;
  virtual char const *name() const = 0;
  friend auto operator<<(std::ostream &os, IntegratorHelper const &obj)
      -> std::ostream & {
    return os << obj.name();
  }
};

#ifdef ESPRESSO_EXTERNAL_FORCES
struct : public IntegratorHelper {
  void set_integrator() const override {
    espresso::system->thermostat->thermo_switch = THERMO_OFF;
    espresso::system->steepest_descent =
        std::make_shared<SteepestDescent>(0., 0.01, 100.);
    espresso::system->propagation->set_integ_switch(
        INTEG_METHOD_STEEPEST_DESCENT);
  }
  void set_particle_properties(int pid) const override {
    set_particle_property(pid, &Particle::ext_force,
                          Utils::Vector3d{20., 0., 0.});
  }
  char const *name() const override { return "SteepestDescent"; }
} steepest_descent;
#endif // ESPRESSO_EXTERNAL_FORCES

struct : public IntegratorHelper {
  void set_integrator() const override {
    espresso::system->thermostat->thermo_switch = THERMO_OFF;
    espresso::system->steepest_descent = std::make_shared<SteepestDescent>();
    espresso::system->propagation->set_integ_switch(INTEG_METHOD_NVT);
  }
  void set_particle_properties(int pid) const override {
    set_particle_v(pid, {20., 0., 0.});
  }
  char const *name() const override { return "VelocityVerlet"; }
} velocity_verlet;

#ifdef ESPRESSO_NPT
struct : public IntegratorHelper {
  void set_integrator() const override {
    auto &npt_iso = espresso::system->thermostat->npt_iso;
    espresso::system->nptiso = std::make_shared<NptIsoParameters>(
        1., 1e16, Utils::Vector<bool, 3>{true, true, true}, true);
    espresso::system->propagation->set_integ_switch(INTEG_METHOD_NPT_ISO_AND);
    espresso::system->thermostat->thermo_switch = THERMO_NPT_ISO;
    espresso::system->thermostat->kT = 1.;
    npt_iso = std::make_shared<IsotropicNptThermostat>();
    npt_iso->rng_initialize(0u);
    npt_iso->gammav = 0.; // disable box volume change
    npt_iso->gamma0 = 0.;
  }
  void set_particle_properties(int pid) const override {
    set_particle_v(pid, {20., 0., 0.});
  }
  char const *name() const override { return "VelocityVerletNpT"; }
} velocity_verlet_npt;
#endif // ESPRESSO_NPT

} // namespace Testing

inline double get_dist_from_last_verlet_update(Particle const &p) {
  return (Utils::Vector3d(p.pos()) -
          Utils::Vector3d(p.pos_at_last_verlet_update()))
      .norm();
}

inline double get_dist_from_pair(Particle const &p1, Particle const &p2) {
  return (Utils::Vector3d(p1.pos()) - Utils::Vector3d(p2.pos())).norm();
}

auto const node_grids = std::vector<Utils::Vector3i>{{4, 1, 1}, {2, 2, 1}};
auto const propagators =
    std::vector<std::reference_wrapper<Testing::IntegratorHelper>>{
        Testing::velocity_verlet,
#ifdef ESPRESSO_NPT
        Testing::velocity_verlet_npt,
#endif // ESPRESSO_NPT
#ifdef ESPRESSO_EXTERNAL_FORCES
        Testing::steepest_descent
#endif // ESPRESSO_EXTERNAL_FORCES
    };

BOOST_TEST_DECORATOR(*utf::precondition(has_4_mpi_ranks) *
                     utf::precondition(has_lj))
BOOST_DATA_TEST_CASE_F(ParticleFactory, verlet_list_update,
                       bdata::make(node_grids) * bdata::make(propagators),
                       node_grid, integration_helper) {
#ifdef ESPRESSO_LENNARD_JONES
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();
  auto const comm = boost::mpi::communicator();
  auto const rank = comm.rank();

  auto const box_l = 8.;
  auto &system = *espresso::system;
  system.set_box_l(Utils::Vector3d::broadcast(box_l));
  ::communicator.set_node_grid(node_grid);
  system.on_node_grid_change();

  // particle properties
  auto const pid1 = 9;
  auto const pid2 = 2;
  auto const pid3 = 5;

  // distance between particles
  auto const dist = 0.2;
  // set up LJ potential
  auto const eps = 1.0;
  auto const sig = 0.05;
  auto const shift = 0.0;
  auto const offset = 0.1;
  auto const min = 0.0;
  auto const r_off = dist - offset;
  auto const cut = r_off + 1e-3;
  LJ_Parameters lj{eps, sig, cut, offset, min, shift};
  system.nonbonded_ias->make_particle_type_exist(1);
  system.nonbonded_ias->get_ia_param(0, 1).lj = lj;
  system.on_non_bonded_ia_change();

  // set up velocity-Verlet integrator
  auto const time_step = 0.01;
  auto const skin = 0.1;
  system.set_time_step(time_step);
  system.cell_structure->set_verlet_skin(skin);
  integration_helper.get().set_integrator();

  // If the Verlet list is not updated, two particles initially placed in
  // different cells will never see each other, even when closer than the
  // interaction cutoff. To reproduce that bug, we check the case where
  // we have 4 cells in one direction and particles in cells #0 and #2
  // (i.e. the non-neighboring cell case) and the case where we have
  // 2 cells in one direction (i.e. the neighboring cell case).
  create_particle({2. - 0.10, 1., 1.}, pid1, 0);
  create_particle({4. + 0.15, 1., 1.}, pid2, 1);
  create_particle({2. - 0.10, 5., 1.}, pid3, 0);
  BOOST_REQUIRE_EQUAL(get_particle_node_parallel(pid1), rank ? -1 : 0);
  BOOST_REQUIRE_EQUAL(get_particle_node_parallel(pid2), rank ? -1 : 2);
  BOOST_REQUIRE_EQUAL(get_particle_node_parallel(pid3),
                      rank ? -1 : ((node_grid[0] == 4) ? 0 : 1));

  // check that particles in different cells will eventually interact during
  // normal integration (the integration helper sets a collision course)
  {
    integration_helper.get().set_particle_properties(pid1);

    // integrate until both particles are closer than cutoff
    {
      system.integrate(11, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const p1_opt = copy_particle_to_head_node(comm, system, pid1);
      auto const p2_opt = copy_particle_to_head_node(comm, system, pid2);
      if (rank == 0) {
        auto const &p1 = *p1_opt;
        auto const &p2 = *p2_opt;
        BOOST_REQUIRE_LT(get_dist_from_pair(p1, p2), cut);
      }
    }

    // check forces and Verlet update
    {
      system.integrate(1, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const p1_opt = copy_particle_to_head_node(comm, system, pid1);
#ifdef ESPRESSO_EXTERNAL_FORCES
      auto const p2_opt = copy_particle_to_head_node(comm, system, pid2);
#endif // ESPRESSO_EXTERNAL_FORCES
      if (rank == 0) {
        auto const &p1 = *p1_opt;
#ifdef ESPRESSO_EXTERNAL_FORCES
        auto const &p2 = *p2_opt;
        BOOST_CHECK_CLOSE(p1.force()[0] - p1.ext_force()[0], 480., 1e-9);
#endif // ESPRESSO_EXTERNAL_FORCES
        BOOST_CHECK_CLOSE(p1.force()[1], 0., tol);
        BOOST_CHECK_CLOSE(p1.force()[2], 0., tol);
#ifdef ESPRESSO_EXTERNAL_FORCES
        BOOST_TEST(Utils::Vector3d(p1.force()) - p1.ext_force() ==
                       -Utils::Vector3d(p2.force()),
                   boost::test_tools::per_element());
#endif // ESPRESSO_EXTERNAL_FORCES
        BOOST_CHECK_LT(get_dist_from_last_verlet_update(p1), skin / 2.);
      }
    }
  }

  // check that particles in different cells will interact when manually
  // placed next to each other (@c set_particle_pos() resorts particles)
  {
    ::set_particle_pos(pid3, {4. + 0.10, 1., 1.0});
    {
      auto const p2_opt = copy_particle_to_head_node(comm, system, pid2);
      auto const p3_opt = copy_particle_to_head_node(comm, system, pid3);
      if (rank == 0) {
        auto const &p2 = *p2_opt;
        auto const &p3 = *p3_opt;
        BOOST_REQUIRE_LT(get_dist_from_pair(p2, p3), cut);
        BOOST_CHECK_GT(get_dist_from_last_verlet_update(p3), skin / 2.);
      }
    }
    {
      system.integrate(0, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const p3_opt = copy_particle_to_head_node(comm, system, pid3);
      if (rank == 0) {
        auto const &p3 = *p3_opt;
        BOOST_CHECK_LT(get_dist_from_last_verlet_update(p3), skin / 2.);
      }
    }
  }
#endif // ESPRESSO_LENNARD_JONES
}

BOOST_AUTO_TEST_SUITE_END()

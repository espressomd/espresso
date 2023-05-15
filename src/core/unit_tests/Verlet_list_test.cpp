/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE Verlet list update test

#include "config/config.hpp"

#ifdef LENNARD_JONES

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;
namespace bdata = boost::unit_test::data;

#include "ParticleFactory.hpp"
#include "particle_management.hpp"

#include "EspressoSystemStandAlone.hpp"
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "integrate.hpp"
#include "integrators/steepest_descent.hpp"
#include "nonbonded_interactions/lj.hpp"
#include "npt.hpp"
#include "particle_node.hpp"
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
static std::unique_ptr<EspressoSystemStandAlone> system;
} // namespace espresso

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

#ifdef EXTERNAL_FORCES
struct : public IntegratorHelper {
  void set_integrator() const override {
    mpi_set_thermo_switch_local(THERMO_OFF);
    register_integrator(SteepestDescentParameters(0., 0.01, 100.));
    set_integ_switch(INTEG_METHOD_STEEPEST_DESCENT);
  }
  void set_particle_properties(int pid) const override {
    set_particle_property(pid, &Particle::ext_force,
                          Utils::Vector3d{20., 0., 0.});
  }
  char const *name() const override { return "SteepestDescent"; }
} steepest_descent;
#endif // EXTERNAL_FORCES

struct : public IntegratorHelper {
  void set_integrator() const override {
    mpi_set_thermo_switch_local(THERMO_OFF);
    set_integ_switch(INTEG_METHOD_NVT);
  }
  void set_particle_properties(int pid) const override {
    set_particle_v(pid, {20., 0., 0.});
  }
  char const *name() const override { return "VelocityVerlet"; }
} velocity_verlet;

#ifdef NPT
struct : public IntegratorHelper {
  void set_integrator() const override {
    ::nptiso = NptIsoParameters(1., 1e9, {true, true, true}, true);
    set_integ_switch(INTEG_METHOD_NPT_ISO);
    mpi_set_temperature_local(1.);
    mpi_npt_iso_set_rng_seed(0);
    mpi_set_thermo_switch_local(thermo_switch | THERMO_NPT_ISO);
    mpi_set_nptiso_gammas_local(0., 0.); // disable box volume change
  }
  void set_particle_properties(int pid) const override {
    set_particle_v(pid, {20., 0., 0.});
  }
  char const *name() const override { return "VelocityVerletNpT"; }
} velocity_verlet_npt;
#endif // NPT

} // namespace Testing

inline double get_dist_from_last_verlet_update(Particle const &p) {
  return (p.pos() - p.pos_at_last_verlet_update()).norm();
}

inline double get_dist_from_pair(Particle const &p1, Particle const &p2) {
  return (p1.pos() - p2.pos()).norm();
}

auto const node_grids = std::vector<Utils::Vector3i>{{4, 1, 1}, {2, 2, 1}};
auto const propagators =
    std::vector<std::reference_wrapper<Testing::IntegratorHelper>>{
        Testing::velocity_verlet,
#ifdef NPT
        Testing::velocity_verlet_npt,
#endif // NPT
#ifdef EXTERNAL_FORCES
        Testing::steepest_descent
#endif // EXTERNAL_FORCES
    };

BOOST_DATA_TEST_CASE_F(ParticleFactory, verlet_list_update,
                       bdata::make(node_grids) * bdata::make(propagators),
                       node_grid, integration_helper) {
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();
  auto const comm = boost::mpi::communicator();
  auto const rank = comm.rank();

  auto const box_l = 8.;
  espresso::system->set_box_l(Utils::Vector3d::broadcast(box_l));
  espresso::system->set_node_grid(node_grid);

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
  make_particle_type_exist(1);
  auto const key = get_ia_param_key(0, 1);
  LJ_Parameters lj{eps, sig, cut, offset, min, shift};
  ::nonbonded_ia_params[key]->lj = lj;
  on_non_bonded_ia_change();

  // set up velocity-Verlet integrator
  auto const time_step = 0.01;
  auto const skin = 0.1;
  espresso::system->set_time_step(time_step);
  espresso::system->set_skin(skin);
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
      integrate(11, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const p1_opt = copy_particle_to_head_node(comm, pid1);
      auto const p2_opt = copy_particle_to_head_node(comm, pid2);
      if (rank == 0) {
        auto const &p1 = *p1_opt;
        auto const &p2 = *p2_opt;
        BOOST_REQUIRE_LT(get_dist_from_pair(p1, p2), cut);
      }
    }

    // check forces and Verlet update
    {
      integrate(1, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const p1_opt = copy_particle_to_head_node(comm, pid1);
#ifdef EXTERNAL_FORCES
      auto const p2_opt = copy_particle_to_head_node(comm, pid2);
#endif // EXTERNAL_FORCES
      if (rank == 0) {
        auto const &p1 = *p1_opt;
#ifdef EXTERNAL_FORCES
        auto const &p2 = *p2_opt;
        BOOST_CHECK_CLOSE(p1.force()[0] - p1.ext_force()[0], 480., 1e-9);
#endif // EXTERNAL_FORCES
        BOOST_CHECK_CLOSE(p1.force()[1], 0., tol);
        BOOST_CHECK_CLOSE(p1.force()[2], 0., tol);
#ifdef EXTERNAL_FORCES
        BOOST_TEST(p1.force() - p1.ext_force() == -p2.force(),
                   boost::test_tools::per_element());
#endif // EXTERNAL_FORCES
        BOOST_CHECK_LT(get_dist_from_last_verlet_update(p1), skin / 2.);
      }
    }
  }

  // check that particles in different cells will interact when manually
  // placed next to each other (@c set_particle_pos() resorts particles)
  {
    ::set_particle_pos(pid3, {4. + 0.10, 1., 1.0});
    {
      auto const p2_opt = copy_particle_to_head_node(comm, pid2);
      auto const p3_opt = copy_particle_to_head_node(comm, pid3);
      if (rank == 0) {
        auto const &p2 = *p2_opt;
        auto const &p3 = *p3_opt;
        BOOST_REQUIRE_LT(get_dist_from_pair(p2, p3), cut);
        BOOST_CHECK_GT(get_dist_from_last_verlet_update(p3), skin / 2.);
      }
    }
    {
      integrate(0, INTEG_REUSE_FORCES_CONDITIONALLY);
      auto const p3_opt = copy_particle_to_head_node(comm, pid3);
      if (rank == 0) {
        auto const &p3 = *p3_opt;
        BOOST_CHECK_LT(get_dist_from_last_verlet_update(p3), skin / 2.);
      }
    }
  }
}

int main(int argc, char **argv) {
  espresso::system = std::make_unique<EspressoSystemStandAlone>(argc, argv);
  // the test case only works for 4 MPI ranks
  boost::mpi::communicator world;
  int error_code = 0;
  if (world.size() == 4) {
    error_code = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  }
  return error_code;
}
#else // ifdef LENNARD_JONES
int main(int argc, char **argv) {}
#endif

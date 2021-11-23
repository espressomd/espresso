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

#define BOOST_TEST_MODULE LB exception mechanism
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"

#include <stdexcept>

BOOST_AUTO_TEST_CASE(exceptions) {
  // getters and setters
  BOOST_CHECK_THROW(lb_lbfluid_get_rng_state(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_rng_state(0u), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_density(-1.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_density(1.), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_density(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_viscosity(-1.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_viscosity(1.), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_viscosity(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_bulk_viscosity(-1.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_bulk_viscosity(1.), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_bulk_viscosity(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_gamma_odd(2.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_gamma_odd({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_gamma_odd(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_gamma_even(2.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_gamma_even({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_gamma_even(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_agrid(-1.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_agrid(1.), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_agrid(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_ext_force_density({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_ext_force_density(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_tau(-1.), std::invalid_argument);
  BOOST_CHECK_THROW(lb_lbfluid_set_tau(1.), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_tau(), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_set_kT({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_kT(), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_get_boundary({}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_set_density({}, {}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_get_density({}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_set_velocity({}, {}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_get_velocity({}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_set_pop({}, {}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_get_pop({}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_get_pressure_tensor({}), std::exception);
  BOOST_CHECK_THROW(lb_lbnode_get_pressure_tensor_neq({}), std::exception);
  // particle coupling and interpolation
  BOOST_CHECK_EQUAL(lb_lbcoupling_get_rng_state(), 0u);
  BOOST_CHECK_THROW(lb_lbfluid_get_interpolated_velocity({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_interpolated_density({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_shape(), std::exception);
  BOOST_CHECK_EQUAL(lb_lbfluid_calc_fluid_momentum(), Utils::Vector3d{});
  BOOST_CHECK_THROW(lb_lbfluid_set_lattice_switch(static_cast<ActiveLB>(100)),
                    std::invalid_argument);
  ::lattice_switch = ActiveLB::CPU;
  mpi_set_interpolation_order_local(InterpolationOrder::quadratic);
  BOOST_CHECK_THROW(lb_lbfluid_get_interpolated_density({}), std::exception);
  BOOST_CHECK_THROW(lb_lbfluid_get_interpolated_velocity({}),
                    std::runtime_error);
  BOOST_CHECK_THROW(lb_lbinterpolation_add_force_density({}, {}),
                    std::runtime_error);
  ::lattice_switch = ActiveLB::GPU;
  BOOST_CHECK_THROW(lb_lbfluid_get_interpolated_density({}), std::exception);
  ::lattice_switch = ActiveLB::NONE;
  mpi_set_interpolation_order_local(InterpolationOrder::linear);
}

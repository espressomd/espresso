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

#include "config.hpp"

#ifdef LB_WALBERLA

#include "FluidNodeWalberla.hpp"

#include "core/MpiCallbacks.hpp"
#include "core/communication.hpp"
#include "core/grid_based_algorithms/lb_interface.hpp"
#include "core/grid_based_algorithms/lb_walberla_instance.hpp"
#include "core/grid_based_algorithms/lb_walberla_interface.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <vector>

namespace Walberla {

static boost::optional<Utils::Vector3d>
mpi_get_node_velocity_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_velocity(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_velocity_local)

static boost::optional<Utils::Vector3d>
mpi_get_node_velocity_at_boundary_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_velocity_at_boundary(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_velocity_at_boundary_local)

static boost::optional<Utils::Vector3d>
mpi_get_node_last_applied_force_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_last_applied_force(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_last_applied_force_local)

static boost::optional<double> mpi_get_node_density_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_density(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_density_local)

static boost::optional<bool>
mpi_get_node_is_boundary_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_is_boundary(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_is_boundary_local)

static boost::optional<Utils::Vector3d>
mpi_get_node_boundary_force_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_boundary_force(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_boundary_force_local)

static boost::optional<std::vector<double>>
mpi_get_node_pop_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_pop(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_pop_local)

static boost::optional<Utils::Vector6d>
mpi_get_node_pressure_tensor_local(Utils::Vector3i ind) {
  return lb_walberla()->get_node_pressure_tensor(ind);
}

REGISTER_CALLBACK_ONE_RANK(mpi_get_node_pressure_tensor_local)

} // namespace Walberla

namespace ScriptInterface::walberla {
FluidNodeWalberla::FluidNodeWalberla() {
  add_parameters(
      {{"_index", AutoParameter::read_only, [this]() { return (m_index); }},
       {"velocity",
        [this](const Variant &v) {
          if (auto lb_fluid = m_lb_fluid.lock()) {
            auto const u = get_value<Utils::Vector3d>(v) * m_conv_velocity;
            lb_fluid->set_node_velocity(m_index, u);
            lb_fluid->ghost_communication();
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
        },
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const result = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_velocity_local, m_index);
            return Variant{result / m_conv_velocity};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"velocity_at_boundary",
        [this](const Variant &v) {
          if (auto lb_fluid = m_lb_fluid.lock()) {
            if (is_none(v)) {
              lb_fluid->remove_node_from_boundary(m_index, true);
              lb_fluid->ghost_communication();
            } else {
              auto const u = get_value<Utils::Vector3d>(v) * m_conv_velocity;
              lb_fluid->set_node_velocity_at_boundary(m_index, u, true);
              lb_fluid->ghost_communication();
            }
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
        },
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const is_boundary = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_is_boundary_local, m_index);
            if (is_boundary) {
              auto const u = ::Communication::mpiCallbacks().call(
                  ::Communication::Result::one_rank,
                  Walberla::mpi_get_node_velocity_at_boundary_local, m_index);
              return Variant{u * (1. / m_conv_velocity)};
            }
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"density",
        [this](const Variant &v) {
          if (auto lb_fluid = m_lb_fluid.lock()) {
            auto const density = get_value<double>(v) * m_conv_dens;
            lb_fluid->set_node_density(m_index, density);
            lb_fluid->ghost_communication();
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
        },
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const result = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_density_local, m_index);
            return Variant{result / m_conv_dens};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"_population",
        [this](const Variant &v) {
          if (auto lb_fluid = m_lb_fluid.lock()) {
            auto const population = get_value<std::vector<double>>(v);
            lb_fluid->set_node_pop(m_index, population);
            lb_fluid->ghost_communication();
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
        },
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const result = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_pop_local, m_index);
            return Variant{result};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"is_boundary", AutoParameter::read_only,
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const result = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_is_boundary_local, m_index);
            return Variant{result};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"boundary_force", AutoParameter::read_only,
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const result = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_boundary_force_local, m_index);
            return Variant{result * (1. / m_conv_force)};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"_pressure_tensor", AutoParameter::read_only,
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            auto lb_fluid = m_lb_fluid.lock();
            assert(lb_fluid.get() == ::lb_walberla().get());
            auto const visc = lb_fluid->get_viscosity();
            auto tensor = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_pressure_tensor_local, m_index);
            Walberla::walberla_off_diagonal_correction(tensor, visc);
            return Variant{tensor * (1. / m_conv_press)};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }},
       {"last_applied_force",
        [this](const Variant &v) {
          if (auto lb_fluid = m_lb_fluid.lock()) {
            auto const f = get_value<Utils::Vector3d>(v) * m_conv_force;
            lb_fluid->set_node_last_applied_force(m_index, f);
            lb_fluid->ghost_communication();
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
        },
        [this]() {
          if (lattice_switch == ActiveLB::WALBERLA) {
            assert(m_lb_fluid.lock() and
                   m_lb_fluid.lock().get() == ::lb_walberla().get());
            auto const result = ::Communication::mpiCallbacks().call(
                ::Communication::Result::one_rank,
                Walberla::mpi_get_node_last_applied_force_local, m_index);
            return Variant{result * (1. / m_conv_force)};
          } else if (context()->is_head_node()) {
            throw NoLBActive();
          }
          return Variant{None{}};
        }}});
}

} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA

/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "lb.hpp"
#include "lb_constants.hpp"
#include "lb_interpolation.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <boost/optional.hpp>

using Utils::get_linear_index;

/* LB CPU callback interface */
namespace detail {

template <typename Kernel>
void lb_set(Utils::Vector3i const &index, Kernel kernel) {
  if (lblattice.is_local(index)) {
    kernel(index);
  }
}

template <typename Kernel>
auto lb_calc(Utils::Vector3i const &index, Kernel kernel) {
  using R = decltype(kernel(index));
  if (lblattice.is_local(index)) {
    return boost::optional<R>(kernel(index));
  }
  return boost::optional<R>();
}

template <typename Kernel>
auto lb_calc_for_pos(Utils::Vector3d const &pos, Kernel kernel) {
  using R = decltype(kernel(pos));
  if (map_position_node_array(pos) == this_node) {
    return boost::optional<R>(kernel(pos));
  }
  return boost::optional<R>();
}

template <class Kernel>
auto lb_calc_fluid_kernel(Utils::Vector3i const &index, Kernel kernel) {
  return lb_calc(index, [&](auto index) {
    auto const linear_index =
        get_linear_index(lblattice.local_index(index), lblattice.halo_grid);
    auto const force_density = lbfields[linear_index].force_density;
    auto const modes = lb_calc_modes(linear_index, lbfluid);
    return kernel(modes, force_density);
  });
}
} // namespace detail

boost::optional<Utils::Vector3d>
mpi_lb_get_interpolated_velocity(Utils::Vector3d const &pos) {
  return detail::lb_calc_for_pos(pos, [&](auto pos) {
    return lb_lbinterpolation_get_interpolated_velocity(pos);
  });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_interpolated_velocity)

boost::optional<double>
mpi_lb_get_interpolated_density(Utils::Vector3d const &pos) {
  return detail::lb_calc_for_pos(pos, [&](auto pos) {
    return lb_lbinterpolation_get_interpolated_density(pos);
  });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_interpolated_density)

auto mpi_lb_get_density(Utils::Vector3i const &index) {
  return detail::lb_calc_fluid_kernel(
      index, [&](auto const &modes, auto const &force_density) {
        return lb_calc_density(modes, lbpar);
      });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_density)

auto mpi_lb_get_populations(Utils::Vector3i const &index) {
  return detail::lb_calc(index, [&](auto index) {
    auto const linear_index =
        get_linear_index(lblattice.local_index(index), lblattice.halo_grid);
    return lb_get_population(linear_index);
  });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_populations)

auto mpi_lb_get_boundary_flag(Utils::Vector3i const &index) {
  return detail::lb_calc(index, [&](auto index) {
#ifdef LB_BOUNDARIES
    auto const linear_index =
        get_linear_index(lblattice.local_index(index), lblattice.halo_grid);
    return lbfields[linear_index].boundary;
#else
    return false;
#endif
  });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_boundary_flag)

void mpi_lb_set_population(Utils::Vector3i const &index,
                           Utils::Vector19d const &population) {
  detail::lb_set(index, [&](auto index) {
    auto const linear_index =
        get_linear_index(lblattice.local_index(index), lblattice.halo_grid);
    lb_set_population(linear_index, population);
  });
}

REGISTER_CALLBACK(mpi_lb_set_population)

void mpi_lb_set_force_density(Utils::Vector3i const &index,
                              Utils::Vector3d const &force_density) {
  detail::lb_set(index, [&](auto index) {
    auto const linear_index =
        get_linear_index(lblattice.local_index(index), lblattice.halo_grid);
    lbfields[linear_index].force_density = force_density;
  });
}

REGISTER_CALLBACK(mpi_lb_set_force_density)

auto mpi_lb_get_momentum_density(Utils::Vector3i const &index) {
  return detail::lb_calc_fluid_kernel(
      index, [&](auto const &modes, auto const &force_density) {
        return lb_calc_momentum_density(modes, force_density);
      });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_momentum_density)

auto mpi_lb_get_pressure_tensor(Utils::Vector3i const &index) {
  return detail::lb_calc_fluid_kernel(
      index, [&](auto const &modes, auto const &force_density) {
        return lb_calc_pressure_tensor(modes, force_density, lbpar);
      });
}

REGISTER_CALLBACK_ONE_RANK(mpi_lb_get_pressure_tensor)

void mpi_bcast_lb_params_local(LBParam field, LB_Parameters const &params) {
  lbpar = params;
  lb_on_param_change(field);
}

REGISTER_CALLBACK(mpi_bcast_lb_params_local)

/** @brief Broadcast a parameter for lattice Boltzmann.
 *  @param[in] field  References the parameter field to be broadcasted.
 *                    The references are defined in lb.hpp
 */
void mpi_bcast_lb_params(LBParam field) {
  mpi_call(mpi_bcast_lb_params_local, field, lbpar);
  lb_on_param_change(field);
}

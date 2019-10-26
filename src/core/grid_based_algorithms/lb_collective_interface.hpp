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
#ifndef LB_COLLECTIVE_INTERFACE_HPP
#define LB_COLLECTIVE_INTERFACE_HPP

#include <boost/optional.hpp>
#include <utils/Vector.hpp>

/* collective getter functions */
boost::optional<Utils::Vector3d>
mpi_lb_get_interpolated_velocity(Utils::Vector3d const &pos);
boost::optional<double> mpi_lb_get_density(Utils::Vector3i const &index);
boost::optional<Utils::Vector19d>
mpi_lb_get_populations(Utils::Vector3i const &index);
boost::optional<int> mpi_lb_get_boundary_flag(Utils::Vector3i const &index);
boost::optional<Utils::Vector3d>
mpi_lb_get_momentum_density(Utils::Vector3i const &index);
boost::optional<Utils::Vector6d>
mpi_lb_get_stress(Utils::Vector3i const &index);

/* collective setter functions */
void mpi_lb_set_population(Utils::Vector3i const &index,
                           Utils::Vector19d const &population);
void mpi_lb_set_force_density(Utils::Vector3i const &index,
                              Utils::Vector3d const &force_density);

/* collective sync functions */
void mpi_bcast_lb_params(LBParam field);

#endif

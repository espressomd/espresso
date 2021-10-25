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
#ifndef LATTICE_INTERPOLATION_HPP
#define LATTICE_INTERPOLATION_HPP

#include <utils/Vector.hpp>

/**
 * @brief Interpolation order for the LB fluid interpolation.
 * @note For the CPU LB only linear interpolation is available.
 */
enum class InterpolationOrder { linear, quadratic };

/**
 * @brief Set the interpolation order for the LB.
 */
void lb_lbinterpolation_set_interpolation_order(
    InterpolationOrder const &interpolation_order);

InterpolationOrder lb_lbinterpolation_get_interpolation_order();
/**
 * @brief Calculates the fluid velocity at a given position of the
 * lattice.
 * @note It can lead to undefined behaviour if the
 * position is not within the local lattice.
 */
const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity(const Utils::Vector3d &p);

/**
 * @brief Calculates the fluid density at a given position of the lattice.
 * @note It can lead to undefined behaviour if the
 * position is not within the local lattice.
 */
double lb_lbinterpolation_get_interpolated_density(const Utils::Vector3d &p);

/**
 * @brief Add a force density to the fluid at the given position.
 */
void lb_lbinterpolation_add_force_density(const Utils::Vector3d &p,
                                          const Utils::Vector3d &force_density);
#endif

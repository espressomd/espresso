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
void mpi_set_interpolation_order_slave(int order, int);

InterpolationOrder lb_lbinterpolation_get_interpolation_order();
/**
 * @brief Calculates the fluid velocity at a given position of the
 * lattice.
 * @note It can lead to undefined behaviour if the
 * position is not within the local lattice. */
const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity(const Utils::Vector3d &p);

/**
 * @brief Calculates the interpolated fluid velocity on the master process.
 * @param pos Position at which the velocity is to be calculated.
 * @retval interpolated fluid velocity.
 */
const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity_global(const Utils::Vector3d &pos);

/**
 * @brief Add a force density to the fluid at the given position.
 */
void lb_lbinterpolation_add_force_density(const Utils::Vector3d &p,
                                          const Utils::Vector3d &force_density);
#endif

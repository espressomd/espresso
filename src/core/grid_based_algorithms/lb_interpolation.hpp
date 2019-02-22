#ifndef LATTICE_INTERPOLATION_HPP
#define LATTICE_INTERPOLATION_HPP

/**
 * @brief Calculates the fluid velocity at a given position of the
 * lattice.
 * @note It can lead to undefined behaviour if the
 * position is not within the local lattice. */
const Vector3d lb_lbinterpolation_get_interpolated_velocity(const Vector3d &p);

/**
 * @brief Calculates the interpolated fluid velocity on the master process.
 * @param pos Position at which the velocity is to be calculated.
 * @retval interpolated fluid velocity.
 */
const Vector3d
lb_lbinterpolation_get_interpolated_velocity_global(const Vector3d &pos);

/**
 * @brief Add a force density to the fluid at the given position.
 */
void lb_lbinterpolation_add_force_density(const Vector3d &p,
                                          const Vector3d &force_density);
#endif

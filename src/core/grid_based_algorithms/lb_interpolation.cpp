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
#include <boost/mpi/collectives.hpp>

#include "communication.hpp"
#include "config.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include <utils/Vector.hpp>

#include "lb.hpp"
#include "lbgpu.hpp"

namespace {
InterpolationOrder interpolation_order = InterpolationOrder::linear;
}

void mpi_set_interpolation_order(InterpolationOrder const &order) {
  interpolation_order = order;
}

REGISTER_CALLBACK(mpi_set_interpolation_order)

void lb_lbinterpolation_set_interpolation_order(
    InterpolationOrder const &order) {
  mpi_call_all(mpi_set_interpolation_order, order);
}

InterpolationOrder lb_lbinterpolation_get_interpolation_order() {
  return interpolation_order;
}

namespace {
template <typename Op>
void lattice_interpolation(Lattice const &lattice, Utils::Vector3d const &pos,
                           Op &&op) {
  Utils::Vector<std::size_t, 8> node_index{};
  Utils::Vector6d delta{};

  /* determine elementary lattice cell surrounding the particle
     and the relative position of the particle in this cell */
  lattice.map_position_to_lattice(pos, node_index, delta);
  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        auto &index = node_index[(z * 2 + y) * 2 + x];
        auto const w = delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2];

        op(index, w);
      }
    }
  }
}

Utils::Vector3d node_u(Lattice::index_t index) {
#ifdef LB_BOUNDARIES
  if (lbfields[index].boundary) {
    return lbfields[index].slip_velocity;
  }
#endif // LB_BOUNDARIES
  auto const modes = lb_calc_modes(index, lbfluid);
  auto const local_density = lbpar.density + modes[0];
  return Utils::Vector3d{modes[1], modes[2], modes[3]} / local_density;
}

} // namespace

const Utils::Vector3d
lb_lbinterpolation_get_interpolated_velocity(const Utils::Vector3d &pos) {
  Utils::Vector3d interpolated_u{};

  /* Calculate fluid velocity at particle's position.
     This is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  lattice_interpolation(lblattice, pos,
                        [&interpolated_u](Lattice::index_t index, double w) {
                          interpolated_u += w * node_u(index);
                        });

  return interpolated_u;
}

void lb_lbinterpolation_add_force_density(
    const Utils::Vector3d &pos, const Utils::Vector3d &force_density) {
  switch (interpolation_order) {
  case (InterpolationOrder::quadratic):
    throw std::runtime_error("The non-linear interpolation scheme is not "
                             "implemented for the CPU LB.");
  case (InterpolationOrder::linear):
    lattice_interpolation(lblattice, pos,
                          [&force_density](Lattice::index_t index, double w) {
                            auto &field = lbfields[index];
                            field.force_density += w * force_density;
                          });
    break;
  }
}

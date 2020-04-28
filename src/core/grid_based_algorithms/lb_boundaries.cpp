/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group,
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
/** \file
 *
 * Boundary conditions for lattice Boltzmann fluid dynamics.
 * Header file for \ref lb_boundaries.hpp.
 *
 */

#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm.hpp>

#include "grid_based_algorithms/lb_boundaries.hpp"

#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lbgpu.hpp"
#include "lbboundaries/LBBoundary.hpp"

#include <utils/index.hpp>

using Utils::get_linear_index;

#include <utils/constants.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

namespace LBBoundaries {

std::vector<std::shared_ptr<LBBoundary>> lbboundaries;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

void add(const std::shared_ptr<LBBoundary> &b) {
  lbboundaries.emplace_back(b);

  on_lbboundary_change();
}

void remove(const std::shared_ptr<LBBoundary> &b) {
  lbboundaries.erase(std::remove(lbboundaries.begin(), lbboundaries.end(), b),
                     lbboundaries.end());

  on_lbboundary_change();
}

void ek_init_boundaries() {
#if defined(CUDA) && defined(EK_BOUNDARIES)
  int number_of_boundnodes = 0;

  std::vector<ekfloat> host_wallcharge_species_density;
  float node_wallcharge = 0.0f;
  int wallcharge_species = -1, charged_boundaries = 0;
  bool node_charged = false;

  for (auto &lbboundary : lbboundaries) {
    lbboundary->set_net_charge(0.0);
  }

  if (ek_initialized) {
    host_wallcharge_species_density.resize(ek_parameters.number_of_nodes);
    for (auto &lbboundary : lbboundaries) {
      if (lbboundary->charge_density() != 0.0) {
        charged_boundaries = 1;
        break;
      }
    }

    for (int n = 0; n < int(ek_parameters.number_of_species); n++)
      if (ek_parameters.valency[n] != 0.0) {
        wallcharge_species = n;
        break;
      }

    ek_gather_wallcharge_species_density(host_wallcharge_species_density.data(),
                                         wallcharge_species);

    if (wallcharge_species == -1 && charged_boundaries) {
      runtimeErrorMsg()
          << "no charged species available to create wall charge\n";
    }

    for (int z = 0; z < int(lbpar_gpu.dim_z); z++) {
      for (int y = 0; y < int(lbpar_gpu.dim_y); y++) {
        for (int x = 0; x < int(lbpar_gpu.dim_x); x++) {
          auto const pos = static_cast<double>(lbpar_gpu.agrid) *
                           (Utils::Vector3d{1. * x, 1. * y, 1. * z} +
                            Utils::Vector3d::broadcast(0.5));
          node_charged = false;
          node_wallcharge = 0.0f;

          std::vector<std::shared_ptr<LBBoundary>> boundaries;
          std::copy_if(lbboundaries.begin(), lbboundaries.end(),
                       std::back_inserter(boundaries), [&pos](auto const lbb) {
                         return not lbb->shape().is_inside(pos);
                       });
          for (auto lbb : boundaries) {
            if ((*lbb).charge_density() != 0.0f) {
              node_charged = true;
              node_wallcharge += (*lbb).charge_density() * ek_parameters.agrid *
                                 ek_parameters.agrid * ek_parameters.agrid;
              (*lbb).set_net_charge(
                  (*lbb).net_charge() +
                  (*lbb).charge_density() * ek_parameters.agrid *
                      ek_parameters.agrid * ek_parameters.agrid);
            }
          }
          if (not boundaries.empty()) {
            number_of_boundnodes++;
          }
          ek_parameters.number_of_boundary_nodes = number_of_boundnodes;

          if (wallcharge_species != -1) {
            if (node_charged)
              host_wallcharge_species_density[ek_parameters.dim_y *
                                                  ek_parameters.dim_x * z +
                                              ek_parameters.dim_x * y + x] =
                  node_wallcharge / ek_parameters.valency[wallcharge_species];
          }
        }
      }
    }
    ek_init_species_density_wallcharge(host_wallcharge_species_density.data(),
                                       wallcharge_species);
  }
#endif
}

/** Initialize boundary conditions for all constraints in the system. */
void lb_init_boundaries() {
  if (lattice_switch == ActiveLB::GPU) {
#if defined(CUDA) && defined(LB_BOUNDARIES_GPU)
    if (this_node != 0) {
      return;
    }
    ek_init_boundaries();
    int number_of_boundnodes = 0;
    std::vector<int> host_boundary_node_list;
    std::vector<int> host_boundary_index_list;
    size_t size_of_index;

    for (int z = 0; z < int(lbpar_gpu.dim_z); z++) {
      for (int y = 0; y < int(lbpar_gpu.dim_y); y++) {
        for (int x = 0; x < int(lbpar_gpu.dim_x); x++) {
          auto const pos = static_cast<double>(lbpar_gpu.agrid) *
                           (Utils::Vector3d{1. * x, 1. * y, 1. * z} +
                            Utils::Vector3d::broadcast(0.5));

          // take last boundary containing the node
          auto const boundary = boost::find_if(
              lbboundaries | boost::adaptors::reversed, [&pos](auto const lbb) {
                return not lbb->shape().is_inside(pos);
              });

          if (boundary != boost::rend(lbboundaries)) {
            size_of_index = (number_of_boundnodes + 1) * sizeof(int);
            host_boundary_node_list.resize(size_of_index);
            host_boundary_index_list.resize(size_of_index);
            host_boundary_node_list[number_of_boundnodes] =
                x + lbpar_gpu.dim_x * y + lbpar_gpu.dim_x * lbpar_gpu.dim_y * z;
            auto const boundary_number =
                std::distance(lbboundaries.begin(), boundary.base()) - 1;
            host_boundary_index_list[number_of_boundnodes] =
                boundary_number + 1;
            number_of_boundnodes++;
          }
        }
      }
    }
    lbpar_gpu.number_of_boundnodes = number_of_boundnodes;
    /* call of cuda fkt */
    std::vector<float> boundary_velocity(3 * (lbboundaries.size() + 1));
    int n = 0;
    for (auto lbb = lbboundaries.begin(); lbb != lbboundaries.end();
         ++lbb, n++) {
      boundary_velocity[3 * n + 0] = (**lbb).velocity()[0];
      boundary_velocity[3 * n + 1] = (**lbb).velocity()[1];
      boundary_velocity[3 * n + 2] = (**lbb).velocity()[2];
    }

    boundary_velocity[3 * lbboundaries.size() + 0] = 0.0f;
    boundary_velocity[3 * lbboundaries.size() + 1] = 0.0f;
    boundary_velocity[3 * lbboundaries.size() + 2] = 0.0f;

    lb_init_boundaries_GPU(lbboundaries.size(), number_of_boundnodes,
                           host_boundary_node_list.data(),
                           host_boundary_index_list.data(),
                           boundary_velocity.data());

#endif /* defined (CUDA) && defined (LB_BOUNDARIES_GPU) */
  } else if (lattice_switch == ActiveLB::CPU) {
#if defined(LB_BOUNDARIES)
    auto const lblattice = lb_lbfluid_get_lattice();

    boost::for_each(lbfields, [](auto &f) { f.boundary = 0; });

    auto const node_pos = calc_node_pos(comm_cart);
    auto const offset = Utils::hadamard_product(node_pos, lblattice.grid);

    for (int z = 0; z < lblattice.grid[2] + 2; z++) {
      for (int y = 0; y < lblattice.grid[1] + 2; y++) {
        for (int x = 0; x < lblattice.grid[0] + 2; x++) {
          auto const pos =
              (offset + Utils::Vector3d{x - 0.5, y - 0.5, z - 0.5}) *
              lblattice.agrid;

          auto const boundary = boost::find_if(
              lbboundaries | boost::adaptors::reversed, [&pos](auto const lbb) {
                return not lbb->shape().is_inside(pos);
              });
          if (boundary != boost::rend(lbboundaries)) {
            auto const index = get_linear_index(x, y, z, lblattice.halo_grid);
            auto const boundary_number =
                std::distance(lbboundaries.begin(), boundary.base()) - 1;
            auto &node = lbfields[index];
            node.boundary = static_cast<int>(boundary_number) + 1;
            node.slip_velocity =
                (*boundary)->velocity() *
                (lb_lbfluid_get_tau() / lb_lbfluid_get_agrid());
          } else {
            lbfields[get_linear_index(x, y, z, lblattice.halo_grid)].boundary =
                0;
          }
        }
      }
    }
#endif
  }
}

Utils::Vector3d lbboundary_get_force(LBBoundary const *lbb) {
  Utils::Vector3d force{};
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  auto const it =
      boost::find_if(lbboundaries, [lbb](std::shared_ptr<LBBoundary> const &i) {
        return i.get() == lbb;
      });
  if (it == lbboundaries.end())
    throw std::runtime_error("You probably tried to get the force of an "
                             "lbboundary that was not added to "
                             "system.lbboundaries.");
  std::vector<double> forces(3 * lbboundaries.size());
  if (lattice_switch == ActiveLB::GPU) {
#if defined(LB_BOUNDARIES_GPU) && defined(CUDA)
    lb_gpu_get_boundary_forces(forces.data());
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
#if defined(LB_BOUNDARIES)
    mpi_gather_stats(8, forces.data(), nullptr, nullptr, nullptr);
#endif
  }
  auto const container_index = std::distance(lbboundaries.begin(), it);
  force[0] = forces[3 * container_index + 0];
  force[1] = forces[3 * container_index + 1];
  force[2] = forces[3 * container_index + 2];
#endif
  return force;
}

#endif /* LB_BOUNDARIES or LB_BOUNDARIES_GPU */

} // namespace LBBoundaries

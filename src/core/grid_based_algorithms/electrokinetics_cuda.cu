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

#include "cuda_wrapper.hpp"

#include "config.hpp"
#ifdef CUDA            /* Terminates at end of file */
#ifdef ELECTROKINETICS /* Terminates at end of file */

#include "cufft_wrapper.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>

#include <utils/memory.hpp>

#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "errorhandling.hpp"
#include "fd-electrostatics.cuh"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "grid_based_algorithms/lbgpu.cuh"
#include "grid_based_algorithms/lbgpu.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

/* TODO: get rid of this code duplication with lb-boundaries.h by solving the
         cuda-mpi incompatibility */

extern ActiveLB lattice_switch;
extern bool ek_initialized;
EK_parameters *lb_ek_parameters_gpu;

// Used to limit register use for the pressure calculation
#define EK_LINK_U00_pressure 0
#define EK_LINK_0U0_pressure 1
#define EK_LINK_00U_pressure 2
#define EK_LINK_D00_pressure 3
#define EK_LINK_0D0_pressure 4
#define EK_LINK_00D_pressure 5

#ifdef EK_BOUNDARIES
void LBBoundaries::lb_init_boundaries();
#endif
/* end of code duplication */

#define PI_FLOAT 3.14159265358979323846f

EK_parameters ek_parameters = {
    // agrid
    -1.0,
    // time_step
    -1.0,
    // lb_density
    -1.0,
    // dim_x
    0,
    // dim_x_padded
    0,
    // dim_y
    0,
    // dim_z
    0,
    // number_of_nodes
    0,
    // viscosity
    -1.0,
    // bulk_viscosity
    -1.0,
    // gamma_odd
    0.0,
    // gamma_even
    0.0,
    // friction
    0.0,
    // T
    -1.0,
    // prefactor
    -1.0,
    // lb_force_density
    {0.0, 0.0, 0.0},
    // number_of_species
    0,
    // reaction_species
    {-1, -1, -1},
    // rho_reactant_reservoir
    -1.0,
    // rho_product0_reservoir
    -1.0,
    // rho_product1_reservoir
    -1.0,
    // reaction_ct_rate
    -1.0,
    // reaction_fraction_0
    -1.0,
    // reaction_fraction_1
    -1.0,
    // mass_reactant
    -1.0,
    // mass_product0
    -1.0,
    // mass_product1
    -1.0,
    // stencil
    0,
    // number_of_boundary_nodes
    -1,
    // fluctuation_amplitude
    -1.0,
    // fluctuation
    false,
    // advection
    true,
    // fluidcoupling_ideal_contribution
    true,
    // es_coupling
    false,
    // charge_potential_buffer
    nullptr,
    // electric_field
    nullptr,
    // charge_potential
    nullptr,
    // j
    nullptr,
    // lb_force_density_previous
    nullptr,
    // rho
    {nullptr},
    // species_index
    {-1},
    // density
    {0.0},
    // D
    {0.0},
    // d
    {0.0},
    // valency
    {0.0},
    // ext_force_density
    {0.0},
    // node_is_catalyst
    nullptr,
};

__device__ __constant__ EK_parameters ek_parameters_gpu[1];
ekfloat *charge_gpu;
EK_parameters *ek_parameters_gpu_pointer;
LB_parameters_gpu *ek_lbparameters_gpu;
CUDA_particle_data *particle_data_gpu;
float *ek_lb_boundary_force;
char *ek_node_is_catalyst;
unsigned int old_number_of_species = 0;
unsigned int old_number_of_boundaries = 0;
Utils::Counter<uint64_t> philox_counter = Utils::Counter<uint64_t>(0);

FdElectrostatics *electrostatics = nullptr;

extern LB_parameters_gpu lbpar_gpu;
extern LB_node_force_density_gpu node_f, node_f_buf;
extern LB_nodes_gpu *current_nodes;
extern EK_parameters *lb_ek_parameters;

LB_rho_v_gpu *ek_lb_device_values;

__device__ cufftReal ek_getNode(unsigned x, unsigned y, unsigned z) {
  auto *field =
      reinterpret_cast<cufftReal *>(ek_parameters_gpu->charge_potential);
  return field[ek_parameters_gpu->dim_y * ek_parameters_gpu->dim_x_padded * z +
               ek_parameters_gpu->dim_x_padded * y + x];
}

__device__ void ek_setNode(unsigned x, unsigned y, unsigned z,
                           cufftReal value) {
  auto *field =
      reinterpret_cast<cufftReal *>(ek_parameters_gpu->charge_potential);
  field[ek_parameters_gpu->dim_y * ek_parameters_gpu->dim_x_padded * z +
        ek_parameters_gpu->dim_x_padded * y + x] = value;
}

__device__ cufftReal ek_getNode(unsigned i) {
  auto const x = i % ek_parameters_gpu->dim_x;
  i /= ek_parameters_gpu->dim_x;
  auto const y = i % ek_parameters_gpu->dim_y;
  auto const z = i / ek_parameters_gpu->dim_y;
  return ek_getNode(x, y, z);
}

__device__ void ek_setNode(unsigned i, cufftReal value) {
  auto const x = i % ek_parameters_gpu->dim_x;
  i /= ek_parameters_gpu->dim_x;
  auto const y = i % ek_parameters_gpu->dim_y;
  auto const z = i / ek_parameters_gpu->dim_y;
  ek_setNode(x, y, z, value);
}

__device__ unsigned int ek_getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x +
         threadIdx.x;
}

__device__ void rhoindex_linear2cartesian(unsigned int index,
                                          unsigned int *coord) {

  coord[0] = index % ek_parameters_gpu->dim_x;
  index /= ek_parameters_gpu->dim_x;
  coord[1] = index % ek_parameters_gpu->dim_y;
  coord[2] = index / ek_parameters_gpu->dim_y;
}

__device__ unsigned int
rhoindex_cartesian2linear(unsigned int x, unsigned int y, unsigned int z) {

  return z * ek_parameters_gpu->dim_y * ek_parameters_gpu->dim_x +
         y * ek_parameters_gpu->dim_x + x;
}

__device__ unsigned int rhoindex_cartesian2linear_padded(unsigned int x,
                                                         unsigned int y,
                                                         unsigned int z) {

  return z * ek_parameters_gpu->dim_y * ek_parameters_gpu->dim_x_padded +
         y * ek_parameters_gpu->dim_x_padded + x;
}

__device__ void jindex_linear2cartesian(unsigned int index, unsigned int *coord,
                                        unsigned int *c) {

  coord[0] = index % ek_parameters_gpu->dim_x;
  index /= ek_parameters_gpu->dim_x;
  coord[1] = index % ek_parameters_gpu->dim_y;
  index /= ek_parameters_gpu->dim_y;
  coord[2] = index % ek_parameters_gpu->dim_z;
  *c = index / ek_parameters_gpu->dim_z;
}

__device__ unsigned int jindex_cartesian2linear(unsigned int x, unsigned int y,
                                                unsigned int z,
                                                unsigned int c) {

  return c * ek_parameters_gpu->number_of_nodes +
         z * ek_parameters_gpu->dim_y * ek_parameters_gpu->dim_x +
         y * ek_parameters_gpu->dim_x + x;
}

// TODO fluxindex fastest running might improve caching
__device__ unsigned int jindex_getByRhoLinear(unsigned int rho_index,
                                              unsigned int c) {

  return c * ek_parameters_gpu->number_of_nodes + rho_index;
}

__device__ void ek_displacement(float *dx, LB_nodes_gpu n,
                                unsigned int node_index,
                                LB_parameters_gpu *ek_lbparameters_gpu) {

  float rho = ek_lbparameters_gpu->rho * ek_lbparameters_gpu->agrid *
              ek_lbparameters_gpu->agrid * ek_lbparameters_gpu->agrid;

  float mode[19];

  for (int i = 0; i < 19; i++) {
    mode[i] = n.vd[i * ek_lbparameters_gpu->number_of_nodes + node_index];
  }

  rho += mode[0] + mode[1] + mode[2] + mode[3] + mode[4] + mode[5] + mode[6] +
         mode[7] + mode[8] + mode[9] + mode[10] + mode[11] + mode[12] +
         mode[13] + mode[14] + mode[15] + mode[16] + mode[17] + mode[18];

  dx[0] = (mode[1] - mode[2]) + (mode[7] - mode[8]) + (mode[9] - mode[10]) +
          (mode[11] - mode[12]) + (mode[13] - mode[14]);

  dx[1] = (mode[3] - mode[4]) + (mode[7] - mode[8]) - (mode[9] - mode[10]) +
          (mode[15] - mode[16]) + (mode[17] - mode[18]);

  dx[2] = (mode[5] - mode[6]) + (mode[11] - mode[12]) - (mode[13] - mode[14]) +
          (mode[15] - mode[16]) - (mode[17] - mode[18]);

  // Velocity requires half the force_density in the previous time step

  dx[0] += 0.5f * ek_parameters_gpu->lb_force_density_previous[node_index];
  dx[1] += 0.5f *
           ek_parameters_gpu
               ->lb_force_density_previous[ek_parameters_gpu->number_of_nodes +
                                           node_index];
  dx[2] +=
      0.5f *
      ek_parameters_gpu
          ->lb_force_density_previous[2 * ek_parameters_gpu->number_of_nodes +
                                      node_index];

  dx[0] *= 1.0f / rho;
  dx[1] *= 1.0f / rho;
  dx[2] *= 1.0f / rho;
}

__device__ void ek_diffusion_migration_lbforce_linkcentered_stencil(
    unsigned int index, unsigned int index_padded,
    unsigned int const *neighborindex, unsigned int const *neighborindex_padded,
    unsigned int species_index, LB_node_force_density_gpu node_f,
    LB_nodes_gpu lb_node) {
  ekfloat flux, force;

  float agrid_inv = 1.0f / ek_parameters_gpu->agrid;
  float sqrt2agrid_inv = 1.0f / (sqrtf(2.0f) * ek_parameters_gpu->agrid);
  float sqrt2_inv = 1.0f / sqrtf(2.0f);
  float twoT_inv = 1.0f / (2.0f * ek_parameters_gpu->T);
  float D_inv = 1.0f / ek_parameters_gpu->D[species_index];
  float force_conv =
      agrid_inv * ek_parameters_gpu->time_step * ek_parameters_gpu->time_step;

  // face in x
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U00]]) *
         agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_U00]]) *
           agrid_inv +
       ek_parameters_gpu->ext_force_density[0][species_index]);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U00]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_U00]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U00)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[index], force * 0.5f);
    atomicAdd(&node_f.force_density[neighborindex[EK_LINK_U00]], force * 0.5f);
  } else {
    force = -1.0f * ek_parameters_gpu->valency[species_index] *
            (((cufftReal *)ek_parameters_gpu
                  ->charge_potential)[neighborindex_padded[EK_LINK_U00]] -
             ((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded]) *
            agrid_inv;

    force *= force_conv;

    atomicAdd(&node_f.force_density[index],
              ek_parameters_gpu->rho[species_index][index] *
                  (force * 0.5f +
                   ek_parameters_gpu->ext_force_density[0][species_index] *
                       force_conv));
  }

  // face in y
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0U0]]) *
         agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_0U0]]) *
           agrid_inv +
       ek_parameters_gpu->ext_force_density[1][species_index]);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0U0]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_0U0]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
              force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_0U0]],
              force * 0.5f);
  } else {
    force = -1.0f * ek_parameters_gpu->valency[species_index] *
            (((cufftReal *)ek_parameters_gpu
                  ->charge_potential)[neighborindex_padded[EK_LINK_0U0]] -
             ((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded]) *
            agrid_inv;

    force *= force_conv;

    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
              ek_parameters_gpu->rho[species_index][index] *
                  (force * 0.5f +
                   ek_parameters_gpu->ext_force_density[1][species_index] *
                       force_conv));

    atomicAdd(
        &node_f.force_density[ek_parameters_gpu->number_of_nodes +
                              neighborindex[EK_LINK_0U0]],
        ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0U0]] *
            force * 0.5f);
  }

  // face in z
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_00U]]) *
         agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_00U]]) *
           agrid_inv +
       ek_parameters_gpu->ext_force_density[2][species_index]);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_00U]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_00U]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_00U)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
        force * 0.5f);
    atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_00U]],
              force * 0.5f);
  } else {
    force = -1.0f * ek_parameters_gpu->valency[species_index] *
            (((cufftReal *)ek_parameters_gpu
                  ->charge_potential)[neighborindex_padded[EK_LINK_00U]] -
             ((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded]) *
            agrid_inv;

    force *= force_conv;

    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
        ek_parameters_gpu->rho[species_index][index] *
            (force * 0.5f +
             ek_parameters_gpu->ext_force_density[2][species_index] *
                 force_conv));

    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                              neighborindex[EK_LINK_00U]],
        ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_00U]] *
            force * 0.5f);
  }

  // edge in z
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UU0]]) *
         sqrt2agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_UU0]]) *
           sqrt2agrid_inv +
       (ek_parameters_gpu->ext_force_density[0][species_index] +
        ek_parameters_gpu->ext_force_density[1][species_index]) *
           sqrt2_inv);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UU0]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_UU0]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[index], force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
              force * 0.5f);
    atomicAdd(&node_f.force_density[neighborindex[EK_LINK_UU0]], force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_UU0]],
              force * 0.5f);
  }

  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UD0]]) *
         sqrt2agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_UD0]]) *
           sqrt2agrid_inv +
       (ek_parameters_gpu->ext_force_density[0][species_index] -
        ek_parameters_gpu->ext_force_density[1][species_index]) *
           sqrt2_inv);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UD0]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_UD0]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[index], force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
              -force * 0.5f);
    atomicAdd(&node_f.force_density[neighborindex[EK_LINK_UD0]], force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_UD0]],
              -force * 0.5f);
  }

  // edge in y
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0U]]) *
         sqrt2agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_U0U]]) *
           sqrt2agrid_inv +
       (ek_parameters_gpu->ext_force_density[0][species_index] +
        ek_parameters_gpu->ext_force_density[2][species_index]) *
           sqrt2_inv);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0U]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_U0U]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[index], force * 0.5f);
    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
        force * 0.5f);
    atomicAdd(&node_f.force_density[neighborindex[EK_LINK_U0U]], force * 0.5f);
    atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_U0U]],
              force * 0.5f);
  }

  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0D]]) *
         sqrt2agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_U0D]]) *
           sqrt2agrid_inv +
       (ek_parameters_gpu->ext_force_density[0][species_index] -
        ek_parameters_gpu->ext_force_density[2][species_index]) *
           sqrt2_inv);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0D]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_U0D]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[index], force * 0.5f);
    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
        -force * 0.5f);
    atomicAdd(&node_f.force_density[neighborindex[EK_LINK_U0D]], force * 0.5f);
    atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_U0D]],
              -force * 0.5f);
  }

  // edge in x
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UU]]) *
         sqrt2agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_0UU]]) *
           sqrt2agrid_inv +
       (ek_parameters_gpu->ext_force_density[1][species_index] +
        ek_parameters_gpu->ext_force_density[2][species_index]) *
           sqrt2_inv);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UU]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_0UU]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
              force * 0.5f);
    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
        force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_0UU]],
              force * 0.5f);
    atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_0UU]],
              force * 0.5f);
  }

  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UD]]) *
         sqrt2agrid_inv;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_0UD]]) *
           sqrt2agrid_inv +
       (ek_parameters_gpu->ext_force_density[1][species_index] -
        ek_parameters_gpu->ext_force_density[2][species_index]) *
           sqrt2_inv);

  flux += force *
          (ek_parameters_gpu->rho[species_index][index] +
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UD]]) *
          twoT_inv;

  flux *= ek_parameters_gpu->d[species_index] * agrid_inv;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_0UD]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
            flux * ek_parameters_gpu->time_step);

  if (ek_parameters_gpu->fluidcoupling_ideal_contribution) {
    force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid * D_inv;
    force *= force_conv;

    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
              force * 0.5f);
    atomicAdd(
        &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
        -force * 0.5f);
    atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_0UD]],
              force * 0.5f);
    atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                    neighborindex[EK_LINK_0UD]],
              -force * 0.5f);
  }
}

__device__ void ek_diffusion_migration_lbforce_nodecentered_stencil(
    unsigned int index, unsigned int index_padded,
    unsigned int const *neighborindex, unsigned int const *neighborindex_padded,
    unsigned int species_index, LB_node_force_density_gpu node_f,
    LB_nodes_gpu lb_node) {
  ekfloat flux, force;

  // face in x
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U00]]) /
         ek_parameters_gpu->agrid;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_U00]]) /
           ek_parameters_gpu->agrid +
       ek_parameters_gpu->ext_force_density[0][species_index]);

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U00]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_U00]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U00)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[index], force / 2.0f);
  atomicAdd(&node_f.force_density[neighborindex[EK_LINK_U00]], force / 2.0f);

  // face in y
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0U0]]) /
         ek_parameters_gpu->agrid;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_0U0]]) /
           ek_parameters_gpu->agrid +
       ek_parameters_gpu->ext_force_density[1][species_index]);

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0U0]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_0U0]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0U0)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
            force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_0U0]],
            force / 2.0f);

  // face in z
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_00U]]) /
         ek_parameters_gpu->agrid;

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_00U]]) /
           ek_parameters_gpu->agrid +
       ek_parameters_gpu->ext_force_density[2][species_index]);

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_00U]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_00U]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_00U)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(
      &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
      force / 2.0f);
  atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_00U]],
            force / 2.0f);

  // edge in z
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UU0]]) /
         (sqrtf(2.0f) * ek_parameters_gpu->agrid);

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_UU0]]) /
           (sqrtf(2.0f) * ek_parameters_gpu->agrid) +
       (ek_parameters_gpu->ext_force_density[0][species_index] +
        ek_parameters_gpu->ext_force_density[1][species_index]) /
           sqrtf(2.0f));

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UU0]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_UU0]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UU0)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[index], force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
            force / 2.0f);
  atomicAdd(&node_f.force_density[neighborindex[EK_LINK_UU0]], force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_UU0]],
            force / 2.0f);

  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UD0]]) /
         (sqrtf(2.0f) * ek_parameters_gpu->agrid);

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_UD0]]) /
           (sqrtf(2.0f) * ek_parameters_gpu->agrid) +
       (ek_parameters_gpu->ext_force_density[0][species_index] -
        ek_parameters_gpu->ext_force_density[1][species_index]) /
           sqrtf(2.0f));

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_UD0]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_UD0]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UD0)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[index], force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
            -force / 2.0f);
  atomicAdd(&node_f.force_density[neighborindex[EK_LINK_UD0]], force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_UD0]],
            -force / 2.0f);

  // edge in y
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0U]]) /
         (sqrtf(2.0f) * ek_parameters_gpu->agrid);

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_U0U]]) /
           (sqrtf(2.0f) * ek_parameters_gpu->agrid) +
       (ek_parameters_gpu->ext_force_density[0][species_index] +
        ek_parameters_gpu->ext_force_density[2][species_index]) /
           sqrtf(2.0f));

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0U]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_U0U]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U0U)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[index], force / 2.0f);
  atomicAdd(
      &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
      force / 2.0f);
  atomicAdd(&node_f.force_density[neighborindex[EK_LINK_U0U]], force / 2.0f);
  atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_U0U]],
            force / 2.0f);

  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0D]]) /
         (sqrtf(2.0f) * ek_parameters_gpu->agrid);

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_U0D]]) /
           (sqrtf(2.0f) * ek_parameters_gpu->agrid) +
       (ek_parameters_gpu->ext_force_density[0][species_index] -
        ek_parameters_gpu->ext_force_density[2][species_index]) /
           sqrtf(2.0f));

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_U0D]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_U0D]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U0D)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[index], force / 2.0f);
  atomicAdd(
      &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
      -force / 2.0f);
  atomicAdd(&node_f.force_density[neighborindex[EK_LINK_U0D]], force / 2.0f);
  atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_U0D]],
            -force / 2.0f);

  // edge in x
  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UU]]) /
         (sqrtf(2.0f) * ek_parameters_gpu->agrid);

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_0UU]]) /
           (sqrtf(2.0f) * ek_parameters_gpu->agrid) +
       (ek_parameters_gpu->ext_force_density[1][species_index] +
        ek_parameters_gpu->ext_force_density[2][species_index]) /
           sqrtf(2.0f));

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UU]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_0UU]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0UU)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
            force / 2.0f);
  atomicAdd(
      &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
      force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_0UU]],
            force / 2.0f);
  atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_0UU]],
            force / 2.0f);

  flux = (ek_parameters_gpu->rho[species_index][index] -
          ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UD]]) /
         (sqrtf(2.0f) * ek_parameters_gpu->agrid);

  force =
      (ek_parameters_gpu->valency[species_index] *
           (((cufftReal *)ek_parameters_gpu->charge_potential)[index_padded] -
            ((cufftReal *)ek_parameters_gpu
                 ->charge_potential)[neighborindex_padded[EK_LINK_0UD]]) /
           (sqrtf(2.0f) * ek_parameters_gpu->agrid) +
       (ek_parameters_gpu->ext_force_density[1][species_index] -
        ek_parameters_gpu->ext_force_density[2][species_index]) /
           sqrtf(2.0f));

  flux +=
      force *
      (static_cast<ekfloat>(force >= 0.0f) *
           ek_parameters_gpu->rho[species_index][index] +
       static_cast<ekfloat>(force < 0.0f) *
           ek_parameters_gpu->rho[species_index][neighborindex[EK_LINK_0UD]]) /
      ek_parameters_gpu->T;

  flux *= ek_parameters_gpu->d[species_index] / ek_parameters_gpu->agrid;

  flux *= static_cast<ekfloat>(!(lb_node.boundary[index] ||
                                 lb_node.boundary[neighborindex[EK_LINK_0UD]]));

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0UD)],
            flux * ek_parameters_gpu->time_step);

  force = flux * ek_parameters_gpu->T * ek_parameters_gpu->agrid /
          ek_parameters_gpu->D[species_index];

  force *= powf(ek_parameters_gpu->agrid, -1) * ek_parameters_gpu->time_step *
           ek_parameters_gpu->time_step;

  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes + index],
            force / 2.0f);
  atomicAdd(
      &node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index],
      -force / 2.0f);
  atomicAdd(&node_f.force_density[ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_0UD]],
            force / 2.0f);
  atomicAdd(&node_f.force_density[2 * ek_parameters_gpu->number_of_nodes +
                                  neighborindex[EK_LINK_0UD]],
            -force / 2.0f);
}

__device__ void
ek_add_advection_to_flux(unsigned int index, unsigned int *neighborindex,
                         unsigned int *coord, unsigned int species_index,
                         LB_node_force_density_gpu node_f, LB_nodes_gpu lb_node,
                         LB_parameters_gpu *ek_lbparameters_gpu) {
  float dx[3];
  int di[3];
  unsigned int node;

  ek_displacement(dx, lb_node, index, ek_lbparameters_gpu);

  di[0] = 1 - signbit(dx[0]);
  di[1] = 1 - signbit(dx[1]);
  di[2] = 1 - signbit(dx[2]);

  dx[0] = fabs(dx[0]);
  dx[1] = fabs(dx[1]);
  dx[2] = fabs(dx[2]);

  unsigned int target_node[3];
  unsigned int target_node_index;
  int not_boundary;

  // face in x
  node = rhoindex_cartesian2linear(
      (coord[0] + di[0] - 1 + ek_parameters_gpu->dim_x) %
          ek_parameters_gpu->dim_x,
      coord[1], coord[2]);

  target_node[0] = (coord[0] + 2 * static_cast<unsigned>(di[0]) - 1 +
                    ek_parameters_gpu->dim_x) %
                   ek_parameters_gpu->dim_x;
  target_node[1] = coord[1];
  target_node[2] = coord[2];
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(node, EK_LINK_U00)],
            (2 * static_cast<ekfloat>(di[0]) - 1) *
                ek_parameters_gpu->rho[species_index][index] * dx[0] *
                (1.0f - dx[1]) * (1.0f - dx[2]) *
                static_cast<ekfloat>(not_boundary));

  // face in y
  node = rhoindex_cartesian2linear(
      coord[0],
      (coord[1] + di[1] - 1 + ek_parameters_gpu->dim_y) %
          ek_parameters_gpu->dim_y,
      coord[2]);

  target_node[0] = coord[0];
  target_node[1] = (coord[1] + 2 * static_cast<unsigned>(di[1]) - 1 +
                    ek_parameters_gpu->dim_y) %
                   ek_parameters_gpu->dim_y;
  target_node[2] = coord[2];
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(node, EK_LINK_0U0)],
            (2 * static_cast<ekfloat>(di[1]) - 1) *
                ek_parameters_gpu->rho[species_index][index] * (1.0f - dx[0]) *
                dx[1] * (1.0f - dx[2]) * static_cast<ekfloat>(not_boundary));

  // face in z
  node = rhoindex_cartesian2linear(
      coord[0], coord[1],
      (coord[2] + di[2] - 1 + ek_parameters_gpu->dim_z) %
          ek_parameters_gpu->dim_z);

  target_node[0] = coord[0];
  target_node[1] = coord[1];
  target_node[2] = (coord[2] + 2 * static_cast<unsigned>(di[2]) - 1 +
                    ek_parameters_gpu->dim_z) %
                   ek_parameters_gpu->dim_z;
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(node, EK_LINK_00U)],
            (2 * static_cast<ekfloat>(di[2]) - 1) *
                ek_parameters_gpu->rho[species_index][index] * (1.0f - dx[0]) *
                (1.0f - dx[1]) * dx[2] * static_cast<ekfloat>(not_boundary));

  // edge in x
  node = rhoindex_cartesian2linear(
      coord[0],
      (coord[1] + di[1] - 1 + ek_parameters_gpu->dim_y) %
          ek_parameters_gpu->dim_y,
      (coord[2] + (1 - di[1]) * (2 * di[2] - 1) + ek_parameters_gpu->dim_z) %
          ek_parameters_gpu->dim_z);

  target_node[0] = coord[0];
  target_node[1] = (coord[1] + 2 * static_cast<unsigned>(di[1]) - 1 +
                    ek_parameters_gpu->dim_y) %
                   ek_parameters_gpu->dim_y;
  target_node[2] = (coord[2] + 2 * static_cast<unsigned>(di[2]) - 1 +
                    ek_parameters_gpu->dim_z) %
                   ek_parameters_gpu->dim_z;
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(
      &ek_parameters_gpu
           ->j[jindex_getByRhoLinear(node, EK_LINK_0UU + (di[1] + di[2] == 1))],
      (2 * static_cast<ekfloat>(di[1]) - 1) *
          ek_parameters_gpu->rho[species_index][index] * (1.0f - dx[0]) *
          dx[1] * dx[2] * static_cast<ekfloat>(not_boundary));

  // edge in y
  node = rhoindex_cartesian2linear(
      (coord[0] + di[0] - 1 + ek_parameters_gpu->dim_x) %
          ek_parameters_gpu->dim_x,
      coord[1],
      (coord[2] + (1 - di[0]) * (2 * di[2] - 1) + ek_parameters_gpu->dim_z) %
          ek_parameters_gpu->dim_z);

  target_node[0] = (coord[0] + 2 * static_cast<unsigned>(di[0]) - 1 +
                    ek_parameters_gpu->dim_x) %
                   ek_parameters_gpu->dim_x;
  target_node[1] = coord[1];
  target_node[2] = (coord[2] + 2 * static_cast<unsigned>(di[2]) - 1 +
                    ek_parameters_gpu->dim_z) %
                   ek_parameters_gpu->dim_z;
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(
      &ek_parameters_gpu
           ->j[jindex_getByRhoLinear(node, EK_LINK_U0U + (di[0] + di[2] == 1))],
      (2 * static_cast<ekfloat>(di[0]) - 1) *
          ek_parameters_gpu->rho[species_index][index] * dx[0] *
          (1.0f - dx[1]) * dx[2] * static_cast<ekfloat>(not_boundary));

  // edge in z
  node = rhoindex_cartesian2linear(
      (coord[0] + di[0] - 1 + ek_parameters_gpu->dim_x) %
          ek_parameters_gpu->dim_x,
      (coord[1] + (1 - di[0]) * (2 * di[1] - 1) + ek_parameters_gpu->dim_y) %
          ek_parameters_gpu->dim_y,
      coord[2]);

  target_node[0] = (coord[0] + 2 * static_cast<unsigned>(di[0]) - 1 +
                    ek_parameters_gpu->dim_x) %
                   ek_parameters_gpu->dim_x;
  target_node[1] = (coord[1] + 2 * static_cast<unsigned>(di[1]) - 1 +
                    ek_parameters_gpu->dim_y) %
                   ek_parameters_gpu->dim_y;
  target_node[2] = coord[2];
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(
      &ek_parameters_gpu
           ->j[jindex_getByRhoLinear(node, EK_LINK_UU0 + (di[0] + di[1] == 1))],
      (2 * static_cast<ekfloat>(di[0]) - 1) *
          ek_parameters_gpu->rho[species_index][index] * dx[0] * dx[1] *
          (1.0f - dx[2]) * static_cast<ekfloat>(not_boundary));

  // corner
  node = rhoindex_cartesian2linear(
      (coord[0] + di[0] - 1 + ek_parameters_gpu->dim_x) %
          ek_parameters_gpu->dim_x,
      (coord[1] + (1 - di[0]) * (2 * di[1] - 1) + ek_parameters_gpu->dim_y) %
          ek_parameters_gpu->dim_y,
      (coord[2] + (1 - di[0]) * (2 * di[2] - 1) + ek_parameters_gpu->dim_z) %
          ek_parameters_gpu->dim_z);

  target_node[0] = (coord[0] + 2 * static_cast<unsigned>(di[0]) - 1 +
                    ek_parameters_gpu->dim_x) %
                   ek_parameters_gpu->dim_x;
  target_node[1] = (coord[1] + 2 * static_cast<unsigned>(di[1]) - 1 +
                    ek_parameters_gpu->dim_y) %
                   ek_parameters_gpu->dim_y;
  target_node[2] = (coord[2] + 2 * static_cast<unsigned>(di[2]) - 1 +
                    ek_parameters_gpu->dim_z) %
                   ek_parameters_gpu->dim_z;
  target_node_index =
      rhoindex_cartesian2linear(target_node[0], target_node[1], target_node[2]);
  not_boundary =
      (lb_node.boundary[index] || lb_node.boundary[target_node_index]) == 0;

  atomicAdd(&ek_parameters_gpu->j[jindex_getByRhoLinear(
                node, (1 - di[0]) * (EK_LINK_UUU + 2 * di[1] + di[2]) +
                          di[0] * (EK_LINK_UDD - 2 * di[1] - di[2]))],
            static_cast<float>(2 * di[0] - 1) *
                ek_parameters_gpu->rho[species_index][index] * dx[0] * dx[1] *
                dx[2] * static_cast<float>(not_boundary));
}

__device__ float4 ek_random_wrapper_philox(unsigned int index,
                                           unsigned int mode,
                                           uint64_t philox_counter) {
  // Split the 64 bit counter into two 32 bit ints.
  auto const philox_counter_hi = static_cast<uint32_t>(philox_counter >> 32);
  auto const philox_counter_low = static_cast<uint32_t>(philox_counter);
  uint4 rnd_ints =
      curand_Philox4x32_10(make_uint4(index, philox_counter_hi, 0, mode),
                           make_uint2(philox_counter_low, 0));
  float4 rnd_floats;
  rnd_floats.w = static_cast<float>(rnd_ints.w) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  rnd_floats.x = static_cast<float>(rnd_ints.x) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  rnd_floats.y = static_cast<float>(rnd_ints.y) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  rnd_floats.z = static_cast<float>(rnd_ints.z) * CURAND_2POW32_INV +
                 (CURAND_2POW32_INV / 2.0f);
  return rnd_floats;
}

__device__ void ek_add_fluctuations_to_flux(unsigned int index,
                                            unsigned int species_index,
                                            unsigned int const *neighborindex,
                                            LB_nodes_gpu lb_node,
                                            uint64_t philox_counter) {
  if (index < ek_parameters_gpu->number_of_nodes) {
    float density = ek_parameters_gpu->rho[species_index][index];
    float *flux = ek_parameters_gpu->j;
    float diffusion = ek_parameters_gpu->D[species_index];
    float time_step = ek_parameters_gpu->time_step;
    float agrid = ek_parameters_gpu->agrid;
    float4 random_floats;
    float random;

#ifdef EK_DEBUG
    float *flux_fluc = ek_parameters_gpu->j_fluc;
#endif
    float fluc = 0.0f;

    for (int i = 0; i < 9; i++) {

      if (i % 4 == 0) {
        random_floats = ek_random_wrapper_philox(index, i + 40, philox_counter);
        random = (random_floats.w - 0.5f) * 2.0f;
      } else if (i % 4 == 1) {
        random = (random_floats.x - 0.5f) * 2.0f;
      } else if (i % 4 == 2) {
        random = (random_floats.y - 0.5f) * 2.0f;
      } else if (i % 4 == 3) {
        random = (random_floats.z - 0.5f) * 2.0f;
      }
      float H = 0.0f;
      float HN = 0.0f;
      float neighbor_density =
          ek_parameters_gpu->rho[species_index][neighborindex[i]];

      H = static_cast<float>(density >= 0.0f) * min(density, 1.0f);
      HN = static_cast<float>(neighbor_density >= 0.0f) *
           min(neighbor_density, 1.0f);

      float average_density = H * HN * (density + neighbor_density) / 2.0f;

      if (i > 2) {
        fluc = 1.0f *
               powf(2.0f * average_density * diffusion * time_step /
                        (agrid * agrid),
                    0.5f) *
               random * ek_parameters_gpu->fluctuation_amplitude / sqrtf(2.0f);
        fluc *= static_cast<float>(
            !(lb_node.boundary[index] || lb_node.boundary[neighborindex[i]]));
#ifdef EK_DEBUG
        flux_fluc[jindex_getByRhoLinear(index, i)] = fluc;
#endif
        flux[jindex_getByRhoLinear(index, i)] += fluc;
      } else {
        fluc = 1.0f *
               powf(2.0f * average_density * diffusion * time_step /
                        (agrid * agrid),
                    0.5f) *
               random * ek_parameters_gpu->fluctuation_amplitude;
        fluc *= static_cast<float>(
            !(lb_node.boundary[index] || lb_node.boundary[neighborindex[i]]));
#ifdef EK_DEBUG
        flux_fluc[jindex_getByRhoLinear(index, i)] = fluc;
#endif
        flux[jindex_getByRhoLinear(index, i)] += fluc;
      }
    }
  }
}

__global__ void ek_calculate_quantities(unsigned int species_index,
                                        LB_nodes_gpu lb_node,
                                        LB_node_force_density_gpu node_f,
                                        LB_parameters_gpu *ek_lbparameters_gpu,
                                        LB_rho_v_gpu *d_v,
                                        uint64_t philox_counter) {

  unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {

    unsigned int coord[3];
    unsigned int neighborindex[9];
    unsigned int neighborindex_padded[9];
    unsigned int index_padded;

    rhoindex_linear2cartesian(index, coord);

    /* Calculate the diffusive fluxes between this node and its neighbors. Only
       the 9 fluxes along the directions of the LB velocities c_i with i odd are
       stored with a node to avoid redundancies. */

    neighborindex[EK_LINK_U00] = rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1], coord[2]);

    neighborindex[EK_LINK_0U0] = rhoindex_cartesian2linear(
        coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2]);

    neighborindex[EK_LINK_00U] = rhoindex_cartesian2linear(
        coord[0], coord[1], (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_UU0] = rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu->dim_x,
        (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2]);

    neighborindex[EK_LINK_UD0] = rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu->dim_x,
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        coord[2]);

    neighborindex[EK_LINK_U0U] = rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1],
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_U0D] = rhoindex_cartesian2linear(
        (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1],
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_0UU] = rhoindex_cartesian2linear(
        coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y,
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_0UD] = rhoindex_cartesian2linear(
        coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y,
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    /* calculate the same indices respecting the FFT padding */

    index_padded =
        rhoindex_cartesian2linear_padded(coord[0], coord[1], coord[2]);

    neighborindex_padded[EK_LINK_U00] = rhoindex_cartesian2linear_padded(
        (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1], coord[2]);

    neighborindex_padded[EK_LINK_0U0] = rhoindex_cartesian2linear_padded(
        coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2]);

    neighborindex_padded[EK_LINK_00U] = rhoindex_cartesian2linear_padded(
        coord[0], coord[1], (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex_padded[EK_LINK_UU0] = rhoindex_cartesian2linear_padded(
        (coord[0] + 1) % ek_parameters_gpu->dim_x,
        (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2]);

    neighborindex_padded[EK_LINK_UD0] = rhoindex_cartesian2linear_padded(
        (coord[0] + 1) % ek_parameters_gpu->dim_x,
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        coord[2]);

    neighborindex_padded[EK_LINK_U0U] = rhoindex_cartesian2linear_padded(
        (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1],
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex_padded[EK_LINK_U0D] = rhoindex_cartesian2linear_padded(
        (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1],
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex_padded[EK_LINK_0UU] = rhoindex_cartesian2linear_padded(
        coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y,
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex_padded[EK_LINK_0UD] = rhoindex_cartesian2linear_padded(
        coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y,
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    /* diffusive contribution to flux and LB force_density*/
    if (ek_parameters_gpu->stencil == 0) // link centered
      ek_diffusion_migration_lbforce_linkcentered_stencil(
          index, index_padded, neighborindex, neighborindex_padded,
          species_index, node_f, lb_node);
    else if (ek_parameters_gpu->stencil == 1) // node centered
      ek_diffusion_migration_lbforce_nodecentered_stencil(
          index, index_padded, neighborindex, neighborindex_padded,
          species_index, node_f, lb_node);

    /* advective contribution to flux */
    if (ek_parameters_gpu->advection)
      ek_add_advection_to_flux(index, neighborindex, coord, species_index,
                               node_f, lb_node, ek_lbparameters_gpu);

    /* fluctuation contribution to flux */
    if (ek_parameters_gpu->fluctuations)
      ek_add_fluctuations_to_flux(index, species_index, neighborindex, lb_node,
                                  philox_counter);
  }
}

__global__ void ek_propagate_densities(unsigned int species_index) {

  unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    unsigned int neighborindex[13];
    unsigned int coord[3];

    rhoindex_linear2cartesian(index, coord);

    /* Indices of the neighbors storing the other half
       of the fluxes associated with this link */
    neighborindex[EK_LINK_D00 - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        coord[1], coord[2]);

    neighborindex[EK_LINK_0D0 - 13] = rhoindex_cartesian2linear(
        coord[0],
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        coord[2]);

    neighborindex[EK_LINK_00D - 13] = rhoindex_cartesian2linear(
        coord[0], coord[1],
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_DD0 - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        coord[2]);

    neighborindex[EK_LINK_DU0 - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2]);

    neighborindex[EK_LINK_D0D - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        coord[1],
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_D0U - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        coord[1], (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_0DD - 13] = rhoindex_cartesian2linear(
        coord[0],
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_0DU - 13] = rhoindex_cartesian2linear(
        coord[0],
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_DDD - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_DDU - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_DUD - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        (coord[1] + 1) % ek_parameters_gpu->dim_y,
        (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

    neighborindex[EK_LINK_DUU - 13] = rhoindex_cartesian2linear(
        (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
        (coord[1] + 1) % ek_parameters_gpu->dim_y,
        (coord[2] + 1) % ek_parameters_gpu->dim_z);

    /* Calculate change of densities due to diffusive fluxes */
    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U00)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_D00 - 13], EK_LINK_U00)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0U0)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_0D0 - 13], EK_LINK_0U0)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_00U)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_00D - 13], EK_LINK_00U)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UU0)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_DD0 - 13], EK_LINK_UU0)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UD0)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_DU0 - 13], EK_LINK_UD0)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U0U)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_D0D - 13], EK_LINK_U0U)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_U0D)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_D0U - 13], EK_LINK_U0D)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0UU)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_0DD - 13], EK_LINK_0UU)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_0UD)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_0DU - 13], EK_LINK_0UD)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UUU)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_DDD - 13], EK_LINK_UUU)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UUD)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_DDU - 13], EK_LINK_UUD)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UDU)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_DUD - 13], EK_LINK_UDU)];

    ek_parameters_gpu->rho[species_index][index] -=
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, EK_LINK_UDD)];
    ek_parameters_gpu->rho[species_index][index] +=
        ek_parameters_gpu->j[jindex_getByRhoLinear(
            neighborindex[EK_LINK_DUU - 13], EK_LINK_UDD)];
  }
}

__global__ void ek_apply_boundaries(unsigned int species_index,
                                    LB_nodes_gpu lbnode,
                                    LB_node_force_density_gpu node_f) {

  unsigned int index = ek_getThreadIndex();
  unsigned int neighborindex[22];
  unsigned int coord[3];

  if (index < ek_parameters_gpu->number_of_nodes) {
    if (lbnode.boundary[index]) {

      rhoindex_linear2cartesian(index, coord);

      /* Indices of the neighbors */
      neighborindex[EK_LINK_D00 - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          coord[1], coord[2]);

      neighborindex[EK_LINK_0D0 - 13] = rhoindex_cartesian2linear(
          coord[0],
          (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
          coord[2]);

      neighborindex[EK_LINK_00D - 13] = rhoindex_cartesian2linear(
          coord[0], coord[1],
          (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_DD0 - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
          coord[2]);

      neighborindex[EK_LINK_DU0 - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2]);

      neighborindex[EK_LINK_D0D - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          coord[1],
          (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_D0U - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          coord[1], (coord[2] + 1) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_0DD - 13] = rhoindex_cartesian2linear(
          coord[0],
          (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
          (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_0DU - 13] = rhoindex_cartesian2linear(
          coord[0],
          (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
          (coord[2] + 1) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_DDD - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
          (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_DDU - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          (coord[1] - 1 + ek_parameters_gpu->dim_y) % ek_parameters_gpu->dim_y,
          (coord[2] + 1) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_DUD - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          (coord[1] + 1) % ek_parameters_gpu->dim_y,
          (coord[2] - 1 + ek_parameters_gpu->dim_z) % ek_parameters_gpu->dim_z);

      neighborindex[EK_LINK_DUU - 13] = rhoindex_cartesian2linear(
          (coord[0] - 1 + ek_parameters_gpu->dim_x) % ek_parameters_gpu->dim_x,
          (coord[1] + 1) % ek_parameters_gpu->dim_y,
          (coord[2] + 1) % ek_parameters_gpu->dim_z);

      /* Clear fluxes on links connecting a boundary node */
      for (int i = 0; i < 13; i++)
        ek_parameters_gpu->j[jindex_getByRhoLinear(index, i)] = 0.0f;

      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_D00 - 13], EK_LINK_U00)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_0D0 - 13], EK_LINK_0U0)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_00D - 13], EK_LINK_00U)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_DD0 - 13], EK_LINK_UU0)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_DU0 - 13], EK_LINK_UD0)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_D0D - 13], EK_LINK_U0U)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_D0U - 13], EK_LINK_U0D)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_0DD - 13], EK_LINK_0UU)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_0DU - 13], EK_LINK_0UD)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_DDD - 13], EK_LINK_UUU)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_DDU - 13], EK_LINK_UUD)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_DUD - 13], EK_LINK_UDU)] = 0.0f;
      ek_parameters_gpu->j[jindex_getByRhoLinear(
          neighborindex[EK_LINK_DUU - 13], EK_LINK_UDD)] = 0.0f;
    }
  }
}

__global__ void ek_clear_fluxes() {
  unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    for (int i = 0; i < 13; i++) {
      ek_parameters_gpu->j[jindex_getByRhoLinear(index, i)] = 0.0f;
#ifdef EK_DEBUG
      ek_parameters_gpu->j_fluc[jindex_getByRhoLinear(index, i)] = 0.0f;
#endif
    }
  }
}

__global__ void ek_init_species_density_homogeneous() {
  unsigned int index = ek_getThreadIndex();
  unsigned int coord[3];

  rhoindex_linear2cartesian(index, coord);

  if (index < ek_parameters_gpu->number_of_nodes) {
    for (int i = 0; i < ek_parameters_gpu->number_of_species; i++) {
      ek_parameters_gpu->rho[i][index] =
          ek_parameters_gpu->density[i] * ek_parameters_gpu->agrid *
          ek_parameters_gpu->agrid * ek_parameters_gpu->agrid;
    }
  }
}

__global__ void ek_gather_species_charge_density() {
  auto const index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    ek_setNode(index, 0.0f);
    cufftReal tmp = 0.0f;
    for (int i = 0; i < ek_parameters_gpu->number_of_species; i++) {
      tmp += ek_parameters_gpu->valency[i] * ek_parameters_gpu->rho[i][index];
    }
    ek_setNode(index, tmp / powf(ek_parameters_gpu->agrid, 3));
  }
}

__global__ void
ek_gather_particle_charge_density(CUDA_particle_data *particle_data,
                                  size_t number_of_particles,
                                  LB_parameters_gpu *ek_lbparameters_gpu) {
  unsigned int index = ek_getThreadIndex();
  int lowernode[3];
  float cellpos[3];
  float gridpos;

  if (index < number_of_particles) {
    gridpos = particle_data[index].p[0] / ek_parameters_gpu->agrid - 0.5f;
    lowernode[0] = (int)floorf(gridpos);
    cellpos[0] = gridpos - static_cast<float>(lowernode[0]);

    gridpos = particle_data[index].p[1] / ek_parameters_gpu->agrid - 0.5f;
    lowernode[1] = (int)floorf(gridpos);
    cellpos[1] = gridpos - static_cast<float>(lowernode[1]);

    gridpos = particle_data[index].p[2] / ek_parameters_gpu->agrid - 0.5f;
    lowernode[2] = (int)floorf(gridpos);
    cellpos[2] = gridpos - static_cast<float>(lowernode[2]);

    lowernode[0] =
        static_cast<int>((lowernode[0] + ek_lbparameters_gpu->dim_x) %
                         ek_lbparameters_gpu->dim_x);
    lowernode[1] =
        static_cast<int>((lowernode[1] + ek_lbparameters_gpu->dim_y) %
                         ek_lbparameters_gpu->dim_y);
    lowernode[2] =
        static_cast<int>((lowernode[2] + ek_lbparameters_gpu->dim_z) %
                         ek_lbparameters_gpu->dim_z);

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  lowernode[0], lowernode[1], lowernode[2])],
              particle_data[index].q * (1 - cellpos[0]) * (1 - cellpos[1]) *
                  (1 - cellpos[2]));

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  (lowernode[0] + 1) % ek_parameters_gpu->dim_x, lowernode[1],
                  lowernode[2])],
              particle_data[index].q * cellpos[0] * (1 - cellpos[1]) *
                  (1 - cellpos[2]));

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  lowernode[0], (lowernode[1] + 1) % ek_parameters_gpu->dim_y,
                  lowernode[2])],
              particle_data[index].q * (1 - cellpos[0]) * cellpos[1] *
                  (1 - cellpos[2]));

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  lowernode[0], lowernode[1],
                  (lowernode[2] + 1) % ek_parameters_gpu->dim_z)],
              particle_data[index].q * (1 - cellpos[0]) * (1 - cellpos[1]) *
                  cellpos[2]);

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  (lowernode[0] + 1) % ek_parameters_gpu->dim_x,
                  (lowernode[1] + 1) % ek_parameters_gpu->dim_y, lowernode[2])],
              particle_data[index].q * cellpos[0] * cellpos[1] *
                  (1 - cellpos[2]));

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  (lowernode[0] + 1) % ek_parameters_gpu->dim_x, lowernode[1],
                  (lowernode[2] + 1) % ek_parameters_gpu->dim_z)],
              particle_data[index].q * cellpos[0] * (1 - cellpos[1]) *
                  cellpos[2]);

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  lowernode[0], (lowernode[1] + 1) % ek_parameters_gpu->dim_y,
                  (lowernode[2] + 1) % ek_parameters_gpu->dim_z)],
              particle_data[index].q * (1 - cellpos[0]) * cellpos[1] *
                  cellpos[2]);

    atomicAdd(&((cufftReal *)ek_parameters_gpu
                    ->charge_potential)[rhoindex_cartesian2linear_padded(
                  (lowernode[0] + 1) % ek_parameters_gpu->dim_x,
                  (lowernode[1] + 1) % ek_parameters_gpu->dim_y,
                  (lowernode[2] + 1) % ek_parameters_gpu->dim_z)],
              particle_data[index].q * cellpos[0] * cellpos[1] * cellpos[2]);
  }
}

__global__ void
ek_spread_particle_force(CUDA_particle_data *particle_data,
                         size_t number_of_particles, float *particle_forces,
                         LB_parameters_gpu *ek_lbparameters_gpu) {

  unsigned int index = ek_getThreadIndex();
  int lowernode[3];
  float cellpos[3];
  float gridpos;

  if (index < number_of_particles) {
    gridpos = particle_data[index].p[0] / ek_parameters_gpu->agrid - 0.5f;
    lowernode[0] = (int)floorf(gridpos);
    cellpos[0] = gridpos - (float)(lowernode[0]);

    gridpos = particle_data[index].p[1] / ek_parameters_gpu->agrid - 0.5f;
    lowernode[1] = (int)floorf(gridpos);
    cellpos[1] = gridpos - (float)(lowernode[1]);

    gridpos = particle_data[index].p[2] / ek_parameters_gpu->agrid - 0.5f;
    lowernode[2] = (int)floorf(gridpos);
    cellpos[2] = gridpos - (float)(lowernode[2]);

    lowernode[0] =
        static_cast<int>((lowernode[0] + ek_lbparameters_gpu->dim_x) %
                         ek_lbparameters_gpu->dim_x);
    lowernode[1] =
        static_cast<int>((lowernode[1] + ek_lbparameters_gpu->dim_y) %
                         ek_lbparameters_gpu->dim_y);
    lowernode[2] =
        static_cast<int>((lowernode[2] + ek_lbparameters_gpu->dim_z) %
                         ek_lbparameters_gpu->dim_z);

    float efield[3] = {0., 0., 0.};
    for (unsigned int dim = 0; dim < 3; ++dim) {
      // 0 0 0
      efield[dim] +=
          ek_parameters_gpu->electric_field[3 * rhoindex_cartesian2linear(
                                                    lowernode[0], lowernode[1],
                                                    lowernode[2]) +
                                            dim] *
          (1 - cellpos[0]) * (1 - cellpos[1]) * (1 - cellpos[2]);

      // 0 0 1
      efield[dim] +=
          ek_parameters_gpu
              ->electric_field[3 * rhoindex_cartesian2linear(
                                       lowernode[0], lowernode[1],
                                       (lowernode[2] + 1) %
                                           ek_lbparameters_gpu->dim_z) +
                               dim] *
          (1 - cellpos[0]) * (1 - cellpos[1]) * cellpos[2];

      // 0 1 0
      efield[dim] +=
          ek_parameters_gpu
              ->electric_field[3 * rhoindex_cartesian2linear(
                                       lowernode[0],
                                       (lowernode[1] + 1) %
                                           ek_lbparameters_gpu->dim_y,
                                       lowernode[2]) +
                               dim] *
          (1 - cellpos[0]) * cellpos[1] * (1 - cellpos[2]);

      // 0 1 1
      efield[dim] +=
          ek_parameters_gpu->electric_field
              [3 * rhoindex_cartesian2linear(
                       lowernode[0],
                       (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                       (lowernode[2] + 1) % ek_lbparameters_gpu->dim_z) +
               dim] *
          (1 - cellpos[0]) * cellpos[1] * cellpos[2];

      // 1 0 0
      efield[dim] +=
          ek_parameters_gpu
              ->electric_field[3 * rhoindex_cartesian2linear(
                                       (lowernode[0] + 1) %
                                           ek_lbparameters_gpu->dim_x,
                                       lowernode[1], lowernode[2]) +
                               dim] *
          cellpos[0] * (1 - cellpos[1]) * (1 - cellpos[2]);

      // 1 0 1
      efield[dim] +=
          ek_parameters_gpu->electric_field
              [3 * rhoindex_cartesian2linear(
                       (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                       lowernode[1],
                       (lowernode[2] + 1) % ek_lbparameters_gpu->dim_z) +
               dim] *
          cellpos[0] * (1 - cellpos[1]) * cellpos[2];

      // 1 1 0
      efield[dim] +=
          ek_parameters_gpu->electric_field
              [3 * rhoindex_cartesian2linear(
                       (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                       (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                       lowernode[2]) +
               dim] *
          cellpos[0] * cellpos[1] * (1 - cellpos[2]);

      // 1 1 1
      efield[dim] +=
          ek_parameters_gpu->electric_field
              [3 * rhoindex_cartesian2linear(
                       (lowernode[0] + 1) % ek_lbparameters_gpu->dim_x,
                       (lowernode[1] + 1) % ek_lbparameters_gpu->dim_y,
                       (lowernode[2] + 1) % ek_lbparameters_gpu->dim_z) +
               dim] *
          cellpos[0] * cellpos[1] * cellpos[2];
    }
    particle_forces[3 * index + 0] += particle_data[index].q * efield[0];
    particle_forces[3 * index + 1] += particle_data[index].q * efield[1];
    particle_forces[3 * index + 2] += particle_data[index].q * efield[2];
  }
}

__global__ void ek_calc_electric_field(const float *potential) {
  unsigned int coord[3];
  const unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    rhoindex_linear2cartesian(index, coord);
    const float agrid_inv = 1.0f / ek_parameters_gpu->agrid;

    ek_parameters_gpu->electric_field[3 * index + 0] =
        -0.5f * agrid_inv *
        (potential[rhoindex_cartesian2linear_padded(
             (coord[0] + 1) % ek_parameters_gpu->dim_x, coord[1], coord[2])] -
         potential[rhoindex_cartesian2linear_padded(
             (coord[0] - 1 + ek_parameters_gpu->dim_x) %
                 ek_parameters_gpu->dim_x,
             coord[1], coord[2])]);
    ek_parameters_gpu->electric_field[3 * index + 1] =
        -0.5f * agrid_inv *
        (potential[rhoindex_cartesian2linear_padded(
             coord[0], (coord[1] + 1) % ek_parameters_gpu->dim_y, coord[2])] -
         potential[rhoindex_cartesian2linear_padded(
             coord[0],
             (coord[1] - 1 + ek_parameters_gpu->dim_y) %
                 ek_parameters_gpu->dim_y,
             coord[2])]);
    ek_parameters_gpu->electric_field[3 * index + 2] =
        -0.5f * agrid_inv *
        (potential[rhoindex_cartesian2linear_padded(
             coord[0], coord[1], (coord[2] + 1) % ek_parameters_gpu->dim_z)] -
         potential[rhoindex_cartesian2linear_padded(
             coord[0], coord[1],
             (coord[2] - 1 + ek_parameters_gpu->dim_z) %
                 ek_parameters_gpu->dim_z)]);
  }
}

__global__ void ek_clear_boundary_densities(LB_nodes_gpu lbnode) {

  unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    if (lbnode.boundary[index]) {
      for (int i = 0; i < ek_parameters_gpu->number_of_species; i++) {
        ek_parameters_gpu->rho[i][index] = 0.0f;
      }
    }
  }
}

__global__ void ek_calculate_system_charge(ekfloat *charge_gpu) {

  unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    for (int i = 0; i < ek_parameters_gpu->number_of_species; i++) {
      atomicAdd(charge_gpu, ek_parameters_gpu->rho[i][index] *
                                ek_parameters_gpu->valency[i]);
    }
  }
}

// TODO delete ?? (it has the previous step setting now)
// This is not compatible with external LB force_densitys!
__global__ void ek_clear_node_force(LB_node_force_density_gpu node_f) {

  unsigned int index = ek_getThreadIndex();

  if (index < ek_parameters_gpu->number_of_nodes) {
    ek_parameters_gpu->lb_force_density_previous[index] =
        node_f.force_density[index];
    ek_parameters_gpu
        ->lb_force_density_previous[ek_parameters_gpu->number_of_nodes +
                                    index] =
        node_f.force_density[ek_parameters_gpu->number_of_nodes + index];
    ek_parameters_gpu
        ->lb_force_density_previous[2 * ek_parameters_gpu->number_of_nodes +
                                    index] =
        node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index];

    node_f.force_density[index] = 0.0f;
    node_f.force_density[ek_parameters_gpu->number_of_nodes + index] = 0.0f;
    node_f.force_density[2 * ek_parameters_gpu->number_of_nodes + index] = 0.0f;
  }
}

void ek_calculate_electrostatic_coupling() {
  const int blocks_per_grid_y = 4;
  const int threads_per_block = 64;

  if ((!ek_parameters.es_coupling) || (!ek_initialized))
    return;

  auto device_particles = gpu_get_particle_pointer();
  auto blocks_per_grid_x = static_cast<int>(
      (device_particles.size() + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_spread_particle_force, dim_grid, threads_per_block,
             device_particles.data(), device_particles.size(),
             gpu_get_particle_force_pointer(), ek_lbparameters_gpu);
}

void ek_integrate_electrostatics() {

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_gather_species_charge_density, dim_grid, threads_per_block);

  if (ek_parameters.es_coupling) {
    cuda_safe_mem(cudaMemcpy(
        ek_parameters.charge_potential_buffer, ek_parameters.charge_potential,
        sizeof(cufftComplex) * ek_parameters.dim_z * ek_parameters.dim_y *
            (ek_parameters.dim_x / 2 + 1),
        cudaMemcpyDeviceToDevice));
    electrostatics->calculatePotential(
        (cufftComplex *)ek_parameters.charge_potential_buffer);
    KERNELCALL(ek_calc_electric_field, dim_grid, threads_per_block,
               ek_parameters.charge_potential_buffer);
  }

  auto device_particles = gpu_get_particle_pointer();
  if (not device_particles
              .empty()) // TODO make it an if number_of_charged_particles != 0
  {
    blocks_per_grid_x = static_cast<int>(
        (device_particles.size() + threads_per_block * blocks_per_grid_y - 1) /
        (threads_per_block * blocks_per_grid_y));
    dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

    particle_data_gpu = device_particles.data();

    KERNELCALL(ek_gather_particle_charge_density, dim_grid, threads_per_block,
               particle_data_gpu, device_particles.size(), ek_lbparameters_gpu);
  }

  electrostatics->calculatePotential();
}

void ek_integrate() {
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  /* Clears the force on the nodes and must be called before fluxes are
     calculated, since in the reaction set up the previous-step LB force is
     added to the flux
     (in ek_calculate_quantities / ek_displacement), which is copied in this
     routine */

  // KERNELCALL( ek_clear_node_force, dim_grid, threads_per_block, node_f );

  /* Integrate diffusion-advection */
  for (int i = 0; i < ek_parameters.number_of_species; i++) {
    KERNELCALL(ek_clear_fluxes, dim_grid, threads_per_block);
    KERNELCALL(ek_calculate_quantities, dim_grid, threads_per_block, i,
               *current_nodes, node_f, ek_lbparameters_gpu, ek_lb_device_values,
               philox_counter.value());

    KERNELCALL(ek_propagate_densities, dim_grid, threads_per_block, i);
  }

  /* Integrate electrostatics */
  ek_integrate_electrostatics();

  /* Integrate Navier-Stokes */
  lb_integrate_GPU();

  philox_counter.increment();
}

#ifdef EK_BOUNDARIES
void ek_gather_wallcharge_species_density(ekfloat *wallcharge_species_density,
                                          int wallcharge_species) {
  if (wallcharge_species != -1) {
    cuda_safe_mem(cudaMemcpy(wallcharge_species_density,
                             ek_parameters.rho[wallcharge_species],
                             ek_parameters.number_of_nodes * sizeof(ekfloat),
                             cudaMemcpyDeviceToHost));
  }
}
void ek_init_species_density_wallcharge(ekfloat *wallcharge_species_density,
                                        int wallcharge_species) {
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_clear_boundary_densities, dim_grid, threads_per_block,
             *current_nodes);

  if (wallcharge_species != -1) {
    cuda_safe_mem(cudaMemcpy(ek_parameters.rho[wallcharge_species],
                             wallcharge_species_density,
                             ek_parameters.number_of_nodes * sizeof(ekfloat),
                             cudaMemcpyHostToDevice));
  }
}
#endif

void ek_init_species(int species) {
  if (!ek_initialized) {
    ek_init();
  }

  if (ek_parameters.species_index[species] == -1) {
    ek_parameters.species_index[species] =
        static_cast<int>(ek_parameters.number_of_species);
    ek_parameters.number_of_species++;

    cuda_safe_mem(cudaMalloc(
        (void **)&ek_parameters.rho[ek_parameters.species_index[species]],
        ek_parameters.number_of_nodes * sizeof(ekfloat)));

    ek_parameters.density[ek_parameters.species_index[species]] = 0.0;
    ek_parameters.D[ek_parameters.species_index[species]] = 0.0;
    ek_parameters.valency[ek_parameters.species_index[species]] = 0.0;
    ek_parameters.ext_force_density[0][ek_parameters.species_index[species]] =
        0.0;
    ek_parameters.ext_force_density[1][ek_parameters.species_index[species]] =
        0.0;
    ek_parameters.ext_force_density[2][ek_parameters.species_index[species]] =
        0.0;
    ek_parameters.d[ek_parameters.species_index[species]] =
        ek_parameters.D[ek_parameters.species_index[species]] /
        (1.0f + 2.0f * sqrt(2.0f));
  }
}

int ek_init() {
  if (ek_parameters.agrid < 0.0 || ek_parameters.viscosity < 0.0 ||
      ek_parameters.T < 0.0 || ek_parameters.prefactor < 0.0) {

    fprintf(stderr, "ERROR: invalid agrid, viscosity, T or prefactor\n");

    return 1;
  }

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x;
  dim3 dim_grid;

  if (!ek_initialized) {
    if (cudaGetSymbolAddress((void **)&ek_parameters_gpu_pointer,
                             HIP_SYMBOL(ek_parameters_gpu)) != cudaSuccess) {
      fprintf(stderr, "ERROR: Fetching constant memory pointer\n");

      return 1;
    }

    for (auto &val : ek_parameters.species_index) {
      val = -1;
    }

    if (lattice_switch != ActiveLB::NONE) {
      fprintf(stderr,
              "ERROR: Electrokinetics automatically initializes the LB on the "
              "GPU and can therefore not be used in conjunction with LB.\n");
      fprintf(stderr, "ERROR: Please run either electrokinetics or LB.\n");

      return 1;
    }

    lattice_switch = ActiveLB::GPU;
    ek_initialized = true;

    lbpar_gpu.agrid = ek_parameters.agrid;
    lbpar_gpu.viscosity = 1.0;      // dummy values (real initialization later)
    lbpar_gpu.bulk_viscosity = 1.0; // dummy values (real initialization later)
    lb_lbcoupling_set_gamma(ek_parameters.friction);

    // Convert the density (given in MD units) to LB units
    lbpar_gpu.rho = (ek_parameters.lb_density < 0.0
                         ? 1.0f
                         : ek_parameters.lb_density * ek_parameters.agrid *
                               ek_parameters.agrid * ek_parameters.agrid);

    lbpar_gpu.is_TRT = true;

    lb_reinit_parameters_gpu();
    lbpar_gpu.viscosity = ek_parameters.viscosity * lbpar_gpu.time_step /
                          (ek_parameters.agrid * ek_parameters.agrid);
    lbpar_gpu.bulk_viscosity = ek_parameters.bulk_viscosity *
                               lbpar_gpu.time_step /
                               (ek_parameters.agrid * ek_parameters.agrid);
    lb_reinit_parameters_gpu();

    lb_init_gpu();

    if (ek_parameters.lb_force_density[0] != 0 ||
        ek_parameters.lb_force_density[1] != 0 ||
        ek_parameters.lb_force_density[2] != 0) {
      lbpar_gpu.external_force_density = 1;
      lbpar_gpu.ext_force_density[0] =
          ek_parameters.lb_force_density[0] * ek_parameters.agrid *
          ek_parameters.agrid * ek_parameters.time_step *
          ek_parameters.time_step;
      lbpar_gpu.ext_force_density[1] =
          ek_parameters.lb_force_density[1] * ek_parameters.agrid *
          ek_parameters.agrid * ek_parameters.time_step *
          ek_parameters.time_step;
      lbpar_gpu.ext_force_density[2] =
          ek_parameters.lb_force_density[2] * ek_parameters.agrid *
          ek_parameters.agrid * ek_parameters.time_step *
          ek_parameters.time_step;
      lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);
    } else {
      lbpar_gpu.external_force_density = 0;
      lbpar_gpu.ext_force_density[0] = 0;
      lbpar_gpu.ext_force_density[1] = 0;
      lbpar_gpu.ext_force_density[2] = 0;
    }

    ek_parameters.dim_x = lbpar_gpu.dim_x;
    ek_parameters.dim_x_padded = (ek_parameters.dim_x / 2 + 1) * 2;
    ek_parameters.dim_y = lbpar_gpu.dim_y;
    ek_parameters.dim_z = lbpar_gpu.dim_z;
    ek_parameters.time_step = lbpar_gpu.time_step;
    ek_parameters.number_of_nodes =
        ek_parameters.dim_x * ek_parameters.dim_y * ek_parameters.dim_z;

    cuda_safe_mem(
        cudaMalloc((void **)&ek_parameters.j,
                   ek_parameters.number_of_nodes * 13 * sizeof(ekfloat)));
#ifdef EK_DEBUG
    cuda_safe_mem(
        cudaMalloc((void **)&ek_parameters.j_fluc,
                   ek_parameters.number_of_nodes * 13 * sizeof(ekfloat)));
#endif

    cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(ek_parameters_gpu),
                                     &ek_parameters, sizeof(EK_parameters)));

    lb_get_para_pointer(&ek_lbparameters_gpu);
    lb_set_ek_pointer(ek_parameters_gpu_pointer);

    cuda_safe_mem(
        cudaMalloc((void **)&ek_parameters.lb_force_density_previous,
                   ek_parameters.number_of_nodes * 3 * sizeof(float)));

    if (ek_parameters.es_coupling) {
      cuda_safe_mem(cudaMalloc((void **)&ek_parameters.charge_potential_buffer,
                               sizeof(cufftComplex) * ek_parameters.dim_z *
                                   ek_parameters.dim_y *
                                   (ek_parameters.dim_x / 2 + 1)));
      cuda_safe_mem(
          cudaMalloc((void **)&ek_parameters.electric_field,
                     ek_parameters.number_of_nodes * 3 * sizeof(float)));
    }

    cuda_safe_mem(cudaMalloc((void **)&charge_gpu, sizeof(ekfloat)));

    lb_get_device_values_pointer(&ek_lb_device_values);

    if (cudaGetLastError() != cudaSuccess) {
      fprintf(stderr, "ERROR: Failed to allocate\n");
      return 1;
    }

    cudaMallocHost((void **)&ek_parameters.node_is_catalyst,
                   sizeof(char) * ek_parameters.dim_z * ek_parameters.dim_y *
                       ek_parameters.dim_x);

    if (cudaGetLastError() != cudaSuccess) {
      fprintf(stderr, "ERROR: Failed to allocate\n");
      return 1;
    }

    // initialize electrostatics
    delete electrostatics;

    FdElectrostatics::InputParameters es_parameters = {
        ek_parameters.prefactor, int(ek_parameters.dim_x),
        int(ek_parameters.dim_y), int(ek_parameters.dim_z),
        ek_parameters.agrid};
    try {
      electrostatics = new FdElectrostatics(es_parameters, stream[0]);
    } catch (std::string e) {
      std::cerr << "Error in initialization of electrokinetics electrostatics "
                   "solver: "
                << e << std::endl;
      return 1;
    }

    ek_parameters.charge_potential = electrostatics->getGrid().grid;
    cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(ek_parameters_gpu),
                                     &ek_parameters, sizeof(EK_parameters)));

    // clear initial LB force and finish up
    blocks_per_grid_x = static_cast<int>(
        (ek_parameters.dim_z * ek_parameters.dim_y * (ek_parameters.dim_x) +
         threads_per_block * blocks_per_grid_y - 1) /
        (threads_per_block * blocks_per_grid_y));
    dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
    KERNELCALL(ek_clear_node_force, dim_grid, threads_per_block, node_f);

    ek_initialized = true;
  } else {
    if (lbpar_gpu.agrid != ek_parameters.agrid ||
        lbpar_gpu.viscosity !=
            ek_parameters.viscosity * ek_parameters.time_step /
                (ek_parameters.agrid * ek_parameters.agrid) ||
        lbpar_gpu.bulk_viscosity !=
            ek_parameters.bulk_viscosity * ek_parameters.time_step /
                (ek_parameters.agrid * ek_parameters.agrid) ||
        lb_lbcoupling_get_gamma() != ek_parameters.friction ||
        lbpar_gpu.rho != ek_parameters.lb_density * ek_parameters.agrid *
                             ek_parameters.agrid * ek_parameters.agrid) {
      fprintf(stderr,
              "ERROR: The LB parameters on the GPU cannot be reinitialized.\n");

      return 1;
    }
    cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(ek_parameters_gpu),
                                     &ek_parameters, sizeof(EK_parameters)));

    blocks_per_grid_x =
        static_cast<int>((ek_parameters.number_of_nodes +
                          threads_per_block * blocks_per_grid_y - 1) /
                         (threads_per_block * blocks_per_grid_y));
    dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

    KERNELCALL(ek_init_species_density_homogeneous, dim_grid,
               threads_per_block);

#ifdef EK_BOUNDARIES
    LBBoundaries::lb_init_boundaries();
    lb_get_boundary_force_pointer(&ek_lb_boundary_force);

    cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(ek_parameters_gpu),
                                     &ek_parameters, sizeof(EK_parameters)));
#endif

    ek_integrate_electrostatics();
  }
  return 0;
}

void lb_set_ek_pointer(EK_parameters *pointeradress) {
  lb_ek_parameters_gpu = pointeradress;
}

unsigned int ek_calculate_boundary_mass() {
  auto *bound_array = (unsigned int *)Utils::malloc(lbpar_gpu.number_of_nodes *
                                                    sizeof(unsigned int));

  lb_get_boundary_flags_GPU(bound_array);

  unsigned int boundary_node_number = 0;

  for (int j = 0; j < ek_parameters.number_of_nodes; j++)
    if (bound_array[j] != 0)
      boundary_node_number++;

  free(bound_array);

  return boundary_node_number;
}

void rhoindex_linear2cartesian_host(unsigned int index, unsigned int *coord) {

  coord[0] = index % ek_parameters.dim_x;
  index /= ek_parameters.dim_x;
  coord[1] = index % ek_parameters.dim_y;
  coord[2] = index / ek_parameters.dim_y;
}

unsigned int jindex_cartesian2linear_host(unsigned int x, unsigned int y,
                                          unsigned int z, unsigned int c) {

  x = (x + ek_parameters.dim_x) %
      ek_parameters
          .dim_x; // this does not happen in the GPU version of this function
  y = (y + ek_parameters.dim_y) % ek_parameters.dim_y;
  z = (z + ek_parameters.dim_z) % ek_parameters.dim_z;

  return c * ek_parameters.number_of_nodes +
         z * ek_parameters.dim_y * ek_parameters.dim_x +
         y * ek_parameters.dim_x + x;
}

unsigned int jindex_getByRhoLinear_host(unsigned int rho_index,
                                        unsigned int c) {

  return c * ek_parameters.number_of_nodes + rho_index;
}

unsigned int rhoindex_cartesian2linear_host(unsigned int x, unsigned int y,
                                            unsigned int z) {

  return z * ek_parameters.dim_y * ek_parameters.dim_x +
         y * ek_parameters.dim_x + x;
}

int ek_lb_print_vtk_velocity(char *filename) {

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  auto *host_values = (LB_rho_v_pi_gpu *)Utils::malloc(
      lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu));
  lb_get_values_GPU(host_values);
  auto const lattice_speed = lbpar_gpu.agrid / lbpar_gpu.tau;
  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
velocity\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\nPOINT_DATA %u\n\
SCALARS velocity float 3\n\
LOOKUP_TABLE default\n",
          lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z,
          lbpar_gpu.agrid * 0.5f, lbpar_gpu.agrid * 0.5f,
          lbpar_gpu.agrid * 0.5f, lbpar_gpu.agrid, lbpar_gpu.agrid,
          lbpar_gpu.agrid, lbpar_gpu.number_of_nodes);

  for (int i = 0; i < lbpar_gpu.number_of_nodes; i++) {
    fprintf(fp, "%e %e %e ", host_values[i].v[0] * lattice_speed,
            host_values[i].v[1] * lattice_speed,
            host_values[i].v[2] * lattice_speed);
  }

  free(host_values);
  fclose(fp);

  return 0;
}

int ek_node_print_velocity(
    int x, int y, int z,
    double *velocity) { // TODO only calculate single node velocity

  auto *host_values = (LB_rho_v_pi_gpu *)Utils::malloc(
      lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu));
  lb_get_values_GPU(host_values);

  auto const i = z * ek_parameters.dim_y * ek_parameters.dim_x +
                 y * ek_parameters.dim_x + x;
  auto const lattice_speed = lbpar_gpu.agrid / lbpar_gpu.tau;

  velocity[0] = host_values[i].v[0] * lattice_speed;
  velocity[1] = host_values[i].v[1] * lattice_speed;
  velocity[2] = host_values[i].v[2] * lattice_speed;

  free(host_values);

  return 0;
}

int ek_lb_print_vtk_density(char *filename) {

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  auto *host_values = (LB_rho_v_pi_gpu *)Utils::malloc(
      lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu));
  lb_get_values_GPU(host_values);

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
density_lb\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS density_lb float 1\n\
LOOKUP_TABLE default\n",
          lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z,
          lbpar_gpu.agrid * 0.5f, lbpar_gpu.agrid * 0.5f,
          lbpar_gpu.agrid * 0.5f, lbpar_gpu.agrid, lbpar_gpu.agrid,
          lbpar_gpu.agrid, lbpar_gpu.number_of_nodes);
  auto const agrid = lb_lbfluid_get_agrid();
  for (int i = 0; i < lbpar_gpu.number_of_nodes; i++) {
    fprintf(fp, "%e ", host_values[i].rho / agrid / agrid / agrid);
  }

  free(host_values);
  fclose(fp);

  return 0;
}

int ek_print_vtk_density(int species, char *filename) {

  if (ek_parameters.species_index[species] == -1) {
    return 1;
  }

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  auto *densities =
      (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * sizeof(ekfloat));

  cuda_safe_mem(cudaMemcpy(
      densities, ek_parameters.rho[ek_parameters.species_index[species]],
      ek_parameters.number_of_nodes * sizeof(ekfloat), cudaMemcpyDeviceToHost));

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
density_%d\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS density_%d float 1\n\
LOOKUP_TABLE default\n",
          species, ek_parameters.dim_x, ek_parameters.dim_y,
          ek_parameters.dim_z, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes, species);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {
    fprintf(fp, "%e\n",
            densities[i] / (ek_parameters.agrid * ek_parameters.agrid *
                            ek_parameters.agrid));
  }

  free(densities);
  fclose(fp);

  return 0;
}

int ek_node_print_density(int species, int x, int y, int z, double *density) {

  if (ek_parameters.species_index[species] == -1) {
    return 1;
  }

  auto *densities =
      (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * sizeof(ekfloat));

  cuda_safe_mem(cudaMemcpy(
      densities, ek_parameters.rho[ek_parameters.species_index[species]],
      ek_parameters.number_of_nodes * sizeof(ekfloat), cudaMemcpyDeviceToHost));

  *density = densities[z * ek_parameters.dim_y * ek_parameters.dim_x +
                       y * ek_parameters.dim_x + x] /
             (ek_parameters.agrid * ek_parameters.agrid * ek_parameters.agrid);

  free(densities);

  return 0;
}

int ek_node_print_flux(int species, int x, int y, int z, double *flux) {

  if (ek_parameters.species_index[species] == -1) {
    return 1;
  }

  ekfloat flux_local_cartesian[3]; // temporary variable for converting fluxes
                                   // into Cartesian coordinates for output
  unsigned int coord[3];

  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  auto *fluxes = (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * 13 *
                                          sizeof(ekfloat));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_clear_fluxes, dim_grid, threads_per_block);
  KERNELCALL(ek_calculate_quantities, dim_grid, threads_per_block,
             ek_parameters.species_index[species], *current_nodes, node_f,
             ek_lbparameters_gpu, ek_lb_device_values, philox_counter.value());
  reset_LB_force_densities_GPU(false);

#ifdef EK_BOUNDARIES
  KERNELCALL(ek_apply_boundaries, dim_grid, threads_per_block,
             ek_parameters.species_index[species], *current_nodes, node_f);
#endif

  cuda_safe_mem(cudaMemcpy(fluxes, ek_parameters.j,
                           ek_parameters.number_of_nodes * 13 * sizeof(ekfloat),
                           cudaMemcpyDeviceToHost));

  auto const i = rhoindex_cartesian2linear_host(coord[0], coord[1], coord[2]);

  flux_local_cartesian[0] =
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U00)];

  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UU0)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UD0)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0U)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0D)];

  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1],
                                                coord[2], EK_LINK_D00 - 13)];

  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] - 1,
                                                coord[2], EK_LINK_DD0 - 13)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] + 1,
                                                coord[2], EK_LINK_DU0 - 13)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1], coord[2] - 1, EK_LINK_D0D - 13)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1], coord[2] + 1, EK_LINK_D0U - 13)];

  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
  flux_local_cartesian[0] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

  flux_local_cartesian[1] =
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0U0)];

  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UU0)];
  flux_local_cartesian[1] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UD0)];
  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UU)];
  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UD)];

  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
  flux_local_cartesian[1] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
  flux_local_cartesian[1] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_cartesian2linear_host(coord[0], coord[1] - 1,
                                                coord[2], EK_LINK_0D0 - 13)];

  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] - 1,
                                                coord[2], EK_LINK_DD0 - 13)];
  flux_local_cartesian[1] -=
      0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] + 1,
                                                coord[2], EK_LINK_DU0 - 13)];
  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0], coord[1] - 1, coord[2] - 1, EK_LINK_0DD - 13)];
  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0], coord[1] - 1, coord[2] + 1, EK_LINK_0DU - 13)];

  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
  flux_local_cartesian[1] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
  flux_local_cartesian[1] -=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
  flux_local_cartesian[1] -=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

  flux_local_cartesian[2] =
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_00U)];

  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0U)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0D)];
  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UU)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UD)];

  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0], coord[1], coord[2] - 1, EK_LINK_00D - 13)];

  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1], coord[2] - 1, EK_LINK_D0D - 13)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1], coord[2] + 1, EK_LINK_D0U - 13)];
  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0], coord[1] - 1, coord[2] - 1, EK_LINK_0DD - 13)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0], coord[1] - 1, coord[2] + 1, EK_LINK_0DU - 13)];

  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
  flux_local_cartesian[2] +=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
  flux_local_cartesian[2] -=
      0.5 * fluxes[jindex_cartesian2linear_host(
                coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

  flux[0] =
      flux_local_cartesian[0] /
      (ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid);
  flux[1] =
      flux_local_cartesian[1] /
      (ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid);
  flux[2] =
      flux_local_cartesian[2] /
      (ek_parameters.time_step * ek_parameters.agrid * ek_parameters.agrid);

  free(fluxes);

  return 0;
}

int ek_node_set_density(int species, int x, int y, int z, double density) {
  if (ek_parameters.species_index[species] != -1) {
    auto index =
        static_cast<int>(z * ek_parameters.dim_y * ek_parameters.dim_x +
                         y * ek_parameters.dim_x + x);
    ekfloat num_particles = density * ek_parameters.agrid *
                            ek_parameters.agrid * ek_parameters.agrid;

    cuda_safe_mem(cudaMemcpy(
        &ek_parameters.rho[ek_parameters.species_index[species]][index],
        &num_particles, sizeof(ekfloat), cudaMemcpyHostToDevice));
  } else
    return 1;

  return 0;
}

int ek_print_vtk_flux(int species, char *filename) {

  if (ek_parameters.species_index[species] == -1) {
    return 1;
  }

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  ekfloat flux_local_cartesian[3]; // temporary variable for converting fluxes
                                   // into Cartesian coordinates for output

  unsigned int coord[3];

  auto *fluxes = (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * 13 *
                                          sizeof(ekfloat));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_clear_fluxes, dim_grid, threads_per_block);
  KERNELCALL(ek_calculate_quantities, dim_grid, threads_per_block,
             ek_parameters.species_index[species], *current_nodes, node_f,
             ek_lbparameters_gpu, ek_lb_device_values, philox_counter.value());
  reset_LB_force_densities_GPU(false);

#ifdef EK_BOUNDARIES
  KERNELCALL(ek_apply_boundaries, dim_grid, threads_per_block,
             ek_parameters.species_index[species], *current_nodes, node_f);
#endif

  cuda_safe_mem(cudaMemcpy(fluxes, ek_parameters.j,
                           ek_parameters.number_of_nodes * 13 * sizeof(ekfloat),
                           cudaMemcpyDeviceToHost));

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
flux_%d\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS flux_%d float 3\n\
LOOKUP_TABLE default\n",
          species, ek_parameters.dim_x, ek_parameters.dim_y,
          ek_parameters.dim_z, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes, species);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {
    rhoindex_linear2cartesian_host(i, coord);

    flux_local_cartesian[0] =
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U00)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UU0)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UD0)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0U)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0D)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1],
                                                  coord[2], EK_LINK_D00 - 13)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] - 1,
                                                  coord[2], EK_LINK_DD0 - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] + 1,
                                                  coord[2], EK_LINK_DU0 - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] - 1, EK_LINK_D0D - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] + 1, EK_LINK_D0U - 13)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

    flux_local_cartesian[1] =
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0U0)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UU0)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UD0)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UU)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UD)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0], coord[1] - 1,
                                                  coord[2], EK_LINK_0D0 - 13)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] - 1,
                                                  coord[2], EK_LINK_DD0 - 13)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] + 1,
                                                  coord[2], EK_LINK_DU0 - 13)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] - 1, EK_LINK_0DD - 13)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] + 1, EK_LINK_0DU - 13)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

    flux_local_cartesian[2] =
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_00U)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0U)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0D)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UU)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UD)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1], coord[2] - 1, EK_LINK_00D - 13)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] - 1, EK_LINK_D0D - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] + 1, EK_LINK_D0U - 13)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] - 1, EK_LINK_0DD - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] + 1, EK_LINK_0DU - 13)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

    fprintf(
        fp, "%e %e %e\n",
        flux_local_cartesian[0] / (ek_parameters.time_step *
                                   ek_parameters.agrid * ek_parameters.agrid),
        flux_local_cartesian[1] / (ek_parameters.time_step *
                                   ek_parameters.agrid * ek_parameters.agrid),
        flux_local_cartesian[2] / (ek_parameters.time_step *
                                   ek_parameters.agrid * ek_parameters.agrid));
  }

  free(fluxes);
  fclose(fp);

  return 0;
}

int ek_print_vtk_flux_fluc(int species, char *filename) {
#ifndef EK_DEBUG
  return 1;
#else
  FILE *fp = fopen(filename, "w");
  ekfloat flux_local_cartesian[3]; // temporary variable for converting fluxes
                                   // into cartesian coordinates for output

  unsigned int coord[3];

  if (fp == nullptr) {
    return 1;
  }

  ekfloat *fluxes = (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes *
                                             13 * sizeof(ekfloat));

  if (ek_parameters.species_index[species] != -1) {
    int threads_per_block = 64;
    int blocks_per_grid_y = 4;
    int blocks_per_grid_x = (ek_parameters.number_of_nodes +
                             threads_per_block * blocks_per_grid_y - 1) /
                            (threads_per_block * blocks_per_grid_y);
    dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

    KERNELCALL(ek_clear_fluxes, dim_grid, threads_per_block);
    KERNELCALL(ek_calculate_quantities, dim_grid, threads_per_block,
               ek_parameters.species_index[species], *current_nodes, node_f,
               ek_lbparameters_gpu, ek_lb_device_values,
               philox_counter.value());
    reset_LB_force_densities_GPU(false);

#ifdef EK_BOUNDARIES
    KERNELCALL(ek_apply_boundaries, dim_grid, threads_per_block,
               ek_parameters.species_index[species], *current_nodes, node_f);
#endif

    cuda_safe_mem(
        cudaMemcpy(fluxes, ek_parameters.j_fluc,
                   ek_parameters.number_of_nodes * 13 * sizeof(ekfloat),
                   cudaMemcpyDeviceToHost));
  } else
    return 1;

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
flux_%d\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS flux_%d float 3\n\
LOOKUP_TABLE default\n",
          species, ek_parameters.dim_x, ek_parameters.dim_y,
          ek_parameters.dim_z, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes, species);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {

    float flux_local_linksum = 0;
    rhoindex_linear2cartesian_host(i, coord);

    flux_local_cartesian[0] =
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U00)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UU0)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UD0)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0U)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0D)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1],
                                                  coord[2], EK_LINK_D00 - 13)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] - 1,
                                                  coord[2], EK_LINK_DD0 - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] + 1,
                                                  coord[2], EK_LINK_DU0 - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] - 1, EK_LINK_D0D - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] + 1, EK_LINK_D0U - 13)];

    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
    flux_local_cartesian[0] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

    flux_local_cartesian[1] =
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0U0)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UU0)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UD0)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UU)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UD)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0], coord[1] - 1,
                                                  coord[2], EK_LINK_0D0 - 13)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] - 1,
                                                  coord[2], EK_LINK_DD0 - 13)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_cartesian2linear_host(coord[0] - 1, coord[1] + 1,
                                                  coord[2], EK_LINK_DU0 - 13)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] - 1, EK_LINK_0DD - 13)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] + 1, EK_LINK_0DU - 13)];

    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
    flux_local_cartesian[1] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
    flux_local_cartesian[1] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

    flux_local_cartesian[2] =
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_00U)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0U)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_U0D)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UU)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_0UD)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUU)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UUD)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDU)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_getByRhoLinear_host(i, EK_LINK_UDD)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1], coord[2] - 1, EK_LINK_00D - 13)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] - 1, EK_LINK_D0D - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1], coord[2] + 1, EK_LINK_D0U - 13)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] - 1, EK_LINK_0DD - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0], coord[1] - 1, coord[2] + 1, EK_LINK_0DU - 13)];

    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] - 1, EK_LINK_DDD - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] - 1, coord[2] + 1, EK_LINK_DDU - 13)];
    flux_local_cartesian[2] +=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] - 1, EK_LINK_DUD - 13)];
    flux_local_cartesian[2] -=
        0.5 * fluxes[jindex_cartesian2linear_host(
                  coord[0] - 1, coord[1] + 1, coord[2] + 1, EK_LINK_DUU - 13)];

    for (int j = 0; j < 13; j++) {
      flux_local_linksum += fluxes[jindex_getByRhoLinear_host(i, j)];
    }

    fprintf(
        fp, "%e %e %e %e\n",
        flux_local_cartesian[0] / (ek_parameters.agrid * ek_parameters.agrid),
        flux_local_cartesian[1] / (ek_parameters.agrid * ek_parameters.agrid),
        flux_local_cartesian[2] / (ek_parameters.agrid * ek_parameters.agrid),
        flux_local_linksum / (ek_parameters.agrid * ek_parameters.agrid));
  }

  free(fluxes);
  fclose(fp);

  return 0;
#endif // EK_DEBUG
}

int ek_print_vtk_flux_link(int species, char *filename) {

  if (ek_parameters.species_index[species] == -1) {
    return 1;
  }

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  unsigned int coord[3];

  auto *fluxes = (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * 13 *
                                          sizeof(ekfloat));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_clear_fluxes, dim_grid, threads_per_block);
  KERNELCALL(ek_calculate_quantities, dim_grid, threads_per_block,
             ek_parameters.species_index[species], *current_nodes, node_f,
             ek_lbparameters_gpu, ek_lb_device_values, philox_counter.value());
  reset_LB_force_densities_GPU(false);

#ifdef EK_BOUNDARIES
  KERNELCALL(ek_apply_boundaries, dim_grid, threads_per_block,
             ek_parameters.species_index[species], *current_nodes, node_f);
#endif

  cuda_safe_mem(cudaMemcpy(fluxes, ek_parameters.j,
                           ek_parameters.number_of_nodes * 13 * sizeof(ekfloat),
                           cudaMemcpyDeviceToHost));

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
flux_%d\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS flux_%d float 3\n\
LOOKUP_TABLE default\n",
          species, ek_parameters.dim_x, ek_parameters.dim_y,
          ek_parameters.dim_z, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.number_of_nodes, species);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {
    rhoindex_linear2cartesian_host(i, coord);

    fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e \n",
            fluxes[jindex_getByRhoLinear_host(i, 0)],
            fluxes[jindex_getByRhoLinear_host(i, 1)],
            fluxes[jindex_getByRhoLinear_host(i, 2)],
            fluxes[jindex_getByRhoLinear_host(i, 3)],
            fluxes[jindex_getByRhoLinear_host(i, 4)],
            fluxes[jindex_getByRhoLinear_host(i, 5)],
            fluxes[jindex_getByRhoLinear_host(i, 6)],
            fluxes[jindex_getByRhoLinear_host(i, 7)],
            fluxes[jindex_getByRhoLinear_host(i, 8)],
            fluxes[jindex_getByRhoLinear_host(i, 9)],
            fluxes[jindex_getByRhoLinear_host(i, 10)],
            fluxes[jindex_getByRhoLinear_host(i, 11)],
            fluxes[jindex_getByRhoLinear_host(i, 12)]);
  }

  free(fluxes);
  fclose(fp);

  return 0;
}

int ek_node_print_potential(int x, int y, int z, double *potential) {
  auto i =
      static_cast<int>(z * ek_parameters.dim_y * ek_parameters.dim_x_padded +
                       y * ek_parameters.dim_x_padded + x);
  float pot;

  cuda_safe_mem(cudaMemcpy(&pot, &ek_parameters.charge_potential[i],
                           1 * sizeof(cufftReal), cudaMemcpyDeviceToHost));

  *potential = pot;
  return 0;
}

int ek_print_vtk_potential(char *filename) {

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  auto *potential =
      (float *)Utils::malloc(ek_parameters.number_of_nodes * sizeof(cufftReal));

  cuda_safe_mem(cudaMemcpy2D(potential, ek_parameters.dim_x * sizeof(cufftReal),
                             ek_parameters.charge_potential,
                             ek_parameters.dim_x_padded * sizeof(cufftReal),
                             ek_parameters.dim_x * sizeof(cufftReal),
                             ek_parameters.dim_z * ek_parameters.dim_y,
                             cudaMemcpyDeviceToHost));

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
potential\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS potential float 1\n\
LOOKUP_TABLE default\n",
          ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.agrid, ek_parameters.number_of_nodes);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {
    fprintf(fp, "%e\n", potential[i]);
  }

  free(potential);
  fclose(fp);

  return 0;
}

int ek_print_vtk_particle_potential(char *filename) {

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  auto *potential =
      (float *)Utils::malloc(ek_parameters.number_of_nodes * sizeof(cufftReal));

  cuda_safe_mem(cudaMemcpy2D(potential, ek_parameters.dim_x * sizeof(cufftReal),
                             ek_parameters.charge_potential_buffer,
                             ek_parameters.dim_x_padded * sizeof(cufftReal),
                             ek_parameters.dim_x * sizeof(cufftReal),
                             ek_parameters.dim_z * ek_parameters.dim_y,
                             cudaMemcpyDeviceToHost));

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
potential\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS potential float 1\n\
LOOKUP_TABLE default\n",
          ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.agrid, ek_parameters.number_of_nodes);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {
    fprintf(fp, "%e\n", potential[i]);
  }

  free(potential);
  fclose(fp);

  return 0;
}

int ek_print_vtk_lbforce_density(char *filename) {
#if !defined(VIRTUAL_SITES_INERTIALESS_TRACERS) && !defined(EK_DEBUG)
  throw std::runtime_error("Please rebuild ESPResSo with EK_DEBUG");
#else

  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  auto *lbforce_density = (lbForceFloat *)Utils::malloc(
      ek_parameters.number_of_nodes * 3 * sizeof(lbForceFloat));

  cuda_safe_mem(
      cudaMemcpy(lbforce_density, node_f.force_density_buf,
                 ek_parameters.number_of_nodes * 3 * sizeof(lbForceFloat),
                 cudaMemcpyDeviceToHost));

  fprintf(fp, "\
# vtk DataFile Version 2.0\n\
lbforce\n\
ASCII\n\
\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %u %u %u\n\
ORIGIN %f %f %f\n\
SPACING %f %f %f\n\
\n\
POINT_DATA %u\n\
SCALARS lbforce float 3\n\
LOOKUP_TABLE default\n",
          ek_parameters.dim_x, ek_parameters.dim_y, ek_parameters.dim_z,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid * 0.5f,
          ek_parameters.agrid * 0.5f, ek_parameters.agrid, ek_parameters.agrid,
          ek_parameters.agrid, ek_parameters.number_of_nodes);

  for (int i = 0; i < ek_parameters.number_of_nodes; i++) {
    fprintf(fp, "%e %e %e\n",
            lbforce_density[i] / (powf(ek_parameters.time_step, 2.0) *
                                  powf(ek_parameters.agrid, 4.0)),
            lbforce_density[i + ek_parameters.number_of_nodes] /
                (powf(ek_parameters.time_step, 2.0) *
                 powf(ek_parameters.agrid, 4.0)),
            lbforce_density[i + 2 * ek_parameters.number_of_nodes] /
                (powf(ek_parameters.time_step, 2.0) *
                 powf(ek_parameters.agrid, 4.0)));
  }

  free(lbforce_density);
  fclose(fp);

  return 0;
#endif
}

void ek_print_parameters() {

  printf("ek_parameters {\n");

  printf("  float agrid = %f;\n", ek_parameters.agrid);
  printf("  float time_step = %f;\n", ek_parameters.time_step);
  printf("  float lb_density = %f;\n", ek_parameters.lb_density);
  printf("  unsigned int dim_x = %d;\n", ek_parameters.dim_x);
  printf("  unsigned int dim_y = %d;\n", ek_parameters.dim_y);
  printf("  unsigned int dim_z = %d;\n", ek_parameters.dim_z);
  printf("  unsigned int number_of_nodes = %d;\n",
         ek_parameters.number_of_nodes);
  printf("  float viscosity = %f;\n", ek_parameters.viscosity);
  printf("  float bulk_viscosity = %f;\n", ek_parameters.bulk_viscosity);
  printf("  float gamma_odd = %f;\n", ek_parameters.gamma_odd);
  printf("  float gamma_even = %f;\n", ek_parameters.gamma_even);
  printf("  float friction = %f;\n", ek_parameters.friction);
  printf("  float T = %f;\n", ek_parameters.T);
  printf("  float prefactor = %f;\n", ek_parameters.prefactor);
  printf("  float lb_force_density[] = {%f, %f, %f};\n",
         ek_parameters.lb_force_density[0], ek_parameters.lb_force_density[1],
         ek_parameters.lb_force_density[2]);
  printf("  unsigned int number_of_species = %d;\n",
         ek_parameters.number_of_species);
  printf("  int reaction_species[] = {%d, %d, %d};\n",
         ek_parameters.reaction_species[0], ek_parameters.reaction_species[1],
         ek_parameters.reaction_species[2]);
  printf("  float rho_reactant_reservoir = %f;\n",
         ek_parameters.rho_reactant_reservoir);
  printf("  float rho_product0_reservoir = %f;\n",
         ek_parameters.rho_product0_reservoir);
  printf("  float rho_product1_reservoir = %f;\n",
         ek_parameters.rho_product1_reservoir);
  printf("  float reaction_ct_rate = %f;\n", ek_parameters.reaction_ct_rate);
  printf("  float reaction_fraction_0 = %f;\n",
         ek_parameters.reaction_fraction_0);
  printf("  float reaction_fraction_1 = %f;\n",
         ek_parameters.reaction_fraction_0);
  printf("  float* j = %p;\n", (void *)ek_parameters.j);

  printf("  float* rho[] = {%p, %p, %p, %p, %p, %p, %p, %p, %p, %p};\n",
         (void *)ek_parameters.rho[0], (void *)ek_parameters.rho[1],
         (void *)ek_parameters.rho[2], (void *)ek_parameters.rho[3],
         (void *)ek_parameters.rho[4], (void *)ek_parameters.rho[5],
         (void *)ek_parameters.rho[6], (void *)ek_parameters.rho[7],
         (void *)ek_parameters.rho[8], (void *)ek_parameters.rho[9]);

  printf("  int species_index[] = {%d, %d, %d, %d, %d, %d, %d, %d, %d, %d};\n",
         ek_parameters.species_index[0], ek_parameters.species_index[1],
         ek_parameters.species_index[2], ek_parameters.species_index[3],
         ek_parameters.species_index[4], ek_parameters.species_index[5],
         ek_parameters.species_index[6], ek_parameters.species_index[7],
         ek_parameters.species_index[8], ek_parameters.species_index[9]);

  printf("  float density = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
         ek_parameters.density[0], ek_parameters.density[1],
         ek_parameters.density[2], ek_parameters.density[3],
         ek_parameters.density[4], ek_parameters.density[5],
         ek_parameters.density[6], ek_parameters.density[7],
         ek_parameters.density[8], ek_parameters.density[9]);

  printf("  float D[] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
         ek_parameters.D[0], ek_parameters.D[1], ek_parameters.D[2],
         ek_parameters.D[3], ek_parameters.D[4], ek_parameters.D[5],
         ek_parameters.D[6], ek_parameters.D[7], ek_parameters.D[8],
         ek_parameters.D[9]);

  printf("  float d[] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
         ek_parameters.d[0], ek_parameters.d[1], ek_parameters.d[2],
         ek_parameters.d[3], ek_parameters.d[4], ek_parameters.d[5],
         ek_parameters.d[6], ek_parameters.d[7], ek_parameters.d[8],
         ek_parameters.d[9]);

  printf("  float valency[] = {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f};\n",
         ek_parameters.valency[0], ek_parameters.valency[1],
         ek_parameters.valency[2], ek_parameters.valency[3],
         ek_parameters.valency[4], ek_parameters.valency[5],
         ek_parameters.valency[6], ek_parameters.valency[7],
         ek_parameters.valency[8], ek_parameters.valency[9]);

  printf("  float ext_force_density[0][] = {%f, %f, %f, %f, %f, %f, %f, %f, "
         "%f, %f};\n",
         ek_parameters.ext_force_density[0][0],
         ek_parameters.ext_force_density[0][1],
         ek_parameters.ext_force_density[0][2],
         ek_parameters.ext_force_density[0][3],
         ek_parameters.ext_force_density[0][4],
         ek_parameters.ext_force_density[0][5],
         ek_parameters.ext_force_density[0][6],
         ek_parameters.ext_force_density[0][7],
         ek_parameters.ext_force_density[0][8],
         ek_parameters.ext_force_density[0][9]);

  printf("  float ext_force_density[1][] = {%f, %f, %f, %f, %f, %f, %f, %f, "
         "%f, %f};\n",
         ek_parameters.ext_force_density[1][0],
         ek_parameters.ext_force_density[1][1],
         ek_parameters.ext_force_density[1][2],
         ek_parameters.ext_force_density[1][3],
         ek_parameters.ext_force_density[1][4],
         ek_parameters.ext_force_density[1][5],
         ek_parameters.ext_force_density[1][6],
         ek_parameters.ext_force_density[1][7],
         ek_parameters.ext_force_density[1][8],
         ek_parameters.ext_force_density[1][9]);

  printf("  float ext_force_density[2][] = {%f, %f, %f, %f, %f, %f, %f, %f, "
         "%f, %f};\n",
         ek_parameters.ext_force_density[2][0],
         ek_parameters.ext_force_density[2][1],
         ek_parameters.ext_force_density[2][2],
         ek_parameters.ext_force_density[2][3],
         ek_parameters.ext_force_density[2][4],
         ek_parameters.ext_force_density[2][5],
         ek_parameters.ext_force_density[2][6],
         ek_parameters.ext_force_density[2][7],
         ek_parameters.ext_force_density[2][8],
         ek_parameters.ext_force_density[2][9]);

  printf("}\n");
}

void ek_print_lbpar() {

  printf("lbpar_gpu {\n");

  printf("    float rho = %f;\n", lbpar_gpu.rho);
  printf("    float mu = %f;\n", lbpar_gpu.mu);
  printf("    float viscosity = %f;\n", lbpar_gpu.viscosity);
  printf("    float gamma_shear = %f;\n", lbpar_gpu.gamma_shear);
  printf("    float gamma_bulk = %f;\n", lbpar_gpu.gamma_bulk);
  printf("    float gamma_odd = %f;\n", lbpar_gpu.gamma_odd);
  printf("    float gamma_even = %f;\n", lbpar_gpu.gamma_even);
  printf("    float agrid = %f;\n", lbpar_gpu.agrid);
  printf("    float tau = %f;\n", lbpar_gpu.tau);
  printf("    float time_step = %f;\n", lbpar_gpu.time_step);
  printf("    float bulk_viscosity = %f;\n", lbpar_gpu.bulk_viscosity);
  printf("    unsigned int dim_x = %d;\n", lbpar_gpu.dim_x);
  printf("    unsigned int dim_y = %d;\n", lbpar_gpu.dim_y);
  printf("    unsigned int dim_z = %d;\n", lbpar_gpu.dim_z);
  printf("    unsigned int number_of_nodes = %d;\n", lbpar_gpu.number_of_nodes);
  printf("    int calc_val = %d;\n", lbpar_gpu.calc_val);
  printf("    int external_force_density = %d;\n",
         lbpar_gpu.external_force_density);
  printf("    float ext_force_density[3] = {%f, %f, %f};\n",
         lbpar_gpu.ext_force_density[0], lbpar_gpu.ext_force_density[1],
         lbpar_gpu.ext_force_density[2]);
  printf("    unsigned int reinit = %d;\n", lbpar_gpu.reinit);

  printf("}\n");
}

int ek_set_agrid(float agrid) {

  ek_parameters.agrid = agrid;
  return 0;
}

int ek_set_lb_density(float lb_density) {

  ek_parameters.lb_density = lb_density;
  return 0;
}

int ek_set_prefactor(float prefactor) {

  ek_parameters.prefactor = prefactor;
  return 0;
}

int ek_set_electrostatics_coupling(bool electrostatics_coupling) {
  ek_parameters.es_coupling = electrostatics_coupling;
  return 0;
}

int ek_set_viscosity(float viscosity) {

  ek_parameters.viscosity = viscosity;
  return 0;
}

int ek_set_friction(float friction) {

  ek_parameters.friction = friction;
  return 0;
}

int ek_set_bulk_viscosity(float bulk_viscosity) {

  ek_parameters.bulk_viscosity = bulk_viscosity;
  return 0;
}

int ek_set_gamma_odd(float gamma_odd) {

  ek_parameters.gamma_odd = gamma_odd;
  return 0;
}

int ek_set_gamma_even(float gamma_even) {

  ek_parameters.gamma_even = gamma_even;
  return 0;
}

int ek_set_stencil(int stencil) {
  if (!ek_parameters.fluidcoupling_ideal_contribution)
    return 1; // combination not implemented

  ek_parameters.stencil = stencil;
  return 0;
}

int ek_set_advection(bool advection) {
  ek_parameters.advection = advection;
  return 0;
}

int ek_set_fluctuations(bool fluctuations) {
  ek_parameters.fluctuations = fluctuations;
  return 0;
}

int ek_set_fluctuation_amplitude(float fluctuation_amplitude) {
  ek_parameters.fluctuation_amplitude = fluctuation_amplitude;
  return 0;
}

int ek_set_fluidcoupling(bool ideal_contribution) {
  if (ek_parameters.stencil != 0)
    return 1; // combination not implemented

  ek_parameters.fluidcoupling_ideal_contribution = ideal_contribution;
  return 0;
}

int ek_set_density(int species, float density) {

  ek_init_species(species);

  ek_parameters.density[ek_parameters.species_index[species]] = density;

  return 0;
}

int ek_set_D(int species, float D) {

  ek_init_species(species);

  ek_parameters.D[ek_parameters.species_index[species]] = D;
  ek_parameters.d[ek_parameters.species_index[species]] =
      D / (1.0f + 2.0f * sqrt(2.0f));

  return 0;
}

int ek_set_T(float T) {

  ek_parameters.T = T;

  return 0;
}

int ek_set_valency(int species, float valency) {

  ek_init_species(species);

  ek_parameters.valency[ek_parameters.species_index[species]] = valency;

  return 0;
}

int ek_set_ext_force_density(int species, float ext_force_density_x,
                             float ext_force_density_y,
                             float ext_force_density_z) {

  ek_init_species(species);

  ek_parameters.ext_force_density[0][ek_parameters.species_index[species]] =
      ext_force_density_x;
  ek_parameters.ext_force_density[1][ek_parameters.species_index[species]] =
      ext_force_density_y;
  ek_parameters.ext_force_density[2][ek_parameters.species_index[species]] =
      ext_force_density_z;

  return 0;
}

struct ek_charge_of_particle {
  __host__ __device__ ekfloat operator()(CUDA_particle_data particle) {
    return particle.q;
  };
};

ekfloat ek_get_particle_charge() {
  auto device_particles = gpu_get_particle_pointer();
  ekfloat particle_charge = thrust::transform_reduce(
      thrust::device_ptr<CUDA_particle_data>(device_particles.begin()),
      thrust::device_ptr<CUDA_particle_data>(device_particles.end()),
      ek_charge_of_particle(), 0.0f, thrust::plus<ekfloat>());
  return particle_charge;
}

ekfloat ek_calculate_net_charge() {
  cuda_safe_mem(cudaMemset(charge_gpu, 0, sizeof(ekfloat)));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  auto blocks_per_grid_x =
      static_cast<int>((ek_parameters.number_of_nodes +
                        threads_per_block * blocks_per_grid_y - 1) /
                       (threads_per_block * blocks_per_grid_y));
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(ek_calculate_system_charge, dim_grid, threads_per_block,
             charge_gpu);

  ekfloat charge;
  cuda_safe_mem(
      cudaMemcpy(&charge, charge_gpu, sizeof(ekfloat), cudaMemcpyDeviceToHost));

  if (ek_parameters.es_coupling)
    charge += ek_get_particle_charge();

  return charge;
}

int ek_neutralize_system(int species) {
  int species_index = ek_parameters.species_index[species];

  if (species_index == -1)
    return 1;

  if (ek_parameters.valency[species_index] == 0.0f)
    return 2;

  ekfloat compensating_species_density = 0.0f;

#ifndef EK_BOUNDARIES
  for (int i = 0; i < ek_parameters.number_of_species; i++)
    compensating_species_density +=
        ek_parameters.density[i] * ek_parameters.valency[i];

  compensating_species_density =
      ek_parameters.density[species_index] -
      compensating_species_density / ek_parameters.valency[species_index];

  if (ek_parameters.es_coupling) {
    ekfloat particle_charge = ek_get_particle_charge();
    compensating_species_density -=
        particle_charge / ek_parameters.valency[species_index] /
        (ek_parameters.agrid * ek_parameters.agrid * ek_parameters.agrid) /
        double(ek_parameters.number_of_nodes);
  }

#else
  ekfloat charge = ek_calculate_net_charge();

  compensating_species_density =
      ek_parameters.density[species_index] -
      (charge / ek_parameters.valency[species_index]) /
          (ek_parameters.agrid * ek_parameters.agrid * ek_parameters.agrid *
           double(ek_parameters.number_of_nodes -
                  ek_parameters.number_of_boundary_nodes));
#endif // EK_BOUNDARIES

  if (compensating_species_density < 0.0f)
    return 3;

  ek_parameters.density[species_index] = compensating_species_density;

  return 0;
}

int ek_save_checkpoint(char *filename, char *lb_filename) {
  std::ofstream fout(filename, std::ofstream::binary);
  auto *densities =
      (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * sizeof(ekfloat));

  for (int i = 0; i < ek_parameters.number_of_species; i++) {
    cuda_safe_mem(cudaMemcpy(densities, ek_parameters.rho[i],
                             ek_parameters.number_of_nodes * sizeof(ekfloat),
                             cudaMemcpyDeviceToHost));

    if (!fout.write((char *)densities,
                    sizeof(ekfloat) * ek_parameters.number_of_nodes)) {
      free(densities);
      fout.close();
      return 1;
    }
  }

  free(densities);
  fout.close();

  lb_lbfluid_save_checkpoint(lb_filename, true);
  return 0;
}

int ek_load_checkpoint(char *filename) {
  std::string fname(filename);
  std::ifstream fin((const char *)(fname + ".ek").c_str(),
                    std::ifstream::binary);
  auto *densities =
      (ekfloat *)Utils::malloc(ek_parameters.number_of_nodes * sizeof(ekfloat));

  for (int i = 0; i < ek_parameters.number_of_species; i++) {
    if (!fin.read((char *)densities,
                  sizeof(ekfloat) * ek_parameters.number_of_nodes)) {
      free(densities);
      fin.close();
      return 1;
    }

    cuda_safe_mem(cudaMemcpy(ek_parameters.rho[i], densities,
                             ek_parameters.number_of_nodes * sizeof(ekfloat),
                             cudaMemcpyHostToDevice));
  }

  free(densities);
  fin.close();

  lb_lbfluid_load_checkpoint((char *)(fname + ".lb").c_str(), true);

  ek_integrate_electrostatics();

  return 0;
}

void ek_set_rng_state(uint64_t counter) {
  if (ek_initialized)
    philox_counter = Utils::Counter<uint64_t>(counter);
}

#endif /* ELECTROKINETICS */

#endif /* CUDA */

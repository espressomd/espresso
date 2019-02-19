/*
   Copyright (C) 2010-2018 The ESPResSo project

   This file is part of ESPResSo.

   ESPResSo is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ESPResSo is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *  %Lattice Boltzmann on GPUs.
 *
 *  The corresponding header file is lbgpu.cuh.
 */

#include "cuda_wrapper.hpp"
#include "curand_wrapper.hpp"

#include "config.hpp"

#ifdef LB_GPU
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "debug.hpp"
#include "errorhandling.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/electrokinetics_pdb_parse.hpp"
#include "grid_based_algorithms/lbgpu.cuh"
#include "grid_based_algorithms/lbgpu.hpp"
#include "utils/Counter.hpp"
#include "utils/Vector.hpp"

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/transform_reduce.h>

#include <cassert>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

/** defining structures residing in global memory */

/** device_rho_v: struct for hydrodynamic fields: this is for internal use
 *  (i.e. stores values in LB units) and should not used for
 *  printing values
 */
static LB_rho_v_gpu *device_rho_v = nullptr;

/** print_rho_v_pi: struct for hydrodynamic fields: this is the interface
 *  and stores values in MD units. It should not used
 *  as an input for any LB calculations. TODO: in the future,
 *  one might want to have several structures for printing
 *  separately rho, v, pi without having to compute/store
 *  the complete set.
 */
static LB_rho_v_pi_gpu *print_rho_v_pi = nullptr;

/** @name structs for velocity densities */
/*@{*/
static LB_nodes_gpu nodes_a = {nullptr, nullptr};
static LB_nodes_gpu nodes_b = {nullptr, nullptr};
;
/** struct for node force density*/

/** struct for node force density */
LB_node_force_density_gpu node_f = {nullptr, nullptr};

static LB_extern_nodeforcedensity_gpu *extern_node_force_densities = nullptr;

#ifdef LB_BOUNDARIES_GPU
static float *lb_boundary_force = nullptr;

static float *lb_boundary_velocity = nullptr;

/** @name pointers for bound index array */
/*@{*/
static int *boundary_node_list;
static int *boundary_index_list;
static size_t size_of_boundindex;
/*@}*/
#endif

EK_parameters *lb_ek_parameters_gpu;

/** @name pointers for additional cuda check flag */
/*@{*/
static int *gpu_check = nullptr;
static int *h_gpu_check = nullptr;
/*@}*/

static unsigned int intflag = 1;
LB_nodes_gpu *current_nodes = nullptr;
/** @name defining size values for allocating global memory */
/*@{*/
static size_t size_of_rho_v;
static size_t size_of_rho_v_pi;
static size_t size_of_extern_node_force_densities;
/*@}*/

/** Parameters residing in constant memory */
__device__ __constant__ LB_parameters_gpu para[1];
static const float c_sound_sq = 1.0f / 3.0f;

/*-------------------------------------------------------*/
/*********************************************************/
/** \name device functions called by kernel functions */
/*********************************************************/
/*-------------------------------------------------------*/

/*-------------------------------------------------------*/

static constexpr float sqrt12 = 3.4641016151377544f;
static Utils::Counter<uint64_t> rng_counter_coupling_gpu;
Utils::Counter<uint64_t> rng_counter_fluid_gpu;
__device__ float4 random_wrapper_philox(unsigned int index, unsigned int mode,
                                        uint64_t philox_counter) {
  // Split the 64 bit counter into two 32 bit ints.
  uint32_t philox_counter_hi = static_cast<uint32_t>(philox_counter >> 32);
  uint32_t philox_counter_low = static_cast<uint32_t>(philox_counter);
  uint4 rnd_ints =
      curand_Philox4x32_10(make_uint4(index, philox_counter_hi, 0, mode),
                           make_uint2(philox_counter_low, 0));
  float4 rnd_floats;
  rnd_floats.w = rnd_ints.w * CURAND_2POW32_INV + (CURAND_2POW32_INV / 2.0f);
  rnd_floats.x = rnd_ints.x * CURAND_2POW32_INV + (CURAND_2POW32_INV / 2.0f);
  rnd_floats.y = rnd_ints.y * CURAND_2POW32_INV + (CURAND_2POW32_INV / 2.0f);
  rnd_floats.z = rnd_ints.z * CURAND_2POW32_INV + (CURAND_2POW32_INV / 2.0f);
  return rnd_floats;
}

/** Transformation from 1d array-index to xyz
 *  @param[in]  index   Node index / thread index
 *  @param[out] xyz     Calculated xyz array
 */
template <typename T> __device__ void index_to_xyz(T index, T *xyz) {
  xyz[0] = index % para->dim_x;
  index /= para->dim_x;
  xyz[1] = index % para->dim_y;
  index /= para->dim_y;
  xyz[2] = index;
}

/** Transformation from xyz to 1d array-index
 *  @param[in] xyz     The xyz array
 */
template <typename T> __device__ T xyz_to_index(T *xyz) {
  T x = (xyz[0] + para->dim_x) % para->dim_x;
  T y = (xyz[1] + para->dim_y) % para->dim_y;
  T z = (xyz[2] + para->dim_z) % para->dim_z;
  return x + para->dim_x * (y + para->dim_y * z);
}

/** Calculate modes from the velocity densities (space-transform)
 *  @param[in]  n_a     Local node residing in array a
 *  @param[in]  index   Node index / thread index
 *  @param[out] mode    Local register values mode
 */
__device__ void calc_m_from_n(LB_nodes_gpu n_a, unsigned int index,
                              float *mode) {
  // The following convention is used:
  // The $\hat{c}_i$ form B. Duenweg's paper are given by:

  /* c_0  = { 0, 0, 0}
     c_1  = { 1, 0, 0}
     c_2  = {-1, 0, 0}
     c_3  = { 0, 1, 0}
     c_4  = { 0,-1, 0}
     c_5  = { 0, 0, 1}
     c_6  = { 0, 0,-1}
     c_7  = { 1, 1, 0}
     c_8  = {-1,-1, 0}
     c_9  = { 1,-1, 0}
     c_10 = {-1, 1, 0}
     c_11 = { 1, 0, 1}
     c_12 = {-1, 0,-1}
     c_13 = { 1, 0,-1}
     c_14 = {-1, 0, 1}
     c_15 = { 0, 1, 1}
     c_16 = { 0,-1,-1}
     c_17 = { 0, 1,-1}
     c_18 = { 0,-1, 1} */

  // The basis vectors (modes) are constructed as follows
  // $m_k = \sum_{i} e_{ki} n_{i}$, where the $e_{ki}$ form a
  // linear transformation (matrix) that is given by

  /* $e{ 0,i} = 1$
     $e{ 1,i} = c_{i,x}$
     $e{ 2,i} = c_{i,y}$
     $e{ 3,i} = c_{i,z}$
     $e{ 4,i} = c_{i}^2 - 1$
     $e{ 5,i} = c_{i,x}^2 - c_{i,y}^2$
     $e{ 6,i} = c_{i}^2 - 3*c_{i,z}^2$
     $e{ 7,i} = c_{i,x}*c_{i,y}$
     $e{ 8,i} = c_{i,x}*c_{i,z}$
     $e{ 9,i} = c_{i,y}*c_{i,z}$
     $e{10,i} = (3*c_{i}^2 - 5)*c_{i,x}$
     $e{11,i} = (3*c_{i}^2 - 5)*c_{i,y}$
     $e{12,i} = (3*c_{i}^2 - 5)*c_{i,z}$
     $e{13,i} = (c_{i,y}^2 - c_{i,z}^2)*c_{i,x}$
     $e{14,i} = (c_{i,x}^2 - c_{i,z}^2)*c_{i,y}$
     $e{15,i} = (c_{i,x}^2 - c_{i,y}^2)*c_{i,z}$
     $e{16,i} = 3*c_{i}^2^2 - 6*c_{i}^2 + 1$
     $e{17,i} = (2*c_{i}^2 - 3)*(c_{i,x}^2 - c_{i,y}^2)$
     $e{18,i} = (2*c_{i}^2 - 3)*(c_{i}^2 - 3*c_{i,z}^2)$ */

  // Such that the transformation matrix is given by

  /* {{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0},
      { 0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1},
      { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1},
      {-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      { 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1},
      { 0, 1, 1, 1, 1,-2,-2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1},
      { 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1},
      { 0,-2, 2, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0},
      { 0, 0, 0,-2, 2, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1},
      { 0, 0, 0, 0, 0,-2, 2, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1},
      { 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,-1, 1,-1, 1},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1,-1, 1, 1,-1},
      { 1,-2,-2,-2,-2,-2,-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      { 0,-1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1},
      { 0,-1,-1,-1,-1, 2, 2, 2, 2, 2, 2,-1,-1,-1,-1,-1,-1,-1,-1}} */

  // With weights

  /* q^{c_{i}} = { 1/3, 1/18, 1/18, 1/18,
                  1/18, 1/18, 1/18, 1/36,
                  1/36, 1/36, 1/36, 1/36,
                  1/36, 1/36, 1/36, 1/36,
                  1/36, 1/36, 1/36 } */

  // Which makes the transformation satisfy the following
  // orthogonality condition:
  // \sum_{i} q^{c_{i}} e_{ki} e_{li} = w_{k} \delta_{kl},
  // where the weights are:

  /* w_{i} = {  1, 1/3, 1/3, 1/3,
              2/3, 4/9, 4/3, 1/9,
              1/9, 1/9, 2/3, 2/3,
              2/3, 2/9, 2/9, 2/9,
                2, 4/9, 4/3 } */

  // mass mode

  mode[0] = n_a.vd[(0) * para->number_of_nodes + index] +
            n_a.vd[(1) * para->number_of_nodes + index] +
            n_a.vd[(2) * para->number_of_nodes + index] +
            n_a.vd[(3) * para->number_of_nodes + index] +
            n_a.vd[(4) * para->number_of_nodes + index] +
            n_a.vd[(5) * para->number_of_nodes + index] +
            n_a.vd[(6) * para->number_of_nodes + index] +
            n_a.vd[(7) * para->number_of_nodes + index] +
            n_a.vd[(8) * para->number_of_nodes + index] +
            n_a.vd[(9) * para->number_of_nodes + index] +
            n_a.vd[(10) * para->number_of_nodes + index] +
            n_a.vd[(11) * para->number_of_nodes + index] +
            n_a.vd[(12) * para->number_of_nodes + index] +
            n_a.vd[(13) * para->number_of_nodes + index] +
            n_a.vd[(14) * para->number_of_nodes + index] +
            n_a.vd[(15) * para->number_of_nodes + index] +
            n_a.vd[(16) * para->number_of_nodes + index] +
            n_a.vd[(17) * para->number_of_nodes + index] +
            n_a.vd[(18) * para->number_of_nodes + index];

  // momentum modes

  mode[1] = (n_a.vd[(1) * para->number_of_nodes + index] -
             n_a.vd[(2) * para->number_of_nodes + index]) +
            (n_a.vd[(7) * para->number_of_nodes + index] -
             n_a.vd[(8) * para->number_of_nodes + index]) +
            (n_a.vd[(9) * para->number_of_nodes + index] -
             n_a.vd[(10) * para->number_of_nodes + index]) +
            (n_a.vd[(11) * para->number_of_nodes + index] -
             n_a.vd[(12) * para->number_of_nodes + index]) +
            (n_a.vd[(13) * para->number_of_nodes + index] -
             n_a.vd[(14) * para->number_of_nodes + index]);

  mode[2] = (n_a.vd[(3) * para->number_of_nodes + index] -
             n_a.vd[(4) * para->number_of_nodes + index]) +
            (n_a.vd[(7) * para->number_of_nodes + index] -
             n_a.vd[(8) * para->number_of_nodes + index]) -
            (n_a.vd[(9) * para->number_of_nodes + index] -
             n_a.vd[(10) * para->number_of_nodes + index]) +
            (n_a.vd[(15) * para->number_of_nodes + index] -
             n_a.vd[(16) * para->number_of_nodes + index]) +
            (n_a.vd[(17) * para->number_of_nodes + index] -
             n_a.vd[(18) * para->number_of_nodes + index]);

  mode[3] = (n_a.vd[(5) * para->number_of_nodes + index] -
             n_a.vd[(6) * para->number_of_nodes + index]) +
            (n_a.vd[(11) * para->number_of_nodes + index] -
             n_a.vd[(12) * para->number_of_nodes + index]) -
            (n_a.vd[(13) * para->number_of_nodes + index] -
             n_a.vd[(14) * para->number_of_nodes + index]) +
            (n_a.vd[(15) * para->number_of_nodes + index] -
             n_a.vd[(16) * para->number_of_nodes + index]) -
            (n_a.vd[(17) * para->number_of_nodes + index] -
             n_a.vd[(18) * para->number_of_nodes + index]);

  // stress modes
  mode[4] = -n_a.vd[(0) * para->number_of_nodes + index] +
            n_a.vd[(7) * para->number_of_nodes + index] +
            n_a.vd[(8) * para->number_of_nodes + index] +
            n_a.vd[(9) * para->number_of_nodes + index] +
            n_a.vd[(10) * para->number_of_nodes + index] +
            n_a.vd[(11) * para->number_of_nodes + index] +
            n_a.vd[(12) * para->number_of_nodes + index] +
            n_a.vd[(13) * para->number_of_nodes + index] +
            n_a.vd[(14) * para->number_of_nodes + index] +
            n_a.vd[(15) * para->number_of_nodes + index] +
            n_a.vd[(16) * para->number_of_nodes + index] +
            n_a.vd[(17) * para->number_of_nodes + index] +
            n_a.vd[(18) * para->number_of_nodes + index];

  mode[5] = (n_a.vd[(1) * para->number_of_nodes + index] +
             n_a.vd[(2) * para->number_of_nodes + index]) -
            (n_a.vd[(3) * para->number_of_nodes + index] +
             n_a.vd[(4) * para->number_of_nodes + index]) +
            (n_a.vd[(11) * para->number_of_nodes + index] +
             n_a.vd[(12) * para->number_of_nodes + index]) +
            (n_a.vd[(13) * para->number_of_nodes + index] +
             n_a.vd[(14) * para->number_of_nodes + index]) -
            (n_a.vd[(15) * para->number_of_nodes + index] +
             n_a.vd[(16) * para->number_of_nodes + index]) -
            (n_a.vd[(17) * para->number_of_nodes + index] +
             n_a.vd[(18) * para->number_of_nodes + index]);

  mode[6] = (n_a.vd[(1) * para->number_of_nodes + index] +
             n_a.vd[(2) * para->number_of_nodes + index]) +
            (n_a.vd[(3) * para->number_of_nodes + index] +
             n_a.vd[(4) * para->number_of_nodes + index]) -
            (n_a.vd[(11) * para->number_of_nodes + index] +
             n_a.vd[(12) * para->number_of_nodes + index]) -
            (n_a.vd[(13) * para->number_of_nodes + index] +
             n_a.vd[(14) * para->number_of_nodes + index]) -
            (n_a.vd[(15) * para->number_of_nodes + index] +
             n_a.vd[(16) * para->number_of_nodes + index]) -
            (n_a.vd[(17) * para->number_of_nodes + index] +
             n_a.vd[(18) * para->number_of_nodes + index]) -
            2.0f * ((n_a.vd[(5) * para->number_of_nodes + index] +
                     n_a.vd[(6) * para->number_of_nodes + index]) -
                    (n_a.vd[(7) * para->number_of_nodes + index] +
                     n_a.vd[(8) * para->number_of_nodes + index]) -
                    (n_a.vd[(9) * para->number_of_nodes + index] +
                     n_a.vd[(10) * para->number_of_nodes + index]));

  mode[7] = (n_a.vd[(7) * para->number_of_nodes + index] +
             n_a.vd[(8) * para->number_of_nodes + index]) -
            (n_a.vd[(9) * para->number_of_nodes + index] +
             n_a.vd[(10) * para->number_of_nodes + index]);

  mode[8] = (n_a.vd[(11) * para->number_of_nodes + index] +
             n_a.vd[(12) * para->number_of_nodes + index]) -
            (n_a.vd[(13) * para->number_of_nodes + index] +
             n_a.vd[(14) * para->number_of_nodes + index]);

  mode[9] = (n_a.vd[(15) * para->number_of_nodes + index] +
             n_a.vd[(16) * para->number_of_nodes + index]) -
            (n_a.vd[(17) * para->number_of_nodes + index] +
             n_a.vd[(18) * para->number_of_nodes + index]);

  // kinetic modes

  mode[10] = -2.0f * (n_a.vd[(1) * para->number_of_nodes + index] -
                      n_a.vd[(2) * para->number_of_nodes + index]) +
             (n_a.vd[(7) * para->number_of_nodes + index] -
              n_a.vd[(8) * para->number_of_nodes + index]) +
             (n_a.vd[(9) * para->number_of_nodes + index] -
              n_a.vd[(10) * para->number_of_nodes + index]) +
             (n_a.vd[(11) * para->number_of_nodes + index] -
              n_a.vd[(12) * para->number_of_nodes + index]) +
             (n_a.vd[(13) * para->number_of_nodes + index] -
              n_a.vd[(14) * para->number_of_nodes + index]);

  mode[11] = -2.0f * (n_a.vd[(3) * para->number_of_nodes + index] -
                      n_a.vd[(4) * para->number_of_nodes + index]) +
             (n_a.vd[(7) * para->number_of_nodes + index] -
              n_a.vd[(8) * para->number_of_nodes + index]) -
             (n_a.vd[(9) * para->number_of_nodes + index] -
              n_a.vd[(10) * para->number_of_nodes + index]) +
             (n_a.vd[(15) * para->number_of_nodes + index] -
              n_a.vd[(16) * para->number_of_nodes + index]) +
             (n_a.vd[(17) * para->number_of_nodes + index] -
              n_a.vd[(18) * para->number_of_nodes + index]);

  mode[12] = -2.0f * (n_a.vd[(5) * para->number_of_nodes + index] -
                      n_a.vd[(6) * para->number_of_nodes + index]) +
             (n_a.vd[(11) * para->number_of_nodes + index] -
              n_a.vd[(12) * para->number_of_nodes + index]) -
             (n_a.vd[(13) * para->number_of_nodes + index] -
              n_a.vd[(14) * para->number_of_nodes + index]) +
             (n_a.vd[(15) * para->number_of_nodes + index] -
              n_a.vd[(16) * para->number_of_nodes + index]) -
             (n_a.vd[(17) * para->number_of_nodes + index] -
              n_a.vd[(18) * para->number_of_nodes + index]);

  mode[13] = (n_a.vd[(7) * para->number_of_nodes + index] -
              n_a.vd[(8) * para->number_of_nodes + index]) +
             (n_a.vd[(9) * para->number_of_nodes + index] -
              n_a.vd[(10) * para->number_of_nodes + index]) -
             (n_a.vd[(11) * para->number_of_nodes + index] -
              n_a.vd[(12) * para->number_of_nodes + index]) -
             (n_a.vd[(13) * para->number_of_nodes + index] -
              n_a.vd[(14) * para->number_of_nodes + index]);

  mode[14] = (n_a.vd[(7) * para->number_of_nodes + index] -
              n_a.vd[(8) * para->number_of_nodes + index]) -
             (n_a.vd[(9) * para->number_of_nodes + index] -
              n_a.vd[(10) * para->number_of_nodes + index]) -
             (n_a.vd[(15) * para->number_of_nodes + index] -
              n_a.vd[(16) * para->number_of_nodes + index]) -
             (n_a.vd[(17) * para->number_of_nodes + index] -
              n_a.vd[(18) * para->number_of_nodes + index]);

  mode[15] = (n_a.vd[(11) * para->number_of_nodes + index] -
              n_a.vd[(12) * para->number_of_nodes + index]) -
             (n_a.vd[(13) * para->number_of_nodes + index] -
              n_a.vd[(14) * para->number_of_nodes + index]) -
             (n_a.vd[(15) * para->number_of_nodes + index] -
              n_a.vd[(16) * para->number_of_nodes + index]) +
             (n_a.vd[(17) * para->number_of_nodes + index] -
              n_a.vd[(18) * para->number_of_nodes + index]);

  mode[16] = n_a.vd[(0) * para->number_of_nodes + index] +
             n_a.vd[(7) * para->number_of_nodes + index] +
             n_a.vd[(8) * para->number_of_nodes + index] +
             n_a.vd[(9) * para->number_of_nodes + index] +
             n_a.vd[(10) * para->number_of_nodes + index] +
             n_a.vd[(11) * para->number_of_nodes + index] +
             n_a.vd[(12) * para->number_of_nodes + index] +
             n_a.vd[(13) * para->number_of_nodes + index] +
             n_a.vd[(14) * para->number_of_nodes + index] +
             n_a.vd[(15) * para->number_of_nodes + index] +
             n_a.vd[(16) * para->number_of_nodes + index] +
             n_a.vd[(17) * para->number_of_nodes + index] +
             n_a.vd[(18) * para->number_of_nodes + index] -
             2.0f * ((n_a.vd[(1) * para->number_of_nodes + index] +
                      n_a.vd[(2) * para->number_of_nodes + index]) +
                     (n_a.vd[(3) * para->number_of_nodes + index] +
                      n_a.vd[(4) * para->number_of_nodes + index]) +
                     (n_a.vd[(5) * para->number_of_nodes + index] +
                      n_a.vd[(6) * para->number_of_nodes + index]));

  mode[17] = -(n_a.vd[(1) * para->number_of_nodes + index] +
               n_a.vd[(2) * para->number_of_nodes + index]) +
             (n_a.vd[(3) * para->number_of_nodes + index] +
              n_a.vd[(4) * para->number_of_nodes + index]) +
             (n_a.vd[(11) * para->number_of_nodes + index] +
              n_a.vd[(12) * para->number_of_nodes + index]) +
             (n_a.vd[(13) * para->number_of_nodes + index] +
              n_a.vd[(14) * para->number_of_nodes + index]) -
             (n_a.vd[(15) * para->number_of_nodes + index] +
              n_a.vd[(16) * para->number_of_nodes + index]) -
             (n_a.vd[(17) * para->number_of_nodes + index] +
              n_a.vd[(18) * para->number_of_nodes + index]);

  mode[18] = -(n_a.vd[(1) * para->number_of_nodes + index] +
               n_a.vd[(2) * para->number_of_nodes + index]) -
             (n_a.vd[(3) * para->number_of_nodes + index] +
              n_a.vd[(4) * para->number_of_nodes + index]) -
             (n_a.vd[(11) * para->number_of_nodes + index] +
              n_a.vd[(12) * para->number_of_nodes + index]) -
             (n_a.vd[(13) * para->number_of_nodes + index] +
              n_a.vd[(14) * para->number_of_nodes + index]) -
             (n_a.vd[(15) * para->number_of_nodes + index] +
              n_a.vd[(16) * para->number_of_nodes + index]) -
             (n_a.vd[(17) * para->number_of_nodes + index] +
              n_a.vd[(18) * para->number_of_nodes + index]) +
             2.0f * ((n_a.vd[(5) * para->number_of_nodes + index] +
                      n_a.vd[(6) * para->number_of_nodes + index]) +
                     (n_a.vd[(7) * para->number_of_nodes + index] +
                      n_a.vd[(8) * para->number_of_nodes + index]) +
                     (n_a.vd[(9) * para->number_of_nodes + index] +
                      n_a.vd[(10) * para->number_of_nodes + index]));
}

__device__ void reset_LB_force_densities(unsigned int index,
                                         LB_node_force_density_gpu node_f,
                                         bool buffer = true) {
#if defined(VIRTUAL_SITES_INERTIALESS_TRACERS) || defined(EK_DEBUG)
  // Store backup of the node forces
  if (buffer) {
    node_f.force_density_buf[0 * para->number_of_nodes + index] =
        node_f.force_density[0 * para->number_of_nodes + index];
    node_f.force_density_buf[1 * para->number_of_nodes + index] =
        node_f.force_density[1 * para->number_of_nodes + index];
    node_f.force_density_buf[2 * para->number_of_nodes + index] =
        node_f.force_density[2 * para->number_of_nodes + index];
  }
#endif

#ifdef EXTERNAL_FORCES
  if (para->external_force_density) {
    node_f.force_density[0 * para->number_of_nodes + index] =
        para->ext_force_density[0];
    node_f.force_density[1 * para->number_of_nodes + index] =
        para->ext_force_density[1];
    node_f.force_density[2 * para->number_of_nodes + index] =
        para->ext_force_density[2];
  } else {
    node_f.force_density[0 * para->number_of_nodes + index] = 0.0f;
    node_f.force_density[1 * para->number_of_nodes + index] = 0.0f;
    node_f.force_density[2 * para->number_of_nodes + index] = 0.0f;
  }
#else
  /* reset force */
  node_f.force_density[0 * para->number_of_nodes + index] = 0.0f;
  node_f.force_density[1 * para->number_of_nodes + index] = 0.0f;
  node_f.force_density[2 * para->number_of_nodes + index] = 0.0f;
#endif
}

__global__ void
reset_LB_force_densities_kernel(LB_node_force_density_gpu node_f,
                                bool buffer = true) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes)
    reset_LB_force_densities(index, node_f, buffer);
}

void reset_LB_force_densities_GPU(bool buffer) {
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_LB_force_densities_kernel, dim_grid, threads_per_block,
             node_f, buffer);
}

/**
 *  @param[in]  mode    Local register values mode
 *  @param[in]  index   Node index / thread index
 *  @param[in]  node_f  Local node force
 *  @param[out] d_v     Local device values
 */
__device__ void update_rho_v(float *mode, unsigned int index,
                             LB_node_force_density_gpu node_f,
                             LB_rho_v_gpu *d_v) {
  float Rho_tot = 0.0f;
  float u_tot[3] = {0.0f, 0.0f, 0.0f};

  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */

  d_v[index].rho = mode[0] + para->rho;
  Rho_tot += mode[0] + para->rho;
  u_tot[0] += mode[1];
  u_tot[1] += mode[2];
  u_tot[2] += mode[3];

  /** if forces are present, the momentum density is redefined to
   * include one half-step of the force action.  See the
   * Chapman-Enskog expansion in [Ladd & Verberg]. */

  u_tot[0] += 0.5f * node_f.force_density[0 * para->number_of_nodes + index];
  u_tot[1] += 0.5f * node_f.force_density[1 * para->number_of_nodes + index];
  u_tot[2] += 0.5f * node_f.force_density[2 * para->number_of_nodes + index];

  u_tot[0] /= Rho_tot;
  u_tot[1] /= Rho_tot;
  u_tot[2] /= Rho_tot;

  d_v[index].v[0] = u_tot[0];
  d_v[index].v[1] = u_tot[1];
  d_v[index].v[2] = u_tot[2];
}

/** lb_relax_modes, means collision update of the modes
 *  @param[in] index     Node index / thread index
 *  @param[in,out] mode  Local register values mode
 *  @param[in] node_f    Local node force
 *  @param[in,out] d_v   Local device values
 */
__device__ void relax_modes(float *mode, unsigned int index,
                            LB_node_force_density_gpu node_f,
                            LB_rho_v_gpu *d_v) {
  float u_tot[3] = {0.0f, 0.0f, 0.0f};

  update_rho_v(mode, index, node_f, d_v);

  u_tot[0] = d_v[index].v[0];
  u_tot[1] = d_v[index].v[1];
  u_tot[2] = d_v[index].v[2];

  float Rho;
  float j[3];
  float modes_from_pi_eq[6];

  Rho = mode[0] + para->rho;
  j[0] = Rho * u_tot[0];
  j[1] = Rho * u_tot[1];
  j[2] = Rho * u_tot[2];

  /** equilibrium part of the stress modes (eq13 schiller) */

  modes_from_pi_eq[0] = ((j[0] * j[0]) + (j[1] * j[1]) + (j[2] * j[2])) / Rho;
  modes_from_pi_eq[1] = ((j[0] * j[0]) - (j[1] * j[1])) / Rho;
  modes_from_pi_eq[2] =
      (((j[0] * j[0]) + (j[1] * j[1]) + (j[2] * j[2])) - 3.0f * (j[2] * j[2])) /
      Rho;
  modes_from_pi_eq[3] = j[0] * j[1] / Rho;
  modes_from_pi_eq[4] = j[0] * j[2] / Rho;
  modes_from_pi_eq[5] = j[1] * j[2] / Rho;

  /** relax the stress modes (eq14 schiller) */

  mode[4] =
      modes_from_pi_eq[0] + para->gamma_bulk * (mode[4] - modes_from_pi_eq[0]);
  mode[5] =
      modes_from_pi_eq[1] + para->gamma_shear * (mode[5] - modes_from_pi_eq[1]);
  mode[6] =
      modes_from_pi_eq[2] + para->gamma_shear * (mode[6] - modes_from_pi_eq[2]);
  mode[7] =
      modes_from_pi_eq[3] + para->gamma_shear * (mode[7] - modes_from_pi_eq[3]);
  mode[8] =
      modes_from_pi_eq[4] + para->gamma_shear * (mode[8] - modes_from_pi_eq[4]);
  mode[9] =
      modes_from_pi_eq[5] + para->gamma_shear * (mode[9] - modes_from_pi_eq[5]);

  /** relax the ghost modes (project them out) */
  /** ghost modes have no equilibrium part due to orthogonality */

  mode[10] = para->gamma_odd * mode[10];
  mode[11] = para->gamma_odd * mode[11];
  mode[12] = para->gamma_odd * mode[12];
  mode[13] = para->gamma_odd * mode[13];
  mode[14] = para->gamma_odd * mode[14];
  mode[15] = para->gamma_odd * mode[15];
  mode[16] = para->gamma_even * mode[16];
  mode[17] = para->gamma_even * mode[17];
  mode[18] = para->gamma_even * mode[18];
}

/** Thermalization of the modes with Gaussian random numbers
 *  @param[in] index     Node index / thread index
 *  @param[in,out] mode  Local register values mode
 *  @param[in] philox_counter
 */
__device__ void thermalize_modes(float *mode, unsigned int index,
                                 uint64_t philox_counter) {
  float Rho;
  float4 random_floats;
  /** mass mode */
  Rho = mode[0] + para->rho;

  /* momentum modes */

  /* stress modes */
  random_floats = random_wrapper_philox(index, 4, philox_counter);
  mode[4] += sqrtf(Rho * (para->mu * (2.0f / 3.0f) *
                          (1.0f - (para->gamma_bulk * para->gamma_bulk)))) *
             (random_floats.w - 0.5f) * sqrt12;
  mode[5] += sqrtf(Rho * (para->mu * (4.0f / 9.0f) *
                          (1.0f - (para->gamma_shear * para->gamma_shear)))) *
             (random_floats.x - 0.5f) * sqrt12;

  mode[6] += sqrtf(Rho * (para->mu * (4.0f / 3.0f) *
                          (1.0f - (para->gamma_shear * para->gamma_shear)))) *
             (random_floats.y - 0.5f) * sqrt12;
  mode[7] += sqrtf(Rho * (para->mu * (1.0f / 9.0f) *
                          (1.0f - (para->gamma_shear * para->gamma_shear)))) *
             (random_floats.z - 0.5f) * sqrt12;

  random_floats = random_wrapper_philox(index, 8, philox_counter);
  mode[8] += sqrtf(Rho * (para->mu * (1.0f / 9.0f) *
                          (1.0f - (para->gamma_shear * para->gamma_shear)))) *
             (random_floats.w - 0.5f) * sqrt12;
  mode[9] += sqrtf(Rho * (para->mu * (1.0f / 9.0f) *
                          (1.0f - (para->gamma_shear * para->gamma_shear)))) *
             (random_floats.x - 0.5f) * sqrt12;

  /* ghost modes */
  mode[10] += sqrtf(Rho * (para->mu * (2.0f / 3.0f) *
                           (1.0f - (para->gamma_odd * para->gamma_odd)))) *
              (random_floats.y - 0.5f) * sqrt12;
  mode[11] += sqrtf(Rho * (para->mu * (2.0f / 3.0f) *
                           (1.0f - (para->gamma_odd * para->gamma_odd)))) *
              (random_floats.z - 0.5f) * sqrt12;

  random_floats = random_wrapper_philox(index, 12, philox_counter);
  mode[12] += sqrtf(Rho * (para->mu * (2.0f / 3.0f) *
                           (1.0f - (para->gamma_odd * para->gamma_odd)))) *
              (random_floats.w - 0.5f) * sqrt12;
  mode[13] += sqrtf(Rho * (para->mu * (2.0f / 9.0f) *
                           (1.0f - (para->gamma_odd * para->gamma_odd)))) *
              (random_floats.x - 0.5f) * sqrt12;

  mode[14] += sqrtf(Rho * (para->mu * (2.0f / 9.0f) *
                           (1.0f - (para->gamma_odd * para->gamma_odd)))) *
              (random_floats.y - 0.5f) * sqrt12;
  mode[15] += sqrtf(Rho * (para->mu * (2.0f / 9.0f) *
                           (1.0f - (para->gamma_odd * para->gamma_odd)))) *
              (random_floats.z - 0.5f) * sqrt12;

  random_floats = random_wrapper_philox(index, 16, philox_counter);
  mode[16] += sqrtf(Rho * (para->mu * (2.0f) *
                           (1.0f - (para->gamma_even * para->gamma_even)))) *
              (random_floats.w - 0.5f) * sqrt12;
  mode[17] += sqrtf(Rho * (para->mu * (4.0f / 9.0f) *
                           (1.0f - (para->gamma_even * para->gamma_even)))) *
              (random_floats.x - 0.5f) * sqrt12;

  mode[18] += sqrtf(Rho * (para->mu * (4.0f / 3.0f) *
                           (1.0f - (para->gamma_even * para->gamma_even)))) *
              (random_floats.y - 0.5f) * sqrt12;
}

/*-------------------------------------------------------*/
/** Normalization of the modes need before back-transformation into velocity
 *  space
 *  @param[in,out] mode  Local register values mode
 */
__device__ void normalize_modes(float *mode) {
  /* normalization factors enter in the back transformation */
  mode[0] *= 1.0f;
  mode[1] *= 3.0f;
  mode[2] *= 3.0f;
  mode[3] *= 3.0f;
  mode[4] *= 3.0f / 2.0f;
  mode[5] *= 9.0f / 4.0f;
  mode[6] *= 3.0f / 4.0f;
  mode[7] *= 9.0f;
  mode[8] *= 9.0f;
  mode[9] *= 9.0f;
  mode[10] *= 3.0f / 2.0f;
  mode[11] *= 3.0f / 2.0f;
  mode[12] *= 3.0f / 2.0f;
  mode[13] *= 9.0f / 2.0f;
  mode[14] *= 9.0f / 2.0f;
  mode[15] *= 9.0f / 2.0f;
  mode[16] *= 1.0f / 2.0f;
  mode[17] *= 9.0f / 4.0f;
  mode[18] *= 3.0f / 4.0f;
}

/*-------------------------------------------------------*/
/** Back-transformation from modespace to densityspace and streaming with
 *  the push method using pbc
 *  @param[in]  index  Node index / thread index
 *  @param[in]  mode   Local register values mode
 *  @param[out] n_b    Local node residing in array b
 */
__device__ void calc_n_from_modes_push(LB_nodes_gpu n_b, float *mode,
                                       unsigned int index) {
  unsigned int xyz[3];
  index_to_xyz(index, xyz);
  unsigned int x = xyz[0];
  unsigned int y = xyz[1];
  unsigned int z = xyz[2];

  n_b.vd[0 * para->number_of_nodes + x + para->dim_x * y +
         para->dim_x * para->dim_y * z] =
      1.0f / 3.0f * (mode[0] - mode[4] + mode[16]);

  n_b.vd[1 * para->number_of_nodes + (x + 1) % para->dim_x + para->dim_x * y +
         para->dim_x * para->dim_y * z] =
      1.0f / 18.0f *
      (mode[0] + mode[1] + mode[5] + mode[6] - mode[17] - mode[18] -
       2.0f * (mode[10] + mode[16]));

  n_b.vd[2 * para->number_of_nodes + (para->dim_x + x - 1) % para->dim_x +
         para->dim_x * y + para->dim_x * para->dim_y * z] =
      1.0f / 18.0f *
      (mode[0] - mode[1] + mode[5] + mode[6] - mode[17] - mode[18] +
       2.0f * (mode[10] - mode[16]));

  n_b.vd[3 * para->number_of_nodes + x + para->dim_x * ((y + 1) % para->dim_y) +
         para->dim_x * para->dim_y * z] =
      1.0f / 18.0f *
      (mode[0] + mode[2] - mode[5] + mode[6] + mode[17] - mode[18] -
       2.0f * (mode[11] + mode[16]));

  n_b.vd[4 * para->number_of_nodes + x +
         para->dim_x * ((para->dim_y + y - 1) % para->dim_y) +
         para->dim_x * para->dim_y * z] =
      1.0f / 18.0f *
      (mode[0] - mode[2] - mode[5] + mode[6] + mode[17] - mode[18] +
       2.0f * (mode[11] - mode[16]));

  n_b.vd[5 * para->number_of_nodes + x + para->dim_x * y +
         para->dim_x * para->dim_y * ((z + 1) % para->dim_z)] =
      1.0f / 18.0f *
      (mode[0] + mode[3] - 2.0f * (mode[6] + mode[12] + mode[16] - mode[18]));

  n_b.vd[6 * para->number_of_nodes + x + para->dim_x * y +
         para->dim_x * para->dim_y * ((para->dim_z + z - 1) % para->dim_z)] =
      1.0f / 18.0f *
      (mode[0] - mode[3] - 2.0f * (mode[6] - mode[12] + mode[16] - mode[18]));

  n_b.vd[7 * para->number_of_nodes + (x + 1) % para->dim_x +
         para->dim_x * ((y + 1) % para->dim_y) +
         para->dim_x * para->dim_y * z] =
      1.0f / 36.0f *
      (mode[0] + mode[1] + mode[2] + mode[4] + 2.0f * mode[6] + mode[7] +
       mode[10] + mode[11] + mode[13] + mode[14] + mode[16] + 2.0f * mode[18]);

  n_b.vd[8 * para->number_of_nodes + (para->dim_x + x - 1) % para->dim_x +
         para->dim_x * ((para->dim_y + y - 1) % para->dim_y) +
         para->dim_x * para->dim_y * z] =
      1.0f / 36.0f *
      (mode[0] - mode[1] - mode[2] + mode[4] + 2.0f * mode[6] + mode[7] -
       mode[10] - mode[11] - mode[13] - mode[14] + mode[16] + 2.0f * mode[18]);

  n_b.vd[9 * para->number_of_nodes + (x + 1) % para->dim_x +
         para->dim_x * ((para->dim_y + y - 1) % para->dim_y) +
         para->dim_x * para->dim_y * z] =
      1.0f / 36.0f *
      (mode[0] + mode[1] - mode[2] + mode[4] + 2.0f * mode[6] - mode[7] +
       mode[10] - mode[11] + mode[13] - mode[14] + mode[16] + 2.0f * mode[18]);

  n_b.vd[10 * para->number_of_nodes + (para->dim_x + x - 1) % para->dim_x +
         para->dim_x * ((y + 1) % para->dim_y) +
         para->dim_x * para->dim_y * z] =
      1.0f / 36.0f *
      (mode[0] - mode[1] + mode[2] + mode[4] + 2.0f * mode[6] - mode[7] -
       mode[10] + mode[11] - mode[13] + mode[14] + mode[16] + 2.0f * mode[18]);

  n_b.vd[11 * para->number_of_nodes + (x + 1) % para->dim_x + para->dim_x * y +
         para->dim_x * para->dim_y * ((z + 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] + mode[1] + mode[3] + mode[4] + mode[5] - mode[6] + mode[8] +
       mode[10] + mode[12] - mode[13] + mode[15] + mode[16] + mode[17] -
       mode[18]);

  n_b.vd[12 * para->number_of_nodes + (para->dim_x + x - 1) % para->dim_x +
         para->dim_x * y +
         para->dim_x * para->dim_y * ((para->dim_z + z - 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] - mode[1] - mode[3] + mode[4] + mode[5] - mode[6] + mode[8] -
       mode[10] - mode[12] + mode[13] - mode[15] + mode[16] + mode[17] -
       mode[18]);

  n_b.vd[13 * para->number_of_nodes + (x + 1) % para->dim_x + para->dim_x * y +
         para->dim_x * para->dim_y * ((para->dim_z + z - 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] + mode[1] - mode[3] + mode[4] + mode[5] - mode[6] - mode[8] +
       mode[10] - mode[12] - mode[13] - mode[15] + mode[16] + mode[17] -
       mode[18]);

  n_b.vd[14 * para->number_of_nodes + (para->dim_x + x - 1) % para->dim_x +
         para->dim_x * y +
         para->dim_x * para->dim_y * ((z + 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] - mode[1] + mode[3] + mode[4] + mode[5] - mode[6] - mode[8] -
       mode[10] + mode[12] + mode[13] + mode[15] + mode[16] + mode[17] -
       mode[18]);

  n_b.vd[15 * para->number_of_nodes + x +
         para->dim_x * ((y + 1) % para->dim_y) +
         para->dim_x * para->dim_y * ((z + 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] + mode[2] + mode[3] + mode[4] - mode[5] - mode[6] + mode[9] +
       mode[11] + mode[12] - mode[14] - mode[15] + mode[16] - mode[17] -
       mode[18]);

  n_b.vd[16 * para->number_of_nodes + x +
         para->dim_x * ((para->dim_y + y - 1) % para->dim_y) +
         para->dim_x * para->dim_y * ((para->dim_z + z - 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] - mode[2] - mode[3] + mode[4] - mode[5] - mode[6] + mode[9] -
       mode[11] - mode[12] + mode[14] + mode[15] + mode[16] - mode[17] -
       mode[18]);

  n_b.vd[17 * para->number_of_nodes + x +
         para->dim_x * ((y + 1) % para->dim_y) +
         para->dim_x * para->dim_y * ((para->dim_z + z - 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] + mode[2] - mode[3] + mode[4] - mode[5] - mode[6] - mode[9] +
       mode[11] - mode[12] - mode[14] + mode[15] + mode[16] - mode[17] -
       mode[18]);

  n_b.vd[18 * para->number_of_nodes + x +
         para->dim_x * ((para->dim_y + y - 1) % para->dim_y) +
         para->dim_x * para->dim_y * ((z + 1) % para->dim_z)] =
      1.0f / 36.0f *
      (mode[0] - mode[2] + mode[3] + mode[4] - mode[5] - mode[6] - mode[9] -
       mode[11] + mode[12] + mode[14] - mode[15] + mode[16] - mode[17] -
       mode[18]);
}

/** Bounce back boundary conditions.
 *
 *  The populations that have propagated into a boundary node
 *  are bounced back to the node they came from. This results
 *  in no slip boundary conditions.
 *
 *  [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 *  @param[in]  index   Node index / thread index
 *  @param[in]  n_curr  Local node receiving the current node field
 *  @param[in]  lb_boundary_velocity  Constant velocity at the boundary,
 *                                    set by the user
 *  @param[out] lb_boundary_force     Force on the boundary nodes
 */
__device__ void bounce_back_boundaries(LB_nodes_gpu n_curr, unsigned int index,
                                       float *lb_boundary_velocity,
                                       float *lb_boundary_force) {
  unsigned int xyz[3];
  int c[3];
  float v[3];
  float shift, weight, pop_to_bounce_back;
  float boundary_force[3] = {0.0f, 0.0f, 0.0f};
  size_t to_index, to_index_x, to_index_y, to_index_z;
  int population, inverse;
  int boundary_index;

  boundary_index = n_curr.boundary[index];
  if (boundary_index != 0) {
    v[0] = lb_boundary_velocity[3 * (boundary_index - 1) + 0];
    v[1] = lb_boundary_velocity[3 * (boundary_index - 1) + 1];
    v[2] = lb_boundary_velocity[3 * (boundary_index - 1) + 2];

    index_to_xyz(index, xyz);

    unsigned int x = xyz[0];
    unsigned int y = xyz[1];
    unsigned int z = xyz[2];

    /** store vd temporary in second lattice to avoid race conditions */

    // TODO : PUT IN EQUILIBRIUM CONTRIBUTION TO THE BOUNCE-BACK DENSITY FOR THE
    // BOUNDARY FORCE
    // TODO : INITIALIZE BOUNDARY FORCE PROPERLY, HAS NONZERO ELEMENTS IN FIRST
    // STEP
    // TODO : SET INTERNAL BOUNDARY NODE VALUES TO ZERO

#define BOUNCEBACK()                                                           \
  shift = 2.0f / para->agrid * para->rho * 3.0f * weight * para->tau *         \
          (v[0] * c[0] + v[1] * c[1] + v[2] * c[2]);                           \
  pop_to_bounce_back = n_curr.vd[population * para->number_of_nodes + index];  \
  to_index_x = (x + c[0] + para->dim_x) % para->dim_x;                         \
  to_index_y = (y + c[1] + para->dim_y) % para->dim_y;                         \
  to_index_z = (z + c[2] + para->dim_z) % para->dim_z;                         \
  to_index = to_index_x + para->dim_x * to_index_y +                           \
             para->dim_x * para->dim_y * to_index_z;                           \
  if (n_curr.boundary[to_index] == 0) {                                        \
    boundary_force[0] += (2.0f * pop_to_bounce_back + shift) * c[0] /          \
                         para->tau / para->tau / para->agrid;                  \
    boundary_force[1] += (2.0f * pop_to_bounce_back + shift) * c[1] /          \
                         para->tau / para->tau / para->agrid;                  \
    boundary_force[2] += (2.0f * pop_to_bounce_back + shift) * c[2] /          \
                         para->tau / para->tau / para->agrid;                  \
    n_curr.vd[inverse * para->number_of_nodes + to_index] =                    \
        pop_to_bounce_back + shift;                                            \
  }

    // the resting population does nothing, i.e., population 0.
    c[0] = 1;
    c[1] = 0;
    c[2] = 0;
    weight = 1. / 18.;
    population = 2;
    inverse = 1;
    BOUNCEBACK();

    c[0] = -1;
    c[1] = 0;
    c[2] = 0;
    weight = 1. / 18.;
    population = 1;
    inverse = 2;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = 1;
    c[2] = 0;
    weight = 1. / 18.;
    population = 4;
    inverse = 3;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = -1;
    c[2] = 0;
    weight = 1. / 18.;
    population = 3;
    inverse = 4;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = 0;
    c[2] = 1;
    weight = 1. / 18.;
    population = 6;
    inverse = 5;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = 0;
    c[2] = -1;
    weight = 1. / 18.;
    population = 5;
    inverse = 6;
    BOUNCEBACK();

    c[0] = 1;
    c[1] = 1;
    c[2] = 0;
    weight = 1. / 36.;
    population = 8;
    inverse = 7;
    BOUNCEBACK();

    c[0] = -1;
    c[1] = -1;
    c[2] = 0;
    weight = 1. / 36.;
    population = 7;
    inverse = 8;
    BOUNCEBACK();

    c[0] = 1;
    c[1] = -1;
    c[2] = 0;
    weight = 1. / 36.;
    population = 10;
    inverse = 9;
    BOUNCEBACK();

    c[0] = -1;
    c[1] = 1;
    c[2] = 0;
    weight = 1. / 36.;
    population = 9;
    inverse = 10;
    BOUNCEBACK();

    c[0] = 1;
    c[1] = 0;
    c[2] = 1;
    weight = 1. / 36.;
    population = 12;
    inverse = 11;
    BOUNCEBACK();

    c[0] = -1;
    c[1] = 0;
    c[2] = -1;
    weight = 1. / 36.;
    population = 11;
    inverse = 12;
    BOUNCEBACK();

    c[0] = 1;
    c[1] = 0;
    c[2] = -1;
    weight = 1. / 36.;
    population = 14;
    inverse = 13;
    BOUNCEBACK();

    c[0] = -1;
    c[1] = 0;
    c[2] = 1;
    weight = 1. / 36.;
    population = 13;
    inverse = 14;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = 1;
    c[2] = 1;
    weight = 1. / 36.;
    population = 16;
    inverse = 15;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = -1;
    c[2] = -1;
    weight = 1. / 36.;
    population = 15;
    inverse = 16;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = 1;
    c[2] = -1;
    weight = 1. / 36.;
    population = 18;
    inverse = 17;
    BOUNCEBACK();

    c[0] = 0;
    c[1] = -1;
    c[2] = 1;
    weight = 1. / 36.;
    population = 17;
    inverse = 18;
    BOUNCEBACK();

    atomicAdd(&lb_boundary_force[3 * (n_curr.boundary[index] - 1) + 0],
              boundary_force[0]);
    atomicAdd(&lb_boundary_force[3 * (n_curr.boundary[index] - 1) + 1],
              boundary_force[1]);
    atomicAdd(&lb_boundary_force[3 * (n_curr.boundary[index] - 1) + 2],
              boundary_force[2]);
  }
}

/** Add external forces within the modespace, needed for particle-interaction
 *  @param[in]     index   Node index / thread index
 *  @param[in,out] mode    Local register values mode
 *  @param[in,out] node_f  Local node force
 *  @param[in]     d_v     Local device values
 */
__device__ void apply_forces(unsigned int index, float *mode,
                             LB_node_force_density_gpu node_f,
                             LB_rho_v_gpu *d_v) {
  float u[3] = {0.0f, 0.0f, 0.0f}, C[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
  /* Note: the values d_v were calculated in relax_modes() */

  u[0] = d_v[index].v[0];
  u[1] = d_v[index].v[1];
  u[2] = d_v[index].v[2];

  C[0] += (1.0f + para->gamma_bulk) * u[0] *
              node_f.force_density[0 * para->number_of_nodes + index] +
          1.0f / 3.0f * (para->gamma_bulk - para->gamma_shear) *
              (u[0] * node_f.force_density[0 * para->number_of_nodes + index] +
               u[1] * node_f.force_density[1 * para->number_of_nodes + index] +
               u[2] * node_f.force_density[2 * para->number_of_nodes + index]);

  C[2] += (1.0f + para->gamma_bulk) * u[1] *
              node_f.force_density[1 * para->number_of_nodes + index] +
          1.0f / 3.0f * (para->gamma_bulk - para->gamma_shear) *
              (u[0] * node_f.force_density[0 * para->number_of_nodes + index] +
               u[1] * node_f.force_density[1 * para->number_of_nodes + index] +
               u[2] * node_f.force_density[2 * para->number_of_nodes + index]);

  C[5] += (1.0f + para->gamma_bulk) * u[2] *
              node_f.force_density[2 * para->number_of_nodes + index] +
          1.0f / 3.0f * (para->gamma_bulk - para->gamma_shear) *
              (u[0] * node_f.force_density[0 * para->number_of_nodes + index] +
               u[1] * node_f.force_density[1 * para->number_of_nodes + index] +
               u[2] * node_f.force_density[2 * para->number_of_nodes + index]);

  C[1] += 1.0f / 2.0f * (1.0f + para->gamma_shear) *
          (u[0] * node_f.force_density[1 * para->number_of_nodes + index] +
           u[1] * node_f.force_density[0 * para->number_of_nodes + index]);

  C[3] += 1.0f / 2.0f * (1.0f + para->gamma_shear) *
          (u[0] * node_f.force_density[2 * para->number_of_nodes + index] +
           u[2] * node_f.force_density[0 * para->number_of_nodes + index]);

  C[4] += 1.0f / 2.0f * (1.0f + para->gamma_shear) *
          (u[1] * node_f.force_density[2 * para->number_of_nodes + index] +
           u[2] * node_f.force_density[1 * para->number_of_nodes + index]);

  /* update momentum modes */
  mode[1] += node_f.force_density[0 * para->number_of_nodes + index];
  mode[2] += node_f.force_density[1 * para->number_of_nodes + index];
  mode[3] += node_f.force_density[2 * para->number_of_nodes + index];

  /* update stress modes */
  mode[4] += C[0] + C[2] + C[5];
  mode[5] += C[0] - C[2];
  mode[6] += C[0] + C[2] - 2.0f * C[5];
  mode[7] += C[1];
  mode[8] += C[3];
  mode[9] += C[4];

  reset_LB_force_densities(index, node_f);
}

/** Calculate hydrodynamic fields in MD units
 *  @param[in]  n_a     Local node residing in array a for boundary flag
 *  @param[out] mode    Local register values mode
 *  @param[out] d_p_v   Local print values
 *  @param[out] d_v     Local device values
 *  @param[in]  node_f  Local node force
 *  @param[in]  index   Node index / thread index
 *  @param[in]  print_index  Node index / thread index
 */
__device__ void
calc_values_in_MD_units(LB_nodes_gpu n_a, float *mode, LB_rho_v_pi_gpu *d_p_v,
                        LB_rho_v_gpu *d_v, LB_node_force_density_gpu node_f,
                        unsigned int index, unsigned int print_index) {
  float j[3];
  float modes_from_pi_eq[6];
  float pi[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

  if (n_a.boundary[index] == 0) {
    /* Ensure we are working with the current values of d_v */

    update_rho_v(mode, index, node_f, d_v);

    d_p_v[print_index].rho =
        d_v[index].rho / para->agrid / para->agrid / para->agrid;

    d_p_v[print_index].v[0] = d_v[index].v[0] * para->agrid / para->tau;
    d_p_v[print_index].v[1] = d_v[index].v[1] * para->agrid / para->tau;
    d_p_v[print_index].v[2] = d_v[index].v[2] * para->agrid / para->tau;
    /* stress calculation */
    float Rho = d_v[index].rho;

    /* note that d_v[index].v[] already includes the 1/2 f term, accounting
     * for the pre- and post-collisional average
     */

    j[0] = Rho * d_v[index].v[0];
    j[1] = Rho * d_v[index].v[1];
    j[2] = Rho * d_v[index].v[2];

    // equilibrium part of the stress modes, which comes from
    // the equality between modes and stress tensor components

    /* m4 = trace(pi) - rho
       m5 = pi_xx - pi_yy
       m6 = trace(pi) - 3 pi_zz
       m7 = pi_xy
       m8 = pi_xz
       m9 = pi_yz */

    // and plugging in the Euler stress for the equilibrium:
    // pi_eq = rho_0*c_s^2*I3 + (j \otimes j)/rho
    // with I3 the 3D identity matrix and
    // rho = \trace(rho_0*c_s^2*I3), which yields

    /* m4_from_pi_eq = j.j
       m5_from_pi_eq = j_x*j_x - j_y*j_y
       m6_from_pi_eq = j.j - 3*j_z*j_z
       m7_from_pi_eq = j_x*j_y
       m8_from_pi_eq = j_x*j_z
       m9_from_pi_eq = j_y*j_z */

    // where the / Rho term has been dropped. We thus obtain:

    modes_from_pi_eq[0] = (j[0] * j[0] + j[1] * j[1] + j[2] * j[2]) / Rho;
    modes_from_pi_eq[1] = (j[0] * j[0] - j[1] * j[1]) / Rho;
    modes_from_pi_eq[2] =
        (j[0] * j[0] + j[1] * j[1] + j[2] * j[2] - 3.0f * j[2] * j[2]) / Rho;
    modes_from_pi_eq[3] = j[0] * j[1] / Rho;
    modes_from_pi_eq[4] = j[0] * j[2] / Rho;
    modes_from_pi_eq[5] = j[1] * j[2] / Rho;

    /* Now we must predict the outcome of the next collision */
    /* We immediately average pre- and post-collision.  */
    /* TODO: need a reference for this.   */

    mode[4] = modes_from_pi_eq[0] + (0.5f + 0.5f * para->gamma_bulk) *
                                        (mode[4] - modes_from_pi_eq[0]);
    mode[5] = modes_from_pi_eq[1] + (0.5f + 0.5f * para->gamma_shear) *
                                        (mode[5] - modes_from_pi_eq[1]);
    mode[6] = modes_from_pi_eq[2] + (0.5f + 0.5f * para->gamma_shear) *
                                        (mode[6] - modes_from_pi_eq[2]);
    mode[7] = modes_from_pi_eq[3] + (0.5f + 0.5f * para->gamma_shear) *
                                        (mode[7] - modes_from_pi_eq[3]);
    mode[8] = modes_from_pi_eq[4] + (0.5f + 0.5f * para->gamma_shear) *
                                        (mode[8] - modes_from_pi_eq[4]);
    mode[9] = modes_from_pi_eq[5] + (0.5f + 0.5f * para->gamma_shear) *
                                        (mode[9] - modes_from_pi_eq[5]);

    // Transform the stress tensor components according to the modes that
    // correspond to those used by U. Schiller. In terms of populations this
    // expression then corresponds exactly to those in Eqs. 116 - 121 in the
    // Duenweg and Ladd paper, when these are written out in populations.
    // But to ensure this, the expression in Schiller's modes has to be
    // different!

    pi[0] +=
        (2.0f * (mode[0] + mode[4]) + mode[6] + 3.0f * mode[5]) / 6.0f; // xx
    pi[1] += mode[7];                                                   // xy
    pi[2] +=
        (2.0f * (mode[0] + mode[4]) + mode[6] - 3.0f * mode[5]) / 6.0f; // yy
    pi[3] += mode[8];                                                   // xz
    pi[4] += mode[9];                                                   // yz
    pi[5] += (mode[0] + mode[4] - mode[6]) / 3.0f;                      // zz

    for (int i = 0; i < 6; i++) {
      d_p_v[print_index].pi[i] = pi[i] / para->tau / para->tau / para->agrid;
    }
  } else {
    d_p_v[print_index].rho = 0.0f;

    for (int i = 0; i < 3; i++)
      d_p_v[print_index].v[i] = 0.0f;

    for (int i = 0; i < 6; i++)
      d_p_v[print_index].pi[i] = 0.0f;
  }
}

/** Calculate hydrodynamic fields in MD units
 *  @param[out] mode_single   Local register values mode
 *  @param[in]  d_v_single    Local device values
 *  @param[out] rho_out       Density
 *  @param[out] j_out         Momentum
 *  @param[out] pi_out        Pressure tensor
 */
__device__ void calc_values_from_m_in_LB_units(float *mode_single,
                                               LB_rho_v_gpu *d_v_single,
                                               float *rho_out, float *j_out,
                                               float *pi_out) {
  float modes_from_pi_eq[6];
  float j[6];
  float Rho;

  // stress calculation

  // Set the rho output value

  Rho = d_v_single->rho;
  *rho_out = d_v_single->rho;

  // note that d_v_single->v[] already includes the 1/2 f term,
  // accounting for the pre- and post-collisional average

  j[0] = Rho * d_v_single->v[0];
  j[1] = Rho * d_v_single->v[1];
  j[2] = Rho * d_v_single->v[2];

  j_out[3] = j[0];
  j_out[3] = j[1];
  j_out[3] = j[2];

  // equilibrium part of the stress modes, which comes from
  // the equality between modes and stress tensor components

  modes_from_pi_eq[0] = (j[0] * j[0] + j[1] * j[1] + j[2] * j[2]) / Rho;
  modes_from_pi_eq[1] = (j[0] * j[0] - j[1] * j[1]) / Rho;
  modes_from_pi_eq[2] =
      (j[0] * j[0] + j[1] * j[1] + j[2] * j[2] - 3.0f * j[2] * j[2]) / Rho;
  modes_from_pi_eq[3] = j[0] * j[1] / Rho;
  modes_from_pi_eq[4] = j[0] * j[2] / Rho;
  modes_from_pi_eq[5] = j[1] * j[2] / Rho;

  // Now we must predict the outcome of the next collision
  // We immediately average pre- and post-collision.

  mode_single[4] =
      modes_from_pi_eq[0] +
      (0.5f + 0.5f * para->gamma_bulk) * (mode_single[4] - modes_from_pi_eq[0]);
  mode_single[5] =
      modes_from_pi_eq[1] + (0.5f + 0.5f * para->gamma_shear) *
                                (mode_single[5] - modes_from_pi_eq[1]);
  mode_single[6] =
      modes_from_pi_eq[2] + (0.5f + 0.5f * para->gamma_shear) *
                                (mode_single[6] - modes_from_pi_eq[2]);
  mode_single[7] =
      modes_from_pi_eq[3] + (0.5f + 0.5f * para->gamma_shear) *
                                (mode_single[7] - modes_from_pi_eq[3]);
  mode_single[8] =
      modes_from_pi_eq[4] + (0.5f + 0.5f * para->gamma_shear) *
                                (mode_single[8] - modes_from_pi_eq[4]);
  mode_single[9] =
      modes_from_pi_eq[5] + (0.5f + 0.5f * para->gamma_shear) *
                                (mode_single[9] - modes_from_pi_eq[5]);

  // Transform the stress tensor components according to the mode_singles.

  pi_out[0] = (2.0f * (mode_single[0] + mode_single[4]) + mode_single[6] +
               3.0f * mode_single[5]) /
              6.0f;           // xx
  pi_out[1] = mode_single[7]; // xy
  pi_out[2] = (2.0f * (mode_single[0] + mode_single[4]) + mode_single[6] -
               3.0f * mode_single[5]) /
              6.0f;                                                      // yy
  pi_out[3] = mode_single[8];                                            // xz
  pi_out[4] = mode_single[9];                                            // yz
  pi_out[5] = (mode_single[0] + mode_single[4] - mode_single[6]) / 3.0f; // zz
}

/**
 *  @param[in]  node_index        Node index around (8) particle
 *  @param[out] mode              Local register values mode
 *  @param[in]  n_a               Local node residing in array a
 */
__device__ void calc_mode(float *mode, LB_nodes_gpu n_a,
                          unsigned int node_index) {
  /* mass mode */
  mode[0] = n_a.vd[0 * para->number_of_nodes + node_index] +
            n_a.vd[1 * para->number_of_nodes + node_index] +
            n_a.vd[2 * para->number_of_nodes + node_index] +
            n_a.vd[3 * para->number_of_nodes + node_index] +
            n_a.vd[4 * para->number_of_nodes + node_index] +
            n_a.vd[5 * para->number_of_nodes + node_index] +
            n_a.vd[6 * para->number_of_nodes + node_index] +
            n_a.vd[7 * para->number_of_nodes + node_index] +
            n_a.vd[8 * para->number_of_nodes + node_index] +
            n_a.vd[9 * para->number_of_nodes + node_index] +
            n_a.vd[10 * para->number_of_nodes + node_index] +
            n_a.vd[11 * para->number_of_nodes + node_index] +
            n_a.vd[12 * para->number_of_nodes + node_index] +
            n_a.vd[13 * para->number_of_nodes + node_index] +
            n_a.vd[14 * para->number_of_nodes + node_index] +
            n_a.vd[15 * para->number_of_nodes + node_index] +
            n_a.vd[16 * para->number_of_nodes + node_index] +
            n_a.vd[17 * para->number_of_nodes + node_index] +
            n_a.vd[18 * para->number_of_nodes + node_index];

  /* momentum modes */
  mode[1] = (n_a.vd[1 * para->number_of_nodes + node_index] -
             n_a.vd[2 * para->number_of_nodes + node_index]) +
            (n_a.vd[7 * para->number_of_nodes + node_index] -
             n_a.vd[8 * para->number_of_nodes + node_index]) +
            (n_a.vd[9 * para->number_of_nodes + node_index] -
             n_a.vd[10 * para->number_of_nodes + node_index]) +
            (n_a.vd[11 * para->number_of_nodes + node_index] -
             n_a.vd[12 * para->number_of_nodes + node_index]) +
            (n_a.vd[13 * para->number_of_nodes + node_index] -
             n_a.vd[14 * para->number_of_nodes + node_index]);

  mode[2] = (n_a.vd[3 * para->number_of_nodes + node_index] -
             n_a.vd[4 * para->number_of_nodes + node_index]) +
            (n_a.vd[7 * para->number_of_nodes + node_index] -
             n_a.vd[8 * para->number_of_nodes + node_index]) -
            (n_a.vd[9 * para->number_of_nodes + node_index] -
             n_a.vd[10 * para->number_of_nodes + node_index]) +
            (n_a.vd[15 * para->number_of_nodes + node_index] -
             n_a.vd[16 * para->number_of_nodes + node_index]) +
            (n_a.vd[17 * para->number_of_nodes + node_index] -
             n_a.vd[18 * para->number_of_nodes + node_index]);

  mode[3] = (n_a.vd[5 * para->number_of_nodes + node_index] -
             n_a.vd[6 * para->number_of_nodes + node_index]) +
            (n_a.vd[11 * para->number_of_nodes + node_index] -
             n_a.vd[12 * para->number_of_nodes + node_index]) -
            (n_a.vd[13 * para->number_of_nodes + node_index] -
             n_a.vd[14 * para->number_of_nodes + node_index]) +
            (n_a.vd[15 * para->number_of_nodes + node_index] -
             n_a.vd[16 * para->number_of_nodes + node_index]) -
            (n_a.vd[17 * para->number_of_nodes + node_index] -
             n_a.vd[18 * para->number_of_nodes + node_index]);
}

/**
 *  (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 *  @param[in]  n_a                Local node residing in array a
 *  @param[out] delta              Weighting of particle position
 *  @param[in]  particle_position  Particle position and velocity
 *  @param[out] node_index         Node index around (8) particle
 *  @param[in]  d_v                Local device values
 *  @param[out] interpolated_u     Interpolated velocity
 */
__device__ __inline__ void
interpolation_three_point_coupling(LB_nodes_gpu n_a, float *particle_position,
                                   unsigned int *node_index, LB_rho_v_gpu *d_v,
                                   float *delta, float *interpolated_u) {
  int my_center[3];
  float temp_delta[27];
  float mode[19];

  /* see Duenweg and Ladd http://arxiv.org/abs/0803.2826 eqn. 301 */
  /* the i index is left node, nearest node, right node */
  for (int i = 0; i < 3; ++i) {
    /* note the -0.5f is to account for the shift of the LB grid relative to
     * the MD */
    float scaledpos = particle_position[i] / para->agrid - 0.5f;
    /* the +0.5 is to turn the floorf into a round function */
    my_center[i] = (int)(floorf(scaledpos + 0.5f));
    scaledpos = scaledpos - 1.0f * my_center[i];
    temp_delta[0 + 3 * i] = (5.0f - 3.0f * abs(scaledpos + 1.0f) -
                             sqrtf(-2.0f + 6.0f * abs(scaledpos + 1.0f) -
                                   3.0f * powf(scaledpos + 1.0f, 2))) /
                            6.0f;
    temp_delta[1 + 3 * i] =
        (1.0f + sqrtf(1.0f - 3.0f * powf(scaledpos, 2))) / 3.0f;
    temp_delta[2 + 3 * i] = (5.0f - 3.0f * abs(scaledpos - 1.0f) -
                             sqrtf(-2.0f + 6.0f * abs(scaledpos - 1.0f) -
                                   3.0f * powf(scaledpos - 1.0f, 2))) /
                            6.0f;
  }

  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        delta[i + 3 * j + 9 * k + 13] =
            temp_delta[i + 1] * temp_delta[3 + j + 1] * temp_delta[6 + k + 1];
      }
    }
  }

  // modulo for negative numbers is strange at best, shift to make sure we are
  // positive
  int x = my_center[0] + para->dim_x;
  int y = my_center[1] + para->dim_y;
  int z = my_center[2] + para->dim_z;
  /* Here we collect the nodes for the three point coupling scheme (27 nodes in
   * 3d) with the analogous numbering scheme of the two point coupling scheme */
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        node_index[i + 3 * j + 9 * k + 13] =
            (x + i + para->dim_x) % para->dim_x +
            para->dim_x * ((y + j + para->dim_y) % para->dim_y) +
            para->dim_x * para->dim_y * ((z + k + para->dim_z) % para->dim_z);
      }
    }
  }

  interpolated_u[0] = 0.0f;
  interpolated_u[1] = 0.0f;
  interpolated_u[2] = 0.0f;
#pragma unroll
  for (int i = 0; i < 27; ++i) {
    float totmass = 0.0f;
    calc_m_from_n(n_a, node_index[i], mode);
    totmass += mode[0] + para->rho;
    /* The boolean expression (n_a.boundary[node_index[i]] == 0) causes boundary
       nodes to couple with velocity 0 to particles. This is necessary, since
       boundary nodes undergo the same LB dynamics as fluid nodes do. The flow
       within the boundaries does not interact with the physical fluid, since
       these populations are overwritten by the bounce back kernel. Particles
       close to walls can couple to this unphysical flow, though.
    */
    interpolated_u[0] +=
        (mode[1] / totmass) * delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[1] +=
        (mode[2] / totmass) * delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] +=
        (mode[3] / totmass) * delta[i] * (n_a.boundary[node_index[i]] == 0);
  }
}

/**
 * (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 * @param[in]  n_a                Local node residing in array a
 * @param[out] delta              Weighting of particle position
 * @param[out] delta_j            Weighting of particle momentum
 * @param[in,out] particle_data   Particle position and velocity
 * @param[in,out] particle_force  Particle force
 * @param[in]  part_index         Particle id / thread id
 * @param[out] node_index         Node index around (8) particle
 * @param[in]  d_v                Local device values
 * @param[in]  flag_cs            Determine if we are at the centre (0,
 *                                typical) or at the source (1, swimmer only)
 * @param[in]  philox_counter
 * @param[in]  friction           Friction constant for the particle coupling
 */
__device__ void calc_viscous_force_three_point_couple(
    LB_nodes_gpu n_a, float *delta, CUDA_particle_data *particle_data,
    float *particle_force, unsigned int part_index, float *delta_j,
    unsigned int *node_index, LB_rho_v_gpu *d_v, int flag_cs,
    uint64_t philox_counter, float friction) {
  float interpolated_u[3];
  float interpolated_rho;
  float viscforce_density[3];

  // Zero out workspace
#pragma unroll
  for (int jj = 0; jj < 3; ++jj) {
    viscforce_density[jj] = 0.0f;
    delta_j[jj] = 0.0f;
  }
  // Zero out only if we are at the centre of the particle <=> flag_cs = 0
  particle_force[3 * part_index + 0] =
      flag_cs * particle_force[3 * part_index + 0];
  particle_force[3 * part_index + 1] =
      flag_cs * particle_force[3 * part_index + 1];
  particle_force[3 * part_index + 2] =
      flag_cs * particle_force[3 * part_index + 2];

  float position[3];
  position[0] = particle_data[part_index].p[0];
  position[1] = particle_data[part_index].p[1];
  position[2] = particle_data[part_index].p[2];

  float velocity[3];
  velocity[0] = particle_data[part_index].v[0];
  velocity[1] = particle_data[part_index].v[1];
  velocity[2] = particle_data[part_index].v[2];

#ifdef ENGINE
  // First calculate interpolated velocity for dipole source,
  // such that we don't overwrite mode, d_v, etc. for the rest of the function
  float direction = float(particle_data[part_index].swim.push_pull) *
                    particle_data[part_index].swim.dipole_length;
  // Extrapolate position by dipole length if we are at the centre of the
  // particle
  position[0] +=
      flag_cs * direction * particle_data[part_index].swim.director[0];
  position[1] +=
      flag_cs * direction * particle_data[part_index].swim.director[1];
  position[2] +=
      flag_cs * direction * particle_data[part_index].swim.director[2];
#endif

  // Do the velocity interpolation
  interpolation_three_point_coupling(n_a, position, node_index, d_v, delta,
                                     interpolated_u);

#ifdef ENGINE
  velocity[0] -= (particle_data[part_index].swim.v_swim) *
                 particle_data[part_index].swim.director[0];
  velocity[1] -= (particle_data[part_index].swim.v_swim) *
                 particle_data[part_index].swim.director[1];
  velocity[2] -= (particle_data[part_index].swim.v_swim) *
                 particle_data[part_index].swim.director[2];

  // The first three components are v_center, the last three v_source
  // Do not use within LB, because these have already been converted back to MD
  // units
  particle_data[part_index].swim.v_cs[0 + 3 * flag_cs] =
      interpolated_u[0] * para->agrid / para->tau;
  particle_data[part_index].swim.v_cs[1 + 3 * flag_cs] =
      interpolated_u[1] * para->agrid / para->tau;
  particle_data[part_index].swim.v_cs[2 + 3 * flag_cs] =
      interpolated_u[2] * para->agrid / para->tau;
#endif

  interpolated_rho = 1.0;

  /** calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  float rhotot = 0;

  rhotot += interpolated_rho;

  /* Viscous force */
  viscforce_density[0] -=
      interpolated_rho * friction *
      (velocity[0] - interpolated_u[0] * para->agrid / para->tau) / rhotot;
  viscforce_density[1] -=
      interpolated_rho * friction *
      (velocity[1] - interpolated_u[1] * para->agrid / para->tau) / rhotot;
  viscforce_density[2] -=
      interpolated_rho * friction *
      (velocity[2] - interpolated_u[2] * para->agrid / para->tau) / rhotot;

#ifdef LB_ELECTROHYDRODYNAMICS
  viscforce_density[0] +=
      interpolated_rho * friction * particle_data[part_index].mu_E[0] / rhotot;
  viscforce_density[1] +=
      interpolated_rho * friction * particle_data[part_index].mu_E[1] / rhotot;
  viscforce_density[2] +=
      interpolated_rho * friction * particle_data[part_index].mu_E[2] / rhotot;
#endif

  /** add stochastic force of zero mean (Ahlrichs, Duenweg equ. 15)*/
  float4 random_floats = random_wrapper_philox(
      particle_data[part_index].identity, LBQ * 32, philox_counter);
  float lb_coupl_pref =
      sqrtf(12.f * 2.f * friction * para->kT / para->time_step);
  viscforce_density[0] += lb_coupl_pref * (random_floats.w - 0.5f);
  viscforce_density[1] += lb_coupl_pref * (random_floats.x - 0.5f);
  viscforce_density[2] += lb_coupl_pref * (random_floats.y - 0.5f);
  /** delta_j for transform momentum transfer to lattice units which is done
    in calc_node_force (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225
    (1999)) */
  // only add to particle_force for particle centre <=> (1-flag_cs) = 1
  particle_force[3 * part_index + 0] += (1 - flag_cs) * viscforce_density[0];
  particle_force[3 * part_index + 1] += (1 - flag_cs) * viscforce_density[1];
  particle_force[3 * part_index + 2] += (1 - flag_cs) * viscforce_density[2];

  // only add to particle_force for particle centre <=> (1-flag_cs) = 1
  delta_j[0] -= (1 - flag_cs) * viscforce_density[0] * para->time_step *
                para->tau / para->agrid;
  delta_j[1] -= (1 - flag_cs) * viscforce_density[1] * para->time_step *
                para->tau / para->agrid;
  delta_j[2] -= (1 - flag_cs) * viscforce_density[2] * para->time_step *
                para->tau / para->agrid;
#ifdef ENGINE
  // add swimming force to source position
  delta_j[0] -= flag_cs * particle_data[part_index].swim.f_swim *
                particle_data[part_index].swim.director[0] * para->time_step *
                para->tau / para->agrid;
  delta_j[1] -= flag_cs * particle_data[part_index].swim.f_swim *
                particle_data[part_index].swim.director[1] * para->time_step *
                para->tau / para->agrid;
  delta_j[2] -= flag_cs * particle_data[part_index].swim.f_swim *
                particle_data[part_index].swim.director[2] * para->time_step *
                para->tau / para->agrid;
#endif
}

/** Calculate the node force caused by the particles, with atomicAdd due to
 *  avoiding race conditions
 *  (Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 *  @param[in]  delta       Weighting of particle position
 *  @param[in]  delta_j     Weighting of particle momentum
 *  @param[in]  node_index  Node index around (8) particle
 *  @param[out] node_f      Pointer to the node force
 */
__device__ void
calc_node_force_three_point_couple(float *delta, float *delta_j,
                                   unsigned int *node_index,
                                   LB_node_force_density_gpu node_f) {
  /* TODO: should the drag depend on the density?? */
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        atomicAdd(&(node_f.force_density[0 * para->number_of_nodes +
                                         node_index[i + 3 * j + 9 * k + 13]]),
                  (delta[i + 3 * j + 9 * k + 13] * delta_j[0]));
        atomicAdd(&(node_f.force_density[1 * para->number_of_nodes +
                                         node_index[i + 3 * j + 9 * k + 13]]),
                  (delta[i + 3 * j + 9 * k + 13] * delta_j[1]));
        atomicAdd(&(node_f.force_density[2 * para->number_of_nodes +
                                         node_index[i + 3 * j + 9 * k + 13]]),
                  (delta[i + 3 * j + 9 * k + 13] * delta_j[2]));
      }
    }
  }
}

/** Calculate temperature of the fluid kernel
 *  @param[out] cpu_jsquared  Result
 *  @param[in]  n_a           Local node residing in array a
 *  @param[out] number_of_non_boundary_nodes  Local node residing in array a
 */
__global__ void temperature(LB_nodes_gpu n_a, float *cpu_jsquared,
                            int *number_of_non_boundary_nodes) {
  float mode[4];
  float jsquared = 0.0f;
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes) {
    if (!n_a.boundary[index]) {
      calc_mode(mode, n_a, index);
      jsquared = mode[1] * mode[1] + mode[2] * mode[2] + mode[3] * mode[3];
      atomicAdd(cpu_jsquared, jsquared);
      atomicAdd(number_of_non_boundary_nodes, 1);
    }
  }
}

/**
 *  (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 *  @param[in]  n_a                Local node residing in array a
 *  @param[in]  particle_position  Particle position
 *  @param[out] node_index         Node index around (8) particle
 *  @param[out] mode               The 19 modes for current lattice point
 *  @param[in]  d_v                Local device values
 *  @param[out] delta              Weighting of particle position
 *  @param[out] interpolated_u     Interpolated velocity
 */
__device__ __inline__ void interpolation_two_point_coupling(
    LB_nodes_gpu n_a, float *particle_position, unsigned int *node_index,
    float *mode, LB_rho_v_gpu *d_v, float *delta, float *interpolated_u) {
  int left_node_index[3];
  float temp_delta[6];
  float temp_delta_half[6];

  // see ahlrichs + duenweg page 8227 equ (10) and (11)
#pragma unroll
  for (int i = 0; i < 3; ++i) {
    float scaledpos = particle_position[i] / para->agrid - 0.5f;
    left_node_index[i] = (int)(floorf(scaledpos));
    temp_delta[3 + i] = scaledpos - left_node_index[i];
    temp_delta[i] = 1.0f - temp_delta[3 + i];
    // further value used for interpolation of fluid velocity at part pos near
    // boundaries
    temp_delta_half[3 + i] = (scaledpos - left_node_index[i]) * 2.0f;
    temp_delta_half[i] = 2.0f - temp_delta_half[3 + i];
  }

  delta[0] = temp_delta[0] * temp_delta[1] * temp_delta[2];
  delta[1] = temp_delta[3] * temp_delta[1] * temp_delta[2];
  delta[2] = temp_delta[0] * temp_delta[4] * temp_delta[2];
  delta[3] = temp_delta[3] * temp_delta[4] * temp_delta[2];
  delta[4] = temp_delta[0] * temp_delta[1] * temp_delta[5];
  delta[5] = temp_delta[3] * temp_delta[1] * temp_delta[5];
  delta[6] = temp_delta[0] * temp_delta[4] * temp_delta[5];
  delta[7] = temp_delta[3] * temp_delta[4] * temp_delta[5];

  // modulo for negative numbers is strange at best, shift to make sure we are
  // positive
  int x = left_node_index[0] + para->dim_x;
  int y = left_node_index[1] + para->dim_y;
  int z = left_node_index[2] + para->dim_z;

  node_index[0] = x % para->dim_x + para->dim_x * (y % para->dim_y) +
                  para->dim_x * para->dim_y * (z % para->dim_z);
  node_index[1] = (x + 1) % para->dim_x + para->dim_x * (y % para->dim_y) +
                  para->dim_x * para->dim_y * (z % para->dim_z);
  node_index[2] = x % para->dim_x + para->dim_x * ((y + 1) % para->dim_y) +
                  para->dim_x * para->dim_y * (z % para->dim_z);
  node_index[3] = (x + 1) % para->dim_x +
                  para->dim_x * ((y + 1) % para->dim_y) +
                  para->dim_x * para->dim_y * (z % para->dim_z);
  node_index[4] = x % para->dim_x + para->dim_x * (y % para->dim_y) +
                  para->dim_x * para->dim_y * ((z + 1) % para->dim_z);
  node_index[5] = (x + 1) % para->dim_x + para->dim_x * (y % para->dim_y) +
                  para->dim_x * para->dim_y * ((z + 1) % para->dim_z);
  node_index[6] = x % para->dim_x + para->dim_x * ((y + 1) % para->dim_y) +
                  para->dim_x * para->dim_y * ((z + 1) % para->dim_z);
  node_index[7] = (x + 1) % para->dim_x +
                  para->dim_x * ((y + 1) % para->dim_y) +
                  para->dim_x * para->dim_y * ((z + 1) % para->dim_z);

  interpolated_u[0] = 0.0f;
  interpolated_u[1] = 0.0f;
  interpolated_u[2] = 0.0f;
#pragma unroll
  for (int i = 0; i < 8; ++i) {
    float totmass = 0.0f;

    calc_m_from_n(n_a, node_index[i], mode);

    totmass += mode[0] + para->rho;

    /* The boolean expression (n_a.boundary[node_index[i]] == 0) causes boundary
       nodes to couple with velocity 0 to particles. This is necessary, since
       boundary nodes undergo the same LB dynamics as fluid nodes do. The flow
       within the boundaries does not interact with the physical fluid, since
       these populations are overwritten by the bounce back kernel. Particles
       close to walls can couple to this unphysical flow, though.
    */
    interpolated_u[0] +=
        (mode[1] / totmass) * delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[1] +=
        (mode[2] / totmass) * delta[i] * (n_a.boundary[node_index[i]] == 0);
    interpolated_u[2] +=
        (mode[3] / totmass) * delta[i] * (n_a.boundary[node_index[i]] == 0);
  }
}

/**
 *  (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 *  @param[in]  n_a                Local node residing in array a
 *  @param[out] delta              Weighting of particle position
 *  @param[out] delta_j            Weighting of particle momentum
 *  @param[in,out] particle_data   Particle position and velocity
 *  @param[in,out] particle_force  Particle force
 *  @param[in]  part_index         Particle id / thread id
 *  @param[out] node_index         Node index around (8) particle
 *  @param[in]  d_v                Local device values
 *  @param[in]  flag_cs            Determine if we are at the centre (0,
 *                                 typical) or at the source (1, swimmer only)
 *  @param[in]  philox_counter
 *  @param[in]  friction           Friction constant for the particle coupling
 */
__device__ void calc_viscous_force(LB_nodes_gpu n_a, float *delta,
                                   CUDA_particle_data *particle_data,
                                   float *particle_force,
                                   unsigned int part_index, float *delta_j,
                                   unsigned int *node_index, LB_rho_v_gpu *d_v,
                                   int flag_cs, uint64_t philox_counter,
                                   float friction) {
  float interpolated_u[3];
  float interpolated_rho;
  float viscforce_density[3];
  float mode[19];
// Zero out workspace
#pragma unroll
  for (int jj = 0; jj < 3; ++jj) {
    viscforce_density[jj] = 0.0f;
    delta_j[jj] = 0.0f;
  }

  // Zero out only if we are at the centre of the particle <=> flag_cs = 0
  particle_force[3 * part_index + 0] =
      flag_cs * particle_force[3 * part_index + 0];
  particle_force[3 * part_index + 1] =
      flag_cs * particle_force[3 * part_index + 1];
  particle_force[3 * part_index + 2] =
      flag_cs * particle_force[3 * part_index + 2];

  float position[3];
  position[0] = particle_data[part_index].p[0];
  position[1] = particle_data[part_index].p[1];
  position[2] = particle_data[part_index].p[2];

  float velocity[3];
  velocity[0] = particle_data[part_index].v[0];
  velocity[1] = particle_data[part_index].v[1];
  velocity[2] = particle_data[part_index].v[2];

#ifdef ENGINE
  // First calculate interpolated velocity for dipole source,
  // such that we don't overwrite mode, d_v, etc. for the rest of the function
  float direction = float(particle_data[part_index].swim.push_pull) *
                    particle_data[part_index].swim.dipole_length;
  // Extrapolate position by dipole length if we are at the centre of the
  // particle
  position[0] +=
      flag_cs * direction * particle_data[part_index].swim.director[0];
  position[1] +=
      flag_cs * direction * particle_data[part_index].swim.director[1];
  position[2] +=
      flag_cs * direction * particle_data[part_index].swim.director[2];
#endif

  // Do the velocity interpolation
  interpolation_two_point_coupling(n_a, position, node_index, mode, d_v, delta,
                                   interpolated_u);

#ifdef ENGINE
  velocity[0] -= particle_data[part_index].swim.v_swim *
                 particle_data[part_index].swim.director[0];
  velocity[1] -= particle_data[part_index].swim.v_swim *
                 particle_data[part_index].swim.director[1];
  velocity[2] -= particle_data[part_index].swim.v_swim *
                 particle_data[part_index].swim.director[2];

  // The first three components are v_center, the last three v_source
  // Do not use within LB, because these have already been converted back to MD
  // units
  particle_data[part_index].swim.v_cs[0 + 3 * flag_cs] =
      interpolated_u[0] * para->agrid / para->tau;
  particle_data[part_index].swim.v_cs[1 + 3 * flag_cs] =
      interpolated_u[1] * para->agrid / para->tau;
  particle_data[part_index].swim.v_cs[2 + 3 * flag_cs] =
      interpolated_u[2] * para->agrid / para->tau;
#endif

  /* for LB we do not reweight the friction force */
  interpolated_rho = 1.0;

  /** calculate viscous force
   * take care to rescale velocities with time_step and transform to MD units
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  float rhotot = 0;

  rhotot += interpolated_rho;

  /* Viscous force */
  viscforce_density[0] -=
      interpolated_rho * friction *
      (velocity[0] - interpolated_u[0] * para->agrid / para->tau) / rhotot;
  viscforce_density[1] -=
      interpolated_rho * friction *
      (velocity[1] - interpolated_u[1] * para->agrid / para->tau) / rhotot;
  viscforce_density[2] -=
      interpolated_rho * friction *
      (velocity[2] - interpolated_u[2] * para->agrid / para->tau) / rhotot;

#ifdef LB_ELECTROHYDRODYNAMICS
  viscforce_density[0] +=
      interpolated_rho * friction * particle_data[part_index].mu_E[0] / rhotot;
  viscforce_density[1] +=
      interpolated_rho * friction * particle_data[part_index].mu_E[1] / rhotot;
  viscforce_density[2] +=
      interpolated_rho * friction * particle_data[part_index].mu_E[2] / rhotot;
#endif

  /** add stochastic force of zero mean (Ahlrichs, Duenweg equ. 15)*/
  float4 random_floats = random_wrapper_philox(
      particle_data[part_index].identity, LBQ * 32, philox_counter);
  /* lb_coupl_pref is stored in MD units (force)
   * Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
   * The factor 12 comes from the fact that we use random numbers
   * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
   * time_step comes from the discretization.
   */
  float lb_coupl_pref =
      sqrtf(12.f * 2.f * friction * para->kT / para->time_step);
  viscforce_density[0] += lb_coupl_pref * (random_floats.w - 0.5f);
  viscforce_density[1] += lb_coupl_pref * (random_floats.x - 0.5f);
  viscforce_density[2] += lb_coupl_pref * (random_floats.y - 0.5f);

  /** delta_j for transform momentum transfer to lattice units which is done
    in calc_node_force (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225
    (1999)) */

  // only add to particle_force for particle centre <=> (1-flag_cs) = 1
  particle_force[3 * part_index + 0] += (1 - flag_cs) * viscforce_density[0];
  particle_force[3 * part_index + 1] += (1 - flag_cs) * viscforce_density[1];
  particle_force[3 * part_index + 2] += (1 - flag_cs) * viscforce_density[2];

  // only add to particle_force for particle centre <=> (1-flag_cs) = 1
  delta_j[0] -= ((1 - flag_cs) * viscforce_density[0]) * para->time_step *
                para->tau / para->agrid;
  delta_j[1] -= ((1 - flag_cs) * viscforce_density[1]) * para->time_step *
                para->tau / para->agrid;
  delta_j[2] -= ((1 - flag_cs) * viscforce_density[2]) * para->time_step *
                para->tau / para->agrid;

#ifdef ENGINE
  // add swimming force to source position
  delta_j[0] -= flag_cs * particle_data[part_index].swim.f_swim *
                particle_data[part_index].swim.director[0] * para->time_step *
                para->tau / para->agrid;
  delta_j[1] -= flag_cs * particle_data[part_index].swim.f_swim *
                particle_data[part_index].swim.director[1] * para->time_step *
                para->tau / para->agrid;
  delta_j[2] -= flag_cs * particle_data[part_index].swim.f_swim *
                particle_data[part_index].swim.director[2] * para->time_step *
                para->tau / para->agrid;
#endif
}

/** Calculate the node force caused by the particles, with atomicAdd due to
 *  avoiding race conditions
 *  (Eq. (14) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
 *  @param[in]  delta              Weighting of particle position
 *  @param[in]  delta_j            Weighting of particle momentum
 *  @param[in]  node_index         Node index around (8) particle
 *  @param[out] node_f             Node force
 */
__device__ void calc_node_force(float *delta, float *delta_j,
                                unsigned int *node_index,
                                LB_node_force_density_gpu node_f) {
  /* TODO: should the drag depend on the density?? */
  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[0]]),
      (delta[0] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[0]]),
      (delta[0] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[0]]),
      (delta[0] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[1]]),
      (delta[1] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[1]]),
      (delta[1] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[1]]),
      (delta[1] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[2]]),
      (delta[2] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[2]]),
      (delta[2] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[2]]),
      (delta[2] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[3]]),
      (delta[3] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[3]]),
      (delta[3] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[3]]),
      (delta[3] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[4]]),
      (delta[4] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[4]]),
      (delta[4] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[4]]),
      (delta[4] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[5]]),
      (delta[5] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[5]]),
      (delta[5] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[5]]),
      (delta[5] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[6]]),
      (delta[6] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[6]]),
      (delta[6] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[6]]),
      (delta[6] * delta_j[2]));

  atomicAdd(
      &(node_f.force_density[(0) * para->number_of_nodes + node_index[7]]),
      (delta[7] * delta_j[0]));
  atomicAdd(
      &(node_f.force_density[(1) * para->number_of_nodes + node_index[7]]),
      (delta[7] * delta_j[1]));
  atomicAdd(
      &(node_f.force_density[(2) * para->number_of_nodes + node_index[7]]),
      (delta[7] * delta_j[2]));
}

/*********************************************************/
/** \name System setup and Kernel functions */
/*********************************************************/

/** Kernel to calculate local populations from hydrodynamic fields.
 *  The mapping is given in terms of the equilibrium distribution.
 *
 *  Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 *  Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 *  @param[out] n_a        %Lattice site
 *  @param[out] gpu_check  Additional check if GPU kernel are executed
 *  @param[out] d_v        Local device values
 *  @param[in]  node_f     Node forces
 */
__global__ void calc_n_from_rho_j_pi(LB_nodes_gpu n_a, LB_rho_v_gpu *d_v,
                                     LB_node_force_density_gpu node_f,
                                     int *gpu_check) {
  /* TODO: this can handle only a uniform density, something similar, but local,
           has to be called every time the fields are set by the user ! */
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;
  if (index < para->number_of_nodes) {
    float mode[19];

    /* default values for fields in lattice units */
    gpu_check[0] = 1;

    float Rho = para->rho;
    float v[3] = {0.0f, 0.0f, 0.0f};
    float pi[6] = {Rho * c_sound_sq, 0.0f, Rho * c_sound_sq, 0.0f, 0.0f,
                   Rho * c_sound_sq};

    float rhoc_sq = Rho * c_sound_sq;
    float avg_rho = para->rho;
    float local_rho, local_j[3], *local_pi, trace;

    local_rho = Rho;

    local_j[0] = Rho * v[0];
    local_j[1] = Rho * v[1];
    local_j[2] = Rho * v[2];

    local_pi = pi;

    // reduce the pressure tensor to the part needed here.

    local_pi[0] -= rhoc_sq;
    local_pi[2] -= rhoc_sq;
    local_pi[5] -= rhoc_sq;

    trace = local_pi[0] + local_pi[2] + local_pi[5];

    float rho_times_coeff;
    float tmp1, tmp2;

    /* update the q=0 sublattice */
    n_a.vd[(0) * para->number_of_nodes + index] =
        1.0f / 3.0f * (local_rho - avg_rho) - 1.0f / 2.0f * trace;

    /* update the q=1 sublattice */
    rho_times_coeff = 1.0f / 18.0f * (local_rho - avg_rho);

    n_a.vd[(1) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 6.0f * local_j[0] + 1.0f / 4.0f * local_pi[0] -
        1.0f / 12.0f * trace;
    n_a.vd[(2) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 6.0f * local_j[0] + 1.0f / 4.0f * local_pi[0] -
        1.0f / 12.0f * trace;
    n_a.vd[(3) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 6.0f * local_j[1] + 1.0f / 4.0f * local_pi[2] -
        1.0f / 12.0f * trace;
    n_a.vd[(4) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 6.0f * local_j[1] + 1.0f / 4.0f * local_pi[2] -
        1.0f / 12.0f * trace;
    n_a.vd[(5) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 6.0f * local_j[2] + 1.0f / 4.0f * local_pi[5] -
        1.0f / 12.0f * trace;
    n_a.vd[(6) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 6.0f * local_j[2] + 1.0f / 4.0f * local_pi[5] -
        1.0f / 12.0f * trace;

    /* update the q=2 sublattice */
    rho_times_coeff = 1.0f / 36.0f * (local_rho - avg_rho);

    tmp1 = local_pi[0] + local_pi[2];
    tmp2 = 2.0f * local_pi[1];
    n_a.vd[(7) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] + local_j[1]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(8) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] + local_j[1]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(9) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] - local_j[1]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(10) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] - local_j[1]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;

    tmp1 = local_pi[0] + local_pi[5];
    tmp2 = 2.0f * local_pi[3];

    n_a.vd[(11) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(12) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(13) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(14) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;

    tmp1 = local_pi[2] + local_pi[5];
    tmp2 = 2.0f * local_pi[4];

    n_a.vd[(15) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[1] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(16) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[1] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(17) * para->number_of_nodes + index] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[1] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(18) * para->number_of_nodes + index] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[1] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;

    calc_m_from_n(n_a, index, mode);
    update_rho_v(mode, index, node_f, d_v);
  }
}

/** Kernel to calculate local populations from hydrodynamic fields
 *  from given flow field velocities. The mapping is given in terms of
 *  the equilibrium distribution.
 *
 *  Eq. (2.15) Ladd, J. Fluid Mech. 271, 295-309 (1994)
 *  Eq. (4) in Berk Usta, Ladd and Butler, JCP 122, 094902 (2005)
 *
 *  @param[out] n_a               Current nodes array (double buffering!)
 *  @param[in]  single_nodeindex  Single node index
 *  @param[in]  velocity          Velocity
 *  @param[out] d_v               Local device values
 *  @param[in]  node_f            Node forces
 */
__global__ void set_u_from_rho_v_pi(LB_nodes_gpu n_a, int single_nodeindex,
                                    float *velocity, LB_rho_v_gpu *d_v,
                                    LB_node_force_density_gpu node_f) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index == 0) {
    float local_rho;
    float local_j[3];
    float local_pi[6];
    float trace, avg_rho;
    float rho_times_coeff;
    float tmp1, tmp2;

    float mode_for_pi[19];
    float rho_from_m;
    float j_from_m[3];
    float pi_from_m[6];

    // Calculate the modes for this node

    calc_m_from_n(n_a, single_nodeindex, mode_for_pi);

    // Reset the d_v

    update_rho_v(mode_for_pi, single_nodeindex, node_f, d_v);

    // Calculate the density, velocity, and pressure tensor
    // in LB unit for this node

    calc_values_from_m_in_LB_units(mode_for_pi, &d_v[single_nodeindex],
                                   &rho_from_m, j_from_m, pi_from_m);

    // Take LB component density and calculate the equilibrium part
    local_rho = rho_from_m;
    avg_rho = para->rho;

    // Take LB component velocity and make it a momentum

    local_j[0] = local_rho * velocity[0];
    local_j[1] = local_rho * velocity[1];
    local_j[2] = local_rho * velocity[2];
    // Take LB component pressure tensor and put in equilibrium

    local_pi[0] = pi_from_m[0];
    local_pi[1] = pi_from_m[1];
    local_pi[2] = pi_from_m[2];
    local_pi[3] = pi_from_m[3];
    local_pi[4] = pi_from_m[4];
    local_pi[5] = pi_from_m[5];

    trace = local_pi[0] + local_pi[2] + local_pi[5];

    // update the q=0 sublattice

    n_a.vd[(0) * para->number_of_nodes + single_nodeindex] =
        1.0f / 3.0f * (local_rho - avg_rho) - 1.0f / 2.0f * trace;

    // update the q=1 sublattice

    rho_times_coeff = 1.0f / 18.0f * (local_rho - avg_rho);

    n_a.vd[(1) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 6.0f * local_j[0] + 1.0f / 4.0f * local_pi[0] -
        1.0f / 12.0f * trace;
    n_a.vd[(2) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 6.0f * local_j[0] + 1.0f / 4.0f * local_pi[0] -
        1.0f / 12.0f * trace;
    n_a.vd[(3) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 6.0f * local_j[1] + 1.0f / 4.0f * local_pi[2] -
        1.0f / 12.0f * trace;
    n_a.vd[(4) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 6.0f * local_j[1] + 1.0f / 4.0f * local_pi[2] -
        1.0f / 12.0f * trace;
    n_a.vd[(5) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 6.0f * local_j[2] + 1.0f / 4.0f * local_pi[5] -
        1.0f / 12.0f * trace;
    n_a.vd[(6) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 6.0f * local_j[2] + 1.0f / 4.0f * local_pi[5] -
        1.0f / 12.0f * trace;

    // update the q=2 sublattice

    rho_times_coeff = 1.0f / 36.0f * (local_rho - avg_rho);

    tmp1 = local_pi[0] + local_pi[2];
    tmp2 = 2.0f * local_pi[1];

    n_a.vd[(7) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] + local_j[1]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(8) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] + local_j[1]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(9) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] - local_j[1]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(10) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] - local_j[1]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;

    tmp1 = local_pi[0] + local_pi[5];
    tmp2 = 2.0f * local_pi[3];

    n_a.vd[(11) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(12) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(13) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[0] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(14) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[0] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;

    tmp1 = local_pi[2] + local_pi[5];
    tmp2 = 2.0f * local_pi[4];

    n_a.vd[(15) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[1] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(16) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[1] + local_j[2]) +
        1.0f / 8.0f * (tmp1 + tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(17) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff + 1.0f / 12.0f * (local_j[1] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;
    n_a.vd[(18) * para->number_of_nodes + single_nodeindex] =
        rho_times_coeff - 1.0f / 12.0f * (local_j[1] - local_j[2]) +
        1.0f / 8.0f * (tmp1 - tmp2) - 1.0f / 24.0f * trace;

    // Calculate the modes for this node

    calc_m_from_n(n_a, single_nodeindex, mode_for_pi);

    // Update the density and velocity field for this mode

    update_rho_v(mode_for_pi, single_nodeindex, node_f, d_v);
  }
}

/** Calculate the mass of the whole fluid kernel
 *  @param[out] sum  Resulting mass
 *  @param[in]  n_a  Local node residing in array a
 */
__global__ void calc_mass(LB_nodes_gpu n_a, float *sum) {
  float mode[4];

  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes) {
    calc_mode(mode, n_a, index);
    float Rho = mode[0] + para->rho;
    atomicAdd(&(sum[0]), Rho);
  }
}

/** (Re-)initialize the node force density / set the external force
 *  density in lb units
 *  @param[out] node_f  Local node force density
 */
__global__ void reinit_node_force(LB_node_force_density_gpu node_f) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes) {
    node_f.force_density[0 * para->number_of_nodes + index] = 0.0f;
    node_f.force_density[1 * para->number_of_nodes + index] = 0.0f;
    node_f.force_density[2 * para->number_of_nodes + index] = 0.0f;
  }
}

/** Set external force on single nodes kernel
 *  @param[in]  n_extern_node_force_densities  Number of nodes
 *  @param[in]  extern_node_force_densities    External node force array
 *  @param[out] node_f                         Node force struct
 */
__global__ void init_extern_node_force_densities(
    int n_extern_node_force_densities,
    LB_extern_nodeforcedensity_gpu *extern_node_force_densities,
    LB_node_force_density_gpu node_f) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;
  float factor = powf(para->agrid, 2) * para->tau * para->tau;
  if (index < n_extern_node_force_densities) {
    node_f.force_density[0 * para->number_of_nodes +
                         extern_node_force_densities[index].index] =
        extern_node_force_densities[index].force_density[0] * factor;
    node_f.force_density[1 * para->number_of_nodes +
                         extern_node_force_densities[index].index] =
        extern_node_force_densities[index].force_density[1] * factor;
    node_f.force_density[2 * para->number_of_nodes +
                         extern_node_force_densities[index].index] =
        extern_node_force_densities[index].force_density[2] * factor;
  }
}

/** Kernel to set the local density
 *
 *  @param[out] n_a              Current nodes array (double buffering!)
 *  @param[in] single_nodeindex  Node to set the velocity for
 *  @param[in] rho               Density to set
 *  @param[in] d_v               Local modes
 */
__global__ void set_rho(LB_nodes_gpu n_a, LB_rho_v_gpu *d_v,
                        int single_nodeindex, float rho) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;
  /*Note: this sets the velocities to zero */
  if (index == 0) {
    float local_rho;

    /** default values for fields in lattice units */
    local_rho = (rho - para->rho);
    d_v[single_nodeindex].rho = rho;

    n_a.vd[0 * para->number_of_nodes + single_nodeindex] =
        1.0f / 3.0f * local_rho;
    n_a.vd[1 * para->number_of_nodes + single_nodeindex] =
        1.0f / 18.0f * local_rho;
    n_a.vd[2 * para->number_of_nodes + single_nodeindex] =
        1.0f / 18.0f * local_rho;
    n_a.vd[3 * para->number_of_nodes + single_nodeindex] =
        1.0f / 18.0f * local_rho;
    n_a.vd[4 * para->number_of_nodes + single_nodeindex] =
        1.0f / 18.0f * local_rho;
    n_a.vd[5 * para->number_of_nodes + single_nodeindex] =
        1.0f / 18.0f * local_rho;
    n_a.vd[6 * para->number_of_nodes + single_nodeindex] =
        1.0f / 18.0f * local_rho;
    n_a.vd[7 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[8 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[9 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[10 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[11 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[12 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[13 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[14 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[15 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[16 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[17 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
    n_a.vd[18 * para->number_of_nodes + single_nodeindex] =
        1.0f / 36.0f * local_rho;
  }
}

/** Set the boundary flag for all boundary nodes
 *  @param[in]  boundary_node_list    Indices of the boundary nodes
 *  @param[in]  boundary_index_list   Flag for the corresponding boundary
 *  @param[in]  number_of_boundnodes  Number of boundary nodes
 *  @param[out] n_a                   Local node residing in array a
 *  @param[out] n_b                   Local node residing in array b
 */
__global__ void init_boundaries(int *boundary_node_list,
                                int *boundary_index_list,
                                int number_of_boundnodes, LB_nodes_gpu n_a,
                                LB_nodes_gpu n_b) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < number_of_boundnodes) {
    n_a.boundary[boundary_node_list[index]] = boundary_index_list[index];
    n_b.boundary[boundary_node_list[index]] = boundary_index_list[index];
  }
}

/** Reset the boundary flag of every node
 *  @param[out] n_a   Local node residing in array a
 *  @param[out] n_b   Local node residing in array b
 */
__global__ void reset_boundaries(LB_nodes_gpu n_a, LB_nodes_gpu n_b) {
  size_t index = blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x +
                 threadIdx.x;
  if (index < para->number_of_nodes)
    n_a.boundary[index] = n_b.boundary[index] = 0;
}

/** Integration step of the LB-fluid-solver
 *  @param[in]     n_a     Local node residing in array a
 *  @param[out]    n_b     Local node residing in array b
 *  @param[in,out] d_v     Local device values
 *  @param[in,out] node_f  Local node force density
 *  @param[in]     ek_parameters_gpu  Parameters for the electrokinetics
 *  @param[in]     philox_counter
 */
__global__ void integrate(LB_nodes_gpu n_a, LB_nodes_gpu n_b, LB_rho_v_gpu *d_v,
                          LB_node_force_density_gpu node_f,
                          EK_parameters *ek_parameters_gpu,
                          unsigned int philox_counter) {
  /*every node is connected to a thread via the index*/
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;
  /*the 19 moments (modes) are only temporary register values */
  float mode[19];

  if (index < para->number_of_nodes) {
    calc_m_from_n(n_a, index, mode);
    relax_modes(mode, index, node_f, d_v);
    if (para->kT > 0.0) {
      thermalize_modes(mode, index, philox_counter);
    }
    apply_forces(index, mode, node_f, d_v);
    normalize_modes(mode);
    calc_n_from_modes_push(n_b, mode, index);
  }
}

/** Particle interaction kernel
 *  @param[in]  n_a                 Local node residing in array a
 *  @param[in,out]  particle_data   Particle position and velocity
 *  @param[in,out]  particle_force  Particle force
 *  @param[out] node_f              Local node force
 *  @param[in]  d_v                 Local device values
 *  @param[in]  couple_virtual
 *  @param[in]  philox_counter
 */
__global__ void
calc_fluid_particle_ia(LB_nodes_gpu n_a, CUDA_particle_data *particle_data,
                       float *particle_force, LB_node_force_density_gpu node_f,
                       LB_rho_v_gpu *d_v, bool couple_virtual,
                       uint64_t philox_counter, float friction) {

  unsigned int part_index = blockIdx.y * gridDim.x * blockDim.x +
                            blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[8];
  float delta[8];
  float delta_j[3];
  if (part_index < para->number_of_particles) {
#if defined(VIRTUAL_SITES)
    if (!particle_data[part_index].is_virtual || couple_virtual)
#endif
    {
      /* force acting on the particle. delta_j will be used later to compute the
       * force that acts back onto the fluid. */
      calc_viscous_force(n_a, delta, particle_data, particle_force, part_index,
                         delta_j, node_index, d_v, 0, philox_counter, friction);
      calc_node_force(delta, delta_j, node_index, node_f);

#ifdef ENGINE
      if (particle_data[part_index].swim.swimming) {
        calc_viscous_force(n_a, delta, particle_data, particle_force,
                           part_index, delta_j, node_index, d_v, 1,
                           philox_counter, friction);
        calc_node_force(delta, delta_j, node_index, node_f);
      }
#endif
    }
  }
}

/** Particle interaction kernel
 *  @param[in]     n_a             Local node residing in array a
 *  @param[in,out] particle_data   Particle position and velocity
 *  @param[in,out] particle_force  Particle force
 *  @param[out]    node_f          Local node force
 *  @param[in]     d_v             Local device values
 *  @param[in]     philox_counter
 */
__global__ void calc_fluid_particle_ia_three_point_couple(
    LB_nodes_gpu n_a, CUDA_particle_data *particle_data, float *particle_force,
    LB_node_force_density_gpu node_f, LB_rho_v_gpu *d_v,
    uint64_t philox_counter, float friction) {
  unsigned int part_index = blockIdx.y * gridDim.x * blockDim.x +
                            blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int node_index[27];
  float delta[27];
  float delta_j[3];
  if (part_index < para->number_of_particles) {
    /**force acting on the particle. delta_j will be used later to compute the
     * force that acts back onto the fluid. */
    calc_viscous_force_three_point_couple(
        n_a, delta, particle_data, particle_force, part_index, delta_j,
        node_index, d_v, 0, philox_counter, friction);
    calc_node_force_three_point_couple(delta, delta_j, node_index, node_f);

#ifdef ENGINE
    if (particle_data[part_index].swim.swimming) {
      calc_viscous_force_three_point_couple(
          n_a, delta, particle_data, particle_force, part_index, delta_j,
          node_index, d_v, 1, philox_counter, friction);
      calc_node_force_three_point_couple(delta, delta_j, node_index, node_f);
    }
#endif
  }
}

#ifdef LB_BOUNDARIES_GPU
/** Bounce back boundary kernel
 *  @param[in]  n_curr  Pointer to local node receiving the current node field
 *  @param[in]  lb_boundary_velocity  Constant velocity at the boundary,
 *                                    set by the user
 *  @param[out] lb_boundary_force     Force on the boundary nodes
 */
__global__ void apply_boundaries(LB_nodes_gpu n_curr,
                                 float *lb_boundary_velocity,
                                 float *lb_boundary_force) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes)
    bounce_back_boundaries(n_curr, index, lb_boundary_velocity,
                           lb_boundary_force);
}

#endif

/** Get physical values of the nodes (density, velocity, ...)
 *  @param[in]  n_a     Local node residing in array a
 *  @param[out] p_v     Local print values
 *  @param[out] d_v     Local device values
 *  @param[in]  node_f  Local node force
 */
__global__ void
get_mesoscopic_values_in_MD_units(LB_nodes_gpu n_a, LB_rho_v_pi_gpu *p_v,
                                  LB_rho_v_gpu *d_v,
                                  LB_node_force_density_gpu node_f) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes) {
    float mode[19];
    calc_m_from_n(n_a, index, mode);
    calc_values_in_MD_units(n_a, mode, p_v, d_v, node_f, index, index);
  }
}

/** Get boundary flags
 *  @param[in]  n_a                 Local node residing in array a
 *  @param[out] device_bound_array  Local device values
 */
__global__ void lb_get_boundaries(LB_nodes_gpu n_a,
                                  unsigned int *device_bound_array) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes)
    device_bound_array[index] = n_a.boundary[index];
}

/** Print single node values kernel
 *  @param[in]  single_nodeindex  Node index
 *  @param[out] d_p_v   Result
 *  @param[in]  n_a     Local node residing in array a
 *  @param[out] d_v     Local device values
 *  @param[in]  node_f  Local node force
 */
__global__ void lb_print_node(int single_nodeindex, LB_rho_v_pi_gpu *d_p_v,
                              LB_nodes_gpu n_a, LB_rho_v_gpu *d_v,
                              LB_node_force_density_gpu node_f) {
  float mode[19];
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index == 0) {
    calc_m_from_n(n_a, single_nodeindex, mode);

    /* the following actually copies rho and v from d_v, and calculates pi */
    calc_values_in_MD_units(n_a, mode, d_p_v, d_v, node_f, single_nodeindex, 0);
  }
}

__global__ void momentum(LB_nodes_gpu n_a, LB_rho_v_gpu *d_v,
                         LB_node_force_density_gpu node_f, float *sum) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index < para->number_of_nodes) {
    float j[3] = {0.0f, 0.0f, 0.0f};
    float mode[4];

    calc_mode(mode, n_a, index);

    j[0] += mode[1] + node_f.force_density[0 * para->number_of_nodes + index];
    j[1] += mode[2] + node_f.force_density[1 * para->number_of_nodes + index];
    j[2] += mode[3] + node_f.force_density[2 * para->number_of_nodes + index];

#ifdef LB_BOUNDARIES_GPU
    if (n_a.boundary[index])
      j[0] = j[1] = j[2] = 0.0f;
#endif

    atomicAdd(&(sum[0]), j[0]);
    atomicAdd(&(sum[1]), j[1]);
    atomicAdd(&(sum[2]), j[2]);
  }
}
__global__ void remove_momentum(LB_nodes_gpu n_a, LB_rho_v_gpu *d_v,
                                LB_node_force_density_gpu node_f, float *sum) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;
  if (index < para->number_of_nodes) {
    node_f.force_density[0 * para->number_of_nodes + index] -=
        sum[0] / para->number_of_nodes;
    node_f.force_density[1 * para->number_of_nodes + index] -=
        sum[1] / para->number_of_nodes;
    node_f.force_density[2 * para->number_of_nodes + index] -=
        sum[2] / para->number_of_nodes;
  }
}

/** Print single node boundary flag
 *  @param[in]  single_nodeindex  Node index
 *  @param[out] device_flag       Result
 *  @param[in]  n_a               Local node residing in array a
 */
__global__ void lb_get_boundary_flag(int single_nodeindex,
                                     unsigned int *device_flag,
                                     LB_nodes_gpu n_a) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

  if (index == 0)
    device_flag[0] = n_a.boundary[single_nodeindex];
}

/**********************************************************************/
/* Host functions to setup and call kernels*/
/**********************************************************************/

void lb_get_para_pointer(LB_parameters_gpu **pointeradress) {
  if (cudaGetSymbolAddress((void **)pointeradress, HIP_SYMBOL(para)) !=
      cudaSuccess) {
    fprintf(stderr,
            "Trouble getting address of LB parameters.\n"); // TODO give proper
                                                            // error message
    errexit();
  }
}

void lb_get_lbpar_pointer(LB_parameters_gpu **pointeradress) {
  *pointeradress = &lbpar_gpu;
}

void lb_get_boundary_force_pointer(float **pointeradress) {
#ifdef LB_BOUNDARIES_GPU
  *pointeradress = lb_boundary_force;
#endif
}

void lb_get_device_values_pointer(LB_rho_v_gpu **pointeradress) {
  *pointeradress = device_rho_v;
}

/** Initialization for the lb gpu fluid called from host
 *  @param lbpar_gpu   Pointer to parameters to setup the lb field
 */
void lb_init_GPU(LB_parameters_gpu *lbpar_gpu) {
#define free_realloc_and_clear(var, size)                                      \
  {                                                                            \
    if ((var) != nullptr)                                                      \
      cuda_safe_mem(cudaFree((var)));                                          \
    cuda_safe_mem(cudaMalloc((void **)&var, size));                            \
    cudaMemset(var, 0, size);                                                  \
  }

  size_of_rho_v = lbpar_gpu->number_of_nodes * sizeof(LB_rho_v_gpu);
  size_of_rho_v_pi = lbpar_gpu->number_of_nodes * sizeof(LB_rho_v_pi_gpu);

  /** Allocate structs in device memory*/
  free_realloc_and_clear(device_rho_v, size_of_rho_v);

  /* TODO: this is a almost a copy of device_rho_v; think about eliminating
   * it, and maybe pi can be added to device_rho_v in this case */
  free_realloc_and_clear(print_rho_v_pi, size_of_rho_v_pi);
  free_realloc_and_clear(nodes_a.vd,
                         lbpar_gpu->number_of_nodes * 19 * sizeof(float));
  free_realloc_and_clear(nodes_b.vd,
                         lbpar_gpu->number_of_nodes * 19 * sizeof(float));
  free_realloc_and_clear(node_f.force_density,
                         lbpar_gpu->number_of_nodes * 3 * sizeof(lbForceFloat));
#if defined(VIRTUAL_SITES_INERTIALESS_TRACERS) || defined(EK_DEBUG)
  free_realloc_and_clear(node_f.force_density_buf,
                         lbpar_gpu->number_of_nodes * 3 * sizeof(lbForceFloat));
#endif
  free_realloc_and_clear(nodes_a.boundary,
                         lbpar_gpu->number_of_nodes * sizeof(unsigned int));
  free_realloc_and_clear(nodes_b.boundary,
                         lbpar_gpu->number_of_nodes * sizeof(unsigned int));

  /*write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(para), lbpar_gpu,
                                   sizeof(LB_parameters_gpu)));

  /*check flag if lb gpu init works*/
  free_realloc_and_clear(gpu_check, sizeof(int));

  if (h_gpu_check != nullptr)
    free(h_gpu_check);

  h_gpu_check = (int *)Utils::malloc(sizeof(int));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_boundaries, dim_grid, threads_per_block, nodes_a, nodes_b);

  /* calc of velocitydensities from given parameters and initialize the
   * Node_Force array with zero */
  KERNELCALL(reinit_node_force, dim_grid, threads_per_block, (node_f));
  KERNELCALL(calc_n_from_rho_j_pi, dim_grid, threads_per_block, nodes_a,
             device_rho_v, node_f, gpu_check);

  intflag = 1;
  current_nodes = &nodes_a;
  h_gpu_check[0] = 0;
  cuda_safe_mem(
      cudaMemcpy(h_gpu_check, gpu_check, sizeof(int), cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();

  if (!h_gpu_check[0]) {
    fprintf(stderr, "initialization of lb gpu code failed! \n");
    errexit();
  }
}

/** Reinitialization for the lb gpu fluid called from host
 *  @param lbpar_gpu   Pointer to parameters to setup the lb field
 */
void lb_reinit_GPU(LB_parameters_gpu *lbpar_gpu) {
  /* write parameters in const memory */
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(para), lbpar_gpu,
                                   sizeof(LB_parameters_gpu)));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  /* calc of velocity densities from given parameters and initialize the
   * Node_Force array with zero */
  KERNELCALL(calc_n_from_rho_j_pi, dim_grid, threads_per_block, nodes_a,
             device_rho_v, node_f, gpu_check);
}

void lb_realloc_particles_GPU_leftovers(LB_parameters_gpu *lbpar_gpu) {
  // copy parameters, especially number of parts to gpu mem
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(para), lbpar_gpu,
                                   sizeof(LB_parameters_gpu)));
}

#ifdef LB_BOUNDARIES_GPU
/** Setup and call boundaries from the host
 *  @param host_n_lb_boundaries        Number of LB boundaries
 *  @param number_of_boundnodes        Number of boundnodes
 *  @param host_boundary_node_list     The indices of the boundary nodes
 *  @param host_boundary_index_list    The flag representing the corresponding
 *                                     boundary
 *  @param host_lb_boundary_velocity   The constant velocity at the boundary,
 *                                     set by the user
 */
void lb_init_boundaries_GPU(int host_n_lb_boundaries, int number_of_boundnodes,
                            int *host_boundary_node_list,
                            int *host_boundary_index_list,
                            float *host_lb_boundary_velocity) {
  if (this_node != 0)
    return;

  size_of_boundindex = number_of_boundnodes * sizeof(int);
  cuda_safe_mem(cudaMalloc((void **)&boundary_node_list, size_of_boundindex));
  cuda_safe_mem(cudaMalloc((void **)&boundary_index_list, size_of_boundindex));
  cuda_safe_mem(cudaMemcpy(boundary_index_list, host_boundary_index_list,
                           size_of_boundindex, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(boundary_node_list, host_boundary_node_list,
                           size_of_boundindex, cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMalloc((void **)&lb_boundary_force,
                           3 * host_n_lb_boundaries * sizeof(float)));
  cuda_safe_mem(cudaMalloc((void **)&lb_boundary_velocity,
                           3 * host_n_lb_boundaries * sizeof(float)));
  cuda_safe_mem(
      cudaMemcpy(lb_boundary_velocity, host_lb_boundary_velocity,
                 3 * LBBoundaries::lbboundaries.size() * sizeof(float),
                 cudaMemcpyHostToDevice));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reset_boundaries, dim_grid, threads_per_block, nodes_a, nodes_b);

  if (LBBoundaries::lbboundaries.size() == 0 && !pdb_boundary_lattice) {
    cudaDeviceSynchronize();
    return;
  }

  if (number_of_boundnodes == 0) {
    fprintf(stderr,
            "WARNING: boundary cmd executed but no boundary node found!\n");
  } else {
    int threads_per_block_bound = 64;
    int blocks_per_grid_bound_y = 4;
    int blocks_per_grid_bound_x =
        (number_of_boundnodes +
         threads_per_block_bound * blocks_per_grid_bound_y - 1) /
        (threads_per_block_bound * blocks_per_grid_bound_y);
    dim3 dim_grid_bound =
        make_uint3(blocks_per_grid_bound_x, blocks_per_grid_bound_y, 1);

    KERNELCALL(init_boundaries, dim_grid_bound, threads_per_block_bound,
               boundary_node_list, boundary_index_list, number_of_boundnodes,
               nodes_a, nodes_b);
  }

  cudaDeviceSynchronize();
}
#endif
/** Setup and call extern single node force initialization from the host
 *  @param lbpar_gpu    Host parameter struct
 */
void lb_reinit_extern_nodeforce_GPU(LB_parameters_gpu *lbpar_gpu) {
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(para), lbpar_gpu,
                                   sizeof(LB_parameters_gpu)));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu->number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(reinit_node_force, dim_grid, threads_per_block, node_f);
}
/** Setup and call extern single node force initialization from the host
 *  @param n_extern_node_force_densities     Number of nodes on which the
 *                                           external force has to be applied
 *  @param host_extern_node_force_densities  Host extern node forces
 *  @param lbpar_gpu                         Host parameter struct
 */
void lb_init_extern_nodeforcedensities_GPU(
    int n_extern_node_force_densities,
    LB_extern_nodeforcedensity_gpu *host_extern_node_force_densities,
    LB_parameters_gpu *lbpar_gpu) {

  size_of_extern_node_force_densities =
      n_extern_node_force_densities * sizeof(LB_extern_nodeforcedensity_gpu);
  cuda_safe_mem(cudaMalloc((void **)&extern_node_force_densities,
                           size_of_extern_node_force_densities));
  cuda_safe_mem(
      cudaMemcpy(extern_node_force_densities, host_extern_node_force_densities,
                 size_of_extern_node_force_densities, cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(para), lbpar_gpu,
                                   sizeof(LB_parameters_gpu)));

  int threads_per_block_exf = 64;
  int blocks_per_grid_exf_y = 4;
  int blocks_per_grid_exf_x =
      (n_extern_node_force_densities +
       threads_per_block_exf * blocks_per_grid_exf_y - 1) /
      (threads_per_block_exf * blocks_per_grid_exf_y);
  dim3 dim_grid_exf =
      make_uint3(blocks_per_grid_exf_x, blocks_per_grid_exf_y, 1);

  KERNELCALL(init_extern_node_force_densities, dim_grid_exf,
             threads_per_block_exf, n_extern_node_force_densities,
             extern_node_force_densities, node_f);
  cudaFree(extern_node_force_densities);
}

/** Setup and call particle kernel from the host */
void lb_calc_particle_lattice_ia_gpu(bool couple_virtual, double friction) {
  if (lbpar_gpu.number_of_particles) {
    /* call of the particle kernel */
    /* values for the particle kernel */
    int threads_per_block_particles = 64;
    int blocks_per_grid_particles_y = 4;
    int blocks_per_grid_particles_x =
        (lbpar_gpu.number_of_particles +
         threads_per_block_particles * blocks_per_grid_particles_y - 1) /
        (threads_per_block_particles * blocks_per_grid_particles_y);
    dim3 dim_grid_particles =
        make_uint3(blocks_per_grid_particles_x, blocks_per_grid_particles_y, 1);

    KERNELCALL(calc_fluid_particle_ia, dim_grid_particles,
               threads_per_block_particles, *current_nodes,
               gpu_get_particle_pointer(), gpu_get_particle_force_pointer(),
               node_f, device_rho_v, couple_virtual,
               rng_counter_coupling_gpu.value(), friction);
    rng_counter_coupling_gpu.increment();
  }
}

/** Setup and call kernel for getting macroscopic fluid values of all nodes
 *  @param host_values   struct to save the gpu values
 */
void lb_get_values_GPU(LB_rho_v_pi_gpu *host_values) {
  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(get_mesoscopic_values_in_MD_units, dim_grid, threads_per_block,
             *current_nodes, print_rho_v_pi, device_rho_v, node_f);
  cuda_safe_mem(cudaMemcpy(host_values, print_rho_v_pi, size_of_rho_v_pi,
                           cudaMemcpyDeviceToHost));
}

/** Get all the boundary flags for all nodes
 *  @param host_bound_array   here go the values of the boundary flag
 */
void lb_get_boundary_flags_GPU(unsigned int *host_bound_array) {
  unsigned int *device_bound_array;
  cuda_safe_mem(cudaMalloc((void **)&device_bound_array,
                           lbpar_gpu.number_of_nodes * sizeof(unsigned int)));
  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(lb_get_boundaries, dim_grid, threads_per_block, *current_nodes,
             device_bound_array);

  cuda_safe_mem(cudaMemcpy(host_bound_array, device_bound_array,
                           lbpar_gpu.number_of_nodes * sizeof(unsigned int),
                           cudaMemcpyDeviceToHost));

  cudaFree(device_bound_array);
}

/** Setup and call kernel for getting macroscopic fluid values of a single
 *  node
 */
void lb_print_node_GPU(int single_nodeindex,
                       LB_rho_v_pi_gpu *host_print_values) {
  LB_rho_v_pi_gpu *device_print_values;
  cuda_safe_mem(
      cudaMalloc((void **)&device_print_values, sizeof(LB_rho_v_pi_gpu)));
  int threads_per_block_print = 1;
  int blocks_per_grid_print_y = 1;
  int blocks_per_grid_print_x = 1;
  dim3 dim_grid_print =
      make_uint3(blocks_per_grid_print_x, blocks_per_grid_print_y, 1);

  KERNELCALL(lb_print_node, dim_grid_print, threads_per_block_print,
             single_nodeindex, device_print_values, *current_nodes,
             device_rho_v, node_f);

  cuda_safe_mem(cudaMemcpy(host_print_values, device_print_values,
                           sizeof(LB_rho_v_pi_gpu), cudaMemcpyDeviceToHost));
  cudaFree(device_print_values);
}

/** Setup and call kernel to calculate the total momentum of the hole fluid
 *  @param mass   value of the mass calculated on the GPU
 */
void lb_calc_fluid_mass_GPU(double *mass) {
  float *tot_mass;
  float cpu_mass = 0.0f;
  cuda_safe_mem(cudaMalloc((void **)&tot_mass, sizeof(float)));
  cuda_safe_mem(
      cudaMemcpy(tot_mass, &cpu_mass, sizeof(float), cudaMemcpyHostToDevice));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(calc_mass, dim_grid, threads_per_block, *current_nodes, tot_mass);

  cuda_safe_mem(
      cudaMemcpy(&cpu_mass, tot_mass, sizeof(float), cudaMemcpyDeviceToHost));

  cudaFree(tot_mass);
  mass[0] = (double)(cpu_mass);
}

/** Setup and call kernel to calculate the total momentum of the whole fluid
 *  @param host_mom   value of the momentum calculated on the GPU
 */
void lb_calc_fluid_momentum_GPU(double *host_mom) {
  float *tot_momentum;
  float host_momentum[3] = {0.0f, 0.0f, 0.0f};
  cuda_safe_mem(cudaMalloc((void **)&tot_momentum, 3 * sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_momentum, host_momentum, 3 * sizeof(float),
                           cudaMemcpyHostToDevice));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(momentum, dim_grid, threads_per_block, *current_nodes,
             device_rho_v, node_f, tot_momentum);

  cuda_safe_mem(cudaMemcpy(host_momentum, tot_momentum, 3 * sizeof(float),
                           cudaMemcpyDeviceToHost));

  cudaFree(tot_momentum);
  host_mom[0] = (double)(host_momentum[0] * lbpar_gpu.agrid / lbpar_gpu.tau);
  host_mom[1] = (double)(host_momentum[1] * lbpar_gpu.agrid / lbpar_gpu.tau);
  host_mom[2] = (double)(host_momentum[2] * lbpar_gpu.agrid / lbpar_gpu.tau);
}

/** Setup and call kernel to remove the net momentum of the whole fluid
 */
void lb_remove_fluid_momentum_GPU(void) {
  float *tot_momentum;
  float host_momentum[3] = {0.0f, 0.0f, 0.0f};
  cuda_safe_mem(cudaMalloc((void **)&tot_momentum, 3 * sizeof(float)));
  cuda_safe_mem(cudaMemcpy(tot_momentum, host_momentum, 3 * sizeof(float),
                           cudaMemcpyHostToDevice));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(momentum, dim_grid, threads_per_block, *current_nodes,
             device_rho_v, node_f, tot_momentum);

  cuda_safe_mem(cudaMemcpy(host_momentum, tot_momentum, 3 * sizeof(float),
                           cudaMemcpyDeviceToHost));

  KERNELCALL(remove_momentum, dim_grid, threads_per_block, *current_nodes,
             device_rho_v, node_f, tot_momentum);

  cudaFree(tot_momentum);
}

/** Setup and call kernel to calculate the temperature of the hole fluid
 *  @param host_temp   value of the temperature calculated on the GPU
 */
void lb_calc_fluid_temperature_GPU(double *host_temp) {
  int host_number_of_non_boundary_nodes = 0;
  int *device_number_of_non_boundary_nodes;
  cuda_safe_mem(
      cudaMalloc((void **)&device_number_of_non_boundary_nodes, sizeof(int)));
  cuda_safe_mem(cudaMemcpy(device_number_of_non_boundary_nodes,
                           &host_number_of_non_boundary_nodes, sizeof(int),
                           cudaMemcpyHostToDevice));

  float host_jsquared = 0.0f;
  float *device_jsquared;
  cuda_safe_mem(cudaMalloc((void **)&device_jsquared, sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_jsquared, &host_jsquared, sizeof(float),
                           cudaMemcpyHostToDevice));

  /* values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(temperature, dim_grid, threads_per_block, *current_nodes,
             device_jsquared, device_number_of_non_boundary_nodes);

  cuda_safe_mem(cudaMemcpy(&host_number_of_non_boundary_nodes,
                           device_number_of_non_boundary_nodes, sizeof(int),
                           cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(&host_jsquared, device_jsquared, sizeof(float),
                           cudaMemcpyDeviceToHost));

  *host_temp = 0;

  *host_temp +=
      (double)(host_jsquared * 1. /
               (3.0f * lbpar_gpu.rho / lbpar_gpu.agrid / lbpar_gpu.agrid /
                lbpar_gpu.agrid * host_number_of_non_boundary_nodes *
                lbpar_gpu.tau * lbpar_gpu.tau * lbpar_gpu.agrid));
}

/** Setup and call kernel for getting macroscopic fluid values of all nodes
 *  @param host_checkpoint_vd         struct to save the gpu populations
 *  @param host_checkpoint_boundary   struct to save the boundary nodes
 *  @param host_checkpoint_force      struct to save the forces on the nodes
 */
void lb_save_checkpoint_GPU(float *host_checkpoint_vd,
                            unsigned int *host_checkpoint_boundary,
                            lbForceFloat *host_checkpoint_force) {
  cuda_safe_mem(cudaMemcpy(host_checkpoint_vd, current_nodes->vd,
                           lbpar_gpu.number_of_nodes * 19 * sizeof(float),
                           cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_boundary, current_nodes->boundary,
                           lbpar_gpu.number_of_nodes * sizeof(unsigned int),
                           cudaMemcpyDeviceToHost));
  cuda_safe_mem(cudaMemcpy(host_checkpoint_force, node_f.force_density,
                           lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat),
                           cudaMemcpyDeviceToHost));
}

/** Setup and call kernel for getting macroscopic fluid values of all nodes
 *  @param host_checkpoint_vd         struct to save the GPU populations
 *  @param host_checkpoint_boundary   struct to save the boundary nodes
 *  @param host_checkpoint_force      struct to save the forces on the nodes
 *  @param host_checkpoint_philox_counter
 */
void lb_load_checkpoint_GPU(float *host_checkpoint_vd,
                            unsigned int *host_checkpoint_boundary,
                            lbForceFloat *host_checkpoint_force) {
  current_nodes = &nodes_a;
  intflag = 1;

  cuda_safe_mem(cudaMemcpy(current_nodes->vd, host_checkpoint_vd,
                           lbpar_gpu.number_of_nodes * 19 * sizeof(float),
                           cudaMemcpyHostToDevice));

  cuda_safe_mem(cudaMemcpy(current_nodes->boundary, host_checkpoint_boundary,
                           lbpar_gpu.number_of_nodes * sizeof(unsigned int),
                           cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(node_f.force_density, host_checkpoint_force,
                           lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat),
                           cudaMemcpyHostToDevice));
}

/** Setup and call kernel to get the boundary flag of a single node
 *  @param single_nodeindex   number of the node to get the flag for
 *  @param host_flag          here goes the value of the boundary flag
 */
void lb_get_boundary_flag_GPU(int single_nodeindex, unsigned int *host_flag) {
  unsigned int *device_flag;
  cuda_safe_mem(cudaMalloc((void **)&device_flag, sizeof(unsigned int)));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag =
      make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);

  KERNELCALL(lb_get_boundary_flag, dim_grid_flag, threads_per_block_flag,
             single_nodeindex, device_flag, *current_nodes);

  cuda_safe_mem(cudaMemcpy(host_flag, device_flag, sizeof(unsigned int),
                           cudaMemcpyDeviceToHost));

  cudaFree(device_flag);
}

/** Set the density at a single node
 *  @param single_nodeindex   the node to set the velocity for
 *  @param host_rho           the density to set
 */
void lb_set_node_rho_GPU(int single_nodeindex, float host_rho) {
  float *device_rho;
  cuda_safe_mem(cudaMalloc((void **)&device_rho, sizeof(float)));
  cuda_safe_mem(
      cudaMemcpy(device_rho, &host_rho, sizeof(float), cudaMemcpyHostToDevice));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag =
      make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);
  KERNELCALL(set_rho, dim_grid_flag, threads_per_block_flag, *current_nodes,
             device_rho_v, single_nodeindex, *device_rho);
  cudaFree(device_rho);
}

/** Set the net velocity at a single node
 *  @param single_nodeindex   the node to set the velocity for
 *  @param host_velocity      the velocity to set
 */
void lb_set_node_velocity_GPU(int single_nodeindex, float *host_velocity) {
  float *device_velocity;
  cuda_safe_mem(cudaMalloc((void **)&device_velocity, 3 * sizeof(float)));
  cuda_safe_mem(cudaMemcpy(device_velocity, host_velocity, 3 * sizeof(float),
                           cudaMemcpyHostToDevice));
  int threads_per_block_flag = 1;
  int blocks_per_grid_flag_y = 1;
  int blocks_per_grid_flag_x = 1;
  dim3 dim_grid_flag =
      make_uint3(blocks_per_grid_flag_x, blocks_per_grid_flag_y, 1);

  KERNELCALL(set_u_from_rho_v_pi, dim_grid_flag, threads_per_block_flag,
             *current_nodes, single_nodeindex, device_velocity, device_rho_v,
             node_f);

  cudaFree(device_velocity);
}

/** Reinitialize parameters
 *  @param lbpar_gpu   struct containing the parameters of the fluid
 */
void reinit_parameters_GPU(LB_parameters_gpu *lbpar_gpu) {
  /*write parameters in const memory*/
  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(para), lbpar_gpu,
                                   sizeof(LB_parameters_gpu)));
}

/**integration kernel for the lb gpu fluid update called from host */
void lb_integrate_GPU() {
  /** values for the kernel call */
  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
#ifdef LB_BOUNDARIES_GPU
  if (LBBoundaries::lbboundaries.size() > 0) {
    cuda_safe_mem(
        cudaMemset(lb_boundary_force, 0,
                   3 * LBBoundaries::lbboundaries.size() * sizeof(float)));
  }
#endif

  /**call of fluid step*/
  if (intflag == 1) {
    KERNELCALL(integrate, dim_grid, threads_per_block, nodes_a, nodes_b,
               device_rho_v, node_f, lb_ek_parameters_gpu,
               rng_counter_fluid_gpu.value());
    current_nodes = &nodes_b;
    intflag = 0;
  } else {
    KERNELCALL(integrate, dim_grid, threads_per_block, nodes_b, nodes_a,
               device_rho_v, node_f, lb_ek_parameters_gpu,
               rng_counter_fluid_gpu.value());
    current_nodes = &nodes_a;
    intflag = 1;
  }

#ifdef LB_BOUNDARIES_GPU
  if (LBBoundaries::lbboundaries.size() > 0) {
    KERNELCALL(apply_boundaries, dim_grid, threads_per_block, *current_nodes,
               lb_boundary_velocity, lb_boundary_force);
  }
#endif
}

void lb_gpu_get_boundary_forces(double *forces) {
#ifdef LB_BOUNDARIES_GPU
  float *temp = (float *)Utils::malloc(3 * LBBoundaries::lbboundaries.size() *
                                       sizeof(float));
  cuda_safe_mem(
      cudaMemcpy(temp, lb_boundary_force,
                 3 * LBBoundaries::lbboundaries.size() * sizeof(float),
                 cudaMemcpyDeviceToHost));

  for (int i = 0; i < 3 * LBBoundaries::lbboundaries.size(); i++) {
    forces[i] = (double)temp[i];
  }
  free(temp);
#endif
}

struct lb_lbfluid_mass_of_particle {
  __device__ float operator()(CUDA_particle_data particle) const {
#ifdef MASS
    return particle.mass;
#else
    return 1.;
#endif
  };
};

void lb_lbfluid_remove_total_momentum() {
  // calculate momentum of fluid and particles
  float total_momentum[3] = {0.0f, 0.0f, 0.0f};
  lb_lbfluid_calc_linear_momentum(total_momentum, /*include_particles*/ 1,
                                  /*include_lbfluid*/ 1);

  thrust::device_ptr<CUDA_particle_data> ptr(gpu_get_particle_pointer());
  float particles_mass = thrust::transform_reduce(
      ptr, ptr + lbpar_gpu.number_of_particles, lb_lbfluid_mass_of_particle(),
      0.0f, thrust::plus<float>());

  // lb_calc_fluid_mass_GPU has to be called with double but we don't
  // want narrowing warnings, that's why we narrow it down by hand.
  double lb_calc_fluid_mass_res;
  lb_calc_fluid_mass_GPU(&lb_calc_fluid_mass_res);
  float fluid_mass = lb_calc_fluid_mass_res;

  /* Momentum fraction of the particles */
  auto const part_frac = particles_mass / (fluid_mass + particles_mass);
  /* Momentum per particle */
  float momentum_particles[3] = {-total_momentum[0] * part_frac,
                                 -total_momentum[1] * part_frac,
                                 -total_momentum[2] * part_frac};

  auto const fluid_frac = fluid_mass / (fluid_mass + particles_mass);
  float momentum_fluid[3] = {-total_momentum[0] * fluid_frac,
                             -total_momentum[1] * fluid_frac,
                             -total_momentum[2] * fluid_frac};

  lb_lbfluid_particles_add_momentum(momentum_particles);
  lb_lbfluid_fluid_add_momentum(momentum_fluid);
}

__global__ void
lb_lbfluid_fluid_add_momentum_kernel(float momentum[3], LB_nodes_gpu n_a,
                                     LB_node_force_density_gpu node_f,
                                     LB_rho_v_gpu *d_v) {
  unsigned int index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int number_of_nodes = para->number_of_nodes;
#ifdef LB_BOUNDARIES_GPU
  number_of_nodes -= para->number_of_boundnodes;
#endif
  if (index < para->number_of_nodes) {
    if (n_a.boundary[index] == 0) {
      float force_factor = powf(para->agrid, 2) * para->tau * para->tau;
      // add force density onto each node (momentum / time_step / Volume)
      node_f.force_density[0 * para->number_of_nodes + index] +=
          momentum[0] / para->tau / (number_of_nodes * powf(para->agrid, 3)) *
          force_factor;
      node_f.force_density[1 * para->number_of_nodes + index] +=
          momentum[1] / para->tau / (number_of_nodes * powf(para->agrid, 3)) *
          force_factor;
      node_f.force_density[2 * para->number_of_nodes + index] +=
          momentum[2] / para->tau / (number_of_nodes * powf(para->agrid, 3)) *
          force_factor;
    }
  }
}

void lb_lbfluid_fluid_add_momentum(float momentum_host[3]) {
  float *momentum_device;
  cuda_safe_mem(cudaMalloc((void **)&momentum_device, 3 * sizeof(float)));
  cuda_safe_mem(cudaMemcpy(momentum_device, momentum_host, 3 * sizeof(float),
                           cudaMemcpyHostToDevice));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (lbpar_gpu.number_of_nodes + threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(lb_lbfluid_fluid_add_momentum_kernel, dim_grid, threads_per_block,
             momentum_device, *current_nodes, node_f, device_rho_v);
}

/** Set the populations of a specific node on the GPU
 *  @param[out] n_a         Local node residing in array a
 *  @param[in]  population  New population
 *  @param[in]  x           x-coordinate of node
 *  @param[in]  y           y-coordinate of node
 *  @param[in]  z           z-coordinate of node
 */
__global__ void lb_lbfluid_set_population_kernel(LB_nodes_gpu n_a,
                                                 float population[LBQ], int x,
                                                 int y, int z) {
  int xyz[3] = {x, y, z};
  int index = xyz_to_index(xyz);

  for (int i = 0; i < LBQ; ++i) {
    n_a.vd[i * para->number_of_nodes + index] = population[i];
  }
}

/** Interface to set the populations of a specific node for the GPU
 *  @param[in] xyz              Node coordinates
 *  @param[in] population_host  Population
 */
void lb_lbfluid_set_population(const Vector3i &xyz,
                               float population_host[LBQ]) {
  float *population_device;
  cuda_safe_mem(cudaMalloc((void **)&population_device, LBQ * sizeof(float)));
  cuda_safe_mem(cudaMemcpy(population_device, population_host,
                           LBQ * sizeof(float), cudaMemcpyHostToDevice));

  dim3 dim_grid = make_uint3(1, 1, 1);
  KERNELCALL(lb_lbfluid_set_population_kernel, dim_grid, 1, *current_nodes,
             population_device, xyz[0], xyz[1], xyz[2]);

  cuda_safe_mem(cudaFree(population_device));
}

/** Get the populations of a specific node on the GPU
 *  @param[in]  n_a         Local node residing in array a
 *  @param[out] population  Population
 *  @param[in]  x           x-coordinate of node
 *  @param[in]  y           y-coordinate of node
 *  @param[in]  z           z-coordinate of node
 */
__global__ void lb_lbfluid_get_population_kernel(LB_nodes_gpu n_a,
                                                 float population[LBQ], int x,
                                                 int y, int z) {
  int xyz[3] = {x, y, z};
  int index = xyz_to_index(xyz);

  for (int i = 0; i < LBQ; ++i) {
    population[i] = n_a.vd[i * para->number_of_nodes + index];
  }
}

/** Interface to get the populations of a specific node for the GPU
 *  @param[in]  xyz              Node coordinates
 *  @param[out] population_host  Population
 */
void lb_lbfluid_get_population(const Vector3i &xyz,
                               float population_host[LBQ]) {
  float *population_device;
  cuda_safe_mem(cudaMalloc((void **)&population_device, LBQ * sizeof(float)));

  dim3 dim_grid = make_uint3(1, 1, 1);
  KERNELCALL(lb_lbfluid_get_population_kernel, dim_grid, 1, *current_nodes,
             population_device, xyz[0], xyz[1], xyz[2]);

  cuda_safe_mem(cudaMemcpy(population_host, population_device,
                           LBQ * sizeof(float), cudaMemcpyDeviceToHost));

  cuda_safe_mem(cudaFree(population_device));
}

struct two_point_interpolation {
  LB_nodes_gpu current_nodes_gpu;
  LB_rho_v_gpu *d_v_gpu;
  two_point_interpolation(LB_nodes_gpu _current_nodes_gpu,
                          LB_rho_v_gpu *_d_v_gpu)
      : current_nodes_gpu(_current_nodes_gpu), d_v_gpu(_d_v_gpu){};
  __device__ float3 operator()(const float3 &position) const {
    unsigned int node_index[8];
    float delta[8];
    float u[3];
    float mode[19];
    float _position[3] = {position.x, position.y, position.z};
    interpolation_two_point_coupling(current_nodes_gpu, _position, node_index,
                                     mode, d_v_gpu, delta, u);
    return make_float3(u[0], u[1], u[2]);
  }
};

void lb_get_interpolated_velocity_gpu(double const *positions,
                                      double *velocities, int length) {
  thrust::host_vector<float3> positions_host(length);
  for (int p = 0; p < 3 * length; p += 3) {
    // Cast double coming from python to float.
    positions_host[p / 3].x = static_cast<float>(positions[p]);
    positions_host[p / 3].y = static_cast<float>(positions[p + 1]);
    positions_host[p / 3].z = static_cast<float>(positions[p + 2]);
  }
  thrust::device_vector<float3> positions_device = positions_host;
  thrust::device_vector<float3> velocities_device(length);
  thrust::transform(positions_device.begin(), positions_device.end(),
                    velocities_device.begin(),
                    two_point_interpolation(*current_nodes, device_rho_v));
  thrust::host_vector<float3> velocities_host = velocities_device;
  int index = 0;
  for (auto v : velocities_host) {
    velocities[index] =
        static_cast<double>(v.x) * lbpar_gpu.agrid / lbpar_gpu.tau;
    velocities[index + 1] =
        static_cast<double>(v.y) * lbpar_gpu.agrid / lbpar_gpu.tau;
    velocities[index + 2] =
        static_cast<double>(v.z) * lbpar_gpu.agrid / lbpar_gpu.tau;
    index += 3;
  }
}

void lb_coupling_set_rng_state_gpu(uint64_t counter) {
  rng_counter_coupling_gpu = Utils::Counter<uint64_t>(counter);
}

void lb_fluid_set_rng_state_gpu(uint64_t counter) {
  rng_counter_fluid_gpu = Utils::Counter<uint64_t>(counter);
}

uint64_t lb_coupling_get_rng_state_gpu() {
  return rng_counter_coupling_gpu.value();
}
uint64_t lb_fluid_get_rng_state_gpu() { return rng_counter_fluid_gpu.value(); }

#endif /* LB_GPU */

/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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
 *
 * Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 * Includes fluctuating LB and coupling to MD particles via frictional
 * momentum transfer.
 *
 */

#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lbgpu.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include <cinttypes>
#include <fstream>

#ifdef LB

#include "cells.hpp"
#include "communication.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lbboundaries.hpp"
#include "halo.hpp"
#include "lb-d3q19.hpp"
#include "thermostat.hpp"
#include "utils/Counter.hpp"
#include "utils/matrix_vector_product.hpp"
#include "utils/u32_to_u64.hpp"
#include "utils/uniform.hpp"
#include "virtual_sites/lb_inertialess_tracers.hpp"

#include <Random123/philox.h>
#include <boost/multi_array.hpp>

#include <cassert>
#include <cstdio>
#include <iostream>
#include <mpi.h>

#include "cuda_interface.hpp"

#ifdef ADDITIONAL_CHECKS
static void lb_check_halo_regions(const LB_Fluid &lbfluid);
#endif // ADDITIONAL_CHECKS

/** Flag indicating momentum exchange between particles and fluid */
int transfer_momentum = 0;

/** Counter for the particle coupling RNG */
namespace {
Utils::Counter<uint64_t> rng_counter_coupling;
Utils::Counter<uint64_t> rng_counter_fluid;

/*
 * @brief Salt for the RNGs
 *
 * This is to avoid correlations between the
 * noise on the particle coupling and the fluid
 * thermalization.
 */
enum class RNGSalt { PARTICLES, FLUID };
} // namespace

/** Struct holding the Lattice Boltzmann parameters */
LB_Parameters lbpar = {
    // rho
    0.0,
    // viscosity
    0.0,
    // bulk_viscosity
    -1.0,
    // agrid
    -1.0,
    // tau
    -1.0,
    // friction
    0.0,
    // ext_force_density
    {0.0, 0.0, 0.0},
    // gamma_odd
    0.,
    // gamma_even
    0.,
    // gamma_shear
    0.,
    // gamma_bulk
    0.,
    // is_TRT
    false,
    // resend_halo
    0};

/** The DnQm model to be used. */
LB_Model<> lbmodel = {::D3Q19::c,
                      ::D3Q19::coefficients,
                      ::D3Q19::w,
                      ::D3Q19::e_ki,
                      ::D3Q19::w_k,
                      ::D3Q19::c_sound_sq,
                      ::D3Q19::e_ki_transposed};

#if (!defined(FLATNOISE) && !defined(GAUSSRANDOMCUT) && !defined(GAUSSRANDOM))
#define FLATNOISE
#endif

/** The underlying lattice structure */
Lattice lblattice;

using LB_FluidData = boost::multi_array<double, 2>;
static LB_FluidData lbfluid_a;
static LB_FluidData lbfluid_b;

/** Pointer to the velocity populations of the fluid.
 * lbfluid contains pre-collision populations, lbfluid_post
 * contains post-collision */
LB_Fluid lbfluid;
LB_Fluid lbfluid_post;

/** Pointer to the hydrodynamic fields of the fluid nodes */
std::vector<LB_FluidNode> lbfields;

/** Communicator for halo exchange between processors */
HaloCommunicator update_halo_comm = {0, nullptr};

/** measures the MD time since the last fluid update */
static double fluidstep = 0.0;

/***********************************************************************/
#endif // LB

#if defined(LB) || defined(LB_GPU)

#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"

/* *********************** C Interface part
 * *************************************/
/* ******************************************************************************/

/*
 * set lattice switch on C-level
 */
int lb_set_lattice_switch(int local_lattice_switch) {
  switch (local_lattice_switch) {
  case 0:
    lattice_switch = LATTICE_OFF;
    mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
    return 0;
  case 1:
    lattice_switch = LATTICE_LB;
    mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
    return 0;
  case 2:
    lattice_switch = LATTICE_LB_GPU;
    mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
    return 0;
  default:
    return 1;
  }
}

#ifdef SHANCHEN
int lb_lbfluid_set_shanchen_coupling(double *p_coupling) {
#ifdef LB_GPU
  int ii, jj, n = 0;
  switch (LB_COMPONENTS) {
  case 1:
    lbpar_gpu.coupling[0] = static_cast<float>(p_coupling[0]);
    lbpar_gpu.coupling[1] = static_cast<float>(p_coupling[1]);
    break;
  default:
    for (ii = 0; ii < LB_COMPONENTS; ii++) {
      for (jj = ii; jj < LB_COMPONENTS; jj++) {
        lbpar_gpu.coupling[LB_COMPONENTS * ii + jj] =
            static_cast<float>(p_coupling[n]);
        lbpar_gpu.coupling[LB_COMPONENTS * jj + ii] =
            static_cast<float>(p_coupling[n]);
        n++;
      }
    }
    break;
  }
  on_lb_params_change_gpu(LBPAR_COUPLING);
#endif // LB_GPU
#ifdef LB
#error not implemented
#endif // LB
  return 0;
}

#ifdef SHANCHEN
int lb_lbfluid_set_mobility(double *p_mobility) {
  int ii;
  for (ii = 0; ii < LB_COMPONENTS - 1; ii++) {
    if (p_mobility[ii] <= 0) {
      return -1;
    }
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.mobility[ii] = static_cast<float>(p_mobility[ii]);
      on_lb_params_change_gpu(LBPAR_MOBILITY);
#endif // LB_GPU
    } else {
#ifdef LB
#error not implemented
#endif // LB
    }
  }
  return 0;
}
#endif /* SHANCHEN */

int affinity_set_params(int part_type_a, int part_type_b, double *affinity) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);
  data->affinity_on = 0;
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (affinity[ii] < 0 || affinity[ii] > 1) {
      return ES_ERROR;
    }
    data->affinity[ii] = affinity[ii];
    if (data->affinity[ii] > 0)
      data->affinity_on = 1;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif // SHANCHEN

int lb_lbfluid_set_density(double *p_dens) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (p_dens[ii] <= 0)
      return -1;
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.rho[ii] = static_cast<float>(p_dens[ii]);
      on_lb_params_change_gpu(LBPAR_DENSITY);
#endif // LB_GPU
    } else {
#ifdef LB
      lbpar.rho = p_dens[ii];
      mpi_bcast_lb_params(LBPAR_DENSITY);
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_set_visc(double *p_visc) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (p_visc[ii] <= 0)
      return -1;
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.viscosity[ii] = (float)p_visc[ii];
      on_lb_params_change_gpu(LBPAR_VISCOSITY);
#endif // LB_GPU
    } else {
#ifdef LB
      lbpar.viscosity = p_visc[ii];
      mpi_bcast_lb_params(LBPAR_VISCOSITY);
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_set_bulk_visc(double *p_bulk_visc) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (p_bulk_visc[ii] <= 0)
      return -1;
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.bulk_viscosity[ii] = (float)p_bulk_visc[ii];
      lbpar_gpu.is_TRT = false;
      on_lb_params_change_gpu(LBPAR_BULKVISC);
#endif // LB_GPU
    } else {
#ifdef LB
      lbpar.bulk_viscosity = p_bulk_visc[ii];
      lbpar.is_TRT = false;
      mpi_bcast_lb_params(LBPAR_BULKVISC);
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_set_gamma_odd(double *p_gamma_odd) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (fabs(p_gamma_odd[ii]) > 1)
      return -1;
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.gamma_odd[ii] = (float)p_gamma_odd[ii];
      lbpar_gpu.is_TRT = false;
      on_lb_params_change_gpu(0);
#endif // LB_GPU
    } else {
#ifdef LB
      lbpar.gamma_odd = *p_gamma_odd;
      lbpar.is_TRT = false;
      mpi_bcast_lb_params(0);
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_set_gamma_even(double *p_gamma_even) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (fabs(p_gamma_even[ii]) > 1)
      return -1;
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.gamma_even[ii] = (float)p_gamma_even[ii];
      lbpar_gpu.is_TRT = false;
      on_lb_params_change_gpu(0);
#endif // LB_GPU
    } else {
#ifdef LB
      lbpar.gamma_even = *p_gamma_even;
      lbpar.is_TRT = false;
      mpi_bcast_lb_params(0);
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_set_friction(double *p_friction) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (p_friction[ii] <= 0)
      return -1;
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      lbpar_gpu.friction[ii] = (float)p_friction[ii];
      on_lb_params_change_gpu(LBPAR_FRICTION);
#endif // LB_GPU
    } else {
#ifdef LB
      lbpar.friction = p_friction[ii];
      mpi_bcast_lb_params(LBPAR_FRICTION);
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_get_friction(double *p_friction) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      p_friction[ii] = (double)lbpar_gpu.friction[ii];
#endif // LB_GPU
    } else {
#ifdef LB
      p_friction[ii] = lbpar.friction;
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_set_couple_flag(int couple_flag) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    if (couple_flag != LB_COUPLE_TWO_POINT &&
        couple_flag != LB_COUPLE_THREE_POINT)
      return -1;
    lbpar_gpu.lb_couple_switch = couple_flag;
#endif // LB_GPU
  } else {
#ifdef LB
    /* Only the two point nearest neighbor coupling is present in the case of
       the cpu, so just throw an error if something else is tried */
    if (couple_flag != LB_COUPLE_TWO_POINT)
      return -1;
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_couple_flag(int *couple_flag) {
  *couple_flag = LB_COUPLE_NULL;
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *couple_flag = lbpar_gpu.lb_couple_switch;
#endif
  } else {
#ifdef LB
    *couple_flag = LB_COUPLE_TWO_POINT;
#endif
  }
  return 0;
}

int lb_lbfluid_set_agrid(double p_agrid) {
  if (p_agrid <= 0)
    return -1;
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.agrid = (float)p_agrid;

    lbpar_gpu.dim_x = static_cast<unsigned int>(rint(box_l[0] / p_agrid));
    lbpar_gpu.dim_y = static_cast<unsigned int>(rint(box_l[1] / p_agrid));
    lbpar_gpu.dim_z = static_cast<unsigned int>(rint(box_l[2] / p_agrid));
    unsigned int tmp[3];
    tmp[0] = lbpar_gpu.dim_x;
    tmp[1] = lbpar_gpu.dim_y;
    tmp[2] = lbpar_gpu.dim_z;
    /* sanity checks */
    for (int dir = 0; dir < 3; dir++) {
      /* check if box_l is compatible with lattice spacing */
      if (fabs(box_l[dir] - tmp[dir] * p_agrid) > ROUND_ERROR_PREC) {
        runtimeErrorMsg() << "Lattice spacing p_agrid= " << p_agrid
                          << " is incompatible with box_l[" << dir
                          << "]=" << box_l[dir] << ", factor=" << tmp[dir]
                          << " err= " << fabs(box_l[dir] - tmp[dir] * p_agrid);
      }
    }
    lbpar_gpu.number_of_nodes =
        lbpar_gpu.dim_x * lbpar_gpu.dim_y * lbpar_gpu.dim_z;
    on_lb_params_change_gpu(LBPAR_AGRID);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.agrid = p_agrid;
    mpi_bcast_lb_params(LBPAR_AGRID);
#endif // LB
  }
  return 0;
}

int lb_lbfluid_set_tau(double p_tau) {
  if (p_tau <= 0)
    return -1;
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.tau = static_cast<float>(p_tau);
    on_lb_params_change_gpu(0);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.tau = p_tau;
    mpi_bcast_lb_params(0);
#endif // LB
  }
  return 0;
}

#ifdef SHANCHEN
int lb_lbfluid_set_remove_momentum(void) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.remove_momentum = 1;
    on_lb_params_change_gpu(0);
#endif // LB_GPU
  } else {
#ifdef LB
    return -1;
#endif // LB
  }
  return 0;
}
#endif // SHANCHEN

int lb_lbfluid_set_ext_force_density(int component, double p_fx, double p_fy,
                                     double p_fz) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    if (lbpar_gpu.tau < 0.0)
      return 2;

    if (lbpar_gpu.rho[component] <= 0.0)
      return 3;

    /* external force density is stored in MD units */
    lbpar_gpu.ext_force_density[3 * component + 0] = static_cast<float>(p_fx);
    lbpar_gpu.ext_force_density[3 * component + 1] = static_cast<float>(p_fy);
    lbpar_gpu.ext_force_density[3 * component + 2] = static_cast<float>(p_fz);
    if (p_fx != 0 || p_fy != 0 || p_fz != 0) {
      lbpar_gpu.external_force_density = 1;
    } else {
      lbpar_gpu.external_force_density = 0;
    }
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);

#endif // LB_GPU
  } else {
#ifdef LB
    if (lbpar.tau < 0.0)
      return 2;

    if (lbpar.rho <= 0.0)
      return 3;

    lbpar.ext_force_density[0] = p_fx;
    lbpar.ext_force_density[1] = p_fy;
    lbpar.ext_force_density[2] = p_fz;
    mpi_bcast_lb_params(LBPAR_EXTFORCE);
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_density(double *p_dens) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      p_dens[ii] = static_cast<double>(lbpar_gpu.rho[ii]);
#endif // LB_GPU
    } else {
#ifdef LB
      p_dens[ii] = lbpar.rho;
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_get_visc(double *p_visc) {
  for (int ii = 0; ii < LB_COMPONENTS; ii++) {
    if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
      p_visc[ii] = static_cast<double>(lbpar_gpu.viscosity[ii]);
#endif // LB_GPU
    } else {
#ifdef LB
      p_visc[ii] = lbpar.viscosity;
#endif // LB
    }
  }
  return 0;
}

int lb_lbfluid_get_bulk_visc(double *p_bulk_visc) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_bulk_visc = lbpar_gpu.bulk_viscosity[0];
#endif // LB_GPU
  } else {
#ifdef LB
    *p_bulk_visc = lbpar.bulk_viscosity;
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_gamma_odd(double *p_gamma_odd) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_gamma_odd = lbpar_gpu.gamma_odd[0];
#endif // LB_GPU
  } else {
#ifdef LB
    *p_gamma_odd = lbpar.gamma_odd;
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_gamma_even(double *p_gamma_even) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_gamma_even = lbpar_gpu.gamma_even[0];
#endif // LB_GPU
  } else {
#ifdef LB
    *p_gamma_even = lbpar.gamma_even;
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_agrid(double *p_agrid) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_agrid = lbpar_gpu.agrid;
#endif // LB_GPU
  } else {
#ifdef LB
    *p_agrid = lbpar.agrid;
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_tau(double *p_tau) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    *p_tau = lbpar_gpu.tau;
#endif // LB_GPU
  } else {
#ifdef LB
    *p_tau = lbpar.tau;
#endif // LB
  }
  return 0;
}

int lb_lbfluid_get_ext_force_density(double *p_f) {
#ifdef SHANCHEN
  fprintf(stderr, "Not implemented yet (%s:%d) ", __FILE__, __LINE__);
  errexit();
#endif // SHANCHEN
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    p_f[0] = lbpar_gpu.ext_force_density[0];
    p_f[1] = lbpar_gpu.ext_force_density[1];
    p_f[2] = lbpar_gpu.ext_force_density[2];
#endif // LB_GPU
  } else {
#ifdef LB
    p_f[0] = lbpar.ext_force_density[0];
    p_f[1] = lbpar.ext_force_density[1];
    p_f[2] = lbpar.ext_force_density[2];
#endif // LB
  }
  return 0;
}

int lb_lbfluid_print_vtk_boundary(char *filename) {
  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    unsigned int *bound_array;
    bound_array = (unsigned int *)Utils::malloc(lbpar_gpu.number_of_nodes *
                                                sizeof(unsigned int));
    lb_get_boundary_flags_GPU(bound_array);

    int j;
    /** print of the calculated phys values */
    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbboundaries\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %u %u %u\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %u\n"
            "SCALARS boundary float 1\nLOOKUP_TABLE default\n",
            lbpar_gpu.dim_x, lbpar_gpu.dim_y, lbpar_gpu.dim_z,
            lbpar_gpu.agrid * 0.5, lbpar_gpu.agrid * 0.5, lbpar_gpu.agrid * 0.5,
            lbpar_gpu.agrid, lbpar_gpu.agrid, lbpar_gpu.agrid,
            lbpar_gpu.number_of_nodes);
    for (j = 0; j < int(lbpar_gpu.number_of_nodes); ++j) {
      /** print of the calculated phys values */
      fprintf(fp, "%d \n", bound_array[j]);
    }
    free(bound_array);
#endif // LB_GPU
  } else {
#ifdef LB
    Vector3i pos;
    int boundary;

    auto const grid_size = lblattice.global_grid;

    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbboundaries\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS boundary float 1\nLOOKUP_TABLE default\n",
            grid_size[0], grid_size[1], grid_size[2], lblattice.agrid[0] * 0.5,
            lblattice.agrid[1] * 0.5, lblattice.agrid[2] * 0.5,
            lblattice.agrid[0], lblattice.agrid[1], lblattice.agrid[2],
            grid_size[0] * grid_size[1] * grid_size[2]);

    for (pos[2] = 0; pos[2] < grid_size[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < grid_size[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < grid_size[0]; pos[0]++) {
          lb_lbnode_get_boundary(pos, &boundary);
          fprintf(fp, "%d \n", boundary);
        }
      }
    }
#endif // LB
  }
  fclose(fp);
  return 0;
}

int lb_lbfluid_print_vtk_velocity(char *filename, std::vector<int> bb1,
                                  std::vector<int> bb2) {
  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  std::vector<int> bb_low;
  std::vector<int> bb_high;

  for (auto val1 = bb1.begin(), val2 = bb2.begin();
       val1 != bb1.end() && val2 != bb2.end(); ++val1, ++val2) {
    if (*val1 == -1 || *val2 == -1) {
      bb_low = {0, 0, 0};
      if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
        bb_high = {static_cast<int>(lbpar_gpu.dim_x) - 1,
                   static_cast<int>(lbpar_gpu.dim_y) - 1,
                   static_cast<int>(lbpar_gpu.dim_z) - 1};
#endif // LB_GPU
      } else {
#ifdef LB
        bb_high = {lblattice.global_grid[0] - 1, lblattice.global_grid[1] - 1,
                   lblattice.global_grid[2] - 1};
#endif // LB
      }
      break;
    }

    bb_low.push_back(std::min(*val1, *val2));
    bb_high.push_back(std::max(*val1, *val2));
  }

  Vector3i pos;
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
    host_values = (LB_rho_v_pi_gpu *)Utils::malloc(size_of_values);
    lb_get_values_GPU(host_values);
    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbfluid_gpu\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS velocity float 3\nLOOKUP_TABLE default\n",
            bb_high[0] - bb_low[0] + 1, bb_high[1] - bb_low[1] + 1,
            bb_high[2] - bb_low[2] + 1, (bb_low[0] + 0.5) * lbpar_gpu.agrid,
            (bb_low[1] + 0.5) * lbpar_gpu.agrid,
            (bb_low[2] + 0.5) * lbpar_gpu.agrid, lbpar_gpu.agrid,
            lbpar_gpu.agrid, lbpar_gpu.agrid,
            (bb_high[0] - bb_low[0] + 1) * (bb_high[1] - bb_low[1] + 1) *
                (bb_high[2] - bb_low[2] + 1));
    for (pos[2] = bb_low[2]; pos[2] <= bb_high[2]; pos[2]++)
      for (pos[1] = bb_low[1]; pos[1] <= bb_high[1]; pos[1]++)
        for (pos[0] = bb_low[0]; pos[0] <= bb_high[0]; pos[0]++) {
          int j = lbpar_gpu.dim_y * lbpar_gpu.dim_x * pos[2] +
                  lbpar_gpu.dim_x * pos[1] + pos[0];
          fprintf(fp, "%f %f %f\n", host_values[j].v[0], host_values[j].v[1],
                  host_values[j].v[2]);
        }
    free(host_values);
#endif // LB_GPU
  } else {
#ifdef LB
    Vector3d u;

    fprintf(fp,
            "# vtk DataFile Version 2.0\nlbfluid_cpu\n"
            "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
            "ORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\n"
            "SCALARS velocity float 3\nLOOKUP_TABLE default\n",
            bb_high[0] - bb_low[0] + 1, bb_high[1] - bb_low[1] + 1,
            bb_high[2] - bb_low[2] + 1, (bb_low[0] + 0.5) * lblattice.agrid[0],
            (bb_low[1] + 0.5) * lblattice.agrid[1],
            (bb_low[2] + 0.5) * lblattice.agrid[2], lblattice.agrid[0],
            lblattice.agrid[1], lblattice.agrid[2],
            (bb_high[0] - bb_low[0] + 1) * (bb_high[1] - bb_low[1] + 1) *
                (bb_high[2] - bb_low[2] + 1));

    for (pos[2] = bb_low[2]; pos[2] <= bb_high[2]; pos[2]++)
      for (pos[1] = bb_low[1]; pos[1] <= bb_high[1]; pos[1]++)
        for (pos[0] = bb_low[0]; pos[0] <= bb_high[0]; pos[0]++) {
          lb_lbnode_get_u(pos, u.data());
          fprintf(fp, "%f %f %f\n", u[0], u[1], u[2]);
        }
#endif // LB
  }
  fclose(fp);

  return 0;
}

int lb_lbfluid_print_boundary(char *filename) {
  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    unsigned int *bound_array;
    bound_array = (unsigned int *)Utils::malloc(lbpar_gpu.number_of_nodes *
                                                sizeof(unsigned int));
    lb_get_boundary_flags_GPU(bound_array);

    Vector3i xyz;
    int j;
    for (j = 0; j < int(lbpar_gpu.number_of_nodes); ++j) {
      xyz[0] = j % lbpar_gpu.dim_x;
      int k = j / lbpar_gpu.dim_x;
      xyz[1] = k % lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;
      /** print of the calculated phys values */
      fprintf(fp, "%f %f %f %u\n", (xyz[0] + 0.5) * lbpar_gpu.agrid,
              (xyz[1] + 0.5) * lbpar_gpu.agrid,
              (xyz[2] + 0.5) * lbpar_gpu.agrid, bound_array[j]);
    }
    free(bound_array);
#endif // LB_GPU
  } else {
#ifdef LB
    Vector3i pos;
    int boundary;
    Vector3i gridsize;

    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];

    for (pos[2] = 0; pos[2] < gridsize[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < gridsize[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
          lb_lbnode_get_boundary(pos, &boundary);
          boundary = (boundary != 0 ? 1 : 0);
          fprintf(fp, "%f %f %f %d\n", (pos[0] + 0.5) * lblattice.agrid[0],
                  (pos[1] + 0.5) * lblattice.agrid[1],
                  (pos[2] + 0.5) * lblattice.agrid[2], boundary);
        }
      }
    }
#endif // LB
  }

  fclose(fp);
  return 0;
}

int lb_lbfluid_print_velocity(char *filename) {
  FILE *fp = fopen(filename, "w");

  if (fp == nullptr) {
    return 1;
  }

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
#ifdef SHANCHEN
    fprintf(stderr, "TODO:adapt for SHANCHEN (%s:%d)\n", __FILE__, __LINE__);
    errexit();
#endif // SHANCHEN
    size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
    host_values = (LB_rho_v_pi_gpu *)Utils::malloc(size_of_values);
    lb_get_values_GPU(host_values);
    Vector3i xyz;
    int j;
    for (j = 0; j < int(lbpar_gpu.number_of_nodes); ++j) {
      xyz[0] = j % lbpar_gpu.dim_x;
      int k = j / lbpar_gpu.dim_x;
      xyz[1] = k % lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;
      /** print of the calculated phys values */
      fprintf(fp, "%f %f %f %f %f %f\n", (xyz[0] + 0.5) * lbpar_gpu.agrid,
              (xyz[1] + 0.5) * lbpar_gpu.agrid,
              (xyz[2] + 0.5) * lbpar_gpu.agrid, host_values[j].v[0],
              host_values[j].v[1], host_values[j].v[2]);
    }
    free(host_values);
#endif // LB_GPU
  } else {
#ifdef LB
    Vector3i pos;
    Vector3d u;
    Vector3i gridsize;

    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];

    for (pos[2] = 0; pos[2] < gridsize[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < gridsize[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
#ifdef SHANCHEN
          fprintf(stderr, "SHANCHEN not implemented for the CPU LB\n", __FILE__,
                  __LINE__);
          errexit();
#endif // SHANCHEN
          lb_lbnode_get_u(pos, u.data());
          fprintf(fp, "%f %f %f %f %f %f\n",
                  (pos[0] + 0.5) * lblattice.agrid[0],
                  (pos[1] + 0.5) * lblattice.agrid[1],
                  (pos[2] + 0.5) * lblattice.agrid[2], u[0], u[1], u[2]);
        }
      }
    }
#endif // LB
  }

  fclose(fp);
  return 0;
}

int lb_lbfluid_save_checkpoint(char *filename, int binary) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float *host_checkpoint_vd =
        (float *)Utils::malloc(lbpar_gpu.number_of_nodes * 19 * sizeof(float));
    unsigned int *host_checkpoint_boundary = (unsigned int *)Utils::malloc(
        lbpar_gpu.number_of_nodes * sizeof(unsigned int));
    lbForceFloat *host_checkpoint_force = (lbForceFloat *)Utils::malloc(
        lbpar_gpu.number_of_nodes * 3 * sizeof(lbForceFloat));
    lb_save_checkpoint_GPU(host_checkpoint_vd, host_checkpoint_boundary,
                           host_checkpoint_force);
    if (!binary) {
      std::fstream cpfile(filename, std::ios::out);
      cpfile << std::fixed;
      cpfile.precision(8);
      for (int n = 0; n < (19 * int(lbpar_gpu.number_of_nodes)); n++) {
        cpfile << host_checkpoint_vd[n] << "\n";
      }
      for (int n = 0; n < int(lbpar_gpu.number_of_nodes); n++) {
        cpfile << host_checkpoint_boundary[n] << "\n";
      }
      for (int n = 0; n < (3 * int(lbpar_gpu.number_of_nodes)); n++) {
        cpfile << host_checkpoint_force[n] << "\n";
      }
      cpfile.close();
    } else {
      std::fstream cpfile(filename, std::ios::out | std::ios::binary);
      cpfile.write(reinterpret_cast<char *>(&host_checkpoint_vd),
                   19 * sizeof(float) * lbpar_gpu.number_of_nodes);
      cpfile.write(reinterpret_cast<char *>(&host_checkpoint_boundary),
                   sizeof(int) * lbpar_gpu.number_of_nodes);
      cpfile.write(reinterpret_cast<char *>(&host_checkpoint_force),
                   3 * sizeof(float) * lbpar_gpu.number_of_nodes);
      cpfile.close();
    }
    free(host_checkpoint_vd);
    free(host_checkpoint_boundary);
    free(host_checkpoint_force);
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    std::fstream cpfile;
    if (binary) {
      cpfile.open(filename, std::ios::out | std::ios::binary);
    } else {
      cpfile.open(filename, std::ios::out);
      cpfile.precision(16);
      cpfile << std::fixed;
    }
    double pop[19];
    Vector3i ind;

    Vector3i gridsize;

    gridsize[0] = box_l[0] / lbpar.agrid;
    gridsize[1] = box_l[1] / lbpar.agrid;
    gridsize[2] = box_l[2] / lbpar.agrid;

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          ind[0] = i;
          ind[1] = j;
          ind[2] = k;
          lb_lbnode_get_pop(ind, pop);
          if (!binary) {
            for (int n = 0; n < 19; n++) {
              cpfile << pop[n];
            }
            cpfile << "\n";
          } else {
            cpfile.write(reinterpret_cast<char *>(&pop[0]),
                         19 * sizeof(double));
          }
        }
      }
    }
    cpfile.close();
#endif // LB
  }
  return ES_OK;
}

int lb_lbfluid_load_checkpoint(char *filename, int binary) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    FILE *cpfile;
    cpfile = fopen(filename, "r");
    if (!cpfile) {
      return ES_ERROR;
    }
    std::vector<float> host_checkpoint_vd(lbpar_gpu.number_of_nodes * 19);
    std::vector<unsigned int> host_checkpoint_boundary(
        lbpar_gpu.number_of_nodes);
    std::vector<lbForceFloat> host_checkpoint_force(lbpar_gpu.number_of_nodes *
                                                    3);
    int res = EOF;
    if (!binary) {
      for (int n = 0; n < (19 * int(lbpar_gpu.number_of_nodes)); n++) {
        res = fscanf(cpfile, "%f", &host_checkpoint_vd[n]);
      }
      for (int n = 0; n < int(lbpar_gpu.number_of_nodes); n++) {
        res = fscanf(cpfile, "%u", &host_checkpoint_boundary[n]);
      }
      for (int n = 0; n < (3 * int(lbpar_gpu.number_of_nodes)); n++) {
        res = fscanf(cpfile, "%f", &host_checkpoint_force[n]);
      }
      if (lbpar_gpu.number_of_nodes && res == EOF)
        throw std::runtime_error("Error while reading LB checkpoint.");
    } else {
      if (fread(host_checkpoint_vd.data(), sizeof(float),
                19 * int(lbpar_gpu.number_of_nodes),
                cpfile) != (unsigned int)(19 * lbpar_gpu.number_of_nodes))
        return ES_ERROR;
      if (fread(host_checkpoint_boundary.data(), sizeof(int),
                int(lbpar_gpu.number_of_nodes),
                cpfile) != (unsigned int)lbpar_gpu.number_of_nodes) {
        fclose(cpfile);
        return ES_ERROR;
      }
      if (fread(host_checkpoint_force.data(), sizeof(lbForceFloat),
                3 * int(lbpar_gpu.number_of_nodes),
                cpfile) != (unsigned int)(3 * lbpar_gpu.number_of_nodes)) {
        fclose(cpfile);
        return ES_ERROR;
      }
    }
    lb_load_checkpoint_GPU(host_checkpoint_vd.data(),
                           host_checkpoint_boundary.data(),
                           host_checkpoint_force.data());
    fclose(cpfile);
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    FILE *cpfile;
    cpfile = fopen(filename, "r");
    if (!cpfile) {
      return ES_ERROR;
    }
    double pop[19];
    Vector3i ind;

    Vector3i gridsize;
    mpi_bcast_lb_params(0);
    gridsize[0] = box_l[0] / lbpar.agrid;
    gridsize[1] = box_l[1] / lbpar.agrid;
    gridsize[2] = box_l[2] / lbpar.agrid;

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          ind[0] = i;
          ind[1] = j;
          ind[2] = k;
          if (!binary) {
            if (fscanf(cpfile,
                       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                       "%lf %lf %lf %lf %lf %lf \n",
                       &pop[0], &pop[1], &pop[2], &pop[3], &pop[4], &pop[5],
                       &pop[6], &pop[7], &pop[8], &pop[9], &pop[10], &pop[11],
                       &pop[12], &pop[13], &pop[14], &pop[15], &pop[16],
                       &pop[17], &pop[18]) != 19) {
              return ES_ERROR;
            }
          } else {
            if (fread(pop, sizeof(double), 19, cpfile) != 19)
              return ES_ERROR;
          }
          lb_lbnode_set_pop(ind, pop);
        }
      }
    }
    fclose(cpfile);
#endif // LB
  } else {
    runtimeErrorMsg() << "To load an LB checkpoint one needs to have already "
                         "initialized the LB fluid with the same grid size.";
    return ES_ERROR;
  }
  return ES_OK;
}

bool lb_lbnode_is_index_valid(const Vector3i &ind) {
  auto within_bounds = [](const Vector3i &ind, const Vector3i &limits) {
    return ind < limits && ind >= Vector3i{};
  };
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return within_bounds(ind, {static_cast<int>(lbpar_gpu.dim_x),
                               static_cast<int>(lbpar_gpu.dim_y),
                               static_cast<int>(lbpar_gpu.dim_z)});
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return within_bounds(ind, lblattice.global_grid);
#endif
  }
  return false;
}

int lb_lbnode_get_rho(const Vector3i &ind, double *p_rho) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    static LB_rho_v_pi_gpu *host_print_values = nullptr;

    if (host_print_values == nullptr)
      host_print_values =
          (LB_rho_v_pi_gpu *)Utils::malloc(sizeof(LB_rho_v_pi_gpu));
    lb_print_node_GPU(single_nodeindex, host_print_values);
    for (int ii = 0; ii < LB_COMPONENTS; ii++) {
      p_rho[ii] = (double)(host_print_values->rho[ii]);
    }
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];
    double rho;
    double j[3];
    double pi[6];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j, pi);
    // unit conversion
    rho *= 1 / lbpar.agrid / lbpar.agrid / lbpar.agrid;
    *p_rho = rho;
#endif // LB
  }
  return 0;
}

int lb_lbnode_get_u(const Vector3i &ind, double *p_u) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    static LB_rho_v_pi_gpu *host_print_values = nullptr;
    if (host_print_values == nullptr)
      host_print_values =
          (LB_rho_v_pi_gpu *)Utils::malloc(sizeof(LB_rho_v_pi_gpu));

    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, host_print_values);

    p_u[0] = (double)(host_print_values[0].v[0]);
    p_u[1] = (double)(host_print_values[0].v[1]);
    p_u[2] = (double)(host_print_values[0].v[2]);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];
    double rho;
    double j[3];
    double pi[6];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j, pi);
    // unit conversion
    p_u[0] = j[0] / rho * lbpar.agrid / lbpar.tau;
    p_u[1] = j[1] / rho * lbpar.agrid / lbpar.tau;
    p_u[2] = j[2] / rho * lbpar.agrid / lbpar.tau;
#endif // LB
  }
  return 0;
}

/** calculates the fluid velocity at a given position of the
 * lattice. Note that it can lead to undefined behavior if the
 * position is not within the local lattice. This version of the function
 * can be called without the position needing to be on the local processor.
 * Note that this gives a slightly different version than the values used to
 * couple to MD beads when near a wall, see
 * lb_lbfluid_get_interpolated_velocity.
 */
int lb_lbfluid_get_interpolated_velocity_global(Vector3d &p, double *v) {
  double local_v[3] = {0, 0, 0},
         delta[6]{}; // velocity field, relative positions to surrounding nodes
  Vector3i ind = {0, 0, 0}, tmpind; // node indices
  int x, y, z;                      // counters

  // convert the position into lower left grid point
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    Lattice::map_position_to_lattice_global(p, ind, delta, lbpar_gpu.agrid);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::map_position_to_lattice_global(p, ind, delta, lbpar.agrid);
#endif // LB
  }

  // set the initial velocity to zero in all directions
  v[0] = 0;
  v[1] = 0;
  v[2] = 0;

  for (z = 0; z < 2; z++) {
    for (y = 0; y < 2; y++) {
      for (x = 0; x < 2; x++) {
        // give the index of the neighbouring nodes
        tmpind[0] = ind[0] + x;
        tmpind[1] = ind[1] + y;
        tmpind[2] = ind[2] + z;

        if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
          if (tmpind[0] == int(lbpar_gpu.dim_x))
            tmpind[0] = 0;
          if (tmpind[1] == int(lbpar_gpu.dim_y))
            tmpind[1] = 0;
          if (tmpind[2] == int(lbpar_gpu.dim_z))
            tmpind[2] = 0;
#endif // LB_GPU
        } else {
#ifdef LB
          if (tmpind[0] == box_l[0] / lbpar.agrid)
            tmpind[0] = 0;
          if (tmpind[1] == box_l[1] / lbpar.agrid)
            tmpind[1] = 0;
          if (tmpind[2] == box_l[2] / lbpar.agrid)
            tmpind[2] = 0;
#endif // LB
        }

        lb_lbnode_get_u(tmpind, local_v);

        v[0] +=
            delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] * local_v[0];
        v[1] +=
            delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] * local_v[1];
        v[2] +=
            delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] * local_v[2];
      }
    }
  }

  return 0;
}

int lb_lbnode_get_pi(const Vector3i &ind, double *p_pi) {
  double p0 = 0;

  lb_lbnode_get_pi_neq(ind, p_pi);

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    for (int ii = 0; ii < LB_COMPONENTS; ii++) {
      p0 += lbpar_gpu.rho[ii] / lbpar_gpu.agrid / lbpar_gpu.tau /
            lbpar_gpu.tau / 3.;
    }
#endif // LB_GPU
  } else {
#ifdef LB
    p0 = lbpar.rho / lbpar.agrid / lbpar.tau / lbpar.tau / 3.;
#endif // LB
  }

  p_pi[0] += p0;
  p_pi[2] += p0;
  p_pi[5] += p0;

  return 0;
}

int lb_lbnode_get_pi_neq(const Vector3i &ind, double *p_pi) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    static LB_rho_v_pi_gpu *host_print_values = nullptr;
    if (host_print_values == nullptr)
      host_print_values =
          (LB_rho_v_pi_gpu *)Utils::malloc(sizeof(LB_rho_v_pi_gpu));

    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, host_print_values);
    for (int i = 0; i < 6; i++) {
      p_pi[i] = host_print_values->pi[i];
    }
    return 0;
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];
    double rho;
    double j[3];
    double pi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j, pi);
    // unit conversion
    p_pi[0] = pi[0] / lbpar.tau / lbpar.tau / lbpar.agrid;
    p_pi[1] = pi[1] / lbpar.tau / lbpar.tau / lbpar.agrid;
    p_pi[2] = pi[2] / lbpar.tau / lbpar.tau / lbpar.agrid;
    p_pi[3] = pi[3] / lbpar.tau / lbpar.tau / lbpar.agrid;
    p_pi[4] = pi[4] / lbpar.tau / lbpar.tau / lbpar.agrid;
    p_pi[5] = pi[5] / lbpar.tau / lbpar.tau / lbpar.agrid;
#endif // LB
  }
  return 0;
}

int lb_lbnode_get_boundary(const Vector3i &ind, int *p_boundary) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    unsigned int host_flag;
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_get_boundary_flag_GPU(single_nodeindex, &host_flag);
    p_boundary[0] = host_flag;
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid_boundary_flag(node, index, p_boundary);
#endif // LB
  }
  return 0;
}

#endif // defined (LB) || defined (LB_GPU)

int lb_lbnode_get_pop(const Vector3i &ind, double *p_pop) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float population[19];

    // c is the LB_COMPONENT for SHANCHEN (not yet interfaced)
    int c = 0;
    lb_lbfluid_get_population(ind, population, c);

    for (int i = 0; i < LBQ; ++i)
      p_pop[i] = population[i];
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_recv_fluid_populations(node, index, p_pop);
#endif // LB
  }
  return 0;
}

int lb_lbnode_set_rho(const Vector3i &ind, double *p_rho) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float host_rho[LB_COMPONENTS];
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    int i;
    for (i = 0; i < LB_COMPONENTS; i++) {
      host_rho[i] = (float)p_rho[i];
    }
    lb_set_node_rho_GPU(single_nodeindex, host_rho);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];
    double rho;
    std::array<double, 3> j;
    std::array<double, 6> pi;

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());
    rho = (*p_rho) * lbpar.agrid * lbpar.agrid * lbpar.agrid;
    mpi_send_fluid(node, index, rho, j, pi);
#endif // LB
  }
  return 0;
}

int lb_lbnode_set_u(const Vector3i &ind, double *u) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float host_velocity[3];
    host_velocity[0] = (float)u[0] * lbpar_gpu.tau / lbpar_gpu.agrid;
    host_velocity[1] = (float)u[1] * lbpar_gpu.tau / lbpar_gpu.agrid;
    host_velocity[2] = (float)u[2] * lbpar_gpu.tau / lbpar_gpu.agrid;
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_set_node_velocity_GPU(single_nodeindex, host_velocity);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];
    double rho;
    std::array<double, 3> j;
    std::array<double, 6> pi;

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    /* transform to lattice units */

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());
    j[0] = rho * u[0] * lbpar.tau / lbpar.agrid;
    j[1] = rho * u[1] * lbpar.tau / lbpar.agrid;
    j[2] = rho * u[2] * lbpar.tau / lbpar.agrid;
    mpi_send_fluid(node, index, rho, j, pi);
#endif // LB
  }
  return 0;
}

int lb_lbnode_set_pop(const Vector3i &ind, double *p_pop) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float population[19];

    for (int i = 0; i < LBQ; ++i)
      population[i] = p_pop[i];

    // c is the LB_COMPONENT for SHANCHEN (not yet interfaced)
    int c = 0;
    lb_lbfluid_set_population(ind, population, c);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_send_fluid_populations(node, index, p_pop);
#endif // LB
  }
  return 0;
}

#ifdef LB
/********************** The Main LB Part *************************************/
/* Halo communication for push scheme */
static void halo_push_communication(LB_Fluid &lbfluid) {
  Lattice::index_t index;
  int x, y, z, count;
  int rnode, snode;
  double *buffer = nullptr, *sbuf = nullptr, *rbuf = nullptr;
  MPI_Status status;

  int yperiod = lblattice.halo_grid[0];
  int zperiod = lblattice.halo_grid[0] * lblattice.halo_grid[1];

  /***************
   * X direction *
   ***************/
  count = 5 * lblattice.halo_grid[1] * lblattice.halo_grid[2];
  sbuf = (double *)Utils::malloc(count * sizeof(double));
  rbuf = (double *)Utils::malloc(count * sizeof(double));

  /* send to right, recv from left i = 1, 7, 9, 11, 13 */
  snode = node_neighbors[1];
  rnode = node_neighbors[0];

  buffer = sbuf;
  index = get_linear_index(lblattice.grid[0] + 1, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      buffer[0] = lbfluid[1][index];
      buffer[1] = lbfluid[7][index];
      buffer[2] = lbfluid[9][index];
      buffer[3] = lbfluid[11][index];
      buffer[4] = lbfluid[13][index];
      buffer += 5;

      index += yperiod;
    }
  }

  if (node_grid[0] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
                 MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);
  } else {
    memmove(rbuf, sbuf, count * sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(1, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      lbfluid[1][index] = buffer[0];
      lbfluid[7][index] = buffer[1];
      lbfluid[9][index] = buffer[2];
      lbfluid[11][index] = buffer[3];
      lbfluid[13][index] = buffer[4];
      buffer += 5;

      index += yperiod;
    }
  }

  /* send to left, recv from right i = 2, 8, 10, 12, 14 */
  snode = node_neighbors[0];
  rnode = node_neighbors[1];

  buffer = sbuf;
  index = get_linear_index(0, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      buffer[0] = lbfluid[2][index];
      buffer[1] = lbfluid[8][index];
      buffer[2] = lbfluid[10][index];
      buffer[3] = lbfluid[12][index];
      buffer[4] = lbfluid[14][index];
      buffer += 5;

      index += yperiod;
    }
  }

  if (node_grid[0] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
                 MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);
  } else {
    memmove(rbuf, sbuf, count * sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(lblattice.grid[0], 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (y = 0; y < lblattice.halo_grid[1]; y++) {
      lbfluid[2][index] = buffer[0];
      lbfluid[8][index] = buffer[1];
      lbfluid[10][index] = buffer[2];
      lbfluid[12][index] = buffer[3];
      lbfluid[14][index] = buffer[4];
      buffer += 5;

      index += yperiod;
    }
  }

  /***************
   * Y direction *
   ***************/
  count = 5 * lblattice.halo_grid[0] * lblattice.halo_grid[2];
  sbuf = Utils::realloc(sbuf, count * sizeof(double));
  rbuf = Utils::realloc(rbuf, count * sizeof(double));

  /* send to right, recv from left i = 3, 7, 10, 15, 17 */
  snode = node_neighbors[3];
  rnode = node_neighbors[2];

  buffer = sbuf;
  index = get_linear_index(0, lblattice.grid[1] + 1, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[3][index];
      buffer[1] = lbfluid[7][index];
      buffer[2] = lbfluid[10][index];
      buffer[3] = lbfluid[15][index];
      buffer[4] = lbfluid[17][index];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  if (node_grid[1] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
                 MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);
  } else {
    memmove(rbuf, sbuf, count * sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0, 1, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[3][index] = buffer[0];
      lbfluid[7][index] = buffer[1];
      lbfluid[10][index] = buffer[2];
      lbfluid[15][index] = buffer[3];
      lbfluid[17][index] = buffer[4];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  /* send to left, recv from right i = 4, 8, 9, 16, 18 */
  snode = node_neighbors[2];
  rnode = node_neighbors[3];

  buffer = sbuf;
  index = get_linear_index(0, 0, 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[4][index];
      buffer[1] = lbfluid[8][index];
      buffer[2] = lbfluid[9][index];
      buffer[3] = lbfluid[16][index];
      buffer[4] = lbfluid[18][index];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  if (node_grid[1] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
                 MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);
  } else {
    memmove(rbuf, sbuf, count * sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0, lblattice.grid[1], 0, lblattice.halo_grid);
  for (z = 0; z < lblattice.halo_grid[2]; z++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[4][index] = buffer[0];
      lbfluid[8][index] = buffer[1];
      lbfluid[9][index] = buffer[2];
      lbfluid[16][index] = buffer[3];
      lbfluid[18][index] = buffer[4];
      buffer += 5;

      ++index;
    }
    index += zperiod - lblattice.halo_grid[0];
  }

  /***************
   * Z direction *
   ***************/
  count = 5 * lblattice.halo_grid[0] * lblattice.halo_grid[1];
  sbuf = Utils::realloc(sbuf, count * sizeof(double));
  rbuf = Utils::realloc(rbuf, count * sizeof(double));

  /* send to right, recv from left i = 5, 11, 14, 15, 18 */
  snode = node_neighbors[5];
  rnode = node_neighbors[4];

  buffer = sbuf;
  index = get_linear_index(0, 0, lblattice.grid[2] + 1, lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[5][index];
      buffer[1] = lbfluid[11][index];
      buffer[2] = lbfluid[14][index];
      buffer[3] = lbfluid[15][index];
      buffer[4] = lbfluid[18][index];
      buffer += 5;

      ++index;
    }
  }

  if (node_grid[2] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
                 MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);
  } else {
    memmove(rbuf, sbuf, count * sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0, 0, 1, lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[5][index] = buffer[0];
      lbfluid[11][index] = buffer[1];
      lbfluid[14][index] = buffer[2];
      lbfluid[15][index] = buffer[3];
      lbfluid[18][index] = buffer[4];
      buffer += 5;

      ++index;
    }
  }

  /* send to left, recv from right i = 6, 12, 13, 16, 17 */
  snode = node_neighbors[4];
  rnode = node_neighbors[5];

  buffer = sbuf;
  index = get_linear_index(0, 0, 0, lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      buffer[0] = lbfluid[6][index];
      buffer[1] = lbfluid[12][index];
      buffer[2] = lbfluid[13][index];
      buffer[3] = lbfluid[16][index];
      buffer[4] = lbfluid[17][index];
      buffer += 5;

      ++index;
    }
  }

  if (node_grid[2] > 1) {
    MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, REQ_HALO_SPREAD, rbuf, count,
                 MPI_DOUBLE, rnode, REQ_HALO_SPREAD, comm_cart, &status);
  } else {
    memmove(rbuf, sbuf, count * sizeof(double));
  }

  buffer = rbuf;
  index = get_linear_index(0, 0, lblattice.grid[2], lblattice.halo_grid);
  for (y = 0; y < lblattice.halo_grid[1]; y++) {
    for (x = 0; x < lblattice.halo_grid[0]; x++) {
      lbfluid[6][index] = buffer[0];
      lbfluid[12][index] = buffer[1];
      lbfluid[13][index] = buffer[2];
      lbfluid[16][index] = buffer[3];
      lbfluid[17][index] = buffer[4];
      buffer += 5;

      ++index;
    }
  }

  free(rbuf);
  free(sbuf);
}

/***********************************************************************/

/** Performs basic sanity checks. */
void lb_sanity_checks() {

  if (lbpar.agrid <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann agrid not set";
  }
  if (lbpar.tau <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann time step not set";
  }
  if (lbpar.rho <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann fluid density not set";
  }
  if (lbpar.viscosity <= 0.0) {
    runtimeErrorMsg() << "Lattice Boltzmann fluid viscosity not set";
  }
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    runtimeErrorMsg() << "LB requires domain-decomposition cellsystem";
  }
  if (skin == 0.0) {
    runtimeErrorMsg() << "LB requires a positive skin";
  }
  if (cell_structure.use_verlet_list && skin >= lbpar.agrid / 2.0) {
    runtimeErrorMsg() << "LB requires either no Verlet lists or that the skin "
                         "of the verlet list to be less than half of "
                         "lattice-Boltzmann grid spacing";
  }
}

void lb_on_integration_start() {
  lb_sanity_checks();

  halo_communication(&update_halo_comm,
                     reinterpret_cast<char *>(lbfluid[0].data()));
}

uint64_t lb_coupling_rng_state_cpu() { return rng_counter_coupling.value(); }

void lb_coupling_set_rng_state_cpu(uint64_t counter) {
  uint32_t high, low;
  std::tie(high, low) = Utils::u64_to_u32(counter);
  mpi_call(mpi_set_lb_coupling_counter, high, low);

  rng_counter_coupling = Utils::Counter<uint64_t>(counter);
}

uint64_t lb_fluid_rng_state_cpu() { return rng_counter_fluid.value(); }

void lb_fluid_set_rng_state_cpu(uint64_t counter) {
  uint32_t high, low;
  std::tie(high, low) = Utils::u64_to_u32(counter);
  mpi_call(mpi_set_lb_fluid_counter, high, low);

  rng_counter_fluid = Utils::Counter<uint64_t>(counter);
}
#endif

uint64_t lb_coupling_rng_state() {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lb_coupling_rng_state_cpu();
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lb_coupling_rng_state_gpu();
#endif
  }
  return {};
}

void lb_coupling_set_rng_state(uint64_t counter) {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_coupling_set_rng_state_cpu(counter);
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_coupling_set_rng_state_gpu(counter);
#endif
  }
}

uint64_t lb_fluid_rng_state() {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lb_fluid_rng_state_cpu();
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lb_fluid_rng_state_gpu();
#endif
  }
  return {};
}

void lb_fluid_set_rng_state(uint64_t counter) {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_fluid_set_rng_state_cpu(counter);
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_fluid_set_rng_state_gpu(counter);
#endif
  }
}
void mpi_set_lb_coupling_counter(int high, int low) {
#ifdef LB
  rng_counter_coupling = Utils::Counter<uint64_t>(Utils::u32_to_u64(
      static_cast<uint32_t>(high), static_cast<uint32_t>(low)));
#endif
}

void mpi_set_lb_fluid_counter(int high, int low) {
#ifdef LB
  rng_counter_fluid = Utils::Counter<uint64_t>(Utils::u32_to_u64(
      static_cast<uint32_t>(high), static_cast<uint32_t>(low)));
#endif
}

#ifdef LB
/***********************************************************************/

/** (Re-)allocate memory for the fluid and initialize pointers. */
static void lb_realloc_fluid() {
  LB_TRACE(printf("reallocating fluid\n"));
  const std::array<int, 2> size = {
      {lbmodel.n_veloc, lblattice.halo_grid_volume}};

  lbfluid_a.resize(size);
  lbfluid_b.resize(size);

  using Utils::Span;
  for (int i = 0; i < size[0]; i++) {
    lbfluid[i] = Span<double>(lbfluid_a[i].origin(), size[1]);
    lbfluid_post[i] = Span<double>(lbfluid_b[i].origin(), size[1]);
  }

  lbfields.resize(lblattice.halo_grid_volume);
}

/** Sets up the structures for exchange of the halo regions.
 *  See also \ref halo.cpp */
static void lb_prepare_communication() {
  int i;
  HaloCommunicator comm = {0, nullptr};

  /* since the data layout is a structure of arrays, we have to
   * generate a communication for this structure: first we generate
   * the communication for one of the arrays (the 0-th velocity
   * population), then we replicate this communication for the other
   * velocity indices by constructing appropriate vector
   * datatypes */

  /* prepare the communication for a single velocity */
  prepare_halo_communication(&comm, &lblattice, FIELDTYPE_DOUBLE, MPI_DOUBLE);

  update_halo_comm.num = comm.num;
  update_halo_comm.halo_info =
      Utils::realloc(update_halo_comm.halo_info, comm.num * sizeof(HaloInfo));

  /* replicate the halo structure */
  for (i = 0; i < comm.num; i++) {
    HaloInfo *hinfo = &(update_halo_comm.halo_info[i]);

    hinfo->source_node = comm.halo_info[i].source_node;
    hinfo->dest_node = comm.halo_info[i].dest_node;
    hinfo->s_offset = comm.halo_info[i].s_offset;
    hinfo->r_offset = comm.halo_info[i].r_offset;
    hinfo->type = comm.halo_info[i].type;

    /* generate the vector datatype for the structure of lattices we
     * have to use hvector here because the extent of the subtypes
     * does not span the full lattice and hence we cannot get the
     * correct vskip out of them */

    MPI_Aint lower;
    MPI_Aint extent;
    MPI_Type_get_extent(MPI_DOUBLE, &lower, &extent);
    MPI_Type_create_hvector(lbmodel.n_veloc, 1,
                            lblattice.halo_grid_volume * extent,
                            comm.halo_info[i].datatype, &hinfo->datatype);
    MPI_Type_commit(&hinfo->datatype);

    halo_create_field_hvector(lbmodel.n_veloc, 1,
                              lblattice.halo_grid_volume * sizeof(double),
                              comm.halo_info[i].fieldtype, &hinfo->fieldtype);
  }

  release_halo_communication(&comm);
}

/** (Re-)initializes the fluid. */
void lb_reinit_parameters() {
  int i;

  if (lbpar.viscosity > 0.0) {
    /* Eq. (80) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar.gamma_shear = 1. - 2. / (6. * lbpar.viscosity + 1.);
  }

  if (lbpar.bulk_viscosity > 0.0) {
    /* Eq. (81) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007). */
    lbpar.gamma_bulk = 1. - 2. / (9. * lbpar.bulk_viscosity + 1.);
  }

  if (lbpar.is_TRT) {
    lbpar.gamma_bulk = lbpar.gamma_shear;
    lbpar.gamma_even = lbpar.gamma_shear;
    lbpar.gamma_odd =
        -(7.0 * lbpar.gamma_even + 1.0) / (lbpar.gamma_even + 7.0);
    // gamma_odd = lbpar.gamma_shear; //uncomment for BGK
  }

  // lbpar.gamma_shear = 0.0; //uncomment for special case of BGK
  // lbpar.gamma_bulk = 0.0;
  // gamma_odd = 0.0;
  // gamma_even = 0.0;

  if (temperature > 0.0) {
    /* fluctuating hydrodynamics ? */
    lbpar.fluct = 1;

    /* Eq. (51) Duenweg, Schiller, Ladd, PRE 76(3):036704 (2007).
     * Note that the modes are not normalized as in the paper here! */
    double mu = temperature / lbmodel.c_sound_sq * lbpar.tau * lbpar.tau /
                (lbpar.agrid * lbpar.agrid);
    // mu *= agrid*agrid*agrid;  // Marcello's conjecture

    for (i = 0; i < 4; i++)
      lbpar.phi[i] = 0.0;
    lbpar.phi[4] =
        sqrt(mu * lbmodel.w_k[4] *
             (1. - Utils::sqr(lbpar.gamma_bulk))); // Utils::sqr(x) == x*x
    for (i = 5; i < 10; i++)
      lbpar.phi[i] =
          sqrt(mu * lbmodel.w_k[i] * (1. - Utils::sqr(lbpar.gamma_shear)));
    for (i = 10; i < 16; i++)
      lbpar.phi[i] =
          sqrt(mu * lbmodel.w_k[i] * (1 - Utils::sqr(lbpar.gamma_odd)));
    for (i = 16; i < 19; i++)
      lbpar.phi[i] =
          sqrt(mu * lbmodel.w_k[i] * (1 - Utils::sqr(lbpar.gamma_even)));

    LB_TRACE(fprintf(
        stderr,
        "%d: lbpar.gamma_shear=%lf lbpar.gamma_bulk=%lf shear_fluct=%lf "
        "bulk_fluct=%lf mu=%lf, bulkvisc=%lf, visc=%lf\n",
        this_node, lbpar.gamma_shear, lbpar.gamma_bulk, lbpar.phi[9],
        lbpar.phi[4], mu, lbpar.bulk_viscosity, lbpar.viscosity));
  } else {
    /* no fluctuations at zero temperature */
    lbpar.fluct = 0;
    for (i = 0; i < lbmodel.n_veloc; i++)
      lbpar.phi[i] = 0.0;
  }
}

/** (Re-)initializes the fluid according to the given value of rho. */
void lb_reinit_fluid() {
  std::fill(lbfields.begin(), lbfields.end(), LB_FluidNode());
  /* default values for fields in lattice units */
  /* here the conversion to lb units is performed */
  double rho = lbpar.rho;
  std::array<double, 3> j = {{0., 0., 0.}};
  std::array<double, 6> pi = {{0., 0., 0., 0., 0., 0.}};

  LB_TRACE(fprintf(stderr,
                   "Initialising the fluid with equilibrium populations\n"););

  for (Lattice::index_t index = 0; index < lblattice.halo_grid_volume;
       ++index) {
    // calculate equilibrium distribution
    lb_calc_n_from_rho_j_pi(index, rho, j, pi);

#ifdef LB_BOUNDARIES
    lbfields[index].boundary = 0;
#endif // LB_BOUNDARIES
  }

#ifdef LB_BOUNDARIES
  LBBoundaries::lb_init_boundaries();
#endif // LB_BOUNDARIES
}

/** Performs a full initialization of
 *  the Lattice Boltzmann system. All derived parameters
 *  and the fluid are reset to their default values. */
void lb_init() {
  LB_TRACE(printf("Begin initialzing fluid on CPU\n"));

  if (lbpar.agrid <= 0.0) {
    runtimeErrorMsg()
        << "Lattice Boltzmann agrid not set when initializing fluid";
  }

  if (check_runtime_errors())
    return;

  double temp_agrid[3];
  double temp_offset[3];
  for (int i = 0; i < 3; i++) {
    temp_agrid[i] = lbpar.agrid;
    temp_offset[i] = 0.5;
  }

  /* initialize the local lattice domain */
  lblattice.init(temp_agrid, temp_offset, 1, 0);

  if (check_runtime_errors())
    return;

  /* allocate memory for data structures */
  lb_realloc_fluid();

  /* prepare the halo communication */
  lb_prepare_communication();

  /* initialize derived parameters */
  lb_reinit_parameters();

  /* setup the initial particle velocity distribution */
  lb_reinit_fluid();

  LB_TRACE(printf("Initialzing fluid on CPU successful\n"));
}

/** Release fluid and communication. */
void lb_release() { release_halo_communication(&update_halo_comm); }

/***********************************************************************/
/** \name Mapping between hydrodynamic fields and particle populations */
/***********************************************************************/
/*@{*/
void lb_calc_n_from_rho_j_pi(const Lattice::index_t index, const double rho,
                             const std::array<double, 3> &j,
                             const std::array<double, 6> &pi) {
  int i;
  double local_rho, local_j[3], local_pi[6], trace;
  local_rho = rho;

  local_j[0] = j[0];
  local_j[1] = j[1];
  local_j[2] = j[2];

  for (i = 0; i < 6; i++)
    local_pi[i] = pi[i];

  trace = local_pi[0] + local_pi[2] + local_pi[5];

  double rho_times_coeff;
  double tmp1, tmp2;

  /* update the q=0 sublattice */
  lbfluid[0][index] = 1. / 3. * (local_rho - lbpar.rho) - 1. / 2. * trace;

  /* update the q=1 sublattice */
  rho_times_coeff = 1. / 18. * (local_rho - lbpar.rho);

  lbfluid[1][index] = rho_times_coeff + 1. / 6. * local_j[0] +
                      1. / 4. * local_pi[0] - 1. / 12. * trace;
  lbfluid[2][index] = rho_times_coeff - 1. / 6. * local_j[0] +
                      1. / 4. * local_pi[0] - 1. / 12. * trace;
  lbfluid[3][index] = rho_times_coeff + 1. / 6. * local_j[1] +
                      1. / 4. * local_pi[2] - 1. / 12. * trace;
  lbfluid[4][index] = rho_times_coeff - 1. / 6. * local_j[1] +
                      1. / 4. * local_pi[2] - 1. / 12. * trace;
  lbfluid[5][index] = rho_times_coeff + 1. / 6. * local_j[2] +
                      1. / 4. * local_pi[5] - 1. / 12. * trace;
  lbfluid[6][index] = rho_times_coeff - 1. / 6. * local_j[2] +
                      1. / 4. * local_pi[5] - 1. / 12. * trace;

  /* update the q=2 sublattice */
  rho_times_coeff = 1. / 36. * (local_rho - lbpar.rho);

  tmp1 = local_pi[0] + local_pi[2];
  tmp2 = 2.0 * local_pi[1];

  lbfluid[7][index] = rho_times_coeff + 1. / 12. * (local_j[0] + local_j[1]) +
                      1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  lbfluid[8][index] = rho_times_coeff - 1. / 12. * (local_j[0] + local_j[1]) +
                      1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  lbfluid[9][index] = rho_times_coeff + 1. / 12. * (local_j[0] - local_j[1]) +
                      1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  lbfluid[10][index] = rho_times_coeff - 1. / 12. * (local_j[0] - local_j[1]) +
                       1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;

  tmp1 = local_pi[0] + local_pi[5];
  tmp2 = 2.0 * local_pi[3];

  lbfluid[11][index] = rho_times_coeff + 1. / 12. * (local_j[0] + local_j[2]) +
                       1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  lbfluid[12][index] = rho_times_coeff - 1. / 12. * (local_j[0] + local_j[2]) +
                       1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  lbfluid[13][index] = rho_times_coeff + 1. / 12. * (local_j[0] - local_j[2]) +
                       1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  lbfluid[14][index] = rho_times_coeff - 1. / 12. * (local_j[0] - local_j[2]) +
                       1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;

  tmp1 = local_pi[2] + local_pi[5];
  tmp2 = 2.0 * local_pi[4];

  lbfluid[15][index] = rho_times_coeff + 1. / 12. * (local_j[1] + local_j[2]) +
                       1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  lbfluid[16][index] = rho_times_coeff - 1. / 12. * (local_j[1] + local_j[2]) +
                       1. / 8. * (tmp1 + tmp2) - 1. / 24. * trace;
  lbfluid[17][index] = rho_times_coeff + 1. / 12. * (local_j[1] - local_j[2]) +
                       1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
  lbfluid[18][index] = rho_times_coeff - 1. / 12. * (local_j[1] - local_j[2]) +
                       1. / 8. * (tmp1 - tmp2) - 1. / 24. * trace;
}

/*@}*/

#include <boost/range/numeric.hpp>

template <typename T, std::size_t N>
std::array<T, N> lb_calc_m_from_n(const std::array<T, N> n) {
  return Utils::matrix_vector_product<T, N, ::D3Q19::e_ki>(n);
}

/** Calculation of hydrodynamic modes */
std::array<double, 19> lb_calc_modes(Lattice::index_t index) {
  std::array<double, 19> n;
  for (int i = 0; i < 19; i++) {
    n[i] = lbfluid[i][index];
  }
  return lb_calc_m_from_n(n);
}

template <typename T>
inline std::array<T, 19> lb_relax_modes(Lattice::index_t index,
                                        const std::array<T, 19> &modes) {
  T rho, j[3], pi_eq[6];

  /* re-construct the real density
   * remember that the populations are stored as differences to their
   * equilibrium value */
  rho = modes[0] + lbpar.rho;

  j[0] = modes[1] + 0.5 * lbfields[index].force_density[0];
  j[1] = modes[2] + 0.5 * lbfields[index].force_density[1];
  j[2] = modes[3] + 0.5 * lbfields[index].force_density[2];

  /* equilibrium part of the stress modes */
  pi_eq[0] = scalar(j, j) / rho;
  pi_eq[1] = (Utils::sqr(j[0]) - Utils::sqr(j[1])) / rho;
  pi_eq[2] = (scalar(j, j) - 3.0 * Utils::sqr(j[2])) / rho;
  pi_eq[3] = j[0] * j[1] / rho;
  pi_eq[4] = j[0] * j[2] / rho;
  pi_eq[5] = j[1] * j[2] / rho;

  return {{modes[0], modes[1], modes[2], modes[3],
           /* relax the stress modes */
           pi_eq[0] + lbpar.gamma_bulk * (modes[4] - pi_eq[0]),
           pi_eq[1] + lbpar.gamma_shear * (modes[5] - pi_eq[1]),
           pi_eq[2] + lbpar.gamma_shear * (modes[6] - pi_eq[2]),
           pi_eq[3] + lbpar.gamma_shear * (modes[7] - pi_eq[3]),
           pi_eq[4] + lbpar.gamma_shear * (modes[8] - pi_eq[4]),
           pi_eq[5] + lbpar.gamma_shear * (modes[9] - pi_eq[5]),
           /* relax the ghost modes (project them out) */
           /* ghost modes have no equilibrium part due to orthogonality */
           lbpar.gamma_odd * modes[10], lbpar.gamma_odd * modes[11],
           lbpar.gamma_odd * modes[12], lbpar.gamma_odd * modes[13],
           lbpar.gamma_odd * modes[14], lbpar.gamma_odd * modes[15],
           lbpar.gamma_even * modes[16], lbpar.gamma_even * modes[17],
           lbpar.gamma_even * modes[18]}};
}

template <typename T>
inline std::array<T, 19>
lb_thermalize_modes(Lattice::index_t index, const r123::Philox4x64::ctr_type &c,
                    const std::array<T, 19> &modes) {
  if (lbpar.fluct) {
    using Utils::uniform;
    using rng_type = r123::Philox4x64;
    using ctr_type = rng_type::ctr_type;

    const T rootrho = std::sqrt(std::fabs(modes[0] + lbpar.rho));
    auto const pref = std::sqrt(12.) * rootrho;

    const ctr_type noise[4] = {
        rng_type{}(c, {{static_cast<uint64_t>(index), 0ul}}),
        rng_type{}(c, {{static_cast<uint64_t>(index), 1ul}}),
        rng_type{}(c, {{static_cast<uint64_t>(index), 2ul}}),
        rng_type{}(c, {{static_cast<uint64_t>(index), 3ul}})};

    auto rng = [&](int i) { return uniform(noise[i / 4][i % 4]); };

    return {/* conserved modes */
            {modes[0], modes[1], modes[2], modes[3],
             /* stress modes */
             modes[4] + pref * lbpar.phi[4] * rng(0),
             modes[5] + pref * lbpar.phi[5] * rng(1),
             modes[6] + pref * lbpar.phi[6] * rng(2),
             modes[7] + pref * lbpar.phi[7] * rng(3),
             modes[8] + pref * lbpar.phi[8] * rng(4),
             modes[9] + pref * lbpar.phi[9] * rng(5),

             /* ghost modes */
             modes[10] + pref * lbpar.phi[10] * rng(6),
             modes[11] + pref * lbpar.phi[11] * rng(7),
             modes[12] + pref * lbpar.phi[12] * rng(8),
             modes[13] + pref * lbpar.phi[13] * rng(9),
             modes[14] + pref * lbpar.phi[14] * rng(10),
             modes[15] + pref * lbpar.phi[15] * rng(11),
             modes[16] + pref * lbpar.phi[16] * rng(12),
             modes[17] + pref * lbpar.phi[17] * rng(13),
             modes[18] + pref * lbpar.phi[18] * rng(14)}};
  } else {
    return modes;
  }
}

template <typename T>
std::array<T, 19> lb_apply_forces(Lattice::index_t index,
                                  const std::array<T, 19> &modes) {
  T rho, u[3], C[6];

  const auto &f = lbfields[index].force_density;

  rho = modes[0] + lbpar.rho;

  /* hydrodynamic momentum density is redefined when external forces present */
  u[0] = (modes[1] + 0.5 * f[0]) / rho;
  u[1] = (modes[2] + 0.5 * f[1]) / rho;
  u[2] = (modes[3] + 0.5 * f[2]) / rho;

  C[0] = (1. + lbpar.gamma_bulk) * u[0] * f[0] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * scalar(u, f);
  C[2] = (1. + lbpar.gamma_bulk) * u[1] * f[1] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * scalar(u, f);
  C[5] = (1. + lbpar.gamma_bulk) * u[2] * f[2] +
         1. / 3. * (lbpar.gamma_bulk - lbpar.gamma_shear) * scalar(u, f);
  C[1] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[0] * f[1] + u[1] * f[0]);
  C[3] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[0] * f[2] + u[2] * f[0]);
  C[4] = 1. / 2. * (1. + lbpar.gamma_shear) * (u[1] * f[2] + u[2] * f[1]);

  return {{modes[0],
           /* update momentum modes */
           modes[1] + f[0], modes[2] + f[1], modes[3] + f[2],
           /* update stress modes */
           modes[4] + C[0] + C[2] + C[5], modes[5] + C[0] - C[2],
           modes[6] + C[0] + C[2] - 2. * C[5], modes[7] + C[1], modes[8] + C[3],
           modes[9] + C[4], modes[10], modes[11], modes[12], modes[13],
           modes[14], modes[15], modes[16], modes[17], modes[18]}};
}

template <typename T>
std::array<T, 19> normalize_modes(const std::array<T, 19> &modes) {
  auto normalized_modes = modes;
  for (int i = 0; i < modes.size(); i++) {
    normalized_modes[i] /= lbmodel.w_k[i];
  }
  return normalized_modes;
}

template <typename T, std::size_t N>
std::array<T, N> lb_calc_n_from_m(const std::array<T, N> &modes) {
  auto ret = Utils::matrix_vector_product<T, N, ::D3Q19::e_ki_transposed>(
      normalize_modes(modes));
  std::transform(ret.begin(), ret.end(), ::D3Q19::w.begin(), ret.begin(),
                 std::multiplies<T>());
  return ret;
}

template <typename T>
inline void lb_calc_n_from_modes_push(LB_Fluid &lbfluid, Lattice::index_t index,
                                      const std::array<T, 19> m) {
  const std::array<int, 3> period = {
      {1, lblattice.halo_grid[0],
       lblattice.halo_grid[0] * lblattice.halo_grid[1]}};
  auto const f = lb_calc_n_from_m(m);
  for (int i = 0; i < 19; i++) {
    auto const next = index + boost::inner_product(period, lbmodel.c[i], 0);
    lbfluid[i][next] = f[i];
  }
}

/* Collisions and streaming (push scheme) */
inline void lb_collide_stream() {
/* loop over all lattice cells (halo excluded) */
#ifdef LB_BOUNDARIES
  for (auto it = LBBoundaries::lbboundaries.begin();
       it != LBBoundaries::lbboundaries.end(); ++it) {
    (**it).reset_force();
  }
#endif // LB_BOUNDARIES

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  // Safeguard the node forces so that we can later use them for the IBM
  // particle update
  // In the following loop the lbfields[XX].force are reset to zero
  // Safeguard the node forces so that we can later use them for the IBM
  // particle update In the following loop the lbfields[XX].force are reset to
  // zero
  for (int i = 0; i < lblattice.halo_grid_volume; ++i) {
    lbfields[i].force_density_buf = lbfields[i].force_density;
  }
#endif

  const r123::Philox4x64::ctr_type c{
      {rng_counter_fluid.value(), static_cast<uint64_t>(RNGSalt::FLUID)}};

  Lattice::index_t index = lblattice.halo_offset;
  for (int z = 1; z <= lblattice.grid[2]; z++) {
    for (int y = 1; y <= lblattice.grid[1]; y++) {
      for (int x = 1; x <= lblattice.grid[0]; x++) {
// as we only want to apply this to non-boundary nodes we can throw out
// the if-clause if we have a non-bounded domain
#ifdef LB_BOUNDARIES
        if (!lbfields[index].boundary)
#endif // LB_BOUNDARIES
        {
          /* calculate modes locally */
          auto const modes = lb_calc_modes(index);

          /* deterministic collisions */
          auto const relaxed_modes = lb_relax_modes(index, modes);

          /* fluctuating hydrodynamics */
          auto const thermalized_modes =
              lb_thermalize_modes(index, c, relaxed_modes);

          /* apply forces */
          auto const modes_with_forces =
              lb_apply_forces(index, thermalized_modes);

          /* reset the force density */
          lbfields[index].force_density = lbpar.ext_force_density;

          /* transform back to populations and streaming */
          lb_calc_n_from_modes_push(lbfluid_post, index, modes_with_forces);
        }

        ++index; /* next node */
      }
      index += 2; /* skip halo region */
    }
    index += 2 * lblattice.halo_grid[0]; /* skip halo region */
  }

  /* exchange halo regions */
  halo_push_communication(lbfluid_post);

#ifdef LB_BOUNDARIES
  /* boundary conditions for links */
  LBBoundaries::lb_bounce_back(lbfluid_post);
#endif // LB_BOUNDARIES

  rng_counter_fluid.increment();

  /* swap the pointers for old and new population fields */
  std::swap(lbfluid, lbfluid_post);

  halo_communication(&update_halo_comm,
                     reinterpret_cast<char *>(lbfluid[0].data()));

#ifdef ADDITIONAL_CHECKS
  lb_check_halo_regions(lbfluid);
#endif
}

/***********************************************************************/
/** \name Update step for the lattice Boltzmann fluid                  */
/***********************************************************************/
/*@{*/
/*@}*/

/** Update the lattice Boltzmann fluid.
 *
 * This function is called from the integrator. Since the time step
 * for the lattice dynamics can be coarser than the MD time step, we
 * monitor the time since the last lattice update.
 */
void lattice_boltzmann_update() {
  int factor = (int)round(lbpar.tau / time_step);

  fluidstep += 1;
  if (fluidstep >= factor) {
    fluidstep = 0;

    lb_collide_stream();
  }
}

namespace {
template <typename Op>
void lattice_interpolation(Lattice const &lattice, Vector3d const &pos,
                           Op &&op) {
  Lattice::index_t node_index[8];
  double delta[6];

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
} // namespace

/***********************************************************************/
/** \name Coupling part */
/***********************************************************************/
/*@{*/

namespace {
/**
 * @brief Add a force to the lattice force density.
 * @param pos Position of the force
 * @param force Force in MD units.
 */
void add_md_force(Vector3d const &pos, Vector3d const &force) {
  /* transform momentum transfer to lattice units
     (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  auto const delta_j = -(time_step * lbpar.tau / lbpar.agrid) * force;

  lattice_interpolation(lblattice, pos,
                        [&delta_j](Lattice::index_t index, double w) {
                          auto &node = lbfields[index];

                          node.force_density[0] += w * delta_j[0];
                          node.force_density[1] += w * delta_j[1];
                          node.force_density[2] += w * delta_j[2];
                        });
}
} // namespace

namespace {
bool in_local_domain(Vector3d const &pos) {
  return (pos[0] >= my_left[0] - 0.5 * lblattice.agrid[0] &&
          pos[0] < my_right[0] + 0.5 * lblattice.agrid[0] &&
          pos[1] >= my_left[1] - 0.5 * lblattice.agrid[1] &&
          pos[1] < my_right[1] + 0.5 * lblattice.agrid[1] &&
          pos[2] >= my_left[2] - 0.5 * lblattice.agrid[2] &&
          pos[2] < my_right[2] + 0.5 * lblattice.agrid[2]);
}

#ifdef ENGINE
void add_swimmer_force(Particle &p) {
  if (p.swim.swimming) {
    // calculate source position
    const double direction = double(p.swim.push_pull) * p.swim.dipole_length;
    auto const director = p.r.calc_director();
    auto const source_position = p.r.p + direction * director;

    if (not in_local_domain(source_position)) {
      return;
    }

    p.swim.v_source = lb_lbfluid_get_interpolated_velocity(source_position);

    add_md_force(source_position, p.swim.f_swim * director);
  }
}
#endif
} // namespace

/** Coupling of a single particle to viscous fluid with Stokesian friction.
 *
 * Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
 *
 * @param p          The coupled particle (Input).
 * @param f_random   Additional force to be included.
 *
 * @return The viscous coupling force plus f_random.
 */
inline Vector3d lb_viscous_coupling(Particle *p, Vector3d const &f_random) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  auto const interpolated_u = lb_lbfluid_get_interpolated_velocity(p->r.p);

  Vector3d v_drift = interpolated_u;
#ifdef ENGINE
  if (p->swim.swimming) {
    v_drift += p->swim.v_swim * p->r.calc_director();
    p->swim.v_center[0] = interpolated_u[0];
    p->swim.v_center[1] = interpolated_u[1];
    p->swim.v_center[2] = interpolated_u[2];
  }
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  v_drift += p->p.mu_E;
#endif

  /* calculate viscous force
   * (Eq. (9) Ahlrichs and Duenweg, JCP 111(17):8225 (1999))
   * */
  auto const force = -lbpar.friction * (p->m.v - v_drift) + f_random;

  add_md_force(p->r.p, force);

  return force;
}

namespace {
Vector3d node_u(Lattice::index_t index) {
#ifdef LB_BOUNDARIES
  if (lbfields[index].boundary) {
    return lbfields[index].slip_velocity;
  }
#endif // LB_BOUNDARIES
  auto const modes = lb_calc_modes(index);
  auto const local_rho = lbpar.rho + modes[0];
  return Vector3d{modes[1], modes[2], modes[3]} / local_rho;
}
} // namespace

/*
 * @brief Interpolate the fluid velocity.
 *
 * @param pos Position
 * @param v Interpolated velocity in MD units.
 */
Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d &pos) {
  Vector3d interpolated_u{};

  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  lattice_interpolation(lblattice, pos,
                        [&interpolated_u](Lattice::index_t index, double w) {
                          interpolated_u += w * node_u(index);
                        });

  return (lbpar.agrid / lbpar.tau) * interpolated_u;
}

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
Vector3d lb_lbfluid_get_interpolated_force(const Vector3d &pos) {
  Vector3d interpolated_f{};
  lattice_interpolation(
      lblattice, pos, [&interpolated_f](Lattice::index_t index, double w) {
#ifdef LB_BOUNDARIES
        if (!lbfields[index].boundary) {
#endif
          auto const modes = lb_calc_modes(index);
          auto const local_rho = modes[0] + lbpar.rho;
          interpolated_f +=
              w / 2 / local_rho *
              (lbfields[index].force_density_buf - lbpar.ext_force_density);
#ifdef LB_BOUNDARIES
        }
#endif
      });
  return interpolated_f;
}
#endif

/** Calculate particle lattice interactions.
 * So far, only viscous coupling with Stokesian friction is
 * implemented.
 * Include all particle-lattice forces in this function.
 * The function is called from \ref force_calc.
 *
 * Parallelizing the fluid particle coupling is not straightforward
 * because drawing of random numbers makes the whole thing nonlocal.
 * One way to do it is to treat every particle only on one node, i.e.
 * the random numbers need not be communicated. The particles that are
 * not fully inside the local lattice are taken into account via their
 * ghost images on the neighbouring nodes. But this requires that the
 * correct values of the surrounding lattice nodes are available on
 * the respective node, which means that we have to communicate the
 * halo regions before treating the ghost particles. Moreover, after
 * determining the ghost couplings, we have to communicate back the
 * halo region such that all local lattice nodes have the correct values.
 * Thus two communication phases are involved which will most likely be
 * the bottleneck of the computation.
 *
 * Another way of dealing with the particle lattice coupling is to
 * treat a particle and all of it's images explicitly. This requires the
 * communication of the random numbers used in the calculation of the
 * coupling force. The problem is now that, if random numbers have to
 * be redrawn, we cannot efficiently determine which particles and which
 * images have to be re-calculated. We therefore go back to the outset
 * and go through the whole system again until no failure occurs during
 * such a sweep. In the worst case, this is very inefficient because
 * many things are recalculated although they actually don't need.
 * But we can assume that this happens extremely rarely and then we have
 * on average only one communication phase for the random numbers, which
 * probably makes this method preferable compared to the above one.
 */
void calc_particle_lattice_ia() {

  if (transfer_momentum) {
    using rng_type = r123::Philox4x64;
    using ctr_type = rng_type::ctr_type;
    using key_type = rng_type::key_type;

    ctr_type c{{rng_counter_coupling.value(),
                static_cast<uint64_t>(RNGSalt::PARTICLES)}};
    rng_counter_coupling.increment();

    /* Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * time_step comes from the discretization.
     */
    auto const noise_amplitude =
        sqrt(12. * 2. * lbpar.friction * temperature / time_step);

    auto f_random = [&c](int id) -> Vector3d {
      key_type k{{static_cast<uint32_t>(id)}};

      auto const noise = rng_type{}(c, k);

      using Utils::uniform;
      return Vector3d{uniform(noise[0]), uniform(noise[1]), uniform(noise[2])} -
             Vector3d::broadcast(0.5);
    };

    /* local cells */
    for (auto &p : local_cells.particles()) {
      if (!p.p.is_virtual || thermo_virtual) {
        auto const force =
            lb_viscous_coupling(&p, noise_amplitude * f_random(p.identity()));
        /* add force to the particle */
        p.f.f += force;
#ifdef ENGINE
        add_swimmer_force(p);
#endif
      }
    }

    /* ghost cells */
    for (auto &p : ghost_cells.particles()) {
      /* for ghost particles we have to check if they lie
       * in the range of the local lattice nodes */
      if (in_local_domain(p.r.p)) {
        if (!p.p.is_virtual || thermo_virtual) {
          lb_viscous_coupling(&p, noise_amplitude * f_random(p.identity()));
#ifdef ENGINE
          add_swimmer_force(p);
#endif
        }
      }
    }
  }
}

/***********************************************************************/

/** Calculate the average density of the fluid in the system.
 * This function has to be called after changing the density of
 * a local lattice site in order to set lbpar.rho consistently. */
void lb_calc_average_rho() {
  Lattice::index_t index;
  double rho, local_rho, sum_rho;

  rho = 0.0;
  local_rho = 0.0;
  index = 0;
  for (int z = 1; z <= lblattice.grid[2]; z++) {
    for (int y = 1; y <= lblattice.grid[1]; y++) {
      for (int x = 1; x <= lblattice.grid[0]; x++) {
        lb_calc_local_rho(index, &rho);
        local_rho += rho;

        index++;
      }
      // skip halo region
      index += 2;
    }
    // skip halo region
    index += 2 * lblattice.halo_grid[0];
  }
  MPI_Allreduce(&rho, &sum_rho, 1, MPI_DOUBLE, MPI_SUM, comm_cart);

  /* calculate average density in MD units */
  // TODO!!!
  lbpar.rho = sum_rho / (box_l[0] * box_l[1] * box_l[2]);
}

/*@}*/

static int compare_buffers(double *buf1, double *buf2, int size) {
  int ret;
  if (memcmp(buf1, buf2, size) != 0) {
    runtimeErrorMsg() << "Halo buffers are not identical";
    ret = 1;
  } else {
    ret = 0;
  }
  return ret;
}

/** Checks consistency of the halo regions (ADDITIONAL_CHECKS)
    This function can be used as an additional check. It test whether the
    halo regions have been exchanged correctly.
*/
void lb_check_halo_regions(const LB_Fluid &lbfluid) {
  Lattice::index_t index;
  int i, x, y, z, s_node, r_node, count = lbmodel.n_veloc;
  double *s_buffer, *r_buffer;
  MPI_Status status[2];

  r_buffer = (double *)Utils::malloc(count * sizeof(double));
  s_buffer = (double *)Utils::malloc(count * sizeof(double));

  if (PERIODIC(0)) {
    for (z = 0; z < lblattice.halo_grid[2]; ++z) {
      for (y = 0; y < lblattice.halo_grid[1]; ++y) {
        index = get_linear_index(0, y, z, lblattice.halo_grid);
        for (i = 0; i < lbmodel.n_veloc; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[1];
        r_node = node_neighbors[0];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index =
              get_linear_index(lblattice.grid[0], y, z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index =
              get_linear_index(lblattice.grid[0], y, z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=" << 0 << " at index=" << index
                      << " y=" << y << " z=" << z << "\n";
          }
        }

        index =
            get_linear_index(lblattice.grid[0] + 1, y, z, lblattice.halo_grid);
        for (i = 0; i < lbmodel.n_veloc; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[0];
        r_node = node_neighbors[1];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index = get_linear_index(1, y, z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index = get_linear_index(1, y, z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=0 at index=" << index
                      << " y=" << y << " z=" << z << "\n";
          }
        }
      }
    }
  }

  if (PERIODIC(1)) {
    for (z = 0; z < lblattice.halo_grid[2]; ++z) {
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index = get_linear_index(x, 0, z, lblattice.halo_grid);
        for (i = 0; i < lbmodel.n_veloc; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[3];
        r_node = node_neighbors[2];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index =
              get_linear_index(x, lblattice.grid[1], z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index =
              get_linear_index(x, lblattice.grid[1], z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=1 at index=" << index
                      << " x=" << x << " z=" << z << "\n";
          }
        }
      }
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index =
            get_linear_index(x, lblattice.grid[1] + 1, z, lblattice.halo_grid);
        for (i = 0; i < lbmodel.n_veloc; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[2];
        r_node = node_neighbors[3];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index = get_linear_index(x, 1, z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index = get_linear_index(x, 1, z, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=1 at index=" << index
                      << " x=" << x << " z=" << z << "\n";
          }
        }
      }
    }
  }

  if (PERIODIC(2)) {
    for (y = 0; y < lblattice.halo_grid[1]; ++y) {
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index = get_linear_index(x, y, 0, lblattice.halo_grid);
        for (i = 0; i < lbmodel.n_veloc; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[5];
        r_node = node_neighbors[4];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index =
              get_linear_index(x, y, lblattice.grid[2], lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index =
              get_linear_index(x, y, lblattice.grid[2], lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=2 at index=" << index
                      << " x=" << x << " y=" << y << " z=" << lblattice.grid[2]
                      << "\n";
          }
        }
      }
    }
    for (y = 0; y < lblattice.halo_grid[1]; ++y) {
      for (x = 0; x < lblattice.halo_grid[0]; ++x) {
        index =
            get_linear_index(x, y, lblattice.grid[2] + 1, lblattice.halo_grid);
        for (i = 0; i < lbmodel.n_veloc; i++)
          s_buffer[i] = lbfluid[i][index];

        s_node = node_neighbors[4];
        r_node = node_neighbors[5];
        if (n_nodes > 1) {
          MPI_Sendrecv(s_buffer, count, MPI_DOUBLE, r_node, REQ_HALO_CHECK,
                       r_buffer, count, MPI_DOUBLE, s_node, REQ_HALO_CHECK,
                       comm_cart, status);
          index = get_linear_index(x, y, 1, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            s_buffer[i] = lbfluid[i][index];
          compare_buffers(s_buffer, r_buffer, count * sizeof(double));
        } else {
          index = get_linear_index(x, y, 1, lblattice.halo_grid);
          for (i = 0; i < lbmodel.n_veloc; i++)
            r_buffer[i] = lbfluid[i][index];
          if (compare_buffers(s_buffer, r_buffer, count * sizeof(double))) {
            std::cerr << "buffers differ in dir=2 at index=" << index
                      << " x=" << x << " y=" << y << "\n";
          }
        }
      }
    }
  }

  free(r_buffer);
  free(s_buffer);
}

void lb_calc_local_fields(Lattice::index_t index, double *rho, double *j,
                          double *pi) {

  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_calc_local_fields in " << __FILE__
                      << __LINE__ << ": CPU LB not switched on.";
    *rho = 0;
    j[0] = j[1] = j[2] = 0;
    pi[0] = pi[1] = pi[2] = pi[3] = pi[4] = pi[5] = 0;
    return;
  }

  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_calc_local_pi in " << __FILE__ << __LINE__
                      << ": CPU LB not switched on.";
    j[0] = j[1] = j[2] = 0;
    return;
  }

#ifdef LB_BOUNDARIES
  if (lbfields[index].boundary) {
    *rho = lbpar.rho;
    j[0] = 0.;
    j[1] = 0.;
    j[2] = 0.;
    if (pi) {
      pi[0] = 0.;
      pi[1] = 0.;
      pi[2] = 0.;
      pi[3] = 0.;
      pi[4] = 0.;
      pi[5] = 0.;
    }
    return;
  }
#endif
  double modes_from_pi_eq[6];
  std::array<double, 19> mode = lb_calc_modes(index);

  *rho = mode[0] + lbpar.rho;

  j[0] = mode[1] + 0.5 * lbfields[index].force_density[0];
  j[1] = mode[2] + 0.5 * lbfields[index].force_density[1];
  j[2] = mode[3] + 0.5 * lbfields[index].force_density[2];

  if (!pi)
    return;

  /* equilibrium part of the stress modes */
  modes_from_pi_eq[0] = scalar(j, j) / *rho;
  modes_from_pi_eq[1] = (Utils::sqr(j[0]) - Utils::sqr(j[1])) / *rho;
  modes_from_pi_eq[2] = (scalar(j, j) - 3.0 * Utils::sqr(j[2])) / *rho;
  modes_from_pi_eq[3] = j[0] * j[1] / *rho;
  modes_from_pi_eq[4] = j[0] * j[2] / *rho;
  modes_from_pi_eq[5] = j[1] * j[2] / *rho;

  /* Now we must predict the outcome of the next collision */
  /* We immediately average pre- and post-collision. */
  mode[4] = modes_from_pi_eq[0] +
            (0.5 + 0.5 * lbpar.gamma_bulk) * (mode[4] - modes_from_pi_eq[0]);
  mode[5] = modes_from_pi_eq[1] +
            (0.5 + 0.5 * lbpar.gamma_shear) * (mode[5] - modes_from_pi_eq[1]);
  mode[6] = modes_from_pi_eq[2] +
            (0.5 + 0.5 * lbpar.gamma_shear) * (mode[6] - modes_from_pi_eq[2]);
  mode[7] = modes_from_pi_eq[3] +
            (0.5 + 0.5 * lbpar.gamma_shear) * (mode[7] - modes_from_pi_eq[3]);
  mode[8] = modes_from_pi_eq[4] +
            (0.5 + 0.5 * lbpar.gamma_shear) * (mode[8] - modes_from_pi_eq[4]);
  mode[9] = modes_from_pi_eq[5] +
            (0.5 + 0.5 * lbpar.gamma_shear) * (mode[9] - modes_from_pi_eq[5]);

  // Transform the stress tensor components according to the modes that
  // correspond to those used by U. Schiller. In terms of populations this
  // expression then corresponds exactly to those in Eqs. 116 - 121 in the
  // Duenweg and Ladd paper, when these are written out in populations.
  // But to ensure this, the expression in Schiller's modes has to be different!

  pi[0] = (2.0 * (mode[0] + mode[4]) + mode[6] + 3.0 * mode[5]) / 6.0; // xx
  pi[1] = mode[7];                                                     // xy
  pi[2] = (2.0 * (mode[0] + mode[4]) + mode[6] - 3.0 * mode[5]) / 6.0; // yy
  pi[3] = mode[8];                                                     // xz
  pi[4] = mode[9];                                                     // yz
  pi[5] = (mode[0] + mode[4] - mode[6]) / 3.0;                         // zz
}

#endif // LB

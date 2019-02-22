#include <fstream>

#include "communication.hpp"
#include "config.hpp"
#include "electrokinetics.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "lb.hpp"
#include "lb_interface.hpp"
#include "lbgpu.hpp"
#include "thermostat.hpp"

#if defined(LB) || defined(LB_GPU)

void lb_lbfluid_update() {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lattice_boltzmann_update();
#endif
  } else if (lattice_switch & LATTICE_LB_GPU and this_node == 0) {
#ifdef LB_GPU
#ifdef ELECTROKINETICS
    if (ek_initialized) {
      ek_integrate();
    } else {
#endif
      lattice_boltzmann_update_gpu();
#ifdef ELECTROKINETICS
    }
#endif
#endif
  }
}

void lb_lbfluid_propagate() {
  lb_lbfluid_update();
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    rng_counter_fluid_gpu.increment();
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    rng_counter_fluid.increment();
#endif
  }
}

void lb_lbfluid_on_integration_start() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_GPU_sanity_checks();
    if (this_node == 0 && lb_reinit_particles_gpu()) {
      lb_realloc_particles_gpu();
      lb_reinit_particles_gpu.validate();
    }
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_sanity_checks();

    halo_communication(&update_halo_comm,
                       reinterpret_cast<char *>(lbfluid[0].data()));
#endif
  }
}

void lb_lbfluid_invalidate_particle_allocation() {
#ifdef LB_GPU
  lb_reinit_particles_gpu.invalidate();
#endif
}

/** (Re-)initialize the fluid. */
void lb_lbfluid_reinit_parameters() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    if (this_node == 0)
      lb_reinit_parameters_gpu();
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_reinit_parameters();
#endif
  }
}

/** (Re-)initialize the fluid according to the value of rho. */

void lb_lbfluid_reinit_fluid() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_reinit_fluid_gpu();
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_reinit_fluid();
#endif
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

/** Perform a full initialization of the lattice Boltzmann system.
 *  All derived parameters and the fluid are reset to their default values.
 */
void lb_lbfluid_init() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_init_gpu();
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_init();
#endif
  }
}

uint64_t lb_lbfluid_get_rng_state() {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lb_fluid_get_rng_state();
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lb_fluid_get_rng_state_gpu();
#endif
  }
  return {};
}

void lb_lbfluid_set_rng_state(uint64_t counter) {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_fluid_set_rng_state(counter);
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_fluid_set_rng_state_gpu(counter);
#endif
  }
}

void lb_lbfluid_set_density(double p_dens) {
  if (p_dens <= 0)
    throw std::invalid_argument("Density has to be > 0.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.rho = static_cast<float>(p_dens);
    lb_lbfluid_on_lb_params_change(LBPAR_DENSITY);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.rho = p_dens;
    mpi_bcast_lb_params(LBPAR_DENSITY);
#endif // LB
  }
}

double lb_lbfluid_get_density() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return static_cast<double>(lbpar_gpu.rho);
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.rho;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbfluid_set_viscosity(double p_visc) {
  if (p_visc <= 0)
    throw std::invalid_argument("Viscosity has to be >0.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.viscosity = static_cast<float>(p_visc);
    lb_lbfluid_on_lb_params_change(LBPAR_VISCOSITY);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.viscosity = p_visc;
    mpi_bcast_lb_params(LBPAR_VISCOSITY);
#endif // LB
  }
}

double lb_lbfluid_get_viscosity() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return static_cast<double>(lbpar_gpu.viscosity);
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.viscosity;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbfluid_set_bulk_viscosity(double p_bulk_visc) {
  if (p_bulk_visc <= 0)
    throw std::invalid_argument("Bulk viscosity has to be >0.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.bulk_viscosity = static_cast<float>(p_bulk_visc);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_on_lb_params_change(LBPAR_BULKVISC);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.bulk_viscosity = p_bulk_visc;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBPAR_BULKVISC);
#endif // LB
  }
}

double lb_lbfluid_get_bulk_viscosity() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.bulk_viscosity;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.bulk_viscosity;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbfluid_set_gamma_odd(double p_gamma_odd) {
  if (fabs(p_gamma_odd) > 1)
    throw std::invalid_argument("Gamma odd has to be <= 1.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.gamma_odd = static_cast<float>(p_gamma_odd);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_on_lb_params_change(0);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.gamma_odd = p_gamma_odd;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(0);
#endif // LB
  }
}

double lb_lbfluid_get_gamma_odd() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.gamma_odd;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.gamma_odd;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbfluid_set_gamma_even(double p_gamma_even) {
  if (fabs(p_gamma_even) > 1)
    throw std::invalid_argument("gamma_even has to be <= 1.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.gamma_even = static_cast<float>(p_gamma_even);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_on_lb_params_change(0);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.gamma_even = p_gamma_even;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(0);
#endif // LB
  }
}

double lb_lbfluid_get_gamma_even() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.gamma_even;
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.gamma_even;
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
  return {};
}

void lb_lbfluid_set_agrid(double agrid) {
  if (agrid <= 0)
    throw std::invalid_argument("agrid has to be > 0.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_set_agrid_gpu(agrid);
    lb_lbfluid_on_lb_params_change(LBPAR_AGRID);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.agrid = agrid;
    mpi_bcast_lb_params(LBPAR_AGRID);
#endif // LB
  }
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.agrid;
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.agrid;
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
  return {};
}

void lb_lbfluid_set_ext_force_density(int component,
                                      const Vector3d &force_density) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    /* external force density is stored in MD units */
    lbpar_gpu.ext_force_density[3 * component + 0] =
        static_cast<float>(force_density[0]);
    lbpar_gpu.ext_force_density[3 * component + 1] =
        static_cast<float>(force_density[1]);
    lbpar_gpu.ext_force_density[3 * component + 2] =
        static_cast<float>(force_density[2]);
    if (force_density[0] != 0 || force_density[1] != 0 ||
        force_density[2] != 0) {
      lbpar_gpu.external_force_density = 1;
    } else {
      lbpar_gpu.external_force_density = 0;
    }
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);

#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.ext_force_density = force_density;
    mpi_bcast_lb_params(LBPAR_EXTFORCE);
#endif // LB
  }
}

const Vector3d lb_lbfluid_get_ext_force_density() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return {{lbpar_gpu.ext_force_density[0], lbpar_gpu.ext_force_density[1],
             lbpar_gpu.ext_force_density[2]}};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.ext_force_density;
#endif // LB
  }
  return {};
}

void lb_lbfluid_set_tau(double p_tau) {
  if (p_tau <= 0)
    throw std::invalid_argument("tau has to be positive.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.tau = static_cast<float>(p_tau);
    lb_lbfluid_on_lb_params_change(0);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.tau = p_tau;
    mpi_bcast_lb_params(0);
#endif // LB
  }
}

double lb_lbfluid_get_tau() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.tau;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.tau;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbfluid_set_lattice_switch(int local_lattice_switch) {
  switch (local_lattice_switch) {
  case 0:
    lattice_switch = LATTICE_OFF;
    mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
    break;
  case 1:
    lattice_switch = LATTICE_LB;
    mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
    break;
  case 2:
    lattice_switch = LATTICE_LB_GPU;
    mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
    break;
  default:
    throw std::invalid_argument("Invalid lattice switch.");
  }
}

void lb_lbfluid_set_kT(double kT) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.kT = kT;
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lbpar.kT = kT;
    mpi_bcast_lb_params(LBPAR_KT);
#endif
  }
}

double lb_lbfluid_get_kT() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return static_cast<double>(lbpar_gpu.kT);
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.kT;
#endif
  }
  return {};
}

void lb_lbfluid_print_vtk_boundary(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
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
          auto boundary = lb_lbnode_get_boundary(pos);
          fprintf(fp, "%d \n", boundary);
        }
      }
    }
#endif // LB
  }
  fclose(fp);
}

void lb_lbfluid_print_vtk_velocity(const std::string &filename,
                                   std::vector<int> bb1, std::vector<int> bb2) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
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
          auto u = lb_lbnode_get_u(pos);
          fprintf(fp, "%f %f %f\n", u[0], u[1], u[2]);
        }
#endif // LB
  }
  fclose(fp);
}

void lb_lbfluid_print_boundary(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
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
    Vector3i gridsize;

    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];

    for (pos[2] = 0; pos[2] < gridsize[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < gridsize[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
          auto boundary = lb_lbnode_get_boundary(pos);
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
}

void lb_lbfluid_print_velocity(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
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
    Vector3i gridsize;

    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];

    for (pos[2] = 0; pos[2] < gridsize[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < gridsize[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
          auto u = lb_lbnode_get_u(pos);
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
}

void lb_lbfluid_save_checkpoint(const std::string &filename, int binary) {
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
          auto pop = lb_lbnode_get_pop(ind);
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
}

void lb_lbfluid_load_checkpoint(const std::string &filename, int binary) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    FILE *cpfile;
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error("Could not open file for reading.");
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
        if (fread(host_checkpoint_boundary.data(), sizeof(int),
                  int(lbpar_gpu.number_of_nodes),
                  cpfile) != (unsigned int)lbpar_gpu.number_of_nodes) {
          fclose(cpfile);
        }
      if (fread(host_checkpoint_force.data(), sizeof(lbForceFloat),
                3 * int(lbpar_gpu.number_of_nodes),
                cpfile) != (unsigned int)(3 * lbpar_gpu.number_of_nodes)) {
        fclose(cpfile);
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
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error("Could not open file for reading.");
    }
    Vector<19, double> pop;
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
            }
          } else {
            if (fread(pop.data(), sizeof(double), 19, cpfile) != 19)
              throw std::runtime_error("Error reading file.");
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
  }
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

double lb_lbnode_get_density(const Vector3i &ind) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    static LB_rho_v_pi_gpu *host_print_values = nullptr;

    if (host_print_values == nullptr)
      host_print_values =
          (LB_rho_v_pi_gpu *)Utils::malloc(sizeof(LB_rho_v_pi_gpu));
    lb_print_node_GPU(single_nodeindex, host_print_values);
    return host_print_values->rho;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;
    double rho;
    double j[3];
    double pi[6];

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_recv_fluid(node, index, &rho, j, pi);
    // unit conversion
    rho *= 1 / lbpar.agrid / lbpar.agrid / lbpar.agrid;
    return rho;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

const Vector3d lb_lbnode_get_u(const Vector3i &ind) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    static LB_rho_v_pi_gpu *host_print_values = nullptr;
    if (host_print_values == nullptr)
      host_print_values =
          (LB_rho_v_pi_gpu *)Utils::malloc(sizeof(LB_rho_v_pi_gpu));

    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_print_node_GPU(single_nodeindex, host_print_values);
    return {{host_print_values->v[0], host_print_values->v[1],
             host_print_values->v[2]}};
#else
    return {};
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;
    auto ind_shifted = ind;
    double rho;
    Vector3d j;
    Vector<6, double> pi;

    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());
    // unit conversion
    return j / rho * lbpar.agrid / lbpar.tau;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

const Vector<6, double> lb_lbnode_get_pi(const Vector3i &ind) {
  double p0 = 0;
  Vector<6, double> p_pi = lb_lbnode_get_pi_neq(ind);

  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    p0 += lbpar_gpu.rho / lbpar_gpu.agrid / lbpar_gpu.tau / lbpar_gpu.tau / 3.;
#endif // LB_GPU
  } else {
#ifdef LB
    p0 = lbpar.rho / lbpar.agrid / lbpar.tau / lbpar.tau / 3.;
#endif // LB
  }

  p_pi[0] += p0;
  p_pi[2] += p0;
  p_pi[5] += p0;

  return p_pi;
}

const Vector<6, double> lb_lbnode_get_pi_neq(const Vector3i &ind) {
  Vector<6, double> p_pi{};
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
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;
    double rho;
    double j[3];
    Vector<6, double> pi{};

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j, pi.data());
    // unit conversion
    p_pi = pi / lbpar.tau / lbpar.tau / lbpar.agrid;
#endif // LB
  }
  return p_pi;
}

int lb_lbnode_get_boundary(const Vector3i &ind) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    unsigned int host_flag;
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_get_boundary_flag_GPU(single_nodeindex, &host_flag);
    return host_flag;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;
    auto ind_shifted = ind;

    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    int p_boundary;
    mpi_recv_fluid_boundary_flag(node, index, &p_boundary);
    return p_boundary;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

const Vector<19, double> lb_lbnode_get_pop(const Vector3i &ind) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float population[19];

    lb_lbfluid_get_population(ind, population);
    Vector<19, double> p_pop;
    for (int i = 0; i < LBQ; ++i)
      p_pop[i] = population[i];
    return p_pop;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;
    auto ind_shifted = ind;

    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    Vector<19, double> p_pop;
    mpi_recv_fluid_populations(node, index, p_pop.data());
    return p_pop;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbnode_set_rho(const Vector3i &ind, double p_rho) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float host_rho;
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    host_rho = static_cast<float>(p_rho);
    lb_set_node_rho_GPU(single_nodeindex, host_rho);
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;
    double rho;
    Vector3d j;
    Vector<6, double> pi;

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());
    mpi_send_fluid(node, index, p_rho, j, pi);
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbnode_set_u(const Vector3i &ind, const Vector3d &u) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float host_velocity[3];
    host_velocity[0] =
        static_cast<float>(u[0]) * lbpar_gpu.tau / lbpar_gpu.agrid;
    host_velocity[1] =
        static_cast<float>(u[1]) * lbpar_gpu.tau / lbpar_gpu.agrid;
    host_velocity[2] =
        static_cast<float>(u[2]) * lbpar_gpu.tau / lbpar_gpu.agrid;
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_set_node_velocity_GPU(single_nodeindex, host_velocity);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::index_t index;
    int node;
    double rho;
    Vector3d j;
    Vector<6, double> pi;

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted);
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
}

void lb_lbnode_set_pop(const Vector3i &ind, const Vector<19, double> &p_pop) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    float population[19];

    for (int i = 0; i < LBQ; ++i)
      population[i] = p_pop[i];

    lb_lbfluid_set_population(ind, population);
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    Lattice::index_t index;
    int node;

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_send_fluid_populations(node, index, p_pop);
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

#ifdef LB
const Lattice &lb_lbfluid_get_lattice() { return lblattice; }
#endif

int lb_lbfluid_get_lattice_switch() { return lattice_switch; }

void lb_lbfluid_on_lb_params_change(int field) {
  if (field == LBPAR_AGRID) {
#ifdef LB
    if (lattice_switch & LATTICE_LB)
      lb_init();
#endif
#ifdef LB_GPU
    if (lattice_switch & LATTICE_LB_GPU)
      lb_init_gpu();
#endif
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    LBBoundaries::lb_init_boundaries();
#endif
  }
  if (field == LBPAR_DENSITY) {
    lb_lbfluid_reinit_fluid();
  }
  lb_lbfluid_reinit_parameters();
}

#endif

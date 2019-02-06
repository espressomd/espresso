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

void lb_update() {
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

void lb_on_integration_start() {
#ifdef LB
  lb_sanity_checks();

  halo_communication(&update_halo_comm,
                     reinterpret_cast<char *>(lbfluid[0].data()));
#endif
}

/** (Re-)initialize the fluid. */
void lb_lbfluid_reinit_parameters() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lb_reinit_parameters_gpu();
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_reinit_parameters();
#endif
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

/** (Re-)initialize the fluid according to the given value of rho. */
void lb_reinit_fluid() {
#ifdef LB
  std::fill(lbfields.begin(), lbfields.end(), LB_FluidNode());
  /* default values for fields in lattice units */
  Vector3d j{};
  Vector<6, double> pi{};

  LB_TRACE(fprintf(stderr,
                   "Initialising the fluid with equilibrium populations\n"););

  for (Lattice::index_t index = 0; index < lblattice.halo_grid_volume;
       ++index) {
    // calculate equilibrium distribution
    lb_calc_n_from_rho_j_pi(index, lbpar.rho, j, pi);

#ifdef LB_BOUNDARIES
    lbfields[index].boundary = 0;
#endif // LB_BOUNDARIES
  }

#ifdef LB_BOUNDARIES
  LBBoundaries::lb_init_boundaries();
#endif // LB_BOUNDARIES
#endif
}

/** Perform a full initialization of the lattice Boltzmann system.
 *  All derived parameters and the fluid are reset to their default values.
 */
void lb_init() {
#ifdef LB
  LB_TRACE(printf("Begin initialzing fluid on CPU\n"));

  if (lbpar.agrid <= 0.0) {
    runtimeErrorMsg()
        << "Lattice Boltzmann agrid not set when initializing fluid";
  }

  if (check_runtime_errors())
    return;

  Vector3d temp_agrid, temp_offset;
  for (int i = 0; i < 3; i++) {
    temp_agrid[i] = lbpar.agrid;
    temp_offset[i] = 0.5;
  }

  /* initialize the local lattice domain */
  lblattice.init(temp_agrid.data(), temp_offset.data(), 1, 0);

  if (check_runtime_errors())
    return;

  /* allocate memory for data structures */
  lb_realloc_fluid();

  /* prepare the halo communication */
  lb_prepare_communication();

  /* initialize derived parameters */
  lb_lbfluid_reinit_parameters();

  /* setup the initial particle velocity distribution */
  lb_reinit_fluid();

  LB_TRACE(printf("Initialzing fluid on CPU successful\n"));
#endif
}

#ifdef LB
int transfer_momentum = 0;
#endif

#ifdef LB_GPU
int transfer_momentum_gpu = 0;
#endif

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

/** Calculate the fluid velocity at a given position of the lattice.
 *  Note that it can lead to undefined behavior if the position is not
 *  within the local lattice. This version of the function can be called
 *  without the position needing to be on the local processor. Note that this
 *  gives a slightly different version than the values used to couple to MD
 *  beads when near a wall, see lb_lbfluid_get_interpolated_velocity.
 */
int lb_lbfluid_get_interpolated_velocity_global(Vector3d &p, double *v) {
  Vector<6, double>
      delta{}; // velocity field, relative positions to surrounding nodes
  Vector3i ind{}, tmpind{}; // node indices
  int x, y, z;              // counters

  // convert the position into lower left grid point
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    Lattice::map_position_to_lattice_global(p, ind, delta.data(),
                                            lbpar_gpu.agrid);
#endif // LB_GPU
  } else {
#ifdef LB
    Lattice::map_position_to_lattice_global(p, ind, delta.data(), lbpar.agrid);
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

        const auto local_v = lb_lbnode_get_u(tmpind);

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

void lb_lbfluid_set_visc(double p_visc) {
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

double lb_lbfluid_get_visc() {
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

void lb_lbfluid_set_bulk_visc(double p_bulk_visc) {
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

double lb_lbfluid_get_bulk_visc() {
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
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.gamma_even;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbfluid_set_friction(double p_friction) {
  if (p_friction <= 0)
    throw std::invalid_argument("friction has to be > 0.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.friction = static_cast<float>(p_friction);
    lb_lbfluid_on_lb_params_change(LBPAR_FRICTION);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.friction = p_friction;
    mpi_bcast_lb_params(LBPAR_FRICTION);
#endif // LB
  }
}

double lb_lbfluid_get_friction() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.friction;
#else
    return {};
#endif // LB_GPU
  } else {
#ifdef LB
    return lbpar.friction;
#else
    return {};
#endif // LB
  }
}

void lb_lbfluid_set_couple_flag(int couple_flag) {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    if (couple_flag != LB_COUPLE_TWO_POINT &&
        couple_flag != LB_COUPLE_THREE_POINT)
      throw std::invalid_argument("Invalid couple flag.");
    lbpar_gpu.lb_couple_switch = couple_flag;
#endif // LB_GPU
  } else {
#ifdef LB
    /* Only the two point nearest neighbor coupling is present in the case of
       the cpu, so just throw an error if something else is tried */
    if (couple_flag != LB_COUPLE_TWO_POINT)
      throw std::invalid_argument("Invalid couple flag.");
#endif // LB
  }
}

int lb_lbfluid_get_couple_flag() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.lb_couple_switch;
#else
    return {};
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return LB_COUPLE_TWO_POINT;
#else
    return {};
#endif
  } else {
    return LB_COUPLE_NULL;
  }
}

void lb_lbfluid_set_agrid(double p_agrid) {
  if (p_agrid <= 0)
    throw std::invalid_argument("agrid has to be > 0.");
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    lbpar_gpu.agrid = static_cast<float>(p_agrid);

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
    lb_lbfluid_on_lb_params_change(LBPAR_AGRID);
#endif // LB_GPU
  } else {
#ifdef LB
    lbpar.agrid = p_agrid;
    mpi_bcast_lb_params(LBPAR_AGRID);
#endif // LB
  }
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lbpar_gpu.agrid;
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.agrid;
#else
    return {};
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
  }
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
#else
    return {};
#endif // LB_GPU
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lbpar.ext_force_density;
#else
    return {};
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

void lb_set_lattice_switch(int local_lattice_switch) {
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

double lb_lbnode_get_rho(const Vector3i &ind) {
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
    Vector3i grid, ind_shifted = ind;
    double rho;
    Vector3d j;
    Vector<6, double> pi;

    node = lblattice.map_lattice_to_node(ind_shifted.data(), grid.data());
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
    int node, grid[3], ind_shifted[3];
    double rho;
    double j[3];
    Vector<6, double> pi{};

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
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
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
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
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
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
    int node, grid[3], ind_shifted[3];
    double rho;
    Vector3d j;
    Vector<6, double> pi;

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
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
    int node, grid[3], ind_shifted[3];
    double rho;
    Vector3d j;
    Vector<6, double> pi;

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
    int node, grid[3], ind_shifted[3];

    ind_shifted[0] = ind[0];
    ind_shifted[1] = ind[1];
    ind_shifted[2] = ind[2];
    node = lblattice.map_lattice_to_node(ind_shifted, grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_send_fluid_populations(node, index, p_pop);
#endif // LB
  } else {
    throw std::runtime_error("LB not activated.");
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

namespace {
#ifdef LB
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
#endif
} // namespace

/*
 * @brief Interpolate the fluid velocity.
 *
 * @param pos Position
 * @param v Interpolated velocity in MD units.
 */
#ifdef LB
const Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d &pos) {
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
#endif

#ifdef LB
void lb_lbfluid_add_force_density(const Vector3d &pos,
                                  const Vector3d &force_density) {
  lattice_interpolation(lblattice, pos,
                        [&force_density](Lattice::index_t index, double w) {
                          auto &node = lbfields[index];

                          node.force_density[0] += w * force_density[0];
                          node.force_density[1] += w * force_density[1];
                          node.force_density[2] += w * force_density[2];
                        });
}
#endif

#ifdef LB
const Lattice &lb_lbfluid_get_lattice() { return lblattice; }
#endif

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
#ifdef LB
    if (lattice_switch & LATTICE_LB)
      lb_reinit_fluid();
#endif
#ifdef LB_GPU
    if (lattice_switch & LATTICE_LB_GPU)
      lb_reinit_fluid_gpu();
#endif
  }
  lb_lbfluid_reinit_parameters();
}

#endif

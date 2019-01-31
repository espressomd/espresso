#include <fstream>

#include "lb_interface.hpp"
#include "lb.hpp"
#include "lbgpu.hpp"
#include "communication.hpp"
#include "global.hpp"
#include "grid.hpp"

#if defined(LB) || defined(LB_GPU)

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

#endif // SHANCHEN

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


/** Calculate the fluid velocity at a given position of the lattice.
 *  Note that it can lead to undefined behavior if the position is not
 *  within the local lattice. This version of the function can be called
 *  without the position needing to be on the local processor. Note that this
 *  gives a slightly different version than the values used to couple to MD
 *  beads when near a wall, see lb_lbfluid_get_interpolated_velocity.
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


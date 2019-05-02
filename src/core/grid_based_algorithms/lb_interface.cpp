#include "lb_interface.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "electrokinetics.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "lb.hpp"
#include "lbgpu.hpp"

#include <utils/index.hpp>
using Utils::get_linear_index;

#include <fstream>

ActiveLB lattice_switch = ActiveLB::NONE;

/* LB CPU callback interface */
namespace {
/** Issue REQ_SEND_FLUID: Send a single lattice site to a processor.
 *  @param node   processor to send to
 *  @param index  index of the lattice site
 *  @param rho    local fluid density
 *  @param j      local fluid velocity
 *  @param pi     local fluid pressure
 */
void mpi_send_fluid(int node, int index, double rho, Utils::Vector3d const &j,
                    Utils::Vector6d const &pi) {
  if (node == this_node) {
    lb_calc_n_from_rho_j_pi(index, rho, j, pi);
  } else if (0 == this_node) {
    mpi_call(mpi_send_fluid, node, index, rho, j, pi);
  }
}

REGISTER_CALLBACK(mpi_send_fluid)

void mpi_recv_fluid_slave(int node, int index) {
  if (node == this_node) {
    double data[10];
    lb_calc_local_fields(index, &data[0], &data[1], &data[4]);
    MPI_Send(data, 10, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
  }
}

REGISTER_CALLBACK(mpi_recv_fluid_slave)

/** Issue REQ_GET_FLUID: Receive a single lattice site from a processor.
 *  @param node   processor to send to
 *  @param index  index of the lattice site
 *  @param rho    local fluid density
 *  @param j      local fluid velocity
 *  @param pi     local fluid pressure
 */
void mpi_recv_fluid(int node, int index, double *rho, double *j, double *pi) {
  if (node == this_node) {
    lb_calc_local_fields(index, rho, j, pi);
  } else {
    double data[10];
    mpi_call(mpi_recv_fluid_slave, node, index);
    MPI_Recv(data, 10, MPI_DOUBLE, node, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
    *rho = data[0];
    j[0] = data[1];
    j[1] = data[2];
    j[2] = data[3];
    pi[0] = data[4];
    pi[1] = data[5];
    pi[2] = data[6];
    pi[3] = data[7];
    pi[4] = data[8];
    pi[5] = data[9];
  }
}

void mpi_send_fluid_populations_slave(int node, int index) {
  if (node == this_node) {
    Utils::Vector19d populations;
    MPI_Recv(populations.data(), 19, MPI_DOUBLE, 0, SOME_TAG, comm_cart,
             MPI_STATUS_IGNORE);
    lb_set_populations(index, populations);
  }
}

REGISTER_CALLBACK(mpi_send_fluid_populations_slave)

/** Issue REQ_SEND_FLUID_POPULATIONS: Send a single lattice site to a processor.
 *  @param node   processor to send to
 *  @param index  index of the lattice site
 *  @param pop    local fluid population
 */
void mpi_send_fluid_populations(int node, int index,
                                const Utils::Vector19d &pop) {
  if (node == this_node) {
    lb_set_populations(index, pop);
  } else {
    mpi_call(mpi_send_fluid_populations_slave, node, index);
    MPI_Send(pop.data(), 19, MPI_DOUBLE, node, SOME_TAG, comm_cart);
  }
}

void mpi_bcast_lb_params_slave(LBParam field, const LB_Parameters &params_) {
  lbpar = params_;
  lb_lbfluid_on_lb_params_change(field);
}

REGISTER_CALLBACK(mpi_bcast_lb_params_slave)

/** @brief Broadcast a parameter for lattice Boltzmann.
 *  @param[in] field  References the parameter field to be broadcasted.
 *                    The references are defined in lb.hpp
 */
void mpi_bcast_lb_params(LBParam field) {
  mpi_call(mpi_bcast_lb_params_slave, field, lbpar);
  lb_lbfluid_on_lb_params_change(field);
}

void mpi_recv_fluid_boundary_flag_slave(int node, int index) {
#ifdef LB_BOUNDARIES
  if (node == this_node) {
    int data;
    lb_local_fields_get_boundary_flag(index, &data);
    MPI_Send(&data, 1, MPI_INT, 0, SOME_TAG, comm_cart);
  }
#endif
}

REGISTER_CALLBACK(mpi_recv_fluid_boundary_flag_slave)

/** Issue REQ_LB_GET_BOUNDARY_FLAG: Receive a single lattice sites boundary
 *  flag from a processor.
 *  @param node      processor to send to
 *  @param index     index of the lattice site
 *  @param boundary  local boundary flag
 */
void mpi_recv_fluid_boundary_flag(int node, int index, int *boundary) {
#ifdef LB_BOUNDARIES
  if (node == this_node) {
    lb_local_fields_get_boundary_flag(index, boundary);
  } else {
    int data = 0;
    mpi_call(mpi_recv_fluid_boundary_flag_slave, node, index);
    MPI_Recv(&data, 1, MPI_INT, node, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
    *boundary = data;
  }
#endif
}

void mpi_recv_fluid_populations_slave(int node, int index) {
  if (node == this_node) {
    double data[19];
    lb_get_populations(index, data);
    MPI_Send(data, 19, MPI_DOUBLE, 0, SOME_TAG, comm_cart);
  }
}

REGISTER_CALLBACK(mpi_recv_fluid_populations_slave)

/** Issue REQ_RECV_FLUID_POPULATIONS: Send a single lattice site to a processor.
 *  @param node   processor to send to
 *  @param index  index of the lattice site
 *  @param pop    local fluid population
 */
void mpi_recv_fluid_populations(int node, int index, double *pop) {
  if (node == this_node) {
    lb_get_populations(index, pop);
  } else {
    mpi_call(mpi_recv_fluid_populations_slave, node, index);
    MPI_Recv(pop, 19, MPI_DOUBLE, node, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  }
}
} // namespace

void lb_lbfluid_update() {
  if (lattice_switch == ActiveLB::CPU) {
    lattice_boltzmann_update();
  } else if (lattice_switch == ActiveLB::GPU and this_node == 0) {
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
  if (lb_lbfluid_get_kT() > 0.0) {
    if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
      rng_counter_fluid_gpu->increment();
#endif
    } else if (lattice_switch == ActiveLB::CPU) {
      rng_counter_fluid->increment();
    }
  }
}

/**
 * @brief Check the boundary velocities.
 * Sanity check if the velocity defined at LB boundaries is within the Mach
 * number limits of the scheme i.e. u < 0.3.
 */
void lb_boundary_mach_check() {
  // Boundary velocities are stored in MD units, therefore we need to scale them
  // in order to get lattice units.
  auto const conv_fac = lb_lbfluid_get_tau() / lb_lbfluid_get_agrid();
  double constexpr mach_limit = 0.3;
  using LBBoundaries::lbboundaries;
  if (std::any_of(lbboundaries.begin(), lbboundaries.end(),
                  [conv_fac, mach_limit](auto const &b) {
                    return (b->velocity() * conv_fac).norm() >= mach_limit;
                  })) {
    runtimeErrorMsg() << "Lattice velocity exceeds the Mach number limit";
  }
}

/**
 * @brief Perform LB parameter and boundary velocity checks.
 */
void lb_lbfluid_sanity_checks() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    if (this_node == 0) {
      lb_GPU_sanity_checks();
      lb_boundary_mach_check();
    }
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_sanity_checks();
    lb_boundary_mach_check();
  }
}

void lb_lbfluid_on_integration_start() {
  lb_lbfluid_sanity_checks();
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    if (this_node == 0 and lb_reinit_particles_gpu()) {
      lb_realloc_particles_gpu();
      lb_reinit_particles_gpu.validate();
    }
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    halo_communication(&update_halo_comm,
                       reinterpret_cast<char *>(lbfluid[0].data()));
  }
}

void lb_lbfluid_invalidate_particle_allocation() {
#ifdef LB_GPU
  lb_reinit_particles_gpu.invalidate();
#endif
}

/** (Re-)initialize the fluid. */
void lb_lbfluid_reinit_parameters() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    if (this_node == 0)
      lb_reinit_parameters_gpu();
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_reinit_parameters();
  }
}

/** (Re-)initialize the fluid according to the value of rho. */

void lb_lbfluid_reinit_fluid() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lb_reinit_fluid_gpu();
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_reinit_fluid();
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

/** Perform a full initialization of the lattice Boltzmann system.
 *  All derived parameters and the fluid are reset to their default values.
 */
void lb_lbfluid_init() {
  if (lattice_switch == ActiveLB::GPU && this_node == 0) {
#ifdef LB_GPU
    lb_init_gpu();
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lb_init();
  }
}

uint64_t lb_lbfluid_get_rng_state() {
  if (lattice_switch == ActiveLB::CPU) {
    return lb_fluid_get_rng_state();
  }
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lb_fluid_get_rng_state_gpu();
#endif
  }
  return {};
}

void lb_lbfluid_set_rng_state(uint64_t counter) {
  if (lattice_switch == ActiveLB::CPU) {
    lb_fluid_set_rng_state(counter);
  } else if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lb_fluid_set_rng_state_gpu(counter);
#endif
  }
}

void lb_lbfluid_set_density(double p_dens) {
  if (p_dens <= 0)
    throw std::invalid_argument("Density has to be > 0.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.rho = static_cast<float>(p_dens);
    lb_lbfluid_on_lb_params_change(LBParam::DENSITY);
#endif // LB_GPU
  } else {
    lbpar.rho = p_dens;
    mpi_bcast_lb_params(LBParam::DENSITY);
  }
}

double lb_lbfluid_get_density() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return static_cast<double>(lbpar_gpu.rho);
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.rho;
  }
  throw std::runtime_error("LB not activated.");
}

void lb_lbfluid_set_viscosity(double p_visc) {
  if (p_visc <= 0)
    throw std::invalid_argument("Viscosity has to be >0.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.viscosity = static_cast<float>(p_visc);
    lb_lbfluid_on_lb_params_change(LBParam::VISCOSITY);
#endif // LB_GPU
  } else {
    lbpar.viscosity = p_visc;
    mpi_bcast_lb_params(LBParam::VISCOSITY);
  }
}

double lb_lbfluid_get_viscosity() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return static_cast<double>(lbpar_gpu.viscosity);
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.viscosity;
  }
  throw std::runtime_error("LB not activated.");
}

void lb_lbfluid_set_bulk_viscosity(double p_bulk_visc) {
  if (p_bulk_visc <= 0)
    throw std::invalid_argument("Bulk viscosity has to be >0.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.bulk_viscosity = static_cast<float>(p_bulk_visc);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_on_lb_params_change(LBParam::BULKVISC);
#endif // LB_GPU
  } else {
    lbpar.bulk_viscosity = p_bulk_visc;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBParam::BULKVISC);
  }
}

double lb_lbfluid_get_bulk_viscosity() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lbpar_gpu.bulk_viscosity;
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.bulk_viscosity;
  }
  throw std::runtime_error("LB not activated.");
}

void lb_lbfluid_set_gamma_odd(double p_gamma_odd) {
  if (fabs(p_gamma_odd) > 1)
    throw std::invalid_argument("Gamma odd has to be <= 1.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.gamma_odd = static_cast<float>(p_gamma_odd);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_on_lb_params_change(LBParam::DENSITY);
#endif // LB_GPU
  } else {
    lbpar.gamma_odd = p_gamma_odd;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBParam::DENSITY);
  }
}

double lb_lbfluid_get_gamma_odd() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lbpar_gpu.gamma_odd;
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.gamma_odd;
  }
  throw std::runtime_error("LB not activated.");
}

void lb_lbfluid_set_gamma_even(double p_gamma_even) {
  if (fabs(p_gamma_even) > 1)
    throw std::invalid_argument("gamma_even has to be <= 1.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.gamma_even = static_cast<float>(p_gamma_even);
    lbpar_gpu.is_TRT = false;
    lb_lbfluid_on_lb_params_change(LBParam::DENSITY);
#endif // LB_GPU
  } else {
    lbpar.gamma_even = p_gamma_even;
    lbpar.is_TRT = false;
    mpi_bcast_lb_params(LBParam::DENSITY);
  }
}

double lb_lbfluid_get_gamma_even() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lbpar_gpu.gamma_even;
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.gamma_even;
  }
  throw std::runtime_error("LB not activated.");

  return {};
}

void lb_lbfluid_set_agrid(double agrid) {
  if (agrid <= 0)
    throw std::invalid_argument("agrid has to be > 0.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lb_set_agrid_gpu(agrid);
    lb_lbfluid_on_lb_params_change(LBParam::AGRID);
#endif // LB_GPU
  } else {
    lbpar.agrid = agrid;
    mpi_bcast_lb_params(LBParam::AGRID);
  }
}

double lb_lbfluid_get_agrid() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lbpar_gpu.agrid;
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.agrid;
  }
  throw std::runtime_error("LB not activated.");

  return {};
}

void lb_lbfluid_set_ext_force_density(const Utils::Vector3d &force_density) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.ext_force_density[0] = static_cast<float>(force_density[0]);
    lbpar_gpu.ext_force_density[1] = static_cast<float>(force_density[1]);
    lbpar_gpu.ext_force_density[2] = static_cast<float>(force_density[2]);
    if (force_density[0] != 0 || force_density[1] != 0 ||
        force_density[2] != 0) {
      lbpar_gpu.external_force_density = 1;
    } else {
      lbpar_gpu.external_force_density = 0;
    }
    lb_reinit_extern_nodeforce_GPU(&lbpar_gpu);

#endif // LB_GPU
  } else {
    lbpar.ext_force_density = force_density;
    mpi_bcast_lb_params(LBParam::EXTFORCE);
  }
}

const Utils::Vector3d lb_lbfluid_get_ext_force_density() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return {{lbpar_gpu.ext_force_density[0], lbpar_gpu.ext_force_density[1],
             lbpar_gpu.ext_force_density[2]}};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.ext_force_density;
  }
  return {};
}

void lb_lbfluid_set_tau(double p_tau) {
  if (p_tau <= 0)
    throw std::invalid_argument("tau has to be positive.");
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.tau = static_cast<float>(p_tau);
    lb_lbfluid_on_lb_params_change(LBParam::DENSITY);
#endif // LB_GPU
  } else {
    lbpar.tau = p_tau;
    mpi_bcast_lb_params(LBParam::DENSITY);
  }
}

double lb_lbfluid_get_tau() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lbpar_gpu.tau;
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.tau;
  }
  throw std::runtime_error("LB not activated.");
}

void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) {
  switch (local_lattice_switch) {
  case ActiveLB::NONE:
  case ActiveLB::CPU:
  case ActiveLB::GPU:
    break;
  default:
    throw std::invalid_argument("Invalid lattice switch.");
  }
  lattice_switch = local_lattice_switch;
  mpi_bcast_parameter(FIELD_LATTICE_SWITCH);
}

void lb_lbfluid_set_kT(double kT) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lbpar_gpu.kT = kT;
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    lbpar.kT = kT;
    mpi_bcast_lb_params(LBParam::KT);
  }
}

double lb_lbfluid_get_kT() {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return static_cast<double>(lbpar_gpu.kT);
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    return lbpar.kT;
  }
  return {};
}

double lb_lbfluid_get_lattice_speed() {
  return lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
}

void lb_lbfluid_print_vtk_boundary(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  if (lattice_switch == ActiveLB::GPU) {
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
    Utils::Vector3i pos;
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
      if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
        bb_high = {static_cast<int>(lbpar_gpu.dim_x) - 1,
                   static_cast<int>(lbpar_gpu.dim_y) - 1,
                   static_cast<int>(lbpar_gpu.dim_z) - 1};
#endif // LB_GPU
      } else {
        bb_high = {lblattice.global_grid[0] - 1, lblattice.global_grid[1] - 1,
                   lblattice.global_grid[2] - 1};
      }
      break;
    }

    bb_low.push_back(std::min(*val1, *val2));
    bb_high.push_back(std::max(*val1, *val2));
  }

  Utils::Vector3i pos;
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
    host_values = (LB_rho_v_pi_gpu *)Utils::malloc(size_of_values);
    lb_get_values_GPU(host_values);
    auto const lattice_speed = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
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
          fprintf(fp, "%f %f %f\n", host_values[j].v[0] * lattice_speed,
                  host_values[j].v[1] * lattice_speed,
                  host_values[j].v[2] * lattice_speed);
        }
    free(host_values);
#endif // LB_GPU
  } else {
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
          auto u = lb_lbnode_get_velocity(pos) * lb_lbfluid_get_lattice_speed();
          fprintf(fp, "%f %f %f\n", u[0], u[1], u[2]);
        }
  }
  fclose(fp);
}

void lb_lbfluid_print_boundary(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    unsigned int *bound_array;
    bound_array = (unsigned int *)Utils::malloc(lbpar_gpu.number_of_nodes *
                                                sizeof(unsigned int));
    lb_get_boundary_flags_GPU(bound_array);

    Utils::Vector3i xyz;
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
    Utils::Vector3i pos;
    Utils::Vector3i gridsize;

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
  }
  fclose(fp);
}

void lb_lbfluid_print_velocity(const std::string &filename) {
  FILE *fp = fopen(filename.c_str(), "w");

  if (fp == nullptr) {
    throw std::runtime_error("Could not open file for writing.");
  }

  auto const lattice_speed = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau();
  auto const agrid = lb_lbfluid_get_agrid();
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    size_t size_of_values = lbpar_gpu.number_of_nodes * sizeof(LB_rho_v_pi_gpu);
    host_values = (LB_rho_v_pi_gpu *)Utils::malloc(size_of_values);
    lb_get_values_GPU(host_values);
    Utils::Vector3i xyz;
    int j;
    for (j = 0; j < int(lbpar_gpu.number_of_nodes); ++j) {
      xyz[0] = j % lbpar_gpu.dim_x;
      int k = j / lbpar_gpu.dim_x;
      xyz[1] = k % lbpar_gpu.dim_y;
      k /= lbpar_gpu.dim_y;
      xyz[2] = k;
      /** print of the calculated phys values */
      fprintf(fp, "%f %f %f %f %f %f\n", (xyz[0] + 0.5) * agrid,
              (xyz[1] + 0.5) * agrid, (xyz[2] + 0.5) * agrid,
              host_values[j].v[0] * lattice_speed,
              host_values[j].v[1] * lattice_speed,
              host_values[j].v[2] * lattice_speed);
    }
    free(host_values);
#endif // LB_GPU
  } else {
    Utils::Vector3i pos;
    Utils::Vector3i gridsize;

    gridsize[0] = box_l[0] / lblattice.agrid[0];
    gridsize[1] = box_l[1] / lblattice.agrid[1];
    gridsize[2] = box_l[2] / lblattice.agrid[2];

    for (pos[2] = 0; pos[2] < gridsize[2]; pos[2]++) {
      for (pos[1] = 0; pos[1] < gridsize[1]; pos[1]++) {
        for (pos[0] = 0; pos[0] < gridsize[0]; pos[0]++) {
          auto const u = lb_lbnode_get_velocity(pos) * lattice_speed;
          fprintf(fp, "%f %f %f %f %f %f\n",
                  (pos[0] + 0.5) * lblattice.agrid[0],
                  (pos[1] + 0.5) * lblattice.agrid[1],
                  (pos[2] + 0.5) * lblattice.agrid[2], u[0], u[1], u[2]);
        }
      }
    }
  }

  fclose(fp);
}

void lb_lbfluid_save_checkpoint(const std::string &filename, int binary) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    auto *host_checkpoint_vd =
        (float *)Utils::malloc(lbpar_gpu.number_of_nodes * 19 * sizeof(float));
    lb_save_checkpoint_GPU(host_checkpoint_vd);
    if (!binary) {
      std::fstream cpfile(filename, std::ios::out);
      cpfile << std::fixed;
      cpfile.precision(8);
      cpfile << lbpar_gpu.dim_x << " ";
      cpfile << lbpar_gpu.dim_y << " ";
      cpfile << lbpar_gpu.dim_z << "\n";
      for (int n = 0; n < (19 * int(lbpar_gpu.number_of_nodes)); n++) {
        cpfile << host_checkpoint_vd[n] << "\n";
      }
      cpfile.close();
    } else {
      std::fstream cpfile(filename, std::ios::out | std::ios::binary);
      cpfile.write(reinterpret_cast<char *>(&lbpar_gpu.dim_x),
                   sizeof(lbpar_gpu.dim_x));
      cpfile.write(reinterpret_cast<char *>(&lbpar_gpu.dim_y),
                   sizeof(lbpar_gpu.dim_y));
      cpfile.write(reinterpret_cast<char *>(&lbpar_gpu.dim_z),
                   sizeof(lbpar_gpu.dim_z));
      cpfile.write(reinterpret_cast<char *>(host_checkpoint_vd),
                   19 * sizeof(float) * lbpar_gpu.number_of_nodes);
      cpfile.close();
    }
    free(host_checkpoint_vd);
#endif // LB_GPU
  } else if (lattice_switch == ActiveLB::CPU) {
    std::fstream cpfile;
    if (binary) {
      cpfile.open(filename, std::ios::out | std::ios::binary);
    } else {
      cpfile.open(filename, std::ios::out);
      cpfile.precision(16);
      cpfile << std::fixed;
    }

    double pop[19];
    Utils::Vector3i ind;
    auto const gridsize = lblattice.global_grid;

    if (!binary) {
      cpfile << gridsize[0] << " " << gridsize[1] << " " << gridsize[2] << "\n";
    } else {
      cpfile.write(reinterpret_cast<const char *>(gridsize.data()),
                   3 * sizeof(gridsize[0]));
    }

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          ind[0] = i;
          ind[1] = j;
          ind[2] = k;
          auto pop = lb_lbnode_get_pop(ind);
          if (!binary) {
            for (int n = 0; n < 19; n++) {
              cpfile << pop[n] << "\n";
            }
          } else {
            cpfile.write(reinterpret_cast<char *>(&pop[0]),
                         19 * sizeof(double));
          }
        }
      }
    }
    cpfile.close();
  }
}

void lb_lbfluid_load_checkpoint(const std::string &filename, int binary) {
  int res;
  std::string err_msg = "Error while reading LB checkpoint: ";
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    FILE *cpfile;
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error(err_msg + "could not open file for reading.");
    }
    std::vector<float> host_checkpoint_vd(lbpar_gpu.number_of_nodes * 19);
    if (!binary) {
      int saved_gridsize[3];
      for (int &n : saved_gridsize) {
        res = fscanf(cpfile, "%i", &n);
        if (res == EOF) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "EOF found.");
        }
        if (res != 1) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "incorrectly formatted data.");
        }
      }
      if (saved_gridsize[0] != lbpar_gpu.dim_x ||
          saved_gridsize[1] != lbpar_gpu.dim_y ||
          saved_gridsize[2] != lbpar_gpu.dim_z) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "grid dimensions mismatch, read [" +
                                 std::to_string(saved_gridsize[0]) + ' ' +
                                 std::to_string(saved_gridsize[1]) + ' ' +
                                 std::to_string(saved_gridsize[2]) +
                                 "], expected [" +
                                 std::to_string(lbpar_gpu.dim_x) + ' ' +
                                 std::to_string(lbpar_gpu.dim_y) + ' ' +
                                 std::to_string(lbpar_gpu.dim_z) + "].");
      }
      for (int n = 0; n < (19 * int(lbpar_gpu.number_of_nodes)); n++) {
        res = fscanf(cpfile, "%f", &host_checkpoint_vd[n]);
        if (res == EOF) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "EOF found.");
        }
        if (res != 1) {
          fclose(cpfile);
          throw std::runtime_error(err_msg + "incorrectly formatted data.");
        }
      }
    } else {
      int saved_gridsize[3];
      if (fread(&saved_gridsize[0], sizeof(int), 3, cpfile) != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
      if (saved_gridsize[0] != lbpar_gpu.dim_x ||
          saved_gridsize[1] != lbpar_gpu.dim_y ||
          saved_gridsize[2] != lbpar_gpu.dim_z) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "grid dimensions mismatch, read [" +
                                 std::to_string(saved_gridsize[0]) + ' ' +
                                 std::to_string(saved_gridsize[1]) + ' ' +
                                 std::to_string(saved_gridsize[2]) +
                                 "], expected [" +
                                 std::to_string(lbpar_gpu.dim_x) + ' ' +
                                 std::to_string(lbpar_gpu.dim_y) + ' ' +
                                 std::to_string(lbpar_gpu.dim_z) + "].");
      }
      if (fread(host_checkpoint_vd.data(), sizeof(float),
                19 * int(lbpar_gpu.number_of_nodes),
                cpfile) != (unsigned int)(19 * lbpar_gpu.number_of_nodes)) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    }
    if (!binary) {
      // skip spaces
      for (int n = 0; n < 2; ++n) {
        res = fgetc(cpfile);
        if (res != (int)' ' && res != (int)'\n')
          break;
      }
    } else {
      res = fgetc(cpfile);
    }
    if (res != EOF) {
      fclose(cpfile);
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
    fclose(cpfile);
    lb_load_checkpoint_GPU(host_checkpoint_vd.data());
#endif // LB_GPU
  } else if (lattice_switch == ActiveLB::CPU) {
    FILE *cpfile;
    cpfile = fopen(filename.c_str(), "r");
    if (!cpfile) {
      throw std::runtime_error(err_msg + "could not open file for reading.");
    }

    Utils::Vector19d pop;
    Utils::Vector3i ind;
    auto const gridsize = lblattice.global_grid;
    int saved_gridsize[3];
    mpi_bcast_lb_params(LBParam::DENSITY);

    if (!binary) {
      res = fscanf(cpfile, "%i %i %i\n", &saved_gridsize[0], &saved_gridsize[1],
                   &saved_gridsize[2]);
      if (res == EOF) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "EOF found.");
      }
      if (res != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    } else {
      if (fread(&saved_gridsize[0], sizeof(int), 3, cpfile) != 3) {
        fclose(cpfile);
        throw std::runtime_error(err_msg + "incorrectly formatted data.");
      }
    }
    if (saved_gridsize[0] != gridsize[0] || saved_gridsize[1] != gridsize[1] ||
        saved_gridsize[2] != gridsize[2]) {
      fclose(cpfile);
      throw std::runtime_error(err_msg + "grid dimensions mismatch, read [" +
                               std::to_string(saved_gridsize[0]) + ' ' +
                               std::to_string(saved_gridsize[1]) + ' ' +
                               std::to_string(saved_gridsize[2]) +
                               "], expected [" + std::to_string(gridsize[0]) +
                               ' ' + std::to_string(gridsize[1]) + ' ' +
                               std::to_string(gridsize[2]) + "].");
    }

    for (int i = 0; i < gridsize[0]; i++) {
      for (int j = 0; j < gridsize[1]; j++) {
        for (int k = 0; k < gridsize[2]; k++) {
          ind[0] = i;
          ind[1] = j;
          ind[2] = k;
          if (!binary) {
            res = fscanf(cpfile,
                         "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                         "%lf %lf %lf %lf %lf %lf \n",
                         &pop[0], &pop[1], &pop[2], &pop[3], &pop[4], &pop[5],
                         &pop[6], &pop[7], &pop[8], &pop[9], &pop[10], &pop[11],
                         &pop[12], &pop[13], &pop[14], &pop[15], &pop[16],
                         &pop[17], &pop[18]);
            if (res == EOF) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "EOF found.");
            }
            if (res != 19) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
          } else {
            if (fread(pop.data(), sizeof(double), 19, cpfile) != 19) {
              fclose(cpfile);
              throw std::runtime_error(err_msg + "incorrectly formatted data.");
            }
          }
          lb_lbnode_set_pop(ind, pop);
        }
      }
    }
    if (!binary) {
      // skip spaces
      for (int n = 0; n < 2; ++n) {
        res = fgetc(cpfile);
        if (res != (int)' ' && res != (int)'\n')
          break;
      }
    } else {
      res = fgetc(cpfile);
    }
    if (res != EOF) {
      fclose(cpfile);
      throw std::runtime_error(err_msg + "extra data found, expected EOF.");
    }
    fclose(cpfile);
  } else {
    throw std::runtime_error(
        "To load an LB checkpoint one needs to have already "
        "initialized the LB fluid with the same grid size.");
  }
}

bool lb_lbnode_is_index_valid(const Utils::Vector3i &ind) {
  auto within_bounds = [](const Utils::Vector3i &ind,
                          const Utils::Vector3i &limits) {
    return ind < limits && ind >= Utils::Vector3i{};
  };
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return within_bounds(ind, {static_cast<int>(lbpar_gpu.dim_x),
                               static_cast<int>(lbpar_gpu.dim_y),
                               static_cast<int>(lbpar_gpu.dim_z)});
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    return within_bounds(ind, lblattice.global_grid);
  }
  return false;
}

double lb_lbnode_get_density(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
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
  }
  if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;
    double rho;
    double j[3];
    double pi[6];

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_recv_fluid(node, index, &rho, j, pi);
    return rho;
  }
  throw std::runtime_error("LB not activated.");
}

const Utils::Vector3d lb_lbnode_get_velocity(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
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
#endif
  }
  if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;
    auto ind_shifted = ind;
    double rho;
    Utils::Vector3d j;
    Utils::Vector6d pi;

    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());
    return j / rho;
  }
  throw std::runtime_error("LB not activated.");

  return {};
}

const Utils::Vector6d lb_lbnode_get_pi(const Utils::Vector3i &ind) {
  Utils::Vector6d p_pi = lb_lbnode_get_pi_neq(ind);

  // Add equilibrium stress to the diagonal (in LB units)
  double const p0 = lb_lbfluid_get_density() * lbmodel.c_sound_sq;

  p_pi[0] += p0;
  p_pi[2] += p0;
  p_pi[5] += p0;

  return p_pi;
}

const Utils::Vector6d lb_lbnode_get_pi_neq(const Utils::Vector3i &ind) {
  Utils::Vector6d p_pi{};
  if (lattice_switch == ActiveLB::GPU) {
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
  } else if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;
    double rho;
    double j[3];
    Utils::Vector6d pi{};

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j, p_pi.data());
  }
  return p_pi;
}

/** calculates the average stress of all nodes by iterating
 * over all nodes and deviding by the number_of_nodes.
 */
const Utils::Vector6d lb_lbfluid_get_stress() {
  Utils::Vector6d p{0, 0, 0, 0, 0, 0};

  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    // Copy observable data from gpu
    std::vector<LB_rho_v_pi_gpu> host_values(lbpar_gpu.number_of_nodes);
    lb_get_values_GPU(host_values.data());
    std::for_each(host_values.begin(), host_values.end(),
                  [&p](LB_rho_v_pi_gpu &v) {
                    for (int i = 0; i < 6; i++)
                      p[i] += v.pi[i];
                  });

    // Normalize
    p *= (1. / lbpar_gpu.number_of_nodes);

    // Add equilibrium stress to the diagonal (in LB units)
    double const p0 = lb_lbfluid_get_density() * lbmodel.c_sound_sq;

    p[0] += p0;
    p[2] += p0;
    p[5] += p0;

#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    for (int i = 0; i < lblattice.global_grid[0]; i++) {
      for (int j = 0; j < lblattice.global_grid[1]; j++) {
        for (int k = 0; k < lblattice.global_grid[2]; k++) {
          const Utils::Vector3i node{{i, j, k}};
          p += lb_lbnode_get_pi(node);
        }
      }
    }

    int const number_of_nodes = lblattice.global_grid[0] *
                                lblattice.global_grid[1] *
                                lblattice.global_grid[2];

    p *= 1. / number_of_nodes;
  } else {
    throw std::runtime_error("LB method called on inactive LB");
  }
  return p;
}

int lb_lbnode_get_boundary(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    unsigned int host_flag;
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_get_boundary_flag_GPU(single_nodeindex, &host_flag);
    return host_flag;
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;
    auto ind_shifted = ind;

    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    int p_boundary;
    mpi_recv_fluid_boundary_flag(node, index, &p_boundary);
    return p_boundary;
  }
  throw std::runtime_error("LB not activated.");
}

const Utils::Vector19d lb_lbnode_get_pop(const Utils::Vector3i &ind) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    float population[19];

    lb_lbfluid_get_population(ind, population);
    Utils::Vector19d p_pop;
    for (int i = 0; i < LBQ; ++i)
      p_pop[i] = population[i];
    return p_pop;
#else
    return {};
#endif // LB_GPU
  }
  if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;
    auto ind_shifted = ind;

    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    Utils::Vector19d p_pop;
    mpi_recv_fluid_populations(node, index, p_pop.data());
    return p_pop;
  }
  throw std::runtime_error("LB not activated.");
}

void lb_lbnode_set_density(const Utils::Vector3i &ind, double p_rho) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    auto const host_rho = static_cast<float>(p_rho);
    lb_set_node_rho_GPU(single_nodeindex, host_rho);
#endif // LB_GPU
  } else if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;
    double rho;
    Utils::Vector3d j;
    Utils::Vector6d pi;

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());
    mpi_send_fluid(node, index, p_rho, j, pi);
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

void lb_lbnode_set_velocity(const Utils::Vector3i &ind,
                            const Utils::Vector3d &u) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    float host_velocity[3];
    host_velocity[0] = static_cast<float>(u[0]);
    host_velocity[1] = static_cast<float>(u[1]);
    host_velocity[2] = static_cast<float>(u[2]);
    int single_nodeindex = ind[0] + ind[1] * lbpar_gpu.dim_x +
                           ind[2] * lbpar_gpu.dim_x * lbpar_gpu.dim_y;
    lb_set_node_velocity_GPU(single_nodeindex, host_velocity);
#endif // LB_GPU
  } else {
    Lattice::index_t index;
    int node;

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);

    double rho;
    Utils::Vector3d j;
    Utils::Vector6d pi;

    mpi_recv_fluid(node, index, &rho, j.data(), pi.data());

    /* transform to lattice units */
    j = rho * u;
    mpi_send_fluid(node, index, rho, j, pi);
  }
}

void lb_lbnode_set_pop(const Utils::Vector3i &ind,
                       const Utils::Vector19d &p_pop) {
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    float population[19];

    for (int i = 0; i < LBQ; ++i)
      population[i] = p_pop[i];

    lb_lbfluid_set_population(ind, population);
#endif // LB_GPU
  } else if (lattice_switch == ActiveLB::CPU) {
    Lattice::index_t index;
    int node;

    auto ind_shifted = ind;
    node = lblattice.map_lattice_to_node(ind_shifted, node_grid);
    index = get_linear_index(ind_shifted[0], ind_shifted[1], ind_shifted[2],
                             lblattice.halo_grid);
    mpi_send_fluid_populations(node, index, p_pop);
  } else {
    throw std::runtime_error("LB not activated.");
  }
}

const Lattice &lb_lbfluid_get_lattice() { return lblattice; }

ActiveLB lb_lbfluid_get_lattice_switch() { return lattice_switch; }

void lb_lbfluid_on_lb_params_change(LBParam field) {
  switch (field) {
  case LBParam::AGRID:
    if (lattice_switch == ActiveLB::CPU)
      lb_init();
#ifdef LB_GPU
    if (lattice_switch == ActiveLB::GPU && this_node == 0)
      lb_init_gpu();
#endif
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    LBBoundaries::lb_init_boundaries();
#endif
    break;
  case LBParam::DENSITY:
    lb_lbfluid_reinit_fluid();
    break;
  case LBParam::VISCOSITY:
  case LBParam::EXTFORCE:
  case LBParam::BULKVISC:
  case LBParam::KT:
    break;
  }
  lb_lbfluid_reinit_parameters();
}

Utils::Vector3d lb_lbfluid_calc_fluid_momentum() {
  Utils::Vector3d fluid_momentum{};
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lb_calc_fluid_momentum_GPU(fluid_momentum.data());
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    mpi_gather_stats(6, fluid_momentum.data(), nullptr, nullptr, nullptr);
  }
  return fluid_momentum;
}

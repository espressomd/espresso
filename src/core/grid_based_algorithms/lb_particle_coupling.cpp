#include <Random123/philox.h>
#include <boost/mpi.hpp>
#include <profiler/profiler.hpp>

#include "cells.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"
#include "lb_interpolation.hpp"
#include "lb_particle_coupling.hpp"
#include "lbgpu.hpp"
#include "random.hpp"

#include "utils/Counter.hpp"
#include "utils/u32_to_u64.hpp"
#include "utils/uniform.hpp"

LB_Particle_Coupling lb_particle_coupling;

void mpi_bcast_lb_particle_coupling_slave(int, int) {
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

void lb_lbcoupling_activate() {
  lb_particle_coupling.couple_to_md = true;
  mpi_bcast_lb_particle_coupling_slave(0, 0);
}

void lb_lbcoupling_deactivate() {
  if (lattice_switch != ActiveLB::NONE && this_node == 0 && n_part) {
    runtimeWarning("Recalculating forces, so the LB coupling forces are not "
                   "included in the particle force the first time step. This "
                   "only matters if it happens frequently during "
                   "sampling.\n");
  }

  lb_particle_coupling.couple_to_md = false;
  mpi_bcast_lb_particle_coupling_slave(0, 0);
}

void lb_lbcoupling_set_gamma(double gamma) {
  lb_particle_coupling.gamma = gamma;
  mpi_bcast_lb_particle_coupling();
}

double lb_lbcoupling_get_gamma() { return lb_particle_coupling.gamma; }

bool lb_lbcoupling_is_seed_required() {
  return not lb_particle_coupling.rng_counter_coupling.is_initialized();
}

uint64_t lb_coupling_get_rng_state_cpu() {
  return lb_particle_coupling.rng_counter_coupling->value();
}

uint64_t lb_lbcoupling_get_rng_state() {
  if (lattice_switch == ActiveLB::CPU) {
#ifdef LB
    return lb_coupling_get_rng_state_cpu();
#endif
  }
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    return lb_coupling_get_rng_state_gpu();
#endif
  }
  return {};
}

void lb_lbcoupling_set_rng_state(uint64_t counter) {
  if (lattice_switch == ActiveLB::CPU) {
#ifdef LB
    lb_particle_coupling.rng_counter_coupling =
        Utils::Counter<uint64_t>(counter);
    mpi_bcast_lb_particle_coupling();
#endif
  } else if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    lb_coupling_set_rng_state_gpu(counter);
#endif
  }
}

#if defined(LB) || defined(LB_GPU)

namespace {
/**
 * @brief Add a force to the lattice force density.
 * @param pos Position of the force
 * @param force Force in MD units.
 */
void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force) {
  /* transform momentum transfer to lattice units
     (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  auto const delta_j = -(time_step / lb_lbfluid_get_lattice_speed()) * force;
  lb_lbinterpolation_add_force_density(pos, delta_j);
}
} // namespace

/** Coupling of a single particle to viscous fluid with Stokesian friction.
 *
 *  Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
 *
 * @param[in,out] p         The coupled particle.
 * @param[in]     f_random  Additional force to be included.
 *
 * @return The viscous coupling force plus f_random.
 */
#ifdef LB
Utils::Vector3d lb_viscous_coupling(Particle *p,
                                    Utils::Vector3d const &f_random) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation
     (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  auto const interpolated_u =
      lb_lbinterpolation_get_interpolated_velocity(p->r.p) *
      lb_lbfluid_get_lattice_speed();

  Utils::Vector3d v_drift = interpolated_u;
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
  auto const force = -lb_lbcoupling_get_gamma() * (p->m.v - v_drift) + f_random;

  add_md_force(p->r.p, force);

  return force;
}
#endif

namespace {
bool in_local_domain(Utils::Vector3d const &pos) {
  auto const lblattice = lb_lbfluid_get_lattice();
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

    p.swim.v_source =
        lb_lbinterpolation_get_interpolated_velocity(source_position) *
        lb_lbfluid_get_lattice_speed();

    add_md_force(source_position, p.swim.f_swim * director);
  }
}
#endif
} // namespace

void lb_lbcoupling_calc_particle_lattice_ia(bool couple_virtual) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
    if (lb_particle_coupling.couple_to_md && this_node == 0) {
      switch (lb_lbinterpolation_get_interpolation_order()) {
      case (InterpolationOrder::linear):
        lb_calc_particle_lattice_ia_gpu<8>(couple_virtual,
                                           lb_lbcoupling_get_gamma());
        break;
      case (InterpolationOrder::quadratic):
        lb_calc_particle_lattice_ia_gpu<27>(couple_virtual,
                                            lb_lbcoupling_get_gamma());
        break;
      }
    }
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
#ifdef LB
    if (lb_particle_coupling.couple_to_md) {
      switch (lb_lbinterpolation_get_interpolation_order()) {
      case (InterpolationOrder::quadratic):
        throw std::runtime_error("The non-linear interpolation scheme is not "
                                 "implemented for the CPU LB.");
      case (InterpolationOrder::linear): {
#ifdef ENGINE
        ghost_communicator(&cell_structure.exchange_ghosts_comm,
                           GHOSTTRANS_SWIMMING);
#endif
        using rng_type = r123::Philox4x64;
        using ctr_type = rng_type::ctr_type;
        using key_type = rng_type::key_type;

        ctr_type c;
        if (lb_lbfluid_get_kT() > 0.0) {
          c = ctr_type{{lb_particle_coupling.rng_counter_coupling->value(),
                        static_cast<uint64_t>(RNGSalt::PARTICLES)}};
        } else {
          c = ctr_type{{0, 0}};
        }

        /* Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
         * The factor 12 comes from the fact that we use random numbers
         * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
         * time_step comes from the discretization.
         */
        auto const noise_amplitude = sqrt(12. * 2. * lb_lbcoupling_get_gamma() *
                                          lb_lbfluid_get_kT() / time_step);
        auto f_random = [&c](int id) -> Utils::Vector3d {
          if (lb_lbfluid_get_kT() > 0.0) {
            key_type k{{static_cast<uint32_t>(id)}};

            auto const noise = rng_type{}(c, k);

            using Utils::uniform;
            return Utils::Vector3d{uniform(noise[0]), uniform(noise[1]),
                                   uniform(noise[2])} -
                   Utils::Vector3d::broadcast(0.5);
          }
          return Utils::Vector3d{};
        };

        /* local cells */
        for (auto &p : local_cells.particles()) {
          if (!p.p.is_virtual or couple_virtual) {
            auto const force = lb_viscous_coupling(
                &p, noise_amplitude * f_random(p.identity()));
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
            if (!p.p.is_virtual || couple_virtual) {
              lb_viscous_coupling(&p, noise_amplitude * f_random(p.identity()));
#ifdef ENGINE
              add_swimmer_force(p);
#endif
            }
          }
        }
        break;
      }
      }
#endif
    }
  }
}

void lb_lbcoupling_propagate() {
  if (lb_lbfluid_get_kT() > 0.0) {
    if (lattice_switch == ActiveLB::CPU) {
#ifdef LB
      lb_particle_coupling.rng_counter_coupling->increment();
#endif
    } else if (lattice_switch == ActiveLB::GPU) {
#ifdef LB_GPU
      rng_counter_coupling_gpu->increment();
#endif
    }
  }
}
#endif

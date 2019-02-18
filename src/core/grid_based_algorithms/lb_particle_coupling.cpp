#include <Random123/philox.h>
#include <boost/mpi.hpp>
#include <profiler/profiler.hpp>

#include "cells.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "lattice.hpp"
#include "lb_interface.hpp"
#include "lb_particle_coupling.hpp"
#include "lbgpu.hpp"
#include "thermostat.hpp"
#include "utils/Counter.hpp"
#include "utils/u32_to_u64.hpp"
#include "utils/uniform.hpp"

LB_Particle_Coupling lb_particle_coupling;

void mpi_set_lb_coupling_counter_slave(int high, int low) {
  lb_particle_coupling.rng_counter_coupling =
      Utils::Counter<uint64_t>(Utils::u32_to_u64(static_cast<uint32_t>(high),
                                                 static_cast<uint32_t>(low)));
}

void mpi_bcast_lb_particle_coupling_slave(int, int) {
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

void lb_lbcoupling_activate() {
  lb_particle_coupling.couple_to_md = true;
  if (this_node == 0)
    mpi_bcast_lb_particle_coupling();
}

void lb_lbcoupling_deactivate() {
  lb_particle_coupling.couple_to_md = false;
  if (this_node == 0)
    mpi_bcast_lb_particle_coupling();
}

void lb_lbcoupling_set_friction(double friction) {
  lb_particle_coupling.friction = friction;
  mpi_bcast_lb_particle_coupling();
}

double lb_lbcoupling_get_friction() { return lb_particle_coupling.friction; }

uint64_t lb_coupling_get_rng_state_cpu() {
  return lb_particle_coupling.rng_counter_coupling.value();
}

uint64_t lb_lbcoupling_get_rng_state() {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    return lb_coupling_get_rng_state_cpu();
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    return lb_coupling_get_rng_state_gpu();
#endif
  }
  return {};
}

void lb_lbcoupling_set_rng_state(uint64_t counter) {
  if (lattice_switch & LATTICE_LB) {
#ifdef LB
    lb_particle_coupling.rng_counter_coupling =
        Utils::Counter<uint64_t>(counter);
    mpi_bcast_lb_particle_coupling();
#endif
  } else if (lattice_switch & LATTICE_LB_GPU) {
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
void add_md_force(Vector3d const &pos, Vector3d const &force) {
  /* transform momentum transfer to lattice units
     (Eq. (12) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */

  const auto agrid = lb_lbfluid_get_agrid();
  auto const delta_j = -(time_step * lb_lbfluid_get_tau() / agrid) * force;
  lb_lbfluid_add_force_density(pos, delta_j);
}
} // namespace

/** Coupling of a single particle to viscous fluid with Stokesian friction.
 *
 *  Section II.C. Ahlrichs and Duenweg, JCP 111(17):8225 (1999)
 *
 * @param p          The coupled particle (Input).
 * @param f_random   Additional force to be included.
 *
 * @return The viscous coupling force plus f_random.
 */
#ifdef LB
Vector3d lb_viscous_coupling(Particle *p, Vector3d const &f_random) {
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
  auto const force =
      -lb_lbcoupling_get_friction() * (p->m.v - v_drift) + f_random;

  add_md_force(p->r.p, force);

  return force;
}
#endif

namespace {
bool in_local_domain(Vector3d const &pos) {
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

    p.swim.v_source = lb_lbfluid_get_interpolated_velocity(source_position);

    add_md_force(source_position, p.swim.f_swim * director);
  }
}
#endif
} // namespace

void lb_lbcoupling_calc_particle_lattice_ia(bool couple_virtual) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  if (lattice_switch & LATTICE_LB_GPU) {
#ifdef LB_GPU
    if (lb_particle_coupling.couple_to_md && this_node == 0)
      lb_calc_particle_lattice_ia_gpu(couple_virtual,
                                      lb_lbcoupling_get_friction());
#endif
  } else if (lattice_switch & LATTICE_LB) {
#ifdef LB
    if (lb_particle_coupling.couple_to_md) {
      using rng_type = r123::Philox4x64;
      using ctr_type = rng_type::ctr_type;
      using key_type = rng_type::key_type;

      ctr_type c{{lb_particle_coupling.rng_counter_coupling.value(),
                  static_cast<uint64_t>(RNGSalt::PARTICLES)}};

      /* Eq. (16) Ahlrichs and Duenweg, JCP 111(17):8225 (1999).
       * The factor 12 comes from the fact that we use random numbers
       * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
       * time_step comes from the discretization.
       */
      auto const noise_amplitude =
          sqrt(12. * 2. * lb_lbcoupling_get_friction() * lb_lbfluid_get_kT() /
               time_step);
      auto f_random = [&c](int id) -> Vector3d {
        if (lb_lbfluid_get_kT() > 0.0) {
          key_type k{{static_cast<uint32_t>(id)}};

          auto const noise = rng_type{}(c, k);

          using Utils::uniform;
          return Vector3d{uniform(noise[0]), uniform(noise[1]),
                          uniform(noise[2])} -
                 Vector3d::broadcast(0.5);
        } else {
          return Vector3d{};
        }
      };

      /* local cells */
      for (auto &p : local_cells.particles()) {
        if (!p.p.is_virtual or thermo_virtual or
            (p.p.is_virtual && couple_virtual)) {
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
#endif
  }
}

void lb_lbcoupling_propagate() {
  lb_particle_coupling.rng_counter_coupling.increment();
}
#endif

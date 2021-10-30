/*
 * Copyright (C) 2010-2020 The ESPResSo project
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

#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS
#include "sd_interface.hpp"

#include "stokesian_dynamics/sd_cpu.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/mpi/scatter_buffer.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
/* type for particle data transfer between nodes */
struct SD_particle_data {
  SD_particle_data() = default;
  explicit SD_particle_data(Particle const &p)
      : type(p.p.type), pos(p.r.p), ext_force(p.f) {}

  int type = 0;

  /* particle position */
  Utils::Vector3d pos = {0., 0., 0.};

  /* external force */
  ParticleForce ext_force;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &type;
    ar &pos;
    ar &ext_force;
  }
};

double sd_viscosity = -1.0;

std::unordered_map<int, double> radius_dict;

double sd_kT = 0.0;

int sd_flags = 0;

/** Buffer that holds the (translational and angular) velocities of the local
 *  particles on each node, used for returning results. */
std::vector<double> v_sd{};

} // namespace

BOOST_IS_BITWISE_SERIALIZABLE(SD_particle_data)

/** Update translational and rotational velocities of all particles. */
void sd_update_locally(ParticleRange const &parts) {
  std::size_t i = 0;

  // Even though on the head node, the v_sd vector is larger than
  // the (local) parts vector, this should still work. Because the local
  // particles correspond to the first 6*n entries in the head node's v_sd
  // (which holds the velocities of ALL particles).

  for (auto &p : parts) {
    // skip virtual particles
    if (p.p.is_virtual) {
      continue;
    }

    // Copy velocities
    p.m.v[0] = v_sd[6 * i + 0];
    p.m.v[1] = v_sd[6 * i + 1];
    p.m.v[2] = v_sd[6 * i + 2];

    p.m.omega[0] = v_sd[6 * i + 3];
    p.m.omega[1] = v_sd[6 * i + 4];
    p.m.omega[2] = v_sd[6 * i + 5];

    i++;
  }
}

void set_sd_viscosity(double eta) {
  if (eta < 0.0) {
    throw std::runtime_error("Viscosity has an invalid value: " +
                             std::to_string(eta));
  }

  sd_viscosity = eta;
}

void set_sd_radius_dict(std::unordered_map<int, double> const &x) {
  /* Check that radii are positive */
  for (auto const &kv : x) {
    if (kv.second < 0.) {
      throw std::runtime_error(
          "Particle radius for type " + std::to_string(kv.first) +
          " has an invalid value: " + std::to_string(kv.second));
    }
  }

  radius_dict = x;
}

void set_sd_kT(double kT) {
  if (kT < 0.0) {
    throw std::runtime_error("kT has an invalid value: " + std::to_string(kT));
  }

  sd_kT = kT;
}

double get_sd_kT() { return sd_kT; }

void set_sd_flags(int flg) { sd_flags = flg; }

void propagate_vel_pos_sd(const ParticleRange &particles,
                          const boost::mpi::communicator &comm,
                          const double time_step) {
  static std::vector<SD_particle_data> parts_buffer{};

  parts_buffer.clear();
  boost::transform(particles, std::back_inserter(parts_buffer),
                   [](auto const &p) { return SD_particle_data(p); });
  Utils::Mpi::gather_buffer(parts_buffer, comm, 0);

  /* Buffer that holds local particle data, and all particles on the head
   * node used for sending particle data to head node. */
  if (comm.rank() == 0) {
    std::size_t n_part = parts_buffer.size();

    static std::vector<double> x_host{};
    static std::vector<double> f_host{};
    static std::vector<double> a_host{};

    x_host.resize(6 * n_part);
    f_host.resize(6 * n_part);
    a_host.resize(n_part);

    std::size_t i = 0;
    for (auto const &p : parts_buffer) {
      x_host[6 * i + 0] = p.pos[0];
      x_host[6 * i + 1] = p.pos[1];
      x_host[6 * i + 2] = p.pos[2];
      // Actual orientation is not needed, just need default.
      x_host[6 * i + 3] = 1;
      x_host[6 * i + 4] = 0;
      x_host[6 * i + 5] = 0;

      f_host[6 * i + 0] = p.ext_force.f[0];
      f_host[6 * i + 1] = p.ext_force.f[1];
      f_host[6 * i + 2] = p.ext_force.f[2];

      f_host[6 * i + 3] = p.ext_force.torque[0];
      f_host[6 * i + 4] = p.ext_force.torque[1];
      f_host[6 * i + 5] = p.ext_force.torque[2];

      double radius = radius_dict[p.type];

      a_host[i] = radius;

      ++i;
    }

    v_sd = sd_cpu(x_host, f_host, a_host, n_part, sd_viscosity,
                  std::sqrt(sd_kT / time_step),
                  static_cast<std::size_t>(stokesian.rng_counter()),
                  static_cast<std::size_t>(stokesian.rng_seed()), sd_flags);
  } else { // if (this_node == 0)
    v_sd.resize(particles.size() * 6);
  } // if (this_node == 0) {...} else

  Utils::Mpi::scatter_buffer(v_sd.data(),
                             static_cast<int>(particles.size() * 6), comm, 0);
  sd_update_locally(particles);
}

#endif // STOKESIAN_DYNAMICS

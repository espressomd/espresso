#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS

#include <iostream>
#include <string>
#include <vector>

#include <boost/config.hpp>

#include "errorhandling.hpp"
#include "integrate.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"

#if defined(BLAS) && defined(LAPACK)
#include "stokesian_dynamics/sd_cpu.hpp"
#endif

#ifdef CUDA
#include "stokesian_dynamics/sd_gpu.hpp"
#endif

#include "stokesian_dynamics/sd_interface.hpp"

namespace {
double sd_viscosity = -1.0;

enum { CPU, GPU, INVALID } device = INVALID;
} // namespace

void set_sd_viscosity(double eta) { sd_viscosity = eta; }

double get_sd_viscosity() { return sd_viscosity; }

void set_sd_device(std::string const &dev) {
  if (dev == "cpu") {
    device = CPU;
  } else if (dev == "gpu") {
    device = GPU;
  } else {
    device = INVALID;
  }
}

std::string get_sd_device() {
  switch (device) {
  case CPU:
    return "cpu";
  case GPU:
    return "gpu";
  default:
    return "invalid";
  }
}

void propagate_vel_pos_sd() {
  if (thermo_switch & THERMO_SD) {
    if (BOOST_UNLIKELY(sd_viscosity < 0.0)) {
      runtimeErrorMsg() << "sd_viscosity has an invalid value: " +
                               std::to_string(sd_viscosity);
      return;
    }

    auto &parts = partCfg();
    std::size_t n_part = parts.size();
    std::vector<int> id(n_part);
    std::vector<double> x_host(6 * n_part);
    std::vector<double> f_host(6 * n_part);
    std::vector<double> a_host(n_part);

    std::size_t i = 0;
    for (auto const &p : parts) {
      id[i] = p.p.is_virtual ? -1 : p.p.identity;
      x_host[6 * i + 0] = p.r.p[0];
      x_host[6 * i + 1] = p.r.p[1];
      x_host[6 * i + 2] = p.r.p[2];
      x_host[6 * i + 3] = 0;
      x_host[6 * i + 4] = 0;
      x_host[6 * i + 5] = 0;

      f_host[6 * i + 0] = p.f.f[0];
      f_host[6 * i + 1] = p.f.f[1];
      f_host[6 * i + 2] = p.f.f[2];
      f_host[6 * i + 3] = p.f.torque[0];
      f_host[6 * i + 4] = p.f.torque[1];
      f_host[6 * i + 5] = p.f.torque[2];

      if (BOOST_UNLIKELY(p.p.radius < 0.0)) {
        runtimeErrorMsg() << "particle " + std::to_string(p.p.identity) +
                                 " has an invalid radius: " +
                                 std::to_string(p.p.radius);
        return;
      }
      a_host[i] = p.p.radius;

      ++i;
    }

    std::vector<double> v_sd;
    switch (device) {

#if defined(BLAS) && defined(LAPACK)
    case CPU:
      v_sd = sd_cpu(x_host, f_host, a_host, n_part, sd_viscosity);
      break;
#endif

#ifdef CUDA
    case GPU:
      v_sd = sd_gpu(x_host, f_host, a_host, n_part, sd_viscosity);
      break;
#endif

    default:
      runtimeErrorMsg()
          << "Invalid device for Stokesian dynamics. Available devices:"
#if defined(BLAS) && defined(LAPACK)
             " cpu"
#endif
#ifdef CUDA
             " gpu"
#endif
          ;
      return;
    }

    for (int i = 0; i < n_part; ++i) {
      // skip virtual particles
      if (id[i] == -1) {
        continue;
      }

      // set new velocities and angular velocities
      double v[] = {v_sd[6 * i + 0], v_sd[6 * i + 1], v_sd[6 * i + 2]};
      set_particle_v(id[i], v);

      Vector3d omega = {v_sd[6 * i + 3], v_sd[6 * i + 4], v_sd[6 * i + 5]};
      set_particle_omega_body(id[i], omega);

      // integrate and set new positions
      double p[] = {x_host[6 * i + 0] + time_step * v[0],
                    x_host[6 * i + 1] + time_step * v[1],
                    x_host[6 * i + 2] + time_step * v[2]};
      place_particle(id[i], p);
    }
  }
}

#endif // STOKESIAN_DYNAMICS

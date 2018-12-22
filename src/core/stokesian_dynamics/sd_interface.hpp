#ifndef STOKESIAN_DYNAMICS_INTERFACE_H
#define STOKESIAN_DYNAMICS_INTERFACE_H

#include <vector>

#include "integrate.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"

#include "stokesian_dynamics/sd_cpu.hpp"
#include "stokesian_dynamics/sd_gpu.hpp"

inline void propagate_vel_pos_sd() {
  if (thermo_switch & THERMO_SD) {
    auto &parts = partCfg();
    std::size_t n_part = parts.size();
    std::vector<double> x_host(6 * n_part);
    std::vector<double> f_host(6 * n_part);

    std::size_t i = 0;
    for (auto const &p : parts) {
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

      ++i;
    }

    std::vector<double> v_sd = sd_cpu(x_host, f_host, n_part);

    /*
    i = 0;
    for (auto &&p : parts) {
      // v_{n+1} = v_n + f(x_n,v_n) dt
      p.m.v[0] += p.f.f[0] * time_step + v_sd[6 * i + 0];
      p.m.v[1] += p.f.f[1] * time_step + v_sd[6 * i + 1];
      p.m.v[2] += p.f.f[2] * time_step + v_sd[6 * i + 2];

      // x_{n+1} = x_n + v_{n+1} dt
      p.r.p[0] += p.m.v[0] * time_step;
      p.r.p[1] += p.m.v[1] * time_step;
      p.r.p[2] += p.m.v[2] * time_step;

      ++i;
    }
    */
  }
}

#endif

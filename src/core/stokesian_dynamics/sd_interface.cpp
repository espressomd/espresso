#include <iostream>
#include <vector>

#include "integrate.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"

#include "stokesian_dynamics/sd_interface.hpp"
#include "stokesian_dynamics/sd_cpu.hpp"
#include "stokesian_dynamics/sd_gpu.hpp"

void propagate_vel_pos_sd() {
  if (thermo_switch & THERMO_SD) {
    auto &parts = partCfg();
    std::size_t n_part = parts.size();
    std::vector<int> id(n_part);
    std::vector<double> x_host(6 * n_part);
    std::vector<double> f_host(6 * n_part);

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

      ++i;
    }

    std::vector<double> v_sd = sd_cpu(x_host, f_host, n_part);

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

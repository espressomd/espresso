#include <iostream>
#include <vector>

#include "sd.hpp"
#include "sd_gpu.hpp"

std::vector<double> sd_gpu(std::vector<double> const &x_host,
                           std::vector<double> const &f_host,
                           std::size_t n_part) {
  double const eta = 1.0;
  sd::solver<policy::device, double> viscous_force{eta, n_part};
  return viscous_force.calc_vel(x_host, f_host);
}

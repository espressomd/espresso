#include <iostream>
#include <vector>

#include "sd.hpp"
#include "sd_cpu.hpp"

std::vector<double> sd_cpu(std::vector<double> x_host,
                           std::vector<double> f_host,
                           std::size_t n_part) {
  double const eta = 1.0;
  sd::solver<policy::host, double> viscous_force{eta, n_part};
  return viscous_force.calc_vel(x_host, f_host);
}

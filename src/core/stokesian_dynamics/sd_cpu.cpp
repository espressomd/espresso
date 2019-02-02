#include <iostream>
#include <vector>

#include "sd.hpp"
#include "sd_cpu.hpp"

std::vector<double> sd_cpu(std::vector<double> const &x_host,
                           std::vector<double> const &f_host,
                           std::vector<double> const &a_host,
                           std::size_t n_part, double eta) {
  sd::solver<policy::host, double> viscous_force{eta, n_part};
  return viscous_force.calc_vel(x_host, f_host, a_host);
}

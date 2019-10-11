#pragma once

#include <vector>

std::vector<double> sd_cpu(std::vector<double> const &x_host,
                           std::vector<double> const &f_host,
                           std::vector<double> const &a_host,
                           std::size_t n_part, double eta, double kT,
                           std::size_t offset, std::size_t seed, int flg);

#pragma once

#include <vector>

std::vector<double> sd_gpu(std::vector<double> const &x_host,
                           std::vector<double> const &f_host,
                           std::size_t n_part);

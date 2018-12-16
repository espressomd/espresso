#pragma once

#include <vector>

std::vector<double> sd_cpu(std::vector<double> x_host,
                           std::vector<double> f_host,
                           std::size_t n_part);

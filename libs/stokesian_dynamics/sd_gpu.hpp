#ifndef SD_GPU_HPP
#define SD_GPU_HPP

#include <vector>

std::vector<double> sd_gpu(std::vector<double> const &x_host,
                           std::vector<double> const &f_host,
                           std::vector<double> const &a_host,
                           std::size_t n_part, double eta, double sqrt_kT_Dt,
                           std::size_t offset, std::size_t seed, int flg);

#endif

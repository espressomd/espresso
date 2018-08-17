#ifndef P3M_ERROR_GPU_HPP
#define P3M_ERROR_GPU_HPP

#include "config.hpp"

#ifdef CUDA
double p3m_k_space_error_gpu(double prefactor, int *mesh, int cao, int npart, double sum_q2, double alpha_L, double *box);
#endif
#endif

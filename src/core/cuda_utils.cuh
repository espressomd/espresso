/*
 * Copyright (C) 2013-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CORE_CUDA_UTILS_CUH
#define CORE_CUDA_UTILS_CUH

#if !defined(__CUDACC__)
#error Do not include CUDA headers in normal C++-code!!!
#endif

#include "cuda_utils.hpp"

#include <cuda.h>

#include <string>

class cuda_runtime_error_cuda : public cuda_runtime_error {
public:
  cuda_runtime_error_cuda(cudaError_t error)
      : cuda_runtime_error(error_message(error)) {}

private:
  std::string error_message(cudaError_t error) {
    const char *cuda_error = cudaGetErrorString(error);
    return std::string("CUDA error: ") + cuda_error;
  }
};

/** Convert CUDA error codes into runtime errors. */
#define CUDA_CHECK(statement)                                                  \
  {                                                                            \
    cudaError_t const error_code = (statement);                                \
    if (error_code != cudaSuccess) {                                           \
      throw cuda_runtime_error_cuda(error_code);                               \
    }                                                                          \
  }

/** CUDA streams for parallel computing on CPU and GPU */
extern cudaStream_t stream[1];

/** In case of error during CUDA memory allocation and memory copy, print
 *  the error message and exit.
 *  @param CU_err cuda error code
 *  @param file  .cu file were the error took place
 *  @param line  line of the file were the error took place
 */
void cuda_safe_mem_exit(cudaError_t CU_err, const char *file,
                        unsigned int line);

/** In case of error during a CUDA operation, print the error message and exit.
 */
void cuda_check_errors_exit(const dim3 &block, const dim3 &grid,
                            const char *function, const char *file,
                            unsigned int line);

#define cuda_safe_mem(a) cuda_safe_mem_exit((a), __FILE__, __LINE__)

#define KERNELCALL_shared(_function, _grid, _block, _stream, ...)              \
  _function<<<_grid, _block, _stream, stream[0]>>>(__VA_ARGS__);               \
  cuda_check_errors_exit(_grid, _block, #_function, __FILE__, __LINE__);

#define KERNELCALL_stream(_function, _grid, _block, _stream, ...)              \
  _function<<<_grid, _block, 0, _stream>>>(__VA_ARGS__);                       \
  cuda_check_errors_exit(_grid, _block, #_function, __FILE__, __LINE__);

#define KERNELCALL(_function, _grid, _block, ...)                              \
  KERNELCALL_shared(_function, _grid, _block, 0, ##__VA_ARGS__)

#endif

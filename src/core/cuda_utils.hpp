/*
  Copyright (C) 2013-2018 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _CUDA_UTILS_HPP
#define _CUDA_UTILS_HPP

#if !defined(__HIPCC__) and !defined(__CUDACC__)
#error Do not include CUDA headers in normal C++-code!!!
#endif

/**cuda streams for parallel computing on cpu and gpu */
extern cudaStream_t stream[1];

/**erroroutput for memory allocation and memory copy
 * @param err cuda error code
 * @param *file .cu file were the error took place
 * @param line line of the file were the error took place
 */

void _cuda_safe_mem(cudaError_t err, const char *file, unsigned int line);

void _cuda_check_errors(const dim3 &block, const dim3 &grid,
                        const char *function, const char *file,
                        unsigned int line);

#define cuda_safe_mem(a) _cuda_safe_mem((a), __FILE__, __LINE__)

#define KERNELCALL_shared(_f, _a, _b, _s, ...)                                 \
  hipLaunchKernelGGL(_f, _a, _b, _s, stream[0], ##__VA_ARGS__);                \
  _cuda_check_errors(_a, _b, #_f, __FILE__, __LINE__);

#define KERNELCALL_stream(_function, _grid, _block, _stream, ...)              \
  hipLaunchKernelGGL(_function, _grid, _block, 0, _stream, ##__VA_ARGS__);     \
  _cuda_check_errors(_grid, _block, #_function, __FILE__, __LINE__);

#define KERNELCALL(_f, _a, _b, ...)                                            \
  KERNELCALL_shared(_f, _a, _b, 0, ##__VA_ARGS__)

#endif

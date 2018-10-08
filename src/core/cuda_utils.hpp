#include "hip/hip_runtime.h"
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

#ifndef __HIPCC__
#error Do not include CUDA headers in normal C++-code!!!
#endif

/**cuda streams for parallel computing on cpu and gpu */
extern hipStream_t stream[1];

extern hipError_t CU_err;
extern hipError_t _err;

/**erroroutput for memory allocation and memory copy
 * @param err cuda error code
 * @param *file .cu file were the error took place
 * @param line line of the file were the error took place
 */

void _cuda_safe_mem(hipError_t err, const char *file, unsigned int line);

void _cuda_check_errors(const dim3 &block, const dim3 &grid,
                        const char *function, const char *file,
                        unsigned int line);

#define cuda_safe_mem(a) _cuda_safe_mem((a), __FILE__, __LINE__)

#define KERNELCALL_shared(_f, _a, _b, _s, _params)                             \
  _f<<<_a, _b, _s, stream[0]>>> _params;                                       \
  _cuda_check_errors(_a, _b, #_f, __FILE__, __LINE__);

#define KERNELCALL_stream(_function, _grid, _block, _stream, _params)          \
  _function<<<_grid, _block, 0, _stream>>> _params;                            \
  _cuda_check_errors(_grid, _block, #_function, __FILE__, __LINE__);

#define KERNELCALL(_f, _a, _b, _params)                                        \
  KERNELCALL_shared(_f, _a, _b, 0, _params)

#endif

/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "errorhandling.hpp"

#include "utils.cuh"

#include <cstdio>

cudaStream_t stream[1];

void cuda_check_errors_exit(const dim3 &block, const dim3 &grid,
                            const char *function, const char *file,
                            unsigned int line) {
  cudaError_t CU_err = cudaGetLastError();
  if (CU_err != cudaSuccess) {
    fprintf(stderr,
            "error \"%s\" calling %s with dim %d %d %d, grid %d %d "
            "%d in %s:%u\n",
            cudaGetErrorString(CU_err), function, block.x, block.y, block.z,
            grid.x, grid.y, grid.z, file, line);
    errexit();
  }
}

void cuda_safe_mem_exit(cudaError_t CU_err, const char *file,
                        unsigned int line) {
  if (CU_err != cudaSuccess) {
    fprintf(stderr, "CUDA Memory error at %s:%u.\n", file, line);
    fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(CU_err));
    if (CU_err == cudaErrorInvalidValue)
      fprintf(stderr, "You may have tried to allocate zero memory at %s:%u.\n",
              file, line);
    errexit();
  } else {
    CU_err = cudaGetLastError();
    if (CU_err != cudaSuccess) {
      fprintf(stderr,
              "Error found during memory operation. Possibly however "
              "from a failed operation before. %s:%u.\n",
              file, line);
      printf("CUDA error: %s\n", cudaGetErrorString(CU_err));
      if (CU_err == cudaErrorInvalidValue)
        fprintf(stderr,
                "You may have tried to allocate zero memory before %s:%u.\n",
                file, line);
      errexit();
    }
  }
}

/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "init.hpp"

#include "config/config.hpp"

#include <cassert>

#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ESPRESSO_FFTW
void fft_on_program_start() {
#ifdef _OPENMP
  int omp_num_threads = 1;
#pragma omp parallel
  {
#pragma omp single
    omp_num_threads = omp_get_num_threads();
  }
  [[maybe_unused]] auto const init_status_success = fftw_init_threads();
  assert(init_status_success);
  fftw_plan_with_nthreads(omp_num_threads);
#endif // _OPENMP
}
#endif // ESPRESSO_FFTW

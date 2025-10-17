/*
 * Copyright (C) 2010-2023 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#pragma once

#ifdef ESPRESSO_LIKWID
#include <likwid.h>
void inline profiling_init() {
#ifdef ESPRESSO_LIKWID
  likwid_markerInit();
#pragma omp parallel
  {
    likwid_markerThreadInit();
  }
#endif
}

void inline profiling_section_begin(const char *name) {
#pragma omp parallel
  {
    likwid_markerStartRegion(name);
#pragma omp barrier
  }
}

void inline profiling_section_end(const char *name) {
#pragma omp parallel
  {
    likwid_markerStopRegion(name);
#pragma omp barrier
  }
}

#define PROFILING_INIT profiling_init()
#define PROFILING_CLOSE likwid_markerClose()
#define PROFILING_SECTION_BEGIN(name) profiling_section_begin(name)
#define PROFILING_SECTION_END(name) profiling_section_end(name)
#define PROFILING_MARK_FUNCTION
#define PROFILING_MARK_LOOP_BEGIN(loop_id, name)
#define PROFILING_MARK_LOOP_ITERATION(loop_id, step)
#define PROFILING_MARK_LOOP_END(loop_id)
#elif defined(ESPRESSO_CALIPER)
#include <caliper/cali.h>
#define PROFILING_INIT
#define PROFILING_CLOSE
#define PROFILING_SECTION_BEGIN(name) CALI_MARK_BEGIN(name)
#define PROFILING_SECTION_END(name) CALI_MARK_END(name)
#define PROFILING_MARK_FUNCTION CALI_CXX_MARK_FUNCTION
#define PROFILING_MARK_LOOP_BEGIN(loop_id, name)                               \
  CALI_CXX_MARK_LOOP_BEGIN(loop_id, name)
#define PROFILING_MARK_LOOP_ITERATION(loop_id, step)                           \
  CALI_CXX_MARK_LOOP_ITERATION(loop_id, step)
#define PROFILING_MARK_LOOP_END(loop_id) CALI_CXX_MARK_LOOP_END(loop_id)
#else
#define PROFILING_INIT
#define PROFILING_CLOSE
#define PROFILING_SECTION_BEGIN(name)
#define PROFILING_SECTION_END(name)
#define PROFILING_MARK_FUNCTION
#define PROFILING_MARK_LOOP_BEGIN(loop_id, name)
#define PROFILING_MARK_LOOP_ITERATION(loop_id, step)
#define PROFILING_MARK_LOOP_END(loop_id)
#endif
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

#ifdef ESPRESSO_CALIPER
#include <caliper/cali.h>
#endif

#ifdef ESPRESSO_LIKWID
#include "likwid.h"
#endif

void inline profiling_init(){
#ifdef ESPRESSO_LIKWID
  likwid_markerInit();
#pragma omp parallel
    {
        likwid_markerThreadInit();
    }
#endif
}
void inline profiling_close(){
#ifdef ESPRESSO_LIKWID
  likwid_markerClose();
#endif
}

void inline profiling_section_begin(const char* name){
#ifdef ESPRESSO_LIKWID
#pragma omp parallel
    {
        likwid_markerStartRegion(name);
    }
#endif
#ifdef ESPRESSO_CALIPER
    CALI_MARK_BEGIN(name);
#endif
}

void inline profiling_section_end(const char* name){
#ifdef ESPRESSO_CALIPER
    CALI_MARK_END(name);
#endif
#ifdef ESPRESSO_LIKWID
#pragma omp parallel
    {
        likwid_markerStopRegion(name);
    }
#endif
}
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

#pragma once

// guard for hdf5.h
#if not defined(_H5public_H)
#ifdef OMPI_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX
#endif
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif
#endif // not defined(_H5public_H)

#if not defined(H5_USE_BOOST)
#define H5_USE_BOOST
#endif

#include <H5public.h>

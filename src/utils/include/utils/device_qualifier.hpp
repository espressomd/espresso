/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef UTILS_DEVICE_QUALIFIER_HPP
#define UTILS_DEVICE_QUALIFIER_HPP

#if defined(__CUDACC__)
#define DEVICE_COMPILER cudacc
#elif defined(__HIPCC__)
#define DEVICE_COMPILER hipcc
#endif

#if defined(DEVICE_COMPILER)
#define DEVICE_THROW(E)
#define DEVICE_QUALIFIER __host__ __device__
#define DEVICE_ASSERT(A) void((A))
#else
#define DEVICE_THROW(E) throw(E)
#define DEVICE_QUALIFIER
#define DEVICE_ASSERT(A) assert((A))
#endif

#endif // ESPRESSO_DEVICE_QUALIFIER_HPP

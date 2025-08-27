/*
 * Copyright (C) 2010-2025 The ESPResSo project
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
#ifndef ESPRESSO_CONFIG_HPP
#define ESPRESSO_CONFIG_HPP

/* Prevent C++ bindings in MPI (there is a DataType called LB in there) */
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#include "config/config-features.hpp"

/**
 * @brief Precision below which a double-precision float is assumed to be zero.
 * Used in comparisons to determine if two floating-point numbers are equal.
 */
inline constexpr auto round_error_prec = 1e-14;

/**
 * @brief Special cutoff value for an inactive interaction.
 * Non-bonded potentials that have this cutoff are never evaluated.
 */
inline constexpr double inactive_cutoff = -1.;

/**
 * @brief Special cutoff value for an inactive bond.
 * Bonds that have this cutoff are never evaluated.
 */
inline constexpr double bonded_inactive_cutoff = -1.;

#endif // ESPRESSO_CONFIG_HPP

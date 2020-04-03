/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef STATISTICS_CHAIN_H
#define STATISTICS_CHAIN_H
/** \file
 *
 *  This file contains the code for statistics on chains.
 */

#include <array>

/** \name Exported Functions */
/************************************************************/
/*@{*/

/**
 * @brief Calculate the end-to-end-distance.
 *
 * Calculates the average end-to-end-distance of a range
 * of monodisperse polymers with continuous ids.
 *
 * @param chain_start The id of the first monomer of the first chain.
 * @param chain_n_chains Number of chains contained in the range.
 * @param chain_length The length of every chain.
 */
std::array<double, 4> calc_re(int chain_start, int chain_n_chains,
                              int chain_length);

/**
 * @brief Calculate the radius of gyration.
 *
 * Calculates the average radius of gyration of a range
 * of monodisperse polymers with continuous ids.
 *
 * @param chain_start The id of the first monomer of the first chain.
 * @param chain_n_chains Number of chains contained in the range.
 * @param chain_length The length of every chain.
 */
std::array<double, 4> calc_rg(int chain_start, int chain_n_chains,
                              int chain_length);

/**
 * @brief Calculate the hydrodynamic radius (ref. Kirkwood-Zimm theory).
 *
 * Calculates the average hydrodynamic radius of a range
 * of monodisperse polymers with continuous ids.
 *
 * @param chain_start The id of the first monomer of the first chain.
 * @param chain_n_chains Number of chains contained in the range.
 * @param chain_length The length of every chain.
 */
std::array<double, 2> calc_rh(int chain_start, int chain_n_chains,
                              int chain_length);
/*@}*/

#endif

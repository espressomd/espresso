/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
/** \file
 *  This file contains the asynchronous MPI communication for interactions.
 */
#ifndef _INTERACTIONS_HPP
#define _INTERACTIONS_HPP

/** Calculate the maximal cutoff of all interactions. */
double maximal_cutoff();

/** Check electrostatic and magnetostatic methods are properly initialized.
 *  @return true if sanity checks failed.
 */
bool long_range_interactions_sanity_checks();

/** Send new IA params.
 *  Also calls \ref on_short_range_ia_change.
 *
 *  Used for both bonded and non-bonded interaction parameters. Therefore
 *  @p i and @p j are used depending on their value:
 *
 *  \param i   particle type for non-bonded interaction parameters /
 *             bonded interaction type number.
 *  \param j   if not negative: particle type for non-bonded interaction
 *             parameters / if negative: flag for bonded interaction
 */
void mpi_bcast_ia_params(int i, int j);

#endif

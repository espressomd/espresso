/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
#ifndef DIAMOND_H
#define DIAMOND_H
/** \file

    This file contains everything needed to create a start-up
    diamond structure-like configuration of (partially charged)
    polymer chains with counterions and salt molecules.

*/

#include "PartCfg.hpp"
#include "particle_data.hpp"
#include "utils/Vector.hpp"

/*************************************************************
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/

/** C implementation of 'counterions \<N_CI\> [options]'.
 *  @param  N_CI        = number of counterions to create
 *  @param  part_id     = particle number of the first counterion (defaults to
 *  'n_total_particles')
 *  @param  mode        = selects setup mode: Self avoiding walk (SAW) or plain
 *  random walk (RW) (defaults to 'SAW')
 *  @param  shield      = shield around each particle another particle's
 *  position may not enter if using SAW (defaults to '0.0')
 *  @param  max_try     = how often a monomer should be reset if current
 *  position collides with a previous particle (defaults to '30000')
 *  @param  val_CI      = valency of the counterions (defaults to '-1.0')
 *  @param  type_CI     = type number of the counterions to be used with "part"
 *  (default to '2')
 *  @return Returns how often the attempt to place a particle failed in the
 *  worst case.
 */
int create_counterions(PartCfg &, int N_CI, int part_id, int mode,
                       double shield, int max_try, double val_CI, int type_CI);

/** C implementation of 'diamond \<a\> \<bond_length\> \<MPC\> [options]' */
int create_diamond(PartCfg &, double a, double bond_length, int MPC, int N_CI,
                   double val_nodes, double val_cM, double val_CI, int cM_dist,
                   int nonet);

#endif

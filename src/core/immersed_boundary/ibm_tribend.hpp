/*
Copyright (C) 2010-2018 The ESPResSo project

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

#ifndef IBM_TRIBEND_H
#define IBM_TRIBEND_H

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "bonded_interactions/bonded_interaction_data.hpp"

// DEBUG stuff
extern double maxBendingForce, maxBendingDist, maxX;

// This function is used to set the parameters
// Also calculates and stores the reference state
int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2,
                          const int ind3, const int ind4, const double kb,
                          const bool flat);
// For reading checkpoints.
// Idea: * parameters are set in the run-continue script
//       * also reference shape is recomputed there
//       * only pass kB value here to check consistency
int IBM_Tribend_ResetParams(const int bond_type, const double kb);

// This function calculates and adds the actual force
void IBM_Tribend_CalcForce(Particle *p1, const int numPartners,
                           Particle **const partners,
                           const Bonded_ia_parameters &iaparams);

#endif

#endif

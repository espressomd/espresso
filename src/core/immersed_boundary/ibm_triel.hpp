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

#ifndef IBM_TRIEL_H
#define IBM_TRIEL_H

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "bonded_interactions/bonded_interaction_data.hpp"

// This function is used to set the parameters
// Also calculates and stores the reference state
int IBM_Triel_SetParams(const int bond_type, const int ind1, const int ind2,
                        const int ind3, const double maxDist,
                        const tElasticLaw elasticLaw, const double k1,
                        const double k2);
// For reading checkpoints.
// Idea: * parameters are set in the run-continue script
//       * also reference shape is recomputed there
//       * only pass two values here to check consistency
int IBM_Triel_ResetParams(const int bond_type, const double k1,
                          const double l0);

// This function calculates and adds the actual force
int IBM_Triel_CalcForce(Particle *p1, Particle *p2, Particle *p3,
                        Bonded_ia_parameters *iaparams);

#endif

#endif

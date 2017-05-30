/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file comfixed.cpp
 *
 *  Implementation of \ref comfixed.hpp
 */
#include "cells.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

#ifdef COMFIXED

int comfixed_set_params(int part_type_a, int part_type_b, int flag) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data)
    return 1;

  if (n_nodes > 1)
    return 2;

  data->COMFIXED_flag = flag;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  for (auto &p : local_cells.particles()) {
    if (p.p.type == part_type_a) {
      for (int j = 0; j < 3; j++) {
        p.m.v[j] = 0.;
        p.f.f[j] = 0.;
      }
    }
  }

  return ES_OK;
}

void calc_comfixed() {
  IA_parameters *ia_params;
  int t0;
  int j;
  double fsum0[3], type_mass;

  for (t0 = 0; t0 < n_particle_types; t0++) {
    ia_params = get_ia_param(t0, t0);
    if (ia_params->COMFIXED_flag == 1) {
      type_mass = 0.0;
      for (j = 0; j < 3; j++) {
        fsum0[j] = 0.;
      }
      for (auto const &p : local_cells.particles()) {
        if (p.p.type == t0) {
          type_mass += (p).p.mass;
          for (j = 0; j < 3; j++) {
            fsum0[j] += p.f.f[j];
          }
        }
      }

      for (auto &p : local_cells.particles()) {
        if (p.p.type == t0) {
          for (j = 0; j < 3; j++) {
            p.f.f[j] -= fsum0[j] / type_mass * (p).p.mass;
          }
        }
      }
    }
  }
}

#endif

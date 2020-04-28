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
#ifndef ESPRESSO_OBSERVABLE_STAT_HPP
#define ESPRESSO_OBSERVABLE_STAT_HPP

#include <vector>

struct Observable_stat {
  /** Status flag for observable calculation. For 'analyze energy': 0
   *  re-initialize observable struct, else everything is fine,
   *  calculation can start. For 'analyze pressure' and 'analyze p_inst':
   *  0 or !(1+v_comp) re-initialize, else all OK.
   */
  int init_status;

  /** Array for observables on each node. */
  std::vector<double> data;

  /** number of Coulomb interactions */
  int n_coulomb;
  /** number of dipolar interactions */
  int n_dipolar;
  /** Number of virtual sites relative (rigid body) contributions */
  int n_virtual_sites;
  /** Number of external field contributions */
  const static int n_external_field = 1;

  /** start of bonded interactions. Right after the special ones */
  double *bonded;
  /** start of observables for non-bonded interactions. */
  double *non_bonded;
  /** start of observables for Coulomb interaction. */
  double *coulomb;
  /** start of observables for Coulomb interaction. */
  double *dipolar;
  /** Start of observables for virtual sites relative (rigid bodies) */
  double *virtual_sites;
  /** Start of observables for external fields */
  double *external_fields;

  /** number of doubles per data item */
  int chunk_size;
};

/** Structure used only in the pressure and stress tensor calculation to
   distinguish
    non-bonded intra- and inter- molecular contributions. */
struct Observable_stat_non_bonded {
  /** Array for observables on each node. */
  std::vector<double> data_nb;

  /** start of observables for non-bonded intramolecular interactions. */
  double *non_bonded_intra;
  /** start of observables for non-bonded intermolecular interactions. */
  double *non_bonded_inter;

  /** number of doubles per data item */
  int chunk_size_nb;
};

#endif // ESPRESSO_OBSERVABLE_STAT_HPP

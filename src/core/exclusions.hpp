/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef ESPRESSO_EXCLUSIONS_HPP
#define ESPRESSO_EXCLUSIONS_HPP

#include "config/config.hpp"

#include "Particle.hpp"

#include <algorithm>

#ifdef EXCLUSIONS

/** Determine if the non-bonded interactions between @p p1 and @p p2 should be
 *  calculated.
 */
inline bool do_nonbonded(Particle const &p1, Particle const &p2) {
  /* check for particle 2 in particle 1's exclusion list. The exclusion list is
   * symmetric, so this is sufficient. */
  return std::none_of(p1.exclusions().begin(), p1.exclusions().end(),
                      [&p2](int id) { return p2.id() == id; });
}

/** Remove exclusion from particle if possible */
void delete_exclusion(Particle &p, int p_id);

/** Insert an exclusion if not already set */
void add_exclusion(Particle &p, int p_id);

#endif // EXCLUSIONS
#endif // ESPRESSO_EXCLUSIONS_HPP

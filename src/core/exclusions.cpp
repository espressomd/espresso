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

#include "config.hpp"

#ifdef EXCLUSIONS

#include "exclusions.hpp"

#include <utils/contains.hpp>

#include <algorithm>

void add_exclusion(Particle *part, int part2) {
  if (Utils::contains(part->exclusions(), part2))
    return;

  part->exclusions().push_back(part2);
}

void delete_exclusion(Particle *part, int part2) {
  auto &el = part->exclusions();

  el.erase(std::remove(el.begin(), el.end(), part2), el.end());
}
#endif

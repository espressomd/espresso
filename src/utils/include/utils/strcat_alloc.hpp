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

#ifndef UTILS_STRCAT_ALLOC_HPP
#define UTILS_STRCAT_ALLOC_HPP

#include "memory.hpp"

#include <cstring>

namespace Utils {
/** Extend a string with another one.
 *  Like strcat, just automatically increases the string space.
 */
inline char *strcat_alloc(char *left, const char *right) {
  if (!right) {
    return left;
  }
  if (!left) {
    return strdup(right);
  }
  size_t newlen = strlen(left) + strlen(right);
  char *res = Utils::realloc(left, newlen + 1);
  strncat(res, right, newlen);
  return res;
}
} // namespace Utils

#endif

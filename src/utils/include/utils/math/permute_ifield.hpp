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
#ifndef UTILS_MATH_PERMUTE_IFIELD_HPP
#define UTILS_MATH_PERMUTE_IFIELD_HPP

namespace Utils {
/** permute an integer array field of size size about permute positions. */
inline void permute_ifield(int *field, int size, int permute) {
  if (permute == 0)
    return;
  if (permute < 0)
    permute = (size + permute);
  while (permute > 0) {
    int tmp = field[0];
    for (int i = 1; i < size; i++)
      field[i - 1] = field[i];
    field[size - 1] = tmp;
    permute--;
  }
}
} // namespace Utils

#endif

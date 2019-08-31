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
#ifndef UTILS_LIST_CONTAINS_HPP
#define UTILS_LIST_CONTAINS_HPP

#include <iterator>

/** @brief Check whether an iterator range contains a value.
 *
 * @param first Beginning of the range
 * @param last End of the range.
 * @param value The value to search for.
 *
 * @return True iff range contains the value.
 *
 * */
namespace Utils {
template <class InputIt, class T>
bool contains(InputIt first, InputIt last, T const &value) {
  return std::any_of(first, last, [value](T const &e) { return e == value; });
}

/** @brief Check whether an range contains a value.
 *
 * @param rng The range to search in.
 * @param value The value to search for.
 *
 * @return True iff range contains the value.
 *
 * */
template <class Range, class T>
bool contains(const Range &rng, T const &value) {
  using std::begin;
  using std::end;

  return contains(begin(rng), end(rng), value);
}
} // namespace Utils

#endif

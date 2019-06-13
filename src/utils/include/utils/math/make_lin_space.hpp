/*
Copyright (C) 2019 The ESPResSo project

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

#ifndef UTILS_MAKE_LIN_SPACE_HPP
#define UTILS_MAKE_LIN_SPACE_HPP

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

namespace Utils {
/**
 * @brief Equally spaced values in interval
 *
 * Returns a range of equally spaced values in
 * the range of start and stop, like numpy.linspace.
 *
 * @tparam T floating point type
 * @param start Start value of the interval
 * @param stop End value of the interval
 * @param number Number of partition points
 * @param endpoint If true, the last point is
 *        stop, otherwise one less.
 * @return Range of equally spaced values
 */
template <class T>
auto make_lin_space(T start, T stop, size_t number, bool endpoint = true) {
  using boost::make_counting_iterator;
  using boost::make_iterator_range;
  using boost::make_transform_iterator;

  auto const dx = (stop - start) / (number - endpoint);
  auto x = [dx, start](size_t i) { return start + i * dx; };

  return make_iterator_range(
      make_transform_iterator(make_counting_iterator(size_t(0)), x),
      make_transform_iterator(make_counting_iterator(number), x));
}
} // namespace Utils

#endif

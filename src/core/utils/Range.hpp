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
#ifndef UTILS_RANGE_HPP
#define UTILS_RANGE_HPP

#include <iterator>
#include <utility>

namespace Utils {

/**
 * @brief Helper class that wraps two iterators into
          a range object.
*/
template <typename Iterator> class Range {
  Iterator m_begin, m_end;

public:
  using iterator = Iterator;
  using value_type = typename std::iterator_traits<Iterator>::value_type;
  using difference_type =
      typename std::iterator_traits<Iterator>::difference_type;

  Range(Iterator begin, Iterator end)
      : m_begin(std::move(begin)), m_end(std::move(end)) {}

  Iterator begin() { return m_begin; }
  Iterator end() { return m_end; }

  bool empty() const { return m_begin == m_end; }
  difference_type size() const {
    using std::distance;
    return distance(m_begin, m_end);
  }

  bool operator==(Range const &rhs) const {
    return (m_begin == rhs.m_begin) && (m_end == rhs.m_end);
  }
};

/**
 * @brief Return a range for a pair of iterators.
 *
 * This is a convinence function so we can template
 * argument deduction to figure out the Range type.
 */
template <typename Iterator>
Range<Iterator> make_range(Iterator &&begin, Iterator &&end) {
  return Range<Iterator>(std::forward<Iterator>(begin),
                         std::forward<Iterator>(end));
}

} /* namespace Utils */

#endif

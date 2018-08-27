/*
  Copyright (C) 2010,2012,2013,2014,2015 The ESPResSo project
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

#ifndef UTILS_BATCH_HPP
#define UTILS_BATCH_HPP

#include <utility>

namespace Utils {
/**
 * @brief Convinience class to chain together multiple
 * Callables, like Funktors or lambdas.
 * Batch(A, B, ...)(args...) calls A(args...), B(args...), ... in order.
 */
template <typename... Callables> class Batch {
public:
  /**
   * @brief This operator has been intentionally left blank.
   */
  template <typename... Args> void operator()(Args &&...) { ; }
};

template <typename Callable, typename... Callables>
class Batch<Callable, Callables...> : public Batch<Callables...> {
public:
  Batch(Callable &&callable, Callables &&... callables)
      : Batch<Callables...>(std::forward<Callables>(callables)...),
        m_callable(std::forward<Callable>(callable)) {}

  template <typename... Args> void operator()(Args &&... args) {
    m_callable(std::forward<Args>(args)...);
    Batch<Callables...>::operator()(std::forward<Args>(args)...);
  }

private:
  Callable m_callable;
};

/**
 * @brief Construct a Batch instance with type deduced
 * from the function parameters.
 */
template <typename... Callables>
Batch<Callables...> make_batch(Callables &&... callables) {
  return Batch<Callables...>(std::forward<Callables>(callables)...);
}
} // namespace Utils

#endif

/*
  Copyright (C) 2017-2018 The ESPResSo project

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

#ifndef CORE_UTILS_PRINT_HPP
#define CORE_UTILS_PRINT_HPP

#include <iostream>

namespace Utils {

/**
 * @brief Python style print function.
 */
template <typename T> void print(T v) { std::cout << v << '\n'; }

/**
 * @brief Python style print function.
 */
template <typename Arg, typename... Args> void print(Arg v, Args... args) {
  std::cout << v << " ";
  print(args...);
}

} /* namespace Utils */

#endif

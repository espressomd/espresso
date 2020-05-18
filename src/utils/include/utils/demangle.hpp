/*
 * Copyright (C) 2019 The ESPResSo project
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
#ifndef UTILS_DEMANGLE_HPP
#define UTILS_DEMANGLE_HPP

#include <boost/core/demangle.hpp>

namespace Utils {
/**
 * @brief Get a human-readable name for a type.
 *
 * Uses boost to demangle the name, for details
 * see documentation for boost::core::demangle.
 *
 * @tparam T type
 * @return name
 */
template <class T> std::string demangle() {
  return boost::core::demangle(typeid(T).name());
}
} // namespace Utils

#endif // UTILS_DEMANGLE_HPP

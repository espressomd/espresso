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
#ifndef CORE_UTILS_MEMORY_HPP
#define CORE_UTILS_MEMORY_HPP

#include <cstddef>
#include <cstdlib>
#include <stdexcept>

namespace Utils {
/** used instead of malloc.
    Makes sure that a zero size allocation returns a nullptr pointer */
inline void *malloc(size_t size) {
  if (size == 0) {
    return nullptr;
  }

  void *p = ::malloc(size);

  if (p == nullptr) {
    throw std::bad_alloc{};
  }
  return p;
}
} // namespace Utils

#endif

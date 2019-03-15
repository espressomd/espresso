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
#ifndef CORE_UTILS_MEMORY_HPP
#define CORE_UTILS_MEMORY_HPP

#include <cstdlib>
#include <new>
#include <stdexcept>

namespace Utils {

/*************************************************************/
/** \name Dynamic memory allocation.                         */
/*************************************************************/
/*@{*/

/* to enable us to make sure that freed pointers are invalidated, we normally
   try to use realloc.
   Unfortunately allocating zero bytes (which should be avoided) actually
   allocates 16 bytes, and
   reallocating to 0 also. To avoid this, we use our own malloc and realloc
   procedures. */

/** used instead of realloc.
    Makes sure that resizing to zero FREEs pointer */
template <typename T> inline T *realloc(T *old, size_t size) {
  if (size == 0) {
    ::free(static_cast<void *>(old));
    return nullptr;
  }

  auto *p = static_cast<T *>(::realloc(static_cast<void *>(old), size));

  if (p == nullptr) {
    throw std::bad_alloc{};
  }
  return p;
}

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

/*@}*/

#endif

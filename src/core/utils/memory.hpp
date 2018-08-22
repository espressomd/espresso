#ifndef CORE_UTILS_MEMORY_HPP
#define CORE_UTILS_MEMORY_HPP

#include <new>
#include <stdexcept>
#include <stdlib.h>

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

  T *p = static_cast<T *>(::realloc(static_cast<void *>(old), size));

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
}

/*@}*/

#endif

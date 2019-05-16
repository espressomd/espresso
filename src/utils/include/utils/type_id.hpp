#ifndef UTILS_TYPE_ID_HPP
#define UTILS_TYPE_ID_HPP

#include <cstdint>

namespace Utils {

using typeid_t = void (*)();

/**
 * @brief Unique identifier for type.
 *
 * The uses the type-dependent address of the
 * function itself to derive a unique identifier
 * for the type. This *should* also work across
 * shared libraries.
 *
 * @tparam T Any type
 * @return A unique integral identifier for T.
 */
template <typename T> typeid_t type_id() {
  return typeid_t(type_id<T>);
}
} // namespace Utils

#endif

/*
  Copyright (C) 2017 The ESPResSo project

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

#ifndef CORE_UTILS_CACHE_HPP
#define CORE_UTILS_CACHE_HPP

#include <memory>
#include <type_traits>
#include <unordered_map>

namespace Utils {
template <typename Key, typename Value, typename Get> class Cache {
  using value_storage_type =
      std::unique_ptr<typename std::add_const<Value>::type>;
  using map_type = std::unordered_map<Key, value_storage_type>;

  static_assert(
      std::is_convertible<value_storage_type, decltype(std::declval<Get>()(
                                                  std::declval<Key>()))>::value,
      "Return value of Get::operator() needs to be convertible to "
      "value_storage_type.");

public:
  using key_type = Key;
  using value_type = const Value *;

public:
  template <typename GetRef>
  explicit Cache(GetRef &&getter) : m_getter(std::forward<GetRef>(getter)) {}

private:
  /** @brief The actual cache, maps id -> Particle copy.
   *
   * This is mutable so that access works on a constant
   * cache. This is fine because it does not change the
   * obsevable behaviour of the class.
  */
  mutable map_type m_cache;
  /** Value getter. */
  Get m_getter;

public:
  /** @brief Clear the cache.
   *
   * This invalidates all pointers into the cache.
   * It is not const on purpose because it can
   * have visible effects, even though m_cache is
   * mutable.
   */
  void invalidate() { m_cache.clear(); }

  /** @brief Query if k is contained in the cache. */
  bool has(Key const &k) const { return m_cache.find(k) != m_cache.end(); }

  /** @brief Get a value.
   *
   * If the value is not cached, it is fetched. If the fetching fails
   * e.g. if Get can not finde the requested element, a nullptr is
   * returned. The value is owned by the Cache and can not be modified.
   */
  Value const *get(Key const &k) const {
    auto const needle = m_cache.find(k);

    if (m_cache.end() != needle) {
      return needle->second.get();
    } else {
      auto const &val = m_cache[k] = m_getter(k);
      return val.get();
    }
  }
};

/**
 * @brief Type deduction helper for Cache.
 */
template <typename Key, typename Value, typename Get>
Cache<Key, Value, Get> make_cache(Get &&getter) {
  return Cache<Key, Value, typename std::remove_reference<Get>::type>(
      std::forward<Get>(getter));
}
}

#endif

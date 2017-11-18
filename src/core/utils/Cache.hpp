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
template <typename Key, typename Value> class Cache {
  using map_type =
      std::unordered_map<Key, typename std::add_const<Value>::type>;

public:
  using key_type = Key;
  using value_type = const Value *;

private:
  /** @brief The actual cache, maps id -> Particle copy.
   *
   * This is mutable so that access works on a constant
   * cache. This is fine because it does not change the
   * obsevable behaviour of the class.
  */
  mutable map_type m_cache;

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

  /** @brief Put a value into the cache. */
  Value const *put(Key const &k, Value &&v) {
    typename map_type::const_iterator it;
    std::tie(it, std::ignore) = m_cache.emplace(k, std::move(v));

    return &(it->second);
  }

  /** @brief Get a value.
   *
   * If the value is not cached, it is fetched. If the fetching fails
   * e.g. if Get can not finde the requested element, a nullptr is
   * returned. The value is owned by the Cache and can not be modified.
   */
  Value const *get(Key const &k) const {
    auto const needle = m_cache.find(k);

    if (m_cache.end() != needle) {
      return &(needle->second);
    } else {
      return nullptr;
    }
  }
};
}

#endif

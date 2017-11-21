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
#include <random>
#include <type_traits>
#include <unordered_map>

namespace Utils {
template <typename Key, typename Value> class Cache {
  using map_type =
      std::unordered_map<Key, typename std::add_const<Value>::type>;

public:
  using key_type = Key;
  using value_type = const Value *;
  using size_type = typename map_type::size_type;

  Cache() {
    /* Default to maximal size of our map. */
    m_max_size = m_cache.max_size();
  }
  explicit Cache(size_type max_size) : m_max_size(max_size), m_rand(max_size) {}

private:
  /** @brief The actual cache, maps Key -> Value.
   */
  map_type m_cache;
  size_type m_max_size;

  /** Cheap linear-congruential PRNG. */
  std::minstd_rand m_rand;

  /** @brief Drop an element.
   *
   * This function drops an element from the cache,
   * it is used when a new item is to be cached, but then
   * maximal size is reached. We use the bucket interface
   * here to be able to pick a random element in constant
   * time.
   */
  void drop_random_element() {
    if (m_cache.empty())
      return;

    auto const bucket_count = m_cache.bucket_count();

    /* Pick a random bucket, this has to terminate because
    * the map is not empty. So there has to be a bucket with
    * at least one element. */
    auto bucket =
        std::uniform_int_distribution<size_type>{0, bucket_count - 1}(m_rand);

    while (0 == m_cache.bucket_size(bucket)) {
      /* Wrap around the end */
      bucket = ((bucket + 1) % bucket_count);
    }

    /* Pick a random elemnt form that bucket. */
    auto const elem_in_bucket = std::uniform_int_distribution<size_type>{
        0, m_cache.bucket_size(bucket) - 1}(m_rand);

    /* Get the element in the bucket */
    auto const drop_key =
        std::next(m_cache.cbegin(bucket), elem_in_bucket)->first;

    /* And drop it. */
    m_cache.erase(drop_key);
  }

public:
  /** @brief Clear the cache.
   *
   * This invalidates all pointers into the cache.
   */
  void invalidate() { m_cache.clear(); }

  /** @brief Query if k is contained in the cache. */
  bool has(Key const &k) const { return m_cache.find(k) != m_cache.end(); }

  /** @brief Number of elements currently cached. */
  size_type size() const { return m_cache.size(); }

  /** @brief Maximal size of the cache. */
  size_type max_size() const { return m_max_size; }

  /** @brief Put a value into the cache.
   *
   * If the value already exists, it is overwritten.
   * When the size of the cache would grow below the
   * maximal size, a random element is removed before
   * putting the new one. */
  template <typename ValueRef> Value const *put(Key const &k, ValueRef &&v) {
    /* If there already is a value for k, overwriting it
     * will not increase the size, so we don't have to
     * make room. */
    if ((m_cache.size() >= m_max_size) && !has(k))
      drop_random_element();

    typename map_type::const_iterator it;
    std::tie(it, std::ignore) = m_cache.emplace(k, std::forward<ValueRef>(v));

    return &(it->second);
  }

  /** @brief Put a range of values into the cache.
   *
   * If the values already exists, it is overwritten.
   * When the size of the cache would grow below the
   * maximal size, a random elements are removed until
   * all of the new values fit. If the given range is
   * larger than max_size(), only the first max_size()
   * elements are put into the caache.
   *
   * @tparam KeyInputIterator iterator of keys, at least InputIterator.
   * @tparam ValueInputIterator iterator of value, at least InputIterator.
   *
   * @returns KeyInputIterator one past the last element that was put
   *          into the cache.
   */
  template <typename KeyInputIterator, typename ValueInputIterator>
  KeyInputIterator put(KeyInputIterator kbegin, KeyInputIterator kend, ValueInputIterator vbegin) {
    auto const range_len = std::distance(kbegin, kend);
    auto const len = (range_len > max_size()) ? max_size() : range_len;
    kend = std::next(kbegin, len);

    /* Make some space. */
    while ((max_size() - size()) < len) {
      drop_random_element();
    }

    while(kbegin != kend) {
      put(*kbegin++, *vbegin++);
    }

    return kend;
  }

  /** @brief Get a value.
   *
   *  The value is owned by the Cache and can not be modified.
   *  Pointers into the cache can be invalidated at any point
   *  and should not be stored beyond the calling function.
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

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
#ifndef ESPRESSO_MEMCPY_ARCHIVE_HPP
#define ESPRESSO_MEMCPY_ARCHIVE_HPP

#include <utils/Span.hpp>

#include <boost/mpl/bool.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/serialization.hpp>

#include <cstring>

namespace Utils {
/** @brief Type trait to indicate that a type is
 *         serializable with a static size, e.g. is
 *         suitable for memcpy serialization. Only
 *         specialize this to std::true_type if it is
 *         guarantueed that serializing this type always
 *         returns the same number of bytes, independent
 *         of object state.
 *
 * @tparam T type under consideration.
 */
template <class T> struct is_statically_serializable : std::false_type {};

namespace detail {
/* Use serialize function only if the type is opt-in but not
 * trivially copyable, in which case memcpy is more efficient. */
template <class T>
using use_serialize =
    std::integral_constant<bool, not std::is_trivially_copyable<T>::value and
                                     is_statically_serializable<T>::value>;

struct SizeArchive {
  using is_saving = boost::mpl::true_;
  using is_loading = boost::mpl::false_;

  size_t size = {};

  template <typename T>
  constexpr auto operator<<(T const &)
      -> std::enable_if_t<std::is_trivially_copyable<T>::value> {
    size += sizeof(T);
  }

  template <class T>
  constexpr auto operator<<(T &v)
      -> std::enable_if_t<detail::use_serialize<T>::value> {
    boost::serialization::serialize(*this, v, 0);
  }

  template <class T> constexpr auto operator<<(boost::optional<T> &) {
    operator<<(bool());
    operator<<(T{});
  }

  template <class T> constexpr SizeArchive &operator&(T &t) {
    operator<<(t);

    return *this;
  }
};

class BasicMemcpyArchive {
  /** Buffer to write to */
  Utils::Span<char> buf;
  /** Current position in the buffer */
  char *insert;

public:
  explicit BasicMemcpyArchive(Utils::Span<char> buf)
      : buf(buf), insert(buf.data()) {}

  size_t bytes_processed() const { return insert - buf.data(); }

  template <typename T>
  auto operator>>(T &value)
      -> std::enable_if_t<std::is_trivially_copyable<T>::value> {
    /* check that there is enough space left in the buffer */
    assert((insert + sizeof(T)) <= buf.end());

    std::memcpy(&value, insert, sizeof(T));
    insert += sizeof(T);
  }

  template <typename T>
  auto operator<<(T const &value)
      -> std::enable_if_t<std::is_trivially_copyable<T>::value> {
    /* check that there is enough space left in the buffer */
    assert((insert + sizeof(T)) <= buf.end());
    std::memcpy(insert, &value, sizeof(T));
    insert += sizeof(T);
  }

  /**
   * @brief Determine the static packing size of a type.
   * @tparam T Type to consider.
   * @return Packed size in bytes.
   */
  template <class T> static size_t packing_size() {
    detail::SizeArchive sa{};
    T t;
    return (sa & t).size;
  }
};
} // namespace detail

/**
 * @brief Archive that deserializes from a buffer via memcpy.
 *
 * Can only process types that have a static serialization size,
 * e.g. that serialize to the same number of bytes independent of
 * the state of the object. This can either be automatically detected
 * for types that are trivially copyable, or by explicitly assuring
 * this by specializing @c is_statically_serializable to std::true_type.
 */
class MemcpyIArchive : public detail::BasicMemcpyArchive {
private:
  using base_type = detail::BasicMemcpyArchive;

public:
  using is_loading = boost::mpl::true_;
  using is_saving = boost::mpl::false_;

  /**
   * @param buf Buffer to read from.
   */
  explicit MemcpyIArchive(Utils::Span<char> buf)
      : detail::BasicMemcpyArchive(buf) {}

  /**
   * @brief Number of bytes read from the buffer.
   * @return Number of bytes read.
   */
  size_t bytes_read() const { return bytes_processed(); }

  /** @copydoc base_type::packing_size */
  using base_type::packing_size;
  using base_type::operator>>;

  template <class T>
  auto operator>>(T &value)
      -> std::enable_if_t<detail::use_serialize<T>::value> {
    boost::serialization::serialize(*this, value, 0);
  }

  template <class T> void operator>>(boost::optional<T> &o) {
    bool initialized{};
    operator>>(initialized);
    T val{};
    operator>>(val);

    if (initialized) {
      o = val;
    } else {
      o = boost::none;
    }
  }

  template <class T> MemcpyIArchive &operator&(T &value) {
    operator>>(value);

    return *this;
  }
};

/**
 * @brief Archive that serializes to a buffer via memcpy.
 *
 * @copydetails MemcpyIArchive
 */
class MemcpyOArchive : public detail::BasicMemcpyArchive {
  using base_type = detail::BasicMemcpyArchive;

public:
  using is_loading = boost::mpl::false_;
  using is_saving = boost::mpl::true_;

  /**
   * @param buf Buffer to write to.
   */
  explicit MemcpyOArchive(Utils::Span<char> buf)
      : detail::BasicMemcpyArchive(buf) {}

  /**
   * @brief Number of bytes written to the buffer.
   * @return Number of bytes written.
   */
  size_t bytes_written() const { return bytes_processed(); }

  /** @copydoc base_type::packing_size */
  using base_type::packing_size;
  using base_type::operator<<;

  template <typename T>
  auto operator<<(T &value)
      -> std::enable_if_t<detail::use_serialize<T>::value> {
    boost::serialization::serialize(*this, value, 0);
  }

  template <class T> void operator<<(boost::optional<T> &o) {
    operator<<(o.is_initialized());
    operator<<(o.value_or(T{}));
  }

  template <class T> MemcpyOArchive &operator&(T &value) {
    operator<<(value);

    return *this;
  }
};
} // namespace Utils

#endif // ESPRESSO_MEMCPY_ARCHIVE_HPP

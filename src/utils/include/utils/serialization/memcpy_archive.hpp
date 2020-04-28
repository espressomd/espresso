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
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>

#include <cstring>

namespace Utils {
/** @brief Type trait to indicate that a type is
 *         serializable with a static size, e.g. is
 *         suitable for memcpy serialization. Only
 *         specialize this to std::true_type if it is
 *         guaranteed that serializing this type always
 *         returns the same number of bytes, independent
 *         of object state.
 *
 * @tparam T type under consideration.
 */
template <class T>
struct is_statically_serializable
    : std::integral_constant<
          bool, std::is_trivially_copyable<T>::value or
                    boost::serialization::is_bitwise_serializable<T>::value> {};

namespace detail {
/* Use memcpy for packing */
template <class T>
using use_memcpy = std::integral_constant<
    bool, std::is_trivially_copyable<T>::value or
              boost::serialization::is_bitwise_serializable<T>::value>;
/* Use serialize function only if the type is opt-in but not
 * trivially copyable, in which case memcpy is more efficient. */
template <class T>
using use_serialize =
    std::integral_constant<bool, not use_memcpy<T>::value and
                                     is_statically_serializable<T>::value>;

template <class Derived> class BasicMemcpyArchive {
  /** Buffer to write to */
  Utils::Span<char> buf;
  /** Current position in the buffer */
  char *insert;

public:
  explicit BasicMemcpyArchive(Utils::Span<char> buf)
      : buf(buf), insert(buf.data()) {}

  size_t get_library_version() const { return 4; }

  size_t bytes_processed() const { return insert - buf.data(); }
  void skip(size_t bytes) {
    assert((insert + bytes) <= buf.end());
    insert += bytes;
  }

private:
  void read(void *data, size_t bytes) {
    /* check that there is enough space left in the buffer */
    assert((insert + bytes) <= buf.end());
    std::memcpy(data, insert, bytes);
    insert += bytes;
  }

  void write(const void *data, size_t bytes) {
    /* check that there is enough space left in the buffer */
    assert((insert + bytes) <= buf.end());
    std::memcpy(insert, data, bytes);
    insert += bytes;
  }

public:
  template <typename T>
  auto operator>>(T &value) -> std::enable_if_t<use_memcpy<T>::value> {
    read(&value, sizeof(T));
  }

  template <typename T>
  auto operator<<(T const &value) -> std::enable_if_t<use_memcpy<T>::value> {
    write(&value, sizeof(T));
  }

private:
  template <class T> void process(T &value) {
    auto const old_pos = insert;
    boost::serialization::serialize_adl(*static_cast<Derived *>(this), value,
                                        4);
    auto const new_pos = insert;
    assert((new_pos - old_pos) <= sizeof(T));

    auto const padding_size = sizeof(T) - (new_pos - old_pos);
    skip(padding_size);
  }

public:
  template <class T>
  auto operator>>(T &value)
      -> std::enable_if_t<detail::use_serialize<T>::value> {
    process(value);
  }

  template <class T>
  auto operator<<(T &value)
      -> std::enable_if_t<detail::use_serialize<T>::value> {
    process(value);
  }

  template <class T> void operator<<(const boost::serialization::nvp<T> &nvp) {
    operator<<(nvp.const_value());
  }

  template <class T> void operator>>(const boost::serialization::nvp<T> &nvp) {
    operator>>(nvp.value());
  }

  /**
   * @brief Determine the static packing size of a type.
   * @tparam T Type to consider.
   * @return Packed size in bytes.
   */
  template <class T> static constexpr size_t packing_size() {
    return sizeof(T);
  }
};
} // namespace detail

/**
 * @brief Archive that deserializes from a buffer via memcpy.
 *
 * Can only process types that have a static serialization size,
 * e.g. that serialize to the same number of bytes independent of
 * the state of the object. This can either be automatically detected
 * for types that are trivially copyable, or by explicitly assured
 * by specializing @c is_statically_serializable to std::true_type.
 */
class MemcpyIArchive : public detail::BasicMemcpyArchive<MemcpyIArchive> {
private:
  using base_type = detail::BasicMemcpyArchive<MemcpyIArchive>;

public:
  using is_loading = boost::mpl::true_;
  using is_saving = boost::mpl::false_;

  /**
   * @param buf Buffer to read from.
   */
  explicit MemcpyIArchive(Utils::Span<char> buf) : base_type(buf) {}

  /**
   * @brief Number of bytes read from the buffer.
   * @return Number of bytes read.
   */
  size_t bytes_read() const { return bytes_processed(); }

  /** @copydoc base_type::packing_size */
  using base_type::packing_size;
  using base_type::operator>>;

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
class MemcpyOArchive : public detail::BasicMemcpyArchive<MemcpyOArchive> {
  using base_type = detail::BasicMemcpyArchive<MemcpyOArchive>;

public:
  using is_loading = boost::mpl::false_;
  using is_saving = boost::mpl::true_;

  /**
   * @param buf Buffer to write to.
   */
  explicit MemcpyOArchive(Utils::Span<char> buf) : base_type(buf) {}

  /**
   * @brief Number of bytes written to the buffer.
   * @return Number of bytes written.
   */
  size_t bytes_written() const { return bytes_processed(); }

  /** @copydoc base_type::packing_size */
  using base_type::packing_size;
  using base_type::operator<<;

  template <class T> MemcpyOArchive &operator&(T &value) {
    operator<<(value);

    return *this;
  }
};
} // namespace Utils

#endif // ESPRESSO_MEMCPY_ARCHIVE_HPP

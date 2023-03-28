/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE serialization of containers
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/* This <boost/serialization/version.hpp> include guards against an issue
 * in boost::serialization from boost 1.74.0 that leads to compiler error
 * "explicit specialization of undeclared template struct 'version'" when
 * including <boost/serialization/optional.hpp>. More details in tickets:
 * https://github.com/boostorg/serialization/issues/210
 * https://github.com/boostorg/serialization/issues/217
 */
#include <boost/serialization/version.hpp>

#include <utils/Array.hpp>
#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/quaternion.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/packed_iarchive.hpp>
#include <boost/mpi/packed_oarchive.hpp>
#include <boost/serialization/array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <sstream>
#include <vector>

UTILS_ARRAY_BOOST_CLASS(std::array, N, object_serializable)
UTILS_ARRAY_BOOST_TRACK(std::array, N, track_never)

namespace boost {
namespace mpi {
using packed_archive = packed_oarchive::buffer_type;
} // namespace mpi
} // namespace boost

/* Data types and reference values used in the testsuite. */
namespace Testing {
using T = int;
constexpr std::array<T, 4> values{{4, 5, 6, 7}};
constexpr auto N = values.size();
using Serial_T = boost::mpi::packed_oarchive::buffer_type::value_type;
} // namespace Testing

/* Convenience functions to handle different container types. */
namespace Testing {

/** Instantiate a container. */
template <template <typename, std::size_t> class Container>
auto make_container() {
  return Container<Testing::T, Testing::N>{};
}

template <template <typename> class Container> auto make_container() {
  return Container<Testing::T>{};
}

/** Initialize an array with reference values. */
template <template <typename, std::size_t...> class A, typename T,
          std::size_t... N>
void constexpr init_container(A<T, N...> &array) {
  for (std::size_t i = 0; i < Testing::N; ++i) {
    array[i] = Testing::values[i];
  }
}

template <typename T, std::size_t N>
void constexpr init_container(Utils::detail::Storage<T, N> &array) {
  for (std::size_t i = 0; i < Testing::N; ++i) {
    array.m_data[i] = Testing::values[i];
  }
}

template <typename T>
void constexpr init_container(Utils::compact_vector<T> &array) {
  array.resize(Testing::N);
  for (std::size_t i = 0; i < Testing::N; ++i) {
    array[i] = Testing::values[i];
  }
}

/** Convert an array to a container type that provides method @c operator[] */
template <template <typename, std::size_t...> class A, typename T,
          std::size_t... N>
auto constexpr testable_container(A<T, N...> const &array) {
  return array;
}

template <typename T, std::size_t N>
auto constexpr testable_container(Utils::detail::Storage<T, N> const &s) {
  return Utils::Array<T, N>{s};
}

template <typename T>
auto constexpr testable_container(Utils::Quaternion<T> const &s) {
  return s.m_data;
}

template <typename T>
auto constexpr testable_container(Utils::compact_vector<T> const &s) {
  assert(s.size() == Testing::N);
  Utils::Array<T, Testing::N> out{};
  for (std::size_t i = 0; i < Testing::N; ++i) {
    out[i] = s[i];
  }
  return out;
}

/** Convert a string buffer to a string vector. */
auto to_vector(std::stringstream &buffer_in) {
  std::vector<std::stringstream::char_type> buffer_out;
  buffer_in.seekg(0);
  while (!buffer_in.eof()) {
    std::stringstream::char_type c = '\0';
    buffer_in.get(c);
    buffer_out.push_back(c);
  }
  return buffer_out;
}

/** Convert a bytestring buffer to a bytestring vector. */
auto to_vector(boost::mpi::packed_archive const &buffer_in) {
  std::vector<Serial_T> buffer_out;
  for (auto const byte : buffer_in) {
    buffer_out.push_back(byte);
  }
  return buffer_out;
}
} // namespace Testing

/** Serialize an array using a string representation.
 *  Check if the deserialized array is equal to the original array and if
 *  the original array remains untouched. Return the string representation
 *  of the serialization buffer.
 */
template <template <typename, std::size_t...> class Container>
auto create_text_archive() {
  std::stringstream buffer;
  auto values_send = Testing::make_container<Container>();
  auto values_recv = Testing::make_container<Container>();
  Testing::init_container(values_send);
  boost::archive::text_oarchive oa{buffer};
  oa << values_send;
  boost::archive::text_iarchive ia{buffer};
  ia >> values_recv;
  auto const test_send = Testing::testable_container(values_send);
  auto const test_recv = Testing::testable_container(values_recv);
  auto const test_ref = Testing::testable_container(Testing::values);
  BOOST_TEST(test_send == test_ref, boost::test_tools::per_element());
  BOOST_TEST(test_recv == test_ref, boost::test_tools::per_element());
  return Testing::to_vector(buffer);
}

/** Serialize an array using a bytestring representation.
 *  Check if the deserialized array is equal to the original array and if
 *  the original array remains untouched. Return the bytestring representation
 *  of the serialization buffer.
 */
template <template <typename, std::size_t...> class Container>
auto create_mpi_archive(boost::mpi::communicator const &comm) {
  boost::mpi::packed_archive buffer;
  auto values_send = Testing::make_container<Container>();
  auto values_recv = Testing::make_container<Container>();
  Testing::init_container(values_send);
  boost::mpi::packed_oarchive oa{comm, buffer};
  oa << values_send;
  boost::mpi::packed_iarchive ia{comm, buffer};
  ia >> values_recv;
  auto const test_send = Testing::testable_container(values_send);
  auto const test_recv = Testing::testable_container(values_recv);
  auto const test_ref = Testing::testable_container(Testing::values);
  BOOST_TEST(test_send == test_ref, boost::test_tools::per_element());
  BOOST_TEST(test_recv == test_ref, boost::test_tools::per_element());
  return Testing::to_vector(buffer);
}

template <std::size_t Length, class InputIt>
auto sorted_view(InputIt const &buffer_it) {
  std::array<Testing::Serial_T, Length> subset;
  std::copy_n(buffer_it, Length, subset.begin());
  std::sort(subset.begin(), subset.end(), std::greater<>());
  return subset;
}

/**
 * @brief Simplistic test for endianness.
 * Replace with @c std::endian once ESPResSo becomes a C++20 project.
 */
bool is_big_endian() {
#ifdef __BYTE_ORDER__
  return __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__;
#else
  return false;
#endif
}

BOOST_AUTO_TEST_CASE(serialization_level_test) {
  boost::mpi::communicator comm;
  auto const buffer = create_mpi_archive<std::array>(comm);

  /* Serialization should produce the following bytestring (little-endian):
   *   3 0 0 0 0 0 0 0
   *   4 0 0 0
   *   5 0 0 0
   *   6 0 0 0
   *   7 0 0 0
   * Since the order is not guaranteed, we extract views of the bytestring and
   * sort them in descending order to easily compare them to reference values.
   */
  auto constexpr metadata_size = sizeof(std::size_t);
  auto constexpr data_size = Testing::N * sizeof(Testing::T);
  auto constexpr buffer_size_ref = metadata_size + data_size;
  BOOST_REQUIRE_EQUAL(buffer.size(), buffer_size_ref);

  // check buffer metadata
  {
    static_assert(Testing::N <= (1ul << (8ul * sizeof(Testing::Serial_T))));
    auto constexpr N = sizeof(std::size_t);
    auto const array_length = sorted_view<N>(buffer.begin());
    // since the array length is smaller than the size of the buffer type,
    // the size appears once in the bytestring and all other bytes are 0
    auto const array_length_ref = std::array<Testing::Serial_T, N>{
        {static_cast<Testing::Serial_T>(Testing::N)}};
    BOOST_TEST(array_length == array_length_ref,
               boost::test_tools::per_element());
  }

  // check buffer data
  for (std::size_t i = 0; i < Testing::N; ++i) {
    auto constexpr N = sizeof(Testing::T);
    auto const offset = sizeof(std::size_t) + i * N;
    auto const array_data =
        sorted_view<N>(buffer.begin() + static_cast<long>(offset));
    auto const array_data_ref = std::array<Testing::Serial_T, N>{
        {static_cast<Testing::Serial_T>(Testing::values[i])}};
    BOOST_TEST(array_data == array_data_ref, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_CASE(mpi_archive_test) {
  boost::mpi::communicator comm;
  auto constexpr data_size = Testing::N * sizeof(Testing::T);
  auto buffer_std_array = create_mpi_archive<std::array>(comm);
  decltype(buffer_std_array) buffer_ref(data_size);
  std::copy_n(buffer_std_array.data() + sizeof(std::size_t), data_size,
              buffer_ref.data());
  auto const buffer_array = create_mpi_archive<Utils::Array>(comm);
  auto const buffer_vector = create_mpi_archive<Utils::Vector>(comm);
  auto const buffer_storage = create_mpi_archive<Utils::detail::Storage>(comm);
  auto const buffer_quat = create_mpi_archive<Utils::Quaternion>(comm);
  auto buffer_cv = create_mpi_archive<Utils::compact_vector>(comm);
  BOOST_TEST(buffer_array == buffer_ref, boost::test_tools::per_element());
  BOOST_TEST(buffer_vector == buffer_ref, boost::test_tools::per_element());
  BOOST_TEST(buffer_storage == buffer_ref, boost::test_tools::per_element());
  BOOST_TEST(buffer_quat == buffer_ref, boost::test_tools::per_element());
  auto const index_lsb = (is_big_endian()) ? 1 : 0;
  auto const index_hsb = (is_big_endian()) ? 0 : 1;
  BOOST_TEST(buffer_cv[index_lsb] == Testing::N);
  BOOST_TEST(buffer_cv[index_hsb] == 0);
  buffer_cv.erase(buffer_cv.begin());
  buffer_cv.erase(buffer_cv.begin());
  BOOST_TEST(buffer_cv == buffer_ref, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(text_archive_test) {
  auto const buffer_array = create_text_archive<Utils::Array>();
  auto const buffer_vector = create_text_archive<Utils::Vector>();
  auto const buffer_storage = create_text_archive<Utils::detail::Storage>();
  auto const buffer_quat = create_text_archive<Utils::Quaternion>();
  auto const buffer_cv = create_text_archive<Utils::compact_vector>();
  BOOST_TEST(buffer_array == buffer_storage, boost::test_tools::per_element());
  BOOST_TEST(buffer_vector == buffer_storage, boost::test_tools::per_element());
  BOOST_TEST(buffer_quat == buffer_storage, boost::test_tools::per_element());
  BOOST_TEST(buffer_cv == buffer_storage, boost::test_tools::per_element());
}

template <template <class> class Trait> void constexpr assert_has_trait() {
  struct S {};
  // arrays of standard types have the trait
  static_assert(Trait<Utils::detail::Storage<int, 2>>::value);
  static_assert(Trait<Utils::Array<int, 2>>::value);
  static_assert(Trait<Utils::Vector<int, 2>>::value);
  static_assert(Trait<Utils::Vector3i>::value);
  static_assert(Trait<Utils::Vector3d>::value);
  static_assert(Trait<Utils::Quaternion<float>>::value);
  BOOST_TEST_PASSPOINT();
  // arrays of non-standard types don't have the trait
  static_assert(!Trait<S>::value);
  static_assert(!Trait<Utils::detail::Storage<S, 2>>::value);
  static_assert(!Trait<Utils::Array<S, 2>>::value);
  static_assert(!Trait<Utils::Vector<S, 2>>::value);
  static_assert(!Trait<Utils::Quaternion<S>>::value);
}

BOOST_AUTO_TEST_CASE(serialization_traits_test) {
  // arrays of MPI datatypes can be serialized as MPI datatypes, but
  // arrays of custom datatypes cannot be serialized as MPI datatypes
  assert_has_trait<boost::mpi::is_mpi_datatype>();
  // arrays of bitwise serializable MPI datatypes are bitwise serializable,
  // but arrays of custom datatypes are not bitwise serializable
  assert_has_trait<boost::serialization::is_bitwise_serializable>();
}

BOOST_AUTO_TEST_CASE(compact_vector_test) {
  using namespace Utils;
  static_assert(boost::mpi::is_mpi_datatype<compact_vector<int>>::value);
  static_assert(sizeof(compact_vector<int>) < sizeof(std::vector<int>));
  BOOST_TEST_PASSPOINT();
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}

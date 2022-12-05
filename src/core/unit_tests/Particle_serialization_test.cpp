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

#define BOOST_TEST_MODULE Particle serialization test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "config/config.hpp"

#include <utils/demangle.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/hana.hpp>
#include <boost/mpl/list.hpp>

#include <cstddef>
#include <regex>
#include <type_traits>
#include <utility>
#include <vector>

namespace traits {

namespace detail {

/**
 * @brief Helper class to detect the presence of a member in a type.
 * If a member exists in both @p T and @p Ref, then it's impossible
 * to access the member in the derived class due to the ambiguity.
 * If the member only exists in @p Ref, then it can be accessed in
 * the derived class. This can be leveraged by SFINAE, even if the
 * member is private. To handle all C++ types, use @c std::is_class
 * in the trait specialization. Only works in a constexpr context!
 * @tparam T     Class to check
 * @tparam Ref   Reference class that contains the member
 */
template <typename T, typename Ref>
struct DetectMember : public T, public Ref {};

struct SerializableClass {
  template <class Archive> void serialize(Archive &, long int);
};

auto hasnt_serialize_method = boost::hana::is_valid(
    [](auto &&x) -> decltype((void)x.serialize(std::declval<int &>(),
                                               std::declval<long int>())) {});

} // namespace detail

/**
 * Does the type contain a <tt>serialize(Archive &, long int)</tt> method.
 */
template <class T, typename Enable = void>
struct has_serialize_method : std::integral_constant<bool, false> {};

template <class T>
struct has_serialize_method<T, typename std::enable_if_t<std::is_class_v<T>>>
    : std::integral_constant<
          bool, !static_cast<bool>(detail::hasnt_serialize_method(
                    detail::DetectMember<T, detail::SerializableClass>{}))> {};

} // namespace traits

/**
 * Recursively walk through the serialization tree of an object and keep
 * track of which serialized types are not matching a specific trait.
 * @tparam Trait    Trait to check for
 */
template <template <class> class Trait> class TraitChecker {
public:
  typedef boost::mpl::bool_<true> is_saving;
  typedef boost::mpl::bool_<false> is_loading;
  typedef std::vector<std::string> buffer_type;

  TraitChecker(buffer_type &buffer) : m_buffer(buffer) {}

  template <class T> auto &operator<<(T &t) {
    if (not Trait<T>::value)
      m_buffer.emplace_back(Utils::demangle<T>());
    else
      recurse(t);
    return *this;
  }

  template <class T> auto &operator&(T &t) { return *this << t; }

  /**
   * Special method that skips the first level in the serialization tree.
   * Allow checking for a trait in members of a type T, even if T itself
   * doesn't have this trait.
   */
  template <class T> auto &operator|(T &t) {
    recurse(t);
    return *this;
  }

  void save_binary(void *, std::size_t) {}

private:
  template <typename T> void recurse(T &t) {
    // standard types are not classes
    if constexpr (std::is_class_v<T>) {
      m_accessor.serialize(*this, t, 0);
    }
  }

  // const objects cannot be processed by boost::serialization
  template <typename T> void recurse(T const &) {}

  // std::vector objects cannot be processed by boost::serialization
  template <typename T> void recurse(std::vector<T> &) {}

  buffer_type &m_buffer;
  boost::serialization::access m_accessor;
};

namespace Testing {

struct NotSerializable {};

class BitwiseSerializable {
  int const a = 1;
  long double b = 2.l;
  int const c = 3;
  long double d = 4.l;

  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int) {
    ar &a &b;
    ar << c << d;
  }
};

class NotBitwiseSerializable {
  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int) {}
};

class MixedSerializable {
  BitwiseSerializable a;
  NotBitwiseSerializable b;

  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int) {
    ar << a << b;
  }
};

} // namespace Testing

BOOST_IS_BITWISE_SERIALIZABLE(Testing::BitwiseSerializable)

void constexpr assert_has_serialize_method() {
  using namespace Testing;
  static_assert(!traits::has_serialize_method<int>::value);
  static_assert(!traits::has_serialize_method<NotSerializable>::value);
  static_assert(traits::has_serialize_method<BitwiseSerializable>::value);
  static_assert(traits::has_serialize_method<MixedSerializable>::value);
}

BOOST_AUTO_TEST_CASE(TraitChecker_test) {
  using boost::serialization::is_bitwise_serializable;
  using Checker = TraitChecker<is_bitwise_serializable>;
  assert_has_serialize_method();
  Checker::buffer_type buffer;
  Checker oa{buffer};
  Testing::BitwiseSerializable serializable;
  oa &serializable;
  BOOST_REQUIRE_EQUAL(buffer.size(), 0);
  Testing::MixedSerializable mixed;
  oa | mixed;
  BOOST_REQUIRE_EQUAL(buffer.size(), 1);
  BOOST_REQUIRE_EQUAL(buffer[0], "Testing::NotBitwiseSerializable");
}

namespace trait_wrappers {
struct is_trivially_copyable {
  template <typename T> using apply = std::is_trivially_copyable<T>;
};
struct is_statically_serializable {
  template <typename T> using apply = Utils::is_statically_serializable<T>;
};
} // namespace trait_wrappers

BOOST_AUTO_TEST_CASE_TEMPLATE(
    particle_member_traits_test, TraitWrapper,
    BOOST_IDENTITY_TYPE(
        (boost::mpl::list<trait_wrappers::is_trivially_copyable,
                          trait_wrappers::is_statically_serializable>))) {

  using Checker = TraitChecker<TraitWrapper::template apply>;

  // struct Particle is not statically serializable
  {
    static_assert(!TraitWrapper::template apply<Particle>::value);
    typename Checker::buffer_type buffer_ref = {"Particle"};
    typename Checker::buffer_type buffer;
    Checker oa{buffer};
    Particle p;
    oa &p;
    BOOST_TEST(buffer == buffer_ref, boost::test_tools::per_element());
  }

  // most members of Particle are statically serializable
  {
    typename Checker::buffer_type buffer_ref = {
        "BondList",
#ifdef EXCLUSIONS
        "Utils::compact_vector<int>",
#endif
    };
    typename Checker::buffer_type buffer;
    Checker oa{buffer};
    Particle p;
    oa | p;
    std::transform(buffer.begin(), buffer.end(), buffer.begin(),
                   [](std::string const &symbol) {
                     return std::regex_replace(symbol, std::regex("std::__1::"),
                                               "std::");
                   });
    BOOST_TEST(buffer == buffer_ref, boost::test_tools::per_element());
  }
}

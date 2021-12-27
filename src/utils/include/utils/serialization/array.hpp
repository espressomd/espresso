/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef UTILS_SERIALIZATION_ARRAY_HPP
#define UTILS_SERIALIZATION_ARRAY_HPP

#include <boost/mpl/bool_fwd.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/integral_c_tag.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/tracking_enum.hpp>

#include <cstddef>

#define UTILS_ARRAY_TEMPLATE_T_N template <typename T, std::size_t N>
#define UTILS_ARRAY_TEMPLATE_T_0 template <typename T>
#define UTILS_ARRAY_CONTAINER_T_N(Container) Container<T, N>
#define UTILS_ARRAY_CONTAINER_T_0(Container) Container<T>

// forward declare
namespace boost {
namespace mpi {
template <typename T> struct is_mpi_datatype;
} /* namespace mpi */
} /* namespace boost */

/**
 * @brief Mark array types as MPI data types.
 * @tparam Container            Template template type of the array
 * @tparam N                    N if @p Container uses std::size_t N, else 0
 */
#define UTILS_ARRAY_BOOST_MPI_T(Container, N)                                  \
  namespace boost {                                                            \
  namespace mpi {                                                              \
  UTILS_ARRAY_TEMPLATE_T_##N struct is_mpi_datatype<                           \
      UTILS_ARRAY_CONTAINER_T_##N(Container)> : public is_mpi_datatype<T> {};  \
  } /* namespace mpi */                                                        \
  } /* namespace boost */

/**
 * @brief Mark array types as MPI bitwise serializable.
 * @tparam Container            Template template type of the array
 * @tparam N                    N if @p Container uses std::size_t N, else 0
 */
#define UTILS_ARRAY_BOOST_BIT_S(Container, N)                                  \
  namespace boost {                                                            \
  namespace serialization {                                                    \
  UTILS_ARRAY_TEMPLATE_T_##N struct is_bitwise_serializable<                   \
      UTILS_ARRAY_CONTAINER_T_##N(Container)>                                  \
      : public is_bitwise_serializable<T> {};                                  \
  } /* namespace serialization */                                              \
  } /* namespace boost */

/**
 * @brief Redefinition of @c BOOST_CLASS_IMPLEMENTATION for array types.
 * @tparam Container            Template template type of the array
 * @tparam N                    N if @p Container uses std::size_t N, else 0
 * @tparam ImplementationLevel  Serialization implementation level
 */
#define UTILS_ARRAY_BOOST_CLASS(Container, N, ImplementationLevel)             \
  namespace boost {                                                            \
  namespace serialization {                                                    \
  UTILS_ARRAY_TEMPLATE_T_##N struct implementation_level_impl<                 \
      const UTILS_ARRAY_CONTAINER_T_##N(Container)> {                          \
    typedef mpl::integral_c_tag tag;                                           \
    typedef mpl::int_<ImplementationLevel> type;                               \
    BOOST_STATIC_CONSTANT(int,                                                 \
                          value = implementation_level_impl::type::value);     \
  };                                                                           \
  } /* namespace serialization */                                              \
  } /* namespace boost */

/**
 * @brief Redefinition of @c BOOST_CLASS_TRACKING for array types.
 * @tparam Container            Template template type of the array
 * @tparam N                    N if @p Container uses std::size_t N, else 0
 * @tparam TrackingLevel        Tracking level
 */
#define UTILS_ARRAY_BOOST_TRACK(Container, N, TrackingLevel)                   \
  namespace boost {                                                            \
  namespace serialization {                                                    \
  UTILS_ARRAY_TEMPLATE_T_##N struct tracking_level<                            \
      UTILS_ARRAY_CONTAINER_T_##N(Container)> {                                \
    typedef mpl::integral_c_tag tag;                                           \
    typedef mpl::int_<TrackingLevel> type;                                     \
    BOOST_STATIC_CONSTANT(int, value = tracking_level::type::value);           \
    BOOST_STATIC_ASSERT(                                                       \
        (mpl::greater<                                                         \
            implementation_level<UTILS_ARRAY_CONTAINER_T_##N(Container)>,      \
            mpl::int_<primitive_type>>::value));                               \
  };                                                                           \
  } /* namespace serialization */                                              \
  } /* namespace boost */

#endif

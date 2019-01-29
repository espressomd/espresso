#ifndef UTILS_VECTOR_TRAITS_HPP
#define UTILS_VECTOR_TRAITS_HPP

#include <cstddef>

#include <boost/type_traits/is_arithmetic.hpp>

#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/tags.hpp>

#include "Vector.hpp"

namespace boost {
namespace geometry {
namespace traits {
namespace detail {

// Create class and specialization to indicate the tag
// for normal cases and the case that the type of the std-array is arithmetic
template <bool> struct Vector_tag { typedef geometry_not_recognized_tag type; };

template <> struct Vector_tag<true> { typedef point_tag type; };

} // namespace detail

// Assign the point-tag, preventing arrays of points getting a point-tag
template <typename CoordinateType, std::size_t DimensionCount>
struct tag<Vector<DimensionCount, CoordinateType>>
    : detail::Vector_tag<boost::is_arithmetic<CoordinateType>::value> {};

template <typename CoordinateType, std::size_t DimensionCount>
struct coordinate_type<Vector<DimensionCount, CoordinateType>> {
  typedef CoordinateType type;
};

template <typename CoordinateType, std::size_t DimensionCount>
struct dimension<Vector<DimensionCount, CoordinateType>>
    : boost::mpl::int_<DimensionCount> {};

template <typename CoordinateType, std::size_t DimensionCount,
          std::size_t Dimension>
struct access<Vector<DimensionCount, CoordinateType>, Dimension> {
  static inline CoordinateType
  get(Vector<DimensionCount, CoordinateType> const &a) {
    return a[Dimension];
  }

  static inline void set(Vector<DimensionCount, CoordinateType> &a,
                         CoordinateType const &value) {
    a[Dimension] = value;
  }
};

template <class T, std::size_t N> struct coordinate_system<Vector<N, T>> {
  using type = cs::cartesian;
};

} // namespace traits
} // namespace geometry
} // namespace boost
#endif
